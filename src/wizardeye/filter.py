# -*- coding: utf-8 -*-

"""Filter alignment BAM utilities for WizardEye.

This Python module provides functions to generate mask using previously generated
cross-mappability tracks and to filter BAM with such masks. It also provides functions
to apply stringency-threshold, to compute raw overlapping-k-mers count per read and
to generate the final filtration report.
"""

from __future__ import annotations
import hashlib
import pyBigWig
import shutil
import subprocess
import tempfile
import numpy as np

from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

from .db import get_tracks, Track
import pysam
from .utils import (
    log,
    from_charlist_to_list,
    validate_bam_compatibility,
    merge_bed_files,
    convert_bigwig_to_bedGraph,
    BWAParameters,
)


# -- Mask creation related functions --


def _filter_bed_by_frequency(bed_path: Path, min_freq: int) -> None:
    """Filter a merged BED file to only keep positions overlapped by at least min_freq tracks.

    The BED file format has each line as: chrom, start, end, track_names
    This function modifies the file to keep only positions with >= min_freq tracks.

    Args:
        bed_path (Path): Path to the BED file to filter.
        min_freq (int): Minimum number of tracks required.
    """
    with open(bed_path, "r") as f:
        lines = f.readlines()

    filtered_lines = []
    for line in lines:
        if line.strip() == "":
            continue
        parts = line.strip().split("\t")
        if len(parts) < 4:
            continue
        track_names = [name for name in parts[3].split(",") if name.strip()]
        if len(track_names) >= min_freq:
            filtered_lines.append(line)

    with open(bed_path, "w") as f:
        f.writelines(filtered_lines)


# -- Mask creation related functions --


def _build_mask_from_track(
    track: Track,
    cross_stringency: float,
    consider_all: bool = False,
    no_cache: bool = False,
) -> Tuple[str, Path]:
    """Compute one track mask.

    Args:
            track (Track): the Track object to process.
            cross_stringency (float): Cross-stringency threshold, must be between 0.0 and 1.0.
            consider_all (bool): If True, all k-mers are considered in the cross-stringency
            no_cache (bool): If True, do not use and generate cached tracks.

    Returns:
            Tuple[str, Path]: A tuple containing the track name and the path to the BED mask.
    """
    track_dir = track.track_dir
    track_name = track.track_name
    source_label = "all" if consider_all else "uniq"

    # Store overlap BED with the track data so lookup/save stay within the track folder.
    if no_cache:
        with tempfile.NamedTemporaryFile(
            prefix=f"wizardeye_{track_name}_", suffix=".bed", delete=False
        ) as tmp_handle:
            per_input_bed = Path(tmp_handle.name)
    else:
        per_input_bed = (
            track_dir / f"mask_{track_name}_s{cross_stringency}_{source_label}.bed"
        )

    if not no_cache and per_input_bed.exists():
        log(
            f"Reusing existing mask for {track_name}... \n\033[0;90m({per_input_bed})\033[0m",
            "I",
        )
        return track_name, per_input_bed

    source_bw_name = "map_all.bw" if consider_all else "map_uniq.bw"
    source_bw = track_dir / source_bw_name
    if not source_bw.exists():
        raise ValueError(
            f"Missing required {source_bw_name} for selected track '{track_name}' in {track_dir}"
        )

    with tempfile.NamedTemporaryFile(
        prefix=f"wizardeye_{track_name}_", suffix=".bg", delete=False
    ) as tmp_bg:
        tmp_input_bg = Path(tmp_bg.name)

    tmp_input_bg = convert_bigwig_to_bedGraph(source_bw, Path(tmp_bg.name))

    per_input_bed = compute_stringency_on_bedGraph(
        input_bg=tmp_input_bg,
        kmer_length=track.identity.parameters.kmer_length,
        offset_step=track.identity.parameters.offset_step,
        cross_stringency=cross_stringency,
        output_bed=per_input_bed,
    )

    if tmp_input_bg.exists():
        tmp_input_bg.unlink()

    log(f"Mask cached for '{track_name}': {per_input_bed}", "S")
    return track.identity.query_species, per_input_bed


def generate_global_mask(
    ref_species: str,
    inputs: List[str],
    kmer_length: int,
    offset_step: int,
    cross_stringency: float,
    consider_all: bool = False,
    output_file: Optional[str] = None,
    bwa_params: Optional[BWAParameters] = None,
    db_root: str = "database",
    no_cache: bool = False,
    n_threads: int = 1,
    min_freq: Optional[int] = None,
) -> Path:
    """Generate a mask based on requested tracks and parameters.

    The merged output is a BED4-formatted file where the fourth column contains
    comma-separated track names contributing to each merged region.

    Selection is performed based on the parameters used for track generations.

        Args:
                ref_species (str): Name of the reference species.
                inputs (List[str]): List of input track names to consider.
                kmer_length (int): Length of the k-mers used during track generations.
                offset_step (int): Step size for the offset used during track generations.
                cross_stringency (float): Cross-stringency threshold, must be between 0.0 and 1.0.
                consider_all (bool): If True, all k-mers are considered in the cross-stringency
            calculation. If False (default), only uniquely aligned k-mers are considered.
                bwa (Optional[BWAParameters]): BWA alignment parameters for track selection.
                db_root (str): Root directory of the database.
                no_cache (bool): If True, do not use and generate cached tracks.
                n_threads (int): Number of threads to use for parallel processing.
                min_freq (Optional[int]): Minimum number of tracks that must overlap a position
                    for it to be included in the mask (frequency filtering). If None, uses standard merging.

        Returns:
                Path: Path to the generated mask.

        Raises:
                ValueError: If any of the input parameters are invalid.
    """
    if kmer_length < 1:
        raise ValueError("kmer_length must be a positive integer")
    if offset_step < 1:
        raise ValueError("offset_step must be a positive integer")
    if cross_stringency <= 0.0:
        raise ValueError("The cross-stringency threshold must be a positive float")
    if n_threads < 1:
        raise ValueError("n_threads must be a positive integer")
    if min_freq is not None and min_freq < 1:
        raise ValueError("min_freq must be a positive integer")

    if cross_stringency > 1.0:
        log("Cross-stringency is superior to 1.0.", "W")

    if min_freq is not None:
        log(
            f"Frequency filtering enabled: positions must be overlapped by at least {min_freq} tracks to be masked.",
            "I",
        )

    requested_inputs = from_charlist_to_list(inputs)
    tracks = get_tracks(
        ref_species=ref_species,
        kmer_length=kmer_length,
        offset_step=offset_step,
        db_root=db_root,
        bwa_params=bwa_params,
    )

    if not tracks:
        raise ValueError(
            f"No track found for reference '{ref_species}' with k={kmer_length} and offset={offset_step}"
        )

    if requested_inputs:
        selected_tracks = []
        selected_track_names = set()
        available_track_names = sorted({t.track_name for t in tracks})
        for requested in requested_inputs:
            matches = [
                track
                for track in tracks
                if requested
                in {track.track_name, track.query_name, track.identity.query_species}
            ]
            if not matches:
                raise ValueError(
                    "Requested input/track not found for this reference/k/offset: "
                    f"{requested}. Available track names: {', '.join(available_track_names)}"
                )
            if len(matches) > 1 and not any(t.track_name == requested for t in matches):
                matching_names = ", ".join(sorted(t.track_name for t in matches))
                raise ValueError(
                    f"Ambiguous input '{requested}' for this reference/k/offset. "
                    f"Use one of these track names instead: {matching_names}"
                )

            # If requested matches an exact track name, keep only this one.
            if any(t.track_name == requested for t in matches):
                matches = [t for t in matches if t.track_name == requested]

            for match in matches:
                track_name = match.track_name
                if track_name in selected_track_names:
                    continue
                selected_tracks.append(match)
                selected_track_names.add(track_name)
    else:
        selected_tracks = tracks

    ref_dir = Path(db_root) / ref_species
    source_label = "all" if consider_all else "uniq"

    # Keep one deterministic merged-mask path per exact track selection/parameters.
    sorted_track_names = sorted({track.track_name for track in selected_tracks})
    track_fingerprint = hashlib.sha1(
        ",".join(sorted_track_names).encode("utf-8")
    ).hexdigest()[:12]

    if output_file:
        merged_mask_bed = Path(output_file)
    elif no_cache:
        with tempfile.NamedTemporaryFile(
            prefix="wizardeye_mask_", suffix=".bed", delete=False
        ) as tmp_handle:
            merged_mask_bed = Path(tmp_handle.name)
    else:
        freq_suffix = f"_f{min_freq}" if min_freq is not None else ""
        merged_mask_bed = ref_dir / (
            f"mask_s{cross_stringency}_k{kmer_length}_o{offset_step}"
            f"_{source_label}_t{len(sorted_track_names)}_{track_fingerprint}{freq_suffix}.bed"
        )

    if not no_cache and merged_mask_bed.exists():
        log(f"Reusing existing mask... \033[0;90m({merged_mask_bed})\033[0m", "I")
        return merged_mask_bed

    if no_cache:
        merged_mask_bed.parent.mkdir(parents=True, exist_ok=True)

    per_track_beds: List[Tuple[str, Path]] = []

    if n_threads > 1 and len(selected_tracks) > 1:
        max_workers = min(n_threads, len(selected_tracks))
        log(
            f"Computing overlap tracks in parallel ({len(selected_tracks)} tracks, {max_workers} workers)...",
            "I",
        )
        with ThreadPoolExecutor(max_workers=max_workers) as pool:
            futures = [
                pool.submit(
                    _build_mask_from_track,
                    track,
                    cross_stringency,
                    consider_all,
                    no_cache,
                )
                for track in selected_tracks
            ]
            for future in as_completed(futures):
                per_track_beds.append(future.result())
    else:
        for track in selected_tracks:
            per_track_beds.append(
                _build_mask_from_track(
                    track=track,
                    cross_stringency=cross_stringency,
                    consider_all=consider_all,
                    no_cache=no_cache,
                )
            )

    try:
        log(f"Merging {len(per_track_beds)} masks into one unique mask...", "I")
        merge_bed_files(per_track_beds, merged_mask_bed)

        if min_freq is not None:
            _filter_bed_by_frequency(merged_mask_bed, min_freq)
    finally:
        if no_cache:
            for _, per_track_bed in per_track_beds:
                try:
                    per_track_bed.unlink(missing_ok=True)
                except TypeError:
                    if per_track_bed.exists():
                        per_track_bed.unlink()
    log(f"Merged BED mask cached to {merged_mask_bed}", "S")
    return merged_mask_bed


def compute_stringency_on_bedGraph(
    input_bg: Path,
    kmer_length: int,
    offset_step: int,
    cross_stringency: float,
    output_bed: Path,
) -> Path:
    """Apply a stringency filter on a bedgraph file using AWK and export the result in the BED format.

    This criterion is defined as the number of overlapping k-mers on a position divided by the size of the k-mers:

    cross_stringency * (kmer_length / offset_step)

    In perfect conditions, a position is overlapped by exactly k k-mers. See more details in README.md.

    Args:
            input_bg (Path): Path to the input bedGraph.
            kmer_length (int): Length of the k-mers used to generate the bedGraph.
            offset_step (int): Step size used to generate the bedGraph.
            cross_stringency (float): Stringency threshold to apply.
            output_bed (Path): Path to the output BED file.

    Returns:
            Path: Path to the output BED file.

    Raises:
            RuntimeError: If bedtools or awk are not found.

    Notes:
            This function will probably be modified to avoid the use of awk in the future.

    See also:
            Heng Li's seqbility tool including the original stringency definition used in here,
                    at https://github.com/lh3/misc/tree/cc0f36a9a19f35765efb9387389d9f3a6756f08f/seq/seqbility.
    """
    bedtools = shutil.which("bedtools")
    awk = shutil.which("awk")
    if not (bedtools and awk):
        raise RuntimeError(
            "bedtools and awk are required to compute mask from map_uniq.bw"
        )

    threshold = cross_stringency * (kmer_length / offset_step)
    output_bed.parent.mkdir(parents=True, exist_ok=True)

    awk_script = f'OFS="\\t" {{ if (($4 + 0) >= {threshold:.12f}) print $1, $2, $3; }}'

    with input_bg.open("r", encoding="utf-8") as bg_in, output_bed.open(
        "w", encoding="utf-8"
    ) as out:
        log(f"{awk} '{awk_script}' {input_bg} | {bedtools} merge -i stdin", "C")
        p1 = subprocess.Popen(
            [awk, awk_script], stdin=bg_in, stdout=subprocess.PIPE, text=True
        )
        p2 = subprocess.Popen(
            [bedtools, "merge", "-i", "stdin"], stdin=p1.stdout, stdout=out, text=True
        )
        if p1.stdout is not None:
            p1.stdout.close()
        rc2 = p2.wait()
        rc1 = p1.wait()
        if rc1 != 0:
            raise subprocess.CalledProcessError(rc1, [awk, awk_script])
        if rc2 != 0:
            raise subprocess.CalledProcessError(rc2, [bedtools, "merge", "-i", "stdin"])
    return output_bed


# -- Filtration related functions --


def filter_bam(
    input_bam: str,
    ref: str,
    db_root: str,
    exclude_tracks: List[str],
    kmer_length: int,
    offset_step: int,
    bwa_params: Optional[BWAParameters] = None,
    stringency: float = 0.99,
    min_freq: Optional[int] = None,
    consider_all: bool = False,
    output_report_tsv: Optional[str] = None,
    export_bam: bool = False,
    output_filtered_bam: Optional[str] = None,
    output_excluded_bam: Optional[str] = None,
) -> Dict[str, object]:
    """Alternative BAM filtering function using pyBigWig directly instead of bedtools mask intersection.

    Inspired by _generate_count_only_report, this function filters reads by checking if the sum of
    k-mers overlapping each read position exceeds the stringency threshold for any track.

    Args:
            input_bam (str): Path to the input BAM file to filter.
            ref (str): Name of the reference species.
            exclude_tracks (List[str]): List of input track names to consider.
            kmer_length (int): Length of the k-mers used during track generations.
            offset_step (int): Step size for the offset used during track generations.
            bwa (Optional[BWAParameters]): BWA alignment parameters for track selection.
            stringency (float): Cross-stringency threshold, must be between 0.0 and 1.0.
            min_freq (Optional[int]): Minimum number of tracks that must overlap a position for it to be masked.
                If None, uses standard logic (any track overlapping).
            consider_all (bool): If True, all k-mers are considered. If False (default), only uniquely aligned k-mers are considered.
            output_report_tsv (Optional[str]): Path to the output TSV report file.
            export_bam (bool): If True, splits the input BAM into filtered and excluded BAM files.
            output_filtered_bam (Optional[str]): Path to the output BAM file containing filtered reads if export_bam is True.
            output_excluded_bam (Optional[str]): Path to the output BAM file containing excluded reads if export_bam is True.
            db_root (str): Root directory of the database.

    Returns:
            Dict[str, object]: Dictionary containing mask path (None for this alternative), filtered/excluded BAM paths, report path, and counts.
    """
    if pysam is None:
        raise RuntimeError("pysam is required for report generation and BAM filtering")

    # Parameters validation
    input_bam_path = Path(input_bam)
    if not input_bam_path.exists() or not input_bam_path.is_file():
        raise FileNotFoundError(f"Input BAM not found: {input_bam_path}")

    if kmer_length < 1:
        raise ValueError("kmer_length must be a positive integer")
    if offset_step < 1:
        raise ValueError("offset_step must be a positive integer")
    if not (0.0 <= stringency <= 1.0):
        raise ValueError("stringency must be between 0.0 and 1.0")
    if min_freq is not None and min_freq < 1:
        raise ValueError("min_freq must be a positive integer")

    normalized_tracks = from_charlist_to_list(exclude_tracks)
    if not normalized_tracks:
        raise ValueError("exclude_tracks cannot be empty")

    # Fail fast on initial-file incompatibility before any processing.
    ref_dir = Path(db_root) / ref
    reference_seq_sizes = ref_dir / f"{ref}.sizes"
    validate_bam_compatibility(
        input_bam_path,
        reference_seq_sizes,
        bwa_params=bwa_params,
    )

    sorted_tracks = sorted(set(normalized_tracks))

    tracks_for_params = get_tracks(
        ref_species=ref,
        kmer_length=kmer_length,
        offset_step=offset_step,
        db_root=db_root,
        bwa_params=bwa_params,
    )
    normalized_requested_tracks = set(sorted_tracks)
    track_to_tags: Dict[str, Set[str]] = {}
    for track in tracks_for_params:
        if track.track_name not in normalized_requested_tracks:
            continue
        tags = set(from_charlist_to_list(track.info.get("tags", []), lowercase=True))
        for key in {
            track.track_name,
            track.identity.query_species,
            track.query_name,
            track.track_name.split("_")[0],
        }:
            track_to_tags.setdefault(key, set()).update(tags)

    report_tsv = (
        Path(output_report_tsv)
        if output_report_tsv
        else _default_output_table(input_bam_path)
    )
    report_tsv.parent.mkdir(parents=True, exist_ok=True)

    selected_tracks = [
        track
        for track in tracks_for_params
        if track.track_name in normalized_requested_tracks
    ]
    if not selected_tracks:
        raise ValueError("No selected tracks available for filtering")

    # Calculate threshold once
    threshold = stringency * (kmer_length / offset_step)

    # Open all BigWig files
    opened_bws = []
    for track in selected_tracks:
        bw_path = track.map_all_bw if consider_all else track.map_uniq_bw
        if not bw_path.exists():
            raise FileNotFoundError(f"Missing BigWig: {bw_path}")
        opened_bws.append(pyBigWig.open(str(bw_path)))

    try:
        read_tracks: Dict[Tuple[str, str, int, int], Set[str]] = {}
        read_tags: Dict[Tuple[str, str, int, int], Set[str]] = {}
        excluded_reads: Set[Tuple[str, str, int, int]] = set()
        n_total_records = 0
        n_mapped_records = 0

        log("Identifying overlapping reads...", "I")

        with pysam.AlignmentFile(str(input_bam_path), "rb") as bam:
            for read in bam.fetch(until_eof=True):
                n_total_records += 1
                read_id = read.query_name or ""

                if not read_id or read.is_unmapped or read.reference_name is None:
                    read_key = (read_id, "", -1, -1)
                    read_tracks[read_key] = set()
                    read_tags[read_key] = set()
                    continue

                n_mapped_records += 1
                chrom = read.reference_name
                start = read.reference_start
                end = read.reference_end

                if start is None or end is None or end <= start:
                    read_key = (read_id, chrom or "", start or -1, end or -1)
                    read_tracks[read_key] = set()
                    read_tags[read_key] = set()
                    continue

                read_key = (read_id, chrom, start, end)
                read_tracks[read_key] = set()
                read_tags[read_key] = set()
                overlapping_tracks = {}
                # Check each track for overlap exceeding threshold (per bp, using max)
                for idx, bw in enumerate(opened_bws):
                    if min_freq is None:
                        track = selected_tracks[idx]
                        stat_res = bw.stats(chrom, start, end, type="max", exact=True)

                        # stat_res is (value,) tuple, we want the max value per bp in the region
                        kmer_max = (
                            stat_res[0] if stat_res and stat_res[0] is not None else 0.0
                        )

                        if kmer_max >= threshold:
                            excluded_reads.add((chrom, read_id, start, end))
                            track_name = track.identity.query_species
                            read_tracks[read_key].add(track_name)
                            read_tags[read_key].update(
                                track_to_tags.get(track_name, set())
                            )
                    else:
                        vals = np.array(bw.values(chrom, start, end), dtype=float)
                        if np.any(np.isnan(vals)):
                            raise ValueError(
                                f"Unexpected NaN values in BigWig for track {selected_tracks[idx].track_name} at {chrom}:{start}-{end}"
                            )
                        overlapping_tracks[idx] = vals >= threshold

                if min_freq is not None:
                    sum_overlapping = np.sum(
                        np.array(
                            [
                                arr
                                for arr in overlapping_tracks.values()
                                if arr is not None and len(arr) > 0
                            ]
                        ),
                        axis=0,
                    )
                    max_freq = np.max(sum_overlapping)
                    if max_freq >= min_freq:
                        print("+1")
                        excluded_reads.add((chrom, read_id, start, end))
                        max_idx = np.argmax(sum_overlapping)
                        for track in overlapping_tracks.keys():
                            if overlapping_tracks[track][max_idx]:
                                track_name = selected_tracks[
                                    int(track)
                                ].identity.query_species
                                read_tracks[read_key].add(track_name)
                                read_tags[read_key].update(
                                    track_to_tags.get(track_name, set())
                                )

        report_path = write_filtration_report(
            output_report_tsv=report_tsv, read_tracks=read_tracks
        )
        n_filtered = max(0, n_mapped_records - len(excluded_reads))

        filtered_bam: Optional[Path] = None
        excluded_bam: Optional[Path] = None

        if export_bam:
            filtered_bam = (
                Path(output_filtered_bam)
                if output_filtered_bam
                else _default_output_bam(input_bam_path, "filtered")
            )
            excluded_bam = (
                Path(output_excluded_bam)
                if output_excluded_bam
                else _default_output_bam(input_bam_path, "excluded")
            )
            log("Filtering BAM from excluded read IDs...", "I")
            filter_bam_from_reads_id(
                input_bam=input_bam_path,
                excluded_reads=excluded_reads,
                output_filtered_bam=filtered_bam,
                output_excluded_bam=excluded_bam,
            )

            for bam_path in (filtered_bam, excluded_bam):
                try:
                    pysam.index(str(bam_path))
                except Exception:
                    log(f"Could not index BAM (possibly unsorted): {bam_path}", "W")

        return {
            "mask": None,
            "filtered_bam": filtered_bam,
            "excluded_bam": excluded_bam,
            "report_tsv": report_path,
            "n_total": n_mapped_records,
            "n_filtered": n_filtered,
            "n_excluded": len(excluded_reads),
            "n_total_records": n_total_records,
        }

    finally:
        for bw in opened_bws:
            bw.close()


def filter_bam_from_reads_id(
    input_bam: Path,
    excluded_reads: Set[Tuple[str, str, int, int]],
    output_filtered_bam: Path,
    output_excluded_bam: Path,
) -> Tuple[int, int, int]:
    """Split BAM in one pysam pass using excluded read IDs with position info.

    Args:
            input_bam: Path to input BAM file.
            excluded_reads: Set of tuples (chrom, read_id, start, stop) to exclude.
            output_filtered_bam: Path to output BAM file with filtered reads.
            output_excluded_bam: Path to output BAM file with excluded reads.

    Returns:
            Tuple[int, int, int]: Number of total records, number of filtered records and number of excluded records."""

    log("Splitting BAM into excluded/filtered files based on filtration...", "I")

    if pysam is None:
        raise RuntimeError("pysam is required to filter BAM from read IDs")

    output_filtered_bam.parent.mkdir(parents=True, exist_ok=True)
    output_excluded_bam.parent.mkdir(parents=True, exist_ok=True)

    n_total_records = 0
    n_excluded_records = 0

    with pysam.AlignmentFile(str(input_bam), "rb") as bam:
        with pysam.AlignmentFile(
            str(output_filtered_bam), "wb", template=bam
        ) as filtered_handle:
            with pysam.AlignmentFile(
                str(output_excluded_bam), "wb", template=bam
            ) as excluded_handle:
                for read in bam.fetch(until_eof=True):
                    n_total_records += 1
                    read_id = read.query_name
                    chrom = read.reference_name
                    start = read.reference_start
                    end = read.reference_end
                    if (
                        read_id
                        and chrom
                        and start is not None
                        and end is not None
                        and (chrom, read_id, start, end) in excluded_reads
                    ):
                        excluded_handle.write(read)
                        n_excluded_records += 1
                    else:
                        filtered_handle.write(read)

    n_filtered_records = n_total_records - n_excluded_records
    return n_total_records, n_filtered_records, n_excluded_records


# --- Count-report related function ---


def count_k_mers_on_bam(
    input_bam: str,
    ref: str,
    db_root: str,
    exclude_tracks: List[str],
    kmer_length: int,
    offset_step: int,
    count_mode: str,
    bwa_params: Optional[BWAParameters] = None,
    consider_all: bool = False,
    output_report_tsv: Optional[str] = None,
    n_threads: int = 1,
) -> Dict[str, object]:
    """Main function to filter a BAM file using several tracks and generate requested outputs.

    Args:
            input_bam (str): Path to the input BAM file to filter.
            ref (str): Name of the reference species.
            exclude_tracks (List[str]): List of input track names to consider.
            kmer_length (int): Length of the k-mers used during track generations.
            offset_step (int): Step size for the offset used during track generations.
            count_mode (str): Type of statistic to compute with count among mean, std, max, min, cov or sum.
            bwa (Optional[BWAParameters]): BWA alignment parameters for track selection.
            consider_all (bool): If True, all k-mers are considered in the cross-stringency
               calculation. If False (default), only uniquely aligned k-mers are considered.
            output_report_tsv (Optional[str]): Path to the output TSV report file.
            db_root (str): Root directory of the database.
            no_cache (bool): If True, do not use and generate cached tracks.
            n_threads (int): Number of threads to use for parallel processing.
    """
    if pysam is None:
        raise RuntimeError("pysam is required for report generation and BAM filtering")

    # Parameters validation
    input_bam_path = Path(input_bam)
    if not input_bam_path.exists() or not input_bam_path.is_file():
        raise FileNotFoundError(f"Input BAM not found: {input_bam_path}")

    if kmer_length < 1:
        raise ValueError("kmer_length must be a positive integer")
    if offset_step < 1:
        raise ValueError("offset_step must be a positive integer")

    normalized_tracks = from_charlist_to_list(exclude_tracks)
    if not normalized_tracks:
        raise ValueError("exclude_tracks cannot be empty")

    # Fail fast on initial-file incompatibility before any mask computation.
    ref_dir = Path(db_root) / ref
    reference_seq_sizes = ref_dir / f"{ref}.sizes"
    validate_bam_compatibility(
        input_bam_path,
        reference_seq_sizes,
        bwa_params=bwa_params,
    )

    sorted_tracks = sorted(set(normalized_tracks))

    tracks_for_params = get_tracks(
        ref_species=ref,
        kmer_length=kmer_length,
        offset_step=offset_step,
        db_root=db_root,
        bwa_params=bwa_params,
    )
    normalized_requested_tracks = set(sorted_tracks)
    track_to_tags: Dict[str, Set[str]] = {}
    for track in tracks_for_params:
        if track.track_name not in normalized_requested_tracks:
            continue
        tags = set(from_charlist_to_list(track.info.get("tags", []), lowercase=True))
        for key in {
            track.track_name,
            track.identity.query_species,
            track.query_name,
            track.track_name.split("_k")[0],
        }:
            track_to_tags.setdefault(key, set()).update(tags)

    report_tsv = (
        Path(output_report_tsv)
        if output_report_tsv
        else _default_output_count_table(input_bam_path)
    )
    report_tsv.parent.mkdir(parents=True, exist_ok=True)

    selected_tracks = [
        track
        for track in tracks_for_params
        if track.track_name in normalized_requested_tracks
    ]
    if not selected_tracks:
        raise ValueError("No selected tracks available to compute count-only report")

    n_total_records, n_mapped_records = _generate_count_only_report(
        input_bam=input_bam_path,
        selected_tracks=selected_tracks,
        consider_all=consider_all,
        count_mode=count_mode,
        output_report_tsv=report_tsv,
    )
    return {
        "report_tsv": report_tsv,
        "n_total_records": n_total_records,
        "n_processed": n_mapped_records,
    }


def _generate_count_only_report(
    input_bam: Path,
    selected_tracks: List[Track],
    count_mode: str,
    consider_all: bool,
    output_report_tsv: Path,
) -> Tuple[int, int]:
    """
    Generate a report summarizing statistics about selected tracks k-mers that can overlap reads found in a BAM file.

        Args:
                input_bam (Path): Path to the input BAM file.
                selected_tracks (List[Any]): List of selected tracks.
                count_mode (str): Type of statistic to compute with count among mean, std, max, min, cov or sum.
                consider_all (bool): Whether to consider all k-mers or only uniquely mapped ones.
                output_report_tsv (Path): Path to the output report TSV file.

        Returns:
                int: the number of read records in the BAM file
                int: the number of processed reads
    """
    track_names = [
        f"{count_mode}_{track.identity.query_species}" for track in selected_tracks
    ]

    opened_bws = []
    for track in selected_tracks:
        bw_path = track.map_all_bw if consider_all else track.map_uniq_bw
        if not bw_path.exists():
            raise FileNotFoundError(f"Missing BigWig: {bw_path}")
        opened_bws.append(pyBigWig.open(str(bw_path)))

    n_total_records = 0
    n_mapped_records = 0

    try:
        with output_report_tsv.open("w", encoding="utf-8") as out_tsv:
            header = ["read_id"] + track_names
            out_tsv.write("\t".join(header) + "\n")

            with pysam.AlignmentFile(str(input_bam), "rb") as bam:
                for read in bam.fetch(until_eof=True):
                    n_total_records += 1

                    read_id = read.query_name or ""
                    if not read_id or read.is_unmapped or read.reference_name is None:
                        continue

                    chrom = read.reference_name
                    start = read.reference_start
                    end = read.reference_end

                    if start is None or end is None or end <= start:
                        continue

                    n_mapped_records += 1
                    row_values = [read_id]

                    for bw in opened_bws:
                        stat_res = bw.stats(
                            chrom, start, end, type=count_mode, exact=True
                        )
                        row_values.append(
                            f"{stat_res[0]:.0f}"
                            if count_mode in ["max", "min", "sum"]
                            else f"{stat_res[0]:.3f}"
                        )

                    out_tsv.write("\t".join(row_values) + "\n")

    finally:
        for bw in opened_bws:
            bw.close()

    return n_total_records, n_mapped_records


# --- Reports and results generation ---


def _default_output_bam(input_bam: Path, suffix: str) -> Path:
    if input_bam.suffix.lower() == ".bam":
        return input_bam.with_suffix(f".{suffix}.bam")
    return Path(f"{str(input_bam)}.{suffix}.bam")


def _default_output_table(input_bam: Path) -> Path:
    if input_bam.suffix.lower() == ".bam":
        return input_bam.with_suffix(".wizardeye.report.tsv")
    return Path(f"{str(input_bam)}.wizardeye.report.tsv")


def _default_output_count_table(input_bam: Path) -> Path:
    if input_bam.suffix.lower() == ".bam":
        return input_bam.with_suffix(".wizardeye.counts.tsv")
    return Path(f"{str(input_bam)}.wizardeye.counts.tsv")


def write_filtration_report(
    output_report_tsv: Path, read_tracks: Dict[Tuple[str, str, int, int], Set[str]]
) -> Path:
    """Write one line per read with exclusion flag, overlapping tracks and tags.

    Args:
            output_report_tsv (Path): Path to the output TSV file.
            read_tracks (Dict[Tuple[str, str, int, int], Set[str]]): Dictionary mapping
                (read_id, chrom, start, end) tuples to sets of overlapping tracks.

    Returns:
            Path: The path to the generated report."""
    output_report_tsv.parent.mkdir(parents=True, exist_ok=True)

    # Check for duplicate read IDs (same read_id with different positions)
    read_id_to_positions: Dict[str, Set[Tuple[str, int, int]]] = {}
    for read_key in read_tracks.keys():
        if isinstance(read_key, tuple) and len(read_key) >= 4:
            rid, chrom, start, end = read_key
            if rid not in read_id_to_positions:
                read_id_to_positions[rid] = set()
            read_id_to_positions[rid].add((chrom, start, end))

    # Log warning if duplicates exist
    duplicate_read_ids = [
        rid for rid, positions in read_id_to_positions.items() if len(positions) > 1
    ]
    if duplicate_read_ids:
        log(
            f"Found {len(duplicate_read_ids)} duplicated read IDs: {duplicate_read_ids}. "
            f"Using ID:chrom:start:end format in report for all reads.",
            "W",
        )

    with output_report_tsv.open("w", encoding="utf-8") as handle:
        handle.write("read_key\tfiltered_out\tassociated_tracks\n")
        for read_key, track_set in read_tracks.items():
            tracks = sorted(track_set)
            excluded = "true" if tracks else "false"
            overlapped = ",".join(tracks) if tracks else ""
            # read_key is (read_id, chrom, start, end) where start is 0-based and end is 0-based exclusive
            if isinstance(read_key, tuple) and len(read_key) >= 4:
                rid, chrom, start, end = read_key
                # Use read_id:chrom:start:end format for ALL reads (1-based positions)
                # start + 1 converts 0-based to 1-based start
                # end is 0-based exclusive, which equals 1-based end (inclusive)
                display_id = f"{rid}:{chrom}:{start + 1}:{end}"
            else:
                display_id = str(read_key)
            # Build line with exactly 3 tab-separated columns
            line_parts = [display_id, excluded, overlapped]
            handle.write("\t".join(line_parts) + "\n")
    return output_report_tsv
