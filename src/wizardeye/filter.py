# -*- coding: utf-8 -*-

from __future__ import annotations

import hashlib
import bisect
import shutil
import subprocess
import tempfile

from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Any, Dict, List, Optional, Set, Tuple

from .db import get_tracks
from .mappability import pysam
from .utils import (
    log,
    from_charlist_to_list,
    validate_initial_bam_reference_compatibility,
    merge_bed_files,
    append_merged_interval,
)

# -- Filtration-mask creation related functions --

def generate_mask(
    ref_species: str,
    inputs: Optional[List[str]],
    kmer_length: int,
    offset_step: int,
    cross_stringency: float,
    consider_all: bool = False,
    n_threads: int = 1,
    db_root: str = "database",
    output_file: Optional[str] = None,
    bwa_missing_prob_err_rate: Optional[float] = None,
    bwa_max_gap_opens: Optional[int] = None,
    bwa_seed_length: Optional[int] = None,
    no_cache: bool = False,
) -> Path:
    """
    Generate overlap BED files for selected inputs and one merged mask BED.

    The merged output is a BED4 where column 4 contains comma-separated
    track names contributing to each merged region.

    Selection is done on:
    - reference
    - selected inputs (or all inputs if omitted)
    - k-mer length
    - offset step
    """
    if kmer_length < 1:
        raise ValueError("kmer_length must be a positive integer")
    if offset_step < 1:
        raise ValueError("offset_step must be a positive integer")
    if not (0.0 <= cross_stringency <= 1.0):
        raise ValueError("cross_stringency must be between 0.0 and 1.0")
    if n_threads < 1:
        raise ValueError("n_threads must be a positive integer")

    requested_inputs = from_charlist_to_list(inputs)
    tracks = get_tracks(
        ref_species=ref_species,
        kmer_length=kmer_length,
        offset_step=offset_step,
        db_root=db_root,
        bwa_missing_prob_err_rate=bwa_missing_prob_err_rate,
        bwa_max_gap_opens=bwa_max_gap_opens,
        bwa_seed_length=bwa_seed_length,
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
                if requested in {track.track_name, track.query_name, track.identity.query_species}
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
    stringency_label = str(cross_stringency)
    source_label = "all" if consider_all else "uniq"

    # Keep one deterministic merged-mask path per exact track selection/parameters.
    sorted_track_names = sorted({track.track_name for track in selected_tracks})
    track_fingerprint = hashlib.sha1(",".join(sorted_track_names).encode("utf-8")).hexdigest()[:12]

    if output_file:
        merged_mask_bed = Path(output_file)
    elif no_cache:
        with tempfile.NamedTemporaryFile(prefix="wizardeye_mask_", suffix=".bed", delete=False) as tmp_handle:
            merged_mask_bed = Path(tmp_handle.name)
    else:
        merged_mask_bed = ref_dir / (
            f"mask_s{stringency_label}_k{kmer_length}_o{offset_step}"
            f"_{source_label}_t{len(sorted_track_names)}_{track_fingerprint}.bed"
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
                    stringency_label,
                    kmer_length,
                    offset_step,
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
            per_track_beds.append(_build_mask_from_track(
                track=track,
                stringency_label=stringency_label,
                kmer_length=kmer_length,
                offset_step=offset_step,
                cross_stringency=cross_stringency,
                consider_all=consider_all,
                no_cache=no_cache,
            ))

    try:
        log(f"Merging {len(per_track_beds)} BED masks into one unique mask...", "I")
        merge_bed_files(per_track_beds, merged_mask_bed)
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

def compute_stringency(
	uniq_track: Path,
	kmer_length: int,
	offset_step: int,
	cross_stringency: float,
	output_bed: Path,
) -> Path:
	"""
	Compute a stringency BED mask from a unique-depth track.

	Parameters
	----------
	uniq_track:
		Path to a unique-depth track. Supported formats: .bw or .bg/.bedgraph
	kmer_length:
		K-mer length used during generation.
	offset_step:
		Sliding window step used during generation.
	cross_stringency:
		Threshold in [0.0, 1.0] applied on unique_depth / (kmer_length / offset_step).
	output_bed:
		Output BED file path where merged intervals are written.
	"""
	uniq_track = Path(uniq_track)
	output_bed = Path(output_bed)

	if kmer_length < 1:
		raise ValueError("kmer_length must be a positive integer")
	if offset_step < 1:
		raise ValueError("offset_step must be a positive integer")
	if not (0.0 <= cross_stringency <= 1.0):
		raise ValueError("cross_stringency must be between 0.0 and 1.0")
	if not uniq_track.exists():
		raise FileNotFoundError(f"Unique track not found: {uniq_track}")

	threshold = cross_stringency * (kmer_length / offset_step)
	tmp_bg_path: Optional[Path] = None

	if uniq_track.suffix.lower() in {".bg", ".bedgraph"}:
		bg_path = uniq_track
	elif uniq_track.suffix.lower() == ".bw":
		tool = shutil.which("bigWigToBedGraph")
		if tool is None:
			raise RuntimeError(
				"bigWigToBedGraph not found in PATH, cannot compute stringency from .bw input"
			)
		with tempfile.NamedTemporaryFile(suffix=".bg", delete=False) as tmp_handle:
			tmp_bg_path = Path(tmp_handle.name)
		bg_path = tmp_bg_path
		log(f"{tool} {str(uniq_track)} {str(bg_path)}", "C")
		subprocess.run([tool, str(uniq_track), str(bg_path)], check=True)
	else:
		raise ValueError(
			f"Unsupported track format for compute_stringency: {uniq_track.suffix}. "
			"Expected .bw, .bg or .bedgraph"
		)

	try:
		selected: List[Tuple[str, int, int]] = []
		with bg_path.open("r", encoding="utf-8") as handle:
			for raw_line in handle:
				line = raw_line.strip()
				if not line:
					continue
				parts = line.split("\t")
				if len(parts) < 4:
					continue
				chrom = parts[0]
				start = int(parts[1])
				end = int(parts[2])
				depth = float(parts[3])
				if end <= start:
					continue
				if depth >= threshold:
					selected.append((chrom, start, end))

		selected.sort(key=lambda interval: (interval[0], interval[1], interval[2]))
		merged: List[Tuple[str, int, int]] = []
		for chrom, start, end in selected:
			append_merged_interval(merged, chrom, start, end)

		output_bed.parent.mkdir(parents=True, exist_ok=True)
		with output_bed.open("w", encoding="utf-8") as handle:
			for chrom, start, end in merged:
				handle.write(f"{chrom}\t{start}\t{end}\n")

		log(
			f"Stringency mask written to: {output_bed} "
			f"(threshold unique-depth >= {threshold:.4f})",
			"S",
		)
		return output_bed
	finally:
		if tmp_bg_path is not None and tmp_bg_path.exists():
			tmp_bg_path.unlink()


def _default_output_bam(input_bam: Path, suffix: str) -> Path:
    if input_bam.suffix.lower() == ".bam":
        return input_bam.with_suffix(f".{suffix}.bam")
    return Path(f"{str(input_bam)}.{suffix}.bam")


def _default_output_table(input_bam: Path) -> Path:
    if input_bam.suffix.lower() == ".bam":
        return input_bam.with_suffix(".wizardeye.report.tsv")
    return Path(f"{str(input_bam)}.wizardeye.report.tsv")

def _build_mask_from_track(
    track: Any,
    stringency_label: str,
    kmer_length: int,
    offset_step: int,
    cross_stringency: float,
    consider_all: bool = False,
    no_cache: bool = False,
) -> Tuple[str, Path]:
    """Compute one track BED mask for a selected input and return (track name, BED path)."""
    track_dir = track.track_dir
    track_name = track.track_name
    source_label = "all" if consider_all else "uniq"

    # Store overlap BED with the track data so lookup/save stay within the track folder.
    if no_cache:
        with tempfile.NamedTemporaryFile(prefix=f"wizardeye_{track_name}_", suffix=".bed", delete=False) as tmp_handle:
            per_input_bed = Path(tmp_handle.name)
    else:
        per_input_bed = track_dir / f"mask_{track_name}_s{stringency_label}_{source_label}.bed"

    if not no_cache and per_input_bed.exists():
        log(f"Reusing existing mask for {track_name}... \n\033[0;90m({per_input_bed})\033[0m", "I")
        return track_name, per_input_bed

    source_bw_name = "map_all.bw" if consider_all else "map_uniq.bw"
    source_bw = track_dir / source_bw_name
    if not source_bw.exists():
        raise ValueError(
            f"Missing required {source_bw_name} for selected track '{track_name}' in {track_dir}"
        )

    compute_mask_bed_from_bigwig(
        source_bw=source_bw,
        kmer_length=kmer_length,
        offset_step=offset_step,
        cross_stringency=cross_stringency,
        output_bed=per_input_bed,
    )

    log(f"BED mask cached for '{track_name}': {per_input_bed}", "S")
    return track_name, per_input_bed

def compute_mask_bed_from_bigwig(
    source_bw: Path,
    kmer_length: int,
    offset_step: int,
    cross_stringency: float,
    output_bed: Path,
) -> Path:
    """Compute one overlap BED directly from BigWig with a fast awk/bedtools path when available."""
    tool = shutil.which("bigWigToBedGraph")
    if tool is None:
        log("bigWigToBedGraph not found in PATH and BAM file is unavailable: cannot compute overlaps.", "FE")

    output_bed.parent.mkdir(parents=True, exist_ok=True)

    with tempfile.NamedTemporaryFile(prefix="wizardeye_", suffix=".bg", delete=False) as tmp_handle:
        tmp_bg_path = Path(tmp_handle.name)

    try:
        log(f"Converting track to bedGraph for overlap computation: {source_bw}", "I")
        log(f"{tool} {str(source_bw)} {str(tmp_bg_path)}", "C")
        subprocess.run([tool, str(source_bw), str(tmp_bg_path)], check=True)

        threshold = cross_stringency * (kmer_length / offset_step)
        bedtools = shutil.which("bedtools")
        awk = shutil.which("awk")
        if not (bedtools and awk):
            raise RuntimeError(
                "bedtools and awk are required to compute mask from map_uniq.bw"
            )

        awk_script = f'OFS="\\t" {{ if (($4 + 0) >= {threshold:.12f}) print $1, $2, $3; }}'
        with tmp_bg_path.open("r", encoding="utf-8") as bg_in, output_bed.open("w", encoding="utf-8") as out:
            log(f"{awk} '{awk_script}' {tmp_bg_path} | {bedtools} merge -i stdin", "C")
            p1 = subprocess.Popen([awk, awk_script], stdin=bg_in, stdout=subprocess.PIPE, text=True)
            p2 = subprocess.Popen([bedtools, "merge", "-i", "stdin"], stdin=p1.stdout, stdout=out, text=True)
            if p1.stdout is not None:
                p1.stdout.close()
            rc2 = p2.wait()
            rc1 = p1.wait()
            if rc1 != 0:
                raise subprocess.CalledProcessError(rc1, [awk, awk_script])
            if rc2 != 0:
                raise subprocess.CalledProcessError(rc2, [bedtools, "merge", "-i", "stdin"])
        return output_bed
    finally:
        if tmp_bg_path.exists():
            tmp_bg_path.unlink()

def _load_mask_from_tracks(
    mask_path: Path,
) -> Dict[str, List[Tuple[int, int, Set[str]]]]:
    """Load BED mask and keep optional track labels from column 4."""
    by_chrom: Dict[str, List[Tuple[int, int, Set[str]]]] = {}
    with mask_path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 3:
                continue
            chrom = parts[0]
            start = int(parts[1])
            end = int(parts[2])
            if end <= start:
                continue
            tracks: Set[str] = set()
            if len(parts) > 3 and parts[3].strip():
                tracks = {track for track in parts[3].split(",") if track}
            by_chrom.setdefault(chrom, []).append((start, end, tracks))

    for chrom, intervals in by_chrom.items():
        intervals.sort(key=lambda item: (item[0], item[1]))
        by_chrom[chrom] = intervals
    return by_chrom


def _collect_overlapping_tracks(
    chrom: str,
    read_start: int,
    read_end: int,
    interval_index_by_chrom: Dict[str, Tuple[List[int], List[Tuple[int, int, Set[str]]]]],
) -> Set[str]:
    """Return all track labels overlapping one read interval on one contig."""
    if read_end <= read_start:
        return set()

    indexed = interval_index_by_chrom.get(chrom)
    if indexed is None:
        return set()

    starts, intervals = indexed
    idx = bisect.bisect_left(starts, read_end)
    overlapped: Set[str] = set()

    probe = idx - 1
    while probe >= 0:
        start, end, tracks = intervals[probe]
        if end <= read_start:
            break
        if start < read_end and end > read_start:
            overlapped.update(tracks)
        probe -= 1

    probe = idx
    while probe < len(intervals):
        start, end, tracks = intervals[probe]
        if start >= read_end:
            break
        if end > read_start:
            overlapped.update(tracks)
        probe += 1

    return overlapped


# -- Main filtration logic --

def filter_bam(
    input_bam: str,
    ref: str,
    db_root: str,
    exclude_tracks: List[str],
    kmer_length: int,
    offset_step: int,
    bwa_missing_prob_err_rate: Optional[float] = None,
    bwa_max_gap_opens: Optional[int] = None,
    bwa_seed_length: Optional[int] = None,
    stringency: float = 0.99,
    consider_all: bool = False,
    no_cache: bool = False,
    n_threads: int = 1,
    output_filtered_bam: Optional[str] = None,
    output_excluded_bam: Optional[str] = None,
    output_report_tsv: Optional[str] = None,
    export_bam: bool = False,
) -> Dict[str, object]:
    """Build a merged mask for selected tracks and filter one input BAM against it.

    If export_bam is enabled, writes two BAM files:
    - filtered: reads that do NOT overlap the generated mask
    - excluded: reads that DO overlap the generated mask
    Otherwise, only the read exclusion report is generated.
    """
    if pysam is None:
        raise RuntimeError("pysam is required for report generation and BAM filtering")

    input_bam_path = Path(input_bam)
    if not input_bam_path.exists() or not input_bam_path.is_file():
        raise FileNotFoundError(f"Input BAM not found: {input_bam_path}")

    if kmer_length < 1:
        raise ValueError("kmer_length must be a positive integer")
    if offset_step < 1:
        raise ValueError("offset_step must be a positive integer")
    if not (0.0 <= stringency <= 1.0):
        raise ValueError("stringency must be between 0.0 and 1.0")

    normalized_tracks = from_charlist_to_list(exclude_tracks)
    if not normalized_tracks:
        raise ValueError("exclude_tracks cannot be empty")

    # Fail fast on initial-file incompatibility before any mask computation.
    ref_dir = Path(db_root) / ref
    reference_seq_sizes = ref_dir / f"{ref}.sizes"
    validate_initial_bam_reference_compatibility(
        input_bam_path,
        reference_seq_sizes,
        bwa_missing_prob_err_rate=bwa_missing_prob_err_rate,
        bwa_max_gap_opens=bwa_max_gap_opens,
        bwa_seed_length=bwa_seed_length,
    )

    sorted_tracks = sorted(set(normalized_tracks))

    merged_mask = generate_mask(
        ref_species=ref,
        inputs=sorted_tracks,
        kmer_length=kmer_length,
        offset_step=offset_step,
        cross_stringency=stringency,
        consider_all=consider_all,
        n_threads=n_threads,
        db_root=db_root,
        output_file=None,
        bwa_missing_prob_err_rate=bwa_missing_prob_err_rate,
        bwa_max_gap_opens=bwa_max_gap_opens,
        bwa_seed_length=bwa_seed_length,
        no_cache=no_cache,
    )

    tracks_for_params = get_tracks(
        ref_species=ref,
        kmer_length=kmer_length,
        offset_step=offset_step,
        db_root=db_root,
        bwa_missing_prob_err_rate=bwa_missing_prob_err_rate,
        bwa_max_gap_opens=bwa_max_gap_opens,
        bwa_seed_length=bwa_seed_length,
    )
    normalized_requested_tracks = set(sorted_tracks)
    track_to_tags: Dict[str, Set[str]] = {}
    for track in tracks_for_params:
        if track.track_name not in normalized_requested_tracks:
            continue
        tags = set(from_charlist_to_list(track.info.get("tags", []), lowercase=True))
        for key in {track.track_name, track.identity.query_species, track.query_name, track.track_name.split("_k")[0]}:
            track_to_tags.setdefault(key, set()).update(tags)

    report_tsv = Path(output_report_tsv) if output_report_tsv else _default_output_table(input_bam_path)
    report_tsv.parent.mkdir(parents=True, exist_ok=True)

    filtered_bam: Optional[Path] = None
    excluded_bam: Optional[Path] = None
    try:
        read_tracks, read_tags, excluded_read_ids, n_total, n_excluded, n_total_records = _filter_bam_with_mask(
            input_bam=input_bam_path,
            mask_path=merged_mask,
            track_to_tags=track_to_tags,
        )
        report_path = write_filtration_report(
            output_report_tsv=report_tsv,
            read_tracks=read_tracks,
            read_tags=read_tags,
        )
        n_filtered = max(0, n_total - n_excluded)

        if export_bam:
            filtered_bam = Path(output_filtered_bam) if output_filtered_bam else _default_output_bam(input_bam_path, "filtered")
            excluded_bam = Path(output_excluded_bam) if output_excluded_bam else _default_output_bam(input_bam_path, "excluded")
            log("Filtering BAM from excluded read IDs...", "I")
            filter_bam_from_reads_id(
                input_bam=input_bam_path,
                excluded_read_ids=excluded_read_ids,
                output_filtered_bam=filtered_bam,
                output_excluded_bam=excluded_bam,
            )

            for bam_path in (filtered_bam, excluded_bam):
                try:
                    pysam.index(str(bam_path))
                except Exception:
                    log(f"Could not index BAM (possibly unsorted): {bam_path}", "W")
    finally:
        if no_cache:
            # Clean up any temporary mask generated when cache is disabled.
            try:
                merged_mask.unlink(missing_ok=True)
            except TypeError:
                if merged_mask.exists():
                    merged_mask.unlink()

    return {
        "mask": merged_mask,
        "filtered_bam": filtered_bam,
        "excluded_bam": excluded_bam,
        "report_tsv": report_path,
        "n_total": n_total,
        "n_filtered": n_filtered,
        "n_excluded": n_excluded,
        "n_total_records": n_total_records,
    }

# -- Report generation related functions --

def write_filtration_report(
    output_report_tsv: Path,
    read_tracks: Dict[str, Set[str]],
    read_tags: Dict[str, Set[str]],
) -> Path:
    """Write one line per read_id with exclusion flag, overlapping tracks and tags."""
    output_report_tsv.parent.mkdir(parents=True, exist_ok=True)
    with output_report_tsv.open("w", encoding="utf-8") as handle:
        handle.write("read_id\texcluded\toverlapped\ttags\n")
        for read_id, track_set in read_tracks.items():
            tracks = sorted(track_set)
            tags = sorted(read_tags.get(read_id, set()))
            excluded = "true" if tracks else "false"
            overlapped = ",".join(tracks)
            joined_tags = ",".join(tags)
            handle.write(f"{read_id}\t{excluded}\t{overlapped}\t{joined_tags}\n")
    return output_report_tsv

# -- Filtration related functions --

def _filter_bam_with_mask(
    input_bam: Path,
    mask_path: Path,
    track_to_tags: Optional[Dict[str, Set[str]]] = None,
) -> Tuple[Dict[str, Set[str]], Dict[str, Set[str]], Set[str], int, int, int]:
    """Compute per-read overlaps/tags and read-id exclusion set in one BAM pass."""
    if pysam is None:
        raise RuntimeError("pysam is required to compute read-level overlap data")

    log("Identifying overlapped reads...", "I")
    track_to_tags = track_to_tags or {}
    intervals_by_chrom = _load_mask_from_tracks(mask_path)
    interval_index_by_chrom: Dict[str, Tuple[List[int], List[Tuple[int, int, Set[str]]]]] = {
        chrom: ([interval[0] for interval in intervals], intervals)
        for chrom, intervals in intervals_by_chrom.items()
    }
    read_tracks: Dict[str, Set[str]] = {}
    read_tags: Dict[str, Set[str]] = {}
    excluded_read_ids: Set[str] = set()
    n_total_records = 0
    n_excluded_records = 0

    track_to_tags_get = track_to_tags.get

    with pysam.AlignmentFile(str(input_bam), "rb") as bam:
        for read in bam.fetch(until_eof=True):
            n_total_records += 1
            read_id = read.query_name
            if not read_id:
                continue

            if read_id not in read_tracks:
                read_tracks[read_id] = set()
                read_tags[read_id] = set()

            if read.is_unmapped or read.reference_name is None:
                continue

            read_start = read.reference_start
            read_end = read.reference_end
            if read_start is None or read_end is None:
                continue

            read_chrom = read.reference_name

            overlapped_tracks = _collect_overlapping_tracks(
                chrom=read_chrom,
                read_start=read_start,
                read_end=read_end,
                interval_index_by_chrom=interval_index_by_chrom,
            )
            if overlapped_tracks:
                excluded_read_ids.add(read_id)
                n_excluded_records += 1
            read_tracks[read_id].update(overlapped_tracks)
            for track_name in overlapped_tracks:
                read_tags[read_id].update(track_to_tags_get(track_name, set()))

            if n_total_records % 100000 == 0:
                log(f"\t\t{n_total_records} reads processed...", "I")

    n_total_reads = len(read_tracks)
    n_excluded_reads = len(excluded_read_ids)
    return (
        read_tracks,
        read_tags,
        excluded_read_ids,
        n_total_reads,
        n_excluded_reads,
        n_total_records,
    )


def filter_bam_from_reads_id(
    input_bam: Path,
    excluded_read_ids: Set[str],
    output_filtered_bam: Path,
    output_excluded_bam: Path,
) -> Tuple[int, int, int]:
    """Split BAM in one pysam pass using excluded read IDs."""

    log(f"Splitting BAM into excluded/filtered files based on filtration...", "I")

    if pysam is None:
        raise RuntimeError("pysam is required to filter BAM from read IDs")

    output_filtered_bam.parent.mkdir(parents=True, exist_ok=True)
    output_excluded_bam.parent.mkdir(parents=True, exist_ok=True)

    n_total_records = 0
    n_excluded_records = 0

    with pysam.AlignmentFile(str(input_bam), "rb") as bam:
        with pysam.AlignmentFile(str(output_filtered_bam), "wb", template=bam) as filtered_handle:
            with pysam.AlignmentFile(str(output_excluded_bam), "wb", template=bam) as excluded_handle:
                for read in bam.fetch(until_eof=True):
                    n_total_records += 1
                    read_id = read.query_name
                    if read_id and read_id in excluded_read_ids:
                        excluded_handle.write(read)
                        n_excluded_records += 1
                    else:
                        filtered_handle.write(read)

    n_filtered_records = n_total_records - n_excluded_records
    return n_total_records, n_filtered_records, n_excluded_records

