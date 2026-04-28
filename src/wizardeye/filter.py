# -*- coding: utf-8 -*-

from __future__ import annotations

import hashlib
import bisect
import shutil
import subprocess
import tempfile
import math

from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Set, Tuple

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


def _default_output_count_table(input_bam: Path) -> Path:
    if input_bam.suffix.lower() == ".bam":
        return input_bam.with_suffix(".wizardeye.counts.tsv")
    return Path(f"{str(input_bam)}.wizardeye.counts.tsv")

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

def _format_count_value(value: float) -> str:
    if math.isfinite(value) and abs(value - round(value)) <= 1e-9:
        return str(int(round(value)))
    return f"{value:.6f}".rstrip("0").rstrip(".")


def write_count_only_report(
    output_report_tsv: Path,
    rows: Iterable[Tuple[str, str, int, int, Dict[str, float]]],
    track_columns: List[str],
) -> Path:
    """Write count-only report with one alignment line and one column per selected track."""
    output_report_tsv.parent.mkdir(parents=True, exist_ok=True)
    with output_report_tsv.open("w", encoding="utf-8") as handle:
        header = ["read_id", "chr", "start", "stop"] + track_columns
        handle.write("\t".join(header) + "\n")
        for read_id, chrom, start, end, counts in rows:
            values = [
                read_id,
                chrom,
                str(start),
                str(end),
            ]
            for track_name in track_columns:
                values.append(_format_count_value(counts.get(track_name, 0.0)))
            handle.write("\t".join(values) + "\n")
    return output_report_tsv


def _generate_count_only_report(
    input_bam: Path,
    selected_tracks: List[Any],
    consider_all: bool,
    output_report_tsv: Path,
    n_threads: int,
) -> Dict[str, int]:
    """Generate per-alignment overlap counts from one selected BigWig source per track.

    This implementation uses:
    1) bigWigToBedGraph on selected source (`map_uniq.bw` or `map_all.bw`)
    2) bedtools unionbedg to align depth columns across tracks
    3) bedtools intersect to map aligned depth segments onto reads
    """
    if pysam is None:
        raise RuntimeError("pysam is required for count-only report generation")

    bigwig_to_bedgraph = shutil.which("bigWigToBedGraph")
    if bigwig_to_bedgraph is None:
        raise RuntimeError("bigWigToBedGraph not found in PATH, cannot produce count-only report")

    bedtools = shutil.which("bedtools")
    if bedtools is None:
        raise RuntimeError("bedtools not found in PATH, cannot produce count-only report")

    source_bw_name = "map_all.bw" if consider_all else "map_uniq.bw"
    log(f"Count-mode enabled. Note that this process may take some time.", "W")

    track_names: List[str] = []
    tmp_bg_files: List[Path] = []
    reads_bed_path: Optional[Path] = None
    tmp_extra_files: List[Path] = []

    def _convert_one_track_to_bg(track: Any) -> Tuple[str, Path]:
        track_name = track.track_name
        source_bw = track.map_all_bw if consider_all else track.map_uniq_bw
        if not source_bw.exists():
            raise FileNotFoundError(f"Missing required file for count-only mode: {source_bw}")

        with tempfile.NamedTemporaryFile(prefix=f"wizardeye_{track_name}_", suffix=".bg", delete=False) as tmp_handle:
            bg_path = Path(tmp_handle.name)

        log(f"{bigwig_to_bedgraph} {source_bw} {bg_path}", "C")
        subprocess.run([bigwig_to_bedgraph, str(source_bw), str(bg_path)], check=True)
        return track_name, bg_path

    try:
        if n_threads > 1 and len(selected_tracks) > 1:
            max_workers = min(n_threads, len(selected_tracks))
            log(
                f"Converting BigWig tracks in parallel for count-only mode ({len(selected_tracks)} tracks, {max_workers} workers)...",
                "I",
            )
            with ThreadPoolExecutor(max_workers=max_workers) as pool:
                futures = [pool.submit(_convert_one_track_to_bg, track) for track in selected_tracks]
                for future in as_completed(futures):
                    track_name, bg_path = future.result()
                    track_names.append(track_name)
                    tmp_bg_files.append(bg_path)
        else:
            for track in selected_tracks:
                track_name, bg_path = _convert_one_track_to_bg(track)
                track_names.append(track_name)
                tmp_bg_files.append(bg_path)

        with tempfile.NamedTemporaryFile(prefix="wizardeye_reads_", suffix=".bed", delete=False) as tmp_reads:
            reads_bed_path = Path(tmp_reads.name)

        records: List[Tuple[int, str, str, int, int]] = []
        read_contigs: Set[str] = set()
        rec_indices_by_contig: Dict[str, List[int]] = {}
        n_total_records = 0
        n_mapped_records = 0
        with reads_bed_path.open("w", encoding="utf-8") as reads_out:
            with pysam.AlignmentFile(str(input_bam), "rb") as bam:
                for read in bam.fetch(until_eof=True):
                    n_total_records += 1
                    read_id = read.query_name or ""
                    if not read_id:
                        continue
                    if read.is_unmapped or read.reference_name is None:
                        continue

                    read_start = read.reference_start
                    read_end = read.reference_end
                    if read_start is None or read_end is None or read_end <= read_start:
                        continue

                    rec_idx = len(records)
                    chrom = read.reference_name
                    read_contigs.add(chrom)
                    records.append((rec_idx, read_id, chrom, read_start, read_end))
                    rec_indices_by_contig.setdefault(chrom, []).append(rec_idx)
                    reads_out.write(f"{chrom}\t{read_start}\t{read_end}\t{rec_idx}\n")
                    n_mapped_records += 1

                    if n_total_records % 100000 == 0:
                        log(f"\t\t{n_total_records} reads processed for count-only report...", "I")

        # Split each bedGraph by contig with awk, restricted to contigs present in reads.
        split_by_track: List[Dict[str, Path]] = [{} for _ in tmp_bg_files]
        sorted_read_contigs = sorted(read_contigs)

        def _split_one_track_by_contig(track_idx: int, bg_path: Path) -> Tuple[int, Dict[str, Path], List[Path]]:
            created_paths: List[Path] = []
            per_contig_paths: Dict[str, Path] = {}

            # Keep only contigs that are present both in reads and in this track.
            awk_unique_contigs = 'BEGIN{FS="\\t"} NF>=4{seen[$1]=1} END{for(c in seen) print c}'
            log(f"awk '{awk_unique_contigs}' {bg_path}", "C")
            contigs_output = subprocess.check_output(
                ["awk", awk_unique_contigs, str(bg_path)],
                text=True,
            )
            track_contigs = {line.strip() for line in contigs_output.splitlines() if line.strip()}
            target_contigs = sorted(track_contigs.intersection(read_contigs))

            with tempfile.NamedTemporaryFile(prefix=f"wizardeye_t{track_idx}_contigs_", suffix=".tsv", delete=False) as tmp_map:
                contig_map_path = Path(tmp_map.name)
            created_paths.append(contig_map_path)

            with contig_map_path.open("w", encoding="utf-8") as map_out:
                for contig in target_contigs:
                    with tempfile.NamedTemporaryFile(
                        prefix=f"wizardeye_t{track_idx}_",
                        suffix=".bg",
                        delete=False,
                    ) as tmp_split:
                        split_path = Path(tmp_split.name)
                    created_paths.append(split_path)
                    per_contig_paths[contig] = split_path
                    map_out.write(f"{contig}\t{split_path}\n")

            if target_contigs:
                awk_script = (
                    "BEGIN{FS=OFS=\"\\t\"} "
                    "NR==FNR{map[$1]=$2; next} "
                    "($1 in map){print >> map[$1]; close(map[$1])}"
                )
                log(
                    f"awk '{awk_script}' {contig_map_path} {bg_path}",
                    "C",
                )
                subprocess.run(
                    ["awk", awk_script, str(contig_map_path), str(bg_path)],
                    check=True,
                )

            return track_idx, per_contig_paths, created_paths

        if n_threads > 1 and len(tmp_bg_files) > 1:
            max_workers = min(n_threads, len(tmp_bg_files))
            log(
                f"Splitting tracks by contig in parallel ({len(tmp_bg_files)} tracks, {max_workers} workers)...",
                "I",
            )
            with ThreadPoolExecutor(max_workers=max_workers) as pool:
                futures = [pool.submit(_split_one_track_by_contig, idx, bg) for idx, bg in enumerate(tmp_bg_files)]
                for future in as_completed(futures):
                    track_idx, per_contig_paths, created_paths = future.result()
                    split_by_track[track_idx] = per_contig_paths
                    tmp_extra_files.extend(created_paths)
        else:
            for idx, bg in enumerate(tmp_bg_files):
                track_idx, per_contig_paths, created_paths = _split_one_track_by_contig(idx, bg)
                split_by_track[track_idx] = per_contig_paths
                tmp_extra_files.extend(created_paths)

        # Keep only contigs that are present in reads and have data in at least one track.
        all_contigs: Set[str] = set()
        for contig in sorted_read_contigs:
            for per_contig_paths in split_by_track:
                split_path = per_contig_paths.get(contig)
                if split_path is not None and split_path.exists() and split_path.stat().st_size > 0:
                    all_contigs.add(contig)
                    break

        if all_contigs:
            with tempfile.NamedTemporaryFile(prefix="wizardeye_empty_", suffix=".bg", delete=False) as tmp_empty:
                empty_bg = Path(tmp_empty.name)
            tmp_extra_files.append(empty_bg)

            contig_union_outputs: Dict[str, Path] = {}

            def _run_unionbedg_for_contig(contig: str) -> Tuple[str, Path]:
                input_paths: List[Path] = []
                for per_contig_paths in split_by_track:
                    input_paths.append(per_contig_paths.get(contig, empty_bg))

                with tempfile.NamedTemporaryFile(
                    prefix=f"wizardeye_union_{contig}_",
                    suffix=".bg",
                    delete=False,
                ) as tmp_contig_union:
                    contig_out = Path(tmp_contig_union.name)
                tmp_extra_files.append(contig_out)

                log(
                    f"{bedtools} unionbedg -i {' '.join(str(path) for path in input_paths)} > {contig_out}",
                    "C",
                )
                with contig_out.open("w", encoding="utf-8") as contig_handle:
                    subprocess.run(
                        [bedtools, "unionbedg", "-i", *[str(path) for path in input_paths]],
                        check=True,
                        stdout=contig_handle,
                    )
                return contig, contig_out

            sorted_contigs = sorted(all_contigs)
            log(
                f"Count-only contig set reduced to {len(sorted_contigs)} contigs present in both reads and tracks.",
                "I",
            )
            if n_threads > 1 and len(sorted_contigs) > 1:
                max_workers = min(n_threads, len(sorted_contigs))
                log(
                    f"Running unionbedg in parallel by contig ({len(sorted_contigs)} contigs, {max_workers} workers)...",
                    "I",
                )
                with ThreadPoolExecutor(max_workers=max_workers) as pool:
                    futures = [pool.submit(_run_unionbedg_for_contig, contig) for contig in sorted_contigs]
                    for future in as_completed(futures):
                        contig, contig_out = future.result()
                        contig_union_outputs[contig] = contig_out
            else:
                for contig in sorted_contigs:
                    contig, contig_out = _run_unionbedg_for_contig(contig)
                    contig_union_outputs[contig] = contig_out

        else:
            sorted_contigs = []

        records_by_idx: Dict[int, Tuple[str, str, int, int]] = {
            rec_idx: (read_id, chrom, read_start, read_end)
            for rec_idx, read_id, chrom, read_start, read_end in records
        }

        # Intersect reads with union bedgraph per contig in parallel to save memory
        sorted_report_contigs = sorted_contigs
        contig_intersect_files: Dict[str, Path] = {}

        def _run_intersect_for_contig(contig: str) -> Tuple[str, Path]:
            """Run bedtools intersect for a specific contig between its reads and union bedgraph."""
            # Create temporary BED file for reads in this contig only
            with tempfile.NamedTemporaryFile(
                prefix=f"wizardeye_reads_{contig}_",
                suffix=".bed",
                delete=False,
            ) as tmp_reads_contig:
                reads_contig_path = Path(tmp_reads_contig.name)
            tmp_extra_files.append(reads_contig_path)

            # Write reads for this contig
            with reads_contig_path.open("w", encoding="utf-8") as reads_contig_out:
                for rec_idx in rec_indices_by_contig.get(contig, []):
                    _, read_id, chrom, read_start, read_end = records[rec_idx]
                    reads_contig_out.write(f"{chrom}\t{read_start}\t{read_end}\t{rec_idx}\n")

            # Get union bedgraph for this contig
            union_contig = contig_union_outputs.get(contig)
            if union_contig is None or not union_contig.exists():
                # Empty intersection file for this contig
                with tempfile.NamedTemporaryFile(
                    prefix=f"wizardeye_intersect_{contig}_",
                    suffix=".tsv",
                    delete=False,
                ) as tmp_contig_intersect:
                    contig_intersect = Path(tmp_contig_intersect.name)
                tmp_extra_files.append(contig_intersect)
                return contig, contig_intersect

            # Run bedtools intersect for this contig
            with tempfile.NamedTemporaryFile(
                prefix=f"wizardeye_intersect_{contig}_",
                suffix=".tsv",
                delete=False,
            ) as tmp_contig_intersect:
                contig_intersect = Path(tmp_contig_intersect.name)
            tmp_extra_files.append(contig_intersect)

            log(
                f"{bedtools} intersect -a {reads_contig_path} -b {union_contig} -wa -wb > {contig_intersect}",
                "C",
            )
            with contig_intersect.open("w", encoding="utf-8") as intersect_out:
                subprocess.run(
                    [bedtools, "intersect", "-a", str(reads_contig_path), "-b", str(union_contig), "-wa", "-wb"],
                    check=True,
                    stdout=intersect_out,
                )

            return contig, contig_intersect

        # Run intersects in parallel by contig
        if n_threads > 1 and len(sorted_report_contigs) > 1:
            max_workers = min(n_threads, len(sorted_report_contigs))
            log(
                f"Running bedtools intersect in parallel by contig ({len(sorted_report_contigs)} contigs, {max_workers} workers)...",
                "I",
            )
            with ThreadPoolExecutor(max_workers=max_workers) as pool:
                futures = [pool.submit(_run_intersect_for_contig, contig) for contig in sorted_report_contigs]
                for future in as_completed(futures):
                    contig, contig_intersect = future.result()
                    contig_intersect_files[contig] = contig_intersect
        else:
            for contig in sorted_report_contigs:
                contig, contig_intersect = _run_intersect_for_contig(contig)
                contig_intersect_files[contig] = contig_intersect

        contig_output_files: Dict[str, Path] = {}

        def _build_contig_table(contig: str) -> Tuple[str, Path]:
            intersect_contig = contig_intersect_files.get(contig)
            rec_indices = rec_indices_by_contig.get(contig, [])
            counts_by_rec: Dict[int, List[float]] = {}

            if intersect_contig is not None and intersect_contig.exists() and intersect_contig.stat().st_size > 0:
                with intersect_contig.open("r", encoding="utf-8") as handle:
                    for raw_line in handle:
                        line = raw_line.strip()
                        if not line:
                            continue
                        parts = line.split("\t")
                        if len(parts) < 7 + len(track_names):
                            continue
                        try:
                            a_start = int(parts[1])
                            a_end = int(parts[2])
                            rec_idx = int(parts[3])
                            b_start = int(parts[5])
                            b_end = int(parts[6])
                        except ValueError:
                            continue

                        overlap_start = max(a_start, b_start)
                        overlap_end = min(a_end, b_end)
                        overlap_len = overlap_end - overlap_start
                        if overlap_len <= 0:
                            continue

                        if rec_idx not in counts_by_rec:
                            counts_by_rec[rec_idx] = [0.0 for _ in track_names]
                        counts_row = counts_by_rec[rec_idx]
                        for col_idx in range(len(track_names)):
                            try:
                                depth = float(parts[7 + col_idx])
                            except ValueError:
                                continue
                            if depth > 0:
                                counts_row[col_idx] += overlap_len * depth

            with tempfile.NamedTemporaryFile(prefix=f"wizardeye_report_{contig}_", suffix=".tsv", delete=False) as tmp_contig_out:
                contig_out = Path(tmp_contig_out.name)
            tmp_extra_files.append(contig_out)

            with contig_out.open("w", encoding="utf-8") as out:
                for rec_idx in rec_indices:
                    rec_meta = records_by_idx.get(rec_idx)
                    if rec_meta is None:
                        continue
                    read_id, chrom, read_start, read_end = rec_meta
                    values = [read_id, chrom, str(read_start), str(read_end)]
                    counts_row = counts_by_rec.get(rec_idx)
                    if counts_row is None:
                        counts_row = [0.0 for _ in track_names]
                    values.extend(_format_count_value(value) for value in counts_row)
                    out.write("\t".join(values) + "\n")

            return contig, contig_out

        if n_threads > 1 and len(sorted_report_contigs) > 1:
            max_workers = min(n_threads, len(sorted_report_contigs))
            log(
                f"Building per-contig count tables in parallel ({len(sorted_report_contigs)} contigs, {max_workers} workers)...",
                "I",
            )
            with ThreadPoolExecutor(max_workers=max_workers) as pool:
                futures = [pool.submit(_build_contig_table, contig) for contig in sorted_report_contigs]
                for future in as_completed(futures):
                    contig, contig_out = future.result()
                    contig_output_files[contig] = contig_out
        else:
            for contig in sorted_report_contigs:
                contig, contig_out = _build_contig_table(contig)
                contig_output_files[contig] = contig_out

        # Merge final report as requested: header + simple cat of per-contig files.
        with tempfile.NamedTemporaryFile(prefix="wizardeye_report_header_", suffix=".tsv", delete=False) as tmp_header:
            header_path = Path(tmp_header.name)
        tmp_extra_files.append(header_path)
        header_line = "\t".join(["read_id", "chr", "start", "stop"] + track_names) + "\n"
        header_path.write_text(header_line, encoding="utf-8")

        ordered_contig_outputs = [contig_output_files[contig] for contig in sorted_report_contigs if contig in contig_output_files]
        cat_inputs = [header_path] + ordered_contig_outputs
        if cat_inputs:
            log(
                f"cat {' '.join(str(path) for path in cat_inputs)} > {output_report_tsv}",
                "C",
            )
            with output_report_tsv.open("w", encoding="utf-8") as final_out:
                subprocess.run(
                    ["cat", *[str(path) for path in cat_inputs]],
                    check=True,
                    stdout=final_out,
                )

        return {
            "n_total_records": n_total_records,
            "n_mapped_records": n_mapped_records,
            "n_rows": n_mapped_records,
        }
    finally:
        for extra_tmp in (reads_bed_path,):
            if extra_tmp is None:
                continue
            try:
                extra_tmp.unlink(missing_ok=True)
            except TypeError:
                if extra_tmp.exists():
                    extra_tmp.unlink()
        for path in tmp_extra_files:
            try:
                path.unlink(missing_ok=True)
            except TypeError:
                if path.exists():
                    path.unlink()
        for path in tmp_bg_files:
            try:
                path.unlink(missing_ok=True)
            except TypeError:
                if path.exists():
                    path.unlink()

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
    count_only: bool = False,
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

    if count_only:
        report_tsv = Path(output_report_tsv) if output_report_tsv else _default_output_count_table(input_bam_path)
    else:
        report_tsv = Path(output_report_tsv) if output_report_tsv else _default_output_table(input_bam_path)
    report_tsv.parent.mkdir(parents=True, exist_ok=True)

    # Count-only mode: produce one report from selected BigWig source, no mask/exclusion split.
    if count_only:
        selected_tracks = [
            track
            for track in tracks_for_params
            if track.track_name in normalized_requested_tracks
        ]
        if not selected_tracks:
            raise ValueError("No selected tracks available to compute count-only report")

        stats = _generate_count_only_report(
            input_bam=input_bam_path,
            selected_tracks=selected_tracks,
            consider_all=consider_all,
            output_report_tsv=report_tsv,
            n_threads=n_threads,
        )
        return {
            "mask": None,
            "filtered_bam": None,
            "excluded_bam": None,
            "report_tsv": report_tsv,
            "n_total": stats["n_mapped_records"],
            "n_filtered": stats["n_mapped_records"],
            "n_excluded": 0,
            "n_total_records": stats["n_total_records"],
            "count_only": True,
        }

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

