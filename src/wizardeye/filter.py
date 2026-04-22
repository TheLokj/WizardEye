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
from .mappability import get_unique_mapping_intervals, pysam
from .utils import log, from_charlist_to_list, write_bed

def _append_merged_interval(
    merged: List[Tuple[str, int, int]],
    chrom: str,
    start: int,
    end: int,
) -> None:
    """Append interval and merge with previous one when overlapping/adjacent."""
    if end <= start:
        return

    if not merged:
        merged.append((chrom, start, end))
        return

    last_chrom, last_start, last_end = merged[-1]
    if last_chrom == chrom and start <= last_end:
        merged[-1] = (last_chrom, last_start, max(last_end, end))
        return

    merged.append((chrom, start, end))


def _compute_overlap_intervals_from_bam(
    bam_file: Path,
    kmer_length: int,
    offset_step: int,
    cross_stringency: float,
) -> List[Tuple[str, int, int]]:
    """Compute stringency-filtered overlap intervals directly from BAM unique mappings."""
    if pysam is None:
        log("pysam is required to compute overlaps from BAM. Install it with: pip install pysam", "FE")

    if cross_stringency == 0.0:
        intervals: List[Tuple[str, int, int]] = []
        with pysam.AlignmentFile(str(bam_file), "rb") as bam:
            for chrom, length in sorted(zip(bam.references, bam.lengths), key=lambda x: x[0]):
                if length and length > 0:
                    intervals.append((chrom, 0, int(length)))
        return intervals

    events_by_chrom: Dict[str, Dict[int, int]] = {}
    for chrom, start, end in get_unique_mapping_intervals(bam_file, kmer_length):
        if end <= start:
            continue
        chrom_events = events_by_chrom.setdefault(chrom, {})
        chrom_events[start] = chrom_events.get(start, 0) + 1
        chrom_events[end] = chrom_events.get(end, 0) - 1

    threshold = cross_stringency * (kmer_length / offset_step)
    merged: List[Tuple[str, int, int]] = []

    for chrom in sorted(events_by_chrom):
        events = events_by_chrom[chrom]
        positions = sorted(events)
        depth = 0
        for idx in range(len(positions) - 1):
            pos = positions[idx]
            next_pos = positions[idx + 1]
            depth += events[pos]
            if next_pos <= pos:
                continue
            if depth >= threshold:
                _append_merged_interval(merged, chrom, pos, next_pos)

    return merged


def _compute_overlap_intervals_from_bigwig(
    uniq_bw: Path,
    kmer_length: int,
    offset_step: int,
    cross_stringency: float,
) -> List[Tuple[str, int, int]]:
    
    """Compute stringency-filtered overlap intervals from mappability_uniq.bw."""
    tool = shutil.which("bigWigToBedGraph")
    if tool is None:
        log("bigWigToBedGraph not found in PATH and BAM file is unavailable: cannot compute overlaps.", "FE")

    with tempfile.NamedTemporaryFile(prefix="wizardeye_", suffix=".bg", delete=False) as tmp_handle:
        tmp_bg_path = Path(tmp_handle.name)

    try:
        log(f"Converting track to bedGraph for overlap computation: {uniq_bw}", "I")
        log(f"{tool} {str(uniq_bw)} {str(tmp_bg_path)}", "C")
        subprocess.run([tool, str(uniq_bw), str(tmp_bg_path)], check=True)

        threshold = cross_stringency * (kmer_length / offset_step)
        merged: List[Tuple[str, int, int]] = []
        with tmp_bg_path.open("r", encoding="utf-8") as handle:
            for line in handle:
                line = line.strip()
                if not line:
                    continue
                parts = line.split("\t")
                if len(parts) < 4:
                    continue
                chrom, start, end, depth = parts[0], int(parts[1]), int(parts[2]), float(parts[3])
                if end <= start:
                    continue
                if depth >= threshold:
                    _append_merged_interval(merged, chrom, start, end)
        return merged
    finally:
        if tmp_bg_path.exists():
            tmp_bg_path.unlink()


def _compute_overlap_bed_from_bigwig(
    uniq_bw: Path,
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
        log(f"Converting track to bedGraph for overlap computation: {uniq_bw}", "I")
        log(f"{tool} {str(uniq_bw)} {str(tmp_bg_path)}", "C")
        subprocess.run([tool, str(uniq_bw), str(tmp_bg_path)], check=True)

        threshold = cross_stringency * (kmer_length / offset_step)
        bedtools = shutil.which("bedtools")
        awk = shutil.which("awk")
        if bedtools and awk:
            awk_script = f'OFS="\\t" {{ if (($4 + 0) >= {threshold:.12f}) print $1, $2, $3; }}'
            with tmp_bg_path.open("r", encoding="utf-8") as bg_in, output_bed.open("w", encoding="utf-8") as out:
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

        # Fallback: Python thresholding/merging when bedtools/awk are unavailable.
        intervals = _compute_overlap_intervals_from_bigwig(
            uniq_bw=uniq_bw,
            kmer_length=kmer_length,
            offset_step=offset_step,
            cross_stringency=cross_stringency,
        )
        write_bed(intervals, output_bed)
        return output_bed
    finally:
        if tmp_bg_path.exists():
            tmp_bg_path.unlink()


def _read_bed_intervals(bed_path: Path) -> List[Tuple[str, int, int]]:
    intervals: List[Tuple[str, int, int]] = []
    with bed_path.open("r", encoding="utf-8") as handle:
        for line in handle:
            parts = line.strip().split("\t")
            if len(parts) < 3:
                continue
            start = int(parts[1])
            end = int(parts[2])
            if end > start:
                intervals.append((parts[0], start, end))
    return intervals


def _merge_bed_files(bed_files: List[Tuple[str, Path]], output_bed: Path) -> Path:
    output_bed.parent.mkdir(parents=True, exist_ok=True)
    if not bed_files:
        output_bed.write_text("", encoding="utf-8")
        return output_bed

    bedtools = shutil.which("bedtools")
    if bedtools:
        with tempfile.NamedTemporaryFile(prefix="wizardeye_", suffix=".bed", delete=False) as tmp_handle:
            tmp_tagged_path = Path(tmp_handle.name)

        try:
            with tmp_tagged_path.open("w", encoding="utf-8") as tagged_out:
                for track_name, bed_path in bed_files:
                    with bed_path.open("r", encoding="utf-8") as handle:
                        for line in handle:
                            parts = line.strip().split("\t")
                            if len(parts) < 3:
                                continue
                            start = int(parts[1])
                            end = int(parts[2])
                            if end <= start:
                                continue
                            tagged_out.write(f"{parts[0]}\t{start}\t{end}\t{track_name.split('_k')[0]}\n")

            with tmp_tagged_path.open("r", encoding="utf-8") as tagged_in:
                p_sort = subprocess.Popen(
                    ["sort", "-k1,1", "-k2,2n"],
                    stdin=tagged_in,
                    stdout=subprocess.PIPE,
                    text=True,
                )
                with output_bed.open("w", encoding="utf-8") as out:
                    p_merge = subprocess.Popen(
                        [bedtools, "merge", "-i", "stdin", "-c", "4", "-o", "distinct"],
                        stdin=p_sort.stdout,
                        stdout=out,
                        text=True,
                    )

                    if p_sort.stdout is not None:
                        p_sort.stdout.close()

                    rc_merge = p_merge.wait()
                    rc_sort = p_sort.wait()
                    if rc_sort != 0:
                        raise subprocess.CalledProcessError(rc_sort, ["sort", "-k1,1", "-k2,2n"])
                    if rc_merge != 0:
                        raise subprocess.CalledProcessError(
                            rc_merge,
                            [bedtools, "merge", "-i", "stdin", "-c", "4", "-o", "distinct"],
                        )
            return output_bed
        finally:
            if tmp_tagged_path.exists():
                tmp_tagged_path.unlink()

    all_intervals: List[Tuple[str, int, int, Set[str]]] = []
    for track_name, bed in bed_files:
        for chrom, start, end in _read_bed_intervals(bed):
            all_intervals.append((chrom, start, end, {track_name}))

    all_intervals.sort(key=lambda x: (x[0], x[1], x[2]))
    merged_all: List[Tuple[str, int, int, Set[str]]] = []
    for chrom, start, end, tracks in all_intervals:
        if not merged_all:
            merged_all.append((chrom, start, end, set(tracks)))
            continue

        last_chrom, last_start, last_end, last_tracks = merged_all[-1]
        if last_chrom == chrom and start <= last_end:
            merged_all[-1] = (last_chrom, last_start, max(last_end, end), last_tracks | tracks)
        else:
            merged_all.append((chrom, start, end, set(tracks)))

    with output_bed.open("w", encoding="utf-8") as out:
        for chrom, start, end, tracks in merged_all:
            joined_tracks = ",".join(sorted(tracks))
            out.write(f"{chrom}\t{start}\t{end}\t{joined_tracks}\n")
    return output_bed


def _build_track_overlaps(
    track: Any,
    ref_dir: Path,
    stringency_label: str,
    kmer_length: int,
    offset_step: int,
    cross_stringency: float,
) -> Tuple[str, Path]:
    """Compute one track overlap BED for a selected input and return (track name, BED path)."""
    track_dir = track.track_dir
    track_name = track.track_name

    # Use canonical track_name to avoid collisions between tracks sharing the same query/input name.
    per_input_bed = ref_dir / f"overlap_{track_name}_s{stringency_label}.bed"
    if per_input_bed.exists():
        log(f"Reusing existing overlap BED for track '{track_name}': {per_input_bed}", "I")
        return track_name, per_input_bed

    uniq_bw = track_dir / "mappability_uniq.bw"
    if uniq_bw.exists():
        _compute_overlap_bed_from_bigwig(
            uniq_bw=uniq_bw,
            kmer_length=kmer_length,
            offset_step=offset_step,
            cross_stringency=cross_stringency,
            output_bed=per_input_bed,
        )
    else:
        bam_candidates = list(track_dir.glob("*.bam"))
        if not bam_candidates:
            raise ValueError(
                f"No mappability_uniq.bw or BAM found for selected track '{track_name}' in {track_dir}"
            )
        input_intervals = _compute_overlap_intervals_from_bam(
            bam_file=bam_candidates[0],
            kmer_length=kmer_length,
            offset_step=offset_step,
            cross_stringency=cross_stringency,
        )
        write_bed(input_intervals, per_input_bed)
        if not input_intervals:
            threshold = cross_stringency * (kmer_length / offset_step)
            log(
                f"No interval passed stringency for track '{track_name}' (threshold unique-depth >= {threshold:.4f}).",
                "W",
            )

    log(f"Overlap BED written for track '{track_name}': {per_input_bed}", "S")
    return track_name, per_input_bed

def generate_mask(
    ref_species: str,
    inputs: Optional[List[str]],
    kmer_length: int,
    offset_step: int,
    cross_stringency: float,
    n_threads: int = 1,
    db_root: str = "database",
    output_file: Optional[str] = None,
    bwa_missing_prob_err_rate: Optional[float] = None,
    bwa_max_gap_opens: Optional[int] = None,
    bwa_seed_length: Optional[int] = None,
) -> Path:
    """
    Generate overlap BED files for selected inputs and one merged all-overlaps BED.

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
                    _build_track_overlaps,
                    track,
                    ref_dir,
                    stringency_label,
                    kmer_length,
                    offset_step,
                    cross_stringency,
                )
                for track in selected_tracks
            ]
            for future in as_completed(futures):
                per_track_beds.append(future.result())
    else:
        for track in selected_tracks:
            per_track_beds.append(_build_track_overlaps(
                track=track,
                ref_dir=ref_dir,
                stringency_label=stringency_label,
                kmer_length=kmer_length,
                offset_step=offset_step,
                cross_stringency=cross_stringency,
            ))

    if output_file:
        all_overlaps_bed = Path(output_file)
    else:
        all_overlaps_bed = ref_dir / f"all_overlaps_mask_{stringency_label}_k{kmer_length}_o{offset_step}.bed"

    _merge_bed_files(per_track_beds, all_overlaps_bed)
    log(f"Merged all-overlaps BED written to: {all_overlaps_bed}", "S")
    return all_overlaps_bed


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
			_append_merged_interval(merged, chrom, start, end)

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


def _normalize_contig_name(contig: str) -> str:
    """Normalize common naming variants (chr prefix and M/MT)."""
    normalized = contig.strip()
    if normalized.startswith("chr"):
        normalized = normalized[3:]
    if normalized == "M":
        return "MT"
    return normalized


def _read_mask_contigs(mask_path: Path) -> Set[str]:
    """Read chromosome names from a BED mask."""
    contigs: Set[str] = set()
    with mask_path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if not parts:
                continue
            chrom = parts[0].strip()
            if chrom:
                contigs.add(chrom)
    return contigs


def _load_mask_with_tracks(mask_path: Path) -> Dict[str, List[Tuple[int, int, Set[str]]]]:
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


def _collect_overlapped_tracks(
    chrom: str,
    read_start: int,
    read_end: int,
    intervals_by_chrom: Dict[str, List[Tuple[int, int, Set[str]]]],
) -> Set[str]:
    """Return all track labels overlapping one read interval on one contig."""
    if read_end <= read_start:
        return set()

    intervals = intervals_by_chrom.get(chrom)
    if not intervals:
        return set()

    starts = [interval[0] for interval in intervals]
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


def _write_read_exclusion_report(
    input_bam: Path,
    mask_path: Path,
    output_report_tsv: Path,
    track_to_tags: Optional[Dict[str, Set[str]]] = None,
) -> Path:
    """Write one line per read_id with exclusion flag, overlapping tracks and tags."""
    if pysam is None:
        raise RuntimeError("pysam is required to generate read-level overlap reports")

    track_to_tags = track_to_tags or {}
    intervals_by_chrom = _load_mask_with_tracks(mask_path)
    read_tracks: Dict[str, Set[str]] = {}
    read_tags: Dict[str, Set[str]] = {}
    read_seen_order: List[str] = []

    with pysam.AlignmentFile(str(input_bam), "rb") as bam:
        for read in bam.fetch(until_eof=True):
            read_id = read.query_name
            if not read_id:
                continue

            if read_id not in read_tracks:
                read_tracks[read_id] = set()
                read_tags[read_id] = set()
                read_seen_order.append(read_id)

            if read.is_unmapped or read.reference_name is None:
                continue

            read_start = read.reference_start
            read_end = read.reference_end
            if read_start is None or read_end is None:
                continue

            overlapped_tracks = _collect_overlapped_tracks(
                chrom=read.reference_name,
                read_start=read_start,
                read_end=read_end,
                intervals_by_chrom=intervals_by_chrom,
            )
            read_tracks[read_id].update(overlapped_tracks)
            for track_name in overlapped_tracks:
                read_tags[read_id].update(track_to_tags.get(track_name, set()))

    output_report_tsv.parent.mkdir(parents=True, exist_ok=True)
    with output_report_tsv.open("w", encoding="utf-8") as handle:
        handle.write("read_id\texcluded\toverlapped\ttags\n")
        for read_id in read_seen_order:
            tracks = sorted(read_tracks.get(read_id, set()))
            tags = sorted(read_tags.get(read_id, set()))
            excluded = "true" if tracks else "false"
            overlapped = ",".join(tracks)
            joined_tags = ",".join(tags)
            handle.write(f"{read_id}\t{excluded}\t{overlapped}\t{joined_tags}\n")

    return output_report_tsv


def _warn_if_contig_names_mismatch(mask_path: Path, bam_path: Path) -> None:
    """Warn when BED/BAM contig naming conventions are likely incompatible."""
    if pysam is None:
        return

    try:
        with pysam.AlignmentFile(str(bam_path), "rb") as bam:
            bam_contigs = set(bam.references)
    except Exception:
        return

    try:
        mask_contigs = _read_mask_contigs(mask_path)
    except Exception:
        return

    if not bam_contigs or not mask_contigs:
        return

    exact_overlap = bam_contigs & mask_contigs
    if exact_overlap:
        return

    bam_norm = {_normalize_contig_name(contig) for contig in bam_contigs}
    mask_norm = {_normalize_contig_name(contig) for contig in mask_contigs}
    if bam_norm & mask_norm:
        sample_bam = ", ".join(sorted(list(bam_contigs))[:5])
        sample_mask = ", ".join(sorted(list(mask_contigs))[:5])
        log(
            "No exact chromosome-name overlap between BAM and mask, but overlap appears after normalization "
            "(likely 'chr' prefix mismatch). Filtering may return excluded=0 unless contig names are harmonized. "
            f"BAM sample: {sample_bam} | Mask sample: {sample_mask}",
            "W",
        )


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
    n_threads: int = 1,
    output_filtered_bam: Optional[str] = None,
    output_excluded_bam: Optional[str] = None,
    output_report_tsv: Optional[str] = None,
) -> Dict[str, object]:
    """Build a merged mask for selected tracks and filter one input BAM against it.

    Output BAM files:
    - filtered: reads that do NOT overlap the generated mask
    - excluded: reads that DO overlap the generated mask
    """
    if shutil.which("samtools") is None:
        raise RuntimeError("samtools not found in PATH, cannot filter BAM")

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

    # Keep the merged mask deterministic for one exact set of tracks and parameters.
    sorted_tracks = sorted(set(normalized_tracks))
    track_fingerprint = hashlib.sha1(",".join(sorted_tracks).encode("utf-8")).hexdigest()[:12]
    ref_dir = Path(db_root) / ref
    mask_path = ref_dir / (
        f"all_overlaps_mask_s{stringency}_k{kmer_length}_o{offset_step}"
        f"_t{len(sorted_tracks)}_{track_fingerprint}.bed"
    )

    merged_mask = generate_mask(
        ref_species=ref,
        inputs=sorted_tracks,
        kmer_length=kmer_length,
        offset_step=offset_step,
        cross_stringency=stringency,
        n_threads=n_threads,
        db_root=db_root,
        output_file=str(mask_path),
        bwa_missing_prob_err_rate=bwa_missing_prob_err_rate,
        bwa_max_gap_opens=bwa_max_gap_opens,
        bwa_seed_length=bwa_seed_length,
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

    _warn_if_contig_names_mismatch(merged_mask, input_bam_path)

    filtered_bam = Path(output_filtered_bam) if output_filtered_bam else _default_output_bam(input_bam_path, "filtered")
    excluded_bam = Path(output_excluded_bam) if output_excluded_bam else _default_output_bam(input_bam_path, "excluded")
    report_tsv = Path(output_report_tsv) if output_report_tsv else _default_output_table(input_bam_path)
    filtered_bam.parent.mkdir(parents=True, exist_ok=True)
    excluded_bam.parent.mkdir(parents=True, exist_ok=True)
    report_tsv.parent.mkdir(parents=True, exist_ok=True)

    cmd = [
        "samtools",
        "view",
        "-h",
        "-b",
        "-L",
        str(merged_mask),
        "-U",
        str(filtered_bam),
        "-o",
        str(excluded_bam),
        str(input_bam_path),
    ]
    log(" ".join(cmd), "C")
    subprocess.run(cmd, check=True)

    log(f"samtools view -c {str(input_bam_path)}", "C")
    n_total = int(subprocess.check_output(["samtools", "view", "-c", str(input_bam_path)], text=True).strip())
    log(f"samtools view -c {str(filtered_bam)}", "C")
    n_filtered = int(subprocess.check_output(["samtools", "view", "-c", str(filtered_bam)], text=True).strip())
    log(f"samtools view -c {str(excluded_bam)}", "C")
    n_excluded = int(subprocess.check_output(["samtools", "view", "-c", str(excluded_bam)], text=True).strip())

    for bam_path in (filtered_bam, excluded_bam):
        try:
            log(f"samtools index {str(bam_path)}", "C")
            subprocess.run(["samtools", "index", str(bam_path)], check=True)
        except subprocess.CalledProcessError:
            log(f"Could not index BAM (possibly unsorted): {bam_path}", "W")

    report_path = _write_read_exclusion_report(
        input_bam=input_bam_path,
        mask_path=merged_mask,
        output_report_tsv=report_tsv,
        track_to_tags=track_to_tags,
    )
    log(f"Read exclusion report written: {report_path}", "S")

    return {
        "mask": merged_mask,
        "filtered_bam": filtered_bam,
        "excluded_bam": excluded_bam,
        "report_tsv": report_path,
        "n_total": n_total,
        "n_filtered": n_filtered,
        "n_excluded": n_excluded,
    }

