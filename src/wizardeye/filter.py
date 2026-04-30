# -*- coding: utf-8 -*-

"""Filter alignment BAM utilities for WizardEye.

This Python module provides functions to generate mask using previously generated
cross-mappability tracks and to filter BAM with such masks. It also provides functions 
to apply stringency-threshold, to compute raw overlapping-k-mers count per read and 
to generate the final filtration report.
"""

from __future__ import annotations
from collections import defaultdict
import hashlib 
import pyBigWig
import shutil
import subprocess
import tempfile

from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

from .db import get_tracks, Track
import pysam
from .utils import (
	run,
	log,
	from_charlist_to_list,
	validate_bam_compatibility,
	merge_bed_files,
	convert_bigwig_to_bedGraph
)

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
		with tempfile.NamedTemporaryFile(prefix=f"wizardeye_{track_name}_", suffix=".bed", delete=False) as tmp_handle:
			per_input_bed = Path(tmp_handle.name)
	else:
		per_input_bed = track_dir / f"mask_{track_name}_s{cross_stringency}_{source_label}.bed"

	if not no_cache and per_input_bed.exists():
		log(f"Reusing existing mask for {track_name}... \n\033[0;90m({per_input_bed})\033[0m", "I")
		return track_name, per_input_bed

	source_bw_name = "map_all.bw" if consider_all else "map_uniq.bw"
	source_bw = track_dir / source_bw_name
	if not source_bw.exists():
		raise ValueError(
			f"Missing required {source_bw_name} for selected track '{track_name}' in {track_dir}"
		)

	with tempfile.NamedTemporaryFile(prefix=f"wizardeye_{track_name}_", suffix=".bg", delete=False) as tmp_bg:
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
	return track_name, per_input_bed

def generate_global_mask(
	ref_species: str,
	inputs: List[str],
	kmer_length: int,
	offset_step: int,
	cross_stringency: float,
	consider_all: bool = False,
	output_file: Optional[str] = None,
	bwa_missing_prob_err_rate: Optional[float] = None,
	bwa_max_gap_opens: Optional[int] = None,
	bwa_seed_length: Optional[int] = None,
	db_root: str = "database",
	no_cache: bool = False,
	n_threads: int = 1,
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
		bwa_missing_prob_err_rate (Optional[float]): BWA missing probability error rate
			used during track generations.
		bwa_max_gap_opens (Optional[int]): Maximum number of gap opens allowed in BWA alignment 
			used during track generations.
		bwa_seed_length (Optional[int]): Seed length for BWA alignment 
			used during track generations.
		db_root (str): Root directory of the database.
		no_cache (bool): If True, do not use and generate cached tracks.
		n_threads (int): Number of threads to use for parallel processing.

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

	if cross_stringency > 1.0:
		log("Cross-stringency is superior to 1.0.", "W")

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
			f"mask_s{cross_stringency}_k{kmer_length}_o{offset_step}"
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
				cross_stringency=cross_stringency,
				consider_all=consider_all,
				no_cache=no_cache,
			))

	try:
		log(f"Merging {len(per_track_beds)} masks into one unique mask...", "I")
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
	
	with input_bg.open("r", encoding="utf-8") as bg_in, output_bed.open("w", encoding="utf-8") as out:
		log(f"{awk} '{awk_script}' {input_bg} | {bedtools} merge -i stdin", "C")
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


# -- Filtration related functions --

def filter_bam(
	input_bam: str,
	ref: str,
	db_root: str,
	exclude_tracks: List[Track],
	kmer_length: int,
	offset_step: int,
	bwa_missing_prob_err_rate: Optional[float] = None,
	bwa_max_gap_opens: Optional[int] = None,
	bwa_seed_length: Optional[int] = None,
	stringency: float = 0.99,
	consider_all: bool = False,
	output_report_tsv: Optional[str] = None,
	export_bam: bool = False,
	output_filtered_bam: Optional[str] = None,
	output_excluded_bam: Optional[str] = None,
	no_cache: bool = False,
	n_threads: int = 1,
) -> Dict[str, object]:
	"""Main function to filter a BAM file using several tracks and generate requested outputs.

	Args:
		input_bam (str): Path to the input BAM file to filter.
		ref (str): Name of the reference species.
		exclude_tracks (List[str]): List of input track names to consider. 
		kmer_length (int): Length of the k-mers used during track generations.
		offset_step (int): Step size for the offset used during track generations.
		cross_stringency (float): Cross-stringency threshold, must be between 0.0 and 1.0.
		consider_all (bool): If True, all k-mers are considered in the cross-stringency
            calculation. If False (default), only uniquely aligned k-mers are considered.
		bwa_missing_prob_err_rate (Optional[float]): BWA missing probability error rate
			used during track generations.
		bwa_max_gap_opens (Optional[int]): Maximum number of gap opens allowed in BWA alignment 
			used during track generations.
		bwa_seed_length (Optional[int]): Seed length for BWA alignment 
			used during track generations.
		output_report_tsv (Optional[str]): Path to the output TSV report file.
		export_bam: If True, splits the input BAM into filtered and excluded BAM files.
		output_filtered_bam (Optional[str]): Path to the output BAM file containing filtered reads if export_bam is True.
		output_excluded_bam (Optional[str]): Path to the output BAM file containing excluded reads if export_bam is True.
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
	if not (0.0 <= stringency <= 1.0):
		raise ValueError("stringency must be between 0.0 and 1.0")

	normalized_tracks = from_charlist_to_list(exclude_tracks)
	if not normalized_tracks:
		raise ValueError("exclude_tracks cannot be empty")

	# Fail fast on initial-file incompatibility before any mask computation.
	ref_dir = Path(db_root) / ref
	reference_seq_sizes = ref_dir / f"{ref}.sizes"
	validate_bam_compatibility(
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

	report_tsv = Path(output_report_tsv) if output_report_tsv else _default_output_table(input_bam_path)
	report_tsv.parent.mkdir(parents=True, exist_ok=True)

	merged_mask = generate_global_mask(
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
		# Clean up any temporary mask generated when cache is disabled.
		if no_cache:
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

def _filter_bam_with_mask(
    input_bam: Path,
    mask_path: Path,
    track_to_tags: Optional[Dict[str, Set[str]]] = None,
) -> Tuple[Dict[str, Set[str]], Dict[str, Set[str]], Set[str], int, int, int]:
    """Filter a BAM file using a mask file and identify overlapping reads.
	
	Args:
		input_bam (Path): Path to the input BAM file.
		mask_path (Path): Path to the mask file.
		track_to_tags (Optional[Dict[str, Set[str]]]): Dictionary mapping track names to sets of tags.

	Returns:
		- Dict[str, Set[str]]: read_tracks mapping read IDs to overlapping track names
		- Dict[str, Set[str]]: read_tags mapping read IDs to overlapping track tags
		- Set[str]: excluded_read_ids set of read IDs that overlap the mask
		- int: n_total total number of distinct read IDs
		- int: n_excluded number of excluded read IDs
		- int: n_total_records total number of input BAM records processed

	Notes:
		This function currently use subprocesses to avoid to parse every .bw in memory.
		This may evolve in future version to use pyBigWig instead. However, it currently
		excludes reads using IDs, which can lead to the exclusion of different reads using
		same IDs. 
	"""

    log("Identifying overlapping reads...", "I")
	
	# Initialize data structures
    track_to_tags = track_to_tags or {}
    read_tracks = defaultdict(set)
    read_tags = defaultdict(set)
    excluded_read_ids: Set[str] = set()
    
    all_read_ids = {}
    n_total_records = 0
    
	# Get reads count
    log(f"samtools view {input_bam}", "C")
    view_cmd = ["samtools", "view", str(input_bam)]
    view_proc = subprocess.Popen(view_cmd, stdout=subprocess.PIPE, text=True)
    for line in view_proc.stdout:
        n_total_records += 1
        read_id = line.split('\t', 1)[0] 
        all_read_ids[read_id] = None
    view_proc.wait()

    for rid in all_read_ids.keys():
        read_tracks[rid] = set()
        read_tags[rid] = set()

	# Get the intersection between the BAM file and the mask file
    log(f"bedtools intersect -abam {input_bam} -b {mask_path} -wa -wb -bed", "C")
    intersect_cmd = [
        "bedtools", "intersect",
        "-abam", str(input_bam),
        "-b", str(mask_path),
        "-wa", "-wb", "-bed"
    ]
    
    process = subprocess.Popen(
        intersect_cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True
    )

	# Process reads 
    try:
        for line in process.stdout:
            parts = line.strip().split("\t")
            if len(parts) >= 16:
                read_id = parts[3]  
                excluded_read_ids.add(read_id)

                tracks_str = parts[-1]
                tracks = tracks_str.split(",") if tracks_str else []

                for track in tracks:
                    if track:
                        read_tracks[read_id].add(track)
                        read_tags[read_id].update(track_to_tags.get(track, set()))
            else:
                raise RuntimeError(f"Unexpected ouput from bedtools: {line}")
    finally:
        process.stdout.close()

    if process.wait() != 0:
        raise subprocess.CalledProcessError(process.returncode, intersect_cmd)

    if len(all_read_ids) != n_total_records:
        log(f"Number of unique IDs in BAM ({len(all_read_ids)}) does not match number of records in BED ({n_total_records}): this could be"\
			" due to duplicates. Currently, WizardEye filter using read-id and will then exclude duplicates if at least one is filtered.", "W")

    return (
        dict(read_tracks),
        dict(read_tags),
        excluded_read_ids,
        len(all_read_ids),   
        len(excluded_read_ids), 
        n_total_records,    
    )

def filter_bam_from_reads_id(
	input_bam: Path,
	excluded_read_ids: Set[str],
	output_filtered_bam: Path,
	output_excluded_bam: Path,
) -> Tuple[int, int, int]:
	"""Split BAM in one pysam pass using excluded read IDs.
	
	Args:
		input_bam: Path to input BAM file.
		excluded_read_ids: Set of read IDs to exclude.
		output_filtered_bam: Path to output BAM file with filtered reads.
		output_excluded_bam: Path to output BAM file with excluded reads.
		
	Returns:
		Tuple[int, int, int]: Number of total records, number of filtered records and number of excluded records."""

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

# --- Count-report related function ---

def count_k_mers_on_bam(
	input_bam: str,
	ref: str,
	db_root: str,
	exclude_tracks: List[Track],
	kmer_length: int,
	offset_step: int,
	count_mode: str,
	bwa_missing_prob_err_rate: Optional[float] = None,
	bwa_max_gap_opens: Optional[int] = None,
	bwa_seed_length: Optional[int] = None,
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
		consider_all (bool): If True, all k-mers are considered in the cross-stringency
            calculation. If False (default), only uniquely aligned k-mers are considered.
		bwa_missing_prob_err_rate (Optional[float]): BWA missing probability error rate
			used during track generations.
		bwa_max_gap_opens (Optional[int]): Maximum number of gap opens allowed in BWA alignment 
			used during track generations.
		bwa_seed_length (Optional[int]): Seed length for BWA alignment 
			used during track generations.
		output_report_tsv (Optional[str]): Path to the output TSV report file.
		export_bam: If True, splits the input BAM into filtered and excluded BAM files.
		output_filtered_bam (Optional[str]): Path to the output BAM file containing filtered reads if export_bam is True.
		output_excluded_bam (Optional[str]): Path to the output BAM file containing excluded reads if export_bam is True.
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

	report_tsv = Path(output_report_tsv) if output_report_tsv else _default_output_count_table(input_bam_path)
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
    track_names = [f"{count_mode}_{track.identity.query_species}" for track in selected_tracks]
    
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
                        stat_res = bw.stats(chrom, start, end, type=count_mode, exact=True)
                        row_values.append(f"{stat_res[0]:.0f}" if count_mode in ['max', 'min', 'sum'] else f"{stat_res[0]:.3f}")

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
	output_report_tsv: Path,
	read_tracks: Dict[str, Set[str]],
	read_tags: Dict[str, Set[str]],
) -> Path:
	"""Write one line per read_id with exclusion flag, overlapping tracks and tags.
	
	Args:
		output_report_tsv (Path): Path to the output TSV file.
		read_tracks (Dict[str, Set[str]]): Dictionary mapping read IDs to sets of overlapping tracks.
		read_tags (Dict[str, Set[str]]): Dictionary mapping read IDs to sets of tags.
		
	Returns:
		Path: The path to the generated report."""
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