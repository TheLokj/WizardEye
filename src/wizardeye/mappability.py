# -*- coding: utf-8 -*-

"""Mappability track generation and export utilities for WizardEye.

This Python module provides functions to create mappability tracks from input FASTA files
by splitting a query sequence into k-mers and aligning them to a reference genome using BWA, 
and exporting the results as BigWig files in order to fill the database. 

This module is the direct Python implementation of the original generate_cross_mappability_filter_bwa.sh script, 
with added features and optimizations.

It includes utilities for handling temporary files, validating inputs, and saving track metadata.

	See Also:
		The original script: https://github.com/TheLokj/generate_cross_mappability_filter/blob/master/BWA/generate_cross_mappability_filter_bwa.sh
"""

import subprocess
import shutil
import yaml
import os
import getpass

from pathlib import Path
from datetime import datetime
from concurrent.futures import ThreadPoolExecutor
from typing import Dict, List, Optional, Tuple

from .version import PACKAGE_VERSION
from .utils import (
	log,
	run,
	file_md5,
	from_charlist_to_list,
	convert_bedgraph_to_bigwig,
	write_seq_sizes_from_bam,
	iterate_mapping_intervals,
	iterate_unique_mapping_intervals,
	Intervals,
	Interval,
)

# -- WizardEye mappability utilities --

def from_alignment_to_bigWig(bam_file: Path, out_dir: Path, kmer_length: int) -> Tuple[Path, Path, Dict[str, int]]:
	"""Convert alignment results to BigWig tracks.
	
	Args:
		bam_file(Path): Path to the BAM file containing all mapped k-mers.
		out_dir(Path): Directory where the output BigWig files and intermediate files will be saved.
		kmer_length(int): Length of the k-mers used for mapping, needed for coverage calculations.

	Returns:
		Tuple[Path, Path, Dict[str, int] :
			- Path: Path to the BigWig file for all mappings (map_all).
			- Path: Path to the BigWig file for uniquely mapping k-mers (map_uniq).
			- Dict[str, int]: dictionary with covered base counts for both masks, e.g. {"map_all": 123456, "map_uniq": 78910}.
	"""

	cov_map_bg = out_dir / "map_all.bg"
	cov_uniq_bg = out_dir / "map_uniq.bg"
	cov_map_bw = out_dir / "map_all.bw"
	cov_uniq_bw = out_dir / "map_uniq.bw"
	seq_sizes = out_dir / "chrom.sizes"

	intervals_all = Intervals()
	for chrom, start, end in iterate_mapping_intervals(bam_file, kmer_length):
		intervals_all.append(Interval(chrom, start, end))
	intervals_all.write_to_bedgraph(cov_map_bg)

	intervals_uniq = Intervals()
	for chrom, start, end in iterate_unique_mapping_intervals(bam_file, kmer_length):
		intervals_uniq.append(Interval(chrom, start, end))
	intervals_uniq.write_to_bedgraph(cov_uniq_bg)

	covered_bp = {
		"map_all": intervals_all.count_covered_bases(),
		"map_uniq": intervals_uniq.count_covered_bases(),
	}
	write_seq_sizes_from_bam(bam_file, seq_sizes)

	convert_bedgraph_to_bigwig(cov_map_bg, seq_sizes, cov_map_bw)
	convert_bedgraph_to_bigwig(cov_uniq_bg, seq_sizes, cov_uniq_bw)

	for tmp_file in (cov_map_bg, cov_uniq_bg, seq_sizes):
		if tmp_file.exists():
			tmp_file.unlink()

	return cov_map_bw, cov_uniq_bw, covered_bp

# --- WizardEye mappability pipeline blocks ---

def align_with_bwa_aln(
	input_fasta: Path,
	reference: Path,
	bwa_missing_prob_err_rate: float,
	bwa_max_gap_opens: int,
	bwa_seed_length: int,
) -> Path:
	"""Align a FASTA file using bwa aln and return the path to the resulting BAM file containing mapped reads.
	
	Args:
		input_fasta(Path): Path to the input FASTA file containing k-mers to align.
		reference(Path): Path to the reference genome FASTA file that has been indexed with BWA.
		bwa_missing_prob_err_rate(float): The -n parameter for bwa aln, controlling the maximum edit distance.
		bwa_max_gap_opens(int): The -o parameter for bwa aln, controlling the maximum number of gap opens.
		bwa_seed_length(int): The -l parameter for bwa aln, controlling the seed length.
		
	Returns:
		Path: the BAM file containing the mapped reads for the input chunk.
	"""
	sai_file = input_fasta.with_suffix(".sai")
	bam_file = input_fasta.with_suffix(".bam")

	with sai_file.open("w", encoding="utf-8") as sai_out:
		run(
			[
				"bwa",
				"aln",
				"-t",
				"1",
				"-n",
				str(bwa_missing_prob_err_rate),
				"-o",
				str(bwa_max_gap_opens),
				"-l",
				str(bwa_seed_length),
				str(reference),
				str(input_fasta),
			],
			check=True,
			stdout=sai_out,
		)

	with bam_file.open("wb") as bam_out:
		log(f"bwa samse -n 1000000 {reference} {sai_file} {input_fasta}", "C")
		bwa_samse = subprocess.Popen(
			[
				"bwa",
				"samse",
				"-n",
				"1000000",
				str(reference),
				str(sai_file),
				str(input_fasta),
			],
			stdout=subprocess.PIPE,
		)
		try:
			run(
				["samtools", "view", "-b", "-F", "4", "-"],
				check=True,
				stdin=bwa_samse.stdout,
				stdout=bam_out,
			)
		finally:
			if bwa_samse.stdout is not None:
				bwa_samse.stdout.close()

		return_code = bwa_samse.wait()
		if return_code != 0:
			raise subprocess.CalledProcessError(return_code, bwa_samse.args)

	for tmp_file in (sai_file, input_fasta):
		if tmp_file.exists():
			tmp_file.unlink()

	return bam_file

def create_mappability_track(
	input_fasta,
	reference,
	track_id,
	kmer_length,
	offset_step,
	chunk_size,
	n_threads,
	bwa_missing_prob_err_rate,
	bwa_max_gap_opens,
	bwa_seed_length,
	db_root="database",
	tags: Optional[List[str]] = None,
	manual_track_id=None,
):
	"""Pipeline to create a mappability track from an input FASTA file by aligning k-mers to a reference genome and exporting BigWig files.
	
	Args:
		input_fasta: Path to the input FASTA file containing sequences to analyze.
		reference: Path to the reference genome FASTA file that has been indexed with BWA.
		track_id: Optional string to use as the track ID in metadata and naming. If not provided, it will be derived from the input FASTA name.
		kmer_length: Length of the k-mers to generate from the input FASTA for alignment.
		offset_step: Step size for the sliding window when generating k-mers. Default is 1 for fully overlapping k-mers.
		chunk_size: Number of k-mers to include in each chunk for parallel alignment. Must be a positive integer.
		n_threads: Number of threads to use for parallel alignment.
		bwa_missing_prob_err_rate: The -n parameter for bwa aln, controlling the maximum edit distance.
		bwa_max_gap_opens: The -o parameter for bwa aln, controlling the maximum number of gap opens.
		bwa_seed_length: The -l parameter for bwa aln, controlling the seed length.
		db_root: Root directory for the database where target-specific directories and track outputs will be saved. Default is "database".
		tags: Optional list of tags to associate with the track in metadata.
		manual_track_id: Optional string to use as the track ID in metadata and naming, overriding automatic derivation.
	Returns:
		Path to the directory where the generated track and its metadata are saved.
	"""

	# Validate inputs and prepare paths
	input_fasta = Path(input_fasta)
	input_target = Path(input_target)
	target_name = input_target.stem
	input_name = input_fasta.stem

	if chunk_size < 1:
		raise ValueError("chunk_size must be a positive integer")
	query_track_id = track_id.strip() if track_id else input_name
	if not query_track_id:
		raise ValueError("track_id cannot be empty")
	
	input_md5 = file_md5(input_fasta)
	track_name = f"{query_track_id}_k{kmer_length}_w{offset_step}_n{float(bwa_missing_prob_err_rate):g}_o{bwa_max_gap_opens}_l{bwa_seed_length}"
	target_dir = Path(db_root) / target_name
	target_dir.mkdir(parents=True, exist_ok=True)
	out_dir = Path(db_root) / target_name / track_name
	out_dir.mkdir(parents=True, exist_ok=True)

	# Persist reference fingerprint once at target level and validate it on subsequent runs.
	target_md5 = file_md5(input_target)
	target_meta_yaml = target_dir / f"{target_name}.yaml"
	target_md5_file = target_dir / f"{target_name}.md5"

	existing_md5 = None
	if target_meta_yaml.exists():
		with open(target_meta_yaml, "r", encoding="utf-8") as f:
			existing_meta = yaml.safe_load(f) or {}
		if isinstance(existing_meta, dict):
			existing_md5 = existing_meta.get("reference_fasta_md5")

	if not existing_md5 and target_md5_file.exists():
		with open(target_md5_file, "r", encoding="utf-8") as f:
			first_token = (f.readline().strip().split() or [None])[0]
			existing_md5 = first_token

	if existing_md5 and existing_md5 != target_md5:
		raise ValueError(
			f"Reference MD5 mismatch for '{target_name}'. "
			f"Stored={existing_md5}, provided={target_md5}. "
			"Use the exact same reference FASTA used to build this target database."
		)

	target_meta_content = {
		"reference_name": target_name,
		"reference_fasta": str(input_target.resolve()),
		"reference_fasta_md5": target_md5,
		"last_updated": datetime.now().isoformat(timespec="seconds"),
	}
	with open(target_meta_yaml, "w", encoding="utf-8") as f:
		yaml.safe_dump(target_meta_content, f, sort_keys=False)
	with open(target_md5_file, "w", encoding="utf-8") as f:
		f.write(f"{target_md5}  {input_target.name}\n")

	target_seq_sizes = target_dir / f"{target_name}.sizes"
	if not target_seq_sizes.exists():
		log(f"Generating sequence sizes file for reference '{target_name}'...", "I")
		write_seq_sizes_from_bam(input_target, target_seq_sizes)

	# Create temporary directory in cluster TMPDIR or /tmp
	tmp_root = Path(os.environ.get("TMPDIR", "/tmp"))
	user = getpass.getuser()
	tmp_dir = tmp_root / f"{user}_wizardeye"
	tmp_dir.mkdir(parents=True, exist_ok=True)

	# Index with BWA if needed
	bwt_file = input_target.with_suffix(input_target.suffix + ".bwt")
	if not bwt_file.exists():
		log(f"BWA index not found for {input_target}, running bwa index...", "I")
		run(["bwa", "index", str(input_target)], check=True)
	else:
		log(f"BWA index found for {input_target}, skipping index.", "I")

	# Generate k-mers, deduplicate identical sequences by default, then split to chunks.
	raw_kmer_fasta = tmp_dir / f"{input_name}_kseq.raw.fasta"
	kmer_fasta = tmp_dir / f"{input_name}_kseq.dedup.fasta"
	seqkit_cmd = [
		"seqkit", "sliding",
		"-W", str(kmer_length),
		"-s", str(offset_step),
		str(input_fasta),
		"-o", str(raw_kmer_fasta)
	]
	log(f"Splitting {input_fasta} into {kmer_length}-mers with a sliding window of {offset_step}...", "I")
	run(seqkit_cmd, check=True)

	rmdup_cmd = [
		"seqkit",
		"rmdup",
		"-s",
		str(raw_kmer_fasta),
		"-o",
		str(kmer_fasta),
	]
	log("Deduplicating k-mers with seqkit rmdup -s...", "I")
	run(rmdup_cmd, check=True)

	log(f"Splitting k-mer FASTA into chunks of {chunk_size} sequences for parallel alignment...", "I")
	chunk_dir = tmp_dir / f"{input_name}_chunks"
	chunk_dir.mkdir(parents=True, exist_ok=True)
	chunk_cmd = [
		"seqkit",
		"split2",
		"-s",
		str(chunk_size),
		"-O",
		str(chunk_dir),
		"--extension",
		".fasta",
		str(kmer_fasta),
	]
	run(chunk_cmd, check=True)

	for tmp_fasta in (raw_kmer_fasta, kmer_fasta):
		if tmp_fasta.exists():
			tmp_fasta.unlink()

	chunk_fastas = sorted(chunk_dir.rglob("*.fasta"))
	if not chunk_fastas:
		raise RuntimeError("No chunk FASTA files were produced by seqkit split2")

	workers = max(1, n_threads)
	if n_threads > 1:
		log(f"Aligning {len(chunk_fastas)} chunks with {workers} parallel workers...", "I")

	chunk_bams: List[Path] = []
	with ThreadPoolExecutor(max_workers=workers) as executor:
		futures = [
			executor.submit(
				align_with_bwa_aln,
				chunk_fasta,
				input_target,
				bwa_missing_prob_err_rate,
				bwa_max_gap_opens,
				bwa_seed_length,
			)
			for chunk_fasta in chunk_fastas
		]
		for future in futures:
			chunk_bams.append(future.result())

	# Merge all chunk-level BAM files into one temporary BAM.
	bam_file = tmp_dir / f"{input_name}_{kmer_length}{f'_s{offset_step}' if offset_step != 1 else ''}_{Path(input_target).stem}.bam"
	
	log(f"Merging {len(chunk_bams)} chunk BAM files into one unique BAM...", "I")
	log(f"samtools cat {' '.join(str(path) for path in chunk_bams)}", "C")
	samtools_cat = subprocess.Popen(
		["samtools", "cat", *[str(path) for path in chunk_bams]],
		stdout=subprocess.PIPE,
	)
	try:
		run(
			["samtools", "sort", "-@", str(workers), "-o", str(bam_file), "-"],
			check=True,
			stdin=samtools_cat.stdout,
		)
	finally:
		if samtools_cat.stdout is not None:
			samtools_cat.stdout.close()

	return_code = samtools_cat.wait()
	if return_code != 0:
		raise subprocess.CalledProcessError(return_code, samtools_cat.args)

	for tmp_bam in chunk_bams:
		if tmp_bam.exists():
			tmp_bam.unlink()

	# Export only final BigWig depth tracks.
	log("Exporting mappability BigWig tracks from temporary BAM...", "I")
	cov_map_bw, cov_uniq_bw, covered_bp = from_alignment_to_bigWig(
		bam_file=bam_file,
		out_dir=out_dir,
		kmer_length=kmer_length,
	)

	# Save parameters
	log("Saving track parameters and metadata...", "I")
	param_yaml = out_dir / "param.yaml"
	param_content = {
		"generation_date": datetime.now().isoformat(timespec="seconds"),
		"wizardeye_version": PACKAGE_VERSION,
		"reference": str(input_target),
		"reference_fasta_md5": target_md5,
		"track_id": str(manual_track_id) if manual_track_id else None,
		"input": str(input_fasta),
		"input_fasta_md5": input_md5,
		"tags": from_charlist_to_list(tags, lowercase=True),
		"kmer_size": kmer_length,
		"sliding_window": offset_step,
		"chunk_size": chunk_size,
		"mapping_tool": "bwa aln",
		"bwa_parameters": {
			"-n": bwa_missing_prob_err_rate,
			"-o": bwa_max_gap_opens,
			"-l": bwa_seed_length,
			"-t": n_threads,
		},
		"covered_bases": {
			"map_all_bp": covered_bp["map_all"],
			"map_uniq_bp": covered_bp["map_uniq"],
		},
	}
	if param_content["track_id"] is None:
		param_content.pop("track_id")
	with open(param_yaml, "w") as f:
		yaml.safe_dump(param_content, f, sort_keys=False)

	log("Cleaning up temporary files and directories...", "I")
	try:
		shutil.rmtree(tmp_dir)
	except Exception as e:
		log(f"Could not remove tmp dir {tmp_dir}: {e}", "WARN")

	# Remove stale files from older runs to keep only BigWig outputs.
	for stale_name in ("n_map.bed", "n_uniq.bed", "cov_map.bg", "cov_uniq.bg", "chrom.sizes"):
		stale_path = out_dir / stale_name
		if stale_path.exists():
			stale_path.unlink()

	# Remove persisted BAM leftovers from older runs: only BigWig and metadata are kept.
	for stale_bam in out_dir.glob("*.bam"):
		stale_bam.unlink()

	log(f"Track saved in: {out_dir}", "S")
	return out_dir