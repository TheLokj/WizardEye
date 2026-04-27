# -*- coding: utf-8 -*-

import subprocess
import shutil
import yaml
import tempfile
import os
import getpass

from pathlib import Path
from datetime import datetime
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
from typing import Dict, Iterable, List, Optional, Tuple

from .version import PACKAGE_VERSION
from .utils import log
from .utils import file_md5
from .utils import (
	count_covered_bases_from_bedgraph,
	from_charlist_to_list,
	write_bedgraph,
	convert_bedgraph_to_bigwig,
	write_seq_sizes_from_bam,
	get_mapping_intervals,
	get_unique_mapping_intervals,
)

try:
	import pysam
except ImportError:
	pysam = None

# -- WizardEye mappability utilities --

def export_mapping_bigwig_files(bam_file: Path, out_dir: Path, kmer_length: int) -> Tuple[Path, Path, Dict[str, int]]:
	"""Export BigWig tracks and return paths with covered-bases metrics for both masks."""
	if pysam is None:
		raise RuntimeError(
			"pysam is required to export BigWig files. Install it with: pip install pysam"
		)

	cov_map_bg = out_dir / "map_all.bg"
	cov_uniq_bg = out_dir / "map_uniq.bg"
	cov_map_bw = out_dir / "map_all.bw"
	cov_uniq_bw = out_dir / "map_uniq.bw"
	seq_sizes = out_dir / "chrom.sizes"

	write_bedgraph(get_mapping_intervals(bam_file, kmer_length), cov_map_bg)
	write_bedgraph(get_unique_mapping_intervals(bam_file, kmer_length), cov_uniq_bg)
	covered_bp = {
		"map_all": count_covered_bases_from_bedgraph(cov_map_bg),
		"map_uniq": count_covered_bases_from_bedgraph(cov_uniq_bg),
	}
	write_seq_sizes_from_bam(bam_file, seq_sizes)

	convert_bedgraph_to_bigwig(cov_map_bg, seq_sizes, cov_map_bw)
	convert_bedgraph_to_bigwig(cov_uniq_bg, seq_sizes, cov_uniq_bw)

	for tmp_file in (cov_map_bg, cov_uniq_bg, seq_sizes):
		if tmp_file.exists():
			tmp_file.unlink()

	return cov_map_bw, cov_uniq_bw, covered_bp

# -- WizardEye mappability pipeline blocks --

def align_chunk_to_bam(
	chunk_fasta: Path,
	input_target: Path,
	bwa_missing_prob_err_rate: float,
	bwa_max_gap_opens: int,
	bwa_seed_length: int,
) -> Path:
	"""Align one chunk FASTA and return the produced BAM path."""
	sai_file = chunk_fasta.with_suffix(".sai")
	bam_file = chunk_fasta.with_suffix(".bam")

	with sai_file.open("w", encoding="utf-8") as sai_out:
		log(f"bwa aln -t 1 -n {bwa_missing_prob_err_rate} -o {bwa_max_gap_opens} -l {bwa_seed_length} {input_target} {chunk_fasta}", "C")
		subprocess.run(
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
				str(input_target),
				str(chunk_fasta),
			],
			check=True,
			stdout=sai_out,
		)

	with bam_file.open("wb") as bam_out:
		log(f"bwa samse -n 1000000 {input_target} {sai_file} {chunk_fasta}", "C")
		bwa_samse = subprocess.Popen(
			[
				"bwa",
				"samse",
				"-n",
				"1000000",
				str(input_target),
				str(sai_file),
				str(chunk_fasta),
			],
			stdout=subprocess.PIPE,
		)
		try:
			log(f"samtools view -b -F 4 -", "C")
			subprocess.run(
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

	for tmp_file in (sai_file, chunk_fasta):
		if tmp_file.exists():
			tmp_file.unlink()

	return bam_file

def create_mappability_track(
	input_fasta,
	input_target,
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
		log(f"bwa index {input_target}", "C")
		subprocess.run(["bwa", "index", str(input_target)], check=True)
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
	log(f"seqkit sliding -W {kmer_length} -s {offset_step} {input_fasta} -o {raw_kmer_fasta}", "C")
	subprocess.run(seqkit_cmd, check=True)

	rmdup_cmd = [
		"seqkit",
		"rmdup",
		"-s",
		str(raw_kmer_fasta),
		"-o",
		str(kmer_fasta),
	]
	log("Deduplicating k-mers with seqkit rmdup -s (default behavior)...", "I")
	log(f"seqkit rmdup -s {raw_kmer_fasta} -o {kmer_fasta}", "C")
	subprocess.run(rmdup_cmd, check=True)

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
	log(f"seqkit split2 -s {chunk_size} -O {chunk_dir} --extension .fasta {kmer_fasta}", "C")
	subprocess.run(chunk_cmd, check=True)

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
				align_chunk_to_bam,
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
		log(f"samtools sort -@ {workers} -o {bam_file} -", "C")
		subprocess.run(
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
	cov_map_bw, cov_uniq_bw, covered_bp = export_mapping_bigwig_files(
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


