# -*- coding: utf-8 -*-

import subprocess
import pysam
import hashlib
import shutil

from typing import Dict, Iterable, List, Optional, Tuple
from collections import defaultdict
from pathlib import Path

# -- Standard utilities --
def log(message: str, type: str, colorful=True):
	"""Log a message with a specific type and optional color formatting."""
	if colorful:
		if type == "I":
			print(f"\033[1;36m[INFO]\033[0m {message}")
		elif type == "S":
			print(f"\033[1;32m[SUCCESS]\033[0m {message}")
		elif type == "W":
			print(f"\033[1;33m[WARN]\033[0m {message}")
		elif type == "E":
			print(f"\033[1;31m[ERROR]\033[0m {message}")
		elif type == "FE":
			print(f"\033[1;31m[ERROR]\033[0m {message}")
			raise RuntimeError(message)
		elif type == "SC":
			print(f"\033[1;90m[CMD]\033[0m {message}")
		elif type == "C":
			print(f"\033[0;90m[CMD] {message}\033[0m")
	else:
		print(f"[{type}] {message}")

def file_md5(file_path: Path, chunk_size: int = 1024 * 1024) -> str:
	"""Return MD5 checksum of a file, reading it in chunks."""
	hasher = hashlib.md5()
	with Path(file_path).open("rb") as handle:
		while True:
			chunk = handle.read(chunk_size)
			if not chunk:
				break
			hasher.update(chunk)
	return hasher.hexdigest()

def from_charlist_to_list(values: Optional[List[str]], lowercase: bool = False) -> List[str]:
	"""Split comma-separated values, trim spaces, drop empties, and deduplicate preserving order."""
	if not values:
		return []

	normalized: List[str] = []
	seen = set()
	for raw_value in values:
		for token in str(raw_value).split(","):
			clean = token.strip()
			if lowercase:
				clean = clean.lower()
			if not clean or clean in seen:
				continue
			normalized.append(clean)
			seen.add(clean)
	return normalized

def query_name_from_param(param_content: Dict) -> str:
    """Extract query/input name from param.yaml, with a safe fallback."""
    input_value = param_content.get("input")
    if input_value:
        return Path(str(input_value)).stem

    track_id_value = param_content.get("track_id")
    if track_id_value:
        return str(track_id_value)

# -- Bioinformatics utilities --

# SAM/BAM fields parsing utilities
def parse_xa_tag(xa_value: str) -> List[Tuple[str, int, str]]:
	"""Parse XA:Z entries into (chrom, pos, cigar) tuples."""
	hits: List[Tuple[str, int, str]] = []
	if not xa_value:
		return hits

	for raw_hit in xa_value.split(";"):
		if not raw_hit:
			continue
		parts = raw_hit.split(",")
		if len(parts) < 3:
			continue
		chrom = parts[0]
		try:
			pos = abs(int(parts[1]))
		except ValueError:
			continue
		cigar = parts[2]
		hits.append((chrom, pos, cigar))

	return hits


def read_bam_sq_lengths(bam_path: Path) -> Dict[str, int]:
    """Read BAM header @SQ entries and return {SN: LN}."""
    header = subprocess.check_output(["samtools", "view", "-H", str(bam_path)], text=True)
    lengths: Dict[str, int] = {}

    for raw_line in header.splitlines():
        if not raw_line.startswith("@SQ"):
            continue
        fields = raw_line.split("\t")
        sn = None
        ln = None
        for field in fields[1:]:
            if field.startswith("SN:"):
                sn = field[3:]
            elif field.startswith("LN:"):
                ln = field[3:]
        if not sn or ln is None:
            continue
        lengths[sn] = int(ln)

    if not lengths:
        raise ValueError(f"No @SQ entries found in BAM header: {bam_path}")

    return lengths


def read_seq_sizes(seq_sizes_path: Path) -> Dict[str, int]:
    """Read sizes file and return {seq: length}."""
    sizes: Dict[str, int] = {}
    with seq_sizes_path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                raise ValueError(f"Invalid sizes line in {seq_sizes_path}: {line}")
            sizes[parts[0]] = int(parts[1])

    if not sizes:
        raise ValueError(f"Reference sizes is empty: {seq_sizes_path}")

    return sizes

def reference_len_from_cigar(cigar: Optional[str], default_k: int) -> int:
	"""Return reference-consuming length from a CIGAR string."""
	if not cigar:
		return default_k

	ref_len = 0
	num = []
	for char in cigar:
		if char.isdigit():
			num.append(char)
			continue
		if char in {"M", "D", "N", "=", "X"}:
			try:
				ref_len += int("".join(num))
			except ValueError:
				pass
		num = []

	return ref_len if ref_len > 0 else default_k


def count_covered_bases_from_bedgraph(bedgraph_path: Path) -> int:
	"""Return the number of covered base pairs from a bedGraph file (depth > 0)."""
	total_bp = 0
	with bedgraph_path.open("r", encoding="utf-8") as handle:
		for line in handle:
			parts = line.strip().split("\t")
			if len(parts) < 4:
				continue
			try:
				start = int(parts[1])
				end = int(parts[2])
				depth = float(parts[3])
			except ValueError:
				continue
			if end > start and depth > 0:
				total_bp += end - start
	return total_bp

# File management utilities
def write_bed(intervals: List[Tuple[str, int, int]], out_path: Path) -> Path:
    """Write BED intervals to out_path."""
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", encoding="utf-8") as handle:
        for chrom, start, end in intervals:
            handle.write(f"{chrom}\t{start}\t{end}\n")
    return out_path

def write_bedgraph(intervals: Iterable[Tuple[str, int, int]], out_path: Path) -> None:
	"""Write depth bedGraph (chrom, start, end, depth) from intervals."""
	events_by_chrom: Dict[str, Dict[int, int]] = defaultdict(lambda: defaultdict(int))

	for chrom, start, end in intervals:
		if end <= start:
			continue
		events_by_chrom[chrom][start] += 1
		events_by_chrom[chrom][end] -= 1

	with out_path.open("w", encoding="utf-8") as handle:
		for chrom in sorted(events_by_chrom):
			events = events_by_chrom[chrom]
			positions = sorted(events)
			depth = 0
			for idx, pos in enumerate(positions):
				depth += events[pos]
				if idx + 1 >= len(positions):
					continue
				next_pos = positions[idx + 1]
				if depth > 0 and next_pos > pos:
					handle.write(f"{chrom}\t{pos}\t{next_pos}\t{depth}\n")

def write_seq_sizes_from_bam(bam_file: Path, out_path: Path) -> None:
	"""Write sequences sizes file required by bedGraphToBigWig."""
	if pysam is None:
		raise RuntimeError(
			"pysam is required to export BigWig files. Install it with: pip install pysam"
		)

	with pysam.AlignmentFile(str(bam_file), "rb") as bam, out_path.open("w", encoding="utf-8") as handle:
		for chrom, length in sorted(zip(bam.references, bam.lengths), key=lambda x: x[0]):
			handle.write(f"{chrom}\t{length}\n")


def write_seq_sizes_from_fasta(reference_fasta: Path, output_sizes: Path) -> None:
	"""Write a chrom.sizes-style file (chrom\tlength) from FASTA using samtools faidx."""
	if shutil.which("samtools") is None:
		raise RuntimeError("samtools not found in PATH, cannot generate reference chrom sizes")

	output_sizes.parent.mkdir(parents=True, exist_ok=True)
	fai_path = Path(f"{reference_fasta}.fai")

	log(f"samtools faidx {reference_fasta}", "C")
	subprocess.run(["samtools", "faidx", str(reference_fasta)], check=True)

	if not fai_path.exists():
		raise RuntimeError(f"samtools faidx did not produce index: {fai_path}")

	with fai_path.open("r", encoding="utf-8") as fai, output_sizes.open("w", encoding="utf-8") as out:
		for raw_line in fai:
			line = raw_line.strip()
			if not line:
				continue
			parts = line.split("\t")
			if len(parts) < 2:
				raise ValueError(f"Invalid FASTA index line in {fai_path}: {line}")
			out.write(f"{parts[0]}\t{parts[1]}\n")


def convert_bedgraph_to_bigwig(bedgraph_path: Path, seq_sizes_path: Path, bigwig_path: Path) -> None:
	"""Convert a bedGraph file to BigWig using bedGraphToBigWig."""
	if shutil.which("bedGraphToBigWig") is None:
		raise RuntimeError("bedGraphToBigWig not found in PATH, cannot produce .bw outputs")

	log(f"bedGraphToBigWig {bedgraph_path} {seq_sizes_path} {bigwig_path}", "C")
	subprocess.run(
		[
			"bedGraphToBigWig",
			str(bedgraph_path),
			str(seq_sizes_path),
			str(bigwig_path),
		],
		check=True,
	)

# BAM processing utilities
def get_mapping_intervals(bam_file: Path, kmer_length: int) -> Iterable[Tuple[str, int, int]]:
	"""Yield intervals for primary alignments and all XA-reported alternative hits."""
	with pysam.AlignmentFile(str(bam_file), "rb") as bam:
		for read in bam.fetch(until_eof=True):
			if read.is_unmapped or read.reference_start < 0:
				continue

			primary_len = read.reference_length if read.reference_length and read.reference_length > 0 else kmer_length
			start = read.reference_start
			yield (read.reference_name, start, start + primary_len)

			if read.has_tag("XA"):
				xa_value = read.get_tag("XA")
				for chrom, pos_1based, cigar in parse_xa_tag(xa_value):
					ref_len = reference_len_from_cigar(cigar, kmer_length)
					alt_start = pos_1based - 1
					yield (chrom, alt_start, alt_start + ref_len)


def get_unique_mapping_intervals(bam_file: Path, kmer_length: int) -> Iterable[Tuple[str, int, int]]:
	"""Yield intervals for mapped reads with MAPQ > 0 and no XA alternatives."""
	with pysam.AlignmentFile(str(bam_file), "rb") as bam:
		for read in bam.fetch(until_eof=True):
			if read.is_unmapped or read.reference_start < 0 or read.mapping_quality <= 0:
				continue
			if read.has_tag("XA") and read.get_tag("XA"):
				continue

			ref_len = read.reference_length if read.reference_length and read.reference_length > 0 else kmer_length
			start = read.reference_start
			yield (read.reference_name, start, start + ref_len)

def validate_initial_bam_reference_compatibility(
	bam_path: Path,
	reference_seq_sizes_path: Path,
) -> None:
	"""Fail fast by checking BAM @SQ SN/LN against reference sequence sizes."""
	if not reference_seq_sizes_path.exists() or not reference_seq_sizes_path.is_file():
		raise FileNotFoundError(
			"Reference sequence sizes not found for compatibility check: "
			f"{reference_seq_sizes_path}. Recreate or re-import this reference."
		)

	bam_sizes = read_bam_sq_lengths(bam_path)
	ref_sizes = read_seq_sizes(reference_seq_sizes_path)

	missing_in_ref = sorted([contig for contig in bam_sizes if contig not in ref_sizes])
	length_mismatches = sorted(
		[
			(contig, bam_sizes[contig], ref_sizes[contig])
			for contig in bam_sizes
			if contig in ref_sizes and bam_sizes[contig] != ref_sizes[contig]
		],
		key=lambda x: x[0],
	)

	if missing_in_ref:
		log(f"Sequences missing in reference: {', '.join(missing_in_ref)}.", "E")

	if length_mismatches:
		mismatch_lengths = ', '.join(
			[f"{contig}(bam={bam_len},ref={ref_len})" for contig, bam_len, ref_len in length_mismatches[:8]]
		)

		log(f"Length mismatches: {mismatch_lengths}.", "E")

	if not missing_in_ref and not length_mismatches:
		return

	raise ValueError(
		"Initial files are incompatible. BAM header @SQ (SN/LN) does not match reference sequence sizes/lengths. \n"
		"Please ensure the reference FASTA and BAM files are compatible and if so, edit the BAM file to match the reference sequences names and lengths."
	)

