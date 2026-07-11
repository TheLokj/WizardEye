# -*- coding: utf-8 -*-

"""Utility functions for WizardEye

This modules contains various utility functions for WizardEye, including subprocess execution with logging,
file hashing, sequence size handling, BAM metadata parsing, and BED file merging.
"""

import subprocess
import pysam
import hashlib
import shutil
import shlex
import sys
import os

from typing import Dict, Iterable, List, Optional, Tuple
from dataclasses import dataclass
from pathlib import Path


def get_bwa_params_hash(bwa_params: "BWAParameters", hash_length: int = 8) -> str:
    """Generate a short hash from BWA parameters for track naming.

    Args:
        bwa_params: BWAParameters object containing alignment parameters.
        hash_length: Length of the hash string to return (default: 8).

    Returns:
        A short hash string representing the BWA parameters.
    """
    # Create a string representation of all relevant BWA parameters
    params_str = f"{bwa_params.missing_prob_err_rate}:{bwa_params.max_gap_opens}:{bwa_params.seed_length}:{bwa_params.all_aln}:{bwa_params.threads}:{bwa_params.r_best_hits}:{bwa_params.samse_n}"

    # Generate MD5 hash and return first N characters
    hash_obj = hashlib.md5(params_str.encode())
    return hash_obj.hexdigest()[:hash_length]


# --- WizardEye development utilities ---


def log(message: str, type: str, colorful=True):
    """Log a message with a type prefix and optional color coding.

    Type can be one of: I (info), S (success), W (warning), E (error), C (subprocess command), SC (subprocess command output).

    Args:
            message (str): The message to log.
            type (str): The type of message, determining prefix and color.
            colorful (bool): Whether to use ANSI color codes in the output.
    """
    colorful = sys.stdout.isatty() and sys.stderr.isatty()
    clear_clr = "\033[0m"

    if type == "I":
        i_clr = "\033[1;36m"
        print(
            f"{i_clr}[INFO]{clear_clr} {message}" if colorful else f"[INFO] {message}"
        )
    elif type == "S":
        s_clr = "\033[1;32m"
        print(
            f"{s_clr}[SUCCESS]{clear_clr} {message}"
            if colorful
            else f"[SUCCESS] {message}"
        )
    elif type == "W":
        w_clr = "\033[1;33m"
        print(
            f"{w_clr}[WARN]{clear_clr} {message}" if colorful else f"[WARN] {message}"
        )
    elif type == "E":
        e_clr = "\033[1;31m"
        print(
            f"{e_clr}[ERROR]{clear_clr} {message}" if colorful else f"[ERROR] {message}"
        )
    elif type == "C":
        c_clr = "\033[0;90m"
        print(f"{c_clr}[CMD] {message} {clear_clr}" if colorful else f"[CMD] {message}")
    elif type == "SC":
        sc_clr = "\033[1;90m"
        print(
            f"{sc_clr}[SUB-CMD] {message} {clear_clr}"
            if colorful
            else f"[SUB-CMD] {message}"
        )
    else:
        print(f"[{type}] {message}")


def run(
    command: List[str], log_output: bool = False, **kwargs
) -> subprocess.CompletedProcess:
    """Execute a command with subprocess.run and log its output.

    Args:
            command (List[str]): Command and arguments.
            log_output (bool): Capture and relay stdout/stderr with log() when possible.
            **kwargs: Forwarded as-is to subprocess.run (e.g. check=True, text=True).

    Returns:
            subprocess.CompletedProcess: The subprocess.run return value.
    """

    def _log_stream(content: Optional[object], prefix: str) -> None:
        if content is None:
            return
        text = str(content).rstrip()
        if not text:
            return
        for line in text.splitlines():
            log(f"[{prefix}] {line}", "SC")

    log(" ".join(shlex.quote(str(arg)) for arg in command), "C")
    run_kwargs = dict(kwargs)

    if log_output:
        # If caller did not specify stream behavior, capture outputs so they can be logged.
        if (
            not run_kwargs.get("capture_output")
            and "stdout" not in run_kwargs
            and "stderr" not in run_kwargs
        ):
            run_kwargs["capture_output"] = True
        # Prefer text output in logs unless caller already requested bytes behavior.
        if "text" not in run_kwargs and "encoding" not in run_kwargs:
            run_kwargs["text"] = True

    try:
        result = subprocess.run(command, **run_kwargs)
    except subprocess.CalledProcessError as exc:
        if log_output:
            _log_stream(exc.stdout, "stdout")
            _log_stream(exc.stderr, "stderr")
        raise

    if log_output:
        _log_stream(result.stdout, "stdout")
        _log_stream(result.stderr, "stderr")

    return result


# --- File, data and string utilities ---


def file_md5(file_path: Path, chunk_size: int = 1024 * 1024) -> str:
    """Return MD5 checksum of a file, reading it in chunks.

    Args:
            file_path (Path): The path to the file to hash.
            chunk_size (int): The size of chunks to read at a time (default: 1MB)."""
    hasher = hashlib.md5()
    with Path(file_path).open("rb") as handle:
        while True:
            chunk = handle.read(chunk_size)
            if not chunk:
                break
            hasher.update(chunk)
    return hasher.hexdigest()


def from_charlist_to_list(
    values: Optional[List[str]], lowercase: bool = False
) -> List[str]:
    """Split comma-separated values, trim spaces, drop empties, and deduplicate preserving order.

    Args:
            values (Optional[List[str]]): A list of strings, each potentially containing comma-separated values.
            lowercase (bool): Whether to convert values to lowercase (default: False).

    Returns:
            List[str]: A list of unique, cleaned values.
    """
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


def get_name_from_param(param_content: Dict) -> str:
    """Extract a name from parameter content, preferring 'input' then 'track_id'.

    Args:
            param_content (Dict): The parameter content dictionary.

    Returns:
            str: The extracted name.
    """
    input_value = param_content.get("input")
    if input_value:
        return Path(str(input_value)).stem

    track_id_value = param_content.get("track_id")
    if track_id_value:
        return str(track_id_value)


# --- Sequence utilities ---


def get_seq_sizes(seq_sizes_path: Path) -> Dict[str, int]:
    """Get sequence sizes from a chrom.sizes-style file (chrom\\tlength).

    Args:
            seq_sizes_path (Path): The path to the sequence sizes file.

    Returns:
            Dict[str, int]: A dictionary mapping sequence names to their lengths.
    """
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


def write_seq_sizes_from_fasta(reference_fasta: Path, output_sizes: Path) -> Path:
    """Write a chrom.sizes-style file (chrom\\tlength) from FASTA using samtools faidx.

    Args:
            reference_fasta (Path): The path to the reference FASTA file.
            output_sizes (Path): The path to write the chrom.sizes file to.

    Returns:
            Path: The path to the generated chrom.sizes file.

    Raises:
            RuntimeError: If samtools is not found in PATH or if faidx does not produce the expected index file.
            ValueError: If the FASTA index file contains invalid lines.
    """
    if shutil.which("samtools") is None:
        raise RuntimeError(
            "samtools not found in PATH, cannot generate reference chrom sizes"
        )

    output_sizes.parent.mkdir(parents=True, exist_ok=True)
    fai_path = Path(f"{reference_fasta}.fai")

    run(["samtools", "faidx", str(reference_fasta)], check=True)

    if not fai_path.exists():
        raise RuntimeError(f"samtools faidx did not produce index: {fai_path}")

    with fai_path.open("r", encoding="utf-8") as fai, output_sizes.open(
        "w", encoding="utf-8"
    ) as out:
        for raw_line in fai:
            line = raw_line.strip()
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                raise ValueError(f"Invalid FASTA index line in {fai_path}: {line}")
            out.write(f"{parts[0]}\t{parts[1]}\n")

    return output_sizes


# --- SAM/BAM utilities ---


def parse_xa_tag(xa_value: str) -> List[Tuple[str, int, str]]:
    """Parse XA:Z entries into (chrom, pos, cigar) tuples.

    Args:
            xa_value (str): The raw string from the XA tag, e.g. "chr1,1000,100M;chr2,2000,150M"

    Returns:
            List[Tuple[str, int, str]]: A list of (chrom, pos, cigar) tuples.
    """
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


def parse_bam_metadata(bam_path: Path) -> Dict[str, object]:
    """Read BAM header metadata: @SQ lengths and latest bwa aln command/options.

    Args:
            bam_path (Path): The path to the BAM file to parse.

    Returns:
            Dict[str, object]: A dictionary containing:
                    - "sq_lengths": Dict[str, int] mapping sequence names to lengths from @SQ.
                    - "bwa_aln_cmd": Optional[str] of the latest bwa aln command found in @PG CL.
                    - "bwa_aln_options": Dict[str, Optional[str]] of extracted bwa aln options (-n, -o, -l) from the command.
    """
    log(f"samtools view -H {bam_path}", "C")
    header = subprocess.check_output(
        ["samtools", "view", "-H", str(bam_path)], text=True
    )

    lengths: Dict[str, int] = {}
    bwa_aln_cmds: List[str] = []

    # Get sequence names and lengths
    for raw_line in header.splitlines():
        if raw_line.startswith("@SQ"):
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
            continue

        if not raw_line.startswith("@PG"):
            continue

        fields = raw_line.split("\t")
        cl = None
        for field in fields[1:]:
            if field.startswith("CL:"):
                cl = field[3:]
        if not cl:
            continue

        try:
            tokens = shlex.split(cl)
        except ValueError:
            tokens = cl.split()

        if len(tokens) < 2:
            continue

        has_bwa_aln = False
        for idx in range(len(tokens) - 1):
            exec_name = Path(tokens[idx]).name
            if "bwa" in exec_name and tokens[idx + 1] == "aln":
                has_bwa_aln = True
                break

        if not has_bwa_aln:
            continue

        bwa_aln_cmds.append(cl)

    if not lengths:
        raise ValueError(f"No @SQ entries found in BAM header: {bam_path}")

    # Get the latest bwa aln command and parse options of interest
    bwa_aln_cmd = bwa_aln_cmds[-1] if bwa_aln_cmds else None
    bwa_aln_options: Dict[str, Optional[str]] = {"-n": None, "-o": None, "-l": None}

    if bwa_aln_cmd:
        try:
            tokens = shlex.split(bwa_aln_cmd)
        except ValueError:
            tokens = bwa_aln_cmd.split()

        aln_start_idx = None
        for idx in range(len(tokens) - 1):
            exec_name = Path(tokens[idx]).name
            if "bwa" in exec_name and tokens[idx + 1] == "aln":
                aln_start_idx = idx + 2
                break

        option_tokens = tokens[aln_start_idx:] if aln_start_idx is not None else tokens

        i = 0
        while i < len(option_tokens):
            token = option_tokens[i]
            for opt in bwa_aln_options:
                if token == opt and i + 1 < len(option_tokens):
                    bwa_aln_options[opt] = option_tokens[i + 1]
                elif token.startswith(opt) and len(token) > len(opt):
                    bwa_aln_options[opt] = token[len(opt) :]
            i += 1

    return {
        "sq_lengths": lengths,
        "bwa_aln_cmd": bwa_aln_cmd,
        "bwa_aln_options": bwa_aln_options,
    }


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


def iterate_mapping_intervals(
    bam_file: Path, kmer_length: int
) -> Iterable[Tuple[str, int, int]]:
    """Iterate mapping intervals from BAM, including primary and XA alternative mappings.

    Args:
            bam_file (Path): The path to the BAM file to parse.
            kmer_length (int): The k-mer length to use as a default interval length when reference length cannot be determined from the read or CIGAR.

    Yields:
            Tuple[str, int, int]: (chrom, start, end) tuples representing mapping intervals, including primary and XA alternatives.
    """
    with pysam.AlignmentFile(str(bam_file), "rb") as bam:
        for read in bam.fetch(until_eof=True):
            if read.is_unmapped or read.reference_start < 0:
                continue

            primary_len = (
                read.reference_length
                if read.reference_length and read.reference_length > 0
                else kmer_length
            )
            start = read.reference_start
            yield (read.reference_name, start, start + primary_len)

            if read.has_tag("XA"):
                xa_value = read.get_tag("XA")
                for chrom, pos_1based, cigar in parse_xa_tag(xa_value):
                    ref_len = reference_len_from_cigar(cigar, kmer_length)
                    alt_start = abs(pos_1based) - 1
                    yield (chrom, alt_start, alt_start + ref_len)


def iterate_unique_mapping_intervals(
    bam_file: Path, kmer_length: int
) -> Iterable[Tuple[str, int, int]]:
    """Iterate intervals for mapped reads with MAPQ > 0 and no XA alternatives.

    Args:
            bam_file (Path): The path to the BAM file to parse.
            kmer_length (int): The k-mer length to use as a default interval length when reference length cannot be determined from the read or CIGAR.

    Yields:
            Tuple[str, int, int]: (chrom, start, end) tuples representing unique mapping intervals.
    """
    with pysam.AlignmentFile(str(bam_file), "rb") as bam:
        for read in bam.fetch(until_eof=True):
            if (
                read.is_unmapped
                or read.reference_start < 0
                or read.mapping_quality <= 0
            ):
                continue
            if read.has_tag("XA") and read.get_tag("XA"):
                continue

            ref_len = (
                read.reference_length
                if read.reference_length and read.reference_length > 0
                else kmer_length
            )
            start = read.reference_start
            yield (read.reference_name, start, start + ref_len)


def validate_bam_compatibility(
    bam_path: Path,
    reference_seq_sizes_path: Path,
    bwa_params: Optional["BWAParameters"] = None,
) -> None:
    """Validate that BAM header @SQ SN/LN entries are compatible with reference sequence sizes, and optionally check for consistency with requested BWA aln parameters.

    Args:
            bam_path (Path): The path to the BAM file to validate.
            reference_seq_sizes_path (Path): The path to the reference sequence sizes file.
            bwa (Optional[BWAParameters]): The expected BWA parameters used during BWA alignment.

    Raises:
            FileNotFoundError: If the reference sequence sizes file does not exist.
            ValueError: If there are incompatibilities between BAM header and reference sequence sizes."""

    if not reference_seq_sizes_path.exists() or not reference_seq_sizes_path.is_file():
        raise FileNotFoundError(
            "Reference sequence sizes not found for compatibility check: "
            f"{reference_seq_sizes_path}. Recreate or re-import this reference."
        )

    log("Validating BAM header metadata against reference and BWA parameters...", "I")
    bam_metadata = parse_bam_metadata(bam_path)
    bam_sizes = bam_metadata["sq_lengths"]
    bwa_cmd = bam_metadata["bwa_aln_cmd"]
    observed = bam_metadata["bwa_aln_options"]

    expected_bwa_params: Dict[str, Optional[object]] = {
        "-n": bwa_params.missing_prob_err_rate if bwa_params else None,
        "-o": bwa_params.max_gap_opens if bwa_params else None,
        "-l": bwa_params.seed_length if bwa_params else None,
    }

    if any(value is not None for value in expected_bwa_params.values()):
        if not bwa_cmd:
            log(
                "No 'bwa aln' command found in BAM header (@PG). Cannot verify requested filtration parameters.",
                "W",
            )
        else:
            mismatches: List[str] = []
            for opt, expected_value in expected_bwa_params.items():
                if expected_value is None:
                    continue
                observed_value = observed[opt]
                if observed_value is None:
                    mismatches.append(f"{opt} missing (expected {expected_value})")
                    continue

                if opt == "-n":
                    try:
                        if abs(float(observed_value) - float(expected_value)) > 1e-12:
                            mismatches.append(
                                f"{opt}={observed_value} (expected {expected_value})"
                            )
                    except ValueError:
                        mismatches.append(
                            f"{opt}={observed_value} (expected {expected_value})"
                        )
                else:
                    try:
                        if int(observed_value) != int(expected_value):
                            mismatches.append(
                                f"{opt}={observed_value} (expected {expected_value})"
                            )
                    except ValueError:
                        mismatches.append(
                            f"{opt}={observed_value} (expected {expected_value})"
                        )

            if mismatches:
                log(
                    "`bwa` command used to compute the alignment differs from parameters used to compute masks: "
                    f"{'; '.join(mismatches)}.\
					\nHeader command: {bwa_cmd}",
                    "W",
                )

    ref_sizes = get_seq_sizes(reference_seq_sizes_path)

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
        mismatch_lengths = ", ".join(
            [
                f"{contig}(bam={bam_len},ref={ref_len})"
                for contig, bam_len, ref_len in length_mismatches[:8]
            ]
        )

        log(f"Length mismatches: {mismatch_lengths}.", "E")

    if not missing_in_ref and not length_mismatches:
        if bwa_cmd and not mismatches:
            log("BAM is compatible with reference and requested BWA parameters.", "I")
        else:
            log("BAM sequences are compatible with reference.", "I")
        return

    raise ValueError(
        "Initial files are incompatible. BAM header @SQ (SN/LN) does not match reference sequence sizes/lengths. \n"
        "Please ensure the reference FASTA and BAM files are compatible and if so, edit the BAM file to match the reference sequences names and lengths."
    )


def write_seq_sizes_from_bam(bam_file: Path, out_path: Path) -> Path:
    """Write sequences sizes file required by bedGraphToBigWig.

    Args:
            bam_file (Path): The path to the BAM file to read sequence names and lengths from.
            out_path (Path): The path to write the chrom.sizes file to.

    Returns:
            Path: The path to the written chrom.sizes file.

    Raises:
            RuntimeError: If pysam is not installed, which is required to read BAM files.
    """
    if pysam is None:
        raise RuntimeError(
            "pysam is required to export BigWig files. Install it with: pip install pysam"
        )

    with pysam.AlignmentFile(str(bam_file), "rb") as bam, out_path.open(
        "w", encoding="utf-8"
    ) as handle:
        for chrom, length in sorted(
            zip(bam.references, bam.lengths), key=lambda x: x[0]
        ):
            handle.write(f"{chrom}\t{length}\n")

    return out_path


def merge_and_sort_bams(
    input_bams: List[Path],
    output_bam: Path,
    n_threads: int = 1,
) -> None:
    """Merge multiple BAM files and sort the result using samtools cat + sort pipeline.

    Runs samtools cat to concatenate BAMs, then pipes through samtools sort.

    Args:
            input_bams (List[Path]): List of BAM files to merge.
            output_bam (Path): Path to write the sorted merged BAM to.
            n_threads (int): Number of threads for samtools sort. Default is 1.

    Raises:
            subprocess.CalledProcessError: If samtools cat or sort fails.
    """
    samtools_cat = subprocess.Popen(
        ["samtools", "cat"] + [str(path) for path in input_bams],
        stdout=subprocess.PIPE,
    )
    try:
        run(
            ["samtools", "sort", "-@", str(n_threads), "-o", str(output_bam), "-"],
            check=True,
            stdin=samtools_cat.stdout,
        )
    finally:
        if samtools_cat.stdout is not None:
            samtools_cat.stdout.close()
        return_code = samtools_cat.wait()
        if return_code != 0:
            raise subprocess.CalledProcessError(return_code, samtools_cat.args)


# --- BED and genomic interval utilities ---


def merge_bed_files(bed_files: List[Tuple[str, Path]], output_bed: Path) -> Path:
    """Merge multiple BED files with track names into a single BED file with merged intervals
    and comma-separated track annotations, using bedtools multiinter.

    Args:
        bed_files: List of (track_name, bed_path) tuples.
        output_bed: Destination BED file path.

    Returns:
        Path to the written output BED file.
    """
    output_bed.parent.mkdir(parents=True, exist_ok=True)

    if not bed_files:
        output_bed.write_text("", encoding="utf-8")
        return output_bed

    bedtools = shutil.which("bedtools")
    if bedtools is None:
        raise RuntimeError("bedtools is required but was not found on PATH. ")

    track_names = [name.split("_k")[0] for name, _ in bed_files]
    bed_paths = [str(p) for _, p in bed_files]

    if len(bed_files) == 1:
        log("Only one track requested.", "W")
        raw_name, bed_path = bed_files[0]
        track_name = raw_name.split("_k")[0]

        cmd = [bedtools, "merge", "-i", str(bed_path)]
        log(" ".join(cmd), "C")
        result = subprocess.run(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            check=True,
        )

        with output_bed.open("w", encoding="utf-8") as out_fh:
            for line in result.stdout.splitlines():
                if not line:
                    continue
                parts = line.split("\t")
                if len(parts) >= 3:
                    chrom, start, end = parts[0], parts[1], parts[2]
                    out_fh.write(f"{chrom}\t{start}\t{end}\t{track_name}\n")

    else:
        cmd = [bedtools, "multiinter", "-i"] + bed_paths + ["-names"] + track_names
        log(" ".join(cmd), "C")

        result = subprocess.run(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            check=True,
        )

        # multiinter output columns:
        #   0: chrom  1: start  2: end  3: num_intersecting  4: track_names (comma-sep)
        with output_bed.open("w", encoding="utf-8") as out_fh:
            for line in result.stdout.splitlines():
                if not line:
                    continue
                parts = line.split("\t")
                if len(parts) < 5:
                    continue
                chrom, start, end, tracks = parts[0], parts[1], parts[2], parts[4]
                seen: set[str] = set()
                deduped = []
                for t in tracks.split(","):
                    if t not in seen:
                        seen.add(t)
                        deduped.append(t)
                out_fh.write(f"{chrom}\t{start}\t{end}\t{','.join(deduped)}\n")

    return output_bed


def sort_bed_file(
    input_path: Path,
    output_path: Optional[Path] = None,
    n_threads: int = 1,
    tmp_dir: Optional[Path] = None,
) -> Path:
    """Sort a BED file by chromosome and position.

    Args:
            input_path (Path): The path to the input BED file to sort.
            output_path (Optional[Path]): The path to write the sorted output to. Defaults to input_path (in-place).
            n_threads (int): Number of threads for parallel sort. Default is 1.
            tmp_dir (Optional[Path]): Temporary directory for sort. Defaults to output_path parent.

    Returns:
            Path: The path to the sorted output file.
    """
    output_path = output_path or input_path
    tmp_dir = tmp_dir or output_path.parent
    run(
        [
            "sort",
            "-k1,1",
            "-k2,2n",
            f"--parallel={n_threads}",
            "-T",
            str(tmp_dir),
            str(input_path),
            "-o",
            str(output_path),
        ],
        check=True,
        env={**os.environ, "LC_ALL": "C"},
    )
    return output_path


def compute_sorted_genome_coverage(
    input_bed: Path,
    seq_sizes: Path,
    output_bg: Path,
    n_threads: int = 1,
    tmp_dir: Optional[Path] = None,
) -> None:
    """Compute genome coverage from BED and write sorted bedGraph output.

    Runs bedtools genomecov with -bga flag and pipes through sort to produce
    a position-sorted bedGraph file.

    Args:
            input_bed (Path): The path to the input BED file.
            seq_sizes (Path): The path to the chrom.sizes file.
            output_bg (Path): The path to write the sorted bedGraph output to.
            n_threads (int): Number of threads for parallel sort. Default is 1.
            tmp_dir (Optional[Path]): Temporary directory for sort. Defaults to output_bg parent.

    Raises:
            subprocess.CalledProcessError: If bedtools genomecov or sort fails.
    """
    tmp_dir = tmp_dir or output_bg.parent
    genomecov_cmd = subprocess.Popen(
        [
            "bedtools",
            "genomecov",
            "-i",
            str(input_bed),
            "-g",
            str(seq_sizes),
            "-bga",
        ],
        stdout=subprocess.PIPE,
    )
    sort_cmd = subprocess.Popen(
        [
            "sort",
            "-k1,1",
            "-k2,2n",
            f"--parallel={n_threads}",
            "-T",
            str(tmp_dir),
        ],
        stdin=genomecov_cmd.stdout,
        stdout=open(output_bg, "w"),
        env={**os.environ, "LC_ALL": "C"},
    )
    genomecov_cmd.stdout.close()
    return_code = sort_cmd.wait()
    if return_code != 0:
        raise subprocess.CalledProcessError(return_code, sort_cmd.args)
    if genomecov_cmd.wait() != 0:
        raise subprocess.CalledProcessError(
            genomecov_cmd.returncode, genomecov_cmd.args
        )


def convert_bedgraph_to_bigwig(
    bedgraph_path: Path, seq_sizes_path: Path, bigwig_path: Path
) -> Path:
    """Convert a bedGraph file to BigWig using bedGraphToBigWig.

    Args:
            bedgraph_path (Path): The path to the input bedGraph file.
            seq_sizes_path (Path): The path to the chrom.sizes file with sequence lengths.
            bigwig_path (Path): The path to write the output BigWig file.

    Returns:
            Path: The path to the generated BigWig file."""
    if shutil.which("bedGraphToBigWig") is None:
        raise RuntimeError(
            "bedGraphToBigWig not found in PATH, cannot produce .bw outputs"
        )

    run(
        ["bedGraphToBigWig", str(bedgraph_path), str(seq_sizes_path), str(bigwig_path)],
        check=True,
    )

    return bigwig_path


def convert_bigwig_to_bedGraph(bigwig_path: Path, bedgraph_path: Path) -> Path:
    """Convert a BigWig file to BedGraph using bigWigToBedGraph.

    Args:
            bigwig_path (Path): The path to the input BigWig file.
            bedgraph_path (Path): The path to write the output BedGraph file.

    Returns:
            Path: The path to the generated BigWig file."""
    if shutil.which("bigWigToBedGraph") is None:
        raise RuntimeError(
            "bigWigToBedGraph not found in PATH, cannot produce .bg outputs"
        )

    run(["bigWigToBedGraph", str(bigwig_path), str(bedgraph_path)], check=True)

    return bedgraph_path


# -- BWA parameters class --


@dataclass
class BWAParameters:
    """BWA alignment parameters.

    Contains all parameters passed to bwa aln command.
    - missing_prob_err_rate: -n parameter (default: 0.01, None means filter by any value)
    - max_gap_opens: -o parameter (default: 2, None means filter by any value)
    - seed_length: -l parameter (default: 16500, None means filter by any value)
    - all_aln: -N parameter (default: False, None means filter by any value)
    - threads: -t parameter (default: 1, None means filter by any value)
    - r_best_hits: -R parameter for bwa aln (default: 30, None means filter by any value)
    - samse_n: -n parameter for bwa samse (default: 2000000000, None means filter by any value)
    """

    missing_prob_err_rate: Optional[float] = None
    max_gap_opens: Optional[int] = None
    seed_length: Optional[int] = None
    all_aln: Optional[bool] = None
    threads: Optional[int] = None
    r_best_hits: Optional[int] = None
    samse_n: Optional[int] = None

    def __post_init__(self) -> None:
        # Apply defaults if None - only for mutable dataclass
        if self.missing_prob_err_rate is None:
            self.missing_prob_err_rate = 0.01
        if self.max_gap_opens is None:
            self.max_gap_opens = 2
        if self.seed_length is None:
            self.seed_length = 16500
        if self.all_aln is None:
            self.all_aln = False
        if self.threads is None:
            self.threads = 1
        if self.r_best_hits is None:
            self.r_best_hits = 30
        if self.samse_n is None:
            self.samse_n = 2000000000
