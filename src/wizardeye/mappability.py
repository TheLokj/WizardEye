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
import sys
import tempfile

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
    write_seq_sizes_from_fasta,
    iterate_mapping_intervals,
    iterate_unique_mapping_intervals,
    sort_bed_file,
    compute_sorted_genome_coverage,
    merge_and_sort_bams,
    BWAParameters,
    get_bwa_params_hash,
)

# --- WizardEye mappability pipeline blocks ---


def _split_fasta_into_kmers(
    input_fasta: Path,
    kmer_length: int,
    offset_step: int,
    chunk_size: int,
    output_dir: Path,
    n_threads: int = 1,
) -> List[Path]:
    """Split FASTA into k-mers, deduplicate, and chunk using seqkit pipeline.

    Runs a 3-stage pipeline via PIPE specific to the mappability workflow:
    1. seqkit sliding: generate k-mers with sliding window
    2. seqkit rmdup: deduplicate sequences
    3. seqkit split2: split into chunks of specified size

    Args:
            input_fasta (Path): Path to the input FASTA file.
            kmer_length (int): Length of k-mers to generate (-W).
            offset_step (int): Step size for sliding window (-s).
            chunk_size (int): Number of sequences per chunk for split2 (-s).
            output_dir (Path): Directory where chunk FASTA files will be written.
            n_threads (int): Number of threads for parallel operations. Default is 1.

    Returns:
            List[Path]: Sorted list of generated chunk FASTA files.

    Raises:
            RuntimeError: If no chunk FASTA files are produced.
            subprocess.CalledProcessError: If any stage of the pipeline fails.
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    p1 = p2 = p3 = None
    try:
        # 1. Sliding window
        p1 = subprocess.Popen(
            [
                "seqkit",
                "sliding",
                "-j",
                str(n_threads),
                "-W",
                str(kmer_length),
                "-s",
                str(offset_step),
                str(input_fasta),
            ],
            stdout=subprocess.PIPE,
            stderr=sys.stderr,
        )

        # 2. Deduplication
        p2 = subprocess.Popen(
            ["seqkit", "rmdup", "-j", str(n_threads), "-s"],
            stdin=p1.stdout,
            stdout=subprocess.PIPE,
            stderr=sys.stderr,
        )

        # 3. Splitting
        p3 = subprocess.Popen(
            [
                "seqkit",
                "split2",
                "-j",
                str(n_threads),
                "-s",
                str(chunk_size),
                "-O",
                str(output_dir),
                "--extension",
                ".fasta",
                "-",
            ],
            stdin=p2.stdout,
            stderr=sys.stderr,
        )

        p1.stdout.close()
        p2.stdout.close()

        p1.wait()
        p2.wait()
        p3.wait()

        for p in (p1, p2, p3):
            if p.returncode != 0:
                raise subprocess.CalledProcessError(p.returncode, p.args)

    finally:
        for p in (p1, p2, p3):
            if p is not None and p.poll() is None:
                p.kill()

    chunk_fastas = sorted(output_dir.rglob("*.fasta"))
    if not chunk_fastas:
        raise RuntimeError("No chunk FASTA files were produced by seqkit split2")

    return chunk_fastas


def _from_alignment_to_bigWig(
    bam_file: Path,
    out_dir: Path,
    kmer_length: int,
    n_threads: int = 1,
    tmp_dir: Optional[Path] = None,
) -> Tuple[Path, Path, Dict[str, int]]:
    """Convert alignment results to BigWig tracks.

    Args:
            bam_file(Path): Path to the BAM file containing all mapped k-mers.
            out_dir(Path): Directory where the output BigWig files and intermediate files will be saved.
            kmer_length(int): Length of the k-mers used for mapping, needed for coverage calculations.
            n_threads(int): Number of threads to use for parallel operations. Default is 1.
            tmp_dir(Optional[Path]): Temporary directory for sort operations. Defaults to out_dir.

    Returns:
            Tuple[Path, Path, Dict[str, int] :
                    - Path: Path to the BigWig file for all mappings (map_all).
                    - Path: Path to the BigWig file for uniquely mapping k-mers (map_uniq).
                    - Dict[str, int]: dictionary with covered base counts for both masks, e.g. {"map_all": 123456, "map_uniq": 78910}.
    """
    if tmp_dir is None:
        tmp_dir = out_dir

    cov_map_bg = out_dir / "map_all.bg"
    cov_uniq_bg = out_dir / "map_uniq.bg"
    cov_map_bw = out_dir / "map_all.bw"
    cov_uniq_bw = out_dir / "map_uniq.bw"
    seq_sizes = out_dir / "chrom.sizes"

    # Write raw intervals to temporary BED files
    n_map_bed = out_dir / "n_map.bed"
    n_uniq_bed = out_dir / "n_uniq.bed"

    with n_map_bed.open("w", encoding="utf-8") as f:
        for chrom, start, end in iterate_mapping_intervals(bam_file, kmer_length):
            f.write(f"{chrom}\t{start}\t{end}\n")

    with n_uniq_bed.open("w", encoding="utf-8") as f:
        for chrom, start, end in iterate_unique_mapping_intervals(
            bam_file, kmer_length
        ):
            f.write(f"{chrom}\t{start}\t{end}\n")

    # Sort BED files by chrom then position (required by bedtools genomecov)
    log("Sorting BED files by chromosome and position...", "I")
    sort_bed_file(n_map_bed, n_map_bed, n_threads, tmp_dir)
    sort_bed_file(n_uniq_bed, n_uniq_bed, n_threads, tmp_dir)

    # Use bedtools genomecov to compute coverage
    write_seq_sizes_from_bam(bam_file, seq_sizes)
    compute_sorted_genome_coverage(n_map_bed, seq_sizes, cov_map_bg, n_threads, tmp_dir)
    compute_sorted_genome_coverage(
        n_uniq_bed, seq_sizes, cov_uniq_bg, n_threads, tmp_dir
    )

    # Count covered bases from bedGraph files
    def count_bp_from_bedgraph(bg_path: Path) -> int:
        total = 0
        with bg_path.open("r", encoding="utf-8") as f:
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) >= 4:
                    start, end = int(parts[1]), int(parts[2])
                    total += end - start
        return total

    covered_bp = {
        "map_all": count_bp_from_bedgraph(cov_map_bg),
        "map_uniq": count_bp_from_bedgraph(cov_uniq_bg),
    }

    convert_bedgraph_to_bigwig(cov_map_bg, seq_sizes, cov_map_bw)
    convert_bedgraph_to_bigwig(cov_uniq_bg, seq_sizes, cov_uniq_bw)

    # Clean up temporary files
    for tmp_file in (n_map_bed, n_uniq_bed, cov_map_bg, cov_uniq_bg, seq_sizes):
        if tmp_file.exists():
            tmp_file.unlink()

    return cov_map_bw, cov_uniq_bw, covered_bp


def _align_with_bwa_aln(
    input_fasta: Path,
    reference: Path,
    bwa_params: BWAParameters = BWAParameters(),
) -> Path:
    """Align a FASTA file using bwa aln and return the path to the resulting BAM file containing mapped reads.

    Args:
            input_fasta(Path): Path to the input FASTA file containing k-mers to align.
            reference(Path): Path to the reference genome FASTA file that has been indexed with BWA.
            bwa(BWAParameters): BWA alignment parameters.

    Returns:
            Path: the BAM file containing the mapped reads for the input chunk.
    """
    sai_file = input_fasta.with_suffix(".sai")
    bam_file = input_fasta.with_suffix(".bam")

    if bwa_params.all_aln:
        log("bwa will compute every alternative mapping.", "I")

    bwa_cmd = [
        "bwa",
        "aln",
        "-t",
        str(bwa_params.threads),
    ]
    if bwa_params.all_aln:
        bwa_cmd.append("-N")
    bwa_cmd.extend(
        [
            "-n",
            str(bwa_params.missing_prob_err_rate),
            "-o",
            str(bwa_params.max_gap_opens),
            "-l",
            str(bwa_params.seed_length),
            "-R",
            str(bwa_params.r_best_hits),
            str(reference),
            str(input_fasta),
        ]
    )

    with sai_file.open("w", encoding="utf-8") as sai_out:
        run(bwa_cmd, check=True, stdout=sai_out)

    with bam_file.open("wb") as bam_out:
        log(
            f"bwa samse -n {str(bwa_params.samse_n)} {reference} {sai_file} {input_fasta}",
            "C",
        )
        bwa_samse = subprocess.Popen(
            [
                "bwa",
                "samse",
                "-n",
                str(bwa_params.samse_n),
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


# --- WizardEye mappability main function ---


def create_mappability_track(
    input_fasta,
    input_target,
    track_id,
    kmer_length,
    offset_step,
    chunk_size,
    n_threads,
    bwa_params: BWAParameters = BWAParameters(),
    db_root="database",
    tags: Optional[List[str]] = None,
    manual_track_id=None,
    tmp_dir_custom: Optional[str] = None,
):
    """Pipeline to create a mappability track from an input FASTA file by aligning k-mers to a reference genome and exporting BigWig files.

    Args:
            input_fasta: Path to the input FASTA file containing sequences to analyze.
            reference: Path to the reference genome FASTA file that has been indexed with BWA.
            track_id: Optional string to use as the track ID in metadata and naming. If not provided, it will be derived from the input FASTA name.
            kmer_length: Length of the k-mers to generate from the input FASTA for alignment.
            offset_step: Step size for the sliding window when generating k-mers. Default is 1 for fully overlapping k-mers.
            chunk_size: Number of k-mers to include in each chunk for parallel alignment. Must be a positive integer.
            n_threads: Number of threads to use for chunk parallel alignment.
            bwa_params: BWA alignment parameters.
            db_root: Root directory for the database where target-specific directories and track outputs will be saved. Default is "database".
            tags: Optional list of tags to associate with the track in metadata.
            manual_track_id: Optional string to use as the track ID in metadata and naming, overriding automatic derivation.
            tmp_dir_custom: Optional custom temporary directory path. If provided, overrides the default TMPDIR=/tmp location.
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
    bwa_hash = get_bwa_params_hash(bwa_params)
    track_name = f"{query_track_id}_k{kmer_length}_w{offset_step}_bwa{bwa_hash}"
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
        write_seq_sizes_from_fasta(input_target, target_seq_sizes)

    # Create temporary directory in cluster TMPDIR or /tmp, or use custom directory if provided
    # Use a unique temp directory per call to avoid conflicts when running multiple aligns in parallel
    if tmp_dir_custom:
        tmp_root = Path(tmp_dir_custom)
        tmp_root.mkdir(parents=True, exist_ok=True)
        tmp_dir = Path(tempfile.mkdtemp(dir=str(tmp_root), prefix="wizardeye_"))
    else:
        tmp_root = Path(os.environ.get("TMPDIR", "/tmp"))
        tmp_dir = Path(tempfile.mkdtemp(dir=str(tmp_root), prefix="wizardeye_"))

    try:
        # Index with BWA if needed
        bwt_file = input_target.with_suffix(input_target.suffix + ".bwt")
        if not bwt_file.exists():
            log(f"BWA index not found for {input_target}, running bwa index...", "I")
            run(["bwa", "index", str(input_target)], check=True)
        else:
            log(f"BWA index found for {input_target}, skipping index.", "I")

        log(
            f"Splitting {input_fasta} into {kmer_length}-mers, deduplicating, and chunking with seqkit pipeline...",
            "I",
        )
        chunk_dir = tmp_dir / f"{input_name}_chunks"
        chunk_fastas = _split_fasta_into_kmers(
            input_fasta=input_fasta,
            kmer_length=kmer_length,
            offset_step=offset_step,
            chunk_size=chunk_size,
            output_dir=chunk_dir,
            n_threads=n_threads,
        )

        workers = max(1, n_threads)
        if n_threads > 1:
            log(
                f"Aligning {len(chunk_fastas)} chunks with {workers} parallel workers...",
                "I",
            )

        chunk_bams: List[Path] = []
        with ThreadPoolExecutor(max_workers=workers) as executor:
            futures = [
                executor.submit(
                    _align_with_bwa_aln,
                    chunk_fasta,
                    input_target,
                    bwa_params,
                )
                for chunk_fasta in chunk_fastas
            ]
            for future in futures:
                chunk_bams.append(future.result())

        # Merge all chunk-level BAM files into one temporary BAM.
        bam_file = (
            tmp_dir
            / f"{input_name}_{kmer_length}{f'_s{offset_step}' if offset_step != 1 else ''}_{Path(input_target).stem}.bam"
        )

        log(f"Merging {len(chunk_bams)} chunk BAM files into one unique BAM...", "I")
        log(f"samtools cat {' '.join(str(path) for path in chunk_bams)}", "C")
        merge_and_sort_bams(chunk_bams, bam_file, workers)

        for tmp_bam in chunk_bams:
            if tmp_bam.exists():
                tmp_bam.unlink()

        # Export only final BigWig depth tracks.
        log("Exporting mappability BigWig tracks from temporary BAM...", "I")
        cov_map_bw, cov_uniq_bw, covered_bp = _from_alignment_to_bigWig(
            bam_file=bam_file,
            out_dir=out_dir,
            kmer_length=kmer_length,
            n_threads=n_threads,
            tmp_dir=tmp_dir,
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
                "-n": bwa_params.missing_prob_err_rate,
                "-o": bwa_params.max_gap_opens,
                "-l": bwa_params.seed_length,
                "-N": bwa_params.all_aln,
                "-t": bwa_params.threads,
                "-R": bwa_params.r_best_hits,
            },
            "bwa_samse_parameters": {
                "-n": bwa_params.samse_n,
            },
            "chunk_parameters": {
                "-j": n_threads,
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

        log(f"Track saved in: {out_dir}", "S")
        return out_dir
    finally:
        log("Cleaning up temporary files and directories...", "I")
        try:
            shutil.rmtree(tmp_dir)
        except Exception as e:
            log(f"Could not remove tmp dir {tmp_dir}: {e}", "WARN")
