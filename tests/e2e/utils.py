"""Utility functions for WizardEye e2e tests.

This module provides shared utility functions used across all end-to-end test files.
"""

import random
from pathlib import Path

import pysam

from . import (
    STANDARD_BWA_HASH,
    STANDARD_KMER_LENGTH,
    STANDARD_OFFSET_STEP,
)


def get_all_fasta_sequences(fasta_path: Path) -> list[tuple[str, str]]:
    """Extract all sequences from a FASTA file.

    Args:
        fasta_path (Path): Path to the FASTA file to parse.

    Returns:
        list[tuple[str, str]]: List of (header, sequence) tuples, one per sequence.
    """
    sequences = []
    current_header = None
    current_sequence = []

    with open(fasta_path, "r") as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if current_header is not None:
                    sequences.append((current_header, "".join(current_sequence)))
                current_header = line[1:].split()[0]
                current_sequence = []
            else:
                current_sequence.append(line)

        if current_header is not None:
            sequences.append((current_header, "".join(current_sequence)))

    return sequences


def extract_random_reads_from_fasta(
    fasta_path: Path,
    num_reads: int,
    read_length: int = STANDARD_KMER_LENGTH,
    seed: int = 42,
) -> list[str]:
    """Extract random reads (subsequences) from a FASTA file.

    Extracts random subsequences of specified length from the FASTA file.
    If a sequence is shorter than read_length, the entire sequence is used.

    Args:
        fasta_path (Path): Path to the FASTA file to extract reads from.
        num_reads (int): Number of reads to extract.
        read_length (int, optional): Length of each read. Defaults to STANDARD_KMER_LENGTH.
        seed (int, optional): Random seed for reproducibility. Defaults to 42.

    Returns:
        list[str]: List of read sequences of length read_length (or shorter if source is shorter).
    """
    random.seed(seed)
    sequences = get_all_fasta_sequences(fasta_path)

    reads = []
    for _ in range(num_reads):
        # Select a random sequence (we only need the sequence, not the header)
        _, seq = random.choice(sequences)

        # If sequence is shorter than read_length, use the whole sequence
        if len(seq) <= read_length:
            reads.append(seq)
        else:
            # Extract a random subsequence of read_length
            start = random.randint(0, len(seq) - read_length)
            reads.append(seq[start : start + read_length])

    return reads


def write_fasta_file(sequences: list[tuple[str, str]], output_path: Path) -> None:
    """Write sequences to a FASTA file.

    Args:
        sequences (list[tuple[str, str]]): List of (header, sequence) tuples to write.
        output_path (Path): Path to write the FASTA file.
    """
    with open(output_path, "w") as f:
        f.writelines(
            f">read_{i}_{header}\n{seq}\n" for i, (header, seq) in enumerate(sequences)
        )


def get_track_name(
    query_fasta: Path,
    kmer_length: int = STANDARD_KMER_LENGTH,
    offset_step: int = STANDARD_OFFSET_STEP,
    bwa_hash: str = STANDARD_BWA_HASH,
) -> str:
    """Generate the WizardEye track name based on FASTA file name and parameters.

    The track name follows the convention: {query_species}_k{kmer_length}_w{offset_step}_bwa{bwa_hash}
    where query_species is the stem of the query FASTA file.

    Args:
        query_fasta (Path): Path to the query FASTA file.
        kmer_length (int, optional): K-mer length used for track generation. Defaults to STANDARD_KMER_LENGTH.
        offset_step (int, optional): Offset/step for sliding window. Defaults to STANDARD_OFFSET_STEP.
        bwa_hash (str, optional): BWA parameters hash. Defaults to STANDARD_BWA_HASH.

    Returns:
        str: The generated track name following WizardEye naming convention.
    """
    query_species = Path(query_fasta).stem
    return f"{query_species}_k{kmer_length}_w{offset_step}_bwa{bwa_hash}"


def count_mapped_reads_in_bam(bam_path: Path) -> int:
    """Count the number of mapped (aligned) reads in a BAM file.

    Only counts reads that are not marked as unmapped in the BAM file.
    This matches WizardEye filter behavior which only processes mapped reads.

    Args:
        bam_path (Path): Path to the BAM file to analyze.

    Returns:
        int: Number of mapped reads in the BAM file.
    """
    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        return sum(1 for read in bam.fetch(until_eof=True) if not read.is_unmapped)


def count_unmapped_reads_in_bam(bam_path: Path) -> int:
    """Count the number of unmapped reads in a BAM file.

    Args:
        bam_path (Path): Path to the BAM file to analyze.

    Returns:
        int: Number of unmapped reads in the BAM file.
    """
    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        return sum(1 for read in bam.fetch(until_eof=True) if read.is_unmapped)


def get_all_fasta_sequence_info(fasta_path):
    """Get all sequence names and lengths from a FASTA file.

    Args:
        fasta_path: Path to the FASTA file to parse.

    Returns:
        List of (sequence_name, sequence_length) tuples, one per sequence.
    """
    sequences = []
    seq_name = None
    seq_len = 0
    with open(fasta_path, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if seq_name is not None:
                    sequences.append((seq_name, seq_len))
                seq_name = line[1:].split()[0]
                seq_len = 0
            else:
                seq_len += len(line)
        if seq_name is not None:
            sequences.append((seq_name, seq_len))
    return sequences
