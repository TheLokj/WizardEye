"""Test script for to control WizardEye exhaustivity.

This script tests different specific cases to check if WizardEye produces expected results
"""

import random
import subprocess
import tempfile
from pathlib import Path

# Constants for test fixtures
PROJECT_ROOT = Path(__file__).parent.parent.parent
FIXTURES_DIR = Path(__file__).parent.parent / "fixtures"
SRC_DIR = PROJECT_ROOT / "src"
HG19_FA = FIXTURES_DIR / "hg19_chr1_25_1kbp.fa"
SUS_SCROFA_FA = FIXTURES_DIR / "sus_scrofa_chr1_25_1kbp.fa"
CANIS_LUPUS_FA = FIXTURES_DIR / "canis_lupus_chr1_25_1kbp.fa"
RATTUS_NORVEGICUS_FA = FIXTURES_DIR / "rattus_norvegicus_chr1_25_1kbp.fa"
SIMULATED_URSUS_BAM = FIXTURES_DIR / "ursus_1000000.uniq.L35MQ25.bam"
SCRIPT_PATH = Path(__file__).parent.parent / "generate_cross_mappability_filter_bwa.sh"

# Standard alignment parameters - reused across tests
STANDARD_KMER_LENGTH = 35
STANDARD_OFFSET_STEP = 1
STANDARD_BWA_MISSING_PROB_ERR_RATE = 0.01
STANDARD_BWA_MAX_GAP_OPENINGS = 2
STANDARD_BWA_SEED_LENGTH = 16500
STANDARD_BWA_R_BEST_HITS = 30
STANDARD_BWA_SAMSE_N = 2000000000
STANDARD_BWA_HASH = "2b5d0c37"  # MD5 hash of "0.01:2:16500:False:1:30:2000000000"
STANDARD_N_THREADS = 1
STANDARD_CHUNK_SIZE = 100000
STANDARD_CROSS_STRINGENCY = 0.99


def test_coverage_repetitive_sequence_with_bN():
    """
    Test that WizardEye refers, for a single k-mer, every alternative mappings even if there is a lot of them and despite their scores.
    """

    kmer_length = 35
    num_kmers = 5_000_000
    random.seed(4815162342)

    # Base motif - a random sequence that will be mutated
    base_motif = "".join(random.choices(["A", "C", "G", "T"], k=kmer_length))
    assert len(base_motif) == kmer_length, (
        f"Base motif length {len(base_motif)} != {kmer_length}"
    )

    with tempfile.TemporaryDirectory(prefix="wizardeye_test_repetitive_") as tmpdir:
        tmpdir = Path(tmpdir)

        # Create reference FASTA: concatenation of N k-mer sequences, each with 0-3 mutations
        # Each k-mer is exactly k bp, non-overlapping
        ref_fasta = tmpdir / "reference.fa"
        ref_sequences = []

        for i in range(num_kmers):
            # Create a variant with mutations and indels
            variant = list(base_motif)
            # Add indels: with some probability, insert or delete a base
            if random.random() < 0.5:  # 50% chance to add an indel
                # Choose position for indel (not at the very end to avoid empty sequences)
                indel_pos = random.randint(6, kmer_length - 6)
                # Randomly choose between insertion and deletion
                if random.random() < 0.5:  # insertion
                    # Insert a random base
                    possible_bases = ["A", "C", "G", "T"]
                    variant.insert(indel_pos, random.choice(possible_bases))
                else:  # deletion
                    # Delete the base at indel_pos (only if variant has enough length)
                    if len(variant) > kmer_length // 2:
                        del variant[indel_pos]
            else:
                # Randomly select 0-3 positions to mutate
                num_mutations = random.randint(0, 3)
                if num_mutations > 0:
                    mutation_positions = random.sample(
                        range(kmer_length), num_mutations
                    )
                    for pos in mutation_positions:
                        # Mutate to a different base
                        current = variant[pos]
                        possible_muts = ["A", "C", "G", "T"]
                        possible_muts.remove(current)
                        variant[pos] = random.choice(possible_muts)
            ref_sequences.append("".join(variant))

        all_contigs = set()

        with open(ref_fasta, "w") as f:
            for i, seq in enumerate(ref_sequences):
                f.write(f">contig_{i}\n{seq}\n")
                all_contigs.add(f"contig_{i}")

        # Create query FASTA: the base motif (should match all variants with <=3 mismatches)
        query_fasta = tmpdir / "query.fa"
        with open(query_fasta, "w") as f:
            f.write(f">base_kmer\n{base_motif}\n")

        # Create database directory
        db_dir = tmpdir / "database"
        db_dir.mkdir()

        # Initialize database
        subprocess.run(
            ["python3", "-m", "wizardeye", "database", "--init", "-d", str(db_dir)],
            check=True,
            capture_output=True,
            text=True,
            env={**subprocess.os.environ, "PYTHONPATH": str(SRC_DIR)},
        )

        # Run wizardeye align with -N to force exhaustive search
        # This should trigger the max_entries limit
        db_path = db_dir / "database"
        ref_stem = ref_fasta.stem

        align_cmd = [
            "python3",
            "-m",
            "wizardeye",
            "align",
            "-i",
            str(query_fasta),
            "-r",
            str(ref_fasta),
            "-k",
            str(kmer_length),
            "-w",
            str(kmer_length),
            "-bn",
            "0.01",
            "-bo",
            "2",
            "-bl",
            "16500",
            "-j",
            "1",
            "-bN",
            "-bR",
            "2000000000",
            "-bsn",
            "2000000000",
            "-d",
            str(db_path),
        ]

        result = subprocess.run(
            align_cmd,
            capture_output=False,
            text=True,
            env={**subprocess.os.environ, "PYTHONPATH": str(SRC_DIR)},
        )

        if result.returncode != 0:
            print(f"wizardeye align failed: {result.stderr}")
            raise RuntimeError(f"Alignment failed with return code {result.returncode}")

        # Find the track directory
        track_pattern = f"query_k{kmer_length}_w{kmer_length}_bwa{STANDARD_BWA_HASH}"
        track_dir = db_path / ref_stem / track_pattern

        if not track_dir.exists():
            # Try to find the actual track directory
            possible_tracks = list(
                (db_path / ref_stem).glob(f"query_k{kmer_length}_w{kmer_length}_*")
            )
            if not possible_tracks:
                raise RuntimeError(f"No track directory found in {db_path / ref_stem}")
            track_dir = possible_tracks[0]

        # Check that map_all.bw was created
        map_all_bw = track_dir / "map_all.bw"
        assert map_all_bw.exists(), f"map_all.bw not found in {track_dir}"

        # Use bigWigToBedGraph to read coverage
        bedgraph_file = tmpdir / "map_all.bg"

        try:
            subprocess.run(
                ["bigWigToBedGraph", str(map_all_bw), str(bedgraph_file)],
                check=True,
                capture_output=True,
                text=True,
            )
        except FileNotFoundError:
            print(
                "Warning: bigWigToBedGraph not found, skipping detailed coverage check"
            )
            return

        with open(bedgraph_file, "r") as f:
            lines = f.readlines()

        if not lines:
            raise RuntimeError("BedGraph file is empty")

        total_depth = 0
        for line in lines:
            if line.startswith("#") or line.startswith("track"):
                continue
            parts = line.strip().split("\t")
            if len(parts) >= 4:
                chrom, _, _, depth = (
                    parts[0],
                    int(parts[1]),
                    int(parts[2]),
                    int(parts[3]),
                )
                total_depth += depth
                if depth > 0 and chrom in all_contigs:
                    all_contigs.remove(chrom)

        print("\nRepetitive sequence test results:")
        print(f"Number of reference k-mers: {num_kmers}")
        print(f"Total depth: {total_depth}")

        assert total_depth >= num_kmers, (
            "Expected depth should be greater than the number of generated k-mers. "
        )

        assert len(all_contigs) == 0, f"{len(all_contigs)} contigs were not mapped."


def test_coverage_repetitive_sequence_with_bR_n_errs():
    """
    Test that WizardEye refers, for a single k-mer, every alternative mappings with a score up to best-1 score even if there is a lot of them.
    This test does not include indels are bwa prefers mismatches to indels.
    """

    # Reproduce the test 3 times to test 0 and 1 errors, 1 and 2 errors and 2 and 3 errors.
    errors = [0, 1, 2]

    kmer_length = 35
    num_kmers = 5_000_000
    random.seed(4815162342)

    # Base motif - a random sequence that will be mutated
    base_motif = "".join(random.choices(["A", "C", "G", "T"], k=kmer_length))
    assert len(base_motif) == kmer_length, (
        f"Base motif length {len(base_motif)} != {kmer_length}"
    )

    for error in errors:
        with tempfile.TemporaryDirectory(prefix="wizardeye_test_repetitive_") as tmpdir:
            tmpdir = Path(tmpdir)

            # Create reference FASTA: concatenation of N k-mer sequences, each with 0-3 mutations
            # Each k-mer is exactly k bp, non-overlapping
            ref_fasta = tmpdir / "reference.fa"
            ref_sequences = []

            for i in range(num_kmers):
                # Create a variant with mutations and indels
                variant = list(base_motif)
                # Only one mutation due to bwa behaviour without -bN (best_score AND best_score+1 mismatch)
                num_mutations = random.randint(error, error + 1)
                if num_mutations > 0:
                    mutation_positions = random.sample(
                        range(kmer_length), num_mutations
                    )
                    for pos in mutation_positions:
                        # Mutate to a different base
                        current = variant[pos]
                        possible_muts = ["A", "C", "G", "T"]
                        possible_muts.remove(current)
                        variant[pos] = random.choice(possible_muts)
                ref_sequences.append("".join(variant))

            all_contigs = set()

            with open(ref_fasta, "w") as f:
                for i, seq in enumerate(ref_sequences):
                    f.write(f">contig_{i}\n{seq}\n")
                    all_contigs.add(f"contig_{i}")

            # Create query FASTA: the base motif (should match all variants with <=3 mismatches)
            query_fasta = tmpdir / "query.fa"
            with open(query_fasta, "w") as f:
                f.write(f">base_kmer\n{base_motif}\n")

            # Create database directory
            db_dir = tmpdir / "database"
            db_dir.mkdir()

            # Initialize database
            subprocess.run(
                ["python3", "-m", "wizardeye", "database", "--init", "-d", str(db_dir)],
                check=True,
                capture_output=True,
                text=True,
                env={**subprocess.os.environ, "PYTHONPATH": str(SRC_DIR)},
            )

            # Run wizardeye align with -N to force exhaustive search
            # This should trigger the max_entries limit
            db_path = db_dir / "database"
            ref_stem = ref_fasta.stem

            align_cmd = [
                "python3",
                "-m",
                "wizardeye",
                "align",
                "-i",
                str(query_fasta),
                "-r",
                str(ref_fasta),
                "-k",
                str(kmer_length),
                "-w",
                str(kmer_length),
                "-bn",
                "0.01",
                "-bo",
                "2",
                "-bl",
                "16500",
                "-j",
                "1",
                "-bR",
                "2000000000",
                "-bsn",
                "2000000000",
                "-d",
                str(db_path),
            ]

            result = subprocess.run(
                align_cmd,
                capture_output=False,
                text=True,
                env={**subprocess.os.environ, "PYTHONPATH": str(SRC_DIR)},
            )

            if result.returncode != 0:
                print(f"wizardeye align failed: {result.stderr}")
                raise RuntimeError(
                    f"Alignment failed with return code {result.returncode}"
                )

            # Find the track directory
            track_pattern = (
                f"query_k{kmer_length}_w{kmer_length}_bwa{STANDARD_BWA_HASH}"
            )
            track_dir = db_path / ref_stem / track_pattern

            if not track_dir.exists():
                # Try to find the actual track directory
                possible_tracks = list(
                    (db_path / ref_stem).glob(f"query_k{kmer_length}_w{kmer_length}_*")
                )
                if not possible_tracks:
                    raise RuntimeError(
                        f"No track directory found in {db_path / ref_stem}"
                    )
                track_dir = possible_tracks[0]

            # Check that map_all.bw was created
            map_all_bw = track_dir / "map_all.bw"
            assert map_all_bw.exists(), f"map_all.bw not found in {track_dir}"

            # Use bigWigToBedGraph to read coverage
            bedgraph_file = tmpdir / "map_all.bg"

            try:
                subprocess.run(
                    ["bigWigToBedGraph", str(map_all_bw), str(bedgraph_file)],
                    check=True,
                    capture_output=True,
                    text=True,
                )
            except FileNotFoundError:
                print(
                    "Warning: bigWigToBedGraph not found, skipping detailed coverage check"
                )
                return

            with open(bedgraph_file, "r") as f:
                lines = f.readlines()

            if not lines:
                raise RuntimeError("BedGraph file is empty")

            total_depth = 0
            for line in lines:
                if line.startswith("#") or line.startswith("track"):
                    continue
                parts = line.strip().split("\t")
                if len(parts) >= 4:
                    chrom, _, _, depth = (
                        parts[0],
                        int(parts[1]),
                        int(parts[2]),
                        int(parts[3]),
                    )
                    total_depth += depth
                    if depth > 0 and chrom in all_contigs:
                        all_contigs.remove(chrom)

            print("\nRepetitive sequence test results:")
            print(f"Number of reference k-mers: {num_kmers}")
            print(f"Total depth: {total_depth}")

            assert total_depth >= num_kmers, (
                f"Expected depth should be equal or greater than the number of generated k-mers [{error}-{error + 1} errors]. "
            )

            assert len(all_contigs) == 0, (
                f"{len(all_contigs)} contigs were not mapped [{error}-{error + 1} errors]."
            )
