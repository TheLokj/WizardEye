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
            ["python3", "-m", "wizardeye", "database", "init", "-d", str(db_dir)],
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
                ["python3", "-m", "wizardeye", "database", "init", "-d", str(db_dir)],
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


def test_filter_duplicate_read_ids():
    """
    Test that filter processes reads with duplicate QNAME independently.

    This test creates a BAM file with reads sharing the same QNAME (read ID) but aligned
    at different positions. Some positions overlap the track mask while others don't.
    It verifies that each read is processed individually during filtering:
    - Reads overlapping masked positions are excluded
    - Reads not overlapping masked positions are kept
    - Duplicate IDs do not cause incorrect consolidation
    - Edge cases: same chr/start/stop but different chr, same chr with different positions, same chr/stop with different start
    """
    kmer_length = 35

    with tempfile.TemporaryDirectory(prefix="wizardeye_test_dup_ids_") as tmpdir:
        tmpdir = Path(tmpdir)

        # Create reference FASTA: 100 bp with first 35 as A, rest as T
        # This ensures the query (35 A's) only matches exactly positions 0-34
        ref_fasta = tmpdir / "reference.fa"
        ref_sequence = "A" * 35 + "C" * 65
        chr2_sequence = "C" * 100
        with open(ref_fasta, "w") as f:
            f.write(f">reference\n{ref_sequence}\n")
            f.write(f">chr2\n{chr2_sequence}\n")

        # Create a query FASTA of 35 A's - will match positions 0-34 of reference
        query_fasta = tmpdir / "query.fa"
        query_sequence = "A" * kmer_length
        with open(query_fasta, "w") as f:
            f.write(f">query1\n{query_sequence}\n")

        # Create SAM file with duplicate read IDs and edge cases
        sam_file = tmpdir / "input.sam"
        read_sequence = "A" * 35  # All A's to match reference A region
        read_quality = "I" * 35
        with open(sam_file, "w") as f:
            # SAM header
            f.write("@HD\tVN:1.6\tSO:unsorted\n")
            f.write("@SQ\tSN:reference\tLN:100\n")
            f.write("@SQ\tSN:chr2\tLN:100\n")
            # Reads with duplicate ID "read1":
            # - position 10: overlaps track (0-34) -> excluded
            # - position 55: in C region (50-99), does NOT overlap track -> kept
            f.write(
                f"read1\t0\treference\t10\t35\t35M\t*\t0\t0\t{read_sequence}\t{read_quality}\n"
            )
            f.write(
                f"read1\t0\treference\t55\t60\t35M\t*\t0\t0\t{read_sequence}\t{read_quality}\n"
            )
            # Reads with duplicate ID "read2":
            # - position 20: overlaps track (0-34) -> excluded
            # - position 52: in C region (50-99), does NOT overlap track -> kept
            f.write(
                f"read2\t0\treference\t20\t35\t35M\t*\t0\t0\t{read_sequence}\t{read_quality}\n"
            )
            f.write(
                f"read2\t0\treference\t52\t60\t35M\t*\t0\t0\t{read_sequence}\t{read_quality}\n"
            )
            # Same start/stop but different chr (both should be processed independently)
            # - read3 on reference:10 overlaps track -> excluded
            # - read3 on chr2:10 does NOT overlap track (chr2 has no A's, no track) -> kept
            f.write(
                f"read3\t0\treference\t10\t35\t35M\t*\t0\t0\t{read_sequence}\t{read_quality}\n"
            )
            f.write(
                f"read3\t0\tchr2\t10\t35\t35M\t*\t0\t0\t{read_sequence}\t{read_quality}\n"
            )
            # Same chr, different positions
            # - read4 at position 10 (covers 10-44) overlaps track -> excluded
            # - read4 at position 60 (covers 60-94) no overlap -> kept
            f.write(
                f"read4\t0\treference\t10\t35\t35M\t*\t0\t0\t{read_sequence}\t{read_quality}\n"
            )
            f.write(
                f"read4\t0\treference\t60\t35\t35M\t*\t0\t0\t{read_sequence}\t{read_quality}\n"
            )
            # Same chr and stop but different start
            # - read5 at position 10 with 35M (10-44) overlaps track -> excluded
            # - read5 at position 11 with 34M (11-44) overlaps track -> excluded
            f.write(
                f"read5\t0\treference\t10\t35\t35M\t*\t0\t0\t{read_sequence}\t{read_quality}\n"
            )
            f.write(
                f"read5\t0\treference\t11\t35\t34M\t*\t0\t0\t{read_sequence[:34]}\t{read_quality[:34]}\n"
            )
            # - read6 at position 65: no overlap -> filtered (reference is 100bp, 65+35=100)
            f.write(
                f"read6\t0\treference\t65\t35\t35M\t*\t0\t0\t{read_sequence}\t{read_quality}\n"
            )

        # Convert SAM to BAM and sort it
        unsorted_bam = tmpdir / "input_unsorted.bam"
        bam_file = tmpdir / "input.bam"
        subprocess.run(
            ["samtools", "view", "-b", str(sam_file), "-o", str(unsorted_bam)],
            check=True,
            capture_output=True,
            text=True,
        )

        # Sort the BAM
        subprocess.run(
            ["samtools", "sort", str(unsorted_bam), "-o", str(bam_file)],
            check=True,
            capture_output=True,
            text=True,
        )

        # Index the BAM
        subprocess.run(
            ["samtools", "index", str(bam_file)],
            check=True,
            capture_output=True,
            text=True,
        )

        # Initialize database
        db_dir = tmpdir / "database"
        db_dir.mkdir()
        subprocess.run(
            ["python3", "-m", "wizardeye", "database", "init", "-d", str(db_dir)],
            check=True,
            capture_output=True,
            text=True,
            env={**subprocess.os.environ, "PYTHONPATH": str(SRC_DIR)},
        )

        # Run wizardeye align to create a track
        # Query (35 A's) will match reference positions 0-34, creating coverage there
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
            str(1),
            "-bn",
            "0.01",
            "-bo",
            "2",
            "-bl",
            "16500",
            "-j",
            "1",
            "-bR",
            "30",
            "-bsn",
            "2000000000",
            "-d",
            str(db_path),
        ]

        result = subprocess.run(
            align_cmd,
            capture_output=True,
            text=True,
            env={**subprocess.os.environ, "PYTHONPATH": str(SRC_DIR)},
        )

        if result.returncode != 0:
            print(f"WizardEye align failed: {result.stderr}")
            raise RuntimeError(f"Alignment failed with return code {result.returncode}")

        # Find the track directory
        track_pattern = f"query_k{kmer_length}_w{1}_bwa{STANDARD_BWA_HASH}"
        track_dir = db_path / ref_stem / track_pattern

        if not track_dir.exists():
            possible_tracks = list(
                (db_path / ref_stem).glob(f"query_k{kmer_length}_w{1}_*")
            )
            if not possible_tracks:
                raise RuntimeError(f"No track directory found in {db_path / ref_stem}")
            track_dir = possible_tracks[0]

        # Verify map_all.bw exists
        map_all_bw = track_dir / "map_all.bw"
        assert map_all_bw.exists(), f"map_all.bw not found in {track_dir}"

        # Run WizardEye filter
        filtered_bam = tmpdir / "filtered.bam"
        excluded_bam = tmpdir / "excluded.bam"
        report_tsv = tmpdir / "report.tsv"

        filter_cmd = [
            "python3",
            "-m",
            "wizardeye",
            "filter",
            "-i",
            str(bam_file),
            "-r",
            "reference",
            "-k",
            str(kmer_length),
            "-w",
            str(1),
            "-bn",
            "0.01",
            "-bo",
            "2",
            "-bl",
            "16500",
            "-bR",
            "30",
            "-bsn",
            "2000000000",
            "-d",
            str(db_path),
            "--export-bam",
            "-o",
            str(filtered_bam),
            "--excluded-output",
            str(excluded_bam),
            "-p",
            "0.01",
            "--report-output",
            str(report_tsv),
            "--exclude-tracks",
            "query",
        ]

        result = subprocess.run(
            filter_cmd,
            capture_output=True,
            text=True,
            env={**subprocess.os.environ, "PYTHONPATH": str(SRC_DIR)},
        )

        if result.returncode != 0:
            print("WizardEye filter failed:")
            print(f"stdout: {result.stdout}")
            print(f"stderr: {result.stderr}")
            raise RuntimeError(f"Filter failed with return code {result.returncode}")

        # Count reads in filtered BAM
        result = subprocess.run(
            ["samtools", "view", "-c", str(filtered_bam)],
            capture_output=True,
            text=True,
            check=True,
        )
        n_filtered = int(result.stdout.strip())

        # Count reads in excluded BAM
        result = subprocess.run(
            ["samtools", "view", "-c", str(excluded_bam)],
            capture_output=True,
            text=True,
            check=True,
        )
        n_excluded = int(result.stdout.strip())

        # Verify BAM outputs exist
        assert filtered_bam.exists(), f"Filtered BAM not found: {filtered_bam}"
        assert excluded_bam.exists(), f"Excluded BAM not found: {excluded_bam}"

        # Verify BAM outputs are valid BAM files
        result = subprocess.run(
            ["samtools", "quickcheck", str(filtered_bam)],
            capture_output=True,
            text=True,
            check=True,
        )
        assert result.returncode == 0, f"Filtered BAM is corrupted: {result.stderr}"

        result = subprocess.run(
            ["samtools", "quickcheck", str(excluded_bam)],
            capture_output=True,
            text=True,
            check=True,
        )
        assert result.returncode == 0, f"Excluded BAM is corrupted: {result.stderr}"

        # Total reads should be 11:
        # - 2x read1 (10, 55)
        # - 2x read2 (20, 52)
        # - 2x read3 (reference:10, chr2:10)
        # - 2x read4 (10, 60)
        # - 2x read5 (10, 11)
        # - 1x read6 (70)
        total_reads = n_filtered + n_excluded
        assert total_reads == 11, (
            f"Expected 11 total reads after filtering, got {total_reads} "
            f"(filtered: {n_filtered}, excluded: {n_excluded}). "
            "This indicates reads were lost during filtering."
        )

        # Verify report file exists
        assert report_tsv.exists(), f"Report TSV not found: {report_tsv}"

        # Verify report format and content
        with open(report_tsv, "r") as f:
            report_lines = f.readlines()

        # Report should have header + 11 data lines (one per read)
        assert len(report_lines) == 12, (
            f"Expected 12 lines in report (header + 11 reads), got {len(report_lines)}"
        )

        # Check header
        header = report_lines[0].strip()
        assert header == "read_key	filtered_out	associated_tracks", (
            f"Unexpected report header: {header}"
        )

        # Expected read identifiers in report
        # read1@55: pysam start=54, end=89 (54+35) -> 55:89
        # read1@10: pysam start=9, end=44 (9+35) -> 10:44
        # read2@52: pysam start=51, end=86 (51+35) -> 52:86
        # read2@20: pysam start=19, end=54 (19+35) -> 20:54
        # read3@chr2:10: pysam start=9, end=44 -> 10:44
        # read3@reference:10: pysam start=9, end=44 -> 10:44
        # read4@60: pysam start=59, end=94 (59+35) -> 60:94
        # read4@10: pysam start=9, end=44 -> 10:44
        # read5@10: pysam start=9, end=44 -> 10:44
        # read5@11: pysam start=10, end=44 (10+34) -> 11:44
        # read6@65: pysam start=64, end=99 (64+35) -> 65:99
        expected_filtered_reads_in_report = {
            "read1:reference:55:89",
            "read2:reference:52:86",
            "read3:chr2:10:44",
            "read4:reference:60:94",
            "read6:reference:65:99",
        }
        expected_excluded_reads_in_report = {
            "read1:reference:10:44",
            "read2:reference:20:54",
            "read3:reference:10:44",
            "read4:reference:10:44",
            "read5:reference:10:44",
            "read5:reference:11:44",
        }

        # Parse report and verify content
        actual_filtered_in_report = set()
        actual_excluded_in_report = set()
        n_excluded_in_report = 0
        n_filtered_in_report = 0

        for i, line in enumerate(report_lines[1:]):
            parts = line.strip().split("\t")
            if len(parts) >= 2:
                read_display_id = parts[0]
                filtered_out = parts[1]
                associated_tracks = parts[2] if len(parts) >= 3 else ""
                if filtered_out == "true":
                    n_excluded_in_report += 1
                    actual_excluded_in_report.add(read_display_id)
                    # Verify excluded reads have associated tracks
                    assert associated_tracks, (
                        f"Excluded read should have associated tracks in line {i + 2}: {line.strip()}"
                    )
                elif filtered_out == "false":
                    n_filtered_in_report += 1
                    actual_filtered_in_report.add(read_display_id)
                else:
                    raise ValueError(
                        f"Unexpected filtered_out value: '{filtered_out}' in line {i + 2}: {line.strip()}"
                    )
            else:
                raise ValueError(
                    f"Line {i + 2} has insufficient columns ({len(parts)} < 2): '{line.strip()}'"
                )

        assert n_excluded_in_report == 6, (
            f"Expected 6 excluded reads in report, got {n_excluded_in_report}"
        )
        assert n_filtered_in_report == 5, (
            f"Expected 5 filtered reads in report, got {n_filtered_in_report}"
        )
        assert actual_filtered_in_report == expected_filtered_reads_in_report, (
            f"Filtered reads in report mismatch.\n"
            f"Expected: {sorted(expected_filtered_reads_in_report)}\n"
            f"Actual: {sorted(actual_filtered_in_report)}"
        )
        assert actual_excluded_in_report == expected_excluded_reads_in_report, (
            f"Excluded reads in report mismatch.\n"
            f"Expected: {sorted(expected_excluded_reads_in_report)}\n"
            f"Actual: {sorted(actual_excluded_in_report)}"
        )

        # With stringency=0, reads overlapping track positions (0-34) should be excluded
        # Track on reference covers positions 0-34 (0-based)
        # Expected results (1-based positions in SAM, converted to 0-based for overlap check):
        # - read1@10: pos 10, 35M covers 10-44 (1-based) = 9-43 (0-based) -> overlap -> excluded
        # - read1@55: pos 55, 35M covers 55-89 = 54-88 -> no overlap -> filtered
        # - read2@20: pos 20, 35M covers 20-54 = 19-53 -> overlap -> excluded
        # - read2@52: pos 52, 35M covers 52-86 = 51-85 -> no overlap -> filtered
        # - read3@reference:10: pos 10, 35M covers 10-44 = 9-43 -> overlap -> excluded
        # - read3@chr2:10: on chr2 (no track) -> filtered
        # - read4@10: pos 10, 35M covers 10-44 = 9-43 -> overlap -> excluded
        # - read4@60: pos 60, 35M covers 60-94 = 59-93 -> no overlap -> filtered
        # - read5@10: pos 10, 35M covers 10-44 = 9-43 -> overlap -> excluded
        # - read5@11: pos 11, 34M covers 11-44 = 10-43 -> overlap -> excluded
        # - read6@70: pos 70, 35M covers 70-104 = 69-103 -> no overlap -> filtered
        # Total: 5 filtered (read1@55, read2@52, read3@chr2:10, read4@60, read6@70), 6 excluded
        assert n_filtered == 5 and n_excluded == 6, (
            f"With stringency=0, expected 5 reads in filtered BAM and 6 in excluded BAM, "
            f"got {n_filtered} filtered and {n_excluded} excluded. "
            "This indicates reads with duplicate IDs or edge cases were not processed independently."
        )

        # Verify exact read content in filtered and excluded BAMs
        # Get reads from filtered BAM
        result = subprocess.run(
            ["samtools", "view", str(filtered_bam)],
            capture_output=True,
            text=True,
            check=True,
        )
        filtered_reads = []
        for line in result.stdout.strip().split("\n"):
            if not line or line.startswith("@"):
                continue
            parts = line.split("\t")
            if len(parts) >= 4:
                read_id = parts[0]
                chrom = parts[2]
                pos = int(parts[3])
                filtered_reads.append((read_id, chrom, pos))

        # Get reads from excluded BAM
        result = subprocess.run(
            ["samtools", "view", str(excluded_bam)],
            capture_output=True,
            text=True,
            check=True,
        )
        excluded_reads_list = []
        for line in result.stdout.strip().split("\n"):
            if not line or line.startswith("@"):
                continue
            parts = line.split("\t")
            if len(parts) >= 4:
                read_id = parts[0]
                chrom = parts[2]
                pos = int(parts[3])
                excluded_reads_list.append((read_id, chrom, pos))

        # Expected filtered reads: read1@55, read2@52, read3@chr2:10, read4@60, read6@65
        expected_filtered = {
            ("read1", "reference", 55),
            ("read2", "reference", 52),
            ("read3", "chr2", 10),
            ("read4", "reference", 60),
            ("read6", "reference", 65),
        }
        # Expected excluded reads: read1@10, read2@20, read3@reference:10, read4@10, read5@10, read5@11
        expected_excluded = {
            ("read1", "reference", 10),
            ("read2", "reference", 20),
            ("read3", "reference", 10),
            ("read4", "reference", 10),
            ("read5", "reference", 10),
            ("read5", "reference", 11),
        }

        actual_filtered = set(filtered_reads)
        actual_excluded = set(excluded_reads_list)

        assert actual_filtered == expected_filtered, (
            f"Filtered BAM content mismatch.\n"
            f"Expected: {sorted(expected_filtered)}\n"
            f"Actual: {sorted(actual_filtered)}"
        )
        assert actual_excluded == expected_excluded, (
            f"Excluded BAM content mismatch.\n"
            f"Expected: {sorted(expected_excluded)}\n"
            f"Actual: {sorted(actual_excluded)}"
        )

        print("\nDuplicate read IDs test results:")
        print("Total reads in input: 11")
        print(f"Reads in filtered BAM: {n_filtered}")
        print(f"Reads in excluded BAM: {n_excluded}")
        print(f"Total reads after filtering: {total_reads}")
        print(f"Report entries (excluded): {n_excluded_in_report}")
        print(f"Report entries (filtered): {n_filtered_in_report}")
