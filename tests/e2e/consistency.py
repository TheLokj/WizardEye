"""Test script for different types of consistency.

This script tests:
    - consistency between the original bash script and the new Python tool,
    - consistency between the outputs of the Python tool across different runs,
    - consistency between different alternative functions of the Python tool and
      classic bioinformatic utilities.

It uses subsequences from hg19 as a reference and subsequences from sus_scrofa,
canis lupus and rattus norvegicus as potential contaminants. Three contaminants are
used to target specific bugs due to successive merge or multiple tracks.
ursus_1000000.uniq.L35MQ25.bam was made using gargammel, a Ursus genome, a Briggs model
and aligning on hg19_chr1_25_1kb.fa.
"""

from __future__ import annotations

import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import pytest

# Import shared constants
from . import (
    CANIS_LUPUS_FA,
    HG19_FA,
    RATTUS_NORVEGICUS_FA,
    SCRIPT_PATH,
    SIMULATED_URSUS_BAM,
    SRC_DIR,
    STANDARD_BWA_HASH,
    STANDARD_BWA_MAX_GAP_OPENINGS,
    STANDARD_BWA_MISSING_PROB_ERR_RATE,
    STANDARD_BWA_SEED_LENGTH,
    STANDARD_CHUNK_SIZE,
    STANDARD_CROSS_STRINGENCY,
    STANDARD_KMER_LENGTH,
    STANDARD_N_THREADS,
    STANDARD_OFFSET_STEP,
    SUS_SCROFA_FA,
)
from .utils import get_all_fasta_sequence_info


@pytest.fixture(scope="module")
def standard_database(tmp_path_factory):
    """Build a WizardEye database with standard tracks once per module.

    Generates mappability tracks for SUS_SCROFA_FA, CANIS_LUPUS_FA, and
    RATTUS_NORVEGICUS_FA against HG19_FA using standard alignment parameters.
    The database is built only once and shared across all tests that need it.

    Returns:
        Path to the database directory (db_root / "database").
    """
    db_root = tmp_path_factory.mktemp("wizardeye_standard_db")
    return generate_standard_database(
        db_root=db_root,
        reference_fasta=HG19_FA,
        kmer_length=STANDARD_KMER_LENGTH,
        offset_step=STANDARD_OFFSET_STEP,
        bwa_missing_prob_err_rate=STANDARD_BWA_MISSING_PROB_ERR_RATE,
        bwa_max_gap_opens=STANDARD_BWA_MAX_GAP_OPENINGS,
        bwa_seed_length=STANDARD_BWA_SEED_LENGTH,
        chunk_size=STANDARD_CHUNK_SIZE,
        n_threads=STANDARD_N_THREADS,
    )


def get_available_cpus() -> int:
    """Get the number of available CPUs on the system."""
    try:
        # Try to get the number of CPUs available to the current process
        cpu_count = len(os.sched_getaffinity(0))
    except (AttributeError, OSError):
        cpu_count = os.cpu_count() or 1
    return cpu_count


def adjust_thread_counts_for_test(
    desired_thread_counts: list[int], min_threads_for_test: int = 2
) -> tuple[list[int], bool, str]:
    """Adjust thread counts based on available CPUs.

    Args:
        desired_thread_counts: List of thread counts to test (e.g., [1, 2, 4, 8, 16]).
        min_threads_for_test: Minimum number of threads needed to run the test
            (default: 2, as testing with only 1 thread doesn't test parallelization).

    Returns:
        Tuple of (adjusted_thread_counts, should_skip, warning_message):
        - adjusted_thread_counts: Filtered list of thread counts that are <= available CPUs.
        - should_skip: True if test should be skipped entirely.
        - warning_message: Warning message to display (empty if no adjustment needed).
    """
    available_cpus = get_available_cpus()

    # Filter to only include thread counts that are <= available CPUs
    adjusted = [t for t in desired_thread_counts if t <= available_cpus]

    warning_message = ""
    should_skip = False

    if not adjusted:
        # No thread counts are feasible
        should_skip = True
        warning_message = (
            f"Skipping parallelization test: desired thread counts {desired_thread_counts} "
            f"all exceed available CPUs ({available_cpus})"
        )
    elif adjusted != desired_thread_counts:
        # Some thread counts were filtered out
        removed = [t for t in desired_thread_counts if t not in adjusted]
        warning_message = (
            f"Reducing thread counts for parallelization test. "
            f"Available CPUs: {available_cpus}. "
            f"Testing with: {adjusted}. "
            f"Skipped: {removed}"
        )
        # Check if we still have enough threads to make the test meaningful
        if max(adjusted) < min_threads_for_test:
            should_skip = True
            warning_message += (
                f" Skipping test entirely as max available threads ({max(adjusted)}) "
                f"is less than minimum required ({min_threads_for_test})"
            )

    return adjusted, should_skip, warning_message


def generate_standard_database(
    db_root: Path,
    reference_fasta: Path = HG19_FA,
    query_fastas: list[Path] | None = None,
    kmer_length: int = STANDARD_KMER_LENGTH,
    offset_step: int = STANDARD_OFFSET_STEP,
    bwa_missing_prob_err_rate: float = STANDARD_BWA_MISSING_PROB_ERR_RATE,
    bwa_max_gap_opens: int = STANDARD_BWA_MAX_GAP_OPENINGS,
    bwa_seed_length: int = STANDARD_BWA_SEED_LENGTH,
    chunk_size: int = STANDARD_CHUNK_SIZE,
    n_threads: int = STANDARD_N_THREADS,
    env: dict | None = None,
) -> Path:
    """Generate a standard WizardEye database with tracks for the given query species.

    This helper creates a database and generates mappability tracks for the specified
    query FASTA files against the reference FASTA. Defaults to creating tracks
    for SUS_SCROFA_FA, CANIS_LUPUS_FA, and RATTUS_NORVEGICUS_FA against HG19_FA.

    Args:
        db_root: Root where to create the database.
        reference_fasta: Reference FASTA file (default: HG19_FA).
        query_fastas: List of query FASTA files. Defaults to the 3 standard species.
        kmer_length: K-mer length for track generation.
        offset_step: Offset/step for sliding window.
        bwa_missing_prob_err_rate: BWA -n parameter.
        bwa_max_gap_opens: BWA -o parameter.
        bwa_seed_length: BWA -l parameter.
        chunk_size: Number of k-mers per chunk.
        n_threads: Number of threads for alignment.
        env: Environment variables for subprocess calls.

    Returns:
        Path to the database root directory.
    """
    if query_fastas is None:
        query_fastas = [SUS_SCROFA_FA, CANIS_LUPUS_FA, RATTUS_NORVEGICUS_FA]

    if env is None:
        env = {**subprocess.os.environ, "PYTHONPATH": str(SRC_DIR)}

    db_path = Path(db_root) / "database"

    # Initialize database if needed
    subprocess.run(
        [sys.executable, "-m", "wizardeye", "database", "init", "-d", str(db_root)],
        check=True,
        capture_output=True,
        text=True,
        encoding="utf-8",
        errors="replace",
        env=env,
    )

    # Generate track for each query FASTA
    for query_fasta in query_fastas:
        # Run WizardEye align
        cmd = [
            sys.executable,
            "-m",
            "wizardeye",
            "align",
            "-i",
            str(query_fasta),
            "-r",
            str(reference_fasta),
            "-k",
            str(kmer_length),
            "-w",
            str(offset_step),
            "-bn",
            str(bwa_missing_prob_err_rate),
            "-bo",
            str(bwa_max_gap_opens),
            "-bl",
            str(bwa_seed_length),
            "-j",
            str(n_threads),
            "-cs",
            str(chunk_size),
            "-d",
            str(db_path),
        ]

        subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            encoding="utf-8",
            errors="replace",
            check=True,
            env=env,
        )

        # Verify track directory was created
        ref_stem = Path(reference_fasta).stem
        query_stem = Path(query_fasta).stem
        track_pattern = (
            f"{query_stem}_k{kmer_length}_w{offset_step}_bwa{STANDARD_BWA_HASH}"
        )
        track_dir = db_path / ref_stem / track_pattern

        if not track_dir.exists():
            raise RuntimeError(f"Track directory not found: {track_dir}")

    return db_path


# --- Helper functions for test setup ---


def compare_bedgraph_files(file1, file2, label):
    """Compare two bedGraph files line by line and report the first difference.

    Args:
        file1: Path to the first bedGraph file.
        file2: Path to the second bedGraph file.
        label: Label used in diagnostic messages to identify the comparison.

    Returns:
        True if the files have the same number of lines and identical content,
        False otherwise.
    """
    with open(file1, "r") as f1, open(file2, "r") as f2:
        lines1 = f1.readlines()
        lines2 = f2.readlines()

    if len(lines1) != len(lines2):
        return False

    for i, (line1, line2) in enumerate(zip(lines1, lines2), 1):
        if line1 != line2:
            return False
    return True


def extract_mapped_read_keys(bam_path):
    """Extract read keys from mapped reads in a BAM file.

    The key format matches the ``read_key`` column produced by
    write_filtration_report: ``query_name:reference_name:start+1:end``.

    Args:
        bam_path: Path to the BAM file to read.

    Returns:
        Set of read key strings for all mapped reads.
    """
    import pysam

    keys = set()
    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        for read in bam.fetch(until_eof=True):
            if read.is_unmapped or read.reference_name is None:
                continue
            start = read.reference_start
            end = read.reference_end
            if start is None or end is None or end <= start:
                continue
            keys.add(f"{read.query_name}:{read.reference_name}:{start + 1}:{end}")
    return keys


def parse_filter_report(report_path):
    """Parse a filter report TSV into a dictionary keyed by read key.

    Args:
        report_path: Path to the filter report TSV file (with a header line).

    Returns:
        Dictionary mapping each read key to a
        (filtered_out, tuple_of_tracks) tuple, where filtered_out is the
        string ``"true"`` or ``"false"`` and tuple_of_tracks lists the
        associated track names.
    """
    rows = {}
    with open(report_path, "r", encoding="utf-8") as handle:
        handle.readline()  # header
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            read_key, filtered_out, tracks_str = parts[0], parts[1], parts[2]
            tracks = tuple(t for t in tracks_str.split(",") if t) if tracks_str else ()
            rows[read_key] = (filtered_out, tracks)
    return rows


# -- Consistency with the script --


def test_align_same_behavior_as_original_script():
    """Test if final alignments are equal between WizardEye align and the original script.

    Files cov_map_* and cov_uniq_* should be equal to map_all and map_uniq.
    Tests all query FASTAs defined at the top of this file against HG19_FA as target.
    """

    # List of query FASTAs to test against HG19_FA as target
    query_fastas = [SUS_SCROFA_FA, CANIS_LUPUS_FA, RATTUS_NORVEGICUS_FA]

    for query_fa in query_fastas:
        # Create temporary directories for outputs
        with (
            tempfile.TemporaryDirectory(
                prefix="wizardeye_test_original_"
            ) as original_tmpdir,
            tempfile.TemporaryDirectory(
                prefix="wizardeye_test_wizardeye_"
            ) as wizardeye_tmpdir,
        ):
            original_output = Path(original_tmpdir) / "output"
            wizardeye_db = Path(wizardeye_tmpdir) / "database"
            wizardeye_db.mkdir(parents=True)

            # Prepare input directory for original script - it needs all FASTAs in one directory
            # The script iterates over *.fa/*.fasta in the directory and skips the target
            original_input_dir = Path(original_tmpdir) / "input"
            original_input_dir.mkdir(parents=True)

            # Copy target and query FASTAs to input directory
            hg19_copy = original_input_dir / HG19_FA.name
            query_copy = original_input_dir / query_fa.name
            shutil.copy(HG19_FA, hg19_copy)
            shutil.copy(query_fa, query_copy)

            # Run original script with hg19 as target, processing all FASTAs in directory
            # The script will process query_fa against hg19 (skipping hg19 itself)
            original_cmd = [
                "bash",
                str(SCRIPT_PATH),
                f"-i={original_input_dir}",
                f"-t={hg19_copy}",
                f"-k={STANDARD_KMER_LENGTH}",
                f"-w={STANDARD_OFFSET_STEP}",
                f"-bn={STANDARD_BWA_MISSING_PROB_ERR_RATE}",
                f"-bo={STANDARD_BWA_MAX_GAP_OPENINGS}",
                f"-bl={STANDARD_BWA_SEED_LENGTH}",
                f"-cs={STANDARD_CHUNK_SIZE}",
                f"-j={STANDARD_N_THREADS}",
                f"-s={STANDARD_CROSS_STRINGENCY}",
                f"-o={original_output}/",
            ]

            subprocess.run(original_cmd, capture_output=True, text=True, check=True)

            # Verify original script created expected outputs in output/tracks/
            tracks_dir = original_output / "tracks"
            assert tracks_dir.exists(), (
                f"Original script tracks directory not found: {tracks_dir}"
            )

            query_basename = query_fa.stem
            original_cov_map_bw = tracks_dir / f"cov_map_{query_basename}.bw"
            original_cov_uniq_bw = tracks_dir / f"cov_uniq_{query_basename}.bw"

            assert original_cov_map_bw.exists(), (
                f"Original script did not create {original_cov_map_bw}"
            )
            assert original_cov_uniq_bw.exists(), (
                f"Original script did not create {original_cov_uniq_bw}"
            )

            # Initialize WizardEye database
            subprocess.run(
                [
                    sys.executable,
                    "-m",
                    "wizardeye",
                    "database",
                    "init",
                    "-d",
                    str(wizardeye_tmpdir),
                ],
                check=True,
                capture_output=True,
                text=True,
                env={**subprocess.os.environ, "PYTHONPATH": str(SRC_DIR)},
            )

            # Run WizardEye align with query_fa against hg19 using same parameters
            wizardeye_cmd = [
                sys.executable,
                "-m",
                "wizardeye",
                "align",
                "-i",
                str(query_fa),
                "-r",
                str(HG19_FA),
                "-k",
                str(STANDARD_KMER_LENGTH),
                "-w",
                str(STANDARD_OFFSET_STEP),
                "-bn",
                str(STANDARD_BWA_MISSING_PROB_ERR_RATE),
                "-bo",
                str(STANDARD_BWA_MAX_GAP_OPENINGS),
                "-bl",
                str(STANDARD_BWA_SEED_LENGTH),
                "-bR",
                str(30),
                "-bsn",
                str(2000000000),
                "-j",
                str(STANDARD_N_THREADS),
                "-cs",
                str(STANDARD_CHUNK_SIZE),
                "-d",
                str(wizardeye_db),
            ]

            subprocess.run(
                wizardeye_cmd,
                capture_output=True,
                text=True,
                check=True,
                env={**subprocess.os.environ, "PYTHONPATH": str(SRC_DIR)},
            )

            # Locate WizardEye output - track dir name uses float :g format for bwa params
            hg19_stem = HG19_FA.stem
            query_stem = query_fa.stem

            track_pattern = f"{query_stem}_k{STANDARD_KMER_LENGTH}_w{STANDARD_OFFSET_STEP}_bwa{STANDARD_BWA_HASH}"
            track_dir = wizardeye_db / hg19_stem / track_pattern

            assert track_dir.exists(), (
                f"wizardeye track directory not found: {track_dir}"
            )

            wizardeye_map_all_bw = track_dir / "map_all.bw"
            wizardeye_map_uniq_bw = track_dir / "map_uniq.bw"

            assert wizardeye_map_all_bw.exists(), (
                f"wizardeye did not create {wizardeye_map_all_bw}"
            )
            assert wizardeye_map_uniq_bw.exists(), (
                f"wizardeye did not create {wizardeye_map_uniq_bw}"
            )

            # Compare outputs by converting BigWig to bedGraph and comparing line by line
            with tempfile.TemporaryDirectory(
                prefix="wizardeye_compare_"
            ) as compare_dir:
                compare_path = Path(compare_dir)

                # Convert both outputs to bedGraph
                orig_map_all_bg = compare_path / "original_map_all.bg"
                orig_map_uniq_bg = compare_path / "original_map_uniq.bg"
                we_map_all_bg = compare_path / "wizardeye_map_all.bg"
                we_map_uniq_bg = compare_path / "wizardeye_map_uniq.bg"

                subprocess.run(
                    [
                        "bigWigToBedGraph",
                        str(original_cov_map_bw),
                        str(orig_map_all_bg),
                    ],
                    check=True,
                )
                subprocess.run(
                    [
                        "bigWigToBedGraph",
                        str(original_cov_uniq_bw),
                        str(orig_map_uniq_bg),
                    ],
                    check=True,
                )
                subprocess.run(
                    ["bigWigToBedGraph", str(wizardeye_map_all_bw), str(we_map_all_bg)],
                    check=True,
                )
                subprocess.run(
                    [
                        "bigWigToBedGraph",
                        str(wizardeye_map_uniq_bw),
                        str(we_map_uniq_bg),
                    ],
                    check=True,
                )

                # Compare: cov_map_* (original) vs map_all (WizardEye) - both represent all mappings
                assert compare_bedgraph_files(
                    orig_map_all_bg, we_map_all_bg, f"map_all ({query_fa.name})"
                ), (
                    f"cov_map_* (original) and map_all.bw (wizardeye) differ for {query_fa.name}"
                )

                # Compare: cov_uniq_* (original) vs map_uniq (WizardEye) - both represent unique mappings
                assert compare_bedgraph_files(
                    orig_map_uniq_bg, we_map_uniq_bg, f"map_uniq ({query_fa.name})"
                ), (
                    f"cov_uniq_* (original) and map_uniq.bw (wizardeye) differ for {query_fa.name}"
                )


def test_export_same_behavior_as_original_script():
    """Test if the exported mask is the same between WizardEye export and the original script.

    Uses common parameters for comparison with the original
    generate_cross_mappability_filter_bwa.sh script.
    """

    # All three query FASTAs to include in the final mask
    query_fastas = [SUS_SCROFA_FA, CANIS_LUPUS_FA, RATTUS_NORVEGICUS_FA]
    query_stems = [f.stem for f in query_fastas]

    with (
        tempfile.TemporaryDirectory(
            prefix="wizardeye_export_original_"
        ) as original_tmpdir,
        tempfile.TemporaryDirectory(
            prefix="wizardeye_export_wizardeye_"
        ) as wizardeye_tmpdir,
    ):
        original_output = Path(original_tmpdir) / "output"
        wizardeye_db = Path(wizardeye_tmpdir) / "database"
        wizardeye_db.mkdir(parents=True)

        # Prepare input directory for original script - copy target and ALL query FASTAs
        original_input_dir = Path(original_tmpdir) / "input"
        original_input_dir.mkdir(parents=True)

        # Copy target FASTA
        hg19_copy = original_input_dir / HG19_FA.name
        shutil.copy(HG19_FA, hg19_copy)

        # Copy all query FASTAs
        for query_fa in query_fastas:
            query_copy = original_input_dir / query_fa.name
            shutil.copy(query_fa, query_copy)

        # Run original script ONCE with hg19 as target - processes all FASTAs in directory
        original_cmd = [
            "bash",
            str(SCRIPT_PATH),
            f"-i={original_input_dir}",
            f"-t={hg19_copy}",
            f"-k={STANDARD_KMER_LENGTH}",
            f"-w={STANDARD_OFFSET_STEP}",
            f"-bn={STANDARD_BWA_MISSING_PROB_ERR_RATE}",
            f"-bo={STANDARD_BWA_MAX_GAP_OPENINGS}",
            f"-bl={STANDARD_BWA_SEED_LENGTH}",
            f"-cs={STANDARD_CHUNK_SIZE}",
            f"-j={STANDARD_N_THREADS}",
            f"-s={STANDARD_CROSS_STRINGENCY}",
            f"-o={original_output}/",
        ]

        subprocess.run(
            original_cmd,
            capture_output=True,
            text=True,
            encoding="utf-8",
            errors="replace",
            check=True,
        )

        # Check that original script created the export mask with all 3 species
        bn_str = f"{float(STANDARD_CROSS_STRINGENCY):g}"
        original_mask = (
            original_output / f"all_overlaps_mask_{bn_str}_k{STANDARD_KMER_LENGTH}.bed"
        )
        assert original_mask.exists(), (
            f"Original script did not create export mask: {original_mask}"
        )

        # Initialize WizardEye database
        subprocess.run(
            [
                sys.executable,
                "-m",
                "wizardeye",
                "database",
                "init",
                "-d",
                str(wizardeye_tmpdir),
            ],
            check=True,
            capture_output=True,
            text=True,
            encoding="utf-8",
            errors="replace",
            env={**subprocess.os.environ, "PYTHONPATH": str(SRC_DIR)},
        )

        # Run WizardEye align for EACH query FASTA against hg19
        for query_fa in query_fastas:
            wizardeye_cmd = [
                sys.executable,
                "-m",
                "wizardeye",
                "align",
                "-i",
                str(query_fa),
                "-r",
                str(HG19_FA),
                "-k",
                str(STANDARD_KMER_LENGTH),
                "-w",
                str(STANDARD_OFFSET_STEP),
                "-bn",
                str(STANDARD_BWA_MISSING_PROB_ERR_RATE),
                "-bo",
                str(STANDARD_BWA_MAX_GAP_OPENINGS),
                "-bl",
                str(STANDARD_BWA_SEED_LENGTH),
                "-j",
                str(STANDARD_N_THREADS),
                "-cs",
                str(STANDARD_CHUNK_SIZE),
                "-d",
                str(wizardeye_db),
            ]

            subprocess.run(
                wizardeye_cmd,
                capture_output=True,
                text=True,
                encoding="utf-8",
                errors="replace",
                check=True,
                env={**subprocess.os.environ, "PYTHONPATH": str(SRC_DIR)},
            )

        # Run WizardEye export ONCE with ALL 3 species in exclude-tracks
        hg19_stem = HG19_FA.stem
        exclude_tracks_str = ",".join(query_stems)

        # Use a specific output file for WizardEye export
        wizardeye_mask = wizardeye_db / "export" / "wizardeye_mask_all_species.bed"
        wizardeye_mask.parent.mkdir(parents=True, exist_ok=True)

        wizardeye_export_cmd = [
            sys.executable,
            "-m",
            "wizardeye",
            "export",
            "-r",
            str(hg19_stem),
            "--exclude-tracks",
            exclude_tracks_str,
            "-k",
            str(STANDARD_KMER_LENGTH),
            "-w",
            str(STANDARD_OFFSET_STEP),
            "-bn",
            str(STANDARD_BWA_MISSING_PROB_ERR_RATE),
            "-bo",
            str(STANDARD_BWA_MAX_GAP_OPENINGS),
            "-bl",
            str(STANDARD_BWA_SEED_LENGTH),
            "-p",
            str(STANDARD_CROSS_STRINGENCY),
            "-o",
            str(wizardeye_mask),
            "-d",
            str(wizardeye_db),
            "--only-unique",
        ]

        subprocess.run(
            wizardeye_export_cmd,
            capture_output=True,
            text=True,
            encoding="utf-8",
            errors="replace",
            check=True,
            env={**subprocess.os.environ, "PYTHONPATH": str(SRC_DIR)},
        )

        assert wizardeye_mask.exists(), (
            f"wizardeye export did not create mask: {wizardeye_mask}"
        )

        # WizardEye export produces per-track intervals, need to merge them globally
        # to match the original script behavior which does a global merge
        wizardeye_mask_merged = (
            wizardeye_db / "export" / "wizardeye_mask_all_species_merged.bed"
        )
        with open(wizardeye_mask_merged, "w") as out_file:
            subprocess.run(
                ["bedtools", "merge", "-i", str(wizardeye_mask)],
                stdout=out_file,
                check=True,
                text=True,
            )

        # Compare the two mask files (only first 3 columns - chrom, start, end)
        # The 4th column differs: script has empty, WizardEye has track names
        with tempfile.TemporaryDirectory(
            prefix="wizardeye_export_compare_"
        ) as compare_dir:
            compare_path = Path(compare_dir)

            # Sort both files for comparison (mask order may differ)
            orig_sorted = compare_path / "original_mask_sorted.bed"
            we_sorted = compare_path / "wizardeye_mask_sorted.bed"

            subprocess.run(
                ["sort", str(original_mask), "-o", str(orig_sorted)], check=True
            )
            subprocess.run(
                ["sort", str(wizardeye_mask_merged), "-o", str(we_sorted)], check=True
            )

            # Compare sorted files by first 3 columns only
            with open(orig_sorted, "r") as f1, open(we_sorted, "r") as f2:
                orig_lines = f1.readlines()
                we_lines = f2.readlines()

            if len(orig_lines) != len(we_lines):
                raise AssertionError("Export masks differ in line count")

            for i, (line1, line2) in enumerate(zip(orig_lines, we_lines), 1):
                # Split by tab and compare only first 3 columns
                parts1 = line1.rstrip().split("\t")
                parts2 = line2.rstrip().split("\t")
                if len(parts1) < 3 or len(parts2) < 3:
                    raise AssertionError(
                        f"Export masks have invalid format at line {i}:\n"
                        f"  Original: {line1.rstrip()}\n"
                        f"  WizardEye: {line2.rstrip()}"
                    )
                if parts1[:3] != parts2[:3]:
                    raise AssertionError(
                        f"Export masks differ at line {i} (first 3 cols):\n"
                        f"  Original: {parts1[:3]}\n"
                        f"  WizardEye: {parts2[:3]}"
                    )


# -- Consistency between executions --


def test_consistency_align_same_launch():
    """Test if the results are the same across 10 align rounds with the same parameters."""

    # List of query FASTAs to test against HG19_FA as target
    query_fastas = [SUS_SCROFA_FA, CANIS_LUPUS_FA, RATTUS_NORVEGICUS_FA]

    # Number of repetitions
    n_repetitions = 10

    for query_fa in query_fastas:
        # Use a single persistent temp directory for all repetitions
        with tempfile.TemporaryDirectory(
            prefix=f"wizardeye_consistency_align_{query_fa.stem}_"
        ) as base_tmpdir:
            base_path = Path(base_tmpdir)

            # Dictionary to store BigWig paths per repetition
            bw_files = {}

            for rep in range(n_repetitions):
                with tempfile.TemporaryDirectory(prefix=f"rep{rep}_") as tmpdir:
                    wizardeye_db = Path(tmpdir) / "database"
                    wizardeye_db.mkdir(parents=True)

                    # Initialize WizardEye database
                    subprocess.run(
                        [
                            sys.executable,
                            "-m",
                            "wizardeye",
                            "database",
                            "init",
                            "-d",
                            str(tmpdir),
                        ],
                        check=True,
                        capture_output=True,
                        text=True,
                        encoding="utf-8",
                        errors="replace",
                        env={**subprocess.os.environ, "PYTHONPATH": str(SRC_DIR)},
                    )

                    # Run WizardEye align
                    wizardeye_cmd = [
                        sys.executable,
                        "-m",
                        "wizardeye",
                        "align",
                        "-i",
                        str(query_fa),
                        "-r",
                        str(HG19_FA),
                        "-k",
                        str(STANDARD_KMER_LENGTH),
                        "-w",
                        str(STANDARD_OFFSET_STEP),
                        "-bn",
                        str(STANDARD_BWA_MISSING_PROB_ERR_RATE),
                        "-bo",
                        str(STANDARD_BWA_MAX_GAP_OPENINGS),
                        "-bl",
                        str(STANDARD_BWA_SEED_LENGTH),
                        "-j",
                        str(STANDARD_N_THREADS),
                        "-cs",
                        str(STANDARD_CHUNK_SIZE),
                        "-d",
                        str(wizardeye_db),
                    ]

                    subprocess.run(
                        wizardeye_cmd,
                        capture_output=True,
                        text=True,
                        encoding="utf-8",
                        errors="replace",
                        check=True,
                        env={**subprocess.os.environ, "PYTHONPATH": str(SRC_DIR)},
                    )

                    # Locate WizardEye output
                    hg19_stem = HG19_FA.stem
                    query_stem = query_fa.stem
                    track_pattern = f"{query_stem}_k{STANDARD_KMER_LENGTH}_w{STANDARD_OFFSET_STEP}_bwa{STANDARD_BWA_HASH}"
                    track_dir = wizardeye_db / hg19_stem / track_pattern

                    assert track_dir.exists(), (
                        f"wizardeye track directory not found: {track_dir}"
                    )

                    wizardeye_map_all_bw = track_dir / "map_all.bw"
                    wizardeye_map_uniq_bw = track_dir / "map_uniq.bw"

                    assert wizardeye_map_all_bw.exists(), (
                        f"wizardeye did not create {wizardeye_map_all_bw}"
                    )
                    assert wizardeye_map_uniq_bw.exists(), (
                        f"wizardeye did not create {wizardeye_map_uniq_bw}"
                    )

                    # Copy BigWig files to persistent base directory before tmpdir is cleaned up
                    dest_dir = base_path / f"rep{rep}"
                    dest_dir.mkdir(parents=True, exist_ok=True)
                    shutil.copy(wizardeye_map_all_bw, dest_dir / "map_all.bw")
                    shutil.copy(wizardeye_map_uniq_bw, dest_dir / "map_uniq.bw")

                    bw_files[rep] = (dest_dir / "map_all.bw", dest_dir / "map_uniq.bw")

            # Compare all outputs with the first one
            first_map_all, first_map_uniq = bw_files[0]

            for rep in range(1, n_repetitions):
                current_map_all, current_map_uniq = bw_files[rep]

                with tempfile.TemporaryDirectory(
                    prefix=f"wizardeye_compare_consistency_{query_fa.stem}_rep{rep}_"
                ) as compare_dir:
                    compare_path = Path(compare_dir)

                    # Convert to bedGraph for comparison
                    first_all_bg = compare_path / "first_map_all.bg"
                    first_uniq_bg = compare_path / "first_map_uniq.bg"
                    current_all_bg = compare_path / "current_map_all.bg"
                    current_uniq_bg = compare_path / "current_map_uniq.bg"

                    subprocess.run(
                        ["bigWigToBedGraph", str(first_map_all), str(first_all_bg)],
                        check=True,
                    )
                    subprocess.run(
                        ["bigWigToBedGraph", str(first_map_uniq), str(first_uniq_bg)],
                        check=True,
                    )
                    subprocess.run(
                        ["bigWigToBedGraph", str(current_map_all), str(current_all_bg)],
                        check=True,
                    )
                    subprocess.run(
                        [
                            "bigWigToBedGraph",
                            str(current_map_uniq),
                            str(current_uniq_bg),
                        ],
                        check=True,
                    )

                    # Compare map_all
                    assert compare_bedgraph_files(
                        first_all_bg,
                        current_all_bg,
                        f"map_all rep0 vs rep{rep} ({query_fa.name})",
                    ), (
                        f"map_all.bw differs between rep 0 and rep {rep} for {query_fa.name}"
                    )

                    # Compare map_uniq
                    assert compare_bedgraph_files(
                        first_uniq_bg,
                        current_uniq_bg,
                        f"map_uniq rep0 vs rep{rep} ({query_fa.name})",
                    ), (
                        f"map_uniq.bw differs between rep 0 and rep {rep} for {query_fa.name}"
                    )


def test_consistency_align_parallelisation():
    """Test if the results are the same with different numbers of cores (1, 2, 4, 8, 16).

    It uses a smaller chunk size to allow the use of different numbers of threads.
    """

    # Use a smaller chunk size for parallelisation testing
    chunk_size = 5000

    # Different thread counts to test
    desired_thread_counts = [1, 2, 4, 8, 16]
    thread_counts, should_skip, warning_msg = adjust_thread_counts_for_test(
        desired_thread_counts, min_threads_for_test=2
    )

    if should_skip:
        pytest.skip(warning_msg)

    # List of query FASTAs to test against HG19_FA as target - use single species for speed
    query_fastas = [SUS_SCROFA_FA, CANIS_LUPUS_FA, RATTUS_NORVEGICUS_FA]

    for query_fa in query_fastas:
        # Use a single persistent temp directory for all thread counts
        with tempfile.TemporaryDirectory(
            prefix=f"wizardeye_consistency_parallel_{query_fa.stem}_"
        ) as base_tmpdir:
            base_path = Path(base_tmpdir)

            # Dictionary to store BigWig paths per thread count
            bw_files = {}

            for n_threads in thread_counts:
                with tempfile.TemporaryDirectory(
                    prefix=f"threads{n_threads}_"
                ) as tmpdir:
                    wizardeye_db = Path(tmpdir) / "database"
                    wizardeye_db.mkdir(parents=True)

                    # Initialize WizardEye database
                    subprocess.run(
                        [
                            sys.executable,
                            "-m",
                            "wizardeye",
                            "database",
                            "init",
                            "-d",
                            str(tmpdir),
                        ],
                        check=True,
                        capture_output=True,
                        text=True,
                        encoding="utf-8",
                        errors="replace",
                        env={**subprocess.os.environ, "PYTHONPATH": str(SRC_DIR)},
                    )

                    # Run WizardEye align with this thread count
                    wizardeye_cmd = [
                        sys.executable,
                        "-m",
                        "wizardeye",
                        "align",
                        "-i",
                        str(query_fa),
                        "-r",
                        str(HG19_FA),
                        "-k",
                        str(STANDARD_KMER_LENGTH),
                        "-w",
                        str(STANDARD_OFFSET_STEP),
                        "-bn",
                        str(STANDARD_BWA_MISSING_PROB_ERR_RATE),
                        "-bo",
                        str(STANDARD_BWA_MAX_GAP_OPENINGS),
                        "-bl",
                        str(STANDARD_BWA_SEED_LENGTH),
                        "-j",
                        str(n_threads),
                        "-cs",
                        str(chunk_size),
                        "-d",
                        str(wizardeye_db),
                    ]

                    subprocess.run(
                        wizardeye_cmd,
                        capture_output=True,
                        text=True,
                        encoding="utf-8",
                        errors="replace",
                        check=True,
                        env={**subprocess.os.environ, "PYTHONPATH": str(SRC_DIR)},
                    )

                    # Locate WizardEye output
                    hg19_stem = HG19_FA.stem
                    query_stem = query_fa.stem
                    track_pattern = f"{query_stem}_k{STANDARD_KMER_LENGTH}_w{STANDARD_OFFSET_STEP}_bwa{STANDARD_BWA_HASH}"
                    track_dir = wizardeye_db / hg19_stem / track_pattern

                    assert track_dir.exists(), (
                        f"wizardeye track directory not found: {track_dir}"
                    )

                    wizardeye_map_all_bw = track_dir / "map_all.bw"
                    wizardeye_map_uniq_bw = track_dir / "map_uniq.bw"

                    assert wizardeye_map_all_bw.exists(), (
                        f"wizardeye did not create {wizardeye_map_all_bw}"
                    )
                    assert wizardeye_map_uniq_bw.exists(), (
                        f"wizardeye did not create {wizardeye_map_uniq_bw}"
                    )

                    # Copy BigWig files to persistent base directory before tmpdir is cleaned up
                    dest_dir = base_path / f"threads{n_threads}"
                    dest_dir.mkdir(parents=True, exist_ok=True)
                    shutil.copy(wizardeye_map_all_bw, dest_dir / "map_all.bw")
                    shutil.copy(wizardeye_map_uniq_bw, dest_dir / "map_uniq.bw")

                    bw_files[n_threads] = (
                        dest_dir / "map_all.bw",
                        dest_dir / "map_uniq.bw",
                    )

            # Compare all outputs with the first one (1 thread)
            first_threads = thread_counts[0]
            first_map_all, first_map_uniq = bw_files[first_threads]

            for n_threads in thread_counts[1:]:
                current_map_all, current_map_uniq = bw_files[n_threads]

                with tempfile.TemporaryDirectory(
                    prefix=f"wizardeye_compare_parallel_threads{n_threads}_"
                ) as compare_dir:
                    compare_path = Path(compare_dir)

                    # Convert to bedGraph for comparison
                    first_all_bg = compare_path / "first_map_all.bg"
                    first_uniq_bg = compare_path / "first_map_uniq.bg"
                    current_all_bg = compare_path / "current_map_all.bg"
                    current_uniq_bg = compare_path / "current_map_uniq.bg"

                    subprocess.run(
                        ["bigWigToBedGraph", str(first_map_all), str(first_all_bg)],
                        check=True,
                    )
                    subprocess.run(
                        ["bigWigToBedGraph", str(first_map_uniq), str(first_uniq_bg)],
                        check=True,
                    )
                    subprocess.run(
                        ["bigWigToBedGraph", str(current_map_all), str(current_all_bg)],
                        check=True,
                    )
                    subprocess.run(
                        [
                            "bigWigToBedGraph",
                            str(current_map_uniq),
                            str(current_uniq_bg),
                        ],
                        check=True,
                    )

                    # Compare map_all
                    assert compare_bedgraph_files(
                        first_all_bg,
                        current_all_bg,
                        f"map_all threads{first_threads} vs threads{n_threads} ({query_fa.name})",
                    ), (
                        f"map_all.bw differs between threads {first_threads} and {n_threads} for {query_fa.name}"
                    )

                    # Compare map_uniq
                    assert compare_bedgraph_files(
                        first_uniq_bg,
                        current_uniq_bg,
                        f"map_uniq threads{first_threads} vs threads{n_threads} ({query_fa.name})",
                    ), (
                        f"map_uniq.bw differs between threads {first_threads} and {n_threads} for {query_fa.name}"
                    )


def test_consistency_count_same_launch(standard_database):
    """Test if the results are the same across 10 count rounds with the same parameters."""

    # Number of repetitions
    n_repetitions = 10

    # Use a single persistent temp directory for all repetitions
    with tempfile.TemporaryDirectory(
        prefix="wizardeye_consistency_count_same_launch_"
    ) as base_tmpdir:
        base_path = Path(base_tmpdir)

        # Dictionary to store report paths per repetition
        report_files = {}

        # Use the shared module-scoped database
        wizardeye_db = standard_database

        # Get track names for exclude-tracks parameter
        query_stems = [
            SUS_SCROFA_FA.stem,
            CANIS_LUPUS_FA.stem,
            RATTUS_NORVEGICUS_FA.stem,
        ]
        exclude_tracks_str = ",".join(query_stems)
        hg19_stem = HG19_FA.stem

        for rep in range(n_repetitions):
            with tempfile.TemporaryDirectory(prefix=f"count_rep{rep}_") as tmpdir:
                output_dir = Path(tmpdir) / "output"
                output_dir.mkdir(parents=True)

                # Run WizardEye count with this repetition
                report_path = output_dir / "count_report.tsv"
                wizardeye_cmd = [
                    sys.executable,
                    "-m",
                    "wizardeye",
                    "count",
                    "-i",
                    str(SIMULATED_URSUS_BAM),
                    "-r",
                    str(hg19_stem),
                    "-k",
                    str(STANDARD_KMER_LENGTH),
                    "-w",
                    str(STANDARD_OFFSET_STEP),
                    "-bn",
                    str(STANDARD_BWA_MISSING_PROB_ERR_RATE),
                    "-bo",
                    str(STANDARD_BWA_MAX_GAP_OPENINGS),
                    "-bl",
                    str(STANDARD_BWA_SEED_LENGTH),
                    "-m",
                    "mean",
                    "--exclude-tracks",
                    exclude_tracks_str,
                    "--report-output",
                    str(report_path),
                    "-d",
                    str(wizardeye_db),
                ]

                subprocess.run(
                    wizardeye_cmd,
                    capture_output=True,
                    text=True,
                    encoding="utf-8",
                    errors="replace",
                    check=True,
                    env={**subprocess.os.environ, "PYTHONPATH": str(SRC_DIR)},
                )

                assert report_path.exists(), (
                    f"wizardeye count did not create report: {report_path}"
                )

                # Copy report to persistent base directory before tmpdir is cleaned up
                dest_dir = base_path / f"rep{rep}"
                dest_dir.mkdir(parents=True, exist_ok=True)
                dest_report = dest_dir / "count_report.tsv"
                shutil.copy(report_path, dest_report)

                report_files[rep] = dest_report

        # Compare all outputs with the first one (rep 0)
        first_report = report_files[0]

        for rep in range(1, n_repetitions):
            current_report = report_files[rep]

            # Compare report files line by line
            with open(first_report, "r") as f1, open(current_report, "r") as f2:
                first_lines = f1.readlines()
                current_lines = f2.readlines()

            if len(first_lines) != len(current_lines):
                raise AssertionError(
                    f"Count reports differ in line count between rep 0 and rep {rep}: "
                    f"{len(first_lines)} vs {len(current_lines)}"
                )

            for i, (line1, line2) in enumerate(zip(first_lines, current_lines), 1):
                if line1 != line2:
                    raise AssertionError(
                        f"Count reports differ at line {i} between rep 0 and rep {rep}:\n"
                        f"  Rep 0: {line1.rstrip()}\n"
                        f"  Rep {rep}: {line2.rstrip()}"
                    )


def test_consistency_count_parallelisation(standard_database):
    """Test if the results are the same with different numbers of cores (1, 2, 4, 8, 16)."""

    # Different thread counts to test for count parallelisation
    desired_thread_counts = [1, 2, 4, 8, 16]
    thread_counts, should_skip, warning_msg = adjust_thread_counts_for_test(
        desired_thread_counts, min_threads_for_test=2
    )

    if should_skip:
        pytest.skip(warning_msg)

    # Use a single persistent temp directory for all thread counts
    with tempfile.TemporaryDirectory(
        prefix="wizardeye_consistency_count_parallel_"
    ) as base_tmpdir:
        base_path = Path(base_tmpdir)

        # Dictionary to store report paths per thread count
        report_files = {}

        # Use the shared module-scoped database
        wizardeye_db = standard_database

        # Get track names for exclude-tracks parameter
        query_stems = [
            SUS_SCROFA_FA.stem,
            CANIS_LUPUS_FA.stem,
            RATTUS_NORVEGICUS_FA.stem,
        ]
        exclude_tracks_str = ",".join(query_stems)
        hg19_stem = HG19_FA.stem

        for n_threads in thread_counts:
            with tempfile.TemporaryDirectory(
                prefix=f"count_threads{n_threads}_"
            ) as tmpdir:
                output_dir = Path(tmpdir) / "output"
                output_dir.mkdir(parents=True)

                # Run WizardEye count with this thread count
                report_path = output_dir / "count_report.tsv"
                wizardeye_cmd = [
                    sys.executable,
                    "-m",
                    "wizardeye",
                    "count",
                    "-i",
                    str(SIMULATED_URSUS_BAM),
                    "-r",
                    str(hg19_stem),
                    "-k",
                    str(STANDARD_KMER_LENGTH),
                    "-w",
                    str(STANDARD_OFFSET_STEP),
                    "-bn",
                    str(STANDARD_BWA_MISSING_PROB_ERR_RATE),
                    "-bo",
                    str(STANDARD_BWA_MAX_GAP_OPENINGS),
                    "-bl",
                    str(STANDARD_BWA_SEED_LENGTH),
                    "-m",
                    "mean",
                    "--exclude-tracks",
                    exclude_tracks_str,
                    "--report-output",
                    str(report_path),
                    "-d",
                    str(wizardeye_db),
                ]

                subprocess.run(
                    wizardeye_cmd,
                    capture_output=True,
                    text=True,
                    encoding="utf-8",
                    errors="replace",
                    check=True,
                    env={**subprocess.os.environ, "PYTHONPATH": str(SRC_DIR)},
                )

                assert report_path.exists(), (
                    f"wizardeye count did not create report: {report_path}"
                )

                # Copy report to persistent base directory before tmpdir is cleaned up
                dest_dir = base_path / f"threads{n_threads}"
                dest_dir.mkdir(parents=True, exist_ok=True)
                dest_report = dest_dir / "count_report.tsv"
                shutil.copy(report_path, dest_report)

                report_files[n_threads] = dest_report

        # Compare all outputs with the first one (1 thread)
        first_threads = thread_counts[0]
        first_report = report_files[first_threads]

        for n_threads in thread_counts[1:]:
            current_report = report_files[n_threads]

            # Compare report files line by line
            with open(first_report, "r") as f1, open(current_report, "r") as f2:
                first_lines = f1.readlines()
                current_lines = f2.readlines()

            if len(first_lines) != len(current_lines):
                raise AssertionError(
                    f"Count reports differ in line count between threads {first_threads} and {n_threads}: "
                    f"{len(first_lines)} vs {len(current_lines)}"
                )

            for i, (line1, line2) in enumerate(zip(first_lines, current_lines), 1):
                if line1 != line2:
                    raise AssertionError(
                        f"Count reports differ at line {i} between threads {first_threads} and {n_threads}:\n"
                        f"  Threads {first_threads}: {line1.rstrip()}\n"
                        f"  Threads {n_threads}: {line2.rstrip()}"
                    )


def test_consistency_filter_same_launch(standard_database):
    """Test if the results are the same across 10 filter rounds with the same parameters."""

    # Number of repetitions
    n_repetitions = 10

    # Use a single persistent temp directory for all repetitions
    with tempfile.TemporaryDirectory(
        prefix="wizardeye_consistency_filter_same_launch_"
    ) as base_tmpdir:
        base_path = Path(base_tmpdir)

        # Dictionary to store report paths per repetition
        report_files = {}

        # Use the shared module-scoped database
        wizardeye_db = standard_database

        # Get track names for exclude-tracks parameter
        query_stems = [
            SUS_SCROFA_FA.stem,
            CANIS_LUPUS_FA.stem,
            RATTUS_NORVEGICUS_FA.stem,
        ]
        exclude_tracks_str = ",".join(query_stems)
        hg19_stem = HG19_FA.stem

        for rep in range(n_repetitions):
            with tempfile.TemporaryDirectory(prefix=f"filter_rep{rep}_") as tmpdir:
                output_dir = Path(tmpdir) / "output"
                output_dir.mkdir(parents=True)

                # Run WizardEye filter with this repetition
                report_path = output_dir / "filter_report.tsv"
                wizardeye_cmd = [
                    sys.executable,
                    "-m",
                    "wizardeye",
                    "filter",
                    "-i",
                    str(SIMULATED_URSUS_BAM),
                    "-r",
                    str(hg19_stem),
                    "-k",
                    str(STANDARD_KMER_LENGTH),
                    "-w",
                    str(STANDARD_OFFSET_STEP),
                    "-bn",
                    str(STANDARD_BWA_MISSING_PROB_ERR_RATE),
                    "-bo",
                    str(STANDARD_BWA_MAX_GAP_OPENINGS),
                    "-bl",
                    str(STANDARD_BWA_SEED_LENGTH),
                    "-p",
                    str(STANDARD_CROSS_STRINGENCY),
                    "--exclude-tracks",
                    exclude_tracks_str,
                    "--report-output",
                    str(report_path),
                    "-d",
                    str(wizardeye_db),
                ]

                subprocess.run(
                    wizardeye_cmd,
                    capture_output=True,
                    text=True,
                    encoding="utf-8",
                    errors="replace",
                    check=True,
                    env={**subprocess.os.environ, "PYTHONPATH": str(SRC_DIR)},
                )

                assert report_path.exists(), (
                    f"wizardeye filter did not create report: {report_path}"
                )

                # Copy report to persistent base directory before tmpdir is cleaned up
                dest_dir = base_path / f"rep{rep}"
                dest_dir.mkdir(parents=True, exist_ok=True)
                dest_report = dest_dir / "filter_report.tsv"
                shutil.copy(report_path, dest_report)

                report_files[rep] = dest_report

        # Compare all outputs with the first one (rep 0)
        first_report = report_files[0]

        for rep in range(1, n_repetitions):
            current_report = report_files[rep]

            # Compare report files line by line
            with open(first_report, "r") as f1, open(current_report, "r") as f2:
                first_lines = f1.readlines()
                current_lines = f2.readlines()

            if len(first_lines) != len(current_lines):
                raise AssertionError(
                    f"Filter reports differ in line count between rep 0 and rep {rep}: "
                    f"{len(first_lines)} vs {len(current_lines)}"
                )

            for i, (line1, line2) in enumerate(zip(first_lines, current_lines), 1):
                if line1 != line2:
                    raise AssertionError(
                        f"Filter reports differ at line {i} between rep 0 and rep {rep}:\n"
                        f"  Rep 0: {line1.rstrip()}\n"
                        f"  Rep {rep}: {line2.rstrip()}"
                    )


def test_consistency_filter_parallelisation(standard_database):
    """Test if the results are the same with different numbers of cores (1, 2, 4, 8, 16).

    It uses a smaller chunk size to allow the use of different numbers of threads.
    """

    # Different thread counts to test for filter parallelisation
    desired_thread_counts = [1, 2, 4, 8, 16]
    thread_counts, should_skip, warning_msg = adjust_thread_counts_for_test(
        desired_thread_counts, min_threads_for_test=2
    )

    if should_skip:
        pytest.skip(warning_msg)

    # Use a single persistent temp directory for all thread counts
    with tempfile.TemporaryDirectory(
        prefix="wizardeye_consistency_filter_parallel_"
    ) as base_tmpdir:
        base_path = Path(base_tmpdir)

        # Dictionary to store report paths per thread count
        report_files = {}

        # Use the shared module-scoped database
        wizardeye_db = standard_database

        # Get track names for exclude-tracks parameter
        query_stems = [
            SUS_SCROFA_FA.stem,
            CANIS_LUPUS_FA.stem,
            RATTUS_NORVEGICUS_FA.stem,
        ]
        exclude_tracks_str = ",".join(query_stems)
        hg19_stem = HG19_FA.stem

        for n_threads in thread_counts:
            with tempfile.TemporaryDirectory(
                prefix=f"filter_threads{n_threads}_"
            ) as tmpdir:
                output_dir = Path(tmpdir) / "output"
                output_dir.mkdir(parents=True)

                # Run WizardEye filter with this thread count
                report_path = output_dir / "filter_report.tsv"
                wizardeye_cmd = [
                    sys.executable,
                    "-m",
                    "wizardeye",
                    "filter",
                    "-i",
                    str(SIMULATED_URSUS_BAM),
                    "-r",
                    str(hg19_stem),
                    "-k",
                    str(STANDARD_KMER_LENGTH),
                    "-w",
                    str(STANDARD_OFFSET_STEP),
                    "-bn",
                    str(STANDARD_BWA_MISSING_PROB_ERR_RATE),
                    "-bo",
                    str(STANDARD_BWA_MAX_GAP_OPENINGS),
                    "-bl",
                    str(STANDARD_BWA_SEED_LENGTH),
                    "-p",
                    str(STANDARD_CROSS_STRINGENCY),
                    "--exclude-tracks",
                    exclude_tracks_str,
                    "--report-output",
                    str(report_path),
                    "-d",
                    str(wizardeye_db),
                ]

                subprocess.run(
                    wizardeye_cmd,
                    capture_output=True,
                    text=True,
                    encoding="utf-8",
                    errors="replace",
                    check=True,
                    env={**subprocess.os.environ, "PYTHONPATH": str(SRC_DIR)},
                )

                assert report_path.exists(), (
                    f"wizardeye filter did not create report: {report_path}"
                )

                # Copy report to persistent base directory before tmpdir is cleaned up
                dest_dir = base_path / f"threads{n_threads}"
                dest_dir.mkdir(parents=True, exist_ok=True)
                dest_report = dest_dir / "filter_report.tsv"
                shutil.copy(report_path, dest_report)

                report_files[n_threads] = dest_report

        # Compare all outputs with the first one (1 thread)
        first_threads = thread_counts[0]
        first_report = report_files[first_threads]

        for n_threads in thread_counts[1:]:
            current_report = report_files[n_threads]

            # Compare report files line by line
            with open(first_report, "r") as f1, open(current_report, "r") as f2:
                first_lines = f1.readlines()
                current_lines = f2.readlines()

            if len(first_lines) != len(current_lines):
                raise AssertionError(
                    f"Filter reports differ in line count between threads {first_threads} and {n_threads}: "
                    f"{len(first_lines)} vs {len(current_lines)}"
                )

            for i, (line1, line2) in enumerate(zip(first_lines, current_lines), 1):
                if line1 != line2:
                    raise AssertionError(
                        f"Filter reports differ at line {i} between threads {first_threads} and {n_threads}:\n"
                        f"  Threads {first_threads}: {line1.rstrip()}\n"
                        f"  Threads {n_threads}: {line2.rstrip()}"
                    )


def test_export_and_filter_same_results(standard_database):
    """Test if filter and export+bed_intersect produce the same final results.

    It also indirectly tests the consistency between the new pyBigWig implementation
    and the original filtering based on bedtools for different stringency values.
    """

    # Different stringency values to test
    stringency_values = [0.01, 0.1, 0.25, 0.5, 0.75, 0.9, 0.99, 1.0]

    # Get track names for exclude-tracks parameter
    query_stems = [SUS_SCROFA_FA.stem, CANIS_LUPUS_FA.stem, RATTUS_NORVEGICUS_FA.stem]
    exclude_tracks_str = ",".join(query_stems)
    hg19_stem = HG19_FA.stem

    # Use the shared module-scoped database
    wizardeye_db = standard_database

    for cross_stringency in stringency_values:
        with tempfile.TemporaryDirectory(
            prefix=f"wizardeye_export_filter_compare_s{cross_stringency}_"
        ) as base_tmpdir:
            # Run WizardEye filter to get TSV report with excluded reads
            filter_tmpdir = tempfile.TemporaryDirectory(
                prefix="filter_", dir=base_tmpdir
            )
            filter_output_dir = Path(filter_tmpdir.name) / "output"
            filter_output_dir.mkdir(parents=True)
            filter_report_path = filter_output_dir / "filter_report.tsv"
            filter_excluded_bam = filter_output_dir / "excluded.bam"

            filter_cmd = [
                sys.executable,
                "-m",
                "wizardeye",
                "filter",
                "-i",
                str(SIMULATED_URSUS_BAM),
                "-r",
                str(hg19_stem),
                "-k",
                str(STANDARD_KMER_LENGTH),
                "-w",
                str(STANDARD_OFFSET_STEP),
                "-bn",
                str(STANDARD_BWA_MISSING_PROB_ERR_RATE),
                "-bo",
                str(STANDARD_BWA_MAX_GAP_OPENINGS),
                "-bl",
                str(STANDARD_BWA_SEED_LENGTH),
                "-p",
                str(cross_stringency),
                "--exclude-tracks",
                exclude_tracks_str,
                "--report-output",
                str(filter_report_path),
                "--excluded-output",
                str(filter_excluded_bam),
                "-d",
                str(wizardeye_db),
            ]

            subprocess.run(
                filter_cmd,
                capture_output=True,
                text=True,
                encoding="utf-8",
                errors="replace",
                check=True,
                env={**subprocess.os.environ, "PYTHONPATH": str(SRC_DIR)},
            )

            assert filter_report_path.exists(), (
                f"wizardeye filter did not create report: {filter_report_path}"
            )
            assert filter_excluded_bam.exists(), (
                f"wizardeye filter did not create excluded BAM: {filter_excluded_bam}"
            )

            # Run WizardEye export to get BED mask
            export_tmpdir = tempfile.TemporaryDirectory(
                prefix="export_", dir=base_tmpdir
            )
            export_output_dir = Path(export_tmpdir.name) / "output"
            export_output_dir.mkdir(parents=True)
            export_mask_path = export_output_dir / "export_mask.bed"

            export_cmd = [
                sys.executable,
                "-m",
                "wizardeye",
                "export",
                "-r",
                str(hg19_stem),
                "-k",
                str(STANDARD_KMER_LENGTH),
                "-w",
                str(STANDARD_OFFSET_STEP),
                "-bn",
                str(STANDARD_BWA_MISSING_PROB_ERR_RATE),
                "-bo",
                str(STANDARD_BWA_MAX_GAP_OPENINGS),
                "-bl",
                str(STANDARD_BWA_SEED_LENGTH),
                "-p",
                str(cross_stringency),
                "--exclude-tracks",
                exclude_tracks_str,
                "-o",
                str(export_mask_path),
                "-d",
                str(wizardeye_db),
            ]

            subprocess.run(
                export_cmd,
                capture_output=True,
                text=True,
                encoding="utf-8",
                errors="replace",
                check=True,
                env={**subprocess.os.environ, "PYTHONPATH": str(SRC_DIR)},
            )

            assert export_mask_path.exists(), (
                f"wizardeye export did not create mask: {export_mask_path}"
            )

            # Use bedtools intersect to find reads overlapping the export mask
            intersect_tmpdir = tempfile.TemporaryDirectory(
                prefix="intersect_", dir=base_tmpdir
            )
            intersect_output_dir = Path(intersect_tmpdir.name) / "output"
            intersect_output_dir.mkdir(parents=True)
            intersect_bam_path = intersect_output_dir / "intersect.bam"

            # Sort the input BAM first for bedtools intersect
            sorted_input_bam = intersect_output_dir / "input_sorted.bam"
            subprocess.run(
                [
                    "samtools",
                    "sort",
                    str(SIMULATED_URSUS_BAM),
                    "-o",
                    str(sorted_input_bam),
                ],
                check=True,
                capture_output=True,
                text=True,
            )

            # Run bedtools intersect to find reads overlapping the mask
            # -abam keeps BAM output format, -wa writes all records from A (the BAM)
            with open(str(intersect_bam_path), "w") as out_file:
                subprocess.run(
                    [
                        "bedtools",
                        "intersect",
                        "-abam",
                        str(sorted_input_bam),
                        "-b",
                        str(export_mask_path),
                        "-wa",
                    ],
                    stdout=out_file,
                    check=True,
                    text=True,
                )

            filter_read_keys = extract_mapped_read_keys(filter_excluded_bam)
            intersect_read_keys = extract_mapped_read_keys(intersect_bam_path)

            if filter_read_keys != intersect_read_keys:
                only_in_filter = filter_read_keys - intersect_read_keys
                only_in_intersect = intersect_read_keys - filter_read_keys
                raise AssertionError(
                    f"Read sets differ for stringency {cross_stringency}:\n"
                    f"  Filtered out using filter: {len(only_in_filter)} reads (e.g., {list(only_in_filter)[:5]})\n"
                    f"  Filtered out using the exported mask: {len(only_in_intersect)} reads (e.g., {list(only_in_intersect)[:5]})"
                )

            # Clean up temp directories
            filter_tmpdir.cleanup()
            export_tmpdir.cleanup()
            intersect_tmpdir.cleanup()


def test_export_and_filter_same_results_with_mf(standard_database):
    """Test filter and export+bed_intersect give the same results with -mf.

    It also indirectly tests the consistency between the new pyBigWig implementation
    and the original filtering based on bedtools for different stringency values.
    """

    # Different stringency values to test
    stringency_values = [0.01, 0.1, 0.25, 0.5, 0.75, 0.9, 0.99, 1.0]
    min_frequencies = [1, 2, 3]

    # Get track names for exclude-tracks parameter
    query_stems = [SUS_SCROFA_FA.stem, CANIS_LUPUS_FA.stem, RATTUS_NORVEGICUS_FA.stem]
    exclude_tracks_str = ",".join(query_stems)
    hg19_stem = HG19_FA.stem

    # Use the shared module-scoped database
    wizardeye_db = standard_database

    for mf in min_frequencies:
        for cross_stringency in stringency_values:
            with tempfile.TemporaryDirectory(
                prefix=f"wizardeye_export_filter_compare_s{cross_stringency}_"
            ) as base_tmpdir:
                # Run WizardEye filter to get TSV report with excluded reads
                filter_tmpdir = tempfile.TemporaryDirectory(
                    prefix="filter_", dir=base_tmpdir
                )
                filter_output_dir = Path(filter_tmpdir.name) / "output"
                filter_output_dir.mkdir(parents=True)
                filter_report_path = filter_output_dir / "filter_report.tsv"
                filter_excluded_bam = filter_output_dir / "excluded.bam"

                filter_cmd = [
                    sys.executable,
                    "-m",
                    "wizardeye",
                    "filter",
                    "-i",
                    str(SIMULATED_URSUS_BAM),
                    "-r",
                    str(hg19_stem),
                    "-k",
                    str(STANDARD_KMER_LENGTH),
                    "-w",
                    str(STANDARD_OFFSET_STEP),
                    "-bn",
                    str(STANDARD_BWA_MISSING_PROB_ERR_RATE),
                    "-bo",
                    str(STANDARD_BWA_MAX_GAP_OPENINGS),
                    "-bl",
                    str(STANDARD_BWA_SEED_LENGTH),
                    "-p",
                    str(cross_stringency),
                    "-mf",
                    str(mf),
                    "--exclude-tracks",
                    exclude_tracks_str,
                    "--report-output",
                    str(filter_report_path),
                    "--excluded-output",
                    str(filter_excluded_bam),
                    "-d",
                    str(wizardeye_db),
                ]

                subprocess.run(
                    filter_cmd,
                    capture_output=False,
                    text=True,
                    encoding="utf-8",
                    errors="replace",
                    check=True,
                    env={**subprocess.os.environ, "PYTHONPATH": str(SRC_DIR)},
                )

                assert filter_report_path.exists(), (
                    f"wizardeye filter did not create report: {filter_report_path}"
                )
                assert filter_excluded_bam.exists(), (
                    f"wizardeye filter did not create excluded BAM: {filter_excluded_bam}"
                )

                # Run WizardEye export to get BED mask
                export_tmpdir = tempfile.TemporaryDirectory(
                    prefix="export_", dir=base_tmpdir
                )
                export_output_dir = Path(export_tmpdir.name) / "output"
                export_output_dir.mkdir(parents=True)
                export_mask_path = export_output_dir / "export_mask.bed"

                export_cmd = [
                    sys.executable,
                    "-m",
                    "wizardeye",
                    "export",
                    "-r",
                    str(hg19_stem),
                    "-k",
                    str(STANDARD_KMER_LENGTH),
                    "-w",
                    str(STANDARD_OFFSET_STEP),
                    "-bn",
                    str(STANDARD_BWA_MISSING_PROB_ERR_RATE),
                    "-bo",
                    str(STANDARD_BWA_MAX_GAP_OPENINGS),
                    "-bl",
                    str(STANDARD_BWA_SEED_LENGTH),
                    "-p",
                    str(cross_stringency),
                    "-mf",
                    str(mf),
                    "--exclude-tracks",
                    exclude_tracks_str,
                    "-o",
                    str(export_mask_path),
                    "-d",
                    str(wizardeye_db),
                ]

                subprocess.run(
                    export_cmd,
                    capture_output=False,
                    text=True,
                    encoding="utf-8",
                    errors="replace",
                    check=True,
                    env={**subprocess.os.environ, "PYTHONPATH": str(SRC_DIR)},
                )

                assert export_mask_path.exists(), (
                    f"wizardeye export did not create mask: {export_mask_path}"
                )

                # Use bedtools intersect to find reads overlapping the export mask
                intersect_tmpdir = tempfile.TemporaryDirectory(
                    prefix="intersect_", dir=base_tmpdir
                )
                intersect_output_dir = Path(intersect_tmpdir.name) / "output"
                intersect_output_dir.mkdir(parents=True)
                intersect_bam_path = intersect_output_dir / "intersect.bam"

                # Sort the input BAM first for bedtools intersect
                sorted_input_bam = intersect_output_dir / "input_sorted.bam"
                subprocess.run(
                    [
                        "samtools",
                        "sort",
                        str(SIMULATED_URSUS_BAM),
                        "-o",
                        str(sorted_input_bam),
                    ],
                    check=True,
                    capture_output=True,
                    text=True,
                )

                # Run bedtools intersect to find reads overlapping the mask
                # -abam keeps BAM output format, -wa writes all records from A (the BAM)
                with open(str(intersect_bam_path), "w") as out_file:
                    subprocess.run(
                        [
                            "bedtools",
                            "intersect",
                            "-abam",
                            str(sorted_input_bam),
                            "-b",
                            str(export_mask_path),
                            "-wa",
                        ],
                        stdout=out_file,
                        check=True,
                        text=True,
                    )

                filter_read_keys = extract_mapped_read_keys(filter_excluded_bam)
                intersect_read_keys = extract_mapped_read_keys(intersect_bam_path)

                if filter_read_keys != intersect_read_keys:
                    only_in_filter = filter_read_keys - intersect_read_keys
                    only_in_intersect = intersect_read_keys - filter_read_keys
                    raise AssertionError(
                        f"Read sets differ for stringency {cross_stringency} and minimum frequency {mf}:\n"
                        f"  Filtered out using filter: {len(only_in_filter)} reads (e.g., {list(only_in_filter)[:5]})\n"
                        f"  Filtered out using the exported mask: {len(only_in_intersect)} reads (e.g., {list(only_in_intersect)[:5]})"
                    )

                # Clean up temp directories
                filter_tmpdir.cleanup()
                export_tmpdir.cleanup()
                intersect_tmpdir.cleanup()


def test_consistency_filter_report_bam_split_and_min_frequency(standard_database):
    """Test filter report matches the excluded/kept BAM split.

    Also verifies that -mf only toggles exclusion while associated_tracks stays
    fixed.

    Verifies, for several stringency and -mf values:
      1. A read flagged ``true`` (filtered_out) is in the excluded BAM and absent from the
         kept BAM; a read flagged ``false`` is in the kept BAM and absent from the
         excluded BAM.
      2. ``associated_tracks`` is identical for every -mf value at a given stringency - only
         the ``filtered_out`` flag varies with -mf.
      3. When the number of associated tracks is strictly below -mf, the report marks the
         read as ``false`` (not excluded).
    """
    stringency_values = [0.01, 0.1, 0.25, 0.5, 0.75, 0.9, 0.99, 1.0]

    min_frequencies: list[int | None] = [None, 1, 2, 3]

    query_stems = [SUS_SCROFA_FA.stem, CANIS_LUPUS_FA.stem, RATTUS_NORVEGICUS_FA.stem]
    exclude_tracks_str = ",".join(query_stems)
    hg19_stem = HG19_FA.stem

    with tempfile.TemporaryDirectory(
        prefix="wizardeye_filter_report_mf_"
    ) as base_tmpdir:
        base_path = Path(base_tmpdir)

        # Use the shared module-scoped database
        wizardeye_db = standard_database

        runs_dir = base_path / "runs"
        runs_dir.mkdir(parents=True, exist_ok=True)

        input_mapped_keys = extract_mapped_read_keys(SIMULATED_URSUS_BAM)

        for cross_stringency in stringency_values:
            reports = {}  # mf -> {read_key: (flag, tracks)}
            bam_keys = {}  # mf -> (kept_keys, excluded_keys)

            for mf in min_frequencies:
                mf_label = "none" if mf is None else str(mf)
                report_path = runs_dir / f"report_s{cross_stringency}_mf{mf_label}.tsv"
                kept_bam = runs_dir / f"kept_s{cross_stringency}_mf{mf_label}.bam"
                excluded_bam = (
                    runs_dir / f"excluded_s{cross_stringency}_mf{mf_label}.bam"
                )

                cmd = [
                    sys.executable,
                    "-m",
                    "wizardeye",
                    "filter",
                    "-i",
                    str(SIMULATED_URSUS_BAM),
                    "-r",
                    str(hg19_stem),
                    "-k",
                    str(STANDARD_KMER_LENGTH),
                    "-w",
                    str(STANDARD_OFFSET_STEP),
                    "-bn",
                    str(STANDARD_BWA_MISSING_PROB_ERR_RATE),
                    "-bo",
                    str(STANDARD_BWA_MAX_GAP_OPENINGS),
                    "-bl",
                    str(STANDARD_BWA_SEED_LENGTH),
                    "-p",
                    str(cross_stringency),
                    "--exclude-tracks",
                    exclude_tracks_str,
                    "--report-output",
                    str(report_path),
                    "--kept-output",
                    str(kept_bam),
                    "--excluded-output",
                    str(excluded_bam),
                    "-d",
                    str(wizardeye_db),
                ]
                if mf is not None:
                    cmd += ["-mf", str(mf)]

                subprocess.run(
                    cmd,
                    capture_output=True,
                    text=True,
                    encoding="utf-8",
                    errors="replace",
                    check=True,
                    env={**subprocess.os.environ, "PYTHONPATH": str(SRC_DIR)},
                )

                assert report_path.exists(), (
                    f"wizardeye filter did not create report: {report_path}"
                )
                assert kept_bam.exists(), (
                    f"wizardeye filter did not create kept BAM: {kept_bam}"
                )
                assert excluded_bam.exists(), (
                    f"wizardeye filter did not create excluded BAM: {excluded_bam}"
                )

                reports[mf] = parse_filter_report(report_path)
                bam_keys[mf] = (
                    extract_mapped_read_keys(kept_bam),
                    extract_mapped_read_keys(excluded_bam),
                )

            # 1. Report flag must match the excluded/kept BAM split.
            for mf in min_frequencies:
                kept_keys, excluded_keys = bam_keys[mf]
                rows = reports[mf]

                assert not (kept_keys & excluded_keys), (
                    f"A mapped read is in both kept and excluded BAM "
                    f"(stringency={cross_stringency}, mf={mf})"
                )
                assert input_mapped_keys == (kept_keys | excluded_keys), (
                    f"Kept + excluded BAM do not cover all mapped input reads "
                    f"(stringency={cross_stringency}, mf={mf})"
                )

                for read_key, (flag, tracks) in rows.items():
                    if read_key not in input_mapped_keys:
                        assert flag == "false", (
                            f"Non-mapped read should be flagged false: {read_key} "
                            f"(stringency={cross_stringency}, mf={mf})"
                        )
                        continue
                    if flag == "true":
                        assert read_key in excluded_keys, (
                            f"Read flagged true not in excluded BAM: {read_key} "
                            f"(stringency={cross_stringency}, mf={mf})"
                        )
                        assert read_key not in kept_keys, (
                            f"Read flagged true also in kept BAM: {read_key} "
                            f"(stringency={cross_stringency}, mf={mf})"
                        )
                        assert tracks, (
                            f"Read flagged true has no associated tracks: {read_key} "
                            f"(stringency={cross_stringency}, mf={mf})"
                        )
                    else:
                        assert read_key in kept_keys, (
                            f"Read flagged false not in kept BAM: {read_key} "
                            f"(stringency={cross_stringency}, mf={mf})"
                        )
                        assert read_key not in excluded_keys, (
                            f"Read flagged false also in excluded BAM: {read_key} "
                            f"(stringency={cross_stringency}, mf={mf})"
                        )

            # 2. associated_tracks must be identical across every -mf value; only the
            #    filtered_out flag may change.
            common_keys = set(reports[min_frequencies[0]].keys())
            for mf in min_frequencies[1:]:
                rows = reports[mf]
                assert set(rows.keys()) == common_keys, (
                    f"Report read keys differ between mf={min_frequencies[0]} and "
                    f"mf={mf} (stringency={cross_stringency})"
                )
                for read_key in common_keys:
                    base_tracks = reports[min_frequencies[0]][read_key][1]
                    cur_tracks = rows[read_key][1]
                    assert cur_tracks == base_tracks, (
                        f"associated_tracks changed with -mf for {read_key}: "
                        f"{base_tracks} -> {cur_tracks} "
                        f"(stringency={cross_stringency}, "
                        f"mf {min_frequencies[0]} -> {mf})"
                    )

            # Monotonicity: as -mf grows, exclusion can only shrink. A read excluded at
            # a higher -mf must also be excluded at every lower -mf, and the no-frequency
            # path must match -mf 1 (both mean "at least one overlapping track").
            mf_values = [1, 2, 3]
            for idx in range(len(mf_values) - 1):
                higher = reports[mf_values[idx + 1]]
                lower = reports[mf_values[idx]]
                for read_key in common_keys:
                    if higher[read_key][0] == "true":
                        assert lower[read_key][0] == "true", (
                            f"Read excluded at mf={mf_values[idx + 1]} but not at "
                            f"mf={mf_values[idx]} (monotonicity broken): {read_key} "
                            f"(stringency={cross_stringency})"
                        )
            for read_key in common_keys:
                assert reports[None][read_key][0] == reports[1][read_key][0], (
                    f"filtered_out differs between mf=None and mf=1 for {read_key} "
                    f"(stringency={cross_stringency})"
                )

            # 3. With -mf set, a read with fewer associated tracks than -mf cannot be
            #    excluded.
            for mf in min_frequencies:
                if mf is None:
                    continue
                rows = reports[mf]
                for read_key, (flag, tracks) in rows.items():
                    if len(tracks) < mf:
                        assert flag == "false", (
                            f"Read with {len(tracks)} associated tracks (< mf={mf}) "
                            f"should be flagged false: {read_key} "
                            f"(stringency={cross_stringency})"
                        )


def test_consistency_align_parallel_same_db():
    """Test if multiple genomes can be aligned simultaneously in parallel on the same database.

    This test verifies that running multiple WizardEye align processes concurrently
    on the same database with different query FASTA files:
    1. Completes without errors (no race conditions or lock issues)
    2. Produces results identical to sequential alignments

    Uses all 3 available query FASTAs (sus_scrofa, canis_lupus, rattus_norvegicus).
    """
    import concurrent.futures

    # List of query FASTAs to test simultaneously
    query_fastas = [SUS_SCROFA_FA, CANIS_LUPUS_FA, RATTUS_NORVEGICUS_FA]

    with tempfile.TemporaryDirectory(
        prefix="wizardeye_consistency_parallel_db_"
    ) as base_tmpdir:
        base_path = Path(base_tmpdir)

        # Create a shared database
        wizardeye_db = base_path / "database"
        wizardeye_db.mkdir(parents=True)

        # Initialize WizardEye database
        subprocess.run(
            [
                sys.executable,
                "-m",
                "wizardeye",
                "database",
                "init",
                "-d",
                str(base_path),
            ],
            check=True,
            capture_output=True,
            text=True,
            encoding="utf-8",
            errors="replace",
            env={**subprocess.os.environ, "PYTHONPATH": str(SRC_DIR)},
        )

        def run_align(query_fasta):
            """Run a single align command and return the output paths."""
            subprocess.run(
                [
                    sys.executable,
                    "-m",
                    "wizardeye",
                    "align",
                    "-i",
                    str(query_fasta),
                    "-r",
                    str(HG19_FA),
                    "-k",
                    str(STANDARD_KMER_LENGTH),
                    "-w",
                    str(STANDARD_OFFSET_STEP),
                    "-bn",
                    str(STANDARD_BWA_MISSING_PROB_ERR_RATE),
                    "-bo",
                    str(STANDARD_BWA_MAX_GAP_OPENINGS),
                    "-bl",
                    str(STANDARD_BWA_SEED_LENGTH),
                    "-j",
                    str(STANDARD_N_THREADS),
                    "-cs",
                    str(STANDARD_CHUNK_SIZE),
                    "-d",
                    str(wizardeye_db),
                ],
                capture_output=True,
                text=True,
                encoding="utf-8",
                errors="replace",
                check=True,
                env={**subprocess.os.environ, "PYTHONPATH": str(SRC_DIR)},
            )

            # Return the track directory path for this query
            hg19_stem = HG19_FA.stem
            query_stem = query_fasta.stem
            track_pattern = f"{query_stem}_k{STANDARD_KMER_LENGTH}_w{STANDARD_OFFSET_STEP}_bwa{STANDARD_BWA_HASH}"
            track_dir = wizardeye_db / hg19_stem / track_pattern

            if not track_dir.exists():
                raise RuntimeError(f"Track directory not found: {track_dir}")

            map_all_bw = track_dir / "map_all.bw"
            map_uniq_bw = track_dir / "map_uniq.bw"

            if not map_all_bw.exists():
                raise RuntimeError(f"map_all.bw not found: {map_all_bw}")
            if not map_uniq_bw.exists():
                raise RuntimeError(f"map_uniq.bw not found: {map_uniq_bw}")

            return (query_fasta, map_all_bw, map_uniq_bw)

        # Run all alignments in parallel
        with concurrent.futures.ThreadPoolExecutor(
            max_workers=len(query_fastas)
        ) as executor:
            future_to_query = {
                executor.submit(run_align, query_fa): query_fa
                for query_fa in query_fastas
            }

            parallel_results = {}
            for future in concurrent.futures.as_completed(future_to_query):
                query_fa = future_to_query[future]
                try:
                    result = future.result()
                    parallel_results[query_fa] = result
                except (
                    RuntimeError,
                    subprocess.CalledProcessError,
                    ValueError,
                    OSError,
                ) as e:
                    raise RuntimeError(
                        f"Parallel alignment failed for {query_fa.name}: {e}"
                    )

        # Verify all alignments completed successfully
        assert len(parallel_results) == len(query_fastas), (
            f"Expected {len(query_fastas)} results, got {len(parallel_results)}"
        )

        # Store parallel results for comparison
        parallel_bw_files = {}
        for query_fa, (_, map_all_bw, map_uniq_bw) in parallel_results.items():
            parallel_bw_files[query_fa] = (map_all_bw, map_uniq_bw)

        # Now run sequential alignments for comparison
        # Use a separate database to ensure isolation
        with tempfile.TemporaryDirectory(
            prefix="wizardeye_consistency_parallel_db_sequential_"
        ) as seq_tmpdir:
            seq_base_path = Path(seq_tmpdir)
            seq_db = seq_base_path / "database"
            seq_db.mkdir(parents=True)

            # Initialize sequential database
            subprocess.run(
                [
                    sys.executable,
                    "-m",
                    "wizardeye",
                    "database",
                    "init",
                    "-d",
                    str(seq_tmpdir),
                ],
                check=True,
                capture_output=True,
                text=True,
                encoding="utf-8",
                errors="replace",
                env={**subprocess.os.environ, "PYTHONPATH": str(SRC_DIR)},
            )

            sequential_bw_files = {}
            for query_fa in query_fastas:
                result = subprocess.run(
                    [
                        sys.executable,
                        "-m",
                        "wizardeye",
                        "align",
                        "-i",
                        str(query_fa),
                        "-r",
                        str(HG19_FA),
                        "-k",
                        str(STANDARD_KMER_LENGTH),
                        "-w",
                        str(STANDARD_OFFSET_STEP),
                        "-bn",
                        str(STANDARD_BWA_MISSING_PROB_ERR_RATE),
                        "-bo",
                        str(STANDARD_BWA_MAX_GAP_OPENINGS),
                        "-bl",
                        str(STANDARD_BWA_SEED_LENGTH),
                        "-j",
                        str(STANDARD_N_THREADS),
                        "-cs",
                        str(STANDARD_CHUNK_SIZE),
                        "-d",
                        str(seq_db),
                    ],
                    capture_output=True,
                    text=True,
                    encoding="utf-8",
                    errors="replace",
                    check=True,
                    env={**subprocess.os.environ, "PYTHONPATH": str(SRC_DIR)},
                )

                # Locate output
                hg19_stem = HG19_FA.stem
                query_stem = query_fa.stem
                track_pattern = f"{query_stem}_k{STANDARD_KMER_LENGTH}_w{STANDARD_OFFSET_STEP}_bwa{STANDARD_BWA_HASH}"
                track_dir = seq_db / hg19_stem / track_pattern

                if not track_dir.exists():
                    raise RuntimeError(
                        f"Sequential track directory not found: {track_dir}"
                    )

                map_all_bw = track_dir / "map_all.bw"
                map_uniq_bw = track_dir / "map_uniq.bw"

                if not map_all_bw.exists():
                    raise RuntimeError(f"Sequential map_all.bw not found: {map_all_bw}")
                if not map_uniq_bw.exists():
                    raise RuntimeError(
                        f"Sequential map_uniq.bw not found: {map_uniq_bw}"
                    )

                sequential_bw_files[query_fa] = (map_all_bw, map_uniq_bw)

            # Compare parallel and sequential results for each query FASTA
            for query_fa in query_fastas:
                par_map_all, par_map_uniq = parallel_bw_files[query_fa]
                seq_map_all, seq_map_uniq = sequential_bw_files[query_fa]

                with tempfile.TemporaryDirectory(
                    prefix=f"wizardeye_compare_parallel_db_{query_fa.stem}_"
                ) as compare_dir:
                    compare_path = Path(compare_dir)

                    # Convert to bedGraph for comparison
                    par_all_bg = compare_path / "parallel_map_all.bg"
                    par_uniq_bg = compare_path / "parallel_map_uniq.bg"
                    seq_all_bg = compare_path / "sequential_map_all.bg"
                    seq_uniq_bg = compare_path / "sequential_map_uniq.bg"

                    subprocess.run(
                        ["bigWigToBedGraph", str(par_map_all), str(par_all_bg)],
                        check=True,
                    )
                    subprocess.run(
                        ["bigWigToBedGraph", str(par_map_uniq), str(par_uniq_bg)],
                        check=True,
                    )
                    subprocess.run(
                        ["bigWigToBedGraph", str(seq_map_all), str(seq_all_bg)],
                        check=True,
                    )
                    subprocess.run(
                        ["bigWigToBedGraph", str(seq_map_uniq), str(seq_uniq_bg)],
                        check=True,
                    )

                    # Compare map_all
                    assert compare_bedgraph_files(
                        par_all_bg,
                        seq_all_bg,
                        f"Parallel vs Sequential map_all ({query_fa.name})",
                    ), (
                        f"map_all.bw differs between parallel and sequential for {query_fa.name}"
                    )

                    # Compare map_uniq
                    assert compare_bedgraph_files(
                        par_uniq_bg,
                        seq_uniq_bg,
                        f"Parallel vs Sequential map_uniq ({query_fa.name})",
                    ), (
                        f"map_uniq.bw differs between parallel and sequential for {query_fa.name}"
                    )


# -- Protection against misuse --


def test_filter_paired_end_vs_single_end(standard_database):
    """Test that WizardEye correctly handles different read types.

    Verifies two behaviors:
    1. Rejects BAM files with paired-end reads (flag 0x1) with an appropriate error.
    2. Accepts and successfully processes BAM files with single-end reads.
    """
    sequences = get_all_fasta_sequence_info(HG19_FA)
    sq_headers = "".join(
        [f"@SQ\tSN:{name}\tLN:{length}\n" for name, length in sequences]
    )

    sam_header = (
        f"@HD\tVN:1.6\tSO:coordinate\n"
        f"{sq_headers}"
        f"@PG\tID:bwa\tPN:bwa\tCL:bwa aln -n 0.01 -o 2 -l 16500 ref.fa reads.fq\n"
    )

    ref_seq_name = sequences[0][0]

    sam_content_paired = sam_header + (
        f"read1\t99\t{ref_seq_name}\t100\t60\t50M\t=\t200\t150\tACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC\tIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\tAS:i:50\n"
        f"read1\t147\t{ref_seq_name}\t200\t60\t50M\t=\t100\t-150\tACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC\tIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\tAS:i:50\n"
    )

    sam_content_single = sam_header + (
        f"read3_se\t0\t{ref_seq_name}\t300\t60\t50M\t*\t0\t0\tTGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATG\tIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\tAS:i:50\n"
        f"read4_se\t0\t{ref_seq_name}\t400\t60\t50M\t*\t0\t0\tTGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATG\tIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\tAS:i:50\n"
    )

    with tempfile.TemporaryDirectory(prefix="wizardeye_reads_test_") as tmpdir:
        tmp_path = Path(tmpdir)

        # Use the shared module-scoped database
        wizardeye_db = standard_database

        query_stems = [
            SUS_SCROFA_FA.stem,
            CANIS_LUPUS_FA.stem,
            RATTUS_NORVEGICUS_FA.stem,
        ]
        exclude_tracks_str = ",".join(query_stems)
        hg19_stem = HG19_FA.stem

        base_cmd = [
            sys.executable,
            "-m",
            "wizardeye",
            "filter",
            "-r",
            str(hg19_stem),
            "-k",
            str(STANDARD_KMER_LENGTH),
            "-w",
            str(STANDARD_OFFSET_STEP),
            "-bn",
            str(STANDARD_BWA_MISSING_PROB_ERR_RATE),
            "-bo",
            str(STANDARD_BWA_MAX_GAP_OPENINGS),
            "-bl",
            str(STANDARD_BWA_SEED_LENGTH),
            "-p",
            str(STANDARD_CROSS_STRINGENCY),
            "--exclude-tracks",
            exclude_tracks_str,
            "-d",
            str(wizardeye_db),
        ]

        # Paired-end
        sam_paired = tmp_path / "paired_end.sam"
        bam_paired = tmp_path / "paired_end.bam"
        sam_paired.write_text(sam_content_paired)
        subprocess.run(
            ["samtools", "view", "-b", str(sam_paired), "-o", str(bam_paired)],
            check=True,
        )

        cmd_paired = base_cmd[:4] + ["-i", str(bam_paired)] + base_cmd[4:]

        result_paired = subprocess.run(
            cmd_paired,
            capture_output=True,
            text=True,
            encoding="utf-8",
            errors="replace",
            check=False,
            env={**subprocess.os.environ, "PYTHONPATH": str(SRC_DIR)},
        )

        assert result_paired.returncode != 0, (
            "WizardEye should have failed with paired-end reads."
        )
        error_output = result_paired.stderr + result_paired.stdout
        assert "paired" in error_output.lower(), (
            f"Expected error about paired-end reads, got: {error_output[:500]}"
        )

        # Single-end
        sam_single = tmp_path / "single_end.sam"
        bam_single = tmp_path / "single_end.bam"
        sam_single.write_text(sam_content_single)
        subprocess.run(
            ["samtools", "view", "-b", str(sam_single), "-o", str(bam_single)],
            check=True,
        )

        cmd_single = base_cmd[:4] + ["-i", str(bam_single)] + base_cmd[4:]

        subprocess.run(
            cmd_single,
            capture_output=True,
            text=True,
            encoding="utf-8",
            errors="replace",
            check=True,
            env={**subprocess.os.environ, "PYTHONPATH": str(SRC_DIR)},
        )


def test_count_paired_end_vs_single_end(standard_database):
    """Test that WizardEye count command correctly handles different read types.

    Verifies two behaviors:
    1. Rejects BAM files with paired-end reads (flag 0x1) with an appropriate error.
    2. Accepts and successfully processes BAM files with single-end reads.
    """
    sequences = get_all_fasta_sequence_info(HG19_FA)
    sq_headers = "".join(
        [f"@SQ\tSN:{name}\tLN:{length}\n" for name, length in sequences]
    )

    sam_header = (
        f"@HD\tVN:1.6\tSO:coordinate\n"
        f"{sq_headers}"
        f"@PG\tID:bwa\tPN:bwa\tCL:bwa aln -n 0.01 -o 2 -l 16500 ref.fa reads.fq\n"
    )

    ref_seq_name = sequences[0][0]

    sam_content_paired = sam_header + (
        f"read1\t99\t{ref_seq_name}\t100\t60\t50M\t=\t200\t150\tACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC\tIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\tAS:i:50\n"
        f"read1\t147\t{ref_seq_name}\t200\t60\t50M\t=\t100\t-150\tACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC\tIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\tAS:i:50\n"
    )

    sam_content_single = sam_header + (
        f"read3_se\t0\t{ref_seq_name}\t300\t60\t50M\t*\t0\t0\tTGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATG\tIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\tAS:i:50\n"
        f"read4_se\t0\t{ref_seq_name}\t400\t60\t50M\t*\t0\t0\tTGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATG\tIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\tAS:i:50\n"
    )

    with tempfile.TemporaryDirectory(prefix="wizardeye_reads_test_") as tmpdir:
        tmp_path = Path(tmpdir)

        # Use the shared module-scoped database
        wizardeye_db = standard_database

        query_stems = [
            SUS_SCROFA_FA.stem,
            CANIS_LUPUS_FA.stem,
            RATTUS_NORVEGICUS_FA.stem,
        ]
        exclude_tracks_str = ",".join(query_stems)
        hg19_stem = HG19_FA.stem

        base_cmd = [
            sys.executable,
            "-m",
            "wizardeye",
            "count",
            "-r",
            str(hg19_stem),
            "-k",
            str(STANDARD_KMER_LENGTH),
            "-w",
            str(STANDARD_OFFSET_STEP),
            "-bn",
            str(STANDARD_BWA_MISSING_PROB_ERR_RATE),
            "-bo",
            str(STANDARD_BWA_MAX_GAP_OPENINGS),
            "-bl",
            str(STANDARD_BWA_SEED_LENGTH),
            "-m",
            "mean",
            "--exclude-tracks",
            exclude_tracks_str,
            "-d",
            str(wizardeye_db),
        ]

        # Paired-end
        sam_paired = tmp_path / "paired_end.sam"
        bam_paired = tmp_path / "paired_end.bam"
        sam_paired.write_text(sam_content_paired)
        subprocess.run(
            ["samtools", "view", "-b", str(sam_paired), "-o", str(bam_paired)],
            check=True,
        )

        cmd_paired = base_cmd[:4] + ["-i", str(bam_paired)] + base_cmd[4:]

        result_paired = subprocess.run(
            cmd_paired,
            capture_output=True,
            text=True,
            encoding="utf-8",
            errors="replace",
            check=False,
            env={**subprocess.os.environ, "PYTHONPATH": str(SRC_DIR)},
        )

        assert result_paired.returncode != 0, (
            "WizardEye count should have failed with paired-end reads."
        )
        error_output = result_paired.stderr + result_paired.stdout
        assert "paired" in error_output.lower(), (
            f"Expected error about paired-end reads, got: {error_output[:500]}"
        )

        # Single-end
        sam_single = tmp_path / "single_end.sam"
        bam_single = tmp_path / "single_end.bam"
        sam_single.write_text(sam_content_single)
        subprocess.run(
            ["samtools", "view", "-b", str(sam_single), "-o", str(bam_single)],
            check=True,
        )

        cmd_single = base_cmd[:4] + ["-i", str(bam_single)] + base_cmd[4:]

        subprocess.run(
            cmd_single,
            capture_output=True,
            text=True,
            encoding="utf-8",
            errors="replace",
            check=True,
            env={**subprocess.os.environ, "PYTHONPATH": str(SRC_DIR)},
        )
