"""Common constants and configuration for WizardEye e2e tests.

This module provides shared constants used across all end-to-end test files.
"""

from pathlib import Path

# Path constants
PROJECT_ROOT = Path(__file__).parent.parent.parent
FIXTURES_DIR = Path(__file__).parent.parent / "fixtures"
SRC_DIR = PROJECT_ROOT / "src"
SCRIPT_PATH = Path(__file__).parent.parent / "generate_cross_mappability_filter_bwa.sh"

# Fixture files
HG19_FA = FIXTURES_DIR / "hg19_chr1_25_1kbp.fa"
SUS_SCROFA_FA = FIXTURES_DIR / "sus_scrofa_chr1_25_1kbp.fa"
CANIS_LUPUS_FA = FIXTURES_DIR / "canis_lupus_chr1_25_1kbp.fa"
RATTUS_NORVEGICUS_FA = FIXTURES_DIR / "rattus_norvegicus_chr1_25_1kbp.fa"
SIMULATED_URSUS_BAM = FIXTURES_DIR / "ursus_1000000.uniq.L35MQ25.bam"

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

# Contaminant FASTA files to test
CONTAMINANT_FAS = [
    SUS_SCROFA_FA,
    CANIS_LUPUS_FA,
    RATTUS_NORVEGICUS_FA,
]
