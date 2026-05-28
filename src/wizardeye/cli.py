# -*- coding: utf-8 -*-

"""Command line interface for WizardEye."""

import os
import subprocess
import time
import typer
import sys
from importlib import metadata
from pathlib import Path
from typing import List, Optional

from .mappability import create_mappability_track
from .filter import generate_global_mask, count_k_mers_on_bam, filter_bam
from .utils import from_charlist_to_list, log, BWAParameters, get_bwa_params_hash
from .version import DISPLAY_VERSION, PACKAGE_VERSION, print_version_message
from .db import (
    import_track,
    init_db,
    clean_db,
    print_full_catalogue,
    update_track_tags,
    valid_database,
    check_track_exists,
    from_tags_get_tracks,
    get_refs,
    get_tracks,
    resolve_requested_track_names,
    migrate_database,
)

app = typer.Typer(
    help="WizardEye: A Python tool to create, manage, and filter by cross-mappability tracks.",
    no_args_is_help=True,
    add_completion=False,
    suggest_commands=False,
)

# Database command group
db_app = typer.Typer(help="Initialize, validate, or inspect tracks in a WizardEye database.")
update_app = typer.Typer(help="Update operations for database tracks.")
db_app.add_typer(update_app, name="update")
app.add_typer(db_app, name="database", no_args_is_help=True)

# -- Helper functions --


def _get_runtime_version() -> str:
    """Get installed package version, or fallback display version in source mode.

    Returns:
            str: The installed package version or the display version.
    """
    try:
        installed_version = metadata.version("wizardeye")
        if installed_version == PACKAGE_VERSION:
            return DISPLAY_VERSION
        return installed_version
    except metadata.PackageNotFoundError:
        return DISPLAY_VERSION


def _version_callback(value: bool) -> None:
    """Callback function to display version information if requested and exit.

    Args:
            value (bool): The version flag value.
    """
    if value:
        typer.echo(f"WizardEye version: {_get_runtime_version()}")
        raise typer.Exit(code=0)


# -- Helper functions --


def _get_valid_tracks_from_exclude_tracks(
    exclude_tracks: Optional[List[str]],
    ref: str,
    kmer_length: Optional[int],
    offset_step: Optional[int],
    db_root: str,
    bwa_params: Optional[BWAParameters] = None,
    context: str = "filtering",
) -> List[str]:
    """Get valid tracks from exclude_tracks parameter."""
    if exclude_tracks is None:
        return []

    requested_tracks = from_charlist_to_list(exclude_tracks)
    available_tracks = get_tracks(
        ref_species=ref,
        kmer_length=kmer_length,
        offset_step=offset_step,
        db_root=db_root,
        bwa_params=bwa_params,
    )
    return resolve_requested_track_names(
        requested_tracks=requested_tracks,
        available_tracks=available_tracks,
        ref=ref,
        context=context,
    )


def _get_valid_tracks_from_exclude_tags(
    exclude_tags: Optional[List[str]],
    ref: str,
    kmer_length: Optional[int],
    offset_step: Optional[int],
    db_root: str,
    bwa_params: Optional[BWAParameters] = None,
) -> List[str]:
    """Get valid tracks from exclude_tags parameter."""
    if exclude_tags is None:
        return []

    requested_tags = from_charlist_to_list(exclude_tags)
    return from_tags_get_tracks(
        ref_species=ref,
        tags=requested_tags,
        kmer_length=kmer_length,
        offset_step=offset_step,
        db_root=db_root,
        bwa_params=bwa_params,
    )


def _request_tracks_from_args(
    ref: str,
    input_bam: str,
    db_root: str,
    exclude_tags: str,
    exclude_tracks: str,
    kmer_length: int,
    offset_step: int,
    only_unique: bool,
    no_cache: bool,
    cross_stringency: float,
    bwa_params: Optional[BWAParameters] = None,
) -> List[str]:
    """Validate CLI arguments and get corresponding tracks from database."""

    # Control input
    if ref is None:
        log("Reference (-r) must be specified.", "E")
        raise typer.Exit(code=1)

    if ref not in get_refs(db_root):
        log(f"Reference '{ref}' not found in database '{db_root}'.", "E")
        raise typer.Exit(code=1)

    if exclude_tags and exclude_tracks:
        log(
            "Cannot use both --exclude-tags and --exclude-tracks at the same time for filtering. Please choose one.",
            "E",
        )
        raise typer.Exit(code=1)

    if exclude_tags is None and exclude_tracks is None:
        log("Please provide either --exclude-tags or --exclude-tracks.", "E")
        raise typer.Exit(code=1)

    if exclude_tracks or exclude_tags:
        missing_params: List[str] = []
        if kmer_length is None:
            missing_params.append("-k")
        if offset_step is None:
            missing_params.append("-w")

        if missing_params:
            log(
                "For safety, all requested track generation parameters are mandatory. "
                f"Missing: {', '.join(missing_params)}",
                "E",
            )
            raise typer.Exit(code=1)

    # Log requested parameters
    colorful = sys.stdout.isatty() and sys.stderr.isatty()
    default = "\033[0;90m(default)\033[0m" if colorful else "(default)"

    log(f"Reference used: {ref}", "I")
    log(f"Input alignment file to process: {input_bam}", "I")
    if bwa_params:
        log(
            f"Associated BWA parameters: -n {bwa_params.missing_prob_err_rate}, -o {bwa_params.max_gap_opens}, -l {bwa_params.seed_length}, -t {bwa_params.threads}",
            "I",
        )
    log(
        f"Requested track generation parameters: -k {kmer_length}, -w {offset_step}",
        "I",
    )
    log(
        f"Mask source mode: {'uniquely aligned k-mers' if only_unique else f'all k-mers {default}'}",
        "I",
    )
    log(
        f"Mask cache mode: {'disabled (--no-cache)' if no_cache else f'enabled {default}'}",
        "I",
    )
    if cross_stringency:
        log(f"Requested mask stringency threshold: {cross_stringency}", "I")

    if exclude_tracks:
        valid_tracks = _get_valid_tracks_from_exclude_tracks(
            exclude_tracks=exclude_tracks,
            ref=ref,
            kmer_length=kmer_length,
            offset_step=offset_step,
            db_root=db_root,
            bwa_params=bwa_params,
            context="filtering",
        )

    if exclude_tags:
        requested_tags = from_charlist_to_list(exclude_tags)
        valid_tracks = _get_valid_tracks_from_exclude_tags(
            exclude_tags=exclude_tags,
            ref=ref,
            kmer_length=kmer_length,
            offset_step=offset_step,
            db_root=db_root,
            bwa_params=bwa_params,
        )

        log(f"Requested tags to filter out: {', '.join(requested_tags)}", "I")

    if not valid_tracks:
        log(
            "No valid tracks were found from tags and parameters. Export cannot proceed.",
            "E",
        )
        raise typer.Exit(code=1)

    excluded_tracks = "\n\t- ".join(valid_tracks)
    log(
        f"Following tracks will be used in the calculation:\n\t- {excluded_tracks}", "I"
    )

    return valid_tracks


# -- CLI commands --


@app.callback(invoke_without_command=True)
def common_options(
    version: Optional[bool] = typer.Option(
        None,
        "--version",
        help="Show WizardEye version and exit.",
        callback=_version_callback,
        is_eager=True,
    ),
):
    """Global CLI options shared by all commands."""
    return


# Database subcommands

@db_app.command(
    help="Initialize a new WizardEye database.",
    no_args_is_help=True,
)
def init(
    db_root: str = typer.Option(
        ..., "-d", "--db-root", help="Path to the database root directory."
    ),
):
    """Initialize a new WizardEye database."""
    print_version_message()
    try:
        init_db(db_root)
    except FileExistsError as e:
        log(str(e), "E")
        raise typer.Exit(code=1)
    raise typer.Exit(code=0)


@db_app.command(
    help="Print the full database catalogue.",
    no_args_is_help=True,
)
def catalogue(
    db_root: str = typer.Option(
        ..., "-d", "--db-root", help="Path to the database root directory."
    ),
):
    """Print the full database catalogue."""
    print_version_message()
    if not valid_database(db_root):
        log(f"To initialize the database, run 'database init {db_root}'", "E")
        raise typer.Exit(code=1)
    print_full_catalogue(db_root)
    raise typer.Exit(code=0)


@db_app.command(
    help="Migrate database tracks from old naming format (sus_scrofa_k35_w1_n1_o2_l1) to new format (_k_w_bwahash).",
    no_args_is_help=True,
)
def migrate(
    db_root: str = typer.Option(
        ..., "-d", "--db-root", help="Path to the database root directory."
    ),
    from_version: str = typer.Option(
        "0.1.2", "--from", help="Source version to migrate from."
    ),
    to_version: str = typer.Option(
        "0.1.3", "--to", help="Target version to migrate to."
    ),
    bwa_r_best_hits: int = typer.Option(
        2147483647,
        "-bR",
        help="bwa aln -R parameter value. When n_best_hits>=-bR, bwa aln do not explore suboptimal hits.",
    ),
    bwa_samse_n: int = typer.Option(
        2147483647,
        "-bsn",
        help="bwa samse -n parameter value. If n_kept_hits>-n, kept hits will NOT be saved in the track.",
    ),
):
    """Migrate database tracks from old naming format to new format with BWA parameters hash."""
    print_version_message()
    
    if not valid_database(db_root):
        log(f"To initialize the database, run 'database init {db_root}'", "E")
        raise typer.Exit(code=1)
    
    log(f"Starting migration from version {from_version} to {to_version}", "I")
    
    if from_version == "0.1.2" and to_version == "0.1.3":
        try:
            result = migrate_database(
                db_root=db_root,
                bwa_r_best_hits=bwa_r_best_hits,
                bwa_samse_n=bwa_samse_n,
            )
        except ValueError as e:
            log(str(e), "E")
            raise typer.Exit(code=1)
    else:
        log("Unsupported version migration.", "E")
        raise typer.Exit(code=1)
    
    # Display warnings
    for warning in result.get("warnings", []):
        log(warning, "W")
    
    # Display errors
    for error in result.get("errors", []):
        log(error, "E")
    
    migrated_count = result.get("migrated_count", 0)
    new_db_path = result.get("new_db_path")
    log(f"Migration completed. {migrated_count} track(s) migrated.", "S")
    log(f"Migrated database created at: {new_db_path}", "S")
    log("The original database has NOT been modified. You can now use the migrated version.", "I")
    
    raise typer.Exit(code=0)


@db_app.command(
    help="Delete all .bed files from the database (asks for confirmation).",
    no_args_is_help=True,
)
def clean(
    db_root: str = typer.Option(
        ..., "-d", "--db-root", help="Path to the database root directory."
    ),
    yes: bool = typer.Option(
        False, "-y", "--yes", help="Skip confirmation prompt."
    ),
):
    """Delete all .bed files from the database."""
    print_version_message()
    if not valid_database(db_root):
        log(f"To initialize the database, run 'database init {db_root}'", "E")
        raise typer.Exit(code=1)

    if not yes:
        log(
            f"About to delete the masks cached in '{db_root}'. This will only impact filtering time.",
            "W",
        )
        answer = typer.prompt("Enter 'yes' to confirm cleanup", default="No")
        if answer.lower() not in ["yes", "y"]:
            log("Cleanup cancelled.", "I")
            raise typer.Exit(code=0)

    clean_db(db_root=db_root)
    raise typer.Exit(code=0)


@update_app.command(
    name="track-tags",
    help="Replace tags for one existing track (requires full track parameters).",
    no_args_is_help=True,
)
def track_tags(
    db_root: str = typer.Option(
        ..., "-d", "--db-root", help="Path to the database root directory."
    ),
    ref: Optional[str] = typer.Option(
        None, "-r", "--ref", help="Reference identifier of the track to update."
    ),
    track: Optional[str] = typer.Option(
        None, "-q", "--track", help="Track query identifier to update."
    ),
    kmer_length: Optional[int] = typer.Option(
        None, "-k", "--kmer-length", help="Track k-mer length."
    ),
    offset_step: Optional[int] = typer.Option(
        None, "-w", "--offset-step", help="Track sliding window/offset step."
    ),
    # BWA aln parameters
    bwa_missing_prob_err_rate: Optional[float] = typer.Option(
        None, "-bn", help="bwa aln -n. Max diff. or missing prob. under 0.02 err rate."
    ),
    bwa_max_gap_opens: Optional[int] = typer.Option(
        None, "-bo", help="bwa aln -o. Maximum number or fraction of gap opens."
    ),
    bwa_seed_length: Optional[int] = typer.Option(
        None, "-bl", help="bwa aln -l. Seed length."
    ),
    bwa_all_aln: Optional[bool] = typer.Option(
        None,
        "-bN",
        help="bwa aln -N. Non-iterative mode: search for all n-difference hits.",
    ),
    bwa_threads: Optional[int] = typer.Option(
        None, "-bj", "--bwa-threads", help="bwa aln -j. Number of threads for bwa aln."
    ),
    bwa_r_best_hits: Optional[int] = typer.Option(
        None,
        "-bR",
        help="bwa aln -R. When n_best_hits>=-bR, bwa aln do not explore suboptimal hits.",
    ),
    bwa_samse_n: Optional[int] = typer.Option(
        None,
        "-bsn",
        help="bwa samse -n. If n_kept_hits>-n, kept hits will NOT be saved in the track.",
    ),
    new_tags: Optional[List[str]] = typer.Option(
        None,
        "-t",
        "--tag",
        "--tags",
        help="Replacement tags for the track (comma-separated).",
    ),
):
    """Replace tags for one existing track."""
    print_version_message()

    if not valid_database(db_root):
        log(f"To initialize the database, run 'database init {db_root}'", "E")
        raise typer.Exit(code=1)

    required_params = {
        "--db-root": db_root,
        "--ref": ref,
        "--track": track,
        "--kmer-length": kmer_length,
        "--offset-step": offset_step,
        "-bn (bwa_missing_prob_err_rate)": bwa_missing_prob_err_rate,
        "-bo (bwa_max_gap_opens)": bwa_max_gap_opens,
        "-bl (bwa_seed_length)": bwa_seed_length,
        "-bN (bwa_all_aln)": bwa_all_aln,
        "-bj (bwa_threads)": bwa_threads,
        "-bR (bwa_r_best_hits)": bwa_r_best_hits,
        "-bsn (bwa_samse_n)": bwa_samse_n,
        "--tags": new_tags,
    }
    missing = [name for name, value in required_params.items() if value is None]
    if missing:
        log(
            "All track parameters are required. Missing: " + ", ".join(missing),
            "E",
        )
        raise typer.Exit(code=1)

    parsed_new_tags = from_charlist_to_list(new_tags, lowercase=True)
    if not parsed_new_tags:
        log("No replacement tags provided. Use --tags with at least one tag.", "E")
        raise typer.Exit(code=1)

    # Create BWAParameters object with all parameters
    bwa_params = BWAParameters(
        missing_prob_err_rate=bwa_missing_prob_err_rate,
        max_gap_opens=bwa_max_gap_opens,
        seed_length=bwa_seed_length,
        all_aln=bwa_all_aln,
        threads=bwa_threads,
        r_best_hits=bwa_r_best_hits,
        samse_n=bwa_samse_n,
    )

    if not check_track_exists(
        ref_species=ref,
        query_species=track,
        kmer_length=kmer_length,
        offset_step=offset_step,
        db_root=db_root,
        bwa_params=bwa_params,
    ):
        log(
            f"Track does not exist for ref='{ref}', track='{track}', "
            f"k={kmer_length}, w={offset_step}, bwa_params={bwa_params}.",
            "E",
        )
        raise typer.Exit(code=1)

    try:
        updated_yaml, old_tags, updated_tags = update_track_tags(
            ref_species=ref,
            query_species=track,
            kmer_length=kmer_length,
            offset_step=offset_step,
            tags=parsed_new_tags,
            db_root=db_root,
            bwa_params=bwa_params,
        )
    except FileNotFoundError as e:
        log(str(e), "E")
        raise typer.Exit(code=1)

    old_tags_str = ", ".join(old_tags) if old_tags else "-"
    new_tags_str = ", ".join(updated_tags) if updated_tags else "-"
    log(f"Track tags updated in: {updated_yaml}", "S")
    log(f"Previous tags: {old_tags_str}", "I")
    log(f"New tags: {new_tags_str}", "I")
    raise typer.Exit(code=0)


# Main alignment command with all parameters for track generation
# Same behavior that generate_cross_mappability_filter_bwa.sh
@app.command(
    help="Generate a cross-mappability track and associated metadata from input FASTA files using specific bwa parameters.",
    no_args_is_help=True,
)
def align(
    # Input parameters
    db_root: str = typer.Option(
        ...,
        "-d",
        "--db-root",
        help="Path to the database root directory where to save the track.",
    ),
    input_fasta: Optional[List[str]] = typer.Option(
        None,
        "-i",
        help="Path to FASTA file to align on the target (query).",
    ),
    input_target: str = typer.Option(
        None, "-r", help="Path to the reference target FASTA."
    ),
    track_id: Optional[str] = typer.Option(
        None,
        "--track_ID",
        "--track-id",
        help="Manual identifier to store in track metadata.",
    ),
    tag: Optional[List[str]] = typer.Option(
        None,
        "--tag",
        "-t",
        help="Tag(s) to associate to the track. Comma-separated (e.g. tag1,tag2).",
    ),
    # Splitting parameters
    kmer_length: int = typer.Option(None, "-k", help="Length of k-mers to produce."),
    offset_step: int = typer.Option(
        1,
        "-w",
        "--offset-step",
        "--sliding-window",
        help="Offset/sliding window step for k-mers.",
    ),
    # Alignment parameters
    bwa_missing_prob_err_rate: float = typer.Option(
        0.01, "-bn", help="bwa aln -n. Max diff. or missing prob. under 0.02 err rate."
    ),
    bwa_max_gap_opens: int = typer.Option(
        2, "-bo", help="bwa aln -o. Maximum number or fraction of gap opens."
    ),
    bwa_seed_length: int = typer.Option(16500, "-bl", help="bwa aln -l. Seed length."),
    bwa_all_aln: bool = typer.Option(
        False,
        "-bN",
        help="bwa aln -N. Non-iterative mode: search for all n-difference hits (please read README.md before using it).",
    ),
    bwa_threads: int = typer.Option(
        1, "-bj", "--bwa-threads", help="bwa aln -j. Number of threads for bwa aln."
    ),
    bwa_r_best_hits: int = typer.Option(
        30,
        "-bR",
        help="bwa aln -R. When n_best_hits>=-bR, bwa aln do not explore suboptimal hits. Processed per score, i.e. does not affect best hits count.",
    ),
    bwa_samse_n: int = typer.Option(
        2000000000,
        "-bsn",
        help="bwa samse -n. If n_kept_hits>-n, kept hits will NOT be saved in the track.",
    ),
    # Parallelisation parameters
    jobs: int = typer.Option(
        1,
        "-j",
        "--jobs",
        help="Number of threads to parallelize bwa instances and data processing.",
    ),
    chunk_size: int = typer.Option(
        2000000,
        "-cs",
        "--chunk-size",
        help="Number of k-mers to align per bwa instances.",
    ),
    tmp_dir: Optional[str] = typer.Option(
        None,
        "--tmp-dir",
        help="Custom temporary directory path. If provided, overrides the default TMPDIR=/tmp location for WizardEye temporary files.",
    ),
    force: bool = typer.Option(
        None, help="Force track creation even if it already exists."
    ),
):
    print_version_message()
    start_time = time.perf_counter()

    # Validate BWA parameters to prevent overflow in bwa aln and bwa samse
    max_bwa_int = 2147483647
    if bwa_r_best_hits >= max_bwa_int:
        log(
            f"Error: -bR value ({bwa_r_best_hits}) must be less than {max_bwa_int} (BWA maximum integer value).",
            "E",
        )
        raise typer.Exit(code=1)
    if bwa_samse_n >= max_bwa_int:
        log(
            f"Error: -bsn value ({bwa_samse_n}) must be less than {max_bwa_int} (BWA maximum integer value).",
            "E",
        )
        raise typer.Exit(code=1)

    if bwa_all_aln:
        log(
            "-N parameter used. Results may still not be exhaustive due to how bwa works. See README for more details.",
            "W",
        )
        log("-R set to 2147483647 to map bwa internal behaviour.", "W")
        bwa_r_best_hits = 2147483647

    if input_fasta and input_target and kmer_length:
        if not valid_database(db_root):
            log(f"To initialize the database, run 'database --init {db_root}'", "E")
            raise typer.Exit(code=1)

        input_fastas = from_charlist_to_list(input_fasta)
        if not input_fastas:
            log("No input FASTA provided after parsing -i values.", "E")
            raise typer.Exit(code=1)

        # Validate thread configuration
        try:
            available_cpus = os.cpu_count() or 1
            total_threads = bwa_threads * jobs
            if total_threads > available_cpus:
                log(
                    f"Error: -bj {bwa_threads} * -j {jobs} = {total_threads} threads "
                    f"exceeds available CPUs ({available_cpus}). "
                    f"Please reduce the number of threads.",
                    "E",
                )
                raise typer.Exit(code=1)
        except Exception as e:
            log(f"Could not determine CPU count: {e}", "W")

        manual_track_id: Optional[str] = None
        if track_id is not None:
            manual_track_id = track_id.strip()
            if not manual_track_id:
                log("--track_ID cannot be empty.", "E")
                raise typer.Exit(code=1)

        bwa_params = BWAParameters(
            missing_prob_err_rate=bwa_missing_prob_err_rate,
            max_gap_opens=bwa_max_gap_opens,
            seed_length=bwa_seed_length,
            all_aln=bwa_all_aln,
            threads=bwa_threads,
            r_best_hits=bwa_r_best_hits,
            samse_n=bwa_samse_n,
        )

        total_inputs = len(input_fastas)
        for idx, one_input_fasta in enumerate(input_fastas, start=1):
            base_query_name = Path(one_input_fasta).stem
            query_track_id = base_query_name

            bwa_hash = get_bwa_params_hash(bwa_params)
            track_name = f"{query_track_id}_k{kmer_length}_w{offset_step}_bwa{bwa_hash}"
            log(
                f"[{idx}/{total_inputs}] Processing input FASTA: {one_input_fasta}", "I"
            )

            if (
                check_track_exists(
                    Path(input_target).stem,
                    query_track_id,
                    kmer_length,
                    offset_step,
                    db_root,
                    bwa_params=bwa_params,
                )
                and not force
            ):
                log(
                    f"Track '{track_name}' already exists for reference '{Path(input_target).stem}', skipping track creation.",
                    "W",
                )
                continue

            try:
                log(
                    f"Creating track for reference '{Path(input_target).stem}' and query '{query_track_id}'...",
                    "I",
                )
                create_mappability_track(
                    input_fasta=one_input_fasta,
                    input_target=input_target,
                    track_id=query_track_id,
                    manual_track_id=manual_track_id,
                    kmer_length=kmer_length,
                    offset_step=offset_step,
                    chunk_size=chunk_size,
                    n_threads=jobs,
                    bwa_params=bwa_params,
                    db_root=db_root,
                    tags=tag,
                    tmp_dir_custom=tmp_dir,
                )
            except subprocess.CalledProcessError as e:
                log(f"Subprocess failed: {e}", "E")
                raise typer.Exit(code=2)
            except ValueError as e:
                log(str(e), "E")
                raise typer.Exit(code=2)

        elapsed = time.perf_counter() - start_time
        log(f"Alignment completed in {elapsed:.2f}s", "S")
        raise typer.Exit(code=0)

    log(
        "align requires a FASTA to split and align on a reference and the associated parameters (-i, -r, -k)",
        "E",
    )
    raise typer.Exit(code=1)


@app.command(
    help="Filter an input BAM using selected cross-mappability tracks and stringency.",
    no_args_is_help=True,
)
def filter(
    # Alignment parameters used to construct the BAM file
    input_bam=typer.Option(
        None, "-i", "--input", help="Path to the BAM file to filter."
    ),
    ref: Optional[str] = typer.Option(None, "-r", help="Reference used for alignment."),
    bwa_missing_prob_err_rate: Optional[float] = typer.Option(
        None, "-bn", help="bwa aln -n used for alignment."
    ),
    bwa_max_gap_opens: Optional[int] = typer.Option(
        None, "-bo", help="bwa aln -o used for alignment."
    ),
    bwa_seed_length: Optional[int] = typer.Option(
        None, "-bl", help="bwa aln -l used for alignment."
    ),
    bwa_all_aln: Optional[bool] = typer.Option(
        None,
        "-bN",
        help="bwa aln -N used for alignment (compute every alternative mapping).",
    ),
    bwa_r_best_hits: int = typer.Option(
        30, "-bR", help="bwa aln -R used for alignment."
    ),
    bwa_samse_n: int = typer.Option(
        2000000000, "-bsn", help="bwa samse -n used for alignment."
    ),
    db_root: str = typer.Option(
        ..., "-d", "--db-root", help="Path to the database root directory."
    ),
    # Filtration parameters
    exclude_tags: Optional[List[str]] = typer.Option(
        None,
        "--exclude-tags",
        help="Tag(s) to filter out. Comma-separated (e.g. tag1,tag2).",
    ),
    exclude_tracks: Optional[List[str]] = typer.Option(
        None,
        "--exclude-tracks",
        help="Track identifier(s) to filter out. Comma-separated (e.g. genius_species1, genius_species2).",
    ),
    cross_stringency: float = typer.Option(
        0.99,
        "-p",
        "--stringency",
        "-rc",
        "--cross-stringency",
        help="Stringency threshold in [0.0, 1.0].",
    ),
    only_unique: bool = typer.Option(
        False,
        "--only-unique",
        "--only_unique",
        help="Consider only unique k-mers (no XA tag & MAPQ>0) area.",
    ),
    output_filtered_bam: Optional[str] = typer.Option(
        None,
        "-o",
        "--output",
        help="Output BAM for reads kept after filtering (non-overlapping mask).",
    ),
    output_excluded_bam: Optional[str] = typer.Option(
        None,
        "--excluded-output",
        help="Output BAM for reads excluded by the generated mask.",
    ),
    output_report_tsv: Optional[str] = typer.Option(
        None,
        "--report-output",
        help="Output TSV report with columns: read_id, excluded, overlapped, tags.",
    ),
    export_bam: bool = typer.Option(
        False,
        "--export-bam",
        help="Write filtered/excluded BAM outputs. By default, only the TSV report is generated.",
    ),
    kmer_length: Optional[int] = typer.Option(
        None,
        "-k",
        help="K-mer length to filter on. Must match track generation parameter.",
    ),
    offset_step: Optional[int] = typer.Option(
        None,
        "-w",
        "--offset-step",
        "--sliding-window",
        help="Offset/sliding window to filter on. Must match track generation parameter.",
    ),
):
    start_time = time.perf_counter()
    print_version_message()

    if not input_bam:
        log("Input BAM (-i/--input) must be specified.", "E")
        raise typer.Exit(code=1)

    bwa_params = BWAParameters(
        missing_prob_err_rate=bwa_missing_prob_err_rate,
        max_gap_opens=bwa_max_gap_opens,
        seed_length=bwa_seed_length,
        all_aln=bwa_all_aln,
        threads=None,
        r_best_hits=bwa_r_best_hits,
        samse_n=bwa_samse_n,
    )
    valid_tracks = _request_tracks_from_args(
        ref,
        input_bam,
        db_root,
        exclude_tags,
        exclude_tracks,
        kmer_length,
        offset_step,
        only_unique,
        None,
        cross_stringency,
        bwa_params,
    )
    print("-" * 80)
    log("Starting filtration...", "I")

    if not input_bam:
        log("Input BAM (-i/--input) must be specified.", "E")
        raise typer.Exit(code=1)

    try:
        filter_result = filter_bam(
            input_bam=input_bam,
            ref=ref,
            db_root=db_root,
            exclude_tracks=valid_tracks,
            kmer_length=kmer_length,
            offset_step=offset_step,
            bwa_params=bwa_params,
            stringency=cross_stringency,
            consider_all=not (only_unique),
            output_filtered_bam=output_filtered_bam,
            output_excluded_bam=output_excluded_bam,
            output_report_tsv=output_report_tsv,
            export_bam=export_bam,
        )
    except FileNotFoundError as e:
        log(str(e), "E")
        raise typer.Exit(code=1)
    except ValueError as e:
        log(str(e), "E")
        raise typer.Exit(code=1)
    except RuntimeError as e:
        log(str(e), "E")
        raise typer.Exit(code=2)
    except subprocess.CalledProcessError as e:
        log(f"Subprocess failed: {e}", "E")
        raise typer.Exit(code=2)

    elapsed = time.perf_counter() - start_time
    print("-" * 80)

    log(
        f"Filtration completed in {elapsed:.2f}s:\n\
		Total reads before filtration: {filter_result['n_total']}\n\
		Total reads after filtration: {filter_result['n_filtered']}\n\
		Total reads excluded: {filter_result['n_excluded']}",
        "S",
    )

    log(f"Read exclusion report saved at {filter_result['report_tsv']}", "S")
    if (
        filter_result["filtered_bam"] is not None
        and filter_result["excluded_bam"] is not None
    ):
        log(f"Filtered BAM (kept reads): {filter_result['filtered_bam']}", "S")
        log(f"Excluded BAM (masked reads): {filter_result['excluded_bam']}", "S")

    log("Thank you for using WizardEye!", "S")

    raise typer.Exit(code=0)


@app.command(
    help="Export a merged BED mask using the same track-selection logic as filter.",
    no_args_is_help=True,
)
def export(
    ref: str = typer.Option(..., "-r", "--ref", help="Reference species name."),
    # Same selector interface as `filter`.
    exclude_tags: Optional[List[str]] = typer.Option(
        None,
        "--exclude-tags",
        "--exclude_tags",
        help="Tag(s) to export from. Comma-separated (e.g. tag1,tag2).",
    ),
    exclude_tracks: Optional[List[str]] = typer.Option(
        None,
        "--exclude-tracks",
        "--exclude_tracks",
        help="Track identifier(s) to export from. Comma-separated (e.g. genius_species1, genius_species2).",
    ),
    kmer_length: Optional[int] = typer.Option(None, "-k", help="K-mer size to target."),
    offset_step: Optional[int] = typer.Option(
        None,
        "-w",
        "--offset-step",
        "--sliding-window",
        help="Offset/sliding window to target.",
    ),
    bwa_missing_prob_err_rate: Optional[float] = typer.Option(
        None, "-bn", help="bwa aln -n used to generate selected tracks."
    ),
    bwa_max_gap_opens: Optional[int] = typer.Option(
        None, "-bo", help="bwa aln -o used to generate selected tracks."
    ),
    bwa_seed_length: Optional[int] = typer.Option(
        None, "-bl", help="bwa aln -l used to generate selected tracks."
    ),
    bwa_all_aln: Optional[bool] = typer.Option(
        None,
        "-bN",
        help="bwa aln -N used to generate selected tracks (compute every alternative mapping).",
    ),
    bwa_r_best_hits: int = typer.Option(
        30, "-bR", help="bwa aln -R used to generate selected tracks."
    ),
    bwa_samse_n: int = typer.Option(
        2000000000, "-bsn", help="bwa samse -n used to generate selected tracks."
    ),
    cross_stringency: float = typer.Option(
        0.99,
        "-p",
        "--stringency",
        "-rc",
        "--cross-stringency",
        help="Stringency threshold in [0.0, 1.0].",
    ),
    only_unique: bool = typer.Option(
        False,
        "--only-unique",
        "--only_unique",
        help="Consider only unique k-mers (no XA tag & MAPQ>0) area.",
    ),
    db_root: str = typer.Option(
        ..., "-d", "--db-root", help="Path to the database root directory."
    ),
    output_bed: Optional[str] = typer.Option(
        None,
        "-o",
        "--output",
        help="Output BED path for the merged all-overlaps mask.",
    ),
    n_threads: int = typer.Option(
        1,
        "-j",
        help="Number of threads for parallel overlap extraction from BigWig tracks.",
    ),
):
    """Export a merged BED mask using the same track-selection logic as `filter`."""
    print_version_message()

    bwa_params = BWAParameters(
        missing_prob_err_rate=bwa_missing_prob_err_rate,
        max_gap_opens=bwa_max_gap_opens,
        seed_length=bwa_seed_length,
        all_aln=bwa_all_aln,
        threads=None,
        r_best_hits=bwa_r_best_hits,
        samse_n=bwa_samse_n,
    )

    valid_tracks = _request_tracks_from_args(
        ref,
        None,
        db_root,
        exclude_tags,
        exclude_tracks,
        kmer_length,
        offset_step,
        only_unique,
        None,
        cross_stringency,
        bwa_params,
    )

    try:
        merged_bed = generate_global_mask(
            ref_species=ref,
            inputs=valid_tracks,
            kmer_length=kmer_length,
            offset_step=offset_step,
            cross_stringency=cross_stringency,
            consider_all=not (only_unique),
            n_threads=n_threads,
            db_root=db_root,
            output_file=output_bed,
            bwa_params=bwa_params,
        )
    except ValueError as e:
        log(str(e), "E")
        raise typer.Exit(code=1)
    except RuntimeError as e:
        log(str(e), "E")
        raise typer.Exit(code=2)
    except subprocess.CalledProcessError as e:
        log(f"Subprocess failed: {e}", "E")
        raise typer.Exit(code=2)

    log(f"Export completed. Merged file: {merged_bed}", "S")


@app.command(
    "import",
    help="Import a track manually by providing BigWig files and generation parameters.",
    no_args_is_help=True,
)
def import_tracks(
    ref: str = typer.Option(
        ..., "-r", "--ref", help="Reference species/target name in the database."
    ),
    query: str = typer.Option(
        ..., "-i", "--input", help="Query/input species name for the track."
    ),
    kmer_length: int = typer.Option(
        ..., "-k", help="K-mer size used to generate the imported track."
    ),
    offset_step: int = typer.Option(
        ...,
        "-w",
        "--offset-step",
        "--sliding-window",
        help="Offset/sliding window used to generate the imported track.",
    ),
    map_all_bw: str = typer.Option(
        ..., "-ma", "--map-all-bw", help="Path to map_all.bw generated externally."
    ),
    map_uniq_bw: str = typer.Option(
        ..., "-mu", "--map-uniq-bw", help="Path to map_uniq.bw generated externally."
    ),
    input_fasta: Optional[str] = typer.Option(
        None,
        "--input-fasta",
        help="Original input FASTA path used to generate the imported track.",
    ),
    reference_fasta: Optional[str] = typer.Option(
        None,
        "--reference-fasta",
        help="Original reference FASTA path used to generate the imported track.",
    ),
    reference_fasta_md5: Optional[str] = typer.Option(
        None,
        "--reference-fasta-md5",
        help="Reference FASTA MD5 used during generation (validated if --reference-fasta is provided).",
    ),
    mapping_tool: str = typer.Option(
        "bwa aln",
        "--mapping-tool",
        help="Mapping tool used to generate imported files (e.g. 'bwa aln').",
    ),
    bwa_missing_prob_err_rate: float = typer.Option(
        0.01,
        "-bn",
        help="bwa aln -n used for alignment.",
    ),
    bwa_max_gap_opens: int = typer.Option(
        2,
        "-bo",
        help="bwa aln -o used for alignment.",
    ),
    bwa_seed_length: int = typer.Option(
        16500,
        "-bl",
        help="bwa aln -l used for alignment.",
    ),
    bwa_all_aln: Optional[bool] = typer.Option(
        None,
        "-bN",
        help="bwa aln -N used for alignment.",
    ),
    bwa_r_best_hits: int = typer.Option(
        30,
        "-bR",
        help="bwa aln -R used for alignment.",
    ),
    bwa_samse_n: int = typer.Option(
        2000000000,
        "-bsn",
        help="bwa samse -n used after alignment.",
    ),
    n_threads: int = typer.Option(
        1,
        "-j",
        help="Thread count used for alignment.",
    ),
    db_root: str = typer.Option(
        ..., "-d", "--db-root", help="Path to the database root directory."
    ),
    tag: Optional[List[str]] = typer.Option(
        None,
        "--tag",
        "-t",
        help="Tag(s) to associate to the imported track. Repeatable and comma-separated.",
    ),
    force: bool = typer.Option(
        False,
        help="Overwrite imported track files/metadata if the track already exists.",
    ),
):
    print_version_message()
    """Import a track manually by providing BigWig files and generation parameters."""
    if not valid_database(db_root):
        log(f"To initialize the database, run 'database --init {db_root}'", "E")
        raise typer.Exit(code=1)

    bwa_params = BWAParameters(
        missing_prob_err_rate=bwa_missing_prob_err_rate,
        max_gap_opens=bwa_max_gap_opens,
        seed_length=bwa_seed_length,
        all_aln=bwa_all_aln,
        r_best_hits=bwa_r_best_hits,
        samse_n=bwa_samse_n,
        threads=1,
    )

    log("Importing track to database...", "I")
    try:
        result = import_track(
            ref_species=ref,
            query_species=query,
            kmer_length=kmer_length,
            offset_step=offset_step,
            map_all_bw=map_all_bw,
            map_uniq_bw=map_uniq_bw,
            db_root=db_root,
            input_fasta=input_fasta,
            reference_fasta=reference_fasta,
            reference_fasta_md5=reference_fasta_md5,
            tags=tag,
            mapping_tool=mapping_tool,
            bwa_params=bwa_params,
            n_threads=n_threads,
            force=force,
        )
    except FileExistsError as e:
        log(str(e), "E")
        raise typer.Exit(code=1)
    except FileNotFoundError as e:
        log(str(e), "E")
        raise typer.Exit(code=1)
    except ValueError as e:
        log(str(e), "E")
        raise typer.Exit(code=1)

    log(f"Track imported in: {result['track_dir']}", "S")
    log(f"param.yaml written to: {result['param_yaml']}", "S")
    log(f"map_all BigWig stored at: {result['map_all_bw']}", "S")
    log(f"map_uniq BigWig stored at: {result['map_uniq_bw']}", "S")
    raise typer.Exit(code=0)


@app.command(
    help="Count for each read the number of overlapping k-mers fron required tracks.",
    no_args_is_help=True,
)
def count(
    # Alignment parameters used to construct the BAM file
    input_bam=typer.Option(
        None, "-i", "--input", help="Path to the BAM file to filter."
    ),
    ref: Optional[str] = typer.Option(None, "-r", help="Reference used for alignment."),
    count_mode: Optional[str] = typer.Option(
        "mean", "-m", "--mode", help="Count mode (mean, std, max, min, cov or sum)"
    ),
    bwa_missing_prob_err_rate: Optional[float] = typer.Option(
        None, "-bn", help="bwa aln -n used for alignment."
    ),
    bwa_max_gap_opens: Optional[int] = typer.Option(
        None, "-bo", help="bwa aln -o used for alignment."
    ),
    bwa_seed_length: Optional[int] = typer.Option(
        None, "-bl", help="bwa aln -l used for alignment."
    ),
    bwa_all_aln: Optional[bool] = typer.Option(
        None,
        "-bN",
        help="bwa aln -N used for alignment (compute every alternative mapping).",
    ),
    bwa_r_best_hits: int = typer.Option(
        30, "-bR", help="bwa aln -R used for alignment."
    ),
    bwa_samse_n: int = typer.Option(
        2000000000, "-bsn", help="bwa samse -n used for alignment."
    ),
    # Track generation parameters
    db_root: str = typer.Option(
        ..., "-d", "--db-root", help="Path to the database root directory."
    ),
    kmer_length: Optional[int] = typer.Option(
        None,
        "-k",
        help="K-mer length to filter on. Must match track generation parameter.",
    ),
    offset_step: Optional[int] = typer.Option(
        None,
        "-w",
        "--offset-step",
        "--sliding-window",
        help="Offset/sliding window to filter on. Must match track generation parameter.",
    ),
    # Filtration parameters
    exclude_tags: Optional[List[str]] = typer.Option(
        None,
        "--exclude-tags",
        help="Tag(s) to filter out. Comma-separated (e.g. tag1,tag2).",
    ),
    exclude_tracks: Optional[List[str]] = typer.Option(
        None,
        "--exclude-tracks",
        help="Track identifier(s) to filter out. Comma-separated (e.g. genius_species1, genius_species2).",
    ),
    only_unique: bool = typer.Option(
        False,
        "--only-unique",
        "--only_unique",
        help="Consider only unique k-mers (no XA tag & MAPQ>0) area.",
    ),
    no_cache: bool = typer.Option(
        False,
        "--no-cache",
        help="Disable mask caching: recompute track masks and avoid writing cache files.",
    ),
    output_report_tsv: Optional[str] = typer.Option(
        None,
        "--report-output",
        help="Output TSV report with columns: read_id, excluded, overlapped, tags.",
    ),
    n_threads: int = typer.Option(
        1,
        "-j",
        help="Number of threads for parallel overlap extraction from BigWig tracks.",
    ),
):

    start_time = time.perf_counter()

    if not input_bam:
        log("Input BAM (-i/--input) must be specified.", "E")
        raise typer.Exit(code=1)

    if count_mode not in ["mean", "std", "max", "min", "cov", "sum"]:
        log(
            f"Invalid count mode {count_mode}. Available: mean, std, max, min, cov or sum."
        )
        raise typer.Exit(1)

    print_version_message()

    bwa_params = BWAParameters(
        missing_prob_err_rate=bwa_missing_prob_err_rate,
        max_gap_opens=bwa_max_gap_opens,
        seed_length=bwa_seed_length,
        all_aln=bwa_all_aln,
        threads=None,
        r_best_hits=bwa_r_best_hits,
        samse_n=bwa_samse_n,
    )

    valid_tracks = _request_tracks_from_args(
        ref,
        input_bam,
        db_root,
        exclude_tags,
        exclude_tracks,
        kmer_length,
        offset_step,
        only_unique,
        no_cache,
        None,
        bwa_params,
    )

    print("-" * 80)
    log(f"Starting computing k-mers {count_mode} per read...", "I")

    try:
        count_result = count_k_mers_on_bam(
            input_bam=input_bam,
            ref=ref,
            db_root=db_root,
            exclude_tracks=valid_tracks,
            kmer_length=kmer_length,
            offset_step=offset_step,
            count_mode=count_mode,
            bwa_params=bwa_params,
            consider_all=not (only_unique),
            n_threads=n_threads,
            output_report_tsv=output_report_tsv,
        )

    except FileNotFoundError as e:
        log(str(e), "E")
        raise typer.Exit(code=1)
    except ValueError as e:
        log(str(e), "E")
        raise typer.Exit(code=1)
    except RuntimeError as e:
        log(str(e), "E")
        raise typer.Exit(code=2)
    except subprocess.CalledProcessError as e:
        log(f"Subprocess failed: {e}", "E")
        raise typer.Exit(code=2)

    elapsed = time.perf_counter() - start_time
    print("-" * 80)

    log(
        f"Count and stats computation completed in {elapsed:.2f}s:\n\
		BAM records processed: {count_result['n_processed']}/{count_result['n_total_records']}\n",
        "S",
    )
    log(f"Count report saved at {count_result['report_tsv']}", "S")
    log("Thank you for using WizardEye!", "S")

    raise typer.Exit(code=0)


# -- Entry point --

if __name__ == "__main__":
    app()
