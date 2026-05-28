# -*- coding: utf-8 -*-

"""Database management for WizardEye tracks.

This module provides functions to initialize and validate the database structure, manage track metadata,
and import externally generated tracks. It defines the Track class to represent tracks with their
parameters and metadata, and includes utilities to query available tracks and display the catalogue of tracks
for each reference species.
"""

import shutil
import yaml

from dataclasses import dataclass, field
from pathlib import Path
from datetime import datetime
from importlib import metadata
from typing import Dict, List, Optional, Tuple, Union

from .utils import (
    log,
    get_name_from_param,
    from_charlist_to_list,
    file_md5,
    BWAParameters,
    get_bwa_params_hash,
)
from .version import PACKAGE_VERSION


# -- Database management functions --
def init_db(base_dir: Union[str, Path] = ".") -> Path:
    """Create /database and its info.yaml file.

    If the database already exists, it will not be modified and the existing info.yaml path will be returned.

    Args:
            base_dir: Base directory where the /database folder will be created. Defaults to current directory.

    Returns:
            Path to the info.yaml file in the initialized database.
    """
    root = Path(base_dir)
    db_dir = root / "database"
    db_dir.mkdir(parents=True, exist_ok=True)

    db_yaml_path = db_dir / "info.yaml"
    now = datetime.now().isoformat(timespec="seconds")
    created_date = now

    if db_yaml_path.exists():
        log(f"Database {db_dir} already exists, skipping creation.", "W")
        return db_yaml_path

    try:
        cross_tool_version = metadata.version("cross_tool")
    except metadata.PackageNotFoundError:
        try:
            cross_tool_version = metadata.version("wizardeye")
        except metadata.PackageNotFoundError:
            cross_tool_version = PACKAGE_VERSION

    db_content = {
        "cross_tool_version": cross_tool_version,
        "created_date": created_date,
        "last_updated": now,
    }

    with db_yaml_path.open("w", encoding="utf-8") as handle:
        yaml.safe_dump(db_content, handle, sort_keys=False)

    log(f"Database initialized: {db_yaml_path}", "S")
    return db_yaml_path


def valid_database(db_root: Union[str, Path]) -> bool:
    """Check if the database root contains a valid info.yaml file."""
    db_root = Path(db_root)
    db_yaml_path = db_root / "info.yaml"
    if not db_yaml_path.exists():
        log(f"Database info.yaml not found at {db_yaml_path}", "E")
        return False
    try:
        with db_yaml_path.open("r", encoding="utf-8") as handle:
            data = yaml.safe_load(handle)
        if not isinstance(data, dict) or "cross_tool_version" not in data:
            log(f"Invalid database info.yaml format at {db_yaml_path}", "E")
            return False
    except Exception as e:
        log(f"Failed to read database info.yaml at {db_yaml_path}: {e}", "E")
        return False
    return True


def clean_db(db_root: Path) -> int:
    """Delete all BED files inside one WizardEye database root."

    Args:
            db_root(Path): Path to the database root.

    Returns:
            int: Number of bed files deleted.
    """
    bed_files = sorted(path for path in Path(db_root).rglob("*.bed") if path.is_file())

    if not bed_files:
        log(f"No BED file found under database'{db_root}'.", "I")
        return 0

    deleted = 0
    total_bytes_deleted = 0
    for bed_path in bed_files:
        total_bytes_deleted += bed_path.stat().st_size
        bed_path.unlink()
        deleted += 1

    log(f"Deleted {deleted} BED file(s) from database '{db_root}'.", "S")
    log(
        f"Disk space saved: {total_bytes_deleted / (1024 * 1024):.2f} Mo "
        f"({total_bytes_deleted} bytes).",
        "S",
    )
    return deleted


# -- Track classes and functions --


@dataclass(frozen=True)
class TrackParameters:
    """Alignment and generation parameters that define one track version."""

    kmer_length: int
    offset_step: int
    bwa_params: BWAParameters = field(default_factory=BWAParameters)

    def __post_init__(self) -> None:
        if self.kmer_length < 1:
            raise ValueError("kmer_length must be a positive integer")
        if self.offset_step < 1:
            raise ValueError("offset_step must be a positive integer")

    @property
    def n_value(self) -> float:
        return float(self.bwa_params.missing_prob_err_rate)


@dataclass(frozen=True)
class TrackIdentity:
    """Track identity independent from file system state."""

    ref_species: str
    query_species: str
    parameters: TrackParameters

    @property
    def name(self) -> str:
        bwa_hash = get_bwa_params_hash(self.parameters.bwa_params)
        return f"{self.query_species}_k{self.parameters.kmer_length}_w{self.parameters.offset_step}_bwa{bwa_hash}"


@dataclass
class Track:
    """Materialized track object with both path and loaded metadata."""

    db_root: Path
    identity: TrackIdentity
    info: Dict

    @property
    def track_dir(self) -> Path:
        return self.db_root / self.identity.ref_species / self.identity.name

    @property
    def track_name(self) -> str:
        return self.identity.name

    @property
    def param_yaml(self) -> Path:
        return self.track_dir / "param.yaml"

    @property
    def map_all_bw(self) -> Path:
        return self.track_dir / "map_all.bw"

    @property
    def map_uniq_bw(self) -> Path:
        return self.track_dir / "map_uniq.bw"

    @property
    def query_name(self) -> str:
        query_name = get_name_from_param(self.info)
        return str(query_name)

    def exists(self) -> bool:
        return self.map_all_bw.exists() and self.map_uniq_bw.exists()

    def to_dict(self) -> Dict:
        return {
            "track_dir": self.track_dir,
            "track_name": self.track_name,
            "query_name": self.query_name,
            "param": self.info,
        }

    @classmethod
    def from_param(
        cls,
        db_root: Union[str, Path],
        ref_species: str,
        query_species: str,
        kmer_length: int,
        offset_step: int,
        bwa_params: BWAParameters = BWAParameters(),
        info: Optional[Dict] = None,
    ) -> "Track":
        params = TrackParameters(
            kmer_length=kmer_length,
            offset_step=offset_step,
            bwa_params=bwa_params,
        )
        identity = TrackIdentity(
            ref_species=ref_species, query_species=query_species, parameters=params
        )
        return cls(db_root=Path(db_root), identity=identity, info=info or {})

    @classmethod
    def from_param_yaml(
        cls, db_root: Union[str, Path], ref_species: str, track_dir: Path
    ) -> Optional["Track"]:
        param_yaml = track_dir / "param.yaml"
        if not param_yaml.exists():
            return None

        with param_yaml.open("r", encoding="utf-8") as handle:
            param_content = yaml.safe_load(handle) or {}
        if not isinstance(param_content, dict):
            return None

        kmer_size = param_content.get("kmer_size")
        sliding_window = param_content.get("sliding_window")
        if kmer_size is None or sliding_window is None:
            return None

        try:
            kmer_size = int(kmer_size)
            sliding_window = int(sliding_window)
        except (TypeError, ValueError):
            return None

        bwa_params = param_content.get("bwa_parameters", {})
        if not isinstance(bwa_params, dict):
            bwa_params = {}

        try:
            bwa_n = float(bwa_params.get("-n", 0.01))
            bwa_o = int(bwa_params.get("-o", 2))
            bwa_l = int(bwa_params.get("-l", 16500))
            bwa_all_aln = bool(bwa_params.get("-N", False))
            bwa_threads = int(bwa_params.get("-t", 1))
            bwa_r = int(bwa_params.get("-R", 30))
            bwa_sn = int(bwa_params.get("-sn", 2000000000))
        except (TypeError, ValueError):
            bwa_n, bwa_o, bwa_l, bwa_all_aln, bwa_threads, bwa_r, bwa_sn = (
                0.01,
                2,
                16500,
                False,
                1,
                30,
                2000000000,
            )

        query_species = track_dir.name.split("_k", 1)[0]
        return cls.from_param(
            db_root=db_root,
            ref_species=ref_species,
            query_species=query_species,
            kmer_length=kmer_size,
            offset_step=sliding_window,
            bwa_params=BWAParameters(
                missing_prob_err_rate=bwa_n,
                max_gap_opens=bwa_o,
                seed_length=bwa_l,
                all_aln=bwa_all_aln,
                threads=bwa_threads,
                r_best_hits=bwa_r,
                samse_n=bwa_sn,
            ),
            info=param_content,
        )

    def load_info(self) -> Dict:
        if not self.param_yaml.exists():
            raise FileNotFoundError(f"Track YAML not found at {self.param_yaml}")

        with self.param_yaml.open("r", encoding="utf-8") as handle:
            content = yaml.safe_load(handle) or {}
        if not isinstance(content, dict):
            content = {}
        self.info = content
        return self.info

    def save_info(self) -> Path:
        with self.param_yaml.open("w", encoding="utf-8") as handle:
            yaml.safe_dump(self.info, handle, sort_keys=False)
        return self.param_yaml

    def update_tags(
        self, tags: Optional[List[str]]
    ) -> Tuple[Path, List[str], List[str]]:
        content = self.load_info()

        new_tags = from_charlist_to_list(tags, lowercase=True)
        existing_tags = content.get("tags", [])
        if not isinstance(existing_tags, list):
            existing_tags = []
        old_tags = from_charlist_to_list(
            [str(tag) for tag in existing_tags], lowercase=True
        )

        content["tags"] = new_tags
        content["last_updated"] = datetime.now().isoformat(timespec="seconds")
        self.info = content
        return self.save_info(), old_tags, new_tags


def check_track_exists(
    ref_species: str,
    query_species: str,
    kmer_length: int,
    offset_step: int,
    db_root: Union[str, Path],
    bwa_params: BWAParameters = BWAParameters(),
) -> bool:
    """Check if a track exists for the given reference/query and k/s parameters."""
    track = Track.from_param(
        db_root=db_root,
        ref_species=ref_species,
        query_species=query_species,
        kmer_length=kmer_length,
        offset_step=offset_step,
        bwa_params=bwa_params,
    )
    return track.exists()


def update_track_tags(
    ref_species: str,
    query_species: str,
    kmer_length: int,
    offset_step: int,
    tags: Optional[List[str]],
    db_root: Union[str, Path],
    bwa_params: BWAParameters = BWAParameters(),
) -> Tuple[Path, List[str], List[str]]:
    """Replace all tags in the track param.yaml and return old/new values."""
    track = Track.from_param(
        db_root=db_root,
        ref_species=ref_species,
        query_species=query_species,
        kmer_length=kmer_length,
        offset_step=offset_step,
        bwa_params=bwa_params,
    )
    return track.update_tags(tags)


def import_track(
    ref_species: str,
    query_species: str,
    kmer_length: int,
    offset_step: int,
    map_all_bw: Union[str, Path],
    map_uniq_bw: Union[str, Path],
    db_root: Union[str, Path],
    input_fasta: Optional[Union[str, Path]] = None,
    reference_fasta: Optional[Union[str, Path]] = None,
    reference_fasta_md5: Optional[str] = None,
    tags: Optional[List[str]] = None,
    mapping_tool: str = "bwa aln",
    bwa_params: BWAParameters = BWAParameters(),
    force: bool = False,
) -> Dict[str, Path]:
    """Import an externally generated track by copying BigWig files and writing param.yaml."""
    db_root = Path(db_root)

    all_bw_src = Path(map_all_bw)
    uniq_bw_src = Path(map_uniq_bw)
    if not all_bw_src.exists() or not all_bw_src.is_file():
        raise FileNotFoundError(f"map_all.bw source file not found: {all_bw_src}")
    if not uniq_bw_src.exists() or not uniq_bw_src.is_file():
        raise FileNotFoundError(f"map_uniq.bw source file not found: {uniq_bw_src}")

    track = Track.from_param(
        db_root=db_root,
        ref_species=ref_species,
        query_species=query_species,
        kmer_length=kmer_length,
        offset_step=offset_step,
        bwa_params=bwa_params,
    )
    track_name = track.track_name
    target_dir = db_root / ref_species
    track_dir = track.track_dir
    target_dir.mkdir(parents=True, exist_ok=True)

    if track_dir.exists() and not force:
        raise FileExistsError(
            f"Track '{track_name}' already exists for reference '{ref_species}'. "
            "Use --force to overwrite imported files and metadata."
        )
    track_dir.mkdir(parents=True, exist_ok=True)

    all_bw_dst = track_dir / "map_all.bw"
    uniq_bw_dst = track_dir / "map_uniq.bw"
    shutil.copy2(all_bw_src, all_bw_dst)
    shutil.copy2(uniq_bw_src, uniq_bw_dst)

    input_path = Path(input_fasta) if input_fasta else None
    reference_path = Path(reference_fasta) if reference_fasta else None

    computed_ref_md5 = reference_fasta_md5
    if reference_path and reference_path.exists() and reference_path.is_file():
        computed_ref_md5 = file_md5(reference_path)
        if reference_fasta_md5 and reference_fasta_md5 != computed_ref_md5:
            raise ValueError(
                "reference_fasta_md5 does not match the provided reference FASTA file. "
                f"provided={reference_fasta_md5}, computed={computed_ref_md5}"
            )

    if reference_path and computed_ref_md5:
        target_meta_yaml = target_dir / f"{ref_species}.yaml"
        target_md5_file = target_dir / f"{ref_species}.md5"
        target_meta_content = {
            "reference_name": ref_species,
            "reference_fasta": str(reference_path.resolve()),
            "reference_fasta_md5": computed_ref_md5,
            "last_updated": datetime.now().isoformat(timespec="seconds"),
        }
        with target_meta_yaml.open("w", encoding="utf-8") as handle:
            yaml.safe_dump(target_meta_content, handle, sort_keys=False)
        with target_md5_file.open("w", encoding="utf-8") as handle:
            handle.write(f"{computed_ref_md5}  {reference_path.name}\n")

    track.info = {
        "generation_date": datetime.now().isoformat(timespec="seconds"),
        "wizardeye_version": PACKAGE_VERSION,
        "reference": str(reference_path) if reference_path else ref_species,
        "reference_fasta_md5": computed_ref_md5,
        "input": str(input_path) if input_path else query_species,
        "tags": from_charlist_to_list(tags, lowercase=True),
        "kmer_size": kmer_length,
        "sliding_window": offset_step,
        "mapping_tool": mapping_tool,
        "bwa_parameters": {
            "-n": bwa_params.missing_prob_err_rate,
            "-o": bwa_params.max_gap_opens,
            "-l": bwa_params.seed_length,
            "-t": bwa_params.threads,
            "-R": bwa_params.r_best_hits,
            "-sn": bwa_params.samse_n,
        },
        "imported_track": True,
        "import_sources": {
            "map_all_bw": str(all_bw_src.resolve()),
            "map_uniq_bw": str(uniq_bw_src.resolve()),
        },
    }

    param_yaml = track.save_info()

    return {
        "track_dir": track_dir,
        "param_yaml": param_yaml,
        "map_all_bw": all_bw_dst,
        "map_uniq_bw": uniq_bw_dst,
    }


def get_tracks(
    ref_species: str,
    kmer_length: int,
    offset_step: int,
    db_root: str,
    bwa_params: Optional[BWAParameters] = None,
) -> List[Track]:
    """Collect tracks for one reference filtered by generation parameters."""
    ref_dir = Path(db_root) / ref_species
    if not ref_dir.exists() or not ref_dir.is_dir():
        raise ValueError(f"Reference '{ref_species}' not found in {db_root}")

    tracks: List[Track] = []
    for track_dir in sorted(
        [p for p in ref_dir.iterdir() if p.is_dir()], key=lambda p: p.name
    ):
        track = Track.from_param_yaml(
            db_root=db_root, ref_species=ref_species, track_dir=track_dir
        )
        if track is None:
            continue

        if (
            kmer_length is not None
            and track.identity.parameters.kmer_length != kmer_length
        ):
            continue

        if (
            offset_step is not None
            and track.identity.parameters.offset_step != offset_step
        ):
            continue

        if bwa_params is not None:
            if (
                track.identity.parameters.bwa_params.missing_prob_err_rate
                != bwa_params.missing_prob_err_rate
            ):
                continue
            if (
                track.identity.parameters.bwa_params.max_gap_opens
                != bwa_params.max_gap_opens
            ):
                continue
            if (
                track.identity.parameters.bwa_params.seed_length
                != bwa_params.seed_length
            ):
                continue
            if track.identity.parameters.bwa_params.all_aln != bwa_params.all_aln:
                continue
            if track.identity.parameters.bwa_params.threads != bwa_params.threads:
                continue
            if (
                track.identity.parameters.bwa_params.r_best_hits
                != bwa_params.r_best_hits
            ):
                continue
            if track.identity.parameters.bwa_params.samse_n != bwa_params.samse_n:
                continue

        tracks.append(track)

    return tracks


def get_corresponding_tracks(
    ref_species: str,
    db_root: Union[str, Path],
    query_species: Optional[str],
    tags: Optional[List[str]],
    kmer_length: Optional[int],
    offset_step: Optional[int],
    bwa_params: Optional[BWAParameters] = None,
) -> Optional[List[Track]]:
    if query_species is None and not tags:
        log(
            "At least one track name or tag must be provided to find corresponding tracks.",
            "E",
        )
        return None

    tracks = get_tracks(
        ref_species=ref_species,
        kmer_length=kmer_length,
        offset_step=offset_step,
        db_root=str(db_root),
        bwa_params=bwa_params,
    )

    if query_species:
        tracks = [t for t in tracks if t.identity.query_species == query_species]

    if tags:
        tracks = [
            t
            for t in tracks
            if any(
                tag in from_charlist_to_list(t.info.get("tags", []), lowercase=True)
                for tag in tags
            )
        ]

    return tracks


def from_tags_get_tracks(
    ref_species: str,
    tags: List[str],
    kmer_length: int,
    offset_step: int,
    db_root: Union[str, Path],
    bwa_params: Optional[BWAParameters] = None,
) -> List[str]:
    """Return track names matching all requested generation parameters and at least one tag."""
    normalized_tags = set(from_charlist_to_list(tags, lowercase=True))
    tracks = get_tracks(
        ref_species=ref_species,
        kmer_length=kmer_length,
        offset_step=offset_step,
        db_root=str(db_root),
        bwa_params=bwa_params,
    )

    matching_track_names: List[str] = []
    for track in tracks:
        track_tags = set(
            from_charlist_to_list(track.info.get("tags", []), lowercase=True)
        )
        if track_tags & normalized_tags:
            matching_track_names.append(track.track_name)

    return sorted(set(matching_track_names))


def get_refs(db_root: str) -> List[str]:
    """List reference species available in the database."""
    db_root = Path(db_root)
    if not db_root.exists() or not db_root.is_dir():
        raise ValueError(f"Database root '{db_root}' does not exist")

    refs = [
        d.name for d in db_root.iterdir() if d.is_dir() and not d.name.startswith(".")
    ]
    return sorted(refs)


def resolve_requested_track_names(
    requested_tracks: List[str],
    available_tracks: List,
    ref: str,
    context: str,
) -> List[str]:
    """Resolve user-provided identifiers to canonical track names.

    Accepted identifiers for one track are:
    - canonical track name (query_k..._w..._n..._o..._l...)
    - query species identifier
    - query/input name from param.yaml
    """
    resolved: List[str] = []
    for requested in requested_tracks:
        matches = [
            track
            for track in available_tracks
            if requested
            in {track.track_name, track.identity.query_species, track.query_name}
        ]
        if not matches:
            log(
                f"Track '{requested}' does not exist for reference '{ref}' with specified parameters and will be ignored in {context}.",
                "W",
            )
            continue

        if len(matches) > 1 and not any(
            track.track_name == requested for track in matches
        ):
            matching_names = ", ".join(sorted(track.track_name for track in matches))
            log(
                f"Track identifier '{requested}' is ambiguous for reference '{ref}' in {context}. "
                f"Use one canonical track name: {matching_names}",
                "E",
            )
            raise ValueError(
                f"Ambiguous track identifier '{requested}' for reference '{ref}' in {context}"
            )

        if any(track.track_name == requested for track in matches):
            matches = [track for track in matches if track.track_name == requested]

        for match in matches:
            if match.track_name not in resolved:
                resolved.append(match.track_name)

    return resolved


# Display functions
def print_available_species(ref_species: str, db_root: Union[str, Path]) -> Path:
    """Print one-row-per-track table with merged tags, k and sliding window values."""
    db_root = Path(db_root)
    ref_dir = db_root / ref_species
    if not ref_dir.exists() or not ref_dir.is_dir():
        log(f"Reference '{ref_species}' not found in {db_root}", "E")
        return

    track_dirs = [
        d for d in ref_dir.iterdir() if d.is_dir() and not d.name.startswith(".")
    ]
    if not track_dirs:
        log(f"No tracks found for reference '{ref_species}' in {ref_dir}", "E")
        return

    grouped_tracks: Dict[str, Dict[str, set]] = {}
    for track_dir in sorted(track_dirs, key=lambda p: p.name):
        track_name = track_dir.name
        logical_track_name = track_name
        param_yaml = track_dir / "param.yaml"
        tags_set = set()
        k_set = set()
        sliding_set = set()

        if param_yaml.exists():
            try:
                with param_yaml.open("r", encoding="utf-8") as handle:
                    content = yaml.safe_load(handle) or {}
                if isinstance(content, dict):
                    query_name = get_name_from_param(content)
                    if query_name:
                        logical_track_name = str(query_name)

                    tags = content.get("tags", [])
                    if isinstance(tags, list):
                        tags_set.update(
                            from_charlist_to_list(
                                [str(tag) for tag in tags], lowercase=True
                            )
                        )
                    k = content.get("kmer_size")
                    s = content.get("sliding_window")
                    if k is not None:
                        k_set.add(str(k))
                    if s is not None:
                        sliding_set.add(str(s))
            except Exception as e:
                log(f"Failed to read {param_yaml}: {e}", "W")

        if logical_track_name not in grouped_tracks:
            grouped_tracks[logical_track_name] = {
                "tags": set(),
                "k": set(),
                "w": set(),
            }

        grouped_tracks[logical_track_name]["tags"].update(tags_set)
        grouped_tracks[logical_track_name]["k"].update(k_set)
        grouped_tracks[logical_track_name]["w"].update(sliding_set)

    rows = []
    for logical_track_name in sorted(grouped_tracks):
        tags_values = sorted(grouped_tracks[logical_track_name]["tags"])
        k_values = sorted(
            grouped_tracks[logical_track_name]["k"],
            key=lambda x: int(x) if x.isdigit() else x,
        )
        w_values = sorted(
            grouped_tracks[logical_track_name]["w"],
            key=lambda x: int(x) if x.isdigit() else x,
        )

        tags_str = ", ".join(tags_values) if tags_values else "-"
        k_value = ", ".join(k_values) if k_values else "-"
        sliding_value = ", ".join(w_values) if w_values else "-"
        rows.append((logical_track_name, tags_str, k_value, sliding_value))

    headers = ("Track name", "Tags", "k", "Offset")
    widths = [len(h) for h in headers]
    for row in rows:
        for idx, value in enumerate(row):
            widths[idx] = max(widths[idx], len(value))

    def _format_row(values: tuple) -> str:
        return " | ".join(str(v).ljust(widths[i]) for i, v in enumerate(values))

    log(f"> {ref_species}", "I")
    print(_format_row(headers))
    print("-+-".join("-" * w for w in widths))
    for row in rows:
        print(_format_row(row))


def print_full_catalogue(db_root: Union[str, Path]) -> None:
    """Print track tables for every reference target found in the database root."""
    db_root = Path(db_root)
    if not db_root.exists() or not db_root.is_dir():
        log(f"Database root '{db_root}' does not exist", "E")
        return

    target_dirs = [
        d for d in db_root.iterdir() if d.is_dir() and not d.name.startswith(".")
    ]
    if not target_dirs:
        log(f"No reference targets found in {db_root}", "E")
        return

    for target_dir in sorted(target_dirs, key=lambda p: p.name):
        print_available_species(target_dir.name, db_root)


# -- Migration functions ---

def _parse_old_track_name(track_dir_name: str) -> Optional[Dict[str, Optional[str]]]:
    """Parse old track directory name format (e.g., sus_scrofa_k35_w1_n1_o2_l1).
    
    Returns dict with keys: query_species, kmer_length, offset_step, 
    bwa_missing_prob_err_rate, bwa_max_gap_opens, bwa_seed_length.
    """
    import re
    
    # Expected format: {query_species}_k{k}_w{w}[_n{n}[_o{o}[_l{l}]]]
    # The query_species can contain underscores, so we need to find the first _k and _w
    # Pattern to match: everything before _k{number}_w{number}
    pattern = r'^(.+)_k(\d+)_w(\d+)(?:_n([\d.]+))?(?:_o(\d+))?(?:_l(\d+))?$'
    match = re.match(pattern, track_dir_name)
    
    if not match:
        return None
    
    query_species = match.group(1)
    kmer_length = match.group(2)
    offset_step = match.group(3)
    bwa_n = match.group(4)
    bwa_o = match.group(5)
    bwa_l = match.group(6)
    
    result: Dict[str, Optional[str]] = {
        "query_species": query_species,
        "kmer_length": kmer_length,
        "offset_step": offset_step,
        "bwa_missing_prob_err_rate": bwa_n,
        "bwa_max_gap_opens": bwa_o,
        "bwa_seed_length": bwa_l,
    }
    
    return result


def migrate_database(
    db_root: Union[str, Path],
    bwa_r_best_hits: int = 2147483647,
    bwa_samse_n: int = 2147483647,
) -> Dict[str, object]:
    """Migrate database tracks from old naming format to new format with bwa hash.
    
    Old format: {query}_k{k}_w{w}_n{n}_o{o}_l{l}/
    New format: {query}_k{k}_w{w}_bwa{hash}/
    
    Creates a copy of the database with _v0_1_3 suffix and applies migration on the copy.
    
    Args:
        db_root: Path to the database root directory.
        bwa_r_best_hits: -bR value for BWA parameters (default: 2147483647).
        bwa_samse_n: -bsn value for BWA parameters (default: 2147483647).
    
    Returns:
        Dict with migration statistics: migrated_count, new_db_path, warnings, errors.
    
    Raises:
        ValueError: If database is invalid or no tracks need migration.
    """
    db_root = Path(db_root)
    
    if not valid_database(db_root):
        raise ValueError(f"Invalid database at {db_root}")
    
    # Create a copy of the database with _v0_1_3 suffix
    new_db_root = db_root.parent / f"{db_root.name}_v0_1_3"
    
    if new_db_root.exists():
        raise ValueError(f"Migrated database already exists at {new_db_root}")
    
    log(f"Creating a copy of the database at {new_db_root}", "I")
    
    # Copy the entire database directory
    import shutil
    shutil.copytree(db_root, new_db_root)
    
    # Get all reference directories in the new copy
    ref_dirs = [
        d for d in new_db_root.iterdir() if d.is_dir() and not d.name.startswith(".")
    ]
    
    if not ref_dirs:
        raise ValueError(f"No reference directories found in {new_db_root}")
    
    migrated_count = 0
    warnings_list: List[str] = []
    errors_list: List[str] = []
    
    # Check if defaults are used
    if bwa_r_best_hits == 2147483647 and bwa_samse_n == 2147483647:
        warnings_list.append(
            "WARNING: Using default -bR and -bsn values (2147483647). "
            "Version 0.1.2 had a major bug that broke tracks and did not reference any alternative alignments. "
            "For more information, see: https://github.com/TheLokj/WizardEye/issues/3"
        )
    
    for ref_dir in sorted(ref_dirs, key=lambda p: p.name):
        ref_name = ref_dir.name
        track_dirs = [
            d for d in ref_dir.iterdir() if d.is_dir() and not d.name.startswith(".")
        ]
        
        for old_track_dir in sorted(track_dirs, key=lambda p: p.name):
            old_name = old_track_dir.name
            
            # Check if already in new format (contains _bwa)
            if "_bwa" in old_name:
                continue
            
            # Parse old format
            parsed = _parse_old_track_name(old_name)
            if parsed is None:
                warnings_list.append(f"Skipping {ref_name}/{old_name}: cannot parse old format")
                continue
            
            try:
                kmer_length = int(parsed["kmer_length"])
                offset_step = int(parsed["offset_step"])
            except (ValueError, TypeError):
                errors_list.append(f"Skipping {ref_name}/{old_name}: invalid kmer_length or offset_step")
                continue
            
            # Extract existing BWA parameters from param.yaml if available
            old_param_yaml = old_track_dir / "param.yaml"
            existing_bwa_n = 0.01
            existing_bwa_o = 2
            existing_bwa_l = 16500
            old_param_content = None
            
            if old_param_yaml.exists():
                try:
                    with old_param_yaml.open("r", encoding="utf-8") as handle:
                        old_param_content = yaml.safe_load(handle) or {}
                    bwa_params = old_param_content.get("bwa_parameters", {})
                    if isinstance(bwa_params, dict):
                        existing_bwa_n = float(bwa_params.get("-n", 0.01))
                        existing_bwa_o = int(bwa_params.get("-o", 2))
                        existing_bwa_l = int(bwa_params.get("-l", 16500))
                except Exception:
                    pass
            
            # Create new BWAParameters with extracted values + provided bR and bns
            new_bwa_params = BWAParameters(
                missing_prob_err_rate=existing_bwa_n,
                max_gap_opens=existing_bwa_o,
                seed_length=existing_bwa_l,
                all_aln=False,
                threads=1,
                r_best_hits=bwa_r_best_hits,
                samse_n=bwa_samse_n,
            )
            
            # Generate new hash and track name
            bwa_hash = get_bwa_params_hash(new_bwa_params)
            query_species = parsed["query_species"]
            new_track_name = f"{query_species}_k{kmer_length}_w{offset_step}_bwa{bwa_hash}"
            
            # Create new track directory path
            new_track_dir = ref_dir / new_track_name
            
            if new_track_dir.exists():
                errors_list.append(f"Skipping {ref_name}/{old_name}: new path {new_track_name} already exists")
                continue
            
            # Rename directory in the copied database
            old_track_dir.rename(new_track_dir)
            
            # Update metadata in param.yaml
            if old_param_content is not None:
                new_param_yaml = new_track_dir / "param.yaml"
                try:
                    param_content = old_param_content.copy() if isinstance(old_param_content, dict) else {}
                    
                    # Update bwa_parameters with new values
                    if "bwa_parameters" not in param_content:
                        param_content["bwa_parameters"] = {}
                    
                    param_content["bwa_parameters"]["-n"] = new_bwa_params.missing_prob_err_rate
                    param_content["bwa_parameters"]["-o"] = new_bwa_params.max_gap_opens
                    param_content["bwa_parameters"]["-l"] = new_bwa_params.seed_length
                    param_content["bwa_parameters"]["-t"] = new_bwa_params.threads
                    param_content["bwa_parameters"]["-R"] = new_bwa_params.r_best_hits
                    param_content["bwa_parameters"]["-sn"] = new_bwa_params.samse_n
                    
                    # Update kmer_size and sliding_window if needed
                    param_content["kmer_size"] = kmer_length
                    param_content["sliding_window"] = offset_step
                    
                    with new_param_yaml.open("w", encoding="utf-8") as handle:
                        yaml.safe_dump(param_content, handle, sort_keys=False)
                except Exception as e:
                    errors_list.append(f"Failed to update metadata for {ref_name}/{new_track_name}: {e}")
            
            migrated_count += 1
            log(f"Migrated {ref_name}/{old_name} -> {new_track_name}", "I")
    
    result = {
        "migrated_count": migrated_count,
        "new_db_path": new_db_root,
        "warnings": warnings_list,
        "errors": errors_list,
    }
    
    return result
