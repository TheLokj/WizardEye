# -*- coding: utf-8 -*-

import sys
import argparse
import shutil
import yaml

from dataclasses import dataclass
from pathlib import Path
from datetime import datetime
from importlib import metadata
from typing import Dict, List, Optional, Tuple, Union

from .utils import log, query_name_from_param
from .utils import file_md5
from .utils import from_charlist_to_list
from .version import PACKAGE_VERSION

# -- Database management functions --
def init_db(base_dir: Union[str, Path] = ".") -> Path:
	"""Create /database and its info.yaml file."""
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

# -- Track classes and functions --
@dataclass(frozen=True)
class TrackParameters:
	"""Alignment and generation parameters that define one track version."""
	kmer_length: int
	offset_step: int
	bwa_missing_prob_err_rate: float = 0.01
	bwa_max_gap_opens: int = 2
	bwa_seed_length: int = 16500

	def __post_init__(self) -> None:
		if self.kmer_length < 1:
			raise ValueError("kmer_length must be a positive integer")
		if self.offset_step < 1:
			raise ValueError("offset_step must be a positive integer")

	@property
	def n_value(self) -> float:
		return float(self.bwa_missing_prob_err_rate)

@dataclass(frozen=True)
class TrackIdentity:
	"""Track identity independent from file system state."""
	ref_species: str
	query_species: str
	parameters: TrackParameters

	@property
	def name(self) -> str:
		return (
			f"{self.query_species}_k{self.parameters.kmer_length}_w{self.parameters.offset_step}"
			f"_n{self.parameters.n_value:g}_o{self.parameters.bwa_max_gap_opens}"
			f"_l{self.parameters.bwa_seed_length}"
		)

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
	def mappability_all_bw(self) -> Path:
		return self.track_dir / "mappability_all.bw"

	@property
	def mappability_uniq_bw(self) -> Path:
		return self.track_dir / "mappability_uniq.bw"

	@property
	def query_name(self) -> str:
		query_name = query_name_from_param(self.info)
		return str(query_name)

	def exists(self) -> bool:
		return self.mappability_all_bw.exists() and self.mappability_uniq_bw.exists()

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
		bwa_missing_prob_err_rate: float = 0.01,
		bwa_max_gap_opens: int = 2,
		bwa_seed_length: int = 16500,
		info: Optional[Dict] = None,
	) -> "Track":
		params = TrackParameters(
			kmer_length=kmer_length,
			offset_step=offset_step,
			bwa_missing_prob_err_rate=bwa_missing_prob_err_rate,
			bwa_max_gap_opens=bwa_max_gap_opens,
			bwa_seed_length=bwa_seed_length,
		)
		identity = TrackIdentity(ref_species=ref_species, query_species=query_species, parameters=params)
		return cls(db_root=Path(db_root), identity=identity, info=info or {})

	@classmethod
	def from_param_yaml(cls, db_root: Union[str, Path], ref_species: str, track_dir: Path) -> Optional["Track"]:
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
		except (TypeError, ValueError):
			bwa_n, bwa_o, bwa_l = 0.01, 2, 16500

		query_species = track_dir.name.split("_k", 1)[0]
		return cls.from_param(
			db_root=db_root,
			ref_species=ref_species,
			query_species=query_species,
			kmer_length=kmer_size,
			offset_step=sliding_window,
			bwa_missing_prob_err_rate=bwa_n,
			bwa_max_gap_opens=bwa_o,
			bwa_seed_length=bwa_l,
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

	def update_tags(self, tags: Optional[List[str]]) -> Tuple[Path, List[str], List[str]]:
		content = self.load_info()

		new_tags = from_charlist_to_list(tags, lowercase=True)
		existing_tags = content.get("tags", [])
		if not isinstance(existing_tags, list):
			existing_tags = []
		old_tags = from_charlist_to_list([str(tag) for tag in existing_tags], lowercase=True)

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
	bwa_missing_prob_err_rate: float = 0.01,
	bwa_max_gap_opens: int = 2,
	bwa_seed_length: int = 16500,
) -> bool:
	"""Check if a track exists for the given reference/query and k/s parameters."""
	track = Track.from_param(
		db_root=db_root,
		ref_species=ref_species,
		query_species=query_species,
		kmer_length=kmer_length,
		offset_step=offset_step,
		bwa_missing_prob_err_rate=bwa_missing_prob_err_rate,
		bwa_max_gap_opens=bwa_max_gap_opens,
		bwa_seed_length=bwa_seed_length,
	)
	return track.exists()

def update_track_tags(
	ref_species: str,
	query_species: str,
	kmer_length: int,
	offset_step: int,
	tags: Optional[List[str]],
	db_root: Union[str, Path],
	bwa_missing_prob_err_rate: float = 0.01,
	bwa_max_gap_opens: int = 2,
	bwa_seed_length: int = 16500,
) -> Tuple[Path, List[str], List[str]]:
	"""Replace all tags in the track param.yaml and return old/new values."""
	track = Track.from_param(
		db_root=db_root,
		ref_species=ref_species,
		query_species=query_species,
		kmer_length=kmer_length,
		offset_step=offset_step,
		bwa_missing_prob_err_rate=bwa_missing_prob_err_rate,
		bwa_max_gap_opens=bwa_max_gap_opens,
		bwa_seed_length=bwa_seed_length,
	)
	return track.update_tags(tags)

def import_track(
	ref_species: str,
	query_species: str,
	kmer_length: int,
	offset_step: int,
	mappability_all_bw: Union[str, Path],
	mappability_uniq_bw: Union[str, Path],
	db_root: Union[str, Path],
	input_fasta: Optional[Union[str, Path]] = None,
	reference_fasta: Optional[Union[str, Path]] = None,
	reference_fasta_md5: Optional[str] = None,
	tags: Optional[List[str]] = None,
	mapping_tool: str = "bwa aln",
	bwa_missing_prob_err_rate: float = 0.01,
	bwa_max_gap_opens: int = 2,
	bwa_seed_length: int = 16500,
	n_threads: int = 1,
	force: bool = False,
) -> Dict[str, Path]:
	"""Import an externally generated track by copying BigWig files and writing param.yaml."""
	db_root = Path(db_root)

	all_bw_src = Path(mappability_all_bw)
	uniq_bw_src = Path(mappability_uniq_bw)
	if not all_bw_src.exists() or not all_bw_src.is_file():
		raise FileNotFoundError(f"mappability_all.bw source file not found: {all_bw_src}")
	if not uniq_bw_src.exists() or not uniq_bw_src.is_file():
		raise FileNotFoundError(f"mappability_uniq.bw source file not found: {uniq_bw_src}")

	track = Track.from_param(
		db_root=db_root,
		ref_species=ref_species,
		query_species=query_species,
		kmer_length=kmer_length,
		offset_step=offset_step,
		bwa_missing_prob_err_rate=bwa_missing_prob_err_rate,
		bwa_max_gap_opens=bwa_max_gap_opens,
		bwa_seed_length=bwa_seed_length,
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

	all_bw_dst = track_dir / "mappability_all.bw"
	uniq_bw_dst = track_dir / "mappability_uniq.bw"
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
		"wizardeye_version": "dev",
		"reference": str(reference_path) if reference_path else ref_species,
		"reference_fasta_md5": computed_ref_md5,
		"input": str(input_path) if input_path else query_species,
		"tags": from_charlist_to_list(tags, lowercase=True),
		"kmer_size": kmer_length,
		"sliding_window": offset_step,
		"mapping_tool": mapping_tool,
		"bwa_parameters": {
			"-n": bwa_missing_prob_err_rate,
			"-o": bwa_max_gap_opens,
			"-l": bwa_seed_length,
			"-t": n_threads,
		},
		"imported_track": True,
		"import_sources": {
			"mappability_all_bw": str(all_bw_src.resolve()),
			"mappability_uniq_bw": str(uniq_bw_src.resolve()),
		},
	}

	param_yaml = track.save_info()

	return {
		"track_dir": track_dir,
		"param_yaml": param_yaml,
		"mappability_all_bw": all_bw_dst,
		"mappability_uniq_bw": uniq_bw_dst,
	}

def get_tracks(
 	ref_species: str,
	kmer_length: int,
	offset_step: int,
	db_root: str,
	bwa_missing_prob_err_rate: Optional[float] = None,
	bwa_max_gap_opens: Optional[int] = None,
	bwa_seed_length: Optional[int] = None,
) -> List[Track]:
	"""Collect tracks for one reference filtered by generation parameters."""
	ref_dir = Path(db_root) / ref_species
	if not ref_dir.exists() or not ref_dir.is_dir():
		raise ValueError(f"Reference '{ref_species}' not found in {db_root}")

	tracks: List[Track] = []
	for track_dir in sorted([p for p in ref_dir.iterdir() if p.is_dir()], key=lambda p: p.name):
		track = Track.from_param_yaml(db_root=db_root, ref_species=ref_species, track_dir=track_dir)
		if track is None:
			continue

		if kmer_length is not None and track.identity.parameters.kmer_length != kmer_length:
			continue

		if offset_step is not None and track.identity.parameters.offset_step != offset_step:
			continue

		if (
			bwa_missing_prob_err_rate is not None
			and track.identity.parameters.bwa_missing_prob_err_rate != bwa_missing_prob_err_rate
		):
			continue

		if (
			bwa_max_gap_opens is not None
			and track.identity.parameters.bwa_max_gap_opens != bwa_max_gap_opens
		):
			continue

		if (
			bwa_seed_length is not None
			and track.identity.parameters.bwa_seed_length != bwa_seed_length
		):
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
	bwa_missing_prob_err_rate: Optional[float] = None,
	bwa_max_gap_opens: Optional[int] = None,
	bwa_seed_length: Optional[int] = None,
) -> Optional[List[Track]]:
	if query_species is None and not tags:
		log("At least one track name or tag must be provided to find corresponding tracks.", "E")
		return None
	
	tracks = get_tracks(
		ref_species=ref_species,
		kmer_length=kmer_length,
		offset_step=offset_step,
		db_root=str(db_root),
		bwa_missing_prob_err_rate=bwa_missing_prob_err_rate,
		bwa_max_gap_opens=bwa_max_gap_opens,
		bwa_seed_length=bwa_seed_length,
	)
	
	if query_species:
		tracks = [t for t in tracks if t.identity.query_species == query_species]

	if tags:
		tracks = [
			t for t in tracks
			if any(tag in from_charlist_to_list(t.info.get("tags", []), lowercase=True) for tag in tags)
		]
	
	return tracks

def from_tags_get_tracks(
	ref_species: str,
	tags: List[str],
	kmer_length: int,
	offset_step: int,
	db_root: Union[str, Path],
	bwa_missing_prob_err_rate: float,
	bwa_max_gap_opens: int,
	bwa_seed_length: int,
) -> List[str]:
	"""Return track names matching all requested generation parameters and at least one tag."""
	normalized_tags = set(from_charlist_to_list(tags, lowercase=True))
	tracks = get_tracks(
		ref_species=ref_species,
		kmer_length=kmer_length,
		offset_step=offset_step,
		db_root=str(db_root),
		bwa_missing_prob_err_rate=bwa_missing_prob_err_rate,
		bwa_max_gap_opens=bwa_max_gap_opens,
		bwa_seed_length=bwa_seed_length,
	)

	matching_track_names: List[str] = []
	for track in tracks:
		track_tags = set(from_charlist_to_list(track.info.get("tags", []), lowercase=True))
		if track_tags & normalized_tags:
			matching_track_names.append(track.track_name)

	return sorted(set(matching_track_names))

def get_refs(db_root: str) -> List[str]:
	"""List reference species available in the database."""
	db_root = Path(db_root)
	if not db_root.exists() or not db_root.is_dir():
		raise ValueError(f"Database root '{db_root}' does not exist")

	refs = [d.name for d in db_root.iterdir() if d.is_dir() and not d.name.startswith('.')]
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
			track for track in available_tracks
			if requested in {track.track_name, track.identity.query_species, track.query_name}
		]
		if not matches:
			log(
				f"Track '{requested}' does not exist for reference '{ref}' with specified parameters and will be ignored in {context}.",
				"W",
			)
			continue

		if len(matches) > 1 and not any(track.track_name == requested for track in matches):
			matching_names = ", ".join(sorted(track.track_name for track in matches))
			log(
				f"Track identifier '{requested}' is ambiguous for reference '{ref}' in {context}. "
				f"Use one canonical track name: {matching_names}",
				"E",
			)
			raise ValueError(f"Ambiguous track identifier '{requested}' for reference '{ref}' in {context}")

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

	track_dirs = [d for d in ref_dir.iterdir() if d.is_dir() and not d.name.startswith('.')]
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
					query_name = query_name_from_param(content)
					if query_name:
						logical_track_name = str(query_name)

					tags = content.get("tags", [])
					if isinstance(tags, list):
						tags_set.update(from_charlist_to_list([str(tag) for tag in tags], lowercase=True))
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
		k_values = sorted(grouped_tracks[logical_track_name]["k"], key=lambda x: int(x) if x.isdigit() else x)
		w_values = sorted(grouped_tracks[logical_track_name]["w"], key=lambda x: int(x) if x.isdigit() else x)

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
	log(_format_row(headers), "I")
	log("-+-".join("-" * w for w in widths), "I")
	for row in rows:
		log(_format_row(row), "I")


def print_full_catalogue(db_root: Union[str, Path]) -> None:
	"""Print track tables for every reference target found in the database root."""
	db_root = Path(db_root)
	if not db_root.exists() or not db_root.is_dir():
		log(f"Database root '{db_root}' does not exist", "E")
		return

	target_dirs = [d for d in db_root.iterdir() if d.is_dir() and not d.name.startswith('.')]
	if not target_dirs:
		log(f"No reference targets found in {db_root}", "E")
		return

	for target_dir in sorted(target_dirs, key=lambda p: p.name):
		print_available_species(target_dir.name, db_root)
