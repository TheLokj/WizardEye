# -*- coding: utf-8 -*-

import subprocess
import typer

from importlib import metadata
from pathlib import Path
from typing import List, Optional

from .mappability import create_mappability_track
from .db import (
	import_track,
	init_db,
	print_full_catalogue,
	update_track_tags,
	valid_database,
	check_track_exists,
	from_tags_get_tracks,
	get_refs,
	get_tracks,
)
from .filter import generate_mask, filter_bam
from .utils import from_charlist_to_list, log
from .version import DISPLAY_VERSION, PACKAGE_VERSION

app = typer.Typer(help="Cross-mappability helper CLI.")


def _resolve_requested_track_names(
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
			raise typer.Exit(code=1)

		if any(track.track_name == requested for track in matches):
			matches = [track for track in matches if track.track_name == requested]

		for match in matches:
			if match.track_name not in resolved:
				resolved.append(match.track_name)

	return resolved


def _get_runtime_version() -> str:
	"""Return installed package version, or fallback display version in source mode."""
	try:
		installed_version = metadata.version("wizardeye")
		if installed_version == PACKAGE_VERSION:
			return DISPLAY_VERSION
		return installed_version
	except metadata.PackageNotFoundError:
		return DISPLAY_VERSION


def _version_callback(value: bool) -> None:
	if value:
		typer.echo(f"WizardEye version: {_get_runtime_version()}")
		raise typer.Exit(code=0)


@app.callback()
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

# -- CLI commands --
@app.command()
def database(
	init: bool = typer.Option(False, help="Path to the database root directory."),
	catalogue: bool = typer.Option(False, "-c", "--catalogue", help="Print the full database catalogue after initialization."),
	db_root: str = typer.Option(..., "-d", "--db-root", help="Path to the database root directory."),
	):
	
	if init:
		try:
			init_db(db_root)
		except FileExistsError as e:
			log(str(e), "E")
			raise typer.Exit(code=1)
	elif catalogue:
		print_full_catalogue(db_root)
		raise typer.Exit(code=0)
	elif not valid_database(db_root):
		log(f"To initialize the database, run 'database --init {db_root}'", "E")
		raise typer.Exit(code=1)
	else:
		log(f"Database at '{db_root}' is valid and ready to use.", "S")
		raise typer.Exit(code=0)

# Main alignment command with all parameters for track generation
# Same behavior that generate_cross_mappability_filter_bwa.sh 
@app.command()
def align(
	# Input parameters
	input_fasta: Optional[List[str]] = typer.Option(
		None,
		"-i",
		help="Path(s) to FASTA file(s) to align on the target (query). Repeat -i and/or use comma-separated values.",
	),
	input_target: str = typer.Option(None, "-r", help="Path to the reference target FASTA."),
	track_id: Optional[str] = typer.Option(
		None,
		"--track_ID",
		"--track-id",
		help="Manual identifier used to name/reference the generated track (defaults to input FASTA stem).",
	),
	tag: Optional[List[str]] = typer.Option(None, "--tag", "-t", help="Tag(s) to associate to the track. Comma-separated (e.g. tag1,tag2)."),

	# Splitting parameters
	kmer_length: int = typer.Option(None, "-k", help="Length of k-mers to produce."),
	offset_step: int = typer.Option(1, "-w", "--offset-step", "--sliding-window", help="Offset/sliding window step for k-mers."),
	
	#  Alignment parameters
	bwa_missing_prob_err_rate: float = typer.Option(0.01, "-bn", help="BWA aln -n."),
	bwa_max_gap_opens: int = typer.Option(2, "-bo", help="BWA aln -o."),
	bwa_seed_length: int = typer.Option(16500, "-bl", help="BWA aln -l."),

	# Parallelisation parameters
	n_threads: int = typer.Option(1, "-j", help="Number of threads for parallelisation."),
	chunk_size: int = typer.Option(
		2000000,
		"-cs",
		"--chunk-size",
		help="Number of k-mers per chunk for parallel alignment.",
	),

	# Other parameters
	db_root: str = typer.Option(..., "-d", "--db-root", help="Path to the database root directory."),
	force: bool = typer.Option(None, help="Force track creation even if it already exists."),
):

	if input_fasta and input_target and kmer_length:
		if not valid_database(db_root):
			log(f"To initialize the database, run 'database --init {db_root}'", "E")
			raise typer.Exit(code=1)

		input_fastas = from_charlist_to_list(input_fasta)
		if not input_fastas:
			log("No input FASTA provided after parsing -i values.", "E")
			raise typer.Exit(code=1)

		manual_track_id: Optional[str] = None
		if track_id is not None:
			manual_track_id = track_id.strip()
			if not manual_track_id:
				log("--track_ID cannot be empty.", "E")
				raise typer.Exit(code=1)
			manual_track_id = manual_track_id.replace("/", "_").replace("\\", "_").replace(" ", "_")

		total_inputs = len(input_fastas)
		for idx, one_input_fasta in enumerate(input_fastas, start=1):
			base_query_name = Path(one_input_fasta).stem
			query_track_id = base_query_name if not manual_track_id else f"{base_query_name}__{manual_track_id}"

			track_name = f"{query_track_id}_k{kmer_length}_w{offset_step}_n{float(bwa_missing_prob_err_rate):g}_o{bwa_max_gap_opens}_l{bwa_seed_length}"
			log(f"[{idx}/{total_inputs}] Processing input FASTA: {one_input_fasta}", "I")

			if check_track_exists(
				Path(input_target).stem,
				query_track_id,
				kmer_length,
				offset_step,
				db_root,
				bwa_missing_prob_err_rate,
				bwa_max_gap_opens,
				bwa_seed_length,
			) and not force:
				if tag:
					try:
						updated_yaml, old_tags, new_tags = update_track_tags(
							ref_species=Path(input_target).stem,
							query_species=query_track_id,
							kmer_length=kmer_length,
							offset_step=offset_step,
							tags=tag,
							db_root=db_root,
							bwa_missing_prob_err_rate=bwa_missing_prob_err_rate,
							bwa_max_gap_opens=bwa_max_gap_opens,
							bwa_seed_length=bwa_seed_length,
						)
						old_tags_str = ", ".join(old_tags) if old_tags else "-"
						new_tags_str = ", ".join(new_tags) if new_tags else "-"
						log(f"Track already exists. Tags updated in: {updated_yaml}", "S")
						log(f"Previous tags: {old_tags_str}", "I")
						log(f"New tags: {new_tags_str}", "I")
					except FileNotFoundError as e:
						log(str(e), "E")
						raise typer.Exit(code=1)
					continue

				log(
					f"Track '{track_name}' already exists for reference '{Path(input_target).stem}', skipping track creation.",
					"W",
				)
				continue

			try:
				log(f"Creating track for reference '{Path(input_target).stem}' and query '{query_track_id}'...", "I")
				create_mappability_track(
					input_fasta=one_input_fasta,
					input_target=input_target,
					track_id=query_track_id,
					manual_track_id=manual_track_id,
					kmer_length=kmer_length,
					offset_step=offset_step,
					chunk_size=chunk_size,
					n_threads=n_threads,
					bwa_missing_prob_err_rate=bwa_missing_prob_err_rate,
					bwa_max_gap_opens=bwa_max_gap_opens,
					bwa_seed_length=bwa_seed_length,
					db_root=db_root,
					tags=tag,
				)
			except subprocess.CalledProcessError as e:
				log(f"Subprocess failed: {e}", "E")
				raise typer.Exit(code=2)
			except ValueError as e:
				log(str(e), "E")
				raise typer.Exit(code=2)

		raise typer.Exit(code=0)

	log(f"align requires BAM generation arguments (-i, -r, -k)", "E")
	raise typer.Exit(code=1)

@app.command()
def filter(
	# Alignment parameters used to construct the BAM file
	input_bam = typer.Option(None, "-i", "--input", help="Path to the BAM file to filter."),
	ref: Optional[str] = typer.Option(None, "-r", help="Reference used for alignment."),
	bwa_missing_prob_err_rate: Optional[float] = typer.Option(None, "-bn", help="BWA aln -n used for alignment."),
	bwa_max_gap_opens: Optional[int] = typer.Option(None, "-bo", help="BWA aln -o used for alignment."),
	bwa_seed_length: Optional[int] = typer.Option(None, "-bl", help="BWA aln -l used for alignment."),
	db_root: str = typer.Option(..., "-d", "--db-root", help="Path to the database root directory."),
	
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
		"-s",
		"--stringency",
		"-rc",
		"--cross-stringency",
		help="Stringency threshold in [0.0, 1.0].",
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
	n_threads: int = typer.Option(
		1,
		"-j",
		help="Number of threads for parallel overlap extraction from BigWig tracks.",
	),

	kmer_length: Optional[int] = typer.Option(None, "-k", help="K-mer length to filter on. Must match track generation parameter."),
	offset_step: Optional[int] = typer.Option(None, "-w", "--offset-step", "--sliding-window", help="Offset/sliding window to filter on. Smaller, more sensitive."),
):
	if ref is None:
		log("Reference (-r) must be specified.", "E")
		raise typer.Exit(code=1)

	if ref not in get_refs(db_root):
		log(f"Reference '{ref}' not found in database '{db_root}'.", "E")
		raise typer.Exit(code=1)

	if exclude_tags and exclude_tracks:
		log("Cannot use both --exclude-tags and --exclude-tracks at the same time for filtering. Please choose one.", "E")
		raise typer.Exit(code=1)

	if exclude_tracks or exclude_tags:
		missing_params: List[str] = []
		if kmer_length is None:
			missing_params.append("-k")
		if offset_step is None:
			missing_params.append("-w")
		if bwa_missing_prob_err_rate is None:
			missing_params.append("-bn")
		if bwa_max_gap_opens is None:
			missing_params.append("-bo")
		if bwa_seed_length is None:
			missing_params.append("-bl")

		if missing_params:
			log(
				"For filter safety, all requested track generation parameters are mandatory. "
				f"Missing: {', '.join(missing_params)}",
				"E",
			)
			raise typer.Exit(code=1)

	if exclude_tracks:
		requested_tracks = from_charlist_to_list(exclude_tracks)
		available_tracks = get_tracks(
			ref_species=ref,
			kmer_length=kmer_length,
			offset_step=offset_step,
			db_root=db_root,
			bwa_missing_prob_err_rate=bwa_missing_prob_err_rate,
			bwa_max_gap_opens=bwa_max_gap_opens,
			bwa_seed_length=bwa_seed_length,
		)
		valid_tracks = _resolve_requested_track_names(
			requested_tracks=requested_tracks,
			available_tracks=available_tracks,
			ref=ref,
			context="filtering",
		)

		if not valid_tracks:
			log("No valid tracks remain after validation. Filtering cannot proceed.", "E")
			raise typer.Exit(code=1)

	if exclude_tags:
		requested_tags = from_charlist_to_list(exclude_tags)

		valid_tracks = from_tags_get_tracks(
			ref_species=ref,
			tags=requested_tags,
			kmer_length=kmer_length,
			offset_step=offset_step,
			db_root=db_root,
			bwa_missing_prob_err_rate=bwa_missing_prob_err_rate,
			bwa_max_gap_opens=bwa_max_gap_opens,
			bwa_seed_length=bwa_seed_length,
		)

		if not valid_tracks:
			log("No valid tracks were found from provided tags and parameters. Filtering cannot proceed.", "E")
			raise typer.Exit(code=1)

		log(f"Specified tags: {', '.join(requested_tags)}", "I")
	
	if exclude_tags is None and exclude_tracks is None:
		log("Please provide either --exclude-tags or --exclude-tracks.", "E")
		raise typer.Exit(code=1)

	log(f"Following tracks will be excluded from filtering: {', '.join(valid_tracks)}", "I")

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
			bwa_missing_prob_err_rate=bwa_missing_prob_err_rate,
			bwa_max_gap_opens=bwa_max_gap_opens,
			bwa_seed_length=bwa_seed_length,
			stringency=cross_stringency,
			n_threads=n_threads,
			output_filtered_bam=output_filtered_bam,
			output_excluded_bam=output_excluded_bam,
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

	log(f"Merged mask generated: {filter_result['mask']}", "S")
	log(f"Filtered BAM (kept reads): {filter_result['filtered_bam']}", "S")
	log(f"Excluded BAM (masked reads): {filter_result['excluded_bam']}", "S")
	log(f"Read exclusion report (TSV): {filter_result['report_tsv']}", "S")
	log(
		"Read counts: "
		f"total={filter_result['n_total']}, kept={filter_result['n_filtered']}, excluded={filter_result['n_excluded']}",
		"I",
	)

	raise typer.Exit(code=0)


@app.command()
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
	offset_step: Optional[int] = typer.Option(None, "-w", "--offset-step", "--sliding-window", help="Offset/sliding window to target."),
	bwa_missing_prob_err_rate: Optional[float] = typer.Option(None, "-bn", help="BWA aln -n used for selected tracks."),
	bwa_max_gap_opens: Optional[int] = typer.Option(None, "-bo", help="BWA aln -o used for selected tracks."),
	bwa_seed_length: Optional[int] = typer.Option(None, "-bl", help="BWA aln -l used for selected tracks."),
	cross_stringency: float = typer.Option(
		0.99,
		"-s",
		"--stringency",
		"-rc",
		"--cross-stringency",
		help="Stringency threshold in [0.0, 1.0].",
	),
	db_root: str = typer.Option(..., "-d", "--db-root", help="Path to the database root directory."),
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
	if not valid_database(db_root):
		log(f"To initialize the database, run 'database --init {db_root}'", "E")
		raise typer.Exit(code=1)

	if ref not in get_refs(db_root):
		log(f"Reference '{ref}' not found in database '{db_root}'.", "E")
		raise typer.Exit(code=1)

	if exclude_tags and exclude_tracks:
		log("Cannot use both --exclude-tags and --exclude-tracks at the same time for export. Please choose one.", "E")
		raise typer.Exit(code=1)

	if exclude_tags is None and exclude_tracks is None:
		log("Please provide either --exclude-tags or --exclude-tracks.", "E")
		raise typer.Exit(code=1)

	missing_params: List[str] = []
	if kmer_length is None:
		missing_params.append("-k")
	if offset_step is None:
		missing_params.append("-w")
	if bwa_missing_prob_err_rate is None:
		missing_params.append("-bn")
	if bwa_max_gap_opens is None:
		missing_params.append("-bo")
	if bwa_seed_length is None:
		missing_params.append("-bl")

	if missing_params:
		log(
			"For export safety, all requested track generation parameters are mandatory. "
			f"Missing: {', '.join(missing_params)}",
			"E",
		)
		raise typer.Exit(code=1)

	if exclude_tracks:
		requested_tracks = from_charlist_to_list(exclude_tracks)
		available_tracks = get_tracks(
			ref_species=ref,
			kmer_length=kmer_length,
			offset_step=offset_step,
			db_root=db_root,
			bwa_missing_prob_err_rate=bwa_missing_prob_err_rate,
			bwa_max_gap_opens=bwa_max_gap_opens,
			bwa_seed_length=bwa_seed_length,
		)
		valid_tracks = _resolve_requested_track_names(
			requested_tracks=requested_tracks,
			available_tracks=available_tracks,
			ref=ref,
			context="export",
		)

		if not valid_tracks:
			log("No valid tracks remain after validation. Export cannot proceed.", "E")
			raise typer.Exit(code=1)

	if exclude_tags:
		requested_tags = from_charlist_to_list(exclude_tags)
		valid_tracks = from_tags_get_tracks(
			ref_species=ref,
			tags=requested_tags,
			kmer_length=kmer_length,
			offset_step=offset_step,
			db_root=db_root,
			bwa_missing_prob_err_rate=bwa_missing_prob_err_rate,
			bwa_max_gap_opens=bwa_max_gap_opens,
			bwa_seed_length=bwa_seed_length,
		)

		if not valid_tracks:
			log("No valid tracks were found from provided tags and parameters. Export cannot proceed.", "E")
			raise typer.Exit(code=1)

		log(f"Specified tags: {', '.join(requested_tags)}", "I")

	log(f"Following tracks will be exported: {', '.join(valid_tracks)}", "I")

	try:
		merged_bed = generate_mask(
			ref_species=ref,
			inputs=valid_tracks,
			kmer_length=kmer_length,
			offset_step=offset_step,
			cross_stringency=cross_stringency,
			n_threads=n_threads,
			db_root=db_root,
			output_file=output_bed,
			bwa_missing_prob_err_rate=bwa_missing_prob_err_rate,
			bwa_max_gap_opens=bwa_max_gap_opens,
			bwa_seed_length=bwa_seed_length,
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
 

@app.command("import")
def import_tracks(
	ref: str = typer.Option(..., "-r", "--ref", help="Reference species/target name in the database."),
	query: str = typer.Option(..., "-i", "--input", help="Query/input species name for the track."),
	kmer_length: int = typer.Option(..., "-k", help="K-mer size used to generate the imported track."),
	offset_step: int = typer.Option(..., "-w", "--offset-step", "--sliding-window", help="Offset/sliding window used to generate the imported track."),
	mappability_all_bw: str = typer.Option(
		..., "-ma", "--mappability-all-bw", help="Path to mappability_all.bw generated externally."
	),
	mappability_uniq_bw: str = typer.Option(
		..., "-mu", "--mappability-uniq-bw", help="Path to mappability_uniq.bw generated externally."
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
		help="BWA aln -n used for generation.",
	),
	bwa_max_gap_opens: int = typer.Option(
		2,
		"-bo",
		help="BWA aln -o used for generation.",
	),
	bwa_seed_length: int = typer.Option(
		16500,
		"-bl",
		help="BWA aln -l used for generation.",
	),
	n_threads: int = typer.Option(
		1,
		"-j",
		help="Thread count used for generation.",
	),
	db_root: str = typer.Option(..., "-d", "--db-root", help="Path to the database root directory."),
	tag: Optional[List[str]] = typer.Option(
		None,
		"--tag",
		"-t",
		help="Tag(s) to associate to the imported track. Repeatable and comma-separated.",
	),
	force: bool = typer.Option(False, help="Overwrite imported track files/metadata if the track already exists."),
):
	"""Import a track manually by providing BigWig files and generation parameters."""
	if not valid_database(db_root):
		log(f"To initialize the database, run 'database --init {db_root}'", "E")
		raise typer.Exit(code=1)


	log("Importing track to database...", "I")
	try:
		result = import_track(
			ref_species=ref,
			query_species=query,
			kmer_length=kmer_length,
			offset_step=offset_step,
			mappability_all_bw=mappability_all_bw,
			mappability_uniq_bw=mappability_uniq_bw,
			db_root=db_root,
			input_fasta=input_fasta,
			reference_fasta=reference_fasta,
			reference_fasta_md5=reference_fasta_md5,
			tags=tag,
			mapping_tool=mapping_tool,
			bwa_missing_prob_err_rate=bwa_missing_prob_err_rate,
			bwa_max_gap_opens=bwa_max_gap_opens,
			bwa_seed_length=bwa_seed_length,
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
	log(f"mappability_all BigWig stored at: {result['mappability_all_bw']}", "S")
	log(f"mappability_uniq BigWig stored at: {result['mappability_uniq_bw']}", "S")
	raise typer.Exit(code=0)

if __name__ == "__main__":
	app()