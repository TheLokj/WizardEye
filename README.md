# WizardEye

![Python](https://img.shields.io/badge/Python->=3.9-green.svg)
[![Version](https://img.shields.io/badge/version-0.1.4-yellow.svg)](https://github.com/TheLokj/WizardEye/releases)
![Beta](https://img.shields.io/badge/beta-orange.svg)
[![Install WizardEye using bioconda](https://img.shields.io/badge/Install%20WizardEye%20using-bioconda-brightblue.svg?style=flat)](https://anaconda.org/channels/bioconda/packages/wizardeye/overview)
[![Python CI](https://github.com/TheLokj/WizardEye/actions/workflows/test.yml/badge.svg)](https://github.com/TheLokj/WizardEye/actions/workflows/test.yml)
![GitHub bugs](https://img.shields.io/github/issues/TheLokj/WizardEye/bug)

WizardEye is a Python tool that filters aligned reads based on the risk of ambiguous alignment to a reference genome.

<img src="wizardeye.png" alt="WizardEye" width="700" />

**WizardEye is currently in beta and under active development. Future updates will bring new features, improve stability, and ensure robustness through comprehensive unit testing.**

To achieve this, WizardEye first identifies all positions in your reference genome that could be targeted by reads from known ambiguous sources using your alignment parameters. For example, it can filter out reads that map to regions conserved between your reference genome and potentially contaminating organisms.

The tool both creates a database of reference genomes and their associated risky sources and filters BAM files using that information. By directly using `bwa`, WizardEye empirically identifies risky reads based on your alignment parameters.

This tool is directly adapted from the script [generate_cross_mappability_filter_bwa.sh](https://github.com/TheLokj/generate_cross_mappability_filter/blob/master/BWA/generate_cross_mappability_filter_bwa.sh).

## Table of Contents

- [How it works](#how-it-works)
- [WizardEye limits](#wizardeye-limits)
- [Installation](#installation)
	- [Dependencies](#dependencies)
- [Usage](#usage)
	- [Create the database](#create-the-database)
	- [Create a new track](#create-a-new-track)
	- [Filter a BAM file](#filter-a-bam-file)
	- [Count and compute statistics](#count-and-compute-statistics)
	- [Export a mask](#export-a-mask)
	- [Import an existing track manually](#import-an-existing-track-manually)
- [Go beyond WizardEye limits](#go-beyond-wizardeye-limits)

## How it works

WizardEye first splits the potential contaminant source into `-k`-mers with a sliding window of `-w`. It then aligns each produced unique k-mer using `bwa aln` with your parameters on the target reference to highlight ambiguous regions of that sequence. As this step can be computationally intensive for a complete genome, it stores the computed cross-mappability track in a database.

You can then use these cross-mappability tracks to filter your alignment. For example, if you are studying the evolution of *Hominidae* and align your reads to the human genome, you can use WizardEye to remove reads that could also come from other mammalian sources, such as hyena or deer, using several cross-mappability tracks generated from non-Hominidae mammalian genomes.

## WizardEye limits

Since WizardEye relies only on alignment and prior knowledge, it is recommended for use when you have an idea of the potential contamination types in your sample. The more exhaustive the database, the more effective (and restrictive) the filtering.

Note that WizardEye only identifies ambiguous regions **that can align**. It **cannot** filter sequences that are very distant from the target. Therefore, it is highly recommended as a complementary step to traditional large-scale filtering. Finally, note that WizardEye filtering can be very strict, leading to the loss of substantial information, as it prioritizes risk avoidance over the statistical power of post-hoc analyses.

Additionally, WizardEye currently only supports single-end reads or merged pairs. Paired-end reads are not supported and will cause the tool to exit with an error.

## Installation

The latest WizardEye version can be installed using conda:

```
conda install -c bioconda wizardeye
```

You can also install WizardEye by cloning this repository and running the following commands from the main folder:

```
python3 -m venv .wizardeye
source .wizardeye/bin/activate
python -m pip install -U pip
python -m pip install -e .
```

### Dependencies

#### External tools

The following command-line tools must be installed and available in your `PATH`:

- `bwa` — for read alignment and reference indexing
- `samtools` — for SAM/BAM file manipulation (conversion, sorting, indexing, concatenation)
- `seqkit` — for k-mer generation, deduplication, and FASTA chunking
- `bedtools` — for interval operations (genomecov, merge)
- `bedGraphToBigWig` — UCSC tool for bedGraph to BigWig conversion (package: `ucsc-bedgraphtobigwig`)
- `bigWigToBedGraph` — UCSC tool for BigWig to bedGraph conversion (package: `ucsc-bigwigtobedgraph`)
- `parallel` — GNU parallel for parallel job execution
- `awk`, `cat`, `sort` — standard Unix utilities

#### Python packages

The following Python packages are automatically installed via pip or conda:

- `typer>=0.9` — for building the command-line interface
- `PyYAML>=6.0` — for YAML configuration file handling
- `pysam>=0.22` — for BAM file parsing and interval extraction
- `pyBigWig>=0.3.25` — for BigWig file parsing and k-mer overlap computation
- `numpy` — for numerical operations

## Usage

### Create the database

To start, create a WizardEye database using the following command:

```
wizardeye database init -d /path/to/database
```

This generates a `/database/` directory compatible with WizardEye. You can move it elsewhere or rename the folder after initialization.

#### About the database

The database is a directory containing subdirectories for targets and processed tracks.

```
database/
├── hg19/                           # reference
│   ├── sus_scrofa_k35_bwa8651fd65/ # a track aligned to the reference
│   │   ├── map_all.bw       	    # all overlapping k-mers
│   │   ├── map_uniq.bw             # unique overlapping k-mers
│   │   └── info.yaml               # information about the track
...
│   ├── ref.md5                     # reference md5
│   └── ref.yaml                    # information about the reference
└── info.yaml                      # information about the database
```

Every subdirectory contains two `BigWig` files that can be opened in a traditional genome browser. The `map_all.bw` file, for every position in the reference genome, represents the number of overlapping k-mers from the risky sequences, while `map_uniq.bw` contains only unique k-mers (i.e., k-mers overlapping only that genomic region). These two files enable stringency computation.

##### Update a track

You can update tags for an existing track using the database command by providing the full track-defining parameters and replacement tags:

```
wizardeye database update track-tags -d /path/to/database \
	-r hg19 -q bos_taurus -k 35 -w 1 -bn 0.01 -bo 2 -bl 16500 -bN false -bj 1 -bR 30 -bsn 2000000000 \
	--tags Mammalia,Ruminantia
```

##### Database cache

Note that masks generated during export are cached in the database to speed up subsequent exports. You can avoid this by specifying `--no-cache` before export. You can delete the cache using the following command:

```
wizardeye database clean -d /path/to/database
```

##### Migrate from version 0.1.2 to 0.1.3

Version 0.1.3 introduced new track naming. To migrate an existing database from an older version to 0.1.3, use the following command:

```
wizardeye database migrate -d /path/to/database --from 0.1.2 --to 0.1.3 -bR 30 -bsn 2000000000
```

### Create a new track

You can compute ambiguous regions of `reference.fa` that can be targeted by reads from `risky.fa` using specific BWA parameters:

```
wizardeye align -i /path/to/risky.fa -r /path/to/reference.fa -d /path/to/database \
					-k 35 -w 1 \
					-bn 0.01 -bo 2 -bl 16500
```

You can provide tags to describe `risky.fa` (e.g., phylogeny and/or environment: `-t Mammalia,Carnivora,Felis,Cave`). This enables direct BAM file filtering based on specific tags. You can manually set a track identifier with `--track_ID` to store it in the database metadata. Note that to prevent misuse, a new track cannot be added if the reference alignment file differs, including in sequence names. This is enforced by an MD5-based check to prevent silently corrupted analyses.

As WizardEye uses `seqkit split2` to distribute alignment computation across threads, several parameters can speed this up:

| Parameters | Default | Definition |
|------------|---|---------|
| `--tmp_dir` | `TMPDIR` | Temporary directory for processing chunks and alignments |
| `--chunk_size` | `2000000` | Number of sequences per chunk |
| `--jobs` | `1` | Number of threads for generating chunks and running parallel bwa instances |
| `--bwa_threads` | `1` | Number of threads per bwa instance |

> [!NOTE]
With default parameters (`-n 0.01 -o 2 -l 16500 -k 35 -w 1`), a cross-mappability track based on hg19 coordinates is generated from a complete mammalian genome in ~48 hours using 128 threads (`--chunk_jobs 128 --bwa_threads 1`). Such computation requires ~200GB of free storage in the temporary directory due to the large number of k-mers generated before alignment.

#### How to deal with exhaustivity

Several bwa parameters control alignment exhaustiveness and, consequently, track exhaustiveness.

| Parameters | Default | Definition |
|------------|---|---------|
| `-bR` | `30` | Number of hits to reach before stopping suboptimal hit search |
| `-bsn` | `2000000000` | Limit of alternative alignments in the SAM file |
| `-bN` | `False` | Should all alternative alignments be identified? |

The `-R` parameter is the most important. It defines the number of hits to reach before stopping the search for suboptimal hits. Hits are processed in score batches; the search only stops if this threshold is reached AND all hits with the same score have been processed. For example, with `-R=30` (the default for bwa and WizardEye), in a scenario with 50 best hits and 100 suboptimal hits, all 50 best hits are included in the final cross-mappability track, while the others are not.

Importantly, in classic `bwa aln` usage, the suboptimal search is also limited by the mismatch difference between the best hit and suboptimal hits. Thus, even with a large `-R`, the search identifies best hits and suboptimal hits with up to 1 mismatch compared to the best ones, regardless of `-n`.

If you want to avoid this behavior and be *truly* exhaustive in your cross-mappability track, you can specify the `-bN` parameter, which is passed as `-N` to `bwa aln`. This parameter forces BWA to save all alternative alignments, regardless of their quality and mismatches (within the `-n` (bwa_missing_prob_err_rate) limit). This parameter allows your cross-mappability track to be exhaustive at the expense of higher coverage. It can be used for extreme filtering but reduces the ability to retain informative reads. **Be careful with this option**: for complete genomes with repetitions and allowed mismatches, it can lead to biased representation due to `bwa aln`'s memory limit (alternative mappings are truncated). See [issue #3](https://github.com/TheLokj/WizardEye/issues/3) for more details.

Finally, the `bwa samse -n` parameter (accessible via `-bsn` in WizardEye) limits the number of alternative alignments retained in the SAM file. If there are more alternative alignments than this limit, bwa removes them all and sets MAPQ to 0. This value is set to a large number in WizardEye to avoid losing information and should not be modified.

> [!WARNING]
Due to BWA's behavior, coverage and intervals are not perfectly deterministic across systems and alignment chunk sizes. In fact, BWA randomly selects the primary alignment among best hits using a random number generator and deletes secondary alignments starting at the exact same position. This can lead to tiny coverage differences and, in some intervals, a 1 bp difference at the end position. However, track generation is consistent and deterministic if WizardEye is run on the same system with the same alignment chunk, regardless of the number of threads used. See [issue #1](https://github.com/TheLokj/WizardEye/issues/1) for more details.

### Filter a BAM file

Once you have a representative database, you can filter reads from your alignment:

```
wizardeye filter -i alignment.bam -r hg19 --exclude-tags Cave -k 35 -w 1 -bn 0.01 -bo 2 -bl 16500 -d /path/to/database
```

This filters out all reads that could be ambiguously aligned to your target if they come from a source with a tag in `--exclude-tags`.

For `filter`, `-r` must be a reference identifier already in your WizardEye database (e.g., `hg19`), and `-d/--db-root` is mandatory.

By default, WizardEye considers all reads that can align to the reference. However, to only consider reads aligning to a single position (`MAPQ>0`, no `XA` tag), specify `--only-unique`.

You can also filter out reads based on specific tracks:

```
wizardeye filter -i alignment.bam -r hg19 --exclude-tracks myotis_alcathoe,ursus_arctos -k 35 -w 1 -bn 0.01 -bo 2 -bl 16500 -d /path/to/database
```

> [!WARNING]
If sequence naming differs between BAM and tracks (e.g., `chr1` vs `1`), filtering stops with an explicit error. Harmonize contig names beforehand.

#### Adjust the filter hardness

##### Stringency

To balance sensitivity and specificity, you can specify a stringency value during filtering. This criterion is defined as follows.

For a position P in the target, there are exactly -k/-w different possible k-mers overlapping this position perfectly. After mapping, n k-mers can overlap the position (maximum -k/-w if no mismatches are allowed, otherwise more). The position is then highlighted as ambiguous if n/(-k/-w) >= -rc; i.e., if the proportion of overlapping k-mers relative to the total possible k-mers exceeds the stringency.

If `--only-unique` is used, the same behavior and formula apply, but only uniquely aligned k-mers are considered (i.e., k-mers without BWA's `XA` tag and with `MAPQ>0`). For example, with one mismatch allowed, `-k=40`, and `-rc=0.25`, if a position is overlapped by 10 exact k-mers and 10 k-mers with one mismatch, the position is retained only if 10 of these k-mers are unique, regardless of exactness. With `-rc=0.50`, the region is highlighted only if all 20 k-mers map uniquely to it.

```
wizardeye filter -i alignment.bam -r hg19 --exclude-tracks myotis_alcathoe,ursus_arctos -k 35 -w 1 -bn 0.01 -bo 2 -bl 16500 -p 0.25 -d /path/to/database --only-unique
```

Note that the ratio is not weighted by depth or mismatches. In another case with a repetitive region, if a position is overlapped by 3000 k-mers, the ratio remains the same.

This behavior aims to reproduce [Heng Li's seqbility](https://github.com/lh3/misc/tree/cc0f36a9a19f35765efb9387389d9f3a6756f08f/seq/seqbility) logic, which is not directly usable in a cross-mappability context.

##### Frequency

Using a frequency threshold, a read is filtered out only if a position is at least overlapped by minimum `-mf` tracks. This threshold is computed *after* the stringency. Then, a read is excluded if one of its position is covered by at least `-mf` tracks that exceeds the stringency threshold at that same position. You can apply such threshold using:

```
wizardeye filter -i alignment.bam -r hg19 --exclude-tracks myotis_alcathoe,ursus_arctos -mf 3 -d /path/to/database
```

Note that the `-mf` parameter is also available in the `export` command to generate masks with frequency filtering applied. This will remove lines with less than `-mf` from the final bed file.

#### Output

WizardEye produces a tabulation-separated report containing, for each read, the decision and overlapping tracks, as follows:

| read_key | filtered_out | associated_tracks |
|---|---|---|
| read_1:chrom1:1:35 | true | sus_scrofa,bos_taurus |
| read_2:chrom1:36:70 | false |  |
| read_3:chrom1:71:105 | true | felis_catus |
| read_4:chrom1:106:140 | true | bos_taurus |
| read_5:chrom1:232:256 | true | bos_taurus,ovis_aries |
| read_5:chrom2:2:36 | false |  |

With `--export-bam`, WizardEye produces two additional files:

- `excluded.bam`: reads excluded by the filtration,
- `filtered.bam`: reads retained by the filtration.

### Count and compute statistics

If you prefer to use your own filter, you can export a per-read report summarizing the number of k-mers overlapping the reference per read interval using the `count` command.

This command accepts the same parameters as `filter`, plus a `--mode` parameter to specify the statistical summary type (`sum`, `max`, `min`, `cov`, `mean`, or `std`). Statistics are computed per interval. For example, if `max` is specified, for each track, the maximum number of overlapping k-mers from that track at any position in the read's interval is reported.

```
wizardeye count -i alignment.bam -r hg19 --exclude-tags Farm -k 35 -w 1 -bn 0.01 -bo 2 -bl 16500 -d /path/to/database -m max
```

In this example, tracks not reported by `filter` with `-r=0.01` have a `0` in their column, meaning the maximum number of overlapping k-mers in the interval is 0.

#### Output

WizardEye produces a tabulation-separated report containing, for each read, the requested statistics:

| read_key | max_sus_scrofa | max_bos_taurus | max_ovis_aries |
|---|---|---|---|
| r1:chr1:2:36 | 0 | 4524 | 0 |
| r2:chr1:37:71 | 422 | 0 | 0 |
| r3:chr1:75:109 | 151 | 1254 | 1515 |
| r4:chr1:106:140 | 0 | 0 | 0 |
| r5:chr1:232:256 | 0 | 0 | 311 |
| r5:chr2:3:37 | 122 | 0 | 111 |

### Export a mask

If you plan to use the same configuration frequently (e.g., in a pipeline), you can export a mask to avoid recomputing it continuously:

```
wizardeye export -r hg19 --exclude-tags Cave -k 35 -w 1 -bn 0.01 -bo 2 -bl 16500 -d /path/to/database -o mask.bed
```

### Import an existing track manually

If you computed a track outside WizardEye, you can import it manually by providing both BigWig files and the parameters used to generate them:

```
wizardeye import -d /path/to/database -r ref -i input -k 35 -w 20 \
	--map-all-bw /path/to/map_all.bw \
	--map-uniq-bw /path/to/map_uniq.bw \
	--reference-fasta /path/to/reference.fa \
	--input-fasta /path/to/input.fa \
	-bn 0.01 -bo 2 -bl 16500 -j 8 \
	-t Mammalia,Carnivora
```

This command creates the target/track directory, copies the two BigWig files as `map_all.bw` and `map_uniq.bw`, and writes a `param.yaml` file with the provided generation metadata.

# Go beyond WizardEye limits

It is recommended to complement your filtering with an evolutionarily-aware method such as Kraken2. This combination is useful for removing both reads from completely different organisms and reads that may be ambiguous between closely related organisms.

*Last documentation update: 0.1.4.*
