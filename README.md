# WizardEye

![Python](https://img.shields.io/badge/Python-green.svg)
![Beta](https://img.shields.io/badge/beta-red.svg)
[![Version](https://img.shields.io/badge/version-0.1.0b0-orange.svg)](https://github.com/TheLokj/WizardEye/releases)
[![Python CI](https://github.com/TheLokj/WizardEye/actions/workflows/test.yml/badge.svg)](https://github.com/TheLokj/WizardEye/actions/workflows/test.yml)
![GitHub bugs](https://img.shields.io/github/issues/TheLokj/WizardEye/bug)

WizardEye is a Python tool that filters aligned reads according to the risk of ambiguous alignment sources on a reference genome. 

<img src="wizardeye.png" alt="WizardEye" width="700" />

**WizardEye is currently in beta and under active development. Future updates will bring new features, improve stability, and ensure robustness through comprehensive unit testing.**

To do this, WizardEye first identifies all positions in your reference genome that can be targeted by reads from known ambiguous sources using your alignment parameters. For example, it can be used to filter out reads in your alignment that map to regions conserved between your reference genome and several potentially contaminating organisms.

The tool is designed both to create a database containing several reference genomes and their associated risky sources, and to filter BAM files according to that stored information. By directly using `bwa`, WizardEye can empirically identify risky reads based on your frequently used alignment parameters.

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

As WizardEye is based only on alignment and prior knowledge, it is recommended to use it when you have an idea of the types of contamination that may be present in the sample. The more exhaustive the database, the more efficient (and restrictive) the filtering.

Please also note that WizardEye is only designed to identify ambiguous regions **that can align**. This tool **cannot** be used to filter out sequences that are very distant from the target. It is therefore highly recommended to use it as a complementary step to traditional large-scale filtering. Finally, note that WizardEye filtering can be very strict and lead to the loss of substantial information, prioritizing risk avoidance over the statistical power of post-hoc analyses.

## Installation

You can install WizardEye by cloning this repository and by running the following commands from the main folder:

```
python3 -m venv .wizardeye
source .wizardeye/bin/activate
python -m pip install -U pip
python -m pip install -e .
```

### Dependencies

WizardEye relies on both Python packages and external command-line tools.

### External tools (called with subprocess)

Required:

- `bwa` (alignment and indexing)
- `samtools` (SAM/BAM conversion, sorting, concatenation)
- `seqkit` (k-mer generation, deduplication, chunking)
- `bedGraphToBigWig` (UCSC tool, package often named `ucsc-bedgraphtobigwig`)
- `bigWigToBedGraph` (UCSC tool, package often named `ucsc-bigwigtobedgraph`)
- `bedtools` (faster interval merging in mask generation)
- standard Unix tools: `awk`, `sort`, `cat` (used in high-performance paths)

### Python packages

Required:
- `typer` (CLI)
- `PyYAML` (database and track metadata YAML files)
- `pysam` (BAM parsing and interval extraction)
- `pyBigWig` (BigWig parsing and overlapping k-mer computation)

## Usage

### Create the database

You should initially create a WizardEye database using the following command:

```
wizardeye database --init -d /path/to/database
```

This will generate a `/database/` directory compatible with WizardEye. Feel free to move it elsewhere or rename the folder after this initialisation. 

#### About the database

The database is a directory containing sub-directories representing targets and already processed tracks.

```
database/
├── hg19/                           # reference
│   ├── sus_scrofa_k35_w1_n1_o2_l1/ # a track divided and aligned on the reference
│   │   ├── map_all.bw      	    # all overlapping k-mers
│   │   ├── map_uniq.bw             # unique overlapping k-mers
│   │   └── info.yaml               # information about the track
...
│   ├── ref.md5                     # reference md5
│   └── ref.yaml                    # information about the reference
└──  info.yaml                      # information about the database
```

Every subdirectory contains two `bigWig` files that can be opened in a traditional genome browser. The `map_all.bw` file represents, for every position of the `reference.fa` genome, the number of overlapping k-mers from the `risky` sequences, while `map_uniq.bw` contains only unique k-mers, i.e., k-mers that overlap only this area of the genome. These two files allow you to compute stringency.

##### Update a track

You can update tags of an existing track from the database command by providing the full track-defining parameters and the replacement tags:

```
wizardeye database --update-track-tags -d /path/to/database \
	-r hg19 --track bos_taurus -k 35 -w 1 -bn 0.01 -bo 2 -bl 16500 \
	--tags Mammalia,Ruminantia
```

##### Database cache

Note that masks generated during exportation are cached in the database, to make next exportation quicker. You can avoid that specifying `--no-cache` before exportation. You can delete the cache using the following command:

```
wizardeye database --clean -d /path/to/database
```

### Create a new track

You can compute ambiguous regions of `reference.fa` that can be targeted, using specific BWA parameters, by reads from `risky.fa` with:

```
wizardeye align -i /path/to/risky.fa -r /path/to/reference.fa -d /path/to/database \ 
					-k k_mer_length -w sliding_windows \
					-bn bwa_missing_prob_err_rate -bo bwa_max_gap_opens -bl bwa_seed_length
```

Note that you can provide tags here to describe `risky.fa`. This can be used to describe the sequence, for example by specifying phylogeny and/or a particular environment: `-t Mammalia,Carnivora,Felis,Cave`. This is useful to filter a `.bam` file directly according to specific tags. You can manually set a track identifier with `--track_ID`, in order to save it in database metadata. Note also that, to prevent misuse, it is not possible to add a new track to the database if the reference alignment file is not exactly the same, including sequence names. This is due to an MD5-based control to avoid silently corrupted analyses.

As WizardEye uses `seqkit split2` to share the computation of the alignment between threads, several parameters can be used to make the computation quicker:

| Parameters | default |  |
|------------|---|---------|
| `--tmp_dir` | `TMPDIR` | Temporary directory used to process chunks and alignments |
| `--chunk_size` | 2000000 | Number of sequences per chunk |
| `--jobs` | 1 | Number of threads used to generate chunks and number of parallel bwa instances processing chunks |
| `--bwa_threads` | 1 | Number of threads allowed per bwa instances |

> [!NOTE]  
Using default parameters (`-n 0.01 -o 2 -l 16500 -k 35 -w 1`), a cross-mappability track based on hg19 coordinates is generated from a complete mammalian genome in ~48 hours using 128 threads (`--chunk_jobs 128 --bwa_threads 1`). Such computation requires ~200GB of free storage in the temporary directory, due to the huge amount of k-mers before alignment.

> [!WARNING]  
Due to BWA behaviour, it is important to note that coverage and intervals are not perfectly determinist between systems and alignment chunk sizes. In fact, BWA is selecting randomly the primary alignment among best hits based on a random number generator and delete secondary alignment starting at the exact same position. This can leads to tiny differences of coverage and, in some intervals, to a difference of 1bp at end position. However, tracks generation is consistent and deterministic if WizardEye is run on the same system with the same alignment chunk, despite the number of threads used. See [issue #1](https://github.com/TheLokj/WizardEye/issues/1) for more details. 

If you want an exhaustive cross-mappability track, you can also specify the `-bN` parameter, which will be passed as `-N` in `bwa aln`. This parameter will force BWA to save every alternative alignments, regardless their qualities and mismatches (within the limit of `-bn bwa_missing_prob_err_rate`). This parameter allows your cross-mappability track to be exhaustive in exchange of highter coverage. It can be used for extreme filtration, but it will reduce your ability to keep informative reads.

### Filter a BAM file

Once you have a representative database, you can use it to filter out reads from your alignment. To do that, you can enter:

```
wizardeye filter -i alignment.bam -r hg19 --exclude-tags Cave -k 35 -w 1 -bn 0.01 -bo 2 -bl 16500 -d /path/to/database
```

This will filter out all reads that can be ambiguously aligned to your target if they come from a source with a tag listed in `--exclude-tags`. 

For `filter`, `-r` must be the reference identifier already present in your WizardEye database (for example `hg19`), and `-d/--db-root` is mandatory.

By default, WizardEye only considers reads that can be aligned with a single position in the reference. However, if you want to exclude reads that can be aligned with multiple positions, and then be even more restrictive, specify the additional parameter `--considere-all`.

You can also filter out reads based on specific tracks:

```
wizardeye filter -i alignment.bam -r hg19 --exclude-tracks myotis_alcathoe,ursus_arctos -d /path/to/database
```

> [!WARNING]  
If sequence naming differs between BAM and tracks (for example `chr1` vs `1`), filtering stops with an explicit error. Harmonize contig names beforehand.

#### Adjust the filter hardness

##### Stringency 

To handle the trade-off between sensitivity and specificity, you can specify a stringency value during filtering. This criterion is defined as described below.

For a position $P$ in the target, there are exactly `-k`/`-w` different possible k-mers that overlap that position perfectly. After mapping, among the `n` k-mers overlapping the position (maximum `-k`/`-w` if no mismatches are allowed, otherwise more), `u` is computed as the number of these k-mers that are unique in the target. The position is then highlighted as ambiguous if `u`/(`-k`/`-w`) >= `-rc`, i.e., if the proportion of unique k-mers overlapping the position compared with the total number of possible k-mers overlapping the position is higher than the stringency.

For example, in a situation with one mismatch allowed, with `-k=40` and `-rc=0.25`, if a position is overlapped by 10 exact k-mers and 10 k-mers with one mismatch, the position is retained only if 10 of these k-mers are unique, regardless of exactness. In a similar case with `-rc=0.50`, the region is highlighted only if all 20 k-mers map uniquely there.

```
wizardeye filter -i alignment.bam -r hg19 --exclude-tracks myotis_alcathoe,ursus_arctos -r 0.25 -d /path/to/database
```

Note that the ratio is not weighted by depth. In another case with a repetitive region, where a position is overlapped by 3000 k-mers, the ratio is still the same: if 40 of these are unique, the position is highlighted.

Unique k-mers are defined as k-mers without BWA `XA` tag and with `MAPQ>0`. 

This behavior aims to reproduce the [Heng Li's seqbility](https://github.com/lh3/misc/tree/cc0f36a9a19f35765efb9387389d9f3a6756f08f/seq/seqbility) logic, which is not directly usable in a cross-mappability context. 

If `--considere-all` is used, the same behavior and formulae are still used but uniqueness is simply no longer required.

##### Frequency

*[Not implemented yet]*

The position will only be masked if at least `-mf` tracks overlap that position.

```
wizardeye filter -i alignment.bam -r hg19 --exclude-tracks myotis_alcathoe,ursus_arctos -f 3 -d /path/to/database
```

If the `-by_tags` parameter is used, frequency is computed per tag rather than per track. A position is then masked only if there is overlap from tracks belonging to the specified tags, where `-f` represents the number of tags that must be represented. For example, here, a read is masked only if it has a position overlapped by a combination of 2 of the specified clades.

```
wizardeye filter -i alignment.bam -r hg19 --exclude-tracks myotis_alcathoe,ursus_arctos -f 2 -by_tags Carnivora,Ruminant -d /path/to/database
```

#### Output

WizardEye produces a tabulation-separated report containing, for each read, the decision and the overlapping tracks and tags, such as follows: 


| read_id                                  | excluded | overlapped | tags |
|------------------------------------------|----------|------------|------|
| read_1								   | true  	  | sus_scrofa,bos_taurus | Suina,Ruminantia |
| read_2								   | false    |            |      |
| read_3								   | true     | felis_catus | Feliformia     |
| read_4								   | true     | bos_taurus | Ruminantia     |
| read_5								   | true    | bos_taurus,ovis_aries           | Ruminantia     |
| read_6								   | false    |          |      |

With `--export-bam`, WizardEye produces two more files:

- a BAM `excluded` with the reads excluded by the filtration,
- a BAM `filtered` with the reads kept by the filtration,

### Count and compute statistics

If you prefer to use your own-made filter, you can export a per-read report which summarizes the number of k-mers which can overlap the reference per read-interval using the `count` command. 

This command take the same parameters as `filter`, plus a parameter `--mode` to specify the type of statistical summary you want betweemn `sum`, `max`, `min`, `cov`, `mean` and `std`. Statistics are computed on interval. For example, if `max` is specified, for each track, the maximum number of overlapping k-mers from this track on a position in the read associated interval is reported. 

```
wizardeye count -i alignment.bam -r hg19 --exclude-tags Farm -k 35 -r 1 -bn 0.01 -bo 2 -bl 16500 -d /path/to/database -m max
```

In this example, track whose are not reported using `filter`+`-r=0.01` should have a `0` in their column, as it means that the maximum number of overlapping k-mers in the interval is 0. 

#### Output

WizardEye produces a tabulation-separated report containing, for each read, the requested : 


| read_id                                  | max_sus_scrofa | max_bos_taurus_max | max_ovis_aries |
|------------------------------------------|----------|------------|------|
| read_1								   | 0 | 4524 | 0 |
| read_2								   | 422 | 0 | 0 |
| read_3								   | 151 | 1254 | 1515 |
| read_4								   | 0 | 0 | 0 |
| read_5								   | 0 | 0 | 311 |
| read_6								   | 122 | 0 | 111 |

### Export a mask

If you plan to use the same configuration often to filter your reads, for example in a pipeline, it is possible to export a mask in order to avoid recomputing it continuously:

```
wizardeye export -r hg19 --exclude-tags Cave -k 35 -w 1 -bn 0.01 -bo 2 -bl 16500 -d /path/to/database -o mask.bed
```

### Import an existing track manually

If you already computed a track outside WizardEye, you can import it manually by providing both BigWig files and the parameters that were used to generate them:

```
wizardeye import -d /path/to/database -r ref -i input -k 35 -w 20 \
	--map-all-bw /path/to/map_all.bw \
	--map-uniq-bw /path/to/map_uniq.bw \
	--reference-fasta /path/to/reference.fa \
	--input-fasta /path/to/input.fa \
	-bn 0.01 -bo 2 -bl 16500 -j 8 \
	-t Mammalia,Carnivora
```
This command creates the target/track directory, copies the two BigWig files as `map_all.bw` and `map_uniq.bw`, and writes a `param.yaml` with the provided generation metadata.

# Go beyond WizardEye limits

It is recommanded to complete your filtration using an evolutionnary-aware method such as Kraken2. This combination is useful to remove both reads that belong to completely different organisms and reads that can be ambiguous between closely related organisms.

*Last update of this documentation: beta-0.1.0.*