# WizardEye 

## Requirements

    biopython, numpy, pandas
  
As WizardEye does not contain the tools, it's also necessary to install them before trying to use it. WizardEye uses the `conda` environment to launch `EukRep`, `Tiara`, `Whokaryote` and `DeepMicroClass`.

The tool `CAT` needs for the moment to be installed in a singularity image (.sif). `CAT` remains unnecessary and its use can be disabled, so you can choose to avoid its installation.

## Installation

WizardEye can be quickly used by cloning the repository with git.

    git clone https://github.com/TheLokj/WizardEye.git
  
## Configuration

Before launching WizardEye, you will need to configure it in the `config.ini` file. This configuration is necessary notably to precise the name of the conda environments where the tools are installed, an email to use the NCBI Entrez services and the details about CAT.

    [user]
    
    mail = example@mail.com
    
    [environment]
    
    EukRep = eukrep
    
    Tiara = tiara
    
    Whokaryote = whokaryote
    
    DeepMicroClass = dmc
    
    [CAT]
    
    CAT_sif = /CAT.sif
    
    CAT_DB = /CAT_database/
    
    CAT_Diamond_DB = /CAT_prepare/Diamond_X.X.X/diamond
    
    CAT_Taxonomy_DB = /CAT_taxonomy/

In the case where you don't plan to use CAT, you can only keep the `[user]` and `[environment]` section.

## Usage

Basic usage :

      main.py -i fasta.fa

Different arguments can be added :

| Shortcut | Options | Description |  
|--|--|--|
| -c | y/n| If yes, run CAT  |  
| -n | x | The number of contigs to analyze  |  
| -mnl | x | The minimum length of the contigs to analyze  |  
| -mxl | x | The maximum length of the contigs to analyze  |  
| -l | 1/2 | The level of prediction precision | 
| -d | y/n | If yes, the tools print their output in the stdout | 

By default, the tool uniformize the result of the tools to make easier the comparaison. It's possible to keep the original results with the level argument `-l 2`.

### About the benchmarking mode

As initially developped to benchmark the tools, some arguments can be added :

| Shortcut | Options | Description |  
|--|--|--|
| -bm | ... | The method to generate the file .tsv containing the truth | 
| -t | ... | the path to the file containing the truth or the kingdom in One-Kingdom mode | 
| -bs | ... | The searched kingdom to calculate the benchmarking scores | 

Three methods can be used in benchmarking mode :

 - `CAMI`: the input fasta was generated with CAMISIM and the real taxonomy can be find in `-t gsa_mapping.tsv`.
 -  `AK` for **A**lready-**K**nown : a file containing the real taxonomy for each contig already exists and is located to`-t path`. It need to respect the format described below.
 - `OK` for **O**ne-**K**ingdom : the inputfasta contain only contigs from the `-t kingdom`.

The `OK` and `CAMI` modes automatically generate a `truth.tsv` with two colums : the contig id and the associated kingdom. Note that in the case where you want to use your own-made `truth.tsv` with the `AK` mode and where your input fasta contain also non-labelled contigs, it's necessary to label them as `unknown`.

| contig ID | division |
|--|--|
| ERZXXXXX.1 | eukaryote |
| ERZXXXXX.23 | unknown |
| ERZXXXXX.3432 | prokaryote |
| ERZXXXXX.444 | prokaryote |
| ERZXXXXX.3 | unknown |

Finally, WizardEye in the benchmarking mode will compares the `predictions.tsv` and the `truth.tsv` and calculates scores where the kingdom of `-s kingdom` is used to calculate the *True Positive*, *False Negative*. These results are saved at `performances.tsv`.

|  | EukRep | Tiara | Whokaryote | Whokaryote+Tiara | DeepMicroClass |
|--|--|--|--|--|--|
| TP | 39 | 40 | 37 | 40 | 37 |
| TN | 158 | 158 | 137 | 154 | 153 |
| FP | 3 | 3 | 24 | 7 | 2 |
| FN | 1 | 0 | 3 | 0 | 9 |

Each generated `truth.tsv`,  `performances.tsv` and tools output can be found in the directory output/`newName`/ where the `newName` is a name build with the original input name, the number of contigs, the minimum selected length, the maximum selected length and the wanted proportion as described below. 

### Use kingdom-specific proportions

Note that you can also add proportion parameters in this mode to select only a precise percentage of each specified kingdom :

| Shortcut | Options | Description |  
|--|--|--|
| -bk | kingdom1:kingdon2 | The name of the kingdoms   | 
| -bp | X:Y | The proportion of each kingdom | 

Each kingdom and proportion need to be separated with `:`.




 

 

