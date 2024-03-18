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

## Usage

Basic usage :

      main.py -i fasta.fa

Different arguments can be added :

| Shortcut | Options | Description |  
|--|--|--|
| -c | y/n| If yes, run CAT  |  
| -n | x | The number of contigs to analyzed  |  
| -ml | x | The minimum length of the contigs to analyzed  |  
| -l | 1/2 | The level of prediction precision | 
| -d | y/n | If yes, the tools print their output in the stdout | 

By default, the tool uniformize the result of the tools to make easier the comparaison. It's possible to keep the original results with the level argument `-l 2`.

### About the benchmarking mode

As initially developped to benchmark the tools, some arguments can be added :

| Shortcut | Options | Description |  
|--|--|--|
| -tf | ... | The method to generate the file .tsv containing the truth | 
| -t | ... | The kingdom, the path or the directory  | 
| -s | ... | The searched kingdom to calculate the benchmarking scores | 

Three methods can be used in benchmarking mode :

 - `OK` for one-kingdom : the analyzed fasta contain only contigs from the `-t kingdom`.
 - `CAMI`: the analyzed fasta was generated with CAMISIM and the real taxonomy can be find in `-t gsa_mapping.tsv`
 - `GSA` : the analyzed fasta is an assembly that allowed the reconstruction of the whole genome of different organisms ; the files containing the GSA of these organisms are localized in the `-t directory`.

Note that in the `GSA` mode, the directory need to be build like that, where `kingdomX` need to be replaced by the real kingdom ("`eukaryote`", "`prokaryote`" or "`virus`") :

 - directory
	 - kingdom1
		 - GSA1.fa
		 - GSA2.fa
	- kingdom2
		- GSA3.fa
		- GSA4.fa
		- GSA5.fa

Keep in mind that in the case of the GSA mode, the assembled genomes need to be from the analyzed fasta and need to share the same contig id.

The benchmarking mode will automatically generate a `truth.tsv` with two colums : the contig id and the associated kingdom. It will then compares the `predictions.tsv` and the `truth.tsv` and calculates scores where the kingdom of `-s kingdom` is used to calculate the *True Positive* and *False Positive*.
