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
		
    Block_size = 4

    Index_chunks = 2 

In the case where you don't plan to use CAT, you can only keep the `[user]` and `[environment]` section.

## Usage

Basic usage :

      main.py basic -i fasta.fa

Different arguments can be added :

| Shortcut | Options | Description |  
|--|--|--|
| -n | x | The number of contigs to analyze  |  
| -mnl | x | The minimum length of the contigs to analyze  |  
| -mxl | x | The maximum length of the contigs to analyze  |  
| -l | 1/2 | The level of prediction precision | 
| -c | y/n| If yes, run CAT  |  
| -cl | y/n | If yes, the tools output are cleaned after the exectution of WizardEye | 
| -d | y/n | If yes, the tools print their output in the stdout | 

By default, the tool uniformize the result of the tools to make easier the comparaison. It's possible to keep the original results with the level argument `-l 2`.

### About the benchmarking mode

#### Basic benchmark mode : compare the results of one launch with the truth

As initially developped to benchmark the tools, some arguments can be added to compare the predictions with the reality :

| Shortcut | Options | Description |  
|--|--|--|
| -t | ... | the path to the file containing the truth or the kingdom in One-Kingdom mode | 
| -tf | ... | The method to generate the file .tsv containing the truth | 
| -pr | ... | The kingdom to use as a positive reference to calculate scores | 

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

Finally, WizardEye in the benchmarking mode will compares the `predictions.tsv` and the `truth.tsv` and calculates scores where the kingdom of `-pr kingdom` is used to calculate the *True Positive*, *False Negative*. These results are saved at `performances.tsv`.

|  | EukRep | Tiara | Whokaryote | Whokaryote+Tiara | DeepMicroClass |
|--|--|--|--|--|--|
| TP | 39 | 40 | 37 | 40 | 37 |
| TN | 158 | 158 | 137 | 154 | 153 |
| FP | 3 | 3 | 24 | 7 | 2 |
| FN | 1 | 0 | 3 | 0 | 9 |

Each generated `truth.tsv`,  `predictions_metrics.tsv` and tools output can be found in the directory output/`newName`/ where the `newName` is a name build with the original input name, the number of contigs, the minimum selected length, the maximum selected length and the wanted proportion as described below. 

##### Use kingdom-specific proportions

Note that you can also add proportion parameters in this mode to select only a precise percentage of each specified kingdom :

| Shortcut | Options | Description |  
|--|--|--|
| -kl | kingdom1:kingdon2 | The name of the wanted kingdoms   | 
| -kp | X:Y | The proportion of each kingdom | 

Each kingdom and proportion need to be separated with `:`.

Note that it's also possible to precise subdivision proportion with `..`, as following :

    -kl eukaryote..fungi.chlorophyta..:prokaryote -kp 30..20.10..:70 

In this case, the subdivision need to be precised in a specific column of the `truth.tsv` :

| contig ID | division | subdivision |
|--|--|--|
| ERZXXXXX.1 | eukaryote | fungi |
| ERZXXXXX.23 | unknown |  |
| ERZXXXXX.3432 | prokaryote |  |
| ERZXXXXX.434 | eukaryote | chlorophyta |
| ERZXXXXX.3 | unknown |  |

#### Expert benchmark mode : launch WizardEye several times adjusting the parameters at each iteration 

    main.py expert -i fasta.fa

Because a benchmark is better with multiple data, WizardEye also provide an expert mode to launch the tools several times. This mode resumes all the predictions and metrics obtained thanks to WizardEye and the tools in several `.tsv` files. 

| Shortcut | Options | Description |  
|--|--|--|
| -na | name | The Benchmark name    | 
| -mo | mode | The Benchmark mode (see below) | 

##### Launching WizardEye with different contig lengths 

The `length` mode allows WizardEye to be run several times, selecting contigs between a minimum length `-mnL` and a maximum length `-mxl`. An additional argument, `-lr`, is required to specify the size of the different intervals.

    python3 main.py expert -i example.fa -na LENGTH -mo length -t truth.tsv -c n -lr 1000 -mnl 3000 -mxl 10000 -n 200 -pr eukaryote

In this example, WizardEye will be launched 7 times, analysing firstly 200 contigs from 3000 pb to 4000 pb, then 200 contigs from 4000 pb to 5000 pb... The positive refence to calculate the scores is `eukaryote`.

##### Launching WizardEye with different kingdom known proportions
 
The `proportion` mode allows WizardEye to be run several times, selecting contigs of `-ekl` kindgoms separated with `:` whose proportions are changing. The proportions are precised with the `-ekp` argument and each iteration is separated with a `,`, the `:` separing once again the kingdom.
 
    python3 main.py expert -i example.fa -na PROP -mo proportion -t truth.tsv -ekl eukaryote:prokaryote -ekp 20,80,50,40:80,20,50,60 -mnl 3000 -mxl 10000 -n 200 -pr eukaryote

In this example, WizardEye will be launched 4 times, analysing firstly 200 contigs between 3000 and 10000pb with 20% of eukaryotic contigs and 80% of prokaryotic contigs, then 200 contigs between 3000 and 10000pb with 80% of eukaryotic contigs and 20% of prokaryotic contigs... The positive refence to calculate the scores is `eukaryote`.

#### Output of the benchmarking mode

The benchmarking mode produces different files summarizing :

- the parameters of the run ;
 - the predictions per tool per contigs and the associated truth ;
 - the accuracy, the F1-score and the MCC score representing the predictions correctness calculated thanks to the positive reference `-pr`  ;
 - the running time (the real and the CPU's ones) per tool.

When WizardEye is launched in the expert mode, each of these files are summarized in bigger files containing these results and the ID of the associated run (basically the length range in `length` mode and the proportion in `proportion` mode). It also produce an empty `.mode` file used by the R function to detect the benchmark nature. 

#### Plot the results with the R script

You can use the script `plot.R` to easily plot the results of an expert benchmark made by WizardEye. 

As it isn't currently directly part of WizardEye, you will need to move by yourself the working directory to the output directory of your choice, by using the R `setwd(path)` command.

Note that this script need the following librairies :

    ggplot2, ggubr, stringr, tidyverse, plyr, reshape2
