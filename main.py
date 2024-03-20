# -*- coding: utf-8 -*-
import os, sys, configparser, argparse
import pandas as pd
import subprocess as sp
from Bio import SeqIO
from modules import fastatools as ft
from modules import truthchecker as tc
from modules import scores as sc

# Parse arguments
parser = argparse.ArgumentParser(description="WizardEye summarize kingdom identification from contigs thanks to different tools specialised in the detection of eukaryotic contigs")
parser.add_argument('-i', '--input', metavar='path', type=str, required=True, help="the input multifasta containing the contigs")
parser.add_argument('-c', '--cat', metavar='y/n', type=str, choices=["y", "n"], help="use or not CAT (pretty long step)")
parser.add_argument('-n', '--ncontigs', metavar="n", type=int, help='the number of contigs to select')
parser.add_argument('-mnl', '--minlength', metavar="n", type=int, default=1, help='the minimun length of the analyzed contigs in pb')
parser.add_argument('-mxl', '--maxlength', metavar="n", type=int, help='the maximum length of the analyzed contigs in pb')
parser.add_argument('-l', '--level', metavar="1/2", type=int, choices=[1, 2], default=1, help='the level of precision of the summary (1 : uniformized : prokaryote, eukaryote, virus ; 2 : raw output of tools)')
parser.add_argument('-d', '--debug', metavar="y/n", type=str, choices=["y", "n"], help='the debug mode, print everything in the stdout')
parser.add_argument('-bm', '--bmode', metavar="format", type=str, choices=["CAMI", "OK", "AK"], help='the benchmarking mode : CAMI, OK (One-Kingdom), AK')
parser.add_argument('-t', '--truth', metavar="...", type=str, help='the path containing the real taxonomy (or the kingdom in One-Kingdom mode)')
parser.add_argument('-bp', '--bproportion', metavar="X:Y", type=str, help='the proportion of each kingdom where')
parser.add_argument('-bk', '--bpropkingdom', metavar="kingdom1:kingdom2:kingdom3", type=str, help='the name of the kingdom use with -bp')
parser.add_argument('-bs', '--bsearch', metavar="kingdom", type=str, help='the searched kingdom to calculate the benchmarking scores')

# Read arguments
args = parser.parse_args()
fastaPath = args.input
ncontigs = args.ncontigs
minLength = args.minlength
maxLength = args.maxlength
precisionLevel = args.level
debugMode = args.debug
benchMode = args.bmode
benchParam = args.truth
benchQuery = args.bsearch
benchProp = args.bproportion
benchKProp = args.bpropkingdom

print("\n▄█     █▄   ▄█   ▄███████▄     ▄████████    ▄████████ ████████▄     ▄████████ ▄██   ▄      ▄████████\n███     ███ ███  ██▀     ▄██   ███    ███   ███    ███ ███   ▀███   ███    ███ ███   ██▄   ███    ███ \n███     ███ ███▌       ▄███▀   ███    ███   ███    ███ ███    ███   ███    █▀  ███▄▄▄███   ███    █▀ \n███     ███ ███▌  ▀█▀▄███▀▄▄   ███    ███  ▄███▄▄▄▄██▀ ███    ███  ▄███▄▄▄     ▀▀▀▀▀▀███  ▄███▄▄▄     \n███     ███ ███▌   ▄███▀   ▀ ▀███████████ ▀▀███▀▀▀▀▀   ███    ███ ▀▀███▀▀▀     ▄██   ███ ▀▀███▀▀▀     \n███     ███ ███  ▄███▀         ███    ███ ▀███████████ ███    ███   ███    █▄  ███   ███   ███    █▄  \n███ ▄█▄ ███ ███  ███▄     ▄█   ███    ███   ███    ███ ███   ▄███   ███    ███ ███   ███   ███    ███ \n▀███▀███▀  █▀    ▀████████▀   ███    █▀    ███    ███ ████████▀    ██████████  ▀█████▀    ██████████\n                                                                                                 v0.0.1 by TheLokj")
print("\nMore information and updates on www.github.com/TheLokj/WizardEye/\n")

# Check if the input and output directories exist
if os.path.exists("input") is False :
    os.mkdir("input")
if os.path.exists("output") is False :
    os.mkdir("output")

# Print the tools output in the stdout or not according to the debug mode
if debugMode == "y":
    toolsOutput = ""
else :
    toolsOutput = "> /dev/null"

# Read config file
if os.path.exists("config.ini") is False :
    print("Error : the configuration file config.ini does not exist.")
    sys.exit()
else :
    config = configparser.ConfigParser()
    config.read("config.ini")

# Format the final name of the fasta
fastaName = f"{os.path.basename(fastaPath).split(".")[0]}_{ncontigs}seq_{minLength}pb"
if maxLength != None :
    fastaName = fastaName[:-2] + f"to{maxLength}pb"
if benchProp != None and benchKProp != None :
    for i in range(len(benchProp.split(":"))) :
        fastaName += "_" + benchProp.split(":")[i] + benchKProp.split(":")[i]

if os.path.exists(f"output/{fastaName}") is False :
    os.mkdir(f"output/{fastaName}")

# If the truth is already known, get the contig ids where the taxonomy is known to only use these ids in the newFasta
if benchMode == "AK":
    print("\nReading real taxonomies...")
    truthPath, knownIds = benchParam, {}
    truth = open(benchParam)
    for line in truth.readlines()[1:] :
        if line.split("\t")[1].replace("\n", "") != "unknown":
            if line.split("\t")[1] not in knownIds.keys() :
                knownIds[line.split("\t")[1]] = []
            knownIds[line.split("\t")[1]].append(line.split("\t")[0])
    truth.close()
    if benchProp != None and benchKProp != None :
        proportions = {}
        for i in range(len(benchKProp.split(":"))) :
            proportions[benchKProp.split(":")[i]] = benchProp.split(":")[i]
        if sum([float(value) for value in proportions.values()]) != float(100):
            print("The sum of the proportions should be 100% !")
            sys.exit()
    else :
        proportions = None
else :
    knownIds, proportions = None,  None

# Make a copy of the multifasta with the selectionned size, n contigs
newFasta  = ft.split_multifasta(fastaPath, fastaName, "input/", nSeq=ncontigs, minLength=minLength,  maxLength=maxLength, knownIds=knownIds, proportions=proportions)

# Format the real taxonomy
if benchMode == "CAMI":
        print("\nReading real taxonomies...")
        mail = config['user']['mail']
        truthPath = tc.get_cami_tax(newFasta, benchParam, f"output/{fastaName}/", mail)
if benchMode == "OK":
    print("\nReading real taxonomies...")
    truthPath = tc.get_oneKingdom_tax(newFasta, benchParam, f"output/{fastaName}/")

# Run tools
print("\nRunning tools...")
if args.cat == "y" :
    if os.path.exists(f"output/{fastaName}/CAT") is False :
        print("Running CAT...")
        os.mkdir(f"output/{fastaName}/CAT")
        os.system(f'singularity run {config["CAT"]["CAT_sif"]} CAT contigs -c {newFasta} -d {config["CAT"]["CAT_DB"]} --path_to_diamond {config["CAT"]["CAT_Diamond_DB"]} -t {config["CAT"]["CAT_Taxonomy_DB"]} --out_prefix "./output/{fastaName}/CAT/out.CAT" --block_size 5 --index_chunks 1 {toolsOutput}')
        os.system(f'singularity run {config["CAT"]["CAT_sif"]} CAT add_fastaNames -i "./output/{fastaName}/CAT/out.CAT.contig2classification.txt" -o "./output/{fastaName}/CAT/out.CAT.lineage" -t {config["CAT"]["CAT_Taxonomy_DB"]} {toolsOutput}')
    else :
        print("CAT predictions found...")
else :
    print("CAT will not be launched.")

if os.path.exists(f"output/{fastaName}/eukrep") is False :
    os.mkdir(f"output/{fastaName}/eukrep")
    print("Running EukRep...")
    os.system(f'eval "$(conda shell.bash hook)" && conda activate {config["environment"]["EukRep"]} && EukRep -i {newFasta} -o output/{fastaName}/eukrep/eukrep.txt {toolsOutput}')
else :
    print("EukRep predictions found...")

if os.path.exists(f"output/{fastaName}/tiara") is False :
    os.mkdir(f"output/{fastaName}/tiara")
    print("Running Tiara...")
    os.system(f'eval "$(conda shell.bash hook)" && conda activate {config["environment"]["Tiara"]} && tiara -i {newFasta} -o output/{fastaName}/tiara/tiara.txt {toolsOutput}')
else :
    print("Tiara predictions found...")

if os.path.exists(f"output/{fastaName}/whokaryoteS") is False :
    print("Running Whokaryote...")
    os.system(f'eval "$(conda shell.bash hook)" && conda activate {config["environment"]["Whokaryote"]} && whokaryote.py --contigs {newFasta} --outdir output/{fastaName}/whokaryoteS --model S --minsize {minLength} {toolsOutput}')
else :
    print("Whokaryote predictions found...")

if os.path.exists(f"output/{fastaName}/whokaryoteT") is False :
    print("Running Whokaryote+Tiara...")
    os.system(f'eval "$(conda shell.bash hook)" && conda activate {config["environment"]["Whokaryote"]} && whokaryote.py --contigs {newFasta} --outdir output/{fastaName}/whokaryoteT --model T --minsize {minLength} {toolsOutput}')
else :
    print("Whokaryote+Tiara predictions found...")

if os.path.exists(f"output/{fastaName}/dmc") is False :
    print("Running DeepMicroClass...")
    os.system(f'eval "$(conda shell.bash hook)" && conda activate {config["environment"]["DeepMicroClass"]} && DeepMicroClass predict -i {newFasta} -o output/{fastaName}/dmc/ {toolsOutput}')
    os.remove("model.ckpt")
else :
    print("DeepMicroClass predictions found...")

# Save the predictions
eukOutput = open(f"output/{fastaName}/eukrep/eukrep.txt", "r")
eukResults = eukOutput.read()
eukOutput.close()
predictions = {}
for record in SeqIO.parse(newFasta, "fasta"):
    predictions[record.id] = {}
    # Checking CAT predictions
    if args.cat == 'y' :
        try:
            predictions[record.id]["CAT"] = sp.check_output([f'grep -P "{record.id+"\t"}" output/{fastaName}/CAT/out.CAT.lineage | grep "(superkingdom)"'],shell=True).decode("UTF-8").split("(superkingdom)")[0].split("\t")[-1]
        except:
            predictions[record.id]["CAT"] = None
    # Checking EukRep predictions
    if str(record.id) in eukResults :
        predictions[record.id]["EukRep"] = 'eukaryote'
    else :
        predictions[record.id]["EukRep"] = 'prokaryote'
    # Checking Tiara, Whokaryote and Whokaryote+Tiara predictions
    predictions[record.id]["Tiara"] = sp.check_output([f'grep -P "{record.id+"\\"+"s"}" output/{fastaName}/tiara/tiara.txt | cut -f2'], shell=True).decode("UTF-8").replace('\n', "")
    predictions[record.id]["Whokaryote"] = sp.check_output([f'grep -P "{record.id+"\\"+"s"}" output/{fastaName}/whokaryoteS/whokaryote_predictions_S.tsv | cut -f2'], shell=True).decode("UTF-8").replace('\n', "")
    predictions[record.id]["Whokaryote+Tiara"] = sp.check_output([f'grep -P "{record.id+"\\"+"s"}" output/{fastaName}/whokaryoteT/whokaryote_predictions_T.tsv | cut -f2'], shell=True).decode("UTF-8").replace('\n', "")
    # Checking DMC predictions
    resultDMC = list(map(float, sp.check_output([f'grep -P "{record.id+"\\"+"s"}" output/{fastaName}/dmc/{fastaName}.fa_pred_one-hot_hybrid.tsv'], shell=True).decode("UTF-8").replace('\n', "").split("\t")[1:]))
    predictions[record.id]["DeepMicroClass"] = ["Eukaryote", "EukaryoteVirus", "Plasmid", "Prokaryote", "ProkaryoteVirus"][resultDMC.index(max(resultDMC))]

# Creation of a DataFrame with the predictions
flatPredictions = {}
for i in range(len(predictions)):
    flatPredictions[list(predictions.keys())[i]] = list(predictions[list(predictions.keys())[i]].values())
dataFrame = pd.DataFrame.from_dict(flatPredictions, orient="index", columns=list(predictions[record.id].keys()))

# Uniformize the predictions
dataFrame.replace("", "not classified", inplace=True)
if precisionLevel == 1 :
    dataFrame.replace("EukaryoteVirus", "virus", inplace=True)
    dataFrame.replace("ProkaryoteVirus", "virus", inplace=True)
    dataFrame.replace("Prokaryota ", "prokaryote", inplace=True)
    dataFrame.replace("Prokaryote", "prokaryote", inplace=True)
    dataFrame.replace("Bacteria ", "prokaryote", inplace=True)
    dataFrame.replace("bacteria", "prokaryote", inplace=True)
    dataFrame.replace("prokarya", "prokaryote", inplace=True)
    dataFrame.replace("archaea", "prokaryote", inplace=True)
    dataFrame.replace("eukarya", "eukaryote", inplace=True)
    dataFrame.replace("Eukaryota ", "eukaryote", inplace=True)
    dataFrame.replace("Eukaryote", "eukaryote", inplace=True)
    dataFrame.replace("Viruses ", "virus", inplace=True)
    dataFrame.replace('Plasmid', 'prokaryote', inplace=True)
elif precisionLevel == 2 :
    pass

dataFrame.to_csv(f"output/{fastaName}/predictions.tsv", sep="\t")
print(f"Predictions of tools output saved at output/{fastaName}/predictions.tsv")

# If a real taxonomy is given, compare the result of the tool with it
if benchParam != None :
    bench = tc.check_truth(f"output/{fastaName}/predictions.tsv", truthPath, benchQuery)
    print("\nCounting table of the predicted and real taxonomy")
    print(bench[0])
    print("\nCounting table of the veracity of classification per tool")
    print(bench[1])
    bench[2].to_csv(f"output/{fastaName}/performances.tsv", sep="\t")
    print(f"\nPredictive performance per tool saved at output/{fastaName}/performances.tsv")
    print(bench[2])
    print("\n")
    print("Scores per tool :")
    for col in bench[2] :
        accuracy = sc.get_accuracy(bench[2][col]['TP'],bench[2][col]['TN'],bench[2][col]['FP'],bench[2][col]['FN'])
        fscore = sc.get_f1score(bench[2][col]['TP'],bench[2][col]['TN'],bench[2][col]['FP'],bench[2][col]['FN'])
        mcc = sc.get_mcc(bench[2][col]['TP'],bench[2][col]['TN'],bench[2][col]['FP'],bench[2][col]['FN'])
        print(f"{col} | Accuracy : {accuracy} | F1-score : {fscore} | MCC : {mcc}")
    print("\n")