# -*- coding: utf-8 -*-
import os, sys, configparser, argparse, shutil, time, resource
import pandas as pd
import subprocess as sp
from Bio import SeqIO
from modules import fastatools as ft
from modules import truthchecker as tc
from modules import scores as sc

# Create a class allowing to deactivate the print in the standard output
class HiddenPrints :
    def __enter__(self):
        self._original_stdout = sys.stdout
        sys.stdout = open(os.devnull, "w")
    def __exit__(self, exc_type, exc_val, exc_tb):
        sys.stdout.close()
        sys.stdout = self._original_stdout

def main(fastaPath, ncontigs = None, minLength = None, maxLength = None, precisionLevel = 1, debugMode = None, truthFormat = None, positiveRef = None, realKingdom = None, kingdomProp = None, kingdomList = None, exName = None, clean=None, cat=None, config=None) :
    """
    This function is the main function of WizardEye. It allows to run a specific
    analysis using different tools to predict the kingdom of contigs contained in
    a multifasta file. Results are saved in the ./output/.../predictions.tsv, with
    some metrics if WizardEye is launched in benchmark mode.
    """

    # Print the tools output in the stdout or not according to the debug mode
    if debugMode == "y":
        toolsOutput = ""
    else :
        toolsOutput = "> /dev/null"

    # If the tool is launched with the expert() function, use the benchmark name
    if exName == None :
        exName = ""
    else :
        exName = "EX" + exName + "_"
    
    # Format the final name of the fasta
    fastaName = f"{exName}{os.path.basename(fastaPath).split(".")[0]}_{ncontigs}seq_{minLength}pb"
    if maxLength != None :
        fastaName = fastaName[:-2] + f"to{maxLength}pb"
    if kingdomProp != None and kingdomList != None :
        fastaName += "_"
        for i in range(len(kingdomProp.split(":"))) :
            fastaName += kingdomProp.split(":")[i] + kingdomList.split(":")[i][:4] + "-"
        fastaName = fastaName[:-1]
    
    if os.path.exists(f"output/{fastaName}") is False :
        os.mkdir(f"output/{fastaName}")
    
    # If the truth is already known, get the contig ids where the taxonomy is known to only use these ids in the newFasta
    if truthFormat == "AK":
        print("\nReading real taxonomies...")
        truthPath, knownIds = realKingdom, {}
        truth = open(realKingdom)
        for line in truth.readlines()[1:] :
            if line.split("\t")[1].replace("\n", "") != "unknown":
                if line.split("\t")[1] not in knownIds.keys() :
                    knownIds[line.split("\t")[1]] = []
                knownIds[line.split("\t")[1]].append(line.split("\t")[0])
        truth.close()
        if kingdomProp != None and kingdomList != None :
            proportions = {}
            for i in range(len(kingdomList.split(":"))) :
                proportions[kingdomList.split(":")[i]] = kingdomProp.split(":")[i]
            if sum([float(value) for value in proportions.values()]) != float(100):
                print("The sum of the proportions should be 100% !")
                sys.exit()
        else :
            proportions = None
        if exName == None :
            shutil.copy(truthPath, f"output/{fastaName}")
    else :
        knownIds, proportions = None,  None
    
    # Make a copy of the multifasta with the selectionned size, n contigs
    newFasta, nFindContig  = ft.split_multifasta(fastaPath, fastaName, "input/", nSeq=ncontigs, minLength=minLength,  maxLength=maxLength, knownIds=knownIds, proportions=proportions, exName=exName)
    
    # Format the real taxonomy for non-already existing truth.tsv
    if truthFormat == "CAMI":
            print("\nReading real taxonomies...")
            mail = config['user']['mail']
            truthPath = tc.get_cami_tax(newFasta, realKingdom, f"output/{fastaName}/", mail)
    if truthFormat == "OK":
        print("\nReading real taxonomies...")
        truthPath = tc.get_oneKingdom_tax(newFasta, realKingdom, f"output/{fastaName}/")
    
    # Run tools and measure the elapsedTimes per tool, minus the time to activate the environment
    elapsedTimes, elapsedProcTimes= [], []
    
    print("\nRunning tools...")
    if cat == "y" :
        if os.path.exists(f"output/{fastaName}/CAT") is False :
            os.mkdir(f"output/{fastaName}/CAT")
            envStartTime, envStartProcTime = time.time(), time.process_time()
            os.system(f'singularity run {config["CAT"]["CAT_sif"]} CAT')
            envActivationTime, envProcActivationTime = time.time() - envStartTime, time.process_time() - envStartProcTime
            print("Running CAT...")
            startTime, startProcTime = time.time(), time.process_time()
            os.system(f'singularity run {config["CAT"]["CAT_sif"]} CAT contigs -c {newFasta} -d {config["CAT"]["CAT_DB"]} --path_to_diamond {config["CAT"]["CAT_Diamond_DB"]} -t {config["CAT"]["CAT_Taxonomy_DB"]} --out_prefix "./output/{fastaName}/CAT/out.CAT" --block_size 4 --index_chunks 2 {toolsOutput}')
            os.system(f'singularity run {config["CAT"]["CAT_sif"]} CAT add_fastaNames -i "./output/{fastaName}/CAT/out.CAT.contig2classification.txt" -o "./output/{fastaName}/CAT/out.CAT.lineage" -t {config["CAT"]["CAT_Taxonomy_DB"]} {toolsOutput}')
            elapsedTimes.append(time.time() - startTime - envActivationTime)
            elapsedProcTimes.append(time.process_time() - startProcTime - envProcActivationTime)
        else :
            print("CAT predictions found...")
    else :
        print("CAT will not be launched.")
    
    if os.path.exists(f"output/{fastaName}/eukrep") is False :
        os.mkdir(f"output/{fastaName}/eukrep")
        envStartTime, envStartProcTime = time.time(), time.process_time()
        os.system(f'eval "$(conda shell.bash hook)" && conda activate {config["environment"]["EukRep"]}')
        envActivationTime, envProcActivationTime = time.time() - envStartTime, time.process_time() - envStartProcTime
        print("Running EukRep...")
        startTime, startProcTime = time.time(), time.process_time()
        os.system(f'eval "$(conda shell.bash hook)" && conda activate {config["environment"]["EukRep"]} && EukRep -i {newFasta} -o output/{fastaName}/eukrep/eukrep.txt {toolsOutput}')
        elapsedTimes.append(time.time() - startTime - envActivationTime)
        elapsedProcTimes.append(time.process_time() - startProcTime - envProcActivationTime)
    else :
        print("EukRep predictions found...")
    
    if os.path.exists(f"output/{fastaName}/tiara") is False :
        os.mkdir(f"output/{fastaName}/tiara")
        envStartTime, envStartProcTime = time.time(), time.process_time()
        os.system(f'eval "$(conda shell.bash hook)" && conda activate {config["environment"]["Tiara"]}')
        envActivationTime, envProcActivationTime = time.time() - envStartTime, time.process_time() - envStartProcTime
        print("Running Tiara...")
        startTime, startProcTime = time.time(), time.process_time()
        os.system(f'eval "$(conda shell.bash hook)" && conda activate {config["environment"]["Tiara"]} && tiara -i {newFasta} -o output/{fastaName}/tiara/tiara.txt {toolsOutput}')
        elapsedTimes.append(time.time() - startTime - envActivationTime)
        elapsedProcTimes.append(time.process_time() - startProcTime - envProcActivationTime)
    else :
        print("Tiara predictions found...")
    
    if os.path.exists(f"output/{fastaName}/whokaryoteS") is False :
        envStartTime, envStartProcTime = time.time(), time.process_time()
        os.system(f'eval "$(conda shell.bash hook)" && conda activate {config["environment"]["Whokaryote"]}')
        envActivationTime, envProcActivationTime = time.time() - envStartTime, time.process_time() - envStartProcTime
        print("Running Whokaryote...")
        startTime, startProcTime = time.time(), time.process_time()
        os.system(f'eval "$(conda shell.bash hook)" && conda activate {config["environment"]["Whokaryote"]} && whokaryote.py --contigs {newFasta} --outdir output/{fastaName}/whokaryoteS --model S --minsize {minLength} {toolsOutput}')
        elapsedTimes.append(time.time() - startTime - envActivationTime)
        elapsedProcTimes.append(time.process_time() - startProcTime - envProcActivationTime)
    else :
        print("Whokaryote predictions found...")
    
    if os.path.exists(f"output/{fastaName}/whokaryoteT") is False :
        envStartTime, envStartProcTime = time.time(), time.process_time()
        os.system(f'eval "$(conda shell.bash hook)" && conda activate {config["environment"]["Whokaryote"]}')
        envActivationTime, envProcActivationTime = time.time() - envStartTime, time.process_time() - envStartProcTime
        print("Running Whokaryote+Tiara...")
        startTime, startProcTime = time.time(), time.process_time()
        os.system(f'eval "$(conda shell.bash hook)" && conda activate {config["environment"]["Whokaryote"]} && whokaryote.py --contigs {newFasta} --outdir output/{fastaName}/whokaryoteT --model T --minsize {minLength} {toolsOutput}')
        elapsedTimes.append(time.time() - startTime - envActivationTime)
        elapsedProcTimes.append(time.process_time() - startProcTime - envProcActivationTime)
    else :
        print("Whokaryote+Tiara predictions found...")
    
    if os.path.exists(f"output/{fastaName}/dmc") is False :
        envStartTime, envStartProcTime = time.time(), time.process_time()
        os.system(f'eval "$(conda shell.bash hook)" && conda activate {config["environment"]["DeepMicroClass"]}')
        envActivationTime, envProcActivationTime = time.time() - envStartTime, time.process_time() - envStartProcTime
        print("Running DeepMicroClass...")
        startTime, startProcTime = time.time(), time.process_time()
        os.system(f'eval "$(conda shell.bash hook)" && conda activate {config["environment"]["DeepMicroClass"]} && DeepMicroClass predict -i {newFasta} -o output/{fastaName}/dmc/ {toolsOutput}')
        elapsedTimes.append(time.time() - startTime - envActivationTime)
        elapsedProcTimes.append(time.process_time() - startProcTime - envProcActivationTime)
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
        if cat == 'y' :
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
    if realKingdom != None :
        bench = tc.check_truth(f"output/{fastaName}/predictions.tsv", truthPath, positiveRef)
        print("\nCounting table of the predicted and real taxonomy")
        print(bench[0])
        print("\nCounting table of the veracity of classification per tool")
        print(bench[1])
        bench[2].to_csv(f"output/{fastaName}/predictions_metrics_{positiveRef}.tsv", sep="\t")
        print(f"\nPredictive performance per tool saved at output/{fastaName}/predictions_metrics_{positiveRef}.tsv")
        print(bench[2])
        print("\nScores per tool :")
        for col in bench[2] :
            accuracy = sc.get_accuracy(bench[2][col]['TP'],bench[2][col]['TN'],bench[2][col]['FP'],bench[2][col]['FN'])
            fscore = sc.get_f1score(bench[2][col]['TP'],bench[2][col]['TN'],bench[2][col]['FP'],bench[2][col]['FN'])
            mcc = sc.get_mcc(bench[2][col]['TP'],bench[2][col]['TN'],bench[2][col]['FP'],bench[2][col]['FN'])
            print(f"{col} | Accuracy : {accuracy} | F1-score : {fscore} | MCC : {mcc}")
    
        if len(elapsedTimes) == len(bench[2].columns) :
            print("\nTimes per tool :")
            times = {"Execution time": {}, "CPU Execution time": {}}
            for i in range(len(list(bench[2].columns))) :
                times["Execution time"][bench[2].columns[i]] = elapsedTimes[i]
                times["CPU Execution time"][bench[2].columns[i]] = elapsedProcTimes[i]
            dfTimes = pd.DataFrame.from_dict(times, orient="index")
            print(dfTimes)
            dfTimes.to_csv(f"output/{fastaName}/performances.tsv", sep="\t")
            print("\nNote : the conda environment activation is substracted from the tool's execution times, which can sometimes result in negative CPU time.\n")
    
        # Print the used parameters in a file
        resume = open(f"output/{fastaName}/params.tsv", "w")
        resume.write(f"inputPath\t{fastaPath}\ntruthPath\t{truthPath}\nnWantedcontigs\t{ncontigs}\nnFoundContigs\t{nFindContig}\nminimum\t{minLength}\nmaximum\t{maxLength}\nproportions\t{kingdomProp}\nof\t{kingdomList}")
        resume.close()
    
    # Clean the output directory
    if clean == "y" :
        print("Cleaning output directory...")
        if cat == "y" :
            shutil.rmtree(f"output/{fastaName}/CAT/")
        shutil.rmtree(f"output/{fastaName}/eukrep/")
        shutil.rmtree(f"output/{fastaName}/tiara/")
        shutil.rmtree(f"output/{fastaName}/whokaryoteS/")
        shutil.rmtree(f"output/{fastaName}/whokaryoteT/")
        shutil.rmtree(f"output/{fastaName}/dmc/")
    
    # If it's a serie of benchmark, move all the output in the serie's directory
    if exName != "" :
        exName = exName[:-1]
        if os.path.exists(f"output/{exName}") is False :
            os.mkdir(f"output/{exName}")
        for dir in os.listdir("output") :
            if exName in dir and exName != dir :
                shutil.move(f"output/{dir}", f"output/{exName}/{dir}")

def expert(mode, fastaPath, name, ncontigs = None, minLength=None, maxLength=None, positiveRef = None, debugMode = None, realKingdom = None, config = None, kingdomProp=None, kingdomList=None) :
    """
    This function allows to run WizardEye several times, modifying the parameters at each
    iteration for benchmarking purposes. Two benchmarking modes are currently available :
    - length, to launch WizardEye on contigs of progressive length
    - proportion, to launch WizardEye on contigs with different proportions of known kingdoms
    """
    # Length mode
    if mode == "length":
        lengthRange = args.lengthrange
        # Set the maximum length and the number of repetition
        if maxLength == None:
            rep = max([len(record) for record in SeqIO.parse(fastaPath, "fasta")]) / lengthRange
        else:
            rep = (maxLength - minLength) / lengthRange

        print(f"WizardEye will be launched {int(rep)} times")

        # Create the ranges
        minLatRep = list(range(minLength, int(minLength + rep * lengthRange), lengthRange))
        maxLatRep = list(range(minLength + lengthRange, int(rep * lengthRange + lengthRange + minLength), lengthRange))

        # Running the tools thanks to WizardEye
        for i in range(int(rep)):
            if os.path.exists(
                    f"output/EX{name}/EX{name}_{os.path.basename(fastaPath).split(".")[0]}_{ncontigs}seq_{minLatRep[i]}to{maxLatRep[i]}pb") is False:
                print(f"[{i + 1}/{int(rep)}] Launching WizardEye for contigs from {minLatRep[i]} to {maxLatRep[i]} pb")
                if debugMode != "y" :
                    with HiddenPrints() :
                        main(fastaPath=fastaPath, ncontigs=ncontigs, minLength=minLatRep[i], maxLength=maxLatRep[i], exName=name, realKingdom=realKingdom, truthFormat="AK", positiveRef=positiveRef, clean="y", config=config, kingdomProp=kingdomProp, kingdomList=kingdomList)
                else :
                    main(fastaPath=fastaPath, ncontigs=ncontigs, minLength=minLatRep[i], maxLength=maxLatRep[i], exName=name, realKingdom=realKingdom, truthFormat="AK", positiveRef=positiveRef, clean="y", config=config, kingdomProp=kingdomProp, kingdomList=kingdomList)
            else:
                print( f"WizardEye was already launched on {fastaPath} for {ncontigs} contigs between {minLatRep[i]} and {maxLatRep[i]} pb")
    # Proportion mode
    elif mode == "proportion":
        propK = kingdomList
        prop = kingdomProp
        if prop == None:
            print("The Taxonomic mode need proportions.")
            sys.exit()
        else:
            dicTax = {}
            tax = propK.split(":")
            for i in range(len(tax)):
                dicTax[tax[i]] = prop.split(":")[i].split(",")
            rep = len(dicTax[tax[0]])

        print(f"WizardEye will be launched {int(rep)} times")
        for i in range(int(rep)):
            iprop = ":".join([dicTax[tax[ntax]][i] for ntax in range(len(tax))])
            fastaName = f"EX{name}_{os.path.basename(fastaPath).split(".")[0]}_{ncontigs}seq_{minLength}pb"
            if maxLength != None:
                fastaName = fastaName[:-2] + f"to{maxLength}pb"
            fastaName += "_"
            for t in range(len(iprop.split(":"))):
                fastaName += iprop.split(":")[t] + propK.split(":")[t][:4] + "-"
            if os.path.exists(f"output/EX{name}/{fastaName[:-1]}") is False:
                print(
                    f"[{i + 1}/{int(rep)}] Launching WizardEye for proportions of {"% & ".join([dicTax[tax[ntax]][i] for ntax in range(len(tax))])}% (respectively for {" & ".join(tax)})")
                if debugMode != "y":
                    with HiddenPrints():
                        main(fastaPath=fastaPath, ncontigs=ncontigs, minLength=minLength, maxLength=maxLength,
                             exName=name, realKingdom=realKingdom, truthFormat="AK", positiveRef=positiveRef, clean="y",
                             kingdomProp=iprop, kingdomList=propK, config=config)
                else:
                    main(fastaPath=fastaPath, ncontigs=ncontigs, minLength=minLength, maxLength=maxLength,
                         exName=name, realKingdom=realKingdom, truthFormat="AK", positiveRef=positiveRef, clean="y",
                         kingdomProp=iprop, kingdomList=propK, config=config)
            else:
                print(
                    f"WizardEye already launched on {fastaPath} on {ncontigs} contigs from proportions of {"% & ".join([dicTax[tax[ntax]][i] for ntax in range(len(tax))])}% (respectively for {" & ".join(tax)})")

    # Save the metrics
    f1score, accuracy, mcc, realTimes, cpuTimes, dfParams = {}, {}, {}, {}, {}, pd.DataFrame()
    for dir in os.listdir(f"output/EX{name}/"):
        label = dir.split("_")[-1]
        if os.path.exists(f"output/EX{name}/{dir}/predictions_metrics_{positiveRef}.tsv") is True:
            df = pd.read_csv(f"output/EX{name}/{dir}/predictions_metrics_{positiveRef}.tsv", sep="\t", index_col=0,
                             header=0)
            dfPerf = pd.read_csv(f"output/EX{name}/{dir}/performances.tsv", sep="\t", index_col=0, header=0)
            params = pd.read_csv(f"output/EX{name}/{dir}/params.tsv", sep="\t", index_col=0, header=None)
            params.loc["ID"] = label
            dfParams = pd.concat([dfParams, params], axis=1)
            for col in df.columns:
                if len(f1score.keys()) != len(df.columns):
                    f1score[col] = {}
                    accuracy[col] = {}
                    mcc[col] = {}
                    realTimes[col] = {}
                    cpuTimes[col] = {}
                f1score[col][f"{label}"] = sc.get_f1score(df[col]['TP'], df[col]['TN'], df[col]['FP'],
                                                              df[col]['FN'])
                accuracy[col][f"{label}"] = sc.get_accuracy(df[col]['TP'], df[col]['TN'], df[col]['FP'],
                                                                df[col]['FN'])
                mcc[col][f"{label}"] = sc.get_mcc(df[col]['TP'], df[col]['TN'], df[col]['FP'], df[col]['FN'])
                realTimes[col][f"{label}"] = dfPerf[col]['Execution time']
                cpuTimes[col][f"{label}"] = dfPerf[col]['CPU Execution time']

    # Save the DataFrame resuming the score per range
    dfFscore = pd.DataFrame.from_dict(f1score)
    dfAccuracy = pd.DataFrame.from_dict(accuracy)
    dfMcc = pd.DataFrame.from_dict(mcc)
    dfrealTimes = pd.DataFrame.from_dict(realTimes)
    dfcpuTimes = pd.DataFrame.from_dict(cpuTimes)
    dfParams.transpose().to_csv(f"output/EX{name}/params.tsv", sep="\t")
    dfFscore.to_csv(f"output/EX{name}/f1score.tsv", sep="\t")
    dfAccuracy.to_csv(f"output/EX{name}/accuracy.tsv", sep="\t")
    dfMcc.to_csv(f"output/EX{name}/mcc.tsv", sep="\t")
    dfrealTimes.to_csv(f"output/EX{name}/realTimes.tsv", sep="\t")
    dfcpuTimes.to_csv(f"output/EX{name}/cpuTimes.tsv", sep="\t")

    # Copy the complete truth file in the directory
    shutil.copy(realKingdom, f"output/EX{name}/truth.tsv")

if __name__ == "__main__" :
    # Parse arguments
    parser = argparse.ArgumentParser(description="WizardEye summarize kingdom identification from contigs thanks to different tools specialised in the detection of eukaryotic contigs",
                                     epilog="More help on the GitHub repository, at www.github.com/TheLokj/WizardEye")
    # Mandatory arguments
    parser.add_argument('mode', choices=["basic", "expert"], help = 'the number of contigs to select')
    parser.add_argument('-i', '--input', metavar='path', type=str, required=True, help="the input multifasta containing the contigs")
    # Contigs options arguments
    parser.add_argument('-n', '--ncontigs', metavar="n", type=int, help='the number of contigs to select')
    parser.add_argument('-mnl', '--minlength', metavar="n", type=int, default=1, help='the minimun length of the analyzed contigs in pb')
    parser.add_argument('-mxl', '--maxlength', metavar="n", type=int,  help='the maximum length of the analyzed contigs in pb')
    # General arguments
    parser.add_argument('-c', '--cat', metavar='y/n', type=str, choices=["y", "n"], help="use or not CAT (pretty long step)")
    parser.add_argument('-l', '--level', metavar="1/2", type=int, choices=[1, 2], default=1, help='the level of precision of the summary (1 : uniformized : prokaryote, eukaryote, virus ; 2 : raw output of tools)')
    parser.add_argument('-d', '--debug', metavar="y/n", type=str, choices=["y", "n"], help='the debug mode, print everything in the stdout')
    parser.add_argument('-cl', '--clean', metavar="y/n", choices=["y", "n"], type=str,help='remove or not the tools output after resume the predictions')
    # Benchmarking specific options
    parser.add_argument('-t', '--truth', metavar="...", type=str, help='(benchmarking) the path containing the real taxonomy (or the kingdom in One-Kingdom mode)')
    parser.add_argument('-tf', '--truthformat', metavar="format", type=str, choices=["AK", "OK", "CAMI"],help='(benchmarking) the truth format : AK (Already-Known), OK (One-Kingdom), CAMI')
    parser.add_argument('-pr', '--positiveref', metavar="kingdom", type=str,help='(benchmarking) the kingdom used as positive for calculating scores')
    parser.add_argument('-kp', '--kingdomprop', metavar="X:Y:Z", type=str, help='(benchmarking, basic mode only) the proportion of each kingdom wanted in the final analysed multifasta (in %%)')
    parser.add_argument('-kl', '--kingdomlist', metavar="kingdom1:kingdom2:kingdom3", type=str, help='(benchmarking, basic mode only) the name of the kingdoms in -kp proportions')
    # Expert mode specific options
    parser.add_argument('-na', '--name', metavar="name", type=str, help='(expert) the benchmark name')
    parser.add_argument('-mo', '--bmode', metavar="mode", choices=["length", "proportion"], type=str, help='(expert) the benchmark mode')
    parser.add_argument('-lr', '--lengthrange', metavar="n", type=int,help="(expert, length mode) The size of the range to increament at each launch")
    parser.add_argument('-ekp', '--exkingdomprop', metavar="n", type=str, help="(expert, proportion mode) the different proportion per repetition (20,80) per taxa (20,80:80,20)  ")
    parser.add_argument('-ekl', '--exkingdomlist', metavar="n", type=str, help="(expert, proportion mode) the taxa separated with : (eukaryote:prokaryote)")
    args = parser.parse_args()

    # Print the name of the tool and its version
    print("\n▄█     █▄   ▄█   ▄███████▄     ▄████████    ▄████████ ████████▄     ▄████████ ▄██   ▄      ▄████████\n███     ███ ███  ██▀     ▄██   ███    ███   ███    ███ ███   ▀███   ███    ███ ███   ██▄   ███    ███ \n███     ███ ███▌       ▄███▀   ███    ███   ███    ███ ███    ███   ███    █▀  ███▄▄▄███   ███    █▀ \n███     ███ ███▌  ▀█▀▄███▀▄▄   ███    ███  ▄███▄▄▄▄██▀ ███    ███  ▄███▄▄▄     ▀▀▀▀▀▀███  ▄███▄▄▄     \n███     ███ ███▌   ▄███▀   ▀ ▀███████████ ▀▀███▀▀▀▀▀   ███    ███ ▀▀███▀▀▀     ▄██   ███ ▀▀███▀▀▀     \n███     ███ ███  ▄███▀         ███    ███ ▀███████████ ███    ███   ███    █▄  ███   ███   ███    █▄  \n███ ▄█▄ ███ ███  ███▄     ▄█   ███    ███   ███    ███ ███   ▄███   ███    ███ ███   ███   ███    ███ \n▀███▀███▀  █▀    ▀████████▀   ███    █▀    ███    ███ ████████▀    ██████████  ▀█████▀    ██████████\n                                                                                                 v0.0.2 by TheLokj")
    print("\nMore information and updates on www.github.com/TheLokj/WizardEye/\n")

    # Read config file
    if os.path.exists("config.ini") is False:
        print("Error : the configuration file config.ini does not exist.")
        sys.exit()
    else:
        config = configparser.ConfigParser()
        config.read("config.ini")

    # Check if the input and output directories exist
    if os.path.exists("input") is False:
        os.mkdir("input")
    if os.path.exists("output") is False:
        os.mkdir("output")

    if args.mode == "basic" :
        # Launch WizardEye in the basic mode, just the main function with specific parameters
        main(fastaPath=args.input, ncontigs=args.ncontigs,
             minLength=args.minlength, maxLength=args.maxlength,
             cat=args.cat, precisionLevel=args.level,
             truthFormat=args.truthformat, realKingdom=args.truth,
             kingdomProp=args.kingdomprop, kingdomList=args.kingdomlist,
             positiveRef=args.positiveref, exName=args.name,
             debugMode=args.debug, clean=args.clean, config=config)
    elif args.mode == "expert" :
        # Launch the benchmarking function that launch ultiple times the main()
        expert(mode=args.bmode, fastaPath=args.input,
               name=args.name, ncontigs = args.ncontigs,
               minLength=args.minlength, maxLength=args.maxlength,
               positiveRef = args.positiveref, realKingdom = args.truth,
               kingdomProp=args.exkingdomprop, kingdomList=args.exkingdomlist,
               debugMode = args.debug, config=config)