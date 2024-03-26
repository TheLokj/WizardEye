# -*- coding: utf-8 -*-
import os, sys
import random
import matplotlib
import matplotlib.pyplot as plt
from Bio import SeqIO
import numpy as np

def split_multifasta(fastaPath, fastaName, outdir, nSeq=None, minLength=None, maxLength=None, knownIds=None, kingdomProp=None, exName=None) :
    """
    This function generates a multifasta from another multifasta containing only
    nSeq sequences of a specified length, as the results of the tools used by WizardEye
    can be length-dependant. It also allow to select only known contig from a list.
    """
    if exName == "" :
        print("Plotting the sequences length histogram...")
        hist_length(fastaPath, fastaName)
    print("Selecting the sequences from the fasta...")
    if minLength == None :
        minLength = 0
    if maxLength == None :
        maxLength = float('inf')
    newFasta, selectedRecs, selectedIds = (f"input/{fastaName}.fa"), [], []
    # To avoid unnecessary step, check if the fasta currently exist
    if os.path.exists(newFasta) is False:
        records = []
        for record in SeqIO.parse(fastaPath, "fasta"):
            records.append(record)
        if nSeq == None :
            nSeq = len(records)
        if kingdomProp == None and isinstance(knownIds, dict) is True :
            knownIds = [ids for kingdoms in knownIds.values() for ids in kingdoms]
        elif kingdomProp != None and isinstance(knownIds, dict) is True :
            nSeqPerKingdom = [round(nSeq*(float(value)/100)) for value in kingdomProp.values()]
        # To avoid a bias in the selection, parse semi-randomly the records except if we want all the seq
        for i in range(len(records)) :
            if nSeq == float('inf') :
                record = records[i]
            else :
                rindex = random.randint(0, len(records)-1)
                record = records[rindex]
            if len(selectedRecs) < nSeq:
                if len(record) >= minLength and len(record) <= maxLength :
                    id = record.id.split(" ")[0]
                    if isinstance(knownIds, list) is True or isinstance(knownIds, dict) is True :
                        if kingdomProp == None :
                            if id in knownIds:
                                selectedRecs.append(record)
                                selectedIds.append(id)
                        else :
                            for iKingdom in range(len(knownIds.keys())) :
                                if id in knownIds[list(knownIds.keys())[iKingdom]] :
                                    if nSeqPerKingdom[iKingdom] != 0 :
                                            selectedRecs.append(record)
                                            selectedIds.append(id)
                                            nSeqPerKingdom[iKingdom] -= 1
                    else :
                        selectedRecs.append(record)
                        selectedIds.append(id)
                records.pop(rindex)
            else:
                break
            # To avoid an analyzis on empty data or partial data, stop the tool if nSeq isn't reach
        if len(selectedRecs) < nSeq :
            print(f"Warning : There are only {len(selectedRecs)} compatible sequences (<{nSeq}) between {minLength} and {maxLength} in {fastaName}")
            if isinstance(knownIds, list) is True or isinstance(knownIds, dict) is True:
                print("This may affect the requested proportions.")
        if len(selectedRecs) == 0 and exName == None:
            print("There isn't contig of this length in the input fasta")
            sys.exit()
        else :
            SeqIO.write(selectedRecs, newFasta, "fasta")
        if maxLength == float("inf") :
            print(f"{len(selectedRecs)} sequences larger or making {minLength} pb selected")
        else :
            print(f"{len(selectedRecs)} sequences between {minLength} and {maxLength} pb selected")
    else :
        print("Fasta found...")
    return newFasta, len(selectedRecs)

def hist_length(fastaPath, fastaName) :
    length = [len(record) for record in SeqIO.parse(fastaPath, "fasta")]
    plt.hist(length, log=True)
    plt.xlabel("Length (pb)")
    plt.ylabel("logCount")
    plt.title(f"Histogram of sequences length from \n{fastaPath.split("/")[-1]}")
    plt.savefig(f"output/{fastaName}/original_length_distribution.png")