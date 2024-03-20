# -*- coding: utf-8 -*-
import os, sys
import random
import matplotlib
import matplotlib.pyplot as plt
from Bio import SeqIO

def split_multifasta(fastaPath, fastaName, outdir, nSeq=None, minLength=None, maxLength=None, knownIds=None, proportions=None) :
    """
    This function generates a multifasta from another multifasta containing only
    nSeq sequences of a specified length, as the results of the tools used by WizardEye
    can be length-dependant. It also allow to select only known contig from a list.
    """
    print("Selecting the sequences from the fasta...")
    # Check if proportions are important
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
        if proportions == None and isinstance(knownIds, dict) is True :
            knownIds = [ids for kingdoms in knownIds.values() for ids in kingdoms]
        elif proportions != None and isinstance(knownIds, dict) is True :
            nSeqPerKingdom = [round(nSeq*(float(value)/100)) for value in proportions.values()]
        # To avoid a bias in the selection, parse semi-randomly the records except if we want all the seq
        for i in range(len(records)) :
            if nSeq == float('inf') :
                record = records[i]
            else :
                record = records[random.randint(0, len(records)-1)]
            if len(selectedRecs) < nSeq:
                if len(record) >= minLength and len(record) <= maxLength :
                    id = record.id.split(" ")[0]
                    # To avoid the selection of an already random selected contig
                    if id not in selectedIds :
                        if isinstance(knownIds, list) is True or isinstance(knownIds, dict) is True :
                            if proportions == None :
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
            else:
                break
        SeqIO.write(selectedRecs, newFasta, "fasta")
        if maxLength == float("inf") :
            print(f"{len(selectedRecs)} sequences larger or making {minLength} pb selected")
        else :
            print(f"{len(selectedRecs)} sequences between {minLength} and {maxLength} pb selected")
    else :
        print("Fasta found...")
    return newFasta