# -*- coding: utf-8 -*-
import sys, os, math
import pandas as pd
from .fastatools import *
from Bio import Entrez
from Bio import SeqIO

def get_cami_tax(fastaPath, truthPath, outdir, mail) :
    """
    This function is used in the CAMI mode. It generates a truth.tsv from the
    gsa_mapping.tsv provided by a dataset generated thanks to CAMISIM containing
    for each record id the associated kingdom.
    """
    Entrez.email = mail

    # Read the fasta
    seqList, taxIds, seq = [], [], {}
    for record in SeqIO.parse(fastaPath, "fasta") :
        seqList.append(record)

    # Check the real taxonomy
    print("Checking for CAMI taxonomies...")
    dfTruth = pd.read_csv(truthPath, sep='\t', index_col=0)
    for i in range(len(seqList)):
        taxIds.append(dfTruth.loc[seqList[i].id]['tax_id'])

    # Obtain the division name
    # As the NCBI as limited to 10000 ids, do multiple request if it's necessary
    request = []
    for i in range(math.ceil(len(seqList)/9999)):
        print(f"Requesting NCBI... ({i+1}/{math.ceil(len(seqList)/9999)})")
        request.extend(Entrez.read(Entrez.efetch(id=taxIds[i*9999:((i+1)*9999)], db='taxonomy', rettype='json')))

    for i in range(len(seqList)):
        seq[seqList[i].id] = {"tax_id": dfTruth.loc[seqList[i].id]['tax_id'], "division": request[i]['Division']}

    # Export the result
    dfTax = pd.DataFrame.from_dict(seq, orient="index")
    dfTax.replace('Bacteria', 'prokaryote', inplace=True)

    dfTax.to_csv(f"{outdir}truth.tsv", sep="\t")
    print(f"Real taxonomy from CAMI saved at {outdir}truth.tsv")

    return f"{outdir}truth.tsv"

def get_oneKingdom_tax(fastaPath, kingdom, outdir, mode="standalone") :
    """
    This function is used in the one-kingdom mode. It generates a truth.tsv containing
    for each record id the provided kingdom.
    """
    # Read the fasta
    seq = {}
    if mode =="standalone" :
        for record in SeqIO.parse(fastaPath, "fasta") :
            seq[record.id] = {"division": kingdom}
    elif mode=="GCA" :
        for record in SeqIO.parse(fastaPath, "fasta") :
            seq[record.description.split('contig: ')[1].split(",")[0]] = {"division": kingdom}

    dfTax = pd.DataFrame.from_dict(seq, orient="index")

    if mode=="standalone" :
        # Export the result
        dfTax.to_csv(f"{outdir}truth.tsv", sep="\t")
        print(f"One-kingdom taxonomy saved at {outdir}truth.tsv")
        return f"{outdir}truth.tsv"
    elif mode=="GCA" :
        return dfTax

def get_gca_tax(directory,  outdir) :
    """
    This function is no longer used
    """
    dfTax = pd.DataFrame()
    for subdir in os.listdir(directory) :
        print(f"Checking the sequences from {subdir}")
        for file in os.listdir(f"{directory}/{subdir}") :
            dfTax = pd.concat([dfTax, get_oneKingdom_tax(f"{directory}/{subdir}/{file}", subdir, None, mode="GCA")])
    dfTax.to_csv(f"{outdir}truth.tsv", sep="\t")
    print(f"GCA known taxonomy saved at {outdir}truth.tsv")
    return list(dfTax.index), f"{outdir}truth.tsv"

def check_truth(checkedPath, truthPath, search=None) :
    """
    This function compares a predictions.tsv produced by WizardEye with a truth.tsv
    containing the record ids and the associated kingdoms. It can also calculate
    scores if a searched kingdom is provided.
    """
    predictions = pd.read_csv(checkedPath, sep='\t', header=0, index_col=0)
    truth = pd.read_csv(truthPath, sep='\t', header=0, index_col=0, dtype=str).loc[predictions.index]
    pd.concat([predictions, truth], axis=1).to_csv(checkedPath, sep="\t")

    sameAsTruth, countTable, sameAsTruthSummary, sameAsSearch, countTFPN = pd.DataFrame(index=predictions.index, columns=predictions.columns), pd.DataFrame(), pd.DataFrame(), pd.DataFrame(), pd.DataFrame(index=['TP', 'TN', 'FP', 'FN'])
    for col in predictions.columns :
        for index in predictions.index :
            sameAsTruth.loc[index, col] = (predictions.loc[index, col] == truth.loc[index, 'division'])
        countTable = pd.concat([countTable, predictions[col].value_counts().rename(index=col)], axis=1)
        if search != None :
            sameAsSearch[col] = (predictions[col] == search)

    if search != None :
        dfTFPN = sameAsTruth.astype(str) + sameAsSearch.astype(str)
        dfTFPN.replace("TrueTrue", "TP", inplace=True)
        dfTFPN.replace("TrueFalse", "TN", inplace=True)
        dfTFPN.replace("FalseTrue", "FP", inplace=True)
        dfTFPN.replace("FalseFalse", "FN", inplace=True)

    countTable = pd.concat([countTable, truth['division'].value_counts().rename(index='Truth')], axis=1)

    for col in predictions.columns :
        sameAsTruthSummary = pd.concat([sameAsTruthSummary, sameAsTruth[col].value_counts().rename(index=col)], axis=1)
        if search != None :
            countTFPN = pd.concat([countTFPN, dfTFPN[col].value_counts().rename(index=col)], axis=1)

    return countTable.fillna(0).astype(int), sameAsTruthSummary.fillna(0).astype(int), countTFPN.fillna(0).astype(int)