#!/usr/bin/env python3
"""
author: Asher Preska Steinberg
merge multiple fasta on sequence identity
"""
import pandas as pd
import glob as bob
import os
import numpy as np
from tqdm import tqdm
import itertools
from numba import njit
import click
from marsilea.upset import Upset, UpsetData
import matplotlib.pyplot as plt

def fasta2df(uniprotfastapath, sample="swissprot"):
    """
    fasta to pandas dataframe
    """
    with open(uniprotfastapath, "r") as file:
        uniprotfasta = file.readlines()
    ### make it a dataframe ...
    status = [] ## canonical, non-canonical, fusion ...
    IDs = []
    seqs = []
    i = 0
    for line in uniprotfasta:
        if line.startswith(">"):
            terms = line.split(" ")
            ID = terms[0]
            _, ID = ID.split(">")
            IDs.append(ID.strip())
        else:
            seqs.append(line.strip())
        status.append(sample)
    uniprotseqdat = pd.DataFrame(zip(IDs, status, seqs), columns = ["protein", "sample", "seq"])
    uniprotseqdat.reset_index(drop=True, inplace=True)
    return uniprotseqdat

def joinset(IDs, sort=False):
    ID = list(set(IDs))
    if sort:
        ID.sort()
    ID = ','.join(ID)
    return ID

def plot_upset(countdat, upset_path):
    upset_data = UpsetData.from_memberships(countdat.conditions.str.split(","))
    us = Upset(upset_data, min_cardinality=15)
    us.render()
    plt.savefig(upset_path, dpi=300, bbox_inches="tight")
    return

@click.command()
@click.option('-i', '--input_csv', required=True, type=click.Path(exists=True),
              help='three column format csv (path: fasta path, name: sample name, condition: condition)')
@click.option('-t', '--info_table', required=False,
              default='info_table.tsv',
              type=click.Path(), help="Path to index tsv for merged protein IDs")
@click.option('-fa','--merged_fasta', required=False,
              type=click.Path(), default='merged.fasta',
              help="Path to merged fasta file")
@click.option('--upset', is_flag=True, default=False, help="plot upset")
@click.option('--upset_path', required=False,
              type=click.Path(), default='upset_plot.svg',
              help="Path to upset plot")
def fasta_merge(input_csv, info_table, merged_fasta, upset, upset_path):
    """
    merge multiple fasta on sequence identity
    Args:
        input_csv: input csv with two columns (path: path to fasta, name: sample name, condition: condition (e.g., tumor, normal))
        info_table: path to index tsv for merged protein IDs (default: info_table.tsv)
        merged_fasta: path to merged fasta file (default: merged.fasta)
        upset: plot upset for different conditions (default: False)
        upset_path: path to upset plot (default: upset_plot.svg)
    Returns:

    """
    # read in metadata
    metadata = pd.read_csv(input_csv)
    # initialize dataframe for merging fasta
    protein_dat = pd.DataFrame()
    print("loading fasta files ...")
    for _, row in metadata.iterrows():
        fasta = row["fasta"]
        sample = row["sample"]
        condition = row["condition"]
        # load in the protein fasta file as well
        seqdat = fasta2df(fasta, sample=sample)
        seqdat["condition"] = condition
        protein_dat = pd.concat([protein_dat, seqdat])
    # perform groupby
    grouped = protein_dat.groupby(by=["seq"])
    # reorganize to some comprehensible output
    data = []
    i = 1
    print("creating final dataframe ...")
    for seq, group in tqdm(grouped):
        # get samples, conditions, protein_ids, number of occurrences
        samples = joinset(group["sample"])
        conditions = joinset(group["condition"], sort=True)
        protein_ids = joinset(group["protein"])
        # store data
        data.append({
            "sequence": seq,
            "Protein_ids": protein_ids,
            "unique_protein_id": f"PG{i}", # give protein a unique identifier
            "samples": samples,
            "conditions": conditions,
            "sample_count": len(group)
        })
        i = i+1
    # write dataframe to tsv
    countdat = pd.DataFrame(data)
    countdat.to_csv(info_table, sep="\t", index=False)
    # write to merged fasta file
    with open(merged_fasta, "w+") as f:
        for _, row in countdat.iterrows():
            id = row["unique_protein_id"]
            seq = row["sequence"]
            f.write(f">{id}\n")
            f.write(f"{seq}\n")
    # plot upset plot:
    if upset:
        plot_upset(countdat, upset_path)
    return
