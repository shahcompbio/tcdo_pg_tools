#!/usr/bin/env python3
"""
author: Asher Preska Steinberg
merge proteomegenerator fasta and results across multiple samples on AA seq identity
"""
import os
import re
import pandas as pd
import click
from tqdm import tqdm
from importlib.resources import files
from marsilea.upset import Upset, UpsetData
import matplotlib.pyplot as plt

def fasta2df(uniprotfastapath, sample="swissprot"):
    """
    fasta to pandas dataframe
    """
    with open(uniprotfastapath, "r") as file:
        uniprotfasta = file.readlines()
    ### make it a dataframe ...
    seqs = []
    data = []
    current_seq = ""
    for line in uniprotfasta:
        if line.startswith(">"):
            # if we already have a sequence, save it
            if current_seq:
                seqs.append(current_seq)
                current_seq = ""
            # now parse the new header
            terms = line.split(" ")
            ID = terms[0]
            _, ID = ID.split(">")
            # extract protein info
            id_terms = ID.split("|")
            if len(id_terms) == 3:
                db, UniqueID, EntryName = id_terms
            else:
                db, UniqueID, EntryName = "", "", ""
            # extract gene name
            geneName = ""
            for term in terms:
                if term.startswith("GN="):
                    _, geneName = term.split("=")
                    break
            data.append({
                "db": db,
                "protein": ID.strip(),
                "accession_number": UniqueID,
                "entry_name": EntryName,
                "gene_name": geneName.strip(),
                "sample": sample,
                "header": line
            })
        else:
            # accumulate sequence lines
            current_seq += line.strip()
    # after the loop, don’t forget the last record
    if current_seq:
        seqs.append(current_seq)
    seqdat = pd.DataFrame(data)
    seqdat["seq"] = seqs
    seqdat.reset_index(drop=True, inplace=True)
    return seqdat

def joinset(IDs, sort=False):
    ID = list(set(IDs))
    if sort:
        ID.sort()
    ID = ','.join(ID)
    return ID

def get_protein_info(group):
    # check if we have a swissprot protein:
    sp_df = group[group["protein"].str.contains("sp|")]
    if len(sp_df) > 0:
        uniprot_id = sp_df.iloc[0]["accession_number"]
        uniprot_entryname = sp_df.iloc[0]["entry_name"]
        gene_name = sp_df.iloc[0]["gene_name"]
        header = sp_df.iloc[0]["header"]
    else:
        uniprot_id = ""
        uniprot_entryname = group.iloc[0]["entry_name"]
        gene_name = group.iloc[0]["gene_name"]
        header = ""
    return uniprot_id, uniprot_entryname, gene_name, header


def plot_upset(countdat, upset_path):
    upset_data = UpsetData.from_memberships(countdat.conditions.str.split(","))
    us = Upset(upset_data, min_cardinality=15)
    us.render()
    plt.savefig(upset_path, dpi=300, bbox_inches="tight")
    return

def merge_proteome(input_csv, info_table, merged_fasta, upset,
                   upset_path, unique_proteins=True, filter="", filter_crap=""):
    """
    merge proteomegenerator fasta/results across multiple samples on AA seq identity
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
        # filter for unique proteins
        # get list of samples where no peptide tsv was provided
        no_quant = []
        if unique_proteins:
            philosopher_path = row["protein_table"]
            if type(philosopher_path) is not float and os.path.exists(philosopher_path):
                philosopher_dat = pd.read_csv(philosopher_path, sep="\t")
                philosopher_dat =  philosopher_dat[ philosopher_dat["Indistinguishable Proteins"].isna()]
                unique_proteins = list(philosopher_dat["Protein"])
                seqdat = seqdat[seqdat["protein"].isin(unique_proteins)]
            else:
                print(f"no protein tsv file provided for sample/condition: {sample}/{condition}")
                no_quant.append(sample)
        # append sample to dataframe
        protein_dat = pd.concat([protein_dat, seqdat])
    # filter out contaminants proteins
    if filter_crap != "":
        crap_dat = fasta2df(filter_crap, sample="crap")
        crap_seqs = list(crap_dat["seq"])
        # filter out contaminants
        print("filtering contaminants")
        protein_dat = protein_dat[~protein_dat["seq"].isin(crap_seqs)]
    # group proteins by sequence
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
        # get protein info
        unique_identifier, uniprot_entryname, gene_name, header = get_protein_info(group)
        if unique_identifier == "":
            unique_identifier = f"PG{i}"
            i += 1
        # store data
        data.append({
            "sequence": seq[0],
            "protein_ids": protein_ids,
            "unique_identifier": unique_identifier, # give protein a unique identifier
            "protein_name": uniprot_entryname,
            "gene_name": gene_name,
            "samples": samples,
            "conditions": conditions,
            "sample_count": len(set(group["sample"])),
            "header": header
        })
        i = i+1
    # write dataframe to tsv
    countdat = pd.DataFrame(data)
    if unique_proteins:
        countdat = countdat[~countdat["samples"].isin(no_quant)]
    # also filter for decoys and contams
    if filter != "":
        prefixes = filter.split(",")
        pattern = "|".join(re.escape(p) for p in prefixes)
        countdat = countdat[~countdat["protein_ids"].str.contains(pattern)]
    # write to merged fasta file
    with open(merged_fasta, "w+") as f:
        for _, row in countdat.iterrows():
            header = row["header"]
            seq = row["sequence"]
            if header == "":
                id = row["unique_identifier"]
                protein = row["protein_name"]
                gene = row["gene_name"]
                header = f">tr|{id}|{protein} PG3 predicted ORF OS=Homo sapiens OX=9606 GN={gene} PE=2\n"
            f.write(header)
            f.write(f"{seq}\n")
    # write dataframe to csv
    countdat = countdat.drop("header", axis=1)
    countdat.to_csv(info_table, sep="\t", index=False)
    # plot upset plot:
    if upset:
        plot_upset(countdat, upset_path)
    return

@click.command()
@click.option('-i', '--input_csv', required=True, type=click.Path(exists=True),
              help='four column csv (fasta: fasta path, '
                   'protein_table: protein.tsv path, '
                   'name: sample name, condition: condition)')
@click.option('-t', '--info_table', required=False,
              default='info_table.tsv',
              show_default=True,
              type=click.Path(), help="Path to index tsv for merged protein IDs")
@click.option('-fa','--merged_fasta', required=False,
              type=click.Path(), default='merged.fasta',
              show_default=True,
              help="Path to merged fasta file")
@click.option('--upset', is_flag=True, default=False, help="plot upset")
@click.option('--upset_path', required=False,
              type=click.Path(), default='upset_plot.svg',
              show_default=True,
              help="Path to upset plot")
@click.option('--filter_by_header',
              default="contam_,rev_,tr|GF",
              show_default=True,
              help="filter out proteins by header prefix (provide comma separated list)")
@click.option('--filter_crap',
              default=files("tcdo_pg_tools").joinpath("250707.crap.fasta"),
              show_default=True,
              help="filter out contaminants")
def merge_pg_results(input_csv, info_table, merged_fasta, upset, upset_path, filter_by_header, filter_crap):
    """
    merge proteomegenerator results across multiple samples on AA seq identity
    """
    return merge_proteome(input_csv, info_table, merged_fasta, upset,
                          upset_path, unique_proteins=True, filter=filter_by_header, filter_crap=filter_crap)

@click.command()
@click.option('-i', '--input_csv', required=True, type=click.Path(exists=True),
              help='three column csv (fasta: fasta path, '
                   'name: sample name, condition: condition)')
@click.option('-t', '--info_table', required=False,
              default='info_table.tsv',
              show_default=True,
              type=click.Path(), help="Path to index tsv for merged protein IDs")
@click.option('-fa','--merged_fasta', required=False,
              type=click.Path(), default='merged.fasta',
              show_default=True, help="Path to merged fasta file")
@click.option('--upset', is_flag=True, default=False, help="plot upset")
@click.option('--upset_path', required=False,
              type=click.Path(), default='upset_plot.svg',
              show_default=True, help="Path to upset plot")
@click.option('--filter_by_header',
              default="contam_,rev_,tr|GF",
              show_default=True,
              help="filter out proteins by header prefix (provide comma separated list)")
@click.option('--filter_crap',
              default=files("tcdo_pg_tools").joinpath("250707.crap.fasta"),
              show_default=True,
              help="filter out contaminants")
def merge_fasta(input_csv, info_table, merged_fasta, upset, upset_path, filter_by_header, filter_crap):
    """
    merge multiple fasta on sequence identity
    """
    return merge_proteome(input_csv, info_table, merged_fasta, upset, upset_path,
                          unique_proteins=False, filter=filter_by_header,filter_crap=filter_crap)

if __name__ == "__main__":
    merge_fasta()  # or merge_fasta() if you're testing that

