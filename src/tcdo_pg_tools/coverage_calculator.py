#!/usr/bin/env python3
"""
author: Asher Preska Steinberg
compute AA seq coverage across multiple enzymes
"""
import pandas as pd
import glob as bob
import os
import numpy as np
from tqdm import tqdm
import itertools
from numba import njit
import click


#@njit
def compute_cov(starts, ends, AAlen):
    AAseq = np.zeros(AAlen)
    for i in np.arange(0, len(starts)):
        start = starts[i] - 1
        end = ends[i]
        AAseq[start:end] = 1
    return np.sum(AAseq) / AAlen * 100, np.sum(AAseq)


def enzyme_tag(combo):
    return ",".join(combo)


def compute_combo_cov(enzymes, proteindat, peptidedat):
    ### compute all enzyme combos

    combos = []
    for i in np.arange(1, len(enzymes) + 1):
        combo = [list(x) for x in itertools.combinations(enzymes, i)]
        combos.extend(combo)

    ### initialize coverage dataframe ....
    cov_df = pd.DataFrame()
    protein_lendat = proteindat[["Protein", "Length"]]
    for i in tqdm(np.arange(0, len(combos))):
        combo = combos[i]
        pepdat1 = peptidedat[peptidedat["enzyme"].isin(combo)]
        ### reduce sizes of these dataframes to the essentials ...
        pepdat1 = pepdat1[["Protein Start", "Protein End", "Protein", "enzyme"]]
        #### merge dataframes ....
        merge_pepdat1 = pd.merge(pepdat1, protein_lendat, on="Protein", how="left")
        grouped = merge_pepdat1.groupby(by="Protein")

        proteins = []
        covs = np.zeros(len(grouped))
        AAs = np.zeros(len(grouped))
        AAlens = np.zeros(len(grouped))
        j = 0
        for protein, group in grouped:
            starts = np.array(group["Protein Start"])
            ends = np.array(group["Protein End"])
            AAlen = np.array(group["Length"])[0]
            cov, AA = compute_cov(starts, ends, int(AAlen))
            proteins.append(protein)
            AAs[j] = AA
            AAlens[j] = AAlen
            covs[j] = cov
            j = j + 1
        cov_df1 = pd.DataFrame(zip(proteins, covs, AAs, AAlens),
                               columns=["Protein", "Coverage", "cum_peptide_length", "protein_length"])
        cov_df1["enzymes"] = enzyme_tag(combo)
        cov_df = pd.concat([cov_df, cov_df1])
    return cov_df

@click.command()
@click.option('--fragpipe_dir', required=True, type=click.Path(exists=True),
              help='fragpipe out dir; MUST have subdirs starting with enzyme (e.g., "trypsin_diaPASEF")')
@click.option('--enzymes', required=True, type=str,
                help='enzyme names, comma separated (e.g., "trypsin,chymotrypsin")')
@click.option("--output_tsv", required=False, default="AA_seq_cov.tsv", help="output tsv filename")
@click.option('--unique_proteins', is_flag=True, default=False, required=False,
              help="just calculate coverage for proteins with unique peptides")
@click.option('--unique_peptides', is_flag=True, default=False, required=False,
              help="just calculate coverage using unique peptides")
def coverage_calculator(fragpipe_dir, enzymes, output_tsv, unique_proteins, unique_peptides):
    """
    calculate wt mean AA coverage for multi-enzyme digests
    """
    # import enzymes ...
    enzymes = enzymes.split(',')
    archive = fragpipe_dir
    # load in protein data
    proteindat = pd.DataFrame()
    for i in np.arange(0, len(enzymes)):
        enzyme = enzymes[i]
        pattern = os.path.join(archive, f"{enzyme}*","**", "protein.tsv")
        paths = bob.glob(pattern, recursive=True)
        assert len(paths) == 1, \
            f"multiple protein.tsv files returned for {enzyme}; fix subdirectory naming scheme"
        path = paths[0]
        df = pd.read_csv(path, sep="\t")
        df["enzyme"] = enzyme
        proteindat = pd.concat([proteindat, df])

    # load in peptide data
    peptidedat = pd.DataFrame()
    for i in np.arange(0, len(enzymes)):
        enzyme = enzymes[i]
        pattern = os.path.join(archive, f"{enzyme}*","**", "peptide.tsv")
        paths = bob.glob(pattern, recursive=True)
        assert len(paths) == 1, \
            f"multiple peptide.tsv files returned for {enzyme}"
        path = paths[0]
        df = pd.read_csv(path, sep="\t")
        df["enzyme"] = enzyme
        peptidedat = pd.concat([peptidedat, df])
    # filter for uniquely detected proteins (if required)
    if unique_proteins:
        proteindat = proteindat[proteindat["Indistinguishable Proteins"].isna()]
        # filter peptide list to just include these proteins
        proteins = set(proteindat["Protein"])
        peptidedat = peptidedat[peptidedat["Protein"].isin(proteins)]
    # filter for uniquely detected peptides (if required)
    if unique_peptides:
        peptidedat = peptidedat[peptidedat["Mapped Proteins"].isna()]
    cov_df = compute_combo_cov(enzymes, proteindat, peptidedat)

    # compute cumulative coverage for each enzyme combo
    canon_cov_df = cov_df
    peplendf = canon_cov_df.pivot_table(columns="enzymes", index="Protein", values="cum_peptide_length")
    peplen_sums = peplendf.sum()
    totlendf = canon_cov_df.pivot_table(columns="enzymes", index="Protein", values="protein_length")
    totlen_sums = totlendf.sum()
    cum_cov = (peplen_sums/totlen_sums)*100
    cum_cov_df = pd.DataFrame(cum_cov)
    cum_cov_df = cum_cov_df.reset_index()
    columns=["enzymes", "wt. mean coverage (%AA)"]
    cum_cov_df.columns = columns
    cum_cov_df.to_csv(output_tsv, index=False, sep="\t")
    return

