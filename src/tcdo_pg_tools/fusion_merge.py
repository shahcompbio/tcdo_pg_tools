#!/usr/bin/env python3
# author: Asher Preska Steinberg
import pandas as pd
import numpy as np
import click

@click.command()
@click.option('-i', '--input_metadata', required=True, type=click.Path(exists=True),
              help='Metadata sheet from isabl_utils with paths to FusionInspector results')
@click.option('-v', '--app_version', required=True, type=str,
                help='Filter for results from a specific Isabl app version')
@click.option('-r', '--result_type', required=True, type=str,
              help='Filter on result type in metadata sheet')
@click.option('-t', '--fusion_table', required=True,
              type=click.Path(), help="Path to output fusion info table")
@click.option('-fa','--output_fasta', required=True, type=click.Path(),
              help="Path to output FASTA file")
def fusion_merge(input_metadata, app_version, result_type, fusion_table, output_fasta):
    """
    Merge fusion calls across multiple samples based on gene symbols and breakpoints.
    """
    result_csv = input_metadata
    app_version = app_version
    result_type = result_type
    fusion_info_path = fusion_table
    output_fasta = output_fasta
    # for testing
    # app_version = "0.0.3" # filter for app version
    # result_type = "fusion_predictions" # name of result type to grab
    # result_csv = "packages/fusion_merge/tests/250304_ctat-lr-fusion.csv" # csv with paths to results
    # outputs
    # fusion_info_path = "packages/fusion_merge/tests/fusion_info.tsv"
    # output_fasta = "packages/fusion_merge/tests/fusion_proteins.fa"
    # get paths
    metadata = pd.read_csv(result_csv)
    paths = metadata[metadata["result_type"] == result_type]
    paths = paths[paths["isabl_version"] == app_version]
    # make a dataframe for fusions
    fusion_df = pd.DataFrame()
    for _, row in paths.iterrows():
        patient_id = row["isabl_patient_id"].split(",")[0]  # in case we did ONT+ILL
        isabl_pk = row["isabl_pk"]  # these will uniquely identify results
        inspectordat = pd.read_csv(row["result_filepath"], sep="\t")
        inspectordat["patient_id"] = patient_id
        inspectordat["isabl_pk"] = isabl_pk
        fusion_df = pd.concat([fusion_df, inspectordat])

    # group by fusion names and predicted ORFs and breakpoints
    fusion_groups = fusion_df.groupby(by=["#FusionName",
                                          "FUSION_TRANSL",
                                          "CDS_LEFT_RANGE",
                                          "CDS_LEFT_ID",
                                          "CDS_RIGHT_ID"])
    # make a dataframe that contains fusion info
    data = []
    i = 1
    for (genes, ORF, cds_lr, cds1_id, cds2_id), group in fusion_groups:
        # get breakpoint in AA position; will be useful when we match to proteomics later
        if not cds_lr == ".":
            AA_brk = np.round(int(cds_lr.split("-")[1]) / 3)
            protein_id = f"GF{i}"
            i = i + 1
        else:
            AA_brk = cds_lr
            protein_id = "."
        data.append({
            '#FusionName': genes,
            'CDS_LEFT_ID': cds1_id,
            'CDS_RIGHT_ID': cds2_id,
            'FUSION_TRANSL': ORF,
            'AA_brk_pos': AA_brk,
            'patients': list(group['patient_id']),
            'isabl_pks': list(group['isabl_pk']),
            'LR_FFPM': list(group['LR_FFPM']),
            'Protein': protein_id
        })
    # make the fusion table
    fusion_info_table = pd.DataFrame(data)
    fusion_info_table.to_csv(fusion_info_path, sep="\t", index=None)
    # now let's write the fasta file
    # drop rows without CDS regions
    fusion_info_table1 = fusion_info_table[fusion_info_table["Protein"].str.startswith("GF")]
    # write output fasta
    with open(output_fasta, "w+") as outfile:
        for _, row in fusion_info_table1.iterrows():
            protein_id = row["Protein"]
            gene = row["#FusionName"]
            ORF_id = row["CDS_LEFT_ID"] + "--" + row["CDS_RIGHT_ID"]
            AAseq = row["FUSION_TRANSL"]
            # write header and AA seq
            header = f">tr|{protein_id}|{ORF_id} PG3 predicted ORF OS=Homo sapiens OX=9606 GN={gene} PE=2\n"
            outfile.write(header)
            outfile.write(f"{AAseq}\n")


