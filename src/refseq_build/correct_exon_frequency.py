import math
import os
import sys

sys.path.append("/home/weber/PycharmProjects/EXOTIC/src")

import pandas as pd
import gzip

pd.options.mode.chained_assignment = None  # default='warn'
import multiprocessing
import parmap
import numpy as np
import collections
from tqdm import tqdm


tqdm.pandas()
from pprint import pprint
import requests
import re
from utils.utils import load_config_file
import json
import subprocess
from pandarallel import pandarallel

pandarallel.initialize(nb_workers=20, progress_bar=True)

## YAML FILES CONFIG
yaml = load_config_file()
exotic_files = yaml

transcript_correction_file = pd.read_csv(exotic_files["GTEX"]["transcript_check_refseq"], compression="gzip", sep="\t")

transcripts_checked = set(transcript_correction_file.RefSeq_mRNAID.values.tolist())

refseq_precomputed = pd.read_csv(exotic_files["EXOTIC"]["refseq_path"], compression="gzip", sep="\t").drop_duplicates(
    subset=["Gene", "ranges", "HGNC"]
)


def compute_mrna_refseq():
    o_file = exotic_files["RefSeq"]["refseq_mrna_cds"]
    if os.path.isfile(o_file) is False:
        mrna_df = pd.read_csv(
            "/home/weber/PycharmProjects/ExoCarto/data/1_interim/RefSeq_GRCh37_complete_new.csv.gz",
            compression="gzip",
            sep="\t",
            # nrows=1000,
        )
        dict_mrna = collections.defaultdict(list)
        current_elem = str()

        l_mrna_cds = list()

        # FIRST PASS / RAW

        for j, row in tqdm(mrna_df.iterrows()):
            if row["Element"] == "gene":
                current_elem_gene = row["Name"]

            elif row["Element"] == "mRNA":
                current_elem_mrna = row["Name"]
            elif row["Element"] == "CDS":
                # dict_mrna[current_elem].append(row["Elem_position_ID"].replace("_", "-"))
                l_mrna_cds.append(
                    {
                        "Gene": current_elem_gene,
                        "mRNA": current_elem_mrna.split(".")[0],
                        "CDS": row["Elem_position_ID"].replace("_", "-"),
                    }
                )

        df = pd.DataFrame(l_mrna_cds)
        df_mrna_cds = df.groupby(["Gene", "CDS"])["mRNA"].apply(list).reset_index()
        df_mrna_cds = df_mrna_cds.rename({"mRNA": "mRNA_exons"}, axis=1)
        df_mrna_cds["mRNA_nb"] = df_mrna_cds["mRNA_exons"].apply(len)

        df_mrna_gene = df.groupby("Gene")["mRNA"].apply(set).reset_index()
        df_mrna_gene = df_mrna_gene.rename({"mRNA": "mRNA_gene"}, axis=1)
        df_mrna_gene["mRNA_nb_total"] = df_mrna_gene["mRNA_gene"].apply(len)

        merge_df_mrna_cds = pd.merge(df_mrna_cds, df_mrna_gene, on="Gene")

        merge_df_mrna_cds["mRNA_gene"] = merge_df_mrna_cds["mRNA_gene"].apply(list)

        merge_df_mrna_cds.to_parquet(o_file)
    else:
        merge_df_mrna_cds = pd.read_parquet(o_file)
    return merge_df_mrna_cds


# GTEx correction
def gtex_correction():
    if os.path.isfile(exotic_files["RefSeq"]["refseq_corrected_by_gtex"]) is False:

        refseq_corrected_transcripts = merge_df_mrna_cds[["Gene", "ranges", "mRNA_exons"]]
        print(refseq_corrected_transcripts.Gene.nunique())

        # EXPLODE FROM LISTS TO ROWS
        refseq_corrected_transcripts = refseq_corrected_transcripts.explode("mRNA_exons")
        refseq_corrected_transcripts["mRNA_exons"] = refseq_corrected_transcripts["mRNA_exons"].apply(lambda r: r.split(".")[0])

        # FILTERING mRNAs which respect GTEx expression cutoffs (reads & TPM)
        refseq_corrected_transcripts = refseq_corrected_transcripts.loc[refseq_corrected_transcripts["mRNA_exons"].isin(transcripts_checked)]

        # GROUPBY
        refseq_corrected_transcripts = refseq_corrected_transcripts.groupby(["Gene", "ranges"])["mRNA_exons"].apply(list)

        refseq_corrected_transcripts = refseq_corrected_transcripts.reset_index()
        refseq_corrected_transcripts = refseq_corrected_transcripts.rename({"mRNA_exons": "new_mRNA_exons"}, axis=1)
        refseq_corrected_transcripts["new_mRNA_nb"] = refseq_corrected_transcripts["new_mRNA_exons"].apply(len)

        # TMP FILE
        refseq_corrected_transcripts.to_csv(
            exotic_files["RefSeq"]["refseq_corrected_by_gtex"],
            compression="gzip",
            sep="\t",
            index=False,
        )

    else:
        refseq_corrected_transcripts = pd.read_csv(exotic_files["RefSeq"]["refseq_corrected_by_gtex"], compression="gzip", sep="\t")
    return refseq_corrected_transcripts


def process_new_file(refseq_corrected_transcripts):
    # COMPUTE mRNAs LIST AT GENE LEVEL BASED ON GROUPBY
    refseq_corrected_transcripts["new_mRNA_exons"] = refseq_corrected_transcripts["new_mRNA_exons"].apply(eval)
    refseq_corrected_transcripts_gene_level = refseq_corrected_transcripts.groupby("Gene")["new_mRNA_exons"].apply(list).reset_index()

    # EXTEND LIST OF LIST INTO ONE LIST TRANSFORMED TO A SET
    refseq_corrected_transcripts_gene_level["new_mRNA_gene"] = refseq_corrected_transcripts_gene_level["new_mRNA_exons"].apply(
        lambda r: list(set([sub_e for e in r for sub_e in e]))
    )

    # MERGE EXON & GENE LEVELS
    refseq_corrected_transcripts = pd.merge(
        refseq_corrected_transcripts, refseq_corrected_transcripts_gene_level[["Gene", "new_mRNA_gene"]], on="Gene"
    )

    # COMPUTE
    refseq_corrected_transcripts["new_mRNA_nb_total"] = refseq_corrected_transcripts["new_mRNA_gene"].apply(len)
    refseq_corrected_transcripts["new_Ratio"] = (
        refseq_corrected_transcripts["new_mRNA_nb"].astype(str) + "/" + refseq_corrected_transcripts["new_mRNA_nb_total"].astype(str)
    )
    refseq_corrected_transcripts["new_Ratio_num"] = refseq_corrected_transcripts["new_Ratio"].apply(eval)
    refseq_corrected_transcripts["new_Ratio_num"] = refseq_corrected_transcripts["new_Ratio_num"].astype(float)

    # CONST & ALT
    refseq_corrected_transcripts.loc[refseq_corrected_transcripts["new_Ratio_num"] < 1, "Const_Alt"] = "Alt"
    refseq_corrected_transcripts.loc[refseq_corrected_transcripts["new_Ratio_num"] == 1, "Const_Alt"] = "Const"

    # ALT BINS
    bins = [0, 0.2, 0.4, 0.6, 0.8, 1]
    labels = bins.copy()
    labels_ratio = [str(round(labels[j], 1)) + " - " + str(round(labels[j + 1], 1)) for j in range(len(labels) - 1)]
    refseq_corrected_transcripts["new_Ratio_num_bins"] = pd.cut(
        refseq_corrected_transcripts["new_Ratio_num"], bins=bins, labels=labels_ratio, include_lowest=True
    )

    # DUPLICATES
    refseq_corrected_transcripts = refseq_corrected_transcripts.drop_duplicates(subset=["Gene", "ranges"])

    # START, STOP, LENGTH
    refseq_corrected_transcripts["Start"] = refseq_corrected_transcripts["ranges"].apply(lambda r: r.split("-")[0])
    refseq_corrected_transcripts["Start"] = refseq_corrected_transcripts["Start"].astype(int)
    refseq_corrected_transcripts["End"] = refseq_corrected_transcripts["ranges"].apply(lambda r: r.split("-")[1])
    refseq_corrected_transcripts["End"] = refseq_corrected_transcripts["End"].astype(int)
    refseq_corrected_transcripts["Length"] = refseq_corrected_transcripts["End"] - refseq_corrected_transcripts["Start"]

    # CDS count
    t = refseq_corrected_transcripts.groupby(["Gene", "new_mRNA_nb_total"])["ranges"].agg("nunique").reset_index()

    # MERGE WITH PREVIOUS
    refseq_corrected_transcripts = pd.merge(
        refseq_corrected_transcripts, t.rename({"ranges": "new_CDS_count"}, axis=1), on=["Gene", "new_mRNA_nb_total"]
    )

    # OUPUTS
    refseq_corrected_transcripts.to_parquet(exotic_files["RefSeq"]["refseq_corrected_lite"], index=False)
    refseq_corrected_transcripts.to_csv(
        exotic_files["RefSeq"]["refseq_corrected_lite"].replace(".parquet", ".csv.gz"),
        compression="gzip",
        sep="\t",
        index=False,
    )


# FIRST PASS
merge_df_mrna_cds = compute_mrna_refseq()
merge_df_mrna_cds = merge_df_mrna_cds.rename({"CDS": "ranges"}, axis=1)

# CORRECT WITH GTEX
refseq_corrected_transcripts = gtex_correction()

#  MODIF & OUTPUT
process_new_file(refseq_corrected_transcripts)