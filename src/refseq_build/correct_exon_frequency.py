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

transcript_correction_file = pd.read_csv(
    exotic_files["GTEX"]["transcript_check_refseq"], compression="gzip", sep="\t"
)

transcripts_checked = set(transcript_correction_file.RefSeq_mRNAID.values.tolist())

refseq_precomputed = pd.read_csv(
    exotic_files["EXOTIC"]["refseq_path"], compression="gzip", sep="\t"
).drop_duplicates(subset=["Gene", "ranges", "HGNC"])


def compute_mrna_refseq():
    o_file = exotic_files["RefSeq"]["refseq_mrna_cds"]
    if os.path.isfile(o_file) is True:
        mrna_df = pd.read_csv(
            "/home/weber/PycharmProjects/ExoCarto/data/1_interim/RefSeq_GRCh37_complete_new.csv.gz",
            compression="gzip",
            sep="\t",
            # nrows=1000,
        )
        dict_mrna = collections.defaultdict(list)
        current_elem = str()

        l_mrna_cds = list()

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
        merge_df_mrna_cds["Ratio"] = (
            merge_df_mrna_cds["mRNA_nb"].astype(str)
            + "/"
            + merge_df_mrna_cds["mRNA_nb_total"].astype(str)
        )

        merge_df_mrna_cds["mRNA_gene"] = merge_df_mrna_cds["mRNA_gene"].apply(list)

        merge_df_mrna_cds.to_parquet(o_file)
    else:
        merge_df_mrna_cds = pd.read_parquet(o_file)
    return merge_df_mrna_cds


merge_df_mrna_cds = compute_mrna_refseq()
merge_df_mrna_cds = merge_df_mrna_cds.rename({"CDS": "ranges"}, axis=1)


if os.path.isfile(exotic_files["RefSeq"]["refseq_corrected_by_gtex"]) is False:

    refseq_corrected_transcripts = merge_df_mrna_cds[["Gene", "ranges", "mRNA_exons"]]
    print(refseq_corrected_transcripts.Gene.nunique())

    # refseq_corrected_transcripts["mRNA_exons"] = refseq_corrected_transcripts[
    # "mRNA_exons"
    # ].apply(eval)
    refseq_corrected_transcripts = refseq_corrected_transcripts.explode("mRNA_exons")
    refseq_corrected_transcripts["mRNA_exons"] = refseq_corrected_transcripts[
        "mRNA_exons"
    ].apply(lambda r: r.split(".")[0])

    refseq_corrected_transcripts = refseq_corrected_transcripts.loc[
        refseq_corrected_transcripts["mRNA_exons"].isin(transcripts_checked)
    ]

    refseq_corrected_transcripts = refseq_corrected_transcripts.groupby(
        ["Gene", "ranges"]
    )["mRNA_exons"].apply(list)

    refseq_corrected_transcripts = refseq_corrected_transcripts.reset_index()
    refseq_corrected_transcripts = refseq_corrected_transcripts.rename(
        {"mRNA_exons": "new_mRNA_exons"}, axis=1
    )
    refseq_corrected_transcripts["new_mRNA_nb"] = refseq_corrected_transcripts[
        "new_mRNA_exons"
    ].apply(len)

    # refseq_corrected_transcripts.to_csv(
    # exotic_files["RefSeq"]["refseq_corrected_by_gtex"],
    # compression="gzip",
    # sep="\t",
    # index=False,
    # )

else:
    refseq_corrected_transcripts = pd.read_csv(
        exotic_files["RefSeq"]["refseq_corrected_by_gtex"], compression="gzip", sep="\t"
    )

refseq_corrected_transcripts_gene_level = (
    refseq_corrected_transcripts.groupby("Gene")["new_mRNA_exons"]
    .apply(list)
    .reset_index()
)

refseq_corrected_transcripts["new_mRNA_exons"] = refseq_corrected_transcripts[
    "new_mRNA_exons"
].apply(eval)

refseq_corrected_transcripts_gene_level[
    "new_mRNA_exons"
] = refseq_corrected_transcripts["new_mRNA_exons"].apply(
    lambda r: list(set([sub_e for e in r for sub_e in e.split(",")]))
)

refseq_corrected_transcripts_gene_level = (
    refseq_corrected_transcripts_gene_level.rename(
        {"new_mRNA_exons": "new_mRNA_gene"}, axis=1
    )
)
refseq_corrected_transcripts = pd.merge(
    refseq_corrected_transcripts, refseq_corrected_transcripts_gene_level, on="Gene"
)

merge_old_new_refseq = pd.merge(
    merge_df_mrna_cds, refseq_corrected_transcripts, on=["Gene", "ranges"]
)


merge_old_new_refseq["Same_nb_mRNA_exons"] = merge_old_new_refseq.parallel_apply(
    lambda r: True if r["mRNA_nb"] == r["new_mRNA_nb"] else False, axis=1
)

merge_old_new_refseq["new_mRNA_nb_total"] = merge_old_new_refseq["new_mRNA_gene"].apply(
    len
)


merge_old_new_refseq["Same_nb_mRNA_gene"] = merge_old_new_refseq.parallel_apply(
    lambda r: True if r["mRNA_nb_total"] == r["new_mRNA_nb_total"] else False, axis=1
)

merge_old_new_refseq["Unchanged"] = merge_old_new_refseq.parallel_apply(
    lambda r: True
    if r["Same_nb_mRNA_exons"] == True and r["Same_nb_mRNA_gene"] == True
    else False,
    axis=1,
)

print(merge_old_new_refseq)

merge_old_new_refseq.to_parquet(exotic_files["RefSeq"]["refseq_old_new_comparison"])
