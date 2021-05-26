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

biomart_grch37 = pd.read_csv(
    exotic_files["EXOTIC"]["biomart_ensembl_hgnc_refseq"],
    sep="\t",
)

biomart_grch37 = biomart_grch37[
    [
        "Transcript stable ID",
        "Transcript name",
        "Transcript type",
        "APPRIS annotation",
        "RefSeq mRNA ID",
        "Transcript support level (TSL)",
        "HGNC ID",
    ]
]

biomart_grch37.columns = [
    "transcript_id",
    "transcript_name",
    "transcript_biotype",
    "APPRIS",
    "RefSeq_mRNAID",
    "TSL",
    "HGNC_ID",
]


def check_file(filepath):
    check_type = "tpm" if "tpm" in filepath else "reads"
    print(check_type)
    o_file = filepath.replace(".txt.gz", "_checked_complete.parquet")
    o_file_lite = filepath.replace(".txt.gz", "_checked_lite.txt.gz")
    if os.path.isfile(o_file) is False:

        df = pd.read_csv(
            filepath,
            compression="gzip",
            sep="\t",
            # nrows=100,
        )

        df["gene_id"] = df["gene_id"].apply(lambda r: r.split(".")[0])

        df = pd.merge(biomart_grch37, df, on="gene_id")

        # df = df[list(df.columns)[:8]]

        cutoff_reads = 6
        cutoff_tpm = 0.1

        def check(r, check_type):
            i = 0
            if check_type == "reads":
                cutoff = cutoff_reads
            elif check_type == "tpm":
                cutoff = cutoff_tpm
            for e in r[7:]:
                if e >= cutoff:
                    i += 1
            return i

        df["check_{}_below_cutoff".format(check_type)] = df.parallel_apply(
            lambda r: check(r, check_type), axis=1
        )
        df["check_{}".format(check_type)] = df[
            "check_{}_below_cutoff".format(check_type)
        ].parallel_apply(lambda r: True if r / 11688 >= 0.2 else False)
        print(df)

        df.to_parquet(o_file, index=False)
        df = df[list(df.columns)[:7] + list(df.columns)[-2:]]
        df.to_csv(o_file_lite, index=False, compression="gzip", sep="\t")
    else:
        df = pd.read_csv(o_file_lite, compression="gzip", sep="\t")

    return df


def merge_check_files():
    if os.path.isfile(exotic_files["GTEX"]["transcript_check"]) is False:

        check_reads_file = check_file(exotic_files["GTEX"]["transcript_reads"])
        check_tpm_file = check_file(exotic_files["GTEX"]["transcript_tpm"])
        merge_check = pd.merge(
            check_reads_file,
            check_tpm_file,
            on=[
                "Chrom",
                "Start",
                "End",
                "Strand",
                "gene_id",
                "gene_symbol",
                "transcript_id",
            ],
        )
        merge_check["check"] = merge_check.apply(
            lambda r: True
            if set([r["check_reads"], r["check_tpm"]]) == {True}
            else False,
            axis=1,
        )
        merge_check.to_csv(
            exotic_files["GTEX"]["transcript_check"],
            compression="gzip",
            sep="\t",
            index=False,
        )
    else:
        merge_check = pd.read_csv(
            exotic_files["GTEX"]["transcript_check"],
            compression="gzip",
            sep="\t",
        )
    return merge_check


merge_check = merge_check_files()
merge_check["transcript_id"] = merge_check["transcript_id"].apply(
    lambda r: r.split(".")[0]
)
merge_check = pd.merge(biomart_grch37, merge_check, on="transcript_id")
merge_check = merge_check.loc[
    (merge_check["check"] == True) & (merge_check["RefSeq_mRNAID"].isna() == False)
]
merge_check.to_csv(
    exotic_files["GTEX"]["transcript_check_refseq"], compression="gzip", sep="\t"
)
print(merge_check)