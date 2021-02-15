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


## YAML FILES CONFIG
yaml = load_config_file()
exotic_files = yaml

directory_sqlseeker2 = exotic_files["QTL"]["sqlseeker2"]
pprint(directory_sqlseeker2)

for k, directory in directory_sqlseeker2.items():
    print(k)
    list_dir = [
        directory + t + "/sqtls-0.05fdr.permuted.tsv.gz"
        for t in os.listdir(directory)
        if "parquet" not in t
    ]

    def mp_tissue_sqlseeker(file, l):
        tissue_dir = file.split("/")[-2]
        df = pd.read_csv(file, compression="gzip", sep="\t")
        df["Tissue"] = tissue_dir
        df["variant_id"] = df["best.snp"].str.replace("_b38", "")
        df[["CHROM", "POS_GRCh38", "REF", "ALT"]] = df["variant_id"].str.split(
            "_", expand=True
        )

        l.append(df)

    m = multiprocessing.Manager()
    l = m.list()

    parmap.starmap(mp_tissue_sqlseeker, list(zip(list_dir)), l, pm_pbar=True)
    df = pd.concat(list(l))
    first_cols = ["CHROM", "POS_GRCh38", "POS_GRCh38", "REF", "ALT", "Tissue"]
    df = df[first_cols + [c for c in list(df.columns) if c not in first_cols]]
    modified_columns = list(df.columns)
    modified_columns[2] = modified_columns[2] + "_copie"
    df.columns = modified_columns
    df.to_csv(
        directory + "Merge_sQTL.csv.gz", index=False, compression="gzip", sep="\t"
    )
    print(df)