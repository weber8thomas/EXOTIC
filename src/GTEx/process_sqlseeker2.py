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
from liftover import ChainFile

converter = ChainFile("/gstock/biolo_datasets/variation/benchmark/Databases/CrossMap/GRCh38_to_GRCh37.chain.gz", "GRCh38", "GRCh37")
# print(converter["1"][103786442])
# exit()

## YAML FILES CONFIG
yaml = load_config_file()
exotic_files = yaml

directory_sqlseeker2 = {d: c for d, c in exotic_files["QTL"]["sqlseeker2"].items() if "8" in d}
pprint(directory_sqlseeker2)

for k, directory in directory_sqlseeker2.items():
    print(k)
    if os.path.isfile(directory + "Merge_sQTL.parquet") is True:
        list_dir = [directory + t + "/sqtls-0.05fdr.permuted.tsv.gz" for t in os.listdir(directory) if "parquet" not in t]

        def apply_liftover(r):
            output = converter[r["CHROM"]][int(r["POS_GRCh38"])]
            try:
                pos_grch37 = output[0][1]
            except IndexError:
                pos_grch37 = np.nan
            return pos_grch37

        def mp_tissue_sqlseeker(file, l):
            tissue_dir = file.split("/")[-2]
            df = pd.read_csv(file, compression="gzip", sep="\t")
            df["Tissue"] = tissue_dir
            df["variant_id"] = df["snpId"].str.replace("_b38", "")
            df[["CHROM", "POS_GRCh38", "REF", "ALT"]] = df["variant_id"].str.split("_", expand=True)
            l.append(df)

            # df["POS_GRCh37"] = df[["CHROM", "POS_GRCh38"]].apply(lambda r: apply_liftover(r), axis=1)
            # print(df.loc[df["POS_GRCh37"].isna() == True])
            # df["POS_GRCh37"] = df["POS_GRCh37"].astype(int)
            # print(df)

        m = multiprocessing.Manager()
        l = m.list()

        # mp_tissue_sqlseeker(list_dir[0], l)
        # exit()
        parmap.starmap(mp_tissue_sqlseeker, list(zip(list_dir)), l, pm_pbar=True)
        df = pd.concat(list(l))
        first_cols = ["CHROM", "POS_GRCh38", "REF", "ALT", "Tissue"]
        df = df[first_cols + [c for c in list(df.columns) if c not in first_cols]]
        # modified_columns = list(df.columns)
        # modified_columns[2] = modified_columns[2] + "_copie"
        # df.columns = modified_columns
        df.to_parquet(directory + "Merge_sQTL.parquet")
    else:
        df_merge = pd.read_csv(directory + "Merge_sQTL.csv.gz", compression="gzip", sep="\t", low_memory=False)
        df_37 = pd.read_csv(directory + "Merge_sQTL_lite_lite_37.csv.gz", compression="gzip", sep="\t", low_memory=False)
        df_37.columns = ["CHROM", "POS_GRCh37", "best.snp"]
        output_df = pd.merge(df_merge, df_37[["POS_GRCh37", "best.snp"]], on="best.snp")
        output_df.to_parquet(directory + "Merge_sQTL_final.parquet")
