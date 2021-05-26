import math
import os
import sys
import pandas as pd

sys.path.append("/home/weber/PycharmProjects/EXOTIC/src")
from statannot import add_stat_annotation

pd.options.mode.chained_assignment = None  # default='warn'
import multiprocessing
import parmap
import numpy as np
import collections
from tqdm import tqdm

tqdm.pandas()
import pandarallel
from pandarallel import pandarallel

pandarallel.initialize(nb_workers=20, progress_bar=True)
# tqdm.pandas()
from pprint import pprint
from scipy.stats import zscore
from scipy import stats

import requests
import re
import json
from utils.utils import *

## YAML FILES CONFIG
yaml = load_config_file(config_file="clean/src/config_clean_clean.yaml")

# JSON
dicts = json.load(open("clean/src/config/EXOTIC_config.json"))


class sQTLs:
    def __init__(self):
        dext = self.load_exotic_biomart(
            output_path=yaml["6_sQTLs"]["TMP"]["dext_with_enst"],
            dext_path=yaml["4_DEXT"]["Final"]["dext"],
            biomart_path=yaml["1_GENOMICS"]["External"]["biomart_37"],
        )
        self.map_sqtl_on_transcripts(
            path=yaml["6_sQTLs"]["Final"]["dext_sqtl_biomart_37"], dext=dext, sqtlseeker_dir=yaml["6_sQTLs"]["External"]["sqtlseeker_dir"]
        )

    def load_exotic_biomart(self, output_path, dext_path, biomart_path):
        print("# Add ENST columns to dext file / File : {}".format(output_path))

        if os.path.isfile(output_path) is False:
            dext = pd.read_parquet(dext_path)
            bins = [round(e, 2) for e in list(np.arange(0, 1.1, 0.1))]
            labels = convert_bins_into_labels(bins)

            dext["dext_down_reversed"] = dext["dext_down"] * (-1)
            dext["dext_bins_down"] = pd.cut(dext["dext_down_reversed"], bins=bins, labels=labels)
            dext["dext_bins_up"] = pd.cut(dext["dext_up"], bins=bins, labels=labels)
            dext["dext_tissues_up"] = dext[dicts["GTEx_tissues_list"] + ["dext_up"]].progress_apply(
                lambda r: [
                    dicts["GTEx_tissues_list"][j] for j, e in enumerate(r) if e == r["dext_up"] and j < len(dicts["GTEx_tissues_list"])
                ],
                axis=1,
            )
            dext["dext_tissues_down"] = dext[dicts["GTEx_tissues_list"] + ["dext_down"]].progress_apply(
                lambda r: [
                    dicts["GTEx_tissues_list"][j] for j, e in enumerate(r) if e == r["dext_down"] and j < len(dicts["GTEx_tissues_list"])
                ],
                axis=1,
            )
            cutoff_down = dext["dext_down_reversed"].quantile(0.95)
            cutoff_up = dext["dext_up"].quantile(0.95)

            # dext.loc[(dext['dext_down_reversed'] >= cutoff_down) | (dext['dext_up'] >= cutoff_up)]
            dext.loc[(dext["dext_down_reversed"] >= cutoff_down), "Down"] = True
            dext.loc[(dext["dext_up"] >= cutoff_up), "Up"] = True
            dext.loc[(dext["dext_up"] >= cutoff_up) & (dext["dext_down_reversed"] >= cutoff_down), ["Down", "Up"]] = True
            dext[["Down", "Up"]] = dext[["Down", "Up"]].fillna(False)

            biomart = pd.read_csv(biomart_path, compression="gzip", sep="\t").rename(
                {"Transcript stable ID": "ENST", "RefSeq mRNA ID": "mRNA_exons"}, axis=1
            )

            dext_exploded = dext.explode("mRNA_exons")
            dext_exploded = pd.merge(dext_exploded, biomart[["ENST", "mRNA_exons"]], on="mRNA_exons")
            dext = pd.merge(dext, dext_exploded[["MAP", "ENST"]].drop_duplicates().groupby("MAP")["ENST"].apply(list), on="MAP")
            dext.to_parquet(output_path, index=False)
        else:
            dext = pd.read_parquet(output_path)
        print(dext)
        return dext

    def map_sqtl_on_transcripts(self, path, dext, sqtlseeker_dir):
        print("# Map sQTL on dext exons / File : {}".format(path))

        if os.path.isfile(path) is False:

            sqtlseeker_listdir = sorted([d for d in os.listdir(sqtlseeker_dir) if "parquet" not in d])

            m = multiprocessing.Manager()
            l_df = m.list()

            dext = dext.explode("ENST")
            parmap.starmap(self.mp_sqtl, list(zip(sqtlseeker_listdir)), sqtlseeker_dir, l_df, dext, pm_pbar=True)

            dext_sqtl = pd.concat(list(l_df)).sort_values(by=["MAP"]).reset_index(drop=True)
            dext_sqtl.to_parquet(path)
            dext_sqtl.head()
        else:
            dext_sqtl = pd.read_parquet(path)
        print(dext_sqtl)
        return dext_sqtl

    @staticmethod
    def mp_sqtl(tissue, sqtlseeker_dir, l_df, dext):
        sqtlseeker_tmp = pd.read_csv(sqtlseeker_dir + tissue + "/sqtls-0.05fdr.permuted.tsv.gz", compression="gzip", sep="\t")
        sqtlseeker_tmp = sqtlseeker_tmp.melt(
            id_vars=[c for c in sqtlseeker_tmp.columns if "tr." not in c],
            value_vars=["tr.first", "tr.second"],
            var_name="tr_type",
            value_name="ENST",
        )
        sqtlseeker_tmp.ENST = sqtlseeker_tmp.ENST.apply(lambda r: r.split(".")[0])
        sqtlseeker_tmp["Tissue"] = tissue
        sqtlseeker_tmp["gene_id"] = sqtlseeker_tmp["geneId"].apply(lambda r: r.split(".")[0])
        sqtlseeker_tmp = sqtlseeker_tmp.loc[sqtlseeker_tmp["md"] >= 0.05]
        merge = pd.merge(dext, sqtlseeker_tmp, on="ENST")
        l_df.append(merge)

        # dext = pd.merge(dext, biomart)
        # dext.head()


if __name__ == "__main__":
    sQTLs()
