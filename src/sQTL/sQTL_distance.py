# General imports
import os
import sys
import pandas as pd

pd.options.mode.chained_assignment = None  # default='warn'
import numpy as np
import subprocess

# Other imports
import multiprocessing
import parmap
import collections
from tqdm import tqdm

tqdm.pandas()
from pandarallel import pandarallel

pandarallel.initialize(nb_workers=60, progress_bar=True)
from pprint import pprint

# Custom utils
sys.path.append("/home/weber/PycharmProjects/EXOTIC/src")
from utils.utils import load_config_file

## YAML FILES CONFIG
yaml = load_config_file(config_file="/home/weber/PycharmProjects/EXOTIC/src/config.yaml")
exotic_files = yaml


# STEP 1 - prepare files


# TODO : add paths
def prepare_files():
    # TODO : paths to change

    exotic_37_38_path = "/gstock/EXOTIC/data/EXOTIC/EXOTIC_37_38.parquet"

    if os.path.isfile(exotic_37_38_path) is False:
        biomart = pd.read_csv("/gstock/EXOTIC/data/OTHERS/biomart_refseq_ensembl_hgnc.txt.gz", compression="gzip", sep="\t")
        biomart_strand = pd.read_csv("/gstock/EXOTIC/data/OTHERS/biomart_with_strand.txt.gz", compression="gzip", sep="\t")
        biomart = pd.merge(
            biomart, biomart_strand[["Gene stable ID", "Transcript stable ID", "Strand"]], on=["Gene stable ID", "Transcript stable ID"]
        )

        exotic = pd.read_parquet("/gstock/EXOTIC/data/EXOTIC/EXOTIC_modified_zscore.parquet")
        exotic["Exon_start"] = exotic.Exon.apply(lambda r: int(r.split("-")[0]))
        exotic["Exon_end"] = exotic.Exon.apply(lambda r: int(r.split("-")[1]))

        exotic = pd.merge(
            exotic,
            biomart[["Gene stable ID", "Chromosome/scaffold name"]]
            .drop_duplicates()
            .rename({"Chromosome/scaffold name": "CHROM", "Gene stable ID": "ensg"}, axis=1),
            on="ensg",
        )
        exotic["CHROM"] = exotic["CHROM"].astype(str)
        exotic[["CHROM", "Exon_start", "Exon_end", "MAP"]].sort_values(by=["CHROM", "Exon_start"]).to_csv(
            "/biolo/ngs/remap/EXOTIC_modified_zscore.bed", header=False, sep="\t", index=False
        )

        # TODO : handle remap_api.pl
        subprocess.call(
            "/usr/local/bin/perl remap_api.pl --mode asm-asm --from GCF_000001405.25 --dest GCF_000001405.26 --annotation /biolo/ngs/remap/EXOTIC_modified_zscore.bed --annot_out /biolo/ngs/remap/EXOTIC_modified_zscore_remap.bed --report_out report --in_format bed --out_format bed",
            shell=True,
        )
        exotic_remap_map = pd.read_csv("/biolo/ngs/remap/EXOTIC_remap.bed", sep="\t", names=["CHROM_38", "Exon_start_38", "Exon_end_38", "MAP"])
        exotic_remap_map = exotic_remap_map.loc[~exotic_remap_map["CHROM_38"].str.contains("HSCHR")]

        exotic_37_38 = pd.merge(exotic, exotic_remap_map, on="MAP")
        exotic_37_38["Gene"] = exotic_37_38.MAP.apply(lambda r: r.split("_")[0])
        exotic_37_38["MAP_38"] = (
            exotic_37_38["Gene"] + "_" + exotic_37_38["Exon_start_38"].astype(str) + "-" + exotic_37_38["Exon_end_38"].astype(str)
        )

        exotic_37_38.to_parquet(exotic_37_38_path)
    else:
        exotic_37_38 = pd.read_parquet(exotic_37_38_path)

    return exotic_37_38


## Merge EXOTIC & sQTL


def merge_exotic_sqtl_files(exotic_37_38):
    exotic_sqtl_38_path = "/gstock/EXOTIC/data/QTL/EXOTIC_sQTL_38_lite.parquet"
    if os.path.isfile(exotic_sqtl_38_path) is False:
        merge_exotic_sqtl = pd.read_parquet("/gstock/EXOTIC/data/QTL/sQTL_ENST_modified_zscore.parquet").reset_index(drop=True)

        merge_exotic_sqtl_38 = pd.merge(merge_exotic_sqtl, exotic_37_38[["MAP", "MAP_38"]], on="MAP")
        merge_exotic_sqtl_38.to_parquet(exotic_sqtl_38_path)
    else:
        merge_exotic_sqtl_38 = pd.read_parquet(exotic_sqtl_38_path)
        # merge_exotic_sqtl_38.sample(frac=1).head(100000).to_parquet(exotic_sqtl_38_path.replace(".parquet", "_lite.parquet"))
    return merge_exotic_sqtl_38


## CHECK EXON POSITION


def check_exon_position(r, nb_bin):
    l = [0] * nb_bin
    for j, b in enumerate(r["Gene_bins"]):
        if j < nb_bin:
            if r["Exon_start_38"] >= b and r["Exon_end_38"] <= r["Gene_bins"][j + 1]:
                l[j] += 1
    if l == [0] * nb_bin:
        for j, b in enumerate(r["Gene_bins"]):
            if j < nb_bin:
                if r["Exon_start_38"] < b < r["Exon_end_38"]:
                    l[j - 1] += 1
                    l[j] += 1
    if r["Strand"] == -1:
        l = l[::-1]

    return l


## COMPUTE EXON POSITION


def compute_exon_position(exotic_37_38, min_max="up", nb_bin=10, exotic_cutoffs=[0.5, 0.8, 0.9]):

    exotic_37_38_tmp = exotic_37_38.copy()
    exotic_37_38_tmp = exotic_37_38_tmp[
        [
            "CHROM",
            "MAP",
            "Exon_start_38",
            "Exon_end_38",
            "Gene start (bp)",
            "Gene end (bp)",
            "Strand",
            "EXOTIC_{}".format(min_max),
            "EXOTIC_bins_{}".format(min_max),
        ]
    ]

    l = list()

    for exotic_cutoff in exotic_cutoffs:
        exotic_37_38_tmp = exotic_37_38_tmp.loc[exotic_37_38_tmp["EXOTIC_{}".format(min_max)] > exotic_cutoff]
        #     exotic_37_38_tmp = exotic_37_38_tmp.head(3000)
        exotic_37_38_tmp["Gene_bin_size"] = (exotic_37_38_tmp["Gene end (bp)"] - exotic_37_38_tmp["Gene start (bp)"]) / nb_bin
        exotic_37_38_tmp["Gene_bins"] = exotic_37_38_tmp.progress_apply(
            lambda r: [r["Gene start (bp)"]] + [int(round((1 + e) * r["Gene_bin_size"], 0) + r["Gene start (bp)"]) for e in range(nb_bin)], axis=1
        )
        exotic_37_38_tmp["Exon_bin"] = exotic_37_38_tmp.progress_apply(lambda r: check_exon_position(r, nb_bin), axis=1)
        tmp_d = {
            k: v
            for k, v in dict(collections.Counter([",".join([str(sub_e) for sub_e in e]) for e in exotic_37_38_tmp.Exon_bin.values.tolist()])).items()
            if v > 10
        }
        tmp_d = [[v if int(e) == 1 else 0 for e in k.split(",")] for k, v in tmp_d.items()]
        tmp_df = pd.DataFrame(tmp_d)
        tmp_df.loc["Total"] = tmp_df.sum(axis=0)
        tmp_df.loc["Ratio"] = 100 * (tmp_df.loc["Total"] / tmp_df.loc["Total"].sum())
        tmp_df.columns = [str(e) for e in range(1, nb_bin + 1, 1)]
        l.append(tmp_df.tail(2))

    concat_df_distribution = pd.concat(l, keys=exotic_cutoffs, axis=0)

    # TODO : CHANGE PATH
    concat_df_distribution.to_excel("/gstock/EXOTIC/data/EXOTIC/EXOTIC_{}_density_exons_{}_bins.xlsx".format(min_max, str(nb_bin)))


def check_position_sqtl(r):

    l = [0] * nb_bin
    for j, b in enumerate(r["Gene_bins"]):
        if j < nb_bin:
            if b <= r["snpId_POS"] < r["Gene_bins"][j + 1]:
                l[j] += 1
    if r["Strand"] == -1:
        l = l[::-1]
    return l


def compute_sqtl_position_for_all():
    sqtl_deduplicated = merge_exotic_sqtl_38[["ensg", "symbol", "MAP_38", "snpId", "Tissue"]].drop_duplicates()
    sqtl_deduplicated["snpId_POS"] = sqtl_deduplicated["snpId"].apply(lambda r: int(r.split("_")[1]))
    sqtl_deduplicated = pd.merge(
        sqtl_deduplicated,
        exotic_37_38_tmp[["Gene_bins", "Exon_bin", "MAP_38", "Strand", "EXOTIC_up", "EXOTIC_down", "EXOTIC_tissues_up", "EXOTIC_tissues_down"]],
        on="MAP_38",
    )
    sqtl_deduplicated["snpId_bin"] = sqtl_deduplicated.parallel_apply(check_position_sqtl, axis=1)


def compute_sqtl_position_match_exotic(match_tissue=True, min_max="up", exotic_cutoffs=[0.5, 0.8, 0.9]):

    tmp_l = list()

    for exotic_cutoff in tqdm(exotic_cutoffs):
        sqtl_deduplicated_up = sqtl_deduplicated.loc[sqtl_deduplicated["EXOTIC_up"] > exotic_cutoff].explode("EXOTIC_tissues_up")
        sqtl_deduplicated_up.loc[sqtl_deduplicated_up["Tissue"] == sqtl_deduplicated_up["EXOTIC_tissues_up"], "Match_tissue"] = True
        sqtl_deduplicated_up.loc[sqtl_deduplicated_up["Tissue"] != sqtl_deduplicated_up["EXOTIC_tissues_up"], "Match_tissue"] = False
        sqtl_deduplicated_up.head()

        tmp_d = {
            k: v
            for k, v in dict(
                collections.Counter(
                    [
                        ",".join([str(sub_e) for sub_e in e])
                        for e in sqtl_deduplicated_up.loc[sqtl_deduplicated_up["Match_tissue"] == match_tissue][["snpId", "snpId_bin"]]
                        .drop_duplicates(subset=["snpId"])
                        .snpId_bin.values.tolist()
                    ]
                )
            ).items()
            if v > 10
        }
        tmp_d = [[v if int(e) == 1 else 0 for e in k.split(",")] for k, v in tmp_d.items()]
        tmp_df = pd.DataFrame(tmp_d)
        tmp_df.loc["Total"] = tmp_df.sum(axis=0)
        tmp_df.loc["Ratio"] = 100 * (tmp_df.loc["Total"] / tmp_df.loc["Total"].sum())
        tmp_df.columns = [str(e) for e in range(1, nb_bin + 1, 1)]

    concat_df_distribution = pd.concat(tmp_l, keys=exotic_cutoffs, axis=0)
    concat_df_distribution.to_excel(
        "/gstock/EXOTIC/data/EXOTIC/sQTL_EXOTIC_{}_match-tissue_{}_tissue_density_exons_{}_bins.xlsx".format(min_max, str(match), str(nb_bin))
    )


if __name__ == "__main__":
    exotic_37_38 = prepare_files()
    compute_exon_position(exotic_37_38, min_max="up", nb_bin=10, exotic_cutoffs=[0.5, 0.8, 0.9])
    print(merge_exotic_sqtl_files(exotic_37_38))