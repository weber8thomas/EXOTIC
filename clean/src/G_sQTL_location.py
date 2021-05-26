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

pandarallel.initialize(nb_workers=20, progress_bar=False)
# tqdm.pandas()
from pprint import pprint
from scipy.stats import zscore
from scipy import stats
import subprocess

import requests
import re
import json
from utils.utils import *

## YAML FILES CONFIG
yaml = load_config_file(config_file="clean/src/config_clean_clean.yaml")

# JSON
dicts = json.load(open("clean/src/config/EXOTIC_config.json"))


# ! Remap gene start End
# ! Remap dext exons start End
# ! Compute bins
# ! sQTL bin location at all thresholds
# ! dext exon bin location at all thresholds


class sQTL_location:
    def __init__(self):
        dext_remap = self.remap_exons(
            dext_path=yaml["4_DEXT"]["Final"]["dext"],
            refseq_genes_coord_path=yaml["1_GENOMICS"]["TMP"]["tmp_refseq_concat_exons_cds_filtered"],
            input_bed=yaml["6_sQTLs"]["TMP"]["dext_bed_37"],
            output_bed=yaml["6_sQTLs"]["TMP"]["dext_bed_38"],
        )

        genes_remap = self.remap_genes(
            refseq_genes_coord_path=yaml["1_GENOMICS"]["TMP"]["tmp_refseq_concat_exons_cds_filtered"],
            input_bed=yaml["6_sQTLs"]["TMP"]["genes_bed_37"],
            output_bed=yaml["6_sQTLs"]["TMP"]["genes_bed_38"],
            dext_path=yaml["4_DEXT"]["Final"]["dext"],
        )
        # genes_remap_with_bins = self.compute_bins(
        # genes_remap=genes_remap,
        # )

        dext_remap_final = self.merge_remap(dext_path=yaml["4_DEXT"]["Final"]["dext"], genes=genes_remap, exons=dext_remap)

        distribution_location_exons = self.compute_exon_bin(
            dext=dext_remap_final, output_path=yaml["6_sQTLs"]["Figures_data"]["exons_location"]
        )

        # distribution_location_all_sqtls = self.compute_all_sqtl_bin(
        #     dext=dext_remap_final,
        #     output_path=yaml["6_sQTLs"]["Figures_data"]["sqtls_location"],
        #     sqtl_path=yaml["6_sQTLs"]["Final"]["dext_sqtl"],
        # )

        distribution_location_all_sqtls = self.compute_sqtl_bin(
            dext=dext_remap_final,
            output_path=yaml["6_sQTLs"]["Figures_data"]["match_sqtls_location"],
            sqtl_path=yaml["6_sQTLs"]["Final"]["dext_sqtl"],
        )

    def remap_genes(self, refseq_genes_coord_path, input_bed, output_bed, dext_path):

        if os.path.isfile(output_bed) is False:

            refseq_with_genes_coord = (
                pd.read_parquet(refseq_genes_coord_path)[["Gene", "Gene_start", "Gene_stop"]].drop_duplicates().sort_values(by="Gene")
            )
            dext = pd.read_parquet(dext_path)
            refseq_with_genes_coord = pd.merge(refseq_with_genes_coord, dext[["CHROM", "Gene"]], on="Gene")[
                ["CHROM", "Gene_start", "Gene_stop", "Gene"]
            ]

            refseq_with_genes_coord.to_csv(input_bed, header=False, sep="\t", index=False)

            self.remap_subprocess_fct(input_bed=input_bed, output_bed=output_bed)

        genes_remap = pd.read_csv(output_bed, sep="\t", names=["CHROM_38", "Gene_start_38", "Gene_end_38", "Gene"])
        # genes_remap = genes_remap.loc[~genes_remap["CHROM_38"].str.contains("HSCHR")].drop_duplicates().reset_index(drop=True)
        genes_remap = genes_remap.drop_duplicates().reset_index(drop=True)
        genes_remap["Gene_length"] = genes_remap["Gene_end_38"] - genes_remap["Gene_start_38"]

        # print(genes_remap)
        return genes_remap

        # genes_remap = pd.read_csv(output_bed, sep="\t", names=["CHROM_38", "Gene_start_38", "Gene_end_38", "Gene"])

        # print(refseq_with_genes_coord)
        # dext = pd.merge(refseq_with_genes_coord, dext, on="Gene")

    def remap_exons(self, dext_path, refseq_genes_coord_path, input_bed, output_bed):

        print("REMAP")

        if os.path.isfile(output_bed) is False:

            dext = pd.read_parquet(dext_path)
            dext[["CHROM", "Start", "End", "MAP"]].sort_values(by=["CHROM", "Start"]).to_csv(input_bed, header=False, sep="\t", index=False)

            self.remap_subprocess_fct(input_bed=input_bed, output_bed=output_bed)

        dext_remap = pd.read_csv(output_bed, sep="\t", names=["CHROM_38", "Exon_start_38", "Exon_end_38", "MAP"])
        dext_remap["CHROM_38"] = dext_remap["CHROM_38"].astype(str)
        # dext_remap = dext_remap.loc[~dext_remap["CHROM_38"].str.contains("HSCHR")].drop_duplicates().reset_index(drop=True)
        dext_remap = dext_remap.drop_duplicates().reset_index(drop=True)
        dext_remap["MAP_38"] = (
            dext_remap["MAP"].apply(lambda r: r.split("_")[0])
            + "_"
            + dext_remap["Exon_start_38"].astype(str)
            + "_"
            + dext_remap["Exon_end_38"].astype(str)
        )
        # print(dext_remap)
        return dext_remap

    @staticmethod
    def remap_subprocess_fct(input_bed, output_bed, perl_location="/usr/local/bin/perl", remap_location="/biolo/ngs/remap/remap_api.pl"):

        subprocess_args = [
            perl_location,
            remap_location,
            "--mode",
            "asm-asm",
            "--from",
            "GCF_000001405.25",
            "--dest",
            "GCF_000001405.26",
            "--annotation",
            input_bed,
            "--annot_out",
            output_bed,
            "--report_out",
            "report",
            "--in_format",
            "bed",
            "--out_format",
            "bed",
        ]
        p1 = subprocess.Popen(subprocess_args)
        p1.wait()

    def compute_bins(self, genes_remap, nb_bins=10):
        # nb_bins = 10
        genes_remap["Gene_bin_size"] = genes_remap["Gene_length"] / nb_bins
        genes_remap["Gene_bin_size"] = genes_remap["Gene_bin_size"]
        genes_remap["Gene_bins"] = genes_remap.apply(
            lambda r: [r["Gene_start_38"] + (round(r["Gene_bin_size"] * j)) for j in list(range(nb_bins + 1))], axis=1
        )
        # pprint(genes_remap.loc[0].to_dict())
        # print(genes_remap)
        return genes_remap

    def merge_remap(self, dext_path, genes, exons):
        dext = pd.read_parquet(dext_path)
        merge = pd.merge(pd.merge(dext, exons.drop(["CHROM_38"], axis=1), on="MAP"), genes, on="Gene")
        return merge

    @staticmethod
    def compute_exon_bin_fct(r, nb_bin):
        l = [0] * nb_bin
        # print(l)
        for j, b in enumerate(r["Gene_bins"]):
            # print(j, b)
            if j < nb_bin:
                if r["Exon_start_38"] >= b and r["Exon_end_38"] <= r["Gene_bins"][j + 1]:
                    l[j] += 1

        #         if j < (nb_bin - 1):
        #             if r['Exon_start_38'] >= b and r['Exon_end_38'] < r['Gene_bins'][j+2]:
        if l == [0] * nb_bin:
            for j, b in enumerate(r["Gene_bins"]):
                if j < nb_bin:
                    if r["Exon_start_38"] < b < r["Exon_end_38"]:
                        l[j - 1] += 1
                        l[j] += 1
        #                     print(j, r['Gene_bins'][j-1], b,  r['Exon_start_38'], r['Exon_end_38'])

        if r["Strand"] == -1:
            # print("Strand")
            # exit()
            l = l[::-1]

        return l

    def compute_exon_bin(self, output_path, dext):
        print("# Compute sQTL location for all bin number given in config / File : {}".format(output_path))

        if os.path.isfile(output_path) is False:

            dext["dext_down_reversed"] = dext["dext_down"] * -1
            dext["Strand"] = dext["Strand"].replace({0: -1})
            # print(dext)
            # print(list(dext.columns))
            # exit()

            l = list()

            dext_percentiles = [0] + list(np.arange(0.5, 0.99, 0.05)) + [0.99]

            for nb_bin in [3, 5, 10, 20]:
                dext = self.compute_bins(dext, nb_bin)

                for min_max in ["up", "down"]:

                    tmp_l = list()
                    tmp_l_cutoffs = list()
                    for prct in dext_percentiles:
                        prct = round(prct, 2)

                        if min_max == "up":

                            cutoff = round(dext["dext_{}".format(min_max)].quantile(prct), 2)
                            dext_tmp_cutoff = dext.loc[dext["dext_{}".format(min_max)] >= cutoff]

                        elif min_max == "down":
                            cutoff = round(dext["dext_{}_reversed".format(min_max)].quantile(prct), 2)
                            dext_tmp_cutoff = dext.loc[dext["dext_{}_reversed".format(min_max)] >= cutoff]

                        print("Nb bin : {}, Up/Down : {}, Top prct diff : {}, Cutoff associated : {}".format(nb_bin, min_max, prct, cutoff))

                        dext_tmp_cutoff["Exon_bin"] = dext_tmp_cutoff.apply(lambda r: self.compute_exon_bin_fct(r, nb_bin), axis=1)
                        print(dext_tmp_cutoff)
                        tmp_d = {
                            k: v
                            for k, v in dict(
                                collections.Counter(
                                    [",".join([str(sub_e) for sub_e in e]) for e in dext_tmp_cutoff.Exon_bin.values.tolist()]
                                )
                            ).items()
                            if v > nb_bin
                        }
                        print(tmp_d)
                        exit()
                        tmp_d = [[v if int(e) == 1 else 0 for e in k.split(",")] for k, v in tmp_d.items()]
                        tmp_df = pd.DataFrame(tmp_d)
                        tmp_df.loc["Total"] = tmp_df.sum(axis=0)
                        tmp_df.loc["Ratio"] = 100 * (tmp_df.loc["Total"] / tmp_df.loc["Total"].sum())
                        tmp_df.columns = [str(e) for e in range(1, nb_bin + 1, 1)]

                        tmp_df["Cutoff"] = cutoff
                        tmp_df["Nb_bin"] = nb_bin
                        tmp_df["Up/Down"] = min_max
                        tmp_l.append(tmp_df.tail(2))

                    concat_df_distribution = pd.concat(tmp_l, keys=dext_percentiles, axis=0)
                    concat_df_distribution = concat_df_distribution.reset_index()
                    concat_df_distribution.columns = ["Prct", "Total/Ratio"] + list(range(1, nb_bin + 1)) + ["Cutoff", "Nb_bin", "Up/Down"]
                    concat_df_distribution = concat_df_distribution.melt(
                        id_vars=["Prct", "Total/Ratio", "Cutoff", "Nb_bin", "Up/Down"],
                        value_vars=list(range(1, nb_bin + 1)),
                        var_name="Bin_num",
                        value_name="value",
                    )
                    l.append(concat_df_distribution)
            complete_distribution_table = pd.concat(l, axis=0)
            complete_distribution_table.to_excel(output_path, index=False)
        else:
            complete_distribution_table = pd.read_excel(output_path)
        print(complete_distribution_table)
        return complete_distribution_table

    @staticmethod
    def compute_sqtl_bin_fct(r, nb_bin):
        l = [0] * nb_bin
        for j, b in enumerate(r["Gene_bins_{}".format(nb_bin)]):
            if j < nb_bin:
                if b <= r["snpId_POS"] < r["Gene_bins_{}".format(nb_bin)][j + 1]:
                    l[j] += 1
        if r["Strand"] == 0:
            # print("OK")
            # exit()
            l = l[::-1]
        return l

    def compute_all_sqtl_bin(self, output_path, sqtl_path, dext):
        print("# Compute sQTL location for all bin number given in config / File : {}".format(output_path))

        if os.path.isfile(output_path) is True:

            # BINS NUMBER TO DISPLAY DISTRIBUTION
            nb_bins = [3, 5, 10, 20]
            # nb_bins = [3, 5]

            l = list()

            # PRECOMPUTE BIN WINDOWS FOR ALL CONFIGS
            l_bins = list()
            for j, nb_bin in tqdm(enumerate(nb_bins), desc="Compute gene bins"):
                dext_genes = dext[["Gene", "Gene_length", "Gene_start_38"]].drop_duplicates().reset_index(drop=True)
                dext_genes = self.compute_bins(dext_genes, nb_bin)
                if j == 0:
                    l_cols = ["Gene", "Gene_bins"]
                else:
                    l_cols = ["Gene_bins"]
                dext_genes = dext_genes[l_cols].rename({"Gene_bins": "Gene_bins_{}".format(str(nb_bin))}, axis=1)
                l_bins.append(dext_genes)
            # SMALL DF WITH ALL CONFIGS FOR EACH GENE
            dext_genes = pd.concat(l_bins, axis=1)
            print("CONFIGS done, reading sQTL file ...")

            # SQTL FILE
            sqtl = pd.read_parquet(sqtl_path)

            sqtl_lite = sqtl[["snpId", "Strand", "Gene", "Tissue"]]
            print(sqtl_lite)
            sqtl_lite["snpIdTissue"] = sqtl_lite["snpId"] + "_" + sqtl_lite["Tissue"]
            sqtl_lite = sqtl_lite.drop_duplicates(subset=["snpIdTissue"])
            print(sqtl_lite)

            # MERGE WITH PREVIOUS & COMPUTE COLUMNS
            print("Merging bins file with sQTL one ...")
            sqtl_map = pd.merge(sqtl, dext_genes, on="Gene")
            sqtl_map["snpId_POS"] = sqtl_map["snpId"].apply(lambda r: int(r.split("_")[1]))
            print(sqtl_map)

            # ITERATE OVER SCENARIOS
            for nb_bin in nb_bins:
                print("Nb bin : {}".format(nb_bin))

                sqtl_map["snpId_bin_{}".format(nb_bin)] = sqtl_map.parallel_apply(
                    lambda r: self.compute_sqtl_bin_fct(r, nb_bin=nb_bin), axis=1
                )

                # COMPUTE FREQUENCY OF EACH BIN POSSIBILITY
                tmp_d = {
                    k: v
                    for k, v in dict(
                        collections.Counter(
                            [",".join([str(sub_e) for sub_e in e]) for e in sqtl_map["snpId_bin_{}".format(nb_bin)].values.tolist()]
                        )
                    ).items()
                    if v > nb_bin
                }

                # REFORMAT & ADD DF TO LIST
                tmp_d = [[v if int(e) == 1 else 0 for e in k.split(",")] for k, v in tmp_d.items()]
                tmp_df = pd.DataFrame(tmp_d)
                tmp_df.loc["Total"] = tmp_df.sum(axis=0)
                tmp_df.loc["Ratio"] = 100 * (tmp_df.loc["Total"] / tmp_df.loc["Total"].sum())
                tmp_df.columns = [str(e) for e in range(1, nb_bin + 1, 1)]
                tmp_df = tmp_df.tail(2).reset_index()

                tmp_df = tmp_df.melt(
                    id_vars=["index"],
                    value_vars=[str(e) for e in list(range(1, nb_bin + 1))],
                    var_name="Bin_num",
                    value_name="value",
                )

                tmp_df.columns = ["Total/Ratio"] + list(tmp_df.columns)[1:]
                tmp_df["Nb_bin"] = nb_bin
                l.append(tmp_df)

            # CONCAT FINAL DF
            concat_df_distribution = pd.concat(l)
            concat_df_distribution.to_excel(output_path, index=False)
        else:
            concat_df_distribution = pd.read_excel(output_path)
        print(concat_df_distribution)
        return concat_df_distribution

    def compute_sqtl_bin(self, output_path, sqtl_path, dext):
        print("# Compute dsQTL & non-dsQTLs location for all bin number given in config / File : {}".format(output_path))

        if os.path.isfile(output_path) is True:
            dext["dext_down_reversed"] = dext["dext_down"] * -1

            dext_percentiles = list(np.arange(0.5, 0.99, 0.05))
            # dext_percentiles = [0.3, 0.5]
            # dext_percentiles = [0]

            # BINS NUMBER TO DISPLAY DISTRIBUTION
            nb_bins = [3, 5, 10, 20]
            # nb_bins = [3, 5]

            l = list()

            # PRECOMPUTE BIN WINDOWS FOR ALL CONFIGS
            l_bins = list()
            for j, nb_bin in tqdm(enumerate(nb_bins), desc="Compute gene bins"):
                dext_genes = dext[["Gene", "Gene_length", "Gene_start_38"]].drop_duplicates().reset_index(drop=True)
                dext_genes = self.compute_bins(dext_genes, nb_bin)
                if j == 0:
                    l_cols = ["Gene", "Gene_bins"]
                else:
                    l_cols = ["Gene_bins"]
                dext_genes = dext_genes[l_cols].rename({"Gene_bins": "Gene_bins_{}".format(str(nb_bin))}, axis=1)
                l_bins.append(dext_genes)
            # SMALL DF WITH ALL CONFIGS FOR EACH GENE
            dext_genes = pd.concat(l_bins, axis=1)
            print("CONFIGS done, reading sQTL file ...")

            # SQTL FILE
            sqtl = pd.read_parquet(sqtl_path)
            sqtl_lite = sqtl[
                ["snpId", "Strand", "MAP", "Gene", "dext_up", "dext_down_reversed", "dext_tissues_up", "dext_tissues_down", "Tissue"]
            ]
            # print(sqtl_lite)

            sqtl_lite["snpIdTissue"] = sqtl_lite["snpId"] + "_" + sqtl_lite["Tissue"]
            sqtl_lite = sqtl_lite.drop_duplicates(subset=["snpIdTissue", "MAP"])
            # print(sqtl_lite)

            # MERGE WITH PREVIOUS & COMPUTE COLUMNS
            print("Merging bins file with sQTL one ...")

            sqtl_map = pd.merge(sqtl, dext_genes, on="Gene")
            sqtl_map["snpId_POS"] = sqtl_map["snpId"].apply(lambda r: int(r.split("_")[1]))

            # print(sqtl_map)
            # exit()

            # ITERATE OVER SCENARIOS
            for nb_bin in nb_bins:

                sqtl_map["snpId_bin_{}".format(nb_bin)] = sqtl_map.parallel_apply(
                    lambda r: self.compute_sqtl_bin_fct(r, nb_bin=nb_bin), axis=1
                )

                # ITERATE OVER OVER&UNDER EXPRESSED
                for min_max in ["up", "down"]:
                    # for min_max in ["down"]:
                    sqtl_map = sqtl_map.explode("dext_tissues_{}".format(min_max))
                    sqtl_map.loc[sqtl_map["Tissue"] == sqtl_map["dext_tissues_{}".format(min_max)], "Match_tissue"] = True
                    sqtl_map["Match_tissue"] = sqtl_map["Match_tissue"].fillna(False)
                    # print(sqtl_map)
                    # print(sqtl_map.Match_tissue.value_counts())

                    tmp_l = list()

                    # ITERATE OVER MATCHED & UNMATCHEd ON TISSUES
                    for match in [True, False]:
                        # print(match)

                        for prct in dext_percentiles:
                            # print(match)

                            prct = round(prct, 2)

                            if min_max == "up":

                                cutoff = round(dext["dext_{}".format(min_max)].quantile(prct), 2)
                                sqtl_map_tmp_cutoff = sqtl_map.loc[
                                    (sqtl_map["dext_{}".format(min_max)] >= cutoff) & (sqtl_map["Match_tissue"] == match)
                                ]

                            elif min_max == "down":
                                cutoff = round(dext["dext_{}_reversed".format(min_max)].quantile(prct), 2)
                                sqtl_map_tmp_cutoff = sqtl_map.loc[
                                    (sqtl_map["dext_{}_reversed".format(min_max)] >= cutoff) & (sqtl_map["Match_tissue"] == match)
                                ]

                            print(
                                "Nb bin : {}, Up/Down : {}, Top prct diff : {}, Cutoff associated : {}, Tissue match : {}".format(
                                    nb_bin, min_max, prct, cutoff, match
                                )
                            )

                            # print(sqtl_map_tmp_cutoff)

                            # COMPUTE FREQUENCY OF EACH BIN POSSIBILITY
                            tmp_d = {
                                k: v
                                for k, v in dict(
                                    collections.Counter(
                                        [
                                            ",".join([str(sub_e) for sub_e in e])
                                            for e in sqtl_map_tmp_cutoff["snpId_bin_{}".format(nb_bin)].values.tolist()
                                        ]
                                    )
                                ).items()
                                if v > nb_bin
                            }
                            # print(tmp_d)

                            # REFORMAT & ADD DF TO LIST
                            tmp_d = [[v if int(e) == 1 else 0 for e in k.split(",")] for k, v in tmp_d.items()]

                            tmp_df = pd.DataFrame(tmp_d)
                            tmp_df.loc["Total"] = tmp_df.sum(axis=0)
                            tmp_df.loc["Ratio"] = 100 * (tmp_df.loc["Total"] / tmp_df.loc["Total"].sum())
                            # print(tmp_df)
                            tmp_df.columns = [str(e) for e in range(1, nb_bin + 1, 1)]
                            tmp_df["Type"] = match
                            tmp_df["Cutoff"] = cutoff
                            tmp_df["Nb_bin"] = nb_bin
                            tmp_df["Up/Down"] = min_max
                            tmp_df["Prct"] = prct
                            # print(tmp_df)

                            tmp_l.append(tmp_df.tail(2))

                    # pprint(tmp_l)

                    # CONCAT TMP DF

                    concat_df_distribution = pd.concat(tmp_l, axis=0)

                    concat_df_distribution = concat_df_distribution.reset_index()
                    concat_df_distribution.columns = (
                        ["Total/Ratio"]
                        + list(range(1, nb_bin + 1))
                        + [
                            "Type",
                            "Cutoff",
                            "Nb_bin",
                            "Up/Down",
                            "Prct",
                        ]
                    )
                    concat_df_distribution = concat_df_distribution.melt(
                        id_vars=[
                            "Total/Ratio",
                            "Type",
                            "Cutoff",
                            "Nb_bin",
                            "Up/Down",
                            "Prct",
                        ],
                        value_vars=list(range(1, nb_bin + 1)),
                        var_name="Bin_num",
                        value_name="value",
                    )

                    l.append(concat_df_distribution)

            # CONCAT FINAL DF
            concat_df_distribution_table = pd.concat(l, axis=0)
            print(concat_df_distribution_table)
            # exit()
            concat_df_distribution_table.to_excel(output_path, index=False)
        else:
            concat_df_distribution_table = pd.read_excel(output_path)

        # print(sqtl)

        # sqtl.sample(frac=1).head(10000).to_parquet(sqtl_path.replace(".parquet", "_lite.parquet"))


if __name__ == "__main__":
    c = sQTL_location()
