from pprint import pprint
from tqdm import tqdm
import _pickle
import collections
import gzip
import json
import math
import multiprocessing
import numpy as np
import os
import pandas as pd

# import parmap
import re
import requests
import subprocess
import sys
import time
import warnings
import math
from pandarallel import pandarallel

pandarallel.initialize(nb_workers=20, progress_bar=False)

pd.options.mode.chained_assignment = None  # default='warn'
sys.path.append("/home/weber/PycharmProjects/EXOTIC/clean/src")
from utils.utils import load_config_file, mkdir, convert_bins_into_labels


# import hail as hl

# hl.init(min_block_size=128)
# os.environ["PYSPARK_SUBMIT_ARGS"] = "--driver-memory 200 pyspark-shell"


tqdm.pandas()

## YAML FILES CONFIG
yaml = load_config_file(config_file="clean/src/config_clean_clean.yaml")

## JSON DICTS CONFIG
dicts = json.load(open("clean/src/config/EXOTIC_config.json"))


class LocationDev:
    def __init__(self):

        self.cluspack("/home/weber/Documents/test.txt", "kmeans", "dpc")

        refseq = self.refseq_stats_distribution_old(refseq_path=yaml["1_GENOMICS"]["Final"]["refseq_cds_with_variable"])
        refseq = self.process_refseq_to_get_complete_exons(
            path="/gstock/EXOTIC/data/GENOMICS/RefSeq_exons_simple_mrnas.parquet", refseq=refseq
        )
        # self.stats_introns(refseq)
        exit()

        # COMPUTE GENES BINS
        refseq_genes_bins = self.compute_genes_bins(
            nb_bins=20,
            refseq_gene_level_path=yaml["1_GENOMICS"]["TMP"]["tmp_refseq_concat_exons_cds_filtered"],
            output_path=yaml["7_LOCATION"]["TMP"]["refseq_genes_bins"],
        )
        print(refseq_genes_bins)
        exit()

        # COMPUTE EXONS BINS
        refseq_exons_bins = self.compute_exons_bins(
            refseq_exon_level_path=yaml["2_EXPRESSION"]["Final"]["refseq_corrected_cds_recomputed"],
            output_path=yaml["7_LOCATION"]["TMP"]["refseq_exons_bins"],
            refseq_genes_bins=refseq_genes_bins,
        )

        refseq_exons_bins["MAP"] = (
            refseq_exons_bins["Gene"] + "_" + refseq_exons_bins["Exon_start"].astype(str) + "-" + refseq_exons_bins["Exon_stop"].astype(str)
        )
        print(refseq_exons_bins)

        # TMP
        refseq_exons_bins = refseq_exons_bins.drop(["Exons_ranges_y"], axis=1).rename({"Exons_ranges_x": "Exons_ranges"}, axis=1)
        print(refseq_exons_bins)

        # COMPUTE INTRONS RANGES
        refseq_exons_bins = refseq_exons_bins.parallel_apply(self.compute_introns_ranges, axis=1)
        print(refseq_exons_bins)

        exit()

        # # CONST & ALT DISTRI
        # self.process_distribution_stats(
        #     refseq_exons_bins, output_path="/gstock/EXOTIC/data/GENOMICS/TMP/distribution_exons_const_alt.xlsx", conditions=["Const_Alt"]
        # )

        # # PREPARE DEXT
        # dext = pd.read_parquet(yaml["4_DEXT"]["Final"]["dext"])
        # dext.loc[(dext["dext_up"] > 0.5), "dext"] = "UP"
        # dext.loc[(dext["dext_down"] < -0.5), "dext"] = "DOWN"
        # # dext["dext"] = dext["dext"].fillna("None")

        # refseq_exons_bins_dext = pd.merge(refseq_exons_bins, dext[["MAP", "dext"]], on="MAP", how="right")
        # # refseq_exons_bins["dext"] = refseq_exons_bins["dext"].fillna("None")
        # refseq_exons_bins_dext = refseq_exons_bins_dext.dropna(subset=["dext"])

        # # DEXT DISTRI
        # self.process_distribution_stats(
        #     refseq_exons_bins_dext, output_path="/gstock/EXOTIC/data/GENOMICS/TMP/distribution_exons_dext.xlsx", conditions=["dext"]
        # )

        # # PATHOGENIC DISTRI

        # clinvar_data = self.read_and_process_clinvar_vcf(
        #     clinvar_vcf_path=yaml["3_EXONS_PROPERTIES"]["External"]["clinvar_latest"],
        #     clinvar_ht=yaml["3_EXONS_PROPERTIES"]["TMP"]["clinvar_tmp_hail_table"],
        #     output_file=yaml["5_PATHOGENICITY"]["TMP"]["clinvar_full"],
        # )

        # clinvar_data = self.prepare_clinvar_for_distribution(
        #     clinvar=clinvar_data, refseq=refseq_exons_bins, output_path="/gstock/EXOTIC/data/VARIATIONS/clinvar_with_distri_tmp.parquet"
        # )

        # clinvar_data["C"] = ""

        # self.process_distribution_stats(
        #     clinvar_data,
        #     output_path="/gstock/EXOTIC/data/GENOMICS/TMP/distribution_clinvar.xlsx",
        #     conditions=["C"],
        #     cols=["snpId_location_Gene", "snpId_location_Exons"],
        # )

        # self.process_distribution_stats(
        #     clinvar_data,
        #     output_path="/gstock/EXOTIC/data/GENOMICS/TMP/distribution_clinvar_mc.xlsx",
        #     conditions=["MC_lite"],
        #     cols=["snpId_location_Gene", "snpId_location_Exons"],
        # )

        sqtl = pd.read_parquet(yaml["6_sQTLs"]["Final"]["dext_sqtl_lite"])
        print(sqtl)
        exit()

    @staticmethod
    def shift_according_to_strand(r, col="Exons_ranges"):
        if r["Strand"] == 0:
            r[col] = r[col][::-1]
        else:
            pass
        return r

    @staticmethod
    def groupby_apply(tmp_gene_df):
        # print(tmp_gene_df)
        exons_ranges = list()
        previous_start, previous_end = 0, 0

        for j, r in tmp_gene_df.iterrows():
            if r["Share"] is True:
                # print(previous_start, previous_end, r["Start"], r["End"])
                if previous_end == int(r["Start"]):
                    # print("OK", previous_start, previous_end, r["Start"], r["End"])
                    exons_ranges.append("{}-{}".format(previous_start, r["End"]))
            elif int(r["ranges"].split("-")[0]) > previous_end:
                # else:
                exons_ranges.append(r["ranges"])

            previous_start = int(r["Start"])
            previous_end = int(r["End"])
        # print(exons_ranges)
        # exit()
        return exons_ranges

    def process_refseq_to_get_complete_exons(self, path, refseq):

        if os.path.isfile(path) is False:
            refseq = refseq.sort_values(by=["Gene", "Start"])

            # refseq = refseq.loc[refseq["Gene"] == "AASDH"]
            # print(refseq)
            # exit()

            l = list()

            refseq = refseq.explode("mRNA_exons")
            refseq = refseq.sort_values(by=["Gene", "mRNA_exons", "Start"])
            refseq_genes = refseq[["Gene", "mRNA_exons"]].drop_duplicates()

            df = pd.DataFrame(refseq.groupby("mRNA_exons").parallel_apply(lambda r: self.groupby_apply(r))).reset_index()
            print(df)
            df.columns = ["mRNA_exons", "Exons_ranges"]
            df = pd.merge(refseq_genes, df, on="mRNA_exons")

            df = df.parallel_apply(self.compute_introns_ranges, axis=1)

            df = pd.merge(
                refseq[["Gene", "mRNA_nb_total", "Strand", "Miso"]].drop_duplicates(),
                df,
                on="Gene",
            )

            print(df)

            df["CDS_count"] = df["Exons_ranges"].str.len()

            df = df.parallel_apply(lambda r: self.shift_according_to_strand(r, col="Exons_ranges"), axis=1)
            df = df.parallel_apply(lambda r: self.shift_according_to_strand(r, col="Introns_ranges"), axis=1)
            print(df)

            df.to_parquet(path)
        else:
            df = pd.read_parquet(path)
        print(df)

        return df

    @staticmethod
    def cluspack(intput_file, cm, nbc):
        output_file_path = intput_file.split(".")[0] + ".clu"

        if os.path.isfile(output_file_path) is False:
            args = [
                "/biolo/cluspack/binNew/cluspack",
                intput_file,
                "-dt=coordinates",
                "-cm={}".format(cm),
                "-nbc={}".format(nbc),
                "-wc",
                "-oclu",
            ]
            p1 = subprocess.Popen(args)
            p1.wait()

        else:
            d_stats = collections.defaultdict(dict)
            output_file = open(output_file_path, "r").readlines()
            nb_clusters = int(output_file[0].strip().split(" : ")[1])
            cluster = 0
            for line in output_file[1:]:
                line = line.strip()
                if line:
                    if "Cluster" in line:
                        cluster = int(line.split(" ; ")[0].replace("Cluster ", ""))
                        size = int(line.split(" ; ")[1].replace("size=", ""))
                        d_stats[cluster]["Size"] = size
                        d_stats[cluster]["Indexes"] = list()
                    else:
                        d_stats[cluster]["Indexes"].append(int(line))
            cluspack_df = pd.DataFrame.from_dict(d_stats).T.reset_index().rename({"index": "cluster"}, axis=1)
            print(cluspack_df)
            exit()
            return cluspack_df

    @staticmethod
    def refseq_stats_distribution_old(refseq_path):
        refseq = pd.read_parquet(refseq_path)
        refseq = refseq.rename({"Exon_start": "Start", "Exon_stop": "End", "mRNA_gene_nb": "mRNA_nb_total"}, axis=1)
        refseq.loc[refseq["mRNA_nb_total"] > 1, "Miso"] = True
        refseq.loc[refseq["mRNA_nb_total"] == 1, "Miso"] = False
        return refseq

    @staticmethod
    def refseq_stats_distribution_new(refseq_path):
        refseq_raw = pd.read_parquet(refseq_path)
        print(refseq_raw)
        exit()
        refseq = refseq_raw.copy()
        refseq.loc[refseq["mRNA_nb_total"] > 1, "Miso"] = True
        refseq.loc[refseq["mRNA_nb_total"] == 1, "Miso"] = False

        return refseq
        # print(refseq)
        # print(refseq.CDS_count.describe())
        # print(refseq.groupby("Miso").CDS_count.describe())
        # print(refseq.Gene.nunique())

        # exit()

    # @staticmethod
    def stats_introns(self, df_raw):
        def shift_according_to_strand(r, col="Exons_ranges"):
            if r["Strand"] == 0:
                r[col] = r[col][::-1]
            else:
                pass
            return r

        df = df_raw.copy()
        df["CDS_count"] = df["Exons_ranges"].str.len()

        df = df.parallel_apply(lambda r: shift_according_to_strand(r, col="Exons_ranges"), axis=1)
        df = df.parallel_apply(lambda r: shift_according_to_strand(r, col="Introns_ranges"), axis=1)

        df.to_parquet("/gstock/EXOTIC/data/GENOMICS/raw_matrix_clustering.parquet")

        # key = 5
        # # key_two = 20
        # df = df.loc[(df["CDS_count"] > key)]
        # print(df)
        # df["Introns_lengths"] = df["Introns_ranges"].apply(lambda r: [int(e.split("-")[1]) - int(e.split("-")[0]) for e in r])

        # df["Introns_lengths_pct"] = df["Introns_lengths"].apply(lambda r: pd.Series(r).rank(pct=True).round(2).values.tolist())

        # df["Introns_lengths_ratio"] = df["Introns_lengths"].apply(
        #     lambda r: pd.Series(pd.Series(r) / pd.Series(r).sum()).round(2).values.tolist()
        # )

        # print(df)
        exit()

        df["Introns_lengths_ratio"] = df["Introns_lengths_ratio"].apply(lambda r: r[:3] + r[-2:])

        df = df.reset_index(drop=True)
        X = pd.DataFrame.from_records(df["Introns_lengths_ratio"])

        print(X)

        from sklearn.cluster import MiniBatchKMeans, KMeans
        from sklearn.decomposition import PCA, IncrementalPCA

        # km = KMeans(n_clusters=3)
        # print(km)
        # km.fit(X.values)
        # cluster_centers = km.cluster_centers_
        # cluster_centers = [[round(sub_e, 2) for sub_e in e] for e in cluster_centers]

        # print(cluster_centers)

        # cluster_map = pd.DataFrame()
        # cluster_map["data_index"] = X.index.values
        # cluster_map["cluster"] = km.labels_
        # print(cluster_map.cluster.value_counts().sort_index())
        # print(cluster_map)
        # print(df.loc[cluster_map.loc[cluster_map["cluster"] == 4].index.tolist()])

        n_components = 2
        pca = PCA(n_components=n_components)
        X_pca = pca.fit_transform(X)

        exit()

    @staticmethod
    def compute_introns_ranges(r):
        l = list()
        # print(r)
        exons = r["Exons_ranges"]
        # print(exons)
        exons_start = list(sorted([e.split("-")[0] for e in exons]))
        exons_end = list(sorted([e.split("-")[1] for e in exons]))
        exons = ["{}-{}".format(s, e) for s, e in zip(exons_start, exons_end)]

        # print(exons)
        # exons = list(sorted([int(sub_e) for e in exons for sub_e in e.split("-")]))
        # print(exons)
        # exons = ["{}-{}".format(e, exons[j + 1]) for j, e in enumerate(exons) if j < len(exons) - 1 if j % 2 == 0]
        # print(exons)
        # print(len(exons))

        for j, e in enumerate(exons):
            if j == 0:
                l.append(int(e.split("-")[1]) + 1)
            elif j > 0 and j < len(exons) - 1:
                l.append(int(e.split("-")[0]) - 1)
                l.append(int(e.split("-")[1]) + 1)
            elif j == len(exons) - 1:
                l.append(int(e.split("-")[0]) - 1)
        # print(l)
        l = ["{}-{}".format(e, l[j + 1]) for j, e in enumerate(l) if j < len(l) - 1 if j % 2 == 0]
        # print(l)
        # exit()
        # introns_start = list(sorted([e.split("-")[0] for e in l]))
        # introns_end = list(sorted([e.split("-")[1] for e in l]))
        # introns = ["{}-{}".format(s, e) for s, e in zip(introns_start, introns_end)]
        r["Introns_ranges"] = l
        # print(l)
        # print(len(l))
        return r

    def compute_sqtl_bin(self, output_path, sqtl_path, dext):
        print("# Compute dsQTL & non-dsQTLs location for all bin number given in config / File : {}".format(output_path))

        if os.path.isfile(output_path) is True:

            dext["dext_down_reversed"] = dext["dext_down"] * -1

            nb_bin = 20

            l = list()

            if os.path.isfile(sqtl_map_path) is False:
                # BINS NUMBER TO DISPLAY DISTRIBUTION

                dext_genes = dext[["Gene", "Gene_length", "Gene_start_38"]].drop_duplicates().reset_index(drop=True)
                dext_genes = self.compute_bins(dext_genes, nb_bin)

                l_cols = ["Gene", "Gene_bins"]

                dext_genes = dext_genes[l_cols].rename({"Gene_bins": "Gene_bins_{}".format(str(nb_bin))}, axis=1)
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

                sqtl_map["snpId_bin_{}".format(nb_bin)] = sqtl_map.parallel_apply(
                    lambda r: self.compute_sqtl_bin_fct(r, nb_bin=nb_bin), axis=1
                )

                sqtl_map.to_parquet(sqtl_map_path)

                # exit()
            else:
                sqtl_map = pd.read_parquet(sqtl_map_path)

            print(sqtl_map)

            # ITERATE OVER SCENARIOS
            for nb_bin in nb_bins:

                # ITERATE OVER OVER&UNDER EXPRESSED
                for min_max in ["up", "down"]:
                    # for min_max in ["down"]:

                    tmp_l = list()

                    # ITERATE OVER MATCHED & UNMATCHEd ON TISSUES
                    for match in [True, False]:
                        # print(match)

                        for prct in dext_percentiles:

                            prct = round(prct, 2)

                            if min_max == "up":

                                # cutoff = round(dext["dext_{}".format(min_max)].quantile(prct), 2)
                                cutoff = 0.5

                                # diff_spec = list(
                                # set(dext.loc[dext["dext_{}".format(min_max)] >= cutoff, "Gene"].unique().tolist()).difference(
                                # set(dext.loc[dext["dext_{}".format(min_max)] < cutoff, "Gene"].unique().tolist())
                                # )
                                # )

                                sqtl_map_tmp_cutoff = sqtl_map.loc[
                                    (sqtl_map["dext_{}".format(min_max)] >= cutoff)
                                    # & (sqtl_map["Match_tissue"] == match)
                                    # & (sqtl_map["Gene"].isin(diff_spec))
                                ]
                                print(sqtl_map_tmp_cutoff)

                                sqtl_map_tmp_cutoff = sqtl_map_tmp_cutoff.explode("dext_tissues_{}".format(min_max))
                                print(sqtl_map_tmp_cutoff)
                                sqtl_map_tmp_cutoff.loc[
                                    sqtl_map_tmp_cutoff["Tissue"] == sqtl_map_tmp_cutoff["dext_tissues_{}".format(min_max)], "Match_tissue"
                                ] = True
                                sqtl_map_tmp_cutoff["Match_tissue"] = sqtl_map_tmp_cutoff["Match_tissue"].fillna(False)

                                sqtl_map_tmp_cutoff = sqtl_map_tmp_cutoff.loc[
                                    # (sqtl_map_tmp_cutoff["dext_{}".format(min_max)] >= cutoff)
                                    (sqtl_map_tmp_cutoff["Match_tissue"] == match)
                                    # & (sqtl_map["Gene"].isin(diff_spec))
                                ]
                                print(sqtl_map_tmp_cutoff)
                                # dext_tmp_cutoff = dext.loc[(dext["dext_{}".format(min_max)] >= cutoff) & (dext["Gene"].isin(diff_spec))]

                            elif min_max == "down":

                                # cutoff = round(dext["dext_{}".format(min_max)].quantile(prct), 2)
                                # cutoff = 0.5

                                # diff_spec = list(
                                # set(dext.loc[dext["dext_{}_reversed".format(min_max)] >= cutoff, "Gene"].unique().tolist()).difference(
                                # set(dext.loc[dext["dext_{}_reversed".format(min_max)] < cutoff, "Gene"].unique().tolist())
                                # )
                                # )

                                sqtl_map_tmp_cutoff = sqtl_map.loc[
                                    (sqtl_map["dext_{}_reversed".format(min_max)] >= cutoff)
                                    # & (sqtl_map["Match_tissue"] == match)
                                    # & (sqtl_map["Gene"].isin(diff_spec))
                                ]

                                sqtl_map_tmp_cutoff = sqtl_map_tmp_cutoff.explode("dext_tissues_{}".format(min_max))

                                sqtl_map_tmp_cutoff.loc[
                                    sqtl_map_tmp_cutoff["Tissue"] == sqtl_map_tmp_cutoff["dext_tissues_{}".format(min_max)], "Match_tissue"
                                ] = True
                                sqtl_map_tmp_cutoff["Match_tissue"] = sqtl_map_tmp_cutoff["Match_tissue"].fillna(False)

                                sqtl_map_tmp_cutoff = sqtl_map_tmp_cutoff.loc[
                                    # (sqtl_map_tmp_cutoff["dext_{}".format(min_max)] >= cutoff)
                                    (sqtl_map_tmp_cutoff["Match_tissue"] == match)
                                    # & (sqtl_map["Gene"].isin(diff_spec))
                                ]
                                # dext_tmp_cutoff = dext.loc[(dext["dext_{}".format(min_max)] >= cutoff) & (dext["Gene"].isin(diff_spec))]

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
                                # if v > nb_bin
                            }
                            print(tmp_d)

                            # REFORMAT & ADD DF TO LIST
                            tmp_d = [[v if int(e) == 1 else 0 for e in k.split(",")] for k, v in tmp_d.items()]

                            tmp_df = pd.DataFrame(tmp_d)
                            print(tmp_df)
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

    @staticmethod
    def read_and_process_clinvar_vcf(clinvar_vcf_path, clinvar_ht, output_file):
        print("# Load ClinVar VCF, build ht + pandas & retrieve corresponding variations / File : {}".format(output_file))

        if os.path.isfile(output_file) is False:

            # HAIL PART

            # hl.import_vcf(clinvar_vcf_path, force_bgz=True).write(clinvar_vcf_path.replace(".vcf.gz", ".ht"))

            clinvar_lite = hl.read_matrix_table("/gstock/biolo_datasets/variation/variation_sets/clinvar/vcf_GRCh37/v2/clinvar_20210123.ht")

            clinvar_lite = clinvar_lite.filter_rows(clinvar_lite.info.CLNVC == "single_nucleotide_variant")
            clinvar_lite = clinvar_lite.select_rows(
                clinvar_lite.rsid,
                clinvar_lite.info.ALLELEID,
                clinvar_lite.info.CLNREVSTAT,
                clinvar_lite.info.CLNSIG,
                clinvar_lite.info.GENEINFO,
                clinvar_lite.info.RS,
                clinvar_lite.info.CLNVI,
                clinvar_lite.info.MC,
            )

            # PANDAS PART

            clinvar_lite_pd_raw = clinvar_lite.make_table().to_pandas()
            clinvar_lite_pd = clinvar_lite_pd_raw.copy()
            clinvar_lite_pd = clinvar_lite_pd.replace(to_replace="None", value=np.nan)

            clinvar_lite_pd = clinvar_lite_pd.dropna(subset=["ALLELEID", "CLNREVSTAT", "CLNSIG", "GENEINFO"], how="any")
            clinvar_lite_pd["CLNREVSTAT"] = clinvar_lite_pd["CLNREVSTAT"].apply(lambda r: ",".join(r))
            clinvar_lite_pd["CLNSIG"] = clinvar_lite_pd["CLNSIG"].apply(lambda r: ",".join(r))
            clinvar_lite_pd["RS_STARS"] = clinvar_lite_pd["CLNREVSTAT"].map(dicts["clinvar_review_status_dict"])
            clinvar_lite_pd = clinvar_lite_pd.loc[
                (clinvar_lite_pd["CLNSIG"].str.contains("athogenic"))
                & (~clinvar_lite_pd["CLNSIG"].str.contains("Conflicting_interpretations_of_pathogenicity"))
            ]
            clinvar_lite_pd = clinvar_lite_pd.loc[~clinvar_lite_pd["CLNREVSTAT"].str.contains("interpret")]

            clinvar_lite_pd.to_parquet(output_file)
        else:
            clinvar_lite_pd = pd.read_parquet(output_file)

        return clinvar_lite_pd

    def prepare_clinvar_for_distribution(self, clinvar, refseq, output_path):

        if os.path.isfile(output_path) is False:

            # CLEAN & FILTER
            clinvar_data = clinvar.dropna(subset=["MC"])
            clinvar_data["Gene"] = clinvar_data["GENEINFO"].apply(lambda r: r.split(":")[0])
            clinvar_data["MC_lite"] = clinvar_data["MC"].apply(lambda r: r[0])
            clinvar_data["MC_lite"] = clinvar_data["MC_lite"].apply(lambda r: r.split("|")[1])
            clinvar_data = clinvar_data.loc[clinvar_data["alleles"].str.len() == 2]

            # FILTER ON MISSENSE
            # clinvar_data = clinvar_data.loc[clinvar_data["MC_lite"].str.contains("missense")]

            # MERGE WITH STRAND &
            clinvar_data = pd.merge(
                refseq[["Gene", "Exons_bins", "Gene_bins", "Strand"]].drop_duplicates(subset=["Gene"]), clinvar_data, on=["Gene"]
            )
            clinvar_data = clinvar_data.parallel_apply(lambda r: self.compute_snv_bin_fct(r, col_cat="Exons_bins"), axis=1)
            clinvar_data = clinvar_data.parallel_apply(lambda r: self.compute_snv_bin_fct(r, col_cat="Gene_bins"), axis=1)
            clinvar_data.to_parquet(output_path)

        else:
            clinvar_data = pd.read_parquet(output_path)
        print(clinvar_data)

        return clinvar_data

    @staticmethod
    def distribution_stats(df, col="Exon_location_Gene"):

        df = df.dropna(subset=[col])
        df_sum = pd.DataFrame.from_records(df[col].values).sum()
        df_sum.index = [e + 1 for e in list(df_sum.index)]
        df_sum = pd.DataFrame(df_sum).reset_index()
        df_sum.columns = ["Bin_num", "Total"]
        df_sum["Ratio"] = 100 * (df_sum["Total"] / df_sum["Total"].sum())
        # print(df_sum)
        return df_sum
        # print(df_sum)
        # exit()
        # df_sum.to_excel("/gstock/EXOTIC/data/VARIATIONS/pathogenic_distribution_missense_exons.xlsx", index=False)

    def process_distribution_stats(self, df, output_path, cols=["Exon_location_Gene", "Exon_location_Exons"], conditions=["Const_Alt"]):

        if os.path.isfile(output_path) is False:

            l = list()
            import itertools

            d = {c: df[c].unique().tolist() for c in conditions}
            combinations = list(itertools.product(*d.values()))
            for combi in combinations:
                tmp_df = df.copy()
                tmp_l = list()
                for col, condition in zip(conditions, combi):
                    tmp_df = tmp_df.loc[tmp_df[col] == condition]
                    tmp_l.append(condition)
                    if col == conditions[-1]:

                        if tmp_df.empty is False:
                            tmp_stats_gene = self.distribution_stats(tmp_df, col=cols[0])
                            tmp_stats_gene["Type"] = "Gene"
                            tmp_stats_gene["Filters"] = "-".join(tmp_l)
                            tmp_stats_exons = self.distribution_stats(tmp_df, col=cols[1])
                            tmp_stats_exons["Type"] = "Exons"
                            tmp_stats_exons["Filters"] = "-".join(tmp_l)
                            l.append(pd.concat([tmp_stats_gene, tmp_stats_exons]))
            concat_final = pd.concat(l)
            if len(conditions) > 1:
                concat_final[conditions] = concat_final["Filters"].str.split("-", expand=True)
            concat_final.to_excel(output_path, index=False)

        else:
            concat_final = pd.read_excel(output_path)
        print(concat_final)

    @staticmethod
    def make_windows(r, nb_bins):
        absolute = r["Exon_all_start"]
        relative = 0
        window = 1

        k = math.ceil(r["Exon_window_bp"] / nb_bins)

        bin_ranges = list()

        for e in r["Exons_ranges"]:
            start = int(e.split("-")[0])
            end = int(e.split("-")[1]) + 1
            for i in range(start, end):
                relative += 1
                if relative == 1:
                    bin_ranges.append(i)
                elif relative % k == 0 and window < nb_bins:
                    window += 1
                    bin_ranges.append(i)
        bin_ranges.append(end - 1)
        # bin_ranges = ["{}-{}".format(e, bin_ranges[j + 1]) for j, e in enumerate(bin_ranges) if j < len(bin_ranges) - 1]
        r["Exons_bins"] = bin_ranges

        return r

    def compute_exons_bins(
        self,
        refseq_exon_level_path,
        output_path,
        refseq_genes_bins,
        nb_bins=20,
    ):

        if os.path.isfile(output_path) is True:

            refseq_raw = pd.read_parquet(refseq_exon_level_path)
            refseq = refseq_raw.groupby("Gene")["ranges"].apply(list).reset_index()
            refseq["ranges"] = refseq["ranges"].apply(lambda r: list(sorted(r)))
            refseq["Exon_window_bp"] = refseq["ranges"].apply(lambda r: sum([int(e.split("-")[1]) - int(e.split("-")[0]) for e in r]))
            refseq["Exon_all_start"] = refseq.ranges.apply(lambda r: r[0].split("-")[0])
            refseq["Exon_all_end"] = refseq.ranges.apply(lambda r: r[-1].split("-")[1])
            refseq[["Exon_all_start", "Exon_all_end"]] = refseq[["Exon_all_start", "Exon_all_end"]].astype(int)
            refseq = refseq.rename({"ranges": "Exons_ranges", "Start": "Exon_start", "End": "Exon_stop"}, axis=1)
            refseq_raw = refseq_raw.rename({"Start": "Exon_start", "End": "Exon_stop"}, axis=1)

            refseq = refseq.parallel_apply(lambda r: self.make_windows(r, nb_bins), axis=1)
            refseq = pd.merge(refseq, refseq_genes_bins, on="Gene")
            refseq = pd.merge(refseq, refseq_raw.drop(["Strand"], axis=1), on="Gene")

            refseq[["Exon_start", "Exon_stop"]] = refseq[["Exon_start", "Exon_stop"]].astype(int)

            refseq = refseq.parallel_apply(lambda r: self.compute_exon_bin_fct(r, nb_bins, col="Gene_bins"), axis=1)
            refseq = refseq.parallel_apply(lambda r: self.compute_exon_bin_fct(r, nb_bins, col="Exons_bins"), axis=1)
            refseq.to_parquet(output_path)
        else:
            refseq = pd.read_parquet(output_path)
        print(refseq)
        return refseq

    # @staticmethod
    def compute_snv_bin_fct(self, r, nb_bin=20, col_snv="locus.position", col_cat="Exons_bins"):
        l = [0] * nb_bin
        for j, b in enumerate(r[col_cat]):
            if j < nb_bin:
                if b <= r[col_snv] < r[col_cat][j + 1]:
                    l[j] += 1
        if r["Strand"] == 0:
            l = l[::-1]
        else:
            l = l
        r["snpId_location_{}".format(col_cat.split("_")[0])] = l
        return r

    @staticmethod
    def compute_genes_bins(refseq_gene_level_path, output_path, nb_bins=10):
        if os.path.isfile(output_path) is True:
            genes_remap_raw = pd.read_parquet(refseq_gene_level_path).sort_values(by=["Gene", "ranges"])
            genes_remap_raw = genes_remap_raw.loc[genes_remap_raw["Exon_type"] == "CDS"].drop_duplicates(subset=["Gene", "ranges"])

            genes_remap = genes_remap_raw[["Gene", "Gene_start", "Gene_stop", "Strand"]].drop_duplicates().sort_values(by="Gene")

            genes_remap["Strand"] = genes_remap["Strand"].replace({"-": 0, "+": 1})

            genes_remap["Gene_length"] = genes_remap["Gene_stop"] - genes_remap["Gene_start"]
            genes_remap["Gene_bin_size"] = genes_remap["Gene_length"] / nb_bins

            genes_remap["Gene_bins"] = genes_remap.apply(
                lambda r: [r["Gene_start"] + (round(r["Gene_bin_size"] * j)) for j in list(range(nb_bins + 1))], axis=1
            )
            genes_remap.to_parquet(output_path)
        else:
            genes_remap = pd.read_parquet(output_path)

        return genes_remap

    @staticmethod
    def compute_exon_bin_fct(r, nb_bin, col="Gene_bins"):

        l = [0] * nb_bin
        for j, b in enumerate(r[col]):
            if j < nb_bin:
                if r["Exon_start"] >= b and r["Exon_stop"] <= r[col][j + 1]:
                    l[j] += 1

        if l == [0] * nb_bin:
            for j, b in enumerate(r[col]):
                if j < nb_bin:
                    if r["Exon_start"] < b < r["Exon_stop"]:
                        l[j - 1] += 1
                        l[j] += 1

        if r["Strand"] == -1:
            l = l[::-1]

        l = [1 if e else 0 for e in l]
        r["Exon_location_{}".format(col.split("_")[0])] = l
        return r


if __name__ == "__main__":
    c = LocationDev()