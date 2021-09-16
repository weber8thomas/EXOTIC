# General imports
import os
import sys
from scipy.stats import fisher_exact
import pandas as pd

pd.options.mode.chained_assignment = None  # default='warn'
import numpy as np
import subprocess

from sklearn.cluster import KMeans

from tqdm import tqdm

# Other imports
tqdm.pandas()

import json

from pprint import pprint

# Custom utils
sys.path.append("/home/weber/PycharmProjects/EXOTIC/src")
from utils.utils import load_config_file

# Figures imports
import matplotlib

from matplotlib.lines import Line2D
import matplotlib.patches as mpatches

import random

import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
from statannot import add_stat_annotation
import matplotlib.font_manager as font_manager

plt.rcParams.update({"font.size": 30})

# Font settings
font_dirs = [
    "/home/weber/Fonts",
]
font_files = font_manager.findSystemFonts(fontpaths=font_dirs)
font_list = font_manager.createFontList(font_files)
font_manager.fontManager.ttflist.extend(font_list)

from matplotlib import rcParams

rcParams["font.family"] = "sans-serif"
rcParams["font.sans-serif"] = ["Arial"]
rcParams["font.weight"] = "light"

## YAML FILES CONFIG
yaml = load_config_file(config_file="/home/weber/PycharmProjects/EXOTIC/clean/src/config_clean_clean.yaml")

dicts = json.load(open("/home/weber/PycharmProjects/EXOTIC/src/EXOTIC_config.json"))

random.seed(1)
colors = sns.color_palette("Paired")
random.shuffle(colors)

conditions = ["Miso", "dext"]


class IntronStats:
    def __init__(self):
        df_raw_raw = self.prepare_file()

        import multiprocessing
        import parmap

        # conditions = ["Miso", "dext", "dext_down", "dext_up"]
        conditions = ["Miso"]
        keys = list(range(3, 31))
        # keys = list(range(3, 31))[:5]
        # signs = [">=", "=="]
        signs = ["=="]

        m = multiprocessing.Manager()

        # l_clusters_df = list()
        # l_fold_df = list()
        cluster_nb = 4
        # seeds = list(range(0, 10))
        seeds = [2]

        self.convert_sign = {">=": "supequal", "==": "equal"}

        # for sign in signs:
        #     for condition in conditions:
        #         df_raw = self.prepare_file_based_on_condition(df_raw, condition)
        #         for intron in tqdm(keys):
        #             X, df, d = self.launch_intron(df_raw, intron, condition, d, sign=sign)
        #             cluster_centers, gb = self.clustering(X, df, cluster_nb, condition, df_raw)
        #             l_clusters_df.append(cluster_centers)
        #             l_fold_df.append(gb)

        for seed in seeds:
            for sign in signs:
                l_clusters_df = m.list()
                l_fold_df = m.list()
                l_df = m.list()
                l_gene_df = m.list()
                for condition in conditions:
                    df_raw = self.prepare_file_based_on_condition(df_raw_raw, condition)
                    print(sign, condition)
                    # print(df_raw_raw[condition].value_counts())
                    # print(df_raw[condition].value_counts())
                    # for intron in tqdm(keys):
                    parmap.starmap(
                        self.test,
                        list(zip(keys)),
                        # intron,
                        df_raw,
                        condition,
                        sign,
                        cluster_nb,
                        l_clusters_df,
                        l_fold_df,
                        l_df,
                        # l_gene_df,
                        df_raw_raw,
                        seed,
                        pm_pbar=True,
                    )
                l_clusters_df = list(l_clusters_df)
                l_fold_df = list(l_fold_df)
                # print(l_clusters_df)
                l_df = list(l_df)

                concat_density, index_df = self.concat_mass_centers(
                    l_clusters_df, conditions, "/gstock/EXOTIC/data/CLUSTERING/tmp.xlsx", keys, cluster_nb
                )
                concat_stats = self.concat_stats(l_fold_df, conditions, "/gstock/EXOTIC/data/CLUSTERING/tmp2.xlsx", keys, cluster_nb)
                map_d = self.dict_new_clusters(index_df, keys, cluster_nb)
                # pprint(map_d)

                concat = pd.concat(l_df)
                concat["cluster_new"] = concat["id"].map(map_d)
                # concat["cluster_new"] = concat["cluster_new"].astype(int)
                concat["Condition"] = concat["id"].apply(lambda r: r.split("-")[0])
                tmp_data_genes_clustering = concat[["Introns_nb", "Miso", "cluster_new", "Gene"]].sort_values(
                    by=["Introns_nb", "Miso", "cluster_new", "Gene"]
                )
                print(tmp_data_genes_clustering)
                tmp_data_genes_clustering.to_excel(
                    "/gstock/EXOTIC/data/CLUSTERING/data_genes_clustering_{}.xlsx".format(self.convert_sign[sign], index=False)
                )
                # exit()
                # print(concat)
                tmp_concat_l = list()
                for condition in concat.Condition.unique():
                    tmp_concat_l.append(concat.loc[(concat["Condition"] == condition) & (concat[condition] == "True")])
                concat = pd.concat(tmp_concat_l)

                # print(concat)
                concat_gb = concat.groupby(["Condition", "cluster_new", "Introns_nb"])["Gene"].apply(list)
                # print(concat_gb)

                data_barplot = self.data_barplot(concat_stats, keys, cluster_nb).reset_index()
                data_profiles = self.data_profiles(map_d, concat_density, keys, cluster_nb)
                # print(data_profiles)
                # exit()
                data_heatmap = self.fold_reordered(map_d, concat_stats, keys, cluster_nb)

                data_barplot.to_excel("/gstock/EXOTIC/data/CLUSTERING/data_barplot_{}.xlsx".format(self.convert_sign[sign]), index=False)
                data_profiles.to_excel("/gstock/EXOTIC/data/CLUSTERING/data_profiles_{}.xlsx".format(self.convert_sign[sign]), index=False)
                data_heatmap.reset_index().to_excel(
                    "/gstock/EXOTIC/data/CLUSTERING/data_heatmap_{}.xlsx".format(self.convert_sign[sign]), index=False
                )
                concat_gb.to_excel("/gstock/EXOTIC/data/CLUSTERING/data_genes_{}.xlsx".format(self.convert_sign[sign]), index=True)
                # print(data_profiles)
                # print(data_heatmap)
                # print(concat_density)
                # print(data_heatmap.fillna(100).reset_index().drop(["Introns_nb"], axis=1).groupby("Condition").describe().T)
                self.launch_heatmap(df_raw, data_heatmap, data_barplot, data_profiles, concat_stats, conditions, keys, sign, seed)

    def test(self, intron, df_raw, condition, sign, cluster_nb, l_clusters_df, l_fold_df, l_df, df_raw_raw, seed):
        X, df = self.launch_intron(df_raw, intron, condition, sign=sign)
        cluster_centers, gb, df = self.clustering(X, df, cluster_nb, intron, condition, df_raw, df_raw_raw, seed)
        # print(intron, df.shape, df[condition].value_counts().to_dict(), gb, cluster_centers)

        l_clusters_df.append(cluster_centers)
        l_fold_df.append(gb)
        l_df.append(df)

    @staticmethod
    def prepare_file():

        df_raw_raw = pd.read_parquet("/gstock/EXOTIC/data/GENOMICS/RefSeq_exons_simple.parquet")
        df_raw = df_raw_raw.copy()
        # print(df_raw)
        # exit()
        dext = pd.read_parquet("/gstock/EXOTIC/data/EXOTIC/dext_matrix.parquet")

        df_raw.loc[df_raw["Gene"].isin(dext[(dext["dext_up"] > 0.5) | (dext["dext_down"] < -0.5)].Gene.unique().tolist()), "dext"] = True
        df_raw.loc[df_raw["Gene"].isin(dext[(dext["dext_up"] > 0.5)].Gene.unique().tolist()), "dext_up"] = True
        df_raw.loc[df_raw["Gene"].isin(dext[(dext["dext_down"] < -0.5)].Gene.unique().tolist()), "dext_down"] = True
        df_raw[["dext", "dext_up", "dext_down"]] = df_raw[["dext", "dext_up", "dext_down"]].fillna(False)

        df_raw["Introns_lengths"] = df_raw["Introns_ranges"].apply(
            lambda r: [
                int(e.split("-")[1]) - int(e.split("-")[0])
                for e in r
                if int(e.split("-")[1]) - int(e.split("-")[0]) > 0 and int(e.split("-")[1]) - int(e.split("-")[0]) < 2e6
            ]
        )
        df_raw["Introns_nb"] = df_raw["Introns_lengths"].str.len()

        return df_raw

    @staticmethod
    def prepare_file_based_on_condition(df, condition):
        anti_condition = "dext_up" if condition == "dext_down" else "dext_down"
        if "dext" in condition:
            df = df.loc[df["Miso"] == True]
        if "up" in condition or "down" in condition:
            df = df.loc[df[anti_condition] == False]
        return df

    @staticmethod
    def launch_intron(df_raw, k, condition, sign="=="):
        # FILTER GENES WITH == OR >=
        if sign == "==":
            df = df_raw.loc[(df_raw["Introns_nb"] == k)]
        elif sign == ">=":
            df = df_raw.loc[(df_raw["Introns_nb"] >= k)]

        # print(k, df.shape, df[condition].value_counts().to_dict())

        # CONVERT TO RATIO & SELECT FIRST ONES
        df["Introns_lengths_ratio"] = df["Introns_lengths"].apply(
            lambda r: pd.Series(pd.Series(r) / pd.Series(r).sum()).round(2).values.tolist()
        )
        # df["Introns_lengths_ratio"] = df["Introns_lengths"].apply(lambda r: pd.Series(r).values.tolist())

        df["Introns_lengths_ratio"] = df["Introns_lengths_ratio"].apply(lambda r: r[:k])
        # print(df)

        # REINDEX & COUNT
        df = df.reset_index(drop=True)

        x = "Introns_lengths_ratio"
        X_raw = pd.DataFrame.from_records(df[x])
        return X_raw, df

    @staticmethod
    def clustering(X_raw, df, cluster, intron, condition, df_raw, df_raw_raw, seed):
        # CLUSTERING
        X = X_raw.copy()
        km = KMeans(n_clusters=cluster, random_state=seed)
        km.fit(X.values)
        cluster_centers = km.cluster_centers_
        cluster_centers = [[round(sub_e, 2) for sub_e in e] for e in cluster_centers]
        cluster_centers = pd.DataFrame.from_records(cluster_centers)
        cluster_centers["Condition"] = condition
        cluster_centers["Introns_nb"] = intron

        # l_clusters_df.append(cluster_centers)

        cluster_map = pd.DataFrame()
        cluster_map["data_index"] = X.index.values
        cluster_map["cluster"] = km.labels_
        X["cluster"] = cluster_map.cluster
        X["Miso"] = df["Miso"]
        X[condition] = df[condition]
        X[condition] = X[condition].astype(str)
        X["Strand"] = df["Strand"]
        X["Gene"] = df["Gene"]
        X["Introns_ranges"] = df["Introns_ranges"]
        X["Introns_lengths"] = df["Introns_lengths"]
        X["Introns_nb"] = X["Introns_lengths"].apply(len)
        X["id"] = str(condition) + "-" + X["Introns_nb"].astype(str) + "-" + X["cluster"].astype(str)

        gb_raw = (
            X.groupby(["cluster"])[condition]
            .value_counts()
            .rename("Count")
            .reset_index()
            .pivot(index="cluster", columns=condition, values="Count")
        )
        gb_raw.loc["Sum"] = gb_raw.sum()

        gb = gb_raw.copy()
        gb.columns = [str(c) for c in gb.columns]
        gb = gb.rename({"False": "False_raw", "True": "True_raw"}, axis=1)

        gb["False"] = 100 * (gb["False_raw"] / df_raw[condition].value_counts().loc[False])
        gb["True"] = 100 * (gb["True_raw"] / df_raw[condition].value_counts().loc[True])

        # gb["Ratio"] = gb["True"] / gb["False"]
        gb["Ratio"] = gb["True_raw"] / gb["False_raw"]
        df_raw_condition_count = df_raw[condition].value_counts().to_dict()
        ref_ratio = df_raw_condition_count[True] / df_raw_condition_count[False]
        gb["Ref_Ratio"] = ref_ratio

        # gb["Fold"] = gb["Ratio"] / gb["Ref_Ratio"]
        gb["Fold"] = gb["Ratio"] / gb.loc["Sum"]["Ratio"]
        gb["Fold"] = gb["Fold"] - 1

        gb["Total"] = gb["False"] + gb["True"]

        # l_fold_df.append(gb)

        gb["Introns_nb"] = intron
        gb["Condition"] = condition

        gb = gb.drop("Sum").reset_index(drop=True)

        # print(condition)
        # print(gb)

        return cluster_centers, gb, X

    @staticmethod
    def changing_index(concat):
        # Index order
        def return_index(r):
            r = r.drop(["Introns_nb", "Condition"], axis=1).dropna(axis=1)
            return r.apply(lambda r: [i for i, e in enumerate(r) if e == max(r)][0])

        index_df = (
            concat.groupby("Introns_nb")
            .apply(return_index)
            .reset_index()
            .pivot(index="Introns_nb", columns="level_1", values=0)[[0, 1, 2]]
            .astype(int)
        )
        index_df[3] = index_df.apply(lambda r: [i for i in range(0, 4) if i not in r.values.tolist()][0], axis=1)
        return index_df

    def concat_mass_centers(self, l_clusters_df, conditions, path, keys, cluster_nb):

        import collections

        def correct_clusters(r):
            baseline = list(range(0, cluster_nb))
            if sum(r) != sum(list(range(0, cluster_nb))):
                counter = collections.Counter(r.values.tolist())
                uncorrects = [k for k, v in counter.items() if v != 1]
                new_values = list(sorted([e for e in baseline if e not in r.values.tolist()] + uncorrects))
                r_raw = [e if e not in uncorrects else -1 for e in r.values.tolist()]
                new_r = list()
                i = 0
                for e in r_raw:
                    if e >= 0:
                        new_r.append(e)
                    else:
                        new_r.append(new_values[i])
                        i += 1
                return new_r
            else:
                return r

        if os.path.isfile(path) is True:

            # Mass centers
            l = list()
            j = 0
            # for i, e in enumerate(l_clusters_df):
            #     tmp_df_to_add = pd.DataFrame.from_records(e)
            #     l.append(tmp_df_to_add)
            concat = pd.concat(l_clusters_df)
            # concat["Condition"] = [sub_e for e in conditions for sub_e in [e] * (len(keys) * cluster_nb)]
            concat.index.names = ["Cluster"]
            concat = concat.reset_index()
            concat = concat.sort_values(by=["Condition", "Introns_nb", "Cluster"])
            # print(concat)
            # concat["Introns_nb"] = len(conditions) * [sub_i for i in keys for sub_i in [i] * cluster_nb]
            # concat = concat[["Condition", "Introns_nb"] + list(range(0, 30))].reset_index()
            concat.to_excel(path)

        else:
            concat = pd.read_excel(path)
        index_df = concat.groupby("Condition").apply(self.changing_index)
        index_df = index_df.apply(correct_clusters, axis=1)

        # index_df = (
        # pd.concat([index_df.reset_index()[["Condition", "Introns_nb"]], pd.DataFrame.from_records(index_df.reset_index()[0])], axis=1)
        # .set_index(["Condition", "Introns_nb"])
        # .T
        # )
        # index_df.index.name = "level_1"
        # index_df = index_df.T
        # print(index_df)

        return concat, index_df

    def concat_stats(self, l_fold_df, conditions, path, keys, cluster_nb):
        if os.path.isfile(path) is False:
            concat = pd.concat(l_fold_df)
            # print(concat)
            # concat["Condition"] = [sub_e for e in conditions for sub_e in [e] * (len(keys) * cluster_nb)]
            concat.index.names = ["Cluster"]
            # concat["Introns_nb"] = len(conditions) * [sub_i for i in keys for sub_i in [i] * cluster_nb]
            concat = concat.reset_index()
            # concat = concat[["Cluster", "Fold", "Condition", "Introns_nb"]]
            concat["id"] = concat["Condition"].astype(str) + "-" + concat["Introns_nb"].astype(str) + "-" + concat["Cluster"].astype(str)
            concat = concat.sort_values(by=["Condition", "Introns_nb", "Cluster"])
            # print(concat)
            # exit()
        else:
            pass
        return concat

    @staticmethod
    def dict_new_clusters(index_df, keys, cluster_nb):
        map_df_new_clusters = (
            index_df.reset_index()
            .melt(id_vars=["Introns_nb", "Condition"], value_vars=list(range(0, cluster_nb)))
            .rename({"level_1": "cluster_new", "value": "cluster"}, axis=1)
            .sort_values(by=["Introns_nb", "cluster"])
        )

        map_df_new_clusters["id"] = (
            map_df_new_clusters["Condition"].astype(str)
            + "-"
            + map_df_new_clusters["Introns_nb"].astype(str)
            + "-"
            + map_df_new_clusters["cluster"].astype(str)
        )

        map_d = map_df_new_clusters[["id", "cluster_new"]].set_index("id").to_dict()["cluster_new"]

        return map_d

    def data_profiles(self, map_d, concat, keys, cluster_nb):

        k = 10 if len(keys) >= 10 else max(keys)
        concat = concat.reset_index()
        fillna = concat.fillna(concat.median())[["Cluster", "Condition", "Introns_nb"] + list(range(0, k))].reset_index()
        # fillna = fillna.rename({"index": "cluster"}, axis=1)
        fillna["id"] = fillna["Condition"].astype(str) + "-" + fillna["Introns_nb"].astype(str) + "-" + fillna["Cluster"].astype(str)
        fillna["cluster_new"] = fillna["id"].map(map_d)
        # pd.options.display.max_rows = 120
        # print(fillna)
        # exit()
        # data_profiles = fillna.groupby(["Condition", "cluster_new"]).mean().drop(["Cluster", "Introns_nb"], axis=1)
        data_profiles = fillna.drop(["Cluster", "id"], axis=1)
        data_profiles = data_profiles.melt(id_vars=["Condition", "Introns_nb", "cluster_new"], value_vars=list(range(0, k)))
        data_profiles["value"] = 100 * data_profiles["value"]
        data_profiles["variable"] = data_profiles["variable"] + 1

        return data_profiles

    def fold_reordered(self, map_d, concat, keys, cluster_nb):
        # print(concat)
        concat["cluster_new"] = concat["id"].map(map_d)

        def pivot(r):
            return r.pivot(columns=["cluster_new"], index="Introns_nb", values="Fold")

        # ! WARNING
        # concat = concat.dropna(subset=["cluster_new"])
        # print(concat)

        data = concat.groupby("Condition").apply(pivot)
        # print(data)

        data = data * 100
        return data

    @staticmethod
    def data_barplot(data, keys, cluster_nb):
        data = data[["False_raw", "True_raw", "True", "False", "Condition", "Introns_nb"]]
        data = data.groupby(["Condition", "Introns_nb"]).sum()
        return data

    @staticmethod
    def heatmap(data_heatmap, data_barplot, data_profiles, concat_stats, path, keys, condition, filtering_used, up_lim=200, down_lim=0):
        # print(data_heatmap)
        # print(data_heatmap.reset_index(drop=True).drop(["Condition"], axis=1).set_index(["Introns_nb"]).fillna(100))
        # exit()
        # print(data_barplot)
        # print(data_profiles)
        modif_keys_intron_positions_columns = [e - keys[0] for e in keys[:10]]
        modif_keys_intron_positions_labels = [e - keys[0] + 1 for e in keys[:10]]

        f = plt.figure(constrained_layout=True, figsize=(30, 30))
        gs = f.add_gridspec(nrows=3, ncols=6, width_ratios=[2, 2, 2, 2, 1, 1], height_ratios=[10, 4, 2])
        ax_heatmap = f.add_subplot(gs[0, :-2])
        ax_barplot1 = f.add_subplot(gs[0, -2])
        ax_barplot2 = f.add_subplot(gs[0, -1])

        ax_cluster1 = f.add_subplot(gs[1, 0])
        ax_cluster2 = f.add_subplot(gs[1, 1])
        ax_cluster3 = f.add_subplot(gs[1, 2])
        ax_cluster4 = f.add_subplot(gs[1, 3])

        ax_h_barplot = f.add_subplot(gs[2, :-2])
        # print(data_heatmap)
        # print(data_heatmap.reset_index(drop=True).drop(["Condition"], axis=1).set_index(["Introns_nb"]).fillna(100))
        # exit()

        data_heatmap = (
            data_heatmap.reset_index(drop=True).drop(["Condition"], axis=1).set_index(["Introns_nb"]).fillna(0).T.reset_index(drop=True).T
        )
        print(data_heatmap)
        # exit()

        stats_heatmap_false = (
            concat_stats[["cluster_new", "False_raw", "Introns_nb"]]
            .pivot(index=["Introns_nb"], columns=["cluster_new"], values=["False_raw"])
            .rename({"False_raw": "Count"}, axis=1)
            .T.reset_index(drop=True)
            .T.fillna(0)
            .astype(int)
            .astype(str)
        )

        stats_heatmap_true = (
            concat_stats[["cluster_new", "True_raw", "Introns_nb"]]
            .pivot(index=["Introns_nb"], columns=["cluster_new"], values=["True_raw"])
            .rename({"True_raw": "Count"}, axis=1)
            .T.reset_index(drop=True)
            .T.fillna(0)
            .astype(int)
            .astype(str)
        )

        data_heatmap_to_display = data_heatmap.round(0).applymap(lambda x: "+" + str(int(x)) if x > 0 else str(int(x)))
        stats_heatmap = stats_heatmap_false + " / " + stats_heatmap_true + " / " + data_heatmap_to_display.astype(str) + "%"

        sns.heatmap(
            data=data_heatmap,
            # center=int((up_lim - down_lim) / 2),
            center=0,
            vmin=down_lim,
            vmax=up_lim,
            annot=stats_heatmap.values.tolist(),
            fmt="",
            cmap="coolwarm",
            linecolor="grey",
            lw=0.01,
            ax=ax_heatmap,
            cbar_kws={"shrink": 0.4, "use_gridspec": False, "location": "left", "label": "% enrichment"},
        )
        ax_heatmap.set_xticklabels(["Cluster {}".format(i) for i in range(0, 4)])
        ax_heatmap.set_ylabel("Gene introns nb")
        ax_heatmap.set_title(
            "{} VS {} enrichment by cluster and by group for condition : {} | Filtering : {}".format(
                "True", "False", condition, filtering_used
            )
        )

        sns.barplot(
            data=data_barplot.melt(id_vars=["Introns_nb"], value_vars=["True_raw", "False_raw"]),
            y="Introns_nb",
            x="value",
            hue="variable",
            ax=ax_barplot1,
        )
        ax_barplot1.set_yticklabels([])
        ax_barplot1.set_ylabel("")
        ax_barplot1.set_xlabel("Genes count")
        ax_barplot1.yaxis.set_ticks_position("none")
        ax_barplot1.spines["right"].set_linewidth(0)
        ax_barplot1.spines["top"].set_linewidth(0)
        ax_barplot1.legend().remove()
        ax_barplot1.grid(axis="x")

        sns.barplot(
            data=data_barplot.melt(id_vars=["Introns_nb"], value_vars=["True", "False"]),
            y="Introns_nb",
            x="value",
            hue="variable",
            ax=ax_barplot2,
        )
        ax_barplot2.set_yticklabels([])
        ax_barplot2.set_ylabel("")
        ax_barplot2.set_xlabel("% compare\nto all genes")
        ax_barplot2.yaxis.set_ticks_position("none")
        ax_barplot2.spines["right"].set_linewidth(0)
        ax_barplot2.spines["top"].set_linewidth(0)
        ax_barplot2.grid(axis="x")

        # ax_barplot2.legend().remove()
        max_ylim = 55
        # print(data_profiles.loc[data_profiles["cluster_new"] == 0]["value"].max())
        sns.lineplot(data=data_profiles.loc[data_profiles["cluster_new"] == 0], x="variable", y="value", color="red", ax=ax_cluster1)
        # ax_cluster1.plot(modif_keys_intron_positions_labels, data_profiles[modif_keys_intron_positions_columns].loc[0].values, color="red")
        ax_cluster1.set_xlabel("Intron ordinal position")
        ax_cluster1.set_xticks([i for i in range(1, 11)])
        ax_cluster1.set_ylim(0, max_ylim)
        ax_cluster1.set_ylabel("Intron size\n(% gene body)")
        ax_cluster1.grid(axis="y")

        # ax.yaxis.get_minor_ticks()
        # ax_cluster1.set_xticklabels([str(i) for i in range(1,11)])

        sns.lineplot(data=data_profiles.loc[data_profiles["cluster_new"] == 1], x="variable", y="value", color="red", ax=ax_cluster2)
        # ax_cluster2.plot(modif_keys_intron_positions_labels, data_profiles[modif_keys_intron_positions_columns].loc[1].values, color="red")
        ax_cluster2.set_xlabel("Intron ordinal position")
        ax_cluster2.set_xticks([i for i in range(1, 11)])
        ax_cluster2.set_ylim(0, max_ylim)
        # ax_cluster2.set_yticks([])
        ax_cluster2.grid(axis="y")

        sns.lineplot(data=data_profiles.loc[data_profiles["cluster_new"] == 2], x="variable", y="value", color="red", ax=ax_cluster3)
        # ax_cluster3.plot(modif_keys_intron_positions_labels, data_profiles[modif_keys_intron_positions_columns].loc[2].values, color="red")
        ax_cluster3.set_xlabel("Intron ordinal position")
        ax_cluster3.set_xticks([i for i in range(1, 11)])
        ax_cluster3.set_ylim(0, max_ylim)
        # ax_cluster3.set_yticks([])
        ax_cluster3.grid(axis="y")

        sns.lineplot(data=data_profiles.loc[data_profiles["cluster_new"] == 3], x="variable", y="value", color="red", ax=ax_cluster4)
        # ax_cluster4.plot(modif_keys_intron_positions_labels, data_profiles[modif_keys_intron_positions_columns].loc[3].values, color="red")
        ax_cluster4.set_xlabel("Intron ordinal position")
        ax_cluster4.set_xticks([i for i in range(1, 11)])
        ax_cluster4.set_ylim(0, max_ylim)
        # ax_cluster4.set_yticks([])
        ax_cluster4.grid(axis="y")

        barplot_h_data = (
            concat_stats[["cluster_new", "False_raw", "True_raw"]]
            .groupby("cluster_new")
            .sum()
            .reset_index()
            .rename({"False_raw": "False", "True_raw": "True"}, axis=1)
            .melt(id_vars="cluster_new", value_vars=["True", "False"])
        )

        sns.barplot(data=barplot_h_data, x="cluster_new", y="value", hue="variable", ax=ax_h_barplot)

        f.savefig(path)
        # exit()

    @staticmethod
    def heatmap_two(data_heatmap, data_barplot, data_profiles, concat_stats, path, keys, condition, filtering_used, up_lim=200, down_lim=0):
        # print(data_heatmap)
        # print(data_heatmap.reset_index(drop=True).drop(["Condition"], axis=1).set_index(["Introns_nb"]).fillna(100))
        # exit()
        # print(data_barplot)
        # print(data_profiles)
        modif_keys_intron_positions_columns = [e - keys[0] for e in keys[:10]]
        modif_keys_intron_positions_labels = [e - keys[0] + 1 for e in keys[:10]]

        f = plt.figure(constrained_layout=True, figsize=(30, 30))
        gs = f.add_gridspec(nrows=3, ncols=4, width_ratios=[2, 2, 1, 1], height_ratios=[10, 4, 2])
        ax_heatmap = f.add_subplot(gs[0, :-2])
        ax_barplot1 = f.add_subplot(gs[0, -2])
        ax_barplot2 = f.add_subplot(gs[0, -1])

        ax_cluster1 = f.add_subplot(gs[1, 0])
        ax_cluster2 = f.add_subplot(gs[1, 1])
        # ax_cluster3 = f.add_subplot(gs[1, 2])
        # ax_cluster4 = f.add_subplot(gs[1, 3])

        ax_h_barplot = f.add_subplot(gs[2, :-2])
        # print(data_heatmap)
        # print(data_heatmap.reset_index(drop=True).drop(["Condition"], axis=1).set_index(["Introns_nb"]).fillna(100))
        # exit()

        # data_heatmap = (
        # data_heatmap.reset_index(drop=True).drop(["Condition"], axis=1).set_index(["Introns_nb"]).fillna(0).T.reset_index(drop=True).T
        # )

        concat_stats = concat_stats.rename({"cluster_lite": "cluster_new", "Introns_nb_bins": "Introns_nb"}, axis=1)
        print(concat_stats)
        # data_heatmap = concat_stats.pivot(index="Introns_nb", columns="cluster_new", values="Fold")
        data_heatmap = data_heatmap

        stats_heatmap_false = (
            concat_stats[["cluster_new", "False_raw", "Introns_nb"]]
            .pivot(index=["Introns_nb"], columns=["cluster_new"], values=["False_raw"])
            .rename({"False_raw": "Count"}, axis=1)
            .T.reset_index(drop=True)
            .T.fillna(0)
            .astype(int)
            .astype(str)
        )

        stats_heatmap_true = (
            concat_stats[["cluster_new", "True_raw", "Introns_nb"]]
            .pivot(index=["Introns_nb"], columns=["cluster_new"], values=["True_raw"])
            .rename({"True_raw": "Count"}, axis=1)
            .T.reset_index(drop=True)
            .T.fillna(0)
            .astype(int)
            .astype(str)
        )

        data_heatmap_to_display = data_heatmap.round(0).applymap(lambda x: "+" + str(int(x)) if x > 0 else str(int(x)))
        stats_heatmap = stats_heatmap_false + " / " + stats_heatmap_true + " / " + data_heatmap_to_display.astype(str) + "%"
        print(stats_heatmap)

        sns.heatmap(
            data=data_heatmap,
            # center=int((up_lim - down_lim) / 2),
            center=0,
            vmin=down_lim,
            vmax=up_lim,
            annot=stats_heatmap.values.tolist(),
            fmt="",
            cmap="coolwarm",
            linecolor="grey",
            lw=0.01,
            ax=ax_heatmap,
            cbar_kws={"shrink": 0.4, "use_gridspec": False, "location": "left", "label": "% enrichment"},
        )
        ax_heatmap.set_xticklabels(["Cluster {}".format(i) for i in range(0, 4)])
        ax_heatmap.set_ylabel("Gene introns nb")
        ax_heatmap.set_title(
            "{} VS {} enrichment by cluster and by group for condition : {} | Filtering : {}".format(
                "True", "False", condition, filtering_used
            )
        )

        sns.barplot(
            data=data_barplot.melt(id_vars=["Introns_nb"], value_vars=["True_raw", "False_raw"]),
            y="Introns_nb",
            x="value",
            hue="variable",
            ax=ax_barplot1,
        )
        ax_barplot1.set_yticklabels([])
        ax_barplot1.set_ylabel("")
        ax_barplot1.set_xlabel("Genes count")
        ax_barplot1.yaxis.set_ticks_position("none")
        ax_barplot1.spines["right"].set_linewidth(0)
        ax_barplot1.spines["top"].set_linewidth(0)
        ax_barplot1.legend().remove()
        ax_barplot1.grid(axis="x")

        sns.barplot(
            data=data_barplot.melt(id_vars=["Introns_nb"], value_vars=["True", "False"]),
            y="Introns_nb",
            x="value",
            hue="variable",
            ax=ax_barplot2,
        )
        ax_barplot2.set_yticklabels([])
        ax_barplot2.set_ylabel("")
        ax_barplot2.set_xlabel("% compare\nto all genes")
        ax_barplot2.yaxis.set_ticks_position("none")
        ax_barplot2.spines["right"].set_linewidth(0)
        ax_barplot2.spines["top"].set_linewidth(0)
        ax_barplot2.grid(axis="x")
        # ax_barplot2.legend().remove()

        print(data_profiles)

        max_ylim = 55
        # print(data_profiles.loc[data_profiles["cluster_new"] == 0]["value"].max())
        sns.lineplot(data=data_profiles.loc[data_profiles["cluster_new"] == 0], x="variable", y="value", color="red", ax=ax_cluster1)
        sns.lineplot(data=data_profiles.loc[data_profiles["cluster_new"] == 1], x="variable", y="value", color="green", ax=ax_cluster1)
        sns.lineplot(data=data_profiles.loc[data_profiles["cluster_new"] == 2], x="variable", y="value", color="blue", ax=ax_cluster1)

        # ax_cluster1.plot(modif_keys_intron_positions_labels, data_profiles[modif_keys_intron_positions_columns].loc[0].values, color="red")
        ax_cluster1.set_xlabel("Intron ordinal position")
        ax_cluster1.set_xticks([i for i in range(1, 11)])
        ax_cluster1.set_ylim(0, max_ylim)
        ax_cluster1.set_ylabel("Intron size\n(% gene body)")
        ax_cluster1.grid(axis="y")
        # ax.yaxis.get_minor_ticks()
        # ax_cluster1.set_xticklabels([str(i) for i in range(1,11)])

        # ax_cluster2.plot(modif_keys_intron_positions_labels, data_profiles[modif_keys_intron_positions_columns].loc[1].values, color="red")
        # ax_cluster2.set_xlabel("Intron ordinal position")
        # ax_cluster2.set_xticks([i for i in range(1, 11)])
        # ax_cluster2.set_ylim(0, max_ylim)
        # ax_cluster2.set_yticks([])

        # ax_cluster3.plot(modif_keys_intron_positions_labels, data_profiles[modif_keys_intron_positions_columns].loc[2].values, color="red")
        # ax_cluster3.set_xlabel("Intron ordinal position")
        # ax_cluster3.set_xticks([i for i in range(1, 11)])
        # ax_cluster3.set_ylim(0, max_ylim)
        # ax_cluster3.set_yticks([])

        sns.lineplot(data=data_profiles.loc[data_profiles["cluster_new"] == 3], x="variable", y="value", color="black", ax=ax_cluster2)
        # ax_cluster4.plot(modif_keys_intron_positions_labels, data_profiles[modif_keys_intron_positions_columns].loc[3].values, color="red")
        ax_cluster2.set_xlabel("Intron ordinal position")
        ax_cluster2.set_xticks([i for i in range(1, 11)])
        ax_cluster2.set_ylim(0, max_ylim)
        ax_cluster2.set_ylabel("")
        ax_cluster2.grid(axis="y")
        # ax_cluster2.set_yticks([])

        barplot_h_data = (
            concat_stats[["cluster_new", "False_raw", "True_raw"]]
            .groupby("cluster_new")
            .sum()
            .reset_index()
            .rename({"False_raw": "False", "True_raw": "True"}, axis=1)
            .melt(id_vars="cluster_new", value_vars=["True", "False"])
        )

        sns.barplot(data=barplot_h_data, x="cluster_new", y="value", hue="variable", ax=ax_h_barplot)

        f.savefig(path)
        # exit()

    @staticmethod
    def heatmap_three(
        data_heatmap, data_barplot, data_profiles, concat_stats, path, keys, condition, filtering_used, up_lim=200, down_lim=0
    ):
        # print(data_heatmap)
        # print(data_heatmap.reset_index(drop=True).drop(["Condition"], axis=1).set_index(["Introns_nb"]).fillna(100))
        # exit()
        # print(data_barplot)
        # print(data_profiles)

        color_miso = "#C43032"
        color_siso = "#7B9FF9"
        palette = [color_miso, color_siso]

        modif_keys_intron_positions_columns = [e - keys[0] for e in keys[:10]]
        modif_keys_intron_positions_labels = [e - keys[0] + 1 for e in keys[:10]]

        f = plt.figure(constrained_layout=True, figsize=(40, 30))
        # f = plt.figure(constrained_layout=True, figsize=(40, 40))
        gs = f.add_gridspec(nrows=4, ncols=6, width_ratios=[1, 1, 1, 1, 0.5, 0.5], height_ratios=[5, 3, 1.5, 1.5])
        gs.update(hspace=0.05)

        ax_heatmap = f.add_subplot(gs[0, :-2])
        ax_barplot1 = f.add_subplot(gs[0, -2])
        ax_barplot2 = f.add_subplot(gs[0, -1])

        ax_cluster1 = f.add_subplot(gs[1, 0])
        ax_cluster2 = f.add_subplot(gs[1, 1])
        ax_cluster3 = f.add_subplot(gs[1, 2])
        ax_cluster4 = f.add_subplot(gs[1, 3])
        ax_legend = f.add_subplot(gs[1, 4])

        ax_h_barplot1 = f.add_subplot(gs[2, :-2])
        ax_h_barplot2 = f.add_subplot(gs[3, :-2])
        # print(data_heatmap)
        # print(data_heatmap.reset_index(drop=True).drop(["Condition"], axis=1).set_index(["Introns_nb"]).fillna(100))
        # exit()

        # data_heatmap = (
        # data_heatmap.reset_index(drop=True).drop(["Condition"], axis=1).set_index(["Introns_nb"]).fillna(0).T.reset_index(drop=True).T
        # )

        concat_stats = concat_stats.rename({"cluster_lite": "cluster_new", "Introns_nb_bins": "Introns_nb"}, axis=1)
        print(concat_stats)
        # data_heatmap = concat_stats.pivot(index="Introns_nb", columns="cluster_new", values="Fold")
        data_heatmap = data_heatmap

        stats_heatmap_false = (
            concat_stats[["cluster_new", "False_raw", "Introns_nb"]]
            .pivot(index=["Introns_nb"], columns=["cluster_new"], values=["False_raw"])
            .rename({"False_raw": "Count"}, axis=1)
            .T.reset_index(drop=True)
            .T.fillna(0)
            .astype(int)
            .astype(str)
        )

        stats_heatmap_true = (
            concat_stats[["cluster_new", "True_raw", "Introns_nb"]]
            .pivot(index=["Introns_nb"], columns=["cluster_new"], values=["True_raw"])
            .rename({"True_raw": "Count"}, axis=1)
            .T.reset_index(drop=True)
            .T.fillna(0)
            .astype(int)
            .astype(str)
        )

        data_heatmap_to_display = data_heatmap.round(0).applymap(lambda x: "+" + str(int(x)) if x > 0 else str(int(x)))
        stats_heatmap = stats_heatmap_false + " / " + stats_heatmap_true + " / " + data_heatmap_to_display.astype(str) + "%"
        print(stats_heatmap)
        # stats_heatmap = data_heatmap_to_display.astype(str) + "%"
        # print(stats_heatmap)

        sns.heatmap(
            data=data_heatmap,
            # center=int((up_lim - down_lim) / 2),
            center=0,
            vmin=down_lim,
            vmax=up_lim,
            annot=stats_heatmap.values.tolist(),
            fmt="",
            cmap="coolwarm",
            linecolor="grey",
            lw=0.01,
            ax=ax_heatmap,
            cbar_kws={"shrink": 0.4, "use_gridspec": False, "location": "left", "label": "% enrichment"},
        )
        ax_heatmap.set_xticklabels(["Cluster {}".format(i) for i in range(1, 5)])
        ax_heatmap.set_xlabel("")
        ax_heatmap.set_ylabel("Gene introns nb")
        ax_heatmap.set_title(
            "{} VS {} enrichment by cluster and by group for condition : {} | Filtering : {}".format(
                "True", "False", condition, filtering_used
            )
        )
        print(data_barplot.melt(id_vars=["Introns_nb"], value_vars=["True_raw", "False_raw"]))
        sns.barplot(
            data=data_barplot.melt(id_vars=["Introns_nb"], value_vars=["True_raw", "False_raw"]),
            y="Introns_nb",
            x="value",
            hue="variable",
            ax=ax_barplot1,
            palette=palette,
        )
        ax_barplot1.set_yticklabels([])
        ax_barplot1.set_ylabel("")
        ax_barplot1.set_xlabel("Genes count")
        ax_barplot1.yaxis.set_ticks_position("none")
        ax_barplot1.spines["right"].set_linewidth(0)
        ax_barplot1.spines["top"].set_linewidth(0)
        ax_barplot1.legend().remove()
        ax_barplot1.grid(axis="x")
        ax_barplot1.set_axisbelow(True)

        print(data_barplot.melt(id_vars=["Introns_nb"], value_vars=["True", "False"]))
        sns.barplot(
            data=data_barplot.melt(id_vars=["Introns_nb"], value_vars=["True", "False"]),
            y="Introns_nb",
            x="value",
            hue="variable",
            ax=ax_barplot2,
            palette=palette,
        )
        ax_barplot2.set_yticklabels([])
        ax_barplot2.set_ylabel("")
        ax_barplot2.set_xlabel("% compare\nto all genes")
        ax_barplot2.yaxis.set_ticks_position("none")
        ax_barplot2.spines["right"].set_linewidth(0)
        ax_barplot2.spines["top"].set_linewidth(0)
        ax_barplot2.grid(axis="x")
        ax_barplot2.legend().remove()
        ax_barplot2.set_axisbelow(True)

        print(data_profiles)

        max_ylim = 55
        # max_ylim = data_profiles.loc[data_profiles["cluster_new"] == 0].value.max()
        # print(data_profiles.loc[data_profiles["cluster_new"] == 0]["value"].max())
        sns.lineplot(data=data_profiles.loc[data_profiles["cluster_new"] == 0], x="variable", y="value", color="red", ax=ax_cluster1)

        # ax_cluster1.plot(modif_keys_intron_positions_labels, data_profiles[modif_keys_intron_positions_columns].loc[0].values, color="red")
        ax_cluster1.set_xlabel("Intron ordinal position")
        ax_cluster1.set_ylabel("")
        ax_cluster1.set_xticks([i for i in range(1, 11)])
        ax_cluster1.set_ylim(0, max_ylim)
        ax_cluster1.set_ylabel("Intron size\n(% gene body)")
        ax_cluster1.grid(axis="y")
        # ax.yaxis.get_minor_ticks()
        # ax_cluster1.set_xticklabels([str(i) for i in range(1,11)])

        # ax_cluster2.plot(modif_keys_intron_positions_labels, data_profiles[modif_keys_intron_positions_columns].loc[1].values, color="red")
        sns.lineplot(data=data_profiles.loc[data_profiles["cluster_new"] == 1], x="variable", y="value", color="green", ax=ax_cluster2)
        ax_cluster2.set_xlabel("Intron ordinal position")
        ax_cluster2.set_ylabel("")
        ax_cluster2.set_xticks([i for i in range(1, 11)])
        ax_cluster2.set_ylim(0, max_ylim)
        # ax_cluster2.set_yticks([])
        ax_cluster2.grid(axis="y")

        # ax_cluster3.plot(modif_keys_intron_positions_labels, data_profiles[modif_keys_intron_positions_columns].loc[2].values, color="red")
        sns.lineplot(data=data_profiles.loc[data_profiles["cluster_new"] == 2], x="variable", y="value", color="blue", ax=ax_cluster3)
        ax_cluster3.set_xlabel("Intron ordinal position")
        ax_cluster3.set_ylabel("")
        ax_cluster3.set_xticks([i for i in range(1, 11)])
        ax_cluster3.set_ylim(0, max_ylim)
        # ax_cluster3.set_yticks([])
        ax_cluster3.grid(axis="y")

        sns.lineplot(data=data_profiles.loc[data_profiles["cluster_new"] == 3], x="variable", y="value", color="black", ax=ax_cluster4)
        # ax_cluster4.plot(modif_keys_intron_positions_labels, data_profiles[modif_keys_intron_positions_columns].loc[3].values, color="red")
        ax_cluster4.set_xlabel("Intron ordinal position")
        ax_cluster4.set_ylabel("")
        ax_cluster4.set_xticks([i for i in range(1, 11)])
        ax_cluster4.set_ylim(0, max_ylim)
        ax_cluster4.set_ylabel("")
        ax_cluster4.grid(axis="y")
        # ax_cluster2.set_yticks([])

        import matplotlib.patches as mpatches

        miso_patch = mpatches.Patch(color=color_miso, label="Multi-isoforms")
        siso_patch = mpatches.Patch(color=color_siso, label="Single-isoform")
        ax_legend.spines["top"].set_linewidth(0)
        ax_legend.spines["bottom"].set_linewidth(0)
        ax_legend.spines["right"].set_linewidth(0)
        ax_legend.spines["left"].set_linewidth(0)
        ax_legend.yaxis.set_ticks_position("none")
        ax_legend.xaxis.set_ticks_position("none")
        ax_legend.legend(handles=[miso_patch, siso_patch], loc="upper center", fontsize="x-large")
        ax_legend.set_yticklabels([])
        ax_legend.set_xticklabels([])

        barplot_h_data = (
            concat_stats[["cluster_new", "False_raw", "True_raw"]]
            .groupby("cluster_new")
            .sum()
            .reset_index()
            .rename({"False_raw": "False", "True_raw": "True"}, axis=1)
            .melt(id_vars="cluster_new", value_vars=["True", "False"])
        )
        print(barplot_h_data)

        sns.barplot(data=barplot_h_data, x="cluster_new", y="value", hue="variable", ax=ax_h_barplot1, palette=palette)
        ax_h_barplot1.grid(axis="y", zorder=0)
        ax_h_barplot1.set_axisbelow(True)
        ax_h_barplot1.legend().remove()
        ax_h_barplot1.set_xticklabels(["Cluster {}".format(i) for i in range(1, 5)])
        ax_h_barplot1.set_xlabel("")
        ax_h_barplot1.set_ylabel("Genes count")
        ax_h_barplot1.spines["right"].set_linewidth(0)
        ax_h_barplot1.spines["top"].set_linewidth(0)

        def barplot_h_fct(r):
            r["ratio"] = 100 * (r["value"] / r["value"].sum())
            return r

        barplot_h_data_ratio = barplot_h_data.groupby("variable").apply(barplot_h_fct)
        print(barplot_h_data_ratio)

        sns.barplot(data=barplot_h_data_ratio, x="cluster_new", y="ratio", hue="variable", ax=ax_h_barplot2, palette=palette)
        ax_h_barplot2.grid(axis="y", zorder=0)
        ax_h_barplot2.set_axisbelow(True)
        ax_h_barplot2.legend().remove()
        ax_h_barplot2.set_xticklabels(["Cluster {}".format(i) for i in range(1, 5)])
        ax_h_barplot2.set_xlabel("")
        ax_h_barplot2.set_ylabel("% compare\nto all genes")
        ax_h_barplot2.spines["right"].set_linewidth(0)
        ax_h_barplot2.spines["top"].set_linewidth(0)

        f.savefig(path)
        # exit()

    def launch_heatmap(
        self,
        df_raw,
        data_heatmap,
        data_barplot,
        data_profiles,
        concat_stats,
        conditions,
        keys,
        sign,
        seed,
        # down_lim=[-50, -100, -200, -500],
        # down_lim=[-50],
        down_lim=[-25],
        # up_lim=[50, 100, 200, 500],
        # up_lim=[50],
        up_lim=[25],
    ):
        data_heatmap = data_heatmap.reset_index()
        data_barplot = data_barplot.reset_index()
        data_barplot["Introns_nb"] = data_barplot["Introns_nb"].astype(str)
        data_profiles = data_profiles.reset_index()
        # print(data_heatmap)
        # print(data_barplot)
        # print(data_profiles)

        path = "/gstock/EXOTIC/data/CLUSTERING/{}/heatmap_test11_{}_{}_{}_{}_{}.png"

        cluster_lite = {
            0: 0,
            1: 0,
            2: 0,
            3: 1,
        }

        for up, down in zip(up_lim, down_lim):
            for condition in conditions:

                def test_stats(r, stats_raw):
                    # print(r)
                    # print(r[["cluster_lite", "True_raw", "False_raw"]].T)
                    l_or = list()
                    l_pvalue = list()
                    for cluster in r["cluster_lite"].unique().tolist():
                        # print(cluster)
                        true_condition = int(r.loc[r["cluster_lite"] == cluster]["True_raw"].values[0])
                        false_condition = int(r.loc[r["cluster_lite"] == cluster]["False_raw"].values[0])
                        true_not_condition = int(r.loc[r["cluster_lite"] != cluster]["True_raw"].sum())
                        false_not_condition = int(r.loc[r["cluster_lite"] != cluster]["False_raw"].sum())
                        # print([[true_condition, true_not_condition], [false_condition, false_not_condition]])
                        odd_ratio, pvalue = fisher_exact([[true_condition, true_not_condition], [false_condition, false_not_condition]])
                        l_or.append(odd_ratio)
                        l_pvalue.append(pvalue)
                        # print(odd_ratio, pvalue)

                    r["OR"] = l_or
                    r["p_value"] = l_pvalue
                    r["Ratio"] = r["True_raw"] / r["False_raw"]
                    r["Fold"] = 100 * ((r["Ratio"] / (r["True_raw"].sum() / r["False_raw"].sum())) - 1)

                    r["True"] = 100 * (r["True_raw"] / stats_raw["True_raw"])
                    r["False"] = 100 * (r["False_raw"] / stats_raw["False_raw"])
                    return r

                print(up, down, condition)
                data_heatmap_condition = data_heatmap.loc[data_heatmap["Condition"] == condition].reset_index(drop=True)
                data_barplot_condition = data_barplot.loc[data_barplot["Condition"] == condition].reset_index(drop=True)
                data_profiles_condition = data_profiles.loc[data_profiles["Condition"] == condition].reset_index(drop=True)
                data_heatmap_stats_condition = concat_stats.loc[concat_stats["Condition"] == condition].reset_index(drop=True)

                # print(data_heatmap_condition)
                # print(data_barplot_condition)
                # print(data_profiles_condition)
                # print(data_heatmap_stats_condition)
                # print(df_raw[condition].value_counts().to_dict())

                # data_heatmap_condition["cluster_lite"] = data_heatmap_condition["Cluster"].map(cluster_lite)
                # data_barplot_condition["cluster_lite"] = data_barplot_condition["cluster_new"].map(cluster_lite)
                # data_profiles_condition["cluster_lite"] = data_profiles_condition["cluster_new"].map(cluster_lite)
                # data_heatmap_stats_condition["cluster_lite"] = data_heatmap_stats_condition["cluster_new"].map(cluster_lite)
                # data_profiles_condition = data_profiles_condition.rename({"cluster_new": "cluster_lite"}, axis=1)
                data_heatmap_stats_condition = data_heatmap_stats_condition.rename({"cluster_new": "cluster_lite"}, axis=1)

                ####

                data_heatmap_stats_condition = data_heatmap_stats_condition[
                    ["Introns_nb", "cluster_lite", "False_raw", "True_raw"]
                ].sort_values(by=["Introns_nb", "cluster_lite"])

                bins = [2, 12, 21, 30]
                # bins = [2, 10, 15, 22, 30]
                # bins = [2, 16, 30]
                labels = [" - ".join([str(int(e) + 1), str(int(bins[j + 1]))]) for j, e in enumerate(bins) if j < len(bins) - 1]
                data_heatmap_stats_condition["Introns_nb_bins"] = pd.cut(
                    data_heatmap_stats_condition["Introns_nb"],
                    bins,
                    labels=labels,
                    include_lowest=False,
                )

                pd.options.display.max_rows = 200
                # print(data_heatmap_stats_condition)
                data_heatmap_stats_condition = data_heatmap_stats_condition.rename(
                    {"Introns_nb_bins": "Introns_nb", "Introns_nb": "Introns_nb_raw"}, axis=1
                )

                data_heatmap_stats_condition = data_heatmap_stats_condition.groupby(["Introns_nb", "cluster_lite"]).sum().reset_index()
                data_heatmap_stats_condition = data_heatmap_stats_condition.groupby("Introns_nb").apply(
                    lambda r: test_stats(r, data_heatmap_stats_condition[["False_raw", "True_raw"]].sum().T.to_dict())
                )

                data_heatmap_condition = data_heatmap_stats_condition.pivot(index="Introns_nb", columns="cluster_lite", values="Fold")
                data_heatmap_stats_condition["Fold"] = data_heatmap_stats_condition["Fold"].fillna(0)
                data_barplot_condition = (
                    data_heatmap_stats_condition[["Introns_nb", "True_raw", "False_raw", "True", "False"]]
                    .groupby("Introns_nb")
                    .sum()
                    .reset_index()
                )
                data_barplot_condition["Introns_nb"] = data_barplot_condition["Introns_nb"].astype(str)
                # print(data_heatmap_stats_condition)
                # print(data_heatmap_condition)
                # print(data_barplot_condition)

                ####

                # data_heatmap_stats_condition_gb_introns = data_heatmap_stats_condition.groupby(["Introns_nb"]).sum().reset_index()
                # data_heatmap_stats_condition_gb_introns["cluster_lite"] = "Sum"
                # print(data_heatmap_stats_condition_gb_cluster)
                # print(data_heatmap_stats_condition_gb_introns)
                # print(
                #     pd.concat([data_heatmap_stats_condition_gb_cluster, data_heatmap_stats_condition_gb_introns], axis=0).sort_values(
                #         by=["Introns_nb", "cluster_lite"]
                #     )
                # )
                # print(data_heatmap_condition)
                # print(data_barplot_condition)
                # print(data_profiles_condition)
                # print(data_heatmap_stats_condition)
                # exit()
                # print(data_profiles_condition)
                self.heatmap_three(
                    data_heatmap_condition,
                    data_barplot_condition,
                    data_profiles_condition,
                    data_heatmap_stats_condition,
                    path.format(self.convert_sign[sign], condition, self.convert_sign[sign], str(up), str(down), str(seed)),
                    keys,
                    condition,
                    sign,
                    up,
                    down,
                )


if __name__ == "__main__":
    c = IntronStats()

# def figure_complete():
#     f, ax = plt.subplots(nrows=len(keys), ncols=3, figsize=(25,5*len(keys)))

#     # ! TO COMPLETE

#     legend_lineplot_handles = list()
#     legend_lineplot_labels = list()
#     for (j, c), color in zip(enumerate(cluster_centers), colors):
#         ax[i][0].plot([str(e) for e in range(1, len(c)+1)], [100*e for e in c], color=color, lw=3)
#         ax[i][0].set_ylabel('Intron size (% gene body)')
#         ax[i][0].set_xlabel('Intron ordinal position')
#         ax[i][0].spines['top'].set_linewidth(0)
#         ax[i][0].spines['right'].set_linewidth(0)

#         legend_lineplot_handles.append(Line2D([0], [0], color=color, lw=3))
#         legend_lineplot_labels.append( "Cluster {}".format(str(j)))

#         if j == cluster-1:
#             ax[i][0].legend(legend_lineplot_handles, legend_lineplot_labels)

#         ax[i][1].bar((j+1) - 0.25/2, gb.loc[j]['True'], color = color, width = 0.25, hatch="o", edgecolor='black')
#         ax[i][1].bar((j+1) + 0.25/2, gb.loc[j]['False'], color = color, width = 0.25, hatch="//", edgecolor='black')
#         ax[i][1].set_ylabel('Genes ratio %')

#         ax[i][1].spines['top'].set_linewidth(0)
#         ax[i][1].spines['right'].set_linewidth(0)

#         ax[i][1].legend(
#             [
#                 mpatches.Patch(facecolor="white", hatch="o",  edgecolor='black'),
#                 mpatches.Patch(facecolor="white", hatch="//",  edgecolor='black')
#             ],
#             [
#                 "True ({})".format(gb_raw.loc['Sum']['True']), "False ({})".format(gb_raw.loc['Sum']['False']),
#             ]
#         )

#         markerline, stemlines, baseline = ax[i][2].stem(["Cluster {}".format(j)], [gb.loc[j]['Fold']], bottom=1, linefmt = 'grey', use_line_collection=True)
#         markerline.set_markerfacecolor(color)
#         ax[i][2].axhline(1, color='grey', alpha=0.7, ls='--')


#         ax[i][2].set_ylabel('Normalized fold (True VS False)')
#         ax[i][2].spines['top'].set_linewidth(0)
#         ax[i][2].spines['right'].set_linewidth(0)
#         ax[i][2].set_ylim(gb['Fold'].min() - 0.2, gb['Fold'].max() + 0.2)


# for condition in conditions:


#     d = dict()
#     d[False] = 0
#     d[True] = 0


#     keys = list(range(3,31))
#     cluster = 4


#     f, ax = plt.subplots(nrows=len(keys), ncols=3, figsize=(25,5*len(keys)))

#     print(df_raw[condition].value_counts())


#     l_fold_df = list()
#     l_clusters_df = list()

#     from tqdm import tqdm

#     for i, k in tqdm(enumerate(keys)):


#         df = df_raw.loc[(df_raw["Introns_nb"] == k)]

#         df["Introns_lengths_ratio"] = df["Introns_lengths"].apply(lambda r: pd.Series(pd.Series(r) / pd.Series(r).sum()).round(2).values.tolist())
#         df["Introns_lengths_ratio"] = df["Introns_lengths_ratio"].apply(lambda r: r[:k])

#         df = df.reset_index(drop=True)
#         d[False] += df[condition].value_counts()[False]
#         d[True] += df[condition].value_counts()[True]

#         x = 'Introns_lengths_ratio'
#         X_raw = pd.DataFrame.from_records(df[x])


#         X = X_raw.copy()

#         km = KMeans(n_clusters=cluster)
#         km.fit(X.values)
#         cluster_centers = km.cluster_centers_
#         cluster_centers = [[round(sub_e, 2) for sub_e in e] for e in cluster_centers]
#         l_clusters_df.append(cluster_centers)


#         cluster_map = pd.DataFrame()
#         cluster_map["data_index"] = X.index.values
#         cluster_map["cluster"] = km.labels_


#         X['cluster'] = cluster_map.cluster

#         X['Miso'] = df['Miso']
#         X[condition] = df[condition]
#         X[condition] = X[condition].astype(str)
#         X['Strand'] = df['Strand']
#         X['Gene'] = df['Gene']
#         X['Introns_ranges'] = df['Introns_ranges']
#         X['Introns_lengths'] = df['Introns_lengths']


#         gb_raw = X.groupby(['cluster'])[condition].value_counts().rename('Count').reset_index().pivot(index='cluster', columns=condition, values='Count')
#         gb_raw.loc['Sum'] = gb_raw.sum()

#         gb = gb_raw.copy()
#         gb.columns = [str(c) for c in gb.columns]


#         gb = gb.rename({'False' : 'False_raw', 'True' : 'True_raw'}, axis=1)


#         gb['False'] = 100 * (gb['False_raw'] / df_raw[condition].value_counts().loc[False])
#         gb['True'] = 100 * (gb['True_raw'] / df_raw[condition].value_counts().loc[True])

#         gb['Ratio'] = gb['True'] / gb['False']
#         gb['Fold'] = gb['Ratio'] / gb.loc['Sum']['Ratio']

#         gb['Total'] = gb['False'] + gb['True']

#         l_fold_df.append(gb)

#         gb = gb.drop('Sum').reset_index(drop=True)

#         legend_lineplot_handles = list()
#         legend_lineplot_labels = list()
#         for (j, c), color in zip(enumerate(cluster_centers), colors):
#             ax[i][0].plot([str(e) for e in range(1, len(c)+1)], [100*e for e in c], color=color, lw=3)
#             ax[i][0].set_ylabel('Intron size (% gene body)')
#             ax[i][0].set_xlabel('Intron ordinal position')
#             ax[i][0].spines['top'].set_linewidth(0)
#             ax[i][0].spines['right'].set_linewidth(0)

#             legend_lineplot_handles.append(Line2D([0], [0], color=color, lw=3))
#             legend_lineplot_labels.append( "Cluster {}".format(str(j)))

#             if j == cluster-1:
#                 ax[i][0].legend(legend_lineplot_handles, legend_lineplot_labels)

#             ax[i][1].bar((j+1) - 0.25/2, gb.loc[j]['True'], color = color, width = 0.25, hatch="o", edgecolor='black')
#             ax[i][1].bar((j+1) + 0.25/2, gb.loc[j]['False'], color = color, width = 0.25, hatch="//", edgecolor='black')
#             ax[i][1].set_ylabel('Genes ratio %')

#             ax[i][1].spines['top'].set_linewidth(0)
#             ax[i][1].spines['right'].set_linewidth(0)

#             ax[i][1].legend(
#                 [
#                     mpatches.Patch(facecolor="white", hatch="o",  edgecolor='black'),
#                     mpatches.Patch(facecolor="white", hatch="//",  edgecolor='black')
#                 ],
#                 [
#                     "True ({})".format(gb_raw.loc['Sum']['True']), "False ({})".format(gb_raw.loc['Sum']['False']),
#                 ]
#             )

#             markerline, stemlines, baseline = ax[i][2].stem(["Cluster {}".format(j)], [gb.loc[j]['Fold']], bottom=1, linefmt = 'grey', use_line_collection=True)
#             markerline.set_markerfacecolor(color)
#             ax[i][2].axhline(1, color='grey', alpha=0.7, ls='--')


#             ax[i][2].set_ylabel('Normalized fold (True VS False)')
#             ax[i][2].spines['top'].set_linewidth(0)
#             ax[i][2].spines['right'].set_linewidth(0)
#             ax[i][2].set_ylim(gb['Fold'].min() - 0.2, gb['Fold'].max() + 0.2)


#     fold_df = pd.concat([e.Fold.drop('Sum').sort_index() for e in l_fold_df], axis=1).T.reset_index(drop=True).reset_index().rename({'index' : 'Introns_nb'}, axis=1).fillna(1)
#     fold_df['Introns_nb'] = fold_df['Introns_nb'] + 3
#     fold_df = fold_df.set_index('Introns_nb')

#     data = pd.concat([fold_df, index_df.apply(list, axis=1).rename('Order')], axis=1)
#     data = pd.DataFrame(data.apply(lambda r: pd.Series([r[i] for i in r['Order'] if i < 4]), axis=1))
#     data = data * 100
#     data.to_excel('/gstock/EXOTIC/data/CLUSTERING/fold_enrichment_{}_{}.xlsx'.format(condition, "supequal"))


#     # Mass centers
#     l = list()
#     for i, e in enumerate(l_clusters_df):
#         tmp_df_to_add = pd.DataFrame.from_records(e)
#         tmp_df_to_add['Introns_nb'] = 3 + i
#         l.append(tmp_df_to_add)
#     concat = pd.concat(l)
#     concat = concat[['Introns_nb'] + list(range(0,30))]
#     concat.index.names = ["Cluster"]
#     concat.index = [e for e in list(concat.index)]
#     concat.to_excel('/gstock/EXOTIC/data/CLUSTERING/mass_center_matrix_{}_{}.xlsx'.format(condition, "supequal"))

#     f.savefig('/gstock/EXOTIC/data/CLUSTERING/vertical_figure_{}_{}.png'.format(condition, "supequal"))


#     # Index order
#     def return_index(r):
#         r = r.drop(['Introns_nb'], axis=1).dropna(axis=1)
#         return r.apply(lambda r: [i+1 for i,e in enumerate(r) if e == max(r)][0])

#     index_df = concat.groupby('Introns_nb').apply(return_index).reset_index().pivot(index='Introns_nb', columns='level_1', values=0)[[0,1,2]].astype(int)
#     index_df[3] = index_df.apply(lambda r: [i for i in range(1,5) if i not in r.values.tolist()][0], axis=1)
#     index_df = index_df - 1
#     # index_df.columns = [int(e) +1 for e in index_df.columns]


#     map_df_new_clusters = index_df.reset_index().melt(id_vars=['Introns_nb'], value_vars=list(range(0,4))).rename({'level_1' : 'cluster_new', 'value' : 'cluster'}, axis=1).sort_values(by=['Introns_nb', 'cluster'])
#     map_df_new_clusters['id'] = map_df_new_clusters['Introns_nb'].astype(str) + '-' + map_df_new_clusters['cluster'].astype(str)
#     map_d = map_df_new_clusters[['id', 'cluster_new']].set_index('id').to_dict()['cluster_new']

#     fillna = concat.fillna(concat.median())[['Introns_nb'] + list(range(0,10))].reset_index()
#     fillna = fillna.rename({'index' : 'cluster'}, axis=1)

#     fillna['id'] = fillna['Introns_nb'].astype(str) + '-' + fillna['cluster'].astype(str)
#     fillna['cluster_new'] = fillna['id'].map(map_d)
#     data_profiles = fillna.groupby('cluster_new').mean().drop(['cluster', 'Introns_nb'], axis=1)
#     data_profiles = 100 * data_profiles


#     fig3 = plt.figure(constrained_layout=True, figsize=(30,15))
#     gs = fig3.add_gridspec(nrows=2, ncols=6, width_ratios = [2,2,2,2,1,1], height_ratios = [10,1])
#     ax_heatmap = fig3.add_subplot(gs[0, : -2])
#     ax_barplot1 = fig3.add_subplot(gs[0, -2])
#     ax_barplot2 = fig3.add_subplot(gs[0, -1])

#     ax_cluster1 = fig3.add_subplot(gs[1, 0])
#     ax_cluster2 = fig3.add_subplot(gs[1, 1])
#     ax_cluster3 = fig3.add_subplot(gs[1, 2])
#     ax_cluster4 = fig3.add_subplot(gs[1, 3])


#     sns.heatmap(data=data, center=100, vmin=0, vmax=200, cmap='coolwarm', linecolor='grey', lw=0.01, ax=ax_heatmap, cbar_kws={"shrink": .4, "use_gridspec" : False, "location" : "left", 'label' : '% enrichment'})
#     ax_heatmap.set_xticklabels(['Cluster {}'.format(i) for i in range(0,4)])
#     ax_heatmap.set_ylabel('Gene introns nb')
#     ax_heatmap.set_title('{} VS {} enrichment by cluster and by group'.format('True', 'False'))


#     concat_gb = pd.concat(l_fold_df).reset_index()
#     concat_gb['Introns_nb'] = [sub_e for e in list(range(3,31)) for sub_e in 5*[e]]
#     concat_gb_raw = concat_gb.loc[concat_gb['cluster'] == 'Sum'].melt(id_vars='Introns_nb', value_vars=['False_raw', 'True_raw'])
#     concat_gb_raw['variable'] = concat_gb_raw['variable'].replace({'False_raw' : 'False', 'True_raw' : 'True'})
#     concat_gb_raw['Introns_nb'] = concat_gb_raw['Introns_nb'].astype(str)

#     concat_gb_ratio = concat_gb.loc[concat_gb['cluster'] == 'Sum'].melt(id_vars='Introns_nb', value_vars=['False', 'True'])
#     concat_gb_ratio['variable'] = concat_gb_ratio['variable'].replace({'False_raw' : 'False', 'True_raw' : 'True'})
#     concat_gb_ratio['Introns_nb'] = concat_gb_ratio['Introns_nb'].astype(str)


#     sns.barplot(data=concat_gb_raw, y='Introns_nb', x='value', hue='variable', ax=ax_barplot1)
#     ax_barplot1.set_yticklabels([])
#     ax_barplot1.set_ylabel('')
#     ax_barplot1.set_xlabel('Genes count')
#     ax_barplot1.yaxis.set_ticks_position('none')
#     ax_barplot1.spines['right'].set_linewidth(0)
#     ax_barplot1.spines['top'].set_linewidth(0)


#     sns.barplot(data=concat_gb_ratio, y='Introns_nb', x='value', hue='variable', ax=ax_barplot2)
#     ax_barplot2.set_yticklabels([])
#     ax_barplot2.set_ylabel('')
#     ax_barplot2.set_xlabel('% compare\nto all genes')
#     ax_barplot2.yaxis.set_ticks_position('none')
#     ax_barplot2.spines['right'].set_linewidth(0)
#     ax_barplot2.spines['top'].set_linewidth(0)
#     ax_barplot2.legend().remove()


#     map_df_new_clusters = index_df.reset_index().melt(id_vars=['Introns_nb'], value_vars=list(range(0,4))).rename({'level_1' : 'cluster_new', 'value' : 'cluster'}, axis=1).sort_values(by=['Introns_nb', 'cluster'])
#     map_df_new_clusters['id'] = map_df_new_clusters['Introns_nb'].astype(str) + '-' + map_df_new_clusters['cluster'].astype(str)
#     map_d = map_df_new_clusters[['id', 'cluster_new']].set_index('id').to_dict()['cluster_new']

#     fillna = concat.fillna(concat.median())[['Introns_nb'] + list(range(0,10))].reset_index()
#     fillna = fillna.rename({'index' : 'cluster'}, axis=1)

#     fillna['id'] = fillna['Introns_nb'].astype(str) + '-' + fillna['cluster'].astype(str)
#     fillna['cluster_new'] = fillna['id'].map(map_d)
#     data_profiles = fillna.groupby('cluster_new').mean().drop(['cluster', 'Introns_nb'], axis=1)
#     data_profiles = 100 * data_profiles


#     ax_cluster1.plot(list(range(1,11)), data_profiles.loc[0].values, color='red')
#     ax_cluster1.set_xlabel('Intron ordinal position')
#     ax_cluster1.set_xticks([i for i in range(1,11)])
#     ax_cluster1.set_ylim(0,50)
#     ax_cluster1.set_ylabel('Intron size\n(% gene body)')
#     # ax.yaxis.get_minor_ticks()
#     # ax_cluster1.set_xticklabels([str(i) for i in range(1,11)])

#     ax_cluster2.plot(list(range(1,11)), data_profiles.loc[1].values, color='red')
#     ax_cluster2.set_xlabel('Intron ordinal position')
#     ax_cluster2.set_xticks([i for i in range(1,11)])
#     ax_cluster2.set_yticks([])

#     ax_cluster3.plot(list(range(1,11)), data_profiles.loc[2].values, color='red')
#     ax_cluster3.set_xlabel('Intron ordinal position')
#     ax_cluster3.set_xticks([i for i in range(1,11)])
#     ax_cluster3.set_ylim(0,50)
#     ax_cluster3.set_yticks([])

#     ax_cluster4.plot(list(range(1,11)), data_profiles.loc[3].values, color='red')
#     ax_cluster4.set_xlabel('Intron ordinal position')
#     ax_cluster4.set_xticks([i for i in range(1,11)])
#     ax_cluster4.set_ylim(0,50)
#     ax_cluster4.set_yticks([])
#     fig3.savefig('/gstock/EXOTIC/data/CLUSTERING/heatmap_figure_{}_{}.png'.format(condition, "supequal"))
