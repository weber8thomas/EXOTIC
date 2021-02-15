import math
import os
import sys
import pandas as pd
import gzip

pd.options.mode.chained_assignment = None  # default='warn'
import multiprocessing
import parmap
import numpy as np
import collections
from tqdm import tqdm
from pandarallel import pandarallel

pandarallel.initialize(nb_workers=20, progress_bar=True)
tqdm.pandas()
from pprint import pprint
from scipy.stats import zscore
import requests
import re
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.colors as mc
from matplotlib.colors import LogNorm
from src.utils.utils import load_config_file
import json
import subprocess


class Compute_GTEx_profile:
    def __init__(self):
        gtex_tpm_stats_path = "/gstock/biolo_datasets/variation/benchmark/Databases/GTEX/V7/GTEx_TPM_STMSD_stats.parquet"
        v_profile_path = "/gstock/biolo_datasets/variation/benchmark/Databases/GTEX/V7/GTEx_TPM_V_profile.parquet"
        h_profile_path = "/gstock/biolo_datasets/variation/benchmark/Databases/GTEX/V7/GTEx_TPM_H_profile.parquet"

        refseq_x_vcf_path = "/home/weber/PycharmProjects/ExoCarto/data/2_processed/refseq_vcf_pext_phylocsf.parquet"
        refseq_x_vcf = pd.read_parquet(refseq_x_vcf_path)
        disease_genes = refseq_x_vcf["CCRS_Gene"].unique()

        if (
            os.path.isfile(v_profile_path) is False
            and os.path.isfile(h_profile_path) is False
        ):

            gtex_tpm_stats = self.load_gtex_tpm_stats(gtex_tpm_stats_path)

            pivot_table = self.pivot_table(gtex_tpm_stats)

            h_profile = self.compute_horizontal_profile(pivot_table)
            v_profile = self.compute_vertical_profile(pivot_table)
            h_profile.columns = ["ENSG", "symbol"] + list(h_profile.columns)[2:]
            v_profile.columns = ["ENSG", "symbol"] + list(h_profile.columns)[2:]

            h_profile = h_profile.dropna()
            h_profile = self.check_horizontal_profile(h_profile)
            v_profile = v_profile.dropna()
            v_profile["H_profile"] = h_profile["H_profile"]
            v_profile["H_tissues"] = h_profile["H_tissues"]

            v_profile.to_parquet(v_profile_path)
            h_profile.to_parquet(h_profile_path)

        else:

            v_profile = pd.read_parquet(v_profile_path)

        self.v_profile = v_profile

        # print(self.h_profile)
        # print(self.v_profile)
        # exit()

        # self.v_profile = self.v_profile.loc[self.v_profile['symbol'].isin(disease_genes)]
        # test_df = self.v_profile.loc[self.v_profile['H_profile'] == 'All'].set_index(['ENSG', 'symbol']).drop(['H_profile'], axis=1)
        # print(test_df)
        # test_df = (test_df >= 0.99)
        # pprint(test_df.apply(pd.value_counts).to_dict())
        # print(self.v_profile.loc[self.v_profile['H_profile'] == 'Specific'].head(10))
        # print(self.h_profile.loc[self.h_profile['H_profile'] == 'Specific'].head(10))

    @staticmethod
    def load_gtex_tpm_stats(filepath):
        df = pd.read_parquet(filepath)
        df = df.reset_index()
        return df[["Name", "Description", "SMTSD", "TPM_value_50%"]]

    @staticmethod
    def pivot_table(df):
        return pd.pivot_table(
            df,
            index=["Name", "Description"],
            columns=["SMTSD"],
            values=["TPM_value_50%"],
        )

    @staticmethod
    def compute_horizontal_profile(df):
        df.columns = df.columns.droplevel(0)

        new_df = df.apply(zscore, axis=1)
        new_df = pd.DataFrame.from_dict(dict(zip(new_df.index, new_df.values))).T
        new_df.columns = df.columns
        new_df = new_df.reset_index()
        return new_df

    @staticmethod
    def compute_vertical_profile(df):
        new_df = df.apply(lambda c: c.rank(pct=True), axis=0)
        new_df = new_df.reset_index()
        return new_df

    @staticmethod
    def check_horizontal_profile(df):
        z_score = 2.58

        def apply_specific(r):
            l = list()
            for t, e in r.items():
                if e > z_score:
                    l.append(True)
                else:
                    l.append(False)
            if True in l:
                return "Specific"
            else:
                return "All"

        def apply_tissues(r):
            l = list()
            for t, e in r.items():
                if e > z_score:
                    l.append(t)
            return l

        df["H_profile"] = df.apply(lambda r: apply_specific(r[2:]), axis=1)
        df["H_tissues"] = df.apply(lambda r: apply_tissues(r[2:-1]), axis=1)
        return df


class EXOTIC:
    def __init__(self, gtex_profile, cpus, omim_detailed):
        ## INIT SETTINGS
        self.cpus = cpus
        self.GTEx_v_profile = gtex_profile
        self.GTEx_v_profile = self.GTEx_v_profile.rename(
            {"Cells - Transformed fibroblasts": "Cells - Cultured fibroblasts"}, axis=1
        )
        self.omim_detailed = omim_detailed

        ## YAML FILES CONFIG
        yaml = load_config_file()
        self.exotic_files = yaml["EXOTIC"]
        self.output_dir = "/gstock/EXOTIC/data/"
        # pprint(exotic_files)

        ## JSON DICTS CONFIG
        self.dicts = json.load(open("src/EXOTIC_config.json"))

        # pprint(file_dict)

        omim, disease_genes_multi_iso, clinvar, pubmed = self.load_clinvar_omim()

        ## GTEx - OMIM correlation gene level
        # TODO: REDO OMIM GENE LEVEL
        # self.compare_to_omim_basic_gene_level(
        # self.GTEx_v_profile, omim, disease_genes_multi_iso
        # )
        # exit()

        if os.path.isfile(self.exotic_files["exotic_gtex_comparison_path"]) is True:
            intermediate_file = self.build_file_exotic_score()
            print(intermediate_file)

            exotic_ok = self.extract_tissues_ok(intermediate_file)
            print(exotic_ok)

            exotic_gtex = self.compare_exons_gtex_profile(exotic_ok)
            print(exotic_gtex)
            exit()

        else:
            exotic_gtex = pd.read_parquet(
                self.exotic_files["exotic_gtex_comparison_path"]
            )

        exotic_gtex = exotic_gtex.loc[
            exotic_gtex["symbol"].isin(disease_genes_multi_iso.Name.values)
        ]
        gene_exon_df = self.compare_to_omim_basic_exon_level(
            exotic_gtex, omim, disease_genes_multi_iso
        )
        # exit()
        # gene_exon_df = self.compare_to_omim_basic_gene_exon_level(
        #     exotic_gtex, omim, disease_genes_multi_iso,
        # )
        # self.mapping_clinvar(clinvar, exotic_gtex, pubmed)
        exit()

    # TODO : return files
    def load_genomic_expression(
        self,
    ):
        l_to_return = [
            e
            for e in [
                refseq_gene_level,
                refseq_detailed,
                pext,
                omim,
                clinvar,
                clinvar_pmid,
            ]
            if e is True
        ]
        variable_names = [k for k, v in locals().items() if v in l_to_return]

        ## REFSEQ GENE LEVEL
        full_human_genes = pd.read_csv(
            self.exotic_files["refseq_gene_level_path"],
            compression="gzip",
            sep="\t",
            low_memory=False,
        ).rename({"MIM": "OMIM"}, axis=1)
        full_human_genes = full_human_genes.loc[full_human_genes["mRNA_nb"] > 1]

        ## REFSEQ DETAILED
        refseq = pd.read_csv(
            self.exotic_files["refseq_path"],
            compression="gzip",
            sep="\t",
            low_memory=False,
        )
        refseq = refseq.loc[refseq["Gene"].isin(full_human_genes.Name.values)]
        refseq = refseq.dropna(subset=["HGNC"])
        refseq["HGNC"] = refseq["HGNC"].astype(int)
        refseq["Ratio_num"] = refseq["Ratio"].apply(eval)

        ## pext
        pext_refseq = pd.read_parquet(self.exotic_files["pext_refseq_path"])

        d_to_return = {n: e for n, e in zip(variable_names, l_to_return)}

        return d_to_return

    def load_clinvar_omim(
        self,
    ):
        ## CLINVAR
        clinvar = pd.read_parquet(self.exotic_files["clinvar_file_path"])
        clinvar["POS"] = clinvar["VAR_ID"].apply(lambda r: r.split("_")[1])
        clinvar["POS"] = clinvar["POS"].astype(int)
        clinvar["MC"] = clinvar["MC"].apply(lambda r: r.split(",")[0])
        pathogenic_clinvar = clinvar.loc[
            (clinvar["Status"] == "Pathogenic")
            & (~clinvar["Real_Status"].str.contains("onflicting") == True)
            & (clinvar["RS_STARS"] >= 0)
        ]
        disease_genes = pathogenic_clinvar.loc[
            pathogenic_clinvar["RS_STARS"] > 0, "GENE"
        ].unique()
        # disease_genes = pathogenic_clinvar.loc[pathogenic_clinvar['RS_STARS'] > 0, "GENE"].unique()

        ## REFSEQ GENE LEVEL
        full_human_genes = pd.read_csv(
            self.exotic_files["refseq_gene_level_path"],
            compression="gzip",
            sep="\t",
            low_memory=False,
        ).rename({"MIM": "OMIM"}, axis=1)
        full_human_genes = full_human_genes.loc[full_human_genes["mRNA_nb"] > 1]

        ## DISEASE GENES
        disease_genes_multi_iso = (
            full_human_genes.loc[
                full_human_genes["Name"].isin(pathogenic_clinvar.GENE.unique()),
                ["OMIM", "Name"],
            ]
            .dropna()
            .drop_duplicates()
        )

        ## OMIM
        omim = pd.read_csv(
            self.exotic_files["omim_detailed"],
            compression="gzip",
            sep="\t",
        )

        ## ADD HGNC TO OMIM
        omim = pd.merge(omim, full_human_genes[["Name", "OMIM", "HGNC"]], on="OMIM")
        omim = omim.dropna(subset=list(omim.columns[6:-2]), how="all")

        omim = omim.loc[omim["Name"].isin(disease_genes_multi_iso.Name.unique())]

        clinvar_pmid_mapping = pd.read_csv(
            self.exotic_files["clinvar_pmid_mapping_path"], sep="\t"
        )

        return omim, disease_genes_multi_iso, pathogenic_clinvar, clinvar_pmid_mapping

    ## TEST

    def build_file_exotic_score(self):
        if os.path.isfile(self.exotic_files["exotic_path"]) is False:
            # READ FILES
            ## REFSEQ GENE LEVEL
            full_human_genes = pd.read_csv(
                self.exotic_files["refseq_gene_level_path"],
                compression="gzip",
                sep="\t",
                low_memory=False,
            ).rename({"MIM": "OMIM"}, axis=1)
            full_human_genes = full_human_genes.loc[full_human_genes["mRNA_nb"] > 1]

            ## REFSEQ DETAILED
            refseq = pd.read_csv(
                self.exotic_files["refseq_path"],
                compression="gzip",
                sep="\t",
                low_memory=False,
            )
            refseq = refseq.loc[refseq["Gene"].isin(full_human_genes.Name.values)]
            refseq = refseq.dropna(subset=["HGNC"])
            refseq["HGNC"] = refseq["HGNC"].astype(int)
            refseq["Ratio_num"] = refseq["Ratio"].apply(eval)

            ## pext
            pext_refseq = pd.read_parquet(self.exotic_files["pext_refseq_path"])

            # print(pext_refseq.loc[pext_refseq["symbol"] == "A1CF"])

            # SELECT MULTI-ISOFORM GENES
            # TODO : get from refseq
            multi_mrna_genes = refseq.loc[refseq["mRNA_nb"] > 1, "HGNC"].unique()

            # FILTER PEXT TO GET ONLY MULTI-ISO GENES
            pext_refseq["HGNC ID"] = pext_refseq["HGNC ID"].astype(int)
            pext_refseq = pext_refseq.loc[pext_refseq["HGNC ID"].isin(multi_mrna_genes)]
            # pext_refseq = pext_refseq.head(1000)

            # print(pext_refseq)

            # ADD COLUMSN TO PEXT
            pext_refseq = pd.merge(
                pext_refseq,
                refseq[["HGNC", "ranges", "Ratio_num", "mRNA_nb"]],
                left_on=["HGNC ID", "Exon"],
                right_on=["HGNC", "ranges"],
            )
            pext_refseq["MAP"] = pext_refseq["symbol"] + "_" + pext_refseq["Exon"]

            # DROP DUPLICATES
            pext_refseq = (
                pext_refseq.drop_duplicates()
                .sort_values(by="MAP")
                .reset_index(drop=True)
            )

            # RENAME COLUMNS & REORDER
            pext_refseq = pext_refseq.rename(self.dicts["convert_tissue_dict"], axis=1)

            pext_refseq = pext_refseq[
                [
                    "symbol",
                    "ensg",
                    "HGNC",
                    "Exon",
                    "Ratio_num",
                    "mRNA_nb",
                    "MAP",
                    "mean_proportion",
                    "Adipose - Subcutaneous",
                    "Adipose - Visceral (Omentum)",
                    "Adrenal Gland",
                    "Artery - Aorta",
                    "Artery - Coronary",
                    "Artery - Tibial",
                    "Bladder",
                    "Brain - Amygdala",
                    "Brain - Anterior cingulate cortex (BA24)",
                    "Brain - Caudate (basal ganglia)",
                    "Brain - Cerebellar Hemisphere",
                    "Brain - Cerebellum",
                    "Brain - Cortex",
                    "Brain - Frontal Cortex (BA9)",
                    "Brain - Hippocampus",
                    "Brain - Hypothalamus",
                    "Brain - Nucleus accumbens (basal ganglia)",
                    "Brain - Putamen (basal ganglia)",
                    "Brain - Spinal cord (cervical c-1)",
                    "Brain - Substantia nigra",
                    "Breast - Mammary Tissue",
                    "Cells - Cultured fibroblasts",
                    "Cells - EBV-transformed lymphocytes",
                    "Cervix - Ectocervix",
                    "Cervix - Endocervix",
                    "Colon - Sigmoid",
                    "Colon - Transverse",
                    "Esophagus - Gastroesophageal Junction",
                    "Esophagus - Mucosa",
                    "Esophagus - Muscularis",
                    "Fallopian Tube",
                    "Heart - Atrial Appendage",
                    "Heart - Left Ventricle",
                    "Kidney - Cortex",
                    "Liver",
                    "Lung",
                    "Minor Salivary Gland",
                    "Muscle - Skeletal",
                    "Nerve - Tibial",
                    "Ovary",
                    "Pancreas",
                    "Pituitary",
                    "Prostate",
                    "Skin - Not Sun Exposed (Suprapubic)",
                    "Skin - Sun Exposed (Lower leg)",
                    "Small Intestine - Terminal Ileum",
                    "Spleen",
                    "Stomach",
                    "Testis",
                    "Thyroid",
                    "Uterus",
                    "Vagina",
                    "Whole Blood",
                ]
            ]

            # ALTERNATIVE EXONS
            pext_refseq = pext_refseq.loc[pext_refseq["Ratio_num"] < 1].reset_index(
                drop=True
            )

            # TODO : CHECK HOW TO HANDLE UNEXPRESSED AND FULLY EXPRESSED (0 & 1 MEAN PROP)
            ### WARNING FLAG
            # UNEXPRESSED AND FULLY EXPRESSED EXONS
            pext_refseq = pext_refseq.loc[
                (pext_refseq["mean_proportion"] > 0)
                & (pext_refseq["mean_proportion"] < 1)
            ].reset_index(drop=True)

            # HANDLE NAN VALUES
            # pext_refseq[pext_refseq.columns[8:]] = pext_refseq[pext_refseq.columns[8:]].replace(0, np.nan)
            # print(pext_refseq)
            # print(pext_refseq.shape)

            # print(np.count_nonzero(pext_refseq.values == 0))
            # print(pext_refseq[pext_refseq == 0])

            # print(pext_refseq[pext_refseq.columns[8:]].isna().sum().sum())
            pext_refseq = pext_refseq.dropna(
                subset=list(pext_refseq.columns[8:]), how="all"
            )
            # print(pext_refseq.shape)
            # print(np.count_nonzero(pext_refseq.values == 0))
            # print(pext_refseq[pext_refseq.columns[8:]].isna().sum().sum())

            # exit()

            # pext_refseq[pext_refseq.columns[7:]] = pext_refseq[pext_refseq.columns[7:]].fillna(0)

            # EXOTIC SIGMOID FUNCTION APPLIED TO ZSCORE
            def transform(r):
                def sigmoid(x):
                    return 1 / (1 + np.exp(-x))

                vfunc = np.vectorize(sigmoid)
                # print(r)
                if r.isnull().values.any():
                    # print(r)
                    # z_zero = zscore(r.fillna(0))
                    # r_zero = pd.Series(vfunc(z_zero), index=r.index)

                    # z = zscore(r.dropna())
                    # r_na = pd.Series(vfunc(z), index=r.dropna().index)
                    # print(r_zero)
                    # print(r_na)
                    # pd.concat([r, r_zero, r_na], axis=1)
                    return pd.Series(vfunc(zscore(r.dropna())), index=r.dropna().index)

                else:
                    return pd.Series(vfunc(zscore(r)), index=r.index)

            # BUILD MATRIX WITH Z-SCORE + SIGMOID
            # test = z_df = pext_refseq[pext_refseq.columns[8:]].apply(lambda r: zscore(r), axis=1)
            # print(test)
            z_df = pext_refseq[pext_refseq.columns[8:]].apply(
                lambda r: transform(r), axis=1
            )
            # print(z_df)
            # exit()
            z_df = pd.DataFrame.from_dict(dict(zip(z_df.index, z_df.values))).T
            z_df.columns = list(pext_refseq.columns[8:])

            # print(pext_refseq)
            # print(z_df)
            # exit()

            # TMP CONCAT EXOTIC SCORE & PEXT VALUES
            ### WARNING FLAG
            # TODO : modify to delete pext part after defining thresholds
            concat_df = pd.concat(
                [
                    pext_refseq[pext_refseq.columns[:8]].reset_index(drop=True),
                    z_df.add_suffix("_exotic").reset_index(drop=True),
                    pext_refseq[pext_refseq.columns[8:]]
                    .add_suffix("_pextvalue")
                    .reset_index(drop=True),
                ],
                axis=1,
            )

            # REMOVE DUPLICATES
            concat_df = (
                concat_df.drop_duplicates().sort_values(by="MAP").reset_index(drop=True)
            )

            # OUTPUT INTERMEDIATE FILE
            concat_df.to_parquet(self.exotic_files["exotic_path"])
            concat_df.to_excel(
                self.exotic_files["exotic_path"].replace("parquet", "xlsx")
            )
        else:
            concat_df = pd.read_parquet(self.exotic_files["exotic_path"])

        return concat_df

    @staticmethod
    def show_values_on_bars(axs, i=0, fontsize=10, rotation=0):
        def _show_on_single_plot(ax):
            for p in ax.patches:
                print(p)
                _x = p.get_x() + p.get_width() / 2
                _y = p.get_y() + p.get_height()
                if i == 0:
                    value = "{:.0f}".format(p.get_height())
                if i == 2:
                    value = "{:.2f}".format(p.get_height())
                ax.text(
                    _x, _y, value, ha="center", fontsize=fontsize, rotation=rotation
                )

                if i == 3:
                    value = "{:.3f}".format(p.get_height())
                ax.text(
                    _x, _y, value, ha="center", fontsize=fontsize, rotation=rotation
                )

        if isinstance(axs, np.ndarray):
            for idx, ax in np.ndenumerate(axs):
                _show_on_single_plot(ax)
        else:
            _show_on_single_plot(axs)

    def exotic_filtering(self, intermediate_file):

        # FILTER FULLY EXPRESSED & UNEXPRESSED
        intermediate_file = intermediate_file.loc[
            (intermediate_file["mean_proportion"] > 0)
            & (intermediate_file["mean_proportion"] < 1)
            # & (intermediate_file["Ratio_num"] <= 0.5)
        ].reset_index(drop=True)

        # TMP
        # intermediate_file = intermediate_file.tail(1000)

        # FILTER FILES
        exotic = intermediate_file.filter(regex="exotic")
        pext = intermediate_file.filter(regex="pext")

        # RETRIEVING COLUMNS
        exotic["MAP"] = intermediate_file["MAP"]
        exotic["Ratio_num"] = intermediate_file["Ratio_num"]
        exotic["mRNA_nb"] = intermediate_file["mRNA_nb"]
        exotic["symbol"] = intermediate_file["symbol"]

        # ORDER COLUMNS
        exotic = exotic[
            ["symbol", "MAP", "Ratio_num", "mRNA_nb"] + list(exotic.columns[:-4])
        ]

        exotic["tissues_bronze"] = (
            exotic.filter(regex="exotic")
            .fillna(0)
            .apply(
                lambda r: [
                    exotic.columns[j].replace("_exotic", "")
                    for j, c in enumerate(r)
                    if c >= 0.90 and c < 0.95
                ],
                axis=1,
            )
        )
        exotic["tissues_silver"] = (
            exotic.filter(regex="exotic")
            .fillna(0)
            .apply(
                lambda r: [
                    exotic.columns[j].replace("_exotic", "")
                    for j, c in enumerate(r)
                    if c >= 0.95 and c < 0.99
                ],
                axis=1,
            )
        )
        exotic["tissues_gold"] = (
            exotic.filter(regex="exotic")
            .fillna(0)
            .apply(
                lambda r: [
                    exotic.columns[j].replace("_exotic", "")
                    for j, c in enumerate(r)
                    if c >= 0.99
                ],
                axis=1,
            )
        )
        exotic["pext_OK"] = (
            pext.filter(regex="pext")
            .fillna(0)
            .apply(
                lambda r: [
                    pext.columns[j].replace("_pextvalue", "")
                    for j, c in enumerate(r)
                    if c > 0.1
                ],
                axis=1,
            )
        )

        exotic["bronze"] = exotic.apply(
            lambda r: list(set(r["tissues_bronze"]).intersection(set(r["pext_OK"]))),
            axis=1,
        )
        exotic["silver"] = exotic.apply(
            lambda r: list(set(r["tissues_silver"]).intersection(set(r["pext_OK"]))),
            axis=1,
        )
        exotic["gold"] = exotic.apply(
            lambda r: list(set(r["tissues_gold"]).intersection(set(r["pext_OK"]))),
            axis=1,
        )
        exotic["All_thresh"] = exotic.apply(
            lambda r: r["bronze"] + r["silver"] + r["gold"], axis=1
        )

        exotic = exotic.loc[exotic["All_thresh"].str.len() > 0].reset_index(drop=True)

        print(exotic)
        exit()

    def exotic_score_analysis(self, intermediate_file):

        # sns.set_style("white")
        sns.set_context("paper")

        intermediate_file = intermediate_file.drop_duplicates()
        print(intermediate_file)
        exit()

        # # pprint(list(intermediate_file.columns))
        # test = intermediate_file[["symbol", "MAP"]].groupby("symbol").nunique().drop(["symbol"], axis=1).reset_index()
        # print(test)
        # sns.set(font_scale=2)
        # sns.set_style
        # f, ax = plt.subplots(figsize=(15, 15))
        # v = sns.boxplot(data=test, y="MAP", ax=ax)
        # # self.show_values_on_bars(v, 0)
        # # v.set_xticklabels(v.get_xticklabels(), rotation=90)
        # plt.xlabel("")
        # ax.grid(True, axis="x")
        # ax.set_ylim(0, 15)
        # f.tight_layout(rect=[0, 0.05, 1, 1])
        # f.savefig("data/paper/1_EXOTIC/violin_alt_cds.png")

        intermediate_file = intermediate_file.loc[
            (intermediate_file["mean_proportion"] > 0)
            & (intermediate_file["mean_proportion"] < 1)
            # & (intermediate_file["Ratio_num"] <= 0.5)
        ].reset_index(drop=True)
        # intermediate_file = intermediate_file.loc[(intermediate_file["Ratio_num"] > 0.5)].reset_index(drop=True)

        # intermediate_file = intermediate_file.tail(1000)

        pext = intermediate_file.filter(regex="pext")

        exotic = intermediate_file.filter(regex="exotic")
        pext["Max"] = pext.parallel_apply(
            lambda r: max(r.dropna()) if r.dropna().shape[0] > 0 else np.nan, axis=1
        )
        pext["MAP"] = intermediate_file["MAP"]
        pext["symbol"] = intermediate_file["symbol"]

        exotic["Max"] = exotic.parallel_apply(
            lambda r: max(r.dropna()) if r.dropna().shape[0] > 0 else np.nan, axis=1
        )
        exotic["Tissues_Max"] = exotic.parallel_apply(
            lambda r: [
                list(exotic.columns)[j].replace("_exotic", "")
                for j, e in enumerate(r[:-1])
                if e == r["Max"]
            ],
            axis=1,
        )
        exotic["Tissues_Max_organs"] = exotic["Tissues_Max"].apply(
            lambda r: list(set([e.split(" - ")[0] for e in r]))
        )
        exotic["pext_max"] = pext["Max"]

        exotic["MAP"] = intermediate_file["MAP"]
        exotic["Ratio_num"] = intermediate_file["Ratio_num"]
        exotic["mRNA_nb"] = intermediate_file["mRNA_nb"]
        exotic["symbol"] = intermediate_file["symbol"]

        # # exotic = exotic.loc[(exotic["Max"] > 0.99)]

        # # print(exotic[["MAP", "Nerve - Tibial_exotic"]])
        # # print(pext.loc[exotic.index])

        # # exotic["tissues_no_threshold"] = (
        # #     exotic.filter(regex="exotic")
        # #     .fillna(0)
        # #     .apply(
        # #         lambda r: [exotic.columns[j].replace("_exotic", "") for j, c in enumerate(r) if c >= 0.90 and c < 0.95], axis=1
        # #     )
        # # )

        exotic["tissues_OK_bronze"] = (
            exotic.filter(regex="exotic")
            .fillna(0)
            .apply(
                lambda r: [
                    exotic.columns[j].replace("_exotic", "")
                    for j, c in enumerate(r)
                    if c >= 0.90 and c < 0.95
                ],
                axis=1,
            )
        )
        exotic["tissues_OK_silver"] = (
            exotic.filter(regex="exotic")
            .fillna(0)
            .apply(
                lambda r: [
                    exotic.columns[j].replace("_exotic", "")
                    for j, c in enumerate(r)
                    if c >= 0.95 and c < 0.99
                ],
                axis=1,
            )
        )
        exotic["tissues_OK_gold"] = (
            exotic.filter(regex="exotic")
            .fillna(0)
            .apply(
                lambda r: [
                    exotic.columns[j].replace("_exotic", "")
                    for j, c in enumerate(r)
                    if c >= 0.99
                ],
                axis=1,
            )
        )
        exotic["pext_OK"] = (
            pext.filter(regex="pext")
            .fillna(0)
            .apply(
                lambda r: [
                    pext.columns[j].replace("_pextvalue", "")
                    for j, c in enumerate(r)
                    if c > 0.1
                ],
                axis=1,
            )
        )
        # exotic["No_threshold"] = exotic.apply(
        #     lambda r: list(set(r["tissues_no_threshold"]).intersection(set(r["pext_OK"]))), axis=1
        # )
        exotic["OK_bronze"] = exotic.apply(
            lambda r: list(set(r["tissues_OK_bronze"]).intersection(set(r["pext_OK"]))),
            axis=1,
        )
        exotic["OK_silver"] = exotic.apply(
            lambda r: list(set(r["tissues_OK_silver"]).intersection(set(r["pext_OK"]))),
            axis=1,
        )
        exotic["OK_gold"] = exotic.apply(
            lambda r: list(set(r["tissues_OK_gold"]).intersection(set(r["pext_OK"]))),
            axis=1,
        )
        exotic["OK"] = exotic.apply(
            lambda r: r["OK_bronze"] + r["OK_silver"] + r["OK_gold"], axis=1
        )

        # print(
        #     pd.concat(
        #         [
        #             # exotic["No_threshold"].explode().value_counts(),
        #             exotic["OK_bronze"].explode().value_counts(),
        #             exotic["OK_silver"].explode().value_counts(),
        #             exotic["OK_gold"].explode().value_counts(),
        #         ],
        #         axis=1,
        #     )
        # )

        t = exotic[["MAP", "symbol", "OK"]]
        t = t.loc[t["OK"].str.len() > 0].explode("OK")
        t = t.groupby(["OK"]).nunique().drop(["OK"], axis=1)
        t.columns = ["Nb of exons", "Nb of genes"]

        sns.set_context("poster")
        f, ax = plt.subplots(figsize=(35, 15))
        v = t.plot.bar(ax=ax)
        self.show_values_on_bars(v, 0)
        v.set_xticklabels(v.get_xticklabels(), rotation=90)
        plt.xlabel("")
        ax.grid(True, axis="y")
        f.tight_layout(rect=[0, 0.05, 1, 1])
        f.savefig("data/paper/1_EXOTIC/barplot_tissues_genes_exons.png")

        exit()

        # bins = [0, 0.2, 0.4, 0.6, 0.8, 1]
        # labels = bins.copy()
        # labels_ratio = [str(round(labels[j], 1)) + " - " + str(round(labels[j + 1], 1)) for j in range(len(labels) - 1)]
        # exotic["Ratio_num_bins"] = pd.cut(exotic["Ratio_num"], bins=bins, labels=labels_ratio, include_lowest=True)

        # bins = [0, 0.2, 0.4, 0.6, 0.8, 1]
        # labels = bins.copy()
        # labels_ratio = [str(round(labels[j], 1)) + " - " + str(round(labels[j + 1], 1)) for j in range(len(labels) - 1)]
        # exotic["pext_max_bins"] = pd.cut(exotic["pext_max"], bins=bins, labels=labels_ratio, include_lowest=True)

        # bins = [0, 0.9, 0.95, 0.99, 1]
        # labels = bins.copy()
        # labels_ratio = [str(round(labels[j], 3)) + " - " + str(round(labels[j + 1], 3)) for j in range(len(labels) - 1)]
        # exotic["Max_bins"] = pd.cut(exotic["Max"], bins=bins, labels=labels_ratio, include_lowest=True)

        ## BARPLOT EXOTIC COUNT CUTOFFS

        print(exotic)

        # tmp_barplot = exotic.groupby("Max_bins").size().reset_index()
        # tmp_barplot.columns = ["variable", "value"]
        # tmp_barplot["variable"] = tmp_barplot["variable"].cat.rename_categories(
        #     {"0 - 0.9": "No threshold", "0.9 - 0.95": "Bronze", "0.95 - 0.99": "Silver", "0.99 - 1": "Gold"}
        # )
        # tmp_barplot["variable"] = tmp_barplot["variable"].cat.remove_categories(["0 - 0.8"])
        # tmp_barplot = tmp_barplot.loc[tmp_barplot["variable"].isna() == False]
        # tmp_barplot["variable"] = tmp_barplot["variable"].cat.add_categories(["Total"])
        # tmp_barplot.loc[0] = ["Total", exotic.shape[0]]
        # tmp_barplot.sort_index(inplace=True)
        # print(tmp_barplot)

        # sns.set(font_scale=2)
        # sns.set_context("poster")
        # f, ax = plt.subplots(figsize=(10, 10))
        # v = sns.barplot(data=tmp_barplot, x="variable", y="value", ax=ax, palette=["#2d3436", "#A77044", "#D7D7D7", "#FEE101"])
        # self.show_values_on_bars(v, 0, 20)
        # v.set_xticklabels(v.get_xticklabels(), rotation=45)
        # plt.xlabel("")
        # plt.ylabel("Nb of exons")
        # ax.grid(True, axis="y")
        # f.tight_layout(rect=[0, 0.05, 1, 1])
        # f.savefig("data/paper/1_EXOTIC/barplot_test_exotic_bins.png")
        # exit()

        # print(exotic)
        # print(exotic.groupby("Max_bins").size())
        # exit()

        # bronze = exotic["OK_bronze"].explode().dropna().value_counts()
        # silver = exotic["OK_silver"].explode().dropna().value_counts()
        # gold = exotic["OK_gold"].explode().dropna().value_counts()

        # tmp_barplot_tissues = (
        #     pd.concat([bronze, silver, gold], axis=1)
        #     .reset_index()
        #     .melt(id_vars="index", value_vars=["OK_bronze", "OK_silver", "OK_gold"])
        # )
        # tmp_barplot_tissues.columns = ["Tissue", "Threshold", "Nb of exons"]
        # tmp_barplot_tissues["Threshold"] = tmp_barplot_tissues["Threshold"].str.replace("OK_", "")
        # print(tmp_barplot_tissues)
        # exit()

        # sns.set_context("poster")
        # f, ax = plt.subplots(figsize=(35, 15))
        # v = sns.barplot(
        #     data=tmp_barplot_tissues,
        #     x="Tissue",
        #     y="Nb of exons",
        #     hue="Threshold",
        #     ax=ax,
        #     palette=["#A77044", "#D7D7D7", "#FEE101"],
        # )
        # self.show_values_on_bars(v, 0, 12, 90)
        # v.set_xticklabels(v.get_xticklabels(), rotation=90)
        # plt.xlabel("")
        # ax.grid(True, axis="y")
        # f.tight_layout(rect=[0, 0.05, 1, 1])
        # f.savefig("data/paper/1_EXOTIC/barplot_test_exotic_bins_tissues_05.png")
        ## HEATMAP

        # exotic["OK_organs"] = exotic["OK"].apply(lambda r: list(set([e.split(" - ")[0] for e in r])))
        # exotic["Len_tissues"] = exotic["OK"].apply(lambda r: len(r))

        # t = exotic.explode("OK_organs")
        # l = exotic["Tissues_Max_organs"].values.tolist()
        # l = [list(set([e.split(" - ")[0] for e in r])) for r in l if len(r) > 1]
        # l = [e for e in l if len(e) <= 2]
        # l = [[e[0], e[0]] if len(e) == 1 else e for e in l]
        # new_l = list()
        # for e in l:
        #     if e[0] == e[1]:
        #         new_l.append({"Tissue_0": e[0], "Tissue_1": e[1]})
        #     else:
        #         new_l.append({"Tissue_0": e[0], "Tissue_1": e[1]})
        #         # new_l.append({"Tissue_0": e[1], "Tissue_1": e[0]})

        # tmp_df_heatmap = pd.DataFrame(new_l).groupby(["Tissue_0", "Tissue_1"]).size()
        # tmp_df_heatmap = tmp_df_heatmap.reset_index()
        # tmp_df_heatmap.columns = ["Tissue_0", "Tissue_1", "value"]
        # print(tmp_df_heatmap.sort_values("value", ascending=False))

        def median_pext(r):
            try:
                return np.median(
                    [
                        e
                        for e in r["pext_values"]
                        if math.isnan(e) is False and e != r["pext_max"]
                    ]
                )
            except:
                print(r)
                exit()

        exotic["pext_max"] = pext["Max"]
        exotic["pext_values"] = pext[pext.columns[:-3]].values.tolist()
        exotic["pext_median_except_max"] = exotic.apply(
            lambda r: median_pext(r), axis=1
        )
        exotic["pext_shift"] = exotic["pext_max"] - exotic["pext_median_except_max"]
        # print(exotic)
        # tmp_barplot = exotic[
        #     [
        #         "Ratio_num_bins",
        #         "pext_max_bins",
        #         "Max_bins",
        #         "pext_max",
        #         "Max",
        #     ]
        # ]
        # tmp_barplot = tmp_barplot.rename({"Ratio_num_bins": "Exon freq"}, axis=1)
        # ].melt(
        #     id_vars=[
        #         "Ratio_num_bins",
        #         "pext_max_bins",
        #         "Max_bins",
        #         "pext_max",
        #     ],
        #     value_vars=[
        #         "pext_median_except_max",
        #         "pext_shift",
        #         # "pext_max",
        #     ],
        # )
        # CLEANING
        exotic = exotic.dropna(subset=["Max"])
        exotic = exotic.loc[exotic["pext_max"] > 0.1]

        # sns.set_context("poster")
        # f, ax = plt.subplots(figsize=(15, 15))
        # grid = sns.FacetGrid(
        #     tmp_barplot,
        #     row="Exon freq",
        #     hue="Exon freq",
        #     # row="Ratio_num_bins",
        #     palette="coolwarm",
        #     # margin_titles=True,
        #     height=4,
        #     aspect=3,
        # )
        # grid.map(sns.violinplot, "Max", linewitdh=2)
        # # self.show_values_on_bars(v, 0)
        # # v.set_xticklabels(v.get_xticklabels(), rotation=90)
        # # plt.xlabel("")
        # # ax.grid(True, axis="x")
        # for x in grid.axes.flat:
        #     x.grid(True, axis="x")
        #     x.set_xlim(0.5, 1)
        # grid.fig.tight_layout()
        # grid.fig.savefig("data/paper/1_EXOTIC/facetgrid_test_ratio_exotic.png")

        # exit()
        # # exotic["MAP"] = intermediate_file["MAP"]

        # # Binning EXOTIC
        # # bins = [0, 0.5, 0.8, 0.9, 0.95, 0.99, 0.995, 0.999, 1]
        # bins = [0, 0.9, 0.95, 0.99, 1]
        # labels = bins.copy()
        # labels_ratio = [str(round(labels[j], 3)) + " - " + str(round(labels[j + 1], 3)) for j in range(len(labels) - 1)]
        # exotic["Max_bins"] = pd.cut(exotic["Max"], bins=bins, labels=labels_ratio, include_lowest=True)

        ## BARPLOT TISSUES CUTOFFS

        # tmp_barplot_tissues = (
        #     exotic[["Max", "Tissues_Max", "pext_max", "Max_bins"]]
        #     .explode("Tissues_Max")
        #     .groupby(["Tissues_Max", "Max_bins"])
        #     .size()
        #     .reset_index()
        # )
        # tmp_barplot_tissues.columns = ["Tissue", "Bin", "value"]
        # tmp_barplot_tissues["Bin"] = tmp_barplot_tissues["Bin"].cat.remove_categories(["0 - 0.9"])

        # # print(tmp_barplot_tissues)

        # sns.set(font_scale=2)
        # f, ax = plt.subplots(figsize=(35, 15))
        # v = sns.barplot(data=tmp_barplot_tissues, x="Tissue", y="value", hue="Bin", ax=ax)
        # self.show_values_on_bars(v, 0)
        # v.set_xticklabels(v.get_xticklabels(), rotation=90)
        # plt.xlabel("")
        # ax.grid(True, axis="x")
        # f.tight_layout(rect=[0, 0.05, 1, 1])
        # f.savefig("data/paper/1_EXOTIC/barplot_test_exotic_bins_tissues_05_01.png")
        # exit()

        ## BARPLOT TISSUES TOTAL

        # tmp_barplot_tissues = (
        #     exotic[["Tissues_Max", "symbol", "MAP"]].explode("Tissues_Max").groupby(["Tissues_Max", "symbol"]).nunique()
        # )
        # # tmp_barplot_tissues["Ratio"] = tmp_barplot_tissues["MAP"] / tmp_barplot_tissues["symbol"]
        # tmp_barplot_tissues.columns = ["Remove", "Genes", "Exons"]
        # tmp_barplot_tissues = tmp_barplot_tissues.drop(["Remove", "Genes"], axis=1).reset_index()
        # print(tmp_barplot_tissues)

        # f, ax = plt.subplots(figsize=(30, 15))
        # v = sns.violinplot(data=tmp_barplot_tissues, x="Tissues_Max", y="Exons", linewidth=2)
        # self.show_values_on_bars(v, 2)
        # v.set_xticklabels(v.get_xticklabels(), rotation=90)
        # plt.xlabel("")
        # ax.grid(True, axis="x")
        # ax.set_ylim(0, 10)

        # f.tight_layout(rect=[0, 0.05, 1, 1])
        # f.savefig("data/paper/1_EXOTIC/violinplot_test_exotic_tissues_genes_exon.png")
        # exit()

        # tmp_barplot_tissues = (
        #     tmp_barplot_tissues.drop(["Remove"], axis=1).reset_index().melt(id_vars="Tissues_Max", value_vars=["Genes", "Exons"])
        # )

        # mask = tmp_barplot_tissues.variable.isin(["Genes", "Exons"])
        # scale = float(tmp_barplot_tissues[~mask].value.mean() / tmp_barplot_tissues[mask].value.mean())
        # tmp_barplot_tissues.loc[mask, "value"] = tmp_barplot_tissues.loc[mask, "value"] * scale
        # print(mask)

        # tmp_barplot_tissues.columns = ["Tissue", "value"]
        # tmp_barplot_tissues["Bin"] = tmp_barplot_tissues["Bin"].cat.remove_categories(["0 - 0.9"])

        # print(tmp_barplot_tissues)

        # sns.set(font_scale=2)
        # f, ax1 = plt.subplots(figsize=(30, 15))
        # v = sns.barplot(data=tmp_barplot_tissues, x="Tissues_Max", y="value", hue="variable", ax=ax1)

        # # Create a second y-axis with the scaled ticks
        # # ax1.set_ylabel("X and Y")
        # # ax2 = ax1.twinx()

        # # Ensure ticks occur at the same positions, then modify labels
        # # ax2.set_ylim(ax1.get_ylim())
        # # ax2.set_yticklabels(np.round(ax1.get_yticks() / scale, 1))
        # # ax2.set_ylabel("A and B")

        # self.show_values_on_bars(v, 0)
        # v.set_xticklabels(v.get_xticklabels(), rotation=90)
        # plt.xlabel("")

        # ax1.grid(True, axis="x")
        # f.tight_layout(rect=[0, 0.05, 1, 1])
        # f.savefig("data/paper/1_EXOTIC/barplot_test_exotic_tissues.png")

        cutoffs = [0.9, 0.95, 0.99]
        l_df = list()
        for j, c in enumerate(cutoffs):
            print(j, len(cutoffs), c)
            if j < (len(cutoffs) - 1):
                tmp = exotic.loc[
                    (exotic["Max"] >= c) & (exotic["Max"] < cutoffs[j + 1]),
                    ["Max", "pext_max", "pext_median_except_max", "pext_shift"],
                ]
            else:
                tmp = exotic.loc[
                    exotic["Max"] >= c,
                    ["Max", "pext_max", "pext_median_except_max", "pext_shift"],
                ]
            tmp["Threshold"] = c
            l_df.append(tmp)
        tmp_df = pd.concat(l_df)
        tmp_df = tmp_df.rename(
            {
                "pext_max": "Max expression ratio",
                "pext_median_except_max": "Residual median expression ratio",
                "pext_shift": "Gap",
            },
            axis=1,
        )
        print(tmp_df)
        tmp_df = tmp_df.melt(
            id_vars="Threshold",
            value_vars=[
                "Max expression ratio",
                "Residual median expression ratio",
                "Gap",
            ],
        )
        print(tmp_df)

        ## BARPLOT PEXT SHIFT
        print(tmp_df.groupby(["Threshold", "variable"]).describe())

        sns.set_context("poster")
        f, ax = plt.subplots(figsize=(15, 15))
        v = sns.barplot(
            data=tmp_df,
            x="variable",
            y="value",
            hue="Threshold",
            ax=ax,
            palette=["#A77044", "#D7D7D7", "#FEE101"],
        )
        v.set_xticklabels(v.get_xticklabels(), rotation=45)
        plt.xlabel("")
        plt.ylabel("")
        ax.grid(True, axis="y")
        f.tight_layout(rect=[0, 0.05, 1, 1])
        f.savefig("data/paper/1_EXOTIC/barplot_test_pext_shift.png")
        exit()

        # exotic["Min"] = exotic.apply(lambda r: min(r), axis=1)
        # exotic["Median"] = exotic.apply(lambda r: np.median(r), axis=1)
        # exotic["Mean"] = exotic.apply(lambda r: np.mean(r), axis=1)

        # plot_basic_stats_exotic = exotic[["Min", "Max", "Median", "Mean"]].melt(value_vars=["Min", "Max", "Median", "Mean"])

        # f, ax = plt.subplots(figsize=(10, 10))
        # v = sns.violinplot(data=plot_basic_stats_exotic, x="variable", y="value", ax=ax, linewidth=3)
        # # v.set_xticklabels(v.get_xticklabels(), rotation=90)
        # # plt.xlabel("")
        # ax.grid(True, axis="y")
        # # f.tight_layout(rect=[0, 0.05, 1, 1])
        # f.savefig("data/paper/1_EXOTIC/violin_stats.png")

        # print(exotic)
        # print(exotic.describe())

        # print(intermediate_file.loc[intermediate_file["Bladder_exotic"] > 0.99])
        # print(intermediate_file.loc[intermediate_file["Bladder_exotic"] > 0.99].filter(regex="exotic"))

        exotic = exotic.melt(id_vars=["MAP"], value_vars=exotic.columns[:-1])
        exotic["variable"] = exotic["variable"].str.replace("_exotic", "")
        pext = pext.melt(id_vars=["MAP", "symbol"], value_vars=pext.columns[:-2])
        pext["variable"] = pext["variable"].str.replace("_pextvalue", "")

        def apply_pext_comparison(r):
            median_pext_value_other_tissues = float(
                (
                    intermediate_file.loc[intermediate_file["MAP"] == r["MAP"]]
                    .filter(regex="pext")
                    .T.dropna()
                    .drop(r["variable"] + "_pextvalue")
                    .median()
                )
            )
            fold = float(r["value"] / median_pext_value_other_tissues)
            r["median_others"] = median_pext_value_other_tissues
            r["fold"] = fold
            return r

        l = list()
        l_pext_values = list()
        values = [0.8, 0.9, 0.95, 0.99]
        # values = [0.999]
        for j, value in enumerate(values):
            if j < len(values) - 1:
                tmp_exotic = exotic.loc[
                    (exotic["value"] >= value) & (exotic["value"] < values[j + 1])
                ]
            else:
                tmp_exotic = exotic.loc[(exotic["value"] >= value)]
            tmp_pext = pext.loc[tmp_exotic.index]

            # print(tmp_pext.parallel_apply(apply_pext_comparison, axis=1))
            tmp_pext["hue"] = value
            # print(tmp_pext)
            # exit

            # tmp_pext_pivot = pd.pivot(tmp_pext, index="MAP", columns="variable", values="value")
            # tmp_pext_others = intermediate_file.loc[intermediate_file["MAP"].isin(list(tmp_pext.MAP.unique()))]
            # print(tmp_pext_pivot)
            # print(tmp_pext_others)
            # exit()
            # low = tmp_pext.loc[tmp_pext["value"] < 0.1]
            # for i, row in low.iterrows():
            # print(row)
            # print(pext.loc[pext["symbol"] == row["symbol"]].pivot(index="MAP", columns="variable", values="value"))
            # print(pext.loc[pext["symbol"].isin(low.symbol.unique().tolist())])
            # exit()
            # print(tmp_exotic)
            # print(tmp_pext)
            l_pext_values.append(tmp_pext)
            tmp_exotic = tmp_exotic.groupby("variable").size()
            l.append(tmp_exotic)
        plot_exotic = pd.concat(l, axis=1).reset_index()
        plot_exotic.columns = ["Tissue"] + values
        plot_exotic = plot_exotic.melt(id_vars="Tissue", value_vars=values)
        # plot_exotic = plot_exotic.loc[~plot_exotic["Tissue"].isin(["Min", "Max", "Median", "Mean"])]
        plot_exotic_pext = pd.concat(l_pext_values, axis=0)

        plt.style.use("ggplot")

        sns.set(font_scale=2)

        f, ax = plt.subplots(figsize=(35, 15))
        v = sns.barplot(data=plot_exotic, x="variable", y="value", ax=ax)
        v.set_xticklabels(v.get_xticklabels(), rotation=90)
        plt.xlabel("")
        ax.grid(True, axis="x")
        f.tight_layout(rect=[0, 0.05, 1, 1])
        f.savefig("data/paper/1_EXOTIC/barplot_test.png")

        sns.set(font_scale=1)

        f, ax = plt.subplots(figsize=(45, 15))
        v = sns.boxplot(data=plot_exotic_pext, x="hue", y="value", ax=ax)
        v.set_xticklabels(v.get_xticklabels(), rotation=90)
        plt.xlabel("")
        ax.grid(True, axis="x")
        # f.tight_layout(rect=[0, 0.05, 1, 1])
        f.savefig("data/paper/1_EXOTIC/violin_pext.png")

    def heatmap_associations_between_tissues(self, concat_df, glob=False):

        l = concat_df["Tissues"].values.tolist()

        if glob is True:
            l = [list(set([e.split(" - ")[0] for e in r])) for r in l if len(r) > 1]
        else:
            l = [list(set([e for e in r])) for r in l if len(r) > 1]

        l = [e for e in l if len(e) <= 2]
        l = [[e[0], e[0]] if len(e) == 1 else e for e in l]
        new_l = list()
        for e in l:
            if e[0] == e[1]:
                new_l.append({"Tissue_0": e[0], "Tissue_1": e[1]})
            else:
                new_l.append({"Tissue_0": e[0], "Tissue_1": e[1]})
                new_l.append({"Tissue_0": e[1], "Tissue_1": e[0]})

        tmp_df_heatmap = pd.DataFrame(new_l).groupby(["Tissue_0", "Tissue_1"]).size()
        if glob is True:
            t = "global"
        else:
            t = "detailed"
        tmp_df_heatmap.to_excel(
            self.output_dir
            + self.subdirs_output_dir[1]
            + "table_heatmap_{}.xlsx".format(t)
        )
        # print(tmp_df_heatmap)
        tmp_df_heatmap = tmp_df_heatmap.reset_index()

        tmp_df_heatmap.columns = ["Tissue_0", "Tissue_1", "value"]
        tmp_df_heatmap["value"] = tmp_df_heatmap["value"].astype(int)
        tmp_df_heatmap = tmp_df_heatmap.pivot("Tissue_0", "Tissue_1", "value")
        tmp_df_heatmap = tmp_df_heatmap.fillna(0)
        tmp_df_heatmap = tmp_df_heatmap.astype(int)
        if glob is True:
            gtex_cols = [
                "Adipose",
                "Adrenal Gland",
                "Artery",
                "Bladder",
                "Brain",
                "Breast",
                "Cells",
                "Cervix",
                "Colon",
                "Esophagus",
                "Fallopian Tube",
                "Heart",
                "Kidney",
                "Liver",
                "Lung",
                "Minor Salivary Gland",
                "Muscle",
                "Nerve",
                "Ovary",
                "Pancreas",
                "Pituitary",
                "Prostate",
                "Skin",
                "Small Intestine",
                "Spleen",
                "Stomach",
                "Testis",
                "Thyroid",
                "Uterus",
                "Vagina",
                "Whole Blood",
            ]
        else:
            gtex_cols = [
                "Adipose - Subcutaneous",
                "Adipose - Visceral (Omentum)",
                "Adrenal Gland",
                "Artery - Aorta",
                "Artery - Coronary",
                "Artery - Tibial",
                "Bladder",
                "Brain - Amygdala",
                "Brain - Anterior cingulate cortex (BA24)",
                "Brain - Caudate (basal ganglia)",
                "Brain - Cerebellar Hemisphere",
                "Brain - Cerebellum",
                "Brain - Cortex",
                "Brain - Frontal Cortex (BA9)",
                "Brain - Hippocampus",
                "Brain - Hypothalamus",
                "Brain - Nucleus accumbens (basal ganglia)",
                "Brain - Putamen (basal ganglia)",
                "Brain - Spinal cord (cervical c-1)",
                "Brain - Substantia nigra",
                "Breast - Mammary Tissue",
                "Cells - Cultured fibroblasts",
                "Cells - EBV-transformed lymphocytes",
                "Cervix - Ectocervix",
                "Cervix - Endocervix",
                "Colon - Sigmoid",
                "Colon - Transverse",
                "Esophagus - Gastroesophageal Junction",
                "Esophagus - Mucosa",
                "Esophagus - Muscularis",
                "Fallopian Tube",
                "Heart - Atrial Appendage",
                "Heart - Left Ventricle",
                "Kidney - Cortex",
                "Kidney - Medulla",
                "Liver",
                "Lung",
                "Minor Salivary Gland",
                "Muscle - Skeletal",
                "Nerve - Tibial",
                "Ovary",
                "Pancreas",
                "Pituitary",
                "Prostate",
                "Skin - Sun Exposed (Lower leg)",
                "Skin - Not Sun Exposed (Suprapubic)",
                "Small Intestine - Terminal Ileum",
                "Spleen",
                "Stomach",
                "Testis",
                "Thyroid",
                "Uterus",
                "Vagina",
                "Whole Blood",
            ]

        for col in gtex_cols:
            if col not in tmp_df_heatmap.columns:
                tmp_df_heatmap[col] = 0
            if col not in tmp_df_heatmap.index:
                tmp_df_heatmap.loc[col] = 0

        tmp_df_heatmap = tmp_df_heatmap[gtex_cols].sort_index()
        print(tmp_df_heatmap)

        f, ax = plt.subplots(figsize=(15, 15))
        if glob is True:
            boundaries = [0, 1, 5, 10, 20, 50, 100, 300]
            cmap = [
                "#ecf0f1",
                "#f1c40f",
                "#f39c12",
                "#e67e22",
                "#d35400",
                "#e74c3c",
                "#c0392b",
            ]
        else:
            boundaries = [0, 1, 5, 10, 20, 50, 100]
            cmap = ["#ecf0f1", "#f1c40f", "#f39c12", "#e67e22", "#d35400", "#c0392b"]
        norm = matplotlib.colors.BoundaryNorm(boundaries, len(cmap))
        vmin, vmax = tmp_df_heatmap.values.min(), tmp_df_heatmap.values.max()

        hmp = sns.heatmap(
            tmp_df_heatmap,
            annot=True,
            linewidth=0.5,
            fmt="d",
            ax=ax,
            cmap=cmap,
            norm=norm,
            vmin=vmin,
            vmax=vmax,
            cbar_kws={"ticks": boundaries},
        )

        if glob is False:
            hmp.hlines(
                [
                    2,
                    3,
                    6,
                    7,
                    20,
                    21,
                    23,
                    25,
                    27,
                    30,
                    31,
                    33,
                    35,
                    36,
                    37,
                    38,
                    39,
                    40,
                    41,
                    42,
                    43,
                    44,
                    46,
                    47,
                    48,
                    49,
                    50,
                    51,
                    52,
                    53,
                    54,
                ],
                *ax.get_xlim()
            )
            hmp.vlines(
                [
                    2,
                    3,
                    6,
                    7,
                    20,
                    21,
                    23,
                    25,
                    27,
                    30,
                    31,
                    33,
                    35,
                    36,
                    37,
                    38,
                    39,
                    40,
                    41,
                    42,
                    43,
                    44,
                    46,
                    47,
                    48,
                    49,
                    50,
                    51,
                    52,
                    53,
                    54,
                ],
                *ax.get_ylim()
            )

        f.tight_layout()
        if glob is True:
            f.savefig(
                self.output_dir + self.subdirs_output_dir[1] + "heatmap_global.png"
            )
        else:
            f.savefig(
                self.output_dir + self.subdirs_output_dir[1] + "heatmap_detailed.png"
            )

    def barplot_gene_exon_nb(self, concat_df):
        for glob in [True, False]:
            df = concat_df[["symbol", "MAP", "Medium", "High", "Very-High", "Tissues"]]
            if glob is True:

                df["Supertissues"] = df.Tissues.apply(
                    lambda r: list([e.split(" - ")[0] for e in list(r)])
                )
            else:
                df["Supertissues"] = df["Tissues"]

            df = df[["symbol", "MAP", "Supertissues"]].explode("Supertissues")
            if glob is True:
                tmp_df = pd.concat(
                    [
                        df.groupby("Supertissues").nunique(),
                        df.groupby("Supertissues").size(),
                    ],
                    axis=1,
                )
            else:
                tmp_df = df.groupby("Supertissues").nunique()

            tmp_df = tmp_df.drop("Supertissues", axis=1)
            tmp_df = tmp_df.rename(
                {"symbol": "Genes nb", "MAP": "Exons nb", 0: "Associations nb"}, axis=1
            )

            if glob is True:
                t = "global"
            else:
                t = "detailed"
            tmp_df.to_excel(
                self.output_dir
                + self.subdirs_output_dir[1]
                + "table_barplot_{}.xlsx".format(t)
            )
            plt.style.use("ggplot")
            sns.set_context("paper")
            f, ax = plt.subplots(figsize=(15, 7))
            if glob is True:
                tmp_df[["Genes nb", "Exons nb", "Associations nb"]].plot.bar(ax=ax)

            else:
                tmp_df[["Genes nb", "Exons nb"]].plot.bar(ax=ax)

            for p in ax.patches:
                ax.annotate(
                    str(p.get_height()),
                    (p.get_x() - 0.001, p.get_height() + 1.2),
                    fontsize=6,
                    rotation=90,
                )
            ax.grid(axis="x")
            plt.xlabel("")
            if glob is True:
                plt.suptitle(
                    "Distribution of genes, exons and total associations (for tissue with more than one subtissue) by tissue"
                )
                plt.legend(ncol=3)
            else:
                plt.suptitle("Distribution of genes, exons by subtissue")
                plt.legend(ncol=2)
            f.tight_layout(rect=[0, 0, 1, 0.95])
            f.savefig(
                self.output_dir
                + self.subdirs_output_dir[1]
                + "barplot_tissues_{}.png".format(t),
                dpi=300,
            )

    def barplot_stringency_nb(
        self,
        concat_df,
    ):
        for glob in [True, False]:
            df = concat_df[["symbol", "MAP", "Medium", "High", "Very-High", "Tissues"]]
            print(df)
            if glob is True:
                for col in ["Medium", "High", "Very-High"]:
                    df[col] = df[col].apply(
                        lambda r: list([e.split(" - ")[0] for e in list(r)])
                    )
            else:
                pass
            df = df.melt(
                id_vars=["symbol", "MAP"], value_vars=["Medium", "High", "Very-High"]
            )
            df = df.explode("value")
            df = (
                df[["variable", "value"]]
                .groupby(["variable", "value"])
                .size()
                .reset_index()
                .pivot(index="value", columns="variable", values=0)
            )
            df = df.fillna(0)
            df = df.astype(int)

            if glob is True:
                t = "global"
            else:
                t = "detailed"

            df = df[["Medium", "High", "Very-High"]]

            df.to_excel(
                self.output_dir
                + self.subdirs_output_dir[1]
                + "table_barplot_stringency_{}.xlsx".format(t)
            )

            plt.style.use("ggplot")
            sns.set_context("paper")
            f, ax = plt.subplots(figsize=(15, 7))
            df.plot.bar(ax=ax)
            for p in ax.patches:
                ax.annotate(
                    str(p.get_height()),
                    (p.get_x() - 0.002, p.get_height() + 1.2),
                    fontsize=6,
                    rotation=90,
                )
            ax.grid(axis="x")
            plt.xlabel("")
            if glob is True:
                plt.suptitle(
                    "Distribution of genes, exons and total associations (for tissue with more than one subtissue) by tissue"
                )

            else:
                plt.suptitle("Distribution of genes, exons by subtissue")
            plt.legend(ncol=3)
            f.tight_layout(rect=[0, 0, 1, 0.95])
            f.savefig(
                self.output_dir
                + self.subdirs_output_dir[1]
                + "barplot_tissues_stringency_{}.png".format(t),
                dpi=300,
            )

    def run_config(
        self, config, intermediate_file, output_dir, disease_genes_multi_iso=list()
    ):
        for rare in [True]:
            # m = multiprocessing.Manager()
            # l = list()
            # parmap.starmap(self.mp_pext_search_exon_of_interest, list(zip(tmp_genes)),
            #    intermediate_file, l, config, rare, pm_pbar=True, pm_processes=self.cpus)
            concat_df = self.mp_pext_search_exon_of_interest(
                intermediate_file, config, rare, output_dir, disease_genes_multi_iso
            )
            # self.heatmap_associations_between_tissues(concat_df)
            # self.barplot_stringency_nb(concat_df)
            # self.barplot_gene_exon_nb(concat_df)
            # exit()

            concat_df = self.compare_exons_gtex_profile(concat_df, exon_usage=True)

            # print(concat_df)
            # print(concat_df.MAP.nunique())
            # print(concat_df.symbol.nunique())
            # concat_df = concat_df.loc[concat_df["symbol"].isin(genes_cardio)]
            # concat_df.to_excel(output_dir + "genes_pext.xlsx")
            return concat_df

    @staticmethod
    def mp_pext_search_exon_of_interest(
        df, config_search, rare, output_dir, disease_genes_multi_iso=list()
    ):
        # pd.options.display.max_columns = 10
        # pd.reset_option('all')

        print("H genes", str(df.symbol.nunique()))
        print("H exons", str(df.MAP.nunique()))

        convert_tissue_dict = {
            "Adipose_Subcutaneous": "Adipose - Subcutaneous",
            "Adipose_Visceral_Omentum_": "Adipose - Visceral (Omentum)",
            "AdrenalGland": "Adrenal Gland",
            "Artery_Aorta": "Artery - Aorta",
            "Artery_Coronary": "Artery - Coronary",
            "Artery_Tibial": "Artery - Tibial",
            "Bladder": "Bladder",
            "Brain_Amygdala": "Brain - Amygdala",
            "Brain_Anteriorcingulatecortex_BA24_": "Brain - Anterior cingulate cortex (BA24)",
            "Brain_Caudate_basalganglia_": "Brain - Caudate (basal ganglia)",
            "Brain_CerebellarHemisphere": "Brain - Cerebellar Hemisphere",
            "Brain_Cerebellum": "Brain - Cerebellum",
            "Brain_Cortex": "Brain - Cortex",
            "Brain_FrontalCortex_BA9_": "Brain - Frontal Cortex (BA9)",
            "Brain_Hippocampus": "Brain - Hippocampus",
            "Brain_Hypothalamus": "Brain - Hypothalamus",
            "Brain_Nucleusaccumbens_basalganglia_": "Brain - Nucleus accumbens (basal ganglia)",
            "Brain_Putamen_basalganglia_": "Brain - Putamen (basal ganglia)",
            "Brain_Spinalcord_cervicalc_1_": "Brain - Spinal cord (cervical c-1)",
            "Brain_Substantianigra": "Brain - Substantia nigra",
            "Breast_MammaryTissue": "Breast - Mammary Tissue",
            "Cells_EBV_transformedlymphocytes": "Cells - EBV-transformed lymphocytes",
            "Cervix_Ectocervix": "Cervix - Ectocervix",
            "Cervix_Endocervix": "Cervix - Endocervix",
            "Colon_Sigmoid": "Colon - Sigmoid",
            "Colon_Transverse": "Colon - Transverse",
            "Esophagus_GastroesophagealJunction": "Esophagus - Gastroesophageal Junction",
            "Esophagus_Mucosa": "Esophagus - Mucosa",
            "Esophagus_Muscularis": "Esophagus - Muscularis",
            "FallopianTube": "Fallopian Tube",
            "Heart_AtrialAppendage": "Heart - Atrial Appendage",
            "Heart_LeftVentricle": "Heart - Left Ventricle",
            "Kidney_Cortex": "Kidney - Cortex",
            "Liver": "Liver",
            "Lung": "Lung",
            "MinorSalivaryGland": "Minor Salivary Gland",
            "Muscle_Skeletal": "Muscle - Skeletal",
            "Nerve_Tibial": "Nerve - Tibial",
            "Ovary": "Ovary",
            "Pancreas": "Pancreas",
            "Pituitary": "Pituitary",
            "Prostate": "Prostate",
            "Skin_NotSunExposed_Suprapubic_": "Skin - Not Sun Exposed (Suprapubic)",
            "Skin_SunExposed_Lowerleg_": "Skin - Sun Exposed (Lower leg)",
            "SmallIntestine_TerminalIleum": "Small Intestine - Terminal Ileum",
            "Spleen": "Spleen",
            "Stomach": "Stomach",
            "Testis": "Testis",
            "Thyroid": "Thyroid",
            "Uterus": "Uterus",
            "Vagina": "Vagina",
            "WholeBlood": "Whole Blood",
        }

        df = df.sort_values(by=["MAP", "variable"])

        # RARE EXONS
        gene_df = df.loc[
            (df["Ratio_num"] <= config_search["Presence"]["Exon_frequency_cutoff"])
        ].reset_index(drop=True)

        # print(gene_df)
        print(gene_df.symbol.nunique())
        print(gene_df.MAP.nunique())
        # WEAK GLOBAL EXPRESSION
        gene_df = gene_df.loc[
            (
                gene_df["mean_proportion"]
                < config_search["Presence"]["Global_max_expression"]
            )
        ].reset_index(drop=True)
        # print(gene_df)
        print(gene_df.symbol.nunique())
        print(gene_df.MAP.nunique())

        gene_df = gene_df.rename(convert_tissue_dict, axis=1)
        gene_df = gene_df.loc[gene_df["variable"] == "PEXT_value"]
        z_df = gene_df[gene_df.columns[8:]].apply(zscore, axis=1)
        z_df = pd.DataFrame.from_dict(dict(zip(z_df.index, z_df.values))).T
        z_df.columns = list(gene_df.columns[8:])
        # z_df = pd.concat([gene_df[gene_df.columns[:8]], z_df], axis=1)
        z_df2 = pd.concat(
            [
                gene_df[["MAP"]].reset_index(drop=True),
                z_df.add_suffix("_zscore").reset_index(drop=True),
                gene_df[gene_df.columns[8:]]
                .add_suffix("_pextvalue")
                .reset_index(drop=True),
            ],
            axis=1,
        )

        def apply_tissues(r, z_score_list):
            l = dict()
            r_filter_zscore = r.filter(regex="zscore")
            r_filter_pextvalue = r.filter(regex="pextvalue")
            values = dict()

            for i in ["Medium", "High", "Very-High"]:
                r[i] = []
            for (t1, z), (t2, p) in zip(
                r_filter_zscore.items(), r_filter_pextvalue.items()
            ):
                t = t1.replace("_zscore", "")
                if z_score_list[0] <= z < z_score_list[1] and p > 0.1:
                    r["Medium"].append(t)
                    values[t1] = z
                    values[t2] = p
                elif z_score_list[1] <= z < z_score_list[2] and p > 0.1:
                    r["High"].append(t)
                    values[t1] = z
                    values[t2] = p
                elif z >= z_score_list[2] and p > 0.1:
                    r["Very-High"].append(t)
                    values[t1] = z
                    values[t2] = p

            check = False
            Tissues = list()

            total_cat = list()

            for t, e in r[-3:].items():
                Tissues += e
                if e:
                    check = True
                    total_cat.append(t)
            if check is True:
                r["Tissue-specific"] = True
            else:
                r["Tissue-specific"] = False
            r["Tissues"] = Tissues
            r["Total_cat"] = total_cat
            r["Values"] = values

            return r

        z_score_list = [
            2.3263,
            3.0902,
            3.719,
        ]
        # z_df["H_profile"] = z_df.apply(lambda r: apply_specific(r[8:]), axis=1)
        # z_df = z_df.head(20).apply(lambda r: apply_tissues(r, z_score_list), axis=1)
        z_df = z_df2.apply(lambda r: apply_tissues(r, z_score_list), axis=1)
        # print(z_df)
        # exit()

        z_df = z_df[z_df.columns.drop(list(z_df.filter(regex="pext|zscore")))]
        z_df = z_df.loc[z_df["Tissue-specific"] == True]

        # for col in z_df.columns[1:4]:
        #     print(col, z_df.loc[z_df[col].str.len() != 0].shape[0], collections.Counter([e for r in z_df[col] for e in r]))
        # print("\n")

        # for col in z_df.columns[4:7]:
        #     print(z_df[col].value_counts())
        #     print("\n")

        # pprint(collections.Counter([e for r in z_df["Tissues"] for e in r]))
        # pprint(collections.Counter([e for r in z_df["Total_cat"] for e in r]))
        # pprint(collections.Counter(["_".join(list(set([e.split(" - ")[0] for e in r]))) for r in z_df["Tissues"] if len(r) > 1]))
        # pprint(collections.Counter(["_".join(list(set([e for e in r]))) for r in z_df["Tissues"] if len(r) > 1]))
        # print("\n")
        # print(np.mean([len(r) for r in z_df["Tissues"]]))
        # print(np.std([len(r) for r in z_df["Tissues"]]))
        # print(np.median([len(r) for r in z_df["Tissues"]]))

        # z_df_tmp = z_df.copy()
        # z_df_tmp["Total_cat"] = z_df_tmp["Total_cat"].astype(str)
        # print(z_df_tmp)
        # # exit()

        # for col in ["Low", "Medium", "High"]:
        #     print(col)
        #     pprint(collections.Counter(["_".join(list(set([e.split(" - ")[0] for e in r]))) for r in z_df[col]]))
        #     pprint(collections.Counter(["_".join(list(set([e for e in r]))) for r in z_df[col]]))
        #     pprint(collections.Counter(["_".join(list(set([e.split(" - ")[0] for e in r]))) for r in z_df[col] if len(r) > 1]))
        #     pprint(collections.Counter(["_".join(list(set([e for e in r]))) for r in z_df[col] if len(r) > 1]))
        #     c_z = collections.defaultdict(list)
        #     c_p = collections.defaultdict(list)
        #     print(z_df_tmp.loc[(z_df_tmp["Total_cat"] == "['{}']".format(col))])
        #     [
        #         c_z[k.replace("_zscore", "")].append(v)
        #         for r in z_df_tmp.loc[(z_df_tmp["Total_cat"] == "['{}']".format(col)), "Values"].values
        #         for k, v in r.items()
        #         if "zscore" in k
        #     ]
        #     [
        #         c_p[k.replace("_pextvalue", "")].append(v)
        #         for r in z_df_tmp.loc[(z_df_tmp["Total_cat"] == "['{}']".format(col)), "Values"].values
        #         for k, v in r.items()
        #         if "pext" in k
        #     ]
        #     tmp_l = list()
        #     for t1, t2 in zip(c_z, c_p):
        #         tmp_l.append(
        #             {
        #                 "Tissue": t1,
        #                 "Mean_z": np.mean(c_z[t1]),
        #                 "Std_z": np.std(c_z[t1]),
        #                 "Median_z": np.median(c_z[t1]),
        #                 "Mean_p": np.mean(c_p[t1]),
        #                 "Std_p": np.std(c_p[t1]),
        #                 "Median_p": np.median(c_p[t1]),
        #                 "Len": len(c_z[t1]),
        #             }
        #         )
        #     print_df = pd.DataFrame(tmp_l)
        #     print_df.loc["Mean"] = print_df.mean(axis=0)
        #     print_df.loc["Std"] = print_df.std(axis=0)
        #     print_df.loc["Median"] = print_df.median(axis=0)
        #     print(print_df.sort_values(by="Tissue"))
        #     print(print_df.sort_values(by="Mean_z", ascending=False))
        #     print(print_df.sort_values(by="Mean_p", ascending=False))
        #     print("\n")

        z_df[["symbol", "Exon"]] = z_df.MAP.str.split("_", expand=True)
        print("genes", str(z_df.symbol.nunique()))
        print("exons", str(z_df.MAP.nunique()))
        exit()
        # print(disease_genes_multi_iso)
        # print(z_df)
        # print(z_df.symbol.nunique())
        # print(z_df.MAP.nunique())

        # print(z_df.loc[z_df["symbol"].isin(disease_genes_multi_iso.Name.unique().tolist())])
        # print(z_df.loc[z_df["symbol"].isin(disease_genes_multi_iso.Name.unique().tolist()), "symbol"].nunique())
        # print(z_df.loc[z_df["symbol"].isin(disease_genes_multi_iso.Name.unique().tolist()), "MAP"].nunique())

        # print(z_df.loc[~z_df["symbol"].isin(disease_genes_multi_iso.Name.unique().tolist())])
        # print(z_df.loc[~z_df["symbol"].isin(disease_genes_multi_iso.Name.unique().tolist()), "symbol"].nunique())
        # print(z_df.loc[~z_df["symbol"].isin(disease_genes_multi_iso.Name.unique().tolist()), "MAP"].nunique())

        # exit()
        return z_df

        # z_df.to_excel(output_dir + "exons_zscore.xlsx", index=False)
        # gene_df.to_excel(output_dir + "exons_pext_value_absence.xlsx", index=False)

        # check_pext_value = [
        #     True
        #     if True in [True if tissue > config_search["Presence"]["Specific_min_expression"] else False for tissue in row]
        #     else False
        #     for row in gene_df.loc[gene_df["variable"] == "PEXT_value", gene_df.columns[8:]].values
        # ]

        # tissue_pext_value = [
        #     [c for c, value in row.items() if value > config_search["Presence"]["Specific_min_expression"]]
        #     for i, row in gene_df.loc[gene_df["variable"] == "PEXT_value", gene_df.columns[8:]].iterrows()
        # ]

        # check_pext_pct = [
        #     True
        #     if True in [True if tissue > config_search["Presence"]["Percentile_cutoff"] else False for tissue in row]
        #     else False
        #     for row in gene_df.loc[gene_df["variable"] == "PEXT_pct", gene_df.columns[8:]].values
        # ]

        # tissue_pext_pct = [
        #     [c for c, value in row.items() if value > config_search["Presence"]["Percentile_cutoff"]]
        #     for i, row in gene_df.loc[gene_df["variable"] == "PEXT_pct", gene_df.columns[8:]].iterrows()
        # ]

        # gene_df.loc[gene_df["variable"] == "PEXT_value", "Check"] = check_pext_value
        # gene_df.loc[gene_df["variable"] == "PEXT_value", "Tissues"] = tissue_pext_value

        # gene_df.loc[gene_df["variable"] == "PEXT_pct", "Check"] = check_pext_pct
        # gene_df.loc[gene_df["variable"] == "PEXT_pct", "Tissues"] = tissue_pext_pct

        # return gene_df

    def compare_exons_gtex_profile(self, df, exon_usage=False):
        def gtex_rank_v_profile(r):
            tissues_rank = self.GTEx_v_profile.loc[
                self.GTEx_v_profile["symbol"] == r["symbol"], r["OK"]
            ].values.tolist()
            return tissues_rank[0]

        def apply_outer_tissues_ok_pext_gtex(r):
            pext_tissues = [str(e) for e in r["OK"]]
            try:
                gtex_tissues = list(r["H_tissues"])
            except TypeError:
                gtex_tissues = []

            if gtex_tissues:
                try:
                    return list(set(pext_tissues).difference(set(gtex_tissues)))
                except TypeError:
                    print(pext_tissues, gtex_tissues)
            else:
                return pext_tissues

        merge_df = pd.merge(
            df, self.GTEx_v_profile[["symbol", "H_profile", "H_tissues"]], on="symbol"
        )

        merge_df["Tissues_rank"] = merge_df.apply(gtex_rank_v_profile, axis=1)

        # merge_df["PEXT_specific_tissues"] = merge_df.apply(apply_outer_tissues_ok_pext_gtex, axis=1)
        # merge_df["PEXT_specific_tissues"] = merge_df.apply(apply_outer_tissues_ok_pext_gtex, axis=1)

        # merge_df = merge_df.dropna(subset=["PEXT_specific_tissues"])

        # merge_df = merge_df[merge_df["PEXT_specific_tissues"].map(lambda d: len(d)) > 0]

        # print(merge_df.symbol.nunique())
        # print(merge_df.MAP.nunique())

        merge_df["Supertissues_gene"] = merge_df.loc[
            merge_df["H_profile"] == "Specific"
        ].H_tissues.apply(lambda r: list(set([e.split(" - ")[0] for e in list(r)])))
        merge_df["Supertissues_gene"] = merge_df["Supertissues_gene"].apply(
            lambda d: d if isinstance(d, list) else []
        )
        merge_df["Supertissues_exon"] = merge_df.OK.apply(
            lambda r: list(set([e.split(" - ")[0] for e in list(r)]))
        )

        merge_df["Corrected_Supertissues"] = merge_df.apply(
            lambda r: list(
                set(r["Supertissues_exon"]).difference(set(r["Supertissues_gene"]))
            ),
            axis=1,
        )
        merge_df["Corrected_Tissues"] = merge_df.loc[
            merge_df["H_profile"] == "Specific"
        ].apply(
            lambda r: [
                e for e in r["OK"] if e.split(" - ")[0] in r["Corrected_Supertissues"]
            ],
            axis=1,
        )
        merge_df.loc[merge_df["H_profile"] == "All", "Corrected_Tissues"] = merge_df[
            "OK"
        ]
        # merge_df = merge_df.loc[merge_df["Corrected_Tissues"].str.len() > 0]

        merge_df.to_parquet(self.exotic_files["exotic_gtex_comparison_path"])

        # print(merge_df[["symbol", "H_profile", "H_tissues", "Tissues", "Corrected_Tissues"]])
        # print(merge_df.symbol.nunique())
        # print(merge_df.MAP.nunique())

        # out = merge_df.loc[merge_df["Corrected_Tissues"].str.len() == 0]
        # print(out[["symbol", "H_profile", "H_tissues", "Tissues", "Corrected_Tissues"]])
        # print(out.symbol.nunique())
        # print(out.MAP.nunique())

        # tmp = merge_df.loc[merge_df["Corrected_Tissues"].str.len() > 0]
        # print(tmp[["symbol", "H_profile", "H_tissues", "Tissues", "Corrected_Tissues"]])
        # print(tmp.symbol.nunique())
        # print(tmp.MAP.nunique())
        # exit()

        # if exon_usage is True:
        #     print("Compare EXOTIC identified tissues & GTEx tissues")
        #     merge_df = merge_df.loc[merge_df["Corrected_Tissues"].str.len() > 0]
        #     merge_df.loc[merge_df["H_profile"] == "Specific"]

        # print(merge_df.loc[merge_df['H_profile'] == 'Specific'].symbol.nunique())
        # print(merge_df.loc[merge_df['H_profile'] == 'Specific'].MAP.nunique())

        # print(merge_df.loc[merge_df['H_profile'] == 'All'])
        # print(merge_df.loc[merge_df['H_profile'] == 'All'].symbol.nunique())
        # print(merge_df.loc[merge_df['H_profile'] == 'All'].MAP.nunique())
        # exit()

        # print(merge_df)
        # print(merge_df.symbol.nunique())
        # print(merge_df.MAP.nunique())
        return merge_df

    def compare_to_omim_basic_gene_level(
        self,
        complete_df,
        omim,
        d_genes,
        omim_detailed=True,
        exon_usage=False,
    ):

        ## GTEx nbs
        # df = complete_df.loc[complete_df["H_profile"].isin(condition), ["symbol", "H_tissues"]]

        complete_df = complete_df.loc[complete_df["symbol"].isin(d_genes.Name.unique())]

        df = complete_df.loc[
            complete_df["H_profile"] == "Specific", ["symbol", "H_tissues", "H_profile"]
        ]
        df["Supertissues"] = df["H_tissues"].apply(
            lambda r: list(set([e.split(" - ")[0] for e in r]))
        )

        # print(df[df['H_tissues'].map(lambda d: len(d)) > 0])

        ## OMIM nbs
        # print(omim.OMIM.nunique())
        # print(omim.loc[~omim['OMIM'].isin(
        #     d_genes.OMIM.values.tolist()), 'OMIM'].nunique())
        # print(omim.loc[omim['OMIM'].isin(
        #     d_genes.OMIM.values.tolist()), 'OMIM's].nunique())
        # omim = omim.loc[omim["OMIM"].isin(d_genes.OMIM.values.tolist())]

        mapping_omim_gtex = self.dicts["mapping_omim_gtex_detailed"]

        omim = omim.where(pd.notnull(omim), None)

        # TEST DUPLICATED PHENO OMIM
        # duplicated = omim.loc[omim.duplicated(keep=False, subset=['OMIM', 'Pheno_OMIM', 'Pheno_prefered_title'])]
        # duplicated.apply(omim_compare_pheno_name_to_prefered_title, axis=1)

        omim = omim.loc[
            ~omim.duplicated(
                keep="last", subset=["OMIM", "Pheno_OMIM", "Pheno_prefered_title"]
            )
        ]

        # df = df.loc[df['variable'] == 'PEXT_value']
        # df = df.loc[df['symbol'] == 'BIN1']
        omim_associations = list()

        for i, row in df.iterrows():
            gene_omim_tmp = omim.loc[omim["Name"] == row["symbol"]]
            # print(row)
            # print(gene_omim_tmp)
            # exit()
            if gene_omim_tmp.empty is False:
                for j, pheno in gene_omim_tmp.iterrows():
                    omim_tmp = pheno.to_dict()

                    # print(omim_tmp)
                    if type(row["H_tissues"]) is np.ndarray:
                        row["H_tissues"] = row["H_tissues"].tolist()
                    elif type(row["H_tissues"]) is str:
                        row["H_tissues"] = eval(row["H_tissues"])

                    if row["H_tissues"] is not None:

                        for t in row["H_tissues"]:

                            if t in mapping_omim_gtex:

                                for omim_body_part in mapping_omim_gtex[t]:
                                    # if omim_tmp[omim_body_part]:
                                    # print(omim_body_part in omim_tmp)
                                    if (
                                        omim_body_part in omim_tmp
                                        and omim_tmp[omim_body_part] is not None
                                    ):
                                        # print(omim_tmp[omim_body_part])
                                        omim_associations.append(
                                            {
                                                "symbol": row["symbol"],
                                                # 'MAP': row['MAP'],
                                                "OMIM_pheno_gene_level": pheno[
                                                    "Pheno_OMIM"
                                                ],
                                                "OMIM_pheno_name_gene_level": pheno[
                                                    "Pheno_Name"
                                                ],
                                                "Tissue_gene_level": t,
                                                "Supertissue": t.split(" - ")[0],
                                                "OMIM_body_part_gene_level": omim_body_part,
                                                # "OMIM_details_gene_level": [
                                                # e.split(" {")[0][1:] for e in omim_tmp[omim_body_part].split("\n"),
                                                "OMIM_details_gene_level": [
                                                    e.split(" {")[0]
                                                    for e in eval(
                                                        omim_tmp[omim_body_part]
                                                    )[0].splitlines()
                                                ],
                                            }
                                        )

        omim_associations = pd.DataFrame(omim_associations)

        # omim_associations = omim_associations.loc[omim_associations['Supertissue'] != 'Brain']

        # print(df)
        # print(df.symbol.nunique())

        # omim_associations[["symbol", "Tissue_gene_level", "OMIM_body_part_gene_level"]].groupby("symbol").nunique().to_excel(
        #     self.output_dir + self.subdirs_output_dir[subdir] + "table_gene_level_tissues_omim_gene_contribution.xlsx"
        # )

        # omim_associations[["symbol", "Tissue_gene_level"]].drop_duplicates(subset=["symbol", "Tissue_gene_level"]).groupby(
        #     "Tissue_gene_level"
        # ).size().to_excel(self.output_dir + self.subdirs_output_dir[subdir] + "table_gene_level_tissues_detailed.xlsx")

        # omim_associations[["symbol", "Tissue_gene_level", "Supertissue"]].drop_duplicates(
        #     subset=["symbol", "Tissue_gene_level"]
        # ).groupby("Supertissue").size().to_excel(
        #     self.output_dir + self.subdirs_output_dir[subdir] + "table_gene_level_tissues_global.xlsx"
        # )

        # omim_associations[["symbol", "Tissue_gene_level", "OMIM_body_part_gene_level"]].groupby(
        #     "OMIM_body_part_gene_level"
        # ).size().to_excel(self.output_dir + self.subdirs_output_dir[subdir] + "table_gene_level_OMIM.xlsx")

        if exon_usage is False:
            print(self.output_dir + "PHENOTYPES/")
            # omim_associations.to_parquet(self.output_dir + 'PHENOTYPES/' + "matrix_omim_gene_level_wt_brain.xlsx", index=False)
            omim_associations.to_excel(
                self.output_dir
                + "PHENOTYPES/"
                + "matrix_omim_gene_level_corrected_neurologic.xlsx",
                index=False,
            )
            exit()

        elif exon_usage is True:

            omim_associations = (
                pd.merge(complete_df, omim_associations, on="symbol")
                .drop_duplicates(
                    subset=[
                        "MAP",
                        # "symbol",
                        "OMIM_pheno_name_gene_level",
                        "Tissue_gene_level",
                        "OMIM_body_part_gene_level",
                    ]
                )
                .reset_index(drop=True)
            )

        # omim_associations.to_excel(
        # self.output_dir + "PHENOTYPES/" + "matrix_omim_exon_level.xlsx"
        # )
        # tp_exon_genes = df.symbol.unique().tolist()
        # omim_associated_genes = omim_associations.symbol.unique().tolist()
        # df = complete_df.loc[(~complete_df["symbol"].isin(omim_associated_genes)) & (complete_df["H_profile"] == "All")]
        return omim_associations

        # print(omim_associations)
        # print(omim_associations.symbol.nunique())
        # print(omim_associations.MAP.nunique())
        # print(omim_associations.loc[omim_associations['H_profile'] == 'Specific'].drop_duplicates())
        # return omim_associations

        # print(row)
        # print(row['symbol'],
        #       row['PEXT_specific_tissues'])
        # print(omim_body_part)
        # print(omim_tmp[omim_body_part])
        # print('\n')

    def extract_tissues_ok(self, df):

        print(df)

        intermediate_file = df.drop_duplicates()
        intermediate_file = intermediate_file.loc[
            (intermediate_file["mean_proportion"] > 0)
            & (intermediate_file["mean_proportion"] < 1)
            # & (intermediate_file["mean_proportion"] < 0.1)
            # & (intermediate_file["Ratio_num"] <= 0.5)
        ].reset_index(drop=True)

        print(intermediate_file)

        pext = intermediate_file.filter(regex="pext")

        exotic = intermediate_file.filter(regex="exotic")
        pext["Max"] = pext.progress_apply(
            lambda r: max(r.dropna()) if r.dropna().shape[0] > 0 else np.nan, axis=1
        )
        pext["MAP"] = intermediate_file["MAP"]
        pext["symbol"] = intermediate_file["symbol"]

        exotic["Max"] = exotic.progress_apply(
            lambda r: max(r.dropna()) if r.dropna().shape[0] > 0 else np.nan, axis=1
        )
        exotic["Tissues_Max"] = exotic.progress_apply(
            lambda r: [
                list(exotic.columns)[j].replace("_exotic", "")
                for j, e in enumerate(r[:-1])
                if e == r["Max"]
            ],
            axis=1,
        )
        exotic["Tissues_Max_organs"] = exotic["Tissues_Max"].apply(
            lambda r: list(set([e.split(" - ")[0] for e in r]))
        )
        exotic["pext_max"] = pext["Max"]

        exotic["MAP"] = intermediate_file["MAP"]
        exotic["Ratio_num"] = intermediate_file["Ratio_num"]
        exotic["mRNA_nb"] = intermediate_file["mRNA_nb"]
        exotic["symbol"] = intermediate_file["symbol"]
        exotic["mean_proportion"] = intermediate_file["mean_proportion"]

        exotic["tissues_OK_bronze"] = (
            exotic.filter(regex="exotic")
            .fillna(0)
            .apply(
                lambda r: [
                    exotic.columns[j].replace("_exotic", "")
                    for j, c in enumerate(r)
                    if c >= 0.90 and c < 0.95
                ],
                axis=1,
            )
        )
        exotic["tissues_OK_silver"] = (
            exotic.filter(regex="exotic")
            .fillna(0)
            .apply(
                lambda r: [
                    exotic.columns[j].replace("_exotic", "")
                    for j, c in enumerate(r)
                    if c >= 0.95 and c < 0.99
                ],
                axis=1,
            )
        )
        exotic["tissues_OK_gold"] = (
            exotic.filter(regex="exotic")
            .fillna(0)
            .apply(
                lambda r: [
                    exotic.columns[j].replace("_exotic", "")
                    for j, c in enumerate(r)
                    if c >= 0.99
                ],
                axis=1,
            )
        )
        exotic["pext_OK"] = (
            pext.filter(regex="pext")
            .fillna(0)
            .apply(
                lambda r: [
                    pext.columns[j].replace("_pextvalue", "")
                    for j, c in enumerate(r)
                    if c > 0.1
                ],
                axis=1,
            )
        )
        # exotic["No_threshold"] = exotic.apply(
        #     lambda r: list(set(r["tissues_no_threshold"]).intersection(set(r["pext_OK"]))), axis=1
        # )
        exotic["OK_bronze"] = exotic.apply(
            lambda r: list(set(r["tissues_OK_bronze"]).intersection(set(r["pext_OK"]))),
            axis=1,
        )
        exotic["OK_silver"] = exotic.apply(
            lambda r: list(set(r["tissues_OK_silver"]).intersection(set(r["pext_OK"]))),
            axis=1,
        )
        exotic["OK_gold"] = exotic.apply(
            lambda r: list(set(r["tissues_OK_gold"]).intersection(set(r["pext_OK"]))),
            axis=1,
        )
        exotic["OK"] = exotic.apply(
            lambda r: r["OK_bronze"] + r["OK_silver"] + r["OK_gold"], axis=1
        )

        exotic = exotic[
            [
                "MAP",
                "Ratio_num",
                "mRNA_nb",
                "symbol",
                "Adipose - Subcutaneous_exotic",
                "Adipose - Visceral (Omentum)_exotic",
                "Adrenal Gland_exotic",
                "Artery - Aorta_exotic",
                "Artery - Coronary_exotic",
                "Artery - Tibial_exotic",
                "Bladder_exotic",
                "Brain - Amygdala_exotic",
                "Brain - Anterior cingulate cortex (BA24)_exotic",
                "Brain - Caudate (basal ganglia)_exotic",
                "Brain - Cerebellar Hemisphere_exotic",
                "Brain - Cerebellum_exotic",
                "Brain - Cortex_exotic",
                "Brain - Frontal Cortex (BA9)_exotic",
                "Brain - Hippocampus_exotic",
                "Brain - Hypothalamus_exotic",
                "Brain - Nucleus accumbens (basal ganglia)_exotic",
                "Brain - Putamen (basal ganglia)_exotic",
                "Brain - Spinal cord (cervical c-1)_exotic",
                "Brain - Substantia nigra_exotic",
                "Breast - Mammary Tissue_exotic",
                "Cells - Cultured fibroblasts_exotic",
                "Cells - EBV-transformed lymphocytes_exotic",
                "Cervix - Ectocervix_exotic",
                "Cervix - Endocervix_exotic",
                "Colon - Sigmoid_exotic",
                "Colon - Transverse_exotic",
                "Esophagus - Gastroesophageal Junction_exotic",
                "Esophagus - Mucosa_exotic",
                "Esophagus - Muscularis_exotic",
                "Fallopian Tube_exotic",
                "Heart - Atrial Appendage_exotic",
                "Heart - Left Ventricle_exotic",
                "Kidney - Cortex_exotic",
                "Liver_exotic",
                "Lung_exotic",
                "Minor Salivary Gland_exotic",
                "Muscle - Skeletal_exotic",
                "Nerve - Tibial_exotic",
                "Ovary_exotic",
                "Pancreas_exotic",
                "Pituitary_exotic",
                "Prostate_exotic",
                "Skin - Not Sun Exposed (Suprapubic)_exotic",
                "Skin - Sun Exposed (Lower leg)_exotic",
                "Small Intestine - Terminal Ileum_exotic",
                "Spleen_exotic",
                "Stomach_exotic",
                "Testis_exotic",
                "Thyroid_exotic",
                "Uterus_exotic",
                "Vagina_exotic",
                "Whole Blood_exotic",
                "mean_proportion",
                "Max",
                "Tissues_Max",
                "Tissues_Max_organs",
                "pext_max",
                "tissues_OK_bronze",
                "tissues_OK_silver",
                "tissues_OK_gold",
                "pext_OK",
                "OK_bronze",
                "OK_silver",
                "OK_gold",
                "OK",
            ]
        ]

        exotic.to_parquet(self.exotic_files["exotic_ok_path"])
        return exotic

    def compare_to_omim_basic_gene_exon_level(
        self, df, omim, d_genes, omim_detailed=True
    ):

        df = df.loc[df["Corrected_Tissues"].str.len() > 0]
        df = df.loc[df["OK"].str.len() > 0]

        gene_level_df = self.compare_to_omim_basic_gene_level(
            df,
            omim,
            d_genes,
            omim_detailed=True,
            exon_usage=True,
        )

        print(gene_level_df.symbol.nunique())
        print(gene_level_df.MAP.nunique())
        exit()

        mapping_omim_gtex = self.dicts["mapping_omim_gtex_detailed"]

        omim = omim.where(pd.notnull(omim), None)

        omim = omim.loc[
            ~omim.duplicated(
                keep="last", subset=["OMIM", "Pheno_OMIM", "Pheno_prefered_title"]
            )
        ]

        # print(gene_level_df.symbol.nunique())
        # print(gene_level_df.MAP.nunique())

        omim_associations = list()
        for i, row in gene_level_df.iterrows():
            gene_omim_tmp = omim.loc[omim["Name"] == row["symbol"]]
            if gene_omim_tmp.empty is False:
                for j, pheno in gene_omim_tmp.iterrows():
                    omim_tmp = pheno.to_dict()
                    if type(row["Corrected_Tissues"]) is str:
                        row["Corrected_Tissues"] = eval(row["Corrected_Tissues"])
                    # print(row["Corrected_Tissues"], type(row["Corrected_Tissues"]))
                    # exit()

                    for t in row["Corrected_Tissues"]:
                        if t in mapping_omim_gtex:
                            for omim_body_part in mapping_omim_gtex[t]:
                                if (
                                    omim_body_part in omim_tmp
                                    and omim_tmp[omim_body_part] is not None
                                ):
                                    threshold = [
                                        e.replace("OK_", "")
                                        for e in ["OK_bronze", "OK_silver", "OK_gold"]
                                        if t in row[e]
                                    ][0]
                                    omim_associations.append(
                                        {
                                            "symbol": row["symbol"],
                                            "MAP": row["MAP"],
                                            "Corrected_Tissues": row["OK"],
                                            "H_tissues": row["H_tissues"],
                                            "Threshold": threshold,
                                            "EXOTIC_value": row[t + "_exotic"],
                                            "OMIM_pheno_exon_level": pheno[
                                                "Pheno_OMIM"
                                            ],
                                            "OMIM_pheno_name_exon_level": pheno[
                                                "Pheno_Name"
                                            ],
                                            "Tissue_exon_level": t,
                                            "OMIM_body_part_exon_level": omim_body_part,
                                            "OMIM_details_exon_level": [
                                                e.split(" {")[0]
                                                for e in eval(omim_tmp[omim_body_part])[
                                                    0
                                                ].splitlines()
                                            ],
                                        }
                                    )
        omim_associations = pd.DataFrame(omim_associations)
        # print(omim_associations)
        # print(omim_associations.MAP.nunique())
        # print(omim_associations.symbol.nunique())

        def compare_omim_gene_exon_level(r):
            r_gene = r["OMIM_details_gene_level"]
            r_exon = r["OMIM_details_exon_level"]
            r["OMIM_details_common"] = list(set(r_gene).intersection(set(r_exon)))
            r["OMIM_details_spec_gene"] = list(set(r_gene).difference(set(r_exon)))
            r["OMIM_details_spec_exon"] = list(set(r_exon).difference(set(r_gene)))
            return r

        # output_df = (
        #     pd.merge(df, omim_associations, on=["symbol", "MAP"])
        #     .drop_duplicates(subset=["MAP", "OMIM_pheno_exon_level", "OMIM_body_part_exon_level", "OMIM_body_part_exon_level"])
        #     .reset_index(drop=True)
        #     .apply(compare_omim_gene_exon_level, axis=1)
        # )
        # output_df = output_df.loc[
        #     (output_df["OMIM_details_spec_gene"].str.len() > 0) & (output_df["OMIM_details_spec_exon"].str.len() > 0)
        # ]

        output_df = (
            pd.merge(gene_level_df, omim_associations, on=["symbol", "MAP"])
            .drop_duplicates(
                subset=[
                    "MAP",
                    "OMIM_pheno_exon_level",
                    "OMIM_body_part_exon_level",
                    "OMIM_body_part_exon_level",
                ]
            )
            .reset_index(drop=True)
        )

        output_df = output_df[
            [c for c in list(output_df.columns) if "_exotic" not in c]
        ]

        print(output_df.loc[output_df["symbol"] == "TBC1D24"])

        # present_in_omim = gene_level_df.loc[gene_level_df["symbol"].isin(omim.Name.values.tolist())]
        # print(present_in_omim.symbol.nunique())
        # print(present_in_omim.MAP.nunique())
        # print(present_in_omim[["symbol", "H_profile", "H_tissues", "Tissues"]].drop_duplicates(subset=["symbol", "Tissue_exon_level", 'OMIM_pheno_exon_level']))
        # print("\n")

        print(output_df)
        print(output_df.symbol.nunique())
        print(output_df.MAP.nunique())
        print(output_df.loc[output_df["Threshold"] == "bronze", "symbol"].nunique())
        print(output_df.loc[output_df["Threshold"] == "bronze", "MAP"].nunique())
        print(output_df.loc[output_df["Threshold"] == "silver", "symbol"].nunique())
        print(output_df.loc[output_df["Threshold"] == "silver", "MAP"].nunique())
        print(output_df.loc[output_df["Threshold"] == "gold", "symbol"].nunique())
        print(output_df.loc[output_df["Threshold"] == "gold", "MAP"].nunique())

        output_df.to_excel(
            self.output_dir + "PHENOTYPES/" + "table_OMIM_gene_exon_level_matrix.xlsx",
            index=False,
        )

        # print(
        #     output_df[
        #         [
        #             "symbol",
        #             "MAP",
        #             "H_profile",
        #             "H_tissues_x",
        #             "OK",
        #             "Tissue_exon_level",
        #             "OMIM_pheno_exon_level",
        #         ]
        #     ]
        #     .drop_duplicates(
        #         subset=["MAP", "Tissue_exon_level", "OMIM_pheno_exon_level"]
        #     )
        #     .shape
        # )
        # print("\n")

        # undupl = output_df[
        #     [
        #         "symbol",
        #         "MAP",
        #         "H_profile",
        #         "H_tissues_x",
        #         "OK",
        #         "Tissue_exon_level",
        #         "OMIM_pheno_exon_level",
        #         "OMIM_body_part_exon_level",
        #     ]
        # ].drop_duplicates(subset=["MAP", "Tissue_exon_level", "OMIM_pheno_exon_level"])
        # print(undupl)

        # undupl[["Tissue_exon_level", "OMIM_body_part_exon_level"]].groupby(
        #     "Tissue_exon_level"
        # ).size().to_excel(
        #     self.output_dir + "PHENOTYPES/" + "table_OMIM_exon_level_tissues.xlsx"
        # )
        # undupl[["Tissue_exon_level", "OMIM_body_part_exon_level"]].groupby(
        #     "OMIM_body_part_exon_level"
        # ).size().to_excel(
        #     self.output_dir + "PHENOTYPES/" + "table_OMIM_exon_level_BP.xlsx"
        # )

        # undupl = output_df[
        #     [
        #         "symbol",
        #         "MAP",
        #         "H_profile",
        #         "H_tissues_x",
        #         "OK",
        #         "Tissue_gene_level",
        #         "OMIM_pheno_gene_level",
        #         "OMIM_body_part_gene_level",
        #     ]
        # ].drop_duplicates(subset=["MAP", "Tissue_gene_level", "OMIM_pheno_gene_level"])

        # undupl[["Tissue_gene_level", "OMIM_body_part_gene_level"]].groupby(
        #     "Tissue_gene_level"
        # ).size().to_excel(
        #     self.output_dir + "PHENOTYPES/" + "table_OMIM_gene_exon_level_tissues.xlsx"
        # )
        # undupl[["Tissue_gene_level", "OMIM_body_part_gene_level"]].groupby(
        #     "OMIM_body_part_gene_level"
        # ).size().to_excel(
        #     self.output_dir + "PHENOTYPES/" + "table_OMIM_gene_exon_level_BP.xlsx"
        # )

        # output_df.to_excel(output_dir + "genes_exon_omim_detailed3.xlsx", index=False)
        # exit()

        # print(omim_associations.Gene.nunique())
        # print(omim_associations.MAP.nunique())
        # print(omim_associations)
        # exit()
        # print(omim_associations[["OMIM_pheno", "Gene"]].drop_duplicates())
        # print(omim_associations.MAP.nunique())
        # print(
        # pd.DataFrame(
        # collections.Counter(omim_associations.OMIM_body_part.values.tolist()).items(), columns=["Tissue", "Value"]
        # )
        # )

        return output_df

    def compare_to_omim_basic_exon_level(self, df, omim, d_genes, omim_detailed=True):

        # print(df)
        # print(df.symbol.nunique())
        # print(df.MAP.nunique())

        df = df.loc[df["OK"].str.len() > 0]

        # df = df.loc[df["mean_proportion"] > 0.1]

        # print(df)
        # print(df.symbol.nunique())
        # print(df.MAP.nunique())
        # print(df.loc[df["OK_bronze"].str.len() > 0, "symbol"].nunique())
        # print(df.loc[df["OK_silver"].str.len() > 0, "symbol"].nunique())
        # print(df.loc[df["OK_gold"].str.len() > 0, "symbol"].nunique())
        # print(df.loc[df["OK_bronze"].str.len() > 0, "MAP"].nunique())
        # print(df.loc[df["OK_silver"].str.len() > 0, "MAP"].nunique())
        # print(df.loc[df["OK_gold"].str.len() > 0, "MAP"].nunique())
        # exit()

        mapping_omim_gtex = self.dicts["mapping_omim_gtex_detailed"]

        omim = omim.where(pd.notnull(omim), None)

        omim = omim.loc[
            ~omim.duplicated(
                keep="last", subset=["OMIM", "Pheno_OMIM", "Pheno_prefered_title"]
            )
        ]

        omim_associations = list()
        for i, row in df.iterrows():
            gene_omim_tmp = omim.loc[omim["Name"] == row["symbol"]]
            if gene_omim_tmp.empty is False:
                for j, pheno in gene_omim_tmp.iterrows():
                    omim_tmp = pheno.to_dict()
                    if type(row["OK"]) is str:
                        row["OK"] = eval(row["OK"])
                    # print(row["Corrected_Tissues"], type(row["Corrected_Tissues"]))
                    # exit()

                    for t in row["OK"]:
                        if t in mapping_omim_gtex:
                            for omim_body_part in mapping_omim_gtex[t]:
                                if (
                                    omim_body_part in omim_tmp
                                    and omim_tmp[omim_body_part] is not None
                                ):
                                    threshold = [
                                        e.replace("OK_", "")
                                        for e in ["OK_bronze", "OK_silver", "OK_gold"]
                                        if t in row[e]
                                    ][0]

                                    omim_associations.append(
                                        {
                                            "symbol": row["symbol"],
                                            "MAP": row["MAP"],
                                            "Corrected_Tissues": row["OK"],
                                            # "H_tissues": row["H_tissues"],
                                            "OMIM_pheno_exon_level": pheno[
                                                "Pheno_OMIM"
                                            ],
                                            "Threshold": threshold,
                                            "EXOTIC_value": row[t + "_exotic"],
                                            "OMIM_pheno_name_exon_level": pheno[
                                                "Pheno_Name"
                                            ],
                                            "Tissue_exon_level": t,
                                            "OMIM_body_part_exon_level": omim_body_part,
                                            "OMIM_details_exon_level": [
                                                e.split(" {")[0]
                                                for e in eval(omim_tmp[omim_body_part])[
                                                    0
                                                ].splitlines()
                                            ],
                                        }
                                    )
        omim_associations = pd.DataFrame(omim_associations)

        print(omim_associations)
        print(omim_associations.symbol.nunique())
        print(omim_associations.MAP.nunique())
        print(
            omim_associations.loc[
                omim_associations["Threshold"] == "bronze", "symbol"
            ].nunique()
        )
        print(
            omim_associations.loc[
                omim_associations["Threshold"] == "silver", "symbol"
            ].nunique()
        )
        print(
            omim_associations.loc[
                omim_associations["Threshold"] == "gold", "symbol"
            ].nunique()
        )
        print(
            omim_associations.loc[
                omim_associations["Threshold"] == "bronze", "MAP"
            ].nunique()
        )
        print(
            omim_associations.loc[
                omim_associations["Threshold"] == "silver", "MAP"
            ].nunique()
        )
        print(
            omim_associations.loc[
                omim_associations["Threshold"] == "gold", "MAP"
            ].nunique()
        )
        omim_associations.to_excel(
            self.output_dir
            + "PHENOTYPES/"
            + "table_OMIM_exon_level_matrix_corrected.xlsx",
            index=False,
        )

        def compare_omim_gene_exon_level(r):
            r_gene = r["OMIM_details_gene_level"]
            r_exon = r["OMIM_details_exon_level"]
            r["OMIM_details_common"] = list(set(r_gene).intersection(set(r_exon)))
            r["OMIM_details_spec_gene"] = list(set(r_gene).difference(set(r_exon)))
            r["OMIM_details_spec_exon"] = list(set(r_exon).difference(set(r_gene)))
            return r

        # output_df = (
        #     pd.merge(df, omim_associations, on=["symbol", "MAP"])
        #     .drop_duplicates(subset=["MAP", "OMIM_pheno_exon_level", "OMIM_body_part_exon_level", "OMIM_body_part_exon_level"])
        #     .reset_index(drop=True)
        #     .apply(compare_omim_gene_exon_level, axis=1)
        # )
        # output_df = output_df.loc[
        #     (output_df["OMIM_details_spec_gene"].str.len() > 0) & (output_df["OMIM_details_spec_exon"].str.len() > 0)
        # ]

        # output_df = (
        #     pd.merge(gene_level_df, omim_associations, on=["symbol", "MAP"])
        #     .drop_duplicates(
        #         subset=[
        #             "MAP",
        #             "OMIM_pheno_exon_level",
        #             "OMIM_body_part_exon_level",
        #             "OMIM_body_part_exon_level",
        #         ]
        #     )
        #     .reset_index(drop=True)
        # )

        print(omim_associations)

        # present_in_omim = gene_level_df.loc[gene_level_df["symbol"].isin(omim.Name.values.tolist())]
        # print(present_in_omim.symbol.nunique())
        # print(present_in_omim.MAP.nunique())
        # print(present_in_omim[["symbol", "H_profile", "H_tissues", "Tissues"]].drop_duplicates(subset=["symbol", "Tissue_exon_level", 'OMIM_pheno_exon_level']))
        # print("\n")

        print(omim_associations)
        print(omim_associations.symbol.nunique())
        print(omim_associations.MAP.nunique())
        # print(
        #     output_df[
        #         [
        #             "symbol",
        #             "MAP",
        #             "H_profile",
        #             "H_tissues_x",
        #             "OK",
        #             "Tissue_exon_level",
        #             "OMIM_pheno_exon_level",
        #         ]
        #     ]
        #     .drop_duplicates(
        #         subset=["MAP", "Tissue_exon_level", "OMIM_pheno_exon_level"]
        #     )
        #     .shape
        # )
        # print("\n")

        # undupl = output_df[
        #     [
        #         "symbol",
        #         "MAP",
        #         "H_profile",
        #         "H_tissues_x",
        #         "OK",
        #         "Tissue_exon_level",
        #         "OMIM_pheno_exon_level",
        #         "OMIM_body_part_exon_level",
        #     ]
        # ].drop_duplicates(subset=["MAP", "Tissue_exon_level", "OMIM_pheno_exon_level"])
        # print(undupl)
        # output_df.to_excel(
        #     self.output_dir + "PHENOTYPES/" + "table_OMIM_gene_exon_level_matrix.xlsx",
        #     index=False,
        # )
        # undupl[["Tissue_exon_level", "OMIM_body_part_exon_level"]].groupby(
        #     "Tissue_exon_level"
        # ).size().to_excel(
        #     self.output_dir + "PHENOTYPES/" + "table_OMIM_exon_level_tissues.xlsx"
        # )
        # undupl[["Tissue_exon_level", "OMIM_body_part_exon_level"]].groupby(
        #     "OMIM_body_part_exon_level"
        # ).size().to_excel(
        #     self.output_dir + "PHENOTYPES/" + "table_OMIM_exon_level_BP.xlsx"
        # )

        # undupl = output_df[
        #     [
        #         "symbol",
        #         "MAP",
        #         "H_profile",
        #         "H_tissues_x",
        #         "OK",
        #         "Tissue_gene_level",
        #         "OMIM_pheno_gene_level",
        #         "OMIM_body_part_gene_level",
        #     ]
        # ].drop_duplicates(subset=["MAP", "Tissue_gene_level", "OMIM_pheno_gene_level"])

        # undupl[["Tissue_gene_level", "OMIM_body_part_gene_level"]].groupby(
        #     "Tissue_gene_level"
        # ).size().to_excel(
        #     self.output_dir + "PHENOTYPES/" + "table_OMIM_gene_exon_level_tissues.xlsx"
        # )
        # undupl[["Tissue_gene_level", "OMIM_body_part_gene_level"]].groupby(
        #     "OMIM_body_part_gene_level"
        # ).size().to_excel(
        #     self.output_dir + "PHENOTYPES/" + "table_OMIM_gene_exon_level_BP.xlsx"
        # )

        # output_df.to_excel(output_dir + "genes_exon_omim_detailed3.xlsx", index=False)
        # exit()

        # print(omim_associations.Gene.nunique())
        # print(omim_associations.MAP.nunique())
        # print(omim_associations)
        # exit()
        # print(omim_associations[["OMIM_pheno", "Gene"]].drop_duplicates())
        # print(omim_associations.MAP.nunique())
        # print(
        # pd.DataFrame(
        # collections.Counter(omim_associations.OMIM_body_part.values.tolist()).items(), columns=["Tissue", "Value"]
        # )
        # )
        # print(pd.DataFrame(collections.Counter(omim_associations.Tissue.values.tolist()).items(), columns=["Tissue", "Value"]))

        return omim_associations

    def mapping_clinvar(self, clinvar, df, pubmed):
        # clinvar[["CHROM", "POS", "REF", "ALT"]] = clinvar["VAR_ID"].str.split(
        #     "_", expand=True
        # )
        # clinvar["ALT"] = clinvar["ALT"].astype(str)
        # clinvar["ALT"] = clinvar["ALT"].str.replace("\['", "")
        # clinvar["ALT"] = clinvar["ALT"].str.replace("'\]", "")
        # clinvar = clinvar.loc[
        #     (clinvar["REF"].str.len() == 1) & (clinvar["ALT"].str.len() == 1)
        # ]

        df = df.loc[df["OK_gold"].str.len() > 0]

        print(clinvar)
        print(clinvar.loc[(clinvar["RS_STARS"] > 0)])
        exit()

        print(df)
        print(df.columns)

        print(pubmed)

        df = df.rename({"Gene": "symbol"}, axis=1)

        df["Exon"] = df["MAP"].apply(lambda r: r.split("_")[1])
        pubmed = pubmed.rename({"#AlleleID": "alleleid"}, axis=1)
        c = dict()
        for j, row in tqdm(df.iterrows()):
            cl = clinvar.loc[
                (clinvar["GENE"] == row["symbol"])
                & (clinvar["POS"] >= int(row["Exon"].split("-")[0]))
                & (clinvar["POS"] <= int(row["Exon"].split("-")[1]))
            ]
            if cl.empty is False:
                for i, clinvar_row in cl.iterrows():
                    c[clinvar_row["rs"]] = row["MAP"]
        clinvar["MAP"] = clinvar["rs"].map(c)
        df = pd.merge(
            df,
            clinvar[
                [
                    "rs",
                    "Status",
                    "Real_Status",
                    "RS_STARS",
                    "CLNREVSTAT",
                    "HPO",
                    "MC",
                    "MAP",
                    "alleleid",
                ]
            ],
            on="MAP",
        )
        lite_df = df

        # lite_df["PEXT_specific_tissues"] = lite_df["PEXT_specific_tissues"].astype(str)
        # lite_df["H_tissues"] = lite_df["H_tissues"].astype(str)
        # lite_df["alleleid"] = lite_df["alleleid"].astype(str)

        lite_df = lite_df.drop_duplicates(subset=["MAP", "rs"])

        # pd.options.display.max_colwidth = 50
        # pd.options.display.max_rows = 100

        lite_df = lite_df.sort_values(by="MAP")

        # print(lite_df.symbol.nunique())
        # print(lite_df.MAP.nunique())
        # print(lite_df.rs.nunique())
        # print(lite_df[lite_df["HPO"].map(lambda d: len(d)) > 0])
        # # print(lite_df)
        # print(lite_df[["OMIM_pheno"]].drop_duplicates())

        # # lite_df.to_excel("/home/weber/PycharmProjects/ExoCarto/data/clean/3_phenotypes/protocole_variants.xlsx")
        # print(pd.DataFrame(collections.Counter(lite_df.OMIM_body_part.values.tolist()).items(), columns=["Tissue", "Value"]))
        # print(pd.DataFrame(collections.Counter(lite_df.Tissue.values.tolist()).items(), columns=["Tissue", "Value"]))
        # print(pd.DataFrame(collections.Counter(lite_df.RS_STARS.values.tolist()).items(), columns=["RS", "Value"]))

        # exit()

        # lite_df = lite_df.loc[~lite_df["symbol"].isin(["ATM", "DMD"])]

        print(lite_df)
        print(lite_df.symbol.nunique())
        print(lite_df.MAP.nunique())
        # lite_df.to_excel(
        # self.output_dir + "VARIATIONS/" + "genes_clinvar_exotic_complete.xlsx",
        # index=False,
        # )
        # exit()

        lite_df = lite_df.dropna(subset=["alleleid"])
        lite_df = pd.merge(lite_df, pubmed[["alleleid", "citation_id"]], on="alleleid")
        lite_df_nbk = lite_df.loc[lite_df["citation_id"].str.contains("NBK")]
        lite_df = lite_df.loc[~lite_df["citation_id"].str.contains("NBK")]

        print(lite_df)
        print(lite_df.symbol.nunique())
        print(lite_df.MAP.nunique())

        lite_df.to_excel(
            self.output_dir + "VARIATIONS/" + "genes_pubmed_test.xlsx", index=False
        )
        exit()

        for j, r in lite_df.iterrows():
            pprint(r.to_dict())
            pm_pmc = self.convert_pmid_to_pmcid(r["citation_id"])
            # print(pm_pmc, pm_pmc.shape)
            # if pm_pmc.shape[0] > 0:
            #     pubtator_content = self.pubtator_api(pm_pmc)
            #     if pubtator_content:
            #         self.pubtator_analysis(pubtator_content, r)

        exit()

        lite_df_rs = (
            lite_df[["rs", "citation_id"]]
            .groupby("rs")
            .agg({"citation_id": ",".join})
            .reset_index()
        )
        lite_df_rs = pd.Series(
            lite_df_rs.citation_id.values, index=lite_df_rs.rs
        ).to_dict()
        lite_df["citation_id"] = lite_df["rs"].map(lite_df_rs)
        lite_df = lite_df.drop_duplicates(subset=["alleleid", "MAP", "citation_id"])
        lite_df = lite_df.dropna(subset=["citation_id"])
        # lite_df = lite_df.groupby(['MAP', 'PEXT_specific_tissues', 'H_tissues']).agg({'alleleid' : ','.join})

        # lite_df = lite_df.loc[lite_df['symbol'] == 'BIN1']

        print(lite_df.symbol.nunique())
        print(lite_df.MAP.nunique())
        print(lite_df.rs.nunique())
        print(
            len(
                set(
                    [
                        c
                        for e in lite_df.citation_id.values.tolist()
                        for c in e.split(",")
                    ]
                )
            )
        )
        print(lite_df[["alleleid", "MAP", "citation_id"]])
        lite_df.to_excel(output_dir + "protocole_variants_publi.xlsx")
        exit()

        for j, r in lite_df.iterrows():
            pm_pmc = self.convert_pmid_to_pmcid(r["citation_id"])
            # print(pm_pmc, pm_pmc.shape)
            if pm_pmc.shape[0] > 0:
                pubtator_content = self.pubtator_api(pm_pmc)
                if pubtator_content:
                    self.pubtator_analysis(pubtator_content, r)
        # 	exit()

        # pprint(lite_df.to_dict())

    def convert_pmid_to_pmcid(self, pmid):
        # df = pd.read_csv(
        # "https://www.ncbi.nlm.nih.gov/pmc/utils/idconv/v1.0/?tool=my_tool&email=my_email@example.com&ids={}&format=csv".format(
        # str(pmid)
        # ),
        # sep=",",
        # )
        # return df
        url = "https://www.ncbi.nlm.nih.gov/pmc/utils/idconv/v1.0/?tool=my_tool&email=my_email@example.com&ids={}&format=json".format(
            pmid
        )
        r = requests.get(url)
        j = r.json()
        pubmed_id = j["records"][0]["pmid"]
        if "pmcid" in j["records"][0]:
            pubmed_id = j["records"][0]["pmcid"]
        pubtator_cat = "pmc" if pubmed_id[0] == "P" else "pm"
        if pubtator_cat == "pmc":
            pubtator_url = "https://www.ncbi.nlm.nih.gov/research/pubtator-api/publications/export/biocjson?{}ids={}".format(
                pubtator_cat, pubmed_id
            )
            pubtator_json = requests.get(pubtator_url).json()
            pprint(pubtator_json)
            exit()

    @staticmethod
    def pubtator_api_rxiv(df_pm):
        # try:
        print(df_pm)
        df_pm["PMCID"] = df_pm["PMCID"].astype(str)
        df_pm["PMID"] = df_pm["PMID"].astype(str)
        pmc_ids = ",".join(df_pm.loc[df_pm["PMCID"] != "nan", "PMCID"].values.tolist())
        pm_ids = ",".join(df_pm.loc[df_pm["PMCID"] == "nan", "PMID"].values.tolist())
        empty = list()
        l_json = list()
        if pmc_ids:
            pm_pmc = "pmc"
            for tmp_id in pmc_ids.split(","):
                # print(tmp_id, pmc_ids)

                url = "https://www.ncbi.nlm.nih.gov/research/pubtator-api/publications/export/biocjson?{}ids={}".format(
                    str(pm_pmc), str(tmp_id)
                )
                r = requests.get(url)
                if r.text:
                    l_json.append(r.json())
                else:
                    empty.append(tmp_id)
        if pm_ids:
            pm_pmc = "pm"
            for tmp_id in pm_ids.split(","):
                # print(tmp_id, pm_ids)

                url = "https://www.ncbi.nlm.nih.gov/research/pubtator-api/publications/export/biocjson?{}ids={}".format(
                    str(pm_pmc), str(tmp_id)
                )
                r = requests.get(url)
                if r.text:
                    l_json.append(r.json())
                else:
                    empty.append(tmp_id)

        # print(len(l_json))
        # print(empty)
        return l_json
        # pprint(final_json)
        # exit()
        # except TypeError:
        # return None

    @staticmethod
    def pubtator_analysis(content, row):
        print(row)
        pext_tissues = eval(row["PEXT_specific_tissues"])
        tissues_decomposed = [
            sub_w.lower()
            for tissue in pext_tissues
            for w in tissue.split(" - ")
            for sub_w in w.split(" ")
        ]
        for elem in content:
            try:
                for passage in elem["passages"]:
                    for annotation in passage["annotations"]:
                        for tissue in tissues_decomposed:
                            # print(tissue, annotation['text'])
                            if tissue in annotation["text"].lower():
                                print(passage["text"])
                                pprint(annotation)
                                print(tissues_decomposed, tissue)
            except TypeError:
                print(elem)
        print("\n\n")

    @staticmethod
    def litvar_api(rsid):
        url = "https://www.ncbi.nlm.nih.gov/research/bionlp/litvar/api/v1/entity/litvar/rs{}%23%23".format(
            str(rsid)
        )
        r = requests.get(url)
        if r.text:
            print(rsid, r)
            pprint(r.json())

    def mapping_hpo(self, hpo, df_check_exon):

        hpo_pext_selected = (
            pd.merge(
                df_check_exon,
                hpo,
                left_on=["HGNC ID", "Exon"],
                right_on=["RefSeq_HGNC", "RefSeq_ranges"],
            )
            .sort_values(["symbol", "Exon", "variable"])
            .drop_duplicates()
        )

        hpo_pext_selected["PEXT_spec_tissues"] = hpo_pext_selected.apply(
            lambda r: [
                hpo_pext_selected.columns[j + 5]
                for j, col in enumerate(r[5:-12])
                if col == "True"
            ],
            axis=1,
        )

        hpo_pext_selected["HPO_groupped"] = hpo_pext_selected["HPO_groupped"].fillna(
            "None"
        )
        hpo_pext_selected = hpo_pext_selected[
            [
                "HGNC ID",
                "symbol",
                "ensg",
                "Exon",
                "HPO_POSITION",
                "HPO_ontology_accessions",
                "HPO_groupped",
                "PEXT_spec_tissues",
            ]
        ]
        hpo_pext_selected = hpo_pext_selected.reset_index(drop=True)


class Variants_CCRS:
    def __init__(
        self,
    ):
        # YAML FILES CONFIG
        yaml = load_config_file()
        self.exotic_files = yaml["EXOTIC"]

        self.tmp_ccrs_path = "/gstock/EXOTIC/data/VARIATIONS/CCRS_modified.parquet"
        self.tmp_variants_path = (
            "/gstock/EXOTIC/data/VARIATIONS/ClinVar_gnomAD_all_2021.parquet"
        )
        self.tmp_clinvar_ccrs_path = (
            "/gstock/EXOTIC/data/VARIATIONS/ClinVar_CCRS_gnomAD_2021.parquet"
        )
        self.tmp_phylocsf_path = (
            "/gstock/EXOTIC/data/CONSERVATION/phyloCSF_modified_2021.parquet"
        )
        self.gnomad_pathogenic_genes_path = "/gstock/biolo_datasets/variation/variation_sets/gnomAD/EXOME_MISTIC/RAW/gnomad_clinvar_genes.csv.gz"
        self.gnomad_pathogenic_genes = self.read_gnomad_file()

        # self.process_ccrs()
        self.process_clinvar()
        # self.merge_clinvar_ccrs()
        # self.process_phylocsf()
        exit()

        # TEST WITHOUT 0
        merge_df = merge_df.loc[merge_df["CCRS_CCR_percentile"] > 0]

        # MELT
        # merge_df = merge_df.melt(id_vars='Const_Alt', value_vars='CCRS_CCR_percentile')
        print(merge_df)

    def read_gnomad_file(self):
        gnomad_clinvar_genes = pd.read_csv(
            self.gnomad_pathogenic_genes_path,
            compression="gzip",
            sep="\t",
            low_memory=False,
        )
        # gnomad_clinvar_genes[
        #     ["CHROM", "POS", "REF", "ALT"]
        # ] = gnomad_clinvar_genes.ID.str.split("_", expand=True)
        # gnomad_clinvar_genes = gnomad_clinvar_genes.head(10000)
        gnomad_clinvar_genes["POS"] = gnomad_clinvar_genes["POS"].astype(int)
        gnomad_clinvar_genes["Genes"] = gnomad_clinvar_genes["Genes"].apply(eval)
        gnomad_clinvar_genes["ALT"] = gnomad_clinvar_genes["ALT"].apply(eval)
        gnomad_clinvar_genes = gnomad_clinvar_genes.explode("Genes")
        gnomad_clinvar_genes = gnomad_clinvar_genes.explode("ALT")
        gnomad_clinvar_genes = gnomad_clinvar_genes.loc[
            (gnomad_clinvar_genes["REF"].str.len() < 50)
            & (gnomad_clinvar_genes["ALT"].str.len() < 50)
        ]
        gnomad_clinvar_genes = gnomad_clinvar_genes.rename({"Genes": "Gene"}, axis=1)

        # gnomad_clinvar_genes["Status"] = "Benign"
        # gnomad_clinvar_genes = gnomad_clinvar_gen es.loc[
        #     (gnomad_clinvar_genes["ALT"].str.len() == 1)
        #     & (gnomad_clinvar_genes["REF"].str.len() == 1)
        # ]
        return gnomad_clinvar_genes

    def process_ccrs(self):
        if os.path.isfile(self.tmp_ccrs_path) is True:

            ## REFSEQ DETAILED
            refseq = pd.read_csv(
                self.exotic_files["refseq_path"],
                compression="gzip",
                sep="\t",
                low_memory=False,
            ).sort_values(by=["Gene", "ranges"])
            refseq = refseq.dropna(subset=["HGNC"])
            refseq["HGNC"] = refseq["HGNC"].astype(int)
            refseq["Ratio_num"] = refseq["Ratio"].apply(eval)
            refseq = refseq.loc[refseq["mRNA_nb"] > 1]
            # refseq = refseq.head(1000)

            # CCRS
            ccrs = pd.read_parquet(self.exotic_files["refseq_ccrs"])

            l_ccrs = list()
            for gene in tqdm(refseq.Gene.unique()):
                refseq_gene = refseq.loc[refseq["Gene"] == gene]
                ccrs_gene = ccrs.loc[ccrs["CCRS_Gene"] == gene]
                l_ccrs.append(self.mapping_ccrs(refseq_gene, ccrs_gene))

            print(pd.concat(l_ccrs, axis=0))
            exit()
            concat_ccrs = pd.concat(list(l_ccrs), axis=0)
            concat_ccrs = concat_ccrs.rename({"CCRS_Gene": "Gene"}, axis=1)
            merge_df = pd.merge(
                concat_ccrs,
                refseq[["Ratio_num", "ranges", "Gene"]],
                on=["Gene", "ranges"],
            )
            merge_df.loc[merge_df["Ratio_num"] == 1, "Const_Alt"] = "Const"
            merge_df.loc[merge_df["Ratio_num"] < 1, "Const_Alt"] = "Alt"
            bins = [0, 0.2, 0.4, 0.6, 0.8, 1]
            labels = bins.copy()
            labels_ratio = [
                str(round(labels[j], 1)) + " - " + str(round(labels[j + 1], 1))
                for j in range(len(labels) - 1)
            ]
            merge_df["Ratio_num_bins"] = pd.cut(
                merge_df["Ratio_num"],
                bins=bins,
                labels=labels_ratio,
                include_lowest=True,
            )
            merge_df.to_parquet(tmp_ccrs_path, index=False)

        else:
            merge_df = pd.read_parquet(self.tmp_ccrs_path)
        return merge_df

    def process_clinvar(self):
        if os.path.isfile(self.tmp_variants_path) is True:

            ## REFSEQ DETAILED
            refseq = pd.read_csv(
                self.exotic_files["refseq_path"],
                compression="gzip",
                sep="\t",
                low_memory=False,
            ).sort_values(by=["Gene", "ranges"])
            refseq = refseq.dropna(subset=["HGNC"])
            refseq["HGNC"] = refseq["HGNC"].astype(int)
            refseq["Ratio_num"] = refseq["Ratio"].apply(eval)
            refseq = refseq.loc[refseq["mRNA_nb"] > 1]
            # refseq = refseq.head(10000)

            ## CLINVAR
            clinvar = pd.read_parquet(self.exotic_files["clinvar_file_path"])
            # clinvar["POS"] = clinvar["VAR_ID"].apply(lambda r: r.split("_")[1])
            clinvar["POS"] = clinvar["POS"].astype(int)
            clinvar["MC"] = clinvar["MC"].apply(lambda r: r.split(",")[0])
            validated_clinvar = clinvar.loc[
                (clinvar["Status"] == "Pathogenic")
                & (~clinvar["Real_Status"].str.contains("onflicting") == True)
                & (clinvar["RS_STARS"] >= 1)
            ]
            # validated_clinvar = clinvar.loc[(~clinvar["Real_Status"].str.contains("onflicting") == True) & (clinvar["RS_STARS"] >= 1)]

            # validated_clinvar[
            # ["CHROM", "POS", "REF", "ALT"]
            # ] = validated_clinvar.VAR_ID.str.split("_", expand=True)
            validated_clinvar["ALT"] = validated_clinvar["ALT"].astype(str)
            validated_clinvar["ALT"] = validated_clinvar["ALT"].str.replace("\['", "")
            validated_clinvar["ALT"] = validated_clinvar["ALT"].str.replace("']", "")
            # validated_clinvar["POS"] = validated_clinvar["POS"].astype(int)
            validated_clinvar["ID"] = (
                validated_clinvar["CHROM"].astype(str)
                + "_"
                + validated_clinvar["POS"].astype(str)
                + "_"
                + validated_clinvar["REF"].astype(str)
                + "_"
                + validated_clinvar["ALT"].astype(str)
            )
            validated_clinvar = validated_clinvar.rename({"GENE": "Gene"}, axis=1)

            validated_clinvar = validated_clinvar.loc[
                (validated_clinvar["ALT"].str.len() == 1)
                & (validated_clinvar["REF"].str.len() == 1)
            ]
            validated_clinvar = (
                pd.merge(
                    validated_clinvar,
                    self.gnomad_pathogenic_genes[["Gene", "ID"]],
                    on=["Gene", "ID"],
                    how="outer",
                    indicator=True,
                )
                .query('_merge=="left_only"')
                .drop(["_merge"], axis=1)
            )

            l_variants = list()

            for gene in tqdm(refseq.Gene.unique()):
                refseq_gene = refseq.loc[refseq["Gene"] == gene]
                clinvar_gene = validated_clinvar.loc[validated_clinvar["Gene"] == gene]
                gnomad_gene = self.gnomad_pathogenic_genes.loc[
                    self.gnomad_pathogenic_genes["Gene"] == gene
                ]
                # if clinvar_gene.empty is False and gnomad_gene.empty is False:

                if clinvar_gene.empty is False:
                    clinvar_map = self.mapping_clinvar(refseq_gene, clinvar_gene)
                else:
                    clinvar_map = pd.DataFrame()

                if gnomad_gene.empty is False:
                    gnomad_map = self.mapping_clinvar(refseq_gene, gnomad_gene)
                else:
                    gnomad_map = pd.DataFrame()

                l_variants.append(pd.concat([clinvar_map, gnomad_map]))

                # if pd.concat(l_variants_clinvar).empty is False:
                # print(concat)
                # print(concat.loc[concat['Status'] == 'Benign'])
                # print(concat.loc[concat['Status'] == 'Pathogenic'])
                # exit()

            concat_variants = pd.concat(list(l_variants), axis=0)
            merge_df = pd.merge(
                concat_variants,
                refseq[["Ratio_num", "ranges", "Gene"]],
                on=["Gene", "ranges"],
            )
            merge_df.loc[merge_df["Ratio_num"] == 1, "Const_Alt"] = "Const"
            merge_df.loc[merge_df["Ratio_num"] < 1, "Const_Alt"] = "Alt"

            bins = [0, 0.2, 0.4, 0.6, 0.8, 1]
            labels = bins.copy()
            labels_ratio = [
                str(round(labels[j], 1)) + " - " + str(round(labels[j + 1], 1))
                for j in range(len(labels) - 1)
            ]
            merge_df["Ratio_num_bins"] = pd.cut(
                merge_df["Ratio_num"],
                bins=bins,
                labels=labels_ratio,
                include_lowest=True,
            )
            print(merge_df)

            merge_df.to_parquet(self.tmp_variants_path, index=False)

        else:
            merge_df = pd.read_parquet(self.tmp_variants_path)
        return merge_df

    def process_phylocsf(self):
        if os.path.isfile(self.tmp_phylocsf_path) is False:

            ## REFSEQ DETAILED
            refseq = pd.read_csv(
                self.exotic_files["refseq_path"],
                compression="gzip",
                sep="\t",
                low_memory=False,
            ).sort_values(by=["Gene", "ranges"])
            refseq = refseq.dropna(subset=["HGNC"])
            refseq["HGNC"] = refseq["HGNC"].astype(int)
            refseq["Ratio_num"] = refseq["Ratio"].apply(eval)
            refseq = refseq.loc[refseq["mRNA_nb"] > 1]
            # refseq = refseq.head(10000)

            # phyloCSF
            phyloCSF = pd.read_parquet(self.exotic_files["phylocsf_path"])

            l_phyloCSF = list()
            for gene in tqdm(refseq.Gene.unique()):
                refseq_gene = refseq.loc[refseq["Gene"] == gene]
                phyloCSF_gene = phyloCSF.loc[phyloCSF["Gene"] == gene]
                l_phyloCSF.append(self.mapping_phyloCSF(refseq_gene, phyloCSF_gene))

            # print(pd.concat(l_phyloCSF, axis=0))
            # exit()
            concat_phyloCSF = pd.concat(list(l_phyloCSF), axis=0)
            # concat_phyloCSF = concat_phyloCSF.rename({'CCRS_Gene' : 'Gene'}, axis=1)
            merge_df = pd.merge(
                concat_phyloCSF,
                refseq[["Ratio_num", "ranges", "Gene"]],
                on=["Gene", "ranges"],
            )
            merge_df.loc[merge_df["Ratio_num"] == 1, "Const_Alt"] = "Const"
            merge_df.loc[merge_df["Ratio_num"] < 1, "Const_Alt"] = "Alt"
            bins = [0, 0.2, 0.4, 0.6, 0.8, 1]
            labels = bins.copy()
            labels_ratio = [
                str(round(labels[j], 1)) + " - " + str(round(labels[j + 1], 1))
                for j in range(len(labels) - 1)
            ]
            merge_df["Ratio_num_bins"] = pd.cut(
                merge_df["Ratio_num"],
                bins=bins,
                labels=labels_ratio,
                include_lowest=True,
            )
            merge_df.to_parquet(self.tmp_phylocsf_path, index=False)

        else:
            merge_df = pd.read_parquet(self.tmp_phylocsf_path)
        return merge_df

    def merge_clinvar_ccrs(
        self,
    ):
        if os.path.isfile(self.tmp_clinvar_ccrs_path) is False:

            # CLINVAR
            clinvar = pd.read_parquet(self.tmp_variants_path)
            clinvar = clinvar.loc[clinvar["Status"] != "Other"]
            # clinvar = clinvar.head(100)

            # CCRS
            ccrs = pd.read_parquet(self.tmp_ccrs_path)
            bins = [0, 20, 80, 90, 95, 99, 100]
            labels = bins.copy()
            labels_ratio = [
                str(round(labels[j], 1)) + " - " + str(round(labels[j + 1], 1))
                for j in range(len(labels) - 1)
            ]
            ccrs["CCRS_bins"] = pd.cut(
                ccrs["CCRS_CCR_percentile"],
                bins=bins,
                labels=labels_ratio,
                include_lowest=True,
            )

            m = multiprocessing.Manager()
            l_ccrs = m.list()
            parmap.starmap(
                self.mapping_ccrs_clinvar,
                list(zip(clinvar.Gene.unique().tolist())),
                clinvar,
                ccrs,
                l_ccrs,
                pm_pbar=True,
            )
            # for gene in tqdm(clinvar.Gene.unique()):
            # clinvar_gene = clinvar.loc[clinvar['Gene'] == gene]
            # ccrs_gene = ccrs.loc[ccrs['Gene'] == gene]
            # self.mapping_ccrs_clinvar(clinvar_gene, ccrs_gene)

            concat_ccrs = pd.concat(list(l_ccrs), axis=0)
            # print(concat_ccrs)
            # exit()
            concat_ccrs.to_parquet(self.tmp_clinvar_ccrs_path, index=False)
        else:
            concat_ccrs = pd.read_parquet(self.tmp_clinvar_ccrs_path)
        return concat_ccrs

    def mp_mapping(self, gene, l_ccrs):
        refseq_gene = refseq.loc[refseq["Gene"] == gene]
        clinvar_gene = pathogenic_clinvar.loc[pathogenic_clinvar["GENE"] == gene]
        ccrs_gene = ccrs.loc[ccrs["CCRS_Gene"] == gene]
        # if ccrs_gene.empty is False:
        # print(ccrs_gene)
        tmp = self.mapping(refseq_gene, clinvar_gene, ccrs_gene)
        l_ccrs.append(tmp)

    def mapping_ccrs(self, refseq_gene, ccrs_gene):
        def map_ccrs(r, tmp_l):
            s = ccrs_gene["CCRS_Start"].between(
                int(r.split("-")[0]), int(r.split("-")[1]), inclusive=True
            )
            e = ccrs_gene["CCRS_End"].between(
                int(r.split("-")[0]), int(r.split("-")[1]), inclusive=True
            )
            concat = pd.concat([s, e], axis=1)
            concat = concat.loc[
                (concat["CCRS_Start"] == True) & (concat["CCRS_End"] == True)
            ]
            concat = ccrs_gene.loc[concat.index]
            concat["ranges"] = r
            tmp_l.append(concat)

        if refseq_gene.empty is False and ccrs_gene.empty is False:
            tmp_l = list()
            refseq_gene["ranges"].apply(lambda r: map_ccrs(r, tmp_l))
            return pd.concat(tmp_l)
        else:
            return pd.DataFrame()

    def mapping_phyloCSF(self, refseq_gene, phyloCSF_gene):
        def map_phyloCSF(r_raw, tmp_l):
            r = r_raw["ranges"]
            # print(r_raw)
            # print(phyloCSF_gene)
            s = phyloCSF_gene["Start"].between(
                int(r.split("-")[0]), int(r.split("-")[1]), inclusive=True
            )
            e = phyloCSF_gene["End"].between(
                int(r.split("-")[0]), int(r.split("-")[1]), inclusive=True
            )
            concat = pd.concat([s, e], axis=1)
            concat = concat.loc[(concat["Start"] == True) & (concat["End"] == True)]
            concat = phyloCSF_gene.loc[concat.index]
            concat["ranges"] = r
            tmp_l.append(concat)

        if refseq_gene.empty is False and phyloCSF_gene.empty is False:
            tmp_l = list()
            refseq_gene.apply(lambda r: map_phyloCSF(r, tmp_l), axis=1)
            return pd.concat(tmp_l)
        else:
            return pd.DataFrame()

    def mapping_ccrs_clinvar(self, gene, clinvar, ccrs, l_ccrs):
        def map_ccrs_clinvar(r, tmp_l):
            ranges = r["CCRS_ranges"]
            m = clinvar_gene["POS"].between(
                int(ranges.split("-")[0]), int(ranges.split("-")[1]), inclusive=True
            )
            if True in m.values.tolist():
                variants_between = clinvar_gene.loc[m.where(m == True).dropna().index]
                variants_between["CCRS_ranges"] = ranges
                variants_between["CCRS_bins"] = r["CCRS_bins"]
                variants_between["CCRS_CCR_percentile"] = r["CCRS_CCR_percentile"]
                tmp_l.append(variants_between)

        clinvar_gene = clinvar.loc[clinvar["Gene"] == gene]
        ccrs_gene = ccrs.loc[ccrs["Gene"] == gene]
        if clinvar_gene.empty is False and ccrs_gene.empty is False:
            tmp_l = list()
            ccrs_gene.apply(lambda r: map_ccrs_clinvar(r, tmp_l), axis=1)
            if tmp_l:
                l_ccrs.append(pd.concat(tmp_l))
            else:
                return pd.DataFrame()
        else:
            return pd.DataFrame()

    def mapping_clinvar(self, refseq_gene, clinvar_gene):
        def map_variant(r, tmp_l):
            m = clinvar_gene["POS"].between(
                int(r.split("-")[0]), int(r.split("-")[1]), inclusive=True
            )
            if True in m.values.tolist():
                variants_between = clinvar_gene.loc[m.where(m == True).dropna().index]
                variants_between["ranges"] = r
                tmp_l.append(variants_between)

        if refseq_gene.empty is False and clinvar_gene.empty is False:
            tmp_l = list()
            refseq_gene["ranges"].apply(lambda r: map_variant(r, tmp_l))
            if tmp_l:
                return pd.concat(tmp_l)
            else:
                return pd.DataFrame()

        else:
            return pd.DataFrame()


class GWAS:
    def __init__(self):
        yaml = load_config_file()
        self.exotic_files = yaml["EXOTIC"]
        self.gwas_files = yaml["GWAS"]

        gwas = self.load_gwas()
        variants = self.load_variants()
        variants = self.filter_variants_EXOTIC(variants)

        self.compare_variants_gwas(gwas, variants)

    def compare_variants_gwas(self, gwas, variants):
        pd.options.display.max_rows = 100
        merge = pd.merge(gwas, variants, left_on="SNPS", right_on="rs")
        print(merge)
        merge["P-VALUE"] = merge["P-VALUE"].astype(float)
        merge = merge.loc[merge["P-VALUE"] < 5 * 10e-8]
        print(merge)
        merge = merge.loc[merge["Parent term"] != "Other measurement"]
        print(merge)
        # merge = merge.loc[merge["Status"] == "Pathogenic"]
        # print(merge)
        merge = merge.loc[merge["OK_bronze"].str.len() > 0]
        print(merge)
        print(
            merge[
                [
                    "DISEASE/TRAIT",
                    "OK_bronze",
                    "Max",
                    "symbol",
                    "ranges",
                    "Ratio_num",
                    "Real_Status",
                    "VAR_ID",
                    "P-VALUE",
                ]
            ].drop_duplicates(subset=["DISEASE/TRAIT", "ranges"])
        )
        print(merge.symbol.nunique())
        print(merge.symbol.unique())

    def load_variants(self):

        variants = pd.read_parquet(self.exotic_files["clinvar_gnomad_pathogenic"])
        variants.loc[variants["Status"] == "Pathogenic", "rs"] = (
            "rs" + variants.loc[variants["Status"] == "Pathogenic", "rs"]
        )
        return variants

    def filter_variants_EXOTIC(self, variants):
        exotic = pd.read_parquet(self.exotic_files["exotic_ok_path"])
        exotic = exotic.loc[exotic["OK"].str.len() > 0]
        exotic["ranges"] = exotic["MAP"].apply(lambda r: r.split("_")[1])
        exotic_clinvar = pd.merge(
            exotic[
                [
                    "symbol",
                    "ranges",
                    "mRNA_nb",
                    "OK_bronze",
                    "OK_silver",
                    "OK_gold",
                    "OK",
                    "Max",
                    "pext_max",
                    "mean_proportion",
                ]
            ],
            variants,
            on="ranges",
        )
        return exotic_clinvar

    def load_gwas(self):
        gwas_associations = pd.read_csv(self.gwas_files["associations"], sep="\t")
        # gwas_studies = pd.read_csv(self.gwas_files['GWAS']['studies'], sep='\t')
        gwas_efo = pd.read_csv(self.gwas_files["efo"], sep="\t")
        gwas_efo_lite = gwas_efo[
            ["EFO term", "EFO URI", "Parent term", "Parent URI"]
        ].drop_duplicates()

        # print(gwas_associations)
        # print(gwas_efo)
        # print(gwas_efo_lite)
        gwas_return = pd.merge(
            gwas_associations,
            gwas_efo_lite,
            left_on="MAPPED_TRAIT_URI",
            right_on="EFO URI",
        )
        return gwas_return


class QTL:
    def __init__(self):

        ## YAML FILES CONFIG
        yaml = load_config_file()
        self.exotic_files = yaml

        self.gstock = "/gstock/biolo_datasets/variation/benchmark/Databases/"
        dict_qtl_files = self.files_qtl_gstock()

        #  for l in [dict_qtl_files["ieQTL_files_raw"], dict_qtl_files["isQTL_files_raw"]]:
        # pbar = tqdm(l)
        # for file in pbar:
        # pbar.set_description("Processing %s" % file.split("/")[-1])
        # self.process_qtl_grch37(file)

        merge_exotic_biomart = self.biomart_exotic()
        merge_exotic_biomart["Gene_ranges"] = (
            merge_exotic_biomart["start"].astype(str)
            + "-"
            + merge_exotic_biomart["end"].astype(str)
        )
        exotic_with_introns = self.retrieve_intron_refseq(merge_exotic_biomart)

        for l in [dict_qtl_files["sQTL_files_37"]]:
            # for intron_bool in [True, False]:
            t = l[0].split("/")[-2].split("_")[-1]
            concat_df_tissues = self.init_mp_merge_qtl_tissues(l, type_qtl=t)

            # if intron_bool is True:
            self.init_mp_mapping_qtl_exotic(
                type_qtl=t,
                qtl_files_37=l,
                exotic_file=merge_exotic_biomart.drop_duplicates(
                    subset=["Gene stable ID", "Gene name", "Gene_ranges"]
                ),
                condition="Gene",
            )
            # else:
            #     self.init_mp_mapping_qtl_exotic(
            #         type_qtl=t,
            #         qtl_files_37=l,
            #         exotic_file=merge_exotic_biomart,
            #         intron=intron_bool,
            #     )

            # self.init_mp_mapping_qtl_exotic(
            #     type_qtl=t,
            #     qtl_files_37=l,
            #     exotic_file=exotic_with_introns,
            #     intron=True,
            # )

        exit()

        # print(len(set(exotic_genes).intersection(set(qtl_genes))))
        # print(len(exotic_genes))
        # print(len(qtl_genes))

        # for i, gene in enumerate(common_genes):
        # if i == 5:
        # break

        # self.convert_to_vcf(concat_df_tissues, type_qtl="sQTL")

    def files_qtl_gstock(self):
        eqtl_files_raw = [
            self.gstock + "GTEX/V8/GTEx_Analysis_v8_eQTL/" + e
            for e in os.listdir(self.gstock + "GTEX/V8/GTEx_Analysis_v8_eQTL/")
            if "signif" in e
        ]
        sqtl_files_raw = [
            self.gstock + "GTEX/V8/GTEx_Analysis_v8_sQTL/" + e
            for e in os.listdir(self.gstock + "GTEX/V8/GTEx_Analysis_v8_sQTL/")
            if "signif" in e
        ]

        ieqtl_files_raw = [
            self.gstock + "GTEX/V8/Cell_type/GTEx_Analysis_v8_ieQTL/" + e
            for e in os.listdir(
                self.gstock + "GTEX/V8/Cell_type/GTEx_Analysis_v8_ieQTL/"
            )
        ]
        isqtl_files_raw = [
            self.gstock + "GTEX/V8/Cell_type/GTEx_Analysis_v8_isQTL/" + e
            for e in os.listdir(
                self.gstock + "GTEX/V8/Cell_type/GTEx_Analysis_v8_isQTL/"
            )
        ]

        eqtl_files_37 = [
            self.gstock + "GTEX/V8/GTEx_Analysis_v8_eQTL/" + e
            for e in os.listdir(self.gstock + "GTEX/V8/GTEx_Analysis_v8_eQTL/")
            if "final" in e
        ]

        sqtl_files_37 = [
            self.gstock + "GTEX/V8/GTEx_Analysis_v8_sQTL/" + e
            for e in os.listdir(self.gstock + "GTEX/V8/GTEx_Analysis_v8_sQTL/")
            if "final" in e
        ]

        ieqtl_files_37 = [
            self.gstock + "GTEX/V8/Cell_type/GTEx_Analysis_v8_ieQTL/" + e
            for e in os.listdir(
                self.gstock + "GTEX/V8/Cell_type/GTEx_Analysis_v8_ieQTL/"
            )
            if "final" in e
        ]
        isqtl_files_37 = [
            self.gstock + "GTEX/V8/Cell_type/GTEx_Analysis_v8_isQTL/" + e
            for e in os.listdir(
                self.gstock + "GTEX/V8/Cell_type/GTEx_Analysis_v8_isQTL/"
            )
            if "final" in e
        ]

        d = {
            "eQTL_files_raw": eqtl_files_raw,
            "sQTL_files_raw": sqtl_files_raw,
            "eQTL_files_37": eqtl_files_37,
            "sQTL_files_37": sqtl_files_37,
            "ieQTL_files_raw": ieqtl_files_raw,
            "isQTL_files_raw": isqtl_files_raw,
            "ieQTL_files_37": ieqtl_files_37,
            "isQTL_files_37": isqtl_files_37,
        }
        return d

    def init_mp_mapping_qtl_exotic(
        self, type_qtl, qtl_files_37, exotic_file, condition="Exon"
    ):
        type_region = condition

        if (
            os.path.isfile(
                "/gstock/EXOTIC/data/QTL/EXOTIC_{}_{}.parquet".format(
                    type_qtl, type_region
                )
            )
            is False
        ):

            concat_df_tissues = self.init_mp_merge_qtl_tissues(
                qtl_files_37, type_qtl=type_qtl
            )

            concat_df_tissues = concat_df_tissues.drop_duplicates()
            if type_qtl == "sQTL":
                concat_df_tissues["gene_id"] = concat_df_tissues["phenotype_id"].apply(
                    lambda r: r.split(":")[-1].split(".")[0]
                )

            concat_df_tissues["gene_id"] = concat_df_tissues["gene_id"].apply(
                lambda r: r.split(".")[0]
            )

            # print(concat_df_tissues)
            # exit()
            exotic_file = exotic_file.drop_duplicates(subset=["symbol", "MAP"])

            exotic_genes = list(exotic_file["Gene stable ID"].unique())
            qtl_genes = [
                e.split(".")[0] for e in list(concat_df_tissues.gene_id.unique())
            ]
            common_genes = list(sorted(set(exotic_genes).intersection(set(qtl_genes))))
            l = list()
            for i, gene in tqdm(enumerate(common_genes)):
                self.mp_mapping_qtl_exotic(
                    gene, exotic_file, concat_df_tissues, l, condition
                )

            concat_df = pd.concat(l)
            print(concat_df)

            concat_df.to_parquet(
                "/gstock/EXOTIC/data/QTL/EXOTIC_{}_{}.parquet".format(
                    type_qtl, type_region
                )
            )
        else:
            concat_df = pd.read_parquet(
                "/gstock/EXOTIC/data/QTL/EXOTIC_{}_{}.parquet".format(
                    type_qtl, type_region
                )
            )

        print(concat_df)
        print(concat_df.gene_id.nunique())
        print(condition)
        print(list(concat_df.columns))
        if condition == "Exon":
            print(concat_df.ranges.nunique())
        elif condition == "Intron":
            print(concat_df.Introns_ranges.nunique())
        elif condition == "Gene":
            print(concat_df.Gene_ranges.nunique())
        print(concat_df.variant_id.nunique())

        concat_df["OK"] = concat_df["OK"].apply(lambda r: r.split(","))

        concat_df["Check_tissue_QTL_EXOTIC"] = concat_df.apply(
            lambda r: True if r["Tissue"] in r["OK"] else False,
            axis=1,
        )
        test = concat_df.loc[concat_df["Check_tissue_QTL_EXOTIC"] == True]
        print(test)
        print(test.gene_id.nunique())
        if condition == "Exon":
            print(concat_df.ranges.nunique())
        elif condition == "Intron":
            print(test.Introns_ranges.nunique())
        elif condition == "Gene":
            print(test.Gene_ranges.nunique())
        print(test.variant_id.nunique())
        test.to_excel(
            "/gstock/EXOTIC/data/QTL/EXOTIC_{}_{}.xlsx".format(type_qtl, type_region),
            index=False,
        )

    def mp_mapping_qtl_exotic(
        self, gene, merge_exotic_biomart, concat_df_tissues, l, condition
    ):
        exotic_gene = merge_exotic_biomart.loc[
            merge_exotic_biomart["Gene stable ID"] == gene
        ]
        qtl_gene = concat_df_tissues.loc[concat_df_tissues["gene_id"] == gene]
        mapping_gene = self.map_qtl_to_exotic(exotic_gene, qtl_gene, condition)
        if mapping_gene.empty is False:
            l.append(mapping_gene)

    @staticmethod
    def map_qtl_to_exotic(exotic_gene, qtl_gene, condition):
        def map_variant(row, tmp_l):
            if condition == "Exon":
                r = row["Exon"]
            elif condition == "Intron":
                r = row["Introns_ranges"]
            elif condition == "Gene":
                r = row["Gene_ranges"]
            m = qtl_gene["POS_GRCh37"].between(
                int(r.split("-")[0]) - 5000, int(r.split("-")[1]) + 5000, inclusive=True
            )
            if True in m.values.tolist():

                variants_between = qtl_gene.loc[m.where(m == True).dropna().index]
                if condition == "Exon":
                    variants_between["ranges"] = r
                elif condition == "Intron":
                    variants_between["Introns_ranges"] = r
                elif condition == "Gene":
                    variants_between["Gene_ranges"] = r

                for col in ["OK_bronze", "OK_silver", "OK_gold", "OK"]:
                    variants_between[col] = ",".join(list(row[col]))
                tmp_l.append(variants_between)

        if exotic_gene.empty is False and qtl_gene.empty is False:
            tmp_l = list()
            exotic_gene.apply(lambda r: map_variant(r, tmp_l), axis=1)
            if tmp_l:
                return pd.concat(tmp_l)
            else:
                return pd.DataFrame()

        else:
            return pd.DataFrame()

    @staticmethod
    def compute_intron(df):
        l_introns = list()

        previous_start = 0
        previous_end = 0
        previous_intron = ""

        for j, row in df.iterrows():
            if j == 0:
                previous_end = row["End"]
                previous_start = row["Start"]

                l_introns.append("")

            elif j > 0 and j < len(df.index):
                intron = "-".join([str(previous_end + 1), str(row["Start"] - 1)])

                if row["Start"] == previous_end:
                    l_introns.append(previous_intron)

                elif row["Start"] != previous_end:
                    l_introns.append(intron)

                previous_end = row["End"]
                previous_start = row["Start"]
                previous_intron = intron

        df["Intron"] = l_introns

        return df

    def retrieve_intron_refseq(self, merge_exotic_biomart):

        if (
            os.path.isfile("/gstock/EXOTIC/data/EXOTIC/EXOTIC_with_introns.parquet")
            is False
        ):

            refseq = pd.read_csv(
                self.exotic_files["CCRS_analysis"]["RefSeq_transformed"],
                compression="gzip",
                sep="\t",
            )
            # print(refseq)

            d = dict()
            for j, gene in tqdm(enumerate(merge_exotic_biomart["Gene name"].unique())):
                # if j == 50:
                # break
                refseq_gene = refseq.loc[refseq["Gene"] == gene].reset_index(drop=True)
                exotic_gene = merge_exotic_biomart.loc[
                    merge_exotic_biomart["Gene name"] == gene
                ]
                refseq_gene = self.compute_intron(refseq_gene)

                exotic_in_refseq = refseq_gene.loc[
                    refseq_gene["ranges"].isin(exotic_gene["Exon"].values.tolist())
                ]
                # print(exotic_in_refseq)
                # print(refseq_gene)
                # intron_to_check_exotic = list()
                for exon_index in exotic_in_refseq.index:
                    if exon_index == 0:
                        intron_list_extended = [
                            refseq_gene.loc[exon_index + 1]["Intron"]
                        ]
                        d[
                            gene
                            + "_"
                            + str(exotic_in_refseq.loc[exon_index]["Start"])
                            + "-"
                            + str(exotic_in_refseq.loc[exon_index]["End"])
                        ] = [e for e in list(set(intron_list_extended)) if e != ""]

                    elif exon_index > 0 and exon_index < len(refseq_gene.index):
                        intron_list_extended = refseq_gene.loc[
                            exon_index : exon_index + 1
                        ]["Intron"].values.tolist()
                        d[
                            gene
                            + "_"
                            + str(exotic_in_refseq.loc[exon_index]["Start"])
                            + "-"
                            + str(exotic_in_refseq.loc[exon_index]["End"])
                        ] = [e for e in list(set(intron_list_extended)) if e != ""]

                    elif exon_index == len(refseq_gene.index):
                        intron_list_extended = [refseq_gene.loc[exon_index]["Intron"]]
                        d[
                            gene
                            + "_"
                            + str(exotic_in_refseq.loc[exon_index]["Start"])
                            + "-"
                            + str(exotic_in_refseq.loc[exon_index]["End"])
                        ] = [e for e in list(set(intron_list_extended)) if e != ""]

                # introns_full_exotic = [
                # e for e in sorted(list(set(intron_to_check_exotic))) if e != ""
                # ]
                # d[list(exotic_gene["Gene name"].unique())[0]] = introns_full_exotic

            merge_exotic_biomart["Introns_ranges"] = merge_exotic_biomart["MAP"].map(d)
            exotic_with_introns = merge_exotic_biomart.loc[
                merge_exotic_biomart["Introns_ranges"].str.len() > 0
            ]

            exotic_with_introns = exotic_with_introns.explode("Introns_ranges")
            exotic_with_introns.to_parquet(
                "/gstock/EXOTIC/data/EXOTIC/EXOTIC_with_introns.parquet"
            )
        else:
            exotic_with_introns = pd.read_parquet(
                "/gstock/EXOTIC/data/EXOTIC/EXOTIC_with_introns.parquet"
            )
        return exotic_with_introns

    def biomart_exotic(self):
        if (
            os.path.isfile("/gstock/EXOTIC/data/EXOTIC/EXOTIC_biomart_for_QTL.parquet")
            is False
        ):
            exotic_ok = pd.read_parquet(self.exotic_files["EXOTIC"]["exotic_ok_path"])
            exotic_ok = exotic_ok.loc[exotic_ok["OK"].str.len() > 0]

            exotic_ok["Exon"] = exotic_ok["MAP"].apply(lambda r: r.split("_")[1])

            exotic_ok["Start"] = exotic_ok["Exon"].apply(lambda r: r.split("-")[0])
            exotic_ok["Start"] = exotic_ok["Start"].astype(int)

            exotic_ok["End"] = exotic_ok["Exon"].apply(lambda r: r.split("-")[1])
            exotic_ok["End"] = exotic_ok["End"].astype(int)

            print(exotic_ok)

            biomart = pd.read_csv(
                self.exotic_files["EXOTIC"]["biomart_ensembl_hgnc"],
                sep="\t",
                compression="gzip",
            )

            biomart["Strand"] = biomart["Strand"].replace({1: "+", -1: "-"})
            biomart = biomart.rename(
                {
                    "Strand": "strand",
                    "Chromosome/scaffold name": "chrom",
                    "Gene start (bp)": "start",
                    "Gene end (bp)": "end",
                },
                axis=1,
            )
            biomart["start"] = biomart["start"].astype(int)
            biomart["end"] = biomart["end"].astype(int)

            biomart_to_convert = biomart[
                [
                    "chrom",
                    "start",
                    "end",
                    "strand",
                    "Gene stable ID",
                    "Gene name",
                ]
            ]

            biomart_to_convert = biomart_to_convert.sort_values(by=["chrom", "start"])

            biomart_to_convert.to_csv(
                self.exotic_files["EXOTIC"]["biomart_ensembl_hgnc"].replace(
                    ".txt.gz", "_to_convert.bed3"
                ),
                # compression="gzip",
                sep="\t",
                index=False,
                # header=False,
            )

            biomart_grch37 = pd.read_csv(
                self.exotic_files["EXOTIC"]["biomart_ensembl_hgnc"].replace(
                    ".txt.gz", "_GRCh37.bed3"
                ),
                sep="\t",
            )
            biomart_grch37.columns = biomart_to_convert.columns

            merge_exotic_biomart = pd.merge(
                biomart_grch37, exotic_ok, left_on="Gene name", right_on="symbol"
            )

            merge_exotic_biomart.to_parquet(
                "/gstock/EXOTIC/data/EXOTIC/EXOTIC_biomart_for_QTL.parquet"
            )
        else:
            merge_exotic_biomart = pd.read_parquet(
                "/gstock/EXOTIC/data/EXOTIC/EXOTIC_biomart_for_QTL.parquet"
            )

        return merge_exotic_biomart

    def handle_duplicated_genes(self, biomart, merge_exotic_biomart):
        def check_bounds(r):
            if r["Start"] >= r["start"] and r["End"] <= r["end"]:
                return True
            else:
                return False

        merge_exotic_biomart_genes = list(merge_exotic_biomart.symbol.unique())

        duplicated_genes = biomart.loc[
            (biomart["Gene name"].isin(merge_exotic_biomart_genes))
            & (biomart["Gene name"].duplicated() == True),
        ]
        print(duplicated_genes)
        exit()

        merge_exotic_biomart_duplicated = merge_exotic_biomart.loc[
            merge_exotic_biomart["symbol"].isin(
                duplicated_genes["Gene name"].unique().tolist()
            )
        ]

        merge_exotic_biomart_duplicated[
            "Check_exon_bounds"
        ] = merge_exotic_biomart_duplicated.apply(check_bounds, axis=1)
        print(merge_exotic_biomart_duplicated)
        print(
            merge_exotic_biomart_duplicated.loc[
                merge_exotic_biomart_duplicated["Check_exon_bounds"] == False
            ]
        )

        print(
            merge_exotic_biomart_duplicated.loc[
                merge_exotic_biomart_duplicated["Check_exon_bounds"] == False, "symbol"
            ].nunique()
        )
        print(
            merge_exotic_biomart_duplicated.loc[
                merge_exotic_biomart_duplicated["Gene name"] == "MCTS1"
            ]
        )

        print(biomart.loc[biomart["Gene name"] == "MCTS1"])
        exit()

        for gene in duplicated_exotic.symbol.unique():
            print(gene)
            duplicated_exotic_gene = duplicated_exotic.loc[
                duplicated_exotic["symbol"] == gene
            ]
            duplicated_biomart_gene = biomart.loc[biomart["Gene name"] == gene]
            print(duplicated_exotic_gene)
            print(duplicated_biomart_gene)
            min_exon = duplicated_exotic_gene["Start"].min()
            max_exon = duplicated_exotic_gene["End"].max()
            print(min_exon)
            print(max_exon)
            exit()

    def init_mp_merge_qtl_tissues(
        self,
        files_37,
        type_qtl="eQTL",
    ):
        if "i" in type_qtl:
            o = (
                self.gstock
                + "GTEX/V8/Cell_type/GTEx_Analysis_v8_{}/".format(type_qtl)
                + "Merge_tissues_{}_signif_variant_GRCh37.parquet".format(type_qtl)
            )
        else:
            o = (
                self.gstock
                + "GTEX/V8/GTEx_Analysis_v8_{}/".format(type_qtl)
                + "Merge_tissues_{}_signif_variant_GRCh37.parquet".format(type_qtl)
            )
        if os.path.isfile(o) is False:
            m = multiprocessing.Manager()
            l = m.list()
            parmap.starmap(
                self.mp_merge_qtl_tissues, list(zip(files_37)), l, pm_pbar=True
            )

            concat_tissues_df = pd.concat(l)

            concat_tissues_df["CHROM"] = concat_tissues_df["CHROM"].astype(str)

            concat_tissues_df.to_parquet(o)
        else:
            concat_tissues_df = pd.read_parquet(o)
        return concat_tissues_df

    def mp_merge_qtl_tissues(self, file, l):
        df = pd.read_csv(
            file,
            compression="gzip",
            sep="\t",
            # nrows=100000,
            low_memory=False,
        )
        df = df.drop(["POS_GRCh37.1"], axis=1)
        df["Tissue"] = file.split(".")[0].split("/")[-1]
        df["Cell_type"] = file.split(".")[1].split("/")[-1]
        l.append(df)
        # return df

    def convert_to_vcf(self, df, type_qtl="eQTL"):
        if type_qtl == "eQTL":
            cols = [
                "Tissue",
                "ID_GRCh38",
                "gene_id",
                "tss_distance",
                "ma_samples",
                "ma_count",
                "maf",
                "pval_nominal",
                "slope",
                "slope_se",
                "pval_nominal_threshold",
                "min_pval_nominal",
                "pval_beta",
            ]
        elif type_qtl == "sQTL":
            cols = [
                "Tissue",
                "ID_GRCh38",
                "phenotype_id",
                "tss_distance",
                "ma_samples",
                "ma_count",
                "maf",
                "pval_nominal",
                "slope",
                "slope_se",
                "pval_nominal_threshold",
                "min_pval_nominal",
                "pval_beta",
            ]

        for col in cols:
            df[col] = col + "=" + df[col].astype(str)

        df["INFO"] = df[cols].apply(
            lambda row: ";".join(row.values.astype(str)), axis=1
        )

        df = df.rename({"POS_GRCh37": "POS", "CHROM": "#CHROM"}, axis=1)

        df["ID"] = "."
        df["QUAL"] = 40
        df["FILTER"] = "."

        df = df[["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"]]
        # df = pd.concat(
        #     [
        #         pd.DataFrame(
        #             [
        #                 "##fileformat=VCFv4.1",
        #                 "##fileDate=20090805",
        #                 "##source=myImputationProgramV3.1",
        #                 "##reference=file:///seq/references/",
        #             ]
        #         ),
        #         df,
        #     ],
        #     axis=0,
        # )

        header = "##fileformat=VCFv4.1\n##fileDate=20090805\n##source=myImputationProgramV3.1\n##reference=file:///seq/references/\n"

        print(df)

        output_VCF = "GTEx_V8_{}_GRCh37_transformed.vcf.gz".format(type_qtl)
        with gzip.open(
            self.gstock + "GTEX/V8/GTEx_Analysis_v8_{}/".format(type_qtl) + output_VCF,
            "wb",
        ) as vcf:
            vcf.write(header.encode())

        df.to_csv(
            self.gstock + "GTEX/V8/GTEx_Analysis_v8_{}/".format(type_qtl) + output_VCF,
            sep="\t",
            mode="a",
            index=False,
            compression="gzip",
        )

    def process_qtl_grch37(self, file):

        df = pd.read_csv(
            file,
            compression="gzip",
            sep="\t",
        )

        df["ID_GRCh38"] = df.variant_id.str.replace("_b38", "")
        df[["CHROM", "POS", "REF", "ALT"]] = df["ID_GRCh38"].str.split("_", expand=True)
        df["CHROM"] = df["CHROM"].str.replace("chr", "")

        df_output = df[
            [
                "CHROM",
                "POS",
                "POS",
                "REF",
                "ALT",
                "variant_id",
            ]
        ]

        if os.path.isfile(file.replace(".txt.gz", "_modified_38.txt.gz")) is False:

            df_output.to_csv(
                file.replace(".txt.gz", "_modified_38.txt.gz"),
                sep="\t",
                index=False,
                header=False,
            )

        if os.path.isfile(file.replace(".txt.gz", "_modified_37.txt.gz")) is False:
            subprocess.call(
                [
                    "/home/weber/.conda/envs/ExoCarto/bin/python",
                    "/home/weber/.conda/envs/ExoCarto/bin/CrossMap.py",
                    "bed",
                    self.gstock + "CrossMap/GRCh38_to_GRCh37.chain.gz",
                    file.replace(".txt.gz", "_modified_38.txt.gz"),
                    file.replace(".txt.gz", "_modified_37.txt"),
                ],
                shell=False,
            )

            subprocess.call(
                [
                    "bgzip",
                    file.replace(".txt.gz", "_modified_37.txt"),
                ],
                shell=False,
            )

        if os.path.isfile(file.replace(".txt.gz", "_final_37.txt.gz")) is False:
            lite_df_37 = pd.read_csv(
                file.replace(".txt.gz", "_modified_37.txt" + ".gz"),
                sep="\t",
                compression="gzip",
            )

            lite_df_37.columns = df_output.columns
            lite_df_37 = lite_df_37.rename({"POS": "POS_GRCh37"}, axis=1)
            lite_df_37 = lite_df_37[["variant_id", "CHROM", "POS_GRCh37", "REF", "ALT"]]

            new_df = pd.merge(
                df.drop(["CHROM", "POS", "REF", "ALT"], axis=1),
                lite_df_37,
                on="variant_id",
            )
            if "eqtl" in file:

                new_df = new_df[
                    [
                        "ID_GRCh38",
                        "CHROM",
                        "POS_GRCh37",
                        "REF",
                        "ALT",
                        "variant_id",
                        "gene_id",
                        "tss_distance",
                        "ma_samples",
                        "ma_count",
                        "maf",
                        "pval_nominal",
                        "slope",
                        "slope_se",
                        "pval_nominal_threshold",
                        "min_pval_nominal",
                        "pval_beta",
                    ]
                ]
            elif "sqtl" in file:
                new_df = new_df[
                    [
                        "ID_GRCh38",
                        "CHROM",
                        "POS_GRCh37",
                        "REF",
                        "ALT",
                        "variant_id",
                        "phenotype_id",
                        "tss_distance",
                        "ma_samples",
                        "ma_count",
                        "maf",
                        "pval_nominal",
                        "slope",
                        "slope_se",
                        "pval_nominal_threshold",
                        "min_pval_nominal",
                        "pval_beta",
                    ]
                ]
            elif "ieQTL" in file or "isQTL" in file:
                new_df = new_df[
                    [
                        "ID_GRCh38",
                        "CHROM",
                        "POS_GRCh37",
                        "REF",
                        "ALT",
                        "variant_id",
                        "gene_id",
                        "gene_name",
                        "biotype",
                        "phenotype_id",
                        "tss_distance",
                        "maf",
                        "ma_samples",
                        "ma_count",
                        "pval_g",
                        "b_g",
                        "b_g_se",
                        "pval_i",
                        "b_i",
                        "b_i_se",
                        "pval_gi",
                        "b_gi",
                        "b_gi_se",
                        "pval_emt",
                        "tests_emt",
                        "pval_adj_bh",
                    ]
                ]
            new_df.to_csv(
                file.replace(".txt.gz", "_final_GRCh37.txt.gz"),
                compression="gzip",
                sep="\t",
                index=False,
            )
        # exit()


if __name__ == "__main__":

    # GTEx_TPM = Compute_GTEx_profile()
    # SearchClass = EXOTIC(gtex_profile=GTEx_TPM.v_profile, cpus=20, omim_detailed=True)
    # SearchClass = EXOTIC(gtex_profile="", cpus=2, omim_detailed=True)
    # V = Variants_CCRS()
    # G = GWAS()
    Q = QTL()

""" 
    pprint(SearchClass.concat_df.loc[SearchClass.concat_df['symbol'] == 'CACNB2'].to_dict())
    pprint(GTEx_TPM.h_profile.loc[GTEx_TPM.h_profile['symbol'] == 'CACNB2'].to_dict())
    pprint(GTEx_TPM.v_profile.loc[GTEx_TPM.v_profile['symbol'] == 'CACNB2'].to_dict())

    total_genes = output_df.loc[output_df["check_exon"]
                                == True, "symbol"].nunique()
    d_genes = len(
        set(
            output_df.loc[
                (output_df["check_exon"] == True)
                & (output_df["HGNC ID"].isin(disease_genes))
            ]
            .drop_duplicates()["symbol"]
            .tolist()
        )
    )
    h_genes = len(
        set(
            output_df.loc[
                (output_df["check_exon"] == True)
                & (~output_df["HGNC ID"].isin(disease_genes))
            ]
            .drop_duplicates()["symbol"]
            .tolist()
        )
    )

    output_df.to_csv('/home/weber/PycharmProjects/ExoCarto/data/test/test_pext_frequent.csv.gz', compression='gzip', sep='\t', index=False)
    output_df.to_parquet('/home/weber/PycharmProjects/ExoCarto/data/test/test_pext.parquet', index=False)

    hpo_pext_selected
    hpo_genes = len(dict(collections.Counter(hpo_pext_selected.symbol.values)))
    results.append(
    	{
    					**{
    						'Mode': mode,
    						'Total genes': total_genes,
    						'D genes': d_genes,
    						'H genes': h_genes,
    						'ClinVar genes': clinvar_genes,
    						'HPO genes': hpo_genes,
    					},
    					**config_search[mode]
    	}

    )

    results_df = pd.DataFrame(results)
    print(results_df)
    results_df.to_excel('/home/weber/PycharmProjects/ExoCarto/data/test/test_stats.xlsx')
 """
