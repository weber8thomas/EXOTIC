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

import json
from pandarallel import pandarallel

pandarallel.initialize(nb_workers=60, progress_bar=True)
from pprint import pprint

# Custom utils
sys.path.append("/home/weber/PycharmProjects/EXOTIC/src")
from utils.utils import load_config_file

## YAML FILES CONFIG
yaml = load_config_file(config_file="/home/weber/PycharmProjects/EXOTIC/src/config.yaml")


def load_omim(omim_path, biomart_omim_path):
    omim = pd.read_csv(
        omim_path,
        compression="gzip",
        sep="\t",
    )

    # print("Total : ", omim.OMIM.nunique())
    omim = omim.dropna(subset=list(omim.columns[6:-2]), how="all")
    # print("Dropna on all cols : " omim.OMIM.nunique())

    biomart_omim = pd.read_csv(biomart_omim_path, sep="\t", compression="gzip").dropna(subset=["MIM gene accession"])
    biomart_omim["MIM gene accession"] = biomart_omim["MIM gene accession"].astype(int)
    biomart_omim = biomart_omim.rename({"MIM gene accession": "OMIM", "Gene stable ID": "ensg", "Gene name": "Name"}, axis=1)

    ## ADD Gene TO OMIM
    omim = pd.merge(biomart_omim[["ensg", "Name", "OMIM"]], omim, on="OMIM")
    # print("Merge BIOMART : " omim.OMIM.nunique())

    # MELT
    omim = omim.melt(id_vars=list(omim.columns)[:7], value_vars=list(omim.columns)[7:], var_name="OMIM_BP", value_name="OMIM_BP_phenotypes").dropna()

    return omim


def load_exotic(exotic_path, min_max, cutoff=0.8):
    exotic = pd.read_parquet(exotic_path)
    exotic = exotic.loc[exotic["EXOTIC_{}".format(min_max)] >= cutoff]
    if min_max == "up":
        exotic["EXOTIC_tissues_above_cutoff_{}".format(min_max)] = exotic[exotic.columns[9:-6]].apply(
            lambda r: [exotic.columns[9:-6][j] for j, e in enumerate(r) if e >= cutoff], axis=1
        )
    else:
        exotic["EXOTIC_tissues_above_cutoff_{}".format(min_max)] = exotic[exotic.columns[9:-6]].apply(
            lambda r: [exotic.columns[9:-6][j] for j, e in enumerate(r) if (1 - e) >= cutoff], axis=1
        )
    # print(exotic.columns[9:-6])
    # print(exotic)
    # exit()
    exotic = exotic.explode("EXOTIC_tissues_above_cutoff_{}".format(min_max))
    return exotic


def compare_to_omim_basic_exon_level(exotic, min_max, omim, omim_detailed=True):
    print(min_max)
    print(exotic.symbol.nunique())
    print(exotic.MAP.nunique())

    mapping_omim_gtex = dicts["mapping_omim_gtex_detailed"]

    omim = omim.where(pd.notnull(omim), None)

    omim = omim.loc[~omim.duplicated(keep="last", subset=["OMIM", "Pheno_OMIM", "Pheno_prefered_title"])]

    print(omim)

    print(exotic)

    merge = pd.merge(omim, exotic, on="ensg")
    merge["EXOTIC_tissue_BP"] = merge["EXOTIC_tissues_above_cutoff_{}".format(min_max)].map(mapping_omim_gtex)
    merge = merge.explode("EXOTIC_tissue_BP")
    merge = merge[
        [
            "ensg",
            "OMIM",
            "OMIM_BP",
            "OMIM_BP_phenotypes",
            "EXOTIC_tissues_above_cutoff_{}".format(min_max),
            "EXOTIC_{}".format(min_max),
            "EXOTIC_bins_{}".format(min_max),
            "EXOTIC_tissue_BP",
        ]
    ]

    print(merge)
    print(merge.loc[merge["OMIM_BP"] == merge["EXOTIC_tissue_BP"]])

    exit()

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
                            if omim_body_part in omim_tmp and omim_tmp[omim_body_part] is not None:
                                threshold = [e.replace("OK_", "") for e in ["OK_bronze", "OK_silver", "OK_gold"] if t in row[e]][0]

                                omim_associations.append(
                                    {
                                        "symbol": row["symbol"],
                                        "MAP": row["MAP"],
                                        "Corrected_Tissues": row["OK"],
                                        # "H_tissues": row["H_tissues"],
                                        "OMIM_pheno_exon_level": pheno["Pheno_OMIM"],
                                        "Threshold": threshold,
                                        "EXOTIC_value": row[t + "_exotic"],
                                        "OMIM_pheno_name_exon_level": pheno["Pheno_Name"],
                                        "Tissue_exon_level": t,
                                        "OMIM_body_part_exon_level": omim_body_part,
                                        "OMIM_details_exon_level": [e.split(" {")[0] for e in eval(omim_tmp[omim_body_part])[0].splitlines()],
                                    }
                                )
    omim_associations = pd.DataFrame(omim_associations)

    print(omim_associations)
    print(omim_associations.symbol.nunique())
    print(omim_associations.MAP.nunique())
    print(omim_associations.loc[omim_associations["Threshold"] == "bronze", "symbol"].nunique())
    print(omim_associations.loc[omim_associations["Threshold"] == "silver", "symbol"].nunique())
    print(omim_associations.loc[omim_associations["Threshold"] == "gold", "symbol"].nunique())
    print(omim_associations.loc[omim_associations["Threshold"] == "bronze", "MAP"].nunique())
    print(omim_associations.loc[omim_associations["Threshold"] == "silver", "MAP"].nunique())
    print(omim_associations.loc[omim_associations["Threshold"] == "gold", "MAP"].nunique())
    omim_associations.to_excel(
        output_dir + "PHENOTYPES/" + "table_OMIM_exon_level_matrix_corrected.xlsx",
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
    #     output_dir + "PHENOTYPES/" + "table_OMIM_gene_exon_level_matrix.xlsx",
    #     index=False,
    # )
    # undupl[["Tissue_exon_level", "OMIM_body_part_exon_level"]].groupby(
    #     "Tissue_exon_level"
    # ).size().to_excel(
    #     output_dir + "PHENOTYPES/" + "table_OMIM_exon_level_tissues.xlsx"
    # )
    # undupl[["Tissue_exon_level", "OMIM_body_part_exon_level"]].groupby(
    #     "OMIM_body_part_exon_level"
    # ).size().to_excel(
    #     output_dir + "PHENOTYPES/" + "table_OMIM_exon_level_BP.xlsx"
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
    #     output_dir + "PHENOTYPES/" + "table_OMIM_gene_exon_level_tissues.xlsx"
    # )
    # undupl[["Tissue_gene_level", "OMIM_body_part_gene_level"]].groupby(
    #     "OMIM_body_part_gene_level"
    # ).size().to_excel(
    #     output_dir + "PHENOTYPES/" + "table_OMIM_gene_exon_level_BP.xlsx"
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


if __name__ == "__main__":

    ## YAML FILES CONFIG
    files = load_config_file(config_file="src/config_clean.yaml")
    output_dir = "/gstock/EXOTIC/data/"
    # pprint(exotic_files)

    ## JSON DICTS CONFIG
    dicts = json.load(open("src/EXOTIC_config.json"))

    print(files["EXOTIC"]["exotic_modified_zscore"])
    exotic_up = load_exotic(files["EXOTIC"]["exotic_modified_zscore"], min_max="up", cutoff=0.7)
    exotic_down = load_exotic(files["EXOTIC"]["exotic_modified_zscore"], min_max="down", cutoff=0.7)
    print(exotic_down)

    omim = load_omim(files["EXOTIC"]["omim_detailed"], files["BIOMART"]["biomart_omim"])

    # print(omim)

    compare_to_omim_basic_exon_level(exotic_up, min_max="up", omim=omim)
