# IMPORTS

import math
import os
import sys
import pandas as pd

# sys.path.append("/home/weber/PycharmProjects/EXOTIC/src")
sys.path.append("/home/weber/PycharmProjects/EXOTIC/clean/src")

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
import seaborn as sns
import json
from utils import utils


## YAML FILES CONFIG
yaml = utils.load_config_file(config_file="/home/weber/PycharmProjects/EXOTIC/clean/src/config_clean_clean.yaml")


class ProcessRefSeq:
    def __init__(self, path):
        """[Main function to launch steps]

        Arguments:
            path {[str]} -- [Output file path]

        Returns:
            [pd.DataFrame] -- [Final processed refseq dataframe]
        """

        if os.path.isfile(path) is False:
            print("# Files don't exist ☒")

            # * 0 Load raw file
            refseq_grch37_gff = yaml["1_GENOMICS"]["External"]["raw_refseq"]
            refseq_df = self.load_refseq(refseq_grch37_gff)

            utils.mkdir(os.path.dirname(yaml["1_GENOMICS"]["TMP"]["tmp_refseq_chroms"]))

            # * 1 Build tmp files by category
            refseq_df_chroms = self.refseq_chroms_fct(yaml["1_GENOMICS"]["TMP"]["tmp_refseq_chroms"], refseq_df)
            refseq_df_pc_genes = self.refseq_pc_genes_fct(yaml["1_GENOMICS"]["TMP"]["tmp_refseq_pc_genes"], refseq_df)
            refseq_df_mrnas = self.refseq_mrnas_fct(yaml["1_GENOMICS"]["TMP"]["tmp_refseq_mrnas"], refseq_df)
            refseq_df_exons = self.refseq_exons_fct(yaml["1_GENOMICS"]["TMP"]["tmp_refseq_exons"], refseq_df)
            refseq_df_cds = self.refseq_cds_fct(yaml["1_GENOMICS"]["TMP"]["tmp_refseq_cds"], refseq_df)

            # * 2 Select exons(with UTRs) & coding exons
            self.select_exons(
                yaml["1_GENOMICS"]["TMP"]["tmp_refseq_exons_filtered"],
                refseq_df_chroms,
                refseq_df_mrnas,
                refseq_df_pc_genes,
                refseq_df_exons,
                "Exon",
            )
            self.select_exons(
                yaml["1_GENOMICS"]["TMP"]["tmp_refseq_cds_filtered"],
                refseq_df_chroms,
                refseq_df_mrnas,
                refseq_df_pc_genes,
                refseq_df_cds,
                "CDS",
            )

            # * 3 Concat tmp files
            concat_exons_cds = self.concat_exons_cds_fct(
                yaml["1_GENOMICS"]["TMP"]["tmp_refseq_concat_exons_cds_filtered"],
                yaml["1_GENOMICS"]["TMP"]["tmp_refseq_exons_filtered"],
                yaml["1_GENOMICS"]["TMP"]["tmp_refseq_cds_filtered"],
            )

            # * 4 Groupby and produce final output file
            self.groupby_mrnas = self.groupby_mrnas_fct(yaml["1_GENOMICS"]["Final"]["refseq_processed"], concat_exons_cds)

            # * 5 Groupby and produce final output file
            self.final_df = self.get_shared_exons(yaml["1_GENOMICS"]["Final"]["refseq_cds_with_variable"], self.groupby_mrnas)

            # * 6 Compute CDS combination and retrieve specificies linked to UTRs
            utrs = self.compute_utrs(yaml["1_GENOMICS"]["Final"]["refseq_miso_utrs"], self.groupby_mrnas)

        else:
            print("# Files exist ✓, Loading ... ")

            # print(pd.read_parquet(yaml["1_GENOMICS"]["TMP"]["tmp_refseq_pc_genes"]))
            # exit()

            # # * 3 Concat tmp files
            # concat_exons_cds = self.concat_exons_cds_fct(
            #     yaml["1_GENOMICS"]["TMP"]["tmp_refseq_concat_exons_cds_filtered"],
            #     yaml["1_GENOMICS"]["TMP"]["tmp_refseq_exons_filtered"],
            #     yaml["1_GENOMICS"]["TMP"]["tmp_refseq_cds_filtered"],
            # )

            # # * 4 Groupby and produce final output file
            # self.groupby_mrnas = self.groupby_mrnas_fct(yaml["1_GENOMICS"]["Final"]["refseq_processed"], concat_exons_cds)
            # print(self.groupby_mrnas)
            # exit()

            self.final_df = pd.read_parquet(yaml["1_GENOMICS"]["Final"]["refseq_cds_with_variable"])
            utrs = pd.read_parquet(yaml["1_GENOMICS"]["Final"]["refseq_miso_utrs"])
            print(self.final_df)
            print(utrs)

    @staticmethod
    def load_refseq(path):
        """[Load RefSeq GFF file]

        Arguments:
            path {[str]} -- [Path to the GFF RefSeq file]

        Returns:
            [pd.DataFrame] -- [RefSeq GFF turned into pandas dataframe]
        """
        print("### Load RefSeq / File = {}".format(path))

        refseq_df = pd.read_csv(
            path,
            compression="gzip",
            sep="\t",
            skiprows=9,
            nrows=10000,
            names=["NC", "RefSeq_validation", "Region_type", "Start", "End", "Score", "Strand", "Phase", "Attributes"],
        )
        refseq_df = refseq_df.dropna(subset=["Start", "End"])
        refseq_df["Start"] = refseq_df["Start"].astype(int)
        refseq_df["End"] = refseq_df["End"].astype(int)
        return refseq_df

    @staticmethod
    def refseq_chroms_fct(path, refseq_df):
        """[Extract chromosomes from RefSeq GFF]

        Arguments:
            path {[str]} -- [Tmp output file]
            refseq_df {[pd.DataFrame]} -- [RefSeq GFF turned into pandas dataframe]

        Returns:
            [pd.DataFrame] -- [RefSeq chromosomes into pandas dataframe]
        """
        print("### Build temp file (chroms part) / File = {}".format(path))

        if os.path.isfile(path) is False:
            print("# Files don't exist ☒")

            refseq_df_chroms = refseq_df.loc[refseq_df["Region_type"] == "region"]
            index_list = list(refseq_df_chroms.index)

            chroms = [(i, index_list[j + 1] - 1) for j, i in enumerate(index_list) if j < (len(index_list) - 1)]
            refseq_df_chroms = refseq_df_chroms.loc[
                (refseq_df_chroms["NC"].str.contains("NC")) & (refseq_df_chroms["RefSeq_validation"] == "RefSeq")
            ]
            refseq_df_chroms.to_parquet(path)
        else:
            print("# Files exist ✓, Loading ... ")
            refseq_df_chroms = pd.read_parquet(path)
        return refseq_df_chroms

    @staticmethod
    def refseq_pc_genes_fct(path, refseq_df):
        """[Extract protein coding genes from RefSeq GFF]

        Arguments:
            path {[str]} -- [Tmp output file]
            refseq_df {[pd.DataFrame]} -- [RefSeq GFF turned into pandas dataframe]

        Returns:
            [pd.DataFrame] -- [RefSeq protein coding genes into pandas dataframe]
        """
        print("### Build temp file (protein coding genes part) / File = {}".format(path))

        if os.path.isfile(path) is False:
            print("# Files don't exist ☒")

            refseq_df_pc_genes = refseq_df.loc[
                (refseq_df["Attributes"].str.contains("gene_biotype=protein_coding")) & (refseq_df["NC"].str.contains("NC_"))
            ]
            refseq_df_pc_genes.to_parquet(path)
        else:
            print("# Files exist ✓, Loading ... ")
            refseq_df_pc_genes = pd.read_parquet(path)
        return refseq_df_pc_genes

    @staticmethod
    def refseq_mrnas_fct(path, refseq_df):
        """[Extract mRNAs from RefSeq GFF]

        Arguments:
            path {[str]} -- [Tmp output file]
            refseq_df {[pd.DataFrame]} -- [RefSeq GFF turned into pandas dataframe]

        Returns:
            [pd.DataFrame] -- [RefSeq mRNAs into pandas dataframe]
        """
        print("### Build temp file (mRNAs part) / File = {}".format(path))

        if os.path.isfile(path) is False:
            print("# Files don't exist ☒")

            refseq_df_mrna = refseq_df.loc[
                (refseq_df["Attributes"].str.contains("NM_")) & (refseq_df["Region_type"] == "mRNA") & (refseq_df["NC"].str.contains("NC_"))
            ]
            refseq_df_mrna.to_parquet(path)
        else:
            print("# Files exist ✓, Loading ... ")
            refseq_df_mrna = pd.read_parquet(path)
        return refseq_df_mrna

    @staticmethod
    def refseq_exons_fct(path, refseq_df):
        """[Extract exons including UTRs from RefSeq GFF]

        Arguments:
            path {[str]} -- [Tmp output file]
            refseq_df {[pd.DataFrame]} -- [RefSeq GFF turned into pandas dataframe]

        Returns:
            [pd.DataFrame] -- [RefSeq exons including UTRs into pandas dataframe]
        """
        print("### Build temp file (Exons (with UTRs) part) / File = {}".format(path))

        if os.path.isfile(path) is False:
            print("# Files don't exist ☒")

            refseq_df_exons = refseq_df.loc[
                (refseq_df["Attributes"].str.contains("exon-NM"))
                & (refseq_df["Region_type"] == "exon")
                & (refseq_df["NC"].str.contains("NC_"))
            ]
            refseq_df_exons.to_parquet(path)
        else:
            print("# Files exist ✓, Loading ... ")
            refseq_df_exons = pd.read_parquet(path)
        return refseq_df_exons

    @staticmethod
    def refseq_cds_fct(path, refseq_df):
        """[Extract coding exons (CDS) from RefSeq GFF]

        Arguments:
            path {[str]} -- [Tmp output file]
            refseq_df {[pd.DataFrame]} -- [RefSeq GFF turned into pandas dataframe]

        Returns:
            [pd.DataFrame] -- [RefSeq coding exons (CDS) into pandas dataframe]
        """
        print("### Build temp file (coding exons part) / File = {}".format(path))

        if os.path.isfile(path) is False:
            print("# Files don't exist ☒")

            refseq_df_cds = refseq_df.loc[
                (refseq_df["Attributes"].str.contains("NP_")) & (refseq_df["Region_type"] == "CDS") & (refseq_df["NC"].str.contains("NC_"))
            ]
            refseq_df_cds.to_parquet(path)
        else:
            print("# Files exist ✓, Loading ... ")
            refseq_df_cds = pd.read_parquet(path)
        return refseq_df_cds

    @staticmethod
    def mp_build_list(t, l, chroms, mrnas, pc_genes, exons, type_exon):
        """[MP function / Extract Exons, CDS & corresponding mRNAs from protein coding genes]

        Arguments:
            t {[tuple]} -- [Enumerate index [0], dataframe index [1]]
            l {[list]} -- [Empty list to fill in order to build the final dataframe]
        """
        # ENUMERATE INDEX & GENE INDEX
        j, index = t[0], t[1]

        if j < len(list(pc_genes.index)) - 1:

            # GET START & STOP INDEX
            start_index, stop_index = index, list(pc_genes.index)[j + 1]

            # GET GENE
            current_gene = [e.replace("Name=", "") for e in pc_genes.loc[start_index]["Attributes"].split(";") if "Name" in e][0]
            strand = pc_genes.loc[start_index]["Strand"]

            # ITERATE TO FIND CORRESPONDING EXONS
            for e, cds in enumerate(list(exons.index)):
                if cds > start_index and cds < stop_index:
                    # print(e, cds, start_index, stop_index)
                    # print(exons.loc[cds])
                    current_gene_cds = [e.replace("gene=", "") for e in exons.loc[cds]["Attributes"].split(";") if "gene" in e][0]
                    parent_mrna = [e.replace("Parent=rna-", "") for e in exons.loc[cds]["Attributes"].split(";") if "Parent" in e][0].split(
                        "."
                    )[0]
                    mrna_index = mrnas.loc[mrnas["Attributes"].str.contains(parent_mrna)].index[0]

                    # COMPARE GENE NAMES
                    if current_gene == current_gene_cds:
                        l.append(
                            {
                                "Gene": current_gene,
                                "Gene_start": pc_genes.loc[start_index]["Start"],
                                "Gene_stop": pc_genes.loc[start_index]["End"],
                                "Exon_start": exons.loc[cds]["Start"],
                                "Exon_stop": exons.loc[cds]["End"],
                                "mRNA": parent_mrna,
                                "mRNA_start": mrnas.loc[mrna_index]["Start"],
                                "mRNA_stop": mrnas.loc[mrna_index]["End"],
                                "Exon_type": type_exon,
                                "Strand": strand,
                            }
                        )

    def select_exons(self, path, chroms, mrnas, pc_genes, exons, type_exon):
        """[Main function to launch MP]

        Arguments:
            path {[str]} -- [description]
            chroms {[pd.DataFrame]} -- [description]
            mrnas {[pd.DataFrame]} -- [description]
            pc_genes {[pd.DataFrame]} -- [description]
            exons {[pd.DataFrame]} -- [description]
            type_exon {[str]} -- [description]

        Returns:
            [pd.DataFrame] -- [Processed dataframe with filtered exons]
        """
        print("### Retrieve exons & coding exons corresponding to genes / File = {}".format(path))

        if os.path.isfile(path) is False:
            print("# Files don't exist ☒")

            # CREATE LIST
            l = list()

            # MP FCT
            m = multiprocessing.Manager()
            l = m.list()
            parmap.starmap(
                self.mp_build_list, list(zip(enumerate(list(pc_genes.index)))), l, chroms, mrnas, pc_genes, exons, type_exon, pm_pbar=True
            )

            df = pd.DataFrame(list(l)).sort_values(by=["Gene", "Exon_start", "Exon_stop"]).drop_duplicates()
            df["Length"] = df["Exon_stop"] - df["Exon_start"]
            df = df.sort_values(by=["Gene", "Exon_start", "Exon_stop", "Length"], ascending=[True, True, True, False])
            # df = df.drop_duplicates(subset=['Gene', 'Exon_start'], keep='last').drop_duplicates(subset=['Gene', 'Exon_stop'], keep='first')
            df.to_parquet(path)
        else:
            print("# Files exist ✓, Loading ... ")
            return pd.read_parquet(path)

    @staticmethod
    def concat_exons_cds_fct(path_tmp_exons_cds, path_tmp_exons, path_tmp_cds):
        """[Concat exons[1] & CDS[2] tmp files]

        Arguments:
            path_tmp_exons_cds {[str]} -- [Output concatenated file]
            path_tmp_exons {[str]} -- [Exons tmp file]
            path_tmp_cds {[str]} -- [Coding exons (CDS) tmp file]

        Returns:
            [pd.DataFrame] -- [Concatenated dataframe]
        """
        if os.path.isfile(path_tmp_exons_cds) is False:
            concat_df_exons_cds = (
                pd.concat([pd.read_parquet(path_tmp_exons), pd.read_parquet(path_tmp_cds)])
                .sort_values(by=["Gene", "Exon_start", "Exon_stop"])
                .reset_index(drop=True)
            )
            concat_df_exons_cds["ranges"] = (
                concat_df_exons_cds["Exon_start"].astype(str) + "-" + concat_df_exons_cds["Exon_stop"].astype(str)
            )
            concat_df_exons_cds.to_parquet(path_tmp_exons_cds)

        else:
            print("# Files exist ✓, Loading ... ")
            concat_df_exons_cds = pd.read_parquet(path_tmp_exons_cds)
        return concat_df_exons_cds

    @staticmethod
    def groupby_mrnas_fct(path, concat_exons_cds):
        """[summary]

        Arguments:
            path {[str]} -- [Output path]
            concat_exons_cds {[pd.DataFrame]} -- [Tmp concatenated file]

        Returns:
            [pd.DataFrame] -- [Output final processed RefSeq file (one line / exon)]
        """
        print("### Groupby mRNAs / File = {}".format(path))

        if os.path.isfile(path) is False:
            print("# Files don't exist ☒")

            # Merge Exon info & mRNA list / exon
            df = (
                pd.merge(
                    concat_exons_cds[["Gene", "ranges", "Exon_type", "mRNA"]]
                    .groupby(["Gene", "Exon_type", "ranges"])["mRNA"]
                    .apply(list)
                    .reset_index(),
                    concat_exons_cds[["Gene", "ranges", "Exon_start", "Exon_stop", "Length", "Exon_type", "Strand"]].drop_duplicates(),
                    on=["Gene", "ranges", "Exon_type"],
                )
                .sort_values(by=["Gene", "ranges"])
                .reset_index(drop=True)
            )

            # Compute columns
            df["mRNA_exons_nb"] = df["mRNA"].apply(len)
            df["Strand"] = df["Strand"].replace({"-": 0, "+": 1})

            # DF with mRNA by gene
            df_gene_mrna = df[["Gene", "mRNA"]].groupby("Gene")["mRNA"].apply(list).reset_index()
            df_gene_mrna["mRNA_gene"] = df_gene_mrna["mRNA"].apply(lambda r: list(set([sub_e for e in r for sub_e in e])))
            df_gene_mrna = df_gene_mrna.drop(["mRNA"], axis=1)
            df_gene_mrna["mRNA_gene_nb"] = df_gene_mrna["mRNA_gene"].apply(len)

            # Merge previous & compute columns
            df = pd.merge(df, df_gene_mrna, on="Gene")
            df = df.rename({"mRNA": "mRNA_exons"}, axis=1)
            df["Ratio"] = df["mRNA_exons_nb"].astype(str) + "/" + df["mRNA_gene_nb"].astype(str)
            df["Ratio_num"] = df["Ratio"].apply(eval)
            df.loc[df["Ratio_num"] < 1, "Const_Alt"] = "Alt"
            df.loc[df["Ratio_num"] == 1, "Const_Alt"] = "Const"
            df.to_parquet(path)
        else:
            print("# Files exist ✓, Loading ... ")
            df = pd.read_parquet(path)
        return df

    @staticmethod
    def process_variable_exons(gene, complete_list, df):
        """[Retrieve overlapping exons to detect modifiable/variable coding regions]

        Arguments:
            gene {[str]} -- [Gene symbol]
            complete_list {[list]} -- [Empty list to fill]
            df {[pd.DataFrame]} -- [Raw df to process]

        Returns:
            [list] -- [Filled list]
        """

        for gene in tqdm(df.Gene.unique().tolist()):
            l = list()
            l_start = list()
            l_stop = list()
            df_gene = df.loc[df["Gene"] == gene]
            df_gene[["Exon_start", "Exon_stop"]] = df_gene.ranges.str.split("-", expand=True)
            df_gene["Share"] = False

            if len(df_gene.Exon_start.unique()) != df_gene.shape[0]:
                previous_start, previous_stop, previous_mrnas = 0, 0, []
                for j, r in df_gene.iterrows():
                    start, stop, mrnas = r["Exon_start"], r["Exon_stop"], r["mRNA_exons"]
                    if previous_start == start:

                        l_start.append(start)

                        first_exon = "{}-{}".format(start, previous_stop)
                        first_mrnas = list(previous_mrnas) + list(mrnas)
                        second_exon = "{}-{}".format(previous_stop, stop)
                        second_mrnas = mrnas
                        l.append(
                            {
                                "Gene": r["Gene"],
                                "ranges": first_exon,
                                "mRNA_exons": first_mrnas,
                                "mRNA_exons_nb": len(first_mrnas),
                                "mRNA_gene": r["mRNA_gene"],
                                "mRNA_gene_nb": r["mRNA_gene_nb"],
                                "Ratio": "{}/{}".format(len(first_mrnas), r["mRNA_gene_nb"]),
                                "Exon_start": start,
                                "Exon_stop": previous_stop,
                                "Share": True,
                                "Strand": r["Strand"],
                                "Exon_type": "CDS",
                            }
                        )
                        l.append(
                            {
                                "Gene": r["Gene"],
                                "ranges": second_exon,
                                "mRNA_exons": second_mrnas,
                                "mRNA_exons_nb": len(second_mrnas),
                                "mRNA_gene": r["mRNA_gene"],
                                "mRNA_gene_nb": r["mRNA_gene_nb"],
                                "Ratio": "{}/{}".format(len(second_mrnas), r["mRNA_gene_nb"]),
                                "Exon_start": previous_stop,
                                "Exon_stop": stop,
                                "Share": True,
                                "Strand": r["Strand"],
                                "Exon_type": "CDS",
                            }
                        )

                    if previous_stop == stop:

                        l_stop.append(stop)

                        first_exon = "{}-{}".format(previous_start, start)
                        first_mrnas = list(mrnas)
                        second_exon = "{}-{}".format(start, stop)
                        second_mrnas = list(previous_mrnas) + list(mrnas)
                        l.append(
                            {
                                "Gene": r["Gene"],
                                "ranges": first_exon,
                                "mRNA_exons": first_mrnas,
                                "mRNA_exons_nb": len(first_mrnas),
                                "mRNA_gene": r["mRNA_gene"],
                                "mRNA_gene_nb": r["mRNA_gene_nb"],
                                "Ratio": "{}/{}".format(len(first_mrnas), r["mRNA_gene_nb"]),
                                "Exon_start": previous_start,
                                "Exon_stop": start,
                                "Share": True,
                                "Strand": r["Strand"],
                                "Exon_type": "CDS",
                            }
                        )
                        l.append(
                            {
                                "Gene": r["Gene"],
                                "ranges": second_exon,
                                "mRNA_exons": second_mrnas,
                                "mRNA_exons_nb": len(second_mrnas),
                                "mRNA_gene": r["mRNA_gene"],
                                "mRNA_gene_nb": r["mRNA_gene_nb"],
                                "Ratio": "{}/{}".format(len(second_mrnas), r["mRNA_gene_nb"]),
                                "Exon_start": start,
                                "Exon_stop": stop,
                                "Share": True,
                                "Strand": r["Strand"],
                                "Exon_type": "CDS",
                            }
                        )

                    previous_start = start
                    previous_stop = stop
                    previous_mrnas = mrnas

            complete_list.append(
                pd.concat(
                    [df_gene.loc[(~df_gene["Exon_start"].isin(l_start)) & (~df_gene["Exon_stop"].isin(l_stop))], pd.DataFrame(l)], axis=0
                ).sort_values(by="ranges")
            )
        return complete_list

    def get_shared_exons(self, path, df):
        print("### Process variable regions / File = {}".format(path))

        if os.path.isfile(path) is False:
            print("# Files don't exist ☒")

            complete_list = list()
            df = df.loc[df["Exon_type"] == "CDS"]
            # parmap.starmap(self.mp_process_variable, list(zip(df.Gene.unique().tolist())), complete_list, df, pm_pbar=True)
            complete_list = self.process_variable_exons(df.Gene.unique().tolist(), complete_list, df)
            df_with_variable = pd.concat(list(complete_list)).reset_index(drop=True)
            df_with_variable["Length"] = df_with_variable["Exon_stop"].astype(int) - df_with_variable["Exon_start"].astype(int)
            df_with_variable["Ratio_num"] = df_with_variable["Ratio"].apply(eval)
            df_with_variable.loc[df_with_variable["Ratio_num"] < 1, "Const_Alt"] = "Alt"
            df_with_variable.loc[df_with_variable["Ratio_num"] == 1, "Const_Alt"] = "Const"
            df_with_variable = pd.merge(
                df_with_variable,
                df_with_variable.groupby("Gene")["ranges"].count().reset_index().rename({"ranges": "CDS_count"}, axis=1),
                on="Gene",
            )
            df_with_variable.to_parquet(path)
        else:
            print("# Files exist ✓, Loading ... ")
            df_with_variable = pd.read_parquet(path)
        return df_with_variable

    def counter(self, r, concat_cds_exon_join):

        three_prime_counter = 0
        five_prime_counter = 0

        counter_dict = dict(collections.Counter(r))
        for k, v in counter_dict.items():
            if v > 1:
                exons_including_utrs = concat_cds_exon_join.loc[concat_cds_exon_join["CDS_ranges"] == k, "Exons_ranges"].tolist()
                exons_including_utrs_split = list(sorted(set([sub_e for e in exons_including_utrs for sub_e in e.split(",")])))
                strand_l = list(set(concat_cds_exon_join.loc[concat_cds_exon_join["CDS_ranges"] == k, "Strand"].values))

                gene = list(set(concat_cds_exon_join.loc[concat_cds_exon_join["CDS_ranges"] == k, "Gene"].values))

                if len(strand_l) > 1:
                    strand = 1
                elif len(strand_l) == 0:
                    strand = 1
                else:
                    strand = strand_l[0]

                tmp_l = [
                    k
                    for k, v in dict(collections.Counter([sub_e for elem in exons_including_utrs for sub_e in elem.split(",")])).items()
                    if v < len(exons_including_utrs)
                ]

                for utr in tmp_l:
                    index = exons_including_utrs_split.index(utr)
                    exons_nb = len(exons_including_utrs_split)
                    if index > round(exons_nb / 2, 0):
                        check = True
                    else:
                        check = False

                    if check is True and int(strand) == -1:
                        check = False
                    elif check is False and int(strand) == -1:
                        check = True

                    t = "3_prime_UTR" if check is True else "5_prime_UTR"
                    if t == "3_prime_UTR":
                        three_prime_counter += 1
                    elif t == "5_prime_UTR":
                        five_prime_counter += 1
        return pd.Series([five_prime_counter, three_prime_counter, len(list(counter_dict.values()))])

    def compute_utrs(self, path, complete_df):
        print("### Process UTRs / File = {}".format(path))

        if os.path.isfile(path) is False:
            print("# Files don't exist ☒")

            df = pd.DataFrame(complete_df.groupby("Gene")["ranges"].apply(lambda r: list(r))).reset_index()
            mrna_count = complete_df.loc[complete_df["mRNA_gene_nb"] > 1, ["Gene", "mRNA_gene_nb"]].drop_duplicates()
            print(mrna_count)

            tmp_df = complete_df.loc[
                (complete_df["Gene"].isin(mrna_count.Gene.values.tolist())) & (complete_df["Exon_type"] == "CDS")
            ].explode("mRNA_exons")
            cds_join = pd.DataFrame(tmp_df.groupby(["Gene", "mRNA_exons"])["ranges"].apply(",".join))
            cds_join.columns = ["CDS_ranges"]
            print(cds_join)

            tmp_df = complete_df.loc[
                (complete_df["Gene"].isin(mrna_count.Gene.values.tolist())) & (complete_df["Exon_type"] == "Exon")
            ].explode("mRNA_exons")
            exon_join = pd.DataFrame(tmp_df.groupby(["Gene", "mRNA_exons"])["ranges"].apply(",".join))
            exon_join.columns = ["Exons_ranges"]

            concat_cds_exon_join = pd.concat([exon_join, cds_join], axis=1)
            concat_cds_exon_join = pd.merge(concat_cds_exon_join, complete_df[["Gene", "Strand"]].drop_duplicates(), on="Gene")

            print(concat_cds_exon_join)

            df = pd.DataFrame(cds_join.groupby("Gene")["CDS_ranges"].apply(lambda r: list(r))).reset_index()
            print(df)
            df[["5_prime_modif", "3_prime_modif", "Nb_combi"]] = df.CDS_ranges.progress_apply(
                lambda r: self.counter(r, concat_cds_exon_join)
            )
            print(df)
            df.to_parquet(path)
        else:
            print("# Files exist ✓, Loading ... ")
            df = pd.read_parquet(path)
        return df


if __name__ == "__main__":

    c = ProcessRefSeq(yaml["1_GENOMICS"]["Final"]["refseq_cds_with_variable"])
    # print(c.groupby_mrnas)
