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
from pandarallel import pandarallel

pandarallel.initialize(nb_workers=20, progress_bar=True)

## YAML FILES CONFIG
yaml = load_config_file(config_file="/home/weber/PycharmProjects/EXOTIC/clean/src/config_clean_clean.yaml")


class CorrectExpression:
    def __init__(self, path):

        if os.path.isfile(path) is False:
            # * 1 Correction by expression
            gtex_corrected = self.launch(yaml["2_EXPRESSION"]["Final"]["transcript_check_refseq"])

            # * 2 Load RefSeq & remove genes with transcript identical on CDS
            refseq = pd.read_parquet("/gstock/EXOTIC/data/GENOMICS/refseq_37_processed_cds_variable.parquet")
            utrs = pd.read_parquet("/gstock/EXOTIC/data/GENOMICS/refseq_37_processed_utrs_analysis.parquet")
            refseq = refseq.loc[~refseq["Gene"].isin(utrs.loc[utrs["Nb_combi"] == 1])]

            # * 3 Remove unexpressed transcripts & recompute frequencies & usage
            refseq_correct_mrnas = self.correct_refseq_by_expression(
                yaml["2_EXPRESSION"]["Final"]["refseq_corrected_cds_with_variable"], refseq, gtex_corrected
            )
            self.final_corrected_refseq = self.process_corrected_file(path, refseq_correct_mrnas)
        else:
            self.final_corrected_refseq = pd.read_parquet(path)
            print(self.final_corrected_refseq)

    def launch(self, path):

        print("### Launch pipeline correction method / File = {}".format(path))

        if os.path.isfile(path) is False:
            print("# Files don't exist ☒")
            # MAIN FCT
            merge_check = self.merge_check_files(
                path_gtex_reads=yaml["2_EXPRESSION"]["External"]["transcript_reads"],
                path_gtex_tpm=yaml["2_EXPRESSION"]["External"]["transcript_tpm"],
                biomart_path=yaml["1_GENOMICS"]["External"]["biomart"],
            )

            # ADD COLUMNS
            merge_check["transcript_id"] = merge_check["transcript_id"].apply(lambda r: r.split(".")[0])
            merge_check = pd.merge(biomart_grch37, merge_check, on="transcript_id")

            # RETRIEVED TRUE & IN REFSEQ
            merge_check = merge_check.loc[(merge_check["check"] == True) & (merge_check["RefSeq_mRNAID"].isna() == False)]

            # OUTPUT
            merge_check.to_csv(path, compression="gzip", sep="\t")
        else:
            print("# Files exist ✓, Loading ... ")

            merge_check = pd.read_csv(path, compression="gzip", sep="\t")
        return merge_check

    @staticmethod
    def load_biomart(path):
        print("### Load BIOMART file / File = {}".format(path))

        biomart_grch37 = pd.read_csv(path, sep="\t", compression="gzip")

        biomart_grch37 = biomart_grch37[
            [
                "Transcript stable ID",
                "Transcript type",
                "APPRIS annotation",
                "RefSeq mRNA ID",
                "Transcript support level (TSL)",
            ]
        ]

        biomart_grch37.columns = [
            "transcript_id",
            "transcript_biotype",
            "APPRIS",
            "RefSeq_mRNAID",
            "TSL",
        ]

        return biomart_grch37

    def check_file(self, filepath, biomart_path):

        print("--- Checking GTEx raw files ---")

        # CHECK IF PROCESSING EITHER TPM OR READS
        check_type = "tpm" if "tpm" in filepath else "reads"
        print(check_type)

        o_file = filepath.replace(".txt.gz", "_checked_complete.parquet")
        o_file_lite = filepath.replace(".txt.gz", "_checked_lite.txt.gz")
        if os.path.isfile(o_file) is False:
            print("# Files don't exist ☒")

            # LOAD FILE
            df = pd.read_csv(
                filepath,
                compression="gzip",
                sep="\t",
                # nrows=100,
            )

            df["gene_id"] = df["gene_id"].apply(lambda r: r.split(".")[0])

            # MERGE WITH BIOMART TO HAVE GENE SYMBOLS
            biomart_grch37 = self.load_biomart(biomart_path)
            df = pd.merge(biomart_grch37, df, on="gene_id")

            # CUTOFFS
            cutoff_reads = 6
            cutoff_tpm = 0.1

            # CHECK FUNCTION
            @staticmethod
            def check(r, check_type):
                i = 0
                if check_type == "reads":
                    cutoff = cutoff_reads
                elif check_type == "tpm":
                    cutoff = cutoff_tpm
                for e in r[7:]:
                    if e >= cutoff:
                        i += 1
                return i

            # PARALLEL APPLY TO CHECK CUTOFF
            df["check_{}_below_cutoff".format(check_type)] = df.parallel_apply(lambda r: check(r, check_type), axis=1)
            df["check_{}".format(check_type)] = df["check_{}_below_cutoff".format(check_type)].parallel_apply(
                lambda r: True if r / 11688 >= 0.2 else False
            )
            print(df)

            # OUTPUT COMPLETE & LITE
            df.to_parquet(o_file, index=False)
            df = df[list(df.columns)[:7] + list(df.columns)[-2:]]
            df.to_csv(o_file_lite, index=False, compression="gzip", sep="\t")
        else:

            print("# Files exist ✓")

            df = pd.read_csv(o_file_lite, compression="gzip", sep="\t")

        return df

    def merge_check_files(self, path_gtex_reads, path_gtex_tpm, biomart_path):
        if os.path.isfile(yaml["2_EXPRESSION"]["Final"]["transcript_check"]) is False:

            # CHECK GTEx FILES TO RETRIEVE TRANSCRIPTS PASSING CUTOFFS
            check_reads_file = self.check_file(path_gtex_reads)
            check_tpm_file = self.check_file(path_gtex_tpm)

            # MERGE FILES
            merge_check = pd.merge(
                check_reads_file,
                check_tpm_file,
                on=[
                    "Chrom",
                    "Start",
                    "End",
                    "Strand",
                    "gene_id",
                    "gene_symbol",
                    "transcript_id",
                ],
            )

            # COMPUTE COLUMN TO RETRIEVE TRANSCRIPTS PASSING BOTH CONDITIONS
            merge_check["check"] = merge_check.apply(
                lambda r: True if set([r["check_reads"], r["check_tpm"]]) == {True} else False,
                axis=1,
            )

            # OUTPUT
            merge_check.to_csv(
                yaml["2_EXPRESSION"]["Final"]["transcript_check"],
                compression="gzip",
                sep="\t",
                index=False,
            )
        else:
            merge_check = pd.read_csv(
                yaml["2_EXPRESSION"]["Final"]["transcript_check"],
                compression="gzip",
                sep="\t",
            )
        return merge_check

    @staticmethod
    def correct_refseq_by_expression(path, raw_df, transcripts_checked):
        print("### Remove unexpressed mRNAs / File = {}".format(path))
        if os.path.isfile(path) is False:
            print("# Files don't exist ☒")

            # EXPLODE FROM LISTS TO ROWS
            df = raw_df.explode("mRNA_exons")
            df["mRNA_exons"] = df["mRNA_exons"].apply(lambda r: r.split(".")[0])

            # FILTERING mRNAs which respect GTEx expression cutoffs (reads & TPM)
            df = df.loc[df["mRNA_exons"].isin(transcripts_checked.RefSeq_mRNAID.unique().tolist())]

            # GROUPBY
            df = df.groupby(["Gene", "ranges"])["mRNA_exons"].apply(list)

            df = df.reset_index()
            # df = df.rename({"mRNA_exons": "mRNA_exons"}, axis=1)
            df["mRNA_nb"] = df["mRNA_exons"].apply(len)

            # MERGE WITH PREVIOUS INFO
            df = pd.merge(df, raw_df[["Gene", "ranges", "Share", "Strand"]], on=["Gene", "ranges"])

            # TMP FILE
            df.to_parquet(path)

        else:
            print("# Files exist ✓, Loading ... ")

            df = pd.read_parquet(path)
        return df

    @staticmethod
    def process_corrected_file(path, df):
        print("### Correct RefSeq frequency & exons / File = {}".format(path))

        if os.path.isfile(path) is False:
            print("# Files don't exist ☒")

            # COMPUTE mRNAs LIST AT GENE LEVEL BASED ON GROUPBY
            # df["mRNA_exons"] = df["mRNA_exons"].apply(eval)
            df_gene_level = df.groupby("Gene")["mRNA_exons"].apply(list).reset_index()

            # EXTEND LIST OF LIST INTO ONE LIST TRANSFORMED TO A SET
            df_gene_level["mRNA_gene"] = df_gene_level["mRNA_exons"].apply(lambda r: list(set([sub_e for e in r for sub_e in e])))

            # MERGE EXON & GENE LEVELS
            df_corrected = pd.merge(df, df_gene_level[["Gene", "mRNA_gene"]], on="Gene")

            # COMPUTE
            df_corrected["mRNA_exons"] = df_corrected["mRNA_exons"].apply(lambda r: list(set(r)))
            df_corrected["mRNA_nb"] = df_corrected["mRNA_exons"].str.len()
            df_corrected["mRNA_nb_total"] = df_corrected["mRNA_gene"].apply(len)
            df_corrected["Ratio"] = df_corrected["mRNA_nb"].astype(str) + "/" + df_corrected["mRNA_nb_total"].astype(str)
            df_corrected["Ratio_num"] = df_corrected["Ratio"].apply(eval)
            df_corrected["Ratio_num"] = df_corrected["Ratio_num"].astype(float)

            # CONST & ALT
            df_corrected.loc[df_corrected["Ratio_num"] < 1, "Const_Alt"] = "Alt"
            df_corrected.loc[df_corrected["Ratio_num"] == 1, "Const_Alt"] = "Const"

            # ALT BINS
            bins = [0, 0.2, 0.4, 0.6, 0.8, 1]
            labels = bins.copy()
            labels_ratio = [str(round(labels[j], 1)) + " - " + str(round(labels[j + 1], 1)) for j in range(len(labels) - 1)]
            df_corrected["Ratio_num_bins"] = pd.cut(df_corrected["Ratio_num"], bins=bins, labels=labels_ratio, include_lowest=True)
            df_corrected.loc[
                (df_corrected["Const_Alt"] == "Alt") & (df_corrected["Ratio_num_bins"] == "0.8 - 1"), "Ratio_num_bins_update"
            ] = "0.8 - 0.99"
            df_corrected.loc[
                (df_corrected["Const_Alt"] == "Alt") & (df_corrected["Ratio_num_bins"] != "0.8 - 1"), "Ratio_num_bins_update"
            ] = df_corrected.loc[(df_corrected["Const_Alt"] == "Alt") & (df_corrected["Ratio_num_bins"] != "0.8 - 1")]["Ratio_num_bins"]
            df_corrected.loc[
                (df_corrected["Const_Alt"] == "Const") & (df_corrected["Ratio_num_bins"] == "0.8 - 1"), "Ratio_num_bins_update"
            ] = "1.00"
            df_corrected["MAP"] = df_corrected["Gene"] + "_" + df_corrected["ranges"]

            # DUPLICATES
            df_corrected = df_corrected.drop_duplicates(subset=["Gene", "ranges"])

            # START, STOP, LENGTH
            df_corrected["Start"] = df_corrected["ranges"].apply(lambda r: r.split("-")[0])
            df_corrected["Start"] = df_corrected["Start"].astype(int)
            df_corrected["End"] = df_corrected["ranges"].apply(lambda r: r.split("-")[1])
            df_corrected["End"] = df_corrected["End"].astype(int)
            df_corrected["Length"] = df_corrected["End"] - df_corrected["Start"]

            # CDS count
            t = df_corrected.groupby(["Gene", "mRNA_nb_total"])["ranges"].agg("nunique").reset_index()

            # MERGE WITH PREVIOUS
            df_corrected = pd.merge(df_corrected, t.rename({"ranges": "CDS_count"}, axis=1), on=["Gene", "mRNA_nb_total"])

            # OUPUTS
            df_corrected.to_parquet(path, index=False)

        else:
            print("# Files exist ✓, Loading ... ")
            df_corrected = pd.read_parquet(path)
        return df_corrected


if __name__ == "__main__":
    c = CorrectExpression(yaml["2_EXPRESSION"]["Final"]["refseq_corrected_cds_recomputed"])
    df = c.final_corrected_refseq
    print(df.loc[df["Share"] == True])
