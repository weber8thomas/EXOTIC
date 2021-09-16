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

import parmap
import re
import requests
import subprocess
import sys
import time
import warnings
import math
from pandarallel import pandarallel

pandarallel.initialize(nb_workers=20, progress_bar=True)

pd.options.mode.chained_assignment = None  # default='warn'
sys.path.append("/home/weber/PycharmProjects/EXOTIC")
from utils.utils import load_config_file, mkdir, convert_bins_into_labels


# import hail as hl

# hl.init(min_block_size=128)
# os.environ["PYSPARK_SUBMIT_ARGS"] = "--driver-memory 200 pyspark-shell"


tqdm.pandas()

## YAML FILES CONFIG
# yaml = load_config_file(config_file="clean/src/config_clean_clean.yaml")
yaml = load_config_file(config_file="/home/weber/PycharmProjects/gene_isoforms/src/config/config_files.yaml")

## JSON DICTS CONFIG
dicts = json.load(open("clean/src/config/EXOTIC_config.json"))


class Pathogenicity:
    def __init__(self):
        # self.turn_refseq_miso_exons_into_bed(
        #     "/gstock/EXOTIC/data/GENOMICS/refseq_37_processed_cds_variable.parquet", "/gstock/EXOTIC/data/GENOMICS/refseq_complete.bed"
        # )

        # refseq = self.exon_bins(
        #     yaml["1_GENOMICS"]["Final"]["refseq_cds_with_variable"], "/gstock/EXOTIC/data/GENOMICS/refeq_exons_bins.parquet"
        # )

        # # gnomad_complete = self.read_gnomad_ht_file(
        # #     "/gstock/EXOTIC/data/GENOMICS/refseq_complete.bed",
        # #     yaml["3_EXONS_PROPERTIES"]["External"]["gnomad_2_1_1"],
        # #     "/gstock/EXOTIC/data/VARIATIONS/gnomad_full_expressed_genes.tsv.gz",
        # # )
        # # print(gnomad_complete)

        # # clinvar_omim = self.load_clinvar(clinvar_path=yaml["3_EXONS_PROPERTIES"]["Final"]["refseq_clinvar_hail_retrieved"])
        # clinvar_data = self.read_and_process_clinvar_vcf(
        #     clinvar_vcf_path=yaml["3_EXONS_PROPERTIES"]["External"]["clinvar_latest"],
        #     output_file=yaml["5_PATHOGENICITY"]["TMP"]["clinvar_full"],
        # )
        # print(clinvar_data)

        # self.distribution_pathogenic_variations(yaml["1_GENOMICS"]["TMP"]["tmp_refseq_concat_exons_cds_filtered"], clinvar_data, refseq)
        # # self.distribution_benign_variations(yaml["1_GENOMICS"]["TMP"]["tmp_refseq_concat_exons_cds_filtered"], gnomad_complete)

        # exit()

        # # * REFSEQ OMIM
        # refseq_omim = self.load_refseq_pc_genes(path=yaml["1_GENOMICS"]["TMP"]["tmp_refseq_pc_genes"])

        # # * DEXT
        # dext_omim = self.load_dext(yaml["4_DEXT"]["Final"]["dext"], refseq_omim, yaml["5_PATHOGENICITY"]["TMP"]["dext_omim"])

        # # * DL OMIM entries
        # # self.download_omim(dext_omim, output_dir=yaml["5_PATHOGENICITY"]["TMP"]["omim_api_dump_directory"])
        clinvar_lite = pd.read_parquet("/gstock/EXOTIC/data/VARIATIONS/clinvar_20210123_lite_table.parquet")
        clinvar_path = clinvar_lite.loc[(clinvar_lite["Status"] == "Pathogenic") & (~clinvar_lite["CLNREVSTAT"].str.contains("onflicting"))]
        clinvar_genes_path = clinvar_path.GENE.unique().tolist()

        genes_path_processed = yaml["1_GENOMICS"]["Final"]["refseq_genes_processed_miso_siso"]
        genes = pd.read_parquet(genes_path_processed)
        genes.loc[genes["Gene"].isin(clinvar_genes_path), "D/H"] = "Disease"
        genes.loc[~genes["Gene"].isin(clinvar_genes_path), "D/H"] = "Healthy"

        genes = genes.loc[genes["D/H"] == "Disease"]

        genes["OMIM"] = genes["Attributes"].apply(lambda r: [e for e in r.split(";") if "MIM" in e])
        genes_omim = genes.loc[genes["OMIM"].str.len() > 0]
        genes_omim["OMIM"] = genes_omim["OMIM"].apply(lambda r: [e for e in r[0].split(",") if "MIM:" in e])
        genes_omim = genes_omim.loc[genes_omim["OMIM"].str.len() > 0]
        genes_omim["OMIM"] = genes_omim["OMIM"].apply(lambda r: r[0].replace("MIM:", ""))

        # * BUILD OMIM DF
        omim = self.build_omim(
            path=yaml["3_DISEASES"]["TMP"]["omim_genes_pheno"],
            omim=genes_omim,
        )
        exit()

        # * MELT & PROCESSED OMIM
        omim, omim_matrix = self.processed_omim(omim)

        # * SEARCH SPECIFIC PHENOS
        omim_specific_pheno = self.search_specific_pheno(omim, yaml["5_PATHOGENICITY"]["TMP"]["omim_specific_pheno"])

        # * MERGE DEXT & OMIM SPECIFIC
        dext_omim_specific = self.merge_dext_omim_specific(dext_omim, omim_specific_pheno)

        # * RETRIEVE GENES through clinvar with variation referenced
        clinvar_omim = self.load_clinvar(clinvar_path=yaml["3_EXONS_PROPERTIES"]["Final"]["refseq_clinvar_hail_retrieved"])

        # * RETRIEVE OMIM PHENO FROM VARIATIONS
        clinvar_omim_processed = self.retrieve_omim_variants_info(
            clinvar_omim, path=yaml["5_PATHOGENICITY"]["TMP"]["clinvar_omim_processed"]
        )

        # ? TMP

        clinvar_omim_processed = clinvar_omim_processed.explode("PHENOS_OMIM")
        clinvar_omim_processed["OMIM"] = clinvar_omim_processed["OMIM"].astype(str)
        clinvar_omim_processed = clinvar_omim_processed.rename({"PHENOS_OMIM": "Pheno_OMIM"}, axis=1)

        print(omim_matrix)
        print(clinvar_omim_processed)
        omim_matrix_sum = omim_matrix[["OMIM", "Pheno_OMIM", "Matrix_vector_sum"]].sort_values(by=["OMIM", "Pheno_OMIM"]).drop_duplicates()
        omim_matrix_count = (
            omim_matrix[["OMIM", "Pheno_OMIM", "Matrix_vector_sum"]]
            .sort_values(by=["OMIM", "Pheno_OMIM"])
            .drop_duplicates()
            .groupby("OMIM")["Pheno_OMIM"]
            .nunique()
            .reset_index()
        )
        omim_matrix_count.columns = ["OMIM", "Count"]
        print(omim_matrix_sum)
        print(omim_matrix_count)
        omim_matrix_merge = pd.merge(omim_matrix_sum, omim_matrix_count, on="OMIM")
        omim_matrix_merge = omim_matrix_merge.loc[omim_matrix_merge["Count"] > 1]
        print(omim_matrix_merge)
        print(omim_matrix_merge.OMIM.nunique())
        print(omim_matrix_merge.Pheno_OMIM.nunique())

        print(clinvar_omim_processed)

        clinvar_matrix = pd.merge(omim_matrix_merge, clinvar_omim_processed, on=["OMIM", "Pheno_OMIM"])
        clinvar_matrix.groupby("OMIM")
        exit()

        pd.options.display.max_rows = 100

        # * REFSEQ
        refseq_const = self.load_refseq_final(path=yaml["1_GENOMICS"]["Final"]["refseq_cds_with_variable"])
        refseq_const = pd.merge(clinvar_omim, refseq_const, on=["Gene", "MAP"])
        refseq_const = refseq_const.loc[refseq_const["Const_Alt"] == "Const"]
        print(refseq_const.loc[refseq_const["Gene"] == "LMNA"][["MAP", "snpId", "OMIM_variant_id"]])
        exit()

        # * MERGE DEXT OMIM & CLINVAR VARIATIONS
        merge_dext_omim_clinvar = self.merge_dext_omim_clinvar(dext_omim_specific, clinvar_omim)
        print(merge_dext_omim_clinvar[["MAP", "OMIM_variant_id", "dext_tissues_up", "gtex_omim_up", "OMIM_BP"]])
        exit()

        # merge_dext_omim_clinvar.loc[merge_dext_omim_clinvar['Gene'].isin()]
        print(refseq_const.loc[refseq_const["Gene"].isin(merge_dext_omim_clinvar.Gene.unique().tolist())])
        genes = sorted(
            list(
                set(
                    refseq_const.loc[refseq_const["Gene"].isin(merge_dext_omim_clinvar.Gene.unique().tolist())].Gene.unique().tolist()
                ).intersection(set(merge_dext_omim_clinvar.Gene.unique().tolist()))
            )
        )

        print(
            refseq_const.loc[refseq_const["Gene"].isin(genes)][["Gene", "MAP", "Share", "snpId", "OMIM_variant_id"]].sort_values(
                by=["MAP", "snpId"]
            )
        )

    @staticmethod
    def exon_bins(refseq_path, output_path):
        def make_windows(r):
            nb_bins = 20

            absolute = r["Exon_start"]
            relative = 0
            window = 1

            k = math.ceil(r["Exon_window_bp"] / nb_bins)

            bin_ranges = list()

            for e in r["ranges"]:
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

        if os.path.isfile(output_path) is True:

            refseq_raw = pd.read_parquet(refseq_path)
            refseq = refseq_raw.groupby("Gene")["ranges"].apply(list).reset_index()
            refseq["ranges"] = refseq["ranges"].apply(lambda r: list(sorted(r)))
            refseq["Exon_window_bp"] = refseq["ranges"].apply(lambda r: sum([int(e.split("-")[1]) - int(e.split("-")[0]) for e in r]))
            refseq["Exon_start"] = refseq.ranges.apply(lambda r: r[0].split("-")[0])
            refseq["Exon_end"] = refseq.ranges.apply(lambda r: r[-1].split("-")[1])
            refseq[["Exon_start", "Exon_end"]] = refseq[["Exon_start", "Exon_end"]].astype(int)
            refseq = refseq.parallel_apply(make_windows, axis=1)
            refseq = pd.merge(refseq_raw[["Gene", "Strand"]], refseq, on="Gene")
            refseq.to_parquet(output_path)
        else:
            refseq = pd.read_parquet(output_path)

        return refseq

    @staticmethod
    def read_and_process_clinvar_vcf(clinvar_vcf_path, output_file):
        print("# Load ClinVar VCF, build ht + pandas & retrieve corresponding variations / File : {}".format(output_file))

        if os.path.isfile(output_file) is False:

            clinvar_review_status = {
                "practice_guideline": 4,
                "reviewed_by_expert_panel": 3,
                "criteria_provided,_multiple_submitters,_no_conflicts": 2,
                "criteria_provided,_conflicting_interpretations": 1,
                "criteria_provided,_single_submitter": 1,
                "no_assertion_for_the_individual_variant": 0,
                "no_assertion_criteria_provided": 0,
                "no_assertion_provided": 0,
            }

            # bed_file = hl.import_bed(bed_path, reference_genome="GRCh37")

            # hl.import_vcf(clinvar_vcf_path, force_bgz=True).write(clinvar_vcf_path.replace(".vcf.gz", ".ht"))
            clinvar_lite = hl.read_matrix_table("/gstock/biolo_datasets/variation/variation_sets/clinvar/vcf_GRCh37/v2/clinvar_20210123.ht")

            # clinvar_lite = clinvar.filter_rows(hl.is_defined(bed_file[clinvar.locus]))
            # print(clinvar_lite.count())

            clinvar_lite = clinvar_lite.filter_rows(clinvar_lite.info.CLNVC == "single_nucleotide_variant")
            # print(clinvar_lite.count())
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

            # clinvar_lite = clinvar_lite.annotate_rows(
            # MAP=bed_file[clinvar_lite.locus].target,
            #     Gene=hl.str(bed_file[clinvar_lite.locus].target).split("_")[0],
            # )

            clinvar_lite_pd_raw = clinvar_lite.make_table().to_pandas()
            clinvar_lite_pd = clinvar_lite_pd_raw.copy()
            clinvar_lite_pd = clinvar_lite_pd.replace(to_replace="None", value=np.nan)

            # print("0", clinvar_lite_pd.shape[0], clinvar_lite_pd.Gene.nunique())
            clinvar_lite_pd = clinvar_lite_pd.dropna(subset=["ALLELEID", "CLNREVSTAT", "CLNSIG", "GENEINFO"], how="any")
            # print("1", clinvar_lite_pd.shape[0], clinvar_lite_pd.Gene.nunique())
            clinvar_lite_pd["CLNREVSTAT"] = clinvar_lite_pd["CLNREVSTAT"].apply(lambda r: ",".join(r))
            clinvar_lite_pd["CLNSIG"] = clinvar_lite_pd["CLNSIG"].apply(lambda r: ",".join(r))
            clinvar_lite_pd["RS_STARS"] = clinvar_lite_pd["CLNREVSTAT"].map(clinvar_review_status)
            clinvar_lite_pd = clinvar_lite_pd.loc[
                (clinvar_lite_pd["CLNSIG"].str.contains("athogenic"))
                & (~clinvar_lite_pd["CLNSIG"].str.contains("Conflicting_interpretations_of_pathogenicity"))
            ]
            # print("2", clinvar_lite_pd.shape[0], clinvar_lite_pd.Gene.nunique())
            clinvar_lite_pd = clinvar_lite_pd.loc[~clinvar_lite_pd["CLNREVSTAT"].str.contains("interpret")]
            # print("3", clinvar_lite_pd.shape[0], clinvar_lite_pd.Gene.nunique())

            # clinvar_lite_pd["snpId"] = (
            #     clinvar_lite_pd["locus.contig"].astype(str)
            #     + "_"
            #     + clinvar_lite_pd["locus.position"].astype(str)
            #     + "_"
            #     + clinvar_lite_pd["alleles"].apply(lambda r: r[0])
            #     + "_"
            #     + clinvar_lite_pd["alleles"].apply(lambda r: r[1])
            # )
            # gnomad = pd.read_parquet("/gstock/EXOTIC/data/EXON_PROPERTIES/refseq_miso_gnomad_variations.parquet")
            # clinvar_lite_pd = clinvar_lite_pd.loc[~clinvar_lite_pd["snpId"].isin(gnomad.id.values.tolist())]

            clinvar_lite_pd.to_parquet(output_file)
        else:
            clinvar_lite_pd = pd.read_parquet(output_file)

        # print(clinvar_lite_pd)
        return clinvar_lite_pd

    @staticmethod
    def turn_refseq_miso_exons_into_bed(refseq_path, bed_path):
        print("# Turn RefSeq into BED / File : {}".format(bed_path))
        if os.path.isfile(bed_path) is False:
            df = pd.read_parquet(refseq_path)
            # df = df.loc[(df["mRNA_nb_total"] > 1) & (df["CDS_count"] > 1)]

            biomart = pd.read_csv(yaml["1_GENOMICS"]["External"]["biomart"], compression="gzip", sep="\t")
            df = pd.merge(
                biomart[["Gene name", "Chromosome/scaffold name"]]
                .drop_duplicates()
                .rename({"Gene name": "Gene", "Chromosome/scaffold name": "CHROM"}, axis=1),
                df,
                on=["Gene"],
            )
            df["CHROM"] = df["CHROM"].astype(str)
            df["MAP"] = df["Gene"] + "_" + df["ranges"]

            df = df.loc[~df["CHROM"].str.contains("HSCHR|HG")]
            df = df[["CHROM", "Exon_start", "Exon_stop", "MAP"]].sort_values(by=["CHROM", "Exon_start"])
            df[["Exon_start", "Exon_stop"]] = df[["Exon_start", "Exon_stop"]].astype(int)
            df["Exon_start"] = df["Exon_start"] - 1
            df["Exon_stop"] = df["Exon_stop"] - 1
            # bed_path = yaml['3_EXONS_PROPERTIES']['TMP']['refseq_bed_gnomad']
            mkdir(os.path.dirname(bed_path))
            df.to_csv(bed_path, header=False, sep="\t", index=False)
        else:
            pass

    @staticmethod
    def read_gnomad_ht_file(bed_path, gnomad_path_ht, output_file):
        print("# Load gnomAD hail table & retrieve corresponding variations / File : {}".format(output_file.replace(".tsv.gz", ".parquet")))

        if os.path.isfile(output_file.replace(".tsv.gz", ".parquet")) is False:
            import hail as hl

            hl.init(min_block_size=128)
            os.environ["PYSPARK_SUBMIT_ARGS"] = "--driver-memory 200 pyspark-shell"

            data = hl.read_table(gnomad_path_ht)

            bed_file = hl.import_bed(bed_path, reference_genome="GRCh37")
            # print(bed_file.head(10).show())

            data_lite = data.filter(hl.is_defined(bed_file[data.locus]))
            data_lite = data_lite.filter(data_lite.variant_type == "snv")
            data_lite = data_lite.filter(data_lite.freq[0]["AC"] > 0)
            # data_lite['AF'] = data_lite.freq[0]['AC']
            data_lite = data_lite.annotate(
                MAP=bed_file[data_lite.locus].target,
                AF=data_lite.freq[0]["AF"],
                id=hl.str(data_lite.vep.id).replace("/", "_"),
                Gene=hl.str(bed_file[data_lite.locus].target).split("_")[0],
            )

            data_lite = data_lite.select(
                data_lite.AF,
                data_lite.rsid,
                data_lite.variant_type,
                data_lite.vep.end,
                data_lite.id,
                data_lite.vep.most_severe_consequence,
                data_lite.MAP,
                data_lite.Gene,
            )
            data_lite.export(output_file)
            data_lite = pd.read_csv(output_file, compression="gzip", sep="\t")
            data_lite.to_parquet(output_file.replace(".tsv.gz", ".parquet"))

        else:
            # data_lite = pd.read_csv(output_file, compression="gzip", sep="\t")

            data_lite = pd.read_parquet(output_file.replace(".tsv.gz", ".parquet"))
            # data_lite.to_parquet(output_file.replace(".tsv.gz", ".parquet"))
        # print(data_lite)
        return data_lite

        # data_lite = data_lite.key_by('id')

    # @staticmethod
    def distribution_benign_variations(self, refseq_path, df):

        df = df.loc[df["most_severe_consequence"] == "missense_variant"]

        # df = df.head(1000)
        df["end"] = df["end"].astype(int)

        refseq_with_genes_coord = (
            pd.read_parquet(refseq_path)[["Gene", "Gene_start", "Gene_stop", "Strand"]].drop_duplicates().sort_values(by="Gene")
        )
        refseq_with_genes_coord["Strand"] = refseq_with_genes_coord["Strand"].replace({"-": 0, "+": 1})

        refseq_with_genes_coord = self.compute_bins(refseq_with_genes_coord, nb_bins=[20])

        df = pd.merge(refseq_with_genes_coord, df, on=["Gene"])

        nb_bin = 20
        df["snpId_bin_{}".format(nb_bin)] = df.parallel_apply(lambda r: self.compute_snv_bin_fct(r, nb_bin, col="end"), axis=1)
        df = df.dropna(subset=["snpId_bin_{}".format(nb_bin)])

        df_sum = pd.DataFrame.from_records(df["snpId_bin_{}".format(nb_bin)].values).sum()
        df_sum.index = [e + 1 for e in list(df_sum.index)]
        df_sum = pd.DataFrame(df_sum).reset_index()
        df_sum.columns = ["Bin_num", "Total"]
        df_sum["Ratio"] = 100 * (df_sum["Total"] / df_sum["Total"].sum())
        print(df_sum)
        # exit()
        df_sum.to_excel("/gstock/EXOTIC/data/VARIATIONS/pathogenic_distribution_missense.xlsx")

    # @staticmethod
    def distribution_pathogenic_variations(self, refseq_path, clinvar, refseq):
        print(clinvar)
        # print(clinvar.CLNSIG.unique())
        # print(clinvar.CLNREVSTAT.unique())
        # exit()

        clinvar_data = clinvar.dropna(subset=["MC"])
        # exit()
        clinvar_data["Gene"] = clinvar_data["GENEINFO"].apply(lambda r: r.split(":")[0])
        clinvar_data["MC_lite"] = clinvar_data["MC"].apply(lambda r: r[0])
        clinvar_data["MC_lite"] = clinvar_data["MC_lite"].apply(lambda r: r.split("|")[1])
        # print(clinvar_data)
        # print(clinvar_data.MC_lite.unique())
        # exit()

        clinvar_data = clinvar_data.loc[clinvar_data["alleles"].str.len() == 2]

        # print(clinvar_data)
        clinvar_data = clinvar_data.loc[clinvar_data["MC_lite"].str.contains("missense")]
        # print(clinvar_data)

        # refseq_with_genes_coord = (
        #     pd.read_parquet(refseq_path)[["Gene", "Gene_start", "Gene_stop", "Strand"]].drop_duplicates().sort_values(by="Gene")
        # )
        # refseq_with_genes_coord["Strand"] = refseq_with_genes_coord["Strand"].replace({"-": 0, "+": 1})

        # refseq_with_genes_coord = refseq_with_genes_coord.loc[refseq_with_genes_coord["Gene"].isin(clinvar_data.Gene.unique().tolist())]

        # refseq_with_genes_coord = self.compute_bins(refseq_with_genes_coord, nb_bins=[20])

        # exit()
        # print(refseq_with_genes_coord)
        # exit()

        clinvar_omim = pd.merge(refseq, clinvar_data, on=["Gene"])
        # print(clinvar_omim)
        # exit()

        nb_bin = 20
        # for nb_bin in nb_bins:
        # clinvar_omim["snpId_bin_{}".format(nb_bin)] = clinvar_omim.head(10).apply(lambda r: self.compute_snv_bin_fct(r, nb_bin), axis=1)
        # clinvar_omim["snpId_bin_{}".format(nb_bin)] = clinvar_omim.parallel_apply(lambda r: self.compute_snv_bin_fct(r, nb_bin), axis=1)
        clinvar_omim["snpId_bin"] = clinvar_omim.parallel_apply(lambda r: self.compute_snv_bin_fct_dev(r), axis=1)

        clinvar_omim = clinvar_omim.dropna(subset=["snpId_bin"])
        df_sum = pd.DataFrame.from_records(clinvar_omim["snpId_bin"].values).sum()
        df_sum.index = [e + 1 for e in list(df_sum.index)]
        df_sum = pd.DataFrame(df_sum).reset_index()
        df_sum.columns = ["Bin_num", "Total"]
        df_sum["Ratio"] = 100 * (df_sum["Total"] / df_sum["Total"].sum())
        # print(df_sum)
        # exit()
        df_sum.to_excel("/gstock/EXOTIC/data/VARIATIONS/pathogenic_distribution_missense_exons.xlsx", index=False)

    @staticmethod
    def load_clinvar(clinvar_path):
        clinvar = pd.read_parquet(clinvar_path)
        print(clinvar)
        # print(clinvar.Gene.nunique())
        # exit()
        # clinvar = clinvar.dropna(subset=["CLNVI"])
        # clinvar["OMIM_variant_id"] = clinvar["CLNVI"].apply(lambda r: [e for e in r if "OMIM" in e])
        # clinvar["OMIM_variant_id"] = clinvar["OMIM_variant_id"].apply(
        #     lambda r: [sub_e.replace("OMIM_Allelic_Variant:", "") for e in r for sub_e in e.split("|") if "OMIM" in sub_e]
        # )
        # clinvar = clinvar.loc[clinvar["OMIM_variant_id"].str.len() > 0]
        # clinvar["OMIM_variant_nb"] = clinvar["OMIM_variant_id"].apply(lambda r: int(r[0].split(".")[1]))
        # clinvar["OMIM"] = clinvar["OMIM_variant_id"].apply(lambda r: int(r[0].split(".")[0]))
        # print(clinvar)
        return clinvar

    @staticmethod
    def load_refseq_pc_genes(path):
        refseq_pc_genes = pd.read_parquet(path)
        refseq_pc_genes["Dbxref"] = refseq_pc_genes["Attributes"].apply(lambda r: [e for e in r.split(";") if "Dbxref" in e][0])
        refseq_pc_genes["Gene"] = refseq_pc_genes["Attributes"].apply(
            lambda r: [e.replace("gene=", "") for e in r.split(";") if "gene=" in e][0]
        )
        refseq_pc_genes["HGNC"] = refseq_pc_genes["Dbxref"].apply(
            lambda r: [e.replace("HGNC:HGNC:", "") for e in r.split(",") if "HGNC" in e]
        )
        refseq_pc_genes["HGNC"] = refseq_pc_genes["HGNC"].apply(lambda r: r[0] if r else np.nan)
        # refseq_pc_genes["HGNC"] = refseq_pc_genes["HGNC"].astype(int)
        refseq_pc_genes["OMIM"] = refseq_pc_genes["Dbxref"].apply(lambda r: [e.replace("MIM:", "") for e in r.split(",") if "MIM" in e])
        refseq_pc_genes["OMIM"] = refseq_pc_genes["OMIM"].apply(lambda r: r[0] if r else np.nan)
        # refseq_pc_genes["OMIM"] = refseq_pc_genes["OMIM"].astype(int)
        refseq_pc_genes["GeneID"] = refseq_pc_genes["Dbxref"].apply(
            lambda r: [e.replace("Dbxref=GeneID:", "") for e in r.split(",") if "GeneID" in e]
        )
        refseq_pc_genes["GeneID"] = refseq_pc_genes["GeneID"].apply(lambda r: r[0] if r else np.nan)

        return refseq_pc_genes

    @staticmethod
    def load_refseq_final(path):
        refseq = pd.read_parquet(path)
        refseq["MAP"] = refseq["Gene"] + "_" + refseq["ranges"]
        return refseq

    @staticmethod
    def load_dext(dext_path, refseq_pc_genes, output_path):
        if os.path.isfile(output_path) is False:

            bins = [round(e, 2) for e in list(np.arange(0, 1.1, 0.1))]
            labels = convert_bins_into_labels(bins)

            dext = pd.read_parquet(dext_path)
            dext = pd.merge(dext, refseq_pc_genes[["Gene", "OMIM"]], on="Gene")
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
            dext.to_parquet(output_path)
        else:
            dext = pd.read_parquet(output_path)

        return dext

    @staticmethod
    def download_omim(df, output_dir):
        print("# Download OMIM API entries into pickle files / File : {}".format(output_dir))

        check_downloaded_omim = [str(e.replace(".pkl", "")) for e in os.listdir(output_dir)]

        j = 0

        omim_api = "https://api.omim.org/api/entry?mimNumber={}&include={}&format=json&apiKey={}"
        api_key = "UMjcVQ3aRlykG5I-NlxZrA"

        stats = collections.defaultdict(int)
        # stats["Genes_to_dl"] = list()
        # stats["Phenos_to_dl"] = list()
        df = df.dropna(subset=["OMIM"])
        genes = sorted(df.OMIM.unique().tolist())

        for gene in tqdm(genes):
            # print(gene)

            # LOCAL VERSION
            if str(gene) not in check_downloaded_omim:
                print("DL", gene)
                json = requests.get(omim_api.format(str(gene), "all", api_key)).json()
                _pickle.dump(json, open(output_dir + "{}.pkl".format(gene), "wb"))
                j += 1
                stats["Genes_to_dl_nb"] += 1
                # stats["Genes_to_dl"].append(gene)
                if j % 20 == 0:
                    print("Sleep ... / Gene level")
                    time.sleep(5)

            else:
                # print("Already dl")
                json = _pickle.load(open(output_dir + "{}.pkl".format(str(gene)), "rb"))
                stats["Genes_already_dl"] += 1

            entry = dict(json["omim"]["entryList"][0]["entry"])

            if "geneMap" in entry:
                if "phenotypeMapList" in entry["geneMap"]:
                    phenotypes = json["omim"]["entryList"][0]["entry"]["geneMap"]["phenotypeMapList"]
                    for pheno in phenotypes:
                        if "phenotypeMimNumber" in pheno["phenotypeMap"]:
                            phenotype_mim_number = pheno["phenotypeMap"]["phenotypeMimNumber"]
                            if str(phenotype_mim_number) not in check_downloaded_omim:
                                print("DL", gene, phenotype_mim_number)
                                second_request = requests.get(omim_api.format(str(phenotype_mim_number), "all", api_key)).json()
                                _pickle.dump(second_request, open(output_dir + "{}.pkl".format(phenotype_mim_number), "wb"))
                                j += 1
                                stats["Phenos_to_dl"] += 1
                                if j % 20 == 0:
                                    time.sleep(5)
                                    print("Sleep ... / Pheno level")
                            else:
                                stats["Phenos_already_dl"] += 1
                                pass
                                # second_request = _pickle.load(open(output_dir + "{}.pkl".format(phenotype_mim_number), "rb"))
                                # pheno_already_download.append(int(phenotype_mim_number))

    @staticmethod
    def build_omim_mp(gene, l, i):
        output_dir = "/gstock/biolo_datasets/variation/benchmark/Databases/OMIM/JSON_API/"

        omim_global_dict = dicts["omim_global_dict"]
        try:

            json = _pickle.load(open(output_dir + "{}.pkl".format(str(gene)), "rb"))

            entry = dict(json["omim"]["entryList"][0]["entry"])
            if "geneMap" in entry:
                if "phenotypeMapList" in entry["geneMap"]:
                    phenotypes = json["omim"]["entryList"][0]["entry"]["geneMap"]["phenotypeMapList"]
                    for pheno in phenotypes:

                        if "phenotypeMimNumber" in pheno["phenotypeMap"]:
                            phenotype_mim_number = pheno["phenotypeMap"]["phenotypeMimNumber"]
                            second_request = _pickle.load(open(output_dir + "{}.pkl".format(phenotype_mim_number), "rb"))

                            for key in second_request["omim"]["entryList"][0]:
                                # print(key)
                                prefered_title = second_request["omim"]["entryList"][0]["entry"]["titles"]["preferredTitle"]

                                if "clinicalSynopsis" in dict(second_request["omim"]["entryList"][0]["entry"]).keys():
                                    clinical_synopsis = dict(second_request["omim"]["entryList"][0]["entry"]["clinicalSynopsis"])

                                    # # FIRST LEVEL - GLOBAL CLINICAL SYNOPSIS BODY PART OMIM
                                    new_dict = collections.defaultdict(list)
                                    for k, v in clinical_synopsis.items():

                                        if "Exists" not in k and k in omim_global_dict:
                                            new_dict[k].append(v)

                                    new_dict["OMIM"] = gene
                                    new_dict["Pheno_OMIM"] = phenotype_mim_number
                                    new_dict["Pheno_Name"] = pheno["phenotypeMap"]["phenotype"]
                                    new_dict["Pheno_prefered_title"] = prefered_title

                                    if "referenceList" in dict(second_request["omim"]["entryList"][0]["entry"]).keys():
                                        pmid_list = [
                                            ref["reference"]["pubmedID"]
                                            for ref in second_request["omim"]["entryList"][0]["entry"]["referenceList"]
                                            if "pubmedID" in ref["reference"]
                                        ]
                                        new_dict["PMID"] = pmid_list

                                    l.append(new_dict)
        except FileNotFoundError:
            # print(gene)
            i += 1

    def build_omim(self, path, omim):
        if os.path.isfile(path) is True:
            genes = sorted(omim.OMIM.dropna().unique().tolist())
            # exit()

            m = multiprocessing.Manager()
            l = m.list()
            i = 0
            parmap.starmap(self.build_omim_mp, list(zip(genes)), l, i, pm_pbar=True)
            print(i)
            df = pd.DataFrame(list(l))
            df = df[
                ["OMIM", "Pheno_OMIM", "Pheno_Name", "Pheno_prefered_title", "PMID"]
                + [c for c in df.columns if c not in ["OMIM", "Pheno_OMIM", "Pheno_Name", "Pheno_prefered_title", "PMID"]]
            ]
            df.to_parquet(path)

            # print(df)
            # exit()

        else:
            df = pd.read_parquet(path)
        print(df)
        # print(df)
        print(df.OMIM.nunique())
        # print(omim)
        # print(omim.OMIM.nunique())
        # omim = omim[["Gene", "OMIM"]].drop_duplicates()
        # print(omim.loc[~omim["OMIM"].isin(df.OMIM.unique().tolist())].dropna())
        # exit()

        # print(df)
        return df

    @staticmethod
    def processed_omim(omim):
        omim = omim.dropna(subset=list(omim.columns[6:]), how="all")

        omim_matrix = omim[omim.columns[6:]].apply(lambda r: [1 if e else 0 for e in r], axis=1)
        omim_matrix.name = "Matrix_vector"
        omim_matrix = pd.concat([omim[omim.columns[:5]], omim_matrix], axis=1)
        omim_matrix["Matrix_vector_sum"] = omim_matrix["Matrix_vector"].apply(sum)

        omim = omim.melt(
            id_vars=list(omim.columns)[:6], value_vars=list(omim.columns)[6:], var_name="OMIM_BP", value_name="OMIM_BP_phenotypes"
        ).dropna()

        return omim, omim_matrix

    @staticmethod
    def search_specific_pheno(omim, output_path):

        if os.path.isfile(output_path) is False:

            l = list()

            for gene in tqdm(omim.OMIM.unique()):

                gene_omim = omim.loc[omim["OMIM"] == gene]
                if gene_omim.Pheno_OMIM.nunique() > 1:
                    # print(gene, gene_omim.shape[0])
                    gene_omim_bp = (
                        gene_omim[["OMIM", "Pheno_OMIM", "OMIM_BP"]].groupby(["OMIM", "Pheno_OMIM"])["OMIM_BP"].apply(set).reset_index()
                    )
                    # print(gene_omim_bp)
                    nb_pheno = gene_omim_bp.Pheno_OMIM.nunique()
                    counter_pheno = collections.Counter([sub_e for e in gene_omim_bp["OMIM_BP"].values.tolist() for sub_e in e])
                    # print(counter_pheno)
                    cutoff = 1 if nb_pheno < 4 else math.floor(nb_pheno / 2)
                    counter_pheno_rare = [k for k, v in counter_pheno.items() if v <= cutoff]
                    # print(counter_pheno_rare)
                    #         print(nb_pheno, cutoff)
                    #         print(counter_pheno_rare)
                    l.append(gene_omim.loc[gene_omim["OMIM_BP"].isin(counter_pheno_rare)])

            #         for pheno_rare in counter_pheno_rare:
            #             print(gene_omim.loc[gene_omim['OMIM_BP'] == pheno_rare, ['Pheno_OMIM', 'Pheno_prefered_title', 'OMIM_BP', 'OMIM_BP_phenotypes']].values.tolist())
            omim_multi_pheno_filterred = pd.concat(l)
            # print(omim_multi_pheno_filterred)
            # print(omim_multi_pheno_filterred.Gene.nunique())
            # print(omim_multi_pheno_filterred.OMIM.nunique())
            # print(omim_multi_pheno_filterred.Pheno_OMIM.nunique())
            omim_multi_pheno_filterred.to_parquet(output_path)
        else:
            omim_multi_pheno_filterred = pd.read_parquet(output_path)

        return omim_multi_pheno_filterred

    @staticmethod
    def merge_dext_omim_specific(dext, omim):
        dext_omim_specific = pd.merge(dext, omim, on="OMIM")

        dext_omim_specific = dext_omim_specific.explode("dext_tissues_up").explode("dext_tissues_down")

        dext_omim_specific["gtex_omim_up"] = dext_omim_specific["dext_tissues_up"].map(dicts["mapping_omim_gtex_neurologic"])
        dext_omim_specific["gtex_omim_down"] = dext_omim_specific["dext_tissues_down"].map(dicts["mapping_omim_gtex_neurologic"])
        dext_omim_specific = dext_omim_specific.explode("gtex_omim_up").explode("gtex_omim_down")
        dext_omim_specific.loc[dext_omim_specific["gtex_omim_up"] == dext_omim_specific["OMIM_BP"], "GTEx_OMIM_match_up"] = True
        dext_omim_specific.loc[dext_omim_specific["gtex_omim_down"] == dext_omim_specific["OMIM_BP"], "GTEx_OMIM_match_down"] = True
        dext_omim_specific[["GTEx_OMIM_match_up", "GTEx_OMIM_match_down"]] = dext_omim_specific[
            ["GTEx_OMIM_match_up", "GTEx_OMIM_match_down"]
        ].fillna(False)
        # print(dext_omim_specific.GTEx_OMIM_match_up.value_counts())
        # print(dext_omim_specific.GTEx_OMIM_match_down.value_counts())

        dext_omim_specific = dext_omim_specific.loc[
            # (dext_omim_specific["GTEx_OMIM_match_up"] == True) | (dext_omim_specific["GTEx_OMIM_match_down"] == True)
            (dext_omim_specific["GTEx_OMIM_match_up"] == True)
        ]
        # print(dext_omim_specific.Gene.nunique())
        # print(dext_omim_specific.OMIM.nunique())
        # print(dext_omim_specific.Pheno_OMIM.nunique())
        return dext_omim_specific

    @staticmethod
    def merge_dext_omim_clinvar(dext_omim, clinvar_omim):
        clinvar_omim["OMIM"] = clinvar_omim["OMIM"].astype(str)
        # print(dext_omim)
        # print(clinvar_omim)

        dext_omim_specific_true = dext_omim.loc[(dext_omim["GTEx_OMIM_match_up"] == True) | (dext_omim["GTEx_OMIM_match_down"] == True)]
        # print(dext_omim_specific_true)

        dext_omim_specific_true = pd.merge(dext_omim_specific_true, clinvar_omim, on=["OMIM", "MAP", "Gene"])
        # print(dext_omim_specific_true)
        # print(dext_omim_specific_true.OMIM.nunique())
        # print(dext_omim_specific_true.Pheno_OMIM.nunique())
        # print(dext_omim_specific_true.MAP.nunique())
        # print(dext_omim_specific_true.snpId.nunique())
        # print(list(dext_omim_specific_true.columns))

        dext_omim_specific_true = dext_omim_specific_true.loc[
            # (dext_omim_specific_true["dext_up"] >= 0.5) | (dext_omim_specific_true["dext_down_reversed"] >= 0.5)
            (dext_omim_specific_true["dext_up"] >= 0.5)
        ]
        # print(dext_omim_specific_true)
        # print(dext_omim_specific_true.OMIM.nunique())
        # print(dext_omim_specific_true.Pheno_OMIM.nunique())
        # print(dext_omim_specific_true.MAP.nunique())
        # print(dext_omim_specific_true.snpId.nunique())
        # print(dext_omim_specific_true.Gene.unique().tolist())
        # print(dext_omim_specific_true.Gene.unique())
        # print(list(dext_omim_specific_true.columns))
        # dext_omim_specific_true.to_excel(yaml["5_PATHOGENICITY"]["TMP"]["dext_omim_clinvar"])
        return dext_omim_specific_true

    @staticmethod
    def retrieve_omim_variants_info_fct(r):
        def forbid(r):
            forbidden = [":", ".", "rs", "ss", "dbSNP"]
            check = True if set([True if e not in r else False for e in forbidden]) == {True} else False
            return check

        if r["OMIM"]:
            omim_gene_json = _pickle.load(
                open("/gstock/biolo_datasets/variation/benchmark/Databases/OMIM/JSON_API/{}.pkl".format(str(r["OMIM"])), "rb")
            )

            if "allelicVariantList" in omim_gene_json["omim"]["entryList"][0]["entry"]:
                # print(omim_gene_json["omim"]["entryList"][0]["entry"]["allelicVariantList"])
                # exit()
                # if int(r["OMIM_variant_nb"] - 1) in omim_gene_json["omim"]["entryList"][0]["entry"]["allelicVariantList"]:
                # print(omim_gene_json)
                # exit()
                # print(r["OMIM_variant_nb"])
                # pprint(omim_gene_json["omim"]["entryList"][0]["entry"]["allelicVariantList"])

                if len(omim_gene_json["omim"]["entryList"][0]["entry"]["allelicVariantList"]) >= r["OMIM_variant_nb"]:
                    if "allelicVariant" in omim_gene_json["omim"]["entryList"][0]["entry"]["allelicVariantList"][r["OMIM_variant_nb"] - 1]:
                        variant_json = omim_gene_json["omim"]["entryList"][0]["entry"]["allelicVariantList"][r["OMIM_variant_nb"] - 1][
                            "allelicVariant"
                        ]
                        #             print(variant_json)
                        # pprint(variant_json)

                        # try:
                        if "text" in variant_json:

                            if variant_json["text"].startswith("See {") is False:

                                pheno_ids = [publi for publi in re.findall("\{.*?\}", variant_json["text"]) if forbid(publi) is True]
                                # print(pheno_ids)
                                #             pheno_id = pheno_id[0]
                                pheno_ids = [pheno.replace("{", "").replace("}", "") for pheno in pheno_ids]
                                pheno_ids = [int(pheno) for pheno in pheno_ids if pheno.isdecimal()]
                                publi_ids = [publi for publi in re.findall("\{.*?\}", variant_json["text"]) if ":" in publi]
                                publi_ids = [
                                    sub_publi for publi in publi_ids for sub_publi in publi.split(":")[0].replace("{", "").split(",")
                                ]
                                #             print(r, pheno_ids)

                                #             if pheno_id == r['Pheno_OMIM']:
                                pmids = []
                                try:
                                    pmids = [
                                        omim_gene_json["omim"]["entryList"][0]["entry"]["referenceList"][
                                            int(publi.split(":")[0].replace("{", "")) - 1
                                        ]["reference"]["pubmedID"]
                                        for publi in publi_ids
                                    ]

                                except:
                                    #                     pprint([omim_gene_json['omim']['entryList'][0]['entry']['referenceList'][int(publi.split(':')[0].replace('{', '')) - 1]['reference'] for publi in publi_ids])
                                    pass
                                r["PMIDS_OMIM"] = pmids
                                r["PHENOS_OMIM"] = pheno_ids

                            # except KeyError:
                            # print(variant_json)
                            # exit()
                            # r["PMIDS_OMIM"] = ""
                            # r["PHENOS_OMIM"] = ""
        return r

    def retrieve_omim_variants_info(self, omim_raw, path):
        if os.path.isfile(path) is False:

            omim = omim_raw.dropna(subset=["OMIM", "OMIM_variant_nb"])
            omim[["OMIM", "OMIM_variant_nb"]] = omim[["OMIM", "OMIM_variant_nb"]].astype(int)
            omim = omim.progress_apply(self.retrieve_omim_variants_info_fct, axis=1)
            omim = omim[list(omim_raw.columns) + ["PMIDS_OMIM", "PHENOS_OMIM"]]
            omim.to_parquet(path)
        else:
            omim = pd.read_parquet(path)
        print(omim)
        return omim

    # @staticmethod
    def compute_snv_bin_fct(self, r, nb_bin, col="locus.position"):
        l = [0] * nb_bin
        for j, b in enumerate(r["Gene_bins_{}".format(nb_bin)]):
            if j < nb_bin:
                if b <= r[col] < r["Gene_bins_{}".format(nb_bin)][j + 1]:
                    l[j] += 1
        if r["Strand"] == 0:
            # print("OK")
            # exit()
            # print(l)a
            l = l[::-1]
            # print(l)
            # exit()
        else:
            l = l
        # print(r)
        # print(l)
        # exit()
        return l

    # @staticmethod
    def compute_snv_bin_fct_dev(self, r, nb_bin=20, col="locus.position"):
        l = [0] * nb_bin
        for j, b in enumerate(r["Exons_bins"]):
            if j < nb_bin:
                if b <= r[col] < r["Exons_bins"][j + 1]:
                    l[j] += 1
        if r["Strand"] == 0:
            l = l[::-1]
        else:
            l = l
        return l

    @staticmethod
    def compute_bins(genes_remap, nb_bins=[10]):
        genes_remap["Gene_length"] = genes_remap["Gene_stop"] - genes_remap["Gene_start"]
        for nb in nb_bins:
            genes_remap["Gene_bin_size"] = genes_remap["Gene_length"] / nb

            genes_remap["Gene_bins_{}".format(nb)] = genes_remap.parallel_apply(
                lambda r: [r["Gene_start"] + (round(r["Gene_bin_size"] * j)) for j in list(range(nb + 1))], axis=1
            )
        return genes_remap

    @staticmethod
    def compute_bins_dev(genes_remap, nb_bins=[10]):
        print(genes_remap)
        exit()
        genes_remap["Gene_length"] = genes_remap["Gene_stop"] - genes_remap["Gene_start"]
        for nb in nb_bins:
            genes_remap["Gene_bin_size"] = genes_remap["Gene_length"] / nb

            genes_remap["Gene_bins_{}".format(nb)] = genes_remap.parallel_apply(
                lambda r: [r["Gene_start"] + (round(r["Gene_bin_size"] * j)) for j in list(range(nb + 1))], axis=1
            )
        return genes_remap


if __name__ == "__main__":
    c = Pathogenicity()