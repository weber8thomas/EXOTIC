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
from utils.utils import load_config_file, mkdir
import json
import subprocess
import hail as hl

hl.init(min_block_size=128)
os.environ["PYSPARK_SUBMIT_ARGS"] = "--driver-memory 200 pyspark-shell"


from pandarallel import pandarallel

pandarallel.initialize(nb_workers=20, progress_bar=True)

## YAML FILES CONFIG
yaml = load_config_file(config_file="/home/weber/PycharmProjects/EXOTIC/clean/src/config_clean_clean.yaml")


class Variants_CCRS:
    def __init__(
        self,
    ):
        # YAML FILES CONFIG

        self.tmp_ccrs_path = "/gstock/EXOTIC/data/VARIATIONS/CCRS_modified_refseq_corrected.parquet"
        self.tmp_variants_path = "/gstock/EXOTIC/data/VARIATIONS/ClinVar_gnomAD_all_2021_refseq_corrected.parquet"
        self.tmp_clinvar_ccrs_path = "/gstock/EXOTIC/data/VARIATIONS/ClinVar_CCRS_gnomAD_2021_refseq_corrected.parquet"
        self.tmp_phylocsf_path = "/gstock/EXOTIC/data/CONSERVATION/phyloCSF_modified_2021_refseq_corrected.parquet"
        self.gnomad_pathogenic_genes_path = (
            "/gstock/biolo_datasets/variation/variation_sets/gnomAD/EXOME_MISTIC/RAW/gnomad_clinvar_genes.csv.gz"
        )
        # self.gnomad_pathogenic_genes = self.read_gnomad_file()

        # refseq_bed = self.turn_refseq_miso_exons_into_bed(
        #     refseq_path=yaml["2_EXPRESSION"]["Final"]["refseq_corrected_cds_recomputed"],
        #     bed_path=yaml["3_EXONS_PROPERTIES"]["TMP"]["refseq_bed"],
        # )

        # gnomad_data = self.read_gnomad_ht_file(
        #     bed_path=yaml["3_EXONS_PROPERTIES"]["TMP"]["refseq_bed"],
        #     gnomad_path_ht=yaml["3_EXONS_PROPERTIES"]["External"]["gnomad_2_1_1"],
        #     output_file=yaml["3_EXONS_PROPERTIES"]["Final"]["refseq_gnomad_hail_retrieved"],
        # )
        # print(gnomad_data)

        # clinvar_data = self.read_and_process_clinvar_vcf(
        # bed_path=yaml["3_EXONS_PROPERTIES"]["TMP"]["refseq_bed"],
        # clinvar_vcf_path=yaml["3_EXONS_PROPERTIES"]["External"]["clinvar_latest"],
        # output_file=yaml["3_EXONS_PROPERTIES"]["Final"]["refseq_clinvar_hail_retrieved"],
        # )

        # print(clinvar_data)

        # self.read_and_process_ccr(
        #     ccr_autosome=yaml["3_EXONS_PROPERTIES"]["External"]["ccr_autosome"],
        #     ccr_x=yaml["3_EXONS_PROPERTIES"]["External"]["ccr_x"],
        #     ouput_file_ccr_concat=yaml["3_EXONS_PROPERTIES"]["External"]["ccr_complete"],
        #     bed_path=yaml["3_EXONS_PROPERTIES"]["TMP"]["refseq_bed"],
        #     output_file_start=yaml["3_EXONS_PROPERTIES"]["TMP"]["refseq_ccr_start"],
        #     output_file_end=yaml["3_EXONS_PROPERTIES"]["TMP"]["refseq_ccr_end"],
        #     output_file_start_end_concat=yaml["3_EXONS_PROPERTIES"]["Final"]["refseq_ccr_start_end"],
        # )

        self.read_and_process_phylocsf(
            phylocsf=yaml["3_EXONS_PROPERTIES"]["External"]["phylocsf_raw"],
            phylocsf_lite=yaml["3_EXONS_PROPERTIES"]["TMP"]["phylocsf_lite"],
            bed_path=yaml["3_EXONS_PROPERTIES"]["TMP"]["refseq_bed"],
            output_file_start=yaml["3_EXONS_PROPERTIES"]["TMP"]["refseq_phylocsf_start"],
            output_file_end=yaml["3_EXONS_PROPERTIES"]["TMP"]["refseq_phylocsf_end"],
            output_file_start_end_concat=yaml["3_EXONS_PROPERTIES"]["Final"]["refseq_phylocsf"],
        )
        # self.process_ccrs()
        # self.process_clinvar()
        # self.merge_clinvar_ccrs()
        # self.process_phylocsf()
        exit()

        # TEST WITHOUT 0
        merge_df = merge_df.loc[merge_df["CCRS_CCR_percentile"] > 0]

        # MELT
        # merge_df = merge_df.melt(id_vars='Const_Alt', value_vars='CCRS_CCR_percentile')
        print(merge_df)

    @staticmethod
    def turn_refseq_miso_exons_into_bed(refseq_path, bed_path):
        print("# Turn RefSeq into BED / File : {}".format(bed_path))
        if os.path.isfile(bed_path) is False:
            df = pd.read_parquet(refseq_path)
            df = df.loc[(df["mRNA_nb_total"] > 1) & (df["CDS_count"] > 1)]

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
            df = df[["CHROM", "Start", "End", "MAP"]].sort_values(by=["CHROM", "Start"])
            df["Start"] = df["Start"] - 1
            df["End"] = df["End"] - 1
            # bed_path = yaml['3_EXONS_PROPERTIES']['TMP']['refseq_bed_gnomad']
            mkdir(os.path.dirname(bed_path))
            df.to_csv(bed_path, header=False, sep="\t", index=False)
        else:
            pass

    @staticmethod
    def read_gnomad_ht_file(bed_path, gnomad_path_ht, output_file):
        print("# Load gnomAD hail table & retrieve corresponding variations / File : {}".format(output_file.replace(".tsv.gz", ".parquet")))

        if os.path.isfile(output_file.replace(".tsv.gz", ".parquet")) is False:
            # import hail as hl

            # hl.init(min_block_size=128)
            # os.environ["PYSPARK_SUBMIT_ARGS"] = "--driver-memory 200 pyspark-shell"

            data = hl.read_table(gnomad_path_ht)

            bed_file = hl.import_bed(bed_path, reference_genome="GRCh37")
            print(bed_file.head(10).show())

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
        return data_lite

        # data_lite = data_lite.key_by('id')

    @staticmethod
    def read_and_process_clinvar_vcf(bed_path, clinvar_vcf_path, output_file):
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

            bed_file = hl.import_bed(bed_path, reference_genome="GRCh37")

            hl.import_vcf(clinvar_vcf_path, force_bgz=True).write(clinvar_vcf_path.replace(".vcf.gz", ".ht"), overwrite=True)
            clinvar = hl.read_matrix_table("/gstock/biolo_datasets/variation/variation_sets/clinvar/vcf_GRCh37/v2/clinvar_20210123.ht")
            clinvar_lite = clinvar.filter_rows(hl.is_defined(bed_file[clinvar.locus]))
            clinvar_lite = clinvar_lite.filter_rows(clinvar_lite.info.CLNVC == "single_nucleotide_variant")
            clinvar_lite = clinvar_lite.select_rows(
                clinvar_lite.rsid,
                clinvar_lite.info.ALLELEID,
                clinvar_lite.info.CLNREVSTAT,
                clinvar_lite.info.CLNSIG,
                clinvar_lite.info.GENEINFO,
                clinvar_lite.info.RS,
            )

            clinvar_lite = clinvar_lite.annotate_rows(
                MAP=bed_file[clinvar_lite.locus].target,
                Gene=hl.str(bed_file[clinvar_lite.locus].target).split("_")[0],
            )

            clinvar_lite_pd_raw = clinvar_lite.make_table().to_pandas()
            clinvar_lite_pd = clinvar_lite_pd_raw.copy()
            clinvar_lite_pd = clinvar_lite_pd.replace(to_replace="None", value=np.nan)
            clinvar_lite_pd = clinvar_lite_pd.dropna(subset=["ALLELEID", "CLNREVSTAT", "CLNSIG", "GENEINFO"], how="any")
            clinvar_lite_pd["CLNREVSTAT"] = clinvar_lite_pd["CLNREVSTAT"].apply(lambda r: ",".join(r))
            clinvar_lite_pd["CLNSIG"] = clinvar_lite_pd["CLNSIG"].apply(lambda r: ",".join(r))
            clinvar_lite_pd["RS_STARS"] = clinvar_lite_pd["CLNREVSTAT"].map(clinvar_review_status)
            clinvar_lite_pd = clinvar_lite_pd.loc[
                (clinvar_lite_pd["CLNSIG"].str.contains("athogenic"))
                & (~clinvar_lite_pd["CLNSIG"].str.contains("Conflicting_interpretations_of_pathogenicity"))
            ]
            clinvar_lite_pd = clinvar_lite_pd.loc[~clinvar_lite_pd["CLNREVSTAT"].str.contains("interpret")]
            clinvar_lite_pd.to_parquet(output_file)
        else:
            clinvar_lite_pd = pd.read_parquet(output_file)
        return clinvar_lite_pd

    def read_and_process_ccr(
        self, ccr_autosome, ccr_x, ouput_file_ccr_concat, bed_path, output_file_start, output_file_end, output_file_start_end_concat
    ):
        if os.path.isfile(ouput_file_ccr_concat) is False:
            cmd = [
                "zcat",
                "{}".format(),
                "{}".format(),
                "|",
                "sed",
                "'s/^chr//g'",
                "|",
                "/biolo/ngs/tabix/bgzip",
                ">",
                "{}".format(yaml["3_EXONS_PROPERTIES"]["External"]["ccr_complete"]),
            ]

            ps = subprocess.Popen(" ".join(cmd), shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
            output = ps.communicate()[0]
        else:
            print("Already done")
            pass

        if os.path.isfile(output_file_start) is False:
            pass

            bed_file = hl.import_bed(bed_path, reference_genome="GRCh37")

            ccrs_table = hl.import_bed(
                ouput_file_ccr_concat,
                reference_genome="GRCh37",
                force_bgz=True,
            )
            ccrs_table = ccrs_table.annotate(
                start=ccrs_table.interval.start,
                end=ccrs_table.interval.end,
            )
            ccrs_table_start = ccrs_table.filter(hl.is_defined(bed_file[ccrs_table.start]))
            ccrs_table_start = ccrs_table_start.annotate(
                MAP=bed_file[ccrs_table_start.start].target,
                Gene=hl.str(bed_file[ccrs_table_start.start].target).split("_")[0],
            )
            ccrs_table_start.export(output_file_start)
            ccrs_table_start.to_parquet(output_file_start.replace(".tsv.gz", ".parquet"))

        else:
            ccrs_table_start = pd.read_parquet(output_file_start.replace(".tsv.gz", ".parquet"))
            print(ccrs_table_start)

        if os.path.isfile(output_file_end) is False:

            bed_file = hl.import_bed(bed_path, reference_genome="GRCh37")

            ccrs_table = hl.import_bed(
                ouput_file_ccr_concat,
                reference_genome="GRCh37",
                force_bgz=True,
            )
            ccrs_table = ccrs_table.annotate(
                start=ccrs_table.interval.start,
                end=ccrs_table.interval.end,
            )
            ccrs_table_end = ccrs_table.filter(hl.is_defined(bed_file[ccrs_table.end]))
            ccrs_table_end = ccrs_table_end.annotate(
                MAP=bed_file[ccrs_table_end.end].target,
                Gene=hl.str(bed_file[ccrs_table_end.end].target).split("_")[0],
            )
            ccrs_table_end.export(output_file_end)
            ccrs_table_end = pd.read_csv(output_file_end, compression="gzip", sep="\t")
            ccrs_table_end.to_parquet(output_file_end.replace(".tsv.gz", ".parquet"))

        else:
            # ccrs_table_end = pd.read_csv(output_file_end, compression="gzip", sep="\t")
            # ccrs_table_end.to_parquet(output_file_end.replace(".tsv.gz", ".parquet"))

            ccrs_table_end = pd.read_parquet(output_file_end.replace(".tsv.gz", ".parquet"))
            print(ccrs_table_end)

        if os.path.isfile(output_file_start_end_concat) is False:
            concat_ccrs_start_end = (
                pd.concat([ccrs_table_start, ccrs_table_end]).sort_values(by="interval").drop_duplicates(subset="interval")
            )
            concat_ccrs_start_end.to_parquet(output_file_start_end_concat)
        else:
            concat_ccrs_start_end = pd.read_parquet(output_file_start_end_concat)

    def read_and_process_phylocsf(
        self, phylocsf, phylocsf_lite, bed_path, output_file_start, output_file_end, output_file_start_end_concat
    ):
        if os.path.isfile(phylocsf_lite) is False:
            print("# Process raw phyloCSF data")

            phylocsf_raw = pd.read_csv(phylocsf, compression="gzip", sep="\t").rename(
                {"chromosome_name": "CHROM", "start_coordinate": "Start", "end_coordinate": "End", "max_score": "phyloCSF_score"}, axis=1
            )
            phylocsf_lite_df = phylocsf_raw[["CHROM", "Start", "End", "phyloCSF_score"]]
            phylocsf_lite_df["CHROM"] = phylocsf_lite_df["CHROM"].str.replace("chr", "")
            phylocsf_lite_df = phylocsf_lite_df.drop_duplicates().sort_values(by=["CHROM", "Start"])
            phylocsf_lite_df.to_csv(phylocsf_lite, header=False, sep="\t", index=False)
            print(phylocsf_raw)
            print(phylocsf_lite_df)

        else:
            print("# Process raw phyloCSF data / Already done")
            pass

        if os.path.isfile(output_file_start) is False:
            print("# Process phyloCSF data with hail to map on RefSeq (Start)")

            bed_file = hl.import_bed(bed_path, reference_genome="GRCh37")

            phylocsf_table = hl.import_bed(
                phylocsf_lite,
                reference_genome="GRCh37",
                # force_bgz=True,
            )

            phylocsf_table = phylocsf_table.annotate(
                start=phylocsf_table.interval.start,
                end=phylocsf_table.interval.end,
            )
            phylocsf_table_start = phylocsf_table.filter(hl.is_defined(bed_file[phylocsf_table.start]))
            phylocsf_table_start = phylocsf_table_start.annotate(
                MAP=bed_file[phylocsf_table_start.start].target,
                Gene=hl.str(bed_file[phylocsf_table_start.start].target).split("_")[0],
            )
            phylocsf_table_start.export(output_file_start)
            phylocsf_table_start = pd.read_csv(output_file_start, compression="gzip", sep="\t")
            phylocsf_table_start.to_parquet(output_file_start.replace(".tsv.gz", ".parquet"))

        else:
            phylocsf_table_start = pd.read_parquet(output_file_start.replace(".tsv.gz", ".parquet"))
        print(phylocsf_table_start)

        if os.path.isfile(output_file_end) is False:
            print("# Process phyloCSF data with hail to map on RefSeq (End)")

            bed_file = hl.import_bed(bed_path, reference_genome="GRCh37")

            phylocsf_table = hl.import_bed(
                phylocsf_lite,
                reference_genome="GRCh37",
                # force_bgz=True,
            )

            phylocsf_table = phylocsf_table.annotate(
                start=phylocsf_table.interval.start,
                end=phylocsf_table.interval.end,
            )
            phylocsf_table_end = phylocsf_table.filter(hl.is_defined(bed_file[phylocsf_table.end]))
            phylocsf_table_end = phylocsf_table_end.annotate(
                MAP=bed_file[phylocsf_table_end.end].target,
                Gene=hl.str(bed_file[phylocsf_table_end.end].target).split("_")[0],
            )
            phylocsf_table_end.export(output_file_end)
            phylocsf_table_end = pd.read_csv(output_file_end, compression="gzip", sep="\t")
            phylocsf_table_end.to_parquet(output_file_end.replace(".tsv.gz", ".parquet"))

        else:
            phylocsf_table_end = pd.read_parquet(output_file_end.replace(".tsv.gz", ".parquet"))
        print(phylocsf_table_end)

        if os.path.isfile(output_file_start_end_concat) is False:
            print("Concat phyloCSF Start & End")
            concat_phylocsf_start_end = (
                pd.concat([phylocsf_table_start, phylocsf_table_end]).sort_values(by="interval").drop_duplicates(subset="interval")
            )
            concat_phylocsf_start_end.to_parquet(output_file_start_end_concat)
        else:
            concat_phylocsf_start_end = pd.read_parquet(output_file_start_end_concat)
        print(concat_phylocsf_start_end)

    def read_gnomad_file(self):
        gnomad_clinvar_genes = pd.read_csv(
            self.gnomad_pathogenic_genes_path,
            compression="gzip",
            sep="\t",
            low_memory=False,
            # nrows=100000,
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
            (gnomad_clinvar_genes["REF"].str.len() < 50) & (gnomad_clinvar_genes["ALT"].str.len() < 50)
        ]
        gnomad_clinvar_genes = gnomad_clinvar_genes.rename({"Genes": "Gene"}, axis=1)
        gnomad_clinvar_genes["Status"] = "Benign"

        print(gnomad_clinvar_genes)

        # gnomad_clinvar_genes = gnomad_clinvar_gen es.loc[
        #     (gnomad_clinvar_genes["ALT"].str.len() == 1)
        #     & (gnomad_clinvar_genes["REF"].str.len() == 1)
        # ]
        return gnomad_clinvar_genes

    def init_mp_test(self, gene, l_ccrs, refseq, ccrs):
        refseq_gene = refseq.loc[refseq["Gene"] == gene]
        ccrs_gene = ccrs.loc[ccrs["CCRS_Gene"] == gene]
        l_ccrs.append(self.mapping_ccrs(refseq_gene, ccrs_gene))

    def process_ccrs(self):
        if os.path.isfile(self.tmp_ccrs_path) is True:
            print("CCRS")

            ## REFSEQ DETAILED
            refseq = pd.read_csv(
                # self.exotic_files["refseq_path"],
                self.files["RefSeq"]["refseq_corrected_lite"].replace(".parquet", ".csv.gz"),
                compression="gzip",
                sep="\t",
                low_memory=False,
                # nrows=1000,
            ).sort_values(by=["Gene", "ranges"])

            refseq.columns = [c.replace("new_", "") for c in refseq.columns]

            print(refseq)
            print(list(refseq.columns))
            # refseq = refseq.dropna(subset=["HGNC"])
            # refseq["HGNC"] = refseq["HGNC"].astype(int)
            refseq = refseq.loc[(refseq["mRNA_nb_total"] > 1) & (refseq["CDS_count"] > 1)]
            # refseq = refseq.head(1000)

            # CCRS
            ccrs = pd.read_parquet(self.exotic_files["refseq_ccrs"])
            # ccrs = ccrs.head(1000)

            l_ccrs = list()

            # m = multiprocessing.Manager()
            # l_ccrs = m.list()
            # parmap.starmap(self.init_mp_test, list(zip(refseq.Gene.unique()))[:100], l_ccrs, refseq, ccrs, pm_pbar=True)

            for gene in tqdm(refseq.Gene.unique()):
                refseq_gene = refseq.loc[refseq["Gene"] == gene]
                ccrs_gene = ccrs.loc[ccrs["CCRS_Gene"] == gene]
                l_ccrs.append(self.mapping_ccrs(refseq_gene, ccrs_gene))

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
            labels_ratio = [str(round(labels[j], 1)) + " - " + str(round(labels[j + 1], 1)) for j in range(len(labels) - 1)]
            merge_df["Ratio_num_bins"] = pd.cut(
                merge_df["Ratio_num"],
                bins=bins,
                labels=labels_ratio,
                include_lowest=True,
            )
            merge_df.to_parquet(self.tmp_ccrs_path, index=False)

        else:
            merge_df = pd.read_parquet(self.tmp_ccrs_path)
        return merge_df

    def process_clinvar(self):
        if os.path.isfile(self.tmp_variants_path) is True:
            print("Variants")

            ## REFSEQ DETAILED
            refseq = pd.read_csv(
                # self.exotic_files["refseq_path"],
                self.files["RefSeq"]["refseq_corrected_lite"].replace(".parquet", ".csv.gz"),
                compression="gzip",
                sep="\t",
                low_memory=False,
                # nrows=1000,
            ).sort_values(by=["Gene", "ranges"])
            refseq.columns = [c.replace("new_", "") for c in refseq.columns]

            # refseq = refseq.dropna(subset=["HGNC"])
            # refseq["HGNC"] = refseq["HGNC"].astype(int)
            # refseq["Ratio_num"] = refseq["Ratio"].apply(eval)
            refseq = refseq.loc[(refseq["mRNA_nb_total"] > 1) & (refseq["CDS_count"] > 1)]
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

            validated_clinvar = validated_clinvar.loc[(validated_clinvar["ALT"].str.len() == 1) & (validated_clinvar["REF"].str.len() == 1)]
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
                gnomad_gene = self.gnomad_pathogenic_genes.loc[self.gnomad_pathogenic_genes["Gene"] == gene]
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
            labels_ratio = [str(round(labels[j], 1)) + " - " + str(round(labels[j + 1], 1)) for j in range(len(labels) - 1)]
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
        if os.path.isfile(self.tmp_phylocsf_path) is True:
            print("phyloCSF")

            ## REFSEQ DETAILED
            refseq = pd.read_csv(
                # self.exotic_files["refseq_path"],
                self.files["RefSeq"]["refseq_corrected_lite"].replace(".parquet", ".csv.gz"),
                compression="gzip",
                sep="\t",
                low_memory=False,
                # nrows=1000,
            ).sort_values(by=["Gene", "ranges"])

            refseq.columns = [c.replace("new_", "") for c in refseq.columns]

            # refseq = refseq.dropna(subset=["HGNC"])
            # refseq["HGNC"] = refseq["HGNC"].astype(int)
            # refseq["Ratio_num"] = refseq["Ratio"].apply(eval)
            refseq = refseq.loc[(refseq["mRNA_nb_total"] > 1) & (refseq["CDS_count"] > 1)]
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
            labels_ratio = [str(round(labels[j], 1)) + " - " + str(round(labels[j + 1], 1)) for j in range(len(labels) - 1)]
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
            # print(clinvar.Real_Status.value_counts())
            # exit()

            # clinvar = clinvar.head(100)

            # CCRS
            ccrs = pd.read_parquet(self.tmp_ccrs_path)
            bins = [0, 20, 80, 90, 95, 99, 100]
            labels = bins.copy()
            labels_ratio = [str(round(labels[j], 1)) + " - " + str(round(labels[j + 1], 1)) for j in range(len(labels) - 1)]
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
            # clinvar_gene = clinvar.loc[clinvar["Gene"] == gene]
            # ccrs_gene = ccrs.loc[ccrs["Gene"] == gene]
            # l_ccrs.append(self.mapping_ccrs_clinvar(clinvar_gene, ccrs_gene))

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
            s = ccrs_gene["CCRS_Start"].between(int(r.split("-")[0]), int(r.split("-")[1]), inclusive=True)
            e = ccrs_gene["CCRS_End"].between(int(r.split("-")[0]), int(r.split("-")[1]), inclusive=True)
            concat = pd.concat([s, e], axis=1)
            concat = concat.loc[(concat["CCRS_Start"] == True) & (concat["CCRS_End"] == True)]
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
            s = phyloCSF_gene["Start"].between(int(r.split("-")[0]), int(r.split("-")[1]), inclusive=True)
            e = phyloCSF_gene["End"].between(int(r.split("-")[0]), int(r.split("-")[1]), inclusive=True)
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
            m = clinvar_gene["POS"].between(int(ranges.split("-")[0]), int(ranges.split("-")[1]), inclusive=True)
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
                return l_ccrs.append(pd.concat(tmp_l))
            else:
                return pd.DataFrame()
        else:
            return pd.DataFrame()

    def mapping_clinvar(self, refseq_gene, clinvar_gene):
        def map_variant(r, tmp_l):
            m = clinvar_gene["POS"].between(int(r.split("-")[0]), int(r.split("-")[1]), inclusive=True)
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


if __name__ == "__main__":
    c = Variants_CCRS()
