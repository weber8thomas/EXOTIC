import sys

sys.path.append("/home/weber/PycharmProjects/EXOTIC/src")

import math
import os
from pprint import pprint
from scipy import stats
import re
import json
import numpy as np
import collections

import pandas as pd

pd.options.mode.chained_assignment = None  # default='warn'

from tqdm import tqdm

tqdm.pandas()

from utils.utils import load_config_file, mkdir

## YAML FILES CONFIG
yaml = load_config_file(config_file="clean/src/config_clean_clean.yaml")

## JSON DICTS CONFIG
dicts = json.load(open("clean/src/config/EXOTIC_config.json"))


from pandarallel import pandarallel

pandarallel.initialize(nb_workers=20, progress_bar=True)


class HandleExonExpression:
    def __init__(self):

        pext_groupby = self.transform_pext(yaml["4_DEXT"]["TMP"]["pext_groupby_refseq"])
        dext = self.build_file_dext_score(
            pext_file=pext_groupby,
            refseq_path=yaml["2_EXPRESSION"]["Final"]["refseq_corrected_cds_recomputed"],
            biomart_path=yaml["1_GENOMICS"]["External"]["biomart_tmp"],
            output_path=yaml["4_DEXT"]["Final"]["dext"],
        )

        dext = dext.sort_values(by=["Gene", "Start", "End"])
        # self.test(dext)
        # exit()

        self.compare_to_utrs(
            dext,
            yaml["1_GENOMICS"]["Final"]["refseq_miso_utrs"],
            yaml["1_GENOMICS"]["TMP"]["tmp_refseq_exons"],
            yaml["1_GENOMICS"]["TMP"]["tmp_refseq_cds"],
        )

    @staticmethod
    def test(df):
        l = list()
        previous_gene = ""
        for gene in tqdm(df.Gene.unique().tolist()):
            # print(gene)
            start, end, previous_start, previous_end = 0, 0, 0, 0
            tmp_df = df.loc[df["Gene"] == gene]
            for j, r in tmp_df.iterrows():
                if j == 0:
                    previous_start = r["Start"]
                    previous_end = r["End"]
                elif j > 0:
                    start = r["Start"]
                    end = r["End"]
                    if previous_start < end:
                        if previous_end > start:
                            l.append({"G": gene, "PS": previous_start, "PE": previous_end, "S": start, "E": end})

                    previous_start = r["Start"]
                    previous_end = r["End"]
        stats = pd.DataFrame(l)
        print(stats)
        print(stats.G.nunique())

    @staticmethod
    def transform_pext(output_file):
        print("# Process pext raw file & groupby exon / File : {}".format(output_file.replace(".tsv.gz", ".parquet")))

        if os.path.isfile(output_file.replace(".tsv.gz", ".parquet")) is False:
            print("# Files don't exist ☒")

            mkdir(output_file)

            import hail as hl

            hl.init(min_block_size=128)

            os.environ["PYSPARK_SUBMIT_ARGS"] = "--driver-memory 200 pyspark-shell"

            gs = yaml["4_DEXT"]["External"]["pext_hail_table"]
            # gs = "/data/scratch/gnomAD/v3/gnomad.genomes.v3.1.sites.ht/"
            data = hl.read_table(gs)

            refseq_bed_path = yaml["3_EXONS_PROPERTIES"]["TMP"]["refseq_bed"]
            bed_file = hl.import_bed(refseq_bed_path, reference_genome="GRCh37")
            # bed_file = bed_file.head(10)
            # bed_file.show()

            data_lite = data.filter(hl.is_defined(bed_file[data.locus]))
            data_lite = data_lite.annotate(
                MAP=bed_file[data_lite.locus].target,
                Gene=hl.str(bed_file[data_lite.locus].target).split("_")[0],
            )
            # data_lite.show()

            data_groupby = data_lite.group_by("MAP").aggregate(
                Adipose_Subcutaneous=hl.agg.mean(data_lite.Adipose_Subcutaneous),
                Adipose_Visceral_Omentum_=hl.agg.mean(data_lite.Adipose_Visceral_Omentum_),
                AdrenalGland=hl.agg.mean(data_lite.AdrenalGland),
                Artery_Aorta=hl.agg.mean(data_lite.Artery_Aorta),
                Artery_Coronary=hl.agg.mean(data_lite.Artery_Coronary),
                Artery_Tibial=hl.agg.mean(data_lite.Artery_Tibial),
                Bladder=hl.agg.mean(data_lite.Bladder),
                Brain_Amygdala=hl.agg.mean(data_lite.Brain_Amygdala),
                Brain_Anteriorcingulatecortex_BA24_=hl.agg.mean(data_lite.Brain_Anteriorcingulatecortex_BA24_),
                Brain_Caudate_basalganglia_=hl.agg.mean(data_lite.Brain_Caudate_basalganglia_),
                Brain_CerebellarHemisphere=hl.agg.mean(data_lite.Brain_CerebellarHemisphere),
                Brain_Cerebellum=hl.agg.mean(data_lite.Brain_Cerebellum),
                Brain_Cortex=hl.agg.mean(data_lite.Brain_Cortex),
                Brain_FrontalCortex_BA9_=hl.agg.mean(data_lite.Brain_FrontalCortex_BA9_),
                Brain_Hippocampus=hl.agg.mean(data_lite.Brain_Hippocampus),
                Brain_Hypothalamus=hl.agg.mean(data_lite.Brain_Hypothalamus),
                Brain_Nucleusaccumbens_basalganglia_=hl.agg.mean(data_lite.Brain_Nucleusaccumbens_basalganglia_),
                Brain_Putamen_basalganglia_=hl.agg.mean(data_lite.Brain_Putamen_basalganglia_),
                Brain_Spinalcord_cervicalc_1_=hl.agg.mean(data_lite.Brain_Spinalcord_cervicalc_1_),
                Brain_Substantianigra=hl.agg.mean(data_lite.Brain_Substantianigra),
                Breast_MammaryTissue=hl.agg.mean(data_lite.Breast_MammaryTissue),
                Cells_EBV_transformedlymphocytes=hl.agg.mean(data_lite.Cells_EBV_transformedlymphocytes),
                Cells_Transformedfibroblasts=hl.agg.mean(data_lite.Cells_Transformedfibroblasts),
                Cervix_Ectocervix=hl.agg.mean(data_lite.Cervix_Ectocervix),
                Cervix_Endocervix=hl.agg.mean(data_lite.Cervix_Endocervix),
                Colon_Sigmoid=hl.agg.mean(data_lite.Colon_Sigmoid),
                Colon_Transverse=hl.agg.mean(data_lite.Colon_Transverse),
                Esophagus_GastroesophagealJunction=hl.agg.mean(data_lite.Esophagus_GastroesophagealJunction),
                Esophagus_Mucosa=hl.agg.mean(data_lite.Esophagus_Mucosa),
                Esophagus_Muscularis=hl.agg.mean(data_lite.Esophagus_Muscularis),
                FallopianTube=hl.agg.mean(data_lite.FallopianTube),
                Heart_AtrialAppendage=hl.agg.mean(data_lite.Heart_AtrialAppendage),
                Heart_LeftVentricle=hl.agg.mean(data_lite.Heart_LeftVentricle),
                Kidney_Cortex=hl.agg.mean(data_lite.Kidney_Cortex),
                Liver=hl.agg.mean(data_lite.Liver),
                Lung=hl.agg.mean(data_lite.Lung),
                MinorSalivaryGland=hl.agg.mean(data_lite.MinorSalivaryGland),
                Muscle_Skeletal=hl.agg.mean(data_lite.Muscle_Skeletal),
                Nerve_Tibial=hl.agg.mean(data_lite.Nerve_Tibial),
                Ovary=hl.agg.mean(data_lite.Ovary),
                Pancreas=hl.agg.mean(data_lite.Pancreas),
                Pituitary=hl.agg.mean(data_lite.Pituitary),
                Prostate=hl.agg.mean(data_lite.Prostate),
                Skin_NotSunExposed_Suprapubic_=hl.agg.mean(data_lite.Skin_NotSunExposed_Suprapubic_),
                Skin_SunExposed_Lowerleg_=hl.agg.mean(data_lite.Skin_SunExposed_Lowerleg_),
                SmallIntestine_TerminalIleum=hl.agg.mean(data_lite.SmallIntestine_TerminalIleum),
                Spleen=hl.agg.mean(data_lite.Spleen),
                Stomach=hl.agg.mean(data_lite.Stomach),
                Testis=hl.agg.mean(data_lite.Testis),
                Thyroid=hl.agg.mean(data_lite.Thyroid),
                Uterus=hl.agg.mean(data_lite.Uterus),
                Vagina=hl.agg.mean(data_lite.Vagina),
                WholeBlood=hl.agg.mean(data_lite.WholeBlood),
                mean_proportion=hl.agg.mean(data_lite.mean_proportion),
            )
            # table_result.export()
            data_groupby.export(output_file)
            data_groupby = pd.read_csv(output_file, compression="gzip", sep="\t")
            data_groupby.to_parquet(output_file.replace(".tsv.gz", ".parquet"))

        else:
            print("# Files exist ✓, Loading ... ")

            data_groupby = pd.read_parquet(output_file.replace(".tsv.gz", ".parquet"))
        return data_groupby

    @staticmethod
    def main_transform(x):
        x = x.dropna()

        def sigmoid(e):
            return 1 / (1 + math.exp(-e))

        vsigmoid = np.vectorize(sigmoid)

        def norm(e):
            return (2 * ((e - 0) / (1 - 0))) - 1

        vnorm = np.vectorize(norm)

        def remove_value(x, i):
            return [sub_e for j, sub_e in enumerate(x) if j != i]

        def modif_zscore(x):
            return [
                e_z * np.abs(e_x - np.median([e for i, e in enumerate(x) if i != j]))
                for (j, e_x), e_z in zip(enumerate(x), stats.zscore(x))
            ]

        z_score_x = modif_zscore(x)
        z_score_x = vsigmoid(z_score_x)
        z_score_x = vnorm(z_score_x)

        return pd.Series(z_score_x, index=x.index)

    @staticmethod
    def dext_up_down(r, stat):
        try:
            return stat(r.dropna())
        except:
            pass

    def build_file_dext_score(self, pext_file, refseq_path, biomart_path, output_path):
        print("# Build dext file : {}".format(output_path))

        if os.path.isfile(output_path) is False:
            stats_l = list()

            print("# Files don't exist ☒")

            # REFSEQ treatments
            refseq = pd.read_parquet(refseq_path)
            refseq = refseq.loc[(refseq["mRNA_nb_total"] > 1) & (refseq["CDS_count"] > 1) & (refseq["Ratio_num"] < 1)]
            stats_l.append(
                {
                    "Filter": "RAW",
                    "Genes nb": refseq.Gene.nunique(),
                    "Transcripts nb": refseq.explode("mRNA_exons")["mRNA_exons"].nunique(),
                    "cExons nb": refseq.MAP.nunique(),
                }
            )

            # BIOMART treatments
            biomart = pd.read_csv(biomart_path, compression="gzip", sep="\t")

            biomart = (
                biomart[["Gene stable ID", "HGNC ID", "Gene name", "Chromosome/scaffold name"]]
                .rename(
                    {"Gene stable ID": "ENSG", "HGNC ID": "HGNC", "Gene name": "Gene", "Chromosome/scaffold name": "CHROM"},
                    axis=1,
                )
                .drop_duplicates()
            )
            # biomart["HGNC"] = biomart["HGNC"].str.replace("HGNC:", "")

            # MERGE IDS
            refseq = pd.merge(refseq, biomart, on="Gene")
            refseq["CHROM"] = refseq["CHROM"].astype(str)
            refseq = refseq.loc[~refseq["CHROM"].str.contains("HSCHR|HG")]
            stats_l.append(
                {
                    "Filter": "Merge BIOMART",
                    "Genes nb": refseq.Gene.nunique(),
                    "Transcripts nb": refseq.explode("mRNA_exons")["mRNA_exons"].nunique(),
                    "cExons nb": refseq.MAP.nunique(),
                }
            )

            # print(refseq)

            ## pext treatments
            pext_refseq = pd.merge(refseq, pext_file, on="MAP")
            pext_refseq = pext_refseq.rename(dicts["convert_tissue_dict"], axis=1)
            stats_l.append(
                {
                    "Filter": "Merge pext",
                    "Genes nb": pext_refseq.Gene.nunique(),
                    "Transcripts nb": pext_refseq.explode("mRNA_exons")["mRNA_exons"].nunique(),
                    "cExons nb": pext_refseq.MAP.nunique(),
                }
            )

            # DROP DUPLICATES
            pext_refseq = pext_refseq.drop_duplicates(subset=["MAP"]).sort_values(by="MAP").reset_index(drop=True)

            # UNEXPRESSED AND FULLY EXPRESSED EXONS
            pext_refseq = pext_refseq.loc[(pext_refseq["mean_proportion"] < 1)].reset_index(drop=True)
            stats_l.append(
                {
                    "Filter": "Filter out unexpressed & fully expressed",
                    "Genes nb": pext_refseq.Gene.nunique(),
                    "Transcripts nb": pext_refseq.explode("mRNA_exons")["mRNA_exons"].nunique(),
                    "cExons nb": pext_refseq.MAP.nunique(),
                }
            )

            # HANDLE NAN VALUES
            pext_refseq = pext_refseq.dropna(subset=dicts["GTEx_tissues_list"], how="all")
            stats_l.append(
                {
                    "Filter": "Filter out NaN in all tissues",
                    "Genes nb": pext_refseq.Gene.nunique(),
                    "Transcripts nb": pext_refseq.explode("mRNA_exons")["mRNA_exons"].nunique(),
                    "cExons nb": pext_refseq.MAP.nunique(),
                }
            )

            print(pd.DataFrame(stats_l))
            exit()

            # BUILD MATRIX WITH Z-SCORE + SIGMOID
            # z_df = pext_refseq[dicts["GTEx_tissues_list"]].head(100).progress_apply(lambda r: self.main_transform(r), axis=1) # ! TMP
            z_df = pext_refseq[dicts["GTEx_tissues_list"]].progress_apply(lambda r: self.main_transform(r), axis=1)

            # COMPUTE MIN & MAX dext COLS
            z_df_processed = z_df.copy()

            for stat in [min, max]:
                str_stat = stat.__name__

                if str_stat == "min":
                    dext_stat_convert = "down"
                elif str_stat == "max":
                    dext_stat_convert = "up"

                z_df_processed["dext_{}".format(dext_stat_convert)] = z_df.apply(lambda r: self.dext_up_down(r, stat), axis=1)

            concat_df = pd.concat(
                [
                    pext_refseq[[c for c in pext_refseq.columns if c not in dicts["GTEx_tissues_list"]]].reset_index(drop=True),
                    z_df_processed.reset_index(drop=True),
                ],
                axis=1,
            )
            stats_l.append(
                {
                    "Filter": "dext",
                    "Genes nb": concat_df.Gene.nunique(),
                    "Transcripts nb": concat_df.explode("mRNA_exons")["mRNA_exons"].nunique(),
                    "cExons nb": concat_df.MAP.nunique(),
                }
            )
            print(concat_df)

            # concat_df.to_parquet(output_path)
            stats_dext = pd.DataFrame(stats_l)
            print(stats_dext)
            stats_dext.to_excel(yaml["4_DEXT"]["Figures_data"]["dext_stats_test"])
            return concat_df

        else:
            print("# Files exist ✓, Loading ... ")
            concat_df = pd.read_parquet(output_path)

        return concat_df

    @staticmethod
    def compare_to_utrs(dext, utrs_path, exons_path, cds_path):
        # print(dext)
        utrs = pd.read_parquet(utrs_path)
        # print(utrs)
        exons = pd.read_parquet(exons_path)
        exons["Gene"] = exons["Attributes"].apply(lambda r: [e.replace("gene=", "") for e in r.split(";") if "gene=" in e][0])
        exons["mRNA"] = exons["Attributes"].apply(lambda r: [e.replace("Parent=rna-", "") for e in r.split(";") if "Parent=rna" in e][0])
        exons["mRNA"] = exons["mRNA"].apply(lambda r: r.split(".")[0])
        exons = exons.sort_values(by=["Gene", "Start", "End"])
        exons["ranges"] = exons["Start"].astype(str) + "-" + exons["End"].astype(str)

        cds = pd.read_parquet(cds_path)
        cds["Gene"] = cds["Attributes"].apply(lambda r: [e.replace("gene=", "") for e in r.split(";") if "gene=" in e][0])
        cds["mRNA"] = cds["Attributes"].apply(lambda r: [e.replace("Parent=rna-", "") for e in r.split(";") if "Parent=rna" in e][0])
        cds["mRNA"] = cds["mRNA"].apply(lambda r: r.split(".")[0])
        cds = cds.sort_values(by=["Gene", "Start", "End"])
        cds["ranges"] = cds["Start"].astype(str) + "-" + cds["End"].astype(str)

        # print(exons)
        # print(cds)

        exons_gb = (
            exons.groupby(["Gene", "mRNA"])["ranges"]
            .apply(lambda r: sorted(list(r)))
            .reset_index()
            .rename({"ranges": "exons_ranges"}, axis=1)
        )
        cds_gb = (
            cds.groupby(["Gene", "mRNA"])["ranges"].apply(lambda r: sorted(list(r))).reset_index().rename({"ranges": "cds_ranges"}, axis=1)
        )

        # print(exons_gb)
        # print(cds_gb)

        merge = pd.merge(exons_gb, cds_gb, on=["Gene", "mRNA"])
        merge["utrs"] = merge.apply(lambda r: list(sorted(list(set(r["exons_ranges"]).difference(set(r["cds_ranges"]))))), axis=1)
        merge = merge.explode("utrs")

        # dext = dext.loc[(dext["dext_up"] > 0.5) | (dext["dext_down"] < -0.5)]

        dext = dext.explode("mRNA_exons")
        merge_dext = pd.merge(dext[["Gene", "mRNA_exons", "ranges"]].rename({"mRNA_exons": "mRNA"}, axis=1), merge, on=["Gene", "mRNA"])

        # merge_dext = merge_dext.dropna(subset=["utrs"])
        merge["utrs"] = merge["utrs"].astype(str)
        print(merge_dext)

        def categorize_utr(r):
            alt_start = r["ranges"].split("-")[0]
            alt_end = r["ranges"].split("-")[1]
            utr_start = r["utrs"].split("-")[0]
            utr_end = r["utrs"].split("-")[1]
            if utr_start == alt_start and utr_end > alt_end:
                r["utr_bound"] = "3_prime"
                r["3_prime"] = r["ranges"]
            elif utr_end == alt_end and utr_start < alt_start:
                r["utr_bound"] = "5_prime"
                r["5_prime"] = r["ranges"]
            else:
                r["utr_bound"] = "Not"
                r["core"] = r["ranges"]
            return r

        merge_dext = merge_dext.parallel_apply(categorize_utr, axis=1)
        print(merge_dext)
        print(merge_dext["Gene"].nunique())
        print(merge_dext["mRNA"].nunique())

        print(merge_dext[["Gene", "mRNA", "3_prime"]].drop_duplicates().dropna())
        print(merge_dext[["Gene", "mRNA", "5_prime"]].drop_duplicates().dropna())
        print(merge_dext[["Gene", "mRNA", "core"]].drop_duplicates().dropna())

        tmp_gb = merge_dext.groupby(["Gene", "mRNA"])["utr_bound"].value_counts()
        print(tmp_gb)
        print(tmp_gb.groupby("utr_bound").sum())
        tmp_gb.name = "Count"
        print(tmp_gb.reset_index().Gene.nunique())
        print(tmp_gb.reset_index().mRNA.nunique())

        # print(merge_dext["3_prime"].nunique())
        # print(merge_dext["5_prime"].nunique())
        # print(merge_dext["core"].nunique())
        exit()


if __name__ == "__main__":

    c = HandleExonExpression()
