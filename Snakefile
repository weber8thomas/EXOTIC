import os
import sys
print(sys.version)
from src.utils import utils
#
PROJECT_DIR = os.path.dirname(os.path.abspath("../../README.md")) + "/"
configfile: "src/config.yaml"
# print(config["refseq"])

# EXPLOIT RESULTS
rule target:
	input:
		config["CCRS_analysis"]["RefSeq_CCRS"],
		config["CCRS_analysis"]["RefSeq_transformed"],
# #------------------------------------------------------------------------------

# BUILD REFSEQ INIT FILE
rule refseq_get_coding_genes:
	input:
		config["RefSeq"]["raw_file"]
	output:
		config["RefSeq"]["build_file"]
	shell:
		"python src/refseq_build/refseq_get_coding_genes.py {input} {output}"

# # TRANSFORM REFSEQ TO ONE ROW / GENE
rule refseq_transform_file:
	input:
		config["RefSeq"]["build_file"]
	output:
		config["RefSeq"]["df_transposed_gene_file"]
	shell:
		"python src/refseq_build/refseq_transform_file.py {input} {output}"

# # EXTRACT DISEASE & HEALTHY GENES FROM CLINVAR + HGMD
rule extract_disease_genes:
	input:
	     clinvar = config["GeneExonLevel"]["clinvar_file"],
	     hgmd = config["GeneExonLevel"]["hgmd_file"]
	output:
	     config["ExtractDiseaseGenes"]["DISEASE_GENES_PATH"]
	shell:
		"python src/gene_exon_level_analysis/extract_disease_genes.py {input.clinvar} {input.hgmd} {output}"

# PRODUCE REFSEQ STATS
rule refseq_stats:
	input:
	     refseq_gene_file = config["RefSeq"]["df_transposed_gene_file"],
	     disease_genes = config["ExtractDiseaseGenes"]["DISEASE_GENES_PATH"],
	output:
	     corrected_df = config["RefSeqStats"]["CORRECTED_DF_PATH"],
	     corrected_multi_iso_df = config["RefSeqStats"]["CORRECTED_MULTI_ISO_PATH"],
	     corrected_disease_df = config["RefSeqStats"]["CORRECTED_DISEASE_PATH"],
	     corrected_healthy_df = config["RefSeqStats"]["CORRECTED_HEALTHY_PATH"],
	shell:
		"python src/gene_exon_level_analysis/refseq_stats.py {input.refseq_gene_file} {input.disease_genes} {output.corrected_df} {output.corrected_multi_iso_df} {output.corrected_disease_df} {output.corrected_healthy_df}"

# # COMPARE DISEASE & HEALTHY COMPARISON ON GENE AND EXON LEVEL
# rule disease_healthy_comparison:
# 	input:
# 	     disease_df = config["RefSeqStats"]["CORRECTED_DISEASE_PATH"],
# 	     healthy_df = config["RefSeqStats"]["CORRECTED_HEALTHY_PATH"],
# 	output:
# 	     count = config["DiseaseHealthyComparison"]["COUNT_mRNA_CDS_RATIO"],
# 	     ks = config["DiseaseHealthyComparison"]["KS_results"],
# 	shell:
# 		"python src/gene_exon_level_analysis/disease_healthy_comparison.py {input.disease_df} {input.healthy_df} {output.count} {output.ks}"

# CCRS
rule ccrs_analysis:
	input:
	     multi_iso_refseq = config["RefSeqStats"]["CORRECTED_MULTI_ISO_PATH"],
	     biomart = config["GeneExonLevel"]["biomart_file"],
	     ccrs_file = config["GeneExonLevel"]["ccrs_file"],
	output:
	     refseq_x_ccrs = config["CCRS_analysis"]["RefSeq_CCRS"],
	     refseq_transformed = config["CCRS_analysis"]["RefSeq_transformed"],
	shell:
		"python src/nucleotide_level_analysis/ccrs_analysis.py {input.multi_iso_refseq} {input.biomart} {input.ccrs_file} {output.refseq_transformed} {output.refseq_x_ccrs}"


# VCF
# rule vcf_analysis:
# 	input:
# 		vcf_path = config["VCF_analysis"]["VCF_PATH_TEST"],
# 		vcf_benign = config["VCF_analysis"]["VCF_BENIGN_TEST"],
# 		bed = config["VCF_analysis"]["RefSeq_TEST"]
# 	output:
# 		output_visu = config["VCF_analysis"]["Output_visu"],
# 		output = config["VCF_analysis"]["Output"],
# 	shell:
# 		"python src/nucleotide_level_analysis/vcf_analysis.py {input.vcf_path} {input.vcf_benign} {input.bed} {output.output_visu} {output.output}"

