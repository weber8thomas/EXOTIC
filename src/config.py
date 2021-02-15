import os

CONFIG_DIR = dict(
	RefSeq=dict(
		raw_file="data/0_external/GRCh37_latest_genomic.gff.gz",
		build_file="data/1_interim/GRCh37_RefSeq_exon_frequency.csv.gz",
		df_transposed_gene_file="data/2_processed/GRCh37_RefSeq_lite.csv.gz",
	),

	GeneExonLevel=dict(
		clinvar_file="data/0_external/clinvar_vep_vcfanno.vcf.gz",
		hgmd_file="data/0_external/hgmd_vep_vcfanno.vcf.gz",
		ccrs_file="data/0_external/ccrs.autosomes.v2.20180420.bed.gz",
		biomart_file="data/0_external/BIOMART_CHR.txt.gz",
	),

	ExtractDiseaseGenes=dict(
		venn_clinvar_hgmd="visualization/Venn_diagram_pathogenic_genes_ClinVar_HGMD.png",
		disease_genes_pathogenic="data/1_interim/disease_genes.pkl",

	),

	RefSeqStats=dict(
		CORRECTED_DF_PATH="data/2_processed/GRCh37_corrected_refseq.csv.gz",
		CORRECTED_MULTI_ISO_PATH="data/2_processed/GRCh37_corrected_refseq_multi_iso.csv.gz",
		CORRECTED_DISEASE_PATH="data/2_processed/GRCh37_corrected_refseq_disease.csv.gz",
		CORRECTED_HEALTHY_PATH="data/2_processed/GRCh37_corrected_refseq_healthy.csv.gz",
		STATS_PATH="statistics/GRCh37_stats_total.xlsx",
		STATS_MULTI_ISO_PATH="statistics/GRCh37_stats_multi_iso.xlsx",
		STATS_DISEASE_PATH="statistics/GRCh37_stats_multi_iso_disease.xlsx",
		STATS_HEALTHY_PATH="statistics/GRCh37_stats_multi_iso_healthy.xlsx",
		mRNA_nb_counter="statistics/mRNA_stats_counter.xlsx",
		mRNA_nb_metrics="statistics/mRNA_stats_metrics.xlsx",
	),

	DiseaseHealthyComparison=dict(
		COUNT_mRNA_CDS_RATIO='data/2_processed/Count_mRNA_CDS_Ratio_disease_healty.csv.gz',
		KS_results="statistics/KS_results_CDS_mRNA.xlsx",
	),

	CCRS_analysis=dict(
		RefSeq_transformed='data/2_processed/RefSeq_transformed.csv.gz',
		RefSeq_CCRS="data/2_processed/RefSeq_CCRS.csv.gz",
	),
)

PROJECT_DIR = os.path.dirname(os.path.abspath("../../README.md")) + '/'