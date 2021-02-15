import pandas as pd
import vaex
import sys, os
import collections
import numpy as np
from tqdm import tqdm

raw_file = '/gstock/biolo_datasets/variation/benchmark/Databases/GTEX/V7/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm.gct.gz'
tissue_file = pd.read_csv('https://storage.googleapis.com/gtex_analysis_v7/annotations/GTEx_v7_Annotations_SampleAttributesDS.txt', sep='\t')[['SAMPID', 'SMTS', 'SMTSD']]

multi_iso_genes = pd.read_csv('/home/weber/PycharmProjects/ExoCarto/data/2_processed/GRCh37_RefSeq_lite_hgnc.csv.gz', compression='gzip', sep='\t')
multi_iso_genes = multi_iso_genes.loc[multi_iso_genes['mRNA_nb'] > 1, 'Name'].values.tolist()



# tex_tpm_raw = pd.read_csv(raw_file, compression='gzip', sep='\t', nrows=100, skiprows=2)

for i, batch_df in tqdm(enumerate(pd.read_csv(raw_file, compression='gzip', sep='\t', skiprows=2, chunksize=1000))):
    gtex_tpm_raw = batch_df.loc[batch_df['Description'].isin(multi_iso_genes)]
    gtex_tpm_melt = pd.melt(gtex_tpm_raw, id_vars=['Name', 'Description'], value_vars=list(gtex_tpm_raw.columns)[3:], var_name='SAMPID', value_name='TPM_value')
    merge_df = pd.merge(gtex_tpm_melt, tissue_file, on='SAMPID')

    output_name = '/gstock/biolo_datasets/variation/benchmark/Databases/GTEX/V7/GTEx_TPM_melt.csv.gz'
    if i == 0:
        merge_df.to_csv(output_name, mode='w', compression='gzip', sep='\t')
    elif i > 0:
        merge_df.to_csv(output_name, mode='a', compression='gzip', sep='\t', header=False)
