import sys
project = "/home/weber/PycharmProjects/ExoCarto/"
sys.path.insert(0, project)

import warnings
import multiprocessing
import collections
import parmap
import os
import pandas as pd
from pprint import pprint
import pandarallel
import _pickle
from src.utils import utils
from tqdm import tqdm
import numpy as np
import warnings
import multiprocessing
import collections
import sys
import vaex

project = "/home/weber/PycharmProjects/ExoCarto/"
sys.path.insert(0, project)

warnings.simplefilter(action='ignore', category=FutureWarning)
# default='warn'from pprint import pprint
pd.options.mode.chained_assignment = None
# pandarallel.pandarallel.initialize(nb_workers=os.cpu_count(), progress_bar=True)

# CONFIG
config = utils.load_config_file(project + "src/config.yaml")

pathogenic_genes_file = project + \
    config['ExtractDiseaseGenes']['DISEASE_GENES_PATH']
refseq = project + config['RefSeq']['df_transposed_gene_file']

if os.path.isfile(project + 'data/2_processed/temp_file.parquet') is True:
    global_df = vaex.open(
        project + 'data/2_processed/refseq_vcf_pext_phylocsf.parquet')

    global_df = global_df[['RefSeq_HGNC', 'CCRS_Gene',
                           'RefSeq_Start', 'RefSeq_End', 'RefSeq_mRNA_IDS']]

    column_names = list(global_df.columns)
    sets = [global_df._set(column_name) for column_name in column_names]
    counts = [set.count for set in sets]
    set_names = [global_df.add_variable(f'set_{column_name}', set, unique=True) for column_name, set in zip(column_names, sets)]
    expression = global_df[f'_ordinal_values({column_names[0]}, {set_names[0]})'].astype(
        'int64')
    product_count = 1
    for count, set_name, column_name in zip(counts[:-1], set_names[1:], column_names[1:]):
        product_count *= count
        expression = expression + global_df[f'_ordinal_values({column_name}, {set_name})'].astype('int64') * product_count
    global_df['row_id'] = expression
    index = global_df._index('row_id')
    unique_row_ids = global_df.row_id.unique()
    indices = index.map_index(unique_row_ids)
    global_df = global_df.take(indices)
    global_df['RefSeq_mRNA_IDS'] = global_df['RefSeq_mRNA_IDS'].apply(eval)
    global_df = global_df.explode('RefSeq_mRNA_IDS')

    # lite_df = global_df[[
        # 'RefSeq_HGNC', 'CCRS_Gene', 'RefSeq_Start', 'RefSeq_End', 'RefSeq_mRNA_IDS',
    # ]].drop_duplicates()




    # lite_df.to_parquet(project + 'data/2_processed/temp_file.parquet')
else:
    lite_df = vaex.open(project + 'data/2_processed/temp_file.parquet')
# print(lite_df)
print(global_df)
# if os.path.isfile(config['GTEX']['tmp_gtex_mapping']) is False:

#     biomart_file = pd.read_csv(
#         config['GTEX']['biomart_mapping_transcripts'], compression='gzip', sep='\t').dropna()
#     biomart_file['HGNC ID'] = biomart_file['HGNC ID'].str.replace('HGNC:', '')
#     biomart_file['HGNC ID'] = biomart_file['HGNC ID'].astype(int)
#     biomart_file = biomart_file.drop('Transcript stable ID version', axis=1)

#     lite_df['RefSeq_mRNA_IDS'] = lite_df['RefSeq_mRNA_IDS'].apply(eval)
#     lite_df = lite_df.explode('RefSeq_mRNA_IDS')
#     lite_df['RefSeq_mRNA_IDS'] = lite_df['RefSeq_mRNA_IDS'].str.split(
#         '.').str[0]

#     merge_df = pd.merge(lite_df, biomart_file, how='left', left_on=[
#                         'RefSeq_HGNC', 'RefSeq_mRNA_IDS'], right_on=['HGNC ID', 'RefSeq mRNA ID'])
#     merge_df = merge_df.drop(
#         ['HGNC ID', 'RefSeq mRNA ID', 'Gene stable ID version'], axis=1)
#     merge_df = merge_df.dropna()
#     merge_df.to_parquet(config['GTEX']['tmp_gtex_mapping'])

# else:
#     merge_df = pd.read_parquet(config['GTEX']['tmp_gtex_mapping'])

# gtex_df = pd.read_csv(config['GTEX']['transcript_average_file'],
#                       compression='gzip', sep='\t', nrows=1000)
# print(gtex_df)
