import pandas as pd
import vaex
import sys
import os
import collections
import numpy as np
from tqdm import tqdm

pivot_file = '/gstock/biolo_datasets/variation/benchmark/Databases/GTEX/V7/GTEx_TPM_melt.csv.gz'

pivot_df = pd.read_csv(pivot_file, compression='gzip', sep='\t')
pivot_df = pivot_df.drop(pivot_df.columns[0], axis=1)

stats_stmsd = pivot_df.groupby(['Name', 'Description', 'SMTSD']).describe()

new_cols = ['_'.join(e) for e in stats_stmsd.columns]
stats_stmsd.columns = new_cols
# stats_stmsd.to_csv('/gstock/biolo_datasets/variation/benchmark/Databases/GTEX/GTEx_TPM_STMSD_stats.csv.gz', compression='gzip', sep='\t')
stats_stmsd.to_parquet('/gstock/biolo_datasets/variation/benchmark/Databases/GTEX/V7/GTEx_TPM_STMSD_stats.parquet')
