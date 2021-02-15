import sys
project = "/home/weber/PycharmProjects/ExoCarto/"
sys.path.insert(0, project)
# IMPORTS
# import plotly.io as pio
# pio.renderers.default='notebook'
import plotly.express as px
import collections
import multiprocessing
import warnings
import plotly.graph_objects as go
import igv
import matplotlib.pyplot as plt
import parmap 

import multiprocessing
import numpy as np
import ptitprince as pt
import seaborn as sns
from tqdm import tqdm
from src.utils import utils
import _pickle
from pprint import pprint
warnings.simplefilter(action='ignore', category=FutureWarning)
import pandas as pd
pd.options.mode.chained_assignment = None  # default='warn'from pprint import pprint
import os
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from networkx.drawing.nx_agraph import graphviz_layout
import matplotlib.pyplot as plt
import networkx as nx
import obonet

# CONFIG
config = utils.load_config_file(project + "src/config.yaml")

pathogenic_genes_file = project + config['ExtractDiseaseGenes']['DISEASE_GENES_PATH']
refseq = project + config['RefSeq']['df_transposed_gene_file']

tmp_file = project + 'data/2_processed/Merge_HPO_OMIM_ORPHA_ClinVar.parquet'


lite_file = project + 'data/2_processed/refseq_vcf_pext_phylocsf.parquet'
vcf_df = pd.read_parquet(lite_file)

# CONST / ALT
vcf_df.loc[vcf_df['RefSeq_Ratio'] == 1, 'Const_Alt'] = 'Const'
vcf_df.loc[vcf_df['RefSeq_Ratio'] != 1, 'Const_Alt'] = 'Alt'


# PHENOTYPES
clinvar = pd.read_parquet(project + 'data/2_processed/' + os.path.basename(config['ClinVar']['09_20']).split('.')[0] + '_lite_table.parquet')
clinvar['CHROM'], clinvar['POS'], clinvar['REF'], clinvar['ALT'] = clinvar['VAR_ID'].str.split('_').str
clinvar = clinvar.loc[clinvar['RS_STARS'] > 0]
clinvar['MC'] = clinvar['MC'].apply(lambda r: r.split(',')[0])

extract_ccrs = vcf_df[['CCRS_Gene', 'RefSeq_HGNC', 'RefSeq_Chrom', 'RefSeq_ranges', 'RefSeq_Start', 'RefSeq_End', 'RefSeq_Ratio', 'Const_Alt', ]].drop_duplicates()
extract_ccrs = extract_ccrs.loc[extract_ccrs['RefSeq_Chrom'].str.contains('CHR') == False]

columns_to_process = ['POS', 'HPO', 'Status', 'RS_STARS', 'CLNREVSTAT', 'VAR_ID', 'MC', 'Real_Status', 'rs']

def select_gene(gene, l):
    ccrs_ranges = extract_ccrs.loc[extract_ccrs['CCRS_Gene'] == gene, 'RefSeq_ranges'].values
    instanciation_dict = collections.defaultdict(list)
    hpo_gene = clinvar.loc[clinvar['GENE'] == gene, columns_to_process]

    for ccrs_range in ccrs_ranges:
        for position in hpo_gene.values:
            if int(ccrs_range.split('-')[0]) <= int(position[0]) < int(ccrs_range.split('-')[1]):
                # print(gene, ccrs_range, position)
                instanciation_dict[ccrs_range].append(position.tolist())
    l.append(instanciation_dict)

genes = extract_ccrs.CCRS_Gene.unique()

m = multiprocessing.Manager()
l_hpo_dict = m.list()

parmap.starmap(select_gene, list(zip(genes)), l_hpo_dict, pm_pbar=True, pm_processes=30)


final_dict = {k: v for d in list(l_hpo_dict) for k, v in d.items()}
extract_ccrs['MAPPING'] = extract_ccrs['RefSeq_ranges'].map(final_dict)
extract_ccrs = extract_ccrs.dropna(subset=['MAPPING']).reset_index(drop=True)
extract_ccrs = extract_ccrs.explode('MAPPING').reset_index(drop=True)

for col, i in zip(columns_to_process, range(len(columns_to_process))):
    extract_ccrs[col] = extract_ccrs.MAPPING.apply(lambda r: r[i])
                            
extract_ccrs = extract_ccrs.drop('MAPPING', axis=1)
extract_ccrs = extract_ccrs.dropna(subset=['HPO'])
extract_ccrs.HPO = extract_ccrs.HPO.apply(lambda r: [e for e in r if 'HP:' in e])
extract_ccrs = extract_ccrs.explode('HPO')
extract_ccrs = extract_ccrs.dropna(subset=['HPO'])


extract_ccrs = extract_ccrs.loc[extract_ccrs['Status'] == 'Pathogenic'].reset_index(drop=True)

tmp = extract_ccrs[['CCRS_Gene', 'Const_Alt']].groupby('CCRS_Gene')['Const_Alt'].apply(list)
genes_const_alt = [gene for gene, const_alt in tmp.to_dict().items() if len(set(const_alt)) == 2]

extract_ccrs = extract_ccrs.loc[extract_ccrs['CCRS_Gene'].isin(genes_const_alt)]
extract_ccrs['MAP'] = extract_ccrs['CCRS_Gene'] + '_' + extract_ccrs['HPO']



genes_to_phenotypes = pd.read_csv('/gstock/biolo_datasets/variation/benchmark/Databases/HPO/genes_to_phenotype.txt', sep='\t', skiprows=1, names=['gene_id', 'gene_symbol', 'HPO', 'HPO_name', 'Frequency_raw', 'Frequency_HPO', 'Additional_info', 'Source', 'XRef'])
genes_to_phenotypes['MAP'] = genes_to_phenotypes['gene_symbol'] + '_' + genes_to_phenotypes['HPO']


def build_omim_orpha_column(row):
    orpha = np.nan
    omim = np.nan
    if 'OMIM:' in row:
        omim = row
    elif 'ORPHA:' in row:
        orpha = row
    return orpha, omim
    
genes_to_phenotypes['OrphaCode'], genes_to_phenotypes['OMIM'] = zip(*genes_to_phenotypes["XRef"].map(build_omim_orpha_column))

xref_orpha_omim = pd.read_csv('/home/weber/PycharmProjects/Orphanet/Disease_xref_orpha_omim.csv', sep='\t')
xref_orpha_omim = xref_orpha_omim[xref_orpha_omim.columns[1:]]
xref_orpha_omim['OrphaCode'] = xref_orpha_omim['OrphaCode'].astype(str)
xref_orpha_omim['OrphaCode'] = 'ORPHA:' + xref_orpha_omim['OrphaCode']

xref_orpha_omim['OMIM'] = xref_orpha_omim['OMIM'].fillna('[]')

xref_orpha_omim['OMIM'] = xref_orpha_omim['OMIM'].apply(eval)
xref_orpha_omim = xref_orpha_omim.explode('OMIM').reset_index(drop=True)
xref_orpha_omim['OMIM'] = xref_orpha_omim['OMIM'].astype(str)
xref_orpha_omim.loc[xref_orpha_omim['OMIM'].isna() == False, 'OMIM'] = 'OMIM:' + xref_orpha_omim['OMIM']


disease_db = pd.concat(
    [
        pd.merge(genes_to_phenotypes.loc[genes_to_phenotypes['OrphaCode'].isna() == False].drop('OMIM', axis=1), xref_orpha_omim, on='OrphaCode'),
        pd.merge(genes_to_phenotypes.loc[genes_to_phenotypes['OMIM'].isna() == False].drop('OrphaCode', axis=1), xref_orpha_omim, on='OMIM'),
    ]
).drop_duplicates()


pd.options.display.max_columns = 10
merge = pd.merge(extract_ccrs, disease_db, on='MAP', how='left').reset_index(drop=True).drop(['gene_symbol', 'gene_id'], axis=1)
print(merge)
merge.to_parquet(tmp_file)
