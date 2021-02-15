import nltk
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
import pandarallel
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
import vaex

# Read the taxrank ontology
url = 'http://purl.obolibrary.org/obo/hp.obo'
graph = obonet.read_obo(url)

# # Mapping from term ID to name
id_to_name = {id_: data.get('name') for id_, data in graph.nodes(data=True)}
name_to_id = {data['name']: id_ for id_,
              data in graph.nodes(data=True) if 'name' in data}


def remove_common(df, common):
    index_list = list()
    for j, row in df.iterrows():
        for r in ['Specific_Const', 'Specific_Alt']:
            if row[r]:
                for e in row[r]:
                    if e in common:
                        index_list.append(j)

    return df.drop(index_list).reset_index(drop=True)


def compare_phenotypes_const_alt(const, alt, tmp_dict):



    # BUILDING STANDARD VENN DIAGRAM
    common_phenotypes = set(const).intersection(set(alt))
    specific_const = set(const).difference(set(alt))
    specific_alt = set(alt).difference(set(const))


    # BUILDING DICT NETWORK PATHS
    d_graph = collections.defaultdict(dict)
    for s, name in zip([common_phenotypes, specific_const, specific_alt], ['Common', 'Const', 'Alt']):
        for elem in s:
            paths = nx.all_simple_paths(graph, source=name_to_id[elem], target=name_to_id['All'])
            d_graph[name][elem] = [[id_to_name[node] for node in path][:-1] for path in paths]



    specific = list()
    remove_list = list()
    for const_alt in ['Const', 'Alt']:
        for const_pheno in d_graph[const_alt]:
            for const_path in d_graph[const_alt][const_pheno]:

                tmp_specific = list()
                for common_pheno in d_graph['Common']:

                    # COMPUTE ALL INTERSECTION BETWEEN LIST OF PATHS AND RETURN THE LENGTH OF BIGGEST INTERSECTION
                    i_max = max([len(sorted(set(const_path).intersection(set(common_path)), key=lambda x: const_path.index(x))) for common_path in d_graph['Common'][common_pheno]])

                    for common_path in d_graph['Common'][common_pheno]:
                        intersection = sorted(set(const_path).intersection(set(common_path)), key=lambda x: const_path.index(x))
                       
                        # print(const_alt)
                        # print(const_path)
                        # print(common_path)
                        # print(intersection)

                        if intersection:
                            
                            # SHOW ONLY MOST RELEVANT INTERSECTION (USEFUL FOR CHILD/CLOSE PATHS)
                            if len(intersection) == i_max:
                                child = False
                                close = False
                                spec = False
                                
                                # FILTER THE CASES WHERE SPECIFIC ARE HIGHER TERMS COMPARE TO COMMON
                                # if len(intersection) < len(const_path):
                               
                                # CHILD HPO (DIRECT DESCENDANT)
                                if len(intersection) == len(common_path) and len(intersection) > 3:
                                    child = True
                                    close = True
                                    spec = False
                                    remove_list.append(const_pheno)

                                # CHILD HPO (DIRECT ASCENDANT)
                                if len(intersection) == len(const_path) and len(intersection) > 3:
                                    child = True
                                    close = True
                                    spec = False
                                    remove_list.append(const_pheno)

                                # CLOSE HPO (MAX SEPARATION = BY 2 NODES)
                                if (len(const_path) - len(common_path)) <= 2 and len(intersection) >= len(const_path)-2 and len(intersection) > 4:
                                    child = False
                                    close = True
                                    spec = True
                                    remove_list.append(const_pheno)
             

    new_specific_const = [e for e in specific_const if e not in set(remove_list)]
    new_specific_alt = [e for e in specific_alt if e not in set(remove_list)]
    # print(common_phenotypes)
    # print(new_specific_const)
    # print(new_specific_alt)
    return_dict = {
        **tmp_dict,
        **{
            'Pheno_Specific_Const': sorted(new_specific_const),
            'Pheno_Specific_Const_nb': len(sorted(new_specific_const)),
            'Pheno_Specific_Alt': sorted(new_specific_alt),
            'Pheno_Specific_Alt_nb': len(sorted(new_specific_alt)),
            'Pheno_Common': sorted(common_phenotypes),
            'Pheno_Common_nb': len(sorted(common_phenotypes)),
        }
    }
    return return_dict



def detailed_hpo_variants_proportion(info_dict, df):
    new_dict_variants_details = collections.defaultdict(dict)
    
    for col in ['Pheno_Specific_Const', 'Pheno_Specific_Alt', 'Pheno_Common']:
        newcol = col + '_details'
        new_dict_variants_details[newcol]['List_variants'] = list()
        new_dict_variants_details[newcol + '_Total_variants'] = 0
        if info_dict[col]:
            for pheno in info_dict[col]:
                new_dict_variants_details[newcol][pheno] = dict()
                gene = info_dict['CCRS_Gene']
                positions = df.loc[(df['OrphaCode'] == info_dict['OrphaCode']) & (df['HPO_name'] == pheno), 'POS']
                variants = df.loc[(df['CCRS_Gene'] == gene) & (df['POS'].isin(positions)), ['VAR_ID', 'MC', 'Real_Status', 'RS_STARS', 'CLNREVSTAT', 'rs']].drop_duplicates(subset=['VAR_ID', 'rs']).reset_index(drop=True)
                variants = variants.to_dict('records')
                new_dict_variants_details[newcol][pheno]['Phenotypes_variants_details'] = variants
                new_dict_variants_details[newcol][pheno]['Variants_nb'] = len(variants)
                for elem in new_dict_variants_details[newcol][pheno]['Phenotypes_variants_details']:
                    new_dict_variants_details[newcol]['List_variants'].append(elem['VAR_ID'])
            new_dict_variants_details[newcol]['List_variants'] = set(new_dict_variants_details[newcol]['List_variants'])
            new_dict_variants_details[newcol + '_Total_variants'] = len(new_dict_variants_details[newcol]['List_variants'])

    return new_dict_variants_details



# CONFIG
config = utils.load_config_file(project + "src/config.yaml")

pathogenic_genes_file = project + config['ExtractDiseaseGenes']['DISEASE_GENES_PATH']
refseq = project + config['RefSeq']['df_transposed_gene_file']
tmp_file = project + 'data/2_processed/Merge_HPO_OMIM_ORPHA_ClinVar.parquet'


merge = pd.read_parquet(tmp_file)
# hpo = pd.read_parquet(project + 'data/2_processed/HPO_ClinVar.parquet')



i_alt = 0
alt = list()
const = list()
i = 0
i_const = 0
genes_output = list()


merge = merge.replace('Abnormality of the foot', 'Abnormal foot morphology')
merge = merge.replace('Abnormality of abdomen morphology', 'Abnormal abdomen morphology')
merge = merge.replace('Focal seizures, afebril', 'Focal-onset seizure')

# for gene in tqdm(['TNNT2']):
for gene in tqdm(merge.CCRS_Gene.unique()):
    tmp_gene_df = merge.loc[merge['CCRS_Gene'] == gene].dropna(subset=['OrphaCode'])
    diseases = tmp_gene_df.OrphaCode.unique()


    for disease in diseases:
        # if disease == 'ORPHA:75249':
        # print(disease)

        tmp_dict = dict()
        for col in ['CCRS_Gene', 'RefSeq_Chrom', 'Disorder_name', 'Disorder_type', 'Disorder_classification']:
            result = list(
                tmp_gene_df.loc[(tmp_gene_df['OrphaCode'] == disease), col].unique())
            if len(result) == 1:
                result = result[0]
            tmp_dict[col] = result
        tmp_dict['OrphaCode'] = disease
        
        exons_variants = tmp_gene_df.loc[(tmp_gene_df['OrphaCode'] == disease), ['RefSeq_ranges', 'VAR_ID']].drop_duplicates().groupby('RefSeq_ranges')['VAR_ID'].apply(list).to_dict()
        tmp_dict['Exons_&_variants'] = exons_variants

        set_const = set(tmp_gene_df.loc[(tmp_gene_df['OrphaCode'] == disease) & (tmp_gene_df['Const_Alt'] == 'Const'), 'HPO_name'].unique())
        set_const = {e for e in set_const if 'inheritance' not in e}
        set_alt = set(tmp_gene_df.loc[(tmp_gene_df['OrphaCode'] == disease) & (tmp_gene_df['Const_Alt'] == 'Alt'), 'HPO_name'].unique())
        set_alt = {e for e in set_alt if 'inheritance' not in e}

        if set_const and set_alt:
            compare_hpo_dict = compare_phenotypes_const_alt(set_const, set_alt, tmp_dict)
            detailed_dict = detailed_hpo_variants_proportion(compare_hpo_dict, tmp_gene_df)
            # pprint({**tmp_dict, **compare_hpo_dict, **{"Details" : detailed_dict}})
            final_dict = {**tmp_dict, **compare_hpo_dict, **detailed_dict}
            genes_output.append(final_dict)



output = pd.DataFrame(genes_output)
output.to_excel(project + 'data/2_processed/HPO_specific_genes_and_disease_const_alt_analysis_detailed.xlsx')
