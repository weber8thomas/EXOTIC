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


def compare_phenotypes_const_alt(const, alt):

    list_compare_phenotypes = list()

    # LOOP OVER CONSTITUTIVE EXONS PHENOTYPES
    for const_elem in const:
        # LOOP OVER ALTERNATIVE EXONS PHENOTYPES
        for alt_elem in alt:

            # CASE N°1 : IDENTICAL HPO
            if const_elem == alt_elem:
                list_compare_phenotypes.append(
                    {
                        'Identical': True,
                        'Child': False,
                        'Term_Const': const_elem,
                        'Term_Alt': alt_elem,
                        # 'List_HPO_Const' :np.nan,
                        # 'List_HPO_Alt' :np.nan,
                        'Overlap': np.nan,
                        'Nb_of_overlapping_terms': np.nan,
                        'Last_common_term': np.nan,
                        'Deepest_term': np.nan,
                        'Specific_Const': [],
                        'Specific_Alt': [],
                    }
                )

            # TODO : common

            # CASE N°2 : DIFFERENT HPO
            elif const_elem != alt_elem:

                # BUILD PATHS FROM OBO GRAPH
                d = collections.OrderedDict()
                d['Const'] = list()
                d['Alt'] = list()
                for elem, const_alt in zip([const_elem, alt_elem], ['Const', 'Alt']):
                    try:
                        paths = nx.all_simple_paths(
                            graph,
                            source=name_to_id[elem],
                            target=name_to_id['All']
                        )
       
                        for path in paths:
                            d[const_alt].append([id_to_name[node] for node in path][:-1])
                    except KeyError:
                        print(elem)

                l_associations = list()
                last_elem_alt = ""
                # ITER OVER PATH ON CONST EXONS HPO
                for sub_list_const in d['Const']:

                    check = True

                    # ITER OVER PATH ON ALT EXONS HPO
                    for sub_list_alt in d['Alt']:

                        # CASE N°2A : INTERSECTION BETWEEN HPO PATHS
                        if len(set(sub_list_const).intersection(set(sub_list_alt))) == 1:

                            print(const_elem, alt_elem)
                            print(set(sub_list_const))
                            print(set(sub_list_alt))
                            print('\n')

                            if last_elem_alt == sub_list_alt[0]:
                                check = False

                            if check is True:

                                check = False
                                list_compare_phenotypes.append(
                                    {
                                        'Identical': False,
                                        'Child': False,
                                        'Term_Const': const_elem,
                                        'Term_Alt': alt_elem,
                                        # 'List_HPO_Const' : sub_list_const,
                                        # 'List_HPO_Alt' : sub_list_alt,
                                        'Overlap': np.nan,
                                        'Nb_of_overlapping_terms': np.nan,
                                        'Last_common_term': np.nan,
                                        'Deepest_term': np.nan,
                                        'Specific_Const': [const_elem],
                                        'Specific_Alt': [alt_elem],
                                    }
                                )

                        # CASE N°2B : INTERSECTION BETWEEN HPO PATHS
                        elif len(set(sub_list_const).intersection(set(sub_list_alt))) > 1:

                            overlap = sorted(set(sub_list_const).intersection(
                                sub_list_alt), key=lambda x: sub_list_const.index(x))
                            spec_const = sorted(set(sub_list_const).difference(
                                sub_list_alt), key=lambda x: sub_list_const.index(x))
                            spec_alt = sorted(set(sub_list_alt).difference(
                                sub_list_const), key=lambda x: sub_list_alt.index(x))

                            shortest = sub_list_const if len(sub_list_const) < len(sub_list_alt) else sub_list_alt
                            longest = sub_list_const if len(sub_list_const) > len(sub_list_alt) else sub_list_alt

                            ps = nltk.stem.PorterStemmer()
                            last_common_term = overlap[0]
                            last_common_term_nltk = ps.stem(last_common_term.lower()).split(' ')
                            deepest_term = longest[0]
                            deepest_term_nltk = ps.stem(deepest_term.lower()).split(' ')

                            child = True if [e for e in last_common_term_nltk if e in deepest_term_nltk] else False

                            if len(overlap) > 3:
                                if spec_const and spec_alt:
                                    child = False
                                l_associations.append(
                                    {
                                        'Identical': False,
                                        'Child': child,
                                        'Term_Const': const_elem,
                                        'Term_Alt': alt_elem,
                                        # 'List_HPO_Const' : sub_list_const,
                                        # 'List_HPO_Alt' : sub_list_alt,
                                        'Overlap': overlap,
                                        'Nb_of_overlapping_terms': len(overlap),
                                        'Last_common_term': last_common_term,
                                        'Deepest_term': deepest_term,
                                        'Specific_Const': spec_const,
                                        'Specific_Alt': spec_alt,
                                    }
                                )

                        last_elem_alt = alt_elem
                if l_associations:
                    max_level_overlap = max(
                        [d['Nb_of_overlapping_terms'] for d in l_associations])
                    l_associations = [d for d in l_associations if d['Nb_of_overlapping_terms'] == max_level_overlap]
                    list_compare_phenotypes.extend(l_associations)
    output_df = pd.DataFrame(list_compare_phenotypes)

    pd.options.display.max_rows = 1000

    specific_alt = output_df.Specific_Alt.values.tolist()
    specific_alt = ["-".join(e) if type(e) is list else e for e in specific_alt]
    specific_const = output_df.Specific_Const.values.tolist()
    specific_const = ["-".join(e) if type(e)
                    is list else e for e in specific_const]

    common = output_df.loc[output_df['Identical']
                        == True, 'Term_Const'].values.tolist()


    # output_df = remove_common(output_df, common)
    print(output_df.sort_values(by=['Term_Const', 'Term_Alt']).reset_index(drop=True))
    print(common)
    print(output_df.loc[output_df['Identical']== True])
    output_df = output_df.drop(output_df.loc[(output_df['Identical'] == False) & (output_df['Term_Const'].isin(common))].index)
    output_df = output_df.drop(output_df.loc[(output_df['Identical'] == False) & (output_df['Term_Alt'].isin(common))].index)
    print(output_df.sort_values(by=['Term_Const', 'Term_Alt']).reset_index(drop=True))
    exit()
    return output_df.sort_values(by=['Term_Const', 'Term_Alt']).reset_index(drop=True)




# CONFIG
config = utils.load_config_file(project + "src/config.yaml")

pathogenic_genes_file = project + config['ExtractDiseaseGenes']['DISEASE_GENES_PATH']
refseq = project + config['RefSeq']['df_transposed_gene_file']
tmp_file = project + 'data/2_processed/Merge_HPO_OMIM_ORPHA_ClinVar.parquet'


merge = pd.read_parquet(tmp_file)
hpo = pd.read_parquet(project + 'data/2_processed/HPO_ClinVar.parquet')

i_alt = 0
alt = list()
const = list()
i = 0
i_const = 0
genes_output = list()


for gene in tqdm(['BBS4']):
    tmp_gene_df = merge.loc[merge['CCRS_Gene']
                            == gene].dropna(subset=['OrphaCode'])
    diseases = tmp_gene_df.OrphaCode.unique()

    for disease in diseases:

        tmp_dict = dict()
        for col in ['CCRS_Gene', 'RefSeq_Chrom', 'Disorder_name', 'Disorder_type', 'Disorder_classification']:
            result = list(
                tmp_gene_df.loc[(tmp_gene_df['OrphaCode'] == disease), col].unique())
            if len(result) == 1:
                result = result[0]
            tmp_dict[col] = result
        cds_const = list(tmp_gene_df.loc[(tmp_gene_df['OrphaCode'] == disease) & (tmp_gene_df['Const_Alt'] == 'Const'), 'RefSeq_ranges'].unique())
        position_const = list(tmp_gene_df.loc[(tmp_gene_df['OrphaCode'] == disease) & (tmp_gene_df['Const_Alt'] == 'Const'), 'HPO_POSITION'].unique())
        variants_const = hpo.loc[(hpo['POSITION'].isin(position_const)) & (hpo['GENE'] == gene), 'VAR_ID'].unique().tolist()

        cds_alt = list(tmp_gene_df.loc[(tmp_gene_df['OrphaCode'] == disease) & (tmp_gene_df['Const_Alt'] == 'Alt'), 'RefSeq_ranges'].unique())
        position_alt = list(tmp_gene_df.loc[(tmp_gene_df['OrphaCode'] == disease) & (tmp_gene_df['Const_Alt'] == 'Alt'), 'HPO_POSITION'].unique())
        variants_alt = hpo.loc[(hpo['POSITION'].isin(position_alt)) & (hpo['GENE'] == gene), 'VAR_ID'].unique().tolist()

        tmp_dict['CDS_Const'] = cds_const
        tmp_dict['VARIANTS_Const'] = variants_const
        tmp_dict['CDS_Alt'] = cds_alt
        tmp_dict['VARIANTS_Alt'] = variants_alt

        set_const = set(tmp_gene_df.loc[(tmp_gene_df['OrphaCode'] == disease) & (tmp_gene_df['Const_Alt'] == 'Const'), 'HPO_name'].unique())
        set_alt = set(tmp_gene_df.loc[(tmp_gene_df['OrphaCode'] == disease) & (tmp_gene_df['Const_Alt'] == 'Alt'), 'HPO_name'].unique())

        if set_alt and set_const and len(set_alt) > 1 and len(set_const) > 1:
            # print(gene)
            i += 1
            pprint(tmp_dict)
            print(compare_phenotypes_const_alt(set_const, set_alt))



#             specific_const = set_const.difference(set_alt)

#             if specific_const:
#                 i_const += 1
#                 const.append(gene)

#             specific_alt = set_alt.difference(set_const)

#             if specific_alt:
#                 i_alt += 1
#                 alt.append(gene)
#             common = set_const.intersection(set_alt)
#             # print(gene, set_const, set_alt, specific_const, specific_alt)
#             output_dict = {
#                 **{
#                     'Gene': gene,
#                     'Disease': disease,
#                     'Specific_Const': specific_const,
#                     'Specific_Alt': specific_alt,
#                     'Common': common,
#                 }, **tmp_dict
#             }
#             genes_output.append(output_dict)
# output = pd.DataFrame(genes_output)
# # output.to_excel(project + 'data/2_processed/HPO_specific_genes_and_disease_const_alt_analysis.xlsx')
# print(i, len(set(const)), len(set(alt)))












# const = ['Polydactyly', 'Syndactyly', 'Brachydactyly', 'Renal cyst', 'Nyctalopia', 'Retinal degeneration',
#     'External genital hypoplasia', 'Abnormality of the dentition', 'Rod-cone dystrophy', 'Intellectual disability', 'Hypogonadism']
# alt = ['Multicystic kidney dysplasia', 'Postaxial hand polydactyly', 'Brachydactyly',  'Nystagmus', 'Obesity', 'Short neck', 'Short stature', 'Medial flaring of the eyebrow', 'Nephrotic syndrome', 'Neurological speech impairment',
#        'Hypoplasia of penis', 'Downslanted palpebral fissures', 'Cryptorchidism', 'Finger syndactyly', 'Hypertension', 'Low-set, posteriorly rotated ears', 'Hearing impairment', 'Prominent nasal bridge', 'Abnormal electroretinogram']

# const = ['Retinal degeneration']
# alt = ['Rod-cone dystrophy']

