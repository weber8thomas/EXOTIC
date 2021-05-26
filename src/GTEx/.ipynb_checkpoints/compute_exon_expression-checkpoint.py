import sys
project = "/home/weber/PycharmProjects/ExoCarto/"
sys.path.insert(0, project)

import collections
import multiprocessing
import warnings
import parmap 
import multiprocessing
import numpy as np
from tqdm import tqdm
from src.utils import utils
import _pickle
import pandarallel
from pprint import pprint
warnings.simplefilter(action='ignore', category=FutureWarning)
import pandas as pd
pd.options.mode.chained_assignment = None  # default='warn'from pprint import pprint
import os
# pandarallel.pandarallel.initialize(nb_workers=os.cpu_count(), progress_bar=True)

# CONFIG
config = utils.load_config_file(project + "src/config.yaml")

pathogenic_genes_file = project + config['ExtractDiseaseGenes']['DISEASE_GENES_PATH']
refseq = project + config['RefSeq']['df_transposed_gene_file']

rename_GTEX_dict = {
        'Adipose-Subcutaneous': 'GTEX-Adipose-Subcutaneous',
        'Adipose-Visceral(Omentum)': 'GTEX-Adipose-Visceral(Omentum)',
        'AdrenalGland': 'GTEX-AdrenalGland',
        'Artery-Aorta': 'GTEX-Artery-Aorta',
        'Artery-Coronary': 'GTEX-Artery-Coronary',
        'Artery-Tibial': 'GTEX-Artery-Tibial',
        'Bladder': 'GTEX-Bladder',
        'Brain-Amygdala': 'GTEX-Brain-Amygdala',
        'Brain-Anteriorcingulatecortex(BA24)': 'GTEX-Brain-Anteriorcingulatecortex(BA24)',
        'Brain-Caudate(basalganglia)': 'GTEX-Brain-Caudate(basalganglia)',
        'Brain-CerebellarHemisphere': 'GTEX-Brain-CerebellarHemisphere',
        'Brain-Cerebellum': 'GTEX-Brain-Cerebellum',
        'Brain-Cortex': 'GTEX-Brain-Cortex',
        'Brain-FrontalCortex(BA9)': 'GTEX-Brain-FrontalCortex(BA9)',
        'Brain-Hippocampus': 'GTEX-Brain-Hippocampus',
        'Brain-Hypothalamus': 'GTEX-Brain-Hypothalamus',
        'Brain-Nucleusaccumbens(basalganglia)': 'GTEX-Brain-Nucleusaccumbens(basalganglia)',
        'Brain-Putamen(basalganglia)': 'GTEX-Brain-Putamen(basalganglia)',
        'Brain-Spinalcord(cervicalc-1)': 'GTEX-Brain-Spinalcord(cervicalc-1)',
        'Brain-Substantianigra': 'GTEX-Brain-Substantianigra',
        'Breast-MammaryTissue': 'GTEX-Breast-MammaryTissue',
        'Cells-EBV-transformedlymphocytes': 'GTEX-Cells-EBV-transformedlymphocytes',
        'Cells-Transformedfibroblasts': 'GTEX-Cells-Transformedfibroblasts',
        'Cervix-Ectocervix': 'GTEX-Cervix-Ectocervix',
        'Cervix-Endocervix': 'GTEX-Cervix-Endocervix',
        'Colon-Sigmoid': 'GTEX-Colon-Sigmoid',
        'Colon-Transverse': 'GTEX-Colon-Transverse',
        'Esophagus-GastroesophagealJunction': 'GTEX-Esophagus-GastroesophagealJunction',
        'Esophagus-Mucosa': 'GTEX-Esophagus-Mucosa',
        'Esophagus-Muscularis': 'GTEX-Esophagus-Muscularis',
        'FallopianTube': 'GTEX-FallopianTube',
        'Heart-AtrialAppendage': 'GTEX-Heart-AtrialAppendage',
        'Heart-LeftVentricle': 'GTEX-Heart-LeftVentricle',
        'Kidney-Cortex': 'GTEX-Kidney-Cortex',
        'Liver': 'GTEX-Liver',
        'Lung': 'GTEX-Lung',
        'MinorSalivaryGland': 'GTEX-MinorSalivaryGland',
        'Muscle-Skeletal': 'GTEX-Muscle-Skeletal',
        'Nerve-Tibial': 'GTEX-Nerve-Tibial',
        'Ovary': 'GTEX-Ovary',
        'Pancreas': 'GTEX-Pancreas',
        'Pituitary': 'GTEX-Pituitary',
        'Prostate': 'GTEX-Prostate',
        'Skin-NotSunExposed(Suprapubic)': 'GTEX-Skin-NotSunExposed(Suprapubic)',
        'Skin-SunExposed(Lowerleg)': 'GTEX-Skin-SunExposed(Lowerleg)',
        'SmallIntestine-TerminalIleum': 'GTEX-SmallIntestine-TerminalIleum',
        'Spleen': 'GTEX-Spleen',
        'Stomach': 'GTEX-Stomach',
        'Testis': 'GTEX-Testis',
        'Thyroid': 'GTEX-Thyroid',
        'Uterus': 'GTEX-Uterus',
        'Vagina': 'GTEX-Vagina',
        'WholeBlood': 'GTEX-WholeBlood',
    }



if os.path.isfile(project + 'data/2_processed/temp_file.parquet') is False:
    mrna_df = pd.read_csv(
    project + 'data/1_interim/RefSeq_GRCh37_complete_new.csv.gz', compression='gzip', sep='\t')
    list_df = list()

    for j, row in tqdm(mrna_df.iterrows()):
        if row['Element'] == 'gene':
            gene = row['Name']
        elif row['Element'] == 'mRNA':
            mrna = row['Name'].split('.')[0]
        elif row['Element'] == 'CDS':
            exon = row['Elem_position_ID'].replace('_', '-')
            list_df.append(
                {
                    'Gene': gene,
                    'mRNA': mrna,
                    'Exon_Start': int(exon.split('-')[0]),
                    'Exon_End': int(exon.split('-')[1]),
                }
            )

    refseq_df = pd.DataFrame(list_df)
    # print(refseq_df.loc[refseq_df['Gene'] == ''])
    global_df = pd.read_parquet(project + 'data/2_processed/refseq_vcf_pext_phylocsf.parquet')
    lite_df = global_df[[
        'RefSeq_HGNC', 'CCRS_Gene', 'RefSeq_Start', 'RefSeq_End', 'RefSeq_mRNA_IDS',
    ]].drop_duplicates()
    lite_df.to_parquet(project + 'data/2_processed/temp_file.parquet')
else:
    lite_df = pd.read_parquet(project + 'data/2_processed/temp_file.parquet')



if os.path.isfile(config['GTEX']['tmp_gtex_mapping']) is False:
    lite_df = lite_df.drop('RefSeq_mRNA_IDS', axis=1)
    lite_df['RefSeq_Start'] = lite_df['RefSeq_Start'].astype(int)
    lite_df['RefSeq_End'] = lite_df['RefSeq_End'].astype(int)

    pd.options.display.max_rows = 100
    pd.options.display.max_colwidth = 100

    lite_df = pd.merge(refseq_df, lite_df, left_on=['Gene', 'Exon_Start', 'Exon_End'], right_on=[
                    'CCRS_Gene', 'RefSeq_Start', 'RefSeq_End'], how='right')
    print(lite_df)
    biomart_file = pd.read_csv(config['GTEX']['biomart_mapping_transcripts'], compression='gzip', sep='\t').dropna()
    biomart_file['HGNC ID'] = biomart_file['HGNC ID'].str.replace('HGNC:', '')
    biomart_file['HGNC ID'] = biomart_file['HGNC ID'].astype(int)
    biomart_file = biomart_file.drop('Transcript stable ID version', axis=1)

    # lite_df['RefSeq_mRNA_IDS'] = lite_df['RefSeq_mRNA_IDS'].apply(eval)
    # lite_df = lite_df.explode('RefSeq_mRNA_IDS')
    # lite_df['RefSeq_mRNA_IDS'] = lite_df['RefSeq_mRNA_IDS'].str.split('.').str[0]

    merge_df = pd.merge(lite_df, biomart_file, how='left', left_on=['RefSeq_HGNC', 'mRNA'], right_on=['HGNC ID', 'RefSeq mRNA ID'])
    merge_df = merge_df.drop(['HGNC ID', 'RefSeq mRNA ID', 'Gene stable ID version', 'Gene'], axis=1)
    merge_df = merge_df.dropna()
    merge_df = merge_df.reset_index(drop=True)
    merge_df.to_parquet(config['GTEX']['tmp_gtex_mapping'])

else:
    merge_df = pd.read_parquet(config['GTEX']['tmp_gtex_mapping'])


gtex_df = pd.read_csv(config['GTEX']['transcript_average_file'], compression='gzip', sep='\t')

merge_refseq_gtex = pd.merge(merge_df, gtex_df, left_on=['Gene stable ID', 'Transcript stable ID'], right_on=['gene_id', 'transcript_id']).drop(['Gene stable ID', 'Transcript stable ID'], axis=1)

merge_refseq_gtex = merge_refseq_gtex.groupby(['RefSeq_HGNC', 'CCRS_Gene', 'gene_id', 'RefSeq_Start', 'RefSeq_End']).agg({'transcript_id' : ",".join, 'mRNA' : ",".join, 'Adipose-Subcutaneous' : np.mean, 'Adipose-Visceral(Omentum)' : np.mean, 'AdrenalGland' : np.mean, 'Artery-Aorta' : np.mean, 'Artery-Coronary' : np.mean, 'Artery-Tibial' : np.mean, 'Bladder' : np.mean, 'Brain-Amygdala' : np.mean, 'Brain-Anteriorcingulatecortex(BA24)' : np.mean, 'Brain-Caudate(basalganglia)' : np.mean, 'Brain-CerebellarHemisphere' : np.mean, 'Brain-Cerebellum' : np.mean, 'Brain-Cortex' : np.mean, 'Brain-FrontalCortex(BA9)' : np.mean, 'Brain-Hippocampus' : np.mean, 'Brain-Hypothalamus' : np.mean, 'Brain-Nucleusaccumbens(basalganglia)' : np.mean, 'Brain-Putamen(basalganglia)' : np.mean, 'Brain-Spinalcord(cervicalc-1)' : np.mean, 'Brain-Substantianigra' : np.mean, 'Breast-MammaryTissue' : np.mean, 'Cells-EBV-transformedlymphocytes' : np.mean, 'Cells-Transformedfibroblasts' : np.mean, 'Cervix-Ectocervix' : np.mean, 'Cervix-Endocervix' : np.mean, 'Colon-Sigmoid' : np.mean, 'Colon-Transverse' : np.mean, 'Esophagus-GastroesophagealJunction' : np.mean, 'Esophagus-Mucosa' : np.mean, 'Esophagus-Muscularis' : np.mean, 'FallopianTube' : np.mean, 'Heart-AtrialAppendage' : np.mean, 'Heart-LeftVentricle' : np.mean, 'Kidney-Cortex' : np.mean, 'Liver' : np.mean, 'Lung' : np.mean, 'MinorSalivaryGland' : np.mean, 'Muscle-Skeletal' : np.mean, 'Nerve-Tibial' : np.mean, 'Ovary' : np.mean, 'Pancreas' : np.mean, 'Pituitary' : np.mean, 'Prostate' : np.mean, 'Skin-NotSunExposed(Suprapubic)' : np.mean, 'Skin-SunExposed(Lowerleg)' : np.mean, 'SmallIntestine-TerminalIleum' : np.mean, 'Spleen' : np.mean, 'Stomach' : np.mean, 'Testis' : np.mean, 'Thyroid' : np.mean, 'Uterus' : np.mean, 'Vagina' : np.mean, 'WholeBlood' : np.mean, }).reset_index()
merge_refseq_gtex = merge_refseq_gtex.rename(rename_GTEX_dict, axis=1)
merge_refseq_gtex['GTEX-mean_proportion'] = merge_refseq_gtex.filter(regex='GTEX').mean(axis=1)
# merge_refseq_gtex = merge_refseq_gtex.head(1000)
merge_refseq_gtex = pd.melt(merge_refseq_gtex, id_vars=['RefSeq_HGNC', 'CCRS_Gene', 'gene_id', 'RefSeq_Start', 'RefSeq_End', 'transcript_id', 'mRNA', ], value_vars=['GTEX-Adipose-Visceral(Omentum)', 'GTEX-AdrenalGland', 'GTEX-Artery-Aorta', 'GTEX-Artery-Coronary', 'GTEX-Artery-Tibial', 'GTEX-Bladder', 'GTEX-Brain-Amygdala', 'GTEX-Brain-Anteriorcingulatecortex(BA24)', 'GTEX-Brain-Caudate(basalganglia)', 'GTEX-Brain-CerebellarHemisphere', 'GTEX-Brain-Cerebellum', 'GTEX-Brain-Cortex', 'GTEX-Brain-FrontalCortex(BA9)', 'GTEX-Brain-Hippocampus', 'GTEX-Brain-Hypothalamus', 'GTEX-Brain-Nucleusaccumbens(basalganglia)', 'GTEX-Brain-Putamen(basalganglia)', 'GTEX-Brain-Spinalcord(cervicalc-1)', 'GTEX-Brain-Substantianigra', 'GTEX-Breast-MammaryTissue', 'GTEX-Cells-EBV-transformedlymphocytes', 'GTEX-Cells-Transformedfibroblasts', 'GTEX-Cervix-Ectocervix', 'GTEX-Cervix-Endocervix', 'GTEX-Colon-Sigmoid', 'GTEX-Colon-Transverse', 'GTEX-Esophagus-GastroesophagealJunction', 'GTEX-Esophagus-Mucosa', 'GTEX-Esophagus-Muscularis', 'GTEX-FallopianTube', 'GTEX-Heart-AtrialAppendage', 'GTEX-Heart-LeftVentricle', 'GTEX-Kidney-Cortex', 'GTEX-Liver', 'GTEX-Lung', 'GTEX-MinorSalivaryGland', 'GTEX-Muscle-Skeletal', 'GTEX-Nerve-Tibial', 'GTEX-Ovary', 'GTEX-Pancreas', 'GTEX-Pituitary', 'GTEX-Prostate', 'GTEX-Skin-NotSunExposed(Suprapubic)', 'GTEX-Skin-SunExposed(Lowerleg)', 'GTEX-SmallIntestine-TerminalIleum', 'GTEX-Spleen', 'GTEX-Stomach', 'GTEX-Testis', 'GTEX-Thyroid', 'GTEX-Uterus', 'GTEX-Vagina', 'GTEX-WholeBlood', 'GTEX-mean_proportion'])
# print(merge_refseq_gtex.loc[merge_refseq_gtex['CCRS_Gene'] == 'TTN', 'variable'].unique())
# print(sorted(merge_refseq_gtex.CCRS_Gene.unique()))
# exit()

median_expression_genes_tissus = pd.read_parquet(config['GTEX']['stats_gtex_median']).reset_index()
median_expression_genes_tissus['SMTSD'] = 'GTEX-' + median_expression_genes_tissus['SMTSD'].str.replace(' ', '')
print(median_expression_genes_tissus.loc[median_expression_genes_tissus['Description'] == 'TTN'])



# merge_gtex = pd.merge(merge_refseq_gtex, median_expression_genes_tissus[['Description', 'SMTSD', 'TPM_value_50%']], left_on=['CCRS_Gene', 'variable'], right_on=['Description', 'SMTSD'])
# merge_gtex['value_normalized'] = merge_gtex['value'] / merge_gtex['TPM_value_50%']
# pd.options.display.max_rows = 20000
# print(merge_gtex.loc[merge_gtex['CCRS_Gene'] == 'TTN'])



