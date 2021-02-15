import os

import sys
project = "/home/weber/PycharmProjects/ExoCarto/"
sys.path.insert(0, project)



import seaborn as sns
import numpy as np
import collections
import multiprocessing

import yaml
from cyvcf2 import VCF
import matplotlib.pyplot as plt
from pandarallel import pandarallel
from parmap import parmap
import _pickle
from src.utils import utils
pandarallel.initialize(progress_bar=True, nb_workers=multiprocessing.cpu_count())
from pprint import pprint
from tqdm import tqdm
import warnings

warnings.simplefilter(action='ignore', category=FutureWarning)
import pandas as pd

pd.options.mode.chained_assignment = None  # default='warn'
pd.set_option('display.float_format', lambda x: '%.3f' % x)
tqdm.pandas()


class VCF_analysis():

	def __init__(self, vcf_file_pathogenic, vcf_file_benign, refseq_transformed, output_file_visu, output_file, pathogenic_genes_file):
		# VEP ANNOTATION FIELD & SEPARATOR
		self.vep_field = 'CSQ'
		self.vep_separator = '|'
		self.short = False

		# MAP VARIANTS TO BED EXON FILE
		# disease_genes = [int(gene) for gene in _pickle.load(open(pathogenic_genes_file, 'rb')) if gene != '']

		if os.path.isfile(output_file) is False:
			
			bed_x_vcf = self.read_pandas(output_file, full=True)
			# bed_x_vcf = self.read_pandas(output_file, full=False, nrows=10000)
			# bed_x_vcf = bed_x_vcf.loc[bed_x_vcf['RefSeq_HGNC'].isin(disease_genes)]

		elif os.path.isfile(output_file) is True:
			# PREPARE VARIANTS DICTIONNARY & INDEX DICT FOR VEP FIELD
			d_variants, self.index_dict = self.prepare_dict_variants_pathogenic(vcf_file_pathogenic, )
			print(len(d_variants))
			# _pickle.dump(list(d_variants.keys()), open(pathogenic_genes_file, 'wb'))

			# d_variants, self.index_dict = self.prepare_dict_variants_benign(vcf_file_benign, d_variants)
			pprint(collections.Counter([variant_content['Source'] for gene, gene_content in d_variants.items() for variant, variant_content in gene_content.items()]))
			pprint(collections.Counter([variant_content['STATUS'] for gene, gene_content in d_variants.items() for variant, variant_content in gene_content.items()]))
			# INSTANCIATION OF EMPTY LIST TO NOT USE ALREADY MAPPED VARIANTS
			# self.check = list()
			bed_x_vcf = self.bed_x_vcf_compute(refseq_transformed, d_variants, output_file)
			exit()
			# print(bed_x_vcf.dropna(subset=['VAR_ID']))

		bed_x_vcf['HUE'] = bed_x_vcf['STATUS'] + '_' + bed_x_vcf['RefSeq_CDS_representation']
		bed_x_vcf.drop('Variant_mapping_list', axis=1, inplace=True)
		for condition in ['row', 'col', 'full']:
		# for condition in ['full']:
			bed_x_vcf_copy = bed_x_vcf.copy()
			self.prepare_data_for_visualization(bed_x_vcf_copy, output_file_visu, condition=condition)

	@staticmethod
	def prepare_header(vcf, vep_field, vep_separator):

		index_dict = dict()
		if vep_field:
			for h in vcf.header_iter():
				try:
					if h.info()['ID'] == vep_field:
						csq_header = h.info()['Description'].split(vep_separator)
						for elem in csq_header:
							index_dict[elem] = csq_header.index(elem)
				except:
					pass
		return index_dict

	@staticmethod
	def read_pandas(file_path, full=True, nrows=1000):
		if full is True:
			return pd.read_csv(file_path, compression='gzip', sep='\t', low_memory=False)
		if full is False:
			return pd.read_csv(file_path, compression='gzip', sep='\t', low_memory=False, nrows=nrows)

	def prepare_dict_variants_benign(self, vcf_file, d_variants):
		vep_consequences = {
			'HIGH': 3,
			'MODERATE': 2,
			'LOW': 1,
			'MODIFIER': 0,
		}
		changing_cols = ['splice_donor', 'splice_acceptor', 'missense', 'start', 'stop', 'intron', ]

		vcf = VCF(vcf_file)

		index_dict = self.prepare_header(vcf, vep_field=self.vep_field, vep_separator=self.vep_separator)


		for counter, variant in enumerate(tqdm(vcf)):
			if counter == 10000 and self.short is True:
				break

			if len(variant.REF) == 1 and len(variant.ALT[0]) == 1:
				csq = variant.INFO.get(self.vep_field)
				new_csq = list()
				if ',' in csq:
					csq = csq.split(',')
				else:
					csq = [csq]
				max_impact = max(
					[vep_consequences[case.split(self.vep_separator)[index_dict['IMPACT']]] for case in csq])
				if max_impact > 0:
					check_case = False

					for case in csq:
						case = case.split(self.vep_separator)

						if case[index_dict['HGNC_ID']]:
							hgnc = case[index_dict['HGNC_ID']]
							hgnc = int(hgnc)
							if hgnc in d_variants:
								if vep_consequences[case[index_dict['IMPACT']]] > 0:

									for consequence in changing_cols:
										if consequence in case[index_dict['Consequence']]:

											if check_case is False:
												d_variants[int(hgnc)][str(variant.POS) + '_' + str(variant.REF) + '_' + str(variant.ALT[0])] = {'VAR_ID' : str(variant.POS) + '_' + str(variant.REF) + '_' + str(variant.ALT[0]), 'Source' : 'gnomAD', 'STATUS': 'Benign', 'MC': consequence}
												check_case = True

						# ADD DICT TO LIST WITH ALL VARIANT INFO

		return d_variants, index_dict

	# FUNCTION TO RETRIEVE PATHOGENIC WELL REVIEWED VARIANTS FROM CLINVAR & HGMD
	def prepare_dict_variants_pathogenic(self, vcf_file, ):
		vep_consequences = {
			'HIGH': 3,
			'MODERATE': 2,
			'LOW': 1,
			'MODIFIER': 0,
		}

		vcf = VCF(vcf_file)

		# INDEX DICT
		index_dict = self.prepare_header(vcf, vep_field=self.vep_field, vep_separator=self.vep_separator)

		# CLNSIG FOR EACH DB
		db_clnsig = {
			'HGMD': 'DM',
			'ClinVar': 'athogenic',
		}

		# INSTANCIATION
		d_variants = collections.defaultdict(dict)
		i = 0

		# LOOP ON VARIANTS
		for counter, variant in enumerate(tqdm(vcf)):

			# BOOLEAN CHECKER FOR REVIEWING STATUS
			check = False

			# UNCOMMENT IF TEST
			if counter == 50000 and self.short is True:
				break

			# SNV ONLY
			if len(variant.REF) == 1 and len(variant.ALT[0]) == 1:

				# CLINVAR CLNSIG
				if 'CLNSIG=' in str(variant):
					cln_sig = variant.INFO.get('CLNSIG')
					db = 'ClinVar'

				# HGMD CLASS
				elif 'CLASS=' in str(variant):
					cln_sig = variant.INFO.get('CLASS')
					db = 'HGMD'

					# ALL HGMD VARIANTS ARE OK FOR REVIEWING STATUS

				try:
					# IF CLNSIG IS PRESENT
					if cln_sig:

						# IF CLNSIG OK FOR SELECTED DB
						if db == 'HGMD':
							if db_clnsig[db] == cln_sig:
								check = True
								# check = False
						elif db == 'ClinVar':
							if db_clnsig[db] in cln_sig and 'onflict' not in cln_sig:

								# # CLNREVSTAT
								# tmp_stats = variant.INFO['CLNREVSTAT'].split(',')
								# good_criteria = ['criteria_provided', 'reviewed_by_expert_panel',
								# 				'_multiple_submitters', 'practice_guideline']
								# bad_criteria = '_conflicting_interpretations'
								# if set(good_criteria) & set(tmp_stats) and bad_criteria not in tmp_stats:
								check = True
			

					# GOOD REVIEW STATUS
					if check is True:


						# VEP CSQ FIELD
						csq = variant.INFO.get(self.vep_field)
						if ',' in csq:
							csq = csq.split(',')
						else:
							csq = [csq]

						# RETRIEVE MAX IMPACT FOR ALL VEP CASES (EXAMPLE : ONE VARIANT OVERLAP ON TWO GENES)
						max_impact = max(
							[vep_consequences[case.split(self.vep_separator)[index_dict['IMPACT']]] for case in
								csq])
						# print(max_impact)
						# print([case.split(self.vep_separator)[index_dict['IMPACT']]  for case in csq])
						# print([case.split(self.vep_separator)[index_dict['Consequence']]  for case in csq])

						# LOOP ON CASES
						for case in csq:
							case = case.split(self.vep_separator)

							# IF CASE == MAX IMPACT
							if vep_consequences[case[index_dict['IMPACT']]] == max_impact and max_impact > 0:
								hgnc = case[index_dict['HGNC_ID']]

								# IF HGNC IS PRESENT
								# print(hgnc)
								# print(case[index_dict['Consequence']])
								if hgnc:
									hgnc = int(hgnc)

									for consequence in ['splice_donor', 'splice_acceptor', 'missense', 'start', 'stop', 'intron']:
										if consequence in case[index_dict['Consequence']]:
								# 			# ADD DICT TO LIST WITH ALL VARIANT INFO
											d_variants[int(hgnc)][str(variant.POS) + '_' + str(variant.REF) + '_' + str(variant.ALT[0])] = {'VAR_ID' : str(variant.POS) + '_' + str(variant.REF) + '_' + str(variant.ALT[0]), 'Source': db, 'STATUS': "Pathogenic", 'MC' : consequence}
					
				except KeyError:
					print(variant)
					pass

		return d_variants, index_dict

	def bed_x_vcf_compute(self, bed_file, d_variants, output_file):
		rename_cols = {
			"Gene": "CCRS_Gene",
			"Chrom": "RefSeq_Chrom",
			"Start": "RefSeq_Start",
			"End": "RefSeq_End",
			"ranges": "RefSeq_ranges",
			"Length": "RefSeq_Length",
			"GeneID": "RefSeq_GeneID",
			"HGNC": "RefSeq_HGNC",
			"CDS_representation": "RefSeq_CDS_representation",
			"Ratio": "RefSeq_Ratio",
			"Sharing_status": "RefSeq_Sharing_status",
			"Sources": "RefSeq_Sources",
			"Count_CDS_alternative": "RefSeq_Count_CDS_alternative",
			"Count_CDS_constitutive": "RefSeq_Count_CDS_constitutive",
			"mRNA_IDS": "RefSeq_mRNA_IDS",
			"mRNA_nb": "RefSeq_mRNA_nb",
		}

		# MP INSTANCIATION
		m = multiprocessing.Manager()
		l_df = m.list()
		# l_df = list()

		# READ REFSEQ FILE
		bed = self.read_pandas(bed_file, full=True)
		# bed = self.read_pandas(bed_file, full=False, nrows=100000)
		# bed['CCRS_ranges'] = bed['CCRS_ranges'].str.split(',')
		# bed = bed.explode('CCRS_ranges')
		# bed[['CCRS_Start', 'CCRS_End']] = bed['CCRS_ranges'].str.split('-', expand=True)
		# bed['CCRS_Start'] = bed['CCRS_Start'].astype(int)
		# bed['CCRS_End'] = bed['CCRS_End'].astype(int)
		bed.reset_index(drop=True, inplace=True)
		# bed = self.read_pandas(bed_file, full=False, nrows=10000)

		# RENAMING
		bed.rename(rename_cols, axis=1, inplace=True)

		# FILTERING
		bed = bed.loc[bed['RefSeq_mRNA_nb'] > 1]

		# PROCESSING HGNC
		bed = bed.loc[bed['RefSeq_HGNC'].isna() == False]
		bed.RefSeq_HGNC = bed.RefSeq_HGNC.astype(np.int64)

		# RefSeq exon freq
		bed['RefSeq_Ratio'] = bed['RefSeq_Ratio'].apply(lambda r: eval(r))

		# VCF PARSING
		bed_genes = set(bed.RefSeq_HGNC.values.tolist())
		vcf_genes = set(d_variants.keys())
		intersection_genes = bed_genes.intersection(vcf_genes)
		print('BED genes nb : {}'.format(str(len(bed_genes))))
		print('VCF genes nb : {}'.format(str(len(vcf_genes))))
		print('Intersection genes nb : {}'.format(str(len(intersection_genes))))
		# import matplotlib_venn
		# plt.figure()
		# matplotlib_venn.venn2([bed_genes, vcf_genes], ['BED', 'VCF'])
		# plt.savefig('visualization/VENN_BED_VCF.png')
		# plt.close()
		if intersection_genes:
			parmap.starmap(self.mp_process_gene, list(zip(intersection_genes)), bed, d_variants, l_df, pm_pbar=True, pm_processes=multiprocessing.cpu_count())
			# l_df = self.mp_process_gene(intersection_genes, bed, d_variants, l_df)
			bed_x_vcf = pd.concat(list(l_df)).sort_values(by='CCRS_Gene').drop(0, axis=1).reset_index(drop=True)
			bed_x_vcf.drop('Variant_mapping_list', axis=1, inplace=True)
			# bed_x_vcf.to_csv(output_file, compression='gzip', sep='\t', index=False)
			bed_x_vcf.to_parquet(output_file)
			return bed_x_vcf

	# FROM FLOAT TO BINS (EX : [0.803, 0.402, 0.657, 0.743] => [0.8-1, 0.4-0.6, 0.6-0.8, 0.6-0.8]
	@staticmethod
	def binning_ratio(df, condition):
		if condition == 'row':
			bins = [0, 1]
		else:
			bins = np.arange(0, 1.1, 0.2)
		labels = bins.copy()
		labels = [str(round(labels[j], 1)) + ' - ' + str(round(labels[j + 1], 1)) for j in range(len(labels) - 1)]
		df['RefSeq_Ratio'] = pd.cut(df['RefSeq_Ratio'], bins=bins, labels=labels)
		print(len(list(df['RefSeq_ranges'].unique())))
		c = dict(collections.Counter(
			df[['RefSeq_Ratio', 'RefSeq_ranges']].drop_duplicates(subset='RefSeq_ranges', keep='first')[
				'RefSeq_Ratio'].values.tolist()))
		for l in labels:
			if l not in c:
				c[l] = 0
		new_labels = [l + ' ({})'.format(c[l]) for l in labels]
		for l, n in zip(labels, new_labels):
			df['RefSeq_Ratio'] = df['RefSeq_Ratio'].str.replace(l, n)
		#
		df['RefSeq_Ratio'] = df['RefSeq_Ratio'].astype(str)
		# new_labels = labels
		return df, new_labels

	# SAME FOR LENGTH
	@staticmethod
	def binning_length(df, condition):
		df['RefSeq_Length'] = df['RefSeq_End'] - df['RefSeq_Start']
		dict_bins = df['RefSeq_Length'].describe().to_dict()
		tmp_bins = list()
		if condition == 'col':
			tmp_bins = [dict_bins['min'], dict_bins['max']]
		else:
			tmp_bins = [dict_bins['min'], dict_bins['25%'], dict_bins['50%'], dict_bins['75%'], dict_bins['max']]
		bins = list()
		for i, b in enumerate(tmp_bins):
			if i == 0:
				previous_elem = b
				bins.append(b)
			elif i != 0:
				curent_elem = b
				if curent_elem != previous_elem:
					bins.append(b)
				previous_elem = b

		labels = bins.copy()
		labels = [str(round(labels[j], 1)) + ' - ' + str(round(labels[j + 1], 1)) for j in range(len(labels) - 1)]
		df['RefSeq_Length_bins'] = pd.cut(df['RefSeq_Length'], bins=bins, labels=labels)
		print(len(list(df['RefSeq_ranges'].unique())))
		c = dict(collections.Counter(
			df[['RefSeq_Length_bins', 'RefSeq_ranges']].drop_duplicates(subset='RefSeq_ranges', keep='first')[
				'RefSeq_Length_bins'].values.tolist()))
		for l in labels:
			if l not in c:
				c[l] = 0
		new_labels = [l + ' ({})'.format(c[l]) for l in labels]
		for l, n in zip(labels, new_labels):
			df['RefSeq_Length_bins'] = df['RefSeq_Length_bins'].str.replace(l, n)
		# new_labels = labels
		return df, new_labels

	# MP FUNCTION TO PROCESS EACH GENE
	# @utils.profile
	def mp_process_gene(self, gene, ccrs, d_variants, l_ccrs_df):
		
		# for gene in tqdm(genes):
		bed_intersection_vcf_gene = ccrs.loc[ccrs['RefSeq_HGNC'] == gene]
		bed_intersection_vcf_gene.reset_index(drop=True, inplace=True)
		# bed_intersection_vcf_gene[['Variant_mapping_list', 'Variant_count', 'splice_donor', 'splice_acceptor', 'missense', 'start', 'stop', 'intron']] = bed_intersection_vcf_gene.apply(lambda r: self.apply_variants(r, d_variants[gene]), axis=1)
		bed_intersection_vcf_gene['Variant_mapping_list'] = self.apply_variants(bed_intersection_vcf_gene['RefSeq_Start'], bed_intersection_vcf_gene['RefSeq_End'], d_variants[gene])
		# bed_intersection_vcf_gene['Variant_mapping_list'] = self.apply_variants(bed_intersection_vcf_gene['CCRS_Start'].values, bed_intersection_vcf_gene['CCRS_End'].values, d_variants[gene])
		# bed_intersection_vcf_gene['Variant_mapping_list'] = bed_intersection_vcf_gene['Variant_mapping_list'].astype(str)
		# print(bed_intersection_vcf_gene.loc[bed_intersection_vcf_gene['Variant_mapping_list'] != '[]', ['CCRS_ranges', 'Variant_mapping_list']].to_string())
		# exit()
		bed_intersection_vcf_gene = bed_intersection_vcf_gene.explode('Variant_mapping_list')
		# print(bed_intersection_vcf_gene[['CCRS_ranges', 'Variant_mapping_list']])
		# exit()
		bed_intersection_vcf_gene = pd.concat([bed_intersection_vcf_gene, bed_intersection_vcf_gene['Variant_mapping_list'].apply(pd.Series)], axis=1)
		l_ccrs_df.append(bed_intersection_vcf_gene)

		# return l_ccrs_df

	# APPLY FUNCTION TO ADD NEW COLUMNS - NUMBER OF VARIANTS FOR EACH CATEGORY
	# @utils.profile
	def apply_variants(self, refseq_start, refseq_end, tmp):


		# CHAMPS INFOS SELECTIONNES
		selected_infos = ['VAR_ID', 'Source', 'STATUS', 'MC']

		def append_variant(refseq_start_pos, refseq_end_pos):
			tmp_list = list()

			# ON REGARDE TOUT LES VARIANTS POUR VOIR CEUX DEJ
			for variant, variant_content in tmp.items():
				pos = int(variant.split('_')[0])

				# ON VERIFIE QUE LE VARIANT EST DANS L'INTERVAL
				if refseq_start_pos <= pos < refseq_end_pos:

					# SI OK, ON L'AJOUTE A LA LISTE POUR L'ITERATION SUIVANTE
					variant_content = {k: v for k, v in variant_content.items() if k in selected_infos}

					# TMP FOR CSQ
					tmp_list.append(variant_content)

			return tmp_list
		
		l = list(map(append_variant, refseq_start, refseq_end))

		return pd.Series(l)


	def prepare_data_for_visualization(self, df, output_file_visu, condition):
		# df['Norm'] = (1000 * df['Variant_count']) / df['RefSeq_Length']
		# df = df.sort_values(by='Norm', ascending=False)
		# df.to_csv('data/2_processed/KB_ratio.csv.gz', compression='gzip', sep='\t', index=False)
		# print(df)



		# BINNING RATIO & LENGTH
		data_bins_ratio, order_ratio = self.binning_ratio(df, condition)
		data_bins_length, order_length = self.binning_length(data_bins_ratio, condition)

		print(df[['RefSeq_Length_bins', 'RefSeq_Ratio', 'HUE', 'RefSeq_Length']].groupby(['RefSeq_Length_bins', 'RefSeq_Ratio', 'HUE', ]).sum())

		# print(data_bins_length.drop(['VAR_ID', 'RefSeq_mRNA_IDS', 'RefSeq_Sources'],axis=1)
		# print(df[['RefSeq_Length_bins', 'RefSeq_Ratio', 'HUE', 'RefSeq_Length']].groupby(['RefSeq_Length_bins', 'RefSeq_Ratio', 'HUE',]).sum().reset_index())
		# exit()

		count_mc = df[['RefSeq_Length_bins', 'RefSeq_Ratio', 'HUE', 'MC', 'VAR_ID']].groupby(['RefSeq_Length_bins', 'RefSeq_Ratio', 'HUE', 'MC']).count().reset_index()
		sum_length = df[['RefSeq_Length_bins', 'RefSeq_Ratio', 'HUE', 'RefSeq_Length']].groupby(['RefSeq_Length_bins', 'RefSeq_Ratio', 'HUE', ]).sum().reset_index()
		final_data = pd.merge(count_mc, sum_length, on=['RefSeq_Length_bins', 'RefSeq_Ratio', 'HUE', ])

		final_data.rename({'MC': 'variable', 'VAR_ID': 'value'}, axis=1, inplace=True)

		final_data['RefSeq_Length'] = final_data['RefSeq_Length'] / 1000
		final_data['value_kb'] = final_data['value'] / final_data['RefSeq_Length']

		# COPY DF
		data_barplot = final_data.copy()

		# SUM ALL TYPE OF VARIANTS
		# data_barplot['all'] = data_barplot[changing_cols].sum(axis=1)

		# ADD TO LIST
		# changing_cols.append('all')

		# DROP STANDARD LENGTH COL

		# CONVERT TYPE TO CATEGORY TO REORDER FOR BOTH RATIO & LENGTH - SORT COLS
		data_barplot['RefSeq_Ratio'] = data_barplot['RefSeq_Ratio'].astype('category')
		data_barplot['RefSeq_Ratio'].cat.set_categories(order_ratio)
		data_barplot['RefSeq_Length_bins'] = data_barplot['RefSeq_Length_bins'].astype('category')
		data_barplot['RefSeq_Length_bins'].cat.set_categories(order_length)
		data_barplot = data_barplot.sort_values(by=['RefSeq_Ratio', 'RefSeq_Length_bins'])

		# RENAME COLS FOR PLOT
		pivot_table = data_barplot.copy()

		data_barplot.rename({'RefSeq_Length_bins': 'Length', 'RefSeq_Ratio': 'Freq'}, axis=1, inplace=True)
		data_barplot.to_excel(output_file_visu.replace('.xlsx', '_' + condition + '.xlsx'), index=False)

		variant_pivot = pivot_table.pivot_table(index=['RefSeq_Length_bins', 'RefSeq_Ratio', 'HUE'], columns='variable', values='value').reset_index()
		# variant_pivot[['splice_donor', 'splice_acceptor', 'missense', 'start', 'stop', 'intron', ]].fillna(0, inplace=True)
		variant_pivot[['missense', 'start', 'stop', 'intron', ]].fillna(0, inplace=True)
		pivot_length = df[['RefSeq_Length_bins', 'RefSeq_Ratio', 'HUE', 'RefSeq_Length']].groupby(['RefSeq_Length_bins', 'RefSeq_Ratio', 'HUE', ]).sum().reset_index()
		merge_pivot = pd.merge(variant_pivot, pivot_length, on=['RefSeq_Length_bins', 'RefSeq_Ratio', 'HUE'])
		merge_pivot['HUE'] = merge_pivot['HUE'].str.replace('Variable_region', 'Variable-region')
		merge_pivot[['STATUS', 'CDS TYPE']] = merge_pivot['HUE'].str.split('_', expand=True)
		merge_pivot['RefSeq_Length'] = merge_pivot['RefSeq_Length'] / 1000
		merge_pivot.sort_values(by=['RefSeq_Length_bins', 'RefSeq_Ratio', 'HUE'], inplace=True)
		merge_pivot = merge_pivot.drop('HUE', axis=1)
		merge_pivot.to_excel(output_file_visu.replace('.xlsx', '_raw_' + condition + '.xlsx'), index=False)


		data_barplot = data_barplot.drop(['RefSeq_Length'], axis=1)
		# return data_barplot
		for sharey in [True, False]:
			self.visualization_and_stats(data_barplot, order_ratio, order_length, condition, sharey)

	# VISUALIZATION FUNCTION
	def visualization_and_stats(self, data_barplot, order_ratio, order_length, condition, sharey):
		# CATEGORY OF PATHOGENIC STUDIED VARIANTS
		# changing_cols = ['splice_donor', 'splice_acceptor', 'missense', 'start', 'stop', 'intron', ]
		changing_cols = ['missense', 'start', 'stop', 'intron', ]

		# PLOT IN FACETGRID - https://seaborn.pydata.org/generated/seaborn.FacetGrid.html
		sns.set_style('whitegrid')  # SEABORN STYLE
		plt.figure()  # START FIGURE

		# START EMPTY FACETGRID - SELECT COL, ROW, ORDERS & OPTIONS
		g = sns.FacetGrid(data_barplot,
		                  col='Length',
		                  row='Freq',
		                  # hue='RefSeq_CDS_representation',
		                  margin_titles=True,
		                  sharey=sharey,
		                  row_order=order_ratio,
		                  col_order=order_length,
		                  )

		# MAP FOR EACH CELL & PRODUCE A BARPLOT - https://seaborn.pydata.org/generated/seaborn.barplot.html
		g.map(sns.barplot,
		      'variable',  # X axis
		      'value_kb',  # Y AXIS
		      'HUE',  # HUE
		      order=[c for c in changing_cols if c in list(data_barplot.variable.unique())],  # ORDER OF X AXIS
		      hue_order=list(sorted(list(data_barplot.HUE.unique()))),  # HUE ORDER
		      ci=None,  # NO CONFIDENCE INTERVAL (BLACK BAR)
		      palette=['#16a085', '#1abc9c', '#c0392b', '#e74c3c'],  # COLORS
		      )

		# # ROTATE X LABELS
		[plt.setp(ax.get_xticklabels(), rotation=90) for ax in g.axes.flat]

		# # X & Y AXIS LABELS SET TO NONE
		[plt.setp(ax.set_xlabel('')) for ax in g.axes.flat]
		[plt.setp(ax.set_ylabel('')) for ax in g.axes.flat]

		# # GRID AXIS
		[ax.yaxis.grid(True) for ax in g.axes.flat]
		[ax.xaxis.grid(True) for ax in g.axes.flat]

		# plt.legend()

		# OPTIMIZATION OF DISPLAY WITH MATPLOTLIB
		plt.tight_layout()

		# SAVE FIGURE
		plt.savefig("visualization/VCF_facetgrid_{}_{}_TP53.png".format(str(condition), str(sharey)))

		# CLOSE PLOT
		plt.close()




if __name__ == '__main__':
	config = utils.load_config_file()
	VCF_analysis(
		config['GeneExonLevel']['combine_path'],
		config['GeneExonLevel']['gnomad_pathogenic_genes'],
		config['CCRS_analysis']['RefSeq_transformed'],
		config['VCF_analysis']['Output_visu'],
		config['VCF_analysis']['Output'],
		config['ExtractDiseaseGenes']['DISEASE_GENES_PATH']
	)
# VCF_analysis(
# 	sys.argv[1],
# 	sys.argv[2],
# 	sys.argv[3],
# 	sys.argv[4],
# 	sys.argv[5],
# )
