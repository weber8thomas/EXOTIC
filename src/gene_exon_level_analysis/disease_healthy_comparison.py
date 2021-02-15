import collections
import itertools
import os
import sys
import warnings

import matplotlib.pyplot as plt
import numpy as np
import parmap
import seaborn as sns
from scipy import stats

warnings.simplefilter(action='ignore', category=FutureWarning)
import pandas as pd

pd.options.mode.chained_assignment = None  # default='warn'
from pandarallel import pandarallel
pandarallel.initialize(progress_bar=True, nb_workers=1)
from tqdm import tqdm
from src.utils import utils
from src.gene_exon_level_analysis.refseq_stats import RefSeqStats
pd.set_option('display.float_format', lambda x: '%.3f' % x)


class DiseaseHealthyComparison(RefSeqStats):


	def __init__(self, disease_df_path, healthy_df_path, output_count, output_ks):

		self.config = utils.load_config_file()
		self.logger.info('Working on Disease & Healthy genes to compare patterns')
		disease_df = pd.read_csv(disease_df_path, compression='gzip', sep='\t')
		healthy_df = pd.read_csv(healthy_df_path, compression='gzip', sep='\t')

		COUNT_mRNA_CDS_RATIO_df = self.plots_and_histo_analysis(disease_df, healthy_df)

	def plots_and_histo_analysis(self, disease_df, healthy_df):

		# FILTERING REFSEQ DISEASE GENES & MULTI-ISOFORMS
		if os.path.isfile(self.config['DiseaseHealthyComparison']['COUNT_mRNA_CDS_RATIO']) is False:

			if os.path.isfile(self.config['DiseaseHealthyComparison']['KS_results']) is False:
				self.kolmogorov_qqplot(disease_df, healthy_df)
			mrna_list = list(sorted(pd.concat([disease_df['mRNA_nb'], healthy_df['mRNA_nb']]).unique()))
			cds_list = list(sorted(pd.concat([disease_df['CDS_nb'], healthy_df['CDS_nb']]).unique()))
			df_disease_mrna_cds_frequency, df_disease_mrna_cds_genes = self.matrix_iso_exons(disease_df, mrna_list, cds_list, 'DISEASE')
			df_healthy_mrna_cds_frequency, df_healthy_mrna_cds_genes = self.matrix_iso_exons(healthy_df, mrna_list, cds_list, 'HEALTHY')
			COUNT_mRNA_CDS_RATIO_df = pd.concat(
				[df_disease_mrna_cds_frequency, df_disease_mrna_cds_genes, df_healthy_mrna_cds_frequency,
				 df_healthy_mrna_cds_genes], axis=0)
			self.write_df(COUNT_mRNA_CDS_RATIO_df, self.config['DiseaseHealthyComparison']['COUNT_mRNA_CDS_RATIO'])
			self.histo_frequency(COUNT_mRNA_CDS_RATIO_df)
		else:
			COUNT_mRNA_CDS_RATIO_df = self.read_df(self.config['DiseaseHealthyComparison']['COUNT_mRNA_CDS_RATIO'])
		return COUNT_mRNA_CDS_RATIO_df

	# @profile
	def kolmogorov_qqplot(self, df1, df2, plot=False):
		self.logger.info('DHC - KS test & plots')

		df1['STATUS'] = 'DISEASE'
		df2['STATUS'] = 'HEALTHY'
		df_concat = pd.concat([df1, df2], axis=0)
		selected_cols = [c for c in list(df_concat.columns) if 'exon' not in c]
		df_concat = df_concat[selected_cols]
		combinations_cds_mrna = list(
			sorted(itertools.product(list(df_concat.mRNA_nb.unique()), list(df_concat.CDS_nb.unique()))))
		output_list = list()
		LOCAL_DIR_FUNCTION = 'COMBINATION_CDS_MRNA_EXON_POSITION'
		utils.mkdir("data/2_processed/" + "DiseaseHealthyComparison" + '/' + LOCAL_DIR_FUNCTION)
		for combi in tqdm(combinations_cds_mrna, desc='KS computing'):
			# for combi in tqdm(combinations_cds_mrna[:20]):
			tmp_combi_df = df_concat.loc[(df_concat['mRNA_nb'] == combi[0]) & (df_concat['CDS_nb'] == combi[1])]
			if len(list(tmp_combi_df.STATUS.unique())) == 2:
				tmp_dict_array_freq_norm = dict()
				tmp_dict_genes = dict()
				tmp_dict_genes_nb = dict()
				tmp_dict_sizes = dict()
				tmp_dict_constitutive_alternative = dict()

				for status, color, align in zip(['DISEASE', 'HEALTHY'], ["#c0392b", "#3498db"], ['left', 'mid']):
					tmp_combi_df_status = tmp_combi_df.loc[tmp_combi_df['STATUS'] == status]
					tmp_dict_sizes[status] = tmp_combi_df_status.filter(regex='total_size$', axis=1).agg(
						'sum').to_dict()
					tmp_dict_constitutive_alternative[status] = tmp_combi_df_status[
						['Count_CDS_alternative', 'Count_CDS_constitutive']].sum().to_dict()
					tmp_combi_df_status['List_ratio_CDS'] = tmp_combi_df_status.CDS_locations_corrected.apply(
						lambda r: [eval(v['Ratio']) for k, v in eval(r).items()])
					tmp_dict_genes[status] = tmp_combi_df_status['Name'].values.tolist()
					gene_nb = tmp_combi_df_status.shape[0]
					tmp_dict_genes_nb[status] = gene_nb
					list_ratio_cds = [sub_l for sub_l in list(tmp_combi_df_status.List_ratio_CDS.values)]
					tmp_dict_array_freq_norm[status] = [sum / gene_nb for sum in np.sum(list_ratio_cds, 0)]

				tmp_df_grid = pd.DataFrame.from_dict(tmp_dict_array_freq_norm)
				tmp_df_grid['CDS_nb'] = range(1, tmp_df_grid.shape[0] + 1)

				tmp_df_grid = pd.melt(tmp_df_grid, id_vars=['CDS_nb'], value_vars=['DISEASE', 'HEALTHY'],
				                      var_name='STATUS', value_name='VALUE')

				if plot is True:
					g = sns.FacetGrid(tmp_df_grid, hue='STATUS', col='CDS_nb', sharex=False, sharey=False,
					                  palette=["#c0392b", "#3498db"], hue_order=['DISEASE', 'HEALTHY'],
					                  margin_titles=False)
					g.map(plt.hist, "VALUE", bins=np.arange(0, 1.2, 0.05), rwidth=2, align='left', alpha=0.3)
					g.fig.subplots_adjust(top=0.8)
					g.fig.suptitle('mRNA nb : {} - CDS nb : {}'.format(str(int(combi[0])), str(combi[1])),
					               fontsize='x-large')
					g.add_legend()

					# plt.show()
					tmp_dir = self.config['DiseaseHealthyComparison']['FACETGRID'].replace('.png', '_' + str(int(combi[0])) + '_' + str(
						combi[1]) + '.png').split('/')
					tmp_dir.insert(1, LOCAL_DIR_FUNCTION)
					# plt.tight_layout()
					plt.savefig("/".join(tmp_dir), dpi=150)
					plt.close()

				ks = stats.ks_2samp(tmp_dict_array_freq_norm['DISEASE'], tmp_dict_array_freq_norm['HEALTHY'])
				output_list.append({'mRNA_nb': combi[0], 'CDS_nb': combi[1], 'Stat_KS': ks[0], 'pvalue_KS': ks[1],
				                    'DISEASE_mean': np.mean(tmp_dict_array_freq_norm['DISEASE']),
				                    'HEALTHY_mean': np.mean(tmp_dict_array_freq_norm['HEALTHY']),
				                    'DISEASE_median': np.median(tmp_dict_array_freq_norm['DISEASE']),
				                    'HEALTHY_median': np.median(tmp_dict_array_freq_norm['HEALTHY']),
				                    'Compare_medians': float(np.median(tmp_dict_array_freq_norm['DISEASE']) - np.median(tmp_dict_array_freq_norm['HEALTHY'])) ** 2,
				                    'DISEASE_CDS_Unique_total_size': tmp_dict_sizes['DISEASE']['CDS_Unique_total_size'],
				                    'DISEASE_CDS_Variable_region_total_size': tmp_dict_sizes['DISEASE']['CDS_Variable_region_total_size'],
				                    'HEALTHY_CDS_Unique_total_size': tmp_dict_sizes['HEALTHY']['CDS_Unique_total_size'],
				                    'HEALTHY_CDS_Variable_region_total_size': tmp_dict_sizes['HEALTHY']['CDS_Variable_region_total_size'],
				                    'DISEASE_CDS_alternative': tmp_dict_constitutive_alternative['DISEASE']['Count_CDS_alternative'],
				                    'DISEASE_CDS_constitutive': tmp_dict_constitutive_alternative['DISEASE']['Count_CDS_constitutive'],
				                    'HEALTHY_CDS_alternative': tmp_dict_constitutive_alternative['HEALTHY']['Count_CDS_alternative'],
				                    'HEALTHY_CDS_constitutive': tmp_dict_constitutive_alternative['HEALTHY']['Count_CDS_constitutive'],
				                    'DISEASE_CDS_total': tmp_dict_constitutive_alternative['DISEASE']['Count_CDS_alternative'] + tmp_dict_constitutive_alternative['DISEASE']['Count_CDS_constitutive'],
				                    'HEALTHY_CDS_total': tmp_dict_constitutive_alternative['HEALTHY']['Count_CDS_alternative'] + tmp_dict_constitutive_alternative['HEALTHY']['Count_CDS_constitutive'],
				                    'DISEASE_gene_nb': tmp_dict_genes_nb['DISEASE'],
				                    'HEALTHY_gene_nb': tmp_dict_genes_nb['HEALTHY'],
				                    'DISEASE/HEALTHY CDS/Gene': float((tmp_dict_constitutive_alternative['DISEASE']['Count_CDS_alternative'] + tmp_dict_constitutive_alternative['DISEASE']['Count_CDS_constitutive']) / (tmp_dict_constitutive_alternative['HEALTHY']['Count_CDS_alternative'] + tmp_dict_constitutive_alternative['HEALTHY']['Count_CDS_constitutive'])) / float(tmp_dict_genes_nb['DISEASE'] / tmp_dict_genes_nb['HEALTHY']),
				                    'Overlapping_genes': len(set(tmp_dict_genes['DISEASE']).intersection(set(tmp_dict_genes['HEALTHY']))),
				                    'DISEASE_gene_list_and_functions': tmp_dict_genes['DISEASE'],
				                    'HEALTHY_gene_list_and_functions': tmp_dict_genes['HEALTHY'],
				                    # 'DISEASE_array': tmp_dict_array_freq_norm['DISEASE'],
				                    # 'HEALTHY_array': tmp_dict_array_freq_norm['HEALTHY']}
				                    })
		output_df = pd.DataFrame(list(output_list))
		output_df = output_df.sort_values(by='Stat_KS')
		output_df.to_excel(self.config['DiseaseHealthyComparison']['KS_results'], index=False)

	@staticmethod
	def matrix_iso_exons(df, mrna_list, cds_list, source):

		init_dict_cds_ratio = collections.defaultdict(dict)
		init_dict_genes = collections.defaultdict(dict)

		# df.mRNA_nb = df.mRNA_nb.astype('category')
		# df.CDS_nb = df.CDS_nb.astype('category')

		# mRNA count
		for mrna_count in tqdm(mrna_list, 'Building matrix for {} condition'.format(source)):
			cds_specific_count_mrna = df.loc[df['mRNA_nb'] == mrna_count]

			# CDS count
			for cds_nb in cds_list:
				tmp_l_ratio = list()
				cds_specific_count_mrna_count_cds = cds_specific_count_mrna.loc[
					cds_specific_count_mrna['CDS_nb'] == cds_nb]

				# GENE COUNT
				gene_list = list(cds_specific_count_mrna_count_cds['Name'].values)

				cds_locations_corrected = cds_specific_count_mrna_count_cds['CDS_locations_corrected']
				for cds_dict in cds_locations_corrected:
					cds_dict = eval(cds_dict)
					tmp_l_ratio += [eval(cds_info['Ratio']) for cds_positions, cds_info in cds_dict.items()]
				init_dict_cds_ratio[mrna_count][cds_nb] = tmp_l_ratio
				init_dict_genes[mrna_count][cds_nb] = gene_list

		# RETURN DF
		return_df_cds = pd.DataFrame.from_dict(init_dict_cds_ratio, orient='index').reset_index()
		return_df_cds['TYPE'] = 'CDS'
		return_df_cds['SOURCE'] = source
		return_df_genes = pd.DataFrame.from_dict(init_dict_genes, orient='index').reset_index()
		return_df_genes['TYPE'] = 'GENES'
		return_df_genes['SOURCE'] = source
		return return_df_cds, return_df_genes

	def histo_frequency(self, stats_df):
		self.logger.info('DHC - building statistics for all combinations')


		stats_df.rename({'index': 'mRNA_count'}, axis=1, inplace=True)
		stats_df = pd.melt(stats_df, id_vars=['mRNA_count', 'TYPE', 'SOURCE'],
		                   value_vars=[e for e in list(stats_df.columns) if e not in ['mRNA_count', 'TYPE', 'SOURCE']],
		                   value_name='List_Frequency_CDS', var_name='CDS_count')
		stats_df['Frequency'], stats_df['Counter'] = zip(
			*stats_df['List_Frequency_CDS'].parallel_apply(lambda r: self.apply_counter_cds_frequency(r)))
		stats_df.drop('List_Frequency_CDS', axis=1, inplace=True)
		stats_df = self.split_dataframe_rows(stats_df, ['Frequency', 'Counter'])
		stats_df['mRNA_count'] = stats_df['mRNA_count'].astype(float)
		stats_df['mRNA_count'] = stats_df['mRNA_count'].astype(int)
		stats_df['CDS_count'] = stats_df['CDS_count'].astype(float)
		stats_df['CDS_count'] = stats_df['CDS_count'].astype(int)
		stats_df['COMBINATION'] = stats_df['mRNA_count'].astype(str) + '_' + stats_df['CDS_count'].astype(str)
		# stats_df = stats_df.loc[(stats_df['COMBINATION'].isin(list(stats_df['COMBINATION'].unique())[:20]))]
		for combi in list(stats_df.COMBINATION.unique()):
			for source in list(stats_df.SOURCE.unique()):
				genes_nb = stats_df.loc[(stats_df['TYPE'] == 'GENES') & (stats_df['COMBINATION'] == combi) & (stats_df['SOURCE'] == source)].shape[0]
				stats_df.loc[(stats_df['TYPE'] == 'CDS') & (stats_df['COMBINATION'] == combi) & (stats_df['SOURCE'] == source), 'GENE_NB'] = genes_nb
		stats_df = stats_df.loc[stats_df['TYPE'] == 'CDS']
		stats_df['Frequency'] = stats_df['Frequency'].astype(float)
		stats_df['Counter'] = stats_df['Counter'] / stats_df['GENE_NB']
		stats_df.rename({'CDS_count': 'CDS', 'mRNA_count': 'mRNA'}, axis=1, inplace=True)
		stats_df.Frequency = stats_df.Frequency.round(2)

		combinations = list(stats_df['COMBINATION'].unique())
		plt.rcParams["patch.force_edgecolor"] = True

		SUB_OUTPUT_FILES = dict()
		SUB_OUTPUT_FILES['FACETGRID'] = 'Combination_exon_position.png'
		CURRENT_OUTPUT_DIR = self.OUTPUT_DIR + '/3_DiseaseHealthyComparison'
		SUB_OUTPUT_FILES = utils.create_new_sub_dir_and_update_dict(CURRENT_OUTPUT_DIR, 'COMBINATION_CDS_MRNA_EXON_POSITION', SUB_OUTPUT_FILES)
		parmap.starmap(self.parallel_plot_hist, list(zip(combinations)), stats_df, SUB_OUTPUT_FILES, pm_pbar=True, pm_processes=30)

	def parallel_plot_hist(self, combi, stats_df, SUB_OUTPUT_FILES):
		# for j, combi in tqdm(enumerate(combinations)):

		tmp_plot_df = stats_df.loc[stats_df.COMBINATION == combi]
		if len(list(tmp_plot_df.SOURCE.unique())) == 2:
			for source, color, align in zip(['DISEASE', 'HEALTHY'], ["#c0392b", "#3498db"], ['left', 'mid']):
				tmp_plot_df_source = tmp_plot_df.loc[tmp_plot_df.SOURCE == source]
				tmp_distplot_dict = dict()
				for i, line in tmp_plot_df_source.iterrows():
					tmp_distplot_dict[line['Frequency']] = line['Counter']
				val, weight = zip(*[(k, v) for k, v in tmp_distplot_dict.items()])
				plt.hist(val, weights=weight, bins=np.arange(0, 1.2, 0.05), rwidth=2, align='left', color=color,
				         alpha=0.3, label=source)
			plt.legend()
			plt.title('mRNA : {} - CDS : {}'.format(combi.split('_')[0], combi.split('_')[1]))
			plt.tight_layout()
			plt.savefig(SUB_OUTPUT_FILES['FACETGRID'].replace('.png', combi + '.png'), dpi=150)
			plt.close()

	@staticmethod
	def apply_counter_cds_frequency(row):
		tmp_d = dict(collections.Counter(eval(row)))
		return list(tmp_d.keys()), list(tmp_d.values())

	@staticmethod
	def split_dataframe_rows(df, column_selectors):
		# we need to keep track of the ordering of the columns
		def _split_list_to_rows(row, row_accumulator, column_selector):
			split_rows = {}
			max_split = 0
			for column_selector in column_selectors:
				split_row = row[column_selector]
				split_rows[column_selector] = split_row
				if len(split_row) > max_split:
					max_split = len(split_row)

			for i in range(max_split):
				new_row = row.to_dict()
				for column_selector in column_selectors:
					try:
						new_row[column_selector] = split_rows[column_selector].pop(0)
					except IndexError:
						new_row[column_selector] = ''
				row_accumulator.append(new_row)

		new_rows = []
		df.apply(_split_list_to_rows, axis=1, args=(new_rows, column_selectors))
		new_df = pd.DataFrame(new_rows, columns=df.columns)
		return new_df


if __name__ == "__main__":
	DiseaseHealthyComparison(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
