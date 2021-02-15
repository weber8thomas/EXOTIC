# -*- coding: utf-8 -*-
import collections
import multiprocessing
import sys
import pandas as pd
import parmap

from src.refseq_build.refseq_get_coding_genes import RefSeqBuildFileCodingGenes
from src import config

pd.set_option('mode.chained_assignment', None)


class RefSeqCDSAnalysis(RefSeqBuildFileCodingGenes):
	def __init__(self, input_name, output_name):
		self.cpus = multiprocessing.cpu_count() - 4
		raw_df, gene_ids = self.process_df(input_name)
		m = multiprocessing.Manager()
		output_list = m.list()
		parmap.starmap(self.mp_genes_id, list(zip(gene_ids)), raw_df, output_list, pm_pbar=True, pm_processes=self.cpus)
		self.process_list_df(output_name, output_list)

	@staticmethod
	def process_df(filename):
		"""

		Args:
			filename:

		Returns:

		"""
		raw_df_input = pd.read_csv(filename, compression='gzip', sep='\t', low_memory=False)
		gene_ids_list = list(raw_df_input['GeneID'].unique())
		raw_df_input['GeneID'] = raw_df_input['GeneID'].astype('category')
		return raw_df_input, gene_ids_list

	@staticmethod
	def process_list_df(output_name, output_list):
		"""

		Args:
		 output_list:
			output_name:

		Returns:

		"""
		df = pd.concat(list(output_list))
		df.to_csv(output_name, compression='gzip', sep='\t', index=False, chunksize=10000)
		print(df)

	@staticmethod
	def getOverlap(a, b):
		"""

		Args:
			a:
			b:

		Returns:

		"""
		return max(0, min(a[1], b[1]) - max(a[0], b[0]))

	def mp_genes_id(self, gene, raw_df, output_list):
		"""

		Args:
			gene:
			raw_df:
			output_list:

		Returns:

		"""
		df_gene = raw_df.loc[raw_df['GeneID'] == gene]
		gene_row = df_gene.loc[df_gene['Element'] == 'gene']
		mrnas_id = list(df_gene.loc[df_gene['Element'] == 'mRNA']['transcript_id'].unique())
		gene_row['mRNA_IDS'] = str(mrnas_id)
		for elem in ['CDS']:
			tmp_list = raw_df.loc[
				(raw_df['GeneID'] == gene) & (raw_df['Element'] == elem), ['Start', 'End', 'Count', 'mRNA_nb',
				                                                           'Dbxref']]
			tmp_list['Count'] = tmp_list['Count'].astype(str) + '/' + df_gene['mRNA_nb'].astype(str)
			tmp_list = tmp_list.values
			cds = sorted(list(set([tuple(e)[:3] for e in tmp_list])))

		positions_list = list()
		for elem in cds:
			positions_list.append(elem[0])
			positions_list.append(elem[1])
		positions_list = list(sorted(set(positions_list)))

		dico_intersection = collections.defaultdict(list)
		for index, position in enumerate(positions_list):
			if index < len(positions_list) - 1:
				for elem in cds:
					if self.getOverlap([position, positions_list[index + 1]], [elem[0], elem[1]]) > 0:
						dico_intersection[(position, positions_list[index + 1])].append(elem)
		final_dico_intersection = collections.defaultdict(dict)
		list_cds_counter = list()

		for elem in dico_intersection:
			final_dico_intersection[elem]['Sources'] = dico_intersection[elem]
			list_cds_counter.extend(list((elem[0], elem[1])))

			new_ratio = []
			nb_mrna = cds[0][2].split('/')[1]
			[new_ratio.append(float(source[2].split('/')[0])) for source in final_dico_intersection[elem]['Sources']]
			new_ratio = str(sum(new_ratio)) + '/' + str(nb_mrna)
			final_dico_intersection[elem]['Ratio'] = new_ratio

		dict_cds_counter = dict(collections.Counter(list_cds_counter))
		for elem in final_dico_intersection:
			check = False
			for e in elem:
				if check is False:
					if dict_cds_counter[e] == 1:
						final_dico_intersection[elem]['CDS_representation'] = 'Unique'
						final_dico_intersection[elem]['Sharing_status'] = False
					if dict_cds_counter[e] > 1:
						final_dico_intersection[elem]['CDS_representation'] = 'Variable_region'
						final_dico_intersection[elem]['Sharing_status'] = True
						check = True

		gene_row['CDS_locations_corrected'] = str(dict(final_dico_intersection))
		output_list.append(gene_row)


if __name__ == '__main__':
	RefSeqCDSAnalysis(sys.argv[1], sys.argv[2])
