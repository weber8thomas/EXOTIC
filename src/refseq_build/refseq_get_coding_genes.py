# -*- coding: utf-8 -*-
import collections
import gzip
import multiprocessing
import sys

import numpy as np
import pandas as pd
import parmap
from tqdm import tqdm

from src import config

pd.set_option('mode.chained_assignment', None)


class RefSeqBuildFileCodingGenes:

	def __init__(self, filename, output_name_path):
		"""

		Args:
			filename (str): refseq raw compressed file
			output_name_path (str): output name for the produced file
		"""
		self.cpus = multiprocessing.cpu_count() - 4

		self.output_name = output_name_path
		file_input = gzip.open(filename, 'rb')
		data = self.process_refseq_file(file_input)
		df = self.process_df(data)
		self.stats = self.get_coding_genes(df)
		self.return_stats()

	def process_refseq_file(self, file_input):
		"""
		Function to process GFF :
			- Remove header
			- Get each record to process_record function
		Args:
			file_input (file): input file name

		Returns:
			refseq_build (list): list of dictionaries
		"""
		data = list()
		# PROCESS EACH ROW OF REFSEQ FILE
		for j, record in enumerate(tqdm(file_input)):
			record = record.decode()
			if record.startswith('#') is False:
				record = record.strip().split('\t')

				# CALL SECONDARY FUNCTION TO DECCOMPOSE ROW
				output_dict = self.process_record(record)
				data.append(output_dict)
		# UNCOMMENT IF TEST
			if j == 100000:
				break
		return data

	@staticmethod
	def process_record(record):
		"""
		Function to process each line of GFF file
		Args:
			record (list): list corresponding to one line of GFF

		Returns:
			output_dict (dict): dictionary with all required info split correctly
		"""
		attributes = record[-1].split(';')
		tmp_dict_attributes = {elem.split('=')[0]: elem.split('=')[1] for elem in attributes}
		tmp_dict_attributes_dbxref = dict()
		try:
			tmp_dict_attributes_dbxref = {elem.split(':')[0]: elem.split(':')[1] for elem in
			                              tmp_dict_attributes['Dbxref'].split(',') if 'HGNC' not in elem.split(':')[0]}
		except KeyError:
			pass
		tmp_second_dict = {c: record[i] for i, c in
		                   enumerate(['ID', 'Source', 'Element', 'Start', 'End', 'Score', 'Brin', 'Phase'])}
		output_dict = {**tmp_dict_attributes, **tmp_second_dict, **tmp_dict_attributes_dbxref}
		return output_dict

	@staticmethod
	def process_df(data):
		"""
		Build dataframe and order columns
		Args:
			data (list): list of dicts

		Returns:
			df (DataFrame): pandas dataframe
		"""
		df = pd.DataFrame.from_records(data)
		ordered_cols = ['ID', 'Source', 'Element', 'Start', 'End', 'Score', 'Brin', 'Phase'] + list(
			sorted(set(list(df.columns)) - {'ID', 'Source', 'Element', 'Start', 'End', 'Score', 'Brin', 'Phase'}))
		df = df[ordered_cols]
		return df

	# @profile
	def get_coding_genes(self, raw_df):
		"""
		Key function of this script :
			- Get all protein coding genes and remove pseudogenes and other genes which are non coding for proteins
			- Iterate over all gene ids corresponding to previous mentioned requirements
			- For each gene id, add columns corresponding to
				- number of : mRNA, exon and CDS in each gene
				- statistics for each gene : (mean, median, std) size for each mRNA, exon and CDS
		Args:
			raw_df (DataFrame): complete DF with all RefSeq

		Returns:
			d_stats (dict): stats dictionnary

		"""
		# SIMPLE STEP TO GET ONLY PROTEIN CODING GENES
		gene_ids = list(sorted(list(set(raw_df[raw_df['gene_biotype'] == 'protein_coding']['GeneID'].tolist()))))
		df = raw_df[raw_df['GeneID'].isin(gene_ids)]

		# FILTER TO KEEP ONLY NC (COMPLETE GENOMIC SEQUENCES)
		ids_nc = [e for e in list(df.loc[df['Element'] == 'gene', 'ID'].unique()) if 'NC' in e]
		df = df.loc[df['ID'].isin(ids_nc)]

		# FILTER TO REMOVE OTHER GENE BIOTYPE
		others = df.loc[df['gene_biotype'] == 'other']
		ids_other = others['ID'].tolist()
		gene_ids_others = others['GeneID'].tolist()
		df.drop(list(df.loc[(df['ID'].isin(ids_other)) & (df['GeneID'].isin(gene_ids_others))].index), inplace=True)
		tmp = df['Genbank'].dropna().str.contains('NR')
		df.drop(list(tmp[tmp].index), inplace=True)
		df['GeneID'] = df['GeneID'].astype('category')

		# MP
		m = multiprocessing.Manager()
		complete_list_df = m.list()
		global_dict = m.dict()
		parmap.starmap(self.mp_refseq, list(zip(gene_ids)), df, complete_list_df, pm_pbar=True, pm_processes=self.cpus)

		# STATS
		d_stats = collections.defaultdict()
		for k, v in dict(global_dict).items():
			d_stats[k] = self.stats(v, k)

		# CONCATENATION
		complete_df = pd.concat(list(complete_list_df))

		# OUTPUTS
		complete_df.to_csv(self.output_name, compression='gzip', sep='\t', index=False, chunksize=10000)
		return d_stats

	def mp_refseq(self, gene, df, complete_list_df, filter_utr=False):
		"""
		MP function to process df with multiple genes at the same time
		Args:
			gene (str): RefSeq GeneID
			df (DataFrame): complete RefSeq DF
			complete_list_df (list): shared mp output list with sub-DF
			filter_utr (bool): enable or not the filtering of UTR

		Returns:
			None
		"""
		# CHECK IF GENE DF CORRESPONDS TO A PROTEIN CODING GENE
		if set(df[df['GeneID'] == gene]['Element'].unique()) == {'mRNA', 'exon', 'gene', 'CDS'}:
			df_gene = df.loc[(df['GeneID'] == gene) & (~df.Element.isin(['exon']))].reset_index(drop=True)

			# ADD LENGTH AND ID COLUMNS
			df_gene['Length'] = pd.to_numeric(df_gene['End']) - pd.to_numeric(df_gene['Start'])
			df_gene['Elem_position_ID'] = df_gene['Start'] + '_' + df_gene['End']

			# FILTER_UTR FUNCTION
			if filter_utr is True:
				if df_gene.loc[df_gene.Element == 'mRNA', 'Element'].shape[0] > 1:
					self.filter_mrna_utr(df_gene)

			# CHECK IF GENE IS ALWAYS THERE
			if 'gene' in df_gene['Element'].tolist():

				# GET LIST OF ALL PROTEINS PRODUCED
				prot_ids = df_gene[df_gene['Element'] == 'CDS']['Name'].unique().tolist()
				cds_list = list()
				mrna_list = list()

				# ITERATE OVER THE PROTEINS LIST
				for p, prot in enumerate(prot_ids):
					tmp_df_prot = df_gene[df_gene['Name'] == prot]

					# GET ALL CDS, ROW WITH NAME EQUAL TO PROT CORRESPOND TO CDS ROWS
					cds = df_gene[df_gene['Name'] == prot]
					df_gene.loc[list(cds.index), 'CDS_num'] = range(1, cds.shape[0] + 1)
					cds_list += [tuple(e) for e in cds[['Start', 'End']].values]

					# PARENT MRNA
					parent_mrna = tmp_df_prot['Parent'].unique()[0].split('-')[1]
					tmp_df_mrna_exons = df_gene[df_gene['transcript_id'] == parent_mrna]

					# MRNA
					mrna = tmp_df_mrna_exons[tmp_df_mrna_exons['Element'] == 'mRNA']

					# MRNA NUMEROTATION
					df_gene.loc[list(mrna.index), 'mRNA_num'] = p + 1
					mrna_list += [tuple(e) for e in mrna[['Start', 'End']].values]

				# MRNA NB COLUMN
				df_gene.loc[:, 'mRNA_nb'] = len(mrna_list)

				# COUNT OCCURENCE OF EACH CDS
				counter_cds = dict(collections.Counter(cds_list))
				counter_cds = {k[0] + '_' + k[1]: v for k, v in counter_cds.items()}
				total_count = len(counter_cds)
				cds_index = list(df_gene[df_gene['Element'] == 'CDS'].index)
				df_gene.loc[cds_index, 'Count'] = df_gene[df_gene['Element'] == 'CDS']['Elem_position_ID'].map(
					counter_cds)

				# COUNT ALTERNATIVE AND CONSTITUTIVE CDS
				counter_cds = {k: v for k, v in counter_cds.items() if v < len(mrna_list)}
				df_gene['Count_CDS_alternative'] = len(counter_cds)
				df_gene['Count_CDS_constitutive'] = total_count - len(counter_cds)

			# ADD TO SHARED LIST
			complete_list_df.append(df_gene)

	@staticmethod
	def filter_mrna_utr(df_gene):
		"""

		Args:
			df_gene (DataFrame): dataframe corresponding to entries gene, mRNA and CDS to a given gene

		Returns:
			df_gene (DataFrame): filter dataframe with mRNA modifications due to CDS changes
		"""

		nb_mrna_total = df_gene.loc[df_gene.Element == 'mRNA', 'Element'].shape[0]
		previous_element = str()
		final_comparison_list = list()
		tmp_l = list()
		mrna_id = str()
		current_element = str()
		dict_mrna_cds_tmp = collections.defaultdict(list)
		for index, row in df_gene.iterrows():
			if index > 0:
				previous_element = current_element
			current_element = str(row.Element)
			if current_element == 'CDS':
				tmp_l.append(str(row.Start) + '_' + str(row.End))
			if current_element == 'mRNA' and previous_element == 'gene':
				mrna_id = str(row.Start) + '_' + str(row.End)
			if current_element == 'mRNA' and previous_element == 'CDS':
				dict_mrna_cds_tmp[mrna_id] = "_".join(tmp_l)
				mrna_id = str(row.Start) + '_' + str(row.End)
				final_comparison_list.append("_".join(tmp_l))
				tmp_l = list()
			if current_element == 'CDS' and (index + 1) == df_gene.shape[0]:
				dict_mrna_cds_tmp[mrna_id] = "_".join(tmp_l)
				final_comparison_list.append("_".join(tmp_l))
		if len(list(set(final_comparison_list))) != nb_mrna_total:
			tmp_counter_mrna_cds = dict(collections.Counter(final_comparison_list))
			# tmp_counter_mrna_cds = [k for k,v in tmp_counter_mrna_cds.items() if v >= 2]
			final_dict = dict()
			for k, v in dict_mrna_cds_tmp.items():
				if v not in list(final_dict.values()):
					if tmp_counter_mrna_cds[v] >= 2:
						final_dict[k] = v
					else:
						final_dict[k] = v
			keep_mrna_list = list(final_dict.keys())
			df_gene = df_gene.drop(df_gene[(~df_gene.Elem_position_ID.isin(keep_mrna_list)) & (
					df_gene.Element == 'mRNA')].index).drop_duplicates(subset=['Element', 'Elem_position_ID'],
			                                                           keep='first')
		return df_gene

	def stats(self, list_values, elem):
		"""

		Args:
			list_values (list): list of tuples position
			elem: type of genomic element

		Returns:
			dict_stats (dict): stats with all precomputed statistics
		"""
		dict_stats = dict()
		dict_stats['nb'] = len(list_values)
		tmp_list_size = list()
		for e in list_values:
			tmp_list_size.append(int(e[1]) - int(e[0]))
		if tmp_list_size:
			dict_basic_stats = self.basic_stats(tmp_list_size)
			dict_stats = {**dict_stats, **dict_basic_stats}
			dict_stats["Element"] = elem
		return dict_stats

	@staticmethod
	def basic_stats(tmp_list):
		"""

		Args:
			tmp_list (list): list with sizes

		Returns:
			dict_stats (dict): dict with stats on sizes

		"""
		dict_stats = dict()
		dict_stats["mean_size"] = np.mean(tmp_list)
		dict_stats["max_size"] = np.max(tmp_list)
		dict_stats["min_size"] = np.min(tmp_list)
		dict_stats["median_size"] = np.median(tmp_list)
		dict_stats["std_size"] = np.std(tmp_list)
		return dict_stats

	def return_stats(self, ):
		"""
		Basic function to return stats of RefSeq file
		Returns:

		"""
		return self.stats


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
	RefSeqBuildFileCodingGenes(sys.argv[1], sys.argv[2])
