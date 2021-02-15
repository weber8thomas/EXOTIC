import sys
from pprint import pprint
import collections
import os
import warnings


warnings.simplefilter(action='ignore', category=FutureWarning)
import pandas as pd

pd.options.mode.chained_assignment = None  # default='warn'
from pandarallel import pandarallel
pandarallel.initialize(progress_bar=True, nb_workers=12)
from tqdm import tqdm

pd.set_option('display.float_format', lambda x: '%.3f' % x)
from src.utils import utils
from src.gene_exon_level_analysis.gene_exon_level import GeneExonLevel

pd.set_option('display.float_format', lambda x: '%.3f' % x)

class RefSeqStats(GeneExonLevel):

	def __init__(self, ref_file, disease_genes, corrected_df, filter_multi_iso_df, disease_df, healthy_df):

		self.config = utils.load_config_file()
		self.logger.info('Producing RefSeq global statistics')
		self.disease_genes = disease_genes
		self.corrected_df, self.filter_multi_iso_df, self.disease_df, self.healthy_df = self.refseq_correction_and_global_stats(ref_file)

	def refseq_correction_and_global_stats(self, ref_file):
		# CORRECTION OF REFSEQ EXONS
		if os.path.isfile(self.config['RefSeqStats']['CORRECTED_DF_PATH']) is False:
			data = self.read_df(ref_file)
			corrected_df = self.correct_df(data)
			self.write_df(corrected_df, self.config['RefSeqStats']['CORRECTED_DF_PATH'])
			self.stats_df(corrected_df, self.config['RefSeqStats']['STATS_PATH'], subset=True)

			filter_multi_iso_df = self.filter_multi_isoforms(corrected_df)
			self.write_df(filter_multi_iso_df, self.config['RefSeqStats']['CORRECTED_MULTI_ISO_PATH'])
			self.stats_df(filter_multi_iso_df, self.config['RefSeqStats']['STATS_MULTI_ISO_PATH'], subset=True)

			disease_df, healthy_df = self.filter_disease_healthy_genes(filter_multi_iso_df, self.disease_genes)
			self.write_df(disease_df, self.config['RefSeqStats']['CORRECTED_DISEASE_PATH'])
			self.write_df(healthy_df, self.config['RefSeqStats']['CORRECTED_HEALTHY_PATH'])
			self.stats_df(disease_df, self.config['RefSeqStats']['STATS_DISEASE_PATH'], subset=True)
			self.stats_df(healthy_df, self.config['RefSeqStats']['STATS_HEALTHY_PATH'], subset=True)

			self.merge_stats_of_different_conditions()
		else:
			corrected_df = self.read_df(self.config['RefSeqStats']['CORRECTED_DF_PATH'])
			# corrected_df = ""
			filter_multi_iso_df = self.read_df(self.config['RefSeqStats']['CORRECTED_DF_PATH'])
			# filter_multi_iso_df = ""

			disease_df = self.read_df(self.config['RefSeqStats']['CORRECTED_DF_PATH'])
			healthy_df = self.read_df(self.config['RefSeqStats']['CORRECTED_DF_PATH'])
		return corrected_df, filter_multi_iso_df, disease_df, healthy_df

	# CORRECTION OF CDS LOCATIONS
	def correct_df(self, df):
		tqdm.pandas()

		output_df = df[['ID', 'Source', 'Element', 'Start', 'End', 'Score', 'Brin', 'Phase', 'Dbxref', 'GeneID', 'MIM',
		                'Name', 'description', 'gbkey', 'gene', 'gene_biotype', 'gene_synonym', 'Length',
		                'Elem_position_ID',
		                'mRNA_nb', 'Count', 'Count_CDS_alternative', 'Count_CDS_constitutive', 'mRNA_IDS',
		                'CDS_locations_corrected']]

		output_df['HGNC'] = output_df['Dbxref'].apply(lambda r: self.add_hgnc(r))
		output_df['Count'] = output_df['CDS_locations_corrected'].parallel_apply(lambda r: len(eval(r)))
		output_df['CDS_locations_corrected'] = output_df['CDS_locations_corrected'].progress_apply(
			lambda r: self.check_integrity(r))
		# df['CDS_locations_corrected'] = df['CDS_locations_corrected'].parallel_apply(lambda r: self.check_integrity(r))
		output_df = output_df.parallel_apply(lambda r: self.apply_correct_cds(r), axis=1)

		output_df[["CDS_Unique_nb", "CDS_Variable_region_nb", ]] = output_df[
			["CDS_Unique_nb", "CDS_Variable_region_nb", ]].fillna(0)
		output_df["CDS_Unique_nb"] = output_df["CDS_Unique_nb"].astype(int)
		output_df["CDS_Variable_region_nb"] = output_df["CDS_Variable_region_nb"].astype(int)
		output_df['CDS_nb'] = output_df["CDS_Unique_nb"] + output_df["CDS_Variable_region_nb"]
		return output_df

	# ADD HGNC COLUMN
	@staticmethod
	def add_hgnc(row):
		row = [e.replace('HGNC:HGNC:', '') for e in row.split(',') if 'HGNC' in e]
		if row:
			return row[0]
		else:
			return None

	# CHECK INTEGRITY OF REFSEQ FILE
	@staticmethod
	def check_integrity(row):
		row = eval(row)
		for cds, cds_content in row.items():
			# FIRST CONDITION - SHARING STATUS -  IF MORE THAN ONE SOURCE, THEN EXON SHARED
			if len(cds_content['Sources']) == 1:
				cds_content['Sharing_status'] = False
			else:
				cds_content['Sharing_status'] = True

			# # SECOND CONDITION - SHARING STATUS & CDS REPRESENTATION - IF EXON SHARED, THEN CDS VARIABLE CONSTANT
			# if (cds_content['Sharing_status'] is True) and (cds_content['CDS_representation'] == 'Unique'):
			# 	cds_content['CDS_representation'] = 'Variable_constant'

			# THIRD CONDITION - IF CDS REPRESENTATION UNIQUE, THEN EXON NOT SHARED
			if cds_content['CDS_representation'] == 'Unique':
				cds_content['Sharing_status'] = False

		# # FOURTH CONDITION - IF CDS REPRESENTATION DISPENSABLE, THEN EXON SHARED
		# if cds_content['CDS_representation'] == 'Variable_dispensable':
		# 	cds_content['Sharing_status'] = True
		return row

	# CORRECT CDS
	def apply_correct_cds(self, row):
		# print(row)
		size_dict_of_list = collections.defaultdict(list)
		ratio_dict_of_list = collections.defaultdict(list)
		row_CDS_locations_corrected = row['CDS_locations_corrected']
		# row_CDS_locations_corrected = eval(row['CDS_locations_corrected'])
		for k, v in row_CDS_locations_corrected.items():
			size_dict_of_list[v['CDS_representation']].append(k)
			ratio_dict_of_list[v['CDS_representation']].append(eval(v['Ratio']))
		for category in size_dict_of_list:
			pd_series_size = self.exons_stats(size_dict_of_list[category], 'CDS_' + category, 'size')
			pd_series_ratio = self.exons_stats(ratio_dict_of_list[category], 'CDS_' + category, 'ratio')
			row = row.append([pd_series_size, pd_series_ratio])
		return row


	# FILTER REFSEQ WITH MULTI-ISO GENES
	def filter_multi_isoforms(self, df):
		self.logger.info('RSS - Filtering multi isoforms genes')

		if set([os.path.isfile(self.config['RefSeqStats'][file]) for file in ['mRNA_nb_counter', 'mRNA_nb_metrics']]) != {True}:
			# STATS
			df_stats_output = 100 * (
					df.loc[df['mRNA_nb'] >= 2, 'mRNA_nb'].value_counts() / df.loc[df['mRNA_nb'] >= 2].shape[
				0])  # STATS OUTPUT
			df_stats_output.to_excel(self.config['RefSeqStats']['mRNA_nb_counter'])
			self.stats_df(df.loc[df['mRNA_nb'] >= 2, 'mRNA_nb'], self.config['RefSeqStats']['mRNA_nb_metrics'])

		# FILTERING
		corrected_df = df.loc[df['mRNA_nb'] >= 2]
		return corrected_df

	# FILTER REFSEQ WITH DISEASE GENES BY COMPARING HGNC IDS
	@staticmethod
	def filter_disease_healthy_genes(corrected_df, disease_genes):
		corrected_df['HGNC'] = corrected_df['HGNC'].astype(str)
		corrected_df['HGNC'] = corrected_df['HGNC'].str.replace('.0', '')
		disease_df = corrected_df.loc[corrected_df['HGNC'].isin(disease_genes)]
		healthy_df = corrected_df.loc[~corrected_df['HGNC'].isin(disease_genes)]
		return disease_df, healthy_df

	def merge_stats_of_different_conditions(self):
		# FILTERING STATS FILE ONLY
		STATS_FILES = {file:v for file, v in self.config.items() if 'STATS' in file}
		pprint(STATS_FILES)

		# ITERATING
		final_df = list()
		name_list = list()
		for file, path in STATS_FILES.items():
			tmp_stats_df = pd.read_excel(path).T
			for col in tmp_stats_df.columns:
				if col == 'count':
					tmp_stats_df[col] = tmp_stats_df[col].astype(int)
				else:
					tmp_stats_df[col] = tmp_stats_df[col].round(2)

			name = path.split('/')[-1].replace('GRCh37_stats_', '').replace('.xlsx', '')
			name_list.append(name)
			final_df.append(tmp_stats_df)
		cols = list(tmp_stats_df.columns)
		final_df = pd.concat(final_df, axis=1)
		final_df.columns = pd.MultiIndex.from_product([name_list, cols])
		print("/".join(path.split('/')[:-1]))
		final_df.to_excel("/".join(path.split('/')[:-1]) + '/GRCh37_stats_complete.xlsx')


if __name__ == "__main__":
	RefSeqStats(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6],)
