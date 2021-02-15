import sys
import warnings

import numpy as np

warnings.simplefilter(action='ignore', category=FutureWarning)
import pandas as pd

pd.options.mode.chained_assignment = None  # default='warn'
from pandarallel import pandarallel
pandarallel.initialize(progress_bar=True, nb_workers=1)
from src.utils import utils

pd.set_option('display.float_format', lambda x: '%.3f' % x)
from matplotlib_venn import venn2
# import matplotlib.pyplot as plt

# from pprint import pprint

class GeneExonLevel:
	# GLOBAL INSTANCIATIONS

	logger = utils.setup_custom_logger('ExoCarto.txt')
	# pprint(utils.load_config_file())

	def __init__(self):
		self.logger.info('----- GeneExonLevel part analysis starting ... -----')

	@staticmethod
	def exons_stats(list_values, elem, type):
		dict_stats = dict()
		# print(list_values)
		CDS_nb = len(list_values)
		tmp_list_size = list()
		if type == 'size':
			dict_stats[elem + '_nb'] = CDS_nb
			for e in list_values:
				tmp_list_size.append(int(e[1]) - int(e[0]))
		else:
			tmp_list_size = list_values
		if tmp_list_size:
			dict_stats[elem + "_mean_" + type] = np.mean(tmp_list_size)
			dict_stats[elem + "_max_" + type] = np.max(tmp_list_size)
			dict_stats[elem + "_min_" + type] = np.min(tmp_list_size)
			dict_stats[elem + "_median_" + type] = np.median(tmp_list_size)
			dict_stats[elem + "_std_" + type] = np.std(tmp_list_size)
			dict_stats[elem + "_total_" + type] = np.sum(tmp_list_size)

		# print(list_values, len([CDS_nb, CDS_mean_size, CDS_max_size, CDS_min_size, CDS_median_size,
		# CDS_std_size]), CDS_nb, CDS_mean_size, CDS_max_size, CDS_min_size, CDS_median_size, CDS_std_size)
		# print(pd.DataFrame.from_dict(dict_stats, orient='index').T) exit()
		return pd.Series(dict_stats)

	# READ FILE
	@staticmethod
	def read_df(input_file, full=True):
		if full is True:
			df = pd.read_csv(input_file, compression='gzip', sep='\t', low_memory=False)
		else:
			df = pd.read_csv(input_file, compression='gzip', sep='\t', low_memory=False, nrows=3000)
		# df = pd.read_csv(file, compression='gzip', sep='\t', low_memory=False)
		return df

	# WRITE FILE
	@staticmethod
	def write_df(df, output):
		df.to_csv(output, compression='gzip', sep='\t', index=False)

	# PRODUCE STATS AND DUMP IT INTO EXCEL FILE
	@staticmethod
	def stats_df(df, name, subset=False):
		if subset is True:
			df = df[
				['CDS_Unique_max_ratio', 'CDS_Unique_max_size', 'CDS_Unique_mean_ratio', 'CDS_Unique_mean_size',
				 'CDS_Unique_median_ratio', 'CDS_Unique_median_size', 'CDS_Unique_min_ratio', 'CDS_Unique_min_size',
				 'CDS_Unique_nb', 'CDS_Unique_std_ratio', 'CDS_Unique_std_size', 'CDS_Unique_total_ratio',
				 'CDS_Unique_total_size', 'CDS_Variable_region_max_ratio', 'CDS_Variable_region_max_size',
				 'CDS_Variable_region_mean_ratio', 'CDS_Variable_region_mean_size',
				 'CDS_Variable_region_median_ratio', 'CDS_Variable_region_median_size',
				 'CDS_Variable_region_min_ratio', 'CDS_Variable_region_min_size', 'CDS_Variable_region_nb',
				 'CDS_Variable_region_std_ratio', 'CDS_Variable_region_std_size', 'CDS_Variable_region_total_ratio',
				 'CDS_Variable_region_total_size', 'Count', 'Count_CDS_alternative', 'Count_CDS_constitutive',
				 'mRNA_nb', 'CDS_nb']]
		describe = pd.concat([df.describe(), df.agg(['sum'])], axis=0)
		describe.to_excel(name, index=True)

	# PROCESS HEADER OF VCF FILE AND RETURN DICT OF VEP CONSTITUTION FIELD
	@staticmethod
	def get_header(vcf, vep_field, vep_separator):
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
	def dump_counter(counter, output_dir, name):
		if output_dir.endswith('/') is False:
			output_dir = output_dir + '/'
		pd.DataFrame.from_dict(counter, orient='index').reset_index().to_excel(output_dir + name + '.xlsx')

	def venn_gene_sets(self, setA, setB, labels):
		venn2([setA, setB], set_labels=labels)
		# plt.savefig(self.OUTPUT_FILES['VENN_DIAGRAM_CLINVAR_HGMD'])

	@staticmethod
	def getOverlap(a, b):
		"""
		Compute overlap between two intervals
		Args:
			a(tuple): tuple of int
			b(tuple): tuple of int
		Returns:

		"""
		return max(0, min(a[1], b[1]) - max(a[0], b[0]))

if __name__ == "__main__":
	GeneExonLevel(sys.argv[1])