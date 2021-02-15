import collections
import os
import sys
import warnings
import parmap
import multiprocessing

warnings.simplefilter(action='ignore', category=FutureWarning)
import pandas as pd

pd.options.mode.chained_assignment = None  # default='warn'
from pandarallel import pandarallel

pandarallel.initialize(progress_bar=True, nb_workers=1)
from tqdm import tqdm
from src.utils import utils
from src.gene_exon_level_analysis.gene_exon_level import GeneExonLevel

pd.set_option('display.float_format', lambda x: '%.3f' % x)


class CCRS_analysis(GeneExonLevel):

	def __init__(self, multi_iso_refseq_path, BIOMART_FILE, CCRS_FILE):
		self.logger.info('----- NucleotideLevel part analysis starting ... -----')

		self.logger.info('Building CCRS X RefSeq file')
		multi_iso_refseq = pd.read_csv(multi_iso_refseq_path, compression='gzip', sep='\t')
		self.config = utils.load_config_file()
		self.launch_class_and_check_files(multi_iso_refseq, BIOMART_FILE, CCRS_FILE)

	def launch_class_and_check_files(self, multi_iso_refseq, BIOMART_FILE, CCRS_FILE):
		## TRANSFORMATION FROM REFSEQ ONE LINE PER GENE TO ONE LINE PER CDS

		# IF REFSEQ MELT FILE NOT PRODUCED
		if os.path.isfile(self.config['CCRS_analysis']['RefSeq_transformed']) is False:
			# multi_iso_refseq = read_file(file, full=False)

			# ADD HGNC IDS
			# multi_iso_refseq = self.process_refseq_hgnc(multi_iso_refseq)
			# multi_iso_refseq = multi_iso_refseq.loc[multi_iso_refseq['HGNC'].isin(['25284', '24149'])] # TEST ON MULTI ISOFORM GENE
			multi_iso_refseq = self.process_refseq_count(multi_iso_refseq)

			# MELT DF
			df = self.function_transform_dataframe(multi_iso_refseq)

			# JOIN CHROM BASED ON GENE THANKS TO BIOMART
			df = self.join_chromosome(df, BIOMART_FILE, self.config['CCRS_analysis']['RefSeq_transformed'])

		# IF FILE PRODUCED, THEN READ IT
		else:
			df = GeneExonLevel.read_df(self.config['CCRS_analysis']['RefSeq_transformed'])
		# scales = scaling(df)

		## JOIN CCRS DATA TO REFSEQ

		# IF RefSeq X CCRS file not produced
		if os.path.isfile(self.config['CCRS_analysis']['RefSeq_CCRS'].replace('.csv.gz', '1.csv.gz')) is False:

			# Process CCRS file
			ccrs = self.process_ccrs_file(CCRS_FILE)

			# Filter to work on multiple isoforms
			refseq = self.filter_multi_isoforms(df)

			# Join CCRS & RefSeq CDS
			self.get_overlaps_between_files(refseq, ccrs, self.config['CCRS_analysis']['RefSeq_CCRS'])
		else:
			join_dfs = GeneExonLevel.read_df(self.config['CCRS_analysis']['RefSeq_CCRS'], full=False)
			print(join_dfs)

	# FIRST PART - CONVERT REFSEQ LITE FILE TO CDS (ONE ROW PER LINE) COMPLETE FILE

	@staticmethod
	def apply_count(row):
		row = eval(row)
		ratios = dict(collections.Counter([eval(v['Ratio']) for k, v in row.items()]))
		if 1. in list(ratios.keys()):
			const = ratios[1.]
		else:
			const = 0
		if len(list(ratios.keys())) > 1:
			alt = sum([v for k, v in ratios.items() if k < 1.])
		else:
			alt = 0
		return pd.Series([const, alt])

	def process_refseq_count(self, df):
		df[['Count_CDS_constitutive', 'Count_CDS_alternative']] = df['CDS_locations_corrected'].parallel_apply(
			lambda r: self.apply_count(r))
		return df

	def function_transform_dataframe(self, df):
		list_region = list()

		# ITERATE OVER EACH GENE ROW
		for j, row in tqdm(df.iterrows()):
			row_dict = row.to_dict()
			cds_locations_corrected = eval(row_dict['CDS_locations_corrected'])
			row_dict = {k: v for k, v in row.to_dict().items() if k != 'CDS_locations_corrected'}
			row_dict['Gene_Start'] = row_dict.pop('Start')
			row_dict['Gene_End'] = row_dict.pop('End')
			tmp_dict_sub_row = dict()
			for sub_row, sub_row_content in cds_locations_corrected.items():
				tmp_dict_sub_row['Start'] = sub_row[0]
				tmp_dict_sub_row['End'] = sub_row[1]
				list_region.append({**tmp_dict_sub_row, **sub_row_content, **row_dict})

		# BUILD NEW DF FROM LIST OF DICT PRODUCED
		df = pd.DataFrame(list_region)
		df = self.check_integrity(df)
		df.rename({'gene': 'Gene'}, axis=1, inplace=True)
		return df

	@staticmethod
	def check_sharing_status(row):
		if len(row) == 1:
			return False
		else:
			return True

	def check_integrity(self, df):
		df['Sharing_status'] = df['Sources'].apply(self.check_sharing_status)
		df.loc[(df['Sharing_status'] == True) & (
				df['CDS_representation'] == 'Unique'), 'CDS_representation'] = 'Variable_region'
		df.loc[df['CDS_representation'] == 'Unique', 'Sharing_status'] = False
		df.loc[df['CDS_representation'] == 'Variable_region', 'Sharing_status'] = True
		return df

	@staticmethod
	def join_chromosome(df, biomart_file, output_file):
		"""
		Args:
			df:
			biomart_file:
			output_file:

		Returns:

		"""
		# Biomart file with NCBI Gene ID and Chrom names
		biomart = pd.read_csv(biomart_file, sep='\t', compression='gzip')[['NCBI gene ID', 'Chromosome/scaffold name']]
		biomart.columns = ['GeneID', 'Chrom']

		# Join RefSeq and Biomart files to get chrom on RefSeq
		join_dfs = pd.merge(df, biomart, how='inner', on='GeneID')

		# Add Length and Range columns
		join_dfs['Length'] = join_dfs['End'] - join_dfs['Start']
		join_dfs['ranges'] = join_dfs['Start'].astype(str) + '-' + join_dfs['End'].astype(str)

		# Reorder, sort and export
		join_dfs = join_dfs.sort_values(by=['Chrom', 'Start', 'End'])
		join_dfs.to_csv(output_file, compression='gzip', sep='\t', index=False)

		return join_dfs

	# SECOND PART - PREPARE CCRS FILE AND MAP IT ON REFSEQ

	def process_ccrs_file(self, ccrs_file):
		# READ FILE
		ccrs = GeneExonLevel.read_df(ccrs_file, full=False)

		# RENAME AND SELECT COLUMNS
		ccrs.rename({'#chrom': 'Chrom', 'gene': 'Gene', 'start': 'Start', 'end': 'End', 'ccr_pct': 'CCR_percentile'},
		            axis=1, inplace=True)
		ccrs = ccrs[['Chrom', 'Start', 'End', 'CCR_percentile', 'Gene', 'ranges']]
		return ccrs

	def mp_overlap_refseq_ccrs(self, gene, refseq, ccrs, return_list):
		test_list = list()

		# SELECTION OF REFSEQ AND CCRS PART FILE CORRESPONDING TO THE CHOSEN GENE
		tmp_refseq = refseq.loc[refseq['Gene'] == gene]
		tmp_ccrs = ccrs.loc[ccrs['Gene'] == gene]
		tmp_range_comparison = collections.defaultdict(list)

		# ITERATE OVER EACH CDS ENTRY
		for i, row_refseq in tmp_refseq.iterrows():

			# BUILDING TUPLE AND DICT INFO
			refseq_tuple = (row_refseq['Start'], row_refseq['End'])
			coverage_refseq_tuple = int(refseq_tuple[1] - refseq_tuple[0])
			tmp_dict_refseq = row_refseq.to_dict()
			tmp_dict_refseq = {k: v for k, v in tmp_dict_refseq.items() if
			                   k in ["CDS_representation", "HGNC", "ranges", "Ratio", "Chrom", "Start", "End",
			                         "CDS_locations_corrected"]}
			tmp_dict_refseq = {'RefSeq_' + k: v for k, v in tmp_dict_refseq.items()}

			# ITERATE OVER CCRS GENE PART TO COMPARE INTERVALS
			for j, row_ccrs in tmp_ccrs.iterrows():

				# IF RANGE ALREADY USE, NEXT
				# if row_ccrs['ranges'] not in ranges_list_not_to_reparse:  # WARNING IF MP
				tmp_dict_ccrs = row_ccrs.to_dict()
				tmp_dict_ccrs = {'CCRS_' + k: v for k, v in tmp_dict_ccrs.items()}
				ccrs_tuple = (row_ccrs[1], row_ccrs[2])

				# COMPARE INTERVALS CCRS AND REFSEQ
				if GeneExonLevel.getOverlap(refseq_tuple, ccrs_tuple) > 0:
					# ADD NEW ROW TO DF CORRESPONDING TO OVERLAP
					tmp_range_comparison[refseq_tuple].append(ccrs_tuple)
					test_list.append({**tmp_dict_refseq, **tmp_dict_ccrs})

		# BUILD DF FROM LIST OF DICTS
		tmp_df = pd.DataFrame(test_list)

		if tmp_df.empty is False:

			# INSTANCIATION OF SMALL DICTS FOR COVERAGE AND ORIGINS OF GAPS
			dict_coverage = dict()
			dict_previous = dict()

			# FIRST LOOP ON REFSEQ RANGES
			for refseq_tuple, ccrs_tuple_list in tmp_range_comparison.items():

				# COMPUTING COVERAGES
				coverage_refseq_tuple = int(refseq_tuple[1] - refseq_tuple[0])
				coverage_ccrs_tuple = sum([int(t[1] - t[0]) for t in ccrs_tuple_list])

				# IF CCRS COVERAGE STRICTLY BELOW REFSEQ
				if coverage_ccrs_tuple < coverage_refseq_tuple:

					# SECOND LOOP ON TUPLE LIST
					previous = tuple()
					for t in ccrs_tuple_list:

						# FIRST ITERATION
						if not previous:
							previous = t

						# SECOND AND OTHER ITERATIONS
						elif previous:

							# IF LAST CCRS RANGE END != FROM CURRENT CCRS START
							if t[0] != previous[1]:
								# FILL DICTS
								dict_coverage['-'.join([str(e) for e in list(t)])] = t[0] - previous[1]
								dict_previous['-'.join([str(e) for e in list(t)])] = '-'.join(
									[str(e) for e in list(previous)])

							# INCREMENT
							previous = t

			try:
				# NEW COLUMNS WITH COVERAGE AND ORIGINS OF GAPS
				tmp_df['Missing_coverage_bp'] = tmp_df['CCRS_ranges'].map(dict_coverage).fillna(0)
				tmp_df['Gap_with_previous_CCRS'] = tmp_df['CCRS_ranges'].map(dict_previous)
			except KeyError:
				print('Error on gene {}'.format(str(gene)))
				print(tmp_df)

		# RETURN
		return_list.append(tmp_df)

	# FILTER REFSEQ WITH MULTI-ISO GENES
	@staticmethod
	def filter_multi_isoforms(df):
		corrected_df = df.loc[df['mRNA_nb'] >= 2]
		return corrected_df

	def get_overlaps_between_files(self, refseq, ccrs, output_file):
		# CATEGORIZATION OF COLUMNS
		# refseq['Gene'] = refseq['Gene'].astype('category')
		# ccrs['Gene'] = ccrs['Gene'].astype('category')
		refseq['Chrom'] = refseq['Chrom'].astype('str')
		ccrs['Chrom'] = ccrs['Chrom'].astype(str)
		# ccrs['Chrom'] = ccrs['Chrom'].astype('category')
		gene_list_ccrs = set(list(ccrs['Gene'].unique()))
		gene_list_refseq = set(list(refseq['Gene'].unique()))

		# TODO : Check correspondance between REFSEQ AND CCRS files
		# INTERSECTION BETWEEN REFSEQ AND CCRS
		# venn_gene_sets(gene_list_ccrs, gene_list_refseq, ['CCRS', 'RefSeq'])
		intersection_genes_ccrs_refseq = list(set(gene_list_ccrs).intersection(set(gene_list_refseq)))

		# pprint(list(set(gene_list_ccrs)))
		# pprint(list(set(gene_list_refseq)))
		# pprint(intersection_genes_ccrs_refseq)
		# print(len(intersection_genes_ccrs_refseq))

		# TUPLE OF CHROM AND GENES
		df_genes_chr = pd.DataFrame(list(set(list(zip(refseq.Chrom, refseq.Gene)))),
		                            columns=['Chrom', 'Gene']).sort_values(by='Chrom')
		df_genes_chr = df_genes_chr.loc[df_genes_chr['Gene'].isin(intersection_genes_ccrs_refseq)]

		m = multiprocessing.Manager()

		# ITERATE OVER CHROM
		for chrom in tqdm(list(df_genes_chr.Chrom.unique())):
			# if chrom == '13':
			if chrom:
				# 	print(chrom, type(chrom))

				# LIST OF GENES TO ITERATE ON THE GIVEN CHROM
				list_genes = list(df_genes_chr.loc[df_genes_chr['Chrom'] == chrom, 'Gene'].values)

				# SELECTION OF ALL GENES ON REFSEQ AND CCRS ON THIS CHROM
				refseq_tmp_chrom = refseq.loc[refseq['Chrom'] == chrom]
				ccrs_tmp_chrom = ccrs.loc[ccrs['Chrom'] == chrom]

				# MP
				chrom_list_return = m.list()
				parmap.starmap(self.mp_overlap_refseq_ccrs, list(zip(list_genes)), refseq_tmp_chrom, ccrs_tmp_chrom,
				               chrom_list_return, pm_pbar=True, pm_processes=multiprocessing.cpu_count() - 4)

				# JOIN RESULTS
				chrom_list_return = pd.concat(list(chrom_list_return))
				# chrom_list_return = chrom_list_return.drop('CCRS_Chrom', axis=1)

				# DUMP
				if chrom_list_return.empty is False:
					chrom_list_return = chrom_list_return[
						['CCRS_Gene', 'RefSeq_HGNC', 'RefSeq_Chrom', 'RefSeq_ranges', 'RefSeq_Start', 'RefSeq_End',
						 'RefSeq_Ratio', 'RefSeq_CDS_representation', 'CCRS_ranges', 'CCRS_Start', 'CCRS_End',
						 'CCRS_CCR_percentile', 'Missing_coverage_bp', 'Gap_with_previous_CCRS', ]]
					chrom_list_return.sort_values(by=['RefSeq_Chrom', 'RefSeq_Start', 'RefSeq_End']).to_csv(
						output_file.replace('.csv.gz', '_chrom{}'.format(str(chrom)) + '.csv.gz'), compression='gzip',
						sep='\t', index=False)

			# print(refseq.to_string())
			# print(ccrs.to_string())


if __name__ == "__main__":
	CCRS_analysis(sys.argv[1], sys.argv[2], sys.argv[3])
