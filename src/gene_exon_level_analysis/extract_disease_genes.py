import _pickle
import collections
import os
import sys
import warnings
from concurrent.futures import ThreadPoolExecutor

from cyvcf2 import VCF, Writer

warnings.simplefilter(action='ignore', category=FutureWarning)
import pandas as pd

pd.options.mode.chained_assignment = None  # default='warn'
from pandarallel import pandarallel

pandarallel.initialize(progress_bar=True, nb_workers=1)
from tqdm import tqdm

pd.set_option('display.float_format', lambda x: '%.3f' % x)

from src.utils import utils
from src.gene_exon_level_analysis.gene_exon_level import GeneExonLevel


class ExtractDiseaseGenes(GeneExonLevel):

	def __init__(self, CLINVAR_PATH, HGMD_PATH, output):

		self.CLINVAR_PATH = CLINVAR_PATH
		self.HGMD_PATH = HGMD_PATH
		self.logger.info('Extracting disease genes from ClinVar & HGMD')
		self.disease_genes = self.get_disease_genes()
		self.config = utils.load_config_file()

	def get_disease_genes(self):
		# GET DISEASE GENES IN CLINVAR & HGMD
		if os.path.isfile(self.config['ExtractDiseaseGenes']['DISEASE_GENES_PATH']) is False:
			with ThreadPoolExecutor(2) as executor:
				list_disease_genes_clinvar = executor.submit(lambda p: self.get_disease_gene_with_vcf(*p),
				                                             [self.CLINVAR_PATH])
				list_disease_genes_hgmd = executor.submit(lambda p: self.get_disease_gene_with_vcf(*p),
				                                          [self.HGMD_PATH])
			list_disease_genes_clinvar = list_disease_genes_clinvar.result()
			list_disease_genes_hgmd = list_disease_genes_hgmd.result()
			disease_genes = list(set(list_disease_genes_clinvar).union(set(list_disease_genes_hgmd)))
			# self.venn_gene_sets(set(list_disease_genes_clinvar), set(list_disease_genes_hgmd), ['ClinVar', 'HGMD'], )
			with open(self.config['ExtractDiseaseGenes']['DISEASE_GENES_PATH'], "wb") as f:
				_pickle.dump(disease_genes, f)
		else:
			with open(self.config['ExtractDiseaseGenes']['DISEASE_GENES_PATH'], "rb") as f:
				disease_genes = _pickle.load(f)
		return disease_genes

	# FILTER VCF FILE TO GET ONLY DISEASE IMPLICATED GENES
	def get_disease_gene_with_vcf(self, file):
		vcf = VCF(file)
		output_dir = os.path.dirname(self.config['ExtractDiseaseGenes']['DISEASE_GENES_PATH'])
		w = Writer(output_dir + '/' + os.path.basename(file).replace('.vcf.gz', '_pathogenic.vcf.gz'), vcf)

		vep_field = 'CSQ'
		vep_separator = '|'
		list_disease_genes = list()
		list_disease_genes_symbol = list()

		index_dict = self.get_header(vcf, vep_field, vep_separator)
		c = list()
		name = 'clinvar' if 'clinvar' in file else 'hgmd'
		pbar = False if name == 'hgmd' else True
		for counter, variant in enumerate(tqdm(vcf, desc=name)):
			# if counter == 10000:
			# 	break
			if len(variant.REF) == 1 and len(variant.ALT[0]) == 1 and vep_field:
				if 'clinvar' in file:
					cln_sig = str(variant.INFO.get('CLNSIG')).lower()
					c.append(cln_sig)
					for state in ['Pathogenic', ]:
						if cln_sig and state[1:] in cln_sig:
							if 'conflicting' not in cln_sig:
								w.write_record(variant)
								for case in variant.INFO.get('CSQ').split(','):
									case = case.split(vep_separator)
									hgnc = case[index_dict['HGNC_ID']]
									genename = case[index_dict['SYMBOL']]
									list_disease_genes.append(hgnc)
									list_disease_genes_symbol.append(genename)

				if 'hgmd' in file:
					cln_sig = variant.INFO.get('CLASS')
					c.append(cln_sig)
					if cln_sig and cln_sig == 'DM':
						w.write_record(variant)
						for case in variant.INFO.get('CSQ').split(','):
							case = case.split(vep_separator)
							hgnc = case[index_dict['HGNC_ID']]
							genename = case[index_dict['SYMBOL']]
							list_disease_genes.append(hgnc)
							list_disease_genes_symbol.append(genename)

		# self.dump_counter(collections.Counter(c), self.LOCAL_DIR, 'Stats_VCF_CLINSIG_' + name)
		# self.dump_counter(collections.Counter(list_disease_genes), self.LOCAL_DIR, 'Stats_VCF_GENES_HGNC_' + name)
		# self.dump_counter(collections.Counter(list_disease_genes_symbol), self.LOCAL_DIR, 'Stats_VCF_GENES_SYMBOL_' + name)

		return list(set(list_disease_genes))


if __name__ == "__main__":
	ExtractDiseaseGenes(sys.argv[1], sys.argv[2], sys.argv[3])
