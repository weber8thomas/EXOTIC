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
from pprint import pprint
import pandas as pd
import os
import requests
import time

omim_file = pd.read_csv("https://omim.org/static/omim/data/mim2gene.txt", skiprows=4, sep='\t')
omim_full = omim_file[omim_file.columns[0]].values

omim_api = "https://api.omim.org/api/entry?mimNumber={}&include={}&format=json&apiKey=fpUv2M2QTmCGhuPexCO-MA"

output_dir = "/gstock/biolo_datasets/variation/benchmark/Databases/OMIM/JSON_API/"

refseq_x_vcf = pd.read_parquet('/home/weber/PycharmProjects/ExoCarto/data/2_processed/refseq_vcf_pext_phylocsf.parquet')
disease_genes = refseq_x_vcf.CCRS_Gene.unique()

refseq = pd.read_csv('/home/weber/PycharmProjects/ExoCarto/data/2_processed/GRCh37_RefSeq_lite.csv.gz', compression='gzip', sep='\t').dropna(subset=['MIM'])
refseq['MIM'] = refseq['MIM'].astype(int)
refseq = refseq.loc[(refseq['mRNA_nb'] > 1) & (refseq['Name'].isin(disease_genes))]


omim_list_complete = refseq.MIM.values
check_downloaded_omim = [int(e.replace('.pkl', '')) for e in os.listdir(output_dir)]
omim_to_download = [int(e) for e in omim_full if int(e) not in check_downloaded_omim]


def check_downloaded_files(check_downloaded_omim):
	for omim in tqdm(check_downloaded_omim):
		json = _pickle.load(open(output_dir + '{}.pkl'.format(str(omim)), 'rb'))
		if str(json)[:7] != "{'omim'":
			print(json)
			exit()


def download_mp(gene):
	# gene = omim_list[0]
	response = requests.get(omim_api.format(str(gene), "all"))
	json = response.json()
	_pickle.dump(json, open(output_dir + "{}.pkl".format(gene), 'wb'))


# DOWNLOAD PART
j = 0
from tqdm import trange
t = trange(len(omim_to_download), desc='Download OMIM', leave=True)
for j, gene in enumerate(omim_to_download):
	t.set_description('Download OMIM : {} download IDs'.format(str(j)))
	download_mp(gene)
	if j % 20 == 0:
		time.sleep(5)
# parmap.starmap(download_mp, list(zip(omim_list)), pm_pbar=True, pm_processes=60)

# CHECK PART
check_downloaded_files(check_downloaded_omim)


