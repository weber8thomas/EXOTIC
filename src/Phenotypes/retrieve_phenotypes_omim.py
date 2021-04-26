import parmap
import time
import requests
import os
import pandas as pd
from pprint import pprint
import _pickle
from tqdm import tqdm
import numpy as np
import warnings
import multiprocessing
import collections


output_dir = "/gstock/biolo_datasets/variation/benchmark/Databases/OMIM/JSON_API/"
check_downloaded_omim = [int(e.replace(".pkl", "")) for e in os.listdir(output_dir)]
genes_to_download = list()
genes_already_download = list()
pheno_omim_to_download = list()
pheno_already_download = list()
j = 0

d_stats = collections.defaultdict(int)

def forbid(r):
    forbidden = [':', '.', 'rs', 'ss', 'dbSNP']
    check = True if set([True if check not in r else False for e in forbidden]) == {True} else False
    return check 
        
        
def omim_api(gene, l):
    global j
    global d_stats

    omim_api = "https://api.omim.org/api/entry?mimNumber={}&include={}&format=json&apiKey=fpUv2M2QTmCGhuPexCO-MA"
    # API VERSION
    # 	response = requests.get(omim_api.format(str(gene), "all"))
    # 	json = response.json()

    # LOCAL VERSION
    # if gene not in check_downloaded_omim:
    genes_to_download.append(int(gene))
    json = requests.get(omim_api.format(str(gene), "all")).json()
    _pickle.dump(json, open(output_dir + "{}.pkl".format(gene), "wb"))
    j += 1
    if j % 20 == 0:
        print("Sleep ...")
        time.sleep(5)

    # else:
    #     json = _pickle.load(open(output_dir + "{}.pkl".format(str(gene)), "rb"))
    #     genes_already_download.append(int(gene))
    #     d_stats["Genes"] += 1

    entry = dict(json["omim"]["entryList"][0]["entry"])
    if "geneMap" in entry:
        if "phenotypeMapList" in entry["geneMap"]:
            phenotypes = json["omim"]["entryList"][0]["entry"]["geneMap"]["phenotypeMapList"]
            for pheno in phenotypes:
                if "phenotypeMimNumber" in pheno["phenotypeMap"]:
                    phenotype_mim_number = pheno["phenotypeMap"]["phenotypeMimNumber"]
                    if phenotype_mim_number in check_downloaded_omim:
                        pheno_omim_to_download.append(int(phenotype_mim_number))
                        second_request = requests.get(omim_api.format(str(phenotype_mim_number), "all")).json()
                        _pickle.dump(second_request, open(output_dir + "{}.pkl".format(phenotype_mim_number), "wb"))
                        j += 1
                        if j % 20 == 0:
                            time.sleep(5)
                            print("Sleep ...")
                    else:
                        # second_request = _pickle.load(open(output_dir + "{}.pkl".format(phenotype_mim_number), "rb"))
                        # pheno_already_download.append(int(phenotype_mim_number))
                        d_stats["Phenotypes"] += 1

                    # for key in second_request['omim']['entryList'][0]:
                    #     if 'clinicalSynopsis' in dict(second_request['omim']['entryList'][0]['entry']).keys():

                    #         clinical_synopsis = dict(
                    #             second_request['omim']['entryList'][0]['entry']['clinicalSynopsis'])
                    #         old_format = False
                    #         if 'oldFormat' in clinical_synopsis:
                    #             old_format = True
                    #             clinical_synopsis = clinical_synopsis['oldFormat']

                    #         for tissue, content in dict(clinical_synopsis).items():
                    #             pheno_name = content[:content.find('{')]
                    #             xref = content[content.find(
                    #                 '{')+1:-1].split('} {')
                    #             xref = [e for e in xref if 'HPO' in e]
                    #             if xref:
                    #                 xref = xref[0].split(' ')
                    #                 xref = [e for e in xref if 'HP:' in e]
                    #             l.append(
                    #                 {
                    #                     'Gene_OMIM': gene,
                    #                     'Pheno_OMIM': phenotype_mim_number,
                    #                     'Pheno_Name': pheno_name,
                    #                     'Old_Format': old_format,
                    #                     'ClinicalSynopsis': tissue,
                    #                     'HPO': xref,
                    #                 }
                    #             )


# refseq_x_vcf = pd.read_parquet("/home/weber/PycharmProjects/ExoCarto/data/2_processed/refseq_vcf_pext_phylocsf.parquet")
# disease_genes = refseq_x_vcf.CCRS_Gene.unique()

# full_human_genes = pd.read_csv(
#     "/home/weber/PycharmProjects/ExoCarto/data/2_processed/GRCh37_RefSeq_lite_hgnc.csv.gz", compression="gzip", sep="\t"
# )[["Name", "MIM", "HGNC", "mRNA_nb"]].dropna()
# full_human_genes = full_human_genes.loc[full_human_genes["mRNA_nb"] > 1]
# full_human_genes["MIM"] = full_human_genes["MIM"].astype(int)
# full_human_genes["HGNC"] = full_human_genes["HGNC"].astype(int)
# full_human_genes.loc[full_human_genes["Name"].isin(disease_genes), "STATUS"] = "DISEASE"
# full_human_genes.loc[~full_human_genes["Name"].isin(disease_genes), "STATUS"] = "HEALTHY"

# # genes = full_human_genes.MIM.unique()
# genes = full_human_genes.loc[full_human_genes["STATUS"] == "HEALTHY", "MIM"].unique()

biomart = pd.read_csv("/gstock/EXOTIC/data/OTHERS/biomart_omim.txt.gz", compression="gzip", sep="\t")
biomart = biomart.dropna(subset=["MIM gene accession"])
biomart["MIM gene accession"] = biomart["MIM gene accession"].astype(int)
genes = biomart["MIM gene accession"].unique()
genes = [
    "604917",
    "138210",
    "617094",
    "617094",
    "617094",
    "600937",
    "616661",
    "616661",
    "617619",
    "610661",
    "610661",
    "602676",
    "605993",
    "613768",
    "190315",
    "616830",
    "603518",
    "615533",
    "601023",
    "619109",
]

m = multiprocessing.Manager()
# l = m.list()
l = list()
# parmap.starmap(omim_api, list(zip(genes)), l, pm_pbar=True, pm_processes=20)
for i, gene in tqdm(enumerate(genes)):
    omim_api(gene, l)
# 	if i % 100:

pprint(d_stats)
exit()

print(len(set(genes_to_download)), len(set(pheno_omim_to_download)))
print(len(set(genes_already_download)), len(set(pheno_already_download)))
omim_hpo = pd.DataFrame(list(l))
print(omim_hpo)
# omim_hpo.to_parquet(project + 'data/clean/3_phenotypes/OMIM_HPO_API.parquet')
