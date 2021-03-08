import sys, os
from pprint import pprint
import pandas as pd
from cyvcf2 import VCF, Writer
from tqdm import tqdm


def prepare_header(vcf, vep_field, vep_separator):
    index_dict = dict()
    if vep_field:
        for h in vcf.header_iter():
            if "ID" in h.info().keys():
                if h.info()["ID"] == vep_field:
                    csq_header = h.info()["Description"].split(vep_separator)
                    for elem in csq_header:
                        index_dict[elem] = csq_header.index(elem)
    return index_dict


def parsing(vcf, w, gene_list, vep_field, vep_separator):
    for j, record in tqdm(enumerate(vcf)):
        # if j == 10000:
        # break
        if int(record.INFO.get("AC")) > 0:
            csq = record.INFO.get(vep_field).split(",")
            check = False
            for case in csq:
                case = case.split(vep_separator)
                if case[index_dict["SYMBOL"]] in gene_list:
                    if "synonymous" in case[index_dict["Consequence"]]:
                        check = True
            if check is True:
                w.write_record(record)


if __name__ == "__main__":
    vcf_path = "/gstock/biolo_datasets/variation/variation_sets/gnomAD/EXOME_MISTIC/RAW/gnomAD_EXOME.vcf.gz"
    vcf_myo_output = "/gstock/biolo_datasets/Sarah/Synonymes/gnomad_myo_synonymous.vcf.gz"
    myo_genes_path = "/gstock/biolo_datasets/Sarah/MuscleGeneTable2020.csv"

    vep_field = "vep"
    vep_separator = "|"

    gene_list = pd.read_csv(myo_genes_path, sep="\t")["Gene symbol"].values.tolist()
    vcf = VCF(vcf_path)
    w = Writer(vcf_myo_output, vcf)

    index_dict = prepare_header(vcf, vep_field, vep_separator)
    parsing(vcf, w, gene_list, vep_field, vep_separator)
