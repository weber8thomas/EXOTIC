import collections
from cyvcf2 import VCF
import multiprocessing
import warnings
import numpy as np
from tqdm import tqdm
import sys

project = "/home/weber/PycharmProjects/EXOTIC/"
sys.path.insert(0, project)
from src.utils import utils
import _pickle
import pandarallel
from pprint import pprint
import pandas as pd
import os
import parmap

warnings.simplefilter(action="ignore", category=FutureWarning)
# default='warn'from pprint import pprint
pd.options.mode.chained_assignment = None
# pandarallel.pandarallel.initialize(nb_workers=os.cpu_count(), progress_bar=True)

# CONFIG
config = utils.load_config_file(project + "src/config.yaml")

pathogenic_genes_file = project + config["ExtractDiseaseGenes"]["DISEASE_GENES_PATH"]

# PHENOTYPES

# hpo = pd.read_parquet(project + "data/2_processed/Ensembl_phenotypes_GRCh37.parquet")
# hpo = hpo.loc[hpo["description"].str.contains("not specified|HGMD") == False]
# hpo["description"] = hpo["description"].str.lower()
# hpo = hpo.drop_duplicates(subset=["GENE", "CHROM", "POSITION", "description"])


# vcf = VCF(config['GeneExonLevel']['clinvar_file'])

# vep_field = 'CSQ'
# vep_separator = '|'
# list_disease_genes = list()
# list_disease_genes_symbol = list()


def get_header(vcf, vep_field, vep_separator):
    index_dict = dict()
    if vep_field:
        for h in vcf.header_iter():
            try:
                if h.info()["ID"] == vep_field:
                    csq_header = h.info()["Description"].split(vep_separator)
                    for elem in csq_header:
                        index_dict[elem] = csq_header.index(elem)
            except:
                pass
    return index_dict


def prepare_dict_variants_pathogenic(
    vcf_file,
):
    vep_consequences = {
        "HIGH": 3,
        "MODERATE": 2,
        "LOW": 1,
        "MODIFIER": 0,
    }

    clinvar_review_status = {
        "practice_guideline": 4,
        "reviewed_by_expert_panel": 3,
        "criteria_provided,_multiple_submitters,_no_conflicts": 2,
        "criteria_provided,_conflicting_interpretations": 1,
        "criteria_provided,_single_submitter": 1,
        "no_assertion_for_the_individual_variant": 0,
        "no_assertion_criteria_provided": 0,
        "no_assertion_provided": 0,
    }

    vcf = VCF(vcf_file)
    vep_field = "CSQ"
    vep_separator = "|"

    # INDEX DICT
    # index_dict = get_header(vcf, vep_field="CSQ", vep_separator="|")

    # HPO DICT
    hpo_dict = list()

    # CLNSIG FOR EACH DB
    db_clnsig = {
        "HGMD": "DM",
        "ClinVar": "athogenic",
    }

    # INSTANCIATION
    d_variants = collections.defaultdict(dict)
    i = 0

    # LOOP ON VARIANTS
    for counter, variant in enumerate(tqdm(vcf)):

        # BOOLEAN CHECKER FOR REVIEWING STATUS
        check = False

        # UNCOMMENT IF TEST
        # if counter == 10000:
        # break

        # SNV ONLY
        if len(variant.REF) <= 50 and len(variant.ALT[0]) <= 50:

            # CLINVAR CLNSIG
            if "CLNSIG=" in str(variant):
                cln_sig = variant.INFO.get("CLNSIG")
                db = "ClinVar"

            # HGMD CLASS
            elif "CLASS=" in str(variant):
                cln_sig = variant.INFO.get("CLASS")
                db = "HGMD"

                # ALL HGMD VARIANTS ARE OK FOR REVIEWING STATUS

            try:
                # IF CLNSIG IS PRESENT
                if cln_sig:

                    # IF CLNSIG OK FOR SELECTED DB
                    if db == "HGMD":
                        if db_clnsig[db] == cln_sig:
                            check = True
                            # check = False
                    elif db == "ClinVar":
                        if (("athogenic" in cln_sig) or ("enign" in cln_sig)) and "onflict" not in cln_sig:

                            # CLNREVSTAT
                            if variant.INFO["CLNREVSTAT"]:
                                tmp_stats = variant.INFO["CLNREVSTAT"]
                                if clinvar_review_status[tmp_stats] > 0 and "conflicting" not in tmp_stats:
                                    check = True

                # GOOD REVIEW STATUS
                if check is True:

                    # VEP CSQ FIELD
                    csq = variant.INFO.get(vep_field)
                    if "," in csq:
                        csq = csq.split(",")
                    else:
                        csq = [csq]

                    # RETRIEVE MAX IMPACT FOR ALL VEP CASES (EXAMPLE : ONE VARIANT OVERLAP ON TWO GENES)
                    max_impact = max([vep_consequences[case.split(vep_separator)[index_dict["IMPACT"]]] for case in csq])

                    # LOOP ON CASES
                    for case in csq:
                        case = case.split(vep_separator)

                        # IF CASE == MAX IMPACT
                        if vep_consequences[case[index_dict["IMPACT"]]] == max_impact and max_impact > 0:
                            hgnc = case[index_dict['HGNC_ID"']]

                            # IF HGNC IS PRESENT
                            if hgnc:
                                hgnc = int(hgnc)
                                if "HP:" in str(variant):
                                    i += 1
                                    clndisdb = variant.INFO.get("CLNDISDB")
                                    if clndisdb:
                                        hpo = clndisdb.split("|")
                                        for hpo_case in hpo:

                                            for sub_hpo_case in hpo_case.split(","):

                                                if "enign" in cln_sig:
                                                    status = "Benign"

                                                elif "athogenic" in cln_sig:
                                                    status = "Pathogenic"

                                                hpo_dict.append(
                                                    {
                                                        "Status": status,
                                                        "VAR_ID": str(variant.CHROM)
                                                        + "_"
                                                        + str(variant.POS)
                                                        + "_"
                                                        + str(variant.REF)
                                                        + "_"
                                                        + str(variant.ALT[0]),
                                                        "Xref": sub_hpo_case,
                                                        "HGNC": hgnc,
                                                        "CLNREVSTAT": clinvar_review_status[tmp_stats],
                                                    }
                                                )

                                    # for consequence in ['splice_donor', 'splice_acceptor', 'missense', 'start', 'stop', 'intron']:
                                    # if consequence in case[index_dict['Consequence']]:
                                    #                       # ADD DICT TO LIST WITH ALL VARIANT INFO
                                    # d_variants[int(hgnc)][str(variant.POS) + '_' + str(variant.REF) + '_' + str(variant.ALT[0])] = {'VAR_ID': str(
                                    # variant.POS) + '_' + str(variant.REF) + '_' + str(variant.ALT[0]), 'Source': db, 'STATUS': "Pathogenic", 'MC': consequence}
            except KeyError:
                pass
    print(i)
    return hpo_dict


def convert_clinvar_to_table(
    vcf_file,
):
    vep_consequences = {
        "HIGH": 3,
        "MODERATE": 2,
        "LOW": 1,
        "MODIFIER": 0,
    }

    clinvar_review_status = {
        "practice_guideline": 4,
        "reviewed_by_expert_panel": 3,
        "criteria_provided,_multiple_submitters,_no_conflicts": 2,
        "criteria_provided,_conflicting_interpretations": 1,
        "criteria_provided,_single_submitter": 1,
        "no_assertion_for_the_individual_variant": 0,
        "no_assertion_criteria_provided": 0,
        "no_assertion_provided": 0,
    }

    vcf = VCF(vcf_file)
    vep_field = "CSQ"
    vep_separator = "|"

    # INDEX DICT
    # index_dict = get_header(vcf, vep_field="CSQ", vep_separator="|")

    # LIST
    return_list = list()

    # CLNSIG FOR EACH DB
    db_clnsig = {
        "HGMD": "DM",
        "ClinVar": "athogenic",
    }

    # INSTANCIATION
    d_variants = collections.defaultdict(dict)
    i = 0

    # LOOP ON VARIANTS
    for counter, variant in enumerate(tqdm(vcf)):

        # BOOLEAN CHECKER FOR REVIEWING STATUS
        check = False

        # UNCOMMENT IF TEST
        # if counter == 10000:
        # break

        # SNV ONLY
        if len(variant.REF) == 1 and len(variant.ALT[0]) == 1:

            # CLINVAR CLNSIG
            if "CLNSIG=" in str(variant):
                cln_sig = variant.INFO.get("CLNSIG")
                db = "ClinVar"

            # HGMD CLASS
            elif "CLASS=" in str(variant):
                cln_sig = variant.INFO.get("CLASS")
                db = "HGMD"

                # ALL HGMD VARIANTS ARE OK FOR REVIEWING STATUS

            try:
                # IF CLNSIG IS PRESENT
                if cln_sig:

                    # IF CLNSIG OK FOR SELECTED DB
                    if db == "ClinVar":
                        if (("athogenic" in cln_sig) or ("enign" in cln_sig)) and "onflict" not in cln_sig:

                            # CLNREVSTAT
                            if variant.INFO["CLNREVSTAT"]:
                                tmp_stats = variant.INFO["CLNREVSTAT"]
                                if clinvar_review_status[tmp_stats] >= 0 and "conflicting" not in tmp_stats:
                                    check = True

                    # GOOD REVIEW STATUS
                    if check is True:

                        # VEP CSQ FIELD
                        # csq = variant.INFO.get(vep_field)
                        # if ',' in csq:
                        # csq = csq.split(',')
                        # else:
                        # csq = [csq]

                        # RETRIEVE MAX IMPACT FOR ALL VEP CASES (EXAMPLE : ONE VARIANT OVERLAP ON TWO GENES)
                        # max_impact = max([vep_consequences[case.split(vep_separator)[index_dict['IMPACT']]] for case in csq])

                        # LOOP ON CASES
                        # for case in csq:
                        # case = case.split(vep_separator)

                        # IF CASE == MAX IMPACT
                        # if vep_consequences[case[index_dict['IMPACT']]] == max_impact and max_impact > 0:
                        # hgnc = case[index_dict['HGNC_ID"']]

                        rs = variant.INFO.get("RS")

                        if "enign" in cln_sig:
                            status = "Benign"
                        #
                        elif "athogenic" in cln_sig:
                            status = "Pathogenic"

                        return_list.append(
                            {
                                "VAR_ID": str(variant.CHROM) + "_" + str(variant.POS) + "_" + str(variant.REF) + "_" + str(variant.ALT[0]),
                                # 'CSQ' : case[index_dict['Consequence']],
                                "MC": variant.MC.split("|")[1],
                                "Status": status,
                                "Real_Status": cln_sig,
                                # 'HGNC': hgnc,
                                "RS_STARS": clinvar_review_status[tmp_stats],
                                "CLNREVSTAT": tmp_stats,
                                "rs": rs,
                            }
                        )

                        # # IF HGNC IS PRESENT
                        # if hgnc:
                        #     hgnc = int(hgnc)
                        #     if 'HP:' in str(variant):
                        #         i += 1
                        #         clndisdb = variant.INFO.get('CLNDISDB')
                        #         if clndisdb:
                        #             hpo = clndisdb.split('|')
                        #             for hpo_case in hpo:

                        #                 for sub_hpo_case in hpo_case.split(','):

                        #                     if 'enign' in cln_sig:
                        #                         status = 'Benign'

                        #                     elif 'athogenic' in cln_sig:
                        #                         status = 'Pathogenic'

                        #                     return_list.append(
                        #                         {
                        #                             'Status': status,
                        #                             'VAR_ID': str(variant.CHROM) + '_' + str(variant.POS) + '_' + str(variant.REF) + '_' + str(variant.ALT[0]),
                        #                             'Xref': sub_hpo_case,
                        #                             'HGNC': hgnc,
                        #                             'CLNREVSTAT': clinvar_review_status[tmp_stats],
                        #                         }
                        #                     )

                        # for consequence in ['splice_donor', 'splice_acceptor', 'missense', 'start', 'stop', 'intron']:
                        # if consequence in case[index_dict['Consequence']]:
                        #                    # ADD DICT TO LIST WITH ALL VARIANT INFO
                        # d_variants[int(hgnc)][str(variant.POS) + '_' + str(variant.REF) + '_' + str(variant.ALT[0])] = {'VAR_ID': str(
                        # variant.POS) + '_' + str(variant.REF) + '_' + str(variant.ALT[0]), 'Source': db, 'STATUS': "Pathogenic", 'MC': consequence}
            except KeyError:
                pass
    print(i)
    return return_list


def convert_clinvar_to_table_wt_vep(vcf_file, d_stats):

    clinvar_review_status = {
        "practice_guideline": 4,
        "reviewed_by_expert_panel": 3,
        "criteria_provided,_multiple_submitters,_no_conflicts": 2,
        "criteria_provided,_conflicting_interpretations": 1,
        "criteria_provided,_single_submitter": 1,
        "no_assertion_for_the_individual_variant": 0,
        "no_assertion_criteria_provided": 0,
        "no_assertion_provided": 0,
    }

    vcf = VCF(vcf_file)

    # LIST
    return_list = list()

    # CLNSIG FOR EACH DB
    db_clnsig = {
        "HGMD": "DM",
        "ClinVar": "athogenic",
    }

    # INSTANCIATION
    d_variants = collections.defaultdict(dict)
    i = 0

    # LOOP ON VARIANTS
    for counter, variant in enumerate(tqdm(vcf)):
        d_stats["Total"] += 1
        # print(variant)

        # BOOLEAN CHECKER FOR REVIEWING STATUS
        check = False

        # UNCOMMENT IF TEST
        # if counter == 10000:
        # break

        # SNV ONLY
        #         if variant.REF and variant.ALT and len(variant.REF) == 1 and len(variant.ALT[0]) == 1:

        # CLINVAR CLNSIG
        if "CLNSIG=" in str(variant):
            cln_sig = variant.INFO.get("CLNSIG")
            db = "ClinVar"

        # HGMD CLASS
        elif "CLASS=" in str(variant):
            cln_sig = variant.INFO.get("CLASS")
            db = "HGMD"

        # ALL HGMD VARIANTS ARE OK FOR REVIEWING STATUS

        try:
            # IF CLNSIG IS PRESENT
            if cln_sig:
                d_stats["CLNSIG"] += 1

                # IF CLNSIG OK FOR SELECTED DB
                if db == "ClinVar":
                    check = True
                    #                                       if (('athogenic' in cln_sig) or ('enign' in cln_sig)) and 'onflict' not in cln_sig:

                    #                                               # CLNREVSTAT
                    if variant.INFO["CLNREVSTAT"]:
                        tmp_stats = variant.INFO["CLNREVSTAT"]
                    else:
                        tmp_stats = "None"
                #                                                       if clinvar_review_status[tmp_stats] >= 0 and 'conflicting' not in tmp_stats:
                #                                                               check = True

                # GOOD REVIEW STATUS
                if check is True:

                    rs = variant.INFO.get("RS")

                    clnvi = np.nan
                    if variant.INFO.get("CLNVI"):
                        if "OMIM_Allelic_Variant" in variant.INFO.get("CLNVI"):
                            clnvi = [
                                e.replace("OMIM_Allelic_Variant:", "") for e in variant.INFO.get("CLNVI").split("|") if "OMIM_Allelic_Variant" in e
                            ][0]

                    if variant.INFO.get("MC") and variant.INFO.get("GENEINFO"):
                        d_stats["MC_GENEINFO"] += 1

                        alleleid = variant.INFO.get("ALLELEID")

                        mc = variant.INFO.get("MC").split("|")[1]
                        if variant.INFO.get("CLNDISDB"):
                            clndisdb = variant.INFO.get("CLNDISDB").split("|")
                        else:
                            clndisdb = []
                        hpo = [sub_e.replace("Human_Phenotype_Ontology:", "") for e in clndisdb for sub_e in e.split(",") if "HP" in sub_e]
                        geneinfo = variant.INFO.get("GENEINFO").split(":")[0]

                        if "enign" in cln_sig:
                            status = "Benign"
                        #
                        elif "athogenic" in cln_sig:
                            status = "Pathogenic"

                        else:
                            status = "Other"

                        d_stats["FINAL"] += 1
                        return_list.append(
                            {
                                "GENE": geneinfo,
                                "VAR_ID": str(variant.CHROM) + "_" + str(variant.POS) + "_" + str(variant.REF) + "_" + str(variant.ALT),
                                "CHROM": str(variant.CHROM),
                                "POS": variant.POS,
                                "REF": str(variant.REF),
                                "ALT": str(variant.ALT),
                                "MC": mc,
                                "Status": status,
                                "Real_Status": cln_sig,
                                "RS_STARS": clinvar_review_status[tmp_stats],
                                "CLNREVSTAT": tmp_stats,
                                "HPO": hpo,
                                "rs": rs,
                                "alleleid": alleleid,
                                "OMIM_VARIANT_ID": clnvi,
                            }
                        )
        except KeyError:
            pass
    print(i)
    return return_list


d_stats = collections.defaultdict(int)
return_list = convert_clinvar_to_table_wt_vep(config["ClinVar"]["Jan_21"], d_stats)

# print(config['GeneExonLevel']['clinvar_file'], type(config['GeneExonLevel']['clinvar_file']))
df = pd.DataFrame(return_list)
# df["ALT"] = df["ALT"].str.replace("['", "")
# df["ALT"] = df["ALT"].str.replace("']", "")
print(df)


df.to_parquet(config["EXOTIC"]["clinvar_file_path"])
