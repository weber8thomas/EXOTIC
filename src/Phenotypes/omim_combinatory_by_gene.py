import time
import requests
import os
import pandas as pd
from pprint import pprint
import _pickle
from tqdm import tqdm
import numpy as np
import parmap
import warnings
import multiprocessing
import collections
import sys

project = "/home/weber/PycharmProjects/ExoCarto/"
sys.path.insert(0, project)


output_dir = "/gstock/biolo_datasets/variation/benchmark/Databases/OMIM/JSON_API/"


def retrieve_clinical_synopsis_omim(gene, l):
    omim_clinical_synopsis_dict = {
        "abdomen": {
            "abdomenBiliaryTract": None,
            "abdomenExternalFeatures": None,
            "abdomenGastrointestinal": None,
            "abdomenLiver": None,
            "abdomenPancreas": None,
            "abdomenSpleen": None,
        },
        "cardiovascular": {"cardiovascularHeart": None, "cardiovascularVascular": None},
        "chest": {
            "chestBreasts": None,
            "chestDiaphragm": None,
            "chestExternalFeatures": None,
            "chestRibsSternumClaviclesAndScapulae": None,
        },
        "endocrineFeatures": None,
        "genitourinary": {
            "genitourinaryBladder": None,
            "genitourinaryExternalGenitaliaFemale": None,
            "genitourinaryExternalGenitaliaMale": None,
            "genitourinaryInternalGenitaliaFemale": None,
            "genitourinaryInternalGenitaliaMale": None,
            "genitourinaryKidneys": None,
            "genitourinaryUreters": None,
        },
        "growth": {"growthHeight": None, "growthOther": None, "growthWeight": None},
        "headAndNeck": {
            "headAndNeckEars": None,
            "headAndNeckEyes": None,
            "headAndNeckFace": None,
            "headAndNeckHead": None,
            "headAndNeckMouth": None,
            "headAndNeckNeck": None,
            "headAndNeckNose": None,
            "headAndNeckTeeth": None,
        },
        "hematology": None,
        "immunology": None,
        "inheritance": None,
        "laboratoryAbnormalities": None,
        "metabolicFeatures": None,
        "miscellaneous": None,
        "molecularBasis": None,
        "muscleSoftTissue": None,
        "neoplasia": None,
        "neurologic": {
            "neurologicBehavioralPsychiatricManifestations": None,
            "neurologicCentralNervousSystem": None,
            "neurologicPeripheralNervousSystem": None,
        },
        "prenatalManifestations": {
            "prenatalManifestationsAmnioticFluid": None,
            "prenatalManifestationsDelivery": None,
            "prenatalManifestationsMaternal": None,
            "prenatalManifestationsMovement": None,
            "prenatalManifestationsPlacentaAndUmbilicalCord": None,
        },
        "respiratory": {
            "respiratoryAirways": None,
            "respiratoryLarynx": None,
            "respiratoryLung": None,
            "respiratoryNasopharynx": None,
        },
        "skeletal": {
            "skeletalFeet": None,
            "skeletalHands": None,
            "skeletalLimbs": None,
            "skeletalPelvis": None,
            "skeletalSkull": None,
            "skeletalSpine": None,
        },
        "skinNailsHair": {
            "skinNailsHairHair": None,
            "skinNailsHairNails": None,
            "skinNailsHairSkin": None,
            "skinNailsHairSkinElectronMicroscopy": None,
            "skinNailsHairSkinHistology": None,
        },
        "voice": None,
    }

    omim_global_dict = {
        "abdomenBiliaryTract": "abdomen",
        "abdomenExternalFeatures": "abdomen",
        "abdomenGastrointestinal": "abdomen",
        "abdomenLiver": "abdomen",
        "abdomenPancreas": "abdomen",
        "abdomenSpleen": "abdomen",
        "cardiovascularHeart": "cardiovascular",
        "cardiovascularVascular": "cardiovascular",
        "chestBreasts": "chest",
        "chestDiaphragm": "chest",
        "chestExternalFeatures": "chest",
        "chestRibsSternumClaviclesAndScapulae": "chest",
        "endocrineFeatures": "endocrineFeatures",
        "genitourinaryBladder": "genitourinary",
        "genitourinaryExternalGenitaliaFemale": "genitourinary",
        "genitourinaryExternalGenitaliaMale": "genitourinary",
        "genitourinaryInternalGenitaliaFemale": "genitourinary",
        "genitourinaryInternalGenitaliaMale": "genitourinary",
        "genitourinaryKidneys": "genitourinary",
        "genitourinaryUreters": "genitourinary",
        "growthHeight": "growth",
        "growthOther": "growth",
        "growthWeight": "growth",
        "headAndNeckEars": "headAndNeck",
        "headAndNeckEyes": "headAndNeck",
        "headAndNeckFace": "headAndNeck",
        "headAndNeckHead": "headAndNeck",
        "headAndNeckMouth": "headAndNeck",
        "headAndNeckNeck": "headAndNeck",
        "headAndNeckNose": "headAndNeck",
        "headAndNeckTeeth": "headAndNeck",
        "hematology": "hematology",
        "immunology": "immunology",
        "inheritance": "inheritance",
        "laboratoryAbnormalities": "laboratoryAbnormalities",
        "metabolicFeatures": "metabolicFeatures",
        "miscellaneous": "miscellaneous",
        "molecularBasis": "molecularBasis",
        "muscleSoftTissue": "muscleSoftTissue",
        "neoplasia": "neoplasia",
        "neurologicBehavioralPsychiatricManifestations": "neurologic",
        "neurologicCentralNervousSystem": "neurologic",
        "neurologicPeripheralNervousSystem": "neurologic",
        "prenatalManifestationsAmnioticFluid": "prenatalManifestations",
        "prenatalManifestationsDelivery": "prenatalManifestations",
        "prenatalManifestationsMaternal": "prenatalManifestations",
        "prenatalManifestationsMovement": "prenatalManifestations",
        "prenatalManifestationsPlacentaAndUmbilicalCord": "prenatalManifestations",
        "respiratoryAirways": "respiratory",
        "respiratoryLarynx": "respiratory",
        "respiratoryLung": "respiratory",
        "respiratoryNasopharynx": "respiratory",
        "skeletalFeet": "skeletal",
        "skeletalHands": "skeletal",
        "skeletalLimbs": "skeletal",
        "skeletalPelvis": "skeletal",
        "skeletalSkull": "skeletal",
        "skeletalSpine": "skeletal",
        "skinNailsHairHair": "skinNailsHair",
        "skinNailsHairNails": "skinNailsHair",
        "skinNailsHairSkin": "skinNailsHair",
        "skinNailsHairSkinElectronMicroscopy": "skinNailsHair",
        "skinNailsHairSkinHistology": "skinNailsHair",
        "voice": "voice",
        "Abdomen": "abdomen",
        "Cardiac": "cardiovascular",
        "Cardiovascular": "cardiovascular",
        "Cranium": "headAndNeck",
        "Ears": "headAndNeck",
        "Endo": "endocrineFeatures",
        "Endocrine": "endocrineFeatures",
        "Eye": "headAndNeck",
        "Eyes": "headAndNeck",
        "Facies": "headAndNeck",
        "Growth": "growth",
        "Hair": "skinNailsHair",
        "Head": "headAndNeck",
        "Heme": "hematology",
        "Immunol": "immunology",
        "Immunology": "immunology",
        "Inheritance": "inheritance",
        "Joints": "skeletal",
        "Limbs": "skeletal",
        "Liver": "abdomen",
        "Mandible": "headAndNeck",
        "Metabolic": "metabolicFeatures",
        "Misc": "miscellaneous",
        "Mouth": "headAndNeck",
        "Muscle": "muscleSoftTissue",
        "Nails": "skinNailsHair",
        "Neck": "headAndNeck",
        "Neuro": "neurologic",
        "Nose": "headAndNeck",
        "Oncology": "neoplasia",
        "Pulm": "respiratory",
        "Pulmonary": "respiratory",
        "Radiology": "laboratoryAbnormalities",
        "Resp": "respiratory",
        "Respiratory": "respiratory",
        "Skel": "skeletal",
        "Skeletal": "skeletal",
        "Skin": "skinNailsHair",
        "Skull": "skeletal",
        "Spine": "neurologic",
        "Teeth": "headAndNeck",
        "Thorax": "chest",
        "Tongue": "headAndNeck",
        "Vascular": "cardiovascular",
        "Voice": "headAndNeck",
    }

    json = _pickle.load(open(output_dir + "{}.pkl".format(str(gene)), "rb"))

    entry = dict(json["omim"]["entryList"][0]["entry"])
    if "geneMap" in entry:
        if "phenotypeMapList" in entry["geneMap"]:
            phenotypes = json["omim"]["entryList"][0]["entry"]["geneMap"]["phenotypeMapList"]
            for pheno in phenotypes:
                # print(gene, pheno)

                if "phenotypeMimNumber" in pheno["phenotypeMap"]:
                    phenotype_mim_number = pheno["phenotypeMap"]["phenotypeMimNumber"]
                    second_request = _pickle.load(open(output_dir + "{}.pkl".format(phenotype_mim_number), "rb"))

                    for key in second_request["omim"]["entryList"][0]:
                        # print(key)
                        prefered_title = second_request["omim"]["entryList"][0]["entry"]["titles"]["preferredTitle"]

                        if "clinicalSynopsis" in dict(second_request["omim"]["entryList"][0]["entry"]).keys():
                            clinical_synopsis = dict(second_request["omim"]["entryList"][0]["entry"]["clinicalSynopsis"])

                            # # FIRST LEVEL - GLOBAL CLINICAL SYNOPSIS BODY PART OMIM
                            new_dict = collections.defaultdict(list)
                            for k, v in clinical_synopsis.items():
                                # print(k, v)
                                # if "Exists" not in k and k in omim_global_dict:
                                #     new_dict[omim_global_dict[k]].append(v)

                                if "Exists" not in k and k in omim_global_dict:
                                    new_dict[k].append(v)

                                # # OLD FORMAT
                                # if "oldFormat" in clinical_synopsis:
                                #     for k in clinical_synopsis["oldFormat"]:
                                #         if "Exists" not in k and k in omim_global_dict:
                                #             new_dict[omim_global_dict[k]].append(v)

                            new_dict["OMIM"] = gene
                            new_dict["Pheno_OMIM"] = phenotype_mim_number
                            new_dict["Pheno_Name"] = pheno["phenotypeMap"]["phenotype"]
                            new_dict["Pheno_prefered_title"] = prefered_title

                            if "referenceList" in dict(second_request["omim"]["entryList"][0]["entry"]).keys():
                                pmid_list = [
                                    ref["reference"]["pubmedID"]
                                    for ref in second_request["omim"]["entryList"][0]["entry"]["referenceList"]
                                    if "pubmedID" in ref["reference"]
                                ]
                                new_dict["PMID"] = pmid_list

                            l.append(new_dict)


clinvar = pd.read_parquet("/home/weber/PycharmProjects/ExoCarto/data/clean/clinvar_20201010_lite_table.parquet")
disease_genes = clinvar.loc[
    (clinvar["Status"] == "Pathogenic") & (~clinvar["Real_Status"].str.contains("onflict")), "GENE"
].unique()

full_human_genes = pd.read_csv(
    "/home/weber/PycharmProjects/ExoCarto/data/2_processed/GRCh37_RefSeq_lite_hgnc.csv.gz",
    compression="gzip",
    sep="\t",
)[["Name", "MIM", "HGNC", "mRNA_nb"]].dropna()

full_human_genes = full_human_genes.loc[full_human_genes["mRNA_nb"] > 1]
full_human_genes["MIM"] = full_human_genes["MIM"].astype(int)
full_human_genes["HGNC"] = full_human_genes["HGNC"].astype(int)
full_human_genes.loc[full_human_genes["Name"].isin(disease_genes), "STATUS"] = "DISEASE"
full_human_genes.loc[~full_human_genes["Name"].isin(disease_genes), "STATUS"] = "HEALTHY"

genes = full_human_genes["MIM"].unique()

m = multiprocessing.Manager()
l = m.list()

parmap.starmap(retrieve_clinical_synopsis_omim, list(zip(genes)), l, pm_pbar=True)
# for i, gene in tqdm(enumerate(genes)):
# retrieve_clinical_synopsis_omim(gene, l)
# print(old)
# print(set(old_format))
# exit()

df = pd.DataFrame(list(l))
df = df[
    ["OMIM", "Pheno_OMIM", "Pheno_Name", "Pheno_prefered_title", "PMID"]
    + [c for c in df.columns if c not in ["OMIM", "Pheno_OMIM", "Pheno_Name", "Pheno_prefered_title", "PMID"]]
]
print(df)
print(df.OMIM.nunique())

df.to_csv(
    "/home/weber/PycharmProjects/ExoCarto/data/clean/3_phenotypes/omim_genes_bodyparts_detailed.csv.gz",
    compression="gzip",
    sep="\t",
    index=False,
)

# for i, gene in tqdm(enumerate(genes[:3])):
# retrieve_clinical_synopsis_omim(gene)
