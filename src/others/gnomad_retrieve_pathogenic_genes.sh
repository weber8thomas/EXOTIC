for file in 1 10 11 12 13 14 15 16 17 18 19 2 20 21 22 3 4 5 6 7 8 9 X Y; do; 
python ~/PycharmProjects/ExoCarto/src/others/vcf_gnomad_parsing_pathogenic_genes.py ~/PycharmProjects/ExoCarto/data/1_interim/clinvar_hgmd_pathogenic_genes_reviewed.pkl gnomad.exomes.r2.1.1.sites."$file".vcf.bgz | bgzip | tqdm > gnomad.exomes.r2.1.1.sites."$file".pathogenic_reviewed.vcf.gz
done;