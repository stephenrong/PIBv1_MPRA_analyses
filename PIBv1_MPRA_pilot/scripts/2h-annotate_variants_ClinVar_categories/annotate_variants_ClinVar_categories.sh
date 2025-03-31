#!/bin/sh

# bgzip -dk ../../results/1a-preprocess_PIBv1_MPRA_pilot/adaptive_variants.bed.gz
tabix -h ../../../../Datasets/disease_annotations/ClinVar_20230326/data_download/GRCh37/clinvar_20230326.vcf.gz -R ../../results/1a-preprocess_PIBv1_MPRA_pilot/adaptive_variants.bed.gz | bcftools view -i 'TYPE="SNP"' | bgzip >| ../../results/2h-annotate_variants_ClinVar_categories/adaptive_variants-clinvar_20230326.vcf.gz
tabix -p vcf ../../results/2h-annotate_variants_ClinVar_categories/adaptive_variants-clinvar_20230326.vcf.gz
# rm ../../results/1a-preprocess_PIBv1_MPRA_pilot/adaptive_variants.bed
