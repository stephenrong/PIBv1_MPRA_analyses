#!/bin/sh

bgzip -dk ../../results/1a-preprocess_PIBv1_MPRA_pilot/adaptive_variants.bed.gz
tabix -h ../../../../Datasets/human_popgen_genomes/1KGP_phase3_genomes/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5c.20130502.sites.vcf.gz -R ../../results/1a-preprocess_PIBv1_MPRA_pilot/adaptive_variants.bed | bcftools view -i 'TYPE="SNP"' | bgzip >| ../../results/2c-annotate_variants_1KGP_AFs/adaptive_variants-ALL.wgs.phase3_shapeit2_mvncall_integrated_v5c.20130502.sites.vcf.gz
tabix -p vcf ../../results/2c-annotate_variants_1KGP_AFs/adaptive_variants-ALL.wgs.phase3_shapeit2_mvncall_integrated_v5c.20130502.sites.vcf.gz
rm ../../results/1a-preprocess_PIBv1_MPRA_pilot/adaptive_variants.bed
