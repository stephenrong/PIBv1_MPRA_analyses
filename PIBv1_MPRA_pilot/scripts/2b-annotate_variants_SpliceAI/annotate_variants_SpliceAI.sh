#!/bin/sh
# conda activate software-spliceai

# SpliceAI raw
# using spliceai precomputed
bgzip -dk ../../results/1a-preprocess_PIBv1_MPRA_pilot/adaptive_variants.bed.gz
tabix -h ../../../Datasets/variant_effect_predictions/SpliceAI/data_download/spliceai_scores.raw.snv.hg19.vcf.gz -R ../../results/1a-preprocess_PIBv1_MPRA_pilot/adaptive_variants.bed | bgzip >| ../../results/2b-annotate_variants_SpliceAI/adaptive_variants-SpliceAI_raw_hg19_temp.vcf.gz
tabix -p vcf ../../results/2b-annotate_variants_SpliceAI/adaptive_variants-SpliceAI_raw_hg19_temp.vcf.gz
bcftools isec -p ../../results/2b-annotate_variants_SpliceAI/SpliceAI_raw_hg19_temp -Oz ../../results/1a-preprocess_PIBv1_MPRA_pilot/adaptive_variants.vcf.gz ../../results/2b-annotate_variants_SpliceAI/adaptive_variants-SpliceAI_raw_hg19_temp.vcf.gz

mv ../../results/2b-annotate_variants_SpliceAI/SpliceAI_raw_hg19_temp/0003.vcf.gz ../../results/2b-annotate_variants_SpliceAI/adaptive_variants-SpliceAI_raw_hg19.vcf.gz
mv ../../results/2b-annotate_variants_SpliceAI/SpliceAI_raw_hg19_temp/0003.vcf.gz.tbi ../../results/2b-annotate_variants_SpliceAI/adaptive_variants-SpliceAI_raw_hg19.vcf.gz.tbi

rm -rf ../../results/2b-annotate_variants_SpliceAI/SpliceAI_raw_hg19_temp/
rm ../../results/2b-annotate_variants_SpliceAI/adaptive_variants-SpliceAI_raw_hg19_temp.vcf.gz
rm ../../results/2b-annotate_variants_SpliceAI/adaptive_variants-SpliceAI_raw_hg19_temp.vcf.gz.tbi
rm ../../results/1a-preprocess_PIBv1_MPRA_pilot/adaptive_variants.bed

# # SpliceAI masked
# # using spliceai precomputed
# bgzip -dk ../../results/1a-preprocess_PIBv1_MPRA_pilot/adaptive_variants.bed.gz
# tabix -h ../../../Datasets/variant_effect_predictions/SpliceAI/data_download/spliceai_scores.masked.snv.hg19.vcf.gz -R ../../results/1a-preprocess_PIBv1_MPRA_pilot/adaptive_variants.bed | bgzip >| ../../results/2b-annotate_variants_SpliceAI/adaptive_variants-SpliceAI_masked_hg19_temp.vcf.gz
# tabix -p vcf ../../results/2b-annotate_variants_SpliceAI/adaptive_variants-SpliceAI_masked_hg19_temp.vcf.gz
# bcftools isec -p ../../results/2b-annotate_variants_SpliceAI/SpliceAI_masked_hg19_temp -Oz ../../results/1a-preprocess_PIBv1_MPRA_pilot/adaptive_variants.vcf.gz ../../results/2b-annotate_variants_SpliceAI/adaptive_variants-SpliceAI_masked_hg19_temp.vcf.gz

# mv ../../results/2b-annotate_variants_SpliceAI/SpliceAI_masked_hg19_temp/0003.vcf.gz ../../results/2b-annotate_variants_SpliceAI/adaptive_variants-SpliceAI_masked_hg19.vcf.gz
# mv ../../results/2b-annotate_variants_SpliceAI/SpliceAI_masked_hg19_temp/0003.vcf.gz.tbi ../../results/2b-annotate_variants_SpliceAI/adaptive_variants-SpliceAI_masked_hg19.vcf.gz.tbi

# rm -rf ../../results/2b-annotate_variants_SpliceAI/SpliceAI_masked_hg19_temp/
# rm ../../results/2b-annotate_variants_SpliceAI/adaptive_variants-SpliceAI_masked_hg19_temp.vcf.gz
# rm ../../results/2b-annotate_variants_SpliceAI/adaptive_variants-SpliceAI_masked_hg19_temp.vcf.gz.tbi
# rm ../../results/1a-preprocess_PIBv1_MPRA_pilot/adaptive_variants.bed
