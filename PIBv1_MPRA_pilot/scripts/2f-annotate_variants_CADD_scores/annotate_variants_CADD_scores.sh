#!/bin/sh

bgzip -dk ../../results/1a-preprocess_PIBv1_MPRA_pilot/adaptive_variants.bed.gz
tabix -h ../../../../Datasets/variant_effect_predictions/CADD/GRCh37_v1.6/whole_genome_SNVs_inclAnno.tsv.gz -R ../../results/1a-preprocess_PIBv1_MPRA_pilot/adaptive_variants.bed | bgzip >| ../../results/2f-annotate_variants_CADD_scores/adaptive_variants-whole_genome_SNVs_inclAnno.tsv.gz
rm ../../results/1a-preprocess_PIBv1_MPRA_pilot/adaptive_variants.bed
