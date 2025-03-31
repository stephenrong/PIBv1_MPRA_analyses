#!/bin/sh

tabix -h ../../../../Datasets/variant_effect_predictions/CADD/GRCh37_v1.6/whole_genome_SNVs_inclAnno.tsv.gz -R ../../results/1a-preprocess_PIBv1_MPRA_final/introgressed_variants.bed | bgzip >| ../../results/2f-annotate_variants_CADD_scores/introgressed_variants-whole_genome_SNVs_inclAnno.tsv.gz
