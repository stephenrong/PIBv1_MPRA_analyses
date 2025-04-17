#!/bin/sh

# bgzip -dk ../../results/1a-preprocess_PIBv1_MPRA_pilot/alternate_variants_lift38_chr/adaptive_variants_lift38_chr.bed.gz
awk '{print $1, $2, $3}' ../../results/1a-preprocess_PIBv1_MPRA_pilot/alternate_variants_lift38_chr/adaptive_variants_lift38_chr.bed | awk '$4=$1 "_" $2 "_" $3' >| ../../results/2g-annotate_variants_phyloP_scores/adaptive_variants_lift38_chr_for_phyloP.bed
bigWigAverageOverBed ../../../Datasets/phylo_conservation/Zoonomia_2020_human_phyloP_scores/data_download/241-mammalian-2020v2.phylop-Homo_sapiens.bigWig ../../results/2g-annotate_variants_phyloP_scores/adaptive_variants_lift38_chr_for_phyloP.bed ../../results/2g-annotate_variants_phyloP_scores/adaptive_variants_lift38_chr-241-mammalian-2020v2.tab -bedOut=../../results/2g-annotate_variants_phyloP_scores/adaptive_variants_lift38_chr-241-mammalian-2020v2.bed
# rm ../../results/1a-preprocess_PIBv1_MPRA_pilot/alternate_variants_lift38_chr/adaptive_variants_lift38_chr.bed
