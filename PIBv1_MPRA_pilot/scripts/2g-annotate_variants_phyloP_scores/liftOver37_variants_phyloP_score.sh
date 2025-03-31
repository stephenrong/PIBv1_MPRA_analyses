#!/bin/sh

INPUT="../../results/2g-annotate_variants_phyloP_scores/adaptive_variants_lift38_chr-241-mammalian-2020v2.bed"
CHAIN="../../../../Datasets/reference_genomes/liftOver_chains/hg38ToHg19.over.chain"
OUTPUT="../../results/2g-annotate_variants_phyloP_scores/adaptive_variants_lift37_final-241-mammalian-2020v2.bed"

CrossMap.py bed ${CHAIN} ${INPUT} ${OUTPUT}
