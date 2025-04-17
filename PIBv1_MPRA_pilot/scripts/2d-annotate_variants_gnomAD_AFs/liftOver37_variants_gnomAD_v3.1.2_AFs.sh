#!/bin/sh

INPUT="../../results/2d-annotate_variants_gnomAD_AFs/adaptive_variants_lift38_chr-gnomad.genomes.v3.1.2.sites.vcf.gz"
CHAIN="../../../Datasets/reference_genomes/liftOver_chains/hg38ToHg19.over.chain"
GENOME="../../../Datasets/reference_genomes/Homo_sapiens-hg19-GRCh37/UCSC_hg19/hg19.fa.gz"
OUTPUT="../../results/2d-annotate_variants_gnomAD_AFs/adaptive_variants_lift38_chr_lift37-gnomad.genomes.v3.1.2.sites.vcf"

CrossMap.py vcf ${CHAIN} ${INPUT} ${GENOME} ${OUTPUT}
bgzip ${OUTPUT}
tabix -p vcf ${OUTPUT}.gz
