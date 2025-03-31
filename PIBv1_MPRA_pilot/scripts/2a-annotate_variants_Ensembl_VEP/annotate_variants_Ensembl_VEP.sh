#!/bin/sh
conda activate software-ensembl-vep

# using Ensembl VEP cache

# as VCF format
vep -i ../../results/1a-preprocess_PIBv1_MPRA_final/adaptive_variants_clean.vcf.gz --cache --vcf --force_overwrite --assembly GRCh37 --compress_output bgzip -o ../../results/2a-annotate_variants_Ensembl_VEP/adaptive_variants-Ensembl_VEP_nopick.vcf.gz --fields "Uploaded_variation,Location,Allele,Gene,Feature,Feature_type,Consequence,cDNA_position,CDS_position,Protein_position,Amino_acids,Codons,Existing_variation,Extra,IMPACT,SYMBOL,SYMBOL_SOURCE,STRAND,BIOTYPE"
vep -i ../../results/1a-preprocess_PIBv1_MPRA_final/adaptive_variants_clean.vcf.gz --cache --vcf --force_overwrite --pick --assembly GRCh37 --compress_output bgzip -o ../../results/2a-annotate_variants_Ensembl_VEP/adaptive_variants-Ensembl_VEP_pick.vcf.gz --fields "Uploaded_variation,Location,Allele,Gene,Feature,Feature_type,Consequence,cDNA_position,CDS_position,Protein_position,Amino_acids,Codons,Existing_variation,Extra,IMPACT,SYMBOL,SYMBOL_SOURCE,STRAND,BIOTYPE"
