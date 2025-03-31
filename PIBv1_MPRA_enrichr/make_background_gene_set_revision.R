#!/bin/R

# Revision uses background gene set for Enricher analysis

# Load libraries
library(tidyverse)
library(data.table)

# Load introgressed functional annotations
introgressed_func_anno <- as_tibble(fread("../../../TucciMPRAFinal/results/3-merge_all_variant_annotations/introgressed_variants_all_annotations_summ.txt.gz"))

# Get nearest and cre-linked gene lists
background_gene_set_expansive <- unique(c(
	unlist(strsplit(introgressed_func_anno$nearest_gene_summ_Symbol, ",")),
	unlist(strsplit(introgressed_func_anno$ABC_summ_Symbol, ",")), 
	unlist(strsplit(introgressed_func_anno$Roadmap_Epigenomics_summ_Symbol, ",")), 
	unlist(strsplit(introgressed_func_anno$EpiMap_summ_Symbol, ",")))
)

background_gene_set_expansive <- background_gene_set_expansive[!is.na(background_gene_set_expansive)]
write(background_gene_set_expansive, "Gene_set_enrichment_background_expansive_revision/Gene_set_enrichment_background_expansive_revision.txt")
