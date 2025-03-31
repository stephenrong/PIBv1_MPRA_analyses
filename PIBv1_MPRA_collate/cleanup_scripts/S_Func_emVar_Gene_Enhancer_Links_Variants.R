#!/bin/R

library(data.table)
library(tidyverse)

S_Func_emVar_Gene_Enhancer_Links_Variants <- fread("../../PIBv1_MPRA_pilot/results/4-table_emVars_linked_gene_lists/Supplementary_Table_SX-Pairs_of_adaptive_introgressed_MPRA_emVar_variants_and_linked_genes_in_pilot_dataset.txt")
S_Func_emVar_Gene_Enhancer_Links_Variants <- S_Func_emVar_Gene_Enhancer_Links_Variants %>% 
	dplyr::select(-`CADD score`, -`phyloP score`) %>% 
	mutate(`Variant ID` = gsub("chr", "", `Variant ID`)) %>% 
	mutate(`Variant ID` = gsub("-", "_", `Variant ID`))

names(S_Func_emVar_Gene_Enhancer_Links_Variants) <- c("Variant ID", "Gene symbol", "Gene-enhancer count", "Gene-enhancer methods")

fwrite(S_Func_emVar_Gene_Enhancer_Links_Variants, file="../submission_files/S_Func_emVar_Gene_Enhancer_Links_Variants.txt", sep="\t")
