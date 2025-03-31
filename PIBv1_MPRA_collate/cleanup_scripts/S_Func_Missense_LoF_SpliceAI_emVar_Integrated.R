#!/bin/R

library(data.table)
library(tidyverse)

S_Func_Missense_Integrated <- fread("../../PIBv1_MPRA_final/results/4-table_consequence_gene_lists/Supplementary_Table_SX-Corehaps_introgressed_missense_variants_in_final_dataset.txt")
S_Func_Missense_Integrated <- S_Func_Missense_Integrated %>% 
	dplyr::select(-`CADD score`, -`phyloP score`) %>% 
	mutate(`Variant ID` = gsub("chr", "", `Variant ID`)) %>% 
	mutate(`Variant ID` = gsub("-", "_", `Variant ID`))

names(S_Func_Missense_Integrated) <- c("Variant ID", "Ensembl VEP", "Gene symbol")

fwrite(S_Func_Missense_Integrated, file="../submission_files/S_Func_Missense_Integrated.txt", sep="\t")


S_Func_SpliceAI_Integrated <- fread("../../PIBv1_MPRA_final/results/4-table_consequence_gene_lists/Supplementary_Table_SX-Corehaps_introgressed_SpliceAI_variants_in_final_dataset.txt")
S_Func_SpliceAI_Integrated <- S_Func_SpliceAI_Integrated %>% 
	dplyr::select(-`CADD score`, -`phyloP score`) %>% 
	mutate(`Variant ID` = gsub("chr", "", `Variant ID`)) %>% 
	mutate(`Variant ID` = gsub("-", "_", `Variant ID`))

fwrite(S_Func_SpliceAI_Integrated, file="../submission_files/S_Func_SpliceAI_Integrated.txt", sep="\t")


S_Func_LoF_Integrated <- fread("../../PIBv1_MPRA_final/results/4-table_consequence_gene_lists/Supplementary_Table_SX-Corehaps_introgressed_LoF_variants_in_final_dataset.txt")
S_Func_LoF_Integrated <- S_Func_LoF_Integrated %>% 
	dplyr::select(-`CADD score`, -`phyloP score`) %>% 
	mutate(`Variant ID` = gsub("chr", "", `Variant ID`)) %>% 
	mutate(`Variant ID` = gsub("-", "_", `Variant ID`))

fwrite(S_Func_LoF_Integrated, file="../submission_files/S_Func_LoF_Integrated.txt", sep="\t")


S_Func_emVar_Integrated <- fread("../../PIBv1_MPRA_final/results/4-table_emVars_linked_gene_lists/Supplementary_Table_SX-Pairs_of_corehaps_introgressed_MPRA_emVar_variants_and_linked_genes_in_final_dataset.txt")
S_Func_emVar_Integrated <- S_Func_emVar_Integrated %>% 
	dplyr::select(-`CADD score`, -`phyloP score`) %>% 
	mutate(`Variant ID` = gsub("chr", "", `Variant ID`)) %>% 
	mutate(`Variant ID` = gsub("-", "_", `Variant ID`))

names(S_Func_emVar_Integrated) <- c("Variant ID", "Gene symbol", "Gene-enhancer count", "Gene-enhancer methods")

fwrite(S_Func_emVar_Integrated, file="../submission_files/S_Func_emVar_Integrated.txt", sep="\t")
