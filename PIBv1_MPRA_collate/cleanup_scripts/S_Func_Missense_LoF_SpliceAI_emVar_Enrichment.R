#!/bin/R

library(data.table)
library(tidyverse)

S_Func_Missense_LoF_SpliceAI_emVar_Enrichment <- fread("../../PIBv1_MPRA_enrichr/Supplementary_tables_and_figures_revision/Supplementary_X-Gene_set_enrichment-combined_figures-Reactome_Pathways_2024_table.txtSupplementary_X-Gene_set_enrichment-combined_figures-Reactome_2022_table.txt")
S_Func_Missense_LoF_SpliceAI_emVar_Enrichment <- S_Func_Missense_LoF_SpliceAI_emVar_Enrichment %>% 
	dplyr::select(`Set`, `Term`, `Overlap`, `P-value`, `Adjusted P-value`, `Odds Ratio`, `Combined Score`, `Genes`)

fwrite(S_Func_Missense_LoF_SpliceAI_emVar_Enrichment, file="../submission_files/S_Func_Missense_LoF_SpliceAI_emVar_Enrichment.txt", sep="\t")
