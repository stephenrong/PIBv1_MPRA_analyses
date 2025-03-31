#!/bin/R

# load packages
library(tidyverse)
library(data.table)
library(eulerr)
library(ggpubr)
library(scales)

# load bio packages
library(plyranges)
library(rtracklayer)

# load files
Denisovan_corehaps_variants <- as_tibble(fread("../../results/3-merge_all_variant_annotations/corehaps_variants_all_annotations_summ.txt.gz"))

Denisovan_corehaps_variants_temp <- Denisovan_corehaps_variants %>% 
	dplyr::select(starts_with("CorehapID_")) %>% 
	mutate_all(function(x) {grepl("Denisovan", x)}) %>% 
	rowSums()
Denisovan_corehaps_variants <- Denisovan_corehaps_variants[which(Denisovan_corehaps_variants_temp > 0),]

# Get Ensembl VEP missense variant and gene tables
# Denisovan_corehaps_variants_Ensembl_VEP_missense <- Denisovan_corehaps_variants %>% 
# 	filter(Ensembl_VEP_summ_Summary %in% c("missense_variant")) %>% 
# 	filter(!is.na(Ensembl_VEP_summ_Symbol))

Denisovan_corehaps_variants_Ensembl_VEP_missense <- Denisovan_corehaps_variants %>% 
	filter(Ensembl_VEP_summ_Summary %in% c("missense_variant")) %>% 
	filter(!is.na(Ensembl_VEP_summ_Gene)) %>% 
	mutate(Ensembl_VEP_summ_Symbol = ifelse(is.na(Ensembl_VEP_summ_Symbol), paste0(Ensembl_VEP_summ_Gene, "_failed_convert"), Ensembl_VEP_summ_Symbol))

Denisovan_corehaps_variants_Ensembl_VEP_missense_abbrev_variants <- Denisovan_corehaps_variants_Ensembl_VEP_missense %>% 
	dplyr::select(VariantID, Ensembl_VEP_summ_Summary, Ensembl_VEP_summ_Symbol, CADD_summ_Score, phyloP_mam241_summ_Score) %>% 
	mutate(VariantID = paste0("chr", gsub("_", "-", VariantID))) %>% 
	mutate(CADD_summ_Score = signif(CADD_summ_Score, 3)) %>% 
	mutate(phyloP_mam241_summ_Score = signif(phyloP_mam241_summ_Score, 3)) %>% 
	mutate(Ensembl_VEP_summ_Summary = gsub("_", " ", Ensembl_VEP_summ_Summary))

Denisovan_corehaps_variants_Ensembl_VEP_missense_abbrev_genes <- Denisovan_corehaps_variants_Ensembl_VEP_missense_abbrev_variants %>% 
	group_by(Ensembl_VEP_summ_Symbol) %>% summarise(`Missense count` = dplyr::n(), `Missense variant list` = paste(VariantID, collapse=","))

Denisovan_corehaps_variants_Ensembl_VEP_missense_abbrev_genes_collapse <- Denisovan_corehaps_variants_Ensembl_VEP_missense_abbrev_genes %>% 
	group_by(`Missense count`) %>% summarise(`Count` = dplyr::n(), `Genes` = paste(Ensembl_VEP_summ_Symbol, collapse=",")) %>% ungroup() %>% 
	arrange(desc(`Missense count`))

names(Denisovan_corehaps_variants_Ensembl_VEP_missense_abbrev_variants) <- c("Variant ID", "Ensembl VEP", "Gene symbol", "CADD score", "phyloP score")
write_tsv(Denisovan_corehaps_variants_Ensembl_VEP_missense_abbrev_variants, "../../results/4-table_consequence_gene_lists/Supplementary_Table_SX-Denisovan_Corehaps_introgressed_missense_variants_in_final_dataset.txt")
names(Denisovan_corehaps_variants_Ensembl_VEP_missense_abbrev_genes) <- c("Gene symbol", "Missense count", "Variant ID list")
write_tsv(Denisovan_corehaps_variants_Ensembl_VEP_missense_abbrev_genes, "../../results/4-table_consequence_gene_lists/Supplementary_Table_SX-Genes_with_Denisovan_corehaps_introgressed_missense_variants_in_final_dataset.txt")
names(Denisovan_corehaps_variants_Ensembl_VEP_missense_abbrev_genes_collapse) <- c("Missense count", "Gene count", "Gene symbol list")
write_tsv(Denisovan_corehaps_variants_Ensembl_VEP_missense_abbrev_genes_collapse, "../../results/4-table_consequence_gene_lists/Supplementary_Table_SX-Gene_lists_by_count_of_Denisovan_corehaps_introgressed_missense_variants_in_final_dataset.txt")


# Get Ensembl VEP LoF variant and gene tables
Denisovan_corehaps_variants_Ensembl_VEP_LoF <- Denisovan_corehaps_variants %>% 
	filter(Ensembl_VEP_summ_Summary %in% c("splice_acceptor_variant", "splice_donor_variant", "start_lost", "stop_gained", "stop_lost")) %>% 
	filter(!is.na(Ensembl_VEP_summ_Symbol))

Denisovan_corehaps_variants_Ensembl_VEP_LoF_abbrev_variants <- Denisovan_corehaps_variants_Ensembl_VEP_LoF %>% 
	dplyr::select(VariantID, Ensembl_VEP_summ_Summary, Ensembl_VEP_summ_Symbol, CADD_summ_Score, phyloP_mam241_summ_Score) %>% 
	mutate(VariantID = paste0("chr", gsub("_", "-", VariantID))) %>% 
	mutate(CADD_summ_Score = signif(CADD_summ_Score, 3)) %>% 
	mutate(phyloP_mam241_summ_Score = signif(phyloP_mam241_summ_Score, 3)) %>% 
	mutate(Ensembl_VEP_summ_Summary = gsub("_", " ", Ensembl_VEP_summ_Summary))

Denisovan_corehaps_variants_Ensembl_VEP_LoF_abbrev_genes <- Denisovan_corehaps_variants_Ensembl_VEP_LoF_abbrev_variants %>% 
	group_by(Ensembl_VEP_summ_Symbol, Ensembl_VEP_summ_Summary) %>% summarise(`LoF count` = dplyr::n(), `LoF variant list` = paste(VariantID, collapse=","))

Denisovan_corehaps_variants_Ensembl_VEP_LoF_abbrev_genes_collapse <- Denisovan_corehaps_variants_Ensembl_VEP_LoF_abbrev_genes %>% 
	group_by(Ensembl_VEP_summ_Summary, `LoF count`) %>% summarise(`Count` = dplyr::n(), `Genes` = paste(Ensembl_VEP_summ_Symbol, collapse=",")) %>% ungroup() %>% 
	arrange(desc(`LoF count`))

names(Denisovan_corehaps_variants_Ensembl_VEP_LoF_abbrev_variants) <- c("Variant ID", "Ensembl VEP", "Gene symbol", "CADD score", "phyloP score")
write_tsv(Denisovan_corehaps_variants_Ensembl_VEP_LoF_abbrev_variants, "../../results/4-table_consequence_gene_lists/Supplementary_Table_SX-Denisovan_Corehaps_introgressed_LoF_variants_in_final_dataset.txt")
names(Denisovan_corehaps_variants_Ensembl_VEP_LoF_abbrev_genes) <- c("Gene symbol", "LoF type", "LoF count", "Variant ID list")
write_tsv(Denisovan_corehaps_variants_Ensembl_VEP_LoF_abbrev_genes, "../../results/4-table_consequence_gene_lists/Supplementary_Table_SX-Genes_with_Denisovan_corehaps_introgressed_LoF_variants_in_final_dataset.txt")
names(Denisovan_corehaps_variants_Ensembl_VEP_LoF_abbrev_genes_collapse) <- c("LoF type", "LoF count", "Gene count", "Gene symbol list")
write_tsv(Denisovan_corehaps_variants_Ensembl_VEP_LoF_abbrev_genes_collapse, "../../results/4-table_consequence_gene_lists/Supplementary_Table_SX-Gene_lists_by_count_of_Denisovan_corehaps_introgressed_LoF_variants_in_final_dataset.txt")


# Get SpliceAI variant and gene tables
# Denisovan_corehaps_variants_SpliceAI <- Denisovan_corehaps_variants %>% 
# 	filter(!is.na(SpliceAI_raw_summ_Summary)) %>% 
# 	filter(!is.na(SpliceAI_raw_summ_Symbol))

Denisovan_corehaps_variants_SpliceAI <- Denisovan_corehaps_variants %>% 
	filter(!is.na(SpliceAI_raw_summ_Summary)) %>% 
	filter(!is.na(SpliceAI_raw_summ_Gene)) %>% 
	mutate(SpliceAI_raw_summ_Symbol = ifelse(is.na(SpliceAI_raw_summ_Symbol), paste0(SpliceAI_raw_summ_Gene, "_failed_convert"), SpliceAI_raw_summ_Symbol))

Denisovan_corehaps_variants_SpliceAI_abbrev_variants <- Denisovan_corehaps_variants_SpliceAI %>% 
	dplyr::select(VariantID, SpliceAI_raw_summ_Summary, SpliceAI_raw_summ_Score, SpliceAI_raw_summ_Description, SpliceAI_raw_summ_Symbol, CADD_summ_Score, phyloP_mam241_summ_Score) %>% 
	mutate(VariantID = paste0("chr", gsub("_", "-", VariantID))) %>% 
	mutate(CADD_summ_Score = signif(CADD_summ_Score, 3)) %>% 
	mutate(phyloP_mam241_summ_Score = signif(phyloP_mam241_summ_Score, 3)) %>% 
	mutate(SpliceAI_raw_summ_Score = signif(SpliceAI_raw_summ_Score, 2))

Denisovan_corehaps_variants_SpliceAI_abbrev_variants_temp <- Denisovan_corehaps_variants_SpliceAI_abbrev_variants %>% 
	group_by(SpliceAI_raw_summ_Symbol) %>% arrange(desc(SpliceAI_raw_summ_Score)) %>% 
	mutate(SpliceAI_raw_summ_Symbol = paste0(SpliceAI_raw_summ_Symbol, paste0(rep("*", dplyr::n()-1), collapse=""))) %>% 
	slice(1) %>% ungroup()

Denisovan_corehaps_variants_SpliceAI_abbrev_genes <- Denisovan_corehaps_variants_SpliceAI_abbrev_variants_temp %>% 
	group_by(SpliceAI_raw_summ_Symbol, SpliceAI_raw_summ_Summary) %>% summarise(`SpliceAI variant list` = paste(VariantID, collapse=","))

Denisovan_corehaps_variants_SpliceAI_abbrev_genes_collapse <- Denisovan_corehaps_variants_SpliceAI_abbrev_genes %>% 
	group_by(SpliceAI_raw_summ_Summary) %>% summarise(`Count` = dplyr::n(), `Genes` = paste(SpliceAI_raw_summ_Symbol, collapse=",")) %>% ungroup()

names(Denisovan_corehaps_variants_SpliceAI_abbrev_variants) <- c("Variant ID", "SpliceAI category", "SpliceAI max score", "SpliceAI prediction", "Gene symbol", "CADD score", "phyloP score")
write_tsv(Denisovan_corehaps_variants_SpliceAI_abbrev_variants, "../../results/4-table_consequence_gene_lists/Supplementary_Table_SX-Denisovan_Corehaps_introgressed_SpliceAI_variants_in_final_dataset.txt")
names(Denisovan_corehaps_variants_SpliceAI_abbrev_genes) <- c("Gene symbol", "SpliceAI category", "Variant ID list")
write_tsv(Denisovan_corehaps_variants_SpliceAI_abbrev_genes, "../../results/4-table_consequence_gene_lists/Supplementary_Table_SX-Genes_with_Denisovan_corehaps_introgressed_SpliceAI_variants_in_final_dataset.txt")
names(Denisovan_corehaps_variants_SpliceAI_abbrev_genes_collapse) <- c("SpliceAI category count", "Gene count", "Gene symbol list")
write_tsv(Denisovan_corehaps_variants_SpliceAI_abbrev_genes_collapse, "../../results/4-table_consequence_gene_lists/Supplementary_Table_SX-Gene_lists_by_count_of_Denisovan_corehaps_introgressed_SpliceAI_variants_in_final_dataset.txt")
