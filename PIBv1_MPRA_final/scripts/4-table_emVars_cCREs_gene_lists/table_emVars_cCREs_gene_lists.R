#!/bin/R

# load packages
library(tidyverse)
library(data.table)

# load bio packages
library(plyranges)

# load files
mpra_emVars_variants_ENCODE_cCREs_UKBB_finemap_pip0.1 <- as_tibble(fread("../../results/3-overlap_emVars_GWAS_cCREs/mpra_emVars_variants_ENCODE_cCREs_UKBB_finemap_pip0.1.tsv.gz"))
mpra_emVars_variants_ENCODE_cCREs_BBJ_finemap_pip0.1 <- as_tibble(fread("../../results/3-overlap_emVars_GWAS_cCREs/mpra_emVars_variants_ENCODE_cCREs_BBJ_finemap_pip0.1.tsv.gz"))

mpra_emVars_variants_UKBB_finemap_pip0.1 <- as_tibble(fread("../../results/3-overlap_emVars_GWAS_cCREs/mpra_emVars_variants_UKBB_finemap_pip0.1.tsv.gz"))
mpra_emVars_variants_BBJ_finemap_pip0.1 <- as_tibble(fread("../../results/3-overlap_emVars_GWAS_cCREs/mpra_emVars_variants_BBJ_finemap_pip0.1.tsv.gz"))

# clean up
mpra_emVars_variants_ENCODE_cCREs_UKBB_finemap_pip0.1_temp <- mpra_emVars_variants_ENCODE_cCREs_UKBB_finemap_pip0.1 %>% 
	dplyr::select(MPRA_VariantID, MPRA_K562_summ_Summary, MPRA_K562_summ_Description_skew_log10FDR, MPRA_K562_summ_Description_skew_log2FC, MPRA_Jurkat_summ_Summary, MPRA_Jurkat_summ_Description_skew_log10FDR, MPRA_Jurkat_summ_Description_skew_log2FC, ENCODE_cCREs_orig_class, ENCODE_cCREs_orig_element, UKBB_finemap_VariantID, UKBB_finemap_orig_method, UKBB_finemap_orig_trait, UKBB_finemap_orig_pip, ABC_summ_Score, ABC_summ_Symbol, ABC_summ_Biosample, ABC_summ_Category, EpiMap_summ_Score, EpiMap_summ_Symbol, EpiMap_summ_Group, Roadmap_Epigenomics_summ_Score, Roadmap_Epigenomics_summ_Symbol, Roadmap_Epigenomics_summ_Biosample, Roadmap_Epigenomics_summ_Category, nearest_gene_summ_Symbol) %>% 
	mutate(MPRA_VariantID = paste0("chr", gsub("_", "-", MPRA_VariantID))) %>% 
	mutate(MPRA_K562_summ_Summary = gsub("_", " ", MPRA_K562_summ_Summary)) %>% 
	mutate(MPRA_K562_summ_Description_skew_log10FDR = ifelse(!is.na(MPRA_K562_summ_Summary), signif(MPRA_K562_summ_Description_skew_log10FDR, 3), NA)) %>% 
	mutate(MPRA_K562_summ_Description_skew_log2FC = ifelse(!is.na(MPRA_K562_summ_Summary), signif(MPRA_K562_summ_Description_skew_log2FC, 3), NA)) %>% 
	mutate(MPRA_Jurkat_summ_Summary = gsub("_", " ", MPRA_Jurkat_summ_Summary)) %>% 
	mutate(MPRA_Jurkat_summ_Description_skew_log10FDR = ifelse(!is.na(MPRA_Jurkat_summ_Summary), signif(MPRA_Jurkat_summ_Description_skew_log10FDR, 3), NA)) %>% 
	mutate(MPRA_Jurkat_summ_Description_skew_log2FC =ifelse(!is.na(MPRA_Jurkat_summ_Summary),  signif(MPRA_Jurkat_summ_Description_skew_log2FC, 3), NA)) %>% 
	mutate(UKBB_finemap_VariantID = paste0("chr", gsub("_", "-", UKBB_finemap_VariantID)))
mpra_emVars_variants_ENCODE_cCREs_BBJ_finemap_pip0.1_temp <- mpra_emVars_variants_ENCODE_cCREs_BBJ_finemap_pip0.1 %>% 
	dplyr::select(MPRA_VariantID, MPRA_K562_summ_Summary, MPRA_K562_summ_Description_skew_log10FDR, MPRA_K562_summ_Description_skew_log2FC, MPRA_Jurkat_summ_Summary, MPRA_Jurkat_summ_Description_skew_log10FDR, MPRA_Jurkat_summ_Description_skew_log2FC, ENCODE_cCREs_orig_class, ENCODE_cCREs_orig_element, BBJ_finemap_VariantID, BBJ_finemap_orig_method, BBJ_finemap_orig_trait, BBJ_finemap_orig_pip, ABC_summ_Score, ABC_summ_Symbol, ABC_summ_Biosample, ABC_summ_Category, EpiMap_summ_Score, EpiMap_summ_Symbol, EpiMap_summ_Group, Roadmap_Epigenomics_summ_Score, Roadmap_Epigenomics_summ_Symbol, Roadmap_Epigenomics_summ_Biosample, Roadmap_Epigenomics_summ_Category, nearest_gene_summ_Symbol) %>% 
	mutate(MPRA_VariantID = paste0("chr", gsub("_", "-", MPRA_VariantID))) %>% 
	mutate(MPRA_K562_summ_Summary = gsub("_", " ", MPRA_K562_summ_Summary)) %>% 
	mutate(MPRA_K562_summ_Description_skew_log10FDR = ifelse(!is.na(MPRA_K562_summ_Summary), signif(MPRA_K562_summ_Description_skew_log10FDR, 3), NA)) %>% 
	mutate(MPRA_K562_summ_Description_skew_log2FC = ifelse(!is.na(MPRA_K562_summ_Summary), signif(MPRA_K562_summ_Description_skew_log2FC, 3), NA)) %>% 
	mutate(MPRA_Jurkat_summ_Summary = gsub("_", " ", MPRA_Jurkat_summ_Summary)) %>% 
	mutate(MPRA_Jurkat_summ_Description_skew_log10FDR = ifelse(!is.na(MPRA_Jurkat_summ_Summary), signif(MPRA_Jurkat_summ_Description_skew_log10FDR, 3), NA)) %>% 
	mutate(MPRA_Jurkat_summ_Description_skew_log2FC =ifelse(!is.na(MPRA_Jurkat_summ_Summary),  signif(MPRA_Jurkat_summ_Description_skew_log2FC, 3), NA)) %>% 
	mutate(BBJ_finemap_VariantID = paste0("chr", gsub("_", "-", BBJ_finemap_VariantID)))

mpra_emVars_variants_UKBB_finemap_pip0.1_temp <- mpra_emVars_variants_UKBB_finemap_pip0.1 %>% 
	dplyr::select(MPRA_VariantID, MPRA_K562_summ_Summary, MPRA_K562_summ_Description_skew_log10FDR, MPRA_K562_summ_Description_skew_log2FC, MPRA_Jurkat_summ_Summary, MPRA_Jurkat_summ_Description_skew_log10FDR, MPRA_Jurkat_summ_Description_skew_log2FC, UKBB_finemap_VariantID, UKBB_finemap_orig_method, UKBB_finemap_orig_trait, UKBB_finemap_orig_pip, ABC_summ_Score, ABC_summ_Symbol, ABC_summ_Biosample, ABC_summ_Category, EpiMap_summ_Score, EpiMap_summ_Symbol, EpiMap_summ_Group, Roadmap_Epigenomics_summ_Score, Roadmap_Epigenomics_summ_Symbol, Roadmap_Epigenomics_summ_Biosample, Roadmap_Epigenomics_summ_Category, nearest_gene_summ_Symbol) %>% 
	mutate(MPRA_VariantID = paste0("chr", gsub("_", "-", MPRA_VariantID))) %>% 
	mutate(MPRA_K562_summ_Summary = gsub("_", " ", MPRA_K562_summ_Summary)) %>% 
	mutate(MPRA_K562_summ_Description_skew_log10FDR = ifelse(!is.na(MPRA_K562_summ_Summary), signif(MPRA_K562_summ_Description_skew_log10FDR, 3), NA)) %>% 
	mutate(MPRA_K562_summ_Description_skew_log2FC = ifelse(!is.na(MPRA_K562_summ_Summary), signif(MPRA_K562_summ_Description_skew_log2FC, 3), NA)) %>% 
	mutate(MPRA_Jurkat_summ_Summary = gsub("_", " ", MPRA_Jurkat_summ_Summary)) %>% 
	mutate(MPRA_Jurkat_summ_Description_skew_log10FDR = ifelse(!is.na(MPRA_Jurkat_summ_Summary), signif(MPRA_Jurkat_summ_Description_skew_log10FDR, 3), NA)) %>% 
	mutate(MPRA_Jurkat_summ_Description_skew_log2FC =ifelse(!is.na(MPRA_Jurkat_summ_Summary),  signif(MPRA_Jurkat_summ_Description_skew_log2FC, 3), NA)) %>% 
	mutate(UKBB_finemap_VariantID = paste0("chr", gsub("_", "-", UKBB_finemap_VariantID)))
mpra_emVars_variants_BBJ_finemap_pip0.1_temp <- mpra_emVars_variants_BBJ_finemap_pip0.1 %>% 
	dplyr::select(MPRA_VariantID, MPRA_K562_summ_Summary, MPRA_K562_summ_Description_skew_log10FDR, MPRA_K562_summ_Description_skew_log2FC, MPRA_Jurkat_summ_Summary, MPRA_Jurkat_summ_Description_skew_log10FDR, MPRA_Jurkat_summ_Description_skew_log2FC, BBJ_finemap_VariantID, BBJ_finemap_orig_method, BBJ_finemap_orig_trait, BBJ_finemap_orig_pip, ABC_summ_Score, ABC_summ_Symbol, ABC_summ_Biosample, ABC_summ_Category, EpiMap_summ_Score, EpiMap_summ_Symbol, EpiMap_summ_Group, Roadmap_Epigenomics_summ_Score, Roadmap_Epigenomics_summ_Symbol, Roadmap_Epigenomics_summ_Biosample, Roadmap_Epigenomics_summ_Category, nearest_gene_summ_Symbol) %>% 
	mutate(MPRA_VariantID = paste0("chr", gsub("_", "-", MPRA_VariantID))) %>% 
	mutate(MPRA_K562_summ_Summary = gsub("_", " ", MPRA_K562_summ_Summary)) %>% 
	mutate(MPRA_K562_summ_Description_skew_log10FDR = ifelse(!is.na(MPRA_K562_summ_Summary), signif(MPRA_K562_summ_Description_skew_log10FDR, 3), NA)) %>% 
	mutate(MPRA_K562_summ_Description_skew_log2FC = ifelse(!is.na(MPRA_K562_summ_Summary), signif(MPRA_K562_summ_Description_skew_log2FC, 3), NA)) %>% 
	mutate(MPRA_Jurkat_summ_Summary = gsub("_", " ", MPRA_Jurkat_summ_Summary)) %>% 
	mutate(MPRA_Jurkat_summ_Description_skew_log10FDR = ifelse(!is.na(MPRA_Jurkat_summ_Summary), signif(MPRA_Jurkat_summ_Description_skew_log10FDR, 3), NA)) %>% 
	mutate(MPRA_Jurkat_summ_Description_skew_log2FC =ifelse(!is.na(MPRA_Jurkat_summ_Summary),  signif(MPRA_Jurkat_summ_Description_skew_log2FC, 3), NA)) %>% 
	mutate(BBJ_finemap_VariantID = paste0("chr", gsub("_", "-", BBJ_finemap_VariantID)))

# rename cols
mpra_emVars_variants_ENCODE_cCREs_UKBB_finemap_pip0.1_named <- mpra_emVars_variants_ENCODE_cCREs_UKBB_finemap_pip0.1_temp
names(mpra_emVars_variants_ENCODE_cCREs_UKBB_finemap_pip0.1_named) <- c("MPRA Variant ID", "MPRA K562", "MPRA K562 log10FDR", "MPRA K562 log2FC", "MPRA Jurkat", "MPRA Jurkat log10FDR", "MPRA Jurkat log2FC", "ENCODE cCRE category", "ENCODE cCRE coordinates", "UKBB Variant ID", "UKBB method", "UKBB trait", "UKBB PIP", "ABC score list", "ABC gene symbol list", "ABC biosample list", "ABC category list", "EpiMap score list", "EpiMap gene symbol list", "EpiMap group list", "Roadmap Epigenomics score list", "Roadmap Epigenomics gene symbol list", "Roadmap Epigenomics biosample list", "Roadmap Epigenomics category list", "Nearest gene symbol list")
mpra_emVars_variants_ENCODE_cCREs_BBJ_finemap_pip0.1_named <- mpra_emVars_variants_ENCODE_cCREs_BBJ_finemap_pip0.1_temp
names(mpra_emVars_variants_ENCODE_cCREs_BBJ_finemap_pip0.1_named) <- c("MPRA Variant ID", "MPRA K562", "MPRA K562 log10FDR", "MPRA K562 log2FC", "MPRA Jurkat", "MPRA Jurkat log10FDR", "MPRA Jurkat log2FC", "ENCODE cCRE category", "ENCODE cCRE coordinates", "BBJ Variant ID", "BBJ method", "BBJ trait", "BBJ PIP", "ABC score list", "ABC gene symbol list", "ABC biosample list", "ABC category list", "EpiMap score list", "EpiMap gene symbol list", "EpiMap group list", "Roadmap Epigenomics score list", "Roadmap Epigenomics gene symbol list", "Roadmap Epigenomics biosample list", "Roadmap Epigenomics category list", "Nearest gene symbol list")
mpra_emVars_variants_UKBB_finemap_pip0.1_named <- mpra_emVars_variants_UKBB_finemap_pip0.1_temp
names(mpra_emVars_variants_UKBB_finemap_pip0.1_named) <- c("MPRA Variant ID", "MPRA K562", "MPRA K562 log10FDR", "MPRA K562 log2FC", "MPRA Jurkat", "MPRA Jurkat log10FDR", "MPRA Jurkat log2FC", "UKBB Variant ID", "UKBB method", "UKBB trait", "UKBB PIP", "ABC score list", "ABC gene symbol list", "ABC biosample list", "ABC category list", "EpiMap score list", "EpiMap gene symbol list", "EpiMap group list", "Roadmap Epigenomics score list", "Roadmap Epigenomics gene symbol list", "Roadmap Epigenomics biosample list", "Roadmap Epigenomics category list", "Nearest gene symbol list")
mpra_emVars_variants_BBJ_finemap_pip0.1_named <- mpra_emVars_variants_BBJ_finemap_pip0.1_temp
names(mpra_emVars_variants_BBJ_finemap_pip0.1_named) <- c("MPRA Variant ID", "MPRA K562", "MPRA K562 log10FDR", "MPRA K562 log2FC", "MPRA Jurkat", "MPRA Jurkat log10FDR", "MPRA Jurkat log2FC", "BBJ Variant ID", "BBJ method", "BBJ trait", "BBJ PIP", "ABC score list", "ABC gene symbol list", "ABC biosample list", "ABC category list", "EpiMap score list", "EpiMap gene symbol list", "EpiMap group list", "Roadmap Epigenomics score list", "Roadmap Epigenomics gene symbol list", "Roadmap Epigenomics biosample list", "Roadmap Epigenomics category list", "Nearest gene symbol list")

# save supp tables
write_tsv(mpra_emVars_variants_ENCODE_cCREs_UKBB_finemap_pip0.1_named, "../../results/4-table_emVars_cCREs_gene_lists/Supplementary_Table_SX-Triplets_of_adaptive_introgressed_MPRA_emVar_variants_overlapping_ENCODE_cCRE_with_UKBB_fine-mapped_variant_PIP_0.1_in_pilot_dataset.txt")
write_tsv(mpra_emVars_variants_ENCODE_cCREs_BBJ_finemap_pip0.1_named, "../../results/4-table_emVars_cCREs_gene_lists/Supplementary_Table_SX-Triplets_of_adaptive_introgressed_MPRA_emVar_variants_overlapping_ENCODE_cCRE_with_BBJ_fine-mapped_variant_PIP_0.1_in_pilot_dataset.txt")
write_tsv(mpra_emVars_variants_UKBB_finemap_pip0.1_named, "../../results/4-table_emVars_cCREs_gene_lists/Supplementary_Table_SX-Pairs_of_adaptive_introgressed_MPRA_emVar_variants_overlapping_UKBB_fine-mapped_variant_PIP_0.1_in_pilot_dataset.txt")
write_tsv(mpra_emVars_variants_BBJ_finemap_pip0.1_named, "../../results/4-table_emVars_cCREs_gene_lists/Supplementary_Table_SX-Pairs_of_adaptive_introgressed_MPRA_emVar_variants_overlapping_BBJ_fine-mapped_variant_PIP_0.1_in_pilot_dataset.txt")
