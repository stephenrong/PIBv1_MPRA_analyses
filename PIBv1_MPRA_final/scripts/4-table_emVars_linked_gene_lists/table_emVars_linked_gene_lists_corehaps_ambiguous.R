#!/bin/R

# load packages
library(tidyverse)
library(data.table)

# load bio packages
library(plyranges)

# load files
Ambiguous_corehaps_variants <- as_tibble(fread("../../results/3-merge_all_variant_annotations/corehaps_variants_all_annotations_summ.txt.gz"))

Ambiguous_corehaps_variants_temp <- Ambiguous_corehaps_variants %>% 
	dplyr::select(starts_with("CorehapID_")) %>% 
	mutate_all(function(x) {grepl("Ambiguous", x)}) %>% 
	rowSums()
Ambiguous_corehaps_variants <- Ambiguous_corehaps_variants[which(Ambiguous_corehaps_variants_temp > 0),]

# Look up each variant and their associated and genes split by method
Ambiguous_corehaps_variants_emVars <- Ambiguous_corehaps_variants %>% 
	filter(!is.na(MPRA_K562_summ_Summary) | !is.na(MPRA_Jurkat_summ_Summary))

Ambiguous_corehaps_variants_emVars_linked <- Ambiguous_corehaps_variants_emVars %>% 
	filter(!is.na(ABC_summ_Summary) | !is.na(EpiMap_summ_Summary) | !is.na(Roadmap_Epigenomics_summ_Summary)) %>% 
	mutate(CADD_summ_Score = signif(CADD_summ_Score, 3)) %>% 
	mutate(phyloP_mam241_summ_Score = signif(phyloP_mam241_summ_Score, 3))

Ambiguous_corehaps_variants_emVars_linked_abbrev_variants <- Ambiguous_corehaps_variants_emVars_linked %>% 
	dplyr::select(VariantID, MPRA_K562_summ_Summary, MPRA_K562_summ_Description_skew_log10FDR, MPRA_K562_summ_Description_skew_log2FC, MPRA_Jurkat_summ_Summary, MPRA_Jurkat_summ_Description_skew_log10FDR, MPRA_Jurkat_summ_Description_skew_log2FC, ABC_summ_Score, ABC_summ_Symbol, ABC_summ_Biosample, ABC_summ_Category, EpiMap_summ_Score, EpiMap_summ_Symbol, EpiMap_summ_Group, Roadmap_Epigenomics_summ_Score, Roadmap_Epigenomics_summ_Symbol, Roadmap_Epigenomics_summ_Biosample, Roadmap_Epigenomics_summ_Category, CADD_summ_Score, phyloP_mam241_summ_Score) %>% 
	mutate(VariantID = paste0("chr", gsub("_", "-", VariantID))) %>% 
	mutate(MPRA_K562_summ_Summary = gsub("_", " ", MPRA_K562_summ_Summary)) %>% 
	mutate(MPRA_K562_summ_Description_skew_log10FDR = ifelse(!is.na(MPRA_K562_summ_Summary), signif(MPRA_K562_summ_Description_skew_log10FDR, 3), NA)) %>% 
	mutate(MPRA_K562_summ_Description_skew_log2FC = ifelse(!is.na(MPRA_K562_summ_Summary), signif(MPRA_K562_summ_Description_skew_log2FC, 3), NA)) %>% 
	mutate(MPRA_Jurkat_summ_Summary = gsub("_", " ", MPRA_Jurkat_summ_Summary)) %>% 
	mutate(MPRA_Jurkat_summ_Description_skew_log10FDR = ifelse(!is.na(MPRA_Jurkat_summ_Summary), signif(MPRA_Jurkat_summ_Description_skew_log10FDR, 3), NA)) %>% 
	mutate(MPRA_Jurkat_summ_Description_skew_log2FC = ifelse(!is.na(MPRA_Jurkat_summ_Summary), signif(MPRA_Jurkat_summ_Description_skew_log2FC, 3), NA))

Ambiguous_corehaps_variants_emVars_linked_abbrev_variants_named <- Ambiguous_corehaps_variants_emVars_linked_abbrev_variants
names(Ambiguous_corehaps_variants_emVars_linked_abbrev_variants_named) <- c("Variant ID", "MPRA K562", "MPRA K562 log10FDR", "MPRA K562 log2FC", "MPRA Jurkat", "MPRA Jurkat log10FDR", "MPRA Jurkat log2FC", "ABC score list", "ABC gene symbol list", "ABC biosample list", "ABC category list", "EpiMap score list", "EpiMap gene symbol list", "EpiMap group list", "Roadmap Epigenomics score list", "Roadmap Epigenomics gene symbol list", "Roadmap Epigenomics biosample list", "Roadmap Epigenomics category list", "CADD score", "phyloP score")
write_tsv(Ambiguous_corehaps_variants_emVars_linked_abbrev_variants_named, "../../results/4-table_emVars_linked_gene_lists/Supplementary_Table_SX-Ambiguous_Corehaps_introgressed_MPRA_emVar_variants_with_gene_enhancer_links_in_final_dataset.txt")

unique_entries_string <- function(string) {
	paste(unique(strsplit(string, ",")[[1]]), collapse=",")
}

# Look up each variant and their associated methods by method
Ambiguous_corehaps_variants_emVars_linked_variant_gene_pairs <- Ambiguous_corehaps_variants_emVars_linked %>% 
	dplyr::select(VariantID, ABC_summ_Symbol, EpiMap_summ_Symbol, Roadmap_Epigenomics_summ_Symbol, CADD_summ_Score, phyloP_mam241_summ_Score) %>% 
	# split genes per method into separate rows
	separate_longer_delim(ABC_summ_Symbol, delim = ",") %>% 
	separate_longer_delim(EpiMap_summ_Symbol, delim = ",") %>% 
	separate_longer_delim(Roadmap_Epigenomics_summ_Symbol, delim = ",") %>% 
	# pivot and count gene-enhancer linkage methods
	pivot_longer(cols=c(ABC_summ_Symbol, EpiMap_summ_Symbol, Roadmap_Epigenomics_summ_Symbol), names_to="Method", values_to="Gene") %>% 
	filter(Gene != "NA") %>% 
	# remove duplicate variant to gene pairs
	unique() %>% 
	# group by genes and collapse methods
	group_by(VariantID, Gene, CADD_summ_Score, phyloP_mam241_summ_Score) %>% summarise(`Method count` = dplyr::n(), `Method list` = paste(Method, collapse=",")) %>% ungroup() %>% 
	# clean up columns
	mutate(VariantID = paste0("chr", gsub("_", "-", VariantID))) %>% 
	mutate(`Method list` = gsub("ABC_summ_Symbol", "ABC", `Method list`)) %>% 
	mutate(`Method list` = gsub("EpiMap_summ_Symbol", "EpiMap", `Method list`)) %>% 
	mutate(`Method list` = gsub("Roadmap_Epigenomics_summ_Symbol", "Roadmap Epigenomics", `Method list`)) %>% 
	dplyr::select(VariantID, Gene, `Method count`, `Method list`, CADD_summ_Score, phyloP_mam241_summ_Score)

Ambiguous_corehaps_variants_emVars_linked_variant_gene_pairs_named <- Ambiguous_corehaps_variants_emVars_linked_variant_gene_pairs
names(Ambiguous_corehaps_variants_emVars_linked_variant_gene_pairs_named) <- c("Variant ID", "Gene symbol", "Gene-enhancer methods count", "Gene-enhancer methods list", "CADD score", "phyloP score")
write_tsv(Ambiguous_corehaps_variants_emVars_linked_variant_gene_pairs_named, "../../results/4-table_emVars_linked_gene_lists/Supplementary_Table_SX-Pairs_of_Ambiguous_corehaps_introgressed_MPRA_emVar_variants_and_linked_genes_in_final_dataset.txt")

# Look up for a given gene which variants are linked
Ambiguous_corehaps_variants_emVars_linked_genes <- Ambiguous_corehaps_variants_emVars_linked_variant_gene_pairs %>% 
	group_by(`Gene`) %>% summarise(`Linked variant count` = dplyr::n(), `Variant ID list` = paste(`VariantID`, collapse=",")) %>% ungroup() %>% 
	arrange(desc(`Linked variant count`))
Ambiguous_corehaps_variants_emVars_linked_genes_named <- Ambiguous_corehaps_variants_emVars_linked_genes
names(Ambiguous_corehaps_variants_emVars_linked_genes_named) <- c("Gene symbol", "emVar+linked variant count", "Variant ID list")
write_tsv(Ambiguous_corehaps_variants_emVars_linked_genes_named, "../../results/4-table_emVars_linked_gene_lists/Supplementary_Table_SX-Genes_with_Ambiguous_corehaps_introgressed_MPRA_emVars_in_final_dataset.txt")

# Look up which genes have the most links
Ambiguous_corehaps_variants_emVars_linked_genes_collapse <- Ambiguous_corehaps_variants_emVars_linked_genes %>% 
	group_by(`Linked variant count`) %>% summarise(`Count` = dplyr::n(), `Genes` = paste(`Gene`, collapse=",")) %>% ungroup() %>% 
	arrange(desc(`Linked variant count`))
Ambiguous_corehaps_variants_emVars_linked_genes_collapse_named <- Ambiguous_corehaps_variants_emVars_linked_genes_collapse
names(Ambiguous_corehaps_variants_emVars_linked_genes_collapse_named) <- c("emVar+linked variant count", "Gene count", "Gene symbol list")
write_tsv(Ambiguous_corehaps_variants_emVars_linked_genes_collapse_named, "../../results/4-table_emVars_linked_gene_lists/Supplementary_Table_SX-Gene_lists_by_count_of_Ambiguous_corehaps_introgressed_MPRA_emVars_in_final_dataset.txt")
	