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
corehaps_variants <- as_tibble(fread("../../results/3-merge_all_variant_annotations/corehaps_variants_all_annotations_summ.txt.gz"))
corehaps_variants_clinvar_print <- corehaps_variants %>% 
	dplyr::select(VariantID, ClinVar_20230326_summ_Summary) %>% 
	filter(!is.na(ClinVar_20230326_summ_Summary))
names(corehaps_variants_clinvar_print) <- c("Variant ID", "ClinVar 20230326")
write_tsv(corehaps_variants_clinvar_print, "../../results/4-table_corehaps_finemapped_UKBB_BBJ_overlap/Supplementary_Table_SX-Corehaps_annotated_in_ClinVar_in_final_dataset.txt")
