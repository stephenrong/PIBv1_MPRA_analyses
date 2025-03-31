#!/bin/R

# for data analysis
library(tidyverse)
library(data.table)

# for bio data analysis
library(plyranges)
library(rtracklayer)
source("../shared_functions/seqinfo_fix_change.R")

# load variants
adaptive_variants_tb <- as_tibble(readRDS("../../results/1a-preprocess_PIBv1_MPRA_pilot/adaptive_variants.rds"))[,1:10]

# load DHS_index_vocabulary
DHS_index_vocabulary_tb <- as_tibble(fread(paste0("../../../../Datasets/gene_regulation_element_catalogs/DHS_index_vocabulary/data_cleanup/GRCh37_merged/DHS_index_vocabulary.txt.gz")))

# rename component
DHS_index_vocabulary_tb <- DHS_index_vocabulary_tb %>% 
	mutate(component = gsub(" ", "_", component))

# DHS_index_vocabulary names
names(DHS_index_vocabulary_tb)[6:length(DHS_index_vocabulary_tb)] <- paste("DHS_index_vocabulary_orig", names(DHS_index_vocabulary_tb)[6:length(DHS_index_vocabulary_tb)], sep="_")
adaptive_variants_DHS_index_vocabulary_names <- names(DHS_index_vocabulary_tb)[6:length(DHS_index_vocabulary_tb)]

# overlap variants and DHS_index_vocabulary
adaptive_variants_DHS_index_vocabulary <- as_tibble(find_overlaps(GRanges(adaptive_variants_tb), GRanges(DHS_index_vocabulary_tb)))

# DHS_index_vocabulary names
adaptive_variants_DHS_index_vocabulary_names <- names(adaptive_variants_DHS_index_vocabulary)[11:length(adaptive_variants_DHS_index_vocabulary)]

# collapse by component
adaptive_variants_DHS_index_vocabulary_collapse <- adaptive_variants_DHS_index_vocabulary %>% 
	group_by(seqnames, start, end, width, strand, VariantID, VariantCHROM, VariantPOS, VariantREF, VariantALT) %>% 
	dplyr::summarise_at(adaptive_variants_DHS_index_vocabulary_names, function(x) {paste(x, collapse=",")}) %>% 
	ungroup() %>% 
	mutate(DHS_index_vocabulary_summ_Component = DHS_index_vocabulary_orig_component) %>% 
	mutate(DHS_index_vocabulary_summ_Summary = "DHS_index_vocabulary_overlap") %>% 
	dplyr::select(seqnames, start, end, width, strand, VariantID, VariantCHROM, VariantPOS, VariantREF, VariantALT, 
		starts_with("DHS_index_vocabulary_orig"), DHS_index_vocabulary_summ_Summary, DHS_index_vocabulary_summ_Component)

# pivot by component
adaptive_variants_DHS_index_vocabulary_pivot <- adaptive_variants_DHS_index_vocabulary %>% 
	mutate(DHS_index_vocabulary_orig_component_temp = DHS_index_vocabulary_orig_component) %>% 
	pivot_wider(names_from=DHS_index_vocabulary_orig_component, names_sep="-", names_prefix="DHS_index_vocabulary_summ_Component-", values_from=DHS_index_vocabulary_orig_component_temp)

# join together
adaptive_variants_DHS_index_vocabulary <- full_join(adaptive_variants_DHS_index_vocabulary_collapse, adaptive_variants_DHS_index_vocabulary_pivot)

# merge to full variants
adaptive_variants_tb <- as_tibble(fread("../../results/1a-preprocess_PIBv1_MPRA_pilot/adaptive_variants.txt.gz")) %>% 
	dplyr::select(seqnames, start, end, width, strand, VariantID, VariantCHROM, VariantPOS, VariantREF, VariantALT)
names_temp <- intersect(names(adaptive_variants_tb), names(adaptive_variants_DHS_index_vocabulary))
adaptive_variants_DHS_index_vocabulary <- adaptive_variants_tb %>% mutate_at(names_temp, as.character)%>% 
	left_join(adaptive_variants_DHS_index_vocabulary %>% mutate_at(names_temp, as.character))

# save all variants
write_tsv(adaptive_variants_DHS_index_vocabulary, gzfile("../../results/2r-annotate_variants_DHS_index_vocabulary/adaptive_variants_DHS_index_vocabulary.txt.gz"))

# save summ variants
adaptive_variants_DHS_index_vocabulary_summ <- adaptive_variants_DHS_index_vocabulary %>% 
	dplyr::select(seqnames, start, end, width, strand, starts_with("Variant"), contains("_summ_"))
write_tsv(adaptive_variants_DHS_index_vocabulary_summ, gzfile("../../results/2r-annotate_variants_DHS_index_vocabulary/adaptive_variants_DHS_index_vocabulary_summ.txt.gz"))
sort(table(adaptive_variants_DHS_index_vocabulary_summ$DHS_index_vocabulary_summ_Summary))
