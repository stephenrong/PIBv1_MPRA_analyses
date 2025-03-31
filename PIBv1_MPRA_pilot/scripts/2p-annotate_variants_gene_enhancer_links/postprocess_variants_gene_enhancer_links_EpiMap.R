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

# load EpiMap
EpiMap_gene_enhancer_links_tb <- as_tibble(fread(paste0("../../../../Datasets/gene_regulation_gene_enhancer_links/EpiMap_gene_enhancer_links/data_cleanup/GRCh37_merged/EpiMap_gene_enhancer_links.txt.gz")))

# remove ambig TSS
EpiMap_gene_enhancer_links_tb <- EpiMap_gene_enhancer_links_tb %>% 
	filter(!is.na(target_symbol)) %>% 
	filter(!is.na(target_tss))

# EpiMap names
names(EpiMap_gene_enhancer_links_tb)[6:length(EpiMap_gene_enhancer_links_tb)] <- paste("EpiMap_orig", names(EpiMap_gene_enhancer_links_tb)[6:length(EpiMap_gene_enhancer_links_tb)], sep="_")
adaptive_variants_EpiMap_names <- names(EpiMap_gene_enhancer_links_tb)[6:length(EpiMap_gene_enhancer_links_tb)]

# convert to GRange
EpiMap_gene_enhancer_links_tb <- EpiMap_gene_enhancer_links_tb %>% 
	mutate(seqnames = gsub("chr", "", seqnames))

# overlap variants and EpiMap
adaptive_variants_EpiMap <- as_tibble(find_overlaps(GRanges(adaptive_variants_tb), GRanges(EpiMap_gene_enhancer_links_tb)))

# save for interact plot
# 	because want not collapsed
adaptive_variants_EpiMap_interact_all <- adaptive_variants_EpiMap
write_tsv(as_tibble(adaptive_variants_EpiMap_interact_all), gzfile("../../results/2p-annotate_variants_gene_enhancer_links/adaptive_variants_EpiMap_interact_all.txt.gz"))

# get max per EpiMap group
adaptive_variants_EpiMap <- adaptive_variants_EpiMap %>% 
	group_by(seqnames, start, end, width, strand, VariantID, VariantCHROM, VariantPOS, VariantREF, VariantALT, EpiMap_orig_group) %>% 
	arrange(EpiMap_orig_score) %>% 
	slice(n()) %>% 
	ungroup()

# save for interact plot
# 	because want not collapsed
adaptive_variants_EpiMap_interact_max <- adaptive_variants_EpiMap
write_tsv(as_tibble(adaptive_variants_EpiMap_interact_max), gzfile("../../results/2p-annotate_variants_gene_enhancer_links/adaptive_variants_EpiMap_interact_max.txt.gz"))

# EpiMap names
adaptive_variants_EpiMap_names <- setdiff(names(adaptive_variants_EpiMap)[11:length(adaptive_variants_EpiMap)], "EpiMap_orig_group")

# add summary cols
adaptive_variants_EpiMap <- adaptive_variants_EpiMap %>% 
	rowwise() %>% 
	mutate(EpiMap_summ_Score = EpiMap_orig_score) %>% 
	mutate(EpiMap_summ_Ensembl_ID = EpiMap_orig_target_ensembl_id) %>% 
	mutate(EpiMap_summ_Symbol = EpiMap_orig_target_symbol) %>% 
	mutate(EpiMap_summ_Group = EpiMap_orig_group)

# EpiMap names
adaptive_variants_EpiMap_names <- setdiff(names(adaptive_variants_EpiMap)[11:length(adaptive_variants_EpiMap)], "EpiMap_orig_group")

# collapse by group
adaptive_variants_EpiMap_collapse <- adaptive_variants_EpiMap %>% 
	group_by(seqnames, start, end, width, strand, VariantID, VariantCHROM, VariantPOS, VariantREF, VariantALT) %>% 
	dplyr::summarise_at(adaptive_variants_EpiMap_names, function(x) {paste(x, collapse=",")}) %>% 
	ungroup() %>% 
	rowwise() %>% 
	mutate(EpiMap_summ_Count = sum(!is.na(strsplit(EpiMap_orig_target_symbol, ",")[[1]]))) %>% 
	ungroup() %>% 
	mutate(EpiMap_summ_Summary = "EpiMap_gene_enhancer_linked") %>% 
	dplyr::select(seqnames, start, end, width, strand, VariantID, VariantCHROM, VariantPOS, VariantREF, VariantALT, 
		starts_with("EpiMap_orig"), EpiMap_summ_Summary, EpiMap_summ_Count, starts_with("EpiMap_summ"))

# pivot by group
adaptive_variants_EpiMap_pivot <- adaptive_variants_EpiMap %>% 
	pivot_wider(names_from=EpiMap_orig_group, names_sep="-", values_from=adaptive_variants_EpiMap_names)

# join together
adaptive_variants_EpiMap <- full_join(adaptive_variants_EpiMap_collapse, adaptive_variants_EpiMap_pivot)

# merge to full variants
adaptive_variants_tb <- as_tibble(fread("../../results/1a-preprocess_PIBv1_MPRA_pilot/adaptive_variants.txt.gz")) %>% 
	dplyr::select(seqnames, start, end, width, strand, VariantID, VariantCHROM, VariantPOS, VariantREF, VariantALT)
names_temp <- intersect(names(adaptive_variants_tb), names(adaptive_variants_EpiMap))
adaptive_variants_EpiMap <- adaptive_variants_tb %>% mutate_at(names_temp, as.character)%>% 
	left_join(adaptive_variants_EpiMap %>% mutate_at(names_temp, as.character))

# save all variants
write_tsv(adaptive_variants_EpiMap, gzfile("../../results/2p-annotate_variants_gene_enhancer_links/adaptive_variants_EpiMap.txt.gz"))

# save summ variants
adaptive_variants_EpiMap_summ <- adaptive_variants_EpiMap %>% 
	dplyr::select(seqnames, start, end, width, strand, starts_with("Variant"), contains("_summ_"))
write_tsv(adaptive_variants_EpiMap_summ, gzfile("../../results/2p-annotate_variants_gene_enhancer_links/adaptive_variants_EpiMap_summ.txt.gz"))
sort(table(adaptive_variants_EpiMap_summ$EpiMap_summ_Summary))
