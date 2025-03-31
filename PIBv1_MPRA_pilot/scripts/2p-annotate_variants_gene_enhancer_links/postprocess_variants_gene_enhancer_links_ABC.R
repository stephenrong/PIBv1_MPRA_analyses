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

# load ABC
ABC_gene_enhancer_links_tb <- as_tibble(fread(paste0("../../../../Datasets/gene_regulation_gene_enhancer_links/ABC_gene_enhancer_links/data_cleanup/GRCh37_merged/ABC_gene_enhancer_links.txt.gz")))

# remove ambig TSS
ABC_gene_enhancer_links_tb <- ABC_gene_enhancer_links_tb %>% 
	filter(!is.na(target_symbol)) %>% 
	filter(!is.na(target_tss))

# ABC names
names(ABC_gene_enhancer_links_tb)[6:length(ABC_gene_enhancer_links_tb)] <- paste("ABC_orig", names(ABC_gene_enhancer_links_tb)[6:length(ABC_gene_enhancer_links_tb)], sep="_")
adaptive_variants_ABC_names <- names(ABC_gene_enhancer_links_tb)[6:length(ABC_gene_enhancer_links_tb)]

# convert to GRange
ABC_gene_enhancer_links_tb <- ABC_gene_enhancer_links_tb %>% 
	mutate(seqnames = gsub("chr", "", seqnames))

# overlap variants and ABC
adaptive_variants_ABC <- as_tibble(find_overlaps(GRanges(adaptive_variants_tb), GRanges(ABC_gene_enhancer_links_tb)))

# save for interact plot
# 	because want not collapsed
adaptive_variants_ABC_interact_all <- adaptive_variants_ABC
write_tsv(as_tibble(adaptive_variants_ABC_interact_all), gzfile("../../results/2p-annotate_variants_gene_enhancer_links/adaptive_variants_ABC_interact_all.txt.gz"))

# get max per ABC category
adaptive_variants_ABC <- adaptive_variants_ABC %>% 
	group_by(seqnames, start, end, width, strand, VariantID, VariantCHROM, VariantPOS, VariantREF, VariantALT, ABC_orig_category) %>% 
	arrange(ABC_orig_score) %>% 
	slice(n()) %>% 
	ungroup()

# save for interact plot
# 	because want not collapsed
adaptive_variants_ABC_interact_max <- adaptive_variants_ABC
write_tsv(as_tibble(adaptive_variants_ABC_interact_max), gzfile("../../results/2p-annotate_variants_gene_enhancer_links/adaptive_variants_ABC_interact_max.txt.gz"))

# ABC names
adaptive_variants_ABC_names <- setdiff(names(adaptive_variants_ABC)[11:length(adaptive_variants_ABC)], "ABC_orig_category")

# add summary cols
adaptive_variants_ABC <- adaptive_variants_ABC %>% 
	rowwise() %>% 
	mutate(ABC_summ_Score = ABC_orig_score) %>% 
	mutate(ABC_summ_Gene = ABC_orig_target_gene) %>% 
	mutate(ABC_summ_Symbol = ABC_orig_target_symbol) %>% 
	mutate(ABC_summ_Biosample = ABC_orig_biosample) %>% 
	mutate(ABC_summ_Category = ABC_orig_category)

# ABC names
adaptive_variants_ABC_names <- setdiff(names(adaptive_variants_ABC)[11:length(adaptive_variants_ABC)], "ABC_orig_category")

# collapse by category
adaptive_variants_ABC_collapse <- adaptive_variants_ABC %>% 
	group_by(seqnames, start, end, width, strand, VariantID, VariantCHROM, VariantPOS, VariantREF, VariantALT) %>% 
	dplyr::summarise_at(adaptive_variants_ABC_names, function(x) {paste(x, collapse=",")}) %>% 
	ungroup() %>% 
	rowwise() %>% 
	mutate(ABC_summ_Count = sum(!is.na(strsplit(ABC_orig_target_symbol, ",")[[1]]))) %>% 
	ungroup() %>% 
	mutate(ABC_summ_Summary = "ABC_gene_enhancer_linked") %>% 
	dplyr::select(seqnames, start, end, width, strand, VariantID, VariantCHROM, VariantPOS, VariantREF, VariantALT, 
		starts_with("ABC_orig"), ABC_summ_Summary, ABC_summ_Count, starts_with("ABC_summ"))

# pivot by category
adaptive_variants_ABC_pivot <- adaptive_variants_ABC %>% 
	pivot_wider(names_from=ABC_orig_category, names_sep="-", values_from=adaptive_variants_ABC_names)

# join together
adaptive_variants_ABC <- full_join(adaptive_variants_ABC_collapse, adaptive_variants_ABC_pivot)

# merge to full variants
adaptive_variants_tb <- as_tibble(fread("../../results/1a-preprocess_PIBv1_MPRA_pilot/adaptive_variants.txt.gz")) %>% 
	dplyr::select(seqnames, start, end, width, strand, VariantID, VariantCHROM, VariantPOS, VariantREF, VariantALT)
names_temp <- intersect(names(adaptive_variants_tb), names(adaptive_variants_ABC))
adaptive_variants_ABC <- adaptive_variants_tb %>% mutate_at(names_temp, as.character)%>% 
	left_join(adaptive_variants_ABC %>% mutate_at(names_temp, as.character))

# save all variants
write_tsv(adaptive_variants_ABC, gzfile("../../results/2p-annotate_variants_gene_enhancer_links/adaptive_variants_ABC.txt.gz"))

# save summ variants
adaptive_variants_ABC_summ <- adaptive_variants_ABC %>% 
	dplyr::select(seqnames, start, end, width, strand, starts_with("Variant"), contains("_summ_"))
write_tsv(adaptive_variants_ABC_summ, gzfile("../../results/2p-annotate_variants_gene_enhancer_links/adaptive_variants_ABC_summ.txt.gz"))
sort(table(adaptive_variants_ABC_summ$ABC_summ_Summary))
