#!/bin/R

# for data analysis
library(tidyverse)
library(data.table)
library(paletteer)

# for bio data analysis
source("../shared_functions/seqinfo_fix_change.R")
source("../shared_functions/create_ucsc_custom_tracks.R")

# load ABC_max
ABC_max_gene_enhancer_links <- as_tibble(fread("../../results/2p-annotate_variants_gene_enhancer_links/adaptive_variants_ABC_interact_max.txt.gz"))
ABC_max_gene_enhancer_links_abbrev <- ABC_max_gene_enhancer_links %>% 
	GRanges() %>% seqinfo_fix("UCSC", "hg19") %>% as_tibble()

ABC_max_gene_enhancer_links_abbrev <- ABC_max_gene_enhancer_links_abbrev %>% 
	mutate(ABC_orig_category = gsub(" ", "_", ABC_orig_category)) %>% 
	mutate(ABC_orig_category = gsub("&", "and", ABC_orig_category))

# filter to only emVars
adaptive_variants_MPRA_K562_summ <- as_tibble(fread("../../results/2i-annotate_variants_MPRA_pilot/adaptive_variants_MPRA_K562_summ.txt.gz")) %>% 
	filter(!is.na(MPRA_K562_summ_Summary))
adaptive_variants_MPRA_Jurkat_summ <- as_tibble(fread("../../results/2i-annotate_variants_MPRA_pilot/adaptive_variants_MPRA_Jurkat_summ.txt.gz")) %>% 
	filter(!is.na(MPRA_Jurkat_summ_Summary))
adaptive_variants_MPRA_HepG2_summ <- as_tibble(fread("../../results/2i-annotate_variants_MPRA_pilot/adaptive_variants_MPRA_HepG2_summ.txt.gz")) %>% 
	filter(!is.na(MPRA_HepG2_summ_Summary))
ABC_max_gene_enhancer_links_emVar <- ABC_max_gene_enhancer_links_abbrev %>% 
	filter(VariantID %in% c(adaptive_variants_MPRA_K562_summ$VariantID, adaptive_variants_MPRA_Jurkat_summ$VariantID, adaptive_variants_MPRA_HepG2_summ$VariantID))

# create table
custom_palette = as.vector(paletteer_d("RColorBrewer::Set1"))

ABC_max_gene_enhancer_links_emVar <- ABC_max_gene_enhancer_links_emVar %>% 
	rowwise() %>% 
	mutate(
		interact_chrom = seqnames, 
		interact_chromStart = format(min(start, ABC_orig_target_tss)-1, scientific=FALSE),
		interact_chromEnd = format(max(end, ABC_orig_target_tss), scientific=FALSE)
	) %>% 
	ungroup() %>% 
	mutate(
		interact_name = paste(VariantID, ABC_orig_target_symbol, ABC_orig_category, ABC_orig_biosample, sep="-"),
		interact_score = 0,
		interact_value = ABC_orig_score,
		interact_exp = ABC_orig_category,
		interact_color = as.factor(ABC_orig_category), 
		interact_sourceChrom = seqnames,
		interact_sourceStart = format(start-1, scientific=FALSE),
		interact_sourceEnd = format(end, scientific=FALSE),
		interact_sourceName = paste(VariantID, ABC_orig_target_symbol, ABC_orig_category, ABC_orig_biosample, sep="-"),
		interact_sourceStrand = "*",
		interact_targetChrom = seqnames, 
		interact_targetStart = format(ABC_orig_target_tss-1, scientific=FALSE), 
		interact_targetEnd = format(ABC_orig_target_tss, scientific=FALSE),
		interact_targetName = ABC_orig_target_symbol,
		interact_targetStrand = "*" 
	) %>% 
	rowwise() %>% 
	mutate(
		interact_color="255,0,0"
	) %>% 
	ungroup() %>% 
	dplyr::select(starts_with("interact"))

create_ucsc_custom_interact_track(
	input_tb=ABC_max_gene_enhancer_links_emVar,
	output_interact="../../results/2p-annotate_variants_gene_enhancer_links/adaptive_variants_ABC_max_emVar_ucsc_custom_track.interact", 
	name="ABC-Max+emVar gene enhancer links", description="ABC-Max+emVar gene enhancer links", interactDirectional="true", maxHeightPixels="200:100:50", visibility="full")
