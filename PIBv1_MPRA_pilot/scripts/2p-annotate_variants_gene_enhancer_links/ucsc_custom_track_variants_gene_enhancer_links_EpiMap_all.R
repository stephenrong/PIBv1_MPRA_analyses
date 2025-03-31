#!/bin/R

# for data analysis
library(tidyverse)
library(data.table)
library(paletteer)

# for bio data analysis
source("../shared_functions/seqinfo_fix_change.R")
source("../shared_functions/create_ucsc_custom_tracks.R")

# load EpiMap_all
EpiMap_all_gene_enhancer_links <- as_tibble(fread("../../results/2p-annotate_variants_gene_enhancer_links/adaptive_variants_EpiMap_interact_all.txt.gz"))
EpiMap_all_gene_enhancer_links_abbrev <- EpiMap_all_gene_enhancer_links %>% 
	GRanges() %>% seqinfo_fix("UCSC", "hg19") %>% as_tibble()

EpiMap_all_gene_enhancer_links_abbrev <- EpiMap_all_gene_enhancer_links_abbrev %>% 
	mutate(EpiMap_orig_group = gsub(" ", "_", EpiMap_orig_group)) %>% 
	mutate(EpiMap_orig_group = gsub("&", "and", EpiMap_orig_group))

# filter to only emVars
adaptive_variants_MPRA_K562_summ <- as_tibble(fread("../../results/2i-annotate_variants_MPRA_pilot/adaptive_variants_MPRA_K562_summ.txt.gz")) %>% 
	filter(!is.na(MPRA_K562_summ_Summary))
adaptive_variants_MPRA_Jurkat_summ <- as_tibble(fread("../../results/2i-annotate_variants_MPRA_pilot/adaptive_variants_MPRA_Jurkat_summ.txt.gz")) %>% 
	filter(!is.na(MPRA_Jurkat_summ_Summary))
adaptive_variants_MPRA_HepG2_summ <- as_tibble(fread("../../results/2i-annotate_variants_MPRA_pilot/adaptive_variants_MPRA_HepG2_summ.txt.gz")) %>% 
	filter(!is.na(MPRA_HepG2_summ_Summary))
EpiMap_all_gene_enhancer_links_emVar <- EpiMap_all_gene_enhancer_links_abbrev %>% 
	filter(VariantID %in% c(adaptive_variants_MPRA_K562_summ$VariantID, adaptive_variants_MPRA_Jurkat_summ$VariantID, adaptive_variants_MPRA_HepG2_summ$VariantID))

# create table
custom_palette = as.vector(paletteer_d("RColorBrewer::Set1"))

EpiMap_all_gene_enhancer_links_emVar <- EpiMap_all_gene_enhancer_links_emVar %>% 
	rowwise() %>% 
	mutate(
		interact_chrom = seqnames, 
		interact_chromStart = format(min(start, EpiMap_orig_target_tss)-1, scientific=FALSE),
		interact_chromEnd = format(max(end, EpiMap_orig_target_tss), scientific=FALSE)
	) %>% 
	ungroup() %>% 
	mutate(
		interact_name = paste(VariantID, EpiMap_orig_target_symbol, EpiMap_orig_group, sep="-"),
		interact_score = 0,
		interact_value = EpiMap_orig_score,
		interact_exp = EpiMap_orig_group,
		interact_color = as.factor(EpiMap_orig_group), 
		interact_sourceChrom = seqnames,
		interact_sourceStart = format(start-1, scientific=FALSE),
		interact_sourceEnd = format(end, scientific=FALSE),
		interact_sourceName = paste(VariantID, EpiMap_orig_target_symbol, EpiMap_orig_group, sep="-"),
		interact_sourceStrand = "*",
		interact_targetChrom = seqnames, 
		interact_targetStart = format(EpiMap_orig_target_tss-1, scientific=FALSE), 
		interact_targetEnd = format(EpiMap_orig_target_tss, scientific=FALSE),
		interact_targetName = EpiMap_orig_target_symbol,
		interact_targetStrand = "*" 
	) %>% 
	rowwise() %>% 
	mutate(
		interact_color="169,169,169"
	) %>% 
	ungroup() %>% 
	dplyr::select(starts_with("interact"))

create_ucsc_custom_interact_track(
	input_tb=EpiMap_all_gene_enhancer_links_emVar,
	output_interact="../../results/2p-annotate_variants_gene_enhancer_links/adaptive_variants_EpiMap_all_emVar_ucsc_custom_track.interact", 
	name="EpiMap-All+emVar gene enhancer links", description="EpiMap-All+emVar gene enhancer links", interactDirectional="true", maxHeightPixels="200:100:50", visibility="full")
