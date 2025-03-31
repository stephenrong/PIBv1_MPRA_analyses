#!/bin/R

# for data analysis
library(tidyverse)
library(data.table)
library(paletteer)

# for bio data analysis
source("../shared_functions/create_ucsc_custom_tracks.R")

# load gr
adaptive_variants_BBJ_finemap_ucsc_custom_track <- as_tibble(fread("../../results/2m-annotate_variants_UKBB_BBJ_finemap/adaptive_variants_BBJ_finemap_summ.txt.gz"))
adaptive_variants_BBJ_finemap_ucsc_custom_track <- adaptive_variants_BBJ_finemap_ucsc_custom_track %>% 
	filter(!is.na(BBJ_finemap_summ_Summary))

# create bed
custom_palette = as.vector(paletteer_d("RColorBrewer::Set1"))

adaptive_variants_BBJ_finemap_ucsc_custom_track_bed <- adaptive_variants_BBJ_finemap_ucsc_custom_track %>% 
	mutate(
		bed_chrom=seqnames,
		bed_chromStart=format(start-1, scientific=FALSE),
		bed_chromEnd=format(end, scientific=FALSE),
		bed_name=paste(VariantID, BBJ_finemap_summ_Summary, sep="-"),
		bed_score=0, 
		bed_strand=".", 
		bed_thickStart=format(start-1, scientific=FALSE), 
		bed_thickEnd=format(end, scientific=FALSE), 
		bed_itemRgb=as.factor(BBJ_finemap_summ_Summary)
	) %>% 
	rowwise() %>% 
	mutate(
		bed_itemRgb=paste(as.vector(col2rgb(custom_palette[1])), collapse=",")
	) %>% 
	ungroup() %>% 
	dplyr::select(starts_with("bed"))

create_ucsc_custom_bed_track(
	input_tb=adaptive_variants_BBJ_finemap_ucsc_custom_track_bed, 
	output_bed="../../results/2m-annotate_variants_UKBB_BBJ_finemap/adaptive_variants_BBJ_finemap_ucsc_custom_track.bed", 
	name="BBJ fine mapped PIP>=0.1", description="BBJ fine mapped PIP>=0.1", visibility="pack", itemRgb="on")
