#!/bin/R

# for data analysis
library(tidyverse)
library(data.table)
library(paletteer)

# for bio data analysis
source("../shared_functions/create_ucsc_custom_tracks.R")

# load gr
adaptive_tracts_ucsc_custom_track <- as_tibble(fread("../../results/1a-preprocess_PIBv1_MPRA_pilot/adaptive_tracts.txt.gz"))

# create bed
custom_palette = as.vector(paletteer_d("RColorBrewer::Set1"))

adaptive_tracts_ucsc_custom_track_bed <- adaptive_tracts_ucsc_custom_track %>% 
	mutate(
		bed_chrom=seqnames,
		bed_chromStart=format(start-1, scientific=FALSE),
		bed_chromEnd=format(end, scientific=FALSE),
		bed_name=paste(TractID, sep="-"),
		bed_score=0, 
		bed_strand=".", 
		bed_thickStart=format(start-1, scientific=FALSE), 
		bed_thickEnd=format(end, scientific=FALSE), 
		bed_itemRgb=as.factor(Pop)
	) %>% 
	rowwise() %>% 
	mutate(
		bed_itemRgb=paste(as.vector(col2rgb(custom_palette[as.integer(bed_itemRgb)])), collapse=",")
	) %>% 
	ungroup() %>% 
	dplyr::select(starts_with("bed"))

create_ucsc_custom_bed_track(
	input_tb=adaptive_tracts_ucsc_custom_track_bed, 
	output_bed="../../results/1a-preprocess_PIBv1_MPRA_pilot/adaptive_tracts_ucsc_custom_track.bed", 
	name="Adaptive introgressed tracts", description="Adaptive introgressed tracts", visibility="pack", itemRgb="on")
