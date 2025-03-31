#!/bin/R

# for data analysis
library(tidyverse)
library(data.table)
library(paletteer)

# for bio data analysis
source("../shared_functions/create_ucsc_custom_tracks.R")

# load gr
introgressed_tracts_ucsc_custom_track <- as_tibble(fread("../../results/1a-preprocess_PIBv1_MPRA_final/introgressed_tracts.txt.gz"))

# create bed
custom_palette = as.vector(paletteer_c("viridis::viridis", n=length(unique(introgressed_tracts_ucsc_custom_track$Tracts_SprimePopulation))))

introgressed_tracts_ucsc_custom_track_bed <- introgressed_tracts_ucsc_custom_track %>% 
	mutate(
		bed_chrom=seqnames,
		bed_chromStart=format(start-1, scientific=FALSE),
		bed_chromEnd=format(end, scientific=FALSE),
		bed_name=paste(Tracts_TractID, sep="-"),
		bed_score=0, 
		bed_strand=".", 
		bed_thickStart=format(start-1, scientific=FALSE), 
		bed_thickEnd=format(end, scientific=FALSE), 
		bed_itemRgb=as.factor(Tracts_SprimePopulation)
	) %>% 
	rowwise() %>% 
	mutate(
		bed_itemRgb=paste(as.vector(col2rgb(custom_palette[as.integer(bed_itemRgb)])), collapse=",")
	) %>% 
	ungroup() %>% 
	dplyr::select(starts_with("bed"))

create_ucsc_custom_bed_track(
	input_tb=introgressed_tracts_ucsc_custom_track_bed, 
	output_bed="../../results/1a-preprocess_PIBv1_MPRA_final/introgressed_tracts_ucsc_custom_track.bed", 
	name="Archaic introgressed tracts", description="Archaic introgressed tracts", visibility="pack", itemRgb="on")
