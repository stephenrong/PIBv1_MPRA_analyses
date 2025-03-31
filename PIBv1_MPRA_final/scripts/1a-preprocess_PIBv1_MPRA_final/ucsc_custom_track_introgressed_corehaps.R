#!/bin/R

# for data analysis
library(tidyverse)
library(data.table)
library(paletteer)

# for bio data analysis
source("../shared_functions/create_ucsc_custom_tracks.R")

# load gr
introgressed_corehaps_ucsc_custom_track <- as_tibble(fread("../../results/1a-preprocess_PIBv1_MPRA_final/introgressed_corehaps.txt.gz"))

# create bed
custom_palette = as.vector(paletteer_c("viridis::viridis", n=length(unique(introgressed_corehaps_ucsc_custom_track$Corehaps_SprimePopulation))))

introgressed_corehaps_ucsc_custom_track_bed <- introgressed_corehaps_ucsc_custom_track %>% 
	mutate(
		bed_chrom=seqnames,
		bed_chromStart=format(start-1, scientific=FALSE),
		bed_chromEnd=format(end, scientific=FALSE),
		bed_name=paste(Corehaps_TractID, sep="-"),
		bed_score=0, 
		bed_strand=".", 
		bed_thickStart=format(start-1, scientific=FALSE), 
		bed_thickEnd=format(end, scientific=FALSE), 
		bed_itemRgb=as.factor(Corehaps_SprimePopulation)
	) %>% 
	rowwise() %>% 
	mutate(
		bed_itemRgb=paste(as.vector(col2rgb(custom_palette[as.integer(bed_itemRgb)])), collapse=",")
	) %>% 
	ungroup() %>% 
	dplyr::select(starts_with("bed"))

create_ucsc_custom_bed_track(
	input_tb=introgressed_corehaps_ucsc_custom_track_bed, 
	output_bed="../../results/1a-preprocess_PIBv1_MPRA_final/introgressed_corehaps_ucsc_custom_track.bed", 
	name="Archaic introgessed core haplotypes", description="Archaic introgessed core haplotypes", visibility="pack", itemRgb="on")
