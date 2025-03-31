#!/bin/R

# for data analysis
library(tidyverse)
library(data.table)
library(paletteer)

# for bio data analysis
source("../shared_functions/create_ucsc_custom_tracks.R")

# load gr
mpra_ccre_variants_ucsc_custom_track <- as_tibble(fread("../../results/1a-preprocess_PIBv1_MPRA_final/mpra_ccre_variants.txt.gz"))

# create bed
custom_color = c("#31a354")

mpra_ccre_variants_ucsc_custom_track_bed <- mpra_ccre_variants_ucsc_custom_track %>% 
	mutate(
		bed_chrom=seqnames,
		bed_chromStart=format(start-1, scientific=FALSE),
		bed_chromEnd=format(end, scientific=FALSE),
		bed_name=paste(VariantID, sep="-"),
		bed_score=0, 
		bed_strand=".", 
		bed_thickStart=format(start-1, scientific=FALSE), 
		bed_thickEnd=format(end, scientific=FALSE), 
		bed_itemRgb=paste(as.vector(col2rgb(custom_color)), collapse=",")
	) %>% 
	dplyr::select(starts_with("bed"))

create_ucsc_custom_bed_track(
	input_tb=mpra_ccre_variants_ucsc_custom_track_bed, 
	output_bed="../../results/1a-preprocess_PIBv1_MPRA_final/mpra_ccre_variants_ucsc_custom_track.bed", 
	name="Adaptive introgressed MPRA variants", description="Adaptive introgressed MPRA variants", visibility=1, itemRgb="on")
