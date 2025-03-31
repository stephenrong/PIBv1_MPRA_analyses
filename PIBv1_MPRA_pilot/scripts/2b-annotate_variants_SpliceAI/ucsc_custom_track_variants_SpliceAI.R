#!/bin/R

# for data analysis
library(tidyverse)
library(data.table)
library(paletteer)

# for bio data analysis
source("../shared_functions/create_ucsc_custom_tracks.R")

# load gr
adaptive_variants_SpliceAI_raw_ucsc_custom_track <- as_tibble(fread("../../results/2b-annotate_variants_SpliceAI/adaptive_variants_SpliceAI_summ.txt.gz"))
adaptive_variants_SpliceAI_raw_ucsc_custom_track <- adaptive_variants_SpliceAI_raw_ucsc_custom_track %>% 
	filter(!is.na(SpliceAI_raw_summ_Summary))

# create bed
custom_palette = as.vector(paletteer_d("RColorBrewer::Set1"))

adaptive_variants_SpliceAI_raw_ucsc_custom_track_bed <- adaptive_variants_SpliceAI_raw_ucsc_custom_track %>% 
	mutate(
		bed_chrom=seqnames,
		bed_chromStart=format(start-1, scientific=FALSE),
		bed_chromEnd=format(end, scientific=FALSE),
		bed_name=paste(VariantID, SpliceAI_raw_summ_Symbol, SpliceAI_raw_summ_Summary, sep="-"),
		bed_score=0, 
		bed_strand=".", 
		bed_thickStart=format(start-1, scientific=FALSE), 
		bed_thickEnd=format(end, scientific=FALSE), 
		bed_itemRgb=as.factor(SpliceAI_raw_summ_Summary)
	) %>% 
	rowwise() %>% 
	mutate(
		bed_itemRgb=paste(as.vector(col2rgb(custom_palette[as.integer(bed_itemRgb)])), collapse=",")
	) %>% 
	ungroup() %>% 
	dplyr::select(starts_with("bed"))

create_ucsc_custom_bed_track(
	input_tb=adaptive_variants_SpliceAI_raw_ucsc_custom_track_bed, 
	output_bed="../../results/2b-annotate_variants_SpliceAI/adaptive_variants_SpliceAI_raw_ucsc_custom_track.bed", 
	name="SpliceAI max score", description="SpliceAI max score", visibility="pack", itemRgb="on")
