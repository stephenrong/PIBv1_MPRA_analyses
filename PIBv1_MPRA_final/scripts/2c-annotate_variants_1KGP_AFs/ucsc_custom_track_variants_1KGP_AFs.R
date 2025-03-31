#!/bin/R

# for data analysis
library(tidyverse)
library(data.table)
library(paletteer)

# for bio data analysis
source("../shared_functions/create_ucsc_custom_tracks.R")

# load gr
introgressed_variants_1KGP_phase3_AFs_ucsc_custom_track <- as_tibble(fread("../../results/2c-annotate_variants_1KGP_AFs/introgressed_variants_1KGP_phase3_AFs_summ.txt.gz"))
introgressed_variants_1KGP_phase3_AFs_ucsc_custom_track <- introgressed_variants_1KGP_phase3_AFs_ucsc_custom_track %>% 
	filter(!is.na(`1KGP_phase3_summ_AF`))

# create bed
introgressed_variants_1KGP_phase3_AFs_ucsc_custom_track_bed <- introgressed_variants_1KGP_phase3_AFs_ucsc_custom_track %>% 
	mutate(
		bedGraph_chrom=seqnames,
		bedGraph_chromStart=format(start-1, scientific=FALSE),
		bedGraph_chromEnd=format(end, scientific=FALSE), 
		bedGraph_score=`1KGP_phase3_summ_AF`
	) %>% 
	dplyr::select(starts_with("bedGraph"))

create_ucsc_custom_bedGraph_track(
	input_tb=introgressed_variants_1KGP_phase3_AFs_ucsc_custom_track_bed, 
	output_bedGraph="../../results/2c-annotate_variants_1KGP_AFs/introgressed_variants_1KGP_phase3_AFs_ucsc_custom_track.bedGraph", 
	name="1KGP AFs", description="1KGP AFs", visibility="full", alwaysZero="on", viewLimits=0:1, color=paste(as.vector(col2rgb("#000000")), collapse=","))
