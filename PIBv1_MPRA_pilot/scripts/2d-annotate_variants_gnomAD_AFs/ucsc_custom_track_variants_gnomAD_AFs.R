#!/bin/R

# for data analysis
library(tidyverse)
library(data.table)
library(paletteer)

# for bio data analysis
source("../shared_functions/create_ucsc_custom_tracks.R")

# load gr
adaptive_variants_gnomAD_v3.1.2_AFs_ucsc_custom_track <- as_tibble(fread("../../results/2d-annotate_variants_gnomAD_AFs/adaptive_variants_gnomAD_v3.1.2_AFs_summ.txt.gz"))
adaptive_variants_gnomAD_v3.1.2_AFs_ucsc_custom_track <- adaptive_variants_gnomAD_v3.1.2_AFs_ucsc_custom_track %>% 
	filter(!is.na(`gnomAD_v3.1.2_summ_AF`))

adaptive_variants_gnomAD_v3.1.2_AFs_ucsc_custom_track_bed <- adaptive_variants_gnomAD_v3.1.2_AFs_ucsc_custom_track %>% 
	mutate(
		bedGraph_chrom=seqnames,
		bedGraph_chromStart=format(start-1, scientific=FALSE),
		bedGraph_chromEnd=format(end, scientific=FALSE), 
		bedGraph_score=`gnomAD_v3.1.2_summ_AF`
	) %>% 
	dplyr::select(starts_with("bedGraph"))

create_ucsc_custom_bedGraph_track(
	input_tb=adaptive_variants_gnomAD_v3.1.2_AFs_ucsc_custom_track_bed, 
	output_bedGraph="../../results/2d-annotate_variants_gnomAD_AFs/adaptive_variants_gnomAD_v3.1.2_AFs_ucsc_custom_track.bedGraph", 
	name="gnomAD AFs", description="gnomAD AFs", visibility="full", alwaysZero="on", viewLimits=0:1, color=paste(as.vector(col2rgb("#000000")), collapse=","))
