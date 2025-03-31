#!/bin/R

# for data analysis
library(tidyverse)
library(data.table)
library(paletteer)

# for bio data analysis
source("../shared_functions/create_ucsc_custom_tracks.R")

# load gr
introgressed_variants_phyloP_ucsc_custom_track <- as_tibble(fread("../../results/2g-annotate_variants_phyloP_scores/introgressed_variants_phyloP_summ.txt.gz"))

# create bed
custom_palette = as.vector(paletteer_d("RColorBrewer::Set1"))

introgressed_variants_phyloP_ucsc_custom_track_bed <- introgressed_variants_phyloP_ucsc_custom_track %>% 
	mutate(
		bedGraph_chrom=seqnames,
		bedGraph_chromStart=format(start-1, scientific=FALSE),
		bedGraph_chromEnd=format(end, scientific=FALSE), 
		bedGraph_score=`phyloP_mam241_summ_Score`
	) %>% 
	dplyr::select(starts_with("bedGraph"))

create_ucsc_custom_bedGraph_track(
	input_tb=introgressed_variants_phyloP_ucsc_custom_track_bed, 
	output_bedGraph="../../results/2g-annotate_variants_phyloP_scores/introgressed_variants_phyloP_ucsc_custom_track.bedGraph", 
	name="phyloP (Zoonomia mammalian)", description="phyloP (Zoonomia mammalian)", visibility="full", alwaysZero="on", viewLimits=0:1, color=paste(as.vector(col2rgb("#000000")), collapse=","))
