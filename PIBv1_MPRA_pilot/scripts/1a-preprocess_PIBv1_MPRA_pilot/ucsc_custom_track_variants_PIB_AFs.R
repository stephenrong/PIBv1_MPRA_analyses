#!/bin/R

# for data analysis
library(tidyverse)
library(data.table)

# for bio data analysis
source("../shared_functions/create_ucsc_custom_tracks.R")

# load gr
adaptive_variants_PIB_AFs_ucsc_custom_track <- as_tibble(fread("../../results/1a-preprocess_PIBv1_MPRA_pilot/adaptive_variants.txt.gz"))
adaptive_variants_PIB_AFs_ucsc_custom_track_bed <- adaptive_variants_PIB_AFs_ucsc_custom_track %>% 
	mutate(
		bedGraph_chrom=seqnames,
		bedGraph_chromStart=format(start-1, scientific=FALSE),
		bedGraph_chromEnd=format(end, scientific=FALSE), 
		bedGraph_score=`PIB_AF`
	) %>% 
	dplyr::select(starts_with("bedGraph"))

adaptive_variants_Ata_AFs_ucsc_custom_track_bed <- adaptive_variants_PIB_AFs_ucsc_custom_track %>% 
	mutate(
		bedGraph_chrom=seqnames,
		bedGraph_chromStart=format(start-1, scientific=FALSE),
		bedGraph_chromEnd=format(end, scientific=FALSE), 
		bedGraph_score=`Ata_AF`
	) %>% 
	dplyr::select(starts_with("bedGraph"))

adaptive_variants_Baining_AFs_ucsc_custom_track_bed <- adaptive_variants_PIB_AFs_ucsc_custom_track %>% 
	mutate(
		bedGraph_chrom=seqnames,
		bedGraph_chromStart=format(start-1, scientific=FALSE),
		bedGraph_chromEnd=format(end, scientific=FALSE), 
		bedGraph_score=`Baining_AF`
	) %>% 
	dplyr::select(starts_with("bedGraph"))

adaptive_variants_Mamusi_AFs_ucsc_custom_track_bed <- adaptive_variants_PIB_AFs_ucsc_custom_track %>% 
	mutate(
		bedGraph_chrom=seqnames,
		bedGraph_chromStart=format(start-1, scientific=FALSE),
		bedGraph_chromEnd=format(end, scientific=FALSE), 
		bedGraph_score=`Mamusi_AF`
	) %>% 
	dplyr::select(starts_with("bedGraph"))

adaptive_variants_Melamela_AFs_ucsc_custom_track_bed <- adaptive_variants_PIB_AFs_ucsc_custom_track %>% 
	mutate(
		bedGraph_chrom=seqnames,
		bedGraph_chromStart=format(start-1, scientific=FALSE),
		bedGraph_chromEnd=format(end, scientific=FALSE), 
		bedGraph_score=`Melamela_AF`
	) %>% 
	dplyr::select(starts_with("bedGraph"))

adaptive_variants_Lavongai_AFs_ucsc_custom_track_bed <- adaptive_variants_PIB_AFs_ucsc_custom_track %>% 
	mutate(
		bedGraph_chrom=seqnames,
		bedGraph_chromStart=format(start-1, scientific=FALSE),
		bedGraph_chromEnd=format(end, scientific=FALSE), 
		bedGraph_score=`Lavongai_AF`
	) %>% 
	dplyr::select(starts_with("bedGraph"))

create_ucsc_custom_bedGraph_track(
	input_tb=adaptive_variants_PIB_AFs_ucsc_custom_track_bed, 
	output_bedGraph="../../results/1a-preprocess_PIBv1_MPRA_pilot/adaptive_variants_PIB_AFs_ucsc_custom_track.bedGraph", 
	name="PIB AFs", description="PIB AFs", visibility="full", alwaysZero="on", viewLimits=0:1, color=paste(as.vector(col2rgb("#000000")), collapse=","))

create_ucsc_custom_bedGraph_track(
	input_tb=adaptive_variants_Ata_AFs_ucsc_custom_track_bed, 
	output_bedGraph="../../results/1a-preprocess_PIBv1_MPRA_pilot/adaptive_variants_Ata_AFs_ucsc_custom_track.bedGraph", 
	name="Ata AFs", description="Ata AFs", visibility="full", alwaysZero="on", viewLimits=0:1, color=paste(as.vector(col2rgb("#000000")), collapse=","))

create_ucsc_custom_bedGraph_track(
	input_tb=adaptive_variants_Baining_AFs_ucsc_custom_track_bed, 
	output_bedGraph="../../results/1a-preprocess_PIBv1_MPRA_pilot/adaptive_variants_Baining_AFs_ucsc_custom_track.bedGraph", 
	name="Baining AFs", description="Baining AFs", visibility="full", alwaysZero="on", viewLimits=0:1, color=paste(as.vector(col2rgb("#000000")), collapse=","))

create_ucsc_custom_bedGraph_track(
	input_tb=adaptive_variants_Mamusi_AFs_ucsc_custom_track_bed, 
	output_bedGraph="../../results/1a-preprocess_PIBv1_MPRA_pilot/adaptive_variants_Mamusi_AFs_ucsc_custom_track.bedGraph", 
	name="Mamusi AFs", description="Mamusi AFs", visibility="full", alwaysZero="on", viewLimits=0:1, color=paste(as.vector(col2rgb("#000000")), collapse=","))

create_ucsc_custom_bedGraph_track(
	input_tb=adaptive_variants_Melamela_AFs_ucsc_custom_track_bed, 
	output_bedGraph="../../results/1a-preprocess_PIBv1_MPRA_pilot/adaptive_variants_Melamela_AFs_ucsc_custom_track.bedGraph", 
	name="Melamela AFs", description="Melamela AFs", visibility="full", alwaysZero="on", viewLimits=0:1, color=paste(as.vector(col2rgb("#000000")), collapse=","))

create_ucsc_custom_bedGraph_track(
	input_tb=adaptive_variants_Lavongai_AFs_ucsc_custom_track_bed, 
	output_bedGraph="../../results/1a-preprocess_PIBv1_MPRA_pilot/adaptive_variants_Lavongai_AFs_ucsc_custom_track.bedGraph", 
	name="Lavongai AFs", description="Lavongai AFs", visibility="full", alwaysZero="on", viewLimits=0:1, color=paste(as.vector(col2rgb("#000000")), collapse=","))
