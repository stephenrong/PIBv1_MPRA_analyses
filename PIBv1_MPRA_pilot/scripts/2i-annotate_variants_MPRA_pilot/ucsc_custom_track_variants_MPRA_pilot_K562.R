#!/bin/R

# for data analysis
library(tidyverse)
library(data.table)
library(paletteer)

# for bio data analysis
source("../shared_functions/create_ucsc_custom_tracks.R")

# load gr
adaptive_variants_MPRA_K562_ucsc_custom_track <- as_tibble(fread("../../results/2i-annotate_variants_MPRA_pilot/adaptive_variants_MPRA_K562_summ.txt.gz"))
adaptive_variants_MPRA_K562_ucsc_custom_track <- adaptive_variants_MPRA_K562_ucsc_custom_track %>% 
	filter(!is.na(MPRA_K562_summ_Description_activity_log10FDR))

# create bed
custom_palette = as.vector(paletteer_d("RColorBrewer::Set1"))

adaptive_variants_MPRA_K562_variant_ucsc_custom_track <- adaptive_variants_MPRA_K562_ucsc_custom_track %>% 
	filter(!is.na(MPRA_K562_summ_Summary)) %>% 
	mutate(
		bed_chrom=seqnames,
		bed_chromStart=format(start-1, scientific=FALSE),
		bed_chromEnd=format(end, scientific=FALSE),
		bed_name=paste(VariantID, MPRA_K562_summ_Summary, sep="-"),
		bed_score=0, 
		bed_strand=".", 
		bed_thickStart=format(start-1, scientific=FALSE), 
		bed_thickEnd=format(end, scientific=FALSE), 
		bed_itemRgb=as.factor(MPRA_K562_summ_Summary)
	) %>% 
	rowwise() %>% 
	mutate(
		bed_itemRgb=paste(as.vector(col2rgb(custom_palette[as.integer(bed_itemRgb)])), collapse=",")
	) %>% 
	ungroup() %>% 
	filter(!is.na(MPRA_K562_summ_Summary)) %>% 
	dplyr::select(starts_with("bed"))

create_ucsc_custom_bed_track(
	input_tb=adaptive_variants_MPRA_K562_variant_ucsc_custom_track, 
	output_bed="../../results/2i-annotate_variants_MPRA_pilot/adaptive_variants_MPRA_K562_variant_ucsc_custom_track.bed", 
	name="MRPA K562 emVar", description="MRPA K562 emVar", visibility="pack", itemRgb="on")

adaptive_variants_MPRA_K562_logSkew_ucsc_custom_track_bed <- adaptive_variants_MPRA_K562_ucsc_custom_track %>% 
	mutate(
		bedGraph_chrom=seqnames,
		bedGraph_chromStart=format(start-1, scientific=FALSE),
		bedGraph_chromEnd=format(end, scientific=FALSE), 
		bedGraph_score=`MPRA_K562_summ_Description_skew_log2FC`
	) %>% 
	dplyr::select(starts_with("bedGraph"))

create_ucsc_custom_bedGraph_track(
	input_tb=adaptive_variants_MPRA_K562_logSkew_ucsc_custom_track_bed, 
	output_bedGraph="../../results/2i-annotate_variants_MPRA_pilot/adaptive_variants_MPRA_K562_logSkew_ucsc_custom_track.bedGraph", 
	name="MPRA K562 log2Skew", description="MPRA K562 log2Skew", visibility="full", alwaysZero="on", viewLimits=0:1, color=paste(as.vector(col2rgb("#000000")), collapse=","))

adaptive_variants_MPRA_K562_logFDR_ucsc_custom_track_bed <- adaptive_variants_MPRA_K562_ucsc_custom_track %>% 
	mutate(
		bedGraph_chrom=seqnames,
		bedGraph_chromStart=format(start-1, scientific=FALSE),
		bedGraph_chromEnd=format(end, scientific=FALSE), 
		bedGraph_score=`MPRA_K562_summ_Description_skew_log10FDR`
	) %>% 
	dplyr::select(starts_with("bedGraph"))

create_ucsc_custom_bedGraph_track(
	input_tb=adaptive_variants_MPRA_K562_logFDR_ucsc_custom_track_bed, 
	output_bedGraph="../../results/2i-annotate_variants_MPRA_pilot/adaptive_variants_MPRA_K562_logFDR_ucsc_custom_track.bedGraph", 
	name="MPRA K562 log10FDR", description="MPRA K562 log10FDR", visibility="full", alwaysZero="on", viewLimits=0:1, color=paste(as.vector(col2rgb("#000000")), collapse=","))
