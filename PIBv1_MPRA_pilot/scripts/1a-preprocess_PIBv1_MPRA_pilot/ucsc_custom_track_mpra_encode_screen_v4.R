#!/bin/R

# for data analysis
library(tidyverse)
library(data.table)
library(paletteer)

# for bio data analysis
library(rtracklayer)
source("../shared_functions/create_ucsc_custom_tracks.R")

# load gr
mpra_encode_screen_v4_ucsc_custom_track <- as_tibble(fread("../../../Datasets/gene_regulation_element_catalogs/ENCODE_SCREEN_V4_cCREs/data_cleanup/sample_agnostic/ENCODE_SCREEN_V4_cCREs_lift37/GRCh38-cCREs.V4.lift37.nochr.sort.bed.gz", header=FALSE))

# create bed
custom_palette = as.vector(paletteer_d("RColorBrewer::Set3"))

mpra_encode_screen_v4_ucsc_custom_track_bed <- mpra_encode_screen_v4_ucsc_custom_track %>% 
	mutate(
		bed_chrom=V1,
		bed_chromStart=format(V2, scientific=FALSE),
		bed_chromEnd=format(V3, scientific=FALSE),
		bed_name=paste(V6, sep="-"),
		bed_score=0, 
		bed_strand=".", 
		bed_thickStart=format(V2, scientific=FALSE), 
		bed_thickEnd=format(V3, scientific=FALSE), 
		bed_itemRgb=as.factor(V6)
	) %>% 
	rowwise() %>% 
	mutate(
		bed_itemRgb=paste(as.vector(col2rgb(custom_palette[as.integer(bed_itemRgb)])), collapse=",")
	) %>% 
	ungroup() %>% 
	dplyr::select(starts_with("bed"))

create_ucsc_custom_bed_track(
	input_tb=mpra_encode_screen_v4_ucsc_custom_track_bed, 
	output_bed="../../results/1a-preprocess_PIBv1_MPRA_pilot/mpra_encode_screen_v4_ucsc_custom_track.bed", 
	name="ENCODE SCREEN v4 cCREs", description="ENCODE SCREEN v4 cCREs", visibility=4, itemRgb="on")

# split versions
create_ucsc_custom_bed_track(
	input_tb=(mpra_encode_screen_v4_ucsc_custom_track_bed %>% filter(grepl("^PLS", bed_name))), 
	output_bed="../../results/1a-preprocess_PIBv1_MPRA_pilot/mpra_encode_screen_v4_ucsc_custom_track_PLS.bed", 
	name="ENCODE SCREEN v4 cCREs (PLS)", description="ENCODE SCREEN v4 cCREs (PLS)", visibility=4, itemRgb="on")

create_ucsc_custom_bed_track(
	input_tb=(mpra_encode_screen_v4_ucsc_custom_track_bed %>% filter(grepl("^pELS", bed_name))), 
	output_bed="../../results/1a-preprocess_PIBv1_MPRA_pilot/mpra_encode_screen_v4_ucsc_custom_track_pELS.bed", 
	name="ENCODE SCREEN v4 cCREs (pELS)", description="ENCODE SCREEN v4 cCREs (pELS)", visibility=4, itemRgb="on")

create_ucsc_custom_bed_track(
	input_tb=(mpra_encode_screen_v4_ucsc_custom_track_bed %>% filter(grepl("^dELS", bed_name))), 
	output_bed="../../results/1a-preprocess_PIBv1_MPRA_pilot/mpra_encode_screen_v4_ucsc_custom_track_dELS.bed", 
	name="ENCODE SCREEN v4 cCREs (dELS)", description="ENCODE SCREEN v4 cCREs (dELS)", visibility=4, itemRgb="on")

create_ucsc_custom_bed_track(
	input_tb=(mpra_encode_screen_v4_ucsc_custom_track_bed %>% filter(grepl("^CA|^TF", bed_name))), 
	output_bed="../../results/1a-preprocess_PIBv1_MPRA_pilot/mpra_encode_screen_v4_ucsc_custom_track_CA-TF.bed", 
	name="ENCODE SCREEN v4 cCREs (CA/TF)", description="ENCODE SCREEN v4 cCREs (CA/TF)", visibility=1, itemRgb="on")
