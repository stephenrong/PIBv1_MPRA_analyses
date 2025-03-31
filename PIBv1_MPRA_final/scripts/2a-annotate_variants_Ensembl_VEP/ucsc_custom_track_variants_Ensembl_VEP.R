#!/bin/R

# for data analysis
library(tidyverse)
library(data.table)
library(paletteer)

# for bio data analysis
source("../shared_functions/create_ucsc_custom_tracks.R")

# load gr
introgressed_variants_Ensembl_VEP <- as_tibble(fread("../../results/2a-annotate_variants_Ensembl_VEP/introgressed_variants_Ensembl_VEP_summ.txt.gz"))

# create bed
custom_palette = as.vector(paletteer_d("RColorBrewer::Set1"))

introgressed_variants_Ensembl_VEP_ucsc_custom_track <- introgressed_variants_Ensembl_VEP %>% 
	mutate(
		bed_chrom=seqnames,
		bed_chromStart=format(start-1, scientific=FALSE),
		bed_chromEnd=format(end, scientific=FALSE),
		bed_name=paste(VariantID, Ensembl_VEP_summ_Symbol, Ensembl_VEP_summ_Summary, sep="-"),
		bed_score=0, 
		bed_strand=".", 
		bed_thickStart=format(start-1, scientific=FALSE), 
		bed_thickEnd=format(end, scientific=FALSE), 
		bed_itemRgb=as.factor(Ensembl_VEP_summ_Summary)
	) %>% 
	rowwise() %>% 
	mutate(
		bed_itemRgb=paste(as.vector(col2rgb(custom_palette[as.integer(bed_itemRgb)])), collapse=",")
	) %>% 
	ungroup() %>% 
	filter(!is.na(Ensembl_VEP_summ_Summary)) %>% 
	dplyr::select(starts_with("bed"))

create_ucsc_custom_bed_track(
	input_tb=introgressed_variants_Ensembl_VEP_ucsc_custom_track, 
	output_bed="../../results/2a-annotate_variants_Ensembl_VEP/introgressed_variants_Ensembl_VEP_ucsc_custom_track.bed", 
	name="Ensembl VEP Consequence", description="Ensembl VEP Consequence", visibility="pack", itemRgb="on")
