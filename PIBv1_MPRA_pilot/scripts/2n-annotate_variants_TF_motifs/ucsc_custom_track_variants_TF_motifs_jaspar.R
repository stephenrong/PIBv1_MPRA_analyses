#!/bin/R

# for data analysis
library(tidyverse)
library(data.table)
library(paletteer)

# for bio data analysis
source("../shared_functions/create_ucsc_custom_tracks.R")

# load gr
adaptive_variants_TF_motifs_jaspar_ucsc_custom_track <- as_tibble(fread("../../results/2n-annotate_variants_TF_motifs/adaptive_variants_TF_motifs_jaspar_summ.txt.gz"))

# create bed
custom_palette = as.vector(paletteer_d("RColorBrewer::Set1"))

adaptive_variants_TF_motifs_jaspar_ucsc_custom_track <- adaptive_variants_TF_motifs_jaspar_ucsc_custom_track %>% 
	mutate(
		bed_chrom=seqnames,
		bed_chromStart=format(start-1, scientific=FALSE),
		bed_chromEnd=format(end, scientific=FALSE),
		bed_name=paste(VariantID, TF_motifs_jaspar_summ_Description, TF_motifs_jaspar_summ_Summary, sep="-"),
		bed_score=0, 
		bed_strand=".", 
		bed_thickStart=format(start-1, scientific=FALSE), 
		bed_thickEnd=format(end, scientific=FALSE), 
		bed_itemRgb=as.factor(TF_motifs_jaspar_summ_Summary)
	) %>% 
	rowwise() %>% 
	mutate(
		bed_itemRgb=paste(as.vector(col2rgb(custom_palette[as.integer(bed_itemRgb)])), collapse=",")
	) %>% 
	ungroup() %>% 
	filter(!is.na(TF_motifs_jaspar_summ_Summary)) %>% 
	dplyr::select(starts_with("bed"))

create_ucsc_custom_bed_track(
	input_tb=adaptive_variants_TF_motifs_jaspar_ucsc_custom_track, 
	output_bed="../../results/2n-annotate_variants_TF_motifs/adaptive_variants_TF_motifs_jaspar_ucsc_custom_track.bed", 
	name="TF motif disrupting (JASPAR)", description="TF motif disrupting (JASPAR)", visibility="pack", itemRgb="on")
