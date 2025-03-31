#!/bin/R

# for data analysis
library(tidyverse)
library(data.table)
library(paletteer)

# for bio data analysis
source("../shared_functions/create_ucsc_custom_tracks.R")

# load gr
introgressed_variants_archaic_genotypes_ucsc_custom_track <- as_tibble(fread("../../results/2e-annotate_variants_archaic_genotypes/introgressed_variants_archaic_genotypes_summ.txt.gz"))

# create bed
introgressed_variants_Altai_Neanderthal_ucsc_custom_track_bed <- introgressed_variants_archaic_genotypes_ucsc_custom_track %>% 
	mutate(
		bed_chrom=seqnames,
		bed_chromStart=format(start-1, scientific=FALSE),
		bed_chromEnd=format(end, scientific=FALSE),
		bed_name=paste(VariantID, sep="-"),
		bed_score=0, 
		bed_strand=".", 
		bed_thickStart=format(start-1, scientific=FALSE), 
		bed_thickEnd=format(end, scientific=FALSE), 
		bed_itemRgb=ifelse(archaic_genotype_summ_Altai_Neanderthal %in% c("1/1"), "0,0,255", ifelse(archaic_genotype_summ_Altai_Neanderthal %in% c("0/1", "1/0"), "255,0,0", NA))
	) %>% 
	filter(!is.na(bed_itemRgb)) %>% 
	dplyr::select(starts_with("bed"))

introgressed_variants_Altai_Denisovan_ucsc_custom_track_bed <- introgressed_variants_archaic_genotypes_ucsc_custom_track %>% 
	mutate(
		bed_chrom=seqnames,
		bed_chromStart=format(start-1, scientific=FALSE),
		bed_chromEnd=format(end, scientific=FALSE),
		bed_name=paste(VariantID, sep="-"),
		bed_score=0, 
		bed_strand=".", 
		bed_thickStart=format(start-1, scientific=FALSE), 
		bed_thickEnd=format(end, scientific=FALSE), 
		bed_itemRgb=ifelse(archaic_genotype_summ_Altai_Denisovan %in% c("1/1"), "0,0,255", ifelse(archaic_genotype_summ_Altai_Denisovan %in% c("0/1", "1/0"), "255,0,0", NA))
	) %>% 
	filter(!is.na(bed_itemRgb)) %>% 
	dplyr::select(starts_with("bed"))

introgressed_variants_Vindija_Neanderthal_ucsc_custom_track_bed <- introgressed_variants_archaic_genotypes_ucsc_custom_track %>% 
	mutate(
		bed_chrom=seqnames,
		bed_chromStart=format(start-1, scientific=FALSE),
		bed_chromEnd=format(end, scientific=FALSE),
		bed_name=paste(VariantID, sep="-"),
		bed_score=0, 
		bed_strand=".", 
		bed_thickStart=format(start-1, scientific=FALSE), 
		bed_thickEnd=format(end, scientific=FALSE), 
		bed_itemRgb=ifelse(archaic_genotype_summ_Vindija_Neanderthal %in% c("1/1"), "0,0,255", ifelse(archaic_genotype_summ_Vindija_Neanderthal %in% c("0/1", "1/0"), "255,0,0", NA))
	) %>% 
	filter(!is.na(bed_itemRgb)) %>% 
	dplyr::select(starts_with("bed"))

introgressed_variants_Chagyrskaya_Neanderthal_ucsc_custom_track_bed <- introgressed_variants_archaic_genotypes_ucsc_custom_track %>% 
	mutate(
		bed_chrom=seqnames,
		bed_chromStart=format(start-1, scientific=FALSE),
		bed_chromEnd=format(end, scientific=FALSE),
		bed_name=paste(VariantID, sep="-"),
		bed_score=0, 
		bed_strand=".", 
		bed_thickStart=format(start-1, scientific=FALSE), 
		bed_thickEnd=format(end, scientific=FALSE), 
		bed_itemRgb=ifelse(archaic_genotype_summ_Chagyrskaya_Neanderthal %in% c("1/1"), "0,0,255", ifelse(archaic_genotype_summ_Chagyrskaya_Neanderthal %in% c("0/1", "1/0"), "255,0,0", NA))
	) %>% 
	filter(!is.na(bed_itemRgb)) %>% 
	dplyr::select(starts_with("bed"))

create_ucsc_custom_bed_track(
	input_tb=introgressed_variants_Altai_Neanderthal_ucsc_custom_track_bed, 
	output_bed="../../results/2e-annotate_variants_archaic_genotypes/introgressed_variants_Altai_Neanderthal_ucsc_custom_track.bed", 
	name="Matching Altai Neanderthal", description="Matching Altai Neanderthal", visibility="dense", itemRgb="on")

create_ucsc_custom_bed_track(
	input_tb=introgressed_variants_Altai_Denisovan_ucsc_custom_track_bed, 
	output_bed="../../results/2e-annotate_variants_archaic_genotypes/introgressed_variants_Altai_Denisovan_ucsc_custom_track.bed", 
	name="Matching Altai Denisovan", description="Matching Altai Denisovan", visibility="dense", itemRgb="on")

create_ucsc_custom_bed_track(
	input_tb=introgressed_variants_Vindija_Neanderthal_ucsc_custom_track_bed, 
	output_bed="../../results/2e-annotate_variants_archaic_genotypes/introgressed_variants_Vindija_Neanderthal_ucsc_custom_track.bed", 
	name="Matching Vindija Neanderthal", description="Matching Vindija Neanderthal", visibility="dense", itemRgb="on")

create_ucsc_custom_bed_track(
	input_tb=introgressed_variants_Chagyrskaya_Neanderthal_ucsc_custom_track_bed, 
	output_bed="../../results/2e-annotate_variants_archaic_genotypes/introgressed_variants_Chagyrskaya_Neanderthal_ucsc_custom_track.bed", 
	name="Matching Chagyrskaya Neanderthal", description="Matching Chagyrskaya Neanderthal", visibility="dense", itemRgb="on")
