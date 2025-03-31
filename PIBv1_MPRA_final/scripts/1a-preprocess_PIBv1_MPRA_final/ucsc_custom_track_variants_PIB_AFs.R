#!/bin/R

# for data analysis
library(tidyverse)
library(data.table)

# for bio data analysis
source("../shared_functions/create_ucsc_custom_tracks.R")

# load gr
introgressed_variants_PIB_SAFs_ucsc_custom_track <- as_tibble(fread("../../results/1a-preprocess_PIBv1_MPRA_final/introgressed_variants.txt.gz"))

# create bed
introgressed_variants_Ata_SAFs_ucsc_custom_track <- introgressed_variants_PIB_SAFs_ucsc_custom_track %>% 
	mutate(
		bedGraph_chrom=seqnames,
		bedGraph_chromStart=format(start-1, scientific=FALSE),
		bedGraph_chromEnd=format(end, scientific=FALSE), 
		bedGraph_score=`SAF_Ata`
	) %>% 
	dplyr::select(starts_with("bedGraph"))

introgressed_variants_BainingKagat_SAFs_ucsc_custom_track <- introgressed_variants_PIB_SAFs_ucsc_custom_track %>% 
	mutate(
		bedGraph_chrom=seqnames,
		bedGraph_chromStart=format(start-1, scientific=FALSE),
		bedGraph_chromEnd=format(end, scientific=FALSE), 
		bedGraph_score=`SAF_BainingKagat`
	) %>% 
	dplyr::select(starts_with("bedGraph"))

introgressed_variants_BainingMali_SAFs_ucsc_custom_track <- introgressed_variants_PIB_SAFs_ucsc_custom_track %>% 
	mutate(
		bedGraph_chrom=seqnames,
		bedGraph_chromStart=format(start-1, scientific=FALSE),
		bedGraph_chromEnd=format(end, scientific=FALSE), 
		bedGraph_score=`SAF_BainingMali`
	) %>% 
	dplyr::select(starts_with("bedGraph"))

introgressed_variants_BellonaRennell_SAFs_ucsc_custom_track <- introgressed_variants_PIB_SAFs_ucsc_custom_track %>% 
	mutate(
		bedGraph_chrom=seqnames,
		bedGraph_chromStart=format(start-1, scientific=FALSE),
		bedGraph_chromEnd=format(end, scientific=FALSE), 
		bedGraph_score=`SAF_BellonaRennell`
	) %>% 
	dplyr::select(starts_with("bedGraph"))

introgressed_variants_GorokaSepik_SAFs_ucsc_custom_track <- introgressed_variants_PIB_SAFs_ucsc_custom_track %>% 
	mutate(
		bedGraph_chrom=seqnames,
		bedGraph_chromStart=format(start-1, scientific=FALSE),
		bedGraph_chromEnd=format(end, scientific=FALSE), 
		bedGraph_score=`SAF_GorokaSepik`
	) %>% 
	dplyr::select(starts_with("bedGraph"))

introgressed_variants_Kove_SAFs_ucsc_custom_track <- introgressed_variants_PIB_SAFs_ucsc_custom_track %>% 
	mutate(
		bedGraph_chrom=seqnames,
		bedGraph_chromStart=format(start-1, scientific=FALSE),
		bedGraph_chromEnd=format(end, scientific=FALSE), 
		bedGraph_score=`SAF_Kove`
	) %>% 
	dplyr::select(starts_with("bedGraph"))

introgressed_variants_LavongaiMussau_SAFs_ucsc_custom_track <- introgressed_variants_PIB_SAFs_ucsc_custom_track %>% 
	mutate(
		bedGraph_chrom=seqnames,
		bedGraph_chromStart=format(start-1, scientific=FALSE),
		bedGraph_chromEnd=format(end, scientific=FALSE), 
		bedGraph_score=`SAF_LavongaiMussau`
	) %>% 
	dplyr::select(starts_with("bedGraph"))

introgressed_variants_Malaita_SAFs_ucsc_custom_track <- introgressed_variants_PIB_SAFs_ucsc_custom_track %>% 
	mutate(
		bedGraph_chrom=seqnames,
		bedGraph_chromStart=format(start-1, scientific=FALSE),
		bedGraph_chromEnd=format(end, scientific=FALSE), 
		bedGraph_score=`SAF_Malaita`
	) %>% 
	dplyr::select(starts_with("bedGraph"))

introgressed_variants_Mamusi_SAFs_ucsc_custom_track <- introgressed_variants_PIB_SAFs_ucsc_custom_track %>% 
	mutate(
		bedGraph_chrom=seqnames,
		bedGraph_chromStart=format(start-1, scientific=FALSE),
		bedGraph_chromEnd=format(end, scientific=FALSE), 
		bedGraph_score=`SAF_Mamusi`
	) %>% 
	dplyr::select(starts_with("bedGraph"))

introgressed_variants_Melamela_SAFs_ucsc_custom_track <- introgressed_variants_PIB_SAFs_ucsc_custom_track %>% 
	mutate(
		bedGraph_chrom=seqnames,
		bedGraph_chromStart=format(start-1, scientific=FALSE),
		bedGraph_chromEnd=format(end, scientific=FALSE), 
		bedGraph_score=`SAF_Melamela`
	) %>% 
	dplyr::select(starts_with("bedGraph"))

introgressed_variants_NailikNotsiTigak_SAFs_ucsc_custom_track <- introgressed_variants_PIB_SAFs_ucsc_custom_track %>% 
	mutate(
		bedGraph_chrom=seqnames,
		bedGraph_chromStart=format(start-1, scientific=FALSE),
		bedGraph_chromEnd=format(end, scientific=FALSE), 
		bedGraph_score=`SAF_NailikNotsiTigak`
	) %>% 
	dplyr::select(starts_with("bedGraph"))

introgressed_variants_NakanaiMangseng_SAFs_ucsc_custom_track <- introgressed_variants_PIB_SAFs_ucsc_custom_track %>% 
	mutate(
		bedGraph_chrom=seqnames,
		bedGraph_chromStart=format(start-1, scientific=FALSE),
		bedGraph_chromEnd=format(end, scientific=FALSE), 
		bedGraph_score=`SAF_NakanaiMangseng`
	) %>% 
	dplyr::select(starts_with("bedGraph"))

introgressed_variants_Nasioi_SAFs_ucsc_custom_track <- introgressed_variants_PIB_SAFs_ucsc_custom_track %>% 
	mutate(
		bedGraph_chrom=seqnames,
		bedGraph_chromStart=format(start-1, scientific=FALSE),
		bedGraph_chromEnd=format(end, scientific=FALSE), 
		bedGraph_score=`SAF_Nasioi`
	) %>% 
	dplyr::select(starts_with("bedGraph"))

introgressed_variants_SantaCruz_SAFs_ucsc_custom_track <- introgressed_variants_PIB_SAFs_ucsc_custom_track %>% 
	mutate(
		bedGraph_chrom=seqnames,
		bedGraph_chromStart=format(start-1, scientific=FALSE),
		bedGraph_chromEnd=format(end, scientific=FALSE), 
		bedGraph_score=`SAF_SantaCruz`
	) %>% 
	dplyr::select(starts_with("bedGraph"))

introgressed_variants_Saposa_SAFs_ucsc_custom_track <- introgressed_variants_PIB_SAFs_ucsc_custom_track %>% 
	mutate(
		bedGraph_chrom=seqnames,
		bedGraph_chromStart=format(start-1, scientific=FALSE),
		bedGraph_chromEnd=format(end, scientific=FALSE), 
		bedGraph_score=`SAF_Saposa`
	) %>% 
	dplyr::select(starts_with("bedGraph"))

introgressed_variants_Tikopia_SAFs_ucsc_custom_track <- introgressed_variants_PIB_SAFs_ucsc_custom_track %>% 
	mutate(
		bedGraph_chrom=seqnames,
		bedGraph_chromStart=format(start-1, scientific=FALSE),
		bedGraph_chromEnd=format(end, scientific=FALSE), 
		bedGraph_score=`SAF_Tikopia`
	) %>% 
	dplyr::select(starts_with("bedGraph"))

introgressed_variants_VellaLavella_SAFs_ucsc_custom_track <- introgressed_variants_PIB_SAFs_ucsc_custom_track %>% 
	mutate(
		bedGraph_chrom=seqnames,
		bedGraph_chromStart=format(start-1, scientific=FALSE),
		bedGraph_chromEnd=format(end, scientific=FALSE), 
		bedGraph_score=`SAF_VellaLavella`
	) %>% 
	dplyr::select(starts_with("bedGraph"))

create_ucsc_custom_bedGraph_track(
	input_tb=introgressed_variants_Ata_SAFs_ucsc_custom_track, 
	output_bedGraph="../../results/1a-preprocess_PIBv1_MPRA_final/introgressed_variants_Ata_SAFs_ucsc_custom_track.bedGraph",
	name="Ata SAFs", description="Ata SAFs", , visibility="full", alwaysZero="on", viewLimits=0:1, color=paste(as.vector(col2rgb("#000000")), collapse=","))
create_ucsc_custom_bedGraph_track(
	input_tb=introgressed_variants_BainingKagat_SAFs_ucsc_custom_track, 
	output_bedGraph="../../results/1a-preprocess_PIBv1_MPRA_final/introgressed_variants_BainingKagat_SAFs_ucsc_custom_track.bedGraph",
	name="BainingKagat SAFs", description="BainingKagat SAFs", , visibility="full", alwaysZero="on", viewLimits=0:1, color=paste(as.vector(col2rgb("#000000")), collapse=","))
create_ucsc_custom_bedGraph_track(
	input_tb=introgressed_variants_BainingMali_SAFs_ucsc_custom_track, 
	output_bedGraph="../../results/1a-preprocess_PIBv1_MPRA_final/introgressed_variants_BainingMali_SAFs_ucsc_custom_track.bedGraph",
	name="BainingMali SAFs", description="BainingMali SAFs", , visibility="full", alwaysZero="on", viewLimits=0:1, color=paste(as.vector(col2rgb("#000000")), collapse=","))
create_ucsc_custom_bedGraph_track(
	input_tb=introgressed_variants_BellonaRennell_SAFs_ucsc_custom_track, 
	output_bedGraph="../../results/1a-preprocess_PIBv1_MPRA_final/introgressed_variants_BellonaRennell_SAFs_ucsc_custom_track.bedGraph",
	name="BellonaRennell SAFs", description="BellonaRennell SAFs", , visibility="full", alwaysZero="on", viewLimits=0:1, color=paste(as.vector(col2rgb("#000000")), collapse=","))
create_ucsc_custom_bedGraph_track(
	input_tb=introgressed_variants_GorokaSepik_SAFs_ucsc_custom_track, 
	output_bedGraph="../../results/1a-preprocess_PIBv1_MPRA_final/introgressed_variants_GorokaSepik_SAFs_ucsc_custom_track.bedGraph",
	name="GorokaSepik SAFs", description="GorokaSepik SAFs", , visibility="full", alwaysZero="on", viewLimits=0:1, color=paste(as.vector(col2rgb("#000000")), collapse=","))
create_ucsc_custom_bedGraph_track(
	input_tb=introgressed_variants_Kove_SAFs_ucsc_custom_track, 
	output_bedGraph="../../results/1a-preprocess_PIBv1_MPRA_final/introgressed_variants_Kove_SAFs_ucsc_custom_track.bedGraph",
	name="Kove SAFs", description="Kove SAFs", , visibility="full", alwaysZero="on", viewLimits=0:1, color=paste(as.vector(col2rgb("#000000")), collapse=","))
create_ucsc_custom_bedGraph_track(
	input_tb=introgressed_variants_LavongaiMussau_SAFs_ucsc_custom_track, 
	output_bedGraph="../../results/1a-preprocess_PIBv1_MPRA_final/introgressed_variants_LavongaiMussau_SAFs_ucsc_custom_track.bedGraph",
	name="LavongaiMussau SAFs", description="LavongaiMussau SAFs", , visibility="full", alwaysZero="on", viewLimits=0:1, color=paste(as.vector(col2rgb("#000000")), collapse=","))
create_ucsc_custom_bedGraph_track(
	input_tb=introgressed_variants_Malaita_SAFs_ucsc_custom_track, 
	output_bedGraph="../../results/1a-preprocess_PIBv1_MPRA_final/introgressed_variants_Malaita_SAFs_ucsc_custom_track.bedGraph",
	name="Malaita SAFs", description="Malaita SAFs", , visibility="full", alwaysZero="on", viewLimits=0:1, color=paste(as.vector(col2rgb("#000000")), collapse=","))
create_ucsc_custom_bedGraph_track(
	input_tb=introgressed_variants_Mamusi_SAFs_ucsc_custom_track, 
	output_bedGraph="../../results/1a-preprocess_PIBv1_MPRA_final/introgressed_variants_Mamusi_SAFs_ucsc_custom_track.bedGraph",
	name="Mamusi SAFs", description="Mamusi SAFs", , visibility="full", alwaysZero="on", viewLimits=0:1, color=paste(as.vector(col2rgb("#000000")), collapse=","))
create_ucsc_custom_bedGraph_track(
	input_tb=introgressed_variants_Melamela_SAFs_ucsc_custom_track, 
	output_bedGraph="../../results/1a-preprocess_PIBv1_MPRA_final/introgressed_variants_Melamela_SAFs_ucsc_custom_track.bedGraph",
	name="Melamela SAFs", description="Melamela SAFs", , visibility="full", alwaysZero="on", viewLimits=0:1, color=paste(as.vector(col2rgb("#000000")), collapse=","))
create_ucsc_custom_bedGraph_track(
	input_tb=introgressed_variants_NailikNotsiTigak_SAFs_ucsc_custom_track, 
	output_bedGraph="../../results/1a-preprocess_PIBv1_MPRA_final/introgressed_variants_NailikNotsiTigak_SAFs_ucsc_custom_track.bedGraph",
	name="NailikNotsiTigak SAFs", description="NailikNotsiTigak SAFs", , visibility="full", alwaysZero="on", viewLimits=0:1, color=paste(as.vector(col2rgb("#000000")), collapse=","))
create_ucsc_custom_bedGraph_track(
	input_tb=introgressed_variants_NakanaiMangseng_SAFs_ucsc_custom_track, 
	output_bedGraph="../../results/1a-preprocess_PIBv1_MPRA_final/introgressed_variants_NakanaiMangseng_SAFs_ucsc_custom_track.bedGraph",
	name="NakanaiMangseng SAFs", description="NakanaiMangseng SAFs", , visibility="full", alwaysZero="on", viewLimits=0:1, color=paste(as.vector(col2rgb("#000000")), collapse=","))
create_ucsc_custom_bedGraph_track(
	input_tb=introgressed_variants_Nasioi_SAFs_ucsc_custom_track, 
	output_bedGraph="../../results/1a-preprocess_PIBv1_MPRA_final/introgressed_variants_Nasioi_SAFs_ucsc_custom_track.bedGraph",
	name="Nasioi SAFs", description="Nasioi SAFs", , visibility="full", alwaysZero="on", viewLimits=0:1, color=paste(as.vector(col2rgb("#000000")), collapse=","))
create_ucsc_custom_bedGraph_track(
	input_tb=introgressed_variants_SantaCruz_SAFs_ucsc_custom_track, 
	output_bedGraph="../../results/1a-preprocess_PIBv1_MPRA_final/introgressed_variants_SantaCruz_SAFs_ucsc_custom_track.bedGraph",
	name="SantaCruz SAFs", description="SantaCruz SAFs", , visibility="full", alwaysZero="on", viewLimits=0:1, color=paste(as.vector(col2rgb("#000000")), collapse=","))
create_ucsc_custom_bedGraph_track(
	input_tb=introgressed_variants_Saposa_SAFs_ucsc_custom_track, 
	output_bedGraph="../../results/1a-preprocess_PIBv1_MPRA_final/introgressed_variants_Saposa_SAFs_ucsc_custom_track.bedGraph",
	name="Saposa SAFs", description="Saposa SAFs", , visibility="full", alwaysZero="on", viewLimits=0:1, color=paste(as.vector(col2rgb("#000000")), collapse=","))
create_ucsc_custom_bedGraph_track(
	input_tb=introgressed_variants_Tikopia_SAFs_ucsc_custom_track, 
	output_bedGraph="../../results/1a-preprocess_PIBv1_MPRA_final/introgressed_variants_Tikopia_SAFs_ucsc_custom_track.bedGraph",
	name="Tikopia SAFs", description="Tikopia SAFs", , visibility="full", alwaysZero="on", viewLimits=0:1, color=paste(as.vector(col2rgb("#000000")), collapse=","))
create_ucsc_custom_bedGraph_track(
	input_tb=introgressed_variants_VellaLavella_SAFs_ucsc_custom_track, 
	output_bedGraph="../../results/1a-preprocess_PIBv1_MPRA_final/introgressed_variants_VellaLavella_SAFs_ucsc_custom_track.bedGraph",
	name="VellaLavella SAFs", description="VellaLavella SAFs", , visibility="full", alwaysZero="on", viewLimits=0:1, color=paste(as.vector(col2rgb("#000000")), collapse=","))
