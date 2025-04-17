#!/bin/R

# for data analysis
library(tidyverse)
library(data.table)

# for bio data analysis
library(plyranges)
library(rtracklayer)
library(genekitr)
source("../shared_functions/seqinfo_fix_change.R")

# load variants
introgressed_variants_nearest_gene <- as_tibble(readRDS("../../results/1a-preprocess_PIBv1_MPRA_final/introgressed_variants.rds"))[,1:10]

# get gene bodies
ucsc_knownGene_abbrev <- readRDS("../../../Datasets/gene_annotations/get_canonical_genes/data_cleanup/canonical_gene_lift37_gr.rds") %>% 
	mutate(nearest_gene_orig_ensembl_id = tx_ensembl_gene_id, nearest_gene_orig_symbol = tx_hgnc_symbol) %>% 
	filter(!is.na(tx_hgnc_symbol))

# get nearest
overlap_temp <- as_tibble(distanceToNearest(GRanges(introgressed_variants_nearest_gene), GRanges(ucsc_knownGene_abbrev), select="all"))
overlap_temp$VariantID <- introgressed_variants_nearest_gene$VariantID[overlap_temp$queryHits]
overlap_temp$nearest_gene_orig_ensembl_id <- ucsc_knownGene_abbrev$nearest_gene_orig_ensembl_id[overlap_temp$subjectHits]
overlap_temp$nearest_gene_orig_symbol <- ucsc_knownGene_abbrev$nearest_gene_orig_symbol[overlap_temp$subjectHits]
overlap_temp$nearest_gene_orig_distance <- as.integer(overlap_temp$distance)
overlap_temp <- overlap_temp %>% 
	filter(distance <= 5e5) %>% 
	dplyr::select(-queryHits, -subjectHits, -distance)

introgressed_variants_nearest_gene <- introgressed_variants_nearest_gene %>% 
	left_join(overlap_temp)

# save names
introgressed_variants_nearest_gene_names <- names(introgressed_variants_nearest_gene)[11:length(introgressed_variants_nearest_gene)]

# collapse multiple entries
introgressed_variants_nearest_gene <- introgressed_variants_nearest_gene %>% 
	group_by(seqnames, start, end, width, strand, VariantID, VariantCHROM, VariantPOS, VariantREF, VariantALT) %>% 
	dplyr::summarise_at(introgressed_variants_nearest_gene_names, function(x) {paste(x, collapse=",")}) %>% 
	mutate_at(introgressed_variants_nearest_gene_names, function(x) {ifelse(x=="NA", NA, x)}) %>% 
	ungroup()

# add summary cols
introgressed_variants_nearest_gene <- introgressed_variants_nearest_gene %>% 
	rowwise() %>% mutate(nearest_gene_summ_Count = sum(!is.na(strsplit(nearest_gene_orig_symbol, ",")[[1]]))) %>% ungroup() %>% 
	mutate(nearest_gene_summ_Symbol = nearest_gene_orig_symbol)

# merge to full variants
introgressed_variants_tb <- as_tibble(fread("../../results/1a-preprocess_PIBv1_MPRA_final/introgressed_variants.txt.gz")) %>% 
	dplyr::select(seqnames, start, end, width, strand, VariantID, VariantCHROM, VariantPOS, VariantREF, VariantALT)
names_inter <- intersect(names(introgressed_variants_tb), names(introgressed_variants_nearest_gene))
introgressed_variants_nearest_gene <- introgressed_variants_tb %>% mutate_at(names_inter, as.character) %>% 
	left_join(introgressed_variants_nearest_gene %>% mutate_at(names_inter, as.character))

# save all variants
write_tsv(introgressed_variants_nearest_gene, gzfile("../../results/2o-annotate_variants_nearest_gene/introgressed_variants_nearest_gene.txt.gz"))

# save summ variants
introgressed_variants_nearest_gene_summ <- introgressed_variants_nearest_gene %>% 
	dplyr::select(seqnames, start, end, width, strand, starts_with("Variant"), contains("_summ_"))
write_tsv(introgressed_variants_nearest_gene_summ, gzfile("../../results/2o-annotate_variants_nearest_gene/introgressed_variants_nearest_gene_summ.txt.gz"))
sort(table(introgressed_variants_nearest_gene_summ$nearest_gene_summ_Symbol))
