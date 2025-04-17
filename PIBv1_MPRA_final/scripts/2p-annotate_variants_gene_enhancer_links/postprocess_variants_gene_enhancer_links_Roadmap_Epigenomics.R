#!/bin/R

# for data analysis
library(tidyverse)
library(data.table)

# for bio data analysis
library(plyranges)
library(rtracklayer)
source("../shared_functions/seqinfo_fix_change.R")

# load variants
introgressed_variants_tb <- as_tibble(readRDS("../../results/1a-preprocess_PIBv1_MPRA_final/introgressed_variants.rds"))[,1:10]

# load Roadmap_Epigenomics
Roadmap_Epigenomics_gene_enhancer_links_tb <- as_tibble(fread(paste0("../../../Datasets/gene_regulation_gene_enhancer_links/Roadmap_Epigenomics_gene_enhancer_links/data_cleanup/GRCh37_merged/Roadmap_Epigenomics_gene_enhancer_links.txt.gz")))

# remove ambig TSS
Roadmap_Epigenomics_gene_enhancer_links_tb <- Roadmap_Epigenomics_gene_enhancer_links_tb %>% 
	filter(!is.na(target_symbol)) %>% 
	filter(!is.na(target_tss))

# Roadmap_Epigenomics names
names(Roadmap_Epigenomics_gene_enhancer_links_tb)[6:length(Roadmap_Epigenomics_gene_enhancer_links_tb)] <- paste("Roadmap_Epigenomics_orig", names(Roadmap_Epigenomics_gene_enhancer_links_tb)[6:length(Roadmap_Epigenomics_gene_enhancer_links_tb)], sep="_")
introgressed_variants_Roadmap_Epigenomics_names <- names(Roadmap_Epigenomics_gene_enhancer_links_tb)[6:length(Roadmap_Epigenomics_gene_enhancer_links_tb)]

# convert to GRange
Roadmap_Epigenomics_gene_enhancer_links_tb <- Roadmap_Epigenomics_gene_enhancer_links_tb %>% 
	mutate(seqnames = gsub("chr", "", seqnames))

# overlap variants and Roadmap_Epigenomics
introgressed_variants_Roadmap_Epigenomics <- as_tibble(find_overlaps(GRanges(introgressed_variants_tb), GRanges(Roadmap_Epigenomics_gene_enhancer_links_tb)))

# get max per Roadmap_Epigenomics group
introgressed_variants_Roadmap_Epigenomics <- introgressed_variants_Roadmap_Epigenomics %>% 
	group_by(seqnames, start, end, width, strand, VariantID, VariantCHROM, VariantPOS, VariantREF, VariantALT, Roadmap_Epigenomics_orig_group) %>% 
	arrange(Roadmap_Epigenomics_orig_score) %>% 
	slice(n()) %>% 
	ungroup()

# Roadmap_Epigenomics names
introgressed_variants_Roadmap_Epigenomics_names <- setdiff(names(introgressed_variants_Roadmap_Epigenomics)[11:length(introgressed_variants_Roadmap_Epigenomics)], "Roadmap_Epigenomics_orig_group")

# add summary cols
introgressed_variants_Roadmap_Epigenomics <- introgressed_variants_Roadmap_Epigenomics %>% 
	rowwise() %>% 
	mutate(Roadmap_Epigenomics_summ_Score = Roadmap_Epigenomics_orig_score) %>% 
	mutate(Roadmap_Epigenomics_summ_Ensembl_ID = Roadmap_Epigenomics_orig_target_ensembl_id) %>% 
	mutate(Roadmap_Epigenomics_summ_Symbol = Roadmap_Epigenomics_orig_target_symbol) %>% 
	mutate(Roadmap_Epigenomics_summ_Biosample = Roadmap_Epigenomics_orig_epigenome) %>% 
	mutate(Roadmap_Epigenomics_summ_Category = Roadmap_Epigenomics_orig_group)

# save for interact plot
# 	because want not collapsed
introgressed_variants_Roadmap_Epigenomics_interact <- introgressed_variants_Roadmap_Epigenomics
write_tsv(as_tibble(introgressed_variants_Roadmap_Epigenomics_interact), gzfile("../../results/2p-annotate_variants_gene_enhancer_links/introgressed_variants_Roadmap_Epigenomics_interact.txt.gz"))

# Roadmap_Epigenomics names
introgressed_variants_Roadmap_Epigenomics_names <- setdiff(names(introgressed_variants_Roadmap_Epigenomics)[11:length(introgressed_variants_Roadmap_Epigenomics)], "Roadmap_Epigenomics_orig_group")

# collapse by group
introgressed_variants_Roadmap_Epigenomics_collapse <- introgressed_variants_Roadmap_Epigenomics %>% 
	group_by(seqnames, start, end, width, strand, VariantID, VariantCHROM, VariantPOS, VariantREF, VariantALT) %>% 
	dplyr::summarise_at(introgressed_variants_Roadmap_Epigenomics_names, function(x) {paste(x, collapse=",")}) %>% 
	ungroup() %>% 
	rowwise() %>% 
	mutate(Roadmap_Epigenomics_summ_Count = sum(!is.na(strsplit(Roadmap_Epigenomics_orig_target_symbol, ",")[[1]]))) %>% 
	ungroup() %>% 
	mutate(Roadmap_Epigenomics_summ_Summary = "Roadmap_Epigenomics_gene_enhancer_linked") %>% 
	dplyr::select(seqnames, start, end, width, strand, VariantID, VariantCHROM, VariantPOS, VariantREF, VariantALT, 
		starts_with("Roadmap_Epigenomics_orig"), Roadmap_Epigenomics_summ_Summary, Roadmap_Epigenomics_summ_Count, starts_with("Roadmap_Epigenomics_summ"))

# pivot by group
introgressed_variants_Roadmap_Epigenomics_pivot <- introgressed_variants_Roadmap_Epigenomics %>% 
	pivot_wider(names_from=Roadmap_Epigenomics_orig_group, names_sep="-", values_from=introgressed_variants_Roadmap_Epigenomics_names)

# join together
introgressed_variants_Roadmap_Epigenomics <- full_join(introgressed_variants_Roadmap_Epigenomics_collapse, introgressed_variants_Roadmap_Epigenomics_pivot)

# merge to full variants
introgressed_variants_tb <- as_tibble(fread("../../results/1a-preprocess_PIBv1_MPRA_final/introgressed_variants.txt.gz")) %>% 
	dplyr::select(seqnames, start, end, width, strand, VariantID, VariantCHROM, VariantPOS, VariantREF, VariantALT)
names_temp <- intersect(names(introgressed_variants_tb), names(introgressed_variants_Roadmap_Epigenomics))
introgressed_variants_Roadmap_Epigenomics <- introgressed_variants_tb %>% mutate_at(names_temp, as.character)%>% 
	left_join(introgressed_variants_Roadmap_Epigenomics %>% mutate_at(names_temp, as.character))

# save all variants
write_tsv(introgressed_variants_Roadmap_Epigenomics, gzfile("../../results/2p-annotate_variants_gene_enhancer_links/introgressed_variants_Roadmap_Epigenomics.txt.gz"))

# save summ variants
introgressed_variants_Roadmap_Epigenomics_summ <- introgressed_variants_Roadmap_Epigenomics %>% 
	dplyr::select(seqnames, start, end, width, strand, starts_with("Variant"), contains("_summ_"))
write_tsv(introgressed_variants_Roadmap_Epigenomics_summ, gzfile("../../results/2p-annotate_variants_gene_enhancer_links/introgressed_variants_Roadmap_Epigenomics_summ.txt.gz"))
sort(table(introgressed_variants_Roadmap_Epigenomics_summ$Roadmap_Epigenomics_summ_Summary))
