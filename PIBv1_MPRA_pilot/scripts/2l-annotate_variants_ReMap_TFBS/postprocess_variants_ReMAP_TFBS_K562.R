#!/bin/R

# for data analysis
library(tidyverse)
library(data.table)

# for bio data analysis
library(plyranges)
library(rtracklayer)

# load VCF 
adaptive_variants_ReMap_TFBS_K562 <- as_tibble(fread("../../results/1a-preprocess_PIBv1_MPRA_pilot/adaptive_variants.txt.gz")) %>% 
	dplyr::select(seqnames, start, end, width, strand, VariantID, VariantCHROM, VariantPOS, VariantREF, VariantALT)

# load BED
ReMap_TFBS_K562_tb <- as_tibble(fread("../../../../Datasets/gene_regulation_binding_catalogs/ReMap2022_ChIP_seq_homo_sapiens/data_cleanup/nonredundant_peaks/remap2022_nr_macs2_hg19_v1_0.K-562.bed.gz"))
ReMap_TFBS_K562_tb <- ReMap_TFBS_K562_tb[,c(1,2,3,4)]
ReMap_TFBS_K562_tb <- unique(ReMap_TFBS_K562_tb)
names(ReMap_TFBS_K562_tb) <- c("seqnames", "start", "end", "ReMap_TFBS_K562_orig_TFs")
ReMap_TFBS_K562_tb <- ReMap_TFBS_K562_tb %>% 
	mutate(seqnames = gsub("chr", "", seqnames)) %>% 
	mutate(start = start+1) %>% 
	mutate(ReMap_TFBS_K562_orig_TFs = gsub(":.*", "", ReMap_TFBS_K562_orig_TFs)) %>% 
	mutate(ReMap_TFBS_K562_orig_TFs_element = paste(paste(seqnames, start, end, sep="_"), ReMap_TFBS_K562_orig_TFs, sep="-"))

# overlap
overlap_temp <- as_tibble(findOverlaps(GRanges(adaptive_variants_ReMap_TFBS_K562), GRanges(ReMap_TFBS_K562_tb)))
overlap_temp$VariantID <- adaptive_variants_ReMap_TFBS_K562$VariantID[overlap_temp$queryHits]
overlap_temp$ReMap_TFBS_K562_orig_TFs <- ReMap_TFBS_K562_tb$ReMap_TFBS_K562_orig_TFs[overlap_temp$subjectHits]
overlap_temp$ReMap_TFBS_K562_orig_TFs_element <- ReMap_TFBS_K562_tb$ReMap_TFBS_K562_orig_TFs_element[overlap_temp$subjectHits]
overlap_temp <- overlap_temp %>% 
	dplyr::select(-queryHits, -subjectHits)
adaptive_variants_ReMap_TFBS_K562 <- adaptive_variants_ReMap_TFBS_K562 %>% 
	left_join(overlap_temp)

# save names
adaptive_variants_ReMap_TFBS_K562_names <- names(adaptive_variants_ReMap_TFBS_K562)[11:length(adaptive_variants_ReMap_TFBS_K562)]

# collapse multiple entries
adaptive_variants_ReMap_TFBS_K562 <- adaptive_variants_ReMap_TFBS_K562 %>% 
	group_by(seqnames, start, end, width, strand, VariantID, VariantCHROM, VariantPOS, VariantREF, VariantALT) %>% 
	dplyr::summarise_at(adaptive_variants_ReMap_TFBS_K562_names, function(x) {paste(x, collapse=",")}) %>% 
	mutate_at(adaptive_variants_ReMap_TFBS_K562_names, function(x) {ifelse(x=="NA", NA, x)}) %>% 
	ungroup()

# add summary cols
adaptive_variants_ReMap_TFBS_K562 <- adaptive_variants_ReMap_TFBS_K562 %>% 
	rowwise() %>% mutate(ReMap_TFBS_K562_summ_Count = sum(!is.na(strsplit(ReMap_TFBS_K562_orig_TFs_element, ",")[[1]]))) %>% ungroup() %>% 
	mutate(ReMap_TFBS_K562_summ_Description = ReMap_TFBS_K562_orig_TFs) %>% 
	mutate(ReMap_TFBS_K562_summ_Summary = ifelse(ReMap_TFBS_K562_summ_Count > 0, "TF_ReMap_ChIP_peak", NA))

# merge to full variants
adaptive_variants_tb <- as_tibble(fread("../../results/1a-preprocess_PIBv1_MPRA_pilot/adaptive_variants.txt.gz")) %>% 
	dplyr::select(seqnames, start, end, width, strand, VariantID, VariantCHROM, VariantPOS, VariantREF, VariantALT)
names_inter <- intersect(names(adaptive_variants_tb), names(adaptive_variants_ReMap_TFBS_K562))
adaptive_variants_ReMap_TFBS_K562 <- adaptive_variants_tb %>% mutate_at(names_inter, as.character) %>% 
	left_join(adaptive_variants_ReMap_TFBS_K562 %>% mutate_at(names_inter, as.character))

# save all variants
write_tsv(adaptive_variants_ReMap_TFBS_K562, gzfile("../../results/2l-annotate_variants_ReMap_TFBS/adaptive_variants_ReMap_TFBS_K562.txt.gz"))

# save summ variants
adaptive_variants_ReMap_TFBS_K562_summ <- adaptive_variants_ReMap_TFBS_K562 %>% 
	dplyr::select(seqnames, start, end, width, strand, starts_with("Variant"), contains("_summ_"))
write_tsv(adaptive_variants_ReMap_TFBS_K562_summ, gzfile("../../results/2l-annotate_variants_ReMap_TFBS/adaptive_variants_ReMap_TFBS_K562_summ.txt.gz"))
sort(table(adaptive_variants_ReMap_TFBS_K562_summ$ReMap_TFBS_K562_summ_Summary))
