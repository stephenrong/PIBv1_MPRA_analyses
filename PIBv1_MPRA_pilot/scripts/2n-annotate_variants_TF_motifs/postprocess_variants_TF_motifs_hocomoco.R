#!/bin/R

# for data analysis
library(tidyverse)
library(data.table)

# for bio data analysis
library(plyranges)

# load TF motif
adaptive_variants_TF_motifs_hocomoco <- as_tibble(readRDS("../../results/2n-annotate_variants_TF_motifs/adaptive_variants_TF_motifs_hocomoco.rds"))

# TF motif names
names(adaptive_variants_TF_motifs_hocomoco)[9:length(adaptive_variants_TF_motifs_hocomoco)] <- paste("TF_motifs_hocomoco_orig", names(adaptive_variants_TF_motifs_hocomoco)[9:length(adaptive_variants_TF_motifs_hocomoco)], sep="_")
adaptive_variants_TF_motifs_hocomoco_names <- names(adaptive_variants_TF_motifs_hocomoco)[9:length(adaptive_variants_TF_motifs_hocomoco)]

# add GRange and Variant cols
adaptive_variants_TF_motifs_hocomoco <- adaptive_variants_TF_motifs_hocomoco %>% 
	mutate(seqnames = gsub("chr", "", seqnames)) %>% 
	mutate(TF_motifs_hocomoco_orig_strand = strand) %>% 
	mutate(strand = "*") %>% 
	mutate(VariantCHROM = seqnames, VariantPOS = start, VariantREF = REF, VariantALT = ALT) %>% 
	mutate(VariantID = paste(VariantCHROM, VariantPOS, VariantREF, VariantALT, sep="_")) %>% 
	dplyr::select(seqnames, start, end, width, strand, VariantID, VariantCHROM, VariantPOS, VariantREF, VariantALT, TF_motifs_hocomoco_orig_strand, starts_with("TF_motifs_hocomoco_orig"))

# collapse multiple entries
adaptive_variants_TF_motifs_hocomoco <- adaptive_variants_TF_motifs_hocomoco %>% 
	group_by(seqnames, start, end, width, strand, VariantID, VariantCHROM, VariantPOS, VariantREF, VariantALT) %>% 
	arrange(desc(abs(TF_motifs_hocomoco_orig_alleleDiff))) %>% 
	dplyr::summarise_at(adaptive_variants_TF_motifs_hocomoco_names, function(x) {paste(x, collapse=",")}) %>% 
	ungroup()

# add summary cols
adaptive_variants_TF_motifs_hocomoco <- adaptive_variants_TF_motifs_hocomoco %>% 
	rowwise() %>% mutate(TF_motifs_hocomoco_summ_Count = sum(!is.na(strsplit(TF_motifs_hocomoco_orig_geneSymbol, ",")[[1]]))) %>% ungroup() %>% 
	mutate(TF_motifs_hocomoco_summ_Description = TF_motifs_hocomoco_orig_geneSymbol) %>% 
	mutate(TF_motifs_hocomoco_summ_Summary = ifelse(TF_motifs_hocomoco_summ_Count > 0, "TF_motif_disrupting", NA))

# merge to full variants
adaptive_variants_tb <- as_tibble(fread("../../results/1a-preprocess_PIBv1_MPRA_pilot/adaptive_variants.txt.gz")) %>% 
	dplyr::select(seqnames, start, end, width, strand, VariantID, VariantCHROM, VariantPOS, VariantREF, VariantALT)
names_temp <- intersect(names(adaptive_variants_tb), names(adaptive_variants_TF_motifs_hocomoco))
adaptive_variants_TF_motifs_hocomoco <- adaptive_variants_tb %>% mutate_at(names_temp, as.character)%>% 
	left_join(adaptive_variants_TF_motifs_hocomoco %>% mutate_at(names_temp, as.character))

# save all variants
write_tsv(adaptive_variants_TF_motifs_hocomoco, gzfile("../../results/2n-annotate_variants_TF_motifs/adaptive_variants_TF_motifs_hocomoco.txt.gz"))

# save summ variants
adaptive_variants_TF_motifs_hocomoco_summ <- adaptive_variants_TF_motifs_hocomoco %>% 
	dplyr::select(seqnames, start, end, width, strand, starts_with("Variant"), contains("_summ_"))
write_tsv(adaptive_variants_TF_motifs_hocomoco_summ, gzfile("../../results/2n-annotate_variants_TF_motifs/adaptive_variants_TF_motifs_hocomoco_summ.txt.gz"))
sort(table(adaptive_variants_TF_motifs_hocomoco_summ$TF_motifs_hocomoco_summ_Summary))

# visual check
pdf("../../results/2n-annotate_variants_TF_motifs/adaptive_variants_TF_motifs_hocomoco_summ_pie.pdf")
pie(rev(sort(table(adaptive_variants_TF_motifs_hocomoco_summ$TF_motifs_hocomoco_summ_Count))))
dev.off()
