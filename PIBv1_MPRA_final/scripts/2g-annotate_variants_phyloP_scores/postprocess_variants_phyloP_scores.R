#!/bin/R

# for data analysis
library(tidyverse)
library(data.table)

# for bio data analysis
library(plyranges)
library(rtracklayer)
library(vcfR)
source("../shared_functions/seqinfo_fix_change.R")

# load BED
introgressed_variants_phyloP <- import("../../results/2g-annotate_variants_phyloP_scores/introgressed_variants_lift37_final-241-mammalian-2020v2.bed") %>% 
	seqinfo_fix("NCBI", "GRCh37") %>% 
	as_tibble()

# add GRange and Variant cols
introgressed_variants_phyloP <- introgressed_variants_phyloP %>% 
	mutate(VariantCHROM = seqnames, VariantPOS = start) %>% 
	dplyr::select(seqnames, start, end, width, strand, VariantCHROM, VariantPOS, score) 

# add summary cols
introgressed_variants_phyloP <- introgressed_variants_phyloP %>%  
	dplyr::rename(phyloP_mam241_summ_Score = score) %>% 
	mutate(phyloP_mam241_summ_Summary = ifelse(phyloP_mam241_summ_Score >= 2.27, "constrained", NA))

# merge to full variants
introgressed_variants_tb <- as_tibble(fread("../../results/1a-preprocess_PIBv1_MPRA_final/introgressed_variants.txt.gz")) %>% 
	dplyr::select(seqnames, start, end, width, strand, VariantID, VariantCHROM, VariantPOS, VariantREF, VariantALT)
names_temp <- intersect(names(introgressed_variants_tb), names(introgressed_variants_phyloP))
introgressed_variants_phyloP <- introgressed_variants_tb %>% mutate_at(names_temp, as.character) %>% 
	left_join(introgressed_variants_phyloP %>% mutate_at(names_temp, as.character))

# save all variants
write_tsv(introgressed_variants_phyloP, gzfile("../../results/2g-annotate_variants_phyloP_scores/introgressed_variants_phyloP.txt.gz"))

# save summ variants
introgressed_variants_phyloP_summ <- introgressed_variants_phyloP %>% 
	dplyr::select(seqnames, start, end, width, strand, starts_with("Variant"), contains("_summ_"))
write_tsv(introgressed_variants_phyloP_summ, gzfile("../../results/2g-annotate_variants_phyloP_scores/introgressed_variants_phyloP_summ.txt.gz"))
sort(table(introgressed_variants_phyloP_summ$phyloP_mam241_summ_Summary))

# visual check
pdf("../../results/2g-annotate_variants_phyloP_scores/introgressed_variants_phyloP_summ_sort.pdf")
plot(sort(introgressed_variants_phyloP$phyloP_mam241_summ_Score), ylab="phyloP mam241 score")
dev.off()

pdf("../../results/2g-annotate_variants_phyloP_scores/introgressed_variants_phyloP_summ_pie.pdf")
pie(rev(sort(table(introgressed_variants_phyloP_summ$phyloP_mam241_summ_Summary))))
dev.off()
