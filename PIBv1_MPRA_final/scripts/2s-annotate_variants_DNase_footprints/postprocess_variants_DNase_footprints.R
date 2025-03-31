#!/bin/R

# for data analysis
library(tidyverse)
library(data.table)

# for bio data analysis
library(plyranges)
library(rtracklayer)

# load VCF 
introgressed_variants_DNase_footprints <- as_tibble(fread("../../results/1a-preprocess_PIBv1_MPRA_final/introgressed_variants.txt.gz")) %>% 
	dplyr::select(seqnames, start, end, width, strand, VariantID, VariantCHROM, VariantPOS, VariantREF, VariantALT)

# load BED
DNase_footprints_tb <- as_tibble(fread("../../../../Datasets/gene_regulation_binding_catalogs/DNase_digital_genomic_footprinting_consensus/data_cleanup/GRCh37_merged/consensus_footprints_and_motifs_hg38_lift37_GRanges.txt.gz"))
names(DNase_footprints_tb) <- c("seqnames", "start", "end", "width", "strand", "identifier", "mean_signal", "num_samples", "num_fps", "summit", "core_start", "core_end", "motif_clusters")

names(DNase_footprints_tb)[6:13] <- paste("DNase_footprints_orig", names(DNase_footprints_tb)[6:13], sep="_")
DNase_footprints_tb <- DNase_footprints_tb %>% 
	mutate(seqnames = gsub("chr", "", seqnames)) %>% 
	mutate(DNase_footprints_orig_motifs = DNase_footprints_orig_motif_clusters) %>% 
	mutate(DNase_footprints_orig_element = paste(paste(seqnames, start, end, sep="_"), "DNase_footprint", sep="-"))
DNase_footprints_gr <- GRanges(DNase_footprints_tb)

# overlap
overlap_temp <- as_tibble(findOverlaps(GRanges(introgressed_variants_DNase_footprints), DNase_footprints_gr))
overlap_temp$VariantID <- introgressed_variants_DNase_footprints$VariantID[overlap_temp$queryHits]
overlap_temp$DNase_footprints_orig_motifs <- DNase_footprints_gr$DNase_footprints_orig_motifs[overlap_temp$subjectHits]
overlap_temp$DNase_footprints_orig_element <- DNase_footprints_gr$DNase_footprints_orig_element[overlap_temp$subjectHits]
overlap_temp <- overlap_temp %>% 
	dplyr::select(-queryHits, -subjectHits)

introgressed_variants_DNase_footprints <- introgressed_variants_DNase_footprints %>% 
	left_join(overlap_temp)

# save names
introgressed_variants_DNase_footprints_names <- names(introgressed_variants_DNase_footprints)[11:length(introgressed_variants_DNase_footprints)]

# collapse multiple entries
introgressed_variants_DNase_footprints <- introgressed_variants_DNase_footprints %>% 
	group_by(seqnames, start, end, width, strand, VariantID, VariantCHROM, VariantPOS, VariantREF, VariantALT) %>% 
	dplyr::summarise_at(introgressed_variants_DNase_footprints_names, function(x) {paste(x, collapse=",")}) %>% 
	mutate_at(introgressed_variants_DNase_footprints_names, function(x) {ifelse(x=="NA", NA, x)}) %>% 
	ungroup()

# add summary cols
introgressed_variants_DNase_footprints <- introgressed_variants_DNase_footprints %>% 
	rowwise() %>% mutate(DNase_footprints_summ_Count = sum(!is.na(strsplit(DNase_footprints_orig_element, ",")[[1]]))) %>% ungroup() %>% 
	mutate(DNase_footprints_summ_Description = DNase_footprints_orig_motifs) %>% 
	mutate(DNase_footprints_summ_Summary = ifelse(DNase_footprints_summ_Count > 0, "DNase_TF_footprint", NA))

# merge to full variants
introgressed_variants_tb <- as_tibble(fread("../../results/1a-preprocess_PIBv1_MPRA_final/introgressed_variants.txt.gz")) %>% 
	dplyr::select(seqnames, start, end, width, strand, VariantID, VariantCHROM, VariantPOS, VariantREF, VariantALT)
names_inter <- intersect(names(introgressed_variants_tb), names(introgressed_variants_DNase_footprints))
introgressed_variants_DNase_footprints <- introgressed_variants_tb %>% mutate_at(names_inter, as.character) %>% 
	left_join(introgressed_variants_DNase_footprints %>% mutate_at(names_inter, as.character))

# save all variants
write_tsv(introgressed_variants_DNase_footprints, gzfile("../../results/2s-annotate_variants_DNase_footprints/introgressed_variants_DNase_footprints.txt.gz"))

# save summ variants
introgressed_variants_DNase_footprints_summ <- introgressed_variants_DNase_footprints %>% 
	dplyr::select(seqnames, start, end, width, strand, starts_with("Variant"), contains("_summ_"))
write_tsv(introgressed_variants_DNase_footprints_summ, gzfile("../../results/2s-annotate_variants_DNase_footprints/introgressed_variants_DNase_footprints_summ.txt.gz"))
sort(table(introgressed_variants_DNase_footprints_summ$DNase_footprints_summ_Summary))
