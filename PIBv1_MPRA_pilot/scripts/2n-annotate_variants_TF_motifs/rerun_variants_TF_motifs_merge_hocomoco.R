#!/bin/R

# for data analysis
library(tidyverse)
library(data.table)

# for bio data analysis
library(plyranges)

# process TF motif
adaptive_variants_TF_motifs_hocomoco_list <- NULL
for (index in 1:83018) {
	print(index)
	adaptive_variants_TF_motifs_hocomoco_list[[index]] <- readRDS(paste0("../../results/2n-annotate_variants_TF_motifs/motifbreakR_job_outputs/adaptive_variants_TF_motifs_hocomoco-", format(index, scientific=FALSE), ".rds"))
}

adaptive_variants_TF_motifs_hocomoco_list <- adaptive_variants_TF_motifs_hocomoco_list %>% discard(is.null)
adaptive_variants_TF_motifs_hocomoco <- unlist(GRangesList(adaptive_variants_TF_motifs_hocomoco_list))
saveRDS(adaptive_variants_TF_motifs_hocomoco, "../../results/2n-annotate_variants_TF_motifs/adaptive_variants_TF_motifs_hocomoco.rds")
