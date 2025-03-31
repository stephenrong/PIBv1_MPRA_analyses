#!/bin/R

# for data analysis
library(tidyverse)
library(data.table)

# for bio data analysis
library(plyranges)

# process TF motif
introgressed_variants_TF_motifs_hocomoco_list <- NULL
for (index in 1:713470) {
	print(index)
	introgressed_variants_TF_motifs_hocomoco_list[[index]] <- readRDS(paste0("../../results/2n-annotate_variants_TF_motifs/motifbreakR_job_outputs/introgressed_variants_TF_motifs_hocomoco-", format(index, scientific=FALSE), ".rds"))
}

introgressed_variants_TF_motifs_hocomoco_list <- introgressed_variants_TF_motifs_hocomoco_list %>% discard(is.null)
introgressed_variants_TF_motifs_hocomoco <- unlist(GRangesList(introgressed_variants_TF_motifs_hocomoco_list))
saveRDS(introgressed_variants_TF_motifs_hocomoco, "../../results/2n-annotate_variants_TF_motifs/introgressed_variants_TF_motifs_hocomoco.rds")
