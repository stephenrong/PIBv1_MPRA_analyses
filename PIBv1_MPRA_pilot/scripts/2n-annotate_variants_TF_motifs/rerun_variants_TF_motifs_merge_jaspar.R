#!/bin/R

# for data analysis
library(tidyverse)
library(data.table)

# for bio data analysis
library(plyranges)

# process TF motif
adaptive_variants_TF_motifs_jaspar_list <- NULL
for (index in 1:83018) {
	print(index)
	adaptive_variants_TF_motifs_jaspar_list[[index]] <- readRDS(paste0("../../results/2n-annotate_variants_TF_motifs/motifbreakR_job_outputs/adaptive_variants_TF_motifs_jaspar-", format(index, scientific=FALSE), ".rds"))
}

adaptive_variants_TF_motifs_jaspar_list <- adaptive_variants_TF_motifs_jaspar_list %>% discard(is.null)
adaptive_variants_TF_motifs_jaspar <- unlist(GRangesList(adaptive_variants_TF_motifs_jaspar_list))
saveRDS(adaptive_variants_TF_motifs_jaspar, "../../results/2n-annotate_variants_TF_motifs/adaptive_variants_TF_motifs_jaspar.rds")
