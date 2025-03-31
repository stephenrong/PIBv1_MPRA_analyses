#!/bin/R

# for data analysis
library(tidyverse)
library(data.table)

# for bio data analysis
library(plyranges)

# process TF motif
introgressed_variants_TF_motifs_jaspar_list <- NULL
for (index in 1:713470) {
	print(index)
	introgressed_variants_TF_motifs_jaspar_list[[index]] <- readRDS(paste0("../../results/2n-annotate_variants_TF_motifs/motifbreakR_job_outputs/introgressed_variants_TF_motifs_jaspar-", format(index, scientific=FALSE), ".rds"))
}

introgressed_variants_TF_motifs_jaspar_list <- introgressed_variants_TF_motifs_jaspar_list %>% discard(is.null)
introgressed_variants_TF_motifs_jaspar <- unlist(GRangesList(introgressed_variants_TF_motifs_jaspar_list))
saveRDS(introgressed_variants_TF_motifs_jaspar, "../../results/2n-annotate_variants_TF_motifs/introgressed_variants_TF_motifs_jaspar.rds")
