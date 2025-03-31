#!/bin/R

# for data analysis
library(tidyverse)
library(data.table)

# for bio data analysis
library(plyranges)
library(motifbreakR)
library(BSgenome.Hsapiens.UCSC.hg19)
source("../shared_functions/seqinfo_fix_change.R")

# command line
index = commandArgs(trailingOnly=TRUE)
i <- as.integer(index[1])
j <- as.integer(index[2])

# load variants
adaptive_variants <- as_tibble(readRDS("../../results/1a-preprocess_PIBv1_MPRA_pilot/adaptive_variants.rds"))[,1:10]

# run motifbreakR
adaptive_variants_TF_motifs <- adaptive_variants %>% 
	mutate(SNP_id = VariantID, REF = VariantREF, ALT = VariantALT) %>% 
	as_tibble() %>% dplyr::select(seqnames, start, end, width, strand, SNP_id, REF, ALT) %>% GRanges()

adaptive_variants_TF_motifs <- adaptive_variants_TF_motifs %>% 
	seqinfo_fix("UCSC", "hg19")  # necessary to match BSgenome

attr(adaptive_variants_TF_motifs, "genome.package") <- 
	"BSgenome.Hsapiens.UCSC.hg19"  # necessary for motifbreakR input

# get motifbreakR results, hocomoco from Hsapiens (426 models)
print("getting motifbreakR results for hocomoco")
data(hocomoco)
for (k in i:j) {
	print(k)
	adaptive_variants_TF_motifs_hocomoco <- motifbreakR(
		snpList = adaptive_variants_TF_motifs[k,],  # snps.mb.temp,
		filterp = TRUE,
		pwmList = hocomoco, 
		threshold = 1e-4, 
		method = "ic",
		bkg = c(A=0.25, C=0.25, G=0.25, T=0.25),
		BPPARAM = BiocParallel::SerialParam()
	)
	saveRDS(adaptive_variants_TF_motifs_hocomoco, paste0("../../results/2n-annotate_variants_TF_motifs/motifbreakR_job_outputs/adaptive_variants_TF_motifs_hocomoco-", k, ".rds"))
}

# get motifbreakR results, jaspar restricted to Hsapiens (691 models)
print("getting motifbreakR results for jaspar")
jaspar <- subset(subset(MotifDb, dataSource %in% c("jaspar2022")), organism == "Hsapiens")
for (k in i:j) {
	print(k)
	adaptive_variants_TF_motifs_jaspar <- motifbreakR(
		snpList = adaptive_variants_TF_motifs[k,],  # snps.mb.temp,
		filterp = TRUE,
		pwmList = jaspar, 
		threshold = 1e-4, 
		method = "ic",
		bkg = c(A=0.25, C=0.25, G=0.25, T=0.25),
		BPPARAM = BiocParallel::SerialParam()
	)
	saveRDS(adaptive_variants_TF_motifs_jaspar, paste0("../../results/2n-annotate_variants_TF_motifs/motifbreakR_job_outputs/adaptive_variants_TF_motifs_jaspar-", k, ".rds"))
}
