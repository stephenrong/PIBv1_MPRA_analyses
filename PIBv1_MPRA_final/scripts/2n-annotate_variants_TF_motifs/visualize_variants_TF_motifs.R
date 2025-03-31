#!/bin/R

# for data analysis
library(tidyverse)
library(data.table)

# for bio data analysis
library(plyranges)
library(motifbreakR)
library(BSgenome.Hsapiens.UCSC.hg19)
source("../shared_functions/seqinfo_fix_change.R")

# load motifbreakR input
introgressed_variants <- readRDS("../../results/1a-preprocess_PIBv1_MPRA_final/introgressed_variants.rds")

# temporary filters
introgressed_variants <- introgressed_variants %>% 
	filter(VariantID %in% c(
		"12_113346546_C_G", "12_113347017_T_C",  # OAS
		"1_89587575_T_C", "1_89592457_C_T", "1_89592489_T_C", "1_89594982_T_G", "1_89596265_G_A", "1_89599705_A_G",  # GBP
		"6_138020511_C_T", "6_138037108_G_A", "6_138039111_A_T", "6_138226361_A_C", "6_138248445_G_A"
		)
	)

# fix motifbreakR input
introgressed_variants_TF_motifs <- introgressed_variants %>% 
	mutate(SNP_id = VariantID, REF = VariantREF, ALT = VariantALT) %>% 
	as_tibble() %>% dplyr::select(seqnames, start, end, width, strand, SNP_id, REF, ALT) %>% GRanges()

introgressed_variants_TF_motifs <- introgressed_variants_TF_motifs %>% 
	seqinfo_fix("UCSC", "hg19")  # necessary to match BSgenome

attr(introgressed_variants_TF_motifs, "genome.package") <- 
	"BSgenome.Hsapiens.UCSC.hg19"  # necessary for motifbreakR input

# get motifbreakR results, hocomoco
data(hocomoco)
introgressed_variants_TF_motifs_hocomoco <- motifbreakR(
	snpList = introgressed_variants_TF_motifs,  # snps.mb.temp,
	filterp = TRUE,
	pwmList = hocomoco, 
	threshold = 1e-4, 
	method = "ic",
	bkg = c(A=0.25, C=0.25, G=0.25, T=0.25),
	BPPARAM = BiocParallel::SerialParam()
)

for (i in 1:length(introgressed_variants_TF_motifs)) {
	if (introgressed_variants_TF_motifs$SNP_id[i] %in% introgressed_variants_TF_motifs_hocomoco$SNP_id) {
		pdf(paste0("../../results/2n-annotate_variants_TF_motifs/variant_level_visualizations/variant_hocomoco-", introgressed_variants_TF_motifs$SNP_id[i], ".pdf"))  # , width=4.5, height=4.5)
		plotMB(results = introgressed_variants_TF_motifs_hocomoco, rsid = introgressed_variants_TF_motifs$SNP_id[i], altAllele = introgressed_variants_TF_motifs$ALT[i])
		dev.off()
		pdf(paste0("../../results/2n-annotate_variants_TF_motifs/variant_level_visualizations/variant_hocomoco_max-", introgressed_variants_TF_motifs$SNP_id[i], ".pdf"))  # , width=4.5, height=4.5)
		plotMB(results = (introgressed_variants_TF_motifs_hocomoco %>% group_by(SNP_id) %>% filter(abs(alleleDiff) == max(abs(alleleDiff))) %>% ungroup()), rsid = introgressed_variants_TF_motifs$SNP_id[i], altAllele = introgressed_variants_TF_motifs$ALT[i])
		dev.off()
	}
}

# get motifbreakR results, jaspar
introgressed_variants_TF_motifs_results_jaspar <- motifbreakR(
	snpList = introgressed_variants_TF_motifs,  # snps.mb.temp,
	filterp = TRUE,
	pwmList = subset(subset(MotifDb, dataSource %in% c("jaspar2022")), organism == "Hsapiens"), 
	threshold = 1e-4, 
	method = "ic",
	bkg = c(A=0.25, C=0.25, G=0.25, T=0.25),
	BPPARAM = BiocParallel::SerialParam()
)

for (i in 1:length(introgressed_variants_TF_motifs)) {
	if (introgressed_variants_TF_motifs$SNP_id[i] %in% introgressed_variants_TF_motifs_results_jaspar$SNP_id) {
		pdf(paste0("../../results/2n-annotate_variants_TF_motifs/variant_level_visualizations/variant_jaspar-", introgressed_variants_TF_motifs$SNP_id[i], ".pdf"))  # , width=4.5, height=4.5)
		plotMB(results = introgressed_variants_TF_motifs_results_jaspar, rsid = introgressed_variants_TF_motifs$SNP_id[i], altAllele = introgressed_variants_TF_motifs$ALT[i])
		dev.off()
		pdf(paste0("../../results/2n-annotate_variants_TF_motifs/variant_level_visualizations/variant_jaspar_max-", introgressed_variants_TF_motifs$SNP_id[i], ".pdf"))  # , width=4.5, height=4.5)
		plotMB(results = (introgressed_variants_TF_motifs_results_jaspar %>% group_by(SNP_id) %>% filter(abs(alleleDiff) == max(abs(alleleDiff))) %>% ungroup()), rsid = introgressed_variants_TF_motifs$SNP_id[i], altAllele = introgressed_variants_TF_motifs$ALT[i])
		dev.off()
	}
}
