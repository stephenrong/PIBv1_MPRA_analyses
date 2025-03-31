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
adaptive_variants <- readRDS("../../results/1a-preprocess_PIBv1_MPRA_pilot/adaptive_variants.rds")

# temporary filters
adaptive_variants <- adaptive_variants %>% 
	filter(VariantID %in% c(
		"12_113346546_C_G", "12_113347017_T_C", "12_113348921_G_T", "12_113348935_G_C",  # ,  # OAS1
		"1_89587575_T_C", "1_89592457_C_T", "1_89592489_T_C", "1_89594982_T_G", "1_89596265_G_A", "1_89599705_A_G",  # GBP2
		"6_138020511_C_T", "6_138037108_G_A", "6_138039111_A_T", "6_138226361_A_C", "6_138248445_G_A", "6_138786583_G_C", #,  # TNFAIP3
		"8_116327463_A_G", "8_116364887_T_C", "8_116393745_G_A", "8_116394033_G_T",  # TRPS1
		"1_65419060_C_T", "1_65431457_G_A", "1_65475897_G_T", "1_65475965_A_G", "1_65475981_T_A", "1_65468183_T_A", "1_65477791_C_T",  #,  # JAK1
		"9_96576740_C_G",  # BARX1
		"2_161933715_T_C", "2_162016242_G_A", "2_162020918_T_C"  # TANK
		)
	)

# fix motifbreakR input
adaptive_variants_TF_motifs <- adaptive_variants %>% 
	mutate(SNP_id = VariantID, REF = VariantREF, ALT = VariantALT) %>% 
	as_tibble() %>% dplyr::select(seqnames, start, end, width, strand, SNP_id, REF, ALT) %>% GRanges()

adaptive_variants_TF_motifs <- adaptive_variants_TF_motifs %>% 
	seqinfo_fix("UCSC", "hg19")  # necessary to match BSgenome

attr(adaptive_variants_TF_motifs, "genome.package") <- 
	"BSgenome.Hsapiens.UCSC.hg19"  # necessary for motifbreakR input

# get motifbreakR results, hocomoco + jaspar
data(hocomoco)
adaptive_variants_TF_motifs_results_hocomoco_and_jaspar <- motifbreakR(
	snpList = adaptive_variants_TF_motifs,  # snps.mb.temp,
	filterp = TRUE,
	pwmList = c(subset(subset(MotifDb, dataSource %in% c("jaspar2022")), organism == "Hsapiens"), hocomoco), 
	threshold = 1e-4, 
	method = "ic",
	bkg = c(A=0.25, C=0.25, G=0.25, T=0.25),
	BPPARAM = BiocParallel::SerialParam()
)

for (i in 1:length(adaptive_variants_TF_motifs)) {
	if (adaptive_variants_TF_motifs$SNP_id[i] %in% adaptive_variants_TF_motifs_results_hocomoco_and_jaspar$SNP_id) {
		pdf(paste0("../../results/2n-annotate_variants_TF_motifs/variant_level_visualizations/variant_hocomoco_and_jaspar-", adaptive_variants_TF_motifs$SNP_id[i], ".pdf"), width=10, height=10)
		plotMB(results = adaptive_variants_TF_motifs_results_hocomoco_and_jaspar, rsid = adaptive_variants_TF_motifs$SNP_id[i], altAllele = adaptive_variants_TF_motifs$ALT[i])
		dev.off()
		if (length(adaptive_variants_TF_motifs_results_hocomoco_and_jaspar %>% filter(SNP_id == adaptive_variants_TF_motifs$SNP_id[i])) > 5) {
			pdf(paste0("../../results/2n-annotate_variants_TF_motifs/variant_level_visualizations/variant_hocomoco_and_jaspar-", adaptive_variants_TF_motifs$SNP_id[i], ".pdf"), width=30, height=30)
			plotMB(results = adaptive_variants_TF_motifs_results_hocomoco_and_jaspar, rsid = adaptive_variants_TF_motifs$SNP_id[i], altAllele = adaptive_variants_TF_motifs$ALT[i])
			dev.off()
		}
		pdf(paste0("../../results/2n-annotate_variants_TF_motifs/variant_level_visualizations/variant_hocomoco_and_jaspar_max-", adaptive_variants_TF_motifs$SNP_id[i], ".pdf"))
		plotMB(results = (adaptive_variants_TF_motifs_results_hocomoco_and_jaspar %>% group_by(SNP_id) %>% filter(alleleEffectSize == max(alleleEffectSize)) %>% ungroup()), rsid = adaptive_variants_TF_motifs$SNP_id[i], altAllele = adaptive_variants_TF_motifs$ALT[i])
		dev.off()
		pdf(paste0("../../results/2n-annotate_variants_TF_motifs/variant_level_visualizations/variant_hocomoco_and_jaspar_min-", adaptive_variants_TF_motifs$SNP_id[i], ".pdf"))
		plotMB(results = (adaptive_variants_TF_motifs_results_hocomoco_and_jaspar %>% group_by(SNP_id) %>% filter(alleleEffectSize == min(alleleEffectSize)) %>% ungroup()), rsid = adaptive_variants_TF_motifs$SNP_id[i], altAllele = adaptive_variants_TF_motifs$ALT[i])
		dev.off()
	}
}

# visualize
for (i in adaptive_variants_TF_motifs_results_hocomoco_and_jaspar$SNP_id) {
	adaptive_variants_TF_motifs_results_hocomoco_and_jaspar_temp <- adaptive_variants_TF_motifs_results_hocomoco_and_jaspar %>% 
		group_by(SNP_id) %>% filter(alleleEffectSize == max(alleleEffectSize)) %>% ungroup() %>% 
		as_tibble() %>% filter(SNP_id == i)
	ggplot(adaptive_variants_TF_motifs_results_hocomoco_and_jaspar_temp) + 
		geom_bar(aes(1, pctRef), stat="identity", fill="#fee5d9") + 
		geom_bar(aes(2, pctAlt), stat="identity", fill="#fb6a4a") + 
		scale_x_discrete("Allele", labels = c("REF", "ALT"), limits = c(1, 2)) + 
		ylab("TF motif percent of best score") + 
		theme(
			aspect.ratio = 2, 
			panel.grid.major = element_blank(), 
			panel.grid.minor = element_blank(),
			panel.background = element_rect(fill = "transparent", color = NA),  
			plot.background = element_rect(fill = "transparent", color = NA), 
			axis.line = element_line(color = "black")
		)
	ggsave(paste0("../../results/2n-annotate_variants_TF_motifs/variant_level_visualizations/variant_hocomoco_and_jaspar_max_activity-", i, ".pdf"), width=3, height=3)
	adaptive_variants_TF_motifs_results_hocomoco_and_jaspar_temp <- adaptive_variants_TF_motifs_results_hocomoco_and_jaspar %>% 
		group_by(SNP_id) %>% filter(alleleEffectSize == min(alleleEffectSize)) %>% ungroup() %>% 
		as_tibble() %>% filter(SNP_id == i)
	ggplot(adaptive_variants_TF_motifs_results_hocomoco_and_jaspar_temp) + 
		geom_bar(aes(1, pctRef), stat="identity", fill="#fee5d9") + 
		geom_bar(aes(2, pctAlt), stat="identity", fill="#fb6a4a") + 
		scale_x_discrete("Allele", labels = c("REF", "ALT"), limits = c(1, 2)) + 
		ylab("TF motif percent of best score") + 
		theme(
			aspect.ratio = 2, 
			panel.grid.major = element_blank(), 
			panel.grid.minor = element_blank(),
			panel.background = element_rect(fill = "transparent", color = NA),  
			plot.background = element_rect(fill = "transparent", color = NA), 
			axis.line = element_line(color = "black")
		)
	ggsave(paste0("../../results/2n-annotate_variants_TF_motifs/variant_level_visualizations/variant_hocomoco_and_jaspar_min_activity-", i, ".pdf"), width=3, height=3)
}
