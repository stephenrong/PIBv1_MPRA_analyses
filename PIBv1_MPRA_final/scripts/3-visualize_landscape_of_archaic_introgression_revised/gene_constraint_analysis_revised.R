#!/bin/R

library(tidyverse)
library(data.table)
library(genekitr)

# load scores
cutoff <- 0.10

s_het <- read_tsv("../../data/gene_constraint_metrics/media-1.tsv")[c("ensg", "post_mean")]
names(s_het) <- c("tx_ensembl_gene_id", "s_het")
s_het <- s_het %>% filter(!duplicated(tx_ensembl_gene_id))

loeuf <- read_tsv("../../data/gene_constraint_metrics/gnomad.v2.1.1.lof_metrics.by_gene.txt")[c("gene", "oe_lof_upper")]

names(loeuf) <- c("tx_hgnc_symbol", "loeuf")
loeuf <- loeuf %>% filter(!duplicated(tx_hgnc_symbol))

transId_map <- deframe(transId(sort(unique(loeuf$tx_hgnc_symbol)), transTo=c("symbol"), keepNA=TRUE, unique=TRUE, hgVersion="v38"))
loeuf <- loeuf %>% 
	rowwise() %>% mutate(tx_hgnc_symbol = ifelse(tx_hgnc_symbol %in% names(transId_map), transId_map[[tx_hgnc_symbol]], NA)) %>% ungroup()
loeuf <- loeuf %>% filter(!duplicated(tx_hgnc_symbol))

# quantile based cutoffs
s_het_bottom <- quantile(s_het$s_het, cutoff, na.rm=T)
s_het_top <- quantile(s_het$s_het, 1-cutoff, na.rm=T)

loeuf_bottom <- 0.35
loeuf_top <- 1

# load tables
skov_introgressed_deserts_win_ucsc_knownGene <- read_tsv("../../results/3-visualize_landscape_of_archaic_introgression_revised/skov_introgressed_deserts_win_ucsc_knownGene.txt.gz")[c("tx_ensembl_gene_id", "tx_hgnc_symbol")] %>% 
	distinct() %>% left_join(s_het) %>% left_join(loeuf)
skov_introgressed_reduced_ucsc_knownGene <- read_tsv("../../results/3-visualize_landscape_of_archaic_introgression_revised/skov_introgressed_reduced_ucsc_knownGene.txt.gz")[c("tx_ensembl_gene_id", "tx_hgnc_symbol")] %>% 
	distinct() %>% left_join(s_het) %>% left_join(loeuf)
introgressed_deserts_OCN_and_noOCN_ucsc_knownGene <- read_tsv("../../results/3-visualize_landscape_of_archaic_introgression_revised/introgressed_deserts_OCN_and_noOCN_ucsc_knownGene.txt.gz")[c("tx_ensembl_gene_id", "tx_hgnc_symbol")] %>% 
	distinct() %>% left_join(s_het) %>% left_join(loeuf)
introgressed_reduced_OCN_and_noOCN_ucsc_knownGene <- read_tsv("../../results/3-visualize_landscape_of_archaic_introgression_revised/introgressed_reduced_OCN_and_noOCN_ucsc_knownGene.txt.gz")[c("tx_ensembl_gene_id", "tx_hgnc_symbol")] %>% 
	distinct() %>% left_join(s_het) %>% left_join(loeuf)

# print for revision
write_tsv(skov_introgressed_deserts_win_ucsc_knownGene, "../../results/3-visualize_landscape_of_archaic_introgression_revised/skov_introgressed_deserts_win_ucsc_knownGene_s_het_loeuf.txt")
write_tsv(skov_introgressed_reduced_ucsc_knownGene, "../../results/3-visualize_landscape_of_archaic_introgression_revised/skov_introgressed_reduced_ucsc_knownGene_s_het_loeuf.txt")
write_tsv(introgressed_deserts_OCN_and_noOCN_ucsc_knownGene, "../../results/3-visualize_landscape_of_archaic_introgression_revised/introgressed_deserts_OCN_and_noOCN_ucsc_knownGene_s_het_loeuf.txt")
write_tsv(introgressed_reduced_OCN_and_noOCN_ucsc_knownGene, "../../results/3-visualize_landscape_of_archaic_introgression_revised/introgressed_reduced_OCN_and_noOCN_ucsc_knownGene_s_het_loeuf.txt")

# fisher's exact tests
skov_introgressed_deserts_win_ucsc_knownGene_summ <- skov_introgressed_deserts_win_ucsc_knownGene %>% 
	summarise(
		s_het_bottom_n = sum(s_het < s_het_bottom, na.rm=T), 
		s_het_top_n = sum(s_het > s_het_top, na.rm=T), 
		loeuf_bottom_n = sum(loeuf > loeuf_bottom, na.rm=T), 
		loeuf_top_n = sum(loeuf < loeuf_top, na.rm=T)
	)
skov_introgressed_reduced_ucsc_knownGene_summ <- skov_introgressed_reduced_ucsc_knownGene %>% 
	summarise(
		s_het_bottom_n = sum(s_het < s_het_bottom, na.rm=T), 
		s_het_top_n = sum(s_het > s_het_top, na.rm=T), 
		loeuf_bottom_n = sum(loeuf > loeuf_bottom, na.rm=T), 
		loeuf_top_n = sum(loeuf < loeuf_top, na.rm=T)
	)
introgressed_deserts_OCN_and_noOCN_ucsc_knownGene_summ <- introgressed_deserts_OCN_and_noOCN_ucsc_knownGene %>% 
	summarise(
		s_het_bottom_n = sum(s_het < s_het_bottom, na.rm=T), 
		s_het_top_n = sum(s_het > s_het_top, na.rm=T), 
		loeuf_bottom_n = sum(loeuf > loeuf_bottom, na.rm=T), 
		loeuf_top_n = sum(loeuf < loeuf_top, na.rm=T)
	)
introgressed_reduced_OCN_and_noOCN_ucsc_knownGene_summ <- introgressed_reduced_OCN_and_noOCN_ucsc_knownGene %>% 
	summarise(
		s_het_bottom_n = sum(s_het < s_het_bottom, na.rm=T), 
		s_het_top_n = sum(s_het > s_het_top, na.rm=T), 
		loeuf_bottom_n = sum(loeuf > loeuf_bottom, na.rm=T), 
		loeuf_top_n = sum(loeuf < loeuf_top, na.rm=T)
	)

fisher_OCN_noOCN_s_het <- fisher.test(matrix(c(
	introgressed_deserts_OCN_and_noOCN_ucsc_knownGene_summ$s_het_top_n, 
	introgressed_deserts_OCN_and_noOCN_ucsc_knownGene_summ$s_het_bottom_n, 
	introgressed_reduced_OCN_and_noOCN_ucsc_knownGene_summ$s_het_top_n, 
	introgressed_reduced_OCN_and_noOCN_ucsc_knownGene_summ$s_het_bottom_n
), nrow=2))
sink(paste0("../../results/3-visualize_landscape_of_archaic_introgression_revised/fisher_OCN_noOCN_s_het_", cutoff, ".txt"))
fisher_OCN_noOCN_s_het
sink()

fisher_skov_s_het <- fisher.test(matrix(c(
	skov_introgressed_deserts_win_ucsc_knownGene_summ$s_het_top_n, 
	skov_introgressed_deserts_win_ucsc_knownGene_summ$s_het_bottom_n, 
	skov_introgressed_reduced_ucsc_knownGene_summ$s_het_top_n, 
	skov_introgressed_reduced_ucsc_knownGene_summ$s_het_bottom_n
), nrow=2))
sink(paste0("../../results/3-visualize_landscape_of_archaic_introgression_revised/fisher_skov_s_het_", cutoff, ".txt"))
fisher_skov_s_het
sink()

fisher_OCN_noOCN_loeuf <- fisher.test(matrix(c(
	introgressed_deserts_OCN_and_noOCN_ucsc_knownGene_summ$loeuf_top_n, 
	introgressed_deserts_OCN_and_noOCN_ucsc_knownGene_summ$loeuf_bottom_n, 
	introgressed_reduced_OCN_and_noOCN_ucsc_knownGene_summ$loeuf_top_n, 
	introgressed_reduced_OCN_and_noOCN_ucsc_knownGene_summ$loeuf_bottom_n
), nrow=2))
sink(paste0("../../results/3-visualize_landscape_of_archaic_introgression_revised/fisher_OCN_noOCN_loeuf_fixed.txt"))
fisher_OCN_noOCN_loeuf
sink()

fisher_skov_loeuf <- fisher.test(matrix(c(
	skov_introgressed_deserts_win_ucsc_knownGene_summ$loeuf_top_n, 
	skov_introgressed_deserts_win_ucsc_knownGene_summ$loeuf_bottom_n, 
	skov_introgressed_reduced_ucsc_knownGene_summ$loeuf_top_n, 
	skov_introgressed_reduced_ucsc_knownGene_summ$loeuf_bottom_n
), nrow=2))
sink(paste0("../../results/3-visualize_landscape_of_archaic_introgression_revised/fisher_skov_loeuf_fixed.txt"))
fisher_skov_loeuf
sink()
