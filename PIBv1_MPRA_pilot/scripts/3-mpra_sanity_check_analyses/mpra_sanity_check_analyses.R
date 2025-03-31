#!/bin/R

# load packages
library(tidyverse)
library(data.table)
library(eulerr)
library(ggpubr)
library(scales)

# load bio packages
library(plyranges)

# for aesthetics
library(RColorBrewer)
library(hrbrthemes)
library(viridis)
library(ggrepel)

# load files
mpra_ccre_variants <- as_tibble(fread("../../results/3-merge_all_variant_annotations/mpra_ccre_variants_all_annotations_summ.txt.gz"))

# emVars and ENCODE cCREs
ENCODE_cCREs_other_map <- c(
	"CA"="CA",
	"CA-CTCF"="CA",
	"CA-H3K4me3"="CA",
	"CA-TF"="CA",
	"CA-TF,CA"="CA",
	"dELS"="dELS",
	"dELS,CA-CTCF"="dELS",
	"dELS,dELS"="dELS",
	"pELS"="pELS",
	"PLS"="PLS",
	"TF"="TF"
)

mpra_emVars_K562_ENCODE_cCREs <- mpra_ccre_variants %>% 
	filter(!is.na(MPRA_K562_summ_Summary)) %>% 
	mutate(ENCODE_cCREs_summ_Summary = ENCODE_cCREs_other_map[ENCODE_cCREs_summ_Summary]) %>% 
	dplyr::select(ENCODE_cCREs_summ_Summary) %>% table() %>% enframe() %>% 
	mutate(percent = 100*value/sum(value)) %>% 
	mutate(name = factor(name, levels=(c(
		"CA", 
		"TF", 
		"dELS", 
		"pELS",
		"PLS"
	))))

palette <- (c("#eff3ff", "#bdd7e7", "#6baed6", "#3182bd", "#08519c"))
ggplot(mpra_emVars_K562_ENCODE_cCREs) + 
	aes(0.4, value, fill = name, label = name) + 
	geom_bar(
		position = "stack", 
		stat = "identity", 
		color = "white", 
		linewidth = 0.2
	) + 
	coord_flip(clip = "off") + 
	geom_text(
		aes(1, value, color = name, 
			label = name
		), 
		color = "black", 
		position = position_stack(vjust = 0.5),
		hjust = 0.1, 
		vjust = -0.1,
		size = 3,
		angle = 45
	) + 
	geom_text(
		aes(0.4, value, color = name, 
			label = format(value, nsmall=1, big.mark=",")
		), 
		color = "black", 
		position = position_stack(vjust = 0.5),
		hjust = 0.5, 
		vjust = 0.5,
		size = 2,
		angle = 90
	) + 
	scale_fill_manual(name = "", values = palette) + 
	scale_color_manual(name = "", values = palette) + 
	theme_ipsum(base_family = "Arial") + 
	scale_x_discrete(labels = NULL, breaks = NULL) + 
	scale_y_continuous(breaks = NULL) + 
	labs(x = "", y = "") + 
	theme(
		aspect.ratio = 0.33, 
		legend.position = "none", 
		panel.grid.major = element_blank(), 
		panel.grid.minor = element_blank(),
		panel.background = element_rect(fill = "transparent", color = NA), 
		plot.background = element_rect(fill = "transparent", color = NA)
	)
ggsave("../../results/3-mpra_sanity_check_analyses/mpra_emVars_ENCODE_cCREs-K562.pdf", device=cairo_pdf, scale=0.6)

mpra_emVars_Jurkat_ENCODE_cCREs <- mpra_ccre_variants %>% 
	filter(!is.na(MPRA_Jurkat_summ_Summary)) %>% 
	mutate(ENCODE_cCREs_summ_Summary = ENCODE_cCREs_other_map[ENCODE_cCREs_summ_Summary]) %>% 
	dplyr::select(ENCODE_cCREs_summ_Summary) %>% table() %>% enframe() %>% 
	mutate(percent = 100*value/sum(value)) %>% 
	mutate(name = factor(name, levels=(c(
		"CA", 
		"TF", 
		"dELS", 
		"pELS",
		"PLS"
	))))

palette <- (c("#eff3ff", "#bdd7e7", "#6baed6", "#3182bd", "#08519c"))
ggplot(mpra_emVars_Jurkat_ENCODE_cCREs) + 
	aes(0.4, value, fill = name, label = name) + 
	geom_bar(
		position = "stack", 
		stat = "identity", 
		color = "white", 
		linewidth = 0.2
	) + 
	coord_flip(clip = "off") + 
	geom_text(
		aes(1, value, color = name, 
			label = name
		), 
		color = "black", 
		position = position_stack(vjust = 0.5),
		hjust = 0.1, 
		vjust = -0.1,
		size = 3,
		angle = 45
	) + 
	geom_text(
		aes(0.4, value, color = name, 
			label = format(value, nsmall=1, big.mark=",")
		), 
		color = "black", 
		position = position_stack(vjust = 0.5),
		hjust = 0.5, 
		vjust = 0.5,
		size = 2,
		angle = 90
	) + 
	scale_fill_manual(name = "", values = palette) + 
	scale_color_manual(name = "", values = palette) + 
	theme_ipsum(base_family = "Arial") + 
	scale_x_discrete(labels = NULL, breaks = NULL) + 
	scale_y_continuous(breaks = NULL) + 
	labs(x = "", y = "") + 
	theme(
		aspect.ratio = 0.33, 
		legend.position = "none", 
		panel.grid.major = element_blank(), 
		panel.grid.minor = element_blank(),
		panel.background = element_rect(fill = "transparent", color = NA), 
		plot.background = element_rect(fill = "transparent", color = NA)
	)
ggsave("../../results/3-mpra_sanity_check_analyses/mpra_emVars_ENCODE_cCREs-Jurkat.pdf", device=cairo_pdf, scale=0.6)


# overlap between mpras emVars
mpra_emVars_venn <- mpra_ccre_variants %>% 
	filter(!is.na(MPRA_K562_summ_Summary) | !is.na(MPRA_Jurkat_summ_Summary)) %>% 
	mutate(`K562 emVars` = !is.na(MPRA_K562_summ_Summary), `Jurkat emVars` = !is.na(MPRA_Jurkat_summ_Summary)) %>% 
	dplyr::select(`K562 emVars`, `Jurkat emVars`)

pdf("../../results/3-mpra_sanity_check_analyses/mpra_emVars_euler-K562-Jurkat.pdf", height=5, width=5)
plot(euler(mpra_emVars_venn), quantities=TRUE, legend=TRUE)
dev.off()

# correlations between mpra log skews
ggscatter(mpra_ccre_variants %>% filter(!is.na(MPRA_K562_summ_Summary) | !is.na(MPRA_Jurkat_summ_Summary)), 
	x = "MPRA_K562_summ_Description_skew_log2FC", 
	y = "MPRA_Jurkat_summ_Description_skew_log2FC", 
	color = "black", shape = 1, size = 3, 
	add = "reg.line",
	add.params = list(color = "blue", fill = "lightgray"), 
	conf.int = TRUE,
	cor.coef = TRUE,
	cor.coeff.args = list(method = "pearson", label.sep = "\n")
) + 
	geom_hline(yintercept = 0, color = "orange") + 
	geom_vline(xintercept = 0, color = "orange") + 
	theme_bw() + theme(
		aspect.ratio = 1, 
		axis.line = element_line(colour = "black"),
		panel.grid.major = element_blank(),
		panel.grid.minor = element_blank(),
		panel.border = element_blank(),
		panel.background = element_blank()
	) + 
	coord_fixed(ratio = 1) + 
	xlab("K562 log2 Skew") + 
	ylab("Jurkat log2 Skew") + 
	scale_x_continuous(breaks = pretty_breaks(n=5)) + 
	scale_y_continuous(breaks = pretty_breaks(n=5))
ggsave("../../results/3-mpra_sanity_check_analyses/mpra_correlations-K562_log2FC-or-Jurkat_log2FC.pdf", scale=0.6)

# enrichment of active sequences in ENCODE cCRE cell lines
mpra_active_K562_ENCODE_DHS <- mpra_ccre_variants %>% 
	filter(!is.na(ENCODE_cCREs_summ_Summary)) %>% 
	mutate(MPRA_K562_active = !is.na(MPRA_K562_summ_Description_active)) %>% 
	mutate(K562_ENCODE_DHS = !(grepl("Low-DNase", ENCODE_cCREs_K562_summ_Category))) %>% 
	dplyr::select(MPRA_K562_active, K562_ENCODE_DHS)
mpra_active_K562_ENCODE_DHS_tab <- table(mpra_active_K562_ENCODE_DHS)
mpra_active_K562_ENCODE_DHS_fisher <- fisher.test(mpra_active_K562_ENCODE_DHS_tab)
mpra_active_K562_ENCODE_DHS_fisher_tb <- as_tibble(t(c(
	"P-value"=mpra_active_K562_ENCODE_DHS_fisher$p.value, 
	"Fold enrichment"=mpra_active_K562_ENCODE_DHS_fisher$estimate[["odds ratio"]], 
	"Number of active"=sum(mpra_active_K562_ENCODE_DHS$MPRA_K562_active
)))) %>% mutate(`Cell line`="K562")

mpra_active_Jurkat_ENCODE_DHS <- mpra_ccre_variants %>% 
	filter(!is.na(ENCODE_cCREs_summ_Summary)) %>% 
	mutate(MPRA_Jurkat_active = !is.na(MPRA_Jurkat_summ_Description_active)) %>% 
	mutate(Jurkat_ENCODE_DHS = !(grepl("Low-DNase", ENCODE_cCREs_Jurkat_summ_Category))) %>% 
	dplyr::select(MPRA_Jurkat_active, Jurkat_ENCODE_DHS)
mpra_active_Jurkat_ENCODE_DHS_tab <- table(mpra_active_Jurkat_ENCODE_DHS)
mpra_active_Jurkat_ENCODE_DHS_fisher <- fisher.test(mpra_active_Jurkat_ENCODE_DHS_tab)
mpra_active_Jurkat_ENCODE_DHS_fisher_tb <- as_tibble(t(c(
	"P-value"=mpra_active_Jurkat_ENCODE_DHS_fisher$p.value, 
	"Fold enrichment"=mpra_active_Jurkat_ENCODE_DHS_fisher$estimate[["odds ratio"]], 
	"Number of active"=sum(mpra_active_Jurkat_ENCODE_DHS$MPRA_Jurkat_active
)))) %>% mutate(`Cell line`="Jurkat")


mpra_active_ENCODE_DHS_fisher_tb <- bind_rows(
	mpra_active_K562_ENCODE_DHS_fisher_tb, 
	mpra_active_Jurkat_ENCODE_DHS_fisher_tb
) %>% 
	mutate(`Cell line` = factor(`Cell line`, levels = c("K562", "Jurkat")))

ggplot(mpra_active_ENCODE_DHS_fisher_tb) + 
	geom_bar(aes(`Cell line`, `Fold enrichment`), stat="identity", fill="#6baed6") + 
	geom_text(aes(`Cell line`, 0, label=`Number of active`), vjust=1, color="#1b9e77") + 
	geom_text(aes(`Cell line`, `Fold enrichment`, label=round(-log10(`P-value`), 2)), vjust=0, color="#7570b3") + 
	geom_hline(yintercept=1, color="#eeac6f") + 
	theme_bw() + 
	ggtitle("Enrichment of active \n sequences in ENCODE DHS regions") +
	theme(
		aspect.ratio=1.25, 
		panel.grid.major = element_blank(),
		panel.grid.minor = element_blank(), 
		plot.title = element_text(hjust = 0.5)
	)
ggsave("../../results/3-mpra_sanity_check_analyses/mpra_active_ENCODE_DHS_enrichment-K562-Jurkat.pdf", scale=0.6)


# enrichment of emVars in ENCODE cCRE cell lines
mpra_emVars_K562_ENCODE_DHS <- mpra_ccre_variants %>% 
	filter(!is.na(ENCODE_cCREs_summ_Summary)) %>% 
	mutate(MPRA_K562_emVars = !is.na(MPRA_K562_summ_Summary)) %>% 
	mutate(K562_ENCODE_DHS = !(grepl("Low-DNase", ENCODE_cCREs_K562_summ_Category))) %>% 
	dplyr::select(MPRA_K562_emVars, K562_ENCODE_DHS)
mpra_emVars_K562_ENCODE_DHS_tab <- table(mpra_emVars_K562_ENCODE_DHS)
mpra_emVars_K562_ENCODE_DHS_fisher <- fisher.test(mpra_emVars_K562_ENCODE_DHS_tab)
mpra_emVars_K562_ENCODE_DHS_fisher_tb <- as_tibble(t(c(
	"P-value"=mpra_emVars_K562_ENCODE_DHS_fisher$p.value, 
	"Fold enrichment"=mpra_emVars_K562_ENCODE_DHS_fisher$estimate[["odds ratio"]], 
	"Number of emVars"=sum(mpra_emVars_K562_ENCODE_DHS$MPRA_K562_emVars
)))) %>% mutate(`Cell line`="K562")

mpra_emVars_Jurkat_ENCODE_DHS <- mpra_ccre_variants %>% 
	filter(!is.na(ENCODE_cCREs_summ_Summary)) %>% 
	mutate(MPRA_Jurkat_emVars = !is.na(MPRA_Jurkat_summ_Summary)) %>% 
	mutate(Jurkat_ENCODE_DHS = !(grepl("Low-DNase", ENCODE_cCREs_Jurkat_summ_Category))) %>% 
	dplyr::select(MPRA_Jurkat_emVars, Jurkat_ENCODE_DHS)
mpra_emVars_Jurkat_ENCODE_DHS_tab <- table(mpra_emVars_Jurkat_ENCODE_DHS)
mpra_emVars_Jurkat_ENCODE_DHS_fisher <- fisher.test(mpra_emVars_Jurkat_ENCODE_DHS_tab)
mpra_emVars_Jurkat_ENCODE_DHS_fisher_tb <- as_tibble(t(c(
	"P-value"=mpra_emVars_Jurkat_ENCODE_DHS_fisher$p.value, 
	"Fold enrichment"=mpra_emVars_Jurkat_ENCODE_DHS_fisher$estimate[["odds ratio"]], 
	"Number of emVars"=sum(mpra_emVars_Jurkat_ENCODE_DHS$MPRA_Jurkat_emVars
)))) %>% mutate(`Cell line`="Jurkat")


mpra_emVars_ENCODE_DHS_fisher_tb <- bind_rows(
	mpra_emVars_K562_ENCODE_DHS_fisher_tb, 
	mpra_emVars_Jurkat_ENCODE_DHS_fisher_tb
) %>% 
	mutate(`Cell line`=factor(`Cell line`, levels = c("K562", "Jurkat")))

ggplot(mpra_emVars_ENCODE_DHS_fisher_tb) + 
	geom_bar(aes(`Cell line`, `Fold enrichment`), stat="identity", fill="#6baed6") + 
	geom_text(aes(`Cell line`, 0, label=`Number of emVars`), vjust=1, color="#1b9e77") + 
	geom_text(aes(`Cell line`, `Fold enrichment`, label=round(-log10(`P-value`), 2)), vjust=0, color="#7570b3") + 
	geom_hline(yintercept=1, color="#eeac6f") + 
	theme_bw() + 
	ggtitle("Enrichment of emVars \n in ENCODE DHS regions") +
	theme(
		aspect.ratio=1.25, 
		panel.grid.major = element_blank(),
		panel.grid.minor = element_blank(), 
		plot.title = element_text(hjust = 0.5)
	)
ggsave("../../results/3-mpra_sanity_check_analyses/mpra_emVars_ENCODE_DHS_enrichment-K562-Jurkat.pdf", scale=0.6)


# enrichment in TF motif disruptions (HOCOMOCO)
mpra_emVars_K562_TF_motifs_hocomoco <- mpra_ccre_variants %>% 
	filter(!is.na(ENCODE_cCREs_summ_Summary)) %>% 
	mutate(MPRA_K562_emVars = !is.na(MPRA_K562_summ_Summary)) %>% 
	mutate(TF_motifs_hocomoco = !is.na(TF_motifs_hocomoco_summ_Summary)) %>% 
	dplyr::select(MPRA_K562_emVars, TF_motifs_hocomoco)
mpra_emVars_K562_TF_motifs_hocomoco_tab <- table(mpra_emVars_K562_TF_motifs_hocomoco)
mpra_emVars_K562_TF_motifs_hocomoco_fisher <- fisher.test(mpra_emVars_K562_TF_motifs_hocomoco_tab)
mpra_emVars_K562_TF_motifs_hocomoco_fisher_tb <- as_tibble(t(c(
	"P-value"=mpra_emVars_K562_TF_motifs_hocomoco_fisher$p.value, 
	"Fold enrichment"=mpra_emVars_K562_TF_motifs_hocomoco_fisher$estimate[["odds ratio"]], 
	"Number of emVars"=sum(mpra_emVars_K562_TF_motifs_hocomoco$MPRA_K562_emVars
)))) %>% mutate(`Cell line`="K562")

mpra_emVars_Jurkat_TF_motifs_hocomoco <- mpra_ccre_variants %>% 
	filter(!is.na(ENCODE_cCREs_summ_Summary)) %>% 
	mutate(MPRA_Jurkat_emVars = !is.na(MPRA_Jurkat_summ_Summary)) %>% 
	mutate(TF_motifs_hocomoco = !is.na(TF_motifs_hocomoco_summ_Summary)) %>% 
	dplyr::select(MPRA_Jurkat_emVars, TF_motifs_hocomoco)
mpra_emVars_Jurkat_TF_motifs_hocomoco_tab <- table(mpra_emVars_Jurkat_TF_motifs_hocomoco)
mpra_emVars_Jurkat_TF_motifs_hocomoco_fisher <- fisher.test(mpra_emVars_Jurkat_TF_motifs_hocomoco_tab)
mpra_emVars_Jurkat_TF_motifs_hocomoco_fisher_tb <- as_tibble(t(c(
	"P-value"=mpra_emVars_Jurkat_TF_motifs_hocomoco_fisher$p.value, 
	"Fold enrichment"=mpra_emVars_Jurkat_TF_motifs_hocomoco_fisher$estimate[["odds ratio"]], 
	"Number of emVars"=sum(mpra_emVars_Jurkat_TF_motifs_hocomoco$MPRA_Jurkat_emVars
)))) %>% mutate(`Cell line`="Jurkat")


mpra_emVars_TF_motifs_hocomoco_fisher_tb <- bind_rows(
	mpra_emVars_K562_TF_motifs_hocomoco_fisher_tb, 
	mpra_emVars_Jurkat_TF_motifs_hocomoco_fisher_tb
) %>% 
	mutate(`Cell line`=factor(`Cell line`, levels = c("K562", "Jurkat")))

ggplot(mpra_emVars_TF_motifs_hocomoco_fisher_tb) + 
	geom_bar(aes(`Cell line`, `Fold enrichment`), stat="identity", fill="#6baed6") + 
	geom_text(aes(`Cell line`, 0, label=`Number of emVars`), vjust=1, color="#1b9e77") + 
	geom_text(aes(`Cell line`, `Fold enrichment`, label=round(-log10(`P-value`), 2)), vjust=0, color="#7570b3") + 
	geom_hline(yintercept=1, color="#eeac6f") + 
	theme_bw() + 
	ggtitle("Enrichment of emVars with \n TF motif disruptions (HOCOMOCO)") +
	theme(
		aspect.ratio=1.25, 
		panel.grid.major = element_blank(),
		panel.grid.minor = element_blank(), 
		plot.title = element_text(hjust = 0.5)
	)
ggsave("../../results/3-mpra_sanity_check_analyses/mpra_emVars_TF_motifs_hocomoco_enrichment-K562-Jurkat.pdf", scale=0.6)


# enrichment in TF motif disruptions (JASPAR)
mpra_emVars_K562_TF_motifs_jaspar <- mpra_ccre_variants %>% 
	filter(!is.na(ENCODE_cCREs_summ_Summary)) %>% 
	mutate(MPRA_K562_emVars = !is.na(MPRA_K562_summ_Summary)) %>% 
	mutate(TF_motifs_jaspar = !is.na(TF_motifs_jaspar_summ_Summary)) %>% 
	dplyr::select(MPRA_K562_emVars, TF_motifs_jaspar)
mpra_emVars_K562_TF_motifs_jaspar_tab <- table(mpra_emVars_K562_TF_motifs_jaspar)
mpra_emVars_K562_TF_motifs_jaspar_fisher <- fisher.test(mpra_emVars_K562_TF_motifs_jaspar_tab)
mpra_emVars_K562_TF_motifs_jaspar_fisher_tb <- as_tibble(t(c(
	"P-value"=mpra_emVars_K562_TF_motifs_jaspar_fisher$p.value, 
	"Fold enrichment"=mpra_emVars_K562_TF_motifs_jaspar_fisher$estimate[["odds ratio"]], 
	"Number of emVars"=sum(mpra_emVars_K562_TF_motifs_jaspar$MPRA_K562_emVars
)))) %>% mutate(`Cell line`="K562")

mpra_emVars_Jurkat_TF_motifs_jaspar <- mpra_ccre_variants %>% 
	filter(!is.na(ENCODE_cCREs_summ_Summary)) %>% 
	mutate(MPRA_Jurkat_emVars = !is.na(MPRA_Jurkat_summ_Summary)) %>% 
	mutate(TF_motifs_jaspar = !is.na(TF_motifs_jaspar_summ_Summary)) %>% 
	dplyr::select(MPRA_Jurkat_emVars, TF_motifs_jaspar)
mpra_emVars_Jurkat_TF_motifs_jaspar_tab <- table(mpra_emVars_Jurkat_TF_motifs_jaspar)
mpra_emVars_Jurkat_TF_motifs_jaspar_fisher <- fisher.test(mpra_emVars_Jurkat_TF_motifs_jaspar_tab)
mpra_emVars_Jurkat_TF_motifs_jaspar_fisher_tb <- as_tibble(t(c(
	"P-value"=mpra_emVars_Jurkat_TF_motifs_jaspar_fisher$p.value, 
	"Fold enrichment"=mpra_emVars_Jurkat_TF_motifs_jaspar_fisher$estimate[["odds ratio"]], 
	"Number of emVars"=sum(mpra_emVars_Jurkat_TF_motifs_jaspar$MPRA_Jurkat_emVars
)))) %>% mutate(`Cell line`="Jurkat")


mpra_emVars_TF_motifs_jaspar_fisher_tb <- bind_rows(
	mpra_emVars_K562_TF_motifs_jaspar_fisher_tb, 
	mpra_emVars_Jurkat_TF_motifs_jaspar_fisher_tb
) %>% 
	mutate(`Cell line`=factor(`Cell line`, levels = c("K562", "Jurkat")))

ggplot(mpra_emVars_TF_motifs_jaspar_fisher_tb) + 
	geom_bar(aes(`Cell line`, `Fold enrichment`), stat="identity", fill="#6baed6") + 
	geom_text(aes(`Cell line`, 0, label=`Number of emVars`), vjust=1, color="#1b9e77") + 
	geom_text(aes(`Cell line`, `Fold enrichment`, label=round(-log10(`P-value`), 2)), vjust=0, color="#7570b3") + 
	geom_hline(yintercept=1, color="#eeac6f") + 
	theme_bw() + 
	ggtitle("Enrichment of emVars with \n TF motif disruptions (JASPAR)") +
	theme(
		aspect.ratio=1.25, 
		panel.grid.major = element_blank(),
		panel.grid.minor = element_blank(), 
		plot.title = element_text(hjust = 0.5)
	)
ggsave("../../results/3-mpra_sanity_check_analyses/mpra_emVars_TF_motifs_jaspar_enrichment-K562-Jurkat.pdf", scale=0.6)


# enrichment in DNase footprint overlaps
mpra_emVars_K562_DNase_footprints <- mpra_ccre_variants %>% 
	filter(!is.na(ENCODE_cCREs_summ_Summary)) %>% 
	mutate(MPRA_K562_emVars = !is.na(MPRA_K562_summ_Summary)) %>% 
	mutate(DNase_footprints = !is.na(DNase_footprints_summ_Summary)) %>% 
	dplyr::select(MPRA_K562_emVars, DNase_footprints)
mpra_emVars_K562_DNase_footprints_tab <- table(mpra_emVars_K562_DNase_footprints)
mpra_emVars_K562_DNase_footprints_fisher <- fisher.test(mpra_emVars_K562_DNase_footprints_tab)
mpra_emVars_K562_DNase_footprints_fisher_tb <- as_tibble(t(c(
	"P-value"=mpra_emVars_K562_DNase_footprints_fisher$p.value, 
	"Fold enrichment"=mpra_emVars_K562_DNase_footprints_fisher$estimate[["odds ratio"]], 
	"Number of emVars"=sum(mpra_emVars_K562_DNase_footprints$MPRA_K562_emVars
)))) %>% mutate(`Cell line`="K562")

mpra_emVars_Jurkat_DNase_footprints <- mpra_ccre_variants %>% 
	filter(!is.na(ENCODE_cCREs_summ_Summary)) %>% 
	mutate(MPRA_Jurkat_emVars = !is.na(MPRA_Jurkat_summ_Summary)) %>% 
	mutate(DNase_footprints = !is.na(DNase_footprints_summ_Summary)) %>% 
	dplyr::select(MPRA_Jurkat_emVars, DNase_footprints)
mpra_emVars_Jurkat_DNase_footprints_tab <- table(mpra_emVars_Jurkat_DNase_footprints)
mpra_emVars_Jurkat_DNase_footprints_fisher <- fisher.test(mpra_emVars_Jurkat_DNase_footprints_tab)
mpra_emVars_Jurkat_DNase_footprints_fisher_tb <- as_tibble(t(c(
	"P-value"=mpra_emVars_Jurkat_DNase_footprints_fisher$p.value, 
	"Fold enrichment"=mpra_emVars_Jurkat_DNase_footprints_fisher$estimate[["odds ratio"]], 
	"Number of emVars"=sum(mpra_emVars_Jurkat_DNase_footprints$MPRA_Jurkat_emVars
)))) %>% mutate(`Cell line`="Jurkat")


mpra_emVars_DNase_footprints_fisher_tb <- bind_rows(
	mpra_emVars_K562_DNase_footprints_fisher_tb, 
	mpra_emVars_Jurkat_DNase_footprints_fisher_tb
) %>% 
	mutate(`Cell line`=factor(`Cell line`, levels = c("K562", "Jurkat")))

ggplot(mpra_emVars_DNase_footprints_fisher_tb) + 
	geom_bar(aes(`Cell line`, `Fold enrichment`), stat="identity", fill="#6baed6") + 
	geom_text(aes(`Cell line`, 0, label=`Number of emVars`), vjust=1, color="#1b9e77") + 
	geom_text(aes(`Cell line`, `Fold enrichment`, label=round(-log10(`P-value`), 2)), vjust=0, color="#7570b3") + 
	geom_hline(yintercept=1, color="#eeac6f") + 
	theme_bw() + 
	ggtitle("Enrichment of emVars \n in DNase footprints") +
	theme(
		aspect.ratio=1.25, 
		panel.grid.major = element_blank(),
		panel.grid.minor = element_blank(), 
		plot.title = element_text(hjust = 0.5)
	)
ggsave("../../results/3-mpra_sanity_check_analyses/mpra_emVars_DNase_footprints_enrichment-K562-Jurkat.pdf", scale=0.6)


# enrichment in ReMap ChIP-seq peaks
mpra_emVars_K562_ReMap_TFBS <- mpra_ccre_variants %>% 
	filter(!is.na(ENCODE_cCREs_summ_Summary)) %>% 
	mutate(MPRA_K562_emVars = !is.na(MPRA_K562_summ_Summary)) %>% 
	mutate(ReMap_TFBS_K562 = !is.na(ReMap_TFBS_K562_summ_Summary)) %>% 
	dplyr::select(MPRA_K562_emVars, ReMap_TFBS_K562)
mpra_emVars_K562_ReMap_TFBS_tab <- table(mpra_emVars_K562_ReMap_TFBS)
mpra_emVars_K562_ReMap_TFBS_fisher <- fisher.test(mpra_emVars_K562_ReMap_TFBS_tab)
mpra_emVars_K562_ReMap_TFBS_fisher_tb <- as_tibble(t(c(
	"P-value"=mpra_emVars_K562_ReMap_TFBS_fisher$p.value, 
	"Fold enrichment"=mpra_emVars_K562_ReMap_TFBS_fisher$estimate[["odds ratio"]], 
	"Number of emVars"=sum(mpra_emVars_K562_ReMap_TFBS$MPRA_K562_emVars
)))) %>% mutate(`Cell line`="K562")

mpra_emVars_Jurkat_ReMap_TFBS <- mpra_ccre_variants %>% 
	filter(!is.na(ENCODE_cCREs_summ_Summary)) %>% 
	mutate(MPRA_Jurkat_emVars = !is.na(MPRA_Jurkat_summ_Summary)) %>% 
	mutate(ReMap_TFBS_Jurkat = !is.na(ReMap_TFBS_Jurkat_summ_Summary)) %>% 
	dplyr::select(MPRA_Jurkat_emVars, ReMap_TFBS_Jurkat)
mpra_emVars_Jurkat_ReMap_TFBS_tab <- table(mpra_emVars_Jurkat_ReMap_TFBS)
mpra_emVars_Jurkat_ReMap_TFBS_fisher <- fisher.test(mpra_emVars_Jurkat_ReMap_TFBS_tab)
mpra_emVars_Jurkat_ReMap_TFBS_fisher_tb <- as_tibble(t(c(
	"P-value"=mpra_emVars_Jurkat_ReMap_TFBS_fisher$p.value, 
	"Fold enrichment"=mpra_emVars_Jurkat_ReMap_TFBS_fisher$estimate[["odds ratio"]], 
	"Number of emVars"=sum(mpra_emVars_Jurkat_ReMap_TFBS$MPRA_Jurkat_emVars
)))) %>% mutate(`Cell line`="Jurkat")


mpra_emVars_ReMap_TFBS_fisher_tb <- bind_rows(
	mpra_emVars_K562_ReMap_TFBS_fisher_tb, 
	mpra_emVars_Jurkat_ReMap_TFBS_fisher_tb
) %>% 
	mutate(`Cell line`=factor(`Cell line`, levels = c("K562", "Jurkat")))

ggplot(mpra_emVars_ReMap_TFBS_fisher_tb) + 
	geom_bar(aes(`Cell line`, `Fold enrichment`), stat="identity", fill="#6baed6") + 
	geom_text(aes(`Cell line`, 0, label=`Number of emVars`), vjust=1, color="#1b9e77") + 
	geom_text(aes(`Cell line`, `Fold enrichment`, label=round(-log10(`P-value`), 2)), vjust=0, color="#7570b3") + 
	geom_hline(yintercept=1, color="#eeac6f") + 
	theme_bw() + 
	ggtitle("Enrichment of emVars \n in ReMap TF ChIP-seq peaks") +
	theme(
		aspect.ratio=1.25, 
		panel.grid.major = element_blank(),
		panel.grid.minor = element_blank(), 
		plot.title = element_text(hjust = 0.5)
	)
ggsave("../../results/3-mpra_sanity_check_analyses/mpra_emVars_ReMap_TFBS_enrichment-K562-Jurkat.pdf", scale=0.6)


# enrichmet in DHS index peaks
mpra_emVars_K562_DHS_index <- mpra_ccre_variants %>% 
	filter(!is.na(ENCODE_cCREs_summ_Summary)) %>% 
	mutate(MPRA_K562_emVars = !is.na(MPRA_K562_summ_Summary)) %>% 
	mutate(DHS_index = !is.na(DHS_index_vocabulary_summ_Summary)) %>% 
	dplyr::select(MPRA_K562_emVars, DHS_index)
mpra_emVars_K562_DHS_index_vocabulary_tab <- table(mpra_emVars_K562_DHS_index)
mpra_emVars_K562_DHS_index_vocabulary_fisher <- fisher.test(mpra_emVars_K562_DHS_index_vocabulary_tab)
mpra_emVars_K562_DHS_index_vocabulary_fisher_tb <- as_tibble(t(c(
	"P-value"=mpra_emVars_K562_DHS_index_vocabulary_fisher$p.value, 
	"Fold enrichment"=mpra_emVars_K562_DHS_index_vocabulary_fisher$estimate[["odds ratio"]], 
	"Number of emVars"=sum(mpra_emVars_K562_DHS_index$MPRA_K562_emVars
)))) %>% mutate(`Cell line`="K562")

mpra_emVars_Jurkat_DHS_index <- mpra_ccre_variants %>% 
	filter(!is.na(ENCODE_cCREs_summ_Summary)) %>% 
	mutate(MPRA_Jurkat_emVars = !is.na(MPRA_Jurkat_summ_Summary)) %>% 
	mutate(DHS_index = !is.na(DHS_index_vocabulary_summ_Summary)) %>% 
	dplyr::select(MPRA_Jurkat_emVars, DHS_index)
mpra_emVars_Jurkat_DHS_index_vocabulary_tab <- table(mpra_emVars_Jurkat_DHS_index)
mpra_emVars_Jurkat_DHS_index_vocabulary_fisher <- fisher.test(mpra_emVars_Jurkat_DHS_index_vocabulary_tab)
mpra_emVars_Jurkat_DHS_index_vocabulary_fisher_tb <- as_tibble(t(c(
	"P-value"=mpra_emVars_Jurkat_DHS_index_vocabulary_fisher$p.value, 
	"Fold enrichment"=mpra_emVars_Jurkat_DHS_index_vocabulary_fisher$estimate[["odds ratio"]], 
	"Number of emVars"=sum(mpra_emVars_Jurkat_DHS_index$MPRA_Jurkat_emVars
)))) %>% mutate(`Cell line`="Jurkat")


mpra_emVars_DHS_index_vocabulary_fisher_tb <- bind_rows(
	mpra_emVars_K562_DHS_index_vocabulary_fisher_tb, 
	mpra_emVars_Jurkat_DHS_index_vocabulary_fisher_tb
) %>% 
	mutate(`Cell line`=factor(`Cell line`, levels = c("K562", "Jurkat")))

ggplot(mpra_emVars_DHS_index_vocabulary_fisher_tb) + 
	geom_bar(aes(`Cell line`, `Fold enrichment`), stat="identity", fill="#6baed6") + 
	geom_text(aes(`Cell line`, 0, label=`Number of emVars`), vjust=1, color="#1b9e77") + 
	geom_text(aes(`Cell line`, `Fold enrichment`, label=round(-log10(`P-value`), 2)), vjust=0, color="#7570b3") + 
	geom_hline(yintercept=1, color="#eeac6f") + 
	theme_bw() + 
	ggtitle("Enrichment of emVars in \n DHS index peaks") +
	theme(
		aspect.ratio=1.25, 
		panel.grid.major = element_blank(),
		panel.grid.minor = element_blank(), 
		plot.title = element_text(hjust = 0.5)
	)
ggsave("../../results/3-mpra_sanity_check_analyses/mpra_emVars_DHS_index_vocabulary_enrichment-K562-Jurkat.pdf", scale=0.6)


# enrichment in Roadmap Epigenomics merged
mpra_emVars_K562_Roadmap_Epigenomics_merged <- mpra_ccre_variants %>% 
	filter(!is.na(ENCODE_cCREs_summ_Summary)) %>% 
	mutate(MPRA_K562_emVars = !is.na(MPRA_K562_summ_Summary)) %>% 
	mutate(Roadmap_Epigenomics_merged = !is.na(Roadmap_Epigenomics_merged_summ_Summary)) %>% 
	dplyr::select(MPRA_K562_emVars, Roadmap_Epigenomics_merged)
mpra_emVars_K562_Roadmap_Epigenomics_merged_tab <- table(mpra_emVars_K562_Roadmap_Epigenomics_merged)
mpra_emVars_K562_Roadmap_Epigenomics_merged_fisher <- fisher.test(mpra_emVars_K562_Roadmap_Epigenomics_merged_tab)
mpra_emVars_K562_Roadmap_Epigenomics_merged_fisher_tb <- as_tibble(t(c(
	"P-value"=mpra_emVars_K562_Roadmap_Epigenomics_merged_fisher$p.value, 
	"Fold enrichment"=mpra_emVars_K562_Roadmap_Epigenomics_merged_fisher$estimate[["odds ratio"]], 
	"Number of emVars"=sum(mpra_emVars_K562_Roadmap_Epigenomics_merged$MPRA_K562_emVars
)))) %>% mutate(`Cell line`="K562")

mpra_emVars_Jurkat_Roadmap_Epigenomics_merged <- mpra_ccre_variants %>% 
	filter(!is.na(ENCODE_cCREs_summ_Summary)) %>% 
	mutate(MPRA_Jurkat_emVars = !is.na(MPRA_Jurkat_summ_Summary)) %>% 
	mutate(Roadmap_Epigenomics_merged = !is.na(Roadmap_Epigenomics_merged_summ_Summary)) %>% 
	dplyr::select(MPRA_Jurkat_emVars, Roadmap_Epigenomics_merged)
mpra_emVars_Jurkat_Roadmap_Epigenomics_merged_tab <- table(mpra_emVars_Jurkat_Roadmap_Epigenomics_merged)
mpra_emVars_Jurkat_Roadmap_Epigenomics_merged_fisher <- fisher.test(mpra_emVars_Jurkat_Roadmap_Epigenomics_merged_tab)
mpra_emVars_Jurkat_Roadmap_Epigenomics_merged_fisher_tb <- as_tibble(t(c(
	"P-value"=mpra_emVars_Jurkat_Roadmap_Epigenomics_merged_fisher$p.value, 
	"Fold enrichment"=mpra_emVars_Jurkat_Roadmap_Epigenomics_merged_fisher$estimate[["odds ratio"]], 
	"Number of emVars"=sum(mpra_emVars_Jurkat_Roadmap_Epigenomics_merged$MPRA_Jurkat_emVars
)))) %>% mutate(`Cell line`="Jurkat")


mpra_emVars_Roadmap_Epigenomics_merged_fisher_tb <- bind_rows(
	mpra_emVars_K562_Roadmap_Epigenomics_merged_fisher_tb, 
	mpra_emVars_Jurkat_Roadmap_Epigenomics_merged_fisher_tb
) %>% 
	mutate(`Cell line`=factor(`Cell line`, levels = c("K562", "Jurkat")))

ggplot(mpra_emVars_Roadmap_Epigenomics_merged_fisher_tb) + 
	geom_bar(aes(`Cell line`, `Fold enrichment`), stat="identity", fill="#6baed6") + 
	geom_text(aes(`Cell line`, 0, label=`Number of emVars`), vjust=1, color="#1b9e77") + 
	geom_text(aes(`Cell line`, `Fold enrichment`, label=round(-log10(`P-value`), 2)), vjust=0, color="#7570b3") + 
	geom_hline(yintercept=1, color="#eeac6f") + 
	theme_bw() + 
	ggtitle("Enrichment of emVars in \n Roadmap Epigenomics merged") +
	theme(
		aspect.ratio=1.25, 
		panel.grid.major = element_blank(),
		panel.grid.minor = element_blank(), 
		plot.title = element_text(hjust = 0.5)
	)
ggsave("../../results/3-mpra_sanity_check_analyses/mpra_emVars_Roadmap_Epigenomics_merged_enrichment-K562-Jurkat.pdf", scale=0.6)


# enrichment in Roadmap Epigenomics promoters
mpra_emVars_K562_Roadmap_Epigenomics_promoters <- mpra_ccre_variants %>% 
	filter(!is.na(ENCODE_cCREs_summ_Summary)) %>% 
	mutate(MPRA_K562_emVars = !is.na(MPRA_K562_summ_Summary)) %>% 
	mutate(Roadmap_Epigenomics_promoters = !is.na(Roadmap_Epigenomics_promoters_summ_Summary)) %>% 
	dplyr::select(MPRA_K562_emVars, Roadmap_Epigenomics_promoters)
mpra_emVars_K562_Roadmap_Epigenomics_promoters_tab <- table(mpra_emVars_K562_Roadmap_Epigenomics_promoters)
mpra_emVars_K562_Roadmap_Epigenomics_promoters_fisher <- fisher.test(mpra_emVars_K562_Roadmap_Epigenomics_promoters_tab)
mpra_emVars_K562_Roadmap_Epigenomics_promoters_fisher_tb <- as_tibble(t(c(
	"P-value"=mpra_emVars_K562_Roadmap_Epigenomics_promoters_fisher$p.value, 
	"Fold enrichment"=mpra_emVars_K562_Roadmap_Epigenomics_promoters_fisher$estimate[["odds ratio"]], 
	"Number of emVars"=sum(mpra_emVars_K562_Roadmap_Epigenomics_promoters$MPRA_K562_emVars
)))) %>% mutate(`Cell line`="K562")

mpra_emVars_Jurkat_Roadmap_Epigenomics_promoters <- mpra_ccre_variants %>% 
	filter(!is.na(ENCODE_cCREs_summ_Summary)) %>% 
	mutate(MPRA_Jurkat_emVars = !is.na(MPRA_Jurkat_summ_Summary)) %>% 
	mutate(Roadmap_Epigenomics_promoters = !is.na(Roadmap_Epigenomics_promoters_summ_Summary)) %>% 
	dplyr::select(MPRA_Jurkat_emVars, Roadmap_Epigenomics_promoters)
mpra_emVars_Jurkat_Roadmap_Epigenomics_promoters_tab <- table(mpra_emVars_Jurkat_Roadmap_Epigenomics_promoters)
mpra_emVars_Jurkat_Roadmap_Epigenomics_promoters_fisher <- fisher.test(mpra_emVars_Jurkat_Roadmap_Epigenomics_promoters_tab)
mpra_emVars_Jurkat_Roadmap_Epigenomics_promoters_fisher_tb <- as_tibble(t(c(
	"P-value"=mpra_emVars_Jurkat_Roadmap_Epigenomics_promoters_fisher$p.value, 
	"Fold enrichment"=mpra_emVars_Jurkat_Roadmap_Epigenomics_promoters_fisher$estimate[["odds ratio"]], 
	"Number of emVars"=sum(mpra_emVars_Jurkat_Roadmap_Epigenomics_promoters$MPRA_Jurkat_emVars
)))) %>% mutate(`Cell line`="Jurkat")


mpra_emVars_Roadmap_Epigenomics_promoters_fisher_tb <- bind_rows(
	mpra_emVars_K562_Roadmap_Epigenomics_promoters_fisher_tb, 
	mpra_emVars_Jurkat_Roadmap_Epigenomics_promoters_fisher_tb
) %>% 
	mutate(`Cell line`=factor(`Cell line`, levels = c("K562", "Jurkat")))

ggplot(mpra_emVars_Roadmap_Epigenomics_promoters_fisher_tb) + 
	geom_bar(aes(`Cell line`, `Fold enrichment`), stat="identity", fill="#6baed6") + 
	geom_text(aes(`Cell line`, 0, label=`Number of emVars`), vjust=1, color="#1b9e77") + 
	geom_text(aes(`Cell line`, `Fold enrichment`, label=round(-log10(`P-value`), 2)), vjust=0, color="#7570b3") + 
	geom_hline(yintercept=1, color="#eeac6f") + 
	theme_bw() + 
	ggtitle("Enrichment of emVars in \n Roadmap Epigenomics promoters") +
	theme(
		aspect.ratio=1.25, 
		panel.grid.major = element_blank(),
		panel.grid.minor = element_blank(), 
		plot.title = element_text(hjust = 0.5)
	)
ggsave("../../results/3-mpra_sanity_check_analyses/mpra_emVars_Roadmap_Epigenomics_promoters_enrichment-K562-Jurkat.pdf", scale=0.6)


# enrichment in Roadmap Epigenomics enhancers
mpra_emVars_K562_Roadmap_Epigenomics_enhancers <- mpra_ccre_variants %>% 
	filter(!is.na(ENCODE_cCREs_summ_Summary)) %>% 
	mutate(MPRA_K562_emVars = !is.na(MPRA_K562_summ_Summary)) %>% 
	mutate(Roadmap_Epigenomics_enhancers = !is.na(Roadmap_Epigenomics_enhancers_summ_Summary)) %>% 
	dplyr::select(MPRA_K562_emVars, Roadmap_Epigenomics_enhancers)
mpra_emVars_K562_Roadmap_Epigenomics_enhancers_tab <- table(mpra_emVars_K562_Roadmap_Epigenomics_enhancers)
mpra_emVars_K562_Roadmap_Epigenomics_enhancers_fisher <- fisher.test(mpra_emVars_K562_Roadmap_Epigenomics_enhancers_tab)
mpra_emVars_K562_Roadmap_Epigenomics_enhancers_fisher_tb <- as_tibble(t(c(
	"P-value"=mpra_emVars_K562_Roadmap_Epigenomics_enhancers_fisher$p.value, 
	"Fold enrichment"=mpra_emVars_K562_Roadmap_Epigenomics_enhancers_fisher$estimate[["odds ratio"]], 
	"Number of emVars"=sum(mpra_emVars_K562_Roadmap_Epigenomics_enhancers$MPRA_K562_emVars
)))) %>% mutate(`Cell line`="K562")

mpra_emVars_Jurkat_Roadmap_Epigenomics_enhancers <- mpra_ccre_variants %>% 
	filter(!is.na(ENCODE_cCREs_summ_Summary)) %>% 
	mutate(MPRA_Jurkat_emVars = !is.na(MPRA_Jurkat_summ_Summary)) %>% 
	mutate(Roadmap_Epigenomics_enhancers = !is.na(Roadmap_Epigenomics_enhancers_summ_Summary)) %>% 
	dplyr::select(MPRA_Jurkat_emVars, Roadmap_Epigenomics_enhancers)
mpra_emVars_Jurkat_Roadmap_Epigenomics_enhancers_tab <- table(mpra_emVars_Jurkat_Roadmap_Epigenomics_enhancers)
mpra_emVars_Jurkat_Roadmap_Epigenomics_enhancers_fisher <- fisher.test(mpra_emVars_Jurkat_Roadmap_Epigenomics_enhancers_tab)
mpra_emVars_Jurkat_Roadmap_Epigenomics_enhancers_fisher_tb <- as_tibble(t(c(
	"P-value"=mpra_emVars_Jurkat_Roadmap_Epigenomics_enhancers_fisher$p.value, 
	"Fold enrichment"=mpra_emVars_Jurkat_Roadmap_Epigenomics_enhancers_fisher$estimate[["odds ratio"]], 
	"Number of emVars"=sum(mpra_emVars_Jurkat_Roadmap_Epigenomics_enhancers$MPRA_Jurkat_emVars
)))) %>% mutate(`Cell line`="Jurkat")


mpra_emVars_Roadmap_Epigenomics_enhancers_fisher_tb <- bind_rows(
	mpra_emVars_K562_Roadmap_Epigenomics_enhancers_fisher_tb, 
	mpra_emVars_Jurkat_Roadmap_Epigenomics_enhancers_fisher_tb
) %>% 
	mutate(`Cell line`=factor(`Cell line`, levels = c("K562", "Jurkat")))

ggplot(mpra_emVars_Roadmap_Epigenomics_enhancers_fisher_tb) + 
	geom_bar(aes(`Cell line`, `Fold enrichment`), stat="identity", fill="#6baed6") + 
	geom_text(aes(`Cell line`, 0, label=`Number of emVars`), vjust=1, color="#1b9e77") + 
	geom_text(aes(`Cell line`, `Fold enrichment`, label=round(-log10(`P-value`), 2)), vjust=0, color="#7570b3") + 
	geom_hline(yintercept=1, color="#eeac6f") + 
	theme_bw() + 
	ggtitle("Enrichment of emVars in \n Roadmap Epigenomics enhancers") +
	theme(
		aspect.ratio=1.25, 
		panel.grid.major = element_blank(),
		panel.grid.minor = element_blank(), 
		plot.title = element_text(hjust = 0.5)
	)
ggsave("../../results/3-mpra_sanity_check_analyses/mpra_emVars_Roadmap_Epigenomics_enhancers_enrichment-K562-Jurkat.pdf", scale=0.6)


# enrichment in Roadmap Epigenomics dyadic
mpra_emVars_K562_Roadmap_Epigenomics_dyadic <- mpra_ccre_variants %>% 
	filter(!is.na(ENCODE_cCREs_summ_Summary)) %>% 
	mutate(MPRA_K562_emVars = !is.na(MPRA_K562_summ_Summary)) %>% 
	mutate(Roadmap_Epigenomics_dyadic = !is.na(Roadmap_Epigenomics_dyadic_summ_Summary)) %>% 
	dplyr::select(MPRA_K562_emVars, Roadmap_Epigenomics_dyadic)
mpra_emVars_K562_Roadmap_Epigenomics_dyadic_tab <- table(mpra_emVars_K562_Roadmap_Epigenomics_dyadic)
mpra_emVars_K562_Roadmap_Epigenomics_dyadic_fisher <- fisher.test(mpra_emVars_K562_Roadmap_Epigenomics_dyadic_tab)
mpra_emVars_K562_Roadmap_Epigenomics_dyadic_fisher_tb <- as_tibble(t(c(
	"P-value"=mpra_emVars_K562_Roadmap_Epigenomics_dyadic_fisher$p.value, 
	"Fold enrichment"=mpra_emVars_K562_Roadmap_Epigenomics_dyadic_fisher$estimate[["odds ratio"]], 
	"Number of emVars"=sum(mpra_emVars_K562_Roadmap_Epigenomics_dyadic$MPRA_K562_emVars
)))) %>% mutate(`Cell line`="K562")

mpra_emVars_Jurkat_Roadmap_Epigenomics_dyadic <- mpra_ccre_variants %>% 
	filter(!is.na(ENCODE_cCREs_summ_Summary)) %>% 
	mutate(MPRA_Jurkat_emVars = !is.na(MPRA_Jurkat_summ_Summary)) %>% 
	mutate(Roadmap_Epigenomics_dyadic = !is.na(Roadmap_Epigenomics_dyadic_summ_Summary)) %>% 
	dplyr::select(MPRA_Jurkat_emVars, Roadmap_Epigenomics_dyadic)
mpra_emVars_Jurkat_Roadmap_Epigenomics_dyadic_tab <- table(mpra_emVars_Jurkat_Roadmap_Epigenomics_dyadic)
mpra_emVars_Jurkat_Roadmap_Epigenomics_dyadic_fisher <- fisher.test(mpra_emVars_Jurkat_Roadmap_Epigenomics_dyadic_tab)
mpra_emVars_Jurkat_Roadmap_Epigenomics_dyadic_fisher_tb <- as_tibble(t(c(
	"P-value"=mpra_emVars_Jurkat_Roadmap_Epigenomics_dyadic_fisher$p.value, 
	"Fold enrichment"=mpra_emVars_Jurkat_Roadmap_Epigenomics_dyadic_fisher$estimate[["odds ratio"]], 
	"Number of emVars"=sum(mpra_emVars_Jurkat_Roadmap_Epigenomics_dyadic$MPRA_Jurkat_emVars
)))) %>% mutate(`Cell line`="Jurkat")


mpra_emVars_Roadmap_Epigenomics_dyadic_fisher_tb <- bind_rows(
	mpra_emVars_K562_Roadmap_Epigenomics_dyadic_fisher_tb, 
	mpra_emVars_Jurkat_Roadmap_Epigenomics_dyadic_fisher_tb
) %>% 
	mutate(`Cell line`=factor(`Cell line`, levels = c("K562", "Jurkat")))

ggplot(mpra_emVars_Roadmap_Epigenomics_dyadic_fisher_tb) + 
	geom_bar(aes(`Cell line`, `Fold enrichment`), stat="identity", fill="#6baed6") + 
	geom_text(aes(`Cell line`, 0, label=`Number of emVars`), vjust=1, color="#1b9e77") + 
	geom_text(aes(`Cell line`, `Fold enrichment`, label=round(-log10(`P-value`), 2)), vjust=0, color="#7570b3") + 
	geom_hline(yintercept=1, color="#eeac6f") + 
	theme_bw() + 
	ggtitle("Enrichment of emVars in \n Roadmap Epigenomics dyadic") +
	theme(
		aspect.ratio=1.25, 
		panel.grid.major = element_blank(),
		panel.grid.minor = element_blank(), 
		plot.title = element_text(hjust = 0.5)
	)
ggsave("../../results/3-mpra_sanity_check_analyses/mpra_emVars_Roadmap_Epigenomics_dyadic_enrichment-K562-Jurkat.pdf", scale=0.6)


# enrichment in ABC gene-enhancer links
mpra_emVars_K562_ABC_links <- mpra_ccre_variants %>% 
	filter(!is.na(ENCODE_cCREs_summ_Summary)) %>% 
	mutate(MPRA_K562_emVars = !is.na(MPRA_K562_summ_Summary)) %>% 
	mutate(ABC_links = !is.na(ABC_summ_Summary)) %>% 
	dplyr::select(MPRA_K562_emVars, ABC_links)
mpra_emVars_K562_ABC_links_tab <- table(mpra_emVars_K562_ABC_links)
mpra_emVars_K562_ABC_links_fisher <- fisher.test(mpra_emVars_K562_ABC_links_tab)
mpra_emVars_K562_ABC_links_fisher_tb <- as_tibble(t(c(
	"P-value"=mpra_emVars_K562_ABC_links_fisher$p.value, 
	"Fold enrichment"=mpra_emVars_K562_ABC_links_fisher$estimate[["odds ratio"]], 
	"Number of emVars"=sum(mpra_emVars_K562_ABC_links$MPRA_K562_emVars
)))) %>% mutate(`Cell line`="K562")

mpra_emVars_Jurkat_ABC_links <- mpra_ccre_variants %>% 
	filter(!is.na(ENCODE_cCREs_summ_Summary)) %>% 
	mutate(MPRA_Jurkat_emVars = !is.na(MPRA_Jurkat_summ_Summary)) %>% 
	mutate(ABC_links = !is.na(ABC_summ_Summary)) %>% 
	dplyr::select(MPRA_Jurkat_emVars, ABC_links)
mpra_emVars_Jurkat_ABC_links_tab <- table(mpra_emVars_Jurkat_ABC_links)
mpra_emVars_Jurkat_ABC_links_fisher <- fisher.test(mpra_emVars_Jurkat_ABC_links_tab)
mpra_emVars_Jurkat_ABC_links_fisher_tb <- as_tibble(t(c(
	"P-value"=mpra_emVars_Jurkat_ABC_links_fisher$p.value, 
	"Fold enrichment"=mpra_emVars_Jurkat_ABC_links_fisher$estimate[["odds ratio"]], 
	"Number of emVars"=sum(mpra_emVars_Jurkat_ABC_links$MPRA_Jurkat_emVars
)))) %>% mutate(`Cell line`="Jurkat")


mpra_emVars_ABC_links_fisher_tb <- bind_rows(
	mpra_emVars_K562_ABC_links_fisher_tb, 
	mpra_emVars_Jurkat_ABC_links_fisher_tb
) %>% 
	mutate(`Cell line`=factor(`Cell line`, levels = c("K562", "Jurkat")))

ggplot(mpra_emVars_ABC_links_fisher_tb) + 
	geom_bar(aes(`Cell line`, `Fold enrichment`), stat="identity", fill="#6baed6") + 
	geom_text(aes(`Cell line`, 0, label=`Number of emVars`), vjust=1, color="#1b9e77") + 
	geom_text(aes(`Cell line`, `Fold enrichment`, label=round(-log10(`P-value`), 2)), vjust=0, color="#7570b3") + 
	geom_hline(yintercept=1, color="#eeac6f") + 
	theme_bw() + 
	ggtitle("Enrichment of emVars with \n ABC gene-enhancer links") +
	theme(
		aspect.ratio=1.25, 
		panel.grid.major = element_blank(),
		panel.grid.minor = element_blank(), 
		plot.title = element_text(hjust = 0.5)
	)
ggsave("../../results/3-mpra_sanity_check_analyses/mpra_emVars_ABC_links_enrichment-K562-Jurkat.pdf", scale=0.6)


# enrichment in EpiMap gene-enhancer links
mpra_emVars_K562_EpiMap_links <- mpra_ccre_variants %>% 
	filter(!is.na(ENCODE_cCREs_summ_Summary)) %>% 
	mutate(MPRA_K562_emVars = !is.na(MPRA_K562_summ_Summary)) %>% 
	mutate(EpiMap_links = !is.na(EpiMap_summ_Summary)) %>% 
	dplyr::select(MPRA_K562_emVars, EpiMap_links)
mpra_emVars_K562_EpiMap_links_tab <- table(mpra_emVars_K562_EpiMap_links)
mpra_emVars_K562_EpiMap_links_fisher <- fisher.test(mpra_emVars_K562_EpiMap_links_tab)
mpra_emVars_K562_EpiMap_links_fisher_tb <- as_tibble(t(c(
	"P-value"=mpra_emVars_K562_EpiMap_links_fisher$p.value, 
	"Fold enrichment"=mpra_emVars_K562_EpiMap_links_fisher$estimate[["odds ratio"]], 
	"Number of emVars"=sum(mpra_emVars_K562_EpiMap_links$MPRA_K562_emVars
)))) %>% mutate(`Cell line`="K562")

mpra_emVars_Jurkat_EpiMap_links <- mpra_ccre_variants %>% 
	filter(!is.na(ENCODE_cCREs_summ_Summary)) %>% 
	mutate(MPRA_Jurkat_emVars = !is.na(MPRA_Jurkat_summ_Summary)) %>% 
	mutate(EpiMap_links = !is.na(EpiMap_summ_Summary)) %>% 
	dplyr::select(MPRA_Jurkat_emVars, EpiMap_links)
mpra_emVars_Jurkat_EpiMap_links_tab <- table(mpra_emVars_Jurkat_EpiMap_links)
mpra_emVars_Jurkat_EpiMap_links_fisher <- fisher.test(mpra_emVars_Jurkat_EpiMap_links_tab)
mpra_emVars_Jurkat_EpiMap_links_fisher_tb <- as_tibble(t(c(
	"P-value"=mpra_emVars_Jurkat_EpiMap_links_fisher$p.value, 
	"Fold enrichment"=mpra_emVars_Jurkat_EpiMap_links_fisher$estimate[["odds ratio"]], 
	"Number of emVars"=sum(mpra_emVars_Jurkat_EpiMap_links$MPRA_Jurkat_emVars
)))) %>% mutate(`Cell line`="Jurkat")


mpra_emVars_EpiMap_links_fisher_tb <- bind_rows(
	mpra_emVars_K562_EpiMap_links_fisher_tb, 
	mpra_emVars_Jurkat_EpiMap_links_fisher_tb
) %>% 
	mutate(`Cell line`=factor(`Cell line`, levels = c("K562", "Jurkat")))

ggplot(mpra_emVars_EpiMap_links_fisher_tb) + 
	geom_bar(aes(`Cell line`, `Fold enrichment`), stat="identity", fill="#6baed6") + 
	geom_text(aes(`Cell line`, 0, label=`Number of emVars`), vjust=1, color="#1b9e77") + 
	geom_text(aes(`Cell line`, `Fold enrichment`, label=round(-log10(`P-value`), 2)), vjust=0, color="#7570b3") + 
	geom_hline(yintercept=1, color="#eeac6f") + 
	theme_bw() + 
	ggtitle("Enrichment of emVars with \n EpiMap gene-enhancer links") +
	theme(
		aspect.ratio=1.25, 
		panel.grid.major = element_blank(),
		panel.grid.minor = element_blank(), 
		plot.title = element_text(hjust = 0.5)
	)
ggsave("../../results/3-mpra_sanity_check_analyses/mpra_emVars_EpiMap_links_enrichment-K562-Jurkat.pdf", scale=0.6)


# enrichment in Roadmap Epigenomics gene-enhancer links
mpra_emVars_K562_Roadmap_Epigenomics_links <- mpra_ccre_variants %>% 
	filter(!is.na(ENCODE_cCREs_summ_Summary)) %>% 
	mutate(MPRA_K562_emVars = !is.na(MPRA_K562_summ_Summary)) %>% 
	mutate(Roadmap_Epigenomics_links = !is.na(Roadmap_Epigenomics_summ_Summary)) %>% 
	dplyr::select(MPRA_K562_emVars, Roadmap_Epigenomics_links)
mpra_emVars_K562_Roadmap_Epigenomics_links_tab <- table(mpra_emVars_K562_Roadmap_Epigenomics_links)
mpra_emVars_K562_Roadmap_Epigenomics_links_fisher <- fisher.test(mpra_emVars_K562_Roadmap_Epigenomics_links_tab)
mpra_emVars_K562_Roadmap_Epigenomics_links_fisher_tb <- as_tibble(t(c(
	"P-value"=mpra_emVars_K562_Roadmap_Epigenomics_links_fisher$p.value, 
	"Fold enrichment"=mpra_emVars_K562_Roadmap_Epigenomics_links_fisher$estimate[["odds ratio"]], 
	"Number of emVars"=sum(mpra_emVars_K562_Roadmap_Epigenomics_links$MPRA_K562_emVars
)))) %>% mutate(`Cell line`="K562")

mpra_emVars_Jurkat_Roadmap_Epigenomics_links <- mpra_ccre_variants %>% 
	filter(!is.na(ENCODE_cCREs_summ_Summary)) %>% 
	mutate(MPRA_Jurkat_emVars = !is.na(MPRA_Jurkat_summ_Summary)) %>% 
	mutate(Roadmap_Epigenomics_links = !is.na(Roadmap_Epigenomics_summ_Summary)) %>% 
	dplyr::select(MPRA_Jurkat_emVars, Roadmap_Epigenomics_links)
mpra_emVars_Jurkat_Roadmap_Epigenomics_links_tab <- table(mpra_emVars_Jurkat_Roadmap_Epigenomics_links)
mpra_emVars_Jurkat_Roadmap_Epigenomics_links_fisher <- fisher.test(mpra_emVars_Jurkat_Roadmap_Epigenomics_links_tab)
mpra_emVars_Jurkat_Roadmap_Epigenomics_links_fisher_tb <- as_tibble(t(c(
	"P-value"=mpra_emVars_Jurkat_Roadmap_Epigenomics_links_fisher$p.value, 
	"Fold enrichment"=mpra_emVars_Jurkat_Roadmap_Epigenomics_links_fisher$estimate[["odds ratio"]], 
	"Number of emVars"=sum(mpra_emVars_Jurkat_Roadmap_Epigenomics_links$MPRA_Jurkat_emVars
)))) %>% mutate(`Cell line`="Jurkat")


mpra_emVars_Roadmap_Epigenomics_links_fisher_tb <- bind_rows(
	mpra_emVars_K562_Roadmap_Epigenomics_links_fisher_tb, 
	mpra_emVars_Jurkat_Roadmap_Epigenomics_links_fisher_tb
) %>% 
	mutate(`Cell line`=factor(`Cell line`, levels = c("K562", "Jurkat")))

ggplot(mpra_emVars_Roadmap_Epigenomics_links_fisher_tb) + 
	geom_bar(aes(`Cell line`, `Fold enrichment`), stat="identity", fill="#6baed6") + 
	geom_text(aes(`Cell line`, 0, label=`Number of emVars`), vjust=1, color="#1b9e77") + 
	geom_text(aes(`Cell line`, `Fold enrichment`, label=round(-log10(`P-value`), 2)), vjust=0, color="#7570b3") + 
	geom_hline(yintercept=1, color="#eeac6f") + 
	theme_bw() + 
	ggtitle("Enrichment of emVars with Roadmap \n Epigenomics gene-enhancer links") +
	theme(
		aspect.ratio=1.25, 
		panel.grid.major = element_blank(),
		panel.grid.minor = element_blank(), 
		plot.title = element_text(hjust = 0.5)
	)
ggsave("../../results/3-mpra_sanity_check_analyses/mpra_emVars_Roadmap_Epigenomics_links_enrichment-K562-Jurkat.pdf", scale=0.6)
