#!/bin/R

# Revision uses background gene set for Enricher analysis

library(tidyverse)
library(data.table)
library(ggpubr)

main_title <- "Reactome (2024)"
plot_scale = 0.8
n_show_thresh= 15
fdr_show_thresh = 1.0
fdr_threshold = 0.01
ylim_max = 7
min_genes= 2 

coding_splicing_emVars <- as_tibble(fread("Gene_set_enrichment_coding_splicing_emVars_revision/Enrichr_coding_splicing_emVars_revision/Reactome_Pathways_2024_table.txt"))

coding_splicing_emVars <- coding_splicing_emVars %>% 
	filter(`Adjusted P-value` < fdr_show_thresh) %>%  # remove gene sets above fdr show threshold
	rowwise() %>% filter(length(strsplit(Genes, ";")[[1]]) >= min_genes) %>% ungroup() %>%  # remove gene sets with single gene
	mutate(Genes_dup = !duplicated(Genes)) %>% filter(Genes_dup) %>% dplyr::select(-Genes_dup) %>%  # remove if same set of genes but lower p-value
	dplyr::select(Term, Overlap, `P-value`, `Adjusted P-value`, `Odds Ratio`, `Combined Score`, Genes) %>% 
	mutate(Genes = gsub(";", ", ", Genes)) %>% 
	mutate(Term = gsub(" R-HSA-.*", "", Term)) %>% 
	mutate(`P-value` = as.numeric(format(`P-value`, scientific = TRUE, digits = 3))) %>% 
	mutate(`Adjusted P-value` = as.numeric(format(`Adjusted P-value`, scientific = TRUE, digits = 3)))

coding_splicing_emVars <- coding_splicing_emVars %>% 
	arrange(`Adjusted P-value`) %>% 
	mutate(index = row_number()) %>% 
	filter(index %in% c(1:n_show_thresh))

coding_splicing_only <- as_tibble(fread("Gene_set_enrichment_coding_splicing_only_revision/Enrichr_coding_splicing_only_revision/Reactome_Pathways_2024_table.txt"))

coding_splicing_only <- coding_splicing_only %>% 
	filter(`Adjusted P-value` < fdr_show_thresh) %>%  # remove gene sets above fdr show threshold
	rowwise() %>% filter(length(strsplit(Genes, ";")[[1]]) >= min_genes) %>% ungroup() %>%  # remove gene sets with single gene
	mutate(Genes_dup = !duplicated(Genes)) %>% filter(Genes_dup) %>% dplyr::select(-Genes_dup) %>%  # remove if same set of genes but lower p-value
	dplyr::select(Term, Overlap, `P-value`, `Adjusted P-value`, `Odds Ratio`, `Combined Score`, Genes) %>% 
	mutate(Genes = gsub(";", ", ", Genes)) %>% 
	mutate(Term = gsub(" R-HSA-.*", "", Term)) %>% 
	mutate(`P-value` = as.numeric(format(`P-value`, scientific = TRUE, digits = 3))) %>% 
	mutate(`Adjusted P-value` = as.numeric(format(`Adjusted P-value`, scientific = TRUE, digits = 3))) %>% 
	filter(Term %in% coding_splicing_emVars$Term)

combined_figure <- bind_rows(mutate(coding_splicing_only, `Set`="coding + splicing"), mutate(coding_splicing_emVars, `Set`="coding + splicing + emVars"))

index_map <- coding_splicing_emVars[c("Term", "index")] %>% deframe()

combined_figure <- combined_figure %>% 
	mutate(index = index_map[Term])
write_tsv(combined_figure, "Supplementary_tables_and_figures_revision/Supplementary_X-Gene_set_enrichment-combined_figures-Reactome_Pathways_2024_table.txt")

ggplot(combined_figure) +
	geom_linerange(aes(x = reorder(gsub(" R-HSA-.*", "", Term), -index), ymin = 0, ymax = -log10(`P-value`), color = `Set`), 
		position = position_dodge(width = 1)) +
	geom_point(aes(x = reorder(gsub(" R-HSA-.*", "", Term), -index), y = -log10(`P-value`), color = `Set`, alpha=(`Adjusted P-value` < fdr_threshold)),
		position = position_dodge(width = 1), size=4) + 
	coord_flip() +
	xlab("") +
	ylab(expression("-log"["10"] ~ "P-value")) +
	theme_pubr() + 
	theme(aspect.ratio=1) + 
	theme(axis.text.y = element_text(size=10)) + 
	theme(legend.position="bottom") + 
	labs(title = main_title) + 
	ylim(c(0, ylim_max)) + 
	scale_alpha_manual(name = paste0("FDR < ", fdr_threshold), breaks = c(TRUE, FALSE), values = c(1, 0.5))
ggsave("Supplementary_tables_and_figures_revision/Supplementary_X-Gene_set_enrichment-combined_figure-Reactome_Pathways_2024_table.pdf", scale=plot_scale)

# fdr version
combined_figure_replot <- combined_figure %>% 
	mutate(`P-value` = -log10(`P-value`)) %>% 
	mutate(`Adjusted P-value` = -log10(`Adjusted P-value`)) %>% 
	dplyr::select(Term, `Adjusted P-value`, `Set`, index) %>% 
	pivot_wider(names_from = `Set`, values_from = `Adjusted P-value`)

ggplot(combined_figure_replot) + 
	geom_linerange(aes(x = reorder(gsub(" R-HSA-.*", "", Term), -index), ymin = `coding + splicing`, ymax = `coding + splicing + emVars`)) + 
	geom_point(aes(x = reorder(gsub(" R-HSA-.*", "", Term), -index), y = -log10(`Adjusted P-value`), color = `Set`, shape = `Set`), data=combined_figure, size=3) + 
	geom_hline(yintercept=-log10(fdr_threshold), color="grey", linetype="dashed") + 
	coord_flip() + 
	xlab("") +
	ylab(expression("-log"["10"] ~ "FDR")) +
	theme_pubr() + 
	theme(aspect.ratio=1.5) + 
	theme(axis.text.y = element_text(size=10)) + 
	theme(legend.position="bottom") + 
	labs(title = main_title) + 
	ylim(c(0, 4)) + 
	theme(legend.title=element_blank())
ggsave("Supplementary_tables_and_figures_revision/Supplementary_X-Gene_set_enrichment-combined_figure-Reactome_Pathways_2024_table-replot.pdf", scale=plot_scale)
