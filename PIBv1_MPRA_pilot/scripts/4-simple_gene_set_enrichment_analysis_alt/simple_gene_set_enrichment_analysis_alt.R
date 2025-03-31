#!/bin/R

# load packages
library(tidyverse)
library(data.table)

# load bio packages
library(plyranges)

# load enrichments
Enrichr_GO_BP_adaptive_all_genomewide_simple_enrichment <- as_tibble(fread("../../results/4-simple_gene_set_enrichment_analysis_alt/emVar_linked_adaptive_Enrichr/GO_Biological_Process_2021_table-emVar_linked.txt"))

Enrichr_GO_BP_adaptive_class_genomewide_simple_enrichment <- bind_rows(
	Enrichr_GO_BP_adaptive_all_genomewide_simple_enrichment %>% mutate(Classification = "")
) %>% mutate(Classification = factor(Classification, levels = c("")))

# get significant
Enrichr_GO_BP_adaptive_class_genomewide_simple_enrichment <- Enrichr_GO_BP_adaptive_class_genomewide_simple_enrichment %>% 
	mutate(`Signif` = (`Adjusted P-value` < 0.1))

# filter significant
Enrichr_GO_BP_adaptive_class_genomewide_simple_enrichment_signif <- Enrichr_GO_BP_adaptive_class_genomewide_simple_enrichment %>% 
	filter(`Adjusted P-value` < 0.1)  # FDR 10%
Enrichr_GO_BP_adaptive_class_genomewide_simple_enrichment <- Enrichr_GO_BP_adaptive_class_genomewide_simple_enrichment %>% 
	filter(`Term` %in% Enrichr_GO_BP_adaptive_class_genomewide_simple_enrichment_signif$`Term`)

# plot dotplot
ggplot(Enrichr_GO_BP_adaptive_class_genomewide_simple_enrichment, 
	aes(y=reorder(gsub("\\(.*", "", `Term`), log2(`Odds Ratio`)), 
		x=`Classification`,
		fill=log2(`Odds Ratio`), 
		color=`Signif`,
		size=as.integer(gsub("/.*", "", `Overlap`))
	)) + 
	geom_point(
		shape=21,
		stroke=1
	) + 
	scale_fill_gradient2(low="#D7191D", mid="white", high="#2C7BB6") + 
	scale_color_manual(breaks=c(TRUE, FALSE), values=c("black", "white")) + 
	labs(y=NULL, x=NULL, color="FDR < 0.1", fill=bquote('Log'[2]~'fold-change'), size="Genes in set") + 
	theme_bw() + 
	ggtitle("GO Biological Process") + 
	theme(
		axis.text.x=element_text(angle=45, hjust=1), 
		axis.ticks = element_blank(),
		panel.grid.major=element_blank(), 
		panel.grid.minor=element_blank(), 
		aspect.ratio=5, 
		plot.title = element_text(hjust = 0.5)  # , 
	)
ggsave("../../results/4-simple_gene_set_enrichment_analysis_alt/Enrichr_GO_Biological_Process_2021_adaptive_class_genomewide_simple_enrichment_plot.pdf", scale=0.9)


# load enrichments
Enrichr_GO_BP_adaptive_all_genomewide_simple_enrichment <- as_tibble(fread("../../results/4-simple_gene_set_enrichment_analysis_alt/emVar_linked_adaptive_Enrichr/GO_Biological_Process_2023_table-emVar_linked.txt"))

Enrichr_GO_BP_adaptive_class_genomewide_simple_enrichment <- bind_rows(
	Enrichr_GO_BP_adaptive_all_genomewide_simple_enrichment %>% mutate(Classification = "")
) %>% mutate(Classification = factor(Classification, levels = c("")))

# get significant
Enrichr_GO_BP_adaptive_class_genomewide_simple_enrichment <- Enrichr_GO_BP_adaptive_class_genomewide_simple_enrichment %>% 
	mutate(`Signif` = (`Adjusted P-value` < 0.1))

# filter significant
Enrichr_GO_BP_adaptive_class_genomewide_simple_enrichment_signif <- Enrichr_GO_BP_adaptive_class_genomewide_simple_enrichment %>% 
	filter(`Adjusted P-value` < 0.1)  # FDR 10%
Enrichr_GO_BP_adaptive_class_genomewide_simple_enrichment <- Enrichr_GO_BP_adaptive_class_genomewide_simple_enrichment %>% 
	filter(`Term` %in% Enrichr_GO_BP_adaptive_class_genomewide_simple_enrichment_signif$`Term`)

# plot dotplot
ggplot(Enrichr_GO_BP_adaptive_class_genomewide_simple_enrichment, 
	aes(y=reorder(gsub("\\(.*", "", `Term`), log2(`Odds Ratio`)), 
		x=`Classification`,
		fill=log2(`Odds Ratio`), 
		color=`Signif`,
		size=as.integer(gsub("/.*", "", `Overlap`))
	)) + 
	geom_point(
		shape=21,
		stroke=1
	) + 
	scale_fill_gradient2(low="#D7191D", mid="white", high="#2C7BB6") + 
	scale_color_manual(breaks=c(TRUE, FALSE), values=c("black", "white")) + 
	labs(y=NULL, x=NULL, color="FDR < 0.1", fill=bquote('Log'[2]~'fold-change'), size="Genes in set") + 
	theme_bw() + 
	ggtitle("GO Biological Process") + 
	theme(
		axis.text.x=element_text(angle=45, hjust=1), 
		axis.ticks = element_blank(),
		panel.grid.major=element_blank(), 
		panel.grid.minor=element_blank(), 
		aspect.ratio=5, 
		plot.title = element_text(hjust = 0.5)  # , 
	)
ggsave("../../results/4-simple_gene_set_enrichment_analysis_alt/Enrichr_GO_Biological_Process_2023_adaptive_class_genomewide_simple_enrichment_plot.pdf", scale=0.7)


# load enrichments
Enrichr_Reactome_adaptive_all_genomewide_simple_enrichment <- as_tibble(fread("../../results/4-simple_gene_set_enrichment_analysis_alt/emVar_linked_adaptive_Enrichr/Reactome_2022_table-emVar_linked.txt"))

Enrichr_Reactome_adaptive_class_genomewide_simple_enrichment <- bind_rows(
	Enrichr_Reactome_adaptive_all_genomewide_simple_enrichment %>% mutate(Classification = "")
) %>% mutate(Classification = factor(Classification, levels = c("")))

# get significant
Enrichr_Reactome_adaptive_class_genomewide_simple_enrichment <- Enrichr_Reactome_adaptive_class_genomewide_simple_enrichment %>% 
	mutate(`Signif` = (`Adjusted P-value` < 0.1))

# filter significant
Enrichr_Reactome_adaptive_class_genomewide_simple_enrichment_signif <- Enrichr_Reactome_adaptive_class_genomewide_simple_enrichment %>% 
	filter(`Adjusted P-value` < 0.1)  # FDR 10%
Enrichr_Reactome_adaptive_class_genomewide_simple_enrichment <- Enrichr_Reactome_adaptive_class_genomewide_simple_enrichment %>% 
	filter(`Term` %in% Enrichr_Reactome_adaptive_class_genomewide_simple_enrichment_signif$`Term`)

# plot dotplot
ggplot(Enrichr_Reactome_adaptive_class_genomewide_simple_enrichment, 
	aes(y=reorder(gsub("\\R-HSA.*", "", `Term`), log2(`Odds Ratio`)), 
		x=`Classification`,
		fill=log2(`Odds Ratio`), 
		color=`Signif`,
		size=as.integer(gsub("/.*", "", `Overlap`))
	)) + 
	geom_point(
		shape=21,
		stroke=1
	) + 
	scale_fill_gradient2(low="#D7191D", mid="white", high="#2C7BB6") + 
	scale_color_manual(breaks=c(TRUE, FALSE), values=c("black", "white")) + 
	labs(y=NULL, x=NULL, color="FDR < 0.1", fill=bquote('Log'[2]~'fold-change'), size="Genes in set") + 
	theme_bw() + 
	ggtitle("Reactome") + 
	theme(
		axis.text.x=element_text(angle=45, hjust=1), 
		axis.ticks = element_blank(),
		panel.grid.major=element_blank(), 
		panel.grid.minor=element_blank(), 
		aspect.ratio=5, 
		plot.title = element_text(hjust = 0.5)  # , 
	)
ggsave("../../results/4-simple_gene_set_enrichment_analysis_alt/Enrichr_Reactome_2022_adaptive_class_genomewide_simple_enrichment_plot.pdf", scale=0.7)


# helper functions
load_data <- function(filename) {
	# Read all lines
	lines <- readLines(filename)

	# Remove trailing tabs
	lines <- gsub("\t*$", "", lines)
	
	# Find the lines which start with CLUSTER_NAME (our tables)
	table_starts <- grep("^CLUSTER_NAME|^SOURCE", lines)
	
	# Add one more index to represent end of file
	table_starts <- c(table_starts, length(lines) + 1)
	
	# Create empty list to store tables
	tables <- list()
	
	# Loop over table starts and read tables
	for (i in seq_along(table_starts)[-length(table_starts)]) {
		# Get start and end indices for the current table
		start_line <- table_starts[i]
		end_line <- table_starts[i+1] - 1
		
		# Subset lines to those for the current table
		table_lines <- lines[start_line:end_line]
		
		# Create a text connection to read the table from
		table_conn <- textConnection(table_lines)
		
		# Read the table and store in list
		tables[[i]] <- read.table(table_conn, sep="\t", header=TRUE)
		
		# Close the connection
		close(table_conn)
	}
	
	# Return the list of tables
	return(tables)
}

# load enrichments
HumanBase_global_adaptive_all_genomewide_simple_enrichment <- as_tibble(load_data("../../results/4-simple_gene_set_enrichment_analysis_alt/emVar_linked_adaptive_HumanBase/global_1685565106836-emVar_linked/global_1685565106836.tsv")[[2]])

HumanBase_global_adaptive_class_genomewide_simple_enrichment <- bind_rows(
	HumanBase_global_adaptive_all_genomewide_simple_enrichment %>% mutate(Classification = "")
) %>% mutate(Classification = factor(Classification, levels = c("")))

# get significant
HumanBase_global_adaptive_class_genomewide_simple_enrichment <- HumanBase_global_adaptive_class_genomewide_simple_enrichment %>% 
	mutate(`Signif` = (`TERM_Q_VALUE` < 0.05))

# filter significant
HumanBase_global_adaptive_class_genomewide_simple_enrichment_signif <- HumanBase_global_adaptive_class_genomewide_simple_enrichment %>% 
	filter(`TERM_Q_VALUE` < 0.05)	# FDR 10%
HumanBase_global_adaptive_class_genomewide_simple_enrichment <- HumanBase_global_adaptive_class_genomewide_simple_enrichment %>% 
	filter(`TERM_NAME` %in% HumanBase_global_adaptive_class_genomewide_simple_enrichment_signif$`TERM_NAME`) %>% 
	filter(`GENE_COUNT` >= 10)

HumanBase_global_adaptive_class_genomewide_simple_enrichment_order <- HumanBase_global_adaptive_class_genomewide_simple_enrichment %>% 
	group_by(`TERM_NAME`) %>% arrange(-`GENE_COUNT`) %>% slice(1) %>% summarise(Order = as.integer(`Classification`)) %>% ungroup()

HumanBase_global_adaptive_class_genomewide_simple_enrichment <- HumanBase_global_adaptive_class_genomewide_simple_enrichment %>% 
	left_join(HumanBase_global_adaptive_class_genomewide_simple_enrichment_order)

# plot dotplot
ggplot(HumanBase_global_adaptive_class_genomewide_simple_enrichment, 
	aes(y=reorder(`TERM_NAME`, -Order), 
		x=`Classification`,
		fill=-log10(`TERM_Q_VALUE`), 
		color=`Signif`,
		size=`GENE_COUNT`
	)) + 
	geom_point(
		shape=21,
		stroke=1
	) + 
	scale_x_discrete(drop=FALSE) + 
	scale_fill_discrete(drop=FALSE) +
	scale_fill_gradient2(low="#D7191D", mid="white", high="#2C7BB6") + 
	scale_color_manual(breaks=c(TRUE, FALSE), values=c("black", "white")) + 
	labs(y=NULL, x=NULL, color="FDR < 0.05", fill=bquote('-Log'[10]~'FDR'), size="Genes in set") + 
	theme_bw() + 
	ggtitle("HumanBase") + 
	theme(
		axis.text.x=element_text(angle=45, hjust=1), 
		axis.ticks = element_blank(),
		panel.grid.major=element_blank(), 
		panel.grid.minor=element_blank(), 
		aspect.ratio=5, 
		plot.title = element_text(hjust = 0.5)
ggsave("../../results/4-simple_gene_set_enrichment_analysis_alt/HumanBase_global_adaptive_class_genomewide_simple_enrichment_plot.pdf", scale=1)
