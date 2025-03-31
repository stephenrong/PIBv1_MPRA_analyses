#!/bin/R

# load data packages
library(tidyverse)
library(data.table)

# load bio data packages
library(plyranges)
library(rtracklayer)
source("../shared_functions/seqinfo_fix_change.R")

# load gwas catalog
gwas_catalog_associations <- as_tibble(fread("../../../../Datasets/associated_regions_gwas/GWAS_Catalog/data_download/e109_r2023-07-05/gwas_catalog_v1.0.2-associations_e109_r2023-07-05.tsv.gz", quote=""))
gwas_catalog_trait_mappings <- as_tibble(fread("../../../../Datasets/associated_regions_gwas/GWAS_Catalog/data_download/e109_r2023-07-05/gwas_catalog_trait-mappings_r2023-07-05.tsv.gz", quote=""))

gwas_catalog_associations_mapped <- gwas_catalog_associations %>% 
	filter(CHR_POS != "") %>%  # must have a chromosomal position
	filter(!grepl("x|;", CHR_POS)) %>%  # must have single position
	filter(as.numeric(`P-VALUE`) <= 5e-6) %>% 
	mutate(`Risk allele` = gsub(".*-", "", `STRONGEST SNP-RISK ALLELE`)) %>% 
	filter(`Risk allele` %in% c("A", "C", "T", "G"))  # must be SNP

# convert to grange
gwas_catalog_associations_mapped <- gwas_catalog_associations_mapped %>% 
	mutate(seqnames = CHR_ID, start = as.numeric(CHR_POS)) %>% 
	rowwise() %>% mutate(end = start) %>% ungroup() %>% 
	filter(seqnames %in% 1:22)

gwas_catalog_associations_mapped <- GRanges(gwas_catalog_associations_mapped) %>% 
	seqinfo_fix("NCBI", "GRCh38") %>% sort()

# liftOver to hg19
hg38toHg19 <- import.chain("../../../../Datasets/reference_genomes/liftOver_chains/hg38ToHg19.over.chain")

gwas_catalog_associations_mapped_lift37 <- gwas_catalog_associations_mapped %>% 
	seqinfo_fix("UCSC", "hg38") %>% 
	liftOver(hg38toHg19) %>% unlist() %>% 
	seqinfo_fix("NCBI", "GRCh37") %>% sort()

# save txt tables
write_tsv(as_tibble(gwas_catalog_associations_mapped), gzfile("../../results/2t-annotate_variants_GWAS_Catalog/gwas_catalog_associations_mapped.txt.gz"))
write_tsv(as_tibble(gwas_catalog_associations_mapped_lift37), gzfile("../../results/2t-annotate_variants_GWAS_Catalog/gwas_catalog_associations_mapped_lift37.txt.gz"))
saveRDS(gwas_catalog_associations_mapped, "../../results/2t-annotate_variants_GWAS_Catalog/gwas_catalog_associations_mapped.rds")
saveRDS(gwas_catalog_associations_mapped_lift37, "../../results/2t-annotate_variants_GWAS_Catalog/gwas_catalog_associations_mapped_lift37.rds")
