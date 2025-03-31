#!/bin/R

# for data analysis
library(tidyverse)
library(data.table)

# for bio data analysis
library(plyranges)
library(vcfR)

# load TSV
introgressed_variants_CADD <- as_tibble(fread("../../results/2f-annotate_variants_CADD_scores/introgressed_variants-whole_genome_SNVs_inclAnno.tsv.gz"))

# add prefix and save names
names(introgressed_variants_CADD)[5:length(introgressed_variants_CADD)] <- paste("CADD_orig", names(introgressed_variants_CADD)[5:length(introgressed_variants_CADD)], sep="_")
introgressed_variants_CADD_names <- names(introgressed_variants_CADD)[5:length(introgressed_variants_CADD)]

# add GRange and Variant cols
introgressed_variants_CADD <- introgressed_variants_CADD %>% 
	dplyr::rename(Chrom = `#Chrom`) %>% 
	mutate(VariantCHROM = Chrom, VariantPOS = Pos, VariantREF = Ref, VariantALT = Alt) %>% 
	mutate(VariantID = paste(VariantCHROM, VariantPOS, VariantREF, VariantALT, sep="_")) %>% 
	mutate(seqnames=VariantCHROM, start=VariantPOS, end=VariantPOS, width=1, strand="*") %>% 
	dplyr::select(seqnames, start, end, width, strand, VariantID, VariantCHROM, VariantPOS, VariantREF, VariantALT, CADD_orig_PHRED)  # all_of(introgressed_variants_CADD_names))
introgressed_variants_CADD_names <- names(introgressed_variants_CADD)[11:length(introgressed_variants_CADD)]

# prefilter variants
introgressed_variants_tb <- as_tibble(fread("../../results/1a-preprocess_PIBv1_MPRA_final/introgressed_variants.txt.gz")) %>% 
	dplyr::select(seqnames, start, end, width, strand, VariantID, VariantCHROM, VariantPOS, VariantREF, VariantALT)
introgressed_variants_CADD <- introgressed_variants_CADD %>% 
	filter(VariantID %in% introgressed_variants_tb$VariantID)

# collapse multiple entries
introgressed_variants_CADD <- introgressed_variants_CADD %>% 
	group_by(seqnames, start, end, width, strand, VariantID, VariantCHROM, VariantPOS, VariantREF, VariantALT) %>% 
	arrange(desc(CADD_orig_PHRED)) %>% 
	dplyr::summarise_at(introgressed_variants_CADD_names, function(x) {paste(x, collapse=",")}) %>% 
	ungroup()

# add summary cols
introgressed_variants_CADD <- introgressed_variants_CADD %>% 
	mutate(CADD_summ_Score = as.numeric(gsub(",.*", "", CADD_orig_PHRED))) %>% 
	mutate(CADD_summ_Summary = ifelse(CADD_summ_Score>=40, "CADD>=40", ifelse(CADD_summ_Score>=30, "CADD>=30", ifelse(CADD_summ_Score>=20, "CADD>=20", ifelse(CADD_summ_Score>=10, "CADD>=10", NA)))))

# merge to full variants
names_temp <- intersect(names(introgressed_variants_tb), names(introgressed_variants_CADD))
introgressed_variants_CADD <- introgressed_variants_tb %>% mutate_at(names_temp, as.character)%>% 
	left_join(introgressed_variants_CADD %>% mutate_at(names_temp, as.character))

# save all variants
write_tsv(introgressed_variants_CADD, gzfile("../../results/2f-annotate_variants_CADD_scores/introgressed_variants_CADD.txt.gz"))

# save summ variants
introgressed_variants_CADD_summ <- introgressed_variants_CADD %>% 
	dplyr::select(seqnames, start, end, width, strand, starts_with("Variant"), contains("_summ_"))
write_tsv(introgressed_variants_CADD_summ, gzfile("../../results/2f-annotate_variants_CADD_scores/introgressed_variants_CADD_summ.txt.gz"))
sort(table(introgressed_variants_CADD_summ$CADD_summ_Summary))

# visual check
pdf("../../results/2f-annotate_variants_CADD_scores/introgressed_variants_CADD_summ_sort.pdf")
plot(sort(introgressed_variants_CADD$CADD_summ_Score), ylab="CADD PHRED score")
dev.off()

pdf("../../results/2f-annotate_variants_CADD_scores/introgressed_variants_CADD_summ_pie.pdf")
pie(rev(sort(table(introgressed_variants_CADD_summ$CADD_summ_Summary))))
dev.off()
