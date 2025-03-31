#!/bin/R

# for data analysis
library(tidyverse)
library(data.table)

# for bio data analysis
library(plyranges)
library(vcfR)

# load TSV
adaptive_variants_CADD <- as_tibble(fread("../../results/2f-annotate_variants_CADD_scores/adaptive_variants-whole_genome_SNVs_inclAnno.tsv.gz"))

# add prefix and save names
names(adaptive_variants_CADD)[5:length(adaptive_variants_CADD)] <- paste("CADD_orig", names(adaptive_variants_CADD)[5:length(adaptive_variants_CADD)], sep="_")
adaptive_variants_CADD_names <- names(adaptive_variants_CADD)[5:length(adaptive_variants_CADD)]

# add GRange and Variant cols
adaptive_variants_CADD <- adaptive_variants_CADD %>% 
	dplyr::rename(Chrom = `#Chrom`) %>% 
	mutate(VariantCHROM = Chrom, VariantPOS = Pos, VariantREF = Ref, VariantALT = Alt) %>% 
	mutate(VariantID = paste(VariantCHROM, VariantPOS, VariantREF, VariantALT, sep="_")) %>% 
	mutate(seqnames=VariantCHROM, start=VariantPOS, end=VariantPOS, width=1, strand="*") %>% 
	dplyr::select(seqnames, start, end, width, strand, VariantID, VariantCHROM, VariantPOS, VariantREF, VariantALT, CADD_orig_PHRED)  # all_of(adaptive_variants_CADD_names))
adaptive_variants_CADD_names <- names(adaptive_variants_CADD)[11:length(adaptive_variants_CADD)]

# prefilter variants
adaptive_variants_tb <- as_tibble(fread("../../results/1a-preprocess_PIBv1_MPRA_pilot/adaptive_variants.txt.gz")) %>% 
	dplyr::select(seqnames, start, end, width, strand, VariantID, VariantCHROM, VariantPOS, VariantREF, VariantALT)
adaptive_variants_CADD <- adaptive_variants_CADD %>% 
	filter(VariantID %in% adaptive_variants_tb$VariantID)

# collapse multiple entries
adaptive_variants_CADD <- adaptive_variants_CADD %>% 
	group_by(VariantID) %>% 
	arrange(desc(CADD_orig_PHRED)) %>% 
	dplyr::summarise_at(adaptive_variants_CADD_names, function(x) {paste(x, collapse=",")}) %>% 
	ungroup()

# add summary cols
adaptive_variants_CADD <- adaptive_variants_CADD %>% 
	mutate(CADD_summ_Score = as.numeric(gsub(",.*", "", CADD_orig_PHRED))) %>% 
	mutate(CADD_summ_Summary = ifelse(CADD_summ_Score>=30, "CADD>=30", ifelse(CADD_summ_Score>=20, "CADD>=20", ifelse(CADD_summ_Score>=10, "CADD>=10", NA))))

# merge to full variants
names_temp <- intersect(names(adaptive_variants_tb), names(adaptive_variants_CADD))
adaptive_variants_CADD <- adaptive_variants_tb %>% mutate_at(names_temp, as.character)%>% 
	left_join(adaptive_variants_CADD %>% mutate_at(names_temp, as.character))

# save all variants
write_tsv(adaptive_variants_CADD, gzfile("../../results/2f-annotate_variants_CADD_scores/adaptive_variants_CADD.txt.gz"))

# save summ variants
adaptive_variants_CADD_summ <- adaptive_variants_CADD %>% 
	dplyr::select(seqnames, start, end, width, strand, starts_with("Variant"), contains("_summ_"))
write_tsv(adaptive_variants_CADD_summ, gzfile("../../results/2f-annotate_variants_CADD_scores/adaptive_variants_CADD_summ.txt.gz"))
sort(table(adaptive_variants_CADD_summ$CADD_summ_Summary))

# visual check
pdf("../../results/2f-annotate_variants_CADD_scores/adaptive_variants_CADD_summ_sort.pdf")
plot(sort(adaptive_variants_CADD$CADD_summ_Score), ylab="CADD PHRED score")
dev.off()

pdf("../../results/2f-annotate_variants_CADD_scores/adaptive_variants_CADD_summ_pie.pdf")
pie(rev(sort(table(adaptive_variants_CADD_summ$CADD_summ_Summary))))
dev.off()
