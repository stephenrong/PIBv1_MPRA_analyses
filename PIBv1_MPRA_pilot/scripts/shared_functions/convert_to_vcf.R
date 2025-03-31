#!/bin/R

# for data analysis
library(tidyverse)
library(data.table)
library(wrapr)

# for bio data analysis
library(plyranges)
source("../shared_functions/seqinfo_fix_change.R")

# custom file conversions
convert_gr_to_vcf <- function(gr, out_file, clean=FALSE, complete=TRUE, compress=TRUE) {
	# gr must have VariantID, VariantCHROM, VariantPOS, VariantREF, VariantALT columns
	tb <- gr %>% as_tibble() %>% 
		dplyr::select(-seqnames, -start, -end, -width, -strand)
	genome <- seqinfo(gr)@genome[1]

	# vcf_fixed
	vcf_fixed <- gr %>% as_tibble() %>% 
		dplyr::select(seqnames, start, VariantID, VariantREF, VariantALT)
	names(vcf_fixed) <- c("#CHROM", "POS", "ID", "REF", "ALT")
	vcf_fixed <- vcf_fixed %>% mutate(QUAL = ".", FILTER = ".", INFO = ".")

	# vcf_info
	if (complete) {
		vcf_info <- tb %>% 
			dplyr::select(setdiff(names(tb), c(names(vcf_fixed))))	
	}
	vcf_info_temp <- vcf_info
	for (i in names(vcf_info_temp)) {
		let(c(VAR=i), vcf_info_temp <- vcf_info_temp %>% 
			mutate(VAR = paste(i, VAR, sep="=")))
	}
	if (!clean) {
		vcf_fixed$INFO <- vcf_info_temp %>% 
			unite(col=INFO, sep=";") %>% .$INFO
		vcf_fixed <- vcf_fixed %>% 
			arrange(`#CHROM`, POS, REF, ALT) %>% 
			mutate(across(everything(), as.character))
	}

	# vcf_contig
	vcf_contig_temp <- as_tibble(rownames_to_column(as.data.frame(
		seqinfo(gr)), "seqname"))
	vcf_contig_temp <- vcf_contig_temp %>% 
		mutate(seqname = gsub("chr", "", seqname)) %>% 
		filter(seqname %in% vcf_fixed$`#CHROM`)
	vcf_contig_temp <- vcf_contig_temp %>% 
		mutate(contig_meta = paste(
			"##contig=<", 
			"ID=", seqname, 
			",length=", seqlengths, 
			",assembly=", "\"\"", 
			",species=", "\"\"", 
			">", sep=""))
	vcf_contig <- vcf_contig_temp$contig_meta

	# vcf_meta
	# 	map r types to vcf format types
	map <- c("logical"="Integer", "integer"="Integer", "double"="Float","character"="String", 
		"complex"="String", "raw"="String", "list"="String", "factor"="String", 
		"ordered"="String", "Date"="String", "POSIXt"="String", "difftime"="String")
	vcf_meta <- c(
		"##fileformat=VCFv4.2",
		paste("##fileDate=", gsub("-", "", Sys.Date()), sep=""), 
		paste("##source=", "convert_gr_to_vcf", sep=""),
		paste("##reference=", genome, sep=""),
		vcf_contig
	)
	if (!clean) {
		vcf_meta <- c(
			vcf_meta,
			as.vector(sapply(names(vcf_info_temp), function(x) {
				paste("##INFO=<ID=", x, ",Number=1,Type=", 
					map[class(vcf_info_temp[[x]])], 
					",Description=\"", x, "\">", sep="")
			}))
		)
	}

	# # concatenate vcf file
	# vcf <- c(vcf_meta, paste(names(vcf_fixed), collapse="\t"), 
	# 	matrix(apply(vcf_fixed, 1, paste, collapse='\t'), ncol=1))
	# return(vcf)

	# output vcf file
	writeLines(vcf_meta, out_file)
	write_tsv(vcf_fixed, out_file, append=TRUE, col_names=TRUE)
	if (compress) {
		system(paste("bgzip -f ", out_file, sep=""))
		system(paste("tabix -p vcf ", out_file, ".gz", sep=""))
	}
}

convert_tb_to_vcf <- function(tb, style, genome, out_file, complete=TRUE, compress=TRUE) {
	gr <- tb %>% seqinfo_fix(style, genome)
	convert_gr_to_vcf(gr, out_file, complete=complete, compress=TRcompressUE)
}
