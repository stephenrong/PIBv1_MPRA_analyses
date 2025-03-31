#!/bin/R

# command line
library(optparse)

option_list = list(
	make_option(
		c("-i", "--input"), type="character", default=NULL, help="input vcf/bed (un/compressed) [default: %default]", metavar="character"), 
	make_option(
		c("-t", "--target"), type="character", default=NULL, help="target bed (un/compressed) [default: %default]", metavar="character"), 
	make_option(
		c("-b", "--background"), type="character", default=NULL, help="background bed (un/compressed) [default: %default]", metavar="character"), 
	make_option(
		c("-o", "--output"), type="character", default=NULL, help="output file (uncompressed) [default: %default]", metavar="character"),
	make_option(
		c("-z", "--compressed"), action="store_true", default=FALSE, help="logical, compress output file [default: %default]", metavar="character"),
	make_option(
		c("-a", "--autosomes"), action="store_true", default=FALSE, help="logical, use autosomes only [default: %default]", metavar="character")
)

opt <- parse_args(OptionParser(option_list=option_list))
file_input <- opt$input
file_target <- opt$target
file_background <- opt$background
file_output <- opt$output
flag_compressed <- opt$compressed
flag_autosomes <- opt$autosomes

# for data analysis
library(tidyverse)
library(data.table)

# for bio data analysis
library(plyranges)
library(rtracklayer)
library(GenomeInfoDb)
library(vcfR)

# load files
drop_mcols <- function(gr) {
	mcols(gr) <- NULL
	return(gr)
}

seqstyle_change <- function(gr, style) {
	seqlevelsStyle(gr) <- style  # NCBI or UCSC
	return(gr)
}

input_gr <- tryCatch(
	# input is vcf file
	expr = {
		as_tibble(read.vcfR(file_input)@fix) %>% 
			mutate(seqnames=CHROM, start=POS, end=POS, width=1, strand="*") %>% 
			dplyr::select(seqnames, start, end, width, strand) %>% 
			GRanges() %>% 
				drop_mcols() %>% 
				seqstyle_change("NCBI")
	},
	# input is bed file
	error = function(e) {
		import(file_input, format="BED") %>% 
			drop_mcols() %>% 
			seqstyle_change("NCBI")
	}
)

target_gr <- tryCatch(
	# input is bed file
	expr = {
		import(file_target, format="BED") %>% 
			drop_mcols() %>% 
			seqstyle_change("NCBI")
	},
	# input is narrowpeak file
	error = function(e) {
		import(file_target, format="narrowPeak") %>% 
			drop_mcols() %>% 
			seqstyle_change("NCBI")
	}
)

background_gr <- tryCatch(
	# input is bed file
	expr = {
		import(file_background, format="BED") %>% 
			drop_mcols() %>% 
			seqstyle_change("NCBI")
	},
	# input is narrowpeak file
	error = function(e) {
		import(file_background, format="narrowPeak") %>% 
			drop_mcols() %>% 
			seqstyle_change("NCBI")
	}
)

# autosomes
if (flag_autosomes) {
	print("autosomes only")
	input_gr <- input_gr %>% 
		filter(seqnames %in% c(1:22))
	target_gr <- target_gr %>% 
		filter(seqnames %in% c(1:22))
	background_gr <- background_gr %>% 
		filter(seqnames %in% c(1:22))
}

# count overlap
a <- length(findOverlaps(input_gr, target_gr))
b <- length(findOverlaps(input_gr, background_gr))

# count total length
c <- sum(width(target_gr))
d <- sum(width(background_gr))

# calculate enrichment
fisher <- fisher.test(matrix(c(a,b,c-a,d-b), nrow=2))
output_tb <- data.frame(list(
	Fisher_Exact_Test_n1=a, Fisher_Exact_Test_n2=b, 
	Fisher_Exact_Test_n3=c, Fisher_Exact_Test_n4=d, 
	Fisher_Exact_Test_p_value=fisher$p.value[[1]], 
	Fisher_Exact_Test_odds_ratio=fisher$estimate[[1]], 
	Fisher_Exact_Test_LCI=fisher$conf.int[[1]], 
	Fisher_Exact_Test_UCI=fisher$conf.int[[2]])
)

# save output
if (!flag_compressed) {
	write_tsv(output_tb, file_output)
} else {
	write_tsv(output_tb, gzfile(file_output))
}
