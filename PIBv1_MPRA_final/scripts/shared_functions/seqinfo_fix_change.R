#!/bin/R

library(plyranges)
# library(GenomeInfoDb)

# # save cache
# ChromInfoFromNCBI_GRCh37 <- getChromInfoFromNCBI("GRCh37", as.Seqinfo=TRUE)[c(1:22,"X","Y")]
# ChromInfoFromNCBI_GRCh38 <- getChromInfoFromNCBI("GRCh38", as.Seqinfo=TRUE)[c(1:22,"X","Y")]
# ChromInfoFromUCSC_hg19 <- getChromInfoFromUCSC("hg19", as.Seqinfo=TRUE)[paste("chr", c(1:22,"X","Y"), sep="")]
# ChromInfoFromUCSC_hg38 <- getChromInfoFromUCSC("hg38", as.Seqinfo=TRUE)[paste("chr", c(1:22,"X","Y"), sep="")]
# # https://github.com/Bioconductor/GenomeInfoDb/issues/82

# saveRDS(ChromInfoFromNCBI_GRCh37, "../../data/shared_functions/ChromInfoFromNCBI_GRCh37.rds")
# saveRDS(ChromInfoFromNCBI_GRCh38, "../../data/shared_functions/ChromInfoFromNCBI_GRCh38.rds")
# saveRDS(ChromInfoFromUCSC_hg19, "../../data/shared_functions/ChromInfoFromUCSC_hg19.rds")
# saveRDS(ChromInfoFromUCSC_hg38, "../../data/shared_functions/ChromInfoFromUCSC_hg38.rds")

# load cache
ChromInfoFromNCBI_GRCh37 <- readRDS("../../data/shared_functions/ChromInfoFromNCBI_GRCh37.rds")
ChromInfoFromNCBI_GRCh38 <- readRDS("../../data/shared_functions/ChromInfoFromNCBI_GRCh38.rds")
ChromInfoFromUCSC_hg19 <- readRDS("../../data/shared_functions/ChromInfoFromUCSC_hg19.rds")
ChromInfoFromUCSC_hg38 <- readRDS("../../data/shared_functions/ChromInfoFromUCSC_hg38.rds")

seqinfo_fix <- function(gr, style, genome) {
	if (style == "NCBI") {
		if (genome == "GRCh37") {
			seqorig <- ChromInfoFromNCBI_GRCh37
		} else if (genome == "GRCh38") {
			seqorig <- ChromInfoFromNCBI_GRCh38
		}
		gr <- gr %>% seqstyle_change("NCBI")
		seqlevels(gr) <- seqlevels(seqorig)[which(seqlevels(seqorig) %in% seqlevels(gr))]
		seqinfo(gr) <- seqorig[seqlevels(gr)[order(match(seqlevels(gr),seqlevels(seqorig)))]]
	} else if (style == "UCSC") {
		if (genome == "hg19") {
			seqorig <- ChromInfoFromUCSC_hg19
		} else if (genome == "hg38") {
			seqorig <- ChromInfoFromUCSC_hg38
		}
		gr <- gr %>% seqstyle_change("UCSC")
		seqlevels(gr) <- seqlevels(seqorig)[which(seqlevels(seqorig) %in% seqlevels(gr))]
		seqinfo(gr) <- seqorig[seqlevels(gr)[order(match(seqlevels(gr),seqlevels(seqorig)))]]
	}
	return(gr)
}

seqstyle_change <- function(gr, style) {
	seqlevelsStyle(gr) <- style
	return(gr)
}
