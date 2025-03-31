#!/bin/R

# for data analysis
library(tidyverse)
library(data.table)

# custom bed function
create_ucsc_custom_bed_track <- function(input_tb, output_bed, browser_position=NA, name=NA, description=NA, visibility=NA, useScore=NA, itemRgb=NA, colorByStrand=NA) {
	# create header line
	header_opts <- c(name=name, description=description, visibility=visibility, useScore=useScore, itemRgb=itemRgb, colorByStrand=colorByStrand)
	header_opts <- header_opts[which(!is.na(header_opts))]
	header = paste0("track type=bed ", paste(sapply(1:length(header_opts), function(i) {paste0(names(header_opts)[i], "=\"", header_opts[i], "\"")}), collapse=" "))

	# create browser position line
	if (!is.na(browser_position)) {
		browser_position = paste0("browser position ", browser_position)
		write(browser_position, output_bed)
	}

	# save bed
	if (!is.na(browser_position)) {
		write(browser_position, output_bed)
		write(header, output_bed, append=TRUE)
	} else {
		write(header, output_bed)
	}
	write.table(input_tb, output_bed, append=TRUE, quote=FALSE, sep="\t", col.names=FALSE, row.names=FALSE)
}

create_ucsc_custom_bedGraph_track <- function(input_tb, output_bedGraph, browser_position=NA, name=NA, description=NA, visibility=NA, color=NA, altColor=NA, priority=NA, autoScale=NA, alwaysZero=NA, gridDefault=NA, maxHeightPixels=NA, graphType=NA, viewLimits=NA, yLineMark=NA, yLineOnOff=NA, windowingFunction=NA, smoothingWindow=NA) {
	# create header line
	header_opts <- c(name=name, description=description, visibility=visibility, color=color, altColor=altColor, priority=priority, autoScale=autoScale, alwaysZero=alwaysZero, gridDefault=gridDefault, maxHeightPixels=maxHeightPixels, graphType=graphType, viewLimits=viewLimits, yLineMark=yLineMark, yLineOnOff=yLineOnOff, windowingFunction=windowingFunction, smoothingWindow=smoothingWindow)
	header_opts <- header_opts[which(!is.na(header_opts))]
	header = paste0("track type=bedGraph ", paste(sapply(1:length(header_opts), function(i) {paste0(names(header_opts)[i], "=\"", header_opts[i], "\"")}), collapse=" "))

	# create browser position line
	if (!is.na(browser_position)) {
		browser_position = paste0("browser position ", browser_position)
		write(browser_position, output_bedGraph)
	}

	# save bedGraph
	if (!is.na(browser_position)) {
		write(browser_position, output_bedGraph)
		write(header, output_bedGraph, append=TRUE)
	} else {
		write(header, output_bedGraph)
	}
	write.table(input_tb, output_bedGraph, append=TRUE, quote=FALSE, sep="\t", col.names=FALSE, row.names=FALSE)
}

create_ucsc_custom_interact_track <- function(input_tb, output_interact, name, description, interactDirectional, maxHeightPixels, visibility, browser_position=NA) {
	# create header line
	header = paste0("track type=interact name=\"", name, "\" description=\"", description, "\" interactDirectional=", interactDirectional, "\" maxHeightPixels=", maxHeightPixels, "\" visibility=", visibility)

	# create browser position line
	if (!is.na(browser_position)) {
		browser_position = paste0("browser position ", browser_position)
		write(browser_position, output_interact)
	}

	# save interact
	if (!is.na(browser_position)) {
		write(browser_position, output_interact)
		write(header, output_interact, append=TRUE)
	} else {
		write(header, output_interact)
	}
	write.table(input_tb, output_interact, append=TRUE, quote=FALSE, sep="\t", col.names=FALSE, row.names=FALSE)
}

