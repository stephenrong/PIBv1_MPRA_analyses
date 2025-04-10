#!/bin/sh
# conda activate sequencing-qc
fastqc *_001.fastq.gz
multiqc	*_fastqc.zip
