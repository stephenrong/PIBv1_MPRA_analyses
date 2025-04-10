#!/bin/sh
gunzip -cd Tucci-MPRA-Enh-Barcode_S1_merged_R1_001.fastq.gz | head -n 400000 | bgzip >| Tucci-MPRA-Enh-Barcode_S1_merged_R1_001_test.fastq.gz
gunzip -cd Tucci-MPRA-Enh-Barcode_S1_merged_R2_001.fastq.gz | head -n 400000 | bgzip >| Tucci-MPRA-Enh-Barcode_S1_merged_R2_001_test.fastq.gz
