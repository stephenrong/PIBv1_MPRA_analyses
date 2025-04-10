#!/bin/sh
cat Tucci-MPRA-Enh-Barcode_S1_L001_R1_001.fastq.gz Tucci-MPRA-Enh-Barcode_S1_L002_R1_001.fastq.gz Tucci-MPRA-Enh-Barcode_S1_L003_R1_001.fastq.gz Tucci-MPRA-Enh-Barcode_S1_L004_R1_001.fastq.gz >| Tucci-MPRA-Enh-Barcode_S1_merged_R1_001.fastq.gz
cat Tucci-MPRA-Enh-Barcode_S1_L001_R2_001.fastq.gz Tucci-MPRA-Enh-Barcode_S1_L002_R2_001.fastq.gz Tucci-MPRA-Enh-Barcode_S1_L003_R2_001.fastq.gz Tucci-MPRA-Enh-Barcode_S1_L004_R2_001.fastq.gz >| Tucci-MPRA-Enh-Barcode_S1_merged_R2_001.fastq.gz
