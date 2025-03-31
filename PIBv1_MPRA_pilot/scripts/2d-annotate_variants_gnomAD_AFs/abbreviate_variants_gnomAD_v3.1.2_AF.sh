bcftools#!/bin/sh

INPUT="../../results/2d-annotate_variants_gnomAD_AFs/adaptive_variants_lift38_chr_lift37-gnomad.genomes.v3.1.2.sites.vcf.gz"
OUTPUT="../../results/2d-annotate_variants_gnomAD_AFs/adaptive_variants_lift38_chr_lift37_abbrev-gnomad.genomes.v3.1.2.sites.vcf.gz"

bcftools annotate -x "^INFO/AF,^INFO/AC,^INFO/AN,^INFO/AF_oth,^INFO/AC_oth,^INFO/AN_oth,^INFO/AF_ami,^INFO/AC_ami,^INFO/AN_ami,^INFO/AF_sas,^INFO/AC_sas,^INFO/AN_sas,^INFO/AF_fin,^INFO/AC_fin,^INFO/AN_fin,^INFO/AF_eas,^INFO/AC_eas,^INFO/AN_eas,^INFO/AF_amr,^INFO/AC_amr,^INFO/AN_amr,^INFO/AF_afr,^INFO/AC_afr,^INFO/AN_afr,^INFO/AF_mid,^INFO/AC_mid,^INFO/AN_mid,^INFO/AF_asj,^INFO/AC_asj,^INFO/AN_asj,^INFO/AF_nfe,^INFO/AC_nfe,^INFO/AN_nfe" ${INPUT} | bgzip >| ${OUTPUT}
tabix -p vcf ${OUTPUT}
