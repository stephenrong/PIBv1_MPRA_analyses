# Long-term isolation and archaic introgression shape functional genetic variation in Near Oceania

Code used for massively parallel reporter assay (MPRA) and functional genomics analyses. See the main repository for the manuscript: https://github.com/YourePrettyGood/PIBv1_manuscript/.

Contact the corresponding author with questions about these scripts. Last updated April 10th 2025. 

## MPRA analyses scripts here

- PIBv1_MPRA_pilot: Analyses based on the pilot archaic introgression dataset used to create the MPRA. Run PIBv1_MPRA_pilot/scripts/main_runall.sh in order.
- PIBv1_MPRA_final: Analyses based on the final archaic introgression dataset used in the manuscript. Run PIBv1_MPRA_final/scripts/main_runall.sh in order.
- PIBv1_MPRA_enrichr: Enrichr gene set enrichment input gene sets, output pathway enrichments, and scripts for figures.
- PIBv1_MPRA_collate: Scripts for collating MPRA and functional genomic analyses into intermediate files used for supplementary tables.
- PIBv1_MPRA_mpra: MRPA sequencing analysis pipeline. Download from Zenodo (see below).
- ../Datasets: Genomic or functional genomic annotations used in this project. Download from Zenodo (see below).

## Genomic and functional genomic annotations on Zenodo

Download Datasets.zip from Zenodo (...).

Unzip "Datasets" and place in same directory as PIBv1_MPRA_analyses.

## MPRA sequencing analysis on Zenodo

Download PIBv1_MPRA_mpra.zip from Zenodo (...).

Unzip "PIBv1_MPRA_mpra" and replace subdirectory PIBv1_MPRA_analyses/PIBv1_MPRA_mpra.

Cd into PIBv1_MPRA_mpra and git clone https://github.com/tewhey-lab/MPRAmodel and https://github.com/tewhey-lab/MPRA_oligo_barcode_pipeline repositories and follow install dependencies instructions therein.

## Software dependencies
Analyses performed on Yale University HPC as Slurm scripts or in RStudio.

### Software versions
- fastqc=0.11.9
- multiqc=v1.10.1
- htslib=1.21
- vcftools=0.1.16
- bedtools=v2.31.1
- crossmap=0.6.4
- ensembl-vep=105.0
- ucsc-bedtobigbed=366
- ucsc-bigwigaverageoverbed=366
- R=4.4.1 (2024-06-14)
- Rstudio=2024.09.0+375

### R sessionInfo()
```
R version 4.4.1 (2024-06-14)
Platform: aarch64-apple-darwin20
Running under: macOS 15.3.2

Matrix products: default
BLAS:   /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/lib/libRblas.0.dylib 
LAPACK: /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.0

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

time zone: Europe/London
tzcode source: internal

attached base packages:
[1] grid      stats4    stats     graphics  grDevices utils     datasets 
[8] methods   base     

other attached packages:
 [1] knitr_1.48                        biomaRt_2.60.1                   
 [3] AnnotationHub_3.12.0              BiocFileCache_2.12.0             
 [5] dbplyr_2.5.0                      wrapr_2.1.0                      
 [7] viridis_0.6.5                     viridisLite_0.4.2                
 [9] vcfR_1.15.0                       lubridate_1.9.3                  
[11] forcats_1.0.0                     stringr_1.5.1                    
[13] dplyr_1.1.4                       purrr_1.0.2                      
[15] readr_2.1.5                       tidyr_1.3.1                      
[17] tibble_3.2.1                      tidyverse_2.0.0                  
[19] scales_1.3.0                      plyranges_1.24.0                 
[21] paletteer_1.6.0                   optparse_1.7.5                   
[23] motifbreakR_2.18.0                MotifDb_1.46.0                   
[25] hrbrthemes_0.8.7                  ggrepel_0.9.6                    
[27] ggpubr_0.6.0                      ggbio_1.52.0                     
[29] ggplot2_3.5.1                     genekitr_1.2.8                   
[31] eulerr_7.0.2                      data.table_1.16.0                
[33] RColorBrewer_1.1-3                BSgenome.Hsapiens.UCSC.hg19_1.4.3
[35] BSgenome_1.72.0                   rtracklayer_1.64.0               
[37] BiocIO_1.14.0                     Biostrings_2.72.1                
[39] XVector_0.44.0                    GenomicRanges_1.56.1             
[41] GenomeInfoDb_1.40.1               IRanges_2.38.1                   
[43] S4Vectors_0.42.1                  BiocGenerics_0.50.0              

loaded via a namespace (and not attached):
  [1] R.methodsS3_1.8.2           dichromat_2.0-0.1          
  [3] progress_1.2.3              urlchecker_1.0.1           
  [5] nnet_7.3-19                 poweRlaw_0.80.0            
  [7] vctrs_0.6.5                 digest_0.6.37              
  [9] png_0.1-8                   deldir_2.0-4               
 [11] permute_0.9-7               MASS_7.3-61                
 [13] fontLiberation_0.1.0        reshape2_1.4.4             
 [15] httpuv_1.6.15               qvalue_2.36.0              
 [17] withr_3.0.1                 xfun_0.48                  
 [19] ggfun_0.1.6                 ellipsis_0.3.2             
 [21] memoise_2.0.1               clusterProfiler_4.12.6     
 [23] gson_0.1.0                  profvis_0.4.0              
 [25] systemfonts_1.1.0           tidytree_0.4.6             
 [27] gtools_3.9.5                R.oo_1.26.0                
 [29] Formula_1.2-5               GGally_2.2.1               
 [31] prettyunits_1.2.0           rematch2_2.1.2             
 [33] KEGGREST_1.44.1             promises_1.3.0             
 [35] pinfsc50_1.3.0              httr_1.4.7                 
 [37] rstatix_0.7.2               restfulr_0.0.15            
 [39] rstudioapi_0.16.0           UCSC.utils_1.0.0           
 [41] miniUI_0.1.1.1              generics_0.1.3             
 [43] DOSE_3.30.5                 base64enc_0.1-3            
 [45] curl_5.2.3                  zlibbioc_1.50.0            
 [47] ggraph_2.2.1                polyclip_1.10-7            
 [49] GenomeInfoDbData_1.2.12     SparseArray_1.4.8          
 [51] RBGL_1.80.0                 xtable_1.8-4               
 [53] ade4_1.7-22                 pracma_2.4.4               
 [55] evaluate_1.0.0              S4Arrays_1.4.1             
 [57] hms_1.1.3                   colorspace_2.1-1           
 [59] filelock_1.0.3              getopt_1.20.4              
 [61] magrittr_2.0.3              later_1.3.2                
 [63] ggtree_3.12.0               lattice_0.22-6             
 [65] XML_3.99-0.17               triebeard_0.4.1            
 [67] shadowtext_0.1.4            cowplot_1.1.3              
 [69] matrixStats_1.4.1           Hmisc_5.1-3                
 [71] pillar_1.9.0                nlme_3.1-166               
 [73] pwalign_1.0.0               caTools_1.18.3             
 [75] compiler_4.4.1              stringi_1.8.4              
 [77] SummarizedExperiment_1.34.0 devtools_2.4.5             
 [79] GenomicAlignments_1.40.0    plyr_1.8.9                 
 [81] crayon_1.5.3                abind_1.4-8                
 [83] gridGraphics_0.5-1          graphlayouts_1.2.0         
 [85] bit_4.5.0                   geneset_0.2.7              
 [87] fastmatch_1.1-4             codetools_0.2-20           
 [89] biovizBase_1.52.0           mime_0.12                  
 [91] splines_4.4.1               Rcpp_1.0.13                
 [93] europepmc_0.4.3             Rttf2pt1_1.3.12            
 [95] interp_1.1-6                blob_1.2.4                 
 [97] utf8_1.2.4                  BiocVersion_3.19.1         
 [99] seqLogo_1.70.0              AnnotationFilter_1.28.0    
[101] fs_1.6.4                    checkmate_2.3.2            
[103] pkgbuild_1.4.4              Gviz_1.48.0                
[105] openxlsx_4.2.7.1            ggsignif_0.6.4             
[107] ggplotify_0.1.2             Matrix_1.7-0               
[109] tzdb_0.4.0                  tweenr_2.0.3               
[111] pkgconfig_2.0.3             tools_4.4.1                
[113] cachem_1.1.0                RSQLite_2.3.7              
[115] DBI_1.2.3                   splitstackshape_1.4.8      
[117] fastmap_1.2.0               rmarkdown_2.28             
[119] usethis_3.0.0               Rsamtools_2.20.0           
[121] broom_1.0.7                 patchwork_1.3.0            
[123] BiocManager_1.30.25         ggstats_0.7.0              
[125] VariantAnnotation_1.50.0    graph_1.82.0               
[127] carData_3.0-5               rpart_4.1.23               
[129] farver_2.1.2                mgcv_1.9-1                 
[131] tidygraph_1.3.1             scatterpie_0.2.4           
[133] yaml_2.3.10                 latticeExtra_0.6-30        
[135] MatrixGenerics_1.16.0       foreign_0.8-87             
[137] cli_3.6.3                   txdbmaker_1.0.1            
[139] lifecycle_1.0.4             Biobase_2.64.0             
[141] sessioninfo_1.2.2           backports_1.5.0            
[143] BiocParallel_1.38.0         annotate_1.82.0            
[145] timechange_0.3.0            gtable_0.3.5               
[147] rjson_0.2.23                parallel_4.4.1             
[149] ape_5.8                     jsonlite_1.8.9             
[151] TFBSTools_1.42.0            bitops_1.0-9               
[153] bit64_4.5.2                 vegan_2.6-8                
[155] yulab.utils_0.1.7           zip_2.3.1                  
[157] urltools_1.7.3              CNEr_1.40.0                
[159] GOSemSim_2.30.2             R.utils_2.12.3             
[161] lazyeval_0.2.2              shiny_1.9.1                
[163] htmltools_0.5.8.1           enrichplot_1.24.4          
[165] GO.db_3.19.1                rappdirs_0.3.3             
[167] ensembldb_2.28.1            glue_1.8.0                 
[169] TFMPvalue_0.0.9             ggvenn_0.1.10              
[171] httr2_1.0.5                 gdtools_0.4.0              
[173] RCurl_1.98-1.16             treeio_1.28.0              
[175] jpeg_0.1-10                 motifStack_1.48.0          
[177] gridExtra_2.3               igraph_2.0.3               
[179] extrafontdb_1.0             R6_2.5.1                   
[181] GenomicFeatures_1.56.0      cluster_2.1.6              
[183] pkgload_1.4.0               aplot_0.2.3                
[185] DirichletMultinomial_1.46.0 DelayedArray_0.30.1        
[187] tidyselect_1.2.1            ProtGenerics_1.36.0        
[189] htmlTable_2.4.3             ggforce_0.4.2              
[191] xml2_1.3.6                  fontBitstreamVera_0.1.1    
[193] car_3.1-3                   AnnotationDbi_1.66.0       
[195] munsell_0.5.1               fontquiver_0.2.1           
[197] htmlwidgets_1.6.4           fgsea_1.30.0               
[199] rlang_1.1.4                 extrafont_0.19             
[201] remotes_2.5.0               fansi_1.0.6                
[203] OrganismDbi_1.46.0 
```
