# HCI011_RNAseq

This repository contains the data and code for the RNA-seq analysis done in the paper __FOXM1 is a biomarker of resistance to PI3Kα inhibition in ER+ breast cancer that is detectable using metabolic imaging__ (Ros et al. 2020).

Input data includes processed gene counts and sample metadata information, and the gene lists used in the analysis. 

## Requirements
R (>= 3.61) with these packages:

* plyr
* dplyr
* tibble
* data.table
* RColorBrewer
* pheatmap
* svglite (only required if wishing to produce svg figures)
* biomaRt
* org.Hs.eg.db
* edgeR
* scales
* GSEABase
* limma

## Usage

To reproduce the analysis, simply switch to the `HCI011_RNAseq/R` directory, open an R session and type `source("PDX_FOXM1_RNAseq_pub.R")`

Figures and tables will be generated under the `HCI011_RNAseq/output` directory. 
