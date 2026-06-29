# R Scripts

This folder contains all the scripts used to produce the results reported in "**KDM: embedding DNA/RNA motifs and sequences in a shared k-mer space for unified discovery, analysis and binding prediction**".

## Dependencies

First, download the `hg38.2bit` file, which is required to extract genomic sequences from BED files:

`cd KDM/DATA/SUPPORT_FILES/`

`wget https://hgdownload.gi.ucsc.edu/goldenPath/hg38/bigZips/hg38.2bit`

All the scripts require our package **rkdMotifs** to be installed (see folder CODE/R-PACKAGE)

The scripts have a number of other dependencies:

### R packages:

##### From CRAN

- **ROCR**
- **glmnet**
- **dplyr**
- **ggplot2**
- **reshape2**
- **RColorBrewer**
- **ggpubr**
- **doMC**
- **stringr**
- **irr**
- **tidyr**

##### **Fro**m Bioconductor

- **universalmotif**

##### From other repositories

- **SeqGL** (<https://github.com/ManuSetty/SeqGL>)

### Other required non-R tools which need to be installed to execute the scripts:

- **GraphProt** (<https://github.com/dmaticzka/GraphProt>)
- **MEME suite (version 5.4.0)** (<https://meme-suite.org/meme/doc/download.html>)
- **ls-GKM** (<https://github.com/Dongwon-Lee/lsgkm/>)
- **RNAProt** (<https://github.com/BackofenLab/RNAProt>)

## SCRIPTS DESCRIPTION

Once you have downloaded To execute the scripts you need to set your R working directory to the CODE folder of this repository

## R VERSION

All the scripts have been tested on R version 4.0.5 (2021-03-31) -- "Shake and Throw"
