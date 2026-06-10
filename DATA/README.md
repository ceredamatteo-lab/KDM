# Data Folder
This folder contains the data used for the analyses reported in "**KDM: embedding DNA/RNA motifs and sequences in a shared k-mer space for unified discovery, analysis and binding prediction**". 

The sub folders are organized as follows:


## RBP-ENCODE and TF-ENCODE: main datasets
For the two extensive dataset used in the paper we report the bed files containing the regions in bed format.
The sub folders RBP-ENCODE and TF-ENCODE contain a metadata table reporting experiments information. 
For each experiment the BED sub folder contains two bed files: one for the training set and one for the test set.
 
Please note:

1. In the bed name field (fourth column) regions starting with "P" are positive sequences while those starting with "N" are the negative ones.
2. In the training sets the score field (fifth column) contain the group id used for **glmnet cross-validation**.

## RDATA: data in R format used by the scripts:
In this folder a number of Rdata files are collected and needed by the R scripts (see **CODE/R-SCRIPTS**):

ELENCARLI

## SUPPORT-FILES: other data needed to run the scripts
In this folder files in varius formats used by the scripts are collected (see **CODE/R-SCRIPTS**):

ELENCARLI


## **RESULTS**: empty folder, here the scripts will save the results.
Please note that the results produced by the scripts are quite bulking







