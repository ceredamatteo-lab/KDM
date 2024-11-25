Each RBP has a dedicated folder with 4 Rdata objects for each cell_line (K562 and HepG2)
01. **RBPname_CellLine_Genome_exons_info.rds**: Table with exons classification (Enhanced-Silenced-Control) and exons information (exons length and introns length)
02. **RBPname_CellLine_Genome_W.rds** : MRMDs matrix with the 70 MRMs distributions obtained after the factorization step
03. **RBPname_CellLine_Genome_ridge.rds** : this Rdata objects is a list that contains multiple object genreated through the analysis
  3.1 *coeff* : a dataframe with the glm coefficents of each MRMD in each subregion   
  3.2 *sel* : MRMDs with at least one subregion in which the overrepresentation is significant
  3.3 *profile* : score of each MRMD in each position of each subregion on both altered and control regions
  3.4 *test*: overrepresentation results of each MRMD in each subregion
  3.5 *pwm*s: MRMD to PWM convertion of each MRMD
  3.6 *kdmtom*: results for the comparison betweeen each MRMDs and mCross PWMs converted into MRMD
04. **RBPname_CellLine_Genome_recall_1000.rds** : a dataframe with the result obtained on the downsampling analysis. This dataframe contains a column *p* with the percentage of total exons used for the downsampling, a column *n* with the fraction of original MRMDs that are retreived in the W matix obtained after downsamling and the *iter* column that contains iterations steps from 1 to 1000
