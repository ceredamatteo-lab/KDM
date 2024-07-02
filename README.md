# KDM

Nella cartella Scripts/ ho aggiunto i file di R per le prime analisi di KDM sui datasets di GraphProt (https://doi.org/10.1186/gb-2014-15-1-r17)
01. **create_CV.R** : per ogni esperimento nel dataset, crea un Rdata con le sequenze positive + negative divise nei 10 sottogruppi necessati per la cross-validation
02. **windows_CV_new.R** : dato un esperimento e una specifica combinazione di kmer length + numero di gaps + numero di motivi, calcola la performance usando una 10-fold cv
03. **generate_data_motifs.R** : per i dataset di graphprot, partendo dal modello kdm ottenuto con una specifica combinazione di parametri (6-2-70) crea gli Rdata necessari per visualizzare le proprietá di ogni motivo
    - enrichement
    - centrality
    - analisi dei coefficenti nel modello _GLM_
    - motivo convertito in PWM
    - confronto con la reference di _mCross_
