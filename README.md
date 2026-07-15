# KDM

[![DOI](https://img.shields.io/badge/DOI-10.5281%2Fzenodo.21383381-blue)](https://doi.org/10.5281/zenodo.21383381)

**Embedding DNA/RNA motifs and sequences in a shared k-mer space for unified discovery, analysis and binding prediction**

This repository contains the analysis code, precomputed intermediate data, and figures accompanying the KDM paper:

> Fumagalli L.\*, Becchi T.\*, Cereda M.†, Pozzoli U.† — *KDM: embedding DNA/RNA motifs and sequences in a shared k-mer space for unified discovery, analysis and binding prediction.* \*Equal contribution. †Joint senior and corresponding authors.

------------------------------------------------------------------------

## Overview

Motif discovery and binding-site prediction in DNA and RNA are central tasks in regulatory genomics, but the methodological landscape is split between **interpretable but rigid** position weight matrices (PWMs) and **high-performing but opaque** machine-learning models.

**KDM** is a unifying framework in which *both* motifs and sequences are represented as probability distributions over a shared **k-mer dictionary**, embedded via the **Hellinger transformation**. This common geometry enables motif–sequence scoring, motif–motif comparison, de novo discovery, and binding prediction with a single primitive: the **Bhattacharyya coefficient**.

Four tools are instantiated on this representation:

| Tool | Purpose |
|-----------------------------|-------------------------------------------|
| **KDMMap** | Positional enrichment analysis |
| **KDMMatch** | Information-content-aware motif matching / annotation |
| **KDMFind** | Unsupervised de novo motif discovery via projective non-negative matrix factorization (PNMF) |
| **KDM-LRLM** | Binding prediction with Lasso-regularized logistic regression |

------------------------------------------------------------------------

## Repository scope

> **Note.** This repository is the **reproducibility / figure-generation companion** to the paper. It holds the R scripts that produce the manuscript figures together with the precomputed `Rdata` they consume. The core KDM engine itself lives in the separate R package **`rkdMotifs`**, which the discovery/analysis scripts depend on (see [Dependencies](#dependencies)).

```         
KDM/
├── CODE/
   ├── R-SCRIPTS   # Scripts used to produce paper figure
   ├── R-PACKAGE   # Code for installing rkdmotifs     
├── DATA/          # Datasets used to produce paper figure + Rdata generated
├── gpl-3.0.md     # GPL license
└── README.md
```

## Citation

If you use KDM, please cite the paper above. The KDM software package is available at <https://github.com/ceredamatteo-lab/KDM>.

## License

© 2010-2026 Uberto Pozzoli and Matteo Cereda

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

## Contact

- Uberto Pozzoli — [uberto.pozzoli\@lanostrafamiglia.it](mailto:uberto.pozzoli@lanostrafamiglia.it){.email}
- Matteo Cereda — [matteo.cereda\@ifom.eu](mailto:matteo.cereda@ifom.eu){.email}
