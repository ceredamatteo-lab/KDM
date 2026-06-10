# KDM

**Embedding DNA/RNA motifs and sequences in a shared k-mer space for unified discovery, analysis and binding prediction**

This repository contains the analysis code, precomputed intermediate data, and figures accompanying the KDM paper:

> Fumagalli L.\*, Becchi T.\*, Cereda M.†, Pozzoli U.† — *KDM: embedding DNA/RNA motifs and sequences in a shared k-mer space for unified discovery, analysis and binding prediction.*
> \*Equal contribution. †Joint senior and corresponding authors.

---

## Overview

Motif discovery and binding-site prediction in DNA and RNA are central tasks in regulatory genomics, but the methodological landscape is split between **interpretable but rigid** position weight matrices (PWMs) and **high-performing but opaque** machine-learning models.

**KDM** is a unifying framework in which *both* motifs and sequences are represented as probability distributions over a shared **k-mer dictionary**, embedded via the **Hellinger transformation**. This common geometry enables motif–sequence scoring, motif–motif comparison, de novo discovery, and binding prediction with a single primitive: the **Bhattacharyya coefficient**.

Four tools are instantiated on this representation:

| Tool | Purpose |
|------|---------|
| **KDMMap** | Positional enrichment analysis |
| **KDMMatch** | Information-content-aware motif matching / annotation |
| **KDMFind** | Unsupervised de novo motif discovery via projective non-negative matrix factorization (PNMF) |
| **KDM-LRLM** | Binding prediction with Lasso-regularized logistic regression |

---

## Repository scope

> **Note.** This repository is the **reproducibility / figure-generation companion** to the paper. It holds the R scripts that produce the manuscript figures together with the precomputed `Rdata` they consume. The core KDM engine itself lives in the separate R package **`rkdMotifs`**, which the discovery/analysis scripts depend on (see [Dependencies](#dependencies)).

```
KDM/
├── Scripts/        # 25 R scripts: analysis pipelines + figure generation
├── Rdata/          # 42 precomputed .rds / .Rdata intermediate objects
├── Tables/         # Supplementary metadata tables (TSV)
├── Figure/         # 109 generated PDF figures, grouped by dataset/panel
├── LICENSE         # MIT
└── README.md
```

### `Figure/` layout

| Subfolder | Content |
|-----------|---------|
| `Paper_Figure/` | Main manuscript panels (P4–P22, IC plots, correlation distributions) for TF and RBP |
| `ENCODE_Dataset/` | ENCODE eCLIP analyses: binding AUCs, classification, KDMM panels, Jaccard agreement |
| `DNA/` | DNA (TF) classification and KDMM panels |
| `RBP_24/` | Curated 24-RBP benchmark: parameter optimization, GraphProt/RNAProt AUC comparisons, ANNOVAR heatmaps |
| `Splicing_Maps/` | Per-RBP splicing maps (e.g. HNRNPK, U2AF2) |
| `Figure0/` | Introductory / schematic examples |

---

## Dependencies

### Core engine

- **[`rkdMotifs`](https://github.com/ceredamatteo-lab/KDM)** — the KDM R package providing the k-mer / Hellinger representation and the `kdm*` functions (`kdmCentrimo`, `kdmTom`, `kdmGetSequence`, `kdmFeaturesProfileFromPWMSet`, `kdmLoadPWMSet`, …) and the `kdmotifs` S4 class. Required by all *discovery/analysis* scripts (`figureA–E.R`, `P*_RBP.R`, `p15.R`).

### R packages (CRAN / Bioconductor)

```r
# CRAN
install.packages(c(
  "ggplot2", "ggpubr", "dplyr", "reshape2", "scales",
  "RColorBrewer", "viridis", "circlize", "irr", "jaccard", "ggstatsplot"
))

# Bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("universalmotif", "ComplexHeatmap"))
```

---

## Usage

### Reproducing figures from precomputed data (recommended)

Most plotting scripts read intermediate objects already shipped in `Rdata/` and write PDFs into `Figure/`. They use **paths relative to the repository root**, so run R from the top of the repo:

```bash
cd /path/to/KDM
Rscript Scripts/ENCODE_Dataset_AUC.R          # binding-prediction AUC panels
Rscript Scripts/ENCODE_Dataset_Classification.R
Rscript Scripts/RBP24_Dataset_AUC.R           # 24-RBP benchmark panels
Rscript Scripts/ENCODE_KDMM.R                 # KDMM example panels (PRPF8, U2AF2, …)
Rscript Scripts/IC_PWM.R                      # information-content comparisons
```

Several scripts are split into a **compute** stage (calls `rkdMotifs`, regenerates the `.rds`/`.Rdata` and is usually commented out / guarded) and a **plot** stage (reads the saved object). To only regenerate figures, run the plot stage.

### Re-running the full analysis from raw data

The discovery scripts (`figureA.R`–`figureE.R`, `P4*–P22*_RBP.R`, `p15.R`) reference upstream inputs that are **not** bundled here — e.g. `DATASET.Rdata`, `mCrossBase.Rdata`, `ENCODE_eCLIP_DATASET.Rdata`, genome 2bit files, and an internal helper (`funzioni_importanti.rtx`). These contain hard-coded absolute paths from the authors' environment and will need to be adapted to your local layout before they can be executed end-to-end.

---

## Script index

| Script | What it produces |
|--------|------------------|
| `00_plot_functions.R` | Shared plotting/statistics helpers (binomial enrichment tests, etc.) |
| `figureA.R` – `figureE.R` | Core method figures (KDMMap/KDMMatch behaviour, IC-aware scoring, agreement with CentriMo/Tomtom) |
| `IC_PWM.R` | Information-content comparison of DNA (HOCOMOCO) vs RNA (mCrossBase) PWMs |
| `ENCODE_KDMM.R` | KDMM positional-decomposition panels for example experiments |
| `ENCODE_Dataset_AUC.R` | KDM-LRLM binding AUCs; competitor correlation; per-cluster / per-cell breakdowns |
| `ENCODE_Dataset_Classification.R`, `*_new.R` | Motif-annotation classification vs mCross / Tomtom / reference; Jaccard agreement |
| `ENCODE_annovar.R` | ENCODE eCLIP ANNOVAR genomic-region enrichment heatmap |
| `RBP24_Dataset_AUC.R` | 24-RBP benchmark: parameter optimization, GraphProt / RNAProt AUC comparison |
| `HEATMAP_annovar.R`, `Heatmap_KDMM_annovar.R` | ANNOVAR region-enrichment heatmaps across methods/experiments |
| `P4A_P9_RBP.R`, `P4B_RBP.R`, `P6_RBP.R`, `P8_RBP.R`, `P10_RBP.R`, `P15_RBP.R`/`p15.R`, `P16_P17_P18_RBP.R`, `P19_RBP.R`, `P22_RBP.R` | Individual manuscript panels (RBP & TF): discovery, comparison, and prediction analyses |

---

## Data

- `Tables/RBP_24_Info.tsv` — metadata for the curated 24-RBP benchmark (experiment, source paper, CLIP method, cell line, protein, function, genomic-region cluster, peak count).
- `Rdata/` — precomputed objects grouped by analysis: `ENCODE_Dataset/`, `RBP_24/`, `Paper_Figure/`, plus top-level `KDMM_*.rds` example decompositions.

---

## Citation

If you use KDM, please cite the paper above. The KDM software package is available at
<https://github.com/ceredamatteo-lab/KDM>.

## License

Released under the [MIT License](LICENSE). © 2024 Matteo Cereda.

## Contact

- Uberto Pozzoli — <uberto.pozzoli@lanostrafamiglia.it>
- Matteo Cereda — <matteo.cereda@ifom.eu>
