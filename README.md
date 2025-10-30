episcope — Epigenomic Integration and Cell-State–Specific Gene Regulatory Network Analysis
=========================================================================================

[![Version](https://img.shields.io/badge/version-0.9.0-blue.svg?style=plastic)](https://github.com/yaoxiangli/episcope)
[![License](https://img.shields.io/badge/license-MIT-green.svg?style=plastic)](LICENSE)
[![Bioconductor](https://img.shields.io/badge/install%20via-BiocManager-orange.svg?style=plastic)](https://bioconductor.org)
[![Maintenance](https://img.shields.io/badge/Maintained%3F-yes-green.svg?style=plastic)](https://github.com/yaoxiangli/episcope/graphs/commit-activity)
[![publication](https://img.shields.io/badge/Publication-in%20prep-blue.svg?style=plastic)]()

---

## Introduction

**episcope** is an R/Bioconductor-style framework for *integrative, reproducible, and cell-state–specific* analysis of chromatin accessibility and transcriptional regulation.  
It bridges ATAC-seq and RNA-seq data to infer **gene regulatory networks (GRNs)** that are **cell-type– and condition-dependent**, with built-in visualization, delta analysis, and network clustering modules.

<img align="right" width=150 src="man/figures/episcope_logo.png">

episcope provides:

- Unified pipelines for **ATAC-seq**, **RNA-seq**, and **multi-condition delta analyses**
- Integration of **GeneHancer**, **JASPAR**, and **HOCOMOCO** regulatory priors
- Automated construction of **TF→enhancer→gene** regulatory graphs
- Interactive **visNetwork**-based TF-centric and topic-centric subnetworks
- LDA and Louvain clustering of differential GRNs to identify shared regulatory programs
- Fully reproducible R package infrastructure compatible with Bioconductor workflows

While developed for **pancreatic ductal adenocarcinoma (PDAC)** nutrient-stress and EMT studies, episcope is broadly applicable to any context combining ATAC-seq and RNA-seq across multiple perturbations or cell types.

---

## Installation

episcope is written as an R package and can be installed directly from GitHub:

```r
# Using remotes
remotes::install_github("yaoxiangli/episcope")

# or using pak
pak::pak("yaoxiangli/episcope")
```

To install all dependencies (recommended for full functionality):
```r
install.packages(c("visNetwork", "igraph", "ggplot2", "data.table", "BiocManager"))
BiocManager::install(c("DESeq2", "GenomicRanges", "SummarizedExperiment"))
```

---

## Usage Overview

All episcope functions are available directly within R.  
The main workflow follows three stages:

1. **Prepare inputs** — processed ATAC and RNA tables, plus enhancer–gene priors (e.g., GeneHancer ELITE)
2. **Integrate and model** — compute differential accessibility and expression, correlate features, and build GRNs
3. **Visualize and interpret** — render interactive networks or topic-based subnetworks for exploratory analysis

Example:

```r
library(episcope)

# Load ATAC and RNA tables
atac <- episcope::load_links_table("AsPC1_FBS_vs_0pctFBS_links.csv")
rna  <- episcope::load_rna_table("AsPC1_FBS_vs_0pctFBS_rna.csv")

# Build condition-specific GRN
grn <- episcope::compare_links_two_conditions(atac, rna)

# Render TF hub (bi-directional)
episcope::render_tf_hub_network(
  comp_csv = "AsPC1_FBS_vs_0pctFBS_delta_links_filtered_lda_K20.csv",
  input_tf = "CEBPB"
)

# Or Δ-links network
episcope::render_tf_hub_delta_network(
  comp_csv = "AsPC1_FBS_vs_0pctFBS_delta_links_filtered_lda_K20.csv",
  input_tf = "CEBPB"
)
```

---

## Major Modules

* [utils_grn_diff.R](R/utils_grn_diff.R): Compute link deltas across conditions  
* [utils_grn_filter.R](R/utils_grn_filter.R): Filter and normalize link tables  
* [utils_grn_link_network_tf.R](R/utils_grn_link_network_tf.R): Render TF-centric subnetworks  
* [utils_grn_link_network_topic.R](R/utils_grn_link_network_topic.R): Render topic-based subnetworks (LDA clusters)  
* [utils_grn_lda.R](R/utils_grn_lda.R): Topic modeling of differential link matrices  
* [utils_tobias_overview_loader_simple.R](R/utils_tobias_overview_loader_simple.R): Simplified TOBIAS output integration  

---

## Example Workflow

The episcope pipeline typically runs as:

1. **Load preprocessed ATAC and RNA data**
2. **Compute link deltas** with `compare_links_bulk()`
3. **Integrate GeneHancer priors** (TF–enhancer–gene)
4. **Filter edges** by correlation and motif evidence
5. **Cluster GRNs** with Louvain or LDA
6. **Visualize** TF- or topic-centric subnetworks interactively

```
ATAC + RNA  →  Δ-link tables  →  GRN (TF–enhancer–gene)
                           ↘
                       visualization & topic analysis
```

---

## Pipelines

episcope modules can be run independently or as part of a unified analysis pipeline.

**Snakemake pipeline (HPC compatible)**  
Pre-set workflow for bias-correction, footprinting, GRN inference, and visualization:  
👉 [episcope_snakemake](https://github.com/yaoxiangli/episcope_snakemake)

**Nextflow pipeline (cloud-ready)**  
Nextflow version for scalable parallel execution:  
👉 [episcope_nextflow](https://github.com/yaoxiangli/episcope_nextflow)

---

## Help and Documentation

- 📘 [Wiki and Tutorials](https://github.com/yaoxiangli/episcope/wiki)
- ❓ [FAQ](https://github.com/yaoxiangli/episcope/wiki/FAQ)
- 🐛 [Report an Issue](https://github.com/yaoxiangli/episcope/issues)

---

## How to Cite

Li, Y. *et al.* (in preparation).  
**episcope:** A reproducible framework for integrative chromatin and transcriptomic network analysis in PDAC.

---

## License

This project is licensed under the [MIT License](LICENSE).
