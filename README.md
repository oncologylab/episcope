episcope — Integrative Multi-Omics Framework for Condition-Specific Gene Regulatory Networks
===========================================================================================

[![Version](https://img.shields.io/badge/version-1.0.0-blue.svg?style=plastic)](https://github.com/oncologylab/episcope)
[![License](https://img.shields.io/badge/license-MIT-green.svg?style=plastic)](https://github.com/oncologylab/episcope/blob/main/LICENSE)
[![Bioconductor](https://img.shields.io/badge/install%20via-BiocManager-orange.svg?style=plastic)](https://bioconductor.org)
[![Maintenance](https://img.shields.io/badge/Maintained%3F-yes-green.svg?style=plastic)](https://github.com/oncologylab/episcope/graphs/commit-activity)
[![publication](https://img.shields.io/badge/Publication-in%20prep-blue.svg?style=plastic)]()

---

## Introduction

**episcope** is a modular, reproducible framework for building and analyzing **condition-specific gene regulatory networks (GRNs)** through the integration of **chromatin accessibility**, **footprinting**, and **transcriptomic** data.

It extends the ideas of ATAC-seq footprinting (as in TOBIAS) to full multi-omics GRN modeling, enabling users to:

- Load ATAC, footprint, and RNA-seq data from any source  
- Align and correct footprint signals across motifs  
- Correlate TF binding and expression to refine assignments  
- Build regulation priors using GeneHancer, TSS proximity, or Hi-C data  
- Infer data-specific TF→enhancer→gene networks  
- Compare and visualize GRNs across multiple conditions or perturbations  

episcope is data- and cell-type–agnostic: it supports **any organism**, **treatment**, or **perturbation**, including large-scale matched ATAC–RNA datasets.

---

## Installation

episcope is written in R and available on GitHub.

```r
# Using remotes
remotes::install_github("oncologylab/episcope")

# or using pak
pak::pak("oncologylab/episcope")
```

Dependencies (Bioconductor and CRAN):
```r
install.packages(c("visNetwork", "igraph", "ggplot2", "data.table", "BiocManager"))
BiocManager::install(c("DESeq2", "GenomicRanges", "SummarizedExperiment"))
```

---

## Usage Overview

episcope modules are fully interoperable and can be executed independently or as part of a full pipeline.

```r
library(episcope)

# Example: build GRN from matched ATAC–RNA dataset
fp <- load_footprints("footprints_corrected.bw", "peaks.bed")
rna <- load_rna("rna_expression.csv")

fp_aligned <- align_footprints(fp)
fp_corrected <- correct_footprints(fp_aligned)
tf_map <- map_tf_to_footprints(fp_corrected, rna)

priors <- build_regulation_priors(method = "genehancer")
refined1 <- correlate_atac_to_gene(priors, fp, rna)
refined2 <- correlate_fp_to_gene(refined1)

grn <- build_basal_grn(refined2)
validate_grn_perturbation(grn, perturb_db = "perturbdb.sqlite")

light_condition_grn(grn, condition = "GlcLow")
compare_grn(grn, conditionA = "Ctrl", conditionB = "Stress")
```

---

## Modules Overview

episcope provides modular functions for each step, from raw footprint data to condition-specific regulatory networks.

### Data Loading and Preprocessing
- `load_footprints()` imports footprint score tracks and ATAC peaks.  
- `load_rna()` loads RNA-seq quantifications.  
- `align_footprints()` consolidates motif-redundant footprint calls into consensus sites.  
- `correct_footprints()` applies bias correction and depth normalization.  
  → [https://github.com/oncologylab/episcope/wiki/load_footprints](https://github.com/oncologylab/episcope/wiki/load_footprints)

### TF Assignment and Regulation Priors
- `map_tf_to_footprints()` correlates footprint signals to TF expression to assign regulators.  
- `build_regulation_priors()` defines TF–gene edges from GeneHancer, ±30 kb TSS, or Hi-C contact maps.  
  → [https://github.com/oncologylab/episcope/wiki/map_tf_to_footprints](https://github.com/oncologylab/episcope/wiki/map_tf_to_footprints)

### Correlation-Based Refinement
- `correlate_atac_to_gene()` refines enhancer–gene links by ATAC–RNA correlation.  
- `correlate_fp_to_gene()` quantifies TF activity by correlating footprint scores with gene expression.  
  → [https://github.com/oncologylab/episcope/wiki/correlate_fp_to_gene](https://github.com/oncologylab/episcope/wiki/correlate_fp_to_gene)

### Network Assembly and Validation
- `build_basal_grn()` constructs a dataset-specific GRN.  
- `validate_grn_perturbation()` validates TF–gene links using knockout/knockdown datasets or perturbation databases.  
  → [https://github.com/oncologylab/episcope/wiki/build_basal_grn](https://github.com/oncologylab/episcope/wiki/build_basal_grn)

### Condition-Specific and Differential GRNs
- `light_condition_grn()` identifies active regulatory edges within each condition.  
- `compare_grn()` computes differential networks between two conditions.  
- `filter_grn_diff()` removes low-confidence differential edges.  
  → [https://github.com/oncologylab/episcope/wiki/compare_grn](https://github.com/oncologylab/episcope/wiki/compare_grn)

### Clustering, Topic Modeling, and Hub Discovery
- `cluster_grn()` groups GRN edges by shared activity patterns.  
- `topic_model_grn()` performs topic modeling across GRN edge matrices.  
- `find_hub_tfs()` ranks TFs by centrality metrics within each topic.  
  → [https://github.com/oncologylab/episcope/wiki/topic_model_grn](https://github.com/oncologylab/episcope/wiki/topic_model_grn)

### Visualization
- `plot_topic_network()` visualizes topic-level GRNs (pairwise or Δ).  
- `plot_tf_network()` visualizes TF-centric subnetworks for specified contrasts.  
  → [https://github.com/oncologylab/episcope/wiki/plot_tf_network](https://github.com/oncologylab/episcope/wiki/plot_tf_network)

---

## Pipelines

episcope supports multiple orchestrated pipelines:

- **Snakemake** — scalable, HPC-optimized  
- **Nextflow** — reproducible, cloud-ready  
- **R scripts** — local exploratory runs

👉 [https://github.com/oncologylab/episcope_snakemake](https://github.com/oncologylab/episcope_snakemake)  
👉 [https://github.com/oncologylab/episcope_nextflow](https://github.com/oncologylab/episcope_nextflow)

---

## Help and Documentation

- 📘 [https://github.com/oncologylab/episcope/wiki](https://github.com/oncologylab/episcope/wiki)  
- ❓ [https://github.com/oncologylab/episcope/wiki/FAQ](https://github.com/oncologylab/episcope/wiki/FAQ)  
- 🐛 [https://github.com/oncologylab/episcope/issues](https://github.com/oncologylab/episcope/issues)

---

## How to Cite

Li, Y. *et al.* (in preparation).  
**episcope:** An integrative multi-omics framework for condition-specific gene regulatory network analysis.

---

## License

This project is licensed under the [MIT License](https://github.com/oncologylab/episcope/blob/main/LICENSE).

---

