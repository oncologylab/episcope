episcope — Integrative ATAC–RNA Framework for Condition-Specific Gene Regulatory Networks
===========================================================================================

[![Version](https://img.shields.io/badge/version-1.0.0-blue.svg?style=plastic)](https://github.com/oncologylab/episcope)
[![License](https://img.shields.io/badge/license-MIT-green.svg?style=plastic)](https://github.com/oncologylab/episcope/blob/main/LICENSE)
[![Bioconductor](https://img.shields.io/badge/install%20via-BiocManager-orange.svg?style=plastic)](https://bioconductor.org)
[![Maintenance](https://img.shields.io/badge/Maintained%3F-yes-green.svg?style=plastic)](https://github.com/oncologylab/episcope/graphs/commit-activity)
[![publication](https://img.shields.io/badge/Publication-in%20prep-blue.svg?style=plastic)]()

---

## Introduction

**episcope** is a modular, reproducible framework for building and analyzing **condition-specific gene regulatory networks (GRNs)** from **chromatin accessibility (ATAC)** and **transcriptomic (RNA)** data with footprinting support.

It extends classical ATAC-seq footprinting to **integrative GRN modeling**, enabling users to:

- Load ATAC peak sets, **footprint scores** (from mainstream pipelines), and matched RNA-seq quantifications  
- Align motif-redundant footprints to consensus sites and **bias-correct** footprint intensity  
- Correlate TF binding evidence with TF expression to refine footprint→TF assignments  
- Build **regulatory priors** (GeneHancer, TSS ±30 kb, Hi-C enhancer–promoter)  
- Infer data-specific **TF→enhancer→gene** networks; perform condition lighting and differential GRN analysis  
- Visualize topic- and TF-centric subnetworks (pairwise or delta views)-
- *Perturbation datasets (KO/KD, drug, CRISPRi/a) can be imported and used for **validation** of inferred GRNs*

---

## Installation

episcope is an R package available on GitHub.

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

episcope modules are interoperable and can be executed independently or as part of a pipeline.

```r
library(episcope)

# Example: build a GRN from matched ATAC–RNA dataset
fp  <- load_footprints(footprint_bw = "footprints_corrected.bw", peaks_bed = "peaks.bed")
rna <- load_rna("rna_expression.csv")

fp_aligned   <- align_footprints(fp)
fp_corrected <- correct_footprints(fp_aligned)
tf_map       <- map_tf_to_footprints(fp_corrected, rna)

priors   <- build_regulation_priors(method = "genehancer")
refined1 <- correlate_atac_to_gene(priors, fp, rna)
refined2 <- correlate_fp_to_gene(refined1)

grn <- build_basal_grn(refined2)
validate_grn_perturbation(grn, perturb_db = "perturbdb.sqlite")

# Condition lighting and differential analysis
lit_grn <- light_condition_grn(grn, condition = "GlcLow", replicate_policy = "pooled")
dgrn    <- compare_grn(grn, condition_a = "Ctrl", condition_b = "Stress")
dgrn_f  <- filter_grn_diff(dgrn, min_abs_delta = 0.2)
```

---

## Modules Overview

episcope provides modular functions for each step, from footprint data to condition-specific regulatory networks.

### Data Loading and Preprocessing
- `load_footprints()` imports footprint score tracks and ATAC peaks.  
  - **Sources supported:** [**fptools**](https://github.com/oncologylab/fptools) and [**TOBIAS**](https://github.com/loosolab/TOBIAS) outputs.  
- `load_rna()` loads RNA-seq quantifications from mainstream tools (STAR+featureCounts counts, Salmon/Kallisto TPM/counts).  
- `align_footprints()` consolidates motif-redundant footprint calls into consensus sites. *(Per-motif consolidation; **not** cross-sample alignment.)*  
- `correct_footprints()` applies bias correction and depth normalization.  
- `quantile_normalize_fp_unique()` performs **peak-wise quantile normalization** of footprint scores (see code below).  
  → https://github.com/oncologylab/episcope/wiki/Data-Loading-and-Preprocessing

### TF Assignment and Regulation Priors
- `map_tf_to_footprints()` correlates footprint signals with TF expression to assign regulators.  
- `build_regulation_priors()` defines TF–gene edges from **GeneHancer**, **TSS ±30 kb**, or **Hi-C** contact maps.  
  → https://github.com/oncologylab/episcope/wiki/TF-Assignment-and-Regulation-Priors

### Correlation-Based Refinement
- `correlate_atac_to_gene()` refines enhancer–gene links by ATAC–RNA correlation.  
- `correlate_fp_to_gene()` quantifies TF regulatory activity by correlating footprint scores with gene expression.  
  → https://github.com/oncologylab/episcope/wiki/Correlation‑Based-Refinement-of-Regulatory-Links

### Network Assembly, Validation, and Condition Lighting
- `build_basal_grn()` constructs a dataset-specific GRN.  
- `validate_grn_perturbation()` validates TF–gene links via KO/KD public datasets (SQLite + external perturbation DBs).  
- `light_condition_grn()` identifies **active regulatory edges** within each condition (replicate pooling or strict replicate consensus).  
  → https://github.com/oncologylab/episcope/wiki/Network-Assembly-Validation-and-Condition-Lighting

### Differential GRNs, Clustering, Topics, and Hubs
- `compare_grn()` computes differential networks between two conditions.  
- `filter_grn_diff()` filters differential GRNs by score/correlation/evidence.  
- `cluster_grn()` groups edges by shared activity patterns (Louvain / hierarchical).  
- `topic_model_grn()` performs **LDA topic modeling** across GRN edge matrices.  
- `find_hub_tfs()` ranks TFs by centrality metrics (HITS/degree) within each topic.  
  → https://github.com/oncologylab/episcope/wiki/topic_model_grn

### Interactive Gene Regulatory Network Visualization
- `plot_topic_network_pairwise()` — topic-level pairwise visualization (left/right).  
- `plot_topic_network_delta()` — topic-level **delta** visualization (single-panel).  
- `plot_tf_network_pairwise()` — TF-centric pairwise visualization (left/right).  
- `plot_tf_network_delta()` — TF-centric **delta** visualization (single-panel).  
  → https://github.com/oncologylab/episcope/wiki/plot_tf_network

---

## Pipelines

episcope supports multiple orchestrated pipelines:

- **Snakemake** — scalable, HPC-optimized  
- **Nextflow** — reproducible, cloud-ready  
- **R scripts** — local exploratory runs

👉 https://github.com/oncologylab/episcope_snakemake  
👉 https://github.com/oncologylab/episcope_nextflow

---

## Help and Documentation

- 📘 Docs https://github.com/oncologylab/episcope/wiki  
- ❓ FAQ https://github.com/oncologylab/episcope/wiki/FAQ  
- 🐛 Bug report https://github.com/oncologylab/episcope/issues

---

## How to Cite

Li, Y., Yi, C. *et al.* (in preparation).  
**episcope:** Integrative ATAC–RNA framework for condition-specific gene regulatory network analysis.

---

## License

This project is licensed under the [GNU General Public License v3.0](https://github.com/oncologylab/episcope/blob/main/LICENSE.md).
