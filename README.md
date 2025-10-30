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
- Align and correct footprint signals across conditions  
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
epi_links <- episcope::LoadFootprints("footprints_corrected.bw", "peaks.bed")
epi_rna   <- episcope::LoadRNA("rna_expression.csv")

aligned_fp <- episcope::AlignFootprints(epi_links)
tf_corr    <- episcope::CorrelateFootprintsToTFs(aligned_fp, epi_rna)

priors     <- episcope::BuildRegulationPriors(method = "genehancer")
refined1   <- episcope::RefineByATACCorrelation(priors, epi_links, epi_rna)
refined2   <- episcope::RefineByFootprintGeneCorrelation(refined1)

grn        <- episcope::BuildBasalNetwork(refined2)
episcope::ValidateNetwork(grn, perturb_db = "perturbdb.sqlite")

episcope::RenderConditionNetwork(grn, condition = "GlcLow")
episcope::CompareDifferentialNetworks(grn, conditionA="Ctrl", conditionB="Stress")
```

---

## Modules Overview

Each module has its own dedicated **Wiki** page with examples and detailed documentation.

---

### 🧩 **1. Data Loading**
| Module | Description | Wiki |
|--------|--------------|------|
| **LoadFootprints** | Load ATAC-seq peaks and footprint score data from fp-tools or TOBIAS outputs. | [Wiki → LoadFootprints](https://github.com/oncologylab/episcope/wiki/LoadFootprints) |
| **LoadRNA** | Load and normalize matched RNA-seq data (TPM/FPKM/Counts). | [Wiki → LoadRNA](https://github.com/oncologylab/episcope/wiki/LoadRNA) |

---

### ⚙️ **2. Footprint Alignment & Correction**
| Module | Description | Wiki |
|--------|--------------|------|
| **AlignFootprints** | Perform data-driven footprint alignment across samples. | [Wiki → AlignFootprints](https://github.com/oncologylab/episcope/wiki/AlignFootprints) |
| **CorrectFootprints** | Apply signal-based correction to harmonize footprint bias. | [Wiki → CorrectFootprints](https://github.com/oncologylab/episcope/wiki/CorrectFootprints) |

---

### 🔗 **3. TF Assignment via Correlation**
| Module | Description | Wiki |
|--------|--------------|------|
| **CorrelateFootprintsToTFs** | Correlate footprint scores to TF expression to assign footprints to regulators. | [Wiki → CorrelateFootprintsToTFs](https://github.com/oncologylab/episcope/wiki/CorrelateFootprintsToTFs) |

---

### 🧬 **4. Build Regulation Priors**
| Module | Description | Wiki |
|--------|--------------|------|
| **BuildRegulationPriors** | Generate TF–gene priors using GeneHancer, ±30kb TSS, or Hi-C data. | [Wiki → BuildRegulationPriors](https://github.com/oncologylab/episcope/wiki/BuildRegulationPriors) |

---

### 🔄 **5. Regulation Refinement**
| Module | Description | Wiki |
|--------|--------------|------|
| **RefineByATACCorrelation** | Correlate ATAC peak signals with gene expression to refine priors. | [Wiki → RefineByATACCorrelation](https://github.com/oncologylab/episcope/wiki/RefineByATACCorrelation) |
| **RefineByFootprintGeneCorrelation** | Correlate footprint signals with target gene expression to identify active edges. | [Wiki → RefineByFootprintGeneCorrelation](https://github.com/oncologylab/episcope/wiki/RefineByFootprintGeneCorrelation) |

---

### 🧠 **6. Network Construction**
| Module | Description | Wiki |
|--------|--------------|------|
| **BuildBasalNetwork** | Build the condition-independent (basal) GRN from refined data. | [Wiki → BuildBasalNetwork](https://github.com/oncologylab/episcope/wiki/BuildBasalNetwork) |

---

### 🧪 **7. Perturbation Validation**
| Module | Description | Wiki |
|--------|--------------|------|
| **ValidateNetwork** | Validate TF–gene edges using SQLite perturbation DB or PerturbDB resources. | [Wiki → ValidateNetwork](https://github.com/oncologylab/episcope/wiki/ValidateNetwork) |

---

### 🌐 **8. Condition-Specific GRN Lighting**
| Module | Description | Wiki |
|--------|--------------|------|
| **RenderConditionNetwork** | Visualize GRN activity per condition with replicate pooling or strict integration. | [Wiki → RenderConditionNetwork](https://github.com/oncologylab/episcope/wiki/RenderConditionNetwork) |

---

### ⚖️ **9. Differential GRN Analysis**
| Module | Description | Wiki |
|--------|--------------|------|
| **CompareDifferentialNetworks** | Identify edges/nodes differing between conditions. | [Wiki → CompareDifferentialNetworks](https://github.com/oncologylab/episcope/wiki/CompareDifferentialNetworks) |
| **FilterDifferentialNetworks** | Filter differential GRNs based on score or correlation thresholds. | [Wiki → FilterDifferentialNetworks](https://github.com/oncologylab/episcope/wiki/FilterDifferentialNetworks) |

---

### 🧩 **10. Clustering & Topic Modeling**
| Module | Description | Wiki |
|--------|--------------|------|
| **ClusterNetworks** | Cluster edges by regulatory activity using Louvain or hierarchical clustering. | [Wiki → ClusterNetworks](https://github.com/oncologylab/episcope/wiki/ClusterNetworks) |
| **TopicModelNetworks** | Perform LDA-based topic modeling across GRN edge matrices. | [Wiki → TopicModelNetworks](https://github.com/oncologylab/episcope/wiki/TopicModelNetworks) |

---

### 🌟 **11. Hub TF Analysis & Visualization**
| Module | Description | Wiki |
|--------|--------------|------|
| **FindHubTFs** | Identify hub TFs in each topic using HITS or degree centrality. | [Wiki → FindHubTFs](https://github.com/oncologylab/episcope/wiki/FindHubTFs) |
| **RenderTopicNetworks** | Visualize topic-level GRNs (pairwise or delta). | [Wiki → RenderTopicNetworks](https://github.com/oncologylab/episcope/wiki/RenderTopicNetworks) |
| **RenderTFNetworks** | Visualize TF-centric subnetworks (pairwise or delta). | [Wiki → RenderTFNetworks](https://github.com/oncologylab/episcope/wiki/RenderTFNetworks) |

---

## Pipelines

episcope supports multiple orchestrated pipelines:

- **Snakemake** — scalable, HPC-optimized  
- **Nextflow** — reproducible, cloud-ready  
- **R scripts** — local exploratory runs

👉 [episcope_snakemake](https://github.com/oncologylab/episcope_snakemake)  
👉 [episcope_nextflow](https://github.com/oncologylab/episcope_nextflow)

---

## Help and Documentation

- 📘 [Wiki and Tutorials](https://github.com/oncologylab/episcope/wiki)
- ❓ [FAQ](https://github.com/oncologylab/episcope/wiki/FAQ)
- 🐛 [Report Issues](https://github.com/oncologylab/episcope/issues)

---

## How to Cite

Li, Y., Yi, C *et al.* (in preparation).  
**episcope:** An integrative multi-omics framework for condition-specific gene regulatory network analysis.

---

## License

This project is licensed under the [MIT License](https://github.com/oncologylab/episcope/blob/main/LICENSE).
