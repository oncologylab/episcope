# CraftGRN <img src="man/figures/logo.png" align="right" height="139" alt="CraftGRN logo" />

[![Version](https://img.shields.io/badge/version-0.1.2-2C3E50.svg?style=flat-square)](https://github.com/oncologylab/craftgrn)
[![License](https://img.shields.io/badge/license-GPLv3-16A085.svg?style=flat-square)](https://github.com/oncologylab/craftgrn/blob/main/LICENSE.md)
[![Documentation](https://img.shields.io/badge/docs-pkgdown-4DBBD5.svg?style=flat-square)](https://oncologylab.github.io/craftgrn/)
[![pkgdown](https://img.shields.io/github/actions/workflow/status/oncologylab/craftgrn/pkgdown.yaml?branch=main&label=pkgdown&style=flat-square&color=27AE60)](https://github.com/oncologylab/craftgrn/actions/workflows/pkgdown.yaml)
[![Last commit](https://img.shields.io/github/last-commit/oncologylab/craftgrn.svg?style=flat-square&color=9B59B6)](https://github.com/oncologylab/craftgrn/commits/main)
[![Publication](https://img.shields.io/badge/publication-in%20prep-E67E22.svg?style=flat-square)]()

[![ATAC-seq](https://img.shields.io/badge/ATAC--seq-footprints-00A087.svg?style=flat-square)](https://oncologylab.github.io/craftgrn/)
[![RNA-seq](https://img.shields.io/badge/RNA--seq-expression-4DBBD5.svg?style=flat-square)](https://oncologylab.github.io/craftgrn/)
[![Module 1](https://img.shields.io/badge/module%201-TFBS-9B59B6.svg?style=flat-square)](https://oncologylab.github.io/craftgrn/articles/craftgrn.html)
[![Module 2](https://img.shields.io/badge/module%202-TF--target-27AE60.svg?style=flat-square)](https://oncologylab.github.io/craftgrn/articles/craftgrn.html)
[![Module 3](https://img.shields.io/badge/module%203-regulatory%20topics-E74C3C.svg?style=flat-square)](https://oncologylab.github.io/craftgrn/articles/craftgrn.html)
[![GRN](https://img.shields.io/badge/output-condition--specific%20GRNs-2C3E50.svg?style=flat-square)](https://oncologylab.github.io/craftgrn/)

## Introduction

**CraftGRN** is a modular framework for integrating chromatin accessibility
profiles from ATAC-seq with matched RNA-seq expression data to infer
condition-specific transcription factor binding sites and reconstruct dynamic
gene regulatory networks.

CraftGRN helps users:

- Collapse overlapping TF motif footprints into consensus, site- and
  motif-nonredundant footprint clusters.
- Infer condition-specific canonical and non-canonical TF binding sites by
  correlating TF expression with footprint or chromatin accessibility scores.
- Refine TF->TFBS->gene regulatory priors using enhancer-gene maps, genomic
  proximity, or user-supplied chromatin interaction data.
- Extract active regulatory links within each condition and compare links
  between conditions.
- Learn regulatory topics from RNA and footprint signals using topic modeling
  and VAE-based representations.
- Generate summaries and visualizations for topic- and condition-specific
  regulatory programs.

<img src="https://raw.githubusercontent.com/oncologylab/craftgrn/main/figures/pipeline_full.svg" alt="CraftGRN pipeline" width="800">

## Installation

CraftGRN can be installed from GitHub:

```r
# Using remotes
remotes::install_github("oncologylab/craftgrn")

# or using pak
pak::pak("oncologylab/craftgrn")
```

Common CRAN and Bioconductor dependencies can be installed with:

```r
install.packages(c("igraph", "ggplot2", "data.table", "BiocManager"))
BiocManager::install(c("DESeq2", "GenomicRanges", "SummarizedExperiment"))
```

## Pipeline Overview

CraftGRN is organized as a three-module workflow.

### Module 1: Predict TF Binding Sites

Module 1 loads matched ATAC, RNA, metadata, and optional footprint score files,
then prepares a multiomic data object for downstream regulatory analysis.

Primary package functions:

- `load_prep_multiomic_data()` loads, filters, aligns, and prepares multiomic
  inputs from a YAML configuration file.
- `correlate_tf_to_fp()` correlates TF expression with footprint or peak signal
  across matched conditions to nominate condition-specific TF binding sites.

<img src="https://raw.githubusercontent.com/oncologylab/craftgrn/main/figures/module_1.svg" alt="Module 1 workflow" width="800">

### Module 2: Connect TFs to Target Genes

Module 2 links TF binding sites to candidate target genes using enhancer-gene
maps, genomic distance windows, or 3D chromatin interaction priors. Candidate
TF->TFBS->target links are filtered by condition-specific expression, binding,
footprint or peak signal, and cross-condition correlation evidence.

<img src="https://raw.githubusercontent.com/oncologylab/craftgrn/main/figures/module_2.svg" alt="Module 2 workflow" width="800">

### Module 3: Learn Regulatory Topics and Visualize Differential GRNs

Module 3 compares condition-specific regulatory links, builds joint RNA and
footprint document-term matrices, trains topic models, assigns regulatory links
to topics, and summarizes pathway and master TF programs.

Primary package functions:

- `train_topic_models()` trains regulatory topic models across a user-defined
  topic-number grid.
- `extract_regulatory_topics()` assigns links and terms to selected regulatory
  topics.
- `run_diff_grn_pathway_analysis()` and related plotting functions summarize
  pathway programs.
- `run_diff_grn_master_tf_summary()` and related plotting functions summarize
  master TF connectivity.

<img src="https://raw.githubusercontent.com/oncologylab/craftgrn/main/figures/module_3.svg" alt="Module 3 workflow" width="800">

## Get Started

For a module-by-module tutorial, see the
[Get started article](https://oncologylab.github.io/craftgrn/articles/craftgrn.html).

## Documentation

- Website: https://oncologylab.github.io/craftgrn/
- Reference: https://oncologylab.github.io/craftgrn/reference/
- Issues: https://github.com/oncologylab/craftgrn/issues

## Citation

Li, Y., Yi, C. et al. (in preparation).
**CraftGRN:** Integrative ATAC-RNA framework for condition-specific gene
regulatory network analysis.

## License

This project is licensed under the GNU General Public License v3.0.
