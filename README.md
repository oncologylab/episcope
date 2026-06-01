# CraftGRN <img src="man/figures/logo.png" align="right" height="139" alt="CraftGRN logo" />

[![Version](https://img.shields.io/badge/version-0.1.2-2C3E50.svg?style=flat-square)](https://github.com/oncologylab/craftgrn)
[![License](https://img.shields.io/badge/license-GPLv3-16A085.svg?style=flat-square)](https://github.com/oncologylab/craftgrn/blob/main/LICENSE.md)
[![Documentation](https://img.shields.io/badge/docs-pkgdown-4DBBD5.svg?style=flat-square)](https://oncologylab.github.io/craftgrn/)
[![pkgdown](https://img.shields.io/github/actions/workflow/status/oncologylab/craftgrn/pkgdown.yaml?branch=main&label=pkgdown&style=flat-square&color=27AE60)](https://github.com/oncologylab/craftgrn/actions/workflows/pkgdown.yaml)
[![Last commit](https://img.shields.io/github/last-commit/oncologylab/craftgrn.svg?style=flat-square&color=9B59B6)](https://github.com/oncologylab/craftgrn/commits/main)
[![Publication](https://img.shields.io/badge/publication-in%20prep-E67E22.svg?style=flat-square)]()

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

## Demo Data

CraftGRN keeps demo datasets outside the source package so installation remains
small and CRAN-friendly. The package helper reports any configured external
demo bundles:

```r
craftgrn::craftgrn_demo_data_info()
```

No external demo bundle is currently configured. To run your own project, point
CraftGRN at a project-level YAML file:

```r
config <- "project.yaml"

omics <- craftgrn::load_prep_multiomic_data(
  config = config,
  label_col = "strict_match_rna",
  do_preprocess = FALSE,
  verbose = TRUE
)

module1 <- craftgrn::predict_tfbs(
  omics_data = omics,
  out_dir = "predict_tf_binding_sites",
  output_format = "auto",
  write_stats = FALSE,
  verbose = TRUE
)
```

Troubleshooting:

- If `craftgrn_demo_data_info()` returns zero rows, no public demo bundle is
  currently advertised by this package version.
- If paths fail after moving a project folder, keep `project.yaml` in the
  project directory and pass that config path explicitly. A portable project
  config should use `base_dir: "."`.
- If memory is limited, start with `load_prep_multiomic_data()` and Module 1
  before running Module 2.

## Pipeline Overview

CraftGRN is organized as a three-module workflow.

### Module 1: Predict TF Binding Sites

Module 1 loads matched ATAC, RNA, metadata, and optional footprint score files,
then prepares a multiomic data object for downstream regulatory analysis.

Primary package functions:

- `load_prep_multiomic_data()` loads, filters, aligns, and prepares multiomic
  inputs from a YAML configuration file.
- `predict_tfbs()` performs direct-bound footprint filtering and TF binding site
  prediction across matched conditions.
- `build_module1_qc_report()` writes an HTML QC report for run parameters, input gates, canonical support, correlation diagnostics, predicted TFBS chunk integrity, top TFs/FPs, condition support, warning checks, and related Module 1 artifacts. The report uses multiple static plot types, including processing funnels, density curves, scatter summaries, heatmaps, lollipop rank plots, and cumulative curves.

<img src="https://raw.githubusercontent.com/oncologylab/craftgrn/main/figures/module_1.svg" alt="Module 1 workflow" width="800">

### Module 2: Connect TFs to Target Genes

Module 2 links TF binding sites to candidate target genes using enhancer-gene
maps, genomic distance windows, or 3D chromatin interaction priors. Candidate
TF->TFBS->target links are filtered by condition-specific expression, binding,
footprint or peak signal, and cross-condition correlation evidence.

Primary package functions:

- `link_tf_targets()` connects predicted TFBS to target genes using TF-target and FP-target correlations.
- `build_module2_qc_report()` writes an HTML QC report for compact handoff checks, TF-target and FP-target filters, candidate source and distance-to-TSS evidence, final-link integrity, condition activity, warning checks, top TF/target/FP summaries, and related browser reports. The report combines relational flow diagrams, density and cumulative distance plots, scatter summaries, heatmaps, and lollipop rank plots.

<img src="https://raw.githubusercontent.com/oncologylab/craftgrn/main/figures/module_2.svg" alt="Module 2 workflow" width="800">

### Module 3: Learn Regulatory Topics and Visualize Differential GRNs

Module 3 compares condition-specific regulatory links, builds joint RNA and
footprint document-term matrices, trains topic models, assigns regulatory links
to topics, and summarizes pathway and master TF programs.

Primary package functions:

- `train_topic_models()` trains regulatory topic models across a user-defined
  topic-number grid using the native `warp_omp` WarpLDA sampler by default.
  Use `warplda_sampler = "warp_ref"` only when you need a slower sequential
  fixed-seed reference run from the native backend.
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
