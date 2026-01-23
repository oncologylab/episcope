Episcope - Integrative ATAC-RNA Framework for Condition-Specific Gene Regulatory Networks
===========================================================================================

[![Version](https://img.shields.io/badge/version-1.0.0-blue.svg?style=plastic)](https://github.com/oncologylab/episcope)
[![License](https://img.shields.io/badge/license-MIT-green.svg?style=plastic)](https://github.com/oncologylab/episcope/blob/main/LICENSE)
[![Bioconductor](https://img.shields.io/badge/install%20via-BiocManager-orange.svg?style=plastic)](https://bioconductor.org)
[![Maintenance](https://img.shields.io/badge/Maintained%3F-yes-green.svg?style=plastic)](https://github.com/oncologylab/episcope/graphs/commit-activity)
[![publication](https://img.shields.io/badge/Publication-in%20prep-blue.svg?style=plastic)]()

---

## Introduction

**Episcope** is a modular, reproducible framework for building and analyzing **condition-specific gene regulatory networks (GRNs)** from **chromatin accessibility (ATAC)** and **transcriptomic (RNA)** data with footprinting support.

It extends classical ATAC-seq footprinting to **integrative GRN modeling**, enabling users to:

- Load ATAC peak sets, **footprint scores** (from mainstream pipelines), and matched RNA-seq quantifications  
- Align motif-redundant footprints to consensus sites and **bias-correct** footprint intensity  
- Correlate TF binding evidence with TF expression to refine footprint→TF assignments  
- Build **regulatory priors** (GeneHancer, TSS ±30 kb, Hi-C enhancer-promoter)  
- Infer data-specific **TF→enhancer→gene** networks; perform condition lighting and differential GRN analysis  
- Visualize topic- and TF-centric subnetworks (pairwise or delta views)-
- *Perturbation datasets (KO/KD, drug, CRISPRi/a) can be imported and used for **validation** of inferred GRNs*

---

## Installation

Episcope is an R package available on GitHub.

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

Episcope modules are interoperable and can be executed independently or as part
of a pipeline.

```r
library(episcope)
load_episcope_config("episcope_grn.yaml")

# Step 0: load footprints, ATAC, RNA; build strict grn_set
fp_manifest <- load_footprints(root_dir = fp_root_dir, db_name = db, out_dir = fp_out_dir)
fp_aligned <- align_footprints(fp_manifest, output_mode = "distinct")

# Step 1: predict TF binding sites with FP–TF correlations
grn_set <- build_grn_set(...)
grn_set <- grn_add_fp_tf_corr(grn_set, method = "pearson", cores = 20L)
grn_set <- grn_filter_fp_tf_corr(grn_set, r_thr = threshold_fp_tf_corr_r, p_thr = threshold_fp_tf_corr_p)

# Step 2: link TFBS to target genes and light per condition
gh_std <- load_genehancer_panc(file.path("inst", "extdata", "GeneHancer_v5.24_elite_panc.csv"))
fp_res_full_pearson <- correlate_fp_to_genes(grn_set, gh_tbl = gh_std, method = "pearson")
fp_res_full_spearman <- correlate_fp_to_genes(grn_set, gh_tbl = gh_std, method = "spearman")
fp_links_filtered <- filter_links_by_fp_rna_criteria(...)
status_res <- build_link_status_matrix(...)
light_by_condition(...)

# Step 3: differential GRN and topic analysis
run_links_deltas_driver(...)
filter_links_deltas_bulk(...)
run_vae_ctf_multivi(...)
```

---

## Modules Overview

Episcope provides modular functions for each step, from footprint data to condition-specific regulatory networks.

### [**Predict TF binding sites**](https://github.com/oncologylab/episcope/wiki/Predict-TF-binding-sites)
Build condition-aware FP-bound and gene-expression matrices, apply FP-wise QN,
and compute FP–TF correlations (Pearson/Spearman/Kendall; canonical/all modes).
Outputs TF binding probability overviews and optional TFBS BEDs.

### [**Connect TF-occupied enhancers to target genes**](https://github.com/oncologylab/episcope/wiki/Connect-TFs-to-Target-Genes)
Link TFBS to target genes via GeneHancer (or distance/loops), build a
condition-specific link-status matrix, compute TF–gene correlations, and
assemble per-condition GRNs with link_score = r(TF–gene) × FP score.

### [**Build basal GRN & identify active regulatory edges per condition**](https://github.com/oncologylab/episcope/wiki/Build-basal-GRN-&-identify-active-regulatory-edges-per-condition)
Assemble basal TF–peak–gene links and light edges per condition under FP-bound
and gene-expression gating. Emits per-condition link tables and an index.

### [**Perform differential GRN analysis & identify master TFs**](https://github.com/oncologylab/episcope/wiki/Perform-differential-GRN-analysis-&-identify-master-TFs)
Compute per-condition contrasts, filter delta links by expression and delta
thresholds, and run topic modeling with VAE-based workflows.

### [**Generate interactive Topic & TF regulatory hub subnetworks**](https://github.com/oncologylab/episcope/wiki/Interactive-Gene-Regulatory-Network-Visualization)
Generate doc-topic heatmaps, topic delta subnet plots, and pathway summaries
from `topic_models/` outputs.



<!-- ### [**Predict TF binding sites**](https://github.com/oncologylab/episcope/wiki/Predict-TF-binding-sites)
- `load_footprints()` imports footprint score tracks and ATAC peaks.  
  - **Sources supported:** [**fptools**](https://github.com/oncologylab/fptools) and [**TOBIAS**](https://github.com/loosolab/TOBIAS) outputs.  
- `fp_manifest_trim()` 
- `fp_manifest_trim_annots()` 
- `align_footprints()` consolidates motif-redundant footprint calls into consensus sites. *(Per-motif consolidation; **not** cross-sample alignment.)*  
- `qn_footprints()` performs **peak-wise quantile normalization** of footprint scores.  
- `save_footprints()` writes the aligned/normalized footprint data for reuse.
- `filter_footprints()` RNA expression and ATAC peak based footprint filter.
- `process_motifs_tf_corr_in_parallel()` Correlate footprint with predicted binded TFs RNA expression per motif.
- `

- `load_rna()` loads RNA-seq quantifications from mainstream tools (STAR+featureCounts counts, Salmon/Kallisto TPM/counts).  

### [**Connect TF-occupied enhancers to target genes**](https://github.com/oncologylab/episcope/wiki/TF-Assignment-and-Regulation-Priors)
- `map_tf_to_footprints()` correlates footprint signals with TF expression to assign regulators.  
- `build_regulation_priors()` defines TF-gene edges from **GeneHancer**, **TSS ±30 kb**, or **Hi-C** contact maps.  

### [**Connect TF-occupied enhancers to target genes**](https://github.com/oncologylab/episcope/wiki/Correlation‑Based-Refinement-of-Regulatory-Links)
- `correlate_atac_to_gene()` refines enhancer-gene links by ATAC-RNA correlation.  
- `correlate_fp_to_gene()` quantifies TF regulatory activity by correlating footprint scores with gene expression.  

### [**Build basal GRN & identify active regulatory edges per condition**](https://github.com/oncologylab/episcope/wiki/Network-Assembly-Validation-and-Condition-Lighting)
- `build_basal_grn()` constructs a dataset-specific GRN.  
- `validate_grn_perturbation()` validates TF-gene links via KO/KD public datasets (SQLite + external perturbation DBs).  
- `light_condition_grn()` identifies **active regulatory edges** within each condition (replicate pooling or strict replicate consensus).  

### [**Perform differential GRN analysis & identify master TFs**](https://github.com/oncologylab/episcope/wiki/Differential-GRNs-Clustering-Topics-and-Hubs)
- `compare_grn()` computes differential networks between two conditions.  
- `filter_grn_diff()` filters differential GRNs by score/correlation/evidence.  
- `cluster_grn()` groups edges by shared activity patterns (Louvain / hierarchical).  
- `topic_model_grn()` performs **LDA topic modeling** across GRN edge matrices.  
- `find_hub_tfs()` ranks TFs by centrality metrics (HITS/degree) within each topic.  

### [**Generate interactive Topic & TF regulatory hub subnetworks**](https://github.com/oncologylab/episcope/wiki/Interactive-Gene-Regulatory-Network-Visualization)
- `plot_topic_network_pairwise()` - topic-level pairwise visualization (left/right).  
- `plot_topic_network_delta()` - topic-level **delta** visualization (single-panel).  
- `plot_tf_network_pairwise()` - TF-centric pairwise visualization (left/right).  
- `plot_tf_network_delta()` - TF-centric **delta** visualization (single-panel).   -->

---

## Pipelines

episcope supports multiple orchestrated pipelines:

- **Snakemake** - scalable, HPC-optimized  
- **Nextflow** - reproducible, cloud-ready  
- **R scripts** - local exploratory runs

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
**Episcope:** Integrative ATAC-RNA framework for condition-specific gene regulatory network analysis.

---

## License

This project is licensed under the [GNU General Public License v3.0](https://github.com/oncologylab/episcope/blob/main/LICENSE.md).
