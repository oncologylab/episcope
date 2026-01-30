options(future.globals.maxSize = 64 * 1024^3)
fp_res_full_kendall <- correlate_fp_to_genes(
  grn_set             = grn_set,
  atac_gene_corr_kept = atac_res$atac_gene_corr_full,
  fdr                 = threshold_fp_gene_corr_p,
  r_abs_min           = threshold_fp_gene_corr_abs_r,
  method              = "kendall",   # c("pearson", "spearman", "kendall")
  workers             = 20,
  cache_dir           = file.path(base_dir, "cache", "fp_gene_corr"),
  cache_tag           = "nutrient_genehancer_kendall",
  cache_chunk_size    = 5000L,
  cache_verbose       = TRUE
)
