module1_tiny_fixture <- function() {
  conditions <- c("cond_a", "cond_b", "cond_c", "cond_d")

  fp_score_condition_qn <- tibble::tibble(
    peak_ID = c(
      "chr1:100-140",
      "chr1:200-240",
      "chr1:300-340",
      "chr2:100-150",
      "chr2:220-260",
      "chr3:500-550"
    ),
    cond_a = c(1, 1, 4, 8, 2, 5),
    cond_b = c(2, 1, 3, 7, 2, 5),
    cond_c = c(8, 9, 2, 1, 2, 5),
    cond_d = c(9, 8, 1, 1, 2, 5)
  )

  fp_bound_condition <- tibble::tibble(
    peak_ID = fp_score_condition_qn$peak_ID,
    cond_a = c(1L, 1L, 1L, 1L, 0L, 1L),
    cond_b = c(1L, 1L, 1L, 1L, 0L, 1L),
    cond_c = c(1L, 1L, 1L, 0L, 0L, 1L),
    cond_d = c(1L, 1L, 1L, 0L, 0L, 1L)
  )

  fp_annotation <- tibble::tibble(
    fp_peak = fp_score_condition_qn$peak_ID,
    atac_peak = c(
      "chr1:90-160",
      "chr1:180-260",
      "chr1:280-360",
      "chr2:80-170",
      "chr2:200-280",
      "chr3:480-570"
    ),
    motifs = c("MOTIF_TF_A", "MOTIF_TF_B", "MOTIF_TF_C", "MOTIF_TF_D", "MOTIF_TF_A", "MOTIF_TF_E"),
    tfs = c("TF_A", "TF_B", "TF_C", "TF_D", "TF_A", "TF_E")
  )

  rna_condition <- tibble::tibble(
    ensembl_gene_id = paste0("gene", seq_len(5)),
    HGNC = c("TF_A", "TF_B", "TF_C", "TF_D", "TF_NO_PASS"),
    cond_a = c(1, 9, 4, 8, 2),
    cond_b = c(2, 8, 3, 7, 2),
    cond_c = c(8, 2, 2, 1, 2),
    cond_d = c(9, 1, 1, 1, 2)
  )

  rna_expressed <- tibble::tibble(
    ensembl_gene_id = rna_condition$ensembl_gene_id,
    HGNC = rna_condition$HGNC,
    cond_a = c(1L, 1L, 1L, 1L, 1L),
    cond_b = c(1L, 1L, 1L, 1L, 1L),
    cond_c = c(1L, 1L, 1L, 1L, 1L),
    cond_d = c(1L, 1L, 1L, 1L, 1L)
  )

  omics_data <- list(
    fp_score_condition_qn = fp_score_condition_qn,
    fp_bound_condition = fp_bound_condition,
    fp_annotation = fp_annotation,
    rna_condition = rna_condition,
    rna_expressed = rna_expressed,
    tf_list = c("TF_A", "TF_B", "TF_C", "TF_D", "TF_NO_PASS", "TF_ABSENT")
  )

  list(
    omics_data = omics_data,
    expected_multi_tf_fp = "chr1:100-140",
    expected_excluded_fp = "chr2:220-260",
    expressed_tfs = c("TF_A", "TF_B", "TF_C", "TF_D", "TF_NO_PASS")
  )
}
