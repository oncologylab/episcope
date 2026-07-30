test_that("overall pathway CSV is not limited by plot pathway caps", {
  out_dir <- file.path(
    tempdir(),
    paste0("topic-pathway-full-table-", sample.int(1e8, 1L))
  )
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  out_file <- file.path(out_dir, "topic_pathway_enrichment.pdf")
  topic_terms <- data.table::data.table(
    topic = c("Topic1", "Topic1", "Topic1"),
    topic_num = c(1L, 1L, 1L),
    term_id = c("GENE:A", "GENE:B", "GENE:C"),
    score = c(0.8, 0.7, 0.6),
    in_topic = TRUE
  )

  with_mocked_bindings(
    plot_topic_pathway_enrichment_heatmap(
      topic_terms = topic_terms,
      edges_docs = NULL,
      option_label = "opt3_gene_fc_expr",
      out_file = out_file,
      dbs = "Toy_DB",
      pathway_species = "human",
      min_genes = 2L,
      top_n_per_topic = Inf,
      dot_top_n_per_topic = 2L,
      max_pathways = 1L,
      make_heatmap = FALSE,
      make_dotplot = TRUE,
      pathway_backend = "enrichly"
    ),
    .pathway_backend_available = function(pathway_backend) TRUE,
    .run_enrichr_cached = function(...) {
      list(Toy_DB = data.frame(
        Term = c("Term one", "Term two"),
        Adjusted.P.value = c(0.01, 0.02),
        P.value = c(0.001, 0.002),
        Overlap = c("2/10", "2/12"),
        Genes = c("A;B", "B;C"),
        Combined.Score = c(8, 6),
        Odds.Ratio = c(2, 1.5),
        query_size = c(3L, 3L),
        term_size = c(10L, 12L),
        background_size = c(20000L, 20000L),
        check.names = FALSE
      ))
    }
  )

  full <- data.table::fread(file.path(
    out_dir,
    "topic_pathway_enrichment_topic_terms.csv"
  ))
  plotted <- data.table::fread(file.path(
    out_dir,
    "topic_pathway_enrichment_dotplot.csv"
  ))
  expect_equal(sort(full$pathway), c("Toy_DB: Term one", "Toy_DB: Term two"))
  expect_equal(data.table::uniqueN(plotted$pathway), 1L)
})
