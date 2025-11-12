summary_files <- list.files(path = lighting_folder, pattern = "_lda_K20_summary.csv", full.names = TRUE)
summary_file <- summary_files[1]
readr::read_csv(summary_file)
