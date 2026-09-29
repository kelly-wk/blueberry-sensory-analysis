if (!requireNamespace("jsonlite", quietly = TRUE)) stop("Package jsonlite is required")

required <- c(
  "REPORT.md", "DATA_CARD.md", "data/blueberry_2013_sensory.csv",
  "data/source_manifest.json", "results/summary.json",
  "results/pca_loading_bootstrap_intervals.csv",
  "results/cluster_internal_validation.csv",
  "results/cluster_subsample_stability.csv",
  "results/genotype_cv_predictions.csv",
  "results/artifact_manifest.csv",
  "figures/pca_biplot.png", "figures/cluster_diagnostics.png",
  "figures/genotype_cv_confusion.png"
)
stopifnot(all(file.exists(required)), all(file.info(required)$size > 0))

data <- utils::read.csv("data/blueberry_2013_sensory.csv", check.names = FALSE)
stopifnot(nrow(data) == 50L, length(unique(data$sample_id)) == 50L,
          length(unique(data$environment)) == 9L,
          length(unique(data$genotype)) == 6L,
          all(complete.cases(data)))
stopifnot(identical(as.integer(table(data$genotype)[c("Emerald", "Endura", "Farthing",
                                                     "Meadowlark", "Primadonna", "Scintilla")]),
                    c(9L, 9L, 9L, 9L, 8L, 6L)))

summary <- jsonlite::read_json("results/summary.json", simplifyVector = TRUE)
stopifnot(summary$metadata$observations == 50L,
          summary$metadata$environments == 9L,
          summary$metadata$genotypes == 6L,
          summary$pca$block_bootstrap_replicates == 1000L,
          summary$clustering$gap_bootstrap_replicates == 200L,
          summary$clustering$gap_first_se_selected_k_standardized >= 1L,
          summary$clustering$gap_first_se_selected_k_standardized <= 8L,
          summary$clustering$exploratory_forced_k_standardized >= 2L,
          summary$clustering$exploratory_forced_k_standardized <= 8L,
          summary$leave_one_environment_out_genotype_validation$accuracy >= 0,
          summary$leave_one_environment_out_genotype_validation$accuracy <= 1)

boot <- utils::read.csv("results/pca_bootstrap_stability.csv")
cluster_stability <- utils::read.csv("results/cluster_subsample_stability.csv")
predictions <- utils::read.csv("results/genotype_cv_predictions.csv")
manifest <- utils::read.csv("results/artifact_manifest.csv")
stopifnot(nrow(boot) == 1000L,
          nrow(cluster_stability) == 500L,
          nrow(predictions) == 50L,
          all(predictions$true_genotype %in% predictions$predicted_genotype),
          all(nchar(manifest$sha256) == 64L),
          !any(grepl("^/|file://|@", manifest$path, ignore.case = TRUE)))

report <- paste(readLines("REPORT.md", warn = FALSE), collapse = "\n")
stopifnot(grepl("non-causal|no causal|makes no causal claim", report, ignore.case = TRUE),
          grepl("gap statistic", report, ignore.case = TRUE),
          grepl("CC BY 4.0", report, fixed = TRUE))

cat("test_outputs.R: all tests passed\n")
