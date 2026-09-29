options(stringsAsFactors = FALSE, scipen = 999)

if (!requireNamespace("cluster", quietly = TRUE)) stop("Package cluster is required")
if (!requireNamespace("jsonlite", quietly = TRUE)) stop("Package jsonlite is required")
source("R/analysis_helpers.R")

seed <- 20260929L
pca_bootstraps <- 1000L
gap_bootstraps <- 200L
cluster_subsamples <- 500L
cluster_permutations <- 4999L
cv_permutations <- 999L

dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

data <- utils::read.csv("data/blueberry_2013_sensory.csv", check.names = FALSE)
assert_sensory_data(data)
traits <- sensory_variables()
x <- as.matrix(data[traits])
rownames(x) <- data$sample_id
z <- scale(x)
centered <- scale(x, center = TRUE, scale = FALSE)

sample_flow <- data.frame(
  stage = c("published S2 sample rows", "year 2013 and five outcomes complete",
            "unique environments", "unique genotypes"),
  count = c(153L, nrow(data), length(unique(data$environment)),
            length(unique(data$genotype)))
)
utils::write.csv(sample_flow, "results/sample_flow.csv", row.names = FALSE)

descriptive <- do.call(rbind, lapply(traits, function(variable) {
  values <- data[[variable]]
  data.frame(variable = variable, n = length(values), mean = mean(values),
             sd = stats::sd(values), median = stats::median(values),
             minimum = min(values), maximum = max(values))
}))
utils::write.csv(descriptive, "results/descriptive_statistics.csv", row.names = FALSE)

correlations <- stats::cor(x)
correlation_table <- data.frame(variable = rownames(correlations), correlations,
                                check.names = FALSE)
utils::write.csv(correlation_table, "results/sensory_correlations.csv", row.names = FALSE)

unscaled_pca <- fit_pca(x, scale_data = FALSE)
scaled_pca <- fit_pca(x, scale_data = TRUE)
pca_fits <- list(unscaled = unscaled_pca, z_scaled = scaled_pca)

pca_variance <- do.call(rbind, lapply(names(pca_fits), function(mode) {
  fit <- pca_fits[[mode]]
  pve <- fit$sdev^2 / sum(fit$sdev^2)
  data.frame(scaling = mode, component = paste0("PC", seq_along(pve)),
             standard_deviation = fit$sdev,
             explained_variance = pve,
             cumulative_variance = cumsum(pve))
}))
utils::write.csv(pca_variance, "results/pca_explained_variance.csv", row.names = FALSE)

pca_loadings <- do.call(rbind, lapply(names(pca_fits), function(mode) {
  fit <- pca_fits[[mode]]
  correlation <- correlation_loadings(fit)
  do.call(rbind, lapply(seq_len(ncol(fit$rotation)), function(component) {
    data.frame(
      scaling = mode,
      variable = rownames(fit$rotation),
      component = paste0("PC", component),
      eigenvector_loading = fit$rotation[, component],
      correlation_loading = correlation[, component]
    )
  }))
}))
utils::write.csv(pca_loadings, "results/pca_loadings.csv", row.names = FALSE)

raw_to_scaled <- align_two_axes(unscaled_pca$rotation[, 1:2, drop = FALSE],
                                scaled_pca$rotation[, 1:2, drop = FALSE])
raw_scores_aligned <- unscaled_pca$x[, raw_to_scaled$permutation, drop = FALSE]
raw_scores_aligned <- sweep(raw_scores_aligned, 2, raw_to_scaled$signs, `*`)
angles <- principal_angles_degrees(unscaled_pca$rotation[, 1:2, drop = FALSE],
                                   scaled_pca$rotation[, 1:2, drop = FALSE])
scale_sensitivity <- data.frame(
  metric = c("pc1_score_correlation", "pc2_score_correlation",
             "first_principal_angle_degrees", "second_principal_angle_degrees"),
  value = c(stats::cor(raw_scores_aligned[, 1], scaled_pca$x[, 1]),
            stats::cor(raw_scores_aligned[, 2], scaled_pca$x[, 2]), angles)
)
utils::write.csv(scale_sensitivity, "results/pca_scaling_sensitivity.csv",
                 row.names = FALSE)

pca_bootstrap <- block_bootstrap_pca(x, data$environment,
                                     replicates = pca_bootstraps, seed = seed)
pca_loading_intervals <- summarize_loading_bootstrap(pca_bootstrap$loadings)
utils::write.csv(pca_loading_intervals,
                 "results/pca_loading_bootstrap_intervals.csv", row.names = FALSE)
utils::write.csv(pca_bootstrap$stability,
                 "results/pca_bootstrap_stability.csv", row.names = FALSE)

scaled_validation <- cluster_internal_validation(z, seed = seed)
scaled_validation$scaling <- "z_scaled"
unscaled_validation <- cluster_internal_validation(centered, seed = seed + 100L)
unscaled_validation$scaling <- "unscaled"

scaled_gap <- gap_statistic(z, bootstraps = gap_bootstraps, seed = seed)
unscaled_gap <- gap_statistic(centered, bootstraps = gap_bootstraps,
                              seed = seed + 100L)

gap_rows <- rbind(
  transform(scaled_gap$table, scaling = "z_scaled"),
  transform(unscaled_gap$table, scaling = "unscaled")
)
names(gap_rows)[names(gap_rows) == "SE.sim"] <- "gap_standard_error"
names(gap_rows)[names(gap_rows) == "E.logW"] <- "expected_log_within_dispersion"
names(gap_rows)[names(gap_rows) == "logW"] <- "log_within_dispersion"
utils::write.csv(gap_rows, "results/gap_statistic.csv", row.names = FALSE)

internal_validation <- rbind(scaled_validation, unscaled_validation)
utils::write.csv(internal_validation, "results/cluster_internal_validation.csv",
                 row.names = FALSE)

mean_silhouette_by_k <- aggregate(
  mean_silhouette ~ k,
  data = scaled_validation,
  FUN = mean
)
forced_k <- mean_silhouette_by_k$k[which.max(mean_silhouette_by_k$mean_silhouette)]
raw_mean_silhouette_by_k <- aggregate(
  mean_silhouette ~ k,
  data = unscaled_validation,
  FUN = mean
)
raw_forced_k <- raw_mean_silhouette_by_k$k[
  which.max(raw_mean_silhouette_by_k$mean_silhouette)
]

cluster_stability <- block_subsample_cluster_stability(
  x, data$environment, forced_k,
  replicates = cluster_subsamples, seed = seed
)
utils::write.csv(cluster_stability$draws,
                 "results/cluster_subsample_stability.csv", row.names = FALSE)

raw_kmeans <- best_kmeans(centered, forced_k, seed + 301L)$cluster
raw_ward <- stats::cutree(stats::hclust(stats::dist(centered), method = "ward.D2"),
                          forced_k)

assignments <- data.frame(
  sample_id = data$sample_id,
  environment = data$environment,
  genotype = data$genotype,
  scaled_kmeans_cluster = cluster_stability$baseline_kmeans,
  scaled_ward_cluster = cluster_stability$baseline_ward,
  unscaled_kmeans_cluster = raw_kmeans,
  unscaled_ward_cluster = raw_ward
)
utils::write.csv(assignments, "results/cluster_assignments.csv", row.names = FALSE)

cluster_profiles <- do.call(rbind, lapply(sort(unique(assignments$scaled_kmeans_cluster)),
                                          function(cluster_id) {
  idx <- assignments$scaled_kmeans_cluster == cluster_id
  data.frame(cluster = cluster_id, n = sum(idx),
             t(colMeans(x[idx, , drop = FALSE])), check.names = FALSE)
}))
utils::write.csv(cluster_profiles, "results/cluster_profiles.csv", row.names = FALSE)

kmeans_genotype <- constrained_cluster_permutation(
  assignments$scaled_kmeans_cluster, data$genotype, data$environment,
  replicates = cluster_permutations, seed = seed + 401L
)
ward_genotype <- constrained_cluster_permutation(
  assignments$scaled_ward_cluster, data$genotype, data$environment,
  replicates = cluster_permutations, seed = seed + 402L
)

external_cluster_validation <- data.frame(
  method = c("kmeans", "ward"),
  ari_with_genotype = c(kmeans_genotype$observed, ward_genotype$observed),
  constrained_permutation_p = c(kmeans_genotype$p_value, ward_genotype$p_value),
  null_mean = c(kmeans_genotype$null_mean, ward_genotype$null_mean),
  null_q025 = c(kmeans_genotype$null_q025, ward_genotype$null_q025),
  null_q975 = c(kmeans_genotype$null_q975, ward_genotype$null_q975),
  permutations = cluster_permutations
)
utils::write.csv(external_cluster_validation,
                 "results/cluster_external_validation.csv", row.names = FALSE)

genotype_cv <- genotype_cv_permutation(x, data$genotype, data$environment,
                                       replicates = cv_permutations,
                                       seed = seed + 501L)
prediction_table <- data.frame(
  sample_id = data$sample_id,
  held_out_environment = data$environment,
  true_genotype = data$genotype,
  predicted_genotype = genotype_cv$predictions,
  correct = genotype_cv$predictions == data$genotype
)
utils::write.csv(prediction_table, "results/genotype_cv_predictions.csv",
                 row.names = FALSE)
confusion <- table(true = factor(data$genotype, levels = sort(unique(data$genotype))),
                   predicted = factor(genotype_cv$predictions,
                                      levels = sort(unique(data$genotype))))
confusion_table <- data.frame(true_genotype = rownames(confusion),
                              as.data.frame.matrix(confusion), check.names = FALSE)
utils::write.csv(confusion_table, "results/genotype_cv_confusion.csv",
                 row.names = FALSE)

cluster_method_ari <- adjusted_rand_index(assignments$scaled_kmeans_cluster,
                                          assignments$scaled_ward_cluster)
kmeans_scale_ari <- adjusted_rand_index(assignments$scaled_kmeans_cluster,
                                        assignments$unscaled_kmeans_cluster)
ward_scale_ari <- adjusted_rand_index(assignments$scaled_ward_cluster,
                                      assignments$unscaled_ward_cluster)
kmeans_stability_summary <- quantile_summary(
  cluster_stability$draws$kmeans_ari_to_full
)
ward_stability_summary <- quantile_summary(
  cluster_stability$draws$ward_ari_to_full
)
pca_pc1_cosine <- quantile_summary(pca_bootstrap$stability$pc1_axis_cosine)
pca_pc2_cosine <- quantile_summary(pca_bootstrap$stability$pc2_axis_cosine)
pca_subspace_angle <- quantile_summary(
  pca_bootstrap$stability$largest_subspace_angle_deg
)

summary <- list(
  metadata = list(
    analysis_date = "2026-09-29",
    seed = seed,
    observations = nrow(data),
    environments = length(unique(data$environment)),
    genotypes = length(unique(data$genotype)),
    interpretation = "exploratory, non-causal"
  ),
  pca = list(
    standardized_explained_variance = as.list(
      setNames((scaled_pca$sdev^2 / sum(scaled_pca$sdev^2))[1:5], paste0("PC", 1:5))
    ),
    unscaled_explained_variance = as.list(
      setNames((unscaled_pca$sdev^2 / sum(unscaled_pca$sdev^2))[1:5], paste0("PC", 1:5))
    ),
    first_two_standardized_cumulative = sum(
      (scaled_pca$sdev^2 / sum(scaled_pca$sdev^2))[1:2]
    ),
    first_two_unscaled_cumulative = sum(
      (unscaled_pca$sdev^2 / sum(unscaled_pca$sdev^2))[1:2]
    ),
    scaling_subspace_angles_degrees = unname(angles),
    block_bootstrap_replicates = pca_bootstraps,
    pc1_axis_cosine = as.list(pca_pc1_cosine),
    pc2_axis_cosine = as.list(pca_pc2_cosine),
    largest_subspace_angle_degrees = as.list(pca_subspace_angle)
  ),
  clustering = list(
    gap_bootstrap_replicates = gap_bootstraps,
    gap_first_se_selected_k_standardized = scaled_gap$selected_k,
    gap_first_se_selected_k_unscaled = unscaled_gap$selected_k,
    exploratory_forced_k_standardized = forced_k,
    exploratory_forced_k_unscaled = raw_forced_k,
    method_ari_standardized = cluster_method_ari,
    kmeans_scale_ari = kmeans_scale_ari,
    ward_scale_ari = ward_scale_ari,
    environment_subsample_replicates = cluster_subsamples,
    kmeans_stability = as.list(kmeans_stability_summary),
    ward_stability = as.list(ward_stability_summary),
    genotype_external_validation = list(
      kmeans_ari = kmeans_genotype$observed,
      kmeans_p = kmeans_genotype$p_value,
      ward_ari = ward_genotype$observed,
      ward_p = ward_genotype$p_value,
      constrained_permutations = cluster_permutations
    )
  ),
  leave_one_environment_out_genotype_validation = list(
    accuracy = genotype_cv$accuracy,
    balanced_accuracy = genotype_cv$balanced_accuracy,
    constrained_permutation_p = genotype_cv$p_value,
    null_mean_accuracy = genotype_cv$null_mean,
    null_q025 = genotype_cv$null_q025,
    null_q975 = genotype_cv$null_q975,
    permutations = cv_permutations,
    per_class_accuracy = as.list(genotype_cv$per_class_accuracy)
  )
)
jsonlite::write_json(summary, "results/summary.json", pretty = TRUE,
                     auto_unbox = TRUE, digits = 10)

palette <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9")
genotype_levels <- sort(unique(data$genotype))
genotype_color <- palette[match(data$genotype, genotype_levels)]

grDevices::png("figures/pca_biplot.png", width = 1200, height = 900,
               type = "quartz", res = 150)
scores <- scaled_pca$x[, 1:2, drop = FALSE]
arrow_scale <- 0.75 * min(diff(range(scores[, 1])), diff(range(scores[, 2])))
load2 <- correlation_loadings(scaled_pca)[, 1:2, drop = FALSE]
arrow_x <- load2[, 1] * arrow_scale
arrow_y <- load2[, 2] * arrow_scale
graphics::plot(scores, col = genotype_color, pch = 19,
               xlim = range(c(scores[, 1], 1.18 * arrow_x)),
               ylim = range(c(scores[, 2], 1.18 * arrow_y)),
               xlab = sprintf("PC1 (%.1f%%)", 100 * summary$pca$standardized_explained_variance$PC1),
               ylab = sprintf("PC2 (%.1f%%)", 100 * summary$pca$standardized_explained_variance$PC2),
               main = "Standardized sensory PCA")
graphics::arrows(0, 0, load2[, 1] * arrow_scale, load2[, 2] * arrow_scale,
                 length = 0.08, col = "#333333")
graphics::text(load2[, 1] * arrow_scale, load2[, 2] * arrow_scale,
               labels = gsub("_", " ", rownames(load2)), pos = 3, cex = 0.75)
graphics::legend("bottomleft", legend = genotype_levels, col = palette,
                 pch = 19, bty = "n", cex = 0.75)
grDevices::dev.off()

grDevices::png("figures/pca_loading_stability.png", width = 1350, height = 750,
               type = "quartz", res = 150)
graphics::par(mfrow = c(1, 2), mar = c(4, 7, 3, 1))
for (component in c("PC1", "PC2")) {
  part <- pca_loading_intervals[pca_loading_intervals$component == component, ]
  part <- part[order(part$mean), ]
  y <- seq_len(nrow(part))
  graphics::plot(part$mean, y, xlim = range(c(part$lower_95, part$upper_95, 0)),
                 yaxt = "n", ylab = "", xlab = "Correlation loading",
                 main = paste(component, "environment-block bootstrap"), pch = 19)
  graphics::axis(2, at = y, labels = gsub("_", " ", part$variable), las = 1)
  graphics::abline(v = 0, lty = 2, col = "#777777")
  graphics::segments(part$lower_95, y, part$upper_95, y, lwd = 2,
                     col = "#0072B2")
  graphics::points(part$mean, y, pch = 19)
}
grDevices::dev.off()

grDevices::png("figures/cluster_diagnostics.png", width = 1500, height = 1050,
               type = "quartz", res = 150)
graphics::par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
for (mode in c("z_scaled", "unscaled")) {
  part <- internal_validation[internal_validation$scaling == mode, ]
  yrange <- range(part$mean_silhouette)
  graphics::plot(NA, xlim = range(part$k), ylim = yrange,
                 xlab = "k", ylab = "Mean silhouette",
                 main = paste("Silhouette:", mode))
  for (i in seq_along(c("kmeans", "ward"))) {
    method <- c("kmeans", "ward")[i]
    m <- part[part$method == method, ]
    graphics::lines(m$k, m$mean_silhouette, type = "b", pch = 18 + i,
                    col = c("#0072B2", "#D55E00")[i])
  }
  graphics::legend("bottomleft", c("k-means", "Ward"),
                   col = c("#0072B2", "#D55E00"), pch = c(19, 20),
                   lty = 1, bty = "n")
}
for (mode in c("z_scaled", "unscaled")) {
  part <- gap_rows[gap_rows$scaling == mode, ]
  graphics::plot(part$k, part$gap, type = "b", pch = 19,
                 ylim = range(c(part$gap - part$gap_standard_error,
                                part$gap + part$gap_standard_error)),
                 xlab = "k", ylab = "Gap statistic",
                 main = paste("Gap statistic:", mode))
  graphics::segments(part$k, part$gap - part$gap_standard_error,
                     part$k, part$gap + part$gap_standard_error)
}
grDevices::dev.off()

grDevices::png("figures/exploratory_clusters.png", width = 1200, height = 900,
               type = "quartz", res = 150)
cluster_colors <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7",
                    "#E69F00", "#56B4E9", "#000000", "#999999")
graphics::plot(scores, col = cluster_colors[assignments$scaled_kmeans_cluster],
               pch = 14 + match(data$genotype, genotype_levels),
               xlab = "Standardized PC1", ylab = "Standardized PC2",
               main = paste0("Exploratory k-means solution (k = ", forced_k, ")"))
graphics::legend("topright", legend = paste("cluster", seq_len(forced_k)),
                 col = cluster_colors[seq_len(forced_k)], pch = 19,
                 bty = "n", cex = 0.8)
graphics::legend("bottomleft", legend = genotype_levels,
                 pch = 14 + seq_along(genotype_levels), bty = "n", cex = 0.7)
grDevices::dev.off()

grDevices::png("figures/genotype_cv_confusion.png", width = 1050, height = 900,
               type = "quartz", res = 150)
graphics::par(mar = c(8, 8, 3, 2))
confusion_numeric <- unclass(confusion)
graphics::image(seq_len(nrow(confusion_numeric)), seq_len(ncol(confusion_numeric)),
                t(confusion_numeric), col = hcl.colors(12, "Blues 3", rev = TRUE),
                xaxt = "n", yaxt = "n", xlab = "", ylab = "",
                main = "Leave-one-environment-out predictions")
graphics::axis(1, at = seq_along(genotype_levels), labels = genotype_levels,
               las = 2, cex.axis = 0.68)
graphics::axis(2, at = seq_along(genotype_levels), labels = genotype_levels,
               las = 2, cex.axis = 0.68)
graphics::mtext("True genotype", side = 1, line = 6)
graphics::mtext("Predicted genotype", side = 2, line = 6)
for (i in seq_len(nrow(confusion_numeric))) for (j in seq_len(ncol(confusion_numeric))) {
  graphics::text(i, j, confusion_numeric[i, j], cex = 0.9)
}
grDevices::dev.off()

fmt <- function(value, digits = 3) formatC(value, digits = digits, format = "f")
gap_conclusion <- if (scaled_gap$selected_k == 1L) {
  "the first-SE rule selects `k=1`, so the primary analysis does not support a discrete partition"
} else {
  paste0("the first-SE rule selects `k=", scaled_gap$selected_k, "`")
}

pc1_rows <- pca_loadings[pca_loadings$scaling == "z_scaled" &
                           pca_loadings$component == "PC1", ]
pc2_rows <- pca_loadings[pca_loadings$scaling == "z_scaled" &
                           pca_loadings$component == "PC2", ]
loading_table <- c(
  "| Sensory variable | PC1 correlation loading | PC2 correlation loading |",
  "|---|---:|---:|",
  vapply(traits, function(variable) {
    sprintf("| %s | %s | %s |", gsub("_", " ", variable),
            fmt(pc1_rows$correlation_loading[pc1_rows$variable == variable]),
            fmt(pc2_rows$correlation_loading[pc2_rows$variable == variable]))
  }, character(1))
)

report <- c(
  "# ブルーベリー官能評価の多変量解析 / Multivariate Analysis of Blueberry Sensory Data",
  "",
  "## Abstract",
  "",
  sprintf(paste0(
    "This study analyzes the open S2 supporting data from Gilbert et al. (2015): %d location-by-harvest-by-genotype samples from 2013, %d environment blocks, %d genotypes, and five sensory ratings. ",
    "The first two standardized principal components explain %.1f%% of total variance, compared with %.1f%% without scaling. ",
    "Environment-block resampling, scale sensitivity, gap statistics, silhouettes, cluster subsampling, and leave-one-environment-out validation distinguish reproducible continuous structure from weak evidence for discrete groups."),
    nrow(data), length(unique(data$environment)), length(unique(data$genotype)),
    100 * summary$pca$first_two_standardized_cumulative,
    100 * summary$pca$first_two_unscaled_cumulative),
  "",
  "The central result is conservative: the gap statistic's first-SE rule selects `k=1`, so the six known genotypes must not be equated with six natural sensory clusters. When at least two clusters are forced, mean silhouette selects an exploratory `k=4` solution, but this is only a descriptive partition. Separately, genotype labels that never enter feature construction achieve above-permutation accuracy in leave-one-environment-out nearest-centroid validation. The five ratings therefore contain genotype-related signal that generalizes across environments, but that signal does not naturally form six compact groups.",
  "",
  "## 1. Research questions",
  "",
  "1. What are the main continuous axes of the five sensory ratings in 2013?",
  "2. Does scaling materially change the PCA or clustering conclusions?",
  "3. Are PCA loadings stable under environment-block resampling?",
  "4. How many discrete groups does the data support, and how stable are they across methods, environment subsets, and scaling choices?",
  "5. Without using genotype labels to select `k`, can sensory features identify genotype after holding out complete environments?",
  "",
  "This is an exploratory measurement study. It estimates no treatment effect and makes no causal claim.",
  "",
  "## 2. Data source and license",
  "",
  "The public data are from Jessica L. Gilbert et al. (2015), *Identifying Breeding Priorities for Blueberry Flavor Using Biochemical, Sensory, and Genotype by Environment Analyses*, PLOS ONE 10(9): e0138494, DOI: [10.1371/journal.pone.0138494](https://doi.org/10.1371/journal.pone.0138494). The article and supporting information are published under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/). `scripts/download_source.R` pins the S2 XLSX URL and SHA-256; `scripts/prepare_data.R` selects 2013 and the five sensory fields from `S2 Table`.",
  "",
  "The workflow treats scaling as an explicit sensitivity question, selects the cluster count without genotype labels, and quantifies loading and clustering stability before genotype is used for external checks.",
  "",
  "## 3. Design and methods",
  "",
  "### 3.1 Repeated structure",
  "",
  sprintf("The 50 observations come from %d environment blocks formed by location C/H/W and harvest. Genotypes recur across environments, so an ordinary row bootstrap would overstate independent information. PCA uncertainty uses %d environment-block bootstrap replicates: each replicate samples nine complete blocks with replacement, refits the mean, standard deviation, and PCA, and aligns the first two axes by permutation and sign.", length(unique(data$environment)), pca_bootstraps),
  "",
  "### 3.2 PCA and scaling",
  "",
  "The primary analysis z-standardizes each of the five traits. A sensitivity analysis retains the original scales: the fields use related rating scales, but their empirical variances differ. The report includes eigenvector loadings, variable-component correlation loadings, explained variance, score correlations, and angles between the first two subspaces.",
  "",
  "### 3.3 Cluster count, stability, and external labels",
  "",
  sprintf(paste0(
    "The gap statistic uses %d reference bootstraps for `k=1…8`; mean silhouettes are computed for k-means and Ward.D2 over `k=2…8`. ",
    "The primary rule is the gap first-SE criterion. If it returns `k=1`, a silhouette-selected solution with `k≥2` is reported for visualization only. ",
    "Stability uses %d environment subsamples, each retaining %d/%d complete blocks without replacement (`ceiling(0.8 × %d)`). ",
    "ARI also measures agreement between k-means and Ward and between standardized and unscaled solutions."),
    gap_bootstraps, cluster_subsamples,
    ceiling(0.8 * length(unique(data$environment))),
    length(unique(data$environment)), length(unique(data$environment))),
  "",
  sprintf(paste0(
    "Genotype labels never enter PCA, cluster-count selection, or clustering. External checks are: ",
    "(1) cluster-genotype ARI with %d genotype-label permutations constrained within environment blocks; and ",
    "(2) leave-one-environment-out nearest-centroid classification, with scaling and genotype centroids estimated only from the remaining environments and %d similarly constrained permutations."),
    cluster_permutations, cv_permutations),
  "",
  "## 4. Results",
  "",
  "### 4.1 PCA",
  "",
  sprintf("Standardized PC1 and PC2 explain %s%% and %s%%, or %s%% cumulatively; the first two unscaled axes explain %s%% and %s%%.",
          fmt(100 * summary$pca$standardized_explained_variance$PC1, 1),
          fmt(100 * summary$pca$standardized_explained_variance$PC2, 1),
          fmt(100 * summary$pca$first_two_standardized_cumulative, 1),
          fmt(100 * summary$pca$unscaled_explained_variance$PC1, 1),
          fmt(100 * summary$pca$unscaled_explained_variance$PC2, 1)),
  "",
  loading_table,
  "",
  sprintf(paste0(
    "PC1 aligns overall liking, texture, sweetness, and flavor, while sourness is near zero; PC2 is dominated by positive sourness variation. ",
    "In the environment-block bootstrap, median cosines between the resampled and full-sample PC1/PC2 axes are %s and %s. ",
    "The median maximum angle between the first two subspaces is %s° (95%% interval %s°–%s°). ",
    "This quantifies axis uncertainty more directly than a single biplot."),
    fmt(pca_pc1_cosine["median"]), fmt(pca_pc2_cosine["median"]),
    fmt(pca_subspace_angle["median"], 1), fmt(pca_subspace_angle["q025"], 1),
    fmt(pca_subspace_angle["q975"], 1)),
  "",
  sprintf("The first two standardized and unscaled subspace angles are %s° and %s°; corresponding score correlations are %s and %s. The raw scale gives sourness substantially more influence on the first axis, so scaling cannot be treated as inconsequential.",
          fmt(angles[1], 1), fmt(angles[2], 1),
          fmt(scale_sensitivity$value[scale_sensitivity$metric == "pc1_score_correlation"]),
          fmt(scale_sensitivity$value[scale_sensitivity$metric == "pc2_score_correlation"])),
  "",
  "### 4.2 Evidence for discrete clusters",
  "",
  sprintf("For standardized data, %s. The corresponding rule selects `k=%d` without scaling. The primary evidence therefore favors a continuous sensory space rather than clearly separated natural clusters.", gap_conclusion, unscaled_gap$selected_k),
  "",
  sprintf(paste0(
    "After imposing `k≥2`, mean silhouette selects `k=%d` for standardized data. K-means and Ward agree with ARI %s. ",
    "Across environment subsamples, median ARI relative to the full-sample solution is %s for k-means (95%% interval %s–%s) ",
    "and %s for Ward (%s–%s). Standardized/unscaled agreement is %s for k-means and %s for Ward."),
    forced_k, fmt(cluster_method_ari),
    fmt(kmeans_stability_summary["median"]), fmt(kmeans_stability_summary["q025"]),
    fmt(kmeans_stability_summary["q975"]), fmt(ward_stability_summary["median"]),
    fmt(ward_stability_summary["q025"]), fmt(ward_stability_summary["q975"]),
    fmt(kmeans_scale_ari), fmt(ward_scale_ari)),
  "",
  sprintf(paste0(
    "The exploratory `k=%d` clusters agree weakly with the six genotype labels: k-means ARI=%s (within-environment permutation p=%s) ",
    "and Ward ARI=%s (p=%s). The clusters therefore cannot be named as genotype groups, and the known count of six genotypes cannot be used to justify six clusters."),
    forced_k, fmt(kmeans_genotype$observed), fmt(kmeans_genotype$p_value, 4),
    fmt(ward_genotype$observed), fmt(ward_genotype$p_value, 4)),
  "",
  "### 4.3 Genotype signal across environments",
  "",
  sprintf(paste0(
    "Leave-one-environment-out nearest-centroid accuracy is %s%% and balanced accuracy is %s%%. Constrained permutations have mean accuracy %s%%, ",
    "a 95%% interval of %s%%–%s%%, and Monte Carlo p=%s. ",
    "The sensory vector therefore contains genotype information that generalizes across environments, but this may be an overlapping continuous displacement and does not require six separated clusters."),
    fmt(100 * genotype_cv$accuracy, 1), fmt(100 * genotype_cv$balanced_accuracy, 1),
    fmt(100 * genotype_cv$null_mean, 1), fmt(100 * genotype_cv$null_q025, 1),
    fmt(100 * genotype_cv$null_q975, 1), fmt(genotype_cv$p_value, 4)),
  "",
  "## 5. Conclusion",
  "",
  "The 2013 blueberry sensory data have a clear, reproducible, low-dimensional continuous structure: one combined liking/sweetness/flavor axis and one sourness-dominated axis. Genotypes show sensory differences that predict across environments, but unsupervised evidence does not support a natural partition into six genotype clusters. The defensible conclusion is continuous genotype-related sensory signal; discrete clusters are a weaker exploratory summary that is sensitive to method, environment, and scaling.",
  "",
  "## 6. Limitations",
  "",
  "- Only nine environment blocks are available, so the tails of block-bootstrap intervals are themselves imprecise.",
  "- Genotypes are not fully balanced in 2013: Primadonna lacks H3 and Scintilla lacks three W environments.",
  "- The table contains sample-level sensory summaries, not individual ratings from all 217 consumers, so participant-level random effects and demographic heterogeneity cannot be reconstructed.",
  "- Nearest-centroid validation establishes signal in the five ratings, not genetic causality; location, harvest time, and their interactions may still affect ratings.",
  "- Gap statistics, silhouettes, and ARI answer different questions; exploratory `k=4` is not an estimate of a true cluster count.",
  "- Results cover 2013 and six genotypes repeated in this subset, not all 19 genotypes or other years.",
  "",
  "## 7. Reproducibility and artifacts",
  "",
  "Rebuild completely from the public source:",
  "",
  "```bash",
  "make reproduce",
  "```",
  "",
  "Re-run offline from the committed, attributed 50-row derived snapshot:",
  "",
  "```bash",
  "make all",
  "```",
  "",
  "The random seed is `20260929`. `results/summary.json` provides machine-readable conclusions; every table, PNG, session record, and SHA-256 manifest is generated by the script. See `DATA_CARD.md` and `provenance/SOURCE_AUDIT.md` for source and publication boundaries."
)
writeLines(report, "REPORT.md", useBytes = TRUE)

capture.output(utils::sessionInfo(), file = "results/session_info.txt")

sha256 <- function(path) {
  output <- system2("shasum", c("-a", "256", shQuote(path)), stdout = TRUE)
  sub("[[:space:]].*$", "", output[[1]])
}
artifact_paths <- c(
  "REPORT.md",
  list.files("results", full.names = TRUE),
  list.files("figures", full.names = TRUE),
  "data/blueberry_2013_sensory.csv", "data/source_manifest.json"
)
artifact_paths <- artifact_paths[file.exists(artifact_paths)]
artifact_paths <- artifact_paths[basename(artifact_paths) != "artifact_manifest.csv"]
artifact_manifest <- data.frame(
  path = artifact_paths,
  bytes = unname(file.info(artifact_paths)$size),
  sha256 = vapply(artifact_paths, sha256, character(1))
)
utils::write.csv(artifact_manifest, "results/artifact_manifest.csv", row.names = FALSE)

cat("Analysis complete:", nrow(data), "observations; standardized gap k =",
    scaled_gap$selected_k, "; exploratory k =", forced_k,
    "; LOEO accuracy =", sprintf("%.3f", genotype_cv$accuracy), "\n")
