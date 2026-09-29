sensory_variables <- function() {
  c("overall_liking", "texture", "sweetness", "sourness", "flavor")
}

assert_sensory_data <- function(data) {
  required <- c("sample_id", "year", "location", "harvest", "environment",
                "genotype", sensory_variables())
  missing_columns <- setdiff(required, names(data))
  if (length(missing_columns)) {
    stop("Missing columns: ", paste(missing_columns, collapse = ", "))
  }
  if (anyDuplicated(data$sample_id)) stop("sample_id must be unique")
  if (any(!complete.cases(data[required]))) stop("Analysis data contain missing values")
  if (any(!vapply(data[sensory_variables()], is.numeric, logical(1)))) {
    stop("All sensory variables must be numeric")
  }
  if (any(as.matrix(data[sensory_variables()]) < 0)) {
    stop("Sensory values must be non-negative")
  }
  invisible(TRUE)
}

choose_two <- function(x) x * (x - 1) / 2

adjusted_rand_index <- function(labels_a, labels_b) {
  if (length(labels_a) != length(labels_b)) stop("Label vectors have different lengths")
  if (length(labels_a) < 2L) return(NA_real_)
  tab <- table(labels_a, labels_b)
  n <- sum(tab)
  total_pairs <- choose_two(n)
  sum_rows <- sum(choose_two(rowSums(tab)))
  sum_cols <- sum(choose_two(colSums(tab)))
  observed <- sum(choose_two(tab))
  expected <- sum_rows * sum_cols / total_pairs
  denominator <- 0.5 * (sum_rows + sum_cols) - expected
  if (abs(denominator) < .Machine$double.eps) {
    return(if (all(as.integer(factor(labels_a)) == as.integer(factor(labels_b)))) 1 else 0)
  }
  (observed - expected) / denominator
}

fit_pca <- function(x, scale_data = TRUE) {
  if (!is.matrix(x)) x <- as.matrix(x)
  if (any(!is.finite(x))) stop("PCA matrix contains non-finite values")
  if (any(apply(x, 2, stats::sd) == 0)) stop("PCA matrix contains a constant column")
  stats::prcomp(x, center = TRUE, scale. = scale_data)
}

correlation_loadings <- function(fit) {
  sweep(fit$rotation, 2, fit$sdev, `*`)
}

align_two_axes <- function(candidate, reference) {
  if (!all(dim(candidate) == dim(reference)) || ncol(candidate) != 2L) {
    stop("candidate and reference must be equal p x 2 matrices")
  }
  identity_score <- sum(abs(diag(crossprod(candidate, reference))))
  swapped <- candidate[, 2:1, drop = FALSE]
  swap_score <- sum(abs(diag(crossprod(swapped, reference))))
  permutation <- if (swap_score > identity_score) c(2L, 1L) else c(1L, 2L)
  aligned <- candidate[, permutation, drop = FALSE]
  signs <- sign(diag(crossprod(aligned, reference)))
  signs[signs == 0] <- 1
  aligned <- sweep(aligned, 2, signs, `*`)
  list(matrix = aligned, permutation = permutation, signs = signs)
}

principal_angles_degrees <- function(candidate, reference) {
  singular_values <- svd(crossprod(candidate, reference), nu = 0, nv = 0)$d
  singular_values <- pmin(1, pmax(0, singular_values))
  sort(acos(singular_values) * 180 / pi)
}

best_kmeans <- function(x, centers, seed, nstart = 100L) {
  set.seed(seed)
  stats::kmeans(x, centers = centers, nstart = nstart, iter.max = 1000)
}

mean_silhouette <- function(labels, distances) {
  if (length(unique(labels)) < 2L) return(NA_real_)
  mean(cluster::silhouette(labels, distances)[, "sil_width"])
}

cluster_internal_validation <- function(x, k_values = 2:8, seed = 20260929L) {
  distances <- stats::dist(x)
  ward <- stats::hclust(distances, method = "ward.D2")
  rows <- list()
  for (k in k_values) {
    km <- best_kmeans(x, k, seed + k)
    rows[[length(rows) + 1L]] <- data.frame(
      k = k, method = "kmeans",
      mean_silhouette = mean_silhouette(km$cluster, distances)
    )
    rows[[length(rows) + 1L]] <- data.frame(
      k = k, method = "ward",
      mean_silhouette = mean_silhouette(stats::cutree(ward, k), distances)
    )
  }
  do.call(rbind, rows)
}

gap_statistic <- function(x, k_max = 8L, bootstraps = 200L, seed = 20260929L) {
  set.seed(seed)
  gap <- cluster::clusGap(
    x,
    FUNcluster = function(z, k) {
      list(cluster = stats::kmeans(z, centers = k, nstart = 50L,
                                   iter.max = 1000)$cluster)
    },
    K.max = k_max,
    B = bootstraps,
    verbose = FALSE
  )
  table <- as.data.frame(gap$Tab)
  table$k <- seq_len(nrow(table))
  selected <- cluster::maxSE(table$gap, table$SE.sim, method = "firstSEmax")
  list(table = table, selected_k = as.integer(selected))
}

block_bootstrap_pca <- function(x, blocks, replicates = 1000L,
                                seed = 20260929L) {
  reference <- fit_pca(x, scale_data = TRUE)
  reference_rotation <- reference$rotation[, 1:2, drop = FALSE]
  block_levels <- unique(blocks)
  loading_rows <- vector("list", replicates)
  stability_rows <- vector("list", replicates)
  set.seed(seed)
  for (b in seq_len(replicates)) {
    sampled_blocks <- sample(block_levels, length(block_levels), replace = TRUE)
    indices <- unlist(lapply(sampled_blocks, function(z) which(blocks == z)),
                      use.names = FALSE)
    fitted <- fit_pca(x[indices, , drop = FALSE], scale_data = TRUE)
    axis_alignment <- align_two_axes(fitted$rotation[, 1:2, drop = FALSE],
                                     reference_rotation)
    correlations <- correlation_loadings(fitted)[, 1:2, drop = FALSE]
    correlations <- correlations[, axis_alignment$permutation, drop = FALSE]
    correlations <- sweep(correlations, 2, axis_alignment$signs, `*`)
    loading_rows[[b]] <- data.frame(
      replicate = b,
      variable = rep(rownames(correlations), 2L),
      component = rep(c("PC1", "PC2"), each = nrow(correlations)),
      correlation_loading = as.vector(correlations)
    )
    aligned_rotation <- axis_alignment$matrix
    cosines <- diag(crossprod(aligned_rotation, reference_rotation))
    angles <- principal_angles_degrees(fitted$rotation[, 1:2, drop = FALSE],
                                        reference_rotation)
    pve <- fitted$sdev^2 / sum(fitted$sdev^2)
    matched_pve <- pve[axis_alignment$permutation]
    stability_rows[[b]] <- data.frame(
      replicate = b,
      pc1_axis_cosine = cosines[1],
      pc2_axis_cosine = cosines[2],
      pc1_explained = matched_pve[1],
      pc2_explained = matched_pve[2],
      largest_subspace_angle_deg = max(angles)
    )
  }
  list(loadings = do.call(rbind, loading_rows),
       stability = do.call(rbind, stability_rows))
}

summarize_loading_bootstrap <- function(draws) {
  split_draws <- split(draws$correlation_loading,
                       interaction(draws$variable, draws$component, drop = TRUE))
  rows <- lapply(names(split_draws), function(key) {
    parts <- strsplit(key, "\\.", fixed = FALSE)[[1]]
    values <- split_draws[[key]]
    data.frame(
      variable = paste(parts[-length(parts)], collapse = "."),
      component = parts[length(parts)],
      mean = mean(values),
      median = stats::median(values),
      lower_95 = unname(stats::quantile(values, 0.025)),
      upper_95 = unname(stats::quantile(values, 0.975)),
      sign_agreement = max(mean(values >= 0), mean(values <= 0))
    )
  })
  out <- do.call(rbind, rows)
  out[order(out$component, out$variable), ]
}

block_subsample_cluster_stability <- function(x, blocks, k, replicates = 500L,
                                              seed = 20260929L) {
  z <- scale(x)
  baseline_kmeans <- best_kmeans(z, k, seed)
  baseline_ward <- stats::cutree(stats::hclust(stats::dist(z), method = "ward.D2"), k)
  block_levels <- unique(blocks)
  number_blocks <- max(2L, ceiling(0.8 * length(block_levels)))
  rows <- vector("list", replicates)
  set.seed(seed + 1000L)
  for (b in seq_len(replicates)) {
    kept_blocks <- sample(block_levels, number_blocks, replace = FALSE)
    indices <- which(blocks %in% kept_blocks)
    z_sub <- scale(x[indices, , drop = FALSE])
    km <- best_kmeans(z_sub, k, seed + 2000L + b, nstart = 50L)$cluster
    ward <- stats::cutree(stats::hclust(stats::dist(z_sub), method = "ward.D2"), k)
    rows[[b]] <- data.frame(
      replicate = b,
      environments_retained = number_blocks,
      observations_retained = length(indices),
      kmeans_ari_to_full = adjusted_rand_index(km, baseline_kmeans$cluster[indices]),
      ward_ari_to_full = adjusted_rand_index(ward, baseline_ward[indices])
    )
  }
  list(
    draws = do.call(rbind, rows),
    baseline_kmeans = baseline_kmeans$cluster,
    baseline_ward = baseline_ward
  )
}

constrained_cluster_permutation <- function(cluster_labels, genotype, blocks,
                                            replicates = 4999L,
                                            seed = 20260929L) {
  observed <- adjusted_rand_index(cluster_labels, genotype)
  null <- numeric(replicates)
  block_levels <- unique(blocks)
  set.seed(seed)
  for (b in seq_len(replicates)) {
    permuted <- as.character(genotype)
    for (block in block_levels) {
      idx <- which(blocks == block)
      permuted[idx] <- sample(permuted[idx])
    }
    null[b] <- adjusted_rand_index(cluster_labels, permuted)
  }
  list(observed = observed,
       p_value = (1 + sum(null >= observed)) / (replicates + 1),
       null_mean = mean(null),
       null_q025 = unname(stats::quantile(null, 0.025)),
       null_q975 = unname(stats::quantile(null, 0.975)),
       replicates = replicates)
}

leave_one_environment_out <- function(x, genotype, blocks) {
  predictions <- rep(NA_character_, nrow(x))
  block_levels <- unique(blocks)
  for (block in block_levels) {
    test <- blocks == block
    train <- !test
    means <- colMeans(x[train, , drop = FALSE])
    scales <- apply(x[train, , drop = FALSE], 2, stats::sd)
    train_z <- sweep(sweep(x[train, , drop = FALSE], 2, means), 2, scales, `/`)
    test_z <- sweep(sweep(x[test, , drop = FALSE], 2, means), 2, scales, `/`)
    groups <- split(seq_len(sum(train)), genotype[train])
    centroids <- do.call(rbind, lapply(groups, function(idx) {
      colMeans(train_z[idx, , drop = FALSE])
    }))
    distances <- vapply(seq_len(nrow(test_z)), function(i) {
      rowSums((centroids - matrix(test_z[i, ], nrow = nrow(centroids),
                                  ncol = ncol(centroids), byrow = TRUE))^2)
    }, numeric(nrow(centroids)))
    if (is.null(dim(distances))) distances <- matrix(distances, ncol = 1L)
    predictions[test] <- rownames(centroids)[apply(distances, 2, which.min)]
  }
  predictions
}

genotype_cv_permutation <- function(x, genotype, blocks, replicates = 999L,
                                    seed = 20260929L) {
  observed_prediction <- leave_one_environment_out(x, genotype, blocks)
  observed_accuracy <- mean(observed_prediction == genotype)
  per_class <- tapply(observed_prediction == genotype, genotype, mean)
  null <- numeric(replicates)
  block_levels <- unique(blocks)
  set.seed(seed)
  for (b in seq_len(replicates)) {
    permuted <- as.character(genotype)
    for (block in block_levels) {
      idx <- which(blocks == block)
      permuted[idx] <- sample(permuted[idx])
    }
    null_prediction <- leave_one_environment_out(x, permuted, blocks)
    null[b] <- mean(null_prediction == permuted)
  }
  list(
    predictions = observed_prediction,
    accuracy = observed_accuracy,
    balanced_accuracy = mean(per_class),
    per_class_accuracy = per_class,
    p_value = (1 + sum(null >= observed_accuracy)) / (replicates + 1),
    null_mean = mean(null),
    null_q025 = unname(stats::quantile(null, 0.025)),
    null_q975 = unname(stats::quantile(null, 0.975)),
    replicates = replicates
  )
}

quantile_summary <- function(x) {
  c(mean = mean(x), median = stats::median(x),
    q025 = unname(stats::quantile(x, 0.025)),
    q975 = unname(stats::quantile(x, 0.975)))
}

