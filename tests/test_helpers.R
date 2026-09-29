source("R/analysis_helpers.R")

stopifnot(adjusted_rand_index(c(1, 1, 2, 2), c("a", "a", "b", "b")) == 1)
stopifnot(adjusted_rand_index(c(1, 1, 2, 2), c("a", "b", "a", "b")) < 0)

reference <- diag(3)[, 1:2]
candidate <- cbind(-reference[, 2], reference[, 1])
aligned <- align_two_axes(candidate, reference)
stopifnot(max(abs(aligned$matrix - reference)) < 1e-12)

angles <- principal_angles_degrees(reference, reference)
stopifnot(max(abs(angles)) < 1e-7)

fixture <- data.frame(
  sample_id = paste0("s", 1:6), year = 2013,
  location = rep(c("A", "B"), each = 3), harvest = 1,
  environment = rep(c("A1", "B1"), each = 3),
  genotype = rep(c("g1", "g2", "g3"), 2),
  overall_liking = 1:6, texture = 2:7, sweetness = c(1, 3, 2, 6, 5, 7),
  sourness = 7:2, flavor = c(2, 4, 3, 7, 6, 8)
)
stopifnot(isTRUE(assert_sensory_data(fixture)))
pca <- fit_pca(as.matrix(fixture[sensory_variables()]), scale_data = TRUE)
stopifnot(ncol(pca$rotation) == 5L, abs(sum(pca$sdev^2 / sum(pca$sdev^2)) - 1) < 1e-12)

bad <- fixture
bad$flavor[1] <- NA
stopifnot(inherits(try(assert_sensory_data(bad), silent = TRUE), "try-error"))

cat("test_helpers.R: all tests passed\n")
