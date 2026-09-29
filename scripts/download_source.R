args <- commandArgs(trailingOnly = TRUE)
output <- if (length(args)) args[[1]] else "data/raw/pone.0138494.s002.xlsx"

source_url <- paste0(
  "https://journals.plos.org/plosone/article/file?type=supplementary&",
  "id=info:doi/10.1371/journal.pone.0138494.s002"
)
expected_sha256 <- "fade6546de831bf15c5b4ef0fa234b13d25764a599ae87a04f906bbf051acab6"

sha256 <- function(path) {
  result <- system2("shasum", c("-a", "256", shQuote(path)), stdout = TRUE)
  sub("[[:space:]].*$", "", result[[1]])
}

dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
if (!file.exists(output)) {
  temporary <- paste0(output, ".download")
  on.exit(if (file.exists(temporary)) unlink(temporary), add = TRUE)
  utils::download.file(source_url, temporary, mode = "wb", quiet = TRUE)
  if (!identical(sha256(temporary), expected_sha256)) {
    stop("Downloaded S2 workbook failed SHA-256 verification")
  }
  if (!file.rename(temporary, output)) stop("Could not move verified download")
}
if (!identical(sha256(output), expected_sha256)) {
  stop("Existing S2 workbook has the wrong SHA-256; remove it and retry")
}
cat(normalizePath(output), expected_sha256, "\n")

