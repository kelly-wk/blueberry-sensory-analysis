options(stringsAsFactors = FALSE)

if (!requireNamespace("readxl", quietly = TRUE)) stop("Package readxl is required")
if (!requireNamespace("jsonlite", quietly = TRUE)) stop("Package jsonlite is required")

raw_path <- "data/raw/pone.0138494.s002.xlsx"
if (!file.exists(raw_path)) {
  status <- system2("Rscript", c("--vanilla", "scripts/download_source.R", raw_path))
  if (status != 0) stop("Source download failed")
}

raw <- readxl::read_excel(raw_path, sheet = "S2 Table", skip = 3,
                          .name_repair = "unique_quiet")
traits <- c("Overall Liking", "Texture", "Sweetness", "Sourness", "Flavor")
year <- suppressWarnings(as.numeric(raw$Year))
numeric_traits <- lapply(raw[traits], function(x) suppressWarnings(as.numeric(x)))

data <- data.frame(
  sample_id = raw$gID,
  year = year,
  location = raw$Location,
  harvest = suppressWarnings(as.integer(raw[["Harvest#"]])),
  environment = paste0(raw$Location, raw[["Harvest#"]]),
  genotype = raw$Genotype,
  overall_liking = numeric_traits[[1]],
  texture = numeric_traits[[2]],
  sweetness = numeric_traits[[3]],
  sourness = numeric_traits[[4]],
  flavor = numeric_traits[[5]],
  stringsAsFactors = FALSE
)
data <- data[!is.na(data$sample_id) & data$year == 2013 &
               complete.cases(data), ]
data[c("overall_liking", "texture", "sweetness", "sourness", "flavor")] <-
  lapply(data[c("overall_liking", "texture", "sweetness", "sourness", "flavor")],
         function(x) round(x, 1))
data <- data[order(data$location, data$harvest, data$genotype), ]
row.names(data) <- NULL

stopifnot(nrow(data) == 50L,
          length(unique(data$sample_id)) == 50L,
          length(unique(data$environment)) == 9L,
          identical(sort(unique(data$genotype)),
                    c("Emerald", "Endura", "Farthing", "Meadowlark",
                      "Primadonna", "Scintilla")))

dir.create("data", showWarnings = FALSE)
utils::write.csv(data, "data/blueberry_2013_sensory.csv",
                 row.names = FALSE, quote = TRUE, na = "")

manifest <- list(
  source = "Gilbert et al. (2015), Supporting Information S2 Table",
  article_doi = "10.1371/journal.pone.0138494",
  supplement_doi = "10.1371/journal.pone.0138494.s002",
  source_url = paste0(
    "https://journals.plos.org/plosone/article/file?type=supplementary&",
    "id=info:doi/10.1371/journal.pone.0138494.s002"
  ),
  source_sha256 = "fade6546de831bf15c5b4ef0fa234b13d25764a599ae87a04f906bbf051acab6",
  source_bytes = 94712,
  retrieved = "2026-09-29",
  license = "Creative Commons Attribution 4.0 International (CC BY 4.0)",
  transformation = list(
    sheet = "S2 Table",
    header_rows_skipped = 3,
    filters = c("Year == 2013", "non-missing gID", "complete five sensory outcomes"),
    selected_fields = names(data),
    output_rows = nrow(data)
  )
)
jsonlite::write_json(manifest, "data/source_manifest.json",
                     pretty = TRUE, auto_unbox = TRUE)
cat("Prepared", nrow(data), "rows in data/blueberry_2013_sensory.csv\n")

