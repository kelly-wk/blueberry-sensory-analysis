required <- c(readxl = "1.4.5", cluster = "2.1.8.2", jsonlite = "2.0.0")
failures <- character()
for (package in names(required)) {
  if (!requireNamespace(package, quietly = TRUE)) {
    failures <- c(failures, paste(package, "is missing"))
  } else {
    installed <- as.character(utils::packageVersion(package))
    cat(package, installed, "\n")
    if (installed != required[[package]]) {
      failures <- c(failures, paste(package, installed, "!=", required[[package]]))
    }
  }
}
cat("R", as.character(getRversion()), "\n")
if (length(failures)) stop(paste(failures, collapse = "; "))

