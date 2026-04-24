args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args) >= 1) args[[1]] else "."

files <- c(
  list.files(file.path(root, "R"), pattern = "\\.[Rr]$", recursive = TRUE, full.names = TRUE),
  list.files(file.path(root, "tests", "testthat"), pattern = "\\.[Rr]$", recursive = TRUE, full.names = TRUE),
  list.files(file.path(root, "data-raw"), pattern = "\\.[Rr]$", recursive = TRUE, full.names = TRUE)
)
files <- unique(files[file.exists(files)])

if (!length(files)) {
  cat("No .R files found.\n")
  quit(status = 0)
}

trim_preview <- function(x, width = 120) {
  x <- gsub("\t", " ", x, fixed = TRUE)
  if (nchar(x, type = "width") > width) paste0(substr(x, 1, width), "...") else x
}

bad_found <- FALSE
total_bad <- 0L

for (f in files) {
  txt <- readLines(f, warn = FALSE, encoding = "UTF-8")
  bad_idx <- grep("[^\\p{ASCII}]", txt, perl = TRUE)

  if (length(bad_idx)) {
    bad_found <- TRUE
    total_bad <- total_bad + length(bad_idx)
    cat(sprintf("\n%s\n", f))
    for (j in bad_idx) {
      cat(sprintf("  line %d: %s\n", j, trim_preview(txt[j])))
    }
  }
}

if (!bad_found) {
  cat("OK: no non-ASCII characters found in scanned .R files.\n")
  quit(status = 0)
}

cat(sprintf("\nTotal problematic lines: %d\n", total_bad))
#quit(status = 1)
