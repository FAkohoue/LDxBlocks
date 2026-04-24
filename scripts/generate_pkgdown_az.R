# scripts/generate_pkgdown_az.R

namespace_file <- "NAMESPACE"

if (!file.exists(namespace_file)) {
  stop("NAMESPACE not found. Run devtools::document() first.", call. = FALSE)
}

ns <- readLines(namespace_file, warn = FALSE)

exports <- grep("^export\\(", ns, value = TRUE)
exports <- sub("^export\\(", "", exports)
exports <- sub("\\)$", "", exports)
exports <- sort(unique(exports))

# Optional: add datasets manually if you want them in the A-Z index
datasets <- c(
  "ldx_blocks",
  "ldx_blues",
  "ldx_blues_list",
  "ldx_geno",
  "ldx_gwas",
  "ldx_snp_info"
)

items <- sort(unique(c(exports, datasets)))

# Group by first letter
first_letter <- toupper(substr(items, 1L, 1L))
groups <- split(items, first_letter)

# Keep only A-Z letters, in order
letters_az <- LETTERS[LETTERS %in% names(groups)]

cat("# ------------------------------------------------------------\n")
cat("# A-Z function index, ASReml style\n")
cat("# ------------------------------------------------------------\n\n")

for (L in letters_az) {
  cat("  - title: \"", L, "\"\n", sep = "")
  cat("    desc: >\n")
  cat("      Functions and datasets beginning with ", L, ".\n", sep = "")
  cat("    contents:\n")

  for (x in groups[[L]]) {
    cat("      - ", x, "\n", sep = "")
  }

  cat("\n")
}
