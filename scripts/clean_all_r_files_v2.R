args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args) >= 1) args[[1]] else "."
backup <- if (length(args) >= 2) as.logical(args[[2]]) else TRUE

clean_all_r_files <- function(root = ".", backup = TRUE) {
  files <- c(
    list.files(file.path(root, "R"), pattern = "\\.[Rr]$", recursive = TRUE, full.names = TRUE),
    #list.files(file.path(root, "tests", "testthat"), pattern = "\\.[Rr]$", recursive = TRUE, full.names = TRUE),
    list.files(file.path(root, "data-raw"), pattern = "\\.[Rr]$", recursive = TRUE, full.names = TRUE)
  )
  files <- unique(files[file.exists(files)])

  if (!length(files)) {
    message("No .R files found in R/, tests/testthat/, or data-raw/.")
    return(invisible(NULL))
  }

  replacements <- c(
    "\u2014" = "-",       # em dash —
    "\u2013" = "-",       # en dash –
    "\u2212" = "-",       # minus sign −
    "\u2260" = "!=",      # ≠
    "\u2264" = "<=",      # ≤
    "\u2265" = ">=",      # ≥
    "\u2192" = "->",      # →
    "\u2190" = "<-",      # ←
    "\u00D7" = "x",       # ×
    "\u039B" = "Lambda",  # Λ
    "\u03BB" = "lambda",  # λ
    "\u0160" = "S",       # Š
    "\u0161" = "s",       # š
    "\u00E1" = "a",       # á
    "\u00E9" = "e",       # é
    "\u00ED" = "i",       # í
    "\u00F3" = "o",       # ó
    "\u00FA" = "u",       # ú
    "\u2026" = "...",     # …
    "\u2500" = "-",       # ─
    "\u2550" = "=",       # ═
    "\u00D7" = "x"        # ×
  )

  changed_files <- 0L

  for (f in files) {
    txt <- readLines(f, warn = FALSE, encoding = "UTF-8")
    original <- txt

    for (sym in names(replacements)) {
      txt <- gsub(sym, replacements[[sym]], txt, fixed = TRUE)
    }

    if (!identical(txt, original)) {
      if (backup) {
        file.copy(f, paste0(f, ".bak"), overwrite = TRUE)
      }
      writeLines(txt, f, useBytes = TRUE)
      changed_files <- changed_files + 1L
      message("Cleaned: ", f)
    }
  }

  message("\nTotal files modified: ", changed_files)
  invisible(changed_files)
}

clean_all_r_files(root = root, backup = backup)

