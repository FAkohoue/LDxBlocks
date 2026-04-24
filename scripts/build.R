
install.packages("rsvg")
rsvg::rsvg_png(
  "man/figures/logo.svg",
  "man/figures/logo.png",
  width  = 1360,
  height = 560
)

# Then build favicons from the PNG instead
pkgdown::build_favicons(overwrite = TRUE)

options(timeout = 3000)  # 5 minutes
pkgdown::build_favicons(overwrite = TRUE)

# Regenerate favicons from the correct logo
pkgdown::build_favicons(overwrite = TRUE)
pkgdown::build_site()
pkgdown::build_site(override = list(template = list(favicon = FALSE)))

library(magick)

# Load your logo
img <- image_read("man/figures/logo.png")

# Create favicon sizes
sizes <- c(16, 32, 48, 64, 180, 192, 512)

dir.create("pkgdown/assets", recursive = TRUE, showWarnings = FALSE)

for (s in sizes) {
  resized <- image_resize(img, paste0(s, "x", s))
  image_write(resized, paste0("pkgdown/assets/favicon-", s, ".png"))
}

# Create favicon.ico (multi-size)
ico <- image_resize(img, "64x64")
image_write(ico, "pkgdown/assets/favicon.ico")

writeLines(
  c(
    '<link rel="icon" type="image/svg+xml" href="logo.svg">',
    '<link rel="icon" type="image/png" sizes="32x32" href="favicon-32.png">',
    '<link rel="apple-touch-icon" href="favicon-180.png">'
  ),
  "pkgdown/assets/favicon.html"
)


# Create the correct folder
dir.create("pkgdown/favicon", showWarnings = FALSE)

# Copy all favicon files from assets/ to favicon/
file.copy(
  from      = list.files("pkgdown/assets", full.names = TRUE),
  to        = "pkgdown/favicon/",
  overwrite = TRUE
)

list.files("pkgdown/favicon")





######################################################################################################
#######################################################################################################

# =============================================================================
# LDxBlocks - clean build script
# Run SESSION A top-to-bottom, let it restart R, then run SESSION B.
# Both sessions must start with the working directory set to the package root
# (the folder containing DESCRIPTION).
# =============================================================================


# ------------------------------- SESSION A -----------------------------------
# Goal: wipe every stale artifact so the next session starts from zero.

# 1. Remove the installed package so no old DLL can be loaded by mistake.
remove.packages("LDxBlocks")

unlink(file.path(.libPaths()[1], "00LOCK-LDxBlocks"), recursive = TRUE)

# 2. Delete the compiled DLL from the source tree (src/*.so / src/*.dll).
#    Do this BEFORE the restart so there is nothing to unload conflicts with.
#devtools::clean_dll()

# 3. Restart R.  All loaded DLLs are released, file locks are cleared.
.rs.restartR()
# -- after restart, continue in SESSION B --------------------------------------




# ------------------------------- SESSION B -----------------------------------
# Goal: rebuild everything from source in the correct dependency order.

# 1. Regenerate example data (.rda files in data/ and flat files in
#    inst/extdata/).  Must come first because devtools::document() will
#    try to lazy-load data/ when it parses roxygen @examples.
source("data-raw/generate_example_data.R")

# 2. First document pass: generates NAMESPACE and man/ from roxygen tags in
#    R/*.R.  RcppExports.R does not exist yet so its tags are skipped.
devtools::document()

# 3. Regenerate RcppExports.R (in R/) and src/RcppExports.cpp from the
#    [[Rcpp::export]] annotations in src/ld_core.cpp.
#    This is the canonical source of truth for the C++ symbol table.
Rcpp::compileAttributes()

# 4. Second document pass: picks up the roxygen tags now present in
#    the freshly written R/RcppExports.R.
devtools::document()

# 5. Compile C++ and install into the library.
#    upgrade = FALSE  - do not touch other packages.
devtools::install()


Sys.setenv(PATH = paste(
  "C:/Program Files/TASSEL5/jre/bin",
  Sys.getenv("PATH"),
  sep = ";"
))

# 6. Run the test suite.  All C++ symbols are now registered in the
#    installed DLL, so load_all() will find them.
devtools::test()

#devtools::test(filter = "association")

# 7. Full CRAN check (run after tests pass).
devtools::check()


# 9. Build vignettes
#options(pkgdown.internet = FALSE)
library(LDxBlocks)

# Build everything except home, then build home separately
#pkgdown::build_favicons(overwrite = TRUE)

#pkgdown::check_pkgdown()


pkgdown::build_reference()
pkgdown::build_articles()
pkgdown::build_news()

# Build home with network disabled at the curl level
httr2_mock <- function(...) stop("no network", call. = FALSE)
pkgdown::build_home()

#unloadNamespace("OptSLDP")

#pkgdown::clean_site(force = TRUE)
pkgdown::build_site()

# 10. Build package
devtools::build()






# Delete the cached rendered output for this vignette
unlink("docs/articles/LDxBlocks-large-scale.html")
unlink("docs/articles/LDxBlocks-large-scale_files", recursive = TRUE)

# Also clear any pkgdown article cache
unlink("vignettes/cache", recursive = TRUE)
unlink("vignettes/LDxBlocks-large-scale_cache", recursive = TRUE)

# Now rebuild just that article from scratch
pkgdown::build_article("LDxBlocks-large-scale")


file <- "R/haplotype_association.R"

tools::showNonASCIIfile(file)

txt <- readLines(file, encoding = "UTF-8")

txt <- gsub("->", "->", txt, fixed = TRUE)
txt <- gsub("x", "x", txt, fixed = TRUE)
txt <- gsub("-", "-", txt, fixed = TRUE)
txt <- gsub("...", "...", txt, fixed = TRUE)

writeLines(txt, file, useBytes = TRUE)

tools::showNonASCIIfile(file)


txt <- readLines(file, encoding = "UTF-8")

txt <- gsub("-", "-", txt, fixed = TRUE)
txt <- gsub("-", "-", txt, fixed = TRUE)

writeLines(txt, file, useBytes = TRUE)

txt <- readLines(file, encoding = "UTF-8")

txt <- gsub("->", "->", txt, fixed = TRUE)
txt <- gsub("x", "x", txt, fixed = TRUE)
txt <- gsub("-", "-", txt, fixed = TRUE)
txt <- gsub("...", "...", txt, fixed = TRUE)
txt <- gsub("union", "union", txt, fixed = TRUE)
txt <- gsub("<=", "<=", txt, fixed = TRUE)
txt <- gsub(">=", ">=", txt, fixed = TRUE)
writeLines(txt, file, useBytes = TRUE)

tools::showNonASCIIfile(file)




txt <- readLines(file, encoding = "UTF-8")

replacements <- c(
  "->"="->",
  "x"="x",
  "-"="-",
  "-"="-",
  "-"="-",
  "..."="...",
  "union"="union",
  "<="="<=",
  ">="=">=",
  "^3"="^3",
  "Lambda"="Lambda"
)

for (sym in names(replacements)) {
  txt <- gsub(sym, replacements[[sym]], txt, fixed = TRUE)
}

writeLines(txt, file, useBytes = TRUE)

tools::showNonASCIIfile(file)
