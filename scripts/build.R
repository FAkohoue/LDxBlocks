
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
# BiocManager::install("SeqArray")
# install.packages("SeqArray", force = TRUE)

# 0. Recreate NAMESPACE and man/ from roxygen tags

remove.packages("LDxBlocks")
.rs.restartR()

# Recreate NAMESPACE and man/ from roxygen tags
source("data-raw/generate_example_data.R")

# Step 2: now document will work
devtools::document()


# 4. Now compileAttributes will work (DESCRIPTION exists, src/ exists)
Rcpp::compileAttributes()

# 5. Document again to pick up RcppExports.R
devtools::document()


# Step 1: Detach and unload the package completely
try(detach("package:LDxBlocks", unload = TRUE, force = TRUE), silent = TRUE)

# Step 2: Restart R session to release the file lock
.rs.restartR()

# 6. Install
devtools::install()

# 7. Test
devtools::test()


# 8. Check
devtools::check()


# 9. Build vignettes
#options(pkgdown.internet = FALSE)
library(LDxBlocks)
# Build everything except home, then build home separately
#pkgdown::build_favicons(overwrite = TRUE)
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
