## ============================================================
## build_package.R
##
## Developer script for building, documenting, and checking the
## assetIS package locally.  Run from the PARENT directory of
## the assetIS/ folder, e.g.:
##
##   Rscript build_package.R
##
## or interactively from within the package root via devtools.
## ============================================================
rm(list = ls())
WD <- "~/Desktop/RpackagesDev/isASSET"
setwd(WD)

## ---- 1. Install build-time dependencies (run once) ----------
# install.packages(c("devtools", "roxygen2", "pkgdown",
#                    "testthat", "covr", "MASS", "knitr",
#                    "rmarkdown", "ggplot2", "tidyr"))

## ---- 2. Load devtools ----------------------------------------
if (!requireNamespace("devtools", quietly = TRUE))
  stop("Please install devtools: install.packages('devtools')")

pkg <- "assetIS"   # adjust path if running from elsewhere

## ---- 3. Generate roxygen2 documentation ----------------------
message("=== Generating documentation with roxygen2 ===")
devtools::document(pkg = pkg)

## ---- 4. Run tests --------------------------------------------
message("=== Running testthat tests ===")
devtools::test(pkg = pkg)

## ---- 5. R CMD check (no vignettes for speed) -----------------
message("=== Running R CMD check ===")
devtools::check(
  pkg            = pkg,
  document       = FALSE,   # already done above
  vignettes      = FALSE,   # build separately with build_vignettes()
  manual         = FALSE,
  args           = c("--as-cran")
)

## ---- 6. Build and check with vignettes (CRAN-style) ----------
## Uncomment for a full pre-submission check:
# devtools::check(pkg = pkg, vignettes = TRUE, manual = TRUE,
#                 args = "--as-cran")


## step  - Builds the package
# This step will generate a .tar.gz package (a zip file)
setwd(WD)
devtools::build(pkg = pkg) ## default argument is pkg = ".", current working directory


## ---- 7. Install locally --------------------------------------
message("=== Installing package locally ===")
devtools::install(pkg = pkg, upgrade = "never")

## ---- 8. Build pkgdown site (optional) -----------------------
## Uncomment to build the documentation website:
# pkgdown::build_site(pkg = pkg)

message("=== Done ===")
