# This file is the entry point for devtools::test() and R CMD check.
# It ensures the compiled C++ shared library is loaded before any test runs.

library(testthat)
library(LDxBlocks)

test_check("LDxBlocks")
