# Runs the test suite from the repository root: Rscript tests/run.R [filter]
for (path in list.files("R", pattern = "[.]R$", full.names = TRUE)) sys.source(path, envir = globalenv())
args <- commandArgs(trailingOnly = TRUE)
testthat::test_dir("tests/testthat", filter = if (length(args)) args[1L] else NULL,
  env = new.env(parent = globalenv()), stop_on_failure = TRUE)
