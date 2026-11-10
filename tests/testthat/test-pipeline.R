test_that("the pipeline fits each model on the units of the models table", {
  env <- new.env()
  withr::with_dir(repo_root, sys.source("_targets.R", envir = env))
  expected <- do.call(rbind, lapply(models$model, function(m)
    data.frame(model = m, unit = setdiff(model_units(m), "all-data"), stringsAsFactors = FALSE)))
  expect_identical(env$validation_fits, expected)
  expect_true(all(validation_units$unit[validation_units$unit != "all-data"] %in% env$validation_fits$unit))
})
