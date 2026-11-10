# Inputs and the partitions of the validation units.

testthat::test_that("the inputs hold every target with its observation id", {
  skip_without_moqa()
  x <- test_inputs()
  t <- x$targets
  testthat::expect_identical(as.vector(table(t$source)[c("ground_no2", "ground_pm25", "moqa")]),
    c(49L, 43L, 21L))
  testthat::expect_false(anyDuplicated(t$observation_id) > 0L)
  testthat::expect_true(all(grepl("^seuma:P0[0-8]:[0-9]{4}-[0-9]{2}$", t$observation_id[t$source == "ground_no2"])))
  testthat::expect_true(all(grepl("^seuma:fortaleza-.*:[0-9]{4}-[0-9]{2}$", t$observation_id[t$source == "ground_pm25"])))
  testthat::expect_true(all(grepl("^moqa:[A-Z0-9]{6}:[0-9]{4}-[0-9]{2}$", t$observation_id[t$source == "moqa"])))
  # The calibration-overlap months of the background location are not targets.
  testthat::expect_false(any(t$source == "ground_pm25" & t$physical_location_id == background_location &
    t$month %in% c("2024-09", "2024-10", "2024-12", "2025-01")))
  testthat::expect_true(all(t$month >= "2020-01" & t$month <= "2025-12"))
  testthat::expect_true(all(t$observed > 0))
  testthat::expect_identical(c(nrow(x$acag), nrow(x$acag_cells), nrow(x$column), nrow(x$calibration)),
    c(60L, 330L, 1584L, 66L))
  testthat::expect_equal(sum(x$acag_cells$weight), 1, tolerance = 1e-12)
  testthat::expect_identical(x$era5$month, .months("2020-01", "2025-12"))
})

testthat::test_that("the public site tables hold only the reference locations", {
  dir <- file.path(repo_root, "data")
  sites <- utils::read.csv(file.path(dir, "processed", "sites.csv"), stringsAsFactors = FALSE)
  covariates <- utils::read.csv(file.path(dir, "processed", "site-covariates.csv"), stringsAsFactors = FALSE)
  testthat::expect_true(all(sites$source == "seuma"))
  testthat::expect_identical(covariates$location_id, sites$location_id)
  testthat::expect_identical(sort(sites$location_id),
    sort(utils::read.csv(file.path(dir, "decisions", "seuma-site-classes.csv"))$location_id))
  testthat::expect_false(any(grepl("^moqa", list.files(file.path(dir, "decisions")))))
})

testthat::test_that("data_files() lists what the readers open, with or without names", {
  files <- data_files("d")
  testthat::expect_identical(unname(files[c("sites", "moqa_sites")]), c("d/processed/sites.csv", "d/moqa/sites.csv"))
  testthat::expect_identical(.data_paths(unname(files)), files)
  testthat::expect_identical(data_file(unname(files), "boundary"), "d/processed/boundary.gpkg")
  testthat::expect_error(.data_paths(files[-1L]), "data_files")
  testthat::expect_error(read_inputs(file.path(tempdir(), "no-data")), "low-cost network")
  skip_without_moqa()
  testthat::expect_true(all(file.exists(data_files(file.path(repo_root, "data")))))
  x <- test_inputs()
  testthat::expect_identical(x$sites$source, rep(c("moqa", "seuma"), c(14L, 8L)))
  testthat::expect_identical(x$site_covariates$location_id, x$sites$location_id)
})

testthat::test_that("each unit withholds what the paper describes and nothing else enters", {
  skip_without_moqa()
  x <- test_inputs()
  counts <- c("omitted-months" = 8L, "new-traffic-1" = 10L, "new-traffic-2" = 10L,
    "new-traffic-3" = 16L, "withheld-background" = 56L, "year-2023" = 8L, "year-2024" = 14L,
    "year-2025" = 22L, "new-devices-1" = 8L, "new-devices-2" = 3L, "new-devices-3" = 10L,
    "all-data" = 0L)
  for (unit in validation_units$unit) {
    p <- partition(x, unit)
    testthat::expect_identical(nrow(p$held), counts[[unit]], info = unit)
    used <- c(p$train$observation_id, p$moqa$observation_id)
    testthat::expect_false(any(p$held$observation_id %in% used), info = unit)
    testthat::expect_true(all(p$train$source %in% c("ground_no2", "ground_pm25")), info = unit)
    # No held location's low-cost rows, pairs or held-month BAM days enter the fit.
    held_ref <- p$held[p$held$source != "moqa", , drop = FALSE]
    testthat::expect_false(any(p$moqa$physical_location_id %in% held_ref$physical_location_id), info = unit)
    testthat::expect_false(any(paste(p$calibration$location_id, p$calibration$month) %in%
      paste(held_ref$physical_location_id[held_ref$source == "ground_pm25"],
        held_ref$month[held_ref$source == "ground_pm25"])), info = unit)
    if (!is.null(p$origin)) {
      testthat::expect_true(all(c(p$train$month, p$moqa$month, p$calibration$month, p$acag$month,
        p$column$month) < p$origin), info = unit)
      testthat::expect_true(all(substr(p$held$month, 1L, 4L) == substr(p$origin, 1L, 4L)), info = unit)
    }
  }
  # The background location carries the calibration pairs: withholding it removes the network.
  testthat::expect_identical(nrow(partition(x, "withheld-background")$moqa), 0L)
  testthat::expect_identical(partition(x, "year-2025")$origin, "2025-01")
  testthat::expect_error(partition(x, "F1"), "Unknown unit")
})

testthat::test_that("every model names known units and the flags of its family", {
  testthat::expect_false(anyDuplicated(models$model) > 0L)
  for (model in models$model) testthat::expect_true(length(model_units(model)) > 0L, info = model)
  testthat::expect_identical(model_units("final"), validation_units$unit)
  final <- model_spec("final")
  testthat::expect_identical(unlist(final[c("tropomi", "temporal_loading", "sign_constraint",
    "pm25_profile", "deviation", "amplitude", "profile_cor", "free_gain", "no2_cov", "pm25_cov")]),
    c(tropomi = 1L, temporal_loading = 0L, sign_constraint = 1L, pm25_profile = 0L, deviation = 0L,
      amplitude = 0L, profile_cor = 0L, free_gain = 1L, no2_cov = 0L, pm25_cov = 1L))
  # Each round-1 model differs from its reference in one declared column.
  flags <- setdiff(names(models), c("model", "units"))
  base <- models[models$model == "full-no-tropomi", flags]
  for (model in grep("^full-no-tropomi-", models$model, value = TRUE)) {
    changed <- flags[unlist(models[models$model == model, flags]) != unlist(base)]
    testthat::expect_length(changed, 1L)
  }
})
