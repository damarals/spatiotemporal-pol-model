# Paired contrasts, location-by-year block bootstrap, and result tables.

contrast_fixture <- function() {
  base <- data.frame(
    model = "with",
    unit = c("new-traffic-1", "new-traffic-1", "new-traffic-2", "new-traffic-3"),
    observation_id = c("o1", "o2", "o3", "o4"),
    pollutant = "NO2",
    physical_location_id = c("s1", "s1", "s2", "s3"),
    month = c("2022-01", "2022-02", "2022-06", "2023-03"),
    observed = c(15, 18, 12, 20),
    source = "ground_no2",
    device_id = "",
    site_class = c("avenue", "avenue", "bus_terminal", "avenue"),
    predictive_mean = c(14, 17, 13, 19),
    predictive_sd = c(2, 2, 2, 2),
    bias = c(-1, -1, 1, -1),
    squared_error = c(1, 1, 1, 1),
    absolute_error = c(1, 1, 1, 1),
    crps = c(1.0, 1.4, 0.8, 1.2),
    nls = c(1.8, 2.0, 1.6, 1.9),
    pit = c(0.4, 0.5, 0.6, 0.4),
    lower_50 = c(12, 15, 11, 17), upper_50 = c(16, 19, 15, 21),
    coverage_50 = c(TRUE, TRUE, TRUE, TRUE), width_50 = c(4, 4, 4, 4),
    lower_80 = c(11, 14, 10, 16), upper_80 = c(17, 20, 16, 22),
    coverage_80 = c(TRUE, TRUE, TRUE, TRUE), width_80 = c(6, 6, 6, 6),
    lower_95 = c(10, 13, 9, 15), upper_95 = c(18, 21, 17, 23),
    coverage_95 = c(TRUE, TRUE, TRUE, TRUE), width_95 = c(8, 8, 8, 8),
    crps_mcse = c(0.03, 0.03, 0.04, 0.05),
    stringsAsFactors = FALSE)
  without <- transform(base, model = "without",
    crps = c(1.2, 1.8, 1.1, 1.5), crps_mcse = c(0.04, 0.04, 0.03, 0.05))
  list(with = base[score_columns], without = without[score_columns])
}

testthat::test_that("paired_contrast weights locations equally and computes paired MCSE", {
  f <- contrast_fixture()
  c1 <- paired_contrast(f$with, f$without, "new-traffic", "NO2", seed = 101L, B = 500L)
  # Location means: s1 = 1.2 vs 1.5 (d = 0.3), s2 = 0.8 vs 1.1 (d = 0.3), s3 = 1.2 vs 1.5 (d = 0.3).
  testthat::expect_identical(c(c1$n, c1$locations), c(4L, 3L))
  testthat::expect_equal(c1$crps_with, mean(c(1.2, 0.8, 1.2)))
  testthat::expect_equal(c1$crps_without, mean(c(1.5, 1.1, 1.5)))
  testthat::expect_equal(c1$gain, 0.3)
  # Paired MCSE per target: sqrt(0.03^2 + 0.04^2) = 0.05 on o1, o2, o3; sqrt(2)*0.05 on o4.
  w <- c(1 / 6, 1 / 6, 1 / 3, 1 / 3)
  m <- c(0.05, 0.05, 0.05, sqrt(2) * 0.05)
  testthat::expect_equal(c1$mcse, sum(w * m))
  testthat::expect_equal(c1$mcse_indep, sqrt(sum(w^2 * m^2)))
  testthat::expect_identical(c1$locations_agreeing, 3L)
  testthat::expect_true(c1$stable)
  lg <- location_gains(f$with, f$without, "new-traffic", "NO2")
  testthat::expect_equal(as.numeric(lg[c("s1", "s2", "s3")]), c(0.3, 0.3, 0.3))
})

testthat::test_that("pooled_contrast excludes withheld-background NO2 and resamples location-year blocks", {
  f <- contrast_fixture()
  bg <- transform(f$with[1L, ], unit = "withheld-background", observation_id = "bg1",
    physical_location_id = background_location, crps = 3.0)
  bg_wo <- transform(bg, model = "without", crps = 9.0)
  p <- pooled_contrast(rbind(f$with, bg), rbind(f$without, bg_wo), seed = 99L, B = 200L)
  testthat::expect_identical(c(p$n_targets, p$n_blocks), c(4L, 3L))
  testthat::expect_equal(p$gain, mean(f$without$crps - f$with$crps))
  testthat::expect_equal(p$p_with_better, 1)
})

testthat::test_that("pooled_contrast does not depend on the order of the input rows", {
  f <- contrast_fixture()
  f$without$crps <- c(1.2, 1.9, 1.0, 1.6)
  p <- pooled_contrast(f$with, f$without, seed = 7L, B = 300L)
  q <- pooled_contrast(f$with[c(4L, 2L, 3L, 1L), ], f$without[c(3L, 1L, 4L, 2L), ], seed = 7L, B = 300L)
  testthat::expect_identical(p, q)
  testthat::expect_true(p$boot95_lo < p$boot95_hi)
})

testthat::test_that("the location gains table gives the gain and count of each traffic location", {
  f <- contrast_fixture()
  steps <- data.frame(reference = "r", component = "c", change = "removed", model_with = "with",
    model_without = "without", stringsAsFactors = FALSE)
  g <- location_gains_table(rbind(f$with, f$without), steps)
  testthat::expect_identical(g$location_id, c("s1", "s2", "s3"))
  testthat::expect_identical(g$n_targets, c(2L, 1L, 1L))
  testthat::expect_equal(g$gain, c(0.3, 0.3, 0.3))
  testthat::expect_identical(g$site_class, c("avenue", "bus_terminal", "avenue"))
  testthat::expect_identical(unique(g$pollutant), "NO2")
})

testthat::test_that("the full model against the final model is the elimination comparison numbered after its steps", {
  f <- contrast_fixture()
  full <- transform(f$with, model = "full")
  final <- transform(f$without, model = "final")
  x <- full_vs_final_table(rbind(full, final))
  testthat::expect_identical(c(x$model_with, x$model_without), c("full", "final"))
  testthat::expect_identical(x$reference, "full model")
  testthat::expect_identical(x$bootstrap_seed, 20261109L + 100L * 27L + 1L)
  y <- elimination_table(rbind(full, final), full_vs_final_step)
  testthat::expect_identical(x[setdiff(names(x), "bootstrap_seed")], y[setdiff(names(y), "bootstrap_seed")])
  p <- paired_contrast(full, final, "new-traffic", "NO2", seed = x$bootstrap_seed)
  testthat::expect_equal(c(x$gain, x$boot95_lo, x$boot95_hi), round(c(p$gain, p$boot95_lo, p$boot95_hi), 3))
})

testthat::test_that("results/full-vs-final.csv compares the two models on every evaluation group", {
  elim <- utils::read.csv(file.path(repo_root, "results", "elimination.csv"), stringsAsFactors = FALSE)
  x <- utils::read.csv(file.path(repo_root, "results", "full-vs-final.csv"), stringsAsFactors = FALSE)
  testthat::expect_identical(names(x), names(elim))
  testthat::expect_identical(paste(x$group, x$pollutant),
    paste(evaluation_groups$group, evaluation_groups$pollutant))
  testthat::expect_true(all(x$model_with == "full" & x$model_without == "final"))
})

testthat::test_that("elimination_steps and paper table selectors match the repository models", {
  testthat::expect_identical(nrow(elimination_steps), 26L)
  testthat::expect_true(all(elimination_steps$model_with %in% models$model))
  testthat::expect_true(all(elimination_steps$model_without %in% models$model))
  if (file.exists(file.path(repo_root, "results", "elimination.csv"))) {
    elim <- utils::read.csv(file.path(repo_root, "results", "elimination.csv"), stringsAsFactors = FALSE)
    testthat::expect_identical(nrow(elim), 186L)
    t5 <- paper_contrasts_table(elim, "structure")
    t6 <- paper_contrasts_table(elim, "sources")
    testthat::expect_identical(nrow(t5), 14L)
    testthat::expect_identical(nrow(t6), 30L)
  }
})