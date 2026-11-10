# Scores of predictive draws against their definitions.

score_case <- function(S = 4000L) {
  set.seed(11)
  held <- data.frame(observation_id = c("a", "b", "c"), pollutant = "PM25",
    physical_location_id = c("l1", "l1", "l2"), month = c("2024-01", "2024-02", "2024-01"),
    observed = c(10, 14, 9), source = "ground_pm25", device_id = "", site_class = "avenue",
    stringsAsFactors = FALSE)
  mu <- log(c(11, 12, 9.5))
  z <- matrix(stats::rnorm(3 * S), 3)
  draws <- exp(mu + 0.2 * z)
  log_density <- matrix(stats::dnorm(log(held$observed), mu, 0.2, log = TRUE) - log(held$observed), 3, S)
  rownames(draws) <- rownames(log_density) <- held$observation_id
  list(held = held, predicted = list(draws = draws, log_density = log_density))
}

testthat::test_that("each score follows its definition", {
  x <- score_case()
  s <- score_targets(x$held, x$predicted, "final", "new-traffic-1")
  draws <- x$predicted$draws
  testthat::expect_identical(names(s), score_columns)
  testthat::expect_identical(s$observation_id, x$held$observation_id)
  testthat::expect_equal(s$predictive_mean, unname(rowMeans(draws)), tolerance = 1e-12)
  testthat::expect_equal(s$predictive_sd, unname(sqrt(rowMeans((draws - rowMeans(draws))^2))),
    tolerance = 1e-12)
  testthat::expect_equal(s$bias, s$predictive_mean - x$held$observed)
  # CRPS of an empirical law: E|X - y| - E|X - X'| / 2.
  x400 <- score_case(400L)
  s400 <- score_targets(x400$held, x400$predicted, "m", "u")
  crps400 <- vapply(1:3, function(i) mean(abs(x400$predicted$draws[i, ] - x400$held$observed[i])) -
    sum(abs(outer(x400$predicted$draws[i, ], x400$predicted$draws[i, ], "-"))) / (2 * 400^2), numeric(1))
  testthat::expect_equal(s400$crps, crps400, tolerance = 1e-10)
  # NLS: minus the log of the mean density; here every draw has the same density.
  testthat::expect_equal(s$nls, -unname(x$predicted$log_density[, 1]), tolerance = 1e-12)
  testthat::expect_equal(s$pit, unname(rowMeans(draws < x$held$observed)), tolerance = 1e-12)
  # Central intervals: type-1 quantiles of the draws.
  for (level in c(50, 80, 95)) {
    a <- (1 - level / 100) / 2
    q <- unname(t(apply(draws, 1L, stats::quantile, probs = c(a, 1 - a), type = 1, names = FALSE)))
    testthat::expect_equal(s[[paste0("lower_", level)]], q[, 1])
    testthat::expect_equal(s[[paste0("upper_", level)]], q[, 2])
    testthat::expect_identical(s[[paste0("coverage_", level)]],
      x$held$observed >= q[, 1] & x$held$observed <= q[, 2])
  }
  # Monte Carlo SE of the CRPS: SD over 8 batches of 500 draws.
  batch <- vapply(1:8, function(b) scoringRules::crps_sample(x$held$observed,
    draws[, (b - 1) * 500 + 1:500], method = "edf"), numeric(3))
  testthat::expect_equal(s$crps_mcse, apply(batch, 1L, stats::sd) / sqrt(8), tolerance = 1e-10)
  testthat::expect_error(score_targets(x$held[3:1, ], x$predicted, "m", "u"), "held order")
})

testthat::test_that("summaries weight locations equally", {
  x <- score_case()
  s <- score_targets(x$held, x$predicted, "final", "new-traffic-1")
  sm <- summarise_locations(s)
  testthat::expect_identical(c(sm$n_targets, sm$n_locations), c(3L, 2L))
  testthat::expect_equal(sm$crps, mean(c(mean(s$crps[1:2]), s$crps[3])))
  testthat::expect_equal(sm$rmse, sqrt(mean(c(mean(s$squared_error[1:2]), s$squared_error[3]))))
  pooled <- rbind(s, transform(s, unit = "new-traffic-2", observation_id = paste0(s$observation_id, "2")),
    transform(s, unit = "new-traffic-3", observation_id = paste0(s$observation_id, "3")))
  groups <- summarise_units(pooled)
  testthat::expect_identical(groups$group, c("new-traffic", "new-traffic-1", "new-traffic-2", "new-traffic-3"))
  testthat::expect_equal(groups$crps[1], sm$crps)
  testthat::expect_error(summarise_locations(rbind(s, transform(s, pollutant = "NO2"))), "one pollutant")
})
