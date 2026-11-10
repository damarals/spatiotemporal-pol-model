# Predictive checks: TROPOMI statistics by hand and against the dense
# month-centering covariance, invariance to common monthly shifts, moments under
# the column law, the tail summary, and the 55 marginal features.

.test_col_rows <- function(cells, months) {
  grid <- expand.grid(cell_id = cells, month = months, stringsAsFactors = FALSE)
  grid$observation_id <- paste0("tropomi:", grid$month, ":", grid$cell_id)
  grid[, c("observation_id", "cell_id", "month")]
}

.test_dense_z <- function(r, d, group, kappa) {
  z <- numeric(length(r))
  for (g in unique(group)) {
    i <- which(group == g)
    C <- diag(length(i)) - matrix(1 / length(i), length(i), length(i))
    covariance <- C %*% diag(d[i]^2 / kappa, length(i)) %*% t(C)
    z[i] <- drop(C %*% r[i]) / sqrt(diag(covariance))
  }
  z
}

testthat::test_that("column statistics match hand and dense calculations", {
  rows <- .test_col_rows(c("1:1", "1:2", "1:3"), c("2024-01", "2024-02", "2024-03"))
  panel <- .column_panel(rows)
  r <- c(1, 2, 6, 0, 3, 3, 2, -1, 5)
  d <- c(1, 1, 2, 2, 1, 1, 1, 2, 1)
  kappa <- 4
  result <- .column_values(r, rep(0, 9), rep(2, 9), d, kappa, panel)
  testthat::expect_named(result, column_statistics)
  testthat::expect_equal(unname(result["dispersion"]), 89 / 36)
  testthat::expect_equal(unname(result["maximum"]), 3)
  testthat::expect_equal(unname(result["monthly"]), 5)
  z <- .test_dense_z(r, d, rep(1:3, each = 3), kappa)
  testthat::expect_equal(unname(result["contrast_dispersion"]), mean(z^2), tolerance = 1e-12)
  lag <- cbind(c(1, 2, 3, 4, 5, 6), c(4, 5, 6, 7, 8, 9))
  testthat::expect_equal(unname(result["within_cell"]), mean(z[lag[, 1]] * z[lag[, 2]]),
    tolerance = 1e-12)
  neighbor <- cbind(c(1, 2, 4, 5, 7, 8), c(2, 3, 5, 6, 8, 9))
  testthat::expect_equal(unname(result["neighbor"]), mean(z[neighbor[, 1]] * z[neighbor[, 2]]),
    tolerance = 1e-12)
  testthat::expect_equal(z[1:3], c(-2, -1, 3) / sqrt(c(1, 1, 2) / 4), tolerance = 1e-12)
})

testthat::test_that("contrasts ignore what the cells of a month share and keep missing cells out", {
  rows <- .test_col_rows(c("1:1", "1:2", "2:1", "2:2"), c("2024-01", "2024-02", "2024-04"))
  rows <- rows[-6, ]
  panel <- .column_panel(rows)
  testthat::expect_equal(nrow(panel$month_pairs), 1L)
  testthat::expect_equal(nrow(panel$lag_pairs), 3L)
  testthat::expect_equal(nrow(panel$neighbor_pairs), 4L + 2L + 4L)
  r <- c(0.3, -1.2, 0.8, 2.0, -0.4, 1.1, 0.5, -0.7, 0.1, 1.6, 0.9)
  d <- c(0.7, 0.9, 1.1, 0.8, 1.0, 0.6, 0.9, 0.8, 1.2, 0.7, 1.0)
  sd <- sqrt(0.04 + d^2 / 9)
  one <- .column_values(r, rep(0, 11), sd, d, 9, panel)
  shifted <- .column_values(r + c(5, -3, 0.5)[panel$group], rep(0, 11), sd, d, 9, panel)
  contrasts <- c("contrast_dispersion", "within_cell", "neighbor")
  testthat::expect_equal(shifted[contrasts], one[contrasts], tolerance = 1e-12)
  testthat::expect_false(isTRUE(all.equal(shifted["dispersion"], one["dispersion"])))
  z <- .test_dense_z(r, d, panel$group, 9)
  testthat::expect_equal(unname(one["contrast_dispersion"]), mean(z^2), tolerance = 1e-12)
  testthat::expect_error(.column_panel(rows[c(1, 1:11), ]), "unique")
  testthat::expect_error(.column_panel(rows[c(1:4, 9), ]), "one cell")
})

testthat::test_that("statistics have their law moments under the column law", {
  rows <- .test_col_rows(c("1:1", "1:2", "2:1", "2:2"), sprintf("2024-%02d", 1:8))
  panel <- .column_panel(rows)
  d <- rep(c(0.6, 0.9, 1.2, 0.75), 8) * rep(c(1, 1.1, 0.9, 1, 1.2, 0.8, 1, 1.05), each = 4)
  kappa <- 6
  rho <- 0.5
  sigma_q <- 0.7
  sd <- sqrt(sigma_q^2 + d^2 / kappa)
  set.seed(91L)
  draws <- t(vapply(seq_len(20000L), function(s) {
    q <- numeric(8)
    q[1] <- stats::rnorm(1, 0, sigma_q)
    for (t in 2:8) q[t] <- rho * q[t - 1] + stats::rnorm(1, 0, sigma_q * sqrt(1 - rho^2))
    .column_values(q[panel$group] + stats::rnorm(32, 0, d / sqrt(kappa)), rep(0, 32), sd, d,
      kappa, panel)
  }, numeric(6)))
  colnames(draws) <- column_statistics
  testthat::expect_lt(abs(mean(draws[, "dispersion"]) - 1), 0.02)
  testthat::expect_lt(abs(mean(draws[, "monthly"]) - rho * sigma_q^2), 0.02)
  testthat::expect_lt(abs(mean(draws[, "contrast_dispersion"]) - 1), 0.01)
  testthat::expect_lt(abs(mean(draws[, "within_cell"])), 0.01)
  expected <- mean(vapply(seq_len(nrow(panel$neighbor_pairs)), function(k) {
    i <- panel$neighbor_pairs[k, ]
    month <- which(panel$group == panel$group[i[1]])
    C <- diag(4) - matrix(0.25, 4, 4)
    V <- C %*% diag(d[month]^2 / kappa) %*% C
    a <- match(i, month)
    V[a[1], a[2]] / sqrt(V[a[1], a[1]] * V[a[2], a[2]])
  }, numeric(1)))
  testthat::expect_lt(expected, 0)
  testthat::expect_lt(abs(mean(draws[, "neighbor"]) - expected), 0.01)
})

testthat::test_that("the tail summary counts ties, splits chains and flags constant indicators", {
  set.seed(5L)
  observed <- stats::rnorm(200)
  replicated <- observed + stats::rnorm(200, 0.3)
  replicated[1:3] <- observed[1:3]
  chain <- rep(1:2, 100)
  tail <- .column_tail(observed, replicated, chain)
  indicator <- as.numeric(replicated >= observed)
  testthat::expect_equal(tail$count, as.integer(sum(indicator)))
  testthat::expect_equal(tail$tail_probability, mean(indicator))
  mat <- cbind(indicator[chain == 1], indicator[chain == 2])
  testthat::expect_equal(tail$mcse, stats::sd(as.vector(mat)) / sqrt(.ess(.split_chains(mat))))
  testthat::expect_equal(tail$mcse_status, "estimated over chains")
  none <- .column_tail(rep(1, 200), rep(0, 200), chain)
  testthat::expect_equal(none$count, 0L)
  testthat::expect_true(is.na(none$mcse))
  testthat::expect_equal(none$mcse_status, "not estimable: constant indicator")
  testthat::expect_equal(c(none$clopper_pearson_lower, none$clopper_pearson_upper),
    stats::binom.test(0, 200)$conf.int[1:2])
  testthat::expect_error(.column_tail(1:3, 1:3, c(1, 1, 2)), "equal lengths")
})

testthat::test_that("the public checks drop the observed low-cost features only", {
  ppc <- data.frame(source = c("NO2 background", "Low-cost network, all devices",
    "Low-cost network, calibration pairs"), feature = "mean_log", observed = c(1, 2, 3), ppp = 0.5,
    stringsAsFactors = FALSE)
  p <- public_ppc(ppc)
  testthat::expect_identical(p$observed, c(1, NA, NA))
  testthat::expect_identical(p[-3L], ppc[-3L])
})

testthat::test_that("the calibration factor uses the within-month variance of the log BAM values", {
  input <- list(calibration = data.frame(month = c("2024-09", "2024-09", "2024-10", "2024-10", "2024-10"),
    bam = exp(c(1, 2, 0, 1, 2)), stringsAsFactors = FALSE))
  fit <- list(model = "m", unit = "u", draws = cbind("b_M[1]" = c(1.2, 1.6)))
  f <- calibration_factor(fit, input)
  s2 <- (0.5 + 2) / 3
  testthat::expect_equal(f$var_within_log_bam, s2)
  testthat::expect_equal(f$factor_at_mean, exp(1.4 * 0.4 * s2 / 2))
  testthat::expect_equal(f$factor_mean, mean(exp(c(1.2 * 0.2, 1.6 * 0.6) * s2 / 2)))
  testthat::expect_identical(c(f$n_pairs, f$n_months), c(5L, 2L))
})

testthat::test_that("posterior_predictive_checks produces the 55 source features", {
  skip_without_moqa()
  set.seed(20261131)
  x <- test_case("all-data", "final")
  draws <- do.call(rbind, replicate(8L, dense_values(dense_random_point(x$prep$data)), simplify = FALSE))
  fit <- list(prep = x$prep, draws = draws, seed = 126L)
  ppc <- posterior_predictive_checks(fit)
  testthat::expect_identical(nrow(ppc), 55L)
  testthat::expect_identical(
    names(ppc),
    c("source", "feature", "observed", "rep_q05", "rep_median", "rep_q95", "ppp", "flagged")
  )
  testthat::expect_true(all(is.finite(ppc$observed)) && all(is.finite(ppc$ppp)))
  testthat::expect_true(all(ppc$ppp >= 0 & ppc$ppp <= 1))
  cc <- column_checks(fit, x$input)
  testthat::expect_identical(cc$summary$statistic, column_statistics)
  testthat::expect_identical(dim(cc$draws), c(8L, 14L))
  testthat::expect_true(all(cc$summary$count >= 0L & cc$summary$count <= 8L))
})