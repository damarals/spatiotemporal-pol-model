# The program's log density against the dense Normal law of the model, written
# row by row (every TROPOMI row), up to the prior constants that `~` drops; the
# Stan data flags of every model; the ported convergence diagnostics.

check_program <- function(x, points, info) {
  d <- x$prep$data
  stan <- vapply(points, function(p) stan_log_density(d, p), numeric(1))
  dense <- vapply(points, function(p) dense_log_density(x, p), numeric(1))
  dropped <- vapply(points, dense_dropped_constant, numeric(1))
  testthat::expect_true(all(is.finite(c(stan, dense))), info = info)
  testthat::expect_equal(stan - dense, dropped, tolerance = 1e-9, info = info)
  invisible(dense)
}

testthat::test_that("the program is the dense law of the full and final models with TROPOMI", {
  skip_without_moqa()
  skip_without_cmdstan()
  set.seed(20261115)
  for (case in list(c("new-devices-1", "full"), c("new-traffic-1", "full"), c("new-traffic-1", "final"))) {
    x <- test_case(case[1], case[2])
    d <- x$prep$data
    testthat::expect_identical(c(d$use_column, d$n_col, d$n_cells, d$K_col, d$K), c(1L, 1584L, 22L, 4L, 4L))
    testthat::expect_true(all(d$count_col == 22L))
    testthat::expect_gt(sum(d$type == 4L), 0L)
    testthat::expect_gt(d$n_C, 0L)
    testthat::expect_gt(sum(d$type == 1L & d$cls == 3L), 0L)
    # The pseudo-rows follow the ground rows, one per month.
    rows <- x$prep$rows
    testthat::expect_identical(which(rows$type == 5L), x$prep$n_ground + 1:72)
    points <- replicate(3L, dense_random_point(d), simplify = FALSE)
    dense <- check_program(x, points, paste(case, collapse = " "))
    testthat::expect_gt(diff(range(dense)), 10)
  }
})

testthat::test_that("a whole-year unit with months without TROPOMI rows equals the dense law", {
  skip_without_moqa()
  skip_without_cmdstan()
  set.seed(20261123)
  x <- test_case("year-2024", "full")
  d <- x$prep$data
  testthat::expect_identical(c(d$use_column, d$n_col, d$T), c(1L, 1056L, 59L))
  testthat::expect_true(any(d$count_col == 0L))
  check_program(x, replicate(2L, dense_random_point(d), simplify = FALSE), "year-2024")
})

testthat::test_that("every model's flags keep the program equal to the dense law", {
  skip_without_moqa()
  skip_without_cmdstan()
  set.seed(20261116)
  for (model in models$model) {
    m <- model_spec(model)
    unit <- model_units(model)[1L]
    x <- test_case(unit, model)
    d <- x$prep$data
    testthat::expect_identical(c(d$use_column, d$col_temporal, d$sign_constraint, d$field_N, d$field_P,
      d$coupled, d$use_kappa, d$n_g, d$use_nugget, d$acag_seasonal, d$free_bA, d$K, d$quartic),
      c(m$tropomi, m$temporal_loading, m$sign_constraint, m$no2_profile, m$pm25_profile, m$r,
        m$amplitude, m$classes, m$deviation, m$seasonal, m$free_gain, m$knots,
        as.integer(m$priors == "quartic")), info = model)
    testthat::expect_identical(d$couple_space, as.integer(m$profile_cor && m$no2_profile && m$pm25_profile),
      info = model)
    testthat::expect_identical(c(d$Q_N, d$Q_P, d$P_N, d$P_P), c(m$no2_cov, m$pm25_cov, 2L * m$era5,
      2L * m$era5), info = model)
    testthat::expect_identical(any(d$type == 3L), m$satpm25 == 1L, info = model)
    testthat::expect_identical(any(d$type == 4L), m$low_cost == 1L, info = model)
    check_program(x, list(dense_random_point(d)), model)
  }
  testthat::expect_error(test_case("new-devices-1", "reduced-without-low-cost"), "Held low-cost rows")
})

testthat::test_that("the sign constraint only bounds lambda_S", {
  skip_without_moqa()
  skip_without_cmdstan()
  set.seed(20261121)
  x <- test_case("new-traffic-1", "full")
  p <- dense_random_point(x$prep$data)
  free <- x$prep$data
  free$sign_constraint <- 0L
  testthat::expect_equal(stan_log_density(free, p), stan_log_density(x$prep$data, p), tolerance = 1e-12)
  p$lambda_S <- -p$lambda_S
  testthat::expect_true(is.finite(stan_log_density(free, p)))
  testthat::expect_error(stan_log_density(x$prep$data, p), "log_prob failed")
})

testthat::test_that("the quartic scale has the declared tail probability", {
  a <- quartic_scale(c(0.6, 1))
  # (sigma / a)^4 ~ Gamma(1/4, 1): Pr(sigma > U) = Pr(G > (U / a)^4).
  testthat::expect_equal(stats::pgamma((c(0.6, 1) / a)^4, 1 / 4, lower.tail = FALSE), c(0.05, 0.05),
    tolerance = 1e-12)
  density <- function(s, a) 4 / (a * gamma(1 / 4)) * exp(-(s / a)^4)
  testthat::expect_equal(stats::integrate(density, 0, Inf, a = a[1])$value, 1, tolerance = 1e-8)
  testthat::expect_equal(stats::integrate(density, 0.6, Inf, a = a[1])$value, 0.05, tolerance = 1e-8)
})

testthat::test_that("R-hat and bulk ESS follow the rank-normalized definitions", {
  set.seed(3)
  iid <- matrix(stats::rnorm(4000), 1000, 4)
  testthat::expect_lt(abs(rank_rhat(iid) - 1), 0.01)
  testthat::expect_gt(bulk_ess(iid), 3000)
  shifted <- iid + rep(c(0, 0, 0, 3), each = 1000)
  testthat::expect_gt(rank_rhat(shifted), 1.3)
  ar <- matrix(stats::filter(stats::rnorm(4000), 0.9, method = "recursive"), 1000, 4)
  testthat::expect_lt(bulk_ess(ar), 600)
})
