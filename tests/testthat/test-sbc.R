# SBC pieces: the prior draw is the program's parameter set within its supports;
# at a fixed point the simulator's Monte Carlo moments match the dense law of
# helper-dense.R; the rank statistic and the uniformity summary.

testthat::test_that("the prior draw is the program's parameter set within its supports", {
  skip_without_moqa()
  set.seed(20261129)
  x <- test_case("omitted-months", "final")
  d <- x$prep$data
  testthat::expect_identical(
    c(d$use_column, d$col_temporal, d$use_nugget, d$field_P, d$free_bA, d$use_kappa,
      d$couple_space, d$Q_N, d$sign_constraint),
    c(1L, 0L, 0L, 0L, 1L, 0L, 0L, 0L, 1L)
  )
  reference <- dense_random_point(d)
  for (k in 1:20) {
    p <- draw_prior(d)
    testthat::expect_setequal(names(p), names(reference))
    testthat::expect_identical(lengths(p)[names(reference)], lengths(reference))
    testthat::expect_true(all(unlist(p[c(
      "s_N", "s_P", "sig_N", "sig_NT", "sig_P", "sig_A", "sd_N",
      "tau_M", "sig_M", "sig_C", "tau_m", "sigma_q", "sigma_beta"
    )]) > 0))
    testthat::expect_true(p$b_A > 0.5 && p$b_A < 2.5 && p$b_M > 0.5 && p$b_M < 2 && p$lambda_S >= 0)
    testthat::expect_true(all(c(p$rho_N, p$rho_P) > 0 & c(p$rho_N, p$rho_P) < 1) &&
      abs(p$r) < 1 && abs(p$rho_q) < 1)
    testthat::expect_false("lambda_T" %in% names(p))
    testthat::expect_length(p$b_cell, 21L)
    testthat::expect_setequal(names(.flatten(p)), colnames(dense_values(p)))
  }
  draws <- replicate(4000L, draw_prior(d)$sig_N)
  testthat::expect_lt(abs(mean(draws > 1) - 0.05), 0.012)
  unsigned <- d
  unsigned$sign_constraint <- 0L
  signs <- replicate(200L, sign(draw_prior(unsigned)$lambda_S))
  testthat::expect_true(any(signs < 0) && any(signs > 0))
  d1 <- d
  d1$col_temporal <- 1L
  testthat::expect_true("lambda_T" %in% names(draw_prior(d1)))
})

testthat::test_that("the simulator's Monte Carlo moments match the dense law at a fixed point", {
  skip_without_moqa()
  set.seed(20261130)
  x <- test_case("omitted-months", "final")
  d <- x$prep$data
  p <- dense_random_point(d)
  p$lambda_S <- abs(p$lambda_S)
  rows <- dense_train_rows(x)
  law <- dense_law(x, rows, p)
  n <- length(law$mean)
  testthat::expect_identical(n, d$n + d$n_col)
  R <- 400L
  sims <- matrix(NA_real_, R, n)
  cal <- matrix(NA_real_, R, d$n_C)
  for (k in seq_len(R)) {
    s <- simulate_data(d, p)
    sims[k, ] <- c(s$ly, s$y_col)
    cal[k, ] <- s$lmoqa_C
  }
  z_mean <- (colMeans(sims) - law$mean) / sqrt(diag(law$V) / R)
  testthat::expect_lt(max(abs(z_mean)), 5)
  emp <- stats::cov(sims)
  z_cov <- (emp - law$V) / sqrt((outer(diag(law$V), diag(law$V)) + law$V^2) / R)
  upper <- z_cov[upper.tri(z_cov, diag = TRUE)]
  testthat::expect_lt(mean(abs(upper) > 3), 0.01)
  testthat::expect_lt(max(abs(upper)), 8)

  column <- rows$type == 6L
  testthat::expect_true(all(law$V[!column, column] == 0))
  testthat::expect_lt(max(abs(z_cov[!column, column])), 6)

  mean_C <- p$a_M + p$b_M * d$lbam_C
  V_C <- diag(p$sig_C^2, d$n_C) + matrix(p$tau_M^2, d$n_C, d$n_C)
  testthat::expect_lt(max(abs((colMeans(cal) - mean_C) / sqrt(diag(V_C) / R))), 5)
  z_C <- (stats::cov(cal) - V_C) / sqrt((outer(diag(V_C), diag(V_C)) + V_C^2) / R)
  testthat::expect_lt(mean(abs(z_C[upper.tri(z_C, diag = TRUE)]) > 3), 0.02)

  d2 <- .replace_observations(d, simulate_data(d, p))
  testthat::expect_identical(length(d2$ly), d$n)
  testthat::expect_identical(dim(d2$y_col), d$n_col)
  testthat::expect_identical(dim(d2$lmoqa_C), d$n_C)
  testthat::expect_false(identical(d2$ly, d$ly))
})

testthat::test_that("rank bins and the overview count what the text quotes", {
  table <- data.frame(rep = 1:20, divergences = rep(c(0, 3), 10), rank_a = rep(c(0, 99), 10),
    rank_b = 0:19 * 5)
  bins <- sbc_bins(table)
  testthat::expect_identical(bins$parameter, c("a", "b"))
  testthat::expect_identical(unname(unlist(bins[1L, paste0("bin_", 1:10)])), c(10L, rep(0L, 8L), 10L))
  testthat::expect_identical(unname(unlist(bins[2L, paste0("bin_", 1:10)])), rep(2L, 10L))
  testthat::expect_identical(bins$expected, c(2, 2))
  o <- sbc_overview(table, data.frame(p_value = c(0.01, 0.5)))
  testthat::expect_identical(c(o$replications, o$with_divergences, o$parameters, o$p_below_005), c(20L, 10L, 2L, 1L))
  testthat::expect_equal(c(o$expected_p_below_005, o$expected_per_bin, o$coverage90_se), c(0.1, 2, sqrt(0.09 / 20)))
})

testthat::test_that("the rank statistic counts thinned draws below the truth", {
  draws <- cbind(a = seq_len(4000), b = rev(seq_len(4000)))
  testthat::expect_identical(unname(.sbc_rank(draws, c(a = 2000.5, b = 2000.5))), c(50, 49))
  testthat::expect_identical(unname(.sbc_rank(draws, c(a = 0.5, b = 5000))), c(0, 99))
  small <- cbind(a = seq_len(600))
  testthat::expect_identical(unname(.sbc_rank(small, c(a = 300.5))), 50)

  table <- data.frame(divergences = rep(c(0L, 3L), 50L), max_rhat = 1.001, min_ess_bulk = 1000,
    rank_a = 0:99, rank_b = rep(0, 100L), check.names = FALSE)
  s <- summarise_sbc(table)
  testthat::expect_identical(s$parameter, c("a", "b"))
  testthat::expect_equal(s$chisq, c(0, 900))
  testthat::expect_equal(s$p_value, stats::pchisq(c(0, 900), 9, lower.tail = FALSE))
  testthat::expect_equal(s$p_holm, stats::p.adjust(s$p_value, "holm"))
  testthat::expect_equal(s$mean_u, c(0.5, 0.005))
  testthat::expect_equal(s$coverage90, c(0.9, 0))
  testthat::expect_identical(s$max_bin, c(10L, 100L))
  testthat::expect_identical(s$n_clean, c(50L, 50L))
  testthat::expect_equal(s$coverage90_clean, c(mean(seq(0, 98, 2) >= 5 & seq(0, 98, 2) <= 94), 0))
})
