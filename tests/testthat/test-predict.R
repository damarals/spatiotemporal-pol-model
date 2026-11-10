# The predictive of predict_held() (joint_moments() on the TROPOMI pseudo-rows
# and the conditioning) against a brute-force dense conditional on every
# training row, every TROPOMI row included.

testthat::test_that("the predictive is the dense conditional law of the held rows", {
  skip_without_moqa()
  set.seed(20261117)
  seed <- 115L
  cases <- list(c("new-devices-1", "full"), c("new-traffic-1", "full"), c("omitted-months", "full"),
    c("withheld-background", "final"), c("year-2025", "final"), c("new-traffic-2", "final"))
  for (case in cases) {
    info <- paste(case, collapse = " ")
    x <- test_case(case[1], case[2])
    prep <- x$prep
    d <- prep$data
    p <- dense_random_point(d)
    train_rows <- dense_train_rows(x)
    held_rows <- dense_held_rows(x)
    law <- dense_law(x, dense_bind_rows(train_rows, held_rows), p)
    n_dense <- length(train_rows$t)
    o <- seq_len(n_dense)
    h <- n_dense + seq_len(nrow(prep$held))
    testthat::expect_identical(n_dense, prep$n_ground + d$n_col, info = info)
    # The ground and held blocks of joint_moments() are the dense law's.
    m <- joint_moments(.draw(dense_values(p), 1L), prep)
    g <- seq_len(prep$n_ground)
    hp <- prep$n_train + seq_len(nrow(prep$held))
    testthat::expect_equal(m$mean[c(g, hp)], law$mean[c(g, h)], tolerance = 1e-12, info = info)
    testthat::expect_equal(m$V[c(g, hp), c(g, hp)], law$V[c(g, h), c(g, h)], tolerance = 1e-12, info = info)
    # Brute-force conditioning on every training row.
    y <- dense_train_values(x)
    gain <- t(solve(law$V[o, o], law$V[o, h]))
    mu <- law$mean[h] + drop(gain %*% (y - law$mean[o]))
    C <- law$V[h, h] - gain %*% law$V[o, h]
    y_held <- log(prep$held$observed)
    predicted <- predict_held(list(prep = prep, draws = dense_values(p), seed = seed))
    testthat::expect_equal(unname(predicted$log_density[, 1L]),
      stats::dnorm(y_held, mu, sqrt(diag(C)), log = TRUE) - y_held, tolerance = 1e-9, info = info)
    z <- .with_seed(seed + 1L, stats::rnorm(length(h)))
    testthat::expect_equal(unname(predicted$draws[, 1L]),
      exp(mu + drop(t(chol(C + diag(1e-10, length(h)))) %*% z)), tolerance = 1e-9, info = info)
    testthat::expect_gt(max(abs(mu - law$mean[h])), 1e-3)
    if (case[1] == "new-traffic-1") {
      # TROPOMI informs the held NO2 rows through e_N when lambda_T != 0.
      og <- seq_len(prep$n_ground)
      gain_g <- t(solve(law$V[og, og], law$V[og, h]))
      mu_g <- law$mean[h] + drop(gain_g %*% (y[og] - law$mean[og]))
      testthat::expect_gt(max(abs(mu[held_rows$type == 1L] - mu_g[held_rows$type == 1L])), 1e-4)
      testthat::expect_false(any(held_rows$site %in% prep$train_sites))
      testthat::expect_true(any(held_rows$type == 2L & held_rows$nugget))
    }
    if (case[1] == "new-devices-1") {
      testthat::expect_true(all(held_rows$type == 4L))
      testthat::expect_false(any(held_rows$device %in% prep$train_devices))
    }
  }
})

testthat::test_that("prediction leaves the caller's random stream untouched", {
  skip_without_moqa()
  set.seed(5)
  x <- test_case("omitted-months", "final")
  p <- dense_random_point(x$prep$data)
  set.seed(9)
  before <- .Random.seed
  kind <- RNGkind()
  predict_held(list(prep = x$prep, draws = dense_values(p), seed = 111L))
  testthat::expect_identical(.Random.seed, before)
  testthat::expect_identical(RNGkind(), kind)
})
