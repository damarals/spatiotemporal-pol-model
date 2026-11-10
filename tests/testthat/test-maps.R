# Map products: at a location, the map's latent log mean equals the latent part
# of its predictive mean under the dense law; grid covariates are clipped to the
# range of the 22 locations; the products have the declared shape.

testthat::test_that("the map's latent log mean at a station is the latent part of its predictive mean", {
  skip_without_moqa()
  set.seed(20261127)
  x <- test_case("omitted-months", "final-no-tropomi")
  prep <- x$prep
  d <- prep$data
  testthat::expect_identical(
    c(d$use_column, d$use_nugget, d$field_P, d$field_N, d$Q_N, d$Q_P, d$free_bA, d$use_kappa),
    c(0L, 0L, 0L, 1L, 0L, 1L, 1L, 0L)
  )
  p <- dense_random_point(d)
  values <- dense_values(p)
  hq <- background_location
  site_row <- x$input$sites[x$input$sites$location_id == hq, , drop = FALSE]
  cov_row <- x$input$site_covariates[x$input$site_covariates$location_id == hq, , drop = FALSE]
  point_row <- cbind(site_row[, c("location_id", "x_m", "y_m")],
    cov_row[, setdiff(names(cov_row), "location_id")])
  testthat::expect_identical(nrow(point_row), 1L)
  t_N <- 10L
  t_P <- 20L

  train_rows <- dense_train_rows(x)
  x_held <- x
  x_held$prep$held <- data.frame(
    observation_id = c("map:N", "map:P"),
    pollutant = c("NO2", "PM25"),
    physical_location_id = hq,
    month = prep$calendar[c(t_N, t_P)],
    observed = 1,
    source = c("ground_no2", "ground_pm25"),
    device_id = "",
    site_class = "background",
    x_m = site_row$x_m,
    y_m = site_row$y_m,
    stringsAsFactors = FALSE
  )
  held_rows <- dense_held_rows(x_held)
  rows <- dense_bind_rows(train_rows, held_rows)
  law <- dense_law(x, rows, p)
  o <- seq_along(train_rows$t)
  h <- length(o) + 1:2
  y <- dense_train_values(x)
  gain <- t(solve(law$V[o, o], law$V[o, h]))
  mu <- law$mean[h] + drop(gain %*% (y - law$mean[o]))

  latent <- .latent_levels(list(prep = prep, draws = values, seed = -1L))
  testthat::expect_identical(dim(latent$mean), c(1L, 2L * d$T))
  spatial <- .spatial_terms(list(prep = prep, draws = values), point_row, x$input$site_covariates)
  map_N <- latent$mean[1L, paste0("latent:N:", prep$calendar[t_N])] + spatial$N[1L, 1L]
  map_P <- latent$mean[1L, paste0("latent:P:", prep$calendar[t_P])] + spatial$P[1L, 1L]
  testthat::expect_equal(unname(c(map_N, map_P)), unname(mu), tolerance = 1e-9)

  support <- .dense_support(site_row$x_m, site_row$y_m, prep$basis$knots, x$model$width_m)
  testthat::expect_equal(
    spatial$N[1L, 1L],
    log(drop(support %*% exp(dense_log_profiles(d, p)$N))),
    tolerance = 1e-12
  )
  testthat::expect_equal(
    spatial$P[1L, 1L],
    drop(held_rows$x_P[2L, , drop = FALSE] %*% p$theta_P),
    tolerance = 1e-12
  )
  testthat::expect_false(any(spatial$clipped))
  testthat::expect_gt(
    abs(latent$mean[1L, paste0("latent:P:", prep$calendar[t_P])] - law$mean[h[2L]] + spatial$P[1L, 1L]),
    1e-3
  )
  testthat::expect_lt(latent$var[1L, paste0("latent:P:", prep$calendar[t_P])], p$s_P^2)
  testthat::expect_gt(min(latent$var), 0)
})

testthat::test_that("the map summary gives ranges, weighted means, levels and relative widths", {
  grid <- data.frame(cell = 1:3, lon = 0, lat = 0, weight = c(0.5, 0.25, 0.25))
  maps <- list(
    annual = data.frame(cell = rep(1:3, 2L), pollutant = rep(c("NO2", "PM25"), each = 3L),
      mean = c(2, 4, 8, 10, 10, 10), q025 = c(1, 3, 6, 9, 8, 7), q975 = c(3, 5, 10, 11, 12, 13)),
    level = data.frame(pollutant = c("NO2", "NO2", "PM25", "PM25"), level_mean = c(3, 5, 9, 11)),
    clipped = matrix(c(FALSE, TRUE, FALSE), 3L))
  s <- map_summary(maps, grid, "m")
  testthat::expect_identical(s$pollutant, c("NO2", "PM25"))
  testthat::expect_equal(s$cell_min, c(2, 10))
  testthat::expect_equal(s$max_over_min, c(4, 1))
  testthat::expect_equal(s$area_mean, c(0.5 * 2 + 0.25 * 4 + 0.25 * 8, 10))
  testthat::expect_equal(s$latent_level, c(4, 10))
  testthat::expect_equal(s$median_relative_width, c(0.5, 0.4))
  testthat::expect_identical(s$cells_clipped, c(0L, 1L))
})

testthat::test_that("grid covariates are clipped to the 22-location range", {
  skip_without_moqa()
  input <- test_inputs()
  grid <- read_grid(file.path(repo_root, "data"))
  testthat::expect_identical(nrow(grid), 330L)
  sites <- .site_covariates(input$site_covariates)
  for (col in c("log_pop_1km", "log_veh_250", "log_bus_500", "log_mtps_100")) {
    cov <- .grid_covariates(grid, input$site_covariates, col)
    z <- (sites[[col]] - mean(sites[[col]])) / stats::sd(sites[[col]])
    testthat::expect_true(all(cov$X[, 1L] >= min(z) - 1e-12 & cov$X[, 1L] <= max(z) + 1e-12), info = col)
    testthat::expect_identical(cov$range$cells_below + cov$range$cells_above, cov$range$cells_clipped)
    testthat::expect_identical(cov$range$cells_clipped, sum(cov$clipped))
  }
  testthat::expect_identical(sum(.grid_covariates(grid, input$site_covariates, "log_pop_1km")$clipped), 55L)
  testthat::expect_identical(sum(.grid_covariates(grid, input$site_covariates, "log_veh_250")$clipped), 116L)
})

testthat::test_that("the map products have the declared shape and the annual mean is the mean of the months", {
  skip_without_moqa()
  set.seed(20261128)
  x <- test_case("omitted-months", "final")
  prep <- x$prep
  grid <- read_grid(file.path(repo_root, "data"))
  values <- rbind(dense_values(dense_random_point(prep$data)), dense_values(dense_random_point(prep$data)))
  fit <- list(prep = prep, draws = values, seed = 126L)
  products <- compute_maps(fit, x$input, grid)
  testthat::expect_identical(nrow(products$monthly), 330L * prep$data$T * 2L)
  testthat::expect_identical(nrow(products$annual), 660L)
  testthat::expect_identical(nrow(products$level), 2L * prep$data$T)
  testthat::expect_true(all(is.finite(products$monthly$mean)) && all(products$monthly$mean > 0))
  testthat::expect_true(all(products$monthly$q025 <= products$monthly$mean + 1e-9) &&
    all(products$monthly$q975 >= products$monthly$mean - 1e-9))
  for (p in c("NO2", "PM25")) {
    monthly <- products$monthly[products$monthly$pollutant == p, ]
    annual <- products$annual[products$annual$pollutant == p, ]
    by_cell <- tapply(monthly$mean, monthly$cell, mean)
    testthat::expect_equal(as.numeric(by_cell[as.character(annual$cell)]), annual$mean,
      tolerance = 1e-9, info = p)
  }
  no2 <- products$annual[products$annual$pollutant == "NO2", ]
  testthat::expect_gt(diff(range(no2$mean)), 1e-3)
  testthat::expect_identical(products$clipping$covariate, unname(prep$covariates[names(prep$covariates) == "pm25"]))
  testthat::expect_identical(products$clipping$cells_clipped, 3L)
  cells <- map_cells(products, grid)
  testthat::expect_identical(nrow(cells), 660L)
  testthat::expect_identical(sum(cells$clipped), 3L)
  testthat::expect_false(any(cells$clipped[cells$pollutant == "NO2"]))
})