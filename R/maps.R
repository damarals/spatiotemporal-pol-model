# Concentration maps of a fit on the 0.01 degree SatPM2.5 cells. In each
# posterior draw the 2T monthly latent log levels ell_p(t) (background NO2 and
# PM2.5 at no location, without observation noise) are drawn jointly from their
# conditional law given the training rows, and the map value of cell s is
#   C_p(s, t) = exp{ell_p(t) + spatial_p(s)},
#   spatial_N(s) = log R_N(s)                   (the NO2 profile; no class offset),
#   spatial_P(s) = X_P(s)' theta_P + log R_P(s)  (log R_P only with a PM2.5 profile),
# with X_P the fit's PM2.5 covariate at the cell, transformed as at the locations,
# clipped to the range of the 22 locations and standardized with their mean and SD.

#' The map grid: SatPM2.5 cells with their raw covariates.
#'
#' @param dir data directory.
#' @param files data_files(dir), possibly without names.
#' @return data.frame: cell, lon, lat, weight, x_m, y_m and the covariates of
#'   grid-covariates.csv.
read_grid <- function(dir = "data", files = data_files(dir)) {
  files <- .data_paths(files)
  cells <- utils::read.csv(files[["acag_cells"]], stringsAsFactors = FALSE)
  covariates <- utils::read.csv(files[["grid_covariates"]], stringsAsFactors = FALSE)
  if (!identical(cells$cell, covariates$cell))
    stop("acag-cells.csv and grid-covariates.csv must list the same cells in order.", call. = FALSE)
  cbind(cells, covariates[setdiff(names(covariates), "cell")])
}

#' Map products of a fit.
#'
#' @param fit fit_model() or read_fit() (the final model fitted to all data).
#' @param input read_inputs().
#' @param grid read_grid().
#' @return list: annual (cell, pollutant, mean, q025, q975 of the 2020-2025 mean
#'   of the monthly values), monthly (cell, month, pollutant, mean, q025, q975),
#'   level (posterior summaries of exp(ell_p(t)) and of ell_p(t) per month),
#'   clipping (site range of each PM2.5 covariate and the cells outside it) and
#'   clipped (cells x covariates, TRUE where the covariate was clipped).
compute_maps <- function(fit, input, grid) {
  needed <- c("cell", "x_m", "y_m")
  if (!all(needed %in% names(grid)))
    stop("The grid table needs ", paste(needed, collapse = ", "), ".", call. = FALSE)
  latent <- .latent_levels(fit)
  spatial <- .spatial_terms(fit, grid, input$site_covariates)
  products <- .map_products(latent, spatial, fit$prep, grid$cell)
  c(products, list(clipping = spatial$range, clipped = spatial$clipped))
}

#' Cell values of the maps: the 2020-2025 mean per cell and pollutant.
#'
#' @param maps compute_maps().
#' @param grid the grid passed to compute_maps().
#' @return data.frame: cell, lon, lat, weight, pollutant, mean, q025, q975,
#'   clipped (TRUE where a PM2.5 covariate lies outside the range of the 22
#'   locations; FALSE for NO2, which has no covariate on the map).
map_cells <- function(maps, grid) {
  a <- maps$annual
  at <- match(a$cell, grid$cell)
  if (anyNA(at)) stop("The map cells differ from the grid cells.", call. = FALSE)
  clipped <- rowSums(maps$clipped) > 0
  data.frame(cell = a$cell, lon = grid$lon[at], lat = grid$lat[at], weight = grid$weight[at],
    pollutant = a$pollutant, mean = a$mean, q025 = a$q025, q975 = a$q975,
    clipped = a$pollutant == "PM25" & clipped[at], stringsAsFactors = FALSE)
}

#' Summary of a fit's maps per pollutant: range of the cell values of the
#' 2020-2025 mean and the ratio of its ends, area-weighted mean, latent municipal
#' level (mean over the months of the posterior mean of exp(ell_p(t))), median
#' width of the central 95 % interval relative to the mean, and clipped cells.
#'
#' @param maps compute_maps().
#' @param grid the grid passed to compute_maps().
#' @param model label written in the table.
map_summary <- function(maps, grid, model) {
  cells <- map_cells(maps, grid)
  do.call(rbind, lapply(c("NO2", "PM25"), function(p) {
    x <- cells[cells$pollutant == p, , drop = FALSE]
    level <- maps$level$level_mean[maps$level$pollutant == p]
    data.frame(model = model, pollutant = p, cell_min = min(x$mean), cell_max = max(x$mean),
      max_over_min = max(x$mean) / min(x$mean), area_mean = sum(x$weight * x$mean),
      latent_level = mean(level), median_relative_width = stats::median(relative_width(x)),
      cells_clipped = sum(x$clipped), stringsAsFactors = FALSE)
  }))
}

#' Width of the central 95 % interval relative to the mean, per map cell.
relative_width <- function(cells) (cells$q975 - cells$q025) / cells$mean

# Latent rows: ell_N(t) as a background NO2 row and ell_P(t) as a reference
# PM2.5 row, at no location and without observation noise.
.latent_rows <- function(prep) {
  T <- prep$data$T
  data.frame(observation_id = c(paste0("latent:N:", prep$calendar), paste0("latent:P:", prep$calendar)),
    type = rep(c(1L, 2L), each = T), t = rep(seq_len(T), 2L), site = NA_character_,
    class = rep(c(1L, 0L), each = T), device = "", x_m = NA_real_, y_m = NA_real_,
    observed = NA_real_, col = NA_integer_, latent = TRUE, stringsAsFactors = FALSE)
}

# Conditional mean and variance of the latent rows given the training rows, and
# one joint draw of them, per posterior draw (S x 2T matrices).
.latent_levels <- function(fit) {
  prep <- fit$prep
  values <- fit$draws
  o <- seq_len(prep$n_train)
  train <- prep$rows[o, , drop = FALSE]
  train$latent <- FALSE
  latent <- .latent_rows(prep)
  rows <- rbind(train, latent)
  h <- prep$n_train + seq_len(nrow(latent))
  y <- ifelse(train$type == 5L, train$observed, log(train$observed))
  S <- nrow(values)
  mean <- var <- draw <- matrix(NA_real_, S, nrow(latent), dimnames = list(NULL, latent$observation_id))
  .with_mt_seed(fit$seed + 2L, {
    for (j in seq_len(S)) {
      m <- joint_moments(.draw(values, j), prep, rows)
      R <- chol(m$V[o, o])
      alpha <- backsolve(R, forwardsolve(t(R), y - m$mean[o]))
      W <- forwardsolve(t(R), m$V[o, h, drop = FALSE])
      mu <- m$mean[h] + drop(crossprod(m$V[o, h, drop = FALSE], alpha))
      C <- m$V[h, h, drop = FALSE] - crossprod(W)
      mean[j, ] <- mu
      var[j, ] <- diag(C)
      L <- tryCatch(t(chol(C + diag(1e-10, nrow(C)))), error = function(e)
        stop("Conditional covariance of the latent levels is not positive definite in draw ", j, ".",
          call. = FALSE))
      draw[j, ] <- mu + drop(L %*% stats::rnorm(nrow(C)))
    }
  })
  if (any(!is.finite(c(mean, var, draw))) || any(var < -1e-8))
    stop("The latent levels are not finite.", call. = FALSE)
  list(rows = latent, mean = mean, var = var, draw = draw)
}

# Covariates `cols` at the points, transformed as at the locations, clipped to the
# range of the locations and standardized with their mean and SD.
.grid_covariates <- function(points, site_covariates, cols) {
  if (!length(cols))
    return(list(X = matrix(0, nrow(points), 0L), range = NULL, clipped = matrix(FALSE, nrow(points), 0L)))
  sites <- .site_covariates(site_covariates)
  grid <- .site_covariates(within(points, location_id <- NA_character_))
  X <- clipped <- matrix(0, nrow(points), length(cols), dimnames = list(NULL, cols))
  clipped <- clipped > 0
  range <- data.frame(covariate = cols, site_min = NA_real_, site_max = NA_real_, site_mean = NA_real_,
    site_sd = NA_real_, cells_below = NA_integer_, cells_above = NA_integer_, cells_clipped = NA_integer_,
    stringsAsFactors = FALSE)
  for (i in seq_along(cols)) {
    s <- as.numeric(sites[[cols[i]]])
    g <- as.numeric(grid[[cols[i]]])
    if (!length(s) || !length(g) || any(!is.finite(s)) || any(!is.finite(g)))
      stop("Covariate ", cols[i], " is missing or not finite at the locations or cells.", call. = FALSE)
    scale <- stats::sd(s)
    if (!(scale > 0)) stop("Covariate ", cols[i], " has no variation across the locations.", call. = FALSE)
    lo <- min(s)
    hi <- max(s)
    below <- g < lo
    above <- g > hi
    clipped[, i] <- below | above
    X[, i] <- (pmin(pmax(g, lo), hi) - mean(s)) / scale
    range[i, c("site_min", "site_max", "site_mean", "site_sd")] <- c(lo, hi, mean(s), scale)
    range[i, c("cells_below", "cells_above", "cells_clipped")] <- c(sum(below), sum(above), sum(clipped[, i]))
  }
  list(X = X, range = range, clipped = clipped)
}

# Spatial terms of the points per draw (points x S): log R_N(s) for background
# NO2 and X_P(s)' theta_P + log R_P(s) for PM2.5.
.spatial_terms <- function(fit, points, site_covariates) {
  prep <- fit$prep
  values <- fit$draws
  d <- prep$data
  S <- nrow(values)
  n <- nrow(points)
  W <- profile_support(cbind(points$x_m, points$y_m), prep$basis$knots, prep$basis$width_m)
  cov <- .grid_covariates(points, site_covariates, unname(prep$covariates[names(prep$covariates) == "pm25"]))
  if (ncol(cov$X) != d$Q_P) stop("The PM2.5 covariates of the map differ from the fit's.", call. = FALSE)
  spatial_N <- spatial_P <- matrix(0, n, S)
  for (j in seq_len(S)) {
    par <- .draw(values, j)
    if (d$field_N) spatial_N[, j] <- log(drop(W %*% exp(.log_profile(prep$basis, par$z_N, par$sd_N))))
    sp <- rep(0, n)
    if (d$field_P) {
      zc <- if (d$couple_space) par$rho_S * par$z_N + sqrt(1 - par$rho_S^2) * par$z_P else par$z_P
      sp <- log(drop(W %*% exp(.log_profile(prep$basis, zc, par$sd_P))))
    }
    if (d$Q_P > 0L) sp <- sp + drop(cov$X %*% par$theta_P)
    spatial_P[, j] <- sp
  }
  list(N = spatial_N, P = spatial_P, range = cov$range, clipped = cov$clipped)
}

# Posterior mean, 2.5 % and 97.5 % of C_p(s, t) per cell and month, of its
# 2020-2025 mean per cell (the mean within each draw, then over draws), and of
# exp(ell_p(t)) per month.
.map_products <- function(latent, spatial, prep, cells) {
  T <- prep$data$T
  S <- nrow(latent$draw)
  n <- length(cells)
  pollutant <- c(N = "NO2", P = "PM25")
  monthly <- level <- vector("list", 2L * T)
  annual <- vector("list", 2L)
  k <- 0L
  for (p in c("N", "P")) {
    sp <- spatial[[p]]
    acc <- matrix(0, n, S)
    for (t in seq_len(T)) {
      id <- paste0("latent:", p, ":", prep$calendar[t])
      ell <- latent$draw[, id]
      Cm <- exp(sweep(sp, 2L, ell, "+"))
      acc <- acc + Cm / T
      q <- apply(Cm, 1L, stats::quantile, c(0.025, 0.975))
      k <- k + 1L
      monthly[[k]] <- data.frame(cell = cells, month = prep$calendar[t], pollutant = pollutant[[p]],
        mean = rowMeans(Cm), q025 = q[1L, ], q975 = q[2L, ], stringsAsFactors = FALSE)
      lv <- exp(ell)
      level[[k]] <- data.frame(month = prep$calendar[t], pollutant = pollutant[[p]],
        level_mean = mean(lv), level_q025 = unname(stats::quantile(lv, 0.025)),
        level_q975 = unname(stats::quantile(lv, 0.975)), log_level_mean = mean(latent$mean[, id]),
        log_level_sd = sqrt(mean(latent$var[, id]) + stats::var(latent$mean[, id])),
        stringsAsFactors = FALSE)
    }
    q <- apply(acc, 1L, stats::quantile, c(0.025, 0.975))
    annual[[p]] <- data.frame(cell = cells, pollutant = pollutant[[p]], mean = rowMeans(acc),
      q025 = q[1L, ], q975 = q[2L, ], stringsAsFactors = FALSE)
  }
  list(annual = do.call(rbind, unname(annual)), monthly = do.call(rbind, monthly),
    level = do.call(rbind, level))
}

# Evaluates expr with the Mersenne-Twister stream of `seed` (R's default
# generator, which the maps and checks use), then restores the caller's RNG.
.with_mt_seed <- function(seed, expr) {
  kind <- RNGkind()
  had <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  if (had) previous <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
  on.exit({
    RNGkind(kind[1L], kind[2L], kind[3L])
    if (had) assign(".Random.seed", previous, envir = globalenv())
  }, add = TRUE)
  set.seed(seed, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
  expr
}
