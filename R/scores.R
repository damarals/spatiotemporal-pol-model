# Scores of held-out predictive draws on the concentration scale: CRPS, negative
# log score (NLS), bias, squared and absolute error, PIT, central 50/80/95%
# intervals with coverage and width, and the Monte Carlo standard error of the
# CRPS; then means within locations (or devices) with equal weight per location.

# Columns of the per-target score table.
score_columns <- c("model", "unit", "observation_id", "pollutant", "physical_location_id", "month",
  "observed", "source", "device_id", "site_class", "predictive_mean", "predictive_sd", "bias",
  "squared_error", "absolute_error", "crps", "nls", "pit", "lower_50", "upper_50", "coverage_50",
  "width_50", "lower_80", "upper_80", "coverage_80", "width_80", "lower_95", "upper_95",
  "coverage_95", "width_95", "crps_mcse")

# Metrics averaged by summarise_locations().
location_metrics <- c("predictive_mean", "predictive_sd", "bias", "squared_error", "absolute_error",
  "crps", "nls", "pit", "coverage_50", "width_50", "coverage_80", "width_80", "coverage_95",
  "width_95")

#' Per-target scores of one fit's predictions.
#'
#' @param held the held targets (partition()$held or fit$prep$held).
#' @param predicted predict_held(): draws and log_density (targets x draws).
#' @param model,unit labels written in the table.
#' @return data.frame with score_columns, one row per target in held order.
score_targets <- function(held, predicted, model, unit) {
  draws <- predicted$draws
  log_density <- predicted$log_density
  if (!identical(rownames(draws), held$observation_id) || !identical(dim(draws), dim(log_density)))
    stop("Draws and log densities must be targets x draws in held order.", call. = FALSE)
  if (any(!is.finite(draws)) || any(!is.finite(log_density)))
    stop("Draws and log densities must be finite.", call. = FALSE)
  weights <- .equal_weights(ncol(draws))
  observed <- held$observed
  scored <- data.frame(model = model, unit = unit, held[c("observation_id", "pollutant",
    "physical_location_id", "month", "observed", "source", "device_id", "site_class")],
    stringsAsFactors = FALSE)
  for (i in seq_len(nrow(draws))) {
    mean <- stats::weighted.mean(draws[i, ], weights)
    scored$predictive_mean[i] <- mean
    scored$predictive_sd[i] <- sqrt(stats::weighted.mean((draws[i, ] - mean)^2, weights))
  }
  scored$bias <- scored$predictive_mean - observed
  scored$squared_error <- scored$bias^2
  scored$absolute_error <- abs(scored$bias)
  scored$crps <- .crps(observed, draws, weights)
  scored$nls <- .nls(log_density, weights)
  scored$pit <- as.numeric((draws < observed) %*% weights) + 0.5 * as.numeric((draws == observed) %*% weights)
  for (level in c(50L, 80L, 95L)) {
    alpha <- (1 - level / 100) / 2
    bounds <- t(vapply(seq_len(nrow(draws)), function(i)
      .quantile_type1(draws[i, ], weights, c(alpha, 1 - alpha)), numeric(2)))
    scored[[paste0("lower_", level)]] <- bounds[, 1L]
    scored[[paste0("upper_", level)]] <- bounds[, 2L]
    scored[[paste0("coverage_", level)]] <- observed >= bounds[, 1L] & observed <= bounds[, 2L]
    scored[[paste0("width_", level)]] <- bounds[, 2L] - bounds[, 1L]
  }
  scored$crps_mcse <- .crps_mcse(observed, draws, weights)
  rownames(scored) <- NULL
  scored[score_columns]
}

#' Per-target scores of the fit whose chain CSVs are in dir.
score_fit <- function(input, unit, model, dir) {
  fit <- read_fit(input, unit, model, dir)
  score_targets(fit$prep$held, predict_held(fit), model, unit)
}

#' Equal-location summary of one pollutant's scored targets: each metric is
#' averaged within a location (a device for low-cost targets), then over
#' locations; RMSE is the square root of the mean squared error.
summarise_locations <- function(scored) {
  if (length(unique(scored$pollutant)) != 1L)
    stop("Summarise one pollutant at a time.", call. = FALSE)
  locations <- split(seq_len(nrow(scored)), scored$physical_location_id)
  means <- vapply(location_metrics, function(metric)
    mean(vapply(locations, function(rows) mean(scored[[metric]][rows]), numeric(1))), numeric(1))
  data.frame(pollutant = as.character(scored$pollutant[1L]), n_targets = nrow(scored),
    n_locations = length(locations), rmse = sqrt(means[["squared_error"]]), as.list(means),
    row.names = NULL, check.names = FALSE)
}

# Summaries reported per unit group: the pooled folds and each fold.
summary_groups <- data.frame(
  group = c("omitted-months", "new-traffic", "new-traffic-1", "new-traffic-2", "new-traffic-3",
    "withheld-background", "whole-years", "year-2023", "year-2024", "year-2025",
    "new-devices", "new-devices-1", "new-devices-2", "new-devices-3"),
  units = c("omitted-months", "new-traffic-1;new-traffic-2;new-traffic-3", "new-traffic-1",
    "new-traffic-2", "new-traffic-3", "withheld-background", "year-2023;year-2024;year-2025",
    "year-2023", "year-2024", "year-2025", "new-devices-1;new-devices-2;new-devices-3",
    "new-devices-1", "new-devices-2", "new-devices-3"),
  stringsAsFactors = FALSE)

#' Equal-location summaries of a model's scored targets by unit group and
#' pollutant. Reference targets are summarised for the reference units and
#' low-cost targets for the new-device units; a group needs all its units.
summarise_units <- function(scored) {
  rows <- list()
  for (i in seq_len(nrow(summary_groups))) {
    members <- strsplit(summary_groups$units[i], ";", fixed = TRUE)[[1L]]
    if (!all(members %in% scored$unit)) next
    devices <- startsWith(summary_groups$group[i], "new-devices")
    for (pollutant in if (devices) "PM25" else c("NO2", "PM25")) {
      source <- if (devices) "moqa" else if (pollutant == "NO2") "ground_no2" else "ground_pm25"
      selected <- scored[scored$unit %in% members & scored$pollutant == pollutant &
        scored$source == source, , drop = FALSE]
      if (!nrow(selected)) next
      rows[[length(rows) + 1L]] <- cbind(data.frame(model = selected$model[1L],
        group = summary_groups$group[i], units = summary_groups$units[i]),
        summarise_locations(selected))
    }
  }
  do.call(rbind, rows)
}

.equal_weights <- function(n) {
  w <- rep(1, n)
  scaled <- w / max(w)
  scaled / sum(scaled)
}

# CRPS of the weighted empirical predictive distribution.
.crps <- function(observed, draws, weights) {
  scoringRules::crps_sample(y = observed, dat = draws, method = "edf",
    w = matrix(rep(weights, each = nrow(draws)), nrow(draws)))
}

# NLS: minus the log of the weighted mean predictive density, by log-sum-exp.
.nls <- function(log_density, weights) {
  terms <- sweep(log_density, 2L, log(weights), "+")
  apply(terms, 1L, function(x) {
    shift <- max(x)
    -(shift + log(sum(exp(x - shift))))
  })
}

# Monte Carlo standard error of the CRPS from 8 consecutive batches of draws.
.crps_mcse <- function(observed, draws, weights) {
  batches <- 8L
  size <- ncol(draws) %/% batches
  if (size < 2L) stop("Too few draws for the Monte Carlo standard error.", call. = FALSE)
  batch_crps <- vapply(seq_len(batches), function(b) {
    i <- (b - 1L) * size + seq_len(size)
    w <- weights[i] / max(weights[i])
    .crps(observed, draws[, i, drop = FALSE], w / sum(w))
  }, numeric(nrow(draws)))
  batch_crps <- matrix(batch_crps, nrow(draws))
  apply(batch_crps, 1L, stats::sd) / sqrt(batches)
}

# Weighted type-1 quantile (inverse of the empirical distribution function).
.quantile_type1 <- function(x, w, probs) {
  o <- order(x)
  x <- x[o]
  Fx <- cumsum(w[o]) / sum(w[o])
  if (anyDuplicated(x)) {
    last <- rev(!duplicated(rev(x)))
    x <- x[last]
    Fx <- Fx[last]
  }
  if (length(x) == 1L) return(rep(x, length(probs)))
  stats::approx(Fx, x, xout = probs, ties = "ordered", rule = 2, method = "constant", f = 1)$y
}
