# Posterior predictive checks of a fit to all data. Replicates come from
# simulate_data(), the integrated law the sampler targets: at each posterior
# draw the AR(1) residuals, TROPOMI's common error q_t, the low-cost monthly and
# device effects and the noise are drawn afresh. The checks are descriptive.

# The six TROPOMI statistics in the panel order of the figure.
column_statistics <- c("dispersion", "maximum", "monthly", "contrast_dispersion", "within_cell", "neighbor")

#' Marginal posterior predictive check: 55 features of the training data by
#' source (mean, SD, min and max of the log values, mean and max on the
#' concentration scale, lag-1 autocorrelation of the monthly series where there
#' is one; for TROPOMI, the monthly means over the cells), over every fourth
#' posterior draw.
#'
#' @param fit fit_model() or read_fit().
#' @return data.frame: source, feature, observed, rep_q05, rep_median, rep_q95,
#'   ppp (share of replicates at or above the observed value, ties half) and
#'   flagged (ppp < 0.05 or > 0.95).
posterior_predictive_checks <- function(fit) {
  prep <- fit$prep
  d <- prep$data
  values <- fit$draws
  keep <- seq(4L, nrow(values), by = 4L)
  background <- match(background_location, prep$train_sites)
  groups <- list(
    "NO2 background" = which(d$type == 1L & d$cls == 1L),
    "NO2 avenue" = which(d$type == 1L & d$cls == 2L),
    "NO2 bus terminal" = which(d$type == 1L & d$cls == 3L),
    "PM2.5 reference, all locations" = which(d$type == 2L),
    "PM2.5 reference, background location" = which(d$type == 2L & d$site == background),
    "SatPM2.5 municipal mean" = which(d$type == 3L),
    "Low-cost network, all devices" = which(d$type == 4L))
  groups <- groups[lengths(groups) > 0L]
  observed <- .ppc_features(d, groups, d$ly, as.numeric(d$y_col), as.numeric(d$lmoqa_C))
  replicated <- .with_mt_seed(20261108L, t(vapply(keep, function(j) {
    sim <- simulate_data(d, .draw(values, j))
    .ppc_features(d, groups, sim$ly, sim$y_col, sim$lmoqa_C)
  }, numeric(length(observed)))))
  ppp <- vapply(seq_along(observed), function(j) {
    x <- replicated[, j]
    mean(x > observed[j], na.rm = TRUE) + 0.5 * mean(x == observed[j], na.rm = TRUE)
  }, numeric(1))
  q <- t(apply(replicated, 2L, stats::quantile, probs = c(0.05, 0.5, 0.95), na.rm = TRUE))
  parts <- do.call(rbind, strsplit(names(observed), " | ", fixed = TRUE))
  data.frame(source = parts[, 1L], feature = parts[, 2L], observed = unname(observed),
    rep_q05 = q[, 1L], rep_median = q[, 2L], rep_q95 = q[, 3L], ppp = ppp,
    flagged = ppp < 0.05 | ppp > 0.95, row.names = NULL, stringsAsFactors = FALSE)
}

#' TROPOMI column checks over every posterior draw. At draw s the row of cell i
#' and month t has the fitted mean mu_it (without q_t), residual r_it = Y_it -
#' mu_it and SD sd_it = sqrt(sigma_q^2 + d_it^2 / kappa_T); the contrast e_it =
#' r_it - mean_j r_jt over the month's n_t cells removes q_t exactly and has
#' variance (d_it^2 (1 - 2 / n_t) + sum_j d_jt^2 / n_t^2) / kappa_T; z_it is
#' its standardized value. Statistics: mean (r / sd)^2, max |r / sd|, mean
#' m_t m_t+1 over consecutive months (m_t the monthly mean residual), mean z^2,
#' mean z_it z_i,t+1 within a cell, mean z_it z_jt over rook-adjacent cells.
#' Observed and replicated statistics share the draw's fitted mean and SDs.
#'
#' @param fit read_fit() of the final model on all-data (spatial loading only).
#' @param input read_inputs().
#' @return list: summary (per statistic: medians and central 95 % intervals,
#'   type-1 quantiles, of observed and replicated values; count and share of
#'   draws with replicated >= observed; its MCSE over the chains; Clopper-
#'   Pearson interval treating draws as independent) and draws (per draw).
column_checks <- function(fit, input) {
  d <- fit$prep$data
  if (d$use_column != 1L || d$col_temporal != 0L)
    stop("The column checks need TROPOMI with the spatial loading only.", call. = FALSE)
  month <- fit$prep$calendar[d$t_col]
  cell_id <- unique(input$column_cells$cell_id)[d$cell_col]
  panel <- .column_panel(data.frame(observation_id = paste0("tropomi:", month, ":", cell_id),
    cell_id = cell_id, month = month, stringsAsFactors = FALSE))
  values <- fit$draws
  S <- nrow(values)
  chain <- rep(seq_len(sampler$chains), each = S / sampler$chains)
  y <- as.numeric(d$y_col)
  d_obs <- as.numeric(d$d_col)
  observed <- replicated <- matrix(0, S, length(column_statistics), dimnames = list(NULL, column_statistics))
  .with_mt_seed(20261109L, {
    for (s in seq_len(S)) {
      p <- .draw(values, s)
      kappa <- exp(p$log_kappa_T)
      mu <- .column_mean(d, p)
      sd <- sqrt(p$sigma_q^2 + d_obs^2 / kappa)
      y_rep <- simulate_data(d, p)$y_col
      observed[s, ] <- .column_values(y, mu, sd, d_obs, kappa, panel)
      replicated[s, ] <- .column_values(y_rep, mu, sd, d_obs, kappa, panel)
    }
  })
  quantiles <- function(x) .quantile_type1(x, .equal_weights(length(x)), c(0.025, 0.5, 0.975))
  summary <- do.call(rbind, lapply(column_statistics, function(k) {
    o <- quantiles(observed[, k])
    r <- quantiles(replicated[, k])
    cbind(data.frame(statistic = k, observed_median = o[2L], observed_lower = o[1L], observed_upper = o[3L],
      replicated_median = r[2L], replicated_lower = r[1L], replicated_upper = r[3L],
      stringsAsFactors = FALSE), .column_tail(observed[, k], replicated[, k], chain))
  }))
  rownames(summary) <- NULL
  list(summary = summary,
    draws = data.frame(draw = seq_len(S), chain = chain, observed = observed, replicated = replicated))
}

#' The posterior predictive checks without the observed features of the
#' low-cost network (its data are not public); their replicated quantiles and
#' p-values stay.
public_ppc <- function(ppc) {
  ppc$observed[startsWith(ppc$source, "Low-cost network")] <- NA_real_
  ppc
}

#' Factor between the monthly and the daily calibration link of the low-cost
#' network. If the daily log BAM values of a month spread with variance s2
#' around their mean and the daily link is log M = a_M + b_M log BAM, the
#' monthly mean of M is exp(a_M) (monthly mean of BAM)^b_M exp{b_M (b_M - 1) s2 / 2}.
#' s2 is the pooled within-month variance of the log daily BAM values of the
#' calibration pairs.
#'
#' @param fit a fit with b_M (the final model on all data).
#' @param input read_inputs().
#' @return one-row data.frame: n_pairs, n_months, s2, b_M posterior mean, the
#'   factor at that mean, and the posterior mean and 95 % interval of the factor.
calibration_factor <- function(fit, input) {
  cal <- input$calibration
  y <- log(cal$bam)
  n_months <- length(unique(cal$month))
  s2 <- sum(tapply(y, cal$month, function(x) sum((x - mean(x))^2))) / (length(y) - n_months)
  b <- fit$draws[, "b_M[1]"]
  factor <- exp(b * (b - 1) * s2 / 2)
  q <- stats::quantile(factor, c(0.025, 0.975), names = FALSE)
  data.frame(model = fit$model, unit = fit$unit, n_pairs = length(y), n_months = n_months,
    var_within_log_bam = s2, b_M_mean = mean(b), factor_at_mean = exp(mean(b) * (mean(b) - 1) * s2 / 2),
    factor_mean = mean(factor), factor_q025 = q[1L], factor_q975 = q[2L], stringsAsFactors = FALSE)
}

# Features of one data set: log values ly, TROPOMI rows y_col, calibration pairs lmoqa_C.
.ppc_features <- function(d, groups, ly, y_col, lmoqa_C) {
  t_row <- rep(seq_len(d$T), d$count)
  basic <- function(v) c(mean_log = mean(v), sd_log = stats::sd(v), min_log = min(v), max_log = max(v),
    mean_conc = mean(exp(v)), max_conc = max(exp(v)))
  with_lag <- c("NO2 background", "PM2.5 reference, background location", "SatPM2.5 municipal mean")
  out <- c()
  for (k in names(groups)) {
    r <- groups[[k]]
    f <- basic(ly[r])
    if (k %in% with_lag) f <- c(f, lag1_log = .lag1(ly[r], t_row[r]))
    if (k == "Low-cost network, all devices") f <- c(f, lag1_log = .lag1_devices(ly[r], t_row[r], d$dev[r]))
    names(f) <- paste0(k, " | ", names(f))
    out <- c(out, f)
  }
  if (d$n_C > 0L) {
    f <- c(mean_log = mean(lmoqa_C), sd_log = stats::sd(lmoqa_C), min_log = min(lmoqa_C), max_log = max(lmoqa_C))
    names(f) <- paste0("Low-cost network, calibration pairs | ", names(f))
    out <- c(out, f)
  }
  if (d$use_column == 1L) {
    monthly <- tapply(y_col, d$t_col, mean)
    f <- c(mean = mean(monthly), sd = stats::sd(monthly), min = min(monthly), max = max(monthly),
      lag1 = .lag1(as.numeric(monthly), as.integer(names(monthly))))
    names(f) <- paste0("TROPOMI monthly mean (1e-5 mol/m2) | ", names(f))
    out <- c(out, f)
  }
  out
}

# Lag-1 correlation over consecutive months of one series.
.lag1 <- function(y, t) {
  o <- order(t)
  y <- y[o]
  t <- t[o]
  i <- which(diff(t) == 1L)
  if (length(i) < 3L) return(NA_real_)
  stats::cor(y[i], y[i + 1L])
}

# Lag-1 correlation over the consecutive months of each device, pooled.
.lag1_devices <- function(y, t, device) {
  pairs <- do.call(rbind, lapply(split(seq_along(y), device), function(r) {
    o <- r[order(t[r])]
    i <- which(diff(t[o]) == 1L)
    cbind(y[o][i], y[o][i + 1L])
  }))
  if (is.null(pairs) || nrow(pairs) < 3L) return(NA_real_)
  stats::cor(pairs[, 1L], pairs[, 2L])
}

# Month groups, consecutive months, within-cell consecutive pairs and same-month
# rook-neighbor pairs of the column rows (cell ids "row:col").
.column_panel <- function(rows) {
  time <- as.integer(substr(rows$month, 1L, 4L)) * 12L + as.integer(substr(rows$month, 6L, 7L))
  support <- data.frame(observation_id = rows$observation_id, group = rows$cell_id, time = time,
    cell = rows$cell_id, stringsAsFactors = FALSE)
  if (anyNA(time) || anyDuplicated(support$observation_id) || anyDuplicated(support[, c("group", "time")]))
    stop("Column rows need unique IDs and one row per cell and month.", call. = FALSE)
  months <- unique(time)
  group <- match(time, months)
  count <- tabulate(group)
  if (any(count < 2L)) stop("A month with one cell has no contrast.", call. = FALSE)
  following <- match(months + 1L, months)
  list(support = support, group = group, count = count,
    month_pairs = cbind(which(!is.na(following)), following[!is.na(following)]),
    lag_pairs = .successive_pairs(support), neighbor_pairs = .neighbor_pairs(support))
}

# Fitted column mean of every row at point p, without q_t (spatial loading only).
.column_mean <- function(d, p) {
  if (!is.null(p$lambda_T)) stop("The column checks need the spatial loading only.", call. = FALSE)
  mu_cell <- p$alpha_T + drop(d$U_cell %*% p$b_cell)
  if (d$field_N) mu_cell <- mu_cell + p$lambda_S * drop(d$B_col %*% p$z_N)
  mu_cell[d$cell_col] + drop(d$H1[d$t_col, , drop = FALSE] %*% p$gamma_T)
}

# The six statistics of column values y at one draw.
.column_values <- function(y, mu, sd, d_obs, kappa, panel) {
  r <- y - mu
  if (any(!is.finite(r)) || any(!is.finite(sd)) || any(sd <= 0) || length(kappa) != 1L ||
      !is.finite(kappa) || kappa <= 0)
    stop("Column statistics need finite residuals and positive SDs.", call. = FALSE)
  m <- as.numeric(rowsum(r, panel$group, reorder = FALSE)) / panel$count
  n <- panel$count[panel$group]
  s2 <- as.numeric(rowsum(d_obs^2, panel$group, reorder = FALSE))[panel$group]
  z <- (r - m[panel$group]) / sqrt((d_obs^2 * (1 - 2 / n) + s2 / n^2) / kappa)
  product <- function(x, pairs) mean(x[pairs[, 1L]] * x[pairs[, 2L]])
  result <- c(dispersion = mean((r / sd)^2), maximum = max(abs(r / sd)),
    monthly = product(m, panel$month_pairs), contrast_dispersion = mean(z^2),
    within_cell = product(z, panel$lag_pairs), neighbor = product(z, panel$neighbor_pairs))
  if (any(!is.finite(result))) stop("A column statistic is not finite.", call. = FALSE)
  result
}

# Share of draws with replicated >= observed (ties count), its MCSE over the
# chains (draws in chain order; not estimable when the indicator is constant)
# and its Clopper-Pearson interval.
.column_tail <- function(observed, replicated, chain) {
  if (length(observed) != length(replicated) || length(chain) != length(observed) ||
      any(!is.finite(c(observed, replicated))))
    stop("The tail summary needs aligned finite draws.", call. = FALSE)
  indicator <- as.numeric(replicated >= observed)
  count <- sum(indicator)
  S <- length(indicator)
  per_chain <- split(indicator, chain)
  if (length(unique(lengths(per_chain))) != 1L) stop("Chains must have equal lengths.", call. = FALSE)
  by_chain <- do.call(cbind, per_chain)
  constant <- count == 0 || count == S
  interval <- stats::binom.test(count, S)$conf.int
  data.frame(count = as.integer(count), draws = S, tail_probability = count / S,
    mcse = if (constant) NA_real_ else stats::sd(indicator) / sqrt(.ess(.split_chains(by_chain))),
    mcse_status = if (constant) "not estimable: constant indicator" else "estimated over chains",
    clopper_pearson_lower = interval[1L], clopper_pearson_upper = interval[2L], stringsAsFactors = FALSE)
}

# Consecutive-time pairs within each group (cell).
.successive_pairs <- function(support) {
  if (nrow(support) < 2L) return(matrix(integer(), 0L, 2L))
  ord <- order(support$group, support$time)
  i <- utils::head(ord, -1L)
  j <- utils::tail(ord, -1L)
  keep <- support$group[i] == support$group[j] & support$time[j] - support$time[i] == 1L
  cbind(i[keep], j[keep])
}

# Same-month pairs of rook-adjacent cells.
.neighbor_pairs <- function(support) {
  if (anyNA(support$cell) || any(!grepl("^-?[0-9]+:-?[0-9]+$", support$cell)))
    stop("Column neighbor pairs need row:col cell ids.", call. = FALSE)
  cells <- unique(support$cell)
  coordinates <- matrix(as.numeric(do.call(rbind, strsplit(cells, ":", fixed = TRUE))), ncol = 2L)
  distance <- abs(outer(coordinates[, 1L], coordinates[, 1L], "-")) +
    abs(outer(coordinates[, 2L], coordinates[, 2L], "-"))
  neighbors <- which(upper.tri(distance) & distance == 1, arr.ind = TRUE)
  keys <- paste(support$time, support$cell, sep = "|")
  pairs <- matrix(integer(), 0L, 2L)
  for (k in seq_len(nrow(neighbors))) {
    i <- which(support$cell == cells[neighbors[k, 1L]])
    j <- match(paste(support$time[i], cells[neighbors[k, 2L]], sep = "|"), keys)
    keep <- !is.na(j)
    pairs <- rbind(pairs, cbind(i[keep], j[keep]))
  }
  pairs
}
