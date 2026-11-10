# Simulation-based calibration (Talts et al. 2018) of a model on the design of a
# fit (every row, month, location, device and TROPOMI cell it trains on), and the
# forward simulator of the joint law that it and the predictive checks share.
# The simulator restates the law of stan/fusion.stan from the Stan data and a
# parameter point, drawing every quantity the program integrates: the bivariate
# AR(1) residuals, TROPOMI's AR(1) common error q_t, the monthly and device
# effects of the low-cost network, the PM2.5 deviations, the calibration device
# effect and every row's noise. It shares no code with joint_moments().

# Seed of replication `rep` is sbc_seed + rep.
sbc_seed <- 20261108L

# Replications whose sampling was clean: no divergence, R-hat <= 1.01, bulk ESS >= 400.
sbc_clean <- list(max_rhat = 1.01, min_ess_bulk = 400)

#' One draw of every parameter from the prior the program declares for data d:
#' Normal fixed effects, half exponential-power scales ((sigma / a)^4 ~
#' Gamma(1/4, 1); half-Normal when quartic = 0), Uniform(0, 1) persistences and
#' kappa, (x + 1) / 2 ~ Beta(2, 2) correlations, truncated Normals for b_A and
#' b_M, lambda_S half-Normal under the sign constraint, no lambda_T without the
#' temporal loading, and b_cell ~ N(0, sigma_beta^2).
#'
#' @param d Stan data (prepare_fit()$data).
#' @return list of parameter blocks, as .draw() names them.
draw_prior <- function(d) {
  a <- d$prior_scale
  scale <- function(i) if (d$quartic) a[i] * stats::rgamma(1L, shape = 1 / 4, rate = 1)^(1 / 4) else
    abs(stats::rnorm(1L, 0, a[i]))
  truncated <- function(m, s, lo, hi) {
    repeat {
      x <- stats::rnorm(1L, m, s)
      if (x > lo && x < hi) return(x)
    }
  }
  correlation <- function() 2 * stats::rbeta(1L, 2, 2) - 1
  has_a <- any(d$type == 3L)
  has_m <- any(d$type == 4L)
  use_cal <- has_m && d$n_C > 0L
  col <- d$use_column == 1L
  coupled_space <- d$couple_space == 1L && d$field_N == 1L && d$field_P == 1L
  p <- list(b_N = stats::rnorm(1L, log(10), 1), b_P = stats::rnorm(1L, log(10), 1),
    g = stats::rnorm(d$n_g, log(1.8), 0.35), beta_N = stats::rnorm(d$P_N, 0, 0.5),
    beta_P = stats::rnorm(d$P_P, 0, 0.5), h_N = stats::rnorm(d$H, 0, 0.5), h_P = stats::rnorm(d$H, 0, 0.5),
    theta_N = stats::rnorm(d$Q_N, 0, 0.5), theta_P = stats::rnorm(d$Q_P, 0, 0.5),
    d_A = if (has_a) stats::rnorm(1L, 0, 0.3),
    gamma_A = if (has_a && d$acag_seasonal) stats::rnorm(2L, 0, 0.3),
    b_A = if (has_a && d$free_bA) truncated(1, 0.5, 0.5, 2.5),
    s_N = scale(1L), s_P = scale(2L), rho_N = stats::runif(1L), rho_P = stats::runif(1L),
    r = if (d$coupled) correlation(), kappa = if (d$use_kappa) stats::runif(1L),
    sig_N = scale(6L), sig_NT = scale(7L), sig_P = scale(8L), sig_A = if (has_a) scale(9L),
    z_N = if (d$field_N) stats::rnorm(d$K - 1L), z_P = if (d$field_P) stats::rnorm(d$K - 1L),
    sd_N = if (d$field_N) scale(12L), sd_P = if (d$field_P) scale(13L),
    rho_S = if (coupled_space) correlation(), tau_P = if (d$use_nugget) scale(3L),
    a_M = if (has_m) stats::rnorm(1L, 0, 1), b_M = if (has_m) truncated(1, 0.3, 0.5, 2),
    tau_M = if (has_m) scale(4L), sig_M = if (has_m) scale(10L), sig_C = if (use_cal) scale(11L),
    tau_m = if (has_m) scale(5L),
    alpha_T = if (col) stats::rnorm(1L, 0, 3), gamma_T = if (col) stats::rnorm(2L, 0, 0.3),
    lambda_S = if (col && d$field_N) {
      x <- stats::rnorm(1L, 0, 0.5)
      if (d$sign_constraint == 1L) abs(x) else x
    },
    lambda_T = if (col && d$col_temporal == 1L) stats::rnorm(1L, 0, 0.5),
    log_kappa_T = if (col) stats::rnorm(1L, 0, 1), sigma_q = if (col) scale(14L),
    rho_q = if (col) correlation(), sigma_beta = if (col) scale(15L))
  if (col) p$b_cell <- stats::rnorm(d$n_cells - 1L, 0, p$sigma_beta)
  p[!vapply(p, is.null, logical(1))]
}

#' Simulated observations of the law at point p on the design d, with every
#' integrated quantity drawn afresh.
#'
#' @param d Stan data.
#' @param p parameter point (draw_prior() or .draw()).
#' @return list: ly (log values of the reference, SatPM2.5 and low-cost rows),
#'   lmoqa_C (log low-cost values of the calibration pairs), y_col (TROPOMI rows).
simulate_data <- function(d, p) {
  T <- d$T
  t_row <- rep(seq_len(T), d$count)
  has_m <- any(d$type == 4L)
  use_cal <- has_m && d$n_C > 0L
  col <- d$use_column == 1L
  f_N <- drop(d$Z_N %*% p$beta_N) + drop(d$Hm %*% p$h_N)
  f_P <- drop(d$Z_P %*% p$beta_P) + drop(d$Hm %*% p$h_P)
  profile <- function(z, sd) {
    a <- drop(d$U %*% z) * sd
    a - log(sum(d$wbar * exp(a)))
  }
  flat <- rep(-log(sum(d$wbar)), d$K)
  lr_N <- if (d$field_N) profile(p$z_N, p$sd_N) else flat
  lr_P <- flat
  if (d$field_P) {
    zc <- if (d$couple_space && d$field_N) p$rho_S * p$z_N + sqrt(1 - p$rho_S^2) * p$z_P else p$z_P
    lr_P <- profile(zc, p$sd_P)
  }
  u_N <- log(drop(d$W_site %*% exp(lr_N))) + (if (d$Q_N > 0L) drop(d$Xs_N %*% p$theta_N) else 0)
  u_P <- log(drop(d$W_site %*% exp(lr_P))) + (if (d$Q_P > 0L) drop(d$Xs_P %*% p$theta_P) else 0)
  # Bivariate AR(1) residuals with correlated innovations, stationary start.
  rr <- if (d$coupled) p$r else 0
  iN <- p$s_N * sqrt(1 - p$rho_N^2)
  iP <- p$s_P * sqrt(1 - p$rho_P^2)
  c0 <- rr * iN * iP / (1 - p$rho_N * p$rho_P)
  L0 <- t(chol(matrix(c(p$s_N^2, c0, c0, p$s_P^2), 2L)))
  Lw <- t(chol(matrix(c(iN^2, rr * iN * iP, rr * iN * iP, iP^2), 2L)))
  e <- matrix(0, T, 2L)
  e[1L, ] <- drop(L0 %*% stats::rnorm(2L))
  for (t in seq_len(T)[-1L])
    e[t, ] <- c(p$rho_N * e[t - 1L, 1L], p$rho_P * e[t - 1L, 2L]) + drop(Lw %*% stats::rnorm(2L))
  m_t <- if (has_m) stats::rnorm(T, 0, p$tau_m) else numeric(T)
  v_dev <- if (has_m && d$D > 0L) stats::rnorm(d$D, 0, p$tau_M) else numeric()
  n_nug <- if (d$use_nugget && d$Jn > 0L) stats::rnorm(d$Jn, 0, p$tau_P) else numeric()
  nugget_of <- function(site) ifelse(site > 0L & d$nug[pmax(site, 1L)] > 0L,
    n_nug[pmax(d$nug[pmax(site, 1L)], 1L)], 0)
  k <- if (d$use_kappa) p$kappa else 0
  bA <- if (d$free_bA && !is.null(p$b_A)) p$b_A else 1
  g_of <- c(0L, 1L, d$n_g)
  type <- d$type
  cls <- d$cls
  site <- d$site
  mean <- noise <- numeric(d$n)
  i1 <- type == 1L
  traffic <- i1 & cls > 1L
  amp <- ifelse(traffic, k, 1)
  shift <- ifelse(traffic, p$g[pmax(g_of[pmax(cls, 1L)], 1L)], 0)
  mean[i1] <- (p$b_N + u_N[pmax(site, 1L)] + shift + amp * (f_N[t_row] + e[t_row, 1L]))[i1]
  noise[i1] <- ifelse(traffic, p$sig_NT, p$sig_N)[i1]
  i2 <- type == 2L
  mean[i2] <- (p$b_P + u_P[pmax(site, 1L)] + f_P[t_row] + e[t_row, 2L] + nugget_of(site))[i2]
  noise[i2] <- p$sig_P
  i3 <- type == 3L
  if (any(i3)) {
    seasonal <- if (d$acag_seasonal) drop(d$H1 %*% p$gamma_A) else numeric(T)
    mean[i3] <- (p$d_A + seasonal[t_row] + p$b_P + bA * (f_P[t_row] + e[t_row, 2L]))[i3]
    noise[i3] <- p$sig_A
  }
  i4 <- type == 4L
  if (any(i4)) {
    mean[i4] <- (p$a_M + p$b_M * (p$b_P + u_P[pmax(site, 1L)] + f_P[t_row] + e[t_row, 2L] + nugget_of(site)) +
      v_dev[pmax(d$dev, 1L)] + m_t[t_row])[i4]
    noise[i4] <- p$sig_M
  }
  ly <- mean + stats::rnorm(d$n, 0, noise)
  lmoqa_C <- if (use_cal) p$a_M + stats::rnorm(1L, 0, p$tau_M) + p$b_M * d$lbam_C +
    stats::rnorm(d$n_C, 0, p$sig_C) else numeric()
  y_col <- numeric()
  if (col) {
    q <- .ar1_path(T, p$rho_q, p$sigma_q)
    mu_cell <- p$alpha_T + drop(d$U_cell %*% p$b_cell)
    if (d$field_N) mu_cell <- mu_cell + p$lambda_S * drop(d$B_col %*% p$z_N)
    lamT <- if (is.null(p$lambda_T)) 0 else p$lambda_T / p$s_N
    mean_col <- mu_cell[d$cell_col] + drop(d$H1[d$t_col, , drop = FALSE] %*% p$gamma_T) +
      lamT * e[d$t_col, 1L] + q[d$t_col]
    y_col <- mean_col + stats::rnorm(d$n_col, 0, d$d_col / sqrt(exp(p$log_kappa_T)))
  }
  list(ly = ly, lmoqa_C = lmoqa_C, y_col = y_col)
}

#' One SBC replication: a prior draw, data simulated from it on the design of
#' prep, a fit of the program to those data (sampler settings of fit_model())
#' and the rank of each true scalar among 99 thinned posterior draws.
#'
#' @param rep replication number (seed sbc_seed + rep).
#' @param prep prepare_fit() of the design (the final model on all-data).
#' @param dir scratch directory for the CmdStan files; removed at the end.
#' @param exe compiled program (compile_program()).
#' @return list: rep, seed, truth (named), rank (named, 0..99), diagnostics.
sbc_replicate <- function(rep, prep, dir, exe = compile_program()) {
  d <- prep$data
  seed <- sbc_seed + as.integer(rep)
  sim <- .with_seed(seed, {
    p <- draw_prior(d)
    list(truth = .flatten(p), data = .replace_observations(d, simulate_data(d, p)))
  })
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(dir, recursive = TRUE), add = TRUE)
  data_file <- file.path(dir, "data.json")
  write_stan_data(sim$data, data_file)
  started <- Sys.time()
  .run_chains(exe, data_file, dir, seed)
  runtime <- as.numeric(difftime(Sys.time(), started, units = "secs"))
  chains <- .read_chains(dir)
  diagnostics <- .diagnostics(chains)
  diagnostics$runtime_s <- runtime
  list(rep = as.integer(rep), seed = seed, truth = sim$truth, rank = .sbc_rank(chains$draws, sim$truth),
    diagnostics = diagnostics)
}

#' Replications as one table.
#'
#' @param reps list of sbc_replicate() results.
#' @return data.frame: rep, seed, divergences, max_rhat, min_ess_bulk,
#'   runtime_s, treedepth_hits, min_ebfmi, truth_<parameter>, rank_<parameter>.
sbc_table <- function(reps) {
  params <- names(reps[[1L]]$truth)
  out <- do.call(rbind, lapply(reps, function(r) {
    if (!identical(names(r$truth), params) || !identical(names(r$rank), params))
      stop("Replication ", r$rep, " has other parameters.", call. = FALSE)
    dg <- r$diagnostics
    data.frame(rep = r$rep, seed = r$seed, divergences = dg$divergences, max_rhat = dg$max_rhat,
      min_ess_bulk = dg$min_ess_bulk, runtime_s = dg$runtime_s, treedepth_hits = dg$treedepth_hits,
      min_ebfmi = dg$min_ebfmi, t(stats::setNames(r$truth, paste0("truth_", params))),
      t(stats::setNames(r$rank, paste0("rank_", params))), check.names = FALSE, stringsAsFactors = FALSE)
  }))
  rownames(out) <- NULL
  out
}

#' Rank uniformity and coverage per parameter, over all replications and over
#' the cleanly sampled ones (an exploratory subset, selected on the prior draw).
#'
#' @param table sbc_table().
#' @return data.frame per parameter: n, chisq and p_value (chi-square on ten
#'   equal bins of rank %/% 10), p_holm (Holm over parameters), mean_u and sd_u
#'   of (rank + 0.5) / 100 (uniform: 0.5, 0.289), min_bin, max_bin,
#'   coverage90 (share of ranks 5..94: the central 90 % interval covers the
#'   truth), and the same with suffix _clean over the clean replications.
summarise_sbc <- function(table) {
  params <- sub("^rank_", "", grep("^rank_", names(table), value = TRUE))
  clean <- table$divergences == 0 & table$max_rhat <= sbc_clean$max_rhat &
    table$min_ess_bulk >= sbc_clean$min_ess_bulk
  one <- function(keep) {
    out <- do.call(rbind, lapply(params, function(k) {
      rank <- table[[paste0("rank_", k)]][keep]
      counts <- .rank_bins(rank)
      test <- suppressWarnings(stats::chisq.test(counts, p = rep(0.1, 10L)))
      u <- (rank + 0.5) / 100
      data.frame(n = length(rank), chisq = unname(test$statistic), p_value = test$p.value,
        mean_u = mean(u), sd_u = stats::sd(u), min_bin = min(counts), max_bin = max(counts),
        coverage90 = mean(rank >= 5 & rank <= 94))
    }))
    out$p_holm <- stats::p.adjust(out$p_value, "holm")
    out[c("n", "chisq", "p_value", "p_holm", "mean_u", "sd_u", "min_bin", "max_bin", "coverage90")]
  }
  all <- one(rep(TRUE, nrow(table)))
  subset <- one(clean)
  names(subset) <- paste0(names(subset), "_clean")
  data.frame(parameter = params, all, subset, row.names = NULL, stringsAsFactors = FALSE)
}

#' Rank counts per parameter in the ten equal bins of rank %/% 10, over all
#' replications (expected count under uniformity: replications / 10).
sbc_bins <- function(table) {
  params <- sub("^rank_", "", grep("^rank_", names(table), value = TRUE))
  counts <- t(vapply(params, function(k) .rank_bins(table[[paste0("rank_", k)]]), integer(10L)))
  colnames(counts) <- paste0("bin_", 1:10)
  data.frame(parameter = params, counts, expected = nrow(table) / 10, row.names = NULL,
    stringsAsFactors = FALSE)
}

#' Counts the paper quotes about the calibration: replications, replications
#' with at least one divergent transition, parameters, parameters with an
#' unadjusted p-value below 0.05 and the number expected under uniformity
#' (parameters x 0.05), expected count per bin and the Monte Carlo SE of a 90 %
#' coverage over the replications.
sbc_overview <- function(table, summary) {
  n <- nrow(table)
  data.frame(replications = n, with_divergences = sum(table$divergences > 0),
    parameters = nrow(summary), p_below_005 = sum(summary$p_value < 0.05),
    expected_p_below_005 = nrow(summary) * 0.05, expected_per_bin = n / 10,
    coverage90_se = sqrt(0.9 * 0.1 / n), stringsAsFactors = FALSE)
}

.rank_bins <- function(rank) tabulate(rank %/% 10 + 1, nbins = 10L)

# Parameters the program declares as plain reals; every other block is a vector
# or a size-one array, named name[i] in the draws.
.sbc_scalars <- c("b_N", "b_P", "s_N", "s_P", "rho_N", "rho_P", "sig_N", "sig_NT", "sig_P")

# A parameter point flattened with the names CmdStan gives the draws.
.flatten <- function(p) {
  unlist(lapply(names(p), function(k) {
    x <- p[[k]]
    if (!length(x)) return(NULL)
    stats::setNames(x, if (k %in% .sbc_scalars) k else paste0(k, "[", seq_along(x), "]"))
  }))
}

# The data d with its observations replaced by a simulation.
.replace_observations <- function(d, sim) {
  d$ly <- sim$ly
  if (d$n_C > 0L) d$lmoqa_C <- .array(sim$lmoqa_C)
  if (d$use_column == 1L) d$y_col <- .array(sim$y_col)
  d
}

# Rank of each true value among 99 draws thinned evenly from the S draws (every
# 40th of 4000, the last dropped): the number of thinned draws below it, 0..99.
.sbc_rank <- function(draws, truth) {
  S <- nrow(draws)
  keep <- round(seq(S / 100, S - S / 100, length.out = 99L))
  thinned <- draws[keep, names(truth), drop = FALSE]
  vapply(names(truth), function(k) sum(thinned[, k] < truth[[k]]), numeric(1))
}

# A stationary AR(1) path of length T with marginal SD s and persistence rho.
.ar1_path <- function(T, rho, s) {
  e <- numeric(T)
  e[1L] <- stats::rnorm(1L, 0, s)
  innovation <- s * sqrt(1 - rho^2)
  for (t in seq_len(T)[-1L]) e[t] <- rho * e[t - 1L] + stats::rnorm(1L, 0, innovation)
  e
}
