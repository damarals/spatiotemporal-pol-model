# Dense reference of the model's law: every row's mean and the joint covariance
# of the log values (and of the natural TROPOMI values) written directly from
# the law, with every Gaussian effect (temporal components, monthly low-cost
# effect, PM2.5 deviations, device effects, calibration device, TROPOMI's
# common AR(1) error) integrated by hand. TROPOMI enters row by row, one row
# per cell and month, never through the collapsed pseudo-rows, so the
# reference shares no code with the Kalman filter, the collapse,
# joint_moments() or predict_held(). Public stats densities only.

.dense_support <- function(x, y, knots, width_m) {
  w <- exp(-(outer(x, knots[, 1], "-")^2 + outer(y, knots[, 2], "-")^2) / width_m^2)
  w / rowSums(w)
}

.dense_contrasts <- function(K) {
  U <- stats::contr.helmert(K)
  unname(sweep(U, 2L, sqrt(colSums(U^2)), "/"))
}

# Cell means of the standardized NO2 contrast h_S(s) = (w(s) - wbar)' U z_N / r_B,
# averaged point by point over the quadrature of each cell.
dense_h_cell <- function(x, z_N) {
  knots <- x$prep$basis$knots
  width_m <- x$model$width_m
  acag <- x$input$acag_cells
  quadrature <- .column_quadrature(x$input$column_cells)
  U <- .dense_contrasts(nrow(knots))
  W_acag <- .dense_support(acag$x_m, acag$y_m, knots, width_m)
  wbar <- colSums(acag$weight * W_acag)
  B_acag <- sweep(W_acag, 2L, wbar) %*% U
  r_B <- sqrt(sum(acag$weight * rowSums(B_acag^2)))
  W_points <- .dense_support(quadrature$points[, 1], quadrature$points[, 2], knots, width_m)
  h_points <- drop(sweep(W_points, 2L, wbar) %*% U %*% z_N) / r_B
  unname(drop(rowsum(quadrature$weight * h_points, quadrature$cell)))
}

# Training rows in the program's order, read from the Stan data, then one row
# (type 6) per TROPOMI cell and month.
dense_train_rows <- function(x) {
  d <- x$prep$data
  site <- d$site
  at <- pmax(site, 1L)
  w <- d$W_site[at, , drop = FALSE]
  w[site == 0L, ] <- 0
  x_N <- d$Xs_N[at, , drop = FALSE]
  x_P <- d$Xs_P[at, , drop = FALSE]
  x_N[site == 0L, ] <- 0
  x_P[site == 0L, ] <- 0
  rows <- list(t = rep(seq_len(d$T), d$count), type = d$type, cls = d$cls,
    site = ifelse(site > 0L, x$prep$train_sites[at], NA_character_), w = w, x_N = x_N, x_P = x_P,
    device = ifelse(d$type == 4L, x$prep$train_devices[pmax(d$dev, 1L)], ""),
    nugget = d$use_nugget == 1L & site > 0L & d$nug[at] > 0L,
    cell = rep(NA_integer_, d$n), d_sd = rep(NA_real_, d$n))
  if (d$use_column == 1L) {
    m <- d$n_col
    rows <- dense_bind_rows(rows, list(t = as.integer(d$t_col), type = rep(6L, m), cls = rep(0L, m),
      site = rep(NA_character_, m), w = matrix(0, m, d$K), x_N = matrix(0, m, d$Q_N),
      x_P = matrix(0, m, d$Q_P), device = rep("", m), nugget = rep(FALSE, m),
      cell = as.integer(d$cell_col), d_sd = as.numeric(d$d_col)))
  }
  rows
}

dense_train_values <- function(x) c(x$prep$data$ly, as.numeric(x$prep$data$y_col))

# Held rows: support from their coordinates, covariates standardized over every
# location of the table, the deviation at every location with reference PM2.5 rows.
dense_held_rows <- function(x) {
  held <- x$prep$held
  d <- x$prep$data
  table <- .site_covariates(x$input$site_covariates)
  standardized <- function(cols) {
    if (!length(cols)) return(matrix(0, nrow(held), 0L))
    z <- vapply(cols, function(col) {
      v <- as.numeric(table[[col]])
      ((v - mean(v)) / stats::sd(v))[match(held$physical_location_id, table$location_id)]
    }, numeric(nrow(held)))
    matrix(z, nrow(held))
  }
  cov <- x$prep$covariates
  type <- ifelse(held$source == "ground_no2", 1L, ifelse(held$source == "ground_pm25", 2L, 4L))
  train <- dense_train_rows(x)
  bam_sites <- unique(c(train$site[train$type == 2L], held$physical_location_id[type == 2L]))
  list(t = match(held$month, x$prep$calendar), type = type,
    cls = ifelse(type == 1L, match(held$site_class, c("background", "avenue", "bus_terminal")), 0L),
    site = held$physical_location_id,
    w = .dense_support(held$x_m, held$y_m, x$prep$basis$knots, x$model$width_m),
    x_N = standardized(unname(cov[names(cov) == "no2"])),
    x_P = standardized(unname(cov[names(cov) == "pm25"])),
    device = ifelse(type == 4L, held$device_id, ""),
    nugget = d$use_nugget == 1L & held$physical_location_id %in% bam_sites,
    cell = rep(NA_integer_, nrow(held)), d_sd = rep(NA_real_, nrow(held)))
}

dense_bind_rows <- function(a, b) list(t = c(a$t, b$t), type = c(a$type, b$type),
  cls = c(a$cls, b$cls), site = c(a$site, b$site), w = rbind(a$w, b$w), x_N = rbind(a$x_N, b$x_N),
  x_P = rbind(a$x_P, b$x_P), device = c(a$device, b$device), nugget = c(a$nugget, b$nugget),
  cell = c(a$cell, b$cell), d_sd = c(a$d_sd, b$d_sd))

# log r_kp of both profiles: a_p = U zc_p sd_p, r_kp = exp(a_kp) / sum_j wbar_j exp(a_jp).
dense_log_profiles <- function(d, p) {
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
  list(N = lr_N, P = lr_P)
}

# Mean and covariance of the values of `rows` at point p:
#   v_p[t] = Z_p[t] beta_p + Hm[t] h_p + e_p[t],  u_p(s) = x_p(s) theta_p + log sum_k w_k(s) r_kp
#   NO2 background: b_N + u_N + v_N[t];  traffic class c: b_N + u_N + g_c + kappa v_N[t]
#   reference PM2.5: b_P + u_P + v_P[t] + n_P(s);  SatPM2.5: d_A + H1[t] gamma_A + b_P + b_A v_P[t]
#   low-cost device d: a_M + b_M (b_P + u_P + v_P[t] + n_P(s)) + v_d + m_t
#   TROPOMI cell i: alpha_T + H1[t] gamma_T + lambda_S h_cell[i] + beta_cell[i]
#                   + (lambda_T / s_N) e_N[t] + q[t],  beta_cell = U_cell b_cell
#   Cov(e_N[t], e_N[u]) = s_N^2 rho_N^|t-u|, Cov(e_P[t], e_P[u]) = s_P^2 rho_P^|t-u|,
#   Cov(e_N[t], e_P[u]) = c0 rho_N^(t-u) for t >= u and c0 rho_P^(u-t) for u > t,
#   c0 = r s_N s_P sqrt((1 - rho_N^2)(1 - rho_P^2)) / (1 - rho_N rho_P);
#   Cov(q[t], q[u]) = sigma_q^2 rho_q^|t-u|;
#   n_P(s) ~ N(0, tau_P^2), v_d ~ N(0, tau_M^2), m_t ~ N(0, tau_m^2);
#   noise sig per log source, d_it^2 / kappa_T per TROPOMI row.
dense_law <- function(x, rows, p) {
  d <- x$prep$data
  f_N <- drop(d$Z_N %*% p$beta_N) + drop(d$Hm %*% p$h_N)
  f_P <- drop(d$Z_P %*% p$beta_P) + drop(d$Hm %*% p$h_P)
  lr <- dense_log_profiles(d, p)
  t <- rows$t
  type <- rows$type
  n <- length(t)
  u_N <- log(drop(rows$w %*% exp(lr$N))) + drop(rows$x_N %*% p$theta_N)
  u_P <- log(drop(rows$w %*% exp(lr$P))) + drop(rows$x_P %*% p$theta_P)
  u_N[type %in% c(3L, 6L)] <- 0
  u_P[type %in% c(3L, 6L)] <- 0
  traffic <- type == 1L & rows$cls > 1L
  offset <- c(0, p$g[1L], p$g[d$n_g])
  kappa <- if (d$use_kappa) p$kappa else 0
  amp <- ifelse(traffic, kappa, 1)
  bA <- if (d$free_bA) p$b_A else 1
  seasonal <- if (d$acag_seasonal && !is.null(p$gamma_A)) drop(d$H1 %*% p$gamma_A) else rep(0, d$T)
  scalar <- function(v) if (is.null(v)) NA_real_ else v
  column <- type == 6L
  mean_col <- rep(0, n)
  if (any(column)) {
    h_cell <- if (d$field_N) p$lambda_S * dense_h_cell(x, p$z_N) else rep(0, d$n_cells)
    beta_cell <- drop(.dense_contrasts(d$n_cells) %*% p$b_cell)
    mean_col[column] <- p$alpha_T + drop(d$H1[t[column], , drop = FALSE] %*% p$gamma_T) +
      h_cell[rows$cell[column]] + beta_cell[rows$cell[column]]
  }
  mean <- ifelse(type == 1L, p$b_N + u_N + offset[pmax(rows$cls, 1L)] + amp * f_N[t],
    ifelse(type == 2L, p$b_P + u_P + f_P[t],
      ifelse(type == 3L, p$d_A + seasonal[t] + p$b_P + bA * f_P[t],
        ifelse(type == 6L, mean_col, scalar(p$a_M) + scalar(p$b_M) * (p$b_P + u_P + f_P[t])))))
  lambda_T <- if (is.null(p$lambda_T)) 0 else p$lambda_T
  l_N <- ifelse(type == 1L, amp, ifelse(type == 6L, lambda_T / p$s_N, 0))
  l_P <- ifelse(type == 2L, 1, ifelse(type == 3L, bA, ifelse(type == 4L, scalar(p$b_M), 0)))
  lag <- outer(t, t, "-")
  V <- p$s_N^2 * p$rho_N^abs(lag) * tcrossprod(l_N) + p$s_P^2 * p$rho_P^abs(lag) * tcrossprod(l_P)
  if (d$coupled) {
    c0 <- p$r * p$s_N * p$s_P * sqrt((1 - p$rho_N^2) * (1 - p$rho_P^2)) / (1 - p$rho_N * p$rho_P)
    C <- c0 * ifelse(lag >= 0, p$rho_N^pmax(lag, 0), p$rho_P^pmax(-lag, 0))
    cross <- outer(l_N, l_P) * C
    V <- V + cross + t(cross)
  }
  same_site <- outer(rows$site, rows$site, "==")
  same_site[is.na(same_site)] <- FALSE
  if (d$use_nugget) {
    l_nug <- ifelse(rows$nugget & type == 2L, 1, ifelse(rows$nugget & type == 4L, scalar(p$b_M), 0))
    V <- V + p$tau_P^2 * tcrossprod(l_nug) * same_site
  }
  moqa <- type == 4L
  if (any(moqa)) {
    V <- V + p$tau_M^2 * (outer(rows$device, rows$device, "==") & outer(moqa, moqa)) +
      p$tau_m^2 * (outer(t, t, "==") & outer(moqa, moqa))
  }
  if (any(column)) V <- V + p$sigma_q^2 * p$rho_q^abs(lag) * outer(column, column)
  noise2 <- ifelse(type == 1L, ifelse(traffic, p$sig_NT, p$sig_N)^2,
    ifelse(type == 2L, p$sig_P^2, ifelse(type == 3L, scalar(p$sig_A)^2,
      ifelse(type == 6L, rows$d_sd^2 / exp(scalar(p$log_kappa_T)), scalar(p$sig_M)^2))))
  V <- V + diag(noise2, n)
  list(mean = mean, V = V)
}

# Quartic scale a with Pr(sigma > U) = 0.05: (sigma / a)^4 ~ Gamma(1/4, 1).
.dense_quartic_a <- function(U) U / stats::qgamma(0.05, shape = 1 / 4, lower.tail = FALSE)^(1 / 4)

.dense_scale_U <- c(s_N = 0.6, s_P = 0.6, tau_P = 0.6, tau_M = 0.6, tau_m = 0.6, sig_N = 1, sig_NT = 1,
  sig_P = 1, sig_A = 1, sig_M = 1, sig_C = 1, sd_N = 1, sd_P = 0.5, sigma_q = 1, sigma_beta = 1)
.dense_half_normal_sd <- c(s_N = 0.3, s_P = 0.3, tau_P = 0.3, tau_M = 0.5, tau_m = 0.3, sig_N = 0.5,
  sig_NT = 0.5, sig_P = 0.5, sig_A = 0.5, sig_M = 0.5, sig_C = 0.5, sd_N = 0.449258049,
  sd_P = 0.449258049, sigma_q = 0.5, sigma_beta = 0.5)

dense_scale_log_prior <- function(d, p) {
  a <- if (d$quartic) .dense_quartic_a(.dense_scale_U) else .dense_half_normal_sd
  present <- c("s_N", "s_P", "sig_N", "sig_NT", "sig_P", if (!is.null(p$sig_A)) "sig_A",
    if (!is.null(p$sd_N)) "sd_N", if (!is.null(p$sd_P)) "sd_P", if (!is.null(p$tau_P)) "tau_P",
    if (!is.null(p$tau_M)) c("tau_M", "sig_M", "tau_m"), if (!is.null(p$sig_C)) "sig_C",
    if (!is.null(p$sigma_q)) c("sigma_q", "sigma_beta"))
  s <- unlist(p[present])
  a <- a[present]
  if (d$quartic) return(sum(log(4) - log(a) - lgamma(1 / 4) - (s / a)^4))
  sum(log(2) + stats::dnorm(s, 0, a, log = TRUE))
}

# Normal prior terms with a fixed SD (value, mean, SD).
dense_prior_terms <- function(p) {
  term <- function(x, m, s) if (length(x)) cbind(x = x, m = rep(m, length(x)), s = rep(s, length(x)))
  rbind(term(c(p$b_N, p$b_P), log(10), 1), term(p$g, log(1.8), 0.35),
    term(c(p$beta_N, p$beta_P, p$h_N, p$h_P, p$theta_N, p$theta_P), 0, 0.5),
    term(p$d_A, 0, 0.3), term(p$gamma_A, 0, 0.3), term(p$b_A, 1, 0.5),
    term(c(p$z_N, p$z_P), 0, 1), term(p$a_M, 0, 1), term(p$b_M, 1, 0.3),
    term(p$alpha_T, 0, 3), term(p$gamma_T, 0, 0.3), term(c(p$lambda_S, p$lambda_T), 0, 0.5),
    term(p$log_kappa_T, 0, 1))
}

dense_dmvn <- function(x, mu, S) {
  R <- chol(S)
  z <- forwardsolve(t(R), x - mu)
  -0.5 * (length(x) * log(2 * pi) + sum(z^2)) - sum(log(diag(R)))
}

# Complete log density of the training data and the priors (every constant).
dense_log_density <- function(x, p) {
  d <- x$prep$data
  law <- dense_law(x, dense_train_rows(x), p)
  terms <- dense_prior_terms(p)
  value <- sum(stats::dnorm(terms[, "x"], terms[, "m"], terms[, "s"], log = TRUE)) +
    dense_scale_log_prior(d, p) + dense_dmvn(dense_train_values(x), law$mean, law$V)
  for (rho in c(p$r, p$rho_S, p$rho_q)) value <- value + stats::dbeta((rho + 1) / 2, 2, 2, log = TRUE)
  if (!is.null(p$b_cell)) value <- value + sum(stats::dnorm(p$b_cell, 0, p$sigma_beta, log = TRUE))
  if (!is.null(p$sig_C)) {
    value <- value + dense_dmvn(d$lmoqa_C, p$a_M + p$b_M * d$lbam_C,
      diag(p$sig_C^2, d$n_C) + matrix(p$tau_M^2, d$n_C, d$n_C))
  }
  value
}

# The constants that `~` drops from the Normal priors.
dense_dropped_constant <- function(p) {
  sum(log(dense_prior_terms(p)[, "s"]) + log(2 * pi) / 2) + length(p$b_cell) * log(2 * pi) / 2
}

# A random point in the interior of the parameter space of data d.
dense_random_point <- function(d) {
  u <- function(k, lo, hi) stats::runif(k, lo, hi)
  r <- function(k, m, s) stats::rnorm(k, m, s)
  has_a <- any(d$type == 3L)
  has_m <- any(d$type == 4L)
  use_cal <- has_m && d$n_C > 0L
  coupled_space <- d$couple_space == 1L && d$field_N == 1L && d$field_P == 1L
  col <- d$use_column == 1L
  p <- list(b_N = r(1, log(10), 0.3), b_P = r(1, log(10), 0.3), g = r(d$n_g, 0.5, 0.2),
    beta_N = r(d$P_N, 0, 0.3), beta_P = r(d$P_P, 0, 0.3), h_N = r(d$H, 0, 0.3), h_P = r(d$H, 0, 0.3),
    theta_N = r(d$Q_N, 0, 0.3), theta_P = r(d$Q_P, 0, 0.3),
    d_A = if (has_a) r(1, 0, 0.2), gamma_A = if (has_a && d$acag_seasonal) r(2, 0, 0.2),
    b_A = if (has_a && d$free_bA) u(1, 0.6, 1.8),
    s_N = u(1, 0.05, 0.5), s_P = u(1, 0.05, 0.5), rho_N = u(1, 0.1, 0.95), rho_P = u(1, 0.1, 0.95),
    r = if (d$coupled) u(1, -0.8, 0.8), kappa = if (d$use_kappa) u(1, 0.05, 0.95),
    sig_N = u(1, 0.05, 0.5), sig_NT = u(1, 0.05, 0.5), sig_P = u(1, 0.05, 0.5),
    sig_A = if (has_a) u(1, 0.05, 0.5),
    z_N = if (d$field_N) r(d$K - 1L, 0, 1), z_P = if (d$field_P) r(d$K - 1L, 0, 1),
    sd_N = if (d$field_N) u(1, 0.1, 0.8), sd_P = if (d$field_P) u(1, 0.1, 0.8),
    rho_S = if (coupled_space) u(1, -0.8, 0.8), tau_P = if (d$use_nugget) u(1, 0.05, 0.4),
    a_M = if (has_m) r(1, 0, 0.3), b_M = if (has_m) u(1, 0.7, 1.5), tau_M = if (has_m) u(1, 0.05, 0.4),
    sig_M = if (has_m) u(1, 0.05, 0.5), sig_C = if (use_cal) u(1, 0.05, 0.5),
    tau_m = if (has_m) u(1, 0.05, 0.3),
    alpha_T = if (col) r(1, 1.3, 0.3), gamma_T = if (col) r(2, 0, 0.2),
    lambda_S = if (col && d$field_N) { if (d$sign_constraint) u(1, 0.05, 0.8) else u(1, -0.8, 0.8) },
    lambda_T = if (col && d$col_temporal) u(1, -0.8, 0.8),
    log_kappa_T = if (col) r(1, 2, 0.5), sigma_q = if (col) u(1, 0.05, 0.4),
    rho_q = if (col) u(1, -0.8, 0.8), sigma_beta = if (col) u(1, 0.05, 0.4),
    b_cell = if (col) r(d$n_cells - 1L, 0, 0.2))
  p[!vapply(p, is.null, logical(1))]
}

# Parameters the program declares as plain reals; the rest are vectors or
# size-one arrays, named name[i].
.dense_scalars <- c("b_N", "b_P", "s_N", "s_P", "rho_N", "rho_P", "sig_N", "sig_NT", "sig_P")

# The point as one row of draws, with the program's variable names.
dense_values <- function(p) {
  v <- unlist(lapply(names(p), function(k) {
    x <- p[[k]]
    if (!length(x)) return(NULL)
    stats::setNames(x, if (k %in% .dense_scalars) k else paste0(k, "[", seq_along(x), "]"))
  }))
  matrix(v, 1L, dimnames = list(NULL, names(v)))
}

# The program's log density at constrained point p (CmdStan's log_prob method),
# without the Jacobian of the constraining transforms.
stan_log_density <- function(d, p, exe = test_program()) {
  dir <- tempfile("fusion-log-prob-")
  dir.create(dir)
  on.exit(unlink(dir, recursive = TRUE), add = TRUE)
  data_file <- file.path(dir, "data.json")
  point_file <- file.path(dir, "point.json")
  output_file <- file.path(dir, "lp.csv")
  write_stan_data(d, data_file)
  json <- lapply(names(p), function(k) if (k %in% .dense_scalars) jsonlite::unbox(p[[k]]) else p[[k]])
  names(json) <- names(p)
  writeLines(jsonlite::toJSON(json, digits = NA, auto_unbox = FALSE), point_file)
  status <- system2(exe, c("method=log_prob", paste0("constrained_params=", point_file), "jacobian=0",
    "data", paste0("file=", data_file), "output", paste0("file=", output_file), "sig_figs=18"),
    stdout = FALSE, stderr = FALSE)
  if (status != 0L) stop("CmdStan log_prob failed.", call. = FALSE)
  utils::read.csv(output_file, comment.char = "#")[1L, 1L]
}
