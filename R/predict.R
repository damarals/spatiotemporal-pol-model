# Predictive draws and log densities at the held-out targets. In each posterior
# draw the log values of the training and held rows are jointly Normal, with
# every Gaussian effect integrated (temporal components, monthly low-cost
# effect, PM2.5 deviations, device effects, TROPOMI's common error); the held
# values are drawn from their conditional law given the training values, and
# each target's log density is lognormal on the concentration scale.

#' Predict the held targets of a fit.
#'
#' @param fit fit_model() or read_fit().
#' @return list: draws and log_density, targets x draws matrices (row names
#'   the targets' observation ids).
predict_held <- function(fit) {
  prep <- fit$prep
  values <- fit$draws
  held <- prep$held
  S <- nrow(values)
  o <- seq_len(prep$n_train)
  h <- prep$n_train + seq_len(nrow(held))
  rows <- prep$rows[o, , drop = FALSE]
  y <- ifelse(rows$type == 5L, rows$observed, log(rows$observed))
  y_held <- log(held$observed)
  draws <- log_density <- matrix(NA_real_, nrow(held), S, dimnames = list(held$observation_id, NULL))
  .with_seed(fit$seed + 1L, {
    for (j in seq_len(if (nrow(held)) S else 0L)) {
      m <- joint_moments(.draw(values, j), prep)
      R <- chol(m$V[o, o])
      alpha <- backsolve(R, forwardsolve(t(R), y - m$mean[o]))
      W <- forwardsolve(t(R), m$V[o, h, drop = FALSE])
      mu <- m$mean[h] + drop(crossprod(m$V[o, h, drop = FALSE], alpha))
      C <- m$V[h, h, drop = FALSE] - crossprod(W)
      sd <- sqrt(pmax(diag(C), 1e-12))
      L <- tryCatch(t(chol(C + diag(1e-10, nrow(C)))), error = function(e)
        stop("Conditional covariance of the held rows is not positive definite in draw ", j, ".",
          call. = FALSE))
      draws[, j] <- exp(mu + drop(L %*% stats::rnorm(length(h))))
      log_density[, j] <- stats::dnorm(y_held, mu, sd, log = TRUE) - y_held
    }
  })
  if (any(!is.finite(draws)) || any(draws <= 0) || any(!is.finite(log_density)))
    stop("Non-finite or non-positive predictive output.", call. = FALSE)
  list(draws = draws, log_density = log_density)
}

#' Mean and covariance of the joint Normal law of `rows` (default: training then
#' held rows of the fit) at one posterior draw `par` (a .draw() list). Rows with
#' latent = TRUE carry no observation noise.
joint_moments <- function(par, prep, rows = prep$rows) {
  d <- prep$data
  n <- nrow(rows)
  t <- rows$t
  type <- rows$type
  f_N <- drop(d$Z_N %*% par$beta_N) + drop(d$Hm %*% par$h_N)
  f_P <- drop(d$Z_P %*% par$beta_P) + drop(d$Hm %*% par$h_P)
  lr_N <- if (d$field_N) .log_profile(prep$basis, par$z_N, par$sd_N) else rep(0, prep$basis$K)
  lr_P <- rep(0, prep$basis$K)
  if (d$field_P) {
    zc <- if (d$couple_space) par$rho_S * par$z_N + sqrt(1 - par$rho_S^2) * par$z_P else par$z_P
    lr_P <- .log_profile(prep$basis, zc, par$sd_P)
  }
  at <- match(rows$site, prep$sites)
  has_site <- !is.na(at)
  first <- ifelse(has_site, at, 1L)
  W <- prep$W_all[first, , drop = FALSE]
  u_N <- ifelse(has_site, log(drop(W %*% exp(lr_N))) +
    drop(prep$Xs_N[first, , drop = FALSE] %*% par$theta_N), 0)
  u_P <- ifelse(has_site, log(drop(W %*% exp(lr_P))) +
    drop(prep$Xs_P[first, , drop = FALSE] %*% par$theta_P), 0)
  k <- if (d$use_kappa) par$kappa else 0
  bA <- if (d$free_bA && !is.null(par$b_A)) par$b_A else 1
  g_of <- c(0L, 1L, d$n_g)
  traffic <- type == 1L & rows$class > 1L
  shift <- ifelse(traffic, par$g[pmax(g_of[pmax(rows$class, 1L)], 1L)], 0)
  amp <- ifelse(traffic, k, 1)
  mean <- numeric(n)
  mean[type == 1L] <- (par$b_N + u_N + shift + amp * f_N[t])[type == 1L]
  mean[type == 2L] <- (par$b_P + u_P + f_P[t])[type == 2L]
  if (any(type == 3L)) {
    seasonal <- if (d$acag_seasonal) drop(d$H1 %*% par$gamma_A) else rep(0, d$T)
    mean[type == 3L] <- (par$d_A + par$b_P + bA * f_P[t] + seasonal[t])[type == 3L]
  }
  if (any(type == 4L)) mean[type == 4L] <- (par$a_M + par$b_M * (par$b_P + u_P + f_P[t]))[type == 4L]
  # TROPOMI pseudo-rows (type 5): the precision-weighted mean of the month's cells.
  column <- type == 5L
  if (any(column)) {
    ps <- prep$pseudo
    at_col <- rows$col[column]
    mean[column] <- par$alpha_T + drop(d$H1[t[column], , drop = FALSE] %*% par$gamma_T) +
      drop(ps$Gbar[at_col, , drop = FALSE] %*% par$b_cell)
    if (d$field_N) mean[column] <- mean[column] + par$lambda_S * drop(ps$Bbar[at_col, , drop = FALSE] %*% par$z_N)
  }
  lambda_T <- if (is.null(par$lambda_T)) 0 else par$lambda_T
  l_N <- ifelse(type == 1L, amp, ifelse(column, lambda_T / par$s_N, 0))
  l_P <- ifelse(type == 2L, 1, ifelse(type == 3L, bA, ifelse(type == 4L, par$b_M, 0)))
  lag <- outer(t, t, "-")
  V <- tcrossprod(l_N) * par$s_N^2 * par$rho_N^abs(lag) + tcrossprod(l_P) * par$s_P^2 * par$rho_P^abs(lag)
  if (d$coupled) {
    c0 <- par$r * par$s_N * par$s_P * sqrt((1 - par$rho_N^2) * (1 - par$rho_P^2)) /
      (1 - par$rho_N * par$rho_P)
    cross <- outer(l_N, l_P) * c0 * ifelse(lag >= 0, par$rho_N^pmax(lag, 0), par$rho_P^pmax(-lag, 0))
    V <- V + cross + t(cross)
  }
  same_site <- outer(rows$site, rows$site, "==")
  same_site[is.na(same_site)] <- FALSE
  if (d$use_nugget) {
    at_nugget <- rows$site %in% prep$nugget_sites
    l_nug <- ifelse(at_nugget & type == 2L, 1, ifelse(at_nugget & type == 4L, par$b_M, 0))
    V <- V + par$tau_P^2 * tcrossprod(l_nug) * same_site
  }
  if (any(type == 4L)) {
    moqa <- type == 4L
    same_device <- outer(rows$device, rows$device, "==") & outer(moqa, moqa)
    V <- V + par$tau_M^2 * same_device + par$tau_m^2 * (outer(t, t, "==") & outer(moqa, moqa))
  }
  scalar <- function(v) if (is.null(v)) NA_real_ else v
  column_sd <- rep(NA_real_, n)
  if (any(column)) {
    V <- V + par$sigma_q^2 * par$rho_q^abs(lag) * outer(column, column)
    column_sd[column] <- 1 / sqrt(exp(par$log_kappa_T) * prep$pseudo$S[rows$col[column]])
  }
  noise <- ifelse(type == 1L, ifelse(traffic, par$sig_NT, par$sig_N),
    ifelse(type == 2L, par$sig_P, ifelse(type == 3L, scalar(par$sig_A),
      ifelse(column, column_sd, scalar(par$sig_M)))))
  if (!is.null(rows$latent)) noise[rows$latent %in% TRUE] <- 0
  diag(V) <- diag(V) + noise^2
  list(mean = mean, V = V)
}

# Parameters of draw j by their program names; an absent block is NULL or empty.
.draw <- function(values, j) {
  columns <- colnames(values)
  get <- function(name) if (name %in% columns) unname(values[j, name]) else NULL
  vec <- function(prefix) {
    cols <- grep(paste0("^", prefix, "\\[[0-9]+\\]$"), columns, value = TRUE)
    if (!length(cols)) return(numeric())
    unname(values[j, cols[order(as.integer(sub(".*\\[([0-9]+)\\]$", "\\1", cols)))]])
  }
  list(b_N = get("b_N"), b_P = get("b_P"), g = vec("g"), beta_N = vec("beta_N"),
    beta_P = vec("beta_P"), h_N = vec("h_N"), h_P = vec("h_P"), theta_N = vec("theta_N"),
    theta_P = vec("theta_P"), d_A = get("d_A[1]"), gamma_A = vec("gamma_A"), b_A = get("b_A[1]"),
    s_N = get("s_N"), s_P = get("s_P"), rho_N = get("rho_N"), rho_P = get("rho_P"), r = get("r[1]"),
    kappa = get("kappa[1]"), sig_N = get("sig_N"), sig_NT = get("sig_NT"), sig_P = get("sig_P"),
    sig_A = get("sig_A[1]"), z_N = vec("z_N"), z_P = vec("z_P"), sd_N = get("sd_N[1]"),
    sd_P = get("sd_P[1]"), rho_S = get("rho_S[1]"), tau_P = get("tau_P[1]"), a_M = get("a_M[1]"),
    b_M = get("b_M[1]"), tau_M = get("tau_M[1]"), sig_M = get("sig_M[1]"), sig_C = get("sig_C[1]"),
    tau_m = get("tau_m[1]"), alpha_T = get("alpha_T[1]"), gamma_T = vec("gamma_T"),
    lambda_S = get("lambda_S[1]"), lambda_T = get("lambda_T[1]"), log_kappa_T = get("log_kappa_T[1]"),
    sigma_q = get("sigma_q[1]"), rho_q = get("rho_q[1]"), sigma_beta = get("sigma_beta[1]"),
    b_cell = vec("b_cell"))
}

# Evaluates expr with the L'Ecuyer-CMRG stream of `seed`, then restores the caller's RNG.
.with_seed <- function(seed, expr) {
  kind <- RNGkind()
  had <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  if (had) previous <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
  on.exit({
    RNGkind(kind[1L], kind[2L], kind[3L])
    if (had) assign(".Random.seed", previous, envir = globalenv())
  }, add = TRUE)
  set.seed(seed, kind = "L'Ecuyer-CMRG", normal.kind = "Inversion")
  expr
}
