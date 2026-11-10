# Fit of one model on one unit: the Stan data of stan/fusion.stan, compilation
# with CmdStan's make into work/stan/, four NUTS chains, and the draws read back
# from CmdStan's CSV files with their diagnostics.

# Scale parameters in the order of prior_scale in the program, their tail
# thresholds U with Pr(sigma > U) = 0.05 (quartic priors) and the SDs of the
# half-normal sensitivity priors.
scale_names <- c("s_N", "s_P", "tau_P", "tau_M", "tau_m", "sig_N", "sig_NT", "sig_P", "sig_A",
  "sig_M", "sig_C", "sd_N", "sd_P", "sigma_q", "sigma_beta")
scale_tail_U <- c(s_N = 0.6, s_P = 0.6, tau_P = 0.6, tau_M = 0.6, tau_m = 0.6, sig_N = 1, sig_NT = 1,
  sig_P = 1, sig_A = 1, sig_M = 1, sig_C = 1, sd_N = 1, sd_P = 0.5, sigma_q = 1, sigma_beta = 1)
half_normal_sd <- c(s_N = 0.3, s_P = 0.3, tau_P = 0.3, tau_M = 0.5, tau_m = 0.3, sig_N = 0.5,
  sig_NT = 0.5, sig_P = 0.5, sig_A = 0.5, sig_M = 0.5, sig_C = 0.5, sd_N = 0.449258049,
  sd_P = 0.449258049, sigma_q = 0.5, sigma_beta = 0.5)

# Candidate site covariates (transformed) and ERA5 drivers of each pollutant.
no2_covariates <- c("log_bus_250", "log_bus_500", "log_mtps_100", "log_mtps_500", "log_veh_250",
  "log_dist_mtp", "log_pop_1km")
pm25_covariates <- c(no2_covariates, "acag_log_ratio_municipal")
era5_no2 <- c("log_vent", "t2m")
era5_pm25 <- c("log_vent", "precip")

# Sampler settings of every fit.
sampler <- list(chains = 4L, warmup = 1000L, draws_per_chain = 1000L, adapt_delta = 0.95,
  max_treedepth = 10L)

#' Fit a model on a unit and read the draws back.
#'
#' @param input read_inputs().
#' @param unit a name in validation_units$unit.
#' @param model a name in models$model.
#' @param dir output directory (data.json and the chain CSVs).
#' @param exe compiled program (compile_program()).
#' @return fit: list(model, unit, seed, prep, draws, diagnostics).
fit_model <- function(input, unit, model, dir, exe = compile_program()) {
  part <- partition(input, unit)
  prep <- prepare_fit(input, part, model)
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  write_stan_data(prep$data, file.path(dir, "data.json"))
  started <- Sys.time()
  .run_chains(exe, file.path(dir, "data.json"), dir, part$seed)
  fit <- read_fit(input, unit, model, dir, prep = prep)
  fit$diagnostics$runtime_s <- as.numeric(difftime(Sys.time(), started, units = "secs"))
  fit
}

#' Fit a model on a unit into work/fits/<model>/<unit>/.
#'
#' @return paths of data.json and the chain CSVs.
run_fit <- function(input, unit, model, exe) {
  dir <- file.path("work", "fits", model, unit)
  fit_model(input, unit, model, dir, exe)
  file.path(dir, c("data.json", paste0("chain-", seq_len(sampler$chains), ".csv")))
}

#' A fit from the chain CSVs already in dir (data rebuilt from the inputs).
read_fit <- function(input, unit, model, dir, prep = NULL) {
  part <- partition(input, unit)
  if (is.null(prep)) prep <- prepare_fit(input, part, model)
  chains <- .read_chains(dir)
  list(model = model, unit = unit, seed = part$seed, prep = prep, draws = chains$draws,
    diagnostics = .diagnostics(chains))
}

#' Stan data and the design the prediction reuses.
#'
#' @param input read_inputs().
#' @param part partition().
#' @param model a name in models$model.
#' @return list: data (the Stan data list) and the rows, sites, supports,
#'   covariates, basis and column pseudo-rows of the joint law.
prepare_fit <- function(input, part, model) {
  m <- model_spec(model)
  ground <- c("ground_no2", "ground_pm25")
  held <- part$held
  covariates <- .site_covariates(input$site_covariates)
  cols_P <- if (m$pm25_cov) .select_covariate(part$train, covariates, "ground_pm25", pm25_covariates)
    else character()
  cols_N <- if (m$no2_cov) .select_covariate(part$train, covariates, "ground_no2", no2_covariates)
    else character()
  train <- part$train
  calibration <- if (m$low_cost) part$calibration else part$calibration[0L, , drop = FALSE]
  use_moqa <- m$low_cost == 1L && nrow(part$moqa) > 0L
  if (any(held$source == "moqa") && !use_moqa)
    stop("Held low-cost rows need the low-cost network in the model.", call. = FALSE)
  if (use_moqa) train <- rbind(train, part$moqa)
  acag <- if (m$satpm25) part$acag else part$acag[0L, , drop = FALSE]
  column_months <- if (m$tropomi) part$column$month else character()
  calendar <- .months("2020-01", max(c(train$month, held$month, acag$month, column_months)))
  harmonic <- .harmonic(calendar)
  reference_end <- max(c(train$month[train$source %in% ground], acag$month))
  Z_N <- .era5_matrix(input$era5, if (m$era5) era5_no2 else character(), calendar, reference_end)
  Z_P <- .era5_matrix(input$era5, if (m$era5) era5_pm25 else character(), calendar, reference_end)

  rows_train <- .law_rows(train, calendar)
  if (nrow(acag)) rows_train <- rbind(rows_train, data.frame(observation_id = acag$observation_id,
    type = 3L, t = match(acag$month, calendar), site = NA_character_, class = 0L, device = "",
    x_m = NA_real_, y_m = NA_real_, observed = acag$value, stringsAsFactors = FALSE))
  # Training rows in month order, as the Kalman filter reads them.
  rows_train <- rows_train[order(rows_train$t, rows_train$type), , drop = FALSE]
  rows_train$col <- rep(NA_integer_, nrow(rows_train))
  n_ground <- nrow(rows_train)
  basis <- profile_basis(input$acag_cells, m$knots, m$width_m)
  pseudo <- NULL
  if (m$tropomi) {
    column <- column_design(part$column, input$column_cells, input$acag_cells, calendar, basis)
    pseudo <- column_pseudo_rows(column)
    rows_train <- rbind(rows_train, data.frame(observation_id = paste0("column:", calendar[pseudo$t]),
      type = 5L, t = pseudo$t, site = NA_character_, class = 0L, device = "", x_m = NA_real_,
      y_m = NA_real_, observed = pseudo$ybar, col = seq_along(pseudo$t), stringsAsFactors = FALSE))
  }
  rows_held <- .law_rows(held, calendar)
  rows_held$col <- rep(NA_integer_, nrow(rows_held))
  rows <- rbind(rows_train, rows_held)
  n_train <- nrow(rows_train)
  # Sites of every row with a point support, training sites first.
  with_site <- !is.na(rows$site)
  train_sites <- unique(rows$site[seq_len(n_train)][with_site[seq_len(n_train)]])
  sites <- c(train_sites, setdiff(unique(rows$site[with_site]), train_sites))
  xy <- unique(rows[with_site, c("site", "x_m", "y_m")])
  if (anyDuplicated(xy$site) || any(!is.finite(c(xy$x_m, xy$y_m))))
    stop("Every site needs one finite pair of coordinates.", call. = FALSE)
  xy <- xy[match(sites, xy$site), , drop = FALSE]
  W_all <- profile_support(cbind(xy$x_m, xy$y_m), basis$knots, basis$width_m)
  rownames(W_all) <- sites
  Xs_N <- .covariate_matrix(covariates, cols_N, sites)
  Xs_P <- .covariate_matrix(covariates, cols_P, sites)
  # Reference sites, which carry the PM2.5 deviation: sites with BAM rows.
  nugget_sites <- if (m$deviation) unique(rows$site[rows$type == 2L]) else character()
  nug_train <- intersect(train_sites, unique(rows_train$site[rows_train$type == 2L]))
  nug_train <- nug_train[nug_train %in% nugget_sites]
  nug <- ifelse(train_sites %in% nug_train, match(train_sites, nug_train), 0L)
  train_devices <- unique(rows_train$device[rows_train$type == 4L])
  ground_rows <- rows_train[seq_len(n_ground), , drop = FALSE]
  count <- tabulate(ground_rows$t, nbins = length(calendar))
  n_C <- if (use_moqa) nrow(calibration) else 0L
  column_data <- if (m$tropomi) .column_stan_data(column) else list(n_cells = 1L, K_col = basis$K,
    B_col = matrix(0, 1L, basis$K - 1L), U_cell = matrix(0, 1L, 0L), n_col = 0L,
    y_col = .array(numeric()), d_col = .array(numeric()), cell_col = .array(integer()),
    t_col = .array(integer()), start_col = .array(rep(1L, length(calendar))),
    count_col = .array(rep(0L, length(calendar))))
  prior_scale <- if (m$priors == "quartic") quartic_scale(scale_tail_U) else half_normal_sd
  data <- c(list(T = length(calendar), P_N = ncol(Z_N), P_P = ncol(Z_P), Z_N = Z_N, Z_P = Z_P,
    H = ncol(harmonic), Hm = harmonic, H1 = harmonic,
    K = basis$K, U = basis$U, wbar = basis$wbar,
    field_N = m$no2_profile, field_P = m$pm25_profile,
    couple_space = as.integer(m$profile_cor && m$no2_profile && m$pm25_profile),
    J = length(train_sites), W_site = W_all[train_sites, , drop = FALSE],
    Q_N = ncol(Xs_N), Xs_N = Xs_N[train_sites, , drop = FALSE],
    Q_P = ncol(Xs_P), Xs_P = Xs_P[train_sites, , drop = FALSE],
    use_nugget = m$deviation, Jn = length(nug_train), nug = .array(as.integer(nug)),
    coupled = m$r, use_kappa = m$amplitude, n_g = m$classes, acag_seasonal = m$seasonal,
    free_bA = m$free_gain,
    n = n_ground, ly = log(ground_rows$observed), type = ground_rows$type, cls = ground_rows$class,
    site = ifelse(is.na(ground_rows$site), 0L, match(ground_rows$site, train_sites)),
    start = pmax(1L, cumsum(count) - count + 1L), count = count,
    D = length(train_devices),
    dev = ifelse(ground_rows$type == 4L, match(ground_rows$device, train_devices), 0L),
    n_C = n_C,
    lmoqa_C = .array(if (n_C) log(calibration$moqa) else numeric()),
    lbam_C = .array(if (n_C) log(calibration$bam) else numeric()),
    use_column = m$tropomi, col_temporal = m$temporal_loading,
    sign_constraint = m$sign_constraint),
    column_data, list(quartic = as.integer(m$priors == "quartic"), prior_scale = unname(prior_scale)))
  list(data = data, calendar = calendar, rows = rows, n_train = n_train, n_ground = n_ground,
    held = held, sites = sites, train_sites = train_sites, W_all = W_all, Xs_N = Xs_N, Xs_P = Xs_P,
    nugget_sites = nugget_sites, train_devices = train_devices, basis = basis, pseudo = pseudo,
    use_moqa = use_moqa, covariates = c(no2 = cols_N, pm25 = cols_P))
}

#' Scale a of the half exponential-power prior with Pr(sigma > U) = tail:
#' (sigma / a)^k ~ Gamma(1 / k, 1).
quartic_scale <- function(U, tail = 0.05, k = 4) {
  U / stats::qgamma(tail, shape = 1 / k, lower.tail = FALSE)^(1 / k)
}

#' Write the Stan data as CmdStan JSON (full precision).
write_stan_data <- function(data, path) {
  jsonlite::write_json(data, path = path, auto_unbox = TRUE, factor = "integer",
    always_decimal = FALSE, digits = NA, pretty = TRUE)
  invisible(path)
}

#' Compile stan/fusion.stan with CmdStan's make into dir; make rebuilds only when
#' the program changed. Returns the executable's path.
compile_program <- function(dir = "work/stan", program = "stan/fusion.stan") {
  cmdstan <- Sys.getenv("CMDSTAN")
  if (!nzchar(cmdstan)) stop("Set CMDSTAN to the CmdStan directory.", call. = FALSE)
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  exe <- file.path(normalizePath(dir), "fusion")
  if (!file.copy(program, paste0(exe, ".stan"), overwrite = TRUE, copy.date = TRUE))
    stop("Cannot copy ", program, " to ", dir, call. = FALSE)
  log <- paste0(exe, "-compile.log")
  path <- Sys.getenv("PATH")
  on.exit(Sys.setenv(PATH = path), add = TRUE)
  # The system assembler and linker first.
  Sys.setenv(PATH = paste(unique(c("/usr/bin", "/bin", strsplit(path, .Platform$path.sep,
    fixed = TRUE)[[1L]])), collapse = .Platform$path.sep))
  status <- system2("make", c("-C", shQuote(cmdstan), "STANCFLAGS=--O1", shQuote(exe)),
    stdout = log, stderr = log)
  if (status != 0L)
    stop("Compilation failed:\n", paste(utils::tail(readLines(log), 30L), collapse = "\n"), call. = FALSE)
  exe
}

# Four chains in parallel, each with the unit's seed and its own chain id.
.run_chains <- function(exe, data_file, dir, seed) {
  status <- parallel::mclapply(seq_len(sampler$chains), function(chain) {
    system2(exe, c(paste0("id=", chain), "random", paste0("seed=", seed), "data",
      paste0("file=", data_file), "output", paste0("file=", file.path(dir, paste0("chain-", chain, ".csv"))),
      "refresh=0", "method=sample", paste0("num_samples=", sampler$draws_per_chain),
      paste0("num_warmup=", sampler$warmup), "save_warmup=0", "algorithm=hmc", "engine=nuts",
      paste0("max_depth=", sampler$max_treedepth), "adapt", paste0("delta=", sampler$adapt_delta)),
      stdout = file.path(dir, paste0("chain-", chain, ".log")),
      stderr = file.path(dir, paste0("chain-", chain, ".log")))
  }, mc.cores = sampler$chains, mc.preschedule = FALSE)
  failed <- !vapply(status, identical, logical(1), 0L)
  if (any(failed)) stop("Chains failed: ", paste(which(failed), collapse = ", "), " (see ", dir, ").",
    call. = FALSE)
  invisible(dir)
}

# Draws of the chain CSVs of dir, stacked chain by chain in chain-id order.
.read_chains <- function(dir) {
  files <- list.files(dir, pattern = "[.]csv$", full.names = TRUE)
  header <- lapply(files, function(f) readLines(f, n = 60L))
  sampled <- vapply(header, function(h) any(grepl("^# method = sample", h)), logical(1))
  files <- files[sampled]
  id <- vapply(header[sampled], function(h) as.integer(sub("^# id = ([0-9]+).*$", "\\1",
    grep("^# id = ", h, value = TRUE)[1L])), integer(1))
  if (!identical(sort(id), seq_len(sampler$chains)))
    stop("Expected chains 1..", sampler$chains, " in ", dir, call. = FALSE)
  files <- files[order(id)]
  chains <- lapply(files, function(f) {
    x <- utils::read.csv(f, comment.char = "#", check.names = FALSE)
    names(x) <- sub("[.]([0-9]+)$", "[\\1]", names(x))
    as.matrix(x)
  })
  if (any(vapply(chains, nrow, integer(1)) != sampler$draws_per_chain))
    stop("Each chain must hold ", sampler$draws_per_chain, " draws.", call. = FALSE)
  all <- do.call(rbind, chains)
  diag_cols <- grepl("__$", colnames(all)) & colnames(all) != "lp__"
  list(draws = all[, !diag_cols, drop = FALSE], sampler = all[, diag_cols, drop = FALSE],
    chain = rep(seq_along(chains), each = sampler$draws_per_chain))
}

# Divergences, treedepth hits, the largest rank-normalized split R-hat and the
# smallest bulk ESS over every parameter, and the smallest E-BFMI of a chain.
.diagnostics <- function(chains) {
  by_chain <- function(x) matrix(x, ncol = sampler$chains)
  params <- setdiff(colnames(chains$draws), "lp__")
  rhat <- vapply(params, function(p) rank_rhat(by_chain(chains$draws[, p])), numeric(1))
  ess <- vapply(params, function(p) bulk_ess(by_chain(chains$draws[, p])), numeric(1))
  energy <- by_chain(chains$sampler[, "energy__"])
  ebfmi <- apply(energy, 2L, function(x) (sum(diff(x)^2) / length(x)) / stats::var(x))
  list(divergences = sum(chains$sampler[, "divergent__"]), max_rhat = max(rhat, na.rm = TRUE),
    min_ess_bulk = min(ess, na.rm = TRUE),
    treedepth_hits = sum(chains$sampler[, "treedepth__"] >= sampler$max_treedepth),
    min_ebfmi = min(ebfmi))
}

#' Rank-normalized split R-hat (Vehtari et al. 2021), the larger of the bulk
#' and the folded (tail) value; x is iterations x chains.
rank_rhat <- function(x) {
  max(.rhat(.z_scale(.split_chains(x))), .rhat(.z_scale(.split_chains(abs(x - stats::median(x))))))
}

#' Bulk effective sample size (Vehtari et al. 2021); x is iterations x chains.
bulk_ess <- function(x) .ess(.z_scale(.split_chains(x)))

.split_chains <- function(x) {
  half <- nrow(x) / 2
  cbind(x[seq_len(floor(half)), , drop = FALSE], x[ceiling(half + 1):nrow(x), , drop = FALSE])
}

.z_scale <- function(x) {
  r <- rank(x, ties.method = "average")
  array(stats::qnorm((r - 3 / 8) / (length(r) - 2 * 3 / 8 + 1)), dim = dim(x))
}

.is_constant <- function(x) anyNA(x) || any(!is.finite(x)) || abs(max(x) - min(x)) < .Machine$double.eps

.rhat <- function(x) {
  if (.is_constant(x)) return(NA_real_)
  n <- nrow(x)
  chain_mean <- colMeans(x)
  chain_var <- apply(x, 2L, stats::var)
  sqrt((n * stats::var(chain_mean) / mean(chain_var) + n - 1) / n)
}

.ess <- function(x) {
  chains <- ncol(x)
  n <- nrow(x)
  if (n < 3L || .is_constant(x)) return(NA_real_)
  acov <- apply(x, 2L, .autocovariance)
  acov_means <- rowMeans(acov)
  mean_var <- acov_means[1] * n / (n - 1)
  var_plus <- mean_var * (n - 1) / n
  if (chains > 1) var_plus <- var_plus + stats::var(colMeans(x))
  rho <- rep.int(0, n)
  t <- 0
  rho_even <- 1
  rho[t + 1] <- rho_even
  rho_odd <- 1 - (mean_var - acov_means[t + 2]) / var_plus
  rho[t + 2] <- rho_odd
  while (t < nrow(acov) - 5 && !is.nan(rho_even + rho_odd) && (rho_even + rho_odd > 0)) {
    t <- t + 2
    rho_even <- 1 - (mean_var - acov_means[t + 1]) / var_plus
    rho_odd <- 1 - (mean_var - acov_means[t + 2]) / var_plus
    if ((rho_even + rho_odd) >= 0) {
      rho[t + 1] <- rho_even
      rho[t + 2] <- rho_odd
    }
  }
  max_t <- t
  if (rho_even > 0) rho[max_t + 1] <- rho_even
  # Geyer's initial monotone sequence.
  t <- 0
  while (t <= max_t - 4) {
    t <- t + 2
    if (rho[t + 1] + rho[t + 2] > rho[t - 1] + rho[t]) {
      rho[t + 1] <- (rho[t - 1] + rho[t]) / 2
      rho[t + 2] <- rho[t + 1]
    }
  }
  ess <- chains * n
  tau <- max(-1 + 2 * sum(rho[seq_len(max_t)]) + rho[max_t + 1], 1 / log10(ess))
  ess / tau
}

.autocovariance <- function(x) {
  n <- length(x)
  var_x <- stats::var(x)
  if (var_x == 0) return(rep(0, n))
  yc <- c(x - mean(x), rep.int(0, 2 * stats::nextn(n) - n))
  ac <- Re(stats::fft(abs(stats::fft(yc))^2, inverse = TRUE)[seq_len(n)])
  ac / ac[1] * var_x * (n - 1) / n
}

# Rows of the joint law from targets: type 1 NO2, 2 BAM, 4 low-cost.
.law_rows <- function(targets, calendar) {
  type <- ifelse(targets$source == "ground_no2", 1L, ifelse(targets$source == "ground_pm25", 2L, 4L))
  class <- ifelse(type == 1L, match(targets$site_class, site_classes), 0L)
  if (any(type == 1L & is.na(class)))
    stop("NO2 rows need site classes background, avenue or bus_terminal.", call. = FALSE)
  data.frame(observation_id = targets$observation_id, type = as.integer(type),
    t = match(targets$month, calendar), site = targets$physical_location_id, class = as.integer(class),
    device = ifelse(type == 4L, targets$device_id, ""), x_m = targets$x_m, y_m = targets$y_m,
    observed = targets$observed, stringsAsFactors = FALSE)
}

# Transformed site covariates: log1p of counts, lengths, distance and population.
.site_covariates <- function(raw) {
  data.frame(location_id = raw$location_id, log_bus_250 = log1p(raw$n_bus_stops_250m),
    log_bus_500 = log1p(raw$n_bus_stops_500m), log_mtps_100 = log1p(raw$len_mtps_100m),
    log_mtps_500 = log1p(raw$len_mtps_500m), log_veh_250 = log1p(raw$len_vehicular_250m),
    log_dist_mtp = log1p(raw$dist_mtp_road_m), log_pop_1km = log1p(raw$pop_1km),
    acag_log_ratio_municipal = raw$acag_log_ratio_municipal, stringsAsFactors = FALSE)
}

# In-fold selection: the candidate with the largest absolute correlation with the
# training-location means of the log observations of `source`.
.select_covariate <- function(train, covariates, source, candidates) {
  rows <- train[train$source == source, , drop = FALSE]
  means <- tapply(log(rows$observed), rows$physical_location_id, mean)
  at <- match(names(means), covariates$location_id)
  r <- vapply(candidates, function(col) {
    x <- covariates[[col]][at]
    if (stats::sd(x) == 0) return(0)
    stats::cor(x, as.numeric(means))
  }, numeric(1))
  candidates[which.max(abs(r))]
}

# Covariates standardized over every location of the table [sites x cols].
.covariate_matrix <- function(covariates, cols, sites) {
  if (!length(cols)) return(matrix(0, length(sites), 0L, dimnames = list(sites, NULL)))
  at <- match(sites, covariates$location_id)
  if (anyNA(at)) stop("Site covariates lack sites: ", paste(sites[is.na(at)], collapse = ", "),
    call. = FALSE)
  X <- vapply(cols, function(col) {
    x <- as.numeric(covariates[[col]])
    (x[at] - mean(x)) / stats::sd(x)
  }, numeric(length(sites)))
  matrix(X, nrow = length(sites), dimnames = list(sites, cols))
}

# ERA5 drivers [T x cols] standardized over 2020-01..reference_end, so a whole-year
# prediction never uses later weather in the scaling.
.era5_matrix <- function(era5, cols, calendar, reference_end) {
  if (!length(cols)) return(matrix(0, length(calendar), 0L))
  at <- match(calendar, era5$month)
  if (anyNA(at)) stop("ERA5 lacks months of the model calendar.", call. = FALSE)
  reference <- era5$month >= "2020-01" & era5$month <= reference_end
  Z <- vapply(cols, function(col) {
    x <- as.numeric(era5[[col]])
    (x[at] - mean(x[reference])) / stats::sd(x[reference])
  }, numeric(length(calendar)))
  matrix(Z, nrow = length(calendar), dimnames = list(calendar, cols))
}

# Annual harmonic pair (sin, cos) of each calendar month.
.harmonic <- function(calendar) {
  month <- as.integer(substr(calendar, 6L, 7L))
  cbind(sin(2 * pi * month / 12), cos(2 * pi * month / 12))
}

.months <- function(from, to) {
  format(seq(as.Date(paste0(from, "-01")), as.Date(paste0(to, "-01")), by = "month"), "%Y-%m")
}

.array <- function(x) {
  dim(x) <- length(x)
  x
}
