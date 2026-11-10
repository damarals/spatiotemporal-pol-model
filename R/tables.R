# Result tables of the paper: Table 4 (scores by validation unit), Tables 5-6
# (key backward-elimination contrasts), posterior summaries of the fit to all
# data, and MCMC sampler diagnostics of every fit.

# Rows of Table 4, in the order printed in the paper.
table4_rows <- data.frame(
  pollutant = c(rep("NO2", 4L), rep("PM25", 5L)),
  group = c("omitted-months", "new-traffic", "withheld-background", "whole-years",
    "omitted-months", "new-traffic", "withheld-background", "new-devices", "whole-years"),
  unit = c("Omitted months", "New traffic locations", "Withheld background location", "Whole years",
    "Omitted months", "New traffic locations", "Withheld background location", "New devices",
    "Whole years"),
  stringsAsFactors = FALSE)

#' Table 4 of the paper: predictive scores of the full model without TROPOMI,
#' the final model and the class and pooled benchmarks by validation unit.
scores_table <- function(scores, bench, round_digits = TRUE) {
  if (is.list(bench) && !is.data.frame(bench)) bench <- bench$scores
  sum_final <- summarise_units(scores[scores$model == "final", , drop = FALSE])
  sum_notrop <- summarise_units(scores[scores$model == "full-no-tropomi", , drop = FALSE])
  loc_mean <- function(rows) {
    if (!nrow(rows)) NA_real_ else mean(tapply(rows$crps, rows$physical_location_id, mean))
  }
  rows <- lapply(seq_len(nrow(table4_rows)), function(i) {
    pol <- table4_rows$pollutant[i]
    grp <- table4_rows$group[i]
    sf <- sum_final[sum_final$group == grp & sum_final$pollutant == pol, , drop = FALSE]
    sn <- sum_notrop[sum_notrop$group == grp & sum_notrop$pollutant == pol, , drop = FALSE]
    bc <- loc_mean(.group_rows(bench[bench$benchmark == "class", , drop = FALSE], grp, pol))
    bp <- loc_mean(.group_rows(bench[bench$benchmark == "pooled", , drop = FALSE], grp, pol))
    data.frame(unit = table4_rows$unit[i], group = grp, pollutant = pol,
      n = sf$n_targets, locations = sf$n_locations,
      crps_full_no_tropomi = sn$crps, crps_final = sf$crps,
      crps_class_benchmark = bc, crps_pooled_benchmark = bp,
      skill_class_pct = 100 * (1 - sf$crps / bc),
      skill_pooled_pct = 100 * (1 - sf$crps / bp),
      bias = sf$bias, rmse = sf$rmse, nls = sf$nls,
      coverage_50_pct = 100 * sf$coverage_50, width_50 = sf$width_50,
      coverage_80_pct = 100 * sf$coverage_80, width_80 = sf$width_80,
      coverage_95_pct = 100 * sf$coverage_95, width_95 = sf$width_95,
      stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  if (!round_digits) return(out)
  num3 <- c("crps_full_no_tropomi", "crps_final", "crps_class_benchmark",
    "crps_pooled_benchmark", "nls")
  num2 <- c("skill_class_pct", "skill_pooled_pct", "bias", "rmse")
  num1 <- c("coverage_50_pct", "width_50", "coverage_80_pct", "width_80",
    "coverage_95_pct", "width_95")
  out[num3] <- lapply(out[num3], round, 3)
  out[num2] <- lapply(out[num2], round, 2)
  out[num1] <- lapply(out[num1], round, 1)
  out
}

# Rows selected from elimination_table() for Tables 5 and 6 of the paper.
table5_rows <- utils::read.csv(stringsAsFactors = FALSE, text = "
section,component,reference,group,pollutant
Removed,PM2.5 spatial profile,reduced model,new-traffic,PM25
Removed,PM2.5 spatial profile,reduced model,new-devices,PM25
Removed,PM2.5 deviation,reduced model,new-traffic,PM25
Removed,PM2.5 deviation,reduced model,new-devices,PM25
Removed,Amplitude kappa,full model without TROPOMI,new-traffic,NO2
Removed,Profile correlation,full model without TROPOMI,new-traffic,NO2
Kept,Innovation correlation r,full model without TROPOMI,omitted-months,NO2
Kept,Innovation correlation r,full model without TROPOMI,omitted-months,PM25
Kept,Innovation correlation r,full model without TROPOMI,whole-years,NO2
Kept,Innovation correlation r,full model without TROPOMI,new-devices,PM25
Kept,Two traffic classes,reduced model,new-traffic,NO2
Kept,Two traffic classes,final model,new-traffic,NO2
Kept,NO2 profile,full model without TROPOMI,new-traffic,NO2
Kept,Four profile knots,final model,new-traffic,NO2
")

table6_rows <- utils::read.csv(stringsAsFactors = FALSE, text = "
section,component,reference,group,pollutant
Removed or changed,SatPM2.5 gain fixed at one,reduced model,new-traffic,PM25
Removed or changed,SatPM2.5 gain fixed at one,reduced model,withheld-background,PM25
Removed or changed,SatPM2.5 gain fixed at one,reduced model,omitted-months,PM25
Removed or changed,\"TROPOMI, both loadings\",full model without TROPOMI,new-traffic,NO2
Removed or changed,\"TROPOMI, both loadings\",full model without TROPOMI,new-traffic,PM25
Removed or changed,\"TROPOMI, both loadings\",full model without TROPOMI,omitted-months,NO2
Removed or changed,\"TROPOMI, both loadings\",full model without TROPOMI,omitted-months,PM25
Removed or changed,\"TROPOMI, both loadings\",full model without TROPOMI,whole-years,NO2
Removed or changed,\"TROPOMI, both loadings\",full model without TROPOMI,whole-years,PM25
Removed or changed,\"TROPOMI, both loadings\",full model without TROPOMI,new-devices,PM25
Removed or changed,NO2 site covariate,reduced model,new-traffic,NO2
Kept,ERA5 drivers,full model without TROPOMI,omitted-months,NO2
Kept,ERA5 drivers,full model without TROPOMI,whole-years,NO2
Kept,ERA5 drivers,full model without TROPOMI,whole-years,PM25
Kept,ERA5 drivers,full model without TROPOMI,new-traffic,PM25
Kept,SatPM2.5,full model without TROPOMI,new-traffic,PM25
Kept,SatPM2.5,full model without TROPOMI,whole-years,PM25
Kept,SatPM2.5,full model without TROPOMI,omitted-months,PM25
Kept,SatPM2.5,full model without TROPOMI,omitted-months,NO2
Kept,SatPM2.5 seasonal factor,full model without TROPOMI,withheld-background,PM25
Kept,SatPM2.5 seasonal factor,full model without TROPOMI,whole-years,PM25
Kept,SatPM2.5 seasonal factor,full model without TROPOMI,omitted-months,PM25
Kept,SatPM2.5 seasonal factor,full model without TROPOMI,new-traffic,PM25
Kept,Low-cost network,reduced model,new-traffic,PM25
Kept,Low-cost network,final model without TROPOMI,new-traffic,PM25
Kept,Low-cost network,final model without TROPOMI,omitted-months,PM25
Kept,\"TROPOMI, spatial loading\",final model,new-traffic,NO2
Kept,PM2.5 site covariate,full model without TROPOMI,new-traffic,PM25
Kept,PM2.5 site covariate,full model without TROPOMI,withheld-background,PM25
Kept,PM2.5 site covariate,full model without TROPOMI,new-devices,PM25
")

#' Select the rows of Table 5 (structure and coupling) or Table 6 (sources)
#' from the 186-row elimination table.
paper_contrasts_table <- function(elimination, which = c("structure", "sources")) {
  which <- match.arg(which)
  spec <- if (which == "structure") table5_rows else table6_rows
  key <- function(x) paste(x$component, x$reference, x$group, x$pollutant, sep = "|")
  at <- match(key(spec), key(elimination))
  if (anyNA(at)) stop("Missing rows in elimination table for ", which, ".", call. = FALSE)
  cbind(section = spec$section, elimination[at, , drop = FALSE], row.names = NULL)
}

#' Posterior summary table (41 structural and observation parameters) from the
#' fit of the final model to all data (2020-2025).
posterior_table <- function(fit, chains = 4L) {
  draws <- if (is.matrix(fit)) fit else fit$draws
  v <- colnames(draws)
  keep <- v[!grepl("^(lp__|z_|zc|b_cell|lr_|u_|f_|a[.\\[]|mu_|log_lik|Pm|r_col)", v)]
  by_chain <- function(x) matrix(x, ncol = chains)
  rows <- lapply(keep, function(p) {
    x <- draws[, p]
    s <- stats::sd(x)
    if (!is.finite(s) || s == 0) return(NULL)
    q <- stats::quantile(x, c(0.025, 0.975), names = FALSE)
    data.frame(variable = sub("[.]([0-9]+)$", "[\\1]", p), mean = mean(x), sd = s,
      `2.5%` = q[1L], `97.5%` = q[2L], rhat = rank_rhat(by_chain(x)),
      ess_bulk = bulk_ess(by_chain(x)), check.names = FALSE, stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

# Columns of the per-target table in results/.
target_columns <- c("model", "unit", "observation_id", "pollutant", "physical_location_id", "month",
  "observed", "predictive_mean", "predictive_sd", "crps", "crps_mcse", "nls", "pit", "coverage_50",
  "coverage_80", "coverage_95")

#' Per-target scores of every fit at the reference stations. The low-cost
#' targets are left out: with the predictive distribution, their scores would
#' give away the network's monthly values, which are not public.
target_table <- function(scores) {
  x <- scores[scores$source != "moqa", target_columns, drop = FALSE]
  x <- x[order(match(x$model, models$model), match(x$unit, validation_units$unit)), , drop = FALSE]
  rownames(x) <- NULL
  x
}

#' Coverage of the central 50, 80 and 95 % predictive intervals over the pooled
#' evaluation targets (pooled_targets()), for each model fitted on all eleven
#' validation units.
coverage_table <- function(scores) {
  units <- validation_units$unit[validation_units$unit != "all-data"]
  rows <- lapply(models$model, function(m) {
    x <- scores[scores$model == m, , drop = FALSE]
    if (!all(units %in% x$unit)) return(NULL)
    x <- pooled_targets(x)
    data.frame(model = m, n_targets = nrow(x), coverage_50 = mean(x$coverage_50),
      coverage_80 = mean(x$coverage_80), coverage_95 = mean(x$coverage_95), stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

#' Scores at the background location: each whole year and the withheld
#' background location, for one model and the class and pooled benchmarks.
#'
#' @return data.frame per unit and pollutant: n_targets, mean_observed,
#'   mean_predictive, crps (model), crps_class, crps_pooled and class_pool (the
#'   training values of the class benchmark: the location's own history; NA at
#'   the withheld background location, which has none).
background_table <- function(scores, bench, model = "final") {
  if (is.list(bench) && !is.data.frame(bench)) bench <- bench$scores
  rows <- list()
  for (u in c("year-2023", "year-2024", "year-2025", "withheld-background")) {
    for (pollutant in c("NO2", "PM25")) {
      x <- .group_rows(scores[scores$model == model, , drop = FALSE], u, pollutant)
      kind <- function(k) .group_rows(bench[bench$benchmark == k, , drop = FALSE], u, pollutant)
      class <- kind("class")
      pooled <- kind("pooled")
      if (!nrow(x) || !nrow(pooled)) stop("No ", pollutant, " targets in ", u, ".", call. = FALSE)
      pool <- unique(class$n_training)
      if (length(pool) > 1L) stop("The class pool differs between targets in ", u, ".", call. = FALSE)
      rows[[length(rows) + 1L]] <- data.frame(model = model, unit = u, pollutant = pollutant,
        n_targets = nrow(x), mean_observed = mean(x$observed), mean_predictive = mean(x$predictive_mean),
        crps = mean(x$crps), crps_class = if (nrow(class)) mean(class$crps) else NA_real_,
        crps_pooled = mean(pooled$crps), class_pool = if (length(pool)) pool else NA_integer_,
        stringsAsFactors = FALSE)
    }
  }
  do.call(rbind, rows)
}

# Transformed parameters whose posterior the paper quotes.
derived_parameters <- data.frame(
  variable = c("exp(g[1])", "exp(g[2])", "exp(d_A[1])", "kappa_T", "z_N[1]"),
  draw = c("g[1]", "g[2]", "d_A[1]", "log_kappa_T[1]", "z_N[1]"),
  transform = c("exp", "exp", "exp", "exp", "identity"), stringsAsFactors = FALSE)

#' Posterior mean, median and 2.5, 5, 95 and 97.5 % quantiles of the derived
#' parameters (traffic class factors, SatPM2.5 factor, TROPOMI precision kappa_T
#' and the first NO2 profile coordinate) present in each fit.
derived_table <- function(fits) {
  rows <- list()
  for (f in fits) {
    for (i in seq_len(nrow(derived_parameters))) {
      v <- derived_parameters$draw[i]
      if (!v %in% colnames(f$draws)) next
      x <- f$draws[, v]
      if (derived_parameters$transform[i] == "exp") x <- exp(x)
      q <- stats::quantile(x, c(0.025, 0.05, 0.5, 0.95, 0.975), names = FALSE)
      rows[[length(rows) + 1L]] <- data.frame(model = f$model, unit = f$unit,
        variable = derived_parameters$variable[i], mean = mean(x), median = q[3L], q025 = q[1L],
        q05 = q[2L], q95 = q[4L], q975 = q[5L], stringsAsFactors = FALSE)
    }
  }
  do.call(rbind, rows)
}

#' Write a result table with write.csv (15 significant digits) and return its path.
write_result <- function(x, path) {
  utils::write.csv(x, path, row.names = FALSE)
  path
}

#' One-row-per-fit MCMC diagnostics table: model, unit, divergences, max R-hat,
#' min bulk ESS, treedepth hits, min E-BFMI. Accepts a list of fit objects or a
#' data.frame of diagnostics.
diagnostics_table <- function(fits) {
  cols <- c("model", "unit", "divergences", "max_rhat", "min_ess_bulk", "treedepth_hits",
    "min_ebfmi")
  if (is.data.frame(fits)) {
    if (!all(cols %in% names(fits)))
      stop("Diagnostics table needs columns: ", paste(cols, collapse = ", "), call. = FALSE)
    out <- fits[cols]
  } else {
    rows <- lapply(fits, function(f) {
      d <- f$diagnostics
      data.frame(model = f$model, unit = f$unit,
        divergences = as.integer(d$divergences),
        max_rhat = as.numeric(d$max_rhat),
        min_ess_bulk = as.numeric(d$min_ess_bulk),
        treedepth_hits = as.integer(d$treedepth_hits),
        min_ebfmi = as.numeric(d$min_ebfmi), stringsAsFactors = FALSE)
    })
    out <- do.call(rbind, rows)
  }
  out <- out[order(match(out$model, models$model), match(out$unit, validation_units$unit)), ,
    drop = FALSE]
  rownames(out) <- NULL
  out
}