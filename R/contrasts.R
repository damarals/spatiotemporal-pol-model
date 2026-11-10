# Paired comparisons of model variants: the 26 steps (183 comparisons) of the
# backward elimination and sensitivity analyses, per-location breakdowns, and
# the location-by-year block bootstrap over all 133 evaluation targets.

# The 26 steps of the backward elimination and sensitivity analyses, in order.
elimination_steps <- utils::read.csv(stringsAsFactors = FALSE, text = "
reference,component,change,model_with,model_without
full model without TROPOMI,Low-cost network,removed,full-no-tropomi,full-no-tropomi-without-low-cost
full model without TROPOMI,PM2.5 deviation,removed,full-no-tropomi,full-no-tropomi-without-deviation
full model without TROPOMI,PM2.5 spatial profile,removed,full-no-tropomi,full-no-tropomi-without-pm25-profile
full model without TROPOMI,NO2 site covariate,removed,full-no-tropomi,full-no-tropomi-without-no2-cov
full model without TROPOMI,NO2 profile,removed,full-no-tropomi,full-no-tropomi-without-no2-profile
full model without TROPOMI,Amplitude kappa,removed,full-no-tropomi,full-no-tropomi-without-amplitude
full model without TROPOMI,Two traffic classes,removed,full-no-tropomi,full-no-tropomi-one-class
full model without TROPOMI,Innovation correlation r,removed,full-no-tropomi,full-no-tropomi-without-r
full model without TROPOMI,ERA5 drivers,removed,full-no-tropomi,full-no-tropomi-without-era5
full model without TROPOMI,PM2.5 site covariate,removed,full-no-tropomi,full-no-tropomi-without-pm25-cov
full model without TROPOMI,SatPM2.5 gain fixed at one,switched,full-no-tropomi,full-no-tropomi-free-gain
full model without TROPOMI,SatPM2.5 seasonal factor,removed,full-no-tropomi,full-no-tropomi-without-seasonal
full model without TROPOMI,SatPM2.5,removed,full-no-tropomi,full-no-tropomi-without-satpm25
full model without TROPOMI,Profile correlation,removed,full-no-tropomi,full-no-tropomi-without-profile-cor
full model without TROPOMI,\"TROPOMI, both loadings\",removed,full,full-no-tropomi
reduced model,Low-cost network,removed,reduced,reduced-without-low-cost
reduced model,PM2.5 deviation,removed,reduced,reduced-without-deviation
reduced model,PM2.5 spatial profile,removed,reduced,reduced-without-pm25-profile
reduced model,NO2 site covariate,removed,reduced,reduced-without-no2-cov
reduced model,Two traffic classes,removed,reduced,reduced-one-class
reduced model,SatPM2.5 gain fixed at one,switched,reduced,reduced-free-gain
final model without TROPOMI,Low-cost network,removed,final-no-tropomi,final-no-tropomi-without-low-cost
final model,\"TROPOMI, spatial loading\",removed,final,final-no-tropomi
final model,Two traffic classes,removed,final,final-one-class
final model,Four profile knots,switched,final,final-eight-knots
final model,Quartic-tail scale priors,switched,final,final-half-normal
")

#' Paired comparison of two models on one unit group and pollutant: equal
#' weight per location, with paired Monte Carlo SE (correlated bound and
#' independent-error version) and uncertainty across locations.
paired_contrast <- function(scores_with, scores_without, group, pollutant, seed, B = 4000L) {
  A <- .group_rows(scores_with, group, pollutant)
  W <- .group_rows(scores_without, group, pollutant)
  if (!nrow(A) || !nrow(W)) return(NULL)
  key <- function(x) paste(x$unit, x$observation_id)
  keep <- intersect(key(A), key(W))
  if (!length(keep)) return(NULL)
  A <- A[match(keep, key(A)), , drop = FALSE]
  W <- W[match(keep, key(W)), , drop = FALSE]
  if (max(abs(A$observed - W$observed)) > 1e-6)
    stop("Observed values differ between ", A$model[1L], " and ", W$model[1L], " in ", group, " ",
      pollutant, ".", call. = FALSE)
  loc <- A$physical_location_id
  w <- 1 / (length(unique(loc)) * as.numeric(table(loc)[loc]))
  d <- W$crps - A$crps
  m <- sqrt(A$crps_mcse^2 + W$crps_mcse^2)
  sa <- summarise_locations(A)
  sw <- summarise_locations(W)
  u <- summarise_differences(d, loc, seed, B = B)
  data.frame(group = group, pollutant = pollutant, n = nrow(A), locations = sa$n_locations,
    crps_with = sa$crps, crps_without = sw$crps, gain = sum(w * d),
    mcse = sum(w * m), mcse_indep = sqrt(sum(w^2 * m^2)),
    se_between = u$se_between, se_basis = u$se_basis, locations_agreeing = u$agree,
    drop1_min = u$drop1_min, drop1_max = u$drop1_max, stable = u$stable,
    boot95_lo = u$boot_lo, boot95_hi = u$boot_hi, stringsAsFactors = FALSE)
}

#' The pooled evaluation targets: every target except the prior-dominated NO2 at
#' the withheld background location.
pooled_targets <- function(scores) {
  scores[!(scores$unit == "withheld-background" & scores$pollutant == "NO2"), , drop = FALSE]
}

#' Mean paired CRPS gain by location at the new traffic locations, for every step
#' of the backward elimination and sensitivity analyses and both pollutants.
#'
#' @param scores per-target scores of every fit.
#' @return data.frame: the step (reference, component, change, model_with,
#'   model_without), pollutant, location_id, site_class, n_targets and gain (CRPS
#'   without minus CRPS with).
location_gains_table <- function(scores, steps = elimination_steps) {
  rows <- list()
  key <- function(x) paste(x$unit, x$observation_id)
  for (i in seq_len(nrow(steps))) for (pollutant in c("NO2", "PM25")) {
    sa <- scores[scores$model == steps$model_with[i], , drop = FALSE]
    sb <- scores[scores$model == steps$model_without[i], , drop = FALSE]
    A <- .group_rows(sa, "new-traffic", pollutant)
    B <- .group_rows(sb, "new-traffic", pollutant)
    A <- A[key(A) %in% key(B), , drop = FALSE]
    if (!nrow(A)) next
    gain <- location_gains(sa, sb, "new-traffic", pollutant)
    location <- names(gain)
    rows[[length(rows) + 1L]] <- data.frame(steps[i, c("reference", "component", "change",
      "model_with", "model_without")], pollutant = pollutant, location_id = location,
      site_class = A$site_class[match(location, A$physical_location_id)],
      n_targets = as.integer(table(A$physical_location_id)[location]), gain = as.numeric(gain),
      row.names = NULL, stringsAsFactors = FALSE)
  }
  do.call(rbind, rows)
}

#' Mean paired CRPS gain (CRPS without minus CRPS with) by physical location.
location_gains <- function(scores_with, scores_without, group, pollutant) {
  A <- .group_rows(scores_with, group, pollutant)
  W <- .group_rows(scores_without, group, pollutant)
  key <- function(x) paste(x$unit, x$observation_id)
  keep <- intersect(key(A), key(W))
  A <- A[match(keep, key(A)), , drop = FALSE]
  W <- W[match(keep, key(W)), , drop = FALSE]
  tapply(W$crps - A$crps, A$physical_location_id, mean)
}

#' All 183 paired comparisons of the backward elimination and sensitivity
#' analyses across the five validation unit groups.
elimination_table <- function(scores, steps = elimination_steps, round_digits = TRUE) {
  rows <- list()
  for (i in seq_len(nrow(steps))) {
    sa <- scores[scores$model == steps$model_with[i], , drop = FALSE]
    sb <- scores[scores$model == steps$model_without[i], , drop = FALSE]
    j <- 0L
    for (u in seq_len(nrow(evaluation_groups))) {
      grp <- evaluation_groups$group[u]
      pol <- evaluation_groups$pollutant[u]
      # First check if common rows exist so j increments exactly as in the paper kit.
      A <- .group_rows(sa, grp, pol)
      B <- .group_rows(sb, grp, pol)
      if (!nrow(A) || !nrow(B)) next
      key <- function(x) paste(x$unit, x$observation_id)
      if (!length(intersect(key(A), key(B)))) next
      j <- j + 1L
      seed <- 20261109L + 100L * i + j
      row <- paired_contrast(sa, sb, grp, pol, seed)
      rows[[length(rows) + 1L]] <- data.frame(
        reference = steps$reference[i], component = steps$component[i], change = steps$change[i],
        model_with = steps$model_with[i], model_without = steps$model_without[i],
        unit = evaluation_groups$unit_label[u], group = grp, pollutant = pol,
        n = row$n, locations = row$locations, crps_with = row$crps_with,
        crps_without = row$crps_without, gain = row$gain, mcse = row$mcse,
        se_between = row$se_between, se_basis = row$se_basis,
        locations_agreeing = row$locations_agreeing, drop1_min = row$drop1_min,
        drop1_max = row$drop1_max, stable = row$stable,
        boot95_lo = row$boot95_lo, boot95_hi = row$boot95_hi, bootstrap_seed = ifelse(
          row$se_basis == "locations", seed, NA_integer_), stringsAsFactors = FALSE)
    }
  }
  out <- do.call(rbind, rows)
  if (!round_digits) return(out)
  num3 <- c("crps_with", "crps_without", "gain", "mcse", "se_between", "drop1_min", "drop1_max",
    "boot95_lo", "boot95_hi")
  out[num3] <- lapply(out[num3], round, 3)
  out
}

#' Pooled comparison over all common targets except the prior-dominated NO2 at
#' the withheld background location, with a location-by-year block bootstrap.
#' Targets and blocks are sorted (unit, observation id; block label) before the
#' resampling, so the result does not depend on the order of the input rows.
pooled_contrast <- function(scores_with, scores_without, seed = 20261109L + 9001L, B = 4000L,
    drop_withheld_background_no2 = TRUE) {
  pick <- function(x) if (drop_withheld_background_no2) pooled_targets(x) else x
  A <- pick(scores_with)
  R <- pick(scores_without)
  key <- function(x) paste(x$unit, x$observation_id)
  common <- sort(intersect(key(A), key(R)), method = "radix")
  if (!length(common)) stop("No common targets for pooled contrast.", call. = FALSE)
  A <- A[match(common, key(A)), , drop = FALSE]
  R <- R[match(common, key(R)), , drop = FALSE]
  d <- R$crps - A$crps
  block <- paste(A$physical_location_id, substr(A$month, 1L, 4L))
  blocks <- sort(unique(block), method = "radix")
  set.seed(seed)
  boot <- replicate(B, {
    chosen <- sample(blocks, length(blocks), replace = TRUE)
    mean(unlist(lapply(chosen, function(b) d[block == b])))
  })
  data.frame(model_with = A$model[1L], model_without = R$model[1L],
    n_targets = length(d), n_blocks = length(blocks), gain = mean(d),
    boot95_lo = stats::quantile(boot, 0.025, names = FALSE),
    boot95_hi = stats::quantile(boot, 0.975, names = FALSE),
    p_with_better = mean(boot > 0), bootstrap_seed = seed, stringsAsFactors = FALSE)
}