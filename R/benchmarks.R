# Empirical benchmarks (class and pooled) and benchmark comparisons by
# validation unit.

benchmark_units <- c("new-traffic-1", "new-traffic-2", "new-traffic-3", "omitted-months",
  "year-2023", "year-2024", "year-2025", "new-devices-1", "new-devices-2", "new-devices-3",
  "withheld-background")

benchmark_fields <- c("observation_id", "pollutant", "physical_location_id", "month", "observed",
  "source", "site_class")

# Unit groups evaluated in Tables 4-6 and in the benchmark comparisons.
evaluation_groups <- data.frame(
  group = c("new-traffic", "new-traffic", "withheld-background", "withheld-background",
    "omitted-months", "omitted-months", "whole-years", "whole-years", "new-devices"),
  unit_label = c("New traffic locations", "New traffic locations",
    "Withheld background location", "Withheld background location",
    "Omitted months", "Omitted months", "Whole years", "Whole years", "New devices"),
  pollutant = c("NO2", "PM25", "NO2", "PM25", "NO2", "PM25", "NO2", "PM25", "PM25"),
  stringsAsFactors = FALSE)

#' Training values one held target is scored against. Every pool keeps the
#' target's source. The class pool keeps its site class (at the background
#' location, its own training history; every low-cost device is of class
#' "sensor"). At a reference traffic target, the pooled benchmark pools every
#' training traffic location and leaves the background location out for both
#' pollutants; pooled_with_background keeps the background location in the
#' pool for both pollutants, and pooled_without_background leaves the background
#' location out even when the target is at the background location.
benchmark_pool <- function(part, target,
    kind = c("class", "pooled", "pooled_with_background", "pooled_without_background")) {
  kind <- match.arg(kind)
  train <- rbind(part$train, part$moqa)
  same <- train[train$source == target$source, , drop = FALSE]
  if (kind == "class") return(same$observed[same$site_class == target$site_class])
  traffic <- target$physical_location_id != background_location
  without_bg <- switch(kind,
    pooled = target$source %in% c("ground_no2", "ground_pm25") && traffic,
    pooled_with_background = FALSE,
    pooled_without_background = TRUE)
  if (without_bg) same$observed[same$physical_location_id != background_location] else same$observed
}

#' CRPS of the empirical distribution of training_values at each observed value.
empirical_crps <- function(training_values, observed) {
  if (!is.numeric(training_values) || !length(training_values) || any(!is.finite(training_values)))
    stop("Training values must be finite and non-empty.", call. = FALSE)
  if (!is.numeric(observed) || !length(observed) || any(!is.finite(observed)))
    stop("Observed values must be finite and non-empty.", call. = FALSE)
  vapply(observed, function(y)
    scoringRules::crps_sample(y = y, dat = training_values, method = "edf"), numeric(1))
}

#' Per-target benchmark scores across validation units. Targets whose pool is
#' empty (the class benchmark at the withheld background location) are returned
#' in `unsupported`.
score_benchmarks <- function(input, units = benchmark_units, variants = TRUE) {
  scores <- unsupported <- list()
  for (u in units) {
    part <- partition(input, u)
    for (i in seq_len(nrow(part$held))) {
      target <- part$held[i, , drop = FALSE]
      ground <- target$source %in% c("ground_no2", "ground_pm25")
      at_bg <- target$physical_location_id == background_location
      kinds <- c("class", "pooled",
        if (variants && ground && !at_bg) "pooled_with_background",
        if (variants && target$source == "ground_no2" && at_bg) "pooled_without_background")
      for (kind in kinds) {
        pool <- benchmark_pool(part, target, kind)
        row <- data.frame(unit = u, benchmark = kind, target[benchmark_fields],
          n_training = length(pool), stringsAsFactors = FALSE, row.names = NULL)
        if (length(pool)) {
          row$crps <- empirical_crps(pool, target$observed)
          scores[[length(scores) + 1L]] <- row
        } else {
          row$reason <- "no_training_support"
          unsupported[[length(unsupported) + 1L]] <- row
        }
      }
    }
  }
  list(scores = do.call(rbind, scores), unsupported = do.call(rbind, unsupported))
}

#' Uncertainty of an equal-location mean difference across held-out locations:
#' SD of location means over sqrt(L), locations agreeing in sign, one-location
#' omissions (stable when a strict majority agrees and every omission keeps the
#' sign), and 95% percentile bootstrap interval over B location resamples. With
#' one location only, the SE over targets is given.
summarise_differences <- function(d, location, seed, B = 4000L) {
  dl <- tapply(d, location, mean)
  L <- length(dl)
  effect <- mean(dl)
  if (L == 1L)
    return(data.frame(effect = effect, se_between = stats::sd(d) / sqrt(length(d)),
      se_basis = "targets", agree = NA_integer_, drop1_min = NA_real_, drop1_max = NA_real_,
      stable = NA, boot_lo = NA_real_, boot_hi = NA_real_, stringsAsFactors = FALSE))
  drop1 <- vapply(seq_len(L), function(l) mean(dl[-l]), 0)
  set.seed(seed)
  boot <- replicate(B, mean(sample(dl, L, replace = TRUE)))
  agree <- sum(sign(dl) == sign(effect))
  data.frame(effect = effect, se_between = stats::sd(dl) / sqrt(L), se_basis = "locations",
    agree = agree, drop1_min = min(drop1), drop1_max = max(drop1),
    stable = agree > L / 2 && all(sign(drop1) == sign(effect)),
    boot_lo = stats::quantile(boot, 0.025, names = FALSE),
    boot_hi = stats::quantile(boot, 0.975, names = FALSE), stringsAsFactors = FALSE)
}

#' Compare one benchmark against one model on one unit group and pollutant.
benchmark_unit <- function(bench, model_scores, group, pollutant, kind, seed) {
  B <- .group_rows(bench[bench$benchmark == kind, , drop = FALSE], group, pollutant)
  M <- .group_rows(model_scores, group, pollutant)
  if (!nrow(B) || !nrow(M)) return(NULL)
  key <- function(x) paste(x$unit, x$observation_id)
  at <- match(key(B), key(M))
  if (anyNA(at) || nrow(B) != nrow(M))
    stop(M$model[1L], " and the ", kind, " benchmark differ in their targets of ", group, " ",
      pollutant, ".", call. = FALSE)
  M <- M[at, , drop = FALSE]
  if (max(abs(B$observed - M$observed)) > 1e-6)
    stop("Observed values differ between ", M$model[1L], " and the benchmark in ", group, " ",
      pollutant, ".", call. = FALSE)
  location_mean <- function(x) mean(tapply(x, B$physical_location_id, mean))
  reference <- location_mean(B$crps)
  crps <- location_mean(M$crps)
  cbind(data.frame(group = group, pollutant = pollutant, benchmark = kind, model = M$model[1L],
    n = nrow(B), locations = length(unique(B$physical_location_id)), benchmark_crps = reference,
    model_crps = crps, skill = 1 - crps / reference, stringsAsFactors = FALSE),
    summarise_differences(B$crps - M$crps, B$physical_location_id, seed))
}

#' Compare benchmarks against reference models across the evaluation unit groups.
compare_benchmarks <- function(bench, scores,
    models = c("final", "final-no-tropomi", "full-no-tropomi", "reduced", "full"),
    variants = FALSE, round_digits = TRUE) {
  if (is.list(bench) && !is.data.frame(bench)) bench <- bench$scores
  kinds <- c("class", "pooled",
    if (variants) c("pooled_without_background", "pooled_with_background"))
  rows <- list()
  for (m in seq_along(models)) {
    mod <- models[m]
    ms <- scores[scores$model == mod, , drop = FALSE]
    for (u in seq_len(nrow(evaluation_groups))) {
      grp <- evaluation_groups$group[u]
      pol <- evaluation_groups$pollutant[u]
      for (k in seq_along(kinds)) {
        kind <- kinds[k]
        if (kind == "pooled_without_background" &&
            !(grp %in% c("omitted-months", "whole-years") && pol == "NO2")) next
        if (kind == "pooled_with_background" && grp != "new-traffic") next
        seed <- 20261109L + 11000L + 100L * m + 10L * u + k
        row <- benchmark_unit(bench, ms, grp, pol, kind, seed)
        if (is.null(row)) next
        row$unit <- evaluation_groups$unit_label[u]
        row$bootstrap_seed <- ifelse(row$se_basis == "locations", seed, NA_integer_)
        rows[[length(rows) + 1L]] <- row
      }
    }
  }
  res <- do.call(rbind, rows)
  out <- data.frame(unit = res$unit, group = res$group, pollutant = res$pollutant,
    benchmark = res$benchmark, model = res$model, n = res$n, locations = res$locations,
    benchmark_crps = res$benchmark_crps, model_crps = res$model_crps,
    skill_pct = 100 * res$skill, difference = res$effect, se_between = res$se_between,
    se_basis = res$se_basis, locations_agreeing = res$agree, drop1_min = res$drop1_min,
    drop1_max = res$drop1_max, stable = res$stable, boot95_lo = res$boot_lo,
    boot95_hi = res$boot_hi, bootstrap_seed = res$bootstrap_seed, stringsAsFactors = FALSE)
  if (!round_digits) return(out)
  num3 <- c("benchmark_crps", "model_crps", "difference", "se_between", "drop1_min", "drop1_max",
    "boot95_lo", "boot95_hi")
  out[num3] <- lapply(out[num3], round, 3)
  out$skill_pct <- round(out$skill_pct, 2)
  out
}

# Rows of `x` belonging to a summary group (or single validation unit) and pollutant.
.group_rows <- function(x, group, pollutant) {
  members <- if (group %in% summary_groups$group)
    strsplit(summary_groups$units[summary_groups$group == group], ";", fixed = TRUE)[[1L]]
  else group
  if (!all(members %in% x$unit)) return(x[0L, , drop = FALSE])
  devices <- startsWith(group, "new-devices")
  source <- if (devices) "moqa" else if (pollutant == "NO2") "ground_no2" else "ground_pm25"
  x[x$unit %in% members & x$pollutant == pollutant & x$source == source, , drop = FALSE]
}