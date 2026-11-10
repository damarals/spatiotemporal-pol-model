# Quality control and monthly means of the reference station (CETESB 2020 coverage rules).

# Value as stored in the model inputs: CSV text written by data.table::fwrite (15 significant
# digits, rounded on its own decimal expansion, which can differ from sprintf in the last digit).
stored_value <- function(x) {
  path <- tempfile(fileext = ".csv")
  on.exit(unlink(path))
  data.table::fwrite(list(x = x), path, na = "")
  as.numeric(readLines(path)[-1])
}

# Hourly values -> daily and monthly means. A day needs at least 16 valid hours; a month needs
# valid days on at least two thirds of its calendar days, rounded up; each mean is the plain
# mean of the level below. hours: key, day (YYYY-MM-DD), value, valid. grouped = TRUE takes
# data.table's grouped mean (one long-double pass, as the reference PM2.5 tables were built)
# instead of base::mean. Value: one row per key and month with valid_days, required_days and
# value (NA when coverage is insufficient).
coverage_means <- function(hours, grouped = FALSE) {
  h <- data.table::as.data.table(hours)[, .(key, day, value, valid)]
  if (grouped) {
    days <- h[valid == TRUE, .(valid_hours = .N, value = mean(value)), by = .(key, day)]
    days <- days[valid_hours >= 16L][, month := substr(day, 1, 7)]
    months <- days[, .(valid_days = .N, value = mean(value)), by = .(key, month)]
    months[, required_days := vapply(month, function(m) as.integer(ceiling(2L * days_in_month(m) / 3L)),
                                     integer(1))]
    months[valid_days < required_days, value := NA_real_]
    return(months[, .(key, month, valid_days, required_days, value)])
  }
  days <- h[, {
    n <- sum(valid)
    list(valid_hours = n, value = if (n >= 16L) mean(value[valid]) else NA_real_)
  }, by = .(key, day)]
  days[, `:=`(valid_day = valid_hours >= 16L, month = substr(day, 1, 7))]
  months <- days[, {
    required <- as.integer(ceiling(2L * days_in_month(month) / 3L))
    n <- sum(valid_day)
    list(valid_days = n, required_days = required,
         value = if (n >= required) mean(value[valid_day]) else NA_real_)
  }, by = .(key, month)]
  months[]
}

days_in_month <- function(month) {
  start <- as.Date(paste0(month[1], "-01"))
  as.integer(seq(start, by = "month", length.out = 2L)[2] - start)
}

# NO2 monthly means. samples: read_seuma_samples(); deployments: data/decisions/seuma-deployments.csv.
# Value: data.frame location_id, station_id, month, value (ug/m3) for the months with coverage.
seuma_no2_monthly <- function(samples, deployments) {
  x <- data.table::copy(samples)
  x <- x[format(timestamp, "%Y-%m", tz = "America/Fortaleza") == month]
  x[, valid := .no2_valid(x)]
  x[, no2_ugm3 := data.table::fifelse(no2_unit == "ppb", no2 * 46.0055 / 24.45, no2)]
  x <- .seuma_collapse_duplicates(x)
  x[, `:=`(key = paste(station_id, month), day = format(timestamp, "%Y-%m-%d", tz = "America/Fortaleza"),
           hour = format(timestamp, "%Y-%m-%dT%H", tz = "America/Fortaleza"))]
  # Every workbook has one sample per hour, so the 75 % sample rule asks for one valid sample.
  hours <- x[, list(day = day[1], valid = any(valid),
                    value = if (any(valid)) mean(no2_ugm3[valid]) else NA_real_), by = .(key, hour)]
  m <- coverage_means(hours)[!is.na(value)]
  m[, station_id := substr(key, 1, 3)]
  out <- data.frame(location_id = .seuma_location(m$station_id, deployments), station_id = m$station_id,
                    month = m$month, value = m$value, stringsAsFactors = FALSE)
  out[order(out$station_id, out$month, method = "radix"), ]
}

# Invalid NO2 samples: missing or negative NO2; all eight channels exactly zero (logger failure);
# NO2 and NOx both zero. A zero NO2 with NOx positive or missing stays valid.
.no2_valid <- function(x) {
  p <- as.matrix(x[, seuma_pollutants, with = FALSE])
  missing <- !is.finite(x$no2)
  all_zero <- rowSums(is.finite(p) & p == 0) == length(seuma_pollutants)
  zero_nox <- is.finite(x$nox) & x$nox == 0
  !missing & !(x$no2 < 0) & !all_zero & !(x$no2 == 0 & zero_nox)
}

# Exact duplicate samples (same station and timestamp, same values) count once; when duplicates
# differ, the valid one is kept and two different valid ones stop the run.
.seuma_collapse_duplicates <- function(x) {
  x[, dup := .N > 1L, by = .(station_id, timestamp)]
  if (!any(x$dup)) return(x[, dup := NULL][])
  sig <- do.call(paste, c(lapply(seuma_pollutants, function(p) format(x[[p]], digits = 17)), sep = "|"))
  x[, sig := sig]
  x[, keep := TRUE]
  x[dup == TRUE, keep := {
    if (data.table::uniqueN(sig) == 1L) seq_len(.N) == 1L else {
      if (data.table::uniqueN(sig[valid]) > 1L) stop("conflicting valid duplicate samples")
      first <- if (any(valid)) which(valid)[1] else 1L
      seq_len(.N) == first
    }
  }, by = .(station_id, timestamp)]
  x[keep == TRUE][, c("dup", "sig", "keep") := NULL][]
}

.seuma_location <- function(station_id, deployments) {
  i <- match(station_id, deployments$station_id)
  if (anyNA(i)) stop("unknown SEUMA station: ", station_id[is.na(i)][1])
  deployments$location_id[i]
}

# Hourly PM2.5 of the reviewed workbook-months. An hour is valid when it has a finite value and
# no flag from the operator records (data/decisions/seuma-pm25-flagged-hours.csv); zeros stay.
# Value: data.table location_id, station_id, month, day, hour (YYYY-MM-DD HH), value, valid.
seuma_pm25_hours <- function(samples, deployments, months, flagged) {
  x <- data.table::copy(samples)[, .(station_id, month, timestamp, pm25)]
  x[, timestamp_local := format(timestamp, "%Y-%m-%d %H:%M:%S", tz = "America/Fortaleza")]
  f <- data.table::as.data.table(flagged)
  if (anyDuplicated(f[, .(station_id, timestamp_local)])) stop("duplicate flagged PM2.5 hour")
  hit <- f[x, on = .(station_id, timestamp_local), which = TRUE]
  if (length(setdiff(seq_len(nrow(f)), hit))) stop("a flagged PM2.5 hour matches no exported row")
  x[, flagged := !is.na(hit)]
  reviewed <- months[months$decision == "reviewed", ]
  x <- x[paste(station_id, month) %in% paste(reviewed$station_id, reviewed$month)]
  if (!all(grepl("^[0-9]{4}-[0-9]{2}-[0-9]{2} [0-9]{2}:00:(00|10)$", x$timestamp_local))) {
    stop("PM2.5 rows must be hourly logger records")
  }
  x[, `:=`(location_id = .seuma_location(station_id, deployments), hour = substr(timestamp_local, 1, 13),
           value = stored_value(pm25))]
  x[, dup := .N > 1L, by = .(location_id, hour)]
  if (any(x[dup == TRUE, .(n = data.table::uniqueN(.SD)), by = .(location_id, hour),
            .SDcols = c("timestamp_local", "value", "flagged")]$n > 1L)) {
    stop("conflicting duplicate PM2.5 hours")
  }
  x <- unique(x, by = c("location_id", "hour"))
  x[, valid := !flagged & is.finite(value)]
  x[, day := substr(hour, 1, 10)]
  x[order(location_id, hour), .(location_id, station_id, month, day, hour, value, valid)]
}

# PM2.5 monthly means of the reviewed months with coverage. Value: data.frame location_id, month, value.
seuma_pm25_monthly <- function(hours) {
  h <- data.table::copy(hours)[, key := location_id]
  m <- coverage_means(h, grouped = TRUE)[!is.na(value)]
  out <- data.frame(location_id = m$key, month = m$month, value = stored_value(m$value),
                    stringsAsFactors = FALSE)
  out[order(out$location_id, out$month, method = "radix"), ]
}
