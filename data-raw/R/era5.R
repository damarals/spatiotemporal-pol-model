# Open-Meteo ERA5 hourly meteorology (Fortaleza grid cell -3.75, -38.50),
# GAM (mgcv::bam) boundary-layer-height imputation for 2024-01..2024-06,
# and monthly aggregation over the model calendar (2020-01..2025-12).
#
# Note: The 8 hourly JSON files in data-raw/era5/ include 2026 hours up to
# 2026-10-31T23:00 because the GAM boundary-layer-height imputation model was
# fit on all complete hourly rows across 2019-12..2026-10 before subsetting the
# monthly output to the 72-month model calendar (2020-01..2025-12).

.era5_hourly_vars <- c(
  "temperature_2m", "relative_humidity_2m", "precipitation", "rain",
  "wind_speed_10m", "wind_speed_100m", "wind_direction_10m",
  "boundary_layer_height", "surface_pressure", "shortwave_radiation",
  "cloud_cover"
)

read_era5_hourly <- function(raw_dir = "data-raw/era5") {
  files <- sort(list.files(raw_dir, pattern = "^era5_hourly_.*[.]json$", full.names = TRUE))
  if (!length(files)) {
    stop("no ERA5 hourly JSON files found in ", raw_dir, call. = FALSE)
  }
  read_one <- function(f) {
    j <- jsonlite::fromJSON(f, simplifyVector = TRUE)
    if (!identical(j$latitude, -3.75) || !identical(j$longitude, -38.5) ||
        !identical(j$timezone, "America/Fortaleza")) {
      stop("unexpected ERA5 grid cell or timezone in ", f, call. = FALSE)
    }
    h <- j$hourly
    d <- data.frame(time = h$time, stringsAsFactors = FALSE)
    for (v in .era5_hourly_vars) d[[v]] <- as.numeric(h[[v]])
    d
  }
  hr <- do.call(rbind, lapply(files, read_one))
  if (anyDuplicated(hr$time)) {
    stop("duplicate timestamps in ERA5 hourly JSON files", call. = FALSE)
  }
  hr <- hr[order(hr$time), ]
  nonblh <- setdiff(.era5_hourly_vars, "boundary_layer_height")
  last_complete <- max(which(stats::complete.cases(hr[, nonblh])))
  hr <- hr[seq_len(last_complete), ]
  if (!all(stats::complete.cases(hr[, nonblh]))) {
    stop("missing non-BLH hourly ERA5 values inside complete window", call. = FALSE)
  }
  tt <- as.POSIXct(hr$time, format = "%Y-%m-%dT%H:%M", tz = "UTC")
  if (!all(diff(as.numeric(tt)) == 3600)) {
    stop("ERA5 hourly series is not strictly consecutive", call. = FALSE)
  }
  hr$hour <- as.integer(substr(hr$time, 12L, 13L))
  hr$doy <- as.integer(format(as.Date(substr(hr$time, 1L, 10L)), "%j"))
  rownames(hr) <- NULL
  hr
}

impute_era5_blh <- function(hr) {
  gap <- is.na(hr$boundary_layer_height)
  if (!any(gap)) {
    hr$blh_imputed <- FALSE
    return(hr)
  }
  blh_form <- log(boundary_layer_height) ~
    s(hour, bs = "cc", k = 12) +
    s(wind_speed_10m) +
    s(wind_speed_100m) +
    s(temperature_2m) +
    s(relative_humidity_2m) +
    s(shortwave_radiation) +
    s(cloud_cover) +
    s(precipitation, k = 5) +
    s(surface_pressure) +
    ti(hour, wind_speed_10m, bs = c("cc", "tp")) +
    ti(hour, temperature_2m, bs = c("cc", "tp")) +
    ti(hour, relative_humidity_2m, bs = c("cc", "tp")) +
    s(doy, bs = "cc", k = 8)

  fit <- mgcv::bam(
    blh_form,
    data = hr[!gap, ],
    discrete = TRUE,
    knots = list(hour = c(0, 24), doy = c(0, 366))
  )
  smear <- mean(exp(stats::residuals(fit)))
  pred <- as.numeric(exp(stats::predict(fit, newdata = hr[gap, ]))) * smear
  hr$blh_imputed <- gap
  hr$boundary_layer_height[gap] <- pred
  if (!all(stats::complete.cases(hr[, .era5_hourly_vars]))) {
    stop("ERA5 hourly series still contains missing values after BLH imputation", call. = FALSE)
  }
  hr
}

aggregate_era5_monthly <- function(hr,
                                   calendar = format(
                                     seq(as.Date("2020-01-01"), as.Date("2025-12-01"), by = "month"),
                                     "%Y-%m"
                                   )) {
  hr$month <- substr(hr$time, 1L, 7L)
  hr$vent_h <- hr$wind_speed_10m * hr$boundary_layer_height
  sp <- split(hr, hr$month)
  mo <- do.call(rbind, lapply(names(sp), function(m) {
    z <- sp[[m]]
    vent <- mean(z$vent_h)
    data.frame(
      month = m,
      n_hours = nrow(z),
      log_vent = log(vent),
      t2m = mean(z$temperature_2m),
      precip = sum(z$precipitation),
      stringsAsFactors = FALSE
    )
  }))
  days_in <- function(m) {
    d1 <- as.Date(paste0(m, "-01"))
    as.integer(seq(d1, by = "month", length.out = 2L)[2L] - d1)
  }
  expected <- vapply(mo$month, days_in, integer(1L)) * 24L
  mo <- mo[mo$n_hours == expected & mo$month %in% calendar, c("month", "log_vent", "t2m", "precip")]
  rownames(mo) <- NULL

  # In the original pipeline, met-build.R wrote met-monthly.csv via write.csv()
  # before export-current-inputs.R read it back with read.csv() and wrote
  # era5-monthly.csv with write_table().
  tmp <- tempfile(fileext = ".csv")
  on.exit(unlink(tmp, force = TRUE), add = TRUE)
  utils::write.csv(mo, tmp, row.names = FALSE)
  out <- utils::read.csv(tmp, stringsAsFactors = FALSE)
  if (!identical(out$month, calendar)) {
    stop("ERA5 monthly table does not match the 72-month model calendar", call. = FALSE)
  }
  out
}

compute_era5_monthly <- function(raw_dir = "data-raw/era5") {
  hr <- read_era5_hourly(raw_dir)
  hr <- impute_era5_blh(hr)
  aggregate_era5_monthly(hr)
}