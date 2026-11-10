# Low-cost PM2.5 network (MoQA): device-months, field months, calibration pairs and sites.
# The export and the operational workbook belong to the Fortaleza municipal government and are
# managed by TRAMA, UFC (see data/moqa/README.md).

# Readings of the listed devices. Value: data.table device_id, label (export Timestamp), pm25,
# in export order, with only well-formed timestamp labels.
read_moqa_readings <- function(path, devices) {
  x <- data.table::fread(path, select = c("moqaID", "Timestamp", "pm25"),
                         colClasses = list(character = c("moqaID", "Timestamp")),
                         na.strings = NULL, strip.white = FALSE, showProgress = FALSE)
  x <- x[moqaID %in% devices]
  data.table::setnames(x, c("moqaID", "Timestamp"), c("device_id", "label"))
  x[.moqa_valid_labels(label)]
}

.moqa_valid_labels <- function(labels) {
  dates <- unique(substr(labels, 1, 10))
  parsed <- as.Date(dates, format = "%Y-%m-%d")
  good <- dates[!is.na(parsed) & format(parsed, "%Y-%m-%d") == dates]
  !is.na(labels) & substr(labels, 1, 10) %in% good &
    grepl("^[0-9]{4}-[0-9]{2}-[0-9]{2}T([01][0-9]|2[0-3]):[0-5][0-9]:[0-5][0-9]Z$", labels)
}

# Device registry of the operational workbook (sheet "Rede de monitoramento").
# Value: data.table device_id, latitude, longitude (EPSG:4326).
read_moqa_registry <- function(path) {
  r <- as.data.frame(readxl::read_excel(path, sheet = "Rede de monitoramento", range = "A1:L32",
                                        col_names = FALSE, col_types = "text", trim_ws = FALSE,
                                        .name_repair = "minimal"))
  if (!identical(trimws(unname(unlist(r[1, 1:7]))),
                 c("Ponto", "MoQA_id", "Zona", "Latitude", "Longitude", "Bairro", "Referência"))) {
    stop("unexpected layout of the MoQA registry sheet")
  }
  r <- data.table::data.table(device_id = r[2:32, 2], latitude = as.numeric(r[2:32, 4]),
                              longitude = as.numeric(r[2:32, 5]))
  if (anyNA(r$device_id) || anyDuplicated(r$device_id)) stop("MoQA registry device IDs must be unique")
  r
}

# Hourly, daily and monthly means of every operational period (data/moqa/moqa-periods.csv).
# The export labels are local civil time despite their "Z", except where a period declares a clock
# shift. An hour is valid with at least 75 % of the readings its cadence implies; days and months
# follow coverage_means(). A month leaves when a hand decision names it (moqa-month-decisions.csv)
# or when more than zero_share_max of its readings are exactly zero.
# Value: list(hours = device_id, period, day, hour, value, valid;
#             months = device_id, location_id, software, month, valid_days, required_days,
#                      zero_share, status, value).
moqa_device_months <- function(readings, periods, decisions, zero_share_max) {
  .moqa_check_periods(periods)
  hours <- months <- vector("list", nrow(periods))
  for (i in seq_len(nrow(periods))) {
    p <- periods[i, ]
    start <- as.numeric(as.POSIXct(p$start_date, tz = "UTC"))
    end <- as.numeric(as.POSIXct(p$end_date_exclusive, tz = "UTC"))
    z <- readings[device_id == p$device_id]
    z[, seconds := as.numeric(as.POSIXct(label, format = "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")) +
        p$clock_shift_h * 3600]
    z <- z[seconds >= start & seconds < end]
    if (any(!is.finite(z$pm25))) stop("non-finite PM2.5 in period ", p$period_id)
    if (anyDuplicated(unique(z[, .(label, pm25)]), by = "label")) {
      stop("conflicting readings at one timestamp in period ", p$period_id)
    }
    z <- unique(z, by = "label")
    z[, hour := format(as.POSIXct(seconds, origin = "1970-01-01", tz = "UTC"), "%Y-%m-%d %H")]
    threshold <- ceiling(0.75 * 3600 / p$cadence_s)
    h <- z[, .(n = .N, value = mean(pm25), zeros = sum(pm25 == 0)), by = hour]
    h[, `:=`(valid = n >= threshold, day = substr(hour, 1, 10), key = p$period_id)]
    m <- coverage_means(h)
    zero <- h[, .(zero_share = sum(zeros) / sum(n)), by = .(month = substr(day, 1, 7))]
    m <- merge(m, zero, by = "month", sort = TRUE)
    m[, `:=`(device_id = p$device_id, location_id = p$location_id, software = p$software)]
    hours[[i]] <- h[, .(device_id = p$device_id, period = p$period_id, day, hour, value, valid)]
    months[[i]] <- m
  }
  months <- data.table::rbindlist(months)
  if (anyDuplicated(months, by = c("device_id", "month"))) stop("device-month in two periods")
  months[, status := data.table::fifelse(is.na(value), "insufficient_coverage", "admitted")]
  d <- data.table::as.data.table(decisions)
  i <- match(paste(d$device_id, d$month), paste(months$device_id, months$month))
  if (anyNA(i)) stop("a MoQA month decision names a month without readings")
  months[i, status := d$decision]
  months[status == "admitted" & zero_share > zero_share_max, status := "zero_share"]
  months <- months[, .(device_id, location_id, software, month, valid_days, required_days,
                       zero_share, status, value)]
  list(hours = data.table::rbindlist(hours), months = months[order(device_id, month)])
}

.moqa_check_periods <- function(periods) {
  cadences <- list(pre_v4_5 = c(27, 28), intermediate = c(79, 80, 81), v4_5 = 35)
  ok <- mapply(function(s, c) c %in% cadences[[s]], periods$software, periods$cadence_s)
  if (!all(ok)) stop("software regime disagrees with the cadence of its period")
  if (anyDuplicated(periods$period_id)) stop("duplicate MoQA period")
  for (id in unique(periods$device_id)) {
    q <- periods[periods$device_id == id, ]
    q <- q[order(q$start_date), ]
    if (any(q$start_date >= q$end_date_exclusive) ||
        any(utils::head(q$end_date_exclusive, -1) > utils::tail(q$start_date, -1))) {
      stop("overlapping or empty MoQA periods for ", id)
    }
  }
  invisible(TRUE)
}

# Field device-months: admitted, in the given software regime, away from the reference locations.
# Value: data.frame device_id, location_id, month, value (sensor scale, ug/m3).
moqa_field_months <- function(months, reference_locations, regime) {
  m <- months[status == "admitted" & software == regime & !location_id %in% reference_locations]
  data.frame(device_id = m$device_id, location_id = m$location_id, month = m$month,
             value = stored_value(m$value), stringsAsFactors = FALSE)
}

# Daily calibration pairs of the device colocated with the reference station, in its admitted months
# of the given regime: means of the sensor and of the reference over the hours valid for both, on
# days with at least 16 such hours. bam_hours: seuma_pm25_hours().
# Value: data.frame device_id, location_id, day, common_hours, moqa, bam.
moqa_calibration_pairs <- function(device_months, reference_locations, regime, bam_hours) {
  m <- device_months$months[status == "admitted" & software == regime &
                              location_id %in% reference_locations]
  if (data.table::uniqueN(m$device_id) != 1L) stop("expected one colocated device")
  h <- device_months$hours[device_id == m$device_id[1] & valid == TRUE & substr(day, 1, 7) %in% m$month]
  b <- data.table::as.data.table(bam_hours)[valid == TRUE & location_id == m$location_id[1]]
  common <- merge(h[, .(device_id, day, hour, moqa = value)], b[, .(hour, bam = value)], by = "hour")
  # Means of the stored hourly values (15 significant digits), summed to the nearest double.
  p <- common[, .(common_hours = .N, moqa = sum(stored_value(moqa)) / .N, bam = sum(bam) / .N),
              by = .(device_id, day)][common_hours >= 16L]
  data.frame(device_id = p$device_id, location_id = m$location_id[1], day = p$day,
             common_hours = p$common_hours, moqa = p$moqa, bam = p$bam,
             stringsAsFactors = FALSE)[order(p$day), ]
}

# Low-cost site rows: devices with at least one admitted month in any regime, away from the
# reference locations, at their registry coordinates in EPSG:31984.
# Value: data.frame location_id, source, device_id, x_m, y_m ordered by device_id.
moqa_sites <- function(months, registry, reference_locations) {
  d <- unique(months[status == "admitted" & !location_id %in% reference_locations,
                     .(device_id, location_id)])
  if (anyDuplicated(d$device_id)) stop("a MoQA device has two locations")
  d <- d[order(device_id, method = "radix")]
  r <- registry[match(d$device_id, registry$device_id)]
  if (anyNA(r$latitude) || anyNA(r$longitude)) stop("MoQA device without registry coordinates")
  xy <- project_31984(r$longitude, r$latitude)
  data.frame(location_id = d$location_id, source = "moqa", device_id = d$device_id,
             x_m = xy[, 1], y_m = xy[, 2], stringsAsFactors = FALSE)
}
