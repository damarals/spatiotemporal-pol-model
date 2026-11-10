root <- repo_root
for (f in c("write.R", "seuma-qc.R", "seuma-sites.R", "moqa.R")) source(file.path(root, "data-raw/R", f), local = TRUE)
suppressPackageStartupMessages(library(data.table))

test_that("the export reader keeps listed devices with well-formed timestamp labels", {
  x <- read_moqa_readings(testthat::test_path("fixtures", "moqa-export.csv"), "TEST01")
  expect_identical(x$label, c("2024-09-01T00:00:05Z", "2024-09-01T00:00:40Z", "2024-09-01T00:00:40Z"))
  expect_identical(x$pm25, c(4.5, 5, 5))
})

# Synthetic readings: `per_hour` readings 45 s apart in each listed hour, value = day of month.
readings_of <- function(device, days, hours = 0:23, per_hour = 40, month = "2024-09", value = NULL) {
  g <- CJ(day = days, hour = hours, k = seq_len(per_hour) - 1L)
  g[, label := sprintf("%s-%02dT%02d:%02d:%02dZ", month, day, hour, (k * 45L) %/% 60L, (k * 45L) %% 60L)]
  data.table(device_id = device, label = g$label, pm25 = if (is.null(value)) as.numeric(g$day) else value)
}
period_of <- function(device, cadence = 79, software = "intermediate", shift = 0, location = "site") {
  data.frame(period_id = paste0(device, "-p"), device_id = device, location_id = location,
             start_date = "2024-09-01", end_date_exclusive = "2024-10-01", cadence_s = cadence,
             software = software, clock_shift_h = shift, stringsAsFactors = FALSE)
}
no_decisions <- data.frame(device_id = character(), month = character(), decision = character())

test_that("an hour needs 75 % of the readings its cadence implies; days and months follow the coverage rule", {
  # 79 s cadence: 3600 / 79 = 45.6 readings per hour, 35 needed.
  r <- rbind(readings_of("A", 1:20, per_hour = 35), readings_of("A", 21, per_hour = 34),
             readings_of("A", 22, hours = 0:14, per_hour = 35))
  dm <- moqa_device_months(r, period_of("A"), no_decisions, 0.01)
  expect_identical(nrow(dm$hours[valid == TRUE]), 20L * 24L + 15L)
  m <- dm$months
  expect_identical(c(m$valid_days, m$required_days), c(20L, 20L))
  expect_identical(m$status, "admitted")
  expect_equal(m$value, mean(1:20), tolerance = 1e-14)
  short <- moqa_device_months(r[!grepl("-20T", label)], period_of("A"), no_decisions, 0.01)$months
  expect_identical(short$status, "insufficient_coverage")
})

test_that("a clock shift moves readings across hours, days and the period edge", {
  r <- readings_of("A", 1:21, per_hour = 35)
  m <- moqa_device_months(r, period_of("A", shift = -3), no_decisions, 0.01)
  expect_identical(min(m$hours$hour), "2024-09-01 00")
  expect_identical(m$months$valid_days, 21L)
})

test_that("zero shares and hand decisions withhold months; decisions name existing months", {
  r <- readings_of("A", 1:21, per_hour = 35)
  r[seq(1, .N, by = 50), pm25 := 0]
  expect_identical(moqa_device_months(r, period_of("A"), no_decisions, 0.01)$months$status, "zero_share")
  expect_identical(moqa_device_months(r, period_of("A"), no_decisions, 0.03)$months$status, "admitted")
  d <- data.frame(device_id = "A", month = "2024-09", decision = "unresolved_sensor_validity")
  expect_identical(moqa_device_months(r, period_of("A"), d, 0.03)$months$status, "unresolved_sensor_validity")
  d$month <- "2024-10"
  expect_error(moqa_device_months(r, period_of("A"), d, 0.03), "without readings")
})

test_that("periods stop on duplicates with different values, cadence mismatch and overlap", {
  r <- readings_of("A", 1, per_hour = 35)
  expect_identical(nrow(moqa_device_months(rbind(r, r[1]), period_of("A"), no_decisions, 1)$hours), 24L)
  expect_error(moqa_device_months(rbind(r, r[1][, pm25 := 99]), period_of("A"), no_decisions, 1), "conflicting")
  expect_error(moqa_device_months(r, period_of("A", cadence = 35), no_decisions, 1), "cadence")
  two <- rbind(period_of("A"), transform(period_of("A"), period_id = "A-q", start_date = "2024-09-15"))
  expect_error(moqa_device_months(r, two, no_decisions, 1), "overlapping")
})

months_of <- function(...) {
  data.table(device_id = c("C", "F", "F", "G"), location_id = c("ref", "f", "f", "g"),
             software = c("v4_5", "v4_5", "v4_5", "intermediate"), month = "2024-09",
             valid_days = 20L, required_days = 20L, zero_share = 0,
             status = c("admitted", "admitted", "zero_share", "admitted"), value = c(1, 2.123456789012345678, 3, 4))[...]
}

test_that("field months keep admitted months of one regime away from the reference locations", {
  f <- moqa_field_months(months_of(), "ref", "v4_5")
  expect_identical(f$device_id, "F")
  expect_identical(f$value, 2.12345678901235)
})

test_that("calibration pairs average both instruments over their common valid hours", {
  day <- function(d, h) sprintf("2024-09-%02d %02d", d, h)
  hours <- data.table(device_id = "C", period = "C-p", day = rep(c("2024-09-01", "2024-09-02"), each = 24),
                      hour = c(day(1, 0:23), day(2, 0:23)), value = c(1:24, 1:24) + 0.5, valid = TRUE)
  hours[c(1:4, 25:28), valid := FALSE]
  bam <- data.table(location_id = "ref", hour = c(day(1, 0:23), day(2, 0:23)), value = c(1:24, 1:24) * 2,
                    valid = c(rep(TRUE, 24), rep(c(TRUE, FALSE), 12)))
  bam[c(5:8, 29:32), valid := c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE)]
  dm <- list(months = months_of(), hours = hours)
  p <- moqa_calibration_pairs(dm, "ref", "v4_5", bam)
  expect_identical(p$day, "2024-09-01")
  expect_identical(p$common_hours, 20L)
  expect_identical(p$moqa, sum((5:24) + 0.5) / 20)
  expect_identical(p$bam, sum((5:24) * 2) / 20)
  expect_identical(p$location_id, "ref")
  expect_error(moqa_calibration_pairs(dm, c("ref", "f"), "v4_5", bam), "one colocated device")
})

test_that("low-cost sites are the admitted devices away from the reference locations", {
  registry <- data.table(device_id = c("G", "F", "C"), latitude = c(0, -3.7, -3.8), longitude = c(-39, -38.5, -38.6))
  s <- moqa_sites(months_of(), registry, "ref")
  expect_identical(s$device_id, c("F", "G"))
  expect_identical(s$location_id, c("f", "g"))
  expect_equal(c(s$x_m[2], s$y_m[2]), c(500000, 10000000), tolerance = 1e-9)
  expect_error(moqa_sites(months_of(), registry[device_id != "G"], "ref"), "registry")
})

test_that("the low-cost tables have the shape of the data contract", {
  dir <- file.path(root, "data/moqa")
  skip_without_moqa()
  m <- read_table(file.path(dir, "moqa-monthly.csv"))
  p <- read_table(file.path(dir, "calibration-pairs.csv"))
  expect_identical(names(m), c("device_id", "location_id", "month", "value"))
  expect_identical(names(p), c("device_id", "location_id", "day", "common_hours", "moqa", "bam"))
  expect_identical(c(nrow(m), nrow(p)), c(21L, 66L))
  expect_true(all(p$common_hours >= 16L))
  sites <- read_table(file.path(dir, "sites.csv"))
  expect_true(all(sites$source == "moqa"))
  expect_true(all(m$location_id %in% sites$location_id))
  covariates <- read_table(file.path(dir, "site-covariates.csv"))
  expect_identical(covariates$location_id, sites$location_id)
  expect_true(all(p$location_id %in% read_table(file.path(root, "data/decisions/seuma-locations.csv"))$location_id))
})
