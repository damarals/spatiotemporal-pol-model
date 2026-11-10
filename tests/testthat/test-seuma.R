root <- repo_root
for (f in c("write.R", "seuma-read.R", "seuma-qc.R", "seuma-sites.R")) source(file.path(root, "data-raw/R", f), local = TRUE)
suppressPackageStartupMessages(library(data.table))
fixture <- function(name) testthat::test_path("fixtures", name)
plain <- function(df) { rownames(df) <- NULL; df }

expect_parsed <- function(x, name) {
  e <- read_table(fixture(paste0(name, "-expected.csv")))
  expect_identical(format(x$timestamp, "%Y-%m-%d %H:%M:%S", tz = "America/Fortaleza"), e$timestamp)
  for (p in seuma_pollutants) expect_identical(x[[p]], as.double(e[[p]]), info = p)
  expect_identical(attr(x, "no2_unit"), unique(e$no2_unit))
}

test_that("legacy workbook excerpt parses as the reference parser read the full workbook", {
  x <- read_seuma_workbook(fixture("seuma-legacy.xlsx"), "P00", "2021-01", "legacy")
  expect_parsed(x, "seuma-legacy")
  expect_identical(attr(x, "pm25_unit"), "ug/m3")
})

test_that("modern workbook excerpt parses and skips the export summary rows", {
  x <- read_seuma_workbook(fixture("seuma-modern.xlsx"), "P08", "2025-12", "modern")
  expect_parsed(x, "seuma-modern")
})

test_that("a header without units needs a documented unit decision", {
  path <- fixture("seuma-modern-no-unit.xlsx")
  expect_error(read_seuma_workbook(path, "P01", "2021-11", "modern"), "without a decision")
  x <- read_seuma_workbook(path, "P01", "2021-11", "modern", unit_basis = "operator report")
  expect_parsed(x, "seuma-modern-no-unit")
})

test_that("Excel serial datetimes round to the second and carry over midnight", {
  t <- .seuma_excel_time(c(44805.000115740739, 44805 + 86399.9 / 86400))
  expect_identical(format(t, "%Y-%m-%d %H:%M:%S", tz = "America/Fortaleza"),
                   c("2022-09-01 00:00:10", "2022-09-02 00:00:00"))
  expect_error(.seuma_day_time("2021-02", 30, 0), "outside")
})

hours_of <- function(group, day, n_valid, value = 1) {
  data.table(day = day, value = c(rep(value, n_valid), rep(99, 24 - n_valid)),
             valid = c(rep(TRUE, n_valid), rep(FALSE, 24 - n_valid)))[, key := group][]
}

test_that("a day needs 16 valid hours and a month two thirds of its days, rounded up", {
  feb <- rbindlist(c(lapply(1:19, function(d) hours_of("a", sprintf("2023-02-%02d", d), 16, d)),
                     list(hours_of("a", "2023-02-20", 15, 100))))
  apr <- rbindlist(lapply(1:19, function(d) hours_of("b", sprintf("2023-04-%02d", d), 24, d)))
  for (grouped in c(FALSE, TRUE)) {
    m <- coverage_means(rbind(feb, apr), grouped = grouped)
    expect_identical(m[key == "a"]$required_days, 19L)
    expect_identical(m[key == "a"]$valid_days, 19L)
    expect_equal(m[key == "a"]$value, mean(1:19), tolerance = 1e-14)
    expect_identical(m[key == "b"]$required_days, 20L)
    expect_true(is.na(m[key == "b"]$value))
  }
  expect_identical(as.integer(ceiling(2 * c(28, 29, 30, 31) / 3)), c(19L, 20L, 20L, 21L))
})

samples_of <- function(no2_values, nox_values = 1, others = 1, unit = "ug/m3") {
  x <- data.table(station_id = "P01", month = "2022-01", row = seq_along(no2_values),
                  timestamp = as.POSIXct("2022-01-01", tz = "America/Fortaleza") + 3600 * (seq_along(no2_values) - 1))
  for (p in seuma_pollutants) set(x, j = p, value = rep(others, length(no2_values)))
  x[, `:=`(no2 = no2_values, nox = nox_values, no2_unit = unit, pm25_unit = "ug/m3")]
  x[]
}

test_that("NO2 sample QC removes missing, negative, all-zero and zero NO2 with zero NOx", {
  x <- samples_of(c(NA, -1, 0, 0, 0, 0, 5), c(1, 1, 0, 0, 3, NA, 1))
  for (p in setdiff(seuma_pollutants, c("no2", "nox"))) set(x, i = 3L, j = p, value = 0)
  expect_identical(.no2_valid(x), c(FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE))
})

test_that("ppb NO2 converts with 46.0055 g/mol and 24.45 L/mol", {
  dep <- data.frame(station_id = "P01", location_id = "loc")
  x <- rbindlist(lapply(1:31, function(d) {
    s <- samples_of(rep(10, 24), unit = "ppb")
    s[, timestamp := timestamp + (d - 1) * 86400]
  }))
  out <- seuma_no2_monthly(x, dep)
  expect_identical(out$value, 10 * 46.0055 / 24.45)
})

test_that("exact duplicate samples count once and conflicting valid duplicates stop", {
  x <- samples_of(c(5, 5, 7))
  x[2, timestamp := x$timestamp[1]]
  x[, valid := .no2_valid(x)]
  expect_identical(nrow(.seuma_collapse_duplicates(copy(x))), 2L)
  y <- copy(x)[2, no2 := 6]
  expect_error(.seuma_collapse_duplicates(y), "conflicting")
  z <- copy(x)[2, `:=`(no2 = -1)][, valid := .no2_valid(.SD)]
  expect_identical(.seuma_collapse_duplicates(z)$no2, c(5, 7))
})

test_that("flagged PM2.5 hours are invalid, zeros stay valid and unmatched flags stop", {
  x <- samples_of(rep(1, 4))
  x[, timestamp := timestamp + 10]
  x[, pm25 := c(0, 12, NA, 8)]
  dep <- data.frame(station_id = "P01", location_id = "loc")
  months <- data.frame(station_id = "P01", month = "2022-01", decision = "reviewed")
  flag <- data.frame(station_id = "P01", timestamp_local = "2022-01-01 01:00:10")
  h <- seuma_pm25_hours(x, dep, months, flag)
  expect_identical(h$valid, c(TRUE, FALSE, FALSE, TRUE))
  expect_identical(h$hour, sprintf("2022-01-01 %02d", 0:3))
  flag$timestamp_local <- "2022-01-01 01:00:00"
  expect_error(seuma_pm25_hours(x, dep, months, flag), "matches no exported row")
  y <- rbind(x, x[2][, pm25 := 13])
  expect_error(seuma_pm25_hours(y, dep, months, flag[0, ]), "conflicting duplicate")
  months$decision <- "withheld"
  expect_identical(nrow(seuma_pm25_hours(x, dep, months, flag[0, ])), 0L)
})

test_that("stored values follow the 15-digit text written by data.table::fwrite", {
  expect_identical(stored_value(c(8.2200684693877548, NA, 2)), c(8.22006846938776, NA, 2))
  expect_false(identical(stored_value(8.2200684693877548), as.numeric(sprintf("%.15g", 8.2200684693877548))))
})

test_that("reference sites project to EPSG:31984", {
  xy <- project_31984(-39, 0)
  expect_equal(unname(xy[1, ]), c(500000, 10000000), tolerance = 1e-9)
  s <- seuma_sites(data.frame(location_id = c("b", "a"), latitude = c(-3.7, -3.8), longitude = c(-38.5, -38.6)))
  expect_identical(s$location_id, c("a", "b"))
  expect_true(all(is.na(s$device_id)) && all(s$source == "seuma"))
})

test_that("decision files agree with each other", {
  d <- function(name) read_table(file.path(root, "data/decisions", name))
  w <- seuma_workbooks(file.path(root, "data-raw/sources/seuma.csv"), file.path(root, "data/decisions/seuma-workbooks.csv"))
  expect_identical(nrow(w), 62L)
  expect_true(all(grepl("^[0-9a-f]{64}$", w$sha256)))
  dep <- d("seuma-deployments.csv")
  expect_identical(dep$station_id, sprintf("P%02d", 0:8))
  span <- aggregate(month ~ station_id, w, range)
  expect_identical(dep$first_month, span$month[, 1])
  expect_identical(dep$last_month, span$month[, 2])
  loc <- d("seuma-locations.csv")
  expect_setequal(unique(dep$location_id), loc$location_id)
  expect_setequal(d("seuma-site-classes.csv")$location_id, loc$location_id)
  pm <- d("seuma-pm25-months.csv")
  expect_identical(sort(paste(pm$station_id, pm$month)), sort(paste(w$station_id, w$month)))
  expect_identical(as.vector(table(pm$decision)[c("reviewed", "withheld")]), c(48L, 14L))
  fl <- d("seuma-pm25-flagged-hours.csv")
  expect_identical(nrow(fl), 958L)
  expect_true(all(fl$station_id %in% dep$station_id))
})

test_that("the workbooks rebuild the committed reference tables", {
  raw <- file.path(root, "work/raw/seuma")
  testthat::skip_if_not(dir.exists(raw) && length(list.files(raw, "[.]xlsx$")) == 62L, "SEUMA workbooks not downloaded")
  old <- setwd(root)
  on.exit(setwd(old))
  w <- seuma_workbooks()
  verify_seuma(w)
  dep <- read_table("data/decisions/seuma-deployments.csv")
  s <- read_seuma_samples(w)
  expect_identical(plain(seuma_no2_monthly(s, dep)), read_table("data/processed/seuma-no2-monthly.csv"))
  h <- seuma_pm25_hours(s, dep, read_table("data/decisions/seuma-pm25-months.csv"),
                        read_table("data/decisions/seuma-pm25-flagged-hours.csv"))
  expect_identical(plain(seuma_pm25_monthly(h)), read_table("data/processed/seuma-pm25-monthly.csv"))
})
