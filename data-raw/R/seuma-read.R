# Reading the reference station workbooks published by the Fortaleza municipal government (SEUMA).
# Two layouts exist: "legacy" (2020-2021, sheet "Estação", DIA/HORA columns, gases in ppb) and
# "modern" (from 2021-11, one sheet, Excel serial datetime, ug/m3).

seuma_pollutants <- c("o3", "co", "no", "no2", "nox", "so2", "pm10", "pm25")

# Workbook table: data-raw/sources/seuma.csv joined with data/decisions/seuma-workbooks.csv.
# Value: one row per workbook with file, url, sha256, bytes, station_id, month, layout, unit_basis.
seuma_workbooks <- function(sources = "data-raw/sources/seuma.csv",
                            decisions = "data/decisions/seuma-workbooks.csv") {
  s <- read_table(sources)
  d <- read_table(decisions)
  if (anyDuplicated(s$file) || !setequal(s$file, d$file) || nrow(s) != nrow(d)) {
    stop("SEUMA sources and workbook decisions list different files")
  }
  w <- merge(s, d, by = "file", sort = FALSE)
  if (anyDuplicated(paste(w$station_id, w$month))) stop("one workbook per station and month")
  if (!all(w$layout %in% c("legacy", "modern"))) stop("unknown SEUMA layout")
  w[order(w$station_id, w$month, method = "radix"), ]
}

# Downloads every listed workbook from the municipal website into dir, then verifies them.
download_seuma <- function(workbooks, dir = "work/raw/seuma") {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  for (i in seq_len(nrow(workbooks))) {
    status <- utils::download.file(workbooks$url[i], file.path(dir, workbooks$file[i]),
                                   mode = "wb", quiet = TRUE)
    if (status != 0) stop("download failed: ", workbooks$url[i])
  }
  verify_seuma(workbooks, dir)
}

# Stops unless every listed workbook is in dir with the recorded size and SHA-256.
verify_seuma <- function(workbooks, dir = "work/raw/seuma") {
  path <- file.path(dir, workbooks$file)
  missing <- !file.exists(path)
  if (any(missing)) stop("SEUMA workbook missing: ", path[missing][1], " (run download_seuma())")
  size <- file.size(path)
  sha <- vapply(path, digest::digest, character(1), algo = "sha256", file = TRUE)
  bad <- size != workbooks$bytes | sha != workbooks$sha256
  if (any(bad)) stop("SEUMA workbook differs from its recorded SHA-256: ", path[bad][1])
  invisible(path)
}

# Parses one workbook. Value: data.table with station_id, month, row, timestamp (POSIXct,
# America/Fortaleza) and the eight pollutant columns as exported; attributes no2_unit, pm25_unit.
read_seuma_workbook <- function(path, station_id, month, layout, unit_basis = "header") {
  x <- switch(layout,
    legacy = .read_seuma_legacy(path, month),
    modern = .read_seuma_modern(path),
    stop("unknown SEUMA layout: ", layout))
  units <- c(no2 = .seuma_unit(x$units[["no2"]], layout == "legacy", unit_basis, path),
             pm25 = .seuma_unit(x$units[["pm25"]], FALSE, unit_basis, path))
  out <- data.table::data.table(station_id = station_id, month = month, row = x$rows,
                                timestamp = x$timestamp)
  for (p in seuma_pollutants) data.table::set(out, j = p, value = x$values[, p])
  data.table::setattr(out, "no2_unit", units[["no2"]])
  data.table::setattr(out, "pm25_unit", units[["pm25"]])
  out
}

# Reads all workbooks. Value: one data.table of exported rows, with columns no2_unit and pm25_unit.
read_seuma_samples <- function(workbooks, dir = "work/raw/seuma") {
  data.table::rbindlist(lapply(seq_len(nrow(workbooks)), function(i) {
    w <- workbooks[i, ]
    x <- read_seuma_workbook(file.path(dir, w$file), w$station_id, w$month, w$layout, w$unit_basis)
    x[, `:=`(no2_unit = attr(x, "no2_unit"), pm25_unit = attr(x, "pm25_unit"))]
  }))
}

.seuma_sheet <- function(path, sheet) {
  x <- readxl::read_excel(path, sheet = sheet, col_names = FALSE, col_types = "text",
                          .name_repair = "minimal")
  if (!nrow(x) || !ncol(x)) stop("empty SEUMA workbook: ", path)
  m <- unname(trimws(as.matrix(x)))
  m[!is.na(m) & !nzchar(m)] <- NA_character_
  width <- ncol(m)
  if (width < 18) m <- cbind(m, matrix(NA_character_, nrow(m), 18 - width))
  attr(m, "width") <- width
  m
}

.seuma_number <- function(text) suppressWarnings(as.numeric(text))

.seuma_anchor <- function(hit, what) {
  row <- which(hit)
  if (length(row) != 1L) stop(what, " header is missing or ambiguous")
  row
}

# Header unit, or the unit fixed by a documented decision when the header declares none.
.seuma_unit <- function(header, ppb, unit_basis, path) {
  expected <- if (ppb) "ppb" else "ug/m3"
  unit <- if (is.na(header)) NA_character_ else
    switch(header, ppb = "ppb", "µg/m³" = , "ug/m³" = , "µg/m3" = , "ug/m3" = "ug/m3", "unknown")
  if (is.na(unit)) {
    if (identical(unit_basis, "header")) stop("undeclared unit without a decision: ", path)
    return(expected)
  }
  if (unit != expected) stop("unexpected unit '", header, "' in ", path)
  unit
}

.read_seuma_legacy <- function(path, month) {
  sheets <- readxl::excel_sheets(path)
  if (sum(sheets == "Estação") != 1L) stop("legacy workbook needs sheet Estação: ", path)
  m <- .seuma_sheet(path, "Estação")
  h <- .seuma_anchor(m[, 1] %in% "DIA" & m[, 2] %in% "HORA" & m[, 7] %in% c("NO2", "NO₂"), "legacy")
  if (!identical(unname(m[h, 4:11]), c("O₃", "CO", "NO", "NO₂", "NOx", "SO₂", "PM10", "PM2,5")) ||
      !identical(m[h + 1, 1], "UNIDADES")) stop("legacy header does not match: ", path)
  rows <- seq.int(h + 2L, nrow(m))
  day <- .seuma_number(m[rows, 1])
  fraction <- .seuma_number(m[rows, 2])
  ok <- is.finite(day) & is.finite(fraction)
  .seuma_skipped_rows(m, rows, ok, "legacy", path)
  cols <- c(o3 = 4, co = 5, no = 6, no2 = 7, nox = 8, so2 = 9, pm10 = 10, pm25 = 11)
  list(rows = rows[ok], timestamp = .seuma_day_time(month, day[ok], fraction[ok]),
       values = .seuma_values(m, rows[ok], cols), units = c(no2 = m[h + 1, 7], pm25 = m[h + 1, 11]))
}

.seuma_modern_headers <- list(
  so2 = c("SO2", "SO2 (µg/m³)"), no = c("NO", "NO (µg/m³)"),
  no2 = c("NO2", "NO₂", "NO2 (µg/m³)", "NO2  (µg/m³)", "NO₂ (µg/m³)"),
  nox = c("NOx", "NOx (µg/m³)"), o3 = c("O3", "O₃", "O3 (µg/m³)", "O3  (µg/m³)"),
  pm10 = c("MP10", "MP10 (µg/m³)", "MP10  (µg/m³)"),
  pm25 = c("MP2,5", "MP2,5 (µg/m³)", "MP2,5  (µg/m³)"), co = c("CO", "CO (ppm)"))

.read_seuma_modern <- function(path) {
  if (length(readxl::excel_sheets(path)) != 1L) stop("modern workbook needs one sheet: ", path)
  m <- .seuma_sheet(path, 1L)
  h <- .seuma_anchor(m[, 1] %in% "MÊS" & m[, 4] %in% "Dióxido de Nitrogênio", "modern") + 1L
  ok_header <- mapply(function(v, choices) v %in% choices, m[h, 2:9], .seuma_modern_headers)
  if (!all(ok_header)) stop("modern header does not match: ", path)
  unit_of <- function(header) if (grepl("(µg/m³)", header, fixed = TRUE)) "µg/m³" else NA_character_
  rows <- seq.int(h + 1L, nrow(m))
  serial <- .seuma_number(m[rows, 1])
  ok <- is.finite(serial)
  .seuma_skipped_rows(m, rows, ok, "modern", path)
  cols <- c(o3 = 6, co = 9, no = 3, no2 = 4, nox = 5, so2 = 2, pm10 = 7, pm25 = 8)
  list(rows = rows[ok], timestamp = .seuma_excel_time(serial[ok]),
       values = .seuma_values(m, rows[ok], cols),
       units = c(no2 = unit_of(m[h, 4]), pm25 = unit_of(m[h, 8])))
}

.seuma_values <- function(m, rows, cols) {
  v <- vapply(cols, function(j) .seuma_number(m[rows, j]), numeric(length(rows)))
  if (!is.matrix(v)) v <- matrix(v, nrow = length(rows), dimnames = list(NULL, names(cols)))
  v[, seuma_pollutants, drop = FALSE]
}

# Rows without a timestamp may only be export summaries, empty rows or the legacy notes block.
.seuma_skipped_rows <- function(m, rows, ok, layout, path) {
  labels <- m[rows, 1]
  summary <- tolower(labels) %in% switch(layout,
    legacy = c("falha operação", "falha de operação", "força maior", "manutenção", "calibração",
               "dados válidos", "eficiência do dia", "dados para o mês:", "eficiência:", "observações:"),
    modern = c("média", "máximo", "mínimo", "minímo", "soma", "soma (acumulado)"))
  if (layout == "legacy") summary <- summary & is.na(m[rows, 2])
  last <- min(if (layout == "legacy") 18L else 16L, attr(m, "width"))
  empty <- is.na(m[rows, seq_len(last), drop = FALSE])
  blank <- rowSums(!empty) == 0
  notes <- rep(FALSE, length(rows))
  if (layout == "legacy") {
    notes <- cumsum(labels %in% "OBSERVAÇÕES:") > 0 & rowSums(!empty[, -1, drop = FALSE]) == 0
  }
  bad <- !ok & !summary & !blank & !notes
  if (any(bad)) stop("unsupported row ", rows[bad][1], " in ", path)
  invisible(TRUE)
}

.seuma_midnight <- function(dates) {
  as.POSIXct(paste(format(dates, "%Y-%m-%d"), "00:00:00"), format = "%Y-%m-%d %H:%M:%S",
             tz = "America/Fortaleza")
}

.seuma_excel_time <- function(serial) {
  whole <- floor(serial)
  dates <- as.Date(whole, origin = "1899-12-30")
  seconds <- round((serial - whole) * 86400)
  carry <- seconds >= 86400
  dates[carry] <- dates[carry] + 1L
  seconds[carry] <- seconds[carry] - 86400
  .seuma_midnight(dates) + seconds
}

.seuma_day_time <- function(month, day, fraction) {
  if (any(day != floor(day) | day < 1 | day > 31 | fraction < 0 | fraction >= 1)) {
    stop("invalid DIA/HORA in legacy workbook of ", month)
  }
  dates <- as.Date(sprintf("%s-%02d", month, as.integer(day)), format = "%Y-%m-%d")
  if (anyNA(dates) || any(format(dates, "%Y-%m") != month)) stop("legacy date outside ", month)
  .seuma_midnight(dates) + round(fraction * 86400)
}
