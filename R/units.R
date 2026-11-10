# Validation units and the fit to all data: what each withholds, its seed, and
# the partition of the inputs into training data and held-out targets.

background_location <- "fortaleza-seuma-hq-cajazeiras"

# held: locations or devices withheld (";"-separated), or the year predicted.
validation_units <- data.frame(
  unit = c("omitted-months", "new-traffic-1", "new-traffic-2", "new-traffic-3",
    "withheld-background", "year-2023", "year-2024", "year-2025",
    "new-devices-1", "new-devices-2", "new-devices-3", "all-data"),
  kind = c("months", "locations", "locations", "locations", "locations", "year", "year", "year",
    "devices", "devices", "devices", "none"),
  held = c("",
    "fortaleza-dhpp-fatima;fortaleza-terminal-parangaba",
    "fortaleza-bombeiros-moura-brasil;fortaleza-terminal-papicu",
    "fortaleza-sda-sao-gerardo;fortaleza-terminal-siqueira;fortaleza-terminal-messejana",
    background_location, "2023", "2024", "2025",
    "EDA13C;P1C1BJ;DD56F0", "DD2BA4;X3TYO5;Z3S59I", "H45UY3;DD2BDC;M1SSJN", ""),
  seed = c(111L, 101L, 102L, 103L, 104L, 123L, 124L, 125L, 107L, 108L, 109L, 126L),
  stringsAsFactors = FALSE)

# Omitted months: six NO2 and two PM2.5 months of 2024 at the background location.
omitted_no2_months <- c("2024-01", "2024-02", "2024-03", "2024-04", "2024-05", "2024-07")
omitted_pm25_months <- c("2024-08", "2024-11")

#' Partition the inputs for one unit.
#'
#' @param input value of read_inputs().
#' @param unit a name in validation_units$unit.
#' @return list: unit, seed, origin (first predicted month of a whole year, else
#'   NULL), train (reference rows), held (targets), moqa (low-cost training rows),
#'   calibration (pairs the fit may use), acag and column (rows available).
partition <- function(input, unit) {
  u <- validation_units[validation_units$unit == unit, , drop = FALSE]
  if (nrow(u) != 1L) stop("Unknown unit: ", unit, call. = FALSE)
  targets <- input$targets
  calibration <- input$calibration
  acag <- input$acag
  column <- input$column
  ids <- strsplit(u$held, ";", fixed = TRUE)[[1L]]
  origin <- NULL
  held <- switch(u$kind,
    months = targets$physical_location_id == background_location &
      ((targets$source == "ground_no2" & targets$month %in% omitted_no2_months) |
       (targets$source == "ground_pm25" & targets$month %in% omitted_pm25_months)),
    locations = targets$physical_location_id %in% ids,
    devices = targets$source == "moqa" & targets$device_id %in% ids,
    year = targets$physical_location_id == background_location & targets$source != "moqa" &
      substr(targets$month, 1L, 4L) == u$held,
    none = rep(FALSE, nrow(targets)))
  if (u$kind == "locations")
    calibration <- calibration[!calibration$location_id %in% ids, , drop = FALSE]
  if (u$kind == "year") {
    # Every source loses its rows from the predicted year on.
    origin <- paste0(u$held, "-01")
    keep <- targets$month < origin | held
    targets <- targets[keep, , drop = FALSE]
    held <- held[keep]
    calibration <- calibration[calibration$month < origin, , drop = FALSE]
    acag <- acag[acag$month < origin, , drop = FALSE]
    column <- column[column$month < origin, , drop = FALSE]
  }
  train <- targets[!held & targets$source != "moqa", , drop = FALSE]
  held <- targets[held, , drop = FALSE]
  # No daily BAM value of a held month enters the calibration.
  held_bam <- held[held$source == "ground_pm25", , drop = FALSE]
  calibration <- calibration[!paste(calibration$location_id, calibration$month) %in%
    paste(held_bam$physical_location_id, held_bam$month), , drop = FALSE]
  # Low-cost rows enter only where calibration pairs exist.
  moqa <- targets[targets$source == "moqa" & !targets$observation_id %in% held$observation_id &
    !targets$physical_location_id %in% held$physical_location_id[held$source != "moqa"], ,
    drop = FALSE]
  if (!nrow(calibration)) moqa <- moqa[0L, , drop = FALSE]
  if (any(held$source == "moqa") && !nrow(moqa))
    stop("Held low-cost rows need low-cost training rows and calibration pairs.", call. = FALSE)
  if (any(calibration$device_id %in% c(moqa$device_id, held$device_id)))
    stop("Calibration and deployment devices must be distinct.", call. = FALSE)
  strip <- function(x) {
    rownames(x) <- NULL
    x
  }
  list(unit = unit, seed = u$seed, origin = origin, train = strip(train), held = strip(held),
    moqa = strip(moqa), calibration = strip(calibration), acag = strip(acag), column = strip(column))
}
