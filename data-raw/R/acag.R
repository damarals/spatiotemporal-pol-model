# SatPM2.5 (ACAG V6.GL.03) monthly 0.01-degree South America NetCDF acquisition,
# municipal cell-area weighting, and monthly aggregation over Fortaleza.

.acag_source <- "https://www.satpm.org/v6-gl-03"
.acag_months <- format(
  seq(as.Date("2020-01-01"), as.Date("2024-12-01"), by = "month"),
  "%Y-%m"
)

.acag_product_name <- function(month) {
  compact <- gsub("-", "", month, fixed = TRUE)
  paste0("V6GL03.CNNPM25.SA.", compact, "-", compact, ".nc")
}

.acag_s3_url <- function(file_or_month) {
  fname <- if (grepl("\\.nc$", file_or_month)) file_or_month else .acag_product_name(file_or_month)
  paste0("https://satpmdata.s3.amazonaws.com/V6GL03/CNNPM25/Monthly/SA/", fname)
}

acag_sources <- function(path = "data-raw/sources/acag.csv") {
  read.csv(path, stringsAsFactors = FALSE)
}

download_acag <- function(sources = acag_sources(), dest_dir = "work/raw/acag") {
  dir.create(dest_dir, recursive = TRUE, showWarnings = FALSE)
  paths <- character(nrow(sources))
  for (i in seq_len(nrow(sources))) {
    row <- sources[i, ]
    dest <- file.path(dest_dir, row$file)
    if (file.exists(dest)) {
      actual_sha <- digest::digest(file = dest, algo = "sha256")
      actual_bytes <- file.info(dest)$size
      if (identical(actual_sha, row$sha256) && identical(as.numeric(actual_bytes), as.numeric(row$bytes))) {
        paths[i] <- dest
        next
      }
    }
    req <- httr2::request(row$url) |>
      httr2::req_retry(max_tries = 5L, retry_on_failure = TRUE) |>
      httr2::req_timeout(seconds = 600)
    tmp <- paste0(dest, ".part")
    on.exit(if (file.exists(tmp)) unlink(tmp, force = TRUE), add = TRUE)
    httr2::req_perform(req, path = tmp)
    actual_bytes <- file.info(tmp)$size
    if (!identical(as.numeric(actual_bytes), as.numeric(row$bytes))) {
      stop("downloaded ACAG byte count mismatch for ", row$file, call. = FALSE)
    }
    actual_sha <- digest::digest(file = tmp, algo = "sha256")
    if (!identical(actual_sha, row$sha256)) {
      stop("downloaded ACAG SHA-256 mismatch for ", row$file, call. = FALSE)
    }
    if (!file.rename(tmp, dest)) {
      stop("failed to move downloaded ACAG file to ", dest, call. = FALSE)
    }
    paths[i] <- dest
  }
  invisible(paths)
}

default_acag_contract <- function() {
  list(
    schema_version = 1L,
    contract_id = "fortaleza_pm25_acag_observations",
    status = "development",
    region = "fortaleza",
    pollutant = "pm25",
    source = "acag",
    product = list(
      version = "V6.GL.03",
      url = .acag_source,
      variable = "PM25",
      unit = "ug/m3",
      resolution_degrees = 0.01,
      crs_epsg = 4326L,
      calendar_start = "2020-01",
      calendar_end = "2024-12",
      months = 60L
    ),
    aggregation = list(
      quantity = "municipal_area_weighted_arithmetic_mean",
      support_type = "municipal_area_mean",
      weight = "exact_spherical_cell_intersection_area",
      intersection = "exactextractr_coverage_area",
      spatial_denominator = "municipal_area_ibge_2022",
      boundary_ibge_code = "2304400",
      boundary_reference_year = 2022L
    ),
    coverage = list(
      minimum = 1,
      tolerance = 1e-6,
      renormalize_incomplete = FALSE,
      incomplete_status = "blocked_incomplete_coverage",
      approved_status = "approved"
    ),
    checks = list(
      area_cross_check_relative_tolerance = 0.005
    )
  )
}

.acag_monthly_columns <- c(
  "month", "status", "value_ug_m3", "coverage", "n_cells", "n_valid_cells",
  "intersected_area_m2", "municipal_area_m2", "boundary_area_m2",
  "time_start", "time_end", "payload_sha256", "product_version"
)

.acag_cell_columns <- c(
  "month", "cell", "lon", "lat", "intersected_area_m2", "weight", "value"
)

.acag_observation_columns <- c(
  "observation_id", "source", "source_product", "station_id",
  "physical_location_id", "month", "time_start", "time_end", "value",
  "unit", "support_type", "uncertainty", "qc_status", "qc_reason_ids",
  "location_version", "product_version", "geometry"
)

.acag_month_period <- function(month) {
  if (!is.character(month) || length(month) != 1L || is.na(month) || !month %in% .acag_months) {
    stop("ACAG month is outside the frozen calendar: ", month, call. = FALSE)
  }
  start <- as.Date(paste0(month, "-01"))
  end <- seq(start, by = "month", length.out = 2L)[[2L]]
  c(
    start = paste0(format(start, "%Y-%m-%d"), "T00:00:00Z"),
    end = paste0(format(end, "%Y-%m-%d"), "T00:00:00Z")
  )
}

.validate_acag_boundary <- function(boundary, contract = default_acag_contract()) {
  if (!inherits(boundary, "sf") || nrow(boundary) != 1L ||
      !"code_muni" %in% names(boundary) ||
      !identical(as.character(boundary$code_muni[[1L]]), contract$aggregation$boundary_ibge_code)) {
    stop("ACAG requires one Fortaleza boundary with IBGE code 2304400", call. = FALSE)
  }
  geometry <- sf::st_geometry(boundary)
  if (!length(geometry) || is.na(sf::st_crs(geometry)) || any(sf::st_is_empty(geometry))) {
    stop("Fortaleza boundary requires non-empty geometry and CRS", call. = FALSE)
  }
  sf::st_union(sf::st_transform(sf::st_make_valid(geometry), 4326))
}

.acag_boundary_sha256 <- function(boundary) {
  digest::digest(
    sf::st_as_binary(sf::st_geometry(boundary), EWKB = TRUE)[[1L]],
    algo = "sha256",
    serialize = FALSE
  )
}

.acag_netcdf_attribute <- function(netcdf, variable, name) {
  attribute <- ncdf4::ncatt_get(netcdf, variable, name)
  if (!isTRUE(attribute$hasatt)) return(NULL)
  attribute$value
}

.validate_acag_netcdf_metadata <- function(path, entry, contract = default_acag_contract()) {
  netcdf <- tryCatch(
    ncdf4::nc_open(path),
    error = function(error) {
      stop("ACAG NetCDF cannot be opened: ", conditionMessage(error), call. = FALSE)
    }
  )
  on.exit(ncdf4::nc_close(netcdf), add = TRUE)

  variable <- netcdf$var[[contract$product$variable]]
  if (is.null(variable)) {
    stop("ACAG NetCDF omits variable PM25", call. = FALSE)
  }
  if (!identical(variable$units, contract$product$unit)) {
    stop("ACAG PM25 unit must equal ug/m3", call. = FALSE)
  }
  dimension_names <- vapply(variable$dim, `[[`, character(1), "name")
  dimension_units <- vapply(variable$dim, `[[`, character(1), "units")
  if (!identical(dimension_names, c("lon", "lat")) ||
      !identical(dimension_units, c("degrees_east", "degrees_north"))) {
    stop("ACAG CRS must equal EPSG:4326", call. = FALSE)
  }
  valid_coordinate_grid <- vapply(variable$dim, function(dimension) {
    coordinates <- as.numeric(dimension$vals)
    if (length(coordinates) < 2L || any(!is.finite(coordinates))) return(FALSE)
    steps <- diff(coordinates)
    effective_resolution <- (coordinates[[length(coordinates)]] - coordinates[[1L]]) / (length(coordinates) - 1L)
    expected <- coordinates[[1L]] + (seq_along(coordinates) - 1L) * effective_resolution
    all(steps > 0) &&
      abs(effective_resolution - contract$product$resolution_degrees) <= 1e-6 &&
      max(abs(coordinates - expected)) <= 1e-5
  }, logical(1))
  if (any(!valid_coordinate_grid)) {
    stop("ACAG resolution must equal 0.01 degrees", call. = FALSE)
  }
  timecoverage <- .acag_netcdf_attribute(netcdf, 0, "TIMECOVERAGE")
  expected_timecoverage <- gsub("-", "", entry$month, fixed = TRUE)
  if (!identical(as.character(timecoverage), expected_timecoverage)) {
    stop("ACAG TIMECOVERAGE disagrees with manifest month", call. = FALSE)
  }
  product_version <- .acag_netcdf_attribute(netcdf, 0, "PRODUCT_VERSION")
  if (!is.null(product_version) && !identical(as.character(product_version), entry$product_version)) {
    stop("ACAG PRODUCT_VERSION disagrees with manifest version", call. = FALSE)
  }
  invisible(TRUE)
}

.acag_raster <- function(path, boundary, contract = default_acag_contract()) {
  raster <- tryCatch(
    terra::rast(path, subds = contract$product$variable, drivers = "NETCDF"),
    error = function(error) {
      stop("ACAG PM25 raster cannot be opened: ", conditionMessage(error), call. = FALSE)
    }
  )
  if (terra::nlyr(raster) != 1L) {
    stop("ACAG PM25 raster must contain exactly one layer", call. = FALSE)
  }
  resolution <- terra::res(raster)
  if (length(resolution) != 2L || any(!is.finite(resolution)) ||
      any(abs(resolution - contract$product$resolution_degrees) > 1e-6)) {
    stop("ACAG resolution must equal 0.01 degrees", call. = FALSE)
  }
  if (!isTRUE(sf::st_crs(terra::crs(raster)) == sf::st_crs(contract$product$crs_epsg))) {
    stop("ACAG CRS must equal EPSG:4326", call. = FALSE)
  }
  raster_extent <- terra::ext(raster)
  boundary_extent <- sf::st_bbox(boundary)
  tolerance <- 1e-8
  contains <- raster_extent$xmin <= boundary_extent[["xmin"]] + tolerance &&
    raster_extent$xmax >= boundary_extent[["xmax"]] - tolerance &&
    raster_extent$ymin <= boundary_extent[["ymin"]] + tolerance &&
    raster_extent$ymax >= boundary_extent[["ymax"]] - tolerance
  if (!isTRUE(contains)) {
    stop("ACAG raster extent does not contain Fortaleza", call. = FALSE)
  }
  raster
}

.acag_grid_identity <- function(raster) {
  list(
    nrow = terra::nrow(raster),
    ncol = terra::ncol(raster),
    extent = unname(as.vector(terra::ext(raster))),
    resolution = unname(terra::res(raster))
  )
}

.read_acag_month_payload <- function(files, index, boundary, contract = default_acag_contract(), root = ".") {
  row <- files[index, ]
  month <- as.character(row$month[[1L]])
  if (!month %in% .acag_months) {
    stop("ACAG month is outside the frozen calendar", call. = FALSE)
  }
  period <- .acag_month_period(month)
  entry <- list(
    id = if ("id" %in% names(row)) as.character(row$id[[1L]]) else paste0("acag:", month),
    path = as.character(row$path[[1L]]),
    sha256 = as.character(row$sha256[[1L]]),
    bytes = as.numeric(row$bytes[[1L]]),
    product_identity = if ("product_identity" %in% names(row)) as.character(row$product_identity[[1L]]) else basename(as.character(row$path[[1L]])),
    product_version = if ("product_version" %in% names(row)) as.character(row$product_version[[1L]]) else "V6.GL.03",
    period_start = if ("period_start" %in% names(row)) as.character(row$period_start[[1L]]) else period[["start"]],
    period_end = if ("period_end" %in% names(row)) as.character(row$period_end[[1L]]) else period[["end"]],
    month = month
  )
  path <- if (grepl("^(/|[A-Za-z]:)", entry$path)) entry$path else file.path(root, entry$path)
  if (!file.exists(path)) {
    stop("missing ACAG file: ", path, call. = FALSE)
  }
  payload_sha256 <- digest::digest(path, algo = "sha256", file = TRUE)
  if (!identical(entry$sha256, payload_sha256)) {
    stop("ACAG monthly checksum mismatch for ", entry$id, call. = FALSE)
  }
  if (!identical(as.numeric(entry$bytes), as.numeric(unname(file.info(path)$size)))) {
    stop("ACAG monthly byte count mismatch for ", entry$id, call. = FALSE)
  }
  .validate_acag_netcdf_metadata(path, entry, contract)
  list(
    entry = entry,
    payload_sha256 = payload_sha256,
    raster = .acag_raster(path, boundary, contract)
  )
}

compute_acag_boundary_weights <- function(files, boundary, contract = default_acag_contract(), root = ".") {
  files <- data.table::as.data.table(files)
  if (!nrow(files)) {
    stop("ACAG weights require verified monthly inputs", call. = FALSE)
  }
  for (m in as.character(files$month)) {
    if (!m %in% .acag_months) stop("ACAG month is outside the frozen calendar", call. = FALSE)
  }
  boundary <- .validate_acag_boundary(boundary, contract)
  payload <- .read_acag_month_payload(files, 1L, boundary, contract, root)
  cells <- exactextractr::exact_extract(
    payload$raster,
    boundary,
    coverage_area = TRUE,
    include_cell = TRUE,
    include_xy = TRUE,
    progress = FALSE
  )[[1L]]
  if (!is.data.frame(cells) || !all(c("cell", "x", "y", "coverage_area") %in% names(cells))) {
    stop("ACAG cell extraction returned an unexpected structure", call. = FALSE)
  }
  if (!nrow(cells) || anyDuplicated(cells$cell)) {
    stop("ACAG cells must intersect Fortaleza exactly once", call. = FALSE)
  }
  area <- as.numeric(cells$coverage_area)
  if (any(!is.finite(area)) || any(area <= 0)) {
    stop("ACAG intersected cell areas must be finite and positive", call. = FALSE)
  }
  municipal_area_m2 <- sum(area)
  boundary_area_m2 <- as.numeric(sf::st_area(boundary))
  relative_difference <- abs(municipal_area_m2 - boundary_area_m2) / boundary_area_m2
  if (!is.finite(relative_difference) ||
      relative_difference > contract$checks$area_cross_check_relative_tolerance) {
    stop("ACAG intersected area disagrees with the municipal boundary area", call. = FALSE)
  }
  weights <- data.table::data.table(
    cell = as.integer(cells$cell),
    lon = as.numeric(cells$x),
    lat = as.numeric(cells$y),
    intersected_area_m2 = area,
    weight = area / municipal_area_m2
  )
  data.table::setorder(weights, cell)
  list(
    cells = weights[],
    municipal_area_m2 = municipal_area_m2,
    boundary_area_m2 = boundary_area_m2,
    grid = .acag_grid_identity(payload$raster),
    boundary_sha256 = .acag_boundary_sha256(boundary)
  )
}

.validate_acag_weights <- function(weights, boundary, contract = default_acag_contract()) {
  if (!is.list(weights) || !identical(
    names(weights),
    c("cells", "municipal_area_m2", "boundary_area_m2", "grid", "boundary_sha256")
  )) {
    stop("ACAG boundary weights use an invalid structure", call. = FALSE)
  }
  cells <- weights$cells
  if (!data.table::is.data.table(cells) || !nrow(cells) ||
      !identical(names(cells), c("cell", "lon", "lat", "intersected_area_m2", "weight"))) {
    stop("ACAG boundary weights use an unexpected cell table", call. = FALSE)
  }
  if (!identical(weights$boundary_sha256, .acag_boundary_sha256(boundary))) {
    stop("ACAG boundary weights were built for another boundary", call. = FALSE)
  }
  if (abs(sum(cells$intersected_area_m2) - weights$municipal_area_m2) > 1e-6 * weights$municipal_area_m2 ||
      abs(sum(cells$weight) - 1) > 1e-9) {
    stop("ACAG cell weights do not sum to the municipal area", call. = FALSE)
  }
  invisible(TRUE)
}

.aggregate_acag_month <- function(files, index, boundary, contract, weights, root) {
  payload <- .read_acag_month_payload(files, index, boundary, contract, root)
  if (!identical(.acag_grid_identity(payload$raster), weights$grid)) {
    stop("ACAG month grid disagrees with the weighted municipal grid", call. = FALSE)
  }
  cells <- weights$cells
  values <- terra::extract(payload$raster, cells$cell)
  if (!is.data.frame(values) || ncol(values) != 1L || nrow(values) != nrow(cells)) {
    stop("ACAG value extraction returned an unexpected structure", call. = FALSE)
  }
  values <- as.numeric(values[[1L]])
  valid <- is.finite(values)
  coverage <- sum(cells$intersected_area_m2[valid]) / weights$municipal_area_m2
  blocked <- coverage < contract$coverage$minimum - contract$coverage$tolerance
  value <- if (blocked) {
    NA_real_
  } else {
    sum(cells$intersected_area_m2[valid] * values[valid]) / weights$municipal_area_m2
  }
  if (!blocked && (!is.finite(value) || value <= 0)) {
    stop("ACAG municipal mean must be finite and positive for ", payload$entry$month, call. = FALSE)
  }
  entry <- payload$entry
  list(
    monthly = data.table::data.table(
      month = entry$month,
      status = if (blocked) contract$coverage$incomplete_status else contract$coverage$approved_status,
      value_ug_m3 = value,
      coverage = coverage,
      n_cells = nrow(cells),
      n_valid_cells = sum(valid),
      intersected_area_m2 = sum(cells$intersected_area_m2[valid]),
      municipal_area_m2 = weights$municipal_area_m2,
      boundary_area_m2 = weights$boundary_area_m2,
      time_start = entry$period_start,
      time_end = entry$period_end,
      payload_sha256 = payload$payload_sha256,
      product_version = entry$product_version
    ),
    cells = data.table::data.table(
      month = entry$month,
      cell = cells$cell,
      lon = cells$lon,
      lat = cells$lat,
      intersected_area_m2 = cells$intersected_area_m2,
      weight = cells$weight,
      value = values
    ),
    product_identity = entry$product_identity
  )
}

aggregate_acag_municipal <- function(files, boundary, contract = default_acag_contract(), weights, root = ".") {
  files <- data.table::as.data.table(files)
  if (!nrow(files)) {
    stop("ACAG aggregation requires verified monthly inputs", call. = FALSE)
  }
  boundary <- .validate_acag_boundary(boundary, contract)
  .validate_acag_weights(weights, boundary, contract)
  months <- lapply(
    seq_len(nrow(files)),
    .aggregate_acag_month,
    files = files,
    boundary = boundary,
    contract = contract,
    weights = weights,
    root = root
  )
  monthly <- data.table::rbindlist(lapply(months, `[[`, "monthly"))
  cells <- data.table::rbindlist(lapply(months, `[[`, "cells"))
  if (!identical(names(monthly), .acag_monthly_columns) ||
      !identical(names(cells), .acag_cell_columns)) {
    stop("ACAG outputs drifted from the declared column contract", call. = FALSE)
  }
  if (anyDuplicated(monthly$month)) {
    stop("ACAG months must be unique", call. = FALSE)
  }
  observations <- sf::st_sf(
    data.frame(
      observation_id = paste0("acag:", monthly$month),
      source = contract$source,
      source_product = vapply(months, `[[`, character(1), "product_identity"),
      station_id = NA_character_,
      physical_location_id = NA_character_,
      month = monthly$month,
      time_start = monthly$time_start,
      time_end = monthly$time_end,
      value = monthly$value_ug_m3,
      unit = contract$product$unit,
      support_type = contract$aggregation$support_type,
      uncertainty = NA_real_,
      qc_status = monthly$status,
      qc_reason_ids = data.table::fifelse(
        monthly$status == contract$coverage$incomplete_status,
        "incomplete_municipal_coverage",
        NA_character_
      ),
      location_version = paste0(
        "ibge-", contract$aggregation$boundary_reference_year, "-",
        contract$aggregation$boundary_ibge_code, ":",
        weights$boundary_sha256
      ),
      product_version = monthly$product_version,
      stringsAsFactors = FALSE
    ),
    geometry = sf::st_sfc(rep(list(boundary[[1L]]), nrow(monthly)), crs = sf::st_crs(boundary))
  )
  observations <- observations[, .acag_observation_columns]
  list(observations = observations, monthly = monthly[], cells = cells[])
}

write_acag_municipal_report <- function(aggregate, root = ".") {
  if (!is.list(aggregate) || !identical(names(aggregate), c("observations", "monthly", "cells"))) {
    stop("ACAG aggregate has invalid fields", call. = FALSE)
  }
  relative <- file.path(
    "data/processed/fortaleza/acag",
    c("monthly-observations.csv", "cell-weights.csv")
  )
  outputs <- file.path(root, relative)
  dir.create(dirname(outputs[[1L]]), recursive = TRUE, showWarnings = FALSE)
  data.table::fwrite(aggregate$monthly, outputs[[1L]])
  data.table::fwrite(aggregate$cells, outputs[[2L]])
  outputs
}

.acag_fwrite_roundtrip <- function(x) {
  f <- tempfile(fileext = ".csv")
  on.exit(unlink(f, force = TRUE), add = TRUE)
  data.table::fwrite(data.table::data.table(x = x), f)
  data.table::fread(f)$x
}

compute_acag_tables <- function(raw_dir = "work/raw/acag",
                                sources_path = "data-raw/sources/acag.csv",
                                boundary_path = "data/processed/boundary.gpkg") {
  src <- acag_sources(sources_path)
  months <- sub("^V6GL03\\.CNNPM25\\.SA\\.([0-9]{4})([0-9]{2})-.*$", "\\1-\\2", src$file)
  files <- data.table::data.table(
    id = paste0("acag:", months),
    path = file.path(raw_dir, src$file),
    sha256 = src$sha256,
    bytes = src$bytes,
    product_identity = src$file,
    product_version = "V6.GL.03",
    month = months
  )
  boundary <- sf::st_read(boundary_path, quiet = TRUE)
  contract <- default_acag_contract()
  weights <- compute_acag_boundary_weights(files, boundary, contract, root = ".")
  agg <- aggregate_acag_municipal(files, boundary, contract, weights, root = ".")

  # In the original pipeline, write_acag_municipal_report() serialized the
  # monthly and cell tables with data.table::fwrite() before downstream model
  # preparation read them back and projected (lon, lat) to EPSG:31984 (x_m, y_m).
  monthly <- data.frame(
    month = as.character(agg$monthly$month),
    value = .acag_fwrite_roundtrip(as.numeric(agg$monthly$value_ug_m3)),
    stringsAsFactors = FALSE
  )
  w_cells <- weights$cells
  lon_rt <- .acag_fwrite_roundtrip(as.numeric(w_cells$lon))
  lat_rt <- .acag_fwrite_roundtrip(as.numeric(w_cells$lat))
  weight_rt <- .acag_fwrite_roundtrip(as.numeric(w_cells$weight))
  pts <- sf::st_transform(
    sf::st_as_sf(data.frame(lon = lon_rt, lat = lat_rt), coords = c("lon", "lat"), crs = 4326),
    31984
  )
  xy <- sf::st_coordinates(pts)
  cells_df <- data.frame(
    cell = as.integer(w_cells$cell),
    lon = lon_rt,
    lat = lat_rt,
    weight = weight_rt,
    x_m = as.numeric(xy[, 1L]),
    y_m = as.numeric(xy[, 2L]),
    stringsAsFactors = FALSE
  )
  list(monthly = monthly, cells = cells_df, aggregate = agg, weights = weights)
}