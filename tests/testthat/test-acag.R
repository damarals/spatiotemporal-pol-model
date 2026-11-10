source(file.path(repo_root, "data-raw/R/write.R"), local = TRUE)
source(file.path(repo_root, "data-raw/R/ibge.R"), local = TRUE)
source(file.path(repo_root, "data-raw/R/tropomi.R"), local = TRUE)
source(file.path(repo_root, "data-raw/R/acag.R"), local = TRUE)
source(file.path(repo_root, "data-raw/R/era5.R"), local = TRUE)
source(file.path(repo_root, "data-raw/R/covariates.R"), local = TRUE)

acag_fixture_contract <- function() {
  default_acag_contract()
}

# The synthetic grid keeps one value per column of cells, so the area-weighted
# municipal mean is the column mean weighted by the covered fraction of each
# column and stays analytic even though the cell areas vary with latitude.
acag_fixture_column_values <- c(10, 20, 30, 80, 999, 999)
acag_fixture_expected_mean <- (0.5 * 10 + 20 + 30 + 0.5 * 80) / 3

.write_acag_fixture_netcdf <- function(path, column_values, variable, unit,
                                      resolution_degrees, include_crs,
                                      netcdf_month, product_version,
                                      missing_column) {
  longitude <- (seq_len(6L) - 0.5) * resolution_degrees
  latitude <- (seq_len(6L) - 1.5) * resolution_degrees
  coordinate_names <- if (include_crs) c("lon", "lat") else c("x", "y")
  coordinate_units <- if (include_crs) {
    c("degrees_east", "degrees_north")
  } else {
    c("", "")
  }
  dimensions <- list(
    ncdf4::ncdim_def(coordinate_names[[1L]], coordinate_units[[1L]], longitude),
    ncdf4::ncdim_def(coordinate_names[[2L]], coordinate_units[[2L]], latitude)
  )
  definition <- ncdf4::ncvar_def(
    variable, unit, dimensions,
    missval = -9999,
    longname = "Synthetic ACAG monthly PM2.5",
    prec = "float"
  )
  netcdf <- ncdf4::nc_create(path, definition, force_v4 = TRUE)
  values <- array(
    rep(column_values, times = length(latitude)),
    dim = c(length(longitude), length(latitude))
  )
  if (!is.null(missing_column)) values[missing_column, ] <- NA_real_
  ncdf4::ncvar_put(netcdf, variable, values)
  ncdf4::ncatt_put(
    netcdf, 0, "TIMECOVERAGE", gsub("-", "", netcdf_month, fixed = TRUE)
  )
  if (!is.null(product_version)) {
    ncdf4::ncatt_put(netcdf, 0, "PRODUCT_VERSION", product_version)
  }
  ncdf4::nc_close(netcdf)
  invisible(path)
}

make_acag_fixture <- function(months = c("2020-01", "2020-02"),
                              column_values = acag_fixture_column_values,
                              variable = "PM25",
                              unit = "ug/m3",
                              resolution_degrees = 0.01,
                              include_crs = TRUE,
                              netcdf_months = months,
                              product_version = "V6.GL.03",
                              missing_column = NULL,
                              missing_month = months[[length(months)]]) {
  root <- tempfile("acag-fixture-")
  directory <- file.path(root, "data/raw/fortaleza/acag-v6gl03")
  dir.create(directory, recursive = TRUE)

  entries <- lapply(seq_along(months), function(index) {
    month <- months[[index]]
    product_name <- .acag_product_name(month)
    relative <- file.path("data/raw/fortaleza/acag-v6gl03", product_name)
    full_path <- file.path(root, relative)
    .write_acag_fixture_netcdf(
      path = full_path,
      column_values = column_values,
      variable = variable,
      unit = unit,
      resolution_degrees = resolution_degrees,
      include_crs = include_crs,
      netcdf_month = netcdf_months[[index]],
      product_version = product_version,
      missing_column = if (identical(month, missing_month)) {
        missing_column
      } else {
        NULL
      }
    )
    start <- as.Date(paste0(month, "-01"))
    end <- seq(start, by = "month", length.out = 2L)[[2L]]
    data.table::data.table(
      manifest_id = "acag_manifest",
      id = paste0("acag:", month),
      path = relative,
      sha256 = digest::digest(file = full_path, algo = "sha256"),
      bytes = as.numeric(file.info(full_path)$size),
      source = .acag_source,
      product_identity = product_name,
      product_version = "V6.GL.03",
      period_start = paste0(format(start, "%Y-%m-%d"), "T00:00:00Z"),
      period_end = paste0(format(end, "%Y-%m-%d"), "T00:00:00Z"),
      month = month
    )
  })

  files <- data.table::rbindlist(entries)

  # The polygon halves the first and last covered column of cells and the first
  # and last row, so partial intersections carry a fractional weight.
  ring <- matrix(
    c(
      0.005, -0.005,
      0.035, -0.005,
      0.035, 0.045,
      0.005, 0.045,
      0.005, -0.005
    ),
    ncol = 2L,
    byrow = TRUE
  )
  boundary <- sf::st_sf(
    code_muni = "2304400",
    geometry = sf::st_sfc(sf::st_polygon(list(ring)), crs = 4674)
  )

  list(
    root = root,
    files = files,
    boundary = boundary,
    contract = acag_fixture_contract()
  )
}

acag_fixture_weights <- function(fixture) {
  compute_acag_boundary_weights(
    fixture$files, fixture$boundary, fixture$contract, fixture$root
  )
}

acag_fixture_aggregate <- function(fixture, weights = acag_fixture_weights(fixture)) {
  aggregate_acag_municipal(
    fixture$files, fixture$boundary, fixture$contract, weights, fixture$root
  )
}

testthat::test_that("the ACAG observation contract declares the municipal support", {
  contract <- acag_fixture_contract()

  testthat::expect_identical(contract$source, "acag")
  testthat::expect_identical(contract$product$unit, "ug/m3")
  testthat::expect_identical(contract$product$months, 60L)
  testthat::expect_identical(
    contract$aggregation$quantity,
    "municipal_area_weighted_arithmetic_mean"
  )
  testthat::expect_identical(
    contract$aggregation$support_type, "municipal_area_mean"
  )
  testthat::expect_identical(
    contract$aggregation$spatial_denominator, "municipal_area_ibge_2022"
  )
  testthat::expect_false(contract$coverage$renormalize_incomplete)
  testthat::expect_identical(
    contract$coverage$incomplete_status, "blocked_incomplete_coverage"
  )
})

testthat::test_that("partial cells carry the fraction of their intersected area", {
  fixture <- make_acag_fixture()
  weights <- acag_fixture_weights(fixture)
  cells <- weights$cells

  testthat::expect_identical(nrow(cells), 24L)
  testthat::expect_equal(sum(cells$weight), 1, tolerance = 1e-12)
  testthat::expect_equal(
    sum(cells$intersected_area_m2),
    weights$municipal_area_m2,
    tolerance = 1e-12
  )

  rows <- split(cells, round(cells$lat, 6L))
  testthat::expect_length(rows, 6L)
  for (row in rows) {
    row <- row[order(row$lon)]
    testthat::expect_identical(nrow(row), 4L)
    testthat::expect_equal(
      row$intersected_area_m2[[1L]] / row$intersected_area_m2[[2L]],
      0.5,
      tolerance = 1e-6
    )
    testthat::expect_equal(
      row$intersected_area_m2[[4L]] / row$intersected_area_m2[[3L]],
      0.5,
      tolerance = 1e-6
    )
    testthat::expect_equal(
      row$intersected_area_m2[[2L]],
      row$intersected_area_m2[[3L]],
      tolerance = 1e-6
    )
  }
})

testthat::test_that("the municipal value is the area-weighted arithmetic mean", {
  fixture <- make_acag_fixture()
  aggregate <- acag_fixture_aggregate(fixture)
  monthly <- aggregate$monthly

  testthat::expect_identical(monthly$month, c("2020-01", "2020-02"))
  testthat::expect_equal(
    monthly$value_ug_m3,
    rep(acag_fixture_expected_mean, 2L),
    tolerance = 1e-6
  )
  testthat::expect_false(isTRUE(all.equal(
    monthly$value_ug_m3[[1L]],
    mean(acag_fixture_column_values[1:4])
  )))
  testthat::expect_equal(monthly$coverage, rep(1, 2L), tolerance = 1e-12)
  testthat::expect_identical(monthly$n_cells, rep(24L, 2L))
  testthat::expect_identical(monthly$n_valid_cells, rep(24L, 2L))
  testthat::expect_identical(monthly$status, rep("approved", 2L))
})

testthat::test_that("a constant ACAG field is preserved exactly", {
  fixture <- make_acag_fixture(column_values = rep(7, 6L))
  aggregate <- acag_fixture_aggregate(fixture)

  testthat::expect_equal(
    aggregate$monthly$value_ug_m3, rep(7, 2L), tolerance = 1e-10
  )
})

testthat::test_that("incomplete coverage blocks the month instead of renormalizing", {
  fixture <- make_acag_fixture(missing_column = 2L)
  aggregate <- acag_fixture_aggregate(fixture)
  monthly <- aggregate$monthly

  testthat::expect_identical(
    monthly$status,
    c("approved", "blocked_incomplete_coverage")
  )
  testthat::expect_equal(monthly$coverage[[2L]], 2 / 3, tolerance = 1e-6)
  testthat::expect_identical(monthly$n_valid_cells, c(24L, 18L))
  testthat::expect_true(is.na(monthly$value_ug_m3[[2L]]))
  testthat::expect_equal(
    monthly$value_ug_m3[[1L]], acag_fixture_expected_mean, tolerance = 1e-6
  )
  renormalized <- (0.5 * 10 + 30 + 0.5 * 80) / 2
  testthat::expect_false(isTRUE(all.equal(
    monthly$coverage[[2L]] * renormalized, monthly$value_ug_m3[[2L]]
  )))
})

testthat::test_that("ACAG observations follow the canonical observation schema", {
  fixture <- make_acag_fixture(missing_column = 2L)
  aggregate <- acag_fixture_aggregate(fixture)
  observations <- aggregate$observations

  testthat::expect_s3_class(observations, "sf")
  testthat::expect_identical(names(observations), .acag_observation_columns)
  testthat::expect_identical(nrow(observations), 2L)
  testthat::expect_identical(
    observations$observation_id, c("acag:2020-01", "acag:2020-02")
  )
  testthat::expect_identical(observations$source, rep("acag", 2L))
  testthat::expect_identical(
    observations$support_type, rep("municipal_area_mean", 2L)
  )
  testthat::expect_identical(observations$unit, rep("ug/m3", 2L))
  testthat::expect_identical(observations$product_version, rep("V6.GL.03", 2L))
  testthat::expect_identical(
    observations$time_start, c("2020-01-01T00:00:00Z", "2020-02-01T00:00:00Z")
  )
  testthat::expect_identical(
    observations$time_end, c("2020-02-01T00:00:00Z", "2020-03-01T00:00:00Z")
  )
  testthat::expect_identical(
    observations$qc_status, c("approved", "blocked_incomplete_coverage")
  )
  testthat::expect_identical(
    observations$qc_reason_ids, c(NA_character_, "incomplete_municipal_coverage")
  )
  testthat::expect_true(all(is.na(observations$uncertainty)))
  testthat::expect_match(
    observations$location_version, "^ibge-2022-2304400:[0-9a-f]{64}$"
  )
  testthat::expect_identical(sf::st_crs(observations)$epsg, 4326L)
})

testthat::test_that("the traceable tables carry one row per month and cell", {
  fixture <- make_acag_fixture()
  aggregate <- acag_fixture_aggregate(fixture)

  testthat::expect_identical(names(aggregate$monthly), .acag_monthly_columns)
  testthat::expect_identical(names(aggregate$cells), .acag_cell_columns)
  testthat::expect_identical(nrow(aggregate$cells), 48L)
  testthat::expect_identical(
    sort(unique(aggregate$cells$month)), c("2020-01", "2020-02")
  )

  outputs <- write_acag_municipal_report(aggregate, fixture$root)
  testthat::expect_identical(basename(outputs), c(
    "monthly-observations.csv", "cell-weights.csv"
  ))
  testthat::expect_true(all(file.exists(outputs)))
  written <- data.table::fread(outputs[[1L]])
  testthat::expect_identical(names(written), .acag_monthly_columns)
  testthat::expect_identical(nrow(written), 2L)
})

testthat::test_that("ACAG NetCDF scientific metadata fail closed", {
  missing_variable <- make_acag_fixture(variable = "OTHER")
  testthat::expect_error(
    acag_fixture_weights(missing_variable),
    "ACAG NetCDF omits variable PM25"
  )

  wrong_unit <- make_acag_fixture(unit = "mg/m3")
  testthat::expect_error(
    acag_fixture_weights(wrong_unit),
    "ACAG PM25 unit must equal ug/m3"
  )

  wrong_resolution <- make_acag_fixture(resolution_degrees = 0.02)
  testthat::expect_error(
    acag_fixture_weights(wrong_resolution),
    "ACAG resolution must equal 0.01 degrees"
  )

  missing_crs <- make_acag_fixture(include_crs = FALSE)
  testthat::expect_error(
    acag_fixture_weights(missing_crs),
    "ACAG CRS must equal EPSG:4326"
  )

  wrong_month <- make_acag_fixture(netcdf_months = c("2020-03", "2020-02"))
  testthat::expect_error(
    acag_fixture_weights(wrong_month),
    "ACAG TIMECOVERAGE disagrees with manifest month"
  )

  wrong_version <- make_acag_fixture(product_version = "WRONG.VERSION")
  testthat::expect_error(
    acag_fixture_weights(wrong_version),
    "ACAG PRODUCT_VERSION disagrees with manifest version"
  )

  without_version <- make_acag_fixture(
    column_values = rep(7, 6L), product_version = NULL
  )
  testthat::expect_equal(
    acag_fixture_aggregate(without_version)$monthly$value_ug_m3,
    rep(7, 2L),
    tolerance = 1e-10
  )
})

testthat::test_that("ACAG payload and support fail closed", {
  fixture <- make_acag_fixture()
  tampered <- fixture
  tampered$files <- data.table::copy(fixture$files)
  data.table::set(tampered$files, 1L, "sha256", strrep("0", 64L))
  testthat::expect_error(
    acag_fixture_weights(tampered),
    "ACAG monthly checksum mismatch"
  )

  outside_ring <- matrix(
    c(0.005, 80, 0.035, 80, 0.035, 80.05, 0.005, 80.05, 0.005, 80),
    ncol = 2L,
    byrow = TRUE
  )
  outside <- fixture
  outside$boundary <- sf::st_sf(
    code_muni = "2304400",
    geometry = sf::st_sfc(sf::st_polygon(list(outside_ring)), crs = 4674)
  )
  testthat::expect_error(
    acag_fixture_weights(outside),
    "ACAG raster extent does not contain Fortaleza"
  )

  other_municipality <- fixture
  other_municipality$boundary <- sf::st_sf(
    code_muni = "3550308",
    geometry = sf::st_geometry(fixture$boundary)
  )
  testthat::expect_error(
    acag_fixture_weights(other_municipality),
    "ACAG requires one Fortaleza boundary with IBGE code 2304400"
  )

  weights <- acag_fixture_weights(fixture)
  moved <- weights
  moved$boundary_sha256 <- strrep("0", 64L)
  testthat::expect_error(
    acag_fixture_aggregate(fixture, moved),
    "ACAG boundary weights were built for another boundary"
  )
})

testthat::test_that("ACAG months are restricted to the frozen calendar", {
  testthat::expect_error(
    acag_fixture_weights(make_acag_fixture(months = c("2025-01", "2025-02"))),
    "ACAG month is outside the frozen calendar"
  )
  testthat::expect_identical(length(.acag_months), 60L)
  testthat::expect_identical(.acag_months[[1L]], "2020-01")
  testthat::expect_identical(.acag_months[[60L]], "2024-12")
})

testthat::test_that("IBGE, TROPOMI, ERA5, covariates and source manifests satisfy contracts", {
  # Source manifests
  src_acag <- acag_sources(file.path(repo_root, "data-raw/sources/acag.csv"))
  testthat::expect_identical(nrow(src_acag), 60L)
  testthat::expect_identical(names(src_acag), c("file", "url", "sha256", "bytes"))
  testthat::expect_true(all(grepl("^https://satpmdata[.]s3[.]amazonaws[.]com/", src_acag$url)))

  src_tropomi <- tropomi_sources(file.path(repo_root, "data-raw/sources/tropomi.csv"))
  testthat::expect_identical(nrow(src_tropomi), 2182L)
  testthat::expect_identical(names(src_tropomi), c("file", "url", "sha256", "bytes"))

  src_wp <- worldpop_sources(file.path(repo_root, "data-raw/sources/worldpop.csv"))
  testthat::expect_identical(nrow(src_wp), 1L)
  testthat::expect_identical(names(src_wp), c("file", "url", "sha256", "bytes"))

  src_osm <- osm_sources(file.path(repo_root, "data-raw/sources/osm.csv"))
  testthat::expect_identical(nrow(src_osm), 19L)
  testthat::expect_identical(names(src_osm), c("file", "url", "sha256", "bytes", "timestamp"))

  # IBGE boundary & basemap comparison
  b_path <- file.path(repo_root, "data/processed/boundary.gpkg")
  b_sf <- sf::st_read(b_path, quiet = TRUE)
  cmp_self <- compare_boundary(b_sf, b_sf)
  testthat::expect_true(cmp_self$ewkb_identical)
  testthat::expect_equal(cmp_self$sym_diff_km2, 0, tolerance = 1e-12)

  bm_path <- file.path(repo_root, "data/processed/basemap.gpkg")
  testthat::expect_true(file.exists(bm_path))
  bm_layers <- sf::st_layers(bm_path)$name
  testthat::expect_true(all(c("municipalities", "urban_areas") %in% bm_layers))

  # TROPOMI product parser, grid & monthly extraction
  parsed <- .parse_tropomi_product_name(
    "00000000-0000-0000-0000-000000000001",
    src_tropomi$file[[1L]]
  )
  testthat::expect_identical(parsed$timeliness, "RPRO")
  testthat::expect_identical(parsed$month, "2020-01")

  boundary <- sf::st_read(b_path, quiet = TRUE)
  grid <- .aligned_harp_grid(boundary, step = 0.035)
  testthat::expect_identical(length(grid$lat_edges), 8L)
  testthat::expect_identical(length(grid$lon_edges), 9L)
  ops <- .harp_operations(grid)
  testthat::expect_match(ops$per_product, "tropospheric_NO2_column_number_density_validity>75", fixed = TRUE)

  trop <- compute_tropomi_tables(
    nc_dir = file.path(repo_root, "data/processed/tropomi"),
    boundary_path = b_path
  )
  ref_tm <- read_table(file.path(repo_root, "data/processed/tropomi-monthly.csv"))
  ref_tc <- read_table(file.path(repo_root, "data/processed/tropomi-cells.csv"))
  testthat::expect_identical(trop$monthly, ref_tm)
  testthat::expect_identical(trop$cells, ref_tc)

  # ERA5 hourly reader
  hr <- read_era5_hourly(file.path(repo_root, "data-raw/era5"))
  testthat::expect_true(nrow(hr) > 60000L)
  testthat::expect_identical(sum(is.na(hr$boundary_layer_height)), 4368L)

  # Covariate ACAG log ratio helper
  synth_cells <- data.frame(
    month = rep(c("2020-01", "2020-02"), each = 2L),
    cell = rep(c(1L, 2L), times = 2L),
    lon = rep(c(-38.55, -38.54), times = 2L),
    lat = rep(c(-3.75, -3.75), times = 2L),
    weight = rep(c(0.5, 0.5), times = 2L),
    value = c(10, 30, 10, 30)
  )
  lr <- compute_acag_cell_log_ratio(synth_cells)
  testthat::expect_equal(lr$mun_period, 20, tolerance = 1e-12)
  pts_synth <- data.frame(uid = c("a", "b"), lon = c(-38.55, -38.54), lat = c(-3.75, -3.75))
  cov_lr <- .acag_cov_at_points(pts_synth, lr)
  testthat::expect_equal(cov_lr$acag_log_ratio_municipal, log(c(10 / 20, 30 / 20)), tolerance = 1e-12)
})