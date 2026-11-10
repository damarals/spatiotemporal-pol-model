# Static spatial covariates at the monitoring sites (site-covariates.csv)
# and the 330 ACAG grid cells (grid-covariates.csv) from OpenStreetMap
# timestamped Overpass XML extracts, WorldPop 2020 UN-adjusted constrained
# population raster, and SatPM2.5 (ACAG V6.GL.03) cell-to-municipal log ratios.

.covariate_columns <- c(
  "n_bus_stops_250m",
  "n_bus_stops_500m",
  "len_mtps_100m",
  "len_mtps_500m",
  "len_vehicular_250m",
  "dist_mtp_road_m",
  "pop_1km",
  "acag_log_ratio_municipal"
)

.vehicular_highways <- c(
  "motorway", "trunk", "primary", "secondary", "tertiary",
  "motorway_link", "trunk_link", "primary_link", "secondary_link", "tertiary_link",
  "unclassified", "residential", "living_street", "service", "road", "busway", "track"
)

worldpop_sources <- function(path = "data-raw/sources/worldpop.csv") {
  read.csv(path, stringsAsFactors = FALSE)
}

osm_sources <- function(path = "data-raw/sources/osm.csv") {
  read.csv(path, stringsAsFactors = FALSE)
}

.as_projected_points <- function(df, id_col) {
  df$uid <- as.character(df[[id_col]])
  pts <- sf::st_as_sf(df, coords = c("x_m", "y_m"), crs = 31984, remove = FALSE)
  ll <- sf::st_coordinates(sf::st_transform(pts, 4326))
  pts$lon <- ll[, 1L]
  pts$lat <- ll[, 2L]
  pts
}

compute_acag_cell_log_ratio <- function(acag_cells_long) {
  # In the original scripts, cell-weights.csv was written by data.table::fwrite()
  # and read back via utils::read.csv() before computing cell and municipal means.
  f <- tempfile(fileext = ".csv")
  on.exit(unlink(f, force = TRUE), add = TRUE)
  data.table::fwrite(data.table::as.data.table(acag_cells_long), f)
  cw <- utils::read.csv(f, stringsAsFactors = FALSE)
  acell <- stats::aggregate(value ~ cell + lon + lat, data = cw, FUN = mean)
  mun_month <- vapply(
    split(cw, cw$month),
    function(d) sum(d$weight * d$value) / sum(d$weight),
    numeric(1L)
  )
  list(acell = acell, mun_period = mean(mun_month))
}

.acag_cov_at_points <- function(pts, acag_summary) {
  acell <- acag_summary$acell
  mun_period <- acag_summary$mun_period
  nearest_cell <- vapply(seq_len(nrow(pts)), function(i) {
    d <- pmax(abs(acell$lon - pts$lon[i]), abs(acell$lat - pts$lat[i]))
    j <- which.min(d)
    if (d[j] > 0.0050001) stop("point outside ACAG cell set: ", pts$uid[i], call. = FALSE)
    j
  }, 1L)
  cell_mean <- acell$value[nearest_cell]
  data.frame(
    uid = pts$uid,
    acag_log_ratio_municipal = log(cell_mean / mun_period),
    stringsAsFactors = FALSE
  )
}

.roads_cov_at_points <- function(pts, maj, loc, stops) {
  old_s2 <- sf::sf_use_s2()
  sf::sf_use_s2(FALSE)
  on.exit(sf::sf_use_s2(old_s2), add = TRUE)

  out <- data.frame(uid = pts$uid, stringsAsFactors = FALSE)
  mtp <- maj[maj$highway %in% c("motorway", "trunk", "primary"), ]
  dm <- sf::st_distance(pts, mtp)
  out$dist_mtp_road_m <- apply(dm, 1L, min)

  mtps <- maj[maj$highway %in% c(
    "motorway", "trunk", "primary", "secondary",
    "motorway_link", "trunk_link", "primary_link", "secondary_link"
  ), ]
  len_in <- function(lines, r) {
    buf <- sf::st_buffer(pts[, "uid"], r)
    ii <- suppressWarnings(sf::st_intersection(buf, lines[, "highway"]))
    o <- tapply(as.numeric(sf::st_length(ii)), ii$uid, sum)
    v <- unname(o[out$uid])
    v[is.na(v)] <- 0
    v
  }
  out$len_mtps_100m <- len_in(mtps, 100)
  out$len_mtps_500m <- len_in(mtps, 500)
  out$len_vehicular_250m <- len_in(loc[loc$highway %in% .vehicular_highways, ], 250)

  ds <- sf::st_distance(pts, stops)
  out$n_bus_stops_250m <- rowSums(units::drop_units(ds) <= 250)
  out$n_bus_stops_500m <- rowSums(units::drop_units(ds) <= 500)
  out
}

.pop_cov_at_points <- function(pts, pop_raster, zero_if_na = FALSE) {
  old_s2 <- sf::sf_use_s2()
  sf::sf_use_s2(FALSE)
  on.exit(sf::sf_use_s2(old_s2), add = TRUE)

  buf <- terra::vect(sf::st_transform(sf::st_buffer(pts[, "uid"], 1000), 4326))
  e <- terra::extract(pop_raster, buf, fun = sum, na.rm = TRUE, exact = TRUE)
  pop_1km <- as.numeric(e[, 2L])
  if (zero_if_na) pop_1km[is.na(pop_1km)] <- 0
  data.frame(uid = pts$uid, pop_1km = pop_1km, stringsAsFactors = FALSE)
}

.read_osm_tiles <- function(tile_files, layer, key) {
  parts <- lapply(tile_files, function(f) {
    tryCatch(sf::st_read(f, layer = layer, quiet = TRUE), error = function(e) NULL)
  })
  parts <- parts[!vapply(parts, is.null, TRUE)]
  parts <- parts[vapply(parts, nrow, 1L) > 0L]
  d <- do.call(rbind, parts)
  k <- do.call(paste, c(lapply(key, function(col) d[[col]]), sep = "|"))
  d <- d[!duplicated(k), ]
  sf::st_transform(d, 31984)
}

.signif_csv_roundtrip <- function(df) {
  num <- vapply(df, is.numeric, TRUE)
  df[num] <- lapply(df[num], function(v) signif(v, 7L))
  f <- tempfile(fileext = ".csv")
  on.exit(unlink(f, force = TRUE), add = TRUE)
  utils::write.csv(df, f, row.names = FALSE)
  utils::read.csv(f, stringsAsFactors = FALSE)
}

# sites_paths: site tables stacked in this order (the low-cost sites in data/moqa/, when
# present, then the reference locations); the site rows keep that order.
compute_covariate_tables <- function(sites_paths = "data/processed/sites.csv",
                                     acag_cells_path = "data/processed/acag-cells.csv",
                                     acag_cells_long = NULL,
                                     acag_raw_dir = "work/raw/acag",
                                     osm_dir = "work/raw/osm",
                                     worldpop_path = "work/raw/worldpop/bra_ppp_2020_UNadj_constrained.tif") {
  if (is.null(acag_cells_long)) {
    if (!exists("compute_acag_tables", mode = "function")) {
      source("data-raw/R/acag.R")
    }
    acag_cells_long <- compute_acag_tables(raw_dir = acag_raw_dir)$aggregate$cells
  }
  acag_summary <- compute_acag_cell_log_ratio(acag_cells_long)

  old_s2 <- sf::sf_use_s2()
  suppressMessages(sf::sf_use_s2(FALSE))
  on.exit(suppressMessages(sf::sf_use_s2(old_s2)), add = TRUE)

  sites_df <- do.call(rbind, lapply(sites_paths, read_table))
  cells_df <- read_table(acag_cells_path)
  pts_sites <- .as_projected_points(sites_df, "location_id")
  pts_cells <- .as_projected_points(cells_df, "cell")

  f1 <- list.files(osm_dir, "^q1-major-roads-.*\\.osm$", full.names = TRUE)
  f2 <- list.files(osm_dir, "^q2-site-buffers-.*\\.osm$", full.names = TRUE)
  f3 <- list.files(osm_dir, "^q3-bus-stops-.*\\.osm$", full.names = TRUE)
  f5 <- sort(list.files(file.path(osm_dir, "q5-tiles"), "^tile-[0-9]+-[0-9]+\\.osm$", full.names = TRUE))
  if (length(f1) != 1L || length(f2) != 1L || length(f3) != 1L || length(f5) != 16L) {
    stop("missing OSM extracts in ", osm_dir, call. = FALSE)
  }

  maj <- sf::st_transform(sf::st_read(f1, layer = "lines", quiet = TRUE), 31984)
  loc2 <- sf::st_transform(sf::st_read(f2, layer = "lines", quiet = TRUE), 31984)
  opt2 <- sf::st_transform(sf::st_read(f2, layer = "points", quiet = TRUE), 31984)
  stops2 <- opt2[!is.na(opt2$highway) & opt2$highway == "bus_stop", ]

  opt3 <- sf::st_transform(sf::st_read(f3, layer = "points", quiet = TRUE), 31984)
  stops3 <- opt3[!is.na(opt3$highway) & opt3$highway == "bus_stop", ]
  loc5 <- .read_osm_tiles(f5, "lines", "osm_id")

  pop_raster <- terra::crop(
    terra::rast(worldpop_path),
    terra::ext(-38.70, -38.35, -3.95, -3.62)
  )

  # Monitoring sites (using q1 + q2)
  s_acag <- .acag_cov_at_points(pts_sites, acag_summary)
  s_roads <- .roads_cov_at_points(pts_sites, maj, loc2, stops2)
  s_pop <- .pop_cov_at_points(pts_sites, pop_raster, zero_if_na = FALSE)
  s_merged <- Reduce(function(u, v) merge(u, v, by = "uid", sort = FALSE), list(s_acag, s_roads, s_pop))
  s_merged <- s_merged[match(pts_sites$uid, s_merged$uid), ]
  s_round <- .signif_csv_roundtrip(s_merged[, .covariate_columns])
  site_cov <- data.frame(
    location_id = sites_df$location_id,
    s_round[, .covariate_columns],
    stringsAsFactors = FALSE
  )
  rownames(site_cov) <- NULL

  # 330 ACAG cells (using q1 + q5 + q3)
  g_acag <- .acag_cov_at_points(pts_cells, acag_summary)
  g_roads <- .roads_cov_at_points(pts_cells, maj, loc5, stops3)
  g_pop <- .pop_cov_at_points(pts_cells, pop_raster, zero_if_na = TRUE)
  g_merged <- Reduce(function(u, v) merge(u, v, by = "uid", sort = FALSE), list(g_acag, g_roads, g_pop))
  g_merged <- g_merged[match(pts_cells$uid, g_merged$uid), ]
  g_round <- .signif_csv_roundtrip(g_merged[, .covariate_columns])
  grid_cov <- data.frame(
    cell = as.integer(cells_df$cell),
    g_round[, .covariate_columns],
    stringsAsFactors = FALSE
  )
  rownames(grid_cov) <- NULL

  list(sites = site_cov, grid = grid_cov)
}