# TROPOMI Sentinel-5P L2 NO2 catalogue query, HARP conversion, selective GDAL
# staging, and extraction of monthly cell means and corners.
#
# Full L2 acquisition (2,182 orbits, ~1.01 TiB) requires CDSE S3 credentials
# (CDSE_S3_ACCESS_KEY, CDSE_S3_SECRET_KEY) and HARP 1.29 executables
# (harpconvert, harpmerge). Selective staging additionally uses GDAL
# gdalmdimtranslate over HTTP byte-range requests against the public AWS mirror
# (meeo-s5p.s3.amazonaws.com) or CDSE S3.
#
# For normal reproduction without downloading 1 TiB of L2 orbits, the 72
# monthly HARP NetCDF grids (201 KB total) are committed under
# data/processed/tropomi/, and compute_tropomi_tables() extracts the 1,584
# monthly cell observations and 88 projected cell corners directly from them.

.cdse_catalogue_url <- "https://catalogue.dataspace.copernicus.eu/odata/v1/Products"
.tropomi_product_prefixes <- c(
  RPRO = "S5P_RPRO_L2__NO2____",
  OFFL = "S5P_OFFL_L2__NO2____"
)
.tropomi_timeliness_precedence <- c("RPRO", "OFFL")
.tropomi_name_pattern <- paste0(
  "^S5P_(RPRO|OFFL)_L2__NO2____",
  "([0-9]{8}T[0-9]{6})_",
  "([0-9]{8}T[0-9]{6})_",
  "([0-9]{5})_([0-9]{2})_([0-9]{6})_",
  "([0-9]{8}T[0-9]{6})[.]nc$"
)
.tropomi_months <- format(
  seq(as.Date("2020-01-01"), as.Date("2025-12-01"), by = "month"),
  "%Y-%m"
)
.tropomi_validity_threshold <- 75L

tropomi_sources <- function(path = "data-raw/sources/tropomi.csv") {
  read.csv(path, stringsAsFactors = FALSE)
}

.format_tropomi_time <- function(value, label) {
  parsed <- as.POSIXct(value, format = "%Y%m%dT%H%M%S", tz = "UTC")
  if (is.na(parsed)) {
    stop("invalid ", label, " in TROPOMI product name", call. = FALSE)
  }
  format(parsed, "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
}

.parse_tropomi_product_name <- function(id = NA_character_, name) {
  if (length(name) != 1L || !is.character(name) || is.na(name)) {
    stop("TROPOMI product Name must be a single string", call. = FALSE)
  }
  match <- regmatches(name, regexec(.tropomi_name_pattern, name, perl = TRUE))[[1L]]
  if (length(match) != 8L) {
    stop("invalid Sentinel-5P L2 NO2 product name: ", name, call. = FALSE)
  }
  sensing_start <- .format_tropomi_time(match[[3L]], "sensing start")
  sensing_end <- .format_tropomi_time(match[[4L]], "sensing end")
  if (as.POSIXct(sensing_start, format = "%Y-%m-%dT%H:%M:%SZ", tz = "UTC") >=
      as.POSIXct(sensing_end, format = "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")) {
    stop("TROPOMI sensing start must precede sensing end", call. = FALSE)
  }
  data.frame(
    Id = as.character(id),
    Name = name,
    timeliness = match[[2L]],
    sensing_start = sensing_start,
    sensing_end = sensing_end,
    orbit = match[[5L]],
    collection = match[[6L]],
    processing_baseline = match[[7L]],
    creation_time = .format_tropomi_time(match[[8L]], "creation time"),
    logical_identity = paste0("S5P_L2__NO2____", match[[3L]], "_", match[[4L]], "_", match[[5L]]),
    month = substr(sensing_start, 1L, 7L),
    stringsAsFactors = FALSE
  )
}

.validate_cdse_products <- function(products) {
  products <- data.table::as.data.table(products)
  if ("content_length" %in% names(products) && !"ContentLength" %in% names(products)) {
    products[, ContentLength := content_length]
  }
  if (!all(c("Id", "Name") %in% names(products))) {
    stop("CDSE products must contain Id and Name", call. = FALSE)
  }
  if (!nrow(products)) {
    return(data.table::data.table(
      Id = character(), Name = character(), timeliness = character(),
      sensing_start = character(), sensing_end = character(), orbit = character(),
      collection = character(), processing_baseline = character(),
      creation_time = character(), logical_identity = character(), month = character()
    ))
  }
  parsed <- data.table::rbindlist(lapply(seq_len(nrow(products)), function(i) {
    data.table::as.data.table(.parse_tropomi_product_name(products$Id[[i]], products$Name[[i]]))
  }))
  if ("ContentLength" %in% names(products)) {
    parsed[, content_length := products$ContentLength]
  }
  identities <- unique(parsed$logical_identity)
  selected <- lapply(identities, function(identity) {
    candidates <- parsed[logical_identity == identity]
    ranks <- match(candidates$timeliness, .tropomi_timeliness_precedence)
    preferred <- candidates[timeliness == candidates$timeliness[[which.min(ranks)]]]
    if (nrow(preferred) != 1L) {
      stop("ambiguous processing baselines for ", identity, call. = FALSE)
    }
    preferred
  })
  selected <- data.table::rbindlist(selected)
  data.table::setorder(selected, sensing_start, sensing_end, orbit)
  selected[]
}

.tropomi_month_bounds <- function(month) {
  start <- as.Date(paste0(month, "-01"))
  next_month <- seq(start, by = "month", length.out = 2L)[[2L]]
  c(
    start = paste0(format(start - 1L, "%Y-%m-%d"), "T00:00:00.000Z"),
    end = paste0(format(next_month + 1L, "%Y-%m-%d"), "T00:00:00.000Z")
  )
}

.cdse_month_filter <- function(month, boundary_wkt) {
  bounds <- .tropomi_month_bounds(month)
  classes <- paste(sprintf("contains(Name,'%s')", unname(.tropomi_product_prefixes)), collapse = " or ")
  paste0(
    "Collection/Name eq 'SENTINEL-5P' and (", classes, ")",
    " and ContentDate/Start ge ", bounds[["start"]],
    " and ContentDate/Start lt ", bounds[["end"]],
    " and OData.CSC.Intersects(area=geography'SRID=4326;", boundary_wkt, "')"
  )
}

.query_cdse_month <- function(month, boundary_wkt, perform = httr2::req_perform) {
  filter <- .cdse_month_filter(month, boundary_wkt)
  query <- list(
    `$filter` = filter,
    `$select` = "Id,Name,ContentLength",
    `$orderby` = "ContentDate/Start asc,Name asc",
    `$top` = 1000L
  )
  req <- do.call(httr2::req_url_query, c(list(httr2::request(.cdse_catalogue_url)), query))
  rows <- list()
  repeat {
    resp <- perform(req)
    httr2::resp_check_status(resp)
    payload <- httr2::resp_body_json(resp, simplifyVector = FALSE)
    if (length(payload$value)) rows <- c(rows, lapply(payload$value, as.list))
    next_link <- payload[["@odata.nextLink"]]
    if (is.null(next_link)) next_link <- payload[["@OData.nextLink"]]
    if (is.null(next_link)) break
    req <- httr2::request(next_link)
  }
  if (!length(rows)) return(data.table::data.table())
  products <- data.table::rbindlist(rows, use.names = TRUE, fill = TRUE)
  validated <- .validate_cdse_products(products[, .(Id, Name, ContentLength)])
  validated[validated$month == month]
}

.cdse_s3_url <- function(product_name, sensing_start = NULL) {
  if (is.null(sensing_start)) {
    sensing_start <- .parse_tropomi_product_name(name = product_name)$sensing_start
  }
  sensing_day <- chartr("-", "/", substr(sensing_start, 1L, 10L))
  paste0(
    "https://eodata.dataspace.copernicus.eu/eodata/Sentinel-5P/TROPOMI/L2__NO2___/",
    sensing_day, "/", product_name
  )
}

.meeo_tropomi_url <- function(product_name) {
  parsed <- .parse_tropomi_product_name(name = product_name)
  sensing_day <- chartr("-", "/", substr(parsed$sensing_start, 1L, 10L))
  paste0(
    "https://meeo-s5p.s3.amazonaws.com/",
    parsed$timeliness, "/L2__NO2___/", sensing_day, "/", product_name
  )
}

download_cdse_product <- function(product_name, destination, expected_sha256 = NULL, expected_bytes = NULL) {
  access_key <- Sys.getenv("CDSE_S3_ACCESS_KEY", unset = "")
  secret_key <- Sys.getenv("CDSE_S3_SECRET_KEY", unset = "")
  if (!nzchar(access_key) || !nzchar(secret_key)) {
    stop("CDSE S3 credentials missing; set CDSE_S3_ACCESS_KEY and CDSE_S3_SECRET_KEY", call. = FALSE)
  }
  url <- .cdse_s3_url(product_name)
  req <- httr2::request(url) |>
    httr2::req_headers(`x-amz-content-sha256` = digest::digest("", algo = "sha256", serialize = FALSE)) |>
    httr2::req_auth_aws_v4(
      aws_access_key_id = access_key,
      aws_secret_access_key = secret_key,
      aws_service = "s3",
      aws_region = "us-east-1"
    ) |>
    httr2::req_retry(max_tries = 5L, retry_on_failure = TRUE) |>
    httr2::req_timeout(seconds = 7200)
  dir.create(dirname(destination), recursive = TRUE, showWarnings = FALSE)
  httr2::req_perform(req, path = destination)
  if (!is.null(expected_bytes)) {
    actual_bytes <- file.info(destination)$size
    if (!identical(as.numeric(actual_bytes), as.numeric(expected_bytes))) {
      stop("downloaded TROPOMI byte count mismatch for ", product_name, call. = FALSE)
    }
  }
  if (!is.null(expected_sha256)) {
    actual_sha <- digest::digest(file = destination, algo = "sha256")
    if (!identical(actual_sha, expected_sha256)) {
      stop("downloaded TROPOMI SHA-256 mismatch for ", product_name, call. = FALSE)
    }
  }
  invisible(destination)
}

.assert_fortaleza_boundary <- function(boundary) {
  if (!inherits(boundary, c("sf", "sfc"))) {
    stop("Fortaleza boundary must be an sf geometry", call. = FALSE)
  }
  geometry <- if (inherits(boundary, "sf")) sf::st_geometry(boundary) else boundary
  if (!length(geometry) || is.na(sf::st_crs(geometry))) {
    stop("Fortaleza boundary requires geometry and CRS", call. = FALSE)
  }
  geometry <- sf::st_transform(sf::st_make_valid(geometry), 4326)
  sf::st_union(geometry)
}

.aligned_harp_grid <- function(boundary, step = 0.035) {
  geometry <- .assert_fortaleza_boundary(boundary)
  bbox <- sf::st_bbox(geometry)
  tolerance <- 1e-10
  lower_index <- floor(c(bbox[["ymin"]], bbox[["xmin"]]) / step + tolerance)
  upper_index <- ceiling(c(bbox[["ymax"]], bbox[["xmax"]]) / step - tolerance)
  list(
    origin_degrees = 0L,
    step_degrees = step,
    lat_edges = round(seq.int(lower_index[[1L]], upper_index[[1L]]) * step, digits = 12L),
    lon_edges = round(seq.int(lower_index[[2L]], upper_index[[2L]]) * step, digits = 12L)
  )
}

.harp_operations <- function(grid) {
  format_edges <- function(edges) paste(sprintf("%.3f", edges), collapse = ",")
  variables <- paste0(
    "datetime_start,latitude,longitude,latitude_bounds,longitude_bounds,",
    "tropospheric_NO2_column_number_density,",
    "tropospheric_NO2_column_number_density_uncertainty"
  )
  per_product <- paste0(
    'set("propagate_uncertainty","correlated");',
    "tropospheric_NO2_column_number_density_validity>", .tropomi_validity_threshold, ";",
    "keep(", variables, ");",
    "bin_spatial((", format_edges(grid$lat_edges), "),(", format_edges(grid$lon_edges), "))"
  )
  list(
    per_product = per_product,
    monthly = 'set("propagate_uncertainty","correlated");bin()'
  )
}

.locate_tropomi_window <- function(latitude, longitude, qa_value, grid) {
  latitude_limits <- range(grid$lat_edges)
  longitude_limits <- range(grid$lon_edges)
  locator_lat <- latitude_limits + c(-1, 1)
  locator_lon <- longitude_limits + c(-1, 1)
  hit <- is.finite(latitude) & is.finite(longitude) &
    latitude >= locator_lat[[1L]] & latitude <= locator_lat[[2L]] &
    longitude >= locator_lon[[1L]] & longitude <= locator_lon[[2L]]
  indices <- which(hit, arr.ind = TRUE)
  if (nrow(indices)) {
    scanline <- range(indices[, "col"])
  } else {
    finite <- is.finite(latitude) & is.finite(longitude)
    survives <- finite & is.finite(qa_value) & round(100 * qa_value) > .tropomi_validity_threshold
    candidates <- if (any(survives)) survives else finite
    distance <- pmax(latitude_limits[[1L]] - latitude, 0, latitude - latitude_limits[[2L]])^2 +
      pmax(longitude_limits[[1L]] - longitude, 0, longitude - longitude_limits[[2L]])^2
    distance[!candidates] <- Inf
    nearest <- arrayInd(which.min(distance), dim(latitude))
    scanline <- rep(nearest[[2L]], 2L)
  }
  list(
    scanline = as.integer(c(max(1L, scanline[[1L]] - 1L), min(ncol(latitude), scanline[[2L]] + 1L)) - 1L),
    ground_pixel = as.integer(c(0L, nrow(latitude) - 1L))
  )
}

.read_harp_month <- function(path, grid) {
  netcdf <- ncdf4::nc_open(path)
  on.exit(ncdf4::nc_close(netcdf), add = TRUE)
  conventions <- ncdf4::ncatt_get(netcdf, 0, "Conventions")$value
  if (!is.character(conventions) || !grepl("HARP-1.0", conventions, fixed = TRUE)) {
    stop("NetCDF must declare HARP-1.0 conventions: ", path, call. = FALSE)
  }
  n_lat <- length(grid$lat_edges) - 1L
  n_lon <- length(grid$lon_edges) - 1L
  lat_bounds_raw <- ncdf4::ncvar_get(netcdf, "latitude_bounds", collapse_degen = FALSE)
  lon_bounds_raw <- ncdf4::ncvar_get(netcdf, "longitude_bounds", collapse_degen = FALSE)
  latitude_bounds <- t(matrix(lat_bounds_raw, nrow = 2L, ncol = n_lat))
  longitude_bounds <- t(matrix(lon_bounds_raw, nrow = 2L, ncol = n_lon))
  latitude <- rowMeans(latitude_bounds)
  longitude <- rowMeans(longitude_bounds)
  value_name <- "tropospheric_NO2_column_number_density"
  uncertainty_name <- paste0(value_name, "_uncertainty")
  values_raw <- ncdf4::ncvar_get(netcdf, value_name, collapse_degen = FALSE)
  uncertainties_raw <- ncdf4::ncvar_get(netcdf, uncertainty_name, collapse_degen = FALSE)
  values <- aperm(values_raw, c(3L, 2L, 1L))
  uncertainties <- aperm(uncertainties_raw, c(3L, 2L, 1L))
  list(
    axes = list(
      latitude = latitude,
      longitude = longitude,
      latitude_bounds = latitude_bounds,
      longitude_bounds = longitude_bounds
    ),
    values = values,
    uncertainties = uncertainties,
    conventions = conventions
  )
}

.harp_cells <- function(harp, boundary_4326) {
  cells <- data.table::CJ(
    lat_index = seq_along(harp$axes$latitude),
    lon_index = seq_along(harp$axes$longitude)
  )
  indices <- cbind(rep(1L, nrow(cells)), cells$lat_index, cells$lon_index)
  cells[, value := harp$values[indices]]
  cells[, uncertainty := harp$uncertainties[indices]]
  cells <- cells[!is.na(value)]
  if (!nrow(cells) || any(!is.finite(cells$value)) || any(!is.finite(cells$uncertainty)) || any(cells$uncertainty <= 0)) {
    stop("TROPOMI month has invalid or non-positive cells", call. = FALSE)
  }
  centers <- sf::st_sfc(lapply(seq_len(nrow(cells)), function(i) {
    sf::st_point(c(
      harp$axes$longitude[[cells$lon_index[[i]]]],
      harp$axes$latitude[[cells$lat_index[[i]]]]
    ))
  }), crs = 4326)
  inside <- lengths(sf::st_within(centers, boundary_4326)) == 1L
  cells[inside]
}

compute_tropomi_tables <- function(nc_dir = "data/processed/tropomi",
                                   boundary_path = "data/processed/boundary.gpkg",
                                   step = 0.035,
                                   unit_molm2 = 1e-5) {
  boundary <- sf::st_read(boundary_path, quiet = TRUE)
  boundary_4326 <- .assert_fortaleza_boundary(boundary)
  grid <- .aligned_harp_grid(boundary_4326, step = step)
  files <- sort(list.files(nc_dir, pattern = "^tropomi-no2-[0-9]{4}-[0-9]{2}\\.nc$", full.names = TRUE))
  if (length(files) != 72L) {
    stop("expected 72 monthly TROPOMI NetCDF files in ", nc_dir, ", found ", length(files), call. = FALSE)
  }
  monthly_list <- lapply(files, function(f) {
    m <- sub("^tropomi-no2-([0-9]{4}-[0-9]{2})\\.nc$", "\\1", basename(f))
    h <- .read_harp_month(f, grid)
    c <- .harp_cells(h, boundary_4326)
    u <- as.numeric(c$uncertainty)
    # In the model input preparation, Gaussian observation precision was scaled
    # by (1 / u^2) * unit_molm2^2 and converted back to internal SD via 1 / sqrt(scale).
    data.frame(
      cell_id = paste0(c$lat_index - 1L, ":", c$lon_index - 1L),
      month = m,
      value = as.numeric(c$value) / unit_molm2,
      sd = 1 / sqrt((1 / u^2) * unit_molm2^2),
      stringsAsFactors = FALSE
    )
  })
  monthly <- do.call(rbind, monthly_list)
  rownames(monthly) <- NULL

  cell_ids <- unique(monthly$cell_id)
  corners <- do.call(rbind, lapply(cell_ids, function(cid) {
    idx <- as.integer(strsplit(cid, ":", fixed = TRUE)[[1L]]) + 1L
    lat <- grid$lat_edges[idx[1L] + 0:1]
    lon <- grid$lon_edges[idx[2L] + 0:1]
    geographic <- data.frame(
      lon = lon[c(1L, 2L, 2L, 1L)],
      lat = lat[c(1L, 1L, 2L, 2L)]
    )
    projected <- sf::st_coordinates(sf::st_transform(
      sf::st_as_sf(geographic, coords = c("lon", "lat"), crs = 4326),
      31984
    ))
    data.frame(
      cell_id = cid,
      corner = c("SW", "SE", "NE", "NW"),
      x_m = as.numeric(projected[, 1L]),
      y_m = as.numeric(projected[, 2L]),
      stringsAsFactors = FALSE
    )
  }))
  rownames(corners) <- NULL
  list(monthly = monthly, cells = corners)
}