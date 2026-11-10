# Municipal boundary of Fortaleza (IBGE 2022) and study-area basemap layers (geobr).

# Fetches the unsimplified municipal polygon of Fortaleza in EPSG:4674.
fetch_boundary <- function(code_muni = 2304400L, year = 2022L,
                           reader = geobr::read_municipality) {
  b <- reader(code_muni = code_muni, year = year, simplified = FALSE, showProgress = FALSE)
  if (!inherits(b, "sf") || nrow(b) != 1L ||
      !identical(as.character(b$code_muni[[1L]]), as.character(code_muni))) {
    stop("geobr did not return a single polygon for municipality ", code_muni, call. = FALSE)
  }
  b <- sf::st_make_valid(sf::st_transform(b, 4674))
  types <- as.character(sf::st_geometry_type(b))
  if (any(sf::st_is_empty(b)) || !all(sf::st_is_valid(b)) ||
      !all(types %in% c("POLYGON", "MULTIPOLYGON"))) {
    stop("invalid municipal boundary geometry", call. = FALSE)
  }
  b
}

# Writes the municipal boundary GeoPackage with layer "fortaleza".
write_boundary <- function(boundary, path = "data/processed/boundary.gpkg") {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  if (file.exists(path)) unlink(path)
  sf::st_write(boundary, path, layer = "fortaleza", quiet = TRUE, append = FALSE)
  invisible(path)
}

# Compares two municipal boundary sf objects on geometry, vertex count and EPSG:31984 area.
compare_boundary <- function(new, reference) {
  g_new <- sf::st_geometry(sf::st_transform(new, 4674))
  g_ref <- sf::st_geometry(sf::st_transform(reference, 4674))
  wkb_same <- identical(sf::st_as_binary(g_new, EWKB = TRUE)[[1L]],
                        sf::st_as_binary(g_ref, EWKB = TRUE)[[1L]])
  p_new <- sf::st_transform(g_new, 31984)
  p_ref <- sf::st_transform(g_ref, 31984)
  area_new <- as.numeric(sf::st_area(p_new))
  area_ref <- as.numeric(sf::st_area(p_ref))
  sym_diff <- if (wkb_same) 0 else sum(as.numeric(sf::st_area(sf::st_sym_difference(p_new, p_ref))))
  data.frame(
    ewkb_identical = wkb_same,
    vertices_new = nrow(sf::st_coordinates(g_new)),
    vertices_ref = nrow(sf::st_coordinates(g_ref)),
    area_new_km2 = area_new / 1e6,
    area_ref_km2 = area_ref / 1e6,
    sym_diff_km2 = sym_diff / 1e6,
    stringsAsFactors = FALSE
  )
}

# Fetches the Ceara 2022 municipalities and 2015 urbanised areas intersecting Ceara in EPSG:4674.
fetch_basemap <- function(code_state = 23L, muni_year = 2022L, urban_year = 2015L,
                          read_muni = geobr::read_municipality,
                          read_urban = geobr::read_urban_area) {
  munis <- read_muni(code_muni = code_state, year = muni_year, showProgress = FALSE)
  if (!inherits(munis, "sf") || !nrow(munis)) stop("failed to fetch state municipalities", call. = FALSE)
  urban <- read_urban(year = urban_year, showProgress = FALSE)
  if (!inherits(urban, "sf") || !nrow(urban)) stop("failed to fetch urban areas", call. = FALSE)
  munis <- sf::st_make_valid(sf::st_transform(munis, 4674))
  urban <- sf::st_make_valid(sf::st_transform(urban, 4674))
  urban <- urban[munis, , op = sf::st_intersects]
  if (!nrow(urban)) stop("no urban areas intersect the state municipalities", call. = FALSE)
  list(municipalities = munis, urban_areas = urban)
}

# Writes the basemap GeoPackage with layers "municipalities" and "urban_areas".
write_basemap <- function(basemap, path = "data/processed/basemap.gpkg") {
  if (!is.list(basemap) || !all(c("municipalities", "urban_areas") %in% names(basemap))) {
    stop("basemap must contain municipalities and urban_areas", call. = FALSE)
  }
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  if (file.exists(path)) unlink(path)
  sf::st_write(basemap$municipalities, path, layer = "municipalities", quiet = TRUE, append = FALSE)
  sf::st_write(basemap$urban_areas, path, layer = "urban_areas", quiet = TRUE, append = TRUE)
  invisible(path)
}