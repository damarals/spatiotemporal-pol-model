# Reference station locations and the projected site table.

# Projects longitude/latitude (EPSG:4326) to SIRGAS 2000 / UTM 24S (EPSG:31984). Value: matrix X, Y.
project_31984 <- function(longitude, latitude) {
  points <- sf::st_as_sf(data.frame(longitude = longitude, latitude = latitude),
                         coords = c("longitude", "latitude"), crs = 4326)
  sf::st_coordinates(sf::st_transform(points, 31984))
}

# Reference site rows (one per location of data/decisions/seuma-locations.csv) ordered by location_id.
# Value: data.frame location_id, source, device_id (NA), x_m, y_m.
seuma_sites <- function(locations) {
  l <- locations[order(locations$location_id, method = "radix"), ]
  xy <- project_31984(l$longitude, l$latitude)
  data.frame(location_id = l$location_id, source = "seuma", device_id = NA_character_,
             x_m = xy[, 1], y_m = xy[, 2], stringsAsFactors = FALSE)
}
