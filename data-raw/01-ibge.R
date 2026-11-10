# IBGE municipal boundary of Fortaleza (2022) and study-area basemap layers (Ceara 2022, urban 2015).
# Downloads them with geobr into work/tmp/ibge/, compares them with the committed
# data/processed/boundary.gpkg and data/processed/basemap.gpkg and reports whether they are
# identical. The committed files are kept either way: a GeoPackage written again differs in its
# timestamps, so the basemap is compared layer by layer (attributes and geometry).
# Run from the repository root.
for (f in c("write.R", "ibge.R")) source(file.path("data-raw/R", f))

stage_dir <- "work/tmp/ibge"
dir.create(stage_dir, recursive = TRUE, showWarnings = FALSE)

committed <- "data/processed/boundary.gpkg"
ref <- sf::st_read(committed, layer = "fortaleza", quiet = TRUE)
new_b <- fetch_boundary()
cmp <- compare_boundary(new_b, ref)
staged_b <- write_boundary(new_b, file.path(stage_dir, "boundary.gpkg"))
same_bytes <- identical(readBin(staged_b, "raw", file.size(staged_b)),
                        readBin(committed, "raw", file.size(committed)))
cat(sprintf("boundary.gpkg: byte-identical=%s, ewkb_identical=%s, vertices new=%d ref=%d, area_km2 new=%.4f ref=%.4f, sym_diff_km2=%.4f (bytes new=%d ref=%d); keeping committed\n",
            same_bytes, cmp$ewkb_identical, cmp$vertices_new, cmp$vertices_ref,
            cmp$area_new_km2, cmp$area_ref_km2, cmp$sym_diff_km2,
            as.integer(file.size(staged_b)), as.integer(file.size(committed))))

committed_bm <- "data/processed/basemap.gpkg"
staged_bm <- write_basemap(fetch_basemap(), file.path(stage_dir, "basemap.gpkg"))
same_layer <- function(layer) {
  a <- sf::st_read(staged_bm, layer = layer, quiet = TRUE)
  b <- sf::st_read(committed_bm, layer = layer, quiet = TRUE)
  identical(sf::st_drop_geometry(a), sf::st_drop_geometry(b)) &&
    identical(sf::st_as_binary(sf::st_geometry(a)), sf::st_as_binary(sf::st_geometry(b)))
}
cat(sprintf("basemap.gpkg: municipalities identical=%s, urban_areas identical=%s; keeping committed\n",
            same_layer("municipalities"), same_layer("urban_areas")))
