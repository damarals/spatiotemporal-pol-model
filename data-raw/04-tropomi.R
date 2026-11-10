# Build data/processed/tropomi-monthly.csv (1,584 rows) and
# data/processed/tropomi-cells.csv (88 rows) from the 72 monthly HARP NetCDF
# files in data/processed/tropomi/ and the municipal polygon in
# data/processed/boundary.gpkg.
#
# Full L2 acquisition (1.01 TiB, 2,182 Sentinel-5P L2 orbits) is documented in
# data-raw/R/tropomi.R and data-raw/sources/tropomi.csv, and is not re-run here.
# Run from the repository root.

source("data-raw/R/write.R")
source("data-raw/R/tropomi.R")

stage_dir <- "work/tmp/tropomi"
dir.create(stage_dir, recursive = TRUE, showWarnings = FALSE)

tables <- compute_tropomi_tables(
  nc_dir = "data/processed/tropomi",
  boundary_path = "data/processed/boundary.gpkg"
)

stage_and_verify <- function(df, name) {
  staged <- file.path(stage_dir, name)
  committed <- file.path("data/processed", name)
  write_table(df, staged)
  raw_staged <- readBin(staged, "raw", file.info(staged)$size)
  raw_committed <- readBin(committed, "raw", file.info(committed)$size)
  if (identical(raw_staged, raw_committed)) {
    write_table(df, committed)
    cat(sprintf("%s: identical (%d rows, %d bytes)\n", name, nrow(df), length(raw_staged)))
  } else {
    cat(sprintf("%s: DIFFERS (staged %d bytes vs committed %d bytes); keeping committed\n",
                name, length(raw_staged), length(raw_committed)))
  }
}

stage_and_verify(tables$monthly, "tropomi-monthly.csv")
stage_and_verify(tables$cells, "tropomi-cells.csv")