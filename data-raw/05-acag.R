# Build data/processed/acag-monthly.csv (60 rows) and
# data/processed/acag-cells.csv (330 rows) from the 60 monthly ACAG V6.GL.03
# NetCDF files in work/raw/acag/ (or downloaded on demand from the anonymous
# S3 mirror listed in data-raw/sources/acag.csv) and the municipal polygon in
# data/processed/boundary.gpkg.
# Run from the repository root.

source("data-raw/R/write.R")
source("data-raw/R/acag.R")

raw_dir <- Sys.getenv(
  "ACAG_RAW_DIR",
  unset = "work/raw/acag"
)
src <- acag_sources("data-raw/sources/acag.csv")
missing_files <- !file.exists(file.path(raw_dir, src$file))
if (any(missing_files)) {
  download_acag(src[missing_files, , drop = FALSE], dest_dir = raw_dir)
}

stage_dir <- "work/tmp/acag"
dir.create(stage_dir, recursive = TRUE, showWarnings = FALSE)

tables <- compute_acag_tables(
  raw_dir = raw_dir,
  sources_path = "data-raw/sources/acag.csv",
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

stage_and_verify(tables$monthly, "acag-monthly.csv")
stage_and_verify(tables$cells, "acag-cells.csv")