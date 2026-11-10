# Build data/processed/site-covariates.csv (the 8 reference locations),
# data/moqa/site-covariates.csv (the 14 low-cost sites, when data/moqa/sites.csv
# is present) and data/processed/grid-covariates.csv (330 rows) from the
# timestamped Overpass OSM extracts in work/raw/osm/, the WorldPop 2020
# constrained 100 m raster in work/raw/worldpop/, and the ACAG V6.GL.03 cell
# values in work/raw/acag/ (or ACAG_RAW_DIR).
# Run from the repository root.

source("data-raw/R/write.R")
source("data-raw/R/acag.R")
source("data-raw/R/covariates.R")

stage_dir <- "work/tmp/covariates"
dir.create(stage_dir, recursive = TRUE, showWarnings = FALSE)

moqa_sites <- "data/moqa/sites.csv"
with_moqa <- file.exists(moqa_sites)
tables <- compute_covariate_tables(
  sites_paths = c(if (with_moqa) moqa_sites, "data/processed/sites.csv"),
  acag_cells_path = "data/processed/acag-cells.csv",
  acag_raw_dir = Sys.getenv("ACAG_RAW_DIR", unset = "work/raw/acag"),
  osm_dir = Sys.getenv("OSM_RAW_DIR", unset = "work/raw/osm"),
  worldpop_path = Sys.getenv(
    "WORLDPOP_RAW_PATH",
    unset = "work/raw/worldpop/bra_ppp_2020_UNadj_constrained.tif"
  )
)

stage_and_verify <- function(df, dir, name) {
  staged <- file.path(stage_dir, paste0(basename(dir), "-", name))
  committed <- file.path(dir, name)
  write_table(df, staged)
  raw_staged <- readBin(staged, "raw", file.info(staged)$size)
  raw_committed <- readBin(committed, "raw", file.info(committed)$size)
  if (identical(raw_staged, raw_committed)) {
    write_table(df, committed)
    cat(sprintf("%s: identical (%d rows, %d bytes)\n", committed, nrow(df), length(raw_staged)))
  } else {
    cat(sprintf("%s: DIFFERS (staged %d bytes vs committed %d bytes); keeping committed\n",
                committed, length(raw_staged), length(raw_committed)))
  }
}

public <- tables$sites$location_id %in% read_table("data/processed/sites.csv")$location_id
stage_and_verify(tables$sites[public, ], "data/processed", "site-covariates.csv")
if (with_moqa) stage_and_verify(tables$sites[!public, ], "data/moqa", "site-covariates.csv")
stage_and_verify(tables$grid, "data/processed", "grid-covariates.csv")
