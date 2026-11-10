# Build data/processed/era5-monthly.csv (72 rows) from the 8 Open-Meteo ERA5
# hourly JSON files in data-raw/era5/, imputing the missing 2024-01..2024-06
# boundary-layer height with mgcv::bam and aggregating to the model calendar
# (2020-01..2025-12).
# Run from the repository root.

source("data-raw/R/write.R")
source("data-raw/R/era5.R")

stage_dir <- "work/tmp/era5"
dir.create(stage_dir, recursive = TRUE, showWarnings = FALSE)

era5 <- compute_era5_monthly("data-raw/era5")

staged <- file.path(stage_dir, "era5-monthly.csv")
committed <- "data/processed/era5-monthly.csv"
write_table(era5, staged)

raw_staged <- readBin(staged, "raw", file.info(staged)$size)
raw_committed <- readBin(committed, "raw", file.info(committed)$size)
if (identical(raw_staged, raw_committed)) {
  write_table(era5, committed)
  cat(sprintf("era5-monthly.csv: identical (%d rows, %d bytes)\n", nrow(era5), length(raw_staged)))
} else {
  cat(sprintf("era5-monthly.csv: DIFFERS (staged %d bytes vs committed %d bytes); keeping committed\n",
              length(raw_staged), length(raw_committed)))
}