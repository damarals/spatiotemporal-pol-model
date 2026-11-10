# Low-cost network (MoQA): field device-months, calibration pairs and the site table.
# Usage: Rscript data-raw/03-moqa.R [output root, default "."]
# Reads work/raw/moqa/ (export and operational workbook, available from TRAMA on request) and
# work/raw/seuma/; applies the decisions in data/moqa/; writes data/moqa/moqa-monthly.csv,
# data/moqa/calibration-pairs.csv, the low-cost sites in data/moqa/sites.csv and the reference
# locations in data/processed/sites.csv.
for (f in c("write.R", "seuma-read.R", "seuma-qc.R", "seuma-sites.R", "moqa.R")) {
  source(file.path("data-raw/R", f))
}
suppressPackageStartupMessages(library(data.table))
out <- c(commandArgs(trailingOnly = TRUE), ".")[1]
raw <- data.frame(
  path = c("work/raw/moqa/moqa_data.csv", "work/raw/moqa/moqa_sensor_metadata.xlsx"),
  sha256 = c("51659871eb71c090b885bc9be7b52d3a7e00b8b6b08f1a97017cc01ec842fa6e",
             "8330577a81000c64b43f26acc27b1d08487ec53c166bea27cbf5a702f07eb87d"))
for (i in seq_len(nrow(raw))) {
  if (digest::digest(raw$path[i], algo = "sha256", file = TRUE) != raw$sha256[i]) {
    stop("unexpected content of ", raw$path[i])
  }
}
# The field months and the calibration pairs share one sensor software regime.
regime <- "v4_5"

periods <- read_table("data/moqa/moqa-periods.csv")
qc <- read_table("data/moqa/moqa-qc.csv")
registry <- read_moqa_registry(raw$path[2])
if (!all(periods$device_id %in% registry$device_id)) stop("period device missing from the registry")
readings <- read_moqa_readings(raw$path[1], unique(periods$device_id))
dm <- moqa_device_months(readings, periods, read_table("data/moqa/moqa-month-decisions.csv"),
                         qc$value[qc$parameter == "zero_share_max"])
rm(readings)

locations <- read_table("data/decisions/seuma-locations.csv")
deployments <- read_table("data/decisions/seuma-deployments.csv")
workbooks <- seuma_workbooks()
verify_seuma(workbooks)
bam <- seuma_pm25_hours(read_seuma_samples(workbooks), deployments,
                        read_table("data/decisions/seuma-pm25-months.csv"),
                        read_table("data/decisions/seuma-pm25-flagged-hours.csv"))

field <- moqa_field_months(dm$months, locations$location_id, regime)
pairs <- moqa_calibration_pairs(dm, locations$location_id, regime, bam)
sites <- moqa_sites(dm$months, registry, locations$location_id)

dir.create(file.path(out, "data/moqa"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out, "data/processed"), recursive = TRUE, showWarnings = FALSE)
write_table(field, file.path(out, "data/moqa/moqa-monthly.csv"))
write_table(pairs, file.path(out, "data/moqa/calibration-pairs.csv"))
write_table(sites, file.path(out, "data/moqa/sites.csv"))
write_table(seuma_sites(locations), file.path(out, "data/processed/sites.csv"))
cat(nrow(field), "field device-months,", nrow(pairs), "calibration pairs,", nrow(sites), "low-cost sites\n")
