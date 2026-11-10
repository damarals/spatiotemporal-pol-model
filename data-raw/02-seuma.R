# Reference station (SEUMA): monthly NO2 and PM2.5 means from the public workbooks.
# Usage: Rscript data-raw/02-seuma.R [--download] [output root, default "."]
# Reads work/raw/seuma/ (filled by --download) and writes
# data/processed/seuma-no2-monthly.csv and data/processed/seuma-pm25-monthly.csv.
for (f in c("write.R", "seuma-read.R", "seuma-qc.R")) source(file.path("data-raw/R", f))
suppressPackageStartupMessages(library(data.table))
args <- commandArgs(trailingOnly = TRUE)
out <- c(setdiff(args, "--download"), ".")[1]

workbooks <- seuma_workbooks()
if ("--download" %in% args) download_seuma(workbooks) else verify_seuma(workbooks)
deployments <- read_table("data/decisions/seuma-deployments.csv")
samples <- read_seuma_samples(workbooks)

no2 <- seuma_no2_monthly(samples, deployments)
hours <- seuma_pm25_hours(samples, deployments, read_table("data/decisions/seuma-pm25-months.csv"),
                          read_table("data/decisions/seuma-pm25-flagged-hours.csv"))
pm25 <- seuma_pm25_monthly(hours)

dir.create(file.path(out, "data/processed"), recursive = TRUE, showWarnings = FALSE)
write_table(no2, file.path(out, "data/processed/seuma-no2-monthly.csv"))
write_table(pm25, file.path(out, "data/processed/seuma-pm25-monthly.csv"))
cat(nrow(no2), "NO2 and", nrow(pm25), "PM2.5 location-months\n")
