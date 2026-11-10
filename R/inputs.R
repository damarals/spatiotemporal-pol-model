# Model inputs read from data/: reference stations, low-cost network, SatPM2.5,
# TROPOMI, ERA5 and site covariates, with the observation ids of the targets.

source_names <- c("ground_no2", "ground_pm25", "moqa")
site_classes <- c("background", "avenue", "bus_terminal")

# Reference PM2.5 months that follow the low-cost rows in the target table. The
# published fits used this order; it fixes the order of the held rows at the
# withheld background location, and so their joint predictive draws.
late_reference_months <- c("2024-08", "2024-11")

#' Files the pipeline reads, in a fixed order; the names give their role. The
#' low-cost network files (moqa/) are not in git.
#'
#' @param dir data directory.
#' @return named character vector of paths.
data_files <- function(dir = "data") {
  files <- c(sites = "processed/sites.csv", moqa_sites = "moqa/sites.csv",
    classes = "decisions/seuma-site-classes.csv", no2 = "processed/seuma-no2-monthly.csv",
    pm25 = "processed/seuma-pm25-monthly.csv", overlap = "decisions/calibration-overlap.csv",
    moqa = "moqa/moqa-monthly.csv", pairs = "moqa/calibration-pairs.csv",
    acag = "processed/acag-monthly.csv", acag_cells = "processed/acag-cells.csv",
    column = "processed/tropomi-monthly.csv", column_cells = "processed/tropomi-cells.csv",
    era5 = "processed/era5-monthly.csv", site_covariates = "processed/site-covariates.csv",
    moqa_covariates = "moqa/site-covariates.csv", grid_covariates = "processed/grid-covariates.csv",
    deployments = "decisions/seuma-deployments.csv", boundary = "processed/boundary.gpkg",
    basemap = "processed/basemap.gpkg")
  stats::setNames(file.path(dir, files), names(files))
}

#' Read every table the model needs.
#'
#' @param dir data directory with processed/, decisions/ and moqa/.
#' @param files data_files(dir), possibly without names (as a targets file target).
#' @return list: targets (one row per monthly reference or low-cost mean, in
#'   target order), acag (municipal SatPM2.5), acag_cells, column (TROPOMI cell
#'   months), column_cells (corners), calibration (daily low-cost/BAM pairs),
#'   era5, site_covariates and sites (low-cost sites first, then the reference
#'   locations).
read_inputs <- function(dir = "data", files = data_files(dir)) {
  files <- .data_paths(files)
  read <- function(name) .read_csv(files[[name]])
  moqa_files <- files[c("moqa_sites", "moqa", "pairs", "moqa_covariates")]
  if (!all(file.exists(moqa_files)))
    stop("The low-cost network files are missing from ", dirname(moqa_files[[1L]]),
      " (available from TRAMA on request).", call. = FALSE)
  sites <- .stack(read("moqa_sites"), read("sites"))
  classes <- read("classes")
  no2 <- read("no2")
  pm25 <- read("pm25")
  overlap <- read("overlap")
  moqa <- read("moqa")
  pairs <- read("pairs")

  overlap <- overlap[overlap$source == "seuma_pm25", , drop = FALSE]
  pm25 <- pm25[!paste(pm25$location_id, pm25$month) %in% paste(overlap$location_id, overlap$month), ,
    drop = FALSE]
  late <- pm25$location_id == "fortaleza-seuma-hq-cajazeiras" & pm25$month %in% late_reference_months

  at <- function(location) {
    i <- match(location, sites$location_id)
    if (anyNA(i)) stop("Locations without coordinates: ", paste(unique(location[is.na(i)]),
      collapse = ", "), call. = FALSE)
    i
  }
  class_of <- function(location) {
    class <- classes$class[match(location, classes$location_id)]
    if (anyNA(class) || !all(class %in% site_classes))
      stop("Reference locations need a site class in seuma-site-classes.csv.", call. = FALSE)
    class
  }
  rows <- function(id, pollutant, location, month, value, source, device, class)
    data.frame(observation_id = id, pollutant = pollutant, physical_location_id = location,
      month = month, observed = value, source = source, device_id = device, site_class = class,
      x_m = sites$x_m[at(location)], y_m = sites$y_m[at(location)], stringsAsFactors = FALSE)
  reference <- function(x) rows(paste0("seuma:", x$location_id, ":", x$month), "PM25", x$location_id,
    x$month, x$value, "ground_pm25", "", class_of(x$location_id))
  device_sites <- sites$location_id[match(moqa$device_id, sites$device_id)]
  if (!identical(device_sites, moqa$location_id))
    stop("Each low-cost device must sit at the location sites.csv gives it.", call. = FALSE)
  targets <- rbind(
    rows(paste("seuma", no2$station_id, no2$month, sep = ":"), "NO2", no2$location_id, no2$month,
      no2$value, "ground_no2", "", class_of(no2$location_id)),
    reference(pm25[!late, , drop = FALSE]),
    rows(paste0("moqa:", moqa$device_id, ":", moqa$month), "PM25", moqa$location_id, moqa$month,
      moqa$value, "moqa", moqa$device_id, "sensor"),
    reference(pm25[late, , drop = FALSE]))
  rownames(targets) <- NULL
  if (anyDuplicated(targets$observation_id)) stop("Duplicate target ids.", call. = FALSE)
  if (any(!(targets$observed > 0))) stop("Concentrations must be positive.", call. = FALSE)

  acag <- read("acag")
  column <- read("column")
  list(targets = targets,
    acag = data.frame(observation_id = paste0("acag:", acag$month), month = acag$month,
      value = acag$value, stringsAsFactors = FALSE),
    acag_cells = read("acag_cells"),
    column = data.frame(observation_id = paste0("tropomi:", column$month, ":", column$cell_id),
      cell_id = column$cell_id, month = column$month, value = column$value, sd = column$sd,
      stringsAsFactors = FALSE),
    column_cells = read("column_cells"),
    calibration = data.frame(device_id = pairs$device_id, location_id = pairs$location_id,
      day = pairs$day, month = substr(pairs$day, 1L, 7L), moqa = pairs$moqa, bam = pairs$bam,
      stringsAsFactors = FALSE),
    era5 = read("era5"),
    site_covariates = .stack(read("moqa_covariates"), read("site_covariates")),
    sites = sites)
}

#' One path of data_files() by its role.
#'
#' @param files data_files(), possibly without names.
#' @param role a name of data_files().
data_file <- function(files, role) .data_paths(files)[[role]]

# Names data_files() gives, restored by position when a targets file target dropped them.
.data_paths <- function(files) {
  roles <- names(data_files())
  if (length(files) != length(roles)) stop("Expected the ", length(roles), " paths of data_files().",
    call. = FALSE)
  if (is.null(names(files))) names(files) <- roles
  if (!identical(names(files), roles)) stop("The paths must follow data_files().", call. = FALSE)
  files
}

# A CSV table with its identifier and month columns read as text.
.read_csv <- function(path) {
  text <- intersect(c("month", "device_id", "cell_id", "day"), names(utils::read.csv(path, nrows = 1L)))
  utils::read.csv(path, stringsAsFactors = FALSE,
    colClasses = stats::setNames(rep("character", length(text)), text))
}

# Rows of the low-cost part followed by the public part, as one table.
.stack <- function(moqa, public) {
  out <- rbind(moqa, public)
  rownames(out) <- NULL
  out
}
