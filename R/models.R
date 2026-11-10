# Every model the paper reports: the full model, the backward elimination
# (round 1 on the full model without TROPOMI, round 2 on the reduced model,
# then the final model without TROPOMI), the final model and its sensitivity
# analyses, with the units each was fitted on.
#
# Columns (1 = component in the model):
#   tropomi          TROPOMI column
#   temporal_loading TROPOMI loading lambda_T on the NO2 temporal component
#   sign_constraint  lambda_S > 0 (TROPOMI spatial loading)
#   no2_profile, pm25_profile   spatial profiles
#   profile_cor      correlation rho_S between the profiles' coordinates
#   r                innovation correlation between the pollutants
#   amplitude        kappa: traffic NO2 follows the temporal component
#   classes          traffic classes with their own offset (2) or one (1)
#   deviation        PM2.5 deviation at each reference location
#   low_cost         low-cost network and its calibration
#   era5, satpm25    ERA5 drivers; SatPM2.5 municipal mean
#   seasonal         seasonal factor of SatPM2.5
#   free_gain        SatPM2.5 gain b_A estimated (0: fixed at one)
#   no2_cov, pm25_cov  site covariate selected in each fit (0: none)
#   knots, width_m   profile basis
#   priors           scale priors: quartic tails or half-normal
#   units            unit groups fitted (see model_units())
models <- utils::read.table(header = TRUE, stringsAsFactors = FALSE, text = "
model                                     tropomi temporal_loading sign_constraint no2_profile pm25_profile profile_cor r amplitude classes deviation low_cost era5 satpm25 seasonal free_gain no2_cov pm25_cov knots width_m priors      units
full                                      1       1                1               1           1            1           1 1         2       1         1        1    1       1        0         1       1        4     8000    quartic     omitted;traffic;background;year-2023;year-2024
full-unsigned                             1       1                0               1           1            1           1 1         2       1         1        1    1       1        0         1       1        4     8000    quartic     traffic
full-no-tropomi                           0       0                0               1           1            1           1 1         2       1         1        1    1       1        0         1       1        4     8000    quartic     omitted;traffic;background;years;devices
full-no-tropomi-without-r                 0       0                0               1           1            1           0 1         2       1         1        1    1       1        0         1       1        4     8000    quartic     omitted;traffic;years;devices
full-no-tropomi-without-no2-profile       0       0                0               0           1            1           1 1         2       1         1        1    1       1        0         1       1        4     8000    quartic     omitted;traffic;background;years
full-no-tropomi-without-low-cost          0       0                0               1           1            1           1 1         2       1         0        1    1       1        0         1       1        4     8000    quartic     omitted;traffic;year-2025
full-no-tropomi-without-pm25-profile      0       0                0               1           0            1           1 1         2       1         1        1    1       1        0         1       1        4     8000    quartic     omitted;traffic;background;years;devices
full-no-tropomi-without-deviation         0       0                0               1           1            1           1 1         2       0         1        1    1       1        0         1       1        4     8000    quartic     omitted;traffic;background;years;devices
full-no-tropomi-without-era5              0       0                0               1           1            1           1 1         2       1         1        0    1       1        0         1       1        4     8000    quartic     omitted;traffic;background;years
full-no-tropomi-free-gain                 0       0                0               1           1            1           1 1         2       1         1        1    1       1        1         1       1        4     8000    quartic     omitted;traffic;background;years
full-no-tropomi-without-seasonal          0       0                0               1           1            1           1 1         2       1         1        1    1       0        0         1       1        4     8000    quartic     omitted;traffic;background;years
full-no-tropomi-without-satpm25           0       0                0               1           1            1           1 1         2       1         1        1    0       1        0         1       1        4     8000    quartic     omitted;traffic;background;years
full-no-tropomi-without-amplitude         0       0                0               1           1            1           1 0         2       1         1        1    1       1        0         1       1        4     8000    quartic     omitted;traffic;background;years
full-no-tropomi-one-class                 0       0                0               1           1            1           1 1         1       1         1        1    1       1        0         1       1        4     8000    quartic     omitted;traffic;background;years
full-no-tropomi-without-profile-cor       0       0                0               1           1            0           1 1         2       1         1        1    1       1        0         1       1        4     8000    quartic     traffic;background
full-no-tropomi-without-no2-cov           0       0                0               1           1            1           1 1         2       1         1        1    1       1        0         0       1        4     8000    quartic     omitted;traffic;background;years
full-no-tropomi-without-pm25-cov          0       0                0               1           1            1           1 1         2       1         1        1    1       1        0         1       0        4     8000    quartic     traffic;background;devices
reduced                                   0       0                0               1           1            0           1 0         2       1         1        1    1       1        0         1       1        4     8000    quartic     omitted;traffic;background;years;devices
reduced-free-gain                         0       0                0               1           1            0           1 0         2       1         1        1    1       1        1         1       1        4     8000    quartic     omitted;traffic;background;years
reduced-without-no2-cov                   0       0                0               1           1            0           1 0         2       1         1        1    1       1        0         0       1        4     8000    quartic     omitted;traffic;background;years
reduced-without-pm25-profile              0       0                0               1           0            0           1 0         2       1         1        1    1       1        0         1       1        4     8000    quartic     omitted;traffic;background;years;devices
reduced-without-low-cost                  0       0                0               1           1            0           1 0         2       1         0        1    1       1        0         1       1        4     8000    quartic     omitted;traffic;year-2025
reduced-without-deviation                 0       0                0               1           1            0           1 0         2       0         1        1    1       1        0         1       1        4     8000    quartic     omitted;traffic;background;years;devices
reduced-one-class                         0       0                0               1           1            0           1 0         1       1         1        1    1       1        0         1       1        4     8000    quartic     omitted;traffic;background;years
final-no-tropomi                          0       0                0               1           0            0           1 0         2       0         1        1    1       1        1         0       1        4     8000    quartic     omitted;traffic;background;years;devices;all
final-no-tropomi-without-low-cost         0       0                0               1           0            0           1 0         2       0         0        1    1       1        1         0       1        4     8000    quartic     omitted;traffic;year-2025
final                                     1       0                1               1           0            0           1 0         2       0         1        1    1       1        1         0       1        4     8000    quartic     omitted;traffic;background;years;devices;all
final-one-class                           1       0                1               1           0            0           1 0         1       0         1        1    1       1        1         0       1        4     8000    quartic     omitted;traffic
final-eight-knots                         1       0                1               1           0            0           1 0         2       0         1        1    1       1        1         0       1        8     4000    quartic     omitted;traffic;background
final-half-normal                         1       0                1               1           0            0           1 0         2       0         1        1    1       1        1         0       1        4     8000    half_normal omitted;traffic;background
")

#' Model row by name.
model_spec <- function(model) {
  m <- models[models$model == model, , drop = FALSE]
  if (nrow(m) != 1L) stop("Unknown model: ", model, call. = FALSE)
  as.list(m)
}

#' Units a model is fitted on, in validation_units order.
model_units <- function(model) {
  groups <- strsplit(model_spec(model)$units, ";", fixed = TRUE)[[1L]]
  members <- list(omitted = "omitted-months",
    traffic = c("new-traffic-1", "new-traffic-2", "new-traffic-3"),
    background = "withheld-background", years = c("year-2023", "year-2024", "year-2025"),
    "year-2023" = "year-2023", "year-2024" = "year-2024", "year-2025" = "year-2025",
    devices = c("new-devices-1", "new-devices-2", "new-devices-3"), all = "all-data")
  unknown <- setdiff(groups, names(members))
  if (length(unknown)) stop("Unknown unit groups: ", paste(unknown, collapse = ", "), call. = FALSE)
  intersect(validation_units$unit, unlist(members[groups]))
}
