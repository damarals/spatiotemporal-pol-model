# Data

The model reads only this folder: the files listed by `data_files()` in `R/inputs.R`, through
`read_inputs()` and `read_grid()` (`R/maps.R`). The scripts in `data-raw/` build `processed/` and
`moqa/` from the original sources. Months are `YYYY-MM`; coordinates `x_m`, `y_m` are metres in
SIRGAS 2000 / UTM zone 24S (EPSG:31984). Every double in the CSV files reads back exactly as
written.

## `processed/`

| File | Rows | Columns | Units and notes | Built by |
|---|---|---|---|---|
| `sites.csv` | 8 | location_id, source (`seuma`), device_id (NA), x_m, y_m | the reference locations; the model reads the 14 low-cost sites of `moqa/sites.csv` first, then these rows, and this order fixes its site index | `03-moqa.R` |
| `seuma-no2-monthly.csv` | 49 | location_id, station_id (P00-P08), month, value | reference NO2, µg/m3 | `02-seuma.R` |
| `seuma-pm25-monthly.csv` | 47 | location_id, month, value | reference PM2.5, µg/m3; the 4 months in `decisions/calibration-overlap.csv` are left out of the model | `02-seuma.R` |
| `acag-monthly.csv` | 60 | month, value | SatPM2.5 municipal area-weighted mean, µg/m3, 2020-01 to 2024-12 | `05-acag.R` |
| `acag-cells.csv` | 330 | cell, lon, lat, weight, x_m, y_m | 0.01 degree SatPM2.5 cells intersecting the municipality, centres (WGS 84 and EPSG:31984) and area weights (sum 1); the map grid | `05-acag.R` |
| `tropomi-monthly.csv` | 1584 | cell_id, month, value, sd | TROPOMI tropospheric NO2 column monthly mean and its standard deviation, 1e-5 mol/m2, 22 cells x 72 months | `04-tropomi.R` |
| `tropomi-cells.csv` | 88 | cell_id, corner, x_m, y_m | corners of the 22 TROPOMI grid cells | `04-tropomi.R` |
| `tropomi/` | 72 files | | monthly gridded TROPOMI NO2 (HARP NetCDF), input of `04-tropomi.R` | `data-raw/R/tropomi.R` |
| `era5-monthly.csv` | 72 | month, log_vent, t2m, precip | 2020-01 to 2025-12: log of the mean hourly ventilation (10 m wind speed x boundary-layer height, m2/s), mean 2 m temperature (°C), total precipitation (mm) | `06-era5.R` |
| `site-covariates.csv` | 8 | location_id and the 8 covariates below | raw values at the reference locations (the low-cost sites are in `moqa/site-covariates.csv`); the model transforms them | `07-covariates.R` |
| `grid-covariates.csv` | 330 | cell and the 8 covariates below | cells as in `acag-cells.csv` | `07-covariates.R` |
| `boundary.gpkg` | 1 | layer `fortaleza` | municipal boundary of Fortaleza, IBGE 2022, EPSG:4674; sets the SatPM2.5 and TROPOMI cells and clips the maps | `01-ibge.R` (checks it) |
| `basemap.gpkg` | 184 + 307 | layers `municipalities`, `urban_areas` | municipalities of Ceará (IBGE 2022) and urbanized areas (IBGE 2015), EPSG:4674; Figure 1 only | `01-ibge.R` |

Covariates: `n_bus_stops_250m`, `n_bus_stops_500m` (bus stops within 250 m and 500 m),
`len_mtps_100m`, `len_mtps_500m` (length of motorway, trunk, primary and secondary roads within
100 m and 500 m, m), `len_vehicular_250m` (length of vehicular roads within 250 m, m),
`dist_mtp_road_m` (distance to the nearest motorway, trunk or primary road, m), `pop_1km`
(population within 1 km) and `acag_log_ratio_municipal` (log of the SatPM2.5 cell mean over the
municipal mean).

## `decisions/`

Choices made by hand, read by `data-raw/` or by the model. Each has its reason or evidence in a
column.

| File | Rows | Contents |
|---|---|---|
| `seuma-site-classes.csv` | 8 | class of each reference location: `background`, `avenue` or `bus_terminal` |
| `calibration-overlap.csv` | 4 | reference PM2.5 months at the background location left out of the monthly likelihood (2024-09, 2024-10, 2024-12, 2025-01) |
| `seuma-deployments.csv` | 9 | deployments P00-P08 of the mobile reference station: location, first and last month |
| `seuma-locations.csv` | 8 | name, address and coordinates of each reference location |
| `seuma-workbooks.csv` | 62 | layout and concentration unit of each public monthly workbook |
| `seuma-pm25-months.csv` | 62 | reference PM2.5 months reviewed or withheld |
| `seuma-pm25-flagged-hours.csv` | 958 | PM2.5 hours flagged in the operator's monthly reports |

The decisions about the low-cost network are with its data in `moqa/`.

## `moqa/`

Everything derived from the low-cost network: the monthly means, the calibration pairs, the 14
low-cost sites and their covariates, the operating periods and the quality-control decisions.
They are not in git; see `moqa/README.md`.
