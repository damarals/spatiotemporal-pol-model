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

## Sources and terms

The MIT License of the repository covers the code only. The files in this folder, the ERA5 files
in `data-raw/era5/` and the results derived from them in `results/` keep the terms of their
sources:

- **SEUMA** (Secretaria Municipal do Urbanismo e Meio Ambiente, Prefeitura de Fortaleza): hourly
  NO2 and PM2.5 of the mobile reference station, from the public monthly workbooks
  (`processed/seuma-*.csv`).
- **MoQA / TRAMA**: low-cost network data of the Prefeitura de Fortaleza, managed by TRAMA,
  Universidade Federal do Ceará; available on request, not redistributed here (`moqa/`).
- **TROPOMI**: Sentinel-5P TROPOMI L2 NO2 (https://doi.org/10.5270/S5P-9bnp8q8), processed with
  HARP (`processed/tropomi*`). Contains modified Copernicus Sentinel data (2020-2025).
- **ERA5** (https://doi.org/10.24381/cds.adbb2d47): obtained through the
  [Open-Meteo](https://open-meteo.com/) historical weather API and aggregated to monthly values
  (`processed/era5-monthly.csv`, `data-raw/era5/`). Contains modified Copernicus Climate Change
  Service information (2026). Open-Meteo data are under
  [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/).
- **Copernicus notice**: neither the European Commission nor ECMWF is responsible for any use that
  may be made of the Copernicus information or data it contains.
- **SatPM2.5** (ACAG, Atmospheric Composition Analysis Group, Washington University in St. Louis):
  monthly surface PM2.5 estimates, version V6.GL.03 (`processed/acag-*.csv`). Cite Shen, S., Li,
  C., van Donkelaar, A., Jacobs, N., Wang, C., Martin, R. V. (2024). Enhancing global estimation
  of fine particulate matter concentrations by including geophysical a priori information in deep
  learning. ACS ES&T Air 1(5), 332-345. https://doi.org/10.1021/acsestair.3c00054.
- **OpenStreetMap**: © OpenStreetMap contributors. The site and grid covariates in
  `processed/site-covariates.csv` and `processed/grid-covariates.csv` are derived from
  OpenStreetMap extracts and are available under the
  [Open Database License (ODbL) 1.0](https://opendatacommons.org/licenses/odbl/1-0/), not under
  the MIT License of the code.
- **WorldPop** (www.worldpop.org): Brazil 2020 constrained population counts, UN adjusted, 100 m,
  under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/). The covariate `pop_1km` is
  derived from it (population summed within 1 km). Cite Bondarenko, M., Kerr, D., Sorichetta, A.,
  Tatem, A. J. (2020). Census/projection-disaggregated gridded population datasets, adjusted to
  match the corresponding UNPD 2020 estimates, for 183 countries in 2020 using Built-Settlement
  Growth Model (BSGM) outputs. WorldPop, University of Southampton.
  https://doi.org/10.5258/SOTON/WP00685.
- **IBGE** (Instituto Brasileiro de Geografia e Estatística): municipal boundaries (2022) and
  urbanized areas (2015), obtained through the `geobr` package (`processed/boundary.gpkg`,
  `processed/basemap.gpkg`). Source: IBGE.
