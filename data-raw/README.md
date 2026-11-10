# Building `data/` from the original sources

The processed data are already in `data/`, so the pipeline does not need these scripts. They
document and repeat how `data/` was built. Run them by hand, once, from the repository root and in
order (`03-moqa.R` needs the reference workbooks of `02-seuma.R`; `07-covariates.R` needs the
output of `05-acag.R`):

```sh
Rscript data-raw/01-ibge.R
Rscript data-raw/02-seuma.R --download
...
Rscript data-raw/07-covariates.R
```

Raw downloads go to `work/raw/<source>/` (not in git). The functions are in `data-raw/R/`; the URL,
size and SHA-256 of each download are in `data-raw/sources/`. Credentials are read from
environment variables and never stored in the repository. Scripts 01 and 04-07 first write to
`work/tmp/<source>/`, compare with the file in `data/` and report whether they are identical; they
replace a committed table only by an identical one, so a new download never changes `data/`
silently. Scripts 02 and 03 write `data/` directly.
Times were measured on an 8-core machine with the raw files already in place.

| Script | Reads or downloads | Credentials | Writes | Time |
|---|---|---|---|---|
| `01-ibge.R` | Municipal boundary of Fortaleza (IBGE 2022), municipalities of Ceará (2022) and urbanized areas (2015), downloaded with `geobr` | none | nothing: compares the download with `data/processed/boundary.gpkg` (bytes and geometry) and `data/processed/basemap.gpkg` (attributes and geometry of each layer) and keeps the committed files | download time |
| `02-seuma.R` | 62 public monthly workbooks of the mobile reference station (Prefeitura de Fortaleza site, `sources/seuma.csv`) into `work/raw/seuma/`; `--download` fetches them, otherwise they are only verified | none | `data/processed/seuma-no2-monthly.csv`, `data/processed/seuma-pm25-monthly.csv` | 25 s |
| `03-moqa.R` | Low-cost network export `moqa_data.csv` (1.2 GB) and operational workbook `moqa_sensor_metadata.xlsx` in `work/raw/moqa/`, available from TRAMA on request (sizes and SHA-256 in `data/moqa/README.md`); `work/raw/seuma/` | none | `data/moqa/moqa-monthly.csv`, `data/moqa/calibration-pairs.csv`, `data/moqa/sites.csv`, `data/processed/sites.csv` | 1.5 min, 1.9 GB of memory |
| `04-tropomi.R` | The 72 monthly NetCDF grids in `data/processed/tropomi/` and `data/processed/boundary.gpkg` | none | `data/processed/tropomi-monthly.csv`, `data/processed/tropomi-cells.csv` | 2 s |
| `05-acag.R` | 60 monthly SatPM2.5 V6.GL.03 NetCDF files (3.1 GB) in `work/raw/acag/` (or `ACAG_RAW_DIR`), downloaded from the public bucket in `sources/acag.csv` when missing | none | `data/processed/acag-monthly.csv`, `data/processed/acag-cells.csv` | 8 s, plus the download |
| `06-era5.R` | 8 hourly ERA5 files from the Open-Meteo archive API in `data-raw/era5/` (in git; URLs and SHA-256 in `data-raw/era5/provenance.tsv`) | none | `data/processed/era5-monthly.csv` | 4 s |
| `07-covariates.R` | 19 OpenStreetMap extracts in `work/raw/osm/` (or `OSM_RAW_DIR`), the WorldPop 2020 raster in `work/raw/worldpop/` (or `WORLDPOP_RAW_PATH`; URL in `sources/worldpop.csv`), `data/processed/acag-cells.csv` with the SatPM2.5 files in `work/raw/acag/` (or `ACAG_RAW_DIR`), and the site tables (`data/moqa/sites.csv` when present) | none | `data/processed/site-covariates.csv`, `data/moqa/site-covariates.csv` (when `data/moqa/sites.csv` is present), `data/processed/grid-covariates.csv` | 40 s |

Notes:

- **TROPOMI.** `04-tropomi.R` starts from the 72 monthly grids in `data/processed/tropomi/`.
  The step from the 2,182 Sentinel-5P L2 NO2 orbits of 2020-2025 (1.01 TiB, listed with their
  SHA-256 in `sources/tropomi.csv`) to those grids is not automated: `data-raw/R/tropomi.R`
  queries the catalogue, downloads and verifies the orbits and writes the HARP operations, but
  does not call HARP. The grids were made with HARP 1.29 in two steps. For each orbit, the
  scanlines within 1 degree of the grid (plus one of margin) of the variables time, delta_time,
  latitude, longitude, latitude_bounds, longitude_bounds, qa_value,
  nitrogendioxide_tropospheric_column and its precision were read with `gdalmdimtranslate`
  (HTTP range requests) into a small HARP-convention NetCDF with validity = round(100 qa_value)
  (`.locate_tropomi_window()` gives the scanlines) and converted with
  `harpconvert --no-history --operations '<per-orbit operations>'`, where the per-orbit
  operations are
  `set("propagate_uncertainty","correlated");tropospheric_NO2_column_number_density_validity>75;keep(datetime_start,latitude,longitude,latitude_bounds,longitude_bounds,tropospheric_NO2_column_number_density,tropospheric_NO2_column_number_density_uncertainty);bin_spatial((-3.920,-3.885,...,-3.675),(-38.675,-38.640,...,-38.395))`
  (a 0.035 degree grid aligned at 0 degrees around the municipal polygon,
  `.aligned_harp_grid()`; orbits left empty by the filter are skipped). The orbits of each month
  were then combined with
  `harpmerge --no-history --post-operations 'set("propagate_uncertainty","correlated");bin()'`.
  Both operation strings are written by `.harp_operations()`. Downloading the orbits from the
  Copernicus Data Space Ecosystem needs an account and its S3 keys in `CDSE_S3_ACCESS_KEY` and
  `CDSE_S3_SECRET_KEY`.
- **ERA5.** The hourly files run to 2026-10, because the model that fills the missing
  boundary-layer height (2024-01 to 2024-06) was fitted on every available hour. The model uses
  only 2020-01 to 2025-12.
- **OpenStreetMap.** The queries are in `data-raw/osm/`; `sources/osm.csv` gives the hash and
  database timestamp of each extract used. Overpass serves the current database, so a new
  download gives slightly different covariates. The extracts are of two dates. The covariates of
  the 22 monitoring locations use the major roads (q1, database of 2026-11-07) and the
  site-buffer extract (q2), whose database timestamp is 2026-06-01: the mirror that served it
  was five months behind. The covariates of the 330 map cells use q1, the bus stops (q3,
  2026-11-08) and the 16 road tiles (q5, 2026-11-08). The model estimates the PM2.5 covariate
  coefficient on the number of bus stops within 500 m at the locations (June 2026 database) and
  the PM2.5 map applies it to the counts at the cells (November 2026 database), so a bus stop
  added or removed in between shifts the map at that cell but not the fit. The q2 query file
  keeps only the 8 reference locations; the query over all 22 sites, with the low-cost
  coordinates, is with the low-cost network files (`data/moqa/osm-q2-site-buffers.txt`).
- **Boundary.** `data/processed/boundary.gpkg` is the polygon `geobr` served when the data were
  built (release v1.7.0, 336.8 km2 in EPSG:31984). The current IBGE polygon is clipped at the
  coastline (312.1 km2); `01-ibge.R` reports the difference and keeps the committed file, which
  defines the SatPM2.5 and TROPOMI cells of the model.
