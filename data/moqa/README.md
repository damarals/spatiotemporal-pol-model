# Low-cost network data (not in this repository)

The low-cost PM2.5 network data belong to the Fortaleza municipal government and are managed by
TRAMA, Universidade Federal do Ceará (UFC). They are not redistributed here; they are available from
TRAMA on request. Put the files below in this folder to run the analysis; the tests that need
them skip when the folder holds only this file.

## Files

`moqa-monthly.csv`: 21 field device-months.

| Column | Meaning |
|---|---|
| device_id | device identifier of the network |
| location_id | installation site (a row of `sites.csv`) |
| month | YYYY-MM |
| value | monthly mean PM2.5 on the uncalibrated sensor scale, ug/m3 |

`calibration-pairs.csv`: 66 daily pairs of the device colocated with the reference station at its
background location (September and October 2024, January 2025).

| Column | Meaning |
|---|---|
| device_id | colocated device |
| location_id | reference location |
| day | YYYY-MM-DD |
| common_hours | hours valid for both instruments that day (at least 16) |
| moqa | mean of the sensor over the common hours, ug/m3 (sensor scale) |
| bam | mean of the reference PM2.5 (beta attenuation monitor) over the same hours, ug/m3 |

`sites.csv`: the 14 low-cost sites, with the columns of `data/processed/sites.csv` (location_id,
source `moqa`, device_id, x_m, y_m in EPSG:31984). The model reads these rows before the 8
reference locations; the order fixes its site index.

`site-covariates.csv`: the covariates of the 14 low-cost sites, with the columns of
`data/processed/site-covariates.csv` (OpenStreetMap-derived, ODbL 1.0).

`moqa-periods.csv` (21 operating periods: site, cadence, software regime, clock shift),
`moqa-month-decisions.csv` (5 device-months excluded by hand) and `moqa-qc.csv` (the automatic
quality-control threshold): the decisions `data-raw/03-moqa.R` applies.

`osm-q2-site-buffers.txt`: the Overpass query of the OpenStreetMap site extract over all 22
sites (`data-raw/osm/q2-site-buffers.txt` keeps the 8 reference locations).

## Rebuilding them

`Rscript data-raw/03-moqa.R` rebuilds `moqa-monthly.csv`, `calibration-pairs.csv` and `sites.csv` (and
`data/processed/sites.csv`) from the network export and its operational workbook, placed in
`work/raw/moqa/`:

| File | Bytes | SHA-256 |
|---|---|---|
| `moqa_data.csv` | 1,163,839,584 | `51659871eb71c090b885bc9be7b52d3a7e00b8b6b08f1a97017cc01ec842fa6e` |
| `moqa_sensor_metadata.xlsx` | 166,384 | `8330577a81000c64b43f26acc27b1d08487ec53c166bea27cbf5a702f07eb87d` |

It also needs the reference station workbooks in `work/raw/seuma/` (`Rscript data-raw/02-seuma.R --download`).
`Rscript data-raw/07-covariates.R` rebuilds `site-covariates.csv`.
