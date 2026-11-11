# More data sources do not guarantee better urban air pollution predictions

Code, processed data and results of the paper by Daniel de Amaral da Silva, Yves Rybarczyk, Bruno
Vieira Bertoncini and Danielo G. Gomes, submitted to *Environmental Modelling & Software*.

The paper links each data source to latent concentration fields through its own observation model
and tests each model component by backward elimination on held-out predictions, for monthly
NO<sub>2</sub> and PM<sub>2.5</sub> in Fortaleza, Brazil, from 2020 to 2025. At new traffic
locations, the continuous ranked probability score falls from 2.08 to 1.65 µg/m³ for
NO<sub>2</sub> and from 1.28 to 0.98 µg/m³ for PM<sub>2.5</sub>. The main limit is the monitoring
design: one background location and seven traffic locations observed for one to three months each.

The repository holds the Stan program ([`stan/fusion.stan`](stan/fusion.stan)), the `targets`
pipeline ([`_targets.R`](_targets.R)), the functions it calls ([`R/`](R)), the processed data
([`data/`](data)), the scripts that build them ([`data-raw/`](data-raw)) and the results
([`results/`](results)). Each model is a set of data flags in [`R/models.R`](R/models.R), and the
validation units are in [`R/units.R`](R/units.R). It is the analysis code for one city, not a
general-purpose package.

**The low-cost network data are not redistributed, and the pipeline does not run without them**
(see [Data](#data)).

## What you can run

| Goal | Needs | Command | Time |
|---|---|---|---|
| Read the result tables and figures | nothing | open [`results/`](results) | none |
| Run the tests | R packages | `Rscript tests/run.R` | about 1 min |
| Rebuild the public tables | R packages, downloads for some | scripts in [`data-raw/`](data-raw/README.md), such as `Rscript data-raw/06-era5.R` | seconds to minutes, plus the downloads |
| Refit every model and rebuild `results/` | R packages, CmdStan, `data/moqa/` | `Rscript -e 'targets::tar_make()'` | about 20 h on 8 cores, the sum of steps timed separately |

The tests should end with `FAIL 0`; those that need `data/moqa/` skip without it.

Table 4 is `results/scores.csv`, except its Full column, which is in `results/full-vs-final.csv`.
Tables 5 and 6 are rows of `results/elimination.csv`, which holds all 186 comparisons. Figures 1 to
4 are in `results/figures/`.

Seeds are fixed, but Stan does not promise identical draws across compilers or hardware, so expect
small Monte Carlo differences on another machine.

## Install

Tested on Ubuntu 22.04 under WSL 2 with an 8-core AMD Ryzen 7 5700X, 15 GB of memory, g++ 11.4,
R 4.6.1 and CmdStan 2.37.0. On Windows use WSL, because the fits run their chains with
`parallel::mclapply()`. Each fit uses four CPU cores.

```sh
# System libraries (Ubuntu)
sudo apt install g++ make cmake gdal-bin libgdal-dev libgeos-dev libproj-dev libnetcdf-dev \
  libudunits2-dev libcurl4-openssl-dev libssl-dev libxml2-dev libfontconfig1-dev libpng-dev \
  libglpk-dev libabsl-dev

# R packages at the versions in renv.lock (needs R 4.6.1), about 10 min from source
MAKEFLAGS=-j8 Rscript -e 'renv::restore(prompt = FALSE)'

# CmdStan 2.37.0, needed only to fit the models
curl -L https://github.com/stan-dev/cmdstan/releases/download/v2.37.0/cmdstan-2.37.0.tar.gz \
  | tar -xz -C "$HOME"
make -C "$HOME/cmdstan-2.37.0" -j4 build
export CMDSTAN="$HOME/cmdstan-2.37.0"
```

## Data

[`data/README.md`](data/README.md) describes each processed file and lists the sources and their
terms. [`data-raw/README.md`](data-raw/README.md) documents how the files were built. The MIT
License covers the code only; the data keep the terms of their sources, including the Copernicus
notices, the Open Database License for the OpenStreetMap covariates and CC BY 4.0 for WorldPop and
Open-Meteo.

The low-cost PM<sub>2.5</sub> network (MoQA) data belong to the Fortaleza municipal government and
are managed by TRAMA (Laboratório de Transportes e Meio Ambiente, Universidade Federal do Ceará).
They are available from TRAMA on request. Put the files listed in
[`data/moqa/README.md`](data/moqa/README.md) in `data/moqa/`; without them the pipeline stops at
its first target, `data_paths`.

---

This study was financed in part by the Coordenação de Aperfeiçoamento de Pessoal de Nível Superior
– Brasil (CAPES) – Finance Code 001.
