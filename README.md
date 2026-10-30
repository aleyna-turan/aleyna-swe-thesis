# Bachelor Thesis: Evaluation of mHM-based Snow Water Equivalent Simulations for Czechia

## Project summary

This repository contains the data-processing pipeline and analysis code for a bachelor
thesis comparing **mHM-simulated snow water equivalent (SWE)** against **ground-based
SWE observations from CHMI (Czech Hydrometeorological Institute)** stations.

The mHM simulations are gridded, based on E-OBS meteorological forcing, covering
1950–2025 (extended in later updates) at daily time step and 0.125° x 0.125° spatial
resolution, cropped to the Czechia domain. Full pan-European simulations exist upstream
of this cropped file.

## Workflow overview

Scripts are numbered in pipeline order and all configuration (periods,
thresholds, criteria) sits at the top of each file.

| # | Script | What it does |
|---|---|---|
| 01 | `01_explore_netcdf.py` | Inspect the mHM NetCDF output: grid, time range, variables, long-term average pattern |
| 02 | `02_station_metadata.py` | Build the station metadata table from CHMI SWE data |
| 03 | `03_extract_modelled_swe.py` | Extract modelled SWE at each station's grid cell, on observation dates only |
| 04 | `04_bias_analysis.R` | First bias and evaluation indicators (superseded by 07; output in `results/*/archive/`) |
| 05 | — | Does not exist. A first-pass availability check, superseded by 06 before it was saved. The numbering gap is intentional. |
| 06 | `06_station_qc_selection.py` | Screen stations for regular measurement over each candidate period |
| 07 | `07_seasonal_statistics.py` | Seasonal metrics, model evaluation, and Theil–Sen/Mann–Kendall trends per hydrological year |
| 08 | `08_duration_bias_diagnosis.py` | Threshold sensitivity and date-by-date disagreement analysis |

Results feed into the thesis text, which is written separately in Overleaf
(LaTeX). This repo is the reproducible analysis backbone, not the thesis text.

### Conventions

- **Hydrological year**: 1 Nov to 31 Oct, named after the year it ends in
  (Czech / CHMI convention). HY 2025 = 1 Nov 2024 to 31 Oct 2025.
- **Period labels** are calendar spans: `1976-2025` = HY 1977–2025, 49 years.
- **Snow threshold**: SWE ≥ 5 mm, justified by `results/figures/04_threshold_tests/`.
- **Sampling**: SWE is measured weekly, so timing metrics have ~7-day resolution.
  CHMI moved the measurement day from Friday to Monday during 1970; scripts
  homogenise this per year. See `docs/notes.md`.

## Data sources

- **mHM simulations**: provided by supervisor, based on E-OBS forcing, NetCDF format.
  Not included in this repo (too large, not ours to publish) — see `data/README.md`.
- **CHMI station SWE data**: obtained from the CHMI portal (ČHMÚ). Exact access route
  (public portal vs. request) still being confirmed — noted in `docs/notes.md`.
- **Future extension**: ISIMIP3b-forced mHM SWE simulations (1850–2100, historical +
  3 SSP scenarios x 5 GCMs) — not yet in scope, noted for later.

## Repository structure

```
environment/   Setup notes and Python package requirements
data/          Local data folders (raw data itself is gitignored)
scripts/       Analysis scripts, numbered in pipeline order
results/
  figures/     FIGURES.md documents every figure: what, how, why
  tables/      TABLES.md documents every table and which script writes it
docs/          Running notes/log, newest entry last
```

Both `results/figures/` and `results/tables/` are grouped by the **question
each output answers**, not by study period:

| Folder | Question | Script |
|---|---|---|
| `00_inputs/` | the data everything is built from (tables only) | 02, 03 |
| `01_data_screening/` | which stations measured regularly | 06 |
| `02_model_evaluation/` | how well mHM reproduces the observations | 07 |
| `03_trends/` | is the snow changing, and does the model capture it | 07 |
| `04_threshold_tests/` | is the duration result real or an artefact | 08 |
| `archive/` | superseded output, kept so older notes stay traceable | — |

Start with **`results/figures/FIGURES.md`** — it explains every figure,
including how each plotted quantity is computed and what to watch out for.

## Status

October 2026. The processing pipeline is complete and reproduces identically on
two machines. Current state:

- **Study period**: 1976–2025 recommended (235 stations, 49 hydrological years);
  awaiting confirmation with the supervisor.
- **Main result**: mHM's peak SWE bias becomes more negative with elevation,
  −2.29 mm per 100 m (Kendall tau = −0.182, p = 3.3e-05, 235 stations), and the
  relationship holds when the three stations above 1000 m are removed.
- **Trends**: peak SWE −3.25 mm/decade and snow season length −7.78 days/decade
  over 1976–2025, both robust to the snow threshold; mHM reproduces them.
- **Known limitation**: only 7 selected stations lie above 800 m and 3 above
  1000 m, in every candidate period. This is a network limitation.

See `docs/notes.md` for the full log, including two corrections to earlier
findings.
