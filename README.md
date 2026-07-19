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

1. **Inspect the mHM NetCDF output** (grid, time range, variable names, long-term
   average pattern) — `scripts/01_explore_netcdf.py`
2. **Build a station metadata table** from CHMI SWE station data (ID, lat/lon,
   observation period, missing values, etc.) — `scripts/02_station_metadata.py`
3. **Extract the modelled SWE time series** at each station's grid cell —
   `scripts/03_extract_modelled_swe.py`
4. **Compute bias and other evaluation indicators** between observed and modelled
   SWE per station — `scripts/04_bias_analysis.py`
5. Results (tables, figures) feed into the actual thesis, which is written separately
   in Overleaf (LaTeX) — this repo is the reproducible analysis backbone, not the
   thesis text itself.

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
data/          Local data folders (raw data itself is gitignored, not committed)
scripts/       Analysis scripts, numbered in pipeline order
results/       Generated tables and figures
docs/          Running notes/log (meeting notes, decisions, email excerpts)
```

## Status

Project just started (July 2026). Supervisor meeting planned for early September or
sooner with intermediate results. See `docs/notes.md` for the up-to-date log.
