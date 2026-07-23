# Project notes / running log

Informal log of decisions, meeting notes, and reasoning — kept here so it's easy to
pull from when writing the methodology section in Overleaf later.

---

## 2026-07-19 — Project kickoff

**From supervisor's email:**
- Will provide mHM-based gridded SWE simulations, E-OBS forcing, 1950–2025, daily,
  0.125° x 0.125° resolution.
- Currently cropping full European simulations down to Czechia only.
- Sample file location: SharePoint (`czuvpraze-my.sharepoint.com`), folders named by
  the date data was shared.
- **Next step (this repo's current focus):**
  1. Explore CHMI portal for station-based SWE time series.
  2. Build a station metadata table: ID, lat, lon, start/end of observation period,
     number of missing values, + other relevant info.
  3. Extract the corresponding modelled SWE time series from mHM for each station.
  4. Calculate bias and other useful indicators.
- Next meeting: early September, or sooner if there are intermediate results/questions.
- Keep supervisor posted via email or Teams.
- After first results: start official assignment paperwork in UIS.
- **Possible future extension:** mHM SWE simulations forced with ISIMIP3b climate
  projections (1850–2100, historical 1950–2014 + 3 SSP scenarios x 5 GCMs from 2015).
  Supervisor still needs to confirm whether SWE was actually written out in those runs.

**File received:**
- `snowpack_czechia_19500101_20260531.nc` — 83MB.
- Supervisor's note: `cdo timavg` long-term average clearly shows terrain pattern —
  highlands at Czech borders (green/blue), lowlands (red), Alps visible as dark blue.

**Open question — CHMI station data access:**
- CHMI publishes SWE ("vodní hodnota sněhu") station data and interpolated maps on
  their portal (chmi.cz), measured weekly (Mondays, 7:00 CET) at stations with snow
  depth > 3 cm. Element code in their CLIDATA database: `SVH`.
- Not yet clear whether historical per-station daily/weekly time series are
  downloadable directly from the public portal, or whether we need to request them
  from CHMI, or whether the supervisor already has access to hand off.
  **Action: confirm data access route.**

---

## 2026-07-19 — Found CHMI open data source

- Confirmed the mHM file supervisor sent (`snowpack_czechia_19500101_20260531.nc`,
  in `data/raw/`) is the model side; supervisor did not separately provide CHMI
  station data — this needs to be obtained from CHMI directly.
- Found CHMI's free, open data portal: `https://opendata.chmi.cz`. No login needed,
  Creative Commons BY 4.0 license.
- Daily station SWE data lives at:
  `https://opendata.chmi.cz/meteorology/climate/historical_csv/data/daily/snow/`
  - Filename pattern: `dly-{WSI}-SVH.csv`, one file per station, full period of
    record in a single file (one row per day since ~1863 for some stations).
  - Element `SVH` = vodní hodnota sněhu = snow water equivalent (mm).
  - Columns: `STATION, ELEMENT, TIMEFUNC, DT, VALUE, FLAG, QUALITY`
  - `QUALITY` code meaning (from `meta4.csv`): 0=Good, 1=Suspect, 2=Poor,
    3=Estimated, 4=**Missing**, 5=Unknown.
  - SWE is only measured Mondays (and only when snow depth > 3cm), so most days
    in a year legitimately have `QUALITY=4` (missing) — this is expected, not
    a data problem.
  - Found **1147 stations** with SVH files.
- Station metadata (name, coordinates, elevation) lives at:
  `https://opendata.chmi.cz/meteorology/climate/historical_csv/metadata/meta1.csv`
  - Columns: `WSI, GH_ID, BEGIN_DATE, END_DATE, FULL_NAME, GEOGR1 (=lon),
    GEOGR2 (=lat), ELEVATION`
  - **Important**: a single station can have multiple rows if its exact
    coordinates/metadata changed over time (e.g. station 11406 "Cheb" has 4+
    historical entries). Script picks the row overlapping the mHM period
    (1950–2026).
- Wrote `scripts/02_station_metadata.py` to download all 1147 SVH files and
  build `results/tables/station_metadata.csv` (ID, lat, lon, elevation,
  observation start/end based on actual non-missing values, missing count).

---

## 2026-07-20 — Bias analysis complete (Step 5), core assignment done

- Decided to do the statistics/analysis step in **R** instead of Python (more
  comfortable working in R). Data hand-off between languages is trivial since
  everything is plain CSV (`modelled_vs_observed_swe.csv.gz`) - no need to
  redo the NetCDF/extraction work in R.
- R environment: R 4.4.1 + RStudio, packages `dplyr`, `readr`, `ggplot2`.
- Wrote `scripts/04_bias_analysis.R`:
  - For each station: bias (mean of modelled - observed), RMSE, correlation,
    n_obs, joined with metadata (lat/lon/elevation).
  - Output: `results/tables/station_bias_metrics.csv` (1131 of 1147 stations
    got valid metrics; 16 dropped/undefined correlation due to zero-variance
    edge cases - e.g. stations with almost no real snow events recorded).
  - Also made a diagnostic plot: bias vs. elevation
    (`results/figures/bias_vs_elevation.png`).

**Key results:**
- Overall: mean bias ~0.2mm (negligible), mean RMSE ~11.3mm, mean
  correlation ~0.83 (strong) - model tracks stations well overall.
- **Important pattern found**: bias vs. elevation plot shows bias stays near
  zero up to ~700-800m, then becomes increasingly **negative** above that
  (model **under-predicts** snow at high-elevation/mountain stations,
  reaching -100 to -180mm bias at the highest stations, ~1300-1500m).
  - Likely explanation: grid cell size (~10-14km at 0.125deg) averages
    together valley + peak elevations, so the model's effective elevation
    for a mountain grid cell is lower than a specific high-altitude station
    sitting at/near a peak - and snow accumulation is highly elevation-
    sensitive. This is a known general limitation of gridded models at this
    resolution in mountainous terrain, worth discussing in the thesis.
  - This nuances the supervisor's original highlands-vs-lowlands framing:
    the model gets the *spatial pattern* of where snow is right, but may be
    getting the *magnitude* systematically wrong at the highest elevations
    specifically.

**Status check against supervisor's original email:**
- [x] Station metadata table
- [x] Extract modelled SWE at each station
- [x] Calculate bias + other indicators
- [ ] Meeting in early September (or sooner - now have real intermediate
      results worth sharing, should message supervisor)
- [ ] UIS assignment paperwork (comes after supervisor discussion)
- [ ] ISIMIP3b SSP extension - explicitly a "possible future extension" in
      the email, gated on supervisor confirming SWE was written out in those
      runs. Not started, intentionally deferred.

**Next planned steps:** possibly break down bias by elevation band into a
clean summary table for the thesis; message supervisor with intermediate
results ahead of September meeting.

---

## Template for future entries
