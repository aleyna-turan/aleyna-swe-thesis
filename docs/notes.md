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

## 2026-09-17 — Study period selection + seasonal statistics (Steps 6–7)

Goal: decide which observation period to use (1956/1966/1976/1986 → 2025) by
first checking which stations actually measured *regularly*, then computing
seasonal SWE statistics for every period so they can be compared.

### Two data quirks found (both matter for the methodology section)

**1. CHMI changed the weekly measurement day from Friday to Monday during 1970.**
- Pre-1970: 92–94% of SVH observations fall on Fridays.
- From the 1970/71 season: 97–99% fall on Mondays.
- Found this the hard way — a first version of the screening homogenised to
  "Mondays only" and returned **zero** stations for 1956 and 1966.
- Fix: for each season/hydrological year, keep only that year's modal weekday
  across the network (`FREQUENCY_HOMOGENISATION = "dominant"`). Logged per
  season in `results/tables/season_measurement_weekday.csv`.
- Why this matters beyond the switch: some station-eras reported more often
  than weekly (up to 166 obs in a single Nov–Apr window vs ~26 expected).
  Without homogenisation those stations and decades are over-weighted in
  every mean, bias and RMSE.

**2. Snow outside Nov–Apr is real and non-trivial.**
- 1239 non-zero SWE observations fall in May–Oct (961 in October, 263 in May),
  at 430 stations, max 686 mm.
- So all seasonal metrics are computed over the **full hydrological year**,
  not a Nov–Apr window. Nov–Apr is used only for the availability screening.

### Hydrological year convention

Czech HY = 1 Nov – 31 Oct, named after the year it **ends** in
(HY 2025 = 1 Nov 2024 → 31 Oct 2025). Every table from script 07 carries a
`hydro_year` column in this convention. Period labels keep the calendar span:
"1976-2025" = 1 Nov 1976 → 31 Oct 2025 = HY 1977–2025 (49 hydrological years).
Timing metrics use day-of-hydrological-year (DOHY, 1 = 1 Nov) so dates around
New Year don't wrap.

### Step 6 — station screening (`scripts/06_station_qc_selection.py`)

Coverage alone is not enough: a station measuring every second year, or one
that switched sampling frequency mid-record, can still show "90% coverage".
Six criteria, all must pass:

| Criterion | Setting | Catches |
|---|---|---|
| coverage | ≥ 0.90 of seasons active | too few active seasons |
| max gap | ≤ 3 consecutive dead seasons | "every second year" stations |
| starts early / ends late | active within 2 seasons of each end | records covering only one end of the period |
| median obs | ≥ 15 per active season | technically active but sparse seasons |
| CV of obs | ≤ 0.60 | wildly varying measurement intensity |

A season counts as active with ≥ 10 observations. Rejection reasons are written
per station to `qc_diagnostics_<period>.csv`.

**Result:**

| Period | HY | Selected | Median elev. | Max elev. | >800 m | >1000 m |
|---|---|---|---|---|---|---|
| 1956–2025 | 1957–2025 | 62 | 420 m | **780 m** | **0** | **0** |
| 1966–2025 | 1967–2025 | 187 | 453 m | 1322 m | 5 | 3 |
| 1976–2025 | 1977–2025 | 235 | 458 m | 1322 m | 7 | 3 |
| 1986–2025 | 1987–2025 | 303 | 477 m | 1322 m | 13 | 4 |

### Step 7 — seasonal statistics (`scripts/07_seasonal_statistics.py`)

Per station × hydrological year, for observed and mHM side by side: mean SWE
(whole HY / DJFM window / snow-days only), peak SWE, peak timing, snow onset,
melt-out, season length, weeks with snow, SWE-days. Then per station:
bias/MAE/RMSE/correlation, plus Theil–Sen slope with Mann–Kendall test for
every metric.

**Period comparison:**

| Period | Stations | Station-years | Obs peak SWE | Peak bias | r | Peak trend (mm/dec) | % sig. declining |
|---|---|---|---|---|---|---|---|
| 1956–2025 | 62 | 4168 | 40.8 mm | +3.3 | 0.78 | −0.43 | 9.7% |
| 1966–2025 | 187 | 10820 | 41.3 mm | −3.3 | 0.80 | −1.17 | 20.9% |
| 1976–2025 | 235 | 11276 | 43.6 mm | −4.8 | 0.81 | −3.25 | 40.4% |
| 1986–2025 | 303 | 11597 | 45.2 mm | −5.1 | 0.83 | −2.75 | 19.8% |

**Trends, 1976–2025, 235 stations (median per decade):**

| Metric | Observed | Sig. neg / pos | mHM | Slope agreement r |
|---|---|---|---|---|
| Peak SWE | −3.25 mm | 95 / 1 | −3.04 mm | 0.67 |
| Snow duration | −9.33 d | 158 / 0 | −8.24 d | 0.51 |
| Snow onset | +3.85 d (later) | 0 / 71 | +3.23 d | 0.54 |
| Melt-out | −6.36 d (earlier) | 126 / 0 | −5.26 d | 0.39 |
| Weeks with snow | −1.00 | 179 / 1 | −1.11 | 0.80 |

Coherent signal: less peak snow, later onset, earlier melt, ~9 fewer snow days
per decade — and mHM reproduces the direction and roughly the magnitude of all
of it.

**New finding — mHM overestimates snow season length by +13 to +25 days** in
every period and across all elevation bands. Larger and more systematic than
the peak SWE bias found in July. Probable cause: the model carries small
continuous SWE values at the shoulders of the season where observations record
none. Still to test: whether raising `SNOW_THRESHOLD_MM` from 1 to 5 mm
removes it (threshold artefact) or not (real model behaviour).

**Elevation bias (1976–2025)** — the July pattern holds and sharpens:

| Band | Stations | Obs peak SWE | Peak SWE bias |
|---|---|---|---|
| 0–400 m | 83 | 24.5 mm | +1.3 mm |
| 400–600 m | 113 | 36.8 mm | −3.3 mm |
| 600–800 m | 32 | 85.4 mm | −21.9 mm |
| 800–1000 m | 4 | 122.9 mm | +25.7 mm |
| >1000 m | 3 | 261.4 mm | −91.3 mm |

Caveat to state explicitly in the thesis: the two highest bands rest on 4 and 3
stations. The mountain under-prediction is a real and physically plausible
result, but the sample supporting it is small in **every** candidate period —
this is a network limitation, not something a different period choice fixes.

### Period decision

Leaning **1976–2025** as the main analysis period: 235 stations × 49 HY =
11,276 station-years, the strongest trend signal (40% of stations
significantly declining), and the same mountain coverage as 1966. 1956–2025
is the one to avoid as a primary period — 62 stations, none above 780 m, and
its weak trend largely reflects which stations survived screening rather than
the longer window. Possible structure: 1976–2025 as main period, 1956–2025 as
a supporting long-record lowland subset.

**Still to do:** confirm period choice with supervisor; test the snow-duration
bias against a 5 mm threshold; check whether observed and modelled curves
diverge over time in `annual_maxswe_1976-2025.png`.

---

---

## 2026-10-06 — Snow-duration bias diagnosed; threshold changed to 5 mm (Step 8)

**Headline: the "+13 to +25 day snow duration overestimate" reported in the
2026-09-17 entry was largely a threshold artefact. That entry is corrected
below.** This is why the test was run before writing anything up.

### The question

Script 07 (at `SNOW_THRESHOLD_MM = 1.0`) found mHM simulating snow seasons
13–25 days longer than observed, in every period and every elevation band.
That was larger and more systematic than the peak-SWE bias, so it needed to be
established as real model behaviour before it could be claimed as a result.

Two candidate explanations:
1. mHM genuinely holds snow too long — a real model error.
2. mHM carries small SWE values at the shoulders of the season that the 1 mm
   threshold counts as "snow cover" but the CHMI measurement never registers —
   an artefact of where the threshold was drawn.

### How it was tested — `scripts/08_duration_bias_diagnosis.py`

**Part A, threshold sensitivity.** Recompute onset, melt-out, duration and
weeks-with-snow at 0.5, 1, 2, 5, 10 and 20 mm, for all four periods. A
threshold artefact shrinks towards zero as the threshold rises; a real model
error persists.

**Part B, date-by-date mechanism.** Classify every observation date into
agree_snow / agree_bare / false_snow (model says snow, observation doesn't) /
missed_snow (observation says snow, model doesn't), and look at both the
magnitude of the modelled values in the false_snow cases and where in the
season they fall.

### Result — it is a threshold artefact

Snow duration bias, 1976–2025 (same pattern in all four periods):

| Threshold | Obs duration | mHM duration | Bias | Onset bias | Melt-out bias |
|---|---|---|---|---|---|
| 0.5 mm | 75.5 d | 98.2 d | **+22.7 d** | −10.0 d | +12.7 d |
| 1 mm | 75.2 d | 88.6 d | **+13.3 d** | −7.1 d | +6.2 d |
| 2 mm | 74.1 d | 80.0 d | +5.9 d | −4.0 d | +1.9 d |
| **5 mm** | 69.4 d | 68.8 d | **−0.7 d** | −1.4 d | −2.1 d |
| 10 mm | 60.6 d | 58.9 d | −1.7 d | −0.9 d | −2.6 d |
| 20 mm | 51.2 d | 51.3 d | +0.1 d | −0.2 d | −0.1 d |

The bias collapses to essentially zero at 5 mm and stays there.

Part B confirms the mechanism directly. At 1 mm, 3.38% of all observation
dates are false_snow against only 0.75% missed_snow, and the **median modelled
value on false_snow dates is 5.12 mm** — sitting right at the boundary. They
are strongly seasonal: 12.4% of December dates and 8.3% of January dates,
falling to 0.1% in May and zero through summer.

So mHM lays down a thin simulated snow cover in early winter that the CHMI
protocol does not register (SVH is measured when snow depth exceeds 3 cm), and
that thin layer inflates the apparent season length at low thresholds.

### Decision

`SNOW_THRESHOLD_MM` changed from 1.0 to **5.0** in `scripts/07_seasonal_statistics.py`
and all tables regenerated. 5 mm is the lowest threshold at which the
comparison is stable, and it is physically closer to what the measurement
protocol actually captures.

Unaffected by the change (no threshold used): mean SWE, DJFM mean, peak SWE,
SWE-days, and all peak-SWE biases and trends.
Affected: onset, melt-out, duration, weeks with snow, mean SWE on snow days.

### What survives, and is worth writing up

**1. Trends are robust to the threshold choice.** Snow duration trend went from
−9.33 d/decade (1 mm, 158/235 stations significant) to −7.78 d/decade (5 mm,
119/235 significant). Magnitude softened, signal unchanged. Worth stating
explicitly in the results — it pre-empts the obvious question.

Trends, 1976–2025, 235 stations, median per decade, α = 0.05:

| Metric | Observed | Sig. declining | mHM | Slope agreement r |
|---|---|---|---|---|
| Peak SWE | −3.25 mm | 95 | −3.04 mm | 0.67 |
| Snow duration | −7.78 d | 119 | −8.51 d | 0.32 |
| Weeks with snow | −0.80 | 156 | −0.96 | 0.76 |
| SWE-days | −116.9 mm·d | 143 | −121.2 | 0.76 |
| Melt-out | −5.50 d earlier | 92 | −4.74 d | 0.17 |
| Snow onset | +2.93 d later | 31 sig. later | +4.00 d | 0.46 |
| Mean SWE (HY) | −0.38 mm | 152 | −0.40 mm | 0.77 |

**2. The elevation result is now coherent across two independent variables.**

| Band | Stations | Obs peak SWE | Peak SWE bias | Obs duration | Duration bias |
|---|---|---|---|---|---|
| 0–400 m | 83 | 24.5 mm | +1.3 mm | 49.6 d | +2.8 d |
| 400–600 m | 113 | 36.8 mm | −3.3 mm | 62.8 d | +0.1 d |
| 600–800 m | 32 | 85.4 mm | −21.9 mm | 102.4 d | −8.3 d |
| 800–1000 m | 4 | 122.9 mm | +25.7 mm | 125.1 d | −7.2 d |
| >1000 m | 3 | 261.4 mm | −91.3 mm | 194.0 d | −26.9 d |

In the lowlands mHM is near-perfect on both peak magnitude and season length.
Above 600 m it under-predicts both — less snow **and** a shorter season. That
is physically coherent as one mechanism with two symptoms: insufficient
accumulation at elevation, therefore earlier melt-out.

Caveats to state explicitly: the top two bands rest on 4 and 3 stations, and
this is a network limitation present in **every** candidate period, not
something a different period choice fixes. The +25.7 mm peak bias at
800–1000 m is the odd one out and may be driven by a single station — still to
check in `model_eval_1976-2025.csv`.

### Open items

- Check whether one station drives the 800–1000 m peak-SWE anomaly.
- Confirm the 1976–2025 period choice with supervisor.
- Consider whether the shoulder-season mismatch deserves its own short results
  subsection as a measurement-definition issue rather than a model error.

---

---

## 2026-10-30 — Repository reorganised, every figure documented

The results folders had grown to 16 figures and 48 tables, all flat, with
names that gave the period but not the question. Reorganised both by the
question each output answers, so the structure matches how the thesis is
argued rather than how the scripts happened to run.

```
results/figures/ and results/tables/
  00_inputs/            station metadata, the paired SWE table   (tables only)
  01_data_screening/    which stations measured regularly        (script 06)
  02_model_evaluation/  how well mHM reproduces observations     (script 07)
  03_trends/            is the snow changing                     (script 07)
  04_threshold_tests/   is the duration result an artefact       (script 08)
  archive/              superseded July output, kept for the record
```

Scripts 06, 07 and 08 were edited to write into these folders, so the
structure is permanent rather than something that has to be re-tidied after
every run. All three were rerun and reproduce the previous numbers exactly
(62 / 187 / 235 / 303 selected stations; peak SWE bias −4.84 mm and r = 0.807
for 1976–2025).

**Every figure now carries its own explanation inside the image**: a title
saying what it shows, a subtitle saying how the plotted quantity was computed,
a legend, and a footnote on how to read it and what the limitations are. The
figures are meant to be understandable without the thesis text.
`results/figures/FIGURES.md` and `results/tables/TABLES.md` document each
output in full — what, how, why, and what to watch out for.

### Figure changes worth noting

- **Overlapping lines were unreadable.** In the annual series figure the
  observed and modelled lines coincide exactly in many years, and with two
  equal-weight lines the observed one simply vanished under the model — it
  looked like missing data. Observed is now drawn as a thick pale band with
  the model thin and dashed on top, so agreement is visible as a dashed line
  inside the band. The caption states this explicitly, and reports how many
  years are exact ties (15 of 39 for 1986–2025).
- **Why exact ties happen at all**: snow season length is measured between
  weekly sampling dates, so it can only be a multiple of 7 days. The flat
  steps in that panel are a property of the sampling, not of the model.
- **`annual_maxswe_*` renamed to `annual_series_*`** — the figure shows peak
  SWE *and* snow duration, and the old name mentioned only one.
- **New figure `bias_vs_elevation_<period>.png`**, generated by script 07.
  This is the headline result and until now it existed only as a July R plot
  made before screening and hydrological-year handling existed. The new
  version carries the Theil–Sen fit, the Mann–Kendall statistics, and names
  the two extreme stations (Bedřichov, Lysá hora) rather than leaving them as
  anonymous points.

### Also fixed

- The September and October notes entries had never been merged into this
  file — they were sitting loose in `docs/` as separate files. Merged;
  originals moved to `docs/archive/`.
- `scripts/05_station_availability.py` does not exist and never did. It was a
  first-pass availability check that was superseded by script 06 before being
  saved, which is why the numbering jumps 04 → 06. Recorded in the README so
  the gap is not mistaken for a lost file.

**Still open:** confirm the 1976–2025 period with the supervisor; check
whether observations-per-season declines with elevation among the selected
stations; decide on the degree-day benchmark experiment.

---

## Template for future entries
