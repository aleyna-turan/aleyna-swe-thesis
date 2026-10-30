# Tables — what each one contains and which script writes it

Grouped by the question they answer, matching `results/figures/`. Conventions
(hydrological year, period labels, 5 mm snow threshold, weekly sampling,
model sampled on observation dates) are documented in
[`../figures/FIGURES.md`](../figures/FIGURES.md).

**Which tables to cite.** Most files here are working output. The ones that
carry results are marked **cite** below; the rest exist so the analysis can be
audited station by station.

---

## 00_inputs/ — the data everything else is built from

| File | Written by | Contents |
|---|---|---|
| `station_metadata.csv` | `02_station_metadata.py` | 1,147 CHMI stations: id, name, lat/lon, elevation, record start and end, valid and missing record counts. |
| `modelled_vs_observed_swe.csv` | `03_extract_modelled_swe.py` | The paired dataset: one row per station per date with a real observation, carrying observed and modelled SWE. ~1.98 M rows. Gitignored (74 MB). |
| `modelled_vs_observed_swe.csv.gz` | compressed copy | The version committed to the repo (9.8 MB). Scripts accept either. |

**Important.** Modelled SWE was extracted **only** on dates where the station
has a real observation. There are therefore no gaps in this table — a missing
observation produced no row at all, rather than a row with a missing value.
This is why screening counts rows per season instead of looking for NaNs, and
why the comparison is fair: the model is only ever evaluated where the
observation exists.

---

## 01_data_screening/ — which stations measured regularly

Written by `scripts/06_station_qc_selection.py`.

| File | Contents |
|---|---|
| **`qc_selection_overview.csv`** | **cite** — final station count per period, with elevation distribution and worst coverage and gap in the selected set. |
| **`qc_rejection_summary.csv`** | **cite** — how many candidate stations failed on each criterion, per period. |
| `qc_selected_<period>.csv` | The station list used by scripts 07 and 08 for that period. |
| `qc_diagnostics_<period>.csv` | Every candidate station with all diagnostics, one boolean column per criterion, and a `fail_reasons` string. Use this to look up why any individual station was dropped. |
| `season_presence_<period>.csv` | Station × season matrix of observation counts for the selected stations. The source of the availability heatmap. |
| `season_measurement_weekday.csv` | The dominant measurement weekday for each season across the network. **This is where the 1970 Friday-to-Monday switch is visible.** |

The six screening criteria, all of which must pass: coverage ≥ 90% of seasons
active; no gap longer than 3 consecutive inactive seasons; active within 2
seasons of each end of the period; median ≥ 15 observations per active season;
coefficient of variation of observations per season ≤ 0.60. A season counts as
active with at least 10 observations.

---

## 02_model_evaluation/ — how well mHM reproduces the observations

Written by `scripts/07_seasonal_statistics.py`.

| File | Contents |
|---|---|
| **`period_comparison.csv`** | **cite** — one row per candidate period: station counts, station-years, observed climatology, model bias and correlation, and headline trends. This is the table the period choice was made from. |
| **`elevation_bands_<period>.csv`** | **cite with care** — metrics and bias averaged within elevation bands. The bands above 800 m hold only 4 and 3 stations and are dominated by individual outliers; prefer the station-level regression. |
| `annual_series_<period>.csv` | Network median and mean of every metric, per hydrological year. Source of the annual series figure. |
| `model_eval_<period>.csv` | Per station: bias, MAE, RMSE, correlation and observed/modelled means for each of the nine metrics. |
| `seasonal_metrics_<period>.csv` | The raw per-station, per-hydrological-year metrics everything else is aggregated from. Gitignored (2.5 MB each); regenerates in about a minute. |

The nine metrics, computed for observed and modelled side by side: mean SWE
over the hydrological year; mean SWE over a fixed 1 Dec–31 Mar window; mean SWE
on snow days only; peak SWE; day-of-hydrological-year of the peak; snow onset;
melt-out; snow season length; weeks with snow; and SWE-days (integrated
storage). Timing metrics use day-of-hydrological-year with 1 = 1 November, so
dates around New Year do not wrap.

---

## 03_trends/ — is the snow changing

Written by `scripts/07_seasonal_statistics.py`.

| File | Contents |
|---|---|
| **`trend_summary_<period>.csv`** | **cite** — per metric: median observed and modelled slope per decade, counts of significantly increasing and decreasing stations, and the correlation between observed and modelled station slopes. |
| `trends_<period>.csv` | Per station per metric: Theil–Sen slope, Kendall's tau, p-value and the number of years used, for both series. |

Slopes are Theil–Sen; significance is Mann–Kendall at alpha = 0.05. Stations
need at least 80% of the period's hydrological years to enter the trend
analysis.

---

## 04_threshold_tests/ — is the duration result an artefact

Written by `scripts/08_duration_bias_diagnosis.py`.

| File | Contents |
|---|---|
| **`duration_threshold_sensitivity.csv`** | **cite** — duration, onset and melt-out bias at six thresholds from 0.5 to 20 mm, for all four periods. The evidence that the 1 mm duration bias was an artefact. |
| `duration_bias_by_elevation.csv` | The same sensitivity, broken down by elevation band, for the main period. |
| `snow_contingency_<period>.csv` | Counts and median values for agree-snow, agree-bare, false-snow and missed-snow dates, at each threshold. |
| `false_snow_by_month_<period>.csv` | The same classification aggregated by calendar month. Source of the seasonality figure. |

---

## archive/ — superseded

| File | Origin | Why superseded |
|---|---|---|
| `station_bias_metrics.csv` | `scripts/04_bias_analysis.R`, July 2026 | Per-station bias over all 1,147 stations with no screening, no hydrological-year handling and no weekday homogenisation. Replaced by `02_model_evaluation/model_eval_<period>.csv`. |
