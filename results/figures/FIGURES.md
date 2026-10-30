# Figures — what each one shows, how it was computed, and why it exists

Figures are grouped by the **question they answer**, not by study period. Each
group is produced by one script, and each figure carries its own title, method
subtitle, legend and footnote inside the image, so it can be read without this
file. This file adds the reasoning that does not fit in a caption.

Conventions used everywhere:

| Convention | Value |
|---|---|
| Hydrological year | 1 Nov to 31 Oct, named after the year it **ends** in (HY 2025 = 1 Nov 2024 to 31 Oct 2025). Czech / CHMI convention. |
| Period label | Calendar span. `1976-2025` = 1 Nov 1976 to 31 Oct 2025 = HY 1977–2025, 49 hydrological years. |
| Snow threshold | SWE ≥ 5 mm counts as snow cover (see `04_threshold_tests/`). |
| Sampling | SWE is measured **weekly**, so all timing metrics have ~7-day resolution. |
| Model values | Sampled on the **same dates** as the observations, never from the full daily model output. |
| Observed colour | Blue, drawn thick and pale, underneath. |
| Model colour | Orange, drawn thin and dashed, on top. |

Drawing the model on top of a thicker observed band is deliberate: where the
two series agree exactly, the dashed line sits visibly inside the band instead
of one line disappearing under the other. **A single visible line means
agreement, not missing data.**

---

## 01_data_screening/ — which stations can be trusted?

Produced by `scripts/06_station_qc_selection.py`.

### `availability_<period>.png`

**Shows.** A heatmap of the measurement record of every station that *passed*
screening for that period. One row per station (sorted by record length), one
column per snow season. Colour is the number of observations recorded in that
season.

**Computed.** Observation dates are taken from the paired SWE table, restricted
to Nov–Apr, and assigned to a snow season (Nov/Dec to that year, Jan–Apr to the
previous). Sampling frequency is homogenised first: for each season, only that
season's dominant measurement weekday across the whole network is kept. The
cell value is then the count of surviving observations. A full weekly season is
about 26.

**Why.** To show that the surviving network really is continuous, and to make
the tolerated gaps visible rather than hidden behind a pass/fail flag. Dark
streaks are seasons that fell below the activity threshold but were within the
allowed gap length. The columns around 1994–2000 are visibly thinner across
many stations, which is a network-wide effect and not a per-station fault.

**Watch out.** Only *selected* stations are shown, so the figure cannot be used
to judge how many stations were rejected — that is in
`tables/01_data_screening/qc_rejection_summary.csv`.

---

## 02_model_evaluation/ — how well does mHM reproduce the observations?

Produced by `scripts/07_seasonal_statistics.py`.

### `annual_series_<period>.png`

**Shows.** Two stacked panels against hydrological year: network median peak
SWE (top), and network median snow season length (bottom), observed and
modelled.

**Computed.** For each station and hydrological year: peak SWE is the largest
SWE value in that year; snow season length is the number of days from the first
to the last observation with SWE ≥ 5 mm. The network value is the **median**
across all selected stations for that year — median, not mean, so that a few
very snowy mountain stations do not dominate.

**Why.** It is the overall sanity check: does the model track year-to-year
variability, and does it drift relative to the observations over time? It also
shows the big snow years (HY 2006) and the snow-free ones (HY 2007, where the
median station had no snow season at all at the 5 mm threshold) as shared
features of both series.

**Watch out.** Snow season length can only take values that are **multiples of
7 days**, because SWE is sampled weekly. Exact ties between observed and
modelled are therefore common — 15 of 39 years in the 1986–2025 figure. The
flat steps and the frequent perfect overlap are a property of the sampling, not
of the model.

### `bias_vs_elevation_<period>.png`

**Shows.** One point per station: the station's elevation against the model's
mean peak SWE bias, with a Theil–Sen fit and the Mann–Kendall test result. This
is the **headline figure of the thesis**.

**Computed.** Bias at a station = mean over all its hydrological years of
(modelled peak SWE − observed peak SWE). Negative means the model holds less
snow than was measured. The line is a Theil–Sen fit: the median of the slopes
of all station pairs. The test is Mann–Kendall (Kendall's tau).

**Why.** Theil–Sen and Mann–Kendall are used rather than ordinary least squares
because a few isolated summit stations have extreme bias and would dominate a
least-squares fit, which squares residuals. The median of pairwise slopes is
unaffected by how far away those stations sit. The two most extreme stations are
labelled by name so they are not anonymous points: for 1976–2025 these are
Bedřichov and Lysá hora.

**Watch out.** The relationship is modest in strength (tau ≈ −0.18) but highly
significant across 235 stations, and it **survives removing the three stations
above 1000 m**. Report the station-level fit, not the elevation-band means — the
bands above 800 m contain only 4 and 3 stations and are dominated by individual
outliers.

---

## 03_trends/ — is the snow changing, and does the model capture it?

Produced by `scripts/07_seasonal_statistics.py`.

### `trend_vs_elevation_<period>.png`

**Shows.** One point per station: elevation against that station's observed
peak SWE trend in mm per decade, split by whether the trend is statistically
significant.

**Computed.** For each station, the Theil–Sen slope of its annual peak SWE
against hydrological year, multiplied by 10 to give mm per decade. Significance
is Mann–Kendall at alpha = 0.05. Only stations with at least 80% of the
period's hydrological years present are included, so a short or gappy record
cannot produce a spurious trend.

**Why.** It separates two different statements that are easy to confuse: *snow
is declining* (most points below zero) and *the decline is stronger at
altitude* (points spreading further below zero towards the right). It also
shows honestly how many stations fail to reach significance.

**Watch out.** This figure is about the **observations**, not the model. The
comparison of observed and modelled trends is in
`tables/03_trends/trend_summary_<period>.csv`, including the correlation
between the two sets of station slopes.

---

## 04_threshold_tests/ — is the duration result real or an artefact?

Produced by `scripts/08_duration_bias_diagnosis.py`.

### `duration_threshold_<period>.png`

**Shows.** The model's snow-duration bias recomputed at thresholds of 0.5, 1,
2, 5, 10 and 20 mm SWE, for all four candidate periods.

**Computed.** At each threshold, snow onset, melt-out and season length are
recalculated from scratch for every station-year, for both series. The plotted
value is the mean of (modelled length − observed length) across all
station-years.

**Why.** This figure exists because an earlier version of the analysis, run at
a 1 mm threshold, appeared to show mHM simulating snow seasons 13 to 25 days
too long — which would have been a major finding. A bias that shrinks towards
zero as the threshold rises is an **artefact** of where the threshold was
drawn; a bias that stays flat would be a real model error. The bias collapses
to about −0.7 days at 5 mm and stays there, so the finding was not real, and
the analysis threshold was changed to 5 mm.

**Watch out.** This is the figure that justifies `SNOW_THRESHOLD_MM = 5.0` in
script 07. Cite it whenever the threshold is stated.

### `false_snow_seasonality_<period>.png`

**Shows.** For each month in hydrological-year order, the share of observation
dates where the model reports snow and the station does not ("false snow"), and
the reverse ("missed snow").

**Computed.** Every observation date is classified at a 1 mm threshold into
agree-snow, agree-bare, false-snow or missed-snow, then aggregated by calendar
month.

**Why.** It identifies the *mechanism* behind the threshold artefact rather
than just demonstrating that it exists. False snow outnumbers missed snow by
more than four to one, the median modelled value on false-snow dates is about
5 mm — right at the boundary — and it is concentrated in early winter (12.4% of
December dates, 8.3% of January). That is the signature of the model laying
down a thin snow cover that the CHMI protocol does not register, since SVH is
measured when snow depth exceeds 3 cm.

---

## archive/ — superseded, kept for the record

Not regenerated by any current script. Kept so that earlier results in
`docs/notes.md` can still be traced, and clearly separated so they are never
mistaken for current output.

| File | Origin | Why superseded |
|---|---|---|
| `bias_vs_elevation_2026-07_R_superseded.png` | `scripts/04_bias_analysis.R`, July 2026 | Used all 1,147 stations with no screening, no hydrological-year handling and no weekday homogenisation. Replaced by `02_model_evaluation/bias_vs_elevation_<period>.png`. |
| `long_term_average_swe.png` | `scripts/01_explore_netcdf.py`, July 2026 | Exploratory map of the raw model grid, made before the station comparison existed. Still useful as a picture of the model domain. |
