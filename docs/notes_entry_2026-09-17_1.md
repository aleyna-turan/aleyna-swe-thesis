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
