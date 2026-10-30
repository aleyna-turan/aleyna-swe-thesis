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
