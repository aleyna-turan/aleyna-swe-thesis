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

## Template for future entries

```
## YYYY-MM-DD — short title

- What was done
- What was decided and why
- Open questions / blockers
```
