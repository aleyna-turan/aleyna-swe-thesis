# Data folder

**Nothing in `raw/` or `stations/` is committed to git** (see `.gitignore`) — these
folders exist locally on your machine only, so the repo stays small and doesn't
redistribute your supervisor's unpublished data or CHMI station data without
permission.

## `data/raw/`

Place the mHM NetCDF simulation file(s) here, e.g.:

```
data/raw/snowpack_czechia_19500101_20260531.nc
```

Source: supervisor's SharePoint (`czuvpraze-my.sharepoint.com`), folder named by the
date it was shared.

- Coverage: 1950–2025+ (see filename for exact end date), daily time step
- Resolution: 0.125° x 0.125°, cropped to Czechia
- Forcing: E-OBS meteorological data

## `data/stations/`

Place CHMI station SWE data here once obtained (format TBD — see `docs/notes.md`
for progress on this).

Expected eventual contents:
- Raw per-station SWE time series files (however CHMI provides them — CSV, or via
  a query tool)
- The generated station metadata table (this one *is* small and could optionally be
  tracked in `results/tables/` instead once finalized — see step 2 of the pipeline)

## Why raw data isn't in git

- The mHM file is tens of MB and will grow as the time series extends — not something
  git handles efficiently.
- It's the supervisor's unpublished research output — shouldn't be redistributed via
  a repo without explicit permission, even a private one.
