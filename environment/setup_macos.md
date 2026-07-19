# Environment setup (macOS)

This project assumes you already have **Homebrew** and **mHM** installed. This file
covers the extra tools needed for the data analysis side (NetCDF handling, station
data processing, plotting).

## 1. Check what you already have

Run these in Terminal first — no harm in running them, they just print info:

```bash
brew --version
cdo --version
python3 --version
python3 -m pip list | grep -Ei "xarray|netCDF4|pandas|numpy|matplotlib|cartopy"
```

- If `cdo --version` prints a version → CDO is installed, skip step 2.
- If the `pip list` command prints nothing → none of those Python packages are
  installed yet, go to step 3.
- If it prints some but not all of them → only install the missing ones in step 3.

## 2. Install CDO (Climate Data Operators) via Homebrew, if missing

```bash
brew install cdo
```

CDO is what your supervisor used for `cdo timavg` — we'll use it for grid inspection,
long-term averages, and possibly point extraction at station coordinates.

## 3. Set up a Python virtual environment

Keeping this project's Python packages isolated avoids conflicts with other projects:

```bash
cd /path/to/aleyna-swe-thesis
python3 -m venv .venv
source .venv/bin/activate
pip install -r environment/requirements.txt
```

Whenever you come back to work on this project, re-activate with:

```bash
source .venv/bin/activate
```

(You'll see `(.venv)` appear in your terminal prompt when it's active.)

## Notes

- `cdo_refcard.pdf` (CDO reference card) is kept in `environment/` for quick command
  lookup while writing scripts.
- If `mHM` needs to be re-run (not just its outputs analyzed), that has its own
  build/run instructions from your supervisor — not duplicated here yet.
