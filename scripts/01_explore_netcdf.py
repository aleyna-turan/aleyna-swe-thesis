"""
01_explore_netcdf.py

Inspects the mHM SWE NetCDF file: variable names, dimensions, grid extent,
time range, and produces a long-term average map (like the supervisor's
`cdo timavg` example) as a sanity check.

Usage:
    source .venv/bin/activate
    python scripts/01_explore_netcdf.py

Input:
    data/raw/snowpack_czechia_19500101_20260531.nc

Output:
    results/figures/long_term_average_swe.png
    Prints a summary to the terminal.
"""

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import xarray as xr

REPO_ROOT = Path(__file__).resolve().parent.parent
RAW_DATA_DIR = REPO_ROOT / "data" / "raw"
FIG_OUT = REPO_ROOT / "results" / "figures" / "long_term_average_swe.png"


def find_nc_file() -> Path:
    """Supervisor names folders by the date data was shared, so the file
    may be directly in data/raw/ or nested in a dated subfolder
    (e.g. data/raw/20260622/...). Search recursively, and if there are
    several, pick the one with the largest file size (fullest time range)."""
    candidates = sorted(RAW_DATA_DIR.rglob("snowpack_czechia_*.nc"))
    candidates = [c for c in candidates if "timavg" not in c.name]
    if not candidates:
        print(f"[ERROR] No 'snowpack_czechia_*.nc' file found anywhere under {RAW_DATA_DIR}")
        print("Make sure the mHM NetCDF file has been placed in data/raw/ "
              "(in a dated subfolder is fine).")
        sys.exit(1)
    if len(candidates) > 1:
        print(f"[NOTE] Found {len(candidates)} matching files, using the "
              f"most recently modified one:")
        for c in candidates:
            print(f"  - {c.relative_to(REPO_ROOT)}")
        candidates.sort(key=lambda p: p.stat().st_mtime)
    return candidates[-1]


def main():
    NC_FILE = find_nc_file()

    print(f"Opening {NC_FILE.name} ...\n")
    ds = xr.open_dataset(NC_FILE)

    print("=" * 60)
    print("DATASET OVERVIEW")
    print("=" * 60)
    print(ds)
    print()

    print("=" * 60)
    print("VARIABLES")
    print("=" * 60)
    for name, var in ds.data_vars.items():
        units = var.attrs.get("units", "n/a")
        long_name = var.attrs.get("long_name", "n/a")
        print(f"- {name}: dims={var.dims}, shape={var.shape}, "
              f"units={units}, long_name={long_name}")
    print()

    print("=" * 60)
    print("COORDINATES / GRID")
    print("=" * 60)
    for coord_name in ds.coords:
        coord = ds.coords[coord_name]
        print(f"- {coord_name}: dtype={coord.dtype}, size={coord.size}")

    lat_name = next((c for c in ds.coords if c.lower() in
                      ("lat", "latitude", "y")), None)
    lon_name = next((c for c in ds.coords if c.lower() in
                      ("lon", "longitude", "x")), None)

    if lat_name and lon_name:
        lats = ds.coords[lat_name].values
        lons = ds.coords[lon_name].values
        print(f"\nLatitude range:  {lats.min():.4f} to {lats.max():.4f}")
        print(f"Longitude range: {lons.min():.4f} to {lons.max():.4f}")
        if len(lats) > 1:
            print(f"Latitude resolution:  {abs(lats[1] - lats[0]):.4f} deg")
        if len(lons) > 1:
            print(f"Longitude resolution: {abs(lons[1] - lons[0]):.4f} deg")
    else:
        print("\n[WARNING] Could not auto-detect lat/lon coordinate names.")
        print("Coordinate names found:", list(ds.coords))
    print()

    print("=" * 60)
    print("TIME RANGE")
    print("=" * 60)
    time_name = next((c for c in ds.coords if c.lower() in
                       ("time", "t")), None)
    if time_name:
        times = ds.coords[time_name].values
        print(f"Start: {times.min()}")
        print(f"End:   {times.max()}")
        print(f"Number of timesteps: {len(times)}")

        time_index = ds.indexes[time_name] if hasattr(ds, "indexes") else None
        if time_index is not None:
            diffs = np.diff(time_index.values).astype("timedelta64[D]")
            gap_days = diffs[diffs != np.timedelta64(1, "D")]
            if len(gap_days) > 0:
                print(f"[NOTE] Found {len(gap_days)} non-1-day gaps in the "
                      f"time series (this may be expected, e.g. leap days "
                      f"handled differently, or real missing periods).")
            else:
                print("Time series is a continuous daily sequence (no gaps).")
    else:
        print("[WARNING] Could not auto-detect a time coordinate.")
    print()

    print("=" * 60)
    print("LONG-TERM AVERAGE MAP")
    print("=" * 60)

    candidate_vars = [
        name for name, var in ds.data_vars.items()
        if not name.endswith("_bnds") and set(var.dims) >= {time_name, lat_name, lon_name}
    ]
    if not candidate_vars:
        print("[ERROR] Could not find a (time, lat, lon) data variable.")
        sys.exit(1)
    main_var = candidate_vars[0]
    print(f"Using variable: {main_var}")

    if time_name:
        avg = ds[main_var].mean(dim=time_name, skipna=True)
    else:
        avg = ds[main_var]

    FIG_OUT.parent.mkdir(parents=True, exist_ok=True)
    from matplotlib.colors import LogNorm

    plt.figure(figsize=(8, 6))
    vmin = max(float(avg.min()), 0.1)
    vmax = float(avg.max())
    avg.plot(cmap="jet_r", norm=LogNorm(vmin=vmin, vmax=vmax))
    plt.title(f"Long-term average {main_var} (1950-2026)")
    plt.tight_layout()
    plt.savefig(FIG_OUT, dpi=150)
    print(f"\nSaved long-term average map to: {FIG_OUT}")

    ds.close()


if __name__ == "__main__":
    main()
