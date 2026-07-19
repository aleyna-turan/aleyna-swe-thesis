"""
03_extract_modelled_swe.py

For each station, finds the nearest mHM grid cell and extracts the
modelled SWE value ONLY on the dates that station actually has an
observation (not every day in 1950-2026) - this is all we need for the
bias comparison in the next step, and keeps the output a manageable size.

Usage:
    source .venv/bin/activate
    python scripts/03_extract_modelled_swe.py

Inputs:
    data/raw/**/snowpack_czechia_*.nc      (the mHM file, found automatically)
    results/tables/station_metadata.csv    (from 02_station_metadata.py)
    data/stations/svh_raw/*.csv            (raw CHMI station files, for their
                                             actual observation dates)

Outputs:
    results/tables/modelled_vs_observed_swe.csv
        columns: station_id, date, observed_swe_mm, modelled_swe_mm
    results/tables/stations_outside_grid.csv
        stations whose coordinates fall outside the mHM grid (excluded)
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

REPO_ROOT = Path(__file__).resolve().parent.parent
RAW_DATA_DIR = REPO_ROOT / "data" / "raw"
STATION_METADATA = REPO_ROOT / "results" / "tables" / "station_metadata.csv"
SVH_RAW_DIR = REPO_ROOT / "data" / "stations" / "svh_raw"
OUT_COMBINED = REPO_ROOT / "results" / "tables" / "modelled_vs_observed_swe.csv"
OUT_EXCLUDED = REPO_ROOT / "results" / "tables" / "stations_outside_grid.csv"


def find_nc_file() -> Path:
    candidates = sorted(RAW_DATA_DIR.rglob("snowpack_czechia_*.nc"))
    candidates = [c for c in candidates if "timavg" not in c.name]
    if not candidates:
        print(f"[ERROR] No 'snowpack_czechia_*.nc' file found under {RAW_DATA_DIR}")
        sys.exit(1)
    candidates.sort(key=lambda p: p.stat().st_mtime)
    return candidates[-1]


def load_station_observations(station_id: str) -> pd.DataFrame:
    """Load a station's actual SWE observations (date + value), keeping
    only rows with a real (non-missing) value."""
    svh_path = SVH_RAW_DIR / f"dly-{station_id}-SVH.csv"
    if not svh_path.exists():
        return pd.DataFrame(columns=["date", "observed_swe_mm"])

    df = pd.read_csv(svh_path)
    df.columns = [c.strip() for c in df.columns]
    df["DT"] = pd.to_datetime(df["DT"], utc=True).dt.tz_localize(None)

    is_missing = df["QUALITY"].eq(4.0) | df["VALUE"].isna()
    df = df.loc[~is_missing, ["DT", "VALUE"]]
    df = df.rename(columns={"DT": "date", "VALUE": "observed_swe_mm"})
    df["date"] = df["date"].dt.normalize()
    return df


def main():
    nc_file = find_nc_file()
    print(f"Opening {nc_file.name} ...")
    ds = xr.open_dataset(nc_file)
    model_times = pd.to_datetime(ds["time"].values).normalize()

    lat_min, lat_max = float(ds["lat"].min()), float(ds["lat"].max())
    lon_min, lon_max = float(ds["lon"].min()), float(ds["lon"].max())
    print(f"Grid extent: lat [{lat_min:.4f}, {lat_max:.4f}], "
          f"lon [{lon_min:.4f}, {lon_max:.4f}]")

    print(f"\nLoading station metadata from {STATION_METADATA.name} ...")
    stations = pd.read_csv(STATION_METADATA)
    print(f"Total stations: {len(stations)}")

    inside_mask = (
        (stations["latitude"] >= lat_min) & (stations["latitude"] <= lat_max) &
        (stations["longitude"] >= lon_min) & (stations["longitude"] <= lon_max)
    )
    stations_inside = stations[inside_mask].copy()
    stations_outside = stations[~inside_mask].copy()

    print(f"\nStations inside mHM grid domain:  {len(stations_inside)}")
    print(f"Stations OUTSIDE mHM grid domain: {len(stations_outside)}")

    if len(stations_outside) > 0:
        OUT_EXCLUDED.parent.mkdir(parents=True, exist_ok=True)
        stations_outside.to_csv(OUT_EXCLUDED, index=False)
        print(f"  -> saved list to {OUT_EXCLUDED}")
        print("  These stations are excluded from extraction below.")

    if len(stations_inside) == 0:
        print("[ERROR] No stations fall inside the grid domain, nothing to extract.")
        sys.exit(1)

    print(f"\nExtracting modelled SWE for {len(stations_inside)} stations "
          f"(nearest grid cell, only on dates with real observations)...")

    all_rows = []
    n_no_obs = 0
    for i, (_, row) in enumerate(stations_inside.iterrows(), start=1):
        station_id = row["station_id"]
        lat = row["latitude"]
        lon = row["longitude"]

        obs = load_station_observations(station_id)
        if obs.empty:
            n_no_obs += 1
            continue

        series = ds["snowpack"].sel(lat=lat, lon=lon, method="nearest")
        series_df = pd.DataFrame({
            "date": model_times,
            "modelled_swe_mm": series.values,
        })

        merged = obs.merge(series_df, on="date", how="left")
        merged.insert(0, "station_id", station_id)
        all_rows.append(merged)

        if i % 100 == 0 or i == len(stations_inside):
            print(f"  processed {i}/{len(stations_inside)}")

    if n_no_obs > 0:
        print(f"\n[NOTE] {n_no_obs} stations had no usable SVH file/observations, skipped.")

    result = pd.concat(all_rows, ignore_index=True)
    OUT_COMBINED.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(OUT_COMBINED, index=False)

    print(f"\nDone. Combined observed+modelled table saved to: {OUT_COMBINED}")
    print(f"Total rows: {len(result)}")
    print(f"Stations included: {result['station_id'].nunique()}")

    ds.close()


if __name__ == "__main__":
    main()
