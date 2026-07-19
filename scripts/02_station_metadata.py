"""
02_station_metadata.py

Downloads CHMI SVH (snow water equivalent) station data from the open data
portal and builds a station metadata table:
  ID, name, latitude, longitude, elevation,
  observation start, observation end, number of missing values, total records

Data source: https://opendata.chmi.cz/meteorology/climate/historical_csv/
Free / open data, Creative Commons BY 4.0 license.

Usage:
    source .venv/bin/activate
    python scripts/02_station_metadata.py

Inputs:
    data/stations/svh_filenames.txt   (list of SVH filenames, one per line)
    -> if this doesn't exist yet, the script fetches the listing itself.

Outputs:
    data/stations/svh_raw/*.csv        (raw downloaded per-station files)
    data/stations/meta1.csv            (raw CHMI station metadata, cached)
    results/tables/station_metadata.csv  (the final table we asked for)
"""

import re
import sys
from pathlib import Path

import pandas as pd
import requests

# --- Config ---------------------------------------------------------------

BASE_URL = "https://opendata.chmi.cz/meteorology/climate/historical_csv"
SNOW_DIR_URL = f"{BASE_URL}/data/daily/snow/"
META1_URL = f"{BASE_URL}/metadata/meta1.csv"

# Period covered by our mHM simulations - used to pick the right metadata
# row for stations whose coordinates changed over time.
MHM_PERIOD_START = pd.Timestamp("1950-01-01", tz="UTC")
MHM_PERIOD_END = pd.Timestamp("2026-05-31", tz="UTC")

REPO_ROOT = Path(__file__).resolve().parent.parent
STATIONS_DIR = REPO_ROOT / "data" / "stations"
SVH_RAW_DIR = STATIONS_DIR / "svh_raw"
RESULTS_TABLE = REPO_ROOT / "results" / "tables" / "station_metadata.csv"

STATIONS_DIR.mkdir(parents=True, exist_ok=True)
SVH_RAW_DIR.mkdir(parents=True, exist_ok=True)
RESULTS_TABLE.parent.mkdir(parents=True, exist_ok=True)


# --- Step 1: get the list of SVH files -------------------------------------

def get_svh_filenames() -> list[str]:
    cache_file = STATIONS_DIR / "svh_filenames.txt"
    if cache_file.exists():
        print(f"Using cached filename list: {cache_file}")
        return cache_file.read_text().splitlines()

    print("Fetching SVH file listing from CHMI...")
    resp = requests.get(SNOW_DIR_URL, timeout=60)
    resp.raise_for_status()
    filenames = sorted(set(re.findall(r'dly-[^"]*SVH\.csv', resp.text)))
    cache_file.write_text("\n".join(filenames))
    print(f"Found {len(filenames)} SVH files. Saved list to {cache_file}")
    return filenames


# --- Step 2: download each SVH file (skip if already downloaded) -----------

def download_svh_files(filenames: list[str]) -> None:
    total = len(filenames)
    for i, fname in enumerate(filenames, start=1):
        out_path = SVH_RAW_DIR / fname
        if out_path.exists():
            continue  # already downloaded, skip
        url = SNOW_DIR_URL + fname
        try:
            resp = requests.get(url, timeout=60)
            resp.raise_for_status()
            out_path.write_bytes(resp.content)
        except requests.RequestException as e:
            print(f"  [WARNING] failed to download {fname}: {e}")
            continue
        if i % 50 == 0 or i == total:
            print(f"  downloaded {i}/{total}")


# --- Step 3: download / load station metadata (meta1.csv) ------------------

def load_station_metadata() -> pd.DataFrame:
    meta1_path = STATIONS_DIR / "meta1.csv"
    if not meta1_path.exists():
        print("Downloading meta1.csv (station metadata)...")
        resp = requests.get(META1_URL, timeout=60)
        resp.raise_for_status()
        meta1_path.write_bytes(resp.content)

    df = pd.read_csv(meta1_path)
    df.columns = [c.strip() for c in df.columns]
    df["BEGIN_DATE"] = pd.to_datetime(df["BEGIN_DATE"], utc=True)
    df["END_DATE"] = pd.to_datetime(df["END_DATE"], utc=True)
    return df


def pick_relevant_row(station_rows: pd.DataFrame) -> pd.Series:
    """A station can have several metadata rows (location changed over
    time). Pick the row whose validity period overlaps our mHM period the
    most; if several overlap, prefer the most recent one."""
    overlapping = station_rows[
        (station_rows["BEGIN_DATE"] <= MHM_PERIOD_END)
        & (station_rows["END_DATE"] >= MHM_PERIOD_START)
    ]
    candidates = overlapping if len(overlapping) > 0 else station_rows
    return candidates.sort_values("BEGIN_DATE").iloc[-1]


# --- Step 4: parse one station's SVH file -----------------------------------

def summarize_svh_file(path: Path) -> dict:
    df = pd.read_csv(path)
    df.columns = [c.strip() for c in df.columns]
    df["DT"] = pd.to_datetime(df["DT"], utc=True)

    # A record counts as "missing" if QUALITY == 4 (Missing) or VALUE is empty
    is_missing = df["QUALITY"].eq(4.0) | df["VALUE"].isna()
    has_value = ~is_missing

    if has_value.any():
        obs_start = df.loc[has_value, "DT"].min()
        obs_end = df.loc[has_value, "DT"].max()
    else:
        obs_start = pd.NaT
        obs_end = pd.NaT

    return {
        "total_records": len(df),
        "n_missing": int(is_missing.sum()),
        "n_valid": int(has_value.sum()),
        "obs_start": obs_start,
        "obs_end": obs_end,
    }


# --- Main pipeline ----------------------------------------------------------

def main():
    filenames = get_svh_filenames()
    print(f"\nDownloading {len(filenames)} SVH files (skips ones already present)...")
    download_svh_files(filenames)

    print("\nLoading station metadata (meta1.csv)...")
    meta1 = load_station_metadata()

    print("\nBuilding station metadata table...")
    rows = []
    for i, fname in enumerate(filenames, start=1):
        wsi = fname.removeprefix("dly-").removesuffix("-SVH.csv")

        station_rows = meta1[meta1["WSI"] == wsi]
        if station_rows.empty:
            print(f"  [WARNING] no metadata found for {wsi}, skipping")
            continue
        meta_row = pick_relevant_row(station_rows)

        svh_path = SVH_RAW_DIR / fname
        if not svh_path.exists():
            print(f"  [WARNING] SVH file missing on disk for {wsi}, skipping")
            continue
        summary = summarize_svh_file(svh_path)

        rows.append({
            "station_id": wsi,
            "station_name": meta_row["FULL_NAME"],
            "latitude": meta_row["GEOGR2"],
            "longitude": meta_row["GEOGR1"],
            "elevation_m": meta_row["ELEVATION"],
            "obs_start": summary["obs_start"],
            "obs_end": summary["obs_end"],
            "n_valid_records": summary["n_valid"],
            "n_missing_records": summary["n_missing"],
            "total_records": summary["total_records"],
        })

        if i % 100 == 0 or i == len(filenames):
            print(f"  processed {i}/{len(filenames)}")

    result = pd.DataFrame(rows)
    result.to_csv(RESULTS_TABLE, index=False)
    print(f"\nDone. Station metadata table saved to: {RESULTS_TABLE}")
    print(f"Total stations in table: {len(result)}")


if __name__ == "__main__":
    sys.exit(main())
