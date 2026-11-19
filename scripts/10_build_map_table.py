"""
10_build_map_table.py

Builds one small CSV holding everything a station map needs, so that the
Colab notebook only has to read a single file instead of three.

Read-only with respect to the analysis: it reads existing tables and writes
one new file. Nothing else is touched.

Output
------
results/tables/05_maps/map_stations.csv   (all 1,147 stations, ~100 KB)

Columns
-------
station_id, station_name, latitude, longitude, elevation_m
selected              True for the 235 stations that passed screening for
                      1976-2025, False for the rest. The map draws the
                      unselected ones in muted grey.
n_years               hydrological years contributing, selected stations only
mean_peak_swe_mm      station mean of annual peak SWE, observed
mean_snow_duration_d  station mean of annual snow season length, observed
bias_peak_swe_mm      modelled minus observed peak SWE, the headline metric
notable               short reason string for the handful of stations worth
                      labelling on the map, blank otherwise

A note on mean_snow_duration_d: under the current definition a snow-free
hydrological year has no duration value and is dropped, so this is the mean
season length across the years that HAD a season. If the definition changes
to count snow-free years as zero, rerun this script and the column updates.

Run:
    cd /Users/aleynaturan/Desktop/a/aleyna-swe-thesis
    source .venv/bin/activate
    python scripts/10_build_map_table.py
"""

from pathlib import Path

import pandas as pd

PERIOD = "1976-2025"

# How many stations to flag for labelling on the map.
N_LABEL_HIGHEST = 5      # highest elevation among the selected stations
N_LABEL_BIAS = 3         # largest absolute model bias

REPO_ROOT = Path(__file__).resolve().parent.parent
TABLES = REPO_ROOT / "results" / "tables"
META_FILE = TABLES / "00_inputs" / "station_metadata.csv"
SELECTED_FILE = TABLES / "01_data_screening" / f"qc_selected_{PERIOD}.csv"
EVAL_FILE = TABLES / "02_model_evaluation" / f"model_eval_{PERIOD}.csv"
OUT_DIR = TABLES / "05_maps"
OUT_FILE = OUT_DIR / "map_stations.csv"

META_COLS = ["station_id", "station_name", "latitude", "longitude", "elevation_m"]
EVAL_COLS = {
    "n_years": "n_years",
    "max_swe_obs_mean": "mean_peak_swe_mm",
    "snow_duration_obs_mean": "mean_snow_duration_d",
    "max_swe_bias": "bias_peak_swe_mm",
}


def need(path: Path) -> Path:
    if not path.exists():
        raise SystemExit(f"Missing: {path}\nRun the earlier scripts first.")
    return path


def main():
    meta = pd.read_csv(need(META_FILE))
    missing = [c for c in META_COLS if c not in meta.columns]
    if missing:
        raise SystemExit(f"{META_FILE.name} is missing columns: {missing}\n"
                         f"It has: {list(meta.columns)}")
    out = meta[META_COLS].copy()
    print(f"station metadata      : {len(out):>5} stations")

    sel = set(pd.read_csv(need(SELECTED_FILE))["station_id"])
    out["selected"] = out["station_id"].isin(sel)
    print(f"selected for {PERIOD} : {int(out['selected'].sum()):>5} stations")

    ev = pd.read_csv(need(EVAL_FILE))
    have = {src: dst for src, dst in EVAL_COLS.items() if src in ev.columns}
    absent = set(EVAL_COLS) - set(have)
    if absent:
        print(f"  [note] not found in {EVAL_FILE.name}, left blank: {sorted(absent)}")
    out = out.merge(ev[["station_id"] + list(have)].rename(columns=have),
                    on="station_id", how="left")

    # --- which stations to label -----------------------------------------
    # The map gets crowded, so only a handful are labelled: the highest
    # selected stations, and those with the largest model bias. Both groups
    # are the ones a reader is likely to ask about.
    out["notable"] = ""
    pool = out[out["selected"]]

    if "elevation_m" in pool:
        for sid in pool.nlargest(N_LABEL_HIGHEST, "elevation_m")["station_id"]:
            out.loc[out["station_id"] == sid, "notable"] = "highest"

    if "bias_peak_swe_mm" in pool.columns:
        extreme = (pool.dropna(subset=["bias_peak_swe_mm"])
                       .assign(abs_bias=lambda d: d["bias_peak_swe_mm"].abs())
                       .nlargest(N_LABEL_BIAS, "abs_bias"))
        for sid in extreme["station_id"]:
            cur = out.loc[out["station_id"] == sid, "notable"].iloc[0]
            out.loc[out["station_id"] == sid, "notable"] = (
                "highest, extreme bias" if cur else "extreme bias")

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    out.to_csv(OUT_FILE, index=False)

    # --- report ------------------------------------------------------------
    print(f"\nwritten: {OUT_FILE}")
    print(f"  {len(out)} rows, {out['selected'].sum()} selected")
    print(f"  longitude {out['longitude'].min():.3f} to {out['longitude'].max():.3f}")
    print(f"  latitude  {out['latitude'].min():.3f} to {out['latitude'].max():.3f}")
    print(f"  elevation {out['elevation_m'].min():.0f} to {out['elevation_m'].max():.0f} m")

    lab = out[out["notable"] != ""]
    if len(lab):
        print(f"\nflagged for labelling ({len(lab)}):")
        cols = [c for c in ["station_name", "elevation_m", "bias_peak_swe_mm",
                            "mean_snow_duration_d", "notable"] if c in lab.columns]
        print(lab.sort_values("elevation_m", ascending=False)[cols].to_string(index=False))

    for c in ("mean_peak_swe_mm", "mean_snow_duration_d", "bias_peak_swe_mm"):
        if c in out.columns:
            v = out[c].dropna()
            if len(v):
                print(f"\n{c}: {len(v)} values, "
                      f"{v.min():.1f} to {v.max():.1f}, median {v.median():.1f}")


if __name__ == "__main__":
    main()
