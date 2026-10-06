"""
08_duration_bias_diagnosis.py

Step 8: is mHM's snow-season-length overestimate (+13 to +25 days, script 07)
a real model behaviour, or an artefact of the 1 mm snow threshold?

Two independent tests:

PART A - threshold sensitivity
    Recompute onset, melt-out, duration and weeks-with-snow at several SWE
    thresholds. If the bias is a threshold artefact it shrinks towards zero as
    the threshold rises. If it is real model behaviour it persists.

PART B - mechanism, date by date
    Every observation date is one of four cases:
        agree_snow   obs >= T and mod >= T
        agree_bare   obs <  T and mod <  T
        false_snow   obs <  T and mod >= T   <- model invents snow
        missed_snow  obs >= T and mod <  T   <- model misses snow
    If the duration bias comes from the model carrying small SWE at the
    shoulders of the season, false_snow should (a) dominate missed_snow,
    (b) carry small modelled values, and (c) cluster in the accumulation and
    melt phases rather than mid-winter.

Run:
    source .venv/bin/activate
    python scripts/08_duration_bias_diagnosis.py

Inputs : results/tables/modelled_vs_observed_swe.csv[.gz]
         results/tables/station_metadata.csv
         results/tables/qc_selected_<period>.csv       (script 06)
Outputs: results/tables/duration_threshold_sensitivity.csv
         results/tables/duration_bias_by_elevation.csv
         results/tables/snow_contingency_<period>.csv
         results/tables/false_snow_by_month_<period>.csv
         results/figures/duration_threshold_<period>.png
         results/figures/false_snow_seasonality_<period>.png

Conventions follow script 07: Czech hydrological year (1 Nov - 31 Oct, named
after the year it ends in), day-of-hydrological-year with 1 = 1 Nov, and the
same dominant-weekday sampling homogenisation.
"""

from pathlib import Path

import numpy as np
import pandas as pd

# --- Config ---------------------------------------------------------------

PERIODS = [(1956, 2025), (1966, 2025), (1976, 2025), (1986, 2025)]

# The period used for the detailed mechanism analysis (Part B) and figures.
MAIN_PERIOD = (1976, 2025)

THRESHOLDS_MM = [0.5, 1.0, 2.0, 5.0, 10.0, 20.0]

FREQUENCY_HOMOGENISATION = "dominant"   # as in scripts 06 and 07
MIN_OBS_PER_HY = 10

ELEVATION_BANDS = [(0, 400), (400, 600), (600, 800), (800, 1000), (1000, 9999)]

REPO_ROOT = Path(__file__).resolve().parent.parent
TABLES = REPO_ROOT / "results" / "tables"
FIGURES = REPO_ROOT / "results" / "figures"
META_FILE = TABLES / "station_metadata.csv"

MAKE_FIGURES = True


def period_label(start: int, end: int) -> str:
    return f"{start}-{end}"


def hydro_years(start: int, end: int) -> tuple:
    return start + 1, end


def find_obs_file() -> Path:
    for name in ("modelled_vs_observed_swe.csv", "modelled_vs_observed_swe.csv.gz"):
        p = TABLES / name
        if p.exists():
            return p
    raise FileNotFoundError(f"No observation table in {TABLES}. Run script 03.")


def load_observations(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, parse_dates=["date"])
    df["hydro_year"] = np.where(df["date"].dt.month >= 11,
                                df["date"].dt.year + 1,
                                df["date"].dt.year)
    hy_start = pd.to_datetime(dict(year=df["hydro_year"] - 1, month=11, day=1))
    df["dohy"] = (df["date"] - hy_start).dt.days + 1
    df = df.dropna(subset=["modelled_swe_mm"])

    if FREQUENCY_HOMOGENISATION != "none":
        dow = df["date"].dt.dayofweek
        if FREQUENCY_HOMOGENISATION == "monday":
            df = df[dow == 0]
        elif FREQUENCY_HOMOGENISATION == "dominant":
            modal = (df.assign(dow=dow).groupby("hydro_year")["dow"]
                       .agg(lambda x: x.value_counts().idxmax()))
            df = df[dow == df["hydro_year"].map(modal)]
        else:
            raise ValueError(FREQUENCY_HOMOGENISATION)
    return df


# --- Part A: threshold sensitivity -----------------------------------------

def duration_metrics(df: pd.DataFrame, thr: float) -> pd.DataFrame:
    """Onset, melt-out, duration and weeks-with-snow per station-HY, for both
    series, at one threshold."""
    out = []
    for (sid, hy), g in df.groupby(["station_id", "hydro_year"], sort=False):
        if len(g) < MIN_OBS_PER_HY:
            continue
        g = g.sort_values("dohy")
        dohy = g["dohy"].to_numpy()
        row = {"station_id": sid, "hydro_year": int(hy)}
        for src, col in (("obs", "observed_swe_mm"), ("mod", "modelled_swe_mm")):
            snow = g[col].to_numpy() >= thr
            if snow.any():
                first, last = dohy[snow].min(), dohy[snow].max()
                row[f"onset_{src}"] = float(first)
                row[f"meltout_{src}"] = float(last)
                row[f"duration_{src}"] = float(last - first)
                row[f"weeks_{src}"] = int(snow.sum())
            else:
                row[f"onset_{src}"] = np.nan
                row[f"meltout_{src}"] = np.nan
                row[f"duration_{src}"] = np.nan
                row[f"weeks_{src}"] = 0
        out.append(row)
    return pd.DataFrame(out)


def threshold_sensitivity(df: pd.DataFrame, meta: pd.DataFrame) -> tuple:
    rows, elev_rows = [], []
    elev = meta.set_index("station_id")["elevation_m"]

    for start, end in PERIODS:
        label = period_label(start, end)
        sel_path = TABLES / f"qc_selected_{label}.csv"
        if not sel_path.exists():
            print(f"  [skip] {sel_path.name} missing - run script 06")
            continue
        sel = pd.read_csv(sel_path)["station_id"].tolist()
        hy0, hy1 = hydro_years(start, end)
        sub = df[df["station_id"].isin(sel) & df["hydro_year"].between(hy0, hy1)]

        for thr in THRESHOLDS_MM:
            dm = duration_metrics(sub, thr)
            ok = dm["duration_obs"].notna() & dm["duration_mod"].notna()
            d = dm[ok]
            rows.append({
                "period": label,
                "threshold_mm": thr,
                "n_station_years": int(len(d)),
                "duration_obs_d": round(float(d["duration_obs"].mean()), 1),
                "duration_mod_d": round(float(d["duration_mod"].mean()), 1),
                "duration_bias_d": round(float((d["duration_mod"] - d["duration_obs"]).mean()), 2),
                "onset_bias_d": round(float((d["onset_mod"] - d["onset_obs"]).mean()), 2),
                "meltout_bias_d": round(float((d["meltout_mod"] - d["meltout_obs"]).mean()), 2),
                "weeks_obs": round(float(d["weeks_obs"].mean()), 2),
                "weeks_mod": round(float(d["weeks_mod"].mean()), 2),
                "weeks_bias": round(float((d["weeks_mod"] - d["weeks_obs"]).mean()), 2),
                # snow-free years the model turns into snow years
                "hy_obs_snowfree_mod_snow": int(
                    (dm["duration_obs"].isna() & dm["duration_mod"].notna()).sum()),
                "hy_obs_snow_mod_snowfree": int(
                    (dm["duration_obs"].notna() & dm["duration_mod"].isna()).sum()),
            })

            if (start, end) == MAIN_PERIOD:
                dd = d.assign(elevation_m=d["station_id"].map(elev))
                dd = dd.assign(bias=dd["duration_mod"] - dd["duration_obs"])
                for lo, hi in ELEVATION_BANDS:
                    b = dd[(dd["elevation_m"] >= lo) & (dd["elevation_m"] < hi)]
                    if b.empty:
                        continue
                    elev_rows.append({
                        "period": label,
                        "threshold_mm": thr,
                        "elev_band": f"{lo}-{hi if hi < 9999 else '+'} m",
                        "n_stations": int(b["station_id"].nunique()),
                        "duration_obs_d": round(float(b["duration_obs"].mean()), 1),
                        "duration_bias_d": round(float(b["bias"].mean()), 2),
                    })

    return pd.DataFrame(rows), pd.DataFrame(elev_rows)


# --- Part B: date-by-date mechanism ----------------------------------------

def contingency(df: pd.DataFrame, thr: float) -> pd.DataFrame:
    o = df["observed_swe_mm"].to_numpy() >= thr
    m = df["modelled_swe_mm"].to_numpy() >= thr
    case = np.select(
        [o & m, ~o & ~m, ~o & m, o & ~m],
        ["agree_snow", "agree_bare", "false_snow", "missed_snow"],
        default="unclassified",
    )
    d = df.assign(case=case)

    rows = []
    n = len(d)
    for c in ("agree_snow", "agree_bare", "false_snow", "missed_snow"):
        s = d[d["case"] == c]
        rows.append({
            "threshold_mm": thr,
            "case": c,
            "n_dates": int(len(s)),
            "pct_of_all": round(100 * len(s) / n, 2),
            "median_obs_mm": round(float(s["observed_swe_mm"].median()), 2) if len(s) else np.nan,
            "median_mod_mm": round(float(s["modelled_swe_mm"].median()), 2) if len(s) else np.nan,
            "p90_mod_mm": round(float(s["modelled_swe_mm"].quantile(0.90)), 2) if len(s) else np.nan,
        })
    return pd.DataFrame(rows), d


def false_snow_by_month(d: pd.DataFrame) -> pd.DataFrame:
    d = d.assign(month=d["date"].dt.month)
    tot = d.groupby("month").size().rename("n_dates")
    fs = d[d["case"] == "false_snow"].groupby("month").size().rename("n_false_snow")
    ms = d[d["case"] == "missed_snow"].groupby("month").size().rename("n_missed_snow")
    out = pd.concat([tot, fs, ms], axis=1).fillna(0).reset_index()
    out["pct_false_snow"] = (100 * out["n_false_snow"] / out["n_dates"]).round(2)
    out["pct_missed_snow"] = (100 * out["n_missed_snow"] / out["n_dates"]).round(2)
    order = [11, 12, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
    out["hy_order"] = out["month"].map({m: i for i, m in enumerate(order)})
    return out.sort_values("hy_order").drop(columns="hy_order")


# --- Figures ----------------------------------------------------------------

def make_figures(sens: pd.DataFrame, fsm: pd.DataFrame, label: str) -> None:
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("  [skip figures] matplotlib not installed")
        return
    FIGURES.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(8, 5))
    for period, g in sens.groupby("period"):
        ax.plot(g["threshold_mm"], g["duration_bias_d"], marker="o", label=period)
    ax.axhline(0, color="k", lw=0.8)
    ax.set_xscale("log")
    ax.set_xticks(THRESHOLDS_MM)
    ax.set_xticklabels([str(t) for t in THRESHOLDS_MM])
    ax.set_xlabel("Snow threshold (mm SWE)")
    ax.set_ylabel("mHM snow duration bias (days)")
    ax.set_title("Does the snow-duration bias survive a higher threshold?")
    ax.legend()
    fig.tight_layout()
    fig.savefig(FIGURES / f"duration_threshold_{label}.png", dpi=150)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(9, 5))
    x = np.arange(len(fsm))
    ax.bar(x - 0.2, fsm["pct_false_snow"], width=0.4, label="false snow (model only)")
    ax.bar(x + 0.2, fsm["pct_missed_snow"], width=0.4, label="missed snow (obs only)")
    ax.set_xticks(x)
    ax.set_xticklabels(["Nov", "Dec", "Jan", "Feb", "Mar", "Apr",
                        "May", "Jun", "Jul", "Aug", "Sep", "Oct"])
    ax.set_ylabel("% of observation dates in month")
    ax.set_xlabel("Month (hydrological year order)")
    ax.set_title(f"Where in the season the model and observations disagree, {label}")
    ax.legend()
    fig.tight_layout()
    fig.savefig(FIGURES / f"false_snow_seasonality_{label}.png", dpi=150)
    plt.close(fig)
    print(f"  figures -> duration_threshold_{label}.png, "
          f"false_snow_seasonality_{label}.png")


# --- Main -------------------------------------------------------------------

def main():
    TABLES.mkdir(parents=True, exist_ok=True)
    obs_file = find_obs_file()
    print(f"Reading {obs_file.name} ...")
    df = load_observations(obs_file)
    meta = pd.read_csv(META_FILE, usecols=["station_id", "elevation_m"])
    print(f"  {len(df)} paired records after homogenisation")

    print("\n--- PART A: threshold sensitivity ---")
    sens, elev = threshold_sensitivity(df, meta)
    sens.to_csv(TABLES / "duration_threshold_sensitivity.csv", index=False)
    elev.to_csv(TABLES / "duration_bias_by_elevation.csv", index=False)
    with pd.option_context("display.width", 200, "display.max_columns", 30):
        print(sens[["period", "threshold_mm", "duration_obs_d", "duration_mod_d",
                    "duration_bias_d", "onset_bias_d", "meltout_bias_d",
                    "weeks_bias"]].to_string(index=False))

    label = period_label(*MAIN_PERIOD)
    print(f"\n--- PART B: mechanism, {label} ---")
    sel = pd.read_csv(TABLES / f"qc_selected_{label}.csv")["station_id"].tolist()
    hy0, hy1 = hydro_years(*MAIN_PERIOD)
    sub = df[df["station_id"].isin(sel) & df["hydro_year"].between(hy0, hy1)]

    cont_all = []
    for thr in THRESHOLDS_MM:
        c, d = contingency(sub, thr)
        cont_all.append(c)
        if thr == 1.0:
            fsm = false_snow_by_month(d)
    cont = pd.concat(cont_all, ignore_index=True)
    cont.to_csv(TABLES / f"snow_contingency_{label}.csv", index=False)
    fsm.to_csv(TABLES / f"false_snow_by_month_{label}.csv", index=False)

    print(cont[cont["threshold_mm"].isin([1.0, 5.0, 20.0])].to_string(index=False))
    print("\nDisagreement by month (threshold 1 mm):")
    print(fsm[["month", "n_dates", "pct_false_snow", "pct_missed_snow"]].to_string(index=False))

    if MAKE_FIGURES:
        make_figures(sens, fsm, label)

    print(f"\nTables written to {TABLES}")


if __name__ == "__main__":
    main()
