"""
07_seasonal_statistics.py

Step 7: seasonal SWE statistics, model evaluation and trends, computed
separately for EACH candidate study period so the periods can be compared
before one is chosen.

Run:
    source .venv/bin/activate
    python scripts/07_seasonal_statistics.py


HYDROLOGICAL YEAR - read this before using any output
-----------------------------------------------------
In Czechia the hydrological year runs 1 November - 31 October and is named
after the calendar year in which it ENDS (CHMI / Czech meteorological
dictionary convention):

    HY 2025 = 1 Nov 2024 ... 31 Oct 2025

Every table produced here is indexed by `hydro_year` in that convention.

The period labels keep the calendar span you specified:

    "1956-2025"  =  1 Nov 1956 ... 31 Oct 2025  =  HY 1957 ... HY 2025  (69 HY)
    "1966-2025"  =  1 Nov 1966 ... 31 Oct 2025  =  HY 1967 ... HY 2025  (59 HY)
    "1976-2025"  =  1 Nov 1976 ... 31 Oct 2025  =  HY 1977 ... HY 2025  (49 HY)
    "1986-2025"  =  1 Nov 1986 ... 31 Oct 2025  =  HY 1987 ... HY 2025  (39 HY)

These are the same windows script 06 screened (it labelled them by season
start year; hydro_year = season + 1), so the selected station sets carry
over unchanged.

Statistics are computed over the FULL hydrological year, not a Nov-Apr
window. The data contain 1239 non-zero SWE observations outside Nov-Apr
(961 in October, 263 in May, up to 686 mm) at 430 stations, so a Nov-Apr
cut would truncate snow onset and melt-out - worst at exactly the
high-elevation stations. Nov-Apr is used only for availability screening
in script 06.

Day-of-hydrological-year (DOHY): 1 = 1 November. Used for all timing
metrics so that dates around New Year do not wrap.


PAIRED COMPARISON - an important limitation
-------------------------------------------
Script 03 extracted modelled SWE only on dates where a station has a real
observation. So modelled metrics here are computed from the SAME weekly
sampling as the observations, not from the full daily model output. This is
the fair comparison (like with like) but it means:
  - peak timing has ~7-day resolution and the true modelled peak may fall
    between two sampling dates;
  - season length is a weekly-resolution estimate, reported both as a
    date span and as a count of weeks with snow.
If you later want metrics from the full daily model series, script 03 has
to be rerun without the date-matching step.


Inputs
------
results/tables/modelled_vs_observed_swe.csv[.gz]   (script 03)
results/tables/station_metadata.csv                (script 02)
results/tables/qc_selected_<period>.csv            (script 06)

Outputs (results/tables/ and results/figures/)
----------------------------------------------
seasonal_metrics_<period>.csv    station x hydro_year, observed and modelled
annual_series_<period>.csv       network-median series per hydro_year
model_eval_<period>.csv          per-station bias/RMSE/correlation per metric
model_eval_summary_<period>.csv  network summary of the above
trends_<period>.csv              per-station Theil-Sen slope + Mann-Kendall p
trend_summary_<period>.csv       network trend summary, observed vs modelled
period_comparison.csv            one row per period - the table to choose from
elevation_bands_<period>.csv     metrics and bias by elevation band
figures/annual_maxswe_<period>.png
figures/trend_vs_elevation_<period>.png
"""

from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

# --- Config ---------------------------------------------------------------

# (period label start, period label end) as calendar years, matching script 06
PERIODS = [(1956, 2025), (1966, 2025), (1976, 2025), (1986, 2025)]

# Sampling-frequency homogenisation, IDENTICAL to script 06. CHMI moved the
# weekly measurement day from Friday to Monday during 1970, and some
# station-eras reported more often than weekly. Without this, stations and
# decades with denser sampling are over-weighted in every mean below.
# "dominant" keeps, for each hydrological year, only that year's modal
# weekday across the network. Options: "dominant" | "monday" | "none".
FREQUENCY_HOMOGENISATION = "dominant"

# SWE at or above this counts as "snow present" (mm water equivalent).
# 0.0 would let trace/rounding noise define season length.
SNOW_THRESHOLD_MM = 1.0

# A hydrological year is usable for a station only if it has at least this
# many observations (weekly sampling -> a full HY holds ~52, a full snow
# season ~26).
MIN_OBS_PER_HY = 10

# A station enters the trend analysis only with at least this fraction of
# the period's hydrological years present.
MIN_YEARS_FRACTION_FOR_TREND = 0.80

ALPHA = 0.05  # significance level for Mann-Kendall

ELEVATION_BANDS = [(0, 400), (400, 600), (600, 800), (800, 1000), (1000, 9999)]

# Metrics carried through model evaluation and trend analysis
METRICS = [
    "mean_swe",        # mean over the whole hydrological year, zeros included
    "mean_swe_djfm",   # mean over 1 Dec - 31 Mar only, a fixed comparable window
    "mean_swe_snow",   # mean over observations with snow only
    "max_swe",         # peak SWE
    "dohy_max",        # timing of peak, days since 1 Nov
    "dohy_first_snow", # snow onset
    "dohy_last_snow",  # melt-out
    "snow_duration",   # last_snow - first_snow, days
    "n_snow_weeks",    # observations with snow (~weeks with snow)
    "swe_days",        # sum(SWE) * 7, approximate mm-days of storage
]

REPO_ROOT = Path(__file__).resolve().parent.parent
TABLES = REPO_ROOT / "results" / "tables"
FIGURES = REPO_ROOT / "results" / "figures"
META_FILE = TABLES / "station_metadata.csv"

MAKE_FIGURES = True


def period_label(start: int, end: int) -> str:
    return f"{start}-{end}"


def hydro_years(start: int, end: int) -> np.ndarray:
    """Period 1 Nov <start> .. 31 Oct <end>  ->  HY start+1 .. HY end."""
    return np.arange(start + 1, end + 1)


def find_obs_file() -> Path:
    for name in ("modelled_vs_observed_swe.csv", "modelled_vs_observed_swe.csv.gz"):
        p = TABLES / name
        if p.exists():
            return p
    raise FileNotFoundError(
        f"No observation table in {TABLES}. Run scripts/03_extract_modelled_swe.py "
        f"or pull the gzipped copy from the repo."
    )


# --- Step 1: load and add hydrological-year columns -------------------------

def load_observations(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, parse_dates=["date"])

    # Czech hydrological year: Nov and Dec belong to the NEXT year's HY.
    df["hydro_year"] = np.where(df["date"].dt.month >= 11,
                                df["date"].dt.year + 1,
                                df["date"].dt.year)

    # Day of hydrological year, 1 = 1 November.
    hy_start = pd.to_datetime(dict(year=df["hydro_year"] - 1, month=11, day=1))
    df["dohy"] = (df["date"] - hy_start).dt.days + 1

    df = df.dropna(subset=["modelled_swe_mm"])

    if FREQUENCY_HOMOGENISATION != "none":
        n_all = len(df)
        dow = df["date"].dt.dayofweek
        if FREQUENCY_HOMOGENISATION == "monday":
            df = df[dow == 0]
        elif FREQUENCY_HOMOGENISATION == "dominant":
            modal = (df.assign(dow=dow)
                       .groupby("hydro_year")["dow"]
                       .agg(lambda x: x.value_counts().idxmax()))
            df = df[dow == df["hydro_year"].map(modal)]
        else:
            raise ValueError(FREQUENCY_HOMOGENISATION)
        print(f"  homogenisation '{FREQUENCY_HOMOGENISATION}': kept {len(df)} "
              f"of {n_all} records ({100 * len(df) / n_all:.1f}%)")

    return df


# --- Step 2: per station x hydrological year metrics ------------------------

def metrics_for_group(swe: np.ndarray, dohy: np.ndarray,
                      in_djfm: np.ndarray) -> dict:
    """All seasonal metrics for one station-HY, from one SWE series.

    Timing metrics use the FULL hydrological year so October onset and May
    melt-out are not truncated. Magnitude means are given both over the
    whole year and over a fixed 1 Dec - 31 Mar window, which is the
    comparable one if sampling differs between stations."""
    snow = swe >= SNOW_THRESHOLD_MM
    out = {
        "mean_swe": float(np.mean(swe)),
        "mean_swe_djfm": float(np.mean(swe[in_djfm])) if in_djfm.any() else np.nan,
        "max_swe": float(np.max(swe)),
        "n_snow_weeks": int(snow.sum()),
        "swe_days": float(np.sum(swe) * 7.0),
    }
    if snow.any():
        out["mean_swe_snow"] = float(np.mean(swe[snow]))
        out["dohy_max"] = float(dohy[int(np.argmax(swe))])
        out["dohy_first_snow"] = float(dohy[snow].min())
        out["dohy_last_snow"] = float(dohy[snow].max())
        out["snow_duration"] = out["dohy_last_snow"] - out["dohy_first_snow"]
    else:
        # snow-free hydrological year: amounts are 0, timings undefined
        for k in ("mean_swe_snow", "dohy_max", "dohy_first_snow",
                  "dohy_last_snow", "snow_duration"):
            out[k] = np.nan
    return out


def seasonal_metrics(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (sid, hy), g in df.groupby(["station_id", "hydro_year"], sort=False):
        if len(g) < MIN_OBS_PER_HY:
            continue
        g = g.sort_values("dohy")
        dohy = g["dohy"].to_numpy()
        in_djfm = g["date"].dt.month.isin([12, 1, 2, 3]).to_numpy()
        row = {"station_id": sid, "hydro_year": int(hy), "n_obs": len(g)}
        for src, col in (("obs", "observed_swe_mm"), ("mod", "modelled_swe_mm")):
            for k, v in metrics_for_group(g[col].to_numpy(), dohy, in_djfm).items():
                row[f"{k}_{src}"] = v
        rows.append(row)
    return pd.DataFrame(rows)


# --- Step 3: model evaluation, per station, per metric ----------------------

def model_evaluation(sm: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for sid, g in sm.groupby("station_id"):
        row = {"station_id": sid, "n_years": len(g)}
        for m in METRICS:
            o, p = g[f"{m}_obs"], g[f"{m}_mod"]
            ok = o.notna() & p.notna()
            n = int(ok.sum())
            row[f"{m}_n"] = n
            if n >= 3:
                d = (p[ok] - o[ok])
                row[f"{m}_bias"] = float(d.mean())
                row[f"{m}_mae"] = float(d.abs().mean())
                row[f"{m}_rmse"] = float(np.sqrt((d ** 2).mean()))
                row[f"{m}_obs_mean"] = float(o[ok].mean())
                row[f"{m}_mod_mean"] = float(p[ok].mean())
                if o[ok].std() > 0 and p[ok].std() > 0:
                    row[f"{m}_r"] = float(np.corrcoef(o[ok], p[ok])[0, 1])
                else:
                    row[f"{m}_r"] = np.nan
            else:
                for suffix in ("bias", "mae", "rmse", "obs_mean", "mod_mean", "r"):
                    row[f"{m}_{suffix}"] = np.nan
        rows.append(row)
    return pd.DataFrame(rows)


# --- Step 4: trends ---------------------------------------------------------

def theil_sen_mk(years: np.ndarray, values: np.ndarray) -> dict:
    """Theil-Sen slope (robust, no normality assumption) + Mann-Kendall test
    via Kendall's tau. Both are standard for hydrological trend work and
    tolerate the non-normal, zero-inflated distributions here."""
    ok = ~np.isnan(values)
    if ok.sum() < 10:
        return {"slope": np.nan, "intercept": np.nan, "tau": np.nan,
                "p_value": np.nan, "n": int(ok.sum())}
    y, v = years[ok], values[ok]
    slope, intercept, _, _ = stats.theilslopes(v, y, 0.95)
    tau, p = stats.kendalltau(y, v)
    return {"slope": float(slope), "intercept": float(intercept),
            "tau": float(tau), "p_value": float(p), "n": int(ok.sum())}


def trends(sm: pd.DataFrame, hys: np.ndarray) -> pd.DataFrame:
    min_years = MIN_YEARS_FRACTION_FOR_TREND * len(hys)
    rows = []
    for sid, g in sm.groupby("station_id"):
        if len(g) < min_years:
            continue
        g = g.sort_values("hydro_year")
        yrs = g["hydro_year"].to_numpy(dtype=float)
        row = {"station_id": sid, "n_years": len(g)}
        for m in METRICS:
            for src in ("obs", "mod"):
                r = theil_sen_mk(yrs, g[f"{m}_{src}"].to_numpy(dtype=float))
                row[f"{m}_{src}_slope"] = r["slope"]      # units per year
                row[f"{m}_{src}_p"] = r["p_value"]
                row[f"{m}_{src}_tau"] = r["tau"]
                row[f"{m}_{src}_n"] = r["n"]
        rows.append(row)
    return pd.DataFrame(rows)


def trend_summary(tr: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for m in METRICS:
        row = {"metric": m, "n_stations": int(len(tr))}
        for src in ("obs", "mod"):
            s = tr[f"{m}_{src}_slope"]
            p = tr[f"{m}_{src}_p"]
            sig = p < ALPHA
            row[f"{src}_median_slope_per_decade"] = round(float(s.median() * 10), 3)
            row[f"{src}_n_sig"] = int(sig.sum())
            row[f"{src}_n_sig_neg"] = int((sig & (s < 0)).sum())
            row[f"{src}_n_sig_pos"] = int((sig & (s > 0)).sum())
            row[f"{src}_pct_sig_neg"] = round(100 * float((sig & (s < 0)).mean()), 1)
        # does the model reproduce the observed station-to-station trend pattern?
        ok = tr[f"{m}_obs_slope"].notna() & tr[f"{m}_mod_slope"].notna()
        row["slope_agreement_r"] = (
            round(float(np.corrcoef(tr.loc[ok, f"{m}_obs_slope"],
                                    tr.loc[ok, f"{m}_mod_slope"])[0, 1]), 3)
            if ok.sum() >= 3 else np.nan
        )
        rows.append(row)
    return pd.DataFrame(rows)


# --- Step 5: network-level annual series and elevation bands ----------------

def annual_series(sm: pd.DataFrame) -> pd.DataFrame:
    agg = {}
    for m in METRICS:
        agg[f"{m}_obs_median"] = (f"{m}_obs", "median")
        agg[f"{m}_mod_median"] = (f"{m}_mod", "median")
        agg[f"{m}_obs_mean"] = (f"{m}_obs", "mean")
        agg[f"{m}_mod_mean"] = (f"{m}_mod", "mean")
    out = sm.groupby("hydro_year").agg(n_stations=("station_id", "nunique"), **agg)
    return out.reset_index()


def elevation_bands(ev: pd.DataFrame, meta: pd.DataFrame) -> pd.DataFrame:
    d = ev.merge(meta[["station_id", "elevation_m"]], on="station_id", how="left")
    rows = []
    for lo, hi in ELEVATION_BANDS:
        sub = d[(d["elevation_m"] >= lo) & (d["elevation_m"] < hi)]
        if sub.empty:
            continue
        row = {"elev_band": f"{lo}-{hi if hi < 9999 else '+'} m",
               "n_stations": int(len(sub))}
        for m in ("max_swe", "mean_swe", "snow_duration", "dohy_max"):
            row[f"{m}_obs_mean"] = round(float(sub[f"{m}_obs_mean"].mean()), 2)
            row[f"{m}_bias"] = round(float(sub[f"{m}_bias"].mean()), 2)
            row[f"{m}_rmse"] = round(float(sub[f"{m}_rmse"].mean()), 2)
            row[f"{m}_r"] = round(float(sub[f"{m}_r"].mean()), 3)
        rows.append(row)
    return pd.DataFrame(rows)


# --- Step 6: figures --------------------------------------------------------

def make_figures(ann: pd.DataFrame, tr: pd.DataFrame, meta: pd.DataFrame,
                 label: str) -> None:
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("  [skip figures] matplotlib not installed")
        return

    FIGURES.mkdir(parents=True, exist_ok=True)

    fig, axes = plt.subplots(2, 1, figsize=(10, 7), sharex=True)
    axes[0].plot(ann["hydro_year"], ann["max_swe_obs_median"], label="observed", lw=1.6)
    axes[0].plot(ann["hydro_year"], ann["max_swe_mod_median"], label="mHM", lw=1.6)
    axes[0].set_ylabel("Network median peak SWE (mm)")
    axes[0].legend()
    axes[0].set_title(f"Peak SWE and snow duration, {label} (hydrological years)")
    axes[1].plot(ann["hydro_year"], ann["snow_duration_obs_median"], label="observed", lw=1.6)
    axes[1].plot(ann["hydro_year"], ann["snow_duration_mod_median"], label="mHM", lw=1.6)
    axes[1].set_ylabel("Network median snow duration (days)")
    axes[1].set_xlabel("Hydrological year (ends 31 Oct)")
    axes[1].legend()
    fig.tight_layout()
    fig.savefig(FIGURES / f"annual_maxswe_{label}.png", dpi=150)
    plt.close(fig)

    d = tr.merge(meta[["station_id", "elevation_m"]], on="station_id", how="left")
    fig, ax = plt.subplots(figsize=(8, 5))
    sig = d["max_swe_obs_p"] < ALPHA
    ax.scatter(d.loc[~sig, "elevation_m"], d.loc[~sig, "max_swe_obs_slope"] * 10,
               s=16, alpha=0.5, label=f"p >= {ALPHA}")
    ax.scatter(d.loc[sig, "elevation_m"], d.loc[sig, "max_swe_obs_slope"] * 10,
               s=22, alpha=0.9, label=f"p < {ALPHA}")
    ax.axhline(0, lw=0.8, color="k")
    ax.set_xlabel("Station elevation (m)")
    ax.set_ylabel("Observed peak SWE trend (mm per decade)")
    ax.set_title(f"Peak SWE trend vs elevation, {label}")
    ax.legend()
    fig.tight_layout()
    fig.savefig(FIGURES / f"trend_vs_elevation_{label}.png", dpi=150)
    plt.close(fig)
    print(f"  figures -> annual_maxswe_{label}.png, trend_vs_elevation_{label}.png")


# --- Main -------------------------------------------------------------------

def main():
    TABLES.mkdir(parents=True, exist_ok=True)
    obs_file = find_obs_file()
    print(f"Reading {obs_file.name} ...")
    df = load_observations(obs_file)
    print(f"  {len(df)} paired records, {df['station_id'].nunique()} stations, "
          f"HY {df['hydro_year'].min()}-{df['hydro_year'].max()}")

    meta = pd.read_csv(META_FILE, usecols=["station_id", "station_name", "latitude",
                                           "longitude", "elevation_m"])

    comparison = []

    for start, end in PERIODS:
        label = period_label(start, end)
        hys = hydro_years(start, end)
        print(f"\n=== {label}  (HY {hys[0]}-{hys[-1]}, {len(hys)} hydrological years) ===")

        sel_path = TABLES / f"qc_selected_{label}.csv"
        if not sel_path.exists():
            print(f"  [skip] {sel_path.name} not found - run script 06 first")
            continue
        selected = pd.read_csv(sel_path)["station_id"].tolist()
        print(f"  selected stations: {len(selected)}")

        sub = df[df["station_id"].isin(selected)
                 & df["hydro_year"].between(hys[0], hys[-1])]

        sm = seasonal_metrics(sub)
        ev = model_evaluation(sm)
        tr = trends(sm, hys)
        ts = trend_summary(tr)
        ann = annual_series(sm)
        eb = elevation_bands(ev, meta)

        sm.merge(meta, on="station_id", how="left") \
          .to_csv(TABLES / f"seasonal_metrics_{label}.csv", index=False)
        ev.merge(meta, on="station_id", how="left") \
          .to_csv(TABLES / f"model_eval_{label}.csv", index=False)
        tr.merge(meta, on="station_id", how="left") \
          .to_csv(TABLES / f"trends_{label}.csv", index=False)
        ts.to_csv(TABLES / f"trend_summary_{label}.csv", index=False)
        ann.to_csv(TABLES / f"annual_series_{label}.csv", index=False)
        eb.to_csv(TABLES / f"elevation_bands_{label}.csv", index=False)

        ev_num = ev.select_dtypes(include=[np.number])
        comparison.append({
            "period": label,
            "hy_first": int(hys[0]), "hy_last": int(hys[-1]),
            "n_hydro_years": len(hys),
            "n_stations": len(selected),
            "n_station_years": int(len(sm)),
            "obs_mean_swe_mm": round(float(sm["mean_swe_obs"].mean()), 2),
            "obs_max_swe_mm": round(float(sm["max_swe_obs"].mean()), 2),
            "obs_dohy_max": round(float(sm["dohy_max_obs"].mean()), 1),
            "obs_snow_duration_d": round(float(sm["snow_duration_obs"].mean()), 1),
            "bias_max_swe_mm": round(float(ev_num["max_swe_bias"].mean()), 2),
            "rmse_max_swe_mm": round(float(ev_num["max_swe_rmse"].mean()), 2),
            "r_max_swe": round(float(ev_num["max_swe_r"].mean()), 3),
            "bias_dohy_max_d": round(float(ev_num["dohy_max_bias"].mean()), 2),
            "bias_duration_d": round(float(ev_num["snow_duration_bias"].mean()), 2),
            "n_trend_stations": int(len(tr)),
            "maxswe_trend_obs_dec": float(
                ts.loc[ts["metric"] == "max_swe", "obs_median_slope_per_decade"].iloc[0]),
            "maxswe_trend_mod_dec": float(
                ts.loc[ts["metric"] == "max_swe", "mod_median_slope_per_decade"].iloc[0]),
            "maxswe_pct_sig_neg_obs": float(
                ts.loc[ts["metric"] == "max_swe", "obs_pct_sig_neg"].iloc[0]),
            "duration_trend_obs_dec": float(
                ts.loc[ts["metric"] == "snow_duration", "obs_median_slope_per_decade"].iloc[0]),
        })

        print(f"  station-years: {len(sm)}   trend stations: {len(tr)}")
        print(f"  observed mean peak SWE: {sm['max_swe_obs'].mean():.1f} mm, "
              f"mean duration: {sm['snow_duration_obs'].mean():.0f} d")
        print(f"  mHM peak SWE bias: {ev_num['max_swe_bias'].mean():+.1f} mm, "
              f"r = {ev_num['max_swe_r'].mean():.2f}")

        if MAKE_FIGURES:
            make_figures(ann, tr, meta, label)

    comp = pd.DataFrame(comparison)
    comp.to_csv(TABLES / "period_comparison.csv", index=False)
    print("\n=== Period comparison (this is the table to choose from) ===")
    with pd.option_context("display.width", 200, "display.max_columns", 50):
        print(comp.to_string(index=False))
    print(f"\nAll tables written to {TABLES}")
    print(f"Settings: snow_threshold={SNOW_THRESHOLD_MM} mm, "
          f"min_obs_per_hy={MIN_OBS_PER_HY}, "
          f"min_years_fraction_for_trend={MIN_YEARS_FRACTION_FOR_TREND}, "
          f"alpha={ALPHA}")


if __name__ == "__main__":
    main()
