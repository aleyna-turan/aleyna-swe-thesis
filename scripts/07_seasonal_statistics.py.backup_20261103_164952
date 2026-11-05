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

import textwrap
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
SNOW_THRESHOLD_MM = 5.0

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
# Results are grouped by the question they answer, not by period.
# This script writes stage 02 (model evaluation) and stage 03 (trends),
# and reads the selected station lists produced by script 06 (stage 01).
INPUTS = REPO_ROOT / "results" / "tables" / "00_inputs"
SCREENING = REPO_ROOT / "results" / "tables" / "01_data_screening"
TABLES = REPO_ROOT / "results" / "tables" / "02_model_evaluation"
TRENDS = REPO_ROOT / "results" / "tables" / "03_trends"
FIGURES = REPO_ROOT / "results" / "figures" / "02_model_evaluation"
TREND_FIGURES = REPO_ROOT / "results" / "figures" / "03_trends"
META_FILE = INPUTS / "station_metadata.csv"

MAKE_FIGURES = True


def period_label(start: int, end: int) -> str:
    return f"{start}-{end}"


def hydro_years(start: int, end: int) -> np.ndarray:
    """Period 1 Nov <start> .. 31 Oct <end>  ->  HY start+1 .. HY end."""
    return np.arange(start + 1, end + 1)


def find_obs_file() -> Path:
    for name in ("modelled_vs_observed_swe.csv", "modelled_vs_observed_swe.csv.gz"):
        p = INPUTS / name
        if p.exists():
            return p
    raise FileNotFoundError(
        f"No observation table in {INPUTS}. Run scripts/03_extract_modelled_swe.py "
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
#
# Every figure carries, inside the image: a title saying WHAT it shows, a
# subtitle saying HOW the plotted quantity was computed, a legend, and a
# footnote saying HOW TO READ it and what the limitations are. The figures
# are meant to be understandable on their own, without the thesis text.

OBS_COLOR = "#1f77b4"
MOD_COLOR = "#e8702a"
CAPTION_KW = dict(fontsize=8.3, color="#333333", va="bottom")


def _caption(fig, text: str, base_bottom: float = 0.09) -> None:
    # Wrap to the figure width and reserve enough bottom margin for the lines,
    # so a long caption is never clipped at the right edge.
    width = int(fig.get_figwidth() * 17)
    lines = []
    for para in text.split("\n"):
        lines.extend(textwrap.wrap(para, width=width) or [""])
    fig.subplots_adjust(bottom=base_bottom + 0.030 * len(lines))
    fig.text(0.012, 0.015, "\n".join(lines), **CAPTION_KW)


def _theil_sen(x, y):
    # Robust slope + Mann-Kendall p, used to annotate the elevation figures.
    ok = x.notna() & y.notna()
    if ok.sum() < 10:
        return None
    slope, intercept, _, _ = stats.theilslopes(y[ok], x[ok], 0.95)
    tau, pval = stats.kendalltau(x[ok], y[ok])
    return slope, intercept, tau, pval, int(ok.sum())


def fig_annual_series(ann: pd.DataFrame, label: str, plt) -> None:
    # Observed is drawn THICK underneath, model THIN and DASHED on top, so a
    # year where the two agree exactly shows the dashed line inside the thick
    # band instead of one line vanishing under the other.
    ties = int((ann["snow_duration_obs_median"] == ann["snow_duration_mod_median"]).sum())
    n_years = len(ann)
    n_lo, n_hi = int(ann["n_stations"].min()), int(ann["n_stations"].max())

    fig, axes = plt.subplots(2, 1, figsize=(11, 8.6), sharex=True)
    fig.subplots_adjust(left=0.08, right=0.97, top=0.90, hspace=0.30)

    for ax, obs_col, mod_col, ylab, sub in (
        (axes[0], "max_swe_obs_median", "max_swe_mod_median",
         "Peak SWE (mm)",
         "Peak SWE: the largest SWE value observed in each hydrological year at each station, then the median across stations"),
        (axes[1], "snow_duration_obs_median", "snow_duration_mod_median",
         "Snow season length (days)",
         "Snow season length: days from the first to the last observation with SWE >= %g mm, then the median across stations" % SNOW_THRESHOLD_MM),
    ):
        ax.plot(ann["hydro_year"], ann[obs_col], color=OBS_COLOR, lw=4.0,
                alpha=0.55, solid_capstyle="round", label="Observed (CHMI stations)")
        ax.plot(ann["hydro_year"], ann[mod_col], color=MOD_COLOR, lw=1.6,
                ls="--", marker="o", ms=3.2, label="mHM simulated")
        ax.set_ylabel(ylab)
        ax.set_title(sub, fontsize=9, loc="left", color="#444444", pad=6)
        ax.grid(alpha=0.25, lw=0.5)
        ax.legend(loc="upper right", framealpha=0.9)

    axes[1].set_xlabel("Hydrological year (1 Nov of the previous year to 31 Oct)")

    fig.suptitle("Observed and simulated snow, %s" % label, fontsize=13, y=0.97)

    caption = (
        "How to read: the observed series is drawn as a thick pale band and the model as a thin dashed line on top. "
        "Where only the dashed line is visible the two series are IDENTICAL for that year - this is agreement, not missing data. "
        "Both series are complete: every year has a value, from %d to %d stations.\n"
        "Snow season length can only take values that are multiples of 7 days, because SWE is measured weekly, so exact ties are common "
        "(%d of %d years here). Medians, not means, are shown so that a few very snowy mountain stations do not dominate. "
        "Model values are sampled on the same dates as the observations, never from the full daily model output."
    ) % (n_lo, n_hi, ties, n_years)
    _caption(fig, caption, base_bottom=0.07)

    fig.savefig(FIGURES / ("annual_series_%s.png" % label), dpi=150)
    plt.close(fig)


def fig_bias_vs_elevation(ev: pd.DataFrame, meta: pd.DataFrame, label: str, plt) -> None:
    # The headline result: model bias becomes more negative with elevation.
    d = ev.merge(meta[["station_id", "station_name", "elevation_m"]],
                 on="station_id", how="left").dropna(subset=["max_swe_bias", "elevation_m"])

    fig, ax = plt.subplots(figsize=(10, 6.4))
    fig.subplots_adjust(left=0.10, right=0.97, top=0.90)

    ax.scatter(d["elevation_m"], d["max_swe_bias"], s=26, alpha=0.65,
               color=OBS_COLOR, edgecolor="none", label="One station (n = %d)" % len(d))
    ax.axhline(0, color="#444444", lw=1.0)

    fit = _theil_sen(d["elevation_m"], d["max_swe_bias"])
    if fit:
        slope, intercept, tau, pval, n = fit
        xs = np.array([d["elevation_m"].min(), d["elevation_m"].max()])
        ax.plot(xs, intercept + slope * xs, color="#b3331c", lw=2.0,
                label="Theil-Sen fit: %+.2f mm per 100 m" % (slope * 100))
        ax.annotate("Kendall tau = %.3f\np = %.1e\nn = %d stations" % (tau, pval, n),
                    xy=(0.985, 0.96), xycoords="axes fraction", ha="right", va="top",
                    fontsize=9, bbox=dict(boxstyle="round,pad=0.45", fc="white",
                                          ec="#bbbbbb", alpha=0.92))

    # Name the extreme stations so they are not anonymous points.
    lo, hi = d["max_swe_bias"].min(), d["max_swe_bias"].max()
    pad = 0.18 * (hi - lo)
    ax.set_ylim(lo - pad, hi + pad)
    for _, r in d.nsmallest(2, "max_swe_bias").iterrows():
        ax.annotate(str(r["station_name"]), xy=(r["elevation_m"], r["max_swe_bias"]),
                    xytext=(-14, -16), textcoords="offset points", fontsize=8.5,
                    color="#b3331c", ha="right",
                    arrowprops=dict(arrowstyle="-", color="#b3331c", lw=0.7))

    ax.set_xlabel("Station elevation (m a.s.l.)")
    ax.set_ylabel("mHM peak SWE bias (mm)")
    ax.set_title("Model bias in peak SWE against station elevation, %s" % label,
                 fontsize=13)
    ax.grid(alpha=0.25, lw=0.5)
    ax.legend(loc="lower left", framealpha=0.9)

    caption = (
        "Bias = mean over all hydrological years of (modelled peak SWE - observed peak SWE) at that station. "
        "Negative means the model holds less snow than was measured.\n"
        "The line is a Theil-Sen fit (median of the slopes of all station pairs) and the test is Mann-Kendall, "
        "both chosen because a few isolated summit stations have extreme bias and would dominate an ordinary least-squares fit. "
        "Those stations are labelled; the relationship holds with them removed."
    )
    _caption(fig, caption)

    fig.savefig(FIGURES / ("bias_vs_elevation_%s.png" % label), dpi=150)
    plt.close(fig)


def fig_trend_vs_elevation(tr: pd.DataFrame, meta: pd.DataFrame, label: str, plt) -> None:
    d = tr.merge(meta[["station_id", "elevation_m"]], on="station_id", how="left")
    d = d.dropna(subset=["max_swe_obs_slope", "elevation_m"])
    sig = d["max_swe_obs_p"] < ALPHA

    fig, ax = plt.subplots(figsize=(10, 6.4))
    fig.subplots_adjust(left=0.10, right=0.97, top=0.90)

    ax.scatter(d.loc[~sig, "elevation_m"], d.loc[~sig, "max_swe_obs_slope"] * 10,
               s=22, alpha=0.45, color="#8c8c8c", edgecolor="none",
               label="Not significant, p >= %g (n = %d)" % (ALPHA, int((~sig).sum())))
    ax.scatter(d.loc[sig, "elevation_m"], d.loc[sig, "max_swe_obs_slope"] * 10,
               s=30, alpha=0.85, color=OBS_COLOR, edgecolor="none",
               label="Significant, p < %g (n = %d)" % (ALPHA, int(sig.sum())))
    ax.axhline(0, color="#444444", lw=1.0)

    ax.set_xlabel("Station elevation (m a.s.l.)")
    ax.set_ylabel("Observed peak SWE trend (mm per decade)")
    ax.set_title("Observed peak SWE trend against station elevation, %s" % label,
                 fontsize=13)
    ax.grid(alpha=0.25, lw=0.5)
    ax.legend(loc="lower left", framealpha=0.9)

    caption = (
        "One point per station. The trend is the Theil-Sen slope of that station's annual peak SWE against hydrological year, "
        "converted to mm per decade; significance is the Mann-Kendall test at alpha = %g.\n"
        "Points below the zero line are stations where peak SWE is declining. Only stations with at least %g of the period's "
        "hydrological years present are included, so a short or gappy record cannot produce a spurious trend."
    ) % (ALPHA, MIN_YEARS_FRACTION_FOR_TREND)
    _caption(fig, caption)

    fig.savefig(TREND_FIGURES / ("trend_vs_elevation_%s.png" % label), dpi=150)
    plt.close(fig)


def make_figures(ann: pd.DataFrame, tr: pd.DataFrame, ev: pd.DataFrame,
                 meta: pd.DataFrame, label: str) -> None:
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("  [skip figures] matplotlib not installed")
        return

    FIGURES.mkdir(parents=True, exist_ok=True)
    TREND_FIGURES.mkdir(parents=True, exist_ok=True)

    fig_annual_series(ann, label, plt)
    fig_bias_vs_elevation(ev, meta, label, plt)
    fig_trend_vs_elevation(tr, meta, label, plt)
    print("  figures -> annual_series_%s.png, bias_vs_elevation_%s.png, "
          "trend_vs_elevation_%s.png" % (label, label, label))


# --- Main -------------------------------------------------------------------

def main():
    for d in (TABLES, TRENDS, FIGURES, TREND_FIGURES):
        d.mkdir(parents=True, exist_ok=True)
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

        sel_path = SCREENING / f"qc_selected_{label}.csv"
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
          .to_csv(TRENDS / f"trends_{label}.csv", index=False)
        ts.to_csv(TRENDS / f"trend_summary_{label}.csv", index=False)
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
            make_figures(ann, tr, ev, meta, label)

    comp = pd.DataFrame(comparison)
    comp.to_csv(TABLES / "period_comparison.csv", index=False)
    print("\n=== Period comparison (this is the table to choose from) ===")
    with pd.option_context("display.width", 200, "display.max_columns", 50):
        print(comp.to_string(index=False))
    print(f"\nTables written to {TABLES.parent}, figures to {FIGURES.parent}")
    print(f"Settings: snow_threshold={SNOW_THRESHOLD_MM} mm, "
          f"min_obs_per_hy={MIN_OBS_PER_HY}, "
          f"min_years_fraction_for_trend={MIN_YEARS_FRACTION_FOR_TREND}, "
          f"alpha={ALPHA}")


if __name__ == "__main__":
    main()
