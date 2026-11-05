"""
compare_duration_definitions.py

Standalone. Reads the pipeline's inputs, writes its own report, figure and
table, and CHANGES NOTHING in scripts 06-08 or their outputs. Nothing has to
be decided before running it.

The question
------------
snow_duration is the number of days from the first observation with
SWE >= 5 mm to the last one, within a hydrological year.

If a station records no observation at or above the threshold all year, there
is no first and no last. Two ways to record that:

  (A) MISSING  - what the pipeline does now. The station-year is dropped from
      every median and from the station's trend fit.

  (B) ZERO     - the snow season lasted zero days. That is a measurement, not
      a gap: the station reported all winter and never reached 5 mm.

Why it is not cosmetic
----------------------
Under (A) the median in a mild year is taken only over stations that still had
a season, so it answers "when there was a snow season, how long was it?" -
a number that can stay flat while snow seasons become rare. Under (B) it
answers "how long was the snow season, counting zero when there was none?",
which captures seasons getting shorter AND seasons disappearing.

The same applies to the trends: Theil-Sen fits only on years with a value, so
under (A) a station that loses its snow season contributes nothing for those
years rather than contributing a zero.

Note the asymmetry this script assumes for (B): counts and amounts take zero,
dates stay missing. A snow-free year has duration 0, but onset, melt-out and
peak timing remain undefined, because there is no date when there is no snow.
mean_swe, n_snow_weeks and swe_days already take zero in the pipeline, so (B)
would make snow_duration consistent with them.

Outputs
-------
results/DURATION_DEFINITION_COMPARISON.txt
results/tables/02_model_evaluation/snow_season_fraction_<period>.csv
results/figures/02_model_evaluation/duration_definition_comparison_<period>.png

Run:
    cd /Users/aleynaturan/Desktop/a/aleyna-swe-thesis
    source .venv/bin/activate
    python scripts/compare_duration_definitions.py
"""

from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

# --- Settings, matching scripts 06-08 -------------------------------------

PERIOD = "1976-2025"
HY_FIRST, HY_LAST = 1977, 2025

SNOW_THRESHOLD_MM = 5.0
FREQUENCY_HOMOGENISATION = "dominant"
MIN_OBS_PER_HY = 10
MIN_YEARS_FRACTION_FOR_TREND = 0.80
ALPHA = 0.05

REPO_ROOT = Path(__file__).resolve().parent.parent
TABLES = REPO_ROOT / "results" / "tables"
INPUTS = TABLES / "00_inputs"
SCREENING = TABLES / "01_data_screening"
EVAL = TABLES / "02_model_evaluation"
FIGURES = REPO_ROOT / "results" / "figures" / "02_model_evaluation"
REPORT = REPO_ROOT / "results" / "DURATION_DEFINITION_COMPARISON.txt"

_lines = []


def say(t=""):
    print(t)
    _lines.append(t)


def head(t):
    say()
    say("=" * 78)
    say(t)
    say("=" * 78)
    say()


# --- Load ------------------------------------------------------------------

def find_obs_file() -> Path:
    for n in ("modelled_vs_observed_swe.csv", "modelled_vs_observed_swe.csv.gz"):
        p = INPUTS / n
        if p.exists():
            return p
    raise FileNotFoundError(f"No paired table in {INPUTS}")


def load() -> pd.DataFrame:
    df = pd.read_csv(find_obs_file(), parse_dates=["date"])
    df["hydro_year"] = np.where(df["date"].dt.month >= 11,
                                df["date"].dt.year + 1, df["date"].dt.year)
    hy_start = pd.to_datetime(dict(year=df["hydro_year"] - 1, month=11, day=1))
    df["dohy"] = (df["date"] - hy_start).dt.days + 1
    df = df.dropna(subset=["modelled_swe_mm"])

    # Same weekday homogenisation as the pipeline: CHMI moved the weekly
    # measurement day from Friday to Monday during 1970, and some station-eras
    # reported more often than weekly. Keeping each year's modal weekday makes
    # every station-year comparable.
    if FREQUENCY_HOMOGENISATION == "dominant":
        dow = df["date"].dt.dayofweek
        modal = (df.assign(dow=dow).groupby("hydro_year")["dow"]
                   .agg(lambda x: x.value_counts().idxmax()))
        df = df[dow == df["hydro_year"].map(modal)]

    sel = set(pd.read_csv(SCREENING / f"qc_selected_{PERIOD}.csv")["station_id"])
    return df[df["station_id"].isin(sel)
              & df["hydro_year"].between(HY_FIRST, HY_LAST)], len(sel)


# --- Duration under both definitions ---------------------------------------

def durations(df: pd.DataFrame) -> pd.DataFrame:
    """One row per station-year, duration under both definitions."""
    rows = []
    for (sid, hy), g in df.groupby(["station_id", "hydro_year"], sort=False):
        if len(g) < MIN_OBS_PER_HY:
            continue
        g = g.sort_values("dohy")
        row = {"station_id": sid, "hydro_year": int(hy)}
        for src, col in (("obs", "observed_swe_mm"), ("mod", "modelled_swe_mm")):
            snow = g[col].to_numpy() >= SNOW_THRESHOLD_MM
            if snow.any():
                d = g["dohy"].to_numpy()[snow]
                v = float(d.max() - d.min())
                row[f"{src}_nan"] = v
                row[f"{src}_zero"] = v
                row[f"{src}_hassnow"] = 1
            else:
                row[f"{src}_nan"] = np.nan     # definition A
                row[f"{src}_zero"] = 0.0       # definition B
                row[f"{src}_hassnow"] = 0
        rows.append(row)
    return pd.DataFrame(rows)


def theil_sen(years, values):
    """Robust slope plus Mann-Kendall p.

    Theil-Sen is the MEDIAN of the slopes of all pairs of points, so a few
    extreme years cannot drag the fit. Mann-Kendall (Kendall's tau) tests only
    whether the ranking relationship exists, assuming nothing about the shape
    of the distribution - which matters here because duration is bounded below
    at zero and strongly skewed in mild years.
    """
    ok = ~np.isnan(values)
    if ok.sum() < 10:
        return np.nan, np.nan, int(ok.sum())
    slope, _, _, _ = stats.theilslopes(values[ok], years[ok], 0.95)
    _, p = stats.kendalltau(years[ok], values[ok])
    return float(slope), float(p), int(ok.sum())


def trends(dur: pd.DataFrame, col: str, n_years: int) -> pd.DataFrame:
    min_years = MIN_YEARS_FRACTION_FOR_TREND * n_years
    rows = []
    for sid, g in dur.groupby("station_id"):
        g = g.sort_values("hydro_year")
        v = g[col].to_numpy(dtype=float)
        if (~np.isnan(v)).sum() < min_years:
            continue          # station does not qualify under this definition
        slope, p, n = theil_sen(g["hydro_year"].to_numpy(dtype=float), v)
        rows.append({"station_id": sid, "slope": slope, "p": p, "n_years": n})
    return pd.DataFrame(rows)


# --- Main ------------------------------------------------------------------

def main():
    say("SNOW DURATION: MISSING versus ZERO FOR A SNOW-FREE YEAR")
    say(f"Period {PERIOD} (hydrological years {HY_FIRST}-{HY_LAST}), "
        f"threshold {SNOW_THRESHOLD_MM:g} mm")
    say()
    say("Definition A (current) : a snow-free year has NO duration value.")
    say("Definition B (proposed): a snow-free year has duration 0.")

    df, n_sel = load()
    dur = durations(df)
    n_years = HY_LAST - HY_FIRST + 1
    say()
    say(f"{n_sel} selected stations, {len(dur):,} station-years with at least "
        f"{MIN_OBS_PER_HY} observations.")

    # ---- 1. how often is a year snow-free -------------------------------
    head("1. HOW OFTEN IS A STATION-YEAR SNOW-FREE?")
    say("This is the quantity that definition A discards.")
    say()
    for src, name in (("obs", "Observed"), ("mod", "mHM")):
        n0 = int((dur[f"{src}_hassnow"] == 0).sum())
        say(f"  {name:9s}: {n0:>6,} of {len(dur):,} station-years "
            f"({100 * n0 / len(dur):.1f}%) had no snow season")

    frac = (dur.groupby("hydro_year")
              .agg(n_stations=("station_id", "size"),
                   obs_with_snow=("obs_hassnow", "sum"),
                   mod_with_snow=("mod_hassnow", "sum"))
              .reset_index())
    frac["obs_pct"] = 100 * frac["obs_with_snow"] / frac["n_stations"]
    frac["mod_pct"] = 100 * frac["mod_with_snow"] / frac["n_stations"]

    say()
    say("  Years where the fewest stations had a snow season (observed):")
    say()
    say(f"  {'hydro year':>11} {'stations':>9} {'with snow':>10} {'%':>7}")
    for _, r in frac.nsmallest(6, "obs_pct").iterrows():
        say(f"  {int(r['hydro_year']):>11} {int(r['n_stations']):>9} "
            f"{int(r['obs_with_snow']):>10} {r['obs_pct']:>6.1f}%")

    sl, p, n = theil_sen(frac["hydro_year"].to_numpy(dtype=float),
                         frac["obs_pct"].to_numpy(dtype=float))
    say()
    say(f"  Trend in the share of stations with a snow season: "
        f"{sl * 10:+.2f} percentage points per decade "
        f"(Mann-Kendall p = {p:.2e}, n = {n} years)")
    say("  This series is new information either way: it separates 'seasons")
    say("  are getting shorter' from 'seasons are happening less often'.")

    # ---- 2. the median in mild years --------------------------------------
    head("2. THE NETWORK MEDIAN, SIDE BY SIDE")
    say("Under A the median is taken only over stations that had a season.")
    say("Under B every selected station counts, with zero where there was none.")
    say()
    med = (dur.groupby("hydro_year")
             .agg(n=("station_id", "size"),
                  n_with_snow=("obs_hassnow", "sum"),
                  A=("obs_nan", "median"),
                  B=("obs_zero", "median"))
             .reset_index())
    say(f"  {'hydro year':>11} {'stations':>9} {'with snow':>10} "
        f"{'A (days)':>10} {'B (days)':>10} {'diff':>8}")
    for _, r in med.nsmallest(8, "n_with_snow").iterrows():
        a, b = r["A"], r["B"]
        a_txt = "n/a" if pd.isna(a) else f"{a:.0f}"
        d_txt = "" if pd.isna(a) else f"{b - a:+.0f}"
        say(f"  {int(r['hydro_year']):>11} {int(r['n']):>9} "
            f"{int(r['n_with_snow']):>10} {a_txt:>10} {b:>10.0f} {d_txt:>8}")
    say()
    say(f"  Mean over all {len(med)} years:  A = {med['A'].mean():.1f} days, "
        f"B = {med['B'].mean():.1f} days")

    # ---- 3. the trends ----------------------------------------------------
    head("3. THE TREND, SIDE BY SIDE")
    say("Per station: Theil-Sen slope of duration against hydrological year,")
    say("Mann-Kendall significance at alpha = %.2f. A station enters the fit" % ALPHA)
    say(f"only with at least {MIN_YEARS_FRACTION_FOR_TREND:.0%} of the period's "
        f"{n_years} years present - which itself differs between the two")
    say("definitions, because under A a snow-free year is absent.")
    say()
    say(f"  {'':24} {'A (missing)':>14} {'B (zero)':>14}")
    out = {}
    for src, name in (("obs", "Observed"), ("mod", "mHM")):
        row = {}
        for tag in ("nan", "zero"):
            t = trends(dur, f"{src}_{tag}", n_years)
            sig = t["p"] < ALPHA
            row[tag] = dict(
                n_stations=len(t),
                median_slope=float(t["slope"].median() * 10) if len(t) else np.nan,
                n_sig_neg=int((sig & (t["slope"] < 0)).sum()),
                pct_sig_neg=100 * float((sig & (t["slope"] < 0)).mean()) if len(t) else np.nan,
            )
        out[src] = row
        say(f"  {name} - stations in fit  "
            f"{row['nan']['n_stations']:>14} {row['zero']['n_stations']:>14}")
        say(f"  {name} - median slope/dec "
            f"{row['nan']['median_slope']:>14.2f} {row['zero']['median_slope']:>14.2f}")
        say(f"  {name} - sig. declining   "
            f"{row['nan']['n_sig_neg']:>14} {row['zero']['n_sig_neg']:>14}")
        say(f"  {name} - % sig. declining "
            f"{row['nan']['pct_sig_neg']:>13.1f}% {row['zero']['pct_sig_neg']:>13.1f}%")
        say()

    # ---- 4. what it means --------------------------------------------------
    head("4. READING THIS")
    o = out["obs"]
    d_slope = o["zero"]["median_slope"] - o["nan"]["median_slope"]
    say("Observed median trend moves from %.2f to %.2f days per decade (%+.2f)"
        % (o["nan"]["median_slope"], o["zero"]["median_slope"], d_slope))
    say("when snow-free years are counted as zero rather than dropped.")
    say()
    say("Definition A answers: when a snow season occurred, how long was it?")
    say("Definition B answers: how long was the snow season, counting zero")
    say("                      when there was none?")
    say()
    say("Both are legitimate. B is the one that registers a station losing its")
    say("snow season entirely, which under A simply leaves the sample. For a")
    say("thesis about long-term snow change that disappearance is arguably the")
    say("change of interest, and it is currently recorded as missing data.")
    say()
    say("Whichever is chosen must be stated in the methods. Under A the number")
    say("of contributing stations varies year to year, and a reader will")
    say("otherwise assume it is constant.")

    # ---- outputs -----------------------------------------------------------
    EVAL.mkdir(parents=True, exist_ok=True)
    frac.to_csv(EVAL / f"snow_season_fraction_{PERIOD}.csv", index=False)
    make_figure(med, frac, dur, n_years)

    REPORT.write_text("\n".join(_lines) + "\n")
    print(f"\nReport  : {REPORT}")
    print(f"Table   : {EVAL / f'snow_season_fraction_{PERIOD}.csv'}")


def make_figure(med, frac, dur, n_years):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("[skip figure] matplotlib not installed")
        return
    import textwrap

    FIGURES.mkdir(parents=True, exist_ok=True)
    A_COLOR, B_COLOR = "#8c8c8c", "#1f77b4"

    fig, axes = plt.subplots(3, 1, figsize=(11, 11.4))
    fig.subplots_adjust(left=0.085, right=0.97, top=0.93, hspace=0.42)

    # panel 1 - the median under both definitions
    ax = axes[0]
    # B thick and pale underneath, A thin and dashed on top, so that years
    # where the two definitions agree show the dashed line INSIDE the band
    # rather than one series vanishing under the other. The two diverge only
    # in years when some station was snow-free.
    ax.plot(med["hydro_year"], med["B"], color=B_COLOR, lw=4.0, alpha=0.55,
            solid_capstyle="round", label="B: snow-free year counts as 0 days")
    ax.plot(med["hydro_year"], med["A"], color=A_COLOR, lw=1.6, ls="--",
            marker="o", ms=3.0, label="A: snow-free year dropped (current)")
    ax.set_ylabel("Median snow season (days)")
    ax.set_title("Network median snow season length under each definition, observed",
                 fontsize=10, loc="left")
    ax.grid(alpha=0.25, lw=0.5)
    ax.legend(fontsize=8.5, framealpha=0.9)

    # panel 2 - the share of stations that had a season
    ax = axes[1]
    ax.plot(frac["hydro_year"], frac["obs_pct"], color="#1f77b4", lw=2.2,
            label="Observed")
    ax.plot(frac["hydro_year"], frac["mod_pct"], color="#e8702a", lw=1.5,
            ls="--", marker="o", ms=3.0, label="mHM")
    ax.set_ylim(0, 102)
    ax.set_ylabel("Stations with a snow season (%)")
    ax.set_title("Share of stations recording any snow season - the information "
                 "definition A discards", fontsize=10, loc="left")
    ax.grid(alpha=0.25, lw=0.5)
    ax.legend(fontsize=8.5, framealpha=0.9)

    # panel 3 - per-station trend distributions
    ax = axes[2]
    ta = trends(dur, "obs_nan", n_years)["slope"] * 10
    tb = trends(dur, "obs_zero", n_years)["slope"] * 10
    lo = float(min(ta.min(), tb.min()))
    hi = float(max(ta.max(), tb.max()))
    bins = np.linspace(lo, hi, 36)
    ax.hist(ta, bins=bins, color=A_COLOR, alpha=0.55,
            label="A (n = %d stations)" % len(ta))
    ax.hist(tb, bins=bins, color=B_COLOR, alpha=0.55,
            label="B (n = %d stations)" % len(tb))
    ax.axvline(0, color="#444444", lw=1.0)
    ax.axvline(float(ta.median()), color=A_COLOR, lw=2.0, ls="--")
    ax.axvline(float(tb.median()), color=B_COLOR, lw=2.0, ls="--")
    ax.set_xlabel("Per-station trend in snow season length (days per decade)")
    ax.set_ylabel("Stations")
    ax.set_title("Distribution of per-station trends; dashed lines are the medians",
                 fontsize=10, loc="left")
    ax.grid(alpha=0.25, lw=0.5)
    ax.legend(fontsize=8.5, framealpha=0.9)

    fig.suptitle("Does a snow-free winter count as zero days or as missing data?",
                 fontsize=13, y=0.975)

    caption = (
        "snow_duration is the number of days from the first to the last observation with SWE >= %g mm. "
        "When a station never reaches that threshold in a hydrological year there is no first and no last, "
        "and the value must either be dropped (A) or recorded as zero (B). In the top panel B is the thick pale "
        "line and A the dashed one on top, so where only the dashed line shows the two definitions agree - that "
        "is a year in which every station had a snow season.\n"
        "Dropping it means the median answers 'when there was a snow season, how long was it?', and a station "
        "that loses its snow season leaves the sample instead of contributing a zero. Recording zero means the "
        "median answers 'how long was the snow season, counting zero when there was none?'.\n"
        "Trends are Theil-Sen slopes with Mann-Kendall significance, both robust to the skew that duration has "
        "in mild years. Under either definition a station needs %.0f%% of the period's years present to enter "
        "the fit, which is itself easier to satisfy under B."
    ) % (SNOW_THRESHOLD_MM, 100 * MIN_YEARS_FRACTION_FOR_TREND)

    width = int(fig.get_figwidth() * 15)
    lines = []
    for para in caption.split("\n"):
        lines.extend(textwrap.wrap(para, width=width) or [""])
    fig.subplots_adjust(bottom=0.055 + 0.021 * len(lines))
    fig.text(0.012, 0.012, "\n".join(lines), fontsize=8.3, color="#333333",
             va="bottom")

    out = FIGURES / f"duration_definition_comparison_{PERIOD}.png"
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"Figure  : {out}")


if __name__ == "__main__":
    main()
