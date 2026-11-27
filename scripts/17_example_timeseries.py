"""
17_example_timeseries.py

Picks one station at random from each elevation band and draws the model
against the observations for it, three ways.

Why these three figures
-----------------------
Everything else in the thesis is an aggregate: medians, bands, regressions
across 235 stations. A reader who has never seen a single station cannot
judge any of it. These figures go first in the results chapter and show
what the underlying data actually looks like.

  1. annual peak SWE per hydrological year, observed and modelled, with
     both Theil-Sen trends. Shows year-to-year variability, whether the
     model follows good and bad winters, the bias, and the trend, in one
     picture.

  2. the mean seasonal cycle, averaged over all years in weekly bins.
     Shows whether the model builds snow up and melts it out at the right
     times, which an annual maximum cannot reveal.

  3. one real winter per station, the year closest to that station's
     median peak. An average hides how ragged a single season is.

Random selection
----------------
One station is drawn at random from each band, using the seed set below.
The same seed always gives the same three stations, and the choice is
written to a CSV so the thesis can state which stations they are and that
they were not hand-picked to flatter the model.

Stations are drawn from those with at least MIN_YEARS hydrological years
present, so an example figure is not built on a gappy record. Screening
already demanded 90 percent coverage, so this rules out almost nothing, but
it is stated here rather than left to chance.

Note on sampling: the modelled series was extracted only on dates that
carry a real observation, so both series have the same, roughly weekly,
sampling. Nothing is interpolated anywhere in these figures.

Outputs
-------
results/figures/02_model_evaluation/example_annual_peak_swe.png
results/figures/02_model_evaluation/example_seasonal_cycle.png
results/figures/02_model_evaluation/example_single_winter.png
results/tables/02_model_evaluation/example_stations.csv

Run:
    cd /Users/aleynaturan/Desktop/a/aleyna-swe-thesis
    source .venv/bin/activate
    python scripts/17_example_timeseries.py
"""

from pathlib import Path
import textwrap

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import theilslopes

# ============================== settings ==============================

SEED = 20261126           # change this to draw a different three stations
PERIOD = "1976-2025"
MIN_YEARS = 40            # of the 49 hydrological years in the period

BAND_EDGES = [400, 700]
BAND_NAMES = ["Lowland (<= 400 m)", "Upland (400-700 m)", "Mountain (> 700 m)"]

# Same two colours as the band figures, so observed and modelled mean the
# same thing in every figure of the thesis. ColorBrewer Dark2.
C_OBS, C_MOD = "#1b9e77", "#d95f02"

# ======================================================================

REPO_ROOT = Path(__file__).resolve().parent.parent
TAB = REPO_ROOT / "results" / "tables"
METRICS = TAB / "02_model_evaluation" / f"seasonal_metrics_{PERIOD}.csv"
DAILY = TAB / "00_inputs" / "modelled_vs_observed_swe.csv"
DAILY_GZ = TAB / "00_inputs" / "modelled_vs_observed_swe.csv.gz"
FIG_DIR = REPO_ROOT / "results" / "figures" / "02_model_evaluation"
OUT_PICK = TAB / "02_model_evaluation" / "example_stations.csv"

RULE = "-" * 68


def banner(t):
    print(f"\n{RULE}\n{t}\n{RULE}")


def caption(fig, text):
    width = int(fig.get_figwidth() * 15)
    wrapped = "\n".join(textwrap.fill(p, width) for p in text.strip().split("\n"))
    fig.text(0.01, -0.01, wrapped, ha="left", va="top", fontsize=7.5, color="#333333")


def band_of(elev):
    return int(np.searchsorted(BAND_EDGES, elev, side="left"))


def sen(x, y):
    """Theil-Sen slope per decade, or None if there is too little to fit."""
    ok = np.isfinite(x) & np.isfinite(y)
    if ok.sum() < 10:
        return None
    s, i, _, _ = theilslopes(np.asarray(y)[ok], np.asarray(x)[ok])
    return s, i


# ------------------------------------------------------------- selection

def choose(metrics):
    """One station per band, drawn at random with the fixed seed."""
    per = (metrics.groupby("station_id")
                  .agg(station_name=("station_name", "first"),
                       elevation_m=("elevation_m", "first"),
                       n_years=("hydro_year", "nunique"))
                  .reset_index())
    per["band"] = per["elevation_m"].apply(band_of)
    eligible = per[per["n_years"] >= MIN_YEARS]

    print(f"stations in the metrics table : {len(per)}")
    print(f"with at least {MIN_YEARS} years      : {len(eligible)}")

    rng = np.random.default_rng(SEED)
    picked = []
    for b in range(3):
        pool = eligible[eligible["band"] == b]
        if pool.empty:
            print(f"  [warning] no eligible station in {BAND_NAMES[b]}")
            continue
        row = pool.iloc[int(rng.integers(len(pool)))].copy()
        row["band_name"] = BAND_NAMES[b]
        row["n_eligible_in_band"] = len(pool)
        picked.append(row)
        print(f"  {BAND_NAMES[b]:<22} -> {row.station_name} "
              f"({row.elevation_m:.0f} m, {row.n_years} years, "
              f"drawn from {len(pool)})")
    return pd.DataFrame(picked)


# --------------------------------------------------------------- fig 1

def fig_annual(metrics, picks):
    fig, axes = plt.subplots(len(picks), 1, figsize=(10, 3.1 * len(picks)),
                             sharex=True)
    axes = np.atleast_1d(axes)

    for ax, (_, p) in zip(axes, picks.iterrows()):
        d = (metrics[metrics.station_id == p.station_id]
             .sort_values("hydro_year"))
        x = d["hydro_year"].values

        for col, colour, label in ((("max_swe_obs"), C_OBS, "Observed"),
                                   (("max_swe_mod"), C_MOD, "mHM")):
            y = d[col].values
            ax.plot(x, y, color=colour, lw=1.3, marker="o", ms=3.2,
                    alpha=.9, label=label, zorder=3)
            f = sen(x, y)
            if f:
                s, i = f
                ax.plot(x, i + s * x, color=colour, lw=1.6, ls="--",
                        alpha=.85, zorder=4)
                ax.plot([], [], " ",
                        label=f"   trend {s * 10:+.1f} mm/decade")

        ax.set_title(f"{p.station_name} — {p.elevation_m:.0f} m — "
                     f"{p.band_name.split(' (')[0]}", fontsize=11, loc="left")
        ax.set_ylabel("Peak SWE (mm)")
        ax.grid(alpha=.18, lw=.4, ls=":")
        # "best" lets matplotlib find the emptiest corner of each panel,
        # which differs between a lowland and a mountain station
        ax.legend(fontsize=8, ncol=2, framealpha=.92, loc="best")

    axes[-1].set_xlabel("Hydrological year")
    fig.suptitle("Annual peak snow water equivalent at three example stations, "
                 f"{PERIOD}", fontsize=12.5, y=0.998)
    caption(fig, f"""
One station drawn at random from each elevation band (seed {SEED}), from stations with at least {MIN_YEARS} of the period's hydrological years. Each point is that year's highest recorded snow water equivalent; the solid line joins them and the dashed line is a Theil-Sen trend. The hydrological year runs from 1 November to 31 October and is named after the year it ends in.
The modelled series was sampled on the same dates as the observations, so neither series is interpolated and the two are directly comparable.
""")
    fig.tight_layout(rect=[0, 0, 1, 0.985])
    return fig


# --------------------------------------------------------------- fig 2

def fig_cycle(daily, picks):
    fig, axes = plt.subplots(1, len(picks), figsize=(4.4 * len(picks), 4.2))
    axes = np.atleast_1d(axes)

    for ax, (_, p) in zip(axes, picks.iterrows()):
        d = daily[daily.station_id == p.station_id].copy()
        # weekly bins through the hydrological year: observations are about
        # weekly, so a per-day average would be spiky rather than informative
        d["week"] = (d["dohy"] - 1) // 7
        g = d.groupby("week").agg(obs=("observed_swe_mm", "mean"),
                                  mod=("modelled_swe_mm", "mean"),
                                  n=("observed_swe_mm", "size")).reset_index()
        g = g[g["n"] >= 5]
        xs = g["week"] * 7 + 4

        ax.fill_between(xs, 0, g["obs"], color=C_OBS, alpha=.18, zorder=2)
        ax.plot(xs, g["obs"], color=C_OBS, lw=2.0, label="Observed", zorder=3)
        ax.plot(xs, g["mod"], color=C_MOD, lw=1.8, ls="--", label="mHM", zorder=4)

        # month ticks along the hydrological year, 1 Nov = day 1
        ticks = [1, 31, 62, 92, 123, 154, 182, 213]
        names = ["Nov", "Dec", "Jan", "Feb", "Mar", "Apr", "May", "Jun"]
        ax.set_xticks(ticks)
        ax.set_xticklabels(names, fontsize=9)
        ax.set_xlim(0, 230)
        ax.set_title(f"{p.station_name}\n{p.elevation_m:.0f} m", fontsize=10.5)
        ax.grid(alpha=.18, lw=.4, ls=":")
        ax.legend(fontsize=8.5, framealpha=.92)

    axes[0].set_ylabel("Mean snow water equivalent (mm)")
    fig.suptitle(f"Mean seasonal cycle, averaged over {PERIOD}", fontsize=12.5)
    caption(fig, """
Average snow water equivalent through the winter, in seven-day bins, over every year in the period. Bins holding fewer than five measurements are left out. Shading marks the observed curve.
This is where accumulation and melt timing show up: a model that peaks too early, or melts out too soon, does so visibly here even when its annual maximum is right. Note that averaging over many years smooths each individual winter, so the curves are flatter than any single season.
""")
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    return fig


# --------------------------------------------------------------- fig 3

def fig_winter(daily, metrics, picks):
    fig, axes = plt.subplots(1, len(picks), figsize=(4.4 * len(picks), 4.2))
    axes = np.atleast_1d(axes)

    for ax, (_, p) in zip(axes, picks.iterrows()):
        m = metrics[(metrics.station_id == p.station_id)
                    & metrics.max_swe_obs.notna()]
        if m.empty:
            continue
        # the year closest to this station's median peak: a typical winter
        # rather than a record one
        target = m["max_swe_obs"].median()
        hy = int(m.iloc[(m["max_swe_obs"] - target).abs().argsort().iloc[0]]["hydro_year"])

        d = daily[(daily.station_id == p.station_id)
                  & (daily.hydro_year == hy)].sort_values("dohy")

        ax.plot(d["dohy"], d["observed_swe_mm"], color=C_OBS, lw=1.6,
                marker="o", ms=3.5, label="Observed", zorder=3)
        ax.plot(d["dohy"], d["modelled_swe_mm"], color=C_MOD, lw=1.5, ls="--",
                marker="s", ms=3.0, label="mHM", zorder=4)

        ticks = [1, 31, 62, 92, 123, 154, 182, 213]
        names = ["Nov", "Dec", "Jan", "Feb", "Mar", "Apr", "May", "Jun"]
        ax.set_xticks(ticks)
        ax.set_xticklabels(names, fontsize=9)
        ax.set_xlim(0, 230)
        ax.set_title(f"{p.station_name}\n{p.elevation_m:.0f} m — "
                     f"winter {hy - 1}/{str(hy)[2:]}", fontsize=10.5)
        ax.grid(alpha=.18, lw=.4, ls=":")
        ax.legend(fontsize=8.5, framealpha=.92)

    axes[0].set_ylabel("Snow water equivalent (mm)")
    fig.suptitle("One typical winter at each station", fontsize=12.5)
    caption(fig, """
A single hydrological year at each station: the year whose peak was closest to that station's own median, so it is a representative winter rather than an extreme one. Markers are the actual measurement dates, roughly weekly, and the modelled values are taken on those same dates.
An average season is smooth; a real one is not. This figure shows the melt and refreeze episodes within a winter that the mean seasonal cycle averages away.
""")
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    return fig


# ----------------------------------------------------------------- main

def main():
    if not METRICS.exists():
        raise SystemExit(f"Missing: {METRICS}")

    banner(f"CHOOSING STATIONS  (seed {SEED})")
    metrics = pd.read_csv(METRICS)
    picks = choose(metrics)
    if picks.empty:
        raise SystemExit("No stations could be chosen.")

    src = DAILY if DAILY.exists() else DAILY_GZ
    if not src.exists():
        raise SystemExit(f"Missing: {DAILY} (or the .gz version)")

    banner("READING THE PAIRED DAILY DATA")
    print(f"{src.name} ... (large file, this takes a moment)")
    daily = pd.read_csv(src, parse_dates=["date"])
    daily = daily[daily.station_id.isin(picks.station_id)].copy()
    print(f"rows for the three chosen stations: {len(daily):,}")

    # hydrological year and day within it: 1 November is day 1, so a date
    # near New Year does not wrap around to the start of the series
    y, mth = daily["date"].dt.year, daily["date"].dt.month
    daily["hydro_year"] = np.where(mth >= 11, y + 1, y)
    start = pd.to_datetime(dict(year=daily["hydro_year"] - 1, month=11, day=1))
    daily["dohy"] = (daily["date"] - start).dt.days + 1

    FIG_DIR.mkdir(parents=True, exist_ok=True)
    OUT_PICK.parent.mkdir(parents=True, exist_ok=True)

    banner("FIGURES")
    for name, fn in (("example_annual_peak_swe", lambda: fig_annual(metrics, picks)),
                     ("example_seasonal_cycle", lambda: fig_cycle(daily, picks)),
                     ("example_single_winter", lambda: fig_winter(daily, metrics, picks))):
        fig = fn()
        out = FIG_DIR / f"{name}.png"
        fig.savefig(out, dpi=200, bbox_inches="tight")
        plt.close(fig)
        print(f"  {out.relative_to(REPO_ROOT)}")

    keep = picks.copy()
    keep.insert(0, "seed", SEED)
    keep.insert(1, "min_years_required", MIN_YEARS)
    keep.to_csv(OUT_PICK, index=False)
    print(f"  {OUT_PICK.relative_to(REPO_ROOT)}")

    banner("FOR THE METHODS SECTION")
    print(f"Random seed {SEED}, one station per elevation band, drawn from")
    print(f"stations with at least {MIN_YEARS} hydrological years present:")
    for _, p in picks.iterrows():
        print(f"  {p.band_name:<22} {p.station_name} ({p.elevation_m:.0f} m), "
              f"one of {p.n_eligible_in_band} eligible")


if __name__ == "__main__":
    main()
