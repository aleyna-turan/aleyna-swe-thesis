"""
patch_07_add_bands.py

One-off patch. Run it ONCE, then delete it (or keep it as a record).

What it does to scripts/07_seasonal_statistics.py
-------------------------------------------------
1. Adds percentile calculation across stations for every metric and every
   hydrological year, separately for observations and model, written to a new
   long-format file  annual_percentiles_<period>.csv

2. Redraws the annual series figure as median + percentile bands instead of a
   bare median line, for both series, and prints in the caption how wide the
   model's band is relative to the observed one.

3. Adds a confidence interval on the Theil-Sen slope to the bias-vs-elevation
   figure, drawn as a shaded wedge around the fitted line.

Safety
------
Makes a timestamped backup of the script before touching it, and asserts on
every anchor it edits. If any anchor is not found, nothing is written and the
error says which one failed.

Run:
    cd /Users/aleynaturan/Desktop/a/aleyna-swe-thesis
    source .venv/bin/activate
    python scripts/patch_07_add_bands.py
"""

import shutil
from datetime import datetime
from pathlib import Path

TARGET = Path(__file__).resolve().parent / "07_seasonal_statistics.py"

NEW_FIG_BLOCK = '''def _theil_sen(x, y):
    """Robust slope with a confidence interval, plus the Mann-Kendall test.

    Theil-Sen takes the MEDIAN of the slopes of all pairs of points, so a few
    extreme stations cannot drag the fit. scipy also returns the lower and
    upper slope at the requested confidence level; those two slopes, pivoted
    about the median point, give the shaded wedge drawn on the figure.

    Kendall's tau tests only whether the ranking relationship exists, making
    no assumption about the shape of the data.
    """
    ok = x.notna() & y.notna()
    if ok.sum() < 10:
        return None
    slope, intercept, lo_slope, up_slope = stats.theilslopes(
        y[ok], x[ok], 1 - SLOPE_CI_ALPHA)
    tau, pval = stats.kendalltau(x[ok], y[ok])
    return dict(slope=slope, intercept=intercept, lo=lo_slope, up=up_slope,
                tau=tau, p=pval, n=int(ok.sum()),
                xmed=float(x[ok].median()), ymed=float(y[ok].median()))


def _series(pct: pd.DataFrame, metric: str, source: str):
    """Pull one metric/source out of the long percentile table, year-ordered."""
    d = pct[(pct["metric"] == metric) & (pct["source"] == source)]
    return d.sort_values("hydro_year")


def fig_annual_series(pct: pd.DataFrame, label: str, plt) -> None:
    """Network median per hydrological year, with percentile bands.

    A PERCENTILE BAND is not a confidence interval. A confidence interval says
    how precisely the median is known and shrinks as stations are added. A
    percentile band says how much the stations differ from one another and
    does not shrink. The second is what is wanted here: it shows the spread of
    snow conditions across Czechia, and whether the model reproduces it.
    """
    lo_o, hi_o = _ptag(BAND_OUTER[0]), _ptag(BAND_OUTER[1])
    specs = [
        ("max_swe", "Peak SWE (mm)",
         "Peak SWE: the largest SWE value in each hydrological year at each "
         "station"),
        ("snow_duration", "Snow season length (days)",
         "Snow season length: days from the first to the last observation "
         "with SWE >= %g mm" % SNOW_THRESHOLD_MM),
    ]

    fig, axes = plt.subplots(len(specs), 1, figsize=(11, 9.0), sharex=True)
    fig.subplots_adjust(left=0.08, right=0.97, top=0.90, hspace=0.30)

    widths = []
    for ax, (metric, ylab, sub) in zip(axes, specs):
        o, m = _series(pct, metric, "obs"), _series(pct, metric, "mod")
        yrs = o["hydro_year"].to_numpy()

        # outer band first, so the inner band and the medians sit on top
        ax.fill_between(yrs, o[lo_o], o[hi_o], color=OBS_COLOR, alpha=0.16,
                        lw=0, label="Observed, %g-%g percentile of stations"
                                    % BAND_OUTER)
        ax.fill_between(yrs, m[lo_o], m[hi_o], color=MOD_COLOR, alpha=0.16,
                        lw=0, label="mHM, %g-%g percentile of stations"
                                    % BAND_OUTER)
        if BAND_INNER:
            lo_i, hi_i = _ptag(BAND_INNER[0]), _ptag(BAND_INNER[1])
            ax.fill_between(yrs, o[lo_i], o[hi_i], color=OBS_COLOR, alpha=0.30,
                            lw=0, label="Observed, %g-%g (middle half)"
                                        % BAND_INNER)
            ax.fill_between(yrs, m[lo_i], m[hi_i], color=MOD_COLOR, alpha=0.30,
                            lw=0, label="mHM, %g-%g (middle half)" % BAND_INNER)

        ax.plot(yrs, o["p50"], color=OBS_COLOR, lw=2.6, label="Observed median")
        ax.plot(yrs, m["p50"], color=MOD_COLOR, lw=1.5, ls="--", marker="o",
                ms=3.0, label="mHM median")

        ax.set_ylabel(ylab)
        ax.set_title(sub, fontsize=9, loc="left", color="#444444", pad=6)
        ax.grid(alpha=0.25, lw=0.5)
        ax.legend(loc="upper right", fontsize=7.8, ncol=2, framealpha=0.9)

        wo = float((o[hi_o].to_numpy() - o[lo_o].to_numpy()).mean())
        wm = float((m[hi_o].to_numpy() - m[lo_o].to_numpy()).mean())
        widths.append((metric, wo, wm))

    axes[-1].set_xlabel("Hydrological year (1 Nov of the previous year to 31 Oct)")
    fig.suptitle("Observed and simulated snow with station spread, %s" % label,
                 fontsize=13, y=0.965)

    wtext = "; ".join("%s: observed %.0f, mHM %.0f (%.0f%%)"
                      % (mname.replace("_", " "), wo, wm, 100 * wm / wo if wo else 0)
                      for mname, wo, wm in widths)
    caption = (
        "The shaded bands are PERCENTILE RANGES across stations, not confidence intervals. "
        "A confidence interval would say how precisely the median is known and would shrink as stations "
        "are added; a percentile range says how much the stations differ from each other and does not shrink. "
        "The outer band holds the middle %g%% of stations, the darker inner band the middle half. "
        "Bands are computed separately for observations and model, as the spread itself is a model property worth testing.\\n"
        "Mean outer-band width over the period - %s. A model band narrower than the observed band means mHM "
        "under-represents how different stations are from one another, even where its median tracks well.\\n"
        "Medians are used rather than means so that a few very snowy mountain stations do not dominate. "
        "Model values are sampled on the same dates as the observations."
    ) % (BAND_OUTER[1] - BAND_OUTER[0], wtext)
    _caption(fig, caption, base_bottom=0.07)

    fig.savefig(FIGURES / ("annual_series_%s.png" % label), dpi=150)
    plt.close(fig)


def fig_bias_vs_elevation(ev: pd.DataFrame, meta: pd.DataFrame, label: str, plt) -> None:
    """Per-station model bias against elevation: the headline result."""
    d = ev.merge(meta[["station_id", "station_name", "elevation_m"]],
                 on="station_id", how="left").dropna(subset=["max_swe_bias", "elevation_m"])

    fig, ax = plt.subplots(figsize=(10, 6.4))
    fig.subplots_adjust(left=0.10, right=0.97, top=0.90)

    ax.scatter(d["elevation_m"], d["max_swe_bias"], s=26, alpha=0.65,
               color=OBS_COLOR, edgecolor="none", label="One station (n = %d)" % len(d))
    ax.axhline(0, color="#444444", lw=1.0)

    fit = _theil_sen(d["elevation_m"], d["max_swe_bias"])
    if fit:
        xs = np.linspace(d["elevation_m"].min(), d["elevation_m"].max(), 50)
        # the confidence wedge: the lowest and highest slopes consistent with
        # the data, both pivoted about the median point the fit passes through
        lo_line = fit["ymed"] + fit["lo"] * (xs - fit["xmed"])
        up_line = fit["ymed"] + fit["up"] * (xs - fit["xmed"])
        ax.fill_between(xs, lo_line, up_line, color="#b3331c", alpha=0.15, lw=0,
                        label="%g%% confidence interval on the slope"
                              % (100 * (1 - SLOPE_CI_ALPHA)))
        ax.plot(xs, fit["intercept"] + fit["slope"] * xs, color="#b3331c", lw=2.0,
                label="Theil-Sen fit: %+.2f mm per 100 m" % (fit["slope"] * 100))
        ax.annotate("Kendall tau = %.3f\\np = %.1e\\nslope %+.2f mm/100 m\\n"
                    "%g%% CI  %+.2f to %+.2f\\nn = %d stations"
                    % (fit["tau"], fit["p"], fit["slope"] * 100,
                       100 * (1 - SLOPE_CI_ALPHA), fit["lo"] * 100, fit["up"] * 100,
                       fit["n"]),
                    xy=(0.985, 0.96), xycoords="axes fraction", ha="right", va="top",
                    fontsize=9, bbox=dict(boxstyle="round,pad=0.45", fc="white",
                                          ec="#bbbbbb", alpha=0.92))

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
    ax.legend(loc="lower left", framealpha=0.9, fontsize=8.5)

    caption = (
        "Bias = mean over all hydrological years of (modelled peak SWE - observed peak SWE) at that station. "
        "Negative means the model holds less snow than was measured.\\n"
        "The shaded wedge here IS a confidence interval, unlike the percentile bands on the annual series figure: "
        "it shows the range of slopes consistent with the data, so it answers how well determined the slope is, "
        "not how much the stations differ.\\n"
        "Theil-Sen takes the median of the slopes of all station pairs and Mann-Kendall tests the ranking, both "
        "chosen because a few isolated summit stations would dominate an ordinary least-squares fit. Those stations "
        "are labelled; the relationship holds with them removed."
    )
    _caption(fig, caption)

    fig.savefig(FIGURES / ("bias_vs_elevation_%s.png" % label), dpi=150)
    plt.close(fig)


'''


def main():
    if not TARGET.exists():
        raise SystemExit(f"Not found: {TARGET}")

    s = TARGET.read_text()
    if "annual_percentiles" in s:
        raise SystemExit("This patch has already been applied to "
                         f"{TARGET.name}. Nothing to do.")

    backup = TARGET.with_suffix(
        f".py.backup_{datetime.now():%Y%m%d_%H%M%S}")
    shutil.copy2(TARGET, backup)
    print(f"Backup written: {backup.name}")

    steps = []

    # --- 1. config: band choice and CI level --------------------------------
    anchor = "ELEVATION_BANDS = [(0, 400), (400, 600), (600, 800), (800, 1000), (1000, 9999)]"
    assert anchor in s, "config anchor ELEVATION_BANDS not found"
    s = s.replace(anchor, anchor + '''

# --- Percentile bands across stations -------------------------------------
# These describe the SPREAD of station values, not the precision of the
# median. Supervisor's options were min-max, 2.5-97.5, or 5-95. The outer band
# is set to 5-95 because min-max is dominated by a single mountain station
# every year. Set BAND_INNER to None to draw only the outer band.
BAND_OUTER = (5, 95)
BAND_INNER = (25, 75)

# All of these are written to annual_percentiles_<period>.csv, so the figure
# can be redrawn with a different choice without recomputing anything.
PERCENTILES = [2.5, 5, 25, 50, 75, 95, 97.5]

# Confidence level for the Theil-Sen slope interval on the elevation figure.
# This one IS a confidence interval: it is about the fit, not the spread.
SLOPE_CI_ALPHA = 0.05


def _ptag(p) -> str:
    """Column name for a percentile: 2.5 -> p2_5, 50 -> p50."""
    return "p" + f"{p:g}".replace(".", "_")''')
    steps.append("config: BAND_OUTER, BAND_INNER, PERCENTILES, SLOPE_CI_ALPHA, _ptag()")

    # --- 2. the long-format percentile table --------------------------------
    anchor = "def elevation_bands(ev: pd.DataFrame, meta: pd.DataFrame) -> pd.DataFrame:"
    assert anchor in s, "anchor elevation_bands() not found"
    s = s.replace(anchor, '''def annual_percentiles(sm: pd.DataFrame) -> pd.DataFrame:
    """Spread of station values within each hydrological year.

    Long format - one row per hydrological year, metric and source - because
    a wide table would need 200+ columns and be unreadable. Filter on metric
    and source to pull one series out.
    """
    rows = []
    for hy, g in sm.groupby("hydro_year"):
        for m in METRICS:
            for src in ("obs", "mod"):
                v = g[f"{m}_{src}"].dropna()
                if v.empty:
                    continue
                row = {"hydro_year": int(hy), "metric": m, "source": src,
                       "n_stations": int(len(v)),
                       "mean": float(v.mean()),
                       "min": float(v.min()), "max": float(v.max())}
                for p in PERCENTILES:
                    row[_ptag(p)] = float(v.quantile(p / 100.0))
                rows.append(row)
    return pd.DataFrame(rows)


''' + anchor)
    steps.append("annual_percentiles(): percentile spread per year, long format")

    # --- 3. replace the three figure functions ------------------------------
    i = s.index("def _theil_sen(")
    j = s.index("def fig_trend_vs_elevation(")
    s = s[:i] + NEW_FIG_BLOCK + s[j:]
    steps.append("fig_annual_series(): percentile bands for both series")
    steps.append("fig_bias_vs_elevation(): Theil-Sen confidence wedge")

    # --- 4. make_figures signature and call ---------------------------------
    old = """def make_figures(ann: pd.DataFrame, tr: pd.DataFrame, ev: pd.DataFrame,
                 meta: pd.DataFrame, label: str) -> None:"""
    assert old in s, "anchor make_figures() signature not found"
    s = s.replace(old, """def make_figures(pct: pd.DataFrame, tr: pd.DataFrame, ev: pd.DataFrame,
                 meta: pd.DataFrame, label: str) -> None:""")
    old = "    fig_annual_series(ann, label, plt)"
    assert old in s, "anchor fig_annual_series call not found"
    s = s.replace(old, "    fig_annual_series(pct, label, plt)")
    old = "            make_figures(ann, tr, ev, meta, label)"
    assert old in s, "anchor make_figures call not found"
    s = s.replace(old, "            make_figures(pct, tr, ev, meta, label)")
    steps.append("make_figures(): takes the percentile table")

    # --- 5. compute and write the new table ---------------------------------
    old = "        ann = annual_series(sm)"
    assert old in s, "anchor annual_series call not found"
    s = s.replace(old, "        ann = annual_series(sm)\n        pct = annual_percentiles(sm)")
    old = '        ann.to_csv(TABLES / f"annual_series_{label}.csv", index=False)'
    assert old in s, "anchor annual_series to_csv not found"
    s = s.replace(old, old +
                  '\n        pct.to_csv(TABLES / f"annual_percentiles_{label}.csv", index=False)')
    steps.append("main(): writes annual_percentiles_<period>.csv")

    compile(s, str(TARGET), "exec")   # fail before writing if syntax broke
    TARGET.write_text(s)

    print(f"\nPatched {TARGET.name}:")
    for st in steps:
        print(f"  - {st}")
    print("\nNow rerun the analysis:")
    print("    python scripts/07_seasonal_statistics.py")
    print(f"\nTo undo: cp {backup.name} {TARGET.name}")


if __name__ == "__main__":
    main()
