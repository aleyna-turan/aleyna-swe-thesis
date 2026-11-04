"""
patch_07_split_spread.py

One-off patch, run ONCE after patch_07_add_bands.py.

Why
---
Putting the percentile bands on the same axes as the medians made the annual
series figure unreadable: the between-station spread (lowland stations with
almost no snow, Lysa hora with 360 mm) is far larger than the year-to-year
variation of the median, so the bands pushed the y-axis from ~100 mm to 350 mm
and squashed the median lines into the bottom of the panel. Four overlapping
translucent fills also blurred into one colour.

That is not a plotting accident, it is a property of the data worth stating in
the thesis: the differences BETWEEN stations are much larger than the changes
BETWEEN years. Two questions, two figures.

What it does to scripts/07_seasonal_statistics.py
-------------------------------------------------
1. fig_annual_series()  -> back to medians only, obs and model overlaid.
   Answers: does the model track the network from year to year?

2. fig_annual_spread()  -> NEW. A 2x2 grid: peak SWE and snow duration down
   the rows, observations and mHM across the columns, shared y-axis per row so
   the two band shapes can be compared directly. Each panel is annotated with
   its own mean band width.
   Answers: does the model reproduce how much the stations differ?

3. Caption wrapping tightened in scripts 06 and 07 so long captions stop
   clipping at the right edge.

Safety
------
Timestamped backups before editing, assertions on every anchor. If an anchor
is missing nothing is written and the error names it.

Run:
    cd /Users/aleynaturan/Desktop/a/aleyna-swe-thesis
    source .venv/bin/activate
    python scripts/patch_07_split_spread.py
"""

import shutil
from datetime import datetime
from pathlib import Path

HERE = Path(__file__).resolve().parent
TARGET = HERE / "07_seasonal_statistics.py"
ALSO_WRAP = [HERE / "06_station_qc_selection.py"]

NEW_BLOCK = '''def fig_annual_series(pct: pd.DataFrame, label: str, plt) -> None:
    """Network median per hydrological year, medians only.

    Why no spread band here: the spread between stations is far larger than
    the variation of the median between years, so a band drawn on these axes
    would swamp the lines it is meant to annotate. Spread has its own figure,
    annual_spread_<period>.png.

    Observed is drawn THICK and pale underneath, the model THIN and dashed on
    top, so a year where the two agree exactly shows the dashed line inside
    the band rather than one series vanishing under the other.
    """
    specs = [
        ("max_swe", "Peak SWE (mm)",
         "Peak SWE: the largest SWE value in each hydrological year at each "
         "station, then the median across stations"),
        ("snow_duration", "Snow season length (days)",
         "Snow season length: days from the first to the last observation with "
         "SWE >= %g mm, then the median across stations" % SNOW_THRESHOLD_MM),
    ]

    fig, axes = plt.subplots(len(specs), 1, figsize=(11, 8.6), sharex=True)
    fig.subplots_adjust(left=0.08, right=0.97, top=0.90, hspace=0.30)

    ties = n_years = n_lo = n_hi = 0
    for ax, (metric, ylab, sub) in zip(axes, specs):
        o, m = _series(pct, metric, "obs"), _series(pct, metric, "mod")
        yrs = o["hydro_year"].to_numpy()
        ax.plot(yrs, o["p50"], color=OBS_COLOR, lw=4.0, alpha=0.55,
                solid_capstyle="round", label="Observed (CHMI stations)")
        ax.plot(yrs, m["p50"], color=MOD_COLOR, lw=1.6, ls="--", marker="o",
                ms=3.2, label="mHM simulated")
        ax.set_ylabel(ylab)
        ax.set_title(sub, fontsize=9, loc="left", color="#444444", pad=6)
        ax.grid(alpha=0.25, lw=0.5)
        ax.legend(loc="upper right", framealpha=0.9)
        if metric == "snow_duration":
            ties = int((o["p50"].to_numpy() == m["p50"].to_numpy()).sum())
            n_years = len(yrs)
            n_lo = int(o["n_stations"].min())
            n_hi = int(o["n_stations"].max())

    axes[-1].set_xlabel("Hydrological year (1 Nov of the previous year to 31 Oct)")
    fig.suptitle("Observed and simulated snow, %s" % label, fontsize=13, y=0.965)

    caption = (
        "How to read: the observed series is a thick pale band, the model a thin dashed line on top. "
        "Where only the dashed line shows, the two are IDENTICAL for that year - agreement, not missing data. "
        "Both series are complete: every year has a value, from %d to %d stations.\\n"
        "Snow season length can only take values that are multiples of 7 days, because SWE is measured weekly, "
        "so exact ties are common (%d of %d years here). Medians, not means, so that a few very snowy mountain "
        "stations do not dominate. Model values are sampled on the same dates as the observations.\\n"
        "Station-to-station spread is shown separately in annual_spread_%s.png."
    ) % (n_lo, n_hi, ties, n_years, label)
    _caption(fig, caption, base_bottom=0.07)

    fig.savefig(FIGURES / ("annual_series_%s.png" % label), dpi=150)
    plt.close(fig)


def fig_annual_spread(pct: pd.DataFrame, label: str, plt) -> None:
    """How much the stations differ, observations against model.

    PERCENTILE RANGE, NOT CONFIDENCE INTERVAL. A confidence interval says how
    precisely the median is known and narrows as stations are added. A
    percentile range says how far apart the stations are and does not narrow.
    The second is what is wanted: whether mHM reproduces the real variety of
    snow conditions across Czechia, not just the typical station.

    Observations and model are drawn in separate panels with a shared y-axis
    per row, because four overlapping translucent bands on one panel are
    impossible to read.
    """
    lo_o, hi_o = _ptag(BAND_OUTER[0]), _ptag(BAND_OUTER[1])
    lo_i, hi_i = _ptag(BAND_INNER[0]), _ptag(BAND_INNER[1])

    specs = [("max_swe", "Peak SWE (mm)"),
             ("snow_duration", "Snow season length (days)")]
    sources = [("obs", "Observed (CHMI stations)", OBS_COLOR),
               ("mod", "mHM simulated", MOD_COLOR)]

    fig, axes = plt.subplots(len(specs), 2, figsize=(12.5, 8.8), sharex=True)
    fig.subplots_adjust(left=0.075, right=0.98, top=0.90, hspace=0.30, wspace=0.14)

    widths = {}
    for row, (metric, ylab) in enumerate(specs):
        series = {src: _series(pct, metric, src) for src, _, _ in sources}
        ymax = max(float(series[s][hi_o].max()) for s in series) * 1.06

        for col, (src, src_label, color) in enumerate(sources):
            ax = axes[row, col]
            d = series[src]
            yrs = d["hydro_year"].to_numpy()

            ax.fill_between(yrs, d[lo_o], d[hi_o], color=color, alpha=0.20,
                            lw=0, label="%g-%g percentile of stations" % BAND_OUTER)
            ax.fill_between(yrs, d[lo_i], d[hi_i], color=color, alpha=0.42,
                            lw=0, label="%g-%g, the middle half" % BAND_INNER)
            ax.plot(yrs, d["p50"], color=color, lw=2.0, label="Median station")

            w = float((d[hi_o].to_numpy() - d[lo_o].to_numpy()).mean())
            widths[(metric, src)] = w
            ax.annotate("mean %g-%g width: %.0f"
                        % (BAND_OUTER[0], BAND_OUTER[1], w),
                        xy=(0.98, 0.95), xycoords="axes fraction",
                        ha="right", va="top", fontsize=8.5,
                        bbox=dict(boxstyle="round,pad=0.35", fc="white",
                                  ec="#cccccc", alpha=0.9))

            ax.set_ylim(0, ymax)
            ax.grid(alpha=0.25, lw=0.5)
            ax.set_title(src_label, fontsize=10, pad=6)
            if col == 0:
                ax.set_ylabel(ylab)
            else:
                ax.tick_params(labelleft=False)
            if row == 0:
                ax.legend(loc="upper left", fontsize=7.8, framealpha=0.9)

    for col in range(2):
        axes[-1, col].set_xlabel("Hydrological year (ends 31 Oct)")

    fig.suptitle("Spread of station values, observations against model, %s" % label,
                 fontsize=13, y=0.965)

    ratios = "; ".join(
        "%s observed %.0f, mHM %.0f (%.0f%%)"
        % (m.replace("_", " "), widths[(m, "obs")], widths[(m, "mod")],
           100 * widths[(m, "mod")] / widths[(m, "obs")] if widths[(m, "obs")] else 0)
        for m, _ in specs)

    caption = (
        "Each panel shows, for every hydrological year, the median station and the range the stations span. "
        "These are PERCENTILE RANGES, not confidence intervals: they describe how far apart the stations are, "
        "and do not narrow as stations are added.\\n"
        "Mean outer-band width - %s. A model band narrower than the observed band means mHM compresses the "
        "difference between snowy and bare stations, under-representing the variety of conditions even where "
        "its median tracks well.\\n"
        "The bands are wide because the network genuinely spans lowland stations with almost no snow and "
        "mountain stations with several hundred mm. That between-station spread is much larger than the "
        "year-to-year movement of the median, which is why the two are drawn as separate figures."
    ) % ratios
    _caption(fig, caption, base_bottom=0.06)

    fig.savefig(FIGURES / ("annual_spread_%s.png" % label), dpi=150)
    plt.close(fig)


'''


def main():
    if not TARGET.exists():
        raise SystemExit(f"Not found: {TARGET}")
    s = TARGET.read_text()

    if "fig_annual_spread" in s:
        raise SystemExit(f"Already applied to {TARGET.name}. Nothing to do.")
    if "annual_percentiles" not in s:
        raise SystemExit("Run patch_07_add_bands.py first — this patch builds "
                         "on it.")

    stamp = f"{datetime.now():%Y%m%d_%H%M%S}"
    backup = TARGET.with_suffix(f".py.backup_{stamp}")
    shutil.copy2(TARGET, backup)
    print(f"Backup written: {backup.name}")

    steps = []

    # --- 1. replace fig_annual_series, add fig_annual_spread ---------------
    i = s.index("def fig_annual_series(")
    j = s.index("def fig_bias_vs_elevation(")
    s = s[:i] + NEW_BLOCK + s[j:]
    steps.append("fig_annual_series(): medians only again, readable")
    steps.append("fig_annual_spread(): NEW 2x2 obs-vs-model band comparison")

    # --- 2. call the new figure --------------------------------------------
    old = "    fig_annual_series(pct, label, plt)"
    assert old in s, "anchor fig_annual_series call not found"
    s = s.replace(old, old + "\n    fig_annual_spread(pct, label, plt)")

    old = ('    print("  figures -> annual_series_%s.png, bias_vs_elevation_%s.png, "\n'
           '          "trend_vs_elevation_%s.png" % (label, label, label))')
    assert old in s, "anchor figure print statement not found"
    s = s.replace(old,
                  '    print("  figures -> annual_series_%s.png, annual_spread_%s.png,"\n'
                  '          " bias_vs_elevation_%s.png, trend_vs_elevation_%s.png"\n'
                  '          % (label, label, label, label))')
    steps.append("make_figures(): draws and reports the new figure")

    # --- 3. caption wrapping ------------------------------------------------
    assert "fig.get_figwidth() * 17" in s, "anchor caption wrap width not found in 07"
    s = s.replace("fig.get_figwidth() * 17", "fig.get_figwidth() * 15")
    steps.append("caption wrap tightened (07) so long lines stop clipping")

    compile(s, str(TARGET), "exec")
    TARGET.write_text(s)

    for extra in ALSO_WRAP:
        if extra.exists():
            t = extra.read_text()
            if "fig.get_figwidth() * 17" in t:
                shutil.copy2(extra, extra.with_suffix(f".py.backup_{stamp}"))
                extra.write_text(t.replace("fig.get_figwidth() * 17",
                                           "fig.get_figwidth() * 15"))
                steps.append(f"caption wrap tightened ({extra.name})")

    print(f"\nPatched {TARGET.name}:")
    for st in steps:
        print(f"  - {st}")
    print("\nNow rerun:")
    print("    python scripts/07_seasonal_statistics.py")
    print("    python scripts/06_station_qc_selection.py   # for the wrap fix")
    print(f"\nTo undo: cp {backup.name} {TARGET.name}")


if __name__ == "__main__":
    main()
