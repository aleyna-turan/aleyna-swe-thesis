"""
12_elevation_bands.py

Splits the stations into three elevation bands and produces:

  figures
    trend_vs_elevation_<metric>.png   per-station trend against elevation,
                                      coloured by band, with a Theil-Sen line
                                      fitted across all stations and the
                                      median level of each band marked
    band_trends_obs_vs_mod.png        four metrics, three bands, observed
                                      against modelled, as box plots

  table
    band_trends_<period>.csv          n, median trend, quartiles and the
                                      count of significant stations, for
                                      every band and metric, observed and
                                      modelled

Reads only results/tables/03_trends/trends_<period>.csv, which already
carries elevation_m and station_name. Writes figures and one table. It does
not touch any existing analysis output.

Why fixed band edges and not terciles
-------------------------------------
Terciles of this network cut at 395 m and 518 m, so the "high" band would
start at 518 m. That number describes where CHMI put its stations, not
where Czech highlands begin, and it moves whenever the network changes.
Fixed edges at 400 m and 700 m give 83 / 130 / 22 stations, the names mean
something, and the cuts stay put. Switch to terciles in one line below if
your supervisor prefers the balanced version.

Why Theil-Sen for the line across elevation
-------------------------------------------
Each dot is already a Theil-Sen trend over 49 years. Fitting a line across
elevation is a second regression, and the elevation axis is badly skewed:
213 of 235 stations sit below 700 m, and a single station sits at 1322 m.
Ordinary least squares minimises squared distance, so those few far-right
points would drag the line. Theil-Sen takes the median of the slopes of all
pairs of points, so it ignores up to about 29 percent of the data being
extreme before it breaks down. Same method, same reasoning, as everywhere
else in this pipeline.

Run:
    cd /Users/aleynaturan/Desktop/a/aleyna-swe-thesis
    source .venv/bin/activate
    python scripts/12_elevation_bands.py
"""

from pathlib import Path
import textwrap

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from scipy.stats import theilslopes, kendalltau

# ============================== settings ==============================

PERIOD = "1976-2025"
ALPHA = 0.05                 # significance level, as everywhere else

# The band edges, in metres. Everything at or below the first edge is the
# lowest band; everything above the last edge is the highest.
BAND_EDGES = [400, 700]
BAND_NAMES = ["Lowland (<= 400 m)", "Upland (400-700 m)", "Mountain (> 700 m)"]

# The balanced alternative, if you decide to switch after the meeting:
# BAND_EDGES = [395, 518]
# BAND_NAMES = ["Lower third", "Middle third", "Upper third"]

# ColorBrewer YlGnBu, three classes taken from the four-class ramp so the
# palest one is dropped and the points stay visible on white. ColorBrewer
# flags this scheme as colourblind safe.
BAND_COLOURS = ["#a1dab4", "#41b6c4", "#225ea8"]

# Observed against modelled, used in the box figure. These two carry the
# comparison the figure exists for, so they are different hues rather than
# two shades of one. ColorBrewer Dark2, which ColorBrewer flags as
# colourblind safe. The bands are named on the axis there, so colour does
# not have to carry them as well.
SRC_COLOURS = {"obs": "#1b9e77", "mod": "#d95f02"}

# The individual station dots sit on top of the boxes, so they use a darker
# shade of the same two hues: dark enough to read against the box fill, and
# still obviously the same colour family as the box they belong to.
SRC_DOT_COLOURS = {"obs": "#0b5d45", "mod": "#8a3d01"}

# Metrics to plot. key -> (nice name, unit after conversion to per decade)
METRICS = {
    "max_swe":         ("Peak SWE",          "mm per decade"),
    "snow_duration":   ("Snow season length", "days per decade"),
    "dohy_first_snow": ("Season onset",      "days per decade"),
    "dohy_last_snow":  ("Melt-out date",     "days per decade"),
}

# The stored per-station slopes are per year; the summary table reports per
# decade. The script verifies this against the summary rather than trusting
# it, and will tell you if the check fails.
SLOPE_TO_DECADE = 10.0

# ======================================================================

REPO_ROOT = Path(__file__).resolve().parent.parent
TRENDS_FILE = REPO_ROOT / "results" / "tables" / "03_trends" / f"trends_{PERIOD}.csv"
SUMMARY_FILE = REPO_ROOT / "results" / "tables" / "03_trends" / f"trend_summary_{PERIOD}.csv"
FIG_DIR = REPO_ROOT / "results" / "figures" / "03_trends"
TAB_DIR = REPO_ROOT / "results" / "tables" / "03_trends"
OUT_TABLE = TAB_DIR / f"band_trends_{PERIOD}.csv"


def caption(fig, text):
    """Put an explanatory note under the figure, wrapped to its width."""
    width = int(fig.get_figwidth() * 15)
    wrapped = "\n".join(textwrap.fill(p, width) for p in text.strip().split("\n"))
    fig.text(0.01, -0.01, wrapped, ha="left", va="top", fontsize=7.5,
             color="#333333")


def assign_band(elev):
    """Return the band index 0, 1 or 2 for an elevation in metres."""
    return int(np.searchsorted(BAND_EDGES, elev, side="left"))


def verify_units(df):
    """Cross-check the per-year assumption against the summary table.

    If the summary says the median peak SWE trend is -3.25 mm per decade and
    our per-station slopes times ten give the same thing, the assumption
    holds. If they disagree by more than a little, something is wrong and it
    is better to stop than to publish a figure that is out by a factor.
    """
    if not SUMMARY_FILE.exists():
        print(f"[note] {SUMMARY_FILE.name} not found, unit check skipped.")
        return
    summ = pd.read_csv(SUMMARY_FILE)
    if "metric" not in summ.columns or "obs_median_slope_per_decade" not in summ.columns:
        print("[note] summary table has unexpected columns, unit check skipped.")
        return

    print("\nunit check (per-station median x 10 against the summary table):")
    ok = True
    for key in METRICS:
        col = f"{key}_obs_slope"
        row = summ[summ["metric"] == key]
        if col not in df.columns or row.empty:
            continue
        ours = df[col].median() * SLOPE_TO_DECADE
        theirs = float(row["obs_median_slope_per_decade"].iloc[0])
        gap = abs(ours - theirs)
        tol = max(0.05, abs(theirs) * 0.02)
        mark = "ok " if gap <= tol else "MISMATCH"
        if gap > tol:
            ok = False
        print(f"  {mark} {key:<16} ours {ours:>8.2f}   summary {theirs:>8.2f}")
    if not ok:
        raise SystemExit(
            "\nThe per-station slopes do not scale to the summary values.\n"
            "Do not trust the figures. Send me this output and I will fix the\n"
            "conversion before anything is plotted."
        )


def fit_line(x, y):
    """Theil-Sen slope of y against x, with its CI, plus Kendall's tau."""
    good = np.isfinite(x) & np.isfinite(y)
    x, y = np.asarray(x)[good], np.asarray(y)[good]
    if len(x) < 10:
        return None
    slope, intercept, lo, hi = theilslopes(y, x, alpha=0.95)
    tau, p = kendalltau(x, y)
    return dict(slope=slope, intercept=intercept, lo=lo, hi=hi,
                tau=tau, p=p, n=len(x))


# ------------------------------------------------------------------ figure 1

def plot_trend_vs_elevation(df, key, nice, unit):
    col_slope, col_p = f"{key}_obs_slope", f"{key}_obs_p"
    if col_slope not in df.columns:
        print(f"  [skip] {key}: no column {col_slope}")
        return None

    d = df.dropna(subset=[col_slope, "elevation_m"]).copy()
    d["trend"] = d[col_slope] * SLOPE_TO_DECADE
    d["sig"] = d[col_p] < ALPHA if col_p in d.columns else False

    fig, ax = plt.subplots(figsize=(9, 6))
    ax.axhline(0, color="k", lw=.9, zorder=2)

    # the band edges, so the reader can see where the groups split
    for e in BAND_EDGES:
        ax.axvline(e, color="#999999", lw=.8, ls=":", zorder=1)

    # points, coloured by band, solid if the station's own trend is significant
    for b in range(3):
        sub = d[d["band"] == b]
        c = BAND_COLOURS[b]
        ns, s = sub[~sub["sig"]], sub[sub["sig"]]
        ax.scatter(ns.elevation_m, ns.trend, s=22, c=c, alpha=.40,
                   edgecolor="none", zorder=3)
        ax.scatter(s.elevation_m, s.trend, s=40, c=c, edgecolor="#222222",
                   linewidth=.5, zorder=4)

    # the median trend of each band, drawn across the band it belongs to
    lo_x, hi_x = d.elevation_m.min(), d.elevation_m.max()
    spans = [(lo_x, BAND_EDGES[0]), (BAND_EDGES[0], BAND_EDGES[1]),
             (BAND_EDGES[1], hi_x)]
    for b in range(3):
        sub = d[d["band"] == b]
        if sub.empty:
            continue
        m = sub["trend"].median()
        a, z = spans[b]
        ax.plot([a, z], [m, m], color=BAND_COLOURS[b], lw=3.4, solid_capstyle="butt",
                zorder=5, path_effects=None)
        ax.plot([a, z], [m, m], color="#222222", lw=.7, zorder=6)

    # the Theil-Sen line across all stations
    fit = fit_line(d.elevation_m.values, d["trend"].values)
    if fit:
        xs = np.array([lo_x, hi_x])
        ax.plot(xs, fit["intercept"] + fit["slope"] * xs,
                color="#c0392b", lw=1.9, zorder=7,
                label="Theil-Sen fit")
        per100 = fit["slope"] * 100
        lo100, hi100 = fit["lo"] * 100, fit["hi"] * 100
        txt = (f"Theil-Sen: {per100:+.2f} per 100 m\n"
               f"95% CI {lo100:+.2f} to {hi100:+.2f}\n"
               f"tau = {fit['tau']:+.3f},  p = {fit['p']:.2g}\n"
               f"n = {fit['n']} stations")
        ax.text(.985, .035, txt, transform=ax.transAxes, ha="right", va="bottom",
                fontsize=8.5, family="monospace",
                bbox=dict(boxstyle="round,pad=0.45", fc="white",
                          ec="#bbbbbb", alpha=.93), zorder=8)

    handles = []
    for b in range(3):
        n = int((d["band"] == b).sum())
        nsig = int(((d["band"] == b) & d["sig"]).sum())
        handles.append(Line2D([], [], marker="o", ls="", mfc=BAND_COLOURS[b],
                              mec="#222222", mew=.5, ms=8,
                              label=f"{BAND_NAMES[b]}  (n = {n}, {nsig} significant)"))
    handles.append(Line2D([], [], marker="o", ls="", mfc="#999999", mec="none",
                          ms=6, alpha=.45, label=f"Faded = not significant (p >= {ALPHA})"))
    handles.append(Line2D([], [], color="#c0392b", lw=1.9, label="Theil-Sen fit"))
    handles.append(Line2D([], [], color="#555555", lw=3.0, label="Band median"))
    ax.legend(handles=handles, loc="upper right", fontsize=8, framealpha=.93)

    ax.set_xlabel("Station elevation (m a.s.l.)")
    ax.set_ylabel(f"Observed {nice.lower()} trend ({unit})")
    ax.set_title(f"Observed {nice.lower()} trend against station elevation, {PERIOD}",
                 fontsize=12.5, pad=10)
    ax.grid(alpha=.18, lw=.4, ls=":")

    caption(fig, f"""
One point per station. The vertical position is that station's Theil-Sen trend in {unit}; the horizontal position is its elevation. Solid points with a dark edge are stations whose own trend is significant by the Mann-Kendall test at alpha = {ALPHA}; faded points are not.
Colour gives the elevation band. The thick horizontal segment in each band is that band's median trend, drawn across the elevation range the band covers. The red line is a Theil-Sen fit across all stations, which measures whether the rate of change itself depends on elevation; the median of pairwise slopes is used rather than least squares because 213 of 235 stations lie below 700 m and a least-squares line would be pulled by the few high ones.
Only stations with at least 0.8 of the period's hydrological years present are included, so a short or gappy record cannot produce a spurious trend.
""")
    fig.tight_layout()
    return fig


# ------------------------------------------------------------------ figure 2

def plot_bands_obs_vs_mod(df):
    keys = [k for k in METRICS if f"{k}_obs_slope" in df.columns]
    if not keys:
        return None
    n = len(keys)
    rows = (n + 1) // 2
    fig, axes = plt.subplots(rows, 2, figsize=(12, 4.6 * rows))
    axes = np.atleast_1d(axes).ravel()

    for ax, key in zip(axes, keys):
        nice, unit = METRICS[key]
        ax.axhline(0, color="k", lw=.9, zorder=1)
        positions, data, colours, dot_colours, ticks = [], [], [], [], []

        for b in range(3):
            sub = df[df["band"] == b]
            for j, src in enumerate(("obs", "mod")):
                col = f"{key}_{src}_slope"
                if col not in df.columns:
                    continue
                v = sub[col].dropna().values * SLOPE_TO_DECADE
                if len(v) == 0:
                    continue
                pos = b * 3 + j
                positions.append(pos)
                data.append(v)
                colours.append(SRC_COLOURS[src])
                dot_colours.append(SRC_DOT_COLOURS[src])
            ticks.append((b * 3 + 0.5,
                          f"{BAND_NAMES[b].split(' (')[0]}\nn = {len(sub)}"))

        bp = ax.boxplot(data, positions=positions, widths=.72, patch_artist=True,
                        showfliers=False, medianprops=dict(color="#222222", lw=1.6),
                        whiskerprops=dict(color="#555555", lw=.9),
                        capprops=dict(color="#555555", lw=.9))
        # Colour says observed or modelled. The modelled box is hatched as
        # well, so the pair is still readable if the figure is printed in
        # greyscale or photocopied.
        for patch, c, pos in zip(bp["boxes"], colours, positions):
            patch.set_facecolor(c)
            patch.set_edgecolor("#222222")
            patch.set_linewidth(.9)
            patch.set_alpha(.88)
            if pos % 3 == 1:                      # modelled
                patch.set_hatch("////")

        # The individual stations, one dot each, jittered sideways so they do
        # not land on top of one another. Coloured by source like the boxes,
        # so a dot can always be traced to observed or modelled.
        rng = np.random.default_rng(0)
        for pos, v, dc in zip(positions, data, dot_colours):
            ax.scatter(pos + rng.uniform(-.17, .17, len(v)), v, s=7,
                       c=dc, alpha=.55, linewidths=0, zorder=3)

        ax.set_xticks([t[0] for t in ticks])
        ax.set_xticklabels([t[1] for t in ticks])
        ax.set_ylabel(unit)
        ax.set_title(nice, fontsize=11.5)
        ax.grid(axis="y", alpha=.18, lw=.4, ls=":")

    for ax in axes[len(keys):]:
        ax.set_visible(False)

    # One legend for the whole figure, placed under the title so it never
    # lands on top of a box. The swatches carry the real styles.
    legend_handles = [
        Patch(facecolor=SRC_COLOURS["obs"], edgecolor="#222222", lw=.9,
              alpha=.88, label="Observed"),
        Patch(facecolor=SRC_COLOURS["mod"], edgecolor="#222222", lw=.9,
              alpha=.88, hatch="////", label="mHM modelled"),
    ]
    fig.legend(handles=legend_handles, loc="upper center",
               bbox_to_anchor=(0.5, 0.955), ncol=2, fontsize=10,
               frameon=True, framealpha=.95)

    fig.suptitle(f"Per-station trends by elevation band, observed against modelled, {PERIOD}",
                 fontsize=13, y=0.995)
    caption(fig, f"""
Each box holds the per-station Theil-Sen trends of the stations in that band: the box spans the interquartile range, the line inside is the median, the whiskers reach 1.5 times the interquartile range. The small dots are the individual stations, one per station, jittered sideways so they do not overlap; they carry the same colour as the box they belong to. Outliers are not drawn separately because the dots already show them.
In every pair the left green box is the observations and the right orange hatched box is the mHM simulation. A model that reproduced the observed change would put the orange box on top of the green one, in every band. Where the orange box is narrower, the model is producing less station-to-station variety than really exists; where it sits lower, the model is changing faster than the observations do.
The mountain band holds {int((df['band'] == 2).sum())} stations, so its quartiles are less stable than the other two. That is a limitation of the station network, not of the method.
""")
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    return fig


# ------------------------------------------------------------------ the table

def build_table(df):
    rows = []
    for key, (nice, unit) in METRICS.items():
        for b in range(3):
            sub = df[df["band"] == b]
            for src in ("obs", "mod"):
                cs, cp = f"{key}_{src}_slope", f"{key}_{src}_p"
                if cs not in df.columns:
                    continue
                v = sub[cs].dropna() * SLOPE_TO_DECADE
                if v.empty:
                    continue
                sig = sub[cp] < ALPHA if cp in sub.columns else pd.Series(dtype=bool)
                rows.append({
                    "metric": key,
                    "metric_name": nice,
                    "unit": unit,
                    "band": BAND_NAMES[b],
                    "band_from_m": (0 if b == 0 else BAND_EDGES[b - 1]),
                    "band_to_m": (BAND_EDGES[b] if b < len(BAND_EDGES) else np.nan),
                    "source": "observed" if src == "obs" else "modelled",
                    "n_stations": len(v),
                    "median_trend": round(float(v.median()), 3),
                    "q25": round(float(v.quantile(.25)), 3),
                    "q75": round(float(v.quantile(.75)), 3),
                    "n_significant": int(sig.sum()) if len(sig) else 0,
                    "n_significant_negative": int((sig & (sub[cs] < 0)).sum()) if len(sig) else 0,
                })
    return pd.DataFrame(rows)


def main():
    if not TRENDS_FILE.exists():
        raise SystemExit(f"Missing: {TRENDS_FILE}\nRun scripts/07_seasonal_statistics.py first.")

    df = pd.read_csv(TRENDS_FILE)
    df = df.dropna(subset=["elevation_m"]).copy()
    df["band"] = df["elevation_m"].apply(assign_band)

    print(f"stations: {len(df)}")
    for b in range(3):
        n = int((df["band"] == b).sum())
        sub = df[df["band"] == b]
        rng = f"{sub.elevation_m.min():.0f} to {sub.elevation_m.max():.0f} m" if n else "-"
        print(f"  {BAND_NAMES[b]:<22} n = {n:>3}   {rng}")

    verify_units(df)

    FIG_DIR.mkdir(parents=True, exist_ok=True)
    TAB_DIR.mkdir(parents=True, exist_ok=True)

    print("\nfigures:")
    for key, (nice, unit) in METRICS.items():
        fig = plot_trend_vs_elevation(df, key, nice, unit)
        if fig is None:
            continue
        out = FIG_DIR / f"trend_vs_elevation_{key}.png"
        fig.savefig(out, dpi=200, bbox_inches="tight")
        plt.close(fig)
        print(f"  {out.relative_to(REPO_ROOT)}")

    fig = plot_bands_obs_vs_mod(df)
    if fig is not None:
        out = FIG_DIR / "band_trends_obs_vs_mod.png"
        fig.savefig(out, dpi=200, bbox_inches="tight")
        plt.close(fig)
        print(f"  {out.relative_to(REPO_ROOT)}")

    tab = build_table(df)
    tab.to_csv(OUT_TABLE, index=False)
    print(f"\ntable: {OUT_TABLE.relative_to(REPO_ROOT)}  ({len(tab)} rows)")

    print("\nobserved median trend by band:")
    for key, (nice, unit) in METRICS.items():
        t = tab[(tab.metric == key) & (tab.source == "observed")]
        if t.empty:
            continue
        print(f"\n  {nice} ({unit})")
        for _, r in t.iterrows():
            print(f"    {r['band']:<22} n = {r['n_stations']:>3}   "
                  f"median {r['median_trend']:>8.2f}   "
                  f"IQR {r['q25']:>7.2f} to {r['q75']:>7.2f}   "
                  f"{r['n_significant']} significant")


if __name__ == "__main__":
    main()
