"""
06_station_qc_selection.py

Step 6: strict selection of stations for each candidate study period.

Coverage alone is not enough. Three very different stations can all show
"90% coverage":
  (a) measures every Monday, every season, for 50 years          -> good
  (b) measures densely but only every second season              -> bad
  (c) reported DAILY until the 1970s then switched to WEEKLY     -> dangerous:
      it contributes ~6x more rows in its early years, so any pooled mean,
      bias or RMSE is silently weighted towards that station and that era.

This script therefore (1) optionally homogenises the sampling frequency by
keeping Mondays only, (2) computes per-station diagnostics, (3) applies an
explicit multi-criteria filter, and (4) writes out WHY each station was
rejected, so nothing is dropped invisibly.

Run:
    source .venv/bin/activate
    python scripts/06_station_qc_selection.py

Inputs (already in the repo)
----------------------------
results/tables/modelled_vs_observed_swe.csv  (or .csv.gz - either is accepted)
results/tables/station_metadata.csv

Outputs (results/tables/ and results/figures/)
----------------------------------------------
qc_diagnostics_<period>.csv    every candidate station + all diagnostics + verdict
qc_selected_<period>.csv       the stations that PASS (this is the study set)
qc_rejection_summary.csv       how many stations failed on each criterion
qc_selection_overview.csv      final station count per period
season_presence_<period>.csv   station x season matrix of obs counts (eyeball it)
station_availability.png       availability heatmap (optional, needs matplotlib)
"""

from pathlib import Path

import numpy as np
import pandas as pd

# --- Config ---------------------------------------------------------------

# (first season, last season), labelled by season START year.
# Season 2024 = Nov 2024 - Apr 2025, the last COMPLETE season in the record
# (the data stops 2025-12-29), so "1956-2025" = seasons 1956..2024.
PERIODS = [(1956, 2024), (1966, 2024), (1976, 2024), (1986, 2024)]

SEASON_MONTHS = [11, 12, 1, 2, 3, 4]

# Homogenise sampling frequency.
#
# IMPORTANT: CHMI changed the weekly measurement day during 1970. Before the
# switch 92-94% of SWE observations fall on FRIDAYS; from the 1970/71 season
# onward 97-99% fall on MONDAYS. Filtering on "Mondays only" therefore
# deletes the whole pre-1970 record and makes the long periods look empty.
#
# "dominant" handles this: for each SEASON, the modal weekday across the
# whole network is found, and only observations on that weekday are kept.
# This gives every station-season the same nominal sampling (max ~26 obs),
# removes the daily-vs-weekly weighting problem, and survives the 1970
# switch automatically.
#
# Options: "dominant" (recommended) | "monday" | "none"
FREQUENCY_HOMOGENISATION = "dominant"

# --- Selection criteria (these are the knobs to argue about) --------------

MIN_OBS_PER_SEASON = 10     # a season counts as "active" at or above this
MIN_COVERAGE = 0.90         # fraction of the period's seasons that must be active
MAX_GAP_SEASONS = 3         # longest run of consecutive inactive seasons allowed
EDGE_TOLERANCE = 2          # must be active within this many seasons of each end
MIN_MEDIAN_OBS = 15         # median obs across its active seasons
MAX_CV_OBS = 0.60           # coeff. of variation of obs/season; catches stations
                            # whose measurement intensity swings wildly

REPO_ROOT = Path(__file__).resolve().parent.parent
TABLES = REPO_ROOT / "results" / "tables"
FIGURES = REPO_ROOT / "results" / "figures"
META_FILE = TABLES / "station_metadata.csv"


def find_obs_file() -> Path:
    """03_extract_modelled_swe.py writes an uncompressed .csv locally;
    the copy committed to the repo is gzipped. Accept either."""
    for name in ("modelled_vs_observed_swe.csv",
                 "modelled_vs_observed_swe.csv.gz"):
        p = TABLES / name
        if p.exists():
            return p
    raise FileNotFoundError(
        f"No observation table found in {TABLES}.\n"
        f"Expected 'modelled_vs_observed_swe.csv' or '.csv.gz'.\n"
        f"Run scripts/03_extract_modelled_swe.py first, or pull the gzipped "
        f"copy from the repo."
    )


MAKE_FIGURE = True


def period_label(start: int, end: int) -> str:
    return f"{start}-{end + 1}"


# --- Step 1: observations -> station x season counts ------------------------

DAY_NAMES = ["Mon", "Tue", "Wed", "Thu", "Fri", "Sat", "Sun"]


def season_counts(obs_path: Path) -> pd.DataFrame:
    df = pd.read_csv(obs_path, usecols=["station_id", "date"], parse_dates=["date"])

    df = df[df["date"].dt.month.isin(SEASON_MONTHS)].copy()
    df["season"] = np.where(df["date"].dt.month >= 11,
                            df["date"].dt.year,
                            df["date"].dt.year - 1)
    df["dow"] = df["date"].dt.dayofweek

    n_all = len(df)
    if FREQUENCY_HOMOGENISATION == "monday":
        df = df[df["dow"] == 0]
    elif FREQUENCY_HOMOGENISATION == "dominant":
        # modal weekday of each season, across the whole network
        modal = df.groupby("season")["dow"].agg(lambda s: s.value_counts().idxmax())
        modal.rename("season_dow").to_frame().assign(
            weekday=lambda x: [DAY_NAMES[d] for d in x["season_dow"]]
        ).to_csv(TABLES / "season_measurement_weekday.csv")
        df = df[df["dow"] == df["season"].map(modal)]
    elif FREQUENCY_HOMOGENISATION != "none":
        raise ValueError(f"unknown FREQUENCY_HOMOGENISATION: "
                         f"{FREQUENCY_HOMOGENISATION}")

    if FREQUENCY_HOMOGENISATION != "none":
        print(f"  homogenisation '{FREQUENCY_HOMOGENISATION}': kept {len(df)} "
              f"of {n_all} winter records ({100 * len(df) / n_all:.1f}%)")

    return (df.groupby(["station_id", "season"])
              .size()
              .rename("n_obs")
              .reset_index())


# --- Step 2: per-station diagnostics for one period -------------------------

def longest_gap(active_seasons: np.ndarray, start: int, end: int) -> int:
    """Longest run of consecutive INACTIVE seasons inside [start, end],
    including the stretches before the first and after the last active
    season."""
    all_seasons = np.arange(start, end + 1)
    is_active = np.isin(all_seasons, active_seasons)
    longest = run = 0
    for a in is_active:
        run = 0 if a else run + 1
        longest = max(longest, run)
    return int(longest)


def diagnose(counts: pd.DataFrame, start: int, end: int) -> pd.DataFrame:
    """One row per station that has ANY data in the period, with every
    diagnostic and a pass/fail verdict."""
    n_seasons = end - start + 1
    window = counts[(counts["season"] >= start) & (counts["season"] <= end)]

    rows = []
    for station_id, grp in window.groupby("station_id"):
        active = grp[grp["n_obs"] >= MIN_OBS_PER_SEASON]
        n_active = len(active)

        if n_active == 0:
            rows.append({
                "station_id": station_id,
                "n_active_seasons": 0, "coverage": 0.0,
                "first_active": np.nan, "last_active": np.nan,
                "max_gap_seasons": n_seasons,
                "median_obs": np.nan, "min_obs": np.nan, "max_obs": np.nan,
                "cv_obs": np.nan,
                "n_seasons_any_obs": int(len(grp)),
            })
            continue

        obs = active["n_obs"].to_numpy()
        rows.append({
            "station_id": station_id,
            "n_active_seasons": n_active,
            "coverage": n_active / n_seasons,
            "first_active": int(active["season"].min()),
            "last_active": int(active["season"].max()),
            "max_gap_seasons": longest_gap(active["season"].to_numpy(), start, end),
            "median_obs": float(np.median(obs)),
            "min_obs": int(obs.min()),
            "max_obs": int(obs.max()),
            "cv_obs": float(obs.std(ddof=0) / obs.mean()) if obs.mean() > 0 else np.nan,
            "n_seasons_any_obs": int(len(grp)),
        })

    d = pd.DataFrame(rows)

    # --- criteria, each as its own boolean column so failures are traceable
    d["ok_coverage"] = d["coverage"] >= MIN_COVERAGE
    d["ok_gap"] = d["max_gap_seasons"] <= MAX_GAP_SEASONS
    d["ok_starts_early"] = d["first_active"] <= start + EDGE_TOLERANCE
    d["ok_ends_late"] = d["last_active"] >= end - EDGE_TOLERANCE
    d["ok_density"] = d["median_obs"] >= MIN_MEDIAN_OBS
    d["ok_stability"] = d["cv_obs"] <= MAX_CV_OBS

    crit = ["ok_coverage", "ok_gap", "ok_starts_early",
            "ok_ends_late", "ok_density", "ok_stability"]
    d[crit] = d[crit].fillna(False)
    d["selected"] = d[crit].all(axis=1)

    d["fail_reasons"] = d.apply(
        lambda r: "" if r["selected"] else
        ";".join(c.removeprefix("ok_") for c in crit if not r[c]),
        axis=1,
    )

    d["period"] = period_label(start, end)
    d["n_seasons_in_period"] = n_seasons
    return d


# --- Step 3: figure ---------------------------------------------------------

def availability_figure(counts: pd.DataFrame, selected_ids: list, start: int,
                        end: int, out_path: Path) -> None:
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("  [skip figure] matplotlib not installed")
        return

    seasons = np.arange(start, end + 1)
    sub = counts[counts["station_id"].isin(selected_ids)]
    mat = (sub.pivot(index="station_id", columns="season", values="n_obs")
              .reindex(columns=seasons))
    mat = mat.loc[mat.notna().sum(axis=1).sort_values(ascending=False).index]

    fig, ax = plt.subplots(figsize=(11, max(4, len(mat) * 0.04)))
    im = ax.imshow(mat.to_numpy(), aspect="auto", interpolation="nearest",
                   cmap="viridis", vmin=0, vmax=26)
    ax.set_xticks(np.arange(0, len(seasons), 5))
    ax.set_xticklabels(seasons[::5], rotation=90, fontsize=7)
    ax.set_yticks([])
    ax.set_xlabel("Snow season (start year)")
    ax.set_ylabel(f"Selected stations (n = {len(mat)})")
    ax.set_title(f"Observations per season, selected stations, "
                 f"{period_label(start, end)}")
    fig.colorbar(im, ax=ax, label="observations in season")
    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print(f"  figure -> {out_path}")


# --- Main -------------------------------------------------------------------

def main():
    TABLES.mkdir(parents=True, exist_ok=True)
    obs_file = find_obs_file()
    print(f"Reading {obs_file.name} ...")
    counts = season_counts(obs_file)
    print(f"  {counts['station_id'].nunique()} stations, "
          f"{len(counts)} station-seasons")
    print(f"  obs per station-season: median {counts['n_obs'].median():.0f}, "
          f"min {counts['n_obs'].min()}, max {counts['n_obs'].max()}")

    meta = pd.read_csv(META_FILE, usecols=["station_id", "station_name", "latitude",
                                           "longitude", "elevation_m"])

    TABLES.mkdir(parents=True, exist_ok=True)
    overview, rejection = [], []

    for start, end in PERIODS:
        label = period_label(start, end)
        print(f"\n=== {label} ({end - start + 1} seasons) ===")

        d = diagnose(counts, start, end).merge(meta, on="station_id", how="left")
        sel = d[d["selected"]].copy()

        d.sort_values(["selected", "coverage"], ascending=[False, False]) \
         .to_csv(TABLES / f"qc_diagnostics_{label}.csv", index=False)
        sel.sort_values("elevation_m", ascending=False) \
           .to_csv(TABLES / f"qc_selected_{label}.csv", index=False)

        presence = (counts[(counts["season"] >= start) & (counts["season"] <= end)]
                    .pivot(index="station_id", columns="season", values="n_obs")
                    .reindex(columns=np.arange(start, end + 1)))
        presence.loc[presence.index.isin(sel["station_id"])] \
                .to_csv(TABLES / f"season_presence_{label}.csv")

        crit = ["coverage", "gap", "starts_early", "ends_late", "density", "stability"]
        rej = {"period": label, "candidates": int(len(d)), "selected": int(len(sel))}
        for c in crit:
            rej[f"failed_{c}"] = int((~d[f"ok_{c}"]).sum())
        rejection.append(rej)

        e = sel["elevation_m"].dropna()
        overview.append({
            "period": label,
            "n_seasons": end - start + 1,
            "n_selected": int(len(sel)),
            "elev_median_m": round(float(e.median()), 0) if len(e) else np.nan,
            "elev_max_m": round(float(e.max()), 0) if len(e) else np.nan,
            "n_above_800m": int((e > 800).sum()),
            "n_above_1000m": int((e > 1000).sum()),
            "median_coverage": round(float(sel["coverage"].median()), 3) if len(sel) else np.nan,
            "worst_coverage": round(float(sel["coverage"].min()), 3) if len(sel) else np.nan,
            "worst_gap": int(sel["max_gap_seasons"].max()) if len(sel) else 0,
            "median_obs_per_season": round(float(sel["median_obs"].median()), 1) if len(sel) else np.nan,
        })

        print(f"  candidates with any data : {len(d)}")
        print(f"  SELECTED                 : {len(sel)}")
        if len(sel):
            print(f"  worst coverage in set    : {sel['coverage'].min():.3f}")
            print(f"  worst gap in set         : {int(sel['max_gap_seasons'].max())} seasons")

        if MAKE_FIGURE and len(sel):
            availability_figure(counts, sel["station_id"].tolist(), start, end,
                                FIGURES / f"availability_{label}.png")

    ov = pd.DataFrame(overview)
    rj = pd.DataFrame(rejection)
    ov.to_csv(TABLES / "qc_selection_overview.csv", index=False)
    rj.to_csv(TABLES / "qc_rejection_summary.csv", index=False)

    print("\n=== Final selection ===")
    print(ov.to_string(index=False))
    print("\n=== Why stations were rejected (not mutually exclusive) ===")
    print(rj.to_string(index=False))
    print(f"\nSettings: homogenisation={FREQUENCY_HOMOGENISATION}, "
          f"min_obs_per_season={MIN_OBS_PER_SEASON}, min_coverage={MIN_COVERAGE}, "
          f"max_gap={MAX_GAP_SEASONS}, edge_tol={EDGE_TOLERANCE}, "
          f"min_median_obs={MIN_MEDIAN_OBS}, max_cv={MAX_CV_OBS}")


if __name__ == "__main__":
    main()
