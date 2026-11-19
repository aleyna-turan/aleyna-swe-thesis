"""
11_band_options.py

Read-only. Prints the elevation distribution of the selected stations and
shows how many stations each candidate banding scheme would put in each
band, so the choice of bands is made from the real numbers.

Writes nothing. Changes nothing. Safe to run as many times as you like.

Why this exists
---------------
Your supervisor asked for 2 or 3 elevation bands. "Equally divided" can
mean equal numbers of stations or equal spans of metres, and with a station
network that is concentrated in the lowlands these give very different
answers. A band holding five stations cannot support a trend estimate, so
we look before we choose.

It also lists what is already in results/tables/03_trends, so the next
script knows which per-station trend table to read.

Run:
    cd /Users/aleynaturan/Desktop/a/aleyna-swe-thesis
    source .venv/bin/activate
    python scripts/11_band_options.py
"""

from pathlib import Path

import numpy as np
import pandas as pd

PERIOD = "1976-2025"

REPO_ROOT = Path(__file__).resolve().parent.parent
TABLES = REPO_ROOT / "results" / "tables"
MAP_FILE = TABLES / "05_maps" / "map_stations.csv"
TRENDS_DIR = TABLES / "03_trends"

# Candidate fixed cuts, in metres. Round numbers that a reader can defend
# as lowland / upland / mountain in a central European context.
FIXED_2 = [600]
FIXED_3 = [400, 700]

RULE = "-" * 68


def banner(text):
    print(f"\n{RULE}\n{text}\n{RULE}")


def counts_for(edges, elev, labels=None):
    """Return a small table of band, range and station count."""
    lo, hi = elev.min(), elev.max()
    bounds = [lo - 1] + list(edges) + [hi + 1]
    rows = []
    for i in range(len(bounds) - 1):
        a, b = bounds[i], bounds[i + 1]
        sel = elev[(elev > a) & (elev <= b)]
        name = labels[i] if labels else f"band {i + 1}"
        rows.append({
            "band": name,
            "from_m": max(a, lo),
            "to_m": min(b, hi),
            "n": len(sel),
            "median_m": sel.median() if len(sel) else np.nan,
        })
    return pd.DataFrame(rows)


def show(title, table, note=""):
    print(f"\n{title}")
    if note:
        print(f"  {note}")
    out = table.copy()
    out["from_m"] = out["from_m"].round(0)
    out["to_m"] = out["to_m"].round(0)
    out["median_m"] = out["median_m"].round(0)
    print(out.to_string(index=False))
    smallest = int(table["n"].min())
    if smallest < 20:
        print(f"  [warning] smallest band has {smallest} stations. "
              f"A trend from that few is unstable.")


def main():
    if not MAP_FILE.exists():
        raise SystemExit(f"Missing: {MAP_FILE}\nRun scripts/10_build_map_table.py first.")

    st = pd.read_csv(MAP_FILE)
    sel = st[st["selected"]].dropna(subset=["elevation_m"]).copy()
    elev = sel["elevation_m"]

    banner(f"ELEVATION OF THE SELECTED STATIONS, {PERIOD}")
    print(f"stations: {len(sel)}")
    print(f"range   : {elev.min():.0f} to {elev.max():.0f} m")
    print(f"mean    : {elev.mean():.0f} m        median: {elev.median():.0f} m")

    print("\ndeciles (10% of the stations lie below each of these):")
    for q in range(10, 100, 10):
        print(f"  {q:>3}% : {np.percentile(elev, q):>6.0f} m")

    print("\nhow many stations above each round height:")
    for h in (400, 500, 600, 700, 800, 900, 1000, 1200):
        n = int((elev > h).sum())
        print(f"  > {h:>4} m : {n:>4} stations  ({100 * n / len(elev):>4.1f} %)")

    # ---------------------------------------------------------------- 2 bands
    banner("OPTION A - TWO BANDS")

    med = elev.median()
    show("A1. equal count (split at the median)",
         counts_for([med], elev, ["lower half", "upper half"]),
         f"cut at {med:.0f} m")

    mid = (elev.min() + elev.max()) / 2
    show("A2. equal width (split halfway up the elevation range)",
         counts_for([mid], elev, ["lower half", "upper half"]),
         f"cut at {mid:.0f} m")

    show("A3. fixed cut",
         counts_for(FIXED_2, elev, ["lowland", "upland"]),
         f"cut at {FIXED_2[0]} m")

    # -------------------------------------------------------------- 3 bands
    banner("OPTION B - THREE BANDS")

    t1, t2 = np.percentile(elev, [100 / 3, 200 / 3])
    show("B1. equal count (terciles)",
         counts_for([t1, t2], elev, ["lower third", "middle third", "upper third"]),
         f"cuts at {t1:.0f} and {t2:.0f} m")

    w1 = elev.min() + (elev.max() - elev.min()) / 3
    w2 = elev.min() + 2 * (elev.max() - elev.min()) / 3
    show("B2. equal width",
         counts_for([w1, w2], elev, ["low", "middle", "high"]),
         f"cuts at {w1:.0f} and {w2:.0f} m")

    show("B3. fixed cuts",
         counts_for(FIXED_3, elev, ["lowland", "upland", "mountain"]),
         f"cuts at {FIXED_3[0]} and {FIXED_3[1]} m")

    # ------------------------------------------------------- what trends exist
    banner("WHAT IS ALREADY IN results/tables/03_trends")
    if not TRENDS_DIR.exists():
        print("That folder does not exist.")
    else:
        files = sorted(TRENDS_DIR.glob("*.csv"))
        if not files:
            print("No CSV files in there.")
        for f in files:
            try:
                head = pd.read_csv(f, nrows=1)
                n = sum(1 for _ in open(f)) - 1
                print(f"\n{f.name}   ({n} rows)")
                print(f"  columns: {list(head.columns)}")
            except Exception as e:
                print(f"\n{f.name}   [could not read: {e}]")

    banner("NEXT")
    print("Paste this whole output back. We choose the banding from these")
    print("numbers, then I write the plotting script.")


if __name__ == "__main__":
    main()
