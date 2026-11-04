"""
diagnose_blank_modelled.py

Read-only diagnostic. Changes nothing, writes one report.

Why this exists
---------------
The verification pass (99_verify_pipeline.py) found 104,299 rows in the paired
table where an observed SWE value is present but the modelled value is blank.
This script finds out WHERE those rows are, so we can say whether they affect
any result or not.

Two possible sources, with very different consequences:

  1. BENIGN - observations predating the model.
     The mHM file covers 1950-01-01 to 2026-05-31. CHMI observations in the
     table start in 1895. Every observation before 1950 therefore has no
     modelled counterpart. These rows are harmless: every study period starts
     in 1956, and scripts 07 and 08 both drop rows with no modelled value
     before computing anything.

  2. NOT BENIGN - stations whose grid cell is masked in the model.
     Script 03 selected stations by a lat/lon bounding box. A box around
     Czechia includes corners outside it, where the model grid may carry no
     data. Such a station returns NaN on every date, including dates inside
     the study period.

     This matters because of an asymmetry in the pipeline:
       - script 06 (screening) reads only station_id and date, so it counts
         observations and never sees whether the model has a value;
       - script 07 (metrics) drops rows with no modelled value.
     A station could therefore pass screening on its observation record while
     contributing fewer usable station-years than screening implied. No plot
     would reveal this.

The question this script answers
--------------------------------
Are there any blank-modelled rows inside the chosen study period
(1 Nov 1976 to 31 Oct 2025, hydrological years 1977-2025), at stations that
passed screening? If the answer is none, the failure is explained and nothing
needs changing. If the answer is some, we need to decide what to do about it.

Run:
    cd /Users/aleynaturan/Desktop/a/aleyna-swe-thesis
    source .venv/bin/activate
    python scripts/diagnose_blank_modelled.py
"""

from pathlib import Path

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parent.parent
TABLES = REPO_ROOT / "results" / "tables"
INPUTS = TABLES / "00_inputs"
SCREENING = TABLES / "01_data_screening"
REPORT = REPO_ROOT / "results" / "BLANK_MODELLED_DIAGNOSIS.txt"

MAIN_PERIOD = "1976-2025"
HY_FIRST, HY_LAST = 1977, 2025          # the hydrological years in that period
MODEL_START = pd.Timestamp("1950-01-01")  # from the mHM filename

_lines = []


def say(text=""):
    print(text)
    _lines.append(text)


def head(title):
    say()
    say("=" * 78)
    say(title)
    say("=" * 78)
    say()


def find_obs_file() -> Path:
    for name in ("modelled_vs_observed_swe.csv", "modelled_vs_observed_swe.csv.gz"):
        p = INPUTS / name
        if p.exists():
            return p
    raise FileNotFoundError(f"No paired table in {INPUTS}")


def main():
    src = find_obs_file()
    say("BLANK MODELLED VALUES — DIAGNOSIS")
    say(f"Reading {src.name} ...")
    df = pd.read_csv(src, parse_dates=["date"])

    df["blank"] = df["modelled_swe_mm"].isna()
    df["hydro_year"] = np.where(df["date"].dt.month >= 11,
                                df["date"].dt.year + 1, df["date"].dt.year)

    n_all, n_blank = len(df), int(df["blank"].sum())
    say(f"{n_all:,} rows total, {n_blank:,} with a blank modelled value "
        f"({100 * n_blank / n_all:.2f}%)")

    # ---- 1. are they simply before the model starts? ---------------------
    head("1. ARE THEY BEFORE THE MODEL PERIOD?")
    say(f"The mHM file starts {MODEL_START.date()}. Observations earlier than")
    say("that cannot have a modelled counterpart. If every blank row is here,")
    say("the explanation is benign.")
    say()
    before = df["blank"] & (df["date"] < MODEL_START)
    after = df["blank"] & (df["date"] >= MODEL_START)
    say(f"  blank AND before {MODEL_START.date()} : {int(before.sum()):>8,}")
    say(f"  blank AND after  {MODEL_START.date()} : {int(after.sum()):>8,}   <-- the ones that matter")
    say()
    if int(after.sum()) == 0:
        say("  VERDICT: every blank is an observation predating the model.")
        say("  Nothing downstream is affected.")
    else:
        say("  VERDICT: some blanks fall inside the model period. These need")
        say("  explaining — see sections 3 and 4.")

    # ---- 2. distribution over time ---------------------------------------
    head("2. WHEN DO THE BLANKS OCCUR?")
    say("Blank rows per decade, so the pattern is visible rather than assumed.")
    say()
    dec = (df[df["blank"]]["date"].dt.year // 10 * 10).value_counts().sort_index()
    tot = (df["date"].dt.year // 10 * 10).value_counts().sort_index()
    say(f"  {'decade':>8}  {'blank':>9}  {'all rows':>9}  {'% blank':>8}")
    for d in tot.index:
        b = int(dec.get(d, 0))
        say(f"  {d:>8}  {b:>9,}  {int(tot[d]):>9,}  {100 * b / tot[d]:>7.1f}%")

    # ---- 3. inside the chosen study period -------------------------------
    head(f"3. BLANKS INSIDE THE STUDY PERIOD ({MAIN_PERIOD}, HY {HY_FIRST}-{HY_LAST})")
    say("This is the question that decides whether anything must change.")
    say()
    inper = df[df["hydro_year"].between(HY_FIRST, HY_LAST)]
    n_in, n_in_blank = len(inper), int(inper["blank"].sum())
    say(f"  rows inside the period        : {n_in:>9,}")
    say(f"  of which blank modelled       : {n_in_blank:>9,}  "
        f"({100 * n_in_blank / n_in:.3f}%)")

    sel_f = SCREENING / f"qc_selected_{MAIN_PERIOD}.csv"
    if sel_f.exists():
        sel = set(pd.read_csv(sel_f)["station_id"])
        selper = inper[inper["station_id"].isin(sel)]
        n_sel, n_sel_blank = len(selper), int(selper["blank"].sum())
        say()
        say(f"  restricted to the {len(sel)} SELECTED stations:")
        say(f"    rows                        : {n_sel:>9,}")
        say(f"    of which blank modelled     : {n_sel_blank:>9,}  "
            f"({100 * n_sel_blank / n_sel:.3f}%)")
        say()
        if n_sel_blank == 0:
            say("    VERDICT: no selected station loses a single row inside the")
            say("    study period. Screening and metrics see the same data, and")
            say("    no result is affected.")
        else:
            say("    VERDICT: some selected stations lose rows that screening")
            say("    counted. Per-station detail below.")
            per = (selper.groupby("station_id")
                   .agg(rows=("blank", "size"), blank=("blank", "sum")))
            per["pct"] = 100 * per["blank"] / per["rows"]
            per = per[per["blank"] > 0].sort_values("blank", ascending=False)
            meta_f = INPUTS / "station_metadata.csv"
            if meta_f.exists():
                meta = pd.read_csv(meta_f, usecols=["station_id", "station_name",
                                                    "elevation_m"])
                per = per.merge(meta, on="station_id", how="left")
            say()
            say(per.head(25).to_string(index=False))
    else:
        say(f"  [skip] {sel_f.name} not found")

    # ---- 4. stations that are blank everywhere ---------------------------
    head("4. ARE ANY STATIONS BLANK ON EVERY DATE?")
    say("A station whose grid cell is masked in the model would return no")
    say("value on any date. Such a station should never have entered the")
    say("analysis at all.")
    say()
    per_st = df.groupby("station_id")["blank"].agg(["size", "sum"])
    allblank = per_st[per_st["sum"] == per_st["size"]]
    say(f"  stations with NO modelled value at all : {len(allblank)}")
    if len(allblank):
        meta_f = INPUTS / "station_metadata.csv"
        show = allblank.reset_index()
        if meta_f.exists():
            meta = pd.read_csv(meta_f, usecols=["station_id", "station_name",
                                                "latitude", "longitude",
                                                "elevation_m"])
            show = show.merge(meta, on="station_id", how="left")
        say()
        say(show.to_string(index=False))
        if sel_f.exists():
            overlap = set(allblank.index) & sel
            say()
            say(f"  of these, {len(overlap)} passed screening for {MAIN_PERIOD}"
                + (" <-- PROBLEM" if overlap else " (none, so no effect)"))

    # ---- 5. what scripts 07 and 08 actually do ---------------------------
    head("5. WHAT THE ANALYSIS SCRIPTS DO WITH THESE ROWS")
    say("Scripts 07 and 08 both call dropna(subset=['modelled_swe_mm']) when")
    say("they load the table, so blank rows are removed before any metric is")
    say("computed. Script 06 does not, because it reads only station_id and")
    say("date — it is measuring whether the STATION reported, which is a")
    say("question about the observations alone.")
    say()
    say("That difference is only a problem if blanks fall inside a study")
    say("period at a selected station. Section 3 above answers that.")

    REPORT.write_text("\n".join(_lines) + "\n")
    print(f"\nReport written to {REPORT}")


if __name__ == "__main__":
    main()
