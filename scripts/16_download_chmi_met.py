"""
16_download_chmi_met.py

Downloads daily temperature and precipitation from the CHMI open data portal
for the same stations you already have snow water equivalent for.

Why
---
Both the station records and the mHM output currently hold snow only. To
show that the 400-700 m band is where winters now sit around freezing, and
to separate "less precipitation fell" from "it fell as rain", you need
temperature and precipitation. CHMI publishes both openly, in exactly the
same format as the SVH files your supervisor gave you, so this is a
download rather than a request.

What it fetches
---------------
    T    daily mean air temperature
    SRA  daily precipitation

from  opendata.chmi.cz/meteorology/climate/historical_csv/data/daily/

Where it puts them
------------------
    data/stations/t_raw/dly-<station>-T.csv
    data/stations/sra_raw/dly-<station>-SRA.csv

Same naming as svh_raw, so the reading code you already have needs only the
folder and element name changed.

How much
--------
SVH is winter-only and sparse. T and SRA are daily all year round, so they
are much larger. All 1147 stations would be roughly 3 GB. By default this
takes only the 235 stations that passed screening, roughly 700 MB. Set
STATIONS = "all" below if you want the lot.

Safe to interrupt
-----------------
Already-downloaded files are skipped, so if it stops halfway, or your
connection drops, just run it again and it picks up where it left off.

Before the full run it probes three stations to confirm the URL pattern,
and stops with a clear message if the portal's layout has changed, rather
than firing off two thousand failing requests.

Terms of use: CHMI data is open, and their licence asks that the source be
credited. Cite it in the thesis as Czech Hydrometeorological Institute
(CHMI) open data.

Run:
    cd /Users/aleynaturan/Desktop/a/aleyna-swe-thesis
    source .venv/bin/activate
    python scripts/16_download_chmi_met.py
"""

from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen
import re
import time

import pandas as pd

# ============================== settings ==============================

ELEMENTS = ["T", "SRA"]          # mean temperature, precipitation
STATIONS = "selected"            # "selected" (235) or "all" (1147)
PERIOD = "1976-2025"

BASE = ("https://opendata.chmi.cz/meteorology/climate/"
        "historical_csv/data/daily")

# The portal files daily data by category, not by element code. The full
# list of category folders is: air_pressure, cloudiness, evaporation,
# humidity, precipitation, radiation, snow, soil_temperature, sunshine,
# temperature, wind. Your SVH files came out of snow/.
FOLDER = {
    "T": "temperature", "TMA": "temperature", "TMI": "temperature",
    "SRA": "precipitation",
    "SVH": "snow", "SCE": "snow", "SNO": "snow",
    "H": "humidity", "F": "wind", "P": "air_pressure",
}
# Add "TMA" and "TMI" to ELEMENTS above if you want daily max and min
# temperature as well; everything else in this script already handles them.

WORKERS = 5                      # polite: a handful at a time, not hundreds
PAUSE = 0.15                     # seconds between requests per worker
TIMEOUT = 60
RETRIES = 3

# ======================================================================

REPO_ROOT = Path(__file__).resolve().parent.parent
STATIONS_DIR = REPO_ROOT / "data" / "stations"
SVH_DIR = STATIONS_DIR / "svh_raw"
MAP_FILE = REPO_ROOT / "results" / "tables" / "05_maps" / "map_stations.csv"
LOG_FILE = STATIONS_DIR / "download_log.csv"

UA = {"User-Agent": "bachelor-thesis-data-collection/1.0"}
RULE = "-" * 66


def banner(t):
    print(f"\n{RULE}\n{t}\n{RULE}")


def wsi_from_filename(name):
    """dly-0-203-0-42109026001-SVH.csv  ->  0-203-0-42109026001"""
    m = re.match(r"dly-(.+)-SVH\.csv$", name)
    return m.group(1) if m else None


def urls_for(wsi, element):
    """Category folder first, then two fallbacks in case the portal is
    reorganised, before a station is called missing."""
    fn = f"dly-{wsi}-{element}.csv"
    urls = []
    if element in FOLDER:
        urls.append(f"{BASE}/{FOLDER[element]}/{fn}")
    urls += [f"{BASE}/{element}/{fn}", f"{BASE}/{fn}"]
    return urls


def fetch(url):
    for attempt in range(RETRIES):
        try:
            with urlopen(Request(url, headers=UA), timeout=TIMEOUT) as r:
                return r.read()
        except HTTPError as e:
            if e.code == 404:
                return None              # this station has no such element
            if attempt == RETRIES - 1:
                raise
        except URLError:
            if attempt == RETRIES - 1:
                raise
        time.sleep(1.5 * (attempt + 1))
    return None


def get_one(wsi, element, out_dir):
    """Returns (status, bytes). Status: saved | skipped | missing | error."""
    out = out_dir / f"dly-{wsi}-{element}.csv"
    if out.exists() and out.stat().st_size > 0:
        return "skipped", out.stat().st_size
    for url in urls_for(wsi, element):
        try:
            data = fetch(url)
        except Exception as e:
            return f"error: {type(e).__name__}", 0
        if data:
            out.write_bytes(data)
            time.sleep(PAUSE)
            return "saved", len(data)
    return "missing", 0


def pick_stations():
    """Station codes to download, taken from the SVH files you already have."""
    files = sorted(SVH_DIR.glob("dly-*-SVH.csv"))
    if not files:
        raise SystemExit(f"No SVH files in {SVH_DIR}")
    all_wsi = [w for w in (wsi_from_filename(f.name) for f in files) if w]
    print(f"stations with SVH data: {len(all_wsi)}")

    if STATIONS == "all":
        print("STATIONS = 'all', taking every one of them")
        return all_wsi

    if not MAP_FILE.exists():
        print(f"[note] {MAP_FILE.name} not found, falling back to all stations")
        return all_wsi

    sel = set(pd.read_csv(MAP_FILE).query("selected")["station_id"].astype(str))
    print(f"stations selected for {PERIOD}: {len(sel)}")

    # The SVH filename carries the WMO station identifier; the analysis
    # tables carry CHMI's own station_id. The STATION column inside each
    # file links the two, so one header row per file is read to build it.
    print("matching filenames to station ids ...")
    keep = []
    for f in files:
        w = wsi_from_filename(f.name)
        if not w:
            continue
        try:
            sid = str(pd.read_csv(f, nrows=1)["STATION"].iloc[0])
        except Exception:
            continue
        if sid in sel:
            keep.append(w)

    if not keep:
        print("[note] no filename matched a selected station_id, so the two")
        print("       use different coding. Falling back to all stations.")
        return all_wsi

    print(f"matched: {len(keep)} of {len(sel)}")
    return keep


def main():
    banner("WHICH STATIONS")
    wsis = pick_stations()

    banner("CHECKING THE PORTAL")
    probe = wsis[:3]
    ok = False
    for w in probe:
        for url in urls_for(w, ELEMENTS[0]):
            try:
                hit = fetch(url) is not None
            except Exception as e:
                print(f"  {type(e).__name__} on {url}")
                continue
            print(f"  {'found   ' if hit else 'not here'}  {url}")
            ok = ok or hit
    if not ok:
        raise SystemExit(
            "\nNone of the probe URLs worked.\n"
            "Either the portal layout changed or there is no network here.\n"
            "Send me this output and I will fix the path before you run the\n"
            "full download.")
    print("\npattern confirmed, starting")

    rows = []
    for element in ELEMENTS:
        out_dir = STATIONS_DIR / f"{element.lower()}_raw"
        out_dir.mkdir(parents=True, exist_ok=True)
        banner(f"{element}  ->  {out_dir.relative_to(REPO_ROOT)}")

        done = total_bytes = 0
        counts = {}
        with ThreadPoolExecutor(max_workers=WORKERS) as pool:
            futs = {pool.submit(get_one, w, element, out_dir): w for w in wsis}
            for fut in as_completed(futs):
                w = futs[fut]
                try:
                    status, nbytes = fut.result()
                except Exception as e:
                    status, nbytes = f"error: {type(e).__name__}", 0
                counts[status.split(":")[0]] = counts.get(status.split(":")[0], 0) + 1
                total_bytes += nbytes
                rows.append({"element": element, "station": w, "status": status,
                             "bytes": nbytes})
                done += 1
                if done % 25 == 0 or done == len(wsis):
                    print(f"  {done:>5} / {len(wsis)}   "
                          f"{total_bytes / 1048576:>7.0f} MB   {counts}")

        print(f"\n  {element}: {counts}")

    LOG_FILE.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(LOG_FILE, index=False)

    banner("DONE")
    df = pd.DataFrame(rows)
    for element in ELEMENTS:
        sub = df[df.element == element]
        got = int((sub.status.isin(["saved", "skipped"])).sum())
        print(f"{element:>4}: {got} of {len(sub)} stations, "
              f"{sub.bytes.sum() / 1048576:.0f} MB")
    miss = df[df.status == "missing"]
    if len(miss):
        print(f"\n{len(miss)} station-element combinations not published "
              f"(that station does not record it)")
    err = df[df.status.str.startswith("error")]
    if len(err):
        print(f"{len(err)} failed with an error. Run the script again, it "
              f"skips what is already downloaded.")
    print(f"\nlog: {LOG_FILE.relative_to(REPO_ROOT)}")


if __name__ == "__main__":
    main()
