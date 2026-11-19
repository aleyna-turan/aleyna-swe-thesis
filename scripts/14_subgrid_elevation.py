"""
14_subgrid_elevation.py

Tests whether mHM's peak SWE bias is explained by the difference between a
station's own elevation and the mean elevation of the 12 km model cell it
sits in.

The idea
--------
mHM runs on 0.125 degree cells, about 9 km east-west by 14 km north-south at
Czech latitudes. Each cell has one elevation: the average of everything
inside it. A station, though, sits at one point in that cell, and that point
is almost never at the average height.

    dz = station elevation - mean elevation of its mHM cell

A station on a summit has a large positive dz. The model is simulating a
much lower, warmer place than where the measurement was taken, so it should
produce too little snow there: a negative bias. A station in a valley floor
has a negative dz, the model simulates a higher, colder place, and should
produce too much snow: a positive bias.

So the prediction, written down before looking:

    bias should become more negative as dz becomes more positive.

You already know bias becomes more negative with station elevation
(-2.29 mm per 100 m, p = 3.3e-5). That is the symptom. dz would be the
cause, because high stations are disproportionately on ridges and summits
and therefore have large positive dz.

The regression that settles it is the third one below: fit bias against
station elevation and dz together. If dz keeps its coefficient while station
elevation collapses toward zero, the subgrid explanation wins. If station
elevation survives, something else is driving the bias (precipitation
gradients, snow physics, gauge undercatch) and that is a different and
equally real finding.

The grids
---------
The DEM is on a regular lat/lon grid of 0.03125 degrees. The mHM grid is
0.125 degrees. 0.125 / 0.03125 = 4 exactly, and the cell edges line up, so
every model cell contains exactly 4 x 4 = 16 DEM cells and none straddles a
boundary. No reprojection, no interpolation, no area weighting: a plain
block mean is exact here. The script checks this rather than assuming it.

Sixteen samples is a thin estimate of a cell's mean elevation, and thinnest
where the terrain is roughest, which is where it matters most. The script
reports the within-cell standard deviation so the thinness is visible.

Outputs
-------
results/tables/02_model_evaluation/subgrid_elevation_<period>.csv
results/figures/02_model_evaluation/bias_vs_subgrid_elevation.png

Run:
    cd /Users/aleynaturan/Desktop/a/aleyna-swe-thesis
    source .venv/bin/activate
    python scripts/14_subgrid_elevation.py
"""

from pathlib import Path
import textwrap

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr
from matplotlib.lines import Line2D
from scipy.stats import theilslopes, kendalltau, t as student_t

# ============================== settings ==============================

PERIOD = "1976-2025"

# Same bands as script 12, so the two figures can be read together.
BAND_EDGES = [400, 700]
BAND_NAMES = ["Lowland (<= 400 m)", "Upland (400-700 m)", "Mountain (> 700 m)"]
BAND_COLOURS = ["#a1dab4", "#41b6c4", "#225ea8"]

DEM_VAR = "dem"
DEM_LAT, DEM_LON = "latitude", "longitude"
MHM_LAT, MHM_LON = "lat", "lon"

# A cell needs at least this many valid DEM pixels for its mean to be used.
MIN_DEM_PIXELS = 8          # out of the 16 that should be there

# ======================================================================

REPO_ROOT = Path(__file__).resolve().parent.parent
DATA_DIR = REPO_ROOT / "data"
DEM_FILE = DATA_DIR / "dem_hicam_3km_eu.nc"
MAP_FILE = REPO_ROOT / "results" / "tables" / "05_maps" / "map_stations.csv"
TAB_DIR = REPO_ROOT / "results" / "tables" / "02_model_evaluation"
FIG_DIR = REPO_ROOT / "results" / "figures" / "02_model_evaluation"
OUT_TABLE = TAB_DIR / f"subgrid_elevation_{PERIOD}.csv"
OUT_FIG = FIG_DIR / "bias_vs_subgrid_elevation.png"

RULE = "-" * 70


def banner(text):
    print(f"\n{RULE}\n{text}\n{RULE}")


def find_model_file():
    """Prefer the tiny time-average file; only the grid is needed from it."""
    cands = sorted(DATA_DIR.rglob("snowpack_czechia*timavg*.nc"))
    cands += sorted(DATA_DIR.rglob("snowpack_czechia*.nc"))
    for c in cands:
        return c
    raise SystemExit(f"No snowpack NetCDF found under {DATA_DIR}")


def edges_from(ds, name):
    """Ascending cell edges for one axis, from its bounds if they exist."""
    bnds = f"{name}_bnds"
    if bnds in ds.variables:
        arr = np.asarray(ds[bnds].values, dtype=float)
        return np.unique(np.concatenate([arr[:, 0], arr[:, 1]]))
    centres = np.sort(np.asarray(ds[name].values, dtype=float))
    step = float(np.median(np.diff(centres)))
    return np.append(centres - step / 2, centres[-1] + step / 2)


def cell_index(edges, values):
    """Which cell each value falls in. -1 means outside the grid."""
    idx = np.searchsorted(edges, values, side="right") - 1
    idx[(values < edges[0]) | (values >= edges[-1])] = -1
    return idx


def ols(y, X, names):
    """Least squares with standard errors, t statistics and p values.

    Written out rather than imported so the script has no dependency on
    statsmodels, and so every number in the table can be traced.
    """
    n, k = X.shape
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    resid = y - X @ beta
    dof = n - k
    sigma2 = float(resid @ resid) / dof
    XtX_inv = np.linalg.inv(X.T @ X)
    se = np.sqrt(np.diag(XtX_inv) * sigma2)
    tstat = beta / se
    p = 2 * (1 - student_t.cdf(np.abs(tstat), dof))
    ss_tot = float(((y - y.mean()) ** 2).sum())
    r2 = 1 - float(resid @ resid) / ss_tot if ss_tot > 0 else np.nan
    return pd.DataFrame({"term": names, "coef": beta, "std_err": se,
                         "t": tstat, "p": p}), r2, n


def report_model(title, table, r2, n, scale_note=""):
    print(f"\n{title}   (n = {n}, R2 = {r2:.3f})")
    if scale_note:
        print(f"  {scale_note}")
    for _, r in table.iterrows():
        star = "*" if r["p"] < 0.05 else " "
        # A p value printed as 0 is floating point underflow, not a real zero.
        # Report it as a bound so it is never written up as "p = 0".
        pstr = "< 1e-16" if r["p"] <= 0 else f"{r['p']:.3g}"
        print(f"  {star} {r['term']:<22} {r['coef']:>9.3f}  "
              f"+/- {r['std_err']:>7.3f}   t = {r['t']:>7.2f}   p = {pstr}")


def caption(fig, text):
    width = int(fig.get_figwidth() * 15)
    wrapped = "\n".join(textwrap.fill(p, width) for p in text.strip().split("\n"))
    fig.text(0.01, -0.01, wrapped, ha="left", va="top", fontsize=7.5, color="#333333")


# ----------------------------------------------------------------------

def main():
    for f in (DEM_FILE, MAP_FILE):
        if not f.exists():
            raise SystemExit(f"Missing: {f}")

    model_file = find_model_file()
    banner("GRIDS")
    print(f"DEM   : {DEM_FILE.relative_to(REPO_ROOT)}")
    print(f"model : {model_file.relative_to(REPO_ROOT)}")

    # --- the model grid --------------------------------------------------
    with xr.open_dataset(model_file, decode_times=False) as mds:
        lat_edges = edges_from(mds, MHM_LAT)
        lon_edges = edges_from(mds, MHM_LON)
    nlat, nlon = len(lat_edges) - 1, len(lon_edges) - 1
    mhm_step_lat = float(np.median(np.diff(lat_edges)))
    mhm_step_lon = float(np.median(np.diff(lon_edges)))
    print(f"\nmodel grid : {nlat} x {nlon} cells")
    print(f"  longitude {lon_edges[0]:.4f} to {lon_edges[-1]:.4f}, step {mhm_step_lon:.5f} deg")
    print(f"  latitude  {lat_edges[0]:.4f} to {lat_edges[-1]:.4f}, step {mhm_step_lat:.5f} deg")

    # --- the DEM ---------------------------------------------------------
    with xr.open_dataset(DEM_FILE, decode_times=False) as dds:
        dlat = np.asarray(dds[DEM_LAT].values, dtype=float)
        dlon = np.asarray(dds[DEM_LON].values, dtype=float)
        dem_step = float(np.median(np.abs(np.diff(dlat))))
        # only the part of the DEM that can fall inside the model grid
        keep_lat = (dlat >= lat_edges[0] - dem_step) & (dlat <= lat_edges[-1] + dem_step)
        keep_lon = (dlon >= lon_edges[0] - dem_step) & (dlon <= lon_edges[-1] + dem_step)
        sub = dds[DEM_VAR].isel({DEM_LAT: np.where(keep_lat)[0],
                                 DEM_LON: np.where(keep_lon)[0]})
        z = np.asarray(sub.values, dtype=float)
        zlat = np.asarray(sub[DEM_LAT].values, dtype=float)
        zlon = np.asarray(sub[DEM_LON].values, dtype=float)

    print(f"\nDEM grid   : step {dem_step:.5f} deg")
    ratio = mhm_step_lon / dem_step
    nested = abs(ratio - round(ratio)) < 1e-6
    print(f"  model step / DEM step = {ratio:.4f}  "
          f"({'exact nesting, ' + str(round(ratio)) + ' x ' + str(round(ratio)) + ' DEM cells per model cell' if nested else 'NOT an integer'})")
    print(f"  DEM window over the model domain: {z.shape[0]} x {z.shape[1]} pixels, "
          f"{100 * np.isfinite(z).mean():.1f}% valid")
    if not nested:
        print("  [note] the grids do not nest exactly. The pixel-centre assignment")
        print("         used below is still correct, just not a clean block mean.")

    # --- aggregate the DEM onto the model grid ---------------------------
    ilat = cell_index(lat_edges, zlat)
    ilon = cell_index(lon_edges, zlon)
    LAT, LON = np.meshgrid(ilat, ilon, indexing="ij")
    good = np.isfinite(z) & (LAT >= 0) & (LON >= 0)
    ids = (LAT[good] * nlon + LON[good]).astype(np.int64)
    vals = z[good]

    ncells = nlat * nlon
    cnt = np.bincount(ids, minlength=ncells).astype(float)
    s1 = np.bincount(ids, weights=vals, minlength=ncells)
    s2 = np.bincount(ids, weights=vals ** 2, minlength=ncells)
    with np.errstate(invalid="ignore", divide="ignore"):
        cell_mean = s1 / cnt
        cell_var = np.clip(s2 / cnt - cell_mean ** 2, 0, None)
    cell_std = np.sqrt(cell_var)
    cell_mean[cnt < MIN_DEM_PIXELS] = np.nan
    cell_std[cnt < MIN_DEM_PIXELS] = np.nan

    filled = int(np.isfinite(cell_mean).sum())
    print(f"\ncells with a usable mean elevation: {filled} of {ncells}")
    print(f"  DEM pixels per cell: median {np.median(cnt[cnt > 0]):.0f}, "
          f"min {cnt[cnt > 0].min():.0f}, max {cnt.max():.0f}")

    # --- stations --------------------------------------------------------
    banner("STATIONS")
    st = pd.read_csv(MAP_FILE)
    st = st[st["selected"]].dropna(subset=["latitude", "longitude", "elevation_m"]).copy()
    print(f"selected stations: {len(st)}")

    s_ilat = cell_index(lat_edges, st["latitude"].values)
    s_ilon = cell_index(lon_edges, st["longitude"].values)
    inside = (s_ilat >= 0) & (s_ilon >= 0)
    if (~inside).sum():
        print(f"  [warning] {int((~inside).sum())} station(s) fall outside the model grid")
    sid = np.where(inside, s_ilat * nlon + s_ilon, -1)

    st["cell_row"] = np.where(inside, s_ilat, -1)
    st["cell_col"] = np.where(inside, s_ilon, -1)
    st["grid_elevation_m"] = [cell_mean[i] if i >= 0 else np.nan for i in sid]
    st["grid_elevation_std_m"] = [cell_std[i] if i >= 0 else np.nan for i in sid]
    st["dem_pixels_in_cell"] = [cnt[i] if i >= 0 else np.nan for i in sid]
    st["dz_m"] = st["elevation_m"] - st["grid_elevation_m"]

    d = st.dropna(subset=["dz_m", "bias_peak_swe_mm"]).copy()
    print(f"stations with both dz and a bias value: {len(d)}")

    share = d.groupby(["cell_row", "cell_col"]).size()
    print(f"model cells holding at least one station: {len(share)}")
    print(f"  stations per occupied cell: median {share.median():.0f}, max {share.max()}")

    print(f"\ndz = station elevation - grid mean elevation")
    print(f"  range  {d.dz_m.min():>8.1f} to {d.dz_m.max():>8.1f} m")
    print(f"  mean   {d.dz_m.mean():>8.1f} m     median {d.dz_m.median():>8.1f} m")
    print(f"  above the grid mean: {int((d.dz_m > 0).sum())} stations, "
          f"below: {int((d.dz_m < 0).sum())}")
    print(f"\nwithin-cell elevation spread (how rough the terrain is):")
    print(f"  median {d.grid_elevation_std_m.median():.0f} m, "
          f"max {d.grid_elevation_std_m.max():.0f} m")

    r_ed = np.corrcoef(d.elevation_m, d.dz_m)[0, 1]
    print(f"\ncorrelation between station elevation and dz: r = {r_ed:+.3f}")
    print("  (they are related, which is exactly why the joint regression below")
    print("   is needed to say which of the two carries the bias signal)")

    # --- the named stations ---------------------------------------------
    if "notable" in d.columns:
        lab = d[d.notable.fillna("") != ""].sort_values("elevation_m", ascending=False)
        if len(lab):
            banner("THE NAMED STATIONS")
            print(f"{'station':<24}{'elev':>7}{'grid':>8}{'dz':>8}"
                  f"{'rough':>8}{'bias':>9}")
            print(f"{'':<24}{'m':>7}{'m':>8}{'m':>8}{'m':>8}{'mm':>9}")
            for _, r in lab.iterrows():
                print(f"{str(r.station_name)[:23]:<24}{r.elevation_m:>7.0f}"
                      f"{r.grid_elevation_m:>8.0f}{r.dz_m:>8.0f}"
                      f"{r.grid_elevation_std_m:>8.0f}{r.bias_peak_swe_mm:>9.1f}")

    # --- the statistics --------------------------------------------------
    banner("DOES BIAS SCALE WITH dz?")

    ts = theilslopes(d.bias_peak_swe_mm.values, d.dz_m.values, alpha=0.95)
    tau, p_tau = kendalltau(d.dz_m.values, d.bias_peak_swe_mm.values)
    print(f"Theil-Sen, bias against dz:")
    print(f"  {ts[0] * 100:+.3f} mm per 100 m of dz")
    print(f"  95% CI {ts[2] * 100:+.3f} to {ts[3] * 100:+.3f}")
    print(f"  Kendall tau = {tau:+.3f}, p = {p_tau:.3g}")
    print(f"\n  prediction was a NEGATIVE slope. "
          f"Observed slope is {'NEGATIVE, as predicted' if ts[0] < 0 else 'POSITIVE, against the prediction'}.")

    y = d.bias_peak_swe_mm.values.astype(float)
    e = d.elevation_m.values.astype(float) / 100.0      # per 100 m
    g = d.dz_m.values.astype(float) / 100.0             # per 100 m
    one = np.ones_like(y)

    banner("WHICH ONE CARRIES THE SIGNAL?")
    print("Coefficients are mm of peak SWE bias per 100 m. A star marks p < 0.05.")

    t1, r2_1, n1 = ols(y, np.column_stack([one, e]), ["intercept", "station elevation"])
    report_model("model 1: bias ~ station elevation", t1, r2_1, n1)

    t2, r2_2, n2 = ols(y, np.column_stack([one, g]), ["intercept", "dz"])
    report_model("model 2: bias ~ dz", t2, r2_2, n2)

    t3, r2_3, n3 = ols(y, np.column_stack([one, e, g]),
                       ["intercept", "station elevation", "dz"])
    report_model("model 3: bias ~ station elevation + dz", t3, r2_3, n3)

    c1 = float(t1.loc[t1.term == "station elevation", "coef"].iloc[0])
    c3 = float(t3.loc[t3.term == "station elevation", "coef"].iloc[0])
    p3_dz = float(t3.loc[t3.term == "dz", "p"].iloc[0])
    p3_el = float(t3.loc[t3.term == "station elevation", "p"].iloc[0])

    print("\nhow to read model 3:")
    print(f"  station elevation on its own : {c1:+.2f} mm per 100 m")
    print(f"  station elevation beside dz  : {c3:+.2f} mm per 100 m "
          f"(p = {p3_el:.3g})")
    shrink = (1 - abs(c3) / abs(c1)) * 100 if c1 != 0 else np.nan
    print(f"  it shrank by {shrink:.0f} percent once dz was included")
    if p3_dz < 0.05 and p3_el >= 0.05:
        print("  -> dz carries the signal and station elevation does not survive.")
        print("     The subgrid explanation is supported.")
    elif p3_dz < 0.05 and p3_el < 0.05:
        print("  -> both survive. Subgrid elevation explains part of the bias,")
        print("     but something else contributes as well.")
    elif p3_dz >= 0.05:
        print("  -> dz does not survive. Subgrid elevation is not the explanation,")
        print("     and the elevation-bias relationship needs another one.")

    print("\nhow much each one adds, which the p values do not show:")
    print(f"  station elevation alone      R2 = {r2_1:.3f}")
    print(f"  dz alone                     R2 = {r2_2:.3f}")
    print(f"  both together                R2 = {r2_3:.3f}")
    print(f"  adding dz to elevation       R2 {r2_1:.3f} -> {r2_3:.3f}  "
          f"(+{100 * (r2_3 - r2_1):.0f} points)")
    print(f"  adding elevation to dz       R2 {r2_2:.3f} -> {r2_3:.3f}  "
          f"(+{100 * (r2_3 - r2_2):.0f} points)")

    # --- does it survive without the most extreme station? ----------------
    # One station with a far larger dz than any other has heavy leverage on a
    # least-squares fit, so the whole conclusion could rest on it. This drops
    # it and refits. That is a data reason for removal, not a cosmetic one:
    # the point is to find out whether the result depends on a single row.
    banner("DOES IT SURVIVE WITHOUT THE MOST EXTREME STATION?")
    k = d["dz_m"].abs().idxmax()
    print(f"dropping {d.loc[k, 'station_name']}: "
          f"dz = {d.loc[k, 'dz_m']:+.0f} m, bias = {d.loc[k, 'bias_peak_swe_mm']:+.1f} mm")
    others = d["dz_m"].abs().drop(index=k)
    print(f"  next largest |dz| in the network is {others.max():.0f} m, "
          f"so this station sits far outside the rest")

    dd = d.drop(index=k)
    ts_d = theilslopes(dd.bias_peak_swe_mm.values, dd.dz_m.values, alpha=0.95)
    tau_d, p_tau_d = kendalltau(dd.dz_m.values, dd.bias_peak_swe_mm.values)
    print(f"\nTheil-Sen, bias against dz, without it:")
    print(f"  {ts_d[0] * 100:+.3f} mm per 100 m   "
          f"(with it: {ts[0] * 100:+.3f})")
    print(f"  95% CI {ts_d[2] * 100:+.3f} to {ts_d[3] * 100:+.3f}")
    print(f"  Kendall tau = {tau_d:+.3f}, p = {p_tau_d:.3g}")

    yd = dd.bias_peak_swe_mm.values.astype(float)
    ed = dd.elevation_m.values.astype(float) / 100.0
    gd = dd.dz_m.values.astype(float) / 100.0
    od = np.ones_like(yd)
    t3d, r2_3d, n3d = ols(yd, np.column_stack([od, ed, gd]),
                          ["intercept", "station elevation", "dz"])
    report_model("model 3 refitted without it: bias ~ station elevation + dz",
                 t3d, r2_3d, n3d)

    c3d_dz = float(t3d.loc[t3d.term == "dz", "coef"].iloc[0])
    c3_dz = float(t3.loc[t3.term == "dz", "coef"].iloc[0])
    p3d_dz = float(t3d.loc[t3d.term == "dz", "p"].iloc[0])
    print(f"\n  dz coefficient with it    : {c3_dz:+.2f} mm per 100 m")
    print(f"  dz coefficient without it : {c3d_dz:+.2f} mm per 100 m")
    if p3d_dz < 0.05 and np.sign(c3d_dz) == np.sign(c3_dz):
        print("  -> the relationship holds without it. The conclusion does not")
        print("     rest on one station.")
    else:
        print("  -> the relationship does NOT hold without it. The conclusion")
        print("     rested on that single station and must be stated as such.")

    # --- figure ----------------------------------------------------------
    d["band"] = np.searchsorted(BAND_EDGES, d["elevation_m"].values, side="left")

    fig, ax = plt.subplots(figsize=(9.5, 6.4))
    ax.axhline(0, color="k", lw=.9, zorder=2)
    ax.axvline(0, color="#999999", lw=.9, ls=":", zorder=1)

    for b in range(3):
        sub_b = d[d["band"] == b]
        ax.scatter(sub_b.dz_m, sub_b.bias_peak_swe_mm, s=34, c=BAND_COLOURS[b],
                   edgecolor="#222222", linewidth=.4, alpha=.9, zorder=3)

    xs = np.array([d.dz_m.min(), d.dz_m.max()])
    ax.plot(xs, ts[1] + ts[0] * xs, color="#c0392b", lw=1.9, zorder=5)

    # Pad the limits so that no station sits on the frame. Without this the
    # most extreme station, which is the one a reader most wants to see,
    # lands in the corner and its label falls outside the axes.
    xlo, xhi = float(d.dz_m.min()), float(d.dz_m.max())
    ylo, yhi = float(d.bias_peak_swe_mm.min()), float(d.bias_peak_swe_mm.max())
    xpad, ypad = (xhi - xlo) * 0.09, (yhi - ylo) * 0.09
    ax.set_xlim(xlo - xpad, xhi + xpad)
    ax.set_ylim(ylo - ypad, yhi + ypad)

    # Bottom left. The bottom right is where the largest positive dz lands,
    # which is the single station a reader most wants to see.
    ax.text(.015, .035,
            f"Theil-Sen: {ts[0] * 100:+.2f} mm per 100 m of dz\n"
            f"95% CI {ts[2] * 100:+.2f} to {ts[3] * 100:+.2f}\n"
            f"tau = {tau:+.3f},  p = {p_tau:.2g}\n"
            f"n = {len(d)} stations",
            transform=ax.transAxes, ha="left", va="bottom", fontsize=8.5,
            family="monospace",
            bbox=dict(boxstyle="round,pad=0.45", fc="white", ec="#bbbbbb",
                      alpha=.93), zorder=8)

    if "notable" in d.columns:
        lab = d[d.notable.fillna("") != ""]
        xmid = sum(ax.get_xlim()) / 2
        ymid = sum(ax.get_ylim()) / 2
        for i, (_, r) in enumerate(lab.iterrows()):
            # Push every label toward the middle of the plot, so one near an
            # edge is never written off the figure.
            to_left = r.dz_m > xmid
            ox, ha = (-20, "right") if to_left else (20, "left")
            oy = 15 if r.bias_peak_swe_mm < ymid else -15
            va = "bottom" if oy > 0 else "top"
            ax.annotate(str(r.station_name), xy=(r.dz_m, r.bias_peak_swe_mm),
                        xytext=(ox, oy), textcoords="offset points",
                        fontsize=8, ha=ha, va=va,
                        bbox=dict(boxstyle="round,pad=0.25", fc="white",
                                  ec="#999999", lw=.5, alpha=.92),
                        arrowprops=dict(arrowstyle="-", lw=.7, color="#444444",
                                        shrinkA=0, shrinkB=3),
                        zorder=7)

    handles = [Line2D([], [], marker="o", ls="", mfc=BAND_COLOURS[b],
                      mec="#222222", mew=.4, ms=8,
                      label=f"{BAND_NAMES[b]}  (n = {int((d['band'] == b).sum())})")
               for b in range(3)]
    handles.append(Line2D([], [], color="#c0392b", lw=1.9, label="Theil-Sen fit"))
    ax.legend(handles=handles, loc="upper right", fontsize=8.5, framealpha=.93)

    ax.set_xlabel("dz  =  station elevation minus mHM cell mean elevation  (m)")
    ax.set_ylabel("mHM peak SWE bias, modelled minus observed  (mm)")
    ax.set_title("Model bias against subgrid elevation difference, "
                 f"{PERIOD}", fontsize=12.5, pad=10)
    ax.grid(alpha=.18, lw=.4, ls=":")

    caption(fig, f"""
One point per station. The horizontal axis is dz, how far the station sits above or below the mean elevation of the 0.125 degree mHM cell that contains it; the cell mean is the average of the {int(round(ratio)) if nested else 'contained'} by {int(round(ratio)) if nested else ''} DEM cells inside it, from the 0.03125 degree HICAM DEM. A station to the right of the dotted line is higher than the terrain the model resolves there; a station to the left is lower.
The vertical axis is the model's peak SWE bias, modelled minus observed, so points below the solid line are stations where mHM produces too little snow.
If subgrid elevation drives the bias, points should run from upper left to lower right: stations the model places too low should be short of snow. The red line is a Theil-Sen fit, used rather than least squares because dz has a long positive tail at the mountain stations. Colour gives the elevation band, so the relationship between the two explanations stays visible.
Each model cell is averaged from only {int(round(ratio ** 2)) if nested else 'a few'} DEM pixels, so dz carries real uncertainty, largest where terrain is roughest.
""")
    fig.tight_layout()

    TAB_DIR.mkdir(parents=True, exist_ok=True)
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_FIG, dpi=200, bbox_inches="tight")
    plt.close(fig)

    keep = ["station_id", "station_name", "latitude", "longitude", "elevation_m",
            "grid_elevation_m", "grid_elevation_std_m", "dz_m",
            "dem_pixels_in_cell", "cell_row", "cell_col",
            "bias_peak_swe_mm", "notable"]
    st[[c for c in keep if c in st.columns]].to_csv(OUT_TABLE, index=False)

    banner("WRITTEN")
    print(f"  {OUT_TABLE.relative_to(REPO_ROOT)}")
    print(f"  {OUT_FIG.relative_to(REPO_ROOT)}")


if __name__ == "__main__":
    main()
