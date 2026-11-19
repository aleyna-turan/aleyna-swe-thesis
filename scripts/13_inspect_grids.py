"""
13_inspect_grids.py

Read-only. Prints the structure of every NetCDF file in data/ so that the
subgrid-elevation script can be written against the real variable names and
the real grid, instead of guessing at them.

Writes nothing. Changes nothing.

What it is looking for
----------------------
1. In dem_hicam_3km_eu.nc   : the elevation variable, its coordinates, the
                              cell size, and the area it covers.
2. In the mHM output file   : the grid the model runs on, so each station
                              can be matched to the cell it falls in.

Why this has to happen first
----------------------------
The test is: for each station, work out the mean elevation of the 12 km mHM
cell it sits in, subtract it from the station's own elevation, and see
whether the model's bias scales with that difference. Doing that needs the
exact coordinate names and cell edges of both files. Those cannot be
guessed; a wrong guess produces a script that fails on your machine or,
worse, one that silently matches stations to the wrong cells.

Run:
    cd /Users/aleynaturan/Desktop/a/aleyna-swe-thesis
    source .venv/bin/activate
    python scripts/13_inspect_grids.py
"""

from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parent.parent
DATA_DIR = REPO_ROOT / "data"

RULE = "=" * 72
THIN = "-" * 72


def human(n):
    for unit in ("B", "KB", "MB", "GB"):
        if n < 1024 or unit == "GB":
            return f"{n:.1f} {unit}"
        n /= 1024


def describe_axis(name, values):
    """Print the span, size and spacing of one coordinate axis."""
    v = np.asarray(values).ravel()
    if v.size == 0:
        print(f"    {name:<14} empty")
        return
    try:
        v = v.astype(float)
    except (TypeError, ValueError):
        print(f"    {name:<14} size {v.size}, first {v[0]!r}, last {v[-1]!r}")
        return
    line = f"    {name:<14} size {v.size:>6}   {v.min():>12.4f} to {v.max():>12.4f}"
    if v.size > 1:
        d = np.diff(v)
        step = np.median(d)
        even = "even" if np.allclose(d, step, rtol=1e-3, atol=1e-9) else "UNEVEN"
        line += f"   step {step:>9.5f} ({even})"
    print(line)


def inspect_with_xarray(path):
    import xarray as xr
    with xr.open_dataset(path, decode_times=False) as ds:
        print("\n  dimensions:")
        for k, v in ds.sizes.items():
            print(f"    {k:<14} {v}")

        print("\n  coordinates:")
        if not ds.coords:
            print("    (none declared as coordinates)")
        for k in ds.coords:
            describe_axis(k, ds[k].values)
            units = ds[k].attrs.get("units")
            if units:
                print(f"      units: {units}")

        print("\n  data variables:")
        for k in ds.data_vars:
            da = ds[k]
            attrs = da.attrs
            bits = [f"dims {tuple(da.dims)}", f"dtype {da.dtype}"]
            if "units" in attrs:
                bits.append(f"units {attrs['units']}")
            print(f"    {k:<20} {'  '.join(bits)}")
            for a in ("long_name", "standard_name", "grid_mapping",
                      "_FillValue", "missing_value", "scale_factor"):
                if a in attrs:
                    print(f"      {a}: {attrs[a]}")
            # a value range, for the small variables only
            if da.size and da.size < 60_000_000 and np.issubdtype(da.dtype, np.number):
                try:
                    arr = da.values
                    finite = arr[np.isfinite(arr)]
                    if finite.size:
                        print(f"      values: {finite.min():.2f} to {finite.max():.2f}, "
                              f"{100 * finite.size / arr.size:.1f}% finite")
                except Exception as e:
                    print(f"      [could not read values: {e}]")

        if ds.attrs:
            print("\n  global attributes:")
            for k, v in list(ds.attrs.items())[:15]:
                s = str(v)
                print(f"    {k}: {s[:110]}{'...' if len(s) > 110 else ''}")


def inspect_with_netcdf4(path):
    from netCDF4 import Dataset
    with Dataset(path) as ds:
        print("\n  dimensions:")
        for k, v in ds.dimensions.items():
            print(f"    {k:<14} {len(v)}")
        print("\n  variables:")
        for k, v in ds.variables.items():
            print(f"    {k:<20} dims {v.dimensions}  dtype {v.dtype}")
            for a in v.ncattrs():
                if a in ("units", "long_name", "standard_name", "grid_mapping"):
                    print(f"      {a}: {v.getncattr(a)}")
            if v.ndim == 1 and len(v) > 0:
                describe_axis(k, v[:])


def main():
    if not DATA_DIR.exists():
        raise SystemExit(f"No data folder at {DATA_DIR}")

    files = sorted(list(DATA_DIR.rglob("*.nc")) + list(DATA_DIR.rglob("*.nc4")))
    if not files:
        print(f"No NetCDF files found under {DATA_DIR}")
        print("Files that are there:")
        for f in sorted(DATA_DIR.iterdir())[:40]:
            print(f"  {f.name}")
        return

    print(f"{len(files)} NetCDF file(s) under {DATA_DIR.relative_to(REPO_ROOT)}\n")

    for path in files:
        print(RULE)
        print(f"{path.relative_to(REPO_ROOT)}   ({human(path.stat().st_size)})")
        print(RULE)
        try:
            inspect_with_xarray(path)
        except ImportError:
            print("  [xarray not installed, falling back to netCDF4]")
            try:
                inspect_with_netcdf4(path)
            except Exception as e:
                print(f"  [could not open: {e}]")
        except Exception as e:
            print(f"  [xarray could not open it: {e}]")
            try:
                inspect_with_netcdf4(path)
            except Exception as e2:
                print(f"  [netCDF4 could not open it either: {e2}]")
        print()

    print(THIN)
    print("Paste this whole output back. From it I can write the subgrid")
    print("script against the real variable names and the real grid.")
    print(THIN)


if __name__ == "__main__":
    main()
