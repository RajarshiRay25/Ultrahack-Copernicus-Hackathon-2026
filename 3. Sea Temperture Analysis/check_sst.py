"""Check that the downloaded SST file is what humboldt_sst_analysis.py expects.

Run:
    python check_sst.py            (looks for sst_humboldt.nc in this folder)
    python check_sst.py otra_ruta.nc

It prints an inventory of the file and a GO / NO-GO verdict with the reasons.
"""

import os
import sys

import numpy as np
import netCDF4 as nc

# Same analysis box as humboldt_sst_analysis.py
LAT_MIN, LAT_MAX = -30.2, -28.5
LON_MIN, LON_MAX = -73.0, -71.0
MIN_MONTHS = 24

# Grid centres fall on multiples of 0.05 deg, so a request that ends exactly on a
# box edge can come back up to half a cell short. Allow one cell of tolerance.
TOLERANCE = 0.06

PATH = sys.argv[1] if len(sys.argv) > 1 else "sst_humboldt.nc"

print("=" * 60)
if not os.path.exists(PATH):
    sys.exit(f"NO-GO: file not found: {os.path.abspath(PATH)}")

ds = nc.Dataset(PATH)
print(f"File    : {PATH} ({os.path.getsize(PATH) / 1e6:.1f} MB)")
print(f"Variables: {', '.join(ds.variables)}")

name = next((v for v in ("analysed_sst", "sea_surface_temperature", "sst")
             if v in ds.variables), None)
if name is None:
    sys.exit(f"NO-GO: no SST variable. Found: {list(ds.variables)}")

var = ds.variables[name]
lat = np.asarray(ds.variables["latitude"][:], dtype=float)
lon = np.asarray(ds.variables["longitude"][:], dtype=float)
tvar = ds.variables["time"]
dates = nc.num2date(tvar[:], tvar.units)
months = {(d.year, d.month) for d in dates}

print(f"SST     : {name} ({getattr(var, 'units', '?')}), shape {var.shape}")
print(f"Dates   : {len(dates)} steps, {dates[0]:%Y-%m-%d} to {dates[-1]:%Y-%m-%d}, "
      f"{len(months)} distinct months")
print(f"Lat     : {lat.min():.3f} to {lat.max():.3f} (step {abs(lat[1] - lat[0]):.4f})")
print(f"Lon     : {lon.min():.3f} to {lon.max():.3f} (step {abs(lon[1] - lon[0]):.4f})")

# A sample of the field, to confirm units and that the ocean is there
i = len(dates) // 2
slab = np.ma.filled(np.ma.asarray(var[i], dtype=float), np.nan)
if str(getattr(var, "units", "")).strip().upper() in ("K", "KELVIN"):
    slab = slab - 273.15
finite = np.isfinite(slab)
print(f"Sample  : {slab[finite].min():.2f} to {slab[finite].max():.2f} degC, "
      f"{100 * (~finite).mean():.1f}% without data (land)")
print("=" * 60)

checks = [
    ("covers latitude  -30.2 to -28.5",
     lat.min() <= LAT_MIN + TOLERANCE and lat.max() >= LAT_MAX - TOLERANCE),
    ("covers longitude -73.0 to -71.0",
     lon.min() <= LON_MIN + TOLERANCE and lon.max() >= LON_MAX - TOLERANCE),
    (f"at least {MIN_MONTHS} distinct months", len(months) >= MIN_MONTHS),
    ("resolution near 0.05 deg", 0.02 <= abs(lat[1] - lat[0]) <= 0.12),
    ("SST values in a plausible range", finite.any() and 0 < slab[finite].max() < 40),
]
print(f"(coverage tolerance: {TOLERANCE} deg, about one grid cell)")
for label, ok in checks:
    print(f"  [{'OK ' if ok else 'FAIL'}] {label}")

if all(ok for _, ok in checks):
    print("\nGO: the file is ready. Run humboldt_sst_analysis.py")
else:
    print("\nNO-GO: fix the failed checks above before running the analysis"
          " (usually the download box or the date range).")
