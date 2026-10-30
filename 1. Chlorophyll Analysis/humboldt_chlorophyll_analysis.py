"""
Humboldt Archipelago - Where is the ocean most productive?
==========================================================

Copernicus Marine Service Hackathon 2026

QUESTION
    Chile's Humboldt Archipelago multiple-use marine protected area (created 2023,
    ~5,750 km2) has no management plan yet. Which parts of it are the most
    biologically productive, and therefore most important to manage first?

ECOLOGICAL REASONING
    Coastal upwelling brings nutrient-rich deep water to the surface ->
    phytoplankton bloom (seen from space as chlorophyll-a) ->
    zooplankton and small fish (anchovy, sardine) ->
    top predators (Humboldt penguin, seabirds, sea lions, dolphins, whales).
    So chlorophyll-a is a measurable proxy for where the food web is fed.

DATA
    1. Chlorophyll-a: Copernicus Marine product OCEANCOLOUR_GLO_BGC_L4_MY_009_104
       dataset cmems_obs-oc_glo_bgc-plankton_my_l4-multi-4km_P1M
       (GlobColour, satellite, 4 km, monthly), variable CHL [mg/m3]
    2. Protected-area boundary: KMZ from Chile's Ministry of Environment (SIMBIO),
       https://simbio.mma.gob.cl/CbaAP/Details/2082

METHOD (each step is a function below)
    Step 1  Read the official boundary polygon from the KMZ file.
    Step 2  Read chlorophyll, cut it to a box around the area, and average over
            all months in the file (one month or many years both work).
    Step 3  Keep only the satellite pixels whose centre falls inside the boundary.
    Step 4  Measure each pixel's distance to the coast and compare chlorophyll
            in distance bands (coastal strip vs offshore).
    Step 5  Draw the map and save a small table of results.
    Step 6  Seasonal cycle: for each calendar month, average chlorophyll in the
            coastal strip and offshore across all years (needs a multi-year file).

HOW TO RUN
    pip install netCDF4 shapely matplotlib numpy scipy
    python humboldt_chlorophyll_analysis.py  <boundary.kmz>  <chlorophyll.nc>
"""

import sys
import re
import zipfile
import csv

import numpy as np
import netCDF4 as nc
import shapely
from shapely.geometry import Polygon
from shapely.ops import unary_union
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from scipy.ndimage import distance_transform_edt

# Box around the protected area (degrees). Slightly larger than the area
# so the map shows context outside the boundary.
LAT_MIN, LAT_MAX = -30.2, -28.5
LON_MIN, LON_MAX = -72.5, -71.0

# Distance-from-coast bands in km (distance from each pixel to the nearest land pixel)
BANDS = [
    ("Coastal strip (0-10 km)", 0, 10),
    ("Inner zone (10-25 km)",   10, 25),
    ("Middle zone (25-50 km)",  25, 50),
    ("Offshore (over 50 km)",   50, 1e9),
]


# ---------------------------------------------------------------- Step 1
def read_boundary(kmz_path):
    """Read every polygon in the KMZ and merge them into one shape."""
    kml = zipfile.ZipFile(kmz_path).read("doc.kml").decode("utf-8")

    def ring(text):
        coords = re.search(r"<coordinates>(.*?)</coordinates>", text, re.S).group(1)
        return [tuple(map(float, c.split(",")[:2])) for c in coords.split()]

    polygons = []
    for block in re.findall(r"<Polygon>(.*?)</Polygon>", kml, re.S):
        outer = ring(re.search(r"<outerBoundaryIs>(.*?)</outerBoundaryIs>", block, re.S).group(1))
        holes = [ring(h) for h in re.findall(r"<innerBoundaryIs>(.*?)</innerBoundaryIs>", block, re.S)]
        polygons.append(Polygon(outer, holes).buffer(0))  # buffer(0) repairs small geometry errors

    area = unary_union(polygons)
    # Rough area check against the official 574,873 ha (~5,750 km2)
    km2 = area.area * 111.32**2 * np.cos(np.radians(-29.3))
    print(f"[Step 1] Boundary read: {len(polygons)} polygons, ~{km2:,.0f} km2 "
          f"(official: ~5,749 km2)")
    return area


# ---------------------------------------------------------------- Step 2
def read_chlorophyll(nc_path):
    """Cut chlorophyll to the study box and average it over time."""
    ds = nc.Dataset(nc_path)
    lat = ds["latitude"][:]
    lon = ds["longitude"][:]
    lat_sel = (lat >= LAT_MIN) & (lat <= LAT_MAX)
    lon_sel = (lon >= LON_MIN) & (lon <= LON_MAX)

    t = ds["time"]
    dates = nc.num2date(t[:], t.units)
    print(f"[Step 2] Months in file: {len(dates)} ({dates[0]:%Y-%m} to {dates[-1]:%Y-%m})")

    # Read month by month so large files don't fill the memory
    total = None
    count = None
    for i in range(len(dates)):
        month = ds["CHL"][i][np.ix_(lat_sel, lon_sel)]
        month = np.ma.filled(month.astype(float), np.nan)
        valid = np.isfinite(month)
        if total is None:
            total = np.zeros_like(month)
            count = np.zeros_like(month)
        total[valid] += month[valid]
        count[valid] += 1

    with np.errstate(invalid="ignore", divide="ignore"):
        mean_chl = np.where(count > 0, total / count, np.nan)

    period = (f"{dates[0]:%B %Y}" if len(dates) == 1
              else f"{dates[0]:%b %Y} - {dates[-1]:%b %Y} average")
    return lat[lat_sel], lon[lon_sel], mean_chl, period


# ---------------------------------------------------------------- Step 3
def pixels_inside(area, lats, lons):
    """True for each satellite pixel whose centre is inside the protected area."""
    LON, LAT = np.meshgrid(lons, lats)
    inside = shapely.contains_xy(area, LON, LAT)
    print(f"[Step 3] Pixels inside the protected area: {inside.sum()}")
    return LON, LAT, inside


# ---------------------------------------------------------------- Step 4
def distance_to_coast(chl, lats, lons):
    """Distance (km) from each ocean pixel to the nearest land pixel.
    Land = pixels with no chlorophyll value in any month."""
    land = ~np.isfinite(chl)
    km_per_deg_lat = 111.32
    km_per_deg_lon = 111.32 * np.cos(np.radians(np.mean(lats)))
    dy = abs(lats[1] - lats[0]) * km_per_deg_lat
    dx = abs(lons[1] - lons[0]) * km_per_deg_lon
    return distance_transform_edt(~land, sampling=(dy, dx))


def compare_zones(chl, dist_km, inside):
    """Median chlorophyll in each distance band (median resists outliers)."""
    rows = []
    all_vals = chl[inside & np.isfinite(chl)]
    rows.append(("Whole protected area", float(np.median(all_vals)), int(all_vals.size)))
    for name, near, far in BANDS:
        sel = inside & (dist_km > near) & (dist_km <= far) & np.isfinite(chl)
        rows.append((name, float(np.median(chl[sel])), int(sel.sum())))

    print("[Step 4] Median chlorophyll-a (mg/m3):")
    for name, value, n in rows:
        print(f"         {name:28s} {value:5.2f}   ({n} pixels)")
    ratio = rows[1][1] / rows[-1][1]
    print(f"         Coastal strip is {ratio:.1f}x the offshore value")
    return rows, ratio


# ---------------------------------------------------------------- Step 5
def draw_map(area, lats, lons, chl, period, out_png):
    fig, ax = plt.subplots(figsize=(8, 9), dpi=150)
    ax.set_facecolor("#d9d4c7")  # grey = land or no data
    mesh = ax.pcolormesh(lons, lats, chl, cmap="viridis",
                         norm=LogNorm(vmin=0.3, vmax=10), shading="nearest")
    bar = fig.colorbar(mesh, ax=ax, shrink=0.8)
    bar.set_ticks([0.3, 0.5, 1, 2, 5, 10])
    bar.set_ticklabels(["0.3", "0.5", "1", "2", "5", "10"])
    bar.set_label("Chlorophyll-a (mg/m³, log scale)")

    parts = area.geoms if area.geom_type == "MultiPolygon" else [area]
    for p in parts:
        x, y = p.exterior.xy
        ax.plot(x, y, color="white", lw=2.2)
        ax.plot(x, y, color="crimson", lw=1.2)

    # Approximate island positions, added by hand for orientation (not from the data)
    islands = {"Isla Chañaral": (-71.58, -29.03),
               "Isla Damas": (-71.52, -29.23),
               "Isla Choros": (-71.54, -29.26)}
    for name, (x, y) in islands.items():
        ax.plot(x, y, "o", ms=4, mfc="white", mec="black")
        ax.annotate(name, (x, y), xytext=(6, 0), textcoords="offset points", fontsize=8,
                    bbox=dict(boxstyle="round,pad=0.2", fc="white", alpha=0.8, lw=0))

    ax.set_xlim(LON_MIN, LON_MAX)
    ax.set_ylim(LAT_MIN, LAT_MAX)
    ax.set_aspect(1 / np.cos(np.radians(29.3)))
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    ax.set_title(f"Chlorophyll-a, {period}\n"
                 "Humboldt Archipelago protected area (red outline)", fontsize=11)
    ax.text(0.01, 0.01,
            "Data: Copernicus Marine, OCEANCOLOUR_GLO_BGC_L4_MY_009_104 (GlobColour, 4 km)\n"
            "Boundary: MMA Chile / SIMBIO. Grey = land or no data. Island positions approximate.",
            transform=ax.transAxes, fontsize=6, va="bottom",
            bbox=dict(fc="white", alpha=0.8, lw=0))
    plt.tight_layout()
    plt.savefig(out_png)
    print(f"[Step 5] Map saved: {out_png}")


def save_table(rows, out_csv):
    with open(out_csv, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["zone", "median_chlorophyll_mg_m3", "pixels"])
        for name, value, n in rows:
            w.writerow([name, round(value, 3), n])
    print(f"[Step 5] Table saved: {out_csv}")


# ---------------------------------------------------------------- Step 6
def seasonal_cycle(nc_path, lats, lons, dist_km, inside, out_png, out_csv):
    """Median chlorophyll per month in the coastal strip and offshore,
    then averaged by calendar month across all years."""
    ds = nc.Dataset(nc_path)
    lat = ds["latitude"][:]
    lon = ds["longitude"][:]
    lat_sel = (lat >= LAT_MIN) & (lat <= LAT_MAX)
    lon_sel = (lon >= LON_MIN) & (lon <= LON_MAX)
    t = ds["time"]
    dates = nc.num2date(t[:], t.units)
    if len(dates) < 12:
        print("[Step 6] Skipped: needs at least 12 months of data")
        return

    zones = {"Coastal strip (0-10 km)": inside & (dist_km <= 10),
             "Offshore (over 50 km)": inside & (dist_km > 50)}
    by_month = {z: {m: [] for m in range(1, 13)} for z in zones}
    for i, d in enumerate(dates):
        month = np.ma.filled(ds["CHL"][i][np.ix_(lat_sel, lon_sel)].astype(float), np.nan)
        for z, mask in zones.items():
            vals = month[mask & np.isfinite(month)]
            if vals.size:
                by_month[z][d.month].append(np.median(vals))

    names = ["Jan", "Feb", "Mar", "Apr", "May", "Jun",
             "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]
    fig, ax = plt.subplots(figsize=(8, 4.5), dpi=150)
    colours = {"Coastal strip (0-10 km)": "#1b9e77", "Offshore (over 50 km)": "#3b4cc0"}
    rows = []
    for z in zones:
        mean = np.array([np.mean(by_month[z][m]) for m in range(1, 13)])
        lo = np.array([np.min(by_month[z][m]) for m in range(1, 13)])
        hi = np.array([np.max(by_month[z][m]) for m in range(1, 13)])
        ax.plot(names, mean, "o-", color=colours[z], label=z)
        ax.fill_between(names, lo, hi, color=colours[z], alpha=0.15)
        for m in range(12):
            rows.append([z, names[m], round(mean[m], 3), round(lo[m], 3),
                         round(hi[m], 3), len(by_month[z][m + 1])])
    ax.set_ylabel("Chlorophyll-a (mg/m³)")
    ax.set_title(f"Seasonal cycle of chlorophyll-a, {dates[0]:%Y}-{dates[-1]:%Y}\n"
                 "Line = average year; shading = lowest to highest year", fontsize=11)
    ax.legend(frameon=False)
    ax.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(out_png)

    with open(out_csv, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["zone", "month", "mean_mg_m3", "lowest_year", "highest_year", "years"])
        w.writerows(rows)
    print(f"[Step 6] Seasonal chart saved: {out_png}")
    for z in zones:
        mean = [np.mean(by_month[z][m]) for m in range(1, 13)]
        print(f"         {z}: peak {names[int(np.argmax(mean))]} ({max(mean):.2f}), "
              f"lowest {names[int(np.argmin(mean))]} ({min(mean):.2f})")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print(__doc__)
        sys.exit(1)
    kmz_path, nc_path = sys.argv[1], sys.argv[2]

    area = read_boundary(kmz_path)
    lats, lons, chl, period = read_chlorophyll(nc_path)
    LON, LAT, inside = pixels_inside(area, lats, lons)
    dist_km = distance_to_coast(chl, lats, lons)
    rows, ratio = compare_zones(chl, dist_km, inside)
    draw_map(area, lats, lons, chl, period, "humboldt_chlorophyll_map.png")
    save_table(rows, "humboldt_chlorophyll_by_zone.csv")
    seasonal_cycle(nc_path, lats, lons, dist_km, inside,
                   "humboldt_chlorophyll_seasonal.png", "humboldt_chlorophyll_seasonal.csv")
