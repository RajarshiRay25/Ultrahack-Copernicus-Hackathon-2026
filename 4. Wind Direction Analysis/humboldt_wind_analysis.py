"""
Humboldt Archipelago - Is the wind driving the upwelling?
=========================================================

Copernicus Marine Service Hackathon 2026

QUESTION
    Winds blowing toward the north, along this north-south coast, push surface
    water offshore (Earth's rotation turns it to the left in the Southern
    Hemisphere). Deep, cold, nutrient-rich water rises to replace it. Do the
    winds here blow this way, how strong are they, and in which months and
    latitudes are they strongest?

DATA
    1. Wind: Copernicus Marine product WIND_GLO_PHY_CLIMATE_L4_MY_012_003
       (monthly, 0.25 degree), variables eastward_wind and northward_wind [m/s]
    2. Protected-area boundary: KMZ from Chile's Ministry of Environment (SIMBIO),
       https://simbio.mma.gob.cl/CbaAP/Details/2082

WIND CATEGORIES (northward wind, i.e. the part that drives upwelling)
    Below 0 m/s   "Downwelling"  wind toward the south, pushes water toward the coast
    0 to 5 m/s    "Weak"         little upwelling
    5 to 7 m/s    "Moderate"     steady upwelling
    7 m/s or more "Strong"       intense upwelling
    These thresholds are a transparent, adjustable choice (see CATEGORIES below),
    not an official standard. Changing them changes the labels, not the data.

METHOD
    Step 1  Read the official boundary polygon from the KMZ file.
    Step 2  Read wind, cut to a box around the area.
    Step 3  Location map: where the wind squares are, coloured by category.
            If the file has only one column of squares (a north-south line),
            the analysis runs along that line ("transect mode").
    Step 4  Latitude x month chart: wind category for each band and month.
    Step 5  Seasonal cycle of the northward wind at the protected area's latitudes.

HOW TO RUN
    pip install netCDF4 shapely matplotlib numpy
    python humboldt_wind_analysis.py <boundary.kmz> <wind.nc>
"""

import sys
import re
import csv
import zipfile
import warnings

import numpy as np
import netCDF4 as nc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm
from matplotlib.patches import Patch
from shapely.geometry import box, Polygon
from shapely.ops import unary_union

warnings.simplefilter("ignore", RuntimeWarning)  # land squares have no wind values

LAT_MIN, LAT_MAX = -31.0, -27.5
LON_MIN, LON_MAX = -74.0, -71.0
MONTHS = ["Jan", "Feb", "Mar", "Apr", "May", "Jun",
          "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]

CATEGORIES = [  # (name, lower limit m/s, upper limit m/s, colour)
    ("Downwelling", -99, 0, "#4575b4"),
    ("Weak", 0, 5, "#fee090"),
    ("Moderate", 5, 7, "#fc8d59"),
    ("Strong", 7, 99, "#d73027"),
]
CAT_CMAP = ListedColormap([c[3] for c in CATEGORIES])
CAT_NORM = BoundaryNorm([CATEGORIES[0][1]] + [c[2] for c in CATEGORIES], CAT_CMAP.N)


def category(v):
    for name, lo, hi, _ in CATEGORIES:
        if lo <= v < hi:
            return name
    return "n/a"


def find_name(ds, must_contain, must_not_contain=()):
    """Find a variable by part of its name (names differ slightly between products)."""
    for name in ds.variables:
        low = name.lower()
        if all(w in low for w in must_contain) and not any(w in low for w in must_not_contain):
            return name
    raise KeyError(f"No variable containing {must_contain} in file: {list(ds.variables)}")


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
        polygons.append(Polygon(outer, holes).buffer(0))  # buffer(0) repairs small errors
    area = unary_union(polygons)
    km2 = area.area * 111.32**2 * np.cos(np.radians(-29.3))
    print(f"[Step 1] Boundary read: {len(polygons)} polygons, ~{km2:,.0f} km2 "
          f"(official: ~5,749 km2)")
    return area


# ---------------------------------------------------------------- Step 2
def read_wind(nc_path):
    ds = nc.Dataset(nc_path)
    lat_name = "latitude" if "latitude" in ds.variables else "lat"
    lon_name = "longitude" if "longitude" in ds.variables else "lon"
    u_name = find_name(ds, ["eastward", "wind"], ["stress", "error", "sdd", "uncertainty"])
    v_name = find_name(ds, ["northward", "wind"], ["stress", "error", "sdd", "uncertainty"])
    lat = ds[lat_name][:]
    lon = ds[lon_name][:]
    lat_sel = (lat >= LAT_MIN) & (lat <= LAT_MAX)
    lon_sel = (lon >= LON_MIN - 0.5) & (lon <= LON_MAX)
    t = ds["time"]
    dates = nc.num2date(t[:], t.units)
    u = np.ma.filled(ds[u_name][:][:, lat_sel][:, :, lon_sel].astype(float), np.nan)
    v = np.ma.filled(ds[v_name][:][:, lat_sel][:, :, lon_sel].astype(float), np.nan)
    lats, lons = np.array(lat[lat_sel]), np.array(lon[lon_sel])
    print(f"[Step 2] Wind variables: {u_name}, {v_name}")
    print(f"[Step 2] Months: {len(dates)} ({dates[0]:%Y-%m} to {dates[-1]:%Y-%m}); "
          f"squares: {len(lats)} latitudes x {len(lons)} longitudes")
    return lats, lons, dates, u, v


def area_rows(area, lats, lons):
    """Which wind squares to use for the protected area.
    Normal mode: squares that touch the area. Transect mode (one column of
    squares that does not touch the area): squares at the same latitudes."""
    half = abs(lats[1] - lats[0]) / 2
    touch = np.zeros((len(lats), len(lons)), dtype=bool)
    for i, la in enumerate(lats):
        for j, lo in enumerate(lons):
            touch[i, j] = box(lo - half, la - half, lo + half, la + half).intersects(area)
    if touch.any():
        return touch, False
    s, n = area.bounds[1], area.bounds[3]
    rows = (lats >= s - half) & (lats <= n + half)
    sel = np.zeros_like(touch)
    sel[rows, :] = True
    km = (area.bounds[0] - lons.max()) * 111.32 * np.cos(np.radians(29.3))
    print(f"[Step 3] TRANSECT MODE: wind data is one line of squares about {km:.0f} km "
          f"west of the protected area; using the {rows.sum()} squares at the area's latitudes")
    return sel, True


# ---------------------------------------------------------------- Step 3
def draw_location_map(area, lats, lons, v_mean, transect, out_png):
    fig, ax = plt.subplots(figsize=(7, 8), dpi=150)
    parts = area.geoms if area.geom_type == "MultiPolygon" else [area]
    for p in parts:
        x, y = p.exterior.xy
        ax.fill(x, y, color="#fde0dd", zorder=1)
        ax.plot(x, y, color="crimson", lw=1.2, zorder=2)
    half = abs(lats[1] - lats[0]) / 2
    for i, la in enumerate(lats):
        for j, lo in enumerate(lons):
            val = v_mean[i, j]
            if not np.isfinite(val):
                continue
            colour = CAT_CMAP(CAT_NORM(val))
            ax.add_patch(plt.Rectangle((lo - half, la - half), 2 * half, 2 * half,
                                       fc=colour, ec="white", lw=0.8, zorder=3))
            ax.text(lo + half + 0.03, la, f"{val:.1f} m/s", ha="left", va="center",
                    fontsize=6, zorder=4)
    LON, LAT = np.meshgrid(lons, lats)
    ax.quiver(LON, LAT, np.nanmean(U_ALL, 0), np.nanmean(V_ALL, 0), color="black",
              scale=150, width=0.004, zorder=5)
    ax.annotate("Protected area", xy=(-71.7, -29.3), fontsize=9, color="crimson", ha="center")
    if transect:
        ax.annotate("", xy=(-72.05, -29.3), xytext=(lons.max() + 0.2, -29.3),
                    arrowprops=dict(arrowstyle="<->", color="grey"))
        km = (area.bounds[0] - lons.max()) * 111.32 * np.cos(np.radians(29.3))
        ax.text((lons.max() - 72.05) / 2 + 0.1, -29.22, f"~{km:.0f} km", fontsize=8,
                color="grey", ha="center")
    ax.set_xlim(LON_MIN, LON_MAX)
    ax.set_ylim(LAT_MIN, LAT_MAX)
    ax.set_aspect(1 / np.cos(np.radians(29.3)))
    ax.set_facecolor("#eef5fb")
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    note = ("Wind data in this file covers only the coloured squares"
            if transect else "Arrows = wind direction")
    ax.set_title("Average northward (upwelling) wind by wind square, 2015-2026\n" + note,
                 fontsize=10)
    labels = {"Downwelling": "Downwelling (below 0 m/s)", "Weak": "Weak (0-5 m/s)",
              "Moderate": "Moderate (5-7 m/s)", "Strong": "Strong (7 m/s or more)"}
    ax.legend(handles=[Patch(fc=c[3], label=labels[c[0]]) for c in CATEGORIES],
              loc="lower right", fontsize=7, title="Upwelling wind category", title_fontsize=8)
    ax.text(0.01, 0.01, "Data: Copernicus Marine WIND_GLO_PHY_CLIMATE_L4_MY_012_003\n"
            "Boundary: MMA Chile / SIMBIO", transform=ax.transAxes, fontsize=6)
    plt.tight_layout()
    plt.savefig(out_png)
    print(f"[Step 3] Location map saved: {out_png}")


# ---------------------------------------------------------------- Step 4
def latitude_month_chart(lats, dates, v, area, out_png, out_csv):
    """Average northward wind for every latitude band and calendar month."""
    line = np.nanmean(v, axis=2)  # average across longitudes -> time x latitude
    grid = np.array([[np.nanmean([line[i, r] for i, d in enumerate(dates) if d.month == m])
                      for m in range(1, 13)] for r in range(len(lats))])
    fig, ax = plt.subplots(figsize=(8, 5.5), dpi=150)
    ax.pcolormesh(np.arange(12), lats, grid, cmap=CAT_CMAP, norm=CAT_NORM, shading="nearest")
    for r, la in enumerate(lats):
        for m in range(12):
            ax.text(m, la, f"{grid[r, m]:.1f}", ha="center", va="center", fontsize=6)
    s, n = area.bounds[1], area.bounds[3]
    ax.axhspan(s, n, fill=False, ec="crimson", lw=2)
    ax.text(11.6, (s + n) / 2, "protected\narea\nlatitudes", color="crimson", fontsize=7,
            va="center")
    ax.set_xticks(range(12))
    ax.set_xticklabels(MONTHS)
    ax.set_xlim(-0.5, 12.6)
    ax.set_ylabel("Latitude")
    ax.set_title("Upwelling wind category by latitude and month (2015-2026 average)\n"
                 "Numbers = northward wind in m/s", fontsize=10)
    ax.legend(handles=[Patch(fc=c[3], label=c[0]) for c in CATEGORIES], loc="upper center",
              bbox_to_anchor=(0.5, -0.08), ncol=4, fontsize=8, frameon=False)
    plt.tight_layout()
    plt.savefig(out_png)
    with open(out_csv, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["latitude"] + [f"{m}_m_s" for m in MONTHS] + [f"{m}_category" for m in MONTHS])
        for r, la in enumerate(lats):
            w.writerow([la] + [round(x, 2) for x in grid[r]] + [category(x) for x in grid[r]])
    print(f"[Step 4] Latitude x month chart saved: {out_png}")
    return grid


# ---------------------------------------------------------------- Step 5
def seasonal_cycle(dates, u, v, sel, out_png, out_csv):
    series_v = np.array([np.nanmean(v[i][sel]) for i in range(len(dates))])
    series_u = np.array([np.nanmean(u[i][sel]) for i in range(len(dates))])
    by_m = {m: series_v[[d.month == m for d in dates]] for m in range(1, 13)}
    v_mean = np.array([by_m[m].mean() for m in range(1, 13)])
    v_lo = np.array([by_m[m].min() for m in range(1, 13)])
    v_hi = np.array([by_m[m].max() for m in range(1, 13)])
    u_mean = np.array([series_u[[d.month == m for d in dates]].mean() for m in range(1, 13)])

    fig, ax = plt.subplots(figsize=(8, 4.2), dpi=150)
    for name, lo, hi, colour in CATEGORIES:
        ax.axhspan(max(lo, -1), min(hi, 10), color=colour, alpha=0.18, lw=0)
        if 0 <= lo < 9:
            ax.text(11.45, lo + 0.1, name, fontsize=7, color="dimgrey")
    ax.plot(MONTHS, v_mean, "o-", color="black", label="Average year")
    ax.fill_between(MONTHS, v_lo, v_hi, color="grey", alpha=0.25, label="Weakest to strongest year")
    ax.set_ylim(min(0, v_lo.min() - 0.5), max(9, v_hi.max() + 0.5))
    ax.set_ylabel("Northward wind (m/s)")
    ax.set_title(f"Upwelling wind at the protected area's latitudes, "
                 f"{dates[0]:%Y}-{dates[-1]:%Y}\n"
                 "Above 0 = wind toward the north (pushes surface water offshore)", fontsize=10)
    ax.legend(loc="lower left", fontsize=8, frameon=False)
    ax.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(out_png)

    with open(out_csv, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["month", "northward_wind_mean_m_s", "category", "weakest_year",
                    "strongest_year", "eastward_wind_mean_m_s", "years"])
        for m in range(12):
            w.writerow([MONTHS[m], round(v_mean[m], 2), category(v_mean[m]), round(v_lo[m], 2),
                        round(v_hi[m], 2), round(u_mean[m], 2), len(by_m[m + 1])])

    print(f"[Step 5] Seasonal chart saved: {out_png}")
    for m in range(12):
        print(f"         {MONTHS[m]}: {v_mean[m]:.2f} m/s  {category(v_mean[m])}")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print(__doc__)
        sys.exit(1)
    kmz_path, nc_path = sys.argv[1:3]

    area = read_boundary(kmz_path)
    lats, lons, dates, U_ALL, V_ALL = read_wind(nc_path)
    v_map = np.nanmean(V_ALL, axis=0)
    sel, transect = area_rows(area, lats, lons)
    draw_location_map(area, lats, lons, v_map, transect, "humboldt_wind_map.png")
    latitude_month_chart(lats, dates, V_ALL, area, "humboldt_wind_latitude_month.png",
                         "humboldt_wind_latitude_month.csv")
    seasonal_cycle(dates, U_ALL, V_ALL, sel, "humboldt_wind_seasonal.png",
                   "humboldt_wind_seasonal.csv")
