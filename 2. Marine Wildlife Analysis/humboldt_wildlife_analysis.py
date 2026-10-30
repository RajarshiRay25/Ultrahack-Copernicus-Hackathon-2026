"""
Humboldt Archipelago - Do the animals gather where the food is?
===============================================================

Copernicus Marine Service Hackathon 2026

QUESTION
    Chlorophyll (the base of the food web) is highest in the coastal strip and
    peaks in spring. Are penguins, sea lions, dolphins and whales also
    concentrated in the coastal strip, and when are they seen?

DATA
    1. Animal sightings: GBIF occurrence download (Simple CSV, tab-separated),
       species: Humboldt penguin, South American sea lion, bottlenose dolphin,
       blue whale, fin whale, humpback whale.
    2. Chlorophyll and boundary: same files as humboldt_chlorophyll_analysis.py
       (this script reuses its functions, so keep both files in one folder).

METHOD
    Step 1  Read sightings; keep only "PRESENT" records with coordinates, 2015 onward
            (to match the chlorophyll period).
    Step 2  Give each sighting the distance-to-coast of the nearest satellite pixel
            and check whether it is inside the protected area.
    Step 3  Sightings per 100 km2 in each distance zone (dividing by area is
            essential, because the offshore zones are much bigger).
    Step 4  Map: sightings on top of the 11-year chlorophyll map.
    Step 5  Seasonal chart: sightings per month vs coastal chlorophyll per month.

IMPORTANT LIMITATION
    Sightings show where and when PEOPLE looked, not only where animals are.
    Most observers are on the coast, near the islands, and in the summer holidays.

HOW TO RUN
    python humboldt_wildlife_analysis.py <boundary.kmz> <chlorophyll.nc> <gbif.csv>
"""

import sys
import csv
from collections import Counter, defaultdict

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import shapely

import humboldt_chlorophyll_analysis as chl_tools

START_YEAR = 2015

GROUPS = {
    "Spheniscus humboldti": ("Humboldt penguin", "#111111", "o"),
    "Otaria byronia": ("Sea lion", "#e66101", "s"),
    "Otaria flavescens": ("Sea lion", "#e66101", "s"),
    "Tursiops truncatus": ("Bottlenose dolphin", "#d01c8b", "^"),
    "Balaenoptera musculus": ("Whales", "#ffffff", "D"),
    "Balaenoptera physalus": ("Whales", "#ffffff", "D"),
    "Megaptera novaeangliae": ("Whales", "#ffffff", "D"),
}

ZONES = [
    ("Coast & islands (0-10 km)", -1, 10),
    ("Inner zone (10-25 km)", 10, 25),
    ("Middle zone (25-50 km)", 25, 50),
    ("Offshore (over 50 km)", 50, 1e9),
]


# ---------------------------------------------------------------- Step 1
def read_sightings(csv_path):
    rows = []
    with open(csv_path, encoding="utf-8") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if r["occurrenceStatus"] != "PRESENT":
                continue  # "ABSENT" means the animal was looked for and not found
            if not r["decimalLatitude"] or not r["year"] or int(r["year"]) < START_YEAR:
                continue
            if r["species"] not in GROUPS:
                continue
            rows.append({
                "group": GROUPS[r["species"]][0],
                "lat": float(r["decimalLatitude"]),
                "lon": float(r["decimalLongitude"]),
                "month": int(r["month"]) if r["month"] else None,
                "source": r["institutionCode"],
            })
    print(f"[Step 1] Sightings kept ({START_YEAR} onward, present, with coordinates): {len(rows)}")
    for g, n in Counter(r["group"] for r in rows).most_common():
        print(f"         {g:20s} {n}")
    return rows


# ---------------------------------------------------------------- Step 2
def assign_zones(rows, area, lats, lons, dist_km):
    for r in rows:
        i = int(np.argmin(np.abs(lats - r["lat"])))
        j = int(np.argmin(np.abs(lons - r["lon"])))
        r["dist_km"] = float(dist_km[i, j])  # 0 = on a land/island pixel (shore observation)
        r["inside"] = bool(shapely.contains_xy(area, r["lon"], r["lat"])) or (
            r["dist_km"] == 0 and area.distance(shapely.Point(r["lon"], r["lat"])) < 0.02)
    n_in = sum(r["inside"] for r in rows)
    print(f"[Step 2] Sightings inside or on the edge of the protected area: {n_in}")
    return rows


# ---------------------------------------------------------------- Step 3
def sightings_per_zone(rows, inside_px, dist_km, lats, lons, out_csv):
    km_lat = abs(lats[1] - lats[0]) * 111.32
    km_lon = abs(lons[1] - lons[0]) * 111.32 * np.cos(np.radians(np.mean(lats)))
    pixel_km2 = km_lat * km_lon
    table = []
    print("[Step 3] Sightings per 100 km2 inside the protected area:")
    for name, near, far in ZONES:
        if near < 0:
            # coastal zone: ocean pixels within 10 km (land pixels add no sea area)
            n_px = int((inside_px & (dist_km > 0) & (dist_km <= far)).sum())
        else:
            n_px = int((inside_px & (dist_km > near) & (dist_km <= far)).sum())
        area_km2 = n_px * pixel_km2
        n = sum(1 for r in rows if r["inside"] and near < r["dist_km"] <= far)
        density = 100 * n / area_km2 if area_km2 else float("nan")
        table.append((name, n, round(area_km2), density))
        print(f"         {name:28s} {n:5d} sightings  {area_km2:6.0f} km2  {density:6.2f}")
    with open(out_csv, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["zone", "sightings", "area_km2", "sightings_per_100km2"])
        for name, n, a, d in table:
            w.writerow([name, n, a, round(d, 2)])
    return table


# ---------------------------------------------------------------- Step 4
def draw_map(area, lats, lons, chl, rows, out_png):
    fig, ax = plt.subplots(figsize=(8, 9), dpi=150)
    ax.set_facecolor("#d9d4c7")
    mesh = ax.pcolormesh(lons, lats, chl, cmap="viridis",
                         norm=LogNorm(vmin=0.3, vmax=10), shading="nearest", alpha=0.85)
    bar = fig.colorbar(mesh, ax=ax, shrink=0.8)
    bar.set_ticks([0.3, 0.5, 1, 2, 5, 10])
    bar.set_ticklabels(["0.3", "0.5", "1", "2", "5", "10"])
    bar.set_label("Chlorophyll-a, 11-year average (mg/m³)")
    parts = area.geoms if area.geom_type == "MultiPolygon" else [area]
    for p in parts:
        x, y = p.exterior.xy
        ax.plot(x, y, color="white", lw=2.2)
        ax.plot(x, y, color="crimson", lw=1.2)
    done = set()
    for label, colour, marker in GROUPS.values():
        if label in done:
            continue
        done.add(label)
        pts = [r for r in rows if r["group"] == label]
        ax.scatter([r["lon"] for r in pts], [r["lat"] for r in pts], s=14, c=colour,
                   marker=marker, edgecolors="black", linewidths=0.4, alpha=0.8,
                   label=f"{label} ({len(pts)})")
    ax.legend(loc="lower left", fontsize=8, framealpha=0.9)
    ax.set_xlim(chl_tools.LON_MIN, chl_tools.LON_MAX)
    ax.set_ylim(chl_tools.LAT_MIN, chl_tools.LAT_MAX)
    ax.set_aspect(1 / np.cos(np.radians(29.3)))
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    ax.set_title(f"Animal sightings ({START_YEAR}-2026) on the chlorophyll map\n"
                 "Humboldt Archipelago protected area (red outline)", fontsize=11)
    ax.text(0.99, 0.01, "Sightings: GBIF (eBird, iNaturalist, Happywhale and others)\n"
            "Chlorophyll: Copernicus Marine OCEANCOLOUR_GLO_BGC_L4_MY_009_104",
            transform=ax.transAxes, fontsize=6, ha="right", va="bottom",
            bbox=dict(fc="white", alpha=0.8, lw=0))
    plt.tight_layout()
    plt.savefig(out_png)
    print(f"[Step 4] Map saved: {out_png}")


# ---------------------------------------------------------------- Step 5
def seasonal_comparison(rows, seasonal_csv, out_png):
    """Sightings per month (all years together) next to coastal chlorophyll."""
    chl_by_month = {}
    with open(seasonal_csv, encoding="utf-8") as f:
        for r in csv.DictReader(f):
            if r["zone"].startswith("Coastal"):
                chl_by_month[r["month"]] = float(r["mean_mg_m3"])
    names = ["Jan", "Feb", "Mar", "Apr", "May", "Jun",
             "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]
    counts = defaultdict(Counter)
    for r in rows:
        if r["month"]:
            counts[r["group"]][r["month"]] += 1

    fig, axes = plt.subplots(2, 1, figsize=(8, 7), dpi=150, sharex=True,
                             gridspec_kw={"height_ratios": [1, 2]})
    axes[0].plot(names, [chl_by_month[m] for m in names], "o-", color="#1b9e77")
    axes[0].set_ylabel("Chlorophyll-a\n(mg/m³)")
    axes[0].set_title("Food supply (coastal chlorophyll) and animal sightings by month", fontsize=11)
    axes[0].grid(alpha=0.3)
    style = {"Humboldt penguin": "#111111", "Sea lion": "#e66101",
             "Bottlenose dolphin": "#d01c8b", "Whales": "#2c7fb8"}
    for g, colour in style.items():
        axes[1].plot(names, [counts[g][m] for m in range(1, 13)], "o-", color=colour, label=g)
    axes[1].set_ylabel("Number of sightings")
    axes[1].legend(frameon=False, fontsize=8)
    axes[1].grid(alpha=0.3)
    axes[1].axvspan(-0.4, 1.5, color="grey", alpha=0.12)
    axes[1].axvspan(10.6, 11.4, color="grey", alpha=0.12)
    axes[1].text(0.0, axes[1].get_ylim()[1] * 0.92, "summer holidays\n(more observers)",
                 fontsize=7, color="dimgrey")
    plt.tight_layout()
    plt.savefig(out_png)
    print(f"[Step 5] Seasonal comparison saved: {out_png}")
    for g in style:
        peak = max(range(1, 13), key=lambda m: counts[g][m])
        print(f"         {g:20s} most sightings in {names[peak - 1]}")


if __name__ == "__main__":
    if len(sys.argv) != 4:
        print(__doc__)
        sys.exit(1)
    kmz_path, nc_path, gbif_path = sys.argv[1:4]

    area = chl_tools.read_boundary(kmz_path)
    lats, lons, chl, period = chl_tools.read_chlorophyll(nc_path)
    LON, LAT, inside_px = chl_tools.pixels_inside(area, lats, lons)
    dist_km = chl_tools.distance_to_coast(chl, lats, lons)

    rows = read_sightings(gbif_path)
    rows = assign_zones(rows, area, lats, lons, dist_km)
    sightings_per_zone(rows, inside_px, dist_km, lats, lons, "humboldt_wildlife_by_zone.csv")
    draw_map(area, lats, lons, chl, rows, "humboldt_wildlife_map.png")
    seasonal_comparison(rows, "humboldt_chlorophyll_seasonal.csv", "humboldt_wildlife_seasonal.png")
