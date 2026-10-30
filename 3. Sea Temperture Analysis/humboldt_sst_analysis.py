"""
Humboldt Archipelago - is the water cold, and where?
====================================================

Copernicus Marine Service Hackathon 2026

QUESTION
    Inside the Humboldt Archipelago multiple-use marine protected area, which
    sectors show the strongest upwelling signature (cold water close to the
    coast)? That is where the food web is being fed, so it is where limited
    management resources should go first.

ECOLOGICAL REASONING
    The alongshore wind drives an offshore Ekman transport; the surface water
    carried offshore is replaced by cold, nutrient-rich water from below.
    A cold anomaly near the coast, relative to open water at the same latitude,
    is therefore the fingerprint of upwelling.

DATA
    1. Sea surface temperature: Copernicus Marine product
       SST_GLO_SST_L4_REP_OBSERVATIONS_010_011 (OSTIA, satellite, 0.05 deg,
       daily, foundation SST), variable analysed_sst [K].
    2. Protected-area boundary: KMZ from Chile's Ministry of Environment (SIMBIO),
       https://simbio.mma.gob.cl/CbaAP/Details/2082

METHOD (each step is a function below)
    Step 1  Read the official boundary polygon from the KMZ file.
    Step 2  Read the SST grid, convert kelvin -> degC, build a land mask.
    Step 3  Measure the distance from every water pixel to the nearest land pixel.
    Step 4  Keep the pixels whose centre falls inside the boundary; define
            sectors by latitude.
    Step 5  One streaming pass over every date: monthly mean per pixel, then
            medians per zone -> seasonal cycle, coastal-offshore gradient series,
            and the time-mean fields for the maps.
    Step 6  Draw the maps and charts, save the tables.

WHY THE OFFSHORE REFERENCE STARTS AT 50 km
    The protected area is only about 55 km wide. Restricting the offshore band to
    pixels inside the polygon leaves a one-cell sliver at the western edge, which
    is not open water. Here the reference band is 50-150 km from the coast, which
    lies outside the boundary almost everywhere, so it is a real open-ocean
    comparison. Use the same definition for chlorophyll to keep the two datasets
    comparable.

HOW TO RUN
    pip install netCDF4 shapely matplotlib numpy scipy
    python humboldt_sst_analysis.py  <boundary.kmz>  <sst.nc>
"""

import sys
import re
import zipfile
import csv

MISSING_PACKAGES_HELP = """
Missing dependencies: {missing}

Install them in the interpreter you are running this script with:

    python -m pip install netCDF4 shapely scipy matplotlib numpy

With conda:

    conda install -c conda-forge netcdf4 shapely scipy matplotlib

In VS Code, check that the selected interpreter is the one where you installed
them: Ctrl+Shift+P -> "Python: Select Interpreter".
"""

try:
    import numpy as np
    import netCDF4 as nc
    import shapely
    from shapely.geometry import Polygon
    from shapely.ops import unary_union
    import matplotlib
    from matplotlib.colors import Normalize
    from scipy.ndimage import distance_transform_edt
except ImportError as exc:
    raise SystemExit(MISSING_PACKAGES_HELP.format(missing=exc.name))

matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Analysis box (degrees). Extends further west than the chlorophyll box so the
# 50-150 km offshore reference band is fully covered.
LAT_MIN, LAT_MAX = -30.2, -28.5
LON_MIN, LON_MAX = -73.0, -71.0

# Maps are drawn on a narrower window: beyond ~150 km the offshore reference band
# ends, so the far west would show a contrast that is not the one we compare.
VIEW_LON_MIN = -72.35

# Distance-from-coast bands (km to the nearest land pixel).
# Last field: inside the protected boundary only? The offshore reference must be
# measured outside the boundary, otherwise the band degenerates into the sliver
# along the western edge of the area.
BANDS = [
    ("Coastal strip (0-10 km)", 0, 10, True),
    ("Inner zone (10-25 km)",   10, 25, True),
    ("Middle zone (25-50 km)",  25, 50, True),
    ("Offshore reference (50-150 km)", 50, 150, False),
]

# The coastal band used for the upwelling index
COAST_KM = 10
OFFSHORE_KM = (50, 150)

# Sectors of the protected area, split by latitude (approximate landmark in brackets)
SECTORS = [
    ("Sector 1: 28.89-29.10 S (Isla Chanaral)",   -29.10, -28.88),
    ("Sector 2: 29.10-29.35 S (Punta de Choros)", -29.35, -29.10),
    ("Sector 3: 29.35-29.55 S (Chanaral de Aceituno)", -29.55, -29.35),
    ("Sector 4: 29.55-29.76 S (southern sector)", -29.76, -29.55),
]

SST_VARIANTS = ("analysed_sst", "sea_surface_temperature", "sst")


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
    # Rough area check against the official 574,873 ha (~5,749 km2)
    km2 = area.area * 111.32**2 * np.cos(np.radians(-29.3))
    print(f"[Step 1] Boundary read: {len(polygons)} polygons, ~{km2:,.0f} km2 "
          f"(official: ~5,749 km2)")
    return area


# ---------------------------------------------------------------- Step 2
def sea_surface_temperature_name(ds):
    """Find the SST variable, whatever the file calls it."""
    for name in SST_VARIANTS:
        if name in ds.variables:
            return name
    raise SystemExit(f"[Step 2] No SST variable found. Variables: {list(ds.variables)}")


def to_celsius(values, sst_var):
    """Convert to degrees Celsius and blank impossible values."""
    arr = np.ma.filled(np.ma.asarray(values, dtype=float), np.nan)
    units = str(getattr(sst_var, "units", "K")).strip()
    if units.upper() in ("K", "KELVIN", "DEGREES_KELVIN", "DEG_K"):
        arr = arr - 273.15
    arr[(arr < -2.5) | (arr > 40)] = np.nan
    return arr


def axis_window(values, low, high, axis_name):
    """Coordinate values inside low..high, plus the slice that selects them in the file.

    Returns (ascending values, (start, stop), flipped). If the file stores the axis
    descending, the slice still indexes the file as stored and the caller has to
    reverse the slab it reads.
    """
    flipped = values[0] > values[-1]
    ascending = values[::-1] if flipped else values
    idx = np.where((ascending >= low) & (ascending <= high))[0]
    if idx.size == 0:
        raise SystemExit(f"[Step 2] The file does not cover the analysis box on {axis_name}: "
                         f"file range {values.min():.3f}..{values.max():.3f}, "
                         f"needed {low}..{high}")
    if flipped:
        start, stop = values.size - 1 - idx[-1], values.size - idx[0]
    else:
        start, stop = idx[0], idx[-1] + 1
    return ascending[idx], (start, stop), flipped


def read_grid(nc_path):
    """Coordinates, dates and a land mask sampled from the data.

    OSTIA is an L4 product with no gaps over the ocean, so a pixel that is NaN is
    land or sea ice. We sample a few dates instead of one, so a single missing day
    cannot be mistaken for land.
    """
    ds = nc.Dataset(nc_path)
    name = sea_surface_temperature_name(ds)
    var = ds.variables[name]
    lat_raw = np.asarray(ds.variables["latitude"][:], dtype=float)
    lon_raw = np.asarray(ds.variables["longitude"][:], dtype=float)
    lats, lat_slice, flip_lat = axis_window(lat_raw, LAT_MIN, LAT_MAX, "latitude")
    lons, lon_slice, flip_lon = axis_window(lon_raw, LON_MIN, LON_MAX, "longitude")

    tvar = ds.variables["time"]
    dates = nc.num2date(tvar[:], tvar.units)
    n = len(dates)
    print(f"[Step 2] SST variable '{name}' ({getattr(var, 'units', '?')}), "
          f"{n} dates ({dates[0]:%Y-%m-%d} to {dates[-1]:%Y-%m-%d})")

    # Sample up to 5 dates spread over the record to build the land mask
    sample = np.unique(np.linspace(0, n - 1, min(5, n)).astype(int))
    land = None
    for i in sample:
        slab = to_celsius(var[i, lat_slice[0]:lat_slice[1], lon_slice[0]:lon_slice[1]], var)
        if flip_lat:
            slab = slab[::-1]
        if flip_lon:
            slab = slab[:, ::-1]
        step_land = ~np.isfinite(slab)
        land = step_land if land is None else (land & step_land)
    print(f"[Step 2] Land/no-data pixels in the box: {int(land.sum())} of {land.size}")
    ds.close()
    return lats, lons, dates, land, lat_slice, lon_slice, name, flip_lat, flip_lon


# ---------------------------------------------------------------- Step 3
def distance_to_coast(land, lats, lons):
    """Distance (km) from each water pixel to the nearest land pixel."""
    km_per_deg_lat = 111.32
    km_per_deg_lon = 111.32 * np.cos(np.radians(np.mean(lats)))
    dy = abs(lats[1] - lats[0]) * km_per_deg_lat
    dx = abs(lons[1] - lons[0]) * km_per_deg_lon
    return distance_transform_edt(~land, sampling=(dy, dx))


# ---------------------------------------------------------------- Step 4
def pixels_inside(area, lats, lons):
    """True for each pixel whose centre is inside the protected area."""
    LON, LAT = np.meshgrid(lons, lats)
    inside = shapely.contains_xy(area, LON, LAT)
    print(f"[Step 4] Pixels inside the protected area: {inside.sum()}")
    return LON, LAT, inside


def zone_masks(dist_km, inside, sea):
    """Boolean masks for the distance bands and for the coastal/offshore pair."""
    masks = {}
    for name, near, far, inside_only in BANDS:
        band = sea & (dist_km > near) & (dist_km <= far)
        masks[name] = (band & inside) if inside_only else band
    masks["Coastal 0-10 km"] = masks["Coastal strip (0-10 km)"]
    masks["Offshore 50-150 km"] = masks["Offshore reference (50-150 km)"]
    return masks


# ---------------------------------------------------------------- Step 5
def stream_stats(nc_path, dist_km, inside, land, slices,
                 sst_name, flipped_lat, flipped_lon):
    """One pass over every date.

    Returns the time-mean SST field, one (date, coastal, offshore) row per month,
    and the per-band median of the time mean.
    """
    ds = nc.Dataset(nc_path)
    var = ds.variables[sst_name]
    lat_slice = slice(*slices[0])
    lon_slice = slice(*slices[1])
    tvar = ds.variables["time"]
    dates = nc.num2date(tvar[:], tvar.units)
    sea = ~land

    masks = zone_masks(dist_km, inside, sea)
    coastal = masks["Coastal 0-10 km"]
    offshore = masks["Offshore 50-150 km"]
    print(f"[Step 5] Cells: coastal {int(coastal.sum())}, offshore {int(offshore.sum())}, "
          f"inside the boundary {int((inside & sea).sum())}")

    total = np.zeros(land.shape, dtype=float)
    count = np.zeros(land.shape, dtype=float)
    monthly = []
    key, m_sum, m_cnt = None, None, None

    def flush(k, s, c):
        if s is None:
            return
        with np.errstate(invalid="ignore", divide="ignore"):
            field = np.where(c > 0, s / c, np.nan)
        cv = field[coastal]
        ov = field[offshore]
        if np.isfinite(cv).any() and np.isfinite(ov).any():
            monthly.append((k, float(np.nanmedian(cv)), float(np.nanmedian(ov))))

    for i, d in enumerate(dates):
        slab = to_celsius(var[i, lat_slice, lon_slice], var)
        if flipped_lat:
            slab = slab[::-1]
        if flipped_lon:
            slab = slab[:, ::-1]
        valid = np.isfinite(slab)
        total[valid] += slab[valid]
        count[valid] += 1

        k = (d.year, d.month)
        if m_sum is None or k != key:
            flush(key, m_sum, m_cnt)
            key, m_sum, m_cnt = k, np.zeros(land.shape), np.zeros(land.shape)
        m_sum[valid] += slab[valid]
        m_cnt[valid] += 1
        if (i + 1) % 365 == 0:
            print(f"         {i + 1}/{len(dates)} dates processed")
    flush(key, m_sum, m_cnt)

    with np.errstate(invalid="ignore", divide="ignore"):
        mean_sst = np.where(count > 0, total / count, np.nan)

    band_stats = []
    whole = mean_sst[inside & sea]
    band_stats.append(("Whole protected area", float(np.nanmedian(whole)), int(np.isfinite(whole).sum())))
    for name, *_ in BANDS:
        vals = mean_sst[masks[name]]
        band_stats.append((name, float(np.nanmedian(vals)), int(np.isfinite(vals).sum())))

    print("[Step 5] Median sea surface temperature (degC), time mean:")
    for name, value, n in band_stats:
        print(f"         {name:34s} {value:5.2f}   ({n} pixels)")

    ds.close()
    return mean_sst, monthly, band_stats


def compare_zones(band_stats):
    """Coastal minus offshore difference from the band table."""
    labelled = {name: value for name, value, _ in band_stats}
    coastal = labelled["Coastal strip (0-10 km)"]
    offshore = labelled["Offshore reference (50-150 km)"]
    print(f"         Coastal water is {offshore - coastal:.2f} degC colder than offshore")
    return offshore - coastal


def sector_stats(mean_sst, dist_km, inside, lats, land):
    """Median SST and coastal-offshore gradient for each latitude sector."""
    sea = ~land
    lon_shape = mean_sst.shape[1]
    LON, LAT = np.meshgrid(np.arange(lon_shape), lats)
    rows = []
    for name, lat_lo, lat_hi in SECTORS:
        band = inside & sea & (LAT >= lat_lo) & (LAT < lat_hi)
        coastal = band & (dist_km <= COAST_KM)
        offshore = sea & (dist_km > OFFSHORE_KM[0]) & (dist_km <= OFFSHORE_KM[1]) & \
                   (LAT >= lat_lo) & (LAT < lat_hi)
        if not np.isfinite(mean_sst[band]).any() or not np.isfinite(mean_sst[coastal]).any():
            continue
        sst = float(np.nanmedian(mean_sst[band]))
        grad = float(np.nanmedian(mean_sst[offshore]) - np.nanmedian(mean_sst[coastal]))
        rows.append((name, sst, grad, int(np.isfinite(mean_sst[band]).sum())))

    print("[Step 5] By sector (median SST, coastal minus offshore gradient):")
    for name, sst, grad, n in rows:
        print(f"         {name:50s} {sst:5.2f} degC   gradient {-grad:5.2f} degC  ({n} px)")
    return rows


# ---------------------------------------------------------------- Step 6
def draw_polygon(ax, area):
    parts = area.geoms if area.geom_type == "MultiPolygon" else [area]
    for p in parts:
        x, y = p.exterior.xy
        ax.plot(x, y, color="white", lw=2.2)
        ax.plot(x, y, color="crimson", lw=1.2)


def annotate_islands(ax):
    """Approximate island positions, added by hand for orientation (not from the data)."""
    islands = {"Isla Chanaral": (-71.58, -29.03),
               "Isla Damas": (-71.52, -29.23),
               "Isla Choros": (-71.54, -29.26)}
    for name, (x, y) in islands.items():
        ax.plot(x, y, "o", ms=4, mfc="white", mec="black")
        ax.annotate(name, (x, y), xytext=(6, 0), textcoords="offset points", fontsize=8,
                    bbox=dict(boxstyle="round,pad=0.2", fc="white", alpha=0.8, lw=0))


def draw_map(area, lats, lons, sst, period, out_png):
    fig, ax = plt.subplots(figsize=(8, 9), dpi=150)
    ax.set_facecolor("#d9d4c7")  # grey = land or no data
    mesh = ax.pcolormesh(lons, lats, sst, cmap="RdYlBu_r", shading="nearest")
    bar = fig.colorbar(mesh, ax=ax, shrink=0.8)
    bar.set_label("Sea surface temperature (Â°C)")

    draw_polygon(ax, area)
    annotate_islands(ax)

    ax.set_xlim(VIEW_LON_MIN, LON_MAX)
    ax.set_ylim(LAT_MIN, LAT_MAX)
    ax.set_aspect(1 / np.cos(np.radians(29.3)))
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    ax.set_title(f"Sea surface temperature, {period}\n"
                 "Humboldt Archipelago protected area (red outline)", fontsize=11)
    ax.text(0.01, 0.01,
            "Data: Copernicus Marine, SST_GLO_SST_L4_REP_OBSERVATIONS_010_011 (OSTIA, 5 km)\n"
            "Boundary: MMA Chile / SIMBIO. Grey = land or no data. Island positions approximate.",
            transform=ax.transAxes, fontsize=6, va="bottom",
            bbox=dict(fc="white", alpha=0.8, lw=0))
    plt.tight_layout()
    plt.savefig(out_png)
    print(f"[Step 6] Map saved: {out_png}")


def draw_cold_anomaly(area, lats, lons, sst, dist_km, period, out_png):
    """SST minus the offshore reference at the same latitude.

    Removes the north-south temperature gradient, so what is left is the
    cross-shore contrast: the cold band along the coast is the upwelling.
    """
    reference_band = (dist_km > OFFSHORE_KM[0]) & (dist_km <= OFFSHORE_KM[1])
    anomaly = np.full_like(sst, np.nan)
    for i in range(sst.shape[0]):
        row = sst[i, :]
        ref_cells = row[reference_band[i, :]]
        if np.isfinite(ref_cells).sum() > 5:
            reference = np.nanmedian(ref_cells)
            if np.isfinite(reference):
                anomaly[i, :] = row - reference

    fig, ax = plt.subplots(figsize=(8, 9), dpi=150)
    ax.set_facecolor("#d9d4c7")
    mesh = ax.pcolormesh(lons, lats, anomaly, cmap="RdBu_r",
                         norm=Normalize(vmin=-3, vmax=3), shading="nearest")
    bar = fig.colorbar(mesh, ax=ax, shrink=0.8)
    bar.set_label("Temperature minus offshore reference at same latitude (Â°C)")

    draw_polygon(ax, area)
    annotate_islands(ax)

    ax.set_xlim(VIEW_LON_MIN, LON_MAX)
    ax.set_ylim(LAT_MIN, LAT_MAX)
    ax.set_aspect(1 / np.cos(np.radians(29.3)))
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    ax.set_title(f"Upwelling signature, {period}\n"
                 "Blue = colder than open water at the same latitude", fontsize=11)
    ax.text(0.01, 0.01,
            "Data: Copernicus Marine, SST_GLO_SST_L4_REP_OBSERVATIONS_010_011 (OSTIA, 5 km)\n"
            "Reference: median of the 50-150 km band in the same row. Source data in kelvin (CF).",
            transform=ax.transAxes, fontsize=6, va="bottom",
            bbox=dict(fc="white", alpha=0.8, lw=0))
    plt.tight_layout()
    plt.savefig(out_png)
    print(f"[Step 6] Upwelling map saved: {out_png}")


def draw_seasonal(monthly, out_png, out_csv):
    """Average year: one value per calendar month, shading = lowest to highest year."""
    if len(monthly) < 24:
        print("[Step 6] Seasonal chart skipped: needs at least 24 months")
        return
    names = ["Jan", "Feb", "Mar", "Apr", "May", "Jun",
             "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]
    fig, ax = plt.subplots(figsize=(8, 4.5), dpi=150)
    colours = {"Coastal strip (0-10 km)": "#1b9e77", "Offshore (50-150 km)": "#3b4cc0"}
    series = {"Coastal strip (0-10 km)": {}, "Offshore (50-150 km)": {}}
    for (year, month), coastal, offshore in monthly:
        series["Coastal strip (0-10 km)"].setdefault(month, []).append(coastal)
        series["Offshore (50-150 km)"].setdefault(month, []).append(offshore)

    rows = []
    for zone, by_month in series.items():
        missing = [m for m in range(1, 13) if not by_month.get(m)]
        if missing:
            print(f"[Step 6] Seasonal chart: {zone} has no data for months {missing}, skipped")
            continue
        mean = np.array([np.mean(by_month[m]) for m in range(1, 13)])
        lo = np.array([np.min(by_month[m]) for m in range(1, 13)])
        hi = np.array([np.max(by_month[m]) for m in range(1, 13)])
        ax.plot(names, mean, "o-", color=colours[zone], label=zone)
        ax.fill_between(names, lo, hi, color=colours[zone], alpha=0.15)
        for m in range(12):
            rows.append([zone, names[m], round(mean[m], 3), round(lo[m], 3),
                         round(hi[m], 3), len(by_month[m + 1])])
        swing = mean.max() - mean.min()
        print(f"         {zone}: warmest {names[int(np.argmax(mean))]} ({mean.max():.2f} degC), "
              f"coldest {names[int(np.argmin(mean))]} ({mean.min():.2f} degC), swing {swing:.2f} degC")

    ax.set_ylabel("Sea surface temperature (Â°C)")
    years = sorted({key[0] for key, _, _ in monthly})
    ax.set_title(f"Seasonal cycle of sea surface temperature, {years[0]}-{years[-1]}\n"
                 "Line = average year; shading = lowest to highest year", fontsize=11)
    ax.legend(frameon=False)
    ax.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(out_png)
    with open(out_csv, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["zone", "month", "mean_degc", "coldest_year", "warmest_year", "years"])
        w.writerows(rows)
    print(f"[Step 6] Seasonal chart saved: {out_png}")


def draw_gradient(monthly, out_png, out_csv):
    """Monthly coastal-minus-offshore difference: the upwelling index."""
    import datetime as dt
    dates = [dt.date(key[0], key[1], 1) for key, _, _ in monthly]
    coastal = np.array([c for _, c, _ in monthly])
    offshore = np.array([o for _, _, o in monthly])
    gradient = coastal - offshore

    fig, ax = plt.subplots(figsize=(9, 4), dpi=150)
    ax.axhline(0, color="black", lw=0.8)
    ax.plot(dates, gradient, "o-", ms=2.5, lw=0.8, color="#1b9e77", alpha=0.45,
            label="Monthly")
    if gradient.size >= 12:
        smooth = np.convolve(gradient, np.ones(12) / 12, mode="valid")
        ax.plot(dates[11:], smooth, lw=2.2, color="#0b4f3c", label="12-month mean")
    ax.set_ylabel("Coastal minus offshore SST (Â°C)")
    ax.set_title("Upwelling index: colder coastal water when the wind lifts deep water\n"
                 "Negative values = upwelling is active", fontsize=11)
    ax.legend(frameon=False)
    ax.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(out_png)
    with open(out_csv, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["date", "coastal_degc", "offshore_degc", "coastal_minus_offshore_degc"])
        for d, c, o, g in zip(dates, coastal, offshore, gradient):
            w.writerow([d.isoformat(), round(c, 3), round(o, 3), round(g, 3)])
    print(f"[Step 6] Gradient chart saved: {out_png}")


def save_table(rows, out_csv, value_column="median_sst_degc"):
    with open(out_csv, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["zone", value_column, "pixels"])
        for name, value, n in rows:
            w.writerow([name, round(value, 3), n])
    print(f"[Step 6] Table saved: {out_csv}")


def save_sectors(rows, out_csv):
    with open(out_csv, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["sector", "median_sst_degc", "coastal_minus_offshore_degc", "pixels"])
        for name, sst, grad, n in rows:
            w.writerow([name, round(sst, 3), round(-grad, 3), n])
    print(f"[Step 6] Sector table saved: {out_csv}")


# ---------------------------------------------------------------- main
if __name__ == "__main__":
    if len(sys.argv) != 3:
        print(__doc__)
        sys.exit(1)
    kmz_path, nc_path = sys.argv[1], sys.argv[2]

    area = read_boundary(kmz_path)
    lats, lons, dates, land, lat_slice, lon_slice, sst_name, flip_lat, flip_lon = read_grid(nc_path)
    dist_km = distance_to_coast(land, lats, lons)
    LON, LAT, inside = pixels_inside(area, lats, lons)

    mean_sst, monthly, band_stats = stream_stats(
        nc_path, dist_km, inside, land,
        (lat_slice, lon_slice), sst_name, flip_lat, flip_lon)
    compare_zones(band_stats)
    sectors = sector_stats(mean_sst, dist_km, inside, lats, land)

    period = (f"{dates[0]:%B %Y}" if len(dates) == 1
              else f"{dates[0]:%b %Y} - {dates[-1]:%b %Y} average")
    draw_map(area, lats, lons, mean_sst, period, "humboldt_sst_map.png")
    draw_cold_anomaly(area, lats, lons, mean_sst, dist_km, period, "humboldt_sst_upwelling.png")
    draw_seasonal(monthly, "humboldt_sst_seasonal.png", "humboldt_sst_seasonal.csv")
    draw_gradient(monthly, "humboldt_sst_gradient.png", "humboldt_sst_monthly.csv")
    save_table(band_stats, "humboldt_sst_by_zone.csv")
    save_sectors(sectors, "humboldt_sst_by_sector.csv")


