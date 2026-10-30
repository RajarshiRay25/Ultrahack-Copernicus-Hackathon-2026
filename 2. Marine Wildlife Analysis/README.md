# Humboldt Archipelago: where is the ocean most productive?

Copernicus Marine Service Hackathon 2026

## The problem

In 2023 Chile created the Humboldt Archipelago Multiple-Use Marine Protected Area (about 5,750 km², Atacama and Coquimbo regions). It protects native Humboldt Current species, island ecosystems, coastal upwelling and sustainable artisanal fishing. It has no management plan yet, so managers still have to decide where to focus limited resources. Because many groups use this area (fishers, tourism, industry, conservation), that decision needs evidence anyone can check.

## The ecological reasoning

Coastal winds push surface water offshore, and cold, nutrient-rich water rises from below (upwelling). Phytoplankton bloom on those nutrients, which satellites detect as chlorophyll-a. Phytoplankton feed zooplankton and small fish such as anchovy, which in turn feed Humboldt penguins, seabirds, sea lions, dolphins and whales. Chlorophyll-a is therefore a measurable proxy for where the base of the food web is concentrated.

## Data

| Data | Source |
|---|---|
| Chlorophyll-a (CHL, mg/m³), 4 km, monthly | Copernicus Marine, OCEANCOLOUR_GLO_BGC_L4_MY_009_104, dataset `cmems_obs-oc_glo_bgc-plankton_my_l4-multi-4km_P1M` |
| Animal sightings (CSV) | GBIF occurrence download (eBird, iNaturalist, Happywhale, Centro Ballena Azul and others) |
| Protected-area boundary (KMZ) | Ministry of Environment of Chile, SIMBIO: https://simbio.mma.gob.cl/CbaAP/Details/2082 |

## Method

1. Read the official boundary polygon.
2. Cut chlorophyll to a box around the area and average it over all months in the file.
3. Keep only the satellite pixels inside the boundary.
4. Measure each pixel's distance to the coast and compare chlorophyll across distance bands.
5. Produce a map and a results table.
6. Build the seasonal cycle: average chlorophyll for each calendar month across all years, coastal strip vs offshore.

Every step is a separate, commented function in `humboldt_chlorophyll_analysis.py` and `humboldt_wildlife_analysis.py`, so anyone can check or rerun it.

```
pip install netCDF4 shapely matplotlib numpy scipy
python humboldt_chlorophyll_analysis.py boundary.kmz chlorophyll.nc
python humboldt_wildlife_analysis.py boundary.kmz chlorophyll.nc gbif_sightings.csv
```

## Results (September 2015 – September 2026, 133 months)

### 1. Where: the coastal strip is the productive core

| Zone | Median chlorophyll-a (mg/m³), 11-year average |
|---|---|
| Whole protected area | 0.94 |
| Coastal strip (0–10 km) | 2.95 |
| Inner zone (10–25 km) | 1.38 |
| Middle zone (25–50 km) | 0.84 |
| Offshore (over 50 km) | 0.68 |

Averaged over 11 years, the coastal strip has about 4.3 times more chlorophyll than the offshore zone. The richest water lies along the coast from Isla Damas and Isla Choros southward. The pattern is consistent over the long term, not a single-month event.

Figure: `humboldt_chlorophyll_map.png`

### 2. When: productivity peaks in austral spring

In the coastal strip, chlorophyll is lowest in winter (about 1.4 mg/m³ in May–July) and peaks in spring (about 6.0 mg/m³ in October), roughly a fourfold seasonal swing. Offshore stays low all year (0.4–1.2 mg/m³). Year-to-year variation is large near the coast, so single months can mislead and multi-year data is needed.

Figure: `humboldt_chlorophyll_seasonal.png`

### 3. Who: animal sightings concentrate in the productive coastal strip

GBIF sightings (2015–2026) of Humboldt penguin, South American sea lion, bottlenose dolphin, and blue, fin and humpback whales: 1,740 records, 636 inside or on the edge of the protected area.

| Zone (inside the protected area) | Sightings | Area (km²) | Sightings per 100 km² |
|---|---|---|---|
| Coast & islands (0–10 km) | 518 | 806 | 64.2 |
| Inner zone (10–25 km) | 107 | 1,613 | 6.6 |
| Middle zone (25–50 km) | 10 | 2,476 | 0.4 |
| Offshore (over 50 km) | 1 | 1,013 | 0.1 |

Sightings follow the same coast-to-offshore pattern as chlorophyll, much more steeply. Whale sightings cluster around Isla Chañaral and Isla Choros and Damas. By month, sightings of all groups peak in January–February, two to four months after the October chlorophyll peak.

Figures: `humboldt_wildlife_map.png`, `humboldt_wildlife_seasonal.png`

**Caution:** sightings record where and when people looked. Most observers are on the coast, on boat tours near the islands, and in the summer holidays (December–February), so these numbers overstate the coastal and summer pattern. The time lag between the chlorophyll peak and the predator peak is biologically plausible (energy takes time to move up the food web) but cannot be separated from observer effort with this data.

## How the chlorophyll and animal data connect

The two datasets have different shapes but share two things: **location** and **time**.

| | Chlorophyll (Copernicus) | Animal sightings (GBIF) |
|---|---|---|
| Shape | A grid of 4 × 4 km squares covering the sea | A list of points, one row per sighting |
| Location | Latitude and longitude of each square | Latitude and longitude of each sighting |
| Time | One value per square for every month, 2015–2026 | Date of each sighting |
| What it tells us | How much food (phytoplankton) is in each square | Which animal was seen, where and when |

**Linking by location.** Each sighting is matched to the chlorophyll square it falls in. Each square already has a measured distance to the coast, so every sighting inherits that distance and is placed in the same zones used for chlorophyll (0–10 km, 10–25 km, 25–50 km, over 50 km). Sightings per zone are divided by the zone's sea area, because the offshore zones are much larger. The chlorophyll table and the sightings table then use identical zones and can be compared directly. The map shows the same link visually, with sightings drawn on top of the 11-year chlorophyll average.

**Linking by time.** Each sighting's date gives its calendar month. Sightings are counted per month across all years and placed under the monthly chlorophyll cycle on a shared January–December axis.

**What the comparison shows.**

- *Space:* both chlorophyll and sightings decrease from the coast to offshore. Chlorophyll is about 4 times higher in the coastal strip than offshore; sightings per 100 km² are hundreds of times higher. Animals are concentrated where the food supply is highest.
- *Time:* chlorophyll peaks in October, and sightings peak in January–February. This two-to-four-month delay fits the time energy needs to move up the food chain (phytoplankton → zooplankton and small fish → penguins, sea lions, dolphins, whales), but it also overlaps with the summer holiday season, when the most people are watching.

**What the link does and does not prove.** The data shows that animals and food are found in the same places, which supports the ecological chain described above. It does not prove the animals are there *because* of the food, and the strength of the coastal pattern in sightings is inflated by where observers go. A systematic survey (fixed boat transects covering all zones in all seasons) would be needed to measure animal distribution without observer bias.

## Limitations

- Satellites see only the surface layer; deep-water habitats, which make up most of the area, are not represented by this measure.
- Animal sightings are opportunistic (citizen science and tour boats), not systematic surveys, so they reflect observer effort as well as animal distribution.
- Distance to coast is computed on the 4 km satellite grid, so it is approximate.

## Next steps

- Sea surface temperature (Copernicus OSTIA) to map upwelling directly.
- Fishing activity (Global Fishing Watch) as a pressure layer.
- Combine wind, temperature, chlorophyll and wildlife into one transparent priority map.
- Apply the same method to other marine protected areas, including in Europe.
