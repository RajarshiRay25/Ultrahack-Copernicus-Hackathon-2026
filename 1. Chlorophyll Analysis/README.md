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
| Protected-area boundary (KMZ) | Ministry of Environment of Chile, SIMBIO: https://simbio.mma.gob.cl/CbaAP/Details/2082 |

## Method

1. Read the official boundary polygon.
2. Cut chlorophyll to a box around the area and average it over all months in the file.
3. Keep only the satellite pixels inside the boundary.
4. Measure each pixel's distance to the coast and compare chlorophyll across distance bands.
5. Produce a map and a results table.
6. Build the seasonal cycle: average chlorophyll for each calendar month across all years, coastal strip vs offshore.

Every step is a separate, commented function in `humboldt_chlorophyll_analysis.py`, so anyone can check or rerun it.

```
pip install netCDF4 shapely matplotlib numpy scipy
python humboldt_chlorophyll_analysis.py boundary.kmz chlorophyll.nc
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

## Limitations

- Satellites see only the surface layer; deep-water habitats, which make up most of the area, are not represented by this measure.
- Chlorophyll measures food supply, not animals directly. Species observations (for example from OBIS or GBIF) would strengthen the link.
- Distance to coast is computed on the 4 km satellite grid, so it is approximate.

## Next steps

- Sea surface temperature (Copernicus OSTIA) to map upwelling directly.
- Overlay of species observations and fishing activity.
- Apply the same method to other marine protected areas, including in Europe.
