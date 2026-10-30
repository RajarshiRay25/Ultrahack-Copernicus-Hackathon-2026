# 🐧 Humboldt Archipelago: Following the Ocean's Engine from Wind to Wildlife
## Members and Contributors

**Rajarshi Ray - Tampere University**  || 
**Yunseon Lee -  University of Helsinki** || 
**Rodrigo Andrés Jara Godoy -  Change of North Chile**

### Demo - https://claude.ai/artifact/TvYFjbv2wNTrNroCiFCrQC#plan 
**Team Ocean Star ⭐ | Copernicus Marine Service Hackathon 2026**

> Using free satellite data from Copernicus Marine, we show **where** and **when** Chile's Humboldt Archipelago marine protected area is most alive, and **why**. The result is simple, transparent evidence that managers can use to decide where and when to focus protection.

---

## 📌 In short

- **Wind** blows north along the coast all year, strongest in **September–November**.
- That wind brings **cold, nutrient-rich water** up from the deep near the coast (upwelling): the coast is **about 1 °C colder than the open ocean in 137 of 138 months**.
- Tiny sea plants (phytoplankton) feed on those nutrients: **chlorophyll is about 4 times higher near the coast** than offshore, and peaks in **October**.
- **Penguins, sea lions, dolphins and whales** are seen mostly in the same coastal strip.
- **Conclusion:** the coastal strip around the islands is the heart of the protected area, and spring is its most important season.

---

## 🌊 The problem

In 2023, Chile created the **Humboldt Archipelago Multiple-Use Marine Protected Area**: about **5,750 km²** of ocean in the Atacama and Coquimbo regions. "Multiple-use" means people can still fish and work there; it is not a no-go zone.

The area protects:
1. native species of the Humboldt Current and the places where they feed and travel,
2. islands and islets,
3. the **upwelling** that feeds the whole food web,
4. sustainable small-scale (artisanal) fishing.

The area **has no management plan yet**. Nobody can protect 5,750 km² equally, so managers must decide where and when to focus. Many groups have interests here (fishers, tourism, industry, conservation), so that decision needs **evidence anyone can check**.

**Our question:** which parts of the protected area, and which times of year, matter most for sea life, and why?

---

## 🔗 How this ocean works (the chain we follow)

```
🌬️ WIND             blows north along the coast, pushes surface water out to sea
      ↓
🌡️ COLD WATER       rises from the deep to replace it, bringing nutrients (upwelling)
      ↓
🌿 PHYTOPLANKTON    tiny sea plants feed on the nutrients and grow (seen as chlorophyll)
      ↓
🐟 SMALL FISH       anchovies and sardines eat the plankton
      ↓
🐧🐋 BIG ANIMALS     penguins, sea lions, dolphins and whales eat the fish
```

Each analysis in this project measures one link in this chain.

| Analysis | Link in the chain | Question it answers |
|---|---|---|
| 🌬️ Wind | Cause | What drives the upwelling, and when is it strongest? |
| 🌡️ Sea temperature | Evidence | Where does the cold deep water reach the surface? |
| 🌿 Chlorophyll | Effect | Where and when is there the most food at the base of the food web? |
| 🐧 Marine life | Who benefits | Are the animals found where the food is? |

---

## 📊 Data

All data is free and public.

| Data | Source | Details |
|---|---|---|
| Chlorophyll-a | Copernicus Marine `OCEANCOLOUR_GLO_BGC_L4_MY_009_104`, dataset `cmems_obs-oc_glo_bgc-plankton_my_l4-multi-4km_P1M` | Satellite, 4 km squares, monthly, Sep 2015 – Sep 2026 |
| Sea surface temperature | Copernicus Marine `SST_GLO_SST_L4_REP_OBSERVATIONS_010_011` (OSTIA) | Satellite, about 5 km squares, averaged to monthly, Jan 2015 – Jun 2026 |
| Wind | Copernicus Marine `WIND_GLO_PHY_CLIMATE_L4_MY_012_003` | Satellite + model, 25 km squares, monthly, Jan 2015 – Jun 2026 |
| Animal sightings | [GBIF](https://www.gbif.org) occurrence download | eBird, iNaturalist, Happywhale, Centro Ballena Azul and others |
| Protected-area boundary | Ministry of Environment of Chile, [SIMBIO](https://simbio.mma.gob.cl/CbaAP/Details/2082) | Official KMZ file |

Study box for downloads: latitude **−30.2 to −28.5**, longitude **−72.5 to −71.0** (wind: −31.0 to −27.5 and −74.0 to −71.0, because wind squares are larger).

---

## 🔍 Results

### 1. 🌬️ Wind: the engine runs all year and peaks in spring

The **northward** part of the wind is the "upwelling wind": on this north-south coast, it pushes surface water offshore. We sorted it into simple categories:

| Category | Northward wind | Meaning |
|---|---|---|
| Downwelling | below 0 m/s | no upwelling |
| Weak | 0–5 m/s | little upwelling |
| Moderate | 5–7 m/s | steady upwelling |
| Strong | 7 m/s or more | intense upwelling |

*These thresholds are our own transparent choice, written at the top of the code so anyone can see or change them.*

| Season | Northward wind | Category |
|---|---|---|
| May–June | about 5.5 m/s | Moderate (calmest) |
| December–April | 6.4–6.9 m/s | Moderate |
| **September–November** | **7.0–7.7 m/s (peak in October)** | **Strong** |

- The wind blows north **in every month**, so upwelling never stops on average.
- It is **strongest in October**, the same month as the chlorophyll peak, and **calmest in May–June**, when chlorophyll is lowest (seasonal match: correlation r = 0.91).
- Some years are much weaker: in the worst autumn–winter months, the wind fell to 1–2 m/s.

![Wind season](humboldt_wind_seasonal.png)
![Wind category by latitude and month](humboldt_wind_latitude_month.png)

### 2. 🌡️ Sea temperature: cold water rises near the coast

**What we look for:** if upwelling happens, the water near the coast should be **colder** than the water offshore, which is the opposite of normal, since shallow water usually warms faster. The temperature difference between coastal and offshore water is our **upwelling indicator**: the more negative the difference, the colder the coast compared with the open ocean.

**By distance from the coast** (median, January 2015 – June 2026):

| Zone | Sea temperature (°C) |
|---|---|
| **Coastal strip (0–10 km)** | **14.80** |
| Inner zone (10–25 km) | 14.82 |
| Middle zone (25–50 km) | 15.10 |
| Offshore reference (50–150 km) | 15.82 |
| **Difference (coast − offshore)** | **−1.0** |

**By sector, north to south** (inside the protected area):

| Sector | Sea temperature (°C) | Coast − offshore (°C) |
|---|---|---|
| 1: Isla Chañaral (28.89–29.10° S) | 15.18 | **−1.62** (strongest) |
| 2: Punta de Choros (29.10–29.35° S) | 15.08 | −1.23 |
| 3: Chañaral de Aceituno (29.35–29.55° S) | 14.94 | −0.97 |
| 4: Southern sector (29.55–29.76° S) | 14.91 | −0.85 |

**By season** (average year):

| Months | Coastal strip (°C) | Offshore (°C) | Coast − offshore (°C) |
|---|---|---|---|
| Jan–Feb (summer) | 17.2–17.6 | 18.0–18.7 | −0.8 to −1.2 |
| **Mar–May (autumn)** | 14.0–16.5 | 15.9–18.2 | **−1.7 to −2.1** (largest) |
| Jun–Jul (winter) | 13.2–13.5 | 14.0–14.9 | −0.8 to −1.4 |
| **Aug–Oct (late winter–spring)** | **13.1–13.5 (coldest)** | 13.6–14.0 | −0.5 |
| Nov–Dec | 14.6–15.7 | 15.1–16.4 | −0.5 to −0.7 |

**Key findings**

- **The coast is colder than the open ocean almost all the time:** in **137 of 138 months**, by about **1 °C** on average. This is direct evidence that upwelling brings cold deep water to the surface all year.
- **The coastal water is coldest in August–October** (about 13.1 °C), the same period when the upwelling wind becomes strong and chlorophyll rises to its October peak.
- **The coast-offshore contrast is largest in autumn (March–May), not in spring.** In spring the whole region is cold, so the difference is small even though upwelling is strong; in summer and autumn the sun warms the open ocean while the upwelled coastal water stays cold, so the difference stands out more. The contrast therefore shows *that* upwelling is happening, but it is not a direct measure of *how strong* it is.
- **The northern sector (Isla Chañaral) shows the strongest cold signal**, getting weaker toward the south. Chlorophyll, however, is highest from Isla Damas and Isla Choros southward. Cold water and plankton do not line up exactly: phytoplankton need time to grow after the nutrients arrive, and calmer bays can hold the water long enough for blooms to build up.
- **Year to year:** 2023, an El Niño year, had warmer coastal water (about 15.0 °C on average) than 2019–2022 (14.2–14.7 °C).

Data tables: `humboldt_sst_by_zone.csv`, `humboldt_sst_by_sector.csv`, `humboldt_sst_seasonal.csv`, `humboldt_sst_monthly.csv`

<!-- Add the temperature figures here, for example:
![Temperature map](humboldt_sst_map.png)
![Temperature season](humboldt_sst_seasonal.png) -->

### 3. 🌿 Chlorophyll: the coastal strip is the productive core

11-year average (September 2015 – September 2026, 133 months), median inside the protected area:

| Zone | Chlorophyll-a (mg/m³) |
|---|---|
| **Coastal strip (0–10 km)** | **2.95** |
| Inner zone (10–25 km) | 1.38 |
| Middle zone (25–50 km) | 0.84 |
| Offshore (over 50 km) | 0.68 |

- The coastal strip has **about 4.3 times more** chlorophyll than offshore.
- The richest water runs along the coast from **Isla Damas and Isla Choros southward**.
- **Season:** lowest in winter (about 1.4 mg/m³ in May–July), highest in spring (**about 6.0 mg/m³ in October**). Offshore stays low all year.

![Chlorophyll map](humboldt_chlorophyll_map.png)
![Chlorophyll season](humboldt_chlorophyll_seasonal.png)

### 4. 🐧 Marine life: animals are found where the food is

We used **1,740 sightings** (2015–2026) of Humboldt penguins, South American sea lions, bottlenose dolphins, and blue, fin and humpback whales.

| Zone (inside the protected area) | Sightings | Sightings per 100 km² |
|---|---|---|
| **Coast & islands (0–10 km)** | **518** | **64.2** |
| Inner zone (10–25 km) | 107 | 6.6 |
| Middle zone (25–50 km) | 10 | 0.4 |
| Offshore (over 50 km) | 1 | 0.1 |

- Sightings follow the **same coast-to-offshore pattern** as chlorophyll.
- Whales gather around **Isla Chañaral** and **Isla Choros and Damas**.
- Sightings peak in **January–February**, two to four months after the October chlorophyll peak. This fits the time energy takes to move up the food chain, but it is also the summer holiday season, when the most people are watching.

![Wildlife map](humboldt_wildlife_map.png)
![Wildlife season](humboldt_wildlife_seasonal.png)

---

## 🧩 How the analyses connect

All datasets share **location** (latitude and longitude) and **time** (dates), so they can be placed on the same map and the same calendar.

- **Same zones:** every satellite square gets its distance from the coast. Animal sightings are matched to the square they fall in, so chlorophyll, temperature and sightings are all compared in the same zones (0–10 km, 10–25 km, 25–50 km, over 50 km).
- **Same calendar:** every dataset is averaged by month (January to December), so we can see which comes first.

**The full picture:**

| | Peak season | Where it is strongest |
|---|---|---|
| 🌬️ Wind | October | Along the whole coast |
| 🌡️ Cold water | Coldest coast Aug–Oct; largest coast-offshore contrast Mar–May | Coastal strip, strongest near Isla Chañaral |
| 🌿 Chlorophyll | October | Coastal strip, near the islands |
| 🐧 Animal sightings | January–February | Coastal strip, near the islands |

---

## 🎯 What this means for managing the protected area

This project supports the hackathon goal of a **resilient, sustainable and well-governed ocean**:

1. **Focus on the coastal strip around the islands.** It holds most of the food and most of the wildlife, while making up only a small part of the area.
2. **Follow the calendar, not just the map.** September–November is when the ecosystem's engine runs hardest. Extra care in this season (for example around fishing pressure and disturbance near breeding colonies) protects the most life while leaving the rest of the year open for normal use.
3. **Use the wind as an early warning.** Copernicus wind data is free and updated regularly. A weak spring wind means less food for wildlife in the following months, giving managers time to act, for example by easing fishing pressure on anchovy in a poor year.
4. **Keep decisions transparent.** All steps, thresholds and data sources are written in the code and this README, so fishers, government and conservation groups can see exactly how each conclusion was reached.
5. **The method travels.** The same code can be used for any marine protected area, including in Europe, by changing the boundary file and the study box.

---

## ⚠️ Limitations

- **Satellites only see the surface.** Deep-sea habitats, which make up most of the area, are not measured.
- **Animal sightings show where people looked,** not only where animals are. Most observers are on the coast, on boat tours near the islands, and out in summer. Some sightings even appear on land because people record their own position.
- **The wind data is from about 180 km offshore** (one line of squares), not inside the area itself. Coastal winds may differ in detail.
- **The wind–chlorophyll match is clear for the average season but weak month by month.** Other factors, such as El Niño, also affect individual years.
- **Distance to the coast** is measured on the 4 km satellite grid, so it is approximate.
- **The wind categories** are our own choice, not an official standard.
- **The temperature contrast (coast − offshore)** shows that upwelling is happening but is affected by how much the sun warms the open ocean, so it is not a direct measure of upwelling strength.

---

## 🚀 Next steps

- Download wind for the full study box to map wind over the protected area itself.
- Add fishing activity (Global Fishing Watch) as a pressure layer.
- Combine wind, temperature, chlorophyll and wildlife into **one transparent priority map**.
- Use systematic survey data instead of opportunistic sightings.
- Test the method on a European marine protected area.

---

## 💻 How to run the code

**Requirements:** Python 3 and these libraries:

```bash
pip install netCDF4 shapely matplotlib numpy scipy
```

**Run each analysis** (keep all scripts in the same folder; the wildlife script reuses the chlorophyll script):

```bash
python humboldt_chlorophyll_analysis.py  boundary.kmz  chlorophyll.nc
python humboldt_wildlife_analysis.py     boundary.kmz  chlorophyll.nc  gbif_sightings.csv
python humboldt_wind_analysis.py         boundary.kmz  wind.nc
```

Each script prints its steps as it runs and saves its maps (PNG) and tables (CSV) in the current folder. Run the chlorophyll script before the wildlife script, because the wildlife script uses its seasonal table.

### Repository contents

| File | What it is |
|---|---|
| `humboldt_chlorophyll_analysis.py` | Chlorophyll: boundary, zones, map, seasonal cycle |
| `humboldt_wildlife_analysis.py` | Animal sightings compared with chlorophyll |
| `humboldt_wind_analysis.py` | Wind categories by latitude and month, seasonal cycle |
| _temperature script_ | _to be added: sea temperature by zone, sector and season_ |
| `*.png` | Maps and charts |
| `*.csv` | Result tables |

### How the code is organised

Every script follows the same simple pattern:

1. **Description** at the top: the question, the data and the steps, in plain English.
2. **Settings:** study box, zones and categories, all in one place.
3. **One function per step:** read the boundary → read the data → cut to the area → compare zones → draw maps and charts → save tables.
4. **Start button** at the bottom (`if __name__ == "__main__":`), which runs the steps in order.

---


## 🙏 Acknowledgements

- Copernicus Marine Service for free ocean data
- GBIF and all its contributors (eBird, iNaturalist, Happywhale, Centro Ballena Azul and others)
- Ministry of Environment of Chile (SIMBIO) for the official boundary
