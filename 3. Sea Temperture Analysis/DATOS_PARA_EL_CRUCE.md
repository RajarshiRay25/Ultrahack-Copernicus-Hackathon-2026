# Humboldt Archipelago — paquete de datos para el cruce de variables

Copernicus Marine Service Hackathon 2026 · equipo Ocean Protector

Este documento define el contrato de datos para cruzar **viento → temperatura → clorofila**
sin que cada persona invente su propia versión de las zonas, las unidades o la llave temporal.

---

## 1. Qué se comparte

### Núcleo (imprescindible)

| Archivo | Tamaño | Contenido |
|---|---|---|
| `humboldt_amp.kmz` | 700 KB | Polígono oficial del AMP (SIMBIO 2082). Fuente única de la geometría |
| `humboldt_sst_analysis.py` | 24 KB | Script que produce las series de SST |
| `humboldt_sst_monthly.csv` | 4,5 KB | **Serie mensual de SST** — la tabla del cruce |
| `humboldt_sst_by_zone.csv` | 0,2 KB | Medianas por banda de distancia a la costa |
| `humboldt_sst_by_sector.csv` | 0,3 KB | Medianas y gradiente por sector del AMP |
| `humboldt_sst_seasonal.csv` | 1,3 KB | Ciclo estacional por mes calendario |
| `humboldt_chlorophyll_seasonal.csv` | 1,2 KB | Ciclo estacional de clorofila (del análisis previo) |
| `humboldt_chlorophyll_by_zone.csv` | 0,2 KB | Medianas de clorofila por banda |
| `humboldt_sst_*.png` | ~460 KB | Figuras listas para el deck |

### Solo si alguien necesita recalcular desde el dato crudo

| Archivo | Tamaño | Contenido |
|---|---|---|
| `sst_humboldt.nc` | 141 MB | SST diaria, 2015-01-01 a 2026-06-30, 0.05° |
| `check_sst.py` | 3 KB | Valida que el `.nc` cubra la caja y las fechas |
| `download_sst.py`, `download_sst_direct.py` | 12 KB | Reproducen la descarga (con cuenta / sin cuenta) |

### No compartir

- Credenciales, `~/.copernicusmarine/`, variables de entorno.
- `.venv/` y `.vscode/` (rutas absolutas de una máquina).

---

## 2. Llave del cruce

**Todas las series mensuales deben usar `date` en formato `YYYY-MM-01`** (primer día del mes, ISO).

| Archivo | Columna llave | Ejemplo |
|---|---|---|
| `humboldt_sst_monthly.csv` | `date` | `2015-01-01` |
| clorofila mensual (por generar) | `date` | `2015-01-01` |
| viento mensual (por generar) | `date` | `2015-01-01` |

Los CSV estacionales (`seasonal.csv`) usan `month` con el nombre del mes (`Jan`…`Dec`) y **no
sirven para el cruce temporal**: son para gráficos de ciclo anual.

### Lo que falta generar

Ni la clorofila ni el viento tienen todavía serie mensual por fecha. Se necesitan con esta forma:

```
date,zone,value
```

o, mejor, una columna por zona en un solo archivo por variable:

```
date,coastal,offshore
```

Regla: **mismas zonas, mismo orden, mismo rango de fechas**.

---

## 3. Definición de zonas (no modificar)

Distancia de cada píxel al píxel de tierra más cercano, en km:

| Banda | Distancia | Dentro del AMP |
|---|---|---|
| Coastal strip | 0–10 km | Sí |
| Inner zone | 10–25 km | Sí |
| Middle zone | 25–50 km | Sí |
| Offshore reference | 50–150 km | **No** |

La referencia offshore se mide **fuera** del polígono a propósito: el AMP tiene ~55 km de ancho,
así que "offshore dentro del AMP" degenera en una franja de una celda en el borde occidental.

Sectores del AMP, por latitud:

| Sector | Latitudes | Referencia aproximada |
|---|---|---|
| 1 | 28,88–29,10°S | Isla Chañaral |
| 2 | 29,10–29,35°S | Punta de Choros / Isla Damas |
| 3 | 29,35–29,55°S | Chañaral de Aceituno |
| 4 | 29,55–29,76°S | sector sur |

---

## 4. Unidades y convenciones

- **SST:** grados Celsius en todas las tablas. El NetCDF viene en kelvin (`analysed_sst`, offset 273,15).
- **Clorofila:** mg/m³.
- **Viento:** m/s, y para el índice de afloramiento interesa la componente **norte** (v > 0 = viento hacia el ecuador = favorable).
- **Gradiente:** `coastal − offshore`. Negativo = agua costera más fría = afloramiento activo.
- Mes en UTC.

---

## 5. Cómo cruzar (pandas)

```python
import pandas as pd

sst = pd.read_csv("humboldt_sst_monthly.csv", parse_dates=["date"])
chl = pd.read_csv("humboldt_chlorophyll_monthly.csv", parse_dates=["date"])
wind = pd.read_csv("humboldt_wind_monthly.csv", parse_dates=["date"])

tabla = (sst.merge(chl, on="date", how="inner", suffixes=("_sst", "_chl"))
            .merge(wind, on="date", how="inner"))
tabla["month"] = tabla["date"].dt.month

print(tabla.corr(numeric_only=True))
```

Para el desfase, correlaciona la serie desfasada:

```python
for lag in range(0, 4):
    r = tabla["coastal_degc"].shift(lag).corr(tabla["coastal_chl"])
    print(f"lag {lag} meses: r = {r:.2f}")
```

---

## 6. Resultados ya medidos (para no repetir trabajo)

- SST media del AMP 2015–2026: **14,99 °C**; franja costera **14,80 °C**; referencia offshore **15,82 °C**.
- El agua costera es **1,02 °C más fría** que el océano abierto (media de 11 años).
- Clorofila costera ≈ **4,3×** la offshore (análisis previo, 0,9 mg/m³ de media en el AMP).
- SST costera: mínimo en **septiembre** (13,1 °C), máximo en **febrero** (17,6 °C).
- Clorofila costera: máximo en **octubre** (6,0 mg/m³) → desfase de ~1 mes respecto del mínimo térmico.
- Ranking de sectores por gradiente térmico: **norte más fuerte** (−1,62 °C en Isla Chañaral) que sur (−0,85 °C).

---

## 7. Advertencias para el análisis

1. **El 1 °C de contraste térmico es menor que el 4,3× de clorofila.** Es esperable: OSTIA es un L4 suavizado, la franja de 10 km son ~2 celdas de 5,5 km, y el afloramiento es episódico mientras esto es una media de 11 años.
2. **La franja costera tiene solo 25 celdas** dentro del AMP. Cualquier estadístico por sector descansa en pocas celdas: reportar n y probar la sensibilidad con banda 0–20 km.
3. **Resoluciones distintas:** SST 0,05° (~5,5 km), clorofila 4 km, viento 0,25° (~28 km). El viento **no puede** tener franja costera de 0–10 km: úsalo como índice regional.
4. **Media móvil de 12 meses rezagada**, no centrada: al anotar años de El Niño, desplazar la etiqueta.
5. El producto SST termina el **30-06-2026**; julio en adelante requiere el producto NRT y no está incluido.

---

## 8. Procedencia

| Dato | Producto | Variable | Periodo |
|---|---|---|---|
| SST | `SST_GLO_SST_L4_REP_OBSERVATIONS_010_011` (OSTIA, 0.05°, diario) | `analysed_sst` (K) | 2015-01-01 – 2026-06-30 |
| Clorofila | `OCEANCOLOUR_GLO_BGC_L4_MY_009_104` (GlobColour, 4 km, mensual) | `CHL` (mg/m³) | 2015-09 – 2026-09 |
| Polígono | SIMBIO, registro 2082, Decreto 31/2023 MMA | — | 2023 |

Procesado con `humboldt_sst_analysis.py`; área del polígono verificada en 5.828 km²
contra 5.749 km² oficiales (+1,4% por proyección).
