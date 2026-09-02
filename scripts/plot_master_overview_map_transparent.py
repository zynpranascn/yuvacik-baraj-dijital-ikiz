"""
master_overview_map.png ile AYNI (guncel, dogru) veriyi kullanir - ayni havza
sinirlari, alt-havzalar, dere agi, dolin noktalari, istasyonlar - ama eski
("saydam") gorsel stiliyle cizer: alt-havza dolgu renkleri daha seffaf, DEM
zemin (hillshade/terrain) daha belirgin gorunur, boylece renklerin altindan
arazi rolyefi okunabilir.

Girdi: DEM_UTM35N.tif, subbasins_by_station.shp, streams_by_station.shp,
       basin_official.shp, reservoir_official.shp,
       hidro_meteoroloji/processed/doline_points_enriched.csv,
       hidro_meteoroloji/processed/station_locations.csv
Cikti: dem_processed/master_overview_map_transparent.png

Calistirma: proje kok dizininden
    .venv\\Scripts\\python.exe scripts\\plot_master_overview_map_transparent.py
"""
import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import rasterio

OUT_DIR = "dem_processed"
PROCESSED_DIR = "hidro_meteoroloji/processed"

COLORS = {"FP1": "#e6194B", "FP2": "#3cb44b", "FP3": "#4363d8", "RG6": "#f58231"}

with rasterio.open(f"{OUT_DIR}/DEM_UTM35N.tif") as src:
    dem = src.read(1)
    dem = np.ma.masked_equal(dem, src.nodata)
    extent = [src.bounds.left, src.bounds.right, src.bounds.bottom, src.bounds.top]

subbasins = gpd.read_file(f"{OUT_DIR}/subbasins_by_station.shp")
streams = gpd.read_file(f"{OUT_DIR}/streams_by_station.shp")
basin_official = gpd.read_file(f"{OUT_DIR}/basin_official.shp")
reservoir = gpd.read_file(f"{OUT_DIR}/reservoir_official.shp")
dolines = pd.read_csv(f"{PROCESSED_DIR}/doline_points_enriched.csv")
stations = pd.read_csv(f"{PROCESSED_DIR}/station_locations.csv")

priority_order = ["FP1", "FP2", "FP3", "RG6"]

fig, ax = plt.subplots(figsize=(14, 11))
# DEM zemini eski versiyondaki gibi belirgin (dusuk alpha yerine yuksek gorunurluk)
ax.imshow(dem, cmap="gist_earth", extent=extent, alpha=0.9)

unassigned_streams = streams[streams.station_id.isna()]
unassigned_streams.plot(ax=ax, color="gray", linewidth=0.3, alpha=0.5, zorder=2)

for sid in priority_order:
    sub = subbasins[subbasins.station_id == sid]
    area = sub.iloc[0]["area_km2"] if len(sub) else 0
    # Seffaf dolgu: arazi rolyefi renklerin altindan gorunsun
    sub.plot(ax=ax, facecolor=COLORS[sid], edgecolor="black", linewidth=1.5, alpha=0.45,
              label=f"{sid} ({area:.1f} km2)", zorder=3)

basin_official.boundary.plot(ax=ax, color="black", linewidth=2, linestyle=(0, (3, 2)),
                              label="Resmi havza siniri", zorder=5)
reservoir.plot(ax=ax, facecolor="#4cc9f0", edgecolor="#1f78b4", linewidth=1, alpha=0.85, zorder=4)

for sid in priority_order:
    sub = streams[streams.station_id == sid]
    sub.plot(ax=ax, color="black", linewidth=0.8, alpha=0.7, zorder=4)

# dolinler: kind'a gore renk (subbasin renkleriyle ayni), karst override -> ucgen
for sid in priority_order:
    sub = dolines[dolines.station_id == sid]
    normal = sub[sub.assign_method != "karst_magara_override"]
    karst = sub[sub.assign_method == "karst_magara_override"]
    ax.scatter(normal.x_utm35n, normal.y_utm35n, color=COLORS[sid], edgecolor="black",
               marker="o", s=45, zorder=6)
    ax.scatter(karst.x_utm35n, karst.y_utm35n, color=COLORS[sid], edgecolor="black",
               marker="^", s=70, zorder=6, label="Dolin (karst/magara override)" if sid == "RG6" else None)
unassigned = dolines[dolines.station_id.isna()]
ax.scatter(unassigned.x_utm35n, unassigned.y_utm35n, color="lightgray", edgecolor="black",
           marker="o", s=45, zorder=6, label="Dolin - atanmamis (harici)")

stations_plot = stations[stations.station_id.isin(
    priority_order + ["RG7", "RG8", "RG9", "RG10", "RG11", "RG12"])]
ax.scatter(stations_plot["x_utm35n"], stations_plot["y_utm35n"], c="white", s=60,
           edgecolor="black", marker="s", zorder=7)
for _, row in stations_plot.iterrows():
    ax.annotate(row["station_id"], (row["x_utm35n"], row["y_utm35n"]), fontsize=8,
                fontweight="bold", xytext=(3, 3), textcoords="offset points", zorder=7)

ax.set_title("Yuvacik Havzasi - Tam Ozet Haritasi (guncel veri, saydam stil)\n"
              "3 Debi Istasyonu (FP1-FP3) + RG6 Baraj Yerel Havzasi + Dere Agi + Dolinler")
ax.set_xlabel("UTM 35N Easting (m)")
ax.set_ylabel("UTM 35N Northing (m)")
ax.legend(loc="lower left", fontsize=7, framealpha=0.9)
ax.set_aspect("equal")
plt.tight_layout()

out_png = f"{OUT_DIR}/master_overview_map_transparent.png"
plt.savefig(out_png, dpi=150)
print(f"Yazildi: {out_png}")
