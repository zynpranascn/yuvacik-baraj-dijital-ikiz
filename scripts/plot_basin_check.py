"""Havza delineasyonu gorsel dogrulama: DEM + havza siniri + akis agi + istasyonlar + calisma_alani."""
import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import rasterio

with rasterio.open("dem_processed/DEM_UTM35N.tif") as src:
    dem = src.read(1)
    dem = np.ma.masked_equal(dem, src.nodata)
    extent = [src.bounds.left, src.bounds.right, src.bounds.bottom, src.bounds.top]

basin = gpd.read_file("dem_processed/basin_delineated.shp")
boundary = gpd.read_file("../../dem/dem/dem.gdb", layer="calisma_alani").to_crs("EPSG:32635")
ref_map = gpd.read_file("dem_processed/basin_delineated_reference_map.shp")
basin_official = gpd.read_file("dem_processed/basin_official.shp")
reservoir_official = gpd.read_file("dem_processed/reservoir_official.shp")
streams_official = gpd.read_file("dem_processed/streams_official.shp")
streams = gpd.read_file("dem_processed/stream_network.geojson")
streams = streams.set_crs("EPSG:32635", allow_override=True)  # geojson CRS metadata eksik
stations = pd.read_csv("hidro_meteoroloji/processed/station_locations.csv")

fig, ax = plt.subplots(figsize=(12, 10))
ax.imshow(dem, cmap="terrain", extent=extent, alpha=0.8)
boundary.boundary.plot(ax=ax, color="black", linewidth=1.5, linestyle="--", label="calisma_alani (kaynak sinir)")
basin.boundary.plot(ax=ax, color="red", linewidth=2, label="Hesaplanan havza (D8, basin_delineated)")
ref_map.boundary.plot(ax=ax, color="magenta", linewidth=2, linestyle=":", label="Referans harita (JPEG digitize)")
basin_official.boundary.plot(ax=ax, color="lime", linewidth=2, linestyle="-.", label="Resmi havza siniri (hava_dis_sinir_geo.tif)")
reservoir_official.plot(ax=ax, color="deepskyblue", zorder=6, label="Rezervuar (resmi)")
streams.plot(ax=ax, color="blue", linewidth=0.5, alpha=0.7)
ax.scatter(stations["x_utm35n"], stations["y_utm35n"], c="orange", s=60, edgecolor="black", zorder=5, label="Istasyonlar")
for _, row in stations.iterrows():
    ax.annotate(row["station_id"], (row["x_utm35n"], row["y_utm35n"]), fontsize=8, xytext=(3, 3), textcoords="offset points")

ax.set_title("Yuvacik Havzasi - Hesaplanan Havza Siniri vs calisma_alani")
ax.legend(loc="lower left")
ax.set_xlabel("UTM 35N Easting (m)")
ax.set_ylabel("UTM 35N Northing (m)")
ax.set_aspect("equal")
plt.tight_layout()
plt.savefig("dem_processed/basin_validation_map.png", dpi=150)
print("Kaydedildi: dem_processed/basin_validation_map.png")
