"""
Dolin (karst obruk) noktalarini DEM + resmi havza siniri + akis agi +
istasyonlarla ayni haritada birlestirir, ve her dolin noktasinin havza
icinde/disinda oldugunu + en yakin istasyon/dere mesafesini hesaplar.

basin_official.shp kullanilir (dem_flow_pipeline.py / extract_official_basin.py
sonrasi projede nihai referans olarak kabul edilen havza siniri, bkz.
scripts/station_fusion.py'deki ayni tercih).

Girdi:
    - dem_processed/DEM_UTM35N.tif
    - dem_processed/basin_official.shp
    - dem_processed/stream_network.geojson
    - raw_data/hidro_meteoroloji/processed/station_locations.csv
    - raw_data/hidro_meteoroloji/processed/doline_points.shp

Cikti:
    - dem_processed/dolines_basin_overlay.png   (gorsel)
    - raw_data/hidro_meteoroloji/processed/doline_points_enriched.csv
      (id, x, y, inside_basin, dist_to_nearest_stream_m, nearest_station_id, dist_to_station_m)

Calistirma: proje kok dizininden
    .venv\\Scripts\\python.exe scripts\\plot_dolines_on_basin.py
"""
import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import rasterio

DEM_PATH = "dem_processed/DEM_UTM35N.tif"
BASIN_PATH = "dem_processed/basin_official.shp"
STREAMS_PATH = "dem_processed/stream_network.geojson"
STATIONS_PATH = "raw_data/hidro_meteoroloji/processed/station_locations.csv"
DOLINES_PATH = "raw_data/hidro_meteoroloji/processed/doline_points.shp"
OUT_DIR = "dem_processed"
PROCESSED_DIR = "raw_data/hidro_meteoroloji/processed"


def main():
    with rasterio.open(DEM_PATH) as src:
        dem = src.read(1)
        dem = np.ma.masked_equal(dem, src.nodata)
        extent = [src.bounds.left, src.bounds.right, src.bounds.bottom, src.bounds.top]

    basin = gpd.read_file(BASIN_PATH)
    streams = gpd.read_file(STREAMS_PATH)
    streams = streams.set_crs("EPSG:32635", allow_override=True)
    stations = pd.read_csv(STATIONS_PATH)
    dolines = gpd.read_file(DOLINES_PATH)

    print(f"Dolin noktasi: {len(dolines)}")
    print(f"Havza (resmi referans): {basin.geometry.area.iloc[0] / 1e6:.2f} km2")

    # --- Havza icinde mi? ---
    basin_geom = basin.geometry.iloc[0]
    dolines["inside_basin"] = dolines.geometry.within(basin_geom)
    n_inside = dolines["inside_basin"].sum()
    print(f"Havza icinde: {n_inside} / {len(dolines)}  (disinda: {len(dolines) - n_inside})")

    # --- En yakin dere hattina mesafe (karstik sizinti/baglanti riski icin faydali) ---
    streams_union = streams.geometry.union_all()
    dolines["dist_to_nearest_stream_m"] = dolines.geometry.apply(
        lambda p: round(p.distance(streams_union), 1)
    )

    # --- En yakin istasyon ---
    station_xy = stations[["x_utm35n", "y_utm35n"]].to_numpy()
    nearest_station, nearest_dist = [], []
    for pt in dolines.geometry:
        d = np.linalg.norm(station_xy - np.array([pt.x, pt.y]), axis=1)
        i = d.argmin()
        nearest_station.append(stations.iloc[i]["station_id"])
        nearest_dist.append(round(float(d[i]), 1))
    dolines["nearest_station_id"] = nearest_station
    dolines["dist_to_station_m"] = nearest_dist

    out_csv = dolines.drop(columns="geometry")
    out_csv.to_csv(f"{PROCESSED_DIR}/doline_points_enriched.csv", index=False)
    print(f"Yazildi: {PROCESSED_DIR}/doline_points_enriched.csv")

    # --- Harita ---
    fig, ax = plt.subplots(figsize=(14, 11))
    ax.imshow(dem, cmap="terrain", extent=extent, alpha=0.85)
    basin.boundary.plot(ax=ax, color="lime", linewidth=2, linestyle="-.",
                         label=f"Resmi havza siniri ({basin_geom.area/1e6:.1f} km2)")
    streams.plot(ax=ax, color="blue", linewidth=0.6, alpha=0.7, label="Dere agi (D8)")

    inside = dolines[dolines["inside_basin"]]
    outside = dolines[~dolines["inside_basin"]]
    ax.scatter(inside.geometry.x, inside.geometry.y, s=70, c="red", edgecolor="black",
               zorder=6, label=f"Dolin - havza ici (n={len(inside)})")
    ax.scatter(outside.geometry.x, outside.geometry.y, s=70, c="orange", edgecolor="black",
               marker="^", zorder=6, label=f"Dolin - havza disi (n={len(outside)})")

    ax.scatter(stations["x_utm35n"], stations["y_utm35n"], c="white", s=50,
               edgecolor="black", marker="s", zorder=5, label="Istasyonlar")
    for _, row in stations.iterrows():
        ax.annotate(row["station_id"], (row["x_utm35n"], row["y_utm35n"]), fontsize=7,
                    xytext=(3, 3), textcoords="offset points")

    ax.set_title("Yuvacik Havzasi - Dolin (karst obruk) noktalari, havza siniri, dere agi, istasyonlar")
    ax.set_xlabel("UTM 35N Easting (m)")
    ax.set_ylabel("UTM 35N Northing (m)")
    ax.legend(loc="lower left", fontsize=8)
    ax.set_aspect("equal")
    plt.tight_layout()
    out_png = f"{OUT_DIR}/dolines_basin_overlay.png"
    plt.savefig(out_png, dpi=150)
    print(f"Yazildi: {out_png}")


if __name__ == "__main__":
    main()
