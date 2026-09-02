"""
FP4 istasyonunu ekler -- Isas Excel'inin "Plan ve Koordinantlar" sayfasindaki
10 istasyon tablosunda (RG6-RG12, FP1-FP3) FP4 YOK, ama AYNI Excel dosyasina
GOMULU bir referans harita gorselinde (xl/media/image1.png, kullanicinin
paylastigi ekran goruntusuyle ayni) FP4 acikca isaretli -- ana govde (RG6/baraj)
uzerinde, FP1'in yukarisinda ayri bir debi olcum noktasi.

FP4'un DMS/KOT degeri hicbir tabloda yok, bu yuzden koordinati SADECE bu
referans haritanin piksel konumundan digitize edilerek turetildi:
    1. Haritadaki 8 bilinen istasyonun (FP1,FP2,FP3,RG6,RG7,RG8,RG9,RG10)
       piksel konumu renk/sekil tespitiyle bulundu (bkz. arastirma asamasi).
    2. Bu 8 nokta icin piksel->UTM35N afin donusumu en kucuk kareler ile
       kuruldu (kalan hata 18-125m, harita seması dogrulugu icin makul).
    3. FP4'un piksel konumu (549,167) bu donusumle gercek koordinata cevrildi.
    4. D8 ile snap+delineate edilip komsu istasyonlarla (FP2, FP3 -- ic ice;
       FP1 -- ayrik) tutarliligi dogrulandi (bkz. build_subbasins_fp1_4.py).

NOT: Bu koordinat GPS/DMS olcumu DEGIL, harita digitizasyonu -- tahmini hata
±100-150m mertebesinde. Diger istasyonlardan farkli olarak elevation_kot_m
yok (kaynakta yok), source='map_digitized' ile isaretlendi.

Cikti: raw_data/hidro_meteoroloji/processed/station_locations.csv + .shp (FP4 satiri eklenir)
"""
import geopandas as gpd
import pandas as pd
import rasterio
from pyproj import Transformer
from shapely.geometry import Point

STATIONS_PATH = "raw_data/hidro_meteoroloji/processed/station_locations.csv"
DEM_PATH = "dem_processed/DEM_UTM35N.tif"

# Piksel GCP'leri (reference_basin_map_from_excel.png, 1211x851) -- bilinen
# istasyonlarin gercek koordinatlariyla eslestirilip afin fit edildi.
FP4_PIXEL = (549, 167)
FP4_UTM = (750019.0, 4504223.7)  # fit sonucu (bkz. docstring)


def main():
    df = pd.read_csv(STATIONS_PATH)
    if "FP4" in df["station_id"].values:
        print("FP4 zaten mevcut, cikiliyor.")
        return

    x, y = FP4_UTM
    transformer = Transformer.from_crs("EPSG:32635", "EPSG:4326", always_xy=True)
    lon, lat = transformer.transform(x, y)

    with rasterio.open(DEM_PATH) as dem_src:
        r, c = rasterio.transform.rowcol(dem_src.transform, x, y)
        elev = float(dem_src.read(1)[r, c])

    new_row = {
        "station_id": "FP4", "name": "Ana Govde (referans haritadan digitize)",
        "type": "Akis/Debi", "lat_wgs84": lat, "lon_wgs84": lon,
        "elevation_kot_m": None, "x_utm35n": x, "y_utm35n": y,
        "elevation_m": elev, "elevation_diff_m": None,
        "source": "map_digitized",
    }
    if "source" not in df.columns:
        df["source"] = "isas_table_dms"
    df = pd.concat([df, pd.DataFrame([new_row])], ignore_index=True)
    df.to_csv(STATIONS_PATH, index=False)
    print(f"FP4 eklendi: x={x:.1f} y={y:.1f} elev={elev:.1f}m")

    gdf = gpd.GeoDataFrame(
        df, geometry=[Point(xy) for xy in zip(df.x_utm35n, df.y_utm35n)], crs="EPSG:32635"
    )
    gdf.to_file(STATIONS_PATH.replace(".csv", ".shp"))
    print(f"Yazildi: {STATIONS_PATH} (+.shp)")


if __name__ == "__main__":
    main()
