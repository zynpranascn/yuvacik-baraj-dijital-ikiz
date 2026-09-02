"""
Istasyon koordinat tablosunu (Isas_Barajlar Excel, 'Plan ve Koordinantlar' sayfasi)
temiz bir CSV + shapefile'a donusturur.

Kaynaktaki UTM WGS84 zon degerleri (34T/35T/36T) tutarsiz oldugu icin
kullanilmiyor; tum koordinatlar DMS lat/lon'dan EPSG:32635'e (UTM 35N)
pyproj ile yeniden projekte ediliyor -> proje genelindeki DEM/raster CRS'i
ile birebir tutarli olsun diye.

Yukseklik: kaynaktaki KOT degeri (elevation_kot_m, muhtemelen elle/GPS ile
olculmus) ile dem_processed/DEM_UTM35N.tif'ten (kontur haritasindan uretilen,
dogrulanmis gercek DEM) ayni koordinatta okunan deger (elevation_m)
karsilastirilir. Aralarinda daglik istasyonlarda 25-47m'ye varan fark var --
muhtemelen istasyon GPS koordinatinin egimli arazide yatay hatasi. Sonraki
hesaplarda (lapse rate vb.) DEM_UTM35N.tif'ten okunan elevation_m esas alinir,
KOT sadece referans icin elevation_kot_m olarak saklanir.
"""
import pandas as pd
import geopandas as gpd
import rasterio
from shapely.geometry import Point
from pyproj import Transformer

# Kaynak: raw_data/hidro_meteoroloji/raw/İSAŞ_Barajlar Günlük Veri_TEMMUZ_2026.xlsx
# "Plan ve Koordinantlar" sayfasi, YUVACIK BARAJI istasyon tablosu (elle dogrulanmis)
stations_raw = [
    # id,   name,        type,   lat_d, lat_m, lat_s,  lon_d, lon_m, lon_s, kot_m
    ("RG6",  "Dolusavak",   "Meteoroloji", 40, 40, 27.4, 29, 58, 19.3, 178),
    ("RG7",  "Tepecik",     "Yagis/Sicaklik", 40, 37, 38.1, 29, 59, 25.4, 700),
    ("RG8",  "Aytepe",      "Yagis/Sicaklik", 40, 36,  2.4, 29, 56,  8.4, 953),
    ("RG9",  "Kartepe",     "Yagis/Sicaklik", 40, 39, 21.0, 30,  5, 44.0, 1340),
    ("RG10", "Cilekli",     "Yagis/Sicaklik", 40, 32, 30.1, 30,  2, 38.7, 805),
    ("RG11", "Kazandere",   "Yagis/Sicaklik", 40, 37, 12.2, 29, 57,  8.4, 732),
    ("RG12", "Haci Osman",  "Yagis/Sicaklik", 40, 33,  1.0, 29, 49,  8.0, 865),
    ("FP1",  "Kirazdere",   "Akis/Debi",   40, 38, 33.1, 29, 56, 38.8, 185),
    ("FP2",  "Kazandere",   "Akis/Debi",   40, 38, 22.6, 29, 57, 37.9, 180),
    ("FP3",  "Serindere",   "Akis/Debi",   40, 38, 48.5, 30,  1,  1.1, 284),
]

def dms_to_dd(d, m, s):
    return d + m / 60 + s / 3600

rows = []
for sid, name, typ, la_d, la_m, la_s, lo_d, lo_m, lo_s, kot in stations_raw:
    lat = dms_to_dd(la_d, la_m, la_s)
    lon = dms_to_dd(lo_d, lo_m, lo_s)
    rows.append({"station_id": sid, "name": name, "type": typ,
                 "lat_wgs84": lat, "lon_wgs84": lon, "elevation_kot_m": kot})

df = pd.DataFrame(rows)

# WGS84 (EPSG:4326) -> UTM 35N (EPSG:32635)
transformer = Transformer.from_crs("EPSG:4326", "EPSG:32635", always_xy=True)
df["x_utm35n"], df["y_utm35n"] = transformer.transform(df["lon_wgs84"].values, df["lat_wgs84"].values)

# Yukseklik icin DEM_UTM35N.tif'ten okunan degeri esas al (kontur haritasindan
# uretilmis, dogrulanmis; istasyon GPS/KOT degerinden daha guvenilir)
with rasterio.open("dem_processed/DEM_UTM35N.tif") as dem_src:
    rows_idx, cols_idx = rasterio.transform.rowcol(
        dem_src.transform, df["x_utm35n"].values, df["y_utm35n"].values
    )
    band = dem_src.read(1)
    df["elevation_m"] = band[rows_idx, cols_idx]
df["elevation_diff_m"] = (df["elevation_m"] - df["elevation_kot_m"]).round(1)

out_dir = "raw_data/hidro_meteoroloji/processed"
import os
os.makedirs(out_dir, exist_ok=True)

csv_path = f"{out_dir}/station_locations.csv"
df.to_csv(csv_path, index=False)
print(f"CSV yazildi: {csv_path}")
print(df.to_string(index=False))

gdf = gpd.GeoDataFrame(
    df, geometry=[Point(xy) for xy in zip(df.x_utm35n, df.y_utm35n)], crs="EPSG:32635"
)
shp_path = f"{out_dir}/station_locations.shp"
gdf.to_file(shp_path)
print(f"\nShapefile yazildi: {shp_path}")
