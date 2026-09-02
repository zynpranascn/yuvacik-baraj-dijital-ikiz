"""
Dolin (karst obruk) nokta konumlarini raw_data/hidro_meteoroloji/raw/dolin/dolin_geo.tif
render'indan cikartir.

Kaynak dosyanin (diger *_geo.tif render'leri gibi) gomulu CRS'i yok ama .tfw'si
projedeki paylasilan 30m/piksel ızgarayla (ayni origin: 730131.57, 4513138.38)
birebir ayni -- bu yuzden EPSG:32635 varsayimi guvenli (bkz. extract_official_basin.py,
georeference_reference_map.py'deki daha once dogrulanan ayni-izgara tespiti).

Yontem: gri-tonlamali hillshade zemin uzerine cizilmis kucuk siyah (RGB<60)
sembolleri tespit et. Ayirt edilmesi gereken iki "gurultu" kaynagi var:
  1. Hillshade'in kendi golge bolgeleri de koyu -- ama bunlar TEK, cok buyuk
     (binlerce piksel) bir baglantili bilesen olusturuyor; kucuk boyut esigiyle
     (<500px) disarida birakiliyor.
  2. Harita cercevesindeki koordinat "tick" isaretleri de siyah, kenarlara
     yakin (<30px) ince/uzun bilesenler -- kenar payi (margin) ile eleniyor.
Kalan kucuk bilesenler kumeleniyor (6px yaricap) cunku tek bir obruk sembolu
anti-aliasing nedeniyle birden fazla ayri kucuk siyah parcaya bolunebiliyor.

NOT: Bu otomatik/piksel-tabanli bir cikarim -- kaynakta okunabilir bir lejant/
öznitelik tablosu olmadigi icin nokta sayisi ve konumlari yaklasiktir (+-1
piksel = 30m). dolin_validation_map.png ile gorsel dogrulama onerilir.

Cikti: raw_data/hidro_meteoroloji/processed/
    - doline_points.csv    (id, x_utm35n, y_utm35n)
    - doline_points.shp    (+ .dbf/.shx/.prj)
    - dolin_validation_map.png (dolin_geo.tif uzerine cikarilan noktalar)

Calistirma: proje kok dizininden
    .venv\\Scripts\\python.exe scripts\\extract_doline_points.py
"""
import os

import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import rasterio
from scipy.ndimage import center_of_mass, label
from shapely.geometry import Point

SRC_PATH = "raw_data/hidro_meteoroloji/raw/dolin/dolin_geo.tif"
OUT_DIR = "raw_data/hidro_meteoroloji/processed"
CRS = "EPSG:32635"

DARK_THRESHOLD = 60      # RGB < bu deger -> "siyah sembol" adayi
MAX_SYMBOL_PX = 500      # bundan buyuk bilesenler hillshade golgesi (havza cok daha buyuk)
BORDER_MARGIN_PX = 100   # harita cercevesi/koordinat tick alani (sol/sag kenar
                          # etiketleri 30px'te yetersiz kaldi -- gercek dolin
                          # kumeleri kenardan >450px iceride, guvenli bir pay)
CLUSTER_RADIUS_PX = 6    # ayni sembolun parcalanmis siyah bloklarini birlestirmek icin


def main():
    os.makedirs(OUT_DIR, exist_ok=True)

    with rasterio.open(SRC_PATH) as src:
        arr = src.read()
        transform = src.transform
        h, w = src.height, src.width

    r, g, b = arr[0].astype(int), arr[1].astype(int), arr[2].astype(int)
    dark = (r < DARK_THRESHOLD) & (g < DARK_THRESHOLD) & (b < DARK_THRESHOLD)

    labeled, n = label(dark)
    sizes = np.bincount(labeled.ravel())
    sizes[0] = 0
    centers = center_of_mass(dark, labeled, range(1, n + 1))

    candidates = []
    for i, (cy, cx) in enumerate(centers, start=1):
        if sizes[i] > MAX_SYMBOL_PX:
            continue
        if cy < BORDER_MARGIN_PX or cy > h - BORDER_MARGIN_PX or \
           cx < BORDER_MARGIN_PX or cx > w - BORDER_MARGIN_PX:
            continue
        candidates.append((cy, cx))
    print(f"Aday siyah sembol parcasi: {len(candidates)} (havza-govde golgesi ve "
          f"cerceve tick'leri elendikten sonra)")

    # Ayni sembolun parcalanmis parcalarini kumele (basit tek-gecis, mesafe<esik)
    candidates = np.array(candidates)
    used = np.zeros(len(candidates), dtype=bool)
    pixel_points = []
    for i in range(len(candidates)):
        if used[i]:
            continue
        d = np.linalg.norm(candidates - candidates[i], axis=1)
        group = d < CLUSTER_RADIUS_PX
        used |= group
        pixel_points.append(candidates[group].mean(axis=0))
    print(f"Kumelenmis dolin noktasi: {len(pixel_points)}")

    rows = []
    xs, ys = [], []
    for idx, (py, px) in enumerate(pixel_points, start=1):
        x, y = transform * (px, py)
        xs.append(x)
        ys.append(y)
        rows.append({"id": idx, "x_utm35n": x, "y_utm35n": y})

    gdf = gpd.GeoDataFrame(
        rows, geometry=[Point(x, y) for x, y in zip(xs, ys)], crs=CRS
    )
    gdf.drop(columns="geometry").to_csv(f"{OUT_DIR}/doline_points.csv", index=False)
    gdf.to_file(f"{OUT_DIR}/doline_points.shp")
    print(f"Yazildi: {OUT_DIR}/doline_points.csv, {OUT_DIR}/doline_points.shp")

    # --- Gorsel dogrulama ---
    fig, ax = plt.subplots(figsize=(14, 10))
    img = np.transpose(np.clip(arr, 0, 255).astype("uint8"), (1, 2, 0))
    extent = [transform.c, transform.c + w * transform.a,
              transform.f + h * transform.e, transform.f]
    ax.imshow(img, extent=extent)
    ax.scatter(xs, ys, s=60, facecolors="none", edgecolors="red", linewidths=1.5,
               label=f"Cikarilan dolin noktalari (n={len(xs)})")
    ax.set_title("Dolin noktalari - otomatik piksel cikarimi (dogrulama gerekli)")
    ax.set_xlabel("UTM 35N Easting (m)")
    ax.set_ylabel("UTM 35N Northing (m)")
    ax.legend(loc="lower left")
    ax.set_aspect("equal")
    plt.tight_layout()
    plt.savefig(f"{OUT_DIR}/dolin_validation_map.png", dpi=150)
    print(f"Yazildi: {OUT_DIR}/dolin_validation_map.png (gorsel dogrulama icin)")


if __name__ == "__main__":
    main()
