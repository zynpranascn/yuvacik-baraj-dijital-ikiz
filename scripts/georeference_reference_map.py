"""
Referans harita gorselini (DemWhatsApp Image) georeferanslar ve kirmizi
"Basin Boundary" cizgisini EPSG:32635'e digitize eder.

DEM tabanli analizden (basin_delineated.shp) VAZGECMIYORUZ -- bu script sadece
ek/karsilastirma katmani olarak basin_delineated_reference_map.shp uretir.

Georeferans yontemi: gorseldeki lat/lon eksen etiketlerinin piksel merkezleri
tespit edilip (bkz. tick-detection asamasi), pikselden lon/lat'e dogrusal
(affine) donusum kuruluyor. Boylece harita "ImageJ ile elle tiklama" gibi
gurultulu bir yontem yerine, etiketlerin gercek piksel konumlarina dayanarak
saglam bir sekilde georeferanslaniyor.

Calistirma: proje kok dizininden
    .venv\\Scripts\\python.exe scripts\\georeference_reference_map.py
"""
import geopandas as gpd
import numpy as np
from PIL import Image
from pyproj import Transformer
from scipy.ndimage import binary_dilation, binary_erosion, binary_fill_holes, label
from shapely.geometry import Polygon, LineString
from skimage import measure

IMG_PATH = "raw_data/hidro_meteoroloji/raw/DemWhatsApp Image 2026-08-10 at 12.05.49.jpeg"
OUT_DIR = "dem_processed"

# Tick-etiket piksel merkezleri (onceki adimda dark-pixel analiziyle bulundu)
LON_PIXELS = [180.5, 424.5, 665.0, 909.0]
LON_DEGS = [29 + 50/60, 29 + 55/60, 30.0, 30 + 5/60]
LAT_PIXELS = [70.5, 390.5, 710.5]
LAT_DEGS = [40 + 40/60, 40 + 35/60, 40 + 30/60]

# Legend kutusunun yaklasik piksel siniri (kirmizi 'Basin Boundary' ornek
# cizgisini havza siniri sanip karistirmamak icin bu bolge maskeleniyor)
LEGEND_BBOX = (870, 605, 1116, 780)  # x0,y0,x1,y1


def fit_linear(px, deg):
    a, b = np.polyfit(px, deg, 1)
    return a, b


def main():
    img = Image.open(IMG_PATH)
    arr = np.array(img)
    h, w = arr.shape[:2]
    print(f"Gorsel boyutu: {w}x{h}")

    a_lon, b_lon = fit_linear(LON_PIXELS, LON_DEGS)
    a_lat, b_lat = fit_linear(LAT_PIXELS, LAT_DEGS)
    print(f"lon(x) = {a_lon:.8f}*x + {b_lon:.4f}")
    print(f"lat(y) = {a_lat:.8f}*y + {b_lat:.4f}")

    def pixel_to_lonlat(x, y):
        return a_lon * x + b_lon, a_lat * y + b_lat

    # --- Kirmizi "Basin Boundary" cizgisini tespit et ---
    r, g, b = arr[:, :, 0].astype(int), arr[:, :, 1].astype(int), arr[:, :, 2].astype(int)
    red_mask = (r > 170) & (g < 90) & (b < 90)

    lx0, ly0, lx1, ly1 = LEGEND_BBOX
    red_mask[ly0:ly1, lx0:lx1] = False

    print(f"Kirmizi piksel sayisi: {red_mask.sum()}")

    # Ince kirmizi cizgideki JPEG sikistirma kaynakli bosluklari kapatmak icin
    # giderek artan miktarda genislet (dilate), her denemede doldurmayi (fill)
    # dene; basarili olunca ayni miktarda asindirarak (erode) siniri geri
    # yaklasik orijinaline getir.
    filled_mask = None
    for n_iter in [2, 4, 6, 8, 10, 12, 16, 20, 25, 30]:
        dilated = binary_dilation(red_mask, iterations=n_iter)
        filled = binary_fill_holes(dilated)
        interior = filled & ~dilated
        labeled, n = label(interior)
        if n == 0:
            continue
        sizes = np.bincount(labeled.ravel())
        sizes[0] = 0
        if sizes.max() < 5000:  # cok kucuk bir bosluk kapatilmis olabilir, gercek havza degil
            continue
        largest_label = sizes.argmax()
        component = (labeled == largest_label) | dilated
        filled_mask = binary_erosion(component, iterations=n_iter)
        print(f"  Basarili: dilation={n_iter} piksel, ic alan={(labeled == largest_label).sum()} piksel")
        break
    if filled_mask is None:
        raise RuntimeError("Kirmizi cizgi hicbir dilation seviyesinde kapali alan olusturmadi")

    # --- Kontur (ordered boundary path) cikar -- artik doldurulmus/dolu alan ---
    contours = measure.find_contours(filled_mask.astype(float), level=0.5)
    contours.sort(key=len, reverse=True)
    main_contour = contours[0]  # (row, col) ciftleri, sirali
    print(f"En buyuk kontur: {len(main_contour)} nokta")

    # (row, col) -> (x_pixel, y_pixel) -> lon/lat -> EPSG:32635
    transformer = Transformer.from_crs("EPSG:4326", "EPSG:32635", always_xy=True)
    lonlat = [pixel_to_lonlat(col, row) for row, col in main_contour]
    lons, lats = zip(*lonlat)
    xs, ys = transformer.transform(np.array(lons), np.array(lats))
    coords = list(zip(xs, ys))

    line = LineString(coords)
    poly = Polygon(coords)
    if not poly.is_valid:
        poly = poly.buffer(0)
    area_km2 = poly.area / 1_000_000
    print(f"Digitize edilen sinir alani: {area_km2:.2f} km2")

    gdf = gpd.GeoDataFrame({"id": [1], "area_km2": [round(area_km2, 2)]}, geometry=[poly], crs="EPSG:32635")
    gdf.to_file(f"{OUT_DIR}/basin_delineated_reference_map.shp")
    print(f"Yazildi: {OUT_DIR}/basin_delineated_reference_map.shp")

    line_gdf = gpd.GeoDataFrame({"id": [1]}, geometry=[line], crs="EPSG:32635")
    line_gdf.to_file(f"{OUT_DIR}/basin_boundary_reference_line.shp")
    print(f"Yazildi: {OUT_DIR}/basin_boundary_reference_line.shp")


if __name__ == "__main__":
    main()
