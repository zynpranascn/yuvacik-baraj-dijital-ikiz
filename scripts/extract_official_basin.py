"""
raw_data/geotiff_renders/hava_dis_sinir_geo.tif (= "havza dis sinir", gercek GeoTIFF, EPSG:32635 ile
zaten georeferansli) dosyasindan havza sinirini, raw_data/geotiff_renders/rezervuar_geo.tif'ten
rezervuar poligonunu, raw_data/geotiff_renders/dere_geo.tif'ten dere agini dogrudan gercek piksel
koordinatlarindan cikartir. JPEG digitize etmekten cok daha kesin, cunku bu
dosyalarin CRS/transform'u dosyanin kendi icinde tanimli.

Cikti: dem_processed/
    - basin_official.shp       (havza siniri, raw_data/geotiff_renders/hava_dis_sinir_geo.tif'ten)
    - reservoir_official.shp   (raw_data/geotiff_renders/rezervuar_geo.tif'ten)
    - streams_official.shp     (raw_data/geotiff_renders/dere_geo.tif'ten)
"""
import geopandas as gpd
import numpy as np
import rasterio
from rasterio.features import shapes
from scipy.ndimage import binary_dilation, binary_erosion, binary_fill_holes, label
from shapely.geometry import shape

OUT_DIR = "dem_processed"


def extract_line_feature_as_polygon(path, out_path, color_check):
    """Bir cizgi katmanini (RGB render) tespit edip en buyuk kapali alanı
    poligon olarak cikartir (havza siniri gibi kapali egriler icin)."""
    with rasterio.open(path) as src:
        bands = [src.read(b) for b in range(1, src.count + 1)]
        transform = src.transform
        crs = src.crs
        nodata = src.nodata

    mask = color_check(bands)
    print(f"  {path}: {mask.sum()} eslesen piksel")

    filled_mask = None
    for n_iter in [1, 2, 3, 5, 8, 12, 16, 20]:
        dilated = binary_dilation(mask, iterations=n_iter) if n_iter > 0 else mask
        filled = binary_fill_holes(dilated)
        interior = filled & ~dilated
        labeled, n = label(interior)
        if n == 0:
            continue
        sizes = np.bincount(labeled.ravel())
        sizes[0] = 0
        if sizes.max() < 1000:
            continue
        largest_label = sizes.argmax()
        component = (labeled == largest_label) | dilated
        filled_mask = binary_erosion(component, iterations=n_iter) if n_iter > 0 else component
        print(f"  dilation={n_iter} ile basarili, ic alan piksel={( labeled == largest_label).sum()}")
        break
    if filled_mask is None:
        print("  UYARI: kapali alan bulunamadi, sadece cizgi maskesi kullanilacak")
        filled_mask = mask

    polygons = [
        shape(geom) for geom, val in shapes(filled_mask.astype("uint8"), mask=filled_mask, transform=transform)
        if val == 1
    ]
    polygons.sort(key=lambda p: p.area, reverse=True)
    poly = polygons[0]
    gdf = gpd.GeoDataFrame({"id": [1], "area_km2": [round(poly.area / 1e6, 2)]}, geometry=[poly], crs=crs)
    gdf.to_file(out_path)
    print(f"  Yazildi: {out_path}, alan={poly.area/1e6:.2f} km2")
    return poly


def extract_simple_mask_as_polygon(path, out_path, expected_xy=None, max_reasonable_km2=50):
    """NOT: bu GeoTIFF'lerin piksel icerigi cerceve/eksen/etiket gibi harita
    duzeni ogelerini de iceriyor (duz numpy->PNG onizlemesinde gorundu). Bu
    yuzden en buyuk poligonu degil, makul boyuttaki (max_reasonable_km2 altinda,
    VE rasterin tam sinirini kaplamayan -- cerceve degil) ve varsa beklenen
    konuma (expected_xy) en yakin poligonu seciyoruz."""
    with rasterio.open(path) as src:
        band = src.read(1)
        transform = src.transform
        crs = src.crs
        full_bounds = src.bounds
    mask = (band == 0)
    print(f"  {path}: {mask.sum()} dolu piksel")
    polygons = [
        shape(geom) for geom, val in shapes(mask.astype("uint8"), mask=mask, transform=transform)
        if val == 1
    ]
    full_w = full_bounds.right - full_bounds.left
    full_h = full_bounds.top - full_bounds.bottom
    candidates = []
    for p in polygons:
        minx, miny, maxx, maxy = p.bounds
        is_frame = (maxx - minx) > 0.7 * full_w and (maxy - miny) > 0.7 * full_h
        if is_frame:
            continue
        if p.area / 1e6 < max_reasonable_km2:
            candidates.append(p)
    if not candidates:
        candidates = polygons
    if expected_xy is not None:
        ex, ey = expected_xy
        candidates.sort(key=lambda p: p.centroid.distance(shape({"type": "Point", "coordinates": (ex, ey)})))
    else:
        candidates.sort(key=lambda p: p.area, reverse=True)
    poly = candidates[0]
    gdf = gpd.GeoDataFrame({"id": [1], "area_km2": [round(poly.area / 1e6, 2)]}, geometry=[poly], crs=crs)
    gdf.to_file(out_path)
    print(f"  Yazildi: {out_path}, alan={poly.area/1e6:.2f} km2")
    return poly


def main():
    print("Havza siniri cikartiliyor (raw_data/geotiff_renders/hava_dis_sinir_geo.tif)...")
    def basin_color(bands):
        r, g, b = bands
        return (r > 150) & (g < 120) & (b < 120)
    extract_line_feature_as_polygon("raw_data/geotiff_renders/hava_dis_sinir_geo.tif", f"{OUT_DIR}/basin_official.shp", basin_color)

    print("\nRezervuar cikartiliyor (raw_data/geotiff_renders/rezervuar_geo.tif)...")
    # RG6 (Dolusavak) barajin hemen yaninda -- beklenen konum olarak kullaniliyor
    extract_simple_mask_as_polygon(
        "raw_data/geotiff_renders/rezervuar_geo.tif", f"{OUT_DIR}/reservoir_official.shp",
        expected_xy=(751195.67, 4506848.0), max_reasonable_km2=20,
    )

    print("\nDere agi cikartiliyor (raw_data/geotiff_renders/dere_geo.tif)...")
    with rasterio.open("raw_data/geotiff_renders/dere_geo.tif") as src:
        bands = [src.read(b) for b in range(1, src.count + 1)]
        transform = src.transform
        crs = src.crs
    r, g, b = bands
    stream_mask = (b > 150) & (r < 150)
    from rasterio.features import shapes as rshapes
    geoms = [shape(geom) for geom, val in rshapes(stream_mask.astype("uint8"), mask=stream_mask, transform=transform) if val == 1]
    gdf = gpd.GeoDataFrame({"id": range(len(geoms))}, geometry=geoms, crs=crs)
    gdf.to_file(f"{OUT_DIR}/streams_official.shp")
    print(f"  Yazildi: {OUT_DIR}/streams_official.shp, {len(geoms)} parca")


if __name__ == "__main__":
    main()
