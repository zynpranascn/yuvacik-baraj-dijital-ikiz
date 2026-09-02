"""
Faz 1 - Asama 2: DEM Isleme Mimarisi (D8 akis yonu, akis biriktirme, havza delineasyonu)

richdem yerine pysheds kullanilir (Windows'ta derleyici gerektirmez, ayni D8
mantigini uygular -- bkz. ortam kurulum notlari).

Girdi: dem_processed/DEM_UTM35N.tif (kontur verisinden uretilen, dogrulanmis gercek DEM)

Cikti: dem_processed/
    - flow_direction_D8.tif   (pysheds/ESRI D8 kodlamasi: 1,2,4,8,16,32,64,128)
    - flow_accumulation.tif   (uint32, hucre sayisi)
    - basin_delineated.shp    (+ .dbf/.shx/.prj)
    - stream_network.geojson  (akis agi, threshold=1000 hucre)
    - validation_report.json

Calistirma: proje kok dizininden
    .venv\\Scripts\\python.exe scripts\\dem_flow_pipeline.py
"""
import json
import os

import geopandas as gpd
import numpy as np
import pandas as pd

# pysheds 0.5, NumPy 2.0'da kaldirilan np.in1d'i cagiriyor (yerine np.isin geldi).
# NumPy'i dusurmek diger paketlerle (scipy/numba/rasterio) ABI riski tasidigindan,
# burada hafif bir uyumluluk yamasi kullaniliyor.
if not hasattr(np, "in1d"):
    np.in1d = np.isin

import rasterio
from pysheds.grid import Grid
from rasterio.features import shapes
from shapely.geometry import shape
from shapely.ops import unary_union

DEM_PATH = "dem_processed/DEM_UTM35N.tif"
BOUNDARY_GDB = "../../dem/dem/dem.gdb"
OUT_DIR = "dem_processed"
STREAM_THRESHOLD = 1000  # hucre -- roadmap 2.2 adim 5
DIRMAP = (64, 128, 1, 2, 4, 8, 16, 32)  # pysheds/ESRI D8 kodlamasi


def resolve_invalid_directions(fdir_arr, dem_arr, dirmap, nodata_dem):
    """DIRMAP disi (gecersiz) yonlu hucreleri, konditone DEM'deki en dik inis
    komsusuna gore duzeltir. Bunlari 0 (nodata) yapmak yerine gercek bir yon
    atamak, akis agindaki (ve dolayisiyla havza delineasyonundaki) yapay
    kopukluklari onler -- bkz. bati kolu baglanti teshisi.

    NOT: pysheds'in resolve_flats() adimi NoData hucrelerini bile hafifce
    kaydiriyor (ör. -9999.0 -> -9998.99998), bu yuzden nodata karsilastirmasi
    tam esitlik (==) yerine bir esik degeriyle yapiliyor. Gercek yukseklikler
    (25-1575m) hep pozitif oldugundan -1000 guvenli bir esik.
    """
    NODATA_THRESHOLD = -1000.0

    # dirmap sirasi: N, NE, E, SE, S, SW, W, NW (pysheds/ESRI konvansiyonu)
    offsets = [(-1, 0), (-1, 1), (0, 1), (1, 1), (1, 0), (1, -1), (0, -1), (-1, -1)]
    code_offset = list(zip(dirmap, offsets))

    invalid_mask = ~np.isin(fdir_arr, dirmap)
    rows, cols = np.where(invalid_mask)
    n_rows, n_cols = dem_arr.shape
    fixed = 0
    for r, c in zip(rows, cols):
        z0 = dem_arr[r, c]
        if z0 < NODATA_THRESHOLD:
            continue
        best_code, best_drop = 0, -np.inf
        for code, (dr, dc) in code_offset:
            nr, nc = r + dr, c + dc
            if not (0 <= nr < n_rows and 0 <= nc < n_cols):
                continue
            zn = dem_arr[nr, nc]
            if zn < NODATA_THRESHOLD:
                continue
            drop = z0 - zn
            if drop > best_drop:
                best_drop, best_code = drop, code
        if best_drop > 0:
            fdir_arr[r, c] = best_code
            fixed += 1
    return fixed, len(rows) - fixed


def snap_outlet(acc_arr, transform, approx_x, approx_y, margin=5, win=15):
    """Yaklasik bir nokta civarinda en yuksek akis biriktirmeli hucreyi bulur
    (pour-point snapping), rasterin kenarindan en az `margin` piksel icerde
    kalacak sekilde."""
    approx_row, approx_col = rasterio.transform.rowcol(transform, approx_x, approx_y)
    r0, r1 = max(margin, approx_row - win), min(acc_arr.shape[0] - margin, approx_row + win)
    c0, c1 = max(margin, approx_col - win), min(acc_arr.shape[1] - margin, approx_col + win)
    window = acc_arr[r0:r1, c0:c1]
    local_r, local_c = np.unravel_index(np.argmax(window), window.shape)
    row, col = r0 + local_r, c0 + local_c
    x, y = transform * (col, row)
    return row, col, x, y


def delineate_polygon(grid, fdir, dirmap, x, y):
    """Verilen pour point'ten catchment hesaplar, en buyuk parcayi poligon
    olarak dondurur."""
    catchment = grid.catchment(x=x, y=y, fdir=fdir, dirmap=dirmap, xytype="coordinate")
    catchment_arr = np.array(catchment).astype("uint8")
    polygons = [
        shape(geom) for geom, val in shapes(catchment_arr, mask=catchment_arr.astype(bool), transform=grid.affine)
        if val == 1
    ]
    polygons.sort(key=lambda p: p.area, reverse=True)
    return polygons[0]


def write_geotiff(path, array, transform, crs, dtype, nodata):
    profile = {
        "driver": "GTiff", "dtype": dtype, "nodata": nodata,
        "width": array.shape[1], "height": array.shape[0], "count": 1,
        "crs": crs, "transform": transform, "compress": "lzw",
    }
    with rasterio.open(path, "w", **profile) as dst:
        dst.write(array.astype(dtype), 1)


def main():
    print("DEM yukleniyor...")
    grid = Grid.from_raster(DEM_PATH)
    dem = grid.read_raster(DEM_PATH)
    print(f"  Boyut: {dem.shape}, CRS: {grid.crs}")

    print("\nCukur/depresyon giderme (pit fill -> depression fill -> flat resolve)...")
    pit_filled = grid.fill_pits(dem)
    flooded = grid.fill_depressions(pit_filled)
    inflated = grid.resolve_flats(flooded)

    print("D8 akis yonu hesaplaniyor...")
    fdir = grid.flowdir(inflated, dirmap=DIRMAP)

    # DIRMAP disindaki degerler (NoData sinirinda cozulememis hucreler, ör. -1/-2)
    # ONCEDEN 0 (nodata) yapiliyordu -- bu, akis agini kopararak havzanin bir
    # kismini (bati kolu) yanlislikla disarida birakiyordu (D8 yol izleme ile
    # dogrulandi). Bunun yerine konditone DEM'e gore en dik inis komsusuna
    # yonlendirerek gercek bir yon atiyoruz, boylece akis sureklilikte kalir.
    # BU DUZELTME ACCUMULATION HESAPLANMADAN ONCE yapilmali.
    dem_for_resolve = np.array(inflated)
    fdir_arr = np.array(fdir)
    fixed, unresolved = resolve_invalid_directions(fdir_arr, dem_for_resolve, DIRMAP, grid.nodata)
    fdir[:, :] = fdir_arr
    print(f"  {fixed} gecersiz yonlu hucre en dik inise gore duzeltildi, "
          f"{unresolved} hucre cozulemedi (nodata komsu / DEM nodata)")
    if unresolved > 0:
        fdir[~np.isin(np.array(fdir), DIRMAP)] = 0

    print("Akis biriktirme (flow accumulation) hesaplaniyor...")
    acc = grid.accumulation(fdir, dirmap=DIRMAP)

    write_geotiff(f"{OUT_DIR}/flow_direction_D8.tif", np.array(fdir).astype("uint8"), grid.affine, grid.crs, "uint8", 0)
    write_geotiff(f"{OUT_DIR}/flow_accumulation.tif", np.array(acc), grid.affine, grid.crs, "uint32", 0)
    print(f"  Yazildi: flow_direction_D8.tif, flow_accumulation.tif")
    print(f"  Flow dir degerleri (temizlenmis): {sorted(np.unique(np.array(fdir)))}")
    print(f"  Flow accum araligi: {int(np.array(acc).min())} - {int(np.array(acc).max())}")

    print("\nHavza outlet'i (pour point) tespit ediliyor: RG6 (Dolusavak/baraj) istasyonuna en yakin akis hattina yapistiriliyor...")
    acc_arr = np.array(acc)
    stations = pd.read_csv("hidro_meteoroloji/processed/station_locations.csv")
    dam_row = stations[stations["station_id"] == "RG6"].iloc[0]

    outlet_row, outlet_col, outlet_x, outlet_y = snap_outlet(
        acc_arr, grid.affine, dam_row["x_utm35n"], dam_row["y_utm35n"]
    )
    print(f"  Yapistirilan outlet: ({outlet_row},{outlet_col}) -> koordinat ({outlet_x:.1f}, {outlet_y:.1f})")
    print(f"  Outlet akis biriktirme degeri: {int(acc_arr[outlet_row, outlet_col])} hucre")

    print("\nAna (dogal, tek-cikisli) havza delineasyonu hesaplaniyor...")
    natural_basin_geom = delineate_polygon(grid, fdir, DIRMAP, outlet_x, outlet_y)
    natural_area_km2 = natural_basin_geom.area / 1_000_000
    gpd.GeoDataFrame({"id": [1]}, geometry=[natural_basin_geom], crs=grid.crs).to_file(
        f"{OUT_DIR}/basin_delineated_natural.shp"
    )
    print(f"  Dogal havza alani: {natural_area_km2:.2f} km2 (referans: basin_delineated_natural.shp)")

    # --- Havza disinda kalan istasyonlar icin alt-havzalari hesapla ve birlestir ---
    # Yuvacik, dogal havzasi disindan da (Sapanca vb.) tunelle su transferi aliyor
    # (bkz. 'Sapanca'dan Temin Edilen Su' sutunu). Bu istasyonlarin bulundugu
    # yerel alt-havzalari ayri ayri delineate edip ana havzayla birlestirerek
    # toplam katkı alanini (referans haritadaki 'Basin Boundary'ye daha yakin)
    # yaklasik olarak temsil ediyoruz.
    print("\nHavza disinda kalan istasyonlar tespit ediliyor...")
    outside_rows = []
    for _, srow in stations.iterrows():
        pt_inside = natural_basin_geom.contains(
            shape({"type": "Point", "coordinates": (srow["x_utm35n"], srow["y_utm35n"])})
        )
        if not pt_inside and srow["station_id"] != "RG6":
            outside_rows.append(srow)
    print(f"  Havza disinda kalan istasyon sayisi: {len(outside_rows)}: "
          f"{[r['station_id'] for r in outside_rows]}")

    # Resmi havza siniri (hava_dis_sinir_geo.tif'ten cikarilan) -- her alt-havza
    # adayini buna karsi DOGRULUYORUZ. Sadece D8 hesaplamasi (spekulasyon) ile
    # eklemek yerine, gercek resmi veriyle kesisimi olculuyor; kesisimi dusuk
    # olan (orn. yanlis yone akan) adaylar reddediliyor.
    official_path = f"{OUT_DIR}/basin_official.shp"
    official_geom = None
    if os.path.exists(official_path):
        official_geom = gpd.read_file(official_path).geometry.iloc[0]
        print(f"  Resmi havza siniri yuklendi (dogrulama icin): {official_geom.area/1e6:.2f} km2")
    else:
        print("  UYARI: basin_official.shp bulunamadi, alt-havzalar dogrulanamayacak "
              "(scripts/extract_official_basin.py once calistirilmali)")

    OVERLAP_THRESHOLD = 0.5  # resmi sinirla en az %50 kesisim sart

    sub_polygons = []
    sub_details = []
    for srow in outside_rows:
        try:
            sub_row, sub_col, sub_x, sub_y = snap_outlet(
                acc_arr, grid.affine, srow["x_utm35n"], srow["y_utm35n"], win=15
            )
            sub_geom = delineate_polygon(grid, fdir, DIRMAP, sub_x, sub_y)
            sub_area = sub_geom.area / 1_000_000

            if official_geom is not None:
                overlap_km2 = sub_geom.intersection(official_geom).area / 1_000_000
                overlap_frac = overlap_km2 / sub_area if sub_area > 0 else 0
                accepted = overlap_frac >= OVERLAP_THRESHOLD
            else:
                overlap_frac = None
                accepted = True  # dogrulama verisi yoksa eskisi gibi kabul et

            status = "KABUL" if accepted else "REDDEDILDI (resmi sinirla ortusmuyor)"
            ov_str = f", resmi kesisim=%{overlap_frac*100:.1f}" if overlap_frac is not None else ""
            print(f"  {srow['station_id']} ({srow['name']}): alt-havza {sub_area:.2f} km2{ov_str} -> {status}")

            sub_details.append({
                "station_id": srow["station_id"], "name": srow["name"],
                "sub_catchment_area_km2": round(sub_area, 2),
                "official_overlap_pct": round(overlap_frac * 100, 1) if overlap_frac is not None else None,
                "accepted": accepted,
            })
            if accepted:
                sub_polygons.append(sub_geom)
        except Exception as e:
            print(f"  {srow['station_id']}: delineasyon basarisiz ({e}), atlaniyor")

    print("\nAna havza + dogrulanmis alt-havzalar birlestiriliyor (union)...")
    combined_geom = unary_union([natural_basin_geom] + sub_polygons)
    basin_gdf = gpd.GeoDataFrame({"id": [1]}, geometry=[combined_geom], crs=grid.crs)
    basin_gdf.to_file(f"{OUT_DIR}/basin_delineated.shp")
    basin_area_km2 = combined_geom.area / 1_000_000
    print(f"  Birlesik (toplam katki) havza alani: {basin_area_km2:.2f} km2")
    print(f"  Yazildi: {OUT_DIR}/basin_delineated.shp")

    print("\nAkis agi (stream network) cikartiliyor...")
    branches = grid.extract_river_network(fdir, acc > STREAM_THRESHOLD, dirmap=DIRMAP)
    with open(f"{OUT_DIR}/stream_network.geojson", "w", encoding="utf-8") as f:
        json.dump(branches, f)
    n_branches = len(branches.get("features", []))
    print(f"  {n_branches} kol, yazildi: {OUT_DIR}/stream_network.geojson")

    print("\n--- Dogrulama: calisma_alani (kaynak havza siniri) ile karsilastirma ---")
    boundary = gpd.read_file(BOUNDARY_GDB, layer="calisma_alani").to_crs(grid.crs)
    ref_area_km2 = boundary.geometry.area.iloc[0] / 1_000_000
    area_diff_pct = (basin_area_km2 - ref_area_km2) / ref_area_km2 * 100
    print(f"  Referans (calisma_alani) alani: {ref_area_km2:.2f} km2")
    print(f"  Hesaplanan havza alani:         {basin_area_km2:.2f} km2")
    print(f"  Fark: {area_diff_pct:+.1f}%")

    dem_arr = np.array(dem)
    dem_valid = dem_arr[dem_arr > -1000]
    report = {
        "dem": {
            "shape": list(dem.shape),
            "elevation_min_m": float(dem_valid.min()),
            "elevation_max_m": float(dem_valid.max()),
        },
        "flow_direction": {"unique_values": [int(v) for v in sorted(np.unique(np.array(fdir)))]},
        "flow_accumulation": {
            "min": int(acc_arr.min()), "max": int(acc_arr.max()),
        },
        "outlet": {"x": float(outlet_x), "y": float(outlet_y), "accumulation": int(acc_arr[outlet_row, outlet_col])},
        "basin": {
            "combined_area_km2": round(basin_area_km2, 2),
            "natural_single_outlet_area_km2": round(natural_area_km2, 2),
            "sub_catchments": sub_details,
            "reference_calisma_alani_area_km2": round(ref_area_km2, 2),
            "area_diff_pct": round(area_diff_pct, 1),
            "note": (
                "Dogal (tek cikisli, RG6/baraj outlet'inden D8 ile hesaplanan) havza "
                f"{natural_area_km2:.2f} km2 -- havza disi istasyonlarin (RG8/RG9/RG12/FP1) "
                "bulundugu bolgeleri dislar. Her biri icin ayri D8 alt-havzasi hesaplanip, "
                "SADECE resmi havza siniriyla (basin_official.shp, hava_dis_sinir_geo.tif'ten) "
                "en az %50 kesisenler ana havzaya birlestirildi (union). Ornegin RG9'un "
                "alt-havzasi resmi sinirla %0 kesisiyordu (o bolge dogal olarak kuzeye, "
                "haritanin disina akiyor -- D8 yol izlemeyle dogrulandi) ve REDDEDILDI; "
                "RG8 (%100), RG12 (%100), FP1 (%98.4) kabul edildi. Detaylar icin "
                "sub_catchments listesindeki 'accepted'/'official_overlap_pct' alanlarina "
                "bakin. 'basin_delineated.shp' artik bu dogrulanmis birlesik alani temsil "
                "ediyor; sadece dogal havza icin 'basin_delineated_natural.shp' ayrica "
                "saklandi. calisma_alani (~1043 km2) havza siniri degil, cok daha genis bir "
                "etut alani (basin_official ile karistirilmamali)."
            ),
        },
        "stream_network": {"threshold_cells": STREAM_THRESHOLD, "n_branches": n_branches},
    }
    with open(f"{OUT_DIR}/validation_report.json", "w", encoding="utf-8") as f:
        json.dump(report, f, indent=2, ensure_ascii=False)
    print(f"\nRapor yazildi: {OUT_DIR}/validation_report.json")


if __name__ == "__main__":
    main()
