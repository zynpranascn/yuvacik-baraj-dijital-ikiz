"""
Faz 1 - Gercek DEM uretimi: kontur cizgilerinden (dem.gdb/izhops) TIN-benzeri
Delaunay lineer interpolasyonla raster DEM insa eder, ve ayni bolgeyi
dem_geo.tif'in gri tonlarindan tersine-cevirerek uretilen yaklasik DEM ile
karsilastirir.

Kaynak: ../../dem/dem/dem.gdb  (katmanlar: izhops [kontur], calisma_alani [sinir])
        CRS: ED_1950_Lambert_Conformal_Conic -> EPSG:32635'e cevrilir

Cikti: dem_processed/
    - DEM_UTM35N.tif                 (konturdan, gercek yukseklik, float32)
    - DEM_from_render_reversed.tif   (renkten tersine cevirme, sadece karsilastirma icin)
    - dem_comparison_report.json

Calistirma: proje kok dizininden
    .venv\\Scripts\\python.exe scripts\\build_dem_from_contours.py
"""
import json
import os

import geopandas as gpd
import numpy as np
import rasterio
from rasterio.features import geometry_mask
from rasterio.transform import from_origin
from scipy.interpolate import LinearNDInterpolator

GDB_PATH = "../../dem/dem/dem.gdb"
RENDER_TIF = "dem_geo.tif"
OUT_DIR = "dem_processed"
RESOLUTION = 20  # metre -- 25m kontur araligina uygun, makul islem suresi
TARGET_CRS = "EPSG:32635"
HOLDOUT_FRACTION = 0.1  # dogrulama icin ayrilan kontur cizgisi orani


def load_contours_and_boundary():
    contours = gpd.read_file(GDB_PATH, layer="izhops").to_crs(TARGET_CRS)
    boundary = gpd.read_file(GDB_PATH, layer="calisma_alani").to_crs(TARGET_CRS)
    return contours, boundary


def lines_to_points(contours_gdf):
    xs, ys, zs = [], [], []
    for elev, geom in zip(contours_gdf["Contour"], contours_gdf.geometry):
        parts = geom.geoms if geom.geom_type == "MultiLineString" else [geom]
        for line in parts:
            for x, y in line.coords:
                xs.append(x)
                ys.append(y)
                zs.append(elev)
    return np.array(xs), np.array(ys), np.array(zs)


def build_grid(bounds, resolution):
    left, bottom, right, top = bounds
    width = int(np.ceil((right - left) / resolution))
    height = int(np.ceil((top - bottom) / resolution))
    transform = from_origin(left, top, resolution, resolution)
    xs = left + (np.arange(width) + 0.5) * resolution
    ys = top - (np.arange(height) + 0.5) * resolution
    grid_x, grid_y = np.meshgrid(xs, ys)
    return grid_x, grid_y, width, height, transform


def write_geotiff(path, array, transform, nodata):
    profile = {
        "driver": "GTiff",
        "dtype": "float32",
        "nodata": nodata,
        "width": array.shape[1],
        "height": array.shape[0],
        "count": 1,
        "crs": TARGET_CRS,
        "transform": transform,
        "compress": "lzw",
    }
    with rasterio.open(path, "w", **profile) as dst:
        dst.write(array.astype("float32"), 1)


def interpolate_dem(points_xyz, boundary, resolution):
    """Kontur noktalarindan (x,y,z) Delaunay lineer interpolasyonla maskelenmis
    DEM raster'i uretir. Donus: (dem_masked, transform, width, height, NODATA)"""
    x, y, z = points_xyz
    interp = LinearNDInterpolator(np.column_stack([x, y]), z)

    bounds = boundary.total_bounds
    grid_x, grid_y, width, height, transform = build_grid(bounds, resolution)

    dem = interp(grid_x, grid_y)

    nan_mask = np.isnan(dem)
    if nan_mask.any():
        from scipy.interpolate import NearestNDInterpolator
        nn = NearestNDInterpolator(np.column_stack([x, y]), z)
        dem[nan_mask] = nn(grid_x[nan_mask], grid_y[nan_mask])
        print(f"  {nan_mask.sum()} piksel (convex hull disi) nearest-neighbor ile dolduruldu")

    mask = geometry_mask(boundary.geometry, out_shape=(height, width), transform=transform, invert=True)
    NODATA = -9999.0
    dem_masked = np.where(mask, dem, NODATA).astype("float32")
    return dem_masked, transform, width, height, NODATA


def main():
    os.makedirs(OUT_DIR, exist_ok=True)

    print("Konturlar ve sinir yukleniyor...")
    contours, boundary = load_contours_and_boundary()
    print(f"  {len(contours)} kontur cizgisi, yukseklik araligi "
          f"{contours['Contour'].min()}-{contours['Contour'].max()}m")

    # Holdout: dogrulama icin bazi kontur CIZGILERINI (nokta degil) disarida birak.
    # Bu SADECE durust bir hata raporu uretmek icin -- nihai DEM asagida TUM
    # konturlarla (holdout dahil) ayrica insa edilir.
    rng = np.random.default_rng(42)
    n_holdout = max(1, int(len(contours) * HOLDOUT_FRACTION))
    holdout_idx = rng.choice(contours.index, size=n_holdout, replace=False)
    train_contours = contours.drop(index=holdout_idx)
    holdout_contours = contours.loc[holdout_idx]
    print(f"  Egitim: {len(train_contours)} cizgi | Holdout (dogrulama): {len(holdout_contours)} cizgi")

    print("\n[1/2] Dogrulama DEM'i (egitim seti, holdout haric) insa ediliyor...")
    x, y, z = lines_to_points(train_contours)
    print(f"  {len(x)} nokta (egitim seti)")
    print("  Delaunay ucgenleme + lineer interpolator kuruluyor (biraz surebilir)...")
    val_dem, val_transform, val_w, val_h, NODATA = interpolate_dem((x, y, z), boundary, RESOLUTION)
    val_path = f"{OUT_DIR}/_validation_only_dem.tif"
    write_geotiff(val_path, val_dem, val_transform, NODATA)
    print(f"  Raster boyutu: {val_w} x {val_h} piksel ({RESOLUTION}m cozunurluk)")

    # --- Karsilastirma: renkten tersine cevirme ---
    print("\n--- Renkten tersine cevirme (karsilastirma icin) ---")
    with rasterio.open(RENDER_TIF) as src:
        gray = src.read(1).astype("float64")
        render_transform = src.transform
        render_nodata = src.nodata
        rows, cols = rasterio.transform.rowcol(render_transform, x, y)
        rows = np.clip(rows, 0, gray.shape[0] - 1)
        cols = np.clip(cols, 0, gray.shape[1] - 1)
        sampled_gray = gray[rows, cols]
        valid_sample = sampled_gray != render_nodata

        # gray = a*elev + b  (en kucuk kareler)
        A = np.vstack([z[valid_sample], np.ones(valid_sample.sum())]).T
        a, b = np.linalg.lstsq(A, sampled_gray[valid_sample], rcond=None)[0]
        print(f"  Uydurulan dogrusal iliski: gray = {a:.4f} * elev + {b:.2f}")

        elev_est = (gray - b) / a
        elev_est[gray == render_nodata] = NODATA
        write_geotiff(f"{OUT_DIR}/DEM_from_render_reversed.tif", elev_est, render_transform, NODATA)
        print(f"  Yazildi: {OUT_DIR}/DEM_from_render_reversed.tif")

    # --- Dogrulama: holdout kontur noktalarinda iki yontemi karsilastir ---
    print("\n--- Holdout dogrulama (egitimde kullanilmayan kontur noktalari) ---")
    hx, hy, hz = lines_to_points(holdout_contours)

    with rasterio.open(val_path) as dem_src:
        h_rows, h_cols = rasterio.transform.rowcol(dem_src.transform, hx, hy)
        h_rows = np.clip(h_rows, 0, dem_src.height - 1)
        h_cols = np.clip(h_cols, 0, dem_src.width - 1)
        dem_band = dem_src.read(1)
        pred_contour_dem = dem_band[h_rows, h_cols]

    with rasterio.open(RENDER_TIF) as src:
        gray_h_rows, gray_h_cols = rasterio.transform.rowcol(src.transform, hx, hy)
        gray_h_rows = np.clip(gray_h_rows, 0, src.height - 1)
        gray_h_cols = np.clip(gray_h_cols, 0, src.width - 1)
        gray_h = src.read(1).astype("float64")[gray_h_rows, gray_h_cols]
        pred_color_dem = (gray_h - b) / a

    valid_h = (pred_contour_dem != NODATA)
    err_contour = pred_contour_dem[valid_h] - hz[valid_h]
    err_color = pred_color_dem[valid_h] - hz[valid_h]

    report = {
        "resolution_m": RESOLUTION,
        "n_holdout_points": int(valid_h.sum()),
        "contour_based_dem": {
            "MAE_m": float(np.mean(np.abs(err_contour))),
            "RMSE_m": float(np.sqrt(np.mean(err_contour ** 2))),
            "max_abs_error_m": float(np.max(np.abs(err_contour))),
        },
        "color_reversed_dem": {
            "MAE_m": float(np.mean(np.abs(err_color))),
            "RMSE_m": float(np.sqrt(np.mean(err_color ** 2))),
            "max_abs_error_m": float(np.max(np.abs(err_color))),
            "fit_gray_to_elev": {"a": float(a), "b": float(b)},
        },
    }
    with open(f"{OUT_DIR}/dem_comparison_report.json", "w", encoding="utf-8") as f:
        json.dump(report, f, indent=2, ensure_ascii=False)

    print(f"\nKontur-DEM  -> MAE={report['contour_based_dem']['MAE_m']:.2f}m  "
          f"RMSE={report['contour_based_dem']['RMSE_m']:.2f}m  "
          f"max={report['contour_based_dem']['max_abs_error_m']:.2f}m")
    print(f"Renk-DEM    -> MAE={report['color_reversed_dem']['MAE_m']:.2f}m  "
          f"RMSE={report['color_reversed_dem']['RMSE_m']:.2f}m  "
          f"max={report['color_reversed_dem']['max_abs_error_m']:.2f}m")
    print(f"\nRapor yazildi: {OUT_DIR}/dem_comparison_report.json")

    # --- Nihai uretim DEM'i: TUM konturlarla (holdout dahil) yeniden insa et ---
    print("\n[2/2] Nihai DEM tum kontur verisiyle (1199 cizgi) insa ediliyor...")
    x_all, y_all, z_all = lines_to_points(contours)
    print(f"  {len(x_all)} nokta (tam veri seti)")
    final_dem, final_transform, final_w, final_h, _ = interpolate_dem(
        (x_all, y_all, z_all), boundary, RESOLUTION
    )
    out_path = f"{OUT_DIR}/DEM_UTM35N.tif"
    write_geotiff(out_path, final_dem, final_transform, NODATA)
    valid = final_dem[final_dem != NODATA]
    print(f"  Nihai DEM yazildi: {out_path}")
    print(f"  Yukseklik araligi: {valid.min():.1f} - {valid.max():.1f} m, ortalama {valid.mean():.1f} m")

    os.remove(val_path)
    print(f"  ({val_path} silindi -- sadece dogrulama icin gecici dosyaydi)")


if __name__ == "__main__":
    main()
