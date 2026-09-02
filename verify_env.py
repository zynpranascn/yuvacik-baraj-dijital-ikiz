"""Faz 1 ortam dogrulama scripti."""
import importlib
import sys

packages = [
    "numpy", "pandas", "scipy", "rasterio", "fiona", "shapely",
    "geopandas", "pyproj", "pysheds", "sklearn", "statsmodels",
    "openpyxl", "pyarrow", "matplotlib", "psycopg2", "sqlalchemy",
]

print(f"Python: {sys.version}\n")
ok = True
for pkg in packages:
    try:
        mod = importlib.import_module(pkg)
        ver = getattr(mod, "__version__", "?")
        print(f"[OK]   {pkg:<15} {ver}")
    except Exception as e:
        ok = False
        print(f"[FAIL] {pkg:<15} {e}")

print("\n--- DEM dosyasi okuma testi (dem_processed/DEM_UTM35N.tif) ---")
import rasterio
with rasterio.open("dem_processed/DEM_UTM35N.tif") as src:
    print("CRS:      ", src.crs)
    print("Boyut:    ", src.width, "x", src.height)
    print("Bounds:   ", src.bounds)
    print("NoData:   ", src.nodata)
    band = src.read(1)
    print("Min/Max:  ", float(band.min()), float(band.max()))

print("\nSonuc:", "TUM PAKETLER HAZIR" if ok else "BAZI PAKETLER EKSIK")
