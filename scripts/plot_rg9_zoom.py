"""RG9 civarinda DEM + akis biriktirme + havza sinirlari + akis izini yakinlastirip gosterir."""
import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import rasterio

DIRMAP = (64, 128, 1, 2, 4, 8, 16, 32)
OFFSETS = {64: (-1,0), 128: (-1,1), 1: (0,1), 2: (1,1), 4: (1,0), 8: (1,-1), 16: (0,-1), 32: (-1,-1)}

with rasterio.open("dem_processed/DEM_UTM35N.tif") as src:
    dem = src.read(1)
    dem = np.ma.masked_equal(dem, src.nodata)
    extent = [src.bounds.left, src.bounds.right, src.bounds.bottom, src.bounds.top]

with rasterio.open("dem_processed/flow_direction_D8.tif") as src:
    fdir = src.read(1)
    transform = src.transform
    fshape = fdir.shape

basin = gpd.read_file("dem_processed/basin_delineated.shp")
official = gpd.read_file("dem_processed/basin_official.shp")

def trace(x, y, max_steps=500):
    row, col = rasterio.transform.rowcol(transform, x, y)
    pts = [(x, y)]
    for _ in range(max_steps):
        if not (0 <= row < fshape[0] and 0 <= col < fshape[1]):
            break
        d = fdir[row, col]
        if d == 0 or d not in OFFSETS:
            break
        dr, dc = OFFSETS[d]
        row, col = row + dr, col + dc
        px, py = rasterio.transform.xy(transform, row, col)
        pts.append((px, py))
    return pts

rg9 = (761709.917347, 4505161.0)
paths = [trace(*rg9), trace(761709.9, 4504661.0), trace(761209.9, 4505161.0), trace(760709.9, 4504161.0)]

fig, ax = plt.subplots(figsize=(10, 9))
ax.imshow(dem, cmap="terrain", extent=extent, alpha=0.85)
basin.boundary.plot(ax=ax, color="red", linewidth=2, label="Hesaplanan havza (D8)")
official.boundary.plot(ax=ax, color="lime", linewidth=2, linestyle="-.", label="Resmi havza siniri")

for i, path in enumerate(paths):
    xs, ys = zip(*path)
    ax.plot(xs, ys, color="black", linewidth=1.5, alpha=0.8, label="Akis izi (RG9 -> ?)" if i == 0 else None)
    ax.plot(xs[0], ys[0], marker="o", color="yellow", markeredgecolor="black", markersize=8, zorder=6)

ax.annotate("RG9", rg9, fontsize=11, fontweight="bold", xytext=(8, 8), textcoords="offset points")
ax.set_xlim(755000, 772461)
ax.set_ylim(4498000, 4513138)
ax.set_title("RG9 Civari: Akis izleri kuzeye (haritanin disina) cikiyor")
ax.legend(loc="lower right")
ax.set_xlabel("UTM 35N Easting (m)")
ax.set_ylabel("UTM 35N Northing (m)")
ax.set_aspect("equal")
plt.tight_layout()
plt.savefig("dem_processed/rg9_flow_trace.png", dpi=150)
print("Kaydedildi: dem_processed/rg9_flow_trace.png")
