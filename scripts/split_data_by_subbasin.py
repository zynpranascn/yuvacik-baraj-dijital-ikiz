"""
Faz 1 - Alt-havzalara gore fiziksel veri bolme.

assign_streams_dolines_to_stations.py'nin ciktisi olan 4 alt-havza catchment'i
(subbasins_by_station.shp: FP1, FP2, FP3, RG6 -- FP4 YOK/kaldirildi, RG8 ve
RG12 artik ayri alt-havza degil, ikisi de FP1'e union edildi, bkz.
assign_streams_dolines_to_stations.py "REVIZYON NOTU" + RG8_MERGED_INTO) icin
HER BIRINE ait ayri bir klasor uretir; DEM/akis rasterlerini o alt-havzaya
kirpar, dere agini ve dolin noktalarini filtreler, ve o alt-havza sinirlari
icinde kalan meteoroloji istasyonlarinin hydro_met_clean.csv sutunlarini ayirir.

Hangi alt-havzada hangi meteoroloji istasyonu var (mekansal kesisim ile
dogrulandi, bkz. arastirma asamasi):
    FP1 -> RG8 (Aytepe yagis/sicaklik, artik FP1 icinde) -- RG12'nin de kendi
           met sutunu yok, sadece RG8'in verisi kullanilir
    FP2 -> RG11 (Kazandere yagis/sicaklik, FP2'nin yaninda)
    RG6 -> RG7 (Tepecik) + barajin kendi rezervuar/su dengesi sutunlari
           (RG6 = baraj, bu yuzden reservoir_level_m, inflow_m3, vb. de dahil;
           RG6 alani artik FP2+FP3'un GECTIGI + kendi yerel havzasi)
    FP3 -> RG10 (Cilekli)

Cikti: dem_processed/subbasins/<station_id>/
    - boundary.shp              (o alt-havzanin siniri)
    - dem.tif                   (DEM_UTM35N.tif kirpilmis)
    - flow_direction_D8.tif     (kirpilmis)
    - flow_accumulation.tif     (kirpilmis)
    - streams.shp               (o alt-havzaya ait dere segmentleri)
    - dolines.csv / dolines.shp (o alt-havzaya ait dolin noktalari)
    - stations.csv              (o alt-havza sinirlari icindeki istasyonlar)
    - hydro_met.csv             (varsa: ilgili istasyon(lar)in gunluk serisi)
    - metadata.json             (ozet: alan, istasyon sayisi, dere/dolin sayisi)

Calistirma: proje kok dizininden
    .venv\\Scripts\\python.exe scripts\\split_data_by_subbasin.py
"""
import json
import os

import geopandas as gpd
import pandas as pd
import rasterio
from rasterio.mask import mask
from shapely.geometry import Point

OUT_ROOT = "dem_processed/subbasins"
DEM_PROCESSED = "dem_processed"
PROCESSED_DIR = "hidro_meteoroloji/processed"

# Alt-havza -> o alt-havzada fiziksel olarak bulunan meteoroloji istasyonu
# (mekansal kesisimle dogrulandi). None = o alt-havzada RG-tipi istasyon yok.
MET_STATION_BY_SUBBASIN = {
    "FP2": "kazandere",
    "RG6": "tepecik", "FP1": "aytepe", "FP3": "cilekli",
}
# Baraj/rezervuar sutunlari -- RG6 (baraj) hangi alt-havzanin sinirlari
# icindeyse (mekansal kesisimle belirlenir, hardcoded degil) o alt-havzaya
# dahil edilir. RG6 zaten kendi alt-havzasi oldugu icin bu her zaman RG6.
RESERVOIR_COLS = [
    "reservoir_level_m", "useful_volume_mcm", "fill_pct", "outflow_treatment_m3",
    "inflow_m3", "outflow_spillway_m3", "sapanca_transfer_m3", "evaporation_m3",
    "inflow_q_m3s",
]


def clip_raster(src_path, geom, out_path):
    with rasterio.open(src_path) as src:
        out_arr, out_transform = mask(src, [geom], crop=True, nodata=src.nodata)
        profile = src.profile.copy()
        profile.update({
            "height": out_arr.shape[1], "width": out_arr.shape[2], "transform": out_transform,
        })
    with rasterio.open(out_path, "w", **profile) as dst:
        dst.write(out_arr)


def main():
    subs = gpd.read_file(f"{DEM_PROCESSED}/subbasins_by_station.shp")
    streams = gpd.read_file(f"{DEM_PROCESSED}/streams_by_station.shp")
    dolines = gpd.read_file(f"{PROCESSED_DIR}/doline_points_enriched.shp")
    stations = pd.read_csv(f"{PROCESSED_DIR}/station_locations.csv")
    hydro_met = pd.read_csv(f"{PROCESSED_DIR}/hydro_met_clean.csv")

    os.makedirs(OUT_ROOT, exist_ok=True)

    for _, srow in subs.iterrows():
        sid = srow["station_id"]
        name = srow["name"]
        geom = srow.geometry
        out_dir = f"{OUT_ROOT}/{sid}"
        os.makedirs(out_dir, exist_ok=True)
        print(f"\n=== {sid} ({name}) -> {out_dir} ===")

        # --- Sinir ---
        gpd.GeoDataFrame([srow.drop("geometry")], geometry=[geom], crs=subs.crs).to_file(
            f"{out_dir}/boundary.shp"
        )

        # --- Rasterlar ---
        for raster_name, out_name in [
            ("DEM_UTM35N.tif", "dem.tif"),
            ("flow_direction_D8.tif", "flow_direction_D8.tif"),
            ("flow_accumulation.tif", "flow_accumulation.tif"),
        ]:
            try:
                clip_raster(f"{DEM_PROCESSED}/{raster_name}", geom, f"{out_dir}/{out_name}")
            except ValueError as e:
                print(f"  UYARI: {raster_name} kirpma basarisiz ({e}) -- atlaniyor")
        print(f"  Rasterlar kirpildi: dem.tif, flow_direction_D8.tif, flow_accumulation.tif")

        # --- Dere agi ---
        sub_streams = streams[streams["station_id"] == sid]
        if len(sub_streams) > 0:
            sub_streams.to_file(f"{out_dir}/streams.shp")
        print(f"  Dere segmenti: {len(sub_streams)}")

        # --- Dolinler ---
        sub_dolines = dolines[dolines["station_id"] == sid]
        if len(sub_dolines) > 0:
            sub_dolines.to_file(f"{out_dir}/dolines.shp")
            sub_dolines.drop(columns="geometry").to_csv(f"{out_dir}/dolines.csv", index=False)
        print(f"  Dolin noktasi: {len(sub_dolines)}")

        # --- Icindeki istasyonlar ---
        inside_mask = stations.apply(
            lambda r: geom.contains(Point(r["x_utm35n"], r["y_utm35n"])), axis=1
        )
        sub_stations = stations[inside_mask]
        sub_stations.to_csv(f"{out_dir}/stations.csv", index=False)
        print(f"  Icindeki istasyonlar: {list(sub_stations['station_id'])}")

        # --- Hidro-meteoroloji zaman serisi (varsa) ---
        met_key = MET_STATION_BY_SUBBASIN.get(sid)
        cols = ["date", "qa_flag", "notes"]
        if met_key:
            cols += [c for c in hydro_met.columns if c.startswith(f"precip_{met_key}")
                     or c.startswith(f"temp_{met_key}") or c.startswith(f"snow_{met_key}")]
        if "RG6" in list(sub_stations["station_id"]):
            cols += RESERVOIR_COLS
        cols = [c for c in cols if c in hydro_met.columns]
        has_met_data = len(cols) > 3  # date/qa_flag/notes disinda bir sey var mi
        if has_met_data:
            hydro_met[cols].to_csv(f"{out_dir}/hydro_met.csv", index=False)
            print(f"  Hidro-met sutunlari: {[c for c in cols if c not in ('date','qa_flag','notes')]}")
        else:
            print("  Hidro-met verisi yok (bu alt-havzada olculen istasyon sutunu bulunmuyor)")

        # --- Metadata ---
        meta = {
            "station_id": sid, "name": name, "kind": srow["kind"],
            "area_km2": float(srow["area_km2"]),
            "n_stream_segments": int(len(sub_streams)),
            "n_dolines": int(len(sub_dolines)),
            "stations_inside": list(sub_stations["station_id"]),
            "has_met_timeseries": has_met_data,
            "met_source_station": met_key,
        }
        with open(f"{out_dir}/metadata.json", "w", encoding="utf-8") as f:
            json.dump(meta, f, indent=2, ensure_ascii=False)

    print(f"\nTum alt-havzalar yazildi: {OUT_ROOT}/<station_id>/")


if __name__ == "__main__":
    main()
