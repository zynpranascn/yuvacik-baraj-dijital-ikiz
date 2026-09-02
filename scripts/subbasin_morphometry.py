"""
Faz 1 - Alt-havza/dere analizini derinlestirme: her alt-havza (FP1-FP3, RG6)
icin morfometrik parametreler + drenaj yogunlugu + (varsa) yagis istatistigi.
(RG8 ve RG12 artik ayri alt-havza degil, ikisi de FP1'e union edildi.)

Hesaplanan parametreler (standart havza morfometrisi, orn. Horton/Schumm/
Strahler literaturu):
    - Alan, cevre uzunlugu
    - Yukseklik: min/mean/max, rolyef (max-min)
    - Egim: ortalama egim (derece), DEM'den np.gradient ile
    - Kompaktlik orani (Gravelius Kc): P / (2*sqrt(pi*A)) -- 1.0 = daire (en
      "kompakt"), buyudukce havza uzuyor/dallaniyor -> taskinin zirveye
      ulasmasi daha yavas/yayili olur
    - Uzama orani (Schumm Re): (2/Lb) * sqrt(A/pi) -- 1.0 = daire, kuculdukce
      havza uzar (Lb = havzanin en uzun ekseni)
    - Drenaj yogunlugu Dd = toplam dere uzunlugu / alan (km/km2) -- yuksek
      Dd = hizli/sivri taskin tepkisi, dusuk Dd = yavas/emici (gecirgen,
      ormanli) havza

Girdi: dem_processed/subbasins/<sid>/{dem.tif, streams.shp, boundary.shp}
       hidro_meteoroloji/processed/hydro_met_clean.csv

Cikti: dem_processed/subbasin_morphometry.csv (+ konsola siralanmis rapor)

Calistirma: proje kok dizininden
    .venv\\Scripts\\python.exe scripts\\subbasin_morphometry.py
"""
import numpy as np
import geopandas as gpd
import pandas as pd
import rasterio

SUB_DIR = "dem_processed/subbasins"
SUBBASIN_IDS = ["FP1", "FP2", "FP3", "RG6"]
NAME_MAP = {"FP1": "Kirazdere", "FP2": "Kazandere", "FP3": "Serindere",
            "RG6": "Baraj Yerel Havzasi"}
PRECIP_COL = {"RG6": "precip_tepecik_mm", "FP2": "precip_kazandere_mm",
              "FP3": "precip_cilekli_mm", "FP1": "precip_aytepe_mm"}


def slope_stats(dem_path):
    with rasterio.open(dem_path) as src:
        dem = src.read(1).astype(float)
        nodata = src.nodata
        res = src.res[0]  # metre/piksel (kare piksel varsayimi)
    dem[dem == nodata] = np.nan
    dy, dx = np.gradient(dem, res)
    slope_deg = np.degrees(np.arctan(np.sqrt(dx**2 + dy**2)))
    valid = ~np.isnan(slope_deg)
    elev_valid = dem[~np.isnan(dem)]
    return {
        "elev_min_m": float(np.nanmin(dem)), "elev_mean_m": float(np.nanmean(dem)),
        "elev_max_m": float(np.nanmax(dem)),
        "slope_mean_deg": float(np.nanmean(slope_deg[valid])),
        "slope_max_deg": float(np.nanmax(slope_deg[valid])),
    }


def main():
    hydro_met = pd.read_csv("hidro_meteoroloji/processed/hydro_met_clean.csv")

    rows = []
    for sid in SUBBASIN_IDS:
        boundary = gpd.read_file(f"{SUB_DIR}/{sid}/boundary.shp").geometry.iloc[0]
        area_km2 = boundary.area / 1e6
        perimeter_km = boundary.length / 1000

        minx, miny, maxx, maxy = boundary.bounds
        Lb_km = max(maxx - minx, maxy - miny) / 1000  # kaba en-uzun-eksen tahmini

        gravelius_kc = perimeter_km / (2 * np.sqrt(np.pi * area_km2))
        schumm_re = (2 / Lb_km) * np.sqrt(area_km2 / np.pi) if Lb_km > 0 else np.nan

        try:
            streams = gpd.read_file(f"{SUB_DIR}/{sid}/streams.shp")
            stream_length_km = streams.geometry.length.sum() / 1000
        except Exception:
            stream_length_km = 0.0
        drainage_density = stream_length_km / area_km2 if area_km2 > 0 else np.nan

        sstats = slope_stats(f"{SUB_DIR}/{sid}/dem.tif")

        precip_col = PRECIP_COL.get(sid)
        mean_annual_precip_mm = None
        if precip_col and precip_col in hydro_met.columns:
            daily_mean = hydro_met[precip_col].mean()
            mean_annual_precip_mm = round(daily_mean * 365.25, 1)

        rows.append({
            "station_id": sid, "dere_adi": NAME_MAP[sid],
            "alan_km2": round(area_km2, 2), "cevre_km": round(perimeter_km, 2),
            "drenaj_yogunlugu_km_km2": round(drainage_density, 3),
            "kompaktlik_Kc": round(gravelius_kc, 3),
            "uzama_orani_Re": round(schumm_re, 3),
            "yukseklik_min_m": round(sstats["elev_min_m"], 1),
            "yukseklik_ort_m": round(sstats["elev_mean_m"], 1),
            "yukseklik_max_m": round(sstats["elev_max_m"], 1),
            "rolyef_m": round(sstats["elev_max_m"] - sstats["elev_min_m"], 1),
            "egim_ort_derece": round(sstats["slope_mean_deg"], 1),
            "egim_max_derece": round(sstats["slope_max_deg"], 1),
            "yillik_ort_yagis_mm": mean_annual_precip_mm,
        })

    df = pd.DataFrame(rows).sort_values("alan_km2", ascending=False).reset_index(drop=True)
    df.to_csv("dem_processed/subbasin_morphometry.csv", index=False)
    print(f"Yazildi: dem_processed/subbasin_morphometry.csv\n")

    print("=== Alan (uzundan kisaya) ===")
    print(df[["dere_adi", "alan_km2", "cevre_km"]].to_string(index=False))

    print("\n=== Drenaj yogunlugu (yuksekten alcaga -- hizli/sivri tepki riski) ===")
    print(df.sort_values("drenaj_yogunlugu_km_km2", ascending=False)
          [["dere_adi", "drenaj_yogunlugu_km_km2"]].to_string(index=False))

    print("\n=== Ortalama egim (yuksekten alcaga -- erozyon/hizli akis riski) ===")
    print(df.sort_values("egim_ort_derece", ascending=False)
          [["dere_adi", "egim_ort_derece", "egim_max_derece"]].to_string(index=False))

    print("\n=== Rolyef (yukseklik farki, yuksekten alcaga) ===")
    print(df.sort_values("rolyef_m", ascending=False)
          [["dere_adi", "yukseklik_min_m", "yukseklik_max_m", "rolyef_m"]].to_string(index=False))

    print("\n=== Kompaktlik orani (Kc, 1.0'a yakin=dairesel/hizli tepki, "
          "buyudukce=uzun/yayili tepki) ===")
    print(df.sort_values("kompaktlik_Kc")
          [["dere_adi", "kompaktlik_Kc", "uzama_orani_Re"]].to_string(index=False))

    have_precip = df[df.yillik_ort_yagis_mm.notna()]
    if len(have_precip):
        print("\n=== Yillik ortalama yagis (mm, sadece kendi met istasyonu olan alt-havzalar) ===")
        print(have_precip.sort_values("yillik_ort_yagis_mm", ascending=False)
              [["dere_adi", "yillik_ort_yagis_mm"]].to_string(index=False))


if __name__ == "__main__":
    main()
