"""
Faz 1 - Asama B: Istasyon Entegrasyonu & Veri Fusion

Roadmap 4.2/4.3 kararlari:
    - Agirliklandirma: IDW (p=2) + yukseklik duzeltmesi (Thiessen degil)
    - Sicaklik lapse rate: sabit -6C/km (standart atmosferik)
    - Yagis: ek yukseklik duzeltmesi yok (bolgesel olcum gerekir, yok -> duz IDW)
    - Referans ET: elimizde zaten olculmus "Buharlasma" (evaporation_m3) var,
      Hargreaves/Thornthwaite hesaplamaya gerek yok

Girdi:
    - hidro_meteoroloji/processed/hydro_met_clean.csv  (Asama A ciktisi, wide format)
    - hidro_meteoroloji/processed/station_locations.csv (DEM-elevation ile)
    - dem_processed/DEM_UTM35N.tif, basin_official.shp (resmi havza siniri --
      istasyon-havza iliskisi / IDW hedef noktasi icin nihai referans olarak
      kabul edildi, bkz. dem_flow_pipeline.py notu; basin_delineated.shp D8'in
      kendi dogrulama surecinde kullandigi ara sonuc, burada KULLANILMIYOR)

Cikti: hidro_meteoroloji/processed/
    - fused_daily_parameters.csv   (date, Q_fused, P_fused, T_fused, E_fused)
    - hydro_met_timeseries_long.csv (date, station_id, parameter, value, qa_flag -- DB'ye hazir)
    - qa_per_station.csv           (istasyon basina veri kalitesi ozeti)

Calistirma: proje kok dizininden
    .venv\\Scripts\\python.exe scripts\\station_fusion.py
"""
import geopandas as gpd
import numpy as np
import pandas as pd
import rasterio

HYDRO_MET_PATH = "hidro_meteoroloji/processed/hydro_met_clean.csv"
STATIONS_PATH = "hidro_meteoroloji/processed/station_locations.csv"
DEM_PATH = "dem_processed/DEM_UTM35N.tif"
BASIN_PATH = "dem_processed/basin_official.shp"
OUT_DIR = "hidro_meteoroloji/processed"

STATIONS = ["tepecik", "aytepe", "kartepe", "cilekli", "kazandere"]
# sutun-adi kisaltmasi -> station_locations.csv'deki station_id.
# NOT: "Kazandere" adinda IKI istasyon var (RG11: yagis/sicaklik, FP2: akis/debi).
# Isimle eslestirme bunlari karistirir (FP2, RG11'in ustune yazar) -- station_id
# kullanmak zorunlu.
STATION_ID_MAP = {
    "tepecik": "RG7", "aytepe": "RG8", "kartepe": "RG9",
    "cilekli": "RG10", "kazandere": "RG11",
}
IDW_POWER = 2
LAPSE_RATE = -6 / 1000  # C / m


def idw_weights(station_xy, target_xy, power=2):
    dist = np.linalg.norm(station_xy - target_xy, axis=1)
    dist = np.where(dist == 0, 1e-6, dist)  # hedef tam istasyonun ustundeyse
    w = 1 / (dist ** power)
    return w / w.sum()


def column_qa_flag(notes, col):
    """hydro_met_clean.csv'deki tek gunluk 'qa_flag' sutunu TUM sutunlarin en
    kotu durumunu tasir (bkz. hydro_met_processing.py); sutun-basina gercek
    ciddiyeti yalnizca 'notes' metnindeki {col}_... etiketinden cikarabiliriz.
    Roadmap 3.2/adim4 olcegiyle tutarli: 0=ok, 1=warning (duzeltildi/tahmini),
    2=reject (istatistiksel aykiri VEYA >30 gun bosluk, deger yok/supheli)."""
    has_gap = notes.str.contains(f"{col}_left_missing_gap", na=False, regex=False)
    has_outlier = notes.str.contains(f"{col}_statistical_outlier", na=False, regex=False)
    has_warning = notes.str.contains(col, na=False, regex=False)
    flag = pd.Series(0, index=notes.index)
    flag[has_warning] = 1
    flag[has_outlier | has_gap] = 2
    return flag


def main():
    print("Veriler yukleniyor...")
    df = pd.read_csv(HYDRO_MET_PATH, parse_dates=["date"])
    stations = pd.read_csv(STATIONS_PATH)
    stations_by_id = {row["station_id"]: row for _, row in stations.iterrows()}

    basin = gpd.read_file(BASIN_PATH)
    target_point = basin.geometry.iloc[0].centroid
    target_xy = np.array([target_point.x, target_point.y])

    with rasterio.open(DEM_PATH) as dem_src:
        r, c = rasterio.transform.rowcol(dem_src.transform, target_point.x, target_point.y)
        target_elev = float(dem_src.read(1)[r, c])
    print(f"  Hedef nokta (havza merkezi): ({target_xy[0]:.1f}, {target_xy[1]:.1f}), "
          f"yukseklik: {target_elev:.1f}m")

    # Istasyon koordinat/yukseklik dizileri (5 met istasyonu, sabit sira)
    station_xy = np.array([
        [stations_by_id[STATION_ID_MAP[s]]["x_utm35n"], stations_by_id[STATION_ID_MAP[s]]["y_utm35n"]]
        for s in STATIONS
    ])
    station_elev = np.array([stations_by_id[STATION_ID_MAP[s]]["elevation_m"] for s in STATIONS])

    weights = idw_weights(station_xy, target_xy, power=IDW_POWER)
    print(f"  IDW agirliklari: {dict(zip(STATIONS, weights.round(3)))}  (toplam={weights.sum():.4f})")

    lapse_correction = LAPSE_RATE * (target_elev - station_elev)
    print(f"  Lapse-rate duzeltmesi (C): {dict(zip(STATIONS, lapse_correction.round(2)))}")

    # --- P_fused: duz IDW (yukseklik duzeltmesi yok -- roadmap 4.3) ---
    precip_cols = [f"precip_{s}_mm" for s in STATIONS]
    P_fused = df[precip_cols].to_numpy() @ weights

    # --- T_fused: once her istasyonu hedef yuksekliğe tasi (lapse), sonra IDW ---
    temp_cols = [f"temp_{s}_c" for s in STATIONS]
    temp_adjusted = df[temp_cols].to_numpy() + lapse_correction  # broadcast
    T_fused = temp_adjusted @ weights

    # --- Q_fused, E_fused: zaten havza-butunlesik olcum, istasyon fuzyonu gerekmez ---
    Q_fused = df["inflow_q_m3s"]
    E_fused = df["evaporation_m3"]

    fused = pd.DataFrame({
        "date": df["date"], "Q_fused_m3s": Q_fused, "P_fused_mm": P_fused,
        "T_fused_c": T_fused, "E_fused_m3": E_fused,
    })
    fused.to_csv(f"{OUT_DIR}/fused_daily_parameters.csv", index=False)
    print(f"\nYazildi: {OUT_DIR}/fused_daily_parameters.csv ({len(fused)} satir)")

    # --- Lapse-rate R^2: gunluk istasyon sicakliklari ile yukseklik arasindaki iliski ---
    sample_days = df.sample(min(200, len(df)), random_state=42)
    r2s = []
    for _, row in sample_days.iterrows():
        temps = row[temp_cols].to_numpy(dtype=float)
        if np.isnan(temps).any():
            continue
        corr = np.corrcoef(station_elev, temps)[0, 1]
        r2s.append(corr ** 2)
    print(f"  Sicaklik-yukseklik iliskisi (ortalama R^2, {len(r2s)} gun ornegi): {np.mean(r2s):.3f}")

    # --- Long format: istasyon parametreleri (DB'ye yuklemeye hazir) ---
    print("\nLong-format tablo olusturuluyor...")
    long_frames = []
    for s in STATIONS:
        sid = STATION_ID_MAP[s]
        for param, col in [("precip_mm", f"precip_{s}_mm"), ("temp_c", f"temp_{s}_c")]:
            if col not in df.columns:
                continue
            sub = df[["date", col]].rename(columns={col: "value"})
            sub["station_id"] = sid
            sub["parameter"] = param
            sub["qa_flag"] = column_qa_flag(df["notes"], col) if "notes" in df.columns else 0
            long_frames.append(sub[["date", "station_id", "parameter", "value", "qa_flag"]])
    long_df = pd.concat(long_frames, ignore_index=True)
    long_df.to_csv(f"{OUT_DIR}/hydro_met_timeseries_long.csv", index=False)
    print(f"  Yazildi: {OUT_DIR}/hydro_met_timeseries_long.csv ({len(long_df)} satir)")

    # --- Istasyon basina QA ozeti ---
    print("\nIstasyon basina QA ozeti hesaplaniyor...")
    qa_rows = []
    for s in STATIONS:
        for param, col in [("precip_mm", f"precip_{s}_mm"), ("temp_c", f"temp_{s}_c"), ("snow_cm", f"snow_{s}_cm")]:
            if col not in df.columns:
                continue
            col_flag = column_qa_flag(df["notes"], col) if "notes" in df.columns else pd.Series(0, index=df.index)
            n_warning = int((col_flag == 1).sum())
            n_reject = int((col_flag == 2).sum())
            flagged = n_warning + n_reject
            qa_rows.append({
                "station_id": STATION_ID_MAP[s], "parameter": param,
                "n_days": len(df), "n_flagged": flagged,
                "n_warning": n_warning, "n_reject": n_reject,
                "flagged_pct": round(100 * flagged / len(df), 2),
            })
    qa_df = pd.DataFrame(qa_rows)
    qa_df.to_csv(f"{OUT_DIR}/qa_per_station.csv", index=False)
    print(f"  Yazildi: {OUT_DIR}/qa_per_station.csv")
    print(qa_df.to_string(index=False))


if __name__ == "__main__":
    main()
