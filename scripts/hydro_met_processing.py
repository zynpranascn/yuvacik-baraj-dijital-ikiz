"""
Faz 1 - Asama A: Hidro-Meteoroloji Veri Temizleme (wide format)

Kaynak: hidro_meteoroloji/raw/Yuvacik_Hidro-Meteoroloji_ Gunluk 21.06.2023.xlsx
        sayfa: "duzenli veri" (2006-01-01 -> 2023-12-31, 5 istasyon)

Bu asama veriyi oldugu gibi (istasyon basina sutun, "wide format") temizler.
Istasyon fuzyonu + long-format donusumu Asama B'de (station_fusion.py) yapilacak.

Calistirma: proje kok dizininden
    .venv\\Scripts\\python.exe scripts\\hydro_met_processing.py
"""
import json
import os

import numpy as np
import pandas as pd
from scipy.stats import zscore

RAW_PATH = "hidro_meteoroloji/raw/Yuvacık_Hidro-Meteoroloji_ Günlük 21.06.2023.xlsx"
SHEET = "düzenli veri"
OUT_DIR = "hidro_meteoroloji/processed"

STATIONS = ["tepecik", "aytepe", "kartepe", "cilekli", "kazandere"]

COLUMN_MAP = {
    "TARİH": "date",
    "Baraj Göl Seviyesi\n (m DSG)": "reservoir_level_m",
    "Faydalı Hacim \n(mM3)": "useful_volume_mcm",
    "Yüzde Doluluk \n(%)": "fill_pct",
    "Arıtma Tesisinden\n Çıkan Su\n(m3/gün)": "outflow_treatment_m3",
    "Baraj Gölüne \nGiren Su\n(m3/gün)": "inflow_m3",
    "Deşarj \n(Radyal Kapak)\n(m3/gün)": "outflow_spillway_m3",
    "Sapanca'dan Temin Edilen Su \n(m3/gün)": "sapanca_transfer_m3",
    "Yağış 1\nTepecik": "precip_tepecik_mm",
    "Yağış 2\nAytepe": "precip_aytepe_mm",
    "Yağış 3\nKartepe": "precip_kartepe_mm",
    "Yağış 4\nÇilekli": "precip_cilekli_mm",
    "Yağış 5\nKazandere": "precip_kazandere_mm",
    "Sıcaklık 1\nTepecik": "temp_tepecik_c",
    "Sıcaklık 2\nAytepe": "temp_aytepe_c",
    "Sıcaklık 3\nKartepe": "temp_kartepe_c",
    "Sıcaklık 4\nÇilekli": "temp_cilekli_c",
    "Sıcaklık 5\nKazandere": "temp_kazandere_c",
    "Kar Yüksekliği 1\nTepecik": "snow_tepecik_cm",
    "Kar Yüksekliği 2\nAytepe": "snow_aytepe_cm",
    "Kar Yüksekliği 3\nKartepe": "snow_kartepe_cm",
    "Kar Yüksekliği 4\nÇilekli": "snow_cilekli_cm",
    # Kaynakta birim belirtilmemis; degerler (0-14,152) mm/gun olarak fiziksel
    # olarak imkansiz (bir golun gunde metrelerce buharlasmasi anlamina gelir).
    # Diger su dengesi sutunlariyla ayni olcekte (m3/gun, hacimsel kayip) --
    # o sekilde ele aliniyor.
    "Buharlaşma": "evaporation_m3",
}

# Fiziksel gecerlilik araliklari: (min, max) - roadmap 3.2 adim 4
RANGE_RULES = {
    "reservoir_level_m": (0, 250),
    "fill_pct": (0, 1.01),
    "inflow_m3": (0, 5_000_000),
    "outflow_treatment_m3": (0, 2_000_000),
    "outflow_spillway_m3": (0, 5_000_000),
    "sapanca_transfer_m3": (0, 3_000_000),
    **{f"precip_{s}_mm": (0, 150) for s in STATIONS},
    **{f"temp_{s}_c": (-20, 50) for s in STATIONS},
    **{f"snow_{s}_cm": (0, 300) for s in STATIONS[:4]},  # kazandere'de kar yuksekligi yok
    "evaporation_m3": (0, 50_000),
}

# Yagis/kar dogas geregi sifir-agirlikli ve saga carpik; IQR/z-score gibi
# simetrik istatistiksel testler yanlislikla normal yagis gunlerini "aykiri"
# olarak isaretler. Bu sutunlarda yalnizca fiziksel aralik kontrolu yapilir.
SKIP_STATISTICAL_OUTLIER = {f"precip_{s}_mm" for s in STATIONS} | {
    f"snow_{s}_cm" for s in STATIONS[:4]
}

# Roadmap 3.3: istasyon olcum sutunlari icin imputasyon uygulanir (su dengesi
# hacim sutunlarina degil -- onlarda eksik/imkansiz deger operasyonel bir
# olayi (transfer yapilmadi vb.) yansitiyor olabilir, tahmin etmek yanlis olur)
IMPUTE_COLS = (
    [f"precip_{s}_mm" for s in STATIONS]
    + [f"temp_{s}_c" for s in STATIONS]
    + [f"snow_{s}_cm" for s in STATIONS[:4]]
    + ["evaporation_m3"]
)


def load_raw():
    df = pd.read_excel(RAW_PATH, sheet_name=SHEET, header=1)
    df = df.rename(columns=COLUMN_MAP)
    df = df.dropna(subset=["date"]).copy()
    df["date"] = pd.to_datetime(df["date"])
    df = df.sort_values("date").drop_duplicates(subset="date").reset_index(drop=True)
    # Debi proxy (m3/s) -- Baraj Golune Giren Su gunluk hacimden
    df["inflow_q_m3s"] = df["inflow_m3"] / 86400
    return df


def flag_out_of_range(df):
    """Fiziksel olarak imkansiz degerleri isaretler VE NaN'a cevirir
    (ör. -170C sicaklik) -- bunlar asagida impute_gaps() ile doldurulur."""
    df["qa_flag"] = 0
    notes = [[] for _ in range(len(df))]
    for col, (lo, hi) in RANGE_RULES.items():
        if col not in df.columns:
            continue
        bad = df[col].notna() & ((df[col] < lo) | (df[col] > hi))
        for idx in df.index[bad]:
            notes[idx].append(f"{col}_out_of_range")
        df.loc[bad, "qa_flag"] = df.loc[bad, "qa_flag"].clip(lower=1)
        if col in IMPUTE_COLS:
            df.loc[bad, col] = np.nan
    return notes


def impute_gaps(df, notes):
    """Roadmap 3.3 imputasyon stratejisi:
    <=3 gun: linear interpolation | 3-30 gun: aylik climatology | >30 gun: birak
    """
    for col in IMPUTE_COLS:
        if col not in df.columns:
            continue
        # Climatology, doldurmadan ONCE mevcut gercek degerlerden hesaplanir
        climatology = df.groupby(df["date"].dt.month)[col].transform("mean")

        missing = df[col].isna()
        if not missing.any():
            continue
        group_id = (missing != missing.shift()).cumsum()
        for gid, group in df[missing].groupby(group_id[missing]):
            idx = group.index
            duration = len(idx)
            if duration <= 3:
                df.loc[idx, col] = df[col].interpolate(method="linear").loc[idx]
                tag = "imputed_linear"
            elif duration <= 30:
                df.loc[idx, col] = climatology.loc[idx]
                tag = "imputed_climatology"
            else:
                tag = f"left_missing_gap_{duration}d"
            for i in idx:
                notes[i].append(f"{col}_{tag}")
    return notes


def flag_outliers(df, notes):
    cols = [c for c in RANGE_RULES if c in df.columns and c not in SKIP_STATISTICAL_OUTLIER]
    for col in cols:
        series = df[col].dropna()
        if len(series) < 30:
            continue
        q1, q3 = series.quantile([0.25, 0.75])
        iqr = q3 - q1
        if iqr == 0:
            continue
        iqr_mask = (df[col] < q1 - 3 * iqr) | (df[col] > q3 + 3 * iqr)
        z = pd.Series(np.nan, index=df.index)
        z.loc[series.index] = zscore(series)
        z_mask = z.abs() > 3
        outlier_mask = (iqr_mask | z_mask) & df[col].notna()
        for idx in df.index[outlier_mask]:
            notes[idx].append(f"{col}_statistical_outlier")
        df.loc[outlier_mask, "qa_flag"] = 2
    return notes


def internal_consistency_checks(df, notes):
    for s in STATIONS:
        col = f"temp_{s}_c"
        if col not in df.columns:
            continue
        dT = df[col].diff().abs()
        anomaly = (dT > 15).fillna(False)  # gunluk >15C sicraması supheli
        for idx in df.index[anomaly]:
            notes[idx].append(f"{col}_sudden_jump")
    return notes


def missing_data_report(df):
    return pd.DataFrame({
        "column": df.columns,
        "missing_count": df.isna().sum().values,
        "missing_pct": (df.isna().mean() * 100).round(2).values,
    })


def time_series_statistics(df):
    stats = {}
    cols = list(RANGE_RULES) + ["inflow_q_m3s"]
    for col in cols:
        if col not in df.columns:
            continue
        s = df[col]
        stats[col] = {
            "mean": float(s.mean()) if s.notna().any() else None,
            "std": float(s.std()) if s.notna().any() else None,
            "min": float(s.min()) if s.notna().any() else None,
            "max": float(s.max()) if s.notna().any() else None,
            "median": float(s.median()) if s.notna().any() else None,
            "monthly_mean": {
                str(m): float(v) for m, v in df.groupby(df["date"].dt.month)[col].mean().items()
            },
        }
    return stats


def main():
    os.makedirs(OUT_DIR, exist_ok=True)

    df = load_raw()
    print(f"Yuklendi: {len(df)} satir, {df['date'].min().date()} -> {df['date'].max().date()}")

    notes = flag_out_of_range(df)
    notes = flag_outliers(df, notes)
    notes = internal_consistency_checks(df, notes)
    notes = impute_gaps(df, notes)
    df["notes"] = [";".join(n) if n else "" for n in notes]

    imputed = df["notes"].str.contains("_imputed_", na=False).sum()
    print(f"Imputasyon yapilan gun sayisi (herhangi bir sutunda): {imputed}")

    miss_report = missing_data_report(df)
    miss_report.to_csv(f"{OUT_DIR}/missing_data_report.csv", index=False)
    print("\nEksik veri raporu (yalnizca >0 olanlar):")
    print(miss_report[miss_report.missing_count > 0].to_string(index=False))

    stats = time_series_statistics(df)
    with open(f"{OUT_DIR}/timeseries_statistics.json", "w", encoding="utf-8") as f:
        json.dump(stats, f, indent=2, ensure_ascii=False)

    qa_summary = df[["date", "qa_flag", "notes"]][df["qa_flag"] > 0]
    qa_summary.to_csv(f"{OUT_DIR}/qa_flags.csv", index=False)
    print(f"\nQA flag > 0 olan gun sayisi: {len(qa_summary)} / {len(df)} (%{100*len(qa_summary)/len(df):.1f})")

    df.to_csv(f"{OUT_DIR}/hydro_met_clean.csv", index=False)
    df.to_parquet(f"{OUT_DIR}/hydro_met_clean.parquet", index=False)
    print(f"\nCikti yazildi: {OUT_DIR}/hydro_met_clean.csv (+ .parquet)")


if __name__ == "__main__":
    main()
