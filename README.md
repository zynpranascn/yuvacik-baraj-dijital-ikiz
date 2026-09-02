# Yuvacık Barajı Dijital İkizi — Faz 1: Veri Hazırlığı

Yuvacık Barajı için DEM, 20 yıllık hidro-meteoroloji verisi, istasyon konumları
ve dolin (karst obruk) verisini birleştirip; havzayı gerçek debi ölçüm
noktalarına (FP1, FP2, FP3) + barajın (RG6) yerel havzasına göre alt-havzalara
ayıran, gelecekteki kuraklık/dolum tahmini modelinin veri temelini kuran çalışma.

## Kurulum

```
git clone <bu-repo-URL'si>
cd GEOTİFF_EPSG_32635
python -m venv .venv
.venv\Scripts\activate          # Windows
pip install -r requirements.txt
python verify_env.py            # paketleri + DEM okumayı doğrular
```

Python 3.11 ile test edildi. `dem.tif`, `dem_geo.tif`, `dere.tif`,
`hillshade.tif`/`hillshade_geo.tif`, `hava_dis_sinir.tif`, `istasyon.tif`/
`istasyon_geo.tif`, `rezervuar.tif` (~600MB, hiçbir scriptte okunmuyor, salt
görsel/tarihsel) repoya dahil **değildir** (`.gitignore`). `hava_dis_sinir_geo.tif`,
`rezervuar_geo.tif`, `dere_geo.tif` ise **dahildir** — `extract_official_basin.py`
(pipeline adım 3) bunları gerçekten girdi olarak okuyor. Pipeline'ı adım 1'den
(kontur → DEM) çalıştırmak isterseniz `../../dem/dem/dem.gdb` (Esri File
Geodatabase) ve karşılaştırma için `dem_geo.tif` ayrıca gerekir, bu repo dışında
tutulur; aksi halde `dem_processed/DEM_UTM35N.tif` zaten üretilmiş halde repoda
mevcuttur ve adım 2'den devam edilebilir (bkz. "Pipeline" bölümü).

## İçindekiler
1. [Hızlı bakış — hangi dosya ne işe yarıyor](#hızlı-bakış)
2. [Dizin yapısı](#dizin-yapısı)
3. [Veri işleme hattı (pipeline) — çalıştırma sırası](#pipeline)
4. [Alt-havza modeli (FP1–FP3 + RG6)](#alt-havza-modeli)
5. [Bilinen sınırlamalar / belgelenmiş override'lar](#bilinen-sınırlamalar)

---

## Hızlı bakış

En sık ihtiyaç duyulacak dosyalar:

| İhtiyaç | Dosya |
|---|---|
| Havza + alt-havza + dere + dolin + istasyon + rezervuar — tek harita | `dem_processed/master_overview_map.png` |
| Alt-havza sınırları (GIS) | `dem_processed/subbasins_by_station.shp` |
| Alt-havza morfometrisi (alan, eğim, drenaj yoğunluğu, yağış) | `dem_processed/subbasin_morphometry.csv` |
| Dere ağı, alt-havzaya atanmış | `dem_processed/streams_by_station.shp` |
| Dere kollarının uzunluğu (sıralı) | `dem_processed/stream_branch_lengths_sorted.csv` |
| Dolin noktaları + hangi dereye/alt-havzaya bağlı | `hidro_meteoroloji/processed/doline_points_enriched.csv` |
| Temizlenmiş 18 yıllık hidro-meteoroloji verisi | `hidro_meteoroloji/processed/hydro_met_clean.csv` |
| İstasyon konumları (FP1–FP3, RG6–RG12) | `hidro_meteoroloji/processed/station_locations.csv` |
| Her alt-havza için ayrı klasör (DEM+dere+dolin+istasyon+iklim) | `dem_processed/subbasins/<FP1\|FP2\|FP3\|RG6>/` |
| Akış yönlendirme (FP1/FP2/FP3 → RG6/baraj) | `dem_processed/subbasin_flow_routing.csv` |

---

## Dizin yapısı

```
GEOTİFF_EPSG_32635/                    <- proje kök dizini
│
├── README.md                          <- bu dosya
├── requirements.txt                   <- Python bağımlılıkları
├── verify_env.py                      <- ortam kurulum doğrulama scripti
├── .venv/                             <- Python sanal ortamı
│
├── dem.tif, dere.tif, hillshade.tif,          ┐
│   hava_dis_sinir.tif, istasyon.tif,          │  HAM RENDER dosyaları (kaynak
│   rezervuar.tif  (+ *_geo.tif/.tfw            │  GIS export'undan) — CRS'siz
│   varyantları)                                │  veya piksel-tabanlı; DEM_UTM35N.tif
│   dömüşüm.prj / dömüşüm.txt                  ┘  ÜRETİLDİKTEN sonra sadece
│                                                  görsel referans/karşılaştırma
│                                                  için kullanıldı, artık aktif
│                                                  pipeline'da GİRDİ DEĞİLLER.
│
├── scripts/                           <- TÜM kalıcı Python pipeline kodu
│   (bkz. "Pipeline" bölümü — çalıştırma sırasıyla)
│
├── hidro_meteoroloji/
│   ├── raw/                           <- HAM kaynak veri (dokunulmamış)
│   │   ├── İSAŞ_Barajlar Günlük Veri_TEMMUZ_2026.xlsx   (istasyon koordinat tablosu)
│   │   ├── Yuvacık_Hidro-Meteoroloji_..._2023.xlsx      (18 yıllık günlük ölçüm)
│   │   ├── DemWhatsApp Image....jpeg                    (referans havza görseli)
│   │   ├── istasyon konumları/                          (Thiessen poligon görselleri)
│   │   └── dolin/dolin_geo.tif                          (dolin/karst render'ı)
│   │
│   └── processed/                     <- Temizlenmiş / türetilmiş çıktılar
│       ├── hydro_met_clean.csv(+.parquet)   (Aşama A çıktısı — 6574 gün, QA'lı)
│       ├── missing_data_report.csv, qa_flags.csv, timeseries_statistics.json
│       ├── station_locations.csv(+.shp)     (FP1–FP3, RG6–RG12 gerçek koordinat — FP4 yok, bkz. Bilinen sınırlamalar)
│       ├── fused_daily_parameters.csv       (IDW+lapse-rate ile havza-geneli Q/P/T/E)
│       ├── hydro_met_timeseries_long.csv    (DB'ye yüklemeye hazır uzun format)
│       ├── qa_per_station.csv
│       ├── doline_points.csv(+.shp)         (ham dolin çıkarımı, 53 nokta)
│       └── doline_points_enriched.csv(+.shp) (+ alt-havza, dere, karst override bilgisi)
│
├── dem_processed/                     <- TÜM DEM/GIS türetilmiş çıktılar
│   ├── DEM_UTM35N.tif                 <- GERÇEK DEM (konturdan üretildi, doğrulandı)
│   ├── flow_direction_D8.tif, flow_accumulation.tif
│   ├── basin_official.shp             <- REFERANS havza sınırı (260.7 km², 4 bağımsız
│   │                                      yöntemle doğrulandı — bkz. validation_report.json)
│   ├── basin_delineated*.shp          <- D8 ara-sonuçları (dogal/birlesik havza)
│   ├── reservoir_official.shp         <- baraj gölü poligonu (1.54 km²)
│   ├── streams_official.shp           <- resmi kaynak dere ağı (render'dan)
│   ├── stream_network.geojson         <- D8-türetilmiş tam dere ağı (1427 segment)
│   │
│   ├── subbasins_by_station.shp       <- ★ GÜNCEL alt-havza sınırları (FP1–FP3, RG6)
│   ├── streams_by_station.shp         <- ★ dere ağı, alt-havzaya atanmış (RG12 bağlantısı dahil)
│   ├── streams_rg12_connector.shp     <- RG12↔FP1 bağlantısının kaynağı (artık streams_by_station.shp'ye gömülü)
│   ├── subbasin_flow_routing.csv      <- ★ FP1/FP2/FP3 → RG6 (baraj) akış yönlendirme tablosu
│   ├── subbasin_morphometry.csv       <- ★ alan/eğim/drenaj yoğunluğu/yağış tablosu
│   ├── stream_branch_lengths_sorted.csv <- ★ her dere kolunun uzunluğu (sıralı)
│   ├── reference_subbasins_from_map.shp <- referans Excel haritasından dijitize 4 bölge
│   │                                        (çapraz doğrulama için, bkz. aşağı)
│   │
│   ├── master_overview_map.png        <- ★ TEK KAPSAMLI ÖZET HARİTA
│   ├── doline_stream_connectivity_map.png <- dolin↔dere bağlantı haritası
│   ├── subbasins_reference_comparison.png <- D8 vs referans harita karşılaştırması
│   ├── subbasins_by_station_map.png, dolines_basin_overlay.png,
│   │   basin_validation_map.png, rg9_flow_trace.png  <- ara-aşama doğrulama görselleri
│   │
│   ├── dem_comparison_report.json     <- kontur-DEM vs renk-DEM doğruluk karşılaştırması
│   ├── validation_report.json         <- D8 pipeline'ının kendi doğrulama raporu
│   │                                      (NOT: buradaki havza alanı rakamları ESKİ/ara-
│   │                                      sonuç, basin_official.shp'yi esas alın)
│   │
│   └── subbasins/                     <- ★ HER ALT-HAVZA İÇİN AYRI, KENDİ KENDİNE
│       ├── FP1/  (Kirazdere, 83.9 km², RG8+RG12 dahil)  YETERLİ KLASÖR: boundary.shp,
│       ├── FP2/  (Kazandere, 23.2 km²)                   dem.tif, flow_direction_D8.tif,
│       ├── FP3/  (Serindere, 122.2 km²)                  flow_accumulation.tif, streams.shp,
│       └── RG6/  (Baraj yerel havzası, 31.4 km²)         dolines.csv/.shp, stations.csv,
│                                                          hydro_met.csv (varsa), metadata.json
│
└── (üst dizinde, proje dışı) ../../dem/dem/dem.gdb, ../../dem/tin/
                                        <- Esri File Geodatabase — DEM_UTM35N.tif'in
                                           ÜRETİLDİĞİ orijinal kontur (izhops) verisi.
                                           Sadece build_dem_from_contours.py tarafından
                                           bir kereliğine okunur, artık pipeline'da
                                           tekrar gerekmez.
```

`★` = en sık kullanılacak / en güncel dosyalar.

---

## Pipeline

Script'ler **çalıştırma sırasına göre** (`scripts/`, hepsi proje kök dizininden
`.venv\Scripts\python.exe scripts\<isim>.py` ile çalıştırılır):

| # | Script | Girdi | Çıktı |
|---|---|---|---|
| 1 | `build_dem_from_contours.py` | `../../dem/dem/dem.gdb` (kontur) | `DEM_UTM35N.tif` |
| 2 | `dem_flow_pipeline.py` | `DEM_UTM35N.tif` | `flow_direction_D8.tif`, `flow_accumulation.tif`, `basin_delineated*.shp` |
| 3 | `extract_official_basin.py` | `hava_dis_sinir_geo.tif`, `rezervuar_geo.tif`, `dere_geo.tif` | `basin_official.shp`, `reservoir_official.shp`, `streams_official.shp` |
| 4 | `georeference_reference_map.py` | WhatsApp referans görseli | `basin_delineated_reference_map.shp` (çapraz doğrulama) |
| 5 | `hydro_met_processing.py` | Ham Excel (hidro-met) | `hydro_met_clean.csv` |
| 6 | `build_station_locations.py` | İSAŞ Excel (koordinat tablosu) | `station_locations.csv` (RG6–RG12, FP1–FP3) |
| 7 | `station_fusion.py` | `hydro_met_clean.csv` + `station_locations.csv` | `fused_daily_parameters.csv`, `hydro_met_timeseries_long.csv` |
| 8 | `extract_doline_points.py` | `dolin_geo.tif` | `doline_points.csv/.shp` |
| 9 | `assign_streams_dolines_to_stations.py` | yukarıdakilerin hepsi | **`subbasins_by_station.shp`, `streams_by_station.shp`, `doline_points_enriched.csv`, `subbasin_flow_routing.csv`** |
| 10 | `split_data_by_subbasin.py` | #9 çıktıları | `dem_processed/subbasins/<sid>/` klasörleri |
| 11 | `subbasin_morphometry.py` | #10 çıktıları | `subbasin_morphometry.csv` |

Sadece görselleştirme/doğrulama amaçlı (pipeline'ın parçası değil, istendiğinde
tekrar çalıştırılır): `plot_basin_check.py`, `plot_dolines_on_basin.py`,
`plot_rg9_zoom.py`, `plot_master_overview_map.py` (★ `master_overview_map.png`'yi
üretir — adım 9/10 sonrası bir değişiklik olduğunda tekrar çalıştırın).

**Artık kullanılmıyor (superseded):** `add_fp4_station.py` — FP4'ü referans
haritadan digitize edip 4. bir debi istasyonu olarak ekliyordu. Kullanıcı
kararıyla FP4 kaldırıldı (bkz. "Bilinen sınırlamalar") — script dosyası
tarihsel referans için duruyor ama artık çalıştırılmıyor, `station_locations.csv`
FP4 satırını içermiyor.

**En sık yapılacak değişiklik:** Alt-havza mantığında bir düzeltme gerekirse,
sadece **adım 9, 10, 11**'i bu sırayla yeniden çalıştırmak yeterli — DEM/havza
sınırı/hidro-met temizliği (1-8) değişmediği sürece tekrar çalıştırılmaz.

---

## Alt-havza modeli

Havza, **3 gerçek debi istasyonuna** (FP1, FP2, FP3) + barajın (RG6) kendi
yerel havzasına göre bölünmüştür. Her alt-havza resmi referans sınırına
(`basin_official.shp`) kırpılmıştır, bu yüzden toplam alan resmi sınırla
**birebir** eşleşir:

| Kod | Dere/bölge adı | Alan | Çevre | Tip |
|---|---|---|---|---|
| **FP1** | Kirazdere | 83.95 km² | 75.00 km | debi istasyonu (RG8 + RG12 dahil, aşağı bkz.) |
| **FP2** | Kazandere | 23.19 km² | 31.06 km | debi istasyonu |
| **FP3** | Serindere | 122.19 km² | 80.94 km | debi istasyonu |
| **RG6** | Baraj Yerel Havzası | 31.36 km² | 37.44 km | baraj_yerel_havza (istasyon değil — kendi yerel akışı + Kartepe karst dolinleri dahil) |

Not: RG8 (Aytepe) ve RG12 (Haci Osman) artık ayrı alt-havza değil — ikisi de
FP1'e union edildi (bkz. Bilinen sınırlamalar). İstasyonların kendileri hâlâ
`subbasins/FP1/stations.csv`'de ayrıca listelenir, sadece kendi poligonları yok.

**Akış yönlendirmesi** (`dem_processed/subbasin_flow_routing.csv`): FP1, FP2 ve
FP3 **birbirinden bağımsız 3 ayrı koldur** — her biri kendi alt-havzasından
akar, yukarıda birbirleriyle BİRLEŞMEZLER. Tek ortak nokta **barajın kendisi**:
- **FP2** (Kazandere) → kendi (yeşil) alt-havzasından akar, doğrudan baraja ulaşır
- **FP3** (Serindere) → kendi (mavi) alt-havzasından akar, RG6 yerel havzasından
  (sarı bölge) **geçerek** baraja ulaşır
- **FP1** (Kirazdere) → ayrı/harici bir sistemle (muhtemelen tünel transferi)
  baraja ulaşır

Üçü de sonunda barajda "buluşur" ama bu, tek bir ortak yukarı-akış kanalını
paylaştıkları anlamına gelmez — su dengesi/modelleme açısından üç bağımsız
girdi olarak ele alınmalıdır.

**Toplam: 260.69 km² — resmi referansla (`basin_official.shp`) tam eşleşiyor**
(fark: 0.0000000 km²). Alt-havza morfometrisi (eğim, rölyef, drenaj yoğunluğu,
yağış) için bkz. `dem_processed/subbasin_morphometry.csv`.

## Literatür referansları (Ayfer hoca makalesi)

**Kaynak:** Özdemir, A. (2021). İklim Değişikliğinin Havza Ölçeğinde Akım ve Sediman
Miktarına Etkilerinin Değerlendirilmesi: Yuvacık Baraj Gölü Havzası. *Jeoloji
Mühendisliği Dergisi*, 45(1), 129-153. DOI: 10.24232/jmd.941528

`hidro_meteoroloji/raw/ayfer_hoca_referans/` altında makalenin tamamı
(`ayfer_hoca_makale_kaynak.pdf`, 25 sayfa) ve içindeki 5 şeklin tamamı 400 DPI'da,
tam başlık/lejant/alıntı metniyle kesilmiş olarak duruyor:

| Dosya | İçerik (makaledeki adı) | Durum |
|---|---|---|
| `sekil1_calisma_alani_ve_dereler.png` | Şekil 1 — Çalışma alanı + Türkiye lokasyon haritası; alt havza haritası (Kirazdere/Serindere/Kazandere/Ara, isimli dereler, RG1–RG12, FP1–FP4, iklim grid) | **Aktif** — `reference_subbasins_from_map.shp` dijitizasyonunun kaynağı; FP1=Kirazdere, FP2=Kazandere, FP3=Serindere adlandırmasını doğruluyor |
| `sekil4_composite_a_e_numarali_alt_havza.png` | Şekil 4 — a: SAM/DEM, b: **numaralı** alt havza haritası (63 alt havza, mağara/kaptaj/keson kuyu + kirlilik noktaları), c: eğim, d: toprak, e: arazi kullanımı | **Aktif (b paneli)** — `hidro_meteoroloji/processed/doline_points_enriched.csv`'deki `karst_magara_override` kararının literatür dayanağı; **Aktif (a paneli, destekleyici)** — DEM yükseklik aralığı (140–1530 m) `DEM_UTM35N.tif` ile (118.8–1525 m) örtüşüyor; c/d/e panelleri **pasif**, ileride eğim/toprak/arazi örtüsü fazında kullanılacak |
| `sekil3_gozlenen_aylik_akim_degerleri.png` | Şekil 3 — FP1 Kirazdere, FP2 Kazandere, FP3 Serindere için 2006-2015 gözlenen aylık ortalama akım (m³/s) | **Aktif** — tam da sizin izlediğiniz 3 istasyon (FP1/FP2/FP3); `station_locations.csv`/`hydro_met_clean.csv`'deki günlük akım verisiyle doğrudan çapraz doğrulanabilir |
| `sekil5_gozlenen_vs_modellenen_akim.png` | Şekil 5 — FP1/FP2/FP3 için gözlenen vs SWAT-modellenen akım (R², 1:1 çizgisi) | Pasif (referans) — SWAT modellemesi sizin Faz 1 kapsamınızda değil, ama gelecekte bir hidrolojik model kurarsanız kalibrasyon karşılaştırması için örnek teşkil eder |
| `sekil2_aylik_yagis_grafigi.png` | Şekil 2 — RG1–RG12 istasyonlarında ölçülen ortalama aylık yağış toplamı (mm) | **Pasif** — `hydro_met_clean.csv`/`fused_daily_parameters.csv`'de günlük çözünürlükte veri zaten var; detaylı veri gelince kendi üretilecek eşdeğer grafikle karşılaştırılacak |

### Alan karşılaştırması (makale Çizelge 1 vs. sizin `subbasin_morphometry.csv`)

Makalenin Çizelge 1'i, README'deki "Bilinen sınırlamalar" bölümünde kaynağı belirsiz
olarak not edilen **~258 km²** rakamının tam kaynağıdır:

| Alt havza | Makale (km²) | Sizin veri (km²) | Fark |
|---|---|---|---|
| Kirazdere (FP1) | 79.54 | 83.95 | +4.41 (+5.5%) |
| Kazandere (FP2) | 23.10 | 23.19 | +0.09 (+0.4%) |
| Serindere (FP3) | 120.53 | 122.19 | +1.66 (+1.4%) |
| Ara Havza / RG6 | 34.69 | 31.36 | −3.33 (−9.6%) |
| **Toplam** | **257.86** | **260.69** | +2.83 (+1.1%) |

Toplamdaki ~%1'lik fark makalenin kendi doğrulama aralığıyla (~%0.1) uyumlu
büyüklükte değil ama makul bir delineasyon-yöntemi farkı; Ara Havza/RG6'daki
%9.6'lık fark en dikkat çekici olanı — muhtemelen RG8/RG12'nin FP1'e union
edilmesi (bkz. aşağı) ile makalenin "Ara Havza" tanımının tam örtüşmemesinden
kaynaklanıyor, ileride incelenebilir.

---

## Bilinen sınırlamalar

Şeffaflık için: bazı alt-havza sınırları **saf D8 hesabından değil**, belgelenmiş
override'lardan geliyor (hepsi `scripts/assign_streams_dolines_to_stations.py`
içinde yorum satırlarıyla gerekçelendirilmiş):

- **RG12 → FP1**: DEM'in bu köşesinde kontur kapsamı dışı düz-plato artefaktı
  yüzünden D8 akışı yanlış yöne gidiyor. Referans harita + saha bilgisi temel
  alınarak FP1'e union edildi. Dere bağlantısı için `streams_rg12_connector.shp`
  (DEM en-kısa-yol tahmini, gerçek D8 değil) ayrıca üretildi — bu artık
  **2 segment** içeriyor: (1) RG12 istasyon noktası → FP1 ana ağı (orijinal
  bağlantı), (2) RG12 istasyon noktası → en yakın D8-türetilmiş dere segmenti
  (Hekim Dere devamı), aralarındaki **344 m'lik görsel kopukluğu** kapatan ek
  segment (2026-08-21'de eklendi, master haritada fark edilen bir boşluktu).
  İkisi de `source=manual_dem_leastcost` etiketiyle `streams_by_station.shp`'ye
  gömülü.
- **Kartepe/RG9 dolinleri → RG6**: 10 dolin coğrafi olarak resmi havza sınırının
  DIŞINDA (`inside_basin_official=False`), ama karst/mağara bağlantısıyla suyun
  baraja ulaştığı bilgisine dayanarak RG6'ya (baraj yerel havzası) atandı.
  Haritada **üçgen sembolle** (normal dolinlerden ayrı) işaretli, veri
  tablosunda `assign_method=karst_magara_override` ile etiketli (toplam 53
  dolinin 43'ü `d8_surface`, 10'u `karst_magara_override`).
- **FP4 kaldırıldı**: Önceki bir revizyonda "FP4" adında 4. bir debi istasyonu
  vardı (referans haritadan digitize edilmişti). Kullanıcı kararıyla kaldırıldı
  — o coğrafi bölge hâlâ duruyor ama artık bağımsız bir istasyon değil,
  **RG6'nın (barajın) yerel havzası**: FP2 ve FP3'ün suyu buradan geçerek baraja
  ulaşıyor (D8 ile doğrulandı — RG6'nın ham catchment'i FP2+FP3'ü %100
  içeriyordu). `scripts/add_fp4_station.py` artık çalıştırılmıyor (tarihsel
  referans için duruyor). **Not:** Özdemir (2021) Şekil 1'de FP4 gerçekten
  görünüyor ama bağımsız bir akım istasyonu değil, "ISU Ara Depo" (baraja
  gelen suyun %2'sini/3.5 milyon m³'ünü sağlayan bir ara depolama noktası) —
  yani makale de FP4'ü Kirazdere/Kazandere/Serindere ile aynı kategoride
  görmüyor, bu da FP4'ü ayrı bir istasyon olarak tutmama kararını destekliyor.
- **Resmi sınıra kırpma**: D8-türetilmiş alt-havza poligonları başlangıçta
  `basin_official.shp` ile birebir örtüşmüyordu (63 küçük parçada, toplam ~5 km²
  taşma vardı — haritada renkli dolgunun kesikli referans çizgisinin dışına
  taştığı yerler olarak görünüyordu). Her alt-havza artık `basin_official.shp`
  ile kesiştirilerek (intersection) kırpılıyor, bu yüzden toplam alan resmi
  referansla tam eşleşiyor (260.69 km², fark yok).
- **RG8 → FP1**: RG8'in kendi D8 catchment'i (0.11 km²) flow-path izlemeyle
  FP1'in ana dere ağına doğru bağlandığı doğrulandı (RG12'deki gibi bir DEM
  artefaktı YOK) — ama ayrı bir kategori olarak tutulunca haritada dere ağına
  hiç bağlanmamış izole bir "ada" gibi görünüyordu (kullanıcı referans harita
  ile tutarsızlık olarak fark etti). RG12 gibi FP1'e union edildi; istasyonun
  kendisi hâlâ FP1'in "içindeki istasyonlar" listesinde ayrıca görünür.
- **`validation_report.json`** güncel değil — eski (RG6-merkezli) havza
  rakamlarını taşıyor, `basin_official.shp`/`subbasins_by_station.shp`'yi esas alın.
- **Toplam alan ~258 km² değil, 260.69 km² — KAYNAK BULUNDU**: Kullanıcının
  literatürden hatırladığı ~258 km² rakamının kaynağı artık netleşti: Özdemir
  (2021) Çizelge 1'de "Yuvacık Barajı Toplam Drenaj Alanı: **257.86 km²**"
  olarak veriliyor (bkz. "Literatür referansları" bölümü, alan karşılaştırma
  tablosu). Mevcut 260.69 km² dört BAĞIMSIZ yöntemle çapraz doğrulandı
  (CRS/dönüşüm kontrolü, piksel-renk analizi, bağlı-bileşen analizi, D8 fiziksel
  hesap — hepsi 260.6–260.8 km² aralığında, ~%0.1 sapmayla birbirini
  destekliyor), bu yüzden DEĞİŞTİRİLMEDİ — makaleyle ~%1.1'lik fark muhtemelen
  farklı delineasyon yöntemi/DEM çözünürlüğü/outlet tanımından kaynaklanıyor
  (alt-havza bazında en büyük sapma Ara Havza/RG6'da, ~%9.6).
