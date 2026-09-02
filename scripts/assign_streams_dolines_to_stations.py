"""
Faz 1 - Alt-havzalari GERCEK DEBI ISTASYONU (FP1/FP2/FP3) catchment'lerine +
RG6 (baraj) yerel havzasina gore tanimlar, sonra TUM dere agi segmentlerini
ve dolin noktalarini bagli olduklari havzaya atar.

REVIZYON NOTU (kullanici karari): Onceki versiyonda "FP4" adinda 4. bir
debi istasyonu vardi (referans haritadan digitize edilmisti, bkz. eski
scripts/add_fp4_station.py). Kullanici bu noktayi KALDIRMAYA karar verdi --
FP4 hicbir resmi tabloda yok (sadece haritada, digitize edilmis), gercek bir
olcum noktasi olarak guvenilir degil. YERINE: o cografi bolge (eskiden
"FP4" diye adlandirilan sari/turuncu alan) hala 4. bir alan olarak duruyor,
ama artik kendi basina bir "istasyon" degil -- RG6'nin (barajin) YEREL
HAVZASI olarak taniminiyor: FP2 VE FP3'un suyu bu bolgeden GECEREK baraja
dokuluyor (D8 ile dogrulandi, asagida), + bu bolgenin kendi yerel akisi.

    3 gercek debi istasyonu: FP1 (Kirazdere), FP2 (Kazandere), FP3 (Serindere).
    RG6 (Dolusavak/baraj) -- rezervuarin kendisi, nihai TOPLAMA noktasi.

Hiyerarsi / akis yonlendirme (D8 ile dogrulandi, bkz. arastirma asamasi):
    FP1 -> RG6:  66.91 km2 -- RG6'nin dogal havzasina (natural, 174.66 km2)
         D8 ile BAGLI DEGIL (ayri/harici sistem, muhtemelen tunel transferi)
         -- referans haritada da FP1'in bolgesi (mor/pembe) digerlerinden
         tamamen ayri/dokunmuyor, bu bagimsizligi teyit ediyor. Yine de
         suyu baraja ulasiyor (bkz. FLOW_ROUTING).
    FP2 -> RG6 (yerel havza uzerinden):  23.41 km2 -- natural'in icinde,
         FP3 ile kesismiyor (kardes kollar).
    FP3 -> RG6 (yerel havza uzerinden):  120.31 km2 -- natural'in icinde,
         FP2 ile kesismiyor.
    RG6 yerel havzasi (eskiden "FP4"): natural - FP2 - FP3 = 30.94 km2 --
         (RG6'nin D8-snap noktasindan hesaplanan ham deger FP2 VE FP3'u
         TAMAMEN icine aliyordu -- %100 ic ice -- bu, FP2/FP3'un suyunun
         gercekten bu bolgeden GECEREK baraja ulastiginin D8 kaniti).
         FP2+FP3 cikarildiktan sonra kalan bu 30.94 km2, bu bolgenin
         KENDINE OZGU (yerel, baska hicbir istasyondan gecmeyen) katkisi.

FLOW_ROUTING (asagida, ayrica dosyaya da yazilir --
subbasin_flow_routing.csv): FP1, FP2 ve FP3'un UCU de nihayetinde RG6'ya
(baraj golune) dokulur. FP2 ve FP3 bunu RG6 yerel havzasindan GECEREK yapar
(D8 ile dogrulanmis fiziksel gecis); FP1 ayri/harici bir sistem uzerinden
(muhtemelen tunel) dogrudan baraja ulasir.

RG8 (Aytepe, 0.11 km2): D8 ile FP1'in ana agina duzgun bagli (flow-path
izlemeyle dogrulandi, artefakt yok) ama kendi basina ayri bir kategori olarak
tutulunca haritada dere agina baglanmamis izole bir "ada" gibi gorunuyordu --
kullanici bunu tutarsizlik olarak isaretledi. Bu yuzden RG12 gibi FP1'e union
edildi (bkz. RG8_MERGED_INTO). Istasyonun kendisi hala FP1'in "icindeki
istasyonlar" listesinde ayrica gorunur, sadece kendi poligonu yok artik.

RG12 (Haci Osman, 0.74 km2) ise D8 ile FP1'e (veya baska hicbir seye)
BAGLANMIYOR -- bkz. RG12_MERGED_INTO yorumu (asagida): DEM'in bu kosesinde
kontur kapsami disi bir "duz plato" artefakti akisi yanlis yone gonderiyor.
Kullanicinin referans haritasi + saha bilgisiyle dogrulanan gercek baglanti
(FP1) icin RG12'nin catchment'i BURADA FP1'e union edilir -- yani RG12
ARTIK ayri bir alt-havza olarak CIKMAZ, FP1'in bir parcasi olarak gorunur
(istasyon noktasi haritada hala ayrica etiketlenir).

Cikti: dem_processed/
    - subbasins_by_station.shp    (station_id, name, kind, area_km2)
    - streams_by_station.shp      (her dere segmentine station_id atanmis)
    - subbasin_flow_routing.csv   (FP1/FP2/FP3 -> RG6 akis yonlendirme tablosu)
    - subbasins_by_station_map.png
  hidro_meteoroloji/processed/
    - doline_points_enriched.csv (station_id / kind sutunlari eklendi, uzerine yazilir)

Calistirma: proje kok dizininden
    .venv\\Scripts\\python.exe scripts\\assign_streams_dolines_to_stations.py
"""
import json

import numpy as np

if not hasattr(np, "in1d"):
    np.in1d = np.isin

import geopandas as gpd
import matplotlib.pyplot as plt
import pandas as pd
import rasterio
from pysheds.grid import Grid
from rasterio.features import shapes
from shapely.geometry import shape

DIRMAP = (64, 128, 1, 2, 4, 8, 16, 32)
OUT_DIR = "dem_processed"
PROCESSED_DIR = "hidro_meteoroloji/processed"

# Debi olculen (Akis/Debi) istasyonlar -- dere/dolin atamasinin ANA kriteri.
# NOT: FP4 kaldirildi (bkz. modul docstring, "REVIZYON NOTU") -- artik sadece
# 3 gercek debi istasyonu var, 4. alan (RG6) asagida ayrica hesaplaniyor.
DISCHARGE_STATIONS = ["FP1", "FP2", "FP3"]
MASTER_OUTLET = "RG6"  # baraj -- FP2/FP3'un GECTIGI + kendi yerel havzasi

# Bilgi amacli, debi olculmeyen istasyonlar -- artik BOS: RG8 de RG12 gibi
# FP1'e union edildi (bkz. RG8_MERGED_INTO, asagida) cunku RG8'in kendi D8
# catchment'i (0.11 km2) FP1'in ortasinda dere agina baglanmayan IZOLE bir
# "ada" gibi goruntulendiginden kullanici tarafindan tutarsizlik olarak
# isaretlendi (referans haritada RG8 boyle bir bosluk/ada olarak gorunmuyor).
INFO_ONLY_STATIONS = []

# RG12'nin kendi D8 catchment'i (0.74 km2) FP1'e (ana govde) BAGLANMIYOR --
# flow-path izlemeyle dogrulandi: RG12'den baslayan akis, DEM'in bu kosesinde
# ~1.2km boyunca TAM DUZ (sabit 850.0m, tam kontur degeri) bir "plato"
# uzerinden batiya/harita disina gidiyor. Bu, gercek arazi degil -- kontur
# ucgenlemesinin convex hull'u disinda kalan bir cep, nearest-neighbor
# doldurmayla (build_dem_from_contours.py) TEK bir komsu kontur noktasinin
# degerine sabitlenmis (bkz. arastirma asamasi). Yani DEM bu bolgede akis
# yonunu guvenilir sekilde cozemiyor.
# Kullanicinin referans haritasi (xl/media/image1.png) VE saha bilgisi
# ("RG12 zaten dogrudan FP1'e gidiyor") RG12'nin FP1 sistemine bagli
# oldugunu acikca gosteriyor -- bu yuzden RG12'nin catchment'i BURADA,
# belgelenmis bir override ile FP1'e union ediliyor (asagida, main() icinde).
RG12_MERGED_INTO = "FP1"

# RG8'in kendi D8 catchment'i (0.11 km2) FLOW-PATH IZLEMEYLE dogrulandi ki
# gercekten FP1'in ana agina baglaniyor (281 adimda accum=171441'e ulasiyor,
# RG12'deki gibi bir "duz plato" artefakti YOK). Yani bu bir D8 hatasi degil
# -- sadece cift-sayimi onlemek icin FP1'den cikarilip ayri "bilgi amacli"
# bir kategori olarak tutulmustu. Ama bu, haritada RG8'i dere agina hic
# baglanmamis, kirmizi bolgenin ortasinda IZOLE bir "ada" gibi gosteriyordu
# -- kullanici bunu referans haritayla tutarsizlik olarak isaretledi. Cozum:
# RG8'i de RG12 gibi FP1'e union et (D8 zaten dogru, sadece kategori/gorsel
# ayrimini kaldiriyoruz).
RG8_MERGED_INTO = "FP1"


def snap_outlet(acc_arr, transform, approx_x, approx_y, margin=5, win=15):
    approx_row, approx_col = rasterio.transform.rowcol(transform, approx_x, approx_y)
    r0, r1 = max(margin, approx_row - win), min(acc_arr.shape[0] - margin, approx_row + win)
    c0, c1 = max(margin, approx_col - win), min(acc_arr.shape[1] - margin, approx_col + win)
    window = acc_arr[r0:r1, c0:c1]
    local_r, local_c = np.unravel_index(np.argmax(window), window.shape)
    row, col = r0 + local_r, c0 + local_c
    x, y = transform * (col, row)
    return x, y


def delineate(grid, fdir, x, y):
    catchment = grid.catchment(x=x, y=y, fdir=fdir, dirmap=DIRMAP, xytype="coordinate")
    arr = np.array(catchment).astype("uint8")
    polys = [
        shape(g) for g, v in shapes(arr, mask=arr.astype(bool), transform=grid.affine)
        if v == 1
    ]
    polys.sort(key=lambda p: p.area, reverse=True)
    return polys[0]


def main():
    print("Flow direction/accumulation yukleniyor...")
    grid = Grid.from_raster(f"{OUT_DIR}/flow_direction_D8.tif")
    fdir = grid.read_raster(f"{OUT_DIR}/flow_direction_D8.tif")
    with rasterio.open(f"{OUT_DIR}/flow_accumulation.tif") as src:
        acc_arr = src.read(1)

    stations = pd.read_csv(f"{PROCESSED_DIR}/station_locations.csv")
    natural = gpd.read_file(f"{OUT_DIR}/basin_delineated_natural.shp").geometry.iloc[0]

    print("Istasyon catchment'leri hesaplaniyor (FP1, FP2, FP3, RG8, RG12)...")
    station_polys = {}
    for sid in DISCHARGE_STATIONS + INFO_ONLY_STATIONS:
        row = stations[stations.station_id == sid].iloc[0]
        x, y = snap_outlet(acc_arr, grid.affine, row["x_utm35n"], row["y_utm35n"])
        poly = delineate(grid, fdir, x, y)
        station_polys[sid] = poly
        print(f"  {sid} ({row['name']}): {poly.area/1e6:.2f} km2 (ham -- ic ice olabilir)")

    # --- RG6 (baraj) yerel havzasi: FP2 + FP3'un GECIP baraja ulastigi bolge ---
    # RG6'nin D8-snap noktasindan hesaplanan ham catchment (natural, 174.66 km2)
    # FP2 VE FP3'u TAMAMEN icine aliyor (%100 ic ice) -- bu, FP2/FP3'un
    # suyunun fiziksel olarak bu bolgeden GECEREK baraja ulastiginin D8
    # kaniti (bkz. modul docstring, "REVIZYON NOTU"). FP2+FP3 cikarildiktan
    # sonra kalan "natural - FP2 - FP3", RG6'nin KENDINE OZGU (yerel) katki
    # alani -- eskiden "FP4" diye adlandirilan, artik istasyon degil, bolge.
    station_polys[MASTER_OUTLET] = natural.difference(
        station_polys["FP2"].union(station_polys["FP3"])
    )
    print(f"  {MASTER_OUTLET} (baraj yerel havzasi = natural - FP2 - FP3): "
          f"{station_polys[MASTER_OUTLET].area/1e6:.2f} km2")

    # --- RG12 -> FP1'e belgelenmis override union (bkz. RG12_MERGED_INTO yorumu) ---
    rg12_row = stations[stations.station_id == "RG12"].iloc[0]
    rg12_x, rg12_y = snap_outlet(acc_arr, grid.affine, rg12_row["x_utm35n"], rg12_row["y_utm35n"])
    rg12_poly = delineate(grid, fdir, rg12_x, rg12_y)
    dist = station_polys[RG12_MERGED_INTO].distance(rg12_poly)
    print(f"  RG12 kendi D8 catchment'i: {rg12_poly.area/1e6:.2f} km2, "
          f"{RG12_MERGED_INTO}'e olan mesafe: {dist:.0f}m (D8 ile fiziksel baglanti YOK -- "
          f"duz DEM plato artefakti, bkz. modul docstring) -- union ile birlestiriliyor")
    station_polys[RG12_MERGED_INTO] = station_polys[RG12_MERGED_INTO].union(rg12_poly)

    # --- RG8 -> FP1'e belgelenmis union (bkz. RG8_MERGED_INTO yorumu) ---
    rg8_row = stations[stations.station_id == "RG8"].iloc[0]
    rg8_x, rg8_y = snap_outlet(acc_arr, grid.affine, rg8_row["x_utm35n"], rg8_row["y_utm35n"])
    rg8_poly = delineate(grid, fdir, rg8_x, rg8_y)
    print(f"  RG8 kendi D8 catchment'i: {rg8_poly.area/1e6:.2f} km2 -- flow-path izlemeyle "
          f"FP1 agina dogru sekilde baglandigi dogrulandi (RG12'deki gibi bir artefakt YOK), "
          f"gorsel/kategori tutarliligi icin {RG8_MERGED_INTO}'e union ediliyor")
    station_polys[RG8_MERGED_INTO] = station_polys[RG8_MERGED_INTO].union(rg8_poly)

    # --- Oncelik sirasi: en KUCUK/en spesifik catchment kazanir (ic ice olma ihtimaline karsi) ---
    priority_order = sorted(
        DISCHARGE_STATIONS + [MASTER_OUTLET] + INFO_ONLY_STATIONS,
        key=lambda sid: station_polys[sid].area,
    )
    print(f"\nOncelik sirasi (kucukten buyuge): {priority_order}")

    # --- Ic ice gecmis (nested) catchment'leri, buyukten cikararak ayristir ---
    # NOT: assign_point() zaten "en kucuk kazanir" ile dogru calisir (nesting
    # olsa bile), ama SAKLANAN poligonlar (subbasins_by_station.shp) ic ice
    # kalirsa alan toplamlari yanlis olur (ayni piksel iki poligonda sayilir).
    # Kontrol: her cift arasinda >0 kesisim var mi -- varsa buyuk olandan
    # kucuk olani cikar (kucuk zaten kendi orijinal halinde kalir).
    for i in range(len(priority_order)):
        for j in range(i + 1, len(priority_order)):
            small, big = priority_order[i], priority_order[j]
            inter = station_polys[small].intersection(station_polys[big]).area
            if inter > 1.0:  # 1 m2 tolerans (float/rasterizasyon gurultusu)
                print(f"  UYARI: {small} ({station_polys[small].area/1e6:.2f} km2) "
                      f"{big} icinde ic ice ({inter/1e6:.4f} km2 kesisim) -- "
                      f"{big}'den cikarilarak ayristiriliyor")
                station_polys[big] = station_polys[big].difference(station_polys[small])

    # --- Resmi havza sinirindaki bosluklari doldur ---
    # D8-turetilmis 5 catchment resmi siniri (basin_official, 260.7 km2) tam
    # kaplamiyor -- aralarinda (cogu D8/resmi sinir kucuk uyumsuzluklarindan
    # kaynaklanan) bosluklar var; en buyugu FP1-RG12 arasinda 14.47 km2
    # (kullanicinin fark ettigi "kirmizi alaninin ici bos" sorunu -- ayni
    # DEM-kosesi/duz-plato sinirlamasi RG12 icin de gecerli, bkz. yukarida).
    # Her boslugu EN UZUN ORTAK SINIRA sahip oldugu komsu alt-havzaya
    # birlestiriyoruz (en mantikli komsuluk varsayimi); hicbir alt-havzaya
    # dokunmayan (yalitilmis) parcalar HARICI olarak birakilir.
    basin_official_geom = gpd.read_file(f"{OUT_DIR}/basin_official.shp").geometry.iloc[0]
    covered = None
    for sid in priority_order:
        covered = station_polys[sid] if covered is None else covered.union(station_polys[sid])
    gaps = basin_official_geom.difference(covered)
    gap_parts = list(gaps.geoms) if gaps.geom_type == "MultiPolygon" else [gaps]
    print(f"\nResmi sinir icinde doldurulacak bosluk sayisi: {len(gap_parts)} "
          f"(toplam {gaps.area/1e6:.2f} km2)")
    for part in gap_parts:
        if part.area < 1:  # <1 m2 rasterizasyon gurultusu, atla
            continue
        best_sid, best_shared = None, 0.0
        for sid in priority_order:
            shared = part.boundary.intersection(station_polys[sid].boundary).length
            if shared > best_shared:
                best_sid, best_shared = sid, shared
        if best_sid is not None and best_shared > 0:
            if part.area > 100_000:  # >0.1 km2 olanlari raporla, kucukleri sessizce ekle
                print(f"  Bosluk ({part.area/1e6:.2f} km2, merkez={part.centroid}) "
                      f"-> {best_sid}'e eklendi (ortak sinir={best_shared:.0f}m)")
            station_polys[best_sid] = station_polys[best_sid].union(part)
        else:
            print(f"  UYARI: bosluk ({part.area/1e6:.2f} km2) hicbir alt-havzaya "
                  f"dokunmuyor, harici birakildi")

    # --- Resmi sinirin DISINA tasan kucuk parcalari kirp ---
    # D8 fiziksel hesabi ile elle dijitalleştirilmis basin_official.shp
    # birebir ortusmuyor -- 63 kucuk parcada (toplam ~5 km2, en buyugu 0.8 km2)
    # D8 poligonlari resmi sinirin disina tasiyordu. Haritada bu, renkli
    # dolgunun kesikli referans cizgisinin disina cikip "oturmamis" gorunmesine
    # neden oluyordu (kullanici ekran goruntusuyle fark etti). basin_official
    # dort bagimsiz yontemle dogrulanmis nihai referans oldugu icin (bkz.
    # README "Alt-havza modeli"), her alt-havzayi onunla kesistirerek disariya
    # tasan kucuk parcalari temizliyoruz. Test edildi: hicbir alt-havza
    # parcalanmiyor (hepsi tek Polygon kaliyor).
    for sid in priority_order:
        clipped = station_polys[sid].intersection(basin_official_geom)
        if clipped.geom_type == "MultiPolygon":
            parts = sorted(clipped.geoms, key=lambda p: -p.area)
            print(f"  UYARI: {sid} kirpma sonrasi {len(parts)} parcaya bolundu, "
                  f"en buyugu tutuluyor ({parts[0].area/1e6:.4f} km2, "
                  f"digerleri toplam {sum(p.area for p in parts[1:])/1e6:.4f} km2 atildi)")
            clipped = parts[0]
        diff = (station_polys[sid].area - clipped.area) / 1e6
        if diff > 0.001:
            print(f"  {sid}: resmi sinir disina tasan {diff:.3f} km2 kirpildi")
        station_polys[sid] = clipped

    kind_map = {sid: ("debi_istasyonu" if sid in DISCHARGE_STATIONS
                       else "baraj_yerel_havza" if sid == MASTER_OUTLET
                       else "harici_yagis_istasyonu")
                for sid in priority_order}

    def assign_point(pt):
        for sid in priority_order:
            if station_polys[sid].contains(pt):
                return sid
        return None

    # --- Alt-havza shapefile'i ---
    name_map = dict(zip(stations.station_id, stations.name))
    rows = [{"station_id": sid, "name": name_map.get(sid, sid), "kind": kind_map[sid],
             "area_km2": round(station_polys[sid].area / 1e6, 2)} for sid in priority_order]
    gdf_sub = gpd.GeoDataFrame(rows, geometry=[station_polys[s] for s in priority_order], crs="EPSG:32635")
    gdf_sub.to_file(f"{OUT_DIR}/subbasins_by_station.shp")
    print(f"\nYazildi: {OUT_DIR}/subbasins_by_station.shp")
    print(gdf_sub.drop(columns='geometry').to_string(index=False))

    # --- Dogrulama: artik hicbir cift ortusmuyor mu? ---
    max_overlap = 0.0
    for i in range(len(priority_order)):
        for j in range(i + 1, len(priority_order)):
            a, b = priority_order[i], priority_order[j]
            ov = station_polys[a].intersection(station_polys[b]).area / 1e6
            max_overlap = max(max_overlap, ov)
    print(f"Dogrulama: en buyuk kalan ikili kesisim = {max_overlap:.6f} km2 "
          f"({'TEMIZ' if max_overlap < 0.001 else 'UYARI: hala ortusme var!'})")

    # --- Dere agini istasyona ata (segment merkez noktasindan) ---
    print("\nDere agi segmentleri istasyonlara atanıyor...")
    with open(f"{OUT_DIR}/stream_network.geojson", encoding="utf-8") as f:
        streams = gpd.GeoDataFrame.from_features(json.load(f)["features"], crs="EPSG:32635")
    streams["station_id"] = streams.geometry.apply(lambda g: assign_point(g.centroid))
    streams["kind"] = streams["station_id"].map(kind_map)
    streams["source"] = "d8"
    print(streams["station_id"].value_counts(dropna=False))

    # --- RG12-FP1 baglantisini ANA DERE AGINA katistir ---
    # Onceden ayri bir dosyada (streams_rg12_connector.shp) ve haritada kesikli
    # cizgiyle tutuluyordu -- kullanici bunun "yagis-akis-dere iliskisi" icin
    # gercek bir dere kolu gibi ISLEV GORMESI gerektigini belirtti (kesikli
    # cizgi kopuk/belirsiz izlenimi veriyordu). Artik ana agin bir PARCASI:
    # ayni sekilde (duz cizgi) cizilir, station_id=FP1 ile aynı sekilde
    # islenir -- ama source='dem_leastcost' etiketiyle NEREDEN geldigi
    # (gercek D8 degil, DEM yukseklik-maliyetli en-kisa-yol) izlenebilir kalir.
    connector = gpd.read_file(f"{OUT_DIR}/streams_rg12_connector.shp")
    connector["kind"] = connector["station_id"].map(kind_map)
    connector = connector[["station_id", "kind", "source", "geometry"]]
    streams = pd.concat([streams, connector], ignore_index=True)
    print(f"  RG12-FP1 baglantisi ana dere agina eklendi (source=dem_leastcost)")

    streams.to_file(f"{OUT_DIR}/streams_by_station.shp")
    print(f"Yazildi: {OUT_DIR}/streams_by_station.shp")

    # --- Dolin noktalarini istasyona ata ---
    print("\nDolin noktalari istasyonlara atanıyor...")
    dolines = gpd.read_file(f"{PROCESSED_DIR}/doline_points.shp")
    dolines["station_id"] = dolines.geometry.apply(assign_point)
    dolines["kind"] = dolines["station_id"].map(kind_map)
    dolines["assign_method"] = "d8_surface"

    # --- Kartepe/RG9 karst override ---
    # RG9 bolgesi D8 yuzey akisiyla resmi havzaya baglanmiyor (dem_flow_pipeline.py:
    # sub-catchment kabul testinde %0 kesisim, "kuzeye/harita disina akiyor" diye
    # REDDEDILMISTI). AMA kullanicinin saha bilgisine gore buradaki dolinler
    # (karstik obruklar) yeraltindaki magara sistemlerine aciliyor, ve o su
    # yeralti yoluyla (yuzey D8 topografyasini ATLAYARAK) baraja ulasiyor.
    # Bu yuzden SADECE bu dolin noktalarini RG6'ya (baraj yerel havzasi) atiyoruz
    # -- cevredeki genel yamac yuzey akisi (D8'in geri kalani) hala "harici/RG9"
    # olarak kaliyor, cunku iddia spesifik olarak dolin/magara mekanizmasiyla
    # ilgili, tum RG9 alaninin yuzey akisiyla degil.
    KARST_OVERRIDE_STATION, KARST_OVERRIDE_RADIUS_M, KARST_OVERRIDE_TARGET = "RG9", 3000, MASTER_OUTLET
    rg9_row = stations[stations.station_id == KARST_OVERRIDE_STATION].iloc[0]
    from shapely.geometry import Point as _Point
    rg9_pt = _Point(rg9_row["x_utm35n"], rg9_row["y_utm35n"])
    karst_mask = dolines["station_id"].isna() & (dolines.geometry.distance(rg9_pt) < KARST_OVERRIDE_RADIUS_M)
    print(f"  Karst override (RG9 magara baglantisi, D8 DEGIL): "
          f"{karst_mask.sum()} dolin -> {KARST_OVERRIDE_TARGET}'e atandi")
    dolines.loc[karst_mask, "station_id"] = KARST_OVERRIDE_TARGET
    dolines.loc[karst_mask, "kind"] = "baraj_yerel_havza"
    dolines.loc[karst_mask, "assign_method"] = "karst_magara_override"

    print(dolines["station_id"].value_counts(dropna=False))

    # --- Her dolini kendi alt-havzasindaki EN YAKIN dere segmentine bagla ---
    # ("dolinlerin baglantili oldugu/aktigi dereler" -- kullanici istegi).
    # Karst-override dolinleri (Kartepe/RG9) icin bu mesafe buyuk cikar
    # (~7-8km) -- bu BEKLENEN: yuzeyde gorunen en yakin RG6 deresi bu, ama
    # gercek baglanti yeralti magara sistemi uzerinden, yuzey mesafesi degil.
    name_map = {"FP1": "Kirazdere", "FP2": "Kazandere", "FP3": "Serindere",
                "RG6": "Baraj Yerel Havzasi", "RG8": "Aytepe"}
    streams["seg_id"] = streams.index
    nearest_seg, nearest_dist, dere_adi = [], [], []
    for _, d in dolines.iterrows():
        sub_streams = streams[streams.station_id == d.station_id]
        if len(sub_streams) == 0:
            nearest_seg.append(None); nearest_dist.append(None); dere_adi.append(None)
            continue
        dists = sub_streams.geometry.distance(d.geometry)
        idx = dists.idxmin()
        nearest_seg.append(int(sub_streams.loc[idx, "seg_id"]))
        nearest_dist.append(round(float(dists.loc[idx]), 1))
        dere_adi.append(name_map.get(d.station_id, d.station_id))
    dolines["dere_adi"] = dere_adi
    dolines["nearest_seg_id"] = nearest_seg
    dolines["dist_to_stream_m"] = nearest_dist

    # Resmi havza sinirinin (basin_official) COGRAFI olarak icinde mi --
    # karst-override dolinleri (Kartepe/RG9) burada False cikar: yuzeysel
    # olarak havza disinda kalirlar, station_id atamasi (RG6) SADECE
    # yeralti/magara baglantisi yuzunden -- bu sutun o ayrimi acikca isaretler.
    basin_official_geom = gpd.read_file(f"{OUT_DIR}/basin_official.shp").geometry.iloc[0]
    dolines["inside_basin_official"] = dolines.geometry.within(basin_official_geom)

    out_csv = dolines.drop(columns="geometry")
    out_csv.to_csv(f"{PROCESSED_DIR}/doline_points_enriched.csv", index=False)
    dolines.to_file(f"{PROCESSED_DIR}/doline_points_enriched.shp")
    print(f"Yazildi: {PROCESSED_DIR}/doline_points_enriched.csv (+.shp)")

    # --- Akis yonlendirme (flow routing) tablosu ---
    # Kullanici istegi: "FP1 FP2 FP3'un baraj golune akitildigini de
    # tanimlamak gerekiyor". Kullanicinin duzeltmesi (onemli): FP1, FP2, FP3
    # BIRBIRINDEN BAGIMSIZ 3 AYRI KOL -- her biri kendi yatagindan akar,
    # BIRBIRIYLE YUKARIDA BIRLESMEZLER. Tek ortak nokta BARAJIN KENDISI
    # (RG6/rezervuar). (NOT: D8 iz surme FP2/FP3'un baraja ~3km kala ayni
    # hucrelerden gectigini gosteriyordu -- ama bu, su dengesi/modelleme
    # amaci icin YANLIS bir izlenim veriyordu; kullanicinin sadelestirilmis
    # modeli esas alindi: her kol bagimsiz, sadece barajda "bulusuyorlar".)
    routing_rows = [
        {"from_station": "FP1", "to_station": "RG6", "mechanism": "harici_sistem_muhtemelen_tunel",
         "note": "Kirazdere -- ayri/harici sistem (muhtemelen tunel), dogrudan baraja"},
        {"from_station": "FP2", "to_station": "RG6", "mechanism": "bagimsiz_kol",
         "note": "Kazandere -- kendi (yesil) alt-havzasindan akar, FP3 ile yukarida birlesmez, baraja ulasir"},
        {"from_station": "FP3", "to_station": "RG6", "mechanism": "bagimsiz_kol_yerel_havzadan_gecis",
         "note": "Serindere -- kendi (mavi) alt-havzasindan akar, RG6 yerel havzasindan (sari bolge) gecerek baraja ulasir"},
        {"from_station": "RG6_yerel_havza", "to_station": "RG6", "mechanism": "d8_yuzey_yerel",
         "note": "Kendi yerel akisi, dogrudan baraja"},
    ]
    pd.DataFrame(routing_rows).to_csv(f"{OUT_DIR}/subbasin_flow_routing.csv", index=False)
    print(f"Yazildi: {OUT_DIR}/subbasin_flow_routing.csv "
          f"(FP1/FP2/FP3 -> RG6, bagimsiz kollar, sadece barajda bulusuyor)")

    # --- Gorsel dogrulama ---
    with rasterio.open(f"{OUT_DIR}/DEM_UTM35N.tif") as src:
        dem = src.read(1)
        dem = np.ma.masked_equal(dem, src.nodata)
        extent = [src.bounds.left, src.bounds.right, src.bounds.bottom, src.bounds.top]

    basin_official = gpd.read_file(f"{OUT_DIR}/basin_official.shp")
    colors = {"FP1": "#e6194B", "FP2": "#3cb44b", "FP3": "#4363d8",
              "RG6": "#f58231", "RG8": "#911eb4"}

    fig, ax = plt.subplots(figsize=(14, 11))
    ax.imshow(dem, cmap="Greys", extent=extent, alpha=0.5)
    for sid in priority_order:
        gpd.GeoSeries([station_polys[sid]], crs="EPSG:32635").plot(
            ax=ax, facecolor=colors[sid], edgecolor="black", linewidth=1, alpha=0.55,
            label=f"{sid} {name_map.get(sid,'')} ({station_polys[sid].area/1e6:.1f} km2, {kind_map[sid]})"
        )
    basin_official.boundary.plot(ax=ax, color="black", linewidth=1.5, linestyle="--",
                                  label="Resmi havza siniri (260.7 km2)")
    for sid in priority_order:
        sub = streams[streams.station_id == sid]
        sub.plot(ax=ax, color=colors[sid], linewidth=0.8)
    unassigned_streams = streams[streams.station_id.isna()]
    unassigned_streams.plot(ax=ax, color="lightgray", linewidth=0.4)

    for sid in priority_order:
        sub = dolines[dolines.station_id == sid]
        ax.scatter(sub.geometry.x, sub.geometry.y, color=colors[sid], edgecolor="black",
                   marker="o", s=45, zorder=6)
    unassigned_dolines = dolines[dolines.station_id.isna()]
    ax.scatter(unassigned_dolines.geometry.x, unassigned_dolines.geometry.y, color="lightgray",
               edgecolor="black", marker="o", s=45, zorder=6, label="Dolin - atanmamis (harici)")

    stations_plot = stations[stations.station_id.isin(priority_order + ["RG6", "RG7", "RG9", "RG10", "RG11", "RG12"])]
    ax.scatter(stations_plot["x_utm35n"], stations_plot["y_utm35n"], c="white", s=50,
               edgecolor="black", marker="s", zorder=7)
    for _, row in stations_plot.iterrows():
        ax.annotate(row["station_id"], (row["x_utm35n"], row["y_utm35n"]), fontsize=7,
                    xytext=(3, 3), textcoords="offset points", zorder=7)

    ax.set_title("Yuvacik Havzasi - Dere agi ve dolinler, debi istasyonu catchment'lerine gore")
    ax.set_xlabel("UTM 35N Easting (m)")
    ax.set_ylabel("UTM 35N Northing (m)")
    ax.legend(loc="lower left", fontsize=7)
    ax.set_aspect("equal")
    plt.tight_layout()
    out_png = f"{OUT_DIR}/subbasins_by_station_map.png"
    plt.savefig(out_png, dpi=150)
    print(f"\nYazildi: {out_png}")


if __name__ == "__main__":
    main()
