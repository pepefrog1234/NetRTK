# nrtk — 開源網路 RTK（虛擬參考站）系統

依《VRS 網路 RTK 系統設計規格書》（2026-10-06）實作的網路 RTK 伺服器原型。以
[RTKLIB-EX v2.5.1](https://github.com/rtklibexplorer/RTKLIB) 為底層（git submodule，
**未修改上游原始碼**），C11/C++17 撰寫：

- **libnrtk**：資料擷取、站網分格、網路解算（基站間模糊度固定與殘差萃取）、
  大氣誤差內插（LIM/LSM）、虛擬觀測量合成、RTCM 3 輸出、NTRIP 前端、監控
- **nrtk-gen**：離線 VRS 產生器（RINEX 3.04 + RTCM 3 + 診斷 CSV），支援留一法
- **nrtk-svr**：即時伺服器（NTRIP v1/v2、rover GGA、每連線一個 VRS session、
  `/status` JSON、`/metrics` Prometheus）
- **nrtk-net**：網路解算 CLI（固定率報表、逐衛星殘差）
- **nrtk-sim**：多站 GNSS 觀測量模擬器（含真值），供開發與回歸測試
- **nrtk-tools**、`tools/*.py`：留一法驗證、觀測量層比對、即時測試用 NTRIP client

> 規格書要求避免以 VRS 商標命名產品，專案名稱暫定為 `nrtk`（network RTK）。

## 目前狀態

| 里程碑 | 內容 | 狀態 |
| --- | --- | --- |
| M0 | RTKLIB-EX submodule、CMake、CI 骨架 | 完成（真實站網資料以香港 SatRef 公開資料代替自架站，見下） |
| M1 | 網路解算（EWL/WL 分層、部分固定、殘差輸出） | 完成；兩種引擎可選，`native` 為預設 |
| M2 | 殘差分離、LIM/LSM、品質指標、cellmgr | 完成 |
| M3 | vrsgen、ANTEX/潮汐、RINEX/RTCM 輸出、MSM 往返測試 | 完成 |
| M4 | 驗證戰役 | 合成資料 + 香港真實資料（1 s 兩小時、5 s 全日）已跑，`tools/run_validation.sh` 可重現；結果與未達標項目見 `docs/VALIDATION.md` |
| M5/M6 | 即時核心、NTRIP 前端、監控 | 完成並以重播資料測試：age of corrections p95 0.28 s、連線到首筆 < 0.1 s、100 client 時單核 CPU 6.5%、單站中斷不斷線、即時輸出與離線輸出一致 |
| M7 | 實地外業、正式發布 | 未做（需自架站與實體 rover） |

詳細設計見 [`docs/ARCHITECTURE.md`](docs/ARCHITECTURE.md)，驗證結果與和規格書的差異見
[`docs/VALIDATION.md`](docs/VALIDATION.md)。

## 建置

需求：CMake ≥ 3.20、C/C++17 編譯器、網路（首次設定時以 FetchContent 取得
yaml-cpp 0.9.0 與 GoogleTest 1.17.0；系統已安裝時直接使用）。Python 3（numpy、
matplotlib 可選）用於驗證腳本。

```bash
git submodule update --init
```

```bash
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=RelWithDebInfo
```

```bash
cmake --build build
```

```bash
cd build && ctest --output-on-failure
```

產出在 `build/bin/`：`nrtk-net nrtk-gen nrtk-svr nrtk-sim nrtk-tools nrtk_tests`，以及
RTKLIB-EX 的 `rnx2rtkp convbin str2str rtkrcv pos2kml`（驗證用，與 libnrtk 使用相同的
`NFREQ/NEXOBS` 編譯定義）。

## 快速開始

### 1. 合成資料（不需要任何外部資料）

```bash
./build/bin/nrtk-sim conf/sim/tw_quiet.yaml
```

模擬北台灣 5 站（站距約 40 km，其中一站 600 m 高、一站自由振盪時鐘、兩種天線），
輸出 RINEX 3.04、混合導航檔、合成 ANTEX 與真值。接著：

```bash
./build/bin/nrtk-net -c conf/sim/tw_quiet_network.yaml -o results/net --sat-csv
```

```bash
./build/bin/nrtk-gen -c conf/sim/tw_quiet_network.yaml
```

```bash
python3 tools/loo_validate.py -c conf/sim/tw_quiet_network.yaml -j LOO_TWE5
```

`loo_validate.py` 以 RTKLIB `rnx2rtkp` 對留一站做三組解（對 VRS、對最近實體站、對主站），
加上 `nrtk-tools obscmp` 的觀測量層比對，並逐項對照規格書第 12 章的驗收指標。整個驗證戰役
（合成兩組、香港 1 s 與 5 s）：

```bash
tools/run_validation.sh sim hk1s hk5s
```

### 2. 真實資料（香港 SatRef）

測試期間無法取得台灣的公開密集站網，改用香港地政總署 SatRef（低緯度、多 GNSS、
1 s 公開下載）。資料授權僅限個人非商業內部使用，**不納入版本控制**；設定檔在
`conf/hk/`，資料放在 `data/real/`，以下列指令下載（需要 `crx2rnx`）：

```bash
tools/hkdata/fetch_hk.sh 2026 270
```

接著：

```bash
./build/bin/nrtk-gen -c conf/hk/hk_1s_sparse.yaml
```

### 3. 即時伺服器

```bash
./build/bin/nrtk-svr -c conf/realtime_replay.yaml
```

`stream:` 支援 `ntrip://user:pass@host:port/mount`、`tcpcli://host:port`、`serial://`、
`file://`（RTCM 3），以及測試用的 `rinex://檔案`（依牆鐘時間重播）。rover 端：

```bash
python3 tools/ntrip_client.py --port 2101 --mount VRS_MSM7 --llh 24.78 121.15 250 --seconds 60
```

狀態頁：`curl http://localhost:2101/status`、`curl http://localhost:2101/metrics`。

## 設定檔

所有鍵與預設值見 [`conf/example.yaml`](conf/example.yaml)（規格書第 8 章的範例可直接使用）。

## 目錄

```
apps/        nrtk-sim nrtk-net nrtk-gen nrtk-svr nrtk-tools
include/nrtk 公開標頭（C API；C++ 管線/即時/NTRIP 為 .hpp）
src/         common ingest netpos atmos cellmgr vrsgen rtcm_out config pipeline ntrip_front
tests/       GoogleTest 單元測試、回歸測試（tests/regression）
tools/       驗證與分析腳本（Python）
conf/        範例與測試設定（sim/、hk/）
third_party/RTKLIB  RTKLIB-EX v2.5.1（submodule，未修改）
```

## 授權

BSD-2-Clause（與 RTKLIB-EX 相同）。`src/common/gmf.c` 為 IERS Conventions GMF 的 C 移植，
保留其授權聲明並註明修改處。測試用真實資料不隨專案散布。
