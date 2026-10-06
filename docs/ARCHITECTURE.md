# 架構與演算法

本文件說明 nrtk 的模組、資料流與演算法，並逐項對應規格書章節。程式碼註解為英文，
符號約定見 `include/nrtk/common.h` 開頭。

## 1. 資料流

```
參考站 (RTCM3 MSM / RINEX)
   │ ingest      RTKLIB stream.c + rtcm3.c（即時）、rnxctr 逐曆元讀 RINEX（離線，常數記憶體）
   ▼
EpochSync        依 GPS 時間收齊同一曆元；第一站到達後最多等 sync_window_s（0.5 s）
   ▼
prep             訊號規劃（處理頻率 → slot 0..2）、天線 PCO/PCV 與 ARP 偏心移除（零天線、
   │             marker）、以已知座標估接收機鐘差、光行時幾何距離（含 Sagnac、固體潮）
   ▼
netpos ×N        每條基線（所有 ≤ max_baseline_km 的站對）：滑移偵測、EWL/WL、L1 模糊度
   │ (thread pool) LAMBDA、獨立模型驗證、fix-and-hold、雙差殘差 → ∇ΔI₁ / ∇ΔD
   ▼
network closure  三角閉合差檢核（整數一致性），不一致的 (基線, 衛星) 剔除並重新固定
   ▼
NetEpoch ──► EpochBus（最新 + 歷史，即 corrections cache）
                │
                ▼ 每個 NTRIP 連線一條執行緒
cellmgr ──► atmos ──► vrsgen ──► rtcm_out ──► NTRIP（v1/v2）
（Delaunay、主/輔站）（LIM/LSM、對流層模型）（合成觀測量）（MSM4/7、1006、1033）
```

離線（`nrtk-gen`）與即時（`nrtk-svr`）共用 `Network`（`src/pipeline/network.cpp`）與
`VrsStream`，差別只在輸入（`OfflineReader` vs `RealtimeCore` 的串流/重播執行緒）與
輸出（`VrsWriter` 檔案 vs `VrsSession` 的 NTRIP 連線）。

## 2. 模組對應

| 規格書模組 | 實作 | 語言 | 重用的 RTKLIB-EX |
| --- | --- | --- | --- |
| ingest | `src/ingest/ingest.c` | C | stream.c、rtcm3.c、rinex.c（`rnxctr`、`readrnxt`） |
| epoch_sync | `src/pipeline/realtime.cpp`（`EpochSync`） | C++ | — |
| cellmgr | `src/cellmgr/cellmgr.cpp` | C++ | — |
| netpos | `src/netpos/netpos.c` | C | rtkpos.c（`rtklib` 引擎）、lambda.c、`filter()`、ephemeris.c |
| atmos | `src/atmos/atmos.cpp` | C++ | — |
| vrsgen | `src/vrsgen/vrsgen.c` | C | ephemeris（`satpos`）、tides.c（`tidedisp`） |
| rtcm_out | `src/rtcm_out/rtcm_out.c` | C | rtcm3e.c（MSM、1033、星曆）；1006 自行編碼 |
| ntrip_front | `src/ntrip_front/*.cpp` | C++ | — |
| monitor | `RealtimeCore::status_json/metrics` | C++ | — |
| config | `src/config/config.cpp`（yaml-cpp） | C++ | — |
| （另）prep、antenna、tropo | `src/common/*.c` | C | — |

## 3. RTKLIB-EX 的納入方式

`cmake/rtklib.cmake` 直接以 submodule 的原始碼建立靜態庫 `rtklib_ex`，**不**
`add_subdirectory` 上游的 CMakeLists：上游以 `add_definitions` 只在自己的目錄定義
`NFREQ=3 NEXOBS=3 ENA*`，外部程式含入 `rtklib.h` 會看到不同的 `obsd_t/rtcm_t` 佈局；
上游也會把產物寫進 submodule 原始碼樹。這裡把同一組定義設為 PUBLIC，所有使用者以
相同佈局編譯。上游原始碼不做任何修改；需要的修正都在 libnrtk 側處理：

- 接收機 ANTEX：RTKLIB 只讀 GPS 頻率、只有仰角 PCV。`src/common/antenna.c` 讀入每個系統
  /頻率區塊（G01…C07、J05）與方位角 PCV，觀測量在 prep 階段即改正到零天線 marker，
  RTKLIB 濾波器不再套用天線模型。
- MSM 編碼：RTKLIB 不分割超過 64 cells 的訊息，且第 64 個 cell 的資料會被丟掉
  （`gen_msm_sig` 中 `cell>=64` 的 off-by-one）。`rtcm_out` 依衛星分組，每則
  nsat×nsig ≤ 63，並設定 multiple message bit。
- 1006：RTKLIB 固定寫入 DF141=0（實體站）、ITRF 年 0；VRS 需 DF141=1，故自行編碼。

## 4. 網路解算（第 5 章）

輸入為兩站的 prep 結果（同一曆元）。兩種浮點/整數引擎共用下列各層：

1. **滑移偵測**：LLI、訊號碼改變、資料中斷 > 10 s、幾何無關組合跳動（0.05 m，隨間隔放寬）、
   MW 偏離平均值。偵測到的滑移以 LLI 位元注入 RTKLIB 引擎的觀測量。
2. **EWL/WL**：單差 MW 時間平均（預設至少 30 s），對**所有已固定衛星**求加權平均的浮點寬巷
   後取整（不依賴單一參考衛星，避免單顆衛星多路徑偏差污染整個系統）；已固定者持續複核。
3. **L1（窄巷）固定**：候選衛星（寬巷已固定、鎖定 ≥ min_lock、仰角 ≥ AR 遮蔽角）依仰角排序；
   以已持有的整數（RTKLIB 引擎另加寬巷）對浮點解做條件化後，取 DD L1 模糊度子集呼叫
   RTKLIB `lambda()`；條件標準差 > 0.15 周者不參與；ratio ≥ 3；失敗時逐次剔除最低仰角
   （部分固定，最多 max_drop 顆）。首次固定需 ≥ min_fix 個 DD。
4. **獨立模型驗證**（不使用浮點濾波器狀態）：固定後的非頻散項 D 與 Saastamoinen+GMF 雙差
   模型之差，須能以單一相對濕延遲 dzwd（強健加權最小平方、先驗 5 cm、上限 12 cm）解釋到
   max_res（6 cm）以內。L1 錯一周使 D 偏 λ_NL（~11 cm），寬巷錯一周偏 λ_NL·f₂/(f₁−f₂)
   （32–38 cm）；後者會判定寬巷可疑並重新平均。已持有的衛星持續監看：
   - 每曆元：連續 5 曆元超過 2×max_res 即解除固定（寬巷錯、未偵測的滑移）；
   - 平均值：殘差的 120 s 指數平均在仰角 ≥ 25° 時超過 6.5 cm（約 0.6 λ_NL）即解除固定。
     低仰角時多路徑與對流層殘差可達 ±10 cm，與 λ_NL 同量級，錯一周的 L1 整數可能通過
     首次驗證；衛星升高後這類錯誤會以穩定的 ~10 cm 偏差出現（香港實測：HKTK 的 G27 在 11°
     固定錯誤，之後以 −10 cm 持有數小時，每曆元門檻 12 cm 永遠抓不到）；
   - 錨定衛星的多數夥伴失敗時改判錨定衛星。
   解除固定時同時清除該衛星的浮點模糊度（native：清除狀態；rtklib：歸零 bias 狀態使 rtkpos
   重新初始化）。fix-and-hold 已把浮點值釘在錯誤整數上，不清除的話下一次搜尋只會重複同一個錯誤。
5. **fix-and-hold**：僅把**驗證過**的 L1/L2（/L5）整數以虛擬觀測量（變異數 0.001 周²）
   每曆元施加於浮點濾波器；MW 寬巷只用於條件化搜尋，不寫入濾波器。
6. **殘差萃取**：δⱼ = ∇Δ(Φⱼ − ρ) − λⱼ∇ΔNⱼ，以 f₁/f₂ 分離
   ∇ΔI₁ = (δ₁−δ₂)/(γ−1)、∇ΔD = (γδ₁−δ₂)/(γ−1)（換算到 GPS L1 頻率的延遲），
   另輸出雙差碼殘差、模型對流層與固定狀態。連續曆元的 δ₁/δ₂ 跳動 > 7 cm 時解除固定。

整數以「單差整數 + 每系統任意基準」記帳：DD = 兩顆衛星單差整數之差，參考衛星可任意切換
（每系統取仰角最高的已固定衛星）而不需重置濾波器（第 5 章第 6 點）。

### 兩種浮點引擎（`solver.engine`）

- **`rtklib`**：規格書的首選設計。RTKLIB-EX `rtkpos()`，`PMODE_FIXED`、`IONOOPT_EST`
  （每衛星單差電離層）、`TROPOPT_EST[G]`、`ARMODE_OFF`。主站 ZWD 每曆元釘在模型值、只估輔站
  ZWD（短基線上只有相對延遲可觀測，否則共模 ZTD 會漂移到數公尺，實測見 VALIDATION）。
- **`native`**（預設）：規格書第 5 章的備案「自寫 netpos 濾波器」。狀態為一個相對 ZWD 與每衛星
  一個浮點 L1 模糊度，觀測量為扣除已固定寬巷項的單差無電離層相位
  y = SD[Φ_IF − ρ − T_model] − c·f₂/(f₁²−f₂²)·N_WL = λ_NL·b + h·dzwd + 接收機項，
  每系統對樞紐衛星做雙差。碼與電離層不進入 L1 浮點解，因此碼多路徑與電離層活動不會使其偏移。
  香港真實資料上兩者的比較見 VALIDATION，該結果即本專案對「未決事項：模糊度固定引擎」的回答。

## 5. 網路閉合差檢核（新增）

同一曆元中三條基線 i→j→k→i 的雙差殘差相加時，每個觀測量各出現兩次且正負相反，所以
**正確固定時閉合差恆為 0**（與雜訊、大氣、軌道無關）；錯誤的整數組合留下 λ 的整數倍。
各基線各自的參考衛星只造成對所有衛星相同的常數，因此對每個三角形、每個系統取各衛星閉合差
的中位數，偏離 1 cm 以上者為不一致。某 (基線, 衛星) 在其參與的所有三角形都不一致、沒有任何
一致者時即剔除該固定並要求該基線重新解算（`solver.closure_check`，預設開啟）。實測（香港，
全日）抓到的錯誤以「寬巷錯 1 周 + L1 補償」為主，例如 Galileo (E1, E5b) = (+3, +4) 周：
此組合只改變無電離層延遲 −2.9 cm，任何以對流層模型為基礎的驗證都看不出來，但會讓 ∇ΔI₁
錯約 0.6 m。

限制：閉合差只說明一組整數不一致，不說明哪一條錯。投票以多數為準，所以當某站自身資料使
它的多數基線都固定到同一個錯誤整數時（彼此一致），反而是唯一正確的那條被剔除。上一節的
平均殘差監看會在衛星升高後清掉那些錯誤固定，兩者合用；全日 A/B（VALIDATION §3）中合用
的結果在三個留一站都優於只用其一。

## 6. 大氣內插（第 6 章）

`nrtk_atmos_interp()`：每個系統選一顆在最多條 (主站, 輔站) 基線上皆已固定的共同參考衛星，
把各基線值改參考到它（並依基線方向變號）。每顆衛星：

- 電離層：∇ΔI₁ 在主站為原點的水平面（北、東，km）做平面擬合 a·N + b·E（過原點）。
- 非頻散：先扣 ∇ΔT_model（Saastamoinen + GMF，各站實際高度與仰角），殘餘量做同樣的擬合，
  在 VRS 位置加回以 VRS 高度與仰角計算的 ∇ΔT_model。
- n = 2 條基線為 LIM（恰定），n ≥ 3 為 LSM（輸出擬合殘差 RMS），`interp.method` 可強制 LIM。
- 品質指標：每衛星參與基線數、擬合殘差 RMS、是否外插（VRS 不在所用站點的凸包內）；
  每曆元電離層活躍度（各基線 60 s 內 ∇ΔI₁ 變化率最大值，m/min）、對流層殘餘量級。
  電離層活躍度 > 0.1 m/min、可用衛星 < 6、格網外插時標記降級（第 7 章）。

實測發現：三條輔站基線時 LSM 只有 1 個自由度，擬合殘差 RMS 與真實內插誤差幾乎無關
（相關係數 −0.02，合成資料），因此不以它剔除衛星，只作為輸出的品質指標。

## 7. 虛擬觀測量（第 7 章）

`nrtk_vrs_generate()`：以主站的 prep 觀測量（零天線 marker、未動鐘差/模糊度/滑移）為底，

- 真實收訊時間 t − dtr_M，對 VRS（固體潮位移後的位置）做光行時疊代算 ρ_V；
- Lⱼ += (ρ_V − ρ_M + ΔD − μⱼΔI)/λⱼ、Pⱼ += ρ_V − ρ_M + ΔD + μⱼΔI、
  Dⱼ −= (ρ̇_V − ρ̇_M)/λⱼ（μⱼ = (f_L1/fⱼ)²）；
- 單差改正 = 內插雙差 + 每系統偏移量；參考衛星改變時偏移量取延續衛星的平均差，使改正連續；
- lock time 歸零（LLI 滑移位元）的條件：新衛星／新 VRS、網路重新固定、參與內插的基線集合
  改變、套用於各處理頻率的改正 (D − μⱼI) 單曆元跳動超過 3 cm；主站本身的滑移原樣透傳；
- 只輸出有內插值且 VRS 仰角 ≥ 10° 的衛星；`pass_unfixed` 可改為透傳並每曆元歸零 lock；
- 有內插改正的衛星少於 5 顆（網路剛啟動、大範圍失鎖）時改為最近站透傳，回到 8 顆以上才
  恢復 VRS（遲滯，避免頻繁切換身分）；判斷在合成之前，兩種模式不會互相重置 lock 狀態；
- 輸出訊號依訊號規劃（GPS L1C/A、L2W、L2C、L5；Galileo E1、E5a、E5b；BDS B1I、B3I、B1C、B2a；
  QZSS 同 GPS），SNR 沿用主站。

降級：剩餘健康站不足 3 或距最近站 > 70 km → 最近站透傳（station ID 與 1006 改為該站，lock 歸零）。

## 8. 即時核心（第 11 章）

- 每站一條輸入執行緒（RTKLIB stream，斷線 30 s 重連）或一條重播執行緒（`rinex://`）。
  星曆由任一站的 1019/1042/1044/1046 合併；處理執行緒每曆元取寫時複製的快照，同一曆元所有
  計算（基線、VRS）使用同一份星曆。
- `EpochSync`：同一 GPS 時間收齊或等滿視窗即釋出；已釋出時間之後才到的資料計為 late 並丟棄；
  10 s 內有資料的站才列入「應到」。
- 站健康度：有資料、延遲 < 2 s、其基線 60 s 固定率最大值 ≥ 0.9（加 0.15 的遲滯）；啟動前
  180 s 視為健康；離線過的站須回到 0.9 才重新納入（第 11 章）。健康改變即重新 Delaunay 三角化。
- 每個 NTRIP 連線：第一筆有效 GGA（品質 ≥ 1）建立 VRS（1 m 格網）；離 VRS > 5 km 或離開原三角形
  時重建（station ID + 1、重送 1006/1033、lock 歸零）；1006/1033 每 10 曆元一次；連續 10 s
  無法產生（VRS 與最近站都沒有）即關閉連線；15 s 內沒有有效 GGA 由 NTRIP 前端關閉。
- 單站中斷：該站缺席的曆元，若它是主站即由同三角形另一頂點接任主站（lock 全部歸零，
  station ID 不變），服務不中斷；10 s 後判為離線、重新三角化，session 在下一筆 GGA 時發現
  三角形改變並重建 VRS。實測見 VALIDATION §5。
- NTRIP v1 回覆 `ICY 200 OK\r\n`（不加空行；RTKLIB 的 client 只移除這一行），v2 為 HTTP/1.1
  chunked。連線執行緒使用 8 MB 堆疊（macOS 次執行緒預設只有 512 KB）。
- `/status`（JSON）：站、基線、session、處理時間、延遲 p95、session 事件計數；
  `/metrics`（Prometheus）。log 可輸出 JSON lines（`server.log_json`）。

## 9. 已知限制

- 只處理靜態基站；1 Hz 輸出（第 14 章）。
- GLONASS 未納入（規格書列為選配）。
- BDS 需要雙頻 B1I/B3I；香港 GR50 接收機只有 B1I，因此實測未含 BDS。合成資料的 BDS 已驗證到
  觀測量層；RTKLIB-EX 的 rover 在 `l1+l2` 下 BDS 用 B1I+B2b，與 B3I/B2a 無法配對，所以 rover
  解算的驗證不含 BDS（見 VALIDATION）。
- 即時重播模式以資料時間作時鐘；真實串流模式需主機 NTP 校時。
