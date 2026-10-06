#!/bin/sh
# Download the Hong Kong SatRef data used by conf/hk/*.yaml into data/real/.
#
# The data are published by the Lands Department, HKSAR (https://www.geodetic.gov.hk)
# for personal, non-commercial internal use; redistribution needs their consent.
# Keep them out of version control (data/real/ is in .gitignore).
#
# Needs curl, gzip, python3 and crx2rnx (Hatanaka decompression, RNXCMP:
# https://terras.gsi.go.jp/ja/crx2rnx.html) in PATH, or CRX2RNX=/path/to/crx2rnx.
#
#   tools/hkdata/fetch_hk.sh [YYYY DDD]        (default 2026 270 = 2026-09-27)
#
# Result (names as in conf/hk/*.yaml):
#   data/real/hk/YYYY_DDD/5s/<STA>00HKG_R_YYYYDDD0000_01D_05S_MO.rnx   (24 hourly files merged)
#   data/real/hk/YYYY_DDD/1s/<STA>00HKG_R_YYYYDDD0400_02H_01S_MO.rnx   (04-06 UTC)
#   data/real/hk/YYYY_DDD/30s/<STA>00HKG_R_YYYYDDD0000_01D_30S_MO.rnx
#   data/real/nav/BRDC00IGS_R_YYYYDDD0000_01D_MN.rnx                  (IGS combined broadcast, BKG)
#   data/real/atx/igs20.atx
set -e
Y=${1:-2026}; D=${2:-270}
STATIONS="hksl hkkt hktk hkws hkoh hklm hkst hksc"
HOURS1S="04 05"
CRX2RNX=${CRX2RNX:-crx2rnx}
HERE=$(cd "$(dirname "$0")" && pwd)
ROOT=$(cd "$HERE/../.." && pwd)
OUT=$ROOT/data/real
TMP=$OUT/hk/${Y}_${D}/download
mkdir -p "$TMP" "$OUT/nav" "$OUT/atx" "$OUT/hk/${Y}_${D}/1s" "$OUT/hk/${Y}_${D}/5s" "$OUT/hk/${Y}_${D}/30s"
command -v "$CRX2RNX" >/dev/null || { echo "crx2rnx not found (set CRX2RNX=...)" >&2; exit 1; }

fetch() {   # url dest
    [ -s "$2" ] && gzip -t "$2" 2>/dev/null && return 0
    curl -fsS -m 300 --retry 3 -o "$2.part" "$1" && gzip -t "$2.part" && mv "$2.part" "$2"
}
crx() {     # file.crx.gz -> file.rnx
    gzip -dc "$1" | "$CRX2RNX" - > "${1%.crx.gz}.rnx" || [ $? -eq 2 ]    # 2: warnings only
}

for s in $STATIONS; do
    S=$(echo $s | tr a-z A-Z)
    base=https://rinex.geodetic.gov.hk/rinex3/$Y/$D/$s
    echo "== $S"
    f=${S}00HKG_R_${Y}${D}0000_01D_30S_MO.crx.gz
    fetch "$base/30s/$f" "$TMP/$f" && crx "$TMP/$f" && mv "$TMP/${f%.crx.gz}.rnx" "$OUT/hk/${Y}_${D}/30s/"
    parts=""
    for h in 00 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20 21 22 23; do
        f=${S}00HKG_R_${Y}${D}${h}00_01H_05S_MO.crx.gz
        fetch "$base/5s/$f" "$TMP/$f" && crx "$TMP/$f" && parts="$parts $TMP/${f%.crx.gz}.rnx"
    done
    python3 "$HERE/merge_rnx3.py" "$OUT/hk/${Y}_${D}/5s/${S}00HKG_R_${Y}${D}0000_01D_05S_MO.rnx" $parts
    parts=""
    for h in $HOURS1S; do
        f=${S}00HKG_R_${Y}${D}${h}00_01H_01S_MO.crx.gz
        fetch "$base/1s/$f" "$TMP/$f" && crx "$TMP/$f" && parts="$parts $TMP/${f%.crx.gz}.rnx"
    done
    python3 "$HERE/merge_rnx3.py" "$OUT/hk/${Y}_${D}/1s/${S}00HKG_R_${Y}${D}0400_02H_01S_MO.rnx" $parts
done

f=BRDC00IGS_R_${Y}${D}0000_01D_MN.rnx.gz
fetch "https://igs.bkg.bund.de/root_ftp/IGS/BRDC/$Y/$D/$f" "$TMP/$f"
gzip -dc "$TMP/$f" > "$OUT/nav/${f%.gz}"
[ -s "$OUT/atx/igs20.atx" ] || curl -fsS -m 300 -o "$OUT/atx/igs20.atx" https://files.igs.org/pub/station/general/igs20.atx
echo "done: $OUT (the download/ folder can be deleted)"
