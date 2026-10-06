#!/bin/bash
# Extract selected stations of a CWA GDMS archive (Taiwan continuous GNSS
# stations, RINEX 2.11 daily files ssssDDD0.YYo.gz, 30 s) into data/real/tw/.
#
#   tools/twdata/extract_gdms.sh GDMSdata.tgz 2026 kyi2 hsin fusn ...
#
# writes data/real/tw/YYYY_DDD/ssssDDD0.YYo.gz for every day found in the
# archive. Broadcast ephemerides are not part of the archive; fetch the IGS
# daily files, e.g. for DOY 271:
#   curl -o data/real/nav/BRDC00IGS_R_20262710000_01D_MN.rnx.gz \
#     https://igs.bkg.bund.de/root_ftp/IGS/BRDC/2026/271/BRDC00IGS_R_20262710000_01D_MN.rnx.gz
# Station coordinates: tools/twdata/station_coords.py (the headers are approximate).
set -eu
cd "$(dirname "$0")/../.."
tgz=$1; year=$2; shift 2
yy=${year:2:2}
tmp=$(mktemp -d)
pat=()
for s in "$@"; do pat+=("GDMSdata/${s}[0-9][0-9][0-9]0.${yy}o.gz"); done
tar xzf "$tgz" -C "$tmp" "${pat[@]}" 2>/dev/null || true
n=0
for f in "$tmp"/GDMSdata/*.gz; do
    [ -e "$f" ] || continue
    b=$(basename "$f")
    doy=${b:4:3}
    mkdir -p "data/real/tw/${year}_${doy}"
    mv "$f" "data/real/tw/${year}_${doy}/"
    n=$((n+1))
done
rm -rf "$tmp"
echo "$n files -> data/real/tw/${year}_DDD/"
