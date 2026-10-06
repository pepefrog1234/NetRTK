#!/bin/bash
# Validation campaign of docs/VALIDATION.md (spec ch.12). Run from the repository root
# after building (build/bin). Synthetic data are generated if missing; the Hong Kong
# SatRef data must be downloaded to data/real/ first (see docs/VALIDATION.md).
#
#   tools/run_validation.sh [sim] [hk1s] [hk5s]      (default: all)
#
# Per case: network solution (nrtk-net), VRS generation (nrtk-gen), leave-one-out
# rover solutions (tools/loo_validate.py). Summaries: results/validation/<case>.txt
set -u
cd "$(dirname "$0")/.."
BIN=build/bin
OUT=results/validation
mkdir -p "$OUT"
cases=("$@")
[ ${#cases[@]} -eq 0 ] && cases=(sim hk1s hk5s)

run_case() {   # name conf systems reset job...
    local name=$1 conf=$2 sys=$3 reset=$4
    shift 4
    local log=$OUT/${name}.txt
    echo "######## $name ($conf)" | tee "$log"
    $BIN/nrtk-net -c "$conf" -o "results/${name}_net" -v 1 2>&1 | tail -n 25 | tee -a "$log"
    $BIN/nrtk-gen -c "$conf" -v 1 2>&1 | tail -n 6 | tee -a "$log"
    for job in "$@"; do
        python3 tools/loo_validate.py -c "$conf" -j "$job" --systems "$sys" --reset "$reset" \
            -o "results/${name}_${job}" 2>&1 | tee -a "$log"
    done
}

for c in "${cases[@]}"; do
    case $c in
    sim)
        for s in tw_quiet tw_day; do
            [ -f data/sim/$s/brdc.nav ] || $BIN/nrtk-sim conf/sim/$s.yaml
            # rover GPS+Galileo+QZSS: RTKLIB-EX pairs BDS B1I with B2b in l1+l2 mode,
            # the simulated BDS signals are B1I/B3I/B1C/B2a (BDS checked at the
            # observation layer only)
            run_case sim_$s conf/sim/${s}_network.yaml GEJ 600 LOO_TWE5
        done ;;
    # rover GPS+QZSS: BDS is single-frequency on the HK receivers and RTKLIB-EX
    # fixes Galileo wrongly on this data (docs/VALIDATION.md)
    hk1s) run_case hk1s conf/hk/hk_1s_sparse.yaml GJ 600 LOO_HKSC LOO_HKKT LOO_HKST ;;
    hk5s) run_case hk5s conf/hk/hk_5s_day_sparse.yaml GJ 1800 LOO_HKSC LOO_HKKT LOO_HKST ;;
    *) echo "unknown case $c" >&2; exit 2 ;;
    esac
done
