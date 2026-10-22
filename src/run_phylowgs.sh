#!/usr/bin/env bash
# PhyloWGS. Default config/params.toml mode is smoke (-B 5 -s 5 -i 50, one chain).
# Set phylowgs.mode = "full" for burn-in 1000, 2500 MCMC samples, and 1000 MH iterations.
# A full chain on this sample is many hours. Smoke output is not comparable to it.
set -euo pipefail

ROOT=$(cd "$(dirname "$0")/.." && pwd)
cd "$ROOT"
# shellcheck disable=SC1090
eval "$(python3 "$ROOT/src/read_params.py" --format shell)"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"
export TMPDIR="${TMPDIR:-/tmp/workshop}"
export PYTHON2="${PYTHON2:-/opt/python27/bin/python}"
export PHYLOWGS_DIR="${PHYLOWGS_DIR:-/opt/phylowgs}"
mkdir -p "$TMPDIR"

if [[ ! -x "$PYTHON2" ]]; then
    echo "Python 2 interpreter not found at $PYTHON2. Run this inside the workshop image." >&2
    exit 1
fi
if [[ ! -f "$PHYLOWGS_DIR/evolve.py" || ! -f "$PHYLOWGS_DIR/mh.o" ]]; then
    echo "PhyloWGS checkout or mh.o is missing under $PHYLOWGS_DIR." >&2
    exit 1
fi

ssm="$ROOT/data/prepared/inputs/phylowgs/ssm_data.txt"
cnv="$ROOT/data/prepared/inputs/phylowgs/cnv_data.txt"
output_dir="$ROOT/results/phylowgs/${PHYLOWGS_MODE}"
if [[ "$PHYLOWGS_MODE" == "smoke" ]]; then
    rm -rf "$output_dir"
fi
mkdir -p "$output_dir"

ulimit -s unlimited || echo "WARNING: could not set unlimited stack size; PhyloWGS may crash." >&2

start=$(date -u +%Y-%m-%dT%H:%M:%SZ)
SECONDS=0
(
    cd "$output_dir"
    "$PYTHON2" "$PHYLOWGS_DIR/evolve.py" \
        -B "$PHYLOWGS_BURNIN" \
        -s "$PHYLOWGS_MCMC" \
        -i "$PHYLOWGS_MH" \
        --random-seed "$SEED" \
        "$ssm" \
        "$cnv"
)
"$PYTHON2" "$PHYLOWGS_DIR/write_results.py" --include-ssm-names \
    "$SAMPLE_ID" \
    "$output_dir/trees.zip" \
    "$output_dir/summ.json" \
    "$output_dir/mutlist.json" \
    "$output_dir/mutass.json"
printf '%s\n' "$PHYLOWGS_MODE" > "$output_dir/mcmc_mode.txt"
end=$(date -u +%Y-%m-%dT%H:%M:%SZ)
cat > "$output_dir/runtime.txt" <<EOF
tool=phylowgs
mode=${PHYLOWGS_MODE}
burnin=${PHYLOWGS_BURNIN}
mcmc=${PHYLOWGS_MCMC}
mh=${PHYLOWGS_MH}
seed=${SEED}
start=${start}
end=${end}
elapsed_sec=${SECONDS}
EOF
echo "PhyloWGS ${PHYLOWGS_MODE} finished in ${SECONDS}s."
