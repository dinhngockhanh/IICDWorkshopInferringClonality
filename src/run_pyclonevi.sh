#!/usr/bin/env bash
# PyClone-VI on multiplicity-1 mutations, including non-diploid copy number.
# Cluster timing command: fit -c 10 -d beta-binomial -r 10, then write-results-file.
set -euo pipefail

ROOT=$(cd "$(dirname "$0")/.." && pwd)
cd "$ROOT"
# shellcheck disable=SC1090
eval "$(python3 "$ROOT/src/read_params.py" --format shell)"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"

PYCLONE_VI="${PYCLONE_VI:-pyclone-vi}"
if ! command -v "$PYCLONE_VI" >/dev/null 2>&1; then
    echo "pyclone-vi is not on PATH. Run this inside the workshop image." >&2
    exit 1
fi

input_path="$ROOT/data/prepared/inputs/pyclonevi/${SAMPLE_ID}.tsv"
output_dir="$ROOT/results/pyclonevi"
mkdir -p "$output_dir"
if [[ ! -f "$input_path" ]]; then
    echo "Missing PyClone-VI input: $input_path" >&2
    exit 1
fi

start=$(date -u +%Y-%m-%dT%H:%M:%SZ)
SECONDS=0
"$PYCLONE_VI" fit \
    -i "$input_path" \
    -o "$output_dir/${SAMPLE_ID}.h5" \
    -c "$PYCLONE_CLUSTERS" \
    -d "$PYCLONE_DENSITY" \
    -r "$PYCLONE_RESTARTS"
"$PYCLONE_VI" write-results-file \
    -i "$output_dir/${SAMPLE_ID}.h5" \
    -o "$output_dir/${SAMPLE_ID}_results.tsv"
end=$(date -u +%Y-%m-%dT%H:%M:%SZ)
cat > "$output_dir/runtime.txt" <<EOF
tool=pyclonevi
start=${start}
end=${end}
elapsed_sec=${SECONDS}
EOF
echo "PyClone-VI finished in ${SECONDS}s."
