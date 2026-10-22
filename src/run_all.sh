#!/usr/bin/env bash
# Prepare S8 inputs, run the five tools, and extract an assignment table after each one.
# DECODE is last because the default particle schedule can take hours.
# A tool failure is recorded and the remaining tools still run.
set -uo pipefail

ROOT=$(cd "$(dirname "$0")/.." && pwd)
cd "$ROOT"
mkdir -p "$ROOT/results"
rm -f "$ROOT/results/failures.txt"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"
export RGL_USE_NULL="${RGL_USE_NULL:-TRUE}"

python3 "$ROOT/src/00_prepare.py" || exit 1

failed=0
run_tool() {
    local name="$1"
    shift
    echo ""
    echo "=== ${name} ==="
    if "$@"; then
        if python3 "$ROOT/src/extract_results.py" --tool "$name"; then
            echo "Extracted ${name}"
        else
            echo "EXTRACT FAILED ${name}" | tee -a "$ROOT/results/failures.txt"
            failed=1
        fi
    else
        echo "RUN FAILED ${name}" | tee -a "$ROOT/results/failures.txt"
        failed=1
    fi
}

run_tool pyclonevi bash "$ROOT/src/run_pyclonevi.sh"
run_tool mobster Rscript "$ROOT/src/run_mobster.R"
run_tool sciclone Rscript "$ROOT/src/run_sciclone.R"
run_tool phylowgs bash "$ROOT/src/run_phylowgs.sh"
run_tool decode Rscript "$ROOT/src/run_decode.R"

echo ""
python3 "$ROOT/src/score.py"

if [[ "$failed" -ne 0 ]]; then
    echo "One or more tools failed. See results/failures.txt." >&2
    exit 1
fi
echo "All requested tools finished."
