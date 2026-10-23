set -euo pipefail

###======================== Parameter Config ============================
SAMPLE_ID="S8-noXY"
# Upper bound on the number of clusters. Empty clusters are not reported.
CLUSTERS=10
# Binomial density for read counts. The other common choice is beta-binomial.
DENSITY="binomial"
# Independent restarts. The best restart is kept.
RESTARTS=10

###======================== Paths =======================================
ROOT=$(cd "$(dirname "$0")/.." && pwd)
cd "$ROOT"
# One OpenMP thread, matching the single-core timing run.
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"

PYCLONE_VI="${PYCLONE_VI:-pyclone-vi}"
if ! command -v "$PYCLONE_VI" >/dev/null 2>&1; then
    echo "pyclone-vi is not on PATH. Run this inside the workshop image." >&2
    exit 1
fi

# Tumour_content is the sample purity and is an input, not a result.
input_path="$ROOT/data/PyCloneVI_mutation_table.tsv"
output_dir="$ROOT/results/pyclonevi"
mkdir -p "$output_dir"
if [[ ! -f "$input_path" ]]; then
    echo "Missing PyClone-VI input: $input_path" >&2
    exit 1
fi

###======================== Run PyClone-VI ==============================
# fit writes an HDF5 file. write-results-file turns it into a TSV with
# cluster_id and cellular_prevalence for every input mutation.
start=$(date -u +%Y-%m-%dT%H:%M:%SZ)
SECONDS=0
"$PYCLONE_VI" fit \
    -i "$input_path" \
    -o "$output_dir/${SAMPLE_ID}.h5" \
    -c "$CLUSTERS" \
    -d "$DENSITY" \
    -r "$RESTARTS"
"$PYCLONE_VI" write-results-file \
    -i "$output_dir/${SAMPLE_ID}.h5" \
    -o "$output_dir/${SAMPLE_ID}_results.tsv"
end=$(date -u +%Y-%m-%dT%H:%M:%SZ)

###======================== Save runtime ================================
cat > "$output_dir/runtime.txt" <<EOF
tool=pyclonevi
start=${start}
end=${end}
elapsed_sec=${SECONDS}
EOF

###======================== Extract clones ==============================
# clones.tsv columns: purity, n_clones, clone, ccf.
# purity is NA because tumour_content was supplied as input.
# n_clones is the number of distinct cluster_id values.
# PyClone-VI reports cellular_prevalence conditional on that purity, so
# ccf = cellular_prevalence * truth purity from data/sample.tsv.
extract_clones() {
    python3 - "$output_dir/${SAMPLE_ID}_results.tsv" "$output_dir/clones.tsv" "$ROOT/data/sample.tsv" <<'PY'
import csv
import sys

result_path, out_path, sample_path = sys.argv[1], sys.argv[2], sys.argv[3]

with open(sample_path, newline="") as handle:
    sample = next(csv.DictReader(handle, delimiter="\t"))
truth_purity = float(sample["purity"])
if not (truth_purity > 0):
    sys.exit("data/sample.tsv purity is not positive.")

# Collect the cellular-prevalence values seen inside each cluster.
grouped = {}
with open(result_path, newline="") as handle:
    for row in csv.DictReader(handle, delimiter="\t"):
        grouped.setdefault(row["cluster_id"], set()).add(row["cellular_prevalence"])

def cluster_key(value):
    # Numeric ids sort as 0, 1, 2 rather than "0", "1", "10".
    return int(value) if value.isdigit() else value

rows = []
for cluster_id in sorted(grouped, key=cluster_key):
    prevalences = grouped[cluster_id]
    ccf = float(next(iter(prevalences))) * truth_purity
    rows.append((cluster_id, f"{ccf:.8g}"))

with open(out_path, "w", newline="") as handle:
    writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
    writer.writerow(["purity", "n_clones", "clone", "ccf"])
    if not rows:
        writer.writerow(["NA", 0, "", ""])
    else:
        for clone, ccf in rows:
            writer.writerow(["NA", len(rows), clone, ccf])
print(f"Wrote {out_path}")
PY
}
extract_clones

###======================== Plot VAF ====================================
# PyClone-VI does not return a parametric SFS density. This histogram is the
# observed VAF, stacked by assigned cluster. Cluster 0 is the bottom of each
# bar and higher cluster ids sit above it.
plot_sfs() {
    local pdf_path="$output_dir/${SAMPLE_ID}_sfs.pdf"
    if ! Rscript - "$input_path" "$output_dir/${SAMPLE_ID}_results.tsv" "$pdf_path" "$SAMPLE_ID" <<'RS'
args <- commandArgs(trailingOnly = TRUE)
input_path <- args[[1]]
result_path <- args[[2]]
pdf_path <- args[[3]]
sample_id <- args[[4]]

suppressPackageStartupMessages(library(ggplot2))
counts <- read.delim(input_path, stringsAsFactors = FALSE)
fits <- read.delim(result_path, stringsAsFactors = FALSE)
# Keep a mutation only when both the read counts and the cluster assignment exist.
merged <- merge(counts, fits, by = "mutation_id", all = FALSE)
depth <- merged$ref_counts + merged$alt_counts
merged <- merged[is.finite(depth) & depth > 0, , drop = FALSE]
if (nrow(merged) == 0L) stop("No PyClone-VI mutations left to plot.")
merged$vaf <- merged$alt_counts / depth
# Numeric ids sort as 0, 1, 2 rather than the order they appear in the file.
# position_stack(reverse = TRUE) puts the first level at the bottom.
merged$cluster <- factor(merged$cluster_id, levels = sort(unique(merged$cluster_id)))

plot <- ggplot(merged, aes(x = vaf, fill = cluster)) +
    geom_histogram(position = position_stack(reverse = TRUE), bins = 80, linewidth = 0) +
    labs(
        title = paste0(sample_id, " PyClone-VI VAF by cluster"),
        subtitle = "Observed VAF stacked by assigned cluster, cluster 0 at the bottom. PyClone-VI does not return a fitted SFS density.",
        x = "Variant allele frequency",
        y = "Mutations",
        fill = "Cluster"
    ) +
    guides(fill = guide_legend(reverse = TRUE)) +
    theme_bw()
ggsave(pdf_path, plot, width = 10, height = 6)
cat("Wrote", pdf_path, "\n")
RS
    then
        echo "WARNING: PyClone-VI SFS plot failed." >&2
    fi
}
plot_sfs

echo "PyClone-VI finished in ${SECONDS}s."
