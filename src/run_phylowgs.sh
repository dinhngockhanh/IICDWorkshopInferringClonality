set -euo pipefail

###======================== Parameter Config ============================
SAMPLE_ID="S8-noXY"
SEED=1001
# Four chains, seeds 1001–1004. Used only when MODE=full.
SEEDS=(1001 1002 1003 1004)
# smoke: one short chain, checks that the install and inputs run.
# full: four chains at the length used for the finished S8-noXY run.
# Change this one word to switch. Smoke output is not comparable to full.
MODE="smoke"
if [[ "$MODE" == "smoke" ]]; then
    # -B burn-in, -s MCMC samples, -i Metropolis-Hastings iterations.
    BURNIN=5
    MCMC=5
    MH=50
    CHAINS=1
elif [[ "$MODE" == "full" ]]; then
    # Four chains. On this sample the wall time is days.
    BURNIN=1000
    MCMC=2500
    MH=5000
    CHAINS=4
else
    echo "MODE must be smoke or full, got ${MODE}" >&2
    exit 1
fi

###======================== Runtime environment =========================
ROOT=$(cd "$(dirname "$0")/.." && pwd)
cd "$ROOT"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"
export TMPDIR="${TMPDIR:-/tmp/workshop}"
# PhyloWGS is a Python 2.7 program. The image builds it at these paths.
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
if [[ "$MODE" == "full" && ! -f "$PHYLOWGS_DIR/multievolve.py" ]]; then
    echo "multievolve.py is missing under $PHYLOWGS_DIR." >&2
    exit 1
fi

###======================== Inputs and output directory =================
# Simple somatic mutations. Columns: id, gene, a, d, mu_r, mu_v.
ssm="$ROOT/data/PhyloWGS_ssm.txt"
# Copy-number records. cell_prev here is the CN fraction, not sample purity.
cnv="$ROOT/data/PhyloWGS_cnv.txt"
output_dir="$ROOT/results/phylowgs/${MODE}"
# A new smoke run replaces the previous smoke directory. A full run is kept
# so an interrupted chain can resume from state.last.pickle.
if [[ "$MODE" == "smoke" ]]; then
    rm -rf "$output_dir"
fi
mkdir -p "$output_dir"
if [[ "$MODE" == "full" ]]; then
    cat > "$output_dir/params.txt" <<EOF
sample=${SAMPLE_ID}
chains=${CHAINS}
burnin=${BURNIN}
mcmc=${MCMC}
mh=${MH}
seeds=$(IFS=,; echo "${SEEDS[*]}")
chain_inclusion_factor=1.1
min_ssms=0.01
include_ssm_names=true
EOF
fi

# evolve.py needs a large stack. A failed ulimit is a warning, not a hard stop.
ulimit -s unlimited || echo "WARNING: could not set unlimited stack size; PhyloWGS may crash." >&2

###======================== Run PhyloWGS ================================
start=$(date -u +%Y-%m-%dT%H:%M:%SZ)
SECONDS=0
skipped=0
# summ.json is gzip-compressed despite the name. mutass.json is a zip.
if [[ "$MODE" == "full" && -s "$output_dir/summ.json" && -s "$output_dir/mutlist.json" && -s "$output_dir/mutass.json" ]]; then
    echo "PhyloWGS JSON already exists in $output_dir. Skipping the chains."
    skipped=1
else
    if [[ "$MODE" == "full" ]]; then
        mkdir -p "$output_dir/tmp"
        # Four chains. A chain that already has state.last.pickle resumes.
        TMPDIR="$output_dir/tmp" "$PYTHON2" "$PHYLOWGS_DIR/multievolve.py" \
            -n "$CHAINS" \
            -r "${SEEDS[@]}" \
            --ssms "$ssm" \
            --cnvs "$cnv" \
            -O "$output_dir" \
            -B "$BURNIN" \
            -s "$MCMC" \
            -i "$MH"
    else
        (
            cd "$output_dir"
            "$PYTHON2" "$PHYLOWGS_DIR/evolve.py" \
                -B "$BURNIN" \
                -s "$MCMC" \
                -i "$MH" \
                --random-seed "$SEED" \
                "$ssm" \
                "$cnv"
        )
    fi
    if [[ ! -s "$output_dir/trees.zip" ]]; then
        echo "trees.zip is missing. The chains may still be incomplete." >&2
        exit 1
    fi
    "$PYTHON2" "$PHYLOWGS_DIR/write_results.py" --include-ssm-names \
        "$SAMPLE_ID" \
        "$output_dir/trees.zip" \
        "$output_dir/summ.json" \
        "$output_dir/mutlist.json" \
        "$output_dir/mutass.json"
    printf '%s\n' "$MODE" > "$output_dir/mcmc_mode.txt"
    end=$(date -u +%Y-%m-%dT%H:%M:%SZ)
    if [[ "$MODE" == "full" ]]; then
        seed_text=$(IFS=,; echo "${SEEDS[*]}")
    else
        seed_text=$SEED
    fi
    cat > "$output_dir/runtime.txt" <<EOF
tool=phylowgs
mode=${MODE}
chains=${CHAINS}
burnin=${BURNIN}
mcmc=${MCMC}
mh=${MH}
seed=${SEED}
seeds=${seed_text}
start=${start}
end=${end}
elapsed_sec=${SECONDS}
EOF
fi

###======================== Extract clones ==============================
# clones.tsv columns: purity, n_clones, clone, ccf.
# The tree with the largest tree_densities value is used.
# Population 0 is the normal root and has no mutations, so it is not a clone.
# Each remaining population is one clone. Its ccf is cellular_prevalence.
# purity is the cellular prevalence of the cancer clone directly under the
# normal root. If that root has several children, their prevalences are added.
extract_clones() {
    python3 - "$output_dir/summ.json" "$ROOT/results/phylowgs/clones.tsv" <<'PY'
import csv
import gzip
import json
import sys

summary_path, out_path = sys.argv[1], sys.argv[2]
with gzip.open(summary_path, "rt") as handle:
    summary = json.load(handle)

# Prefer the posterior mode. Fall back to the first tree if densities are absent.
densities = summary.get("tree_densities") or {}
trees = summary.get("trees") or {}
if densities:
    tree_id = str(max(densities, key=lambda key: float(densities[key])))
elif trees:
    tree_id = str(next(iter(trees)))
else:
    sys.exit("PhyloWGS summary contains no trees.")

def lookup(mapping, key):
    # JSON object keys are strings. Tree structure lists may store child ids
    # as ints. Accept either form, including "1" versus 1.
    if key in mapping:
        return mapping[key]
    text = str(key)
    if text in mapping:
        return mapping[text]
    try:
        number = int(float(text))
    except ValueError:
        number = None
    if number is not None and number in mapping:
        return mapping[number]
    if number is not None and str(number) in mapping:
        return mapping[str(number)]
    return None

tree = lookup(trees, tree_id)
if tree is None:
    sys.exit(f"Selected PhyloWGS tree {tree_id} is not in summ.json.")
populations = tree["populations"]
structure = tree.get("structure") or {}

def prevalence(population_id):
    # One sample, so cellular_prevalence is a one-element list.
    population = lookup(populations, population_id)
    if population is None:
        sys.exit(f"Population {population_id} is missing from the selected PhyloWGS tree.")
    value = population.get("cellular_prevalence")
    if isinstance(value, list):
        if len(value) != 1:
            sys.exit("PhyloWGS cellular prevalence is not a single-sample value.")
        value = value[0]
    if value is None:
        sys.exit("PhyloWGS population has no cellular prevalence.")
    return float(value)

# Drop the normal root. Remaining ids are the cancer populations.
cancer_ids = sorted((str(population_id) for population_id in populations if str(population_id) != "0"), key=int)
with open(out_path, "w", newline="") as handle:
    writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
    writer.writerow(["purity", "n_clones", "clone", "ccf"])
    if not cancer_ids:
        writer.writerow(["NA", 0, "", ""])
    else:
        children = lookup(structure, "0") or []
        if not children:
            sys.exit("PhyloWGS tree has cancer populations but no clone under the normal root.")
        purity = sum(prevalence(child) for child in children)
        purity_text = f"{purity:.8g}"
        for population_id in cancer_ids:
            writer.writerow([purity_text, len(cancer_ids), population_id, f"{prevalence(population_id):.8g}"])
print(f"Wrote {out_path} from tree {tree_id}")
PY
}
extract_clones

###======================== Plot selected tree ==========================
# PhyloWGS has no fitted-SFS PDF. The official witness/ viewer is a browser
# app and is not used here. This figure is the same tree as clones.tsv:
# the one with the largest tree_densities value. Population 0 is the normal
# root. Labels are cellular prevalence. The title records MODE so a smoke
# chain is not read as a full-chain result. A plot failure must not discard
# the fit.
plot_tree() {
    local pdf_path="$output_dir/${SAMPLE_ID}_tree.pdf"
    if ! Rscript - "$output_dir/summ.json" "$pdf_path" "$MODE" "$SAMPLE_ID" <<'RS'
args <- commandArgs(trailingOnly = TRUE)
summary_path <- args[[1]]
pdf_path <- args[[2]]
mode <- args[[3]]
sample_id <- args[[4]]

suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(jsonlite))

summary <- jsonlite::fromJSON(gzfile(summary_path), simplifyVector = FALSE)
densities <- summary$tree_densities
trees <- summary$trees
if (is.null(trees) || !length(trees)) stop("PhyloWGS summary contains no trees.")
if (!is.null(densities) && length(densities)) {
    scores <- vapply(densities, function(value) as.numeric(value[[1]]), numeric(1))
    tree_id <- names(scores)[[which.max(scores)]]
} else {
    tree_id <- names(trees)[[1]]
}
tree <- trees[[tree_id]]
if (is.null(tree)) stop("Selected PhyloWGS tree ", tree_id, " is not in summ.json.")

as_ids <- function(value) {
    if (is.null(value) || !length(value)) return(character())
    as.character(unlist(value, use.names = FALSE))
}
children <- list()
structure <- tree$structure
if (!is.null(structure)) {
    for (parent in names(structure)) {
        children[[as.character(parent)]] <- as_ids(structure[[parent]])
    }
}
populations <- tree$populations
if (is.null(populations) || !length(populations)) stop("Selected PhyloWGS tree has no populations.")

prevalence <- function(population_id) {
    population <- populations[[as.character(population_id)]]
    if (is.null(population)) stop("Population ", population_id, " is missing from the selected PhyloWGS tree.")
    value <- population$cellular_prevalence
    if (is.list(value)) value <- unlist(value, use.names = FALSE)
    if (length(value) != 1L || !is.finite(value[[1]])) {
        stop("PhyloWGS cellular prevalence is not a single finite value.")
    }
    as.numeric(value[[1]])
}

# Leaves take the next x position. An internal node sits at the mean of its children.
# y is minus the depth, so the normal root is at the top.
state <- new.env(parent = emptyenv())
state$x <- 1
state$seen <- character()
state$rows <- list()
place <- function(node, depth) {
    node <- as.character(node)
    if (node %in% state$seen) stop("Cycle at PhyloWGS population ", node)
    state$seen <- c(state$seen, node)
    kids <- children[[node]]
    if (is.null(kids) || !length(kids)) {
        x <- state$x
        state$x <- state$x + 1
    } else {
        xs <- vapply(kids, function(kid) place(kid, depth + 1), numeric(1))
        x <- mean(xs)
    }
    state$rows[[length(state$rows) + 1L]] <- data.frame(
        id = node, x = x, y = -depth, stringsAsFactors = FALSE
    )
    x
}

ids <- as.character(names(populations))
child_ids <- unique(unlist(children, use.names = FALSE))
roots <- ids[!ids %in% child_ids]
if (!length(roots)) roots <- ids[[1]]
for (root in roots) place(root, 0)
for (node in setdiff(ids, state$seen)) place(node, 0)

nodes <- do.call(rbind, state$rows)
nodes$ccf <- vapply(nodes$id, prevalence, numeric(1))
nodes$label <- ifelse(
    nodes$id == "0",
    sprintf("normal\n%.3f", nodes$ccf),
    sprintf("pop %s\n%.3f", nodes$id, nodes$ccf)
)

edges <- data.frame(
    x = numeric(), y = numeric(), xend = numeric(), yend = numeric(),
    stringsAsFactors = FALSE
)
for (parent in names(children)) {
    parent_row <- nodes[nodes$id == parent, , drop = FALSE]
    if (!nrow(parent_row)) next
    for (kid in children[[parent]]) {
        child_row <- nodes[nodes$id == kid, , drop = FALSE]
        if (!nrow(child_row)) next
        edges <- rbind(edges, data.frame(
            x = parent_row$x[[1]], y = parent_row$y[[1]],
            xend = child_row$x[[1]], yend = child_row$y[[1]],
            stringsAsFactors = FALSE
        ))
    }
}

plot <- ggplot() +
    geom_segment(
        data = edges,
        aes(x = x, y = y, xend = xend, yend = yend),
        linewidth = 0.6
    ) +
    geom_label(data = nodes, aes(x = x, y = y, label = label), size = 4, linewidth = 0.3) +
    labs(
        title = paste0(sample_id, " PhyloWGS ", mode, " tree ", tree_id),
        subtitle = paste(
            "Posterior mode by tree density. Population 0 is the normal root.",
            "Labels are cellular prevalence. A smoke chain is not comparable to a full chain.",
            sep = "\n"
        )
    ) +
    theme_void(base_size = 14) +
    theme(
        plot.title = element_text(face = "bold"),
        plot.subtitle = element_text(color = "grey30"),
        plot.margin = margin(16, 16, 16, 16)
    )
n_leaves <- max(1, sum(vapply(ids, function(node) {
    kids <- children[[node]]
    is.null(kids) || !length(kids)
}, logical(1))))
depth <- max(0, -min(nodes$y))
ggsave(
    pdf_path, plot,
    width = max(9, 2.2 * n_leaves),
    height = max(4.5, 1.6 * (depth + 1))
)
cat("Wrote", pdf_path, "from tree", tree_id, "\n")
RS
    then
        echo "WARNING: PhyloWGS tree plot failed." >&2
    fi
}
plot_tree

if [[ "$skipped" == "1" ]]; then
    echo "PhyloWGS ${MODE} JSON was already present. Clones and the tree plot were refreshed."
else
    echo "PhyloWGS ${MODE} finished in ${SECONDS}s."
fi
if [[ "$MODE" == "full" ]]; then
    echo "Build the one-page view with: python3 \"$ROOT/src/make_phylowgs_html.py\" \"$output_dir\""
fi
