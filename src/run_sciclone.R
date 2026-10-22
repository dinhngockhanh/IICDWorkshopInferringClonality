#!/usr/bin/env Rscript

# sciClone on the same diploid subset as MOBSTER.
# Copy-number calls and the LOH bed are passed, as in the single-sample cluster script.
# minimumDepth = 10 removes shallow mutations inside sciClone rather than in the input file.

source_file <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- dirname(normalizePath(sub("^--file=", "", source_file[[1]]), mustWork = TRUE))
source(file.path(script_dir, "load_params.R"))

root <- workshop_root()
params <- load_params(root)
sample_id <- params$SAMPLE_ID
set.seed(as.integer(params$SEED))

input_dir <- file.path(root, "data", "prepared", "inputs", "sciclone")
output_dir <- file.path(root, "results", "sciclone")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
vaf_path <- file.path(input_dir, "vafs.tsv")
cn_path <- file.path(input_dir, "cn.tsv")
loh_path <- file.path(input_dir, "exclude.loh.bed")
for (path in c(vaf_path, cn_path, loh_path)) {
    if (!file.exists(path)) stop("Missing sciClone input: ", path)
}

if (!requireNamespace("sciClone", quietly = TRUE)) {
    stop("R package sciClone is not installed. Run this inside the workshop image.")
}
suppressPackageStartupMessages(library(sciClone))

vafs <- read.delim(vaf_path, stringsAsFactors = FALSE)
copy_number <- read.delim(cn_path, stringsAsFactors = FALSE)
loh <- read.delim(loh_path, header = FALSE, stringsAsFactors = FALSE)
cat("sciClone mutations:", nrow(vafs), " segments:", nrow(copy_number), " LOH regions:", nrow(loh), "\n")

started <- Sys.time()
fit <- sciClone(
    vafs = vafs,
    copyNumberCalls = copy_number,
    sampleNames = sample_id,
    regionsToExclude = loh,
    minimumDepth = as.integer(params$SCICLONE_MINIMUM_DEPTH),
    doClustering = TRUE,
    doClusteringAlongMargins = as_flag(params$SCICLONE_ALONG_MARGINS)
)
finished <- Sys.time()
if (is.null(fit)) {
    stop("sciClone returned NULL. No copy-number neutral sites were clustered.")
}

saveRDS(fit, file.path(output_dir, paste0(sample_id, "_sciclone.rds")))
cluster_path <- file.path(output_dir, paste0(sample_id, "_clusters.tsv"))
writeClusterTable(fit, cluster_path)
plot_path <- file.path(output_dir, paste0(sample_id, "_clusters_1d.pdf"))
tryCatch(
    sc.plot1d(fit, plot_path),
    error = function(err) cat("WARNING: sciClone plot failed:", conditionMessage(err), "\n")
)
write_runtime(file.path(output_dir, "runtime.txt"), "sciclone", started, finished)
cat("sciClone cluster means:\n")
print(fit@clust$cluster.means)
cat("sciClone finished.\n")
