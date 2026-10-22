#!/usr/bin/env Rscript

# DECODE on diploid S8 mutations.
# Matches the DREAM pipeline call: neutral_tail = FALSE, compute_parallel = FALSE,
# and the package default ABCSMCDRF_nParticles = rep(1000, 10).
# Multiplicity is already in the mutation table, so this script does not call GRITIC.

source_file <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- dirname(normalizePath(sub("^--file=", "", source_file[[1]]), mustWork = TRUE))
source(file.path(script_dir, "load_params.R"))

root <- workshop_root()
params <- load_params(root)
sample_id <- params$SAMPLE_ID
set.seed(as.integer(params$SEED))

input_path <- file.path(root, "data", "prepared", "inputs", "decode", "mutations.tsv")
output_dir <- file.path(root, "results", "decode")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
if (!file.exists(input_path)) stop("Missing DECODE input: ", input_path)

if (!requireNamespace("DECODE", quietly = TRUE)) {
    stop("R package DECODE is not installed. Run this inside the workshop image.")
}
suppressPackageStartupMessages(library(DECODE))

mutations <- read.delim(input_path, stringsAsFactors = FALSE)
cat("DECODE mutations:", nrow(mutations), "\n")

started <- Sys.time()
decode_result <- DECODE(
    sample_id = sample_id,
    mutation_table = mutations[, c("Ref_count", "Alt_count")],
    neutral_tail = as_flag(params$DECODE_NEUTRAL_TAIL),
    compute_parallel = as_flag(params$DECODE_COMPUTE_PARALLEL)
)
finished <- Sys.time()

saveRDS(decode_result, file.path(output_dir, paste0(sample_id, "_DECODE_result.rds")))
write_runtime(file.path(output_dir, "runtime.txt"), "decode", started, finished)
cat(
    "DECODE finished. best_N_clusters=", decode_result$best_N_clusters,
    " best_with_tail=", decode_result$best_with_tail, "\n",
    sep = ""
)
