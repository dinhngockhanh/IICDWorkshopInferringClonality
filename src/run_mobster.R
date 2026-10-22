#!/usr/bin/env Rscript

# MOBSTER on the diploid multiplicity-1 subset.
# K = 1:3 and tail = TRUE/FALSE are the package defaults used in the cluster timing run.
# parallel = FALSE keeps the fit on one core.
# mobster_fit resets the RNG from its seed argument, so the configured seed is passed in.

source_file <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- dirname(normalizePath(sub("^--file=", "", source_file[[1]]), mustWork = TRUE))
source(file.path(script_dir, "load_params.R"))

root <- workshop_root()
params <- load_params(root)
sample_id <- params$SAMPLE_ID
seed <- as.integer(params$SEED)
set.seed(seed)

input_path <- file.path(root, "data", "prepared", "inputs", "mobster", paste0(sample_id, ".vaf.csv"))
output_dir <- file.path(root, "results", "mobster")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
if (!file.exists(input_path)) stop("Missing MOBSTER input: ", input_path)

if (!requireNamespace("mobster", quietly = TRUE)) {
    stop("R package mobster is not installed. Run this inside the workshop image.")
}
suppressPackageStartupMessages(library(mobster))

mut <- read.csv(input_path, stringsAsFactors = FALSE)
keep <- !is.na(mut$vaf) & mut$vaf > 0 & mut$vaf < 1
mut <- mut[keep, , drop = FALSE]
write.table(
    mut[, "mutation_id", drop = FALSE],
    file.path(output_dir, "used_mutation_ids.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE
)
cat("MOBSTER mutations:", nrow(mut), "\n")

started <- Sys.time()
fit <- mobster_fit(
    data.frame(mutation_id = mut$mutation_id, VAF = mut$vaf, stringsAsFactors = FALSE),
    K = seq.int(as.integer(params$MOBSTER_K_MIN), as.integer(params$MOBSTER_K_MAX)),
    tail = c(TRUE, FALSE),
    parallel = as_flag(params$MOBSTER_PARALLEL),
    seed = seed,
    description = sample_id
)
finished <- Sys.time()

saveRDS(fit, file.path(output_dir, paste0(sample_id, "_mobster_fit.rds")))
pdf_path <- file.path(output_dir, paste0(sample_id, "_mobster_fit.pdf"))
tryCatch(
    {
        pdf(pdf_path, width = 10, height = 8)
        print(plot(fit$best))
        dev.off()
    },
    error = function(err) {
        if (dev.cur() != 1) dev.off()
        cat("WARNING: MOBSTER plot failed:", conditionMessage(err), "\n")
    }
)
write_runtime(file.path(output_dir, "runtime.txt"), "mobster", started, finished)
cat("MOBSTER clusters:\n")
print(fit$best$Clusters)
cat("MOBSTER finished.\n")
