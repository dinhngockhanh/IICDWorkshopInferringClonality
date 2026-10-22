#!/usr/bin/env Rscript

# Fit MOBSTER on diploid multiplicity-1 mutations, save the fit, then write clones.tsv.
# MOBSTER clusters VAF. It does not estimate purity and does not read purity.

###======= Parameter Config =======================================
sample_id <- "S8-noXY"
# Passed into mobster_fit(seed = ...). That function resets its own RNG,
# so set.seed alone would not control the fit.
seed <- 1001L
# FALSE keeps the fit on one core.
parallel <- FALSE
# Try 1, 2, and 3 Beta clusters. MOBSTER picks the best K.
k_min <- 1L
k_max <- 3L
# Fit both with and without a power-law tail, then keep the better model.
tail <- c(TRUE, FALSE)

###======= Paths ==================================================
file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
root <- normalizePath(file.path(dirname(sub("^--file=", "", file_arg[[1]])), ".."))
set.seed(seed)

# Diploid, multiplicity 1, VAF already in (0, 1): mutation_id, vaf, alt_count, depth.
input_path <- file.path(root, "data", "MOBSTER_mutation_table.csv")
output_dir <- file.path(root, "results", "mobster")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
if (!file.exists(input_path)) stop("Missing MOBSTER input: ", input_path)

if (!requireNamespace("mobster", quietly = TRUE)) {
    stop("R package mobster is not installed. Run this inside the workshop image.")
}
suppressPackageStartupMessages(library(mobster))

###======= Load input =============================================
mut <- read.csv(input_path, stringsAsFactors = FALSE)
# Drop any row whose VAF is missing or sits on the 0/1 boundary.
keep <- !is.na(mut$vaf) & mut$vaf > 0 & mut$vaf < 1
mut <- mut[keep, , drop = FALSE]
write.table(
    mut[, "mutation_id", drop = FALSE],
    file.path(output_dir, "used_mutation_ids.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE
)
cat("MOBSTER mutations:", nrow(mut), "\n")

###======= Run MOBSTER ============================================
started <- Sys.time()
fit <- mobster_fit(
    data.frame(mutation_id = mut$mutation_id, VAF = mut$vaf, stringsAsFactors = FALSE),
    K = seq.int(k_min, k_max),
    tail = tail,
    parallel = parallel,
    seed = seed,
    description = sample_id
)
finished <- Sys.time()

###======= Save results ===========================================
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
elapsed <- as.numeric(difftime(finished, started, units = "secs"))
writeLines(
    c(
        "tool=mobster",
        paste0("start=", format(started, "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")),
        paste0("end=", format(finished, "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")),
        paste0("elapsed_sec=", format(round(elapsed, 2), nsmall = 2))
    ),
    file.path(output_dir, "runtime.txt")
)

###======= Extract clones =========================================
# clones.tsv columns: purity, n_clones, clone, ccf.
# A cluster whose name contains "tail" is sequencing noise, not a clone.
# These mutations are diploid with multiplicity 1, so
#   ccf = 2 * mean(VAF) / sample_purity.
# sample_purity is the cellularity in data/sample.tsv. It is only used for
# that conversion. The purity column written below is the largest clone CCF.
extract_clones <- function() {
    best <- if (!is.null(fit$best)) fit$best else fit
    data_frame <- as.data.frame(best$data)
    if (!("cluster" %in% colnames(data_frame))) {
        stop("MOBSTER best$data has no cluster column: ", paste(colnames(data_frame), collapse = ", "))
    }
    if (!("VAF" %in% colnames(data_frame))) stop("MOBSTER best$data has no VAF column.")

    sample <- read.delim(file.path(root, "data", "sample.tsv"), stringsAsFactors = FALSE)
    sample_purity <- sample$purity[[1]]
    if (!is.finite(sample_purity) || sample_purity <= 0) stop("data/sample.tsv purity is not positive.")

    labels <- as.character(data_frame$cluster)
    means <- tapply(data_frame$VAF, labels, mean)
    clone_ids <- names(means)[!grepl("tail", names(means), ignore.case = TRUE)]
    # Order C1, C2, ... by the number in the label.
    clone_ids <- clone_ids[order(suppressWarnings(as.integer(sub("^[^0-9]*", "", clone_ids))), clone_ids)]
    ccfs <- 2 * as.numeric(means[clone_ids]) / sample_purity
    if (length(ccfs) > 0L && any(!is.finite(ccfs))) stop("Clone CCF is not finite.")
    purity <- if (length(ccfs) == 0L) NA_real_ else max(ccfs)
    if (length(clone_ids) == 0L) {
        clones <- data.frame(purity = "NA", n_clones = 0L, clone = "", ccf = "", stringsAsFactors = FALSE)
    } else {
        clones <- data.frame(
            purity = sprintf("%.8g", purity),
            n_clones = length(clone_ids),
            clone = clone_ids,
            ccf = sprintf("%.8g", ccfs),
            stringsAsFactors = FALSE
        )
    }
    path <- file.path(output_dir, "clones.tsv")
    write.table(clones, path, sep = "\t", quote = FALSE, row.names = FALSE)
    cat("Wrote", path, "\n")
}
extract_clones()

cat("MOBSTER clusters:\n")
print(fit$best$Clusters)
cat("MOBSTER finished.\n")
