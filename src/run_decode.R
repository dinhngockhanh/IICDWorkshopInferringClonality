###======================== Parameter Configuration =====================

sample_id <- "S8-noXY"
seed <- 1001L
# FALSE fits clusters only. TRUE would add a neutral-tail component.
neutral_tail <- FALSE
# FALSE keeps the fit on one core.
compute_parallel <- FALSE

###======================== Paths =======================================

file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
root <- normalizePath(file.path(dirname(sub("^--file=", "", file_arg[[1]])), ".."))
set.seed(seed)

# Diploid mutations only. Columns used below are Ref_count and Alt_count.
input_path <- file.path(root, "data", "DECODE_mutation_table.tsv")
output_dir <- file.path(root, "results", "decode")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
if (!file.exists(input_path)) stop("Missing DECODE input: ", input_path)

if (!requireNamespace("DECODE", quietly = TRUE)) {
    stop("R package DECODE is not installed. Run this inside the workshop image.")
}
suppressPackageStartupMessages(library(DECODE))

###======================== Load input ==================================

mutations <- read.delim(input_path, stringsAsFactors = FALSE)
cat("DECODE mutations:", nrow(mutations), "\n")

###======================== Run DECODE ==================================

started <- Sys.time()
decode_result <- DECODE(
    sample_id = sample_id,
    mutation_table = mutations[, c("Ref_count", "Alt_count")],
    neutral_tail = neutral_tail,
    compute_parallel = compute_parallel
)
finished <- Sys.time()

###======================== Save results ================================

saveRDS(decode_result, file.path(output_dir, paste0(sample_id, "_DECODE_result.rds")))
elapsed <- as.numeric(difftime(finished, started, units = "secs"))
writeLines(
    c(
        "tool=decode",
        paste0("start=", format(started, "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")),
        paste0("end=", format(finished, "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")),
        paste0("elapsed_sec=", format(round(elapsed, 2), nsmall = 2))
    ),
    file.path(output_dir, "runtime.txt")
)

###======================== Plot fits ===================================

model_selection_prefix <- file.path(output_dir, paste0(sample_id, "_DECODE_model_selection"))
tryCatch(
    {
        DECODE_plot_model_selection(
            decode_result,
            filename = model_selection_prefix,
            filetype = "pdf"
        )
        cat("Wrote", paste0(model_selection_prefix, ".pdf"), "\n")
    },
    error = function(err) {
        if (dev.cur() != 1) dev.off()
        cat("WARNING: DECODE model-selection plot failed:", conditionMessage(err), "\n")
    }
)

sfs_prefix <- file.path(output_dir, paste0(sample_id, "_DECODE_SFS"))
tryCatch(
    {
        DECODE_plot_SFS(
            decode_result,
            filename = sfs_prefix,
            filetype = "pdf"
        )
        cat("Wrote", paste0(sfs_prefix, ".pdf"), "\n")
    },
    error = function(err) {
        if (dev.cur() != 1) dev.off()
        cat("WARNING: DECODE SFS plot failed:", conditionMessage(err), "\n")
    }
)

###======================== Extract clones ==============================
# clones.tsv columns: purity, n_clones, clone, ccf.
# Cluster 1 is the clonal cluster. For a diploid mutation its VAF is purity / 2,
# so purity = 2 * mean(Cluster_1_freq).
# Clone i has ccf = 2 * mean(Cluster_i_freq).
# n_clones is best_N_clusters. A tail, if the fit selected one, is not a clone.

extract_clones <- function() {
    # Use the tail or no-tail fit that DECODE itself selected.
    if (isTRUE(decode_result$best_with_tail)) {
        fits <- decode_result$fits_with_tail$all_fits
    } else {
        fits <- decode_result$fits_without_tail$all_fits
    }
    n_clones <- as.integer(decode_result$best_N_clusters)
    key <- paste0(n_clones, "_clusters")
    fit_results <- fits[[key]]
    clonal <- fit_results$parameters[["Cluster_1_freq"]]
    clonal_vaf <- mean(as.numeric(unlist(clonal)))
    
    ccfs <- numeric(n_clones)
    for (cluster_id in seq_len(n_clones)) {
        column <- fit_results$parameters[[paste0("Cluster_", cluster_id, "_freq")]]
        cluster_vaf <- mean(as.numeric(unlist(column)))
        ccfs[[cluster_id]] <- 2 * cluster_vaf
    }
    clones <- data.frame(
        purity = sprintf("%.8g", 2 * clonal_vaf),
        n_clones = n_clones,
        clone = seq_len(n_clones),
        ccf = sprintf("%.8g", ccfs),
        stringsAsFactors = FALSE
    )
    path <- file.path(output_dir, "clones.tsv")
    write.table(clones, path, sep = "\t", quote = FALSE, row.names = FALSE)
    cat("Wrote", path, "\n")
}
extract_clones()

###======================== Assign mutations ============================
# The fit above used diploid mutations only. Assignment covers every mutation,
# including those that are not diploid. DECODE_mutation_assignment needs
# normalized_VAF: the VAF the mutation would have if it were diploid with
# multiplicity 1.

assign_mutations <- function() {
    if (isTRUE(decode_result$best_with_tail)) {
        fits <- decode_result$fits_with_tail$all_fits
    } else {
        fits <- decode_result$fits_without_tail$all_fits
    }
    key <- paste0(as.integer(decode_result$best_N_clusters), "_clusters")
    clonal <- fits[[key]]$parameters[["Cluster_1_freq"]]
    clonal_vaf <- mean(as.numeric(unlist(clonal)))
    purity <- 2 * clonal_vaf

    assignment_input <- file.path(root, "data", "DECODE_mutation_table_with_multiplicity.csv")
    if (!file.exists(assignment_input)) stop("Missing DECODE assignment input: ", assignment_input)

    mutations_all <- read.csv(assignment_input, stringsAsFactors = FALSE)
    depth <- mutations_all$Ref_count + mutations_all$Alt_count

    observed_vaf <- mutations_all$Alt_count / depth
    total_cn <- mutations_all$major_cn + mutations_all$minor_cn
    mutations_all$normalized_VAF <- observed_vaf * (total_cn * purity + 2 * (1 - purity)) / (2 * mutations_all$multiplicity)

    assigned <- DECODE_mutation_assignment(
        DECODE_result = decode_result,
        mutation_table = mutations_all
    )
    path <- file.path(output_dir, paste0(sample_id, "_mutation_assignment.csv"))
    write.table(assigned, path, sep = ",", quote = FALSE, row.names = FALSE)
    cat("Wrote", path, "\n")
}
assign_mutations()

cat(
    "DECODE finished. best_N_clusters=", decode_result$best_N_clusters,
    " best_with_tail=", decode_result$best_with_tail, "\n",
    sep = ""
)
