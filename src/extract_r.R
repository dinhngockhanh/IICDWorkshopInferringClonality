#!/usr/bin/env Rscript

# Write per-mutation rows for the R tools.
# Columns: mutation_id, cluster, ccf, note
# extract_results.py joins these onto the full mutation table.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: Rscript extract_r.R <decode|mobster|sciclone>")
tool <- args[[1]]

source_file <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- dirname(normalizePath(sub("^--file=", "", source_file[[1]]), mustWork = TRUE))
source(file.path(script_dir, "load_params.R"))
root <- workshop_root()
params <- load_params(root)
sample_id <- params$SAMPLE_ID

read_purity <- function() {
    sample <- read.delim(file.path(root, "data", "prepared", "sample.tsv"), stringsAsFactors = FALSE)
    sample$purity[[1]]
}

write_components <- function(path, frame) {
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    write.table(frame, path, sep = "\t", quote = FALSE, row.names = FALSE)
}

format_ccf <- function(values) {
    ifelse(is.na(values), "", sprintf("%.8g", values))
}

extract_decode <- function() {
    fit_path <- file.path(root, "results", "decode", paste0(sample_id, "_DECODE_result.rds"))
    if (!file.exists(fit_path)) stop("Missing DECODE result: ", fit_path)
    if (!requireNamespace("DECODE", quietly = TRUE)) {
        stop("R package DECODE is not installed.")
    }
    suppressPackageStartupMessages(library(DECODE))
    decode_result <- readRDS(fit_path)

    if (isTRUE(decode_result$best_with_tail)) {
        fits <- decode_result$fits_with_tail$all_fits
    } else {
        fits <- decode_result$fits_without_tail$all_fits
    }
    key <- paste0(decode_result$best_N_clusters, "_clusters")
    if (!(key %in% names(fits))) {
        stop("DECODE fit ", key, " not found. Available: ", paste(names(fits), collapse = ", "))
    }
    fit_results <- fits[[key]]

    clonal_vaf <- mean(as.numeric(unlist(fit_results$parameters[["Cluster_1_freq"]])))
    purity <- 2 * clonal_vaf
    writeLines(sprintf("%.8g", purity), file.path(root, "results", "decode", "purity.txt"))

    sfs_components <- list()
    if (isTRUE(decode_result$best_with_tail) && "SFS_validation_tail" %in% names(fit_results)) {
        sfs_components[["tail"]] <- colMeans(fit_results[["SFS_validation_tail"]], na.rm = TRUE)
    }
    cluster_names <- grep("^SFS_validation_cluster_", names(fit_results), value = TRUE)
    for (cluster_name in cluster_names) {
        component <- sub("^SFS_validation_cluster_", "cluster_", cluster_name)
        sfs_components[[component]] <- colMeans(fit_results[[cluster_name]], na.rm = TRUE)
    }
    if (length(sfs_components) == 0L) stop("DECODE fit has no SFS validation components.")

    sfs <- as.data.frame(do.call(cbind, sfs_components), check.names = FALSE)
    colnames(sfs) <- names(sfs_components)
    row_sums <- rowSums(sfs, na.rm = TRUE)
    sfs_norm <- sfs
    for (i in seq_len(nrow(sfs))) {
        if (!is.na(row_sums[[i]]) && row_sums[[i]] > 0) {
            sfs_norm[i, ] <- sfs[i, ] / row_sums[[i]]
        } else {
            sfs_norm[i, ] <- 0
        }
    }
    sfs_norm[is.na(sfs_norm)] <- 0
    sfs_norm[!is.finite(as.matrix(sfs_norm))] <- 0
    write.table(
        sfs_norm,
        file.path(root, "results", "decode", "sfs_proportions.tsv"),
        sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA
    )

    mutations <- read.delim(file.path(root, "data", "prepared", "mutations.tsv"), stringsAsFactors = FALSE)
    multiplicity <- as.numeric(mutations$multiplicity)
    total_cn <- as.numeric(mutations$major_cn) + as.numeric(mutations$minor_cn)
    denom <- (multiplicity * purity) / (total_cn * purity + 2 * (1 - purity))
    adjust <- (purity / 2) / denom
    adjust[is.na(denom) | denom == 0 | total_cn == 0] <- NA
    adjusted_vaf <- (mutations$alt_count / mutations$depth) * adjust

    n_bins <- nrow(sfs_norm)
    component_names <- colnames(sfs_norm)
    boundaries <- seq(0, 1, length.out = n_bins + 1)
    assignment <- rep(NA_integer_, nrow(mutations))
    for (i in seq_len(nrow(mutations))) {
        value <- adjusted_vaf[[i]]
        if (is.na(value) || value < 0 || value > 1) next
        bin_idx <- findInterval(value, boundaries, rightmost.closed = TRUE)
        if (bin_idx == 0) bin_idx <- 1
        if (bin_idx > n_bins) bin_idx <- n_bins
        proportions <- as.numeric(sfs_norm[bin_idx, ])
        cluster_idx <- which(component_names != "tail")
        if (length(cluster_idx) == 0L) next
        chosen <- cluster_idx[[which.max(proportions[cluster_idx])]]
        assignment[[i]] <- as.integer(sub("^cluster_", "", component_names[[chosen]]))
    }

    max_cluster <- if (all(is.na(assignment))) 0L else max(assignment, na.rm = TRUE)
    mean_vaf <- rep(NA_real_, max_cluster)
    for (cluster_id in seq_len(max_cluster)) {
        column <- fit_results$parameters[[paste0("Cluster_", cluster_id, "_freq")]]
        if (!is.null(column)) mean_vaf[[cluster_id]] <- mean(as.numeric(unlist(column)))
    }
    ccf <- rep(NA_real_, nrow(mutations))
    known <- !is.na(assignment) & assignment >= 1 & assignment <= length(mean_vaf)
    ccf[known] <- mean_vaf[assignment[known]] / clonal_vaf

    if (max_cluster > 0L) {
        summary_frame <- data.frame(
            cluster = seq_len(max_cluster),
            mean_vaf = mean_vaf,
            ccf = mean_vaf / clonal_vaf,
            stringsAsFactors = FALSE
        )
        write.table(
            summary_frame,
            file.path(root, "results", "decode", "cluster_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE
        )
    }

    keep <- !is.na(assignment)
    write_components(
        file.path(root, "results", "decode", "components.tsv"),
        data.frame(
            mutation_id = mutations$mutation_id[keep],
            cluster = as.character(assignment[keep]),
            ccf = format_ccf(ccf[keep]),
            note = "",
            stringsAsFactors = FALSE
        )
    )
}

extract_mobster <- function() {
    fit_path <- file.path(root, "results", "mobster", paste0(sample_id, "_mobster_fit.rds"))
    if (!file.exists(fit_path)) stop("Missing MOBSTER result: ", fit_path)
    if (!requireNamespace("mobster", quietly = TRUE)) stop("R package mobster is not installed.")
    suppressPackageStartupMessages(library(mobster))
    fit <- readRDS(fit_path)
    best <- if (!is.null(fit$best)) fit$best else fit
    data_frame <- as.data.frame(best$data)
    if (!("cluster" %in% colnames(data_frame))) {
        stop("MOBSTER best$data has no cluster column: ", paste(colnames(data_frame), collapse = ", "))
    }
    if (!("VAF" %in% colnames(data_frame))) stop("MOBSTER best$data has no VAF column.")
    if (!("mutation_id" %in% colnames(data_frame))) {
        id_path <- file.path(root, "results", "mobster", "used_mutation_ids.tsv")
        ids <- read.delim(id_path, stringsAsFactors = FALSE)
        if (nrow(ids) != nrow(data_frame)) {
            stop("MOBSTER output dropped mutation_id and the row count does not match used_mutation_ids.tsv.")
        }
        data_frame$mutation_id <- ids$mutation_id
    }

    purity <- read_purity()
    is_tail <- grepl("tail", as.character(data_frame$cluster), ignore.case = TRUE)
    means <- tapply(data_frame$VAF, as.character(data_frame$cluster), mean)
    ccf <- ifelse(
        is_tail,
        NA_real_,
        2 * as.numeric(means[as.character(data_frame$cluster)]) / purity
    )
    write_components(
        file.path(root, "results", "mobster", "components.tsv"),
        data.frame(
            mutation_id = data_frame$mutation_id,
            cluster = as.character(data_frame$cluster),
            ccf = format_ccf(ccf),
            note = ifelse(is_tail, "tail", ""),
            stringsAsFactors = FALSE
        )
    )
}

cluster_mean_vector <- function(means) {
    # sciClone stores cluster.means as samples x clusters, on the VAF fraction scale.
    mat <- as.matrix(means)
    if (nrow(mat) == 1L) {
        values <- as.numeric(mat[1, ])
    } else if (ncol(mat) == 1L) {
        values <- as.numeric(mat[, 1])
    } else {
        values <- as.numeric(mat[1, ])
    }
    names(values) <- as.character(seq_along(values))
    values
}

pick_col <- function(candidates, available) {
    found <- intersect(candidates, available)
    if (length(found) == 0L) NA_character_ else found[[1]]
}

extract_sciclone <- function() {
    fit_path <- file.path(root, "results", "sciclone", paste0(sample_id, "_sciclone.rds"))
    table_path <- file.path(root, "results", "sciclone", paste0(sample_id, "_clusters.tsv"))
    if (!file.exists(fit_path)) stop("Missing sciClone result: ", fit_path)
    if (!file.exists(table_path)) stop("Missing sciClone cluster table: ", table_path)
    if (!requireNamespace("sciClone", quietly = TRUE)) stop("R package sciClone is not installed.")
    suppressPackageStartupMessages(library(sciClone))
    fit <- readRDS(fit_path)
    clusters <- read.delim(table_path, stringsAsFactors = FALSE, check.names = FALSE)
    names(clusters) <- tolower(names(clusters))
    chrom_col <- pick_col(c("chr", "chromosome", "chrom"), names(clusters))
    pos_col <- pick_col(c("st", "pos", "position", "start"), names(clusters))
    cluster_col <- pick_col(c("cluster", "clust"), names(clusters))
    if (is.na(chrom_col) || is.na(pos_col) || is.na(cluster_col)) {
        stop("Unrecognized sciClone cluster columns: ", paste(names(clusters), collapse = ", "))
    }
    clusters <- clusters[!is.na(clusters[[cluster_col]]), , drop = FALSE]

    purity <- read_purity()
    means <- cluster_mean_vector(fit@clust$cluster.means)
    cluster_ids <- as.character(clusters[[cluster_col]])
    if (!all(cluster_ids %in% names(means))) {
        stop(
            "sciClone cluster ids are not in cluster.means. ids=",
            paste(unique(cluster_ids), collapse = ","),
            " means=",
            paste(names(means), collapse = ",")
        )
    }
    ccf <- 2 * means[cluster_ids] / purity
    mutation_id <- paste(clusters[[chrom_col]], clusters[[pos_col]], sep = ":")
    write_components(
        file.path(root, "results", "sciclone", "components.tsv"),
        data.frame(
            mutation_id = mutation_id,
            cluster = cluster_ids,
            ccf = format_ccf(ccf),
            note = "",
            stringsAsFactors = FALSE
        )
    )
}

if (tool == "decode") {
    extract_decode()
} else if (tool == "mobster") {
    extract_mobster()
} else if (tool == "sciclone") {
    extract_sciclone()
} else {
    stop("Unknown tool: ", tool)
}
cat("Wrote", tool, "components.\n")
