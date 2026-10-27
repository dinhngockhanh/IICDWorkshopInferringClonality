###======================== Parameter Config ============================
sample_id <- "S8-noXY"
seed <- 1001L
# Mutations with depth below this value are dropped by sciClone during the fit.
minimum_depth <- 10L
# FALSE clusters on VAF only, not along the margins of the VAF distribution.
do_clustering_along_margins <- FALSE

###======================== Paths =======================================
file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
root <- normalizePath(file.path(dirname(sub("^--file=", "", file_arg[[1]])), ".."))
set.seed(seed)

output_dir <- file.path(root, "results", "sciclone")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
# VAF is stored as a percent from 0 to 100. Chromosomes have no "chr" prefix.
vaf_path <- file.path(root, "data", "sciClone_mutation_table.tsv")
# Copy-number calls: chr, start, stop, segment_mean = major_cn + minor_cn.
cn_path <- file.path(root, "data", "sciClone_cn.tsv")
# LOH intervals to exclude: chrom, start, end, with no header.
loh_path <- file.path(root, "data", "sciClone_exclude_loh.bed")
for (path in c(vaf_path, cn_path, loh_path)) {
    if (!file.exists(path)) stop("Missing sciClone input: ", path)
}

if (!requireNamespace("sciClone", quietly = TRUE)) {
    stop("R package sciClone is not installed. Run this inside the workshop image.")
}
suppressPackageStartupMessages(library(sciClone))

###======================== Load input ==================================
vafs <- read.delim(vaf_path, stringsAsFactors = FALSE)
copy_number <- read.delim(cn_path, stringsAsFactors = FALSE)
loh <- read.delim(loh_path, header = FALSE, stringsAsFactors = FALSE)
cat("sciClone mutations:", nrow(vafs), " segments:", nrow(copy_number), " LOH regions:", nrow(loh), "\n")

###======================== Run sciClone ================================
started <- Sys.time()
fit <- sciClone(
    vafs = vafs,
    copyNumberCalls = copy_number,
    sampleNames = sample_id,
    regionsToExclude = loh,
    minimumDepth = minimum_depth,
    doClustering = TRUE,
    doClusteringAlongMargins = do_clustering_along_margins
)
finished <- Sys.time()
if (is.null(fit)) {
    stop("sciClone returned NULL. No copy-number neutral sites were clustered.")
}

###======================== Save results ================================
saveRDS(fit, file.path(output_dir, paste0(sample_id, "_sciclone.rds")))
cluster_path <- file.path(output_dir, paste0(sample_id, "_clusters.tsv"))
writeClusterTable(fit, cluster_path)
plot_path <- file.path(output_dir, paste0(sample_id, "_clusters_1d.pdf"))
tryCatch(
    sc.plot1d(fit, plot_path),
    error = function(err) cat("WARNING: sciClone plot failed:", conditionMessage(err), "\n")
)

###======================== Plot fitted VAF =============================
# The raw curve integrates to about 50 on [0, 1]. Scale it to a density, then
# to expected counts in a 1% bin, and lay it over the observed histogram.
# Cluster 1 is the bottom layer. Mutations dropped for low depth are omitted.
# A plot failure must not discard the fit.
plot_sfs <- function() {
    pdf_path <- file.path(output_dir, paste0(sample_id, "_sfs.pdf"))
    tryCatch({
        if (!requireNamespace("ggplot2", quietly = TRUE)) {
            stop("R package ggplot2 is not installed.")
        }
        suppressPackageStartupMessages(library(ggplot2))
        clust <- fit@clust
        curves <- clust$individual.fits.y
        if (!identical(clust$cluster.method, "bmm") || is.null(curves) || is.null(clust$fit.x) || is.null(clust$fit.y)) {
            stop("sciClone fit has no BMM density curves.")
        }
        x <- as.numeric(clust$fit.x) / 100
        y_total_raw <- as.numeric(clust$fit.y)
        y_raw <- lapply(curves, as.numeric)
        if (length(x) == 0L || any(vapply(y_raw, length, integer(1)) != length(x)) || length(y_total_raw) != length(x)) {
            stop("sciClone BMM curve lengths do not match.")
        }
        # Grid is 999 points on (0, 1), so the step is 1 / (n + 1).
        dx <- 1 / (length(y_total_raw) + 1)
        integral <- sum(y_total_raw) * dx
        if (!is.finite(integral) || integral <= 1e-9) stop("sciClone BMM density does not integrate.")
        scale <- 1 / integral

        merged <- fit@vafs.merged
        ref_name <- paste0(sample_id, ".ref")
        var_name <- paste0(sample_id, ".var")
        depth <- merged[[ref_name]] + merged[[var_name]]
        keep <- is.finite(depth) & depth > 0 & !is.na(merged$cluster)
        vaf <- merged[[var_name]][keep] / depth[keep]
        vaf <- vaf[is.finite(vaf) & vaf > 0 & vaf <= 1]
        if (length(vaf) == 0L) stop("No sciClone mutations left to plot.")
        bin_width <- 0.01
        breaks <- seq(0, 1, by = bin_width)
        bin_idx <- pmin(as.integer(floor(vaf / bin_width)), length(breaks) - 2L) + 1L
        counts <- tabulate(bin_idx, nbins = length(breaks) - 1L)
        hist_data <- data.frame(
            xmin = breaks[-length(breaks)],
            xmax = breaks[-1],
            count = counts
        )
        to_counts <- function(y) y * scale * bin_width * length(vaf)
        y_total <- to_counts(y_total_raw)
        ribbons <- vector("list", length(y_raw))
        y_bottom <- rep(0, length(x))
        for (k in seq_along(y_raw)) {
            y_top <- y_bottom + to_counts(y_raw[[k]])
            ribbons[[k]] <- data.frame(
                vaf = x, ymin = y_bottom, ymax = y_top,
                cluster = factor(k, levels = seq_along(y_raw)),
                stringsAsFactors = FALSE
            )
            y_bottom <- y_top
        }
        ribbon <- do.call(rbind, ribbons)
        total <- data.frame(vaf = x, count = y_total)
        # Default ggplot2 discrete fills, in cluster-id order, at 75% opacity.
        # That is the same scale run_pyclonevi.sh uses: #F8766D, #00BA38, #619CFF.
        fills <- grDevices::adjustcolor(
            scales::hue_pal()(length(y_raw)),
            alpha.f = 0.75
        )

        plot <- ggplot() +
            geom_rect(
                data = hist_data,
                aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = count),
                fill = "#999999", color = NA
            ) +
            geom_ribbon(
                data = ribbon,
                aes(x = vaf, ymin = ymin, ymax = ymax, fill = cluster),
                color = NA
            ) +
            geom_line(data = total, aes(x = vaf, y = count), color = "grey20", linewidth = 0.6) +
            scale_fill_manual(values = fills, drop = FALSE) +
            scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
            scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
            labs(
                title = paste0(sample_id, " SciClone VAF fit"),
                subtitle = "Grey bars are observed counts in 1% bins. Colored areas are BMM components; the line is their sum.",
                x = "Variant allele frequency",
                y = "Mutations",
                fill = "Cluster"
            ) +
            theme_bw()
        ggsave(pdf_path, plot, width = 10, height = 6)
        cat("Wrote", pdf_path, "\n")
    }, error = function(err) {
        cat("WARNING: sciClone VAF plot failed:", conditionMessage(err), "\n")
    })
}
plot_sfs()
elapsed <- as.numeric(difftime(finished, started, units = "secs"))
writeLines(
    c(
        "tool=sciclone",
        paste0("start=", format(started, "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")),
        paste0("end=", format(finished, "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")),
        paste0("elapsed_sec=", format(round(elapsed, 2), nsmall = 2))
    ),
    file.path(output_dir, "runtime.txt")
)

###======================== Extract clones ==============================
# clones.tsv columns: purity, n_clones, clone, ccf.
# fit@clust$cluster.means is samples x clusters, on the VAF fraction scale (0 to 1),
# even though the input file stores VAF as a percent.
# This workshop has one sample, so the matrix has one row and one column per clone.
# For diploid multiplicity-1 mutations,
#   ccf = 2 * cluster_mean.
# The purity column is the largest clone CCF. Mutations dropped for low depth are not listed.
extract_clones <- function() {
    mat <- as.matrix(fit@clust$cluster.means)
    if (nrow(mat) == 1L) {
        means <- as.numeric(mat[1, ])
    } else if (ncol(mat) == 1L) {
        means <- as.numeric(mat[, 1])
    } else {
        means <- as.numeric(mat[1, ])
    }
    ccfs <- 2 * means
    if (length(ccfs) > 0L && any(!is.finite(ccfs))) stop("Clone CCF is not finite.")
    purity <- if (length(ccfs) == 0L) NA_real_ else max(ccfs)
    if (length(means) == 0L) {
        clones <- data.frame(purity = "NA", n_clones = 0L, clone = "", ccf = "", stringsAsFactors = FALSE)
    } else {
        clones <- data.frame(
            purity = sprintf("%.8g", purity),
            n_clones = length(means),
            clone = seq_along(means),
            ccf = sprintf("%.8g", ccfs),
            stringsAsFactors = FALSE
        )
    }
    path <- file.path(output_dir, "clones.tsv")
    write.table(clones, path, sep = "\t", quote = FALSE, row.names = FALSE)
    cat("Wrote", path, "\n")
}
extract_clones()

cat("sciClone cluster means:\n")
print(fit@clust$cluster.means)
cat("sciClone finished.\n")
