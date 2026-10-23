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
