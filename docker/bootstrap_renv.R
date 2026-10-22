# Build renv.lock inside rocker/r-ver:4.5.1.
# The project library stays in /tmp so the mounted workshop repo only receives the lock.
# Run from docker/Dockerfile.renv-bootstrap. Students do not need to re-run this.

options(repos = c(CRAN = "https://cloud.r-project.org"))
options(configure.args = list(rgl = "--disable-opengl"))
options(install.opts = list(rgl = "--configure-args='--disable-opengl'"))
options(Ncpus = min(4L, max(1L, parallel::detectCores(logical = TRUE))))
options(timeout = 600)
Sys.setenv(RENV_CONFIG_SANDBOX_ENABLED = "FALSE")
Sys.setenv(RENV_CONFIG_PPM_ENABLED = "FALSE")
# One install job avoids a dependent starting before its dependencies are visible.
Sys.setenv(RENV_CONFIG_INSTALL_JOBS = "1")
Sys.setenv(RENV_CONFIG_INSTALL_STAGED = "FALSE")
Sys.setenv(RENV_CONFIG_INSTALL_TRANSACTIONAL = "FALSE")
# Install GitHub remotes ourselves. Leaving this on also installs MOBSTER's
# Suggests remote CNAqc, which the fit does not load.
Sys.setenv(RENV_CONFIG_INSTALL_REMOTES = "FALSE")
Sys.setenv(RENV_CONFIG_INSTALL_VERBOSE = "TRUE")

if (!requireNamespace("renv", quietly = TRUE)) {
    install.packages("renv")
}

project_dir <- "/tmp/renvproj"
dir.create(project_dir, recursive = TRUE, showWarnings = FALSE)
renv::init(
    project = project_dir,
    bare = TRUE,
    bioconductor = TRUE,
    force = TRUE
)

# Download GitHub archives directly. renv treats codeload URLs as GitHub remotes
# and then calls the rate-limited API. A local tarball skips that.
# dependencies = FALSE keeps genome/bmm from being replaced by the CRAN package
# of the same name, which imports brms.
strong <- c("Depends", "Imports", "LinkingTo")

stamp_github <- function(pkg, user, repo, sha) {
    lib <- renv::paths$library(project = project_dir)
    desc_path <- file.path(lib, pkg, "DESCRIPTION")
    if (!file.exists(desc_path)) {
        stop("Missing DESCRIPTION for ", pkg)
    }
    desc <- read.dcf(desc_path)
    keep <- !startsWith(colnames(desc), "Remote") & !startsWith(colnames(desc), "Github")
    desc <- desc[, keep, drop = FALSE]
    extra <- matrix(
        c("github", "api.github.com", user, repo, sha, sha,
          "api.github.com", user, repo, sha, sha),
        nrow = 1,
        dimnames = list(NULL, c(
            "RemoteType", "RemoteHost", "RemoteUsername", "RemoteRepo", "RemoteRef", "RemoteSha",
            "GithubHost", "GithubUsername", "GithubRepo", "GithubRef", "GithubSHA1"
        ))
    )
    write.dcf(cbind(desc, extra), desc_path)
}

install_github_tarball <- function(pkg, user, repo, sha) {
    url <- sprintf("https://codeload.github.com/%s/%s/tar.gz/%s", user, repo, sha)
    dest <- file.path("/tmp", paste0(pkg, "-", sha, ".tar.gz"))
    message("Downloading ", url)
    status <- download.file(url, dest, mode = "wb", quiet = FALSE)
    if (status != 0) {
        stop("Download failed: ", url)
    }
    # An unused field name: NULL and FALSE both fall back to the default
    # dependency set, which resolves CRAN's bmm and pulls in brms.
    pkg_root <- tempfile("pkg-")
    dir.create(pkg_root)
    untar(dest, exdir = pkg_root)
    pkg_dir <- list.dirs(pkg_root, full.names = TRUE, recursive = FALSE)[[1]]
    unlink(file.path(pkg_dir, c("vignettes", "Meta", "docs")), recursive = TRUE)
    renv::install(pkg_dir, prompt = FALSE, dependencies = "Config/Needs/none", verbose = TRUE)
    stamp_github(pkg, user, repo, sha)
}

use_bioconductor_repos <- function() {
    repos <- c(
        CRAN = "https://cloud.r-project.org",
        BioCsoft = "https://bioconductor.org/packages/3.22/bioc",
        BioCann = "https://bioconductor.org/packages/3.22/data/annotation",
        BioCexp = "https://bioconductor.org/packages/3.22/data/experiment"
    )
    options(repos = repos)
    repos
}

cran_and_bioc <- c(
    "bioc::limma",
    "bioc::IRanges",
    "plotrix",
    "TeachingDemos",
    "RColorBrewer",
    "remotes"
)
message("Installing CRAN and Bioconductor dependencies")
renv::install(cran_and_bioc, prompt = FALSE)

message("Installing MKmisc (CRAN, then archive if it has been removed)")
mkmisc_status <- tryCatch(
    {
        renv::install("MKmisc", prompt = FALSE)
        TRUE
    },
    error = function(err) {
        message("CRAN MKmisc failed: ", conditionMessage(err))
        FALSE
    }
)
if (!isTRUE(mkmisc_status)) {
    renv::install(
        "https://cran.r-project.org/src/contrib/Archive/MKmisc/MKmisc_1.9.tar.gz",
        prompt = FALSE
    )
}

message("Installing rgl without OpenGL")
renv::install("rgl", prompt = FALSE)

message("Installing sciClone CRAN dependencies without resolving the bmm name")
renv::install(
    c("dplyr", "ggplot2", "tibble", "tidyr", "stringr"),
    prompt = FALSE,
    dependencies = strong
)
install_github_tarball("bmm", "genome", "bmm", "ae560f70b0d40e5bf138921bcfa9e795d7ea6824")
install_github_tarball("sciClone", "genome", "sciclone", "1b3d055b0fe1bd523189a29926791e912577a7e0")
bmm_imports <- packageDescription("bmm", lib.loc = renv::paths$library(project = project_dir))$Imports
if (!is.null(bmm_imports) && grepl("brms", bmm_imports)) {
    stop("Installed bmm imports brms. That is the CRAN package, not genome/bmm.")
}

message("Installing Bioconductor packages required by dndscv")
use_bioconductor_repos()
renv::install(
    c(
        "bioc::GenomeInfoDbData",
        "bioc::GenomeInfoDb",
        "bioc::XVector",
        "bioc::Rhtslib",
        "bioc::GenomicRanges",
        "bioc::Biostrings",
        "bioc::Rsamtools",
        "bioc::rtracklayer",
        "bioc::BiocParallel"
    ),
    prompt = FALSE,
    dependencies = c("Depends", "Imports", "LinkingTo")
)

message("Installing CRAN dependencies of MOBSTER, ctree, easypar, and DECODE")
renv::install(
    c(
        "cli", "clisymbols", "cowplot", "crayon", "data.table", "doParallel", "dplyr",
        "drf", "entropy", "foreach", "ggpubr", "ggplot2", "ggraph", "ggrepel",
        "ggthemes", "gridExtra", "igraph", "magrittr", "markdown", "MASS", "matrixcalc",
        "pbapply", "plotrix", "plyr", "poilog", "prettydoc", "progress", "RColorBrewer",
        "readr", "reshape2", "sads", "seqinr", "shadowtext", "stringr", "tibble",
        "tidygraph", "tidyverse", "tidyr", "truncnorm"
    ),
    prompt = FALSE,
    dependencies = strong
)
install_github_tarball("pio", "caravagn", "pio", "2cfe575ba8c41905fc43540fd489e79b1bac0d5e")
install_github_tarball("easypar", "caravagnalab", "easypar", "226e2162ee6eeac44175ab6fb33baf5d5410b6ad")
install_github_tarball("dndscv", "im3sanger", "dndscv", "43c5e2f1beaaf3ac26df9c046b281ceb0edbfff7")
install_github_tarball("ctree", "caravagnalab", "ctree", "d012b580882a26364c7fe804f25bc1f58edfa4c1")
install_github_tarball("mobster", "caravagnalab", "mobster", "85c898f087b46e15f79144bf57bdc019678e2481")
install_github_tarball("DECODE", "dinhngockhanh", "DECODE", "654cbcf8aded8ab6b0322e5762b68c49031fa61d")

lock_path <- "/out/renv.lock"
renv::snapshot(
    project = project_dir,
    lockfile = lock_path,
    type = "all",
    prompt = FALSE
)
message("Wrote ", lock_path)
lib_dir <- renv::paths$library(project = project_dir)
.libPaths(c(lib_dir, .libPaths()))
for (pkg in c("DECODE", "mobster", "sciClone", "bmm")) {
    message(pkg, " ", as.character(packageVersion(pkg)))
}
