# Install the locked R library into the image site library.
options(repos = c(CRAN = "https://cloud.r-project.org"))
options(configure.args = list(rgl = "--disable-opengl"))
options(install.opts = list(
    rgl = "--configure-args='--disable-opengl'"
))
options(Ncpus = min(4L, max(1L, parallel::detectCores(logical = TRUE))))
options(timeout = 600)
Sys.setenv(RENV_CONFIG_SANDBOX_ENABLED = "FALSE")
Sys.setenv(RENV_CONFIG_INSTALL_REMOTES = "FALSE")
Sys.setenv(RENV_CONFIG_INSTALL_JOBS = "2")
Sys.setenv(RENV_CONFIG_PPM_ENABLED = "FALSE")
Sys.setenv(RENV_CONFIG_INSTALL_VERBOSE = "TRUE")
Sys.setenv(RGL_USE_NULL = "TRUE")

# renv 1.3.1 treats any archive that contains Meta/ as a binary package.
# easypar commits Meta/vignette.rds and has no Meta/package.rds, so restore
# would unpack the git tree and then fail while loading it.
patch_renv_package_type <- function() {
    ns <- asNamespace("renv")
    if (!exists("renv_package_type", envir = ns, inherits = FALSE)) {
        stop("renv_package_type is missing; this workaround needs renv 1.3.1")
    }
    orig <- get("renv_package_type", envir = ns)
    patched <- function(path, quiet = FALSE, default = "source") {
        type <- orig(path, quiet = TRUE, default = default)
        if (!identical(type, "binary")) {
            return(type)
        }
        files <- if (dir.exists(path)) {
            list.files(path, recursive = TRUE, all.files = TRUE)
        } else {
            tryCatch(utils::untar(path, list = TRUE), error = function(err) character())
        }
        if (!any(grepl("Meta/package\\.rds$", files))) {
            return("source")
        }
        type
    }
    unlockBinding("renv_package_type", ns)
    assign("renv_package_type", patched, envir = ns)
    lockBinding("renv_package_type", ns)
}

patch_renv_package_type()

library_dir <- "/usr/local/lib/R/site-library"
dir.create(library_dir, recursive = TRUE, showWarnings = FALSE)
renv::restore(
    lockfile = "/opt/workshop/renv.lock",
    library = library_dir,
    prompt = FALSE
)

stopifnot(requireNamespace("DECODE", quietly = TRUE))
stopifnot(requireNamespace("mobster", quietly = TRUE))
stopifnot(requireNamespace("sciClone", quietly = TRUE))
stopifnot(requireNamespace("easypar", quietly = TRUE))
message("R packages restored")
