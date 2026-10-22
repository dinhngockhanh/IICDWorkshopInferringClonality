# Read config/params.toml through src/read_params.py.
# Rscript entry points live in src/, so the workshop root is one directory up.

workshop_root <- function() {
    file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
    if (length(file_arg) >= 1L) {
        script_dir <- dirname(normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = TRUE))
        return(normalizePath(file.path(script_dir, "..")))
    }
    normalizePath(getwd())
}

load_params <- function(root) {
    python <- Sys.which("python3")
    if (!nzchar(python)) {
        stop("python3 is required to read config/params.toml")
    }
    lines <- system2(
        python,
        c(file.path(root, "src", "read_params.py"), "--format", "shell"),
        stdout = TRUE
    )
    status <- attr(lines, "status")
    if (!is.null(status) && status != 0) {
        stop("src/read_params.py failed")
    }
    params <- list()
    for (line in lines) {
        parts <- strsplit(line, "=", fixed = TRUE)[[1]]
        key <- parts[[1]]
        raw <- paste(parts[-1], collapse = "=")
        raw <- gsub("^'|'$", "", raw)
        params[[key]] <- raw
    }
    params
}

as_flag <- function(value) {
    toupper(value) %in% c("TRUE", "T", "1")
}

write_runtime <- function(path, tool, start_time, end_time) {
    elapsed <- as.numeric(difftime(end_time, start_time, units = "secs"))
    writeLines(
        c(
            paste0("tool=", tool),
            paste0("start=", format(start_time, "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")),
            paste0("end=", format(end_time, "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")),
            paste0("elapsed_sec=", format(round(elapsed, 2), nsmall = 2))
        ),
        path
    )
}
