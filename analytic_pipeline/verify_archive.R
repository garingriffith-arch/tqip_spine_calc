# Verify public SPINE-TRACT source-archive integrity without ACS data.
# Run from the repository root:
#   Rscript analytic_pipeline/verify_archive.R

manifest_path <- file.path("analytic_pipeline", "PIPELINE_MANIFEST.csv")
archive_path <- file.path("analytic_pipeline", "SOURCE_ARCHIVE_MANIFEST.csv")

if (!file.exists(manifest_path) || !file.exists(archive_path)) {
  stop("Missing pipeline/source archive manifest.", call. = FALSE)
}

pipeline <- read.csv(manifest_path, stringsAsFactors = FALSE,
                     check.names = FALSE)
archive <- read.csv(archive_path, stringsAsFactors = FALSE,
                    check.names = FALSE)

expected <- pipeline[pipeline$step != "PENDING", , drop = FALSE]
if (nrow(expected) != 20L || nrow(archive) != 20L) {
  stop("Expected 20 archived implementation stages.", call. = FALSE)
}

required <- c("step", "entrypoint", "implementation", "git_blob_sha1")
if (!all(required %in% names(archive))) {
  stop("Source archive manifest has missing required columns.", call. = FALSE)
}

if (anyDuplicated(archive$step) || anyDuplicated(archive$implementation)) {
  stop("Duplicate source archive entry.", call. = FALSE)
}

idx <- match(expected$step, archive$step)
if (anyNA(idx)) {
  stop("Pipeline and source manifests have different stage identifiers.",
       call. = FALSE)
}

a <- archive[idx, , drop = FALSE]
if (!identical(expected$entrypoint, a$entrypoint) ||
    !identical(expected$implementation, a$implementation)) {
  stop("Pipeline and source manifests disagree on file paths.",
       call. = FALSE)
}

if (!all(grepl("^[0-9a-f]{40}$", a$git_blob_sha1))) {
  stop("Invalid Git blob identifier in source archive manifest.",
       call. = FALSE)
}

entry_paths <- file.path("analytic_pipeline", a$entrypoint)
source_paths <- file.path("analytic_pipeline", a$implementation)
missing <- c(entry_paths[!file.exists(entry_paths)],
             source_paths[!file.exists(source_paths)])

if (length(missing)) {
  stop("Missing archived file(s):\n", paste(missing, collapse = "\n"),
       call. = FALSE)
}

# Git's blob identifiers are recorded for provenance. This checker verifies
# manifest consistency and file presence; it does not execute models.
cat(sprintf("PASS: %d numbered entrypoints and %d implementations present.\n",
            length(entry_paths), length(source_paths)))
cat("Not checked: licensed data availability, runtime dependencies, numerical reproduction, deployment artifacts.\n")
