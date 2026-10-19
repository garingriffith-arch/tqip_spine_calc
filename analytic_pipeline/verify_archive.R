# Verify the public source archive without reading licensed ACS data.
# Run from the repository root:
#   Rscript analytic_pipeline/verify_archive.R

root <- getwd()
base <- file.path(root, "analytic_pipeline")
pipeline <- read.csv(file.path(base, "PIPELINE_MANIFEST.csv"),
                     stringsAsFactors = FALSE, check.names = FALSE)
archive <- read.csv(file.path(base, "SOURCE_ARCHIVE_MANIFEST.csv"),
                    stringsAsFactors = FALSE, check.names = FALSE)

required <- c("step", "entrypoint", "implementation", "git_blob_sha1")
stopifnot(all(required %in% names(archive)),
          all(c("step", "entrypoint", "implementation") %in% names(pipeline)))
stages <- pipeline[pipeline$step != "PENDING", , drop = FALSE]
if (nrow(stages) != 23L || nrow(archive) != 23L) {
  stop("Expected 23 archived sources (steps 01–21, 31 and 32).", call. = FALSE)
}
if (anyDuplicated(stages$step) || anyDuplicated(archive$step) ||
    anyDuplicated(archive$implementation)) {
  stop("Duplicate source-archive step or implementation.", call. = FALSE)
}
idx <- match(stages$step, archive$step)
if (anyNA(idx)) stop("Pipeline and source manifests disagree.", call. = FALSE)
archive <- archive[idx, , drop = FALSE]
if (!identical(stages$entrypoint, archive$entrypoint) ||
    !identical(stages$implementation, archive$implementation)) {
  stop("Manifest paths disagree.", call. = FALSE)
}
if (any(!grepl("^[0-9a-f]{40}$", archive$git_blob_sha1))) {
  stop("Invalid Git blob SHA-1 in source manifest.", call. = FALSE)
}
entry <- file.path(base, archive$entrypoint)
impl <- file.path(base, archive$implementation)
missing <- c(entry[!file.exists(entry)], impl[!file.exists(impl)])
if (length(missing)) {
  stop("Missing file(s):\n", paste(missing, collapse = "\n"), call. = FALSE)
}
# An entrypoint must refer to the implementation recorded in both manifests.
for (i in seq_along(entry)) {
  lines <- readLines(entry[i], warn = FALSE)
  if (!any(grepl(basename(impl[i]), lines, fixed = TRUE))) {
    stop("Entrypoint does not reference its recorded implementation: ", entry[i],
         call. = FALSE)
  }
}

# The historical deployment manifest field named source_file_md5 is a frozen
# DATASET checksum, not the Stage 21 R script MD5. Never compare it to the script.
stage21_path <- file.path(base, "implementation",
 "21_build_SPINE_TRACT_deployment_2020_2024_FIXED_20261019.R")
stage21_md5 <- unname(tools::md5sum(stage21_path))
if (!identical(stage21_md5, "e2e81bb9c273b5f182517601796a511f")) {
  stop("Recovered Stage 21 source MD5 differs from archived file.", call. = FALSE)
}

# Require Git and verify file content, not just manifest shape.
git <- Sys.which("git")
if (!nzchar(git)) stop("Git is required for source SHA verification.", call. = FALSE)
for (i in seq_along(impl)) {
  rel <- file.path("analytic_pipeline", archive$implementation[i])
  result <- suppressWarnings(system2(git,
             args = c("hash-object", shQuote(rel)),
             stdout = TRUE, stderr = TRUE))
  code <- attr(result, "status")
  if ((!is.null(code) && code != 0L) ||
      length(result) != 1L ||
      !identical(trimws(result), archive$git_blob_sha1[i])) {
    stop("Source checksum mismatch: ", rel, call. = FALSE)
  }
}
pending <- pipeline[pipeline$step == "PENDING", , drop = FALSE]
cat(sprintf("PASS: %d entrypoints and source SHA-1 values verified.\n",
            nrow(archive)))
if (nrow(pending)) {
  cat("Not archived: ", paste(pending$entrypoint, collapse = ", "), "\n",
      sep = "")
}
cat("Not assessed: ACS data, local environment, numerical reproduction, model training.\n")
