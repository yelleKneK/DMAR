# Make an immutable development build for testing dependent packages.
#
# Usage, from the package root:
#   Rscript tools/dev_build.R
#   Rscript tools/dev_build.R --commit=<rev> --dest=<directory>
#
# The destination defaults to the checkout's parent directory. Each artifact
# is named DMAR_<version>_dev-<commit>.tar.gz with a checksum sidecar. It is
# built from git archive using the commit's vignette preparation machinery.
# Existing artifacts are never overwritten. This command runs no package
# checks and does not certify an artifact for release.
#
# Install development builds into isolated libraries. tools/ is excluded
# from the source package.

args <- commandArgs(trailingOnly = TRUE)
known <- c("--commit", "--dest")
usage <- "usage: Rscript tools/dev_build.R [--commit=<rev>] [--dest=DIR]"
bad <- args[!sub("=.*", "", args) %in% known | !grepl("=.", args)]
if (length(bad)) { cat("unknown or incomplete argument(s):", bad, "\n", usage, "\n"); quit(status = 2L) }
flag <- function(name, default) {
  v <- sub(paste0("^", name, "="), "", grep(paste0("^", name, "="), args, value = TRUE))
  if (length(v)) v[length(v)] else default
}
rev <- flag("--commit", NA_character_)
dest <- flag("--dest", normalizePath("..", mustWork = TRUE))
refuse <- function(...) { cat("REFUSED: ", ..., "\n", sep = ""); quit(status = 1L) }

if (!file.exists("DESCRIPTION") || !file.exists(".git")) {
  cat("Run from the root of the package's git checkout.\n"); quit(status = 2L)
}
if (!dir.exists(dest)) { cat("the destination", dest, "does not exist\n"); quit(status = 2L) }
dest <- normalizePath(dest)
git <- function(...) suppressWarnings(system2("git", c(...), stdout = TRUE, stderr = FALSE))

## ---- The commit -------------------------------------------------------------------------
if (is.na(rev)) {
  dirty <- git("status", "--porcelain", "--untracked-files=no")
  if (length(dirty) || !is.null(attr(dirty, "status")))
    refuse("the working tree has uncommitted changes to tracked files, which a build named for HEAD would not ",
           "contain:\n", paste0("  ", trimws(head(dirty, 20L)), collapse = "\n"),
           "\nCommit them, or name the commit to build as committed with --commit=<rev>.")
  rev <- "HEAD"
}
sha <- git("rev-parse", "--verify", "--quiet", paste0(rev, "^{commit}"))
if (length(sha) != 1L || !nzchar(sha)) { cat("no commit", rev, "in this repository\n"); quit(status = 2L) }
short <- substr(sha, 1L, 8L)
subject <- git("log", "-1", "--format=%s", sha)
desc <- read.dcf(textConnection(git("show", paste0(sha, ":DESCRIPTION"))), fields = c("Package", "Version"))
pkg_name <- unname(desc[1L, "Package"]); version <- unname(desc[1L, "Version"])
name <- sprintf("%s_%s_dev-%s.tar.gz", pkg_name, version, short)
final <- file.path(dest, name)
sidecar <- paste0(final, ".sha256")

## ---- Never over an existing build ---------------------------------------------------------
taken <- function() {
  if (file.exists(final))
    refuse(sprintf("%s already exists (sha256 %s, %.1f MB, modified %s). A development build is never rebuilt in place: use it as it is (shasum -a 256 -c %s, run in %s).",
                   final, unname(tools::sha256sum(final)), file.size(final) / 2^20,
                   format(file.mtime(final), "%Y-%m-%d %H:%M"), basename(sidecar), dest))
  if (file.exists(sidecar))
    refuse(sprintf("%s exists without its tarball; it is left as it is (move it away to build this commit)", sidecar))
}
taken()

## ---- Build, as the gate builds -------------------------------------------------------------
work <- tempfile("dmar_dev_build_")
build_dir <- file.path(work, pkg_name)
dir.create(build_dir, recursive = TRUE)   # under tempdir(), which R removes when this script ends
st <- system(sprintf("git archive %s | tar -x -C %s", sha, shQuote(build_dir)))
if (st != 0L) refuse("git archive of ", short, " failed")
twins <- sort(sub("[.]Rmd$", "", list.files(file.path(build_dir, "vignettes"), pattern = "[.]Rmd$")))
# A commit from before 2026-11-13 has no tools/asis_vignettes.R; its gate let
# R CMD build render the vignettes, and so does this script (the SIEVER pins
# of early November are such commits).
asis_script <- file.path(build_dir, "tools", "asis_vignettes.R")
asis <- if (file.exists(asis_script)) tryCatch({
  source(asis_script, local = TRUE)
  asis_vignettes(build_dir)
}, error = function(e) e)
if (inherits(asis, "error")) refuse("tools/asis_vignettes.R of ", short, " failed: ", conditionMessage(asis))
if (!is.null(asis) && !setequal(asis$name, twins))
  refuse("vignette twins without prebuilt HTML: ", paste(setdiff(twins, asis$name), collapse = ", "))
log <- file.path(work, "build.log")
st <- withr::with_dir(work, system2("R", c("CMD", "build", pkg_name), stdout = log, stderr = log))
built <- file.path(work, sprintf("%s_%s.tar.gz", pkg_name, version))
if (st != 0L || !file.exists(built))
  refuse("R CMD build failed:\n", paste0("  ", tail(readLines(log, warn = FALSE), 8L), collapse = "\n"))
sha256 <- unname(tools::sha256sum(built))

## ---- Deliver under a temporary name, then link to the final name if it is still free ------
partial <- file.path(dest, sprintf(".%s.partial-%d", name, Sys.getpid()))
if (!isTRUE(file.copy(built, partial, overwrite = FALSE))) refuse("the copy to ", partial, " failed")
if (!identical(unname(tools::sha256sum(partial)), sha256)) { unlink(partial); refuse("the copy at ", partial, " does not have the built sha256") }
taken()   # another build of this commit may have finished meanwhile
linked <- suppressWarnings(file.link(partial, final))
unlink(partial)
if (!isTRUE(linked)) refuse("could not create ", final, " (a build of this commit may have been delivered meanwhile)")
writeLines(paste0(sha256, "  ", name), sidecar)
ok <- identical(unname(tools::sha256sum(final)), sha256) &&
  identical(readLines(sidecar, warn = FALSE), paste0(sha256, "  ", name))
if (!ok) refuse(final, " or its sidecar does not carry the built sha256; do not use it")

cat(sprintf("DEVELOPMENT BUILD, NOT CERTIFIED: %s %s at commit %s (%s)\n", pkg_name, version, short, subject),
    sprintf("  tarball: %s (%.1f MB)\n", final, file.size(final) / 2^20),
    sprintf("  sha256:  %s\n", sha256),
    sprintf("  sidecar: %s (shasum -a 256 -c reads it)\n", sidecar),
    if (!is.null(asis))
      sprintf(paste0("Built as tools/release_gate.R builds (git archive of the commit, vignettes prebuilt by the commit's\n",
                     "tools/asis_vignettes.R, %d of %d, then R CMD build), but nothing was checked: no R CMD check, no test,\n"),
              nrow(asis), length(twins))
    else paste0("Built as the gate of this commit built (git archive of the commit, then R CMD build, which rendered its\n",
                sprintf("%d vignettes; the commit predates tools/asis_vignettes.R), but nothing was checked: no R CMD check, no test,\n",
                        length(twins))),
    "no gate. It is for testing dependent packages against this commit; it is never uploaded to CRAN. Only a\n",
    sprintf("full gate run certifies a tarball, as %s_%s.tar.gz with its record in Quality Control/%s.\n",
            pkg_name, version, pkg_name), sep = "")
quit(status = 0L)
