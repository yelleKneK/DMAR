# Live cross-package oracle checks for DMAR.
#
# Every numerical anchor that the shipped test suite pins as a constant was
# captured from a live call to another package (MBESS, metafor, semTools,
# emmeans, multcomp, irr, irrCAC, psych, gsl, BayesFactor, mirt, sirt, lme4
# companions, and others since removed from Suggests). The pinned constants
# are the operative guard on CRAN and for users; this script re-runs the
# LIVE comparisons on the maintainer's machine, where those packages exist,
# so drift between DMAR and its oracles is caught at release time without
# any of them being a dependency.
#
# Usage, from the package root, before each release:
#   Rscript tools/oracle_checks.R
# Each snippet file under tools/oracle_snippets/ is self-contained: it
# installs nothing, states which oracle package and version it was written
# against, and ends every block in stopifnot(). A missing oracle package
# skips its file with a message rather than failing, so the script degrades
# gracefully on machines that lack one.
#
# This directory is excluded from the tarball by .Rbuildignore (^tools$).

files <- sort(list.files(file.path("tools", "oracle_snippets"),
                         pattern = "[.]R$", full.names = TRUE))
if (!length(files)) stop("run from the package root; no snippets found")
ok <- 0L; skipped <- 0L; failed <- character()
for (f in files) {
  message("== ", basename(f))
  res <- tryCatch({ source(f, local = new.env()); "ok" },
                  error = function(e) conditionMessage(e))
  if (identical(res, "ok")) { ok <- ok + 1L; next }
  if (grepl("there is no package called", res)) {
    skipped <- skipped + 1L
    message("   skipped (oracle not installed): ", sub(".*called ", "", res))
  } else {
    failed <- c(failed, paste0(basename(f), ": ", res))
  }
}
message(sprintf("oracle checks: %d ok, %d skipped, %d failed",
                ok, skipped, length(failed)))
if (length(failed)) stop(paste(failed, collapse = "\n"), call. = FALSE)
