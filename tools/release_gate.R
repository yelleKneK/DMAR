# The release gate for DMAR: the check CRAN runs, run here first, in full.
#
# No tarball leaves this machine without a full-mode run of this script
# printing ALL GATES PASSED and writing its record (with the tarball's
# sha256) to the quality control folder. The gate exists because two
# CRAN returns (2026-07-31 check time, 2026-08-19 LaTeX in a help page)
# were both detectable locally and both slipped past checks run with a
# step disabled. Nothing here is optional and nothing is skipped for
# speed in full mode.
#
# Usage, from the package root:
#   Rscript tools/release_gate.R            # full mode, before any tarball
#   Rscript tools/release_gate.R --quick    # after any documentation edit
#
# Quick mode runs the static invariants, regenerates the Rd files and
# refuses stale documentation, renders the whole PDF manual through
# LaTeX (the step that catches Rd markup inside math), and runs the
# spelling and pkgdown gates. Full mode adds: a clean working tree, the
# tarball built from the committed tree, R CMD check --as-cran with both
# manuals and the remote URL checks, the complete local test suite, the
# oracle harness, and a sweep of the tarball and git history for
# anything that should not be in them. tools/ is excluded from the
# tarball by .Rbuildignore.

args  <- commandArgs(trailingOnly = TRUE)
quick <- "--quick" %in% args
dest  <- sub("^--dest=", "", grep("^--dest=", args, value = TRUE))
if (!length(dest)) dest <- path.expand("~/Desktop")
qc_dir <- "/Users/kkelley/Dropbox/R/Quality Control/DMAR"

root <- normalizePath(".")
if (!file.exists(file.path(root, "DESCRIPTION"))) stop("Run from the package root.", call. = FALSE)
pkg <- read.dcf("DESCRIPTION", fields = c("Package", "Version"))
tarball_name <- sprintf("%s_%s.tar.gz", pkg[, "Package"], pkg[, "Version"])

results <- list()
gate <- function(name, ok, detail = character()) {
  status <- if (isTRUE(ok)) "PASS" else "FAIL"
  cat(sprintf("[%s] %s\n", status, name))
  for (d in detail) cat("       ", d, "\n", sep = "")
  results[[length(results) + 1L]] <<- list(name = name, ok = isTRUE(ok), detail = detail)
  invisible(isTRUE(ok))
}
run <- function(cmd, args, log) {
  status <- system2(cmd, args, stdout = log, stderr = log)
  list(status = status, lines = readLines(log, warn = FALSE))
}
work <- file.path(tempdir(), "release_gate")
dir.create(work, showWarnings = FALSE, recursive = TRUE)
t0 <- Sys.time()

## ---- 1. Static invariants (the recorded release_checks.R set) ----------
inv_log <- file.path(work, "release_checks.log")
inv <- run("Rscript", "tools/release_checks.R", inv_log)
inv_detail <- character()
if (inv$status != 0L) {
  keep <- FALSE
  for (l in inv$lines) {
    if (grepl("^\\[FAIL\\]", l)) { keep <- TRUE; inv_detail <- c(inv_detail, l); next }
    if (grepl("^\\[PASS\\]", l)) { keep <- FALSE; next }
    if (keep && nzchar(trimws(l))) inv_detail <- c(inv_detail, trimws(l))
  }
}
gate("static invariants (tools/release_checks.R)", inv$status == 0L, inv_detail)

## ---- 2. Rd markup inside math: the 2026-08-19 return ------------------
# Parsed, not grepped: the markup and the \eqn{ it sits in can be on
# different lines (drinks_trial.Rd was exactly that case).
math_markup <- "\\\\(code|emph|link|strong|pkg|bold|var|samp|file|env|dQuote|sQuote|href|url)\\{"
math_nodes <- function(rd) {
  out <- character()
  walk <- function(x) {
    tag <- attr(x, "Rd_tag")
    if (!is.null(tag) && tag %in% c("\\eqn", "\\deqn")) out <<- c(out, paste(unlist(x), collapse = ""))
    if (is.list(x)) for (el in x) walk(el)
  }
  walk(rd); out
}
rd_db <- tools::Rd_db(dir = ".")
hits <- names(Filter(function(rd) any(grepl(math_markup, math_nodes(rd))), rd_db))
gate("no Rd markup inside \\eqn{} or \\deqn{} (parsed, all pages)", length(hits) == 0L, hits)

## ---- 3. Documentation current ------------------------------------------
doc_files <- c(list.files("man", pattern = "[.]Rd$", full.names = TRUE), "NAMESPACE")
before <- tools::md5sum(doc_files)
suppressMessages(roxygen2::roxygenise(roclets = c("rd", "namespace")))
doc_files <- c(list.files("man", pattern = "[.]Rd$", full.names = TRUE), "NAMESPACE")
after <- tools::md5sum(doc_files)
stale <- c(setdiff(names(after), names(before)), setdiff(names(before), names(after)),
           names(before)[names(before) %in% names(after) & before != after[names(before)]])
gate("documentation current (roxygen regenerates nothing)", length(stale) == 0L, stale)

## ---- 4. The PDF manual renders through LaTeX -----------------------------
pdf_log <- file.path(work, "Rd2pdf.log")
pdf_out <- file.path(work, "manual.pdf")
unlink(pdf_out)
pdf <- run("R", c("CMD", "Rd2pdf", "--no-preview", "--force",
                  paste0("--output=", pdf_out), "."), pdf_log)
latex_err <- grep("^! |LaTeX errors|Missing \\$ inserted", pdf$lines, value = TRUE)
gate("PDF manual renders (R CMD Rd2pdf, full LaTeX pass)",
     pdf$status == 0L && file.exists(pdf_out) && length(latex_err) == 0L,
     c(if (pdf$status != 0L) sprintf("Rd2pdf exit status %d", pdf$status), head(latex_err, 10)))

## ---- 5. Spelling and pkgdown --------------------------------------------
sp <- spelling::spell_check_package(".")
gate("spelling (inst/WORDLIST)", nrow(sp) == 0L,
     if (nrow(sp)) sprintf("%s: %s", sp$word, vapply(sp$found, paste, "", collapse = ", ")))
pk <- tryCatch({ pkgdown::check_pkgdown("."); TRUE }, error = function(e) conditionMessage(e))
gate("pkgdown reference index complete", isTRUE(pk), if (!isTRUE(pk)) pk)

if (quick) {
  ok_all <- all(vapply(results, `[[`, logical(1), "ok"))
  cat(sprintf("\nquick mode: %d gates, %d failed, %.0f s\n", length(results),
              sum(!vapply(results, `[[`, logical(1), "ok")),
              as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  if (!ok_all) quit(status = 1L)
  cat("QUICK GATES PASSED (not a release: run full mode before any tarball leaves)\n")
  quit(status = 0L)
}

## ---- 6. Clean tree; the tarball is the committed tree ------------------
dirty <- system2("git", c("status", "--porcelain"), stdout = TRUE)
gate("working tree clean (the gate certifies a commit, not a draft)", length(dirty) == 0L, dirty)
sha <- system2("git", c("rev-parse", "--short", "HEAD"), stdout = TRUE)
src <- file.path(work, "src"); unlink(src, recursive = TRUE); dir.create(file.path(src, pkg[, "Package"]), recursive = TRUE)
system(sprintf("git archive HEAD | tar -x -C %s", shQuote(file.path(src, pkg[, "Package"]))))
old <- setwd(src)
bld <- run("R", c("CMD", "build", pkg[, "Package"]), file.path(work, "build.log"))
setwd(old)
tarball <- file.path(src, tarball_name)
gate(sprintf("tarball built from commit %s", sha), bld$status == 0L && file.exists(tarball),
     if (bld$status != 0L) tail(bld$lines, 5))

## ---- 7. R CMD check --as-cran, manuals on, remote URL checks on --------
chk_log <- file.path(work, "check.log")
old <- setwd(src)
Sys.setenv(`_R_CHECK_CRAN_INCOMING_REMOTE_` = "true", `_R_CHECK_FORCE_SUGGESTS_` = "true")
chk <- run("R", c("CMD", "check", "--as-cran", tarball_name), chk_log)
setwd(old)
log <- chk$lines
status_line <- grep("^Status:", log, value = TRUE)
# NOTEs allowed for a new submission, and nothing else.
note_heads <- grep("\\.\\.\\. .*NOTE$", log, value = TRUE)
allowed_note <- function(h) grepl("CRAN incoming feasibility", h) ||
  (grepl("HTML version of manual", h) && any(grepl("recent enough HTML Tidy", log)))
bad_notes <- note_heads[!vapply(note_heads, allowed_note, logical(1))]
feas <- log[seq(from = max(1, grep("CRAN incoming feasibility", log)[1]),
                length.out = 25)]
feas_extra <- grep("invalid file URIs|Non-FOSS|Insufficient package version|Size of tarball|Days since last update",
                   feas, value = TRUE)
# URL problems: take the whole block, then re-probe any 5xx (a gateway
# timeout at the checker's moment is environmental); only a reproducible
# problem fails the gate, and every flagged URL is recorded either way.
url_start <- grep("invalid URLs", log)
url_block <- character(); url_bad <- character()
if (length(url_start)) {
  k <- url_start[1] + 1L
  while (k <= length(log) && grepl("^\\s", log[k])) { url_block <- c(url_block, trimws(log[k])); k <- k + 1L }
  urls <- sub("^URL: ", "", grep("^URL: ", url_block, value = TRUE))
  statuses <- sub("^Status: ", "", grep("^Status: ", url_block, value = TRUE))
  for (i in seq_along(urls)) {
    st <- if (i <= length(statuses)) statuses[i] else "?"
    u <- sub(" \\(moved to .*\\)$", "", urls[i])
    if (grepl("^5", st)) {
      code <- tryCatch(system2("curl", c("-s", "-o", "/dev/null", "-w", "%{http_code}", "-L", "--max-time", "20", shQuote(u)), stdout = TRUE), error = function(e) "000")
      if (identical(code, "200")) { cat("        transient ", st, " at check time; re-probe returned 200: ", u, "\n", sep = ""); next }
    }
    url_bad <- c(url_bad, sprintf("%s [%s]", urls[i], st))
  }
}
gate("R CMD check --as-cran (manuals built): no ERROR, no WARNING",
     chk$status == 0L && !any(grepl("ERROR|WARNING", status_line)),
     c(status_line, grep("WARNING$|ERROR$", log, value = TRUE)))
gate("R CMD check: no NOTE beyond new-submission and spelling",
     length(bad_notes) == 0L && length(feas_extra) == 0L && length(url_bad) == 0L,
     c(bad_notes, feas_extra, if (length(url_bad)) c("reproducible URL problems:", url_bad)))
gate("PDF manual built by check", any(grepl("checking PDF version of manual \\.\\.\\..*OK", log)),
     grep("PDF version of manual", log, value = TRUE))
gate("HTML manual: no math rendering problems", !any(grepl("math rendering problems", log)),
     grep("math rendering", log, value = TRUE))
if (any(grepl("recent enough HTML Tidy", log)))
  cat("        note: HTML validation skipped locally (old tidy); win-builder validates it\n")
timings <- grep("checking (examples|tests|re-building of vignette outputs) \\.\\.\\.", log, value = TRUE)
for (t in timings) cat("        ", t, "\n", sep = "")

## ---- 8. The full local test suite, Monte Carlo blocks included ----------
Sys.setenv(NOT_CRAN = "true")
tst <- tryCatch(as.data.frame(testthat::test_local(".", reporter = "silent", stop_on_failure = FALSE)),
                error = function(e) NULL)
tst_ok <- !is.null(tst) && sum(tst$failed) == 0L && sum(tst$error) == 0L && sum(tst$warning) == 0L
gate(sprintf("full local test suite: %s passed, %s failed, %s errors, %s warnings, %s skipped",
             if (is.null(tst)) "?" else sum(tst$passed), if (is.null(tst)) "?" else sum(tst$failed),
             if (is.null(tst)) "?" else sum(tst$error), if (is.null(tst)) "?" else sum(tst$warning),
             if (is.null(tst)) "?" else sum(tst$skipped)), tst_ok)

## ---- 9. Oracle harness, against the tarball just checked ----------------
# The snippets call DMAR:: explicitly, so they must see this tarball, not
# whatever version happens to be installed in the user library.
gate_lib <- file.path(work, "lib"); unlink(gate_lib, recursive = TRUE); dir.create(gate_lib)
inst <- run("R", c("CMD", "INSTALL", "--no-docs", "--no-multiarch", paste0("--library=", gate_lib), tarball),
            file.path(work, "install.log"))
gate("tarball installs into a scratch library for the oracle run", inst$status == 0L,
     if (inst$status != 0L) tail(inst$lines, 8))
orc_log <- file.path(work, "oracle.log")
old_libs <- Sys.getenv("R_LIBS")
Sys.setenv(R_LIBS = paste(c(gate_lib, .libPaths()), collapse = .Platform$path.sep))
orc <- run("Rscript", "tools/oracle_checks.R", orc_log)
Sys.setenv(R_LIBS = old_libs)
gate("oracle harness (live cross-package anchors)", orc$status == 0L,
     grep("oracle checks:|^Error", orc$lines, value = TRUE))

## ---- 10. Nothing in the tarball or the history that should not be ------
# The pattern is assembled from character codes so this file never carries
# the literal strings it polices.
banned <- paste(c(intToUtf8(c(99,108,97,117,100,101)), intToUtf8(c(97,110,116,104,114,111,112,105,99)),
                  intToUtf8(c(99,111,45,97,117,116,104,111,114,101,100))), collapse = "|")
names_hit <- system(sprintf("tar -tzf %s | grep -icE '%s'", shQuote(tarball), banned), intern = TRUE)
content_hit <- system(sprintf("tar -xzOf %s 2>/dev/null | grep -icE '%s'", shQuote(tarball), banned), intern = TRUE)
log_hit <- system(sprintf("git log --format=%%B | grep -icE '%s'", banned), intern = TRUE)
# grep -c exits 1 on a zero count, and system(intern = TRUE) then tags the
# "0" with a status attribute; compare the numbers, not the objects.
counts <- suppressWarnings(as.integer(c(names_hit, content_hit, log_hit)))
gate("tarball and git history carry no authorship residue",
     length(counts) == 3L && !anyNA(counts) && sum(counts) == 0L,
     sprintf("file names %s, file content %s, commit messages %s", names_hit, content_hit, log_hit))
surfaces <- unique(sub("/.*", "", sub("^[^/]+/", "", system(sprintf("tar -tzf %s", shQuote(tarball)), intern = TRUE))))
gate("tarball holds only the standard surfaces", all(surfaces %in% c("DESCRIPTION", "NAMESPACE", "NEWS.md", "README.md", "R", "build", "data", "inst", "man", "tests", "vignettes", "LICENSE", "LICENSE.md")),
     setdiff(surfaces, c("DESCRIPTION", "NAMESPACE", "NEWS.md", "README.md", "R", "build", "data", "inst", "man", "tests", "vignettes", "LICENSE", "LICENSE.md")))

## ---- Record and hand-off ---------------------------------------------------
ok_all <- all(vapply(results, `[[`, logical(1), "ok"))
sha256 <- unname(tools::sha256sum(tarball))
elapsed <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
stamp <- format(Sys.time(), "%Y-%m-%d_%H%M")
record <- c(
  sprintf("# %s release gate, %s, commit %s", pkg[, "Package"], format(Sys.time(), "%Y-%m-%d %H:%M"), sha),
  "", sprintf("Result: %s (%d gates, %.1f minutes)", if (ok_all) "ALL GATES PASSED" else "FAILED",
              length(results), elapsed), "",
  sprintf("Tarball: %s", tarball_name), sprintf("sha256: %s", sha256), "",
  unlist(lapply(results, function(r) c(sprintf("- [%s] %s", if (r$ok) "PASS" else "FAIL", r$name),
                                       if (length(r$detail)) paste0("    ", r$detail)))),
  "", "Check timings:", paste0("    ", timings), "",
  paste0("Check status: ", status_line))
rec_path <- file.path(qc_dir, sprintf("%s_release_gate_%s_%s.md", pkg[, "Package"], stamp, sha))
writeLines(record, rec_path)
cat(sprintf("\nrecord: %s\n", rec_path))
if (ok_all) {
  file.copy(tarball, file.path(dest, tarball_name), overwrite = TRUE)
  cat(sprintf("tarball: %s\nsha256:  %s\n", file.path(dest, tarball_name), sha256))
  cat(sprintf("\nALL GATES PASSED (%d gates, %.1f min)\n", length(results), elapsed))
  quit(status = 0L)
}
cat(sprintf("\nGATE FAILED: %d of %d gates failed; no tarball handed off\n",
            sum(!vapply(results, `[[`, logical(1), "ok")), length(results)))
quit(status = 1L)
