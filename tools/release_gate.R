# Release validation for DMAR, run from the package root.
#
# Usage:
#   caffeinate -i Rscript tools/release_gate.R            # full release validation
#   Rscript tools/release_gate.R --quick                 # development checks
# Options:
#   --dest=DIR         certified-artifact directory (default: parent checkout directory)
#   --replace-stale    explicitly permit replacement of an existing versioned artifact
#
# Full mode builds from the committed tree, checks both manuals, examples,
# vignettes, the complete test suite, numerical oracles and dependent contracts.
# Quick mode checks source invariants, documentation, examples and detector
# controls. Full mode requires a clean tree and an awake machine on macOS.
# Reports go to DMAR_QC_DIR when set, otherwise to Quality Control/DMAR
# beside the checkout. tools/ is excluded from the source package.

args  <- commandArgs(trailingOnly = TRUE)
# A mistyped flag stops here rather than starting a full run without it.
known_flags <- c("--quick", "--dest", "--replace-stale")
unknown <- args[!sub("=.*", "", args) %in% known_flags | (grepl("^--dest", args) & !grepl("^--dest=.", args))]
if (length(unknown)) {
  cat("unknown argument(s):", unknown, "\n",
      "usage: Rscript tools/release_gate.R [--quick] [--dest=DIR] [--replace-stale]\n")
  quit(status = 2L)
}
quick <- "--quick" %in% args
replace_stale <- "--replace-stale" %in% args
dest  <- sub("^--dest=", "", grep("^--dest=", args, value = TRUE))
# Default destinations are relative to the package checkout and remain configurable.
if (!length(dest)) dest <- normalizePath("..", mustWork = TRUE)
qc_dir <- Sys.getenv("DMAR_QC_DIR",
                     unset = file.path(normalizePath("..", mustWork = TRUE), "Quality Control", "DMAR"))

root <- normalizePath(".")
if (!file.exists(file.path(root, "DESCRIPTION"))) stop("Run from the package root.", call. = FALSE)
pkg <- read.dcf("DESCRIPTION", fields = c("Package", "Version", "Date"))
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

## ---- Calibration: the accepted release, measured here and on the farm ----
# Local timings of the published 1.0.0 tree on this machine (2026-09-22,
# awake, examples and the CRAN test subset exactly as the check runs
# them; the network-bound incoming feasibility step is excluded because
# its time is the network's), and the farm's total check times for the
# same release from tools/cran_status.R. Full mode fails a check whose
# local time grows past the ratio that projects the r-devel Linux farm
# total to the incoming pretest's 10 minute line (section 7a), or any
# timed step over its budget (Linux runs about 1.8 times, Intel macOS
# about 3 times the local wall clock). Re-measure and update this block
# when a new version is published.
accepted <- list(
  measured     = "2026-09-22",
  total_s      = 185,    # R CMD check wall clock less the feasibility step
  step_s       = c(examples = 62, tests = 12, vignettes = 15, manual = 17),
  # The Linux and macOS totals were read on 2026-09-22. The r-devel-windows
  # total is 1.0.0's on CRAN's check results page as read on 2026-11-15
  # (R-devel of that week); it was added then, after the 1.1.0 pretest
  # returned the overall check time from that flavor. At the 1.1.0 check
  # of 2026-11-15 (266 s locally) it projects 1022 s, and the pretest
  # measured 17 minutes.
  farm_total_s = c(`r-devel-linux-x86_64-fedora-gcc` = 403,
                   `r-release-macos-x86_64` = 681,
                   `r-devel-windows-x86_64` = 711)
)
# Project the part of a Windows check not explained by examples, tests
# or vignettes. These timings compare the 2026-11-15 incoming pretest with
# the local check of the same artifact. The 17-minute rounded total gives
# a floor estimate with about 30 seconds of resolution. The projection
# is informational; it does not change the pass/fail decision.
windows_pretest <- list(
  measured     = "2026-11-15",
  flavor       = "r-devel-windows-x86_64",
  overall_s    = 17 * 60,
  win_s        = c(r_code = 238, rd = 27, examples = 279, tests = 25, pdf_manual = 120,
                   html_manual = 121, feasibility = 38),
  local_cpu_s     = c(r_code = 51, rd = 6, examples = 78, tests = 6, vignettes = 6, pdf_manual = 32,
                      html_manual = 1, install = 25, feasibility = 15),
  local_elapsed_s = c(r_code = 53, rd = 6, examples = 79, tests = 6, vignettes = 8, pdf_manual = 33,
                      html_manual = 1, install = 26, feasibility = 53),
  local_total_s   = 266   # the whole check (319 s) less the feasibility step (53 s)
)
incoming_limit_s <- 600
# The r-devel-windows lines of section 7a, from a local check's elapsed
# step times (r_code, rd, examples, tests, vignettes, pdf_manual,
# html_manual) and its whole time less the feasibility step. The total is
# projected by the whole-check ratio, as for the other flavors; the floor
# step by step, as described above. Section 5d runs this on the calibration
# check itself as a control.
.gate_windows_projection <- function(elapsed, check_total_s, cal = windows_pretest, limit = incoming_limit_s) {
  flavor <- cal$flavor
  total <- accepted$farm_total_s[[flavor]] * check_total_s / accepted$total_s
  need <- c("r_code", "rd", "examples", "tests", "vignettes", "pdf_manual", "html_manual")
  miss <- need[!need %in% names(elapsed) | is.na(elapsed[need])]
  out <- list(total = total, over_total = isTRUE(total > limit), floor = NA_real_, over_floor = NA)
  info_total <- if (out$over_total)
    sprintf(paste("information, not a failure: the projected %s total %.0f s is over the %d s line, so the",
                  "incoming pretest will report the overall check time (the %s pretest reported %.0f minutes)"),
            flavor, total, limit, cal$measured, cal$overall_s / 60)
  if (length(miss) || is.na(check_total_s)) {
    out$lines <- c(sprintf("%s floor not computed: this check's log has no time for %s", flavor,
                           paste(c(miss, if (is.na(check_total_s)) "the whole check"), collapse = ", ")), info_total)
    return(out)
  }
  le <- cal$local_elapsed_s; ws <- cal$win_s
  floor_cal <- cal$overall_s - ws[["examples"]] - ws[["tests"]]
  rest_win <- floor_cal - sum(ws[c("r_code", "rd", "pdf_manual", "html_manual", "feasibility")])
  rest_local <- max(check_total_s - sum(elapsed[need]), 0)
  part <- c(r_code = ws[["r_code"]] * elapsed[["r_code"]] / le[["r_code"]],
            rd = ws[["rd"]] * elapsed[["rd"]] / le[["rd"]],
            pdf_manual = ws[["pdf_manual"]] * elapsed[["pdf_manual"]] / le[["pdf_manual"]],
            html_manual = ws[["html_manual"]] * elapsed[["pdf_manual"]] / le[["pdf_manual"]],
            feasibility = ws[["feasibility"]],
            rest = rest_win * rest_local / (cal$local_total_s - sum(le[need])))
  out$floor <- sum(part)
  out$over_floor <- isTRUE(out$floor > limit)
  out$lines <- c(
    sprintf(paste("%s floor (the check less its examples, tests, and vignettes), projected step by step from the",
                  "%s pretest: %.0f s (R code %.0f, Rd files %.0f, PDF manual %.0f, HTML manual %.0f, incoming",
                  "feasibility %.0f as measured, the rest %.0f); examples %.0f s and tests %.0f s come on top"),
            flavor, cal$measured, out$floor, part[["r_code"]], part[["rd"]], part[["pdf_manual"]],
            part[["html_manual"]], part[["feasibility"]], part[["rest"]],
            ws[["examples"]] * elapsed[["examples"]] / le[["examples"]], ws[["tests"]] * elapsed[["tests"]] / le[["tests"]]),
    info_total,
    if (out$over_floor)
      sprintf(paste("information, not a failure: the projected %s floor %.0f s is over the %d s line by itself, so",
                    "trimming examples and tests cannot clear it; the maintainer can",
                    "request manual inspection with the measured timings"), flavor, out$floor, limit))
  out
}
# Current budgets allow a projected 720-second Linux check and 60 seconds
# for CRAN-path tests. The incoming 600-second line is also printed so a
# maintainer can assess check time before submission.
gate_limit_s <- 720
tests_budget_s <- 60

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

## ---- 5a. Commit messages carry no authorship residue --------------------
# Section 10's sweep of every commit message reachable from HEAD, here too
# because it costs a second and a full run failed on it alone (2026-11-11,
# 9f97925: three commit message lines). The pattern is assembled from
# character codes so this file never carries the literal strings it
# polices; section 10 uses the same pattern and the same function.
.gate_banned <- paste(c(intToUtf8(c(99,108,97,117,100,101)), intToUtf8(c(97,110,116,104,114,111,112,105,99)),
                        intToUtf8(c(99,111,45,97,117,116,104,111,114,101,100))), collapse = "|")
# The count grep -c prints for the commit messages of the repository at
# `dir`, as a string; callers compare the parsed number (grep -c exits 1 on
# a zero count, which system() reports as a warning and a status attribute).
.gate_log_residue <- function(dir = ".") {
  suppressWarnings(system(sprintf("git -C %s log --format=%%B | grep -icE '%s'", shQuote(dir), .gate_banned),
                          intern = TRUE))
}
log_res <- .gate_log_residue(".")
n_log_res <- suppressWarnings(as.integer(log_res))
n_commits <- length(suppressWarnings(system2("git", c("rev-list", "HEAD"), stdout = TRUE, stderr = FALSE)))
log_offenders <- if (!identical(n_log_res, 0L))
  suppressWarnings(system2("git", c("log", "-i", "-E", shQuote(paste0("--grep=", .gate_banned)), "--format='%h %s'"),
                           stdout = TRUE, stderr = FALSE))
gate(sprintf("commit messages carry no authorship residue (%d commits from HEAD; section 10's sweep)", n_commits),
     length(n_log_res) == 1L && identical(n_log_res, 0L) && n_commits > 0L,
     if (!identical(n_log_res, 0L))
       c(sprintf("commit message lines that match: %s", paste(log_res, collapse = " ")),
         paste0("  ", head(log_offenders, 10L)), "reword them (git rebase) before any push or release"))

## ---- 5b. The oracle snippets call what the namespace has ------------------
# On 2026-11-11 a full run failed after half an hour because a snippet still
# called ci_nc_t(t_value = ), an argument the morning's renames had replaced
# with ncp (record DMAR_release_gate_2026-11-11_0420_b29bc3a.md). Here each
# snippet is parsed, not run: every DMAR::f it calls must be an export of
# the NAMESPACE that section 3 just regenerated, every DMAR::x it names
# without calling must be an export or a dataset, every DMAR:::f must be
# defined in R/, and every argument it names must be a formal of the
# function (exactly, or as the unique start of one, as R matches it). A
# function that passes ... on to another (uniroot(), lavaan) can take
# names that are not its formals, so those are checked only against the
# names it refuses: the string literals its body, or a helper it hands
# ... to, tests with %in% against names(list(...)) or ...names(), which is
# how DMAR refuses a 1.0.x argument name. This check finds the 2026-11-11
# failure in that commit's snippets.
.gate_is_empty_arg <- function(x) is.symbol(x) && !nzchar(as.character(x))
# Whether a function body hands ... on to another call. names(list(...)),
# ...length(), ...names(), and missing() read the dots without passing them
# on; any other use, list(...) on its own included, counts as passing them.
.gate_dots_forwarded <- function(body) {
  found <- FALSE
  walk <- function(x) {
    if (found || !is.call(x)) return(invisible())
    h <- if (is.name(x[[1L]])) as.character(x[[1L]]) else ""
    if (h %in% c("...length", "...names", "missing")) return(invisible())
    if (h == "names" && length(x) == 2L && identical(x[[2L]], quote(list(...)))) return(invisible())
    for (k in seq_along(x)[-1L]) {
      if (.gate_is_empty_arg(x[[k]])) next
      if (identical(x[[k]], as.name("..."))) { found <<- TRUE; return(invisible()) }
      walk(x[[k]])
    }
    if (is.call(x[[1L]])) walk(x[[1L]])
    invisible()
  }
  walk(body)
  found
}
# The names a function body refuses in ..., and the functions it hands ... to.
.gate_dots_refused <- function(body) {
  holders <- character(); refused <- character(); helpers <- character()
  is_dots_names <- function(x) identical(x, quote(...names())) || identical(x, quote(names(list(...))))
  strings <- function(x) if (is.character(x)) x else if (is.call(x) && identical(x[[1L]], as.name("c")))
    unlist(lapply(as.list(x)[-1L], function(a) if (is.character(a)) a)) else character()
  walk <- function(x) {
    if (!is.call(x)) return(invisible())
    h <- if (is.name(x[[1L]])) as.character(x[[1L]]) else ""
    if (h %in% c("<-", "=") && length(x) == 3L && is.name(x[[2L]]) && is_dots_names(x[[3L]]))
      holders <<- c(holders, as.character(x[[2L]]))
    if (h == "%in%" && length(x) == 3L &&
        (is_dots_names(x[[3L]]) || (is.name(x[[3L]]) && as.character(x[[3L]]) %in% holders)))
      refused <<- c(refused, strings(x[[2L]]))
    if (nzchar(h) && any(vapply(seq_along(x)[-1L], function(k) identical(x[[k]], as.name("...")), logical(1))))
      helpers <<- c(helpers, h)
    for (k in seq_along(x)) if (!.gate_is_empty_arg(x[[k]])) walk(x[[k]])
    invisible()
  }
  walk(body)
  list(refused = unique(refused), helpers = unique(helpers))
}
# What the source tree at `root` defines: the exports in NAMESPACE, the
# datasets in data/, every name assigned at the top level of R/, and for
# each function so defined its formals and what it does with ... . The
# function expressions are evaluated in the base environment, which makes
# closures and evaluates no default and no body.
.gate_namespace_facts <- function(root) {
  exports <- unlist(lapply(parse(file.path(root, "NAMESPACE"), keep.source = FALSE), function(e)
    if (is.call(e) && identical(e[[1L]], as.name("export"))) vapply(as.list(e)[-1L], as.character, "")))
  datasets <- unlist(lapply(list.files(file.path(root, "data"), pattern = "[.](rda|RData)$", full.names = TRUE),
                            function(f) load(f, envir = new.env())))
  defined <- character(); fmls <- list()
  for (f in list.files(file.path(root, "R"), pattern = "[.][Rr]$", full.names = TRUE)) {
    for (e in parse(f, keep.source = FALSE)) {
      if (!(is.call(e) && is.name(e[[1L]]) && as.character(e[[1L]]) %in% c("<-", "=") && length(e) == 3L)) next
      nm <- if (is.name(e[[2L]])) as.character(e[[2L]]) else if (is.character(e[[2L]])) e[[2L]] else next
      defined <- c(defined, nm)
      if (is.call(e[[3L]]) && identical(e[[3L]][[1L]], as.name("function"))) {
        fn <- eval(e[[3L]], baseenv())
        fmls[[nm]] <- c(list(formals = names(formals(fn)), dots_forwarded = .gate_dots_forwarded(body(fn))),
                        .gate_dots_refused(body(fn)))
      }
    }
  }
  list(exports = exports, datasets = datasets, defined = unique(defined), formals = fmls)
}
.gate_snippet_problems <- function(files, facts) {
  problems <- character(); n_calls <- 0L; n_args <- 0L; n_passed_on <- 0L
  for (f in files) {
    exprs <- tryCatch(parse(f, keep.source = FALSE), error = function(e) e)
    if (inherits(exprs, "error")) {
      problems <- c(problems, sprintf("%s does not parse: %s", basename(f), conditionMessage(exprs))); next
    }
    ns_name <- function(x) if (is.call(x) && length(x) == 3L && is.name(x[[1L]]) &&
                               as.character(x[[1L]]) %in% c("::", ":::") && identical(as.character(x[[2L]]), "DMAR"))
      c(op = as.character(x[[1L]]), name = as.character(x[[3L]]))
    check_name <- function(op, name, called) {
      ok <- if (op == ":::") name %in% facts$defined else name %in% c(facts$exports, if (!called) facts$datasets)
      if (!ok) problems <<- c(problems, sprintf("%s: DMAR%s%s is not %s", basename(f), op, name,
                                                if (op == ":::") "defined in R/" else if (called) "an export"
                                                else "an export or a dataset"))
      ok
    }
    walk <- function(x) {
      if (!is.call(x)) return(invisible())
      head <- ns_name(x[[1L]])
      if (!is.null(head)) {
        n_calls <<- n_calls + 1L
        fm <- if (check_name(head[["op"]], head[["name"]], TRUE)) facts$formals[[head[["name"]]]]
        given <- names(x)[-1L]
        given <- given[!is.na(given) & nzchar(given)]
        n_args <<- n_args + length(given)
        named <- setdiff(fm$formals, "...")
        for (a in given) {
          if (is.null(fm) || a %in% named) next
          if (fm$dots_forwarded) {
            refused <- unique(c(fm$refused, unlist(lapply(fm$helpers, function(h) facts$formals[[h]]$refused))))
            if (a %in% refused) problems <<- c(problems, sprintf("%s: DMAR%s%s() refuses the argument '%s'",
                                                                 basename(f), head[["op"]], head[["name"]], a))
            else n_passed_on <<- n_passed_on + 1L
          } else if (sum(startsWith(named, a)) != 1L) {
            problems <<- c(problems, sprintf("%s: DMAR%s%s() has no argument '%s'", basename(f), head[["op"]],
                                             head[["name"]], a))
          }
        }
      } else {
        own <- ns_name(x)
        if (!is.null(own)) { check_name(own[["op"]], own[["name"]], FALSE); return(invisible()) }
        walk(x[[1L]])
      }
      for (k in seq_along(x)[-1L]) if (!.gate_is_empty_arg(x[[k]])) walk(x[[k]])
      invisible()
    }
    for (e in exprs) walk(e)
  }
  list(problems = unique(problems), n_calls = n_calls, n_args = n_args, n_passed_on = n_passed_on)
}
snippet_files <- sort(list.files(file.path("tools", "oracle_snippets"), pattern = "[.]R$", full.names = TRUE))
ns_facts <- tryCatch(.gate_namespace_facts("."), error = function(e) e)
snip <- if (inherits(ns_facts, "error")) NULL else .gate_snippet_problems(snippet_files, ns_facts)
gate(if (is.null(snip)) "oracle snippets call only what the namespace has (the namespace could not be read)"
     else sprintf(paste("oracle snippets call only what the namespace has (%d files, %d DMAR calls, %d named arguments;",
                        "%d passed on through ... and checked against the names refused there)"),
                  length(snippet_files), snip$n_calls, snip$n_args, snip$n_passed_on),
     !is.null(snip) && length(snippet_files) > 0L && snip$n_calls > 0L && !length(snip$problems),
     c(if (inherits(ns_facts, "error")) conditionMessage(ns_facts), snip$problems))

## ---- 5c. Every example runs from the source tree ---------------------------------
# R CMD check runs every example, but only in full mode; a rename that
# misses an example (the rename protocol runs them by hand) was otherwise
# found half an hour into a full run. The runner below does what
# devtools::run_examples(run_donttest = FALSE, document = FALSE) does
# (load_all() with the exports only, then pkgload::run_example() on each
# help page, sorted without regard to case, each in a fresh environment),
# in a separate R process started in a temporary directory, and goes on
# past a page that fails, so that every failing page is named. On
# 2026-11-15 it took 36 s for 360 pages, and with it and sections 5a, 5b,
# and 5d the quick gate took 114 s against 77 s before (both timed in the
# shared R test slot); a check that added more than about 90 s would have
# been put behind a flag.
.gate_examples_runner <- c(
  "args <- commandArgs(trailingOnly = TRUE)",
  "pkg <- args[1L]",
  "suppressMessages(pkgload::load_all(pkg, reset = TRUE, export_all = FALSE, helpers = FALSE, quiet = TRUE))",
  "files <- if (length(args) > 1L) args[-1L] else list.files(file.path(pkg, 'man'), pattern = '[.][Rr]d$', full.names = TRUE)",
  "files <- files[order(tolower(basename(files)))]",
  "failed <- 0L",
  "for (f in files) {",
  "  r <- tryCatch({ pkgload::run_example(f, run_donttest = FALSE, run_dontrun = FALSE, quiet = TRUE); NULL },",
  "                error = function(e) conditionMessage(e))",
  "  if (!is.null(r)) { failed <- failed + 1L",
  "    cat(sprintf('EXAMPLE ERROR %s: %s\\n', basename(f), gsub('[[:space:]]+', ' ', r))) }",
  "}",
  "cat(sprintf('examples run: %d pages, %d with an error\\n', length(files), failed))")
.gate_run_examples <- function(pkg, rd = character(), log) {
  pkg <- normalizePath(pkg); rd <- normalizePath(rd); log <- normalizePath(log, mustWork = FALSE)
  wd <- tempfile("gate_examples_"); dir.create(wd)
  script <- file.path(wd, "run_examples.R")
  writeLines(.gate_examples_runner, script)
  old <- setwd(wd); on.exit(setwd(old), add = TRUE)
  t1 <- Sys.time()
  st <- system2("Rscript", shQuote(c(script, pkg, rd)), stdout = log, stderr = log,
                env = "PROCESSX_NOTIFY_OLD_SIGCHLD=true")
  lines <- readLines(log, warn = FALSE)
  summary <- grep("^examples run: [0-9]+ pages, [0-9]+ with an error$", lines, value = TRUE)
  counts <- if (length(summary) == 1L) as.integer(regmatches(summary, gregexpr("[0-9]+", summary))[[1L]]) else integer()
  list(status = st, pages = if (length(counts) == 2L) counts[1L] else NA_integer_,
       failed = if (length(counts) == 2L) counts[2L] else NA_integer_,
       errors = grep("^EXAMPLE ERROR ", lines, value = TRUE), tail = tail(lines, 5L),
       seconds = as.numeric(difftime(Sys.time(), t1, units = "secs")))
}
n_rd <- length(list.files("man", pattern = "[.][Rr]d$"))
ex <- .gate_run_examples(".", log = file.path(work, "examples.log"))
gate(sprintf("every example runs from the source tree (%s of %d pages, %s with an error, %.0f s)",
             if (is.na(ex$pages)) "?" else ex$pages, n_rd, if (is.na(ex$failed)) "?" else ex$failed, ex$seconds),
     isTRUE(ex$status == 0L) && identical(ex$pages, n_rd) && identical(ex$failed, 0L),
     c(sub("^EXAMPLE ERROR ", "", head(ex$errors, 20L)),
       if (length(ex$errors) > 20L) sprintf("and %d more", length(ex$errors) - 20L),
       if (is.na(ex$pages)) c("the runner printed no summary line; its last lines:", paste0("  ", ex$tail))))

## ---- 5d. Positive controls ---------------------------------------------------
# Each detector above, and the r-devel-windows projection that section 7a
# prints, is shown a known-bad input, through the function the gate itself
# calls, and must flag it; a detector that passes its planted defect fails
# here, so none passes vacuously.
ctrl_dir <- file.path(work, "controls"); unlink(ctrl_dir, recursive = TRUE); dir.create(ctrl_dir)
ctrl <- list()
# (a) A repository with one clean commit message and then one that carries
# a banned word.
ctrl$commit_messages <- local({
  repo <- file.path(ctrl_dir, "repo"); dir.create(repo)
  g <- function(...) suppressWarnings(system2("git", c("-C", shQuote(repo), "-c", "user.name=gate",
                                                       "-c", "user.email=gate@localhost", "-c", "commit.gpgsign=false",
                                                       "-c", "core.hooksPath=/dev/null", ...), stdout = TRUE, stderr = TRUE))
  g("init", "-q")
  g("commit", "-q", "--allow-empty", "-m", shQuote("A clean message"))
  clean <- suppressWarnings(as.integer(.gate_log_residue(repo)))
  g("commit", "-q", "--allow-empty", "-m",
    shQuote(paste("A planted message naming", intToUtf8(c(65,110,116,104,114,111,112,105,99)))))
  bad <- suppressWarnings(as.integer(.gate_log_residue(repo)))
  list(ok = identical(clean, 0L) && identical(bad, 1L),
       why = sprintf("a clean history counted %s, the planted one %s (want 0 and 1)", clean, bad))
})
# (b) A snippet that calls a missing function, a missing argument, a refused
# argument, a missing internal, and a missing dataset, beside a call that is
# right.
ctrl$oracle_snippets <- local({
  if (inherits(ns_facts, "error")) return(list(ok = FALSE, why = "the namespace could not be read"))
  plain <- Filter(function(n) { fm <- ns_facts$formals[[n]]
    !is.null(fm) && !fm$dots_forwarded && length(setdiff(fm$formals, "...")) > 0L }, sort(ns_facts$exports))
  refusing <- Filter(function(n) { fm <- ns_facts$formals[[n]]; !is.null(fm) && fm$dots_forwarded && length(fm$refused) },
                     sort(ns_facts$exports))
  fn <- plain[1L]; arg <- setdiff(ns_facts$formals[[fn]]$formals, "...")[1L]
  lines <- c(sprintf("DMAR::%s(%s = 1)", fn, arg), "DMAR::no_such_function_planted(1)",
             sprintf("DMAR::%s(planted_argument = 1)", fn), "DMAR:::.no_such_internal_planted()",
             "x <- DMAR::no_such_dataset_planted",
             if (length(refusing)) sprintf("DMAR::%s(%s = 1)", refusing[1L], ns_facts$formals[[refusing[1L]]]$refused[1L]))
  f <- file.path(ctrl_dir, "planted_snippet.R"); writeLines(lines, f)
  got <- .gate_snippet_problems(f, ns_facts)$problems
  want <- c("no_such_function_planted", "'planted_argument'", ".no_such_internal_planted", "no_such_dataset_planted",
            if (length(refusing)) sprintf("%s() refuses", refusing[1L]))
  found <- vapply(want, function(w) any(grepl(w, got, fixed = TRUE)), logical(1))
  list(ok = all(found) && length(got) == length(want),
       why = sprintf("%d of %d planted defects flagged, %d flags in all%s", sum(found), length(want), length(got),
                  if (length(refusing)) "" else "; no export refuses a name in ..., so that case was not planted"))
})
# (c) A help page whose example stops.
ctrl$examples <- local({
  rd <- file.path(ctrl_dir, "planted_example.Rd")
  writeLines(c("\\name{planted_example}", "\\alias{planted_example}", "\\title{A planted failing example}",
               "\\description{A planted failing example.}", "\\examples{", "stop(\"planted failure\")", "}"), rd)
  r <- .gate_run_examples(".", rd = rd, log = file.path(ctrl_dir, "examples.log"))
  list(ok = identical(r$pages, 1L) && identical(r$failed, 1L) &&
         any(grepl("planted_example.Rd: planted failure", r$errors, fixed = TRUE)),
       why = sprintf("the planted page ran as %s page with %s error", r$pages, r$failed))
})
# (d) The r-devel-windows lines of section 7a, on the local check the
# floor is calibrated on (they must give back the pretest's floor and say
# both things) and on a check of half its size (they must say neither).
ctrl$windows_projection <- local({
  at <- .gate_windows_projection(windows_pretest$local_elapsed_s, windows_pretest$local_total_s)
  half <- .gate_windows_projection(windows_pretest$local_elapsed_s / 2, windows_pretest$local_total_s / 2)
  floor_cal <- windows_pretest$overall_s - sum(windows_pretest$win_s[c("examples", "tests")])
  list(ok = isTRUE(abs(at$floor - floor_cal) < 0.5) && isTRUE(at$over_total) && isTRUE(at$over_floor) &&
         identical(half$over_total, FALSE) && identical(half$over_floor, FALSE),
       why = sprintf("the 703125fb check projects a floor of %.0f s (pretest %.0f s) and a total of %.0f s; half that check %.0f s and %.0f s",
                     at$floor, floor_cal, at$total, half$floor, half$total))
})
ctrl_failed <- names(ctrl)[!vapply(ctrl, function(x) isTRUE(x$ok), logical(1))]
gate(sprintf("positive controls: %d detectors each flag a planted defect (%s)", length(ctrl), paste(names(ctrl), collapse = ", ")),
     !length(ctrl_failed),
     unlist(lapply(names(ctrl), function(n) sprintf("%s %s: %s", if (n %in% ctrl_failed) "FAILED" else "ok", n, ctrl[[n]]$why))))

## ---- Tarballs of this version already at the hand-off ------------------------
# The file at the hand-off path is what goes to CRAN and what the dependent
# packages were told to test against; between 2026-11-13 and 2026-11-14 it
# was rewritten in place three times under SIEVER and SMSD. Full mode fails
# at once (section 6a) when a tarball of this version is already at the
# hand-off directory or in the package root, unless --replace-stale is
# given, and runs on so the record shows everything else; quick mode lists
# them. A fresh build can never match an earlier one byte for byte (R CMD
# build stamps the time into DESCRIPTION), so an existing file is always one
# the run would replace. Development builds for the dependents come from
# tools/dev_build.R, under names of their own.
.gate_stale_tarballs <- function(dirs, name) {
  p <- unique(file.path(normalizePath(dirs[dir.exists(dirs)]), name))
  p[file.exists(p)]
}
.gate_describe_file <- function(p) sprintf("%s (%.1f KB, modified %s, sha256 %s)", p, file.info(p)$size / 1024,
                                           format(file.info(p)$mtime, "%Y-%m-%d %H:%M"), unname(tools::sha256sum(p)))
stale_at_start <- .gate_stale_tarballs(c(dest, root), tarball_name)

if (quick) {
  if (length(stale_at_start))
    cat(sprintf("[INFO] already present (full mode fails on it unless --replace-stale is given): %s\n",
                .gate_describe_file(stale_at_start)), sep = "")
  ok_all <- all(vapply(results, `[[`, logical(1), "ok"))
  cat(sprintf("\nquick mode: %d gates, %d failed, %.0f s\n", length(results),
              sum(!vapply(results, `[[`, logical(1), "ok")),
              as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  if (!ok_all) quit(status = 1L)
  cat("QUICK GATES PASSED (not a release: run full mode before any tarball leaves)\n")
  quit(status = 0L)
}

## ---- 6a. Awake, and a version and date that can be released ------------
if (Sys.info()[["sysname"]] == "Darwin") {
  assertions <- tryCatch(system2("pmset", c("-g", "assertions"), stdout = TRUE), error = function(e) character())
  awake <- any(grepl("caffeinate", assertions, fixed = TRUE, useBytes = TRUE))
  gate("machine kept awake for the timing gates (caffeinate assertion present)", awake,
       if (!awake) "start with: caffeinate -i Rscript tools/release_gate.R")
  if (!awake) { cat("\nGATE FAILED before the check ran; nothing recorded\n"); quit(status = 1L) }
}
gate(if (length(stale_at_start) && replace_stale)
       sprintf("%d stale %s reported (--replace-stale given)", length(stale_at_start), tarball_name)
     else sprintf("no %s already at the hand-off directory or in the package root", tarball_name),
     !length(stale_at_start) || replace_stale,
     c(if (length(stale_at_start)) .gate_describe_file(stale_at_start),
       if (length(stale_at_start) && replace_stale)
         c("none is deleted; on a pass the certified tarball overwrites the one at the hand-off path, if any,",
           "and any other stays where it is, for you to delete")
       else if (length(stale_at_start))
         c("an uncertified file with the certified file's name invites the wrong upload: delete or move it, or rerun with --replace-stale",
           "(archive a superseded CRAN artifact before replacing it)")))
cran_ver <- tryCatch({
  ap <- utils::available.packages(repos = "https://cloud.r-project.org")
  if (pkg[, "Package"] %in% rownames(ap)) unname(ap[pkg[, "Package"], "Version"]) else NA_character_
}, error = function(e) NA_character_)
this_ver <- unname(pkg[, "Version"])
ver_ok <- !is.na(cran_ver) && package_version(this_ver) > package_version(cran_ver)
gate(sprintf("version %s is above the version on CRAN (%s)", this_ver,
             if (is.na(cran_ver)) "lookup failed" else cran_ver), ver_ok,
     if (is.na(cran_ver)) "CRAN could not be reached (the URL checks need the network too)"
     else if (!ver_ok) "bump Version in DESCRIPTION and open its NEWS heading before certifying")
desc_date <- tryCatch(as.Date(unname(pkg[, "Date"])), error = function(e) as.Date(NA))
date_age <- as.numeric(Sys.Date() - desc_date)
gate(sprintf("DESCRIPTION Date %s is the certification day (within three days)", format(desc_date)),
     !is.na(date_age) && date_age >= 0 && date_age <= 3,
     if (is.na(date_age)) "DESCRIPTION has no parseable Date field"
     else "date the release the day it is certified; CRAN notes a Date over a month old")

## ---- 6. Clean tree; the tarball is the committed tree ------------------
dirty <- system2("git", c("status", "--porcelain"), stdout = TRUE)
gate("working tree clean (the gate certifies a commit, not a draft)", length(dirty) == 0L, dirty)
sha <- system2("git", c("rev-parse", "--short", "HEAD"), stdout = TRUE)
src <- file.path(work, "src"); unlink(src, recursive = TRUE); dir.create(file.path(src, pkg[, "Package"]), recursive = TRUE)
system(sprintf("git archive HEAD | tar -x -C %s", shQuote(file.path(src, pkg[, "Package"]))))
# The vignettes go to CRAN as prebuilt HTML (2026-11-13). R CMD check
# rendered every precomputed twin through rmarkdown and pandoc on CRAN's
# clock (18 s of the 1.1.0 candidate's local check) for HTML that is the
# same at every check. In this copy of the commit only, the commit's own
# tools/asis_vignettes.R renders each twin by the call R CMD build makes,
# writes an R.rsp::asis entry beside the HTML, and removes the twin and its
# figures, so the check copies the HTML instead of rendering it. Each
# vignette's executable source, its .Rmd.orig, stays and ships beside the
# .asis entry. The rendered output exists only here and in the tarball; it
# is never committed. The twins and the checksums of their .Rmd.orig
# sources as committed are read first, for the vignette gate in section 10.
build_dir <- file.path(src, pkg[, "Package"])
archived_twins <- sort(sub("[.]Rmd$", "", list.files(file.path(build_dir, "vignettes"), pattern = "[.]Rmd$")))
orig_md5 <- tools::md5sum(file.path(build_dir, "vignettes", paste0(archived_twins, ".Rmd.orig")))
names(orig_md5) <- archived_twins
asis <- tryCatch({
  source(file.path(build_dir, "tools", "asis_vignettes.R"), local = TRUE)
  asis_vignettes(build_dir)
}, error = function(e) e)
asis_ok <- !inherits(asis, "error") && length(asis$name) == length(archived_twins) &&
  setequal(asis$name, archived_twins)
vig_render_s <- if (inherits(asis, "error")) NA_real_ else attr(asis, "render_s")
gate(sprintf("vignettes prebuilt as HTML in the build copy (%s of %d twins rendered in %s s, off the check's clock)",
             if (inherits(asis, "error")) 0L else nrow(asis), length(archived_twins),
             if (is.na(vig_render_s)) "?" else sprintf("%.0f", vig_render_s)),
     asis_ok && length(archived_twins) > 0L,
     if (inherits(asis, "error")) conditionMessage(asis)
     else if (!asis_ok) c("twins without HTML:", paste0("  ", setdiff(archived_twins, asis$name))))
old <- setwd(src)
bld <- run("R", c("CMD", "build", pkg[, "Package"]), file.path(work, "build.log"))
setwd(old)
tarball <- file.path(src, tarball_name)
gate(sprintf("tarball built from commit %s", sha), bld$status == 0L && file.exists(tarball),
     if (bld$status != 0L) tail(bld$lines, 5))

## ---- 7. R CMD check --as-cran, manuals on, remote URL checks on --------
chk_log <- file.path(work, "check.log")
old <- setwd(src)
Sys.setenv(`_R_CHECK_CRAN_INCOMING_REMOTE_` = "true", `_R_CHECK_FORCE_SUGGESTS_` = "true",
           `_R_CHECK_TIMINGS_` = "0")   # every step reports [cpu/elapsed], not only the slow ones
# The check runs the CRAN path whatever the shell that starts the gate
# has set: NOT_CRAN unset, so tests/testthat.R runs the curated subset
# and skip_on_cran() skips; at most two testthat workers, the most CRAN's
# policy lets a check use (the subset itself runs in one process); and
# one thread for the native libraries. These hold
# for the check's process only, and the caller's values are back when it
# returns.
check_env <- c(NOT_CRAN = NA_character_, TESTTHAT_CPUS = "2",
               OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
               MKL_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1")
chk_t0 <- Sys.time()
chk <- withr::with_envvar(check_env,
                          run("R", c("CMD", "check", "--as-cran", tarball_name), chk_log))
chk_wall_s <- as.numeric(difftime(Sys.time(), chk_t0, units = "secs"))
setwd(old)
log <- chk$lines
status_line <- grep("^Status:", log, value = TRUE)
# The NOTEs the gate allows, and nothing else:
# - the CRAN incoming feasibility NOTE when every item in it is one of
#   these: the Maintainer line, "New submission", possibly misspelled
#   words in DESCRIPTION, and invalid URLs whose 5xx status a re-probe
#   clears (a gateway timeout at the checker's moment is environmental).
#   Any other item fails the gate and is named: a URL that stays bad (one
#   or several), an invalid DOI, a Title not in title case, the size of
#   the tarball, days since the last update, a URL check that could not
#   run. R prints the items as paragraphs separated by blank lines
#   (tools:::format.check_package_CRAN_incoming), and the gate reads
#   them against this list, so an item it does not know fails rather
#   than passes;
# - the HTML manual NOTE whose whole text is R's message that the local
#   HTML Tidy is too old to validate the pages (win-builder validates
#   them).
# A NOTE's result is on its step's line or, when the step printed
# output first (a network warning, the tests' Running lines), on a line
# of its own, which belongs to the nearest step above it. Each NOTE is
# listed with its text under the gate's line, so a PASS shows what it
# tolerated.
.gate_reprobe <- function(u) {
  tryCatch(system2("curl", c("-s", "-o", "/dev/null", "-w", "%{http_code}", "-L",
                             "--max-time", "20", shQuote(u)), stdout = TRUE),
           error = function(e) "000")
}
.gate_feasibility_items <- function(text, reprobe = .gate_reprobe) {
  blank <- !nzchar(trimws(text))
  pars <- split(text[!blank], cumsum(blank)[!blank])
  flags <- character(); info <- character()
  for (p in pars) {
    first <- trimws(p[1L]); rest <- p[-1L]
    if (length(p) == 1L && (grepl("^Maintainer: ", first) || first == "New submission")) next
    if (first == "Possibly misspelled words in DESCRIPTION:" &&
        all(grepl("^\\s+\\S.* \\([0-9:, ]+\\)$", rest))) next
    if (grepl("^Found the following \\(possibly\\) invalid URLs?:$", first)) {
      # Entry lines are URL, From, Status, and Message, and continuations
      # of From and Message indented six spaces or more; R appends other
      # URL complaints to the same paragraph (a CRAN URL not in canonical
      # form, a publisher URL that should be a DOI), and those are flagged.
      odd <- rest[!grepl("^\\s*(URL|From|Status|Message): |^\\s{6,}\\S", rest)]
      flags <- c(flags, sprintf("in the URL block: %s", trimws(odd)))
      at <- grep("^\\s*URL: ", rest)
      ends <- c(at[-1L] - 1L, length(rest))
      for (k in seq_along(at)) {
        entry <- rest[at[k]:ends[k]]
        url <- sub("^\\s*URL: ", "", entry[1L])
        st <- sub("^\\s*Status: ", "", grep("^\\s*Status: ", entry, value = TRUE)[1L])
        if (is.na(st)) st <- "?"
        if (grepl("^5", st) && identical(reprobe(sub(" \\(moved to .*\\)$", "", url)), "200")) {
          info <- c(info, sprintf("transient %s at check time; a re-probe returned 200: %s", st, url))
        } else {
          flags <- c(flags, sprintf("invalid URL: %s [%s]", url, st))
        }
      }
      next
    }
    flags <- c(flags, first)
  }
  list(flags = flags, info = info)
}
.gate_old_tidy_only <- function(text) {
  s <- paste(trimws(text[nzchar(trimws(text))]), collapse = " ")
  grepl(paste0("^Skipping checking HTML validation: '[^']*' doesn't look like recent enough HTML Tidy\\.",
               "( Please obtain a recent version of HTML Tidy by downloading a binary release",
               " or compiling the source code from <https://www\\.html-tidy\\.org/>\\.)?$"), s)
}
.gate_check_notes <- function(log, reprobe = .gate_reprobe) {
  step_at <- grep("^\\* ", log)
  stop_at <- c(step_at, grep("^Status:", log), length(log) + 1L)
  note_at <- grep("^\\* .*\\.\\.\\..* NOTE$|^\\s*(\\[[^]]*\\] )?NOTE$", log)
  lapply(note_at, function(r) {
    above <- step_at[step_at <= r]
    h <- if (length(above)) max(above) else r
    e <- min(stop_at[stop_at > r]) - 1L
    text <- if (e > r) log[(r + 1L):e] else character()
    step <- if (h == r) trimws(log[r]) else
      paste(sub("^(\\* .*?\\.\\.\\.).*$", "\\1", log[h], perl = TRUE), trimws(log[r]))
    pre <- if (h < r) c(sub("^\\* .*?\\.\\.\\.", "", log[h], perl = TRUE),
                        if (r - h > 1L) log[(h + 1L):(r - 1L)])
    if (grepl("CRAN incoming feasibility", step)) {
      fi <- .gate_feasibility_items(text, reprobe)
      allowed <- length(fi$flags) == 0L
    } else if (grepl("HTML version of manual", step)) {
      allowed <- .gate_old_tidy_only(text)
      fi <- list(flags = if (!allowed) "its text is more than the message that the local HTML Tidy is too old",
                 info = character())
    } else {
      allowed <- FALSE
      fi <- list(flags = character(), info = character())
    }
    shown <- trimws(c(pre, text))
    list(step = step, text = shown[nzchar(shown)], flags = fi$flags, info = fi$info,
         allowed = allowed)
  })
}
.gate_note_lines <- function(n) {
  c(paste0(if (n$allowed) "allowed: " else "not allowed: ", n$step),
    if (length(n$flags)) paste0("  flagged: ", n$flags),
    if (length(n$info)) paste0("  ", n$info),
    if (length(n$text)) paste0("  ", n$text))
}
notes <- .gate_check_notes(log)
note_ok <- vapply(notes, `[[`, logical(1), "allowed")
# A NOTE that the Status line counts but the log does not show as a
# result fails the gate rather than passing unlisted.
status_notes <- regmatches(status_line, regexpr("[0-9]+ NOTE", status_line))
status_notes <- if (length(status_notes)) as.integer(sub(" NOTE", "", status_notes[1L])) else 0L
unclassified <- status_notes - length(notes)
# A check that did not run to its end passes neither gate below. Given a
# tarball that does not exist (the build failed), R CMD check warns,
# prints no Status line, and exits 0.
chk_done <- length(status_line) == 1L
gate("R CMD check --as-cran (manuals built): no ERROR, no WARNING",
     chk$status == 0L && chk_done && !any(grepl("ERROR|WARNING", status_line)),
     c(if (!chk_done) c("the check printed no Status line: it did not run to its end",
                        head(log[nzchar(trimws(log))], 5L)),
       status_line, grep("WARNING$|ERROR$", log, value = TRUE)))
gate(paste("R CMD check: no NOTE beyond incoming feasibility (new submission, spelling in",
           "DESCRIPTION, a URL server error that a re-probe clears) and HTML validation",
           "skipped by an old local HTML Tidy"),
     chk_done && all(note_ok) && unclassified <= 0L,
     c(if (!chk_done) "the check printed no Status line, so its NOTEs cannot be counted",
       unlist(lapply(notes[!note_ok], .gate_note_lines)),
       if (unclassified > 0L)
         sprintf("not allowed: the Status line counts %d NOTE(s) and the log shows %d to classify",
                 status_notes, length(notes)),
       unlist(lapply(notes[note_ok], .gate_note_lines))))
gate("PDF manual built by check", any(grepl("checking PDF version of manual \\.\\.\\..*OK", log)),
     grep("PDF version of manual", log, value = TRUE))
gate("HTML manual: no math rendering problems", !any(grepl("math rendering problems", log)),
     grep("math rendering", log, value = TRUE))
if (any(grepl("recent enough HTML Tidy", log)))
  cat("        note: HTML validation skipped locally (old tidy); win-builder validates it\n")
timings <- grep("\\[[0-9]+s/[0-9]+s\\]", log, value = TRUE)
for (t in grep("checking (examples|re-building of vignette outputs|PDF version)|Running .testthat", timings, value = TRUE))
  cat("        ", t, "\n", sep = "")

# Per-page example times from the check just run. CRAN flags a page over
# 5 s on its farm, and win-builder has run this package's examples at about
# three times the local elapsed time (122 s against 42 s, 2026-08-17), so
# the local ceiling is 1.7 s per page. Pages over 1 s are listed so a creep
# toward the ceiling is visible before it arrives.
ex_timings <- file.path(src, paste0(pkg[, "Package"], ".Rcheck"), paste0(pkg[, "Package"], "-Ex.timings"))
if (file.exists(ex_timings)) {
  et <- utils::read.table(ex_timings, header = TRUE)
  et$max <- pmax(et$user + et$system, et$elapsed)
  et <- et[order(-et$max), ]
  over <- et[et$max > 1.7, ]
  near <- et[et$max > 1.0 & et$max <= 1.7, ]
  gate(sprintf("no example page over 1.7 s locally (slowest %.2f s: %s; %d pages, %.0f s in all)",
               et$max[1], et$name[1], nrow(et), sum(et$elapsed)),
       nrow(over) == 0L,
       c(sprintf("%s %.2f s", over$name, over$max),
         if (nrow(near)) c("pages over 1 s:", sprintf("%s %.2f s", near$name, near$max))))
} else {
  gate("example timings file present after the check", FALSE, ex_timings)
}

## ---- 7a. The check as a whole, against the accepted release ---------------
# One [cpu/elapsed] pair per timed step, parsed from the log. The overall
# checktime is what the 2026-09-05 pretest returned the package for, so
# it is gated here in full, not only per page: the local wall clock less
# the network-bound feasibility step against the accepted release's, and
# each timed step against its own budget.
#
# Each timed step is compared with its accepted baseline and noise
# allowance. Whole-check and test budgets are set above. Network-bound
# feasibility time is kept separate from computation time.
step_secs <- function(pattern) {
  l <- grep(pattern, timings, value = TRUE)[1]
  if (is.na(l)) return(c(cpu = NA_real_, elapsed = NA_real_))
  v <- as.numeric(strsplit(gsub("[][s]", "", regmatches(l, regexpr("\\[[0-9]+s/[0-9]+s\\]", l))), "/")[[1]])
  c(cpu = v[1], elapsed = v[2])
}
steps <- rbind(examples  = step_secs("checking examples"),
               tests     = step_secs("Running .testthat[.]R."),
               vignettes = step_secs("re-building of vignette outputs"),
               manual    = step_secs("PDF version of manual"))
feas_s <- step_secs("CRAN incoming feasibility")[["elapsed"]]
# When the feasibility step prints a warning (a CRAN address that timed
# out), its [cpu/elapsed] pair is on a later line of the log, not on the
# line that names the step, and the pattern above finds nothing; the time
# of the network was then counted as the check's (2026-11-11: 226 s of a
# 60 s timeout and its retries, a whole check of 460 s against 234 s).
# Read the first pair between the step's line and the next step's.
if (is.na(feas_s)) {
  f0 <- grep("checking CRAN incoming feasibility", log)[1]
  if (!is.na(f0)) {
    nxt <- grep("^\\* checking ", log)
    f1 <- min(c(nxt[nxt > f0], length(log) + 1L)) - 1L
    pair <- regmatches(log[f0:f1], regexpr("\\[[0-9]+s/[0-9]+s\\]", log[f0:f1]))[1]
    if (!is.na(pair)) feas_s <- as.numeric(strsplit(gsub("[][s]", "", pair), "/")[[1]])[2]
  }
}
check_total_s <- chk_wall_s - (if (is.na(feas_s)) 0 else feas_s)
ratio <- check_total_s / accepted$total_s
ratio_max <- gate_limit_s /
  accepted$farm_total_s[["r-devel-linux-x86_64-fedora-gcc"]]
budget <- accepted$step_s * ratio_max + 5
budget[["tests"]] <- max(budget[["tests"]], tests_budget_s)
over_step <- !is.na(steps[, "elapsed"]) & steps[, "elapsed"] > budget[rownames(steps)]
# Since 2026-11-13 the vignette step copies the prebuilt HTML (section 6)
# instead of rendering the twins: 21 R processes that each load R.rsp and
# copy one file, 6 s of CPU and 8 to 9 s elapsed against 14 to 15 s of CPU
# and 18 to 19 s elapsed for the render (two checks of each, back to back,
# on 2026-11-13; the accepted 15 s is 1.0.0's render). Its CPU time is held
# to vignettes_cpu_s, so a check that renders again fails here as well as
# in the vignette gate of section 10. CPU time rather than elapsed, because
# the load of a shared machine moves elapsed time far more than CPU time.
vignettes_cpu_s <- 10
vignettes_cpu_over <- isTRUE(steps["vignettes", "cpu"] > vignettes_cpu_s)
projected <- round(accepted$farm_total_s * ratio)
# The r-devel-windows floor and the two lines of information about it
# (calibration block; information, never a failure).
win_elapsed <- c(r_code = step_secs("checking R code for possible problems")[["elapsed"]],
                 rd = step_secs("checking Rd files")[["elapsed"]],
                 examples = steps["examples", "elapsed"], tests = steps["tests", "elapsed"],
                 vignettes = steps["vignettes", "elapsed"], pdf_manual = steps["manual", "elapsed"],
                 html_manual = step_secs("checking HTML version of manual")[["elapsed"]])
win <- tryCatch(.gate_windows_projection(win_elapsed, check_total_s),
                error = function(e) list(lines = paste("r-devel-windows projection not computed:", conditionMessage(e))))
time_lines <- c(
  sprintf("%-10s %4.0f s elapsed (%3.0f s CPU); accepted %3d s, budget %3.0f s%s",
          rownames(steps), steps[, "elapsed"], steps[, "cpu"], accepted$step_s[rownames(steps)],
          budget[rownames(steps)], ifelse(over_step, "  <- over budget", "")),
  sprintf("vignettes %.0f s CPU against %d s for copying the prebuilt HTML%s; the HTML was rendered by the gate in %s s",
          steps["vignettes", "cpu"], vignettes_cpu_s, if (vignettes_cpu_over) "  <- over budget" else "",
          if (is.na(vig_render_s)) "?" else sprintf("%.0f", vig_render_s)),
  sprintf("whole check %.0f s locally (feasibility step %.0f s excluded); accepted %d s on %s",
          check_total_s, if (is.na(feas_s)) 0 else feas_s, accepted$total_s, accepted$measured),
  sprintf("projected farm totals at this ratio: %s",
          paste(sprintf("%s %d s", names(projected), projected), collapse = "; ")),
  sprintf("incoming pretest line %d s; slowest farm flavor for 1.0.0 %d s (%s)",
          incoming_limit_s, max(accepted$farm_total_s), names(which.max(accepted$farm_total_s))),
  if (projected[["r-devel-linux-x86_64-fedora-gcc"]] > incoming_limit_s)
    sprintf("projected r-devel Linux total %d s is over the %d s line the 2026-09-05 pretest enforced: trim before submitting",
            projected[["r-devel-linux-x86_64-fedora-gcc"]], incoming_limit_s),
  win$lines)
gate(sprintf("check time within ratio %.2f of the accepted release (%.0f s against %d s, ratio %.2f)",
             ratio_max, check_total_s, accepted$total_s, ratio),
     ratio <= ratio_max && !any(over_step) && !vignettes_cpu_over, time_lines)

# A run that slept or shared the machine shows elapsed time far above CPU
# time; the numbers above are then meaningless and the run is repeated.
if (file.exists(ex_timings)) {
  cpu_ex <- sum(et$user + et$system); el_ex <- sum(et$elapsed)
  gate(sprintf("examples ran awake and unshared (%.0f s elapsed against %.0f s CPU)", el_ex, cpu_ex),
       el_ex <= 1.5 * cpu_ex + 5,
       if (el_ex > 1.5 * cpu_ex + 5) "elapsed far above CPU: the machine slept or was busy; rerun the gate")
}

# The CRAN path of the test suite is the curated subset in tests/testthat.R.
# Its summary line is read from the check's own run, so a filter that
# silently matched fewer files (a renamed test) shows up as a count below
# the floor, and its time is held to a budget. The floor was 2000 while the
# subset ran 110 files. On 2026-11-13 the subset was cut to 16 files and 266
# expectations, and the anchors of the closed-form sample size planners
# then took it to 40 files and 527 expectations, the promised detectors 102
# of them. The floor is that count, not a round number under it: a
# top-level skip_on_cran() placed above an anchor still reports one skip,
# so only the count shows a file that lost its anchor. An added expectation
# passes; a lost one fails until the floor is lowered on purpose.
cran_tests_floor <- 527L
rout <- file.path(src, paste0(pkg[, "Package"], ".Rcheck"), "tests", "testthat.Rout")
if (!file.exists(rout)) rout <- paste0(rout, ".fail")
rout_summary <- if (file.exists(rout)) grep("^\\[ FAIL", readLines(rout, warn = FALSE), value = TRUE) else character()
counts <- if (length(rout_summary)) as.integer(regmatches(rout_summary[length(rout_summary)],
                                                          gregexpr("[0-9]+", rout_summary[length(rout_summary)]))[[1]]) else integer()
tests_s <- steps["tests", "elapsed"]
tests_ok <- length(counts) == 4L && counts[1] == 0L && counts[2] == 0L && counts[4] >= cran_tests_floor &&
  !is.na(tests_s) && tests_s <= tests_budget_s
gate(sprintf("CRAN-path tests (curated subset): %s passed, %s failed, %s warnings, %s skipped in %s s",
             if (length(counts) == 4L) counts[4] else "?", if (length(counts) == 4L) counts[1] else "?",
             if (length(counts) == 4L) counts[2] else "?", if (length(counts) == 4L) counts[3] else "?",
             if (is.na(tests_s)) "?" else sprintf("%.0f", tests_s)), tests_ok,
     c(if (length(counts) != 4L) "no testthat summary line in the check's testthat.Rout",
       if (length(counts) == 4L && counts[4] < cran_tests_floor)
         sprintf("under %d expectations: the cran_files filter is losing files", cran_tests_floor),
       if (!is.na(tests_s) && tests_s > tests_budget_s) sprintf("over the %.0f s budget for the CRAN path", tests_budget_s)))

## ---- 7b. References against their registered metadata (information) ----
# tools/doi_audit.py compares every \\doi{} reference in the roxygen
# sources with the metadata Crossref holds for it (first author, year,
# title, volume, pages) and reports the references that disagree. This is
# printed and recorded, and it never fails the gate: a disagreement is a
# prompt to open the article (registries hold errors too), and the audit
# needs the network. Set DOI_AUDIT_CONTACT to send Crossref a contact
# address with the requests; none is sent otherwise.
doi_lines <- tryCatch(
  suppressWarnings(system2("python3", c("tools/doi_audit.py", "--quiet"),
                           stdout = TRUE, stderr = FALSE)),
  error = function(e) character())
doi_lines <- grep("^DOI audit:", doi_lines, value = TRUE)
if (length(doi_lines) == 0L) doi_lines <- "DOI audit: not run (python3 or tools/doi_audit.py unavailable)"
cat(sprintf("[INFO] %s\n", doi_lines), sep = "")

## ---- 8. The full local test suite, Monte Carlo blocks included ----------
# The suite passes when it ran: at least one expectation passed, every
# total is a finite count, nothing failed, errored, warned, or was
# skipped, and every test file ran at least one test. A skip is a test
# that did not run on the machine certifying the release (a Suggests
# package missing, say), and an empty or unreadable result is not a
# pass. A skip at the top level of a file (eleven files open with
# skip_if_not_installed() for lavaan or boot) leaves no row in the
# result at all, so the files are counted as well as the tests.
.gate_full_suite_ok <- function(x) {
  required <- c("passed", "failed", "error", "warning", "skipped")
  if (!is.data.frame(x) || nrow(x) == 0L || !all(required %in% names(x))) return(FALSE)
  if (!all(vapply(x[required], function(v) is.numeric(v) || is.logical(v), logical(1))))
    return(FALSE)
  totals <- vapply(x[required], function(v) as.numeric(sum(v)), numeric(1))
  all(is.finite(totals)) && totals[["passed"]] > 0 &&
    all(totals[c("failed", "error", "warning", "skipped")] == 0)
}
# The test files (testthat's own pattern) that have no row in the result.
.gate_suite_files_missing <- function(x, files) {
  ran <- if (is.data.frame(x) && "file" %in% names(x)) unique(as.character(x$file)) else character()
  setdiff(files, ran)
}
suite_files <- list.files(file.path("tests", "testthat"), pattern = "^test.*[.][rR]$")
Sys.setenv(NOT_CRAN = "true")
tst <- tryCatch(as.data.frame(testthat::test_local(".", reporter = "silent", stop_on_failure = FALSE)),
                error = function(e) e)
tst_err <- if (inherits(tst, "error")) conditionMessage(tst)
if (inherits(tst, "error")) tst <- NULL
tst_missing <- if (is.null(tst)) character() else .gate_suite_files_missing(tst, suite_files)
tst_ok <- .gate_full_suite_ok(tst) && length(suite_files) > 0L && length(tst_missing) == 0L
tst_total <- function(col) {
  v <- tryCatch(sum(tst[[col]]), error = function(e) NA)
  if (is.null(tst) || length(v) != 1L || is.na(v)) "?" else v
}
tst_bad <- tryCatch({
  b <- tst$failed > 0 | tst$error | tst$warning > 0 | tst$skipped
  b[is.na(b)] <- TRUE
  # Code outside test_that() (a stop() at the top of a file) has no test name.
  unique(sprintf("%s: %s", tst$file[b],
                 ifelse(is.na(tst$test[b]), "(code outside test_that())", tst$test[b])))
}, error = function(e) character())
gate(sprintf("full local test suite: %s passed, %s failed, %s errors, %s warnings, %s skipped; %d of %d files ran",
             tst_total("passed"), tst_total("failed"), tst_total("error"),
             tst_total("warning"), tst_total("skipped"),
             if (is.null(tst)) 0L else length(suite_files) - length(tst_missing), length(suite_files)),
     tst_ok,
     c(if (!is.null(tst_err)) paste("test_local() stopped:", tst_err),
       head(tst_bad, 20L),
       if (length(tst_bad) > 20L) sprintf("and %d more", length(tst_bad) - 20L),
       if (length(tst_missing))
         c("files that ran no test (a skip at the top of the file, or no test_that() in it):",
           paste0("  ", head(tst_missing, 20L)),
           if (length(tst_missing) > 20L) sprintf("  and %d more", length(tst_missing) - 20L)),
       if (!tst_ok && is.null(tst_err) && !length(tst_bad) && !length(tst_missing))
         "no expectation passed, no test file was found, or the result had no counts to read"))

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
# The harness skips a snippet whose oracle package is not installed and
# still exits 0, so its exit status alone would pass a run that compared
# nothing. The gate reads the one summary line the harness prints and
# requires every snippet to have run: at least one ok, none skipped, and
# none failed.
.gate_oracles_ok <- function(status, lines) {
  summary <- grep("^oracle checks: [0-9]+ ok, [0-9]+ skipped, [0-9]+ failed$",
                  lines, value = TRUE)
  if (!isTRUE(status == 0L) || length(summary) != 1L) return(FALSE)
  counts <- suppressWarnings(as.integer(regmatches(summary, gregexpr("[0-9]+", summary))[[1L]]))
  isTRUE(length(counts) == 3L && counts[1L] > 0L && counts[2L] == 0L && counts[3L] == 0L)
}
gate("oracle harness (live cross-package anchors, none skipped)",
     .gate_oracles_ok(orc$status, orc$lines),
     grep("oracle checks:|^Error|skipped \\(oracle not installed\\)", orc$lines, value = TRUE))

## ---- 9a. What DMAR's dependents rely on, against the tarball --------------
# tools/revdep_contract.R reads SIEVER's and SMSD's sources for every DMAR
# object they reach, compares each object's formals (or a dataset's
# columns) with its committed snapshot, and runs the canonical calls the
# snapshot records, against the DMAR just installed in the scratch library.
# Its positive control: a copy of the snapshot in which the first recorded
# formal of one function is renamed must exit 1 and name that function. A
# dependent whose checkout is missing makes the script exit 3, and this
# step fails: what that dependent relies on now was not read, and nothing
# in full mode is skipped.
.gate_revdep_summary <- function(lines) {
  s <- grep("^revdep contract: [0-9]+ of [0-9]+ dependents read \\(.*\\); [0-9]+ differences?$", lines, value = TRUE)
  if (length(s) != 1L) return(NULL)
  n <- as.integer(regmatches(sub("\\(.*\\)", "", s), gregexpr("[0-9]+", sub("\\(.*\\)", "", s)))[[1L]])
  list(read = n[1L], of = n[2L], differences = n[3L], line = s)
}
.gate_revdep_step <- function(lib, work, snapshot = file.path("tools", "revdep_contract.json")) {
  rc <- run("Rscript", c("tools/revdep_contract.R", paste0("--lib=", lib), paste0("--snapshot=", snapshot)),
            file.path(work, "revdep_contract.log"))
  s <- .gate_revdep_summary(rc$lines)
  planted <- file.path(work, "revdep_contract_planted.json")
  target <- tryCatch({
    snap <- jsonlite::read_json(snapshot, simplifyVector = FALSE)
    fns <- names(Filter(function(o) identical(o$kind, "function") && length(o$formals), snap$objects))
    t <- fns[1L]
    names(snap$objects[[t]]$formals)[1L] <- paste0(names(snap$objects[[t]]$formals)[1L], "_planted")
    jsonlite::write_json(snap, planted, pretty = TRUE, auto_unbox = TRUE, digits = NA, null = "null")
    t
  }, error = function(e) NA_character_)
  pc <- if (!is.na(target)) run("Rscript", c("tools/revdep_contract.R", paste0("--lib=", lib), paste0("--snapshot=", planted)),
                                file.path(work, "revdep_contract_planted.log"))
  ctrl_ok <- !is.na(target) && identical(as.integer(pc$status), 1L) &&
    any(grepl(sprintf("^\\[DIFF\\] %s \\(.*\\): formals differ", target), pc$lines))
  ok <- identical(as.integer(rc$status), 0L) && !is.null(s) && s$read == s$of && s$differences == 0L && ctrl_ok
  list(ok = ok,
       name = sprintf("contract with DMAR's dependents (tools/revdep_contract.R, against the tarball): %s; positive control %s",
                      if (is.null(s)) "no summary line" else sub("^revdep contract: ", "", s$line),
                      if (ctrl_ok) "flagged" else "NOT flagged"),
       detail = c(grep("^\\[(OK|DIFF|NOT RUN)\\]", rc$lines, value = TRUE),
                  if (is.null(s)) c(sprintf("exit status %s; the last lines:", rc$status), paste0("  ", tail(rc$lines, 5L))),
                  if (!is.null(s) && s$read < s$of) "a dependent was not read: in full mode that is a failure, not a skip",
                  sprintf("positive control: the snapshot with the first formal of %s renamed %s", target,
                          if (ctrl_ok) "was flagged (exit 1)" else sprintf("was NOT flagged (exit %s)", if (is.null(pc)) "?" else pc$status))))
}
rcs <- .gate_revdep_step(gate_lib, work)
gate(rcs$name, rcs$ok, rcs$detail)

## ---- 10. Nothing in the tarball or the history that should not be ------
# The pattern (.gate_banned, section 5a) is assembled from character codes
# so this file never carries the literal strings it polices; the sweep of
# the commit messages is section 5a's function.
banned <- .gate_banned
names_hit <- system(sprintf("tar -tzf %s | grep -icE '%s'", shQuote(tarball), banned), intern = TRUE)
content_hit <- system(sprintf("tar -xzOf %s 2>/dev/null | grep -icE '%s'", shQuote(tarball), banned), intern = TRUE)
log_hit <- .gate_log_residue(".")
# grep -c exits 1 on a zero count, and system(intern = TRUE) then tags the
# "0" with a status attribute; compare the numbers, not the objects.
counts <- suppressWarnings(as.integer(c(names_hit, content_hit, log_hit)))
gate("tarball and git history carry no authorship residue",
     length(counts) == 3L && !anyNA(counts) && sum(counts) == 0L,
     sprintf("file names %s, file content %s, commit messages %s", names_hit, content_hit, log_hit))
entries <- system(sprintf("tar -tzf %s", shQuote(tarball)), intern = TRUE)
surfaces <- unique(sub("/.*", "", sub("^[^/]+/", "", entries)))
gate("tarball holds only the standard surfaces", all(surfaces %in% c("DESCRIPTION", "NAMESPACE", "NEWS.md", "README.md", "R", "build", "data", "inst", "man", "tests", "vignettes", "LICENSE", "LICENSE.md")),
     setdiff(surfaces, c("DESCRIPTION", "NAMESPACE", "NEWS.md", "README.md", "R", "build", "data", "inst", "man", "tests", "vignettes", "LICENSE", "LICENSE.md")))
# The vignettes in the tarball are the commit's twins, prebuilt (section 6).
# vignettes/ holds, for each twin, one <name>.html.asis and one
# <name>.Rmd.orig and nothing else, so no twin or figure is left for the
# check to render and every vignette ships its executable source; each
# shipped .Rmd.orig is byte for byte the committed one; inst/doc holds
# the .asis file and <name>.html; each HTML is byte for byte the one
# rendered from its twin; each HTML still carries the record, which
# tools/precompute_vignettes.R writes into the twin, of the MD5 of the
# .Rmd.orig the twin was knitted from, and that MD5 is the committed
# .Rmd.orig's (tools/release_checks.R checks the record on the twins; this
# checks it on what ships); and the index R CMD build wrote, which
# vignette() and browseVignettes() read, lists every twin with its title
# and nothing else, so no .Rmd.orig was taken for a vignette.
pfx <- paste0(pkg[, "Package"], "/")
in_dir <- function(d) {
  hits <- grep(sprintf("^%s%s/.*[^/]$", pfx, d), entries, value = TRUE)  # files at any depth
  sort(substring(hits, nchar(pfx) + nchar(d) + 2L))
}
vig_diff <- function(where, got, want)
  c(if (length(x <- setdiff(got, want))) sprintf("%s holds %s, which is neither a prebuilt vignette nor its source", where, x),
    if (length(x <- setdiff(want, got))) sprintf("%s lacks %s", where, x))
want_asis <- paste0(archived_twins, ".html.asis")
want_orig <- paste0(archived_twins, ".Rmd.orig")
vig_bad <- c(vig_diff("vignettes/", in_dir("vignettes"), c(want_asis, want_orig)),
             vig_diff("inst/doc/", in_dir("inst/doc"), c(want_asis, paste0(archived_twins, ".html"))))
vig_ex <- file.path(work, "tarball_vignettes"); unlink(vig_ex, recursive = TRUE); dir.create(vig_ex)
html_in <- paste0(pfx, "inst/doc/", archived_twins, ".html")
orig_in <- paste0(pfx, "vignettes/", want_orig)
index_in <- paste0(pfx, "build/vignette.rds")
present <- intersect(c(html_in, orig_in, index_in), entries)
if (length(present)) utils::untar(tarball, files = present, exdir = vig_ex)
for (k in seq_along(archived_twins)) {
  v <- archived_twins[k]
  o <- file.path(vig_ex, orig_in[k])
  if (file.exists(o) && !identical(unname(tools::md5sum(o)), unname(orig_md5[[v]])))
    vig_bad <- c(vig_bad, sprintf("vignettes/%s.Rmd.orig is not the committed source", v))
  f <- file.path(vig_ex, html_in[k])
  if (!file.exists(f)) next   # named above
  if (!asis_ok || !identical(unname(tools::md5sum(f)), asis$html_md5[match(v, asis$name)]))
    vig_bad <- c(vig_bad, sprintf("inst/doc/%s.html is not the HTML rendered from its twin", v))
  rec <- sprintf("precomputed by tools/precompute_vignettes.R from %s.Rmd.orig, md5 %s", v, orig_md5[[v]])
  if (!any(grepl(rec, readLines(f, warn = FALSE), fixed = TRUE)))
    vig_bad <- c(vig_bad, sprintf("inst/doc/%s.html does not record the MD5 of %s.Rmd.orig as committed", v, v))
}
index_file <- file.path(vig_ex, index_in)
if (file.exists(index_file)) {
  vix <- readRDS(index_file)
  at <- match(want_asis, vix$File)
  if (anyNA(at) || nrow(vix) != length(want_asis) ||
      !identical(unname(vix$PDF[at]), paste0(archived_twins, ".html")) ||
      (asis_ok && !identical(unname(vix$Title[at]), asis$title[match(archived_twins, asis$name)])))
    vig_bad <- c(vig_bad, "build/vignette.rds does not list every twin's HTML with the twin's title")
} else {
  vig_bad <- c(vig_bad, "the tarball has no build/vignette.rds, so vignette() would list nothing")
}
gate(sprintf("tarball vignettes: the commit's %d twins as prebuilt HTML with their .Rmd.orig sources as committed (HTML rendered from each twin, MD5 record of each .Rmd.orig, every title indexed)",
             length(archived_twins)),
     length(archived_twins) > 0L && length(vig_bad) == 0L, vig_bad)

## ---- Record and hand-off ---------------------------------------------------
ok_all <- all(vapply(results, `[[`, logical(1), "ok"))
sha256 <- unname(tools::sha256sum(tarball))
elapsed <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
stamp <- format(Sys.time(), "%Y-%m-%d_%H%M")
record <- c(
  sprintf("# %s release gate, %s, commit %s", pkg[, "Package"], format(Sys.time(), "%Y-%m-%d %H:%M"), sha),
  "", sprintf("Result: %s (%d gates, %.1f minutes)", if (ok_all) "ALL GATES PASSED" else "FAILED",
              length(results), elapsed), "",
  sprintf("Tarball: %s", tarball_name), sprintf("sha256: %s", sha256),
  sprintf("Tarballs of this version present at the start: %s%s",
          if (length(stale_at_start)) paste(stale_at_start, collapse = ", ") else "none",
          if (length(stale_at_start) && replace_stale) " (--replace-stale given)" else ""), "",
  unlist(lapply(results, function(r) c(sprintf("- [%s] %s", if (r$ok) "PASS" else "FAIL", r$name),
                                       if (length(r$detail)) paste0("    ", r$detail)))),
  "", "References (information only):", paste0("    ", doi_lines),
  "", "Check timings:", paste0("    ", timings), "",
  "Check time against the accepted release:", paste0("    ", time_lines), "",
  paste0("Check status: ", status_line))
rec_path <- file.path(qc_dir, sprintf("%s_release_gate_%s_%s.md", pkg[, "Package"], stamp, sha))
writeLines(record, rec_path)
cat(sprintf("\nrecord: %s\n", rec_path))
# Delivery is verified, not assumed. A sidecar an earlier run left at the
# destination is removed first (emptied, where the folder allows no
# removal), so that a sidecar there afterward vouches only for a tarball
# this run checked; the copy must exist and carry the sha256 of the
# tarball the gate checked; and a new sidecar in the form shasum -a 256 -c
# reads (the hash, two spaces, the file name) goes beside it. The record
# is written first, so no tarball leaves without one, and is written
# again with the delivery's outcome: a failed delivery changes its Result
# line and fails the run although every gate passed.
.gate_deliver <- function(tarball, dest, sha256, replace_stale = FALSE) {
  delivered <- file.path(dest, basename(tarball))
  sidecar <- paste0(delivered, ".sha256")
  quietly <- function(expr) tryCatch(expr, warning = function(w) conditionMessage(w),
                                     error = function(e) conditionMessage(e))
  # A file at the hand-off path whose sha256 is not the checked tarball's is
  # overwritten only with --replace-stale (it may be the file that went to
  # CRAN, or the one the dependents tested); neither it nor its sidecar is
  # touched otherwise.
  if (file.exists(delivered)) {
    had <- unname(tools::sha256sum(delivered))
    if (!identical(had, sha256) && !isTRUE(replace_stale))
      return(list(ok = FALSE, detail = c(
        sprintf("%s is already there with sha256 %s, not the checked tarball's %s; it was not overwritten",
                delivered, had, sha256),
        "move it away (archive a superseded CRAN artifact before replacing it), or rerun with --replace-stale")))
  }
  if (file.exists(sidecar)) {
    cleared <- quietly(file.remove(sidecar))
    if (file.exists(sidecar))
      cleared <- quietly({ writeLines(character(), sidecar); isTRUE(file.size(sidecar) == 0) })
    if (!isTRUE(cleared))
      return(list(ok = FALSE, detail = c(
        sprintf("the sidecar %s from an earlier run could be neither removed nor emptied; it does not vouch for this run's tarball",
                sidecar),
        if (is.character(cleared)) cleared)))
  }
  copied <- quietly(file.copy(tarball, delivered, overwrite = TRUE))
  if (!isTRUE(copied) || !file.exists(delivered))
    return(list(ok = FALSE, detail = c(sprintf("the copy to %s failed", delivered),
                                       if (is.character(copied)) copied)))
  got <- unname(tools::sha256sum(delivered))
  if (!identical(got, sha256))
    return(list(ok = FALSE, detail = sprintf("%s has sha256 %s, not the checked tarball's %s",
                                             delivered, got, sha256)))
  line <- paste0(sha256, "  ", basename(tarball))
  wrote <- quietly({ writeLines(line, sidecar); TRUE })
  back <- if (isTRUE(wrote)) quietly(readLines(sidecar, warn = FALSE))
  if (!isTRUE(wrote) || !identical(back, line))
    return(list(ok = FALSE, detail = c(sprintf("the checksum sidecar %s was not written", sidecar),
                                       if (is.character(wrote)) wrote)))
  list(ok = TRUE, detail = c(sprintf("tarball: %s (sha256 matches the checked tarball)", delivered),
                             sprintf("sidecar: %s", sidecar)))
}
if (ok_all) {
  dlv <- .gate_deliver(tarball, dest, sha256, replace_stale)
  if (!dlv$ok) {
    dlv$detail <- c(dlv$detail, sprintf("any %s at %s is not the certified tarball; check it against the sha256 on the record",
                                        tarball_name, dest))
    record[grep("^Result: ", record)[1L]] <- sprintf(
      "Result: ALL GATES PASSED, DELIVERY FAILED; the hand-off is not certified (%d gates, %.1f minutes)",
      length(results), elapsed)
  }
  writeLines(c(record, "", if (dlv$ok) "Delivery: verified" else "Delivery: failed; the hand-off is not certified",
               paste0("    ", dlv$detail)), rec_path)
  if (!dlv$ok) {
    cat(sprintf("\nDELIVERY FAILED after all %d gates passed; the hand-off is not certified\n",
                length(results)))
    for (d in dlv$detail) cat("        ", d, "\n", sep = "")
    quit(status = 1L)
  }
  cat(sprintf("tarball: %s\nsha256:  %s\nsidecar: %s.sha256\n", file.path(dest, tarball_name), sha256,
              file.path(dest, tarball_name)))
  cat(sprintf("\nALL GATES PASSED (%d gates, %.1f min)\n", length(results), elapsed))
  cat("next: (1) R CMD INSTALL the certified tarball so the console matches the\n",
      "submission; (2) a win-builder r-devel run closes the HTML Tidy gap this\n",
      "machine cannot check; (3) upload exactly the file above, after\n",
      sprintf("shasum -a 256 -c %s.sha256, run in its folder, passes and its\n", tarball_name),
      "sha256 matches this record.\n", sep = "")
  quit(status = 0L)
}
cat(sprintf("\nGATE FAILED: %d of %d gates failed; no tarball handed off\n",
            sum(!vapply(results, `[[`, logical(1), "ok")), length(results)))
quit(status = 1L)
