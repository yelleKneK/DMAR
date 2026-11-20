# Pre-release quality control for DMAR.
#
# Runs every invariant the project conventions state plus the ones learned during the
# 1.0.0 check time work, mechanically, and prints one PASS/FAIL line per
# check with the offending files when a check fails. Exit status is nonzero
# on any failure, so this can gate a release script.
#
# Usage, from the package root:
#   Rscript tools/release_checks.R
#
# This file is not shipped: tools/ is excluded by .Rbuildignore.

root <- normalizePath(".")
if (!file.exists(file.path(root, "DESCRIPTION"))) {
  stop("Run from the package root.", call. = FALSE)
}

failures <- 0L
check <- function(name, ok, detail = character()) {
  status <- if (ok) "PASS" else "FAIL"
  cat(sprintf("[%s] %s\n", status, name))
  if (!ok) {
    failures <<- failures + 1L
    for (d in detail) cat("       ", d, "\n", sep = "")
  }
  invisible(ok)
}

r_files    <- list.files("R", pattern = "[.]R$", full.names = TRUE)
rd_files   <- list.files("man", pattern = "[.]Rd$", full.names = TRUE)
vig_rmd    <- list.files("vignettes", pattern = "[.]Rmd$", full.names = TRUE)
vig_orig   <- list.files("vignettes", pattern = "[.]Rmd[.]orig$", full.names = TRUE)
test_files <- list.files("tests/testthat", pattern = "^test-.*[.]R$", full.names = TRUE)

grep_files <- function(pattern, files, ...) {
  hits <- character()
  for (f in files) {
    lines <- readLines(f, warn = FALSE)
    m <- grep(pattern, lines, ...)
    if (length(m)) hits <- c(hits, sprintf("%s:%d", f, m))
  }
  hits
}

## ---- Check time invariants (the 1.0.0 rejection class) ----

# One \donttest{} anywhere makes --as-cran run the whole example corpus
# twice; the package ships none, and no \dontrun{} either (the house idiom
# for "shown but not run" is a plain comment).
h <- grep_files("\\\\donttest|\\\\dontrun", c(r_files, rd_files))
check("no \\donttest or \\dontrun anywhere", length(h) == 0, h)

# No commented-out code in any @examples block (the 2026-09-04 CRAN return:
# "Some code lines in examples are commented out. Please never do that.").
# Every contiguous window of comment lines inside an @examples block is
# parsed with the markers stripped; a window that parses to a call or an
# assignment is commented-out code. Prose does not parse and is never
# flagged. The same detector runs on the built Rd files in
# tests/testthat/test-rd_hygiene.R; this one reads the roxygen sources so it
# fires before documentation is regenerated.
is_code <- function(txt) {
  parsed <- tryCatch(parse(text = txt, keep.source = FALSE), error = function(e) NULL)
  if (is.null(parsed) || !length(parsed)) return(FALSE)
  any(vapply(parsed, function(e) is.call(e) && !identical(e[[1L]], as.name("~")), logical(1)))
}
commented_hits <- character()
for (f in r_files) {
  lines <- readLines(f, warn = FALSE)
  for (s in grep("^#' @examples", lines)) {
    i <- s + 1L; ex <- character(); idx <- integer()
    while (i <= length(lines) && grepl("^#'", lines[i]) && !grepl("^#' @[a-zA-Z]", lines[i])) {
      ex <- c(ex, sub("^#' ?", "", lines[i])); idx <- c(idx, i); i <- i + 1L
    }
    cm <- grep("^\\s*#", ex)
    if (!length(cm)) next
    for (run in split(cm, cumsum(c(1L, diff(cm) != 1L)))) {
      txt <- sub("^\\s*#+\\s?", "", ex[run]); n <- length(run); flagged <- integer()
      for (a in seq_len(n)) for (b in a:n) {
        w <- paste(txt[a:b], collapse = "\n")
        if (nzchar(trimws(w)) && is_code(w)) flagged <- c(flagged, run[a:b])
      }
      if (length(flagged)) commented_hits <- c(commented_hits, sprintf("%s:%d", f, idx[sort(unique(flagged))]))
    }
  }
}
check("no commented-out code in any @examples block", length(commented_hits) == 0, commented_hits)

# The run-time policies from the same review, grepped in the sources with
# comment and roxygen lines excluded: nothing reads or writes .Random.seed or
# the global environment (seeds go through .dmar_local_seed()), no <<-, no
# options(warn = -1), and no function that writes a file carries a default
# path (the sensitivity family writes only when `filename` is supplied and
# has no `save` switch).
code_lines <- function(f) {
  lines <- readLines(f, warn = FALSE)
  lines[!grepl("^\\s*#", lines)]
}
grep_code <- function(pattern) {
  hits <- character()
  for (f in r_files) {
    m <- grep(pattern, code_lines(f))
    if (length(m)) hits <- c(hits, sprintf("%s (%d line%s)", f, length(m), if (length(m) == 1) "" else "s"))
  }
  hits
}
h <- grep_code("\\.Random\\.seed|\\.GlobalEnv|globalenv\\(\\)")
check("no function reads or writes .Random.seed or the global environment", length(h) == 0, h)
h <- grep_code("<<-")
check("no <<- anywhere in R/", length(h) == 0, h)
# The option is `warn` itself, followed by "=" (`warning.length`, which
# .no_max_refuse() raises while it prints its lavaan code, is another).
h <- grep_code("options\\(\\s*warn\\s*=")
check("no options(warn = ...) anywhere in R/", length(h) == 0, h)
h <- grep_code("(filename|file|path)\\s*=\\s*[\"']")
check("no default file path in any function signature or call", length(h) == 0, h)
h <- grep_code("^\\s*(save\\s*=\\s*(TRUE|FALSE)|.*_sensitivity\\s*<-\\s*function\\(.*\\bsave\\s*=)")
check("no `save` switch in the sensitivity family", length(h) == 0, h)

## ---- DESCRIPTION, NEWS, and the CRAN test subset ----

# The title is the one CRAN granted by email (K. Lauseker, 2026-09-05):
# the expansion of the name with the acronym in parentheses. The
# description does not say "in R" (the 2026-09-04 review). Both are
# checked as strings so a well-meant edit cannot reopen the exchange.
desc <- read.dcf("DESCRIPTION")
title_ok <- identical(unname(desc[, "Title"]), "Design, Measurement, and Analysis in R (DMAR)")
check("DESCRIPTION Title is the granted form", title_ok, unname(desc[, "Title"]))

# Versions are written x.y.z (1.0.0 styling, decided 2026-09-24). A fourth
# component is never used: R orders x.y.z.0 as the same version as x.y.z,
# so CRAN refuses it as an insufficient version, and a development suffix
# is not part of this package's practice.
ver_parts <- unlist(package_version(unname(desc[, "Version"])))
check(sprintf("DESCRIPTION Version %s has exactly three components", desc[, "Version"]),
      length(ver_parts) == 3L)

# CITATION.cff (read by GitHub and by Zenodo for the release DOI) carries the
# released version and its CRAN publication date, both set by hand at each
# release; the version must track DESCRIPTION.
if (file.exists("CITATION.cff")) {
  cff <- yaml::read_yaml("CITATION.cff")
  check(sprintf("CITATION.cff version matches DESCRIPTION (%s)", cff$version),
        identical(as.character(cff$version), unname(desc[, "Version"])))
} else check("CITATION.cff present at the package root", FALSE)
desc_text <- gsub("\\s+", " ", unname(desc[, "Description"]))
h <- grepl("\\bin R\\b", desc_text)
check("DESCRIPTION Description does not say \"in R\"", !h,
      if (h) regmatches(desc_text, regexpr(".{0,40}\\bin R\\b.{0,40}", desc_text)))

# NEWS.md leads with the version in DESCRIPTION, so a release cannot ship
# with last release's notes on top.
news_head <- grep("^# ", readLines("NEWS.md", warn = FALSE), value = TRUE)[1]
news_ver <- sub("^# [A-Za-z]+ ", "", news_head)
check(sprintf("NEWS.md leads with version %s", desc[, "Version"]),
      identical(news_ver, unname(desc[, "Version"])), news_head)

# The curated CRAN test subset (tests/testthat.R, cran_files): every name
# in it must be a file that exists, because a renamed or deleted test file
# drops out of the CRAN path silently (test_check's filter matches
# nothing and reports nothing), and the mechanical detectors the review
# rounds promised must be in it. Files outside the subset are listed for
# information, so a new test file is placed deliberately.
tt_exprs <- parse("tests/testthat.R", keep.source = FALSE)
cf_expr <- Filter(function(e) is.call(e) && identical(e[[1L]], as.name("<-")) &&
                    identical(e[[2L]], as.name("cran_files")), tt_exprs)
if (length(cf_expr) == 1L) {
  cran_files <- eval(cf_expr[[1L]][[3L]], envir = baseenv())
  present <- file.exists(file.path("tests/testthat", paste0("test-", cran_files, ".R")))
  check("every name in the CRAN test subset is an existing test file", all(present),
        cran_files[!present])
  promised <- c("rd_hygiene", "ss_sensitivity_filename")
  check("the promised detectors are in the CRAN test subset", all(promised %in% cran_files),
        setdiff(promised, cran_files))
  all_tests <- sub("^test-(.*)[.]R$", "\\1", basename(test_files))
  outside <- setdiff(all_tests, cran_files)
  cat(sprintf("       CRAN subset: %d of %d test files; %d run locally only\n",
              length(cran_files), length(all_tests), length(outside)))
} else {
  check("tests/testthat.R defines cran_files once", FALSE)
}

## ---- Vignette invariants ----

# Since 2026-09-07 every vignette is a precomputed twin of a .Rmd.orig
# source (the 2026-09-05 pretest's overall checktime NOTE). A .Rmd with no
# .orig is a live vignette that the check would knit on CRAN's clock.
no_orig <- vig_rmd[!file.exists(paste0(vig_rmd, ".orig"))]
check("every vignette is a precomputed twin (has a .Rmd.orig source)",
      length(no_orig) == 0, no_orig)

# Every precomputed twin was knitted from its .orig as it now stands. A
# twin knitted from an earlier .orig means someone edited the source and
# forgot tools/precompute_vignettes.R, so the shipped document silently
# shows stale output. The script ends each twin with a comment recording
# the MD5 checksum of the .orig it knitted, and the check compares that
# record with the .orig on disk. It replaced (2026-11-01) a comparison of
# modification times, which a fresh checkout or a synced folder reorders
# in either direction, so that it could fail a current twin or pass a
# stale one; the checksum misreports neither. A twin knitted any other way
# carries no record and fails, the safe direction.
twin_record <- paste0("<!-- precomputed by tools/precompute_vignettes.R ",
                      "from ([^ ]+), md5 ([0-9a-f]{32})")
stale <- character()
for (o in vig_orig) {
  twin <- sub("[.]orig$", "", o)
  if (!file.exists(twin)) { stale <- c(stale, paste(twin, "missing")); next }
  rec <- grep(twin_record, readLines(twin, warn = FALSE), value = TRUE)
  if (!length(rec)) {
    stale <- c(stale, paste(twin, "carries no record of the source it was",
                            "knitted from; re-knit it with",
                            "tools/precompute_vignettes.R"))
    next
  }
  rec <- rec[length(rec)]
  rec_src <- sub(paste0(".*", twin_record, ".*"), "\\1", rec)
  rec_md5 <- sub(paste0(".*", twin_record, ".*"), "\\2", rec)
  if (!identical(rec_src, basename(o)) ||
      !identical(rec_md5, unname(tools::md5sum(o)))) {
    stale <- c(stale, paste(twin, "was knitted from another version of",
                            basename(o)))
  }
}
check("precomputed vignette twins current (knitted from the .orig as it stands)",
      length(stale) == 0, stale)

# A precomputed twin must contain no executable chunks.
live <- character()
for (o in vig_orig) {
  twin <- sub("[.]orig$", "", o)
  if (file.exists(twin) &&
      any(grepl("^```\\{r", readLines(twin, warn = FALSE)))) {
    live <- c(live, twin)
  }
}
check("precomputed twins carry no live code", length(live) == 0, live)

# Every figure a vignette references exists, and no orphan figure ships.
fig_missing <- character(); fig_orphan <- character()
referenced <- character()
for (v in vig_rmd) {
  lines <- readLines(v, warn = FALSE)
  refs <- regmatches(lines, gregexpr('[A-Za-z0-9_-]+-fig-[A-Za-z0-9_-]+[.]png', lines))
  refs <- unique(unlist(refs))
  referenced <- c(referenced, refs)
  for (r in refs) {
    if (!file.exists(file.path("vignettes", r))) fig_missing <- c(fig_missing, paste(v, "references", r))
  }
}
shipped_figs <- basename(list.files("vignettes", pattern = "-fig-.*[.]png$"))
fig_orphan <- setdiff(shipped_figs, unique(referenced))
check("vignette figures: none missing", length(fig_missing) == 0, fig_missing)
check("vignette figures: none orphaned", length(fig_orphan) == 0,
      paste("vignettes/", fig_orphan, sep = ""))

# No shipped surface points at a vignette that is not in vignettes/.
vig_names <- sub("[.]Rmd$", "", basename(vig_rmd))
dangling <- character()
h <- grep_files('vignette\\("([^"]+)"', c(r_files, vig_rmd))
for (hit in h) {
  parts <- strsplit(hit, ":")[[1]]
  line <- readLines(parts[1], warn = FALSE)[as.integer(parts[2])]
  m <- regmatches(line, gregexpr('vignette\\("([^"]+)"', line))[[1]]
  for (call in m) {
    nm <- sub('vignette\\("([^"]+)".*', "\\1", call)
    if (!nm %in% vig_names) dangling <- c(dangling, paste0(hit, " -> ", nm))
  }
}
check("no vignette() call names a vignette that does not ship", length(dangling) == 0, dangling)

## ---- Convention pre-release greps ----

# Full precision: no round()/signif() feeding a return value. Reported for
# eyes, not auto-failed, because display-layer and example uses are fine;
# fail only on the known-bad pattern round(...fitMeasures...).
h <- grep_files("round\\(.*fitMeasures", r_files)
check("no rounded fitMeasures stored (the cfa_1 defect)", length(h) == 0, h)

# Seed discipline: no function formal defaults a seed to a number.
seed_hits <- character()
for (f in r_files) {
  lines <- readLines(f, warn = FALSE)
  m <- grep("^[^#]*function\\(.*seed *= *[0-9]", lines)
  if (length(m)) seed_hits <- c(seed_hits, sprintf("%s:%d", f, m))
}
check("no baked-in default seed in any signature", length(seed_hits) == 0, seed_hits)

# Example hygiene: no requireNamespace guard inside an @examples block.
guard_hits <- character()
for (f in r_files) {
  lines <- readLines(f, warn = FALSE)
  in_ex <- FALSE
  for (i in seq_along(lines)) {
    l <- lines[i]
    if (grepl("^#' @examples", l)) { in_ex <- TRUE; next }
    if (in_ex && !grepl("^#'", l)) in_ex <- FALSE
    if (in_ex && grepl("requireNamespace", l) && grepl("^#' [^#]", l)) {
      guard_hits <- c(guard_hits, sprintf("%s:%d", f, i))
    }
  }
}
check("no requireNamespace guard executing in @examples", length(guard_hits) == 0, guard_hits)

# Citation currency. The bare words appear in ordinary prose ("happen to
# appear"), so match only the parenthesized citation-status forms.
h <- grep_files("\\(submitted\\)|\\(under review\\)|\\(in preparation\\)|\\(in press\\)|\\(forthcoming\\)|\\(to appear\\)",
                c(r_files, vig_rmd), ignore.case = TRUE)
check("no '(submitted)' style citation status left in shipped text",
      length(h) == 0, h)

# Corrected citations stay corrected. Each row below is a citation form
# that was found wrong and corrected: the work, the wrong form in words, the
# form that replaced it, the commit or commits that made the correction, a
# Perl pattern for the wrong form, and examples of the text the pattern
# must find (wrong) and must pass (right). The check fails when a wrong form
# appears in the R sources (roxygen, other comments, and code, so a message
# string too), in man/ as it stands before roxygen runs below, in a vignette
# source (.Rmd.orig), in the tests (test names, comments, and snapshots),
# in DESCRIPTION, inst/CITATION, or inst/WORDLIST, or in README.Rmd and the
# README.md knitted from it, so that a page copied from an older one, or an
# entry retyped from memory, cannot bring it back. NEWS.md is not read,
# since it describes the corrections in the words of the forms they
# replaced, and neither are the knitted vignette twins, which the twin
# checksum check above ties to the sources read here. Before the scan, each
# pattern is run on its own examples, so a row whose pattern misses its
# wrong form or flags a right one fails here rather than passing silently.
#
# Each file's lines are trimmed and joined before matching, so a form broken
# across lines is still found (a form split between two strings of a
# message is not). A blank line ends a paragraph, and so, in R code, does a
# change from roxygen to another comment or to code. The rows that pair an
# author or a title with a DOI, page range, or equation match only within
# one paragraph, that is, within one reference entry, so the same DOI in
# another work's entry is not flagged; this holds because the package
# separates reference entries by a blank line, and the DOI rows also stop
# at the first DOI after the author's name. A title pattern starts with
# (?i) because a source may set the title in title case. The seed list is
# the corrections made from 2026-10-01 to 2026-11-05; a later correction
# adds a row, with its examples.
#
# Some corrections of the same period changed what a work is cited for
# rather than how it is cited. The work stays correctly cited for other
# claims, so a pattern cannot tell the corrected use from the right ones,
# and these are left out of the table on purpose: the close fit bound
# credited to Shi, Maydeu-Olivares, and DiStefano (2018) rather than to
# Ximénez, Maydeu-Olivares, Shi, and Revuelta (2022) (a270075b, 3ff9b923);
# the IFI credited to Bollen's (1989) article rather than his book
# (3278bf06); the 90 percent RMSEA interval traced to Steiger and Lind
# (1980), with Browne and Cudeck (1993) for its details, and Browne (1984)
# cited for Proposition 4, which holds at any consistent estimates, rather
# than for a need for scale equivariant estimates (bf4dc1dc); Zhang and
# Savalei (2023) credited with the RMSEA and the CFI only (46f3737d); and
# the two-stage comments citing the residual based statistic as Eq. 8 of
# the published Savalei and Bentler (2009), not Eq. 9 of their working
# paper, and saying that the article cites Yuan and Bentler (1998), not Hu,
# Bentler, and Kano (1992), for the statistic's behavior in small samples
# (16387b17); and the plan of a sample size on Fisher's Z scale on
# ?ss_aipe_partial_r credited to Bonett and Wright (2000, Section 3), who
# give it for a simple correlation, rather than to Bonett (2008), which is
# about standardized linear contrasts of means (51b3b98, 2026-11-11).
wrong_citation <- function(work, was, now, commit, pattern, wrong, right) {
  list(work = work, was = was, now = now, commit = commit,
       pattern = pattern, wrong = wrong, right = right)
}
cite_ifi <- "(?i:\\bIFI\\b|\\bincremental\\s+fit\\s+index\\b)"
cite_bollen_eq5 <- paste0(
  "(?:\\bBollen(?:'s|\u2019s)?(?:,\\s+|\\s+\\()1989(?:\\)|,\\s+p\\.\\s*\\d+)?,?",
  "\\s+(?:Eq\\.|Equation)\\s*5(?!\\.?\\d)|\\b(?:Eq\\.|Equation)\\s*5(?!\\.?\\d)",
  "\\s+(?:of|in)\\s+Bollen(?:'s|\u2019s)?(?:,\\s+|\\s+\\()1989\\b)")
# Up to 200 characters within one clause (no semicolon, and no period
# followed by a capital) that do not pass another mention of Bollen, so the
# IFI and an Eq. 5 of Bollen (1989) are tied only when no other citation of
# Bollen stands between them.
cite_gap <- "(?:(?!\\bBollen\\b|;|\\.\\s+[A-Z])[^\\n]){0,200}?"
wrong_citations <- list(
  wrong_citation("Kelley, Maxwell, and Rausch (2003)",
                 "the journal as \"Evaluation and the Health Professions\"",
                 "\"Evaluation & the Health Professions\"", "8ed74773",
                 "(?i)\\bEvaluation\\s+and\\s+the\\s+Health\\s+Professions\\b",
                 wrong = c("\\emph{Evaluation and the Health Professions, 26}(3)",
                           "*Evaluation and the Health Professions*"),
                 right = c("\\emph{Evaluation & the Health Professions, 26}(3)",
                           "*Evaluation & the Health Professions*")),
  wrong_citation("Kelley, Maxwell, and Rausch (2003)",
                 "the subtitle as \"Delineating methods of sample size planning\"",
                 "\"Delineating methods of sample-size planning\"", "8ed74773",
                 "(?i)\\bDelineating\\s+methods\\s+of\\s+sample\\s+size\\s+planning\\b",
                 wrong = c("Delineating methods of sample size planning.",
                           "Delineating Methods of Sample Size Planning"),
                 right = c("Delineating methods of sample-size planning.",
                           "Delineating Methods of Sample-Size Planning")),
  wrong_citation("Steiger and Fouladi (1997)",
                 "the title as ending \"the evaluation of statistical methods\"",
                 "\"... the evaluation of statistical models\"", "232d57c5",
                 "(?i)\\binterval\\s+estimation\\s+and\\s+the\\s+evaluation\\s+of\\s+statistical\\s+methods\\b",
                 wrong = "Noncentrality interval estimation and the evaluation of statistical methods.",
                 right = "Noncentrality interval estimation and the evaluation of statistical models."),
  wrong_citation("Bargmann (1970)", "the name spelled with one n",
                 "\"Bargmann\"",
                 "36f62b84, after d8c132c1 removed most one-n mentions",
                 "\\bBargman\\b",
                 wrong = c("Bargman (1970)", "Bargman's (1970) example",
                           "(Bargman 1970 / Steiger 2004 example)", "Bargman's"),
                 right = c("Bargmann (1970)", "Bargmann's (1970) example",
                           "Bargmann's")),
  wrong_citation("Kiers (1994)",
                 "the DOI 10.1007/BF02294393, which is Atlas and Overall (1994)",
                 "10.1007/BF02294392", "8384662d",
                 "\\bKiers\\b(?:(?!\\b10\\.\\d{4,9}/)[^\\n]){0,300}?\\b10\\.1007/BF02294393\\b",
                 wrong = "Kiers, H. A. L. (1994). Simplimax: Oblique rotation to an optimal target with simple structure. \\emph{Psychometrika, 59}(4), 567--579. \\doi{10.1007/BF02294393}",
                 right = c("Kiers, H. A. L. (1994). Simplimax: Oblique rotation to an optimal target with simple structure. \\emph{Psychometrika, 59}(4), 567--579. \\doi{10.1007/BF02294392}",
                           "Kiers (1994) \\doi{10.1007/BF02294392} Atlas and Overall (1994) \\doi{10.1007/BF02294393}")),
  wrong_citation("Saunders (1961)",
                 "the DOI 10.1007/BF02289799, which is Feldt (1961)",
                 "10.1007/BF02289800", "8384662d",
                 "\\bSaunders\\b(?:(?!\\b10\\.\\d{4,9}/)[^\\n]){0,300}?\\b10\\.1007/BF02289799\\b",
                 wrong = "Saunders, D. R. (1961). The rationale for an \"oblimax\" method of transformation in factor analysis. \\emph{Psychometrika, 26}(3), 317--324. \\doi{10.1007/BF02289799}",
                 right = c("Saunders, D. R. (1961). The rationale for an \"oblimax\" method of transformation in factor analysis. \\emph{Psychometrika, 26}(3), 317--324. \\doi{10.1007/BF02289800}",
                           "Saunders (1961) \\doi{10.1007/BF02289800} Feldt (1961) \\doi{10.1007/BF02289799}")),
  wrong_citation("Fieller (1941), the insulin study",
                 "the year 1940 in the citation",
                 "1941, the year Crossref gives (prose may say that Fieller, 1954, cites it as 1940)",
                 "5d2acb3f",
                 "\\bFieller(?:'s|\u2019s)?(?:,\\s+E\\.\\s*C\\.)?(?:,\\s*|\\s+\\(?)1940[a-z]?\\b",
                 wrong = c("(Fieller, 1940, 1954)", "Fieller, E. C. (1940).",
                           "Fieller (1940)", "(Fieller 1940)", "Fieller's (1940)",
                           "(Fieller, 1940a)"),
                 right = c("(Fieller, 1941, 1954)", "Fieller, E. C. (1941).",
                           "which Fieller (1954) cites with the year 1940",
                           "Fieller's (1954, p. 180) example")),
  wrong_citation("Fieller (1941), the insulin study", "the pages 1--64",
                 "1--54", "5d2acb3f",
                 "(?i)\\bstandardization\\s+of\\s+insulin\\b[^\\n]{0,200}?\\b1\\s*(?:--|-|\u2013)\\s*64\\b",
                 wrong = "The biological standardization of insulin. \\emph{Supplement to the Journal of the Royal Statistical Society, 7}(1), 1--64.",
                 right = "The biological standardization of insulin. \\emph{Supplement to the Journal of the Royal Statistical Society, 7}(1), 1--54."),
  wrong_citation("Bollen (1989)",
                 "Eq. 5 for the incremental fit index",
                 "Eq. 3", "bf4dc1dc",
                 paste0(cite_ifi, cite_gap, cite_bollen_eq5, "|",
                        cite_bollen_eq5, cite_gap, cite_ifi),
                 wrong = c("\\eqn{IFI = (T_b - T)/(T_b - df)} (Bollen, 1989, Eq. 5).",
                           "The IFI is Eq. 5 of Bollen (1989).",
                           "Bollen's (1989) Eq. 5 gives the IFI.",
                           "The IFI is Bollen's (1989) Eq. 5.",
                           "Incremental Fit Index (Bollen, 1989, Eq. 5)",
                           "the IFI (Bollen, 1989, p. 306, Eq. 5)"),
                 right = c("\\eqn{IFI = (T_b - T)/(T_b - df)} (Bollen, 1989, Eq. 3).",
                           "The IFI (Bollen, 1989, Eq. 3); the other index (Bollen, 1989, Eq. 5) differs.",
                           "The IFI (Bollen, 1989, Eq. 5.2)",
                           "the incremental fit index is due to Bollen (1989).",
                           "Bollen, K. A. (1989). A new incremental fit index for general structural equation models.")),
  wrong_citation("MBESS (Kelley, 2007a, 2007b), as DMAR describes its relation to it",
                 "\"reimplementation and expansion\", the wording the decision of 2026-07-30 retired",
                 "\"a more general and greatly expanded reimagining\"",
                 "1a035fe8",
                 "(?i)\\breimplementation\\s+and\\s+expansion\\b",
                 wrong = c("DMAR is a reimplementation and expansion of MBESS",
                           "## Why a Reimplementation and Expansion"),
                 right = c("DMAR is a more general and greatly expanded reimagining of MBESS",
                           "This is the reimplementation of \\code{var.ete()} from \\pkg{MBESS}")),
  # The 2012 Wiley book is the third edition, as its catalog record (ISBN
  # 978-0-470-17896-6) and the Crossref record of its review in
  # International Statistical Review (doi:10.1111/insr.12020_20) have it.
  # The edition must follow the 2012 within one citation: a year in between
  # starts another work's citation.
  wrong_citation("Rencher and Christensen (2012)",
                 "the edition as \"4th ed.\"",
                 "\"3rd ed.\"", "2c253feb",
                 paste0("\\bRencher\\b(?:,\\s+A\\.\\s*C\\.)?,?\\s+(?:&|\\\\&|and)\\s+",
                        "Christensen\\b(?:'s|\u2019s)?(?:,\\s+W\\.\\s*F\\.)?(?:,\\s*|\\s+\\(?)",
                        "2012[a-z]?\\b(?:(?!\\b(?:1[5-9]|20)\\d{2}[a-z]?\\b)[^\\n]){0,200}?",
                        "(?i:\\b(?:4th|fourth)\\s+ed(?:ition)?\\b)"),
                 wrong = c("Rencher, A. C., & Christensen, W. F. (2012). \\emph{Methods of multivariate analysis} (4th ed.). Wiley.",
                           "Rencher and Christensen (2012, 4th ed.)",
                           "(Rencher & Christensen, 2012, fourth edition)",
                           "Rencher and Christensen's (2012) Methods of Multivariate Analysis, 4th Edition"),
                 right = c("Rencher, A. C., & Christensen, W. F. (2012). \\emph{Methods of multivariate analysis} (3rd ed.). Wiley. \\doi{10.1002/9781118391686}",
                           "Rencher and Christensen (2012, 3rd ed.) and a later book (2020, 4th ed.)",
                           "Rencher, A. C. (2002). \\emph{Methods of multivariate analysis} (2nd ed.). Wiley.")),
  # The 2003 Wiley book is the third edition, as the publisher's page (ISBN
  # 978-0-471-36091-9) and catalog records of the print edition have it;
  # Crossref holds no record of it. The edition must follow the 2003 within
  # one citation: a year in between starts another work's citation.
  wrong_citation("Anderson (2003)",
                 "the edition as \"4th ed.\"",
                 "\"3rd ed.\"", "cb80989f",
                 paste0("\\bAnderson\\b(?:'s|\u2019s)?(?:,\\s+T\\.\\s*W\\.)?(?:,\\s*|\\s+\\(?)",
                        "2003[a-z]?\\b(?:(?!\\b(?:1[5-9]|20)\\d{2}[a-z]?\\b)[^\\n]){0,200}?",
                        "(?i:\\b(?:4th|fourth)\\s+ed(?:ition)?\\b)"),
                 wrong = c("Anderson, T. W. (2003). \\emph{An introduction to multivariate statistical analysis} (4th ed.). Wiley. (See Chapter 13.)",
                           "Anderson, T. W. (2003). \\emph{An Introduction to Multivariate Statistical Analysis} (4th ed.). Wiley.",
                           "Anderson (2003, 4th ed.)",
                           "(Anderson, 2003, fourth edition)",
                           "Anderson's (2003) An Introduction to Multivariate Statistical Analysis, 4th Edition"),
                 right = c("Anderson, T. W. (2003). \\emph{An introduction to multivariate statistical analysis} (3rd ed.). Wiley. (See Chapter 13.)",
                           "Anderson, T. W. (2003). \\emph{An introduction to multivariate statistical analysis} (3rd ed.), Sections 4.2 and 4.3. Wiley.",
                           "Anderson (2003, 3rd ed.) and a later book (2020, 4th ed.)",
                           "(Anderson, 2003, Theorem 13.3.1; Muirhead, 1982, Section 9.7), and Kline (2016, 4th ed.)",
                           "Anderson, T. W. (1984). \\emph{An introduction to multivariate statistical analysis} (2nd ed.). Wiley.")),
  # The three rows that follow are the corrections of the 2026-11-09
  # assessment, each checked against the Crossref record of its DOI.
  wrong_citation("Owen (1968)",
                 "the DOI 10.1080/00401706.1968.10490590, which is a letter to the editor (Rothman, 1968)",
                 "10.1080/00401706.1968.10490594", "01ab7a91",
                 "\\bOwen\\b(?:(?!\\b10\\.\\d{4,9}/)[^\\n]){0,300}?\\b10\\.1080/00401706\\.1968\\.10490590\\b",
                 wrong = "Owen, D. B. (1968). A survey of properties and applications of the noncentral t-distribution. \\emph{Technometrics, 10}(3), 445--478. \\doi{10.1080/00401706.1968.10490590}",
                 right = c("Owen, D. B. (1968). A survey of properties and applications of the noncentral t-distribution. \\emph{Technometrics, 10}(3), 445--478. \\doi{10.1080/00401706.1968.10490594}",
                           "Owen (1968) \\doi{10.1080/00401706.1968.10490594} Rothman (1968) \\doi{10.1080/00401706.1968.10490590}")),
  # PubMed (PMID 2348380) gives the same first page as Crossref. The pages
  # must follow the (1990) within one entry, before any DOI.
  wrong_citation("Phillips (1990)", "the first page as 139",
                 "137--144", "01ab7a91",
                 paste0("\\bPhillips\\b(?:(?!\\b(?:1[5-9]|20)\\d{2}[a-z]?\\)\\.)[^\\n]){0,40}?",
                        "\\(1990\\)(?:(?!\\b10\\.\\d{4,9}/)[^\\n]){0,250}?",
                        "\\b139\\s*(?:--|-|\u2013)\\s*144\\b"),
                 wrong = "Phillips, K. F. (1990). Power of the two one-sided tests procedure in bioequivalence. \\emph{Journal of Pharmacokinetics and Biopharmaceutics, 18}(2), 139--144. \\doi{10.1007/BF01063556}",
                 right = "Phillips, K. F. (1990). Power of the two one-sided tests procedure in bioequivalence. \\emph{Journal of Pharmacokinetics and Biopharmaceutics, 18}(2), 137--144. \\doi{10.1007/BF01063556}"),
  wrong_citation("Hofmann (1977)", "the pages 58--66",
                 "96(1), 103--110", "138a2836",
                 paste0("\\bHofmann\\b[^\\n]{0,40}?\\(1977\\)(?:(?!\\b10\\.\\d{4,9}/)[^\\n]){0,200}?",
                        "\\b58\\s*(?:--|-|\u2013)\\s*66\\b"),
                 wrong = "Hofmann, R. J. (1977). Indices descriptive of factor complexity. \\emph{The Journal of General Psychology, 96}, 58--66.",
                 right = "Hofmann, R. J. (1977). Indices descriptive of factor complexity. \\emph{The Journal of General Psychology, 96}(1), 103--110. \\doi{10.1080/00221309.1977.9920803}"),
  # DMAR 1.0.0 gave the expected value of the sample coefficient of
  # variation as an expansion from Sokal and Rohlf (1995); it is Expression 3
  # of Sokal and Braumann (1980, p. 51). No page cites the 1995 book now, so
  # the row matches the citation in either form.
  wrong_citation("Sokal and Rohlf (1995)",
                 "cited as the source of the expected value of the sample coefficient of variation",
                 "Sokal and Braumann (1980, p. 51, Expression 3)", "870c11b0",
                 "\\bSokal\\b[^\\n]{0,15}?(?:&|\\band\\b)\\s+Rohlf\\b[^\\n]{0,15}?\\b1995\\b",
                 wrong = c("(Sokal & Rohlf, 1995)",
                           "Sokal, R. R., & Rohlf, F. J. (1995). \\emph{Biometry: The principles and practice of statistics in biological research} (3rd ed.). W. H. Freeman."),
                 right = c("Expression 3 of Sokal and Braumann (1980, p. 51)",
                           "Sokal, R. R., & Braumann, C. A. (1980). Significance tests for coefficients of variation and variability profiles. \\emph{Systematic Zoology, 29}(1), 50--66. \\doi{10.1093/sysbio/29.1.50}")),
  # In the three rows that follow, the pages must follow the volume directly,
  # with no issue between them, in a help source ("23}, 226") or in Markdown
  # ("23*, 226").
  wrong_citation("Kelley, Darku, and Chattopadhyay (2018)",
                 "the volume without its issue number",
                 "23(2), 226--243", "c0651375",
                 "\\bPsychological\\s+Methods,\\s*23[}*]?,\\s*226\\s*(?:--|-|\u2013)\\s*243\\b",
                 wrong = "Kelley, K., Darku, F. B., & Chattopadhyay, B. (2018). Accuracy in parameter estimation for a general class of effect sizes: A sequential approach. \\emph{Psychological Methods, 23}, 226--243. \\doi{10.1037/met0000127}",
                 right = "Kelley, K., Darku, F. B., & Chattopadhyay, B. (2018). Accuracy in parameter estimation for a general class of effect sizes: A sequential approach. \\emph{Psychological Methods, 23}(2), 226--243. \\doi{10.1037/met0000127}"),
  wrong_citation("Lachowicz, Preacher, and Kelley (2018)",
                 "the volume without its issue number",
                 "23(2), 244--261", "19bd5df4",
                 "\\bPsychological\\s+Methods,\\s*23[}*]?,\\s*244\\s*(?:--|-|\u2013)\\s*261\\b",
                 wrong = "Lachowicz, M. J., Preacher, K. J., & Kelley, K. (2018). A novel measure of effect size for mediation analysis. \\emph{Psychological Methods, 23}, 244--261. \\doi{10.1037/met0000165}",
                 right = "Lachowicz, M. J., Preacher, K. J., & Kelley, K. (2018). A novel measure of effect size for mediation analysis. \\emph{Psychological Methods, 23}(2), 244--261. \\doi{10.1037/met0000165}"),
  wrong_citation("Tofighi and Kelley (2020)",
                 "the volume without its issue number",
                 "25(4), 496--515", "8aa2d881",
                 "\\bPsychological\\s+Methods,\\s*25[}*]?,\\s*496\\s*(?:--|-|\u2013)\\s*515\\b",
                 wrong = "Tofighi, D., & Kelley, K. (2020). Improved inference in mediation analysis: Introducing the model-based constrained optimization procedure. \\emph{Psychological Methods, 25}, 496--515. \\doi{10.1037/met0000259}",
                 right = "Tofighi, D., & Kelley, K. (2020). Improved inference in mediation analysis: Introducing the model-based constrained optimization procedure. \\emph{Psychological Methods, 25}(4), 496--515. \\doi{10.1037/met0000259}")
)
selftest_hits <- character()
for (w in wrong_citations) {
  row <- sprintf("%s, %s (%s)", w$work, w$was, w$commit)
  if (!length(w$wrong) || !length(w$right)) {
    selftest_hits <- c(selftest_hits, sprintf("%s: the row needs wrong and right examples", row))
    next
  }
  missed  <- w$wrong[!grepl(w$pattern, w$wrong, perl = TRUE)]
  flagged <- w$right[grepl(w$pattern, w$right, perl = TRUE)]
  selftest_hits <- c(selftest_hits,
                     sprintf("%s: misses \"%s\"", row, missed),
                     sprintf("%s: flags \"%s\"", row, flagged))
}
check(sprintf("each known wrong citation pattern finds its wrong examples and passes its right ones (%d rows)",
              length(wrong_citations)),
      length(selftest_hits) == 0, selftest_hits)
# The joined text of a file and the position at which each line starts in
# it. Every line break stays one character (a space inside a paragraph, a
# newline between paragraphs), so a match's position still gives its line.
citation_text <- function(f) {
  lines <- readLines(f, warn = FALSE, encoding = "UTF-8")
  n <- length(lines)
  if (!n) return(list(text = "", starts = integer()))
  kind <- rep("text", n)
  if (grepl("[.]R$", f)) {
    kind <- ifelse(grepl("^\\s*#'", lines), "roxygen",
                   ifelse(grepl("^\\s*#", lines), "comment", "code"))
    lines <- sub("^\\s*#+'?", "", lines)
  }
  lines <- trimws(lines)
  same <- nzchar(lines[-n]) & nzchar(lines[-1L]) & kind[-n] == kind[-1L]
  txt <- paste0(lines, c(ifelse(same, " ", "\n"), ""), collapse = "")
  list(text = txt, starts = cumsum(c(1L, nchar(lines) + 1L))[seq_len(n)])
}
citation_files <- c(r_files, rd_files, vig_orig,
                    list.files("tests", pattern = "[.](R|md)$", recursive = TRUE,
                               full.names = TRUE),
                    Filter(file.exists, c("DESCRIPTION", "inst/CITATION",
                                          "inst/WORDLIST", "README.Rmd", "README.md")))
wrong_hits <- character()
for (f in citation_files) {
  ft <- citation_text(f)
  for (w in wrong_citations) {
    m <- gregexpr(w$pattern, ft$text, perl = TRUE)[[1]]
    if (m[1] == -1L) next
    ends <- m + attr(m, "match.length") - 1L
    wrong_hits <- c(wrong_hits, sprintf(
      "%s:%d: %s, %s; corrected to %s in %s",
      f, findInterval(ends, ft$starts), w$work, w$was, w$now, w$commit))
  }
}
check(sprintf("no citation in a form already corrected (%d known wrong forms)",
              length(wrong_citations)),
      length(wrong_hits) == 0, wrong_hits)

# House style: no em dash in package text.
h <- grep_files("—", c(r_files, vig_rmd, "DESCRIPTION", "NEWS.md"))
check("no em dashes", length(h) == 0, h)

## ---- Structural invariants ----

# Every export that can draw random numbers takes seed = NULL (the
# package's development notes, "Reproducibility seed"; the 1.1.0 sweep added it to the sensitivity
# family, procrustes_phi(), ss_aipe_equivalence_r(), and the
# simulate_* generators). An export counts when its body, or the body of a
# non-exported function of the package that it reaches by calls, calls one
# of the generators below (bare or with a namespace prefix), and also when
# it is named *_sensitivity or simulate_*, because four sensitivity
# functions draw only through another export. The source is parsed, not
# loaded, so the check sees what is on disk. A function that reaches a
# generator only through another export (a model function that offers a
# bootstrap through an export that has its own seed, for example) is not
# listed.
rng_calls <- c("rnorm", "runif", "rbinom", "rchisq", "rt", "rgamma", "rbeta",
               "rexp", "rpois", "rmultinom", "rlogis", "rweibull", "rcauchy",
               "rf", "rgeom", "rhyper", "rnbinom", "rsignrank", "rwilcox",
               "rlnorm", "sample", "sample.int", "rWishart", "mvrnorm",
               "rmvnorm")
pkg_fns <- list()
for (f in r_files) {
  for (e in parse(f, keep.source = FALSE)) {
    if (is.call(e) && (identical(e[[1L]], as.name("<-")) ||
                       identical(e[[1L]], as.name("="))) &&
        is.name(e[[2L]]) && is.call(e[[3L]]) &&
        identical(e[[3L]][[1L]], as.name("function")))
      pkg_fns[[as.character(e[[2L]])]] <- e[[3L]]
  }
}
calls_in <- function(x) {
  found <- new.env(parent = emptyenv())
  walk <- function(e) {
    if (is.call(e)) {
      h <- e[[1L]]
      if (is.name(h)) {
        found[[as.character(h)]] <- TRUE
      } else if (is.call(h) && is.name(h[[1L]]) &&
                 as.character(h[[1L]]) %in% c("::", ":::")) {
        found[[as.character(h[[3L]])]] <- TRUE
      }
      for (a in as.list(e)) if (!missing(a)) walk(a)
    } else if (is.pairlist(e) || is.list(e)) {
      for (a in e) if (!missing(a)) walk(a)
    }
  }
  walk(x)
  ls(found, all.names = TRUE)
}
fn_calls <- lapply(pkg_fns, calls_in)
fn_draws <- vapply(fn_calls, function(k) any(k %in% rng_calls), logical(1))
ns_lines <- readLines("NAMESPACE", warn = FALSE)
exported <- sub("^export\\((.*)\\)$", "\\1",
                grep("^export\\(", ns_lines, value = TRUE))
exported <- gsub("\"", "", exported)
can_draw <- function(f) {
  seen <- character()
  todo <- f
  while (length(todo)) {
    g <- todo[1L]
    todo <- todo[-1L]
    if (g %in% seen) next
    seen <- c(seen, g)
    if (isTRUE(fn_draws[[g]])) return(TRUE)
    nxt <- intersect(fn_calls[[g]], names(pkg_fns))
    todo <- c(todo, setdiff(nxt[!(nxt %in% exported)], seen))
  }
  FALSE
}
stochastic <- Filter(function(f) can_draw(f) ||
                       grepl("_sensitivity$|^simulate_", f),
                     intersect(exported, names(pkg_fns)))
no_seed <- Filter(function(f) {
  fm <- pkg_fns[[f]][[2L]]
  !("seed" %in% names(fm) && is.null(fm[["seed"]]))
}, stochastic)
check(sprintf("every export that draws random numbers has seed = NULL (%d exports)",
              length(stochastic)),
      length(no_seed) == 0, no_seed)

# man/ is generated: regenerating roxygen must be a no-op. Checked by hash
# rather than mtime so it is robust to touch.
tmp <- tempfile(); dir.create(tmp)
before <- vapply(rd_files, function(f) unname(tools::md5sum(f)), character(1))
suppressMessages(roxygen2::roxygenise(roclets = "rd"))
after <- vapply(rd_files, function(f) unname(tools::md5sum(f)), character(1))
drifted <- rd_files[before != after]
check("man/ in sync with roxygen sources", length(drifted) == 0, drifted)

# Rd line widths (added 2026-11-11). R CMD check gives a NOTE for a usage
# line over 90 characters or an examples line over 100, because the PDF
# manual truncates them; the release gate of 2026-11-11 failed on eight
# example lines that the renamed arguments had lengthened. The check's own
# function is run here on man/ so the quick gate sees it. Its limits are
# five above the widths it reports (the Rd indentation).
wide <- utils::capture.output(print(tools:::.check_Rd_line_widths(".", limit = c(usage = 95, examples = 105))))
wide <- wide[nzchar(trimws(wide)) & wide != "NULL"]
check("no usage line over 90 or examples line over 100 characters in any help page",
      length(wide) == 0, wide)

## ---- Citation invariants ----

# Every author-year mention in a help page's text must have an entry in that
# page's \references (or \source) section. The scan reads the Rd files the
# sync check above has just regenerated and parses them, so a mention split
# across source lines is still one mention. Three forms are read: "Name
# (1970)" and "Name and Name (1970)"; "(Name, 1970)", "(e.g., Name and Name,
# 1970)", "(Name et al., 1970)", wherever they stand inside a parenthesis;
# and "Name 1970" inside a parenthesis. A
# mention passes when its first surname and its year appear in the same
# reference paragraph (possessives and year suffixes such as 2007a are
# normalized; a hyphenated pair such as Hakstian-Whalen is also tried by its
# first name). Life spans ("1784--1846") are not citations. An exception
# goes in the table below with its reason.
citation_allow <- data.frame(
  page = character(), name = character(), year = character(),
  why = character(), stringsAsFactors = FALSE)
`%||%` <- function(a, b) if (is.null(a)) b else a
rd_flatten <- function(x) {
  if (is.character(x)) return(paste(x, collapse = ""))
  paste(vapply(x, rd_flatten, character(1)), collapse = "")
}
citation_scan <- function(rd_file) {
  rd <- tools::parse_Rd(rd_file)
  tags <- vapply(rd, function(s) attr(s, "Rd_tag") %||% "", character(1))
  ref_tags <- c("\\references", "\\source")
  skip <- c(ref_tags, "\\usage", "\\name", "\\alias", "\\keyword", "\\concept",
            "\\docType", "\\encoding", "\\author", "COMMENT", "TEXT")
  body <- vapply(rd[!tags %in% skip], rd_flatten, character(1))
  ex <- tags == "\\examples"
  body[ex[!tags %in% skip]] <- gsub("(^|\n)\\s*#+", "\\1",
                                    body[ex[!tags %in% skip]])
  body <- gsub("\\s+", " ", paste(body, collapse = " "))
  refs <- paste(vapply(rd[tags %in% ref_tags], rd_flatten, character(1)),
                collapse = "\n\n")
  paras <- gsub("\\s+", " ", strsplit(refs, "\n[[:space:]]*\n")[[1]])
  nm  <- paste0("(?:(?<![[:alpha:]])(?:von|van|de|der|du|ten|da|le|la)\\s+)?",
                "[A-Z][[:alpha:]\u00c0-\u017f'\u2019-]+")
  sep <- "(?:,\\s*&\\s*|,\\s+and\\s+|,\\s*|\\s+and\\s+|\\s*&\\s*)"
  chain <- sprintf("((?:%s%s)*%s)(?:\\s+et\\s+al\\.)?", nm, sep, nm)
  yr <- "((?:1[6-9]|20)\\d{2})[a-z]?(?!\\s*(?:--|\u2013|-|to)\\s*\\d{4})(?![\\d.])"
  narrative <- paste0(chain, "(?:'s|\u2019s|')?\\s+\\(", yr)
  in_paren  <- c(paste0(chain, ",\\s+", yr),
                 paste0("(", nm, ")(?:'s|\u2019s)?\\s+", yr))
  groups <- regmatches(body, gregexpr("\\([^()]*\\)", body, perl = TRUE))[[1]]
  hits <- list()
  m <- regmatches(body, gregexpr(narrative, body, perl = TRUE))[[1]]
  for (h in m) hits[[length(hits) + 1]] <- regmatches(h, regexec(narrative, h, perl = TRUE))[[1]]
  for (g in groups) for (p in in_paren) {
    m <- regmatches(g, gregexpr(p, g, perl = TRUE))[[1]]
    for (h in m) hits[[length(hits) + 1]] <- regmatches(h, regexec(p, h, perl = TRUE))[[1]]
  }
  stop_words <- c("In", "See", "Since", "From", "By", "Version", "Table",
                  "Figure", "Section", "Chapter", "Equation", "Exercise",
                  "Appendix", "Volume", "January", "February", "March",
                  "April", "May", "June", "July", "August", "September",
                  "October", "November", "December", "R", "DMAR", "CRAN")
  out <- character()
  for (g in hits) {
    who <- trimws(strsplit(g[2], sep, perl = TRUE)[[1]][1])
    who <- sub("(?:'s|\u2019s|')$", "", who, perl = TRUE)
    year <- g[3]
    if (who %in% stop_words) next
    cands <- unique(c(who, sub("-.*$", "", who)))
    ok <- any(vapply(paras, function(pp)
      any(vapply(cands, function(cn) grepl(cn, pp, fixed = TRUE), logical(1))) &&
        grepl(sprintf("\\b%s[a-z]?\\b", year), pp, perl = TRUE), logical(1)))
    if (!ok) out <- c(out, sprintf("%s (%s)", who, year))
  }
  unique(out)
}
cite_hits <- character()
for (f in list.files("man", pattern = "[.]Rd$", full.names = TRUE)) {
  page <- sub("[.]Rd$", "", basename(f))
  miss <- citation_scan(f)
  allowed <- paste(citation_allow$name[citation_allow$page == page],
                   citation_allow$year[citation_allow$page == page])
  miss <- miss[!sub("[ ]\\((\\d{4})\\)$", " \\1", miss) %in% allowed]
  if (length(miss)) cite_hits <- c(cite_hits, sprintf("%s: %s", page, paste(miss, collapse = "; ")))
}
check("every author-year mention in a help page has an entry on that page",
      length(cite_hits) == 0, cite_hits)

# Tarball hygiene: build one and confirm nothing ships that should not.
cat("       building the tarball to inspect it (about a minute)...\n")
bld <- tempfile("bld"); dir.create(bld)
old_wd <- setwd(bld)
ok_build <- system2("R", c("CMD", "build", "--no-manual", "--no-build-vignettes", root),
                    stdout = FALSE, stderr = FALSE) == 0
setwd(old_wd)
if (ok_build) {
  tb <- list.files(bld, pattern = "[.]tar[.]gz$", full.names = TRUE)[1]
  contents <- untar(tb, list = TRUE)
  # Root-level .md files other than NEWS, README, and LICENSE are working
  # files (the instruction file, decision notes) and must not ship.
  bad <- c(grep("_problems|[.]orig$|[.]Rmd[.]orig|tools/|dev/|[.]Rcheck|cran-comments",
                contents, value = TRUE),
           grep("^[^/]+/(?!NEWS[.]md$|README[.]md$|LICENSE[.]md$)[^/]+[.]md$",
                contents, value = TRUE, perl = TRUE))
  check("tarball contains no working files", length(bad) == 0, bad)
  # Every top-level file and directory of the tarball is one R CMD check
  # accepts in a source package: the list in tools:::.check_packages() (R
  # 4.6), less tools/, which DMAR keeps out of the tarball. Anything else
  # draws the "Non-standard files/directories found at top level" NOTE.
  # The grep above names only the working paths it knows, which let a
  # methodology/ directory through (the 1.1.0 review); this check needs
  # no list of what to keep out.
  standard_top <- c(
    "DESCRIPTION", "INDEX", "LICENCE", "LICENSE", "LICENCE.note",
    "LICENSE.note", "MD5", "NAMESPACE", "NEWS", "PORTING", "COPYING",
    "COPYING.LIB", "GPL-2", "GPL-3", "BUGS", "Bugs", "ChangeLog",
    "Changelog", "CHANGELOG", "CHANGES", "Changes", "INSTALL", "README",
    "THANKS", "TODO", "ToDo", "SHA256", "SHA256.sig", "INSTALL.windows",
    "README.md", "NEWS.md", "configure", "configure.win", "cleanup",
    "cleanup.win", "configure.ucrt", "cleanup.ucrt", "configure.ac",
    "configure.in", "datafiles", "R", "data", "demo", "exec", "inst", "man",
    "po", "src", "tests", "vignettes", "build", ".aspell", "upstream",
    "java", "noweb")
  entries <- grep("^[^/]+/[^/]+", contents, value = TRUE)
  top <- unique(sub("^[^/]+/([^/]+).*$", "\\1", entries))
  is_dir <- top %in% sub("^[^/]+/([^/]+)/.*$", "\\1",
                         grep("^[^/]+/[^/]+/", entries, value = TRUE))
  nonstandard <- top[!top %in% standard_top]
  check("tarball top level holds only what R CMD check accepts",
        length(nonstandard) == 0,
        paste0(nonstandard, ifelse(is_dir[match(nonstandard, top)], "/", "")))
  # This build omits rendered vignettes, so its size is a floor rather
  # than the full delivered artifact's size. The current floor tripwire
  # is 7.1 MiB; the release gate separately inspects the delivered build.
  sz <- file.info(tb)$size / 1048576
  check(sprintf("tarball floor %.2f Mb (CRAN guideline 5 Mb; tripwire 7.1 Mb)", sz), sz < 7.1)
} else {
  check("R CMD build succeeds", FALSE, "build failed; run it by hand for the error")
}

cat(sprintf("\n%d failure(s).\n", failures))
if (failures > 0) quit(status = 1L)
