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

# No example runs a bootstrap confidence interval. (The randomization tests
# resample by definition and are exempt.)
boot_hits <- character()
for (f in setdiff(r_files, c("R/randomization_test.R", "R/randomization_test_paired.R"))) {
  lines <- readLines(f, warn = FALSE)
  in_ex <- FALSE
  for (i in seq_along(lines)) {
    l <- lines[i]
    if (grepl("^#' @examples", l)) { in_ex <- TRUE; next }
    if (in_ex && !grepl("^#'", l)) in_ex <- FALSE
    if (in_ex && grepl("^#' [^#]", l) &&
        grepl('ci_method *= *"(percentile|bca|boot|bootstrap)"|boot *= *TRUE|method *= *"bootstrap"', l)) {
      boot_hits <- c(boot_hits, sprintf("%s:%d", f, i))
    }
  }
}
check("no example executes a bootstrap interval", length(boot_hits) == 0, boot_hits)

## ---- Vignette invariants ----

# Every precomputed twin is newer than its source. A twin older than its
# .orig means someone edited the source and forgot tools/precompute_vignettes.R,
# so the shipped document silently shows stale output.
stale <- character()
for (o in vig_orig) {
  twin <- sub("[.]orig$", "", o)
  if (!file.exists(twin)) { stale <- c(stale, paste(twin, "missing")); next }
  if (file.mtime(o) > file.mtime(twin)) stale <- c(stale, paste(twin, "older than its .orig"))
}
check("precomputed vignette twins current", length(stale) == 0, stale)

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

# House style: no em dash in package text.
h <- grep_files("—", c(r_files, vig_rmd, "DESCRIPTION", "NEWS.md"))
check("no em dashes", length(h) == 0, h)

## ---- Structural invariants ----

# man/ is generated: regenerating roxygen must be a no-op. Checked by hash
# rather than mtime so it is robust to touch.
tmp <- tempfile(); dir.create(tmp)
before <- vapply(rd_files, function(f) unname(tools::md5sum(f)), character(1))
suppressMessages(roxygen2::roxygenise(roclets = "rd"))
after <- vapply(rd_files, function(f) unname(tools::md5sum(f)), character(1))
drifted <- rd_files[before != after]
check("man/ in sync with roxygen sources", length(drifted) == 0, drifted)

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
  bad <- grep("_problems|[.]orig$|[.]Rmd[.]orig|tools/|dev/|[.]Rcheck|[A-Z]{4,}[.]md|cran-comments",
              contents, value = TRUE)
  check("tarball contains no working files", length(bad) == 0, bad)
  # Note: built with --no-build-vignettes, so this size understates the real
  # tarball, which adds the rendered vignette HTML. The figure is a floor.
  sz <- file.info(tb)$size / 1048576
  check(sprintf("tarball floor %.2f Mb below the 5 Mb guideline", sz), sz < 5)
} else {
  check("R CMD build succeeds", FALSE, "build failed; run it by hand for the error")
}

cat(sprintf("\n%d failure(s).\n", failures))
if (failures > 0) quit(status = 1L)
