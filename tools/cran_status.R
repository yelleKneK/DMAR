# Post-release monitoring for DMAR: what CRAN's check farm says today.
#
# Once a package is on CRAN the farm re-checks it on every flavor, and a
# result other than OK starts a clock that ends in archival. This script
# reads the package's check results page and its landing page, prints one
# line per flavor, and exits nonzero when anything needs attention: a
# status other than OK, an "Additional issues" section (valgrind, noLD,
# and the like), a flavor whose total check time has grown by more than
# a sixth over the value recorded at acceptance, or a flavor over the
# ten minute line that the incoming pretest enforces. Run it after a
# release lands and whenever a CRAN email arrives:
#
#   Rscript tools/cran_status.R
#
# The acceptance values below are the farm's numbers for 1.0.0 on
# 2026-09-22, the day after publication; update them when a new version
# is published and its results have settled (a week is enough).
#
# This file is not shipped: tools/ is excluded by .Rbuildignore.

pkg <- "DMAR"
accepted <- c(
  `r-devel-linux-x86_64-fedora-clang` = 398.5,
  `r-devel-linux-x86_64-fedora-gcc`   = 403.1,
  `r-release-macos-x86_64`            = 681.0,
  `r-oldrel-macos-x86_64`             = 632.0
)
incoming_limit <- 600   # seconds; the overall checktime line of the incoming pretest

fetch <- function(url) {
  txt <- tryCatch(readLines(url, warn = FALSE), error = function(e) NULL)
  if (is.null(txt)) stop("could not read ", url, call. = FALSE)
  paste(txt, collapse = " ")
}

results_url <- sprintf("https://cran.r-project.org/web/checks/check_results_%s.html", pkg)
landing_url <- sprintf("https://cran.r-project.org/package=%s", pkg)
html <- fetch(results_url)
landing <- fetch(landing_url)

published <- regmatches(landing, regexpr("Published:</td>\\s*<td>\\s*[0-9-]+", landing))
published <- sub(".*<td>\\s*", "", published)
version_cran <- regmatches(landing, regexpr("Version:</td>\\s*<td>\\s*[0-9.-]+", landing))
version_cran <- sub(".*<td>\\s*", "", version_cran)
cat(sprintf("%s %s on CRAN, published %s\n\n", pkg, version_cran, published))

rows <- strsplit(html, "<tr>", fixed = TRUE)[[1]]
rows <- rows[grepl("check_flavors.html#", rows, fixed = TRUE)]
cell <- function(row, pattern) {
  m <- regmatches(row, regexpr(pattern, row, perl = TRUE))
  if (!length(m)) return(NA_character_)
  trimws(gsub("<[^>]+>", "", sub("^[^>]*>", "", m)))
}
tab <- do.call(rbind, lapply(rows, function(r) {
  flavor <- sub(".*check_flavors.html#([^\"]+)\".*", "\\1", r)
  nums <- regmatches(r, gregexpr("class=\"r\">\\s*[0-9.]*\\s*<", r))[[1]]
  nums <- as.numeric(trimws(gsub("[^0-9.]", "", nums)))
  if (length(nums) < 3L) nums <- c(nums, rep(NA_real_, 3L - length(nums)))
  status <- sub(".*class=\"check_([a-z]+)\".*", "\\1", r)
  version <- cell(r, "</a>\\s*</td>\\s*<td>[^<]*<")
  version <- sub(".*<td>\\s*", "", regmatches(r, regexpr("</a>\\s*</td>\\s*<td>\\s*[0-9.-]+", r)))
  data.frame(flavor = flavor, version = version, t_install = nums[1], t_check = nums[2],
             t_total = nums[3], status = toupper(status), stringsAsFactors = FALSE)
}))
if (is.null(tab) || !nrow(tab)) stop("no result rows found on ", results_url, call. = FALSE)

problems <- character()
for (i in seq_len(nrow(tab))) {
  f <- tab$flavor[i]; tt <- tab$t_total[i]
  base <- if (f %in% names(accepted)) accepted[[f]] else NA_real_
  growth <- if (!is.na(base) && !is.na(tt)) tt / base - 1 else NA_real_
  flag <- character()
  if (tab$status[i] != "OK") flag <- c(flag, tab$status[i])
  if (!is.na(growth) && growth > 1 / 6) flag <- c(flag, sprintf("%+.0f%% over acceptance", 100 * growth))
  if (!is.na(tt) && tt > incoming_limit && (is.na(base) || tt > base))
    flag <- c(flag, sprintf("over the %d s incoming line", incoming_limit))
  cat(sprintf("%-36s %-8s %7s s  %-4s %s\n", f, tab$version[i],
              if (is.na(tt)) "" else sprintf("%.0f", tt), tab$status[i],
              if (length(flag)) paste("<-", paste(flag, collapse = "; ")) else ""))
  if (length(flag)) problems <- c(problems, paste(f, paste(flag, collapse = "; "), sep = ": "))
}
if (grepl("Additional issues", html, fixed = TRUE)) {
  problems <- c(problems, "the results page lists Additional issues (valgrind, noLD, or similar); read it")
}
if (any(tab$version != version_cran)) {
  cat("\nnote: some flavors still show an earlier version; the farm has not caught up.\n")
}
cat("\n")
if (length(problems)) {
  cat("ATTENTION NEEDED\n"); for (p in problems) cat("  ", p, "\n", sep = "")
  cat(sprintf("details: %s\n", results_url))
  quit(status = 1L)
}
cat(sprintf("all %d flavors OK; slowest total %.0f s (%s); incoming line %d s\n",
            nrow(tab), max(tab$t_total, na.rm = TRUE), tab$flavor[which.max(tab$t_total)], incoming_limit))
