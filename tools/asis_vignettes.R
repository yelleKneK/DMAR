# Prebuilt HTML vignettes for the tarball that goes to CRAN.
#
# Every vignette is a precomputed twin (vignettes/<name>.Rmd, the static
# knit of its .Rmd.orig), so R CMD check runs no vignette computation. It
# still renders each twin through rmarkdown and pandoc in the step
# "checking re-building of vignette outputs", which took 18 s of the 1.1.0
# candidate's local check (win-builder has run this package's checks at
# about three times the local time) for HTML that does not change from one
# check to the next. This script moves that rendering off CRAN's clock.
# Given a copy of the package (the gate's git archive of the commit, never
# the working tree), it
#
#   1. renders every twin exactly as R CMD build would, through the same
#      call R CMD build makes (tools::buildVignettes() in a fresh R
#      process, which hands each twin to knitr's rmarkdown engine and so
#      to rmarkdown::render() with the html_vignette format its YAML
#      names), giving a self-contained <name>.html;
#   2. writes beside it <name>.html.asis, which carries the twin's
#      %\VignetteIndexEntry, %\VignetteKeyword, and %\VignetteEncoding
#      lines and %\VignetteEngine{R.rsp::asis};
#   3. removes the twin and its <name>-fig-*.png figures from the copy
#      (the figures are inside the HTML, and the .Rmd.orig regenerates
#      the twin and its figures through tools/precompute_vignettes.R); and
#   4. drops the two lines of .Rbuildignore that keep vignettes/*.html and
#      vignettes/*.Rmd.orig out of a tarball, so that R CMD build sees the
#      HTML and ships each vignette's executable source.
#
# R CMD build then moves each <name>.html to inst/doc, and R CMD check
# copies it rather than rendering it. The tarball's vignettes/ holds, for
# each vignette, <name>.html.asis and <name>.Rmd.orig and nothing else.
# CRAN's policy asks a package to provide its source, and the .Rmd.orig is
# the source of the vignette: tools/precompute_vignettes.R knits it to the
# twin, from which this script renders the HTML. No vignette engine claims
# a file whose name ends in .orig, so R CMD build and R CMD check ship it
# as a plain file and neither runs nor renders it. All 21 vignettes still
# ship, with the same titles and the same HTML. The source tree keeps the
# twins as they are, and its .Rbuildignore still leaves out the HTML and
# the .Rmd.orig sources, so devtools::build(), devtools::build_vignettes(),
# and pkgdown work as before; rendered output is never committed.
# DESCRIPTION lists R.rsp in Suggests and VignetteBuilder for the asis
# engine, so R.rsp must be installed to build the package either way.
#
# tools/release_gate.R runs this on its archive of the commit before
# R CMD build. To build a tarball by hand the same way, from the package
# root:
#
#   rm -rf /tmp/dmar_build && mkdir -p /tmp/dmar_build/DMAR
#   git archive HEAD | tar -x -C /tmp/dmar_build/DMAR
#   Rscript tools/asis_vignettes.R /tmp/dmar_build/DMAR
#   cd /tmp/dmar_build && R CMD build DMAR
#
# This file is not shipped: tools/ is excluded by .Rbuildignore.

asis_vignettes <- function(pkg_dir) {
  pkg_dir <- normalizePath(pkg_dir, mustWork = TRUE)
  vig_dir <- file.path(pkg_dir, "vignettes")
  if (!file.exists(file.path(pkg_dir, "DESCRIPTION")) || !dir.exists(vig_dir)) {
    stop(pkg_dir, " is not a package source with a vignettes/ directory.", call. = FALSE)
  }
  # The rendered HTML is a build product and is never committed, so this
  # refuses a git checkout or worktree; it is meant for an archive of one.
  if (file.exists(file.path(pkg_dir, ".git"))) {
    stop(pkg_dir, " is a git checkout. Run this on a copy made by git archive, ",
         "not on the working tree.", call. = FALSE)
  }
  builders <- read.dcf(file.path(pkg_dir, "DESCRIPTION"), fields = "VignetteBuilder")[1L, 1L]
  if (is.na(builders) || !grepl("\\bR[.]rsp\\b", builders)) {
    stop("DESCRIPTION does not list R.rsp in VignetteBuilder, so R CMD build ",
         "cannot read the .html.asis files.", call. = FALSE)
  }
  if (!requireNamespace("R.rsp", quietly = TRUE)) {
    stop("R.rsp is not installed; install.packages(\"R.rsp\").", call. = FALSE)
  }

  names_v <- sort(sub("[.]Rmd$", "", list.files(vig_dir, pattern = "[.]Rmd$")))
  if (!length(names_v)) stop("No vignettes/*.Rmd to render in ", pkg_dir, call. = FALSE)
  twins <- paste0(names_v, ".Rmd")
  # The executable source of each twin stays in the copy and ships.
  origs <- paste0(names_v, ".Rmd.orig")
  if (!all(file.exists(file.path(vig_dir, origs)))) {
    stop("No executable source for: ",
         paste(twins[!file.exists(file.path(vig_dir, origs))], collapse = ", "), call. = FALSE)
  }

  # The vignette metadata of each twin, read from its YAML header before
  # the twin is removed.
  meta <- lapply(twins, function(f) {
    lines <- readLines(file.path(vig_dir, f), warn = FALSE, encoding = "UTF-8")
    fences <- which(trimws(lines) == "---")
    if (length(fences) < 2L || fences[1L] != 1L) {
      stop(f, " has no YAML header.", call. = FALSE)
    }
    header <- trimws(lines[seq_len(fences[2L])])
    pick <- function(field) grep(sprintf("^%%\\\\Vignette%s\\{", field), header, value = TRUE)
    engine <- pick("Engine")
    if (!identical(engine, "%\\VignetteEngine{knitr::rmarkdown}")) {
      stop(f, " is not a knitr::rmarkdown vignette (", paste(engine, collapse = ", "),
           "); this script renders only those.", call. = FALSE)
    }
    if (!any(grepl("rmarkdown::html_vignette", header, fixed = TRUE))) {
      stop(f, " does not name rmarkdown::html_vignette as its output.", call. = FALSE)
    }
    index <- pick("IndexEntry")
    if (length(index) != 1L) stop(f, " needs exactly one %\\VignetteIndexEntry.", call. = FALSE)
    list(index = index, keyword = pick("Keyword"), encoding = pick("Encoding"),
         title = sub("^%\\\\VignetteIndexEntry\\{(.*)\\}$", "\\1", index))
  })

  # 1. Render, as R CMD build renders: the same call, in an R process of
  # its own started the same way, from the package directory.
  rscript <- file.path(R.home("bin"), "Rscript")
  t0 <- Sys.time()
  out <- withr::with_dir(pkg_dir, system2(
    rscript,
    c("--vanilla", "--default-packages=", "-e",
      shQuote("tools::buildVignettes(dir = '.', tangle = FALSE)")),
    stdout = TRUE, stderr = TRUE))
  render_s <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  status <- attr(out, "status")
  if (!is.null(status) && status != 0L) {
    stop("Rendering the vignettes failed (exit status ", status, "):\n",
         paste(out, collapse = "\n"), call. = FALSE)
  }
  html <- file.path(vig_dir, paste0(names_v, ".html"))
  if (!all(file.exists(html))) {
    stop("No HTML was rendered for: ", paste(names_v[!file.exists(html)], collapse = ", "),
         call. = FALSE)
  }
  # The figures are about to be removed, so the HTML must carry them
  # inside it; a reference to a figure file would be a broken image.
  for (i in seq_along(html)) {
    body <- readLines(html[i], warn = FALSE, encoding = "UTF-8")
    if (any(grepl(sprintf("(src|href)=[\"']%s-fig-", names_v[i]), body))) {
      stop(basename(html[i]), " links to a figure file rather than embedding it.",
           call. = FALSE)
    }
  }

  # 2. The .asis file beside each HTML, and 3. the twin and its figures out.
  for (i in seq_along(names_v)) {
    m <- meta[[i]]
    writeLines(c(m$index, "%\\VignetteEngine{R.rsp::asis}", m$keyword, m$encoding),
               paste0(html[i], ".asis"), useBytes = TRUE)
    files <- list.files(vig_dir)
    figs <- files[startsWith(files, paste0(names_v[i], "-fig-")) & endsWith(files, ".png")]
    unlink(file.path(vig_dir, c(twins[i], figs)))
  }

  # 4. Let the HTML and the executable sources into the tarball, and
  # confirm that no pattern left in .Rbuildignore excludes them (R CMD
  # build matches them with perl = TRUE and ignore.case = TRUE against
  # paths relative to the package root).
  rbi <- file.path(pkg_dir, ".Rbuildignore")
  if (file.exists(rbi)) {
    pats <- readLines(rbi, warn = FALSE)
    pats <- pats[!pats %in% c("^vignettes/.*\\.html$", "^vignettes/.*\\.Rmd\\.orig$")]
    writeLines(pats, rbi)
    shipped <- paste0("vignettes/", c(paste0(names_v, ".html"), paste0(names_v, ".html.asis"), origs))
    live <- pats[nzchar(trimws(pats))]
    hit <- shipped[vapply(shipped, function(p)
      any(vapply(live, grepl, logical(1), x = p, perl = TRUE, ignore.case = TRUE)), logical(1))]
    if (length(hit)) {
      stop(".Rbuildignore still excludes ", paste(hit, collapse = ", "), call. = FALSE)
    }
  }

  structure(data.frame(
    name = names_v,
    title = vapply(meta, `[[`, "", "title"),
    html_md5 = unname(tools::md5sum(html)),
    html_bytes = file.size(html),
    stringsAsFactors = FALSE
  ), render_s = render_s)
}

# Rscript tools/asis_vignettes.R <package copy>
if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) != 1L) {
    stop("Usage: Rscript tools/asis_vignettes.R <copy of the package made by git archive>",
         call. = FALSE)
  }
  res <- asis_vignettes(args[1L])
  cat(sprintf("%d vignettes prebuilt as HTML in %.0f s (%.1f MB of HTML) in %s\n",
              nrow(res), attr(res, "render_s"), sum(res$html_bytes) / 1048576,
              file.path(normalizePath(args[1L]), "vignettes")))
}
