# Regenerate the precomputed vignettes. The executable sources are the
# vignettes/*.Rmd.orig files (tracked in git, and shipped in the release
# tarball beside the prebuilt HTML; see tools/asis_vignettes.R); the
# vignettes/*.Rmd are their statically knitted twins, from which
# tools/asis_vignettes.R renders that HTML, so R CMD check re-runs no
# vignette computation. Re-run this script from
# the package root after editing any .Rmd.orig, then commit both files
# and any *-fig-*.png the knit writes.
#
# Usage, from the package root:
#   Rscript tools/precompute_vignettes.R               # every vignette
#   Rscript tools/precompute_vignettes.R mlmr reliability
# Named vignettes (with or without .Rmd or .Rmd.orig) are knitted alone;
# composite_sem_planning takes about an hour, so knit only what changed.
#
# Each twin ends with a comment that records the MD5 checksum of the
# .Rmd.orig it was knitted from. tools/release_checks.R compares that
# checksum with the .Rmd.orig on disk, so a twin is current exactly when
# it was knitted from the source as it now stands; file modification
# times, which a checkout or a synced folder can reorder, play no part.
# A twin knitted any other way carries no record and fails that check.
#
# After the knit the figures are optimized without loss: every
# vignettes/*.png and man/figures/logo.png is passed to oxipng
# (-o 4 --strip safe) or, when that is not installed, to optipng (-o5).
# Both rewrite a PNG only into another PNG with the same pixels (a figure
# with 256 colors or fewer becomes an indexed-color file; the others are
# compressed again); neither changes a pixel, and lossy quantization is
# never used. The tool, its version, and the bytes before and after are
# printed. When neither tool is installed and python3 can import Pillow,
# tools/optimize_png.py tries other encodings of each figure (rows filtered
# in other ways, and for 256 colors or fewer an indexed-color file holding
# exactly those colors) and keeps the smallest only when it is smaller and
# decodes, with Pillow and with the script's own decoder, to the same
# values in the same channels as the original, with the iCCP, eXIf, and
# pHYs chunks copied byte for byte. A figure it cannot rewrite that way
# (16-bit samples, an APNG, a gAMA, cHRM, sRGB, or sBIT chunk, and the
# other cases its header lists) is left as it is. It prints a line for
# each figure and one with the totals (added 2026-11-12: the unoptimized
# figures, stored twice in the tarball, in vignettes/ and inside the HTML
# vignettes, had put it past the 1e7 bytes at which R CMD check --as-cran
# notes its size). When none of the three is available the step prints a
# message and is skipped, and the figures are as knitr wrote them. The
# tarball floor of tools/release_checks.R cannot see that note: it is built
# without the vignettes, so it leaves out the HTML vignettes and the
# figures inside them. The gate on it is tools/release_gate.R, which builds
# the full tarball, runs R CMD check --as-cran on it, and fails when the
# check reports the size of the tarball.
stamp_prefix <- "<!-- precomputed by tools/precompute_vignettes.R from "

optimize_figures <- function() {
  files <- c(list.files("vignettes", pattern = "[.]png$", full.names = TRUE),
             if (file.exists("man/figures/logo.png")) "man/figures/logo.png")
  if (length(files) == 0L) return(invisible(NULL))
  tool <- unname(Sys.which(c("oxipng", "optipng")))
  tool <- tool[nzchar(tool)][1L]
  if (is.na(tool)) {
    # The fallback: Pillow, through tools/optimize_png.py, when python3 can
    # import it (it may be in the user's site library, so python3 runs
    # without -I here; the script reads only the package's own figures).
    python <- unname(Sys.which("python3"))
    pillow <- nzchar(python) && identical(suppressWarnings(system2(
      python, c("-c", shQuote("import PIL")), stdout = FALSE, stderr = FALSE)), 0L)
    if (pillow) {
      out <- system2(python, c(shQuote("tools/optimize_png.py"), shQuote(files)),
                     stdout = TRUE, stderr = TRUE)
      status <- attr(out, "status")
      if (!is.null(status) && status != 0L) {
        stop("tools/optimize_png.py failed (exit status ", status, "):\n",
             paste(out, collapse = "\n"), call. = FALSE)
      }
      message(paste(c("PNG figures passed to tools/optimize_png.py:", out),
                    collapse = "\n"))
      return(invisible(NULL))
    }
    message("PNG figures not optimized: neither oxipng, optipng, nor ",
            "python3 with Pillow is available (", length(files), " files, ",
            sum(file.size(files)), " bytes, left as written).")
    return(invisible(NULL))
  }
  version <- tryCatch(system2(tool, "--version", stdout = TRUE,
                              stderr = TRUE)[1L],
                      error = function(e) "version not reported")
  before <- sum(file.size(files))
  flags <- if (basename(tool) == "oxipng") {
    c("-o", "4", "--strip", "safe", "--quiet")
  } else {
    c("-o5", "-quiet")
  }
  status <- system2(tool, c(flags, shQuote(files)))
  if (!identical(status, 0L)) {
    stop("The PNG optimizer ", tool, " failed (exit status ", status, ").",
         call. = FALSE)
  }
  after <- sum(file.size(files))
  message(sprintf("PNG figures optimized without loss by %s (%s): %d files, %d bytes before, %d after.",
                  basename(tool), version, length(files), before, after))
  invisible(NULL)
}

args <- commandArgs(trailingOnly = TRUE)
origs <- list.files("vignettes", pattern = "[.]Rmd[.]orig$")
if (length(args)) {
  wanted <- paste0(sub("[.]Rmd([.]orig)?$", "", basename(args)), ".Rmd.orig")
  unknown <- setdiff(wanted, origs)
  if (length(unknown)) {
    stop("No such vignette source: ", paste(unknown, collapse = ", "),
         call. = FALSE)
  }
  origs <- wanted
}

# Several vignettes are knitted each in an R process of its own, so that
# each prints what knitting it alone prints. In one process, the options a
# vignette sets (missing_data and sem_fit_measures set the width to 85)
# and the namespaces it loads (whose S3 registrations print a message
# when a later vignette's package overwrites one) reached every vignette
# knitted after it.
if (length(origs) > 1L) {
  script <- sub("^--file=", "",
                grep("^--file=", commandArgs(FALSE), value = TRUE)[1L])
  rscript <- file.path(R.home("bin"), "Rscript")
  for (orig in origs) {
    status <- system2(rscript, c(shQuote(script), shQuote(orig)),
                      env = "DMAR_PRECOMPUTE_CHILD=1")
    if (!identical(status, 0L)) {
      stop("Knitting ", orig, " failed (exit status ", status, ").",
           call. = FALSE)
    }
  }
  # Once, after every vignette has been knitted.
  optimize_figures()
  quit(save = "no")
}

devtools::load_all(".", quiet = TRUE)
owd <- setwd("vignettes")
for (orig in origs) {
  v <- sub("[.]Rmd[.]orig$", "", orig)
  twin <- paste0(v, ".Rmd")
  # The checksum is taken before the knit, of the bytes the knit reads.
  md5 <- unname(tools::md5sum(orig))
  knitr::opts_chunk$set(fig.path = paste0(v, "-fig-"))
  knitr::knit(orig, output = twin, quiet = TRUE,
              envir = new.env(parent = globalenv()))
  cat("\n", stamp_prefix, orig, ", md5 ", md5,
      ". Edit the .Rmd.orig and re-knit; do not edit this file. -->\n",
      file = twin, sep = "", append = TRUE)
  message(v, " precomputed")
}
setwd(owd)
# A vignette knitted on its own; a process started by the loop above leaves
# the figures to that loop.
if (!nzchar(Sys.getenv("DMAR_PRECOMPUTE_CHILD"))) optimize_figures()
