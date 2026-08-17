# Regenerate the precomputed vignettes. The executable sources are the
# vignettes/*.Rmd.orig files (tracked in git, excluded from the tarball);
# the shipped vignettes/*.Rmd are their statically knitted twins, so
# R CMD check re-runs no vignette computation. Re-run this script from
# the package root after editing any .Rmd.orig, then commit both files
# and any *-fig-*.png the knit writes.
devtools::load_all(".", quiet = TRUE)
owd <- setwd("vignettes")
on.exit(setwd(owd))
for (orig in list.files(".", pattern = "[.]Rmd[.]orig$")) {
  v <- sub("[.]Rmd[.]orig$", "", orig)
  knitr::opts_chunk$set(fig.path = paste0(v, "-fig-"))
  knitr::knit(orig, output = paste0(v, ".Rmd"), quiet = TRUE,
              envir = new.env(parent = globalenv()))
  message(v, " precomputed")
}
