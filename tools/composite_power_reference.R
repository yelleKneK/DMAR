# Regenerate the reference values in vignettes/composite_power_ancova.Rmd.
#
# That vignette runs its simulations at reps = 25 so it knits in under a
# second, which is far too few to estimate a power. Every simulated quantity
# is therefore reported beside the same quantity at reps = 10000, and those
# reference values are carried in the vignette as literals in the
# `reference-values` chunk, because a 10,000-replication sweep takes minutes
# and cannot run while the document is being knitted.
#
# Set DMAR_REF_REPS to a small number to smoke test the plumbing without
# waiting for the full sweep, for example DMAR_REF_REPS=50 Rscript <this file>.
#
# Run this script from the package root after changing the population, the
# design, the contrasts, or anything else the simulations depend on, and paste
# the printed block into that chunk. If the vignette's numbers and this
# script's output ever disagree, the vignette is stale.
#
# This file is not shipped: tools/ is excluded by .Rbuildignore.

devtools::load_all(".", quiet = TRUE)
stopifnot(requireNamespace("knitr", quietly = TRUE))

R <- as.integer(Sys.getenv("DMAR_REF_REPS", "10000"))
vignette_path <- "vignettes/composite_power_ancova.Rmd"

# Take the population, the contrasts, and the two simulation functions from the
# vignette itself rather than restating them here, so the reference values
# cannot drift from the document they annotate. Chunks are evaluated one at a
# time because purl() emits the display-only chunks out of dependency order.
# purl() resolves each chunk's eval= option in the global environment, so the
# flags the vignette branches on have to exist there BEFORE purling. Setting
# them afterwards, or only inside the evaluation environment, silently drops
# every `eval = has_ggplot2` chunk from the extracted script, including the one
# that defines the regression coefficients the three way simulation closes over.
has_ggplot2   <- requireNamespace("ggplot2", quietly = TRUE)
has_patchwork <- requireNamespace("patchwork", quietly = TRUE)
if (has_ggplot2) library(ggplot2)

purled <- tempfile(fileext = ".R")
knitr::purl(vignette_path, output = purled, quiet = TRUE)
env <- new.env(parent = globalenv())
for (e in parse(purled)) {
  try(suppressWarnings(suppressMessages(eval(e, env))), silent = TRUE)
}

need <- c("composite_power", "simulate_threeway_power", "mu", "sigma2",
          "targets", "n_design", "N_3way", "gamma_small",
          "b0", "bA", "bB", "bAB")
missing <- need[!vapply(need, exists, logical(1), envir = env)]
if (length(missing)) {
  stop("the vignette did not define: ", paste(missing, collapse = ", "),
       ". A chunk it depends on failed to evaluate.", call. = FALSE)
}
for (nm in need) assign(nm, get(nm, envir = env))

message("Running ", R, " replications per cell. This takes a few minutes.")

mc      <- composite_power(n_design, mu, sigma2, reps = R)
n_seq   <- seq(20, 40, by = 2)
comp    <- vapply(n_seq, function(nn)
             composite_power(nn, mu, sigma2, reps = R)$composite, numeric(1))
tw      <- simulate_threeway_power(N_3way, gamma = gamma_small, reps = R)
N_try   <- c(N_3way, 420, 440, 460)
tw_by_N <- vapply(N_try, function(NN)
             simulate_threeway_power(NN, gamma = gamma_small, reps = R)[["three_way"]],
             numeric(1))

fmt <- function(x) paste(sprintf("%.4f", x), collapse = ", ")
cat("\nPaste the following into the `reference-values` chunk of\n",
    vignette_path, ":\n\n", sep = "")
cat("ref <- list(\n")
cat("  marginal   = c(", fmt(mc$marginal), "),\n", sep = "")
cat("  composite  = ", fmt(mc$composite), ",\n", sep = "")
cat("  product    = ", fmt(mc$product), ",\n", sep = "")
cat("  n_seq      = seq(20, 40, by = 2),\n")
cat("  comp_seq   = c(", fmt(comp), "),\n", sep = "")
cat("  n_comp80   = ", n_seq[which(comp >= 0.80)[1]], ",\n", sep = "")
cat("  threeway   = c(three_way = ", fmt(tw[["three_way"]]),
    ", two_way = ", fmt(tw[["two_way"]]), "),\n", sep = "")
cat("  N_try      = c(", paste(N_try, collapse = ", "), "),\n", sep = "")
cat("  tw_by_N    = c(", fmt(tw_by_N), "),\n", sep = "")
cat("  N_mc80     = ", N_try[which(tw_by_N >= 0.80)[1]], "\n", sep = "")
cat(")\n")
