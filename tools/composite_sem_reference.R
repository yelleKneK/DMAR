# Regenerate the reference values in vignettes/composite_sem_planning.Rmd.orig.
#
# That vignette runs every planner at G = 25 so it knits in about twenty
# seconds. Twenty five replications is far too few to plan against: the same
# first call returns a necessary N near 76 at G = 25 and near 93 at G = 200.
# Every planner result in the vignette is therefore reported beside the same
# call at G = 10000, and those reference values are carried as literals in the
# `reference-values` chunk.
#
# The eight calls are run in parallel, one per core. Expect the better part of
# an hour even so; the latent growth searches are the slow ones.
#
# Run this from the package root after changing either population model, either
# analysis model, the parameter sets, or the width and power targets, and paste
# the printed block into that chunk. If the vignette's numbers and this script's
# output disagree, the vignette is stale.
#
# This file is not shipped: tools/ is excluded by .Rbuildignore.

devtools::load_all(".", quiet = TRUE)

G_REF <- 10000

pop_med <- "
  x ~~ 1*x
  m ~ 0.4*x
  m ~~ 0.84*m
  y ~ 0.35*m + 0.15*x
  y ~~ 0.813*y
"
med_model <- "
  m ~ a*x
  y ~ b*m + cp*x
  ab := a*b
"
pop_lgm <- "
  i =~ 1*t1 + 1*t2 + 1*t3 + 1*t4
  s =~ 0*t1 + 1*t2 + 2*t3 + 3*t4
  i ~~ 1*i
  s ~~ 0.2*s
  i ~~ -0.15*s
  t1 ~~ 0.5*t1; t2 ~~ 0.5*t2; t3 ~~ 0.5*t3; t4 ~~ 0.5*t4
  t1 ~ 0*1; t2 ~ 0*1; t3 ~ 0*1; t4 ~ 0*1
  i ~ 5*1
  s ~ 0.3*1
"
lgm_model <- "
  i =~ 1*t1 + 1*t2 + 1*t3 + 1*t4
  s =~ 0*t1 + 1*t2 + 2*t3 + 3*t4
  i ~~ cov_is*s
  t1 ~ 0*1; t2 ~ 0*1; t3 ~ 0*1; t4 ~ 0*1
  i ~ 1
  s ~ mu_s*1
"

# One entry per displayed result in the vignette, named for the object it
# annotates there.
jobs <- list(
  med_at_100 = function() ss_power_composite_sem(
    model = med_model, pop_model = pop_med, parameters = c("a", "b", "ab"),
    N = 100, G = G_REF, seed = 113),
  med_plan = function() ss_power_composite_sem(
    model = med_model, pop_model = pop_med, parameters = c("a", "b", "ab"),
    desired_power = 0.80, G = G_REF, seed = 113),
  med_plan_cp = function() ss_power_composite_sem(
    model = med_model, pop_model = pop_med, parameters = c("a", "b", "cp", "ab"),
    desired_power = 0.80, G = G_REF, seed = 113),
  med_aipe = function() ss_aipe_composite_sem(
    model = med_model, pop_model = pop_med, parameters = c("a", "b", "ab"),
    desired_width = c(a = 0.25, b = 0.25, ab = 0.15), G = G_REF, seed = 113),
  med_aipe_80 = function() ss_aipe_composite_sem(
    model = med_model, pop_model = pop_med, parameters = c("a", "b", "ab"),
    desired_width = c(a = 0.25, b = 0.25, ab = 0.15), assurance = 0.80,
    G = G_REF, seed = 113),
  lgm_at_150 = function() ss_power_composite_sem(
    model = lgm_model, pop_model = pop_lgm, parameters = c("mu_s", "cov_is"),
    N = 150, G = G_REF, seed = 113),
  lgm_plan = function() ss_power_composite_sem(
    model = lgm_model, pop_model = pop_lgm, parameters = c("mu_s", "cov_is"),
    desired_power = 0.80, G = G_REF, seed = 113),
  lgm_aipe = function() ss_aipe_composite_sem(
    model = lgm_model, pop_model = pop_lgm, parameters = c("mu_s", "cov_is"),
    desired_width = c(mu_s = 0.15, cov_is = 0.25), assurance = 0.80,
    G = G_REF, seed = 113)
)

message("Running ", length(jobs), " planner calls at G = ", G_REF,
        " on ", min(length(jobs), parallel::detectCores() - 1), " cores.")

res <- parallel::mclapply(
  jobs, function(f) tryCatch(f(), error = function(e) conditionMessage(e)),
  mc.cores = min(length(jobs), max(1L, parallel::detectCores() - 1L))
)

failed <- names(res)[vapply(res, is.character, logical(1))]
if (length(failed)) {
  stop("these calls failed: ", paste(failed, collapse = ", "), "\n",
       paste(unlist(res[failed]), collapse = "\n"), call. = FALSE)
}

# The vignette reads single rows out of each table, so carry the whole table as
# a named numeric vector rather than guessing now which rows it will want.
as_vec <- function(x) stats::setNames(x$value, x$term)

cat("\nPaste the following into the `reference-values` chunk of\n",
    "vignettes/composite_sem_planning.Rmd.orig:\n\n", sep = "")
cat("ref <- list(\n")
for (i in seq_along(res)) {
  v <- as_vec(res[[i]])
  cat("  ", names(res)[i], " = c(",
      paste(sprintf('%s = %s', names(v), format(round(v, 4), trim = TRUE)),
            collapse = ", "),
      ")", if (i < length(res)) "," else "", "\n", sep = "")
}
cat(")\n")

saveRDS(res, file.path(tempdir(), "composite_sem_reference.rds"))
message("Raw results also saved to ",
        file.path(tempdir(), "composite_sem_reference.rds"))
