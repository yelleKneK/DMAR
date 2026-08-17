# Live metafor comparisons removed from the test suite on 2026-08-09.
# Each block recomputes the DMAR value and the metafor oracle and asserts
# agreement at the tolerance the pinned test uses.

## from tests/testthat/test-meta_r.R
local({
  r <- c(.28, .35, .22, .40, .31)
  n <- c(120, 85, 200, 60, 150)
  ours <- DMAR::meta_r(r, n, hartung_knapp = FALSE)
  v <- function(t) ours$value[ours$term == t]
  theirs <- metafor::rma(yi = atanh(r), vi = 1 / (n - 3), method = "REML")
  stopifnot(
    isTRUE(all.equal(v("estimate"), tanh(as.numeric(theirs$beta)),
                     tolerance = 1e-6)),
    isTRUE(all.equal(v("lower_limit"), tanh(as.numeric(theirs$ci.lb)),
                     tolerance = 1e-6)),
    isTRUE(all.equal(v("upper_limit"), tanh(as.numeric(theirs$ci.ub)),
                     tolerance = 1e-6)),
    isTRUE(all.equal(v("tau2"), theirs$tau2, tolerance = 1e-5))
  )
})

## from tests/testthat/test-meta_smd.R
local({
  d  <- DMAR::teacher_expectancy$d
  ne <- DMAR::teacher_expectancy$n_experimental
  nc <- DMAR::teacher_expectancy$n_control
  ours <- DMAR::meta_smd(smd = d, n_1 = ne, n_2 = nc, hartung_knapp = FALSE)
  # metafor's SMD route uses the approximate J = 1 - 3/(4m - 1); feed rma()
  # the exact gamma-function J that DMAR uses and allow the tiny gap.
  df <- ne + nc - 2
  J  <- exp(lgamma(df / 2) - log(sqrt(df / 2)) - lgamma((df - 1) / 2))
  g  <- J * d
  vg <- (ne + nc) / (ne * nc) + g^2 / (2 * (ne + nc))
  theirs <- metafor::rma(yi = g, vi = vg, method = "REML")
  stopifnot(
    isTRUE(all.equal(ours$value[ours$term == "estimate"],
                     as.numeric(theirs$beta), tolerance = 1e-3)),
    abs(ours$value[ours$term == "tau2"] - theirs$tau2) < 1e-4
  )
})

## from tests/testthat/test-meta_es.R
local({
  # meta_es() agrees with metafor across estimators (REML, DL, PM, FE).
  set.seed(113)
  k  <- 12
  vi <- runif(k, 0.01, 0.08)
  yi <- 0.30 + rnorm(k, 0, 0.15) + rnorm(k, 0, sqrt(vi))
  for (m in c("reml", "dl", "pm", "fe")) {
    ours <- DMAR::meta_es(yi, vi, method = m, hartung_knapp = FALSE)
    v <- function(t) ours$value[ours$term == t]
    theirs <- metafor::rma(yi = yi, vi = vi, method = toupper(m))
    stopifnot(
      isTRUE(all.equal(v("estimate"), as.numeric(theirs$beta),
                       tolerance = 1e-6)),
      isTRUE(all.equal(v("se"), theirs$se, tolerance = 1e-5)),
      isTRUE(all.equal(v("tau2"), theirs$tau2, tolerance = 1e-4)),
      isTRUE(all.equal(v("Q"), theirs$QE, tolerance = 1e-8)),
      m == "fe" || isTRUE(all.equal(v("I2"), theirs$I2, tolerance = 1e-3))
    )
  }
})

## from tests/testthat/test-meta_es.R
local({
  # meta_es() Hartung-Knapp matches metafor test = "knha".
  set.seed(113)
  k  <- 12
  vi <- runif(k, 0.01, 0.08)
  yi <- 0.30 + rnorm(k, 0, 0.15) + rnorm(k, 0, sqrt(vi))
  ours <- DMAR::meta_es(yi, vi, method = "reml", hartung_knapp = TRUE)
  v <- function(t) ours$value[ours$term == t]
  theirs <- metafor::rma(yi = yi, vi = vi, method = "REML", test = "knha")
  stopifnot(
    isTRUE(all.equal(v("estimate"), as.numeric(theirs$beta),
                     tolerance = 1e-6)),
    isTRUE(all.equal(v("se"), theirs$se, tolerance = 1e-6)),
    isTRUE(all.equal(v("t"), as.numeric(theirs$zval), tolerance = 1e-6)),
    isTRUE(all.equal(v("p_value"), as.numeric(theirs$pval),
                     tolerance = 1e-5)),
    isTRUE(all.equal(v("lower_limit"), as.numeric(theirs$ci.lb),
                     tolerance = 1e-6)),
    isTRUE(all.equal(v("upper_limit"), as.numeric(theirs$ci.ub),
                     tolerance = 1e-6))
  )
})

## from tests/testthat/test-meta_es.R
local({
  # meta_es() tau2 Q-profile interval and prediction interval match metafor.
  set.seed(113)
  k  <- 12
  vi <- runif(k, 0.01, 0.08)
  yi <- 0.30 + rnorm(k, 0, 0.15) + rnorm(k, 0, sqrt(vi))
  ours <- DMAR::meta_es(yi, vi, method = "reml", hartung_knapp = FALSE)
  v <- function(t) ours$value[ours$term == t]
  theirs <- metafor::rma(yi = yi, vi = vi, method = "REML")
  ci <- confint(theirs)
  # DMAR's prediction interval follows Higgins, Thompson, and Spiegelhalter
  # (2009): t on k - 2 df. metafor's default uses z, so reconstruct the HTS
  # interval from metafor's own tau2 and SE and compare to that.
  half <- qt(0.975, df = k - 2) * sqrt(theirs$tau2 + theirs$se^2)
  stopifnot(
    isTRUE(all.equal(v("tau2_lower"), ci$random["tau^2", "ci.lb"],
                     tolerance = 1e-4)),
    isTRUE(all.equal(v("tau2_upper"), ci$random["tau^2", "ci.ub"],
                     tolerance = 1e-4)),
    isTRUE(all.equal(v("prediction_lower"), as.numeric(theirs$beta) - half,
                     tolerance = 1e-4)),
    isTRUE(all.equal(v("prediction_upper"), as.numeric(theirs$beta) + half,
                     tolerance = 1e-4))
  )
})

## from tests/testthat/test-meta_es.R
local({
  # meta_es() handles the homogeneous boundary like metafor (teacher
  # expectancy); REML sits at the tau2 = 0 boundary here.
  d  <- DMAR::teacher_expectancy$d
  ne <- DMAR::teacher_expectancy$n_experimental
  nc <- DMAR::teacher_expectancy$n_control
  yi <- d
  vi <- (ne + nc) / (ne * nc) + d^2 / (2 * (ne + nc))
  ours <- DMAR::meta_es(yi, vi, method = "reml", hartung_knapp = FALSE)
  v <- function(t) ours$value[ours$term == t]
  theirs <- metafor::rma(yi = yi, vi = vi, method = "REML")
  stopifnot(
    isTRUE(all.equal(v("estimate"), as.numeric(theirs$beta),
                     tolerance = 1e-3)),
    isTRUE(all.equal(v("Q"), theirs$QE, tolerance = 1e-8))
  )
})
