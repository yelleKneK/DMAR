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
    isTRUE(all.equal(v("ci_lower"), tanh(as.numeric(theirs$ci.lb)),
                     tolerance = 1e-6)),
    isTRUE(all.equal(v("ci_upper"), tanh(as.numeric(theirs$ci.ub)),
                     tolerance = 1e-6)),
    isTRUE(all.equal(v("tau2"), theirs$tau2, tolerance = 1e-5))
  )
})

## from tests/testthat/test-meta_r.R
local({
  # meta_r() with reliabilities: every row against metafor given the
  # corrected Fisher Z values and their delta method variances, computed
  # here by a central difference and not by DMAR (2026-11-11).
  r   <- c(.21, .35, .42, .18, .30, .47)
  n   <- c(48, 120, 75, 300, 60, 150)
  rxx <- c(.72, .85, .78, .90, .65, .80)
  ryy <- c(.80, .70, .88, .75, .82, .90)
  a   <- rxx * ryy
  k   <- length(r)
  h   <- 1e-6
  z   <- atanh(r)
  dz  <- (atanh(tanh(z + h) / sqrt(a)) - atanh(tanh(z - h) / sqrt(a))) /
    (2 * h)
  yi  <- atanh(r / sqrt(a))
  vi  <- dz^2 / (n - 3)
  for (m in c("fe", "dl", "reml", "pm")) for (hk in c(FALSE, TRUE)) {
    if (m == "fe" && hk) next
    ours <- DMAR::meta_r(r, n, reliability_x = rxx, reliability_y = ryy,
                         method = m, hartung_knapp = hk)
    v <- function(t) ours$value[ours$term == t]
    theirs <- metafor::rma(yi = yi, vi = vi, method = toupper(m),
                           test = if (hk) "knha" else "z",
                           control = list(threshold = 1e-12, tol = 1e-14))
    w <- 1 / (vi + theirs$tau2)
    stopifnot(
      isTRUE(all.equal(v("estimate"), tanh(sum(w * yi) / sum(w)),
                       tolerance = 1e-7)),
      isTRUE(all.equal(v("estimate"), tanh(as.numeric(theirs$beta)),
                       tolerance = 1e-7)),
      isTRUE(all.equal(v("se"), theirs$se, tolerance = 1e-6)),
      isTRUE(all.equal(v("test_statistic"), as.numeric(theirs$zval),
                       tolerance = 1e-6)),
      isTRUE(all.equal(v("p_value"), theirs$pval, tolerance = 1e-5)),
      isTRUE(all.equal(c(v("ci_lower"), v("ci_upper")),
                       tanh(c(theirs$ci.lb, theirs$ci.ub)),
                       tolerance = 1e-6)),
      isTRUE(all.equal(v("tau2"), theirs$tau2, tolerance = 1e-6)),
      isTRUE(all.equal(v("Q"), theirs$QE, tolerance = 1e-10)),
      isTRUE(all.equal(v("Q_p"), theirs$QEp, tolerance = 1e-8))
    )
    if (m == "fe") next
    ci <- confint(theirs, control = list(tol = 1e-14))$random
    half <- qt(0.975, df = k - 2) * sqrt(theirs$tau2 + theirs$se^2)
    stopifnot(
      isTRUE(all.equal(v("I2"), theirs$I2, tolerance = 1e-6)),
      isTRUE(all.equal(v("H2"), theirs$H2, tolerance = 1e-6)),
      isTRUE(all.equal(c(v("tau2_ci_lower"), v("tau2_ci_upper")),
                       unname(ci["tau^2", c("ci.lb", "ci.ub")]),
                       tolerance = 1e-9)),
      isTRUE(all.equal(c(v("I2_ci_lower"), v("I2_ci_upper")),
                       unname(ci["I^2(%)", c("ci.lb", "ci.ub")]),
                       tolerance = 1e-9)),
      isTRUE(all.equal(c(v("prediction_lower"), v("prediction_upper")),
                       tanh(as.numeric(theirs$beta) + c(-1, 1) * half),
                       tolerance = 1e-6))
    )
  }
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
    # REML is run to convergence: at metafor's default threshold of 1e-5
    # its tau2 stops short of the maximum in the fourth significant digit
    # on this fixture, which DMAR's estimate does not (2026-11-09). PM is
    # run to convergence as well (tol = 1e-14, 2026-11-11).
    theirs <- metafor::rma(yi = yi, vi = vi, method = toupper(m),
                           control = list(threshold = 1e-12, tol = 1e-14))
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
  theirs <- metafor::rma(yi = yi, vi = vi, method = "REML", test = "knha",
                         control = list(threshold = 1e-12))
  stopifnot(
    isTRUE(all.equal(v("estimate"), as.numeric(theirs$beta),
                     tolerance = 1e-6)),
    isTRUE(all.equal(v("se"), theirs$se, tolerance = 1e-6)),
    isTRUE(all.equal(v("test_statistic"), as.numeric(theirs$zval), tolerance = 1e-6)),
    isTRUE(all.equal(v("p_value"), as.numeric(theirs$pval),
                     tolerance = 1e-5)),
    isTRUE(all.equal(v("ci_lower"), as.numeric(theirs$ci.lb),
                     tolerance = 1e-6)),
    isTRUE(all.equal(v("ci_upper"), as.numeric(theirs$ci.ub),
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
  theirs <- metafor::rma(yi = yi, vi = vi, method = "REML",
                         control = list(threshold = 1e-12))
  ci <- confint(theirs, control = list(tol = 1e-14))
  # DMAR's prediction interval follows Higgins, Thompson, and Spiegelhalter
  # (2009): t on k - 2 df. metafor's default uses z, so reconstruct the HTS
  # interval from metafor's own tau2 and SE and compare to that.
  half <- qt(0.975, df = k - 2) * sqrt(theirs$tau2 + theirs$se^2)
  stopifnot(
    isTRUE(all.equal(v("tau2_ci_lower"), ci$random["tau^2", "ci.lb"],
                     tolerance = 1e-4)),
    isTRUE(all.equal(v("tau2_ci_upper"), ci$random["tau^2", "ci.ub"],
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

## from tests/testthat/test-meta_es.R
local({
  # The Paule-Mandel estimate and the Q-profile interval against metafor
  # with its root searches run to convergence (2026-11-11). metafor's own
  # searches are on tau2 itself, with an upper end of 100 by default, so
  # the comparison is made at the units of the fixtures; the tests check
  # DMAR's results for other units against a bisection oracle.
  set.seed(113)
  vh <- runif(12, 0.01, 0.08)
  yh <- 0.30 + rnorm(12, 0, 0.15) + rnorm(12, 0, sqrt(vh))
  fixtures <- list(list(yi = yh, vi = vh),
                   list(yi = c(.1, .9, -.4, 1.2, .3, -1),
                        vi = c(.001, .01, .05, .2, .5, 1)))
  for (f in fixtures) {
    ours <- DMAR::meta_es(f$yi, f$vi, method = "pm", hartung_knapp = FALSE)
    v <- function(t) ours$value[ours$term == t]
    theirs <- metafor::rma(yi = f$yi, vi = f$vi, method = "PM",
                           control = list(tol = 1e-14))
    ci <- confint(theirs, control = list(tol = 1e-14))$random
    stopifnot(
      isTRUE(all.equal(v("tau2"), theirs$tau2, tolerance = 1e-9)),
      isTRUE(all.equal(v("estimate"), as.numeric(theirs$beta),
                       tolerance = 1e-9)),
      isTRUE(all.equal(v("se"), theirs$se, tolerance = 1e-9)),
      isTRUE(all.equal(v("I2"), theirs$I2, tolerance = 1e-9)),
      isTRUE(all.equal(c(v("tau2_ci_lower"), v("tau2_ci_upper")),
                       unname(ci["tau^2", c("ci.lb", "ci.ub")]),
                       tolerance = 1e-9)),
      isTRUE(all.equal(c(v("I2_ci_lower"), v("I2_ci_upper")),
                       unname(ci["I^2(%)", c("ci.lb", "ci.ub")]),
                       tolerance = 1e-9))
    )
  }
})

## from tests/testthat/test-ci_ratio_of_means.R
local({
  # ci_ratio_of_means(ci_method = "log") is the response ratio of
  # metafor's escalc("ROM") without its small-sample bias correction: the
  # log ratio is yi, and the squared standard error is vi, "LS" for
  # unequal variances and "HO" for pooled ones (metafor 5.0.1).
  x <- c(37.4, 56.9, 53.1, 44.7, 56.1, 52.1, 52.5, 62.1, 58.9, 40.2, 64.6,
         48.2, 58, 54.7, 49.3, 52.7, 72.4, 56.9, 51.9, 65.9, 55.2, 36.1,
         48.2, 43.6)
  y <- c(31.4, 43.2, 47.7, 44.7, 56, 65.8, 66, 46.3, 41, 34.5, 43.2, 32.6,
         44.6, 46.1, 40.2, 54.5, 39, 27.5, 20, 50.5, -0.5, 40.1, 37, 41.9,
         61.3, 21, 35.5, 40.6, 52.8, 50.8)
  for (pooled in c(FALSE, TRUE)) {
    ours <- DMAR::ci_ratio_of_means(x, y, ci_method = "log",
                                    var_equal = pooled)
    v <- function(t) ours$value[ours$term == t]
    theirs <- metafor::escalc("ROM", m1i = mean(x), sd1i = sd(x),
                              n1i = length(x), m2i = mean(y), sd2i = sd(y),
                              n2i = length(y), correct = FALSE,
                              vtype = if (pooled) "HO" else "LS")
    se <- (log(v("ci_upper")) - log(v("ci_lower"))) /
      (2 * stats::qt(0.975, if (pooled) 52 else attr(ours, "df_limits")[1]))
    stopifnot(
      isTRUE(all.equal(log(v("ratio")), as.numeric(theirs$yi),
                       tolerance = 1e-12)),
      isTRUE(all.equal(se^2, as.numeric(theirs$vi), tolerance = 1e-10)),
      # The constants the test pins.
      isTRUE(all.equal(as.numeric(theirs$vi),
                       if (pooled) 0.00474425205250613 else
                         0.00476047839886595, tolerance = 1e-12))
    )
  }
})

## from tests/testthat/test-ci_ratio_of_means_paired.R
local({
  # ci_ratio_of_means_paired(ci_method = "log") is the response ratio for
  # correlated means, metafor's escalc("ROMC") without its bias correction.
  # The clustered data of the test (see tools/oracle_snippets/sandwich.R),
  # used here without the clusters or the baseline.
  set.seed(19540512)
  cluster <- rep(1:12, times = c(4, 5, 3, 6, 4, 5, 4, 3, 6, 5, 4, 4))
  n <- length(cluster)
  shared <- stats::rnorm(12)[cluster]
  z <- round(10 + shared + stats::rnorm(n), 2)
  y <- round(z + 3 + 0.5 * shared + stats::rnorm(n), 2)
  x <- round(z + 2.4 + 0.5 * shared + stats::rnorm(n), 2)
  ours <- DMAR::ci_ratio_of_means_paired(x, y, ci_method = "log")
  v <- function(t) ours$value[ours$term == t]
  theirs <- metafor::escalc("ROMC", m1i = mean(x), sd1i = sd(x),
                            m2i = mean(y), sd2i = sd(y), ni = n,
                            ri = cor(x, y), correct = FALSE)
  se <- (log(v("ci_upper")) - log(v("ci_lower"))) /
    (2 * stats::qt(0.975, n - 1))
  stopifnot(
    isTRUE(all.equal(log(v("ratio")), as.numeric(theirs$yi),
                     tolerance = 1e-12)),
    isTRUE(all.equal(se^2, as.numeric(theirs$vi), tolerance = 1e-10)),
    # The constants the test pins.
    isTRUE(all.equal(c(as.numeric(theirs$yi), as.numeric(theirs$vi)),
                     c(-0.049602638656155065, 0.000231707617997902),
                     tolerance = 1e-12))
  )
})
