v <- function(tab, t) tab$value[tab$term == t]

# (yi, vi) fixture: the teacher expectancy studies in the d metric.
.te_yivi <- function() {
  d  <- teacher_expectancy$d
  ne <- teacher_expectancy$n_experimental
  nc <- teacher_expectancy$n_control
  list(yi = d, vi = (ne + nc) / (ne * nc) + d^2 / (2 * (ne + nc)))
}

test_that("meta_es() returns the documented tidy schema", {
  f <- .te_yivi()
  res <- meta_es(f$yi, f$vi)
  expect_s3_class(res, "dmar_tbl")
  expect_true(all(c("estimate", "se", "test_statistic", "p_value", "ci_lower",
                    "ci_upper", "prediction_lower", "prediction_upper",
                    "tau2", "tau2_ci_lower", "tau2_ci_upper", "tau", "I2",
                    "I2_ci_lower", "I2_ci_upper", "H2",
                    "Q", "Q_df", "Q_p", "n_studies") %in% res$term))
  # Only tau2 and I2 carry interval rows; ?meta_es documents the table
  # accordingly. If tau or H2 limits are ever added, update the help page
  # alongside this expectation.
  expect_false(any(c("tau_lower", "tau_upper", "H2_lower", "H2_upper")
                   %in% res$term))
  expect_identical(attr(res, "method"), "reml")
  expect_true(attr(res, "hartung_knapp"))
  expect_identical(attr(res, "conf_level"), 0.95)
  expect_equal(v(res, "n_studies"), 19)
  expect_equal(v(res, "tau"), sqrt(v(res, "tau2")))
})

# Heterogeneous fixture: true tau well inside the parameter space, so the
# REML optimum is interior and estimator agreement can be held tight.
.het_yivi <- function() {
  set.seed(113)
  k  <- 12
  vi <- runif(k, 0.01, 0.08)
  yi <- 0.30 + rnorm(k, 0, 0.15) + rnorm(k, 0, sqrt(vi))
  list(yi = yi, vi = vi)
}

test_that("meta_es() agrees with metafor across estimators (REML, DL, PM, FE)", {
  f <- .het_yivi()
  # Pinned from metafor::rma(yi, vi, method = toupper(m)) on this fixture
  # (metafor 5.0.1, 2026-08-09); live comparison in tools/oracle_checks.R.
  # The REML row is from metafor's Fisher scoring run to convergence
  # (control = list(threshold = 1e-12), 2026-11-09). At metafor's default
  # threshold of 1e-5 its tau2 is 0.0279436, short of the maximum in the
  # fourth significant digit, as DMAR's was before the tolerance of its
  # own search was tightened. The PM row is from metafor's root search run
  # to convergence (control = list(tol = 1e-14), 2026-11-11): at its
  # default tolerance, the fourth root of the machine epsilon, its tau2 is
  # 0.04310566, short of the root in the sixth significant digit, which
  # the PM search of DMAR, now made on the scale of the data, is not.
  pinned <- list(
    reml = list(estimate = 0.4164533038082895, se = 0.07425684681461567,
                tau2 = 0.0279387014282716, Q = 20.61130014356519,
                I2 = 44.77163976595503),
    dl   = list(estimate = 0.4158314259576409, se = 0.07556718741348845,
                tau2 = 0.03011305513068721, Q = 20.61130014356519,
                I2 = 46.63121722850569),
    pm   = list(estimate = 0.4129886282193178, se = 0.08285465366310646,
                tau2 = 0.04310632185035629, Q = 20.61130014356519,
                I2 = 55.57065485525343),
    fe   = list(estimate = 0.436158095665153, se = 0.05263784460136482,
                tau2 = 0, Q = 20.61130014356519)
  )
  for (m in names(pinned)) {
    ours <- meta_es(f$yi, f$vi, method = m, hartung_knapp = FALSE)
    ref  <- pinned[[m]]
    expect_equal(v(ours, "estimate"), ref$estimate,
                 tolerance = 1e-6, label = paste("estimate", m))
    expect_equal(v(ours, "se"), ref$se, tolerance = 1e-5,
                 label = paste("se", m))
    expect_equal(v(ours, "tau2"), ref$tau2, tolerance = 1e-4,
                 label = paste("tau2", m))
    expect_equal(v(ours, "Q"), ref$Q, tolerance = 1e-8,
                 label = paste("Q", m))
    if (m != "fe") {
      expect_equal(v(ours, "I2"), ref$I2, tolerance = 1e-3,
                   label = paste("I2", m))
    }
  }
})

# On CRAN the file stops here, after the anchor above; the blocks below run
# in the full local suite (see cran_files in tests/testthat.R).
skip_on_cran()

test_that("meta_es() Hartung-Knapp matches metafor test='knha'", {
  f <- .het_yivi()
  ours <- meta_es(f$yi, f$vi, method = "reml", hartung_knapp = TRUE)
  # Pinned from metafor::rma(yi, vi, method = "REML", test = "knha",
  # control = list(threshold = 1e-12)), REML run to convergence
  # (metafor 5.0.1, 2026-11-09); live comparison in tools/oracle_checks.R.
  expect_equal(v(ours, "estimate"), 0.4164533038082895, tolerance = 1e-6)
  expect_equal(v(ours, "se"), 0.08099516083465047, tolerance = 1e-6)
  expect_equal(v(ours, "test_statistic"), 5.141705992268701, tolerance = 1e-6)
  expect_equal(v(ours, "p_value"), 0.0003223777640718407, tolerance = 1e-5)
  expect_equal(v(ours, "ci_lower"), 0.2381841567719883,
               tolerance = 1e-6)
  expect_equal(v(ours, "ci_upper"), 0.5947224508445907,
               tolerance = 1e-6)
})

test_that("meta_es() tau2 Q-profile interval and prediction interval match metafor", {
  f <- .het_yivi()
  ours <- meta_es(f$yi, f$vi, method = "reml", hartung_knapp = FALSE)
  # Pinned from confint(metafor::rma(yi, vi, method = "REML")), the
  # Q-profile interval for tau2, with its root search run to convergence
  # (control = list(tol = 1e-14); metafor 5.0.1, 2026-11-11). At metafor's
  # default tolerance the upper limit is 0.2262314670, short of the root
  # in the sixth significant digit; the pin was that value until the
  # Q-profile search of DMAR was made on the scale of the data.
  expect_equal(v(ours, "tau2_ci_lower"), 0, tolerance = 1e-4)
  expect_equal(v(ours, "tau2_ci_upper"), 0.2262306486670355, tolerance = 1e-8)
  # DMAR's prediction interval follows Higgins, Thompson, and Spiegelhalter
  # (2009): t on k - 2 df. metafor's default uses z, so the pinned values
  # reconstruct the HTS interval beta -/+ qt(.975, k - 2) *
  # sqrt(tau2 + se^2) from metafor's own tau2 and SE, with REML run to
  # convergence (control = list(threshold = 1e-12), 2026-11-09).
  expect_equal(v(ours, "prediction_lower"),
               0.008924449635623333, tolerance = 1e-4)
  expect_equal(v(ours, "prediction_upper"),
               0.8239821579809556, tolerance = 1e-4)
})

test_that("meta_es() handles the homogeneous boundary like metafor (teacher expectancy)", {
  f <- .te_yivi()
  ours <- meta_es(f$yi, f$vi, method = "reml", hartung_knapp = FALSE)
  # REML sits at the tau2 = 0 boundary here; estimates agree to the width
  # of the optimizers' stopping rules.
  # Pinned from metafor::rma(yi, vi, method = "REML") (metafor 5.0.1,
  # 2026-08-09); live comparison in tools/oracle_checks.R.
  expect_lt(v(ours, "tau2"), 1e-4)
  expect_equal(v(ours, "estimate"), 0.05494697909581286, tolerance = 1e-3)
  expect_equal(v(ours, "Q"), 17.0673484133128, tolerance = 1e-8)
})

test_that("meta_es() reproduces the Raudenbush (1984) heterogeneity by hand", {
  # Study-level (18 studies): the paper reports chi-square(17) = 14.65; with
  # the standard large-sample d variances the statistic is 14.9 (his hand
  # computations rounded). The data are homogeneous-looking overall, which
  # is exactly why the focused contrast mattered.
  s  <- teacher_expectancy[-c(4, 5), ]
  d  <- append(s$d, .52, after = 3)
  ne <- append(s$n_experimental, 22, after = 3)
  nc <- append(s$n_control, 22, after = 3)
  vi <- (ne + nc) / (ne * nc) + d^2 / (2 * (ne + nc))
  res <- meta_es(d, vi, method = "fe")
  expect_equal(v(res, "Q"), 14.65, tolerance = 0.025)
  expect_equal(v(res, "Q_df"), 17)
  expect_gt(v(res, "Q_p"), 0.5)
})

test_that("meta_es() validates its arguments", {
  skip_on_cran()  # The extended validation grid runs in the full local suite.
  expect_error(meta_es(c(1, 2), c(0.1, -0.1)), "positive sampling variance")
  expect_error(meta_es(1, 0.1), "two or more")
  expect_error(meta_es(c(1, 2), c(0.1, 0.1), conf_level = 1.2),
               "\\(0, 1\\)")
  expect_error(meta_es(c(1, 2), c(0.1, 0.1), hartung_knapp = NA),
               "TRUE or FALSE")
})

test_that("meta_es() refuses effect sizes or variances that are not finite", {
  skip_on_cran()
  for (m in c("reml", "pm", "dl", "fe")) {
    expect_error(meta_es(c(.1, Inf, .3), rep(.1, 3), method = m),
                 "finite effect sizes", label = m)
    expect_error(meta_es(c(.1, -Inf, .3), rep(.1, 3), method = m),
                 "effect sizes", label = m)
    expect_error(meta_es(c(.1, .2, .3), c(.1, Inf, .1), method = m),
                 "finite sampling variance", label = m)
  }
  expect_error(plot_forest(c(.1, .2, .3), c(.1, Inf, .1)),
               "sampling variance")
})

# Identical effect sizes under the Hartung-Knapp adjustment: the weighted
# residuals are all zero, so the adjusted standard error is zero. The
# numbers are those metafor's test = "knha" returns (se 0, t Inf, p 0, and
# limits at the estimate; metafor 5.0.1); DMAR adds a classed warning.
test_that("identical effect sizes under Hartung-Knapp keep their numbers and warn once", {
  regular <- meta_es(c(.30, .32, .29, .31, .33), c(.04, .05, .03, .06, .04))
  seen <- new.env(parent = emptyenv())
  seen$conditions <- list()
  res <- withCallingHandlers(
    meta_es(rep(.5, 4), rep(.04, 4)),
    warning = function(w) {
      seen$conditions <- c(seen$conditions, list(w))
      invokeRestart("muffleWarning")
    })
  expect_length(seen$conditions, 1L)
  expect_s3_class(seen$conditions[[1]], "dmar_meta_hk_zero_se")
  expect_match(conditionMessage(seen$conditions[[1]]),
               "standard error is zero")
  expect_match(conditionMessage(seen$conditions[[1]]), "no width")
  expect_identical(nrow(res), 20L)
  expect_identical(res$term, regular$term)
  expect_equal(v(res, "estimate"), 0.5)
  expect_identical(v(res, "se"), 0)
  expect_identical(v(res, "test_statistic"), Inf)
  expect_identical(v(res, "p_value"), 0)
  expect_equal(c(v(res, "ci_lower"), v(res, "ci_upper")), c(0.5, 0.5))
  expect_equal(c(v(res, "prediction_lower"), v(res, "prediction_upper")),
               c(0.5, 0.5))
  expect_identical(v(res, "tau2"), 0)
  expect_equal(v(res, "Q"), 0)
  expect_equal(v(res, "n_studies"), 4)
  # The table prints.
  expect_output(print(res), "estimate")
})

test_that("the zero standard error warning covers every identical case and no other", {
  # A common effect size of zero (the statistic is 0 / 0), two studies with
  # unequal variances, and effect sizes equal up to rounding error
  # (0.1 + 0.2 is not 0.3 in floating point; the standard error is 2e-17).
  expect_warning(z <- meta_es(rep(0, 4), rep(.04, 4)),
                 class = "dmar_meta_hk_zero_se")
  expect_identical(v(z, "se"), 0)
  expect_true(is.nan(v(z, "test_statistic")))
  expect_equal(c(v(z, "ci_lower"), v(z, "ci_upper")), c(0, 0))
  expect_warning(meta_es(c(.5, .5), c(.04, .09)),
                 class = "dmar_meta_hk_zero_se")
  expect_warning(r <- meta_es(c(0.1 + 0.2, 0.3, 0.3), rep(.04, 3)),
                 class = "dmar_meta_hk_zero_se")
  expect_lt(v(r, "se"), 1e-15)
  for (m in c("pm", "dl")) {
    expect_warning(meta_es(rep(.5, 4), rep(.04, 4), method = m),
                   class = "dmar_meta_hk_zero_se")
  }
  # Effect sizes that differ in the ninth decimal are not identical.
  expect_no_warning(near <- meta_es(c(.5, .5, .5 + 1e-9), rep(.04, 3)))
  expect_true(all(is.finite(c(v(near, "ci_lower"), v(near, "ci_upper")))))
  expect_gt(v(near, "ci_upper") - v(near, "ci_lower"), 0)
  # Without the adjustment, and for the fixed effect analysis, the interval
  # comes from the sampling variances and nothing is degenerate.
  expect_no_warning(
    plain <- meta_es(rep(.5, 4), rep(.04, 4), hartung_knapp = FALSE))
  expect_equal(c(v(plain, "ci_lower"), v(plain, "ci_upper")),
               0.5 + c(-1, 1) * qnorm(.975) * 0.1)
  expect_no_warning(meta_es(rep(.5, 4), rep(.04, 4), method = "fe"))
})

test_that("the Hartung-Knapp interval is not bounded below by the normal interval", {
  # The weighted residual mean square is 0.0069 here, so the adjusted
  # interval is far narrower than the normal one. Pinned from
  # metafor::rma(yi, vi, method = "REML", test = "knha") (metafor 5.0.1,
  # 2026-11-09).
  yi <- c(.30, .32, .29, .31, .33)
  vi <- c(.04, .05, .03, .06, .04)
  expect_no_warning(hk <- meta_es(yi, vi))
  expect_equal(v(hk, "se"), 0.007561157543, tolerance = 1e-5)
  expect_equal(v(hk, "ci_lower"), 0.287201305594, tolerance = 1e-5)
  expect_equal(v(hk, "ci_upper"), 0.329187583295, tolerance = 1e-5)
  wald <- meta_es(yi, vi, hartung_knapp = FALSE)
  expect_lt(v(hk, "ci_upper") - v(hk, "ci_lower"),
            v(wald, "ci_upper") - v(wald, "ci_lower"))
})

test_that("the REML estimate of tau2 agrees with metafor to five significant digits", {
  # Pinned from metafor::rma(yi, vi, method = "REML", control =
  # list(threshold = 1e-12)), its Fisher scoring run to convergence
  # (metafor 5.0.1, 2026-11-09). With optimize() at its default tolerance
  # the two estimates were 0.0279437 and 0.3438417, off in the fourth and
  # the fifth significant digit.
  f <- .het_yivi()
  expect_equal(v(meta_es(f$yi, f$vi), "tau2"), 0.02793870142981249,
               tolerance = 1e-5)
  set.seed(20261109)
  vi <- runif(8, 0.01, 0.2)
  yi <- rnorm(8, 0.3, sqrt(vi + 0.1))
  expect_equal(v(meta_es(yi, vi), "tau2"), 0.3438464137834456,
               tolerance = 1e-5)
})

test_that("the REML estimate of tau2 does not depend on the units of the effect sizes", {
  # optimize() stops at an absolute step in its argument, so a search on
  # tau2 itself stopped 5e-9 (earlier 3e-5) from the maximum whatever the
  # units, which is the size of tau2 for effect sizes in small units. The
  # search is on tau2 over the scale of the data, and agrees across units
  # to the flatness of the likelihood at its maximum.
  f <- .het_yivi()
  ref <- v(meta_es(f$yi, f$vi), "tau2")
  for (unit in c(1e-3, 1e-6, 1e3)) {
    scaled <- meta_es(f$yi * unit, f$vi * unit^2)
    expect_equal(v(scaled, "tau2"), ref * unit^2, tolerance = 1e-6)
    expect_equal(v(scaled, "estimate"),
                 v(meta_es(f$yi, f$vi), "estimate") * unit, tolerance = 1e-6)
    expect_equal(v(scaled, "I2"), v(meta_es(f$yi, f$vi), "I2"),
                 tolerance = 1e-6)
  }
  # Ten risk differences near .02 with sampling variances near 1e-5. The
  # earlier search returned tau2 = 0 here. Pinned from metafor::rma(yi, vi,
  # method = "REML", control = list(threshold = 1e-12)) (metafor 5.0.1,
  # 2026-11-09).
  set.seed(2718)
  n_1 <- round(runif(10, 4e3, 4e4))
  p <- 0.1
  yi <- rnorm(10, 0.02, sqrt(2 * p * (1 - p) / n_1 + 0.003^2))
  vi <- 2 * p * (1 - p) / n_1
  expect_equal(v(meta_es(yi, vi), "tau2"), 9.07549e-06, tolerance = 1e-5)
})

# Twelve studies, eleven with sampling variances near 1e-3 and one with a
# sampling variance of 1e6, so that the scale of the data, max(v) + var(y),
# is about 1e8 times tau2.
.outlier_yivi <- function() {
  list(yi = c(0.0982, -0.23, -0.123, -0.00669, -0.13, 0.0172, -0.00542,
              -0.168, 0.0938, -0.15, 0.044, 472),
       vi = c(0.00057, 0.00095, 0.00048, 0.00048, 0.00095, 0.00051,
              0.00085, 0.00026, 0.00094, 0.00032, 4e-04, 1e6))
}

# The REML estimate as the root of the restricted score equation,
# -sum(w) + sum(w^2) / sum(w) + sum(w^2 (y - mu)^2) = 0 with
# w = 1 / (v + tau2), found by plain bisection on tau2 until the bracket is
# two adjacent floating-point numbers: a different computation from the
# maximization of the restricted likelihood in .meta_tau2_reml(), with no
# tolerance, in the units of the effect sizes.
.reml_score_bisect <- function(yi, vi) {
  score <- function(t) {
    w <- 1 / (vi + t)
    mu <- sum(w * yi) / sum(w)
    -sum(w) + sum(w^2) / sum(w) + sum(w^2 * (yi - mu)^2)
  }
  if (score(0) <= 0) return(0)
  lo <- 0
  hi <- max(vi)
  while (score(hi) > 0) hi <- 2 * hi
  repeat {
    mid <- (lo + hi) / 2
    if (mid <= lo || mid >= hi) break
    if (score(mid) > 0) lo <- mid else hi <- mid
  }
  (lo + hi) / 2
}

test_that("the REML estimate of tau2 is accurate relative to itself when one sampling variance dwarfs the rest", {
  # optimize() once stopped at a distance fixed relative to the scale of the
  # data, max(v) + var(y), about 1e8 times tau2 here, and the estimate came
  # back as 0.01051 where the root of the score equation is 0.01206, in
  # every unit (an error of 13 percent). The stopping distance is now
  # relative to the estimate. The remaining difference is the flatness of
  # the restricted likelihood at its maximum.
  f <- .outlier_yivi()
  for (unit in c(1e-6, 1, 1e3)) {
    yi <- f$yi * unit
    vi <- f$vi * unit^2
    expect_equal(v(meta_es(yi, vi, method = "reml"), "tau2"),
                 .reml_score_bisect(yi, vi), tolerance = 1e-6,
                 label = paste("unit", unit))
  }
  # And on the heterogeneous fixture, where tau2 is near the scale of the
  # data, the same oracle.
  h <- .het_yivi()
  expect_equal(v(meta_es(h$yi, h$vi), "tau2"), .reml_score_bisect(h$yi, h$vi),
               tolerance = 1e-7)
})

# Independent oracles for the Paule-Mandel estimate and the Q-profile limits
# (finding F06). The generalized Q is written out here, and its root in tau2
# is found by plain bisection, halving until the bracket is two adjacent
# floating-point numbers, from a bracket whose upper end is doubled from the
# largest sampling variance until Q falls below the target. No tolerance
# enters, so the oracle is the root in whatever units the effect sizes are.
.q_gen <- function(yi, vi, tau2) {
  w <- 1 / (vi + tau2)
  sum(w * (yi - sum(w * yi) / sum(w))^2)
}
.q_root_bisect <- function(yi, vi, target) {
  if (.q_gen(yi, vi, 0) <= target) return(0)
  lo <- 0
  hi <- max(vi)
  while (.q_gen(yi, vi, hi) > target) hi <- 2 * hi
  repeat {
    mid <- (lo + hi) / 2
    if (mid <= lo || mid >= hi) break
    if (.q_gen(yi, vi, mid) > target) lo <- mid else hi <- mid
  }
  (lo + hi) / 2
}
.pm_oracle <- function(yi, vi, conf_level = 0.95) {
  k <- length(yi)
  a <- 1 - conf_level
  tau2 <- .q_root_bisect(yi, vi, k - 1)
  lims <- c(.q_root_bisect(yi, vi, qchisq(1 - a / 2, k - 1)),
            .q_root_bisect(yi, vi, qchisq(a / 2, k - 1)))
  w  <- 1 / (vi + tau2)
  mu <- sum(w * yi) / sum(w)
  se <- sqrt(1 / sum(w))
  half <- qt(1 - a / 2, k - 2) * sqrt(tau2 + se^2)
  wf <- 1 / vi
  s2 <- (k - 1) * sum(wf) / (sum(wf)^2 - sum(wf^2))
  list(tau2 = tau2, estimate = mu, se = se,
       prediction = mu + c(-1, 1) * half,
       tau2_ci = lims, I2_ci = 100 * lims / (lims + s2))
}
.pm_fixtures <- function() {
  list(
    heterogeneous = .het_yivi(),
    # Sampling variances over three orders of magnitude, with the lower
    # limit of the tau2 interval above 0.
    wide = list(yi = c(.1, .9, -.4, 1.2, .3, -1),
                vi = c(.001, .01, .05, .2, .5, 1)),
    # Q at tau2 = 0 between the lower chi square quantile and k - 1: the
    # PM estimate and the lower limit are 0, the upper limit is not.
    pm_at_zero = list(yi = c(.1, .3, .2, .4, .25), vi = rep(.04, 5)),
    # Q at tau2 = 0 below the lower quantile: the interval is (0, 0).
    all_at_zero = list(yi = c(.30, .32, .29, .31, .33),
                       vi = c(.04, .05, .03, .06, .04)),
    # One sampling variance more than 1e9 times the others. With the
    # stopping distance of the root search fixed relative to the scale of
    # the data, the PM estimate was off by 1.8e-5 relative to itself here
    # (0.01119676 where the root is 0.01119655) and the limits by 9.3e-6
    # and 2.7e-6.
    outlier_variance = .outlier_yivi())
}

test_that("the Paule-Mandel estimate and the Q-profile limits do not depend on the units of the effect sizes (equal variances)", {
  skip_on_cran()
  # With equal sampling variances v the weights do not depend on tau2, Q at
  # tau2 is SS / (v + tau2) with SS the sum of squares of the effect sizes,
  # and each root is SS / target - v. For y = (-1, 0, 1) and v = .01 the PM
  # estimate is 0.99 and the interval 0.2611 to 39.49. The earlier search,
  # made on tau2 itself with an absolute tolerance, gave 10.1 and the
  # interval (101, 101) in units of 1e-6 (after rescaling), and a pooled
  # standard error of 1.836 in place of 0.5774.
  k  <- 3
  ss <- 2
  v0 <- .01
  tau2 <- ss / (k - 1) - v0
  lims <- ss / qchisq(c(.975, .025), k - 1) - v0
  se   <- sqrt((v0 + tau2) / k)
  half <- qt(.975, k - 2) * sqrt(tau2 + se^2)
  for (unit in c(1e-6, 1e-3, 1, 1e3)) {
    lab <- paste("unit", unit)
    fit <- meta_es(c(-1, 0, 1) * unit, rep(v0, 3) * unit^2, method = "pm",
                   hartung_knapp = FALSE)
    expect_equal(v(fit, "tau2") / unit^2, 0.99, tolerance = 1e-10,
                 label = lab)
    expect_equal(v(fit, "se") / unit, se, tolerance = 1e-10, label = lab)
    expect_equal(c(v(fit, "prediction_lower"),
                   v(fit, "prediction_upper")) / unit,
                 c(-half, half), tolerance = 1e-10, label = lab)
    expect_equal(c(v(fit, "tau2_ci_lower"), v(fit, "tau2_ci_upper")) / unit^2,
                 lims, tolerance = 1e-10, label = lab)
    expect_equal(c(v(fit, "I2_ci_lower"), v(fit, "I2_ci_upper")),
                 100 * lims / (lims + v0), tolerance = 1e-10, label = lab)
    # At the PM root Q equals k - 1, so the Hartung-Knapp q is 1 and the
    # adjusted standard error is the conventional one.
    hk <- meta_es(c(-1, 0, 1) * unit, rep(v0, 3) * unit^2, method = "pm")
    expect_equal(attr(hk, "hartung_knapp_q"), 1, tolerance = 1e-10,
                 label = lab)
    expect_equal(v(hk, "se") / unit, se, tolerance = 1e-10, label = lab)
  }
})

test_that("the Paule-Mandel estimate and the Q-profile limits agree with a bisection oracle in any units, unequal variances and boundaries included", {
  skip_on_cran()  # numerical oracle check, run in the full local suite
  for (nm in names(.pm_fixtures())) {
    f <- .pm_fixtures()[[nm]]
    for (unit in c(1e-6, 1e-3, 1, 1e3)) {
      lab <- paste(nm, "unit", unit)
      yi <- f$yi * unit
      vi <- f$vi * unit^2
      or <- .pm_oracle(yi, vi)
      fit <- meta_es(yi, vi, method = "pm", hartung_knapp = FALSE)
      expect_equal(v(fit, "tau2"), or$tau2, tolerance = 1e-9, label = lab)
      expect_equal(v(fit, "estimate"), or$estimate, tolerance = 1e-9,
                   label = lab)
      expect_equal(v(fit, "se"), or$se, tolerance = 1e-9, label = lab)
      expect_equal(c(v(fit, "prediction_lower"), v(fit, "prediction_upper")),
                   or$prediction, tolerance = 1e-9, label = lab)
      expect_equal(c(v(fit, "tau2_ci_lower"), v(fit, "tau2_ci_upper")),
                   or$tau2_ci, tolerance = 1e-9, label = lab)
      expect_equal(c(v(fit, "I2_ci_lower"), v(fit, "I2_ci_upper")), or$I2_ci,
                   tolerance = 1e-9, label = lab)
      # The Q-profile limits are one computation for every estimator.
      for (m in c("reml", "dl")) {
        other <- meta_es(yi, vi, method = m, hartung_knapp = FALSE)
        expect_identical(c(v(other, "tau2_ci_lower"), v(other, "tau2_ci_upper")),
                         c(v(fit, "tau2_ci_lower"), v(fit, "tau2_ci_upper")),
                         label = paste(lab, m))
      }
    }
  }
  # The boundaries fall where the fixtures place them, in every unit.
  for (unit in c(1e-6, 1, 1e3)) {
    f <- .pm_fixtures()$pm_at_zero
    fit <- meta_es(f$yi * unit, f$vi * unit^2, method = "pm")
    expect_identical(c(v(fit, "tau2"), v(fit, "tau2_ci_lower")), c(0, 0))
    expect_gt(v(fit, "tau2_ci_upper"), 0)
    f <- .pm_fixtures()$all_at_zero
    fit <- meta_es(f$yi * unit, f$vi * unit^2, method = "pm",
                   hartung_knapp = FALSE)
    expect_identical(c(v(fit, "tau2"), v(fit, "tau2_ci_lower"),
                       v(fit, "tau2_ci_upper"), v(fit, "I2_ci_upper")),
                     c(0, 0, 0, 0))
    f <- .pm_fixtures()$wide
    expect_gt(v(meta_es(f$yi * unit, f$vi * unit^2, method = "pm"),
                "tau2_ci_lower"), 0)
  }
})

test_that("the generalized Q equals its target at the Paule-Mandel estimate and at each Q-profile limit", {
  skip_on_cran()
  for (f in .pm_fixtures()) {
    k <- length(f$yi)
    for (unit in c(1e-6, 1, 1e3)) {
      yi <- f$yi * unit
      vi <- f$vi * unit^2
      for (cl in c(.90, .95)) {
        fit <- meta_es(yi, vi, method = "pm", conf_level = cl)
        targets <- c(tau2 = k - 1,
                     tau2_ci_lower = qchisq(1 - (1 - cl) / 2, k - 1),
                     tau2_ci_upper = qchisq((1 - cl) / 2, k - 1))
        for (nm in names(targets)) {
          t2 <- v(fit, nm)
          # A limit or estimate of 0 means Q at 0 is already at or below
          # the target; otherwise Q at the value is the target.
          if (t2 == 0) {
            expect_lte(.q_gen(yi, vi, 0), targets[[nm]])
          } else {
            expect_equal(.q_gen(yi, vi, t2), targets[[nm]], tolerance = 1e-9,
                         label = paste(nm, unit, cl))
          }
        }
      }
    }
  }
})

test_that("a generalized Q at tau2 = 0 within rounding error of its target gives a value near 0, not an error", {
  skip_on_cran()  # numerical oracle check, run in the full local suite
  # The boundary test (Q at 0 against the target) was made on the data as
  # given and the root search on the data over their scale; when Q at 0 was
  # within a few units of rounding error of the target the two could
  # disagree in sign, and uniroot() stopped with "f() values at end points
  # not of opposite sign" (an unclassed error) for every estimator. Both are
  # now made on the same scaled Q. These inputs, found by an independent
  # review, stopped; each is written to 17 significant digits so that it is
  # the same double.
  inputs <- list(
    # The upper limit of the 90 percent interval; the default call stopped.
    list(yi = c(0.45159433562859058, 0.47355991830297534),
         vi = c(0.12139033597351338, 0.0013130322236605126), cl = 0.90,
         at = "tau2_ci_upper"),
    # The PM estimate, k = 2.
    list(yi = c(-0.61121366057262982, -0.35713048494116501),
         vi = c(0.058532045085625671, 0.006026215053344104), cl = 0.95,
         at = "tau2"),
    # The PM estimate, k = 3, effect sizes in units near 1e-6.
    list(yi = c(1.5198822018657498e-06, 4.9521492748799948e-07,
                6.8394261886765083e-07),
         vi = c(5.8619691093186751e-13, 8.3735597050776827e-15,
                1.2150710745890525e-13), cl = 0.95, at = "tau2"))
  for (b in inputs) {
    k <- length(b$yi)
    s <- max(b$vi) + var(b$yi)
    a <- 1 - b$cl
    or <- c(tau2 = .q_root_bisect(b$yi, b$vi, k - 1),
            tau2_ci_lower = .q_root_bisect(b$yi, b$vi, qchisq(1 - a / 2, k - 1)),
            tau2_ci_upper = .q_root_bisect(b$yi, b$vi, qchisq(a / 2, k - 1)))
    for (m in c("reml", "dl", "pm")) {
      fit <- expect_no_error(meta_es(b$yi, b$vi, method = m,
                                     conf_level = b$cl))
      nms <- c("tau2_ci_lower", "tau2_ci_upper", if (m == "pm") "tau2")
      for (nm in nms) {
        # Each value agrees with the bisection oracle; the one on the
        # boundary is within rounding error of 0 on the scale of the data.
        got <- v(fit, nm)
        expect_true(abs(got - or[[nm]]) <= 1e-10 * s + 1e-9 * or[[nm]],
                    label = paste(m, nm))
        if (nm == b$at) {
          expect_true(got >= 0 && got <= 1e-10 * s, label = paste(m, nm))
        }
      }
    }
  }
  # Sets placed on each boundary: the deviations from the weighted mean are
  # rescaled so that Q at 0 equals k - 1 or a chi square quantile, then
  # nudged by a few units of relative rounding error, in three units. Every
  # call returns, the value is within rounding error of 0 on the scale of
  # the data and of the bisection oracle, and a positive value puts Q on
  # its target.
  eps <- .Machine$double.eps
  for (f in .pm_fixtures()[c("heterogeneous", "wide")]) {
    k <- length(f$yi)
    w <- 1 / f$vi
    mu <- sum(w * f$yi) / sum(w)
    Q0 <- sum(w * (f$yi - mu)^2)
    targets <- c(tau2 = k - 1, tau2_ci_lower = qchisq(.975, k - 1),
                 tau2_ci_upper = qchisq(.025, k - 1))
    for (nm in names(targets)) {
      for (j in c(-4, -2, -1, 0, 1, 2, 4)) {
        for (unit in c(1e-6, 1, 1e3)) {
          yi <- (mu + sqrt(targets[[nm]] / Q0) * (1 + j * eps) *
                   (f$yi - mu)) * unit
          vi <- f$vi * unit^2
          s <- max(vi) + var(yi)
          lab <- paste(nm, j, unit)
          fit <- expect_no_error(meta_es(yi, vi, method = "pm",
                                         hartung_knapp = FALSE))
          t2 <- v(fit, nm)
          expect_true(t2 >= 0 && t2 <= 1e-10 * s, label = lab)
          expect_true(abs(t2 - .q_root_bisect(yi, vi, targets[[nm]])) <=
                        1e-10 * s, label = lab)
          if (t2 > 0) {
            expect_equal(.q_gen(yi, vi, t2), targets[[nm]], tolerance = 1e-9,
                         label = lab)
          }
        }
      }
    }
  }
})

test_that("the Paule-Mandel estimate and the Q-profile interval agree with metafor run to convergence", {
  # Pinned from metafor::rma(yi, vi, method = "PM", control = list(tol =
  # 1e-14)) and confint() of the same fit (control = list(tol = 1e-14))
  # (metafor 5.0.1, 2026-11-11); live comparison in
  # tools/oracle_snippets/metafor.R.
  pinned <- list(
    heterogeneous = c(tau2 = 0.04310632185035629,
                      estimate = 0.4129886282193178,
                      se = 0.08285465366310646, I2 = 55.57065485525343,
                      tau2_ci_lower = 0, tau2_ci_upper = 0.2262306486670355,
                      I2_ci_lower = 0, I2_ci_upper = 86.77994414065644),
    wide = c(tau2 = 0.3565975544087453, estimate = 0.301238149150264,
             se = 0.292542276336004, I2 = 94.3060300340828,
             tau2_ci_lower = 0.07858957121971007,
             tau2_ci_upper = 3.522839637964757, I2_ci_lower = 78.495322717237,
             I2_ci_upper = 99.39254375539961))
  for (nm in names(pinned)) {
    f <- .pm_fixtures()[[nm]]
    fit <- meta_es(f$yi, f$vi, method = "pm", hartung_knapp = FALSE)
    for (q in names(pinned[[nm]])) {
      expect_equal(v(fit, q), pinned[[nm]][[q]], tolerance = 1e-9,
                   label = paste(nm, q))
    }
  }
})

test_that("effect sizes that differ in a late digit do not raise the zero SE warning", {
  expect_no_warning(meta_es(c(0.5, 0.5, 0.5 + 1e-13), rep(0.04, 3)))
  expect_warning(meta_es(c(0.5, 0.5, 0.5), rep(0.04, 3)),
                 class = "dmar_meta_hk_zero_se")
})

test_that("meta_es() refuses a non-finite effect size or variance under every estimator, the fixed effect included", {
  for (m in c("reml", "pm", "dl", "fe")) {
    expect_error(meta_es(c(0.1, 0.2, 0.3), c(0.1, Inf, 0.1), method = m),
                 "Supply a finite sampling variance")
    expect_error(meta_es(c(0.1, Inf, 0.3), rep(0.1, 3), method = m),
                 "Supply finite effect sizes")
  }
})

# The warning for identical effect sizes says what the conventional
# analysis gives, in its text and in its fields (finding F12).
test_that("the zero standard error warning states and carries the conventional standard error and interval", {
  skip_on_cran()
  w <- expect_warning(meta_es(rep(.5, 4), rep(.04, 4)),
                      class = "dmar_meta_hk_zero_se")
  expect_identical(w$n_studies, 4L)
  expect_equal(w$estimate, 0.5)
  expect_equal(w$se_conventional, 0.1)
  expect_equal(c(w$lower_conventional, w$upper_conventional),
               0.5 + c(-1, 1) * qnorm(0.975) * 0.1)
  expect_equal(c(w$lower_conventional, w$upper_conventional),
               c(0.3040036, 0.6959964), tolerance = 1e-7)
  # The fields are the interval hartung_knapp = FALSE reports.
  plain <- meta_es(rep(.5, 4), rep(.04, 4), hartung_knapp = FALSE)
  expect_equal(w$se_conventional, v(plain, "se"))
  expect_equal(c(w$lower_conventional, w$upper_conventional),
               c(v(plain, "ci_lower"), v(plain, "ci_upper")))
  msg <- conditionMessage(w)
  expect_match(msg, "The 4 effect sizes are identical", fixed = TRUE)
  expect_match(msg, "the square root of q", fixed = TRUE)
  expect_match(msg, "has no width (0.5, 0.5), nor has the prediction interval",
               fixed = TRUE)
  expect_match(msg, "the t statistic is infinite.", fixed = TRUE)
  expect_match(msg, "This is not evidence that the effect is known exactly",
               fixed = TRUE)
  expect_match(msg,
               "is 0.1, with the 95% normal interval 0.304 to 0.696.",
               fixed = TRUE)
  expect_match(msg, "set hartung_knapp = FALSE", fixed = TRUE)

  # The level of the normal interval is the level of the call.
  w90 <- expect_warning(meta_es(rep(.5, 4), rep(.04, 4), conf_level = 0.90),
                        class = "dmar_meta_hk_zero_se")
  expect_equal(c(w90$lower_conventional, w90$upper_conventional),
               0.5 + c(-1, 1) * qnorm(0.95) * 0.1)
  expect_match(conditionMessage(w90), "with the 90% normal interval",
               fixed = TRUE)

  # A common effect size of zero: the statistic is 0 / 0.
  w0 <- expect_warning(meta_es(rep(0, 4), rep(.04, 4)),
                       class = "dmar_meta_hk_zero_se")
  expect_match(conditionMessage(w0),
               "the t statistic is undefined, since the common effect size is 0",
               fixed = TRUE)
  expect_equal(c(w0$lower_conventional, w0$upper_conventional),
               c(-1, 1) * qnorm(0.975) * 0.1)

  # Two studies: no prediction interval is reported, and none is named.
  w2 <- expect_warning(meta_es(c(.5, .5), c(.04, .09)),
                       class = "dmar_meta_hk_zero_se")
  expect_identical(w2$n_studies, 2L)
  expect_false(grepl("prediction interval", conditionMessage(w2),
                     fixed = TRUE))
  expect_equal(w2$se_conventional, sqrt(1 / (1 / .04 + 1 / .09)))

  # Effect sizes equal only up to rounding error: the standard error is
  # 2e-17, not 0, and the statistic is finite.
  wr <- expect_warning(meta_es(c(0.1 + 0.2, 0.3, 0.3), rep(.04, 3)),
                       class = "dmar_meta_hk_zero_se")
  expect_match(conditionMessage(wr), "is the estimate divided by a rounding error",
               fixed = TRUE)
  expect_equal(wr$se_conventional, sqrt(.04 / 3))
})

test_that("meta_es() records q and the conventional standard error whenever the adjustment is applied", {
  skip_on_cran()
  yi <- c(.30, .32, .29, .31, .33)
  vi <- c(.04, .05, .03, .06, .04)
  # q is far below 1 here and no warning is given: a small q is the method
  # working as defined.
  expect_no_warning(fit <- meta_es(yi, vi))
  expect_identical(nrow(fit), 20L)
  expect_equal(attr(fit, "hartung_knapp_q"), 0.006860532, tolerance = 1e-7)
  expect_equal(attr(fit, "se_conventional"), 0.0912870929, tolerance = 1e-9)
  # The adjusted standard error is sqrt(q) times the conventional one.
  expect_equal(v(fit, "se"),
               sqrt(attr(fit, "hartung_knapp_q")) * attr(fit, "se_conventional"),
               tolerance = 1e-12)
  expect_equal(v(fit, "se"), 0.007561158, tolerance = 1e-6)
  # The conventional standard error is the one hartung_knapp = FALSE reports.
  plain <- meta_es(yi, vi, hartung_knapp = FALSE)
  expect_equal(attr(fit, "se_conventional"), v(plain, "se"))
  expect_null(attr(plain, "hartung_knapp_q"))
  expect_null(attr(plain, "se_conventional"))
  expect_null(attr(meta_es(yi, vi, method = "fe"), "hartung_knapp_q"))
  # With heterogeneity the weights include tau^2.
  f <- .te_yivi()
  te <- meta_es(f$yi, f$vi, method = "dl")
  w <- 1 / (f$vi + v(te, "tau2"))
  mu <- sum(w * f$yi) / sum(w)
  expect_equal(attr(te, "hartung_knapp_q"),
               sum(w * (f$yi - mu)^2) / (length(f$yi) - 1))
  expect_equal(attr(te, "se_conventional"), sqrt(1 / sum(w)))
  # Identical effect sizes: q is 0.
  z <- suppressWarnings(meta_es(rep(.5, 4), rep(.04, 4)))
  expect_identical(attr(z, "hartung_knapp_q"), 0)
  expect_equal(attr(z, "se_conventional"), 0.1)
  expect_identical(nrow(z), 20L)
})

test_that("the adjusted interval is narrower than the normal one exactly when q is below (z / t)^2", {
  skip_on_cran()
  # The values ?meta_es states for a 95 percent interval.
  k <- c(2, 3, 5, 10)
  expect_equal(round((qnorm(.975) / qt(.975, k - 1))^2, 2),
               c(0.02, 0.21, 0.50, 0.75))
  set.seed(1109)
  for (i in 1:40) {
    k_i <- sample(2:10, 1)
    yi <- rnorm(k_i, 0.3, 0.2)
    vi <- runif(k_i, 0.02, 0.08)
    hk <- meta_es(yi, vi, method = "dl")
    plain <- meta_es(yi, vi, method = "dl", hartung_knapp = FALSE)
    narrower <- (v(hk, "ci_upper") - v(hk, "ci_lower")) <
      (v(plain, "ci_upper") - v(plain, "ci_lower"))
    expect_identical(narrower,
                     attr(hk, "hartung_knapp_q") <
                       (qnorm(.975) / qt(.975, k_i - 1))^2)
  }
})

test_that("the test statistic is the row test_statistic, with its reference distribution in the attributes", {
  skip_on_cran()
  # DMAR 1.0.x named the row t under the Hartung-Knapp adjustment and z
  # without it, so the name changed with an argument; under the NEWS rule
  # for such statistics it is test_statistic, and the reference is
  # recorded in test_distribution and test_df and printed under the table.
  # The statistic is estimate / se either way, and the p-value is its
  # two-sided tail area under the recorded reference.
  f <- .te_yivi()
  k <- length(f$yi)
  hk <- meta_es(f$yi, f$vi)
  no <- meta_es(f$yi, f$vi, hartung_knapp = FALSE)
  fe <- meta_es(f$yi, f$vi, method = "fe")
  for (x in list(hk, no, fe)) {
    expect_s3_class(x, c("dmar_meta", "dmar_tbl", "data.frame"), exact = TRUE)
    expect_identical(x$term[3], "test_statistic")
    expect_false(any(c("t", "z") %in% x$term))
    expect_identical(v(x, "test_statistic"), v(x, "estimate") / v(x, "se"))
  }
  expect_identical(attr(hk, "test_distribution"), "t")
  expect_identical(attr(hk, "test_df"), k - 1)
  expect_equal(v(hk, "p_value"),
               2 * pt(-abs(v(hk, "test_statistic")), df = k - 1),
               tolerance = 1e-14)
  for (x in list(no, fe)) {
    expect_identical(attr(x, "test_distribution"), "normal")
    expect_identical(attr(x, "test_df"), NA_real_)
    expect_equal(v(x, "p_value"), 2 * pnorm(-abs(v(x, "test_statistic"))),
                 tolerance = 1e-14)
  }
  expect_output(print(hk), sprintf(
    "test_statistic: referred to Student's t on %d df (Hartung-Knapp adjustment)",
    k - 1), fixed = TRUE)
  expect_output(print(no), "test_statistic: referred to the standard normal",
                fixed = TRUE)
  # tidy() keeps the row, and meta_r() and meta_smd() return the same table.
  expect_true("test_statistic" %in% generics::tidy(hk)$term)
  r <- meta_r(c(.28, .35, .22, .40, .31), c(120, 85, 200, 60, 150))
  expect_identical(r$term, hk$term)
  expect_identical(attr(r, "test_distribution"), "t")
  s <- meta_smd(smd = teacher_expectancy$d,
                n_1 = teacher_expectancy$n_experimental,
                n_2 = teacher_expectancy$n_control, hartung_knapp = FALSE)
  expect_identical(s$term, hk$term)
  expect_identical(attr(s, "test_distribution"), "normal")
})
