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
  expect_true(all(c("estimate", "se", "t", "p_value", "lower_limit",
                    "upper_limit", "prediction_lower", "prediction_upper",
                    "tau2", "tau2_lower", "tau2_upper", "tau", "I2",
                    "I2_lower", "I2_upper", "H2",
                    "Q", "Q_df", "Q_p", "k") %in% res$term))
  # Only tau2 and I2 carry interval rows; ?meta_es documents the table
  # accordingly. If tau or H2 limits are ever added, update the help page
  # alongside this expectation.
  expect_false(any(c("tau_lower", "tau_upper", "H2_lower", "H2_upper")
                   %in% res$term))
  expect_identical(attr(res, "method"), "reml")
  expect_true(attr(res, "hartung_knapp"))
  expect_identical(attr(res, "conf_level"), 0.95)
  expect_equal(v(res, "k"), 19)
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
  pinned <- list(
    reml = list(estimate = 0.4164518538709743, se = 0.0742598072604618,
                tau2 = 0.02794355697336737, Q = 20.61130014356519,
                I2 = 44.77593674399343),
    dl   = list(estimate = 0.4158314259576409, se = 0.07556718741348845,
                tau2 = 0.03011305513068721, Q = 20.61130014356519,
                I2 = 46.63121722850569),
    pm   = list(estimate = 0.4129887431681661, se = 0.08285430501875582,
                tau2 = 0.04310566448959827, Q = 20.61130014356519,
                I2 = 55.57027834057255),
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

test_that("meta_es() Hartung-Knapp matches metafor test='knha'", {
  f <- .het_yivi()
  ours <- meta_es(f$yi, f$vi, method = "reml", hartung_knapp = TRUE)
  # Pinned from metafor::rma(yi, vi, method = "REML", test = "knha")
  # (metafor 5.0.1, 2026-08-09); live comparison in tools/oracle_checks.R.
  expect_equal(v(ours, "estimate"), 0.4164518538709743, tolerance = 1e-6)
  expect_equal(v(ours, "se"), 0.08099591175605304, tolerance = 1e-6)
  expect_equal(v(ours, "t"), 5.141640421621054, tolerance = 1e-6)
  expect_equal(v(ours, "p_value"), 0.0003224106762200866, tolerance = 1e-5)
  expect_equal(v(ours, "lower_limit"), 0.2381810540678096,
               tolerance = 1e-6)
  expect_equal(v(ours, "upper_limit"), 0.5947226536741389,
               tolerance = 1e-6)
})

test_that("meta_es() tau2 Q-profile interval and prediction interval match metafor", {
  f <- .het_yivi()
  ours <- meta_es(f$yi, f$vi, method = "reml", hartung_knapp = FALSE)
  # Pinned from confint(metafor::rma(yi, vi, method = "REML")), the
  # Q-profile interval for tau2 (metafor 5.0.1, 2026-08-09); live
  # comparison in tools/oracle_checks.R.
  expect_equal(v(ours, "tau2_lower"), 0, tolerance = 1e-4)
  expect_equal(v(ours, "tau2_upper"), 0.2262314670033073, tolerance = 1e-4)
  # DMAR's prediction interval follows Higgins, Thompson, and Spiegelhalter
  # (2009): t on k - 2 df. metafor's default uses z, so the pinned values
  # reconstruct the HTS interval beta -/+ qt(.975, k - 2) *
  # sqrt(tau2 + se^2) from metafor's own tau2 and SE.
  expect_equal(v(ours, "prediction_lower"),
               0.008890747225550666, tolerance = 1e-4)
  expect_equal(v(ours, "prediction_upper"),
               0.8240129605163979, tolerance = 1e-4)
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
  expect_error(meta_es(c(1, 2), c(0.1, -0.1)), "positive sampling variance")
  expect_error(meta_es(1, 0.1), "two or more")
  expect_error(meta_es(c(1, 2), c(0.1, 0.1), conf_level = 1.2),
               "\\(0, 1\\)")
  expect_error(meta_es(c(1, 2), c(0.1, 0.1), hartung_knapp = NA),
               "TRUE or FALSE")
})
