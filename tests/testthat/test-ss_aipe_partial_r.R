test_that("ss_aipe_partial_r() plans N = p + ceiling(z^2 V / h^2) on the raw scale", {
  z <- stats::qnorm(.975)
  for (rho in c(.05, .25, .72)) for (p in c(2, 4, 7)) for (w in c(.05, .4)) {
    expect_identical(
      ss_aipe_partial_r(true_cor = rho, p = p, width = w)$value[1],
      p + ceiling(z^2 * (1 - rho^2)^2 / (w / 2)^2))
    # On Fisher's Z scale the plan is never below p + 3.
    expect_gte(ss_aipe_partial_r(true_cor = rho, p = p, width = 1.9,
                                 fisher_z = TRUE)$value[1], p + 3)
  }
})

# On CRAN the file stops here, after the anchor above, which comes first for
# that reason; the blocks below run in the full local suite (see cran_files
# in tests/testthat.R).
skip_on_cran()

test_that("ss_aipe_partial_r() returns documented rows", {
  res <- ss_aipe_partial_r(true_cor = 0.30, p = 3, width = 0.20)
  expect_true(all(c("necessary_N", "expected_width", "true_cor", "p",
                    "width", "conf_level") %in% res$term))
  # The rows carry the values given, and no row has a 1.0.x name.
  expect_identical(res$value[res$term == "true_cor"], 0.30)
  expect_identical(res$value[res$term == "p"], 3)
  expect_false(any(c("rho", "true_rho", "J", "n") %in% res$term))
  n_val <- res$value[res$term == "necessary_N"]
  expect_gt(n_val, 3)   # sanity: more than p = 3
})

test_that("ss_aipe_partial_r() yields expected width <= target", {
  res <- ss_aipe_partial_r(true_cor = 0.30, p = 3, width = 0.20)
  expect_lte(res$value[res$term == "expected_width"], 0.20 + 1e-6)
})

test_that("ss_aipe_partial_r() Fisher's Z route gives a similar sample size", {
  res_raw <- ss_aipe_partial_r(true_cor = 0.30, p = 3, width = 0.20,
                                fisher_z = FALSE)
  res_fz  <- ss_aipe_partial_r(true_cor = 0.30, p = 3, width = 0.20,
                                fisher_z = TRUE)
  # Two routes should be in the same ballpark (within a factor of 2):
  ratio <- res_fz$value[1] / res_raw$value[1]
  expect_gt(ratio, 0.5)
  expect_lt(ratio, 2)
})

# ---- The assurance step ---------------------------------------------------

test_that("ss_aipe_partial_r() with assurance never plans fewer cases, and reports the request and the attained probability", {
  res_50  <- ss_aipe_partial_r(true_cor = 0.30, p = 3, width = 0.20)
  res_80  <- ss_aipe_partial_r(true_cor = 0.30, p = 3, width = 0.20,
                                assurance = 0.80)
  expect_gt(res_80$value[1], res_50$value[1])
  expect_identical(res_50$term, c("necessary_N", "expected_width",
                                  "true_cor", "p",
                                  "width", "conf_level"))
  expect_identical(res_80$term, c(res_50$term, "assurance",
                                  "assurance_attained"))
  expect_equal(res_80$value[res_80$term == "assurance"], 0.80)
  expect_gte(res_80$value[res_80$term == "assurance_attained"], 0.80)
  expect_lte(res_80$value[res_80$term == "expected_width"], 0.20)
  expect_s3_class(res_80, "dmar_tbl")
})

test_that("ss_aipe_partial_r() assurance plans equal the values computed from the distribution of r", {
  N <- function(...) ss_aipe_partial_r(...)$value[1]
  # One variable partialed out (p = 2), width .10, assurance .80. DMAR
  # 1.0.x returned 426 at rho = .70
  # and 1318 at rho = .30.
  expect_equal(N(true_cor = .70, p = 2, width = .10, assurance = .80), 448)
  expect_equal(N(true_cor = .30, p = 2, width = .10, assurance = .80), 1310)
  expect_equal(N(true_cor = .50, p = 2, width = .10, assurance = .80), 915)
  expect_equal(N(true_cor = .85, p = 2, width = .10, assurance = .80), 151)
  expect_equal(N(true_cor = .95, p = 2, width = .10, assurance = .80), 28)
  # Narrow widths, where the noncentrality of the usual form passes 37.62.
  expect_equal(N(true_cor = .95, p = 2, width = .02, assurance = .80), 427)
  expect_equal(N(true_cor = .70, p = 2, width = .04, assurance = .80), 2617)
  # The sign of rho does not matter.
  expect_equal(N(true_cor = -.70, p = 2, width = .10, assurance = .80), 448)
})

test_that("ss_aipe_partial_r() returns the smallest sufficient n at or above the plan without assurance", {
  skip_on_cran()
  z <- stats::qnorm(.975)
  for (fz in c(FALSE, TRUE)) {
    for (rho in c(.1, .3, .7, -.85)) {
      for (J in c(1, 5)) {
        for (a in c(.80, .99)) {
          info <- paste("rho", rho, "J", J, "assurance", a, "fisher_z", fz)
          n_0 <- ss_aipe_partial_r(rho, p = J + 1, width = .3, fisher_z = fz)$value[1]
          res <- ss_aipe_partial_r(rho, p = J + 1, width = .3, fisher_z = fz,
                                   assurance = a)
          n <- res$value[1]
          attained <- res$value[res$term == "assurance_attained"]
          expect_gte(n, n_0, label = info)
          expect_gte(attained, a, label = info)
          expect_equal(attained,
                       DMAR:::.partial_r_width_probability(n, rho, J, .15, z,
                                                           fz),
                       tolerance = 1e-12, info = info)
          if (n > n_0)
            expect_lt(DMAR:::.partial_r_width_probability(n - 1, rho, J, .15,
                                                          z, fz),
                      a, label = info)
        }
      }
    }
  }
})

test_that("ss_aipe_partial_r() adds nothing for assurance at rho = 0, and does not go below the plan without it", {
  # At rho = 0 the interval is widest at r = 0, which the plan without
  # assurance already makes narrow enough, so every sample meets the width.
  for (fz in c(FALSE, TRUE)) {
    n_0 <- ss_aipe_partial_r(true_cor = 0, p = 2, width = .10,
                             fisher_z = fz)$value[1]
    res <- ss_aipe_partial_r(true_cor = 0, p = 2, width = .10, fisher_z = fz,
                             assurance = .80)
    expect_equal(res$value[1], n_0)
    expect_equal(res$value[res$term == "assurance_attained"], 1)
  }
  expect_equal(ss_aipe_partial_r(true_cor = 0, p = 2, width = .10)$value[1], 1539)
  # With a width of .40 the plan is 99. One case below it the probability
  # is .893 (the interval is too wide only when |r| is very near 0), which
  # already exceeds a request of .80, so a search that looked below the
  # plan would return fewer cases than the plan without assurance. The
  # answer stays at the plan.
  expect_equal(ss_aipe_partial_r(true_cor = 0, p = 2, width = .40)$value[1], 99)
  p_below <- DMAR:::.partial_r_width_probability(98, 0, 1, .20,
                                                 stats::qnorm(.975), FALSE)
  expect_equal(p_below, .8934, tolerance = 1e-3)
  expect_equal(ss_aipe_partial_r(true_cor = 0, p = 2, width = .40,
                                 assurance = .80)$value[1], 99)
})

test_that("ss_aipe_partial_r() treats width as a half-width under which_width = 'Lower' in the assurance step too", {
  full  <- ss_aipe_partial_r(true_cor = .5, p = 3, width = .20, assurance = .90)
  lower <- ss_aipe_partial_r(true_cor = .5, p = 3, width = .10, assurance = .90,
                             which_width = "Lower")
  upper <- ss_aipe_partial_r(true_cor = .5, p = 3, width = .10, assurance = .90,
                             which_width = "Upper", fisher_z = FALSE)
  expect_equal(lower$value[1], full$value[1])
  expect_equal(upper$value[1], full$value[1])
  expect_equal(lower$value[lower$term == "assurance_attained"],
               full$value[full$term == "assurance_attained"])
  expect_equal(lower$value[lower$term == "width"], .10)
})

test_that("ss_aipe_partial_r() checks assurance", {
  expect_error(ss_aipe_partial_r(.3, p = 3, .2, assurance = .5), "assurance")
  expect_error(ss_aipe_partial_r(.3, p = 3, .2, assurance = 1), "assurance")
  expect_error(ss_aipe_partial_r(.3, p = 3, .2, assurance = c(.8, .9)),
               "assurance")
  expect_error(ss_aipe_partial_r(.3, p = 3, .2, assurance = "high"), "assurance")
})

test_that(".partial_r_cdf() is the distribution function of a sample correlation", {
  skip_on_cran()
  cdf <- DMAR:::.partial_r_cdf
  # At rho = 0, sqrt(m - 2) r / sqrt(1 - r^2) is Student's t on m - 2 df.
  for (m in c(3, 4, 12, 300)) {
    for (c in c(-.9, -.3, 0, .1, .7)) {
      expect_lt(abs(cdf(c, m, 0) -
                      stats::pt(sqrt(m - 2) * c / sqrt(1 - c^2), m - 2)),
                1e-11, label = paste("difference at m", m, "c", c))
    }
  }
  # The ends, symmetry in the sign of rho, and monotonicity in c.
  expect_identical(cdf(-1, 10, .4), 0)
  expect_identical(cdf(1, 10, .4), 1)
  expect_equal(cdf(.2, 15, .6), 1 - cdf(-.2, 15, -.6), tolerance = 1e-10)
  grid <- seq(-.95, .95, by = .05)
  expect_true(all(diff(vapply(grid, cdf, 0, m = 9, rho = -.5)) > 0))
  # The noncentral t form, P(T <= sqrt(m - 2) a | ncp = theta sqrt(u))
  # integrated over u chi square on m - 1 df, where the noncentrality is
  # moderate and R's pt() is accurate.
  alt <- function(c, m, rho) {
    a <- c / sqrt(1 - c^2)
    theta <- rho / sqrt(1 - rho^2)
    stats::integrate(function(u) vapply(u, function(v)
      stats::pt(sqrt(m - 2) * a, m - 2, ncp = theta * sqrt(v)), 0) *
        stats::dchisq(u, m - 1),
      stats::qchisq(1e-13, m - 1), stats::qchisq(1 - 1e-13, m - 1),
      rel.tol = 1e-10, subdivisions = 2000L)$value
  }
  expect_equal(cdf(.60, 12, .7), alt(.60, 12, .7), tolerance = 1e-7)
  expect_equal(cdf(-.10, 30, -.3), alt(-.10, 30, -.3), tolerance = 1e-7)
  expect_equal(cdf(.45, 60, .5), alt(.45, 60, .5), tolerance = 1e-7)
  # m = 400, rho = .95, where that form is a normal approximation in R: a
  # simulation of 40,000,000 correlations gave .267578 (SE .000070) at
  # c = .947.
  expect_equal(cdf(.947, 400, .95), .267521, tolerance = 2e-6)
  expect_lt(abs(cdf(.947, 400, .95) - .267578), 3 * .000070)
})

test_that("the attained assurance of ss_aipe_partial_r() agrees with simulated regression data", {
  skip_on_cran()
  # Partial correlations from data: Y, X, and J variables partialed out of
  # both by least squares. 20,000 samples; the standard error of a
  # proportion near .8 is .0028.
  check <- function(rho, J, width, a, fz, G = 20000) {
    res <- ss_aipe_partial_r(rho, p = J + 1, width, fisher_z = fz, assurance = a)
    n <- res$value[1]
    z <- stats::qnorm(.975)
    r <- vapply(seq_len(G), function(g) {
      Z <- cbind(1, matrix(stats::rnorm(n * J), n))
      x <- stats::rnorm(n)
      y <- rho * x + sqrt(1 - rho^2) * stats::rnorm(n) + Z[, 2]
      q <- qr(Z)
      stats::cor(qr.resid(q, y), qr.resid(q, x))
    }, numeric(1))
    half <- DMAR:::.partial_r_half_width(r, n, J, z, fz)
    share <- mean(half <= width / 2)
    attained <- res$value[res$term == "assurance_attained"]
    expect_lt(abs(share - attained),
              4 * sqrt(attained * (1 - attained) / G))
  }
  set.seed(20261110)
  check(.70, 2, .25, .80, FALSE)
  check(-.40, 4, .30, .90, TRUE)
})

test_that("ss_aipe_partial_r() rejects bad inputs", {
  expect_error(ss_aipe_partial_r(true_cor = 1.5, p = 3, width = 0.2), "in \\(-1, 1\\)")
  expect_error(ss_aipe_partial_r(true_cor = 1.5, p = 3, width = 0.2),
               "'true_cor' must be a single value in (-1, 1).", fixed = TRUE)
  expect_error(ss_aipe_partial_r(true_cor = 0.3, p = 1, width = 0.2),
               "whole number of at least 2")
  expect_error(ss_aipe_partial_r(true_cor = 0.3, p = 2.5, width = 0.2),
               "whole number of at least 2")
  expect_error(ss_aipe_partial_r(true_cor = 0.3, p = NA_real_, width = 0.2),
               "whole number of at least 2")
  expect_error(ss_aipe_partial_r(true_cor = 0.3, p = 0, width = 0.2),
               "whole number of at least 2")
  expect_error(ss_aipe_partial_r(true_cor = 0.3, width = 0.2),
               "'p' (the number of predictors) is required.", fixed = TRUE)
  expect_error(ss_aipe_partial_r(true_cor = 0.3, p = 3, width = 0),   "positive")
})

# ---- The count in its place, and the 1.0.x names --------------------------

test_that("ss_aipe_partial_r() takes p in its place, by position or by name", {
  named <- ss_aipe_partial_r(true_cor = .3, p = 3, width = .2)
  expect_identical(named$value[1], 322)
  expect_identical(ss_aipe_partial_r(.3, 3, .2), named)
  expect_identical(do.call(ss_aipe_partial_r, list(.3, 3, .2)), named)
  fwd <- function(...) ss_aipe_partial_r(...)
  expect_identical(fwd(.3, 3, .2), named)
})

test_that("ss_aipe_partial_r() does not accept the DMAR 1.0.x names", {
  expect_error(ss_aipe_partial_r(rho = .3, p = 3, width = .2),
               "unused argument")
  expect_error(ss_aipe_partial_r(true_cor = .3, J = 2, width = .2),
               "unused argument")
  expect_identical(names(formals(ss_aipe_partial_r)),
                   c("true_cor", "p", "width", "which_width", "conf_level",
                     "fisher_z", "assurance"))
  expect_error(ss_aipe_partial_r(true_cor = .3, width = .2),
               "'p' (the number of predictors) is required.", fixed = TRUE)
})

test_that("the assurance step on the Fisher scale returns pinned plans", {
  # Each equals an independent calculation (Hotelling's density of r,
  # integrated over the region where the Fisher-scale interval is no wider
  # than the target), which is below the request at N - 1: attained
  # probabilities .90054, .60028, and 1. The half-width solved for must be
  # the Fisher-scale one; with the raw-scale half-width the plans differ.
  plan <- function(...) {
    res <- suppressWarnings(ss_aipe_partial_r(..., fisher_z = TRUE))
    c(res$value[res$term == "necessary_N"],
      res$value[res$term == "assurance_attained"])
  }
  a <- plan(true_cor = .92, p = 2, width = .10, assurance = .90)
  expect_equal(a[1], 70)
  expect_equal(a[2], 0.900544265, tolerance = 1e-8)
  b <- plan(true_cor = -.45, p = 4, width = .18, assurance = .60)
  expect_equal(b[1], 314)
  expect_equal(b[2], 0.6002842655, tolerance = 1e-8)
  cc <- plan(true_cor = .15, p = 9, width = .22, assurance = .99)
  expect_equal(cc[1], 326)
  expect_equal(cc[2], 1, tolerance = 1e-8)
})

test_that("ss_aipe_partial_r() gives no warning, with or without assurance", {
  expect_no_warning(ss_aipe_partial_r(true_cor = .7, p = 2, width = .1, assurance = .8))
  expect_no_warning(ss_aipe_partial_r(true_cor = .7, p = 2, width = .1))
  expect_no_warning(ss_aipe_partial_r(true_cor = .7, p = 2, width = .1,
                                      fisher_z = TRUE))
})
