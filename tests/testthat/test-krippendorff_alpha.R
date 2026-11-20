# The bootstraps in this file are mechanics tests run with a small B; the
# rule on the number of replications is tested in test-dmar_bootstrap.R.
local_small_B()

test_that("krippendorff_alpha() returns documented rows", {
  ratings <- matrix(c(
    1, 2, 3, 3, 2, 1, 4, 1, 2, NA, NA, NA,
    1, 2, 3, 3, 2, 2, 4, 1, 2, 5,  NA, 3,
    NA, 3, 3, 3, 2, 3, 4, 2, 2, 5,  1,  NA,
    1, 2, 3, 3, 2, 4, 4, 1, 2, 5,  1,  NA
  ), nrow = 12, ncol = 4)
  res <- krippendorff_alpha(ratings, level = "nominal", boot = FALSE)
  expect_identical(res$term,
                   c("krippendorff_alpha", "D_observed", "D_expected",
                     "n_pairable", "ci_lower", "ci_upper", "B_used"))
})

test_that("krippendorff_alpha() is 1.0 on perfect agreement", {
  # All four raters give identical scores on every unit.
  R <- matrix(rep(1:5, times = 4), 5, 4)
  res <- krippendorff_alpha(R, boot = FALSE)
  expect_equal(res$value[res$term == "krippendorff_alpha"], 1,
               tolerance = 1e-10)
})

test_that("krippendorff_alpha() interval level returns finite alpha", {
  set.seed(113)
  r1 <- rnorm(30); r2 <- r1 + rnorm(30, 0, 0.2)
  res <- krippendorff_alpha(cbind(r1, r2), level = "interval",
                            boot = FALSE)
  expect_true(is.finite(res$value[1]))
})

test_that("krippendorff_alpha() bootstrap returns CI bounds", {
  skip_on_cran()  # Additional model fits run in the full local suite.
  set.seed(113)
  R <- matrix(sample(1:4, 40, TRUE), 20, 2)
  res <- krippendorff_alpha(R, boot = TRUE, B = 200L)
  expect_true("ci_lower" %in% res$term)
  expect_true("ci_upper" %in% res$term)
})

test_that("krippendorff_alpha() ratio level handles an identical zero pair (HIGH-05)", {
  # An identical zero pair (0, 0) makes the ratio distance ((a - b) / (a + b))^2
  # evaluate 0 / 0 = NaN, which used to crash the comparison. The oracle is
  # irr's independent implementation.
  # Pinned from irr::kripp.alpha (irr 0.85, 2026-08-09); live comparison in
  # tools/oracle_checks.R.
  ratings <- rbind(c(0, 0), c(1, 1), c(0, 1), c(1, 0))
  oracle  <- 0.125
  res     <- krippendorff_alpha(ratings, level = "ratio", boot = FALSE)
  alpha   <- res$value[res$term == "krippendorff_alpha"]
  expect_false(is.nan(alpha))
  expect_equal(alpha, oracle, tolerance = 1e-8)
})

test_that("the bootstrap is opt-in: no analysis runs one unless requested", {
  set.seed(113)
  r1 <- rnorm(20)
  r2 <- r1 + rnorm(20, 0, 0.4)
  state <- .Random.seed
  res <- krippendorff_alpha(cbind(r1, r2), level = "interval")
  # No resample was drawn, so the interval rows are empty.
  expect_identical(.Random.seed, state)
  lim <- res$value[res$term %in% c("ci_lower", "ci_upper", "B_used")]
  expect_length(lim, 3L)
  expect_true(all(is.na(lim)))
})

test_that("krippendorff_alpha() returns one shape whatever boot", {
  skip_on_cran()  # Additional model fits run in the full local suite.
  # Through 1.0.1 the ci_lower, ci_upper, and B_used rows existed
  # only when boot = TRUE, so a script indexing the table, and tidy() and
  # glance(), saw a different shape by argument (the defect htmt() had).
  # The rows are now always present, NA when no bootstrap is run, and the
  # point estimate and disagreement rows are unchanged.
  set.seed(113)
  r1 <- rnorm(30)
  r2 <- r1 + rnorm(30, 0, 0.3)
  R <- cbind(r1, r2)
  plain <- krippendorff_alpha(R, level = "interval")
  boot <- krippendorff_alpha(R, level = "interval", boot = TRUE, B = 200L,
                             seed = 113)
  expect_identical(plain$term, boot$term)
  expect_identical(names(plain), names(boot))
  expect_identical(plain$value[1:4], boot$value[1:4])
  lim <- c("ci_lower", "ci_upper", "B_used")
  expect_true(all(is.na(plain$value[plain$term %in% lim])))
  expect_true(all(is.finite(boot$value[boot$term %in% lim])))
  expect_identical(boot$value[boot$term == "B_used"], 200)
  expect_lt(boot$value[boot$term == "ci_lower"], boot$value[1])
  expect_gt(boot$value[boot$term == "ci_upper"], boot$value[1])

  # The broom verbs answer with the same columns either way.
  expect_identical(names(generics::tidy(plain)), names(generics::tidy(boot)))
  expect_identical(generics::tidy(plain)$term, generics::tidy(boot)$term)
  expect_identical(names(generics::glance(plain)),
                   names(generics::glance(boot)))
  expect_true(is.na(generics::glance(plain)$B_used))

  # The confidence level is recorded, and printed, only beside an interval
  # that was formed.
  expect_null(attr(plain, "conf_level"))
  expect_identical(attr(boot, "conf_level"), 0.95)
  expect_false(any(grepl("Confidence level", capture.output(print(plain)))))
  expect_true(any(grepl("Confidence level: 95%", capture.output(print(boot)),
                        fixed = TRUE)))

  # No sentence is written from limits that were never formed.
  expect_error(results_sentence(plain), "limits are NA")
  expect_match(results_sentence(boot, label = "alpha"),
               "^alpha = [0-9.]+, 95% CI \\[[0-9.]+, [0-9.]+\\]$")
})

test_that("krippendorff_alpha() validates B and reports the resamples it drops", {
  skip_on_cran()  # Additional model fits run in the full local suite.
  # The stage 2b review: B was not validated (B = 10 gave a percentile
  # interval from ten resamples), and a resample on which alpha could not
  # be computed was dropped without a word. B is now a whole number of at
  # least 100, the minimum of the other agreement bootstraps, and dropped
  # resamples are counted in one warning, as the package's bootstrap rule
  # asks.
  ok <- cbind(c(1, 2, 3, 1, 2, 3), c(1, 2, 3, 1, 3, 3))
  for (bad in list(1, 150.5, c(200, 300), NA_real_, "200")) {
    expect_error(krippendorff_alpha(ok, boot = TRUE, B = bad),
                 "single whole number", class = "dmar_count_invalid")
  }
  # B is not checked when no bootstrap is run.
  expect_no_error(krippendorff_alpha(ok, B = 10))

  # Ten of twelve units are rated 1 by both raters and two are rated 2, so
  # a resample that draws neither of those two has every value alike and
  # no expected disagreement, and alpha is undefined there (about one
  # resample in nine).
  R <- cbind(c(rep(1, 10), 2, 2), c(rep(1, 10), 2, 2))
  expect_equal(krippendorff_alpha(R)$value[1], 1)
  w <- testthat::capture_warnings(
    res <- krippendorff_alpha(R, boot = TRUE, B = 300L, seed = 113))
  expect_length(w, 1L)
  used <- res$value[res$term == "B_used"]
  expect_lt(used, 300)
  expect_gte(used, 100)
  expect_match(w, paste0("^krippendorff_alpha\\(\\): ", 300 - used,
                         " of 300 bootstrap replications \\([0-9.]+%\\) ",
                         "left alpha undefined"))
  expect_match(w, paste0("computed from the ", used, " that remained"),
               fixed = TRUE)
  expect_identical(attr(res, "n_dropped"), 300L - as.integer(used))
  expect_true(all(is.finite(res$value[res$term %in%
                                        c("ci_lower", "ci_upper")])))
  # (DMAR 1.0.x stopped when fewer than 100 resamples remained; the
  # limits are now withheld under a classed warning, tested in
  # test-dmar_bootstrap.R.)

  # A coefficient the data leave undefined gets no bootstrap: the interval
  # rows stay NA, no confidence level is recorded, and no resample is
  # drawn.
  same <- cbind(rep(1, 6), rep(1, 6))
  set.seed(113)
  state <- .Random.seed
  und <- krippendorff_alpha(same, boot = TRUE, B = 200L)
  expect_identical(.Random.seed, state)
  expect_true(all(is.na(und$value)[und$term %in%
                                    c("krippendorff_alpha", "ci_lower",
                                      "ci_upper", "B_used")]))
  expect_null(attr(und, "conf_level"))
})
