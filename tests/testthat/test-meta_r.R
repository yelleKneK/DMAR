v <- function(tab, t) tab$value[tab$term == t]

test_that("meta_r() matches metafor's Fisher's Z pipeline, back-transformed", {
  r <- c(.28, .35, .22, .40, .31)
  n <- c(120, 85, 200, 60, 150)
  ours <- meta_r(r, n, hartung_knapp = FALSE)
  # Pinned from metafor::rma(yi = atanh(r), vi = 1 / (n - 3), method = "REML"),
  # back-transformed with tanh (metafor 5.0.1, 2026-08-09); live comparison in
  # tools/oracle_checks.R.
  expect_equal(v(ours, "estimate"), 0.2897158161550648, tolerance = 1e-6)
  expect_equal(v(ours, "lower_limit"), 0.2148407135896157, tolerance = 1e-6)
  expect_equal(v(ours, "upper_limit"), 0.3612051691843595, tolerance = 1e-6)
  expect_equal(v(ours, "tau2"), 0, tolerance = 1e-5)
})

test_that("meta_r() reliability corrections disattenuate before pooling", {
  r <- c(.28, .35, .22, .40, .31)
  n <- c(120, 85, 200, 60, 150)
  plain <- meta_r(r, n, hartung_knapp = FALSE)
  fixed <- meta_r(r, n, reliability_y = 0.80, hartung_knapp = FALSE)
  # Correcting for criterion unreliability raises the pooled validity, and
  # equals pooling the hand-corrected correlations directly.
  expect_gt(v(fixed, "estimate"), v(plain, "estimate"))
  by_hand <- meta_r(r / sqrt(0.80), n, hartung_knapp = FALSE)
  expect_equal(v(fixed, "estimate"), v(by_hand, "estimate"))
  # Per-study reliabilities are honored elementwise.
  rel <- c(.7, .8, .75, .9, .85)
  el <- meta_r(r, n, reliability_x = rel, hartung_knapp = FALSE)
  byh <- meta_r(r / sqrt(rel), n, hartung_knapp = FALSE)
  expect_equal(v(el, "estimate"), v(byh, "estimate"))
})

test_that("meta_r() validates its arguments", {
  r <- c(.28, .35, .22); n <- c(50, 60, 70)
  expect_error(meta_r(c(.5, 1), c(50, 50)), "\\(-1, 1\\)")
  expect_error(meta_r(r, c(50, 60)), "for each study")
  expect_error(meta_r(r, n, reliability_x = c(.8, .9)), "length 1 or one per")
  expect_error(meta_r(c(.9, .95), c(50, 50), reliability_x = 0.5),
               "cannot support")
})
