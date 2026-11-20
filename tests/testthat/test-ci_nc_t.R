# The noncentral t probability against the oracle (.oracle_nct() in
# helper-nc_t_oracle.R), inside and beyond the range of pt().

test_that(".nc_t_cdf() is the noncentral t probability, by the oracle, where pt() is not", {
  cdf <- DMAR:::.nc_t_cdf
  # Inside the range of pt(), against pt() and the oracle.
  for (x in list(c(2, 10, 1), c(-3, 10, 1), c(31.62, 38, 35))) {
    for (lower in c(TRUE, FALSE)) {
      # pt() sums its series to an error of about 1e-12.
      expect_lt(abs(cdf(x[1], x[2], x[3], lower) -
                      stats::pt(x[1], x[2], x[3], lower.tail = lower)), 1e-11)
      expect_lt(abs(cdf(x[1], x[2], x[3], lower) -
                      .oracle_nct(x[1], x[2], x[3], lower)), 1e-12)
    }
  }
  # Beyond it: large noncentralities, one and two degrees of freedom with
  # quantiles far in the tail (where an integral over the chi square in
  # one piece returned 0), and a million degrees of freedom.
  for (x in list(c(50, 1, 40), c(100, 99, 95), c(37, 3, 45), c(-40, 5, -45),
                 c(-318309.9, 1, 0), c(1e6, 1, 3e4), c(5000, 2, 4000),
                 c(1e4, 30, 1e4), c(40, 1e6, 41))) {
    for (lower in c(TRUE, FALSE)) {
      expect_lt(abs(cdf(x[1], x[2], x[3], lower) -
                      .oracle_nct(x[1], x[2], x[3], lower)), 1e-12,
                label = paste(c(x, lower), collapse = " "))
    }
  }
  expect_equal(cdf(-318309.9, 1, 0), 1e-6, tolerance = 1e-6)
  expect_equal(cdf(3, Inf, 1), stats::pnorm(3, 1))
  # P(|T| <= k), the probability the planners use.
  abs_prob <- DMAR:::.nc_t_abs_probability
  expect_lt(abs(abs_prob(42, 1598, 40) -
                  (.oracle_nct(42, 1598, 40) - .oracle_nct(-42, 1598, 40))),
            1e-12)
  expect_identical(abs_prob(0, 10, 1), 0)
})

# On CRAN the file stops here, after the anchor above, which comes first for
# that reason; the blocks below run in the full local suite (see cran_files
# in tests/testthat.R).
skip_on_cran()

test_that("ci_nc_t() returns a data frame with expected columns", {
  result <- ci_nc_t(ncp = 2.83, df = 126, conf_level = 0.95)
  expect_s3_class(result, "data.frame")
  expect_true("ci_lower" %in% result$term)
  expect_true("ci_upper" %in% result$term)
  expect_named(result, c("term", "value", "prob_less", "prob_greater"))
})

test_that("ci_nc_t() ci_lower < ncp < ci_upper", {
  ncp <- 2.83
  result <- ci_nc_t(ncp = ncp, df = 126, conf_level = 0.95)
  ll <- result[result$term == "ci_lower", "value"]
  ul <- result[result$term == "ci_upper", "value"]
  expect_lt(ll, ncp)
  expect_gt(ul, ncp)
})

test_that("ci_nc_t() one-sided upper interval: ci_lower is -Inf", {
  result <- ci_nc_t(ncp = 2, df = 30, alpha_lower = 0, alpha_upper = 0.05, conf_level = NULL)
  ll <- result[result$term == "ci_lower", "value"]
  expect_equal(ll, -Inf)
})

test_that("ci_nc_t() one-sided lower interval: ci_upper is Inf", {
  result <- ci_nc_t(ncp = 2, df = 30, alpha_lower = 0.05, alpha_upper = 0, conf_level = NULL)
  ul <- result[result$term == "ci_upper", "value"]
  expect_equal(ul, Inf)
})

test_that("ci_nc_t() errors when df <= 0", {
  expect_error(ci_nc_t(ncp = 1, df = 0), "degrees of freedom")
})

test_that("ci_nc_t() has one name for the observed t value and refuses t_value", {
  expect_false("t_value" %in% names(formals(ci_nc_t)))
  for (call in list(quote(ci_nc_t(t_value = 2.5, df = 50)),
                    quote(ci_nc_t(ncp = 2.5, t_value = 3, df = 50)),
                    quote(ci_nc_t(df = 50)))) {
    err <- tryCatch(eval(call), error = function(e) e)
    expect_s3_class(err, "dmar_argument_invalid")
    expect_identical(err$arg, "ncp")
    expect_match(conditionMessage(err), "'ncp'", fixed = TRUE)
  }
})

test_that("ci_nc_t() achieves target tail probabilities to high precision", {
  # P(T > ncp | nu_lower) should equal alpha_lower; P(T < ncp | nu_upper)
  # should equal alpha_upper. These are the achieved tail probabilities of
  # the search.
  result <- ci_nc_t(ncp = 2.83, df = 126, conf_level = 0.95)
  alpha <- 0.025
  expect_equal(result[result$term == "ci_lower", "prob_greater"], alpha, tolerance = 1e-8)
  expect_equal(result[result$term == "ci_upper", "prob_less"],    alpha, tolerance = 1e-8)
})

test_that("ci_nc_t() central case is symmetric about zero when ncp = 0", {
  # When ncp = 0 the problem is symmetric: nu_L solves P(T > 0 | nu_L) = alpha
  # and nu_U solves P(T < 0 | nu_U) = alpha, so nu_U = -nu_L by symmetry of
  # the noncentral t about its noncentrality.
  result <- ci_nc_t(ncp = 0, df = 30, conf_level = 0.95)
  ll <- result[result$term == "ci_lower", "value"]
  ul <- result[result$term == "ci_upper", "value"]
  expect_equal(ul, -ll, tolerance = 1e-10)
})

test_that("ci_nc_t() asymmetric alphas use them in the correct direction", {
  result <- ci_nc_t(ncp = 2.83, df = 126, alpha_lower = 0.01, alpha_upper = 0.04, conf_level = NULL)
  expect_equal(result[result$term == "ci_lower", "prob_greater"], 0.01, tolerance = 1e-8)
  expect_equal(result[result$term == "ci_upper", "prob_less"],    0.04, tolerance = 1e-8)
})

test_that("ci_nc_t(verbose = FALSE) drops the probability columns", {
  result <- ci_nc_t(ncp = 2.83, df = 126, conf_level = 0.95, verbose = FALSE)
  expect_named(result, c("term", "value"))
})

test_that("ci_nc_t() rejects mixing conf_level with explicit alphas", {
  expect_error(
    ci_nc_t(ncp = 2.83, df = 126, conf_level = 0.95, alpha_lower = 0.025),
    "Specify either"
  )
})

test_that("ci_nc_t() rejects alphas that sum to 1 or more", {
  expect_error(
    ci_nc_t(ncp = 1, df = 30, conf_level = NULL, alpha_lower = 0.6, alpha_upper = 0.5),
    "must be non-negative and sum to less than 1"
  )
})

# ---- The exact noncentral t where pt() is not accurate ---------------------
#
# The oracles (.oracle_nct(), .oracle_tails()) are in helper-nc_t_oracle.R,
# shared with test-nc_t_internals.R.

test_that("ci_nc_t() inverts the exact noncentral t where a limit leaves the range of pt()", {
  # Two groups of 20 at d = 10, 99 percent: the observed t is 31.6 and the
  # upper limit 41.4. pt() is a normal approximation past 37.62, and the
  # interval of DMAR 1.0.x was [7.003, 12.999], 5.996 wide.
  ci <- expect_no_warning(ci_smd(smd = 10, n_1 = 20, n_2 = 20, conf_level = .99))
  expect_equal(ci$value[3] - ci$value[1], 6.095749, tolerance = 1e-6)
  expect_equal(ci$value[1], 7.003356, tolerance = 1e-6)
  res <- ci_nc_t(ncp = 10 * sqrt(10), df = 38, conf_level = .99)
  expect_equal(unname(.oracle_tails(res, 10 * sqrt(10), 38)), c(.005, .005),
               tolerance = 1e-9)
  expect_equal(res$prob_greater[1], .005, tolerance = 1e-9)
  expect_equal(res$prob_less[2], .005, tolerance = 1e-9)

  # Observed t = 34.5 on 38 df: only the upper limit passes 37.62, and
  # DMAR 1.0.x did not warn (the observed value is below 37.62); its upper
  # limit was 42.27323, the lower 26.51062 as now.
  res <- expect_no_warning(ci_nc_t(ncp = 34.5, df = 38))
  expect_equal(res$value, c(26.5106180994, 42.4474546869), tolerance = 1e-10)
  expect_equal(unname(.oracle_tails(res, 34.5, 38)), c(.025, .025),
               tolerance = 1e-9)

  # More than 2,000 degrees of freedom. At 2,500 the upper limit passes
  # 37.62 (1.0.x: 38.19576); above 400,000 pt() is the normal approximation
  # at every noncentrality (1.0.x: 8.039933020, 11.960056980).
  res <- expect_no_warning(ci_nc_t(ncp = 36, df = 2500))
  expect_equal(res$value, c(33.7977695039, 38.1964348015), tolerance = 1e-10)
  expect_equal(unname(.oracle_tails(res, 36, 2500)), c(.025, .025),
               tolerance = 1e-9)
  res <- expect_no_warning(ci_nc_t(ncp = 10, df = 5e5))
  expect_equal(res$value, c(8.039933020, 11.960056981), tolerance = 1e-9)
  expect_equal(unname(.oracle_tails(res, 10, 5e5)), c(.025, .025),
               tolerance = 1e-9)
  # Between 33 and 37.62 with 100,000 degrees of freedom, where pt() sums
  # its series but stops short of convergence in the upper tail.
  res <- ci_nc_t(ncp = 35, df = 1e5)
  expect_equal(unname(.oracle_tails(res, 35, 1e5)), c(.025, .025),
               tolerance = 1e-9)

  # One degree of freedom, where the normal approximation is far off
  # (DMAR 1.0.x: [-107.157, 257.157] at t = 100 and 99 percent).
  res <- expect_no_warning(ci_nc_t(ncp = 100, df = 1, conf_level = .99))
  expect_equal(res$value, c(0.394261657, 280.717411728), tolerance = 1e-8)
  expect_equal(unname(.oracle_tails(res, 100, 1)), c(.005, .005),
               tolerance = 1e-9)
})

test_that("ci_nc_t() keeps the limits of pt() inside its accurate range", {
  # Pinned from ci_nc_t() at 3978364, which inverted pt() everywhere; these
  # limits are inside the range of .nc_t_pt_limit(), where the change
  # leaves them to rounding (2e-14 over 202 intervals in the evidence).
  cases <- list(
    list(10, 38, .95, c(6.99384410888906, 12.9466944627785)),
    list(30, 1000, .95, c(27.6346913940257, 32.3547212736691)),
    list(-20, 126, .99, c(-24.1460795343988, -15.8663907552701)),
    list(0.5, 1, .90, c(-1.31175739400151, 2.12389616076711)),
    list(-2.5, 3, .95, c(-5.14451386386315, 0.315765123561642)),
    list(25, 2500, .99, c(22.2659577319663, 27.7300861291098)))
  for (x in cases) {
    res <- ci_nc_t(ncp = x[[1]], df = x[[2]], conf_level = x[[3]])
    expect_equal(res$value, x[[4]], tolerance = 1e-12,
                 info = paste(x[[1]], x[[2]], x[[3]]))
  }
})

test_that("the range in which pt() is used is the one the evidence supports", {
  # R/nc_t_internals.R records the grid: within 1e-9 of the integral up to
  # 37.62 with at most 2,000 df and up to 33 with more, through 400,000.
  limit <- DMAR:::.nc_t_pt_limit
  expect_identical(limit(1), 37.62)
  expect_identical(limit(2000), 37.62)
  expect_identical(limit(2001), 33)
  expect_identical(limit(4e5), 33)
  expect_identical(limit(4e5 + 1), -Inf)
  expect_true(DMAR:::.nc_t_pt_accurate(38, c(-37.62, 37.62)))
  expect_false(DMAR:::.nc_t_pt_accurate(5000, 34))
})

test_that("ci_nc_t() asks pt() for the tail that does not warn", {
  # pt() warns that full precision may not have been achieved when the
  # probability it sums directly is within 1e-10 of 1; .nc_t_pt() asks for
  # the other tail. DMAR 1.0.x muffled the warning with suppressWarnings().
  nct_pt <- DMAR:::.nc_t_pt
  expect_warning(stats::pt(30, 10, ncp = 0.5), "full precision")
  expect_no_warning(lower <- nct_pt(30, 10, 0.5, TRUE))
  expect_equal(lower, 1 - stats::pt(30, 10, ncp = 0.5, lower.tail = FALSE))
  expect_warning(stats::pt(-30, 10, ncp = 0.5, lower.tail = FALSE), "full precision")
  expect_no_warning(upper <- nct_pt(-30, 10, 0.5, FALSE))
  expect_equal(upper, 1 - stats::pt(-30, 10, ncp = 0.5))
  expect_no_warning(ci_nc_t(ncp = 0.5, df = 10, conf_level = .999))
  expect_no_warning(ci_nc_t(ncp = 37.62, df = 200))
  expect_no_warning(ci_nc_t(ncp = 50, df = 200))
})
