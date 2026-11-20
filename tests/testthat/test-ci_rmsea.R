test_that("ci_rmsea() returns the expected 3-row tidy data frame", {
  rmsea_pt <- sqrt(30 / (15 * 199))
  result <- ci_rmsea(rmsea = rmsea_pt, df = 15, N = 200, conf_level = 0.95)
  expect_s3_class(result, "data.frame")
  expect_named(result, c("term", "value"))
  expect_equal(result$term, c("ci_lower", "rmsea", "ci_upper"))
})

test_that("ci_rmsea() ci_lower <= point estimate <= ci_upper", {
  rmsea_pt <- sqrt(30 / (15 * 199))
  result <- ci_rmsea(rmsea = rmsea_pt, df = 15, N = 200, conf_level = 0.95)
  ll <- result$value[result$term == "ci_lower"]
  est <- result$value[result$term == "rmsea"]
  ul <- result$value[result$term == "ci_upper"]
  expect_gte(ll, 0)
  expect_lte(ll, est)
  expect_lte(est, ul)
})

test_that("ci_rmsea() clamps ci_lower to 0 for a well-fitting (small chi^2) model", {
  rmsea_pt <- sqrt(5 / (15 * 199))
  expect_message(
    result <- ci_rmsea(rmsea = rmsea_pt, df = 15, N = 200, conf_level = 0.95),
    "lower RMSEA limit is set to 0"
  )
  expect_equal(result$value[result$term == "ci_lower"], 0)
})

test_that("the help page example's 90 percent upper limit sits just above 0.05", {
  # Anchors the ?ci_rmsea example prose: for rmsea = .035, df = 40, N = 425
  # the 90 percent upper limit is 0.052, just above the Browne and Cudeck
  # close fit threshold of 0.05 (an earlier version of the page described it
  # as landing below 0.05). Recomputed independently by inverting the
  # noncentral chi square lower tail with uniroot on pchisq directly.
  suppressMessages(
    res <- ci_rmsea(rmsea = .035, df = 40, N = 425, conf_level = .90)
  )
  ul <- res$value[res$term == "ci_upper"]
  chi_sq <- .035^2 * 40 * 424 + 40
  lam <- uniroot(function(l) pchisq(chi_sq, df = 40, ncp = l) - 0.05,
                 c(0, 500), tol = 1e-10)$root
  expect_equal(ul, sqrt(lam / (40 * 424)), tolerance = 1e-6)
  expect_gt(ul, 0.05)
  expect_lt(ul, 0.053)
})

test_that("sample_size_convention sets the multiplier, and the default is N - 1", {
  T <- 85.30552; df <- 24; N <- 301
  # Default: the earlier results, with the multiplier N - 1, here at the
  # default level of 90 percent, so the upper limit leaves .05 below it.
  a <- suppressMessages(ci_rmsea(sqrt((T - df) / ((N - 1) * df)), df, N))
  b <- suppressMessages(ci_rmsea(sqrt((T - df) / ((N - 1) * df)), df, N,
                                 sample_size_convention = "N-1"))
  expect_identical(a$value, b$value)
  lam <- uniroot(function(l) pchisq(T, df = df, ncp = l) - 0.05,
                 c(0, 500), tol = 1e-12)$root
  expect_equal(a$value[3], sqrt(lam / (df * (N - 1))), tolerance = 1e-6)
  # An RMSEA computed with N, as lavaan computes it under its default
  # likelihood, given with the matching convention, reproduces the interval
  # from the chi square statistic.
  r <- suppressMessages(ci_rmsea(sqrt((T - df) / (N * df)), df, N, conf_level = .90,
                                 sample_size_convention = "N"))
  s <- sem_fit_measures(chisq = T, df = df, N = N, conf_level = .90,
                        sample_size_convention = "N")
  expect_equal(r$value[c(1, 3)],
               s$value[match(c("rmsea_ci_lower", "rmsea_ci_upper"), s$term)],
               tolerance = 1e-7)
  expect_error(ci_rmsea(.05, 24, 301, sample_size_convention = "N-2"), "should be one of")
})

test_that("degenerate degrees of freedom, sample sizes, and levels are refused", {
  # A saturated model: lavaan reports rmsea = 0 on 0 df, and the limits were
  # once Inf, the lower one above the estimate.
  expect_error(ci_rmsea(rmsea = 0, df = 0, N = 301), "greater than 0")
  expect_error(ci_rmsea(rmsea = 0.05, df = 0, N = 200), "saturated")
  # N = 1 under "N-1" makes the multiplier 0 (upper limit Inf).
  expect_error(ci_rmsea(0.05, df = 5, N = 1), "greater than 1")
  expect_error(ci_rmsea(0.05, df = 5, N = 0, sample_size_convention = "N"),
               "greater than 0")
  expect_error(ci_rmsea(0.05, df = 5, N = 200, conf_level = 1.5), "'conf_level'")
  expect_error(ci_rmsea(0.05, df = 5, N = 200, conf_level = 0), "'conf_level'")
  # Fractional degrees of freedom, as of a mean and variance adjusted
  # statistic, and N = 1 under "N" are accepted.
  expect_s3_class(suppressMessages(ci_rmsea(0.05, df = 9.92, N = 301)), "dmar_tbl")
  expect_s3_class(suppressMessages(ci_rmsea(0.05, df = 5, N = 1,
                                            sample_size_convention = "N")),
                  "dmar_tbl")
})

test_that("an RMSEA from a statistic in the millions has limits around it", {
  skip_on_cran()
  # T = 3e6 on 3000 df with N = 1e5: the limits once collapsed to 0.0997,
  # below the estimate 0.09995. The noncentrality interval is close to
  # normal at this size.
  x <- ci_rmsea(rmsea = sqrt((3e6 - 3000) / (1e5 * 3000)), df = 3000, N = 1e5,
                conf_level = 0.90, sample_size_convention = "N")
  lambda <- 3e6 - 3000
  h <- stats::qnorm(0.95) * sqrt(2 * (3000 + 2 * lambda))
  expect_equal(x$value[c(1, 3)], sqrt((lambda + c(-h, h)) / (1e5 * 3000)),
               tolerance = 1e-5)
  expect_lt(x$value[1], x$value[2])
  expect_gt(x$value[3], x$value[2])
})

test_that("the default level is 0.90, and the table states it", {
  # Since DMAR 1.1.0 every RMSEA interval defaults to 90 percent, whose two
  # limits decide the tests of close and not close fit at .05; DMAR 1.0.x
  # defaulted to 95 percent.
  expect_equal(eval(formals(ci_rmsea)$conf_level), 0.90)
  x <- ci_rmsea(rmsea = .035, df = 40, N = 425)
  expect_identical(x$value, ci_rmsea(rmsea = .035, df = 40, N = 425,
                                     conf_level = 0.90)$value)
  expect_equal(attr(x, "conf_level"), 0.90)
  expect_true("Confidence level: 90%" %in% utils::capture.output(print(x)))
  y <- ci_rmsea(rmsea = .035, df = 40, N = 425, conf_level = 0.95)
  expect_true("Confidence level: 95%" %in% utils::capture.output(print(y)))
  expect_lt(y$value[1], x$value[1])
  expect_gt(y$value[3], x$value[3])
  expect_identical(dmar_metadata(x)$rmsea_conf_level, 0.90)
})

test_that("an interval at the default level, or at a level given, gives no warning", {
  expect_no_warning(x <- ci_rmsea(.035, 40, 425))
  expect_equal(attr(x, "conf_level"), 0.90)
  expect_no_warning(ci_rmsea(.035, 40, 425, conf_level = 0.90))
  expect_no_warning(ci_rmsea(.035, 40, 425, conf_level = 0.95))
  expect_no_warning(suppressMessages(
    ci_rmsea(.035, 40, 425, alpha_lower = 0.05, alpha_upper = 0.05)))
  expect_error(ci_rmsea(.035, 0, 425), "saturated")
})
