sem_value <- function(x, term) x$value[match(term, x$term)]

eq_terms <- c("rmsea", "ci_lower", "ci_upper", "rmsea_bound",
              "null_noncentrality", "p_value", "equivalent",
              "chi_square", "df", "N", "n_groups")

test_that("the test of not close fit uses the lower tail at the null noncentrality", {
  x <- equivalence_rmsea(100, 90, 500, 0.05)
  expect_s3_class(x, c("dmar_sem_fit", "dmar_tbl", "data.frame"))
  expect_identical(x$term, eq_terms)
  # The default multiplier is N - 1.
  expect_equal(sem_value(x, "null_noncentrality"), 499 * 90 * 0.05^2)
  expect_equal(sem_value(x, "p_value"), stats::pchisq(100, 90, ncp = 499 * 90 * 0.05^2))
  expect_equal(sem_value(x, "rmsea"), sqrt(10 / (499 * 90)))
  expect_equal(sem_value(x, "equivalent"), 1)
  expect_lt(sem_value(x, "ci_upper"), 0.05)
  y <- equivalence_rmsea(100, 90, 500, 0.05, sample_size_convention = "N")
  expect_equal(sem_value(y, "null_noncentrality"), 500 * 90 * 0.05^2)
  expect_equal(sem_value(y, "p_value"), stats::pchisq(100, 90, ncp = 500 * 90 * 0.05^2))
  z <- equivalence_rmsea(300, 90, 500, 0.05)
  expect_equal(sem_value(z, "equivalent"), 0)
  expect_gt(sem_value(z, "ci_upper"), 0.05)
})

test_that("Yuan et al. (2016, Example 1) is reproduced: T-size RMSEA .083", {
  # T = 28.098 on 23 df with N = 145, the convention N - 1 of the article:
  # delta_t = 22.705 and RMSEA_t = .083.
  x <- equivalence_rmsea(28.098, 23, 145, rmsea_bound = 0.09)
  expect_equal(sem_value(x, "ci_upper"), 0.083, tolerance = 5e-3)
  expect_equal(144 * 23 * sem_value(x, "ci_upper")^2, 22.705, tolerance = 1e-4)
  expect_equal(sem_value(x, "equivalent"), 1)
  expect_equal(sem_value(equivalence_rmsea(28.098, 23, 145, 0.082), "equivalent"), 0)
  expect_equal(sem_value(equivalence_rmsea(28.098, 23, 145, 0.084), "equivalent"), 1)
})

# On CRAN the file stops here, after the anchor above; the blocks below run
# in the full local suite (see cran_files in tests/testthat.R).
skip_on_cran()

test_that("the decision equals the upper limit of the 1 - 2 alpha interval below the bound", {
  skip_on_cran()
  for (alpha in c(0.05, 0.10)) {
    for (chisq_obs in c(5, 20, 35, 50, 80, 110, 160, 230)) {
      for (conv in c("N-1", "N")) {
        z <- equivalence_rmsea(chisq_obs, 50, 400, 0.06, alpha_level = alpha,
                               sample_size_convention = conv)
        expect_equal(sem_value(z, "equivalent"),
                     as.numeric(sem_value(z, "ci_upper") < 0.06))
        expect_equal(attr(z, "conf_level"), 1 - 2 * alpha)
        expect_lte(sem_value(z, "ci_lower"), sem_value(z, "rmsea"))
        expect_true(sem_value(z, "ci_upper") == 0 ||
                      sem_value(z, "rmsea") <= sem_value(z, "ci_upper"))
      }
    }
  }
  # The upper limit is the upper endpoint of ci_rmsea() at 1 - 2 alpha under
  # the same convention.
  x <- equivalence_rmsea(85.30552, 24, 301, 0.08)
  r <- suppressMessages(ci_rmsea(sem_value(x, "rmsea"), 24, 301, conf_level = 0.90))
  expect_equal(sem_value(x, c("ci_lower", "ci_upper")),
               r$value[c(1, 3)], tolerance = 1e-7)
})

test_that("a statistic in the lower alpha tail has an upper limit of 0 and is equivalent", {
  for (chisq_obs in c(0, 10)) {
    x <- equivalence_rmsea(chisq_obs, 24, 301, 0.05)
    expect_equal(sem_value(x, c("ci_lower", "ci_upper")), c(0, 0))
    expect_equal(sem_value(x, "equivalent"), 1)
    expect_match(attr(x, "sem_fit")$boundary, "upper limit is 0")
  }
  # A statistic below df but not in the lower tail keeps a positive limit.
  y <- equivalence_rmsea(20, 24, 301, 0.05)
  expect_equal(sem_value(y, "rmsea"), 0)
  expect_gt(sem_value(y, "ci_upper"), 0)
  expect_lt(sem_value(y, "ci_upper"), 0.05)
  expect_equal(sem_value(y, "equivalent"), 1)
})

test_that("a multiple group statistic carries Steiger's factor", {
  x <- equivalence_rmsea(115.8513, 48, 301, 0.10, n_groups = 2,
                         sample_size_convention = "N")
  expect_equal(sem_value(x, "null_noncentrality"), 301 * 48 * 0.10^2 / 2)
  expect_equal(sem_value(x, "rmsea"), sqrt(2 * (115.8513 - 48) / (301 * 48)))
  # lavaan's values for the two-school model: notclose p at .10 is 0.4290122,
  # and the 90 percent upper limit 0.1195828.
  expect_equal(sem_value(x, "p_value"), 0.4290122, tolerance = 1e-5)
  expect_equal(sem_value(x, "ci_upper"), 0.1195828, tolerance = 1e-5)
  expect_equal(sem_value(x, "equivalent"), 0)
  y <- equivalence_rmsea(115.8513, 48, 301, 0.10, n_groups = 2)
  expect_equal(attr(y, "sem_fit")$multiplier, 299)
})

test_that("the test agrees with lavaan under its normal and Wishart likelihoods", {
  skip_on_cran()
  skip_if_not_installed("lavaan")
  model <- 'visual =~ x1 + x2 + x3
            textual =~ x4 + x5 + x6
            speed =~ x7 + x8 + x9'
  d <- lavaan::HolzingerSwineford1939
  for (lik in c("normal", "wishart")) {
    f <- lavaan::cfa(model, data = d, likelihood = lik)
    fm <- lavaan::fitMeasures(f)
    # Named values straight from fitMeasures(), to check the row names too.
    x <- equivalence_rmsea(fm["chisq"], fm["df"], 301, rmsea_bound = 0.08,
                           sample_size_convention = if (lik == "normal") "N" else "N-1")
    expect_identical(x$term, eq_terms)
    expect_equal(sem_value(x, "p_value"), unname(fm[["rmsea.notclose.pvalue"]]),
                 tolerance = 1e-8)
    expect_equal(sem_value(x, "ci_upper"), unname(fm[["rmsea.ci.upper"]]),
                 tolerance = 1e-6)
    expect_equal(sem_value(x, "rmsea"), unname(fm[["rmsea"]]))
  }
  g <- lavaan::cfa(model, data = d, group = "school")
  fg <- lavaan::fitMeasures(g, fm.args = list(rmsea.h0.notclosefit = 0.10))
  y <- equivalence_rmsea(fg[["chisq"]], fg[["df"]], 301, rmsea_bound = 0.10,
                         n_groups = 2, sample_size_convention = "N")
  expect_equal(sem_value(y, "p_value"), unname(fg[["rmsea.notclose.pvalue"]]),
               tolerance = 1e-8)
  expect_equal(sem_value(y, "ci_upper"), unname(fg[["rmsea.ci.upper"]]),
               tolerance = 1e-6)
})

test_that("the printed table states the verdict, and inputs are validated", {
  skip_on_cran()  # The extended validation grid runs in the full local suite.
  x <- equivalence_rmsea(28.098, 23, 145, rmsea_bound = 0.05)
  out <- utils::capture.output(print(x))
  expect_true(any(grepl("is not rejected", out)))
  expect_true(any(grepl("Confidence level: 90%", out)))
  expect_match(attr(x, "verdict"), "not established")
  expect_error(equivalence_rmsea(100, 0, 500, 0.05), "greater")
  expect_error(equivalence_rmsea(100, 90, 500, 0), "greater")
  expect_error(equivalence_rmsea(100, 90, 500, 0.05, alpha_level = 1), "below")
  expect_error(equivalence_rmsea(100, 90, 500, 0.05, alpha_level = 0.5), "\\(0, 0.5\\)")
  expect_error(equivalence_rmsea(100, 90, 500, 0.05, alpha_level = 0.6), "\\(0, 0.5\\)")
  expect_error(equivalence_rmsea(100, 90, 500.5, 0.05), "integer")
  expect_error(equivalence_rmsea(100, 90, 2, 0.05, n_groups = 2), "exceed")
  expect_error(equivalence_rmsea(100, 90.5, 500, 0.05), "integer")
})

test_that("a statistic in the millions gets the right p-value and verdict", {
  # R's pchisq() returned 0 for P(chi square on 3000 df, noncentrality 3e6,
  # <= 3e6), so a bound of .10 was endorsed; the probability is .193.
  x <- equivalence_rmsea(chisq = 3e6, df = 3000, N = 1e5, rmsea_bound = 0.10,
                         sample_size_convention = "N")
  expect_equal(sem_value(x, "p_value"), 0.1933, tolerance = 1e-3)
  expect_equal(sem_value(x, "equivalent"), 0)
  expect_gt(sem_value(x, "ci_upper"), sem_value(x, "rmsea"))
  skip_on_cran()
  # At 2e6 the p-value .0307 had declared equivalence; it is .0564.
  y <- equivalence_rmsea(chisq = 2e6, df = 3000, N = 1e5, rmsea_bound = 0.08168,
                         sample_size_convention = "N")
  expect_equal(sem_value(y, "p_value"), 0.0564, tolerance = 1e-2)
  expect_equal(sem_value(y, "equivalent"), 0)
})

test_that("the recorded interval names its levels in full whatever getOption(\"digits\") is", {
  # The levels were formatted with getOption("digits") significant
  # digits, so under options(digits = 2) the 97.5% interval of
  # alpha_level = .0125 was named a 98 percent interval.
  withr::local_options(digits = 2)
  z <- equivalence_rmsea(85.30552, 24, 301, 0.08, alpha_level = 0.0125)
  expect_equal(attr(z, "conf_level"), 0.975)
  expect_match(attr(z, "sem_fit")$interval,
               "^two-sided 97\\.5 percent interval; .* one-sided 98\\.75 percent upper bound$")
})
