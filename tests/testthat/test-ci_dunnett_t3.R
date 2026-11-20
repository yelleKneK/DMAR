# Dunnett's T3: separate-variance standard errors, Welch-Satterthwaite
# degrees of freedom, and the studentized maximum modulus (SMM) for the
# a(a - 1)/2 comparisons, at each pair's degrees of freedom.

# The reactivated-story memory scores of Maxwell, Delaney, and Kelley (2027,
# Chapter 12, Exercise 25), percentage correct for three groups of 13 from
# Kroes et al. (2014, Nature Neuroscience, 17, 204--206; doi:10.1038/nn.3609):
# 1 = ECT, tested 24 hours later; 2 = ECT, tested immediately; 3 = no ECT.
# The exercise's answer key reports SPSS's Dunnett's T3 for these scores.
ect <- data.frame(
  y = c(30, 22, 21, 12, 30, 28, 34, 24, 27, 26, 30, 21, 25,
        25, 29, 54, 32, 34, 21, 30, 30, 28, 34, 51, 37, 30,
        60, 36, 32, 60, 57, 73, 74, 44, 25, 33, 36, 40, 67),
  g = factor(rep(1:3, each = 13))
)

# The blood pressure data of Maxwell, Delaney, and Kelley (2027, Chapter 5),
# 24 mild hypertensives in four treatments, as test-ci_games_howell.R uses
# them; the chapter's Table 5.7 prints SPSS's pairwise comparisons for them.
bp <- data.frame(
  y = c(84, 95, 93, 104, 99, 106,     # 1 drug therapy
        81, 84, 92, 101, 80, 108,     # 2 biofeedback
        98, 95, 86, 87, 94, 101,      # 3 diet
        91, 78, 85, 80, 81, 76),      # 4 combination
  g = factor(rep(1:4, each = 6))
)

# Tamhane's T2 for a table's pairs: Welch's t at the Sidak per-comparison
# level, the computation the help page sets beside T3.
t2_pairs <- function(out, n_comparisons, conf_level = .95) {
  cv <- stats::qt((1 + conf_level^(1 / n_comparisons)) / 2, out$df)
  list(p = 1 - (2 * stats::pt(abs(out$t_value), out$df) - 1)^n_comparisons,
       lower = out$mean_difference - cv * out$se,
       upper = out$mean_difference + cv * out$se)
}

test_that("ci_dunnett_t3() reproduces SPSS's Dunnett's T3 for Exercise 25", {
  out <- ci_dunnett_t3(ect$y, ect$g)
  # SPSS prints I - J (1 - 2, 1 - 3, 2 - 3); DMAR reports the later level
  # minus the earlier, so the signs and the roles of the limits flip.
  expect_identical(out$contrast, c("2 - 1", "3 - 1", "3 - 2"))
  expect_equal(round(out$se, 5), c(3.03432, 4.91927, 5.34288))
  expect_equal(round(out$p_adjusted, 3), c(.044, .001, .026))
  expect_equal(round(out$ci_lower, 2), c(0.19, 10.45, 1.60))
  expect_equal(round(out$ci_upper, 2), c(15.96, 36.78, 29.47))
  # The degrees of freedom are those of Welch's t test, unrounded.
  expect_equal(out$df[1],
               as.numeric(stats::t.test(y ~ g, data = droplevels(subset(ect, g != "3")))$parameter),
               tolerance = 1e-10)
})

test_that("the Sidak formula is Tamhane's T2, which T3 never exceeds", {
  # The formula once circulated for these scores as "Tamhane's T3",
  # p = 1 - (2 * pt(|t|, df) - 1)^C with the matching critical value,
  # gives these values; it is Tamhane's T2, and it misses SPSS's T3 in the
  # third decimal through the critical value alone (the standard errors and
  # centers agree). The circulated table printed 29.54 for the last upper
  # limit, 29.53500 at the default tolerance of the uniroot() it used to
  # find the critical value; the critical value itself gives 29.53499.
  out <- ci_dunnett_t3(ect$y, ect$g)
  t2  <- t2_pairs(out, 3)
  expect_equal(round(t2$p, 3), c(.045, .001, .027))
  expect_equal(round(t2$lower, 2), c(0.16, 10.37, 1.54))
  expect_equal(round(t2$upper, 2), c(15.99, 36.86, 29.53))
  expect_true(all(out$p_adjusted < t2$p))
  expect_true(all(out$ci_lower > t2$lower & out$ci_upper < t2$upper))
})

test_that("MDK Table 5.7's SPSS Tamhane rows are T2, and T3 is narrower", {
  # Table 5.7 prints SPSS's Tamhane's T2 for the blood pressure data, in the
  # pair order (1, 2), (1, 3), (1, 4), (2, 3), (2, 4), (3, 4) with I - J
  # differences; ci_dunnett_t3() lists the same pairs as J - I.
  out <- ci_dunnett_t3(bp$y, bp$g)
  t2  <- t2_pairs(out, 6)
  expect_equal(round(out$se, 3), c(5.724, 4.084, 3.957, 5.284, 5.186, 3.288))
  expect_equal(round(t2$p, 3), c(.913, .967, .027, .998, .535, .032))
  expect_equal(round(-t2$upper, 2), c(-13.38, -10.26, 1.65, -21.17, -9.49, 0.91))
  expect_equal(round(-t2$lower, 2), c(25.04, 16.93, 28.35, 16.17, 27.83, 22.42))
  expect_true(all(out$p_adjusted < t2$p))
  expect_true(all(out$ci_upper - out$ci_lower < t2$upper - t2$lower))
})

test_that("with two groups ci_dunnett_t3() is exactly Welch's t test", {
  # The SMM for one comparison is the two-sided t, so the interval and the
  # p-value are those of t.test(var.equal = FALSE).
  d  <- droplevels(subset(bp, g %in% c("1", "4")))
  t3 <- ci_dunnett_t3(d$y, d$g)
  tt <- stats::t.test(y ~ g, data = d, var.equal = FALSE)
  # t.test reports mean_1 - mean_4; ci_dunnett_t3 reports 4 - 1.
  expect_equal(c(t3$ci_lower, t3$ci_upper), -rev(as.numeric(tt$conf.int)),
               tolerance = 1e-12)
  expect_equal(t3$p_adjusted, tt$p.value, tolerance = 1e-12)
  expect_equal(t3$t_value, -as.numeric(tt$statistic), tolerance = 1e-12)
  expect_equal(t3$df, as.numeric(tt$parameter), tolerance = 1e-12)
})

test_that("each half-width is cv_smm()'s critical value times the standard error", {
  skip_on_cran()
  out <- ci_dunnett_t3(bp$y, bp$g, conf_level = .90)
  cv  <- vapply(out$df, function(d)
    cv_smm(alpha_level = .10, df = d, n_comparisons = 6, verbose = FALSE)$value,
    numeric(1))
  expect_equal(out$ci_upper - out$mean_difference, cv * out$se, tolerance = 1e-12)
  expect_equal(out$mean_difference - out$ci_lower, cv * out$se, tolerance = 1e-12)
  expect_equal(attr(out, "conf_level"), .90)
  # A lower family-wise level gives narrower intervals.
  wide <- ci_dunnett_t3(bp$y, bp$g)
  expect_true(all(out$ci_upper - out$ci_lower <
                    wide$ci_upper - wide$ci_lower))
})

test_that("p_adjusted is below alpha when, and only when, the interval excludes zero", {
  skip_on_cran()
  # At its critical value the SMM tail area is alpha to a relative 1e-6 or
  # better, so the test and the interval could disagree only for a |t| that
  # close to the critical value, which none of these is.
  for (d in c(1.3, 4.7, 14, 19.625, 80)) {
    for (m in c(2, 3, 6, 15)) {
      cv <- .smm_upper_cv(.05, d, m)
      expect_equal(.smm_cdf(cv, d, m, lower_tail = FALSE), .05, tolerance = 1e-6)
    }
  }
  for (out in list(ci_dunnett_t3(ect$y, ect$g), ci_dunnett_t3(bp$y, bp$g))) {
    expect_identical(out$p_adjusted < .05,
                     out$ci_lower > 0 | out$ci_upper < 0)
  }
})

test_that("the SMM tail matches an independent integral", {
  # Pinned from a trapezoid rule on 1e6 points of log S from 1e-100 to 1000
  # (the tail_brute() integral of test-cv_smm.R on a wider, finer grid),
  # which is independent of both integrals in R/smm_internals.R; on 2e5
  # points from 1e-12 to 50 no value moves by more than 2e-15. With an
  # identity correlation the multivariate t is the SMM, so mvtnorm::pmvt()
  # computes these tails too, but its randomized lattice rule
  # (GenzBretz(maxpts = 5e6, abseps = 1e-10, releps = 0) after
  # set.seed(113), mvtnorm 1.4.2) misses them by up to 1.3e-7, so it is the
  # live comparison, at 1e-6, in tools/oracle_snippets/pmcmrplus_dunnett_t3.R
  # (run by tools/oracle_checks.R) rather than the source of these values.
  ref <- data.frame(
    q  = c(2.5, 3.2, 4.0733, 2.66, 4.8, 2.9, 3, 6),
    df = c(14, 14, 14, 20, 15, 19, 5, 3),
    m  = c(6, 6, 6, 3, 3, 3, 10, 2),
    p  = c(0.131697456346617, 0.0355629151895955, 0.00654223956427049,
           0.0434797894760439, 0.00069420114873595, 0.026745672793548,
           0.178509488907323, 0.0162722244972054)
  )
  ours <- mapply(.smm_cdf, ref$q, ref$df, ref$m,
                 MoreArgs = list(lower_tail = FALSE))
  expect_lt(max(abs(ours / ref$p - 1)), 1e-9)
})

test_that("the SMM distribution function has its closed forms and its two tails", {
  # One comparison: the two-sided t. Known variance: independent moduli.
  expect_equal(.smm_cdf(2.3, 17, 1, lower_tail = FALSE), 2 * stats::pt(-2.3, 17))
  expect_equal(.smm_cdf(2.3, Inf, 4), (2 * stats::pnorm(2.3) - 1)^4)
  expect_equal(.smm_cdf(2.3, Inf, 4, lower_tail = FALSE),
               1 - (2 * stats::pnorm(2.3) - 1)^4)
  expect_identical(.smm_cdf(0, 12, 3, lower_tail = FALSE), 1)
  expect_identical(.smm_cdf(Inf, 12, 3, lower_tail = FALSE), 0)
  # The tails are integrated in different coordinates; they must sum to one.
  for (q in c(.4, 1.5, 2.7, 4)) {
    expect_equal(.smm_cdf(q, 9.3, 6) + .smm_cdf(q, 9.3, 6, lower_tail = FALSE),
                 1, tolerance = 1e-9)
  }
})

test_that("the SMM tail stays between the single t and Sidak values far out", {
  skip_on_cran()  # The extended grid runs in the full local suite.
  # P(|T| > q) <= P(M > q) <= 1 - (1 - P(|T| > q))^m for every q, df, and m.
  # Far into the tail, where 1 minus the lower tail is zero to machine
  # precision, an integral that missed the mass would return zero or a value
  # orders of magnitude off; the bounds catch both.
  grid <- expand.grid(q = c(.5, 2, 4, 8, 20, 60),
                      df = c(1, 1.3, 2, 4.7, 19.6, 200, 5000),
                      m = c(2, 3, 10, 45))
  p  <- mapply(.smm_cdf, grid$q, grid$df, grid$m,
               MoreArgs = list(lower_tail = FALSE))
  lb <- 2 * stats::pt(-grid$q, grid$df)
  ub <- -expm1(grid$m * log1p(-lb))
  # Where even the single t tail underflows to zero, so may the SMM tail.
  expect_true(all(p[lb > 0] > 0))
  expect_true(all(p >= lb * (1 - 1e-8) & p <= ub * (1 + 1e-8)))
})

test_that("ci_dunnett_t3() agrees with pinned PMCMRplus::dunnettT3Test() values", {
  # Pinned from PMCMRplus::dunnettT3Test() (PMCMRplus 1.9.12, 2026-11-05)
  # after set.seed(113); live comparison in tools/oracle_checks.R. PMCMRplus
  # takes the SMM tail from mvtnorm::pmvt() at its default absolute error of
  # 1e-3 and rounds the Welch degrees of freedom, so its p-values are
  # compared to 1e-3; its statistics are exact. T2 would differ from these by
  # .027 on the first pair of the second design.
  out <- ci_dunnett_t3(ect$y, ect$g)
  expect_equal(out$t_value, c(2.66185566446, 4.80058782773, 2.90825655579),
               tolerance = 1e-10)
  expect_lt(max(abs(out$p_adjusted -
                      c(0.043287292063699, 0.000645511284608, 0.026305023232608))),
            1e-3)
  # Equal sizes and equal variances make every Welch df exactly 14, so the
  # rounding is moot here.
  e <- c(-7, -3, -2, 0, 1, 2, 4, 5)
  d <- data.frame(y = c(10 + e, 14 + e, 22 + e, 25 + e),
                  g = factor(rep(c("a", "b", "c", "d"), each = 8)))
  out <- ci_dunnett_t3(d$y, d$g)
  expect_equal(out$df, rep(14, 6), tolerance = 1e-10)
  expect_equal(out$t_value,
               c(2.03670030886926, 6.11010092660779, 7.63762615825973,
                 4.07340061773852, 5.60092584939047, 1.52752523165195),
               tolerance = 1e-10)
  expect_lt(max(abs(out$p_adjusted -
                      c(2.87497369449e-01, 1.83050360723e-04, 2.42495946335e-05,
                        6.54683327532e-03, 3.39587058900e-04, 5.76636395427e-01))),
            1e-3)
})

test_that("ci_dunnett_t3() returns the documented structure", {
  skip_on_cran()
  out <- ci_dunnett_t3(bp$y, bp$g)
  expect_s3_class(out, "dmar_post_hoc_ci")
  expect_s3_class(out, "dmar_tbl")
  expect_named(out, c("contrast", "mean_difference", "se", "df",
                      "t_value", "ci_lower", "ci_upper",
                      "p_adjusted"))
  expect_equal(nrow(out), 6L)                    # a(a-1)/2 with a = 4
  expect_identical(out$contrast,
                   ci_games_howell(bp$y, bp$g)$contrast)
  expect_equal(out$t_value, out$mean_difference / out$se)
  expect_equal(attr(out, "conf_level"), .95)
  expect_identical(attr(out, "ci_adjustment"), "dunnett_t3")
  expect_identical(attr(out, "p_adjustment"), "dunnett_t3")
  expect_match(paste(utils::capture.output(print(out)), collapse = "\n"),
               "simultaneous interval (Dunnett's T3)", fixed = TRUE)
})

test_that("ci_dunnett_t3() tidies and glances through the shared methods", {
  out <- ci_dunnett_t3(bp$y, bp$g)
  td  <- generics::tidy(out)
  expect_named(td, c("term", "estimate", "ci_lower", "ci_upper",
                     "p_adjusted", "conf_level", "ci_adjustment",
                     "p_adjustment"))
  expect_equal(nrow(td), 6L)
  expect_equal(td$estimate, out$mean_difference)
  expect_equal(generics::glance(out)$n_contrasts, 6L)
  # One row reads as a sentence with its family-wise adjustment.
  expect_match(results_sentence(out[out$contrast == "4 - 1", ]),
               "multiplicity-adjusted p", fixed = TRUE)
})

test_that("ci_dunnett_t3() accepts a fitted model and drops missing values", {
  skip_on_cran()
  a <- ci_dunnett_t3(bp$y, bp$g)
  b <- ci_dunnett_t3(stats::aov(y ~ g, data = bp))
  expect_equal(a$mean_difference, b$mean_difference)
  expect_equal(a$ci_lower, b$ci_lower)
  expect_equal(a$p_adjusted, b$p_adjusted)
  # A missing outcome or group drops that observation, as the complete cases
  # would.
  y <- c(bp$y, NA, 90)
  g <- factor(c(as.character(bp$g), "2", NA), levels = levels(bp$g))
  expect_equal(as.data.frame(ci_dunnett_t3(y, g)), as.data.frame(a))
})

test_that("ci_dunnett_t3() rejects bad input", {
  expect_error(ci_dunnett_t3(bp$y, bp$g, conf_level = 0), "conf_level")
  expect_error(ci_dunnett_t3(bp$y, bp$g, conf_level = 1), "conf_level")
  expect_error(ci_dunnett_t3(bp$y, bp$g, conf_level = NA), "conf_level")
  expect_error(ci_dunnett_t3(bp$y, bp$g[-1]), "same length")
  expect_error(ci_dunnett_t3(bp$y), "numeric vector")
  expect_error(ci_dunnett_t3(c(1, 2, 3), factor(c("a", "a", "a"))), "two groups")
  expect_error(ci_dunnett_t3(c(1, 2), factor(c("a", "b"))), "two observations")
  expect_error(ci_dunnett_t3(stats::lm(y ~ g + I(seq_along(y)), data = bp)),
               "one factor predictor")
})

test_that("ci_dunnett_t3() refuses an infinite outcome value", {
  # An infinite value makes its group's variance NaN, which must not pass for
  # a constant group, whether one group carries one or two do. NaN is a
  # missing value and is dropped, as NA is.
  y <- c(1, 2, 3, Inf, 4, 5, 6, 7, 2, 9, 4, 1)
  g <- factor(rep(c("a", "b", "c"), each = 4))
  expect_error(ci_dunnett_t3(y, g), "infinite values")
  y[5] <- -Inf
  expect_error(ci_dunnett_t3(y, g), "infinite values")
  y <- c(1, 2, 3, 8, 4, 5, 6, 7, 2, 9, 4, 1)
  expect_equal(as.data.frame(ci_dunnett_t3(replace(y, 2, NaN), g)),
               as.data.frame(ci_dunnett_t3(replace(y, 2, NA), g)))
})

test_that("ci_dunnett_t3() errors on a pair of constant groups, naming it", {
  # Both groups essentially constant makes the pairwise standard error zero
  # and the Welch df an undefined 0/0, as in ci_games_howell().
  y <- c(1, 1, 1, 2, 2, 2, 3, 5, 8)
  g <- factor(rep(c("a", "b", "c"), each = 3))
  expect_error(ci_dunnett_t3(y, g), "constant")
  expect_error(ci_dunnett_t3(y, g), "b - a")
})

test_that("ci_dunnett_t3() still works when only one group is constant", {
  # Welch's t handles a single constant group: the standard error is the
  # other group's, and the degrees of freedom are its n - 1.
  y <- c(1, 1, 1, 2, 3, 4)
  g <- factor(rep(c("a", "b"), each = 3))
  out <- ci_dunnett_t3(y, g)
  expect_equal(out$se, sqrt(var(c(2, 3, 4)) / 3), tolerance = 1e-12)
  expect_equal(out$df, 2, tolerance = 1e-12)
})

test_that("ci_dunnett_t3() works at small degrees of freedom and a high level", {
  # Two observations a group give Welch df between 1 and 2, where the
  # chi squared integral alone used to fail at a .99 level; the critical
  # value must still put the tail area at .01.
  y <- c(10, 14, 11, 30, 12, 13.5, 20, 21)
  g <- factor(rep(c("a", "b", "c", "d"), each = 2))
  out <- ci_dunnett_t3(y, g, conf_level = .99)
  expect_true(all(out$df >= 1 & out$df <= 2))
  cv <- (out$ci_upper - out$mean_difference) / out$se
  tails <- mapply(.smm_cdf, cv, out$df, MoreArgs = list(m = 6, lower_tail = FALSE))
  expect_equal(tails, rep(.01, 6), tolerance = 1e-6)
})
