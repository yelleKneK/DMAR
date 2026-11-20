# Live oracle comparisons for ci_dunnett_t3() and the studentized maximum
# modulus (SMM) behind it and cv_smm(), run at release time after DMAR is
# installed, when PMCMRplus (and mvtnorm, which it imports) are installed
# locally. Neither is a DMAR dependency; the shipped tests pin the values.
## Written against PMCMRplus 1.9.12 and mvtnorm 1.4.2.

## from tests/testthat/test-ci_dunnett_t3.R
## ci_dunnett_t3() against PMCMRplus::dunnettT3Test().
local({
  # PMCMRplus takes the SMM tail from mvtnorm::pmvt() at its default absolute
  # error of 1e-3 and rounds the Welch degrees of freedom to integers, so the
  # p-values agree to 1e-3 and the t statistics exactly. Both list the pairs
  # as later level minus earlier, in the same order.
  e <- c(-7, -3, -2, 0, 1, 2, 4, 5)
  equal_df <- data.frame(y = c(10 + e, 14 + e, 22 + e, 25 + e),
                         g = factor(rep(c("a", "b", "c", "d"), each = 8)))
  # Maxwell, Delaney, and Kelley (2027, Chapter 12, Exercise 25), reactivated.
  ect <- data.frame(
    y = c(30, 22, 21, 12, 30, 28, 34, 24, 27, 26, 30, 21, 25,
          25, 29, 54, 32, 34, 21, 30, 30, 28, 34, 51, 37, 30,
          60, 36, 32, 60, 57, 73, 74, 44, 25, 33, 36, 40, 67),
    g = factor(rep(1:3, each = 13)))
  for (d in list(equal_df, ect)) {
    ours <- DMAR::ci_dunnett_t3(d$y, d$g)
    set.seed(113)  # pmvt() is a randomized lattice rule; fix its RNG state
    theirs <- PMCMRplus::dunnettT3Test(d$y, d$g)
    stat <- as.vector(theirs$statistic); stat <- stat[!is.na(stat)]
    p    <- as.vector(theirs$p.value);   p    <- p[!is.na(p)]
    stopifnot(isTRUE(all.equal(ours$t_value, stat, tolerance = 1e-10)),
              max(abs(ours$p_adjusted - p)) < 1e-3)
  }
})

## from tests/testthat/test-ci_dunnett_t3.R
## The SMM tail and critical value against high-precision mvtnorm.
local({
  # With an identity correlation the multivariate t is the SMM. Equal sizes
  # and variances make every Welch df exactly 14, an integer as pmvt() needs.
  e <- c(-7, -3, -2, 0, 1, 2, 4, 5)
  d <- data.frame(y = c(10 + e, 14 + e, 22 + e, 25 + e),
                  g = factor(rep(c("a", "b", "c", "d"), each = 8)))
  ours <- DMAR::ci_dunnett_t3(d$y, d$g)
  # At this precision pmvt() takes about half a minute for the seven calls;
  # a looser rule cannot resolve the 1e-6 agreement.
  precise <- mvtnorm::GenzBretz(maxpts = 5e6, abseps = 1e-10, releps = 0)
  set.seed(113)
  theirs <- vapply(abs(ours$t_value), function(q)
    1 - mvtnorm::pmvt(lower = rep(-q, 6), upper = rep(q, 6), df = 14,
                      corr = diag(6), algorithm = precise)[1], numeric(1))
  stopifnot(max(abs(ours$p_adjusted - theirs)) < 1e-6)
  # The critical value leaves .05 in the tail (checked at the value rather
  # than by mvtnorm::qmvt(), whose search at this precision takes minutes).
  cv <- DMAR::cv_smm(.05, df = 14, n_comparisons = 6, verbose = FALSE)$value
  set.seed(113)
  tail_at_cv <- 1 - mvtnorm::pmvt(lower = rep(-cv, 6), upper = rep(cv, 6),
                                  df = 14, corr = diag(6),
                                  algorithm = precise)[1]
  stopifnot(abs(tail_at_cv - .05) < 1e-6)
})

## from tests/testthat/test-ci_dunnett_t3.R
## The eight SMM tail cells the test pins, against high-precision mvtnorm.
local({
  # The test pins these tails from a deterministic trapezoid rule on log S.
  # pmvt() at this precision, after set.seed(113) and in this order, misses
  # them by up to 1.3e-7, so DMAR's integral and the pinned values are both
  # compared with it at 1e-6. The eight calls take about 20 s.
  ref <- data.frame(
    q  = c(2.5, 3.2, 4.0733, 2.66, 4.8, 2.9, 3, 6),
    df = c(14, 14, 14, 20, 15, 19, 5, 3),
    m  = c(6, 6, 6, 3, 3, 3, 10, 2),
    p  = c(0.131697456346617, 0.0355629151895955, 0.00654223956427049,
           0.0434797894760439, 0.00069420114873595, 0.026745672793548,
           0.178509488907323, 0.0162722244972054))
  precise <- mvtnorm::GenzBretz(maxpts = 5e6, abseps = 1e-10, releps = 0)
  set.seed(113)
  theirs <- mapply(function(q, df, m)
    1 - mvtnorm::pmvt(lower = rep(-q, m), upper = rep(q, m), df = df,
                      corr = diag(m), algorithm = precise)[1],
    ref$q, ref$df, ref$m)
  ours <- mapply(DMAR:::.smm_cdf, ref$q, ref$df, ref$m,
                 MoreArgs = list(lower_tail = FALSE))
  stopifnot(max(abs(theirs - ref$p)) < 1e-6,
            max(abs(theirs - ours)) < 1e-6)
})

## from tests/testthat/test-ci_dunnett_t3.R
## Tamhane's T2, which the help page sets beside T3, is the Sidak formula.
local({
  # PMCMRplus::tamhaneT2Test() adjusts each pair's two-sided Welch p-value by
  # 1 - (1 - p)^C; the help page describes T2 that way, and T3's adjusted
  # p-values never exceed it.
  ect <- data.frame(
    y = c(30, 22, 21, 12, 30, 28, 34, 24, 27, 26, 30, 21, 25,
          25, 29, 54, 32, 34, 21, 30, 30, 28, 34, 51, 37, 30,
          60, 36, 32, 60, 57, 73, 74, 44, 25, 33, 36, 40, 67),
    g = factor(rep(1:3, each = 13)))
  ours <- DMAR::ci_dunnett_t3(ect$y, ect$g)
  t2 <- as.vector(PMCMRplus::tamhaneT2Test(ect$y, ect$g)$p.value)
  t2 <- t2[!is.na(t2)]
  sidak <- 1 - (2 * stats::pt(abs(ours$t_value), ours$df) - 1)^3
  stopifnot(isTRUE(all.equal(t2, sidak, tolerance = 1e-10)),
            all(ours$p_adjusted < t2))
})
