# Worked examples from printed sources, pinned at full precision.
#
# Each test opens with its record: a "Record:" line naming the source, the
# printed page, and the table or equation, and a "Category:" line. Where
# one test checks values of different categories, each value is labeled.
# The categories are kept apart:
#   exact reproduction     DMAR (or, where stated, the arithmetic of the
#                          test) reproduces a value the source states
#                          exactly.
#   rounded comparison     the source prints a rounded value; the full
#                          precision value (DMAR's, or, where no DMAR
#                          function computes it, the test's own value of an
#                          input) is pinned to an independent calculation
#                          and compared with the print only at the printed
#                          precision.
#   constructed illustration
#                          inputs built here to realize what the source
#                          describes (not the source's raw data).
#   adaptation             the source's data or design used for a DMAR
#                          analysis or method that the source does not use;
#                          DMAR's value is pinned to direct algebra or to
#                          a re-implementation, never to a printed value.
#   algorithm oracle       a formula of the source, coded separately, checks
#                          a DMAR calculation or internal.
# A printed value that comes from interpolation in a table or from a chart
# read by eye is not a tight oracle and is not compared.
# Full precision values were computed independently (direct algebra, base R
# distribution functions, and boot) before these tests were written.
#
# Sources:
#   Anderson, T. W., An Introduction to Multivariate Statistical Analysis,
#     Chapter 4 (Sections 4.3 and 4.4), pp. 140, 143, 144, 152.
#   Cohen, J., Statistical Power Analysis for the Behavioral Sciences,
#     Chapters 8 and 9, pp. 277, 410, 411, 414, 434, 439.
#   Grissom, R. J., & Kim, J. J., Effect Sizes for Research: Univariate and
#     Multivariate Applications, Chapter 5, pp. 150, 151, 160.
#   Efron, B., & Tibshirani, R. J., An Introduction to the Bootstrap,
#     Chapter 14, pp. 179, 180, 185, 186.
# The pages were read from photographs; no title or copyright page was
# photographed, so the editions are not confirmed here.

# ---- Anderson: partial and multiple correlation, Hooker's hay data -----

# Correlations of yield (X1), spring rainfall (X2), and accumulated spring
# temperature (X3) over N = 20 years (p. 140, display (16)).
anderson_R <- matrix(c(1, .80, -.40,
                       .80, 1, -.56,
                       -.40, -.56, 1), 3, byrow = TRUE)
anderson_r13_2 <- (-.40 - .80 * -.56) / sqrt((1 - .80^2) * (1 - .56^2))
anderson_r12_3 <- (.80 - -.40 * -.56) / sqrt((1 - .40^2) * (1 - .56^2))

test_that("Anderson p. 140, Eq. (17), and p. 144: the partial correlations", {
  # Record: Anderson, p. 140, display (16) and Eq. (17); p. 144, Example 1.
  # Category: rounded comparison. Anderson prints r_13.2 = 0.097 (p. 140)
  # and r_12.3 = 0.759 (p. 144). The first-order formula and the precision
  # matrix are two independent routes to the full precision values.
  P <- solve(anderson_R)
  expect_equal(anderson_r13_2, 0.0965609099170536, tolerance = 1e-13)
  expect_equal(anderson_r12_3, 0.758567449489217, tolerance = 1e-13)
  expect_equal(-P[1, 3] / sqrt(P[1, 1] * P[3, 3]), anderson_r13_2,
               tolerance = 1e-13)
  expect_equal(-P[1, 2] / sqrt(P[1, 1] * P[2, 2]), anderson_r12_3,
               tolerance = 1e-13)
  expect_identical(round(anderson_r13_2, 3), 0.097)
  expect_identical(round(anderson_r12_3, 3), 0.759)
})

test_that("Anderson pp. 143-144, Theorem 4.3.5: Fisher Z interval for rho_12.3 at effective N 19", {
  # Record: Anderson, p. 143, Theorem 4.3.5; p. 144, Example 1 and Eq. (27).
  # Category: algorithm oracle. Theorem 4.3.5: a partial correlation from N
  # cases with p - q variables partialed out has the distribution of a
  # simple correlation from N - (p - q) cases. With one variable partialed
  # out of 20 cases the effective size is 19, so ci_r() is called with
  # N = 19, and Fisher's Z with the scale sqrt(N - (p - q) - 3) of
  # Eq. (27), coded separately, checks it. That 19 is the effective sample
  # size of the theorem, not the study's sample size, which is 20.
  #
  # Not compared: Anderson's Example 1 reads David's charts at 19 and
  # prints .50 < rho_12.3 < .88. David's charts are read by eye from the
  # exact distribution of r; ci_r() is the Fisher Z interval. The two are
  # different methods, and the chart reading is not a tight oracle for
  # either: the exact-distribution interval, from the density of r
  # integrated numerically and inverted in rho, is [.4533194, .8957938],
  # and the Fisher Z interval is [.4643488, .9019956].
  res <- ci_r(r = anderson_r12_3, N = 19, conf_level = 0.95)
  z <- atanh(anderson_r12_3)
  by_hand <- tanh(z + c(-1, 1) * qnorm(0.975) / sqrt(19 - 3))
  expect_equal(res$value[res$term == "ci_lower"], 0.464348784510494,
               tolerance = 1e-12)
  expect_equal(res$value[res$term == "ci_upper"], 0.901995609828508,
               tolerance = 1e-12)
  expect_equal(res$value[res$term == "ci_lower"], by_hand[1], tolerance = 1e-13)
  expect_equal(res$value[res$term == "ci_upper"], by_hand[2], tolerance = 1e-13)
  expect_equal(res$value[res$term == "r"], anderson_r12_3, tolerance = 1e-14)
})

test_that("Anderson p. 144, Theorem 4.3.5 and Eq. (27): var_partial_r() at the study's N", {
  # Record: Anderson, p. 143, Theorem 4.3.5; p. 144, Eq. (27) and
  # Example 2, which prints the factor sqrt(20 - 1 - 3).
  # Category: algorithm oracle. The scale of Eq. (27), sqrt(N - (p - q) -
  # 3), coded separately as the variance of Fisher's Z, 1 / (20 - 1 - 3) =
  # 1/16, checks var_partial_r() called with the study's N = 20. The page
  # prints the factor, not the variance, so this is not called an exact
  # reproduction. var_partial_r() takes p = 2 (X2, the focal predictor, and
  # X3, the one variable partialed out), and its Fisher Z variance is
  # 1 / (N - p - 2) = 1 / (20 - 2 - 2) = 1/16.
  v <- var_partial_r(r = anderson_r12_3, N = 20, p = 2, fisher_z = TRUE)
  expect_equal(v$value, 1 / (20 - 1 - 3), tolerance = 1e-15)
  expect_equal(v$value, 1 / 16, tolerance = 1e-15)
})

test_that("Anderson p. 144, Example 2, Eq. (27): Fisher's Z test of rho_13.2 = 0", {
  # Record: Anderson, p. 144, Example 2 and Eq. (27).
  # Category: rounded comparison. Anderson computes sqrt(20 - 1 - 3)
  # (0.0973) = 0.3892 from the printed, rounded r_13.2 = .097. With that
  # rounded input convert_r_Z() reproduces the printed statistic at the
  # printed precision; the full precision partial correlation gives
  # 0.3874508.
  z_rounded <- convert_r_Z(0.097)$value
  expect_equal(z_rounded, atanh(0.097), tolerance = 1e-15)
  expect_identical(round(z_rounded, 4), 0.0973)
  expect_identical(round(sqrt(20 - 1 - 3) * z_rounded, 4), 0.3892)
  z_full <- convert_r_Z(anderson_r13_2)$value
  expect_equal(sqrt(16) * z_full, 0.387450846939011, tolerance = 1e-13)
  expect_lt(abs(sqrt(16) * z_full), qnorm(0.975))
})

test_that("Anderson p. 152, Eq. (31) and Theorem 4.4.3: multiple R and the F test", {
  # Record: Anderson, p. 152, Eqs. (31) and (32), Theorem 4.4.3.
  # Category: rounded comparison, for every value: 1 - R^2 (printed 0.357),
  # R (0.802), the F statistic (15.3), and F_{2,17}(.01) (6.11). Eq. (31)
  # prints the two determinants as matrices; their values, .2448 and .6864,
  # are computed here and are not printed. At full precision
  # R^2 = .4416 / .6864 and F = [R^2 / (1 - R^2)] (17 / 2) is exactly 46/3.
  # The error degrees of freedom: Anderson's (N - p) with his p = 3 (the
  # outcome and the two predictors), printed as (20 - 3) / (3 - 1), and
  # DMAR's N - p - 1 with p = 2 predictors both give 17.
  R2 <- 1 - det(anderson_R) / det(anderson_R[2:3, 2:3])
  expect_equal(R2, 0.4416 / 0.6864, tolerance = 1e-13)
  expect_equal(sqrt(R2), 0.802095158542079, tolerance = 1e-13)
  expect_identical(round(1 - R2, 3), 0.357)
  expect_identical(round(sqrt(R2), 3), 0.802)
  F_res <- convert_R2_f(R2 = R2, N = 20, p = 2)
  expect_equal(F_res$value, 46 / 3, tolerance = 1e-12)
  expect_identical(round(F_res$value, 1), 15.3)
  expect_identical(round(qf(0.99, 2, 17), 2), 6.11)
})

# ---- Cohen (1988): f, f^2, and the noncentrality of Chapter 9 ----------
#
# DMAR keeps its own power method: with fixed predictors the noncentrality
# of the test of q predictors added to w others is lambda = N f^2.
# Cohen computes lambda = f^2 (u + v + 1) (Eq. 9.3.1, p. 414), which for an
# increment is f^2 (N - w), and reads power from Tables 9.3.1 and 9.3.2 by
# interpolation. DMAR offers no Cohen-style option (the maintainer's default
# for 1.1.0). Cohen's numbers can be reproduced numerically by an omnibus
# call (p_restricted = 0) at the reduced size u + v + 1; that reduced size
# is a device for reproducing his noncentrality and is not the study's
# sample size.
#
# Two helpers keep the methods apart. dmar_power_01() is DMAR's own power
# at the study's sample size. cohen_convention_power() is the
# numerical-equivalent omnibus call at the reduced size u + v + 1. ncf_power()
# is the power of the noncentral F test at a given noncentrality, coded with
# base R's pf().

dmar_power_01 <- function(...) {
  # DMAR's fixed-predictor power (lambda = N f^2) at alpha .01; the N given
  # is the study's sample size.
  res <- ss_power_R2_change(..., alpha_level = 0.01, random_predictors = FALSE)
  res$value[res$term == "actual_power"]
}
cohen_convention_power <- function(f2, u, v) {
  # Cohen's lambda = f^2 (u + v + 1) through DMAR: an omnibus call with u
  # predictors at specified_N = u + v + 1, so that N f^2 is Cohen's lambda
  # and the error degrees of freedom are v. That specified_N is not the
  # study's sample size.
  res <- ss_power_R2_change(true_f2 = f2, p_full = u, p_restricted = 0,
                            specified_N = u + v + 1, alpha_level = 0.01,
                            random_predictors = FALSE)
  res$value[res$term == "actual_power"]
}
ncf_power <- function(lambda, u, v, a = 0.01) {
  1 - pf(qf(1 - a, u, v), u, v, ncp = lambda)
}

test_that("ss_power_R2_change() has no Cohen-style noncentrality option", {
  # Record: Cohen, p. 414, Eq. (9.3.1) (method).
  # Category: none; a statement of DMAR's method, not a numerical
  # comparison. The maintainer's default for 1.1.0: no option that switches
  # to f^2 (u + v + 1). The distinction is documented instead.
  expect_false(any(grepl("cohen|convention|noncentral",
                         names(formals(ss_power_R2_change)),
                         ignore.case = TRUE)))
})

test_that("Cohen p. 277, Eq. (8.2.8): seven means in Pattern 1 give f = .5 / sqrt(14)", {
  # Record: Cohen, p. 277, Eq. (8.2.8), Pattern 1.
  # Category: rounded comparison. Means -d/2, 0, 0, 0, 0, 0, +d/2 with
  # range d = .5 within-population standard deviations: f = d sqrt(1 / 14),
  # printed .267d and .13. The arguments are given by position so that the
  # call does not depend on the names of cohen_f()'s first two arguments.
  res <- cohen_f(c(-0.25, 0, 0, 0, 0, 0, 0.25), 1)
  expect_equal(res$value, 0.5 / sqrt(14), tolerance = 1e-14)
  expect_equal(res$value, 0.133630620956212, tolerance = 1e-13)
  expect_identical(round(sqrt(1 / 14), 3), 0.267)
  expect_identical(round(res$value, 2), 0.13)
})

test_that("Cohen pp. 410-411, Eqs. (9.2.2)-(9.2.4): f^2 in Cases 0, 1, and 2", {
  # Record: Cohen, p. 410, Eqs. (9.2.2) and (9.2.3); p. 411, Eq. (9.2.4).
  # Categories, value by value:
  #   Case 0, exact reproduction: .20 / (1 - .20) = .25, printed .25
  #   (p. 410), from ss_power_R2_change()'s true_f2.
  #   Case 1, rounded comparison: (.45 - .30) / (1 - .45) = 3/11, printed
  #   .15/.55 = .2727 (p. 410), from ss_power_R2_change()'s true_f2.
  #   Case 2, exact reproduction of the printed arithmetic only: .15 / .40
  #   = .3750 (p. 411), with the error from the larger model,
  #   R^2_Y.A,B,C = .60. No DMAR call is checked: DMAR takes Case 2
  #   through true_f2, with p_full counting the larger model (the help's
  #   "Error from a larger model").
  f2 <- function(...) {
    res <- ss_power_R2_change(..., specified_N = 100, random_predictors = FALSE)
    res$value[res$term == "true_f2"]
  }
  expect_equal(f2(true_R2_full = 0.20, true_R2_restricted = 0, p_full = 1,
                  p_restricted = 0), 0.25, tolerance = 1e-15)
  expect_equal(f2(true_R2_full = 0.45, true_R2_restricted = 0.30, p_full = 2,
                  p_restricted = 1), 3 / 11, tolerance = 1e-15)
  expect_identical(round(f2(true_R2_full = 0.45, true_R2_restricted = 0.30,
                            p_full = 2, p_restricted = 1), 4), 0.2727)
  expect_equal((0.45 - 0.30) / (1 - 0.60), 0.375, tolerance = 1e-15)
})

test_that("Cohen p. 434, Example 9.10: DMAR's power at N = 90 and Cohen's N - w convention", {
  # Record: Cohen, p. 434, Example 9.10; p. 414, Eq. (9.3.1); Table 9.3.1.
  # Design: u = 4 ethnic-group predictors added to w = 5 for age and
  # education, R^2 from .25 to .30, N = 90, alpha .01, v = 90 - 4 - 5 - 1 =
  # 80, f^2 = .05 / .70 = .0714.
  # Categories, value by value:
  #   DMAR's power at the study's N = 90, from lambda = N f^2 = 6.43, is
  #   .2413: adaptation (Cohen's design, a method Cohen does not use),
  #   pinned to pf() at that noncentrality.
  #   Cohen's convention at N = 90, from lambda = .0714 (85) = 6.1, is
  #   .2225: algorithm oracle. Eq. (9.3.1) coded with pf() checks the
  #   numerical-equivalent omnibus call at 85 = u + v + 1.
  #   Not compared: Cohen's printed .24, an interpolation in Table 9.3.1
  #   (.22 at v = 60, .24 at v = 120), is not a tight oracle for either.
  ours <- ss_power_R2_change(true_R2_full = 0.30, true_R2_restricted = 0.25,
                             p_full = 9, p_restricted = 5, specified_N = 90,
                             alpha_level = 0.01, random_predictors = FALSE)
  val <- function(t) ours$value[ours$term == t]
  expect_equal(val("actual_power"), 0.241296492072416, tolerance = 1e-12)
  expect_equal(val("actual_power"), ncf_power(90 * 0.05 / 0.70, 4, 80),
               tolerance = 1e-13)
  expect_equal(val("noncentral_f_parm"), 90 / 14, tolerance = 1e-14)
  expect_equal(val("df_2"), 80)
  expect_equal(ncf_power(85 * 0.05 / 0.70, 4, 80), 0.222508250709684,
               tolerance = 1e-12)
  # 85 = u + v + 1 is not the study's sample size of 90.
  expect_equal(cohen_convention_power(0.05 / 0.70, u = 4, v = 80),
               0.222508250709684, tolerance = 1e-12)
  # The same example with N doubled to 180 (v = 170).
  #   DMAR's power at the study's N = 180 (lambda 12.86), .6168: adaptation.
  #   Cohen's convention (lambda = .0714 (175) = 12.5), .5991: algorithm
  #   oracle, through the omnibus call at 175 = u + v + 1, which is not
  #   the study's sample size of 180.
  #   Not compared: the printed .61, again an interpolation.
  expect_equal(dmar_power_01(true_R2_full = 0.30, true_R2_restricted = 0.25,
                             p_full = 9, p_restricted = 5, specified_N = 180),
               0.616809467370202, tolerance = 1e-12)
  expect_equal(dmar_power_01(true_R2_full = 0.30, true_R2_restricted = 0.25,
                             p_full = 9, p_restricted = 5, specified_N = 180),
               ncf_power(180 * 0.05 / 0.70, 4, 170), tolerance = 1e-13)
  expect_equal(ncf_power(175 * 0.05 / 0.70, 4, 170), 0.599078496755353,
               tolerance = 1e-12)
  # 175 = u + v + 1 is not the study's sample size of 180.
  expect_equal(cohen_convention_power(0.05 / 0.70, u = 4, v = 170),
               0.599078496755353, tolerance = 1e-12)
})

test_that("Cohen p. 439, Example 9.13: Case 2 tests with error from all three groups", {
  # Record: Cohen, p. 439, Example 9.13 (the tests of H.D and of D);
  # Eqs. (9.2.4), (9.3.1), and (9.3.10); Table 9.3.1.
  # Design: N = 200; group D has 5 variables, H has 3, and S has 4, and
  # R^2_Y.D,H,S = .20, so the Model 2 error is .80 on v = 187 df.
  #
  # H.D, Case 2 with set A = D (w = 5), set B = H (u = 3), and set C = S
  # (z = 4): f^2 = .06 / .80 = .075, whose numerator is R^2_Y.D,H -
  # R^2_Y.D; S enters only the error.
  #   Cohen's convention, lambda = .075 (191) = 14.325, gives .7407:
  #   algorithm oracle for the omnibus call at 191 = u + v + 1, which is
  #   not the study's sample size of 200. The printed .74 is read from
  #   Table 9.3.1 by interpolation at the printed lambda 14.3, so it is not
  #   a tight oracle and is not compared.
  #   DMAR's power, from lambda = N f^2 = 15, is .7670: adaptation (the
  #   value test-ss_power_R2_change.R pins to the G*Power 3.1 manual,
  #   Section 14.3.2; that manual was not read for this file). DMAR's call
  #   gives true_f2 = .075, whose denominator is the Model 2 error, with
  #   p_full = 12 (the model that supplies the error) and p_restricted =
  #   12 - 3 (the help's "Error from a larger model").
  expect_equal(ncf_power(0.075 * 191, 3, 187), 0.740721888312373,
               tolerance = 1e-12)
  # 191 = u + v + 1 is not the study's sample size of 200.
  expect_equal(cohen_convention_power(0.075, u = 3, v = 187),
               0.740721888312373, tolerance = 1e-12)
  expect_equal(dmar_power_01(true_f2 = 0.075, p_full = 12,
                             p_restricted = 12 - 3, specified_N = 200),
               0.766989564222325, tolerance = 1e-12)
  expect_equal(dmar_power_01(true_f2 = 0.075, p_full = 12,
                             p_restricted = 12 - 3, specified_N = 200),
               ncf_power(200 * 0.075, 3, 187), tolerance = 1e-13)
  #
  # D, Case 2 with set A empty (w = 0), set B = D (u = 5), and set C = H
  # and S (z = 7): f^2 = .10 / .80 = .125; H and S enter only the error.
  #   Cohen's convention, lambda = .125 (193) = 24.125 (printed 24.1),
  #   gives .9186: algorithm oracle for the omnibus call at 193 = u + v + 1,
  #   which is not the study's sample size of 200. The printed power "92"
  #   (that is, .92) is read from Cohen's power table by interpolation, so it is not
  #   a tight oracle and is not compared.
  #   DMAR's power, from lambda = N f^2 = 25, is .9303: adaptation. DMAR's
  #   call gives true_f2 = .125 with p_full = 12 (the model that supplies
  #   the error) and p_restricted = 12 - 5.
  expect_equal(ncf_power(0.125 * 193, 5, 187), 0.918581137810925,
               tolerance = 1e-12)
  # 193 = u + v + 1 is not the study's sample size of 200.
  expect_equal(cohen_convention_power(0.125, u = 5, v = 187),
               0.918581137810925, tolerance = 1e-12)
  expect_equal(dmar_power_01(true_f2 = 0.125, p_full = 12,
                             p_restricted = 12 - 5, specified_N = 200),
               0.930262941528867, tolerance = 1e-12)
  expect_equal(dmar_power_01(true_f2 = 0.125, p_full = 12,
                             p_restricted = 12 - 5, specified_N = 200),
               ncf_power(200 * 0.125, 5, 187), tolerance = 1e-13)
})

# ---- Grissom and Kim: probability of superiority and the CL ------------

test_that("Grissom and Kim pp. 150-151, Eqs. (5.1)-(5.2): A from wins and ties", {
  # Record: Grissom and Kim, p. 150, Eqs. (5.1)-(5.2); pp. 150-151, ties.
  # Category: constructed illustration. The book describes samples of
  # three, in which sample a wins five of the nine pairings (estimate 5/9 =
  # .56), and then four wins and one tie, with half of the tie to each
  # sample (4.5/9 = .5). The scores below are built to realize those
  # counts; they are not the book's data, which it does not give.
  five_wins <- vargha_delaney_A(c(1, 4, 6), c(2, 3, 5), test = "t")
  expect_equal(sum(outer(c(1, 4, 6), c(2, 3, 5), ">")), 5L)
  expect_equal(five_wins$A, 5 / 9, tolerance = 1e-14)
  expect_identical(round(five_wins$A, 2), 0.56)
  one_tie <- vargha_delaney_A(c(1, 3, 6), c(2, 3, 5), test = "t")
  expect_equal(sum(outer(c(1, 3, 6), c(2, 3, 5), ">")), 4L)
  expect_equal(sum(outer(c(1, 3, 6), c(2, 3, 5), "==")), 1L)
  expect_equal(one_tie$A, 4.5 / 9, tolerance = 1e-14)
  # Cliff's delta on the same scores is 2A - 1, with the tie in p_tied.
  cd <- cliff_delta(c(1, 3, 6), c(2, 3, 5))
  expect_equal(cd$value[cd$term == "cliff_delta"], 0, tolerance = 1e-14)
  expect_equal(cd$value[cd$term == "p_tied"], 1 / 9, tolerance = 1e-14)
  expect_equal(cliff_delta(c(1, 4, 6), c(2, 3, 5))$value[1], 2 * 5 / 9 - 1,
               tolerance = 1e-14)
})

test_that("Grissom and Kim p. 160, Eqs. (5.6)-(5.7): the CL from a standardized mean difference", {
  # Record: Grissom and Kim, p. 160, Eqs. (5.6)-(5.7).
  # Category: rounded comparison. PS = Phi(d / sqrt(2)); the book uses .707
  # for 1 / sqrt(2) and prints .76 (".760") for d = 1 and .28 for d = -.84.
  # cles() uses 1 / sqrt(2) itself, as the book recommends for more
  # accuracy. A discrepancy in the source: for d = 1 the book rounds .707
  # to .71 and reads "approximately .760" below z = .71, but Phi(.71) =
  # .7611; .760 is Phi(.707) = .7602 (and Phi(1 / sqrt(2)) = .76025).
  # Because that print comes from a table lookup at a rounded z, it is not
  # compared; the value for d = 1 is pinned to Phi(1 / sqrt(2)) only.
  minus <- cles(smd = -0.84)
  plus <- cles(smd = 1)
  expect_equal(minus$value[minus$term == "cles"], pnorm(-0.84 / sqrt(2)),
               tolerance = 1e-15)
  expect_equal(minus$value[minus$term == "cles"], 0.276266190786987,
               tolerance = 1e-13)
  expect_equal(plus$value[plus$term == "cles"], 0.760249938906523,
               tolerance = 1e-13)
  expect_identical(round(minus$value[minus$term == "cles"], 2), 0.28)
  expect_identical(round(0.707 * -0.84, 5), -0.59388)
})

# ---- Efron and Tibshirani: the spatial test data, Table 14.1 -----------

efron_A <- c(48, 36, 20, 29, 42, 42, 20, 42, 22, 41, 45, 14, 6,
             0, 33, 28, 34, 4, 32, 24, 47, 41, 24, 26, 30, 41)
efron_B <- c(42, 33, 16, 39, 38, 36, 15, 33, 20, 43, 34, 22, 7,
             15, 34, 29, 41, 13, 38, 25, 27, 41, 28, 14, 28, 40)

# The plug-in variance (Eq. 14.1) and the BCa quantities of pp. 185-186,
# coded from the book: the acceleration of Eq. (14.15), with the
# leave-one-out values centered at their mean; the bias correction of
# Eq. (14.14) for 1116 of 2000 replications below the estimate; and the
# adjusted tail positions of Eq. (14.10) at the 90 percent level.
efron_theta <- function(z) mean((z - mean(z))^2)
efron_jack <- vapply(seq_along(efron_A), function(i) efron_theta(efron_A[-i]),
                     numeric(1))
efron_u <- mean(efron_jack) - efron_jack
efron_a <- sum(efron_u^3) / (6 * sum(efron_u^2)^1.5)
efron_z0 <- qnorm(1116 / 2000)
efron_zz <- qnorm(c(0.05, 0.95))
efron_pos <- pnorm(efron_z0 + (efron_z0 + efron_zz) /
                     (1 - efron_a * (efron_z0 + efron_zz)))

test_that("Efron and Tibshirani p. 179, Eqs. (14.1)-(14.2): plug-in and unbiased variance", {
  # Record: Efron and Tibshirani, p. 179, Eqs. (14.1)-(14.2); p. 180,
  # Table 14.1.
  # Category: rounded comparison of the inputs to the acceleration below.
  # n = 26. The plug-in estimate divides by n (printed 171.5); the usual
  # unbiased estimate divides by n - 1 (printed 178.4).
  plug_in <- efron_theta(efron_A)
  expect_equal(plug_in, (sum(efron_A^2) - 771^2 / 26) / 26, tolerance = 1e-14)
  expect_equal(plug_in, 171.534023668639, tolerance = 1e-13)
  expect_equal(stats::var(efron_A), 178.395384615385, tolerance = 1e-13)
  expect_identical(round(plug_in, 1), 171.5)
  expect_identical(round(stats::var(efron_A), 1), 178.4)
})

test_that("Efron and Tibshirani p. 186, Eq. (14.15): the arithmetic of .ratio_acceleration()", {
  # Record: Efron and Tibshirani, p. 186, Eq. (14.15).
  # Category: rounded comparison of the acceleration (.0612401, printed
  # .061), and algorithm oracle for the arithmetic of DMAR's
  # .ratio_acceleration() only. The centered values u are computed here
  # and passed in, so this checks the ratio of cubes, sum(u^3) /
  # (6 (sum u^2)^1.5), not a DMAR jackknife; DMAR's own jackknife is
  # checked against Eq. (14.15) in the test of .kappa_boot_limits() below.
  # (DMAR's ratio functions center their influence values at the
  # full-sample estimate, as boot's usual jackknife does; for this
  # statistic that gives .0641402, and for the ratio of the means of A and
  # B .0100678 against .0100692.)
  expect_equal(efron_a, 0.0612401198123075, tolerance = 1e-12)
  expect_identical(round(efron_a, 3), 0.061)
  expect_equal(DMAR:::.ratio_acceleration(list(efron_u), 26), efron_a,
               tolerance = 1e-14)
})

test_that("Efron and Tibshirani pp. 185-186, Eqs. (14.10), (14.14), (14.15): DMAR's kappa BCa engine", {
  # Record: Efron and Tibshirani, pp. 185-186, Eqs. (14.10), (14.13),
  # (14.14), and (14.15).
  # Category: algorithm oracle, deterministic. .kappa_boot_limits()
  # (R/kappa_internals.R, used by cohen_kappa() and fleiss_kappa()) computes
  # its own leave-one-out values through the function it is given, centers
  # them at their mean, and applies the BCa adjustment. Given the plug-in
  # variance of A as that function, the replicates 1, ..., 2000, and an
  # estimate with 1116 replicates below it, it must land at the positions
  # of Eq. (14.10) for the book's z0 = qnorm(.558) and a = .0612401
  # (printed (.110, .985) in Eq. 14.13). A type 7 quantile of 1:2000 at
  # position p is 1 + 1999 p, which returns the positions exactly. A
  # jackknife centered at the estimate instead would move the positions by
  # .001.
  #
  # Two estimates: 1116.5, which no replicate equals, and 1117, which one
  # replicate equals. Eq. (14.14) counts the replicates strictly below the
  # estimate, so 1116 count in both cases and the positions are the same.
  for (est in c(1116.5, 1117)) {
    out <- DMAR:::.kappa_boot_limits(as.numeric(1:2000), 2000L, est,
                                     function(i) efron_theta(efron_A[i]),
                                     26L, 0.90, "bca", "test")
    expect_true(out$verdict$formed)
    expect_equal((out$limits - 1) / 1999,
                 c(0.109906834295439, 0.984506078796717), tolerance = 1e-12)
    expect_equal((out$limits - 1) / 1999, efron_pos, tolerance = 1e-12)
  }
})

test_that("Efron and Tibshirani pp. 185-186, Eqs. (14.10), (14.13), (14.14): DMAR's ratio BCa positions", {
  # Record: Efron and Tibshirani, pp. 185-186, Eqs. (14.10), (14.13),
  # and (14.14).
  # Category: algorithm oracle, deterministic. With 1116 of 2000
  # replications below the estimate, z0 = qnorm(.558) = .146 (Eq. 14.14),
  # and with a = .0612401 the 90 percent BCa interval reads its limits at
  # alpha_1 = .1099068 and alpha_2 = .9845061 (Eq. 14.10; printed (.110,
  # .985) in Eq. 14.13, the 220th and 1970th of the 2000 ordered values).
  # DMAR's BCa step for the ratio functions is run here on the replicates
  # 1, ..., 2000 with the book's acceleration supplied, at the estimates
  # 1116.5 (no replicate equals it) and 1117 (one replicate equals it; the
  # strict count of Eq. 14.14 still gives 1116, and a count of
  # replicates at or below the estimate would move the positions by
  # .0004). A type 7 quantile of 1:2000 at position p is 1 + 1999 p, which
  # returns the positions exactly. The adjusted positions need
  # 10 / (1 - .9845061) = 646 replications, fewer than the 2000 requested,
  # so the interval is formed. The comparisons with the printed .146,
  # (.110, .985), and the 220th and 1970th ordered values are rounded
  # comparisons.
  expect_equal(efron_z0, 0.145900420032994, tolerance = 1e-12)
  expect_identical(round(efron_z0, 3), 0.146)
  expect_equal(efron_pos, c(0.109906834295439, 0.984506078796717),
               tolerance = 1e-12)
  rule <- DMAR:::.ratio_boot_rule(0.05, 0.05, 0.90, "bca")
  kept <- DMAR:::.ratio_boot_keep(as.numeric(1:2000), rep(1, 2000), 2000L,
                                  rule, "test")
  for (est in c(1116.5, 1117)) {
    lims <- DMAR:::.ratio_boot_limits(kept, est, 0.05, 0.05, "bca",
                                      acc = function() efron_a)$limits
    expect_equal((lims - 1) / 1999, efron_pos, tolerance = 1e-12)
  }
  expect_identical(round(efron_pos, 3), c(0.110, 0.985))
  expect_identical(round(2000 * round(efron_pos, 3)), c(220, 1970))
})

test_that("Efron and Tibshirani Table 14.1, adapted: paired ratio of the means of A and B", {
  # Record: Efron and Tibshirani, p. 180, Table 14.1 (data only).
  # Category: adaptation. The book's spatial test data used for an analysis
  # the book does not report: the ratio of the mean of test A to the mean
  # of test B over the 26 children, a paired design. This is not a
  # reproduction of any interval in the book (its Table 14.2 intervals are
  # for var(A)).
  #
  # Bootstrap record. Statistic: sum(A*) / sum(B*), the ratio of the means
  # of the resampled pairs (no divisor enters a ratio of means). Resampling
  # unit: the child, (A_i, B_i) kept together. Algorithm: nonparametric
  # bootstrap; B = 2000 replications drawn in one call of
  # sample.int(26, 26 * 2000, replace = TRUE) after the seed; BCa with bias
  # correction qnorm(mean(ratio* < ratio)) and the jackknife acceleration
  # over children centered at the estimate. Seed: 1988 (set for the call
  # and restored afterward). Quantile rule: stats::quantile() type 7 at the
  # adjusted positions (BCa) or at .05 and .95 (percentile). Level: 90
  # percent, the level of the book's Table 14.2.
  #
  # The seeded limits are pinned to a re-implementation that draws the
  # same stream (a check of the arithmetic). The Monte Carlo tolerance is
  # four times the standard deviation of a B = 2000 limit over 400 seeds
  # (BCa .0032 lower and .0037 upper; percentile .0023 and .0027), around
  # the long-run limits (BCa .94533 and 1.11405, where the mean over those
  # 400 seeds and boot::boot.ci() with R = 200,000 agree to within .00002;
  # percentile .94322 and 1.11160, the mean over the 400 seeds).
  bca <- ci_ratio_of_means_paired(efron_A, efron_B, ci_method = "bca",
                                  conf_level = 0.90, B = 2000, seed = 1988)
  pct <- ci_ratio_of_means_paired(efron_A, efron_B, ci_method = "percentile",
                                  conf_level = 0.90, B = 2000, seed = 1988)
  val <- function(x, t) x$value[x$term == t]
  expect_equal(val(bca, "ratio"), 771 / 751, tolerance = 1e-15)
  expect_equal(val(bca, "ratio"), 1.02663115845539, tolerance = 1e-13)
  expect_equal(attr(bca, "B_used"), 2000)
  expect_identical(attr(bca, "bootstrap_status"), "complete")

  withr::with_seed(1988, {
    idx <- sample.int(26L, 26L * 2000L, replace = TRUE)
  })
  rb <- colSums(matrix(efron_A[idx], nrow = 26)) /
    colSums(matrix(efron_B[idx], nrow = 26))
  est <- sum(efron_A) / sum(efron_B)
  loo <- (sum(efron_A) - efron_A) / (sum(efron_B) - efron_B)
  l <- 25 * (est - loo)
  a_hat <- sum(l^3) / (6 * sum(l^2)^1.5)
  z0 <- qnorm(mean(rb < est))
  zz <- qnorm(c(0.05, 0.95))
  pos <- pnorm(z0 + (z0 + zz) / (1 - a_hat * (z0 + zz)))
  expect_equal(c(val(bca, "ci_lower"), val(bca, "ci_upper")),
               stats::quantile(rb, pos, names = FALSE), tolerance = 1e-12)
  expect_equal(c(val(pct, "ci_lower"), val(pct, "ci_upper")),
               stats::quantile(rb, c(0.05, 0.95), names = FALSE),
               tolerance = 1e-12)
  expect_equal(c(val(bca, "ci_lower"), val(bca, "ci_upper")),
               c(0.943897102158344, 1.11428988118469), tolerance = 1e-10)
  expect_equal(c(val(pct, "ci_lower"), val(pct, "ci_upper")),
               c(0.941168585396625, 1.11096651445967), tolerance = 1e-10)
  expect_lt(abs(val(bca, "ci_lower") - 0.94533), 4 * 0.0032)
  expect_lt(abs(val(bca, "ci_upper") - 1.11405), 4 * 0.0037)
  expect_lt(abs(val(pct, "ci_lower") - 0.94322), 4 * 0.0023)
  expect_lt(abs(val(pct, "ci_upper") - 1.11160), 4 * 0.0027)
})
