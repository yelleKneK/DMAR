## The mean and variance of a standard deviation under normality
## (R/c4_internals.R): for W chi square on df degrees of freedom,
## E(S) / sigma = E(sqrt(W / df)) = sqrt(2 / df) Gamma((df + 1)/2) / Gamma(df/2)
## and Var(S) / sigma^2 = 1 - {E(S) / sigma}^2; c4(n) is the mean factor at
## df = n - 1. ss_aipe_mean(), ss_aipe_two_phase(), design_consequences(),
## sd_unbiased(), var_ete(), and expected_r() (through it,
## expected_partial_r()) use them. Both factors are formed from
## L = log{E(S) / sigma}: from df = 100 on its asymptotic expansion, and
## below that the expansion carried down by the exact recurrence
## L(df) = L(df + 2) - log1p(1 / {df (df + 2)}) / 2.
##
## Oracle. The pinned values were computed with Rmpfr at 320 bits, as L,
## exp(L), and -expm1(2 L) with L = log(2 / df) / 2 + lgamma((df + 1) / 2) -
## lgamma(df / 2), and are given to 17 significant digits. The code agrees
## with them to less than one machine epsilon (2.2e-16) at every pin, and
## to a few units in the last place (5.5e-16, relative, the largest found)
## over some 490,000 df from 0.05 to 1e15 (see R/c4_internals.R);
## the relative tolerance of 1e-15 below leaves room for a platform whose
## exp(), expm1(), and log1p() round differently.

pin_df <- c(0.5, 1, 8, 39, 58, 99, 100, 118, 500, 5000, 1e5, 1e6, 1e8, 1e12,
            1e13, 1e14, 1.7e14, 1e15)
pin_mean <- c(0.67597824006728473, 0.79788456080286536, 0.96931069971395408,
              0.99361094283188581, 0.99569914393430970, 0.99747797607126351,
              0.99750316395510509, 0.99788362397817775, 0.99950012531233438,
              0.99995000125031248, 0.99999750000312504, 0.99999975000003125,
              0.99999999750000000, 0.99999999999975000, 0.99999999999997500,
              0.99999999999999750, 0.99999999999999853, 0.99999999999999975)
pin_var <- c(0.54305341895553637, 0.36338022763241866, 0.060436767420044747,
             0.012737294284730944, 0.0085832147684828248,
             0.0050376872527758590, 0.0049874378995547391,
             0.0042282729961787637, 0.00099949950062787165,
             9.9994999500062529e-5, 4.9999874999375004e-6,
             4.9999987499993750e-7, 4.9999999874999999e-9,
             4.9999999999987500e-13, 4.9999999999998750e-14,
             4.9999999999999875e-15, 2.9411764705882310e-15,
             4.9999999999999987e-16)
pin_log <- c(-0.39159439270683678, -0.22579135264472743, -0.031170078945558090,
             -0.0064095545464535145, -0.0043101313511378138,
             -0.0025252095884189043, -0.0024999583383318163,
             -0.0021186187103619850, -0.00049999966666826665,
             -4.9999999666666683e-5, -2.4999999999583333e-6,
             -2.4999999999995833e-7, -2.5000000000000000e-9,
             -2.5000000000000000e-13, -2.5000000000000000e-14,
             -2.5000000000000000e-15, -1.4705882352941176e-15,
             -2.5000000000000000e-16)

# expect_equal() compares on the absolute scale whenever the expected value
# is smaller than the tolerance, as all.equal() does, so the pins of 5e-16
# and -2.5e-16 at df = 1e15, checked with tolerance = 1e-15, would pass for
# any value near 0. Every pin is checked on its relative error instead.
expect_relative <- function(object, expected, tolerance, what) {
  expect_lt(abs(object / expected - 1), tolerance,
            label = paste("the relative error of", what))
}

test_that("the chi mean factor is below 1 and increases with df from 1 to 1e15", {
  # A tenth of a decade apart, consecutive values differ by more than 6e-14
  # up to df = 1e12, hundreds of units in the last place, so the factor
  # must increase strictly there; beyond 1e12 the steps shrink below a unit
  # in the last place, and the factor must not decrease. The difference of
  # two lgamma() values that design_consequences() used through 1.0.x fails
  # both conditions: on this grid it reaches 1 at df = 10^7.6. The ratio
  # formed through lbeta() comes out above 1 at some df from about 1.7e14
  # on.
  df <- 10^seq(0, 15, by = 0.1)
  f <- DMAR:::.chi_mean_factor(df)
  expect_true(all(f > 0 & f < 1))
  expect_true(all(diff(f[df <= 1e12]) > 0))
  expect_true(all(diff(f) >= 0))
  # A hundredth of a decade apart.
  f_fine <- DMAR:::.chi_mean_factor(10^seq(0, 15, by = 0.01))
  expect_true(all(f_fine < 1))
  expect_true(all(diff(f_fine) >= 0))
  # Every integer df up to 1000, where the steps are largest, and df below 1.
  f_int <- DMAR:::.chi_mean_factor(c(seq(0.05, 0.95, by = 0.05), 1:1000))
  expect_true(all(f_int > 0 & f_int < 1))
  expect_true(all(diff(f_int) > 0))
  # Beyond 1e15 the factor is never above 1, and it is exactly 1 only from
  # about df = 4.5e15 on, where the true factor, 1 - 1/(4 df), rounds to 1.
  df_big <- c(1e15, 4e15, 5e15, 1e16, 1e20, 1e100, 1e300, Inf)
  f_big <- DMAR:::.chi_mean_factor(df_big)
  expect_true(all(f_big <= 1))
  expect_true(all(diff(f_big) >= 0))
  expect_lt(DMAR:::.chi_mean_factor(4e15), 1)
  expect_identical(DMAR:::.chi_mean_factor(1e16), 1)
})

test_that("the chi factors are monotone where the recurrence changes its number of steps", {
  # From df = 100 on, L is the expansion; below, the recurrence carries it
  # down from df + 2m, the first such point at or above 100, so m changes
  # at every even df from 2 to 98, and the expansion takes over at 100.
  # Across each of those points, on grids 1e-13 apart, the mean factor and
  # L must not decrease nor the variance factor increase. At df = 100 one
  # step moves the mean factor by about 2.5e-18, far less than a unit in
  # the last place, so a jump of even one unit where two forms meet would
  # fail. The lbeta() ratio, used below 100 df in an earlier form, fails on
  # the grids around 20, 98, and 100 (more than 800 decreases on each), and
  # so does 1 minus its square.
  #
  # From one double to the next the factors are not exactly monotone, here
  # or anywhere a computed value carries a rounding error. Near df = 42 the
  # mean factor changes by about 1e-18 from one double to the next, a
  # hundredth of its unit in the last place, and an error of a fraction of
  # a unit steps it back by one unit at df = 42 + 1.13e-11; L and the
  # variance factor step back by one unit from the double below 62 to 62,
  # where m changes. Only a correctly rounded evaluation would rule this
  # out.
  for (center in seq(2, 100, by = 2)) {
    df <- center + (-1000:1000) * 1e-13
    expect_true(all(diff(DMAR:::.chi_mean_factor(df)) >= 0),
                info = paste("df near", center))
    expect_true(all(diff(DMAR:::.chi_log_mean_factor(df)) >= 0),
                info = paste("df near", center))
    expect_true(all(diff(DMAR:::.chi_var_factor(df)) <= 0),
                info = paste("df near", center))
  }
})

test_that("the chi mean factor is exact at small df and follows its expansion at large df", {
  # Closed forms: E|Z| = sqrt(2 / pi) at df = 1, E(chi_2) / sqrt(2) =
  # sqrt(pi) / 2 at df = 2, and E(chi_3) / sqrt(3) = 2 sqrt(2 / (3 pi)).
  expect_equal(DMAR:::.chi_mean_factor(1), sqrt(2 / pi), tolerance = 1e-15)
  expect_equal(DMAR:::.chi_mean_factor(2), sqrt(pi) / 2, tolerance = 1e-15)
  expect_equal(DMAR:::.chi_mean_factor(3), 2 * sqrt(2 / (3 * pi)),
               tolerance = 1e-15)
  # log{E(S) / sigma} = -1/(4 df) + 1/(24 df^3) - 1/(20 df^5) + O(df^-7),
  # the expansion of log Gamma(x + 1/2) - log Gamma(x) - log(x) / 2 at
  # x = df / 2; the omitted terms are below 1e-21 for df >= 1000. The
  # lgamma() difference is already off by 2e-12 at df = 1e4.
  for (df in 10^(3:15)) {
    expansion <- exp(-1 / (4 * df) + 1 / (24 * df^3) - 1 / (20 * df^5))
    expect_equal(DMAR:::.chi_mean_factor(df), expansion, tolerance = 1e-15,
                 info = paste("df =", df))
  }
  # The lbeta() ratio, right below df = 100 to a few units in the last
  # place (4e-15 at most, at df = 19.95, on a 0.01 grid).
  lbeta_form <- function(df) sqrt(2 / df) * exp(lgamma(0.5) - lbeta(df / 2, 0.5))
  df <- seq(0.05, 99.95, by = 0.05)
  expect_lt(max(abs(DMAR:::.chi_mean_factor(df) / lbeta_form(df) - 1)), 1e-14)
  for (i in seq_along(pin_df)) {
    expect_relative(DMAR:::.chi_mean_factor(pin_df[i]), pin_mean[i], 1e-15,
                    paste("the mean factor at df =", pin_df[i]))
  }
})

test_that("the chi variance factor is in (0, 1), decreases with df, and agrees with 320-bit values", {
  df <- 10^seq(0, 15, by = 0.1)
  v <- DMAR:::.chi_var_factor(df)
  expect_true(all(v > 0 & v < 1))
  expect_true(all(diff(v) < 0))
  v_int <- DMAR:::.chi_var_factor(c(seq(0.05, 0.95, by = 0.05), 1:1000))
  expect_true(all(v_int > 0 & v_int < 1))
  expect_true(all(diff(v_int) < 0))
  # The mean and the variance factor describe one distribution.
  expect_lt(max(abs(DMAR:::.chi_mean_factor(df)^2 + v - 1)), 1e-15)
  # Against the 320-bit values. One minus the square of the lgamma() mean
  # factor, which design_consequences() used through 1.0.x, is off by 2e-12
  # at df = 118; one minus the square of the lbeta() ratio is off by up to
  # 3e-13 below df = 100, and even with the mean factor right to rounding,
  # the subtraction is off by 1e-9 at df = 1e6.
  for (i in seq_along(pin_df)) {
    expect_relative(DMAR:::.chi_var_factor(pin_df[i]), pin_var[i], 1e-15,
                    paste("the variance factor at df =", pin_df[i]))
  }
  # Leading terms: Var(S) / sigma^2 = 1/(2 df) - 1/(8 df^2) + O(df^-3).
  for (big in c(1e8, 1e12, 1e15)) {
    expect_relative(DMAR:::.chi_var_factor(big), 1 / (2 * big) - 1 / (8 * big^2),
                    1e-15, paste("the variance factor at df =", big))
  }
  expect_gt(DMAR:::.chi_var_factor(1e300), 0)
})

test_that("the log of the chi mean factor keeps its relative precision at every df", {
  # expected_r() forms its bias from L itself. At df = 1e15, log() of even
  # a correctly rounded mean factor would keep only one significant digit
  # (-2.2e-16 for -2.5e-16, a relative error of 0.11).
  for (i in seq_along(pin_df)) {
    expect_relative(DMAR:::.chi_log_mean_factor(pin_df[i]), pin_log[i], 1e-15,
                    paste("L at df =", pin_df[i]))
  }
  expect_identical(DMAR:::.chi_mean_factor(pin_df),
                   exp(DMAR:::.chi_log_mean_factor(pin_df)))
})

test_that("the chi factors are vectorized and give NA and NaN where they should", {
  df <- c(3, NA, 250, 0, -1, 40)
  expect_identical(DMAR:::.chi_mean_factor(df)[c(1, 3, 6)],
                   DMAR:::.chi_mean_factor(c(3, 250, 40)))
  expect_true(is.na(DMAR:::.chi_mean_factor(df)[2]))
  expect_true(all(is.nan(DMAR:::.chi_mean_factor(df)[4:5])))
  expect_length(DMAR:::.chi_mean_factor(numeric(0)), 0L)
  expect_identical(DMAR:::.chi_mean_factor(Inf), 1)
  expect_identical(DMAR:::.chi_var_factor(Inf), 0)
})

test_that("c4(n) is the chi mean factor at df = n - 1", {
  n <- c(2:500, 10^(3:12))
  expect_identical(DMAR:::.c4(n), DMAR:::.chi_mean_factor(n - 1))
})

# ---- Hedges's J from the logarithm of the chi mean factor ----------------

# Reference values to 20 digits from a 60-digit evaluation of
# Gamma(df/2) / {sqrt(df/2) Gamma((df - 1)/2)} and of 1/J - 1.
hedges_df  <- c(2, 3, 7.25, 10, 50, 99, 99.5, 100, 101, 300, 1e5, 1e8, 1e12)
hedges_J   <- c(0.56418958354775628695, 0.72360125455826765936,
                0.8922185011154894354, 0.92274560805308710336,
                0.98491194229557260115, 0.99240185111442286017,
                0.99244014507826257364, 0.9924780549814448666,
                0.99255274548746187935, 0.99749756684385860828,
                0.99999249997812492969, 0.99999999249999997812,
                0.99999999999925)
hedges_inv <- c(0.7724538509055160273, 0.38197659788534191706,
                0.12080168562942548323, 0.083722307939143636462,
                0.015319194596484509532, 0.0076563227658682404827,
                0.0076174416756803663704, 0.0075789534900051385945,
                0.007503132248030435602, 0.0025087110378216144082,
                7.5000781258203206007e-6, 7.5000000781250008203e-9,
                7.5000000000078125e-13)

test_that("Hedges's J and 1/J - 1 are right to rounding from df = 2 to 1e12", {
  rel <- function(a, b) abs(a - b) / abs(b)
  expect_lt(max(rel(DMAR:::.hedges_J(hedges_df), hedges_J)), 1e-15)
  expect_lt(max(rel(expm1(-DMAR:::.hedges_log_J(hedges_df)), hedges_inv)),
            1e-15)
  # The difference of lgamma() values used through 1.0.x does not keep
  # 1/J - 1 at a large df: off by a factor of about 3 at df = 1e8.
  old <- function(df) exp(lgamma(df / 2) - lgamma((df - 1) / 2)) / sqrt(df / 2)
  expect_gt(rel(1 / old(1e8) - 1, hedges_inv[12]), 0.5)
})

test_that("Hedges's J equals the gamma ratio at df 2 to 300", {
  df <- 2:300
  by_lgamma <- exp(lgamma(df / 2) - lgamma((df - 1) / 2)) / sqrt(df / 2)
  expect_lt(max(abs(DMAR:::.hedges_J(df) / by_lgamma - 1)), 1e-12)
  # gamma() itself is finite to an argument of about 171; its ratio is
  # right to about 2e-13 there, so the tolerance is that of the ratio.
  small <- 2:340
  by_gamma <- gamma(small / 2) / (sqrt(small / 2) * gamma((small - 1) / 2))
  expect_lt(max(abs(DMAR:::.hedges_J(small) / by_gamma - 1)), 1e-12)
})

test_that("Hedges's J is below 1, increasing, vectorized, and 0 at df = 1", {
  df <- c(2:400, 10^(3:15))
  J <- DMAR:::.hedges_J(df)
  expect_true(all(J < 1 | df >= 1e15))
  expect_true(all(J <= 1))
  expect_true(all(diff(J) >= 0))
  expect_identical(DMAR:::.hedges_J(c(1, NA, 0.5, 10)),
                   c(0, NA, NaN, DMAR:::.hedges_J(10)))
  expect_identical(DMAR:::.hedges_log_J(1), -Inf)
})
