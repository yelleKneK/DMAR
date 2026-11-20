# Differences with exactly the moments of the Bland and Altman (1986)
# example that Carkeet (2015) reanalyzes: n = 17, mean -2.1, SD 38.8.
# scale() centers and scales exactly, so the moments hold to machine
# precision and the published anchors below are exact benchmarks.
carkeet_differences <- function() {
  set.seed(113)
  z <- as.numeric(scale(rnorm(17)))
  z * 38.8 - 2.1
}

test_that("limits_of_agreement() returns mean, SD, and the two LoA with their CIs", {
  set.seed(113)
  a <- rnorm(40, 100, 15)
  b <- a + rnorm(40, 0, 3)
  res <- limits_of_agreement(a, b)
  expect_setequal(res$term,
                   c("mean_difference", "sd_difference",
                     "loa_lower", "loa_lower_ci_lower",
                     "loa_lower_ci_upper",
                     "loa_upper", "loa_upper_ci_lower",
                     "loa_upper_ci_upper"))
})

test_that("limits_of_agreement() 95% LoA = mean +/- 1.96 SD", {
  set.seed(113)
  a <- rnorm(40, 100, 15); b <- a + rnorm(40, 0, 3)
  res <- limits_of_agreement(a, b)
  d_m <- res$value[res$term == "mean_difference"]
  d_s <- res$value[res$term == "sd_difference"]
  expect_equal(res$value[res$term == "loa_upper"], d_m + qnorm(0.975) * d_s,
               tolerance = 1e-10)
  expect_equal(res$value[res$term == "loa_lower"], d_m - qnorm(0.975) * d_s,
               tolerance = 1e-10)
})

test_that("limits_of_agreement() CI bounds bracket the point LoA under both methods", {
  set.seed(113)
  a <- rnorm(40, 100, 15); b <- a + rnorm(40, 0, 3)
  for (m in c("pair", "individual")) {
    res <- expect_no_warning(limits_of_agreement(a, b, method = m))
    loa_l <- res$value[res$term == "loa_lower"]
    loa_u <- res$value[res$term == "loa_upper"]
    expect_lte(res$value[res$term == "loa_lower_ci_lower"], loa_l)
    expect_gte(res$value[res$term == "loa_lower_ci_upper"], loa_l)
    expect_lte(res$value[res$term == "loa_upper_ci_lower"], loa_u)
    expect_gte(res$value[res$term == "loa_upper_ci_upper"], loa_u)
  }
})

test_that("loa() is the sanctioned short alias of limits_of_agreement()", {
  expect_identical(loa, limits_of_agreement)
  set.seed(113)
  a <- rnorm(30, 100, 15); b <- a + rnorm(30, 0, 4)
  expect_identical(loa(a, b), limits_of_agreement(a, b))
})

test_that("the default method is 'pair' and the method is recorded", {
  set.seed(113)
  a <- rnorm(20, 50, 8); b <- a + rnorm(20, 1, 2)
  res_default <- limits_of_agreement(a, b)
  expect_identical(res_default, limits_of_agreement(a, b, method = "pair"))
  expect_identical(attr(res_default, "method"), "pair")
  res_ind <- limits_of_agreement(a, b, method = "individual")
  expect_identical(attr(res_ind, "method"), "individual")
  expect_error(limits_of_agreement(a, b, method = "both"))
})

test_that("limits_of_agreement(method = 'individual') CI limits follow the documented Carkeet noncentral t quantiles", {
  set.seed(113)
  x <- rnorm(10); y <- x + rnorm(10, 0, 2)
  # Through 3338679 qt() warned here (full precision may not have been
  # achieved) for the negative noncentrality, and the warning was muffled.
  res <- expect_no_warning(
    limits_of_agreement(x, y, coverage = 0.90, conf_level = 0.90, method = "individual"))
  v <- function(t) res$value[res$term == t]
  d <- y - x; n <- 10
  d_m <- mean(d); d_s <- sd(d)
  k <- stats::qnorm(0.95); delta <- k * sqrt(n)
  q <- function(p, ncp) suppressWarnings(stats::qt(p, df = n - 1, ncp = ncp))
  expect_equal(v("loa_upper_ci_lower"),
               d_m + d_s / sqrt(n) * q(0.05, delta), tolerance = 1e-10)
  expect_equal(v("loa_upper_ci_upper"),
               d_m + d_s / sqrt(n) * q(0.95, delta), tolerance = 1e-10)
  expect_equal(v("loa_lower_ci_lower"),
               d_m + d_s / sqrt(n) * q(0.05, -delta), tolerance = 1e-10)
  expect_equal(v("loa_lower_ci_upper"),
               d_m + d_s / sqrt(n) * q(0.95, -delta), tolerance = 1e-10)

  # Published anchor: (upper CI limit of the upper LoA - mean) / SD is
  # the one-sided normal tolerance factor, tabled as 2.911 for n = 10,
  # 95% content, 95% confidence (Odeh & Owen, 1980; Hahn & Meeker,
  # 1991). The symmetric central-t form previously printed on the page
  # cannot reproduce this value.
  expect_equal((v("loa_upper_ci_upper") - d_m) / d_s, 2.911,
               tolerance = 1e-3)

  # The exact individual CI is asymmetric about the sample LoA, wider
  # on the side away from the mean difference.
  expect_gt(v("loa_upper_ci_upper") - v("loa_upper"),
            v("loa_upper") - v("loa_upper_ci_lower"))
})

test_that("limits_of_agreement(method = 'individual') reproduces Carkeet's (2015) worked example and Table 1", {
  d <- carkeet_differences()
  res <- limits_of_agreement(rep(0, 17), d, method = "individual")
  v <- function(t) res$value[res$term == t]
  expect_equal(v("mean_difference"), -2.1, tolerance = 1e-12)
  expect_equal(v("sd_difference"), 38.8, tolerance = 1e-12)
  expect_equal(v("loa_upper"), -2.1 + qnorm(0.975) * 38.8, tolerance = 1e-10)

  # Carkeet (2015), pp. e74-e75: the exact CI on the upper LoA is
  # [48.9, 120.0]; on the lower LoA, [-124.2, -53.1].
  expect_equal(v("loa_upper_ci_lower"), 48.9, tolerance = 1e-3)
  expect_equal(v("loa_upper_ci_upper"), 120.0, tolerance = 1e-3)
  expect_equal(v("loa_lower_ci_lower"), -124.2, tolerance = 1e-3)
  expect_equal(v("loa_lower_ci_upper"), -53.1, tolerance = 1e-3)

  # Table 1 coefficients at nu = 16, F = 0.025 / 0.975: 1.3150 / 3.1483.
  expect_equal((v("loa_upper_ci_lower") - v("mean_difference")) / 38.8,
               1.3150, tolerance = 1e-4)
  expect_equal((v("loa_upper_ci_upper") - v("mean_difference")) / 38.8,
               3.1483, tolerance = 1e-4)
})

test_that("limits_of_agreement(method = 'pair') reproduces Carkeet's (2015) worked example and Table 2", {
  d <- carkeet_differences()
  res <- limits_of_agreement(rep(0, 17), d)  # pair is the default
  expect_identical(attr(res, "method"), "pair")
  v <- function(t) res$value[res$term == t]
  d_m <- v("mean_difference")

  # Carkeet (2015), p. e76: the pair bounds are -2.1 +/- 57.81 (inner)
  # and -2.1 +/- 119.60 (outer).
  expect_equal(v("loa_upper_ci_lower") - d_m,   57.81, tolerance = 1e-4)
  expect_equal(v("loa_upper_ci_upper") - d_m,  119.60, tolerance = 1e-4)
  expect_equal(v("loa_lower_ci_upper") - d_m,  -57.81, tolerance = 1e-4)
  expect_equal(v("loa_lower_ci_lower") - d_m, -119.60, tolerance = 1e-4)

  # Table 2 coefficients at nu = 16, F = 0.025 / 0.975: 1.4900 / 3.0824.
  expect_equal((v("loa_upper_ci_lower") - d_m) / 38.8, 1.4900,
               tolerance = 1e-4)
  expect_equal((v("loa_upper_ci_upper") - d_m) / 38.8, 3.0824,
               tolerance = 1e-4)

  # The pair CI bounds are symmetric about the mean difference, not
  # about the sample LoA.
  expect_equal(v("loa_upper_ci_upper") - d_m,
               d_m - v("loa_lower_ci_lower"), tolerance = 1e-10)
  expect_equal(v("loa_upper_ci_lower") - d_m,
               d_m - v("loa_lower_ci_upper"), tolerance = 1e-10)
})

test_that("limits_of_agreement(method = 'individual') is exact past a noncentrality of 37.62", {
  # ncp = 1.96 sqrt(n) passes 37.62 at n = 369. Here n = 1000, and qt()
  # (a normal approximation there; DMAR through 3338679) put the interval
  # on the upper LoA at [6.051969469915, 6.686201336532] and the one on the
  # lower LoA at [-5.681620574917, -5.047388708300], each about 1e-3 too
  # far from the mean difference (the differences have an SD of 3).
  set.seed(1000)
  a <- rnorm(1000, 100, 15); b <- a + rnorm(1000, 0.5, 3)
  res <- expect_no_warning(limits_of_agreement(a, b, method = "individual"))
  v <- function(t) res$value[res$term == t]
  expect_equal(c(v("loa_upper_ci_lower"), v("loa_upper_ci_upper")),
               c(6.051190025909, 6.685130194723), tolerance = 1e-10)
  expect_equal(c(v("loa_lower_ci_lower"), v("loa_lower_ci_upper")),
               c(-5.680549433108, -5.046609264294), tolerance = 1e-10)
  # The tail probabilities at the limits, by the oracle.
  d <- b - a; n <- 1000; ncp <- qnorm(.975) * sqrt(n)
  t_of <- function(limit) (limit - mean(d)) * sqrt(n) / sd(d)
  expect_equal(.oracle_nct(t_of(v("loa_upper_ci_upper")), n - 1, ncp, lower_tail = FALSE),
               .025, tolerance = 1e-9)
  expect_equal(.oracle_nct(t_of(v("loa_upper_ci_lower")), n - 1, ncp), .025, tolerance = 1e-9)
  expect_equal(.oracle_nct(t_of(v("loa_lower_ci_lower")), n - 1, -ncp), .025, tolerance = 1e-9)
})
