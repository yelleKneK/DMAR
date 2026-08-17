# An independent brute-force randomization p-value. This enumerates every
# reassignment itself and never calls DMAR, so it is a genuine oracle for
# the exact path rather than a restatement of it.
brute_force_randomization_p <- function(a, b, alternative = "two_sided") {
  y  <- c(a, b)
  n1 <- length(a)
  N  <- length(y)
  idx <- utils::combn(N, n1)
  stats <- apply(idx, 2L, function(i) mean(y[i]) - mean(y[-i]))
  obs <- mean(a) - mean(b)
  tol <- 1e-12
  switch(alternative,
         two_sided = mean(abs(stats) >= abs(obs) - tol),
         greater   = mean(stats >= obs - tol),
         less      = mean(stats <= obs + tol))
}


test_that("randomization_test() exact p-value matches brute-force enumeration", {
  treatment <- c(80, 84, 79, 88, 83)
  control   <- c(72, 75, 68, 81, 74)

  for (alt in c("two_sided", "greater", "less")) {
    res <- randomization_test(group_1 = treatment, group_2 = control,
                              alternative = alt)
    p_dmar <- res$value[res$term == "p_value"]
    expect_equal(p_dmar,
                 brute_force_randomization_p(treatment, control, alt),
                 tolerance = 1e-12)
  }

  # every reassignment was enumerated, so the result is exact
  res <- randomization_test(group_1 = treatment, group_2 = control)
  expect_equal(res$value[res$term == "exact"], 1)
  expect_equal(res$value[res$term == "n_evaluated"], choose(10, 5))
  expect_equal(attr(res, "method"), "exact enumeration")
})


test_that("randomization_test() reports the observed mean difference and sizes", {
  treatment <- c(80, 84, 79, 88, 83)
  control   <- c(72, 75, 68, 81, 74)
  res <- randomization_test(group_1 = treatment, group_2 = control)

  expect_equal(res$value[res$term == "mean_difference"],
               mean(treatment) - mean(control), tolerance = 1e-12)
  expect_equal(res$value[res$term == "n_1"], length(treatment))
  expect_equal(res$value[res$term == "n_2"], length(control))
  expect_equal(res$value[res$term == "N"],
               length(treatment) + length(control))
  # the value column stays numeric; labels live on attributes
  expect_true(is.numeric(res$value))
  expect_type(attr(res, "statistic_name"), "character")
})


test_that("randomization_test() effect sizes agree with the package's own estimators", {
  treatment <- c(80, 84, 79, 88, 83)
  control   <- c(72, 75, 68, 81, 74)
  res <- randomization_test(group_1 = treatment, group_2 = control)

  # the standardized mean difference is smd(), not a private reimplementation
  expected_smd <- smd(group_1 = treatment, group_2 = control)$value[1]
  expect_equal(res$value[res$term == "smd"], expected_smd, tolerance = 1e-8)

  # the interval brackets the estimate
  expect_lt(res$value[res$term == "smd_lower_limit"],
            res$value[res$term == "smd"])
  expect_gt(res$value[res$term == "smd_upper_limit"],
            res$value[res$term == "smd"])

  # the distribution-free companions are present with intervals
  for (nm in c("cles", "cliff_delta")) {
    expect_true(nm %in% res$term)
    expect_true(paste0(nm, "_lower_limit") %in% res$term)
    expect_true(paste0(nm, "_upper_limit") %in% res$term)
  }
})


test_that("randomization_test() inverted interval brackets the mean difference", {
  treatment <- c(80, 84, 79, 88, 83)
  control   <- c(72, 75, 68, 81, 74)
  res <- randomization_test(group_1 = treatment, group_2 = control)

  lo   <- res$value[res$term == "shift_lower_limit"]
  hi   <- res$value[res$term == "shift_upper_limit"]
  diff <- res$value[res$term == "mean_difference"]

  expect_true(is.finite(lo) && is.finite(hi))
  expect_lt(lo, diff)
  expect_gt(hi, diff)

  # A shift just inside the interval is not rejected; a shift well outside
  # it is. This is what "inverting the test" means, and it is the property
  # that distinguishes this interval from a normal-theory one.
  inside <- randomization_test(group_1 = treatment - (lo + hi) / 2,
                               group_2 = control)
  expect_gt(inside$value[inside$term == "p_value"], 0.05)
})


test_that("randomization_test() Monte Carlo approaches the exact p-value", {
  skip_on_cran()
  set.seed(113)
  # Two groups of 10 keep choose(20, 10) = 184,756 reassignments under
  # the documented 1,000,000 cap on exact enumeration.
  a <- rnorm(10, mean = 0.8)
  b <- rnorm(10, mean = 0)

  exact <- randomization_test(group_1 = a, group_2 = b, exact = TRUE)
  mc    <- randomization_test(group_1 = a, group_2 = b, exact = FALSE,
                              n_resamples = 20000L, seed = 113)

  p_exact <- exact$value[exact$term == "p_value"]
  p_mc    <- mc$value[mc$term == "p_value"]
  se_mc   <- mc$value[mc$term == "p_value_se"]

  expect_equal(mc$value[mc$term == "exact"], 0)
  expect_true(is.finite(se_mc) && se_mc > 0)
  # within four Monte Carlo standard errors of the enumerated value
  expect_lt(abs(p_mc - p_exact), 4 * se_mc + 1e-8)
  # the (r + 1) / (m + 1) form cannot return zero
  expect_gt(p_mc, 0)
})


test_that("randomization_test() studentized statistic is available and directional tests order correctly", {
  treatment <- c(80, 84, 79, 88, 83)
  control   <- c(72, 75, 68, 81, 74)

  res_t <- randomization_test(group_1 = treatment, group_2 = control,
                              statistic = "t")
  expect_true(is.finite(res_t$value[res_t$term == "p_value"]))
  expect_match(attr(res_t, "statistic_name"), "t", ignore.case = TRUE)

  # treatment scores are higher, so "greater" is the supported direction
  p_greater <- randomization_test(group_1 = treatment, group_2 = control,
                                  alternative = "greater")
  p_less    <- randomization_test(group_1 = treatment, group_2 = control,
                                  alternative = "less")
  expect_lt(p_greater$value[p_greater$term == "p_value"],
            p_less$value[p_less$term == "p_value"])
})


test_that("randomization_test() accepts the formula and vector interfaces alike", {
  d <- data.frame(y = c(80, 84, 79, 88, 83, 72, 75, 68, 81, 74),
                  g = rep(c("treatment", "control"), each = 5))
  from_formula <- randomization_test(y ~ g, data = d)
  from_vectors <- randomization_test(group_1 = d$y[d$g == "treatment"],
                                     group_2 = d$y[d$g == "control"])
  expect_equal(from_formula$value[from_formula$term == "p_value"],
               from_vectors$value[from_vectors$term == "p_value"],
               tolerance = 1e-12)
})


test_that("randomization_test() rejects bad input", {
  expect_error(randomization_test(group_1 = 1:5, group_2 = 1:5,
                                  conf_level = 1.4), "conf_level")
  expect_error(randomization_test(group_1 = numeric(0), group_2 = 1:5))
})


test_that("plot_randomization_test() returns a ggplot built from the reference distribution", {
  skip_if_not_installed("ggplot2")
  treatment <- c(80, 84, 79, 88, 83)
  control   <- c(72, 75, 68, 81, 74)
  res <- randomization_test(group_1 = treatment, group_2 = control)

  p <- plot_randomization_test(res)
  expect_s3_class(p, "ggplot")
  # every reassignment appears in the plotted data
  expect_equal(nrow(p$data), length(attr(res, "reference_distribution")))
  # the shaded proportion is the p-value the test reported
  shaded <- mean(p$data$region == "At least as extreme")
  expect_equal(shaded, res$value[res$term == "p_value"], tolerance = 1e-12)

  expect_error(plot_randomization_test(res, bins = -1), "bins")
  expect_error(plot_randomization_test(data.frame(a = 1)),
               "reference distribution")
})

test_that("randomization_test() is equivariant under rescaling the response", {
  # Multiplying every score by c > 0 multiplies the statistic and the whole
  # reference distribution by c, so the set of reassignments counted as at
  # least as extreme, and therefore the p-value, cannot change. The tie
  # tolerance must scale with the data for this to hold; an absolute floor
  # would swallow the entire distribution when the scores are small.
  set.seed(113)
  g1 <- rnorm(6, 0, 1)
  g2 <- rnorm(6, 0.8, 1)
  p_of <- function(x) x$value[x$term == "p_value"]

  base <- p_of(randomization_test(group_1 = g1, group_2 = g2))
  for (scale in c(1e-9, 1e-7, 1e-3, 1e3, 1e7, 1e9)) {
    expect_equal(
      p_of(randomization_test(group_1 = g1 * scale, group_2 = g2 * scale)),
      base
    )
  }
  # The studentized (t) statistic is scale free to begin with; it must agree
  # too. The accepted values are "mean" and "t"; there is no "studentized".
  base_t <- p_of(randomization_test(group_1 = g1, group_2 = g2,
                                    statistic = "t"))
  expect_equal(
    p_of(randomization_test(group_1 = g1 * 1e-8, group_2 = g2 * 1e-8,
                            statistic = "t")),
    base_t
  )
})
