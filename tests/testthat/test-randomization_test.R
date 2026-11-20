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
  skip_on_cran()  # The extended validation grid runs in the full local suite.
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
  expect_lt(res$value[res$term == "smd_ci_lower"],
            res$value[res$term == "smd"])
  expect_gt(res$value[res$term == "smd_ci_upper"],
            res$value[res$term == "smd"])

  # the distribution-free companions are present with intervals
  for (nm in c("cles", "cliff_delta")) {
    expect_true(nm %in% res$term)
    expect_true(paste0(nm, "_ci_lower") %in% res$term)
    expect_true(paste0(nm, "_ci_upper") %in% res$term)
  }
})


test_that("randomization_test() inverted interval brackets the mean difference", {
  treatment <- c(80, 84, 79, 88, 83)
  control   <- c(72, 75, 68, 81, 74)
  res <- randomization_test(group_1 = treatment, group_2 = control)

  lo   <- res$value[res$term == "shift_ci_lower"]
  hi   <- res$value[res$term == "shift_ci_upper"]
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
  skip_on_cran()
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
  skip_on_cran()  # The extended validation grid runs in the full local suite.
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
  skip_on_cran()
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


# ---- Ties, boundary values, and the count of resamples ------------------
#
# The reference values below are exact. For the difference in means on
# scores that are whole multiples of a unit, the reassignments are ordered
# by whole-number arithmetic. For the studentized statistic, the p-values
# in fixtures/randomization_test_rational_oracle.csv were computed outside
# R in rational arithmetic on the decimal values.

rt_value <- function(x, term) x$value[x$term == term]

# Exact p-value of the difference in means for whole-number scores k_1 and
# k_2: n_2 * sum(first group) - n_1 * sum(second group) orders the
# reassignments as the mean difference does and is a whole number.
whole_number_randomization_p <- function(k_1, k_2, alternative) {
  k <- c(k_1, k_2)
  n_1 <- length(k_1)
  n_2 <- length(k_2)
  idx <- utils::combn(length(k), n_1)
  total <- sum(k)
  s_1 <- colSums(matrix(k[idx], nrow = n_1))
  stat <- n_2 * s_1 - n_1 * (total - s_1)
  obs  <- n_2 * sum(k_1) - n_1 * sum(k_2)
  switch(alternative,
         two_sided = mean(abs(stat) >= abs(obs)),
         greater   = mean(stat >= obs),
         less      = mean(stat <= obs))
}


test_that("randomization_test() counts every tie under the studentized statistic with two-valued data", {
  skip_on_cran()
  g_1 <- c(1, 1, 1, 1, 0)
  g_2 <- c(0, 0, 0, 0, 1)
  two <- randomization_test(group_1 = g_1, group_2 = g_2, statistic = "t")
  expect_equal(rt_value(two, "p_value"), 52 / 252, tolerance = 1e-12)
  one <- randomization_test(group_1 = g_1, group_2 = g_2, statistic = "t",
                            alternative = "greater")
  expect_equal(rt_value(one, "p_value"), 26 / 252, tolerance = 1e-12)
  # The reassignments that put the five ones in one group have no
  # within-group spread, and their statistic is infinite.
  expect_true(any(is.infinite(attr(two, "reference_distribution"))))
  expect_identical(sum(attr(two, "reference_extreme")), 52L)
})


test_that("randomization_test() refuses the studentized statistic for constant separated groups", {
  skip_on_cran()
  # Both groups constant at different values: the separate-variances
  # standard error is 0 and the studentized statistic has no value. Round 1
  # of 1.1.0 development reported Inf, and 1.0.0 reported 0 with p = 1.
  pairs <- list(c(1, 2), c(-0.41, 1.99), c(-3.2, -1.4), c(3, 7.62),
                c(0.1, 0.3), c(1.1, 0.2))
  for (ab in pairs) for (n in 3:5) {
    e <- expect_error(
      randomization_test(group_1 = rep(ab[1], n), group_2 = rep(ab[2], n),
                         statistic = "t"),
      class = "dmar_studentized_statistic_undefined")
    expect_s3_class(e, "dmar_refusal")
    expect_identical(e$reason, "constant_groups")
    expect_identical(e$alternative, "two_sided")
    expect_identical(c(e$n_1, e$n_2), c(n, n))
    expect_equal(e$p_value_mean, 2 / choose(2 * n, n), tolerance = 1e-14)
    expect_match(conditionMessage(e), "zero denominator", fixed = TRUE)
    expect_match(conditionMessage(e), "exchanges the two groups",
                 fixed = TRUE)

    # The difference in means is defined for the same data: the table, the
    # warning that the standardized mean difference is undefined, and the
    # smallest p-value the design allows.
    expect_warning(
      res <- withCallingHandlers(
        randomization_test(group_1 = rep(ab[1], n), group_2 = rep(ab[2], n)),
        warning = function(w) {
          # cliff_delta() says that the groups are completely separated;
          # that warning is its own and is not at issue here.
          if (!inherits(w, "dmar_smd_undefined")) invokeRestart("muffleWarning")
        }),
      class = "dmar_smd_undefined")
    direction <- sign(ab[1] - ab[2])
    expect_equal(rt_value(res, "statistic"), ab[1] - ab[2], tolerance = 1e-12)
    expect_equal(rt_value(res, "p_value"), e$p_value_mean, tolerance = 1e-12)
    for (term in c("smd", "smd_ci_lower", "smd_ci_upper",
                   "cles", "cles_ci_lower", "cles_ci_upper"))
      expect_true(is.na(rt_value(res, term)))
    expect_equal(rt_value(res, "cliff_delta"), direction)
    if (n == 4L) {
      # 2 / 70 is below 0.05, so the interval is the one shift at which
      # the two groups coincide.
      expect_equal(rt_value(res, "shift_ci_lower"), ab[1] - ab[2],
                   tolerance = 1e-6)
      expect_equal(rt_value(res, "shift_ci_upper"), ab[1] - ab[2],
                   tolerance = 1e-6)
    }
  }
})


test_that("randomization_test() states the exact p-value of the mean in its refusal, for each alternative and for unequal groups", {
  skip_on_cran()
  refusal <- function(g_1, g_2, alternative) {
    expect_error(randomization_test(group_1 = g_1, group_2 = g_2,
                                    statistic = "t",
                                    alternative = alternative),
                 class = "dmar_studentized_statistic_undefined")
  }
  mean_p <- function(g_1, g_2, alternative) {
    rt_value(suppressWarnings(
      randomization_test(group_1 = g_1, group_2 = g_2,
                         alternative = alternative, shift_ci = FALSE)),
      "p_value")
  }
  designs <- list(list(rep(5, 4), rep(3, 4)), list(rep(3, 4), rep(5, 4)),
                  list(rep(5, 3), rep(3, 5)), list(rep(1, 2), rep(4, 6)),
                  list(rep(0.3, 6), rep(0.1, 2)))
  for (d in designs) for (alternative in c("two_sided", "less", "greater")) {
    e <- refusal(d[[1]], d[[2]], alternative)
    expect_equal(e$p_value_mean, mean_p(d[[1]], d[[2]], alternative),
                 tolerance = 1e-14)
    expect_identical(c(e$n_1, e$n_2), c(length(d[[1]]), length(d[[2]])))
  }
  # 3 against 5: one reassignment in 56 in the direction of the data, and
  # no reassignment exchanges groups of unequal size.
  e <- refusal(rep(5, 3), rep(3, 5), "two_sided")
  expect_identical(e$p_value_mean, 1 / 56)
  expect_match(conditionMessage(e),
               "only the observed one leaves both groups constant",
               fixed = TRUE)
  expect_match(conditionMessage(e), "the exact p-value 1/56 = 0.01786",
               fixed = TRUE)
  e <- refusal(rep(5, 3), rep(3, 5), "less")
  expect_identical(e$p_value_mean, 1)
  expect_match(conditionMessage(e),
               "the exact p-value 1 (less), because every reassignment",
               fixed = TRUE)
  e <- refusal(rep(5, 4), rep(3, 4), "two_sided")
  expect_match(conditionMessage(e),
               "Every score in the first group is 5 and every score in the second is 3,",
               fixed = TRUE)
  expect_match(conditionMessage(e), "of the 70 reassignments of the 8 scores",
               fixed = TRUE)
  expect_match(conditionMessage(e),
               "the exact p-value 2/70 = 0.02857 (two_sided)", fixed = TRUE)
  # Through the formula interface as well.
  d <- data.frame(y = rep(c(5, 3), each = 4), g = rep(c("a", "b"), each = 4))
  expect_error(randomization_test(y ~ g, data = d, statistic = "t"),
               class = "dmar_studentized_statistic_undefined")
})


test_that("randomization_test() refuses the studentized statistic in the Monte Carlo branch without drawing", {
  set.seed(8128)
  before <- .Random.seed
  e <- expect_error(
    randomization_test(group_1 = rep(1, 30), group_2 = rep(3, 30),
                       statistic = "t"),
    class = "dmar_studentized_statistic_undefined")
  expect_identical(.Random.seed, before)
  expect_equal(e$p_value_mean, 2 / choose(60, 30), tolerance = 1e-14)
  expect_match(conditionMessage(e), "of the about 1.18e+17 reassignments",
               fixed = TRUE)
  expect_match(conditionMessage(e),
               paste0("With 60 scores the test is run by Monte Carlo ",
                      "sampling, which estimates this value and cannot ",
                      "report less than 1/(n_resamples + 1) = 0.00009999."),
               fixed = TRUE)
  e <- expect_error(
    randomization_test(group_1 = rep(1, 4), group_2 = rep(3, 4),
                       statistic = "t", exact = FALSE, seed = 5),
    class = "dmar_studentized_statistic_undefined")
  expect_identical(.Random.seed, before)
  expect_match(conditionMessage(e),
               "With exact = FALSE the test is run by Monte Carlo sampling",
               fixed = TRUE)
})


test_that("randomization_test() keeps the infinite ordering of reassignments inside the reference set", {
  skip_on_cran()
  # The observed groups have spread, so the observed statistic is finite;
  # the one reassignment that puts the four ones in the first group, and
  # the one that puts them in the second, have no spread and are ordered
  # as Inf and -Inf, where the difference in means also places them.
  g_1 <- c(1, 1, 1, 0)
  g_2 <- c(0, 0, 0, 1)
  res <- randomization_test(group_1 = g_1, group_2 = g_2, statistic = "t")
  expect_equal(rt_value(res, "statistic"), sqrt(2), tolerance = 1e-12)
  ref <- attr(res, "reference_distribution")
  expect_identical(sum(ref == Inf), 1L)
  expect_identical(sum(ref == -Inf), 1L)
  expect_identical(sum(attr(res, "reference_extreme")), 34L)
  expect_equal(rt_value(res, "p_value"), 34 / 70, tolerance = 1e-12)
  by_mean <- randomization_test(group_1 = g_1, group_2 = g_2)
  expect_equal(rt_value(by_mean, "p_value"), 34 / 70, tolerance = 1e-12)
  # A group that is constant only to within one part in 1e12 is not
  # constant: the statistic is very large and finite.
  near <- suppressWarnings(
    randomization_test(group_1 = c(5, 5, 5, 5 + 1e-12), group_2 = rep(3, 4),
                       statistic = "t"))
  expect_true(is.finite(rt_value(near, "statistic")))
  expect_gt(rt_value(near, "statistic"), 1e12)
})


test_that("randomization_test() gives the observed split and its mirror image the same studentized statistic", {
  # Two tight clusters. The split that separates them and its mirror image
  # have the same absolute statistic, so the two-sided p-value is 2 / 252.
  g_1 <- c(4.800000002, 4.800000000, 4.800000003, 4.800000001, 4.800000002)
  g_2 <- c(8.000000000, 8.000000001, 8.000000003, 8.000000003, 8.000000003)
  # The standardized mean difference of clusters this tight is beyond the
  # range of R's noncentral t, and cliff_delta() notes the separation;
  # neither warning concerns the test.
  res <- suppressWarnings(
    randomization_test(group_1 = g_1, group_2 = g_2, statistic = "t",
                       shift_ci = FALSE))
  expect_equal(rt_value(res, "p_value"), 2 / 252, tolerance = 1e-12)
})


test_that("randomization_test() returns p = 1 when every score is equal", {
  skip_on_cran()
  for (statistic in c("mean", "t")) {
    # Two warnings: the standardized mean difference is undefined, and
    # cliff_delta() says that its limits, both 0, are degenerate.
    expect_warning(
      expect_warning(
        res <- randomization_test(group_1 = rep(3, 4), group_2 = rep(3, 5),
                                  statistic = statistic),
        class = "dmar_smd_undefined"),
      class = "dmar_interval_degenerate")
    # The degenerate interval warning is passed on under this function's
    # name, with the class and the elements cliff_delta() gave it.
    w <- NULL
    withCallingHandlers(
      randomization_test(group_1 = rep(3, 4), group_2 = rep(3, 5),
                         statistic = statistic),
      dmar_interval_degenerate = function(cond) {
        w <<- cond
        invokeRestart("muffleWarning")
      },
      warning = function(cond) invokeRestart("muffleWarning"))
    expect_match(conditionMessage(w), "^randomization_test\\(\\): ")
    expect_identical(w$reason, "no_variation")
    expect_identical(w$n, c(4L, 5L))
    expect_identical(rt_value(res, "statistic"), 0)
    expect_identical(rt_value(res, "p_value"), 1)
    expect_true(is.na(rt_value(res, "smd")))
    expect_identical(rt_value(res, "cliff_delta"), 0)
  }
})


test_that("randomization_test() gives one p-value for both statistics in a balanced design", {
  skip_on_cran()
  # With n_1 = n_2 the separate-variances standard error is the pooled one
  # and the pooled t increases with the mean difference, whatever the two
  # variances are.
  set.seed(20261109)
  for (i in 1:50) {
    n <- sample(3:7, 1)
    g_1 <- rnorm(n, 0, 6)
    g_2 <- rnorm(n, 1, 1)
    alternative <- c("two_sided", "less", "greater")[i %% 3 + 1]
    # cliff_delta() warns when a sample happens to be completely
    # separated; the warning does not concern the test.
    p_mean <- rt_value(suppressWarnings(
      randomization_test(group_1 = g_1, group_2 = g_2,
                         alternative = alternative, shift_ci = FALSE)),
      "p_value")
    p_t <- rt_value(suppressWarnings(
      randomization_test(group_1 = g_1, group_2 = g_2, statistic = "t",
                         alternative = alternative, shift_ci = FALSE)),
      "p_value")
    expect_identical(p_t, p_mean)
  }
})


test_that("randomization_test() studentized p-values equal rational arithmetic on unbalanced designs", {
  skip_on_cran()
  oracle <- utils::read.csv(
    test_path("fixtures", "randomization_test_rational_oracle.csv"),
    stringsAsFactors = FALSE)
  scores <- function(s) as.numeric(strsplit(s, " ", fixed = TRUE)[[1]])
  # Two constant groups at different values have no studentized statistic,
  # and the call is refused. The oracle's value for those rows, the limit
  # that orders the observed assignment as most extreme, is the p-value of
  # the difference in means, which the refusal carries.
  p_dmar <- vapply(seq_len(nrow(oracle)), function(i) {
    res <- tryCatch(
      suppressWarnings(
        randomization_test(group_1 = scores(oracle$group_1[i]),
                           group_2 = scores(oracle$group_2[i]),
                           statistic = "t",
                           alternative = oracle$alternative[i],
                           shift_ci = FALSE)),
      dmar_studentized_statistic_undefined = function(e) e)
    if (inherits(res, "condition")) res$p_value_mean
    else rt_value(res, "p_value")
  }, numeric(1))
  refused <- vapply(seq_len(nrow(oracle)), function(i) {
    inherits(tryCatch(
      suppressWarnings(
        randomization_test(group_1 = scores(oracle$group_1[i]),
                           group_2 = scores(oracle$group_2[i]),
                           statistic = "t",
                           alternative = oracle$alternative[i],
                           shift_ci = FALSE)),
      dmar_studentized_statistic_undefined = function(e) e), "condition")
  }, logical(1))
  constant <- vapply(seq_len(nrow(oracle)), function(i) {
    g_1 <- scores(oracle$group_1[i])
    g_2 <- scores(oracle$group_2[i])
    length(unique(g_1)) == 1L && length(unique(g_2)) == 1L &&
      g_1[1] != g_2[1]
  }, logical(1))
  expect_identical(refused, constant)
  expect_gt(sum(refused), 0)

  # Binary, three-valued, Likert with an offset of 1000.3, three decimals
  # with and without an offset of a million, and constant groups.
  equal <- oracle$scenario != "two_point"
  expect_gt(sum(equal), 100)
  expect_equal(p_dmar[equal], oracle$p_exact[equal], tolerance = 1e-12)

  # Scores within 3e-9 of 4.8 or 8. Reassignments whose statistics differ
  # by less than the rounding bound are counted as tied, so the p-value can
  # be above the exact one. It is not below it. This is a stated limit of
  # the studentized comparison, not an equality.
  near <- !equal
  expect_true(all(p_dmar[near] >= oracle$p_exact[near] - 1e-12))
})


test_that("randomization_test() mean statistic equals whole-number arithmetic for many-digit scores", {
  skip_on_cran()
  set.seed(20261109)
  for (i in 1:30) {
    n_1 <- sample(3:6, 1)
    n_2 <- sample(3:7, 1)
    k_1 <- sample(0:9, n_1, replace = TRUE)
    k_2 <- sample(0:9, n_2, replace = TRUE)
    alternative <- c("two_sided", "less", "greater")[i %% 3 + 1]
    exact_p <- whole_number_randomization_p(k_1, k_2, alternative)

    # Thousandths on an offset of a million: the mean differences are tied
    # in exact arithmetic wherever the whole numbers are, and the rounding
    # of scores this size is far larger than the rounding of their
    # difference.
    res <- randomization_test(group_1 = (1e9 + k_1) / 1000,
                              group_2 = (1e9 + k_2) / 1000,
                              alternative = alternative, shift_ci = FALSE)
    expect_equal(rt_value(res, "p_value"), exact_p, tolerance = 1e-12)

    # One score about 1e9 times the spread of the rest.
    k_1[1] <- 3e9
    exact_p <- whole_number_randomization_p(k_1, k_2, alternative)
    res <- randomization_test(group_1 = k_1, group_2 = k_2,
                              alternative = alternative, shift_ci = FALSE)
    expect_equal(rt_value(res, "p_value"), exact_p, tolerance = 1e-12)
  }
})


test_that("randomization_test() shift interval is the set of shifts the test accepts", {
  skip_on_cran()
  set.seed(113)
  g_1 <- round(rnorm(6, 52, 8), 1)
  g_2 <- round(rnorm(7, 45, 8), 1)
  for (statistic in c("mean", "t")) {
    res <- randomization_test(group_1 = g_1, group_2 = g_2,
                              statistic = statistic)
    lo <- rt_value(res, "shift_ci_lower")
    hi <- rt_value(res, "shift_ci_upper")
    expect_true(is.finite(lo) && is.finite(hi))
    width <- hi - lo
    grid <- c(seq(lo - width, lo - width / 200, length.out = 12),
              seq(lo + width / 200, hi - width / 200, length.out = 25),
              seq(hi + width / 200, hi + width, length.out = 12))
    # Far outside the interval the shifted groups are completely
    # separated, which cliff_delta() warns about.
    p_grid <- vapply(grid, function(shift) {
      rt_value(suppressWarnings(
        randomization_test(group_1 = g_1 - shift, group_2 = g_2,
                           statistic = statistic, shift_ci = FALSE)),
        "p_value")
    }, numeric(1))
    inside <- grid > lo & grid < hi
    expect_true(all(p_grid[inside] > 0.05))
    expect_true(all(p_grid[!inside] <= 0.05))
  }
})


test_that("randomization_test() cannot reject any shift with three scores per group", {
  skip_on_cran()
  # The smallest two-sided p-value is 2 / 20, so the 95 percent interval is
  # the whole line for either statistic.
  for (statistic in c("mean", "t")) {
    res <- suppressWarnings(
      randomization_test(group_1 = c(80, 84, 79), group_2 = c(72, 75, 68),
                         statistic = statistic))
    expect_identical(rt_value(res, "shift_ci_lower"), -Inf)
    expect_identical(rt_value(res, "shift_ci_upper"), Inf)
  }
})


test_that("randomization_test() Monte Carlo draws and counts as DMAR 1.1.0 did for the mean statistic", {
  # The p-value recorded from the reviewed 1.1.0 build (df120a0d) for this
  # seeded call. It pins the order of the draws and shows that data with no
  # ties are counted as before.
  set.seed(113)
  g_1 <- rnorm(14, 0.7)
  g_2 <- rnorm(12)
  res <- randomization_test(group_1 = g_1, group_2 = g_2,
                            n_resamples = 2000, seed = 113)
  expect_equal(rt_value(res, "p_value"), 0.90604697651174415,
               tolerance = 1e-14)
  expect_equal(rt_value(res, "exact"), 0)
})


test_that("randomization_test() Monte Carlo handles infinite studentized statistics", {
  skip_on_cran()
  g_1 <- c(1, 1, 1, 1, 0)
  g_2 <- c(0, 0, 0, 0, 1)
  res <- randomization_test(group_1 = g_1, group_2 = g_2, statistic = "t",
                            exact = FALSE, n_resamples = 20000, seed = 113,
                            shift_ci = FALSE)
  expect_true(any(is.infinite(attr(res, "reference_distribution"))))
  exact_p <- 52 / 252
  mc_se <- sqrt(exact_p * (1 - exact_p) / 20000)
  expect_lt(abs(rt_value(res, "p_value") - exact_p), 3 * mc_se)
})


test_that("randomization_test() requires n_resamples to be a whole number of at least 1", {
  g_1 <- c(80, 84, 79, 88, 83)
  g_2 <- c(72, 75, 68, 81, 74)
  for (bad in list(2.5, 0, -1, NA_real_, Inf, "10", c(10, 20))) {
    expect_error(randomization_test(group_1 = g_1, group_2 = g_2,
                                    exact = FALSE, n_resamples = bad),
                 "'n_resamples'", class = "dmar_count_invalid")
  }
  a <- randomization_test(group_1 = g_1, group_2 = g_2, exact = FALSE,
                          n_resamples = 200, seed = 113)
  b <- randomization_test(group_1 = g_1, group_2 = g_2, exact = FALSE,
                          n_resamples = 200L, seed = 113)
  expect_identical(a$value, b$value)
  expect_identical(rt_value(a, "n_evaluated"), 200)
})


test_that("randomization_test() refuses a response that is not finite", {
  expect_error(randomization_test(group_1 = c(1, 2, Inf), group_2 = c(2, 3, 4)),
               "must be finite")
  expect_error(randomization_test(c(1, 2, 3, -Inf, 5, 6),
                                  rep(c("a", "b"), each = 3)),
               "must be finite")
})


test_that("plot_randomization_test() shades the reassignments the p-value counted", {
  skip_on_cran()
  skip_if_not_installed("ggplot2")
  shaded_share <- function(res) {
    mean(plot_randomization_test(res)$data$region == "At least as extreme")
  }
  treatment <- c(80, 84, 79, 88, 83)
  control   <- c(72, 75, 68, 81, 74)
  designs <- list(
    list(treatment, control),
    # Scores so small that an absolute tolerance would shade everything.
    list(treatment * 1e-13, control * 1e-13),
    # Tied mean differences on an offset, which an absolute tolerance
    # smaller than their rounding would leave unshaded.
    list(1e5 + c(0.1, 0.4, 0.7, 0.2), 1e5 + c(0.3, 0.5, 0.6, 0.9))
  )
  for (d in designs) {
    res <- randomization_test(group_1 = d[[1]], group_2 = d[[2]])
    expect_equal(shaded_share(res), rt_value(res, "p_value"),
                 tolerance = 1e-12)
  }
  # The third design in whole numbers gives the exact p-value.
  res <- randomization_test(group_1 = designs[[3]][[1]],
                            group_2 = designs[[3]][[2]])
  expect_equal(rt_value(res, "p_value"),
               whole_number_randomization_p(c(1, 4, 7, 2), c(3, 5, 6, 9),
                                            "two_sided"),
               tolerance = 1e-12)

  # A result saved before the counted reassignments were stored is shaded
  # by the counting rule.
  res <- randomization_test(group_1 = treatment, group_2 = control)
  attr(res, "reference_extreme") <- NULL
  expect_equal(shaded_share(res), rt_value(res, "p_value"), tolerance = 1e-12)
})


test_that("plot_randomization_test() draws a result with infinite reference values", {
  skip_on_cran()
  skip_if_not_installed("ggplot2")
  res <- randomization_test(group_1 = c(1, 1, 1, 1, 0),
                            group_2 = c(0, 0, 0, 0, 1), statistic = "t")
  n_infinite <- sum(is.infinite(attr(res, "reference_distribution")))
  expect_identical(n_infinite, 2L)
  p <- plot_randomization_test(res)
  expect_equal(nrow(p$data), 252 - n_infinite)
  expect_match(p$labels$subtitle,
               "2 reassignments with an infinite statistic are not drawn")
  expect_no_warning(ggplot2::ggplot_build(p))

  # The observed statistic is finite, so its line is drawn.
  layers <- vapply(p$layers, function(l) class(l$geom)[1], character(1))
  expect_true("GeomVline" %in% layers)
})


test_that("randomization_test() rejects a shift whose p-value equals alpha", {
  # 1 - 0.90 is 0.09999999999999998 in binary arithmetic, below the p-value
  # 2 / 20 that the 0.10 level rejects, and the 90 percent interval of three
  # scores per group was the whole line. The acceptance set, from every
  # breakpoint of the p-value as a function of the shift, is (-5.7, 5.3).
  g1 <- c(-1.5, 3.6, 0)
  g2 <- c(4.2, 3.0, -1.7)
  res <- suppressWarnings(
    randomization_test(group_1 = g1, group_2 = g2, conf_level = 0.90))
  expect_equal(rt_value(res, "shift_ci_lower"), -5.7, tolerance = 1e-5)
  expect_equal(rt_value(res, "shift_ci_upper"), 5.3, tolerance = 1e-5)
  # The shift just outside the lower limit has p = 0.1 and is rejected.
  out <- suppressWarnings(
    randomization_test(group_1 = g1 + 5.71, group_2 = g2, shift_ci = FALSE))
  expect_equal(rt_value(out, "p_value"), 0.1)
  # At 0.95 no shift can be rejected, as before.
  res <- suppressWarnings(
    randomization_test(group_1 = g1, group_2 = g2, conf_level = 0.95))
  expect_identical(rt_value(res, "shift_ci_lower"), -Inf)
  expect_identical(rt_value(res, "shift_ci_upper"), Inf)
})


test_that("randomization_test() studentized statistic does not depend on the unit of the scores", {
  # Squared residuals of scores near 1e-200 underflow and of scores near
  # 1e200 overflow; the statistic was Inf (or 0) with p = 1.
  g1 <- c(3, 5, 9, 4)
  g2 <- c(1, 2, 2, 6, 3)
  ref <- randomization_test(group_1 = g1, group_2 = g2, statistic = "t",
                            shift_ci = FALSE)
  for (unit in c(1e-200, 1e200)) {
    res <- suppressWarnings(
      randomization_test(group_1 = g1 * unit, group_2 = g2 * unit,
                         statistic = "t", shift_ci = FALSE))
    expect_equal(rt_value(res, "statistic"), rt_value(ref, "statistic"),
                 tolerance = 1e-12)
    expect_equal(rt_value(res, "p_value"), rt_value(ref, "p_value"))
  }
})


test_that("randomization_test() studentized p-value is unchanged by an offset the bound can separate", {
  skip_on_cran()
  # Whole-number scores on an offset of 1e9 are exact in binary arithmetic.
  # The rounding bound of the studentized statistic there is of the order of
  # 1e-6, well below the gaps between distinct statistics, so each p-value
  # is that of the scores without the offset. A bound 10,000 times wider
  # merges statistics that differ.
  set.seed(7)
  for (i in 1:20) {
    n_1 <- sample(4:6, 1)
    n_2 <- sample(4:6, 1)
    g1 <- sample(1:30, n_1, replace = TRUE)
    g2 <- sample(1:30, n_2, replace = TRUE) + 2
    plain <- randomization_test(group_1 = g1, group_2 = g2, statistic = "t",
                                shift_ci = FALSE)
    offset <- randomization_test(group_1 = g1 + 1e9, group_2 = g2 + 1e9,
                                 statistic = "t", shift_ci = FALSE)
    expect_equal(rt_value(offset, "p_value"), rt_value(plain, "p_value"),
                 tolerance = 1e-12)
  }
})


test_that("randomization_test() Monte Carlo shift interval uses the add-one p-value", {
  # With 50 draws and alpha = .05 the add-one p-value (1 + r) / 51 exceeds
  # alpha from r = 2 and r / 50 only from r = 3, so the two rules give
  # different limits. The same seed draws the same reassignments, so the
  # p-value at a shift is the p-value of the call with that shift removed
  # from the first group.
  set.seed(113)
  g1 <- round(rnorm(12, 1, 2), 2)
  g2 <- round(rnorm(12, 0, 2), 2)
  res <- randomization_test(group_1 = g1, group_2 = g2, exact = FALSE,
                            n_resamples = 50, seed = 113)
  p_at <- function(shift) {
    rt_value(randomization_test(group_1 = g1 - shift, group_2 = g2,
                                exact = FALSE, n_resamples = 50, seed = 113,
                                shift_ci = FALSE), "p_value")
  }
  lower <- rt_value(res, "shift_ci_lower")
  upper <- rt_value(res, "shift_ci_upper")
  expect_gt(p_at(lower + 1e-3), 0.05)
  expect_lte(p_at(lower - 1e-3), 0.05)
  expect_gt(p_at(upper - 1e-3), 0.05)
  expect_lte(p_at(upper + 1e-3), 0.05)
})


test_that("plot_randomization_test() shades by the stored flags, not by its fallback rule", {
  skip_if_not_installed("ggplot2")
  # One score of 3e9 among single digits: the p-value counts 35 of the 126
  # reassignments, and the fallback rule for a result saved without its
  # flags, whose tolerance is relative to the largest statistic, counts 56.
  res <- randomization_test(group_1 = c(3e9, 4, 7, 2),
                            group_2 = c(3, 5, 6, 9, 1))
  expect_equal(rt_value(res, "p_value"), 35 / 126)
  p <- plot_randomization_test(res)
  expect_identical(sum(p$data$region == "At least as extreme"), 35L)
})


test_that("plot_randomization_test() says when every counted reassignment is infinite", {
  skip_on_cran()
  skip_if_not_installed("ggplot2")
  # Under exact enumeration the observed assignment is among those counted
  # and its statistic is finite, so this arises only in a Monte Carlo
  # sample: one that holds the reassignment with every one in the first
  # group (statistic Inf) and none of the sixteen tied with the observed
  # one. The first seed that draws such a sample of 6 is used.
  g_1 <- c(1, 1, 1, 0)
  g_2 <- c(0, 0, 0, 1)
  res <- NULL
  for (seed in 1:5000) {
    res <- randomization_test(group_1 = g_1, group_2 = g_2, statistic = "t",
                              alternative = "greater", exact = FALSE,
                              n_resamples = 6, seed = seed, shift_ci = FALSE)
    counted <- attr(res, "reference_extreme")
    ref <- attr(res, "reference_distribution")
    if (any(counted) && all(is.infinite(ref[counted]))) break
    res <- NULL
  }
  expect_false(is.null(res))
  n_infinite <- sum(is.infinite(attr(res, "reference_distribution")))
  p <- plot_randomization_test(res)
  lines <- strsplit(p$labels$subtitle, "\n", fixed = TRUE)[[1]]
  expect_length(lines, 2L)
  expect_false(grepl("shaded area", lines[1], fixed = TRUE))
  expect_match(lines[1], "; p = ", fixed = TRUE)
  expect_match(lines[2], "with an infinite statistic", fixed = TRUE)
  expect_match(lines[2], "not drawn", fixed = TRUE)
  expect_match(lines[2], "no area is shaded", fixed = TRUE)
  expect_equal(nrow(p$data), 6 - n_infinite)
  expect_false(any(p$data$region == "At least as extreme"))
  # With finite counted reassignments as well, the first line keeps
  # "shaded area is p".
  res <- randomization_test(group_1 = c(1, 1, 1, 1, 0),
                            group_2 = c(0, 0, 0, 0, 1), statistic = "t")
  lines <- strsplit(plot_randomization_test(res)$labels$subtitle, "\n",
                    fixed = TRUE)[[1]]
  expect_match(lines[1], "shaded area is p", fixed = TRUE)
  expect_length(lines, 2L)
})
