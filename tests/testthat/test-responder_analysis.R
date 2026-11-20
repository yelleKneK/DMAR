.resp_data <- function() {
  set.seed(113)
  list(x = c(rnorm(60, 8, 9), rnorm(60, 13, 9)),
       g = rep(c("control", "treatment"), each = 60))
}

test_that("ci_proportion() Wilson matches the closed form and stays in [0, 1]", {
  res <- ci_proportion(17, 50)
  p <- 17 / 50; z <- qnorm(.975); n <- 50
  den <- 1 + z^2 / n
  ctr <- (p + z^2 / (2 * n)) / den
  hw  <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / den
  v <- function(t) res$value[res$term == t]
  expect_equal(v("ci_lower"), ctr - hw)
  expect_equal(v("ci_upper"), ctr + hw)
  # Boundary count: still inside [0, 1], unlike Wald.
  at0 <- ci_proportion(0, 20)
  expect_gte(at0$value[at0$term == "ci_lower"], 0)
  expect_gt(at0$value[at0$term == "ci_upper"], 0)
  # Wald option matches its textbook form.
  w <- ci_proportion(17, 50, method = "wald")
  expect_equal(w$value[w$term == "ci_lower"], p - z * sqrt(p * (1 - p) / n))
  expect_error(ci_proportion(5, 4), "at least 'successes'")
})

test_that("responder_analysis() counts, intervals, difference, NNT, omnibus", {
  d <- .resp_data()
  res <- responder_analysis(d$x, d$g, threshold = 10)
  expect_s3_class(res, "dmar_tbl")
  expect_identical(res$group,
                   c("control", "treatment", "difference", "nnt", "omnibus"))
  # Counts agree with direct tabulation.
  expect_equal(res$responders[1], sum(d$x[d$g == "control"] >= 10))
  expect_equal(res$responders[2], sum(d$x[d$g == "treatment"] >= 10))
  expect_equal(res$estimate[1], res$responders[1] / 60)
  # Per-group rows reproduce ci_proportion exactly.
  cp <- ci_proportion(res$responders[1], 60)
  expect_equal(res$ci_lower[1], cp$value[cp$term == "ci_lower"])
  # Difference = treatment - control; NNT = 1/|difference|.
  expect_equal(res$estimate[3], res$estimate[2] - res$estimate[1])
  expect_equal(res$estimate[4], 1 / abs(res$estimate[3]))
  # Omnibus equals the uncorrected chi-square test.
  tab <- table(d$g, d$x >= 10)
  om  <- chisq.test(tab, correct = FALSE)
  expect_equal(res$chi_square[5], unname(om$statistic))
  expect_equal(res$p_value[5], om$p.value)
})

test_that("responder_analysis() direction, sweep, and NNT interval logic", {
  d <- .resp_data()
  # 'le' flips the definition of response.
  le <- responder_analysis(d$x, d$g, threshold = 10, direction = "le")
  ge <- responder_analysis(d$x, d$g, threshold = 10)
  expect_equal(le$responders[1] + ge$responders[1], 60)
  # Sweep stacks per-threshold tables with a leading threshold column.
  sw <- responder_analysis(d$x, d$g, threshold = 10, sweep = c(5, 15))
  expect_identical(names(sw)[1], "threshold")
  expect_equal(nrow(sw), 3 * 5)
  # When the difference interval includes zero, the NNT interval is NA.
  set.seed(113)
  null_x <- rnorm(80); null_g <- rep(c("a", "b"), 40)
  nr <- responder_analysis(null_x, null_g, threshold = 0)
  expect_true(is.na(nr$ci_lower[nr$group == "nnt"]))
  expect_error(responder_analysis(d$x, d$g[-1], 10), "one label per")
})

test_that("responder_analysis() carries the confidence set for the NNT and prints it when it is two intervals", {
  d <- .resp_data()
  res <- responder_analysis(d$x, d$g, threshold = 10)
  expect_s3_class(res, "dmar_responder")
  expect_s3_class(res, "dmar_tbl")
  diff_row <- res[res$group == "difference", ]
  # The interval for the difference includes zero: no limits on the nnt
  # row, and the set is two intervals whose finite ends are the
  # reciprocals of the limits of the difference.
  expect_lt(diff_row$ci_lower, 0)
  expect_gt(diff_row$ci_upper, 0)
  expect_identical(attr(res, "nnt_set"), "two_intervals")
  set <- attr(res, "nnt_confidence_set")
  expect_identical(names(set), c("lower", "upper"))
  expect_equal(set$lower, c(-Inf, 1 / diff_row$ci_upper))
  expect_equal(set$upper, c(1 / diff_row$ci_lower, Inf))
  expect_identical(attr(res, "group_labels"), c("control", "treatment"))
  expect_output(print(res),
                paste0("NNT <= -10.806 (control responds more often) or ",
                       "NNT >= 3.966 (treatment responds more often), ",
                       "joined at infinity (Altman, 1998)"),
                fixed = TRUE)
  # The numbers of the worked example on the help page.
  expect_identical(res$responders[1:2], c(30, 35))
  expect_equal(round(res$estimate[1:4], 3), c(0.500, 0.583, 0.083, 12))
  expect_equal(round(res$ci_lower[1:3], 3), c(0.377, 0.457, -0.093))
  expect_equal(round(res$ci_upper[1:3], 3), c(0.623, 0.699, 0.252))
  expect_equal(round(c(res$chi_square[5], res$p_value[5]), 3), c(0.839, 0.360))
  expect_equal(round(c(set$upper[1], set$lower[2]), 2), c(-10.81, 3.97))

  # An interval for the difference that excludes zero: one interval, on the
  # signed scale, whose absolute ends are the limits of the nnt row; no
  # line is printed.
  one <- responder_analysis(d$x, d$g, threshold = 5)
  expect_identical(attr(one, "nnt_set"), "interval")
  set_one <- attr(one, "nnt_confidence_set")
  expect_identical(nrow(set_one), 1L)
  nnt_row <- one[one$group == "nnt", ]
  expect_equal(c(set_one$lower, set_one$upper),
               c(nnt_row$ci_lower, nnt_row$ci_upper))
  expect_equal(round(c(nnt_row$estimate, nnt_row$ci_lower,
                       nnt_row$ci_upper), 2), c(5.45, 3.03, 37.98))
  expect_false(any(grepl("two intervals", capture.output(print(one)))))
  # With the second group responding less often the interval is negative.
  flipped <- responder_analysis(d$x, factor(d$g, levels = c("treatment",
                                                            "control")),
                                threshold = 5)
  set_flip <- attr(flipped, "nnt_confidence_set")
  expect_equal(c(set_flip$lower, set_flip$upper),
               -c(set_one$upper, set_one$lower))
  expect_equal(flipped$ci_lower[flipped$group == "nnt"],
               nnt_row$ci_lower)

  # A sweep: one value of nnt_set per threshold, and a threshold column.
  sw <- responder_analysis(d$x, d$g, threshold = 10, sweep = c(5, 15))
  expect_identical(attr(sw, "nnt_set"),
                   c("two_intervals", "interval", "two_intervals"))
  set_sw <- attr(sw, "nnt_confidence_set")
  expect_identical(names(set_sw), c("threshold", "lower", "upper"))
  expect_identical(set_sw$threshold, c(10, 10, 5, 15, 15))
  expect_equal(set_sw[1:2, c("lower", "upper")], set, ignore_attr = TRUE)
  out <- capture.output(print(sw))
  expect_true(any(grepl("two intervals at 2 of the 3 thresholds", out,
                        fixed = TRUE)))
  expect_true(any(grepl("at threshold 15, NNT <= -9.979", out, fixed = TRUE)))

  # More than two groups: no difference, no NNT, no set.
  three <- responder_analysis(c(1, 2, 3, 4, 5, 6, 7, 8, 9),
                              rep(c("a", "b", "c"), 3), 4)
  expect_identical(three$group, c("a", "b", "c", "omnibus"))
  expect_null(attr(three, "nnt_set"))
  expect_null(attr(three, "nnt_confidence_set"))
})

test_that("responder_analysis() returns what its help page says at each edge", {
  g <- rep(c("a", "b"), each = 10)
  # Equal proportions: difference 0, NNT infinite, no limits, two intervals.
  eq <- responder_analysis(rep(c(0, 1, 0, 1), each = 5), g, 0.5)
  expect_identical(eq$estimate[3:4], c(0, Inf))
  expect_true(all(is.na(c(eq$ci_lower[4], eq$ci_upper[4]))))
  expect_identical(attr(eq, "nnt_set"), "two_intervals")
  expect_identical(eq$chi_square[5], 0)
  # None in one group and all in the other.
  sep <- responder_analysis(rep(c(0, 1), each = 10), g, 0.5)
  expect_identical(sep$estimate[3:4], c(1, 1))
  expect_equal(sep$ci_upper[3], 1)
  expect_equal(sep$ci_lower[1], 0, tolerance = 1e-15)
  expect_equal(sep$ci_upper[2], 1, tolerance = 1e-15)
  expect_equal(sep$ci_lower[4], 1)
  expect_identical(attr(sep, "nnt_set"), "interval")
  # None in either group: the chi square statistic is 0 / 0.
  none <- responder_analysis(rep(0, 20), g, 0.5)
  expect_identical(none$responders[1:2], c(0, 0))
  expect_identical(none$estimate[3], 0)
  expect_true(is.nan(none$chi_square[5]))
  expect_true(is.nan(none$p_value[5]))
  # A group of one.
  single <- responder_analysis(c(1, 0, 0, 1, 1), c("a", "b", "b", "b", "b"),
                               0.5)
  expect_identical(single$n[1:2], c(1, 4))
  expect_equal(round(c(single$ci_lower[1], single$ci_upper[1]), 3),
               c(0.207, 1))
  # An outcome at the threshold is a responder under either direction.
  x <- c(10, 10, 9, 11, 10, 12)
  ge <- responder_analysis(x, rep(c("a", "b"), 3), 10)
  le <- responder_analysis(x, rep(c("a", "b"), 3), 10, direction = "le")
  expect_identical(sum(ge$responders[1:2]) + sum(le$responders[1:2]),
                   6 + sum(x == 10))
  # Missing values are removed before anything is counted.
  with_na <- responder_analysis(c(x, NA, 12), c(rep(c("a", "b"), 3), "a", NA),
                                10)
  expect_identical(with_na$n[1:2], ge$n[1:2])
  expect_error(responder_analysis(c(1, 2, 3, 4), rep("a", 4), 2),
               "at least two groups")
  expect_error(responder_analysis(c(1, 2, 3), c("a", "a", "b"), 2),
               "at least 4 observations")
})
