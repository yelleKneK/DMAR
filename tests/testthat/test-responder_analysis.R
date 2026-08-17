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
  expect_equal(v("lower_limit"), ctr - hw)
  expect_equal(v("upper_limit"), ctr + hw)
  # Boundary count: still inside [0, 1], unlike Wald.
  at0 <- ci_proportion(0, 20)
  expect_gte(at0$value[at0$term == "lower_limit"], 0)
  expect_gt(at0$value[at0$term == "upper_limit"], 0)
  # Wald option matches its textbook form.
  w <- ci_proportion(17, 50, method = "wald")
  expect_equal(w$value[w$term == "lower_limit"], p - z * sqrt(p * (1 - p) / n))
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
  expect_equal(res$lower_limit[1], cp$value[cp$term == "lower_limit"])
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
  expect_true(is.na(nr$lower_limit[nr$group == "nnt"]))
  expect_error(responder_analysis(d$x, d$g[-1], 10), "one label per")
})
