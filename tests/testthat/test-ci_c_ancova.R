test_that("ci_c_ancova returns the per-comparison t interval by default", {
  # Maxwell, Delaney, and Kelley (2027) example: three groups, ten per group.
  # SSwithin_x = 752.5 is the within-groups sum of squares of the pretest
  # covariate, recomputed exactly from the chapter data (shipped as
  # depression_bdi): sum by group of (bdi_pre - group mean)^2. An earlier
  # version of this example carried 313.37, which is a different row of the
  # same decomposition, the within-groups SS the covariate explains in the
  # posttest (SS_within(Post) - SS_error(ANCOVA) = 313.365).
  out <- ci_c_ancova(adj_means = c(7.5, 12, 14), s_ancova = sqrt(29),
                     c_weights = c(1, -1, 0), n_per_group = 10,
                     cov_means = c(17, 17.7, 17.4), SSwithin_x = 752.5)
  expect_s3_class(out, "data.frame")
  expect_equal(out$term, c("ci_lower", "contrast", "ci_upper"))
  # Point estimate is the contrast of adjusted means.
  expect_equal(out$value[out$term == "contrast"], 7.5 - 12)
  # Interval is symmetric about the point estimate (t-based, single contrast).
  psi <- out$value[out$term == "contrast"]
  expect_equal(out$value[out$term == "ci_lower"] - psi,
               -(out$value[out$term == "ci_upper"] - psi))
  # Anchor the limits so the constant cannot drift again.
  expect_equal(out$value[out$term == "ci_lower"], -9.458423, tolerance = 1e-6)
  expect_equal(out$value[out$term == "ci_upper"],  0.458423, tolerance = 1e-6)
})

test_that("the chapter example's SSwithin_x reproduces from the chapter data", {
  # The authoritative source for the 752.5 in the example above: the
  # book's own data, shipped as depression_bdi.
  data("depression_bdi", package = "DMAR", envir = environment())
  d <- depression_bdi
  ss_within_x <- sum(tapply(d$bdi_pre, d$condition,
                            function(x) sum((x - mean(x))^2)))
  expect_equal(ss_within_x, 752.5)
  # And the fingerprint quantities the example quotes.
  fit <- lm(bdi_post ~ bdi_pre + condition, data = d)
  expect_equal(sum(residuals(fit)^2) / fit$df.residual, 29, tolerance = 0.005)
  expect_equal(as.vector(tapply(d$bdi_pre, d$condition, mean)),
               c(17, 17.7, 17.4))
})

test_that("procedure = 'bryant_paulson' forwards to ci_c_ancova_bp", {
  skip_on_cran()  # two Bryant-Paulson root finds; the t interval anchor above runs on CRAN
  adj <- c(7.5, 12, 14)
  s <- sqrt(29)

  via_ci_c_ancova <- ci_c_ancova(adj_means = adj, s_ancova = s, n_per_group = 10,
                                 procedure = "bryant_paulson")
  direct <- ci_c_ancova_bp(adj_means = adj, s_ancova = s, n_per_group = 10)

  expect_equal(via_ci_c_ancova, direct)
  expect_identical(attr(via_ci_c_ancova, "critical_value"),
                   attr(direct, "critical_value"))
})

test_that("procedure = 'bryant_paulson' forwards c_weights and extra arguments", {
  skip_on_cran()  # two Bryant-Paulson root finds; the t interval anchor above runs on CRAN
  adj <- c(3.595, 3.619, 4.102, 4.515, 4.618, 4.876)
  s <- sqrt(0.01326)
  cw <- c(0.5, 0.5, -0.25, -0.25, -0.25, -0.25)

  via_ci_c_ancova <- ci_c_ancova(adj_means = adj, s_ancova = s, c_weights = cw,
                                 n_per_group = 4, procedure = "bryant_paulson",
                                 df = 14, contrast_type = "allowance")
  direct <- ci_c_ancova_bp(adj_means = adj, s_ancova = s, c_weights = cw,
                           n_per_group = 4, df = 14, contrast_type = "allowance")

  expect_equal(via_ci_c_ancova, direct)
})

test_that("procedure = 'bryant_paulson' requires adj_means", {
  expect_error(
    ci_c_ancova(contrast = -4.5, s_ancova = sqrt(29), n_per_group = 10,
                procedure = "bryant_paulson"),
    "requires 'adj_means'"
  )
})

test_that("an unknown procedure is rejected", {
  expect_error(
    ci_c_ancova(adj_means = c(7.5, 12, 14), s_ancova = sqrt(29), n_per_group = 10,
                procedure = "tukey"),
    "should be one of"
  )
})

test_that("ci_c_ancova() requires a positive s_ancova on the contrast and adj_means paths", {
  # Through 1.0.1 s_ancova was never checked: omitted, each limit added psi
  # to a product with NULL, which data.frame() recycled, so the table gave a
  # zero-width interval at the estimate (-4.5, -4.5, -4.5 here) labeled
  # "Confidence level: 95%" and the sentence "psi = -4.50, 95% CI [-4.50,
  # -4.50]"; an s_ancova of 0 gave the same, and a negative one swapped the
  # limits.
  args <- list(adj_means = c(7.5, 12, 14), c_weights = c(1, -1, 0),
               n_per_group = c(10, 10, 10), cov_means = c(5, 6, 7), SSwithin_x = 750)
  expect_error(do.call(ci_c_ancova, args),
               "standard deviation of the errors of the ANCOVA model")
  expect_error(ci_c_ancova(contrast = -4.5, c_weights = c(1, -1, 0), n_per_group = 10,
                           cov_means = c(5, 6, 7), SSwithin_x = 750),
               "standard deviation of the errors of the ANCOVA model")
  for (bad in list(0, -1, NA_real_, Inf, c(1, 2), "1")) {
    expect_error(do.call(ci_c_ancova, c(args, list(s_ancova = bad))),
                 "'s_ancova' must be a single positive, finite number",
                 info = format(bad))
  }
  # The Bryant-Paulson procedure refuses the same values in the same words.
  expect_error(ci_c_ancova(adj_means = c(7.5, 12, 14), n_per_group = 10,
                           procedure = "bryant_paulson"),
               "standard deviation of the errors of the ANCOVA model")
  expect_error(ci_c_ancova(adj_means = c(7.5, 12, 14), s_ancova = NA_real_,
                           n_per_group = 10, procedure = "bryant_paulson"),
               "'s_ancova' must be a single positive, finite number")
  # A valid s_ancova still gives the interval it did.
  ok <- do.call(ci_c_ancova, c(args, list(s_ancova = 2)))
  se <- 2 * sqrt(2 / 10 + (5 - 6)^2 / 750)
  expect_equal(ok$value, -4.5 + c(-1, 0, 1) * stats::qt(.975, 26) * se,
               tolerance = 1e-12)
})
