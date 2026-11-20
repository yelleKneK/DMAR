## ci_c() / ci_sc() / ci_c_ancova() / ci_sc_ancova() -- contrast CIs.

# Common balanced three-group fixture.
.fix <- function() {
  list(
    means     = c(10, 12, 14),
    s_anova   = 2.5,
    c_weights = c(-1, 0, 1),
    n         = c(20, 20, 20),
    N         = 60
  )
}

test_that("ci_c() returns a tidy data.frame bracketing the unstandardized contrast", {
  f <- .fix()
  res <- ci_c(means = f$means, s_anova = f$s_anova, c_weights = f$c_weights,
              n_per_group = f$n, N = f$N, conf_level = 0.95)
  expect_s3_class(res, "data.frame")
  expect_named(res, c("term", "value"))
  expect_setequal(res$term, c("ci_lower", "contrast", "ci_upper"))
  expect_lt(res$value[res$term == "ci_lower"], res$value[res$term == "contrast"])
  expect_gt(res$value[res$term == "ci_upper"], res$value[res$term == "contrast"])
})

test_that("ci_sc() returns a tidy data.frame with a standardized-contrast CI", {
  f <- .fix()
  res <- ci_sc(means = f$means, s_anova = f$s_anova, c_weights = f$c_weights,
               n_per_group = f$n, N = f$N, conf_level = 0.95)
  expect_s3_class(res, "data.frame")
  expect_setequal(res$term, c("ci_lower", "std_contrast", "ci_upper"))
  expect_lt(res$value[res$term == "ci_lower"], res$value[res$term == "std_contrast"])
})

test_that("ci_c_ancova() returns a tidy ANCOVA-contrast CI bracketing psi", {
  res <- ci_c_ancova(adj_means = c(10, 12, 14), s_ancova = 2.3,
                     c_weights = c(-1, 0, 1),
                     n_per_group = c(20, 20, 20),
                     cov_means = c(5, 5.2, 4.9),
                     SSwithin_x = 50,
                     conf_level = 0.95)
  expect_s3_class(res, "data.frame")
  expect_named(res, c("term", "value"))
  expect_setequal(res$term, c("ci_lower", "contrast", "ci_upper"))
  expect_lt(res$value[res$term == "ci_lower"], res$value[res$term == "contrast"])
})

test_that("ci_sc_ancova() returns a tidy standardized ANCOVA-contrast CI", {
  res <- ci_sc_ancova(adj_means = c(10, 12, 14), s_anova = 2.5, s_ancova = 2.3,
                      c_weights = c(-1, 0, 1),
                      n_per_group = c(20, 20, 20),
                      cov_means = c(5, 5.2, 4.9),
                      SSwithin_x = 50,
                      conf_level = 0.95)
  expect_s3_class(res, "data.frame")
  expect_named(res, c("term", "value"))
  expect_setequal(res$term, c("ci_lower", "std_contrast", "ci_upper"))
  expect_type(res$value, "double")
  expect_lt(res$value[res$term == "ci_lower"], res$value[res$term == "std_contrast"])
  expect_identical(attr(res, "standardizer"), "s_ancova")
})

test_that("ci_sc_ancova(standardizer = 's_anova') records the alternative divisor", {
  res <- ci_sc_ancova(adj_means = c(10, 12, 14), s_anova = 2.5, s_ancova = 2.3,
                      c_weights = c(-1, 0, 1),
                      n_per_group = c(20, 20, 20),
                      cov_means = c(5, 5.2, 4.9),
                      SSwithin_x = 50,
                      standardizer = "s_anova",
                      conf_level = 0.95)
  expect_identical(attr(res, "standardizer"), "s_anova")
})

test_that("ci_sc_ancova() requires a positive s_ancova, and s_anova when it standardizes", {
  # Through 1.0.1 s_ancova was checked only when s_anova was the
  # standardizer, and then only for being omitted: the default call without
  # it stopped with "argument is of length zero", and neither standard
  # deviation was refused at 0 or below (an s_ancova of 0 gave an infinite
  # estimate).
  args <- list(adj_means = c(7.5, 12, 14), c_weights = c(0.5, -0.5, 0),
               n_per_group = c(10, 10, 10), cov_means = c(5, 6, 7), SSwithin_x = 750)
  expect_error(do.call(ci_sc_ancova, args),
               "standard deviation of the errors of the ANCOVA model")
  expect_error(do.call(ci_sc_ancova, c(args, list(standardizer = "s_anova",
                                                  s_anova = 3))),
               "standard deviation of the errors of the ANCOVA model")
  expect_error(do.call(ci_sc_ancova, c(args, list(standardizer = "s_anova",
                                                  s_ancova = 2))),
               "standard deviation of the errors of the ANOVA model")
  for (bad in list(0, -1, NA_real_, Inf, c(1, 2), "1")) {
    expect_error(do.call(ci_sc_ancova, c(args, list(s_ancova = bad))),
                 "'s_ancova' must be a single positive, finite number",
                 info = format(bad))
    expect_error(do.call(ci_sc_ancova, c(args, list(s_ancova = 2,
                                                    s_anova = bad,
                                                    standardizer = "s_anova"))),
                 "'s_anova' must be a single positive, finite number",
                 info = format(bad))
  }
  expect_error(do.call(ci_sc_ancova, c(args, list(s_ancova = 2,
                                                  standardizer = "sd"))),
               "either 's_anova' or 's_ancova'")
  # s_anova is not needed when s_ancova standardizes, and a valid call
  # still gives the interval it did.
  ok <- do.call(ci_sc_ancova, c(args, list(s_ancova = 2)))
  expect_equal(ok$value[ok$term == "std_contrast"], -2.25 / 2, tolerance = 1e-12)
})
