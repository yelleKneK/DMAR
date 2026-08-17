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
              n = f$n, N = f$N, conf_level = 0.95)
  expect_s3_class(res, "data.frame")
  expect_named(res, c("term", "value"))
  expect_setequal(res$term, c("lower_limit", "contrast", "upper_limit"))
  expect_lt(res$value[res$term == "lower_limit"], res$value[res$term == "contrast"])
  expect_gt(res$value[res$term == "upper_limit"], res$value[res$term == "contrast"])
})

test_that("ci_sc() returns a tidy data.frame with a standardized-contrast CI", {
  f <- .fix()
  res <- ci_sc(means = f$means, s_anova = f$s_anova, c_weights = f$c_weights,
               n = f$n, N = f$N, conf_level = 0.95)
  expect_s3_class(res, "data.frame")
  expect_setequal(res$term, c("lower_limit", "std_contrast", "upper_limit"))
  expect_lt(res$value[res$term == "lower_limit"], res$value[res$term == "std_contrast"])
})

test_that("ci_c_ancova() returns a tidy ANCOVA-contrast CI bracketing psi", {
  res <- ci_c_ancova(adj_means = c(10, 12, 14), s_ancova = 2.3,
                     c_weights = c(-1, 0, 1),
                     n = c(20, 20, 20),
                     cov_means = c(5, 5.2, 4.9),
                     SSwithin_x = 50,
                     conf_level = 0.95)
  expect_s3_class(res, "data.frame")
  expect_named(res, c("term", "value"))
  expect_setequal(res$term, c("lower_limit", "psi", "upper_limit"))
  expect_lt(res$value[res$term == "lower_limit"], res$value[res$term == "psi"])
})

test_that("ci_sc_ancova() returns a tidy standardized ANCOVA-contrast CI", {
  res <- ci_sc_ancova(adj_means = c(10, 12, 14), s_anova = 2.5, s_ancova = 2.3,
                      c_weights = c(-1, 0, 1),
                      n = c(20, 20, 20),
                      cov_means = c(5, 5.2, 4.9),
                      SSwithin_x = 50,
                      conf_level = 0.95)
  expect_s3_class(res, "data.frame")
  expect_named(res, c("term", "value"))
  expect_setequal(res$term, c("lower_limit", "psi", "upper_limit"))
  expect_type(res$value, "double")
  expect_lt(res$value[res$term == "lower_limit"], res$value[res$term == "psi"])
  expect_identical(attr(res, "standardizer"), "s_ancova")
})

test_that("ci_sc_ancova(standardizer = 's_anova') records the alternative divisor", {
  res <- ci_sc_ancova(adj_means = c(10, 12, 14), s_anova = 2.5, s_ancova = 2.3,
                      c_weights = c(-1, 0, 1),
                      n = c(20, 20, 20),
                      cov_means = c(5, 5.2, 4.9),
                      SSwithin_x = 50,
                      standardizer = "s_anova",
                      conf_level = 0.95)
  expect_identical(attr(res, "standardizer"), "s_anova")
})
