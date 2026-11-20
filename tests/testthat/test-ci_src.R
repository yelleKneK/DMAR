test_that("ci_src() returns a 3-row data frame with the estimate between the limits", {
  result <- ci_src(std_reg_coef_j = 0.4, se_std_reg_coef_j = 0.1, N = 100, p = 5)
  expect_s3_class(result, "data.frame")
  expect_named(result, c("term", "value", "prob_less", "prob_greater"))
  expect_equal(result$term, c("ci_lower", "src", "ci_upper"))
  expect_equal(result$value[result$term == "src"], 0.4)
})

test_that("ci_src() interval brackets the standardized point estimate", {
  beta_j <- 0.4
  result <- ci_src(std_reg_coef_j = beta_j, se_std_reg_coef_j = 0.1, N = 100, p = 5)
  ll <- result$value[result$term == "ci_lower"]
  ul <- result$value[result$term == "ci_upper"]
  expect_lt(ll, beta_j)
  expect_gt(ul, beta_j)
})

# The noncentral t limits of a standardized coefficient (Kelley, 2007),
# coded without the package: the limits for the noncentrality parameter
# are the values at which the observed t is the 1 - alpha/2 and alpha/2
# quantiles, by root finding on stats::pt, and the limits for the
# coefficient are those times the standard error.
.nc_coef_limits <- function(est, se, df, conf_level = .95) {
  t_obs <- est / se
  a <- (1 - conf_level) / 2
  f <- function(ncp, target) suppressWarnings(pt(t_obs, df, ncp)) - target
  lo <- uniroot(f, t_obs + c(-10, 10), target = 1 - a, tol = 1e-12)$root
  hi <- uniroot(f, t_obs + c(-10, 10), target = a, tol = 1e-12)$root
  c(lo, hi) * se
}

test_that("ci_src() takes the coefficient from the squared multiple correlations, as its help says", {
  # Through the first 1.1.0 builds this route stopped with R's "argument is
  # of length zero" at the check of the coefficient's size.
  N <- 100; p <- 3
  R2 <- .30; R2_without <- .20; R2_j <- .10
  beta <- sqrt((R2 - R2_without) / (1 - R2_j))
  se <- sqrt((1 - R2) / ((1 - R2_j) * (N - p - 1)))
  seen <- new.env(parent = emptyenv())
  res <- withCallingHandlers(
    ci_src(N = N, p = p, R2_Y_X = R2, R2_Y_X_without_j = R2_without,
           R2_j_X_without_j = R2_j),
    dmar_result_qualified = function(x) {
      seen$w <- x
      invokeRestart("muffleWarning")
    })
  w <- seen$w
  expect_equal(res$term, c("ci_lower", "src", "ci_upper"))
  expect_equal(res$value[2], beta, tolerance = 1e-12)
  expect_equal(res$value[c(1, 3)], .nc_coef_limits(beta, se, N - p - 1),
               tolerance = 1e-8)
  # The warning that the sign is not identified names this function's
  # arguments, not those of ci_reg_coef().
  expect_s3_class(w, "dmar_result_qualified")
  expect_match(conditionMessage(w), "'std_reg_coef_j'", fixed = TRUE)
  expect_match(conditionMessage(w), "'se_std_reg_coef_j'", fixed = TRUE)
  expect_false(grepl("(^|[^_])reg_coef_j", conditionMessage(w)))
})

test_that("ci_src() takes the coefficient from t and its standard error, with its sign", {
  N <- 120; p <- 4
  R2 <- .25; R2_j <- .15
  se <- sqrt((1 - R2) / ((1 - R2_j) * (N - p - 1)))
  for (t in c(-2.7, 3.1)) {
    by_t_R2 <- ci_src(t_value = t, N = N, p = p, R2_Y_X = R2, R2_j_X_without_j = R2_j)
    expect_equal(by_t_R2$value[2], t * se, tolerance = 1e-12)
    expect_equal(by_t_R2$value[c(1, 3)], .nc_coef_limits(t * se, se, N - p - 1),
                 tolerance = 1e-8)
    by_t_se <- ci_src(t_value = t, se_std_reg_coef_j = se, N = N, p = p)
    expect_equal(by_t_se$value, by_t_R2$value, tolerance = 1e-10)
    by_coef <- ci_src(std_reg_coef_j = t * se, se_std_reg_coef_j = se, N = N, p = p)
    expect_equal(by_coef$value, by_t_R2$value, tolerance = 1e-10)
  }
  # Reversing the sign of the coefficient reflects the interval.
  pos <- ci_src(std_reg_coef_j = .3, se_std_reg_coef_j = .1, N = N, p = p)$value
  neg <- ci_src(std_reg_coef_j = -.3, se_std_reg_coef_j = .1, N = N, p = p)$value
  expect_equal(neg, -rev(pos), tolerance = 1e-10)
})

test_that("ci_src() refuses inputs that determine no interval, by name", {
  e <- tryCatch(ci_src(N = 100, p = 3, R2_Y_X = .3), error = function(e) e)
  expect_s3_class(e, "dmar_argument_invalid")
  expect_identical(e$arg, "std_reg_coef_j")
  e <- tryCatch(ci_src(std_reg_coef_j = .3, N = 100, p = 3), error = function(e) e)
  expect_s3_class(e, "dmar_argument_invalid")
  expect_true("se_std_reg_coef_j" %in% e$arg)
  e <- tryCatch(ci_src(std_reg_coef_j = .3, se_std_reg_coef_j = .1), error = function(e) e)
  expect_s3_class(e, "dmar_argument_invalid")
  expect_identical(e$arg, c("N", "p"))
  e <- tryCatch(ci_src(reg_coef_j = .3, se_reg_coef_j = .1, N = 100, p = 3, sd_X = 1),
                error = function(e) e)
  expect_identical(e$arg, "sd_Y")
})

test_that("ci_src() qualifies a large standardized coefficient of either sign", {
  expect_warning(ci_src(std_reg_coef_j = 1.4, se_std_reg_coef_j = .1, N = 100, p = 3),
                 class = "dmar_result_qualified")
  expect_warning(ci_src(std_reg_coef_j = -1.4, se_std_reg_coef_j = .1, N = 100, p = 3),
                 class = "dmar_result_qualified")
  expect_no_warning(ci_src(std_reg_coef_j = -.9, se_std_reg_coef_j = .1, N = 100, p = 3))
})

test_that("ci_src() from the unstandardized coefficient and the standard deviations", {
  # b = 2, sd_X = 1.5, sd_Y = 6: the standardized coefficient is b sd_X / sd_Y = .5.
  a <- ci_src(reg_coef_j = 2, se_reg_coef_j = .4, sd_Y = 6, sd_X = 1.5, N = 80, p = 2)
  b <- ci_src(std_reg_coef_j = .5, se_std_reg_coef_j = .1, N = 80, p = 2)
  expect_equal(a$value, b$value, tolerance = 1e-12)
})
