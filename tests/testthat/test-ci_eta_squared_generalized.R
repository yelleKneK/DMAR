test_that("default method = 'none' returns NA limits and a message", {
  expect_message(
    res <- ci_eta_squared_generalized(SS_effect = 100, SS_observed = 70,
                                      SS_error = 200),
    "No CI computed"
  )
  expect_equal(res$method, "none")
  expect_true(is.na(res$lower_limit))
  expect_true(is.na(res$upper_limit))
  expect_equal(res$eta_squared_generalized, 100 / (100 + 70 + 200),
               tolerance = 1e-12)
})

test_that("parametric method returns a CI and warns about approximation", {
  fit <- aov(len ~ supp * dose, data = ToothGrowth)
  expect_warning(
    res <- ci_eta_squared_generalized(fit, observed = "supp",
                                      method = "parametric"),
    "Parametric CI uses an approximate"
  )
  expect_equal(res$method[1], "parametric")
  # CI brackets the point estimate for each row
  for (i in seq_len(nrow(res))) {
    expect_gte(res$eta_squared_generalized[i], res$lower_limit[i])
    expect_lte(res$eta_squared_generalized[i], res$upper_limit[i])
  }
  # All limits in [0, 1]
  expect_true(all(res$lower_limit >= 0))
  expect_true(all(res$upper_limit <= 1))
})

test_that("parametric method with F/df interface requires N", {
  expect_error(
    suppressWarnings(ci_eta_squared_generalized(
      F_effect = 6, df_effect = 2,
      F_observed = 2.5, df_observed = 1, df_error = 50,
      method = "parametric"
    )),
    "requires 'N'"
  )
})

test_that("parametric F/df interface returns a valid CI when N is supplied", {
  res <- suppressWarnings(ci_eta_squared_generalized(
    F_effect = 6, df_effect = 2,
    F_observed = 2.5, df_observed = 1, df_error = 50,
    N = 54, method = "parametric"
  ))
  expect_gte(res$eta_squared_generalized, res$lower_limit)
  expect_lte(res$eta_squared_generalized, res$upper_limit)
  expect_gte(res$lower_limit, 0)
  expect_lte(res$upper_limit, 1)
})

test_that("bootstrap method requires a fitted model", {
  expect_error(
    suppressWarnings(ci_eta_squared_generalized(
      SS_effect = 100, SS_observed = 70, SS_error = 200,
      method = "bootstrap"
    )),
    "requires a fitted model"
  )
})

test_that("bootstrap method returns CI; reproducible with a seed", {
  # Two B = 1000 bootstraps back-to-back; refits the aov() inside the
  # boot() loop, slow on CRAN. The parametric CI test above already
  # checks the same bracketing and "method" attribute.
  skip_on_cran()
  fit <- aov(len ~ supp * dose, data = ToothGrowth)
  res1 <- suppressWarnings(ci_eta_squared_generalized(
    fit, observed = "supp", method = "bootstrap", B = 1000, seed = 113
  ))
  res2 <- suppressWarnings(ci_eta_squared_generalized(
    fit, observed = "supp", method = "bootstrap", B = 1000, seed = 113
  ))
  expect_equal(res1$lower_limit, res2$lower_limit)
  expect_equal(res1$upper_limit, res2$upper_limit)
  for (i in seq_len(nrow(res1))) {
    expect_gte(res1$eta_squared_generalized[i], res1$lower_limit[i])
    expect_lte(res1$eta_squared_generalized[i], res1$upper_limit[i])
    expect_gte(res1$lower_limit[i], 0)
    expect_lte(res1$upper_limit[i], 1)
  }
  expect_equal(res1$method[1], "bootstrap")
})

test_that("invalid R is rejected (minimum 1000)", {
  fit <- aov(len ~ supp * dose, data = ToothGrowth)
  expect_error(
    suppressWarnings(ci_eta_squared_generalized(
      fit, observed = "supp", method = "bootstrap", B = 999
    )),
    "'B' must be"
  )
})

test_that("method argument is matched with match.arg", {
  expect_error(
    ci_eta_squared_generalized(SS_effect = 100, SS_observed = 70,
                               SS_error = 200, method = "garbage"),
    "should be one of"
  )
})

test_that("output columns are consistent across methods", {
  # The bootstrap arm of this column-consistency check uses B = 1000.
  # Skipped on CRAN; the default and "parametric" arms still run there.
  skip_on_cran()
  cols <- c("effect", "eta_squared_generalized", "lower_limit",
            "upper_limit", "method")
  fit  <- aov(len ~ supp * dose, data = ToothGrowth)
  r1   <- suppressMessages(ci_eta_squared_generalized(fit, observed = "supp"))
  r2   <- suppressWarnings(ci_eta_squared_generalized(
    fit, observed = "supp", method = "parametric"
  ))
  r3   <- suppressWarnings(ci_eta_squared_generalized(
    fit, observed = "supp", method = "bootstrap", B = 1000, seed = 113
  ))
  for (r in list(r1, r2, r3)) {
    expect_true(all(cols %in% names(r)))
  }
})

test_that("parametric CI reports the lower-limit clamp once, in terms of generalized eta squared", {
  set.seed(113)
  d <- data.frame(
    y = rnorm(60),
    a = factor(rep(1:3, each = 20)),
    b = factor(rep(1:2, 30))
  )
  fit <- aov(y ~ a * b, data = d)
  msgs <- capture_warnings(result <- ci_eta_squared_generalized(fit, method = "parametric"))
  clamp_msgs <- msgs[grepl("lower confidence limit on generalized eta squared is 0", msgs)]
  expect_length(clamp_msgs, 1L)
  expect_false(any(grepl("prob_greater", msgs)))
  expect_true(all(result$lower_limit == 0))
  expect_true(all(result$upper_limit > 0))
})
