test_that("cohen_f() returns a 1-row tidy data frame with term cohen_f", {
  result <- cohen_f(eta_squared = 0.10)
  expect_s3_class(result, "data.frame")
  expect_named(result, c("term", "value"))
  expect_equal(result$term, "cohen_f")
})

test_that("cohen_f() three calling modes agree on a common ground truth", {
  # Construct compatible inputs by hand:
  mu  <- c(94, 91, 92, 83)
  sig2 <- 67.375
  mu_bar     <- mean(mu)
  sigma_m_sq <- mean((mu - mu_bar)^2)
  f_raw   <- sqrt(sigma_m_sq / sig2)

  # Equivalent eta^2 = sigma_m^2 / (sigma_m^2 + sigma^2)
  eta2 <- sigma_m_sq / (sigma_m_sq + sig2)
  f_eta <- sqrt(eta2 / (1 - eta2))

  expect_equal(cohen_f(mu = mu, sigma_squared = sig2)$value, f_raw, tolerance = 1e-12)
  expect_equal(cohen_f(eta_squared = eta2)$value,             f_eta, tolerance = 1e-12)
  expect_equal(cohen_f(sigma_m = sqrt(sigma_m_sq), sigma = sqrt(sig2))$value, f_raw, tolerance = 1e-12)
  # All three forms agree:
  expect_equal(f_raw, f_eta, tolerance = 1e-12)
})

test_that("cohen_f() unequal n weighting matches the SS_between formulation", {
  mu  <- c(94, 91, 92, 83)
  sig2 <- 67.375
  n   <- c(4, 6, 5, 5)
  N   <- sum(n)
  mu_bar <- sum(n * mu) / N
  SS_between <- sum(n * (mu - mu_bar)^2)
  expected_f <- sqrt((SS_between / N) / sig2)
  expect_equal(cohen_f(mu = mu, sigma_squared = sig2, n = n)$value, expected_f,
               tolerance = 1e-12)
})

test_that("cohen_f() rejects ambiguous or under-specified calls", {
  expect_error(cohen_f(),                              "exactly one of")
  expect_error(cohen_f(eta_squared = 0.1, sigma_m = 1, sigma = 1), "exactly one of")
  expect_error(cohen_f(mu = c(1, 2)),                  "exactly one of")
  expect_error(cohen_f(eta_squared = -0.1),            "in \\[0, 1\\)")
  expect_error(cohen_f(eta_squared = 1),               "in \\[0, 1\\)")
  expect_error(cohen_f(sigma_m = 1, sigma = 0),        "must be positive")
})
