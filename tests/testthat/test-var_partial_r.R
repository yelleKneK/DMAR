# Tests for var_partial_r() and var_semipartial_r() against the closed
# forms from Olkin & Finn (1995), Hotelling (1953), and Cohen et al. (2003).

test_that("var_partial_r() matches Olkin-Finn (1995) closed form", {
  r <- 0.35; n <- 80; J <- 2
  expected <- (1 - r^2)^2 / (n - J - 1)
  expect_equal(var_partial_r(r = r, n = n, J = J)$value,
               expected, tolerance = 1e-12)
})

test_that("var_partial_r() Fisher's Z option matches Hotelling (1953) 1/(n-J-3)", {
  expect_equal(var_partial_r(r = 0.35, n = 80, J = 2, fisher_z = TRUE)$value,
               1 / (80 - 2 - 3), tolerance = 1e-12)
})

test_that("var_partial_r() reduces to simple-correlation result at J = 0... not allowed; check J = 1", {
  # At J = 1, partial r asymptotic variance should be (1 - r^2)^2 / (n - 2),
  # one fewer df than the simple-correlation form.
  expect_equal(var_partial_r(r = 0.5, n = 50, J = 1)$value,
               (1 - 0.5^2)^2 / (50 - 1 - 1), tolerance = 1e-12)
})

test_that("var_partial_r() rejects malformed input", {
  expect_error(var_partial_r(r = -1.5, n = 50, J = 1), "r")
  expect_error(var_partial_r(r = 0.5,  n = 3,  J = 2, fisher_z = TRUE),
               "Fisher")
  expect_error(var_partial_r(r = 0.5,  n = 50, J = 0), "J")
})

test_that("var_partial_r() variance scales like 1/(n - J - 1)", {
  v_50  <- var_partial_r(r = 0.4, n = 50,  J = 2)$value
  v_500 <- var_partial_r(r = 0.4, n = 500, J = 2)$value
  # Ratio of denominators: (500 - 3) / (50 - 3) = 497 / 47.
  expect_equal(v_50 / v_500, 497 / 47, tolerance = 1e-10)
})

test_that("var_semipartial_r() Olkin-Finn-style variance matches the closed form", {
  r_sp <- 0.25; n <- 100; J <- 3
  expected <- (1 - r_sp^2)^2 / (n - J - 1)
  expect_equal(var_semipartial_r(r_sp = r_sp, n = n, J = J)$value,
               expected, tolerance = 1e-12)
})

test_that("var_semipartial_r() with R2_full returns Cohen et al. (2003) null variance", {
  r_sp <- 0.25; n <- 100; J <- 3; R2 <- 0.42
  expected <- (1 - R2) / (n - J - 2)
  out <- var_semipartial_r(r_sp = r_sp, n = n, J = J, R2_full = R2)
  expect_equal(out$value, expected, tolerance = 1e-12)
  expect_equal(out$term,  "var_semipartial_r_under_null")
})

test_that("var_semipartial_r() rejects invalid R2_full", {
  expect_error(var_semipartial_r(r_sp = 0.2, n = 50, J = 1, R2_full = 1.1),
               "R2_full")
  expect_error(var_semipartial_r(r_sp = 0.2, n = 50, J = 1, R2_full = -0.1),
               "R2_full")
})
