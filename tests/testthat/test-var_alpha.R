test_that("var_alpha() returns van Zyl and Bonett rows", {
  res <- var_alpha(alpha = 0.80, N = 100, n_items = 10)
  expect_setequal(res$term, c("var_alpha_van_zyl", "var_alpha_bonett"))
  expect_true(all(res$value > 0))
})

test_that("var_alpha() variance shrinks with larger n", {
  small <- var_alpha(0.90, N = 30,  n_items = 5)$value[1]
  large <- var_alpha(0.90, N = 300, n_items = 5)$value[1]
  expect_gt(small, large)
})

test_that("var_alpha() approximations converge as n grows", {
  out <- var_alpha(0.90, N = 1000, n_items = 5)
  expect_equal(out$value[1], out$value[2], tolerance = 1e-3)
})

test_that("var_alpha() rejects bad inputs", {
  expect_error(var_alpha(1.2, 100, 10), "in \\[0, 1\\)")
  expect_error(var_alpha(0.8, 2,   10), ">= 3")
  expect_error(var_alpha(0.8, 100, 1),  ">= 2")
})

test_that("var_alpha() refuses the 1.0.x name n, which R would match to n_items", {
  # (suppressWarnings: under warnPartialMatchArgs R warns of the partial
  # match before the function refuses it.)
  err <- suppressWarnings(tryCatch(var_alpha(0.80, n = 100, 10), error = function(e) e))
  expect_s3_class(err, "dmar_argument_invalid")
  expect_identical(err$arg, "N")
  expect_match(conditionMessage(err), "N = ", fixed = TRUE)
  expect_error(suppressWarnings(do.call(var_alpha, list(0.80, n = 100, 10))),
               class = "dmar_argument_invalid")
  # By position and by full name the call runs.
  expect_equal(var_alpha(0.80, 100, 10)$value,
               var_alpha(alpha = 0.80, N = 100, n_items = 10)$value)
})

test_that("var_alpha() rows are the van Zyl et al. (2000) and Bonett forms", {
  # Oracles, from the articles: van Zyl, Neudecker, and Nel (2000, Eq. 13),
  # lim n Var(R) = 2 k^2 p / (p - 1) with k = 1 - alpha, p items, and
  # n = N - 1 the degrees of freedom of the Wishart matrix A (p. 272); and
  # Bonett (2003, Eqs. 2 and 4, citing Bonett, 2002), the variance of
  # ln(1 - alpha_hat) is 2k / ((k - 1)(n - 2)) with k parts and n
  # respondents, which the delta method carries to the alpha scale by the
  # factor (1 - alpha)^2. DMAR 1.0.x had the two in each other's rows.
  for (a in c(0, .5, .8, .95)) for (N in c(3, 10, 57, 400))
    for (J in c(2, 5, 12)) {
      res <- var_alpha(alpha = a, N = N, n_items = J)
      k <- 1 - a; n <- N - 1
      expect_identical(res$term, c("var_alpha_van_zyl", "var_alpha_bonett"))
      expect_equal(res$value[1], 2 * k^2 * J / ((J - 1) * n),
                   tolerance = 1e-14)
      expect_equal(res$value[2], (1 - a)^2 * 2 * J / ((J - 1) * (N - 2)),
                   tolerance = 1e-14)
    }
})

test_that("var_alpha() van Zyl row equals the general van Zyl variance at a compound symmetric matrix", {
  # The general asymptotic variance of van Zyl et al. (2000, Eq. 21), as
  # reliability_alpha(ci_method = "ml") evaluates it, divided by N - 1,
  # reduces at a compound symmetric covariance matrix to the row
  # var_alpha_van_zyl; the item variance (here 4) does not matter.
  for (J in c(3, 6)) for (rho in c(.2, .5)) {
    S <- matrix(rho, J, J); diag(S) <- 1; S <- 4 * S
    a <- J * rho / (1 + (J - 1) * rho)
    expect_equal(var_alpha(alpha = a, N = 50, n_items = J)$value[1],
                 DMAR:::.alpha_vanzyl_avar(S) / 49, tolerance = 1e-12)
  }
})

test_that("var_alpha() help-page comparison with the exact variance holds", {
  # Oracle for the exact variance under normality and compound symmetry:
  # R = 1 - kF with F ~ F(n(J - 1), n), n = N - 1 (van Zyl et al., 2000,
  # Section 2), so Var(R) = k^2 Var(F), with Var(F) from numerical
  # integration of the F density, independent of the closed form (Eq. 12)
  # the page quotes.
  var_F <- function(d1, d2) {
    m1 <- integrate(function(x) x * df(x, d1, d2), 0, Inf,
                    rel.tol = 1e-12)$value
    m2 <- integrate(function(x) x^2 * df(x, d1, d2), 0, Inf,
                    rel.tol = 1e-12)$value
    m2 - m1^2
  }
  closed <- function(a, N, J) {
    n <- N - 1
    2 * (1 - a)^2 * n * (n * J - 2) / ((J - 1) * (n - 2)^2 * (n - 4))
  }
  for (N in c(6, 10, 30, 100)) for (J in c(2, 5, 20)) {
    a <- .7
    exact <- (1 - a)^2 * var_F((N - 1) * (J - 1), N - 1)
    expect_equal(closed(a, N, J), exact, tolerance = 1e-9)
    v <- var_alpha(alpha = a, N = N, n_items = J)$value
    # Larger than both forms; the N - 2 form is the closer.
    expect_gt(exact, v[2])
    expect_gt(v[2], v[1])
    expect_equal(v[2] / v[1], (N - 1) / (N - 2), tolerance = 1e-14)
  }
  # The ratios the page quotes for J = 5.
  ratio <- function(N) {
    closed(.7, N, 5) / var_alpha(alpha = .7, N = N, n_items = 5)$value[2]
  }
  expect_equal(round(c(ratio(10), ratio(30), ratio(100)), 2),
               c(2.53, 1.27, 1.07))
})

test_that("var_alpha() N for a target variance is one larger with N - 2 than with N - 1", {
  # The smallest N with c / (N - 2) <= v is ceiling(2 + c / v); with
  # c / (N - 1) <= v it is ceiling(1 + c / v), one less, unless the
  # minimum N = 3 binds.
  for (a in c(.6, .8, .9)) for (J in c(4, 10)) for (target in c(.01, .002)) {
    smallest <- function(row) {
      N <- 3
      while (var_alpha(alpha = a, N = N, n_items = J)$value[row] > target)
        N <- N + 1
      N
    }
    expect_identical(smallest(2) - smallest(1), 1)
  }
})

test_that("var_alpha() gives no warning", {
  expect_no_warning(var_alpha(alpha = .8, N = 100, n_items = 10))
})
