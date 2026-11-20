# Tests for var_partial_r() and var_semipartial_r() against their closed
# forms, and for the first-order variance of the semipartial correlation
# against a delta method computed here from the correlation matrix.

test_that("var_partial_r() matches the first-order closed form", {
  r <- 0.35; n <- 80; J <- 2
  expected <- (1 - r^2)^2 / (n - J - 1)
  expect_equal(var_partial_r(r = r, N = n, p = J + 1)$value,
               expected, tolerance = 1e-12)
})

test_that("var_partial_r() Fisher's Z option is 1/(N - p - 2)", {
  expect_equal(var_partial_r(r = 0.35, N = 80, p = 3, fisher_z = TRUE)$value,
               1 / (80 - 3 - 2), tolerance = 1e-12)
})

test_that("var_partial_r() at p = 2 has one fewer df than the simple correlation form", {
  # At p = 2 (one variable partialed out) the partial r asymptotic variance
  # is (1 - r^2)^2 / (N - 2), one fewer df than the simple-correlation form.
  expect_equal(var_partial_r(r = 0.5, N = 50, p = 2)$value,
               (1 - 0.5^2)^2 / (50 - 1 - 1), tolerance = 1e-12)
})

test_that("var_partial_r() rejects malformed input", {
  expect_error(var_partial_r(r = -1.5, N = 50, p = 2), "r")
  expect_error(var_partial_r(r = 0.5,  N = 3,  p = 3, fisher_z = TRUE),
               "Fisher")
  expect_error(var_partial_r(r = 0.5,  N = 50, p = 1),
               "single whole number of at least 2. With p = 1 nothing is partialed out",
               fixed = TRUE)
  expect_error(var_partial_r(r = 0.5,  N = 50, p = 2.5),
               "whole number of at least 2")
  expect_error(var_partial_r(r = 0.5,  N = 50, p = NA_real_),
               "whole number of at least 2")
  expect_error(var_partial_r(r = 0.5,  N = 50, p = 0),
               "whole number of at least 2")
  # The messages state the inequalities in N and p, and the boundaries are
  # where they say: the raw variance needs N > p, Fisher's Z N > p + 2.
  expect_error(var_partial_r(r = 0.5, N = 4, p = 4),
               "Need N > p for the first-order variance", fixed = TRUE)
  expect_no_error(var_partial_r(r = 0.5, N = 5, p = 4))
  expect_error(var_partial_r(r = 0.5, N = 6, p = 4, fisher_z = TRUE),
               "Need N > p + 2 for the Fisher's Z variance.", fixed = TRUE)
  expect_no_error(var_partial_r(r = 0.5, N = 7, p = 4, fisher_z = TRUE))
  expect_error(var_partial_r(r = 0.5, N = -1, p = 2),
               "'N' must be a single positive number.", fixed = TRUE)
})

test_that("var_partial_r() and var_semipartial_r() take p in its place, by position or by name", {
  # p is the number of predictors, an ordinary third argument.
  expect_identical(var_partial_r(.3, 100, 3),
                   var_partial_r(r = .3, N = 100, p = 3))
  expect_identical(var_partial_r(.3, 100), var_partial_r(.3, 100, p = 2))
  expect_identical(var_semipartial_r(.3, 100, 3, R2_full = .4),
                   var_semipartial_r(r_sp = .3, N = 100, p = 3, R2_full = .4))
  expect_identical(do.call(var_partial_r, list(.3, 100, 3)),
                   var_partial_r(r = .3, N = 100, p = 3))
  expect_identical(Vectorize(var_partial_r)(.3, 100, 3)[, 1]$value,
                   var_partial_r(r = .3, N = 100, p = 3)$value)
  # The count is the number of predictors: p - 1 variables are partialed
  # out, so the divisor is N - p.
  expect_equal(var_partial_r(.3, 100, 3)$value, (1 - .3^2)^2 / (100 - 3),
               tolerance = 1e-14)
})

test_that("var_partial_r() and var_semipartial_r() do not accept the DMAR 1.0.x names", {
  expect_error(var_partial_r(r = .3, N = 100, J = 2), "unused argument")
  expect_error(var_partial_r(r = .3, n = 100, p = 3), "unused argument")
  expect_error(var_semipartial_r(r_sp = .3, N = 100, J = 2, R2_full = .4),
               "unused argument")
  expect_error(var_semipartial_r(r_sp = .3, n = 100, p = 3, R2_full = .4),
               "unused argument")
  expect_identical(names(formals(var_partial_r)),
                   c("r", "N", "p", "fisher_z"))
  expect_identical(names(formals(var_semipartial_r)),
                   c("r_sp", "N", "p", "R2_full", "under_null"))
})

test_that("var_semipartial_r(under_null = TRUE) equals r_sp^2 / t^2 from a regression fit", {
  # The identity r_sp = t sqrt((1 - R2_full) / (N - p - 1)) (Aloe & Becker,
  # 2012, Equation 2), from data generated here and fitted by lm().
  set.seed(113)
  for (dims in list(c(60, 2), c(137, 5), c(500, 9))) {
    N <- dims[1]; p <- dims[2]
    Z <- matrix(rnorm(N * (p - 1)), N, p - 1)
    X <- as.numeric(0.4 * Z[, 1] + rnorm(N))
    Y <- as.numeric(0.5 * X + Z %*% rep(0.3, p - 1) + rnorm(N))
    fit <- lm(Y ~ X + Z)
    t_x <- summary(fit)$coefficients["X", "t value"]
    R2 <- summary(fit)$r.squared
    r_sp <- cor(Y, residuals(lm(X ~ Z)))
    expect_equal(r_sp, t_x * sqrt((1 - R2) / (N - p - 1)), tolerance = 1e-10)
    expect_equal(var_semipartial_r(r_sp, N, p = p, R2_full = R2)$value,
                 r_sp^2 / t_x^2, tolerance = 1e-10)
    expect_equal(r_sp^2, R2 - summary(lm(Y ~ Z))$r.squared,
                 tolerance = 1e-10)
  }
})

test_that("var_partial_r() variance scales like 1/(N - p)", {
  v_50  <- var_partial_r(r = 0.4, N = 50,  p = 3)$value
  v_500 <- var_partial_r(r = 0.4, N = 500, p = 3)$value
  # Ratio of denominators: (500 - 3) / (50 - 3) = 497 / 47.
  expect_equal(v_50 / v_500, 497 / 47, tolerance = 1e-10)
})

test_that("var_semipartial_r() without R2_full keeps its value and warns by class", {
  r_sp <- 0.25; n <- 100; J <- 3
  expected <- (1 - r_sp^2)^2 / (n - J - 1)
  expect_warning(
    out <- var_semipartial_r(r_sp = r_sp, N = n, p = J + 1),
    class = "dmar_semipartial_r_no_R2_full")
  expect_equal(out$value, expected, tolerance = 1e-12)
  expect_equal(out$term, "var_semipartial_r")
  # under_null has no effect without R2_full.
  expect_warning(
    out_f <- var_semipartial_r(r_sp = r_sp, N = n, p = J + 1, under_null = FALSE),
    class = "dmar_semipartial_r_no_R2_full")
  expect_identical(out_f$value, out$value)
})

test_that("var_semipartial_r() gives the first-order variance with R2_full and under_null = FALSE", {
  # V = (1 - r^2)^2 - q (1 - 2 r^2), q = R2_full - r^2: at r = .3 and
  # R2_full = .59, .8281 - .50 * .82 = .4181.
  out <- var_semipartial_r(0.3, 200, p = 4, R2_full = 0.59, under_null = FALSE)
  expect_equal(out$value, 0.4181 / 196, tolerance = 1e-12)
  expect_equal(out$term, "var_semipartial_r")

  # With R2_full = r_sp^2 (controls unrelated to the outcome) it is the
  # value returned without R2_full.
  no_R2 <- suppressWarnings(var_semipartial_r(0.3, 200, p = 4))$value
  expect_equal(
    var_semipartial_r(0.3, 200, p = 4, R2_full = 0.09, under_null = FALSE)$value,
    no_R2, tolerance = 1e-12)
  expect_equal(
    var_semipartial_r(0.3, 200, p = 4, R2_full = 0.3^2, under_null = FALSE)$value,
    no_R2, tolerance = 1e-12)

  # At r_sp = 0 it is (1 - R2_full) / (n - J - 1).
  expect_equal(
    var_semipartial_r(0, 500, p = 2, R2_full = 0.64, under_null = FALSE)$value,
    0.36 / 498, tolerance = 1e-12)

  # Above r_sp^2 = 1/2 the variance grows with the R2 of the controls.
  expect_equal(
    var_semipartial_r(0.8, 100, p = 2, R2_full = 0.94, under_null = FALSE)$value,
    0.2136 / 98, tolerance = 1e-12)
})

test_that("var_semipartial_r() first-order variance equals a delta method on the correlation matrix", {
  # Independent route: the sample semipartial correlation is a function of
  # the sample correlations, whose normal-theory covariances are the
  # Pearson-Filon expressions; its gradient is taken by central differences.
  # The three populations have controls correlated with each other, with X,
  # and with Y.
  semipartial <- function(P) {
    z <- 3:nrow(P)
    Pzz_inv <- solve(P[z, z, drop = FALSE])
    drop((P[1, 2] - P[1, z] %*% Pzz_inv %*% P[z, 2]) /
           sqrt(1 - P[2, z] %*% Pzz_inv %*% P[z, 2]))
  }
  R2_full <- function(P) {
    idx <- 2:nrow(P)
    drop(P[1, idx, drop = FALSE] %*% solve(P[idx, idx]) %*% P[idx, 1])
  }
  cov_rr <- function(P, i, j, k, l) {
    0.5 * P[i, j] * P[k, l] *
      (P[i, k]^2 + P[i, l]^2 + P[j, k]^2 + P[j, l]^2) +
      P[i, k] * P[j, l] + P[i, l] * P[j, k] -
      (P[i, j] * P[i, k] * P[i, l] + P[j, i] * P[j, k] * P[j, l] +
         P[k, i] * P[k, j] * P[k, l] + P[l, i] * P[l, j] * P[l, k])
  }
  delta_variance <- function(P) {
    pairs <- which(upper.tri(P), arr.ind = TRUE)
    m <- nrow(pairs); h <- 1e-6
    grad <- vapply(seq_len(m), function(a) {
      up <- P; down <- P
      i <- pairs[a, 1]; j <- pairs[a, 2]
      up[i, j] <- up[j, i] <- P[i, j] + h
      down[i, j] <- down[j, i] <- P[i, j] - h
      (semipartial(up) - semipartial(down)) / (2 * h)
    }, numeric(1))
    Omega <- matrix(0, m, m)
    for (a in seq_len(m)) for (b in seq_len(m))
      Omega[a, b] <- cov_rr(P, pairs[a, 1], pairs[a, 2],
                            pairs[b, 1], pairs[b, 2])
    drop(t(grad) %*% Omega %*% grad)
  }
  P1 <- matrix(c(1, .5, .6,
                 .5, 1, .3,
                 .6, .3, 1), 3, byrow = TRUE)
  P2 <- matrix(c(1, .45, .40, .30,
                 .45, 1, .35, .20,
                 .40, .35, 1, .50,
                 .30, .20, .50, 1), 4, byrow = TRUE)
  P3 <- matrix(c(1, .80, .30, .25, .20,
                 .80, 1, .15, .10, .25,
                 .30, .15, 1, .40, .30,
                 .25, .10, .40, 1, .35,
                 .20, .25, .30, .35, 1), 5, byrow = TRUE)
  n <- 400
  for (P in list(P1, P2, P3)) {
    J <- nrow(P) - 2
    r_sp <- semipartial(P)
    R2 <- R2_full(P)
    dmar <- var_semipartial_r(r_sp, n, p = J + 1, R2_full = R2,
                              under_null = FALSE)$value * (n - J - 1)
    expect_equal(dmar, delta_variance(P), tolerance = 1e-7,
                 info = paste("J =", J))
    # The controls matter: the variance is not (1 - r_sp^2)^2 here.
    expect_gt(abs(dmar - (1 - r_sp^2)^2), 0.01)
  }
})

test_that("var_semipartial_r() refuses an R2_full below r_sp squared, by class", {
  expect_error(var_semipartial_r(0.9, 100, p = 2, R2_full = 0.1),
               class = "dmar_semipartial_r_incoherent")
  expect_error(var_semipartial_r(0.9, 100, p = 2, R2_full = 0.1,
                                 under_null = FALSE),
               class = "dmar_refusal")
  cond <- tryCatch(var_semipartial_r(0.9, 100, p = 2, R2_full = 0.1),
                   dmar_refusal = function(e) e)
  expect_equal(cond$code, "semipartial_r_incoherent")
  expect_match(conditionMessage(cond), "use R2_full = r_sp^2", fixed = TRUE)
  # R2_full equal to r_sp^2, as floating point gives it, is coherent.
  expect_silent(var_semipartial_r(0.8, 100, p = 2, R2_full = 0.64))
})

test_that("var_semipartial_r() checks p, N, and under_null", {
  expect_error(var_semipartial_r(0.3, 50, p = 2.5, R2_full = 0.4),
               "whole number of at least 2")
  expect_error(var_semipartial_r(0.3, 50, p = 1, R2_full = 0.4),
               "whole number of at least 2")
  expect_error(var_semipartial_r(0.3, 50, p = 0, R2_full = 0.4),
               "whole number of at least 2")
  expect_error(var_semipartial_r(0.3, NA_real_, p = 2, R2_full = 0.4),
               "'N' must be a single finite positive number.", fixed = TRUE)
  # The boundaries the messages state: N > p + 1 for the t test's squared
  # standard error, N > p for the asymptotic variance.
  expect_error(var_semipartial_r(0.3, 5, p = 4, R2_full = 0.4),
               "Need N > p + 1 for the squared standard error of the t test.",
               fixed = TRUE)
  expect_no_error(var_semipartial_r(0.3, 6, p = 4, R2_full = 0.4))
  expect_error(var_semipartial_r(0.3, 4, p = 4, R2_full = 0.4,
                                 under_null = FALSE),
               "Need N > p for the asymptotic variance.", fixed = TRUE)
  expect_no_error(var_semipartial_r(0.3, 5, p = 4, R2_full = 0.4,
                                    under_null = FALSE))
  expect_error(var_semipartial_r(0.3, 50, p = 2, R2_full = 0.4,
                                 under_null = NA), "under_null")
  expect_error(var_semipartial_r(0.3, 50, p = 2, R2_full = 0.4,
                                 under_null = "no"), "under_null")
})

test_that("var_semipartial_r() with R2_full returns Cohen et al. (2003) null variance", {
  r_sp <- 0.25; n <- 100; J <- 3; R2 <- 0.42
  expected <- (1 - R2) / (n - J - 2)
  out <- var_semipartial_r(r_sp = r_sp, N = n, p = J + 1, R2_full = R2)
  expect_equal(out$value, expected, tolerance = 1e-12)
  expect_equal(out$term,  "var_semipartial_r_under_null")
})

test_that("var_semipartial_r() rejects invalid R2_full", {
  expect_error(var_semipartial_r(r_sp = 0.2, N = 50, p = 2, R2_full = 1.1),
               "R2_full")
  expect_error(var_semipartial_r(r_sp = 0.2, N = 50, p = 2, R2_full = -0.1),
               "R2_full")
})

test_that("var_semipartial_r() refuses R2_full below r_sp squared and passes a rounded pair", {
  # A pair that differs by 1e-10 or 1e-9 is R2_full = r_sp^2 in floating
  # point: it passes silently, and the variance is the one at q = 0, not
  # one computed from a negative q.
  expect_no_error(var_semipartial_r(.3, 100, p = 2, R2_full = .09 - 1e-10))
  expect_equal(var_semipartial_r(.3, 100, p = 2, R2_full = .09 - 1e-9,
                                 under_null = FALSE)$value,
               (1 - .09)^2 / 98, tolerance = 1e-12)
})

test_that("var_semipartial_r() reads a pair that two-decimal rounding can produce, with a classed warning", {
  # .46 squared is .2116 and a full-model R2 reported as .21 is below it,
  # but .21 + .005 >= (.46 - .005)^2 = .207025, so the pair is consistent
  # with rounding. DMAR 1.0.0 returned .00814433 for this call; the first
  # builds of 1.1.0 refused it.
  w <- NULL
  out <- withCallingHandlers(
    var_semipartial_r(.46, 100, p = 2, R2_full = .21),
    dmar_semipartial_r_rounded = function(cnd) {
      w <<- cnd
      invokeRestart("muffleWarning")
    })
  expect_equal(out$value, 0.00814433, tolerance = 1e-7)
  expect_equal(out$value, (1 - .21) / 97, tolerance = 1e-14)
  expect_s3_class(w, c("dmar_semipartial_r_rounded", "warning", "condition"),
                  exact = TRUE)
  expect_identical(w$R2_full, .21)
  expect_identical(w$r_sp, .46)
  expect_match(conditionMessage(w), "'R2_full' (0.21) is below 'r_sp' squared (0.2116) by no more than rounding to two decimals allows",
               fixed = TRUE)
  # The squared standard error of the t test uses R2_full as given and
  # does not involve r_sp, and the warning says that; it does not say that
  # R2_full is taken as r_sp^2, which would give (1 - .2116) / 97.
  expect_match(conditionMessage(w),
               "The value returned does not depend on r_sp and uses R2_full as given.",
               fixed = TRUE)
  expect_no_match(conditionMessage(w), "(R2_full = r_sp^2)", fixed = TRUE)
  expect_false(isTRUE(all.equal(out$value, (1 - .46^2) / 97,
                                tolerance = 1e-6)))
  # The count by position is p, the number of predictors (the 1.0.x call
  # gave J = 1 there; see the breaking changes in NEWS).
  expect_equal(
    suppressWarnings(var_semipartial_r(.46, 100, 2, R2_full = .21))$value,
    out$value, tolerance = 1e-14)
  # The first-order variance is the one at q = 0, and its warning says so.
  w0 <- NULL
  first_order <- withCallingHandlers(
    var_semipartial_r(.46, 100, p = 2, R2_full = .21, under_null = FALSE),
    dmar_semipartial_r_rounded = function(cnd) {
      w0 <<- cnd
      invokeRestart("muffleWarning")
    })
  expect_equal(first_order$value, (1 - .46^2)^2 / 98, tolerance = 1e-14)
  expect_match(conditionMessage(w0), "(R2_full = r_sp^2).", fixed = TRUE)
  # The boundary of the rule: (|r| - .005)^2 - .005, with r = .46: .202025.
  expect_warning(var_semipartial_r(.46, 100, p = 2, R2_full = .2021),
                 class = "dmar_semipartial_r_rounded")
  expect_error(var_semipartial_r(.46, 100, p = 2, R2_full = .2019),
               class = "dmar_semipartial_r_incoherent")
  # Further apart than rounding allows: refused as before.
  expect_error(var_semipartial_r(.46, 100, p = 2, R2_full = .15),
               class = "dmar_semipartial_r_incoherent")
  # A negative r_sp, and an |r_sp| below .005, where the bound is 0.
  expect_warning(var_semipartial_r(-.46, 100, p = 2, R2_full = .21),
                 class = "dmar_semipartial_r_rounded")
  expect_warning(var_semipartial_r(.004, 100, p = 2, R2_full = 0),
                 class = "dmar_semipartial_r_rounded")
  # The planner and the sensitivity function apply the same rule.
  expect_warning(
    plan <- ss_aipe_semipartial_r(true_r_sp = .46, p = 2, width = .2,
                                  true_R2_full = .21),
    class = "dmar_semipartial_r_rounded")
  expect_equal(plan$value[plan$term == "necessary_N"],
               ss_aipe_semipartial_r(true_r_sp = .46, p = 2, width = .2,
                                     true_R2_full = .46^2)$value[1])
  expect_error(ss_aipe_semipartial_r(true_r_sp = .46, p = 2, width = .2,
                                     true_R2_full = .15),
               class = "dmar_semipartial_r_incoherent")
})

test_that("var_partial_r() names the first-order variance when N is too small for it", {
  expect_error(var_partial_r(r = .3, N = 3, p = 3),
               "Need N > p for the first-order variance of the partial correlation.",
               fixed = TRUE)
})
