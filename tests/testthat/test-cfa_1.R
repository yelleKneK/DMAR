cfa_example_cov <- function() {
  matrix(
    c(1.384, 1.484, 1.988, 2.429, 3.031,
      1.484, 2.756, 2.874, 3.588, 4.390,
      1.988, 2.874, 4.845, 4.894, 6.080,
      2.429, 3.588, 4.894, 6.951, 7.476,
      3.031, 4.390, 6.080, 7.476, 10.313),
    nrow = 5
  )
}

cfa_example_data <- function(N = 300, loadings = c(.5, .6, .65, .7, .8),
                             seed = 113) {
  set.seed(seed)
  J <- length(loadings)
  eta <- rnorm(N)
  X <- sweep(matrix(rep(eta, J), N, J), 2, loadings, `*`) +
    matrix(rnorm(N * J), N, J) %*% diag(sqrt(1 - loadings^2))
  colnames(X) <- paste0("y", seq_len(J))
  X
}

test_that("cfa_1() is exactly the one factor cfa_k() call it wraps", {
  skip_if_not_installed("lavaan")
  X <- cfa_example_data()
  direct <- cfa_k(X, factors = list(f1 = colnames(X)))
  wrapped <- cfa_1(X)
  expect_identical(as.data.frame(wrapped), as.data.frame(direct))
  expect_s3_class(wrapped, "dmar_cfa_k")
})

test_that("cfa_1() default output is the cfa_k() verbose table", {
  skip_if_not_installed("lavaan")

  res <- cfa_1(S = cfa_example_cov(), N = 300)

  expect_s3_class(res, "data.frame")
  expect_named(res, c("syntax", "term", "estimate", "se", "z_value",
                      "p_value", "ci_lower", "ci_upper"))
  expect_true(all(c("lambda_f1_1", "lambda_f1_5", "psi_f1_1",
                    "psi_f1_5", "phi_f1", "omega_f1") %in% res$term))

  fit_rows <- c("chi_square", "df", "p_chi_square", "cfi", "tli", "nnfi",
                "rmsea", "rmsea_ci_lower", "rmsea_ci_upper", "AIC", "BIC",
                "H0", "H1")
  expect_true(all(fit_rows %in% res$term))
  expect_lte(res$estimate[res$term == "rmsea_ci_lower"],
             res$estimate[res$term == "rmsea_ci_upper"])
  expect_identical(attr(res, "fixed_terms"), c("AIC", "BIC", "H0", "H1"))
})

test_that("cfa_1() requires exactly one of data and S", {
  skip_if_not_installed("lavaan")
  X <- cfa_example_data(N = 50)
  expect_error(cfa_1(), "exactly one of 'data'")
  expect_error(cfa_1(X, S = cov(X), N = 50), "exactly one of 'data'")
  expect_error(cfa_1(S = cfa_example_cov()), "'N' must be provided")
  expect_error(cfa_1(S = as.data.frame(cfa_example_cov()), N = 300),
               "symmetric covariance matrix")
})

test_that("cfa_1() auto-names unnamed input and honors 'items'", {
  skip_if_not_installed("lavaan")
  X <- cfa_example_data()

  # Unnamed covariance matrix: auto-named y1..y5, same fit as named.
  named   <- cfa_1(S = cov(X), N = 300)
  unnamed <- cfa_1(S = unname(cov(X)), N = 300)
  expect_identical(as.data.frame(named), as.data.frame(unnamed))

  # Unnamed raw matrix.
  X2 <- X; colnames(X2) <- NULL
  expect_identical(as.data.frame(cfa_1(X2)), as.data.frame(cfa_1(X)))

  # 'items' selects a subset; a non-item column is simply not modeled.
  with_extra <- cbind(as.data.frame(X), id = seq_len(nrow(X)))
  sub <- cfa_1(with_extra, items = paste0("y", 1:4))
  expect_true("lambda_f1_4" %in% sub$term)
  expect_false("lambda_f1_5" %in% sub$term)
  expect_false(any(grepl("y5|id", sub$syntax)))

  expect_error(cfa_1(X, items = c("y1", "y2")), "three or more")
})

test_that("cfa_1() omega agrees with reliability_omega() and the closed form", {
  skip_if_not_installed("lavaan")
  X <- cfa_example_data(N = 200, loadings = seq(0.4, 0.8, length.out = 6))

  res <- cfa_1(X)
  om <- res$estimate[res$term == "omega_f1"]

  rl <- reliability_omega(data = X, ci_method = "none",
                          denominator = "model_implied")
  expect_equal(om, rl$value[rl$term == "estimate"], tolerance = 1e-8)

  # Hand computation from an independently specified lavaan fit.
  q <- ncol(X); vn <- colnames(X)
  load_line <- paste(paste0("a", 1:q, "*", vn), collapse = " + ")
  err_line  <- paste(paste0(vn, " ~~ b", 1:q, "*", vn), collapse = "\n")
  model <- paste0("f1 =~ NA*", vn[1], " + ", load_line,
                  "\nf1 ~~ 1*f1\n", err_line)
  fit <- lavaan::cfa(model, data = as.data.frame(X), se = "none")
  pe <- lavaan::parameterEstimates(fit)
  lams <- pe[pe$op == "=~", "est"]
  psis <- pe[pe$op == "~~" & pe$lhs == pe$rhs & pe$lhs %in% vn, "est"]
  manual <- (sum(lams))^2 / ((sum(lams))^2 + sum(psis))
  expect_equal(om, manual, tolerance = 1e-6)
})

test_that("cfa_1(estimator='MLR') runs on raw data and produces robust SEs", {
  skip_if_not_installed("lavaan")
  X <- cfa_example_data(N = 250)

  ml  <- cfa_1(X, estimator = "ML",  se = "standard")
  mlr <- cfa_1(X, estimator = "MLR", se = "robust.sem")

  expect_equal(ml$estimate[ml$term == "omega_f1"],
               mlr$estimate[mlr$term == "omega_f1"], tolerance = 1e-6)
  expect_true(is.finite(ml$se[ml$term == "omega_f1"]))
  expect_true(is.finite(mlr$se[mlr$term == "omega_f1"]))
})

test_that("cfa_1(missing='ml') handles missing values via FIML on raw data", {
  skip_if_not_installed("lavaan")
  X <- cfa_example_data(N = 200, loadings = rep(0.7, 5))
  missing_idx <- sample.int(length(X), size = round(0.1 * length(X)))
  X[missing_idx] <- NA

  res <- cfa_1(as.data.frame(X), missing = "ml")
  om <- res$estimate[res$term == "omega_f1"]
  expect_true(is.finite(om))
  expect_gt(om, 0)
  expect_lt(om, 1)
})

test_that("cfa_1() rejects bad estimator, missing, and ordered", {
  skip_if_not_installed("lavaan")
  expect_error(cfa_1(S = cfa_example_cov(), N = 300, estimator = 1),
               regexp = "'estimator' must be")
  expect_error(cfa_1(S = cfa_example_cov(), N = 300, missing = TRUE),
               regexp = "'missing' must be")
  expect_error(cfa_1(S = cfa_example_cov(), N = 300, estimator = "DWLS"),
               regexp = "'estimator' must be")
  expect_error(cfa_1(S = cfa_example_cov(), N = 300,
                     missing = "mean.imputation"),
               regexp = "'missing' must be")
  expect_error(cfa_1(cfa_example_data(), ordered = TRUE),
               regexp = "ordered categorical")
})

test_that("cfa_1(output='fit') returns the raw lavaan fit object", {
  skip_if_not_installed("lavaan")
  fit <- cfa_1(S = cfa_example_cov(), N = 300, output = "fit")
  expect_s4_class(fit, "lavaan")
  expect_true(lavaan::lavInspect(fit, "converged"))
})

test_that("missing = 'ml' estimates the mean structure with free intercepts", {
  skip_if_not_installed("lavaan")
  # Regression test: the bare lavaan::lavaan() interface leaves
  # int.ov.free = FALSE, so FIML fixed every item intercept at zero and
  # the loadings absorbed the item means. Items are shifted well away
  # from zero here so the defect could not hide.
  set.seed(113)
  N <- 300; J <- 4
  lam <- c(.8, .7, .6, .7)
  eta <- rnorm(N)
  X <- outer(eta, lam) +
       matrix(rnorm(N * J), N, J) %*% diag(sqrt(1 - lam^2))
  X <- sweep(X, 2, c(5, 10, 3, 7), `+`)
  colnames(X) <- paste0("y", seq_len(J))
  X[sample(N, 90), 2] <- NA
  X[sample(N, 60), 4] <- NA
  d <- as.data.frame(X)

  fit_dmar <- cfa_1(d, estimator = "ML", missing = "ml", output = "fit")
  ref_model <- paste0("f1 =~ NA*y1 + ",
                      paste0("l", seq_len(J), "*y", seq_len(J),
                             collapse = " + "),
                      "\nf1 ~~ 1*f1")
  fit_ref <- lavaan::cfa(ref_model, data = d, missing = "ml")

  expect_equal(unname(lavaan::coef(fit_dmar)[seq_len(J)]),
               unname(lavaan::coef(fit_ref)[seq_len(J)]),
               tolerance = 1e-8)
  expect_equal(as.numeric(lavaan::logLik(fit_dmar)),
               as.numeric(lavaan::logLik(fit_ref)), tolerance = 1e-8)

  # The intercepts are free and recover the item means.
  pe <- lavaan::parameterEstimates(fit_dmar)
  ints <- pe[pe$op == "~1" & pe$lhs %in% paste0("y", seq_len(J)), ]
  expect_true(all(ints$se > 0))
  expect_equal(ints$est, unname(colMeans(d, na.rm = TRUE)),
               tolerance = 0.2)
})
