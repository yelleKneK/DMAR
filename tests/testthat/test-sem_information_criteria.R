ic_value <- function(x, term) x$value[match(term, x$term)]
ic_terms <- c("AIC", "BIC", "SABIC", "HBIC", "IBIC", "SPBIC",
              "log_likelihood", "npar", "N")
hs_three <- "visual =~ x1 + x2 + x3
             textual =~ x4 + x5 + x6
             speed =~ x7 + x8 + x9"

# ---- summary input: closed forms, shape, and arguments ------------------

test_that("the criteria agree with their closed forms", {
  x <- sem_information_criteria(log_likelihood = -120, npar = 2, N = 200,
                                parameter_estimate = c(0.4, 1), prior_mean = c(0, 1),
                                information_observed = diag(c(100, 50)),
                                information_expected = diag(c(80, 40)),
                                measures = c("AIC", "BIC", "SABIC", "HBIC",
                                             "IBIC", "SPBIC"))
  expect_s3_class(x, "dmar_sem_information")
  expect_s3_class(x, "dmar_tbl")
  expect_identical(x$term, ic_terms)
  expect_equal(ic_value(x, "AIC"), 244)
  expect_equal(ic_value(x, "BIC"), 240 + 2 * log(200))
  expect_equal(ic_value(x, "SABIC"), 240 + 2 * log(202 / 24))
  expect_equal(ic_value(x, "HBIC"), 240 + 2 * log(200 / (2 * pi)))
  # IBIC = -2l + d log(N / 2 pi) + log|I_E / N| (Bollen et al., 2012, Eq. 20).
  expect_equal(ic_value(x, "IBIC"), 240 + 2 * log(200 / (2 * pi)) + log(3200 / 200^2))
  # Q = 0.4^2 * 100 = 16 > d = 2, so SPBIC = -2l + d (1 + log(Q / d)).
  expect_equal(ic_value(x, "SPBIC"), 240 + 2 * (1 + log(8)))
  expect_equal(ic_value(x, "log_likelihood"), -120)
  expect_equal(ic_value(x, "npar"), 2)
  expect_equal(ic_value(x, "N"), 200)
  z <- attr(x, "sem_information")
  expect_equal(z$quadratic_form, 16)
  expect_equal(z$prior_mean, c(0, 1))
  expect_false(z$boundary)
  expect_equal(z$logdet_total_expected_information, log(3200))
  expect_identical(z$expected_information, "supplied")
  # Per observation matrices are multiplied by N.
  y <- sem_information_criteria(log_likelihood = -120, npar = 2, N = 200,
                                parameter_estimate = c(0.4, 1), prior_mean = c(0, 1),
                                information_observed = diag(c(0.5, 0.25)),
                                information_expected = diag(c(0.4, 0.2)),
                                information_scale = "per_observation",
                                measures = c("IBIC", "SPBIC"))
  expect_equal(ic_value(y, c("IBIC", "SPBIC")), ic_value(x, c("IBIC", "SPBIC")))
})

test_that("every table has the same rows, with NA for criteria not requested", {
  x <- sem_information_criteria(log_likelihood = -120, npar = 5, N = 200)
  expect_identical(x$term, ic_terms)
  expect_true(all(is.na(ic_value(x, c("IBIC", "SPBIC")))))
  expect_false(anyNA(ic_value(x, c("AIC", "BIC", "SABIC", "HBIC"))))
  y <- sem_information_criteria(log_likelihood = -120, npar = 5, N = 200,
                                measures = "bic")
  expect_identical(y$term, ic_terms)
  expect_equal(ic_value(y, "BIC"), ic_value(x, "BIC"))
  expect_true(all(is.na(ic_value(y, c("AIC", "SABIC", "HBIC", "IBIC", "SPBIC")))))
  expect_identical(attr(x, "fixed_terms"), c(ic_terms[1:6], "log_likelihood"))
  expect_output(print(y), "Not requested: AIC, SABIC, HBIC, IBIC, SPBIC")
})

test_that("SPBIC handles centering, the boundary branch, and linear reparameterization", {
  spbic_of <- function(...) {
    ic_value(sem_information_criteria(log_likelihood = -10, npar = 2, N = 100,
                                      measures = "SPBIC", ...), "SPBIC")
  }
  # Q <= d: the point mass branch, -2l + Q (Bollen et al., 2012, Eq. 15).
  expect_warning(v <- spbic_of(parameter_estimate = c(0, 0), prior_mean = c(0, 0),
                               information_observed = diag(2)), "boundary branch")
  expect_equal(v, 20)
  expect_warning(v <- spbic_of(parameter_estimate = c(1, 1), prior_mean = c(0, 0),
                               information_observed = diag(2)), "boundary branch")
  expect_equal(v, 22)
  # The two branches meet at Q = d.
  expect_warning(at <- spbic_of(parameter_estimate = c(1, 1), prior_mean = c(0, 0),
                                information_observed = diag(2)), "boundary")
  above <- spbic_of(parameter_estimate = c(1, 1) * (1 + 1e-7), prior_mean = c(0, 0),
                    information_observed = diag(2))
  expect_equal(at, above, tolerance = 1e-6)
  x <- spbic_of(parameter_estimate = c(3, 4), prior_mean = c(0, 4), information_observed = diag(2))
  expect_equal(x, spbic_of(parameter_estimate = c(3, 4), zero_parameters = 1,
                           information_observed = diag(2)))
  expect_equal(x, spbic_of(parameter_estimate = c(a = 3, b = 4), zero_parameters = "a",
                           information_observed = diag(2)))
  # Transform the estimates, their prior mean, and the information together.
  A <- matrix(c(2, 0, 1, 3), 2)
  I_new <- t(solve(A)) %*% diag(2) %*% solve(A)
  expect_equal(x, spbic_of(parameter_estimate = as.vector(A %*% c(3, 4)),
                           prior_mean = as.vector(A %*% c(0, 4)),
                           information_observed = I_new))
})

test_that("Bollen et al. (2014) Table 6 is reproduced in all four rows", {
  # Bollen, Harden, Ray, and Zavisca (2014, Table 6, p. 18): five models of
  # Bollen and Curran (2004), N = 3912, each relative to the saturated model
  # (d = 20; 15 for SPBIC). The inputs are the printed T, d, the printed log
  # determinants of the expected information per observation, and the
  # printed quadratic forms; these are arithmetic reconstructions, not refits.
  N <- 3912
  Tml <- c(534.27, 412.93, 26.26, 62.44, 203.33)
  d <- c(14, 10, 17, 14, 13)
  d_spbic <- c(9, 3, 10, 7, 6)
  logdet <- c(102.61, 74.26, 130.18, 102.77, 94.19)
  Q <- c(1137006314, 180774515, 214642815, 251962808, 189536063)
  criteria <- function(ll, d, d_spbic, logdet, Q) {
    x <- sem_information_criteria(
      log_likelihood = ll, npar = d, N = N,
      measures = c("BIC", "HBIC", "IBIC", "SPBIC"),
      parameter_estimate = c(sqrt(Q), rep(0, d - 1)), prior_mean = rep(0, d),
      information_observed = diag(d) / N,
      information_expected = diag(c(exp(logdet), rep(1, d - 1))),
      information_scale = "per_observation", spbic_dimension = d_spbic)
    ic_value(x, c("BIC", "HBIC", "IBIC", "SPBIC"))
  }
  saturated <- criteria(0, 20, 15, 134.89, 195428707)
  relative <- vapply(seq_along(Tml), function(i) {
    criteria(-Tml[i] / 2, d[i], d_spbic[i], logdet[i], Q[i]) - saturated
  }, numeric(4))
  printed <- rbind(BIC = c(484.64, 330.21, 1.44, 12.81, 145.43),
                   HBIC = c(495.67, 348.59, 6.95, 23.83, 158.29),
                   IBIC = c(463.39, 287.95, 2.25, -8.29, 117.59),
                   SPBIC = c(450.42, 208.93, -55.66, -69.51, 52.20))
  # The printed values are rounded to two decimals from rounded inputs.
  expect_lt(max(abs(relative - printed)), 0.015)
  # The table marks the ALT model with rho free as best by SPBIC, but its
  # own SPBIC row is smallest for the ALT model with rho equal.
  expect_identical(unname(apply(relative, 1, which.min)), c(3L, 3L, 4L, 4L))
})

test_that("Bollen et al. (2012) Table 1 agrees except for its double count of log N in IBIC", {
  # Bollen, Ray, Zavisca, and Harden (2012, Table 1, p. 306): three
  # regression models, N = 50, with the printed log likelihoods, parameter
  # counts, quadratic forms, and log determinants of the whole sample
  # information.
  N <- 50
  ll <- c(-87.21, -125.37, -84.80)
  d <- c(2, 2, 3)
  Q <- c(270.74, 21.26, 296.67)
  logdet_total <- c(6.42, 3.92, 9.87)
  printed <- rbind(SPBIC = c(186.24, 257.47, 186.38),
                   IBIC = c(184.99, 258.82, 185.70),
                   HBIC = c(178.57, 254.90, 175.83),
                   BIC = c(182.24, 258.57, 181.34))
  ours <- function(i, scale) {
    x <- sem_information_criteria(
      log_likelihood = ll[i], npar = d[i], N = N,
      measures = c("SPBIC", "IBIC", "HBIC", "BIC"),
      parameter_estimate = c(sqrt(Q[i]), rep(0, d[i] - 1)), prior_mean = rep(0, d[i]),
      information_observed = diag(d[i]) / if (scale == "total") 1 else N,
      information_expected = diag(c(exp(logdet_total[i]), rep(1, d[i] - 1))),
      information_scale = scale)
    ic_value(x, c("SPBIC", "IBIC", "HBIC", "BIC"))
  }
  eq20 <- vapply(1:3, ours, numeric(4), scale = "total")
  # Two decimals from inputs rounded to two decimals: within 0.015.
  expect_lt(max(abs(eq20[-2, ] - printed[-2, ])), 0.015)
  # The table's IBIC adds the whole sample log determinant to d log(N / 2 pi),
  # which exceeds Eq. 20 by d log N.
  expect_lt(max(abs(eq20[2, ] + d * log(N) - printed["IBIC", ])), 0.015)
})

test_that("an input a requested criterion does not use is refused", {
  skip_on_cran()  # The extended validation grid runs in the full local suite.
  base <- list(log_likelihood = -120, npar = 2, N = 200)
  expect_error(do.call(sem_information_criteria, c(base, list(parameter_estimate = c(0.4, 1)))),
               "'parameter_estimate' is used only for SPBIC")
  expect_error(do.call(sem_information_criteria, c(base, list(prior_mean = c(0, 1),
                                                              zero_parameters = 1))),
               "'prior_mean', 'zero_parameters' are used only for SPBIC")
  expect_error(do.call(sem_information_criteria,
                       c(base, list(information_observed = diag(2)))), "only for SPBIC")
  expect_error(do.call(sem_information_criteria, c(base, list(spbic_dimension = 1))),
               "only for SPBIC")
  expect_error(do.call(sem_information_criteria,
                       c(base, list(information_expected = diag(2), measures = "SPBIC",
                                    parameter_estimate = c(1, 1), prior_mean = c(0, 0),
                                    information_observed = diag(2)))),
               "used only for IBIC")
  expect_error(do.call(sem_information_criteria, c(base, list(measures = "IBIC"))),
               "IBIC needs 'information_expected'")
  expect_error(do.call(sem_information_criteria,
                       c(base, list(measures = "SPBIC", parameter_estimate = c(1, 1),
                                    prior_mean = c(0, 0)))),
               "SPBIC needs 'information_observed'")
  expect_error(do.call(sem_information_criteria,
                       c(base, list(measures = "SPBIC", parameter_estimate = c(1, 1),
                                    prior_mean = c(0, 0), spbic_dimension = 3,
                                    information_observed = diag(2)))),
               "must not exceed npar")
  expect_error(sem_information_criteria(log_likelihood = -120, npar = 2),
               "summary input")
})

test_that("invalid inputs fail instead of fabricating a criterion", {
  skip_on_cran()
  base <- list(log_likelihood = -10, npar = 2, N = 100, measures = "SPBIC",
               parameter_estimate = c(1, 2), prior_mean = c(0, 0))
  expect_error(do.call(sem_information_criteria,
                       c(base, list(information_observed = diag(3)))), "2 by 2")
  expect_error(do.call(sem_information_criteria,
                       c(base, list(information_observed = diag(c(1, 0))))),
               "positive definite")
  expect_error(do.call(sem_information_criteria,
                       c(base, list(information_observed = matrix(c(1, 2, 0, 1), 2)))),
               "symmetric")
  expect_error(sem_information_criteria(log_likelihood = -10, npar = 2, N = 100,
                                        measures = "SPBIC", parameter_estimate = c(1, 2),
                                        information_observed = diag(2)),
               "prior_mean")
  expect_error(sem_information_criteria(log_likelihood = -10, npar = 2, N = 100,
                                        measures = "SPBIC", parameter_estimate = c(a = 1, b = 2),
                                        prior_mean = c(b = 0, a = 0),
                                        information_observed = diag(2)), "order")
  expect_error(sem_information_criteria(log_likelihood = Inf, npar = 2, N = 100), "finite")
  expect_error(sem_information_criteria(log_likelihood = -10, npar = 2.5, N = 100), "integer")
  expect_error(sem_information_criteria(log_likelihood = -10, npar = 2, N = 1), "greater")
  expect_error(sem_information_criteria(log_likelihood = -10, npar = 2, N = 100,
                                        measures = "AICc"), "distinct names")
  expect_error(sem_information_criteria(log_likelihood = -10, npar = 2, N = 100,
                                        measures = c("BIC", "bic")), "distinct names")
})

test_that("SABIC warns when its penalty is not positive", {
  expect_warning(x <- sem_information_criteria(log_likelihood = -10, npar = 2, N = 22),
                 "not positive")
  expect_equal(ic_value(x, "SABIC"), 20)
  expect_silent(sem_information_criteria(log_likelihood = -10, npar = 2, N = 22,
                                         measures = "BIC"))
  expect_silent(sem_information_criteria(log_likelihood = -10, npar = 2, N = 23))
})

# ---- fitted lavaan models -------------------------------------------------

test_that("a fitted model gives lavaan's AIC, BIC, and SABIC and Bollen's HBIC", {
  skip_on_cran()  # Additional model fits run in the full local suite.
  skip_if_not_installed("lavaan")
  f <- lavaan::cfa(hs_three, data = lavaan::HolzingerSwineford1939)
  fm <- lavaan::fitMeasures(f, c("aic", "bic", "bic2", "logl", "npar", "ntotal"))
  x <- sem_information_criteria(f)
  expect_equal(ic_value(x, c("AIC", "BIC", "SABIC")), unname(fm[c("aic", "bic", "bic2")]))
  expect_equal(ic_value(x, "HBIC"), -2 * fm[["logl"]] + 21 * log(301 / (2 * pi)))
  expect_equal(ic_value(x, c("log_likelihood", "npar", "N")),
               unname(fm[c("logl", "npar", "ntotal")]))
  z <- attr(x, "sem_information")
  expect_identical(z$likelihood, "the normal log likelihood of the complete data")
  expect_false(z$meanstructure)
  expect_length(z$conditioning, 0L)
  expect_output(print(x), "Log likelihood: the normal log likelihood")
  expect_error(sem_information_criteria(f, N = 301), "not both")
  expect_error(sem_information_criteria(f, information_scale = "total"),
               "remove 'information_scale'")
})

test_that("AIC, BIC, and SABIC agree with lavaan across the fits a user meets", {
  skip_on_cran()
  skip_if_not_installed("lavaan")
  d <- lavaan::HolzingerSwineford1939
  dm <- d
  dm$x1[1:15] <- NA
  dm$x5[16:30] <- NA
  fits <- list(
    mlr = lavaan::cfa(hs_three, data = d, estimator = "MLR"),
    two_groups = lavaan::cfa(hs_three, data = d, group = "school"),
    two_groups_equal = lavaan::cfa(hs_three, data = d, group = "school",
                                   group.equal = "loadings"),
    fiml = lavaan::cfa(hs_three, data = dm, missing = "ml"),
    fiml_mlr = lavaan::cfa(hs_three, data = dm, missing = "ml", estimator = "MLR"),
    listwise = lavaan::cfa(hs_three, data = dm),
    means = lavaan::cfa(hs_three, data = d, meanstructure = TRUE),
    fixed_x = lavaan::sem("visual =~ x1 + x2 + x3; visual ~ ageyr + grade",
                          data = d),
    labels = lavaan::cfa("f =~ x1 + a*x2 + a*x3", data = d),
    explicit = lavaan::cfa("f =~ x1 + a*x2 + b*x3; a == 2*b", data = d),
    efa = lavaan::cfa('efa("b")*f1 + efa("b")*f2 =~ x1 + x2 + x3 + x4 + x5 + x6',
                      data = d),
    moments = lavaan::cfa("f =~ x1 + x2 + x3",
                          sample.cov = stats::cov(d[, c("x1", "x2", "x3")]),
                          sample.nobs = 301))
  for (nm in names(fits)) {
    f <- fits[[nm]]
    fm <- lavaan::fitMeasures(f, c("aic", "bic", "bic2", "npar"))
    x <- sem_information_criteria(f)
    expect_equal(ic_value(x, c("AIC", "BIC", "SABIC")),
                 unname(fm[c("aic", "bic", "bic2")]), label = nm)
    expect_equal(ic_value(x, "npar"), unname(fm[["npar"]]), label = nm)
  }
  # The rotated block has 19 coefficients and 17 independent parameters.
  expect_equal(ic_value(sem_information_criteria(fits$efa), "npar"), 17)
  expect_match(attr(sem_information_criteria(fits$fiml), "sem_information")$likelihood,
               "FIML")
  expect_match(attr(sem_information_criteria(fits$listwise), "sem_information")$likelihood,
               "30 cases removed by listwise deletion")
  expect_identical(attr(sem_information_criteria(fits$fixed_x), "sem_information")$conditioning,
                   c("ageyr", "grade"))
})

test_that("automatic IBIC and SPBIC use lavaan's information at the right scale and prior", {
  skip_on_cran()  # Additional model fits run in the full local suite.
  skip_if_not_installed("lavaan")
  f <- lavaan::cfa(hs_three, data = lavaan::HolzingerSwineford1939,
                   information = "observed")
  fm <- lavaan::fitMeasures(f)
  x <- sem_information_criteria(f, measures = c("AIC", "IBIC", "SPBIC"))
  I_O <- lavaan::lavInspect(f, "information.observed") * 301
  expect_equal(I_O, solve(lavaan::vcov(f)), ignore_attr = TRUE, tolerance = 1e-8)
  z <- attr(x, "sem_information")
  theta <- lavaan::coef(f)
  pt <- lavaan::parTable(f)
  selected <- pt$free[pt$free > 0 & pt$op %in% c("=~", "~")]
  centered <- theta * 0
  centered[selected] <- theta[selected]
  expect_equal(z$quadratic_form, as.numeric(t(centered) %*% I_O %*% centered))
  expect_equal(z$prior_mean[-selected], theta[-selected])
  expect_equal(unname(z$prior_mean[selected]), rep(0, length(selected)))
  I_E <- lavaan::lavInspect(f, "information.expected") * 301
  expect_equal(ic_value(x, "IBIC"),
               -2 * fm[["logl"]] - length(theta) * log(2 * pi) +
                 as.numeric(determinant(I_E, logarithm = TRUE)$modulus))
  expect_identical(z$identification,
                   c(visual = "marker indicator", textual = "marker indicator",
                     speed = "marker indicator"))
  expect_match(z$observed_information, "negative Hessian")
  expect_match(z$expected_information, "fitted moments")
  expect_output(print(x), "IBIC and SPBIC depend on how the latent variables are identified")
  s <- lavaan::cfa(hs_three, data = lavaan::HolzingerSwineford1939, std.lv = TRUE)
  expect_identical(unname(attr(sem_information_criteria(s, measures = "IBIC"),
                               "sem_information")$identification),
                   rep("variance fixed", 3))
})

test_that("IBIC depends on the units and the identification; the others do not", {
  skip_on_cran()
  skip_if_not_installed("lavaan")
  d <- lavaan::HolzingerSwineford1939
  d10 <- d
  d10$x1 <- 10 * d10$x1
  all6 <- c("AIC", "BIC", "SABIC", "HBIC", "IBIC", "SPBIC")
  difference <- function(data, ...) {
    a <- sem_information_criteria(lavaan::cfa(hs_three, data = data, ...),
                                  measures = all6)
    b <- sem_information_criteria(lavaan::cfa(paste(hs_three, "\n x1 ~~ x3"),
                                              data = data, ...), measures = all6)
    stats::setNames(ic_value(b, all6) - ic_value(a, all6), all6)
  }
  original <- difference(d)
  rescaled <- difference(d10)
  # The extra residual covariance of x1 and x3 carries the units of x1 once,
  # so its IBIC difference moves by -2 log 10 and here changes sign.
  expect_equal(rescaled[["IBIC"]] - original[["IBIC"]], -2 * log(10), tolerance = 1e-3)
  expect_gt(original[["IBIC"]], 0)
  expect_lt(rescaled[["IBIC"]], 0)
  expect_equal(rescaled[1:4], original[1:4], tolerance = 1e-6)
  expect_equal(rescaled[["SPBIC"]], original[["SPBIC"]], tolerance = 1e-3)
  # The identification changes IBIC and SPBIC but not the likelihood criteria.
  marker <- sem_information_criteria(lavaan::cfa(hs_three, data = d), measures = all6)
  variance <- sem_information_criteria(lavaan::cfa(hs_three, data = d, std.lv = TRUE),
                                       measures = all6)
  expect_equal(ic_value(marker, all6[1:4]), ic_value(variance, all6[1:4]),
               tolerance = 1e-6)
  expect_gt(abs(ic_value(marker, "IBIC") - ic_value(variance, "IBIC")), 1)
  expect_gt(abs(ic_value(marker, "SPBIC") - ic_value(variance, "SPBIC")), 1)
})

test_that("fits whose log likelihood the criteria are not defined on are refused", {
  skip_on_cran()  # Additional model fits run in the full local suite.
  skip_if_not_installed("lavaan")
  d <- lavaan::HolzingerSwineford1939
  ordered <- d
  for (v in c("x1", "x2", "x3")) ordered[[v]] <- cut(ordered[[v]], 3, labels = FALSE)
  expect_error(sem_information_criteria(
    lavaan::cfa("f =~ x1 + x2 + x3", data = ordered, ordered = c("x1", "x2", "x3"))),
    "ordered categorical indicators")
  expect_error(sem_information_criteria(
    lavaan::cfa("f =~ x1 + x2 + x3", data = d, estimator = "GLS")),
    "estimator = \"GLS\"")
  expect_error(sem_information_criteria(
    lavaan::cfa("f =~ x1 + x2 + x3", data = d, likelihood = "wishart")),
    "likelihood = \"wishart\"")
  weighted <- d
  weighted$w <- rep(c(0.5, 1.5), length.out = nrow(d))
  fw <- lavaan::cfa("f =~ x1 + x2 + x3", data = weighted, sampling.weights = "w")
  expect_error(sem_information_criteria(fw), "sampling weights")
  expect_error(sem_information_criteria(fw, measures = "AIC"), "pseudo log likelihood")
  expect_error(sem_information_criteria(lm(x1 ~ x2, data = d)), "must be a lavaan model")
})

test_that("a two-level model gives AIC only", {
  skip_on_cran()
  skip_if_not_installed("lavaan")
  f <- lavaan::sem("level: 1
                      fw =~ y1 + y2 + y3
                    level: 2
                      fb =~ y1 + y2 + y3",
                   data = lavaan::Demo.twolevel, cluster = "cluster")
  expect_error(sem_information_criteria(f), "two-level model")
  x <- sem_information_criteria(f, measures = "AIC")
  expect_equal(ic_value(x, "AIC"), lavaan::fitMeasures(f, "aic")[[1]])
})

test_that("automatic IBIC and SPBIC are refused where lavaan's matrices are not theirs", {
  skip_on_cran()  # Additional model fits run in the full local suite.
  skip_if_not_installed("lavaan")
  d <- lavaan::HolzingerSwineford1939
  spbic <- function(f) sem_information_criteria(f, measures = "SPBIC")
  ibic <- function(f) sem_information_criteria(f, measures = "IBIC")
  expect_error(spbic(lavaan::cfa(hs_three, data = d, estimator = "MLR")),
               "Automatic SPBIC.*estimator = \"MLR\"")
  expect_error(spbic(lavaan::cfa("f =~ x1 + a*x2 + a*x3", data = d)),
               "Automatic SPBIC.*constraints")
  expect_error(spbic(lavaan::cfa("f =~ x1 + a*x2 + b*x3; a == 2*b", data = d)),
               "constraints")
  expect_error(spbic(lavaan::sem("x1 ~ x2 + x3", data = d)),
               "conditional on the exogenous observed variables")
  incomplete <- d
  incomplete$x1[1:10] <- NA
  expect_error(spbic(lavaan::cfa("f =~ x1 + x2 + x3", data = incomplete, missing = "ml")),
               "missing values")
  expect_error(spbic(lavaan::cfa("f =~ x1 + x2 + x3", data = d, group = "school")),
               "more than one group")
  efa <- lavaan::cfa('efa("b")*f1 + efa("b")*f2 =~ x1 + x2 + x3 + x4 + x5 + x6', data = d)
  expect_error(ibic(efa), "Automatic IBIC.*exploratory")
  expect_error(spbic(efa), "Automatic SPBIC.*exploratory")
  h1 <- lavaan::cfa(hs_three, data = d, observed.information = "h1")
  expect_error(spbic(h1), "observed.information = \"h1\"")
  expect_s3_class(ibic(h1), "dmar_sem_information")
  unstructured <- lavaan::cfa(hs_three, data = d, h1.information = "unstructured")
  expect_error(ibic(unstructured), "h1.information = \"unstructured\"")
  expect_s3_class(spbic(unstructured), "dmar_sem_information")
  # The likelihood criteria remain available for each of these fits.
  expect_s3_class(sem_information_criteria(efa), "dmar_sem_information")
})

test_that("two-stage criteria are Lai's (2021) and say which likelihood they use", {
  skip_on_cran()
  skip_if_not_installed("lavaan")
  d <- lavaan::HolzingerSwineford1939
  set.seed(113)
  for (v in c("x2", "x5", "x8")) d[[v]][sample(301, 40)] <- NA
  three <- lavaan::cfa(hs_three, data = d, missing = "two.stage")
  one <- lavaan::cfa("g =~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9", data = d,
                     missing = "two.stage")
  x3 <- sem_information_criteria(three)
  x1 <- sem_information_criteria(one)
  expect_match(attr(x3, "sem_information")$likelihood, "two-stage")
  # Lai (2021, Eq. 32): AIC_TS = n F_TS + 2q. lavaan's fmin is half the ML
  # discrepancy function, so n F_TS = 2 n fmin. lavaan's two-stage log
  # likelihood differs from -n F_TS / 2 by a constant shared by both models.
  lai <- function(f) 2 * 301 * lavaan::fitMeasures(f, "fmin")[[1]] +
    2 * lavaan::fitMeasures(f, "npar")[[1]]
  expect_equal(ic_value(x3, "AIC") - ic_value(x1, "AIC"), lai(three) - lai(one),
               tolerance = 1e-6)
  expect_error(sem_information_criteria(three, measures = "IBIC"), "missing values")
})
