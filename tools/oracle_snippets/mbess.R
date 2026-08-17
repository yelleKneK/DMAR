# Live oracle comparisons against MBESS, removed from the test suite when
# the anchors were pinned (MBESS 4.9.3, 2026-08-09). Each block is
# self-contained and run at release time when MBESS is installed locally.

## from tests/testthat/test-numerical-correctness.R
local({
  dmar  <- DMAR::ci_R2(R2 = 0.50, N = 100, p = 5, conf_level = 0.95,
                       random_predictors = FALSE)
  mbess <- MBESS::ci.R2(R2 = 0.50, N = 100, K = 5, conf.level = 0.95,
                        Random.Predictors = FALSE)
  stopifnot(
    isTRUE(all.equal(dmar$value[dmar$term == "lower_limit"],
                     mbess$Lower.Conf.Limit.R2, tolerance = 1e-6)),
    isTRUE(all.equal(dmar$value[dmar$term == "upper_limit"],
                     mbess$Upper.Conf.Limit.R2, tolerance = 1e-6))
  )
})

## from tests/testthat/test-numerical-correctness.R
local({
  dmar  <- DMAR::ci_R2(R2 = 0.50, N = 100, p = 5, conf_level = 0.95,
                       random_predictors = TRUE)
  mbess <- MBESS::ci.R2(R2 = 0.50, N = 100, K = 5, conf.level = 0.95,
                        Random.Predictors = TRUE)
  stopifnot(
    isTRUE(all.equal(dmar$value[dmar$term == "lower_limit"],
                     mbess$Lower.Conf.Limit.R2, tolerance = 1e-5)),
    isTRUE(all.equal(dmar$value[dmar$term == "upper_limit"],
                     mbess$Upper.Conf.Limit.R2, tolerance = 1e-5))
  )
})

## from tests/testthat/test-numerical-correctness.R
local({
  dmar  <- DMAR::ci_smd(ncp = 4.0, n_1 = 30, n_2 = 30, conf_level = 0.95)
  mbess <- MBESS::ci.smd(ncp = 4.0, n.1 = 30, n.2 = 30, conf.level = 0.95)
  stopifnot(
    isTRUE(all.equal(dmar$value[dmar$term == "lower_limit"],
                     mbess$Lower.Conf.Limit.smd, tolerance = 1e-6)),
    isTRUE(all.equal(dmar$value[dmar$term == "upper_limit"],
                     mbess$Upper.Conf.Limit.smd, tolerance = 1e-6))
  )
})

## from tests/testthat/test-numerical-correctness.R
local({
  dmar  <- DMAR::conf_limits_ncf(F_value = 6.0, df_1 = 3, df_2 = 50,
                                 conf_level = 0.95)
  mbess <- MBESS::conf.limits.ncf(F.value = 6.0, df.1 = 3, df.2 = 50,
                                  conf.level = 0.95)
  stopifnot(
    isTRUE(all.equal(dmar$value[dmar$term == "lower_limit"],
                     mbess$Lower.Limit, tolerance = 1e-5)),
    isTRUE(all.equal(dmar$value[dmar$term == "upper_limit"],
                     mbess$Upper.Limit, tolerance = 1e-5))
  )
})

## from tests/testthat/test-numerical-correctness.R
local({
  dmar  <- DMAR::conf_limits_nct(t_value = 2.83, df = 126, conf_level = 0.95)
  mbess <- MBESS::conf.limits.nct(t.value = 2.83, df = 126,
                                  conf.level = 0.95)
  stopifnot(
    isTRUE(all.equal(dmar$value[dmar$term == "lower_limit"],
                     mbess$Lower.Limit, tolerance = 1e-5)),
    isTRUE(all.equal(dmar$value[dmar$term == "upper_limit"],
                     mbess$Upper.Limit, tolerance = 1e-5))
  )
})

## from tests/testthat/test-numerical-correctness.R
local({
  dmar  <- suppressWarnings(
    DMAR::ss_aipe_R2(population_R2 = 0.50, conf_level = 0.95,
                     width = 0.10, p = 5)
  )
  mbess <- suppressMessages(MBESS::ss.aipe.R2(
    Population.R2 = 0.50, conf.level = 0.95, width = 0.10, K = 5
  ))
  stopifnot(isTRUE(abs(dmar$value - as.numeric(mbess)) <= 1))
})

## from tests/testthat/test-numerical-correctness.R
local({
  dmar  <- DMAR::ss_aipe_smd(delta = 0.5, conf_level = 0.95, width = 0.4)
  mbess <- MBESS::ss.aipe.smd(delta = 0.5, conf.level = 0.95, width = 0.4)
  n_dmar <- dmar$value[dmar$term == "necessary_n_per_group"]
  stopifnot(isTRUE(abs(n_dmar - as.numeric(mbess)) <= 1))
})

## from tests/testthat/test-cfa_k.R
local({
  S <- matrix(
    c(1.384, 1.484, 1.988, 2.429, 3.031,
      1.484, 2.756, 2.874, 3.588, 4.390,
      1.988, 2.874, 4.845, 4.894, 6.080,
      2.429, 3.588, 4.894, 6.951, 7.476,
      3.031, 4.390, 6.080, 7.476, 10.313),
    nrow = 5
  )
  dimnames(S) <- list(paste0("y", 1:5), paste0("y", 1:5))
  res_k <- DMAR::cfa_k(S = S, factors = list(f = paste0("y", 1:5)), N = 300)
  oracle <- MBESS::ci.reliability(S = unname(S), N = 300, type = "omega",
                                  interval.type = "ml")
  stopifnot(
    isTRUE(all.equal(res_k$estimate[res_k$term == "omega_f"], oracle$est,
                     tolerance = 1e-6)),
    isTRUE(all.equal(res_k$se[res_k$term == "omega_f"], oracle$se,
                     tolerance = 1e-6))
  )
})

## from tests/testthat/test-reliability_omega_categorical.R
local({
  set.seed(113)
  N <- 400
  J <- 5
  loadings <- rep(0.7, J)
  eta <- rnorm(N)
  latent <- outer(eta, loadings) +
            matrix(rnorm(N * J), N, J) %*% diag(sqrt(1 - loadings^2))
  items_categorical <- apply(latent, 2, function(x)
    as.integer(cut(x, breaks = c(-Inf, -1.5, -0.5, 0.5, 1.5, Inf),
                   labels = FALSE)))
  colnames(items_categorical) <- paste0("y", seq_len(J))
  r <- DMAR::reliability_omega_categorical(data = items_categorical,
                                           ci_method = "none")
  dmar <- r$value[r$term == "estimate"]
  # MBESS::ci.reliability() reaches lavaan through do.call("lavaan", ...) with
  # no namespace qualifier, so lavaan must be attached, not merely installed.
  # Attach it for this comparison and detach afterward so the search path is
  # left as it was found.
  attached <- "package:lavaan" %in% search()
  if (!attached) {
    library(lavaan)
    on.exit(detach("package:lavaan", unload = FALSE, character.only = TRUE),
            add = TRUE)
  }
  oracle <- MBESS::ci.reliability(data = items_categorical,
                                  type = "categorical",
                                  interval.type = "none")
  stopifnot(isTRUE(all.equal(dmar, oracle$est, tolerance = 1e-6)))
})

## from tests/testthat/test-reliability_omega.R
local({
  set.seed(113)
  J <- 6
  N <- 250
  loadings <- seq(0.4, 0.8, length.out = J)
  eta <- rnorm(N)
  errors <- matrix(rnorm(N * J), N, J) %*% diag(sqrt(1 - loadings^2))
  items_continuous <- sweep(matrix(rep(eta, J), N, J), 2, loadings, `*`) +
    errors
  colnames(items_continuous) <- paste0("y", seq_len(J))

  r <- DMAR::reliability_omega(data = items_continuous,
                               denominator = "observed", ci_method = "none")
  dmar <- r$value[r$term == "estimate"]
  oracle <- MBESS::ci.reliability(data = items_continuous,
                                  type = "hierarchical",
                                  interval.type = "none")
  stopifnot(isTRUE(all.equal(dmar, oracle$est, tolerance = 1e-6)))
})

## from tests/testthat/test-reliability_omega.R
local({
  set.seed(113)
  J <- 6
  N <- 250
  loadings <- seq(0.4, 0.8, length.out = J)
  eta <- rnorm(N)
  errors <- matrix(rnorm(N * J), N, J) %*% diag(sqrt(1 - loadings^2))
  items_continuous <- sweep(matrix(rep(eta, J), N, J), 2, loadings, `*`) +
    errors
  colnames(items_continuous) <- paste0("y", seq_len(J))

  r <- DMAR::reliability_omega(data = items_continuous,
                               denominator = "model_implied",
                               ci_method = "likelihood")
  lo <- r$value[r$term == "lower_limit"]
  hi <- r$value[r$term == "upper_limit"]
  oracle <- MBESS::ci.reliability(data = items_continuous, type = "omega",
                                  interval.type = "ll")
  stopifnot(
    isTRUE(all.equal(lo, oracle$ci.lower, tolerance = 5e-4)),
    isTRUE(all.equal(hi, oracle$ci.upper, tolerance = 5e-4))
  )
})

## from tests/testthat/test-reliability_alpha_model_implied.R
local({
  set.seed(113)
  J <- 6
  N <- 250
  loadings <- seq(0.4, 0.8, length.out = J)
  eta <- rnorm(N)
  errors <- matrix(rnorm(N * J), N, J) %*% diag(sqrt(1 - loadings^2))
  items_continuous <- sweep(matrix(rep(eta, J), N, J), 2, loadings, `*`) +
    errors
  colnames(items_continuous) <- paste0("y", seq_len(J))

  # Point estimate.
  r <- DMAR::reliability_alpha(estimator = "model_implied",
                               data = items_continuous, ci_method = "none")
  dmar <- r$value[r$term == "estimate"]
  oracle <- MBESS::ci.reliability(data = items_continuous,
                                  type = "alpha-analytic",
                                  interval.type = "none")
  stopifnot(isTRUE(all.equal(dmar, oracle$est, tolerance = 1e-6)))

  # ML Wald interval.
  r <- DMAR::reliability_alpha(estimator = "model_implied",
                               data = items_continuous, ci_method = "ml")
  oracle <- MBESS::ci.reliability(data = items_continuous,
                                  type = "alpha-analytic",
                                  interval.type = "ml")
  stopifnot(
    isTRUE(all.equal(r$value[r$term == "lower_limit"], oracle$ci.lower,
                     tolerance = 2e-3)),
    isTRUE(all.equal(r$value[r$term == "upper_limit"], oracle$ci.upper,
                     tolerance = 2e-3))
  )

  # Likelihood interval.
  r <- DMAR::reliability_alpha(estimator = "model_implied",
                               data = items_continuous,
                               ci_method = "likelihood")
  oracle <- MBESS::ci.reliability(data = items_continuous,
                                  type = "alpha-analytic",
                                  interval.type = "ll")
  stopifnot(
    isTRUE(all.equal(r$value[r$term == "lower_limit"], oracle$ci.lower,
                     tolerance = 5e-4)),
    isTRUE(all.equal(r$value[r$term == "upper_limit"], oracle$ci.upper,
                     tolerance = 5e-4))
  )
})

## from tests/testthat/test-reliability_alpha.R
local({
  set.seed(113)
  J <- 6
  N <- 200
  loadings <- seq(0.4, 0.8, length.out = J)
  eta <- rnorm(N)
  errors <- matrix(rnorm(N * J), N, J) %*% diag(sqrt(1 - loadings^2))
  items_continuous <- sweep(matrix(rep(eta, J), N, J), 2, loadings, `*`) +
    errors
  colnames(items_continuous) <- paste0("y", seq_len(J))

  r <- DMAR::reliability_alpha(data = items_continuous,
                               estimator = "model_implied",
                               ci_method = "likelihood")
  lo <- r$value[r$term == "lower_limit"]
  hi <- r$value[r$term == "upper_limit"]
  # MBESS::ci.reliability() reaches lavaan through do.call("lavaan", ...) with
  # no namespace qualifier, so lavaan must be attached, not merely installed.
  attached <- "package:lavaan" %in% search()
  if (!attached) {
    library(lavaan)
    on.exit(detach("package:lavaan", unload = FALSE, character.only = TRUE),
            add = TRUE)
  }
  oracle <- MBESS::ci.reliability(data = items_continuous, type = "alpha-cfa",
                                  interval.type = "ll")
  stopifnot(
    isTRUE(all.equal(lo, oracle$ci.lower, tolerance = 5e-4)),
    isTRUE(all.equal(hi, oracle$ci.upper, tolerance = 5e-4))
  )
})

## from tests/testthat/test-reliability_alpha.R
local({
  set.seed(113)
  J <- 5; n <- 300
  lam <- c(.8, .7, .6, .5, .4)
  eta <- rnorm(n)
  d <- as.data.frame(outer(eta, lam) +
         matrix(rnorm(n * J), n, J) %*% diag(sqrt(1 - lam^2)))
  names(d) <- paste0("y", seq_len(J))

  a <- DMAR::reliability_alpha(data = d, estimator = "analytic",
                               ci_method = "none")
  m <- DMAR::reliability_alpha(data = d, estimator = "model_implied",
                               ci_method = "none")
  est <- function(x) x$value[x$term == "estimate"]
  # MBESS::ci.reliability() reaches lavaan through do.call("lavaan", ...) with
  # no namespace qualifier, so lavaan must be attached, not merely installed.
  attached <- "package:lavaan" %in% search()
  if (!attached) {
    library(lavaan)
    on.exit(detach("package:lavaan", unload = FALSE, character.only = TRUE),
            add = TRUE)
  }
  stopifnot(
    isTRUE(all.equal(est(a),
                     MBESS::ci.reliability(data = d, type = "alpha",
                                           interval.type = "none")$est,
                     tolerance = 1e-6)),
    isTRUE(all.equal(est(m),
                     MBESS::ci.reliability(data = d, type = "alpha-cfa",
                                           interval.type = "none")$est,
                     tolerance = 1e-6))
  )
})

## from tests/testthat/test-var_ete.R
local({
  grid <- expand.grid(
    covariate_value = c("sample_mean", "sd", "fixed"),
    type = c("sample", "population"),
    stringsAsFactors = FALSE
  )
  mbess_cv <- c(sample_mean = "sample.mean", sd = "SD", fixed = "fixed")

  for (i in seq_len(nrow(grid))) {
    ours <- DMAR::var_ete(
      sigma2 = 150.5, sigma2_Z = 210.7, n_1 = 64, n_2 = 246,
      beta_1 = 0.968, beta_2 = 0.778, mu_Z = 100, fixed_value = 115,
      type = grid$type[i], covariate_value = grid$covariate_value[i])
    oracle <- MBESS::var.ete(
      sigma2 = 150.5, sigmaz2 = 210.7, n1 = 64, n2 = 246,
      beta1 = 0.968, beta2 = 0.778, muz = 100, c = 115,
      type = grid$type[i],
      covariate.value = mbess_cv[[grid$covariate_value[i]]])
    stopifnot(isTRUE(all.equal(ours$value, unname(oracle),
                               tolerance = 1e-12)))
  }
})
