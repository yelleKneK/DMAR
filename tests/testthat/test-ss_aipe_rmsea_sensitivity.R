## ss_aipe_rmsea_sensitivity() carries out a Monte Carlo sensitivity study for
## the RMSEA AIPE planner. It fits the proposed model with lavaan (the SEM
## backend used throughout DMAR); the population RMSEA is read from
## lavaan::fitMeasures(), which is computed reliably for the large-N population
## fit. (Earlier versions used the sem package, whose summary()$RMSEA returned
## NA for this population fit, so the function could not run.) The Monte Carlo
## body is skipped on CRAN to stay within the test-time budget.

test_that("ss_aipe_rmsea_sensitivity() is exported with the documented API", {
  expect_true(is.function(ss_aipe_rmsea_sensitivity))
  expect_true(all(c("width", "model", "true_cov_matrix", "specified_N", "conf_level", "G",
                    "filename") %in%
                    names(formals(ss_aipe_rmsea_sensitivity))))
  expect_false("Sigma" %in% names(formals(ss_aipe_rmsea_sensitivity)))
  # 0.90 since DMAR 1.1.0, the conventional level of the RMSEA interval and
  # the default of ss_aipe_rmsea() and ci_rmsea(); DMAR 1.0.x used 0.95.
  expect_equal(eval(formals(ss_aipe_rmsea_sensitivity)$conf_level), 0.90)
  expect_null(eval(formals(ss_aipe_rmsea_sensitivity)$specified_N))
  # Nothing is written unless the caller names a file; the old `save`
  # switch and its default path in the working directory are gone.
  expect_null(eval(formals(ss_aipe_rmsea_sensitivity)$filename))
  expect_false("save" %in% names(formals(ss_aipe_rmsea_sensitivity)))
})

test_that("ss_aipe_rmsea_sensitivity() recovers the population RMSEA and runs end to end", {
  skip_on_cran()
  skip_if_not_installed("lavaan")
  skip_if_not_installed("MASS")

  # True data generating model: two correlated factors (r = 0.5), loadings 0.7.
  # Build the implied population covariance Sigma = Lambda Phi Lambda' + Psi.
  Lambda <- matrix(0, 6, 2); Lambda[1:3, 1] <- 0.7; Lambda[4:6, 2] <- 0.7
  Phi   <- matrix(c(1, 0.5, 0.5, 1), 2, 2)
  Sigma <- Lambda %*% Phi %*% t(Lambda) + diag(1 - 0.7^2, 6)
  dimnames(Sigma) <- list(paste0("x", 1:6), paste0("x", 1:6))
  proposed <- "g =~ x1 + x2 + x3 + x4 + x5 + x6"

  # The population RMSEA from lavaan is finite (not NA, the sem failure mode)
  # and positive, since the one-factor model is misspecified for two-factor data.
  pop_rmsea <- unname(lavaan::fitMeasures(
    lavaan::sem(proposed, sample.cov = Sigma, sample.nobs = 1e6), "rmsea"))
  expect_true(is.finite(pop_rmsea))
  expect_gt(pop_rmsea, 0)

  set.seed(113)
  res <- ss_aipe_rmsea_sensitivity(width = 0.05, model = proposed, true_cov_matrix = Sigma, G = 10)
  expect_s3_class(res, "dmar_tbl")
  expect_true(all(c("mean_rmsea", "median_rmsea", "sd_rmsea",
                    "mean_ci_width", "median_ci_width", "sd_ci_width",
                    "pct_ci_less_width", "pct_ci_miss_low", "pct_ci_miss_high",
                    "total_type_I_error",
                    "total_N", "df", "true_rmsea", "width",
                    "conf_level") %in% res$term))
  expect_equal(res$value[res$term == "df"], 9)
  # The reported population RMSEA matches the standalone lavaan population fit.
  expect_equal(res$value[res$term == "true_rmsea"], pop_rmsea, tolerance = 1e-6)
  expect_false("suc_rep" %in% res$term)
  expect_gt(attr(res, "simulation")$returned, 0)
})

test_that("ss_aipe_rmsea_sensitivity() requires named Sigma", {
  skip_if_not_installed("lavaan")
  skip_if_not_installed("MASS")
  expect_error(
    ss_aipe_rmsea_sensitivity(width = 0.05, model = "g =~ x1 + x2 + x3",
                              true_cov_matrix = diag(3), G = 2),
    "row and column names"
  )
})

test_that("ss_aipe_rmsea_sensitivity() writes per-replication results only to a named file", {
  skip_on_cran()
  skip_if_not_installed("lavaan")
  skip_if_not_installed("MASS")

  Lambda <- matrix(0, 6, 2); Lambda[1:3, 1] <- 0.7; Lambda[4:6, 2] <- 0.7
  Phi   <- matrix(c(1, 0.5, 0.5, 1), 2, 2)
  Sigma <- Lambda %*% Phi %*% t(Lambda) + diag(1 - 0.7^2, 6)
  dimnames(Sigma) <- list(paste0("x", 1:6), paste0("x", 1:6))
  proposed <- "g =~ x1 + x2 + x3 + x4 + x5 + x6"

  # A malformed filename is refused before any replication runs.
  expect_error(
    ss_aipe_rmsea_sensitivity(width = 0.05, model = proposed, true_cov_matrix = Sigma,
                              specified_N = 200, G = 3, filename = 42),
    "'filename' must be NULL or a single character string", fixed = TRUE)

  out_file <- tempfile(fileext = ".csv")
  on.exit(unlink(out_file), add = TRUE)
  set.seed(113)
  res <- ss_aipe_rmsea_sensitivity(width = 0.05, model = proposed,
                                   true_cov_matrix = Sigma, specified_N = 200, G = 3,
                                   filename = out_file)
  written <- utils::read.csv(out_file)
  # The limits are ci_lower and ci_upper, as in every other file and table
  # (DMAR 1.0.x wrote them as ci_low and ci_up).
  expect_named(written, c("iteration", "rmsea_hat", "ci_lower", "ci_upper", "width"))
  expect_equal(nrow(written), attr(res, "simulation")$returned)
  expect_equal(written$width, written$ci_upper - written$ci_lower)

  # A second run appends to the existing file and says so.
  set.seed(113)
  expect_message(
    ss_aipe_rmsea_sensitivity(width = 0.05, model = proposed, true_cov_matrix = Sigma,
                              specified_N = 200, G = 3, filename = out_file),
    "already exists", fixed = TRUE)
  expect_equal(nrow(utils::read.csv(out_file)), 2L * nrow(written))
})

# The population and the misspecified model of the tests below.
.sens_rmsea_setup <- function() {
  Lambda <- matrix(0, 6, 2); Lambda[1:3, 1] <- 0.7; Lambda[4:6, 2] <- 0.7
  Phi   <- matrix(c(1, 0.5, 0.5, 1), 2, 2)
  Sigma <- Lambda %*% Phi %*% t(Lambda) + diag(1 - 0.7^2, 6)
  dimnames(Sigma) <- list(paste0("x", 1:6), paste0("x", 1:6))
  list(Sigma = Sigma, model = "g =~ x1 + x2 + x3 + x4 + x5 + x6")
}

# lavaan's own fit measures of the samples a study draws: the population
# fit draws nothing, so after the same seed replication g's sample is the
# g-th draw, whether or not that replication converged.
.sens_rmsea_lavaan <- function(model, Sigma, N, G, seed, level = 0.90, ...) {
  set.seed(seed)
  out <- matrix(NA_real_, G, 4L,
                dimnames = list(NULL, c("chisq", "rmsea", "lower", "upper")))
  for (g in seq_len(G)) {
    d <- MASS::mvrnorm(n = N, mu = rep(0, nrow(Sigma)), Sigma = Sigma)
    colnames(d) <- rownames(Sigma)
    fit <- try(suppressWarnings(lavaan::sem(model, data = as.data.frame(d),
                                            ...)), silent = TRUE)
    if (inherits(fit, "try-error") ||
        !isTRUE(lavaan::lavInspect(fit, "converged"))) next
    out[g, ] <- unname(lavaan::fitMeasures(
      fit, c("chisq", "rmsea", "rmsea.ci.lower", "rmsea.ci.upper"),
      fm.args = list(rmsea.ci.level = level)))
  }
  as.data.frame(out)
}

test_that("each replication's interval is lavaan's, formed with lavaan's multiplier", {
  skip_on_cran()
  skip_if_not_installed("lavaan")
  skip_if_not_installed("MASS")
  s <- .sens_rmsea_setup()
  out_file <- tempfile(fileext = ".csv")
  on.exit(unlink(out_file), add = TRUE)
  # Each replication's interval is the one lavaan reports for its sample,
  # at the same level. Through DMAR 1.0.x the interval was formed from
  # lavaan's RMSEA, computed with N, under the multiplier N - 1, and
  # matched neither.
  set.seed(113)
  res <- ss_aipe_rmsea_sensitivity(width = 0.05, model = s$model,
                                   true_cov_matrix = s$Sigma, specified_N = 200, G = 2,
                                   filename = out_file)
  written <- utils::read.csv(out_file)
  lav <- .sens_rmsea_lavaan(s$model, s$Sigma, N = 200, G = 2, seed = 113)
  lav <- lav[written$iteration, ]
  expect_equal(written$rmsea_hat, lav$rmsea, tolerance = 1e-10)
  expect_equal(written$ci_lower, lav$lower, tolerance = 1e-6)
  expect_equal(written$ci_upper, lav$upper, tolerance = 1e-6)
  for (i in seq_len(nrow(written))) {
    ci_old <- ci_rmsea(written$rmsea_hat[i], df = 9, N = 200,
                       conf_level = 0.90)
    expect_gt(abs(ci_old$value[3] - written$ci_upper[i]), 1e-5)
  }
  # The simulated sample is the planner's necessary_N under either
  # likelihood: the size a study would collect.
  pop <- lavaan::sem(s$model, sample.cov = s$Sigma, sample.nobs = 1e6)
  planned <- suppressMessages(ss_aipe_rmsea(
    true_rmsea = unname(lavaan::fitMeasures(pop, "rmsea")), df = 9, width = 0.05,
    conf_level = 0.90))$value
  set.seed(113)
  default_plan <- ss_aipe_rmsea_sensitivity(width = 0.05, model = s$model,
                                            true_cov_matrix = s$Sigma, G = 1)
  expect_equal(default_plan$value[default_plan$term == "total_N"], planned)
  set.seed(113)
  wishart_plan <- suppressMessages(ss_aipe_rmsea_sensitivity(
    width = 0.05, model = s$model, true_cov_matrix = s$Sigma, G = 1,
    likelihood = "wishart"))
  expect_equal(wishart_plan$value[wishart_plan$term == "total_N"], planned)
  # Under likelihood = "wishart", passed through to lavaan, the statistic
  # and the RMSEA are formed with N - 1, and so is the interval.
  set.seed(113)
  suppressMessages(
    ss_aipe_rmsea_sensitivity(width = 0.05, model = s$model, true_cov_matrix = s$Sigma,
                              specified_N = 200, G = 1, filename = out_file,
                              likelihood = "wishart"))
  wishart <- utils::read.csv(out_file)
  last <- wishart[nrow(wishart), ]
  lav_w <- .sens_rmsea_lavaan(s$model, s$Sigma, N = 200, G = 1, seed = 113,
                              likelihood = "wishart")
  expect_equal(c(last$ci_lower, last$ci_upper), c(lav_w$lower, lav_w$upper),
               tolerance = 1e-6)
  ci_w <- ci_rmsea(last$rmsea_hat, df = 9, N = 200, conf_level = 0.90,
                   sample_size_convention = "N-1")
  expect_equal(c(last$ci_lower, last$ci_upper), ci_w$value[c(1, 3)],
               tolerance = 1e-6)
})

test_that("a replication whose statistic is at or below df gets lavaan's interval for that statistic", {
  skip_on_cran()
  skip_if_not_installed("lavaan")
  skip_if_not_installed("MASS")
  # The two-factor population fit by its own model (8 df): a correctly
  # specified model, whose statistics fall below df in about half the
  # replications. lavaan reports an RMSEA of 0 for each of them, from
  # which the statistic could be rebuilt only as df, the interval through
  # DMAR 1.0.x; the interval is formed from the statistic itself, as
  # lavaan forms it. At conf_level = 0.80 the upper limit is 0 for a
  # statistic at or below the 10th percentile of the central chi square on
  # 8 df (3.49), which two of these twelve samples reach, and between that
  # and df it is positive and below the limit rebuilt from an RMSEA of 0.
  s <- .sens_rmsea_setup()
  two <- "f1 =~ x1 + x2 + x3\n f2 =~ x4 + x5 + x6"
  G <- 12
  out_file <- tempfile(fileext = ".csv")
  on.exit(unlink(out_file), add = TRUE)
  set.seed(113)
  res <- ss_aipe_rmsea_sensitivity(width = 0.05, model = two,
                                   true_cov_matrix = s$Sigma, specified_N = 200, G = G,
                                   conf_level = 0.80, filename = out_file)
  written <- utils::read.csv(out_file)
  lav <- .sens_rmsea_lavaan(two, s$Sigma, N = 200, G = G, seed = 113,
                            level = 0.80)
  lav <- lav[written$iteration, ]
  below <- lav$chisq <= 8
  at_floor <- lav$chisq <= stats::qchisq(0.10, 8)
  expect_true(any(below & !at_floor))
  expect_true(any(at_floor))
  expect_true(all(written$rmsea_hat[below] == 0))
  expect_equal(written$ci_lower, lav$lower, tolerance = 1e-6)
  expect_equal(written$ci_upper, lav$upper, tolerance = 1e-6)
  expect_true(all(written$ci_upper[at_floor] == 0))
  # The limit rebuilt from an RMSEA of 0 is above every one of them.
  old_up <- suppressMessages(ci_rmsea(0, df = 8, N = 200, conf_level = 0.80,
                                      sample_size_convention = "N"))$value[3]
  expect_true(all(written$ci_upper[below] < old_up - 1e-4))
  expect_true(all(written$ci_upper[below & !at_floor] > 0))
  v <- function(term) res$value[res$term == term]
  expect_equal(v("true_rmsea"), 0)
  expect_equal(v("mean_ci_width"), mean(lav$upper - lav$lower),
               tolerance = 1e-6)
  expect_equal(v("pct_ci_less_width"), mean(written$width < 0.05))
})

test_that("a study at the default level gives no warning and reports the level", {
  skip_on_cran()
  skip_if_not_installed("lavaan")
  skip_if_not_installed("MASS")
  s <- .sens_rmsea_setup()
  set.seed(113)
  expect_no_warning(ss_aipe_rmsea_sensitivity(
    width = 0.05, model = s$model, true_cov_matrix = s$Sigma, specified_N = 200, G = 1,
    conf_level = 0.90))
  set.seed(113)
  expect_no_warning(res <- ss_aipe_rmsea_sensitivity(
    width = 0.05, model = s$model, true_cov_matrix = s$Sigma, specified_N = 200, G = 1))
  expect_equal(res$value[res$term == "conf_level"], 0.90)
})
