sem_value <- function(x, term) x$value[match(term, x$term)]

residual_terms <- c("srmr", "srmr_se_exact_fit", "p_exact_fit",
                    "unbiased_srmr", "unbiased_srmr_se",
                    "unbiased_srmr_ci_lower", "unbiased_srmr_ci_upper",
                    "close_fit_bound", "p_close_fit", "srmr_mean",
                    "largest_absolute_residual", "N")

hs_model <- 'visual =~ x1 + x2 + x3
             textual =~ x4 + x5 + x6
             speed =~ x7 + x8 + x9'

# An independent computation of every row from Maydeu-Olivares (2017) as
# Maydeu-Olivares, Shi, and Rosseel (2018, Eqs. 16, 20 to 22) give it. The
# standardized residuals e = G vech(S - Sigma) have asymptotic covariance
# Xi = J Gamma J' / N, with J = G (I - Delta (Delta' W Delta)^{-1} Delta' W),
# W the normal theory weight D'(Sigma^{-1} x Sigma^{-1})D / 2 at the model
# implied Sigma (expected information), G the fixed standardizer, and Gamma
# the asymptotic covariance of vech(S): normal theory, or from the fourth
# order moments of the data.
srmr_oracle <- function(fit, X, gamma = c("normal", "adf"), level = 0.95, c0 = 0.06) {
  gamma <- match.arg(gamma)
  v <- colnames(X)
  N <- nrow(X)
  Xc <- scale(X, scale = FALSE)
  S <- crossprod(Xc) / N
  idx <- which(lower.tri(S, diag = TRUE), arr.ind = TRUE)
  vech <- function(A) A[lower.tri(A, diag = TRUE)]
  t <- nrow(idx)
  Sig <- lavaan::lavInspect(fit, "implied")$cov[v, v]
  Si <- solve(Sig)
  W <- G0 <- matrix(0, t, t)
  for (a in seq_len(t)) for (b in seq_len(t)) {
    i <- idx[a, 1]; j <- idx[a, 2]; k <- idx[b, 1]; l <- idx[b, 2]
    W[a, b] <- (Si[i, k] * Si[j, l] + Si[i, l] * Si[j, k]) / 2 *
      (if (i != j) 2 else 1) * (if (k != l) 2 else 1)
    G0[a, b] <- S[i, k] * S[j, l] + S[i, l] * S[j, k]
  }
  Gam <- if (gamma == "normal") G0 else {
    Dm <- t(apply(Xc, 1, function(x) vech(tcrossprod(x))))
    crossprod(scale(Dm, scale = FALSE)) / N
  }
  Del <- lavaan::lavInspect(fit, "delta")
  Q <- diag(t) - Del %*% solve(t(Del) %*% W %*% Del, t(Del) %*% W)
  sds <- sqrt(diag(S))
  gg <- 1 / (sds[idx[, 1]] * sds[idx[, 2]])
  J <- gg * Q
  e <- gg * vech(S - Sig)
  Xi <- J %*% Gam %*% t(J) / N
  tr1 <- sum(diag(Xi)); tr2 <- sum(Xi * Xi); ee <- sum(e^2)
  eXe <- drop(t(e) %*% Xi %*% e)
  k <- 1 - (tr2 + 2 * eXe) / (4 * ee^2)
  u <- sqrt(max(ee - tr1, 0) / t) / k
  se_u <- sqrt((tr2 + 2 * eXe) / (2 * t * ee) / k^2)
  srmr <- sqrt(ee / t)
  se0 <- sqrt(tr2 / (2 * t * tr1))
  e0 <- sqrt(tr1 / t) * (4 * tr1^2 - tr2) / (4 * tr1^2)
  z <- stats::qnorm(1 - (1 - level) / 2)
  c(srmr = srmr, srmr_se_exact_fit = se0,
    p_exact_fit = stats::pnorm(max(srmr - e0, 0) / se0, lower.tail = FALSE),
    unbiased_srmr = u, unbiased_srmr_se = se_u,
    unbiased_srmr_ci_lower = u - z * se_u, unbiased_srmr_ci_upper = u + z * se_u,
    close_fit_bound = c0,
    p_close_fit = stats::pnorm((u - c0) / se_u, lower.tail = FALSE))
}

test_that("every row equals an independent computation, normal theory and ADF", {
  skip_if_not_installed("lavaan")
  d <- lavaan::HolzingerSwineford1939
  X <- as.matrix(d[, paste0("x", 1:9)])
  f <- lavaan::cfa(hs_model, data = d, std.lv = TRUE)
  x <- sem_residual_fit(f, close_fit = 0.06)
  oracle <- srmr_oracle(f, X, "normal", level = 0.90)
  expect_equal(sem_value(x, names(oracle)), unname(oracle), tolerance = 1e-8)
  fm <- lavaan::cfa(hs_model, data = d, std.lv = TRUE, estimator = "MLM")
  y <- sem_residual_fit(fm, conf_level = 0.95, close_fit = 0.06)
  oracle <- srmr_oracle(fm, X, "adf", level = 0.95)
  expect_equal(sem_value(y, names(oracle)), unname(oracle), tolerance = 1e-8)
  # The sample SRMR is also lavaan's covariance SRMR and sem_fit_measures()'.
  expect_equal(sem_value(x, "srmr"),
               unname(lavaan::fitMeasures(f, "srmr_bentler_nomean")))
  expect_equal(sem_value(x, "srmr"), sem_value(sem_fit_measures(f), "srmr"))
})

test_that("the residual covariance follows the fit, and the attribute says which", {
  skip_if_not_installed("lavaan")
  d <- lavaan::HolzingerSwineford1939
  f <- lavaan::cfa(hs_model, data = d, std.lv = TRUE)
  mlm <- lavaan::cfa(hs_model, data = d, std.lv = TRUE, estimator = "MLM")
  mlr <- lavaan::cfa(hs_model, data = d, std.lv = TRUE, estimator = "MLR")
  mlr_e <- lavaan::cfa(hs_model, data = d, std.lv = TRUE, estimator = "MLR",
                       information = "expected")
  expect_equal(attr(sem_residual_fit(f), "sem_fit")$residual_acov, "normal theory")
  expect_match(attr(sem_residual_fit(mlm), "sem_fit")$residual_acov,
               "fourth order moments")
  r <- sem_residual_fit(mlr)
  expect_match(attr(r, "sem_fit")$residual_acov, "sandwich")
  expect_equal(attr(r, "sem_fit")$information, "observed")
  # MLR with expected information reproduces MLM; observed information
  # differs slightly.
  expect_equal(sem_residual_fit(mlr_e)$value, sem_residual_fit(mlm)$value,
               tolerance = 1e-8)
  expect_false(isTRUE(all.equal(r$value, sem_residual_fit(mlm)$value)))
})

test_that("the rows are fixed, with the close fit rows NA unless a bound is given", {
  skip_if_not_installed("lavaan")
  f <- lavaan::cfa(hs_model, data = lavaan::HolzingerSwineford1939, std.lv = TRUE)
  x <- sem_residual_fit(f)
  y <- sem_residual_fit(f, close_fit = 0.06)
  expect_s3_class(x, c("dmar_sem_fit", "dmar_tbl", "data.frame"))
  expect_identical(x$term, residual_terms)
  expect_identical(y$term, residual_terms)
  expect_true(all(is.na(sem_value(x, c("close_fit_bound", "p_close_fit")))))
  expect_equal(sem_value(y, "close_fit_bound"), 0.06)
  expect_true(is.na(sem_value(x, "srmr_mean")))
  expect_equal(sem_value(x, "N"), 301)
  # The default level is .90 (since DMAR 1.1.0, as every RMSEA interval's
  # and lavaan's), which pairs with a test of close fit at .05, and the
  # table states it; the exact fit p-value is at most .5.
  expect_equal(eval(formals(sem_residual_fit)$conf_level), 0.90)
  expect_equal(x$value, sem_residual_fit(f, conf_level = 0.90)$value)
  expect_equal(attr(x, "conf_level"), 0.90)
  expect_true("Confidence level: 90%" %in% utils::capture.output(print(x)))
  w <- sem_residual_fit(f, conf_level = 0.95)
  expect_lt(sem_value(w, "unbiased_srmr_ci_lower"),
            sem_value(x, "unbiased_srmr_ci_lower"))
  expect_true("Confidence level: 95%" %in% utils::capture.output(print(w)))
  expect_lte(sem_value(x, "p_exact_fit"), 0.5)
  r <- attr(x, "standardized_residuals")
  expect_equal(sem_value(x, "largest_absolute_residual"), max(abs(r[lower.tri(r)])))
  expect_equal(sem_value(x, "srmr"), sqrt(mean(r[lower.tri(r, diag = TRUE)]^2)))
  big <- attr(x, "largest_residual")
  expect_equal(abs(big$residual), sem_value(x, "largest_absolute_residual"))
  expect_equal(r[big$item_2, big$item_1], big$residual)
  out <- utils::capture.output(print(x))
  expect_true(any(grepl("Residual covariance: normal theory", out)))
})

test_that("a mean structure leaves the covariance SRMR alone and is reported apart", {
  skip_if_not_installed("lavaan")
  d <- lavaan::HolzingerSwineford1939
  f <- lavaan::cfa(hs_model, data = d, std.lv = TRUE)
  fm <- lavaan::cfa(hs_model, data = d, std.lv = TRUE, meanstructure = TRUE)
  x <- sem_residual_fit(f)
  y <- sem_residual_fit(fm)
  keep <- setdiff(residual_terms, "srmr_mean")
  expect_equal(sem_value(y, keep), sem_value(x, keep), tolerance = 1e-8)
  expect_equal(sem_value(y, "srmr_mean"), 0)
  restricted <- lavaan::cfa(paste(hs_model, "\n x4 ~ a*1\n x5 ~ a*1"), data = d,
                            std.lv = TRUE, meanstructure = TRUE)
  expect_warning(z <- sem_residual_fit(restricted), "restricts the means")
  s <- lavaan::lavResiduals(restricted, type = "cor.bentler", summary = TRUE)$summary
  expect_equal(sem_value(z, "srmr_mean"), as.numeric(s["srmr", "mean"]))
})

test_that("a covariance matrix gives the values the raw data give", {
  skip_if_not_installed("lavaan")
  d <- lavaan::HolzingerSwineford1939
  f <- lavaan::cfa(hs_model, data = d, std.lv = TRUE)
  g <- lavaan::cfa(hs_model, sample.cov = stats::cov(d[, paste0("x", 1:9)]),
                   sample.nobs = nrow(d), std.lv = TRUE)
  expect_equal(sem_residual_fit(g)$value, sem_residual_fit(f)$value,
               tolerance = 1e-6)
})

test_that("fits outside the computations are refused with a message that says why", {
  skip_if_not_installed("lavaan")
  d <- lavaan::HolzingerSwineford1939
  expect_error(sem_residual_fit(lavaan::cfa(hs_model, data = d, group = "school")),
               "single group")
  expect_error(sem_residual_fit(lavaan::cfa(hs_model, data = d, estimator = "GLS")),
               "maximum likelihood")
  expect_error(sem_residual_fit(lavaan::cfa(hs_model, data = d, estimator = "ULS")),
               "maximum likelihood")
  expect_error(sem_residual_fit(lavaan::cfa(hs_model, data = d,
                                            likelihood = "wishart")),
               "wishart")
  fiml_data <- d
  fiml_data$x1[1:10] <- NA
  expect_error(sem_residual_fit(lavaan::cfa(hs_model, data = fiml_data,
                                            missing = "ml")),
               "missing = \"ml\"")
  expect_error(sem_residual_fit(lavaan::sem("x3 ~ x1 + x2", data = d,
                                            conditional.x = TRUE)),
               "conditional.x")
  w <- d
  w$wt <- rep(c(0.5, 1.5), length.out = nrow(w))
  expect_error(sem_residual_fit(lavaan::cfa(hs_model, data = w,
                                            sampling.weights = "wt")),
               "unweighted")
  expect_error(sem_residual_fit(stats::lm(x1 ~ x2, data = d)), "lavaan model")
  f <- lavaan::cfa(hs_model, data = d)
  expect_error(sem_residual_fit(f, conf_level = 1), "below 1")
  expect_error(sem_residual_fit(f, close_fit = 0), "greater")
  skip_on_cran()
  ordered_data <- d
  for (v in paste0("x", 1:9)) {
    ordered_data[[v]] <- cut(d[[v]], stats::quantile(d[[v]], c(0, .25, .5, .75, 1)),
                             include.lowest = TRUE, labels = FALSE)
  }
  expect_error(sem_residual_fit(lavaan::cfa(hs_model, data = ordered_data,
                                            ordered = TRUE)),
               "continuous")
})

test_that("listwise deletion is named in a warning, and a saturated model explains its NA", {
  skip_if_not_installed("lavaan")
  d <- lavaan::HolzingerSwineford1939
  d$x1[1:20] <- NA
  expect_warning(x <- sem_residual_fit(lavaan::cfa(hs_model, data = d)),
                 "removed 20 of 301 cases")
  expect_equal(sem_value(x, "N"), 281)
  expect_equal(attr(x, "sem_fit")$N_original, 301)
  s <- lavaan::cfa("f =~ x1 + x2 + x3", data = lavaan::HolzingerSwineford1939)
  expect_warning(y <- sem_residual_fit(s), "saturated")
  expect_lt(sem_value(y, "srmr"), 1e-6)
  expect_true(all(is.na(sem_value(y, c("srmr_se_exact_fit", "p_exact_fit",
                                       "unbiased_srmr_se", "unbiased_srmr_ci_lower",
                                       "unbiased_srmr_ci_upper")))))
})

test_that("an older lavaan is refused before lavResiduals() is called", {
  expect_error(.sem_require_lavaan_version("sem_residual_fit", lavaan_version = "0.6-19"),
               "0.7-2 or newer")
  expect_true(.sem_require_lavaan_version("sem_residual_fit", lavaan_version = "0.7-2"))
  expect_true(.sem_require_lavaan_version("sem_residual_fit", lavaan_version = "0.10-1"))
})

test_that("a FIML or pairwise fit of complete data is accepted with the listwise values", {
  skip_if_not_installed("lavaan")
  d <- lavaan::HolzingerSwineford1939
  m <- "visual =~ x1 + x2 + x3
        textual =~ x4 + x5 + x6"
  # FIML on complete data is the listwise fit with a mean structure and
  # observed information, FIML's defaults; it was refused as if the data
  # were incomplete.
  ff <- lavaan::cfa(m, data = d, std.lv = TRUE, missing = "ml")
  fl <- lavaan::cfa(m, data = d, std.lv = TRUE, meanstructure = TRUE,
                    information = "observed")
  x <- sem_residual_fit(ff)
  expect_equal(x$value, sem_residual_fit(fl)$value, tolerance = 1e-10)
  expect_match(attr(x, "sem_fit")$footer[1], "observed information")
  fp <- lavaan::cfa(m, data = d, std.lv = TRUE, missing = "pairwise")
  expect_equal(sem_residual_fit(fp)$value,
               sem_residual_fit(lavaan::cfa(m, data = d, std.lv = TRUE))$value,
               tolerance = 1e-10)
  # Missing values, and two-stage estimation, are still refused.
  d$x1[1:10] <- NA
  expect_error(sem_residual_fit(lavaan::cfa(m, data = d, missing = "ml")),
               "with missing values in its analysis data")
  expect_error(sem_residual_fit(lavaan::cfa(m, data = lavaan::HolzingerSwineford1939,
                                            missing = "two.stage")),
               "missing = \"two.stage\"")
})

test_that("a fit with auxiliary variables is refused", {
  skip_if_not_installed("lavaan")
  fit <- suppressMessages(lavaan::cfa(hs_model, data = lavaan::HolzingerSwineford1939,
                                      missing = "fiml", aux = "ageyr"))
  expect_error(sem_residual_fit(fit), class = "dmar_refusal")
  expect_error(sem_residual_fit(fit), "auxiliary variables (ageyr)", fixed = TRUE)
})
