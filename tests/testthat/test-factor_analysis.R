# The bootstraps in this file are mechanics tests run with a small B; the
# rule on the number of replications is tested in test-dmar_bootstrap.R.
local_small_B()

fa_test_items <- function() c("t1_visual_perception", "t2_cubes", "t4_lozenges",
  "t6_paragraph_comprehension", "t7_sentence", "t9_word_meaning")

test_that("EFA reproduces independent ML covariance and rotated inference", {
  skip_if_not_installed("lavaan")
  d <- holzinger_swineford[fa_test_items()]
  x <- fa_exploratory(d, 2, seed = 113)
  fit <- attr(x, "fit")
  independent <- stats::factanal(d, factors = 2, rotation = "none")
  L <- unclass(independent$loadings)
  Rhat <- tcrossprod(L) + diag(independent$uniquenesses)
  expect_equal(stats::cov2cor(lavaan::fitted(fit)$cov), Rhat, tolerance = 3e-5)
  pe <- lavaan::parameterEstimates(fit, level = .95)
  raw <- x[x$metric == "original", ]
  expect_equal(raw$estimate, pe$est)
  expect_equal(raw$ci_lower, pe$ci.lower)
  expect_equal(raw$ci_upper, pe$ci.upper)
  ss <- lavaan::standardizedSolution(fit, level = .95)
  residual <- ss[ss$op == "~~" & ss$lhs == ss$rhs & ss$lhs %in% names(d), ]
  h <- x[x$op == "communality", ]
  expect_equal(h$estimate, 1 - residual$est.std)
  expect_equal(h$ci_lower, 1 - residual$ci.upper)
  expect_equal(h$ci_upper, 1 - residual$ci.lower)
  expect_identical(generics::glance(x)$estimator, "ML")
  expect_identical(generics::tidy(x)$term, x$term)
  expect_output(print(x), "condition on the selected factor count")
  expect_error(fa_exploratory(d, 4), "identified EFA")
  expect_error(fa_exploratory(d, 10), "identified EFA")
  expect_error(fa_exploratory(S = cov(d), N = nrow(d), estimator = "MLR"), "raw data")
})

test_that("EFA keeps raw FIML, robust SEs, ordered methods and caller RNG", {
  skip_if_not_installed("lavaan")
  d <- holzinger_swineford[fa_test_items()]
  d[seq(2, 300, 3), 2] <- NA
  set.seed(42); expected <- runif(3)
  set.seed(42)
  x <- fa_exploratory(d, 2, estimator = "MLR", seed = 113, conf_level = .99)
  expect_equal(runif(3), expected)
  expect_identical(attr(x, "method")$se, "robust.huber.white")
  expect_identical(attr(x, "method")$missing, "fiml")
  expect_equal(attr(x, "method")$N, nrow(d))
  fit <- attr(x, "fit")
  pe <- lavaan::parameterEstimates(fit, level = .99)
  expect_equal(x$se[x$metric == "original"], pe$se)
  expect_equal(x$ci_upper[x$metric == "original"], pe$ci.upper)
  for (j in seq_along(d)) d[[j]] <- as.integer(cut(d[[j]],
    unique(quantile(d[[j]], seq(0, 1, length.out = 5), na.rm = TRUE)), include.lowest = TRUE))
  expect_error(fa_exploratory(d, 2, ordered = TRUE), "Ordered EFA requires")
  o <- fa_exploratory(d, 2, estimator = "WLSMV", ordered = TRUE, missing = "pairwise", seed = 113)
  expect_identical(attr(o, "method")$estimator, "WLSMV")
  expect_identical(attr(o, "method")$missing, "pairwise")
  expect_true(all(is.finite(o$se[o$op == "=~"])))
})

test_that("parallel analysis uses the same one-factor reduction for each matrix", {
  d <- holzinger_swineford[fa_test_items()]
  set.seed(113)
  R <- cor(d)
  one <- stats::factanal(covmat = R, factors = 1, rotation = "none")
  diag(R) <- 1 - one$uniquenesses
  observed <- eigen(R, symmetric = TRUE, only.values = TRUE)$values
  reference <- replicate(12, {
    z <- matrix(rnorm(nrow(d) * ncol(d)), nrow(d), ncol(d))
    R <- cor(z)
    f <- stats::factanal(covmat = R, factors = 1, rotation = "none")
    diag(R) <- 1 - f$uniquenesses
    eigen(R, symmetric = TRUE, only.values = TRUE)$values
  })
  x <- parallel_analysis(d, B = 12, seed = 113)
  expect_equal(x$estimate, observed)
  expect_equal(unname(attr(x, "reference_draws")), t(reference))
  expect_equal(x$reference_quantile, apply(reference, 1L, quantile, .95), ignore_attr = TRUE)
  expect_equal(attr(x, "suggested_factors"), sum(cumprod(observed > x$reference_quantile)))
  # The Monte Carlo limits of the reference quantile are named for it: the
  # binomial order statistics of the 12 draws at the default 95 percent
  # (the upper one is past the last draw, so infinite).
  lo <- stats::qbinom(.025, 12, .95)
  hi <- stats::qbinom(.975, 12, .95) + 1L
  expect_equal(x$reference_quantile_ci_lower,
               apply(reference, 1L, function(v) sort(v)[lo]))
  expect_identical(hi > 12, TRUE)
  expect_identical(x$reference_quantile_ci_upper, rep(Inf, ncol(d)))
  expect_false(any(grepl("_mc_", names(x))))
  expect_equal(generics::glance(x)$attempted, 12)
  expect_output(print(x), "not sampling CIs")
  expect_false("conf_level" %in% names(x))
  set.seed(42); expected <- runif(2)
  set.seed(42); y <- parallel_analysis(d, B = 12, seed = 113)
  expect_equal(runif(2), expected)
  expect_equal(x, y)
})

test_that("parallel reference resampling and ordinal margins are explicit", {
  skip_if_not_installed("lavaan")
  d <- holzinger_swineford[fa_test_items()]
  for (j in seq_along(d)) d[[j]] <- as.integer(cut(d[[j]],
    unique(quantile(d[[j]], seq(0, 1, length.out = 4))), include.lowest = TRUE))
  x <- parallel_analysis(d, ordered = TRUE, missing = "pairwise",
                         reference = "resample", B = 8, seed = 113)
  expect_identical(attr(x, "method")$reference, "resample")
  expect_identical(attr(x, "method")$ordered, names(d))
  expect_equal(attr(x, "counts")$usable, 8)
  # Complete ordered items with 'missing' left out are analyzed listwise
  # (test-parallel_analysis.R); an explicit FIML request is refused.
  expect_error(parallel_analysis(d, ordered = TRUE, missing = "fiml", B = 8),
               "normal-theory FIML")
  expect_error(parallel_analysis(S = cov(d), N = nrow(d), reference = "resample"),
               "continuous normal reference only")
})

test_that("parallel analysis reduces covariance input to the named items", {
  nine <- grep("^t[1-9]_", names(holzinger_swineford), value = TRUE)
  S9 <- cov(holzinger_swineford[nine])
  N <- nrow(holzinger_swineford)
  items <- rev(fa_test_items())
  named <- parallel_analysis(S = S9, N = N, items = items, B = 12, seed = 113)
  reduced <- parallel_analysis(S = S9[items, items], N = N, B = 12, seed = 113)
  expect_identical(nrow(named), 6L)
  expect_equal(named, reduced)
  expect_identical(colnames(attr(named, "reference_draws")), items)
  expect_identical(nrow(parallel_analysis(S = S9, N = N, B = 12, seed = 113)), 9L)
  expect_error(parallel_analysis(S = S9, N = N, items = c(items[1:2], "t99_unknown"), B = 12),
               "at least three distinct items of 'S'. Not items of 'S': 't99_unknown'.",
               fixed = TRUE, class = "dmar_items_invalid")
  # Only the named items are checked: a missing entry and a variance of
  # zero outside them do not refuse the call.
  other <- setdiff(nine, items)
  holed <- S9
  holed[other[1], other[2]] <- holed[other[2], other[1]] <- NA
  holed[other[3], other[3]] <- 0
  expect_error(parallel_analysis(S = holed, N = N, B = 12),
               "'S' must be a finite symmetric covariance matrix", fixed = TRUE)
  expect_equal(parallel_analysis(S = holed, N = N, items = items, B = 12, seed = 113),
               reduced)
})

test_that("parallel analysis counts failed and defined boundary references separately", {
  d <- holzinger_swineford[fa_test_items()]
  state <- new.env(parent = emptyenv()); state$n <- 0L
  reduction <- DMAR:::.pa_reduced_eigen
  fake <- function(R, N) {
    state$n <- state$n + 1L
    if (state$n > 3L) stop("synthetic failed reference")
    reduction(R, N)
  }
  local_mocked_bindings(.pa_reduced_eigen = fake)
  expect_warning(x <- parallel_analysis(d, B = 6, seed = 113), "Thresholds and the factor count are withheld")
  expect_equal(attr(x, "counts")$usable, 2)
  expect_equal(attr(x, "counts")$failed, 4)
  expect_true(all(is.na(x$reference_quantile)))
  expect_true(is.na(attr(x, "suggested_factors")))
})

test_that("SAM inference matches an independent measurement and moment oracle", {
  skip_if_not_installed("lavaan")
  N <- 400L; p <- 8L
  loading <- c(.8, .7, .9, .75, .75, .85, .7, .8)
  L <- matrix(0, p, 2); L[1:4, 1] <- loading[1:4]; L[5:8, 2] <- loading[5:8]
  S <- L %*% matrix(c(1, .45, .45, 1), 2) %*% t(L) + diag(1 - loading^2)
  dimnames(S) <- list(paste0("y", 1:p), paste0("y", 1:p))
  positions <- which(lower.tri(S, diag = TRUE), arr.ind = TRUE)
  independent <- function(S) {
    L <- matrix(0, p, 2); theta <- numeric(p)
    for (b in 1:2) {
      ids <- (4 * (b - 1L) + 1L):(4 * b)
      fa <- stats::factanal(covmat = S[ids, ids], factors = 1, rotation = "none",
        control = list(lower = 1e-6, opt = list(factr = 10, pgtol = 1e-10)))
      original <- sqrt(diag(S)[ids]) * unclass(fa$loadings)[, 1L]
      L[ids, b] <- original / original[1L]
      theta[ids] <- diag(S)[ids] * fa$uniquenesses
    }
    M <- solve(crossprod(L, diag(1 / theta) %*% L), crossprod(L, diag(1 / theta)))
    V <- M %*% (S - diag(theta)) %*% t(M)
    unname(V[2, 1] / V[1, 1])
  }
  J <- vapply(seq_len(nrow(positions)), function(j) {
    i <- positions[j, 1]; k <- positions[j, 2]; h <- 5e-6
    plus <- minus <- S
    plus[i, k] <- plus[k, i] <- S[i, k] + h
    minus[i, k] <- minus[k, i] <- S[i, k] - h
    (independent(plus) - independent(minus)) / (2 * h)
  }, numeric(1))
  Gamma <- outer(seq_len(nrow(positions)), seq_len(nrow(positions)), Vectorize(function(a, b) {
    i <- positions[a, 1]; j <- positions[a, 2]; k <- positions[b, 1]; l <- positions[b, 2]
    S[i, k] * S[j, l] + S[i, l] * S[j, k]
  }))
  normal_var <- drop(J %*% Gamma %*% J) / N
  model <- "f1 =~ y1+y2+y3+y4\nf2 =~ y5+y6+y7+y8\nf2 ~ f1"
  set.seed(113)
  for (shape in c("normal", "skewed")) {
    raw <- matrix(rnorm(N * p), N, p)
    if (shape == "skewed") raw <- exp(raw * .7)
    raw <- scale(raw, center = TRUE, scale = FALSE)
    d <- raw %*% solve(chol(crossprod(raw) / N)) %*% chol(S)
    colnames(d) <- colnames(S)
    moments <- t(vapply(seq_len(N), function(i) (tcrossprod(d[i, ]) - S)[positions], numeric(nrow(positions))))
    hc0 <- sum(drop(moments %*% J)^2) / N^2
    c1 <- N * (N - 1) / ((N - 2) * (N - 3))
    c2 <- N / ((N - 2) * (N - 3))
    for (method in c("ML", "MLR")) for (unbiased in c(TRUE, FALSE)) {
      x <- sem_sam(model, as.data.frame(d), estimator = method, gamma_unbiased = unbiased)
      path <- x[x$lhs == "f2" & x$op == "~", ]
      expected_var <- if (method == "ML") normal_var else if (unbiased) c1 * hc0 - c2 * normal_var else hc0
      expect_equal(path$estimate, independent(S), tolerance = 1e-6)
      expect_equal(path$se / sqrt(expected_var), 1, tolerance = 1e-5)
      expect_equal(path$ci_lower, path$estimate - qnorm(.975) * path$se)
      expect_identical(generics::glance(x)$gamma_unbiased, unbiased)
      expect_true(any(x$stage == "measurement"))
    }
  }
})

test_that("SAM refuses unsupported analyses before changing their method or cases", {
  skip_if_not_installed("lavaan")
  d <- setNames(as.data.frame(matrix(seq_len(1200), 200, 6)), paste0("y", 1:6))
  model <- "f1 =~ y1+y2+y3\nf2 =~ y4+y5+y6\nf2 ~ f1"
  d[1, 1] <- NA
  expect_error(sem_sam(model, d), "No cases were deleted")
  expect_error(sem_sam(paste0(model, "\nf1 ~ f2"), d), "Feedback paths")
  expect_error(sem_sam(paste0(model, "\neffect := 2"), d), "derived parameters")
  expect_error(sem_sam(paste0(model, "\ny1 ~~ y2"), d), "residual covariances")
  expect_error(sem_sam(paste0(model, "\nf2 ~ y1"), d), "Observed covariates")
  # An indicator the data lack is the misnamed item of the factor analysis
  # functions, refused with their class, and named as an indicator: the
  # function has no 'items' argument for the message of .fa_data() to
  # speak of.
  e <- tryCatch(sem_sam(sub("y6", "y7", model, fixed = TRUE), d),
                dmar_items_invalid = function(e) e)
  expect_s3_class(e, c("dmar_items_invalid", "dmar_refusal", "error"))
  expect_identical(e$reason, "unknown")
  expect_identical(e$items, "y7")
  expect_identical(conditionMessage(e), "Indicators not found in 'data': y7.")
  # The same from a matrix, whose column names are read as .fa_data()
  # reads them.
  e <- tryCatch(sem_sam(sub("y6", "y7", model, fixed = TRUE), as.matrix(d)),
                dmar_items_invalid = function(e) e)
  expect_identical(e$items, "y7")
  # Columns of data that share a name are still the data's fault, found by
  # .fa_data() once every indicator is present.
  e <- tryCatch(sem_sam(model, cbind(d, d["y1"])),
                dmar_items_invalid = function(e) e)
  expect_identical(e$reason, "data_names")
  expect_identical(e$items, "y1")
})

test_that("EFA by ULS is the minimum residual solution, fitted to the correlation matrix", {
  skip_if_not_installed("lavaan")
  d <- holzinger_swineford[fa_test_items()]
  x <- fa_exploratory(d, 2, estimator = "ULS", seed = 113)
  expect_identical(attr(x, "method")$estimator, "ULS")
  expect_identical(attr(x, "method")$se, "robust.sem.nt")
  expect_identical(attr(x, "method")$structure, "correlation")
  expect_true(all(is.finite(x$se[x$op == "=~"])))
  # Independent control: the loadings that minimize the sum of squared
  # off-diagonal residuals of the correlation matrix (Harman & Jones, 1966),
  # found by direct minimization over an unrotated loading matrix from
  # several starts. Communalities and the minimized sum are rotation free.
  R <- stats::cor(d); p <- ncol(R); k <- 2L
  ssr <- function(l) {
    E <- R - tcrossprod(matrix(l, p, k)); sum(E[lower.tri(E)]^2)
  }
  best <- NULL
  for (s in 1:12) {
    set.seed(s)
    o <- stats::optim(stats::runif(p * k, -.5, .5), ssr, method = "BFGS",
                      control = list(maxit = 5000, reltol = 1e-15))
    if (is.null(best) || o$value < best$value) best <- o
  }
  L <- matrix(best$par, p, k)
  expect_equal(x$estimate[x$op == "communality"], rowSums(L^2), tolerance = 1e-6)
  Rhat <- stats::cov2cor(lavaan::fitted(attr(x, "fit"))$cov)
  expect_equal(sum(((R - Rhat)[lower.tri(R)])^2), best$value, tolerance = 1e-8)
  # The table carries no warning of its own, and the fit measures attribute
  # keeps only what the ULS statistic supports.
  expect_no_warning(fa_exploratory(d, 2, estimator = "ULS", seed = 113))
  fm <- attr(x, "fit_measures")
  expect_true(is.finite(fm$value[fm$term == "nfi"]))
  expect_true(is.finite(fm$value[fm$term == "srmr"]))
  expect_true(is.na(fm$value[fm$term == "rmsea"]))
  expect_true(is.na(fm$value[fm$term == "cfi"]))
  # The chi_square row is N - 1 times the ULS discrepancy of the correlations,
  # with no p-value: the number lavaan reports, not a chi square.
  expect_equal(fm$value[fm$term == "chi_square"],
               (nrow(d) - 1) * sum(((R - Rhat)[lower.tri(R)])^2), tolerance = 1e-6)
  expect_true(is.na(fm$value[fm$term == "p_chi_square"]))
  expect_output(print(x), "chi_square row is N - 1 times the ULS discrepancy")
  expect_identical(generics::glance(x)$structure, "correlation")
  # Other rotations keep the sandwich standard errors finite, rotation to a
  # partially specified target among them.
  keys <- list(f1 = fa_test_items()[1:3], f2 = fa_test_items()[4:6])
  for (rot in c("varimax", "geomin", "target")) {
    xr <- fa_exploratory(d, 2, estimator = "ULS", rotation = rot,
                         target = if (rot == "target") keys, seed = 113)
    expect_true(all(is.finite(xr$se[xr$op == "=~"])))
    expect_equal(xr$estimate[xr$op == "communality"], rowSums(L^2), tolerance = 1e-6)
  }
  # Ordered items are refused under ULS as under the other estimators.
  expect_error(fa_exploratory(d, 2, estimator = "ULS", ordered = fa_test_items()[1]), "WLSMV")
  expect_output(print(x), "minimum residual")
  # From a supplied covariance matrix the correlation matrix is fitted, the
  # communalities are the same, and the standard errors are NA.
  y <- fa_exploratory(S = stats::cov(d), N = nrow(d), nfactors = 2, estimator = "ULS", seed = 113)
  expect_equal(y$estimate[y$op == "communality"], x$estimate[x$op == "communality"],
               tolerance = 1e-8)
  expect_true(all(is.na(y$se)))
  expect_true(all(is.na(y$ci_lower)))
  expect_identical(generics::glance(y)$se_method, "none")
  expect_output(print(y), "supply raw data for standard errors")
  expect_false(grepl("Intervals condition", attr(y, "analysis_note"), fixed = TRUE))
  # A correlation matrix given as S is the same fit.
  yr <- fa_exploratory(S = stats::cor(d), N = nrow(d), nfactors = 2, estimator = "ULS", seed = 113)
  expect_equal(yr$estimate[yr$op == "communality"], y$estimate[y$op == "communality"],
               tolerance = 1e-10)
})

test_that("EFA by ULS refuses incomplete data rather than deleting or imputing silently", {
  skip_if_not_installed("lavaan")
  d <- holzinger_swineford[fa_test_items()]
  d[seq(2, 300, 3), 2] <- NA
  expect_error(fa_exploratory(d, 2, estimator = "ULS"), "has no likelihood")
  expect_error(fa_exploratory(d, 2, estimator = "ULS", missing = "ml"), "has no likelihood")
  expect_error(fa_exploratory(d, 2, estimator = "ULS", missing = "pairwise"), "pairwise deletion")
  z <- fa_exploratory(d, 2, estimator = "ULS", missing = "listwise", seed = 113)
  expect_equal(attr(z, "method")$N, sum(stats::complete.cases(d)))
  expect_identical(attr(z, "method")$missing, "listwise")
  # Complete data under the default missing setting are accepted as they are.
  expect_no_error(fa_exploratory(d[stats::complete.cases(d), ], 2, estimator = "ULS", seed = 113))
  expect_error(fa_exploratory(S = stats::cov(d, use = "complete.obs"), N = 200, nfactors = 2,
                   estimator = "MLR"), "raw data")
})

test_that("efa_k() is fa_exploratory() under a shorter name", {
  skip_if_not_installed("lavaan")
  d <- holzinger_swineford[fa_test_items()]
  a <- efa_k(d, 2, seed = 113)
  b <- fa_exploratory(d, 2, seed = 113)
  # The lavaan object carries wall-clock timings, so it is compared by its
  # estimates and the table by everything else.
  strip <- function(x) { attr(x, "fit") <- NULL; x }
  expect_identical(strip(a), strip(b))
  expect_identical(lavaan::coef(attr(a, "fit")), lavaan::coef(attr(b, "fit")))
  expect_identical(formals(efa_k), formals(fa_exploratory))
})

test_that("EFA display shows the solution, with complexity and variance accounted for", {
  skip_if_not_installed("lavaan")
  d <- holzinger_swineford[fa_test_items()]
  x <- fa_exploratory(d, 2, seed = 113)
  expect_s3_class(x, "dmar_fa_exploratory")
  expect_s3_class(x, "dmar_factor_analysis")
  # Complexity rows: Hofmann's index from the standardized pattern loadings.
  load <- x[x$op == "=~" & x$metric == "standardized", ]
  L <- matrix(load$estimate, ncol = 2, dimnames = list(unique(load$rhs), unique(load$lhs)))
  cx <- x[x$op == "complexity", ]
  expect_identical(cx$lhs, rownames(L))
  expect_equal(cx$estimate, unname(rowSums(L^2)^2 / rowSums(L^4)))
  expect_true(all(is.na(cx$se)) && all(is.na(cx$ci_lower)))
  expect_true(all(cx$estimate >= 1 & cx$estimate <= 2))
  expect_true("complexity" %in% generics::tidy(x)$op)
  # Variance accounted for: diag(Phi Lambda' Lambda), summing to the
  # communalities, with proportions of the items' and of the common variance.
  phi <- x$estimate[x$op == "~~" & x$metric == "standardized" & x$lhs == "f1" & x$rhs == "f2"]
  Phi <- matrix(c(1, phi, phi, 1), 2)
  va <- attr(x, "variance_accounted")
  expect_identical(va$factor, c("f1", "f2", "total"))
  expect_equal(va$sum_sq_loadings[1:2], unname(diag(Phi %*% crossprod(L))))
  expect_equal(va$sum_sq_loadings[3], sum(x$estimate[x$op == "communality"]))
  expect_equal(va$proportion_variance, va$sum_sq_loadings / 6)
  expect_equal(va$proportion_common, va$sum_sq_loadings / va$sum_sq_loadings[3])
  # The display: header, loadings matrix with intervals, communality,
  # complexity, the variance table, the factor correlations, and the fit.
  local_reproducible_output(width = 200)
  out <- capture.output(print(x))
  expect_match(out, "^Exploratory factor analysis: 2 factors, 6 items, N = 301$", all = FALSE)
  expect_match(out, "rotation: oblimin \\(oblique\\); SE: standard", all = FALSE)
  expect_match(out, "^ *item +f1 +f2 +communality +complexity *$", all = FALSE)
  expect_match(out, "t4_lozenges +0\\.766 \\[0\\.629, 0\\.903\\]", all = FALSE)
  expect_match(out, "Variance accounted for", all = FALSE)
  expect_match(out, "^ *total +3\\.379 +0\\.563 +1\\.000 *$", all = FALSE)
  expect_match(out, "Factor correlations with 95% intervals", all = FALSE)
  expect_match(out, "f1 +1 +0\\.302 \\[0\\.145, 0\\.459\\]", all = FALSE)
  # The fit measures are measured values, so an RMSEA of 0 and a CFI of 1
  # keep their decimals; the degrees of freedom are a count.
  # The RMSEA interval is at rmsea_conf_level, 90 percent by default, and
  # the line names it.
  expect_match(out, "^Fit: chi square 2\\.109, df 4, p = 0\\.7157; RMSEA with 90% interval 0\\.000 \\[0\\.000, 0\\.064\\]; CFI 1\\.000; TLI 1\\.011; SRMR 0\\.008$", all = FALSE)
  one <- capture.output(print(fa_exploratory(d, 1, seed = 113)))
  expect_match(one, "^Fit: chi square 103\\.230, df 9, p < 0\\.0001; RMSEA with 90% interval 0\\.187", all = FALSE)
  # With one factor every complexity is 1, some of them only up to rounding
  # error, and each reads 1.000.
  expect_match(one, "^ *t2_cubes +0\\.212 \\[0\\.097, 0\\.328\\] +0\\.045 \\[-0\\.004, 0\\.094\\] +1\\.000$", all = FALSE)
  expect_false(any(grepl("\\] +1$", one)))
  expect_match(out, "^Version: robust", all = FALSE)
  expect_match(out, "Confidence level: 95%", all = FALSE)
  expect_match(out, "condition on the selected factor count", all = FALSE)
  expect_false(any(grepl("~1", out, fixed = TRUE)))
  # A subset loses the attributes and prints as an ordinary table.
  sub <- capture.output(print(x[x$op == "=~", ]))
  expect_false(any(grepl("Variance accounted for", sub)))
  expect_match(sub, "f1 =~ t1_visual_perception", all = FALSE)
  # knit_print() and as_kable() show the same display.
  skip_if_not_installed("knitr")
  kb <- as_kable(x, format = "pipe")
  expect_match(as.character(kb), "communality", all = FALSE)
  kp <- knitr::knit_print(x)
  expect_match(as.character(kp), "Variance accounted for", all = FALSE)
  expect_match(as.character(kp), "Fit: chi square", all = FALSE)
})

test_that("EFA rotations run in the form named, orthogonal ones with zero factor correlations", {
  skip_if_not_installed("lavaan")
  d <- holzinger_swineford[fa_test_items()]
  rotations <- DMAR:::.efa_rotations()
  # The fit does not depend on the rotation: every criterion reports the
  # unrotated model's measures (lavaan's FIML-corrected indices of a
  # rotated object would differ by rotation).
  unrotated <- attr(fa_exploratory(d, 2, rotation = "none", seed = 113), "fit_measures")
  # The two target rotations (lavaan's "pst") take a target: each item at
  # zero on the factor it is not keyed to.
  keys <- list(f1 = fa_test_items()[1:3], f2 = fa_test_items()[4:6])
  for (r in rotations$name) {
    pst <- rotations$lavaan[rotations$name == r] == "pst"
    x <- fa_exploratory(d, 2, rotation = r, target = if (pst) keys, seed = 113)
    m <- attr(x, "method")
    expect_equal(attr(x, "fit_measures")$value, unrotated$value, info = r)
    expect_identical(m$rotation, r, info = r)
    expect_identical(m$orthogonal, rotations$orthogonal[rotations$name == r], info = r)
    expect_true(all(is.finite(x$se[x$op == "=~"])), info = r)
    phi <- x[x$op == "~~" & x$metric == "standardized" & x$lhs == "f1" & x$rhs == "f2", ]
    if (m$orthogonal) {
      expect_true(nrow(phi) == 0L || isTRUE(all.equal(phi$estimate, 0)), info = r)
    } else {
      expect_true(is.finite(phi$se), info = r)
    }
    expect_identical(lavaan::lavInspect(attr(x, "fit"), "options")$rotation,
                     rotations$lavaan[rotations$name == r], info = r)
  }
  v <- fa_exploratory(d, 2, rotation = "varimax", seed = 113)
  out <- capture.output(print(v))
  expect_match(out, "rotation: varimax \\(orthogonal\\), Kaiser row weights", all = FALSE)
  expect_match(out, "Factor correlations: fixed at zero \\(orthogonal rotation\\)", all = FALSE)
  u <- fa_exploratory(d, 2, rotation = "none", seed = 113)
  expect_match(capture.output(print(u)), "rotation: none \\(unrotated\\)", all = FALSE)
  expect_error(fa_exploratory(d, 2, rotation = "promax"), "arg")
  one <- fa_exploratory(d, 1, seed = 113)
  expect_identical(attr(one, "variance_accounted")$factor, "f1")
  expect_equal(one$estimate[one$op == "complexity"], rep(1, 6))
  expect_false(any(grepl("Factor correlations", capture.output(print(one)))))
})

test_that("EFA reports each factor's determinacy and H with delta method intervals", {
  skip_if_not_installed("lavaan")
  d <- holzinger_swineford[fa_test_items()]
  x <- fa_exploratory(d, 2, seed = 113)
  fit <- attr(x, "fit")
  det <- x[x$op == "determinacy", ]
  H <- x[x$op == "H", ]
  expect_identical(det$lhs, c("f1", "f2"))
  expect_identical(det$rhs, det$lhs)
  expect_identical(H$lhs, c("f1", "f2"))
  expect_true(all(c(det$metric, H$metric) == "standardized"))
  expect_true(all(is.na(c(det$statistic, det$p_value, H$statistic, H$p_value))))
  # Determinacy: the regression score determinacy lavaan and factor_scores()
  # report, the square root of diag(Phi Lambda' R^-1 Lambda Phi) with R the
  # correlations the model implies.
  expect_equal(det$estimate, unname(lavaan::lavInspect(fit, "fs.determinacy")), tolerance = 1e-8)
  fs <- attr(factor_scores(fit), "determinacy")
  expect_equal(det$estimate, fs$determinacy_regression, tolerance = 1e-8)
  load <- x[x$op == "=~" & x$metric == "standardized", ]
  L <- matrix(load$estimate, ncol = 2, dimnames = list(unique(load$rhs), unique(load$lhs)))
  R <- stats::cov2cor(lavaan::lavInspect(fit, "implied")$cov)[rownames(L), rownames(L)]
  phi <- x$estimate[x$op == "~~" & x$metric == "standardized" & x$lhs == "f1" & x$rhs == "f2"]
  Phi <- matrix(c(1, phi, phi, 1), 2)
  expect_equal(det$estimate^2, unname(diag(Phi %*% t(L) %*% solve(R, L) %*% Phi)), tolerance = 1e-8)
  # H: Hancock and Mueller's closed form on each factor's standardized loadings.
  expect_equal(H$estimate, unname(apply(L, 2, function(l) {
    t <- sum(l^2 / (1 - l^2)); t / (1 + t)
  })))
  # Intervals: finite standard errors on the index scale, and limits formed
  # on the logit scale, inside (0, 1) and around the estimate.
  expect_true(all(is.finite(c(det$se, H$se))))
  expect_true(all(c(det$ci_lower, H$ci_lower) > 0 & c(det$ci_upper, H$ci_upper) < 1))
  expect_true(all(det$ci_lower < det$estimate & det$estimate < det$ci_upper))
  expect_true(all(H$ci_lower < H$estimate & H$estimate < H$ci_upper))
  logit <- function(p) log(p / (1 - p))
  expect_equal(logit(H$ci_upper) - logit(H$estimate),
               qnorm(0.975) * H$se / (H$estimate * (1 - H$estimate)))
  expect_equal(logit(det$estimate) - logit(det$ci_lower),
               qnorm(0.975) * det$se / (det$estimate * (1 - det$estimate)))
  # The attribute holds the same values, one row per factor, and tidy() the rows.
  fq <- attr(x, "factor_quality")
  expect_identical(fq$factor, c("f1", "f2"))
  expect_identical(names(fq), c("factor", "determinacy", "determinacy_se",
                                "determinacy_ci_lower", "determinacy_ci_upper",
                                "min_score_correlation", "min_score_correlation_se",
                                "min_score_correlation_ci_lower",
                                "min_score_correlation_ci_upper",
                                "H", "H_se", "H_ci_lower", "H_ci_upper"))
  # 2 rho^2 - 1, the lowest correlation of two equally consistent score
  # sets (Guttman, 1955): its delta method SE, and the determinacy's
  # limits mapped through the increasing map.
  expect_equal(fq$min_score_correlation, 2 * fq$determinacy^2 - 1)
  expect_equal(fq$min_score_correlation_se, 4 * fq$determinacy * fq$determinacy_se)
  expect_equal(fq$min_score_correlation_ci_lower, 2 * fq$determinacy_ci_lower^2 - 1)
  expect_equal(fq$min_score_correlation_ci_upper, 2 * fq$determinacy_ci_upper^2 - 1)
  expect_equal(fq$determinacy, det$estimate)
  expect_equal(fq$determinacy_se, det$se)
  expect_equal(fq$H_ci_lower, H$ci_lower)
  expect_equal(fq$H_ci_upper, H$ci_upper)
  expect_true(all(c("determinacy", "H") %in% generics::tidy(x)$op))
  # With one factor H is the squared determinacy, and its standard error is
  # the one fa_confirmatory() (lavaan's defined parameter, under the same
  # FIML fit with observed information) and reliability_H() (central
  # differences over the free parameters of the same fit) give; the logit
  # limits are reliability_H()'s default limits.
  one <- fa_exploratory(d, 1, seed = 113)
  d1 <- one[one$op == "determinacy", ]
  h1 <- one[one$op == "H", ]
  expect_equal(h1$estimate, d1$estimate^2)
  cfa <- fa_confirmatory(d, factors = list(f1 = names(d)), missing = "fiml")
  expect_equal(h1$estimate, cfa$estimate[cfa$term == "H_f1"], tolerance = 1e-6)
  expect_equal(h1$se, cfa$se[cfa$term == "H_f1"], tolerance = 1e-5)
  expect_equal(h1$ci_lower, cfa$ci_lower[cfa$term == "H_f1"], tolerance = 1e-5)
  expect_equal(h1$ci_upper, cfa$ci_upper[cfa$term == "H_f1"], tolerance = 1e-5)
  rh <- reliability_H(fit = attr(one, "fit"))
  expect_equal(h1$se, rh$value[rh$term == "se"], tolerance = 1e-5)
  expect_equal(h1$ci_lower, rh$value[rh$term == "ci_lower"], tolerance = 1e-5)
  expect_equal(h1$ci_upper, rh$value[rh$term == "ci_upper"], tolerance = 1e-5)
  # An orthogonal rotation: Phi is the identity.
  v <- fa_exploratory(d, 2, rotation = "varimax", seed = 113)
  lv <- v[v$op == "=~" & v$metric == "standardized", ]
  Lv <- matrix(lv$estimate, ncol = 2, dimnames = list(unique(lv$rhs), unique(lv$lhs)))
  Rv <- stats::cov2cor(lavaan::lavInspect(attr(v, "fit"), "implied")$cov)[rownames(Lv), rownames(Lv)]
  expect_equal(v$estimate[v$op == "determinacy"]^2, unname(diag(t(Lv) %*% solve(Rv, Lv))), tolerance = 1e-8)
  expect_true(all(is.finite(v$se[v$op %in% c("determinacy", "H")])))
  # Without standard errors (ULS from a supplied matrix) the estimates stand
  # alone; ULS from raw data carries the sandwich.
  us <- fa_exploratory(S = stats::cov(d), N = nrow(d), nfactors = 2, estimator = "ULS", seed = 113)
  expect_true(all(is.na(us$se[us$op %in% c("determinacy", "H")])))
  expect_true(all(is.na(us$ci_lower[us$op %in% c("determinacy", "H")])))
  expect_true(all(is.finite(us$estimate[us$op %in% c("determinacy", "H")])))
  ur <- fa_exploratory(d, 2, estimator = "ULS", seed = 113)
  expect_true(all(is.finite(ur$se[ur$op %in% c("determinacy", "H")])))
  # The display: a block under the factor correlations, the estimates alone
  # when the fit has no standard errors.
  local_reproducible_output(width = 200)
  out <- capture.output(print(x))
  expect_match(out, "^Factor determinacy and coefficient H with 95% intervals:$", all = FALSE)
  expect_match(out, "^ *factor +determinacy +min_r +H *$", all = FALSE)
  expect_match(out, "^ *f1 +0\\.833 \\[0\\.7[0-9]+, 0\\.8[0-9]+\\] +0\\.38[0-9] \\[-?[0-9.]+, -?[0-9.]+\\] +0\\.678 \\[0\\.5[0-9]+, 0\\.7[0-9]+\\] *$", all = FALSE)
  expect_match(out, "min_r: the lowest correlation", all = FALSE)
  expect_match(out, "limits from the logit scale", all = FALSE)
  expect_true(which(grepl("^Factor determinacy", out)) > which(grepl("^Factor correlations", out)))
  expect_true(which(grepl("^Factor determinacy", out)) < which(grepl("^Fit: ", out)))
  outs <- capture.output(print(us))
  expect_match(outs, "^Factor determinacy and coefficient H:$", all = FALSE)
  expect_false(any(grepl("logit scale", outs)))
  skip_if_not_installed("knitr")
  expect_match(as.character(knitr::knit_print(x)), "Factor determinacy", all = FALSE)
})


test_that("fa_exploratory() takes listwise deletion for ML and its default for complete ordered items", {
  skip_if_not_installed("lavaan")
  d <- holzinger_swineford[c("t1_visual_perception", "t2_cubes", "t4_lozenges",
                             "t6_paragraph_comprehension", "t7_sentence",
                             "t9_word_meaning")]
  d$t2_cubes[1:30] <- NA
  lw <- fa_exploratory(d, 2, missing = "listwise", seed = 113)
  expect_identical(attr(lw, "method")$missing, "listwise")
  expect_identical(attr(lw, "method")$N, sum(stats::complete.cases(d)))
  fi <- fa_exploratory(d, 2, seed = 113)
  expect_identical(attr(fi, "method")$missing, "fiml")
  expect_match(paste(capture.output(print(fi)), collapse = "\n"), "FIML")
  o <- holzinger_swineford[names(d)]
  o[] <- lapply(o, function(v) cut(v, 4, labels = FALSE))
  x <- fa_exploratory(o, 2, estimator = "WLSMV", ordered = TRUE, seed = 113)
  expect_identical(attr(x, "method")$missing, "listwise")
  o$t2_cubes[1:30] <- NA
  expect_error(fa_exploratory(o, 2, estimator = "WLSMV", ordered = TRUE),
               "when item values are missing")
})
