.cf_spatial <- c("t1_visual_perception", "t2_cubes", "t3_paper_form_board",
                 "t4_lozenges", "t25_paper_form_board_r", "t26_flags")

.cf_shape <- c("S", "mu", "N", "N_complete", "patterns", "converged")

# An EM algorithm for the saturated multivariate normal model written out
# in base R, one pattern at a time, run to a far tighter criterion than
# lavaan's default: the independent oracle for cov_fiml().
.cf_em_oracle <- function(x, tol = 1e-12, max_iter = 20000L) {
  x <- as.matrix(x)
  n <- nrow(x)
  p <- ncol(x)
  obs <- !is.na(x)
  key <- apply(obs, 1L, paste, collapse = "")
  groups <- split(seq_len(n), key)
  mu <- colMeans(x, na.rm = TRUE)
  S <- diag(apply(x, 2L, stats::var, na.rm = TRUE), p)
  for (it in seq_len(max_iter)) {
    T1 <- numeric(p)
    T2 <- matrix(0, p, p)
    for (g in groups) {
      o <- obs[g[1L], ]
      m <- !o
      xg <- x[g, , drop = FALSE]
      C <- matrix(0, p, p)
      if (any(m)) {
        B <- S[m, o, drop = FALSE] %*% solve(S[o, o, drop = FALSE])
        centered <- sweep(xg[, o, drop = FALSE], 2L, mu[o])
        xg[, m] <- sweep(centered %*% t(B), 2L, mu[m], `+`)
        C[m, m] <- S[m, m, drop = FALSE] - B %*% S[o, m, drop = FALSE]
      }
      T1 <- T1 + colSums(xg)
      T2 <- T2 + crossprod(xg) + length(g) * C
    }
    mu_new <- T1 / n
    S_new <- T2 / n - tcrossprod(mu_new)
    change <- max(abs(c(mu_new - mu, S_new - S)))
    mu <- mu_new
    S <- S_new
    if (change < tol) break
  }
  dimnames(S) <- list(colnames(x), colnames(x))
  names(mu) <- colnames(x)
  list(S = S, mu = mu, iterations = it)
}

# The six spatial tests with values on the cubes test (t2) also missing at
# random given visual perception (t1), so there are four patterns.
.cf_spatial_mar <- function(seed = 113) {
  d <- holzinger_swineford[.cf_spatial]
  set.seed(seed)
  z <- as.numeric(scale(d$t1_visual_perception))
  d$t2_cubes[stats::runif(nrow(d)) < stats::plogis(-1 + z)] <- NA
  d
}

test_that("complete data give the sample moments with the maximum likelihood divisor", {
  skip_on_cran()
  skip_if_not_installed("lavaan")
  d <- holzinger_swineford[.cf_spatial[1:4]]
  res <- cov_fiml(d)
  expect_type(res, "list")
  expect_identical(names(res), .cf_shape)
  n <- nrow(d)
  expect_equal(res$S, stats::cov(d) * (n - 1) / n, tolerance = 1e-10)
  expect_equal(res$mu, colMeans(d), tolerance = 1e-12)
  expect_identical(res$N, 301L)
  expect_identical(res$N_complete, 301L)
  expect_identical(res$patterns, 1L)
  expect_true(res$converged)
})

test_that("the spatial tests use every pupil, and the moments are lavaan's h1 moments", {
  skip_if_not_installed("lavaan")
  d <- holzinger_swineford[.cf_spatial]
  res <- cov_fiml(d)
  expect_identical(names(res), .cf_shape)
  expect_identical(dimnames(res$S), list(.cf_spatial, .cf_spatial))
  expect_identical(names(res$mu), .cf_spatial)
  expect_true(isSymmetric(res$S))
  expect_true(all(eigen(res$S, only.values = TRUE)$values > 0))
  expect_identical(res$N, 301L)
  expect_identical(res$N_complete, 145L)
  expect_identical(res$patterns, 2L)
  expect_true(res$converged)
  # A test the Pasteur pupils took keeps its full data mean; the tests
  # they did not take borrow strength from the ones they did.
  expect_equal(res$mu[["t1_visual_perception"]],
               mean(d$t1_visual_perception), tolerance = 1e-10)

  fit <- lavaan::cfa("f =~ t1_visual_perception + t2_cubes +
                      t3_paper_form_board + t4_lozenges +
                      t25_paper_form_board_r + t26_flags",
                     data = d, missing = "ml", std.lv = TRUE)
  h1 <- lavaan::lavInspect(fit, "sampstat.h1")
  # The EM runs on standardized variables (so that its criterion does not
  # depend on the units), lavaan's in the raw metric; both stop at
  # lavaan's criterion, so they agree to about seven significant digits
  # rather than to the last one.
  expect_equal(res$S, unclass(h1$cov)[.cf_spatial, .cf_spatial],
               tolerance = 1e-6, ignore_attr = TRUE)
  expect_equal(unname(res$mu), as.numeric(h1$mean[.cf_spatial]),
               tolerance = 1e-6)
})

test_that("the estimates match an EM algorithm written out in base R", {
  skip_on_cran()
  skip_if_not_installed("lavaan")
  d <- .cf_spatial_mar()
  res <- cov_fiml(d)
  oracle <- .cf_em_oracle(d)
  expect_lt(oracle$iterations, 20000L)
  expect_equal(res$S, oracle$S, tolerance = 1e-6)
  expect_equal(res$mu, oracle$mu, tolerance = 1e-6)
  expect_identical(res$patterns, 4L)
  expect_identical(res$N_complete,
                   sum(stats::complete.cases(d)))
})

test_that("the matrix is the reliability family's FIML basis", {
  skip_on_cran()  # Additional model fits run in the full local suite.
  skip_if_not_installed("lavaan")
  d <- .cf_spatial_mar()[.cf_spatial[1:4]]
  res <- cov_fiml(d)
  # .fiml_cov() reads the same EM estimates (since 1.1.0; through 1.0.1 it
  # maximized the likelihood directly and agreed to about 1e-5).
  basis <- .fiml_cov(d, items = names(d))
  expect_true(basis$converged)
  expect_identical(basis$S, res$S)
  J <- ncol(res$S)
  alpha_em <- (J / (J - 1)) * (1 - sum(diag(res$S)) / sum(res$S))
  alpha_fiml <- reliability_alpha(d, missing = "fiml", ci_method = "none")
  expect_equal(alpha_em, alpha_fiml$value[alpha_fiml$term == "estimate"],
               tolerance = 1e-12)
})

test_that("an auxiliary informs the estimates without entering them", {
  skip_on_cran()  # Additional model fits run in the full local suite.
  skip_if_not_installed("lavaan")
  d <- .cf_spatial_mar()[c("t1_visual_perception", "t2_cubes",
                           "t4_lozenges")]
  vars <- c("t2_cubes", "t4_lozenges")
  with_aux <- cov_fiml(d, aux = "t1_visual_perception")
  everything <- cov_fiml(d)
  expect_identical(dimnames(with_aux$S), list(vars, vars))
  expect_identical(names(with_aux$mu), vars)
  # The analysis block of the joint EM estimates is what is returned.
  expect_equal(with_aux$S, everything$S[vars, vars], tolerance = 1e-12)
  expect_equal(with_aux$mu, everything$mu[vars], tolerance = 1e-12)
  # Missingness on t2 depends on t1, so leaving t1 out moves the mean of
  # t2 back toward that of the pupils who took it.
  without_aux <- cov_fiml(d[vars])
  full <- mean(holzinger_swineford$t2_cubes)
  off <- c(with_aux = abs(with_aux$mu[["t2_cubes"]] - full),
           without_aux = abs(without_aux$mu[["t2_cubes"]] - full),
           respondents = abs(mean(d$t2_cubes, na.rm = TRUE) - full))
  expect_identical(order(off), 1:3)
  expect_identical(with_aux$patterns, 2L)
  expect_identical(without_aux$patterns, 2L)
})

test_that("a case observed only on the auxiliaries is used", {
  skip_on_cran()  # Additional model fits run in the full local suite.
  skip_if_not_installed("lavaan")
  set.seed(113)
  n <- 2000
  z <- stats::rnorm(n)
  d <- data.frame(y1 = 0.6 * z + stats::rnorm(n, sd = 0.8),
                  y2 = 0.6 * z + stats::rnorm(n, sd = 0.8),
                  z = z)
  # The whole analysis block goes missing, more often for high z: MAR
  # given z, and only the cases that answered nothing else carry it.
  gone <- stats::runif(n) < stats::plogis(-0.5 + 1.5 * z)
  d[gone, c("y1", "y2")] <- NA
  res <- cov_fiml(d, aux = "z")
  expect_identical(res$N, as.integer(n))
  expect_identical(res$N_complete, as.integer(sum(!gone)))
  fit <- lavaan::lavCor(d, missing = "ml", output = "fit")
  h1 <- lavaan::lavInspect(fit, "sampstat.h1")
  expect_equal(unname(res$mu), as.numeric(h1$mean[c("y1", "y2")]),
               tolerance = 1e-6)
  # The population means are zero; the respondents alone are well below.
  respondents <- colMeans(d[!gone, c("y1", "y2")])
  expect_lt(max(abs(res$mu)), 0.1)
  expect_true(all(respondents < -0.15))
})

test_that("a case with every variable missing is set aside without a warning", {
  skip_on_cran()
  skip_if_not_installed("lavaan")
  d <- holzinger_swineford[.cf_spatial[1:3]]
  d[1:3, ] <- NA
  d[4:6, 1:2] <- NA
  expect_no_warning(res <- cov_fiml(d, aux = "t3_paper_form_board"))
  expect_identical(res$N, 298L)
  expect_identical(res$N_complete, 295L)
  expect_identical(res$patterns, 2L)
  # One analysis variable is allowed, and its unobserved cases drop out.
  one <- cov_fiml(holzinger_swineford["t25_paper_form_board_r"])
  expect_identical(one$N, 145L)
  expect_equal(one$S[1, 1],
               stats::var(stats::na.omit(
                 holzinger_swineford$t25_paper_form_board_r)) * 144 / 145,
               tolerance = 1e-10)
})

test_that("a matrix and awkward column names are accepted", {
  skip_on_cran()  # Additional model fits run in the full local suite.
  skip_if_not_installed("lavaan")
  m <- as.matrix(holzinger_swineford[c("t1_visual_perception", "t2_cubes")])
  m[1:5, 1] <- NA
  named <- cov_fiml(m)
  colnames(m) <- NULL
  unnamed <- cov_fiml(m)
  expect_identical(dimnames(unnamed$S), list(c("V1", "V2"), c("V1", "V2")))
  expect_equal(unnamed$S, named$S, ignore_attr = TRUE, tolerance = 1e-12)
  odd <- data.frame(`item one` = m[, 1], `item+2` = m[, 2],
                    check.names = FALSE)
  res <- cov_fiml(odd)
  expect_identical(names(res$mu), c("item one", "item+2"))
  expect_equal(res$S, named$S, ignore_attr = TRUE, tolerance = 1e-12)
})

test_that("a stop at the iteration limit is reported, not hidden", {
  skip_on_cran()  # Additional model fits run in the full local suite.
  skip_if_not_installed("lavaan")
  x <- as.matrix(holzinger_swineford[grep("^t[0-9]",
                                          names(holzinger_swineford))])
  expect_no_warning(capped <- .cov_fiml_em(x, em_args = list(max_iter = 10L)))
  expect_false(capped$converged)
  expect_true(.cov_fiml_em(x)$converged)

  real_em <- .cov_fiml_em
  local_mocked_bindings(.cov_fiml_em = function(x, em_args = NULL) {
    real_em(x, em_args = list(max_iter = 10L))
  })
  expect_warning(res <- cov_fiml(holzinger_swineford[.cf_spatial]),
                 "iteration limit")
  expect_identical(names(res), .cf_shape)
  expect_false(res$converged)
})

test_that("the estimates do not depend on the units of the variables", {
  skip_on_cran()  # Additional model fits run in the full local suite.
  skip_if_not_installed("lavaan")
  # The reviewer's case: x3 missing at random given x1 for about 65 percent
  # of the cases. lavaan's EM stops when no mean or covariance changes by
  # more than 1e-5, an absolute criterion; run in the raw metric with every
  # variable multiplied by .001 (variances near 1e-6) it stopped after one
  # accelerated step, with converged = TRUE and the correlations of x3
  # at .115, .054, and .170 in place of .265, .185, and .334.
  set.seed(3)
  n <- 200
  R <- matrix(0.5, 4, 4)
  diag(R) <- 1
  x <- as.data.frame(matrix(stats::rnorm(n * 4), n) %*% chol(R))
  names(x) <- c("y", "x1", "x2", "x3")
  x$x3[stats::runif(n) < stats::plogis(1 + 2 * x$x1)] <- NA
  x$x2[sample(n, 3)] <- NA
  base <- cov_fiml(x)
  expect_equal(unname(stats::cov2cor(base$S)[1:3, "x3"]),
               c(0.265, 0.185, 0.334), tolerance = 0.002)
  for (k in c(1e-3, 1e3)) {
    expect_no_warning(scaled <- cov_fiml(x * k))
    expect_true(scaled$converged)
    expect_equal(scaled$S / k^2, base$S, tolerance = 1e-10, info = k)
    expect_equal(scaled$mu / k, base$mu, tolerance = 1e-10, info = k)
  }
  # One variable at a time, with a shift: S becomes D S D and mu becomes
  # a + D mu.
  d <- c(1, 1e-4, 1, 1e4)
  shifted <- cov_fiml(sweep(sweep(x, 2L, d, `*`), 2L, c(0, 5, 0, -3), `+`))
  expect_equal(shifted$S, base$S * outer(d, d), tolerance = 1e-10)
  expect_equal(shifted$mu, base$mu * d + c(0, 5, 0, -3), tolerance = 1e-10)
  # The estimates are the maximum likelihood ones in any units: the
  # base R EM run to a far tighter criterion agrees at .001.
  oracle <- .cf_em_oracle(x * 1e-3, tol = 1e-16)
  small <- cov_fiml(x * 1e-3)
  expect_equal(small$S / 1e-6, oracle$S / 1e-6, tolerance = 1e-6)
  expect_equal(small$mu / 1e-3, oracle$mu / 1e-3, tolerance = 1e-6)
})

test_that("a warning from lavaan's EM names the variables of 'data'", {
  skip_on_cran()  # Additional model fits run in the full local suite.
  skip_if_not_installed("lavaan")
  # lavaan names the columns of the matrix its EM works on V1, V2, ...,
  # and .cov_fiml_em() hands them over as v1, v2, ...
  nm <- c("t1_visual_perception", "t2_cubes", "t4_lozenges")
  expect_identical(
    .cov_fiml_user_names(paste0("the variances of some variables are ",
                                "(near) zero: (\"V1\", \"V3\")"), nm),
    paste0("the variances of some variables are (near) zero: ",
           "(\"t1_visual_perception\", \"t4_lozenges\")"))
  expect_identical(.cov_fiml_user_names("v2 and V12 and HV1", nm),
                   "t2_cubes and V12 and HV1")
  expect_identical(.cov_fiml_user_names("no names here", nm),
                   "no names here")
  real_em <- .cov_fiml_em
  local_mocked_bindings(.cov_fiml_em = function(x, em_args = NULL) {
    warning("lavaan->lav_mvn_mi_h1_est_moments(): (near) zero: ",
            "(\"V1\", \"V2\")", call. = FALSE)
    real_em(x, em_args)
  })
  expect_warning(cov_fiml(holzinger_swineford[nm]),
                 "(\"t1_visual_perception\", \"t2_cubes\")", fixed = TRUE)
})

test_that("cov_fiml() refuses what it cannot estimate", {
  skip_on_cran()  # The extended validation grid runs in the full local suite.
  skip_if_not_installed("lavaan")
  hs <- holzinger_swineford
  two <- hs[c("t1_visual_perception", "t2_cubes")]
  expect_error(cov_fiml(list(a = 1:3)), "data frame or a numeric matrix")
  expect_error(cov_fiml(two[0, ]), "no rows")
  expect_error(cov_fiml(hs[c("t1_visual_perception", "sex")]),
               "not numeric: .*sex")
  expect_error(cov_fiml(data.frame(a = c(TRUE, FALSE, TRUE), b = 1:3)),
               "not numeric")
  expect_error(cov_fiml(data.frame(a = c(1, 1, 1, NA), b = 1:4)),
               "fewer than two distinct observed values.*a")
  expect_error(cov_fiml(data.frame(a = NA_real_, b = 1)),
               "fewer than two distinct")
  expect_error(cov_fiml(data.frame(a = c(1, Inf, 3), b = 1:3)),
               "infinite")
  never <- data.frame(a = c(1, 2, 3, NA, NA, NA, 2, 5),
                      b = c(NA, NA, NA, 1, 3, 2, NA, NA),
                      c = c(1, 3, 2, 5, 4, 6, 2, 3))
  expect_error(cov_fiml(never), "never observed together.*a with b")
  expect_error(cov_fiml(two, aux = "zz"), "not found in 'data'")
  expect_error(cov_fiml(two, aux = c("t2_cubes", "t2_cubes")),
               "duplicated")
  expect_error(cov_fiml(two, aux = 1), "character vector")
  expect_error(cov_fiml(two, aux = names(two)), "no analysis variable")
  expect_error(cov_fiml(data.frame(a = c(1, NA), b = c(NA_real_, NA))),
               "fewer than two distinct")
})

test_that("a pair never observed together is refused in the words the reliability family uses", {
  skip_on_cran()  # Additional model fits run in the full local suite.
  skip_if_not_installed("lavaan")
  # cov_fiml() and the reliability family under FIML refused such a pair
  # with two copies of the check, worded differently; both now call one
  # internal, so the messages are identical.
  never <- data.frame(a = c(1, 2, 3, NA, NA, NA, 2, 5, 4, 1),
                      b = c(NA, NA, NA, 1, 3, 2, NA, NA, NA, NA),
                      c = c(1, 3, 2, 5, 4, 6, 2, 3, 5, 2))
  from_cov <- tryCatch(cov_fiml(never), error = conditionMessage)
  from_alpha <- tryCatch(
    suppressMessages(reliability_alpha(never, missing = "fiml")),
    error = conditionMessage)
  expect_identical(from_cov, from_alpha)
  expect_identical(
    from_cov,
    paste0("These pairs of variables are never observed together, so the ",
           "data carry no information about their covariance, which full ",
           "information maximum likelihood needs: a with b. lavaan only ",
           "warns about such a pair, and its estimate of that covariance ",
           "is arbitrary, as is anything computed from it. Leave out one ",
           "variable of each pair (an auxiliary variable where the pair ",
           "has one, since it only informs the estimates); ",
           "missing_coverage() reports the cases behind every pair."))
  # A pair with an auxiliary variable is named the same way.
  expect_identical(tryCatch(cov_fiml(never, aux = "b"),
                            error = conditionMessage), from_cov)
  # So is a pair observed together only where one of them takes a single
  # value, or in one row only: those rows say nothing about how the two
  # vary together.
  flat <- data.frame(a = c(1, 2, 3, 4, NA, NA, 2, 5),
                     b = c(NA, NA, NA, 1, 3, 2, 3, 1),
                     c = c(1, 3, 2, 5, 4, 6, 2, 3))
  flat$a[7:8] <- 4
  from_cov <- tryCatch(cov_fiml(flat), error = conditionMessage)
  expect_identical(from_cov, tryCatch(
    suppressMessages(reliability_alpha(flat, missing = "fiml")),
    error = conditionMessage))
  expect_identical(
    from_cov,
    paste0("These pairs of variables are observed together only in rows ",
           "where one of them takes a single value, so those rows carry no ",
           "information about how the two vary together, which full ",
           "information maximum likelihood needs to estimate their ",
           "covariance: a with b (a takes one value in the 3 rows that ",
           "observe both). Its estimate would rest on the form of the normal ",
           "model rather than on the data, or not be determined at all, and ",
           "so would anything computed from it. Leave out one variable of ",
           "each pair (an auxiliary variable where the pair has one, since ",
           "it only informs the estimates); missing_coverage() reports the ",
           "cases behind every pair."))
  once <- flat
  once$b[7:8] <- NA
  expect_error(cov_fiml(once), "a with b (one row observes both)",
               fixed = TRUE)
  # The candidates found by the cross products are confirmed on the
  # values: a variable whose values in common differ by little is not
  # refused.
  near <- flat
  near$a[8] <- 4 + 1e-9
  expect_no_error(cov_fiml(near))
  expect_identical(nrow(DMAR:::.single_valued_pairs(near)), 0L)
})
