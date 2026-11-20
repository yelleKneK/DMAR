# The bootstraps in this file are mechanics tests run with a small B; the
# rule on the number of replications is tested in test-dmar_bootstrap.R.
local_small_B()

.bifactor_fit <- function() {
  skip_if_not_installed("lavaan")
  # Simulate clean orthogonal bifactor data so the solution is proper.
  set.seed(11)
  n <- 600
  g <- rnorm(n); f <- list(rnorm(n), rnorm(n), rnorm(n))
  x <- vapply(1:9, function(i) {
    grp <- ceiling(i / 3)
    0.5 * g + 0.5 * f[[grp]] + sqrt(1 - 0.25 - 0.25) * rnorm(n)
  }, numeric(n))
  colnames(x) <- paste0("x", 1:9)
  d <- as.data.frame(x)
  model <- paste(
    "g  =~", paste(colnames(x), collapse = " + "), "\n",
    "f1 =~ x1 + x2 + x3\n",
    "f2 =~ x4 + x5 + x6\n",
    "f3 =~ x7 + x8 + x9")
  lavaan::cfa(model, data = d, orthogonal = TRUE, std.lv = TRUE)
}

test_that("a general loading held at zero has the same eligibility as an omitted loading", {
  fit <- .bifactor_fit()
  d <- as.data.frame(lavaan::lavInspect(fit, "data"))
  group_model <- "f1 =~ x1 + x2 + x3\nf2 =~ x4 + x5 + x6\nf3 =~ x7 + x8 + x9"
  rest <- paste0("x", 2:9, collapse = " + ")
  for (general_model in c(paste("g =~", rest),
                           paste("g =~ 0*x1 +", rest),
                           paste("g =~ a*x1 +", rest, "\na == 0"))) {
    alt <- suppressWarnings(lavaan::cfa(paste(general_model, group_model, sep = "\n"),
      data = d, orthogonal = TRUE, std.lv = TRUE))
    expect_true(lavaan::lavInspect(alt, "converged"))
    expect_error(bifactor_indices(alt, general = "g"), "every item.*zero.*lavaan::cfa")
  }
})

test_that("bifactor_indices() reproduces ECV and omega_H by hand", {
  fit <- .bifactor_fit()
  res <- bifactor_indices(fit)
  expect_s3_class(res, "dmar_tbl")
  expect_equal(nrow(res), 4L)                       # general + 3 group factors
  expect_equal(res$factor[1], "g")
  expect_false(isTRUE(attr(res, "improper")))

  # Hand computation from the standardized loadings.
  std <- lavaan::standardizedSolution(fit)
  lam <- std[std$op == "=~", ]
  items <- unique(lam$rhs)
  lg <- setNames(lam$est.std[lam$lhs == "g"][match(items, lam$rhs[lam$lhs == "g"])], items)
  groups <- c("f1", "f2", "f3"); ls <- setNames(rep(NA_real_, length(items)), items)
  grp <- ls
  for (g in groups) {
    gi <- lam$rhs[lam$lhs == g]
    ls[gi] <- lam$est.std[lam$lhs == g][match(gi, lam$rhs[lam$lhs == g])]
    grp[gi] <- g
  }
  theta <- 1 - lg^2 - ls^2
  rv <- std[std$op == "~~" & std$lhs == std$rhs & std$lhs %in% items, ]
  theta[rv$lhs] <- rv$est.std

  ecv <- sum(lg^2) / (sum(lg^2) + sum(ls^2))
  total_common <- sum(lg)^2 + sum(vapply(groups, function(g) sum(ls[grp == g])^2, numeric(1)))
  total_var <- total_common + sum(theta)
  omega_h <- sum(lg)^2 / total_var

  expect_equal(res$ECV[1], ecv, tolerance = 1e-8)
  expect_equal(res$omega_H[1], omega_h, tolerance = 1e-8)

  # PUC: 9 items, 3 groups of 3 -> within pairs 3*3 = 9, total C(9,2) = 36.
  expect_equal(res$PUC[1], 1 - 9 / 36, tolerance = 1e-10)

  # Group rows carry omega_HS, not omega_H/PUC.
  expect_true(all(is.na(res$omega_H[-1])))
  expect_true(all(is.na(res$PUC[-1])))
  expect_true(all(is.finite(res$omega_HS[-1])))
})

test_that("bifactor_indices() validates and detects the general factor", {
  skip_if_not_installed("lavaan")
  fit1 <- lavaan::cfa(
    "a =~ t1_visual_perception + t2_cubes + t4_lozenges
     b =~ t6_paragraph_comprehension + t7_sentence + t9_word_meaning",
    data = holzinger_swineford, std.lv = TRUE)
  expect_error(bifactor_indices(fit1), "general factor")  # no factor loads all
  expect_error(bifactor_indices(lm(1 ~ 1)), "lavaan fit")
})

test_that("bifactor_indices() keeps the model implied denominators by default", {
  fit <- .bifactor_fit()
  res <- bifactor_indices(fit)
  res_m <- bifactor_indices(fit, denominator = "model_implied")
  expect_equal(attr(res, "denominator"), "model_implied")
  expect_equal(as.data.frame(res), as.data.frame(res_m))
})

test_that("denominator = 'observed' divides by the composite variances estimated from the data", {
  fit <- .bifactor_fit()
  res_m <- bifactor_indices(fit)
  res_o <- bifactor_indices(fit, denominator = "observed")
  expect_equal(attr(res_o, "denominator"), "observed")

  std <- lavaan::standardizedSolution(fit)
  lam <- std[std$op == "=~", ]
  items <- unique(lam$rhs)
  lg <- setNames(lam$est.std[lam$lhs == "g"][match(items, lam$rhs[lam$lhs == "g"])], items)
  R <- stats::cov2cor(unclass(lavaan::lavInspect(fit, "sampstat")$cov))[items, items]
  expect_equal(res_o$omega_H[1], sum(lg)^2 / sum(R), tolerance = 1e-8)
  f1 <- c("x1", "x2", "x3")
  ls1 <- lam$est.std[lam$lhs == "f1"][match(f1, lam$rhs[lam$lhs == "f1"])]
  expect_equal(res_o$omega_HS[res_o$factor == "f1"],
               sum(ls1)^2 / sum(R[f1, f1]), tolerance = 1e-8)
  # Indices without a denominator do not change.
  expect_equal(res_o$ECV, res_m$ECV)
  expect_equal(res_o$PUC, res_m$PUC)
  expect_equal(res_o$H, res_m$H)
  # The model is correctly specified, so the two settings agree closely.
  expect_lt(max(abs(res_o$omega_H[1] - res_m$omega_H[1])), 0.02)
  expect_output(print(res_o), "Denominator: observed composite variance")
})

# Bifactor data with three group factors of k items each, generating
# loadings of .5 on the general and the group factor. `rescov` adds a
# shared residual of variance v to the items it names.
.bifactor_data <- function(n, seed, k = 3, rescov = NULL) {
  set.seed(seed)
  J <- 3 * k
  g <- rnorm(n); f <- list(rnorm(n), rnorm(n), rnorm(n))
  e <- matrix(rnorm(n * J), n, J) * sqrt(0.5)
  if (!is.null(rescov)) {
    s <- rnorm(n) * sqrt(rescov$v)
    e[, rescov$items] <- e[, rescov$items] * sqrt((0.5 - rescov$v) / 0.5) + s
  }
  x <- vapply(seq_len(J), function(i) 0.5 * g + 0.5 * f[[ceiling(i / k)]],
              numeric(n)) + e
  colnames(x) <- paste0("x", seq_len(J))
  as.data.frame(x)
}
.bifactor_model <- function(k = 3) {
  items <- paste0("x", seq_len(3 * k))
  paste0("g =~ ", paste(items, collapse = " + "), "\n",
         paste0("f", 1:3, " =~ ",
                vapply(1:3, function(j)
                  paste(items[(j - 1) * k + seq_len(k)], collapse = " + "),
                  character(1)),
                collapse = "\n"))
}
.std_loadings <- function(fit, factor, items) {
  std <- lavaan::standardizedSolution(fit)
  lam <- std[std$op == "=~" & std$lhs == factor, ]
  lam$est.std[match(items, lam$rhs)]
}

test_that("an orthogonal fit without residual covariances keeps the earlier model implied values", {
  # The model implied denominators are now block sums of the correlation
  # matrix the fit reproduces; without residual covariances that is the
  # earlier formula to machine precision.
  fit <- .bifactor_fit()
  res <- bifactor_indices(fit)
  std <- lavaan::standardizedSolution(fit)
  items <- paste0("x", 1:9)
  lg <- .std_loadings(fit, "g", items)
  grp <- rep(c("f1", "f2", "f3"), each = 3)
  ls <- unlist(lapply(c("f1", "f2", "f3"), function(f)
    .std_loadings(fit, f, items[grp == f])))
  rv <- std[std$op == "~~" & std$lhs == std$rhs & std$lhs %in% items, ]
  theta <- rv$est.std[match(items, rv$lhs)]
  common <- sum(lg)^2 + sum(tapply(ls, grp, sum)^2)
  expect_equal(res$omega[1], common / (common + sum(theta)), tolerance = 1e-12)
  expect_equal(res$omega_H[1], sum(lg)^2 / (common + sum(theta)),
               tolerance = 1e-12)
  for (f in c("f1", "f2", "f3")) {
    sel <- grp == f
    expect_equal(res$omega_HS[res$factor == f],
                 sum(ls[sel])^2 / (sum(lg[sel])^2 + sum(ls[sel])^2 + sum(theta[sel])),
                 tolerance = 1e-12)
  }
})

test_that("a residual covariance the model specifies enters the model implied denominators", {
  # The 1.1.0 review: the denominators left out the residual covariance,
  # so a model that reproduces the data exactly showed a gap between the
  # two settings. Across groups (x2 with x5) it changes omega and
  # omega_H; within a group (x1 with x2, four items per group so the
  # model is identified) it changes that group's omega_HS.
  items <- paste0("x", 1:9)
  d <- .bifactor_data(3000, 11, rescov = list(items = c(2, 5), v = 0.25))
  fit <- lavaan::cfa(paste(.bifactor_model(), "\nx2 ~~ x5"), data = d,
                     orthogonal = TRUE, std.lv = TRUE)
  R_imp <- unclass(lavaan::lavInspect(fit, "cor.ov"))[items, items]
  lg <- .std_loadings(fit, "g", items)
  res <- bifactor_indices(fit)
  expect_equal(res$omega_H[1], sum(lg)^2 / sum(R_imp), tolerance = 1e-10)
  # The model reproduces the data, so the two settings now agree.
  res_o <- bifactor_indices(fit, denominator = "observed")
  expect_lt(abs(res$omega_H[1] - res_o$omega_H[1]), 0.002)
  # The earlier, diagonal-only denominator was visibly higher.
  std <- lavaan::standardizedSolution(fit)
  diag_only <- sum(R_imp) - 2 * std$est.std[std$op == "~~" & std$lhs == "x2" &
                                               std$rhs == "x5"]
  expect_gt(sum(lg)^2 / diag_only - res$omega_H[1], 0.005)

  items12 <- paste0("x", 1:12)
  d12 <- .bifactor_data(3000, 12, k = 4, rescov = list(items = c(1, 2), v = 0.2))
  fit12 <- lavaan::cfa(paste(.bifactor_model(4), "\nx1 ~~ x2"), data = d12,
                       orthogonal = TRUE, std.lv = TRUE)
  R12 <- unclass(lavaan::lavInspect(fit12, "cor.ov"))[items12, items12]
  f1 <- paste0("x", 1:4)
  res12 <- bifactor_indices(fit12)
  expect_equal(res12$omega_HS[res12$factor == "f1"],
               sum(.std_loadings(fit12, "f1", f1))^2 / sum(R12[f1, f1]),
               tolerance = 1e-10)
})

test_that("correlated factors are refused, since the indices assume orthogonal ones", {
  # A covariance fixed away from zero keeps the fit identified while
  # making two group factors correlate; a free one is refused the same
  # way, since the check reads the fitted factor correlations.
  d <- .bifactor_data(600, 11)
  fit <- lavaan::cfa(paste(.bifactor_model(), "\nf1 ~~ 0.3*f2"), data = d,
                     orthogonal = TRUE, std.lv = TRUE)
  expect_error(bifactor_indices(fit),
               "needs orthogonal general and group factors .* correlates f1 and f2 \\(0\\.300\\)")
  expect_error(bifactor_indices(fit, denominator = "observed"),
               "correlates f1 and f2")
})

test_that("a multiple-group fit is refused as average_variance_extracted() refuses it", {
  # Through DMAR 1.0.1 the loadings came from the first group and the
  # residual variances from the last, so no index described either group
  # (a general factor omega_H of .623 for this fit); the first builds of
  # 1.1.0 described the first group with a message. Ruling 13 of
  # 2026-11-11: refused, with the class and wording of
  # average_variance_extracted().
  skip_if_not_installed("lavaan")
  set.seed(11)
  n <- 400
  g <- rnorm(n); f <- list(rnorm(n), rnorm(n), rnorm(n))
  x <- vapply(1:9, function(i)
    0.5 * g + 0.5 * f[[ceiling(i / 3)]] + sqrt(0.5) * rnorm(n), numeric(n))
  colnames(x) <- paste0("x", 1:9)
  d <- data.frame(x, grp = rep(c("a", "b"), each = n / 2))
  fit <- lavaan::cfa(.bifactor_model(), data = d, orthogonal = TRUE,
                     std.lv = TRUE, group = "grp")
  expect_true(lavaan::lavInspect(fit, "converged"))
  ave <- expect_error(average_variance_extracted(fit),
                      class = "dmar_argument_invalid")
  # Nothing is read from the fit before the refusal, under either
  # denominator.
  local_mocked_bindings(standardizedSolution = function(...) stop("read"),
                        .package = "lavaan")
  for (den in c("model_implied", "observed")) {
    e <- expect_error(bifactor_indices(fit, denominator = den),
                      class = "dmar_argument_invalid")
    expect_s3_class(e, "dmar_refusal")
    expect_identical(e$arg, "fit")
    expect_identical(e$groups, c("a", "b"))
    expect_identical(e$groups, ave$groups)
    # The wording of average_variance_extracted(): the same opening, the
    # same reason, and the same route to one group's value.
    opening <- paste0("'fit' is a multiple-group fit (2 groups: a, b). ",
                      "Each group has its own standardized loadings")
    expect_true(startsWith(conditionMessage(e), opening))
    expect_true(startsWith(conditionMessage(ave), opening))
    for (phrase in c("neither their unweighted mean nor a mean weighted by group size",
                     "describes a fit of one group. For one group of this model, ",
                     "fit the model to that group's data alone (which drops any equality constraints across groups).")) {
      expect_match(conditionMessage(e), phrase, fixed = TRUE)
      expect_match(conditionMessage(ave), phrase, fixed = TRUE)
    }
    expect_match(conditionMessage(e), "bifactor_indices() describes",
                 fixed = TRUE)
  }
  # ci_bifactor_indices() refuses the fit as well, before any resample.
  expect_error(ci_bifactor_indices(fit, B = 400, seed = 113),
               class = "dmar_argument_invalid")
})

test_that("each group's indices are those of the model fitted to that group alone", {
  # The route the refusal names. With no equality constraints across
  # groups, the multiple-group fit is the two single-group fits side by
  # side, so the oracle for each group is omega_H formed from that group's
  # own standardized loadings and model implied correlations in the
  # multiple-group fit, (sum of general loadings)^2 / (sum of the
  # correlation matrix). The two groups differ (.621 and .674), which is
  # why no single table describes the fit.
  skip_if_not_installed("lavaan")
  set.seed(11)
  n <- 400
  g <- rnorm(n); f <- list(rnorm(n), rnorm(n), rnorm(n))
  x <- vapply(1:9, function(i)
    0.5 * g + 0.5 * f[[ceiling(i / 3)]] + sqrt(0.5) * rnorm(n), numeric(n))
  colnames(x) <- paste0("x", 1:9)
  d <- data.frame(x, grp = rep(c("a", "b"), each = n / 2))
  fit <- lavaan::cfa(.bifactor_model(), data = d, orthogonal = TRUE,
                     std.lv = TRUE, group = "grp")
  items <- paste0("x", 1:9)
  std <- lavaan::standardizedSolution(fit)
  oracle <- vapply(1:2, function(k) {
    s <- std[std$group == k & std$op == "=~" & std$lhs == "g", ]
    lg <- s$est.std[match(items, s$rhs)]
    R <- unclass(lavaan::lavInspect(fit, "cor.ov")[[k]])[items, items]
    sum(lg)^2 / sum(R)
  }, numeric(1))
  alone <- vapply(c("a", "b"), function(gr) {
    one <- lavaan::cfa(.bifactor_model(), data = d[d$grp == gr, ],
                       orthogonal = TRUE, std.lv = TRUE)
    bifactor_indices(one)$omega_H[1]
  }, numeric(1))
  expect_equal(unname(alone), oracle, tolerance = 1e-5)
  expect_equal(unname(alone), c(.6206, .6744), tolerance = 1e-3)
})

test_that("under missing = 'ml' the observed denominators are the full information estimates", {
  # Missing at random on x2, x5, and x8, driven by x1: the complete-case
  # correlations would give a visibly different denominator.
  items <- paste0("x", 1:9)
  d <- .bifactor_data(600, 11)
  set.seed(11)
  pr <- stats::plogis(-1 + 1.5 * d$x1)
  for (v in c("x2", "x5", "x8")) d[[v]][stats::runif(nrow(d)) < pr] <- NA
  fit <- lavaan::cfa(.bifactor_model(), data = d, orthogonal = TRUE,
                     std.lv = TRUE, missing = "ml")
  res_o <- bifactor_indices(fit, denominator = "observed")
  R <- stats::cov2cor(.fiml_cov(d, items)$S)[items, items]
  lg <- .std_loadings(fit, "g", items)
  f1 <- c("x1", "x2", "x3")
  expect_equal(res_o$omega_H[1], sum(lg)^2 / sum(R), tolerance = 1e-6)
  expect_equal(res_o$omega_HS[res_o$factor == "f1"],
               sum(.std_loadings(fit, "f1", f1))^2 / sum(R[f1, f1]),
               tolerance = 1e-6)
  R_cc <- stats::cor(d[stats::complete.cases(d), items])
  expect_gt(abs(res_o$omega_H[1] - sum(lg)^2 / sum(R_cc)), 0.01)
})

test_that("a fit to ordered items names the latent response metric of its denominators", {
  # lavaan's sample statistics for ordered items are the polychoric
  # correlations, and the loadings are on the latent response metric, so
  # both settings describe the sum of the latent responses, not of the
  # summed categories.
  d <- .bifactor_data(600, 11)
  X <- as.data.frame(lapply(d, function(v) as.integer(cut(v, c(-Inf, -1, 0, 1, Inf)))))
  fit <- lavaan::cfa(.bifactor_model(), data = X, ordered = TRUE,
                     orthogonal = TRUE, std.lv = TRUE)
  res_o <- bifactor_indices(fit, denominator = "observed")
  res_m <- bifactor_indices(fit)
  expect_identical(attr(res_o, "metric"), "latent_response")
  expect_identical(attr(res_m, "metric"), "latent_response")
  items <- paste0("x", 1:9)
  P <- unclass(lavaan::lavCor(X, ordered = names(X)))[items, items]
  lg <- .std_loadings(fit, "g", items)
  expect_equal(res_o$omega_H[1], sum(lg)^2 / sum(P), tolerance = 1e-8)
  expect_output(print(res_o), "Denominator: variance of the sum of the latent responses, from the polychoric correlations",
                fixed = TRUE)
  expect_output(print(res_m), "Denominator: variance of the sum of the latent responses implied by the fitted model",
                fixed = TRUE)
  # A fit to continuous items is on the observed metric.
  expect_identical(attr(bifactor_indices(.bifactor_fit()), "metric"), "observed")
})

# Nine items in three groups of three with the general loadings `lg` and
# group loadings `ls` (each a scalar or one value per item), drawn
# (`n` cases) or, with n = NULL, given as the population correlation
# matrix with a large nominal N, so the fit reproduces it exactly.
.bifactor_loadings_fit <- function(lg, ls, n = NULL, seed = 113) {
  items <- paste0("x", 1:9)
  lg <- rep_len(lg, 9L)
  ls <- rep_len(ls, 9L)
  grp <- rep(1:3, each = 3)
  model <- .bifactor_model()
  if (is.null(n)) {
    L <- cbind(lg, ls * outer(grp, 1:3, "=="))
    S <- tcrossprod(L)
    diag(S) <- 1
    dimnames(S) <- list(items, items)
    return(lavaan::cfa(model, sample.cov = S, sample.nobs = 1e5,
                       orthogonal = TRUE, std.lv = TRUE))
  }
  set.seed(seed)
  g <- rnorm(n); f <- list(rnorm(n), rnorm(n), rnorm(n))
  x <- vapply(1:9, function(i) {
    lg[i] * g + ls[i] * f[[grp[i]]] + sqrt(1 - lg[i]^2 - ls[i]^2) * rnorm(n)
  }, numeric(n))
  colnames(x) <- items
  lavaan::cfa(model, data = as.data.frame(x), orthogonal = TRUE,
              std.lv = TRUE)
}

test_that("a group factor's row gives the general factor's share of its items' common variance", {
  # With general loadings of .7 and group loadings of .3 the general
  # factor's share of a subscale's common variance (.84) and the group
  # factor's own share (.16) are far apart, and the table carries the
  # first; with equal loadings the two coincide and a test cannot tell
  # them apart.
  fit <- .bifactor_loadings_fit(0.7, 0.3, n = 600)
  res <- bifactor_indices(fit)
  items <- paste0("x", 1:9)
  R <- unclass(lavaan::lavInspect(fit, "cor.ov"))[items, items]
  lg <- .std_loadings(fit, "g", items)
  for (j in 1:3) {
    s <- items[(j - 1) * 3 + 1:3]
    lg_s <- lg[(j - 1) * 3 + 1:3]
    ls_s <- .std_loadings(fit, paste0("f", j), s)
    row <- res[res$factor == paste0("f", j), ]
    general_share <- sum(lg_s^2) / sum(lg_s^2 + ls_s^2)
    expect_equal(row$ECV, general_share, tolerance = 1e-12)
    expect_gt(row$ECV - (1 - general_share), 0.5)
    # omega: the subscale total score, both factors counted as true score.
    expect_equal(row$omega, (sum(lg_s)^2 + sum(ls_s)^2) / sum(R[s, s]),
                 tolerance = 1e-12)
    expect_equal(row$omega_HS, sum(ls_s)^2 / sum(R[s, s]), tolerance = 1e-12)
    t_s <- sum(ls_s^2 / (1 - ls_s^2))
    expect_equal(row$H, t_s / (1 + t_s), tolerance = 1e-12)
  }
  foot <- DMAR:::.dmar_tbl_footer(res)
  expect_match(foot[names(foot) == "bifactor_definitions"][1],
               "general factor's share of the common variance of that factor's items")
})

test_that("H is the index of Rodriguez, Reise, and Haviland, and the maximal reliability is beside it", {
  # The index counts the other factors' common variance as uncorrelated
  # error. In the population with every loading .5 it is 3/4 for the
  # general factor, whose maximal reliability is 9/14; with the group
  # loadings .5, -.5, .5 in each subscale the index is still 3/4 but the
  # optimal weights cancel group variance, and the maximal reliability is
  # 21/26.
  same <- bifactor_indices(.bifactor_loadings_fit(0.5, 0.5))
  mixed <- bifactor_indices(.bifactor_loadings_fit(0.5, c(0.5, -0.5, 0.5)))
  expect_equal(same$H[1], 3 / 4, tolerance = 1e-6)
  expect_equal(mixed$H[1], 3 / 4, tolerance = 1e-6)
  max_same <- attr(same, "maximal_reliability")
  max_mixed <- attr(mixed, "maximal_reliability")
  expect_identical(names(max_same), c("g", "f1", "f2", "f3"))
  expect_equal(unname(max_same[1]), 9 / 14, tolerance = 1e-6)
  expect_equal(unname(max_mixed[1]), 21 / 26, tolerance = 1e-6)
  expect_gt(same$H[1], max_same[[1]])
  expect_lt(mixed$H[1], max_mixed[[1]])

  # On a sample it is l' R^{-1} l with the model implied correlations,
  # the square of the regression score determinacy factor_scores()
  # reports, and it does not depend on the denominator.
  fit <- .bifactor_fit()
  res <- bifactor_indices(fit)
  items <- paste0("x", 1:9)
  R <- unclass(lavaan::lavInspect(fit, "cor.ov"))[items, items]
  L <- cbind(g = .std_loadings(fit, "g", items),
             f1 = c(.std_loadings(fit, "f1", items[1:3]), rep(0, 6)),
             f2 = c(rep(0, 3), .std_loadings(fit, "f2", items[4:6]), rep(0, 3)),
             f3 = c(rep(0, 6), .std_loadings(fit, "f3", items[7:9])))
  expect_equal(attr(res, "maximal_reliability"),
               colSums(L * solve(R, L)), tolerance = 1e-10)
  det <- attr(factor_scores(fit), "determinacy")
  expect_equal(unname(attr(res, "maximal_reliability")),
               det$determinacy_regression[match(c("g", "f1", "f2", "f3"),
                                                det$factor)]^2,
               tolerance = 1e-6)
  expect_identical(attr(bifactor_indices(fit, denominator = "observed"),
                        "maximal_reliability"),
                   attr(res, "maximal_reliability"))
  foot <- DMAR:::.dmar_tbl_footer(res)
  expect_match(foot[names(foot) == "bifactor_definitions"][2],
               "Rodriguez, Reise, and Haviland \\(2016\\).*g 0\\.[0-9]{3}, f1")
  expect_output(print(res), "maximal reliability the fitted model implies")
})

test_that("with a strong general factor a group factor's maximal reliability is above its index, though the general factor's is below", {
  # The help page's direction rule holds row by row. On the general
  # factor's row, with positive and equal group loadings, the index is
  # the larger. On a group factor's row the optimal weights remove the
  # general factor through the other subscales' items, so with general
  # loadings of .7 and group loadings of .3 (ECV .84) the maximal
  # reliability is the larger, 0.279 against 0.229; with loadings of .5
  # and .5 the index is the larger on every row.
  strong <- bifactor_indices(.bifactor_loadings_fit(0.7, 0.3))
  even <- bifactor_indices(.bifactor_loadings_fit(0.5, 0.5))
  index_H <- function(l) sum(l^2 / (1 - l^2)) / (1 + sum(l^2 / (1 - l^2)))
  # The population values, from the loadings alone.
  items <- paste0("x", 1:9)
  L <- cbind(0.7, 0.3 * outer(rep(1:3, each = 3), 1:3, "=="))
  R <- tcrossprod(L)
  diag(R) <- 1
  maximal <- colSums(L * solve(R, L))
  expect_equal(strong$H, c(index_H(rep(0.7, 9)), rep(index_H(rep(0.3, 3)), 3)),
               tolerance = 1e-6)
  expect_equal(unname(attr(strong, "maximal_reliability")), maximal,
               tolerance = 1e-6)
  expect_equal(strong$H[1:2], c(0.8963, 0.2288), tolerance = 1e-4)
  expect_equal(unname(maximal[1:2]), c(0.8647, 0.2785), tolerance = 1e-4)
  expect_gt(strong$H[1], attr(strong, "maximal_reliability")[["g"]])
  for (f in c("f1", "f2", "f3")) {
    expect_lt(strong$H[strong$factor == f],
              attr(strong, "maximal_reliability")[[f]])
    expect_gt(even$H[even$factor == f],
              attr(even, "maximal_reliability")[[f]])
  }
  # The group factor's maximal reliability belongs to its regression
  # score from all nine items, which weights the six items of the other
  # subscales negatively; the best weighting of the subscale's own three
  # items reaches only its omega_HS here (equal loadings), 0.125, below
  # the index as well as the maximal reliability.
  w <- solve(R, L[, 2])
  expect_equal(w / max(w), rep(c(1, -0.405), c(3, 6)), tolerance = 1e-3)
  own <- sum(L[1:3, 2] * solve(R[1:3, 1:3], L[1:3, 2]))
  expect_equal(own, 0.125, tolerance = 1e-6)
  expect_equal(strong$omega_HS[strong$factor == "f1"], own, tolerance = 1e-6)
  expect_lt(own, strong$H[strong$factor == "f1"])
  # The better the general factor is measured outside a subscale, the
  # closer the group factor's maximal reliability comes to the index's
  # closed form with the unique variance 1 - .7^2 - .3^2 as each item's
  # error (0.391): with 20 subscales of three it is 0.372.
  G <- 20
  L20 <- cbind(0.7, 0.3 * outer(rep(seq_len(G), each = 3), seq_len(G), "=="))
  R20 <- tcrossprod(L20)
  diag(R20) <- 1
  q <- 3 * 0.3^2 / (1 - 0.7^2 - 0.3^2)
  expect_equal(sum(L20[, 2] * solve(R20, L20[, 2])), 0.3722, tolerance = 1e-3)
  expect_lt(sum(L20[, 2] * solve(R20, L20[, 2])), q / (1 + q))
})

test_that("an item on two group factors is refused, and a loading fixed or constrained to zero does not count", {
  # Through 1.0.1 the group factor listed last took the item, and its
  # other loading dropped out of ECV, omega, omega_HS, H, and PUC: here
  # ECV was .611 where every loading counted gives .594.
  set.seed(113)
  n <- 800
  g <- rnorm(n); f <- list(rnorm(n), rnorm(n), rnorm(n))
  X <- vapply(1:9, function(i) {
    0.6 * g + 0.5 * f[[ceiling(i / 3)]] + sqrt(1 - 0.36 - 0.25) * rnorm(n)
  }, numeric(n))
  X[, 3] <- 0.6 * g + 0.45 * f[[1]] + 0.45 * f[[2]] +
    sqrt(1 - 0.36 - 2 * 0.2025) * rnorm(n)
  colnames(X) <- paste0("x", 1:9)
  X <- as.data.frame(X)
  general <- "g =~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9\n"
  two <- lavaan::cfa(paste0(general, "f1 =~ x1 + x2 + x3\n",
                            "f2 =~ x3 + x4 + x5 + x6\nf3 =~ x7 + x8 + x9"),
                     data = X, orthogonal = TRUE, std.lv = TRUE)
  expect_error(bifactor_indices(two),
               "one group factor at most, but x3 loads on f1 and f2")
  # The same model written with x3's f2 loading fixed at zero is the
  # model without it, and gets its indices.
  zero <- lavaan::cfa(paste0(general, "f1 =~ x1 + x2 + x3\n",
                             "f2 =~ 0*x3 + x4 + x5 + x6\nf3 =~ x7 + x8 + x9"),
                      data = X, orthogonal = TRUE, std.lv = TRUE)
  without <- lavaan::cfa(paste0(general, "f1 =~ x1 + x2 + x3\n",
                                "f2 =~ x4 + x5 + x6\nf3 =~ x7 + x8 + x9"),
                         data = X, orthogonal = TRUE, std.lv = TRUE)
  a <- bifactor_indices(zero)
  b <- bifactor_indices(without)
  expect_equal(as.data.frame(a)[, -1], as.data.frame(b)[, -1],
               tolerance = 1e-5)
  expect_equal(a$PUC[1], 1 - 9 / 36)
  # So is the loading constrained to zero by an equality constraint, which
  # lavaan estimates as a free parameter; read as a loading, it made x3
  # an item of two group factors, and the fit was refused.
  constrained <- lavaan::cfa(
    paste0(general, "f1 =~ x1 + x2 + x3\n",
           "f2 =~ a*x3 + x4 + x5 + x6\nf3 =~ x7 + x8 + x9\na == 0"),
    data = X, orthogonal = TRUE, std.lv = TRUE)
  e <- bifactor_indices(constrained)
  expect_equal(as.data.frame(e)[, -1], as.data.frame(b)[, -1],
               tolerance = 1e-5)
  expect_equal(attr(e, "maximal_reliability"), attr(b, "maximal_reliability"),
               tolerance = 1e-5)
  # An item on no group factor is on the general factor alone (four items
  # per group, so the group factor that loses it is still identified).
  # Through DMAR 1.0.x and 1.1.0 development it was refused; it now counts
  # in the general factor's indices and in PUC, whose pairs with it are
  # uncontaminated, and in no subscale, as fa_bifactor() has it.
  items12 <- paste0("x", 1:12)
  none <- lavaan::cfa(
    paste0("g =~ ", paste(items12, collapse = " + "), "\n",
           "f1 =~ x1 + x2 + x3\nf2 =~ x5 + x6 + x7 + x8\n",
           "f3 =~ x9 + x10 + x11 + x12"),
    data = .bifactor_data(600, 11, k = 4), orthogonal = TRUE,
    std.lv = TRUE)
  alone <- bifactor_indices(none)
  expect_identical(alone$factor, c("g", "f1", "f2", "f3"))
  expect_equal(alone$PUC[1L], 1 - (3 + 6 + 6) / 66)
  std <- lavaan::standardizedSolution(none)
  lg <- std$est.std[std$op == "=~" & std$lhs == "g"]
  ls <- std$est.std[std$op == "=~" & std$lhs != "g"]
  expect_equal(alone$ECV[1L], sum(lg^2) / (sum(lg^2) + sum(ls^2)),
               tolerance = 1e-10)
  expect_true(all(is.finite(attr(alone, "maximal_reliability"))))
  b4 <- fa_bifactor(.bifactor_data(600, 11, k = 4),
                    list(f1 = c("x1", "x2", "x3"),
                         f2 = c("x5", "x6", "x7", "x8"),
                         f3 = c("x9", "x10", "x11", "x12")),
                    general = "g", general_only = "x4",
                    missing = "listwise", output = "measurement")
  expect_equal(b4$estimate[match(paste0("ECV_", alone$factor), b4$term)],
               alone$ECV, tolerance = 1e-5)
  expect_equal(b4$estimate[match(paste0("H_", alone$factor), b4$term)],
               alone$H, tolerance = 1e-5)
  # A fit with one group factor is refused, and the message names the
  # items on no group factor.
  one <- lavaan::cfa(
    paste0("g =~ ", paste(items12, collapse = " + "), "\n",
           "f1 =~ x1 + x2 + x3 + x4"),
    data = .bifactor_data(600, 11, k = 4), orthogonal = TRUE,
    std.lv = TRUE)
  expect_error(bifactor_indices(one),
               "this fit has one \\(f1\\), with items on no group factor \\(x5,")
  # ci_bifactor_indices() refuses the fit with the same message.
  skip_if_not_installed("boot")
  expect_error(ci_bifactor_indices(two, B = 10),
               "x3 loads on f1 and f2")
})

test_that("a fit with one group factor is refused, by bifactor_indices() and ci_bifactor_indices() alike", {
  # A general factor beside one group factor, every item listed under it
  # (the loadings of x4 to x9 fixed at zero), which DMAR 1.0.0 described.
  # It is not a bifactor model with subscales.
  d <- .bifactor_data(600, 11)
  items <- paste0("x", 1:9)
  model <- paste0("g =~ ", paste(items, collapse = " + "), "\n",
                  "f1 =~ x1 + x2 + x3 + ",
                  paste0("0*", items[4:9], collapse = " + "))
  fit <- lavaan::cfa(model, data = d, orthogonal = TRUE, std.lv = TRUE)
  msg <- paste0("this fit has one \\(f1\\)\\. A general factor \\(g\\) ",
                "beside a single group factor is not a bifactor model ",
                "with subscales.*two or more group factors")
  expect_error(bifactor_indices(fit, general = "g"), msg)
  expect_no_warning(expect_error(bifactor_indices(fit, general = "g"), msg))
  skip_if_not_installed("boot")
  expect_error(ci_bifactor_indices(fit, general = "g", B = 10), msg)
})

test_that("H is not defined for a factor with a standardized loading past one", {
  # An improper solution in which x1's general loading is 1.06. Through
  # 1.0.1 H left that loading out and reported the H of the other eight
  # general loadings, 0.841, lower than H at any proper loading of x1
  # from .9 up.
  set.seed(7)
  n <- 150
  lg <- c(0.9, rep(0.6, 8)); ls <- c(0.42, rep(0.5, 8))
  g <- rnorm(n); f <- list(rnorm(n), rnorm(n), rnorm(n))
  X <- vapply(1:9, function(i) {
    lg[i] * g + ls[i] * f[[ceiling(i / 3)]] +
      sqrt(1 - lg[i]^2 - ls[i]^2) * rnorm(n)
  }, numeric(n))
  colnames(X) <- paste0("x", 1:9)
  fit <- suppressWarnings(lavaan::cfa(.bifactor_model(),
                                      data = as.data.frame(X),
                                      orthogonal = TRUE, std.lv = TRUE))
  items <- paste0("x", 1:9)
  l_g <- .std_loadings(fit, "g", items)
  expect_gt(abs(l_g[1]), 1)
  expect_warning(res <- bifactor_indices(fit), "Improper \\(Heywood\\)")
  expect_true(attr(res, "improper"))
  expect_true(is.na(res$H[1]))
  # The group factors' loadings are below one, so their H stands.
  t1 <- sum(.std_loadings(fit, "f1", items[1:3])^2 /
              (1 - .std_loadings(fit, "f1", items[1:3])^2))
  expect_equal(res$H[2], t1 / (1 + t1), tolerance = 1e-12)
  expect_true(all(is.na(attr(res, "maximal_reliability"))))
  foot <- DMAR:::.dmar_tbl_footer(res)
  expect_match(foot[names(foot) == "bifactor_definitions"][2],
               "not defined for an improper solution")
  # Kept as it is, the point estimate gives no interval.
  skip_if_not_installed("boot")
  expect_error(suppressWarnings(
    ci_bifactor_indices(fit, B = 10, improper = "keep")),
    "not defined or not a proportion")
})


test_that("supplied two-stage auxiliary models cannot silently lose auxiliaries in bootstrap refits", {
  fit <- .bifactor_fit()
  d <- as.data.frame(lavaan::lavInspect(fit, "data"))
  set.seed(73)
  d$a <- rnorm(nrow(d))
  model <- paste("g =~", paste0("x", 1:9, collapse = " + "),
    "\nf1 =~ x1 + x2 + x3\nf2 =~ x4 + x5 + x6\nf3 =~ x7 + x8 + x9")
  aux_fit <- suppressWarnings(lavaan::cfa(model, data = d, aux = "a",
    missing = "two.stage", orthogonal = TRUE, std.lv = TRUE))
  expect_null(.fiml_fit_h1_diagnostic(aux_fit))
  expect_true(lavaan::lavInspect(aux_fit, "converged"))
  h_fit <- suppressWarnings(lavaan::cfa(
    "f1 =~ x1 + x2 + x3\nf2 =~ x4 + x5 + x6\nf3 =~ x7 + x8 + x9",
    data = d, aux = "a", missing = "two.stage", std.lv = TRUE))
  for (f in list(
    function() ci_bifactor_indices(aux_fit, B = 10),
    function() reliability_H(h_fit, factor = "f1", ci_method = "percentile", B = 10),
    function() average_variance_extracted(aux_fit, ci_method = "percentile", B = 100))) {
    expect_error(f(), "cannot yet carry separately stored two-stage auxiliary", fixed = TRUE)
  }
})

test_that("a fit whose optimizer did not converge is refused, under either denominator", {
  fit <- .bifactor_fit()
  d <- as.data.frame(lavaan::lavInspect(fit, "data"))
  # Through DMAR 1.0.1 and the first builds of 1.1.0 a fit stopped after one
  # iteration gave a full table of indices; for the nine items of the
  # deep review fixture the general factor's ECV, omega, and H were .558,
  # .851, and .764, against .610, .810, and .743 for the converged fit.
  # lavaan's own warning that the optimizer stopped is the point of it.
  bad <- suppressWarnings(lavaan::cfa(.bifactor_model(), data = d,
    orthogonal = TRUE, std.lv = TRUE, control = list(iter.max = 1L)))
  expect_false(lavaan::lavInspect(bad, "converged"))
  for (den in c("model_implied", "observed")) {
    e <- expect_error(bifactor_indices(bad, denominator = den),
                      class = "dmar_not_converged")
    expect_identical(e$fn, "bifactor_indices()")
    expect_identical(e$iterations, 1L)
  }
  # ci_bifactor_indices() refuses it with the same class and helper.
  skip_if_not_installed("boot")
  e <- expect_error(ci_bifactor_indices(bad, B = 10),
                    class = "dmar_not_converged")
  expect_identical(e$fn, "ci_bifactor_indices()")
  # The converged fit of the same model is read as before: its general
  # factor's ECV is the oracle computed from the standardized loadings
  # of lavaan::lavInspect(fit, "std").
  lam <- lavaan::lavInspect(fit, "std")$lambda
  ecv <- sum(lam[, "g"]^2) / sum(lam^2)
  expect_equal(bifactor_indices(fit)$ECV[1], ecv, tolerance = 1e-10)
})

test_that("a multilevel fit is refused rather than read across its levels", {
  skip_on_cran()  # a two-level bifactor fit
  skip_if_not_installed("lavaan")
  set.seed(113)
  J <- 200L; m <- 8L; cl <- rep(seq_len(J), each = m); n <- J * m
  gw <- rnorm(n); sw <- matrix(rnorm(n * 3), n, 3); gb <- rnorm(J)[cl]
  X <- vapply(1:9, function(j) .5 * gw + .45 * sw[, ceiling(j / 3)] +
                .4 * gb + rnorm(J, 0, .2)[cl] + rnorm(n, 0, .6), numeric(n))
  colnames(X) <- paste0("x", 1:9)
  d <- data.frame(X, cluster = cl)
  # Through DMAR 1.0.1 the table of this fit read the within level's
  # loadings with the between level's residual variances (a general factor
  # omega of .942 where the within level's own solution gives .885), and
  # the first builds of 1.1.0 stopped with R's "incorrect number of
  # dimensions".
  fit <- lavaan::cfa(paste0("level: 1\n", .bifactor_model(),
                            "\nlevel: 2\n g =~ ",
                            paste0("x", 1:9, collapse = " + ")),
                     data = d, cluster = "cluster", std.lv = TRUE,
                     orthogonal = TRUE)
  expect_true(lavaan::lavInspect(fit, "converged"))
  e <- expect_error(bifactor_indices(fit), class = "dmar_argument_invalid")
  expect_identical(e$arg, "fit")
  expect_match(conditionMessage(e), "multilevel fit (2 levels: within, cluster)",
               fixed = TRUE)
})
