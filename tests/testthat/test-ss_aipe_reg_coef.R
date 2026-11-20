test_that("ss_aipe_reg_coef() returns a positive integer sample size", {
  result <- ss_aipe_reg_coef(true_R2_Y_X = 0.3, true_R2_j_X_without_j = 0.2, p = 3, true_reg_coef_j = 0.4, width = 0.3)
  expect_s3_class(result, "data.frame")
  expect_named(result, c("term", "value"))
  expect_gt(result$value, 0)
})

test_that("ss_aipe_reg_coef() requires a larger n for a narrower width", {
  n_wide   <- ss_aipe_reg_coef(true_R2_Y_X = 0.3, true_R2_j_X_without_j = 0.2, p = 3, true_reg_coef_j = 0.4, width = 0.40)$value
  n_narrow <- ss_aipe_reg_coef(true_R2_Y_X = 0.3, true_R2_j_X_without_j = 0.2, p = 3, true_reg_coef_j = 0.4, width = 0.20)$value
  expect_lt(n_wide, n_narrow)
})

# One predictor ------------------------------------------------------------

# The expected-width criterion of ss_aipe_reg_coef() (Kelley & Maxwell,
# 2003, Eq. 1 with the expected squared multiple correlations in place of
# the population ones), coded here without the package: the Gauss
# hypergeometric 2F1(1, 1; c; z) by its series, whose terms are
# k! z^k / (c)_k. The squared multiple correlation of X_j with the other
# predictors is that of a regression on p - 1 predictors (0 with p = 1).
.hyp2F1_11 <- function(c, z) {
  term <- 1
  total <- 1
  k <- 0
  repeat {
    term <- term * (k + 1) / (c + k) * z
    total <- total + term
    k <- k + 1
    if (abs(term) < 1e-16 * abs(total)) break
  }
  total
}
.expected_R2 <- function(R2, N, p) {
  max(0, 1 - (N - p - 1) / (N - 1) * (1 - R2) * .hyp2F1_11((N + 1) / 2, R2))
}
.planner_width <- function(N, R2, R2_j, p, conf_level = .95, sd_ratio = 1,
                           tails = rep((1 - conf_level) / 2, 2)) {
  (qt(1 - tails[1], N - p - 1) + qt(1 - tails[2], N - p - 1)) * sd_ratio *
    sqrt((1 - .expected_R2(R2, N, p)) /
           ((1 - .expected_R2(R2_j, N, p - 1)) * (N - p - 1)))
}

test_that("ss_aipe_reg_coef() plans for one predictor on the covariance structure path as on the direct path", {
  skip_on_cran()
  r <- .4
  for (width in c(.15, .2, .3)) {
    direct <- ss_aipe_reg_coef(true_R2_Y_X = r^2, true_R2_j_X_without_j = 0,
                               p = 1, true_reg_coef_j = r, width = width)$value
    expect_equal(ss_aipe_reg_coef(p = 1, true_cor_XX = matrix(1), true_cor_YX = r,
                                  which_predictor = 1, width = width)$value, direct)
    # p and which_predictor follow from a 1 x 1 matrix.
    expect_equal(ss_aipe_reg_coef(true_cor_XX = matrix(1), true_cor_YX = r,
                                  width = width)$value, direct)
    expect_equal(ss_aipe_reg_coef(p = 1, true_cor_XX = matrix(1), true_cor_YX = r,
                                  which_predictor = 1, true_reg_coef_j = r,
                                  width = width)$value, direct)
    # The answer is the smallest N whose expected width meets the target,
    # by the criterion coded above.
    expect_lte(.planner_width(direct, r^2, 0, 1), width)
    expect_gt(.planner_width(direct - 1, r^2, 0, 1), width)
    # 327 at width .2, the plan NEWS quotes (328 in the development
    # versions, with the expectation of X_j on one predictor).
    if (width == .2) expect_identical(direct, 327)
  }
  # With assurance and on the noncentral path the two routes agree too.
  expect_equal(
    ss_aipe_reg_coef(p = 1, true_cor_XX = matrix(1), true_cor_YX = r,
                     which_predictor = 1, width = .2, assurance = .9)$value,
    ss_aipe_reg_coef(true_R2_Y_X = r^2, true_R2_j_X_without_j = 0, p = 1,
                     true_reg_coef_j = r, width = .2, assurance = .9)$value)
  expect_equal(
    ss_aipe_reg_coef(p = 1, true_cor_XX = matrix(1), true_cor_YX = r,
                     which_predictor = 1, width = .2, noncentral = TRUE)$value,
    ss_aipe_reg_coef(true_R2_Y_X = r^2, true_R2_j_X_without_j = 0, p = 1,
                     true_reg_coef_j = r, width = .2, noncentral = TRUE)$value)
  # A negative correlation plans the same N.
  expect_equal(ss_aipe_reg_coef(true_cor_XX = matrix(1), true_cor_YX = -r,
                                width = .2)$value,
               ss_aipe_reg_coef(true_cor_XX = matrix(1), true_cor_YX = r,
                                width = .2)$value)
})

test_that("ss_aipe_reg_coef() with p = 1 refuses a structure of the wrong size by name", {
  e <- tryCatch(ss_aipe_reg_coef(p = 1, true_cor_XX = diag(2), true_cor_YX = c(.4, .3),
                                 which_predictor = 1, width = .2),
                error = function(e) e)
  expect_s3_class(e, "dmar_argument_invalid")
  expect_identical(e$arg, "true_cor_XX")
  e <- tryCatch(ss_aipe_reg_coef(p = 1, true_cor_XX = matrix(1), true_cor_YX = c(.4, .3),
                                 which_predictor = 1, width = .2),
                error = function(e) e)
  expect_s3_class(e, "dmar_argument_invalid")
  expect_identical(e$arg, "true_cor_YX")
  e <- tryCatch(ss_aipe_reg_coef(p = 1, true_R2_Y_X = .16, true_R2_j_X_without_j = .1,
                                 true_reg_coef_j = .4, width = .2),
                error = function(e) e)
  expect_s3_class(e, "dmar_argument_invalid")
  expect_identical(e$arg, "true_R2_j_X_without_j")
})

# The scale of true_reg_coef_j --------------------------------------------

test_that("ss_aipe_reg_coef() checks true_reg_coef_j on the scale true_sd_Y and true_sd_X describe", {
  # Two predictors with covariances var(Y) = 4, cov(Y, X) = (.8, .6), and
  # unit predictor variances correlated .3. By hand the unstandardized
  # coefficient of X1 is the first element of solve(cov_XX, cov_YX), .681319,
  # and the standardized one is that times sd_X1 / sd_Y = 1 / 2, .340659
  # (Kelley & Maxwell, 2003, p. 316).
  cov_XX <- matrix(c(1, .3, .3, 1), 2)
  cov_YX <- c(.8, .6)
  b <- solve(cov_XX, cov_YX)[1]
  beta <- b / 2
  expect_equal(b, .62 / .91)
  cor_YX <- cov_YX / (2 * 1)
  plan <- function(...) ss_aipe_reg_coef(width = .3, true_cor_XX = cov_XX,
                                         true_cor_YX = cor_YX, which_predictor = 1, ...)
  # The unstandardized coefficient with the standard deviations that make
  # it unstandardized is accepted, and the plan is that of the same call
  # without it.
  expect_equal(plan(true_sd_Y = 2, true_sd_X = 1, true_reg_coef_j = b)$value,
               plan(true_sd_Y = 2, true_sd_X = 1)$value)
  # The standardized value given as the unstandardized one is refused.
  e <- tryCatch(plan(true_sd_Y = 2, true_sd_X = 1, true_reg_coef_j = beta),
                error = function(e) e)
  expect_s3_class(e, "dmar_argument_invalid")
  expect_identical(e$arg, "true_reg_coef_j")
  expect_equal(e$implied, b)
  # With unit standard deviations the coefficient is the standardized one.
  expect_equal(plan(true_reg_coef_j = beta)$value, plan()$value)
  expect_error(plan(true_reg_coef_j = b), class = "dmar_argument_invalid")
  # The widths: the plan for the unstandardized coefficient is the plan for
  # the standardized one at a width scaled by sd_X / sd_Y.
  expect_equal(plan(true_sd_Y = 2, true_sd_X = 1)$value,
               ss_aipe_reg_coef(width = .15, true_cor_XX = cov_XX,
                                true_cor_YX = cor_YX, which_predictor = 1)$value)
})

# ss_aipe_src() in its own names ------------------------------------------

test_that("ss_aipe_src() states the engine's conditions in its own argument names", {
  skip_on_cran()
  e <- tryCatch(ss_aipe_src(true_R2_Y_X = .3, true_R2_j_X_without_j = .2, p = 3,
                            width = .2),
                error = function(e) e)
  expect_s3_class(e, "dmar_argument_invalid")
  expect_match(conditionMessage(e), "'true_std_reg_coef_j' must also be specified",
               fixed = TRUE)
  expect_false(grepl("(^|[^_])true_reg_coef_j", conditionMessage(e)))
  expect_identical(e$arg, "true_std_reg_coef_j")

  cov_XX <- matrix(c(1, .3, .3, 1), 2)
  e <- tryCatch(ss_aipe_src(true_cor_XX = cov_XX, true_cor_YX = c(.4, .3),
                            which_predictor = 1, true_std_reg_coef_j = .9,
                            width = .2),
                error = function(e) e)
  expect_s3_class(e, "dmar_argument_invalid")
  expect_identical(e$arg, "true_std_reg_coef_j")
  expect_match(conditionMessage(e), "'true_std_reg_coef_j' (0.9)", fixed = TRUE)
  expect_false(grepl("(^|[^_])true_reg_coef_j", conditionMessage(e)))
  # The implied standardized coefficient is accepted.
  beta <- solve(cov_XX, c(.4, .3))[1]
  expect_equal(ss_aipe_src(true_cor_XX = cov_XX, true_cor_YX = c(.4, .3),
                           which_predictor = 1, true_std_reg_coef_j = beta,
                           width = .2)$value,
               ss_aipe_src(true_cor_XX = cov_XX, true_cor_YX = c(.4, .3),
                           which_predictor = 1, width = .2)$value)

  # A standard deviation other than 1 is refused by ss_aipe_src() itself,
  # not by the noncentral path of the engine, which named an argument
  # ('noncentral') the user never gave.
  for (nm in c("true_sd_Y", "true_sd_X")) {
    args <- list(true_cor_XX = cov_XX, true_cor_YX = c(.4, .3),
                 which_predictor = 1, width = .2)
    args[[nm]] <- 2
    e <- tryCatch(do.call(ss_aipe_src, args), error = function(e) e)
    expect_s3_class(e, "dmar_argument_invalid")
    expect_identical(e$arg, nm)
    expect_false(grepl("noncentral", conditionMessage(e)))
  }
})

# The assurance step (assurance_method) ---------------------------------

# The probe and the examples of Kelley and Maxwell (2003, pp. 312--313).
.probe_R <- matrix(c(1, .3, .3, 1), 2)
.probe_r <- c(.4, .3)
.exch <- function(p, rxx) { R <- matrix(rxx, p, p); diag(R) <- 1; R }
.km_R3 <- rbind(c(1, .4, .6), c(.4, 1, .05), c(.6, .05, 1))
.km_r3 <- c(.5, .3, .1)

# The probability that the central interval for the coefficient in the
# units of the data is no wider than `width`, coded here without the
# package and without pf(): the residual sum of squares of Y is
# sd_Y^2 (1 - R2) times a chi square on N - p - 1 degrees of freedom,
# independent of the predictors, and the residual sum of squares of X_j on
# the others is sd_X^2 (1 - R2_j) times a chi square V on N - p, so the
# probability is the expectation over V of a chi square probability,
# computed by numerical integration.
.raw_prob_integral <- function(N, p, R2, R2_j, sd_ratio, width,
                               a_lower = .025, a_upper = .025) {
  nu <- N - p - 1
  q <- qt(1 - a_upper, nu) + qt(1 - a_lower, nu)
  ratio <- (1 - R2_j) / (sd_ratio^2 * (1 - R2))
  f <- function(v) pchisq(width^2 * nu * ratio * v / q^2, nu) * dchisq(v, N - p)
  integrate(f, qchisq(1e-15, N - p), qchisq(1 - 1e-15, N - p),
            rel.tol = 1e-12, subdivisions = 2000L)$value
}
.R2s <- function(R, r, j) {
  c(R2 = drop(r %*% solve(R, r)),
    R2_j = if (ncol(R) > 1) 1 - 1 / solve(R)[j, j] else 0)
}
.val <- function(x, term) x$value[x$term == term]
# The chi square step of Kelley and Maxwell (2003, Equation 3), the plan of
# DMAR 1.0.x with assurance for the central interval, which no argument of
# ss_aipe_reg_coef() gives now: the engine's analytic approximation, the
# one ss_aipe_src(noncentral = FALSE) uses at unit standard deviations.
.chisq_plan <- function(...) {
  a <- list(true_R2_Y_X = NULL, true_R2_j_X_without_j = NULL, p = NULL,
            true_reg_coef_j = NULL, which_width = "Full", true_sd_Y = 1,
            true_sd_X = 1, true_cor_XX = NULL, true_cor_YX = NULL,
            which_predictor = NULL, noncentral = FALSE, alpha_lower = NULL,
            alpha_upper = NULL, conf_level = .95, assurance = NULL,
            step = "approximation", G = NULL)
  given <- list(...)
  a[names(given)] <- given
  do.call(DMAR:::.ss_aipe_reg_coef_plan, a)
}

test_that("the analytic assurance step for the central interval returns the smallest N that attains the request", {
  sets <- list(
    list(R = .probe_R, r = .probe_r, j = 1, width = .3, a = .9, N = 187),
    list(R = .exch(5, .4), r = rep(.3, 5), j = 1, width = .2, a = .8, N = 491),
    list(R = .km_R3, r = .km_r3, j = 3, width = .3, a = .9, N = 238),
    list(R = matrix(1), r = .4, j = 1, width = .2, a = .99, N = 410))
  for (s in sets) {
    plan <- ss_aipe_reg_coef(
      width = s$width, true_cor_XX = s$R, true_cor_YX = s$r,
      which_predictor = s$j, assurance = s$a)
    expect_identical(plan$term, c("necessary_N", "assurance", "assurance_attained"))
    expect_identical(attr(plan, "assurance_method"), "analytic")
    expect_identical(attr(plan, "assurance_probability"), "exact")
    N <- .val(plan, "necessary_N")
    expect_equal(N, s$N)
    rr <- .R2s(s$R, s$r, s$j)
    p <- ncol(s$R)
    oracle <- function(n) .raw_prob_integral(n, p, rr[["R2"]], rr[["R2_j"]], 1, s$width)
    # The probability reported is the probability at N, by the integral.
    expect_equal(.val(plan, "assurance_attained"), oracle(N), tolerance = 1e-8)
    # Minimal: the request is met at N and not at N - 1, and N is not
    # below the plan for the expected width.
    expect_gte(oracle(N), s$a)
    expect_lt(oracle(N - 1), s$a)
    expect_gte(N, ss_aipe_reg_coef(width = s$width, true_cor_XX = s$R,
                                   true_cor_YX = s$r, which_predictor = s$j)$value)
  }
  # The probe: .9016 at 187 (the published plan, 182, attains .8616).
  plan <- ss_aipe_reg_coef(width = .3, true_cor_XX = .probe_R,
                           true_cor_YX = .probe_r, which_predictor = 1,
                           assurance = .9)
  expect_equal(.val(plan, "assurance_attained"), .9016055, tolerance = 1e-6)
  expect_equal(.raw_prob_integral(182, 2, .R2s(.probe_R, .probe_r, 1)[["R2"]],
                                  .R2s(.probe_R, .probe_r, 1)[["R2_j"]], 1, .3),
               .8616, tolerance = 1e-4)
  # Unequal tails: the interval has the two t quantiles of its own tails.
  plan <- ss_aipe_reg_coef(
    width = .3, true_cor_XX = .probe_R, true_cor_YX = .probe_r, which_predictor = 1,
    conf_level = NULL, alpha_lower = .01, alpha_upper = .04, assurance = .9)
  rr <- .R2s(.probe_R, .probe_r, 1)
  N <- .val(plan, "necessary_N")
  expect_equal(N, 201)
  expect_gte(.raw_prob_integral(N, 2, rr[["R2"]], rr[["R2_j"]], 1, .3, .01, .04), .9)
  expect_lt(.raw_prob_integral(N - 1, 2, rr[["R2"]], rr[["R2_j"]], 1, .3, .01, .04), .9)
})

test_that("the probability of the central plan matches raw multivariate normal data", {
  skip_on_cran()
  # The probe at its plan, 187: samples of raw multivariate normal data
  # (no DMAR code), least squares, and the central t interval; with tails
  # .01 and .04 the same samples against the F law of that interval.
  plan <- ss_aipe_reg_coef(width = .3, true_cor_XX = .probe_R,
                           true_cor_YX = .probe_r, which_predictor = 1,
                           assurance = .9)
  N <- .val(plan, "necessary_N")
  set.seed(20261111)
  Sigma <- rbind(c(1, .probe_r), cbind(.probe_r, .probe_R))
  Rc <- chol(Sigma)
  G <- 4000L
  w <- rep(NA_real_, G)
  for (g in seq_len(G)) {
    D <- matrix(rnorm(N * 3), N) %*% Rc
    X <- cbind(1, D[, -1])
    XtXi <- tryCatch(solve(crossprod(X)), error = function(e) NULL)
    if (is.null(XtXi)) next
    b <- XtXi %*% crossprod(X, D[, 1])
    w[g] <- sqrt(sum((D[, 1] - X %*% b)^2) / (N - 3) * XtXi[2, 2])
  }
  usable <- sum(is.finite(w))
  expect_identical(usable, G)  # attempted 4000, failed 0
  rr <- .R2s(.probe_R, .probe_r, 1)
  for (tl in list(c(.025, .025), c(.01, .04))) {
    share <- mean((qt(1 - tl[1], N - 3) + qt(1 - tl[2], N - 3)) * w[is.finite(w)] <= .3)
    attained <- .raw_prob_integral(N, 2, rr[["R2"]], rr[["R2_j"]], 1, .3, tl[1], tl[2])
    mcse <- sqrt(attained * (1 - attained) / usable)
    expect_lt(abs(share - attained), 4 * mcse)
  }
})

test_that("with unequal tails the central plans are for the interval ci_reg_coef() reports", {
  skip_on_cran()  # numerical oracle check, run in the full local suite
  # DMAR 1.0.x took both quantiles at the average tail: 157 for the probe
  # with tails .01 and .04, where that interval is .312 wide.
  for (s in list(list(R = .probe_R, r = .probe_r, j = 1, w = .3),
                 list(R = .km_R3, r = .km_r3, j = 2, w = .25)))
    for (tl in list(c(.01, .04), c(.001, .049), c(.1, .05))) for (sd in c(1, 20)) {
      plan <- function(tl, fun = ss_aipe_reg_coef) fun(
        width = s$w * sd, true_cor_XX = s$R, true_cor_YX = s$r, which_predictor = s$j,
        true_sd_Y = sd, conf_level = NULL, alpha_lower = tl[1], alpha_upper = tl[2])
      N <- plan(tl)$value
      rr <- .R2s(s$R, s$r, s$j)
      wd <- function(n) .planner_width(n, rr[["R2"]], rr[["R2_j"]], ncol(s$R), sd_ratio = sd,
                                       tails = tl)
      expect_true(wd(N) <= s$w * sd && wd(N - 1) > s$w * sd)
      expect_identical(plan(rev(tl)), plan(tl))
      expect_identical(plan(tl, ss_aipe_rc)$value, N)
    }
  probe <- function(tl, ...) ss_aipe_reg_coef(
    width = .3, true_cor_XX = .probe_R, true_cor_YX = .probe_r, which_predictor = 1,
    conf_level = NULL, alpha_lower = tl[1], alpha_upper = tl[2], ...)
  expect_equal(c(probe(c(.01, .04))$value, probe(c(.025, .025))$value), c(169, 156))
  expect_no_warning(src <- ss_aipe_src(
    width = .3, true_cor_XX = .probe_R, true_cor_YX = .probe_r, which_predictor = 1,
    noncentral = FALSE, conf_level = NULL, alpha_lower = .04, alpha_upper = .01))
  expect_equal(src$value, 169)
  # The chi square step (Kelley & Maxwell, 2003, Equation 3) at 169 with
  # both quantiles: 194 (DMAR 1.0.x: 182).
  rr <- .R2s(.probe_R, .probe_r, 1)
  K <- (1 - .expected_R2(rr[["R2"]], 169, 2)) / (1 - .expected_R2(rr[["R2_j"]], 169, 1))
  chisq <- function(tl) .chisq_plan(
    width = .3, true_cor_XX = .probe_R, true_cor_YX = .probe_r, which_predictor = 1,
    conf_level = NULL, alpha_lower = tl[1], alpha_upper = tl[2], assurance = .9)
  N <- chisq(c(.01, .04))$value[1]
  expect_equal(c(N, ceiling(((qt(.99, 166) + qt(.96, 166)) / .3)^2 * K *
                              qchisq(.9, 168) / 166 + 3)), c(194, 194))
  # It is the analytic plan of ss_aipe_src(noncentral = FALSE).
  expect_identical(ss_aipe_src(
    width = .3, true_cor_XX = .probe_R, true_cor_YX = .probe_r, which_predictor = 1,
    noncentral = FALSE, conf_level = NULL, alpha_lower = .01, alpha_upper = .04,
    assurance = .9)$value[1], 194)
  # The steps with assurance are the same with the tails exchanged.
  expect_identical(probe(c(.04, .01), assurance = .9),
                   probe(c(.01, .04), assurance = .9))
  expect_identical(chisq(c(.04, .01)), chisq(c(.01, .04)))
  # Noncentral: at a coefficient of 0 the limits are the normal quantiles;
  # 248 (the central plan is 255).
  wd <- function(n) (qnorm(1 - 1e-5) + qnorm(.95001)) *
    sqrt((1 - .expected_R2(.1, n, 2)) / ((1 - .expected_R2(.2, n, 1)) * (n - 3)))
  for (tl in list(c(1e-5, .04999), c(.04999, 1e-5))) {
    a <- list(true_R2_Y_X = .1, true_R2_j_X_without_j = .2, p = 2, width = .4,
              conf_level = NULL, alpha_lower = tl[1], alpha_upper = tl[2])
    N <- do.call(ss_aipe_reg_coef, c(a, true_reg_coef_j = 0, noncentral = TRUE))$value
    expect_true(N == 248 && wd(N) <= .4 && wd(N - 1) > .4)
    expect_identical(do.call(ss_aipe_src, c(a, true_std_reg_coef_j = 0))$value, N)
  }
})

# The smallest N, from p + 2, whose planned width (.planner_width()) meets
# the target, by a scan of every N.
.smallest_N <- function(width, R2, R2_j, p, tails) {
  N <- p + 2
  while (.planner_width(N, R2, R2_j, p, tails = tails) > width) N <- N + 1
  N
}

test_that("with unequal tails the plans take the quantiles of both tails, without a warning", {
  skip_on_cran()
  rr <- .R2s(.probe_R, .probe_r, 1)
  plan <- function(tl, ..., fun = ss_aipe_reg_coef) fun(
    width = .3, true_cor_XX = .probe_R, true_cor_YX = .probe_r, which_predictor = 1,
    conf_level = NULL, alpha_lower = tl[1], alpha_upper = tl[2], ...)
  # The oracle: 156 with both quantiles at the average tail, 169 with the
  # two quantiles of the reported interval.
  expect_identical(.smallest_N(.3, rr[["R2"]], rr[["R2_j"]], 2, rep(.025, 2)), 156)
  expect_identical(.smallest_N(.3, rr[["R2"]], rr[["R2_j"]], 2, c(.01, .04)), 169)
  # Equal tails, and tails whose plans are those at the average tail (.0249
  # and .0251: the oracle gives 156 either way).
  expect_identical(.smallest_N(.3, rr[["R2"]], rr[["R2_j"]], 2, c(.0249, .0251)), 156)
  for (tl in list(c(.025, .025), c(.0249, .0251))) for (a in list(NULL, .9))
    for (m in list(NULL, "analytic")) {
      expect_no_warning(do.call(plan, c(list(tl), assurance = a,
                                        assurance_method = m)))
    }
  # With tails .0249 and .0251 the plan with assurance is 187, where the chi
  # square step at the average tail is 181.
  expect_no_warning(expect_identical(plan(c(.0249, .0251), assurance = .9)$value[1],
                                     187))
  expect_identical(.chisq_plan(width = .3, true_cor_XX = .probe_R, true_cor_YX = .probe_r,
                               which_predictor = 1, conf_level = NULL,
                               alpha_lower = .0249, alpha_upper = .0251,
                               assurance = .9)$value[1], 181)
  # The plan without assurance, in each function's name.
  expect_no_warning(expect_identical(plan(c(.01, .04))$value, 169))
  expect_identical(round(.planner_width(157, rr[["R2"]], rr[["R2_j"]], 2, tails = c(.01, .04)), 3),
                   .311)
  expect_no_warning(expect_identical(plan(c(.04, .01), fun = ss_aipe_rc)$value, 169))
  # The chi square step: 194, and 181 at the average tail (Kelley &
  # Maxwell, 2003, Equation 3, coded here).
  km3 <- function(q, N) ceiling((q / .3)^2 * (1 - .expected_R2(rr[["R2"]], N, 2)) /
    (1 - .expected_R2(rr[["R2_j"]], N, 1)) * qchisq(.9, N - 1) / (N - 3) + 3)
  expect_identical(c(km3(2 * qt(.975, 153), 156), km3(qt(.99, 166) + qt(.96, 166), 169)), c(181, 194))
  # The plan with assurance is the exact plan, 201.
  expect_no_warning(expect_identical(plan(c(.01, .04), assurance = .9)$value[1], 201))
  # The chi square step with the tails given and at the average tail, coded
  # here.
  km3_at <- function(R2, R2_j, p, tl, a) {
    N <- .smallest_N(.4, R2, R2_j, p, tl)
    ceiling(((qt(1 - tl[1], N - p - 1) + qt(1 - tl[2], N - p - 1)) / .4)^2 *
      (1 - .expected_R2(R2, N, p)) / (1 - .expected_R2(R2_j, N, p - 1)) *
      qchisq(a, N - 1) / (N - p - 1) + p + 1)
  }
  wide <- function(R, r, tl, a, ...) ss_aipe_reg_coef(
    width = .4, true_cor_XX = R, true_cor_YX = r, which_predictor = 1, conf_level = NULL,
    alpha_lower = tl[1], alpha_upper = tl[2], assurance = a, ...)
  # The probe with width .4, tails .2 and .1, and an assurance of .6: 31
  # (the probability .6166, and .5777 at 30), where the chi square step is
  # 31 at the average tail and 32 with the tails given.
  expect_identical(c(km3_at(rr[["R2"]], rr[["R2_j"]], 2, c(.15, .15), .6),
                     km3_at(rr[["R2"]], rr[["R2_j"]], 2, c(.2, .1), .6)), c(31, 32))
  expect_equal(vapply(30:31, .raw_prob_integral, 0, p = 2, R2 = rr[["R2"]], R2_j = rr[["R2_j"]],
                      sd_ratio = 1, width = .4, a_lower = .2, a_upper = .1),
               c(.5777, .6166), tolerance = 1e-4)
  expect_no_warning(expect_identical(wide(.probe_R, .probe_r, c(.2, .1), .6)$value[1], 31))
  # One predictor correlated .35 with Y, tails .03 and .07, and an
  # assurance of .7: 73 (.7009, and .6783 at 72), where the chi square step
  # with the tails given is 72 and at the average tail 69.
  r1 <- .R2s(matrix(1), .35, 1)
  expect_identical(c(km3_at(r1[["R2"]], 0, 1, c(.05, .05), .7),
                     km3_at(r1[["R2"]], 0, 1, c(.03, .07), .7)), c(69, 72))
  expect_equal(vapply(72:73, .raw_prob_integral, 0, p = 1, R2 = r1[["R2"]], R2_j = 0,
                      sd_ratio = 1, width = .4, a_lower = .03, a_upper = .07),
               c(.6783, .7009), tolerance = 1e-4)
  expect_no_warning(expect_identical(wide(matrix(1), .35, c(.03, .07), .7)$value[1], 73))
  # assurance_method = "analytic" given is the method of the default.
  expect_no_warning(expect_identical(
    plan(c(.01, .04), assurance = .9, assurance_method = "analytic")$value[1], 201))
  # ss_aipe_src() and the noncentral plans.
  expect_no_warning(plan(c(.01, .04), fun = ss_aipe_src))
  expect_no_warning(plan(c(.01, .04), fun = ss_aipe_src, noncentral = FALSE))
  expect_no_warning(plan(c(.01, .04), fun = ss_aipe_src, noncentral = FALSE,
                         assurance = .9))
  expect_no_warning(plan(c(.01, .04), noncentral = TRUE))
  expect_no_warning(plan(c(.01, .04), noncentral = TRUE, assurance = .9))
})

test_that("with unequal tails ss_aipe_reg_coef() and ss_aipe_rc() plan 169, and 201 with assurance", {
  skip_on_cran()
  for (fun in c("ss_aipe_reg_coef", "ss_aipe_rc")) {
    plan <- function(...) get(fun)(width = .3, true_cor_XX = .probe_R, true_cor_YX = .probe_r,
                                   which_predictor = 1, conf_level = NULL,
                                   alpha_lower = .01, alpha_upper = .04, ...)
    expect_no_warning(expect_identical(plan()$value[1], 169, info = fun))
    expect_no_warning(expect_identical(plan(assurance = .9)$value[1], 201,
                                       info = fun))
  }
})

test_that("reversing the sign of the coefficient gives the same plan", {
  for (a in c(.8, .99)) {
    expect_identical(
      ss_aipe_reg_coef(width = .3, true_cor_XX = .probe_R,
                       true_cor_YX = -.probe_r, which_predictor = 1,
                       assurance = a),
      ss_aipe_reg_coef(width = .3, true_cor_XX = .probe_R,
                       true_cor_YX = .probe_r, which_predictor = 1,
                       assurance = a))
    expect_identical(
      ss_aipe_reg_coef(true_R2_Y_X = .3, true_R2_j_X_without_j = .2,
                       p = 3, true_reg_coef_j = -.4, width = .3,
                       assurance = a),
      ss_aipe_reg_coef(true_R2_Y_X = .3, true_R2_j_X_without_j = .2,
                       p = 3, true_reg_coef_j = .4, width = .3,
                       assurance = a))
  }
})

test_that("reversing the sign of the coefficient gives the same simulated plan", {
  skip_on_cran()  # the simulated step runs in the full local suite
  # The noncentral step simulates; with the same seed the draws are the
  # same, and the plan and its attained probability are identical.
  nc <- function(r, ...) ss_aipe_reg_coef(
    width = .3, true_cor_XX = .km_R3, true_cor_YX = r, which_predictor = 1,
    noncentral = TRUE, assurance = .9, assurance_method = "simulation",
    G = 2000, seed = 3, ...)
  expect_identical(nc(-.km_r3), nc(.km_r3))
  # With unequal tails, reversing the sign exchanges them.
  expect_identical(nc(-.km_r3, conf_level = NULL, alpha_lower = .01, alpha_upper = .04),
                   nc(.km_r3, conf_level = NULL, alpha_lower = .04, alpha_upper = .01))
  src <- function(r, ...) ss_aipe_src(
    width = .3, true_cor_XX = .probe_R, true_cor_YX = r, which_predictor = 1,
    assurance = .85, assurance_method = "simulation", G = 1000, seed = 8, ...)
  expect_identical(src(-.probe_r), src(.probe_r))
  # The central interval: the tails exchanged give the same plan.
  cen <- function(tl) src(.probe_r, noncentral = FALSE, conf_level = NULL,
                          alpha_lower = tl[1], alpha_upper = tl[2])
  expect_identical(cen(c(.04, .01)), cen(c(.01, .04)))
})

test_that("changing the units of Y and the predictor leaves the plan for the coefficient as it was", {
  # The coefficient and the width scale by sd_Y / sd_X = 10 / .5 = 20.
  base <- ss_aipe_reg_coef(width = .3, true_cor_XX = .probe_R,
                           true_cor_YX = .probe_r, which_predictor = 1,
                           assurance = .9)
  for (sd in list(c(10, .5), c(.01, 3), c(1, 7))) {
    ratio <- sd[1] / sd[2]
    units <- ss_aipe_reg_coef(
      width = .3 * ratio, true_cor_XX = .probe_R, true_cor_YX = .probe_r,
      which_predictor = 1, true_sd_Y = sd[1], true_sd_X = sd[2], assurance = .9)
    expect_identical(units$value[1], base$value[1])
    expect_equal(units$value, base$value, tolerance = 1e-12)
    rc <- ss_aipe_rc(
      true_R2_Y_X = .3, true_R2_j_X_without_j = .2, p = 3,
      true_reg_coef_j = .4 * ratio, width = .25 * ratio, true_sd_Y = sd[1],
      true_sd_X = sd[2], assurance = .95)
    rc_base <- ss_aipe_rc(
      true_R2_Y_X = .3, true_R2_j_X_without_j = .2, p = 3,
      true_reg_coef_j = .4, width = .25, assurance = .95)
    expect_identical(rc$value[1], rc_base$value[1])
    expect_equal(rc$value, rc_base$value, tolerance = 1e-12)
  }
})

test_that("the published approximations are those of DMAR 1.0.x from this release's plans", {
  skip_on_cran()
  # DMAR 1.0.0 and MBESS 5.0.1 ss.aipe.reg.coef() agreed on every call
  # below and gave one more case in each (174, 182, 203, 180, 487, 478,
  # 271, 182, 232, 300, 176, 233, 362, 377): they took the expected squared
  # multiple correlation of X_j with the other predictors on p predictors,
  # where the regression has p - 1 ("Earlier versions" in
  # ?ss_aipe_reg_coef). The calls: the probe, the five exchangeable
  # predictors of Kelley and Maxwell (2003, pp. 312--313), their
  # three-predictor structure (p. 313), and one predictor.
  # The central approximation is no longer a plan of ss_aipe_reg_coef()
  # (its analytic plan is exact) and is the engine's step that
  # ss_aipe_src(noncentral = FALSE) uses; the noncentral one is the analytic
  # plan, the default.
  pub <- function(R, r, j, width, a, nc = FALSE) suppressWarnings(
    if (nc) ss_aipe_reg_coef(width = width, true_cor_XX = R, true_cor_YX = r,
                             which_predictor = j, noncentral = TRUE, assurance = a)
    else .chisq_plan(width = width, true_cor_XX = R, true_cor_YX = r,
                     which_predictor = j, assurance = a))
  expect_equal(pub(.probe_R, .probe_r, 1, .3, .8)$value, c(173, .8, NA))
  expect_equal(pub(.probe_R, .probe_r, 1, .3, .9)$value, c(181, .9, NA))
  expect_equal(pub(.probe_R, .probe_r, 1, .3, .99)$value, c(202, .99, NA))
  expect_equal(pub(.probe_R, .probe_r, 1, .3, .9, nc = TRUE)$value, c(179, .9, NA))
  nc_pub <- pub(.probe_R, .probe_r, 1, .3, .9, nc = TRUE)
  expect_identical(attr(nc_pub, "assurance_method"), "analytic")
  expect_identical(attr(nc_pub, "assurance_probability"), "not computed")
  expect_identical(ss_aipe_reg_coef(width = .3, true_cor_XX = .probe_R,
                                    true_cor_YX = .probe_r, which_predictor = 1,
                                    noncentral = TRUE, assurance = .9,
                                    assurance_method = "analytic"), nc_pub)
  expect_match(paste(utils::capture.output(print(nc_pub)), collapse = " "),
               "the plan is the published approximation", fixed = TRUE)
  expect_equal(pub(.exch(5, .4), rep(.3, 5), 1, .2, .8)$value[1], 486)
  expect_equal(pub(.exch(5, .4), rep(.3, 5), 1, .2, .8, nc = TRUE)$value[1], 477)
  expect_equal(vapply(1:3, function(j) pub(.km_R3, .km_r3, j, .3, .9)$value[1], 0),
               c(270, 181, 231))
  expect_equal(vapply(1:3, function(j) pub(.km_R3, .km_r3, j, .3, .9, nc = TRUE)$value[1], 0),
               c(299, 175, 232))
  expect_equal(pub(matrix(1), .4, 1, .2, .9)$value[1], 361)
  expect_equal(pub(matrix(1), .4, 1, .2, .9, nc = TRUE)$value[1], 376)
  # The central approximation is Kelley and Maxwell's (2003) Equation 3,
  # coded here, at the plan without assurance found by a scan of every N.
  eq3 <- function(R, r, j, width, a) {
    rr <- .R2s(R, r, j)
    p <- ncol(R)
    N <- .smallest_N(width, rr[["R2"]], rr[["R2_j"]], p, c(.025, .025))
    ceiling((2 * qt(.975, N - p - 1) / width)^2 * (1 - .expected_R2(rr[["R2"]], N, p)) /
              (1 - .expected_R2(rr[["R2_j"]], N, p - 1)) * qchisq(a, N - 1) /
              (N - p - 1) + p + 1)
  }
  expect_equal(eq3(.probe_R, .probe_r, 1, .3, .9), 181)
  expect_equal(eq3(.exch(5, .4), rep(.3, 5), 1, .2, .8), 486)
  expect_equal(vapply(1:3, function(j) eq3(.km_R3, .km_r3, j, .3, .9), 0), c(270, 181, 231))
  expect_equal(eq3(matrix(1), .4, 1, .2, .9), 361)
  # The wrapper of the standardized coefficient, for both intervals.
  expect_equal(ss_aipe_src(width = .3, true_cor_XX = .probe_R, true_cor_YX = .probe_r,
                           which_predictor = 1, assurance = .9,
                           noncentral = FALSE)$value, c(181, .9, NA))
  expect_equal(ss_aipe_src(width = .3, true_cor_XX = .probe_R, true_cor_YX = .probe_r,
                           which_predictor = 1, assurance = .9)$value, c(179, .9, NA))
  # The examples of Kelley and Maxwell (2008, p. 186), as printed there:
  # ss.aipe.reg.coef(Rho2.Y_X = 0.5, Rho2.j_X.without.j = 0.2, p = 5,
  # b.j = 3, width = 2, sigma.Y = 50, sigma.X = 5, conf.level = 0.95)
  # gives 250, with degree.of.certainty = 0.85 it gives 278, and the
  # standardized analog (b.j = 0.3, width = 0.2, unit standard deviations,
  # Noncentral = TRUE) gives 264; with the expected squared multiple
  # correlation of X_j on p - 1 predictors each is one smaller. (The text
  # gives 50 and 5 as variances; the calls pass them as standard
  # deviations, and so does this test.)
  km08 <- function(...) ss_aipe_reg_coef(
    true_R2_Y_X = .5, true_R2_j_X_without_j = .2, p = 5, ...)
  expect_equal(km08(true_reg_coef_j = 3, width = 2, true_sd_Y = 50,
                    true_sd_X = 5)$value, 249)
  expect_equal(.chisq_plan(true_R2_Y_X = .5, true_R2_j_X_without_j = .2, p = 5,
                           true_reg_coef_j = 3, width = 2, true_sd_Y = 50,
                           true_sd_X = 5, assurance = .85)$value[1], 277)
  expect_equal(km08(true_reg_coef_j = .3, width = .2, noncentral = TRUE)$value, 263)
  # The analytic plan for the first: 280, where the probability is .8531
  # (at 279 it is below .85).
  exact <- km08(true_reg_coef_j = 3, width = 2, true_sd_Y = 50, true_sd_X = 5,
                assurance = .85)
  expect_equal(exact$value[1], 280)
  expect_lt(.raw_prob_integral(278, 5, .5, .2, 10, 2), .85)
  expect_lt(.raw_prob_integral(279, 5, .5, .2, 10, 2), .85)
  expect_equal(exact$value[3], .raw_prob_integral(280, 5, .5, .2, 10, 2),
               tolerance = 1e-8)
  # The published approximation draws nothing.
  set.seed(5); before <- .Random.seed
  invisible(pub(.probe_R, .probe_r, 1, .3, .9, nc = TRUE))
  expect_identical(.Random.seed, before)
})

test_that("the plans without assurance are the smallest N from p + 2 that a scan of every N finds", {
  # A scan of every N from p + 2, the smallest sample with a residual
  # degree of freedom, coded here: .smallest_N() for the central interval,
  # and for the noncentral one the limits of the noncentrality at the
  # t statistic of the expected squared multiple correlations by uniroot()
  # on pt(), times the standard error. MBESS 5.0.1 ss.aipe.reg.coef()
  # (width = .3, the structure of Kelley and Maxwell, 2003, p. 313,
  # degree.of.certainty = .8, the chi square step here) gave
  # c(260, 174, 222), and c(294, 171, 227) with Noncentral = TRUE, one more
  # case in each than DMAR: it took the expected squared multiple
  # correlation of X_j on p predictors and began its search at the larger
  # of 2p + 1 and ten below the normal approximation ("Earlier versions"
  # in ?ss_aipe_reg_coef).
  plan <- function(...) ss_aipe_reg_coef(...)$value[1]
  expect_identical(vapply(1:3, function(j) .chisq_plan(
    width = .3, true_cor_XX = .km_R3, true_cor_YX = .km_r3, which_predictor = j,
    assurance = .8)$value[1], 0), c(259, 173, 221))
  for (j in 1:3) {
    rr <- .R2s(.km_R3, .km_r3, j)
    expect_identical(plan(width = .3, true_cor_XX = .km_R3, true_cor_YX = .km_r3,
                          which_predictor = j),
                     .smallest_N(.3, rr[["R2"]], rr[["R2_j"]], 3, c(.025, .025)))
  }
  # MBESS and DMAR 1.0.x returned 21 here, 2p + 1, the start of their
  # search, where 14 meets the width.
  args <- list(true_R2_Y_X = .8, true_R2_j_X_without_j = 0, p = 10,
               true_reg_coef_j = .3, width = 2)
  expect_identical(do.call(plan, args), 14)
  expect_identical(.smallest_N(2, .8, 0, 10, c(.025, .025)), 14)
  # The floor: a plan of p + 2 when that sample already meets the width.
  expect_identical(plan(true_R2_Y_X = .5, true_R2_j_X_without_j = 0, p = 3,
                        true_reg_coef_j = .3, width = 50), 5)
  expect_identical(plan(true_R2_Y_X = .5, true_R2_j_X_without_j = 0, p = 3,
                        true_reg_coef_j = .3, width = 50, noncentral = TRUE), 5)
  skip_on_cran()
  nc_width <- function(N, R2, R2_j, p, beta, tails) {
    se <- sqrt((1 - .expected_R2(R2, N, p)) /
                 ((1 - .expected_R2(R2_j, N, p - 1)) * (N - p - 1)))
    t <- beta / se
    # pt() warns that full precision may not have been achieved at the ends
    # of the bracket, ten from t, where it is near 0 or 1; at each root
    # found it gives no warning and the target to 1e-12 of itself.
    f <- function(ncp, target) suppressWarnings(pt(t, N - p - 1, ncp)) - target
    lim <- vapply(c(1 - tails[1], tails[2]), function(target)
      uniroot(f, c(t - 10, t + 10), target = target, tol = 1e-12,
              extendInt = "yes")$root, 0)
    (lim[2] - lim[1]) * se
  }
  nc_scan <- function(width, R2, R2_j, p, beta, tails) {
    N <- p + 2
    while (nc_width(N, R2, R2_j, p, beta, tails) > width) N <- N + 1
    N
  }
  for (j in 1:3) {
    rr <- .R2s(.km_R3, .km_r3, j)
    beta <- solve(.km_R3, .km_r3)[j]
    expect_identical(plan(width = .3, true_cor_XX = .km_R3, true_cor_YX = .km_r3,
                          which_predictor = j, noncentral = TRUE),
                     nc_scan(.3, rr[["R2"]], rr[["R2_j"]], 3, beta, c(.025, .025)))
  }
  # DMAR 1.0.x and the first builds of 1.1.0 returned 358 here, from a
  # search that began six below the central plan.
  args <- list(true_R2_Y_X = .1, true_R2_j_X_without_j = .2, p = 2,
               true_reg_coef_j = 0, width = .4, noncentral = TRUE,
               conf_level = NULL, alpha_lower = .0002, alpha_upper = .0002)
  expect_identical(do.call(plan, args), nc_scan(.4, .1, .2, 2, 0, c(.0002, .0002)))
  expect_identical(do.call(plan, args), 355)
})

test_that("the search of the planners finds the first N that meets, from any guess", {
  skip_on_cran()
  # .reg_coef_first_n() steps from the guess by doubling and then bisects;
  # the cost grows with the logarithm of the distance from the guess.
  for (answer in c(4, 5, 37, 1000, 123457)) {
    for (guess in c(1, 4, answer - 1, answer, answer + 1, 3 * answer, 1e7)) {
      count <- new.env()
      count$n <- 0
      meets <- function(n) {
        count$n <- count$n + 1
        n >= answer
      }
      expect_identical(DMAR:::.reg_coef_first_n(meets, lower = 4, guess = guess), answer)
      expect_lte(count$n, 2 * log2(abs(guess - answer) + 2) + 3)
    }
  }
  # Nothing below `lower` is tried, and `lower` is returned when it meets.
  expect_identical(DMAR:::.reg_coef_first_n(function(n) {
    stopifnot(n >= 10)
    TRUE
  }, lower = 10, guess = 500), 10)
})

test_that("with one predictor the squared multiple correlation of X_j is 0", {
  # The regression of X_j on the other p - 1 predictors has none when
  # p = 1, so its expected squared multiple correlation is 0 at every N,
  # and a planning value other than 0 is refused.
  expect_identical(.expected_R2(0, 50, 0), 0)
  e <- tryCatch(ss_aipe_reg_coef(true_R2_Y_X = .16, true_R2_j_X_without_j = .2, p = 1,
                                 true_reg_coef_j = .4, width = .2),
                error = function(e) e)
  expect_s3_class(e, "dmar_argument_invalid")
  expect_identical(e$arg, "true_R2_j_X_without_j")
  N <- ss_aipe_reg_coef(true_R2_Y_X = .16, true_R2_j_X_without_j = 0, p = 1,
                        true_reg_coef_j = .4, width = .2)$value
  width_at <- function(n) 2 * qt(.975, n - 2) * sqrt((1 - .expected_R2(.16, n, 1)) / (n - 2))
  expect_true(width_at(N) <= .2 && width_at(N - 1) > .2)
})

# The noncentral step ----------------------------------------------------

test_that("the draws of the noncentral step are the statistics least squares gives", {
  skip_on_cran()
  # The step draws entries of the triangular factor L of the Wishart matrix
  # of sums of squares and products of (Z_1, ..., Z_{p-1}, U, E), then
  # forms the statistics by closed expressions. Here the same entries are
  # put in L, the matrix A = L L' is formed, mapped to (X_{-j}, X_j, Y) by
  # the population's linear map, and the statistics are found by solve().
  p <- 4; N <- 40; beta <- .35; R2 <- .4; R2_j <- .3; ad <- .05
  set.seed(9)
  d <- DMAR:::.reg_coef_std_draws(6, N, p, beta, R2, R2_j)
  dd <- DMAR:::.reg_coef_std_d(R2, R2_j, beta, ad, p)
  s_u <- sqrt(1 - R2_j); s_e <- sqrt(1 - R2)
  Tm <- rbind(cbind(diag(3), 0, 0),
              c(sqrt(R2_j), 0, 0, s_u, 0),
              c(dd, 0, beta * s_u, s_e))
  # The population: unit variances of X_j and Y, and the two squared
  # multiple correlations and the coefficient given.
  P <- Tm %*% t(Tm)
  expect_equal(diag(P)[4:5], c(1, 1))
  expect_equal(drop(P[5, 1:4] %*% solve(P[1:4, 1:4], P[1:4, 5])), R2)
  expect_equal(drop(P[4, 1:3] %*% solve(P[1:3, 1:3], P[1:3, 4])), R2_j)
  expect_equal(solve(P[1:4, 1:4], P[1:4, 5])[4], beta)
  expect_equal(P[4, 5] - beta * (1 - R2_j), ad)
  one <- function(t) rep(1, length(t))
  w1 <- DMAR:::.reg_coef_std_width(d, ad, one)
  for (g in 1:6) {
    L <- matrix(0, p + 1, p + 1)
    L[1, 1] <- d$l_11[g]
    L[2, 1:2] <- c(d$l_21[g], d$l_22[g])
    L[3, 1:3] <- c(.7, -1.2, 5.3)  # the third Z does not enter the statistics
    L[4, ] <- c(d$l_u1[g], d$l_u2[g], d$l_um[g, 1], d$l_u[g], 0)
    L[5, ] <- c(d$l_e1[g], d$l_e2[g], d$l_em[g, 1], d$l_eu[g], d$l_e[g])
    S <- Tm %*% (L %*% t(L)) %*% t(Tm)
    b <- solve(S[1:4, 1:4], S[1:4, 5])
    rss <- S[5, 5] - sum(S[1:4, 5] * b)
    ssr_j <- 1 / solve(S[1:4, 1:4])[4, 4]
    se_raw <- sqrt(rss / ((N - p - 1) * ssr_j))
    expect_equal(d$t[g], b[4] / se_raw, tolerance = 1e-10)
    expect_equal(d$se_raw[g], se_raw, tolerance = 1e-10)
    expect_equal(d$ss_x[g], S[4, 4], tolerance = 1e-10)
    # The standard error of the standardized coefficient,
    # sqrt{(1 - R2_hat) / [(1 - R2_j_hat)(N - p - 1)]}.
    R2_hat <- 1 - rss / S[5, 5]
    R2_j_hat <- 1 - ssr_j / S[4, 4]
    expect_equal(w1[g], sqrt((1 - R2_hat) / ((1 - R2_j_hat) * (N - p - 1))),
                 tolerance = 1e-10)
  }
})

test_that("the interpolated noncentral distance is that of ci_nc_t()", {
  skip_on_cran()
  set.seed(4)
  d <- DMAR:::.reg_coef_std_draws(5000, 60, 3, .5, .4, .2)
  for (tails in list(c(.025, .025), c(.01, .04))) {
    delta <- DMAR:::.reg_coef_nc_delta(d$t, 56, tails[1], tails[2])
    tt <- c(range(d$t), sample(d$t, 25))
    direct <- vapply(tt, function(v) diff(ci_nc_t(
      ncp = v, df = 56, conf_level = NULL, alpha_lower = tails[1],
      alpha_upper = tails[2])$value), 0)
    expect_lt(max(abs(delta(tt) / direct - 1)), 1e-8)
  }
})

test_that("the noncentral plan attains the request in raw multivariate normal data", {
  skip_on_cran()
  plan <- ss_aipe_reg_coef(
    width = .3, true_cor_XX = .probe_R, true_cor_YX = .probe_r, which_predictor = 1,
    noncentral = TRUE, assurance = .9, assurance_method = "simulation", seed = 11)
  expect_identical(plan$term, c("necessary_N", "assurance", "assurance_attained", "G"))
  expect_identical(attr(plan, "assurance_method"), "simulation")
  expect_identical(attr(plan, "assurance_probability"), "simulation")
  N <- .val(plan, "necessary_N")
  expect_equal(N, 176)
  attained <- .val(plan, "assurance_attained")
  mc_se <- attr(plan, "assurance_mc_se")
  expect_equal(mc_se, sqrt(attained * (1 - attained) / 20000))
  # Raw data (no DMAR code): standardize, least squares, and the
  # noncentral t limits found here by uniroot() on pt().
  nc_width <- function(t, df) {
    hw <- max(10, 8 * sqrt(1 + t^2 / (2 * df)))
    lo <- uniroot(function(nc) pt(t, df, nc, lower.tail = FALSE) - .025,
                  c(t - hw, t), tol = 1e-10, extendInt = "upX")$root
    up <- uniroot(function(nc) pt(t, df, nc) - .025, c(t, t + hw),
                  tol = 1e-10, extendInt = "downX")$root
    up - lo
  }
  set.seed(7)
  Rc <- chol(rbind(c(1, .probe_r), cbind(.probe_r, .probe_R)))
  G <- 2000L
  w <- rep(NA_real_, G)
  for (g in seq_len(G)) {
    D <- scale(matrix(rnorm(N * 3), N) %*% Rc)
    X <- cbind(1, D[, -1])
    XtXi <- solve(crossprod(X))
    b <- XtXi %*% crossprod(X, D[, 1])
    se <- sqrt(sum((D[, 1] - X %*% b)^2) / (N - 3) * XtXi[2, 2])
    w[g] <- tryCatch(nc_width(b[2] / se, N - 3) * se, error = function(e) NA_real_)
  }
  expect_identical(sum(is.finite(w)), G)  # attempted 2000, failed 0
  share <- mean(w <= .3)
  expect_lt(abs(share - attained),
            4 * sqrt(attained * (1 - attained) / G + mc_se^2))
})

test_that("the simulated step uses one set of samples at every N and returns the smallest N their share allows", {
  skip_on_cran()  # the simulated step runs in the full local suite
  # The step draws G samples once (.reg_coef_std_base_draws(), right after
  # the seed is set) and maps them to every N it tries. Rebuilt here from
  # the same seed, the share of those samples is: assurance_attained at the
  # returned N, at least `assurance` there, below it at N - 1 (unless N is
  # the plan without assurance), and never smaller at a larger N.
  sets <- list(
    list(R = .probe_R, r = .probe_r, j = 1, width = .3, a = .9),
    list(R = matrix(c(1, .5, .5, 1), 2), r = c(.7, .2), j = 1, width = .2, a = .9),
    list(R = .km_R3, r = .km_r3, j = 2, width = .3, a = .8),
    list(R = matrix(1), r = .4, j = 1, width = .2, a = .6))
  for (s in sets) {
    p <- ncol(s$R)
    rr <- .R2s(s$R, s$r, s$j)
    beta <- solve(s$R, s$r)[s$j]
    ad <- DMAR:::.reg_coef_std_ad_values(rr[["R2"]], rr[["R2_j"]], beta, p,
                                         s$r[s$j])
    N_0 <- ss_aipe_reg_coef(width = s$width, true_cor_XX = s$R, true_cor_YX = s$r,
                            which_predictor = s$j, noncentral = TRUE)$value
    for (seed in 1:3) {
      plan <- ss_aipe_reg_coef(
        width = s$width, true_cor_XX = s$R, true_cor_YX = s$r, which_predictor = s$j,
        noncentral = TRUE, assurance = s$a, assurance_method = "simulation",
        G = 500, seed = seed)
      N <- .val(plan, "necessary_N")
      attained <- .val(plan, "assurance_attained")
      set.seed(seed)
      base <- DMAR:::.reg_coef_std_base_draws(500, p)
      share <- function(n) min(DMAR:::.reg_coef_std_share(
        DMAR:::.reg_coef_std_draws_at(base, n, beta, rr[["R2"]], rr[["R2_j"]]),
        ad, s$width, .025, .025))
      expect_identical(share(N), attained)
      expect_gte(attained, s$a)
      if (N > N_0) expect_lt(share(N - 1), s$a)
      expect_gte(N, N_0)
      shares <- vapply(seq(N_0, N + 15), share, 0)
      expect_true(all(diff(shares) >= 0))
    }
  }
  # The share at the returned N attains the request.
  plan <- ss_aipe_reg_coef(
    width = .2, true_cor_XX = matrix(1), true_cor_YX = .4, which_predictor = 1,
    noncentral = TRUE, assurance = .6, assurance_method = "simulation",
    G = 500, seed = 1)
  expect_gte(.val(plan, "assurance_attained"), .6)
  # The samples at one N are those of the sampling distribution at that N:
  # the chi entries are chi square quantiles of uniform draws, so their
  # squares have the chi square mean and variance.
  set.seed(12)
  base <- DMAR:::.reg_coef_std_base_draws(20000, 3)
  d <- DMAR:::.reg_coef_std_draws_at(base, 50, .3, .4, .2)
  expect_lt(abs(mean(d$l_u^2) - 47) / sqrt(2 * 47 / 20000), 4)
  expect_lt(abs(mean(d$l_e^2) - 46) / sqrt(2 * 46 / 20000), 4)
  expect_lt(abs(mean(d$l_22^2) - 48) / sqrt(2 * 48 / 20000), 4)
  # The variance of a chi square on 49 degrees of freedom is 98; the
  # standard error of a sample variance of G of them is
  # 98 sqrt{(2 + 12 / 49) / G}.
  expect_lt(abs(var(d$l_11^2) - 98) / (98 * sqrt((2 + 12 / 49) / 20000)), 4)
})

test_that("the squared multiple correlations path takes the worst correlation of Y with the predictor", {
  rr <- .R2s(.km_R3, .km_r3, 1)
  beta <- solve(.km_R3, .km_r3)[1]
  ad_actual <- .km_r3[1] - beta * (1 - rr[["R2_j"]])
  grid <- DMAR:::.reg_coef_std_ad_values(rr[["R2"]], rr[["R2_j"]], beta, 3)
  expect_length(grid, 21L)
  expect_true(min(grid) <= ad_actual && ad_actual <= max(grid))
  # Every value of the grid is a population: |a'd| <= |a| |d|.
  bound <- sqrt(rr[["R2_j"]] * (rr[["R2"]] - beta^2 * (1 - rr[["R2_j"]])))
  expect_equal(max(abs(grid)), bound)
  # With two predictors a'd is one of the two ends.
  rr2 <- .R2s(.probe_R, .probe_r, 1)
  beta2 <- solve(.probe_R, .probe_r)[1]
  ends <- DMAR:::.reg_coef_std_ad_values(rr2[["R2"]], rr2[["R2_j"]], beta2, 2)
  expect_equal(sort(abs(ends)), rep(abs(.probe_r[1] - beta2 * (1 - rr2[["R2_j"]])), 2))
})

test_that("the plan on the squared multiple correlations path guards the worst correlation", {
  skip_on_cran()  # the simulated step runs in the full local suite
  rr <- .R2s(.km_R3, .km_r3, 1)
  beta <- solve(.km_R3, .km_r3)[1]
  # The plan on this path is at least the plan with the correlation known,
  # up to the Monte Carlo error of the two searches.
  r2_path <- ss_aipe_reg_coef(
    true_R2_Y_X = rr[["R2"]], true_R2_j_X_without_j = rr[["R2_j"]], p = 3,
    true_reg_coef_j = beta, width = .3, noncentral = TRUE, assurance = .9,
    assurance_method = "simulation", G = 5000, seed = 2)
  cov_path <- ss_aipe_reg_coef(
    width = .3, true_cor_XX = .km_R3, true_cor_YX = .km_r3, which_predictor = 1,
    noncentral = TRUE, assurance = .9, assurance_method = "simulation",
    G = 5000, seed = 2)
  expect_gte(r2_path$value[1], cov_path$value[1] - 2)
})

test_that("the assurance arguments are checked", {
  skip_on_cran()
  call <- function(...) ss_aipe_reg_coef(width = .3, true_cor_XX = .probe_R,
                                         true_cor_YX = .probe_r, which_predictor = 1, ...)
  expect_identical(eval(formals(ss_aipe_reg_coef)$assurance_method),
                   c("analytic", "simulation"))
  for (bad in list(NA, "yes", c(TRUE, FALSE), 1, "sim")) {
    e <- tryCatch(call(assurance = .9, assurance_method = bad), error = function(e) e)
    expect_s3_class(e, "dmar_argument_invalid")
    expect_identical(e$arg, "assurance_method")
    # The central interval offers "analytic" only, and the message says so.
    expect_match(conditionMessage(e), paste0(
      "the only method ss_aipe_reg_coef() with noncentral = FALSE offers"),
      fixed = TRUE)
  }
  # The formal's default, both methods, is "analytic" for either interval.
  expect_identical(call(assurance = .9,
                        assurance_method = c("analytic", "simulation"))$value,
                   call(assurance = .9)$value)
  # "simulation" is refused for the central interval, whose probability is
  # exact, with or without assurance, and offered for the noncentral one.
  for (a in list(.9, NULL)) {
    e <- tryCatch(call(assurance = a, assurance_method = "simulation"),
                  error = function(e) e)
    expect_s3_class(e, "dmar_argument_invalid")
    expect_identical(e$arg, "assurance_method")
    expect_match(conditionMessage(e), paste0(
      "assurance_method = \"simulation\" is not offered by ss_aipe_reg_coef() ",
      "with noncentral = FALSE: the analytic probability that the interval is ",
      "no wider than 'width' is exact for this interval"), fixed = TRUE)
  }
  e <- tryCatch(ss_aipe_rc(width = .3, true_cor_XX = .probe_R, true_cor_YX = .probe_r,
                           which_predictor = 1, assurance = .9,
                           assurance_method = "simulation"), error = function(e) e)
  expect_s3_class(e, "dmar_argument_invalid")
  expect_match(conditionMessage(e), "is not offered by ss_aipe_rc():", fixed = TRUE)
  expect_identical(formals(ss_aipe_rc)$assurance_method, "analytic")
  # ss_aipe_rc() offers "analytic" only, so both methods together are not
  # its default.
  e <- tryCatch(ss_aipe_rc(width = .3, true_cor_XX = .probe_R, true_cor_YX = .probe_r,
                           which_predictor = 1, assurance = .9,
                           assurance_method = c("analytic", "simulation")),
                error = function(e) e)
  expect_s3_class(e, "dmar_argument_invalid")
  expect_match(conditionMessage(e), "the only method ss_aipe_rc() offers",
               fixed = TRUE)
  expect_identical(call(noncentral = TRUE, assurance_method = "simulation")$value,
                   call(noncentral = TRUE)$value)
  # A limit at infinity is never narrow enough, ci_reg_coef() takes no
  # tail above .5, and two of .5 are a confidence level of 0.
  for (tl in list(c(0, .05), c(.6, .01), c(.5, .5))) {
    e <- tryCatch(call(conf_level = NULL, alpha_lower = tl[1], alpha_upper = tl[2]),
                  error = function(e) e)
    expect_s3_class(e, "dmar_argument_invalid")
    expect_identical(e$arg, c("alpha_lower", "alpha_upper"))
  }
  # G is checked only when the noncentral step simulates.
  expect_error(call(noncentral = TRUE, assurance = .9, G = 99,
                    assurance_method = "simulation"),
               class = "dmar_count_invalid")
  expect_error(ss_aipe_src(width = .3, true_cor_XX = .probe_R, true_cor_YX = .probe_r,
                           which_predictor = 1, assurance = .9, G = 50,
                           assurance_method = "simulation"),
               class = "dmar_count_invalid")
  expect_identical(suppressWarnings(call(assurance = .9, G = 1)),
                   suppressWarnings(call(assurance = .9)))
  expect_identical(call(noncentral = TRUE, G = 1), call(noncentral = TRUE))
  expect_identical(call(noncentral = TRUE, assurance = .9, G = 1),
                   call(noncentral = TRUE, assurance = .9))
  expect_identical(ss_aipe_src(width = .3, true_cor_XX = .probe_R,
                               true_cor_YX = .probe_r, which_predictor = 1,
                               assurance = .9, G = 50)$value[1], 179)
  # Squared multiple correlations no population allows with the coefficient.
  e <- tryCatch(ss_aipe_reg_coef(true_R2_Y_X = .1, true_R2_j_X_without_j = .2, p = 3,
                                 true_reg_coef_j = .5, width = .3, noncentral = TRUE,
                                 assurance = .9, assurance_method = "simulation",
                                 G = 200),
                error = function(e) e)
  expect_s3_class(e, "dmar_argument_invalid")
  expect_true("true_reg_coef_j" %in% e$arg)
})

test_that("the sensitivity functions pass assurance_method to the planner", {
  skip_on_cran()
  sens <- function(f, ...) {
    out <- suppressWarnings(f(true_var_Y = 1, true_cov_YX = .probe_r,
                              true_cov_XX = .probe_R, width = .3, assurance = .9,
                              G = 2, print_iter = FALSE, seed = 1, ...))
    out$value[out$term == "total_N"]
  }
  planned <- function(...) ss_aipe_reg_coef(
    width = .3, true_cor_XX = .probe_R, true_cor_YX = .probe_r,
    which_predictor = 1, assurance = .9, ...)$value[1]
  for (f in list(ss_aipe_reg_coef_sensitivity, ss_aipe_rc_sensitivity)) {
    expect_identical(eval(formals(f)$assurance_method), c("analytic", "simulation"))
    expect_equal(sens(f), planned())
    expect_equal(sens(f), 187)
    expect_equal(sens(f, assurance_method = "analytic"), 187)
    # The noncentral interval: the approximation, 179, by default, and the
    # simulated plan with the planner's G from the sensitivity function's
    # seed.
    expect_equal(sens(f, noncentral = TRUE), 179)
    expect_equal(sens(f, noncentral = TRUE, assurance_method = "simulation"),
                 planned(noncentral = TRUE, assurance_method = "simulation", seed = 1))
    # "simulation" is refused for the central interval, in the function's
    # own name.
    e <- tryCatch(sens(f, assurance_method = "simulation"), error = function(e) e)
    expect_s3_class(e, "dmar_argument_invalid")
    expect_identical(e$arg, "assurance_method")
    expect_match(conditionMessage(e), "with noncentral = FALSE", fixed = TRUE)
    # It names the function's own noncentral interval among those that
    # offer it.
    fn <- if (identical(f, ss_aipe_rc_sensitivity)) "ss_aipe_rc_sensitivity" else
      "ss_aipe_reg_coef_sensitivity"
    expect_match(conditionMessage(e), paste0(
      "Simulation is offered by ", fn, "() with noncentral = TRUE, by ",
      "ss_aipe_sc_ancova(), by ss_aipe_reg_coef() with noncentral = TRUE, ",
      "and by ss_aipe_src()"), fixed = TRUE)
  }
  # The standardized sensitivity function plans the noncentral interval:
  # the approximation, 179, by default, and the simulated plan, near 176.
  expect_identical(eval(formals(ss_aipe_src_sensitivity)$assurance_method),
                   c("analytic", "simulation"))
  expect_equal(sens(ss_aipe_src_sensitivity), 179)
  expect_lte(abs(sens(ss_aipe_src_sensitivity, assurance_method = "simulation") -
                   planned(noncentral = TRUE, assurance_method = "simulation",
                           seed = 1)), 3)
  # With noncentral = FALSE it plans the central interval for the
  # standardized coefficient through ss_aipe_src(): by default the chi
  # square approximation, 181, and with "simulation" the plan from the
  # draws the sensitivity function's seed fixes, 171 (the F law of the
  # interval in the units of the data gives 187).
  expect_equal(sens(ss_aipe_src_sensitivity, noncentral = FALSE), 181)
  src_central <- ss_aipe_src(width = .3, true_cor_XX = .probe_R, true_cor_YX = .probe_r,
                             which_predictor = 1, noncentral = FALSE,
                             assurance = .9, assurance_method = "simulation",
                             seed = 1)$value[1]
  expect_equal(sens(ss_aipe_src_sensitivity, noncentral = FALSE,
                    assurance_method = "simulation"), src_central)
  expect_equal(src_central, 171)
})

test_that("an abbreviated which_width is the value it stands for", {
  skip_on_cran()
  # DMAR 1.0.x checked which_width with char.expand() and planned with the
  # value as given: "L", "Up", and "" gave the plan for the full width
  # (607 here and 626 for ss_aipe_src() then, 606 and 625 now), while
  # "Lower" and "Upper" were refused. An abbreviation of "Lower" or "Upper" is now refused as they
  # are, and "F" and "Fu" are "Full".
  args <- list(true_R2_Y_X = .3, true_R2_j_X_without_j = .2, p = 5,
               true_reg_coef_j = .25, width = .15)
  full <- do.call(ss_aipe_reg_coef, c(args, which_width = "Full"))$value
  expect_identical(full, 606)
  for (w in c("F", "Fu", "Ful")) {
    expect_identical(do.call(ss_aipe_reg_coef, c(args, which_width = w))$value,
                     full)
    expect_identical(do.call(ss_aipe_rc, c(args, which_width = w))$value, full)
  }
  for (w in c("L", "Lo", "Lower", "U", "Up", "Upper")) {
    for (f in list(ss_aipe_reg_coef, ss_aipe_rc)) {
      e <- tryCatch(do.call(f, c(args, which_width = w)), error = function(e) e)
      expect_s3_class(e, "dmar_argument_invalid")
      expect_identical(e$arg, "which_width")
      expect_match(conditionMessage(e), "only the plan for the full width",
                   fixed = TRUE)
    }
    e <- tryCatch(ss_aipe_src(true_R2_Y_X = .3, true_R2_j_X_without_j = .2,
                              p = 5, true_std_reg_coef_j = .25, width = .15,
                              which_width = w),
                  error = function(e) e)
    expect_s3_class(e, "dmar_argument_invalid")
    expect_identical(e$arg, "which_width")
  }
  e <- tryCatch(do.call(ss_aipe_reg_coef, c(args, which_width = "Up")),
                error = function(e) e)
  expect_match(conditionMessage(e),
               "\"Up\" was given, which stands for \"Upper\".", fixed = TRUE)
  expect_identical(ss_aipe_src(true_R2_Y_X = .3, true_R2_j_X_without_j = .2,
                               p = 5, true_std_reg_coef_j = .25, width = .15,
                               which_width = "F")$value, 625)
  # Not an abbreviation of one of the three: refused with the same class,
  # where DMAR 1.0.x planned "" for the full width and stopped on a number
  # or two values with an error of base R.
  for (w in list("", NA_character_, "full", "x", 1, c("Full", "Lower"))) {
    e <- tryCatch(do.call(ss_aipe_reg_coef, c(args, list(which_width = w))),
                  error = function(e) e)
    expect_s3_class(e, "dmar_argument_invalid")
    expect_identical(e$arg, "which_width")
  }
})

test_that("the regression planners and their sensitivity functions plan with assurance without a warning", {
  skip_on_cran()
  R <- matrix(c(1, .3, .3, 1), 2)
  r <- c(.4, .3)
  expect_no_warning(ss_aipe_reg_coef(width = .3, true_cor_XX = R, true_cor_YX = r,
                                     which_predictor = 1, assurance = .9))
  expect_no_warning(ss_aipe_rc(width = .3, true_cor_XX = R, true_cor_YX = r,
                               which_predictor = 1, assurance = .9))
  for (f in list(ss_aipe_reg_coef_sensitivity, ss_aipe_rc_sensitivity)) {
    expect_no_warning(f(true_var_Y = 1, true_cov_YX = r, true_cov_XX = R,
                        width = .3, assurance = .9, G = 2, print_iter = FALSE,
                        seed = 1))
  }
  # Ten exchangeable predictors: the exact plan and the chi square step are
  # both 244.
  R10 <- matrix(.2, 10, 10); diag(R10) <- 1
  expect_no_warning(expect_identical(ss_aipe_reg_coef(
    width = .28, true_cor_XX = R10, true_cor_YX = rep(.2, 10), which_predictor = 3,
    assurance = .9)$value[1], 244))
  expect_identical(ss_aipe_src(
    width = .28, true_cor_XX = R10, true_cor_YX = rep(.2, 10), which_predictor = 3,
    noncentral = FALSE, assurance = .9)$value[1], 244)
})
