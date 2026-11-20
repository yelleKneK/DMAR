# The bootstraps in this file are mechanics tests run with a small B; the
# rule on the number of replications is tested in test-dmar_bootstrap.R.
local_small_B()

## Full information maximum likelihood and auxiliary variables in the
## reliability family: reliability_alpha() and reliability_omega() with
## missing = c("listwise", "fiml") and aux.
##
## The oracle for the saturated correlates model is a lavaan model
## specified from scratch in these tests (not a third-party auxiliary
## variable wrapper, and not MBESS).

# Simulate congeneric items (with nonzero means, so a wrong FIML mean
# structure could not hide) plus an auxiliary correlated with the factor.
.fiml_test_data <- function(N = 300, lam = c(0.8, 0.7, 0.6, 0.75),
                            seed = 4258) {
  set.seed(seed)
  J <- length(lam)
  psi <- 1 - lam^2
  eta <- rnorm(N)
  X <- outer(eta, lam) + matrix(rnorm(N * J), N, J) %*% diag(sqrt(psi)) + 2
  colnames(X) <- paste0("y", seq_len(J))
  data.frame(X, z = eta + rnorm(N, sd = 0.7))
}

# The same design with MAR missingness on y2 driven by the auxiliary
# (plus some MCAR missingness on y4).
.fiml_test_data_mar <- function(N = 300, seed = 4258) {
  d <- .fiml_test_data(N = N, seed = seed)
  p_miss <- stats::plogis(-1 + 1.5 * as.numeric(scale(d$z)))
  d$y2[stats::runif(N) < p_miss] <- NA
  d$y4[stats::runif(N) < 0.15] <- NA
  d
}


## ---------------------------------------------------------------------------
## (1) With complete data, fiml reduces to listwise.
## ---------------------------------------------------------------------------

test_that("with complete data, missing = 'fiml' equals listwise to 1e-8 (alpha, both estimators)", {
  skip_if_not_installed("lavaan")
  d <- .fiml_test_data()[paste0("y", 1:4)]

  a_lw <- reliability_alpha(data = d, ci_method = "none")
  a_fi <- reliability_alpha(data = d, missing = "fiml", ci_method = "none")
  expect_equal(a_fi$value[a_fi$term == "estimate"],
               a_lw$value[a_lw$term == "estimate"], tolerance = 1e-8)

  m_lw <- reliability_alpha(data = d, estimator = "model_implied",
                            ci_method = "ml")
  m_fi <- reliability_alpha(data = d, estimator = "model_implied",
                            missing = "fiml", ci_method = "ml")
  expect_equal(m_fi$value[m_fi$term == "estimate"],
               m_lw$value[m_lw$term == "estimate"], tolerance = 1e-8)
  # The standard errors are not bitwise equal on complete data because
  # the FIML fit uses the observed information matrix where the
  # complete-data ML fit uses the expected one; they agree to about
  # three decimals here.
  expect_equal(m_fi$value[m_fi$term == "se"],
               m_lw$value[m_lw$term == "se"], tolerance = 5e-3)
})

test_that("with complete data, missing = 'fiml' equals listwise to 1e-8 (omega, both denominators)", {
  skip_if_not_installed("lavaan")
  d <- .fiml_test_data()[paste0("y", 1:4)]

  o_lw <- reliability_omega(data = d, denominator = "model_implied",
                            missing = "listwise", ci_method = "ml")
  o_fi <- reliability_omega(data = d, denominator = "model_implied",
                            missing = "fiml", ci_method = "ml")
  expect_equal(o_fi$value[o_fi$term == "estimate"],
               o_lw$value[o_lw$term == "estimate"], tolerance = 1e-8)

  r_lw <- reliability_omega(data = d, denominator = "observed",
                            missing = "listwise", ci_method = "none")
  r_fi <- reliability_omega(data = d, denominator = "observed",
                            missing = "fiml", ci_method = "none")
  expect_equal(r_fi$value[r_fi$term == "estimate"],
               r_lw$value[r_lw$term == "estimate"], tolerance = 1e-8)
})


## ---------------------------------------------------------------------------
## (2) Under MAR missingness driven by the auxiliary, FIML with the
##     auxiliary is less biased than listwise deletion.
## ---------------------------------------------------------------------------

test_that("FIML with the auxiliary is less biased than listwise under MAR", {
  skip_on_cran()
  skip_if_not_installed("lavaan")

  set.seed(4257)
  J <- 4
  lam <- rep(0.7, J)
  psi <- 1 - lam^2
  Sigma_pop <- tcrossprod(lam) + diag(psi)
  alpha_pop <- (J / (J - 1)) * (1 - sum(diag(Sigma_pop)) / sum(Sigma_pop))

  G <- 200
  N <- 200
  bias_lw <- bias_fi <- numeric(G)
  for (g in seq_len(G)) {
    eta <- rnorm(N)
    X <- outer(eta, lam) + matrix(rnorm(N * J), N, J) %*% diag(sqrt(psi))
    colnames(X) <- paste0("y", seq_len(J))
    z <- eta + rnorm(N, sd = 0.5)
    d <- data.frame(X, z = z)
    p_miss <- stats::plogis(-0.6 + 1.6 * as.numeric(scale(z)))
    d$y1[stats::runif(N) < p_miss] <- NA
    d$y3[stats::runif(N) < p_miss] <- NA
    r_lw <- reliability_alpha(data = d[paste0("y", seq_len(J))],
                              missing = "listwise", ci_method = "none")
    r_fi <- reliability_alpha(data = d, aux = "z", ci_method = "none")
    bias_lw[g] <- r_lw$value[r_lw$term == "estimate"] - alpha_pop
    bias_fi[g] <- r_fi$value[r_fi$term == "estimate"] - alpha_pop
  }
  # With seed 4257 and G = 200 replications: listwise mean bias -0.0729
  # (Monte Carlo SE 0.0035), FIML-with-auxiliary mean bias -0.0019
  # (Monte Carlo SE 0.0021), a roughly 38-fold reduction; the comparison
  # below leaves an order of magnitude of headroom over both MC SEs.
  expect_lt(abs(mean(bias_fi)), abs(mean(bias_lw)) / 5)
})


## ---------------------------------------------------------------------------
## (3) Against an independent oracle: a saturated correlates model
##     specified from scratch in lavaan.
## ---------------------------------------------------------------------------

test_that("a hand-built lavaan saturated-correlates model reproduces the estimates", {
  skip_if_not_installed("lavaan")
  d <- .fiml_test_data_mar()
  J <- 4

  # Congeneric measurement model with the auxiliary correlated with
  # every item residual and never loading on the factor, written out by
  # hand through lavaan::cfa() (a different interface from the
  # package's explicit lavaan::lavaan() call).
  oracle_model <- "
    f1 =~ NA*y1 + y1 + y2 + y3 + y4
    f1 ~~ 1*f1
    z ~~ z
    z ~~ y1 + y2 + y3 + y4
  "
  ofit <- lavaan::cfa(oracle_model, data = d, missing = "ml",
                      check.post = FALSE)
  pe <- lavaan::parameterEstimates(ofit)
  lams <- pe$est[pe$op == "=~"]
  psis <- pe$est[pe$op == "~~" & pe$lhs == pe$rhs &
                 pe$lhs %in% paste0("y", seq_len(J))]
  omega_oracle <- sum(lams)^2 / (sum(lams)^2 + sum(psis))

  r <- reliability_omega(data = d, aux = "z",
                         denominator = "model_implied", ci_method = "ml")
  expect_equal(r$value[r$term == "estimate"], omega_oracle,
               tolerance = 1e-6)

  # The analytic path: coefficient alpha applied to the item block of
  # the FIML covariance matrix from a hand-built saturated model over
  # the items plus the auxiliary.
  vars <- c(paste0("y", seq_len(J)), "z")
  sat_model <- paste(
    c(paste(vars, "~~", vars),
      combn(vars, 2, function(v) paste(v[1], "~~", v[2]))),
    collapse = "\n")
  sfit <- lavaan::lavaan(sat_model, data = d, missing = "ml",
                         meanstructure = TRUE, int.ov.free = TRUE)
  Sig <- lavaan::lavInspect(sfit, "implied")$cov[paste0("y", seq_len(J)),
                                                 paste0("y", seq_len(J))]
  alpha_oracle <- (J / (J - 1)) * (1 - sum(diag(Sig)) / sum(Sig))

  ra <- reliability_alpha(data = d, aux = "z", ci_method = "ml")
  expect_equal(ra$value[ra$term == "estimate"], alpha_oracle,
               tolerance = 1e-6)
})


## ---------------------------------------------------------------------------
## (4) A pure-noise auxiliary does not disturb the measurement model.
## ---------------------------------------------------------------------------

test_that("a pure-noise auxiliary leaves the estimate unchanged", {
  skip_if_not_installed("lavaan")
  d <- .fiml_test_data(seed = 20)[paste0("y", 1:4)]
  set.seed(21)
  d$noise <- rnorm(nrow(d))
  items_only <- d[paste0("y", 1:4)]

  # Complete data: the saturated correlates block leaves the profile
  # likelihood of the measurement parameters untouched, so the match is
  # exact up to optimizer tolerance. An auxiliary wrongly entered as an
  # indicator (or correlated with the factor) would break this.
  m_no  <- reliability_alpha(data = items_only, estimator = "model_implied",
                             missing = "fiml", ci_method = "none")
  m_yes <- reliability_alpha(data = d, aux = "noise",
                             estimator = "model_implied",
                             ci_method = "none")
  expect_equal(m_yes$value[m_yes$term == "estimate"],
               m_no$value[m_no$term == "estimate"], tolerance = 1e-6)

  # Incomplete (MCAR) data: the two estimators no longer coincide
  # exactly, but a noise auxiliary carries no information, so the
  # difference is far inside sampling noise (about 6e-5 here).
  set.seed(77)
  d2 <- items_only
  d2$y1[sample(nrow(d2), 45)] <- NA
  d2$y3[sample(nrow(d2), 30)] <- NA
  d2_aux <- d2
  d2_aux$noise <- d$noise
  f_no  <- reliability_alpha(data = d2, missing = "fiml",
                             ci_method = "none")
  f_yes <- reliability_alpha(data = d2_aux, aux = "noise",
                             ci_method = "none")
  expect_equal(f_yes$value[f_yes$term == "estimate"],
               f_no$value[f_no$term == "estimate"], tolerance = 5e-3)
})


## ---------------------------------------------------------------------------
## (5) N and N_complete accounting.
## ---------------------------------------------------------------------------

test_that("N counts cases with an observed item or auxiliary variable; N_complete the complete cases", {
  skip_if_not_installed("lavaan")
  d <- .fiml_test_data(N = 100, seed = 30)
  d$y1[1:20] <- NA                      # 20 partially observed rows
  d[21:25, paste0("y", 1:4)] <- NA      # 5 rows with no observed item,
  d$z[25] <- NA                         # one of them with nothing observed

  pick <- function(r, term) r$value[r$term == term]

  # The four rows observed only on the auxiliary inform the item moments
  # through it and are used (they were set aside through 1.0.1); the row
  # with nothing observed carries nothing.
  r <- reliability_alpha(data = d, aux = "z", ci_method = "none")
  expect_equal(pick(r, "N"), 99)
  expect_equal(pick(r, "N_complete"), 75)

  ro <- reliability_omega(data = d, aux = "z",
                          denominator = "model_implied",
                          ci_method = "none")
  expect_equal(pick(ro, "N"), 99)
  expect_equal(pick(ro, "N_complete"), 75)

  rh <- reliability_H(data = d, aux = "z", ci_method = "none")
  expect_equal(pick(rh, "N"), 99)
  expect_equal(pick(rh, "N_complete"), 75)

  # Without an auxiliary a row with no observed item carries nothing.
  r_no <- reliability_alpha(data = d[paste0("y", 1:4)], missing = "fiml",
                            ci_method = "none")
  expect_equal(pick(r_no, "N"), 95)
  # Two-stage estimation uses the four rows observed only on the
  # auxiliary as FIML does, since its first stage reads them, and says so
  # in a message (through the first versions of 1.1.0 they were set
  # aside, and N was 95).
  expect_message(
    r_ts <- reliability_alpha(data = d, aux = "z", missing = "two_stage",
                              ci_method = "none"),
    "reliability_alpha(): 4 cases are observed on an auxiliary variable",
    fixed = TRUE)
  expect_equal(pick(r_ts, "N"), 99)
  expect_equal(pick(r_ts, "N_complete"), 75)

  # Listwise: N and N_complete coincide, and both count complete cases
  # on the items.
  r_lw <- reliability_alpha(data = d[paste0("y", 1:4)],
                            missing = "listwise", ci_method = "none")
  expect_equal(pick(r_lw, "N"), 75)
  expect_equal(pick(r_lw, "N_complete"), 75)

  # On fully complete data the two are equal under fiml as well.
  d_c <- .fiml_test_data(N = 80, seed = 31)[paste0("y", 1:4)]
  r_c <- reliability_alpha(data = d_c, missing = "fiml",
                           ci_method = "none")
  expect_equal(pick(r_c, "N"), 80)
  expect_equal(pick(r_c, "N_complete"), 80)
})

test_that("a case observed only on an auxiliary variable informs the estimates", {
  skip_if_not_installed("lavaan")
  # The whole item block goes missing at random given the auxiliary, the
  # case through 1.0.1 set aside: those analyses described the cases
  # with an observed item, whose z (and so factor scores) are lower.
  d <- .fiml_test_data(N = 4000, seed = 11)
  lam <- c(0.8, 0.7, 0.6, 0.75)
  Sig <- tcrossprod(lam) + diag(1 - lam^2)
  alpha_pop <- (4 / 3) * (1 - sum(diag(Sig)) / sum(Sig))
  omega_pop <- sum(lam)^2 / sum(Sig)
  items <- paste0("y", 1:4)
  set.seed(12)
  gone <- stats::runif(nrow(d)) < stats::plogis(-0.5 + 2 * as.numeric(
    scale(d$z)))
  d[gone, items] <- NA
  d$y1[stats::runif(nrow(d)) < 0.1] <- NA
  expect_gt(sum(gone), 1500)
  pick <- function(r, term) r$value[r$term == term]

  a <- reliability_alpha(data = d, aux = "z", ci_method = "ml")
  expect_equal(pick(a, "N"), nrow(d))
  # cov_fiml() keeps the same cases, and its matrix is alpha's basis.
  em <- cov_fiml(d, aux = "z")
  expect_identical(em$N, nrow(d))
  expect_equal(pick(a, "estimate"), DMAR:::.alpha_from_S(em$S),
               tolerance = 1e-12)
  # The estimate is near the population value; the analysis of the cases
  # with an observed item, what was computed through 1.0.1, is not (0.762
  # here, eight of its standard errors below).
  expect_lt(abs(pick(a, "estimate") - alpha_pop), 3 * pick(a, "se"))
  a_old <- reliability_alpha(data = d[!gone, ], aux = "z", ci_method = "ml")
  expect_gt(alpha_pop - pick(a_old, "estimate"), 5 * pick(a_old, "se"))

  # Model implied omega: lavaan's own FIML fit of the saturated correlates
  # model keeps the cases observed only on z.
  o <- reliability_omega(data = d, aux = "z",
                         denominator = "model_implied", ci_method = "ml")
  expect_equal(pick(o, "N"), nrow(d))
  ofit <- lavaan::cfa("
    f1 =~ NA*y1 + y1 + y2 + y3 + y4
    f1 ~~ 1*f1
    z ~~ z
    z ~~ y1 + y2 + y3 + y4
  ", data = d, missing = "ml", check.post = FALSE)
  expect_identical(lavaan::lavInspect(ofit, "nobs"), nrow(d))
  pe <- lavaan::parameterEstimates(ofit)
  lams <- pe$est[pe$op == "=~"]
  psis <- pe$est[pe$op == "~~" & pe$lhs == pe$rhs & pe$lhs %in% items]
  expect_equal(pick(o, "estimate"), sum(lams)^2 / (sum(lams)^2 + sum(psis)),
               tolerance = 1e-6)
  expect_lt(abs(pick(o, "estimate") - omega_pop), 3 * pick(o, "se"))
})

test_that("missing and aux travel as attributes", {
  skip_if_not_installed("lavaan")
  d <- .fiml_test_data_mar(N = 150, seed = 32)

  r <- reliability_alpha(data = d, aux = "z", ci_method = "none")
  expect_identical(attr(r, "missing"), "fiml")
  expect_identical(attr(r, "aux"), "z")

  r_lw <- reliability_alpha(data = d[paste0("y", 1:4)],
                            missing = "listwise", ci_method = "none")
  expect_identical(attr(r_lw, "missing"), "listwise")
  expect_null(attr(r_lw, "aux"))
  # Without 'missing' the incomplete items are analyzed by FIML (the
  # default since 1.1.0), and auxiliary variables are not implied.
  r_df <- suppressWarnings(reliability_alpha(data = d[paste0("y", 1:4)],
                                             ci_method = "none"))
  expect_identical(attr(r_df, "missing"), "fiml")
  expect_null(attr(r_df, "aux"))
})


## ---------------------------------------------------------------------------
## (6) Validation paths: every error fires.
## ---------------------------------------------------------------------------

test_that("aux with an explicit missing = 'listwise' is an error, not an override", {
  d <- .fiml_test_data(N = 50, seed = 40)
  expect_error(
    reliability_alpha(data = d, aux = "z", missing = "listwise"),
    "cannot be combined with"
  )
  expect_error(
    reliability_omega(data = d, aux = "z", missing = "listwise"),
    "cannot be combined with"
  )
})

test_that("fiml and aux reject covariance input", {
  d <- .fiml_test_data(N = 50, seed = 41)
  S <- cov(d[paste0("y", 1:4)])
  expect_error(
    reliability_alpha(S = S, N = 50, missing = "fiml"),
    "requires raw 'data'"
  )
  expect_error(
    reliability_alpha(S = S, N = 50, aux = "z"),
    "requires raw 'data'"
  )
  expect_error(
    reliability_omega(S = S, N = 50, missing = "fiml"),
    "requires raw 'data'"
  )
  expect_error(
    reliability_omega(S = S, N = 50, aux = "z"),
    "requires raw 'data'"
  )
})

test_that("aux names are validated against the data", {
  d <- .fiml_test_data(N = 50, seed = 42)
  expect_error(reliability_alpha(data = d, aux = "not_here"),
               "not found in 'data'")
  expect_error(reliability_alpha(data = d, aux = c("z", "z")),
               "duplicated")
  expect_error(reliability_alpha(data = d, aux = 3),
               "character vector")
  # Naming all but one column as auxiliary leaves too few items.
  expect_error(
    reliability_alpha(data = d, aux = c("z", paste0("y", 1:3))),
    "fewer than two item columns"
  )
  # Unnamed matrix input gets default V1, V2, ... names from
  # as.data.frame(), so an aux name that matched the original data frame
  # is reported as not found.
  d_un <- unname(as.matrix(d))
  expect_error(reliability_alpha(data = d_un, aux = "z"),
               "not found in 'data'")
})

test_that("a constant or non-numeric auxiliary is rejected", {
  d <- .fiml_test_data(N = 50, seed = 43)
  d$flat <- 1
  expect_error(reliability_alpha(data = d[c(paste0("y", 1:4), "flat")],
                                 aux = "flat"),
               "constant")
  d$label <- letters[seq_len(nrow(d)) %% 26 + 1]
  expect_error(reliability_alpha(data = d[c(paste0("y", 1:4), "label")],
                                 aux = "label"),
               "not numeric")
})

test_that("complete-data closed forms, ADF, and the profile likelihood refuse fiml", {
  d <- .fiml_test_data_mar(N = 100, seed = 44)
  for (m in c("feldt", "fisher", "bonett", "hakstian_whalen")) {
    expect_error(reliability_alpha(data = d, aux = "z", ci_method = m),
                 "complete-data closed form")
  }
  expect_error(reliability_alpha(data = d, aux = "z", ci_method = "adf"),
               "distribution-free")
  expect_error(reliability_alpha(data = d, aux = "z",
                                 estimator = "model_implied",
                                 ci_method = "likelihood"),
               "not been implemented")
  expect_error(reliability_omega(data = d, aux = "z",
                                 denominator = "model_implied",
                                 ci_method = "bonett"),
               "complete-data closed form")
  expect_error(reliability_omega(data = d, aux = "z",
                                 denominator = "model_implied",
                                 ci_method = "adf"),
               "distribution-free")
  expect_error(reliability_omega(data = d, aux = "z",
                                 denominator = "model_implied",
                                 ci_method = "likelihood"),
               "not been implemented")
})

test_that("the reliability() wrapper forwards missing and aux, and guards the other types", {
  skip_if_not_installed("lavaan")
  d <- .fiml_test_data_mar(N = 150, seed = 45)

  r_direct <- reliability_alpha(data = d, aux = "z", ci_method = "ml")
  r_wrap <- reliability(data = d, type = "alpha", aux = "z",
                        ci_method = "ml")
  expect_equal(r_wrap$value, r_direct$value, tolerance = 1e-12)
  expect_identical(attr(r_wrap, "missing"), "fiml")

  expect_error(reliability(data = d, type = "kr20", missing = "fiml"),
               "apply only to")
  expect_error(reliability(data = d, type = "kr20", aux = "z"),
               "apply only to")
  # Since 1.1.0 categorical omega takes auxiliary variables, under
  # multiple imputation, so the wrapper forwards them and the family
  # function makes its own check.
  expect_error(reliability(data = d, type = "omega_categorical",
                           aux = "z", missing = "pairwise"),
               "use missing = \"mi\"", fixed = TRUE)
})


## ---------------------------------------------------------------------------
## Unit tests for the saturated correlates syntax builder.
## ---------------------------------------------------------------------------

test_that(".aux_saturated_syntax builds the saturated correlates block", {
  syn <- DMAR:::.aux_saturated_syntax(paste0("y", 1:3), paste0("a", 1:2))
  lines <- strsplit(syn, "\n", fixed = TRUE)[[1]]
  expect_true("a1 ~~ a1" %in% lines)
  expect_true("a2 ~~ a2" %in% lines)
  expect_true("a1 ~~ a2" %in% lines)
  expect_true("a1 ~~ y1 + y2 + y3" %in% lines)
  expect_true("a2 ~~ y1 + y2 + y3" %in% lines)
  # The auxiliaries are never indicators of the factor.
  expect_false(any(grepl("=~", lines, fixed = TRUE)))

  syn1 <- DMAR:::.aux_saturated_syntax("y1", "a1")
  lines1 <- strsplit(syn1, "\n", fixed = TRUE)[[1]]
  expect_identical(lines1, c("a1 ~~ a1", "a1 ~~ y1"))
})


## ---------------------------------------------------------------------------
## Bootstrap under fiml resamples the partially observed rows.
## ---------------------------------------------------------------------------

test_that("bootstrap intervals refit by FIML on resampled rows", {
  skip_on_cran()
  skip_if_not_installed("lavaan")
  skip_if_not_installed("boot")
  d <- .fiml_test_data_mar(N = 150, seed = 46)

  r <- reliability_alpha(data = d, aux = "z", ci_method = "percentile",
                         B = 60, seed = 113)
  pick <- function(term) r$value[r$term == term]
  expect_true(is.finite(pick("ci_lower")))
  expect_true(is.finite(pick("ci_upper")))
  expect_lt(pick("ci_lower"), pick("estimate"))
  expect_gt(pick("ci_upper"), pick("estimate"))
  expect_identical(attr(r, "B"), 60L)
  expect_identical(attr(r, "missing"), "fiml")
})


## ---------------------------------------------------------------------------
## The FIML covariance matrix is the EM estimate of the saturated model.
## ---------------------------------------------------------------------------

test_that("the FIML covariance matrix is lavaan's EM estimate, with and without an auxiliary", {
  skip_if_not_installed("lavaan")
  # The Pasteur students took neither t25 nor t26.
  sp <- c("t1_visual_perception", "t2_cubes", "t3_paper_form_board",
          "t4_lozenges", "t25_paper_form_board_r", "t26_flags")
  d <- holzinger_swineford[c(sp, "age_months")]
  # The unrestricted (h1) moments of any FIML fit to the same variables,
  # here a single-factor model, and with the auxiliary as a saturated
  # correlate written out by hand. The EM behind .fiml_cov() runs on the
  # standardized variables (cov_fiml(), so that its stopping rule does not
  # depend on their units) and lavaan's in their raw metric, and both stop
  # at lavaan's criterion, so they agree to about seven digits.
  measurement <- paste("f =~", paste(sp, collapse = " + "))
  fit <- lavaan::cfa(measurement, data = d[sp], missing = "ml")
  h1 <- unclass(lavaan::lavInspect(fit, "sampstat.h1")$cov)[sp, sp]
  fc <- DMAR:::.fiml_cov(d, items = sp)
  expect_true(fc$converged)
  expect_null(fc$fit)
  expect_identical(fc$N, 301L)
  expect_equal(fc$S, h1, tolerance = 1e-6, ignore_attr = TRUE)
  expect_identical(dimnames(fc$S), list(sp, sp))

  fit_aux <- lavaan::cfa(paste0(
    measurement, "\nage_months ~~ age_months\nage_months ~~ ",
    paste(sp, collapse = " + ")), data = d, missing = "ml",
    check.post = FALSE)
  h1_aux <- unclass(lavaan::lavInspect(fit_aux, "sampstat.h1")$cov)[sp, sp]
  fc_aux <- DMAR:::.fiml_cov(d, items = sp, aux = "age_months")
  expect_equal(fc_aux$S, h1_aux, tolerance = 1e-6, ignore_attr = TRUE)

  # The estimate does not depend on the interval: the saturated fit behind
  # the standard error is evaluated at the EM estimates, with no iteration,
  # and its standard error is that of lavaan's own fit of the saturated
  # model.
  a_ml <- reliability_alpha(d[sp], missing = "fiml", ci_method = "ml")
  a_none <- reliability_alpha(d[sp], missing = "fiml", ci_method = "none")
  expect_identical(a_ml$value[a_ml$term == "estimate"],
                   a_none$value[a_none$term == "estimate"])
  fs <- DMAR:::.fiml_cov(d, items = sp, se = "standard")
  expect_identical(lavaan::lavInspect(fs$fit, "iterations"), 0L)
  y <- paste0("y", seq_along(sp))
  expect_equal(unname(unclass(lavaan::lavInspect(fs$fit, "implied")$cov)[y, y]),
               unname(fs$S), tolerance = 1e-14)
  pe <- lavaan::parameterEstimates(fs$fit)
  expect_equal(pe$est[pe$lhs == "alpha"], fs$alpha, tolerance = 1e-12)
  sat <- paste(c(paste(sp, "~~", sp),
                 utils::combn(sp, 2, function(v) paste(v[1], "~~", v[2]))),
               collapse = "\n")
  sfit <- lavaan::lavaan(sat, data = d[sp], missing = "ml",
                         meanstructure = TRUE, int.ov.free = TRUE)
  V <- lavaan::lavInspect(sfit, "vcov")
  S_hat <- unclass(lavaan::lavInspect(sfit, "implied")$cov)[sp, sp]
  # The delta method for alpha written out: its gradient with respect to
  # the variances and covariances.
  J <- length(sp)
  tot <- sum(S_hat)
  tr <- sum(diag(S_hat))
  grad_var <- (J / (J - 1)) * (-1 / tot + tr / tot^2)
  grad_cov <- (J / (J - 1)) * (2 * tr / tot^2)
  pt <- lavaan::parTable(sfit)
  free <- pt[pt$free > 0L, ]
  g <- ifelse(free$op == "~~" & free$lhs == free$rhs, grad_var,
              ifelse(free$op == "~~", grad_cov, 0))
  expect_equal(a_ml$value[a_ml$term == "se"],
               sqrt(as.numeric(t(g) %*% V %*% g)), tolerance = 1e-5)
  expect_identical(a_ml$value[a_ml$term == "se"], fs$alpha_se)

  # Complete data: the EM estimate is the sample covariance matrix with
  # divisor N, so the classical alpha is the listwise one.
  d_c <- holzinger_swineford[sp[1:4]]
  a_fi <- reliability_alpha(d_c, missing = "fiml", ci_method = "none")
  a_lw <- reliability_alpha(d_c, ci_method = "none")
  expect_equal(a_fi$value[a_fi$term == "estimate"],
               a_lw$value[a_lw$term == "estimate"], tolerance = 1e-12)
})

test_that("26 incomplete items are analyzed without a warning, at the maximum", {
  skip_on_cran()
  skip_if_not_installed("lavaan")
  # Every Holzinger and Swineford test; the Pasteur students lack t25 and
  # t26. Through 1.0.1 the saturated model was fitted from lavaan's
  # default starting values, which took 2484 iterations here and stopped
  # with covariances up to 0.0105 from the maximum likelihood values,
  # and on other data passed lavaan's warnings that its initial model
  # implied matrix was not positive definite to the user.
  items <- grep("^t[0-9]+_", names(holzinger_swineford), value = TRUE)
  expect_length(items, 26L)
  d <- holzinger_swineford[items]
  em <- cov_fiml(d)
  expect_no_warning(a <- reliability_alpha(d, missing = "fiml",
                                           ci_method = "ml"))
  expect_equal(a$value[a$term == "estimate"], DMAR:::.alpha_from_S(em$S),
               tolerance = 1e-12)
  expect_true(is.finite(a$value[a$term == "se"]))
  expect_no_warning(o <- suppressMessages(
    reliability_omega(d, missing = "fiml", denominator = "observed")))
  expect_true(is.finite(o$value[o$term == "estimate"]))
})

test_that("an EM stopped at its iteration limit is reported, not used", {
  skip_if_not_installed("lavaan")
  d <- .fiml_test_data_mar()
  items <- paste0("y", 1:4)
  expect_true(DMAR:::.fiml_cov(d, items = items, aux = "z")$converged)
  # Two iterations of lavaan's EM: its last iterate is not the maximum,
  # and the saturated model is reported as not converged, as a fit that
  # stopped short was through 1.0.1.
  em <- DMAR:::.cov_fiml_em
  testthat::local_mocked_bindings(
    .cov_fiml_em = function(x, em_args = NULL) {
      em(x, em_args = list(max_iter = 2L))
    })
  fc <- DMAR:::.fiml_cov(d, items = items, aux = "z", se = "standard")
  expect_false(fc$converged)
  expect_null(fc$S)
  expect_null(fc$fit)
  expect_error(reliability_alpha(data = d, aux = "z", ci_method = "ml"),
               "did not converge")
})


## ---------------------------------------------------------------------------
## A pair of variables never observed together is refused by name.
## ---------------------------------------------------------------------------

test_that("unrestricted-moment routes refuse pairs absent from the data while constrained FIML can remain", {
  skip_if_not_installed("lavaan")
  # y1 is observed in the first half of the cases and y2 in the second, so
  # their covariance is not in the data. Through 1.0.1 lavaan only warned,
  # and analytic alpha came back with a standard error of zero and a
  # zero-width interval around an arbitrary estimate.
  d <- .fiml_test_data()[paste0("y", 1:4)]
  d$y1[1:150] <- NA
  d$y2[151:300] <- NA
  msg <- "never observed together"
  calls <- list(
    function() reliability_alpha(d, missing = "fiml"),
    function() reliability_alpha(d, missing = "two_stage"),
    function() reliability_omega(d, missing = "fiml"),
    function() reliability_omega(d, missing = "two_stage",
                                 denominator = "model_implied"))
  for (f in calls) {
    expect_error(f(), msg)
    expect_error(f(), "y1 with y2", fixed = TRUE)
    expect_error(f(), "missing_coverage()", fixed = TRUE)
  }
  # The constrained models have overlapping indicators and need not
  # estimate the absent unrestricted covariance. Their H1 qualification
  # remains, and each interval describes its own model-implied target.
  h0 <- list(
    function() reliability_alpha(d, missing = "fiml", estimator = "model_implied"),
    function() reliability_omega(d, missing = "fiml", denominator = "model_implied"),
    function() reliability_H(data = d, missing = "fiml"))
  for (f in h0) {
    warnings <- testthat::capture_warnings(result <- f())
    expect_match(warnings, "Unrestricted-H1 diagnostics", fixed = TRUE, all = FALSE)
    expect_true(all(is.finite(result$value[result$term %in%
      c("estimate", "ci_lower", "ci_upper")])))
    expect_match(attr(result, "h1_diagnostic"), "y1 and y2", fixed = TRUE)
  }
  equal <- suppressWarnings(lavaan::cfa(
    "f =~ l*y1 + l*y2 + l*y3 + l*y4", data = d, std.lv = TRUE, missing = "ml"))
  expected <- .alpha_from_S(lavaan::fitted(equal)$cov)
  alpha <- suppressWarnings(h0[[1L]]())
  expect_equal(alpha$value[alpha$term == "estimate"], expected, tolerance = 1e-6)
  expect_error(reliability_alpha(d, missing = "two_stage"),
               "two-stage maximum likelihood needs", fixed = TRUE)
  # An auxiliary variable never observed with an item is refused too.
  d <- .fiml_test_data()
  d$z[seq(1, 300, by = 2)] <- NA
  d$y3[seq(2, 300, by = 2)] <- NA
  expect_error(reliability_alpha(d, aux = "z"), "y3 with z", fixed = TRUE)
  # Listwise deletion is unaffected: it uses the complete cases.
  d <- .fiml_test_data()[paste0("y", 1:4)]
  d$y1[1:10] <- NA
  expect_no_error(reliability_alpha(d, missing = "fiml"))
})

test_that("bootstrap pair screening names lost and single-valued joint observations", {
  skip_if_not_installed("lavaan")
  d <- .fiml_test_data(N = 120)[paste0("y", 1:4)]
  d$y1[62:120] <- NA; d$y2[1:59] <- NA
  for (m in c("fiml", "two_stage")) {
    checked <- DMAR:::.relia_checked_point(function(x) DMAR:::.fiml_cov(
      x, items = names(x), se = "none", missing = m), missing = m,
      items = names(d), fitted = TRUE, composite = TRUE)
    expect_error(checked(d[-c(60, 61), ]), "y1 with y2", fixed = TRUE)
    expect_error(checked(d[-60, ]), "y1 with y2", fixed = TRUE)
    # Saturated-moment estimates need more than this pairwise screen.
    expect_error(reliability_alpha(d, missing = m), class = "dmar_fiml_unbounded")
  }
})


test_that("a pair observed together only where one of them takes a single value is refused, and listwise deletion does not offer FIML for it", {
  skip_if_not_installed("lavaan")
  # Regression test. t2_cubes takes one value in the 90 rows that observe
  # t1_visual_perception and varies only in the 10 that do not, so the
  # rows that observe both say nothing about how the two vary together.
  # Listwise deletion refused t2 as constant in the 90 complete rows and
  # offered FIML, which returned robust omega 1.11 (a Heywood case: t1's
  # error variance -48 against an observed variance of 50), the analytic
  # alpha 0.511 with no standard error, and the model implied omega 0.890
  # [0.281, 1]. Version 1.0.1 returned robust omega 1.06, the analytic
  # alpha 0.113 with no interval, and the same model implied omega.
  items <- c("t1_visual_perception", "t2_cubes", "t4_lozenges",
             "t6_paragraph_comprehension")
  d <- holzinger_swineford[1:100, items]
  d$t2_cubes[1:90] <- 25
  d$t1_visual_perception[91:100] <- NA
  err <- tryCatch(reliability_alpha(d, missing = "listwise"),
                  error = conditionMessage)
  expect_match(err, paste0("the item varies in the rows it removes. Full ",
                           "information maximum likelihood (missing = ",
                           "\"fiml\") cannot help either: it would be ",
                           "refused for the pair t1_visual_perception with ",
                           "t2_cubes (t2_cubes takes one value in the 90 ",
                           "rows that observe both)"), fixed = TRUE)
  expect_false(grepl("uses those rows as well", err, fixed = TRUE))
  msg <- paste0("needs to estimate their covariance: t1_visual_perception ",
                "with t2_cubes (t2_cubes takes one value in the 90 rows ",
                "that observe both).")
  calls <- list(
    function() reliability_alpha(d, missing = "fiml"),
    function() reliability_alpha(d, missing = "fiml", ci_method = "ml"),
    function() reliability_alpha(d, missing = "two_stage"),
    function() reliability_omega(d, missing = "fiml"),
    function() reliability_omega(d, missing = "two_stage"),
    function() cov_fiml(d))
  for (f in calls) {
    expect_error(suppressMessages(f()), msg, fixed = TRUE)
  }
  # With t2's variation spread over rows that observe t1 too, the pair is
  # informed: listwise deletion offers FIML, and FIML computes.
  e <- holzinger_swineford[1:100, items]
  e$t2_cubes[1:90] <- 25
  e$t1_visual_perception[91:95] <- NA
  e$t4_lozenges[96:100] <- NA
  expect_error(reliability_alpha(e, missing = "listwise"),
               "the item varies in the rows it removes; missing = \"fiml\" uses those rows as well.",
               fixed = TRUE)
  a <- reliability_alpha(e, missing = "fiml")
  expect_equal(a$value[a$term == "estimate"], 0.5388108, tolerance = 1e-6)
  expect_true(all(is.finite(a$value[a$term %in% c("ci_lower",
                                                    "ci_upper")])))
  # An auxiliary variable is checked with the items.
  f <- holzinger_swineford[1:100, c(items, "t7_sentence")]
  f$t7_sentence[1:60] <- 20
  f$t2_cubes[61:100] <- NA
  expect_error(reliability_alpha(f, aux = "t7_sentence"),
               "t2_cubes with t7_sentence (t7_sentence takes one value in the 60 rows that observe both)",
               fixed = TRUE)
})


test_that("affine collapse refuses saturated-moment routes without a universal FIML count rule", {
  skip_if_not_installed("lavaan")
  items <- c("t1_visual_perception", "t2_cubes", "t4_lozenges",
             "t6_paragraph_comprehension", "t7_sentence")
  d <- holzinger_swineford[1:5, items]
  d[4, 1] <- NA; d[5, 2] <- NA
  certificate <- DMAR:::.fmi_unbounded_subset(d)
  expect_true(certificate$certified_by_dimension)
  expect_equal(certificate$n, 3)
  for (m in c("fiml", "two_stage")) {
    expect_error(reliability_alpha(d, missing = m), class = "dmar_fiml_unbounded")
    expect_error(reliability_omega(d, missing = m, denominator = "observed"),
                 class = "dmar_fiml_unbounded")
  }
  # These particular constrained fits fail their own convergence checks.
  # They are not refused by the obsolete saturated-model count rule.
  expect_error(suppressWarnings(reliability_alpha(d, missing = "fiml",
    estimator = "model_implied")), "did not converge")
  expect_error(suppressWarnings(reliability_H(data = d, missing = "fiml")),
               "did not converge")
  a <- holzinger_swineford[1:8, c(items, "t9_word_meaning", "t10_addition", "t11_code")]
  a[8, 1] <- NA
  expect_error(reliability_alpha(a, aux = c("t9_word_meaning", "t10_addition", "t11_code")),
               class = "dmar_fiml_unbounded")
  e <- holzinger_swineford[1:12, items]
  e[11, 1] <- NA; e[12, 2] <- NA
  expect_true(is.finite(reliability_alpha(e, missing = "fiml", ci_method = "none")$value[1]))
})


test_that("a refusal of listwise deletion does not offer FIML when FIML would refuse the composite", {
  skip_if_not_installed("lavaan")
  # Regression test. x1 is constant in the 20 complete rows, in which x2 +
  # x3 + x4 is zero, and 30 more rows each miss one of x2, x3, and x4.
  # Listwise deletion refused x1 and offered FIML, which then refused the
  # composite, constant in the 20 rows that observe every item (through
  # the first versions of 1.1.0 FIML returned alpha -0.586).
  set.seed(1)
  a <- matrix(stats::rnorm(40), 20)
  comp <- cbind(5, a, -rowSums(a))
  inc <- matrix(stats::rnorm(120), 30)
  for (k in 1:3) inc[(10 * (k - 1) + 1):(10 * k), k + 1] <- NA
  d <- as.data.frame(rbind(comp, inc))
  names(d) <- paste0("x", 1:4)
  cannot <- paste0("Full information maximum likelihood (missing = ",
                   "\"fiml\") cannot help either: the composite, the sum of ",
                   "the items, does not vary in the 20 rows that observe ",
                   "every item, so the likelihood of the saturated model ",
                   "it estimates first has no maximum.")
  for (f in list(function(m) reliability_alpha(d, missing = m),
                 function(m) reliability_omega(d, missing = m))) {
    err <- tryCatch(f("listwise"), error = conditionMessage)
    expect_match(err, paste0("x1 takes a single value in the 20 complete ",
                             "rows"), fixed = TRUE)
    expect_match(err, paste0("the item varies in the rows it removes. ",
                             cannot), fixed = TRUE)
    expect_false(grepl("uses those rows as well", err, fixed = TRUE))
    expect_error(suppressWarnings(f("fiml")),
                 "does not vary in the 20 rows that observe every item",
                 fixed = TRUE)
  }
  # H weights its items, and FIML does not refuse its composite, so the
  # offer stands.
  expect_error(reliability_H(data = d, missing = "listwise"),
               "missing = \"fiml\" uses those rows as well", fixed = TRUE)
  # The refusal of too few distinct complete rows for a fitted coefficient
  # follows the same rule: three ipsative rows, each given twice, and 30
  # incomplete rows.
  v <- c("t1_visual_perception", "t2_cubes", "t4_lozenges")
  ip <- as.matrix(holzinger_swineford[1:3, v])
  ip <- ip - rowMeans(ip)
  rest <- as.matrix(holzinger_swineford[4:33, v])
  rest[cbind(1:30, rep(1:3, 10))] <- NA
  e <- as.data.frame(rbind(ip, ip, rest))
  err <- tryCatch(reliability_omega(e, missing = "listwise"),
                  error = conditionMessage)
  expect_match(err, "Only 6 rows of 'data' are observed on every item, of which 3 are distinct",
               fixed = TRUE)
  expect_match(err, sub("the 20 rows", "the 6 rows", cannot, fixed = TRUE),
               fixed = TRUE)
  expect_error(suppressWarnings(reliability_omega(e, missing = "fiml")),
               "does not vary in the 6 rows that observe every item",
               fixed = TRUE)
})


test_that("a refusal of listwise deletion offers FIML with two columns of one name only on the condition that they are named apart", {
  skip_if_not_installed("lavaan")
  # Regression test. Listwise deletion without auxiliary variables reads the
  # columns by position and FIML by name, refusing a name given to two
  # columns before anything else; the refusals of listwise deletion offered
  # FIML without a word about the name, and FIML then refused the call.
  hs <- holzinger_swineford
  d <- hs[1:100, c("t1_visual_perception", "t2_cubes", "t4_lozenges",
                   "t6_paragraph_comprehension")]
  d$t2_cubes[1:90] <- 25
  d$t1_visual_perception[91:95] <- NA
  d$t4_lozenges[96:100] <- NA
  names(d)[4] <- "t1_visual_perception"
  # Listwise deletion is asked for by name, since the default analyzes
  # incomplete items by FIML, which refuses the name at once.
  err <- tryCatch(reliability_alpha(d, missing = "listwise"),
                  error = conditionMessage)
  expect_match(err, paste0("missing = \"fiml\" uses those rows as well. ",
                           "It selects the items by name, so first give the ",
                           "columns named 't1_visual_perception' distinct ",
                           "names (or leave out a column that copies ",
                           "another)."), fixed = TRUE)
  expect_error(reliability_alpha(d, missing = "fiml"),
               "The item names must be distinct", fixed = TRUE)
  renamed <- d
  names(renamed)[4] <- "t6_paragraph_comprehension"
  expect_true(is.finite(reliability_alpha(renamed, missing = "fiml",
                                          ci_method = "none")$value[1]))
  # A data obstacle is still reported, with each column of the shared name
  # given its position.
  set.seed(1)
  e <- data.frame(x1 = 5, x2 = stats::rnorm(40), x3 = stats::rnorm(40),
                  x2 = stats::rnorm(40), check.names = FALSE)
  e$x1[31:40] <- stats::rnorm(10)
  e$x3[31:40] <- NA
  err <- tryCatch(reliability_alpha(e, missing = "listwise"),
                  error = conditionMessage)
  expect_match(err, "cannot help either: it would be refused for the pair x1 with x3",
               fixed = TRUE)
  expect_no_match(err, "distinct names", fixed = TRUE)
  # KR-20 offers reliability_alpha() on the same condition.
  set.seed(2)
  k <- as.data.frame(matrix(stats::rbinom(240, 1, .5), 60, 4))
  names(k) <- c("k1", "k2", "k3", "k3")
  k[1:20, 1] <- NA
  k[21:40, 2] <- NA
  k[41:60, 3] <- NA
  expect_error(reliability_kr20(k),
               "It selects the items by name, so first give the columns named 'k3' distinct names",
               fixed = TRUE)
  # The fitted coefficients refuse the name under listwise deletion too, so
  # they now do so first, before the refusal of the rows.
  d2 <- hs[1:8, c("t1_visual_perception", "t2_cubes", "t4_lozenges",
                  "t6_paragraph_comprehension", "t7_sentence")]
  d2[3:8, 1] <- NA
  names(d2)[5] <- "t2_cubes"
  expect_error(reliability_omega(d2, missing = "listwise"),
               "The item names must be distinct; 't2_cubes' names more than one column.",
               fixed = TRUE)
})


test_that("an H1 refusal does not promise a listwise analysis with undefined input", {
  skip_if_not_installed("lavaan")
  items <- c("t1_visual_perception", "t2_cubes", "t4_lozenges",
             "t6_paragraph_comprehension", "t7_sentence")
  d <- holzinger_swineford[1:5, items]
  d[4:5, 1] <- NA; d$t2_cubes[1:3] <- 25
  for (m in c("fiml", "two_stage")) {
    expect_error(reliability_alpha(d, missing = m), "t2_cubes takes one value in the 3 rows", fixed = TRUE)
  }
  expect_error(reliability_alpha(d, missing = "listwise"),
               "cannot help either", fixed = TRUE)
  a <- data.frame(x1 = c(3, 5, 2, NA), x2 = c(4, 1, NA, 6),
                  x3 = c(2, 6, 7, 1), x4 = c(6, 3, 4, 2))
  expect_error(reliability_alpha(a, missing = "fiml"), class = "dmar_fiml_unbounded")
  expect_error(reliability_alpha(a, missing = "listwise"), "does not vary",
               fixed = TRUE)
  two <- holzinger_swineford[1:5, items]
  two[3:5, 1] <- NA
  expect_error(reliability_alpha(two, missing = "fiml"), class = "dmar_fiml_unbounded")
  expect_warning(x <- reliability_alpha(two, ci_method = "feldt",
                                        missing = "listwise"),
    "The estimate of coefficient alpha is negative", fixed = TRUE)
  S <- cov(two[complete.cases(two), ])
  expect_equal(x$value[1], ncol(S) / (ncol(S) - 1) * (1 - sum(diag(S)) / sum(S)))
})


test_that("the case count behind the FIML and two-stage refusal counts distinct cases with an observed item, in the data and in every replication", {
  skip_if_not_installed("lavaan")
  skip_if_not_installed("boot")
  # Regression test. The rule compared the rows read with J + K, so a row
  # observed on an auxiliary variable alone, or a row repeated, lifted
  # the count past the bound without adding a dimension: four pupils of
  # three tests and an auxiliary variable were refused, and a fifth pupil
  # observed on the auxiliary alone gave the analytic alpha 0.418 with a
  # standard error of 0 and the interval [0.418, 0.418] under two-stage
  # estimation and FIML; five pupils of five tests were refused, and the
  # same rows given twice gave alpha 0.156 with a standard error of 0.
  v <- c("t1_visual_perception", "t2_cubes", "t4_lozenges",
         "t5_general_information")
  d <- holzinger_swineford[1:6, v]
  d[5:6, 1:3] <- NA
  for (m in c("two_stage", "fiml")) {
    for (k in 5:6) {
      err <- tryCatch(suppressMessages(
        reliability_alpha(d[1:k, ], aux = "t5_general_information",
                          missing = m, ci_method = "ml")),
        error = conditionMessage)
      expect_match(err, "4 rows observing", fixed = TRUE)
      expect_match(err, "affine-collapse direction", fixed = TRUE)
    }
    if (m == "two_stage") expect_error(suppressMessages(reliability_omega(
      d[1:5, ], aux = "t5_general_information", missing = m,
      denominator = "model_implied", ci_method = "ml")), class = "dmar_fiml_unbounded")
  }
  # One more case with an observed item computes, the case observed on
  # the auxiliary alone with it.
  e <- holzinger_swineford[1:7, v]
  e[7, 1:3] <- NA
  e[2, 1] <- NA
  expect_true(is.finite(suppressMessages(reliability_alpha(
    e, aux = "t5_general_information", missing = "two_stage",
    ci_method = "none"))$value[1]))
  # Repeated rows count once.
  items <- c("t1_visual_perception", "t2_cubes", "t4_lozenges",
             "t6_paragraph_comprehension", "t7_sentence")
  r <- holzinger_swineford[1:5, items]
  r[4, 1] <- NA
  r[5, 2] <- NA
  for (m in c("fiml", "two_stage")) {
    expect_error(reliability_alpha(rbind(r, r), missing = m),
                 class = "dmar_fiml_unstable")
  }
  expect_error(reliability_omega(rbind(r, r), missing = "two_stage"),
               class = "dmar_fiml_unstable")
  # Listwise deletion counts distinct complete rows for a fitted
  # coefficient, which lavaan otherwise refused in its own words.
  dc <- holzinger_swineford[1:5, items]
  expect_error(reliability_omega(rbind(dc, dc)),
               "'data' has only 10 rows, of which 5 are distinct",
               fixed = TRUE)
  # A replication follows the rule: at 12 cases and J + K = 8 most
  # resamples hold 8 or fewer distinct cases, and each of those is set
  # aside and counted, never analyzed.
  it <- c(items, "t9_word_meaning")
  ax <- c("t10_addition", "t11_code")
  b <- holzinger_swineford[1:12, c(it, ax)]
  b$t1_visual_perception[2] <- NA
  b$t7_sentence[8] <- NA
  res <- suppressWarnings(reliability_alpha(b, aux = ax,
                                            ci_method = "percentile",
                                            B = 60, seed = 1))
  bt <- attr(res, "bootstrap")
  set.seed(1)
  distinct <- boot::boot(b, function(x, i) length(unique(i)), R = 60)$t[, 1]
  expect_true(any(distinct <= 8L))
  expect_true(all(bt$status[distinct <= 8L] == "error"))
  expect_true(all(is.na(bt$t[distinct <= 8L])))
  expect_true(any(grepl("distinct cases of this resample have an observed item",
                        names(bt$error_messages), fixed = TRUE)))
  # A replication of a fitted coefficient under listwise deletion with no
  # more distinct rows than items is refused by the same rule.
  chk <- DMAR:::.relia_checked_point(function(x) 1, fitted = TRUE)
  expect_error(chk(rbind(dc, dc)), "Only 5 distinct rows of this resample",
               fixed = TRUE)
  expect_identical(DMAR:::.relia_checked_point(function(x) 1)(rbind(dc, dc)),
                   1)
})


test_that("no bootstrap interval is reported when more than half of the replications have too few distinct cases", {
  skip_if_not_installed("lavaan")
  skip_if_not_installed("boot")
  # Regression test. Once every replication followed the rule on the
  # number of cases, a resample of six pupils of five tests under FIML or
  # two-stage estimation (J + K + 1 cases) was kept only when it held all
  # six, the data reordered (6!/6^6 = 1.5 percent of resamples), so the
  # percentile interval [0.04347722, 0.04347727] had no width and lay
  # above the estimate 0.04347698, from 3 of 200 replications; BCa fell
  # back to the same limits. Before the rule reached the replications, the
  # other resamples were analyzed on singular EM estimates ([0, 0.898]).
  items <- c("t1_visual_perception", "t2_cubes", "t4_lozenges",
             "t6_paragraph_comprehension", "t7_sentence")
  d <- holzinger_swineford[1:6, items]
  S <- cov(d)
  expected <- ncol(S) / (ncol(S) - 1) * (1 - sum(diag(S)) / sum(S))
  for (m in c("fiml", "two_stage")) {
    for (ci in c("percentile", "bca")) {
      w <- testthat::capture_warnings(
        r <- reliability_alpha(d, missing = m, ci_method = ci, B = 200,
                               seed = 1))
      expect_equal(r$value[1], expected, tolerance = 1e-6)
      expect_true(all(is.na(r$value[2:4])))
      bs <- attr(r, "bootstrap")
      expect_identical(bs$returned, 3L)
      expect_gt(bs$too_few_cases, 100L)
      expect_lte(bs$too_few_cases, bs$reasons[["error"]])
      expect_match(w, paste0("had too few distinct cases for the ",
                             "coefficient, more than half of the 200, so ",
                             "no interval is reported"),
                   fixed = TRUE, all = FALSE)
      expect_false(any(grepl("Caution", w, fixed = TRUE)))
      expect_match(DMAR:::.bootstrap_footer_line(r),
                   paste0("; more than half had too few distinct cases, so ",
                          "no interval is reported$"))
    }
  }
  # Listwise deletion holds a fitted coefficient to more distinct rows
  # than items, and six complete pupils of five tests (J + 1) left the
  # omega interval [0.2240232, 0.2240232] (through the first versions of
  # 1.1.0, [0.023, 0.975] from 27 of 200, those lavaan fitted to singular
  # matrices among them).
  dc <- holzinger_swineford[1:6, items]
  for (f in list(
    function() reliability_omega(dc, ci_method = "percentile", B = 200,
                                 seed = 1),
    function() reliability_alpha(dc, estimator = "model_implied",
                                 ci_method = "percentile", B = 200,
                                 seed = 1))) {
    r <- suppressWarnings(f())
    expect_true(is.finite(r$value[1]))
    expect_true(all(is.na(r$value[2:4])))
    expect_gt(attr(r, "bootstrap")$too_few_cases, 100L)
  }
  # With twelve pupils, almost every resample holds more than five
  # distinct ones, and the interval is formed.
  d12 <- holzinger_swineford[1:12, items]
  d12[4, 1] <- NA
  d12[5, 2] <- NA
  r <- suppressWarnings(reliability_alpha(d12, missing = "fiml",
                                          ci_method = "percentile", B = 100,
                                          seed = 1))
  expect_lt(attr(r, "bootstrap")$too_few_cases, 50L)
  expect_true(all(is.finite(r$value[2:4])))
})


test_that("a standard error lavaan could not compute is NA, not zero", {
  skip_if_not_installed("lavaan")
  # Regression test. When the information matrix cannot be inverted,
  # lavaan gives a defined parameter a standard error of zero, which was
  # read as the standard error of alpha, so the Wald interval had no
  # width.
  v <- c("t1_visual_perception", "t2_cubes", "t4_lozenges",
         "t5_general_information")
  d <- holzinger_swineford[1:5, v]
  d[5, 1:3] <- NA
  fc <- suppressWarnings(DMAR:::.fiml_cov(d, items = v[1:3], aux = v[4],
                                          se = "standard"))
  expect_true(fc$converged)
  # This accessor retries the failed inversion and repeats lavaan's
  # singular-information warnings; the null covariance is the assertion.
  expect_null(suppressWarnings(lavaan::lavInspect(fc$fit, "vcov")))
  expect_identical(fc$alpha_se, NA_real_)
  ok <- DMAR:::.fiml_cov(holzinger_swineford[, v[1:3]], items = v[1:3],
                         se = "standard")
  expect_true(is.finite(ok$alpha_se) && ok$alpha_se > 0)
})

