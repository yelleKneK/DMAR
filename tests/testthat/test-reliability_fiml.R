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
                            ci_method = "ml")
  o_fi <- reliability_omega(data = d, denominator = "model_implied",
                            missing = "fiml", ci_method = "ml")
  expect_equal(o_fi$value[o_fi$term == "estimate"],
               o_lw$value[o_lw$term == "estimate"], tolerance = 1e-8)

  r_lw <- reliability_omega(data = d, denominator = "observed",
                            ci_method = "none")
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
                              ci_method = "none")
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

test_that("N counts cases with at least one observed item; N_complete the complete cases", {
  skip_if_not_installed("lavaan")
  d <- .fiml_test_data(N = 100, seed = 30)
  d$y1[1:20] <- NA                      # 20 partially observed rows
  d[21:25, paste0("y", 1:4)] <- NA      # 5 rows with no observed item
  d$z[21] <- 0.5                        # aux observed does not rescue a row

  pick <- function(r, term) r$value[r$term == term]

  r <- reliability_alpha(data = d, aux = "z", ci_method = "none")
  expect_equal(pick(r, "N"), 95)
  expect_equal(pick(r, "N_complete"), 75)

  ro <- reliability_omega(data = d, aux = "z",
                          denominator = "model_implied",
                          ci_method = "none")
  expect_equal(pick(ro, "N"), 95)
  expect_equal(pick(ro, "N_complete"), 75)

  # Listwise: N and N_complete coincide, and both count complete cases
  # on the items.
  r_lw <- reliability_alpha(data = d[paste0("y", 1:4)],
                            ci_method = "none")
  expect_equal(pick(r_lw, "N"), 75)
  expect_equal(pick(r_lw, "N_complete"), 75)

  # On fully complete data the two are equal under fiml as well.
  d_c <- .fiml_test_data(N = 80, seed = 31)[paste0("y", 1:4)]
  r_c <- reliability_alpha(data = d_c, missing = "fiml",
                           ci_method = "none")
  expect_equal(pick(r_c, "N"), 80)
  expect_equal(pick(r_c, "N_complete"), 80)
})

test_that("missing and aux travel as attributes", {
  skip_if_not_installed("lavaan")
  d <- .fiml_test_data_mar(N = 150, seed = 32)

  r <- reliability_alpha(data = d, aux = "z", ci_method = "none")
  expect_identical(attr(r, "missing"), "fiml")
  expect_identical(attr(r, "aux"), "z")

  r_lw <- reliability_alpha(data = d[paste0("y", 1:4)],
                            ci_method = "none")
  expect_identical(attr(r_lw, "missing"), "listwise")
  expect_null(attr(r_lw, "aux"))
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
  expect_error(reliability(data = d, type = "omega_categorical",
                           aux = "z"),
               "apply only to")
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
  expect_true(is.finite(pick("lower_limit")))
  expect_true(is.finite(pick("upper_limit")))
  expect_lt(pick("lower_limit"), pick("estimate"))
  expect_gt(pick("upper_limit"), pick("estimate"))
  expect_identical(attr(r, "B"), 60)
  expect_identical(attr(r, "missing"), "fiml")
})
