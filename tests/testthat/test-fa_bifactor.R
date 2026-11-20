# fa_bifactor(): the bifactor model fit by fa_confirmatory()'s engine.
# Every number is checked against lavaan::cfa(std.lv = TRUE, orthogonal =
# TRUE) on the same items, against bifactor_indices() where the two
# describe the same quantity, and against the algebra of the indices
# written out here.

skip_if_not_installed("lavaan")

hs_groups <- list(
  spatial = c("t1_visual_perception", "t2_cubes", "t3_paper_form_board",
              "t4_lozenges"),
  verbal  = c("t5_general_information", "t6_paragraph_comprehension",
              "t7_sentence", "t8_word_classification", "t9_word_meaning"),
  memory  = c("t14_word_recognition", "t15_number_recognition",
              "t16_figure_recognition", "t17_object_number",
              "t18_number_figure", "t19_figure_word"))

# The lavaan syntax of the same model, the general factor on `items`.
bf_syntax <- function(groups, items = unlist(groups, use.names = FALSE),
                      general = "general") {
  paste(c(paste0(general, " =~ ", paste(items, collapse = " + ")),
          paste0(names(groups), " =~ ",
                 vapply(groups, paste, character(1L), collapse = " + "))),
        collapse = "\n")
}

# The indices of a lavaan bifactor fit, by hand: omega, omega_H, omega_HS
# for sums of the raw items; ECV from the standardized loadings; H with
# every other source of an item's variance as its error; the maximal
# reliability l' Sigma^-1 l; PUC from the structure.
bf_hand <- function(fit, groups, general = "general") {
  pe <- lavaan::parameterEstimates(fit)
  items <- lavaan::lavNames(fit, "ov")
  latent <- c(general, names(groups))
  L <- matrix(0, length(items), length(latent),
              dimnames = list(items, latent))
  ld <- pe[pe$op == "=~", ]
  L[cbind(ld$rhs, ld$lhs)] <- ld$est
  psi_rows <- pe$op == "~~" & pe$lhs == pe$rhs & pe$lhs %in% items
  psi <- stats::setNames(pe$est[psi_rows], pe$lhs[psi_rows])[items]
  v <- rowSums(L^2) + psi
  lg <- L[, general]
  common_s <- vapply(names(groups), function(s) sum(L[groups[[s]], s])^2,
                     numeric(1L))
  total <- sum(lg)^2 + sum(common_s) + sum(psi)
  Sigma <- tcrossprod(L) + diag(psi)
  out <- c(omega = (sum(lg)^2 + sum(common_s)) / total,
           omega_H = sum(lg)^2 / total,
           ECV = sum(lg^2 / v) / sum((rowSums(L^2)) / v),
           H = { t <- sum(lg^2 / (v - lg^2)); t / (1 + t) },
           maximal = sum(lg * solve(Sigma, lg)),
           PUC = 1 - sum(choose(lengths(groups), 2)) /
             choose(length(items), 2))
  for (s in names(groups)) {
    r <- groups[[s]]
    cg <- sum(lg[r])^2
    sub <- cg + common_s[[s]] + sum(psi[r])
    ls <- L[r, s]
    out[paste0(c("omega_", "omega_HS_", "ECV_", "H_", "maximal_"), s)] <- c(
      (cg + common_s[[s]]) / sub, common_s[[s]] / sub,
      sum(lg[r]^2 / v[r]) / sum((lg[r]^2 + ls^2) / v[r]),
      { t <- sum(ls^2 / (v[r] - ls^2)); t / (1 + t) },
      sum(L[, s] * solve(Sigma, L[, s])))
  }
  out
}

est <- function(x, term) x$estimate[match(term, x$term)]

test_that("the estimates and the fit are lavaan's bifactor model, under FIML and listwise deletion", {
  skip_on_cran()
  for (m in c("fiml", "listwise")) {
    x <- fa_bifactor(holzinger_swineford, hs_groups, missing = m)
    ref <- lavaan::cfa(bf_syntax(hs_groups), data = holzinger_swineford,
                       std.lv = TRUE, orthogonal = TRUE,
                       missing = if (m == "fiml") "ml" else "listwise")
    pe <- lavaan::parameterEstimates(ref)
    ld <- pe[pe$op == "=~", ]
    items <- unlist(hs_groups, use.names = FALSE)
    # The general loadings are numbered over the items, a group factor's
    # within it.
    expect_equal(est(x, paste0("lambda_general_", seq_along(items))),
                 ld$est[ld$lhs == "general"], tolerance = 1e-5)
    expect_equal(x$se[match(paste0("lambda_general_", seq_along(items)),
                            x$term)],
                 ld$se[ld$lhs == "general"], tolerance = 1e-4)
    for (s in names(hs_groups)) {
      expect_equal(est(x, paste0("lambda_", s, "_", seq_along(hs_groups[[s]]))),
                   ld$est[ld$lhs == s], tolerance = 1e-5)
      psi <- pe[pe$op == "~~" & pe$lhs == pe$rhs &
                  pe$lhs %in% hs_groups[[s]], ]
      expect_equal(est(x, paste0("psi_", s, "_", seq_along(hs_groups[[s]]))),
                   psi$est, tolerance = 1e-5)
    }
    expect_equal(est(x, "chi_square"),
                 unname(lavaan::fitMeasures(ref, "chisq")), tolerance = 1e-5)
    expect_equal(est(x, "df"), unname(lavaan::fitMeasures(ref, "df")))
    # Intercepts join the table under FIML, as in fa_confirmatory().
    expect_identical(any(startsWith(x$term, "nu_")), m == "fiml")
    expect_identical(attr(x, "missing"), m)
    # The indices are the algebra above, from lavaan's own fit.
    h <- bf_hand(ref, hs_groups)
    expect_equal(est(x, c("omega_general", "omega_H_general", "ECV_general",
                          "H_general", "maximal_reliability_general",
                          "PUC_general")),
                 unname(h[c("omega", "omega_H", "ECV", "H", "maximal",
                            "PUC")]), tolerance = 1e-5)
    for (s in names(hs_groups)) {
      expect_equal(est(x, paste0(c("omega_", "omega_HS_", "ECV_", "H_",
                                   "maximal_reliability_"), s)),
                   unname(h[paste0(c("omega_", "omega_HS_", "ECV_", "H_",
                                     "maximal_"), s)]), tolerance = 1e-5)
    }
  }
})

test_that("the omega rows are the values semTools::compRelSEM() gives, and ECV, H, the maximal reliability, and PUC are bifactor_indices()'s", {
  skip_on_cran()
  x <- fa_bifactor(holzinger_swineford, hs_groups, output = "measurement")
  fit <- fa_bifactor(holzinger_swineford, hs_groups, output = "fit")
  expect_s4_class(fit, "lavaan")
  b <- bifactor_indices(fit)
  latent <- c("general", names(hs_groups))
  expect_equal(est(x, paste0("ECV_", latent)), b$ECV, tolerance = 1e-6)
  expect_equal(est(x, paste0("H_", latent)), b$H, tolerance = 1e-6)
  expect_equal(est(x, paste0("maximal_reliability_", latent)),
               unname(attr(b, "maximal_reliability")), tolerance = 1e-6)
  expect_equal(est(x, "PUC_general"), b$PUC[1L])
  # The omega rows describe the sum of the raw items, bifactor_indices()'s
  # the sum of the standardized items: on these tests, whose variances
  # range from 8.0 to 153.3, omega_H is .552 against .577.
  expect_equal(round(est(x, "omega_H_general"), 3), 0.552)
  expect_equal(round(b$omega_H[1L], 3), 0.577)
  # The standardized solution evaluates the defined indices at the
  # standardized parameters, which gives bifactor_indices()'s omegas.
  std <- fa_bifactor(holzinger_swineford, hs_groups, output = "standardized")
  sv <- function(t) std$est.std[match(t, std$lhs)]
  expect_equal(sv(paste0("omega_", latent)), b$omega, tolerance = 1e-6)
  expect_equal(sv("omega_H_general"), b$omega_H[1L], tolerance = 1e-6)
  expect_equal(sv(paste0("omega_HS_", names(hs_groups))), b$omega_HS[-1L],
               tolerance = 1e-6)
  # Pinned from semTools::compRelSEM(fit, obs.var = FALSE) on this fit,
  # without and with true = list(<factor> = "<factor>") (semTools 0.5.9,
  # lavaan 0.7.3, 2026-11-05): omega for the general factor and the three
  # subscales, then omega_H and the three omega_HS. They are the hand
  # algebra of bf_hand() on lavaan's own fit of the model, which agreed
  # with them to 1e-15, and the live comparison is
  # tools/oracle_snippets/semtools_bifactor.R.
  pinned_omega <- c(0.8755302000, 0.6951105684, 0.9023000439, 0.7283031141)
  pinned_hier <- c(0.5523387266, 0.2660922789, 0.6124102558, 0.4816024383)
  expect_equal(est(x, paste0("omega_", latent)), pinned_omega,
               tolerance = 1e-6)
  expect_equal(est(x, c("omega_H_general",
                        paste0("omega_HS_", names(hs_groups)))),
               pinned_hier, tolerance = 1e-6)
  ref <- lavaan::cfa(bf_syntax(hs_groups), data = holzinger_swineford,
                     std.lv = TRUE, orthogonal = TRUE, missing = "ml")
  h <- bf_hand(ref, hs_groups)
  expect_equal(unname(h[c("omega", paste0("omega_", names(hs_groups)))]),
               pinned_omega, tolerance = 1e-6)
  expect_equal(unname(h[c("omega_H", paste0("omega_HS_", names(hs_groups)))]),
               pinned_hier, tolerance = 1e-6)
})

test_that("the standard errors are the delta method's, and the limits are formed on the logit scale", {
  skip_on_cran()
  x <- fa_bifactor(holzinger_swineford, hs_groups, missing = "listwise")
  fit <- fa_bifactor(holzinger_swineford, hs_groups, missing = "listwise",
                     output = "fit")
  # A delta method written here over lavaan's coef() and vcov(), by label,
  # with numerical derivatives: it reproduces lavaan's standard error of H.
  theta <- lavaan::coef(fit)
  V <- lavaan::vcov(fit)
  items <- unlist(hs_groups, use.names = FALSE)
  latent <- c("general", names(hs_groups))
  # The labels of the loadings, items by factors, and of the error
  # variances, as fa_bifactor() writes them.
  lab_L <- matrix(NA_character_, length(items), length(latent))
  lab_L[, 1L] <- paste0("lambda_general_", seq_along(items))
  lab_psi <- character(length(items))
  for (s in names(hs_groups)) {
    r <- match(hs_groups[[s]], items)
    lab_L[r, match(s, latent)] <- paste0("lambda_", s, "_", seq_along(r))
    lab_psi[r] <- paste0("psi_", s, "_", seq_along(r))
  }
  on <- !is.na(lab_L)
  build <- function(th) {
    L <- matrix(0, length(items), length(latent))
    L[on] <- th[lab_L[on]]
    list(L = L, psi = unname(th[lab_psi]))
  }
  H_general <- function(th) {
    m <- build(th)
    v <- rowSums(m$L^2) + m$psi
    t <- sum(m$L[, 1L]^2 / (v - m$L[, 1L]^2))
    t / (1 + t)
  }
  dH <- vapply(seq_along(theta), function(i) {
    h <- 1e-6 * max(1, abs(theta[i]))
    up <- dn <- theta
    up[i] <- up[i] + h
    dn[i] <- dn[i] - h
    (H_general(up) - H_general(dn)) / (2 * h)
  }, numeric(1L))
  expect_equal(x$se[x$term == "H_general"],
               sqrt(drop(dH %*% V %*% dH)), tolerance = 1e-5)
  # The maximal reliability's standard error by a route independent of the
  # implementation's numerical derivatives over lavaan's parameter table:
  # the analytic gradient of rho_f = l_f' Sigma^-1 l_f with
  # Sigma = L L' + diag(psi). With a = Sigma^-1 l_f, d rho_f / d L_ik is
  # 2 a_i [k = f] - 2 a_i (L' a)_k, and d rho_f / d psi_i is -a_i^2, each
  # placed on the label of its parameter in coef().
  analytic_se <- function(f) {
    m <- build(theta)
    Sigma <- tcrossprod(m$L) + diag(m$psi)
    a <- solve(Sigma, m$L[, f])
    gL <- -2 * outer(a, drop(crossprod(m$L, a)))
    gL[, f] <- gL[, f] + 2 * a
    g <- stats::setNames(numeric(length(theta)), names(theta))
    g[lab_L[on]] <- gL[on]
    g[lab_psi] <- -a^2
    sqrt(drop(g %*% V %*% g))
  }
  for (k in seq_along(latent)) {
    expect_equal(x$se[x$term == paste0("maximal_reliability_", latent[k])],
                 analytic_se(k), tolerance = 1e-5)
  }
  # Limits on the logit scale, mapped back to (0, 1).
  for (t in c("omega_H_general", "ECV_general", "H_general",
              "maximal_reliability_general", "omega_HS_spatial")) {
    r <- x[x$term == t, ]
    eta <- stats::qlogis(r$estimate)
    d <- stats::qnorm(0.975) * r$se / (r$estimate * (1 - r$estimate))
    expect_equal(c(r$ci_lower, r$ci_upper), stats::plogis(eta + c(-d, d)),
                 tolerance = 1e-8)
  }
  # PUC depends on the structure alone: 1 - (6 + 10 + 15) / 105.
  r <- x[x$term == "PUC_general", ]
  expect_equal(r$estimate, 1 - 31 / 105)
  expect_true(all(is.na(unlist(r[c("se", "test_statistic", "p_value", "ci_lower",
                                   "ci_upper")]))))
  # A finite wald_df takes the t reference, for the indices on the logit
  # scale.
  y <- fa_bifactor(holzinger_swineford, hs_groups, missing = "listwise",
                   wald_df = 20, output = "measurement")
  r <- y[y$term == "maximal_reliability_general", ]
  d <- stats::qt(0.975, 20) * r$se / (r$estimate * (1 - r$estimate))
  expect_equal(c(r$ci_lower, r$ci_upper),
               stats::plogis(stats::qlogis(r$estimate) + c(-d, d)),
               tolerance = 1e-8)
  expect_equal(r$p_value, 2 * stats::pt(-abs(r$test_statistic), 20))
  expect_identical(attr(y, "wald_df"), 20)
})

test_that("items on the general factor alone join its rows and no subscale, and bifactor_indices() takes them", {
  skip_on_cran()
  alone <- c("t20_deduction", "t21_numerical_puzzles")
  x <- fa_bifactor(holzinger_swineford, hs_groups, general_only = alone,
                   missing = "listwise")
  items <- c(unlist(hs_groups, use.names = FALSE), alone)
  ref <- lavaan::cfa(bf_syntax(hs_groups, items), data = holzinger_swineford,
                     std.lv = TRUE, orthogonal = TRUE)
  expect_equal(est(x, "chi_square"),
               unname(lavaan::fitMeasures(ref, "chisq")), tolerance = 1e-5)
  # Their error variances are numbered as their general loadings.
  expect_true(all(c("lambda_general_16", "lambda_general_17",
                    "psi_general_16", "psi_general_17") %in% x$term))
  h <- bf_hand(ref, hs_groups)
  expect_equal(est(x, c("omega_general", "omega_H_general", "ECV_general",
                        "H_general", "maximal_reliability_general")),
               unname(h[c("omega", "omega_H", "ECV", "H", "maximal")]),
               tolerance = 1e-5)
  # PUC counts their pairs as uncontaminated: 1 - (6 + 10 + 15) / 136.
  expect_equal(est(x, "PUC_general"), 1 - 31 / 136)
  expect_identical(attr(x, "general_only"), alone)
  # bifactor_indices() took an item on no group factor as an error through
  # DMAR 1.0.x and 1.1.0 development; it now reads it as an item on the
  # general factor alone and agrees.
  b <- bifactor_indices(ref)
  latent <- c("general", names(hs_groups))
  expect_equal(est(x, paste0("ECV_", latent)), b$ECV, tolerance = 1e-5)
  expect_equal(est(x, paste0("H_", latent)), b$H, tolerance = 1e-5)
  expect_equal(est(x, "PUC_general"), b$PUC[1L])
  expect_equal(est(x, paste0("maximal_reliability_", latent)),
               unname(attr(b, "maximal_reliability")), tolerance = 1e-5)
})

test_that("covariance input, robust ML, and incomplete data follow fa_confirmatory()", {
  skip_on_cran()
  items <- unlist(hs_groups, use.names = FALSE)
  S <- stats::cov(holzinger_swineford[items])
  s_fit <- fa_bifactor(S = S, N = 301, factors = hs_groups)
  raw <- fa_bifactor(holzinger_swineford, hs_groups, missing = "listwise")
  keep <- seq_len(match("chi_square", raw$term) - 1L)
  expect_equal(est(s_fit, raw$term[keep]), raw$estimate[keep],
               tolerance = 1e-5)
  expect_equal(est(s_fit, c("chi_square", "df")),
               est(raw, c("chi_square", "df")), tolerance = 1e-5)
  expect_identical(attr(s_fit, "missing"), "none")
  # MLR: the estimates of ML, lavaan's Huber-White standard errors.
  mlr <- fa_bifactor(holzinger_swineford, hs_groups, estimator = "MLR")
  ref <- lavaan::cfa(bf_syntax(hs_groups), data = holzinger_swineford,
                     std.lv = TRUE, orthogonal = TRUE, estimator = "MLR",
                     missing = "ml")
  ld <- lavaan::parameterEstimates(ref)
  ld <- ld[ld$op == "=~" & ld$lhs == "general", ]
  expect_equal(mlr$se[match(paste0("lambda_general_", seq_along(items)),
                            mlr$term)], ld$se, tolerance = 1e-4)
  expect_identical(attr(mlr, "test"), "scaled")
  # Incomplete data: FIML by default, as lavaan's FIML fit, without a
  # warning.
  hs_inc <- holzinger_swineford
  set.seed(113)
  hs_inc$t7_sentence[sample(301, 60)] <- NA
  expect_no_warning(inc <- fa_bifactor(hs_inc, hs_groups))
  ref <- lavaan::cfa(bf_syntax(hs_groups), data = hs_inc, std.lv = TRUE,
                     orthogonal = TRUE, missing = "ml")
  expect_equal(est(inc, "chi_square"),
               unname(lavaan::fitMeasures(ref, "chisq")), tolerance = 1e-5)
  expect_identical(attr(inc, "missing"), "fiml")
  lw <- fa_bifactor(hs_inc, hs_groups, missing = "listwise")
  expect_identical(attr(lw, "missing"), "listwise")
  # The default with a standard error FIML cannot give, on incomplete data,
  # stops and names the options rather than deleting cases.
  expect_error(fa_bifactor(hs_inc, hs_groups, se = "robust.sem"),
               "no cases are deleted unasked")
})

test_that("a covariance matrix with means takes equal intercepts, as lavaan's model with shared labels", {
  skip_on_cran()
  items <- unlist(hs_groups, use.names = FALSE)
  S <- stats::cov(holzinger_swineford[items])
  M <- colMeans(holzinger_swineford[items])
  eq <- c(spatial = TRUE, verbal = FALSE, memory = FALSE)
  x <- fa_bifactor(S = S, N = 301, M = M, factors = hs_groups,
                   equal_intercept = eq)
  model <- paste(c(bf_syntax(hs_groups),
                   paste0(hs_groups$spatial, " ~ a*1")), collapse = "\n")
  ref <- lavaan::cfa(model, sample.cov = S, sample.mean = M,
                     sample.nobs = 301, std.lv = TRUE, orthogonal = TRUE,
                     meanstructure = TRUE)
  expect_equal(est(x, c("chi_square", "df")),
               unname(lavaan::fitMeasures(ref, c("chisq", "df"))),
               tolerance = 1e-5)
  expect_true("nu_spatial" %in% x$term)
  expect_equal(est(x, "nu_spatial"),
               lavaan::coef(ref)[["a"]], tolerance = 1e-5)
  expect_identical(attr(x, "equal")$intercept, "spatial")
  expect_true(attr(x, "has_means"))
  # The same model from the raw data, which are complete.
  raw <- fa_bifactor(holzinger_swineford, hs_groups, missing = "listwise",
                     meanstructure = TRUE, equal_intercept = eq)
  expect_equal(est(x, c("chi_square", "omega_H_general", "nu_spatial")),
               est(raw, c("chi_square", "omega_H_general", "nu_spatial")),
               tolerance = 1e-5)
  # Equal intercepts ask about the mean structure, which S alone lacks.
  expect_error(fa_bifactor(S = S, N = 301, factors = hs_groups,
                           equal_intercept = eq),
               "requires the item means 'M'")
})

# A fixture whose information matrix is singular at the estimates: six
# items in two groups of three, equal error variances in the first group,
# and data with equal general and group loadings there, so that the
# general and group loadings of that group are not separated. lavaan
# inverts the observed information of the FIML fit and the expected
# information of the listwise and covariance matrix fits without refusing.
bf_singular <- function() {
  L <- cbind(rep(0.6, 6), c(0.4, 0.4, 0.4, 0, 0, 0),
             c(0, 0, 0, 0.5, 0.3, 0.45))
  L[4:6, 1L] <- c(0.6, 0.5, 0.7)
  Sigma <- tcrossprod(L)
  diag(Sigma) <- 1
  set.seed(113)
  X <- as.data.frame(MASS::mvrnorm(400, rep(0, 6), Sigma, empirical = TRUE))
  names(X) <- paste0("x", 1:6)
  X
}

test_that("a singular information matrix leaves every standard error and interval NA, the maximal reliability's included, whatever the information", {
  skip_on_cran()
  X <- bf_singular()
  groups <- list(a = paste0("x", 1:3), b = paste0("x", 4:6))
  eq <- c(a = TRUE, b = FALSE)
  calls <- list(
    fiml = function() fa_bifactor(X, groups, equal_error = eq,
                                  output = "measurement"),
    listwise = function() fa_bifactor(X, groups, equal_error = eq,
                                      missing = "listwise",
                                      output = "measurement"),
    S = function() fa_bifactor(S = stats::cov(X), N = 400, factors = groups,
                               equal_error = eq, output = "measurement"))
  for (k in names(calls)) {
    seen <- new.env(parent = emptyenv())
    seen$msgs <- character()
    x <- withCallingHandlers(calls[[k]](), warning = function(w) {
      seen$msgs <- c(seen$msgs, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
    expect_true(any(startsWith(seen$msgs,
                               "The information matrix of this fit is singular")),
                label = k)
    # Through the first version of fa_bifactor() the maximal reliability
    # rows kept a standard error of 0 under FIML, and the listwise and
    # covariance matrix fits gave every row a finite interval.
    expect_true(all(is.na(x$se)), label = k)
    expect_true(all(is.na(x$test_statistic) & is.na(x$p_value)), label = k)
    expect_true(all(is.na(x$ci_lower) & is.na(x$ci_upper)), label = k)
    expect_false(anyNA(x$estimate), label = k)
    expect_equal(est(x, "maximal_reliability_general"), 0.7007, tolerance = 1e-3,
                 label = k)
  }
  # The standardized solution of the listwise fit follows the same rule.
  std <- suppressWarnings(fa_bifactor(X, groups, equal_error = eq,
                                      missing = "listwise",
                                      output = "standardized"))
  expect_true(all(is.na(std$se)))
  # When lavaan itself cannot invert the expected information (covariances
  # all equal, two groups of four), the existing branch gives the same NA
  # rows, the maximal reliability's included.
  Sigma <- matrix(0.3, 8L, 8L)
  diag(Sigma) <- 1
  dimnames(Sigma) <- list(paste0("x", 1:8), paste0("x", 1:8))
  four <- list(a = paste0("x", 1:4), b = paste0("x", 5:8))
  seen <- new.env(parent = emptyenv())
  seen$msgs <- character()
  y <- withCallingHandlers(
    fa_bifactor(S = Sigma, N = 300, factors = four, output = "measurement"),
    warning = function(w) {
      seen$msgs <- c(seen$msgs, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  expect_true(any(startsWith(seen$msgs,
                             "lavaan could not compute the standard errors")))
  expect_true(all(is.na(y$se)) && all(is.na(y$ci_lower)))
  expect_false(anyNA(y$estimate[startsWith(y$term, "maximal_reliability_")]))
})

test_that("a finite wald_df takes the t reference in the standardized output, the indices on the logit scale", {
  skip_on_cran()
  std <- fa_bifactor(holzinger_swineford, hs_groups, missing = "listwise",
                     wald_df = 20, output = "standardized")
  expect_identical(attr(std, "wald_df"), 20)
  q <- stats::qt(0.975, 20)
  ld <- std[std$op == "=~" & std$lhs == "general", ][1L, ]
  expect_equal(c(ld$ci.lower, ld$ci.upper), ld$est.std + c(-q, q) * ld$se,
               tolerance = 1e-8)
  expect_equal(ld$pvalue, 2 * stats::pt(-abs(ld$z), 20), tolerance = 1e-10)
  r <- std[std$lhs == "omega_H_general", ]
  d <- q * r$se / (r$est.std * (1 - r$est.std))
  expect_equal(c(r$ci.lower, r$ci.upper),
               stats::plogis(stats::qlogis(r$est.std) + c(-d, d)),
               tolerance = 1e-8)
})

test_that("the lavaan fit of output = \"fit\" feeds factor_scores(), whose regression determinacy is the square root of the maximal reliability", {
  skip_on_cran()
  fit <- fa_bifactor(holzinger_swineford, hs_groups, output = "fit")
  x <- fa_bifactor(holzinger_swineford, hs_groups, output = "measurement")
  fs <- factor_scores(fit)
  expect_s3_class(fs, "dmar_factor_scores")
  latent <- c("general", names(hs_groups))
  expect_identical(nrow(fs), 301L * length(latent))
  det <- attr(fs, "determinacy")
  expect_equal(det$determinacy_regression[match(latent, det$factor)]^2,
               est(x, paste0("maximal_reliability_", latent)),
               tolerance = 1e-6)
})

test_that("equality constraints are written with shared labels", {
  skip_on_cran()
  x <- fa_bifactor(holzinger_swineford, hs_groups, missing = "listwise",
                   equal_loading = c(general = FALSE, spatial = TRUE,
                                     verbal = FALSE, memory = FALSE),
                   equal_error = c(spatial = FALSE, verbal = FALSE,
                                   memory = TRUE))
  model <- paste(c(
    paste0("general =~ ", paste(unlist(hs_groups), collapse = " + ")),
    paste0("spatial =~ a*", paste(hs_groups$spatial, collapse = " + a*")),
    paste0("verbal =~ ", paste(hs_groups$verbal, collapse = " + ")),
    paste0("memory =~ ", paste(hs_groups$memory, collapse = " + ")),
    paste0(hs_groups$memory, " ~~ e*", hs_groups$memory)), collapse = "\n")
  ref <- lavaan::cfa(model, data = holzinger_swineford, std.lv = TRUE,
                     orthogonal = TRUE)
  expect_equal(est(x, "chi_square"),
               unname(lavaan::fitMeasures(ref, "chisq")), tolerance = 1e-5)
  expect_equal(est(x, "df"), unname(lavaan::fitMeasures(ref, "df")))
  expect_true(all(c("lambda_spatial", "psi_memory") %in% x$term))
  expect_identical(attr(x, "equal")$loading, "spatial")
  expect_identical(attr(x, "equal")$error, "memory")
  h <- bf_hand(ref, hs_groups)
  expect_equal(est(x, c("omega_spatial", "H_spatial", "ECV_memory")),
               unname(h[c("omega_spatial", "H_spatial", "ECV_memory")]),
               tolerance = 1e-5)
  # A two-item group factor is identified with equal loadings on it.
  two <- hs_groups
  two$spatial <- c("t1_visual_perception", "t4_lozenges")
  expect_s3_class(fa_bifactor(holzinger_swineford, two,
                              equal_loading = c(general = FALSE,
                                                spatial = TRUE,
                                                verbal = FALSE,
                                                memory = FALSE),
                              output = "measurement"),
                  "dmar_fa_bifactor")
})

test_that("structures outside the model, not identified, or with names that clash are refused by class, with the reason", {
  hs <- holzinger_swineford
  refused <- function(expr, pattern, reason) {
    e <- tryCatch(expr, error = function(e) e)
    expect_s3_class(e, "dmar_bifactor_unidentified")
    expect_s3_class(e, "dmar_refusal")
    expect_identical(e$code, "bifactor_unidentified")
    expect_identical(e$reason, reason)
    expect_match(conditionMessage(e), pattern)
    invisible(e)
  }
  # One group factor is a choice of scope, not a question of
  # identification: with items on the general factor alone it is
  # identified.
  e <- refused(fa_bifactor(hs, hs_groups["spatial"]),
               "two or more group factors.*names one \\('spatial'\\)",
               "scope")
  expect_match(conditionMessage(e), "rather than a question of identification")
  expect_false(grepl("not identified|unidentified", conditionMessage(e)))
  refused(fa_bifactor(hs, hs_groups["spatial"],
                      general_only = c("t20_deduction", "t21_numerical_puzzles")),
          "names one \\('spatial'\\)", "scope")
  one <- hs_groups
  one$spatial <- "t1_visual_perception"
  e <- refused(fa_bifactor(hs, one), "'spatial' has one item", "unidentified")
  expect_identical(e$factors, "spatial")
  two <- hs_groups
  two$spatial <- c("t1_visual_perception", "t4_lozenges")
  refused(fa_bifactor(hs, two), "'spatial' has two items with free loadings",
          "unidentified")
  both <- hs_groups
  both$verbal <- c(both$verbal, "t2_cubes")
  e <- refused(fa_bifactor(hs, both),
               "t2_cubes is listed under spatial and verbal", "scope")
  expect_identical(e$items, "t2_cubes")
  # An item listed twice under one group factor is a slip in the list, and
  # the message says so rather than calling it an item on two factors.
  twice <- hs_groups
  twice$spatial <- c(twice$spatial, "t2_cubes")
  e <- refused(fa_bifactor(hs, twice),
               "t2_cubes is listed more than once under spatial", "names")
  expect_false(grepl("two group factors", conditionMessage(e)))
  expect_identical(e$items, "t2_cubes")
  refused(fa_bifactor(hs, hs_groups, general = "verbal"),
          "general factor and a group factor are both named 'verbal'",
          "names")
  refused(fa_bifactor(hs, hs_groups, general_only = "t2_cubes"),
          "t2_cubes is named both under a group factor", "names")
  named_item <- hs_groups
  names(named_item)[1L] <- "t20_deduction"
  refused(fa_bifactor(hs, named_item, general_only = "t20_deduction"),
          "A factor and an item are both named 't20_deduction'", "names")
  # Labels made from the names coincide: a general factor 'g' numbers its
  # first loading lambda_g_1, which a group factor 'g_1' with equal
  # loadings also writes.
  clash <- hs_groups
  names(clash)[1L] <- "g_1"
  refused(fa_bifactor(hs, clash, general = "g",
                      equal_loading = c(g = FALSE, g_1 = TRUE,
                                        verbal = FALSE, memory = FALSE)),
          "labels that coincide \\(lambda_g_1\\)", "names")
  # Two groups of three items with nothing on the general factor alone:
  # the general loadings are fixed only up to a scale, and lavaan cannot
  # invert the information matrix of this fit.
  three <- list(spatial = c("t1_visual_perception", "t2_cubes",
                            "t4_lozenges"),
                verbal = c("t6_paragraph_comprehension", "t7_sentence",
                           "t9_word_meaning"))
  e <- refused(fa_bifactor(hs, three), "not identified", "unidentified")
  expect_setequal(e$factors, c("spatial", "verbal"))
  # (and so are two groups of two items with equal loadings)
  pairs <- list(a = three$spatial[1:2], b = three$verbal[1:2])
  refused(fa_bifactor(hs, pairs,
                      equal_loading = c(general = FALSE, a = TRUE, b = TRUE)),
          "not identified", "unidentified")
  # The checks run before anything is fit, so they need no lavaan fit.
  expect_error(fa_bifactor(hs, hs_groups, general = c("a", "b")),
               "one syntactically valid R name")
})

test_that("two group factors fit with a classed warning about empirical underidentification", {
  skip_on_cran()
  hs <- holzinger_swineford
  three <- list(spatial = c("t1_visual_perception", "t2_cubes",
                            "t4_lozenges"),
                verbal = c("t6_paragraph_comprehension", "t7_sentence",
                           "t9_word_meaning"))
  four <- three
  four$spatial <- c(four$spatial, "t3_paper_form_board")
  seen <- new.env(parent = emptyenv())
  x <- withCallingHandlers(
    fa_bifactor(hs, four, output = "measurement"),
    dmar_bifactor_empirical = function(cnd) {
      seen$w <- cnd
      invokeRestart("muffleWarning")
    })
  expect_s3_class(seen$w, "dmar_bifactor_empirical")
  expect_match(conditionMessage(seen$w), "Green & Yang, 2018")
  ref <- lavaan::cfa(bf_syntax(four), data = hs, std.lv = TRUE,
                     orthogonal = TRUE, missing = "ml")
  expect_equal(est(x, "omega_H_general"),
               unname(bf_hand(ref, four)["omega_H"]), tolerance = 1e-5)
  # Equal loadings within a group of three fix the scale.
  expect_warning(
    fa_bifactor(hs, three, output = "measurement",
                equal_loading = c(general = FALSE, spatial = TRUE,
                                  verbal = FALSE)),
    class = "dmar_bifactor_empirical")
  # The identification checks run after the arguments and the input, so a
  # call that then fails gets its error and no warning about a model it
  # does not fit.
  quiet_error <- function(expr, pattern) {
    seen <- new.env(parent = emptyenv())
    seen$w <- list()
    e <- tryCatch(withCallingHandlers(expr, warning = function(w) {
      seen$w <- c(seen$w, list(w))
      invokeRestart("muffleWarning")
    }), error = function(e) e)
    expect_s3_class(e, "error")
    expect_match(conditionMessage(e), pattern)
    expect_length(seen$w, 0L)
    invisible(e)
  }
  quiet_error(fa_bifactor(hs, four, estimator = "WLSMV"),
              "'estimator' must be one of")
  # A misnamed item, on a group factor or on the general factor alone, is
  # refused by the engine of fa_confirmatory() with the class of the factor
  # analysis items.
  items_refusal <- function(e, named) {
    expect_s3_class(e, c("dmar_items_invalid", "dmar_refusal", "error"))
    expect_identical(e$reason, "unknown")
    expect_identical(e$items, named)
  }
  misspelled <- four
  misspelled$verbal[1L] <- "t6_paragraph"
  items_refusal(quiet_error(fa_bifactor(hs, misspelled),
                            "Items not found in the data"), "t6_paragraph")
  items_refusal(quiet_error(fa_bifactor(hs, four, general_only = "t20_deductions"),
                            "Items not found in the data"), "t20_deductions")
  S <- cov(hs[c(unlist(four), "t20_deduction")], use = "complete.obs")
  items_refusal(quiet_error(fa_bifactor(S = S, factors = misspelled, N = nrow(hs)),
                            "Items not found in the covariance matrix"),
                "t6_paragraph")
  quiet_error(fa_bifactor(hs, four, conf_level = 2), "'conf_level' must be")
})

test_that("an improper solution leaves the maximal reliabilities NA, with warnings", {
  skip_on_cran()
  nine <- list(spatial = c("t1_visual_perception", "t2_cubes", "t4_lozenges"),
               verbal = c("t6_paragraph_comprehension", "t7_sentence",
                          "t9_word_meaning"),
               speed = c("t10_addition", "t12_counting_groups_of_dots",
                         "t13_straight_and_curved_capitals"))
  seen <- new.env(parent = emptyenv())
  seen$msgs <- character()
  x <- withCallingHandlers(
    fa_bifactor(holzinger_swineford, nine),
    warning = function(w) {
      seen$msgs <- c(seen$msgs, paste(class(w)[1L], conditionMessage(w)))
      invokeRestart("muffleWarning")
    })
  msgs <- seen$msgs
  expect_true(any(startsWith(msgs, "dmar_heywood_warning")))
  expect_true(any(grepl("Improper solution: the error variance of t1_visual_perception",
                        msgs)))
  imp <- attr(x, "improper")
  expect_identical(names(imp), "t1_visual_perception")
  expect_lt(imp[[1L]], 0)
  m <- x[startsWith(x$term, "maximal_reliability_"), ]
  expect_true(all(is.na(m$estimate)) && all(is.na(m$ci_lower)))
  # The other rows are reported as estimated.
  expect_false(is.na(est(x, "omega_H_general")))
  expect_true(any(grepl("Improper solution", .dmar_tbl_footer(x))))
})

test_that("what is not available is refused with the reason", {
  hs <- holzinger_swineford
  expect_error(fa_bifactor(hs, hs_groups, ordered = TRUE),
               "ordered items are not available.*bifactor_indices\\(\\)")
  expect_error(fa_bifactor(hs, hs_groups, ord = TRUE),
               "ordered items are not available")
  expect_error(fa_bifactor(hs, hs_groups, estimator = "WLSMV"),
               "must be one of .ML., .MLR.\\. fa_bifactor\\(\\) fits")
  expect_error(fa_bifactor(hs, hs_groups, missing = "two_stage"),
               "'missing' must be one of")
  expect_error(fa_bifactor(hs, hs_groups, wald_df = "missing_information"),
               "not available for fa_bifactor\\(\\): the calibration study")
  expect_error(fa_bifactor(hs, hs_groups, group = "school"),
               "'group' is not available in fa_bifactor\\(\\)")
  expect_error(fa_bifactor(hs, hs_groups, equal_error = c(general = TRUE)),
               "named with exactly the group factor names")
  expect_error(fa_bifactor(hs, hs_groups, equal_loading = c(spatial = TRUE)),
               "names of the general factor and the group factors")
})

test_that("the table's class, header, tidy(), glance(), and metadata", {
  skip_on_cran()
  x <- fa_bifactor(holzinger_swineford, hs_groups)
  expect_s3_class(x, c("dmar_fa_bifactor", "dmar_tbl", "data.frame"),
                  exact = TRUE)
  expect_identical(names(x), c("syntax", "term", "estimate", "se", "test_statistic",
                               "p_value", "ci_lower", "ci_upper"))
  # The rows: the parameters, then the indices in order, then the fit rows;
  # the factor covariances, fixed at 0 by the model, are not rows.
  idx <- c("omega_general", "omega_H_general", "ECV_general", "H_general",
           "maximal_reliability_general", "PUC_general",
           "omega_spatial", "omega_HS_spatial", "ECV_spatial", "H_spatial",
           "maximal_reliability_spatial")
  expect_identical(x$term[match("omega_general", x$term) + 0:10], idx)
  expect_identical(utils::tail(x$term, 1L), "H1")
  expect_false(any(x$syntax %in% c("general ~~ spatial", "general ~~ memory",
                                   "spatial ~~ verbal", "verbal ~~ memory")))
  expect_identical(est(x, "phi_general"), 1)
  lines <- .dmar_tbl_footer(x)
  expect_match(lines[names(lines) == "header"][1L],
               "the general factor \\(general\\) loads on all 15 items")
  expect_true(any(grepl("unit-weighted sums of the items as scored", lines)))
  td <- generics::tidy(x)
  expect_false("chi_square" %in% td$term)
  expect_true("maximal_reliability_general" %in% td$term)
  gl <- generics::glance(x)
  expect_identical(nrow(gl), 1L)
  expect_equal(gl$chi_square, est(x, "chi_square"))
  expect_identical(gl$missing, "fiml")
  md <- dmar_metadata(x)
  expect_identical(md$result, "dmar_fa_bifactor")
  expect_identical(md$estimator, "ML")
  expect_identical(attr(x, "omega_metric"), "raw_item_sum")
  sm <- fa_bifactor(holzinger_swineford, hs_groups, output = "summary")
  expect_s3_class(sm, "lavaan.summary")
})

test_that("the engine builds only the structures it knows", {
  # The engine's structure argument names only the structures it builds.
  expect_error(.fa_cfa_engine(data = holzinger_swineford,
                              factors = hs_groups, structure = "other"),
               "Unknown structure 'other'")
})

test_that("the RMSEA interval is at rmsea_conf_level, apart from the intervals of the parameters and indices", {
  skip_if_not_installed("lavaan")
  skip_on_cran()  # three bifactor fits
  # Every RMSEA interval defaults to 90 percent, whose limits decide the
  # tests of close and not close fit at .05; the parameters and indices
  # keep conf_level, 0.95 by default, and neither level reaches the other.
  rows <- c("rmsea_ci_lower", "rmsea_ci_upper", "rmsea_ci_level")
  get <- function(x, t) x$estimate[match(t, x$term)]
  x <- fa_bifactor(holzinger_swineford, hs_groups)
  y <- fa_bifactor(holzinger_swineford, hs_groups, rmsea_conf_level = 0.95)
  z <- fa_bifactor(holzinger_swineford, hs_groups, conf_level = 0.99)
  fit <- fa_bifactor(holzinger_swineford, hs_groups, output = "fit")
  fm <- function(level) {
    unname(lavaan::fitMeasures(fit, c("rmsea.ci.lower", "rmsea.ci.upper"),
                               fm.args = list(rmsea.ci.level = level)))
  }
  expect_equal(eval(formals(fa_bifactor)$rmsea_conf_level), 0.90)
  expect_equal(get(x, rows), c(fm(0.90), 0.90), tolerance = 1e-10)
  expect_equal(get(y, rows), c(fm(0.95), 0.95), tolerance = 1e-10)
  expect_identical(get(z, rows), get(x, rows))
  par <- !x$term %in% c(rows, "rmsea")
  expect_identical(y$ci_lower[par], x$ci_lower[par])
  expect_false(isTRUE(all.equal(z$ci_lower[par], x$ci_lower[par])))
  expect_true("Confidence level: 95% (90% for the RMSEA interval)" %in%
                utils::capture.output(print(x)))
  expect_identical(generics::glance(y)$rmsea_ci_level, 0.95)
  expect_identical(dmar_metadata(x)$rmsea_conf_level, 0.90)
})

# The standardized indices under equality constraints ---------------------
#
# The bifactor indices of the standardized items: a, b, and e are the
# standardized general loadings, group loadings (0 for an item on the
# general factor alone), and error variances, named by item.
bf_std_idx <- function(a, b, e, groups, general = "general") {
  cs <- vapply(names(groups), function(s) sum(b[groups[[s]]])^2, numeric(1L))
  tot <- sum(a)^2 + sum(cs) + sum(e)
  q <- sum(a^2 / (b^2 + e))
  out <- c(omega = (sum(a)^2 + sum(cs)) / tot, omega_H = sum(a)^2 / tot,
           ECV = sum(a^2) / sum(a^2 + b^2), H = q / (1 + q))
  names(out) <- paste0(names(out), "_", general)
  for (s in names(groups)) {
    r <- groups[[s]]
    sub <- sum(a[r])^2 + cs[[s]] + sum(e[r])
    q <- sum(b[r]^2 / (a[r]^2 + e[r]))
    out[paste0(c("omega_", "omega_HS_", "ECV_", "H_"), s)] <- c(
      (sum(a[r])^2 + cs[[s]]) / sub, cs[[s]] / sub,
      sum(a[r]^2) / sum(a[r]^2 + b[r]^2), q / (1 + q))
  }
  out
}

# The same indices from the items' own rows of a standardized table.
bf_std_direct <- function(z, groups, items, general = "general") {
  pick <- function(lhs) {
    r <- z$op == "=~" & z$lhs == lhs
    stats::setNames(z$est.std[r], z$rhs[r])
  }
  b <- stats::setNames(numeric(length(items)), items)
  for (s in names(groups)) {
    v <- pick(s)
    b[names(v)] <- v
  }
  r <- z$op == "~~" & z$lhs == z$rhs & z$lhs %in% items
  bf_std_idx(pick(general)[items], b,
             stats::setNames(z$est.std[r], z$lhs[r])[items], groups, general)
}

# An oracle on another code path: the indices as functions of the free
# parameter vector by matrix algebra on the parameter table (Sigma =
# Lambda Lambda' + Theta, each item standardized by its model implied
# variance), with the delta method standard error from a central
# difference Jacobian and lavaan's vcov().
bf_std_oracle <- function(fit, groups, general = "general") {
  pt <- lavaan::parTable(fit)
  x0 <- lavaan::coef(fit)
  rows <- which(pt$free > 0)
  pos <- if (length(x0) == length(rows)) seq_along(rows) else pt$free[rows]
  items <- lavaan::lavNames(fit, "ov")
  latent <- c(general, names(groups))
  at <- function(x) {
    est <- pt$est
    est[rows] <- x[pos]
    L <- matrix(0, length(items), length(latent),
                dimnames = list(items, latent))
    ld <- pt$op == "=~"
    L[cbind(pt$rhs[ld], pt$lhs[ld])] <- est[ld]
    th <- pt$op == "~~" & pt$lhs == pt$rhs & pt$lhs %in% items
    psi <- stats::setNames(est[th], pt$lhs[th])[items]
    d <- 1 / sqrt(rowSums(L^2) + psi)
    Ls <- L * d
    bf_std_idx(Ls[, general], rowSums(Ls[, names(groups), drop = FALSE]),
               psi * d^2, groups, general)
  }
  g0 <- at(x0)
  J <- vapply(seq_along(x0), function(k) {
    h <- 1e-5 * max(1, abs(x0[k]))
    xp <- x0
    xm <- x0
    xp[k] <- xp[k] + h
    xm[k] <- xm[k] - h
    (at(xp) - at(xm)) / (2 * h)
  }, numeric(length(g0)))
  list(est = g0, se = sqrt(diag(J %*% lavaan::vcov(fit) %*% t(J))))
}

test_that("standardized bifactor indices use each item's own standardized values when items share a label", {
  # Equal loadings on the spatial factor and equal error variances on the
  # memory factor, whose tests' variances differ. lavaan evaluates a
  # defined parameter at the first item's standardized value of a shared
  # label: through 1.1.0 development this fit gave omega_general .9113,
  # omega_HS_spatial .1222, H_spatial .1972, and omega_memory .8995, where
  # the standardized items give .8731, .2195, .4025, and .7233, and
  # lavaan's warning that the values were ambiguous reached the caller.
  eq_l <- c(general = FALSE, spatial = TRUE, verbal = FALSE, memory = FALSE)
  eq_e <- c(spatial = FALSE, verbal = FALSE, memory = TRUE)
  items <- unlist(hs_groups, use.names = FALSE)
  expect_no_warning(z <- fa_bifactor(holzinger_swineford, hs_groups,
                                     equal_loading = eq_l,
                                     equal_error = eq_e,
                                     output = "standardized"))
  fit <- fa_bifactor(holzinger_swineford, hs_groups, equal_loading = eq_l,
                     equal_error = eq_e, output = "fit")
  def <- z[z$op == ":=", ]
  got <- stats::setNames(def$est.std, def$lhs)
  se <- stats::setNames(def$se, def$lhs)
  direct <- bf_std_direct(z, hs_groups, items)
  expect_equal(got[names(direct)], direct, tolerance = 1e-10)
  expect_equal(unname(got[c("omega_general", "omega_HS_spatial", "H_spatial",
                            "omega_memory")]),
               c(0.873133, 0.219503, 0.402452, 0.723332), tolerance = 1e-5)
  o <- bf_std_oracle(fit, hs_groups)
  expect_equal(got[names(o$est)], o$est, tolerance = 1e-10)
  expect_equal(se[names(o$se)], o$se, tolerance = 1e-6)
  # The verbal rows use no shared label and are lavaan's own.
  lav <- suppressWarnings(lavaan::standardizedSolution(fit))
  v_rows <- paste0(c("omega_", "omega_HS_", "ECV_", "H_"), "verbal")
  for (col in c("est.std", "se", "z", "pvalue")) {
    expect_identical(def[[col]][match(v_rows, def$lhs)],
                     lav[[col]][match(v_rows, lav$lhs)])
  }
  # The z, p-value, and logit-scale limits of a recomputed row.
  r <- def[def$lhs == "omega_HS_spatial", ]
  expect_equal(r$z, r$est.std / r$se)
  expect_equal(r$pvalue, 2 * (1 - stats::pnorm(abs(r$z))))
  d <- stats::qnorm(0.975) * r$se / (r$est.std * (1 - r$est.std))
  expect_equal(c(r$ci.lower, r$ci.upper),
               stats::plogis(stats::qlogis(r$est.std) + c(-d, d)),
               tolerance = 1e-12)
  # ECV and H do not depend on the metric of the items, so they are the
  # rows of the unstandardized table, which lavaan reads unambiguously.
  m <- fa_bifactor(holzinger_swineford, hs_groups, equal_loading = eq_l,
                   equal_error = eq_e, output = "measurement")
  latent <- c("general", names(hs_groups))
  for (t in c(paste0("ECV_", latent), paste0("H_", latent))) {
    expect_equal(unname(got[t]), est(m, t), tolerance = 1e-10, label = t)
    expect_equal(unname(se[t]), m$se[m$term == t], tolerance = 1e-6,
                 label = t)
  }
})

test_that("equal general loadings with items on the general factor alone, listwise deletion, and se = \"none\" follow the same rule", {
  skip_on_cran()  # five more fits
  eq_l <- c(general = TRUE, spatial = FALSE, verbal = FALSE, memory = FALSE)
  eq_e <- c(spatial = FALSE, verbal = FALSE, memory = TRUE)
  alone <- c("t20_deduction", "t21_numerical_puzzles")
  items <- c(unlist(hs_groups, use.names = FALSE), alone)
  for (m in c("fiml", "listwise")) {
    expect_no_warning(z <- fa_bifactor(holzinger_swineford, hs_groups,
                                       general_only = alone,
                                       equal_loading = eq_l,
                                       equal_error = eq_e, missing = m,
                                       output = "standardized"))
    fit <- fa_bifactor(holzinger_swineford, hs_groups, general_only = alone,
                       equal_loading = eq_l, equal_error = eq_e,
                       missing = m, output = "fit")
    def <- z[z$op == ":=", ]
    o <- bf_std_oracle(fit, hs_groups)
    expect_equal(stats::setNames(def$est.std, def$lhs)[names(o$est)], o$est,
                 tolerance = 1e-10, label = m)
    expect_equal(stats::setNames(def$se, def$lhs)[names(o$se)], o$se,
                 tolerance = 1e-6, label = m)
    expect_equal(stats::setNames(def$est.std, def$lhs)[names(o$est)],
                 bf_std_direct(z, hs_groups, items)[names(o$est)],
                 tolerance = 1e-10, label = m)
  }
  # Without standard errors the estimates are the same, and the standard
  # error, test, and interval columns are NA.
  expect_no_warning(n <- fa_bifactor(holzinger_swineford, hs_groups,
                                     general_only = alone,
                                     equal_loading = eq_l,
                                     equal_error = eq_e, missing = "listwise",
                                     se = "none", output = "standardized"))
  expect_equal(n$est.std[n$op == ":="], z$est.std[z$op == ":="],
               tolerance = 1e-10)
  expect_true(all(is.na(unlist(n[c("se", "z", "pvalue", "ci.lower",
                                   "ci.upper")]))))
})
