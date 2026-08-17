# =============================================================================
#  Internal helpers for the reliability family
# =============================================================================
#
# This file contains the small, focused helper functions that the exported
# reliability functions
#
#     reliability()            -- general wrapper
#     reliability_alpha()      -- coefficient alpha (Guttman 1945; Cronbach 1951)
#     reliability_kr20()       -- Kuder-Richardson formula 20 (1937)
#     reliability_omega()      -- McDonald's coefficient omega (1999), with
#                                 a model implied or observed total-variance
#                                 denominator (the latter is the coefficient
#                                 Kelley & Pornprasertmanit 2016, Eq. 15,
#                                 call hierarchical omega)
#     reliability_omega_categorical()    -- categorical omega (Green & Yang 2009;
#                                 Kelley & Pornprasertmanit 2016, Eq. 17-18)
#
# share. Nothing in this file is exported. Names are dot-prefixed by
# convention to mark them as internal. There are no help pages for any of
# these symbols; a knowledgeable user can reach them via DMAR:::.name
# for inspection, but they are not part of the package's public API and
# their signatures may change.
#
# All equation numbers below refer to Kelley, K., & Pornprasertmanit, S.
# (2016). Confidence intervals for population reliability coefficients:
# Evaluation of methods, recommendations, and software for composite
# measures. Psychological Methods, 21(1), 69-92.
#
# Layout:
#   (A) Point estimators                .alpha_from_S, .kr20_from_data,
#                                       .omega_fit_cfa, .omega_observed_from_data,
#                                       .omega_c_from_data
#       Categorical-omega support       .bvn_cdf, .ordered_thresholds,
#                                       .polychoric_from_fit
#       Missing data (FIML) support     .aux_saturated_syntax, .fiml_cov
#   (B) Closed-form CI helpers          .ci_feldt, .ci_fisher, .ci_bonett,
#                                       .ci_hakstian_whalen
#   (C) Delta method CI helpers         .ci_alpha_ml, .ci_alpha_adf,
#                                       .ci_omega_delta
#   (D) Logistic transformation         .logistic_ci
#   (E) Bootstrap engine                .bootstrap_ci
#   (F) Tidy result formatter           .relia_result
#   (G) Input-validation helpers        .relia_resolve_inputs
#
# Style: code is written for readability, not for brevity. Each helper
# has a short comment block before it explaining what it computes, the
# inputs and outputs, and the relevant equation in the paper. Trivial
# one-liners are left without further comment.
# =============================================================================


# -----------------------------------------------------------------------------
# (A) Point estimators
# -----------------------------------------------------------------------------

# .alpha_from_S(S)
#
# Sample coefficient alpha computed directly from a covariance matrix of
# items. The closed-form expression is
#
#     alpha_hat = (J / (J - 1)) * (1 - sum(diag(S)) / sum(S))
#
# where J is the number of items, sum(diag(S)) is the sum of item
# variances (the trace of S), and sum(S) is the sum of all entries in S,
# which equals the variance of the composite Y = X_1 + ... + X_J by the
# variance-sum identity.
#
# Inputs:
#   S : a (J x J) symmetric numeric covariance matrix among items.
#
# Returns:
#   A single numeric value, the sample coefficient alpha. The value can
#   in principle fall outside [0, 1] (for example when items are
#   negatively correlated); the wrapper functions clamp the CI bounds,
#   not the point estimate.
.alpha_from_S <- function(S) {
  J <- nrow(S)
  (J / (J - 1)) * (1 - sum(diag(S)) / sum(S))
}


# .kr20_from_data(data)
#
# Sample Kuder-Richardson formula 20 (Kuder & Richardson, 1937) for a
# composite of dichotomous (0/1) items. On binary data the item variance
# equals p_j * q_j with the 1/N (biased) convention, or N/(N-1) * p_j * q_j
# with Bessel correction; the classical KR-20 formula
#
#     KR_20 = (J / (J - 1)) * (1 - sum(p_j * q_j) / s_Y^2)
#
# is algebraically identical to coefficient alpha computed from the item
# covariance matrix, provided that the same bias convention is used in
# both numerator and denominator. R's cov() uses 1/(N-1), so we route
# KR-20 through .alpha_from_S(cov(data)) and the two functions return the
# same number on the same dichotomous data. (Mixing 1/N for p_j q_j with
# 1/(N-1) for the composite variance produces a value that differs from
# alpha by a factor of N/(N-1).)
#
# Inputs:
#   data : a numeric matrix or data frame of 0/1 item scores; rows are
#          respondents, columns are items. Rows containing NAs are
#          listwise-deleted before computation.
#
# Returns:
#   A single numeric value, the sample KR-20.
#
# Errors:
#   If data contains values other than 0, 1, or NA.
.kr20_from_data <- function(data) {
  data <- as.matrix(data)
  if (!all(data %in% c(0, 1, NA))) {
    stop("KR-20 requires items scored 0/1; non-binary values were found.",
         call. = FALSE)
  }
  data <- data[stats::complete.cases(data), , drop = FALSE]
  .alpha_from_S(stats::cov(data))
}


# .omega_fit_cfa(data = NULL, S = NULL, N = NULL,
#                equal_loadings = FALSE,
#                estimator = "ML", se = "standard", missing = "listwise",
#                aux = NULL)
#
# Adapter on top of cfa_1() that fits a single-factor congeneric CFA and
# returns the parameter estimates in the structured-list form the
# reliability family expects. This helper exists because the CFA table
# data-frame outputs are designed for end-user readability, whereas the
# bootstrap engine and other internal callers in
# reliability_internals.R need direct access to numeric vectors of
# loadings and error variances. The adapter calls cfa_1(output = "fit"),
# which returns the raw lavaan fit object, and then extracts everything
# the reliability family needs in one place.
#
# When `aux` names auxiliary variables (columns of `data` that are not
# items), cfa_1() cannot be used, because the saturated correlates block
# is not part of the model cfa_1()/cfa_k() build. The helper then builds
# the same single-factor measurement syntax, appends the saturated
# correlates block from .aux_saturated_syntax(), and calls lavaan
# directly with missing = "ml". Phase 2 of the missing-data work
# (cfa_1() and cfa_k() themselves) can lift this path into cfa_k()
# proper.
#
# Equivalent model syntax to the legacy hand-rolled lavaan call:
#
#     f1 =~ NA*y1 + lambda_1*y1 + lambda_2*y2 + ... + lambda_J*yJ
#     f1 ~~ 1*f1
#     yj ~~ psi_j*yj          (for j = 1..J)
#     loading_sum := lambda_1 + ... + lambda_J
#     error_sum   := psi_1    + ... + psi_J
#     omega       := (loading_sum^2) / ((loading_sum^2) + error_sum)
#
# When equal_loadings = TRUE, all loadings share the single label
# "lambda" (tau-equivalence), implemented via cfa_1's equal_loading
# argument.
#
# Inputs:
#   data           : raw items (N rows by J columns); supply this OR S+N.
#   S, N           : a J x J covariance matrix and a sample size.
#   equal_loadings : if TRUE the fit imposes tau-equivalence; default
#                    FALSE (congeneric model).
#   estimator      : forwarded to cfa_1(); one of "ML" (default), "MLR",
#                    "WLS", "WLSMV", "GLS".
#   se             : forwarded to cfa_1(); typically "standard",
#                    "robust.sem", or "none".
#   missing        : forwarded to cfa_1() and used by lavaan only when
#                    raw data are supplied; "listwise" (default), "ml",
#                    "fiml", or "pairwise".
#   aux            : optional character vector naming auxiliary columns
#                    of `data` (saturated correlates; requires raw data
#                    and missing = "ml"). The remaining columns are the
#                    items.
#
# Returns:
#   A list with elements:
#     converged : logical, whether the fit converged.
#     omega     : numeric, sample coefficient omega.
#     se_omega  : numeric, lavaan's delta method SE for the omega
#                 defined parameter (NA_real_ when se = "none" or the
#                 SE could not be computed).
#     loadings  : numeric vector of length J, the lambda estimates.
#     errors    : numeric vector of length J, the psi^2 estimates.
#     J         : integer, the number of items.
#     fit       : the lavaan fit object (or NULL when the fit failed).
.omega_fit_cfa <- function(data = NULL, S = NULL, N = NULL,
                           equal_loadings = FALSE,
                           estimator = "ML",
                           se = "standard",
                           missing = "listwise",
                           aux = NULL) {
  if (!requireNamespace("lavaan", quietly = TRUE)) {
    stop("Coefficient omega requires the 'lavaan' package. ",
         "Install with install.packages(\"lavaan\").", call. = FALSE)
  }

  if (!is.null(aux)) {
    if (is.null(data)) {
      stop("Auxiliary variables require raw 'data'.", call. = FALSE)
    }
    return(.omega_fit_cfa_aux(data = data, aux = aux,
                              equal_loadings = equal_loadings,
                              estimator = estimator, se = se))
  }

  # Resolve item names (used when extracting estimates) and figure out J.
  if (!is.null(data)) {
    data <- as.data.frame(data)
    J <- ncol(data)
    var_names <- colnames(data)
    if (is.null(var_names)) {
      var_names <- paste0("y", seq_len(J))
      colnames(data) <- var_names
    }
    cfa_input <- data
    cfa_N     <- NULL
  } else {
    if (is.null(S) || is.null(N)) {
      stop("Either 'data' or both 'S' and 'N' must be supplied.",
           call. = FALSE)
    }
    J <- ncol(S)
    var_names <- colnames(S)
    if (is.null(var_names)) {
      var_names <- paste0("y", seq_len(J))
      colnames(S) <- rownames(S) <- var_names
    }
    cfa_input <- S
    cfa_N     <- N
  }

  # cfa_1() (a wrapper over cfa_k() since 2026-07-30) keeps the caller's
  # variable names, so downstream extraction matches against them; the
  # unnamed-input case was named var_names above, matching what the
  # wrapper would auto-assign.
  item_labels <- var_names

  # cfa_1() warns on an improper (Heywood) solution rather than stopping.
  # Catch that warning here instead of letting it out: this helper is also
  # the inner loop of the reliability bootstrap, where one warning per
  # replication would bury the user. The fact travels back on the returned
  # list as `improper`, so the single-fit callers can warn once and the
  # bootstrap can count.
  improper <- FALSE
  fit <- try(
    withCallingHandlers(
      cfa_1(data = if (is.null(cfa_N)) cfa_input else NULL,
            S = if (is.null(cfa_N)) NULL else cfa_input, N = cfa_N,
            equal_loading = equal_loadings,
            equal_error   = FALSE,
            estimator = estimator,
            se        = se,
            missing   = missing,
            output    = "fit"),
      dmar_heywood_warning = function(w) {
        improper <<- TRUE
        invokeRestart("muffleWarning")
      }
    ),
    silent = TRUE
  )

  if (inherits(fit, "try-error") ||
      isFALSE(lavaan::lavInspect(fit, "converged"))) {
    return(list(converged = FALSE, improper = FALSE,
                omega = NA_real_, se_omega = NA_real_,
                loadings = rep(NA_real_, J),
                errors   = rep(NA_real_, J),
                J = J, fit = NULL))
  }

  pe <- lavaan::parameterEstimates(fit)
  om_row   <- pe[pe$lhs == "omega_f1" & pe$op == ":=", , drop = FALSE]
  ld_rows  <- pe[pe$op == "=~", "est"]
  err_rows <- pe[pe$op == "~~" & pe$lhs == pe$rhs &
                 pe$lhs %in% item_labels, "est"]
  se_omega <- if (!is.null(om_row$se) && nrow(om_row) >= 1L) {
    om_row$se[1]
  } else {
    NA_real_
  }

  list(converged = TRUE,
       improper = improper,
       omega = om_row$est[1],
       se_omega = se_omega,
       loadings = ld_rows,
       errors   = err_rows,
       J = J,
       fit = fit)
}


# .omega_fit_cfa_aux(data, aux, equal_loadings, estimator, se)
#
# The auxiliary-variables branch of .omega_fit_cfa(): a single-factor
# model over the item columns of `data` with the auxiliary columns
# entered as saturated correlates, estimated by full information
# maximum likelihood (lavaan missing = "ml"). The measurement syntax
# mirrors cfa_1() exactly (fixed unit factor variance, labeled loadings
# and error variances, the loading_sum / error_sum / omega defined
# parameters), so the estimates and the omega delta method standard
# error are the same quantities the cfa_1() route produces; only the
# saturated correlates block and the FIML mean structure are added.
#
# Inputs:
#   data           : data frame containing item AND auxiliary columns.
#   aux            : character vector naming the auxiliary columns; the
#                    remaining columns are the items, in column order.
#   equal_loadings : TRUE imposes tau-equivalence.
#   estimator      : "ML" or "MLR" (WLS cannot be combined with FIML;
#                    callers validate before reaching this point).
#   se             : lavaan se type ("standard", "robust.sem", "none").
#
# Returns:
#   The same structured list as .omega_fit_cfa().
.omega_fit_cfa_aux <- function(data, aux, equal_loadings = FALSE,
                               estimator = "ML", se = "standard") {
  data <- as.data.frame(data)
  items <- setdiff(colnames(data), aux)
  J <- length(items)
  K <- length(aux)

  # Rename to positional labels, as cfa_1() does, so the syntax never
  # trips over user column names that are not valid lavaan identifiers.
  y_names <- paste0("y", seq_len(J))
  a_names <- paste0("a", seq_len(K))
  d <- data[, c(items, aux), drop = FALSE]
  colnames(d) <- c(y_names, a_names)

  loading_labels <- if (equal_loadings) {
    rep("lambda", J)
  } else {
    paste0("lambda_", seq_len(J))
  }
  error_labels <- paste0("psi_", seq_len(J))

  measurement <- paste0(
    "f1 =~ NA*y1 + ",
    paste(paste0(loading_labels, "*", y_names), collapse = " + "), "\n",
    "f1 ~~ 1*f1\n",
    paste(paste0(y_names, " ~~ ", error_labels, "*", y_names),
          collapse = "\n"), "\n",
    "loading_sum := ", paste(loading_labels, collapse = " + "), "\n",
    "error_sum := ", paste(error_labels, collapse = " + "), "\n",
    "omega := (loading_sum^2) / ((loading_sum^2) + error_sum)\n"
  )
  model_syntax <- paste0(measurement,
                         .aux_saturated_syntax(y_names, a_names), "\n")

  # check.post = FALSE: lavaan's post-estimation check warns that the
  # residual covariance matrix (theta) is not positive definite whenever
  # an auxiliary correlates substantially with the factor, because the
  # saturated correlates parameterization forces the entire item-aux
  # association through the residual covariances (the aux never touches
  # the factor). That is the expected signature of the specification,
  # not a defect of the data or the fit; the joint implied covariance
  # matrix is still proper. The check that matters for reliability, a
  # negative item error variance (Heywood case), is performed explicitly
  # below and travels back on the `improper` flag.
  fit <- try(lavaan::lavaan(model_syntax, data = d,
                            missing = "ml",
                            estimator = estimator,
                            se = se,
                            meanstructure = TRUE,
                            int.ov.free = TRUE,
                            check.post = FALSE),
             silent = TRUE)

  if (inherits(fit, "try-error") ||
      isFALSE(lavaan::lavInspect(fit, "converged"))) {
    return(list(converged = FALSE, improper = FALSE,
                omega = NA_real_, se_omega = NA_real_,
                loadings = rep(NA_real_, J),
                errors   = rep(NA_real_, J),
                J = J, fit = NULL))
  }

  pe <- lavaan::parameterEstimates(fit)
  om_row   <- pe[pe$lhs == "omega" & pe$op == ":=", , drop = FALSE]
  ld_rows  <- pe[pe$op == "=~", "est"]
  err_rows <- pe[pe$op == "~~" & pe$lhs == pe$rhs &
                 pe$lhs %in% y_names, "est"]
  se_omega <- if (!is.null(om_row$se) && nrow(om_row) >= 1L) {
    om_row$se[1]
  } else {
    NA_real_
  }

  list(converged = TRUE,
       improper = any(err_rows < 0),
       omega = om_row$est[1],
       se_omega = se_omega,
       loadings = ld_rows,
       errors   = err_rows,
       J = J,
       fit = fit)
}


# .aux_saturated_syntax(items, aux)
#
# lavaan syntax for the saturated correlates treatment of auxiliary
# variables (Graham, 2003): each auxiliary gets a free variance, the
# auxiliaries covary freely with each other, and each auxiliary covaries
# with every item's residual. The auxiliaries are never indicators of
# the factor and never enter the composite, so the measurement model's
# parameters are unchanged in expectation while full information
# maximum likelihood uses the auxiliaries' information about the
# missing values.
#
# Inputs:
#   items : character vector of item names as they appear in the model
#           syntax (an item's `~~` with an auxiliary is a residual
#           covariance, because the item is an indicator of the factor).
#   aux   : character vector of auxiliary variable names.
#
# Returns:
#   A single string of newline-separated lavaan syntax lines.
.aux_saturated_syntax <- function(items, aux) {
  K <- length(aux)
  lines <- paste0(aux, " ~~ ", aux)
  if (K >= 2L) {
    for (k in seq_len(K - 1L)) {
      lines <- c(lines,
                 paste0(aux[k], " ~~ ",
                        paste(aux[(k + 1L):K], collapse = " + ")))
    }
  }
  for (a in aux) {
    lines <- c(lines, paste0(a, " ~~ ", paste(items, collapse = " + ")))
  }
  paste(lines, collapse = "\n")
}


# .fiml_cov(data, items, aux = NULL, se = "none")
#
# Full information maximum likelihood estimate of the covariance matrix
# among the items, for the analytic (closed form) reliability path under
# missing = "fiml". A saturated model over the items plus any auxiliary
# variables (every variance and covariance free, saturated mean
# structure) is fit with lavaan missing = "ml"; the item block of the
# implied covariance matrix is the FIML Sigma-hat the classical formula
# is then applied to. Including the auxiliaries in the saturated model
# is what lets their information flow into the item block.
#
# The model also defines coefficient alpha of the item block as a `:=`
# parameter,
#
#     alpha := (J/(J-1)) * (1 - trace(Sigma_items) / sum(Sigma_items)),
#
# so that when se = "standard" lavaan's delta method gives a standard
# error for the analytic alpha that comes from the FIML information
# matrix (and therefore reflects the actual missingness pattern, unlike
# a complete-data closed form evaluated at the used N).
#
# Inputs:
#   data  : data frame containing the item (and auxiliary) columns.
#   items : character vector of item column names.
#   aux   : optional character vector of auxiliary column names.
#   se    : "none" (default; point estimates only, the fast path the
#           bootstrap uses) or "standard" (delta method SE for alpha).
#
# Returns:
#   A list with elements:
#     converged : logical.
#     S         : (J x J) FIML covariance matrix of the items, with the
#                 original item names (NULL when the fit failed).
#     N         : number of cases lavaan used.
#     alpha     : the defined-parameter estimate of coefficient alpha.
#     alpha_se  : its delta method SE (NA_real_ when se = "none").
.fiml_cov <- function(data, items, aux = NULL, se = "none") {
  if (!requireNamespace("lavaan", quietly = TRUE)) {
    stop("missing = \"fiml\" requires the 'lavaan' package. ",
         "Install with install.packages(\"lavaan\").", call. = FALSE)
  }
  data <- as.data.frame(data)
  J <- length(items)
  K <- length(aux)
  y_names <- paste0("y", seq_len(J))
  # paste0() pads a zero-length argument to "" by default (recycle0 is
  # FALSE), so the no-auxiliary case must be handled explicitly.
  a_names <- if (K > 0L) paste0("a", seq_len(K)) else character(0)
  d <- data[, c(items, aux), drop = FALSE]
  colnames(d) <- c(y_names, a_names)
  vars <- c(y_names, a_names)
  p <- length(vars)

  # Saturated covariance structure with a label on every parameter; the
  # item-block labels feed the alpha definition.
  lab <- function(i, j) paste0("s", min(i, j), "_", max(i, j))
  lines <- vapply(seq_len(p), function(i) {
    paste0(vars[i], " ~~ ", lab(i, i), "*", vars[i])
  }, character(1))
  for (i in seq_len(p - 1L)) {
    for (j in (i + 1L):p) {
      lines <- c(lines, paste0(vars[i], " ~~ ", lab(i, j), "*", vars[j]))
    }
  }
  trace_terms <- vapply(seq_len(J), function(i) lab(i, i), character(1))
  cov_terms <- character(0)
  for (i in seq_len(J - 1L)) {
    for (j in (i + 1L):J) cov_terms <- c(cov_terms, lab(i, j))
  }
  trace_expr <- paste(trace_terms, collapse = " + ")
  total_expr <- paste0(trace_expr, " + 2*(",
                       paste(cov_terms, collapse = " + "), ")")
  alpha_def <- paste0("alpha := (", J, "/", J - 1L, ") * (1 - (",
                      trace_expr, ") / (", total_expr, "))")
  model_syntax <- paste(c(lines, alpha_def), collapse = "\n")

  fit <- try(lavaan::lavaan(model_syntax, data = d,
                            missing = "ml",
                            estimator = "ML",
                            se = se,
                            meanstructure = TRUE,
                            int.ov.free = TRUE),
             silent = TRUE)
  if (inherits(fit, "try-error") ||
      isFALSE(lavaan::lavInspect(fit, "converged"))) {
    return(list(converged = FALSE, S = NULL, N = NA_integer_,
                alpha = NA_real_, alpha_se = NA_real_))
  }

  Sigma <- lavaan::lavInspect(fit, "implied")$cov
  S_items <- as.matrix(Sigma[y_names, y_names, drop = FALSE])
  dimnames(S_items) <- list(items, items)
  pe <- lavaan::parameterEstimates(fit)
  al_row <- pe[pe$lhs == "alpha" & pe$op == ":=", , drop = FALSE]
  alpha_se <- if (se == "none" || is.null(al_row$se) ||
                  nrow(al_row) < 1L) {
    NA_real_
  } else {
    al_row$se[1]
  }

  list(converged = TRUE,
       S = S_items,
       N = lavaan::lavInspect(fit, "nobs"),
       alpha = al_row$est[1],
       alpha_se = alpha_se)
}


# .omega_observed_from_data(data, estimator = "mlr")
#
# Point estimator behind reliability_omega(denominator = "observed"):
# coefficient omega with the observed (rather than model implied) total
# variance in the denominator, the coefficient Kelley & Pornprasertmanit
# (2016, Eq. 15) call hierarchical omega. From raw item data:
#
#     omega_observed = (sum lambda_j)^2 / var(Y)
#
# The numerator is the squared sum of loadings from a single-factor CFA.
# The denominator is the observed variance of the composite on the same
# maximum likelihood (N-divisor) metric as the fitted loadings, i.e.
# sum(cov(data)) * (n - 1) / n (cov() uses the N - 1 divisor). This keeps
# the ratio consistent and matches MBESS::ci.reliability(type="hierarchical").
# The observed total variance absorbs minor-factor variability into the
# denominator instead of forcing it into the error term, and it is a
# consistent estimate of var(Y) whether or not the single-factor model is
# correctly specified. For perfectly fitting models the two denominators
# agree and the coefficient equals the model implied omega.
#
# Inputs:
#   data      : raw items (N rows by J columns).
#   estimator : lavaan estimator for the single-factor fit; defaults to
#               "mlr" (robust ML), per the paper's recommendation under
#               nonnormality.
#
# Returns:
#   A single numeric value. NA_real_ if the CFA fit did not converge.
.omega_observed_from_data <- function(data, estimator = "MLR") {
  fit <- .omega_fit_cfa(data = data, equal_loadings = FALSE,
                        estimator = estimator, se = "none",
                        missing = "listwise")
  if (!fit$converged) return(NA_real_)
  d_cc <- stats::na.omit(as.matrix(data))
  S <- stats::cov(d_cc)
  n <- nrow(d_cc)
  # Rescale the observed total variance to the fitted loadings' N-divisor
  # (ML) metric so the ratio matches MBESS::ci.reliability(type = "hierarchical").
  # Listwise deletion in the fit keeps numerator and denominator on the
  # same complete-case sample.
  (sum(fit$loadings))^2 / (sum(S) * (n - 1) / n)
}


# .omega_c_from_data(data)
#
# Categorical omega (omega_C; Green & Yang, 2009; Kelley &
# Pornprasertmanit, 2016, Eqs. 17-18) from raw ordered-categorical item
# data. Each observed item X_j is treated as a categorization of a
# latent continuous variable X_j* via thresholds t_{j,c}. Under the
# delta parameterization Var(X_j*) = 1, the model implied covariance of
# (X_j, X_{j'}) at latent correlation r is
#
#     sigma_{j j'}(r) = sum_c sum_{c'} Phi_2(t_{j,c}, t_{j',c'}; r)
#                       - (sum_c Phi_1(t_{j,c})) *
#                         (sum_{c'} Phi_1(t_{j',c'}))
#
# where Phi_1 is the univariate standard normal CDF and Phi_2 is the
# bivariate standard normal CDF at correlation r. Categorical omega is
#
#     omega_C = [ sum_{j, j'} sigma_{j j'}(lambda_j * lambda_{j'}) ]
#               --------------------------------------------------
#               [ sum_{j, j'} sigma_{j j'}(rho_{X_j* X_{j'}*})    ]
#
# with model implied (factor-derived) polychoric correlations in the
# numerator and freely-estimated polychoric correlations in the
# denominator. Treating ordered-categorical items as continuous and
# computing coefficient omega would mis-model the link between latent
# ability and observed responses, especially when categories are few or
# threshold patterns differ across items.
#
# Inputs:
#   data : a data frame or matrix of ordered-categorical item scores
#          (integer-valued). Columns are coerced to ordered factors
#          before fitting.
#
# Returns:
#   A single numeric value, sample omega_C. NA_real_ if the categorical
#   CFA fit did not converge.
#
# Dependencies:
#   lavaan (single-factor categorical CFA), mvtnorm (bivariate normal
#   CDF).
.omega_c_from_data <- function(data) {
  if (!requireNamespace("lavaan", quietly = TRUE)) {
    stop("Categorical omega requires the 'lavaan' package. ",
         "Install with install.packages(\"lavaan\").", call. = FALSE)
  }
  if (!requireNamespace("mvtnorm", quietly = TRUE)) {
    stop("Categorical omega requires the 'mvtnorm' package. ",
         "Install with install.packages(\"mvtnorm\").", call. = FALSE)
  }

  data <- as.data.frame(data)
  J <- ncol(data)
  for (j in seq_len(J)) data[[j]] <- ordered(data[[j]])
  var_names <- paste0("y", seq_len(J))
  colnames(data) <- var_names

  load_line <- paste(paste0("a", seq_len(J), "*", var_names),
                     collapse = " + ")
  model <- paste0("f1 =~ NA*", var_names[1], " + ", load_line,
                  "\nf1 ~~ 1*f1\n")
  fit <- try(lavaan::cfa(model, data = data, se = "none",
                         ordered = var_names),
             silent = TRUE)
  if (inherits(fit, "try-error") ||
      isFALSE(lavaan::lavInspect(fit, "converged"))) {
    return(NA_real_)
  }

  # Extract loadings, factor variance (= 1), thresholds, and a separate
  # estimate of the full polychoric correlation matrix.
  est <- lavaan::lavInspect(fit, "est")
  lambda <- est$lambda
  psi    <- est$psi
  truevar <- lambda %*% psi %*% t(lambda)
  thresh  <- .ordered_thresholds(fit)
  poly    <- .polychoric_from_fit(fit, data, var_names)
  sigma_hat <- lavaan::lavInspect(fit, "implied")$cov
  inv_sd  <- 1 / sqrt(diag(sigma_hat))
  poly_r  <- diag(inv_sd) %*% truevar %*% diag(inv_sd)

  .gy_omega_from_components(thresh, poly_r, poly)
}


# .gy_omega_from_components(thresh, rho_num, rho_den)
#
# The Green and Yang (2009) assembly shared by the single-factor
# categorical omega and the per-factor categorical omega of cfa_k():
# given the item thresholds and two latent-response correlation
# matrices, sum the bivariate normal cell probabilities over all
# threshold pairs to obtain the covariance of the categorical sum
# score attributable to the common factor(s) (numerator, rho_num =
# model implied correlations) and in total (denominator, rho_den =
# saturated polychoric correlations), and return their ratio.
#
# Inputs:
#   thresh  : list of numeric threshold vectors, one per item, in item
#             order.
#   rho_num : (J x J) model implied latent-response correlation matrix
#             (common-factor part, including communalities on the
#             diagonal).
#   rho_den : (J x J) saturated polychoric correlation matrix.
#
# Returns:
#   A single numeric value, the Green and Yang categorical omega.
.gy_omega_from_components <- function(thresh, rho_num, rho_den) {
  J <- length(thresh)
  num <- 0; den <- 0
  for (j in seq_len(J)) {
    for (jp in seq_len(J)) {
      t1 <- thresh[[j]]
      t2 <- thresh[[jp]]
      sum_p2_num <- 0
      sum_p2_den <- 0
      for (c in seq_along(t1)) {
        for (cp in seq_along(t2)) {
          sum_p2_num <- sum_p2_num + .bvn_cdf(t1[c], t2[cp], rho_num[j, jp])
          sum_p2_den <- sum_p2_den + .bvn_cdf(t1[c], t2[cp], rho_den[j, jp])
        }
      }
      p1  <- sum(stats::pnorm(t1))
      p1p <- sum(stats::pnorm(t2))
      num <- num + (sum_p2_num - p1 * p1p)
      den <- den + (sum_p2_den - p1 * p1p)
    }
  }
  num / den
}


# .omega_categorical_from_kfit(fit, items_f)
#
# Green and Yang (2009) categorical omega for the sum score of one
# factor's items, computed from a fitted multi-factor categorical
# model (cfa_k() with ordered items). The numerator correlations are
# the model implied latent-response correlations restricted to the
# factor's items (the common-factor part, Lambda Phi Lambda' scaled to
# the correlation metric); the denominator correlations are the
# saturated polychorics lavaan computed as the categorical sample
# statistics; the thresholds come from the same fit.
#
# Inputs:
#   fit     : a converged lavaan fit with ordered indicators.
#   items_f : character vector of the factor's item names.
#
# Returns:
#   A single numeric value, the factor's categorical sum score omega.
.omega_categorical_from_kfit <- function(fit, items_f) {
  est <- lavaan::lavInspect(fit, "est")
  lambda <- est$lambda
  phi    <- est$psi
  truevar <- lambda %*% phi %*% t(lambda)
  dimnames(truevar) <- list(rownames(lambda), rownames(lambda))
  # Scale the common part to the correlation metric of the latent
  # responses; under the theta parameterization the implied variances
  # are not 1, and this rescaling makes the computation
  # parameterization invariant.
  implied <- lavaan::lavInspect(fit, "implied")$cov
  inv_sd <- 1 / sqrt(diag(implied))
  rho_num <- diag(inv_sd) %*% truevar %*% diag(inv_sd)
  dimnames(rho_num) <- dimnames(truevar)
  samp <- lavaan::lavInspect(fit, "sampstat")
  rho_den <- samp$cov
  # Saturated thresholds on the standardized latent response metric,
  # matching the polychoric correlations in rho_den.
  th_names <- names(samp$th)
  bar_pos <- regexpr("\\|", th_names)
  th_var <- substr(th_names, 1, bar_pos - 1)
  thresh <- split(as.numeric(samp$th),
                  factor(th_var, levels = unique(th_var)))
  .gy_omega_from_components(thresh[items_f],
                            rho_num[items_f, items_f, drop = FALSE],
                            rho_den[items_f, items_f, drop = FALSE])
}


# .bvn_cdf(t1, t2, r)
#
# Bivariate standard-normal cumulative distribution function evaluated
# at the upper limits (t1, t2) with correlation r. Used in the inner
# loop of categorical omega.
#
# Inputs:
#   t1, t2 : numeric scalars (threshold values).
#   r      : numeric in [-1, 1] (correlation).
# Returns:
#   A single numeric value: P(Z_1 <= t1, Z_2 <= t2) where
#   (Z_1, Z_2) ~ N(0, [[1, r], [r, 1]]).
.bvn_cdf <- function(t1, t2, r) {
  Sigma <- matrix(c(1, r, r, 1), 2, 2)
  mvtnorm::pmvnorm(upper = c(t1, t2), mean = c(0, 0), corr = Sigma)[1]
}


# .ordered_thresholds(fit)
#
# Extract per-item threshold vectors from a fitted ordered-categorical
# lavaan CFA. lavaan stores thresholds in a tau matrix whose rownames
# have the form "y1|t1", "y1|t2", "y2|t1", ...; we split those rows into
# a list of vectors indexed by item.
#
# Inputs:
#   fit : a fitted lavaan ordered-categorical CFA.
# Returns:
#   A named list of numeric vectors, one per item, each containing the
#   (C - 1) thresholds for an item with C categories.
.ordered_thresholds <- function(fit) {
  coef_lst <- lavaan::lavInspect(fit, "est")
  tau <- coef_lst$tau
  names_tau <- rownames(tau)
  bar_pos <- regexpr("\\|", names_tau)
  var_for_tau <- substr(names_tau, 1, bar_pos - 1)
  split(as.numeric(tau),
        factor(var_for_tau, levels = unique(var_for_tau)))
}


# .polychoric_from_fit(fit, data, var_names)
#
# Estimate the polychoric correlation matrix of the items by refitting a
# saturated model in which every pair of items is correlated. The theta
# matrix returned by lavaan from this refit is the polychoric correlation
# matrix used in the denominator of categorical omega.
#
# Note on portability: the legacy MBESS implementation extracted the
# CFA function from the previous fit's @call slot and called it via
# do.call(); that fails silently when lavaan is not attached to the
# search path. We call lavaan::cfa() directly here for robustness.
#
# Inputs:
#   fit       : the previous lavaan fit (used only to inherit conventions).
#   data      : data frame of ordered-factor item scores.
#   var_names : character vector of item names.
# Returns:
#   A (J x J) numeric matrix of estimated polychoric correlations.
.polychoric_from_fit <- function(fit, data, var_names) {
  script <- ""
  for (i in 2:length(var_names)) {
    rhs <- paste0(var_names[seq_len(i - 1)], collapse = " + ")
    script <- paste0(script, var_names[i], " ~~ ", rhs, "\n")
  }
  refit <- suppressWarnings(
    lavaan::cfa(script, data = data, ordered = var_names, se = "none")
  )
  theta <- lavaan::lavInspect(refit, "est")$theta
  # The saturated script introduces the items in the order y2, y3, ..., yJ, y1,
  # so lavaan returns theta in that order of appearance. Reorder to var_names so
  # the caller, which indexes positionally, pairs each item's thresholds with
  # its own polychoric correlations rather than a permuted set.
  theta[var_names, var_names, drop = FALSE]
}


# -----------------------------------------------------------------------------
# (B) Closed-form CI helpers
#
# Each helper takes a point estimate, sample size N, number of items J
# (where the formula depends on it), and a confidence level. None of
# them needs a fitted model; they apply mechanically to any quantity on
# the (0, 1) reliability scale.
#
# All helpers return a list with named elements `se`, `lower`, `upper`,
# and `se` is always on the reliability (coefficient) scale. The
# transformation-based helpers (.ci_fisher, .ci_bonett,
# .ci_hakstian_whalen) compute their standard error on the
# transformation scale, so each also returns `se_transformed` (that
# standard error, the quantity the Wald interval is built from) and
# `se_transform_scale` (a label naming the scale); their `se` is the
# delta method back-transform
#
#     se = se_transformed * |d(back-transform)/dz|
#
# evaluated at the estimate (decided 2026-08-12, so the se row reads on
# the same scale as the estimate while simulation users keep the
# transformation-scale value). When a method does not produce an SE
# (.ci_feldt, an F pivot), `se` is NA_real_.
# -----------------------------------------------------------------------------

# .ci_feldt(estimate, N, J, conf_level)
#
# Feldt's (1965) F-distribution interval for the population reliability
# of a homogeneous parallel-items composite. The interval is
#
#     [ 1 - (1 - alpha_hat) * F_{1 - q, N - 1, (N - 1)(J - 1)} ,
#       1 - (1 - alpha_hat) * F_{q,     N - 1, (N - 1)(J - 1)} ]
#
# where q = (1 - C)/2 and C is the confidence level. The interval is
# exact under multivariate normality and parallel items. The method does
# not produce an SE.
.ci_feldt <- function(estimate, N, J, conf_level) {
  q   <- (1 - conf_level) / 2
  df1 <- N - 1
  df2 <- (N - 1) * (J - 1)
  f_upper <- stats::qf(1 - q, df1, df2)   # gives the lower CI limit
  f_lower <- stats::qf(    q, df1, df2)   # gives the upper CI limit
  list(se    = NA_real_,
       lower = 1 - (1 - estimate) * f_upper,
       upper = 1 - (1 - estimate) * f_lower)
}


# .ci_fisher(estimate, N, conf_level)
#
# Fisher's Z transformation of the reliability coefficient (Fisher,
# 1950). The coefficient is treated like a correlation and transformed
# via z = 0.5 * log((1 + r) / (1 - r)) with SE = 1 / sqrt(N - 3). The
# Wald interval is built on the z scale and back-transformed to the
# reliability scale. Kelley & Pornprasertmanit (2016) found this method
# to over-cover and did not recommend it; it is included for parity with
# the legacy MBESS interface and for users who want it. Requires N >= 4:
# the z variance 1/(N - 3) is infinite at N = 3 (the interval would be
# vacuous) and undefined below it. Returns both standard error scales:
# `se_transformed` is the z-scale value the interval is built from, and
# `se` is its delta method back-transform to the coefficient scale (the
# back-transform is r = tanh(z), so |dr/dz| = 1 - r^2 at the estimate).
.ci_fisher <- function(estimate, N, conf_level) {
  if (N < 4) {
    stop("The Fisher confidence interval requires 'N' of at least 4: the ",
         "variance of the transformed coefficient is 1/(N - 3), which is ",
         "infinite at N = 3, where the interval would be vacuous, and ",
         "undefined for smaller 'N'.", call. = FALSE)
  }
  crit <- stats::qnorm(1 - (1 - conf_level) / 2)
  z  <- 0.5 * log((1 + estimate) / (1 - estimate))
  se <- sqrt(1 / (N - 3))
  z_lower <- z - crit * se
  z_upper <- z + crit * se
  list(se    = se * (1 - estimate^2),
       se_transformed     = se,
       se_transform_scale = "fisher_z",
       lower = (exp(2 * z_lower) - 1) / (exp(2 * z_lower) + 1),
       upper = (exp(2 * z_upper) - 1) / (exp(2 * z_upper) + 1))
}


# .ci_bonett(estimate, N, J, conf_level)
#
# Bonett's (2002) log transformation z = log(1 - r). The SE on the
# transformed scale is sqrt(2J / ((J - 1) * (N - 2))); the Wald interval
# is built on the z scale and back-transformed via 1 - exp(z). The
# transformation accounts for the asymmetry of the sampling distribution
# of a reliability coefficient near 1. Performs well under normality in
# Kelley & Pornprasertmanit's (2016) simulations. Returns both standard
# error scales: `se_transformed` is the log-scale value the interval is
# built from, and `se` is its delta method back-transform to the
# coefficient scale (the back-transform is r = 1 - exp(z), so
# |dr/dz| = exp(z) = 1 - r at the estimate).
.ci_bonett <- function(estimate, N, J, conf_level) {
  crit <- stats::qnorm(1 - (1 - conf_level) / 2)
  z  <- log(1 - estimate)
  se <- sqrt(2 * J / ((J - 1) * (N - 2)))
  list(se    = se * (1 - estimate),
       se_transformed     = se,
       se_transform_scale = "log(1-alpha)",
       lower = 1 - exp(z + crit * se),   # large z' -> small reliability
       upper = 1 - exp(z - crit * se))
}


# .ci_hakstian_whalen(estimate, N, J, conf_level)
#
# Hakstian & Whalen's (1976) cube-root transformation z = (1 - r)^(1/3).
# The transformation reduces the skewness of the sampling distribution
# more than Fisher's or Bonett's transformations do. The standard error
# on the transformed scale (Kelley & Pornprasertmanit, 2016, Eq. 29) is
#
#     SE = sqrt( (2 J / (9 * df_d)) * (1 - r)^(2/3)
#                / (1 - 2 / (9 * df_n))^2 )
#
# with df_n = N - 1 and df_d = (N - 1)(J - 1). The back-transformation
# uses an `a` factor (Kelley & Pornprasertmanit, 2016, Eq. 30) that
# corrects for the small-N bias of the cube-root mean. Returns both
# standard error scales: `se_transformed` is the cube-root-scale value
# the interval is built from, and `se` is its delta method
# back-transform to the coefficient scale (the back-transform is
# r = 1 - a^3 z^3, so |dr/dz| = 3 a^3 z^2 = 3 a^3 (1 - r)^(2/3) at the
# estimate).
.ci_hakstian_whalen <- function(estimate, N, J, conf_level) {
  crit <- stats::qnorm(1 - (1 - conf_level) / 2)
  df_n <- N - 1
  df_d <- (N - 1) * (J - 1)
  a_factor <- (1 - 2 / (9 * df_n)) / (1 - 2 / (9 * df_d))
  numerator   <- (2 * J / (9 * df_d)) * (1 - estimate)^(2 / 3)
  denominator <- (1 - 2 / (9 * df_n))^2
  se  <- sqrt(numerator / denominator)
  z   <- (1 - estimate)^(1 / 3)
  z_lower <- z + crit * se   # large z -> small reliability
  z_upper <- z - crit * se
  list(se    = se * 3 * a_factor^3 * (1 - estimate)^(2 / 3),
       se_transformed     = se,
       se_transform_scale = "cube_root",
       lower = 1 - a_factor^3 * z_lower^3,
       upper = 1 - a_factor^3 * z_upper^3)
}


# .ci_reliability_likelihood(data = NULL, S = NULL, N = NULL,
#                            equal_loadings = FALSE, conf_level = 0.95)
#
# Profile likelihood confidence interval for the model implied
# reliability of a single-factor model, the "likelihood" ci_method of
# reliability_omega() (equal_loadings = FALSE, congeneric) and
# reliability_alpha() (equal_loadings = TRUE, tau-equivalent; the
# profile of the CFA-based alpha).
#
# The interval is the set of population values w0 not rejected by the
# likelihood ratio test: fit the single-factor model unconstrained,
# then refit subject to the nonlinear equality constraint
#
#     (sum lambda_j)^2 / ((sum lambda_j)^2 + sum psi_j) == w0
#
# (lavaan handles the constraint directly), and find by root search the
# w0 on each side of the estimate where
# 2 * (logLik_unconstrained - logLik_constrained) equals the
# chi square critical value with one degree of freedom. The interval
# respects [0, 1] by construction and is not forced to be symmetric,
# unlike the Wald forms.
#
# Inputs:
#   data           : raw items (N rows by J columns), or NULL.
#   S, N           : covariance matrix and sample size, used when data
#                    is NULL (maximum likelihood on summary input).
#   equal_loadings : FALSE for the congeneric model (omega); TRUE for
#                    the tau-equivalent model (CFA-based alpha).
#   conf_level     : confidence level.
#
# Returns:
#   A list with `se` (NA; the method has no standard error), `lower`,
#   and `upper`. Bounds are NA with a warning if a constrained fit
#   fails to converge during the search.
.ci_reliability_likelihood <- function(data = NULL, S = NULL, N = NULL,
                                       equal_loadings = FALSE,
                                       conf_level = 0.95) {
  if (!requireNamespace("lavaan", quietly = TRUE)) {
    stop("The likelihood interval requires the 'lavaan' package.",
         call. = FALSE)
  }

  if (!is.null(data)) {
    data <- as.data.frame(data)
    data <- data[stats::complete.cases(data), , drop = FALSE]
    J <- ncol(data)
    colnames(data) <- paste0("y", seq_len(J))
  } else {
    J <- ncol(S)
    dimnames(S) <- list(paste0("y", seq_len(J)), paste0("y", seq_len(J)))
  }
  vars <- paste0("y", seq_len(J))

  lab <- if (equal_loadings) rep("l1", J) else paste0("l", seq_len(J))
  load_line <- paste(paste0(lab, "*", vars), collapse = " + ")
  err_line  <- paste(paste0(vars, " ~~ t", seq_len(J), "*", vars),
                     collapse = "\n")
  lsum <- paste0("(", paste(unique(lab), collapse = " + "), ")")
  if (equal_loadings) lsum <- paste0("(", J, "*l1)")
  tsum <- paste(paste0("t", seq_len(J)), collapse = " + ")
  base_syntax <- paste0("f1 =~ NA*", vars[1], " + ", load_line,
                        "\nf1 ~~ 1*f1\n", err_line)
  rel_expr <- paste0(lsum, "^2 / (", lsum, "^2 + ", tsum, ")")

  fit_lav <- function(syntax, start = NULL) {
    if (!is.null(data)) {
      lavaan::cfa(syntax, data = data, estimator = "ML", se = "none",
                  start = start)
    } else {
      lavaan::cfa(syntax, sample.cov = S, sample.nobs = N,
                  estimator = "ML", se = "none", start = start)
    }
  }

  fit0 <- fit_lav(base_syntax)
  if (!lavaan::lavInspect(fit0, "converged")) {
    warning("The unconstrained single-factor fit did not converge; ",
            "the likelihood interval is unavailable.", call. = FALSE)
    return(list(se = NA_real_, lower = NA_real_, upper = NA_real_))
  }
  ll0 <- as.numeric(lavaan::logLik(fit0))
  est <- lavaan::parameterEstimates(fit0)
  lam_hat <- est[est$op == "=~", "est"]
  psi_hat <- est[est$op == "~~" & est$lhs == est$rhs &
                 est$lhs %in% vars, "est"]
  w_hat <- (sum(lam_hat))^2 / ((sum(lam_hat))^2 + sum(psi_hat))
  crit <- stats::qchisq(conf_level, df = 1)

  # 2 * (ll0 - ll(w0)) - crit: negative inside the interval, positive
  # outside. Root-find on each side of the estimate.
  g <- function(w0) {
    con <- paste0(base_syntax, "\n", rel_expr, " == ", format(w0, digits = 12))
    fc <- tryCatch(fit_lav(con, start = fit0), error = function(e) NULL)
    if (is.null(fc) || !lavaan::lavInspect(fc, "converged")) return(NA_real_)
    2 * (ll0 - as.numeric(lavaan::logLik(fc))) - crit
  }

  find_bound <- function(direction) {
    # Step outward from the estimate until the profile crosses the
    # critical value, then root-find in the bracket.
    eps <- 1e-6
    w_out <- w_hat
    g_out <- -crit
    for (step in c(0.02, 0.05, 0.1, 0.2, 0.4, 0.8)) {
      w_try <- if (direction < 0) max(w_hat - step, eps)
               else min(w_hat + step, 1 - eps)
      g_try <- g(w_try)
      if (is.na(g_try)) return(NA_real_)
      w_out <- w_try
      g_out <- g_try
      if (g_out > 0) break
      if ((direction < 0 && w_try <= eps) ||
          (direction > 0 && w_try >= 1 - eps)) break
    }
    if (g_out <= 0) {
      # The profile never reaches the critical value before the
      # boundary; the bound is the boundary itself.
      return(if (direction < 0) 0 else 1)
    }
    inner <- if (direction < 0) w_out else w_hat
    outer <- if (direction < 0) w_hat else w_out
    root <- tryCatch(
      stats::uniroot(function(w) g(w), lower = inner, upper = outer,
                     f.lower = if (direction < 0) g_out else -crit,
                     f.upper = if (direction < 0) -crit else g_out,
                     tol = 1e-7)$root,
      error = function(e) NA_real_)
    root
  }

  lower <- find_bound(-1)
  upper <- find_bound(+1)
  if (anyNA(c(lower, upper))) {
    warning("A constrained fit failed to converge during the profile ",
            "likelihood search; the affected bound is NA.", call. = FALSE)
  }
  list(se = NA_real_, lower = lower, upper = upper)
}


# -----------------------------------------------------------------------------
# (C) Delta method CI helpers
#
# Build a confidence interval from an explicit standard error. The SE
# either comes from a closed-form formula (.ci_alpha_ml from van Zyl et
# al. 2000, .ci_alpha_adf from Maydeu-Olivares et al. 2007) or from a
# fitted lavaan model (.ci_omega_delta). When logistic = TRUE, the Wald
# CI is constructed on the logit scale via .logistic_ci() and
# back-transformed.
# -----------------------------------------------------------------------------

# .ci_alpha_ml(estimate, S, N, J, conf_level, logistic = FALSE)
#
# Closed-form maximum likelihood standard error for sample coefficient
# alpha, derived under multivariate normality of the items by van Zyl,
# Neudecker, & Nel (2000). The asymptotic variance is
#
#     Var(alpha_hat) = ( J^2 / (J - 1)^2 ) * ( 2 / (1' R 1)^3 )
#                       * [ (1' R 1) * (tr(R^2) + (tr R)^2)
#                           - 2 * (tr R) * (1' R^2 1) ]
#                       / (N - 1)
#
# where R = cor(S) is the item correlation matrix, 1 is a J-vector of
# ones, and tr() denotes the trace. Wald CI on the reliability scale
# unless logistic = TRUE, in which case the CI is computed on the logit
# scale and back-transformed.
.ci_alpha_ml <- function(estimate, S, N, J, conf_level, logistic = FALSE) {
  R <- stats::cov2cor(S)
  one     <- matrix(1, nrow = J, ncol = 1)
  oneRone <- as.numeric(t(one) %*% R %*% one)
  trR     <- sum(diag(R))
  trR2    <- sum(diag(R %*% R))
  oneR2one <- as.numeric(t(one) %*% (R %*% R) %*% one)

  step1 <- J^2 / (J - 1)^2
  g1 <- 2 / oneRone^3
  g2 <- oneRone * (trR2 + trR^2) - 2 * trR * oneR2one
  variance <- (step1 * g1 * g2) / (N - 1)
  se <- sqrt(variance)
  crit <- stats::qnorm(1 - (1 - conf_level) / 2)

  if (logistic) {
    lg <- .logistic_ci(estimate, se, crit)
    list(se = se, lower = lg$lower, upper = lg$upper)
  } else {
    list(se    = se,
         lower = estimate - crit * se,
         upper = estimate + crit * se)
  }
}


# .ci_alpha_adf(estimate, data, conf_level, logistic = FALSE)
#
# Asymptotic distribution-free (ADF) standard error of coefficient alpha
# (Maydeu-Olivares, Coffman, & Hartmann, 2007). Unlike .ci_alpha_ml,
# this method does NOT assume multivariate normality of the items.
# Requires raw data because the variance is computed via a sandwich
# estimator over the per-person deviations from the sample covariance
# matrix.
#
# Algorithm:
#   1. Compute the gradient (vector "est_delta") of alpha with respect
#      to each unique element of the item covariance matrix S, evaluated
#      at the sample S. This gradient has two values: one for diagonal
#      entries (item variances) and one for off-diagonal entries (item
#      covariances), with closed-form expressions in terms of trace(S)
#      and sum(S).
#   2. For each respondent i, form the contribution s_i to the
#      half-vectorized covariance matrix from that person's centered
#      scores. The deviation s_i - vech(S) summed across persons
#      provides an empirical estimate of the asymptotic covariance of
#      vech(S).
#   3. Sandwich the gradient around that empirical covariance to get
#      Var(alpha_hat).
#
# Wald CI on the reliability scale unless logistic = TRUE.
.ci_alpha_adf <- function(estimate, data, conf_level, logistic = FALSE) {
  d <- as.matrix(stats::na.omit(data))
  N <- nrow(d); J <- ncol(d)
  S <- stats::cov(d)
  trS  <- sum(diag(S))
  totS <- sum(S)

  # vech() returns the column-stacked lower triangle of a matrix
  # (including the diagonal). It is the standard "half-vectorization"
  # operator from multivariate statistics.
  vech <- function(M) {
    unlist(lapply(seq_len(ncol(M)), function(i) M[i:nrow(M), i]))
  }
  vec_S <- vech(S)

  # Gradient of alpha with respect to each unique covariance element.
  off_diag <-  (2 * J / (J - 1)) * (trS / totS^2)
  on_diag  <- -(J / (J - 1)) * ((totS - trS) / totS^2)
  delta_mat <- matrix(off_diag, J, J)
  diag(delta_mat) <- on_diag
  est_delta <- vech(delta_mat)

  d_centered <- scale(d, center = TRUE, scale = FALSE)
  accumulator <- 0
  for (i in seq_len(N)) {
    dy  <- as.matrix(d_centered[i, ])
    s_i <- vech(dy %*% t(dy))
    dsi <- s_i - vec_S
    accumulator <- accumulator +
      as.numeric(t(est_delta) %*% dsi %*% t(dsi) %*% est_delta)
  }
  variance <- accumulator / (N * (N - 1))
  se <- sqrt(variance)
  crit <- stats::qnorm(1 - (1 - conf_level) / 2)

  if (logistic) {
    lg <- .logistic_ci(estimate, se, crit)
    list(se = se, lower = lg$lower, upper = lg$upper)
  } else {
    list(se    = se,
         lower = estimate - crit * se,
         upper = estimate + crit * se)
  }
}


# .ci_omega_delta(fit, conf_level, logistic = FALSE)
#
# Wald (or logistic-transformed Wald) confidence interval for coefficient
# omega using the standard error that lavaan reports for the := defined
# parameter "omega" in the model. The estimator (ML, MLR, WLS) and the
# se= argument were chosen at fit time; this helper just builds the CI.
#
# Inputs:
#   fit        : the list returned by .omega_fit_cfa() with $omega and
#                $se_omega populated.
#   conf_level : confidence level.
#   logistic   : if TRUE, build the CI on the logit scale and
#                back-transform.
.ci_omega_delta <- function(fit, conf_level, logistic = FALSE) {
  if (!fit$converged) {
    return(list(se = NA_real_, lower = NA_real_, upper = NA_real_))
  }
  crit <- stats::qnorm(1 - (1 - conf_level) / 2)
  se  <- fit$se_omega
  est <- fit$omega
  if (logistic) {
    lg <- .logistic_ci(est, se, crit)
    list(se = se, lower = lg$lower, upper = lg$upper)
  } else {
    list(se    = se,
         lower = est - crit * se,
         upper = est + crit * se)
  }
}


# -----------------------------------------------------------------------------
# (D) Logistic transformation
# -----------------------------------------------------------------------------

# .logistic_ci(estimate, se, crit)
#
# Browne's (1982) logit transformation applied to a point estimate and
# standard error on the (0, 1) scale. The transformed estimate is the
# logit r* = log(r / (1 - r)), with delta method SE on the transformed
# scale equal to se / (r * (1 - r)). A Wald CI on the logit scale is
# back-transformed via the inverse logit so that the resulting interval
# always falls inside (0, 1) regardless of how close r is to a boundary.
#
# Inputs:
#   estimate : reliability point estimate in (0, 1).
#   se       : standard error on the reliability scale.
#   crit     : critical value (e.g., qnorm(0.975) for a 95% CI).
# Returns:
#   A list with `lower` and `upper` numeric values. Returns NA endpoints
#   if estimate or se is NA, or if estimate is on the boundary {0, 1}.
.logistic_ci <- function(estimate, se, crit) {
  if (is.na(estimate) || is.na(se) || estimate <= 0 || estimate >= 1) {
    return(list(lower = NA_real_, upper = NA_real_))
  }
  log_est <- log(estimate / (1 - estimate))
  log_se  <- se / (estimate * (1 - estimate))
  log_lower <- log_est - crit * log_se
  log_upper <- log_est + crit * log_se
  list(lower = 1 / (1 + exp(-log_lower)),
       upper = 1 / (1 + exp(-log_upper)))
}


# -----------------------------------------------------------------------------
# (E) Bootstrap engine
# -----------------------------------------------------------------------------

# .bootstrap_ci(data, point_fn, B, conf_level, kind, seed)
#
# Nonparametric bootstrap confidence interval for a reliability
# coefficient. Resamples rows of `data` with replacement B times, applies
# the user-supplied `point_fn` to each resample, and constructs a CI of
# the requested kind:
#
#   "bootstrap_se"          : Wald CI with SE = sd of bootstrap estimates.
#   "bootstrap_se_logistic" : same SE, but built on the logit scale.
#   "percentile"            : empirical 100*(1-C)/2 and 100*(1-(1-C)/2)
#                             quantiles of the bootstrap estimates.
#   "bca"                   : bias-corrected and accelerated bootstrap
#                             (delegates to boot::boot.ci with type = "bca").
#
# Seed handling: if seed is not NULL, set.seed(seed) is called before
# resampling so that runs are reproducible, and the caller's RNG state
# is saved on entry and restored on exit (the package-wide seed
# discipline; see the roxygen of every caller).
#
# Inputs:
#   data       : data frame or matrix passed row-wise to boot::boot.
#   point_fn   : function(d) returning a single numeric estimate from a
#                resample. NA_real_ is allowed (it skips that replication).
#   B          : number of bootstrap replications.
#   conf_level : confidence level.
#   kind       : one of "bootstrap_se", "bootstrap_se_logistic",
#                "percentile", "bca".
#   seed       : integer seed or NULL.
# Returns:
#   A list with `se`, `lower`, `upper`.
.bootstrap_ci <- function(data, point_fn, B, conf_level, kind, seed) {
  if (!requireNamespace("boot", quietly = TRUE)) {
    stop("Bootstrap CIs require the 'boot' package. ",
         "Install with install.packages(\"boot\").", call. = FALSE)
  }
  if (!is.null(seed)) {
    # Save and restore the user's RNG state so that supplying a seed
    # for reproducibility does not pollute their global RNG.
    if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      old_seed <- get(".Random.seed", envir = .GlobalEnv,
                      inherits = FALSE)
      on.exit(assign(".Random.seed", old_seed, envir = .GlobalEnv),
              add = TRUE)
    } else {
      on.exit(rm(".Random.seed", envir = .GlobalEnv), add = TRUE)
    }
    set.seed(as.integer(seed))
  }

  statistic <- function(d, i) {
    val <- tryCatch(point_fn(d[i, , drop = FALSE]),
                    error = function(e) NA_real_)
    if (is.null(val) || length(val) != 1L) NA_real_ else as.numeric(val)
  }
  boot_out <- boot::boot(data = data, statistic = statistic,
                         R = B, stype = "i")

  crit <- stats::qnorm(1 - (1 - conf_level) / 2)
  se  <- stats::sd(boot_out$t, na.rm = TRUE)
  est <- boot_out$t0

  if (kind == "bootstrap_se") {
    lower <- est - crit * se
    upper <- est + crit * se
  } else if (kind == "bootstrap_se_logistic") {
    lg <- .logistic_ci(est, se, crit)
    lower <- lg$lower
    upper <- lg$upper
  } else if (kind == "percentile") {
    # boot.ci() names the component "percent" even for type = "perc";
    # use the full name so the access is exact, not a partial match.
    bci <- boot::boot.ci(boot_out, conf = conf_level, type = "perc")$percent
    lower <- bci[4]
    upper <- bci[5]
  } else if (kind == "bca") {
    bci <- boot::boot.ci(boot_out, conf = conf_level, type = "bca")$bca
    lower <- bci[4]
    upper <- bci[5]
  } else {
    stop("Unknown bootstrap kind: ", kind, call. = FALSE)
  }
  list(se = se, lower = lower, upper = upper)
}


# -----------------------------------------------------------------------------
# (F) Tidy result formatter
# -----------------------------------------------------------------------------

# .relia_result(estimate, se, lower, upper, conf_level, N, J,
#               coefficient, ci_method, B = NA_integer_,
#               N_complete = N, missing = "listwise", aux = NULL,
#               se_transformed = NA_real_, se_transform_scale = NULL)
#
# Pack the point estimate, standard error, CI bounds, and bookkeeping
# values into the package-standard data.frame(term, value).
# Reliability is bounded on [0, 1] by definition, so this helper clamps
# the CI bounds before returning. Metadata that does not belong in the
# numeric `value` column (which coefficient was computed, which CI
# method, the missing-data treatment, the auxiliary variables, and the
# number of bootstrap reps when applicable) travels as attributes.
#
# The `se` row is always on the reliability (coefficient) scale. When
# the interval was built on a transformation scale (the Fisher, Bonett,
# and Hakstian-Whalen methods), the caller passes `se_transformed` (the
# transformation-scale standard error) and `se_transform_scale` (its
# label, e.g. "fisher_z"), and an `se_transformed` row is inserted
# after `se` with the label on the `se_transform_scale` attribute, so
# simulation users keep the quantity the interval is built from.
#
# The table always carries both an `N` row (cases the analysis used;
# under missing = "fiml" that is every case with at least one observed
# item) and an `N_complete` row (cases complete on the items, what
# listwise deletion would keep), so the cost of listwise deletion is
# visible at a glance and the returned shape never depends on how the
# function was called. Under listwise deletion, and for covariance
# input, the two are equal.
.relia_result <- function(estimate, se = NA_real_,
                          lower = NA_real_, upper = NA_real_,
                          conf_level, N, J,
                          coefficient, ci_method, B = NA_integer_,
                          N_complete = N, missing = "listwise",
                          aux = NULL,
                          se_transformed = NA_real_,
                          se_transform_scale = NULL) {
  if (!is.na(lower)) lower <- max(lower, 0)
  if (!is.na(upper)) upper <- min(upper, 1)
  terms  <- c("estimate", "se", "lower_limit", "upper_limit",
              "conf_level", "N", "N_complete", "J")
  values <- c(estimate, se, lower, upper, conf_level, N, N_complete, J)
  if (!is.null(se_transform_scale)) {
    terms  <- append(terms,  "se_transformed", after = 2L)
    values <- append(values, se_transformed,   after = 2L)
  }
  out <- data.frame(term = terms, value = values)
  attr(out, "coefficient") <- coefficient
  attr(out, "ci_method")   <- ci_method
  attr(out, "missing")     <- missing
  if (!is.null(se_transform_scale)) {
    attr(out, "se_transform_scale") <- se_transform_scale
  }
  if (!is.null(aux)) attr(out, "aux") <- aux
  if (!is.na(B)) attr(out, "B") <- B
  # Tag with a leading class so the generics::tidy / generics::glance
  # S3 methods dispatch; data.frame behavior is preserved through
  # inheritance.
  class(out) <- c("dmar_reliability", "dmar_tbl", "data.frame")
  out
}


# -----------------------------------------------------------------------------
# (G) Input-validation helpers
# -----------------------------------------------------------------------------

# .relia_resolve_inputs(data, S, N, allow_S = TRUE,
#                       missing = "listwise", aux = NULL)
#
# Validate and normalize the (data, S, N) inputs that the reliability
# functions accept. The five family functions vary in which inputs they
# allow:
#
#   reliability_alpha       : data OR (S, N)
#   reliability_kr20        : data only
#   reliability_omega       : data OR (S, N); with denominator =
#                             "observed", (S, N) supports the point
#                             estimate only (the bootstrap CI needs data)
#   reliability_omega_categorical     : data only
#
# Setting allow_S = FALSE flips this to "data only" and produces a
# clear error if the caller supplied S.
#
# Missing data handling (raw data only; the caller has already resolved
# the aux-implies-fiml logic and rejected fiml/aux with covariance
# input):
#
#   missing = "listwise" : rows with any NAs are listwise-deleted and
#                          the sample covariance matrix S is computed
#                          once for downstream use (the historical
#                          behavior, unchanged).
#   missing = "fiml"     : rows with at least one observed item are
#                          kept; N is that count and no complete-case S
#                          is computed (S is NULL; the estimation paths
#                          go through lavaan with missing = "ml" or
#                          through .fiml_cov()).
#
# When `aux` names auxiliary columns of `data`, the remaining columns
# are the items; the auxiliaries are validated here (they must exist in
# `data`, be numeric, and not be constant). When S is supplied, it is
# checked for symmetry.
#
# Returns:
#   A list with elements (data, S, N, J, items, aux, missing,
#   N_complete), where exactly one of `data` and `S` is non-NULL
#   according to what the caller supplied. `items` names the item
#   columns of `data` (NULL for covariance input or unnamed data);
#   `N_complete` is the number of cases complete on the items (equal
#   to N under listwise deletion and for covariance input).
.relia_resolve_inputs <- function(data, S, N, allow_S = TRUE,
                                  missing = "listwise", aux = NULL) {
  if (is.null(data) && is.null(S)) {
    stop("Either 'data' (a matrix or data frame of item scores) or ",
         "'S' (a covariance matrix) must be supplied.", call. = FALSE)
  }
  if (!is.null(data) && !is.null(S)) {
    stop("Supply 'data' OR 'S', not both.", call. = FALSE)
  }
  if (!is.null(S) && !allow_S) {
    stop("This coefficient requires raw 'data'; a covariance matrix is ",
         "not sufficient.", call. = FALSE)
  }

  if (!is.null(data)) {
    data <- as.data.frame(data)

    if (!is.null(aux)) {
      if (!is.character(aux) || length(aux) < 1L || anyNA(aux)) {
        stop("'aux' must be a character vector of column names in ",
             "'data'.", call. = FALSE)
      }
      if (anyDuplicated(aux)) {
        stop("'aux' contains duplicated names.", call. = FALSE)
      }
      # as.data.frame() has already guaranteed column names (unnamed
      # matrix input gets V1, V2, ...), so a mismatch surfaces through
      # the not-found error below.
      absent <- setdiff(aux, colnames(data))
      if (length(absent) > 0L) {
        stop("Auxiliary variable(s) not found in 'data': ",
             paste(sQuote(absent), collapse = ", "),
             ". Supply 'data' containing both the items and the ",
             "auxiliary columns; the non-auxiliary columns are the ",
             "items.", call. = FALSE)
      }
      items <- setdiff(colnames(data), aux)
      if (length(items) < 2L) {
        stop("After removing the auxiliary column(s), fewer than two ",
             "item columns remain in 'data'. An auxiliary variable is ",
             "not an item; it cannot make up the composite.",
             call. = FALSE)
      }
      for (a in aux) {
        a_col <- data[[a]]
        if (!is.numeric(a_col)) {
          stop("Auxiliary variable ", sQuote(a), " is not numeric.",
               call. = FALSE)
        }
        a_obs <- a_col[!is.na(a_col)]
        if (length(a_obs) < 2L || stats::var(a_obs) == 0) {
          stop("Auxiliary variable ", sQuote(a), " is constant (or has ",
               "fewer than two observed values); a constant carries no ",
               "information about the items or the missingness.",
               call. = FALSE)
        }
      }
    } else {
      items <- colnames(data)
    }

    if (identical(missing, "fiml")) {
      if (is.null(items)) {
        # Unnamed input; the FIML paths need names for the lavaan
        # syntax, so assign the package's positional convention.
        colnames(data) <- paste0("y", seq_len(ncol(data)))
        items <- colnames(data)
      }
      item_cols <- data[, items, drop = FALSE]
      if (ncol(item_cols) < 2L) {
        stop("At least two items are required.", call. = FALSE)
      }
      any_item_observed <- rowSums(!is.na(item_cols)) > 0L
      data <- data[any_item_observed, , drop = FALSE]
      item_cols <- item_cols[any_item_observed, , drop = FALSE]
      if (nrow(data) < 2L) {
        stop("Fewer than two cases have at least one observed item.",
             call. = FALSE)
      }
      list(data = data,
           S    = NULL,
           N    = nrow(data),
           J    = ncol(item_cols),
           items = items,
           aux  = aux,
           missing = "fiml",
           N_complete = sum(stats::complete.cases(item_cols)))
    } else {
      data <- data[stats::complete.cases(data), , drop = FALSE]
      if (ncol(data) < 2L) {
        stop("At least two items are required.", call. = FALSE)
      }
      list(data = data,
           S    = stats::cov(data),
           N    = nrow(data),
           J    = ncol(data),
           items = items,
           aux  = NULL,
           missing = "listwise",
           N_complete = nrow(data))
    }
  } else {
    if (is.null(N)) {
      stop("'N' must be supplied when 'S' is given.", call. = FALSE)
    }
    if (!isSymmetric(unname(S), tol = 1e-6)) {
      stop("'S' must be a symmetric covariance matrix.", call. = FALSE)
    }
    if (ncol(S) < 2L) {
      stop("At least two items are required.", call. = FALSE)
    }
    list(data = NULL, S = S, N = N, J = ncol(S),
         items = colnames(S), aux = NULL, missing = "listwise",
         N_complete = N)
  }
}
