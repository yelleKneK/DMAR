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
#                                 Kelley & Pornprasertmanit 2016, Eq. 17-18),
#                                 with an unrestricted (the default) or
#                                 model implied polychoric denominator
#     reliability_H()          -- coefficient H, the maximal reliability of
#                                 the optimally weighted composite (Hancock
#                                 & Mueller 2001; Gabler & Raykov 2017)
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
#                                       .omega_observed_parts,
#                                       .omega_observed_fiml_parts,
#                                       .omega_c_from_data
#       Improper single-factor fits     .omega_fit_none, .triad_no_solution,
#                                       .named_cov2cor, .fit_item_cov,
#                                       .is_lavaan_improper_notice,
#                                       .omega_fit_verdict,
#                                       .omega_point_value,
#                                       .omega_point_gate,
#                                       .improper_items_phrase,
#                                       .no_solution_message,
#                                       .negative_error_message
#       Categorical-omega support       .omega_c_parts,
#                                       .stop_omega_c_two_items,
#                                       .improper_message,
#                                       .gy_omega_from_components,
#                                       .omega_categorical_from_kfit,
#                                       .kfit_improper_message,
#                                       .gauss_legendre,
#                                       .GAUSS_LEGENDRE, .pbinorm,
#                                       .pbinorm_genz,
#                                       .ordered_thresholds,
#                                       .polychoric_from_fit,
#                                       .undetermined_polychorics,
#                                       .undetermined_message,
#                                       .omega_c_identification,
#                                       .omega_c_unidentified_message,
#                                       .omega_c_undetermined_note
#       Local dependence and labels     .residual_cor_summary,
#                                       .omega_residual_cor,
#                                       .omega_ml_residual_cor,
#                                       .denominator_phrase
#       Missing data (FIML) support     .aux_saturated_syntax, .fiml_cov
#   (B) Closed-form CI helpers          .ci_feldt, .ci_fisher, .ci_bonett,
#                                       .ci_hakstian_whalen
#   (C) Delta method CI helpers         .alpha_vanzyl_avar, .ci_alpha_ml,
#                                       .ci_alpha_adf, .ci_omega_delta
#   (D) Logistic transformation         .logistic_ci
#   (E) Bootstrap engine                .bootstrap_ci, .BOOT_STATUS,
#                                       .boot_statistic, .boot_replicate,
#                                       .boot_classify, .boot_tally,
#                                       .boot_refit_errors,
#                                       .boot_refit_one,
#                                       .boot_jackknife, .boot_rows,
#                                       .boot_with_progress,
#                                       .boot_on_rows, .boot_lapply,
#                                       .boot_worker, .boot_worker_stage,
#                                       .boot_worker_run,
#                                       .boot_block_size, .boot_progress,
#                                       .check_progress,
#                                       .boot_leave_one_out, .boot_caution,
#                                       .check_ncpus,
#                                       .boot_parallel_type, .pct,
#                                       .bootstrap_reason_phrase,
#                                       .bootstrap_footer_line
#   (F) Tidy result formatter           .relia_result
#   (G) Input-validation helpers        .relia_declared_missing,
#                                       .relia_resolve_inputs
#   (H) Coefficient H                   .H_warn_marginal,
#                                       .H_setup, .H_eval,
#                                       .H_est_at, .H_from_fit, .H_replicate,
#                                       .H_from_loadings, .H_refit_fun,
#                                       .H_improper_note,
#                                       .H_footer_lines, .H_finish, .H_wald,
#                                       .H_attributes, .H_fit_data
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
#   A single numeric value, the sample coefficient alpha. The value is
#   at most one but can fall below zero, when the item covariances are
#   negative on balance (most often an item keyed in the opposite
#   direction to the others). .relia_result() clamps each confidence
#   limit into [0, 1] and warns about a negative estimate; the estimate
#   itself is returned as computed. A composite that does not vary, whose
#   variance sum(S) is zero up to rounding error, is refused
#   (.relia_check_composite()): its reliability is not defined, and the
#   ratio above has no finite value. Through version 1.0.1 the ratio was
#   computed, giving -Inf (KR-20 of two rows on which as many items are
#   passed one way as the other) or, where rounding error left sum(S) a
#   tiny number of either sign, an estimate of the order of 1e16 (alpha of
#   ipsative items, whose scores sum to the same total in every row).
#   Every path that computes alpha from a matrix goes through here: the
#   analytic estimate from data, from S, and from the FIML or two-stage
#   moments, KR-20, and every bootstrap replication of them, where the
#   refusal is an error the engine counts and quotes. A matrix the EM
#   algorithm estimated (`em`, the FIML and two-stage moments) is accurate
#   only to the algorithm's convergence, not to rounding error, and is held
#   to the looser tolerance of .relia_composite_varies().
#
# Inputs:
#   S  : the covariance matrix.
#   em : TRUE when S is an EM estimate.
.alpha_from_S <- function(S, em = FALSE) {
  J <- nrow(S)
  .relia_check_composite(S, em = em)
  (J / (J - 1)) * (1 - sum(diag(S)) / sum(S))
}


# .relia_composite_varies(S, em = FALSE)
#
# Whether the composite (the unweighted sum of the items) varies under the
# covariance matrix S: whether its variance, sum(S), is above a tolerance
# times the sum of the item variances. A sample covariance matrix, or one
# the caller supplies, is exact up to rounding error, and the tolerance is
# sqrt(.Machine$double.eps), about 1.5e-8. An estimate of the EM algorithm
# (`em`) is accurate only to its convergence criterion: when the data hold
# the composite constant in every complete row (ipsative items with some
# values missing), the saturated likelihood grows without bound as the
# composite's variance goes to zero, and the EM algorithm drifts there
# and stops with that variance of the order of 1e-8 to 1e-7 times the sum
# of the item variances, at or above the rounding tolerance, so that the
# analytic alpha under FIML came back as -5e7 or -8e7 rather than refused
# (with 10 to 15 percent of the values of three to five ipsative
# holzinger_swineford tests missing). For an EM estimate the tolerance is
# 1e-6: a composite whose variance is a millionth of the sum of its item
# variances has an alpha below -1e6 J/(J - 1), which no composite with a
# reliability gives, and the drift above stays 15 times or more below it.
# A sum that is not finite is left to the caller.
.relia_composite_varies <- function(S, em = FALSE) {
  total <- sum(S)
  tol <- if (em) 1e-6 else sqrt(.Machine$double.eps)
  !is.finite(total) || total > tol * sum(diag(S))
}


# .stop_relia_constant_composite(em = FALSE, rows = NULL, missing = NULL)
#
# The refusal of a composite (the unweighted sum of the items) that does
# not vary: its variance, the sum of the entries of the covariance
# matrix, is zero up to rounding error (or, for an EM estimate, up to the
# precision of the estimate; .relia_composite_varies()). That happens when
# the items' covariances cancel their variances exactly, as with ipsative
# items, whose scores sum to the same total in every row, or two rows on
# which as many items move one way as the other. The reliability of a
# constant is not defined (a ratio of true score variance to observed
# variance with both zero), and coefficient alpha, J/(J - 1) times one
# minus the ratio of the item variances to the composite's, has no finite
# value; the single-factor fits stop on the same data, whose covariance
# matrix is singular. With `rows`, the check was made on the rows that
# observe every item under full information maximum likelihood or
# two-stage estimation (`missing`), where the composite of those rows is
# all the data say about it: when it is constant there, the saturated
# likelihood grows without bound as the estimate of its variance goes to
# zero, so those methods have no estimate to give. The message is one
# sentence with no line break, since the bootstrap engine quotes the
# first line of a replication's error. Called through
# .relia_check_composite().
.stop_relia_constant_composite <- function(em = FALSE, rows = NULL,
                                           missing = NULL) {
  ipsative <- paste0("as when the items are ipsative (their scores sum to ",
                     "the same total in every row)")
  undefined <- paste0("so the reliability of the composite is not defined ",
                      "(alpha, J/(J - 1) times one minus the ratio of the ",
                      "item variances to that variance, has no finite ",
                      "value")
  .dmar_refuse_argument(
    if (!is.null(rows)) {
 paste0("The composite, the sum of the items, does not vary in the ",
        rows, " rows that observe every item: its variance there is ",
        "zero up to rounding error, ", ipsative, ", ", undefined,
        "), and the likelihood of the saturated model that ",
        .relia_incomplete_label(missing), " estimates first grows ",
        "without bound as its estimate of that variance goes to zero, ",
        "so the method has no estimate to give either.")
} else if (em) {
 paste0("The composite, the sum of the items, does not vary: its ",
        "variance, the sum of the entries of the covariance matrix the ",
        "EM algorithm estimates, is below a millionth of the sum of the ",
        "item variances, zero up to the precision of that estimate, ",
        ipsative, ", ", undefined, ", and the single-factor model ",
        "cannot be fit to the singular covariance matrix).")
} else {
 paste0("The composite, the sum of the items, does not vary: its ",
        "variance, the sum of the entries of the covariance matrix, is ",
        "zero up to rounding error, ", ipsative, ", ", undefined,
        ", and the single-factor model cannot be fit to the singular ",
        "covariance matrix).")
}, arg = NULL)
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
    .dmar_refuse_argument(
      "KR-20 requires items scored 0/1; non-binary values were found.",
      arg = NULL)
  }
  data <- data[stats::complete.cases(data), , drop = FALSE]
  .alpha_from_S(stats::cov(data))
}


# .omega_fit_cfa(data = NULL, S = NULL, N = NULL,
#                equal_loadings = FALSE,
#                estimator = "ML", se = "standard", missing = "listwise",
#                aux = NULL, h1_cov = NULL)
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
# items) under FIML, cfa_1() cannot be used, because the saturated
# correlates block is not part of the model cfa_1()/fa_confirmatory() build. The
# helper then builds the same single-factor measurement syntax, appends
# the saturated correlates block from .aux_saturated_syntax(), and calls
# lavaan directly with missing = "ml". Phase 2 of the missing-data work
# (cfa_1() and fa_confirmatory() themselves) can lift this path into fa_confirmatory()
# proper. Under two-stage estimation the auxiliaries need no block in the
# model, since they enter the first stage only, so cfa_1() is used with
# lavaan's aux argument passed through.
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
#   data           : raw items (N rows by J columns); supply this or S+N.
#   S, N           : a J x J covariance matrix and a sample size.
#   equal_loadings : if TRUE the fit imposes tau-equivalence; default
#                    FALSE (congeneric model).
#   estimator      : forwarded to cfa_1(); one of "ML" (default), "MLR",
#                    "WLS", "WLSMV", "GLS".
#   se             : forwarded to cfa_1(); typically "standard",
#                    "robust.sem", or "none".
#   missing        : forwarded to cfa_1() and used by lavaan only when
#                    raw data are supplied; "listwise" (default), "ml",
#                    "fiml", "pairwise", or lavaan's two-stage values
#                    "two.stage" and "robust.two.stage" (see
#                    .relia_lavaan_missing()).
#   aux            : optional character vector naming auxiliary columns
#                    of `data`; requires raw data. Under missing = "ml"
#                    they enter as saturated correlates
#                    (.omega_fit_cfa_aux()); under two-stage estimation
#                    they enter stage one only, through lavaan's aux
#                    argument, and the measurement model is fit as without
#                    them. The remaining columns are the items.
#   h1_cov         : optional; under FIML or two-stage estimation, the
#                    item block of the saturated (EM) covariance matrix
#                    (.fiml_cov()$S, or the stage one moments) when the
#                    caller has already computed it, which the check for
#                    three items with no solution then uses rather than
#                    estimating it again.
#
# Returns:
#   A list with elements:
#     converged   : logical, whether the fit converged.
#     improper    : logical, whether some error variance estimate is
#                   negative (a Heywood case); see .omega_fit_verdict().
#     no_solution : logical, TRUE when three items under the congeneric
#                   model have covariances whose product is not positive,
#                   so no single-factor solution exists and no fit is
#                   attempted (or the fit is discarded); converged is
#                   then FALSE and improper TRUE. See
#                   .triad_no_solution().
#     triad       : the 3 x 3 correlation matrix behind a no_solution
#                   verdict, which the refusal quotes (NULL otherwise).
#     omega       : numeric, sample coefficient omega.
#     se_omega    : numeric, lavaan's delta method SE for the omega
#                   defined parameter (NA_real_ when se = "none" or the
#                   SE could not be computed).
#     loadings    : numeric vector of length J, the lambda estimates.
#     errors      : numeric vector of length J, the psi^2 estimates.
#     items       : character vector of the item names, parallel to
#                   loadings and errors.
#     lavaan_items: the same items as the lavaan fit names them, y1,
#                   y2, ..., in both branches (absent when there is no
#                   fit).
#     J           : integer, the number of items.
#     fit         : the lavaan fit object (or NULL when the fit failed).
#     fit_error   : when lavaan stopped with an error, its message
#                   (absent otherwise; a fit that ran but did not
#                   converge has none).
#     coef_name   : the name of the defined parameter holding omega in
#                   the fit ("omega_f1", or "omega" in the auxiliary
#                   variable branch), which a Wald interval's wald_df
#                   resolution looks up (absent when there is no fit).
#     no_maximum  : when the likelihood of the fit has no maximum (an
#                   error variance that runs off below zero, the
#                   likelihood flat or still rising;
#                   .lavaan_no_maximum()), the refusal's fields with the
#                   caller's names for the items (.no_max_renamed()),
#                   beside no fit, as .omega_fit_no_maximum() records it;
#                   absent otherwise.
.omega_fit_cfa <- function(data = NULL, S = NULL, N = NULL,
                           equal_loadings = FALSE,
                           estimator = "ML",
                           se = "standard",
                           missing = "listwise",
                           aux = NULL,
                           h1_cov = NULL, h1_diagnostic = NULL) {
  if (!requireNamespace("lavaan", quietly = TRUE)) {
    .dmar_refuse(c(
      "Coefficient omega requires the 'lavaan' package. ",
      "Install with install.packages(\"lavaan\")."),
      "package_required", package = "lavaan")
  }

  two_stage <- missing %in% c("two.stage", "robust.two.stage")
  numeric_data <- is.null(data) || if (is.data.frame(data)) {
    all(vapply(data, is.numeric, logical(1L)))
  } else is.numeric(data)
  if (is.null(h1_diagnostic) && !is.null(data) && numeric_data &&
      missing %in% c("ml", "fiml")) {
    h1_diagnostic <- .fiml_data_h1_diagnostic(data)
  }
  # Under FIML and two-stage estimation a pair of variables never observed
  # together has a covariance the data do not determine, and a pair
  # observed together only where one of them takes a single value has one
  # they do not inform; both are refused by name before any fit
  # (.stop_never_observed_together(), .stop_single_valued_together()).
  if (!is.null(data) && two_stage) {
    label <- .relia_incomplete_label(if (two_stage) "two_stage" else "fiml")
    .stop_never_observed_together(data, label)
    .stop_single_valued_together(data, label)
  }
  if (!is.null(aux)) {
    if (is.null(data)) {
      .dmar_refuse_argument(
        "Auxiliary variables require raw 'data'.", arg = "data")
    }
    # The items and the auxiliary variables are selected by name, so a
    # name given to two columns is refused, whichever they are
    # (.stop_relia_duplicated_names()).
    if (anyDuplicated(colnames(data))) {
      .stop_relia_duplicated_names(colnames(data), aux)
    }
    if (!two_stage) {
      return(.omega_fit_cfa_aux(data = data, aux = aux,
                                equal_loadings = equal_loadings,
                                estimator = estimator, se = se,
                                h1_cov = h1_cov, h1_diagnostic = h1_diagnostic))
    }
  }

  # Resolve item names (used when extracting estimates) and figure out J.
  # Under two-stage estimation with auxiliary variables the data hold both,
  # and the items are the columns not named in `aux`.
  if (!is.null(data)) {
    data <- as.data.frame(data)
    if (is.null(colnames(data))) {
      colnames(data) <- paste0("y", seq_len(ncol(data)))
    }
    var_names <- colnames(data)[!colnames(data) %in% aux]
  } else {
    if (is.null(S) || is.null(N)) {
      .dmar_refuse_argument(
        "Either 'data' or both 'S' and 'N' must be supplied.",
        arg = c("data", "S", "N"))
    }
    var_names <- colnames(S)
    if (is.null(var_names)) var_names <- paste0("y", seq_len(ncol(S)))
  }
  J <- length(var_names)
  # The items enter the fit by position, so two items with one name
  # could not be told apart in a message (and a lookup by name would
  # find the first twice); they are refused by name instead.
  if (anyDuplicated(var_names)) .stop_relia_duplicated_names(var_names)

  # The items enter lavaan as y1, ..., yJ and the auxiliary variables as
  # a1, ..., aK, as in .omega_fit_cfa_aux() and .fiml_cov(), so a column
  # name that is not a valid lavaan identifier ("item 1", "1a") cannot
  # break the model syntax. cfa_1() keeps the names it is given, and
  # through DMAR 1.0.1 they were the caller's: lavaan then stopped on such
  # a name, the try() below caught it, and the call reported a model that
  # did not converge. The caller's names come back in `items` and in
  # every message, lavaan's own error among them
  # (.relia_lavaan_names_back()) and the warnings lavaan raises on the way
  # (.relia_with_user_names(), around everything below that reaches
  # lavaan: the saturated estimate of the check for three items, whose EM
  # names the columns v1, v2, ..., and the fit); `lavaan_items` holds the
  # fit's names.
  item_labels <- paste0("y", seq_len(J))
  aux_labels <- if (length(aux)) paste0("a", seq_along(aux)) else NULL
  if (!is.null(data)) {
    cfa_input <- data[, c(var_names, aux), drop = FALSE]
    colnames(cfa_input) <- c(item_labels, aux_labels)
    cfa_N     <- NULL
  } else {
    cfa_input <- S
    dimnames(cfa_input) <- list(item_labels, item_labels)
    cfa_N     <- N
  }
  user_named <- function(expr) {
    .relia_with_user_names(
      expr, c(item_labels, aux_labels, paste0("v", seq_len(J + length(aux)))),
      rep(c(var_names, aux), 2L))
  }

  # Three items under the congeneric model: the model is saturated, and
  # it has a solution only when the product of the three covariances is
  # positive (see .triad_no_solution()). The check is made before
  # fitting, on the covariance matrix the fit would reproduce: the one
  # supplied, that of the complete cases, or under FIML the saturated
  # (EM) estimate, one quick fit of a saturated model. That spares the
  # point estimate, and every bootstrap replication with such a
  # resample, the thousands of iterations lavaan spends before it either
  # stops or reports a solution that has run away (about 4.6 s for a
  # FIML fit to 300 cases). If the saturated fit fails, the check is
  # made on the fit's own sample statistics below.
  triad_checked <- FALSE
  if (!equal_loadings && J == 3L && is.null(h1_diagnostic)) {
    S_check <- user_named(if (!is.null(cfa_N)) {
      cfa_input
    } else if (identical(missing, "listwise")) {
      cc <- stats::na.omit(as.matrix(cfa_input))
      if (nrow(cc) >= 2L) {
        tryCatch(stats::cov(cc), error = function(e) NULL)
      } else {
        NULL
      }
    } else if (missing %in% c("ml", "fiml")) {
      if (!is.null(h1_cov)) {
        h1_cov
      } else {
        fc <- .fiml_cov(cfa_input, items = item_labels, se = "none")
        if (isTRUE(fc$converged)) fc$S else NULL
      }
    } else if (two_stage) {
      # The stage one (EM) moments the two-stage fit will use, over the
      # items and any auxiliaries, by lavaan's own EM step and without a
      # fit (.ts_em_moments()).
      if (!is.null(h1_cov)) {
        h1_cov
      } else {
        em <- tryCatch(.ts_em_moments(cfa_input),
                       error = function(e) NULL)
        if (is.null(em)) NULL else em$cov[seq_len(J), seq_len(J)]
      }
    } else {
      NULL
    })
    if (!is.null(S_check)) {
      triad_checked <- TRUE
      if (.triad_no_solution(S_check)) {
        return(.omega_fit_none(J, var_names, no_solution = TRUE,
                               triad = .named_cov2cor(S_check, var_names)))
      }
    }
  }

  # cfa_1() warns on an improper (Heywood) solution rather than stopping.
  # Catch that warning here instead of letting it out: this helper is also
  # the inner loop of the reliability bootstrap, where one warning per
  # replication would bury the user. The fact travels back on the returned
  # list as `improper`, so the single-fit callers can warn once (or
  # refuse; see .omega_point_gate()) and the bootstrap can count, through
  # the `dmar_reason` its point functions attach (.omega_point_value()).
  # lavaan's own notices of the same fact (its post-estimation check of
  # negative variances, and, since lavaan 0.7, its note that a solution
  # "seems to have run away") are muffled for the same reason, as
  # .omega_c_parts() muffles the first for the categorical fit. The
  # handler records the fact on a small state environment, the package's
  # idiom for a handler that must report back to the function that
  # installed it.
  state <- new.env(parent = emptyenv())
  state$improper <- FALSE
  fit <- tryCatch(user_named(
    withCallingHandlers(
      do.call(cfa_1, c(
        list(data = if (is.null(cfa_N)) cfa_input else NULL,
             items = if (is.null(cfa_N)) item_labels else NULL,
             S = if (is.null(cfa_N)) NULL else cfa_input, N = cfa_N,
             equal_loading = equal_loadings,
             equal_error   = FALSE,
             estimator = estimator,
             se        = se,
             missing   = missing,
             output    = "fit"),
        # lavaan's aux argument, through cfa_1() and fa_confirmatory(): under
        # two-stage estimation the auxiliaries enter stage one only.
        if (!is.null(aux)) list(aux = aux_labels))),
      dmar_heywood_warning = function(w) {
        state$improper <- TRUE
        invokeRestart("muffleWarning")
      },
      dmar_h1_unavailable = function(w) {
        if (inherits(w, "warning")) invokeRestart("muffleWarning")
      },
      warning = function(w) {
        if (.is_lavaan_improper_notice(w)) invokeRestart("muffleWarning")
      },
      # Under two-stage estimation fa_confirmatory() says that it kept the cases
      # observed on an auxiliary variable alone; the family says so once
      # per call (.relia_resolve_inputs()), not once per fit.
      dmar_two_stage_aux_only_message = function(m) {
        invokeRestart("muffleMessage")
      },
      # When fa_confirmatory() could not compute the two-stage fit again with every
      # case its first stage read, its warning speaks of its own table;
      # the family's says what that means for the coefficient.
      dmar_two_stage_n = function(w) {
        .ts_coefficient_n_warning(w$read, w$nobs, w$extra)
        invokeRestart("muffleWarning")
      }
    )),
    error = function(e) e)

  # A fit that stopped with an error is kept apart from one that did not
  # converge. fa_confirmatory() stops on both, with classed errors: a fit that ran
  # and did not converge (dmar_not_converged) is reported as such, with the
  # advice for an incomplete data method (.relia_fit_failure()), and
  # lavaan's own error (dmar_lavaan_error) travels back as `fit_error`,
  # without fa_confirmatory()'s "Model fitting with lavaan failed: " and with the
  # variables given the caller's names, for the callers to quote. Through
  # the first versions of 1.1.0 every error was recorded as lavaan's, so a
  # fit that did not converge was reported as "lavaan stopped with this
  # error: The model did not converge." and lost that advice. Any other
  # error (one of fa_confirmatory()'s own checks) is quoted as DMAR's.
  # A fit whose likelihood has no maximum (the engine's refusal,
  # .lavaan_no_maximum()) is kept apart as well, with the caller's names,
  # for the point estimate to refuse under the name of the function the
  # user called and a bootstrap replication to be counted
  # (.omega_fit_no_maximum()).
  if (inherits(fit, "dmar_no_maximum")) {
    return(.omega_fit_no_maximum(J, var_names, .no_max_renamed(
      .no_max_fields(fit), c(item_labels, aux_labels), c(var_names, aux))))
  }
  if (inherits(fit, "dmar_not_converged")) {
    return(.omega_fit_none(J, var_names))
  }
  if (inherits(fit, "error")) {
    from_lavaan <- inherits(fit, "dmar_lavaan_error")
    return(.omega_fit_none(J, var_names, error = .relia_lavaan_names_back(
      if (from_lavaan) fit$lavaan_message else conditionMessage(fit),
      c(item_labels, aux_labels), c(var_names, aux)),
      from_lavaan = from_lavaan))
  }
  if (isFALSE(lavaan::lavInspect(fit, "converged"))) {
    return(.omega_fit_none(J, var_names))
  }

  if (!equal_loadings && J == 3L && !triad_checked && is.null(h1_diagnostic)) {
    S_h1 <- .fit_item_cov(fit, item_labels)
    if (!is.null(S_h1) && .triad_no_solution(S_h1)) {
      return(.omega_fit_none(J, var_names, no_solution = TRUE,
                             triad = .named_cov2cor(S_h1, var_names)))
    }
  }

  pe <- lavaan::parameterEstimates(fit)
  om_row   <- pe[pe$lhs == "omega_f1" & pe$op == ":=", , drop = FALSE]
  ld_rows  <- pe[pe$op == "=~", "est"]
  err_rows <- pe[pe$op == "~~" & pe$lhs == pe$rhs &
                 pe$lhs %in% item_labels, "est"]
  se_omega <- .lavaan_defined_se(fit, om_row$se)

  list(converged = TRUE,
       h1_diagnostic = h1_diagnostic,
       improper = state$improper,
       no_solution = FALSE,
       triad = NULL,
       omega = om_row$est[1],
       se_omega = se_omega,
       loadings = ld_rows,
       errors   = err_rows,
       items = var_names,
       lavaan_items = item_labels,
       J = J,
       fit = fit,
       coef_name = "omega_f1")
}


# .omega_fit_none(J, items, no_solution = FALSE, triad = NULL,
#                 error = NULL, from_lavaan = TRUE)
#
# The list .omega_fit_cfa() and .omega_fit_cfa_aux() return when there is
# no fit to report: a fit that failed or did not converge, or (with
# no_solution = TRUE) three items whose covariances admit no
# single-factor solution, which is improper by construction. `error` is
# the message of a fit that stopped with an error, returned as
# `fit_error`, which the callers' refusals quote (.relia_fit_failure());
# `from_lavaan` says whether lavaan raised it (FALSE for one of fa_confirmatory()'s
# own checks), returned as `fit_error_lavaan` beside it.
.omega_fit_none <- function(J, items, no_solution = FALSE, triad = NULL,
                            error = NULL, from_lavaan = TRUE) {
  out <- list(converged = FALSE, improper = no_solution,
              no_solution = no_solution, triad = triad,
              omega = NA_real_, se_omega = NA_real_,
              loadings = rep(NA_real_, J),
              errors   = rep(NA_real_, J),
              items = items,
              J = J, fit = NULL)
  if (!is.null(error)) {
    out$fit_error <- error
    out$fit_error_lavaan <- from_lavaan
  }
  out
}


# .omega_fit_no_maximum(J, items, nm)
#
# The list .omega_fit_cfa() and .omega_fit_cfa_aux() return for a fit
# whose likelihood has no maximum (.lavaan_no_maximum(); `nm`, its fields
# with the caller's names, .no_max_renamed()): no fit, as .omega_fit_none()
# records one, with the short reason as `fit_error` (so a caller that
# reads only that, the planners among them, quotes it), and `no_maximum`,
# from which the point estimate's gate raises the refusal under the name
# of the function the user called (.omega_point_gate()) and a bootstrap
# replication raises its short stop, which the engine counts
# (.omega_point_value()).
.omega_fit_no_maximum <- function(J, items, nm) {
  out <- .omega_fit_none(J, items, error = .no_max_short(nm),
                         from_lavaan = FALSE)
  out$no_maximum <- nm
  out
}


# .relia_lavaan_names_back(text, lavaan_names, names)
#
# A lavaan message about a reliability fit with the variables given the
# caller's names. The fits enter the items into lavaan as y1, ..., yJ and
# the auxiliary variables as a1, ..., aK (.omega_fit_cfa(),
# .omega_fit_cfa_aux(), .omega_c_parts()), so a column name that is not
# a valid lavaan identifier cannot break the model syntax, and an error
# lavaan stops with names the variables that way ("Some model variables
# are of type character: \"y2\""). Through the first versions of 1.1.0
# the error was quoted as lavaan gave it (.relia_fit_failure()), which
# left the user to work out which column y2 was. Each name is replaced as
# a whole word, a maximal run of letters, digits, and underscores, so y10
# is never read as y1 followed by 0, and all of them at once, so items
# the caller named y2 and y1, in that order, are not exchanged twice. A
# word that is not one of `lavaan_names` (lavaan's own function names,
# a parameter label) is left as it is.
#
# Inputs:
#   text         : the message (NULL passes through).
#   lavaan_names : the names the fit gave the variables.
#   names        : the caller's names, parallel to `lavaan_names`.
# Returns:
#   The message with each of `lavaan_names` replaced by its counterpart
#   in `names`.
.relia_lavaan_names_back <- function(text, lavaan_names, names) {
  if (is.null(text) || !length(lavaan_names)) return(text)
  lookup <- stats::setNames(as.character(names), lavaan_names)
  words <- gregexpr("[A-Za-z0-9_]+", text, perl = TRUE)
  regmatches(text, words) <- lapply(regmatches(text, words), function(w) {
    hit <- w %in% lavaan_names
    w[hit] <- lookup[w[hit]]
    w
  })
  text
}


# .relia_with_user_names(expr, lavaan_names, names)
#
# Evaluates `expr`, a call that reaches lavaan with the variables under
# the fit's names, and raises each warning lavaan gives on the way again
# with the caller's names (.relia_lavaan_names_back()), muffling the
# original; a warning that names no variable passes as it is. Through
# DMAR 1.0.1 the single-factor fits passed the caller's names to lavaan,
# so its warnings named the columns as the data did ("some observed
# variables are perfectly correlated; ... variables involved are: dep
# wor", the only place lavaan names the items its refusal that follows is
# about); since the items entered the fit as y1, ..., yJ the same
# warning said "y2 y4", and
# the saturated estimate under FIML, whose EM names the columns v1, v2,
# ..., said "v2 v4" (as it had through 1.0.1). The handler maps only the
# names of the call it wraps: a call that already gives the caller's
# names (an inner fit that maps its own) is left outside it, since
# mapping such a message again could exchange a caller's name that
# happens to look like a fit's name ("y2") for another.
#
# Inputs:
#   expr         : the call, evaluated lazily inside the handler.
#   lavaan_names : the names the call gives the variables.
#   names        : the caller's names, parallel to `lavaan_names`.
# Returns:
#   The value of `expr`.
.relia_with_user_names <- function(expr, lavaan_names, names) {
  withCallingHandlers(expr, warning = function(w) {
    .relia_warn_named(w, .relia_lavaan_names_back(conditionMessage(w),
                                                  lavaan_names, names))
  })
}


# .relia_warn_named(w, named)
#
# Inside a calling handler for the warning `w`: raises it again with the
# message `named` (keeping its class and call) and muffles the original,
# when `named` differs from its message; otherwise does nothing, and the
# warning goes on as it is.
.relia_warn_named <- function(w, named) {
  if (identical(named, conditionMessage(w))) return(invisible(NULL))
  w$message <- named
  warning(w)
  invokeRestart("muffleWarning")
}


# .triad_no_solution(S)
#
# Whether three items have no single-factor solution. A single factor
# implies sigma_12 = lambda_1 lambda_2, and so on, so the product of the
# three covariances is (lambda_1 lambda_2 lambda_3)^2, which cannot be
# negative, and the loadings are determined by it, lambda_1^2 =
# sigma_12 sigma_13 / sigma_23 (cyclically). When the sample product is
# negative no real loadings reproduce the covariances, so the likelihood
# has no maximum at any finite point: its supremum is approached as one
# loading grows without bound while its residual variance falls toward
# minus infinity, and an optimizer either fails to converge or stops
# somewhere along that path (the "run away" solution lavaan 0.7 names).
# A product of zero leaves a loading undetermined. Either way the model
# yields no estimate.
#
# Inputs:
#   S : a 3 x 3 covariance (or correlation) matrix.
# Returns:
#   TRUE when the product of the three off-diagonal elements is not
#   positive; FALSE otherwise, or when it is not finite.
.triad_no_solution <- function(S) {
  S <- as.matrix(S)
  p <- S[1L, 2L] * S[1L, 3L] * S[2L, 3L]
  is.finite(p) && p <= 0
}


# .named_cov2cor(S, items)
#
# The correlation matrix of S with the item names as dimnames, for the
# messages that quote it.
.named_cov2cor <- function(S, items) {
  R <- stats::cov2cor(unname(as.matrix(S)))
  dimnames(R) <- list(items, items)
  R
}


# .fit_item_cov(fit, items)
#
# The covariance matrix of the items the fit was estimated from: the
# sample covariance matrix for complete data, and under missing = "ml"
# the saturated (EM) estimate lavaan fits first, which is what
# lavInspect(fit, "sampstat") returns there. NULL when lavaan cannot
# supply it.
.fit_item_cov <- function(fit, items) {
  S <- tryCatch(lavaan::lavInspect(fit, "sampstat")$cov,
                error = function(e) NULL)
  if (is.null(S) || !all(items %in% rownames(S))) return(NULL)
  unclass(S)[items, items, drop = FALSE]
}


# .is_lavaan_improper_notice(w)
#
# TRUE for the two lavaan warnings that announce an improper solution,
# which the reliability fits report through their own `improper` flag:
# the post-estimation check ("some estimated ov variances are negative")
# and the note that the solution "seems to have run away" (lavaan 0.7).
# lavaan wraps its messages to the console width, so the text is matched
# with its white space collapsed.
.is_lavaan_improper_notice <- function(w) {
  msg <- gsub("[[:space:]]+", " ", conditionMessage(w))
  grepl("estimated ov variances are negative", msg, fixed = TRUE) ||
    grepl("seems to have run away", msg, fixed = TRUE)
}


# .omega_fit_cfa_aux(data, aux, equal_loadings, estimator, se,
#                    h1_cov = NULL)
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
# Unlike the cfa_1() route, the syntax defines omega but not H or the
# average variance extracted (fa_confirmatory()'s H_f1 and ave_f1): reliability_H()
# reads H from the fit's loadings and residual variances
# (.H_from_fit()), not from a defined parameter, and nothing here needs
# the AVE, so a caller that wants either from this fit computes it.
#
# Inputs:
#   data           : data frame containing both item and auxiliary columns.
#   aux            : character vector naming the auxiliary columns; the
#                    remaining columns are the items, in column order.
#   equal_loadings : TRUE imposes tau-equivalence.
#   estimator      : "ML" or "MLR" (WLS cannot be combined with FIML;
#                    callers validate before reaching this point).
#   se             : lavaan se type ("standard", "robust.sem", "none").
#   h1_cov         : optional item block of the saturated (EM) covariance
#                    matrix over the items and auxiliaries, as in
#                    .omega_fit_cfa().
#
# Returns:
#   The same structured list as .omega_fit_cfa().
.omega_fit_cfa_aux <- function(data, aux, equal_loadings = FALSE,
                               estimator = "ML", se = "standard",
                               h1_cov = NULL, h1_diagnostic = NULL) {
  data <- as.data.frame(data)
  # The columns are selected by name, so a name given to two of them is
  # refused (setdiff() kept one of each, and the other was dropped).
  if (anyDuplicated(colnames(data))) {
    .stop_relia_duplicated_names(colnames(data), aux)
  }
  items <- colnames(data)[!colnames(data) %in% aux]
  J <- length(items)
  K <- length(aux)

  # Rename to positional labels, as .omega_fit_cfa() and .fiml_cov() do,
  # so the syntax never trips over user column names that are not valid
  # lavaan identifiers.
  y_names <- paste0("y", seq_len(J))
  a_names <- paste0("a", seq_len(K))
  d <- data[, c(items, aux), drop = FALSE]
  colnames(d) <- c(y_names, a_names)
  # lavaan's warnings name the variables as the fit does, and are raised
  # again with the caller's names (.relia_with_user_names()), as in
  # .omega_fit_cfa().
  user_named <- function(expr) {
    .relia_with_user_names(expr, c(y_names, a_names), c(items, aux))
  }

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

  # Three items under the congeneric model have a solution only when the
  # product of their covariances is positive (.triad_no_solution()); as
  # in .omega_fit_cfa(), the check is made before fitting, here on the
  # item block of the saturated (EM) estimate over the items and the
  # auxiliaries, and again on the fit's own sample statistics below if
  # that saturated fit fails.
  triad_checked <- FALSE
  if (!equal_loadings && J == 3L && is.null(h1_diagnostic)) {
    S_check <- if (!is.null(h1_cov)) {
      h1_cov
    } else {
      fc <- user_named(.fiml_cov(d, items = y_names, aux = a_names,
                                 se = "none"))
      if (isTRUE(fc$converged)) fc$S else NULL
    }
    if (!is.null(S_check)) {
      triad_checked <- TRUE
      if (.triad_no_solution(S_check)) {
        return(.omega_fit_none(J, items, no_solution = TRUE,
                               triad = .named_cov2cor(S_check, items)))
      }
    }
  }

  # check.post = FALSE: lavaan's post-estimation check warns that the
  # residual covariance matrix (theta) is not positive definite whenever
  # an auxiliary correlates substantially with the factor, because the
  # saturated correlates parameterization forces the entire item-aux
  # association through the residual covariances (the aux never touches
  # the factor). That is the expected signature of the specification,
  # not a defect of the data or the fit; the joint implied covariance
  # matrix is still proper. The check that matters for reliability, a
  # negative item error variance (Heywood case), is performed explicitly
  # below and travels back on the `improper` flag, so lavaan's note that
  # such a solution has run away is muffled as in .omega_fit_cfa().
  fit_args <- list(model = model_syntax, data = d, missing = "ml",
                   estimator = estimator, se = se, meanstructure = TRUE,
                   int.ov.free = TRUE, check.post = FALSE)
  fit <- try(user_named(withCallingHandlers(
    do.call(lavaan::lavaan, fit_args),
    warning = function(w) {
      if (.is_lavaan_improper_notice(w)) invokeRestart("muffleWarning")
    })), silent = TRUE)

  # lavaan's error, with the caller's names (.relia_lavaan_names_back()).
  if (inherits(fit, "try-error")) {
    return(.omega_fit_none(J, items, error = .relia_lavaan_names_back(
      conditionMessage(attr(fit, "condition")),
      c(y_names, a_names), c(items, aux))))
  }
  if (isFALSE(lavaan::lavInspect(fit, "converged"))) {
    return(.omega_fit_none(J, items))
  }
  # A fit whose likelihood has no maximum, as the CFA engine finds it for
  # the fits of .omega_fit_cfa() (.lavaan_no_maximum()).
  no_max <- .lavaan_no_maximum(fit)
  if (!is.null(no_max)) {
    no_max$syntax <- model_syntax
    no_max$lavaan_call <- .no_max_call_text("lavaan::lavaan", fit_args)
    return(.omega_fit_no_maximum(J, items, .no_max_renamed(
      no_max, c(y_names, a_names), c(items, aux))))
  }

  if (!equal_loadings && J == 3L && !triad_checked && is.null(h1_diagnostic)) {
    S_h1 <- .fit_item_cov(fit, y_names)
    if (!is.null(S_h1) && .triad_no_solution(S_h1)) {
      return(.omega_fit_none(J, items, no_solution = TRUE,
                             triad = .named_cov2cor(S_h1, items)))
    }
  }

  pe <- lavaan::parameterEstimates(fit)
  om_row   <- pe[pe$lhs == "omega" & pe$op == ":=", , drop = FALSE]
  ld_rows  <- pe[pe$op == "=~", "est"]
  err_rows <- pe[pe$op == "~~" & pe$lhs == pe$rhs &
                 pe$lhs %in% y_names, "est"]
  se_omega <- .lavaan_defined_se(fit, om_row$se)

  list(converged = TRUE,
       h1_diagnostic = h1_diagnostic,
       improper = any(err_rows < 0),
       no_solution = FALSE,
       triad = NULL,
       omega = om_row$est[1],
       se_omega = se_omega,
       loadings = ld_rows,
       errors   = err_rows,
       items = items,
       lavaan_items = y_names,
       J = J,
       fit = fit,
       coef_name = "omega")
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


# .stop_never_observed_together() refuses a pair of variables never
# observed together before a fit by full information maximum likelihood
# or two-stage estimation, and .stop_single_valued_together() a pair
# observed together only where one of them takes a single value. They
# live in R/cov_fiml.R, since cov_fiml() gives the same refusals in the
# same words; .omega_fit_cfa() and .fiml_cov() pass them the method's
# name (.relia_incomplete_label()).


# .fiml_cov(data, items, aux = NULL, se = "none",
#           fit = !identical(se, "none"), missing = "fiml")
#
# Full information maximum likelihood estimate of the covariance matrix
# among the items, for the analytic (closed form) reliability path under
# missing = "fiml" (and "two_stage", whose stage one is this estimate).
# The estimate is the item block of the maximum likelihood means and
# covariances of the saturated model over the items plus any auxiliary
# variables (every mean, variance, and covariance free), the estimates the
# EM algorithm converges to (Dempster, Laird, & Rubin, 1977), computed by
# lavaan's EM through .cov_fiml_em() in R/cov_fiml.R, the engine behind
# cov_fiml(). They are the unrestricted (h1) moments lavaan computes for
# every FIML fit to the same variables, lavInspect(fit, "sampstat.h1").
# Including the auxiliaries in the saturated model is what lets their
# information flow into the item block (Graham, 2003). A case with every
# variable missing carries no information and is set aside here, as
# lavaan sets it aside.
#
# Through DMAR 1.0.1 the estimate came from fitting the saturated model
# with lavaan (missing = "ml") from lavaan's default starting values,
# which maximizes the same likelihood directly over its p(p + 3)/2
# parameters. That is slow, and on a flat likelihood it stops short of
# the maximum: at the 26 Holzinger and Swineford tests (the Pasteur
# students missing the last two) it took 2484 iterations, and left
# covariances up to 0.0105 from the maximum likelihood values with the
# log likelihood 3e-7 below its maximum.
#
# When `fit` is TRUE the saturated model is also fitted, with every
# variance and covariance labeled and coefficient alpha of the item block
# defined as a `:=` parameter,
#
#     alpha := (J/(J-1)) * (1 - trace(Sigma_items) / sum(Sigma_items)),
#
# so that with se = "standard" lavaan's delta method gives a standard
# error for the analytic alpha that comes from the FIML information
# matrix (and therefore reflects the actual missingness pattern, unlike a
# complete-data closed form evaluated at the used N). The fit starts at
# the EM estimates with the optimizer switched off (optim.method =
# "none"), so its estimates are the EM estimates exactly and the
# information is evaluated there; the defined parameter then equals the
# alpha computed from S, and the estimate does not depend on whether a
# standard error was asked for. The fit also serves a Wald interval's
# wald_df resolution and the realized fraction of missing information of
# ss_aipe_reliability_sensitivity() (.aipe_realized_fmi()).
#
# Inputs:
#   data  : data frame containing the item (and auxiliary) columns.
#   items : character vector of item column names.
#   aux   : optional character vector of auxiliary column names.
#   se    : "none" (default; point estimates only, the fast path the
#           bootstrap uses) or "standard" (delta method SE for alpha).
#   fit   : whether to return the lavaan fit of the saturated model;
#           TRUE whenever a standard error is asked for.
#   missing : "fiml" (default) or "two_stage", the method the estimate
#           serves, which names it in the refusal of a pair of variables
#           never observed together (.stop_never_observed_together()).
#
# Returns:
#   A list with elements:
#     converged : logical; FALSE when the EM stopped at its iteration
#                 limit or failed, or the fit failed.
#     S         : (J x J) FIML covariance matrix of the items, with the
#                 original item names (NULL when converged is FALSE).
#     N         : number of cases used (those with an observed item or
#                 auxiliary variable).
#     alpha     : coefficient alpha of S.
#     alpha_se  : its delta method SE (NA_real_ when se = "none").
#     fit       : the lavaan fit of the saturated model (NULL when `fit`
#                 is FALSE or the fit failed).
#     coef_name : "alpha", the name of the defined parameter, which a
#                 Wald interval's wald_df resolution looks up in the fit.
.fiml_cov <- function(data, items, aux = NULL, se = "none",
                      fit = !identical(se, "none"), missing = "fiml") {
  if (!requireNamespace("lavaan", quietly = TRUE)) {
    .dmar_refuse(c(
      "missing = \"fiml\" requires the 'lavaan' package. ",
      "Install with install.packages(\"lavaan\")."),
      "package_required", package = "lavaan")
  }
  data <- as.data.frame(data)
  J <- length(items)
  K <- length(aux)
  y_names <- paste0("y", seq_len(J))
  # paste0() pads a zero-length argument to "" by default (recycle0 is
  # FALSE), so the no-auxiliary case must be handled explicitly.
  a_names <- if (K > 0L) paste0("a", seq_len(K)) else character(0)
  d <- data[, c(items, aux), drop = FALSE]
  .stop_never_observed_together(d, .relia_incomplete_label(missing))
  .stop_single_valued_together(d, .relia_incomplete_label(missing))
  colnames(d) <- c(y_names, a_names)
  vars <- c(y_names, a_names)
  p <- length(vars)
  d <- d[rowSums(!is.na(d)) > 0L, , drop = FALSE]
  failed <- list(converged = FALSE, S = NULL, N = nrow(d),
                 alpha = NA_real_, alpha_se = NA_real_, fit = NULL,
                 coef_name = "alpha")
  if (nrow(d) < 2L) return(failed)

  # The warnings lavaan raises about the data (two variables perfectly
  # correlated, variances of very different sizes) reach the user from the
  # EM, with the caller's names: the EM names the columns v1, v2, ...
  # (.cov_fiml_em(), and V1, V2, ... in a message about a matrix), and the
  # fit below y1, ..., yJ and a1, ..., aK, and each warning is raised again
  # with the names of `items` and `aux` (.relia_warn_named()); through
  # DMAR 1.0.1 and the first versions of 1.1.0 they reached the user as
  # lavaan gave them. The fit reads the same data and would raise them a
  # second time, under its own names, so it muffles a warning whose
  # renamed text the EM already gave.
  internal <- c(paste0("v", seq_len(p)), paste0("V", seq_len(p)), vars)
  callers <- rep(c(items, aux), 3L)
  seen <- new.env(parent = emptyenv())
  seen$messages <- character(0)
  em <- tryCatch(
    withCallingHandlers(.cov_fiml_em(d), warning = function(w) {
      named <- .relia_lavaan_names_back(conditionMessage(w), internal,
                                        callers)
      seen$messages <- c(seen$messages, named)
      .relia_warn_named(w, named)
    }),
    error = function(e) NULL)
  if (is.null(em) || !isTRUE(em$converged) || anyNA(em$S)) return(failed)
  S_items <- em$S[y_names, y_names, drop = FALSE]
  dimnames(S_items) <- list(items, items)
  out <- list(converged = TRUE,
              S = S_items,
              N = nrow(d),
              alpha = .alpha_from_S(S_items, em = TRUE),
              alpha_se = NA_real_,
              fit = NULL,
              coef_name = "alpha")
  if (!fit) return(out)

  # Saturated covariance structure with a label on every parameter; the
  # item-block labels feed the alpha definition.
  lab <- function(i, j) paste0("s", min(i, j), "_", max(i, j))
  ij <- which(upper.tri(diag(p), diag = TRUE), arr.ind = TRUE)
  ij <- ij[order(ij[, 1L], ij[, 2L]), , drop = FALSE]
  lines <- paste0(vars[ij[, 1L]], " ~~ ",
                  vapply(seq_len(nrow(ij)), function(k) {
                    lab(ij[k, 1L], ij[k, 2L])
                  }, character(1)),
                  "*", vars[ij[, 2L]])
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
  start <- data.frame(
    lhs = c(vars[ij[, 1L]], vars),
    op  = c(rep("~~", nrow(ij)), rep("~1", p)),
    rhs = c(vars[ij[, 2L]], rep("", p)),
    est = c(em$S[ij], em$mu),
    stringsAsFactors = FALSE)

  # The gradient check is left off: the estimates are the EM's by
  # construction, and the EM met its own convergence criterion.
  sat <- try(withCallingHandlers(
    lavaan::lavaan(model_syntax, data = d,
                   missing = "ml",
                   estimator = "ML",
                   se = se,
                   meanstructure = TRUE,
                   int.ov.free = TRUE,
                   start = start,
                   optim.method = "none",
                   check.gradient = FALSE),
    warning = function(w) {
      named <- .relia_lavaan_names_back(conditionMessage(w), internal,
                                        callers)
      if (named %in% seen$messages) invokeRestart("muffleWarning")
      .relia_warn_named(w, named)
    }), silent = TRUE)
  if (inherits(sat, "try-error") ||
      isFALSE(lavaan::lavInspect(sat, "converged"))) {
    return(failed)
  }

  pe <- lavaan::parameterEstimates(sat)
  al_row <- pe[pe$lhs == "alpha" & pe$op == ":=", , drop = FALSE]
  out$alpha_se <- if (identical(se, "none")) {
    NA_real_
  } else {
    .lavaan_defined_se(sat, al_row$se)
  }
  out$fit <- sat
  out
}


# .lavaan_defined_se(fit, se)
#
# The standard error lavaan reports for a defined (`:=`) parameter, or NA
# when lavaan computed no covariance matrix of the estimates. When the
# information matrix cannot be inverted, lavaan warns that it could not
# compute standard errors, leaves the covariance matrix empty, gives every
# free parameter a standard error of NA, and every row it does not
# estimate a standard error of zero (lav_model_vcov_se() sets the rows
# with free == 0 to zero), a defined parameter among them. Through the
# first versions of 1.1.0 (and DMAR 1.0.1) that zero was read as the
# standard error of alpha or omega, and the Wald interval came back with
# no width: on five pupils of three holzinger_swineford tests with a
# fourth as an auxiliary variable, the analytic alpha under two-stage
# estimation was 0.418 with the interval [0.418, 0.418]. A covariance
# matrix that was computed but is singular is not a failure: the
# equal-loadings fits have one by construction, through their equality
# constraints, and their standard errors are valid.
#
# Inputs:
#   fit : a lavaan fit.
#   se  : the `se` column of the defined parameter's row of
#         lavaan::parameterEstimates(fit) (NULL when there is none).
# Returns:
#   A number, or NA_real_.
.lavaan_defined_se <- function(fit, se) {
  if (is.null(se) || !length(se) || !is.finite(se[1L])) return(NA_real_)
  v <- tryCatch(lavaan::lavInspect(fit, "vcov"), error = function(e) NULL)
  if (is.null(v) || !length(v)) NA_real_ else se[1L]
}


# .omega_observed_from_data(data, estimator = "ML")
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
# The loadings come from a maximum likelihood fit. Robust ML (MLR) has
# the same estimates, to the last bit (it changes only the standard
# errors and the test statistic, which this estimator does not use), so
# through DMAR 1.0.1, when this helper fit by MLR, every bootstrap
# replication of robust omega paid for a scaled test statistic it never
# read: 1.3, 10, and 37 times the time of the ML fit at 9, 24, and 40
# items. The ML fit also serves the residual correlations the result
# reports, so one fit now does both (see .omega_observed_parts()).
#
# Inputs:
#   data      : raw items (N rows by J columns).
#   estimator : lavaan estimator for the single-factor fit; "ML" (the
#               default since DMAR 1.1.0; "MLR" through 1.0.1, with the
#               same estimate).
#
# Returns:
#   A single numeric value in the form a bootstrap point function
#   returns (see .omega_point_value()): the estimate; the estimate with
#   a `dmar_reason` of "improper" when the solution is improper; NA_real_
#   with that reason when the solution leaves no coefficient; NA_real_
#   if the CFA fit did not converge.
.omega_observed_from_data <- function(data, estimator = "ML") {
  p <- .omega_observed_parts(data, estimator = estimator)
  .omega_point_value(p$fit, p$estimate)
}


# .omega_observed_parts(data, estimator = "ML")
#
# The single-factor fit and robust omega computed from it, for the point
# estimate, which needs the fit as well (its verdict and its residual
# correlations). .omega_observed_from_data() is the bootstrap's view of
# the same computation.
#
# Returns:
#   A list with `estimate` (NA_real_ when the fit did not converge) and
#   `fit` (the .omega_fit_cfa() list).
.omega_observed_parts <- function(data, estimator = "ML") {
  fit <- .omega_fit_cfa(data = data, equal_loadings = FALSE,
                        estimator = estimator, se = "none",
                        missing = "listwise")
  if (!fit$converged) return(list(estimate = NA_real_, fit = fit))
  d_cc <- stats::na.omit(as.matrix(data))
  S <- stats::cov(d_cc)
  n <- nrow(d_cc)
  # Rescale the observed total variance to the fitted loadings' N-divisor
  # (ML) metric so the ratio matches MBESS::ci.reliability(type = "hierarchical").
  # Listwise deletion in the fit keeps numerator and denominator on the
  # same complete-case sample.
  list(estimate = (sum(fit$loadings))^2 / (sum(S) * (n - 1) / n),
       fit = fit)
}


# .omega_observed_fiml_parts(data, items, aux = NULL)
#
# Robust omega under full information maximum likelihood: the loadings
# in the numerator from the single-factor fit with missing = "ml" (with
# any auxiliary variables as saturated correlates), and the total
# variance in the denominator from the FIML estimate of the item
# covariance matrix (the saturated model over items plus auxiliaries).
# The FIML covariance matrix is already on the maximum likelihood (N
# divisor) metric, so numerator and denominator are internally
# consistent without the complete-case rescaling. The fit is by ML, for
# the reason .omega_observed_from_data() gives (it was MLR through DMAR
# 1.0.1, with the same estimate). The saturated fit comes first, so that
# with three items the single-factor fit's check for a solution
# (.triad_no_solution()) reads the same covariance matrix rather than
# fitting the saturated model a second time.
#
# Returns:
#   A list with `estimate` (NA_real_ when either fit failed) and `fit`
#   (the .omega_fit_cfa() list).
.omega_observed_fiml_parts <- function(data, items, aux = NULL) {
  fc <- .fiml_cov(data, items = items, aux = aux, se = "none")
  mf <- .omega_fit_cfa(data = data, equal_loadings = FALSE,
                       estimator = "ML", se = "none",
                       missing = "ml", aux = aux,
                       h1_cov = if (isTRUE(fc$converged)) fc$S else NULL)
  if (!mf$converged || !isTRUE(fc$converged)) {
    return(list(estimate = NA_real_, fit = mf))
  }
  list(estimate = (sum(mf$loadings))^2 / sum(fc$S), fit = mf)
}


# .omega_observed_ts_parts(data, items, aux = NULL, estimates_only = FALSE)
#
# Robust omega under two-stage estimation (missing = "two_stage"): the
# loadings in the numerator from the single-factor model fit to the stage
# one (EM) moments, and the total variance in the denominator from those
# same moments, the sum of the EM covariance matrix of the items (with any
# auxiliary variables in stage one, its item block). The stage one
# estimate is the saturated maximum likelihood estimate FIML's denominator
# uses too, and the fit carries it (lavInspect(fit, "sampstat")), so one
# fit gives both pieces. The two-stage and FIML coefficients differ only
# in the loadings. The fit is the normal theory one: the robust variant
# has the same estimates. A bootstrap replication passes estimates_only =
# TRUE, which fits the second stage alone to the EM moments
# (.omega_fit_cfa_ts_estimates()): the coefficient is the same, without
# the two-stage standard errors and test a replication never reads.
#
# Returns:
#   A list with `estimate` (NA_real_ when the fit failed) and `fit` (the
#   .omega_fit_cfa() list).
.omega_observed_ts_parts <- function(data, items, aux = NULL,
                                     estimates_only = FALSE) {
  mf <- if (estimates_only) {
    .omega_fit_cfa_ts_estimates(data, items = items, aux = aux)
  } else {
    .omega_fit_cfa(data = data[, c(items, aux), drop = FALSE],
                   equal_loadings = FALSE, estimator = "ML",
                   se = "none", missing = "two.stage", aux = aux)
  }
  if (!mf$converged) return(list(estimate = NA_real_, fit = mf))
  S <- .fit_item_cov(mf$fit, mf$lavaan_items)
  if (is.null(S)) return(list(estimate = NA_real_, fit = mf))
  list(estimate = (sum(mf$loadings))^2 / sum(S), fit = mf)
}


# .omega_fit_cfa_ts_estimates(data, items, aux = NULL,
#                             equal_loadings = FALSE)
#
# The single-factor fit of a bootstrap replication under two-stage
# estimation, for its estimates alone. lavaan's two-stage fit computes the
# two-stage standard errors and the scaled statistic whatever is asked
# (lav_options_set() replaces se = "none" and test = "none" with the
# two-stage ones, with a warning), and with them the observed and first
# order information of the saturated model and a fit of the baseline
# model, none of which a replication reads. The estimates are those of
# the complete data maximum likelihood fit to the stage one moments, so
# here stage one is lavaan's own EM (.ts_em_moments(), over the items and
# any auxiliary variables) and stage two is .omega_fit_cfa() on the item
# block of its covariance matrix, which keeps that function's check for
# three items with no solution, its verdict on an improper solution, and
# its convergence gate (and whatever else that function does with the
# items' names). The fit to a covariance matrix still computes lavaan's
# complete data statistic, which costs little. The two fits minimize the
# same discrepancy over the same covariance matrix, the free intercepts of
# the two-stage fit reproducing the EM means exactly, so they share their
# minimizer, and where the likelihood has a well defined minimum (proper
# solutions, and improper ones inside the parameter space) their
# estimates agree to rounding error (omega to 1e-15 on the tests' data,
# and on 120 cases with a fifth of the values missing). They do not
# follow the same path to it: the two-stage fit also carries the
# intercepts, and on a small, heavily incomplete resample, whose
# likelihood is nearly flat or has no proper maximum, the two optimizers
# stop at different points. The values then agree only to the
# optimizer's tolerance (6.6e-7 for omega on a resample of 45 cases with
# 30 percent of the values missing, up to 3.9e-5 on a kept improper
# replication), and a solution that runs away can be reported as not
# converged by one fit and as converged and improper by the other, in
# either direction, so the replication is counted as non-finite rather
# than improper or the reverse, and under improper = "keep" the set of
# replications behind the interval, and with it the interval, can
# differ. Neither fit is the better reference there (the two-stage fit
# has been seen to stop at a point that is not a minimum, with a
# discrepancy three times the fast fit's). A replication reclassified so
# is, under either verdict, one the bootstrap counts (as set aside, or as
# improper and kept), so the counts, and the caution that reads them,
# include it, though one such replication among many need not bring the
# caution on. The cost of the
# two-stage standard errors grows
# with the number of moments: on 20 items and 400 cases a replication
# here took about a third of the time of the two-stage fit (0.11 s against
# 0.31 s), while on 9 items, where lavaan's fixed cost per call dominates,
# the two took about the same time.
# The EM's warnings reach the caller as the two-stage fit's did, and the
# second fit's repetition of one of them is muffled.
#
# A pair of items or auxiliary variables that no row of the replication
# observes together is refused first, by name, with the words the
# two-stage fit uses (.stop_never_observed_together()): the second stage
# reads a covariance matrix, so .omega_fit_cfa()'s own check of the data
# is never reached, and lavaan's EM only warns and returns an arbitrary
# value for that covariance. The refusal is an error the bootstrap engine
# counts and quotes, as it counts the same resample under FIML. Through
# the first versions of 1.1.0 such a replication (a resample, or a BCa
# jackknife deletion, that had lost the last case observing a pair) was
# kept as a proper one, built on that arbitrary value. A pair the
# replication observes together only where one of them takes a single
# value is refused next, for the same reason
# (.stop_single_valued_together()).
#
# Inputs:
#   data           : the replication's rows, items and auxiliary columns.
#   items          : the item columns.
#   aux            : the auxiliary columns, or NULL.
#   equal_loadings : as in .omega_fit_cfa().
#
# Returns:
#   The .omega_fit_cfa() list (the .omega_fit_none() list when stage one
#   fails); an error when a pair is never observed together.
.omega_fit_cfa_ts_estimates <- function(data, items, aux = NULL,
                                        equal_loadings = FALSE) {
  data <- as.data.frame(data)
  .stop_never_observed_together(data[, c(items, aux), drop = FALSE],
                                .relia_incomplete_label("two_stage"))
  .stop_single_valued_together(data[, c(items, aux), drop = FALSE],
                               .relia_incomplete_label("two_stage"))
  J <- length(items)
  # The EM names the columns v1, v2, ... (.ts_em_moments()); its warnings
  # are raised again with the caller's names (.relia_warn_named()), and
  # the second fit, which gives the caller's names itself
  # (.omega_fit_cfa()), muffles a warning the EM already gave.
  internal <- paste0("v", seq_len(J + length(aux)))
  seen <- new.env(parent = emptyenv())
  seen$messages <- character(0)
  em <- tryCatch(
    withCallingHandlers(
      .ts_em_moments(data[, c(items, aux), drop = FALSE], quiet = FALSE),
      warning = function(w) {
        named <- .relia_lavaan_names_back(conditionMessage(w), internal,
                                          c(items, aux))
        seen$messages <- c(seen$messages, named)
        .relia_warn_named(w, named)
      }),
    error = function(e) NULL)
  # The cases the first stage read, a row observed on an auxiliary variable
  # alone among them, as in the full fit (.ts_all_rows()); the estimates do
  # not depend on the count, which only undoes lavaan's rescaling below.
  n <- sum(rowSums(!is.na(data[, c(items, aux), drop = FALSE])) > 0L)
  # A stage one that could not be estimated leaves no fit. The refusals
  # above leave no data on which the EM fails, so this is a guard (the
  # tests reach it by simulation); n < 2 cannot occur after them, and is
  # kept for the rescaling by n / (n - 1) below.
  if (is.null(em) || anyNA(em$cov) || n < 2L) {
    return(.omega_fit_none(J, items))
  }
  # The EM estimate has the N divisor, and lavaan multiplies a supplied
  # covariance matrix by (N - 1)/N, so it is passed multiplied by
  # N/(N - 1): the fit then reads the EM estimate itself.
  S <- em$cov[seq_len(J), seq_len(J), drop = FALSE] * n / (n - 1)
  dimnames(S) <- list(items, items)
  withCallingHandlers(
    .omega_fit_cfa(S = S, N = n, equal_loadings = equal_loadings,
                   estimator = "ML", se = "none", missing = "listwise"),
    warning = function(w) {
      if (conditionMessage(w) %in% seen$messages) {
        invokeRestart("muffleWarning")
      }
    })
}


# .relia_lavaan_missing(missing, robust = FALSE)
#
# lavaan's `missing` value for a reliability family `missing` option:
# "listwise" for listwise deletion, "ml" for FIML, and for two-stage
# estimation "two.stage" (normal theory standard errors) or, when `robust`
# (the "mlr" intervals), "robust.two.stage" (standard errors robust to
# nonnormality; Yuan & Bentler, 2000; Savalei & Falk, 2014). The point
# estimates of the two two-stage variants are the same.
.relia_lavaan_missing <- function(missing, robust = FALSE) {
  switch(missing,
         fiml = "ml",
         two_stage = if (robust) "robust.two.stage" else "two.stage",
         "listwise")
}


# .relia_convergence_note(missing)
#
# The sentence a reliability family message about a fit that did not
# converge adds for an incomplete data method ("" under listwise
# deletion).
.relia_convergence_note <- function(missing) {
  if (identical(missing, "fiml")) {
    paste0(" Full information maximum likelihood fails to ",
           "converge more often than complete-case estimation, ",
           "especially with small samples, many auxiliary ",
           "variables, or sparse missingness patterns.")
  } else if (identical(missing, "two_stage")) {
    paste0(" Two-stage estimation can fail at either stage, more often ",
           "with small samples, many auxiliary variables, or sparse ",
           "missingness patterns.")
  } else {
    ""
  }
}


# .relia_fit_failure(fit, model = "single-factor CFA", estimator = NULL,
#                    missing = "listwise")
#
# The refusal when a single-factor reliability fit (the .omega_fit_cfa()
# list) returned no solution: lavaan's own message when the fit stopped
# with an error (which the fit records with the caller's names for the
# variables, .relia_lavaan_names_back()), and otherwise that the fit did
# not converge, with the note for an incomplete data method
# (.relia_convergence_note()).
# Through DMAR 1.0.1 a fit that stopped with an error was reported as one
# that did not converge, which sent the user looking for a problem with
# the sample (on column names lavaan could not parse, for instance). An
# error that is not lavaan's (fit_error_lavaan FALSE, one of fa_confirmatory()'s own
# checks) is quoted without naming lavaan.
.relia_fit_failure <- function(fit, model = "single-factor CFA",
                               estimator = NULL, missing = "listwise") {
  with_est <- if (is.null(estimator)) {
    ""
  } else {
    paste0(" with estimator = '", estimator, "'")
  }
  error <- fit[["fit_error"]]
  if (!is.null(error)) {
    return(paste0("The ", model, " model could not be fit", with_est,
                  if (isFALSE(fit[["fit_error_lavaan"]])) {
                    ": "
                  } else {
                    "; lavaan stopped with this error: "
                  },
                  trimws(gsub("[[:space:]]+", " ", error))))
  }
  paste0("The ", model, " model did not converge", with_est, ".",
         .relia_convergence_note(missing))
}


# .relia_incomplete_label(missing)
#
# The words the reliability family's messages use for an incomplete data
# method.
.relia_incomplete_label <- function(missing) {
  if (identical(missing, "two_stage")) {
    "two-stage maximum likelihood"
  } else {
    "full information maximum likelihood"
  }
}


# .omega_fit_verdict(fit)
#
# What a single-factor reliability fit (the .omega_fit_cfa() list) can
# report, in the terms of McDonald's (2004) account of improper
# solutions, which follows van Driel (1978):
#
#   "proper"         : every error variance estimate is nonnegative.
#   "improper"       : some error variance is negative (a Heywood case),
#                      but the error variances sum to a nonnegative
#                      number, so the model still implies a nonnegative
#                      error variance for the composite and a model
#                      implied omega between zero and one. Such a
#                      solution can reflect sampling error around a
#                      proper structure (McDonald's Type I, a negative
#                      variance near zero) or a misspecified model (his
#                      Type II), which one sample cannot always tell
#                      apart, so the coefficient is reported with a
#                      warning, as through DMAR 1.0.1.
#   "negative_error" : the error variances sum to less than zero, so the
#                      model implies a negative error variance for the
#                      composite and a model implied omega above one,
#                      which no reliability coefficient can take. No
#                      sampling fluctuation near a zero residual
#                      variance does this unless every item is nearly
#                      error free; it is the signature of a solution
#                      that has run away along a ridge of the likelihood
#                      (a misspecified or empirically underidentified
#                      model, McDonald's Types II and III), and there is
#                      no coefficient to report.
#   "no_solution"    : three items whose covariances admit no
#                      single-factor solution (.triad_no_solution()).
#   "not_converged"  : the fit failed or did not converge.
.omega_fit_verdict <- function(fit) {
  if (isTRUE(fit$no_solution)) return("no_solution")
  if (!isTRUE(fit$converged)) return("not_converged")
  if (!isTRUE(fit$improper)) return("proper")
  if (sum(fit$errors) < 0) "negative_error" else "improper"
}


# .omega_point_value(fit, value = fit$omega)
#
# A single-factor coefficient in the form a bootstrap point function
# returns to .bootstrap_ci(): `value` as it is for a proper solution;
# `value` carrying a `dmar_reason` of "improper" for an improper solution
# that still yields a coefficient (kept or set aside by the engine's
# `improper` rule); NA_real_ with that reason when the solution yields no
# coefficient because its error variances sum below zero
# ("negative_error"); NA_real_ with the reason "no_solution" for three
# items with no single-factor solution (both set aside under either rule
# and counted as improper, the second also on its own; see
# .boot_classify()); and NA_real_ when the fit did not converge or the
# coefficient could not be computed. A fit whose likelihood has no maximum
# (`no_maximum`, .omega_fit_no_maximum()) gives no value: the replication
# stops with .no_max_replication()'s short refusal, and the engine counts
# it among those that raised an error, naming the variable whose error
# variance kept falling.
.omega_point_value <- function(fit, value = fit$omega) {
  # A replication whose likelihood has no maximum stops, and the engine
  # counts it among the replications that raised an error, with the short
  # message that names the variable (.no_max_replication()).
  if (!is.null(fit$no_maximum)) .no_max_replication(fit$no_maximum)
  verdict <- .omega_fit_verdict(fit)
  if (verdict == "no_solution") {
    return(structure(NA_real_, dmar_reason = "no_solution"))
  }
  if (verdict == "negative_error") {
    return(structure(NA_real_, dmar_reason = "improper"))
  }
  if (length(value) != 1L || !is.finite(value)) return(NA_real_)
  if (verdict == "improper") {
    return(structure(value, dmar_reason = "improper"))
  }
  value
}


# .omega_point_gate(fit, improper = "keep", bootstrap = FALSE,
#                   model = "single-factor", fn = "reliability_omega")
#
# What the point estimate of reliability_omega() and of the model
# implied reliability_alpha() does with its fit's verdict
# (.omega_fit_verdict()), after the one refusal that comes first: a fit
# whose likelihood has no maximum (`no_maximum`, .omega_fit_no_maximum())
# stops with the refusal of class dmar_no_maximum under the name of the
# calling function, `fn`, with the lavaan code of the single-factor fit
# (.no_max_refuse()).
#
#   "no_solution", "negative_error" : stop, with a message that says why
#       the solution yields no coefficient and what usually causes it.
#       Through DMAR 1.0.1 such a fit was returned with a warning about
#       a Heywood case; on three items whose correlations had a negative
#       product that was an omega of 471 from a loading of 43.6 and a
#       residual variance of -1899, the point at which lavaan's optimizer
#       stopped on a likelihood that has no maximum. McDonald (2004)
#       separates improper estimates that reflect sampling error around a
#       proper structure, which may be reported, from those that reflect
#       a misspecified or underidentified model, which call for a
#       different model rather than a number; a solution that implies a
#       negative error variance for the composite is of the second kind.
#   "improper" : under improper = "drop" with a bootstrap interval,
#       stop: the interval would rest on proper replications only and
#       would not describe this estimate. Otherwise warn, as through
#       1.0.1, and continue.
#   "not_converged", "proper" : nothing; the caller handles a fit that
#       did not converge with its own message.
#
# Returns the verdict invisibly.
.omega_point_gate <- function(fit, improper = "keep", bootstrap = FALSE,
                              model = "single-factor",
                              fn = "reliability_omega") {
  if (!is.null(fit$no_maximum)) .no_max_refuse(fn, fit$no_maximum)
  verdict <- .omega_fit_verdict(fit)
  if (verdict == "no_solution") {
    .dmar_refuse(.no_solution_message(fit$triad), "computation_failed")
  }
  if (verdict == "negative_error") {
    .dmar_refuse(.negative_error_message(fit, model), "computation_failed")
  }
  if (verdict == "improper") {
    neg <- fit$items[fit$errors < 0]
    if (identical(improper, "drop") && isTRUE(bootstrap)) {
      .dmar_refuse_argument(c(
        "The point estimate comes from an improper solution (a ",
        "Heywood case: the error variance of ",
        paste(neg, collapse = " and "), " is negative), and ",
        "improper = \"drop\" sets aside every bootstrap replication ",
        "with such a solution, so the interval would rest on proper ",
        "solutions only and would not describe this estimate. Use ",
        "improper = \"keep\" (the default), under which the ",
        "estimate is returned with a warning, or reconsider the ",
        model, " model for these items."), arg = NULL)
    }
    .dmar_warn(c(
      "Improper solution (Heywood case): ",
      .improper_items_phrase(fit), ", so the model reproduces an ",
      "impossible variance. The error variances still sum to a ",
      "nonnegative number, so the coefficient remains a proportion ",
      "and is returned, but interpret it with care: a negative ",
      "error variance near zero is what sampling error around a ",
      "proper structure produces, and one as large as the item's ",
      "own variance points to a misspecified model (McDonald, ",
      "2004). Consider whether the ", model, " model is ",
      "appropriate for these items."), "heywood_warning")
  }
  invisible(verdict)
}


# .improper_items_phrase(fit)
#
# The phrase naming each item with a negative error variance in an
# improper single-factor fit, with its observed variance when the fit
# can supply it, for the warning .omega_point_gate() gives: "the error
# variance of y1 is -0.85 (its observed variance is 0.53)". The numbers
# are formatted for the message only.
.improper_items_phrase <- function(fit) {
  neg <- which(fit$errors < 0)
  fmt <- function(v) trimws(formatC(v, format = "fg", digits = 3))
  obs <- if (is.null(fit$h1_diagnostic) && !is.null(fit$fit) && !is.null(fit$lavaan_items)) {
    S <- .fit_item_cov(fit$fit, fit$lavaan_items)
    if (is.null(S)) NULL else diag(S)[neg]
  }
  one <- length(neg) == 1L
  paste0("the error ", if (one) "variance" else "variances", " of ",
         paste(fit$items[neg], collapse = " and "),
         if (one) " is " else " are ",
         paste(fmt(fit$errors[neg]), collapse = " and "),
         if (!is.null(obs)) {
           paste0(" (", if (one) "its observed variance is " else
                    "their observed variances are ",
                  paste(fmt(obs), collapse = " and "), ")")
         } else {
           ""
         })
}


# .no_solution_message(triad)
#
# The refusal for three items whose correlations have a product that is
# not positive (.triad_no_solution()), quoting the three correlations.
.no_solution_message <- function(triad) {
  items <- rownames(triad)
  pairs <- rbind(c(1L, 2L), c(1L, 3L), c(2L, 3L))
  r <- triad[pairs]
  quoted <- paste0("r(", items[pairs[, 1]], ", ", items[pairs[, 2]],
                   ") = ", formatC(r, format = "f", digits = 3),
                   collapse = ", ")
  why <- if (prod(r) < 0) {
    paste0("is negative (", quoted, "), but a single factor implies ",
           "that it equals the square of the product of the three ",
           "loadings, which cannot be negative. No real loadings ",
           "reproduce these correlations")
  } else {
    paste0("is zero (", quoted, "), which a single factor implies only ",
           "when a loading is zero, and then the other two loadings are ",
           "not determined")
  }
  paste0("The single-factor model has no solution for these three ",
         "items: the product of their correlations ", why, ", so the ",
         "maximum likelihood estimate does not exist, and the ",
         "coefficient at whatever point an optimizer stops would ",
         "describe the optimizer rather than the data. An item keyed in ",
         "the opposite direction to the others (reverse score it), or ",
         "items that do not share a common factor, are the usual causes.")
}


# .negative_error_message(fit, model = "single-factor")
#
# The refusal for an improper solution whose error variances sum to less
# than zero (.omega_fit_verdict()), naming each item with a negative
# error variance and the coefficient the solution implies. The numbers
# are formatted for the message only.
.negative_error_message <- function(fit, model = "single-factor") {
  neg <- which(fit$errors < 0)
  fmt <- function(v) trimws(formatC(v, format = "fg", digits = 4))
  which_items <- paste0(fit$items[neg], " (", fmt(fit$errors[neg]), ")",
                        collapse = " and ")
  paste0("The ", model, " solution is improper beyond use: the error ",
         "variance of ", which_items, " is negative, and the error ",
         "variances sum to ", fmt(sum(fit$errors)), ", so the model ",
         "implies a negative error variance for the composite and a ",
         "model implied omega of ", fmt(fit$omega), ", which no ",
         "reliability coefficient can take. Sampling error around a ",
         "proper structure gives a negative error variance near zero, ",
         "not this; an improper solution this far from the boundary ",
         "indicates a misspecified or empirically underidentified model ",
         "(McDonald, 2004). Two items that share a specific factor, an ",
         "item keyed in the opposite direction to the others, or items ",
         "that do not share a common factor are the usual causes; the ",
         "item correlations, or cfa_1() on these items, show which.")
}


# .omega_c_from_data(data, denominator = "observed")
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
# numerator. The denominator correlations depend on `denominator`:
#
#   "observed"      : the unrestricted polychoric correlations, each
#                     estimated from its pair's bivariate table (Kelley &
#                     Pornprasertmanit, 2016, Eq. 18), so structure the
#                     single-factor model omits stays in the total
#                     variance. This is the
#                     categorical counterpart of robust omega and the
#                     DMAR default.
#   "model_implied" : the latent-response correlations reproduced by the
#                     fitted single-factor model, as in Green and Yang's
#                     (2009) Eq. 21, whose appendix program takes "the
#                     reproduced (model-based) correlation matrix" as the
#                     input to the denominator.
#
# When the single-factor model reproduces the polychoric correlations
# the two coincide. Treating ordered-categorical items as continuous and
# computing coefficient omega would mis-model the link between latent
# ability and observed responses, especially when categories are few or
# threshold patterns differ across items.
#
# Inputs:
#   data        : a data frame or matrix of ordered-categorical item
#                 scores (integer-valued). Columns are coerced to
#                 ordered factors before fitting.
#   denominator : "observed" (default) or "model_implied".
#   missing     : "listwise" (default) or "pairwise", passed to
#                 .omega_c_parts() (see there).
#
# Returns:
#   A single numeric value, sample omega_C. NA_real_ if the categorical
#   CFA fit stopped with an error or did not converge; NA_real_ carrying
#   a `dmar_reason` attribute when the solution is improper ("improper":
#   a standardized loading above one, which no denominator can use) or
#   when the data do not determine some polychoric correlation
#   ("undetermined_polychoric"): under the observed denominator any
#   such correlation, which it needs, and under the model implied one
#   those that leave the loadings unidentified. .bootstrap_ci() reads
#   that attribute to count the replications it sets aside by cause.
.omega_c_from_data <- function(data, denominator = "observed",
                               missing = "listwise") {
  parts <- .omega_c_parts(data, missing = missing)
  # No fit: a replication whose fit lavaan stopped with an error, or which
  # did not converge, has no value, and is counted among the non-finite
  # ones, as the continuous fit's is (.omega_point_value()).
  if (!isTRUE(parts$converged)) return(NA_real_)
  # Correlations the data do not determine can leave the loadings
  # unidentified (.omega_c_identification()), and then categorical omega
  # under either denominator is whatever the start values led to; the
  # replication is set aside and counted with those that left a
  # correlation undetermined. The check comes before the one for an
  # improper solution, whose loadings are as arbitrary as the rest.
  # Through the first versions of 1.1.0 the model implied denominator
  # kept such a replication.
  if (isFALSE(parts$identified)) {
    return(structure(NA_real_, dmar_reason = "undetermined_polychoric"))
  }
  # An improper solution (a communality above one) has no categorical
  # omega under either denominator: the bivariate normal at a
  # correlation above one does not exist.
  if (length(parts$improper)) {
    return(structure(NA_real_, dmar_reason = "improper"))
  }
  # An observed denominator needs every unrestricted polychoric; a
  # bootstrap replication in which the data do not determine one is
  # set aside rather than filled with an arbitrary value.
  if (!identical(denominator, "model_implied") &&
      nrow(parts$undetermined) > 0L) {
    return(structure(NA_real_, dmar_reason = "undetermined_polychoric"))
  }
  rho_den <- if (identical(denominator, "model_implied")) {
    parts$rho_model
  } else {
    parts$rho_sat
  }
  .gy_omega_from_components(parts$thresh, parts$rho_num, rho_den)
}


# .omega_c_parts(data, missing = "listwise")
#
# The single-factor ordered-categorical fit behind categorical omega,
# returned as the pieces that both denominators and the local
# dependence summary draw on, so the point estimate and the residual
# correlations come from one fit.
#
# Inputs:
#   data    : a data frame or matrix of ordered-categorical item scores
#             (integer-valued).
#   missing : "listwise" (the default; lavaan's own default for ordered
#             indicators, so the fit is the one DMAR made through 1.0.1)
#             or "pairwise", passed to lavaan: each item's thresholds
#             come from its own observed values and each polychoric
#             correlation from the cases that observe both items
#             (missing = "pairwise" of reliability_omega_categorical()).
#             The caller removes the rows with no observed item first,
#             which lavaan would otherwise drop with a warning. Under
#             pairwise deletion lavaan's Gamma matrix accounts for the
#             number of cases behind each polychoric, so
#             .undetermined_polychorics() applies unchanged.
#
# Returns:
#   When there is no fit, a list whose `converged` is FALSE and whose
#   `fit_error`, when lavaan stopped with an error, is that error's
#   message (absent when the fit ran but did not converge), the shape of
#   .omega_fit_none() for the continuous fit, so the callers can quote
#   lavaan's message (.relia_fit_failure()) and keep it apart from a
#   convergence failure. Through the first versions of 1.1.0 (and DMAR
#   1.0.1) the helper returned NULL for both, and a lavaan error reached
#   the user as a model that did not converge. The callers also read a
#   NULL as no fit. Otherwise a list with
#     converged : TRUE.
#     thresh    : list of threshold vectors, one per item, in item order.
#     rho_num   : (J x J) model implied common-factor correlations
#                 (lambda lambda' scaled to the latent-response
#                 correlation metric; communalities on the diagonal).
#     rho_model : (J x J) latent-response correlations reproduced by the
#                 fitted single-factor model (unit diagonal).
#     rho_sat   : (J x J) unrestricted polychoric correlations.
#     undetermined : two-column character matrix of the item pairs whose
#                 polychoric the data do not determine (zero rows when
#                 none); see .undetermined_polychorics().
#     identified : TRUE when the correlations the data do determine
#                 identify the loadings (always, when none is
#                 undetermined), FALSE otherwise; the loadings, the
#                 model implied correlations, and categorical omega under
#                 either denominator are then arbitrary, and the callers
#                 report no value (.omega_c_identification()).
#     identification : the list .omega_c_identification() returned (the
#                 blocks and sides its message names).
#     improper  : the names of the items whose communality exceeds one
#                 (a negative residual variance, a Heywood case;
#                 character(0) when none). Categorical omega is not
#                 defined for such a solution, since its assembly
#                 evaluates the bivariate normal at the communality.
#     loadings  : the standardized loadings, named by item, which the
#                 message about an improper solution quotes.
#     df        : the degrees of freedom of the single-factor
#                 correlation structure, J (J - 3) / 2 (the thresholds
#                 are saturated), so a three-item fit is recognized as
#                 saturated.
#   Each matrix carries the item names of `data` as dimnames.
#
# The single-factor model for ordered items is not identified with two
# items (the one polychoric correlation determines only the product of
# the two loadings), so the helper refuses fewer than three, as a
# backstop behind the check in reliability_omega_categorical().
#
# Dependencies:
#   lavaan (single-factor categorical CFA). The bivariate normal
#   probabilities of the assembly are DMAR's own (.pbinorm()); through
#   version 1.0.1 they came from mvtnorm, which categorical omega
#   therefore also required.
.omega_c_parts <- function(data, missing = "listwise") {
  if (!requireNamespace("lavaan", quietly = TRUE)) {
    .dmar_refuse(c(
      "Categorical omega requires the 'lavaan' package. ",
      "Install with install.packages(\"lavaan\")."),
      "package_required", package = "lavaan")
  }

  data <- as.data.frame(data)
  J <- ncol(data)
  if (J < 3L) .stop_omega_c_two_items()
  item_names <- colnames(data)
  if (is.null(item_names)) item_names <- paste0("y", seq_len(J))
  for (j in seq_len(J)) data[[j]] <- ordered(data[[j]])
  var_names <- paste0("y", seq_len(J))
  colnames(data) <- var_names

  load_line <- paste(paste0("a", seq_len(J), "*", var_names),
                     collapse = " + ")
  model <- paste0("f1 =~ NA*", var_names[1], " + ", load_line,
                  "\nf1 ~~ 1*f1\n")
  # Only the estimates and the sample statistics are used, so neither
  # standard errors nor a test statistic is computed (test = "none"
  # leaves every estimate, the polychorics, and the Gamma matrix as they
  # are, takes less than half the time, and spares a bootstrap the
  # robust test statistic's warning on a resample whose information
  # matrix is singular). lavaan's post-fit check warns when an estimated
  # residual variance is negative; the `improper` element below reports
  # the same condition by item, and the callers act on it (an error
  # naming the item for the point estimate, a counted reason in the
  # bootstrap), so the warning would only repeat once per bootstrap
  # replication. Every other warning is raised again with the caller's
  # names for the items (.relia_with_user_names(); "correlation between
  # variables y2 and y1 is (nearly) 1.0" said y2 and y1 through the first
  # versions of 1.1.0), the item names only, as for the error below.
  fit <- tryCatch(.relia_with_user_names(withCallingHandlers(
    lavaan::cfa(model, data = data, se = "none", test = "none",
                ordered = var_names, missing = missing),
    warning = function(w) {
      if (.is_lavaan_improper_notice(w)) invokeRestart("muffleWarning")
    }), var_names, item_names), error = function(e) e)
  # A fit that lavaan stopped with an error carries the error's message
  # back as `fit_error` (not `error`, which `$` would match partially
  # against a longer name), apart from a fit that ran and did not
  # converge, with the items given the caller's names
  # (.relia_lavaan_names_back()). Only the item names are mapped: a1,
  # ..., aJ here are the loading labels of the syntax, not auxiliary
  # variables.
  if (inherits(fit, "error")) {
    return(list(converged = FALSE, fit_error = .relia_lavaan_names_back(
      conditionMessage(fit), var_names, item_names)))
  }
  if (isFALSE(lavaan::lavInspect(fit, "converged"))) {
    return(list(converged = FALSE))
  }

  # Extract loadings, factor variance (= 1), thresholds, and the
  # unrestricted polychoric correlations from the same fit.
  # Every matrix is indexed by item name, so the pieces stay paired item
  # by item whatever order lavaan stores its observed variables in.
  est <- lavaan::lavInspect(fit, "est")
  lambda <- est$lambda[var_names, , drop = FALSE]
  psi    <- est$psi
  truevar <- lambda %*% psi %*% t(lambda)
  thresh  <- .ordered_thresholds(fit)
  poly    <- .polychoric_from_fit(fit, var_names)
  sigma_hat <- lavaan::lavInspect(fit, "implied")$cov
  sigma_hat <- sigma_hat[var_names, var_names, drop = FALSE]
  inv_sd  <- 1 / sqrt(diag(sigma_hat))
  poly_r  <- diag(inv_sd) %*% truevar %*% diag(inv_sd)
  # The model implied latent-response correlations: under the delta
  # parameterization the implied variances are already one, and the
  # rescaling makes that explicit rather than assumed.
  rho_model <- diag(inv_sd) %*% sigma_hat %*% diag(inv_sd)
  # An improper solution: a communality above one, which under the delta
  # parameterization is a negative residual variance (both are checked,
  # so the theta parameterization is covered as well).
  theta_diag <- diag(as.matrix(est$theta))[var_names]
  improper <- item_names[diag(poly_r) > 1 |
                           (!is.na(theta_diag) & theta_diag < 0)]

  named <- function(m) {
    m <- unname(as.matrix(m))
    dimnames(m) <- list(item_names, item_names)
    m
  }
  undetermined <- .undetermined_polychorics(fit, var_names, item_names)
  identification <- .omega_c_identification(item_names, undetermined)
  list(converged = TRUE,
       thresh = unname(thresh[var_names]),
       rho_num = named(poly_r),
       rho_model = named(rho_model),
       rho_sat = named(poly),
       undetermined = undetermined,
       identified = identification$identified,
       identification = identification,
       improper = improper,
       loadings = stats::setNames(as.numeric(lambda[, 1] * inv_sd),
                                  item_names),
       df = J * (J - 3) / 2)
}


# .stop_omega_c_two_items()
#
# The refusal of a two-item categorical omega, shared by
# reliability_omega_categorical() and the backstop in .omega_c_parts()
# so the two cannot word it differently. It follows the wording of the
# same refusal in reliability_omega().
.stop_omega_c_two_items <- function() {
  .dmar_refuse_argument(c(
    "Categorical omega needs three or more items: the single-factor ",
    "model for ordered items is not identified with two (the one ",
    "polychoric correlation determines only the product of the two ",
    "loadings). For a two-item scale use reliability_alpha(), whose ",
    "tau-equivalent model is just identified."), arg = NULL)
}


# .improper_message(parts)
#
# The error text reliability_omega_categorical() gives for an improper
# single-factor solution: it names each item whose communality exceeds
# one with its standardized loading, so the user sees which item the
# solution could not accommodate. The loadings are formatted here for
# the message only; nothing stored is rounded.
#
# Inputs:
#   parts : the list .omega_c_parts() returned, with `improper` nonempty.
# Returns:
#   A single character string.
.improper_message <- function(parts) {
  items <- parts$improper
  load <- formatC(parts$loadings[items], format = "f", digits = 3)
  which <- if (length(items) == 1L) {
    paste0("the standardized loading of ", items, " is ", load,
           ", so its residual variance is negative")
  } else {
    paste0("the standardized loadings of ",
           paste(items, collapse = " and "), " are ",
           paste(load, collapse = " and "),
           ", so their residual variances are negative")
  }
  paste0("The single-factor solution is improper (a Heywood case): ",
         which, ", and categorical omega is not defined for it. ",
         "Collapsing sparse categories, a larger sample, or reconsidering ",
         "the item are the usual remedies.")
}


# .gy_omega_from_components(thresh, rho_num, rho_den, weights = NULL)
#
# The Green and Yang (2009) assembly shared by the single-factor
# categorical omega and the per-factor categorical omega of fa_confirmatory():
# given the item thresholds and two latent-response correlation
# matrices, sum the bivariate normal cell probabilities over all
# threshold pairs to obtain the covariance of the categorical sum
# score attributable to the common factor(s) (numerator, rho_num =
# model implied correlations) and in total (denominator, rho_den = the
# unrestricted polychoric correlations, or the latent-response
# correlations the fitted model reproduces), and return their ratio.
#
# With the categories of item j scored 0, 1, ..., C_j - 1 (the number of
# thresholds its latent response exceeds), the covariance of the scores
# of items j and j' at latent-response correlation r is
#   sigma_jj'(r) = sum_c sum_c' Phi2(t_jc, t_j'c'; r)
#                  - (sum_c Phi(t_jc)) (sum_c' Phi(t_j'c'))
# (Green & Yang, 2009, Eqs. 13--14), and the reliability of the
# composite sum_j w_j X_j is
#   sum_j sum_j' w_j w_j' sigma_jj'(rho_num[j, j'])
#     / sum_j sum_j' w_j w_j' sigma_jj'(rho_den[j, j'])
# (their Eq. 21 when every w_j is one; Kelley & Pornprasertmanit, 2016,
# Eqs. 17--18). On the diagonal rho_num holds the item's communality and
# rho_den one, where Phi2(t, t'; 1) = Phi(min(t, t')), so sigma_jj(1) is
# the variance of the item's score.
#
# Every cell (a threshold of item j with a threshold of item j') of
# every ordered item pair (j, j') is laid out in one vector, in the
# order the loops of version 1.0.1 visited them, and .pbinorm()
# evaluates the cells in one call per correlation matrix. Those loops
# made one mvtnorm::pmvnorm() call per cell; each cell now agrees with
# that call to 4e-16, and the coefficients differ only by the rounding
# of their sums: by at most 9e-16 in Green and Yang's ten population
# designs (tests/testthat/test-pbinorm.R pins them) and in six to 24
# Holzinger and Swineford tests cut into two to seven categories, under
# both denominators and in fa_confirmatory()'s per-factor omega. The pairs
# (j, j') and (j', j) are both evaluated, as before, although sigma_jj'
# is symmetric: computing each pair once would halve a cost that is
# already a small fraction of the factor model's fit, at the price of a
# further reordering of the sums. At 12 items with 7 categories the
# assembly takes about 0.002 s against 0.37 s for the loops, so a
# bootstrap replication (refit and assembly) takes 0.07 s rather than
# 0.43 s.
#
# Inputs:
#   thresh  : list of numeric threshold vectors, one per item, in item
#             order.
#   rho_num : (J x J) model implied latent-response correlation matrix
#             (common-factor part, including communalities on the
#             diagonal).
#   rho_den : (J x J) latent-response correlation matrix behind the
#             total variance (unit diagonal).
#   weights : NULL (the default) for unit weights, the unweighted sum
#             score whose reliability is categorical omega; or J finite
#             item weights, for the reliability of a weighted composite
#             of the category scores (no exported function passes
#             weights yet).
#
# Returns:
#   A single numeric value, the Green and Yang categorical omega.
.gy_omega_from_components <- function(thresh, rho_num, rho_den,
                                      weights = NULL) {
  J <- length(thresh)
  if (is.null(weights)) {
    weights <- rep(1, J)
  } else if (length(weights) != J || !all(is.finite(weights))) {
    .dmar_refuse_argument(c(
      "The weights of a categorical composite must be ", J,
      " finite numbers, one per item."), arg = NULL)
  }
  n_t <- lengths(thresh)
  j1 <- rep(seq_len(J), each = J)
  j2 <- rep(seq_len(J), times = J)
  # Cell q of a pair pairs threshold q %/% C2 + 1 of item j1 with
  # threshold q %% C2 + 1 of item j2 (C2 the number of thresholds of
  # item j2), so the second item's thresholds vary fastest.
  n_cells <- n_t[j1] * n_t[j2]
  pair <- rep.int(seq_along(j1), n_cells)
  q <- sequence(n_cells) - 1L
  tv <- unlist(thresh, use.names = FALSE)
  start <- cumsum(c(0L, n_t))[seq_len(J)]
  n2 <- n_t[j2][pair]
  t1 <- tv[start[j1][pair] + q %/% n2 + 1L]
  t2 <- tv[start[j2][pair] + q %% n2 + 1L]
  at <- cbind(j1, j2)
  # rowsum() keeps the pairs in the order they first appear, which is
  # the order of the pairs that have cells.
  has <- which(n_cells > 0L)
  joint_num <- joint_den <- numeric(length(j1))
  joint_num[has] <- rowsum(.pbinorm(t1, t2, rho_num[at][pair]), pair,
                           reorder = FALSE)[, 1L]
  joint_den[has] <- rowsum(.pbinorm(t1, t2, rho_den[at][pair]), pair,
                           reorder = FALSE)[, 1L]
  margin <- vapply(thresh, function(t) sum(stats::pnorm(t)), numeric(1L))
  product <- margin[j1] * margin[j2]
  ww <- weights[j1] * weights[j2]
  sum(ww * (joint_num - product)) / sum(ww * (joint_den - product))
}


# .omega_categorical_from_kfit(fit, items_f)
#
# Categorical omega for the sum score of one factor's items (the Green
# and Yang, 2009, assembly with the unrestricted polychoric denominator
# of Kelley and Pornprasertmanit, 2016, Eq. 18), computed from a fitted
# multi-factor categorical model (fa_confirmatory() with ordered items). The numerator correlations are
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
#   When the solution is improper for one of the factor's items (a
#   negative residual variance of its latent response, so a communality
#   above one), categorical omega is not defined, since the assembly
#   evaluates the bivariate normal distribution at the communality: the
#   value is then NA_real_ with a `dmar_reason` of "improper" and an
#   `improper` attribute naming the items with their communalities,
#   which fa_confirmatory() quotes in its warning (.kfit_improper_message()).
#   Through DMAR 1.0.1 such a solution stopped fa_confirmatory() with mvtnorm's
#   "'corr' is not a correlation matrix".
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
  # An improper solution, checked as .omega_c_parts() checks the
  # single-factor fit: a communality above one, or a negative residual
  # variance of the latent response (the theta parameterization fa_confirmatory()
  # uses estimates it directly).
  communality <- diag(rho_num)[items_f]
  theta_diag <- diag(as.matrix(est$theta))
  names(theta_diag) <- rownames(as.matrix(est$theta))
  theta_f <- theta_diag[items_f]
  bad <- items_f[(!is.na(communality) & communality > 1) |
                   (!is.na(theta_f) & theta_f < 0)]
  if (length(bad)) {
    return(structure(NA_real_, dmar_reason = "improper",
                     improper = communality[bad]))
  }
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


# .kfit_improper_message(factor, communality, ordered = TRUE, omega = TRUE)
#
# The warning fa_confirmatory() gives when the solution is improper for some of a
# factor's items (a negative error variance, or for an ordered item a
# negative residual variance of its latent response, so a communality
# above one), which leaves some of the factor's measurement properties
# undefined and their rows NA. It names the rows, each item with its
# communality, and why each property is undefined: the average variance
# extracted averages the communalities, and coefficient H, the largest
# reliability that a weighted sum of the items attains, is at least the
# largest communality (the item weighted alone attains it), so neither is
# a proportion of variance; categorical omega evaluates the bivariate
# normal distribution at each communality; and linear omega is NA only
# when the error variances sum below zero, and otherwise reported. It is
# the counterpart of .improper_message() for the single-factor fit. The
# communalities are formatted for the message only.
#
# Inputs:
#   factor      : the factor's name.
#   communality : named numeric vector, the communalities above one (of
#                 the items whose error variance, or latent response's
#                 residual variance, is negative).
#   ordered     : TRUE for a factor of ordered items, whose omega is
#                 categorical omega.
#   omega       : TRUE when the factor's omega row is NA too (always for
#                 an ordered factor; for a continuous one, when the error
#                 variances sum below zero).
# Returns:
#   A single character string.
.kfit_improper_message <- function(factor, communality, ordered = TRUE,
                                   omega = TRUE) {
  items <- names(communality)
  vals <- formatC(communality, format = "f", digits = 3)
  one <- length(items) == 1L
  which <- if (ordered && one) {
    paste0("the latent response behind ", items, " has a negative ",
           "residual variance (its communality is ", vals, ")")
  } else if (ordered) {
    paste0("the latent responses behind ", paste(items, collapse = " and "),
           " have negative residual variances (their communalities are ",
           paste(vals, collapse = " and "), ")")
  } else if (one) {
    paste0(items, " has a negative error variance (its communality is ",
           vals, ")")
  } else {
    paste0(paste(items, collapse = " and "), " have negative error ",
           "variances (their communalities are ",
           paste(vals, collapse = " and "), ")")
  }
  rows <- paste0(c(if (omega) "omega_", "ave_", "H_"), factor)
  rows <- if (length(rows) == 3L) {
    paste0(rows[1L], ", ", rows[2L], ", and ", rows[3L])
  } else {
    paste(rows, collapse = " and ")
  }
  largest <- if (one) "that communality" else "the largest of them"
  ave_h <- paste0("the average variance extracted averages the ",
                  "communalities, and coefficient H, the largest ",
                  "reliability a weighted sum of the items attains, is at ",
                  "least ", largest, ", since the item weighted alone ",
                  "attains it")
  body <- if (ordered) {
    paste0("Categorical omega is not defined for such a solution, since ",
           "its assembly evaluates the bivariate normal distribution at ",
           "the communality; ", ave_h, ", so neither is a proportion of ",
           "variance either. Collapsing sparse categories, a larger ",
           "sample, or reconsidering the item are the usual remedies.")
  } else if (omega) {
    paste0("The error variances sum to a negative value, so the model ",
           "implies a negative error variance for the unit-weighted sum ",
           "and an omega above one; ", ave_h, ", so neither is a ",
           "proportion of variance either. A larger sample, or ",
           "reconsidering the model for these items, are the usual ",
           "remedies.")
  } else {
    paste0(toupper(substr(ave_h, 1L, 1L)), substring(ave_h, 2L),
           ", so neither is a proportion of variance here. omega_", factor,
           " is reported: the error variances still sum to a positive ",
           "value, so the reliability of the unit-weighted sum remains a ",
           "proportion. A larger sample, or reconsidering the model for ",
           "these items, are the usual remedies.")
  }
  paste0(rows, " are NA: the solution is improper (a Heywood case): ",
         which, ". ", body)
}


# .gauss_legendre(n)
#
# The n-point Gauss-Legendre rule on [-1, 1], by the construction of
# Golub and Welsch (1969): the nodes are the eigenvalues of the
# symmetric tridiagonal Jacobi matrix of the Legendre polynomials, whose
# off-diagonal elements are i / sqrt(4 i^2 - 1) for i = 1, ..., n - 1,
# and each weight is twice the squared first component of its node's
# unit eigenvector. The rule integrates a polynomial of degree 2n - 1
# exactly, which tests/testthat/test-pbinorm.R checks.
#
# Inputs:
#   n : the number of nodes, a whole number of at least two.
# Returns:
#   A list with x (the nodes, in decreasing order) and w (their weights,
#   which sum to two).
.gauss_legendre <- function(n) {
  i <- seq_len(n - 1L)
  off <- i / sqrt(4 * i^2 - 1)
  jacobi <- matrix(0, n, n)
  jacobi[cbind(i, i + 1L)] <- off
  jacobi[cbind(i + 1L, i)] <- off
  e <- eigen(jacobi, symmetric = TRUE)
  list(x = e$values, w = 2 * e$vectors[1L, ]^2)
}


# .GAUSS_LEGENDRE
#
# The 6-, 12-, and 20-point rules .pbinorm() uses, in that order,
# computed once when the package is installed rather than on every
# call.
.GAUSS_LEGENDRE <- lapply(c(6L, 12L, 20L), .gauss_legendre)


# .pbinorm(h, k, r)
#
# The bivariate standard normal distribution function,
# P(Z_1 <= h, Z_2 <= k) for standard normal Z_1 and Z_2 with correlation
# r, vectorized over its three arguments (recycled to a common length).
# Categorical omega evaluates it at every pair of thresholds of every
# pair of items (.gy_omega_from_components()); through version 1.0.1 it
# called mvtnorm::pmvnorm() once per such cell, so categorical omega
# needed a package beyond lavaan and spent most of a bootstrap
# replication in those calls.
#
# The method is that of Drezner and Wesolowsky (1990) with the
# refinements of Genz (2004), the algorithm mvtnorm itself applies to a
# two-dimensional normal probability (Genz's Fortran function MVBVU in
# the code mvtnorm compiles), ported step for step and vectorized. Genz writes the probability as
# the upper orthant P(X > a, Y > b), with a = -h and b = -k.
#   * For |r| < .925 it is Plackett's identity, the integral of the
#     bivariate normal density over the correlation from 0 to r, taken
#     over the arcsine of the correlation so the integrand stays smooth:
#       P = Phi(-a) Phi(-b) + asin(r) / (4 pi) *
#           int_{-1}^{1} exp((s a b - (a^2 + b^2) / 2) / (1 - s^2)) dx,
#     with s = sin(asin(r) (x + 1) / 2), by Gauss-Legendre quadrature
#     with 6, 12, or 20 nodes as |r| passes .3 and .75.
#   * For .925 <= |r| < 1 that integrand is steep, and the integral is
#     taken instead in sqrt(1 - r^2): a series approximation whose
#     integral is known in closed form, plus 20-node quadrature of the
#     difference between the exact integrand and the series. A negative
#     r is reflected to a positive one (b and a b change sign) and the
#     result reflected back.
#   * At |r| = 1 the limit is exact: Phi(min(h, k)) for r = 1, and
#     max(0, Phi(h) + Phi(k) - 1) for r = -1.
# An infinite limit gives the exact marginal or bound (0, 1, Phi(h), or
# Phi(k)). The normal distribution function is R's pnorm(), which is
# also what mvtnorm's Fortran calls.
#
# Accuracy, measured against mvtnorm 1.4.2 (tests/testthat/test-pbinorm.R
# pins its values; tools/oracle_snippets/mvtnorm.R makes the live
# comparison): on 2,873 grid points (h and k from -3 to 3 by .5, and 17
# correlations from -.999 to .999 on both sides of .3, .75, and .925)
# the largest difference is 4.0e-16, and on 20,000 random points with
# 1 - |r| above 1e-6 it is 2.2e-16. Nearer |r| = 1 mvtnorm loses
# accuracy, because it first reduces the nearly singular problem
# through a Cholesky factor (4e-14 at 1 - |r| = 1e-7, 2e-11 at 1e-10),
# and below about 1 - |r| = 5e-11 it returns the limit at |r| = 1, an
# error of 1.3e-6 at h = k = 0 and 1 - |r| = 3e-11. There this port
# agrees with the exact orthant probability at h = k = 0,
# 1/4 + asin(r) / (2 pi), to 1e-16, and elsewhere with numerical
# integration of the defining integral to 1e-15. Categorical omega
# meets such a correlation only for an item whose communality, or a
# pair whose polychoric correlation, is within 1e-7 of one, a solution
# at the edge of a Heywood case. A correlation within
# sqrt(.Machine$double.eps) of [-1, 1] is taken as its bound, as
# mvtnorm accepts it; for one between about 1e-11 and that far above
# 1 mvtnorm returns 0, which is wrong, and this function returns the
# limit.
#
# Inputs:
#   h, k : upper limits (numeric; recycled).
#   r    : correlations in [-1, 1] (recycled).
# Returns:
#   A numeric vector of the probabilities, of the common length. A
#   missing value in any argument, or a correlation outside [-1, 1],
#   is an error, as it was for mvtnorm::pmvnorm(): no caller should
#   produce one (an improper solution is refused before the assembly).
#
# References:
#   Drezner, Z., & Wesolowsky, G. O. (1990). On the computation of the
#     bivariate normal integral. Journal of Statistical Computation and
#     Simulation, 35(1--2), 101--107.
#   Genz, A. (2004). Numerical computation of rectangular bivariate and
#     trivariate normal and t probabilities. Statistics and Computing,
#     14(3), 251--260.
#   Golub, G. H., & Welsch, J. H. (1969). Calculation of Gauss
#     quadrature rules. Mathematics of Computation, 23(106), 221--230.
.pbinorm <- function(h, k, r) {
  # As for R's own distribution functions, an argument of length zero
  # gives a result of length zero.
  if (!length(h) || !length(k) || !length(r)) return(numeric(0))
  n <- max(length(h), length(k), length(r))
  h <- rep_len(as.numeric(h), n)
  k <- rep_len(as.numeric(k), n)
  r <- rep_len(as.numeric(r), n)
  if (anyNA(h) || anyNA(k) || anyNA(r)) {
    .dmar_refuse_argument(c(
      "The bivariate normal distribution function was given a ",
      "missing limit or correlation."), arg = NULL)
  }
  if (any(abs(r) > 1 + sqrt(.Machine$double.eps))) {
    .dmar_refuse_argument(c(
      "The bivariate normal distribution function needs ",
      "correlations in [-1, 1]; it was given ",
      format(r[abs(r) > 1 + sqrt(.Machine$double.eps)][1L],
             digits = 15), "."), arg = NULL)
  }
  r <- pmin(pmax(r, -1), 1)
  out <- numeric(n)
  inf <- is.infinite(h) | is.infinite(k)
  if (any(inf)) {
    hi <- h[inf]
    ki <- k[inf]
    out[inf] <- ifelse(hi == -Inf | ki == -Inf, 0,
                       ifelse(hi == Inf, stats::pnorm(ki),
                              stats::pnorm(hi)))
  }
  fin <- which(!inf)
  # The quadrature forms a matrix of one row per probability and one
  # column per node; blocks of 65,536 probabilities hold each such
  # matrix to about 10 MB however many cells the assembly sends.
  block <- 65536L
  n_blocks <- ceiling(length(fin) / block)
  for (s in seq.int(1L, by = block, length.out = n_blocks)) {
    i <- fin[s:min(length(fin), s + block - 1L)]
    out[i] <- .pbinorm_genz(-h[i], -k[i], r[i])
  }
  out
}


# .pbinorm_genz(a, b, r)
#
# The upper orthant P(X > a, Y > b) for finite a and b and r in
# [-1, 1], vectorized; the computation .pbinorm() describes. The
# variable names follow Genz's Fortran (hk = a b, hs = (a^2 + b^2) / 2,
# asr = asin(r), and as, bs, c, and d in the series for |r| >= .925,
# with c and d spelled cc and dd so as not to mask c()), so the port
# can be read beside it. In the steep branch Genz evaluates the two
# halves of the symmetric 20-node rule in two algebraically equal
# forms, the second with the common exponential factor taken out of
# the difference; the port does the same. The weighted sums over the
# nodes are rowSums() rather than a matrix product, so they do not
# depend on the BLAS R is linked to (rowSums() accumulates in extended
# precision where the platform has it); the two agree with mvtnorm
# equally well here.
.pbinorm_genz <- function(a, b, r) {
  out <- numeric(length(r))
  hk <- a * b
  ar <- abs(r)
  rule <- ifelse(ar < 0.3, 1L, ifelse(ar < 0.75, 2L, 3L))
  gentle <- ar < 0.925
  for (s in 1:3) {
    i <- which(gentle & rule == s)
    if (!length(i)) next
    gl <- .GAUSS_LEGENDRE[[s]]
    hs <- (a[i] * a[i] + b[i] * b[i]) / 2
    asr <- asin(r[i])
    sn <- sin(outer(asr, gl$x + 1) / 2)
    f <- exp((sn * hk[i] - hs) / (1 - sn * sn))
    out[i] <- rowSums(f * rep(gl$w, each = length(i))) * asr / (4 * pi) +
      stats::pnorm(-a[i]) * stats::pnorm(-b[i])
  }
  i <- which(!gentle)
  if (!length(i)) return(out)
  gl <- .GAUSS_LEGENDRE[[3L]]
  half <- gl$x < 0
  x <- gl$x[half]
  w <- gl$w[half]
  neg <- r[i] < 0
  ai <- a[i]
  bi <- ifelse(neg, -b[i], b[i])
  hki <- ifelse(neg, -hk[i], hk[i])
  bvn <- numeric(length(i))
  inner <- which(abs(r[i]) < 1)
  if (length(inner)) {
    ri <- r[i][inner]
    hkj <- hki[inner]
    as <- (1 - ri) * (1 + ri)
    A <- sqrt(as)
    bs <- (ai[inner] - bi[inner])^2
    cc <- (4 - hkj) / 8
    dd <- (12 - hkj) / 16
    v <- A * exp(-(bs / as + hkj) / 2) *
      (1 - cc * (bs - as) * (1 - dd * bs / 5) / 3 + cc * dd * as * as / 5)
    far <- hkj > -160
    if (any(far)) {
      B <- sqrt(bs[far])
      v[far] <- v[far] - exp(-hkj[far] / 2) * sqrt(2 * pi) *
        stats::pnorm(-B / A[far]) * B *
        (1 - cc[far] * bs[far] * (1 - dd[far] * bs[far] / 5) / 3)
    }
    A <- A / 2
    xs <- outer(A, x + 1)^2
    rs <- sqrt(1 - xs)
    g1 <- exp(-bs / (2 * xs) - hkj / (1 + rs)) / rs -
      exp(-(bs / xs + hkj) / 2) * (1 + cc * xs * (1 + dd * xs))
    xs <- outer(as, (-x + 1)^2) / 4
    rs <- sqrt(1 - xs)
    g2 <- exp(-(bs / xs + hkj) / 2) *
      (exp(-hkj * xs / (2 * (1 + rs)^2)) / rs - (1 + cc * xs * (1 + dd * xs)))
    v <- v + A * rowSums((g1 + g2) * rep(w, each = length(inner)))
    bvn[inner] <- -v / (2 * pi)
  }
  pos <- !neg
  bvn[pos] <- bvn[pos] + stats::pnorm(-pmax(ai[pos], bi[pos]))
  bvn[neg] <- -bvn[neg]
  low <- neg & bi > ai & ai < 0
  bvn[low] <- bvn[low] + stats::pnorm(bi[low]) - stats::pnorm(ai[low])
  high <- neg & bi > ai & ai >= 0
  bvn[high] <- bvn[high] + stats::pnorm(-ai[high]) - stats::pnorm(-bi[high])
  out[i] <- bvn
  out
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


# .polychoric_from_fit(fit, var_names)
#
# The unrestricted polychoric correlation matrix of the items, taken
# from the sample statistics of the single-factor fit (the pairwise
# maximum likelihood polychorics lavaan estimates before fitting any
# structure, which lavaan::lavCor() also returns). Before DMAR 1.1.0
# (through the 1.0.0 release and the 1.0.1 code) this matrix came from
# a second fit, a saturated model in which every pair of items was
# correlated (the legacy MBESS approach). A converged
# saturated fit reproduces the sample polychorics exactly, so the two
# agree whenever every polychoric is determined by the data. They part
# ways when one is not: the diagonally weighted fit gives that
# correlation zero weight and returns its starting value, 0, so the
# refit silently removed the pair from the total variance. Reading the
# sample statistics directly removes that failure and a second fit per
# bootstrap replication; .undetermined_polychorics() identifies the
# pairs the data do not determine.
#
# Inputs:
#   fit       : a converged lavaan fit with the items as ordered
#               indicators.
#   var_names : character vector of item names as lavaan knows them.
# Returns:
#   A (J x J) numeric matrix of polychoric correlations, rows and
#   columns in the order of var_names.
.polychoric_from_fit <- function(fit, var_names) {
  poly <- unclass(lavaan::lavInspect(fit, "sampstat")$cov)
  poly[var_names, var_names, drop = FALSE]
}


# .undetermined_polychorics(fit, var_names, labels = var_names)
#
# Item pairs whose polychoric correlation the data do not determine:
# the asymptotic standard error lavaan reports for the sample
# polychoric (from the Gamma matrix of the fit) is not finite or
# exceeds 1, the width of half the parameter space, so the data place
# no useful bound on the correlation. This happens when the bivariate
# table of two items has empty cells in the pattern that two items
# skewed in opposite directions produce (for example, every respondent
# below the top category of one item in the bottom category of the
# other). The single-factor fit gives such a correlation zero weight,
# so the loadings and the model implied correlations rest on the
# correlations that are determined, and remain estimable only when those
# identify them (.omega_c_identification()); the unrestricted
# polychoric, and with it the observed denominator, does not.
#
# Inputs:
#   fit       : a converged lavaan fit with the items as ordered
#               indicators.
#   var_names : character vector of item names as lavaan knows them.
#   labels    : names to report for the items, parallel to var_names.
# Returns:
#   A two-column character matrix with one row per undetermined pair
#   (zero rows when every polychoric is determined or when the Gamma
#   matrix is unavailable).
.undetermined_polychorics <- function(fit, var_names, labels = var_names) {
  none <- matrix(character(0), ncol = 2L)
  gam <- try(lavaan::lavInspect(fit, "gamma"), silent = TRUE)
  if (inherits(gam, "try-error") || is.null(rownames(gam))) return(none)
  n_obs <- lavaan::lavInspect(fit, "nobs")
  se <- sqrt(pmax(diag(gam), 0) / n_obs)
  pairs <- which(upper.tri(diag(length(var_names))), arr.ind = TRUE)
  out <- none
  for (k in seq_len(nrow(pairs))) {
    a <- var_names[pairs[k, 1]]; b <- var_names[pairs[k, 2]]
    nm <- intersect(c(paste0(a, "~~", b), paste0(b, "~~", a)), names(se))
    if (length(nm) && (!is.finite(se[nm[1]]) || se[nm[1]] > 1)) {
      out <- rbind(out, labels[pairs[k, ]])
    }
  }
  out
}


# .undetermined_message(pairs)
#
# The phrase naming the item pairs .undetermined_polychorics() returned,
# shared by the error, warning, and message texts that report them.
#
# Inputs:
#   pairs : two-column character matrix with at least one row.
# Returns:
#   A single character string, e.g. "the polychoric correlation of y1
#   and y3" or "the polychoric correlations of y1 and y3 and of y2 and
#   y4".
.undetermined_message <- function(pairs) {
  each <- paste(pairs[, 1], "and", pairs[, 2])
  if (length(each) == 1L) {
    paste("the polychoric correlation of", each)
  } else {
    paste("the polychoric correlations of",
          paste(each, collapse = " and of "))
  }
}


# .omega_c_identification(items, undetermined)
#
# Whether the polychoric correlations the data determine identify the
# loadings of the single-factor model. The thresholds are saturated, and
# a correlation the data do not determine gets a weight of about 1e-13 in
# the diagonally weighted fit (the inverse of its unbounded asymptotic
# variance), so the loadings answer only to the determined correlations,
# rho_jk = lambda_j lambda_k, one for each edge of the graph whose
# vertices are the items and whose edges are the determined pairs. Taking
# logarithms of the absolute values gives a linear system in
# log |lambda_j| whose matrix is the graph's unsigned incidence matrix,
# of rank J minus the number of bipartite components; the signs solve
# s_j s_k = sign(rho_jk), of rank J minus the number of components. So
# the loadings are identified, up to the reflection of all of them that
# omega does not see, exactly when the graph is connected and not
# bipartite, that is, when it is connected and has a cycle of an odd
# number of items (three items whose three correlations are determined
# is the shortest). A disconnected graph leaves the relative sign of its
# blocks free, and a bipartite block the scale of its two sides (one
# side's loadings multiplied by c and the other's divided by c fit as
# well). In the usual sample every correlation is determined, and the
# complete graph on three or more items passes. The check is of the
# graph only: a determined correlation near zero leaves the loadings
# weakly identified however the graph looks, which it does not test.
#
# On Green and Yang's (2009) mixed-skew design (two items piled into
# their top category and two into their bottom; N = 300), 13 of 60
# samples failed the check. In the 2 with every cross pair undetermined
# (two blocks), lavaan's default start values gave loadings of opposite
# signs and a categorical omega of 0.05 and 0.09, against a population
# value of 0.716. The other 11 were bipartite (one block): in 8 the
# value moved with the start values by up to 0.06, and 3 gave an
# improper solution at the default start values. In the samples that
# passed, every start value gave the same omega to four decimals.
#
# Inputs:
#   items        : the item names.
#   undetermined : the two-column matrix .undetermined_polychorics()
#                  returns.
# Returns:
#   A list with `identified` (TRUE or FALSE), `blocks` (the connected
#   components of the graph, a list of item name vectors in item order),
#   and `sides` (for each block, its two sides when it is bipartite, a
#   single item counting as bipartite, and NULL otherwise).
.omega_c_identification <- function(items, undetermined) {
  J <- length(items)
  adj <- matrix(TRUE, J, J, dimnames = list(items, items))
  diag(adj) <- FALSE
  for (k in seq_len(nrow(undetermined))) {
    adj[undetermined[k, 1L], undetermined[k, 2L]] <- FALSE
    adj[undetermined[k, 2L], undetermined[k, 1L]] <- FALSE
  }
  # Breadth-first search from each item not yet reached, coloring the
  # items alternately; an edge between two items of one color closes a
  # cycle of an odd number of items.
  block <- integer(J)
  color <- integer(J)
  odd <- logical()
  for (s in seq_len(J)) {
    if (block[s] > 0L) next
    b <- length(odd) + 1L
    odd[b] <- FALSE
    block[s] <- b
    queue <- s
    while (length(queue)) {
      v <- queue[1L]
      queue <- queue[-1L]
      for (w in which(adj[v, ])) {
        if (block[w] == 0L) {
          block[w] <- b
          color[w] <- 1L - color[v]
          queue <- c(queue, w)
        } else if (color[w] == color[v]) {
          odd[b] <- TRUE
        }
      }
    }
  }
  blocks <- lapply(seq_along(odd), function(b) items[block == b])
  sides <- lapply(seq_along(odd), function(b) {
    if (odd[b]) return(NULL)
    in_b <- block == b
    list(items[in_b & color == 0L], items[in_b & color == 1L])
  })
  list(identified = length(odd) == 1L && odd[1L], blocks = blocks,
       sides = sides)
}


# .omega_c_unidentified_message(identification, undetermined)
#
# Why categorical omega cannot be computed when the determined
# polychoric correlations do not identify the loadings
# (.omega_c_identification()): which correlations the data do not
# determine, how the determined ones leave the loadings free, and that
# this holds under either denominator, since the numerator takes the
# loadings.
#
# Inputs:
#   identification : the list .omega_c_identification() returned, with
#                    `identified` FALSE.
#   undetermined   : the undetermined pairs (at least one row).
# Returns:
#   A single character string.
.omega_c_unidentified_message <- function(identification, undetermined) {
  set <- function(v) {
    if (length(v) == 1L) v else paste0("{", paste(v, collapse = ", "), "}")
  }
  blocks <- identification$blocks
  sides <- identification$sides
  how <- if (length(blocks) > 1L) {
    alone <- lengths(blocks) == 1L
    free_scale <- !vapply(sides, is.null, logical(1L)) & !alone
    paste0(
      "no correlation the data determine links ",
      paste(vapply(blocks, set, character(1L)), collapse = " to "),
      ", so the relative sign of the loadings of these blocks is ",
      "arbitrary",
      if (any(alone)) {
        paste0(", and the loading of ",
               paste(unlist(blocks[alone]), collapse = " and of "),
               ", which shares a determined correlation with no item, ",
               "is arbitrary altogether")
      },
      if (any(free_scale)) {
        paste0("; within ",
               paste(vapply(blocks[free_scale], set, character(1L)),
                     collapse = " and within "),
               ", whose determined correlations form no cycle of an odd ",
               "number of items, the scale of the loadings is arbitrary ",
               "as well")
      })
  } else {
    s <- sides[[1L]]
    paste0(
      "every correlation the data determine joins an item of ",
      set(s[[1L]]), " to an item of ", set(s[[2L]]), ", so they form no ",
      "cycle of an odd number of items (three items whose three ",
      "correlations are all determined, for example) and fix only ",
      "products of a loading of one set and a loading of the other: the ",
      "loadings of one set multiplied by any number and those of the ",
      "other divided by it fit as well")
  }
  paste0(
    "The single-factor loadings are not identified: the data do not ",
    "determine ", .undetermined_message(undetermined), " (the bivariate ",
    "tables have empty cells, which two items skewed in opposite ",
    "directions produce), a correlation the data do not determine carries ",
    "no weight in the fit, and ", how, ". The loadings, and categorical ",
    "omega under either denominator, would be whatever the start values ",
    "of the fit led to. Respondents in the empty cells, a larger sample ",
    "in practice, are what determine those correlations.")
}


# .omega_c_undetermined_note(denominator)
#
# The sentence a bootstrap warning gives about the replications set
# aside because a resample left a polychoric correlation undetermined,
# for the bootstrap of the point estimate (.bootstrap_ci()) and for
# bootstrap then impute (.boot_mi_notes()). Under the observed
# denominator a replication is set aside for any such correlation; under
# the model implied one only when the correlations left determined do
# not identify the loadings (.omega_c_identification()). Through the
# first versions of 1.1.0 the sentence said that the model implied
# denominator keeps every such replication, which kept replications
# whose loadings were arbitrary.
.omega_c_undetermined_note <- function(denominator) {
  start <- paste0(
    "Resampling thins the rare cells of a bivariate table, so a ",
    "resample can leave a polychoric correlation undetermined that the ",
    "sample determines")
  if (identical(denominator, "model_implied")) {
    paste0(start, ". The model implied denominator needs no unrestricted ",
           "polychoric correlation, but the loadings rest on the ",
           "correlations left determined, which identify them only when ",
           "they connect every item to the others, directly or through ",
           "other items, and include at least one cycle of an odd number ",
           "of items (three items whose three correlations are all ",
           "determined, for example); a replication in which they do not ",
           "has no categorical omega, and the interval conditions on the ",
           "replications in which they do")
  } else {
    paste0(start, ". The interval conditions on the replications ",
           "that redrew the respondents who determine it, not on all of ",
           "them; denominator = \"model_implied\" needs no unrestricted ",
           "polychoric correlation and sets aside only the replications ",
           "whose determined correlations leave the loadings unidentified")
  }
}


# .residual_cor_summary(R, df = NULL, tol = 1e-6)
#
# Local dependence summary for a single-factor reliability fit: the
# residual correlation matrix (observed or unrestricted polychoric
# correlations minus those the fitted single-factor model reproduces)
# and its largest off-diagonal element in absolute value, with the item
# pair it belongs to. A denominator estimated from the data protects
# the denominator of omega against structure the single-factor model
# omits, but not the numerator: a doublet or a method factor among a
# few items raises their loadings, and omega with them, under either
# denominator. The residual correlations are where that structure
# usually shows, so the reliability functions report them beside the
# coefficient. When the structure sits on the items that dominate the
# factor, the single factor absorbs it and the residuals appear among
# the other items instead (a population fit with a doublet on the two
# strongest of six items puts the largest residual on another pair), so
# the largest residual signals local dependence more reliably than it
# locates it. No cutoff is applied; the number is reported, not
# labeled.
#
# A saturated fit is recognized rather than summarized. With three items
# the single-factor model has as many free parameters as correlations
# (df = 0), and a proper solution reproduces every correlation, so its
# residuals are zero up to the optimizer's tolerance and the "largest"
# of them, with the pair it falls on, is noise. When df is zero and no
# residual exceeds `tol` the summary says so (`saturated = TRUE`) and
# names no pair. When df is zero but a residual does exceed `tol` the
# model cannot reproduce the correlations (with three items that
# happens when the product of the three correlations is negative, so no
# real loadings fit them), and the residual is reported as usual. The
# default tolerance sits well above the optimizer noise of a saturated
# fit (about 1e-8 for the continuous fit, 1e-16 for the categorical
# one) and well below any residual worth printing to three decimals.
#
# Inputs:
#   R   : a (J x J) residual correlation matrix with item names as
#         dimnames, or NULL.
#   df  : the degrees of freedom of the fitted correlation structure, or
#         NULL when unknown (no saturation check).
#   tol : the largest absolute residual still read as zero in a
#         saturated fit.
#
# Returns:
#   NULL when R is NULL, has fewer than two items, or has no finite
#   residual; otherwise a list with `matrix` (R, symmetrized),
#   `largest` (the signed residual correlation largest in absolute
#   value; NA for a saturated fit), `items` (the two item names it
#   belongs to; NULL for a saturated fit), and `saturated` (TRUE when
#   df is zero and every residual is zero within `tol`).
.residual_cor_summary <- function(R, df = NULL, tol = 1e-6) {
  if (is.null(R)) return(NULL)
  R <- as.matrix(R)
  if (nrow(R) < 2L || nrow(R) != ncol(R)) return(NULL)
  nm <- rownames(R)
  if (is.null(nm)) nm <- paste0("y", seq_len(nrow(R)))
  R <- (R + t(R)) / 2
  dimnames(R) <- list(nm, nm)
  idx <- which(upper.tri(R), arr.ind = TRUE)
  vals <- R[idx]
  if (!any(is.finite(vals))) return(NULL)
  saturated <- isTRUE(all.equal(as.numeric(df), 0)) &&
    max(abs(vals[is.finite(vals)])) <= tol
  if (saturated) {
    return(list(matrix = R, largest = NA_real_, items = NULL,
                saturated = TRUE))
  }
  k <- which.max(abs(replace(vals, !is.finite(vals), 0)))
  list(matrix = R, largest = vals[k], items = nm[idx[k, ]],
       saturated = FALSE)
}


# .omega_residual_cor(fit, items = NULL, lavaan_items = items)
#
# Residual correlations of a fitted single-factor lavaan model (the
# observed minus the model implied correlations, lavaan's "cor.bollen"
# residuals), restricted to the items so that auxiliary variables
# entered as saturated correlates drop out, named by the user's item
# names, and summarized by .residual_cor_summary() with the fit's
# degrees of freedom. The auxiliary variable branch of .omega_fit_cfa()
# names the items y1, ..., yJ and the auxiliaries a1, ..., aK inside the
# fit (.omega_fit_cfa_aux()), so the items are found by the names the fit
# gives them (`lavaan_items`) and then given the user's names. A lookup by
# the user's names alone matched nothing there, and the summary kept the
# auxiliary rows and named its pair by position ("y2 with y3"). Returns
# NULL rather than an error when the residuals cannot be computed or the
# fit does not hold every item: the summary is a diagnostic beside the
# coefficient and never a reason to withhold it, and a summary over other
# variables than the items would mislead.
#
# Inputs:
#   fit          : a converged lavaan fit (the `fit` element returned by
#                  .omega_fit_cfa()), or NULL.
#   items        : character vector of item names as the user gave them,
#                  or NULL for all observed variables of the fit.
#   lavaan_items : the same items as the fit names them, in the order of
#                  `items` (the `lavaan_items` element of
#                  .omega_fit_cfa()); by default the user's names.
#
# Returns:
#   The .residual_cor_summary() list, or NULL.
.omega_residual_cor <- function(fit, items = NULL, lavaan_items = items) {
  if (is.null(fit)) return(NULL)
  res <- try(suppressWarnings(lavaan::residuals(fit, type = "cor")),
             silent = TRUE)
  if (inherits(res, "try-error") || is.null(res$cov)) return(NULL)
  R <- unclass(res$cov)
  R <- matrix(as.numeric(R), nrow(R), ncol(R),
              dimnames = list(rownames(res$cov), colnames(res$cov)))
  if (!is.null(lavaan_items)) {
    if (!all(lavaan_items %in% rownames(R))) return(NULL)
    R <- R[lavaan_items, lavaan_items, drop = FALSE]
    if (!is.null(items)) dimnames(R) <- list(items, items)
  }
  # Saturated auxiliary correlates add no degrees of freedom, so the
  # fit's df is that of the single-factor structure among the items.
  df <- tryCatch(unname(lavaan::fitMeasures(fit, "df")),
                 error = function(e) NULL)
  .residual_cor_summary(R, df = df)
}


# .omega_ml_residual_cor(inp, missing, fit = NULL, estimator = NULL)
#
# The local dependence summary reliability_omega() attaches as its
# residual_cor attribute, always from a maximum likelihood fit of the
# single-factor model to the cases the coefficient used, so the summary
# does not depend on the interval method. ML and MLR share their point
# estimates, so a main fit by either is reused as it is; any other fit
# (the WLS fit behind ci_method = "adf" or "adf_logistic", whose
# residuals differ from the ML ones and can fall on a different pair),
# or no fit at all (the observed denominator's point estimate does not
# leave one behind), is replaced by one ML fit made here. It runs once,
# never inside the bootstrap. Under every missing data option the summary
# covers the items only, under the user's names: with auxiliary variables
# under FIML the fit renames the items and carries the auxiliaries beside
# them (.omega_fit_cfa_aux()), so the items are found through the fit's
# `lavaan_items`.
#
# Inputs:
#   inp       : the list .relia_resolve_inputs() returned.
#   missing   : "listwise", "fiml", or "two_stage".
#   fit       : the main fit, the list .omega_fit_cfa() returned (the
#               lavaan object in `fit`, the item names in `items` and
#               `lavaan_items`), or NULL.
#   estimator : the estimator of `fit`, or NULL.
#
# Returns:
#   The .omega_residual_cor() list, or NULL when the ML fit did not
#   converge or the residuals cannot be computed.
.omega_ml_residual_cor <- function(inp, missing, fit = NULL,
                                   estimator = NULL) {
  if (!is.null(inp$h1_diagnostic)) return(NULL)
  if (!is.null(fit$fit) && isTRUE(estimator %in% c("ML", "MLR"))) {
    return(.omega_residual_cor(fit$fit, fit$items %||% inp$items,
                               fit$lavaan_items %||% inp$items))
  }
  ml_fit <- .omega_fit_cfa(data = inp$data,
                           S = if (is.null(inp$data)) inp$S else NULL,
                           N = if (is.null(inp$data)) inp$N else NULL,
                           equal_loadings = FALSE,
                           estimator = "ML", se = "none",
                           missing = .relia_lavaan_missing(missing),
                           aux = inp$aux)
  if (!isTRUE(ml_fit$converged)) return(NULL)
  .omega_residual_cor(ml_fit$fit, ml_fit$items %||% inp$items,
                      ml_fit$lavaan_items %||% inp$items)
}


# .denominator_phrase(denominator, coefficient = NULL, metric = NULL)
#
# The words the printed table uses to name how the total variance in the
# denominator of an omega coefficient was estimated, so a reader of the
# output never has to know the default to know which coefficient is on
# the page. Used by .dmar_tbl_footer() for any table that carries a
# `denominator` attribute (reliability_omega(),
# reliability_omega_categorical(), bifactor_indices()). A bifactor fit
# to ordered items carries `metric = "latent_response"`: its loadings
# and both of its denominators are on the metric of the continuous
# latent responses behind the categories, not of the summed categories,
# and the phrase says so.
#
# Inputs:
#   denominator : "observed" or "model_implied".
#   coefficient : the table's `coefficient` attribute, or NULL.
#   metric      : the table's `metric` attribute ("observed" or
#                 "latent_response"), or NULL.
#
# Returns:
#   A single character string.
.denominator_phrase <- function(denominator, coefficient = NULL,
                                metric = NULL) {
  if (identical(coefficient, "omega_categorical")) {
    return(switch(
      denominator,
      observed = paste("variance of the summed categories from the",
                       "unrestricted polychoric correlations",
                       "(Kelley & Pornprasertmanit, 2016)"),
      model_implied = paste("variance of the summed categories from the",
                            "model implied polychoric correlations",
                            "(Green & Yang, 2009)"),
      denominator))
  }
  if (identical(metric, "latent_response")) {
    return(switch(
      denominator,
      observed = paste("variance of the sum of the latent responses,",
                       "from the polychoric correlations"),
      model_implied = paste("variance of the sum of the latent responses",
                            "implied by the fitted model"),
      denominator))
  }
  switch(
    denominator,
    observed = if (identical(coefficient, "omega")) {
      "observed composite variance (robust omega)"
    } else {
      "observed composite variance"
    },
    model_implied = "composite variance implied by the fitted model",
    denominator)
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

# .relia_two_rows_alternative(kr20 = FALSE)
#
# The interval that remains when an interval refuses two rows (N = 2, the
# complete rows under listwise deletion): Feldt's (1965) F interval, whose
# reference, F on N - 1 and (N - 1)(J - 1) degrees of freedom, is exact
# for any N of at least 2 under multivariate normality and parallel items
# (the normal persons by items model with equal error variances), and for
# 0/1 items, which cannot be normal, an approximation to the sampling
# distribution of KR-20 (Feldt, 1965, whose title says so). Two rows give
# a covariance matrix of rank one, every pair of items correlating +1 or
# -1 whatever the items, and alpha (and KR-20) is then a function of the
# difference between the two rows alone; the intervals that read that
# matrix as the population's ("ml", "adf", and their logistic variants),
# and the bootstrap, whose every resample is one row twice or both rows,
# have no width, and Bonett's has an infinite variance. The cube root
# interval of Hakstian and Whalen (1976) is defined at N = 2, but its
# upper limit then lies above one for every data set with alpha below one
# at any confidence level above about .90 (.ci_hakstian_whalen()), so it
# is not named here.
# Shared by .ci_bonett(), .ci_alpha_ml(), and .ci_alpha_adf(), which take
# the sentence as an argument, and by the callers of .bootstrap_ci() that
# have Feldt's interval to offer (reliability_alpha()'s analytic
# estimator under listwise deletion, and reliability_kr20()); any other
# caller's bootstrap names none.
#
# Inputs:
#   kr20 : TRUE when the caller is reliability_kr20().
# Returns:
#   A sentence.
.relia_two_rows_alternative <- function(kr20 = FALSE) {
  if (kr20) {
    paste0("With two rows use ci_method = \"feldt\", whose F reference ",
           "approximates the sampling distribution of KR-20 (Feldt, 1965) ",
           "for any 'N' of at least 2.")
  } else {
    paste0("With two rows use ci_method = \"feldt\", whose F reference is ",
           "exact for any 'N' of at least 2 under multivariate normality ",
           "and parallel items.")
  }
}


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
    .dmar_refuse_argument(c(
      "The Fisher confidence interval requires 'N' of at least 4: the ",
      "variance of the transformed coefficient is 1/(N - 3), which is ",
      "infinite at N = 3, where the interval would be vacuous, and ",
      "undefined for smaller 'N'."), arg = "N")
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
# |dr/dz| = exp(z) = 1 - r at the estimate). Requires N >= 3: the
# variance is infinite at N = 2, where the interval would be vacuous (the
# limits clamped to [0, 1] and an infinite standard error through 1.0.1),
# as Fisher's is at N = 3 (.ci_fisher()).
.ci_bonett <- function(estimate, N, J, conf_level,
                       two_rows = .relia_two_rows_alternative()) {
  if (N < 3) {
    .dmar_refuse_argument(c(
      "The Bonett confidence interval requires 'N' of at least 3: the ",
      "variance of the transformed coefficient is 2J/((J - 1)(N - 2)), ",
      "which is infinite at N = 2, where the interval would be ",
      "vacuous. ", two_rows), arg = "N")
  }
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
# estimate). At N = 2 (df_n = 1) the upper limit on the cube root scale
# is (1 - r)^(1/3) (1 - crit sqrt(2J/(9(J - 1))) / (7/9)), whose second
# factor does not depend on the data and is negative whenever crit
# exceeds (7/9) / sqrt(2J/(9(J - 1))), which is below 1.650 for every J
# (1.167 at J = 2, 1.650 as J grows): at any confidence level above
# about .90 the upper limit lies above one for every data set with alpha
# below one and is clamped to one (at alpha = 1, which two rows can give,
# the standard error is zero and the interval is [1, 1]). The interval is
# defined there, so it is not refused,
# but the refusals of two rows name Feldt's interval alone.
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
    .dmar_refuse(
      "The likelihood interval requires the 'lavaan' package.",
      "package_required", package = "lavaan")
  }

  if (!is.null(data)) {
    data <- as.data.frame(data)
    data <- data[stats::complete.cases(data), , drop = FALSE]
    J <- ncol(data)
    item_names <- colnames(data)
    colnames(data) <- paste0("y", seq_len(J))
  } else {
    J <- ncol(S)
    item_names <- colnames(S)
    dimnames(S) <- list(paste0("y", seq_len(J)), paste0("y", seq_len(J)))
  }
  vars <- paste0("y", seq_len(J))
  if (is.null(item_names)) item_names <- vars

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

  # Only the log likelihood of each fit is used. The test is named so that
  # lavaan (0.7-3 and later) does not also attempt Browne's residual based
  # statistic, which it cannot compute under the nonlinear equality
  # constraint of the profile and would warn about at every step.
  # lavaan's warnings are raised again with the caller's names for the
  # items (.relia_with_user_names()).
  fit_lav <- function(syntax, start = NULL) {
    .relia_with_user_names(if (!is.null(data)) {
      lavaan::cfa(syntax, data = data, estimator = "ML", se = "none",
                  test = "standard", start = start)
    } else {
      lavaan::cfa(syntax, sample.cov = S, sample.nobs = N,
                  estimator = "ML", se = "none", test = "standard",
                  start = start)
    }, vars, item_names)
  }

  fit0 <- fit_lav(base_syntax)
  if (!lavaan::lavInspect(fit0, "converged")) {
    .dmar_warn(c(
      "The unconstrained single-factor fit did not converge; ",
      "the likelihood interval is unavailable."), "result_qualified")
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
    .dmar_warn(c(
      "A constrained fit failed to converge during the profile ",
      "likelihood search; the affected bound is NA."), "result_qualified")
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
# back-transformed. Each takes wald_df, the degrees of freedom of the t
# reference distribution (Inf, the normal reference, by default).
# -----------------------------------------------------------------------------

# The interval methods of reliability_omega() and reliability_alpha() that
# are Wald intervals, the ones whose reference distribution wald_df sets
# (for omega, under the model implied denominator only).
.RELIA_WALD_METHODS <- c("ml", "ml_logistic", "mlr", "mlr_logistic",
                         "adf", "adf_logistic")

# .relia_wald_df(wald_df, fn, fit = NULL, N, J)
#
# The degrees of freedom of a reliability coefficient's Wald interval:
# wald_df itself when it is a number, and under the keyword the rule
# .wald_df_rule() holds for fn (reliability_omega() only: von Hippel's
# max(3, N(1 - FMI) - 1) from the coefficient's own fraction of missing
# information, the variant the calibration study chose; see section (F)
# of R/missing_information_internals.R). `fit` is the list
# .omega_fit_cfa() or .fiml_cov() returned, whose lavaan fit and defined
# parameter name a rule reads, or NULL for a closed form computed from
# complete data. Returns list(df, table): df unnamed, so the interval
# limits stay unnamed, and table the df_missing() table under the keyword
# (NULL otherwise), which the result carries for the line under it.
#
# The fit names the items y1, ..., yJ, the auxiliary variables a1, ...,
# aK, and the factor f1, so the warnings and errors of the computation
# behind the keyword (the fraction of missing information of the fit,
# .fmi_core()) are given again with the caller's names `items` and `aux`
# (.relia_lavaan_names_back(); an error keeps its class), as the family's
# lavaan warnings are (.relia_with_user_names()). The warning that the
# fit is not at the maximum of its likelihood
# (dmar_fmi_not_at_maximum, from .fmi_check_maximum()) is written anew:
# its parameter is named in words ("the loading of q4"), and its remedy,
# to refit with rescale_data = TRUE, is one the caller cannot follow,
# since the function builds the fit; the remedy given is the one the data
# allow, a common rescaling of the items, which leaves omega as it is
# (rescaling items one at a time changes the composite, and omega with
# it), and a rescaling of each auxiliary variable. Through the first
# versions of 1.1.0 the warning reached the user as "the estimate of y4 ~~
# a1" with the refit advice.
.relia_wald_df <- function(wald_df, fn, fit = NULL, N, J, items = NULL,
                           aux = NULL) {
  lavaan_fit <- if (is.null(fit)) NULL else fit$fit
  name <- if (is.null(fit)) "reliability" else fit$coef_name
  rows <- data.frame(lhs = name, op = ":=", rhs = "",
                     stringsAsFactors = FALSE)
  coefficient <- if (identical(fn, "reliability_alpha")) "alpha" else "omega"
  resolve <- function() {
    .wald_df_resolve(wald_df, fn, fit = lavaan_fit, rows = rows,
                     context = list(N = N, n_items = J, omega = name,
                                    terms = coefficient))
  }
  if (is.null(lavaan_fit) || is.null(items)) {
    res <- resolve()
    return(list(df = unname(res$df), table = res$table))
  }
  internal <- c(paste0("y", seq_along(items)),
                if (length(aux)) paste0("a", seq_along(aux)))
  display <- c(items, aux)
  rename <- function(text) .relia_lavaan_names_back(text, internal, display)
  word <- function(v) {
    k <- match(v, internal)
    if (is.na(k)) v else display[k]
  }
  res <- withCallingHandlers(
    tryCatch(resolve(), error = function(e) {
      e$message <- rename(conditionMessage(e))
      e$call <- NULL
      stop(e)
    }),
    dmar_fmi_not_at_maximum = function(w) {
      fmt <- function(v) format(v, digits = 4L)
      what <- if (identical(w$op, "=~")) {
        paste("the loading of", word(w$rhs))
      } else if (identical(w$op, "~~") && identical(w$lhs, w$rhs)) {
        paste(if (w$lhs %in% paste0("y", seq_along(items))) {
          "the error variance of"
        } else {
          "the variance of"
        }, word(w$lhs))
      } else if (identical(w$op, "~~")) {
        paste("the covariance of", word(w$lhs), "and", word(w$rhs))
      } else if (identical(w$op, "~1")) {
        paste("the mean of", word(w$lhs))
      } else {
        rename(trimws(paste(w$lhs, w$op, w$rhs)))
      }
      .dmar_warn(
        paste0(
  "The single-factor fit behind the degrees of freedom of the ",
  coefficient, " interval (wald_df = \"missing_information\") is ",
  "not at the maximum of its likelihood: one Newton step from its ",
  "estimates moves ", what, " from ", fmt(w$from), " to ",
  fmt(w$to), " and raises the log likelihood by ", fmt(w$gain),
  if (!is.null(w$change)) {
    paste0(", changing the fraction of missing information of ",
           rename(w$change$term), " from ", fmt(w$change$from),
           " to ", fmt(w$change$to))
  },
  ", so the ", coefficient, " estimate, its interval, and the ",
  "fraction of missing information behind the degrees of freedom ",
  "are not quite the maximum likelihood values.",
  if (isFALSE(w$units)) {
    paste0(" lavaan's optimizer can stop within its convergence ",
           "tolerance on a likelihood that is flat in some ",
           "direction, as when few values of an item are observed. ",
           "The change above measures the resulting error in the ",
           "fractions of missing information.")
  } else {
    paste0(
      " lavaan can stop short of the maximum for variables ",
      "measured in units far from one another that it does not ",
      "rescale. Divide the items by one common constant, which ",
      "leaves ", coefficient, " as it is",
      if (length(aux)) {
        ", and each auxiliary variable by a constant of its own"
      }, ", so that the standard deviations are near 1, and call ",
      fn, "() again (rescaling the items one at a time would ",
      "change the composite and its ", coefficient, ").")
  }), "result_qualified")
      invokeRestart("muffleWarning")
    },
    warning = function(w) {
      named <- rename(conditionMessage(w))
      if (!identical(named, conditionMessage(w))) {
        .dmar_warn(named, "result_qualified")
        invokeRestart("muffleWarning")
      }
    })
  list(df = unname(res$df), table = res$table)
}

# .alpha_vanzyl_avar(S)
#
# The asymptotic variance of sample coefficient alpha, times N - 1, under
# multivariate normality of the items (van Zyl, Neudecker, & Nel, 2000,
# Section 3; Terry & Kelley, 2012, Eqs. 7 and 8; Kelley &
# Pornprasertmanit, 2016, Eq. 20):
#
#     (J / (J - 1))^2 * ( 2 / (1' S 1)^3 )
#       * [ (1' S 1) * (tr(S^2) + (tr S)^2) - 2 * (tr S) * (1' S^2 1) ]
#
# where S is the item covariance matrix (the population one when planning,
# the sample one when estimating), 1 a J-vector of ones, and tr() the
# trace. The expression is unchanged when S is multiplied by a constant,
# but not when the items are rescaled one by one, so it is evaluated on the
# covariance matrix itself, never on cov2cor(S): on the correlation matrix
# it is the variance of the alpha of standardized items, which differs from
# the variance of alpha of the raw items whenever the item variances
# differ. Before DMAR 1.1.0, .ci_alpha_ml() evaluated it on cov2cor(S), as
# MBESS's ci.reliability() (its internal .interval.type31_32_type1()) and
# ss.aipe.reliability() still do; the standard error was then too small or
# too large for raw alpha, depending on how the item variances lined up
# with the loadings. Shared by .ci_alpha_ml() and the normal theory closed
# form of ss_aipe_reliability() for the tau equivalent model (its factor
# analytic closed form is the delta method variance of the model implied
# coefficient, .ss_aipe_omega_avar(equal_loadings = TRUE)).
.alpha_vanzyl_avar <- function(S) {
  S <- as.matrix(S)
  J <- ncol(S)
  one      <- matrix(1, nrow = J, ncol = 1)
  S2       <- S %*% S
  oneSone  <- as.numeric(t(one) %*% S %*% one)
  trS      <- sum(diag(S))
  trS2     <- sum(diag(S2))
  oneS2one <- as.numeric(t(one) %*% S2 %*% one)
  (J^2 / (J - 1)^2) * (2 / oneSone^3) *
    (oneSone * (trS2 + trS^2) - 2 * trS * oneS2one)
}

# .ci_alpha_ml(estimate, S, N, J, conf_level, logistic = FALSE,
#              wald_df = Inf)
#
# Closed-form maximum likelihood standard error for sample coefficient
# alpha, derived under multivariate normality of the items by van Zyl,
# Neudecker, & Nel (2000):
#
#     Var(alpha_hat) = .alpha_vanzyl_avar(S) / (N - 1),
#
# with S the item covariance matrix (see .alpha_vanzyl_avar() for the
# expression and for why it is not evaluated on the correlation matrix).
# J is kept in the signature for the callers; the expression takes it from
# S. Wald CI on the reliability scale unless logistic = TRUE, in which case
# the CI is computed on the logit scale and back-transformed. The critical
# value is the normal quantile unless wald_df gives the degrees of freedom
# of a t reference (.wald_crit(), R/missing_information_internals.R).
.ci_alpha_ml <- function(estimate, S, N, J, conf_level, logistic = FALSE,
                         wald_df = Inf,
                         two_rows = .relia_two_rows_alternative()) {
  if (N < 3) {
    .dmar_refuse_argument(c(
      "The \"ml\" confidence interval requires 'N' of at least 3: the ",
      "covariance matrix of two rows has rank one (every pair of items ",
      "correlates +1 or -1, whatever the items), and the asymptotic ",
      "variance of alpha evaluated at it is zero for any data, so the ",
      "interval would have no width. ", two_rows), arg = "N")
  }
  variance <- .alpha_vanzyl_avar(S) / (N - 1)
  se <- sqrt(variance)
  crit <- .wald_crit(conf_level, wald_df)

  if (logistic) {
    lg <- .logistic_ci(estimate, se, crit)
    list(se = se, lower = lg$lower, upper = lg$upper)
  } else {
    list(se    = se,
         lower = estimate - crit * se,
         upper = estimate + crit * se)
  }
}


# .ci_alpha_adf(estimate, data, conf_level, logistic = FALSE,
#               wald_df = Inf)
#
# Asymptotic distribution-free (ADF) standard error of coefficient alpha
# (Maydeu-Olivares, Coffman, & Hartmann, 2007). Unlike .ci_alpha_ml,
# this method does not assume multivariate normality of the items.
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
# Wald CI on the reliability scale unless logistic = TRUE, with the
# normal quantile unless wald_df gives a t reference.
.ci_alpha_adf <- function(estimate, data, conf_level, logistic = FALSE,
                          wald_df = Inf,
                          two_rows = .relia_two_rows_alternative()) {
  d <- as.matrix(stats::na.omit(data))
  N <- nrow(d); J <- ncol(d)
  # From two rows each case's centered cross products are the same matrix,
  # half the covariance matrix, so every deviation s_i - vech(S) is a
  # multiple of vech(S). Alpha is homogeneous of degree zero in S, so its
  # gradient is orthogonal to vech(S) (Euler's theorem), and the variance
  # below, a sum of squares of gradient times deviation, is zero for any
  # data. Through version 1.0.1 it was computed: rounding error of either
  # sign under that zero gave NaN (the square root of a negative number)
  # or a standard error near zero.
  if (N < 3L) {
    .dmar_refuse_argument(c(
      "The \"adf\" confidence interval requires 'N' of at least 3: ",
      "from two rows each case's centered cross products are ",
      "proportional to the covariance matrix, and alpha does not change ",
      "when that matrix is rescaled, so the distribution free variance ",
      "of alpha is zero for any data (a computed negative value or NaN ",
      "is rounding error), and the interval would have no width. ",
      two_rows), arg = "N")
  }
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
  variance <- accumulator / (as.numeric(N) * (N - 1))
  se <- sqrt(variance)
  crit <- .wald_crit(conf_level, wald_df)

  if (logistic) {
    lg <- .logistic_ci(estimate, se, crit)
    list(se = se, lower = lg$lower, upper = lg$upper)
  } else {
    list(se    = se,
         lower = estimate - crit * se,
         upper = estimate + crit * se)
  }
}


# .ci_omega_delta(fit, conf_level, logistic = FALSE, wald_df = Inf)
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
#   wald_df    : degrees of freedom of the t reference distribution;
#                Inf (the default) is the normal one.
.ci_omega_delta <- function(fit, conf_level, logistic = FALSE,
                            wald_df = Inf) {
  if (!fit$converged) {
    return(list(se = NA_real_, lower = NA_real_, upper = NA_real_))
  }
  crit <- .wald_crit(conf_level, wald_df)
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
#   crit     : critical value (e.g., qnorm(0.975) for a 95% CI, or the t
#              quantile .wald_crit() gives when wald_df is finite, which
#              is then applied on the logit scale).
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

# .bootstrap_ci(data, point_fn, B, conf_level, kind, seed,
#               ncpus = 1L, improper = "keep", stat_names = NULL,
#               improper_note = NULL, two_rows = NULL,
#               undetermined_note = NULL, progress = NULL, fun = NULL)
#
# Nonparametric bootstrap confidence interval for a reliability
# coefficient: the one engine behind every bootstrap interval of the
# reliability family, and the engine other internal callers build on
# (its interface, documented here, is kept stable for them). Resamples
# rows of `data` with replacement B times, applies `point_fn` to each
# resample, and constructs a CI of the requested kind:
#
#   "bootstrap_se"          : Wald CI with SE = sd of bootstrap estimates.
#   "bootstrap_se_logistic" : same SE, but built on the logit scale.
#   "percentile"            : empirical 100*(1-C)/2 and 100*(1-(1-C)/2)
#                             quantiles of the bootstrap estimates.
#   "bca"                   : bias-corrected and accelerated bootstrap
#                             (boot::boot.ci with type = "bca", with the
#                             influence values computed here; see below).
#
# The point function. point_fn(d) receives a resample (rows of `data`,
# a data frame or matrix as `data` is) and returns a single number:
#
#   a finite number                    : the estimate from a proper
#                                        solution;
#   a finite number with a dmar_reason : the estimate from an improper
#     attribute of "improper"            (Heywood) solution that still
#                                        yields a coefficient, used or
#                                        set aside by the `improper`
#                                        rule;
#   NA_real_ with dmar_reason          : an improper solution that yields
#     "improper"                         no coefficient (categorical
#                                        omega, which is not defined for
#                                        it; a single-factor solution
#                                        whose error variances sum below
#                                        zero), set aside under either
#                                        rule;
#   NA_real_ with dmar_reason          : three items with no single-factor
#     "no_solution"                      solution (.triad_no_solution()),
#                                        set aside under either rule and
#                                        counted as improper and also on
#                                        its own;
#   NA_real_ with dmar_reason          : the resample left a polychoric
#     "undetermined_polychoric"          correlation undetermined, which
#                                        the observed denominator of
#                                        categorical omega needs;
#   NA_real_ (or NaN, Inf) otherwise   : no estimate, as from a fit that
#                                        did not converge;
#   NA_real_ with dmar_reason          : a fit that did not converge,
#     "not_converged"                    counted under that name (a
#                                        point function that does not
#                                        name it reports the fit as
#                                        non-finite);
#   an error                           : caught and counted.
#
# The single-factor point functions produce this form through
# .omega_point_value(), the categorical one through .omega_c_from_data().
#
# Several statistics from one set of replications. A point function may
# instead return a numeric vector of m statistics, with `stat_names`
# naming them (its length is m; NULL, the default, means one statistic),
# as ci_bifactor_indices() does for the bifactor indices. The
# dmar_reason attribute then describes the replication as a whole, and a
# replication with any statistic that is not finite is set aside whole,
# so every interval rests on the same replications and one count
# describes them all. Each interval is formed from its own column of
# boot's replicate matrix (index = j), and the jackknife behind a BCa
# acceleration runs once for all of them. With one statistic every
# result is as it was before `stat_names` existed. One exception: a
# statistic that is not defined for some replications, as the maximal
# reliability of a bifactor model is not for an improper solution,
# which ci_bifactor_indices() keeps for its other indices under
# improper = "keep". The point function marks it with the attribute
# dmar_undefined, a logical vector of length m that is TRUE for each
# statistic it does not define; such a statistic is NA in the
# replication, which is used or set aside on the others alone, and its
# interval rests on the replications that define it (an interval needs
# two of them, and the BCa regression estimate of the acceleration as
# many as there are cases, counted for that statistic). The first
# statistic is never marked, since the counts read the first column.
# A point function must be deterministic given its resample (the
# replications that raised an error are refit once, on the same
# processes, to recover their messages), and with ncpus > 1 it runs in
# forked (multicore) or socket (snow) workers, so it may not rely on
# side effects in the calling session.
#
# Every replication is counted. The statistic handed to boot::boot()
# returns three numbers per replication: the value the interval uses (the
# estimate, or NA when the replication is set aside), a status code (1
# proper, 2 improper, 3 error, 4 non_finite, 5 undetermined, 6
# not_converged; see .BOOT_STATUS), and the estimate itself (NA where
# there is none); with m statistics, the m values, the status code, and
# the m estimates, 2m + 1 numbers in all. The
# counts are read from boot's replicate matrix rather than from counters
# in an environment, which a forked worker updates in its own copy and
# the parent never sees (204 failures counted serially came back as 0
# under parallel = "multicore"). The interval is formed from the first
# column (boot.ci() and empinf() with index = 1), which equals the single
# column the engine handed boot through DMAR 1.0.1: boot draws its whole
# index array in the parent before any replication runs, so neither the
# extra columns nor the number of workers changes a draw or a limit.
# Through 1.0.1 boot::boot() was also called without a `parallel`
# argument, so a session option boot.parallel = "multicore" ran it in
# forked workers and lost every count; the engine now passes `parallel`
# itself.
#
# A replication is set aside, with its reason recorded, when point_fn
# raises an error, returns no finite estimate, or returns an improper
# estimate that the rule sets aside, and the interval is formed from the
# replications that return a value: the standard error is their standard
# deviation and the percentile and BCa limits are their quantiles. A
# single warning states the counts whenever a replication was set aside.
# Since 1.0.1; in 1.0.0 the failed replications were dropped without a
# count.
#
# Improper replications. Under improper = "keep" (the default) the
# estimate of an improper solution that still yields a coefficient is
# used like any other; under "drop" it is set aside and counted. Through
# DMAR 1.0.1 such estimates were kept and not counted. Keeping them is
# the default for the single-factor coefficients because, in a small
# coverage check (model implied one-factor omega, N = 60, 120 data sets
# per condition, B = 300; loadings .95, .70, .60, .50, .40), percentile
# intervals that kept them covered at .958 while intervals that dropped
# them covered at .933 (a Monte Carlo standard error near .02, so
# suggestive rather than decisive): setting aside the replications in
# which an error variance fell below zero (a standardized loading past
# one; the fits are on the items' own metric, where a loading above one
# is routine) conditions the bootstrap distribution on admissibility and
# trims its upper tail. Yuan and Hayashi (2003, Section 6) make the
# parallel argument for test statistics, counting a replication that
# fails to converge because its resample lies far from the model as an
# extreme value rather than discarding it. That check kept every
# improper replication; the engine sets aside, under either rule, those
# whose solution yields no coefficient (.omega_fit_verdict()), which
# were rare in that condition (the largest gap between a kept and a
# dropped limit was 1.5 on the unclamped scale). Two kinds yield none. A
# solution whose error variances sum below zero has run away along a
# ridge of the likelihood, and the omega it implies exceeds one: like
# the improper replications kept, it lies beyond the top of the
# distribution. Three items whose correlations have a product that is
# not positive have no single-factor solution at all
# (.triad_no_solution()), and where such a resample belongs depends on
# the scale it is read on. On the estimator's own scale it continues the
# top of the distribution: as a resample approaches that boundary, one
# loading, and omega with it, grows without bound, and an optimizer left
# to run stops at an omega in the hundreds. On the data's scale it points
# low: one of its correlations is negative or zero, and its coefficient
# alpha is typically among the lowest of the replications. Neither
# placement is principled (Yuan and Hayashi's extreme value belongs to a
# misfit statistic, which has a direction, and omega has none at that
# boundary), so such a replication is counted, also on its own as
# `no_solution`, rather than used. Where such resamples are common the
# percentile lower limit is too high: in a coverage check made for the
# 1.1.0 review (the saturated three-item solution in closed form, N =
# 60, 500 data sets per condition, B = 999), it exceeded the population
# coefficient in 11.6 percent of samples with loadings .70, .60, .30 and
# in 13.1 percent with loadings .775, .645, .155, against the nominal 2.5
# percent, and the upper limit never fell short of it. Ranking these
# replications below every other instead brought the rate to 3.9 and 2.6
# percent only by putting the median lower limit at zero, which is no
# lower limit at all. Much of the excess comes from the estimator
# itself, which runs high near that boundary and is refused outright on
# the samples beyond it: in the same two conditions (4,000 samples each)
# 15 and 31 percent of samples were refused, and the estimate exceeded
# the population coefficient in 72 and 77 percent of the rest. No
# placement of these replications repairs that, and the caution says
# that the lower limit is then likely too high.
#
# The caution. When the replications whose treatment is in question,
# those set aside for any reason and the improper ones kept, number at
# least the tail area (1 - C)/2 of the B, they can decide where a limit
# falls, so the interval depends on how they were treated; the warning
# and the printed line say so and the `caution` element is TRUE. (When
# about half of the replications are improper, as in the second
# condition of the coverage check above, neither rule gives an interval
# worth reporting, and what the result can state is the count.)
#
# Fewer than half returned. An interval is formed
# from the replications that return a value, so it conditions on
# whatever set the others aside. When fewer than half return one, the
# rest are a minority of the resamples, chosen by whether they returned
# a value, and their spread describes those resamples, not the sampling
# distribution of the coefficient; no interval
# is reported then (the standard error and the limits are NA), and the
# warning and the printed line say so. This is the criterion the point
# estimate of categorical omega under multiple imputation applies to its
# completed data sets (Lim & Cheung, 2022; .omega_c_mi_check()), applied
# to the replications. Categorical omega asks for it: its replications
# are set aside for undetermined polychoric correlations, unidentified
# loadings, and improper solutions, all of which a resample meets far
# more often than the sample, and on Green and Yang's mixed-skew design
# at N = 300 a percentile bootstrap of the model implied denominator set
# aside more than half of the replications in 15 of 45 samples whose
# point estimate was identified and proper (median 42%, at most 90%,
# B = 50). The other callers previously kept two returned values, the
# rule through 1.0.1. Since 1.1.0 the half-return safeguard applies to
# every caller by default. AVE retains its additional minimum of 100
# returned estimates. Neither floor guarantees interval coverage.
#
# The count rule. The number .dmar_boot_required() gives for the kind and
# the level (10 beyond each limit for "percentile", 50 beyond each nominal
# limit for "bca", 100 for the two standard error kinds) is held to the
# number of replications REQUESTED, B (.dmar_boot_count_met() in
# R/dmar_bootstrap.R, where the rule is stated once), and so is the count
# of 10 where boot.ci() places a BCa limit. Replications set aside from a
# request that met the count do not withhold the interval by their
# number; the share rules above decide (fewer than half returned, more
# than half with too few distinct cases, the caution), and the warning
# says how many were set aside and, when those that returned a value are
# fewer than the number required, to request more. Through the first
# versions of 1.1.0 the count was applied to the replications that
# returned a value, so one replication set aside at exactly the required
# B withheld every interval (15 of them in ci_bifactor_indices(B = 400)
# with 399 usable). Under options(dmar.bootstrap_exploratory = TRUE) a
# request below the count runs and is marked by a warning of class
# dmar_bootstrap_exploratory. The warning of a run with replications set
# aside has class dmar_bootstrap_replications_dropped when every interval
# was formed and dmar_bootstrap_insufficient when one was withheld.
#
# More than half with too few distinct cases. One cause of setting a
# replication aside is held to the same criterion for every caller: the
# rule on the number of cases (.relia_check_resample(), which raises an
# error of class dmar_too_few_cases, counted in `too_few_cases`). A
# resample repeats rows, and the rule keeps only the resamples that repeat
# few enough of them, so the replications it keeps are chosen by their
# composition, not by a sampling accident of the coefficient: they hold
# more distinct cases than a resample typically does, and their spread is
# smaller than that of the sampling distribution. At the bound itself the
# selection is complete: when the data hold the fewest distinct cases the
# rule allows (p + 1 for a rule asking for more than p), a resample that
# passes holds every distinct case, so without repeated rows it is the data
# reordered and gives the estimate itself (up to the convergence of an EM
# estimate), and a percentile interval has no width, and can exclude the
# estimate (five tests of six pupils with two values missing under FIML,
# the percentile interval [0.04347722, 0.04347727] about the estimate
# 0.04347698, from the 3 of 200 resamples that held all six pupils, 6!/6^6
# = 1.5 percent of them; in the first versions of 1.1.0, which analyzed
# the other resamples on singular EM estimates, [0, 0.898]). When more than
# half of the B are set aside by the rule, no interval is reported (the
# standard error and the limits are NA), and
# the warning and the printed line say why. Without repeated rows that
# happens below about 1.65 p cases: the share of resamples of n cases with
# more than p distinct ones, sum over k > p of choose(n, k) S(n, k) k! /
# n^n with S the Stirling numbers of the second kind, first reaches one
# half at n = 6 for p = 3, 9 for p = 5, 14 for p = 8, 17 for p = 10, and
# 33 for p = 20, and it is at most 2/9 at n = p + 1.
#
# Parallel replications. With ncpus > 1 the replications are spread over
# forked workers (boot's parallel = "multicore") where the platform
# allows it, and over a socket cluster (parallel = "snow") on Windows,
# whose workers load the installed DMAR; the jackknife behind the BCa
# acceleration uses the same workers. Because boot draws its index array
# in the parent and every point function is deterministic, a given seed
# returns the same interval, counts, and replicates for any ncpus.
#
# The BCa acceleration. boot::boot.ci estimates the acceleration from
# empirical influence values L through a = sum(L^3) / (6 sum(L^2)^1.5),
# and chooses how to get L by its own rule (boot::empinf): when at least
# n replications are available it regresses the replicates on the
# resampling frequencies, which costs no refit, and otherwise it runs the
# jackknife, one refit of `point_fn` per case. Under that rule a single
# jackknife refit that fails made every influence value NA in boot's
# own jackknife (the values are centered on the full-data estimate, but
# boot propagates the NA through the whole vector), so `a` was NA and
# boot.ci stopped with "estimated adjustment 'a' is NA"; the 2026-10-14
# report from the Psychological Methods task hit this on every BFI
# construct (n about 1,355, B = 1,000). The repair keeps boot's rule for
# which estimate to use, so a run on which nothing fails returns
# boot's interval to machine precision, and computes the influence
# values here:
#
#   - regression (returned replications >= n): boot::empinf(type = "reg")
#     on the replications that returned a value; if any influence value
#     is not finite (a rank-deficient regression), the jackknife is run
#     instead;
#   - jackknife (otherwise): L_i = (n - 1) (t0 - t_(i)) exactly as
#     boot's usual.jack() forms it, over the leave-one-out replicates
#     that returned a value (improper ones under the same rule as the
#     bootstrap), with the rest counted. The (n - 1) scale cancels in
#     `a`; the centering on t0 rather than on the jackknife mean is
#     boot's convention and is kept so a clean run reproduces boot's
#     acceleration.
#
# When `a` still cannot be estimated as a finite number (fewer than two
# jackknife replicates returned a value, the floor for a skewness, or
# the influence values have no spread), or when the bias correction
# cannot be (every returned replicate on one side of the estimate,
# boot's "estimated adjustment 'w' is infinite"), the percentile
# interval is returned in its place, `ci_method` says so, and the
# warning says why. A BCa request never stops on either adjustment.
#
# Progress and cancellation (since 1.1.0). A caller may pass a callback,
# `progress`, which an application uses to show how far a long bootstrap
# has run and to let its user stop it (a Shiny app, for example, which
# runs the bootstrap in a process its user cannot interrupt directly). It
# is called as progress(done, total), with the replications attempted so
# far and B: after each replication when the replications run in this
# process, and with ncpus > 1 after each block of replications the workers
# complete together (.boot_block_size(): about B/50, rounded up to a
# multiple of ncpus). When the BCa acceleration needs the jackknife (for a
# statistic whose interval is formed and whose returned replicates are not
# all equal, when fewer of its replications than cases returned a value or
# boot::empinf()'s regression fails or gives influence values that are not
# all finite; see `interval` below; it runs at most once in a call), the
# callback is called again after each leave-one-out refit (each block of
# them, of the size the same rule gives for n refits), with done and total
# both B, so a stop is honored there too. If it returns FALSE itself (as
# identical() tests it), the engine stops before the next replication or
# block with an error of class dmar_bootstrap_cancelled carrying done and
# total; any other value is ignored, and an error inside the callback
# stops the call with an error of class dmar_progress_error, the
# callback's condition its parent (.boot_progress()). With a callback the
# engine draws the resamples through boot::boot() exactly as without one
# (boot's own call, with a statistic that computes nothing), reads them
# back with boot::boot.array(), and runs the replications itself, in order
# (.boot_with_progress()); since boot draws every resample before any
# replication and a point function is deterministic, the replicate matrix,
# and so every interval, count, and draw, is identical to the one
# boot::boot() returns. The callback runs outside the random number stream
# (withr::with_preserve_seed()), so even one that draws random numbers
# changes neither the result nor the caller's stream. With progress = NULL
# nothing of this runs, and the engine calls boot::boot() as it always
# has.
#
# Seed handling: if seed is not NULL, .dmar_local_seed() sets it before
# resampling so that runs are reproducible, and the caller's RNG state is
# restored when this function exits (the package-wide seed discipline;
# see R/dmar_seed.R and the roxygen of every caller). This function is
# where the random draws begin and end, so it is the one that owns the
# seed.
#
# Inputs:
#   data       : data frame or matrix passed row-wise to boot::boot.
#   point_fn   : function(d) returning a single number, as above.
#   B          : number of bootstrap replications.
#   conf_level : confidence level.
#   kind       : one of "bootstrap_se", "bootstrap_se_logistic",
#                "percentile", "bca".
#   seed       : integer seed or NULL.
#   ncpus      : number of processes the replications are spread over
#                (1, the default, runs them in this process).
#   improper   : "keep" (default) or "drop", the treatment of the
#                estimates of improper solutions that yield a
#                coefficient.
#   stat_names : NULL (one statistic) or the names of the m statistics
#                a point function returns.
#   improper_note : NULL, or the sentence the warning gives about the
#                replications set aside as improper, for a caller whose
#                solutions are not single-factor fits (the default
#                sentence speaks of an error variance below zero).
#   two_rows   : NULL, or the sentence the refusal of two rows ends with,
#                naming the interval the caller has for them
#                (.relia_two_rows_alternative()); a caller that has none
#                passes NULL, so the refusal never names an interval the
#                call would refuse.
#   undetermined_note : NULL, or the sentence the warning gives about the
#                replications set aside because a polychoric correlation
#                was undetermined, which depends on categorical omega's
#                denominator (.omega_c_undetermined_note(); NULL gives
#                the observed denominator's).
#   progress   : NULL (the default), or a function(done, total), as above.
#   fun        : NULL, or the name of the exported function that called
#                the engine ("reliability_alpha"), without parentheses.
#                Every warning the engine gives then begins with
#                "<fun>(): ", as the warnings of the rule do
#                (.dmar_boot_message()), and carries the name in the
#                condition's `fun` field. Every caller in the package
#                passes it.
# Returns:
#   A list with `se`, `lower`, `upper`, `ci_method` (the kind requested,
#   or the percentile fallback label when a BCa adjustment could not be
#   estimated), each of length m (one for a single statistic), and
#   `bootstrap`, which .relia_result() attaches as the `bootstrap`
#   attribute. Its elements, in this order (attempted through improper,
#   acceleration, acceleration_method, and jackknife have had these
#   meanings since 1.0.1; improper_kept, not_converged, too_few_cases,
#   improper_rule, caution, require_half, ncpus, t, status, and the
#   jackknife's improper_kept, not_converged, and too_few_cases were added
#   in 1.1.0, so the attribute only gained elements):
#     attempted, returned, failed : B, the replications whose value the
#                                   interval uses, and the difference;
#     reasons                     : c(error = , non_finite = ), how the
#                                   failed replications were set aside
#                                   (non_finite: no value the interval
#                                   could use, including the causes
#                                   named in `undetermined` and
#                                   `improper`);
#     error_messages              : the distinct error messages with
#                                   their counts (empty when none);
#     undetermined                : the replications set aside because
#                                   the resample left a polychoric
#                                   correlation undetermined (0 for
#                                   every coefficient but categorical
#                                   omega with the observed denominator);
#     improper                    : the replications from an improper
#                                   solution that were set aside (a
#                                   subset of reasons[["non_finite"]]);
#     improper_kept               : the replications from an improper
#                                   solution whose estimate was used
#                                   (a subset of `returned`; 0 under
#                                   improper = "drop");
#     no_solution                 : the replications set aside because
#                                   three items had no single-factor
#                                   solution (a subset of `improper`;
#                                   since 1.1.0);
#     not_converged               : the replications set aside because
#                                   the fit did not converge, as the
#                                   point function named them (a subset
#                                   of reasons[["non_finite"]]);
#     too_few_cases               : the replications set aside because
#                                   the resample had too few distinct
#                                   cases (an error of class
#                                   dmar_too_few_cases; a subset of
#                                   reasons[["error"]]); with more than
#                                   half of B, no interval is formed;
#     improper_rule               : "keep" or "drop";
#     caution                     : TRUE when failed + improper_kept is
#                                   at least (1 - conf_level)/2 of B;
#     require_half                : TRUE: no interval is formed from
#                                   fewer than half of the B (always
#                                   TRUE; a documented element, kept);
#     required                    : the fewest replications a request for
#                                   an interval of this kind and level
#                                   may ask for (.dmar_boot_required());
#     bootstrap_status            : "complete", "dropped", "exploratory",
#                                   "insufficient_count", or
#                                   "insufficient_share" (the verdict of
#                                   the rule; `status`, below, is the
#                                   fate of each replication);
#     B_used                      : the smallest count any reported
#                                   limit rests on;
#     acceleration                : the BCa acceleration used (NA when
#                                   no BCa interval was formed), one
#                                   per statistic;
#     acceleration_method         : "regression" or "jackknife" (NA when
#                                   no BCa interval was formed), one per
#                                   statistic;
#     jackknife                   : the same attempted / returned /
#                                   failed / reasons / error_messages /
#                                   undetermined / improper /
#                                   improper_kept / no_solution /
#                                   not_converged / too_few_cases for
#                                   the leave-one-out refits behind the
#                                   acceleration (all zero when no
#                                   jackknife ran);
#     ncpus                       : the number of processes used;
#     t                           : the estimate of every replication,
#                                   in order (NA where there is none;
#                                   an improper estimate set aside under
#                                   "drop" is kept here, so the interval
#                                   under the other rule can be formed
#                                   without refitting); with m
#                                   statistics, a B x m matrix whose
#                                   columns are named by stat_names;
#     status                      : a factor of the fate of every
#                                   replication, with levels "proper",
#                                   "improper", "error", "non_finite",
#                                   "undetermined", and "not_converged".
#                                   The interval uses the replications
#                                   whose status is "proper", and under
#                                   "keep" those whose status is
#                                   "improper" and whose t is finite;
#   and `replicate_mean`, beside `bootstrap` rather than in it (since
#   1.1.0): the mean of the replicates whose value an interval uses, one
#   per statistic (NA when no interval is formed), from which a caller
#   can estimate the bias of the sample coefficient (reliability_H()
#   does).
.bootstrap_ci <- function(data, point_fn, B, conf_level, kind, seed,
                          ncpus = 1L, improper = "keep",
                          stat_names = NULL, improper_note = NULL,
                          two_rows = NULL, undetermined_note = NULL,
                          progress = NULL, fun = NULL) {
  if (!requireNamespace("boot", quietly = TRUE)) {
    .dmar_refuse(c(
      "Bootstrap CIs require the 'boot' package. ",
      "Install with install.packages(\"boot\")."),
      "package_required", package = "boot")
  }
  if (!kind %in% c("bootstrap_se", "bootstrap_se_logistic",
                   "percentile", "bca")) {
    .dmar_refuse_argument(c("Unknown bootstrap kind: ", kind), arg = NULL)
  }
  improper <- match.arg(improper, c("keep", "drop"))
  ncpus <- .check_ncpus(ncpus)
  progress <- .check_progress(progress)
  m <- max(1L, length(stat_names))
  n <- nrow(data)
  # Two rows leave nothing to resample: each replication is one row drawn
  # twice, which has no variance, or both rows, which give the estimate
  # itself, so every interval would have no width.
  if (n < 3L) {
    .dmar_refuse_argument(c(
      "A bootstrap interval requires at least 3 rows of data: from ",
      "two, each resample is one row drawn twice, which has no ",
      "variance, or both rows, which give the estimate itself, so the ",
      "interval would have no width.",
      if (!is.null(two_rows)) paste0(" ", two_rows)), arg = NULL)
  }
  .dmar_local_seed(seed)

  par_type <- .boot_parallel_type(ncpus)
  cl <- NULL
  if (par_type == "snow") {
    cl <- parallel::makePSOCKcluster(ncpus)
    on.exit(parallel::stopCluster(cl), add = TRUE)
  }

  statistic <- .boot_statistic(point_fn, improper, m)
  boot_out <- if (is.null(progress)) {
    boot::boot(data = data, statistic = statistic, R = B, stype = "i",
               parallel = par_type, ncpus = ncpus, cl = cl)
  } else {
    .boot_with_progress(data, statistic, B, m, par_type, ncpus, cl,
                        progress)
  }
  tt <- boot_out$t
  code_col <- m + 1L
  counts <- .boot_tally(tt, B)
  err <- which(tt[, code_col] == 3L)
  if (length(err)) {
    idx <- boot::boot.array(boot_out, indices = TRUE)
    refit <- .boot_refit_errors(point_fn, data,
                                lapply(err, function(r) idx[r, ]),
                                par_type, ncpus, cl)
    counts$error_messages <- refit$messages
    counts$too_few_cases <- refit$too_few_cases
  }
  # The jackknife behind a BCa acceleration runs at most once, the first
  # time a statistic needs it, and every statistic reads its leave-one-out
  # replicates; it is kept in a small environment of its own.
  jk <- new.env(parent = emptyenv())
  jk$t <- NULL
  jk$tally <- .boot_tally(matrix(numeric(0), 0L, 2L * m + 1L), 0L)
  t0 <- as.numeric(boot_out$t0[seq_len(m)])
  returned <- counts$returned
  crit <- stats::qnorm(1 - (1 - conf_level) / 2)
  # Fewer than half of the B returning a value leaves no interval (see
  # above), as fewer than two always does, and so does more than half of
  # the B set aside for too few distinct cases (see above).
  short_of_cases <- returned >= 2L && counts$too_few_cases > B / 2
  short_of_half <- returned >= 2L && returned < B / 2 && !short_of_cases
  # The count rule of ?dmar_bootstrap: an interval is formed only when
  # the number of replications requested, B, is the number its kind and
  # level require (10 beyond each limit for "percentile", 50 beyond each
  # nominal limit for "bca", 100 for the two standard error kinds),
  # unless the exploratory option is set, under which the result is
  # marked. The count is not applied to the replications that returned a
  # value (.dmar_boot_count_met()); the share rules above are.
  boot_kind <- .boot_method_kind(kind)
  required <- .dmar_boot_required(boot_kind, tail = (1 - conf_level) / 2)
  count_met <- .dmar_boot_count_met(B, required)
  short_of_count <- returned >= 2L && !short_of_half && !short_of_cases &&
    !count_met
  enough <- returned >= 2L && !short_of_half && !short_of_cases &&
    !short_of_count

  # The interval for statistic j: column j of boot's replicate matrix,
  # over the replications that returned a value for it (all those used,
  # unless the point function marked the statistic not defined for some).
  interval <- function(j) {
    t_boot <- tt[, j]
    est <- t0[j]
    # The replications that returned a value for this statistic, which a
    # point function can mark not defined for some; the interval needs
    # two of them, besides the call's own rule (`enough`).
    returned <- sum(is.finite(t_boot))
    formed <- enough && returned >= 2L && returned >= B / 2
    se <- if (formed) stats::sd(t_boot, na.rm = TRUE) else NA_real_
    out <- list(se = se, lower = NA_real_, upper = NA_real_, method = kind,
                acceleration = NA_real_,
                acceleration_method = NA_character_,
                fallback = NULL,  # why the percentile interval replaced BCa
                constant = FALSE,  # every returned replicate the same
                bca_withheld = NULL)  # the verdict of the position rule
    if (!formed) {
      # Nothing to form an interval from, or fewer than half returned a
      # value; the warning below says so.
      return(out)
    }
    out$constant <- .boot_constant(t_boot)
    if (out$constant && kind %in% c("bootstrap_se",
                                    "bootstrap_se_logistic")) {
      # A standard deviation of zero gives the estimate at both limits, on
      # either scale. The logit interval is formed apart from it, since
      # at a value of one (identical items) the logit is infinite and
      # .logistic_ci() returns NA limits, which came back through version
      # 1.0.1, and in the first versions of 1.1.0 under a warning that the
      # interval had no width. Rounding decided the case in 1.0.1: where
      # the estimate and the replicates came out just below one the logit
      # was finite and the limits one, and otherwise they were NA, so
      # KR-20 of three copies of a 0/1 item returned one with some data
      # and NA with others.
      out$lower <- est
      out$upper <- est
    } else if (kind == "bootstrap_se") {
      out$lower <- est - crit * out$se
      out$upper <- est + crit * out$se
    } else if (kind == "bootstrap_se_logistic") {
      lg <- .logistic_ci(est, out$se, crit)
      out$lower <- lg$lower
      out$upper <- lg$upper
    } else if (kind == "percentile") {
      pc <- .boot_percentile(boot_out, conf_level, j)
      out$lower <- pc[1L]
      out$upper <- pc[2L]
    } else if (kind == "bca" && out$constant) {
      # With every returned replicate equal there is no bias correction to
      # estimate (no replicate falls below the estimate, so w is infinite)
      # and no spread for an acceleration; boot.ci() would print a note
      # and return NULL. The percentile interval, the common value at both
      # ends, is reported in its place.
      out$fallback <- paste0("every returned replicate has the same ",
                             "value, so the bias correction is not ",
                             "defined")
    } else if (kind == "bca") {
      L <- NULL
      if (returned >= n) {
        # boot's regression estimate of the influence values, the one
        # boot.ci would use itself with this many finite replicates.
        L <- tryCatch(boot::empinf(boot_out, type = "reg", index = j),
                      error = function(e) NULL)
        if (!is.null(L) && !all(is.finite(L))) L <- NULL
        if (!is.null(L)) out$acceleration_method <- "regression"
      }
      if (is.null(L)) {
        # The jackknife, one leave-one-out refit per case, with every
        # refit that returns no usable value counted instead of poisoning
        # the vector.
        if (is.null(jk$t)) {
          jt <- .boot_jackknife(point_fn, data, improper, par_type, ncpus,
                                cl, m, progress = progress, B = B)
          tally <- .boot_tally(jt, n)
          jerr <- which(jt[, code_col] == 3L)
          if (length(jerr)) {
            jrefit <- .boot_refit_errors(
              point_fn, data, lapply(jerr, function(k) seq_len(n)[-k]),
              par_type, ncpus, cl)
            tally$error_messages <- jrefit$messages
            tally$too_few_cases <- jrefit$too_few_cases
          }
          jk$t <- jt
          jk$tally <- tally
        }
        fin <- is.finite(jk$t[, j])
        L <- (n - 1) * (est - jk$t[fin, j])
        out$acceleration_method <- "jackknife"
      }
      # The acceleration is a skewness of the influence values, so it
      # needs at least two of them with some spread. boot centers the
      # jackknife values on t0 rather than on their own mean, which is
      # why a single returned value would otherwise give the finite but
      # meaningless extreme of 1/6 in magnitude; the floor is imposed
      # here.
      a <- if (length(L) >= 2L) sum(L^3) / (6 * sum(L^2)^1.5) else NA_real_
      if (is.finite(a)) {
        bci <- tryCatch(
          .boot_ci_endpoints_quiet(
            boot::boot.ci(boot_out, conf = conf_level, type = "bca", L = L,
                          index = j)$bca,
            !is.na(required) && returned < required),
          error = function(e) e)
        if (inherits(bci, "error") || is.null(bci)) {
          # The remaining way bca.ci stops: the bias correction w is
          # infinite because every returned replicate sits on one side
          # of the estimate. Any other failure is named by its own
          # message.
          msg <- if (inherits(bci, "error")) conditionMessage(bci) else ""
          out$fallback <- if (!nzchar(msg) || grepl("'w' is infinite", msg,
                                                    fixed = TRUE)) {
            paste0("BCa bias correction could not be estimated: every ",
                   "returned replicate fell on one side of the estimate")
          } else {
            paste0("boot.ci() stopped: ", trimws(gsub("[[:space:]]+", " ",
                                                      msg)))
          }
          out$acceleration_method <- NA_character_
        } else {
          # The count of 10 beyond each limit, applied where the bias
          # correction and the acceleration place the limits. boot.ci()
          # returns the two order-statistic positions on the (R + 1)
          # scale, R the replicates that returned a value.
          at <- .dmar_boot_bca_positions(
            .dmar_boot_settle(B, returned, kind = "bca",
                              tail = (1 - conf_level) / 2,
                              conf_level = conf_level, fun = "",
                              quantity = "the coefficient",
                              dropped_because = "returned no value",
                              signal = FALSE),
            c(bci[2], bci[3]) / (returned + 1))
          out$acceleration <- a
          if (isTRUE(at$bca_withheld)) {
            out$bca_withheld <- at
          } else {
            out$lower <- bci[4]
            out$upper <- bci[5]
          }
        }
      } else {
        jfail <- jk$tally$failed
        out$fallback <- if (jfail > 0L) {
          sprintf("BCa acceleration could not be estimated: %d jackknife %s failed",
                  jfail, if (jfail == 1L) "replicate" else "replicates")
        } else {
          "BCa acceleration could not be estimated: the influence values have no spread"
        }
        out$acceleration_method <- NA_character_
      }
    }
    if (kind == "bca" && !is.null(out$fallback)) {
      pc <- .boot_percentile(boot_out, conf_level, j)
      out$lower <- pc[1L]
      out$upper <- pc[2L]
      out$method <- sprintf("percentile (%s)", out$fallback)
    }
    out
  }
  res <- lapply(seq_len(m), interval)
  jack <- jk$tally
  pick <- function(what, type) vapply(res, function(r) r[[what]], type)

  # The caution: the replications whose treatment is in question, those
  # set aside and the improper ones kept, reach the tail area beyond a
  # limit (.boot_caution()).
  in_question <- counts$failed + counts$improper_kept
  tail_area <- (1 - conf_level) / 2
  caution <- .boot_caution(in_question, B, conf_level)

  # One warning per call, stating every count the user should know.
  notes <- character()
  target_returned <- colSums(is.finite(tt[, seq_len(m), drop = FALSE]))
  target_short <- target_returned < 2L | target_returned < B / 2
  required_shown <- max(2L, ceiling(B / 2))
  if (enough && any(target_short)) {
    notes <- c(notes, paste0(
      "Intervals withheld for ", paste(stat_names[target_short], collapse = ", "),
      ": their defined estimates number ",
      paste(target_returned[target_short], collapse = ", "), " of ", B,
      ", below the ", required_shown, " (half of those requested) an ",
      "interval is reported from. Other targets retain intervals when eligible"))
  }
  # The sentence that states the count rule when it withholds.
  count_rule <- if (boot_kind == "se") {
    sprintf(paste0("a standard error taken from the replications is ",
                   "reported only from a request of at least %d"),
            required)
  } else {
    sprintf(paste0(
      "a %s%% %sinterval read from the quantiles of the replications is ",
      "reported only when at least 10 replications lie beyond each ",
      "limit%s, which takes a request of %d at this level"),
      .format_conf_pct(conf_level),
      if (boot_kind == "bca") "BCa " else "",
      if (boot_kind == "bca") paste0(
        " where the bias correction and the acceleration place it (50 ",
        "beyond each nominal limit)") else "", required)
  }
  if (counts$failed > 0L) {
    notes <- c(notes, sprintf(
      "%d of %d bootstrap replications %s set aside (%s)%s",
      counts$failed, B,
      if (counts$failed == 1L) "was" else "were",
      .bootstrap_reason_phrase(counts$reasons, counts$error_messages,
                               counts$undetermined, counts$improper,
                               counts$not_converged, counts$no_solution),
      if (enough) {
        sprintf("; the %s computed from the %d that returned a value%s%s",
                if (m == 1L) "interval is" else "eligible intervals are",
                returned,
                if (counts$improper_kept > 0L) {
                  sprintf(paste0(", %d of them from an improper (Heywood) ",
                                 "solution and kept under ",
                                 "improper = \"keep\""),
                          counts$improper_kept)
                } else {
                  ""
                },
                # The count is held to the request; the note says so and
                # says to request more where replications can fail.
                sub("\\.$", "", sub(
                  "^ The", "; the",
                  .dmar_boot_request_note(B, returned, required,
                                          "replications"))))
      } else if (short_of_cases) {
        sprintf(paste0(
          "; %d of them had too few distinct cases for the coefficient, ",
          "more than half of the %d, so no %s reported: the %d that ",
          "returned a value are the resamples that repeat the fewest rows, ",
          "and their spread understates the sampling variability of the ",
          "coefficient (when the data hold the fewest distinct cases the ",
          "rule allows, each of them holds every distinct case, and without ",
          "repeated rows it is the data reordered and gives the estimate ",
          "itself)"),
          counts$too_few_cases, B,
          if (m == 1L) "interval is" else "intervals are", returned)
      } else if (short_of_half) {
        sprintf(paste0(
          "; fewer than half returned a value, so no %s reported: the %d ",
          "that did are a minority of the resamples, chosen by whether ",
          "they returned a value, and their spread describes those ",
          "resamples, not the sampling distribution of the coefficient"),
          if (m == 1L) "interval is" else "intervals are", returned)
      } else if (short_of_count) {
        # B itself is below the requirement (a caller that did not hold
        # B to it).
        sprintf(paste0(
          "; the %d requested are fewer than the %d required at this ",
          "level, so no %s reported: %s. The point estimate is not ",
          "affected"),
          B, required,
          if (m == 1L) "interval is" else "intervals are", count_rule)
      } else {
        sprintf("; fewer than two returned a value, so no %s reported",
                if (m == 1L) "interval is" else "intervals are")
      }))
  } else if (short_of_count) {
    # Nothing was set aside, but B itself is below the requirement (a
    # caller that did not hold B to it).
    notes <- c(notes, sprintf(
      "Only %d bootstrap replications were attempted, fewer than the %d required at this level, so no %s reported: %s",
      B, required, if (m == 1L) "interval is" else "intervals are",
      count_rule))
  } else if (returned < 2L) {
    # No replication was set aside, but there were fewer than two to
    # begin with (a caller that does not hold B to at least 2).
    notes <- c(notes, sprintf(
      "Only %d bootstrap %s attempted, and an interval needs at least two that return a value, so no %s reported",
      B, if (B == 1L) "replication was" else "replications were",
      if (m == 1L) "interval is" else "intervals are"))
  }
  if (jack$failed > 0L) {
    notes <- c(notes, sprintf(
      "%d of %d jackknife replicates behind the BCa acceleration failed (%s)%s",
      jack$failed, jack$attempted,
      .bootstrap_reason_phrase(jack$reasons, jack$error_messages,
                               jack$undetermined, jack$improper,
                               jack$not_converged, jack$no_solution),
      if (any(pick("acceleration_method", character(1)) %in% "jackknife")) {
        sprintf("; the acceleration was estimated from the %d that returned a value",
                jack$returned)
      } else {
        ""
      }))
  }
  if (counts$undetermined > 0L) {
    notes <- c(notes, undetermined_note %||%
                 .omega_c_undetermined_note("observed"))
  }
  # The improper replications set aside: those with a negative error
  # variance (a Heywood case), and, named apart, those of three items with
  # no single-factor solution, for which no loading was estimated and
  # whose place in the distribution is a different question
  # (.BOOT_NO_SOLUTION).
  heywood <- counts$improper - counts$no_solution
  if (counts$improper > 0L && !is.null(improper_note)) {
    notes <- c(notes, improper_note)
  } else {
    if (heywood > 0L) {
      notes <- c(notes, paste0(
        "An improper replication is one in which an error variance fell ",
        "below zero (a standardized loading past one), so setting these ",
        "aside can trim the upper tail of the bootstrap distribution",
        if (improper == "keep") {
          paste0("; under improper = \"keep\" an improper replication is ",
                 "set aside only when its solution yields no coefficient ",
                 "(its error variances sum below zero)")
        } else {
          ""
        }))
    }
    if (counts$no_solution > 0L) {
      notes <- c(notes, paste0(
        "A replication with no single-factor solution is a resample of ",
        "three items whose correlations have a product that is not ",
        "positive, which no loadings reproduce, so it yields no ",
        "coefficient under either rule. Its correlations point to low ",
        "reliability, while the estimator runs high as a resample nears ",
        "that boundary, so it belongs on neither side of the bootstrap ",
        "distribution and is set aside"))
    }
  }
  # The caution speaks of the limits, so it is given only when an interval
  # was formed (`enough`: two or more returned a value, no more than half
  # had too few distinct cases, and at least half
  # returned a value); through the first versions of 1.1.0 it
  # followed "no interval is reported" when fewer than two did. It names an
  # improper solution only when replications from one were kept: analytic
  # alpha and KR-20 fit no model, and theirs are set aside.
  if (caution && enough) {
    notes <- c(notes, sprintf(paste0(
      "Caution: %d of %d replications (%s) %s, at least the %s of the ",
      "bootstrap distribution that lies beyond each limit of a %s%% ",
      "interval, so the limits depend on how those replications were ",
      "treated; report the interval with these counts%s"),
      in_question, B, .pct(in_question / B),
      if (counts$improper_kept > 0L) {
        "came from an improper solution or were set aside"
      } else {
        "were set aside"
      }, .pct(tail_area),
      .format_conf_pct(conf_level),
      if (counts$no_solution > 0L) {
        paste0(". With replications that have no single-factor solution ",
               "among them, the lower limit is likely too high: the ",
               "estimator runs high near the boundary where no solution ",
               "exists, and the resamples beyond it, which point to low ",
               "reliability, are not in the distribution")
      } else {
        ""
      }))
  }
  for (j in seq_len(m)) {
    if (!is.null(res[[j]]$fallback)) {
      notes <- c(notes, sprintf(
        "The BCa interval%s could not be formed (%s); the percentile interval is reported instead, and the ci_method attribute says so",
        if (m > 1L) paste0(" for ", stat_names[j]) else "",
        res[[j]]$fallback))
    }
  }
  # A BCa interval withheld at its adjusted positions (the text of the
  # rule, R/dmar_bootstrap.R, without a function name).
  bca_withheld <- !vapply(res, function(r) is.null(r$bca_withheld),
                          logical(1))
  for (j in which(bca_withheld)) {
    text <- sub("^\\(\\): no BCa interval is reported",
                paste0("No BCa interval is reported",
                       if (m > 1L) paste0(" for ", stat_names[j]) else ""),
                .dmar_boot_message(res[[j]]$bca_withheld))
    notes <- c(notes, sub("\\.$", "", text))
  }
  constant <- pick("constant", logical(1))
  if (any(constant)) {
    notes <- c(notes, sprintf(paste0(
      "Every replication that returned a value gave the same estimate%s, ",
      "so the bootstrap distribution is a single point and the interval ",
      "has no width: the coefficient does not vary from resample to ",
      "resample, which happens when every item is an exact linear ",
      "function of the same variable (copies of one item, up to a ",
      "change of scale and origin)"),
      if (m > 1L) paste0(" of ", paste(stat_names[constant],
                                       collapse = ", ")) else ""))
  }
  # The verdict in the words of the rule. B_used is the smallest count any
  # reported limit rests on (the replications that returned a value, or
  # with several statistics the fewest among those with an interval).
  reported <- !target_short & !bca_withheld
  B_used <- if (enough && any(reported)) {
    as.integer(min(target_returned[reported]))
  } else {
    as.integer(returned)
  }
  boot_status <- if (!enough) {
    if (short_of_half || short_of_cases) "insufficient_share" else
      "insufficient_count"
  } else if (any(target_short) || any(bca_withheld)) {
    if (any(target_returned < B / 2)) "insufficient_share" else
      "insufficient_count"
  } else if (!is.na(required) && B < required) {
    "exploratory"
  } else if (counts$failed > 0L) {
    "dropped"
  } else {
    "complete"
  }
  # The number required that the result records. For a BCa interval
  # withheld at its adjusted positions it is the count at those positions
  # (the largest among the statistics withheld; NA when a limit lies beyond
  # every replication), the number the message states, and not the number
  # the request was held to.
  required_recorded <- required
  if (any(bca_withheld)) {
    at_positions <- vapply(res[bca_withheld], function(r)
      as.numeric(r$bca_withheld$B_required), numeric(1))
    required_recorded <- if (anyNA(at_positions)) NA_integer_ else
      as.integer(max(at_positions))
  }
  fields <- list(B = as.integer(B), B_used = B_used,
                 n_dropped = as.integer(B) - B_used,
                 B_required = required_recorded,
                 share_dropped = (B - B_used) / B,
                 caution = isTRUE(caution && enough), status = boot_status)
  if (any(bca_withheld)) {
    first_withheld <- res[[which(bca_withheld)[1L]]]$bca_withheld
    fields$positions <- first_withheld$positions
    fields$B_suggest <- first_withheld$B_suggest
  }
  # One warning per call. Its class says whether an interval was formed:
  # dmar_bootstrap_insufficient when one was withheld by the share or the
  # count rule, dmar_bootstrap_replications_dropped when replications were
  # set aside and every interval was formed. A note about something else
  # (a BCa interval replaced by the percentile interval, replicates that
  # are all the same) has the class dmar_bootstrap_note and the same
  # fields.
  # Every warning names the function that was called (`fun`) and carries
  # the name as a field.
  head <- if (is.null(fun)) "" else paste0(fun, "(): ")
  if (!is.null(fun)) fields$fun <- fun
  if (length(notes)) {
    text <- paste0(head, paste(notes, collapse = ". "), ".")
    if (startsWith(boot_status, "insufficient")) {
      do.call(.dmar_warn, c(list(text, "bootstrap_insufficient"), fields))
    } else if (counts$failed > 0L) {
      do.call(.dmar_warn, c(list(text, "bootstrap_replications_dropped"),
                            fields))
    } else {
      do.call(.dmar_warn, c(list(text, "bootstrap_note"), fields))
    }
  }
  # Under options(dmar.bootstrap_exploratory = TRUE) a run below the
  # requirement is marked in every call.
  if (identical(boot_status, "exploratory")) {
    marked <- .dmar_boot_settle(
      B, B_used, kind = boot_kind, tail = (1 - conf_level) / 2,
      conf_level = conf_level, fun = "", quantity = "the coefficient",
      dropped_because = "returned no value", signal = FALSE)
    text <- paste0(head, sub("^\\(\\): ", "",
                             .dmar_boot_message(marked, "exploratory")))
    do.call(.dmar_warn, c(list(text, "bootstrap_exploratory"), fields))
  }

  counts$improper_rule <- improper
  counts$caution <- caution
  counts$require_half <- TRUE
  counts$required <- required_recorded
  counts$bca_withheld <- any(bca_withheld)
  counts$bootstrap_status <- boot_status
  counts$B_used <- B_used
  counts$target_returned <- stats::setNames(target_returned, stat_names)
  counts$target_withheld <- stats::setNames(vapply(res, function(x)
    anyNA(c(x$lower, x$upper)), logical(1)), stat_names)
  counts$acceleration <- pick("acceleration", numeric(1))
  counts$acceleration_method <- pick("acceleration_method", character(1))
  counts$jackknife <- jack
  counts$ncpus <- if (par_type == "no") 1L else ncpus
  counts$t <- if (m == 1L) {
    as.numeric(tt[, 3L])
  } else {
    est_cols <- tt[, m + 1L + seq_len(m), drop = FALSE]
    dimnames(est_cols) <- list(NULL, stat_names)
    est_cols
  }
  status_code <- as.integer(tt[, code_col])
  status_code[status_code == .BOOT_NO_SOLUTION] <- 2L
  counts$status <- factor(.BOOT_STATUS[status_code], levels = .BOOT_STATUS)
  # The mean of the replicates whose value an interval uses travels back
  # beside the interval, not in `bootstrap`, so the attribute every
  # member of the family carries is unchanged. When no interval is formed
  # it is NA: the replicates that returned a value then describe a
  # minority of the resamples, or a single one (through the first versions
  # of 1.1.0 it was their mean, from which reliability_H() reported a bias
  # beside an interval it did not report).
  replicate_mean <- vapply(seq_len(m), function(j) {
    values <- tt[, j]
    n <- sum(is.finite(values))
    if (!enough || n < 2L || n < B / 2) return(NA_real_)
    mean(values[is.finite(values)])
  }, numeric(1))
  list(se = pick("se", numeric(1)), lower = pick("lower", numeric(1)),
       upper = pick("upper", numeric(1)),
       ci_method = pick("method", character(1)), bootstrap = counts,
       replicate_mean = replicate_mean)
}


# The fate of a replication, in the order of the status codes the
# bootstrap statistic returns (1 to 6); .bootstrap_ci() reports them as
# the levels of the `status` factor. "not_converged" (since 1.1.0, for
# the bifactor indices) is a fit that did not converge, which a point
# function names with a dmar_reason of "not_converged"; a point function
# that does not name it reports such a fit as "non_finite".
.BOOT_STATUS <- c("proper", "improper", "error", "non_finite",
                  "undetermined", "not_converged")

# The status code of a replication set aside because three items had no
# single-factor solution (a dmar_reason of "no_solution";
# .triad_no_solution()). It is an improper replication, reported under
# the status "improper", and is counted also on its own (`no_solution`),
# so the warning and the printed line can name it apart from a Heywood
# case: no loading was estimated for it, and it says something different
# about where the bootstrap distribution lies. The code is 8 because the
# bootstrap then impute loop of categorical omega gives 7 to its own
# seventh status, "imputation" (R/multiple_imputation_internals.R).
.BOOT_NO_SOLUTION <- 8L


# .boot_statistic(point_fn, improper, m = 1L)
#
# The statistic .bootstrap_ci() hands boot::boot(): it applies point_fn
# to the resample and classifies the result (.boot_classify()), for a
# point function that returns m statistics. Built in its own small
# closure so a socket worker receives point_fn, the rule, and m, not the
# whole frame of .bootstrap_ci() (which holds the cluster).
.boot_statistic <- function(point_fn, improper, m = 1L) {
  force(point_fn)
  force(improper)
  force(m)
  function(d, i) .boot_replicate(point_fn, d[i, , drop = FALSE], improper, m)
}


# .boot_replicate(point_fn, d, improper, m = 1L)
#
# One replication: point_fn on the rows `d`, with an error caught and
# classified like any other outcome.
.boot_replicate <- function(point_fn, d, improper, m = 1L) {
  value <- tryCatch(withCallingHandlers(point_fn(d), warning = function(w) {
    message <- conditionMessage(w)
    if (grepl("lavaan", message, fixed = TRUE) &&
        grepl("zero coverage|less than 10%[[:space:]]+coverage", message)) {
      # The parent fit reports this H1 qualification. Repeating it for
      # each valid constrained-H0 refit adds no information about failure.
      invokeRestart("muffleWarning")
    }
  }), error = function(e) e)
  .boot_classify(value, improper, m)
}


# .boot_classify(val, improper, m = 1L)
#
# The 2m + 1 numbers a replication contributes to boot's replicate
# matrix: c(the m values the intervals use, status code, the m
# estimates), which for a single statistic is c(value, status code,
# estimate). See the point function contract in .bootstrap_ci() and
# .BOOT_STATUS for the codes (and .BOOT_NO_SOLUTION for the code of a
# replication with no single-factor solution, which is NA in the value
# and the estimate as an improper solution with no coefficient is). The
# values are NA when the replication is
# set aside; the estimates are NA when there are none (an error, a
# non-finite return, an improper solution that yields no coefficient,
# a fit that did not converge). A replication is used or set aside as a
# whole: one non-finite statistic sets aside all m, so every interval
# rests on the same replications, except a statistic the point function
# marks as not defined for the replication (the dmar_undefined
# attribute; see .bootstrap_ci()), which is NA and is not consulted.
.boot_classify <- function(val, improper, m = 1L) {
  none <- rep(NA_real_, m)
  if (inherits(val, "error")) return(c(none, 3, none))
  reason <- attr(val, "dmar_reason", exact = TRUE)
  undefined <- attr(val, "dmar_undefined", exact = TRUE)
  est <- if (is.numeric(val) && length(val) == m) {
    as.numeric(val)
  } else {
    none
  }
  needed <- rep(TRUE, m)
  if (is.logical(undefined) && length(undefined) == m &&
      !anyNA(undefined) && !undefined[1L]) {
    needed <- !undefined
    est[undefined] <- NA_real_
  }
  if (!all(is.finite(est[needed]))) {
    code <- if (identical(reason, "undetermined_polychoric")) {
      5
    } else if (identical(reason, "improper")) {
      2
    } else if (identical(reason, "no_solution")) {
      .BOOT_NO_SOLUTION
    } else if (identical(reason, "not_converged")) {
      6
    } else {
      4
    }
    return(c(none, code, none))
  }
  if (identical(reason, "improper")) {
    return(c(if (improper == "keep") est else none, 2, est))
  }
  c(est, 1, est)
}


# .boot_tally(tt, attempted)
#
# The counts of a set of replications from their rows of the replicate
# matrix (.boot_classify()): the bootstrap's or the jackknife's, with
# 2m + 1 columns for m statistics. Error messages, and the errors of too
# few distinct cases among them (`too_few_cases`, since 1.1.0), are filled
# in by the caller (.boot_refit_errors()). `improper` counts the improper
# replications set aside, those with no single-factor solution among
# them, and `no_solution` counts those alone (since 1.1.0).
.boot_tally <- function(tt, attempted) {
  m <- (ncol(tt) - 1L) %/% 2L
  code <- as.integer(tt[, m + 1L])
  used <- is.finite(tt[, 1L])
  error <- sum(code == 3L)
  failed <- sum(!used)
  list(attempted = as.integer(attempted),
       returned = as.integer(sum(used)),
       failed = as.integer(failed),
       reasons = c(error = error, non_finite = failed - error),
       error_messages = integer(),
       undetermined = sum(code == 5L),
       improper = sum(code %in% c(2L, .BOOT_NO_SOLUTION) & !used),
       improper_kept = sum(code == 2L & used),
       no_solution = sum(code == .BOOT_NO_SOLUTION),
       not_converged = sum(code == 6L),
       too_few_cases = 0L)
}


# .boot_refit_errors(point_fn, data, index_sets, par_type = "no",
#                    ncpus = 1L, cl = NULL)
#
# The distinct error messages of the replications that raised an error,
# with their counts, recovered by refitting each of them once on its own
# rows (`index_sets`, a list of row index vectors), and how many of those
# errors were the refusal of a resample with too few distinct cases
# (class dmar_too_few_cases; .relia_check_resample()). A condition caught
# inside a worker process does not travel back to the parent as a status
# code does, so the engine recovers every message this way, whatever
# ncpus is, and neither the warning text nor the count depends on it.
# The refits run on the engine's own processes (par_type, ncpus, and cl,
# as .bootstrap_ci() set them up), through .boot_lapply() without a
# progress callback, as the jackknife's do: one call over all of
# index_sets, with the function sent to each socket worker once. A
# replication raises an error most often when its fit runs off along a
# likelihood with no maximum (class dmar_no_maximum), and such a fit is
# among the slowest lavaan makes, so the refits are not cheap. On a
# bootstrap on eight processes (40 replications of a nine-item bifactor
# model, 6 with no maximum), the refits took 32 s of a 50 s run when they
# ran in this process alone, as they did through the first versions of
# 1.1.0; run on those eight processes, they take about 8 s of 25 s.
# Each refit returns only its message and whether it was a refusal for
# too few cases (.boot_refit_one()), so nothing else crosses back from a
# worker. A worker that died returns no such record (mclapply() gives a
# try-error or NULL in its place), and that replication is refit in this
# process, so the messages and the counts are those of the serial path.
# Warnings of the refit are muffled: the replication itself already gave
# them.
# (Through the first versions of 1.1.0, .boot_refit_messages() returned
# the messages alone.)
#
# Returns:
#   A list with `messages`, a named integer vector (empty when none), and
#   `too_few_cases`, an integer.
.boot_refit_errors <- function(point_fn, data, index_sets, par_type = "no",
                               ncpus = 1L, cl = NULL) {
  one <- .boot_refit_one(point_fn, data)
  res <- .boot_lapply(index_sets, one, par_type, ncpus, cl)
  lost <- which(!vapply(res, function(r) {
    is.list(r) && is.logical(r$error) && length(r$error) == 1L &&
      !is.na(r$error)
  }, logical(1)))
  for (k in lost) res[[k]] <- one(index_sets[[k]])
  res <- res[vapply(res, function(r) r$error, logical(1))]
  if (!length(res)) {
    return(list(messages = integer(), too_few_cases = 0L))
  }
  few <- sum(vapply(res, function(r) r$too_few_cases, logical(1)))
  tab <- table(vapply(res, function(r) r$message, character(1)))
  list(messages = stats::setNames(as.integer(tab), names(tab)),
       too_few_cases = as.integer(few))
}


# .boot_refit_one(point_fn, data)
#
# The function of a row index vector that .boot_refit_errors() hands its
# workers: point_fn refit on those rows of `data`, returning
# list(error = FALSE) when the refit raises no error, and otherwise
# list(error = TRUE, message =, too_few_cases =), the condition's message
# and whether it is of class dmar_too_few_cases. Built, like
# .boot_leave_one_out(), in its own small closure with its arguments
# forced, so a socket worker receives the point function and the data as
# values and nothing of the caller's frame (which holds the cluster).
.boot_refit_one <- function(point_fn, data) {
  force(point_fn)
  force(data)
  function(i) {
    tryCatch({
      suppressWarnings(point_fn(data[i, , drop = FALSE]))
      list(error = FALSE)
    }, error = function(e) {
      list(error = TRUE, message = conditionMessage(e),
           too_few_cases = inherits(e, "dmar_too_few_cases"))
    })
  }
}


# .boot_jackknife(point_fn, data, improper, par_type, ncpus, cl, m = 1L,
#                 progress = NULL, B = NULL)
#
# The leave-one-out replications behind the BCa acceleration, as an
# n x (2m + 1) matrix in the form of boot's replicate matrix, spread
# over the bootstrap's workers when there are any. With a progress
# callback (see .bootstrap_ci()) it is called after each refit, or each
# block of refits on several processes, with done and total both B, the
# bootstrap's count, which the jackknife follows; a FALSE stops it.
.boot_jackknife <- function(point_fn, data, improper, par_type, ncpus, cl,
                            m = 1L, progress = NULL, B = NULL) {
  n <- nrow(data)
  one <- .boot_leave_one_out(point_fn, data, improper, m)
  res <- .boot_lapply(seq_len(n), one, par_type, ncpus, cl,
                      progress = progress, done = function(k) B, total = B)
  .boot_rows(res, m)
}


# .boot_rows(res, m)
#
# A list of replications, each the 2m + 1 numbers .boot_classify()
# returns, as the rows of a matrix in the form of boot's replicate
# matrix. A worker that died returns no such vector (mclapply() gives a
# try-error or NULL in its place), and it is counted as an error.
.boot_rows <- function(res, m) {
  width <- 2L * m + 1L
  res <- lapply(res, function(v) {
    if (is.numeric(v) && length(v) == width) {
      v
    } else {
      c(rep(NA_real_, m), 3, rep(NA_real_, m))
    }
  })
  matrix(unlist(res), ncol = width, byrow = TRUE)
}


# .boot_with_progress(data, statistic, B, m, par_type, ncpus, cl,
#                     progress)
#
# The boot object .bootstrap_ci() works from, when a progress callback was
# given (see .bootstrap_ci()). boot::boot() is called as the engine calls
# it without one, with the statistic replaced by one that computes
# nothing, so it draws the same B resamples from the same random number
# stream and records the state it drew them from; boot::boot.array()
# then returns them, and the replications run here, in order, through
# .boot_lapply(), which calls the callback after each one (or each
# block). The estimate on the data, the replicate matrix, and the
# statistic are then put in the object where boot::boot() puts them, so
# boot::boot.ci(), boot::empinf(), and boot::boot.array() read it as
# they read boot's own. boot computes the statistic on the data after
# drawing the resamples and before the first replication, and it is
# computed here in the same place.
.boot_with_progress <- function(data, statistic, B, m, par_type, ncpus, cl,
                                progress) {
  width <- 2L * m + 1L
  boot_out <- boot::boot(data = data,
                         statistic = function(d, i) rep(NA_real_, width),
                         R = B, stype = "i", parallel = "no")
  idx <- boot::boot.array(boot_out, indices = TRUE)
  rows <- lapply(seq_len(B), function(r) idx[r, ])
  rm(idx)
  boot_out$t0 <- statistic(data, seq_len(nrow(data)))
  res <- .boot_lapply(rows, .boot_on_rows(statistic, data), par_type,
                      ncpus, cl, progress = progress)
  boot_out$t <- .boot_rows(res, m)
  boot_out$statistic <- statistic
  boot_out
}


# .boot_on_rows(statistic, data)
#
# The function of a resample's row indices that .boot_with_progress()
# hands .boot_lapply(): the statistic on those rows of `data`. Built in
# its own small closure with its arguments forced, like .boot_statistic(),
# so a socket worker receives the statistic and the data and nothing of
# the caller's frame (which holds the cluster and every resample).
.boot_on_rows <- function(statistic, data) {
  force(statistic)
  force(data)
  function(i) statistic(data, i)
}


# .boot_lapply(X, fun, par_type, ncpus, cl, progress = NULL,
#              done = function(k) k, total = length(X))
#
# lapply(X, fun) spread over the engine's processes as boot::boot()
# spreads its replications: parallel::mclapply() on forked workers,
# parallel::parLapply() on the socket cluster `cl`, lapply() in this
# process. With progress = NULL that is one call over all of X, as the
# engine has always made it. With a callback the elements run in blocks
# (.boot_block_size(): one at a time in this process), and after each
# block the callback is called (.boot_progress()) with done(k), k the
# elements of X completed so far, and `total`; a FALSE stops the run
# before the next block. Every element is computed by the same function
# on the same input either way, so the results do not depend on the
# blocks.
#
# On a socket cluster each parLapply() call serializes its function, with
# everything its environment holds, to every worker: the data, and for the
# bootstrap then impute loop the B x n resamples and the B seeds as well
# (.mi_replicate_fn()); a forked worker inherits them instead. So that the
# blocks, about fifty calls, do not send them again each time, `fun` is
# sent to each worker once, before the first block
# (.boot_worker_stage()), and a block sends only its own elements of X
# with a function that applies the one staged (.boot_worker_run()). The
# jackknife, which runs on the same cluster after the replications, stages
# its own function in place of theirs, and what was staged goes with the
# cluster, which belongs to the call of the engine and is stopped when it
# returns. Without a callback the one parLapply() call sends `fun` to each
# worker once anyway, and nothing is staged.
.boot_lapply <- function(X, fun, par_type, ncpus, cl, progress = NULL,
                         done = function(k) k, total = length(X)) {
  if (is.null(progress)) {
    return(switch(par_type,
                  multicore = parallel::mclapply(X, fun, mc.cores = ncpus),
                  snow = parallel::parLapply(cl, X, fun),
                  lapply(X, fun)))
  }
  if (par_type == "snow") {
    parallel::clusterCall(cl, .boot_worker_stage, fun)
  }
  run <- function(x) {
    switch(par_type,
           multicore = parallel::mclapply(x, fun, mc.cores = ncpus),
           snow = parallel::parLapply(cl, x, .boot_worker_run),
           lapply(x, fun))
  }
  n <- length(X)
  size <- .boot_block_size(n, if (par_type == "no") 1L else ncpus)
  out <- vector("list", n)
  for (first in seq.int(1L, n, by = size)) {
    k <- first:min(n, first + size - 1L)
    out[k] <- run(X[k])
    .boot_progress(progress, done(max(k)), total)
  }
  out
}


# .boot_worker, .boot_worker_stage(fun), .boot_worker_run(x)
#
# The function a socket worker applies to the elements of each block when
# .boot_lapply() runs in blocks. .boot_worker_stage() is sent to each
# worker once, with `fun`, by parallel::clusterCall(), and keeps `fun` in
# .boot_worker, an environment of the worker's own DMAR namespace made for
# it (a socket worker loads the installed DMAR to run any of the engine's
# functions); .boot_worker_run(), which holds nothing but a reference to
# that namespace, is what each block sends, and applies the function kept
# there to one element. In this process .boot_worker stays empty.
.boot_worker <- new.env(parent = emptyenv())

.boot_worker_stage <- function(fun) {
  assign("fun", fun, envir = .boot_worker)
  invisible(NULL)
}

.boot_worker_run <- function(x) {
  get("fun", envir = .boot_worker, inherits = FALSE)(x)
}


# .boot_block_size(n, ncpus)
#
# How many replications run between two calls of a progress callback: one
# when they run in this process, and on ncpus processes about a fiftieth
# of the n, rounded up to a multiple of ncpus so that every process has
# the same number of them (and at least one), which bounds the cost of
# starting the workers for each block at about fifty starts in all.
.boot_block_size <- function(n, ncpus) {
  if (ncpus <= 1L) return(1L)
  as.integer(ncpus * ceiling(n / (50 * ncpus)))
}


# .boot_progress(progress, done, total)
#
# Calls a progress callback (see .bootstrap_ci()) with the replications
# attempted so far and the number requested, both as integers. It runs
# outside the random number stream (withr::with_preserve_seed()): a
# callback that draws random numbers leaves the stream as it found it.
# An error inside it stops the call with an error of class
# dmar_progress_error, whose message names the callback and which carries
# done, total, and the callback's own condition as `parent`, its class
# intact, so an application can tell a failing callback from a failing
# analysis, and still recognize a condition of its own (a Shiny req(), for
# example) behind it. A return value of FALSE, and only that
# (identical(answer, FALSE)), stops the bootstrap with an error of class
# dmar_bootstrap_cancelled carrying done and total; any other value is
# ignored. Does nothing when progress is NULL.
.boot_progress <- function(progress, done, total) {
  if (is.null(progress)) return(invisible(NULL))
  done <- as.integer(done)
  total <- as.integer(total)
  answer <- withr::with_preserve_seed(tryCatch(
    progress(done, total),
    error = function(e) {
      stop(errorCondition(
        sprintf(paste0("The progress callback failed after %d of %d ",
                       "bootstrap replications: %s"),
                done, total, conditionMessage(e)),
        class = "dmar_progress_error", done = done, total = total,
        parent = e, call = NULL))
    }))
  if (identical(answer, FALSE)) {
    stop(errorCondition(
      sprintf(paste0("The bootstrap was cancelled after %d of %d ",
                     "replications: the progress callback returned FALSE."),
              done, total),
      class = "dmar_bootstrap_cancelled", done = done, total = total,
      call = NULL))
  }
  invisible(NULL)
}


# .check_progress(progress)
#
# Validates the `progress` argument of the bootstrap functions (NULL or
# a function; see .bootstrap_ci()) and returns it.
.check_progress <- function(progress) {
  if (!is.null(progress) && !is.function(progress)) {
    .dmar_refuse_argument(c(
      "'progress' must be NULL or a function of two arguments, the ",
      "bootstrap replications attempted so far and the number ",
      "requested (done and total)."), arg = "progress")
  }
  progress
}


# .boot_constant(t)
#
# Whether the finite replicates in `t` (one column of boot's replicate
# matrix, NA where a replication was set aside) are all the same value,
# to a tolerance of 1e-8 relative to their mean (absolute below a mean
# of one). Every column boot::boot.ci() itself calls constant meets this
# test, whose tolerance is never smaller than boot's.
.boot_constant <- function(t) {
  t <- t[is.finite(t)]
  if (!length(t)) return(FALSE)
  mid <- mean(t)
  all(abs(t - mid) < 1e-8 * max(1, abs(mid)))
}


# .boot_percentile(boot_out, conf_level, j)
#
# The percentile limits of statistic j, c(lower, upper), the
# (1 - C)/2 and (1 + C)/2 quantiles of its returned replicates, from
# boot::boot.ci(type = "perc"). boot.ci() cannot form them when every
# returned replicate is the same value: it prints "All values of t are
# equal" and returns NULL, and when some replications were set aside (NA
# in the column) its test for that case meets if (NA) and stops. Through
# DMAR 1.0.1 either way stopped the call with an error about vapply()
# that named neither. The quantiles of a single point are that point, so
# the limits are then the common value (the bootstrap_se interval of
# such replicates, with a standard deviation of zero, is the same
# point), and .bootstrap_ci() says so in its warning.
.boot_percentile <- function(boot_out, conf_level, j) {
  t <- boot_out$t[, j]
  if (.boot_constant(t)) {
    v <- t[is.finite(t)][1L]
    return(c(v, v))
  }
  # boot.ci() names the component "percent" even for type = "perc";
  # use the full name so the access is exact, not a partial match. The
  # index is required: with more than one column boot.ci() would
  # otherwise read the second as the variance of the first.
  bci <- .boot_ci_endpoints_quiet(
    boot::boot.ci(boot_out, conf = conf_level, type = "perc",
                  index = j)$percent,
    sum(is.finite(t)) <
      .dmar_boot_required("quantile", (1 - conf_level) / 2))
  c(bci[4], bci[5])
}

# boot::boot.ci() warns "extreme order statistics used as endpoints" when
# a limit falls at the smallest or the largest replicate. In an exploratory
# run below the required number (options(dmar.bootstrap_exploratory =
# TRUE); `below` is TRUE) DMAR's own warning of class
# dmar_bootstrap_exploratory already says how few replications each limit
# rests on, so that unclassed warning is not repeated beside it. In every
# other run it passes.
.boot_ci_endpoints_quiet <- function(expr, below) {
  if (!(isTRUE(below) && .dmar_boot_exploratory())) return(expr)
  withCallingHandlers(expr, warning = function(w) {
    if (grepl("extreme order statistics used as endpoints",
              conditionMessage(w), fixed = TRUE)) {
      invokeRestart("muffleWarning")
    }
  })
}


# .boot_leave_one_out(point_fn, data, improper, m = 1L)
#
# The function of k that .boot_jackknife() hands its workers: the
# replication without case k. Built, like .boot_statistic(), in its own
# small closure with its arguments forced, so a socket worker receives
# the point function, the data, the rule, and m as values; an unforced
# argument would travel as a promise whose environment the worker
# cannot evaluate (every leave-one-out refit then failed on the socket
# path).
.boot_leave_one_out <- function(point_fn, data, improper, m = 1L) {
  force(point_fn)
  force(data)
  force(improper)
  force(m)
  function(k) .boot_replicate(point_fn, data[-k, , drop = FALSE], improper, m)
}


# .check_ncpus(ncpus)
#
# Validates the `ncpus` argument of the reliability functions and returns
# it as an integer.
.check_ncpus <- function(ncpus) {
  if (!is.numeric(ncpus) || length(ncpus) != 1L || !is.finite(ncpus) ||
      ncpus < 1 || ncpus %% 1 != 0) {
    .dmar_refuse_argument(c(
      "'ncpus' must be a single whole number of at least 1, the ",
      "number of processes the bootstrap replications are spread ",
      "over."), arg = "ncpus")
  }
  as.integer(ncpus)
}


# .boot_parallel_type(ncpus)
#
# boot's `parallel` setting for a number of processes: "no" for one,
# "multicore" (forked workers) where the platform has fork, and "snow"
# (a socket cluster) on Windows, which has not.
.boot_parallel_type <- function(ncpus) {
  if (ncpus <= 1L) {
    "no"
  } else if (.Platform$OS.type == "windows") {
    "snow"
  } else {
    "multicore"
  }
}


# .pct(p)
#
# A proportion as a percent for a message, to three significant figures
# ("2.5%", "7.03%"). Display only.
.pct <- function(p) {
  paste0(format(100 * p, digits = 3, trim = TRUE), "%")
}


# .bootstrap_reason_phrase(reasons, error_messages, undetermined = 0L,
#                          improper = 0L, not_converged = 0L,
#                          no_solution = 0L, refit_failed = 0L)
#
# The parenthetical that names why replications were set aside, shared by
# the warning in .bootstrap_ci() and the one-line count that
# .bootstrap_footer_line() adds under a table carrying a `bootstrap`
# attribute: "3 raised an error (most often: ...), 2 returned a
# non-finite value". The most frequent error message is quoted, first
# line only, so the user sees the cause without unpacking the attribute.
# The non-finite replications whose cause is named (`undetermined`,
# `improper`, and `not_converged`, all subsets of
# reasons[["non_finite"]]; `improper` counts the improper replications
# set aside, not those kept) are named by that cause, and only the rest
# are described as non-finite. Among the improper ones, those of three
# items with no single-factor solution (`no_solution`, a subset of
# `improper`) are named apart from the Heywood cases, since no loading
# was estimated for them. A `bootstrap` attribute written before
# these counts existed has none of them, and the defaults of zero
# describe it as before. The bootstrap then impute interval of
# categorical omega under missing = "mi" (.omega_c_mi() in
# R/multiple_imputation_internals.R) adds a further reason,
# `imputation`, for the replications whose imputation itself failed; it
# is not part of reasons[["non_finite"]], and every other result's
# `reasons` has no such element, so the phrase is unchanged for them.
#
# `refit_failed`, a further subset of reasons[["non_finite"]], counts the
# replications whose refit failed or did not converge where the
# bootstrap cannot tell the two apart: lavaan's bootstrap
# (lavaan::lavBootstrap(), behind average_variance_extracted()) returns
# NA for either without saying which. They are named for both causes,
# "failed to refit or did not converge", the words of the warning that
# function gives, where "returned a non-finite value" named neither;
# every other bootstrap has none, and its phrase is unchanged.
.bootstrap_reason_phrase <- function(reasons, error_messages = integer(),
                                     undetermined = 0L, improper = 0L,
                                     not_converged = 0L, no_solution = 0L,
                                     refit_failed = 0L) {
  if (is.null(undetermined)) undetermined <- 0L
  if (is.null(improper)) improper <- 0L
  if (is.null(not_converged)) not_converged <- 0L
  if (is.null(no_solution)) no_solution <- 0L
  if (is.null(refit_failed)) refit_failed <- 0L
  parts <- character()
  if (reasons[["error"]] > 0L) {
    msg <- ""
    if (length(error_messages)) {
      top <- names(error_messages)[which.max(error_messages)]
      top <- trimws(gsub("\\s+", " ", strsplit(top, "\n")[[1]][1]))
      # The message's own closing period would end the parenthesis
      # ("refused.)").
      top <- sub("[.]+$", "", top)
      msg <- sprintf(" (most often: %s)", top)
    }
    parts <- c(parts, sprintf("%d raised an error%s", reasons[["error"]], msg))
  }
  if (undetermined > 0L) {
    parts <- c(parts, sprintf(
      "%d left a polychoric correlation undetermined",
      undetermined))
  }
  if (improper - no_solution > 0L) {
    parts <- c(parts, sprintf(
      "%d gave an improper (Heywood) solution",
      improper - no_solution))
  }
  if (no_solution > 0L) {
    parts <- c(parts, sprintf("%d had no single-factor solution",
                              no_solution))
  }
  if (not_converged > 0L) {
    parts <- c(parts, sprintf("%d did not converge", not_converged))
  }
  if (refit_failed > 0L) {
    parts <- c(parts, sprintf("%d failed to refit or did not converge",
                              refit_failed))
  }
  other <- reasons[["non_finite"]] - undetermined - improper -
    not_converged - refit_failed
  if (other > 0L) {
    parts <- c(parts, sprintf("%d returned a non-finite value", other))
  }
  if ("imputation" %in% names(reasons) && reasons[["imputation"]] > 0L) {
    parts <- c(parts, sprintf("%d could not be imputed",
                              reasons[["imputation"]]))
  }
  paste(parts, collapse = ", ")
}


# .bootstrap_footer_line(x)
#
# The line a table carrying a `bootstrap` attribute prints under itself,
# which .dmar_tbl_footer() (R/dmar_output_helpers.R) takes from here so
# the console and a knitted report show the same text: how many
# replications were attempted and returned a value, why the rest were set
# aside, how many of those used came from an improper solution, the
# jackknife behind a BCa acceleration, and the caution when the
# replications in question reach the tail area beyond a limit, with the
# likely direction of the error when some had no single-factor solution
# (see .bootstrap_ci()), or, when no interval was formed (fewer than two
# replications returned a value, more than half had too few distinct
# cases, or fewer than half returned a value), that no
# interval is reported, and no caution about limits that do not
# exist; the caution names an improper solution only when replications
# from one were kept; for a table of several factors that counts its
# improper replications by factor (improper_kept_by_factor), the caution
# gives the share of each factor it holds for. An attribute written
# before the 1.1.0 elements existed prints as it did.
#
# Inputs:
#   x : the table (its `bootstrap` attribute, and its conf_level, read
#       from a `conf_level` attribute or row, for the caution).
# Returns:
#   A single character string, or NULL when x carries no bootstrap.
.bootstrap_footer_line <- function(x) {
  bs <- attr(x, "bootstrap", exact = TRUE)
  if (!is.list(bs) || is.null(bs$attempted)) return(NULL)
  line <- sprintf("Bootstrap: %d replications attempted, %d returned a value",
                  bs$attempted, bs$returned)
  if (isTRUE(bs$failed > 0L)) {
    line <- paste0(line, sprintf(
      " (%d set aside: %s)", bs$failed,
      .bootstrap_reason_phrase(bs$reasons, bs$error_messages,
                               bs$undetermined, bs$improper,
                               bs$not_converged, bs$no_solution,
                               bs[["refit_failed"]])))
  }
  if (isTRUE(bs$improper_kept > 0L)) {
    line <- paste0(line, sprintf(
      "; %d of the %d came from an improper (Heywood) solution and %s kept",
      bs$improper_kept, bs$returned,
      if (bs$improper_kept == 1L) "was" else "were"))
  }
  jk <- bs$jackknife
  if (is.list(jk) && isTRUE(jk$attempted > 0L)) {
    line <- paste0(line, sprintf(
      "; jackknife for the BCa acceleration: %d attempted, %d returned a value",
      jk$attempted, jk$returned))
    if (isTRUE(jk$failed > 0L)) {
      line <- paste0(line, sprintf(
        " (%d set aside: %s)", jk$failed,
        .bootstrap_reason_phrase(jk$reasons, jk$error_messages,
                                 jk$undetermined, jk$improper,
                                 jk$not_converged, jk$no_solution)))
    }
  }
  if (isTRUE(bs$returned >= 2L) &&
      isTRUE(bs$too_few_cases > bs$attempted / 2)) {
    # No interval: more than half had too few distinct cases
    # (.bootstrap_ci()).
    return(paste0(line, "; more than half had too few distinct cases, so ",
                  "no interval is reported"))
  }
  if (isTRUE(bs$interval_withheld)) {
    return(paste0(line, "; at least ", bs$minimum_returned,
                  " returned estimates are required, so no interval is reported"))
  }
  if (isTRUE(bs$require_half) && isTRUE(bs$returned >= 2L) &&
      isTRUE(bs$returned < bs$attempted / 2)) {
    # No interval: fewer than half returned a value (.bootstrap_ci()).
    return(paste0(line, "; fewer than half returned a value, so no ",
                  "interval is reported"))
  }
  if (isTRUE(bs$bca_withheld) && length(bs$target_withheld) <= 1L) {
    # No BCa interval: too few replications beyond a limit where the bias
    # correction and the acceleration place it (the position rule of
    # ?dmar_bootstrap).
    return(paste0(line, if (is.na(bs$required)) {
      paste0("; the bias correction and the acceleration place a limit ",
             "beyond every replication, so no BCa interval is reported")
    } else {
      paste0("; ", bs$required, " are required where the bias correction ",
             "and the acceleration place the limits, so no BCa interval ",
             "is reported")
    }))
  }
  if (identical(bs$bootstrap_status, "insufficient_count") &&
      isTRUE(bs$returned >= 2L) && isTRUE(bs$attempted < bs$required)) {
    # No interval: fewer were requested than the level requires (the
    # count rule of ?dmar_bootstrap, which is held to the request).
    return(paste0(line, "; fewer than the ", bs$required, " required at ",
                  "this level were requested, so no interval is reported"))
  }
  if (isTRUE(bs$returned < 2L) && isTRUE(bs$attempted >= 1L)) {
    # No interval either, and no caution about its limits (.bootstrap_ci()).
    return(paste0(line, "; fewer than two returned a value, so no ",
                  "interval is reported"))
  }
  if (length(bs$target_withheld) > 1L && any(bs$target_withheld)) {
    line <- paste0(line, "; intervals withheld for ",
                   paste(names(bs$target_withheld)[bs$target_withheld], collapse = ", "),
                   "; defined target counts: ",
                   paste(paste0(names(bs$target_returned), " = ", bs$target_returned),
                         collapse = ", "))
  }
  if (isTRUE(bs$caution)) {
    in_question <- bs$failed + (bs$improper_kept %||% 0L)
    conf_level <- attr(x, "conf_level", exact = TRUE)
    if (is.null(conf_level) && "term" %in% names(x) &&
        "conf_level" %in% x$term) {
      conf_level <- x[[2L]][x$term == "conf_level"][1L]
    }
    # A table of several factors with intervals (average_variance_
    # extracted()) counts its improper replications for each factor, and
    # the caution gives each factor its own share, as its warning does:
    # improper_kept counts a replication improper for any factor, so its
    # share was not that of any one factor.
    by_factor <- bs$improper_kept_by_factor
    if (!is.null(bs$factors_with_intervals)) {
      by_factor <- by_factor[intersect(names(by_factor), bs$factors_with_intervals)]
      if (length(by_factor) == 1L) in_question <- bs$failed + unname(by_factor)
    }
    whose <- if (isTRUE(if (length(by_factor)) any(by_factor > 0L) else
      bs$improper_kept > 0L)) "an improper solution"
    if (length(by_factor) >= 2L) {
      n <- bs$failed + by_factor
      hit <- if (is.numeric(conf_level) && length(conf_level) == 1L) {
        vapply(n, .boot_caution, logical(1L), B = bs$attempted,
               conf_level = conf_level)
      } else {
        n > 0L
      }
      if (any(hit)) {
        share <- sprintf("for %s, %s", names(n)[hit],
                         vapply(n[hit] / bs$attempted, .pct, character(1L)))
        share[1L] <- paste(share[1L], "of the replications")
        in_question <- NULL
        whose <- "a solution improper for the factor"
      }
    }
    line <- paste0(line, sprintf(
      "; caution: %s %s%s%s",
      if (is.null(in_question)) {
        paste0(paste(share, collapse = ", and "),
               if (length(share) > 1L) "," else "")
      } else {
        paste(.pct(in_question / bs$attempted), "of the replications")
      },
      if (is.null(whose)) {
        "were set aside"
      } else {
        paste("came from", whose, "or were set aside")
      },
      if (is.numeric(conf_level) && length(conf_level) == 1L) {
        sprintf(", at least the %s beyond each limit",
                .pct((1 - conf_level) / 2))
      } else {
        ""
      },
      if (isTRUE(bs$no_solution > 0L)) {
        "; the lower limit is likely too high"
      } else {
        ""
      }))
  }
  if (identical(bs$bootstrap_status, "exploratory")) {
    line <- paste0(line, "; at least ", bs$required, " are required for a ",
                   "reported interval. Exploratory: below the number ",
                   "required for a reported interval")
  }
  line
}


# -----------------------------------------------------------------------------
# (F) Tidy result formatter
# -----------------------------------------------------------------------------

# .relia_result(estimate, se, lower, upper, conf_level, N, J,
#               coefficient, ci_method, B = NA_integer_,
#               N_complete = N, missing = "listwise", aux = NULL,
#               se_transformed = NA_real_, se_transform_scale = NULL,
#               bootstrap = NULL)
#
# Pack the point estimate, standard error, CI bounds, and bookkeeping
# values into the package-standard data.frame(term, value).
# Reliability is bounded on [0, 1] by definition, so this helper clamps
# each confidence limit into [0, 1] before returning (see the comment in
# the body), and warns when the estimate itself is negative
# (.relia_warn_negative()). Metadata that does not belong in the
# numeric `value` column (which coefficient was computed, which CI
# method, the missing-data treatment, the auxiliary variables, and the
# number of bootstrap reps when applicable) travels as attributes. A
# bootstrap call also passes `bootstrap`, the replication counts from
# .bootstrap_ci(), which become the `bootstrap` attribute (absent for
# every other method, so no existing attribute changes), with the
# number of replications the interval rests on as the `B_used`
# attribute, the effective count every bootstrap in the package reports
# under that name (since 1.1.0).
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
                          se_transform_scale = NULL,
                          bootstrap = NULL) {
  # Each limit is clamped into [0, 1], the range of a reliability
  # coefficient. The clamp is monotone, so the limits stay in order, and
  # any value in [0, 1] that the unclamped interval covers, the clamped
  # one covers too: clamping never lowers the coverage of a parameter
  # that lies in its range. Through DMAR 1.0.1 the lower limit was only
  # raised to 0 and the upper only lowered to 1, so the interval around a
  # negative alpha or KR-20 that lay wholly below zero came back
  # inverted, a lower limit of 0 above a negative upper limit; it is now
  # [0, 0]. NA limits were considered and rejected: when only the lower
  # unclamped limit is negative, the clamped interval still says how
  # high the coefficient can plausibly be, and the warning below says
  # what a negative estimate means.
  raw_upper <- upper
  if (!is.na(lower)) lower <- min(max(lower, 0), 1)
  if (!is.na(upper)) upper <- min(max(upper, 0), 1)
  if (!is.na(estimate) && estimate < 0) {
    .relia_warn_negative(estimate, coefficient, ci_method, lower, upper,
                         raw_upper)
  }
  terms  <- c("estimate", "se", "ci_lower", "ci_upper",
              "conf_level", "N", "N_complete", "n_items")
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
  # A count, stored as an integer like the other counts of ?dmar_bootstrap.
  if (!is.na(B)) attr(out, "B") <- as.integer(B)
  if (!is.null(bootstrap)) {
    attr(out, "bootstrap") <- bootstrap
    if (!is.null(bootstrap$returned)) attr(out, "B_used") <- bootstrap$returned
    # The counts in the names every bootstrap result of the package
    # carries (?dmar_bootstrap); the list above is the fuller record.
    if (!is.null(bootstrap$bootstrap_status)) {
      if (!is.null(bootstrap$B_used)) attr(out, "B_used") <- bootstrap$B_used
      attr(out, "n_dropped") <- as.integer(bootstrap$attempted) -
        as.integer(attr(out, "B_used"))
      attr(out, "B_required") <- bootstrap$required
      attr(out, "bootstrap_status") <- bootstrap$bootstrap_status
    }
  }
  # Tag with a leading class so the generics::tidy / generics::glance
  # S3 methods dispatch; data.frame behavior is preserved through
  # inheritance.
  class(out) <- c("dmar_reliability", "dmar_tbl", "data.frame")
  out
}


# .relia_warn_negative(estimate, coefficient, ci_method, lower, upper,
#                      raw_upper)
#
# The warning, class dmar_negative_reliability, that a reliability
# estimate is negative. Population alpha (and KR-20, its dichotomous
# case) is at most one but has no lower bound: it is negative whenever
# the item covariances are negative on balance, which in practice
# almost always means an item scored in the opposite direction to the
# others. A reliability coefficient cannot be negative, so the estimate
# lies outside the parameter space, and the warning says so, says what
# usually causes it, and says what became of the limits (clamped into
# [0, 1] by .relia_result(); NA for an interval built on the logit
# scale, which a negative estimate does not have). The class lets a
# caller that loops over many data sets (a Monte Carlo study) muffle it
# by name. The estimate is formatted for the message only.
.relia_warn_negative <- function(estimate, coefficient, ci_method, lower,
                                 upper, raw_upper) {
  label <- switch(coefficient %||% "",
                  alpha = "coefficient alpha",
                  kr20 = "KR-20",
                  omega = "coefficient omega",
                  omega_categorical = "categorical omega",
                  H = "coefficient H",
                  "the reliability coefficient")
  limits <- if (is.na(lower) && is.na(upper)) {
    if (grepl("_logistic", ci_method %||% "", fixed = TRUE)) {
      paste0(" An interval on the logit scale is not defined for a ",
             "negative estimate, so the limits are NA.")
    } else {
      ""
    }
  } else if (!is.na(raw_upper) && raw_upper < 0) {
    paste0(" Each confidence limit is clamped into [0, 1], and the ",
           "interval, which lies wholly below zero, is reported as ",
           "[0, 0].")
  } else {
    " Each confidence limit is clamped into [0, 1]."
  }
  .dmar_warn(
    paste0(
      "The estimate of ", label, " is negative (",
      trimws(formatC(estimate, format = "fg", digits = 3)), "), outside ",
      "the range [0, 1] of a reliability coefficient. A negative ",
      "estimate means that the items covary negatively on balance, most ",
      "often because an item is keyed in the opposite direction to the ",
      "others and needs reverse scoring; the item correlations show ",
      "which.", limits),
    "negative_reliability")
}


# -----------------------------------------------------------------------------
# (G) Input-validation helpers
# -----------------------------------------------------------------------------

# .relia_declared_missing(data)
#
# SPSS user-defined missing values made NA. haven::read_sav(user_na =
# TRUE) keeps them as values in a column of class haven_labelled_spss,
# with the values declared missing in its `na_values` attribute (and a
# range in `na_range`); haven's is.na() method reports them as missing,
# and SPSS leaves them out of every analysis, but complete.cases(),
# lavaan, and mice read them as values. Through version 1.0.1 the
# reliability family analyzed them as data: a 5-point item coded 9 for
# "refused" had 9 as a sixth category of categorical omega under every
# missing data treatment, and as a score of 9 in alpha and omega. The
# column is replaced by its plain values with the declared ones set to
# NA, as read_sav() returns them by default (user_na = FALSE), and a
# message names the columns and counts. haven need not be loaded: the
# attributes are read directly. Other columns, and a matrix, are
# returned unchanged. The columns are walked by position, so a data
# frame with two columns of one name (which listwise deletion accepts)
# has the declared values of each made NA, and the message names a
# repeated name with its position. The raters of cohen_kappa() arrive
# as a list of two vectors and are read the same way; icc() and
# fleiss_kappa(), which refuse missing values, then refuse these.
#
# Inputs:
#   data : the data as supplied (a data frame, a list of vectors, or a
#          matrix).
# Returns:
#   `data`, with each haven_labelled_spss column that declares missing
#   values made a plain numeric vector with those values NA.
.relia_declared_missing <- function(data) {
  if (!is.list(data)) return(data)
  hit <- character()
  nms <- names(data) %||% character(length(data))
  for (j in seq_along(data)) {
    x <- data[[j]]
    if (!inherits(x, "haven_labelled_spss")) next
    v <- as.vector(unclass(x))
    na_values <- attr(x, "na_values", exact = TRUE)
    na_range <- attr(x, "na_range", exact = TRUE)
    declared <- !is.na(v) & v %in% na_values
    if (length(na_range) == 2L) {
      declared <- declared | (!is.na(v) & v >= na_range[1L] &
                                v <= na_range[2L])
    }
    if (!any(declared)) next
    v[declared] <- NA
    data[[j]] <- v
    label <- if (!nzchar(nms[j])) {
      sprintf("column %d", j)
    } else if (sum(nms == nms[j]) > 1L) {
      sprintf("%s (column %d)", nms[j], j)
    } else {
      nms[j]
    }
    hit <- c(hit, sprintf("%d in %s", sum(declared), label))
  }
  if (length(hit)) {
    message("SPSS user-defined missing values (haven_labelled_spss) are ",
            "treated as missing, as SPSS and haven's is.na() treat them: ",
            paste(hit, collapse = ", "), ".")
  }
  data
}


# The interval methods of reliability_omega() and reliability_alpha() that
# full information maximum likelihood and two-stage estimation refuse:
# the closed forms, whose standard errors treat every case as complete;
# ADF, which needs complete raw-data fourth moments; and the profile
# likelihood, not implemented for a FIML fit. Each function refuses them
# under an explicit request with its own message; the default of
# `missing` reads them here (.relia_resolve_missing()).
.RELIA_COMPLETE_DATA_METHODS <- c("feldt", "fisher", "bonett",
                                  "hakstian_whalen", "adf",
                                  "adf_logistic", "likelihood")

# .relia_resolve_missing(missing, supplied, data, aux,
#                        unavailable = character(0), fiml_options = NULL)
#
# The missing data treatment of reliability_omega(), reliability_alpha(),
# and reliability_H() with `data` or `S`; reliability() forwards an
# unsupplied `missing` unsupplied, so its calls resolve here too. The
# default of `missing` is "fiml" since DMAR 1.1.0 (DMAR 1.0.x deleted
# incomplete rows listwise), so that the omega of these functions and
# fa_confirmatory()'s agree on incomplete data. It resolves by the rule
# fa_confirmatory() applies to the paths FIML cannot serve
# (.cfa_k_resolve_missing()), applied here to every path:
#
#   - A covariance matrix, and data in which no value is missing, are
#     analyzed as they stand ("listwise", which then deletes nothing).
#     fa_confirmatory() fits complete data by FIML; these functions do
#     not, because their complete-data intervals (Feldt, Fisher, Bonett,
#     Hakstian-Whalen, ADF, and the profile likelihood, refused under
#     FIML) need the complete data analysis, which also keeps every
#     complete data result what DMAR 1.0.x computed. With complete data
#     the FIML estimate is the listwise one (the saturated FIML moments
#     are the sample moments), so the estimate is FIML's, and everything
#     else is what DMAR 1.0.x computed: standard errors from the expected
#     information and the robust.sem sandwich rather than the observed
#     information and the Huber-White sandwich FIML uses (they differ by
#     about one percent or less on 200 simulated cases), Bonett's
#     interval as the analytic alpha's default, and a bootstrap that
#     recomputes the coefficient from each resample's covariance matrix
#     rather than refitting a saturated FIML model to it (for alpha on 200
#     cases, about 70 times as long). A result DMAR 1.0.x computed
#     therefore changes only where a call that left `missing` at its
#     default analyzes incomplete items by FIML.
#   - A row with no observed item does not make the data incomplete:
#     without auxiliary variables listwise deletion and FIML both drop it,
#     so data complete but for such rows (a blank line at the end of a
#     spreadsheet) are analyzed as they stand, as above.
#   - Incomplete data are analyzed by FIML when every option of the call
#     is available under it. When one is not (`unavailable`, the options
#     the caller named that FIML cannot serve), the call stops and names
#     the choices rather than deleting cases unasked: missing =
#     "listwise" for the complete cases, or `fiml_options` under FIML.
#     The callers check first that those options suit the rest of the
#     call (omega's denominator, alpha's estimator), so that neither
#     choice the message names meets another refusal.
#
# An explicit request is returned as given, and the checks that follow
# refuse an impossible one. Auxiliary variables imply FIML whether or not
# a value is missing, as they have since they were introduced. SPSS
# user-defined missing values count as missing, as they do in the
# analysis (.relia_declared_missing(), whose message is left to
# .relia_resolve_inputs() so that it is given once).
#
# Inputs:
#   missing      : the matched value of `missing`.
#   supplied     : TRUE when the caller supplied `missing`.
#   data         : the raw data, or NULL with a covariance matrix.
#   aux          : the auxiliary variable names, or NULL.
#   unavailable  : the options of the call that FIML cannot serve, each
#                  written as the caller would write it (for example
#                  ci_method = "bonett"); empty when there are none.
#   fiml_options : what the call can use under FIML instead, a phrase
#                  that completes the message.
# Returns:
#   "listwise", "fiml", or "two_stage". A default resolves to "fiml" only
#   for incomplete items.
.relia_resolve_missing <- function(missing, supplied, data, aux,
                                   unavailable = character(0),
                                   fiml_options = NULL) {
  if (supplied) return(missing)
  if (!is.null(aux)) return("fiml")
  if (is.null(data)) return("listwise")
  items <- suppressMessages(.relia_declared_missing(as.data.frame(data)))
  observed <- rowSums(!is.na(items)) > 0L
  if (!anyNA(items[observed, , drop = FALSE])) return("listwise")
  if (!length(unavailable)) return("fiml")
  .dmar_refuse_argument(c(
    "Full information maximum likelihood, the default of 'missing', is ",
    "not available with ", paste(unavailable, collapse = " and "),
    ", and the items have missing values, so no cases are deleted ",
    "unasked. Use missing = \"listwise\" to analyze the complete cases",
    if (!is.null(fiml_options)) paste0(", or ", fiml_options),
    "."), arg = "missing")
}

# .relia_resolve_inputs(data, S, N, allow_S = TRUE,
#                       missing = "listwise", aux = NULL, fun = NULL)
#
# Validate and normalize the (data, S, N) inputs that the reliability
# functions accept. The five family functions vary in which inputs they
# allow:
#
#   reliability_alpha       : data or (S, N)
#   reliability_kr20        : data only
#   reliability_omega       : data or (S, N); with denominator =
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
#                          behavior). Fewer than two complete rows, from
#                          which no variance can be computed, stop the
#                          call here, before any fit, with a message that
#                          names the cause and the treatment that uses the
#                          incomplete rows
#                          (.stop_relia_too_few_complete_rows()). Through
#                          DMAR 1.0.1 such data reached lavaan, which
#                          printed its table of the variables and stopped
#                          with an error the callers reported as a model
#                          that did not converge, and the analytic alpha
#                          returned NA with a warning that NaNs were
#                          produced. A coefficient fitted by the
#                          single-factor model (`single_factor` names it,
#                          "omega", "H", or "the model implied alpha";
#                          NULL for alpha and KR-20, which need only the
#                          sums of the matrix) needs more complete rows
#                          than items, since with no more the matrix is
#                          singular and the maximum likelihood fit is not
#                          defined, and is refused before any fit
#                          (.stop_relia_singular_complete_rows()).
#   missing = "fiml"     : rows with at least one observed item, or
#                          with auxiliary variables at least one observed
#                          item or auxiliary variable, are kept; N is that
#                          count and no complete-case S is computed (S is
#                          NULL; the estimation paths go through lavaan
#                          with missing = "ml" or through .fiml_cov()).
#                          reliability_omega_categorical() resolves its
#                          multiple imputation inputs this way, so a row
#                          observed only on an auxiliary is imputed.
#   missing = "two_stage": the same rows as under "fiml", for two-stage
#                          estimation (lavaan missing = "two.stage" or
#                          "robust.two.stage"), whose first stage reads a
#                          row observed on an auxiliary variable alone.
#                          When there are such rows, the message fa_confirmatory()
#                          and measurement_invariance() give
#                          (.ts_aux_only_note() in
#                          R/two_stage_internals.R) counts them, naming
#                          `fun`, the calling function; this helper runs
#                          once per call, before any bootstrap, so it is
#                          given once. Through the first versions of 1.1.0
#                          such a row was set aside.
#
# When `aux` names auxiliary columns of `data`, the remaining columns
# are the items; the auxiliaries are validated here (they must exist in
# `data`, be numeric, and not be constant). Wherever the columns are
# selected by name (under FIML and two-stage estimation, and with
# auxiliary variables), a name given to more than one column is refused
# (.stop_relia_duplicated_names()). When S is supplied, it is checked
# for symmetry.
#
# The items of a coefficient for continuous items (`continuous`, TRUE for
# every caller but reliability_omega_categorical(), whose items are
# ordered categories and are checked there) are checked before anything
# is computed from them, on every path. An item that is not numeric (or
# logical, which stats::cov() reads as 0 and 1) is refused by name, with
# its class; through DMAR 1.0.1 listwise deletion stopped with stats::cov's
# "is.numeric(x) || is.logical(x) is not TRUE" and FIML with a report that
# the saturated model did not converge. An item that takes a single value
# in the rows analyzed (the complete rows under listwise deletion, the
# rows that observe it otherwise), or a zero variance on the diagonal of
# `S`, is refused by name as well (.relia_check_constant_items()); lavaan
# stopped on it with an error that names no variable, after printing its
# table of them, and analytic alpha counted it among the items.
#
# Returns:
#   A list with elements (data, S, N, J, items, aux, missing,
#   N_complete), where exactly one of `data` and `S` is non-NULL
#   according to what the caller supplied. `items` names the item
#   columns of `data` (NULL for covariance input or unnamed data);
#   `N_complete` is the number of cases complete on the items (equal
#   to N under listwise deletion and for covariance input).
.relia_resolve_inputs <- function(data, S, N, allow_S = TRUE,
                                  missing = "listwise", aux = NULL,
                                  fun = NULL, continuous = TRUE,
                                  single_factor = NULL, needs_saturated = TRUE) {
  if (is.null(data) && is.null(S)) {
    .dmar_refuse_argument(c(
      "Either 'data' (a matrix or data frame of item scores) or ",
      "'S' (a covariance matrix) must be supplied."), arg = c("data", "S"))
  }
  if (!is.null(data) && !is.null(S)) {
    .dmar_refuse_argument(
      "Supply either 'data' or 'S', not both.", arg = c("data", "S"))
  }
  if (!is.null(S) && !allow_S) {
    .dmar_refuse_argument(c(
      "This coefficient requires raw 'data'; a covariance matrix is ",
      "not sufficient."), arg = "data")
  }

  if (!is.null(data)) {
    data <- .relia_declared_missing(as.data.frame(data))

    # The incomplete data paths, and every call with auxiliary variables,
    # select the columns by name, and a name given to two columns selects
    # the first of them: the second was dropped (with auxiliary variables,
    # whose items were found by setdiff(), which keeps one of each name,
    # so the refusal below the next block never saw the duplicate) or the
    # first analyzed twice, and the coefficient of the wrong items came
    # back without a word. Such a name is refused, whatever the columns
    # are. Listwise deletion without auxiliary variables reads the columns
    # by position, as it did through 1.0.1 (and the fitted coefficients
    # refuse the name in .omega_fit_cfa(), where a message could not tell
    # the two apart).
    # The fitted coefficients (`single_factor`) refuse such a name under
    # listwise deletion too (in .omega_fit_cfa()), so it is refused here,
    # before the refusals of the rows, which would otherwise come first and
    # offer a method that then refuses the name.
    if ((!is.null(aux) || missing %in% c("fiml", "two_stage") ||
         !is.null(single_factor)) &&
        anyDuplicated(colnames(data))) {
      .stop_relia_duplicated_names(colnames(data), aux)
    }

    if (!is.null(aux)) {
      if (!is.character(aux) || length(aux) < 1L || anyNA(aux)) {
        .dmar_refuse_argument(c(
          "'aux' must be a character vector of column names in ",
          "'data'."), arg = c("aux", "data"))
      }
      if (anyDuplicated(aux)) {
        .dmar_refuse_argument("'aux' contains duplicated names.", arg = "aux")
      }
      # as.data.frame() has already guaranteed column names (unnamed
      # matrix input gets V1, V2, ...), so a mismatch surfaces through
      # the not-found error below.
      absent <- setdiff(aux, colnames(data))
      if (length(absent) > 0L) {
        .dmar_refuse_argument(c(
          "Auxiliary variable(s) not found in 'data': ",
          paste(sQuote(absent), collapse = ", "),
          ". Supply 'data' containing both the items and the ",
          "auxiliary columns; the non-auxiliary columns are the ",
          "items."), arg = "data")
      }
      items <- colnames(data)[!colnames(data) %in% aux]
      if (length(items) < 2L) {
        .dmar_refuse_argument(c(
          "After removing the auxiliary column(s), fewer than two ",
          "item columns remain in 'data'. An auxiliary variable is ",
          "not an item; it cannot make up the composite."), arg = "data")
      }
      for (a in aux) {
        a_col <- data[[a]]
        if (!is.numeric(a_col)) {
          .dmar_refuse_argument(c(
            "Auxiliary variable ", sQuote(a), " is not numeric."), arg = NULL)
        }
        a_obs <- a_col[!is.na(a_col)]
        if (length(a_obs) < 2L || stats::var(a_obs) == 0) {
          .dmar_refuse_argument(c(
            "Auxiliary variable ", sQuote(a), " is constant (or has ",
            "fewer than two observed values); a constant carries no ",
            "information about the items or the missingness."), arg = NULL)
        }
      }
    } else {
      items <- colnames(data)
    }
    # Without auxiliary variables listwise deletion reads the columns by
    # position (two columns may share a name there), so every column is
    # checked as it stands.
    if (continuous) {
      .relia_check_numeric_items(if (is.null(aux)) data else data[items])
    }

    if (missing %in% c("fiml", "two_stage")) {
      if (is.null(items)) {
        # Unnamed input; the FIML paths need names for the lavaan
        # syntax, so assign the package's positional convention.
        colnames(data) <- paste0("y", seq_len(ncol(data)))
        items <- colnames(data)
      }
      # The items are selected by name; a name given to two columns has
      # been refused above.
      item_cols <- data[, items, drop = FALSE]
      if (ncol(item_cols) < 2L) {
        .dmar_refuse_argument("At least two items are required.", arg = NULL)
      }
      # The rows that carry information. With auxiliary variables a row
      # with every item missing still informs the item moments through an
      # auxiliary observed on it, because the saturated correlates model of
      # FIML, and the first stage of two-stage estimation, estimate the
      # joint distribution of the items and the auxiliaries (Graham, 2003;
      # Savalei & Bentler, 2009, whose first stage likelihood, Eq. 10, sums
      # over every case), and lavaan reads such a row. Setting it aside
      # conditions the analysis on having an observed item, which depends on
      # the auxiliaries whenever the whole item block is missing at random
      # given them, and the estimates are then biased: through DMAR 1.0.1
      # the row was set aside under FIML, and through the first versions of
      # 1.1.0 under two-stage estimation. lavaan's two-stage estimator
      # (lavaan 0.7-3) reads the row in its first stage but leaves it out of
      # its N; the family's two-stage fits go through fa_confirmatory(), which
      # computes them again with every case the first stage read
      # (.ts_all_rows()), and the message says so once per call.
      keep <- rowSums(!is.na(item_cols)) > 0L
      on_item <- keep
      with_aux <- missing %in% c("fiml", "two_stage") && !is.null(aux)
      if (with_aux) {
        on_aux <- rowSums(!is.na(data[, aux, drop = FALSE])) > 0L
        # The message is given here, once per call and before any fit, so
        # it says only what holds whatever the fit: the coefficient is
        # estimated from every case read. Its standard error counts them
        # once fa_confirmatory() has computed the two-stage fit again with every
        # case (.ts_all_rows()); when that could not be done, the family
        # says so in a warning of its own (.ts_coefficient_n_warning(),
        # from .omega_fit_cfa()).
        if (missing == "two_stage" && any(on_aux & !keep)) {
          .ts_aux_only_note(if (is.null(fun)) "reliability" else fun,
                            sum(on_aux & !keep), sum(keep | on_aux),
                            sum(keep),
                            what = paste0("the coefficient is estimated ",
                                          "from all ", sum(keep | on_aux),
                                          " cases it reads"))
        }
        keep <- keep | on_aux
      }
      data <- data[keep, , drop = FALSE]
      item_cols <- item_cols[keep, , drop = FALSE]
      on_item <- on_item[keep]
      if (nrow(data) < 2L) {
        .dmar_refuse_argument(
          if (with_aux) {
       "Fewer than two cases have an observed item or auxiliary variable."
     } else {
       "Fewer than two cases have at least one observed item."
     }, arg = NULL)
      }
      h1_diagnostic <- NULL
      if (continuous) {
        .relia_check_constant_items(item_cols, missing)
        # A composite that does not vary has no reliability, under any
        # treatment (.relia_check_composite()). With more complete rows than
        # items it is checked on them, exactly; with fewer, a fitted
        # coefficient checks the EM estimate (the analytic alpha does so in
        # .fiml_cov()). H weights its items, and is left to its fit.
        if (!identical(single_factor, "H")) {
          complete <- stats::complete.cases(item_cols)
          if (sum(complete) > ncol(item_cols)) {
            .relia_check_composite(
              stats::cov(item_cols[complete, , drop = FALSE]),
              rows = sum(complete), missing = missing)
          }
        }
        # Diagnose from every actual likelihood row, including auxiliaries.
        h1_diagnostic <- .fiml_data_h1_diagnostic(data[, c(items, aux), drop = FALSE])
        if (needs_saturated) {
          .stop_never_observed_together(data[, c(items, aux), drop = FALSE],
                                         .relia_incomplete_label(missing))
          .stop_single_valued_together(data[, c(items, aux), drop = FALSE],
                                        .relia_incomplete_label(missing))
          .stop_fiml_unbounded(data[, c(items, aux), drop = FALSE],
                                .relia_incomplete_label(missing))
        }
        if (needs_saturated && !identical(single_factor, "H") && !is.null(single_factor) &&
            sum(stats::complete.cases(item_cols)) <= ncol(item_cols)) {
          suppressWarnings(.fiml_cov(data, items = items, aux = aux,
                                     se = "none", missing = missing))
        }
      }
      list(data = data,
           S    = NULL,
           N    = nrow(data),
           J    = ncol(item_cols),
           items = items,
           aux  = aux,
           missing = missing,
           N_complete = sum(stats::complete.cases(item_cols)),
           h1_diagnostic = h1_diagnostic)
    } else {
      complete <- stats::complete.cases(data)
      if (ncol(data) < 2L) {
        .dmar_refuse_argument("At least two items are required.", arg = NULL)
      }
      if (sum(complete) < 2L) .stop_relia_too_few_complete_rows(data, needs_saturated = needs_saturated)
      all_rows <- data
      data <- data[complete, , drop = FALSE]
      # Distinct rows, since a repeated row adds nothing to the rank of the
      # covariance matrix (.relia_case_count()).
      if (!is.null(single_factor) &&
          sum(!duplicated(data)) <= ncol(data)) {
        .stop_relia_singular_complete_rows(all_rows, single_factor, needs_saturated = needs_saturated)
      }
      if (continuous) {
        .relia_check_constant_items(data, "listwise", all_rows = all_rows,
                                    composite = !identical(single_factor,
                                                           "H"),
                                    needs_saturated = needs_saturated)
        # A composite that does not vary (ipsative items) has no
        # reliability; H weights its items, and is left to its fit.
        if (!identical(single_factor, "H")) {
          .relia_check_composite(stats::cov(data))
        }
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
      .dmar_refuse_argument(
        "'N' must be supplied when 'S' is given.", arg = c("N", "S"))
    }
    if (!isSymmetric(unname(S), tol = 1e-6)) {
      .dmar_refuse_argument(
        "'S' must be a symmetric covariance matrix.", arg = "S")
    }
    if (ncol(S) < 2L) {
      .dmar_refuse_argument("At least two items are required.", arg = NULL)
    }
    .relia_check_constant_items(S = S)
    ev <- .relia_negative_eigenvalue(S)
    if (!is.null(ev)) .stop_relia_not_covariance(ev, S)
    if (!identical(single_factor, "H")) .relia_check_composite(S)
    list(data = NULL, S = S, N = N, J = ncol(S),
         items = colnames(S), aux = NULL, missing = "listwise",
         N_complete = N)
  }
}


# .relia_check_numeric_items(items)
#
# Refuses, by name and class, the item columns that are neither numeric
# nor logical (which stats::cov() and lavaan read as 0 and 1), before
# anything is computed from them, as auxiliary variables are refused in
# .relia_resolve_inputs(). Through DMAR 1.0.1 a character or factor item
# stopped listwise deletion with stats::cov's "is.numeric(x) ||
# is.logical(x) is not TRUE", which names neither the item nor the
# reason, and FIML with a report that the saturated model did not
# converge.
#
# Inputs:
#   items : data frame of the item columns.
# Returns:
#   TRUE, invisibly, or an error.
.relia_check_numeric_items <- function(items) {
  ok <- vapply(items, function(x) is.numeric(x) || is.logical(x),
               logical(1L))
  if (all(ok)) return(invisible(TRUE))
  bad <- names(items)[!ok]
  cls <- vapply(items[!ok], function(x) class(x)[1L], character(1L))
  .dmar_refuse_argument(c(
    "The items must be numeric, and ",
    if (length(bad) == 1L) "this column" else "these columns",
    " of 'data' ", if (length(bad) == 1L) "is" else "are", " not: ",
    paste0(sQuote(bad), " (", cls, ")", collapse = ", "), ". Recode ",
    "an item scored on a numeric scale as numeric; ordered categorical ",
    "items take reliability_omega_categorical()."), arg = "data")
}


# .relia_check_constant_items(item_data = NULL, missing = "listwise",
#                             all_rows = NULL, kr20 = FALSE, S = NULL,
#                             composite = TRUE)
#
# Refuses, by name, an item that takes a single value in the rows
# analyzed (the complete rows under listwise deletion, the rows that
# observe it under FIML and two-stage estimation), or whose variance on
# the diagonal of a supplied covariance matrix `S` is not positive, before
# anything is computed. Such an item adds a constant to the composite,
# which leaves the composite's variance and its reliability as they are,
# but a coefficient computed with it counts it among the J items:
# coefficient alpha (and KR-20) with it is J(J - 2)/(J - 1)^2 times alpha
# without it, the alpha of the same composite (a tenth lower at J = 4,
# since 8/9 of .9 is .8), and the single-factor model behind omega and H
# cannot be fit with it (lavaan stops with "some variables have no values
# (only missings) or no variance", naming none, after printing its table
# of the variables). An item no row observes is named as such. Under
# listwise deletion the message says when the item varies in the rows
# deletion removed, and offers full information maximum likelihood,
# which uses them (through reliability_alpha() for KR-20, kr20 = TRUE).
# reliability_omega_categorical() checks its ordered items itself
# (.check_omega_c_categories()), in the terms of their thresholds.
#
# Inputs:
#   item_data : data frame of the item columns of the rows analyzed (NULL
#               with `S`).
#   missing   : "listwise", "fiml", or "two_stage".
#   all_rows  : under listwise deletion, the item columns before deletion.
#   kr20      : TRUE when the caller is reliability_kr20().
#   S         : a covariance matrix supplied in place of data.
#   composite : FALSE for H, whose composite is weighted, so that full
#               information maximum likelihood does not refuse it for a
#               composite that does not vary (.relia_fiml_obstacle()).
# Returns:
#   TRUE, invisibly, or an error.
.relia_check_constant_items <- function(item_data = NULL,
                                        missing = "listwise",
                                        all_rows = NULL, kr20 = FALSE,
                                        S = NULL, composite = TRUE, needs_saturated = TRUE) {
  consequence <- paste0(
    "An item that does not vary adds a constant to the composite, which ",
    "leaves its variance and its reliability as they are, but a ",
    "coefficient computed with it counts it among the J items (alpha ",
    "with one such item is J(J - 2)/(J - 1)^2 times alpha without it), ",
    "and the single-factor model cannot be fit with it.")
  # "a", "a and b", "a, b, and c".
  listed <- function(x) {
    if (length(x) <= 2L) return(paste(x, collapse = " and "))
    paste0(paste(x[-length(x)], collapse = ", "), ", and ", x[length(x)])
  }
  if (!is.null(S)) {
    S <- as.matrix(S)
    v <- diag(S)
    nm <- colnames(S)
    if (is.null(nm)) nm <- paste("item", seq_along(v))
    bad <- !(is.finite(v) & v > 0)
    if (!any(bad)) return(invisible(TRUE))
    zero <- bad & is.finite(v) & v == 0
    .dmar_refuse_argument(
      if (all(zero[bad])) {
   paste0("The variance of ", listed(nm[bad]),
          " in 'S' is zero. ", consequence, " Leave ",
          if (sum(bad) == 1L) "it" else "them", " out of 'S'.")
 } else {
   paste0("The diagonal of 'S' must hold the item variances, which ",
          "are positive, and ",
          if (sum(bad) == 1L) "the entry for " else "the entries for ",
          listed(nm[bad]), if (sum(bad) == 1L) " is" else " are",
          " not (", paste(trimws(format(v[bad], digits = 3L)),
                          collapse = ", "), ").")
 }, arg = NULL)
  }
  item_data <- as.data.frame(item_data)
  n_obs <- colSums(!is.na(item_data))
  unobserved <- names(item_data)[n_obs == 0L]
  if (length(unobserved)) {
    .dmar_refuse_argument(c(
      "No row observes ", paste(unobserved, collapse = " or "),
      ", so neither its variance nor its covariances with the other ",
      "items can be estimated, under listwise deletion or any other ",
      "treatment. Leave it out of 'data'."), arg = "data")
  }
  # By position, since listwise deletion without auxiliary variables
  # allows two columns one name.
  at <- which(.relia_single_valued(item_data))
  single <- names(item_data)[at]
  if (!length(single)) return(invisible(TRUE))
  who <- if (length(single) == 1L) {
    paste(single, "takes")
  } else {
    paste(listed(single), "each take")
  }
  where <- if (identical(missing, "listwise")) {
    paste0("the ", nrow(item_data), " complete rows")
  } else if (length(single) == 1L) {
    paste0("the ", n_obs[[at]], " ",
           if (n_obs[[at]] == 1L) "row that observes" else
             "rows that observe", " it")
  } else {
    "the rows that observe them"
  }
  varies <- identical(missing, "listwise") && !is.null(all_rows) &&
    !all(.relia_single_valued(as.data.frame(all_rows)[at]))
  .dmar_refuse_argument(c(
    who, " a single value in ", where, ", so ",
    if (length(single) == 1L) "its variance is" else
      "their variances are", " zero. ", consequence,
    if (varies) {
      # FIML is offered only when it could use the rows deletion
      # removed (.relia_fiml_obstacle()): the variation may lie only in
      # rows that miss another item.
      obstacle <- .relia_fiml_obstacle(all_rows, kr20 = kr20,
                                       composite = composite,
                                       needs_saturated = needs_saturated)
      paste0(" Listwise deletion analyzes only the complete rows, and ",
             if (length(single) == 1L) "the item varies" else
               "these items vary",
             " in the rows it removes",
             if (is.null(obstacle)) {
               paste0("; ", .relia_fiml_words(kr20),
                      " uses those rows as well.",
                      .relia_fiml_precondition(colnames(all_rows)))
             } else {
               paste0(". ", obstacle)
             })
    } else {
      paste0(" Leave ", if (length(single) == 1L) "it" else "them",
             " out of 'data'.")
    }), arg = NULL)
}


# .relia_single_valued(x)
#
# Which columns of `x` (a data frame or matrix of numeric or logical
# columns, possibly with NA) take fewer than two distinct observed
# values: TRUE for a column every observed value of which equals the
# first, and for a column with no observed value. Each column is compared
# with its first observed value in one vectorized pass, so the check costs
# little more than reading the data once; it runs on the data and, through
# .relia_check_resample(), on every bootstrap replication, where a call
# to unique() per column would cost several times what the analytic alpha
# of the replication does. The result is by position, since listwise
# deletion allows two columns one name.
#
# Inputs:
#   x : data frame or matrix.
# Returns:
#   A logical vector, one element per column.
.relia_single_valued <- function(x) {
  n <- NROW(x)
  if (is.data.frame(x)) {
    x <- matrix(vapply(x, as.numeric, numeric(n)), nrow = n)
  } else {
    x <- matrix(as.numeric(x), nrow = n)
  }
  p <- ncol(x)
  if (n == 0L) return(rep(TRUE, p))
  obs <- !is.na(x)
  # The row of each column's first observed value (row 1 for a column with
  # none, whose value is then NA and compares as NA throughout).
  first <- max.col(t(obs), ties.method = "first")
  ref <- x[cbind(first, seq_len(p))]
  colSums(x != rep(ref, each = n), na.rm = TRUE) == 0
}


# .relia_check_resample(d, missing = "listwise", items = NULL, aux = NULL,
#                       fitted = FALSE, composite = FALSE, on_fit = NULL)
#
# The refusals of the data, applied to a bootstrap replication as
# .relia_resolve_inputs() applies them to the data: a resample, or a
# jackknife deletion behind the BCa acceleration, is analyzed as data, so
# it follows the same rules, and the replication raises the error, which
# the bootstrap engine counts and quotes, and the caution reports the
# share set aside.
#
# An item (or auxiliary variable) that does not vary. A resample can lose
# every case that gave an item its variation (a rare 0/1 item passed by
# one respondent loses it in about a third of the resamples, since a case
# is left out of a resample with probability (1 - 1/N)^N, near exp(-1)).
# The coefficient is not defined for such an item under the rule the call
# itself follows (.relia_check_constant_items()). Through version 1.0.1
# (and in the first versions of 1.1.0, which checked the data alone) the
# analytic alpha and KR-20 of such a replication counted the item among
# the J items, which is J(J - 2)/(J - 1)^2 times the coefficient of the
# same composite without it, and kept the replication; the fitted
# coefficients had lavaan stop on the item, printing its table of the
# variables to the console for every such replication, and counted it as
# a replication with no value. Under listwise deletion every column is an
# item and every row complete; under full information maximum likelihood
# and two-stage estimation the items and the auxiliary variables are
# checked on the rows that observe them, and a pair that the resample
# observes together only where one of them takes a single value is
# refused by the fit itself (.stop_single_valued_together() in
# R/cov_fiml.R). A column no row of the resample observes is refused
# first, in the data's words ("No row observes"): its variance is not
# identified, not zero.
#
# Too few distinct cases. Under full information maximum likelihood and
# two-stage estimation the resample needs more distinct cases with an
# observed item than the J + K variables of the saturated model
# (.relia_case_count(), .stop_relia_too_few_cases()), and under listwise
# deletion a coefficient fitted by the single-factor model (`fitted`)
# needs more distinct rows than items
# (.stop_relia_singular_complete_rows()). A resample repeats rows, so at
# a sample size just above the rule's bound most resamples fall below it
# (with 12 cases and J + K = 8, three resamples in four hold no more than
# 8 distinct cases). Through the first versions of 1.1.0 the rule was
# applied to the data alone, and such a replication was analyzed: its EM
# estimate was singular or nearly so, and its value entered the interval.
# The refusal is an error of class dmar_too_few_cases
# (.stop_resample_too_few_cases()), and when more than half of the
# replications are set aside by it the engine forms no interval: those it
# keeps are the resamples that repeat the fewest rows, chosen by their
# composition, and at the bound every one of them holds each distinct case
# (.bootstrap_ci()).
#
# The bootstrap of a lavaan fit (reliability_H(fit = ); .H_refit_fun())
# applies the same checks to the model's continuous variables, under the
# rule of the fit's own treatment of incomplete rows (`missing` is
# "fiml" for lavaan's "ml" and "ml.x", "two_stage" for its "two.stage"
# and "robust.two.stage", "listwise" for "listwise", and lavaan's own
# label otherwise), and says so (`on_fit`). A listwise fit, which includes
# a fit to complete data, is refitted to the sample covariance matrix of
# the resample, whose rank is at most the distinct rows less one, so with
# no more distinct rows than variables it is singular; the rule of the
# saturated model applies only under full information maximum likelihood
# and two-stage estimation; and a fit by pairwise deletion, whose
# covariances come from different rows, is held to neither rule, since
# its matrix need not be singular with few rows. reliability_H(fit = )
# does not refuse the fit's data by these rules, so the messages of that
# path describe the resample alone. Through the first versions of 1.1.0
# every fit was checked as a FIML fit, and a replication of a fit to
# complete data was refused in the words of the saturated model that full
# information maximum likelihood estimates first, with "the data would be
# refused for the same reason".
#
# Inputs:
#   d       : the replication's rows.
#   missing : "listwise", "fiml", or "two_stage" (on the path of a fit,
#             also lavaan's label for another treatment).
#   items   : the item columns (NULL: every column, by position).
#   aux     : the auxiliary columns, or NULL.
#   fitted  : TRUE for a coefficient fitted by the single-factor model.
#   composite : TRUE to refuse a composite that does not vary in the
#           complete rows of the resample (.relia_check_composite()), when
#           they number more than the items, for the coefficients whose
#           point function does not apply the check exactly itself (the
#           fitted coefficients, and the analytic alpha under FIML and
#           two-stage estimation, whose EM estimate is held to a looser
#           tolerance; not H, whose composite is weighted).
#   on_fit  : NULL for the data of a call, or, for the bootstrap of a
#             lavaan fit, the words for the variables checked
#             ("variables", or "continuous variables" when the model also
#             has ordered ones).
# Returns:
#   TRUE, invisibly, or an error.
.relia_check_resample <- function(d, missing = "listwise", items = NULL,
                                  aux = NULL, fitted = FALSE,
                                  composite = FALSE, on_fit = NULL) {
  cols <- if (is.null(items)) seq_len(NCOL(d)) else c(items, aux)
  x <- if (is.data.frame(d)) d[cols] else d[, cols, drop = FALSE]
  labels <- colnames(x)
  if (is.null(labels)) labels <- paste("item", seq_len(NCOL(x)))
  fit_path <- !is.null(on_fit)
  # A column no row of the resample observes is named as the data's is
  # (.relia_check_constant_items()): its variance and covariances are not
  # identified, rather than zero. .relia_single_valued() flags it too, and
  # through the first versions of 1.1.0 it was reported as taking a single
  # value, and tallied with the columns that do.
  none <- colSums(!is.na(x)) == 0L
  if (any(none)) {
    nm <- labels[none]
    one <- length(nm) == 1L
    .dmar_refuse_argument(c(
      "No row of this resample observes ",
      if (one) nm else paste(paste(nm[-length(nm)], collapse = ", "),
                             "or", nm[length(nm)]),
      ", so neither ", if (one) "its variance nor its" else
        "their variances nor their", " covariances with the other ",
      if (fit_path || all(nm %in% aux)) "variables" else "items",
      " can be estimated",
      if (fit_path) "." else
        ", and the data would be refused for the same reason."), arg = NULL)
  }
  flat <- .relia_single_valued(x)
  if (any(flat)) {
    nm <- labels[flat]
    one <- length(nm) == 1L
    what <- if (fit_path) {
      "a variable that does not vary is"
    } else if (all(nm %in% aux)) {
      "an auxiliary variable that does not vary is"
    } else {
      "an item that does not vary is"
    }
    .dmar_refuse_argument(c(
      if (one) nm else paste(paste(nm[-length(nm)], collapse = ", "),
                             "and", nm[length(nm)]),
      if (one) " takes" else " each take", " a single value in ",
      if (identical(missing, "listwise")) {
        "this resample"
      } else {
        paste0("the rows of this resample that observe ",
               if (one) "it" else "them")
      },
      ", so ", if (one) "its variance is" else "their variances are",
      " zero, and ", what, " refused."), arg = NULL)
  }
  if (any(missing == c("fiml", "two_stage")) && !is.null(items)) {
    p <- length(items) + length(aux)
    cases <- .relia_case_count(x, items, aux)
    if (cases$n <= p) {
      .stop_resample_too_few_cases(paste0(
        "Only ", cases$n, " distinct cases of this resample have an observed ",
        if (fit_path) "model variable" else "item", ", no more than the ", p,
        " model variables. It does not meet DMAR's distinct-case estimability ",
        "screen for incomplete-data resamples. This conservative screen ",
        "does not establish likelihood nonexistence."))
    }
  } else if (identical(missing, "listwise") && fitted) {
    n <- sum(!duplicated(as.data.frame(x)))
    if (n <= NCOL(x)) {
      .stop_resample_too_few_cases(if (fit_path) {
        paste0("Only ", n, " distinct rows of this resample, no more than ",
               "the ", NCOL(x), " ", on_fit, " of the model, so their sample ",
               "covariance matrix is singular (its rank is at most ",
               n - 1L, ").")
      } else {
        paste0("Only ", n, " distinct rows of this resample, no more than ",
               "the ", NCOL(x), " items, so their covariance matrix is ",
               "singular and the single-factor model cannot be fit to it, ",
               "and the data would be refused for the same reason.")
      })
    }
  }
  if (composite) {
    y <- if (is.null(items)) x else x[items]
    complete <- stats::complete.cases(y)
    if (sum(complete) > NCOL(y)) {
      .relia_check_composite(
        stats::cov(y[complete, , drop = FALSE]),
        rows = if (!identical(missing, "listwise")) sum(complete),
        missing = missing)
    }
  }
  invisible(TRUE)
}


# .stop_resample_too_few_cases(message)
#
# The refusal of a resample with too few distinct cases
# (.relia_check_resample()), an error of class dmar_too_few_cases, which
# the bootstrap engine counts apart (`too_few_cases`): when more than half
# of the replications are set aside for it, no interval is formed, since
# those it keeps are chosen by their composition (see .bootstrap_ci()).
.stop_resample_too_few_cases <- function(message) {
  .dmar_refuse(
    message,
    "too_few_cases")
}


# .relia_checked_point(point_fn, missing = "listwise", items = NULL,
#                      aux = NULL, fitted = FALSE)
#
# A bootstrap point function of the reliability family with the checks of
# .relia_check_resample() in front of it: the function .bootstrap_ci()
# applies to every resample and jackknife deletion. Every continuous
# coefficient's bootstrap goes through it (reliability_alpha(),
# reliability_omega(), reliability_kr20(), and reliability_H() from
# data). Categorical omega uses the same column checks, without the
# continuous fit's case-count or composite checks, including on completed
# data sets under multiple imputation. A lavaan fit's bootstrap in
# reliability_H() checks its variables in .H_refit_fun().
#
# Inputs:
#   point_fn : function(d) returning the coefficient of the rows d.
#   missing, items, aux, fitted, composite : as in
#            .relia_check_resample().
#   needs_saturated : whether the statistic requires unrestricted moments;
#            screen every selected resample, independently of its parent.
# Returns:
#   A function(d).
.relia_checked_point <- function(point_fn, missing = "listwise",
                                 items = NULL, aux = NULL, fitted = FALSE,
                                 composite = FALSE, needs_saturated = FALSE) {
  force(point_fn)
  force(missing)
  force(items)
  force(aux)
  force(fitted)
  force(composite)
  force(needs_saturated)
  function(d) {
    .relia_check_resample(d, missing = missing, items = items, aux = aux,
                          fitted = fitted, composite = composite)
    if (needs_saturated && missing %in% c("fiml", "two_stage")) {
      cols <- if (is.null(items)) names(d) else c(items, aux)
      .stop_never_observed_together(d[, cols, drop = FALSE],
                                     .relia_incomplete_label(missing))
      .stop_single_valued_together(d[, cols, drop = FALSE],
                                    .relia_incomplete_label(missing))
      .stop_fiml_unbounded(d[, cols, drop = FALSE],
                           .relia_incomplete_label(missing))
    }
    point_fn(d)
  }
}


# .relia_check_composite(S, em = FALSE, rows = NULL, missing = NULL)
#
# Refuses, before anything is fitted, a covariance matrix whose composite
# (the unweighted sum of the items) does not vary
# (.relia_composite_varies()), with the message of
# .stop_relia_constant_composite(). It is applied to the complete rows
# under listwise deletion, to a supplied S, and under full information
# maximum likelihood and two-stage estimation to the rows that observe
# every item when they number more than the items (`rows`), and otherwise,
# for the fitted coefficients, to the EM estimate (`em`), in
# .relia_resolve_inputs(), and to every bootstrap replication of the
# composite's coefficients (.relia_check_resample()); .alpha_from_S()
# applies it to the matrix alpha is computed from. Through the first
# versions of 1.1.0 it was applied only to the complete rows under
# listwise deletion and to S: under FIML and two-stage estimation the
# model implied alpha of ipsative items came back as about 1e-13 with an
# interval, two-stage omega as 5.6e-6 or a model that did not converge,
# and with more values missing the analytic alpha as -5e7.
.relia_check_composite <- function(S, em = FALSE, rows = NULL,
                                   missing = NULL) {
  if (!.relia_composite_varies(S, em = em)) {
    # A sum not above the tolerance proves only that the composite has no
    # variance if S is a covariance matrix. A supplied S may not be one
    # (a pairwise-deletion matrix, a typing error), and then a sum below
    # zero is a negative variance, not a zero one, and the matrix need not
    # be singular; it is named as what it is (.relia_negative_eigenvalue()).
    ev <- .relia_negative_eigenvalue(S)
    if (!is.null(ev)) .stop_relia_not_covariance(ev, S)
    .stop_relia_constant_composite(em = em, rows = rows, missing = missing)
  }
  invisible(TRUE)
}


# .relia_negative_eigenvalue(S)
#
# The smallest eigenvalue of S when it is below zero by more than rounding
# error, -sqrt(.Machine$double.eps) times the largest absolute eigenvalue
# (or times one, when they are all smaller), and NULL otherwise: the rule by
# which ss_aipe_reliability() refuses a population that is not a covariance
# matrix (.ss_aipe_check_S() in R/ss_aipe_reliability.R), so the planner
# and the estimation functions draw the line in one place. The covariance
# matrix of ipsative items has a smallest eigenvalue of rounding size
# (about -1e-14 for three holzinger_swineford tests), which is a covariance
# matrix's.
.relia_negative_eigenvalue <- function(S) {
  S <- unname(as.matrix(S))
  # Standardization preserves definiteness under arbitrary item units.
  # A tolerance on the raw largest eigenvalue can conceal an invalid
  # correlation when item variances differ by many orders of magnitude.
  d <- diag(S)
  standardized <- all(is.finite(d)) && all(d > 0)
  if (standardized) {
    scales <- sqrt(d)
    checked <- S / outer(scales, scales)
  } else {
    scale <- max(abs(S))
    if (scale == 0) return(NULL)
    checked <- S / scale
  }
  ev <- eigen(checked, symmetric = TRUE, only.values = TRUE)$values
  if (min(ev) >= -sqrt(.Machine$double.eps) * max(abs(ev))) return(NULL)
  out <- min(ev)
  attr(out, "standardized") <- standardized
  out
}


# .stop_relia_not_covariance(ev, S)
#
# The refusal of a supplied S that is not a covariance matrix, reached when
# its composite variance, sum(S), is not above the tolerance of
# .relia_composite_varies(). Through the first versions of 1.1.0 such a
# matrix was refused with the message of .stop_relia_constant_composite(),
# which says the variance is zero up to rounding error and the matrix
# singular: for three items that correlate -.9 the sum is -2.4 and the
# matrix is indefinite (eigenvalues 1.9, 1.9, and -0.8) and not singular.
# Version 1.0.1 returned alpha 3.38 with the limits [4.31, 1] for it. The
# composite's variance is quoted when it is below zero beyond rounding
# error.
.stop_relia_not_covariance <- function(ev, S) {
  num <- function(x) trimws(formatC(x, digits = 3, format = "fg"))
  total <- sum(S)
  .dmar_refuse_argument(c(
    "'S' is not a covariance matrix: its smallest eigenvalue ",
    if (isTRUE(attr(ev, "standardized"))) "after standardizing item units " else "after common rescaling ", "is ", num(ev), ", so some weighted sum of the items would have a negative ",
    "variance, which no population has",
    if (total < -sqrt(.Machine$double.eps) * sum(abs(diag(S)))) {
      paste0(", and the composite, the sum of the items, would have the ",
             "variance ", num(total), " (the sum of the entries of 'S')")
    }, "."), arg = "S")
}


# .stop_relia_duplicated_names(names, aux = NULL)
#
# The refusal of a name given to more than one column, where the columns
# are selected by name (.relia_resolve_inputs(), .omega_fit_cfa(), and
# .omega_fit_cfa_aux()). A name `aux` lists is named as an auxiliary
# variable's, since the two columns may be an auxiliary variable and an
# item.
#
# Inputs:
#   names : the column names, with at least one duplicated.
#   aux   : the names of the auxiliary variables, or NULL.
# Returns:
#   Nothing; it stops.
.stop_relia_duplicated_names <- function(names, aux = NULL) {
  dup <- unique(names[duplicated(names)])
  in_aux <- dup %in% aux
  .dmar_refuse_argument(c(
    if (any(in_aux)) "The column names of 'data'" else "The item names",
    " must be distinct; ", paste(sQuote(dup), collapse = ", "),
    if (length(dup) == 1L) " names" else " each name",
    " more than one column",
    if (any(in_aux)) {
      paste0(" (", paste(sQuote(dup[in_aux]), collapse = ", "),
             " named in 'aux')")
    },
    "."), arg = NULL)
}


# .stop_relia_too_few_complete_rows(data, kr20 = FALSE)
#
# The refusal when listwise deletion leaves fewer than two rows
# (.relia_resolve_inputs(), and reliability_kr20(), which deletes on its
# own), from which no variance can be computed: lavaan printed its table
# of the variables and stopped on them ("some variables have no values"
# or "only 1 observation"), which the family reported as a model that
# did not converge, and the analytic alpha and KR-20 returned NA with a
# warning that NaNs were produced. The message names the cause and the
# treatment that uses the incomplete rows, as
# reliability_omega_categorical() does for its own items
# (.stop_omega_c_no_complete_rows()). An item no row observes is named
# first, since no treatment can estimate its variance; full information
# maximum likelihood is offered only when every pair of items is
# observed together somewhere, since without that it has no data for
# some covariance either (.stop_never_observed_together() in
# R/cov_fiml.R). With two or more complete rows coefficient alpha and
# KR-20 are computed, since they need only the sums of the covariance
# matrix; the coefficients fitted by the single-factor model need more
# complete rows than items, and with no more they are refused in
# .relia_resolve_inputs() (.stop_relia_singular_complete_rows()).
#
# reliability_kr20() has no `missing` argument, so for it (kr20 = TRUE)
# the message offers full information maximum likelihood through
# reliability_alpha(), which computes the same coefficient (alpha of 0/1
# items is KR-20, and alpha does not depend on the divisor of the
# covariance matrix) and takes missing = "fiml"; the rest of the message
# is the family's.
#
# Inputs:
#   data : the item columns before deletion (listwise deletion has no
#          auxiliary variables).
#   kr20 : TRUE when the caller is reliability_kr20().
# Returns:
#   Nothing; it stops.
.stop_relia_too_few_complete_rows <- function(data, kr20 = FALSE, needs_saturated = TRUE) {
  data <- as.data.frame(data)
  obs <- !is.na(data)
  nm <- colnames(data)
  unobserved <- nm[colSums(obs) == 0L]
  if (length(unobserved)) {
    .dmar_refuse_argument(c(
      "No row observes ", paste(unobserved, collapse = " or "),
      ", so neither its variance nor its covariances with the other ",
      "items can be estimated, under listwise deletion or any other ",
      "treatment. Leave it out of 'data'."), arg = "data")
  }
  n_complete <- sum(rowSums(obs) == ncol(obs))
  remedy <- .relia_incomplete_remedy(data, kr20 = kr20, needs_saturated = needs_saturated)
  # Listwise deletion was the family's default through DMAR 1.0.x, and the
  # message called it that; since 1.1.0 the default analyzes incomplete
  # items by FIML (.relia_resolve_missing()), so the message names no
  # default.
  .dmar_refuse_argument(c(
    if (n_complete == 0L) "No row" else "Only one row", " of 'data' ",
    "is observed on every item, so listwise deletion ",
    if (kr20) "(the only treatment reliability_kr20() offers) ",
    "leaves nothing to analyze. ", remedy), arg = "data")
}


# .relia_incomplete_remedy(data, kr20 = FALSE, composite = TRUE)
#
# The sentence a refusal of listwise deletion ends with: full information
# maximum likelihood, which uses the incomplete rows as well, or, when
# the data give it nothing to estimate from, why it cannot help either
# (.relia_fiml_obstacle()). For reliability_kr20() (kr20 = TRUE), which
# has no `missing` argument, it is offered through reliability_alpha().
# Shared by .stop_relia_too_few_complete_rows() and
# .stop_relia_singular_complete_rows().
#
# Inputs:
#   data : the item columns before deletion.
#   kr20 : TRUE when the caller is reliability_kr20().
#   composite : FALSE for H (see .relia_fiml_obstacle()).
# Returns:
#   A sentence.
.relia_incomplete_remedy <- function(data, kr20 = FALSE, composite = TRUE, needs_saturated = TRUE) {
  obstacle <- .relia_fiml_obstacle(data, kr20 = kr20, composite = composite, needs_saturated = needs_saturated)
  if (!is.null(obstacle)) return(obstacle)
  fiml <- .relia_fiml_words(kr20)
  paste0(if (kr20) {
    paste0("Use ", fiml, ", which computes the same coefficient (alpha ",
           "of 0/1 items is KR-20) from the incomplete rows as well and ",
           "assumes the values missing at random, where listwise deletion ",
           "assumes them missing completely at random.")
  } else {
    paste0("Use ", fiml, ", which uses the incomplete rows as ",
           "well and assumes the values missing at random, where listwise ",
           "deletion assumes them missing completely at random.")
  }, .relia_fiml_precondition(colnames(data)))
}


# .relia_fiml_precondition(names)
#
# The sentence an offer of full information maximum likelihood adds when
# two columns of 'data' share a name, or "": that method selects the items
# by name and refuses such a name before anything else
# (.stop_relia_duplicated_names(), in .relia_resolve_inputs()), while
# listwise deletion without auxiliary variables reads the columns by
# position. The data are not the obstacle (once the columns are named
# apart, FIML uses the rows listwise deletion removed), so the offer stands
# with this condition rather than giving way to "cannot help either".
# Through the first versions of 1.1.0 the offer was made without it, and
# FIML then refused the call on the name.
.relia_fiml_precondition <- function(names) {
  dup <- unique(names[duplicated(names)])
  if (!length(dup)) return("")
  paste0(" It selects the items by name, so first give the columns named ",
         paste(sQuote(dup), collapse = ", "), " distinct names (or leave ",
         "out a column that copies another).")
}


# .relia_fiml_words(kr20 = FALSE)
#
# How a refusal of listwise deletion names full information maximum
# likelihood: missing = "fiml", or, for reliability_kr20(), which has no
# `missing` argument, reliability_alpha() with missing = "fiml" (alpha of
# 0/1 items is KR-20).
.relia_fiml_words <- function(kr20 = FALSE) {
  if (kr20) {
    "reliability_alpha() with missing = \"fiml\""
  } else {
    "missing = \"fiml\""
  }
}


# .relia_fiml_obstacle(data, kr20 = FALSE, composite = TRUE)
#
# A detected obstacle to full information maximum likelihood, or NULL
# when no implemented obstacle was found: the sentence a refusal of
# listwise deletion gives in place of the offer of FIML, so that it never
# offers a method the same data would be refused under. The obstacles
# are the ones the data put in FIML's way, in the order it meets them
# (after a name given to two columns, which it refuses first but which is
# not an obstacle of the data: the offer then carries the condition of
# .relia_fiml_precondition() instead): an item that
# takes a single value in every row that observes it
# (.relia_check_constant_items()), no more cases than items
# (.stop_relia_too_few_cases()), a composite that does not vary in the
# rows that observe every item when they outnumber the items
# (.relia_check_composite(), for every coefficient but H, whose composite
# is weighted: `composite`), a pair of items no row observes together
# (.stop_never_observed_together() in R/cov_fiml.R), and a pair observed
# together only where one of them takes a single value
# (.stop_single_valued_together()). The last is the one the refusal of a
# constant item met when it offered FIML unconditionally: an item
# constant in the complete rows whose variation lies only in rows that
# miss another item was offered FIML, which returned robust omega 1.11
# for the four tests of the example in .single_valued_pairs(). The
# composite is the one it met once FIML refused a composite that does not
# vary: with x1 constant in 20 complete rows in which x2 + x3 + x4 is
# zero, and 30 more rows that each miss one of x2, x3, and x4, the refusal
# of the constant item offered FIML, which then refused the composite
# (through the first versions of 1.1.0 FIML returned alpha -0.586 for
# these data). With no more complete rows than
# items, FIML checks the composite on its EM estimate instead
# (.relia_composite_varies(em = TRUE)), which is not mirrored here, since
# it would fit the EM inside a refusal; on three to five ipsative
# holzinger_swineford tests with one, two, or J complete rows among 301,
# the others each missing one or more tests, the EM left the composite's
# variance at 1e-5 to 1e-2 of the sum of the item variances, above the
# tolerance of 1e-6, so FIML did not refuse them.
#
# Inputs:
#   data : the item columns before deletion (listwise deletion has no
#          auxiliary variables).
#   kr20 : TRUE when the caller is reliability_kr20().
#   composite : FALSE for H, whose composite is weighted.
# Returns:
#   A sentence, or NULL.
.relia_fiml_obstacle <- function(data, kr20 = FALSE, composite = TRUE, needs_saturated = TRUE) {
  data <- as.data.frame(data)
  nm <- colnames(data)
  # Where two columns share a name, each is named with its position, so a
  # pair or item named below is not ambiguous.
  dup <- nm %in% nm[duplicated(nm)]
  nm[dup] <- paste0(nm[dup], " (column ", which(dup), ")")
  lead <- paste0("Full information maximum likelihood (",
                 .relia_fiml_words(kr20), ") cannot help either: ")
  flat <- which(.relia_single_valued(data))
  if (length(flat)) {
    return(paste0(lead, nm[flat[1L]], " takes a single value in every ",
                  "row that observes it, so its variance is zero there ",
                  "too."))
  }
  if (composite) {
    complete <- stats::complete.cases(data)
    k <- sum(complete)
    if (k > ncol(data) &&
        !.relia_composite_varies(stats::cov(data[complete, , drop = FALSE]))) {
      return(paste0(lead, "the composite, the sum of the items, does not ",
                    "vary in the ", k, " rows that observe every item, so ",
                    "the likelihood of the saturated model it estimates ",
                    "first has no maximum."))
    }
  }
  likelihood <- .fiml_likelihood_obstacle(data)
  if (needs_saturated && !is.null(likelihood)) return(paste0(lead, likelihood))
  joint <- crossprod(!is.na(data))
  never <- which(joint == 0 & upper.tri(joint), arr.ind = TRUE)
  if (nrow(never) > 0L) {
    return(paste0(lead, "no row observes both ", nm[never[1L, 1L]],
                  " and ", nm[never[1L, 2L]], ", so the data carry no ",
                  "information about their covariance; ",
                  "missing_coverage() reports the cases behind every ",
                  "pair."))
  }
  pairs <- .single_valued_pairs(data)
  if (nrow(pairs) > 0L) {
    return(paste0(lead, "it would be refused for the pair ",
                  .single_valued_pair_phrase(pairs, nm, 1L), ", whose ",
                  "rows in common carry no information about how the two ",
                  "vary together; missing_coverage() reports the cases ",
                  "behind every pair."))
  }
  NULL
}


# .relia_case_count(d, items = NULL, aux = NULL)
#
# The cases that count toward the rule of .stop_relia_too_few_cases():
# the distinct rows of `d` (over the items and the auxiliary variables)
# among those that observe at least one item. A row observed on an
# auxiliary variable alone does not count. The saturated model factors
# into the distribution of the auxiliary variables and that of the items
# given them, f(a) f(y | a), and such a row enters only the first
# (Anderson, 1957; Little & Rubin, 2019, ch. 7): it informs the means and covariances of the
# auxiliary variables, and through them the item moments, but not the
# regression of the items on the auxiliary variables, whose residual
# covariance matrix the rows that observe an item must make nonsingular.
# In the subspace argument of .stop_relia_too_few_cases(), its missing
# items can be filled in to put it on any hyperplane through the other
# cases that is not parallel to the item axes, so it cannot bound the
# likelihood. A repeated row adds no dimension to the space the cases
# span, and the likelihood of data with a row repeated is the likelihood
# of the distinct rows with that row weighted, whose supremum is the
# same. Through the first versions of 1.1.0 every row read was counted:
# five pupils of three holzinger_swineford tests with a fourth as an
# auxiliary variable were refused, and a sixth pupil observed on the
# auxiliary alone turned the refusal into an analytic alpha of 0.418 with
# the interval [0.418, 0.418] under two-stage estimation (and under FIML);
# the same five pupils with every row given twice were refused once and
# analyzed the second time. Distinct rows give an upper bound on the
# dimension the cases span, not the dimension itself (a row that another
# completes, as (1, NA) is completed by (1, 2), counts as distinct), so
# the rule remains a necessary condition, as before.
#
# Inputs:
#   d     : the rows (a data frame, or a matrix with column names).
#   items : the item columns (NULL: every column).
#   aux   : the auxiliary columns, or NULL.
# Returns:
#   A list with `n`, the distinct rows that observe an item, and `rows`,
#   the rows that observe an item.
.relia_case_count <- function(d, items = NULL, aux = NULL) {
  if (!is.data.frame(d)) d <- as.data.frame(d)
  y <- if (is.null(items)) d else d[items]
  on_item <- rowSums(!is.na(y)) > 0L
  x <- if (is.null(items)) d else d[c(items, aux)]
  x <- x[on_item, , drop = FALSE]
  list(n = sum(!duplicated(x)), rows = nrow(x))
}


# .relia_listwise_obstacle(data)
#
# Why coefficient alpha under listwise deletion would refuse the complete
# rows of `data` (the item columns), or NULL when it would compute: fewer
# than two complete rows (.stop_relia_too_few_complete_rows()), an item
# that takes a single value in them (.relia_check_constant_items()), or a
# composite that does not vary in them (.relia_check_composite()), in the
# order .relia_resolve_inputs() checks them. It is the mirror of
# .relia_fiml_obstacle(): the refusal of full information maximum
# likelihood for too few cases offers listwise deletion only when listwise
# deletion would accept the data. Through the first versions of 1.1.0 the
# offer was made whenever two rows were complete: with t2_cubes constant in
# the three complete rows of five pupils, the FIML refusal named listwise
# deletion, and listwise deletion's refusal said FIML could not help.
#
# Inputs:
#   data : the item columns.
# Returns:
#   A clause, or NULL.
.relia_listwise_obstacle <- function(data) {
  data <- as.data.frame(data)
  cc <- data[stats::complete.cases(data), , drop = FALSE]
  n <- nrow(cc)
  if (n < 2L) {
    return(paste0(if (n == 0L) "no row" else "only one row", " is observed ",
                  "on every item"))
  }
  flat <- names(cc)[.relia_single_valued(cc)]
  if (length(flat)) {
    who <- if (length(flat) == 1L) {
      paste(flat, "takes")
    } else {
      paste(paste(flat[-length(flat)], collapse = ", "), "and",
            flat[length(flat)], "each take")
    }
    return(paste0(who, " a single value in the ", n, " complete rows, so ",
                  if (length(flat) == 1L) "its variance is" else
                    "their variances are", " zero there"))
  }
  if (!.relia_composite_varies(stats::cov(cc))) {
    return(paste0("the composite, the sum of the items, does not vary in the ",
                  n, " complete rows"))
  }
  NULL
}


# .stop_relia_singular_complete_rows(data, coefficient)
#
# The refusal of a coefficient fitted by the maximum likelihood
# single-factor model (robust omega and the model implied omega, the
# model implied alpha, H from data) when listwise deletion leaves no more
# rows than items. N complete rows give a covariance matrix of rank at
# most N - 1, so with N <= J it is singular, and the maximum likelihood
# discrepancy function, log|Sigma| + tr(S Sigma^-1) - log|S| - J, which
# compares the model with the saturated one whose estimate is S, is not
# defined (Joreskog, 1967; Lawley & Maxwell, 1971): lavaan refuses the
# matrix ("sample covariance matrix is not positive-definite"), or, at N
# = J, stops at a point that is not a maximum, and through the first
# versions of 1.1.0 the refusal quoted it, a report that the model did
# not converge, or now and then a coefficient. Coefficient alpha and
# KR-20 need only the sums of the matrix (its trace and the sum of its
# entries), which any two rows give, so they are computed from two or
# more complete rows; this message says so. With incomplete rows, the
# message offers full information maximum likelihood
# (.relia_incomplete_remedy()); with none, more cases.
#
# Inputs:
#   data        : the item columns before deletion.
#   coefficient : the words for the coefficient ("omega", "H", "the
#                 model implied alpha").
# Returns:
#   Nothing; it stops.
.stop_relia_singular_complete_rows <- function(data, coefficient, needs_saturated = TRUE) {
  data <- as.data.frame(data)
  J <- ncol(data)
  complete <- stats::complete.cases(data)
  n_rows <- sum(complete)
  n <- sum(!duplicated(data[complete, , drop = FALSE]))
  incomplete <- n_rows < nrow(data)
  alpha <- if (identical(coefficient, "the model implied alpha")) {
    "The analytic alpha (estimator = \"analytic\")"
  } else {
    "Coefficient alpha (reliability_alpha())"
  }
  .dmar_refuse_argument(c(
    if (incomplete) {
 paste0("Only ", n_rows, " rows of 'data' are observed on every item")
} else {
 paste0("'data' has only ", n_rows, " rows")
}, if (n < n_rows) {
 paste0(", of which ", n, " are distinct (a repeated row adds nothing ",
        "to the rank)")
}, ", no more than the ", J, " items, so the covariance matrix of ",
if (incomplete) "the complete rows" else "the items",
" is singular (its rank is at most ", n - 1L, "), and the maximum ",
"likelihood fit of the single-factor model behind ", coefficient,
" needs it positive definite: lavaan refuses such a matrix or stops at ",
"a point that is not a maximum. ",
if (incomplete) {
 .relia_incomplete_remedy(data, composite = !identical(coefficient, "H"),
                           needs_saturated = needs_saturated)
} else {
 "The model needs more cases than items."
}, " ", alpha, " needs only the sums of the covariance matrix, which ",
"two or more complete rows give."), arg = NULL)
}


# -----------------------------------------------------------------------------
# (H) Coefficient H
# -----------------------------------------------------------------------------

# .H_warn_marginal()
#
# The warning, class dmar_H_marginal_se, that an interval for H was built
# from the marginal standard errors of the loadings alone. Signaled on
# every call on the loadings path with `se_loadings` and an interval,
# because the problem belongs to the call, not to the session: the loading estimates of a
# fitted factor model are correlated, and the diagonal of their
# covariance matrix can understate or overstate the variance of H.
.H_warn_marginal <- function() {
  .dmar_warn(
    paste0(
      "The interval for H uses only the marginal standard errors in ",
      "'se_loadings' and ignores the covariances among the loading ",
      "estimates, which are correlated in a fitted factor model, so it ",
      "can be too narrow or too wide. Supply 'fit' (or 'vcov_loadings', ",
      "the covariance matrix of the standardized loadings) for the delta ",
      "method with the full covariance matrix."),
    "H_marginal_se")
}


# .pt_zero_rows(pt)
#
# Which rows of a lavaan parameter table the model holds at zero: a row
# fixed at zero (free == 0 and est == 0, the row 0*x writes), and a row
# whose label is constrained to the constant zero by an equality
# constraint (a*x with a == 0), directly or through a chain of equalities
# between labels (a == b with b == 0). Either row describes the same model
# as leaving it out, with the chi square and degrees of freedom of the
# model without it; lavaan keeps a labeled row as a free parameter, so
# the second is not caught by the first test. A constraint that is not a
# plain equality between two labels or between a label and a constant
# (a + b == 0, a == 2*b) is not followed, and a row it constrains is read
# as estimated. A constraint can name a row by its label or by the plabel
# lavaan assigns (.p4.), and both are followed. reliability_H()
# (.H_setup()), average_variance_extracted() (.ave_zero_loadings()),
# factor_scores() (.fs_items() for the items of each factor,
# .fs_auxiliary() for the auxiliary variables it sets aside, and
# .fs_zero_held() for the model matrices its scores and reasons are built
# from) read loadings and regressions through
# this, so a model reaches each of them the same way whichever of the
# three spellings (0*x, a*x with a == 0 directly or through a chain, or no
# row) it was written in.
# bifactor_indices() reads through it only the group loadings of an item
# listed under two or more group factors, where the three spellings agree.
# Elsewhere it reads a listed row as a listing: a general loading held at
# zero is a general loading of zero (a model with no such row is refused,
# since its general factor does not load on every item, and PUC, which
# presumes that it does, counts the item's correlations as general ones),
# and an item's only group loading held at zero keeps the item in that
# group factor's subscale, as its help documents (a model with no such
# row is refused, the item loading on no group factor).
#
# Inputs:
#   pt : lavaan::parTable() of a fitted model.
# Returns:
#   A logical vector over the rows of `pt`, FALSE for the rows of the
#   constraints and defined parameters themselves.
.pt_zero_rows <- function(pt) {
  constraint <- pt$op %in% c("==", "<", ">", ":=")
  fixed <- !constraint & pt$free == 0L & pt$est == 0
  label <- if (is.null(pt$label)) rep("", nrow(pt)) else pt$label
  plabel <- if (is.null(pt$plabel)) rep("", nrow(pt)) else pt$plabel
  label[is.na(label)] <- ""
  plabel[is.na(plabel)] <- ""
  ceq <- pt$op == "=="
  if (!any(ceq)) return(fixed)
  lhs <- trimws(pt$lhs[ceq])
  rhs <- trimws(pt$rhs[ceq])
  is_zero <- function(side) {
    v <- suppressWarnings(as.numeric(side))
    !is.na(v) & v == 0
  }
  held <- c(lhs[is_zero(rhs)], rhs[is_zero(lhs)], label[fixed],
            plabel[fixed])
  held <- unique(held[nzchar(held)])
  repeat {
    more <- setdiff(c(lhs[rhs %in% held], rhs[lhs %in% held]), held)
    if (!length(more)) break
    held <- c(held, more)
  }
  constrained <- !constraint &
    ((nzchar(label) & label %in% held) | (nzchar(plabel) & plabel %in% held))
  fixed | constrained
}


# .H_setup(fit, factor = NULL, group = NULL)
#
# Everything coefficient H of one factor needs from a fitted lavaan model,
# checked once. For a factor f with indicators b,
#
#     q = phi_f * lambda_b' Theta_bb^{-1} lambda_b,      H = q / (1 + q),
#
# where lambda_b are the unstandardized loadings, Theta_bb the residual
# covariance block of the indicators, and phi_f the model implied variance
# of the factor, [(I - B)^{-1} Psi (I - B)^{-T}]_ff. With a diagonal
# Theta_bb this is the H of Hancock and Mueller (2001); with residual
# covariances inside the block it is the maximal reliability of Gabler
# and Raykov (2017, Eq. 9), the reliability of the composite weighted by
# Theta_bb^{-1} lambda_b. Rescaling the factor rescales lambda_b and phi_f
# inversely, so H does not depend on how the factor was identified. For
# ordered items under the delta parameterization the residual variance of
# a latent response is not a parameter but whatever makes the latent
# response variance 1 / delta^2, so it is recomputed from the loadings
# and the factor variance at every evaluation (otherwise the gradient
# would miss that dependence).
#
# Only exported lavaan functions are used (lavInspect, lavNames, parTable),
# and only the parameter table rows that enter H are mapped: the loadings
# of the block, the residual (co)variances inside it, the latent
# (co)variances Psi, the regressions among latent variables B (including
# the loadings of a higher-order factor), and the delta scaling factors.
#
# A loading or a regression the model holds at zero (.pt_zero_rows():
# fixed at zero, 0*x, or constrained to zero, a*x with a == 0) is not
# one: the model written with that row is the model written without it,
# and it gets the same H. So an item listed under a factor with such a
# loading is not one of its indicators, and a cross-loading, a
# higher-order loading, or a regression held at zero refuses nothing;
# factor_scores() (.fs_items(), .fs_zero_held()) and
# average_variance_extracted() read such rows the same way, and
# bifactor_indices() does for the group loadings of an item listed under
# two or more group factors only (see .pt_zero_rows()).
#
# Refusals, each a stop() that names the remedy: a fit that is not a
# converged lavaan model, a multilevel fit, several groups without
# `group`, several measured factors without `factor`, a factor measured by
# other factors, an indicator that also loads on another factor or is
# regressed on another variable, fewer than two indicators, and a model
# whose assembled factor variance or indicator covariance block does not
# reproduce lavaan's own (a covariate predicting the factor, which lavaan
# carries outside the latent variable matrices, is what this catches).
#
# Returns a list: the row indices and matrix positions of each piece, the
# parameter table's `est0` and `free` columns and a `key` that identifies
# its rows, `x0` (the free parameter vector, by free number), `used_free`
# (the free parameters that enter H), and the labels `factor`, `items`,
# `g` (group index), `group` (label, or NA for a single-group fit),
# `ordered` (the block's ordered items), `delta_param`.
.H_setup <- function(fit, factor = NULL, group = NULL) {
  if (!inherits(fit, "lavaan")) {
    .dmar_refuse_argument(c(
      "'fit' must be a fitted lavaan model (from lavaan::cfa(), ",
      "lavaan::sem(), or fa_confirmatory(output = \"fit\"))."), arg = "fit")
  }
  if (!isTRUE(lavaan::lavInspect(fit, "converged"))) {
    .dmar_refuse(c(
      "The fitted model did not converge; H is not computed from a ",
      "nonconverged solution."), "not_converged")
  }
  if (lavaan::lavInspect(fit, "nlevels") > 1L) {
    .dmar_refuse_argument(c(
      "H is computed here for single-level models; this fit has ",
      lavaan::lavInspect(fit, "nlevels"), " levels."), arg = NULL)
  }

  ng <- lavaan::lavInspect(fit, "ngroups")
  labels <- lavaan::lavInspect(fit, "group.label")
  if (ng > 1L) {
    if (is.null(group)) {
      .dmar_refuse_argument(c(
        "The fit has ", ng, " groups (", paste(labels, collapse = ", "),
        "); name one with 'group'. H is computed for one group at a ",
        "time, from that group's parameters."), arg = "group")
    }
    g <- if (length(group) == 1L) match(as.character(group), labels) else NA
    if (is.na(g)) {
      .dmar_refuse_argument(c(
        "'group' must be one of the fit's group labels: ",
        paste(labels, collapse = ", "), "."), arg = "group")
    }
    group_label <- labels[g]
  } else {
    if (!is.null(group)) {
      .dmar_refuse_argument(c(
        "'group' applies to a fit with several groups; this fit has ",
        "one."), arg = "group")
    }
    g <- 1L
    group_label <- NA_character_
  }

  pt <- lavaan::parTable(fit)
  if (is.null(pt$group)) pt$group <- rep(1L, nrow(pt))
  lv <- lavaan::lavNames(fit, "lv")
  opts <- lavaan::lavInspect(fit, "options")
  in_g <- pt$group == g
  # Rows fixed or constrained at zero, which describe no loading or
  # regression (see above).
  zero <- .pt_zero_rows(pt)

  # Rows in which a factor is measured by an observed variable.
  ind <- in_g & pt$op == "=~" & !(pt$rhs %in% lv) & !zero
  measured <- unique(pt$lhs[ind])
  if (is.null(factor)) {
    if (length(measured) != 1L) {
      .dmar_refuse_argument(c(
        "The fit has ", length(measured), " factors measured by ",
        "observed items (", paste(measured, collapse = ", "),
        "); name one with 'factor'."), arg = "factor")
    }
    factor <- measured
  } else if (!is.character(factor) || length(factor) != 1L ||
             !factor %in% lv) {
    .dmar_refuse_argument(c(
      "'factor' must name one latent variable of the fit: ",
      paste(lv, collapse = ", "), "."), arg = "factor")
  }
  if (any(in_g & pt$op == "=~" & pt$lhs == factor & pt$rhs %in% lv &
          !zero)) {
    .dmar_refuse_argument(c(
      "'", factor, "' is measured by other factors (a higher-order ",
      "factor). H is defined here for a factor measured by observed ",
      "items; compute it for each first-order factor."), arg = NULL)
  }
  lam_rows <- which(ind & pt$lhs == factor)
  items <- pt$rhs[lam_rows]
  if (length(items) < 2L) {
    .dmar_refuse_argument(c(
      "H needs two or more indicators of '", factor, "'; the fit has ",
      length(items), "."), arg = NULL)
  }
  cross <- unique(items[items %in% pt$rhs[ind & pt$lhs != factor]])
  if (length(cross)) {
    .dmar_refuse_argument(c(
      "Item(s) ", paste(cross, collapse = ", "), " also load on ",
      "another factor. H is defined here for a block of items that load ",
      "on '", factor, "' alone (a loading fixed or constrained to zero ",
      "does not count, so a cross-loading written 0*item, or a*item ",
      "with a == 0, is not one). For a bifactor ",
      "model, bifactor_indices() ",
      "reports the index Rodriguez, Reise, and Haviland (2016) call H, ",
      "which counts the other factors' common variance as uncorrelated ",
      "error and is not this maximal reliability, and its ",
      "\"maximal_reliability\" attribute gives the maximal reliability ",
      "the fitted model implies for each factor; ci_bifactor_indices() ",
      "gives both an interval."), arg = NULL)
  }
  regressed <- unique(items[items %in% pt$lhs[in_g & pt$op == "~" & !zero]])
  if (length(regressed)) {
    .dmar_refuse_argument(c(
      "Item(s) ", paste(regressed, collapse = ", "), " are regressed ",
      "on other variables, so their variance is not the factor's and ",
      "their residuals' alone. Compute H from the measurement model of ",
      "'", factor, "'."), arg = NULL)
  }

  th_rows  <- which(in_g & pt$op == "~~" & pt$lhs %in% items &
                    pt$rhs %in% items)
  psi_rows <- which(in_g & pt$op == "~~" & pt$lhs %in% lv & pt$rhs %in% lv)
  # A higher-order loading (f2 =~ f1) is the regression of f1 on f2.
  bho_rows <- which(in_g & pt$op == "=~" & pt$lhs %in% lv & pt$rhs %in% lv)
  breg_rows <- which(in_g & pt$op == "~" & pt$lhs %in% lv & pt$rhs %in% lv)
  dl_rows  <- which(in_g & pt$op == "~*~" & pt$lhs %in% items)
  ordered_items <- intersect(items, lavaan::lavNames(fit, "ov.ord"))

  s <- list(
    lam_rows = lam_rows,
    th_rows  = th_rows,
    th_idx   = cbind(match(pt$lhs[th_rows], items),
                     match(pt$rhs[th_rows], items)),
    psi_rows = psi_rows,
    psi_idx  = cbind(match(pt$lhs[psi_rows], lv),
                     match(pt$rhs[psi_rows], lv)),
    b_rows   = c(bho_rows, breg_rows),
    b_idx    = rbind(cbind(match(pt$rhs[bho_rows], lv),
                           match(pt$lhs[bho_rows], lv)),
                     cbind(match(pt$lhs[breg_rows], lv),
                           match(pt$rhs[breg_rows], lv))),
    dl_rows  = dl_rows,
    dl_pos   = match(pt$lhs[dl_rows], items),
    ord_pos  = match(ordered_items, items),
    f = match(factor, lv), nlv = length(lv), J = length(items),
    delta_param = length(ordered_items) > 0L &&
      identical(opts$parameterization, "delta"),
    est0 = pt$est,
    free = pt$free,
    key  = paste(pt$lhs, pt$op, pt$rhs, pt$group),
    factor = factor, items = items, g = g, group = group_label,
    ordered = ordered_items)
  npar <- max(c(0L, pt$free))
  s$x0 <- pt$est[match(seq_len(npar), pt$free)]
  used <- c(lam_rows, th_rows, psi_rows, s$b_rows, dl_rows)
  s$used_free <- sort(unique(pt$free[used][pt$free[used] > 0L]))

  # The assembled matrices must reproduce lavaan's own factor variance and
  # indicator covariance block; anything the representation above leaves
  # out (an observed covariate predicting the factor, which lavaan carries
  # outside the latent variable matrices) shows up as a difference.
  e <- .H_eval(s, s$est0)
  cov_lv <- lavaan::lavInspect(fit, "cov.lv")
  if (ng > 1L) cov_lv <- cov_lv[[g]]
  phi_lav <- cov_lv[factor, factor]
  sig_lav <- lavaan::lavInspect(fit, "cov.ov")
  if (ng > 1L) sig_lav <- sig_lav[[g]]
  sig_lav <- sig_lav[items, items, drop = FALSE]
  # With ordered items lavaan reports the latent response covariances
  # scaled by the delta factors, Delta Sigma Delta (unit variances), in
  # either parameterization; continuous items have delta = 1.
  sig <- (e$phi * tcrossprod(e$lambda) + e$theta) * tcrossprod(e$delta)
  if (!isTRUE(abs(e$phi - phi_lav) <= 1e-6 * max(1, abs(phi_lav))) ||
      !isTRUE(max(abs(sig - sig_lav)) <= 1e-6 * max(1, abs(sig_lav)))) {
    .dmar_refuse_argument(c(
      "The variance of '", factor, "' or the covariances of its items ",
      "depend on parts of the model that H does not represent (for ",
      "example an observed covariate that predicts the factor). Fit ",
      "the measurement model of the factor's items on its own and ",
      "compute H from that fit."), arg = NULL)
  }
  s
}


# .H_eval(s, est)
#
# H at the parameter values `est` (a full parameter table `est` column),
# with the pieces it was built from (`delta` is 1 for an item without a
# scaling factor). `proper` is FALSE when the factor
# variance is not positive or the residual covariance block of the
# indicators is not positive definite (an improper, Heywood, solution),
# for which H is not defined; H is then NA.
.H_eval <- function(s, est) {
  lam <- est[s$lam_rows]
  psi <- matrix(0, s$nlv, s$nlv)
  psi[s$psi_idx] <- est[s$psi_rows]
  psi[s$psi_idx[, 2:1, drop = FALSE]] <- est[s$psi_rows]
  bm <- matrix(0, s$nlv, s$nlv)
  if (length(s$b_rows)) bm[s$b_idx] <- est[s$b_rows]
  ib_f <- solve(diag(s$nlv) - bm)[s$f, ]
  phi <- drop(ib_f %*% psi %*% ib_f)
  th <- matrix(0, s$J, s$J)
  th[s$th_idx] <- est[s$th_rows]
  th[s$th_idx[, 2:1, drop = FALSE]] <- est[s$th_rows]
  dlt <- rep(1, s$J)
  dlt[s$dl_pos] <- est[s$dl_rows]
  if (s$delta_param) {
    for (j in s$ord_pos) th[j, j] <- 1 / dlt[j]^2 - lam[j]^2 * phi
  }
  proper <- is.finite(phi) && phi > 0 && all(is.finite(th)) &&
    all(diag(th) > 0) &&
    min(eigen(th, symmetric = TRUE, only.values = TRUE)$values) > 0
  H <- if (proper) {
    q <- phi * sum(lam * solve(th, lam))
    q / (1 + q)
  } else {
    NA_real_
  }
  list(H = H, proper = proper, phi = phi, lambda = lam, theta = th,
       delta = dlt)
}


# .H_est_at(s, x)
#
# The parameter table's `est` column with the free parameters replaced by
# the vector `x` (indexed by free number), for the numerical gradient.
.H_est_at <- function(s, x) {
  est <- s$est0
  fr <- s$free > 0L
  est[fr] <- x[s$free[fr]]
  est
}


# .H_from_fit(fit, factor = NULL, group = NULL, se = TRUE)
#
# Coefficient H of one factor of a fitted lavaan model, with its delta
# method variance over the full covariance matrix of the fit's free
# parameters:
#
#     Var(H-hat) = g' V g,    g = dH/dx at the estimate,   V = vcov(fit),
#
# so the variance carries whatever standard error the fit was estimated
# with (standard, robust.sem, robust.huber.white under FIML, the robust
# WLSMV matrix). The gradient is a central difference of the closed form
# with step 1e-6 * max(1, |x_k|), taken over the free parameters that
# enter H only (`used_free`); every other element of g is zero. An
# improper solution stops with an error of class dmar_H_improper.
#
# The gradient runs over free numbers, and parameters that share a label
# (or a group.equal constraint) share one free number when the fit has
# ceq.simple = TRUE. lavaan's vcov() is then still laid out over the
# unconstrained parameters, one row for every parameter table row with
# free > 0, in table order: K V_r K', where V_r is the covariance matrix
# of the free numbers and K copies each of them into every row that
# carries it. The rows of V that belong to free number k are therefore
# found by matching k among the free numbers of those rows; every row
# that carries k gives the same V_r entries, and the first is taken.
# Under ceq.simple = FALSE (lavaan's default) each row has a free number
# of its own and the match is the identity. Indexing V by free number
# under ceq.simple = TRUE would read the wrong rows for every parameter
# after the first shared label in the table (any group but the first of
# a group.equal fit, or a block with labelled equal loadings) and give a
# standard error that is too small or too large.
#
# Returns a list: H, var (NA when se = FALSE, or when the covariance
# matrix of the estimates cannot be used, with a warning), the setup
# `s`, the standardized loadings of the block (named), the number of residual
# covariances inside the block, `metric` ("observed" or
# "latent_response"), `se_type`, `missing` ("listwise", "fiml",
# "pairwise", or lavaan's own label), `N` (that group's cases), and
# `N_complete` (cases complete on the block's items; N for a fit to
# summary statistics).
.H_from_fit <- function(fit, factor = NULL, group = NULL, se = TRUE) {
  s <- .H_setup(fit, factor = factor, group = group)
  e <- .H_eval(s, s$est0)
  if (!e$proper) {
    .dmar_refuse(
      paste0(
        "The solution is improper (a Heywood case): a residual variance ",
        "of an item of '", s$factor, "' is not positive, or the factor ",
        "variance is not, so H is not defined for it. Reconsider the ",
        "model for these items before reporting their reliability."),
      "H_improper")
  }
  opts <- lavaan::lavInspect(fit, "options")

  var_H <- NA_real_
  # The covariance matrix of the estimates is used only when the factor
  # analysis tables would use it (.cfa_k_vcov_status()): when lavaan
  # could not compute it, or computed it from an information matrix that
  # is singular in the directions the fit's equality constraints leave
  # free, a parameter is not identified at the estimates and g' V g is
  # meaningless (an H of .853 with a standard error of 534 and the
  # interval [0, 1] on the fixture of the test of this guard, where
  # fa_confirmatory() withheld every standard error of the same fit).
  # The estimate stays; the variance is NA, with the tables' warning.
  # Every information is tested, observed or expected, as the tables test
  # it, since a fit passed in may have either and lavaan can invert either
  # when it is singular. Through 1.1.0 development the variance was formed
  # from whatever lavaan returned, and a matrix lavaan could not compute
  # stopped the call with an R error.
  vcov_status <- if (se) .cfa_k_vcov_status(fit)
  if (se && vcov_status %in% c("singular", "failed")) {
    .cfa_k_vcov_warn(vcov_status, failed_rest = c(
      "the standard error and interval of H are NA."))
    se <- FALSE
  }
  if (se) {
    V <- lavaan::vcov(fit)
    k <- s$used_free
    grad <- vapply(k, function(i) {
      h <- 1e-6 * max(1, abs(s$x0[i]))
      xp <- s$x0
      xm <- s$x0
      xp[i] <- xp[i] + h
      xm[i] <- xm[i] - h
      (.H_eval(s, .H_est_at(s, xp))$H -
         .H_eval(s, .H_est_at(s, xm))$H) / (2 * h)
    }, numeric(1))
    row_free <- s$free[s$free > 0L]
    v_rows <- if (nrow(V) == length(row_free)) {
      match(k, row_free)
    } else if (nrow(V) == length(s$x0)) {
      k
    } else {
      .dmar_refuse_argument(c(
        "The covariance matrix of the estimates has ", nrow(V),
        " rows, matching neither the ", length(s$x0), " free ",
        "parameters of the fit nor its ", length(row_free),
        " unconstrained ones, so the delta method standard error of ",
        "H cannot be formed. Use ci_method = \"none\" or a bootstrap ",
        "method."), arg = NULL)
    }
    var_H <- drop(crossprod(grad, V[v_rows, v_rows, drop = FALSE] %*% grad))
  }

  sigma_jj <- e$lambda^2 * e$phi + diag(e$theta)
  loadings_std <- stats::setNames(e$lambda * sqrt(e$phi) / sqrt(sigma_jj),
                                  s$items)
  n_rescov <- sum(e$theta[upper.tri(e$theta)] != 0)

  nobs <- lavaan::lavInspect(fit, "nobs")
  N <- as.numeric(nobs[s$g])
  X <- tryCatch(lavaan::lavInspect(fit, "data"), error = function(err) NULL)
  if (is.list(X) && !is.data.frame(X)) X <- X[[s$g]]
  N_complete <- if (is.null(X) || !length(X) || is.null(colnames(X)) ||
                    !all(s$items %in% colnames(X))) {
    N
  } else {
    sum(stats::complete.cases(X[, s$items, drop = FALSE]))
  }
  missing <- if (opts$missing %in% c("ml", "ml.x", "fiml", "direct")) {
    "fiml"
  } else {
    opts$missing
  }

  list(H = e$H, var = var_H, setup = s,
       loadings = loadings_std, residual_covariances = n_rescov,
       metric = if (length(s$ordered)) "latent_response" else "observed",
       se_type = opts$se, estimator = opts$estimator, missing = missing,
       N = N, N_complete = N_complete)
}


# .H_replicate(s, fit)
#
# H from a refit of the same model (a bootstrap replication), read through
# the setup `s` of the original fit when the refit's parameter table has
# the same rows, which it has whenever the same model was refitted; a
# table that differs gets a setup of its own. Returns H, or NA_real_ with
# the attribute dmar_reason = "improper" for an improper solution, which
# .bootstrap_ci() counts by that cause.
.H_replicate <- function(s, fit) {
  pt <- lavaan::parTable(fit)
  if (is.null(pt$group)) pt$group <- rep(1L, nrow(pt))
  if (!identical(paste(pt$lhs, pt$op, pt$rhs, pt$group), s$key)) {
    s <- .H_setup(fit, factor = s$factor)
    pt$est <- s$est0
  }
  e <- .H_eval(s, pt$est)
  if (!e$proper) return(structure(NA_real_, dmar_reason = "improper"))
  e$H
}


# .H_fit_data(fit)
#
# The rows the bootstrap of the fit path resamples: the fit's own data
# (lavInspect(fit, "data")), as a data frame named by the observed
# variables, without the rows that have no observed value (lavaan ignores
# them under FIML and would warn about them on every replication). Stops
# when the fit carries no raw data (a fit to summary statistics).
.H_fit_data <- function(fit) {
  X <- tryCatch(lavaan::lavInspect(fit, "data"), error = function(err) NULL)
  ov <- lavaan::lavNames(fit, "ov")
  if (is.null(X) || !length(X) || is.null(colnames(X)) ||
      !all(ov %in% colnames(X))) {
    .dmar_refuse_argument(c(
      "A bootstrap needs the raw data the model was fitted to, and this ",
      "fit does not carry them (for example, it was fitted to a ",
      "covariance matrix). Refit with 'data', or use ci_method = ",
      "\"delta_logistic\"."), arg = "data")
  }
  X <- as.data.frame(X[, ov, drop = FALSE])
  X[rowSums(!is.na(X)) > 0L, , drop = FALSE]
}


# .H_refit_fun(fit, s)
#
# The point function of the fit path bootstrap: given resampled rows d, it
# refits the user's model through lavaan's public interface,
# lavaan::lavaan(model = pt, data = d, slotOptions = opts), with the
# fitted parameter table as the model (its estimates as starting values)
# and the fit's own options with standard errors, test statistics, the
# baseline model, and the unrestricted model switched off (none of them
# enters H), and returns H. reliability_H() refits the original rows
# with it before any resample is drawn and stops unless that refit
# reproduces the estimate; it refuses a fit with sampling weights or a
# cluster variable, which a refit of resampled rows would not carry. A
# resample that leaves a response category of an ordered item
# empty cannot estimate that item's thresholds, so it stops with an error
# that says so (the bootstrap engine counts it and quotes the message). A
# resample in which a continuous variable of the model takes a single
# value is refused the same way, by name (.relia_check_resample()):
# lavaan would stop on the variable after printing its table of the
# variables to the console, once for every such replication. So is a
# resample with too few distinct rows for the fit's own treatment of
# incomplete rows (`rule`): a listwise fit, or one to complete data, needs
# more distinct rows than continuous variables, since the sample covariance
# matrix of the resample is otherwise singular, and a FIML or two-stage
# fit a conservative distinct-case screen, separately from the H1
# likelihood diagnostic, which is also applied to each two-stage refit.
# Historical count-policy background: the rule of the
# saturated model (see .relia_check_resample()). A refit whose likelihood
# has no maximum (an error variance running off below zero, the
# likelihood flat or still rising; .lavaan_no_maximum()) stops with
# .no_max_replication()'s short refusal,
# which the engine counts among the replications that raised an error.
.H_refit_fun <- function(fit, s) {
  # The fit's estimates are kept in the est column, which is where lavaan
  # reads the starting values of a parameter table given as the model: it
  # starts from est when the table has one (unless an estimate is missing
  # or a free variance is exactly zero) and computes its own default
  # starting values otherwise, overwriting the start column either way.
  # Through the first versions of 1.1.0 the estimates were moved to start
  # and est was dropped, so every refit iterated from lavaan's defaults.
  pt <- lavaan::parTable(fit)
  pt$se <- NULL
  opts <- lavaan::lavInspect(fit, "options")
  fixed_predictors <- if (isTRUE(opts$fixed.x) || isTRUE(opts$conditional.x)) {
    lavaan::lavNames(fit, "ov.x")
  } else character()
  rule <- switch(opts$missing,
                 ml = , ml.x = "fiml",
                 two.stage = , robust.two.stage = "two_stage",
                 opts$missing)
  opts$se <- "none"
  opts$test <- "none"
  opts$baseline <- FALSE
  opts$h1 <- FALSE
  ord <- lavaan::lavNames(fit, "ov.ord")
  X <- .H_fit_data(fit)
  continuous <- setdiff(colnames(X), ord)
  vars <- if (length(ord)) "continuous variables" else "variables"
  n_cat <- vapply(ord, function(v) length(unique(stats::na.omit(X[[v]]))),
                  integer(1))
  function(d) {
    if (length(ord)) {
      kept <- vapply(ord, function(v) length(unique(stats::na.omit(d[[v]]))),
                     integer(1))
      if (any(kept < n_cat)) {
        .dmar_refuse_argument(c(
          "the resample left a response category of an ordered item ",
          "empty"), arg = NULL)
      }
    }
    if (length(continuous)) {
      .relia_check_resample(d, missing = rule, items = continuous,
                            fitted = TRUE, on_fit = vars)
    }
    if (rule == "two_stage") {
      response_data <- d[, !names(d) %in% fixed_predictors, drop = FALSE]
      .stop_never_observed_together(response_data, "Two-stage coefficient H refit")
      .stop_single_valued_together(response_data, "Two-stage coefficient H refit")
      .stop_fiml_unbounded(response_data, "Two-stage coefficient H refit")
    }
    # lavaan's post-estimation checks warn about a nonconverged or improper
    # solution on a resample; both are counted below instead (a
    # nonconverged refit returns NA, an improper one NA with its reason),
    # so the warnings, one per replication, are noise here.
    refit <- suppressWarnings(
      lavaan::lavaan(model = pt, data = d, slotOptions = opts))
    if (!isTRUE(lavaan::lavInspect(refit, "converged"))) return(NA_real_)
    # A refit whose likelihood has no maximum stops, and the engine counts
    # it among the replications that raised an error
    # (.no_max_replication()).
    no_max <- .lavaan_no_maximum(refit)
    if (!is.null(no_max)) .no_max_replication(no_max)
    .H_replicate(s, refit)
  }
}


# .H_from_loadings(loadings, V = NULL)
#
# H from standardized loadings, the closed form of Hancock and Mueller
# (2001) with uncorrelated residuals (Theta = diag(1 - l^2)):
#
#     s = sum(l^2 / (1 - l^2)),   H = s / (1 + s),
#     dH/dl_i = 2 l_i / ((1 - l_i^2)^2 (1 + s)^2),
#
# with the delta method variance g' V g when V, the covariance matrix of
# the standardized loading estimates (full, or diagonal from marginal
# standard errors), is supplied. Returns list(H, var).
.H_from_loadings <- function(loadings, V = NULL) {
  s <- sum(loadings^2 / (1 - loadings^2))
  H <- s / (1 + s)
  var_H <- if (is.null(V)) {
    NA_real_
  } else {
    g <- 2 * loadings / ((1 - loadings^2)^2 * (1 + s)^2)
    drop(crossprod(g, V %*% g))
  }
  list(H = H, var = var_H)
}


# .H_improper_note()
#
# The sentences the bootstrap engine's warning gives, in place of its
# own, about the replications of coefficient H set aside as improper. H
# is not defined for an improper solution, so reliability_H() passes
# improper = "drop" to .bootstrap_ci(), which records that rule. An item
# with a residual variance of exactly zero (a standardized loading of
# one) leaves H undefined too, which the engine's sentence, about an
# error variance below zero, does not say. The three-item resample whose
# correlations admit no single-factor solution, which the data path
# counts in `no_solution` as reliability_omega() does, is described as
# the engine describes it: H, like omega, runs high as a resample nears
# that boundary, while the correlations of one beyond it point to low
# reliability, so it belongs on neither side of the distribution.
.H_improper_note <- function() {
  paste0(
    "An improper replication is one whose solution gives an item a ",
    "residual variance that is not positive (a standardized loading at or ",
    "past one), or, with three items, one whose correlations admit no ",
    "single-factor solution. H is not defined for either, so such ",
    "replications are always set aside. Setting aside the first kind can ",
    "trim the upper tail of the bootstrap distribution; the second ",
    "belongs on neither side of it, since its correlations point to low ",
    "reliability while the estimator runs high as a resample nears that ",
    "boundary")
}


# .H_footer_lines(x)
#
# The lines a table of coefficient H prints under itself, beside the
# family's own (the confidence level, the bootstrap counts): which factor
# (and group) of a fitted model H describes, whether residual covariances
# inside the block entered the weights, that H on ordered items describes
# the composite of the latent responses, that a standard error built from
# marginal standard errors ignored their covariances, what an interval
# from standardized loadings needs when they came with nothing it could
# be formed from, and the bootstrap estimate of bias. .dmar_tbl_footer() in
# R/dmar_output_helpers.R calls it, so the console and a knitted report
# show the same lines. Returns a named character vector, empty for any
# table that is not coefficient H.
.H_footer_lines <- function(x) {
  lines <- character(0)
  if (!identical(attr(x, "coefficient", exact = TRUE), "H")) return(lines)
  parts <- character(0)
  fac <- attr(x, "factor", exact = TRUE)
  if (is.character(fac) && length(fac) == 1L && !is.na(fac)) {
    grp <- attr(x, "group", exact = TRUE)
    parts <- c(parts, paste0(
      "Factor: ", fac,
      if (is.character(grp) && length(grp) == 1L && !is.na(grp)) {
        paste0(" (group ", grp, ")")
      } else {
        ""
      }))
  }
  k <- attr(x, "residual_covariances", exact = TRUE)
  if (is.numeric(k) && length(k) == 1L && isTRUE(k > 0)) {
    parts <- c(parts, sprintf(paste0(
      "%d residual %s inside the block %s the optimal weights ",
      "(Gabler & Raykov, 2017)"),
      as.integer(k), if (k == 1) "covariance" else "covariances",
      if (k == 1) "enters" else "enter"))
  }
  if (identical(attr(x, "metric", exact = TRUE), "latent_response")) {
    parts <- c(parts, paste0(
      "H describes the optimally weighted composite of the latent ",
      "responses behind the ordered categories, not of the categories"))
  }
  if (identical(attr(x, "vcov_basis", exact = TRUE), "diagonal")) {
    parts <- c(parts, paste0(
      "the standard error uses the marginal standard errors of the ",
      "loadings and ignores their covariances"))
  }
  # Standardized loadings alone carry nothing an interval can be formed
  # from, so a call with them and no covariance matrix returns the
  # estimate alone, and the table says what an interval would need. A
  # call that supplied vcov_loadings or se_loadings and asked for
  # ci_method = "none" has what an interval needs, and gets no such line.
  if (identical(attr(x, "metric", exact = TRUE), "unknown") &&
      identical(attr(x, "ci_method", exact = TRUE), "none") &&
      identical(attr(x, "vcov_supplied", exact = TRUE), "none")) {
    parts <- c(parts, paste0(
      "no interval (ci_method = \"none\"): an interval from standardized ",
      "loadings needs 'vcov_loadings', the covariance matrix of their ",
      "estimates"))
  }
  if (length(parts)) {
    txt <- paste(parts, collapse = "; ")
    lines["H"] <- paste0(toupper(substr(txt, 1L, 1L)), substring(txt, 2L))
  }
  bs <- attr(x, "bootstrap", exact = TRUE)
  if (is.list(bs) && is.numeric(bs$bias) && length(bs$bias) == 1L &&
      is.finite(bs$bias)) {
    lines["H_bias"] <- sprintf(paste0(
      "Bootstrap estimate of the bias of H (mean of the replicates minus ",
      "the estimate): %s"),
      formatC(bs$bias, digits = 3, format = "g"))
  }
  lines
}


# .H_finish(h, ci_method, boot_out, conf_level, B, aux, factor, group)
#
# The reliability family table for coefficient H from a fit-based result
# `h` (.H_from_fit(), with N, N_complete, and missing possibly set by the
# caller) and, for a bootstrap, the list .bootstrap_ci() returned. The
# delta method interval is the Wald interval on the coefficient scale, or
# on the logit scale (Browne, 1982) for the _logistic methods; a
# bootstrap's interval is the engine's, and its `bootstrap` attribute
# gains `bias`, the mean of the returned replicates minus the estimate.
.H_finish <- function(h, ci_method, boot_out, conf_level, B, aux,
                      factor, group, h1_diagnostic = NULL) {
  bootstrap <- NULL
  if (!is.null(boot_out)) {
    se <- boot_out$se
    lim <- list(lower = boot_out$lower, upper = boot_out$upper)
    method <- boot_out$ci_method
    bootstrap <- boot_out$bootstrap
    bootstrap$bias <- boot_out$replicate_mean - h$H
    basis <- "none"
  } else if (ci_method == "none") {
    se <- NA_real_
    lim <- list(lower = NA_real_, upper = NA_real_)
    method <- ci_method
    basis <- "none"
  } else {
    se <- sqrt(h$var)
    lim <- .H_wald(h$H, se, conf_level, endsWith(ci_method, "_logistic"))
    method <- ci_method
    # A variance .H_from_fit() withheld (a covariance matrix of the
    # estimates that cannot be used) formed no standard error.
    basis <- if (is.na(se)) "none" else "full"
  }
  out <- .relia_result(
    estimate = h$H, se = se, lower = lim$lower, upper = lim$upper,
    conf_level = conf_level, N = h$N, J = length(h$loadings),
    coefficient = "H", ci_method = method,
    B = if (!is.null(boot_out)) B else NA_integer_,
    N_complete = h$N_complete, missing = h$missing, aux = aux,
    bootstrap = bootstrap)
  out <- .H_attributes(out, vcov_basis = basis, se_type = h$se_type,
                factor = factor, group = group, metric = h$metric,
                residual_covariances = h$residual_covariances,
                loadings = h$loadings)
  .relia_h1_qualify(out, h1_diagnostic)
}


# .H_wald(estimate, se, conf_level, logistic)
#
# The Wald interval for H, untransformed or on the logit scale through
# .logistic_ci(); NA limits when the standard error is NA.
.H_wald <- function(estimate, se, conf_level, logistic) {
  if (is.na(se)) return(list(lower = NA_real_, upper = NA_real_))
  crit <- stats::qnorm(1 - (1 - conf_level) / 2)
  if (logistic) {
    .logistic_ci(estimate, se, crit)
  } else {
    list(lower = estimate - crit * se, upper = estimate + crit * se)
  }
}


# .H_attributes(out, ...)
#
# Attach the attributes that record how coefficient H was computed (see
# the Value section of ?reliability_H) to a family table. `vcov_supplied`
# ("full", "diagonal", or "none") is what came with standardized
# loadings, set on that path only.
.H_attributes <- function(out, vcov_basis, se_type, factor, group, metric,
                          residual_covariances, loadings,
                          vcov_supplied = NULL) {
  attr(out, "vcov_basis") <- vcov_basis
  if (!is.null(vcov_supplied)) attr(out, "vcov_supplied") <- vcov_supplied
  attr(out, "se_type") <- se_type
  attr(out, "factor") <- factor
  attr(out, "group") <- group
  attr(out, "metric") <- metric
  attr(out, "residual_covariances") <- as.integer(residual_covariances)
  attr(out, "loadings") <- loadings
  out
}
