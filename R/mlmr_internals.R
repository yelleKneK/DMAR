#' Internal Helpers for Mlmr
#'
#' Not exported. The helpers translate the formula and data into a
#' lavaan model and data frame, construct the lavaan syntax, and
#' compute confidence intervals.
#'
#' @keywords internal
#' @name mlmr-internals
NULL


# Guard the lavaan dependency for mlmr() and mlmr_mv(). lavaan is a
# Suggests, so confirm it is installed, then refuse an old lavaan in an
# environment where the CPU count is unknown. lavaan builds its default
# 'ncpus' option from max(1L, parallel::detectCores() - 1L); when
# detectCores() returns NA (some containers, restricted schedulers, and
# CI runners report this) that default becomes NA, and lavaan 0.6-19
# then fails deep in option validation with the opaque message "missing
# value where TRUE/FALSE needed". lavaan 0.7-2 (published 2026-07-16) is
# the first CRAN version reported to resolve this, so the floor is
# enforced only when the CPU count cannot be detected, where an old
# lavaan cannot succeed anyway. The floor and core count are arguments
# so the decision is testable without a live NA-core environment. The
# version comparison goes through as.numeric_version() so "0.6-19" and
# "0.10-1" order numerically rather than lexicographically.
.mlmr_require_lavaan <- function(fn,
                                 n_cores = parallel::detectCores(),
                                 lavaan_version =
                                   utils::packageVersion("lavaan"),
                                 lavaan_floor = "0.7-2") {
  if (!requireNamespace("lavaan", quietly = TRUE)) {
    .dmar_refuse(c(
      "The 'lavaan' package is required for ", fn, "(). ",
      "Install it with install.packages(\"lavaan\")."),
      "package_required", package = "lavaan")
  }
  if (isTRUE(is.na(n_cores)) &&
      as.numeric_version(lavaan_version) < as.numeric_version(lavaan_floor)) {
    .dmar_refuse(c(
      "This session cannot detect the number of CPU cores ",
      "(parallel::detectCores() is NA), and the installed 'lavaan' (",
      as.character(lavaan_version), ") derives its default 'ncpus' ",
      "option from that count. With an undetectable core count that ",
      "option is NA and ", fn, "() fails inside lavaan with the opaque ",
      "message \"missing value where TRUE/FALSE needed\". Upgrade ",
      "lavaan to ", lavaan_floor, " or newer with ",
      "install.packages(\"lavaan\"), or run where ",
      "parallel::detectCores() returns a value."),
      "package_required", package = "lavaan")
  }
  invisible(TRUE)
}


# Expand the user formula into a lavaan-friendly model frame.
#
# Returns a list with:
#   y_name        the internal name for the response ("y")
#   X_names       internal lavaan-safe predictor names ("x_1", ...)
#   X_display     display names from model.matrix() (e.g., "wt",
#                 "groupB", "wt:hp"); same length and order as X_names
#   has_intercept whether the formula has an intercept
#   model_data    data.frame the lavaan fit consumes, columns
#                 c(y_name, X_names), missing values preserved
#   model_frame   the original model.frame (na.action = na.pass)
#   terms         the terms object
#
# The renaming step exists because model.matrix can produce column
# names with characters (colons, parentheses) that are not legal in
# lavaan's syntax. We keep display names separately so the output
# tables read naturally.
mlmr_build_design <- function(formula, data) {
  trms <- stats::terms(formula, data = data)
  mf <- stats::model.frame(trms, data = data,
                           na.action = stats::na.pass)
  has_intercept <- as.logical(attr(trms, "intercept"))

  y_vec <- stats::model.response(mf)
  if (is.null(y_vec)) {
    .dmar_refuse_argument(
      "'formula' must have a response on the left-hand side.", arg = "formula")
  }
  if (!is.numeric(y_vec)) {
    .dmar_refuse_argument(c(
      "The response variable must be numeric. Got class: ",
      paste(class(y_vec), collapse = ", "), "."), arg = NULL)
  }

  X_full <- stats::model.matrix(trms, mf)
  X <- if (has_intercept) X_full[, -1L, drop = FALSE] else X_full
  X_display <- colnames(X)
  X_names <- if (ncol(X) > 0L) paste0("x_", seq_len(ncol(X)))
             else character(0)

  y_name <- "y"
  model_data <- data.frame(y = as.numeric(y_vec))
  names(model_data) <- y_name
  if (length(X_names) > 0L) {
    X_df <- as.data.frame(X)
    names(X_df) <- X_names
    model_data <- cbind(model_data, X_df)
  }

  list(
    y_name = y_name,
    X_names = X_names,
    X_display = X_display,
    has_intercept = has_intercept,
    model_data = model_data,
    model_frame = mf,
    terms = trms
  )
}


# Build the lavaan model syntax string.
#
# The regression line uses labeled coefficients ("b0" for the
# intercept and "b_1", "b_2", ... for the slopes) so the labels can
# be referenced later when imposing equality constraints during the
# profile likelihood CI computation. The residual variance carries a
# label as well ("sigma2_e") for downstream extraction.
mlmr_build_syntax <- function(y_name, X_names, has_intercept) {
  slope_labels <- paste0("b_", seq_along(X_names))

  if (length(X_names) > 0L) {
    slope_terms <- paste0(slope_labels, "*", X_names, collapse = " + ")
  } else {
    slope_terms <- ""
  }

  intercept_piece <- if (has_intercept) "b0*1" else "0*1"
  rhs <- if (nzchar(slope_terms)) paste(intercept_piece, slope_terms,
                                         sep = " + ")
         else intercept_piece
  reg_line <- paste0(y_name, " ~ ", rhs)

  resid_line <- paste0(y_name, " ~~ sigma2_e*", y_name)

  paste(reg_line, resid_line, sep = "\n")
}


# Validate the 'auxiliary' argument and build the renamed auxiliary
# data frame. Shared by mlmr() and mlmr_mv().
#
# auxiliary   character vector of variable names in 'data', or NULL
# data        the user's data.frame
# model_vars  the variables already used by the model (all.vars of the
#             formula: response(s) and predictors). An auxiliary may
#             not double as a modeled variable.
# fixed_x     the fixed_x setting; auxiliaries require fixed_x = FALSE
#             because the saturated correlates approach models the
#             joint distribution.
# missing     the missing data option after match.arg(). Every option
#             but full information maximum likelihood is refused: an
#             auxiliary carries information only through the full
#             information likelihood. lavaan's listwise deletion would
#             instead drop every row on which an auxiliary is missing, so
#             the coefficients would be those of a smaller sample, unlike
#             lm()'s (145 of 301 rows with t26 as the auxiliary), and
#             reliability_alpha() and reliability_H() refuse that
#             combination the same way. Under pairwise deletion
#             ("pairwise", and "available.cases", which gives the same
#             moments for these continuous variables) the moments are
#             computed from each pair of variables on its own, and with
#             the model saturated the coefficients are
#             a function of the moments of the outcomes and predictors
#             alone, so an auxiliary changes nothing: before 1.1.0 the
#             fit was the fit without it (on holzinger_swineford, t6 on
#             t7, t9, and t26 with t13 as the auxiliary, every
#             coefficient and standard error agreed to 1e-7), while the
#             result recorded the auxiliary as used.
#
# Returns NULL when 'auxiliary' is NULL, otherwise a list with:
#   aux_names    lavaan-safe names ("aux_1", "aux_2", ...)
#   aux_display  the original variable names, same order
#   aux_data     data.frame of the auxiliary columns, renamed
mlmr_prepare_auxiliary <- function(auxiliary, data, model_vars,
                                   fixed_x, missing = "fiml") {
  if (is.null(auxiliary)) {
    return(NULL)
  }
  if (!is.character(auxiliary) || length(auxiliary) == 0L ||
      anyNA(auxiliary) || !all(nzchar(auxiliary))) {
    .dmar_refuse_argument(c(
      "'auxiliary' must be NULL or a character vector of variable ",
      "names in 'data'."), arg = c("auxiliary", "data"))
  }
  if (identical(missing, "listwise")) {
    .dmar_refuse_argument(c(
      "'auxiliary' cannot be combined with missing = \"listwise\": an ",
      "auxiliary variable carries information only through full ",
      "information maximum likelihood, and listwise deletion would ",
      "instead drop every row on which an auxiliary is missing. Drop ",
      "'auxiliary' or set missing = \"fiml\"."), arg = "auxiliary")
  }
  if (missing %in% c("pairwise", "available.cases")) {
    .dmar_refuse_argument(c(
      "'auxiliary' cannot be combined with missing = \"", missing,
      "\": an auxiliary variable carries information only through full ",
      "information maximum likelihood, and pairwise deletion computes ",
      "the moments from each pair of variables on its own, so the ",
      "auxiliary would leave every estimate unchanged. Drop ",
      "'auxiliary' or set missing = \"fiml\"."), arg = "auxiliary")
  }
  auxiliary <- unique(auxiliary)
  not_found <- setdiff(auxiliary, names(data))
  if (length(not_found) > 0L) {
    .dmar_refuse_argument(c(
      "auxiliary variable(s) not found in 'data': ",
      paste(not_found, collapse = ", "), "."), arg = "data")
  }
  clash <- intersect(auxiliary, model_vars)
  if (length(clash) > 0L) {
    .dmar_refuse_argument(c(
      "auxiliary variable(s) also appear in 'formula': ",
      paste(clash, collapse = ", "),
      ". An auxiliary variable is correlated with the model, not a ",
      "term in it; remove it from one place or the other."), arg = "formula")
  }
  if (isTRUE(fixed_x)) {
    .dmar_refuse_argument(c(
      "auxiliary variables require 'fixed_x = FALSE' (the default); ",
      "the saturated correlates approach models the joint ",
      "distribution of the predictors, outcome(s), and auxiliaries."),
      arg = NULL)
  }
  aux_cols <- data[, auxiliary, drop = FALSE]
  is_num <- vapply(aux_cols, is.numeric, logical(1L))
  if (!all(is_num)) {
    .dmar_refuse_argument(c(
      "auxiliary variable(s) must be numeric: ",
      paste(auxiliary[!is_num], collapse = ", "),
      ". Encode a factor as numeric indicator(s) before passing it ",
      "as an auxiliary."), arg = NULL)
  }
  aux_names <- paste0("aux_", seq_along(auxiliary))
  aux_data <- as.data.frame(aux_cols)
  names(aux_data) <- aux_names
  list(aux_names = aux_names, aux_display = auxiliary, aux_data = aux_data)
}


# Saturated correlates syntax (Graham, 2003). Each auxiliary variable
# correlates with the residual of every outcome, with every predictor,
# and with every other auxiliary, but is not a predictor of any
# outcome, and the predictors keep their free covariances with each
# other. This recovers the auxiliary's information for the full
# information likelihood under MAR without changing the meaning of the
# focal regression coefficients. The residual covariance with each
# outcome ("<y> ~~ <aux>") is the piece lavaan does not add on its own
# and is what carries the auxiliary's information into the outcome's
# likelihood; relying on lavaan's automatic exogenous covariances
# instead shifts the slope estimates, so every correlate is stated
# explicitly here.
#
# The covariances among the predictors ("<x_j> ~~ <x_k>") are stated
# too. A predictor that covaries with an auxiliary is no longer
# exogenous in lavaan's sense, so lavaan's automatic covariances among
# exogenous observed variables (auto.cov.x) no longer reach it, and
# without these lines every covariance among the predictors is fixed at
# zero. The model is then not saturated (one degree of freedom with two
# predictors), its implied variance of the outcome omits the predictors'
# covariance, and when a predictor is missing the likelihood treats the
# predictors as uncorrelated and the slopes are biased. Versions before
# 1.1.0 omitted these lines. The variances and means of the predictors
# stay free through lavaan's defaults for sem() (auto.var, and
# int.ov.free with the mean structure), so with these lines the model
# over the outcomes, the predictors, and the auxiliaries is saturated.
mlmr_auxiliary_syntax <- function(outcome_names, X_names, aux_names) {
  lines <- character(0)
  for (a in aux_names) {
    for (yk in outcome_names) {
      lines <- c(lines, paste0(yk, " ~~ ", a))
    }
    for (xk in X_names) {
      lines <- c(lines, paste0(a, " ~~ ", xk))
    }
  }
  if (length(aux_names) > 1L) {
    cmb <- utils::combn(aux_names, 2L)
    for (j in seq_len(ncol(cmb))) {
      lines <- c(lines, paste0(cmb[1L, j], " ~~ ", cmb[2L, j]))
    }
  }
  if (length(X_names) > 1L) {
    cmb <- utils::combn(X_names, 2L)
    for (j in seq_len(ncol(cmb))) {
      lines <- c(lines, paste0(cmb[1L, j], " ~~ ", cmb[2L, j]))
    }
  }
  paste(lines, collapse = "\n")
}


# Wald confidence intervals: symmetric, with the normal quantile by
# default and the t quantile on wald_df degrees of freedom (one value, or
# one per coefficient) when a finite wald_df is given (.wald_crit(),
# R/missing_information_internals.R).
mlmr_wald_ci <- function(est, se, conf_level, wald_df = Inf) {
  crit <- .wald_crit(conf_level, wald_df)
  lo <- est - crit * se
  hi <- est + crit * se
  cbind(lower = lo, upper = hi)
}


# Profile likelihood confidence intervals.
#
# For each beta label, the function refits the model with the
# parameter constrained to a candidate value, computes the likelihood
# ratio statistic LR(b) = 2 (logL_full - logL_constrained), and
# searches for the two values of b where LR(b) equals
# qchisq(conf_level, 1). The bracket is initialized at +/- 4 SE on
# each side of the point estimate, with uniroot's extendInt argument
# enlarging the bracket if necessary. Refit failures are treated as
# Inf (i.e., "outside the confidence region"), which is appropriate
# only when failure happens far from the maximum likelihood estimate.
#
# Iterative refits emit one warning per refit when lavaan reports
# convergence issues; these are counted in a small state environment
# that the objective closure updates, and re-emitted as a single
# summary warning at the end so the caller is not flooded.
mlmr_profile_ci <- function(fit, syntax, fit_args, beta_labels,
                            est, se, conf_level) {
  ll_full <- as.numeric(lavaan::logLik(fit))
  crit <- stats::qchisq(conf_level, df = 1)

  ci <- matrix(NA_real_, nrow = length(beta_labels), ncol = 2L,
               dimnames = list(NULL, c("lower", "upper")))

  state <- new.env(parent = emptyenv())
  state$warn_count <- 0L
  refit_with_constraint <- function(label, value) {
    constrained_syntax <- paste0(syntax, "\n", label, " == ",
                                 formatC(value, format = "g",
                                         digits = 15))
    new_args <- fit_args
    new_args$model <- constrained_syntax
    new_fit <- try(suppressWarnings(do.call(lavaan::sem, new_args)),
                   silent = TRUE)
    if (inherits(new_fit, "try-error")) return(NA_real_)
    if (!lavaan::lavInspect(new_fit, "converged")) return(NA_real_)
    as.numeric(lavaan::logLik(new_fit))
  }

  for (i in seq_along(beta_labels)) {
    lab <- beta_labels[i]
    point <- est[i]
    s <- se[i]
    if (!is.finite(point) || !is.finite(s) || s <= 0) next

    obj <- function(b) {
      ll_c <- refit_with_constraint(lab, b)
      if (!is.finite(ll_c)) {
        state$warn_count <- state$warn_count + 1L
        return(Inf)
      }
      2 * (ll_full - ll_c) - crit
    }

    width <- 4 * s
    lower_root <- try(stats::uniroot(obj,
                                     lower = point - width,
                                     upper = point - 1e-8 * max(1, abs(point)),
                                     extendInt = "downX",
                                     tol = 1e-4,
                                     maxiter = 60)$root,
                      silent = TRUE)
    upper_root <- try(stats::uniroot(obj,
                                     lower = point + 1e-8 * max(1, abs(point)),
                                     upper = point + width,
                                     extendInt = "upX",
                                     tol = 1e-4,
                                     maxiter = 60)$root,
                      silent = TRUE)
    if (!inherits(lower_root, "try-error")) ci[i, 1L] <- lower_root
    if (!inherits(upper_root, "try-error")) ci[i, 2L] <- upper_root
  }

  if (state$warn_count > 0L) {
    .dmar_warn(c(
      "mlmr profile CI: ", state$warn_count,
      " constrained refits did not converge and were treated ",
      "as outside the confidence region. CI bounds for the ",
      "affected parameter may be conservative or NA."), "result_qualified")
  }
  na_sides <- which(is.na(ci), arr.ind = TRUE)
  if (nrow(na_sides) > 0L) {
    failed_terms <- unique(na_sides[, 1L])
    failed_labels <- beta_labels[failed_terms]
    .dmar_warn(c(
      "mlmr profile CI: uniroot could not bracket a sign ",
      "change for ", length(failed_labels),
      " parameter(s) (", paste(failed_labels, collapse = ", "),
      "). The corresponding CI bound(s) are reported as NA. ",
      "This often indicates the constrained model is ",
      "numerically indistinguishable from the full model over ",
      "the search bracket, or that the model is unidentified ",
      "in some direction. Consider ci_method = \"wald\" or ",
      "ci_method = \"boot\" for a finite interval."), "result_qualified")
  }
  ci
}


# Per-predictor effect sizes and omnibus likelihood ratio test.
#
# Refits the model once per slope with that slope constrained to
# zero, recording the constrained residual variance and converting to
# the semi-partial R-squared (sr2 = R2_full - R2_reduced) and Cohen's
# f-squared (f2 = sr2 / (1 - R2_full)). The overall LR test is one
# additional refit with all slopes constrained to zero
# simultaneously, giving the omnibus deviance statistic on K degrees
# of freedom (analogous to the lm() F-test of all slopes equal to
# zero).
#
# Returns a list with:
#   sr2     numeric vector, length K, per-slope semi-partial R^2
#   f2      numeric vector, length K, per-slope Cohen's f^2
#   omnibus data.frame with statistic, df, p_value for the LR test
#           of "all slopes = 0"
mlmr_compute_effect_sizes <- function(fit_args, syntax, slope_labels,
                                      sigma2_y, full_R2, ll_full) {
  K <- length(slope_labels)
  sr2 <- rep(NA_real_, K)
  f2  <- rep(NA_real_, K)

  refit_constrained <- function(constraint_lines) {
    constrained_syntax <- paste(syntax,
                                paste(constraint_lines, collapse = "\n"),
                                sep = "\n")
    new_args <- fit_args
    new_args$model <- constrained_syntax
    rfit <- try(suppressWarnings(do.call(lavaan::sem, new_args)),
                silent = TRUE)
    if (inherits(rfit, "try-error")) return(NULL)
    if (!lavaan::lavInspect(rfit, "converged")) return(NULL)
    rfit
  }

  for (j in seq_len(K)) {
    rfit <- refit_constrained(paste0(slope_labels[j], " == 0"))
    if (is.null(rfit)) next
    rpe <- lavaan::parameterEstimates(rfit, ci = FALSE)
    reduced_sigma2_e <- rpe$est[rpe$label == "sigma2_e"][1]
    reduced_R2 <- 1 - reduced_sigma2_e / sigma2_y
    sr2[j] <- full_R2 - reduced_R2
    if (full_R2 < 1) f2[j] <- sr2[j] / (1 - full_R2) else f2[j] <- Inf
  }

  null_fit <- refit_constrained(paste0(slope_labels, " == 0"))
  if (is.null(null_fit)) {
    omnibus <- data.frame(statistic = NA_real_, df = K,
                          p_value = NA_real_, stringsAsFactors = FALSE)
  } else {
    ll_null <- as.numeric(lavaan::logLik(null_fit))
    LR <- 2 * (ll_full - ll_null)
    omnibus <- data.frame(statistic = LR, df = K,
                          p_value = stats::pchisq(LR, df = K,
                                                  lower.tail = FALSE),
                          stringsAsFactors = FALSE)
  }

  list(sr2 = sr2, f2 = f2, omnibus = omnibus)
}


# Bootstrap confidence intervals.
#
# Uses lavaan::lavBootstrap with fun returning the slopes of
# interest. type = "ordinary" resamples the rows of the original
# data; type = "bollen.stine" applies the Bollen and Stine (1992)
# model-based bootstrap.
#
# boot_seed: optional integer for RNG reproducibility. If supplied,
# the seed is set for the duration of this function through
# .dmar_local_seed(), which restores the caller's generator state when
# the function exits so the user's random stream is not perturbed. If
# NULL (the default), successive calls draw fresh resamples from the
# user's current RNG.
#
# Returns the matrix of limits (one row per coefficient, columns "lower"
# and "upper") with the verdict of the rule on the replications as its
# "bootstrap" attribute (.dmar_boot_settle(); `fun` names the exported
# caller in the conditions). A replication with a coefficient that is not
# a finite number, and one lavaan did not return, has no value: it is
# dropped for every coefficient, counted, and reported, so all intervals
# rest on the same resamples. Through DMAR 1.0.x a missing value was
# removed coefficient by coefficient without a count or a warning. When
# the rule withholds the intervals the limits are NA.
mlmr_boot_ci <- function(fit, fit_args, beta_labels, conf_level,
                         B, type, boot_seed = NULL, fun = "mlmr") {
  .dmar_local_seed(boot_seed)

  boot_fun <- function(fit_b) {
    pe <- lavaan::parameterEstimates(fit_b, ci = FALSE)
    vapply(beta_labels, function(lab) {
      row <- pe[pe$label == lab, , drop = FALSE]
      if (nrow(row) != 1L) NA_real_ else row$est[[1L]]
    }, numeric(1))
  }

  # lavaan warns of the number of runs that failed or did not converge.
  # Those runs are the replications counted as dropped below, and the
  # classed warning of the rule reports them with the counts, so lavaan's
  # unclassed count is not repeated beside it; its other warnings pass.
  boot_res <- try(withCallingHandlers(
    lavaan::lavBootstrap(fit, r = B, type = type, fun = boot_fun),
    warning = function(w) {
      if (grepl("bootstrap runs? failed or did not converge",
                conditionMessage(w))) {
        invokeRestart("muffleWarning")
      }
    }), silent = TRUE)
  if (inherits(boot_res, "try-error")) {
    .dmar_refuse(c(
      "Bootstrap failed: ",
      attr(boot_res, "condition")$message), "computation_failed")
  }

  if (is.list(boot_res) && !is.data.frame(boot_res)) {
    boot_mat <- do.call(rbind, boot_res)
  } else if (is.matrix(boot_res)) {
    boot_mat <- boot_res
  } else {
    boot_mat <- as.matrix(boot_res)
  }

  alpha <- (1 - conf_level) / 2
  boot_mat <- matrix(as.numeric(boot_mat), ncol = length(beta_labels))
  keep <- rowSums(!is.finite(boot_mat)) == 0L
  verdict <- .dmar_boot_settle(
    B, sum(keep), tail = alpha, conf_level = conf_level, fun = fun,
    quantity = "some coefficient",
    defined_as = "every coefficient is defined",
    dropped_because = paste0("gave a refit that failed, did not ",
                             "converge, or returned a coefficient that ",
                             "is not a finite number"))
  ci <- matrix(NA_real_, nrow = length(beta_labels), ncol = 2L)
  if (verdict$formed) {
    ci <- t(apply(boot_mat[keep, , drop = FALSE], 2L, stats::quantile,
                  probs = c(alpha, 1 - alpha), names = FALSE))
  }
  dimnames(ci) <- list(NULL, c("lower", "upper"))
  attr(ci, "bootstrap") <- verdict
  ci
}

# The elements a fit object carries for a bootstrap interval: B, B_used,
# n_dropped, B_required, and bootstrap_status, from the verdict
# mlmr_boot_ci() attaches. NULL when no bootstrap was run.
mlmr_boot_counts <- function(verdict) {
  if (is.null(verdict)) return(NULL)
  list(B = verdict$B, B_used = verdict$B_used,
       n_dropped = verdict$n_dropped, B_required = verdict$B_required,
       bootstrap_status = verdict$status)
}

# The line print() and summary() show for a bootstrap interval
# (.dmar_boot_footer()); nothing is written when no bootstrap was run.
mlmr_boot_line <- function(x) {
  if (is.null(x$bootstrap_status)) return(invisible(NULL))
  cat(.dmar_boot_footer(x[c("B", "B_used", "n_dropped", "B_required",
                            "bootstrap_status")]), "\n", sep = "")
  invisible(NULL)
}
