#' Internal Helpers for mlmr_mv
#'
#' Not exported. The helpers translate a multivariate formula into
#' a lavaan-friendly model frame, construct the lavaan syntax for
#' the joint multivariate regression, and assemble per-outcome
#' confidence intervals and effect sizes.
#'
#' @keywords internal
#' @name mlmr_mv-internals
NULL


# Expand a multivariate formula into a lavaan-friendly model frame.
#
# Returns a list with:
#   Y_names       internal lavaan-safe outcome names ("y_1", "y_2", ...)
#   Y_display     display names for the outcomes from the cbind()
#                 call (e.g., "mpg", "disp")
#   X_names       internal lavaan-safe predictor names
#                 ("x_1", "x_2", ...)
#   X_display     display names from model.matrix() (e.g., "wt",
#                 "groupB", "wt:hp")
#   has_intercept whether the formula has an intercept
#   model_data    data.frame the lavaan fit consumes, columns
#                 c(Y_names, X_names), missing values preserved
#   model_frame   the original model.frame (na.action = na.pass)
#   terms         the terms object
mlmr_mv_build_design <- function(formula, data) {
  trms <- stats::terms(formula, data = data)
  mf <- stats::model.frame(trms, data = data,
                           na.action = stats::na.pass)
  has_intercept <- as.logical(attr(trms, "intercept"))

  Y_block <- stats::model.response(mf)
  if (is.null(Y_block)) {
    stop("'formula' must have a response on the left-hand side.",
         call. = FALSE)
  }
  if (!is.matrix(Y_block)) {
    stop("'formula' must use cbind() on the left-hand side for ",
         "multivariate regression, e.g., cbind(y1, y2) ~ x1 + x2. ",
         "For a single outcome, use mlmr() instead.", call. = FALSE)
  }
  J <- ncol(Y_block)
  if (J < 2L) {
    stop("mlmr_mv requires two or more outcomes; got ", J, ".",
         call. = FALSE)
  }
  Y_display <- colnames(Y_block)
  Y_names <- paste0("y_", seq_len(J))

  X_full <- stats::model.matrix(trms, mf)
  X <- if (has_intercept) X_full[, -1L, drop = FALSE] else X_full
  X_display <- colnames(X)
  X_names <- if (ncol(X) > 0L) paste0("x_", seq_len(ncol(X)))
             else character(0)

  model_data <- as.data.frame(Y_block)
  names(model_data) <- Y_names
  if (length(X_names) > 0L) {
    X_df <- as.data.frame(X)
    names(X_df) <- X_names
    model_data <- cbind(model_data, X_df)
  }

  list(
    Y_names = Y_names,
    Y_display = Y_display,
    X_names = X_names,
    X_display = X_display,
    has_intercept = has_intercept,
    model_data = model_data,
    model_frame = mf,
    terms = trms
  )
}


# Build the lavaan model syntax for a multivariate regression.
#
# For outcomes (y_1, y_2) and predictors (x_1, x_2):
#   y_1 ~ b0_y_1*1 + b_y_1_1*x_1 + b_y_1_2*x_2
#   y_2 ~ b0_y_2*1 + b_y_2_1*x_1 + b_y_2_2*x_2
#   y_1 ~~ sigma2_e_y_1*y_1
#   y_2 ~~ sigma2_e_y_2*y_2
#   y_1 ~~ sigma_e_y_1_y_2*y_2
#
# The labels are chosen so that the per-outcome regression
# coefficients (b0_<y>, b_<y>_<k>) can be referenced individually in
# constraints during profile-CI and effect size computation, and so
# that the residual variances and covariances can be extracted by
# name from parameterEstimates().
mlmr_mv_build_syntax <- function(Y_names, X_names, has_intercept) {
  J <- length(Y_names)
  K <- length(X_names)

  reg_lines <- vapply(Y_names, function(y) {
    slope_terms <- if (K > 0L) {
      paste0("b_", y, "_", seq_len(K), "*", X_names,
             collapse = " + ")
    } else {
      ""
    }
    intercept_piece <- if (has_intercept) paste0("b0_", y, "*1") else "0*1"
    rhs <- if (nzchar(slope_terms)) {
      paste(intercept_piece, slope_terms, sep = " + ")
    } else {
      intercept_piece
    }
    paste0(y, " ~ ", rhs)
  }, character(1L))

  resid_var_lines <- vapply(Y_names, function(y) {
    paste0(y, " ~~ sigma2_e_", y, "*", y)
  }, character(1L))

  resid_cov_lines <- character(0)
  if (J >= 2L) {
    for (a in seq_len(J - 1L)) {
      for (b in seq.int(a + 1L, J)) {
        resid_cov_lines <- c(
          resid_cov_lines,
          paste0(Y_names[a], " ~~ sigma_e_", Y_names[a], "_",
                 Y_names[b], "*", Y_names[b])
        )
      }
    }
  }

  paste(c(reg_lines, resid_var_lines, resid_cov_lines),
        collapse = "\n")
}


# Vector of beta labels for one outcome (b0_y, b_y_1, b_y_2, ...).
mlmr_mv_beta_labels <- function(y_name, K, has_intercept) {
  c(if (has_intercept) paste0("b0_", y_name),
    paste0("b_", y_name, "_", seq_len(K)))
}


# Per-outcome confidence intervals. Returns an array of dimension
# (P, 2, J) where P = K + intercept, 2 = (lower, upper), J = outcomes.
mlmr_mv_compute_ci <- function(fit, syntax, fit_args,
                               beta_labels_all, est_mat, se_mat,
                               conf_level, ci_method, B,
                               boot_type, boot_seed) {
  J <- ncol(est_mat)
  P <- nrow(est_mat)
  out <- array(NA_real_, dim = c(P, 2L, J),
               dimnames = list(rownames(est_mat),
                               c("lower", "upper"),
                               colnames(est_mat)))

  if (ci_method == "wald") {
    z_crit <- stats::qnorm(1 - (1 - conf_level) / 2)
    out[, 1L, ] <- est_mat - z_crit * se_mat
    out[, 2L, ] <- est_mat + z_crit * se_mat
    return(out)
  }

  if (ci_method == "profile") {
    # Per-outcome profile CIs reusing mlmr_profile_ci on the joint
    # fit. The helper does not need to know that the syntax has
    # multiple regressions; it constrains one label at a time.
    for (j in seq_len(J)) {
      P_j <- P
      labs_j <- beta_labels_all[((j - 1L) * P_j + 1L):(j * P_j)]
      ci_j <- mlmr_profile_ci(
        fit = fit, syntax = syntax, fit_args = fit_args,
        beta_labels = labs_j,
        est = est_mat[, j], se = se_mat[, j],
        conf_level = conf_level
      )
      out[, , j] <- ci_j
    }
    return(out)
  }

  # ci_method == "boot": one bootstrap run that returns the full
  # vector of slopes (and intercepts) for every outcome.
  boot_ci_long <- mlmr_boot_ci(
    fit = fit, fit_args = fit_args,
    beta_labels = beta_labels_all,
    conf_level = conf_level, R = B, type = boot_type,
    boot_seed = boot_seed
  )
  # Reshape: boot_ci_long has length(beta_labels_all) rows; map back
  # to (P, J) layout.
  for (j in seq_len(J)) {
    idx <- ((j - 1L) * P + 1L):(j * P)
    out[, 1L, j] <- boot_ci_long[idx, 1L]
    out[, 2L, j] <- boot_ci_long[idx, 2L]
  }
  out
}


# Per-outcome effect sizes (semi-partial R^2 and Cohen's f^2) via
# K constrained refits per outcome, fixing one slope at a time to
# zero. Returns a long data.frame:
#   outcome | term | sr2 | f2
mlmr_mv_compute_effect_sizes <- function(fit_args, syntax,
                                         Y_names, Y_display,
                                         X_display, K,
                                         sigma2_y, full_R2) {
  J <- length(Y_names)
  rows <- vector("list", J * K)
  idx <- 1L

  refit_constrained <- function(label) {
    constrained_syntax <- paste(syntax, paste0(label, " == 0"),
                                sep = "\n")
    new_args <- fit_args
    new_args$model <- constrained_syntax
    rfit <- try(suppressWarnings(do.call(lavaan::sem, new_args)),
                silent = TRUE)
    if (inherits(rfit, "try-error")) return(NULL)
    if (!lavaan::lavInspect(rfit, "converged")) return(NULL)
    rfit
  }

  for (j in seq_len(J)) {
    y_j <- Y_names[j]
    for (k in seq_len(K)) {
      lab <- paste0("b_", y_j, "_", k)
      rfit <- refit_constrained(lab)
      sr2 <- f2 <- NA_real_
      if (!is.null(rfit)) {
        rpe <- lavaan::parameterEstimates(rfit, ci = FALSE)
        reduced_sigma2_e <- rpe$est[rpe$label ==
                                    paste0("sigma2_e_", y_j)][1L]
        if (!is.na(reduced_sigma2_e) && !is.na(sigma2_y[j])) {
          reduced_R2 <- 1 - reduced_sigma2_e / sigma2_y[j]
          sr2 <- full_R2[j] - reduced_R2
          if (!is.na(full_R2[j]) && full_R2[j] < 1) {
            f2 <- sr2 / (1 - full_R2[j])
          } else {
            f2 <- Inf
          }
        }
      }
      rows[[idx]] <- data.frame(
        outcome = Y_display[j],
        term    = X_display[k],
        sr2     = sr2,
        f2      = f2,
        stringsAsFactors = FALSE
      )
      idx <- idx + 1L
    }
  }
  do.call(rbind, rows)
}
