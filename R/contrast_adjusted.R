#' Confidence Interval for a Contrast of Covariate-Adjusted Cell Means in a Factorial ANCOVA
#'
#' @description
#' Given a fitted \code{\link[stats]{lm}} or \code{\link[stats]{aov}} object
#' for a factorial analysis of covariance (one or more crossed factors plus one
#' or more covariates) and a numeric contrast vector over the cells of the
#' factorial design, \code{contrast_adjusted()} forms the contrast of the
#' covariate-adjusted cell means,
#' \eqn{\hat{\psi} = \sum_j c_j \, \hat{\bar{Y}}_j}, where each
#' \eqn{\hat{\bar{Y}}_j} is the model's predicted mean for cell \eqn{j}
#' evaluated at the mean of every covariate (the adjusted, or least-squares,
#' cell mean). It returns the point estimate, a \emph{t} confidence interval on
#' the model's residual degrees of freedom, and the accompanying \emph{t}
#' statistic and two-sided \emph{p}-value for \eqn{H_0\!: \psi = 0}.
#'
#' @details
#' The adjusted cell means are the means the ANCOVA actually tests: the
#' predicted outcome for each combination of factor levels, holding every
#' covariate at its sample mean. Writing \eqn{L} for the linear map that sends
#' the model coefficients to that contrast of adjusted means (built by
#' evaluating the model's design matrix at each cell with the covariates set to
#' their means and combining the rows with the contrast weights), the point
#' estimate is \eqn{\hat{\psi} = L' \hat{\beta}} and its standard error is the
#' square root of the quadratic form \eqn{L' \, \mathrm{vcov}(\hat{\beta}) \, L}.
#' The interval is \eqn{\hat{\psi} \pm t_{1 - \alpha/2,\, \nu}\, \mathrm{SE}},
#' with \eqn{\nu} the residual degrees of freedom of the fitted model.
#'
#' Because \eqn{L} is read off the fitted model's own design matrix, the
#' function is agnostic to how the factors are parameterized: a cell-means
#' parameterization (\code{y ~ 0 + cell + x}) and the crossed-factor
#' parameterization (\code{y ~ A * B + x}) give the same contrast estimate and
#' standard error, provided the contrast vector is ordered to match the cells
#' of the reference grid (see \code{Note}).
#'
#' The \emph{t} interval returned here has exact per-comparison coverage for a
#' single contrast chosen in advance. For a family of contrasts examined
#' together, adjust the critical value for multiplicity (for example the
#' Scheffe critical value for the full cell space, \code{\link{cv_scheffe}}, or
#' a Bryant--Paulson simultaneous interval, \code{\link{ci_c_ancova_bp}}).
#'
#' @param model A fitted \code{\link[stats]{lm}} or \code{\link[stats]{aov}}
#'   object for a factorial ANCOVA: one or more crossed factors and one or more
#'   numeric covariates on the right-hand side of the formula.
#' @param contrast A numeric vector of contrast weights, one weight per cell of
#'   the factorial design (the crossing of the model's factors). Its length must
#'   equal the number of cells. The weights typically sum to zero.
#' @param conf_level The confidence level for the interval (default \code{0.95}).
#'
#' @return
#' A five-row \code{dmar_tbl} (a \code{data.frame} with columns \code{term} and
#' \code{value}). The \code{term} values are \code{"contrast"} (the point
#' estimate \eqn{\hat{\psi}} of the contrast of adjusted cell means),
#' \code{"lower_limit"} and \code{"upper_limit"} (the confidence limits),
#' \code{"t"} (the \emph{t} statistic), and \code{"p"} (the two-sided
#' \emph{p}-value). The stored \code{value} column is numeric at full precision.
#'
#' @note
#' The contrast weights are matched to the cells of the reference grid, which is
#' the crossing of the model's factors in the order the factors appear in the
#' model formula, with the first factor varying fastest (the order
#' \code{\link[base]{expand.grid}} produces over the factor levels). For a
#' single factor this is simply the order of its levels. When in doubt, fit the
#' cell-means form \code{y ~ 0 + cell + covariates} with
#' \code{cell = interaction(A, B, ...)} and order the weights to match
#' \code{levels(cell)}.
#'
#' @references
#' Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). \emph{Designing
#'   experiments and analyzing data: A model comparison perspective}
#'   (4th ed.). Routledge. (See Chapter 9 on designs with covariates.)
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso
#' \code{\link{contrast_test}} for contrasts of unadjusted group means in a
#' one-way design; \code{\link{ci_c_ancova}} for a single-covariate ANCOVA
#' contrast from summary statistics.
#'
#' @examples
#' # A 2 x 2 factorial ANCOVA with one covariate.
#' set.seed(113)
#' d <- data.frame(
#'   A = factor(rep(c("a1", "a2"), each = 40)),
#'   B = factor(rep(rep(c("b1", "b2"), each = 20), 2)),
#'   x = rnorm(80)
#' )
#' d$y <- 5 + 2 * (d$A == "a2") + 1.5 * (d$B == "b2") +
#'   0.8 * d$x + rnorm(80)
#' fit <- lm(y ~ A * B + x, data = d)
#'
#' # Cells in reference-grid order: (a1,b1), (a2,b1), (a1,b2), (a2,b2).
#' # Main effect of A, averaged over B: mean(a2 cells) - mean(a1 cells).
#' contrast_adjusted(fit, contrast = c(-0.5, 0.5, -0.5, 0.5))
#'
#' @keywords htest design
#'
#' @family confidence intervals for effect sizes
#'
#' @export
#' @import stats
contrast_adjusted <- function(model, contrast, conf_level = 0.95) {
  if (!inherits(model, c("lm", "aov"))) {
    stop("'model' must be a fitted lm or aov object.")
  }
  if (!is.numeric(conf_level) || length(conf_level) != 1L ||
      conf_level <= 0 || conf_level >= 1) {
    stop("'conf_level' must be a single number strictly between 0 and 1.")
  }
  if (!is.numeric(contrast) || !is.null(dim(contrast))) {
    stop("'contrast' must be a numeric vector of weights, one per cell.")
  }

  # Split the model frame into factors (which define the cells) and numeric
  # covariates (held at their means to obtain adjusted cell means).
  mf  <- stats::model.frame(model)
  trm <- stats::terms(model)
  resp <- attr(trm, "response")
  pred_vars <- setdiff(seq_len(ncol(mf)), resp)
  is_fac <- vapply(pred_vars, function(j) {
    v <- mf[[j]]
    is.factor(v) || is.character(v) || is.logical(v)
  }, logical(1L))
  fac_cols <- pred_vars[is_fac]
  cov_cols <- pred_vars[!is_fac]
  if (length(fac_cols) == 0L) {
    stop("'model' must contain at least one factor to define the cells.")
  }
  if (length(cov_cols) == 0L) {
    stop("'model' must contain at least one covariate (a numeric predictor); ",
         "for a contrast of unadjusted group means use contrast_test().")
  }

  # Reference grid: the crossing of the factor levels, first factor fastest
  # (the order expand.grid() produces). Covariates enter at their sample means.
  fac_levels <- lapply(fac_cols, function(j) {
    v <- mf[[j]]
    if (is.factor(v)) levels(v) else sort(unique(v))
  })
  names(fac_levels) <- names(mf)[fac_cols]
  grid <- expand.grid(fac_levels, KEEP.OUT.ATTRS = FALSE,
                      stringsAsFactors = FALSE)
  # Restore the original factor typing so the design matrix uses the model's
  # contrasts and level ordering.
  for (nm in names(grid)) {
    v <- mf[[nm]]
    if (is.factor(v)) grid[[nm]] <- factor(grid[[nm]], levels = levels(v))
  }
  for (j in cov_cols) {
    nm <- names(mf)[j]
    grid[[nm]] <- mean(mf[[j]])
  }

  n_cells <- nrow(grid)
  if (length(contrast) != n_cells) {
    stop("'contrast' has length ", length(contrast), " but the design has ",
         n_cells, " cells; supply one weight per cell.")
  }

  # Design-matrix rows for the reference grid, then the coefficient-space map
  # L = t(X_cells) %*% contrast so that psi = L' beta and its variance is the
  # quadratic form L' vcov(model) L. Non-estimable coefficients (aliased
  # columns dropped by lm) carry NA; they must not appear in the contrast.
  X_cells <- stats::model.matrix(stats::delete.response(trm), data = grid,
                                 contrasts.arg = model$contrasts)
  beta <- stats::coef(model)
  if (ncol(X_cells) != length(beta)) {
    # Align columns to the fitted coefficient names, guarding against aliasing.
    common <- intersect(colnames(X_cells), names(beta))
    X_cells <- X_cells[, common, drop = FALSE]
    beta <- beta[common]
  }
  keep <- !is.na(beta)
  # Rank-deficient fit: an empty factorial cell aliases a coefficient, which lm
  # drops to NA. The contrast of adjusted means is estimable only when its full
  # coefficient-space image is invariant to that aliasing, that is, the weight
  # the aliased columns carry equals the weight the kept columns already imply.
  # Otherwise lm's arbitrary zero for the aliased coefficient would make the
  # reported value an artifact (emmeans flags such a contrast as non-estimable).
  # The reference grid itself is full rank, so estimability is judged against
  # the fitting design matrix, not against X_cells.
  if (!all(keep)) {
    L_full   <- drop(crossprod(X_cells, contrast))
    X_fit    <- stats::model.matrix(model)[, colnames(X_cells), drop = FALSE]
    aliasing <- qr.solve(X_fit[, keep, drop = FALSE],
                         X_fit[, !keep, drop = FALSE])
    resid_estim <- L_full[!keep] - drop(crossprod(aliasing, L_full[keep]))
    if (max(abs(resid_estim)) > 1e-8 * max(1, max(abs(L_full)))) {
      stop("The requested contrast is not estimable: it places weight on an ",
           "aliased cell of a rank-deficient design, such as an empty ",
           "factorial cell. No covariate-adjusted value exists for this ",
           "contrast; check the design and the contrast weights.",
           call. = FALSE)
    }
  }
  L <- drop(crossprod(X_cells[, keep, drop = FALSE], contrast))
  V <- stats::vcov(model)[keep, keep, drop = FALSE]

  psi <- sum(L * beta[keep])
  se  <- sqrt(drop(t(L) %*% V %*% L))
  nu  <- stats::df.residual(model)
  t_value <- psi / se
  p_value <- 2 * stats::pt(-abs(t_value), df = nu)

  alpha <- 1 - conf_level
  t_crit <- stats::qt(1 - alpha / 2, df = nu)
  lo <- psi - t_crit * se
  hi <- psi + t_crit * se

  out <- data.frame(
    term  = c("contrast", "lower_limit", "upper_limit", "t", "p"),
    value = c(psi, lo, hi, t_value, p_value)
  )
  .as_dmar_tbl(out, conf_level = conf_level, p_terms = "p")
}
