# Covariate-adjusted cell and marginal means from a fitted linear model.
#' Adjusted Cell and Marginal Means From a Fitted Linear Model
#'
#' @description
#' Given a fitted \code{\link[stats]{lm}} or \code{\link[stats]{aov}} object
#' with one or more factors among its predictors, \code{adjusted_means()}
#' returns the means the model actually compares, sometimes called
#' least-squares means or estimated marginal means. By default the table has
#' one row per cell of the crossed factor design, each cell's mean being the
#' model's predicted response at that combination of factor levels with every
#' covariate held at its sample mean (the adjusted cell means of an ANCOVA;
#' for a model without covariates, the model-based cell means). Naming one or
#' more factors in \code{by} instead returns the marginal means of those
#' factors, formed by averaging the cell predictions over the remaining
#' factors with either equal or frequency-proportional weights. Every mean is
#' accompanied by its standard error and a \emph{t} confidence interval on the
#' model's residual degrees of freedom.
#'
#' @details
#' \strong{The reference grid and adjusted cell means.} The reference grid is
#' the crossing of the model's factor levels, enumerated in the order the
#' factors appear in the model formula with the first factor varying fastest
#' (the order \code{\link[base]{expand.grid}} produces). This is the same cell
#' order \code{\link{contrast_adjusted}} expects, so contrast weights can be
#' read off this table row by row. Every covariate enters the grid at its
#' sample mean, and a transformed covariate is evaluated by applying the
#' transformation to the mean of the raw variable: with \code{log(x)} in the
#' formula the grid carries \code{mean(x)} and the model matrix applies
#' \code{log()}, and a \code{poly(x, 2)} basis is evaluated at \eqn{\bar{x}},
#' matching \code{\link[stats]{predict}} on new data at the covariate mean.
#' Writing \eqn{L} for the matrix whose rows are the design-matrix rows of the
#' grid cells, the cell means are \eqn{L \hat{\beta}}, each standard error is
#' the square root of the corresponding diagonal element of
#' \eqn{L \, \mathrm{vcov}(\hat{\beta}) \, L'}, and each interval is the
#' \emph{t} interval on the model's residual degrees of freedom.
#'
#' \strong{Marginal means and the two weightings.} With \code{by}, the cell
#' predictions are averaged over the factors not named there, and the
#' averaging happens in the coefficient map itself: the marginal mean's
#' \eqn{L} row is the weighted average of its cells' rows, so the estimate and
#' the standard error both follow from one linear function of the
#' coefficients. \code{weights = "equal"} weights every combination of the
#' averaged-over factors equally; this is the population marginal mean of
#' Searle, Speed, and Milliken (1980), the mean for a population in which
#' every cell is equally represented regardless of the sample's cell sizes.
#' \code{weights = "proportional"} weights each averaged-over combination by
#' its observed frequency (in a weighted fit, by its total prior weight), so
#' the marginal mean targets a population whose margins are shaped like the
#' sample's. With balanced data the two weightings coincide; with unbalanced
#' data they generally differ, and the choice between them is a substantive
#' question about the population of interest, not a technical one (Maxwell,
#' Delaney, and Kelley, 2027, Chapter 7).
#'
#' \strong{Nonestimable means.} When the fitted design is rank deficient (for
#' example an empty factorial cell), the model has no predicted value for the
#' affected cell, and a marginal mean that averages over such a cell does not
#' exist either. \code{adjusted_means()} refuses with an error naming the
#' affected rows rather than reporting a value contaminated by \code{lm}'s
#' arbitrary zero for the aliased coefficient.
#'
#' \strong{Scope.} The function covers single-stratum \code{lm} and
#' \code{aov} fits with a single response. Multi-stratum \code{aovlist} fits
#' (within-subjects designs fit with an \code{Error()} term) are refused,
#' because a within-subjects marginal mean takes its standard error from the
#' matching error stratum, which this function does not compute. Factors must
#' enter the model as variables in the data, not as conversions inside the
#' formula: \code{y ~ factor(g) + x} is refused, so convert \code{g} in the
#' data first.
#'
#' @param model A fitted \code{\link[stats]{lm}} or \code{\link[stats]{aov}}
#'   object with one or more factors (and optionally covariates) on the
#'   right-hand side of the formula.
#' @param by \code{NULL} (default) for the cell means table, or a character
#'   vector naming one or more of the model's factors for their marginal
#'   means. The output rows cross the named factors in the order given, first
#'   factor varying fastest.
#' @param weights Weighting used to average cell predictions into marginal
#'   means, so it matters only when \code{by} is supplied and the data are
#'   unbalanced. \code{"equal"} (default) weights every combination of the
#'   averaged-over factors equally; \code{"proportional"} weights each
#'   combination by its observed frequency.
#' @param conf_level The confidence level for the intervals (default
#'   \code{0.95}).
#'
#' @return
#' A \code{data.frame} (class \code{dmar_tbl}) with one row per cell of the
#' reference grid or, with \code{by}, one row per combination of the named
#' factors. The leading columns give the factor levels; the numeric columns
#' are \code{estimate} (the adjusted mean), \code{se} (its standard error),
#' and \code{ci_lower} / \code{ci_upper} (the \emph{t} confidence limits).
#' The residual degrees of freedom of the intervals are attached as the
#' \code{df_residual} attribute and, when \code{by} is supplied, the
#' weighting as the \code{weights} attribute. The stored values keep full
#' precision; only the display rounds (see \code{\link{dmar_tbl}}).
#'
#' @references
#' Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). \emph{Designing
#'   experiments and analyzing data: A model comparison perspective}
#'   (4th ed.). Routledge. (See Chapter 7 on nonorthogonal factorial designs
#'   and Chapter 9 on designs with covariates.)
#'
#' Searle, S. R., Speed, F. M., & Milliken, G. A. (1980). Population marginal
#'   means in the linear model: An alternative to least squares means.
#'   \emph{The American Statistician, 34}(4), 216--221.
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso
#' \code{\link{contrast_adjusted}} for a confidence interval on a single
#' contrast of the adjusted cell means; \code{\link{ancova}} for the one-way
#' ANCOVA table; \code{\link{ci_dunnett}} for simultaneous many-to-one
#' comparisons.
#'
#' @examples
#' # 1. Cell means of a 2 x 3 factorial (no covariate): end-of-study IQ in
#' #    the pygmalion expectancy experiment, grades 1 through 3. The grade
#' #    factor is created in the data, not inside the formula.
#' pyg <- subset(pygmalion, grade <= 3)
#' pyg$grade <- factor(pyg$grade)
#' fit <- lm(iq_8 ~ treatment * grade, data = pyg)
#' adjusted_means(fit)
#'
#' # 2. Marginal means of grade, averaging the cell means over treatment.
#' adjusted_means(fit, by = "grade")
#'
#' # 3. An ANCOVA: each adjusted mean holds the covariate, here the pretest
#' #    depression score, at its sample mean.
#' fit_ancova <- lm(bdi_post ~ condition + bdi_pre, data = depression_bdi)
#' adjusted_means(fit_ancova)
#'
#' # 4. With unbalanced cells (few bloomers in every grade) the two
#' #    weightings answer different questions.
#' adjusted_means(fit, by = "treatment")
#' adjusted_means(fit, by = "treatment", weights = "proportional")
#'
#' @keywords models design
#'
#' @family hypothesis tests
#'
#' @export
adjusted_means <- function(model, by = NULL,
                           weights = c("equal", "proportional"),
                           conf_level = 0.95) {
  if (inherits(model, "aovlist")) {
    stop("'model' is a multi-stratum aov fit (a within-subjects design with ",
         "an Error() term). adjusted_means() covers single-stratum lm and ",
         "aov fits only: a within-subjects marginal mean takes its standard ",
         "error from the matching error stratum, which this function does ",
         "not compute.", call. = FALSE)
  }
  if (!inherits(model, c("lm", "aov")) || inherits(model, "glm")) {
    stop("'model' must be a fitted lm or aov object.")
  }
  if (inherits(model, "mlm")) {
    stop("'model' is a multivariate lm fit; adjusted_means() covers a ",
         "single response.")
  }
  weights <- match.arg(weights)
  if (!is.numeric(conf_level) || length(conf_level) != 1L ||
      conf_level <= 0 || conf_level >= 1) {
    stop("'conf_level' must be a single number strictly between 0 and 1.")
  }

  trm <- stats::terms(model)
  if (length(attr(trm, "offset")) || !is.null(model$offset)) {
    stop("'model' contains an offset term, which adjusted_means() does not ",
         "support.")
  }
  dtrm <- stats::delete.response(trm)
  mf   <- stats::model.frame(model)

  # Classify the model's variables. Factors (with characters and logicals
  # treated as factors, as the model matrix treats them) define the cells of
  # the reference grid; numeric variables are covariates, held at their
  # sample means. A transformed covariate (log(x), poly(x, 2)) is evaluated
  # at the mean of the raw variable: the grid carries the raw variable and
  # model.matrix() applies the basis, exactly as predict() would on new data.
  var_exprs <- as.list(attr(dtrm, "variables"))[-1L]
  var_names <- vapply(var_exprs, deparse1, character(1L))
  is_cat <- vapply(var_names, function(nm) {
    v <- mf[[nm]]
    is.factor(v) || is.character(v) || is.logical(v)
  }, logical(1L))

  fac_names <- var_names[is_cat]
  if (length(fac_names) == 0L) {
    stop("'model' must contain at least one factor to define the cells.")
  }
  bad_fac <- !vapply(var_exprs[is_cat], is.symbol, logical(1L))
  if (any(bad_fac)) {
    stop("The categorical predictor '", fac_names[bad_fac][1L], "' is ",
         "created inside the model formula. Convert the variable in the ",
         "data before fitting (for example d$g <- factor(d$g)) so the ",
         "reference grid can be built from the variable itself.",
         call. = FALSE)
  }

  cov_exprs <- var_exprs[!is_cat]
  cov_raw   <- unique(unlist(lapply(cov_exprs, all.vars)))
  clash <- intersect(cov_raw, fac_names)
  if (length(clash)) {
    stop("The variable '", clash[1L], "' enters the model both as a factor ",
         "and inside a numeric term, so no reference grid can hold it at ",
         "both a factor level and a covariate mean.", call. = FALSE)
  }

  # Raw covariate values. When every covariate enters untransformed, the
  # model frame already holds the raw variables. Otherwise recover the data
  # the model was fit to (re-evaluating the data argument of the call, as
  # predict-style machinery does) and align its rows to the model frame by
  # row name, which drops any rows removed by subset= or na.action.
  cov_means <- numeric(0L)
  if (length(cov_raw)) {
    raw <- if (all(vapply(cov_exprs, is.symbol, logical(1L)))) mf else NULL
    if (is.null(raw)) {
      raw <- tryCatch({
        data_call <- stats::getCall(model)$data
        r <- if (is.null(data_call)) stats::get_all_vars(dtrm) else
          stats::get_all_vars(dtrm,
                              data = eval(data_call, envir = environment(trm)))
        r[rownames(mf), , drop = FALSE]
      }, error = function(e) NULL)
      if (is.null(raw) || !all(cov_raw %in% names(raw)) ||
          anyNA(raw[cov_raw])) {
        stop("A covariate enters the model through a transformation, and ",
             "the data used to fit 'model' could not be recovered to ",
             "compute the raw covariate mean. Refit the model with its data ",
             "available, or add the transformed covariate to the data and ",
             "use it untransformed in the formula.", call. = FALSE)
      }
    }
    cov_means <- vapply(cov_raw, function(nm) {
      v <- raw[[nm]]
      if (!is.numeric(v) || !is.null(dim(v))) {
        stop("The covariate '", nm, "' is not a numeric vector; ",
             "adjusted_means() holds numeric covariates at their means.",
             call. = FALSE)
      }
      mean(v)
    }, numeric(1L))
  }

  # Reference grid: the crossing of the factor levels, first factor fastest
  # (the order expand.grid() produces), with the model's own level sets so
  # the design matrix uses the fit's contrasts and level ordering.
  fac_levels <- lapply(fac_names, function(nm) {
    v <- mf[[nm]]
    if (is.factor(v)) levels(v)
    else if (is.logical(v)) c(FALSE, TRUE)
    else if (!is.null(model$xlevels[[nm]])) model$xlevels[[nm]]
    else sort(unique(v))
  })
  names(fac_levels) <- fac_names
  grid <- expand.grid(fac_levels, KEEP.OUT.ATTRS = FALSE,
                      stringsAsFactors = FALSE)
  for (nm in fac_names) {
    v <- mf[[nm]]
    if (is.factor(v)) {
      grid[[nm]] <- factor(grid[[nm]], levels = levels(v),
                           ordered = is.ordered(v))
    }
  }
  for (nm in cov_raw) grid[[nm]] <- cov_means[[nm]]
  n_cells <- nrow(grid)

  # Design-matrix rows for the grid, mirroring predict.lm(): the terms'
  # predvars evaluate any covariate basis (poly, splines) with the
  # coefficients stored at fit time, and xlev / contrasts reproduce the
  # fit's factor coding.
  mf_grid <- stats::model.frame(dtrm, grid, xlev = model$xlevels)
  X_grid  <- stats::model.matrix(dtrm, mf_grid,
                                 contrasts.arg = model$contrasts)
  beta <- stats::coef(model)
  if (!identical(colnames(X_grid), names(beta))) {
    if (!all(names(beta) %in% colnames(X_grid))) {
      stop("The reference grid's design matrix does not carry every fitted ",
           "coefficient; adjusted_means() cannot map the grid onto this ",
           "model.")
    }
    X_grid <- X_grid[, names(beta), drop = FALSE]
  }

  # The linear map L: one row per output mean. For the cell table L is the
  # grid's design matrix itself; for a marginal table each row is the
  # weighted average of its cells' rows, so estimate and standard error both
  # come from the same linear function of the coefficients.
  if (is.null(by)) {
    L <- X_grid
    label_src <- grid[fac_names]
  } else {
    if (!is.character(by) || length(by) < 1L || anyNA(by)) {
      stop("'by' must be a character vector of factor names.")
    }
    by <- unique(by)
    if (!all(by %in% fac_names)) {
      stop("'by' must name factors in the model; the model's factor",
           if (length(fac_names) > 1L) "s are " else " is ",
           paste0("'", fac_names, "'", collapse = ", "), ".")
    }
    over <- setdiff(fac_names, by)
    out_grid <- expand.grid(fac_levels[by], KEEP.OUT.ATTRS = FALSE,
                            stringsAsFactors = FALSE)
    # Cell weights for the averaging. Equal weighting gives every
    # averaged-over combination the same weight (the population marginal
    # means of Searle, Speed, and Milliken, 1980). Proportional weighting
    # gives each averaged-over combination its observed frequency (its total
    # prior weight, in a weighted fit), computed once from the data's
    # margins, so the same weights apply within every level of 'by'.
    if (weights == "proportional" && length(over)) {
      w_obs <- stats::weights(model)
      if (is.null(w_obs)) w_obs <- rep(1, nrow(mf))
      key_obs <- do.call(paste,
                         c(unname(lapply(mf[over], as.character)),
                           list(sep = "\r")))
      freq <- tapply(w_obs, key_obs, sum)
      key_grid <- do.call(paste,
                          c(unname(lapply(grid[over], as.character)),
                            list(sep = "\r")))
      w_cell <- as.numeric(freq[key_grid])
      w_cell[is.na(w_cell)] <- 0
    } else {
      w_cell <- rep(1, n_cells)
    }
    key_by_grid <- do.call(paste,
                           c(unname(lapply(grid[by], as.character)),
                             list(sep = "\r")))
    key_by_out <- do.call(paste,
                          c(unname(lapply(out_grid, as.character)),
                            list(sep = "\r")))
    L <- matrix(0, nrow(out_grid), ncol(X_grid),
                dimnames = list(NULL, colnames(X_grid)))
    for (i in seq_len(nrow(out_grid))) {
      idx <- which(key_by_grid == key_by_out[i])
      w_i <- w_cell[idx]
      if (sum(w_i) <= 0) {
        stop("No observations fall in the margin being averaged over, so ",
             "the proportional marginal mean is undefined; use ",
             "weights = \"equal\" or check the data.", call. = FALSE)
      }
      L[i, ] <- colSums(X_grid[idx, , drop = FALSE] * (w_i / sum(w_i)))
    }
    label_src <- out_grid
  }

  label_df <- as.data.frame(lapply(label_src, as.character),
                            stringsAsFactors = FALSE, check.names = FALSE)
  lab_parts <- Map(function(nm, v) paste0(nm, " = ", v),
                   names(label_df), label_df)
  labels_chr <- do.call(paste, c(unname(lab_parts), list(sep = ", ")))

  # Rank-deficient fit: an empty factorial cell aliases a coefficient, which
  # lm drops to NA. An adjusted mean is estimable only when its
  # coefficient-space row is invariant to that aliasing, that is, the weight
  # its aliased columns carry equals the weight the kept columns already
  # imply; otherwise lm's arbitrary zero for the aliased coefficient would
  # make the reported value an artifact. The reference grid itself is full
  # rank, so estimability is judged against the fitting design matrix.
  keep <- !is.na(beta)
  if (!all(keep)) {
    X_fit <- stats::model.matrix(model)[, names(beta), drop = FALSE]
    aliasing <- qr.solve(X_fit[, keep, drop = FALSE],
                         X_fit[, !keep, drop = FALSE])
    resid_estim <- L[, !keep, drop = FALSE] -
      L[, keep, drop = FALSE] %*% aliasing
    bad <- apply(abs(resid_estim), 1L, max) >
      1e-8 * pmax(1, apply(abs(L), 1L, max))
    if (any(bad)) {
      n_bad <- sum(bad)
      stop("The adjusted mean", if (n_bad > 1L) "s" else "", " for ",
           paste(labels_chr[bad], collapse = "; "),
           if (n_bad > 1L) " are" else " is", " not estimable: the fitted ",
           "design is rank deficient (for example an empty factorial ",
           "cell), so the model cannot predict every cell ",
           if (is.null(by)) "of the grid" else "these means average over",
           ". No adjusted value exists; check the design for empty cells.",
           call. = FALSE)
    }
  }

  L_keep <- L[, keep, drop = FALSE]
  V <- stats::vcov(model, complete = TRUE)[keep, keep, drop = FALSE]
  est <- drop(L_keep %*% beta[keep])
  se  <- sqrt(pmax(0, rowSums((L_keep %*% V) * L_keep)))
  nu  <- stats::df.residual(model)
  if (!is.finite(nu) || nu < 1) {
    stop("'model' has no residual degrees of freedom, so no standard error ",
         "or interval can be formed.")
  }
  t_crit <- stats::qt(1 - (1 - conf_level) / 2, df = nu)

  out <- cbind(
    label_df,
    data.frame(estimate = est, se = se,
               ci_lower = est - t_crit * se,
               ci_upper = est + t_crit * se,
               row.names = NULL)
  )
  rownames(out) <- NULL
  attr(out, "df_residual") <- nu
  if (!is.null(by)) attr(out, "weights") <- weights
  .as_dmar_tbl(out, conf_level = conf_level)
}
