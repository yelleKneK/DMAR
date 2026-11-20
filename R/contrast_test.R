# Contrast tests with optional multiple-comparison adjustment.
#' Tests One or More Contrasts of Group Means in a One-Way Design
#'
#' Given a fitted one-way \code{\link[stats]{aov}} or \code{\link[stats]{lm}}
#' object and a set of contrast weights, computes for every contrast the
#' estimate \eqn{\hat{\psi} = \sum_i c_i \bar{Y}_i}, its standard error,
#' \emph{t}-statistic, degrees of freedom, two-sided \emph{p}-value, and
#' confidence interval. Supports several common multiple-comparison
#' adjustments and either equal-variance (pooled) or Welch-style unequal-variance
#' inference.
#'
#' @param object A fitted \code{\link[stats]{aov}} or \code{\link[stats]{lm}}
#'   object for a one-way design (a single grouping factor on the right-hand
#'   side of the formula).
#' @param contrasts Specification of one or more contrasts. Any of:
#'   \describe{
#'     \item{\code{"pairwise"} (default)}{All pairwise comparisons among the
#'       group means.}
#'     \item{a named list of numeric vectors}{Each vector is one contrast and
#'       its name is used as the row label.}
#'     \item{a numeric matrix}{Each row is one contrast; \code{rownames}, if
#'       present, are used as labels.}
#'     \item{a numeric vector}{Treated as a single contrast.}
#'   }
#'   Each contrast vector must have length equal to the number of groups, and
#'   the weights are typically chosen to sum to zero.
#' @param adjust Multiple-comparison adjustment. One of \code{"none"}
#'   (default), \code{"bonferroni"}, \code{"scheffe"}, \code{"tukey"} (pairwise
#'   contrasts only), or any of the sequential methods supported by
#'   \code{\link[stats]{p.adjust}} (\code{"holm"}, \code{"hochberg"},
#'   \code{"BH"}, \code{"BY"}).
#' @param conf_level Confidence level for the interval (default \code{0.95}).
#' @param var_equal Logical. If \code{TRUE} (default), uses the pooled
#'   error variance \eqn{\mathit{MS}_{\text{error}}} and the residual degrees
#'   of freedom from \code{object}. If \code{FALSE}, uses each group's own
#'   sample variance and a Welch-Satterthwaite approximate \eqn{df} per
#'   contrast.
#'
#' @return A \code{data.frame} with one row per contrast and columns
#'   \code{contrast}, \code{estimate}, \code{se}, \code{t_value}, \code{df},
#'   \code{p_value}, \code{p_adjusted}, \code{ci_lower}, and
#'   \code{ci_upper}. The
#'   adjustment, confidence level, and variance assumption are stored as
#'   \code{attr(*, "adjust")}, \code{attr(*, "conf_level")}, and
#'   \code{attr(*, "var_equal")}. With \code{adjust} \code{"bonferroni"},
#'   \code{"scheffe"}, or \code{"tukey"}, the recorded level is the
#'   simultaneous level of the family of intervals, the convention of the
#'   post hoc intervals and of \code{\link{pairwise_within}} with
#'   \code{bonferroni_ci = TRUE}; under Bonferroni the line beneath the
#'   table also names the level of each interval, \eqn{1 - \alpha / m}.
#'   The table prints through the \code{\link{dmar_tbl}} display layer and works with
#'   \code{\link[generics]{tidy}} and \code{\link[generics]{glance}} (see
#'   \code{\link{dmar_tidiers}}).
#'
#' @details
#' \strong{Test statistic.} For a contrast with weights \eqn{c_1, \ldots, c_k}
#' (\eqn{k} = number of groups), the estimate is
#' \eqn{\hat{\psi} = \sum_i c_i \bar{Y}_i}. Under equal variances, the standard
#' error is \eqn{\sqrt{\mathit{MS}_{\text{error}} \sum_i c_i^2 / n_i}} with
#' \eqn{df = N - k}; under unequal variances, the standard error is
#' \eqn{\sqrt{\sum_i c_i^2 s_i^2 / n_i}} with the Welch-Satterthwaite df,
#' \deqn{df_{\text{Welch}} = \frac{\left(\sum_i c_i^2 s_i^2 / n_i\right)^2}{\sum_i (c_i^2 s_i^2 / n_i)^2 / (n_i - 1)}.}
#' The unadjusted \emph{p}-value is two-sided based on the \emph{t} reference
#' distribution.
#'
#' \strong{Adjustments.} The \code{p_adjusted} and confidence interval critical
#' value are computed as follows.
#' \itemize{
#'   \item \code{"none"}: no adjustment; the CI uses
#'     \eqn{t_{1-\alpha/2,df}}.
#'   \item \code{"bonferroni"}: \eqn{p_{\text{adj}} = \min(1, m\, p)} for
#'     \eqn{m} contrasts, with CI based on
#'     \eqn{t_{1-\alpha/(2m),df}}.
#'   \item \code{"scheffe"}: appropriate for any contrast (or family of
#'     contrasts). \eqn{p_{\text{adj}}} comes from the upper tail of an
#'     \emph{F} reference distribution applied to \eqn{t^2 / (k-1)}, and the
#'     CI uses \eqn{\sqrt{(k-1)\, F_{1-\alpha,\,k-1,df}}}.
#'   \item \code{"tukey"}: requires every contrast to be pairwise. Uses the
#'     studentized range distribution \eqn{Q_{k, df}}, so that
#'     \eqn{p_{\text{adj}} = P(Q_{k, df} > |t|\sqrt{2})} and the CI uses
#'     \eqn{q_{1-\alpha,\,k,df} / \sqrt{2}}. Both come from base R's
#'     \code{\link[stats]{ptukey}} and \code{\link[stats]{qtukey}} where those
#'     are accurate, and from the exact distribution where they are not (at
#'     small df, for example fewer than 9 df with up to 6 groups; from 10000
#'     to \eqn{10^8} df; at confidence levels above .999; with more than 50
#'     groups, below the .75 point of the distribution (.80 with more than
#'     1000 groups); and wherever
#'     \code{qtukey} stops short of its root; with 20 groups and 2 df
#'     \code{qtukey(.99, 20, 2)} is 44.526 where the exact value is 37.943);
#'     \code{\link{cv_tukey_hsd}} describes the computation and the
#'     regions.
#'   \item \code{"holm"}, \code{"hochberg"}, \code{"BH"}, \code{"BY"}:
#'     \code{\link[stats]{p.adjust}} is applied to the unadjusted
#'     \emph{p}-values; the CI uses the unadjusted \eqn{t}-critical value
#'     because these methods do not give simultaneous CIs in closed form.
#' }
#'
#' \strong{Variance assumption with adjustments.} The Tukey and Scheffé
#' procedures assume equal variances; combining them with
#' \code{var_equal = FALSE} is at the user's risk (the resulting Type I error
#' rate is no longer guaranteed). For unequal variances, common alternatives
#' are Games-Howell (Tukey-style) and Brown-Forsythe (Scheffé-style); these
#' are not currently supported here.
#'
#' \strong{Scope.} Only one-way designs are supported in v1 (one outcome,
#' one grouping factor). Multi-way designs throw an informative error.
#'
#' @references
#' Hsu, J. C. (1996). \emph{Multiple comparisons: Theory and methods}.
#' Chapman & Hall.
#'
#' Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027).
#' \emph{Designing experiments and analyzing data: A model comparison
#' perspective} (4th ed.). Routledge.
#'
#' Scheffe, H. (1953). A method for judging all contrasts in the analysis of
#' variance. \emph{Biometrika, 40}, 87--104.
#'
#' Smith, J. E., Meyers, R. J., & Delaney, H. D. (1998). The community
#'   reinforcement approach with homeless alcohol-dependent individuals.
#'   \emph{Journal of Consulting and Clinical Psychology, 66}(3), 541--548.
#'
#' Tukey, J. W. (1953). The problem of multiple comparisons. Unpublished
#' manuscript, Princeton University.
#'
#' @examples
#' # All pairwise comparisons among the three arms of the depression_bdi
#' # treatment study.
#' fit <- aov(bdi_post ~ condition, data = depression_bdi)
#' contrast_test(fit, contrasts = "pairwise")
#'
#' # Custom contrasts with a Tukey-protected family-wise error rate.
#' contrast_test(fit, contrasts = "pairwise", adjust = "tukey")
#'
#' # A user-defined contrast: the SSRI arm vs. the average of the placebo
#' # and wait list arms. With levels ordered ssri, placebo, wait_list, the
#' # weights c(1, -0.5, -0.5) estimate that difference.
#' contrast_test(
#'   fit,
#'   contrasts = list("ssri vs non-drug arms" = c(1, -0.5, -0.5)),
#'   adjust = "scheffe"
#' )
#'
#' # Welch-style inference: the wait list variance is about twice the
#' # SSRI variance, so the pooled error term is worth questioning.
#' contrast_test(fit, contrasts = "pairwise", var_equal = FALSE)
#'
#' # Pairwise treatment comparisons in the Smith, Meyers, and Delaney
#' # (1998) drinking trial, on the normalizing log scale. Each row is
#' # one pairwise contrast of the three treatment means.
#' fit_drinks <- aov(log_drinks ~ treatment, data = drinks_trial)
#' contrast_test(fit_drinks, contrasts = "pairwise")
#'
#' # An a priori contrast: the two active CRA arms (averaged) versus
#' # standard care. With levels ordered Standard, CRA, CRA + Disulfiram,
#' # the weights c(-1, 0.5, 0.5) compare the active arms against Standard.
#' contrast_test(
#'   fit_drinks,
#'   contrasts = list("CRA arms vs Standard" = c(-1, 0.5, 0.5))
#' )
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso \code{\link[stats]{TukeyHSD}}, \code{\link[stats]{pairwise.t.test}},
#'   \code{\link[stats]{p.adjust}}, \code{\link{cv_tukey_hsd}}, and
#'   \code{\link{dmar_tidiers}} for the tidy methods
#'
#' @keywords htest design
#'
#' @family hypothesis tests
#'
#' @export
#' @import stats

contrast_test <- function(
  object,
  contrasts  = "pairwise",
  adjust     = "none",
  conf_level = 0.95,
  var_equal  = TRUE
) {
  .dmar_check_conf_level(conf_level)
  adjust <- match.arg(
    adjust,
    c("none", "bonferroni", "scheffe", "tukey",
      "holm", "hochberg", "BH", "BY")
  )
  if (!is.numeric(conf_level) || length(conf_level) != 1L ||
      conf_level <= 0 || conf_level >= 1) {
    .dmar_refuse_argument(
      "'conf_level' must be a single number strictly between 0 and 1.",
      arg = "conf_level")
  }
  if (!inherits(object, c("aov", "lm"))) {
    .dmar_refuse_argument(
      "'object' must be a fitted aov or lm object.", arg = "object")
  }

  # Extract one-way structure from the model frame.
  mf <- stats::model.frame(object)
  if (ncol(mf) != 2L) {
    .dmar_refuse_argument(c(
      "contrast_test() currently supports one-way designs (one outcome ~ one grouping factor); ",
      "the supplied model has ", ncol(mf) - 1L, " predictors."), arg = NULL)
  }
  y_vec <- mf[[1]]
  g_vec <- mf[[2]]
  if (!is.factor(g_vec)) g_vec <- factor(g_vec)
  g_vec <- droplevels(g_vec)

  levs   <- levels(g_vec)
  k      <- length(levs)
  if (k < 2L) .dmar_refuse_argument(
    "At least two groups are required.", arg = NULL)
  n_per  <- as.integer(table(g_vec))
  means  <- as.numeric(tapply(y_vec, g_vec, mean))
  vars_  <- as.numeric(tapply(y_vec, g_vec, stats::var))

  cm <- .resolve_contrasts(contrasts, levs)
  if (ncol(cm) != k) {
    .dmar_refuse_argument(c(
      "Each contrast must have length equal to the number of groups (",
      k, ")."), arg = NULL)
  }
  # A contrast whose weights are all zero defines none; it gave the
  # interval [0, 0] with a t value and p value of NaN. A set that holds no
  # contrast (a matrix with no row, an empty list) is refused here too.
  .dmar_check_contrast_nonzero(cm, arg = "contrasts", margin = 1L)

  # Warn (don't error) if a contrast row does not sum to zero.
  row_sums <- rowSums(cm)
  if (any(abs(row_sums) > 1e-8)) {
    bad <- rownames(cm)[abs(row_sums) > 1e-8]
    .dmar_warn(c(
      "Contrast(s) [", paste(bad, collapse = ", "),
      "] do not sum to zero; results are still computed but interpret with care."),
      "result_qualified")
  }

  # Pooled-variance ingredients (used when var_equal = TRUE).
  # A group with a single observation has no sample variance (var() gives NA)
  # and adds nothing to the within-group sum of squares, so it is left out of
  # the sum rather than turning the pooled variance into NA.
  ms_error    <- sum(((n_per - 1) * vars_)[n_per > 1]) / (sum(n_per) - k)
  df_residual <- sum(n_per) - k
  m           <- nrow(cm)

  # Per-contrast estimate, SE, t, df.
  estimate <- as.numeric(cm %*% means)
  if (isTRUE(var_equal)) {
    se <- sqrt(ms_error * rowSums(sweep(cm^2, 2L, n_per, "/")))
    df <- rep(df_residual, m)
  } else {
    # term_ij = c_ij^2 * s_j^2 / n_j
    terms <- sweep(cm^2, 2L, vars_ / n_per, "*")
    se <- sqrt(rowSums(terms))
    df <- rowSums(terms)^2 / rowSums(sweep(terms^2, 2L, n_per - 1, "/"))
  }
  t_stat  <- estimate / se
  p_unadj <- 2 * stats::pt(-abs(t_stat), df = df)

  # Adjustment-specific p_adj and CI critical values.
  alpha <- 1 - conf_level
  if (adjust == "none") {
    p_adj  <- p_unadj
    t_crit <- stats::qt(1 - alpha / 2, df = df)
  } else if (adjust == "bonferroni") {
    p_adj  <- pmin(1, p_unadj * m)
    t_crit <- stats::qt(1 - alpha / (2 * m), df = df)
  } else if (adjust == "scheffe") {
    F_obs  <- t_stat^2 / (k - 1)
    p_adj  <- stats::pf(F_obs, df1 = k - 1, df2 = df, lower.tail = FALSE)
    t_crit <- sqrt((k - 1) * stats::qf(1 - alpha, df1 = k - 1, df2 = df))
  } else if (adjust == "tukey") {
    pairwise_ok <- apply(cm, 1L, .is_pairwise_contrast)
    if (!all(pairwise_ok)) {
      .dmar_refuse_argument(c(
        "'adjust = \"tukey\"' requires every contrast to be pairwise ",
        "(exactly one +c and one -c, with the rest zero). Failing rows: ",
        paste(rownames(cm)[!pairwise_ok], collapse = ", "), "."), arg = NULL)
    }
    # stats::ptukey()/qtukey() where they are accurate, the exact distribution
    # elsewhere (R/studentized_range_internals.R).
    p_adj  <- 1 - .sr_cdf(abs(t_stat) * sqrt(2), k, df)
    t_crit <- .sr_quantile(1 - alpha, k, df) / sqrt(2)
  } else {
    # Sequential methods via p.adjust(); no closed-form simultaneous CI, so
    # the CI uses the unadjusted t critical value.
    p_adj  <- stats::p.adjust(p_unadj, method = adjust)
    t_crit <- stats::qt(1 - alpha / 2, df = df)
  }

  out <- data.frame(
    contrast   = rownames(cm),
    estimate   = estimate,
    se         = se,
    t_value    = t_stat,
    df         = df,
    p_value    = p_unadj,
    p_adjusted = p_adj,
    ci_lower   = estimate - t_crit * se,
    ci_upper   = estimate + t_crit * se,
    stringsAsFactors = FALSE,
    row.names  = NULL
  )

  attr(out, "adjust")     <- adjust
  attr(out, "p_adjustment") <- adjust
  attr(out, "ci_adjustment") <- if (adjust %in% c("bonferroni", "tukey", "scheffe")) adjust else "none"
  attr(out, "conf_level") <- conf_level
  attr(out, "var_equal")  <- var_equal
  # The Bonferroni intervals hold jointly at conf_level, the level the
  # table records; the line beneath the table names the level of each
  # interval, as pairwise_within(bonferroni_ci = TRUE) does.
  if (adjust == "bonferroni" && m > 1L) {
    attr(out, "simultaneous") <- list(method = "Bonferroni",
                                      n_comparisons = m,
                                      each = 1 - alpha / m)
  }

  # Route through the dmar_tbl display layer (idempotent), then layer the
  # broom-dispatch subclass ahead of dmar_tbl so print() still falls through
  # to print.dmar_tbl (see R/dmar_tidiers.R).
  if (!inherits(out, "dmar_tbl")) {
    out <- .as_dmar_tbl(out, conf_level = conf_level)
  }
  class(out) <- c("dmar_contrast_test", class(out))
  out
}


# The tidy() and glance() methods for contrast_test(). They are documented
# with the other DMAR tidiers at ?dmar_tidiers (R/dmar_tidiers.R), which
# states the contract these methods share with the confidence interval, post
# hoc, and sample size planning families.

#' @rdname dmar_tidiers
#' @importFrom generics tidy
#' @exportS3Method generics::tidy
tidy.dmar_contrast_test <- function(x, ...) {
  out <- data.frame(
    term       = x$contrast,
    estimate   = x$estimate,
    ci_lower   = x$ci_lower,
    ci_upper   = x$ci_upper,
    statistic  = x$t_value,
    df         = x$df,
    p_value    = x$p_value,
    p_adjusted = x$p_adjusted,
    conf_level = attr(x, "conf_level") %||% NA_real_,
    stringsAsFactors = FALSE,
    row.names  = NULL
  )
  .dmar_tidy_adjustments(out, x)
}

#' @rdname dmar_tidiers
#' @importFrom generics glance
#' @exportS3Method generics::glance
glance.dmar_contrast_test <- function(x, ...) {
  out <- data.frame(
    n_contrasts     = nrow(x),
    adjust          = attr(x, "adjust")    %||% NA_character_,
    var_equal       = attr(x, "var_equal") %||% NA,
    p_adjusted_min  = if (nrow(x)) min(x$p_adjusted) else NA_real_,
    conf_level      = attr(x, "conf_level") %||% NA_real_,
    stringsAsFactors = FALSE,
    row.names       = NULL
  )
  .dmar_tidy_adjustments(out, x)
}


# ----- Internal helpers -----------------------------------------------------

# Resolve the user's 'contrasts' argument into a numeric matrix with row names.
.resolve_contrasts <- function(contrasts, levs) {
  k <- length(levs)

  # "pairwise": all C(k, 2) pairs as level1 - level2.
  if (is.character(contrasts) && length(contrasts) == 1L &&
      identical(contrasts, "pairwise")) {
    pairs <- utils::combn(k, 2L)
    cm <- matrix(0, nrow = ncol(pairs), ncol = k)
    rn <- character(ncol(pairs))
    # Convention matches stats::TukeyHSD(): later level minus earlier level.
    for (i in seq_len(ncol(pairs))) {
      cm[i, pairs[2L, i]] <-  1
      cm[i, pairs[1L, i]] <- -1
      rn[i] <- paste(levs[pairs[2L, i]], "-", levs[pairs[1L, i]])
    }
    rownames(cm) <- rn
    return(cm)
  }

  # A list or matrix that holds no contrast becomes a matrix with no row,
  # which the check of the weights in contrast_test() refuses; naming its
  # rows stopped with a message from R before.
  if (is.list(contrasts)) {
    cm <- do.call(rbind, contrasts)
    if (is.null(cm)) cm <- matrix(numeric(0), nrow = 0L, ncol = k)
    if (is.null(rownames(cm)) && !is.null(names(contrasts))) {
      rownames(cm) <- names(contrasts)
    }
    if (nrow(cm) > 0L && (is.null(rownames(cm)) || any(rownames(cm) == ""))) {
      rownames(cm) <- paste0("C", seq_len(nrow(cm)))
    }
    return(cm)
  }

  if (is.matrix(contrasts)) {
    cm <- contrasts
    if (nrow(cm) > 0L && is.null(rownames(cm))) {
      rownames(cm) <- paste0("C", seq_len(nrow(cm)))
    }
    return(cm)
  }

  if (is.numeric(contrasts) && is.null(dim(contrasts))) {
    cm <- matrix(contrasts, nrow = 1L)
    rownames(cm) <- "C1"
    return(cm)
  }

  .dmar_refuse_argument(
    "'contrasts' must be \"pairwise\", a named list, a matrix, or a numeric vector.",
    arg = "contrasts")
}

# A "pairwise" contrast has exactly one positive weight, one negative weight of
# equal magnitude, and the rest zero. Tolerance allows for floating-point dust.
.is_pairwise_contrast <- function(c_i, tol = 1e-8) {
  nz <- abs(c_i) > tol
  if (sum(nz) != 2L) return(FALSE)
  vals <- c_i[nz]
  isTRUE(abs(vals[1L] + vals[2L]) < tol) &&
    isTRUE(abs(abs(vals[1L]) - abs(vals[2L])) < tol)
}
