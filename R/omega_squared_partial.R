# Partial omega squared (effect size for ANOVA).
#' Partial Omega Squared (Effect Size for ANOVA)
#'
#' Computes the sample \emph{partial} omega squared (\eqn{\omega^2_p}), Hays'
#' (1994) bias-corrected estimator of the proportion of population variance in
#' the dependent variable accounted for by a fixed effect after the variance
#' attributable to the other effects in the model has been removed:
#' \deqn{\hat{\omega}^2_p \;=\;
#'   \frac{\mathit{SS}_{\text{effect}} - df_{\text{effect}} \cdot \mathit{MS}_{\text{error}}}
#'        {\mathit{SS}_{\text{effect}} + (N - df_{\text{effect}}) \cdot \mathit{MS}_{\text{error}}}
#'   \;=\;
#'   \frac{df_{\text{effect}} (F - 1)}{df_{\text{effect}} (F - 1) + N}.}
#' Accepts either the raw ANOVA summary (\emph{F}, effect df, error df, total
#' \emph{N}) or a fitted \code{aov}/\code{lm}/\code{aovlist} object, in which
#' case the function returns one row per effect (with stratum identification
#' for within-subjects fits).
#'
#' This function is the explicitly-named counterpart of
#' \code{\link{omega_squared}}. The two share the same point-estimate formula:
#' in a one-way ANOVA they coincide with the total \eqn{\omega^2}; in a
#' factorial ANOVA both return the per-effect \emph{partial} value computed
#' against the model's residual mean square. \code{omega_squared_partial} is
#' provided so that user code that explicitly intends partial \eqn{\omega^2}
#' carries that meaning in its name, parallel to the
#' \code{\link{eta_squared}} / \code{\link{eta_squared_partial}} pair.
#'
#' @param object Optional. A fitted model object of class
#'   \code{\link[stats]{aov}}, \code{\link[stats]{lm}}, or
#'   \code{aovlist} (multi-stratum aov fit, e.g.\
#'   \code{aov(y ~ A + Error(subject/A), data = d)}). When supplied, the
#'   function loops over the non-\code{Residuals} rows and returns one row
#'   per effect.
#' @param F_value Observed \emph{F}-value from the fixed-effects ANOVA
#'   (ignored if \code{object} is supplied).
#' @param df_effect Numerator degrees of freedom for the effect
#'   (ignored if \code{object} is supplied).
#' @param df_error Error (residual) degrees of freedom
#'   (ignored if \code{object} is supplied).
#' @param N Total sample size (ignored if \code{object} is supplied;
#'   \code{\link[stats]{nobs}(object)} is used instead).
#'
#' @return A \code{data.frame} with one row per effect. The columns
#'   are \code{effect}, \code{omega_squared_partial} (point estimate),
#'   \code{F_value}, \code{df_effect}, \code{df_error}, and \code{N}.
#'   When the raw-argument interface is used, \code{effect} is
#'   \code{"overall"}. Negative point estimates (which occur whenever
#'   \emph{F} < 1) are truncated to zero, matching the convention used by
#'   \code{\link{omega_squared}} and \code{\link{ci_omega_squared}}.
#'
#' @details
#' \strong{Why partial omega squared and not total.} In a one-way ANOVA,
#' partial \eqn{\omega^2} reduces to total \eqn{\omega^2}; in a factorial
#' ANOVA the two diverge. Total \eqn{\omega^2} for an effect divides its
#' variance contribution by the \emph{total} population variance of \eqn{Y},
#' so adding orthogonal factors to a study mechanically shrinks each effect's
#' total \eqn{\omega^2}. Partial \eqn{\omega^2} divides instead by the
#' variance that is left after the other effects in the model have been
#' partialled out, so a given fixed effect's partial \eqn{\omega^2} is
#' (approximately) invariant to whether additional orthogonal factors are
#' present (Olejnik & Algina, 2003; Maxwell, Delaney, & Kelley, 2027,
#' Sections 7.4.4 and 8.4). For that reason, partial \eqn{\omega^2} is the
#' chapter's preferred effect size index when off-factors are "extrinsic"
#' (i.e., would not vary in a hypothetical full replication of the population
#' setup).
#'
#' \strong{Bias correction vs.\ partial eta squared.} \eqn{\hat{\eta}^2_p},
#' the sample partial eta squared, is the proportion of \emph{sample}
#' variance accounted for and is upward-biased as an estimator of the
#' population \eqn{\eta^2_p}. \eqn{\hat{\omega}^2_p} subtracts
#' \eqn{df_{\text{effect}} \cdot \mathit{MS}_{\text{error}}} from the
#' effect's sum of squares and rescales, yielding an estimator of the
#' \emph{population} variance proportion with substantially smaller bias
#' (Hays, 1994; Olejnik & Algina, 2000; Kelley, 2007). Truncation at zero
#' is conventional when the unbiased estimator goes negative because
#' \eqn{\omega^2 \ge 0} by definition.
#'
#' \strong{Hand-in-hand with \code{ci_omega_squared()}.} Pair this function
#' with \code{\link{ci_omega_squared}} when reporting effect sizes:
#' \code{omega_squared_partial()} returns the point estimate(s) and
#' \code{ci_omega_squared()} returns the same point estimate plus its
#' noncentral \emph{F} confidence limits (Steiger, 2004; Kelley, 2007). The
#' columns shared by the two functions are aligned so the outputs compose
#' cleanly with \code{merge()} or a join.
#'
#' \strong{Sums of squares in unbalanced factorial designs.}
#' \code{\link[stats]{anova}()} on an \code{aov}/\code{lm} uses Type I
#' (sequential) sums of squares. For balanced designs all three SS types
#' agree; for unbalanced designs they differ. If Type II or III
#' \emph{F}-values are required, compute them with e.g.\
#' \code{car::Anova(object, type = 3)} and pass the relevant \emph{F} and
#' degrees of freedom into the raw-argument interface.
#'
#' @references
#' Cohen, J. (1973). Eta-squared and partial eta-squared in fixed factor
#'   ANOVA designs. \emph{Educational and Psychological Measurement, 33}(1),
#'   107--112.
#'
#' Hays, W. L. (1994). \emph{Statistics} (5th ed.). Fort Worth, TX:
#'   Harcourt Brace College Publishers.
#'
#' Kelley, K. (2007). Confidence intervals for standardized
#'   effect sizes: Theory, application, and implementation.
#'   \emph{Journal of Statistical Software, 20}(8), 1--24.
#'   \doi{10.18637/jss.v020.i08}
#'
#' Kelley, K., & Preacher, K. J. (2012). On effect size.
#'   \emph{Psychological Methods, 17}, 137--152. \doi{10.1037/a0028086}
#'
#' Keppel, G., & Wickens, T. D. (2004). \emph{Design and analysis: A
#'   researcher's handbook} (4th ed.). Pearson Prentice Hall.
#'
#' Keren, G., & Lewis, C. (1979). Partial omega squared for ANOVA designs.
#'   \emph{Educational and Psychological Measurement, 39}(1), 119--128.
#'
#' Maxwell, S. E., Camp, C. J., & Arvey, R. D. (1981). Measures of strength
#'   of association: A comparative examination. \emph{Journal of Applied
#'   Psychology, 66}(5), 525--534.
#'
#' Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). \emph{Designing
#'   experiments and analyzing data: A model comparison perspective}
#'   (4th ed.). Routledge. (See Chapter 3 on \eqn{\eta^2}, Chapter 7 on
#'   factorial designs, and Chapter 11 on generalized \eqn{\eta^2} for
#'   within-subjects designs.)
#'
#' Olejnik, S., & Algina, J. (2000). Measures of effect size for comparative
#'   studies: Applications, interpretations, and limitations.
#'   \emph{Contemporary Educational Psychology, 25}(3), 241--286.
#'   \doi{10.1006/ceps.2000.1040}
#'
#' Olejnik, S., & Algina, J. (2003). Generalized eta and omega squared
#'   statistics: Measures of effect size for some common research designs.
#'   \emph{Psychological Methods, 8}(4), 434--447.
#'   \doi{10.1037/1082-989X.8.4.434}
#'
#' Steiger, J. H. (2004). Beyond the \emph{F} test: Effect size confidence
#'   intervals and tests of close fit in the analysis of variance and
#'   contrast analysis. \emph{Psychological Methods, 9}(2), 164--182.
#'   \doi{10.1037/1082-989X.9.2.164}
#'
#' @seealso \code{\link{omega_squared}}, \code{\link{ci_omega_squared}},
#'   \code{\link{eta_squared_partial}}, \code{\link{ci_eta_squared_partial}}
#'
#' @examples
#' # 1. Raw-argument interface (Bargman 1970 / Steiger 2004 example):
#' #        five groups of 11, observed F = 11.221.
#' omega_squared_partial(F_value = 11.221, df_effect = 4, df_error = 50, N = 55)
#'
#' # 2. Two-factor ANOVA: partial omega squared per effect on the
#' #        pygmalion data (expectancy treatment x grade, unequal cell
#' #        sizes, N = 310). The treatment is manipulated while grade is
#' #        a measured classification, and the partial value for each
#' #        effect removes the variance the other accounts for. The
#' #        treatment by grade interaction is weak here (F = 1.19), so
#' #        the additive model is used.
#' fit_additive <- aov(iq_8 ~ treatment + factor(grade), data = pygmalion)
#' omega_squared_partial(fit_additive)
#'
#' # 3. omega_squared_partial() and ci_omega_squared() agree on the
#' #        point estimate row-by-row.
#' pt  <- omega_squared_partial(fit_additive)
#' ci  <- ci_omega_squared(fit_additive)
#' pt$omega_squared_partial
#' ci$omega_squared
#'
#' # 4. The named pair: omega_squared() and omega_squared_partial()
#' #        report identical numbers in this design; the only difference
#' #        is the name of the value column, which makes the user's
#' #        intent (partial) explicit.
#' omega_squared(fit_additive)$omega_squared
#' omega_squared_partial(fit_additive)$omega_squared_partial
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @keywords htest design
#'
#' @family effect size estimates
#'
#' @export
#' @import stats

omega_squared_partial <- function(
  object    = NULL,
  F_value   = NULL,
  df_effect = NULL,
  df_error  = NULL,
  N         = NULL
) {
  if (!is.null(object)) {
    return(.omega_squared_partial_from_model(object))
  }
  if (is.null(F_value) || is.null(df_effect) || is.null(df_error) || is.null(N)) {
    stop("Either provide 'object' (a fitted aov or lm) or all of 'F_value', 'df_effect', 'df_error', and 'N'.")
  }
  .omega_squared_partial_one(F_value, df_effect, df_error, N,
                             effect_label = "overall")
}


# Single-effect worker. Returns a 1-row data.frame with the
# Hays (1994) point estimate and the inputs echoed back so the row
# is self-describing and joins cleanly with ci_omega_squared() output.
.omega_squared_partial_one <- function(F_value, df_effect, df_error, N,
                                       effect_label) {
  if (!is.numeric(F_value) || length(F_value) != 1L || F_value <= 0)
    stop("'F_value' must be a single positive number.")
  if (!is.numeric(df_effect) || length(df_effect) != 1L || df_effect <= 0)
    stop("'df_effect' must be a single positive number.")
  if (!is.numeric(df_error) || length(df_error) != 1L || df_error <= 0)
    stop("'df_error' must be a single positive number.")
  if (!is.numeric(N) || length(N) != 1L || N <= df_effect + df_error)
    stop("'N' must be greater than df_effect + df_error.")

  num <- df_effect * (F_value - 1)
  omega_sq_p <- max(0, num / (num + N))

  data.frame(
    effect                = effect_label,
    omega_squared_partial = omega_sq_p,
    F_value               = F_value,
    df_effect             = df_effect,
    df_error              = df_error,
    N                     = N,
    stringsAsFactors      = FALSE,
    row.names             = NULL
  )
}


# Model-based worker: loops over the non-Residuals rows of anova(object).
# Dispatches between flat lm/aov fits and within-subjects aovlist fits
# (returned by stats::aov() with an Error() term).
.omega_squared_partial_from_model <- function(object) {
  if (inherits(object, "aovlist"))
    return(.omega_squared_partial_from_aovlist(object))
  if (!inherits(object, c("aov", "lm"))) {
    stop("'object' must be an aov, lm, or aovlist fit.")
  }
  tbl <- stats::anova(object)
  resid_row <- which(rownames(tbl) == "Residuals")
  if (length(resid_row) == 0L) {
    stop("Could not find a 'Residuals' row in anova(object); is this a standard fixed-effects ANOVA?")
  }
  df_error <- tbl[resid_row, "Df"]
  N <- stats::nobs(object)
  effect_rows <- setdiff(rownames(tbl), "Residuals")

  results <- lapply(effect_rows, function(effect) {
    F_val <- tbl[effect, "F value"]
    df_e  <- tbl[effect, "Df"]
    if (is.na(F_val) || F_val <= 0) return(NULL)
    .omega_squared_partial_one(F_val, df_e, df_error, N,
                               effect_label = effect)
  })
  results <- Filter(Negate(is.null), results)
  if (length(results) == 0L) {
    stop("No effects with a computable F statistic were found in the model.")
  }

  out <- do.call(rbind, results)
  rownames(out) <- NULL
  out
}

# Worker for aovlist (within-subjects / Error() fits): one row per
# effect, using each effect's own stratum residual as its denominator.
.omega_squared_partial_from_aovlist <- function(object) {
  smry <- summary(object)
  # Total N: sum the strata degrees of freedom plus the intercept
  # (the shared aovlist helper defined in eta_squared.R).
  N <- .aovlist_nobs(object)

  rows <- list()
  for (stratum_nm in names(smry)) {
    tbl <- smry[[stratum_nm]][[1L]]
    if (is.null(tbl)) next
    rn  <- trimws(rownames(tbl))
    resid_idx <- which(rn == "Residuals")
    if (length(resid_idx) != 1L) next
    df_error <- tbl[resid_idx, "Df"]
    if (is.na(df_error) || df_error <= 0) next
    for (i in seq_len(nrow(tbl))) {
      if (rn[i] == "Residuals") next
      F_val <- tbl[i, "F value"]
      df_e  <- tbl[i, "Df"]
      if (is.na(F_val) || F_val <= 0) next
      rows[[length(rows) + 1L]] <-
        .omega_squared_partial_one(F_val, df_e, df_error, N,
                                   effect_label = rn[i])
    }
  }
  if (length(rows) == 0L)
    stop("No testable effects found in aovlist object.")
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}
