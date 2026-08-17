#' Omega Squared (Effect Size for ANOVA)
#'
#' Computes the sample omega squared (\eqn{\omega^2}), Hays' (1994)
#' bias-corrected estimator of the proportion of variance in the dependent
#' variable accounted for by a fixed effect. Accepts either the raw ANOVA
#' summary (\emph{F}, the effect and error degrees of freedom, and total
#' \emph{N}) or a fitted \code{\link[stats]{aov}} or \code{\link[stats]{lm}}
#' object, in which case the function returns one row per effect (partial
#' \eqn{\omega^2} in factorial designs).
#'
#' The confidence interval is provided by the separate
#' \code{\link{ci_omega_squared}}, paralleling the existing
#' \code{\link{smd}}/\code{\link{ci_smd}} and
#' \code{\link{eta_squared}}/\code{\link{ci_eta_squared}} pairings.
#'
#' @param object Optional. A fitted \code{\link[stats]{aov}} or
#'   \code{\link[stats]{lm}} object. When supplied, the function loops
#'   over the non-\code{Residuals} rows of \code{\link[stats]{anova}(object)}
#'   and returns one row per effect.
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
#'   are \code{effect}, \code{omega_squared} (point estimate),
#'   \code{F_value}, \code{df_effect}, \code{df_error}, and \code{N}.
#'   When the raw-argument interface is used, \code{effect} is
#'   \code{"overall"}.
#'
#' @details
#' \strong{Point estimate.} The reported value is Hays' (1994) sample
#' omega squared, which for a one-way design is
#' \deqn{\hat{\omega}^2 = \frac{df_{\text{effect}} (F - 1)}{df_{\text{effect}} (F - 1) + N}.}
#' For factorial designs the same formula applied per effect yields
#' \emph{partial} omega squared (Olejnik & Algina, 2003); negative values
#' are truncated to zero. This is the same point-estimate convention used
#' by \code{\link{ci_omega_squared}}, so the two functions agree on the
#' point estimate row by row.
#'
#' \strong{Hand-in-hand with \code{ci_omega_squared()}.} Pair this
#' function with \code{\link{ci_omega_squared}} when reporting effect
#' sizes: \code{omega_squared()} returns the point estimate(s), and
#' \code{ci_omega_squared()} returns the same point estimate plus its
#' noncentrality-based confidence limits (Steiger, 2004). The columns
#' shared by the two functions (\code{effect}, \code{omega_squared},
#' \code{F_value}, \code{df_effect}, \code{df_error}, \code{N}) are
#' aligned so the outputs compose cleanly with \code{merge()} or a join.
#'
#' \strong{Sums of squares in factorial designs.} \code{anova()} on an
#' \code{aov}/\code{lm} uses Type I (sequential) sums of squares. For
#' balanced designs all three types agree; for unbalanced designs they
#' differ. If Type II or III \emph{F}-values are required, compute them
#' with e.g.\ \code{car::Anova(object, type = 3)} and pass the relevant
#' \emph{F} and degrees of freedom into the raw-argument interface.
#'
#' @references
#' Hays, W. L. (1994). \emph{Statistics} (5th ed.). Fort Worth, TX:
#'   Harcourt Brace College Publishers.
#'
#' Kelley, K. (2007). Confidence intervals for standardized effect
#'   sizes: Theory, application, and implementation. \emph{Journal of Statistical
#'   Software, 20}(8), 1--24. \doi{10.18637/jss.v020.i08}
#'
#' Kelley, K., & Preacher, K. J. (2012). On effect size.
#'   \emph{Psychological Methods, 17}, 137--152. \doi{10.1037/a0028086}
#'
#' Keppel, G. (1991). \emph{Design and analysis: A researcher's handbook}
#'   (3rd ed.). Prentice Hall.
#'
#' Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). \emph{Designing
#'   experiments and analyzing data: A model comparison perspective}
#'   (4th ed.). Routledge. (See Chapter 3 on \eqn{\eta^2}, Chapter 7 on
#'   factorial designs, and Chapter 11 on generalized \eqn{\eta^2} for
#'   within-subjects designs.)
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
#' @seealso \code{\link{ci_omega_squared}}, \code{\link{eta_squared}},
#'   \code{\link{ci_eta_squared}}, \code{\link{ci_pvaf}}
#'
#' @examples
#' # 1. Raw-argument interface. Bargman's (1970) 5-group one-way ANOVA,
#' #        also used in Venables (1975), Fleishman (1980), and Steiger (2004):
#' #        11 subjects per group, observed F = 11.221.
#' omega_squared(F_value = 11.221, df_effect = 4, df_error = 50, N = 55)
#'
#' # 2. One way ANOVA from a fitted model (depression_bdi: three
#' #        treatment arms, 10 per arm, N = 30).
#' fit_one <- aov(bdi_post ~ condition, data = depression_bdi)
#' omega_squared(fit_one)
#'
#' # 3. Two-factor ANOVA: partial omega squared per effect for the
#' #        manipulated expectancy treatment and the measured grade
#' #        classification (pygmalion data, unequal cell sizes,
#' #        N = 310). The treatment by grade interaction is weak here
#' #        (F = 1.19), so the additive model is used.
#' fit_additive <- aov(iq_8 ~ treatment + factor(grade), data = pygmalion)
#' omega_squared(fit_additive)
#'
#' # 4. omega_squared() and ci_omega_squared() compose: the point
#' #        estimates agree row-by-row.
#' pt  <- omega_squared(fit_additive)
#' ci  <- ci_omega_squared(fit_additive)
#' merge(pt, ci, by = c("effect", "omega_squared",
#'                      "F_value", "df_effect", "df_error", "N"))
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @keywords htest design
#'
#' @family effect size estimates
#'
#' @export
#' @import stats

omega_squared <- function(
  object    = NULL,
  F_value   = NULL,
  df_effect = NULL,
  df_error  = NULL,
  N         = NULL
) {
  if (!is.null(object)) {
    return(.omega_squared_from_model(object))
  }
  if (is.null(F_value) || is.null(df_effect) || is.null(df_error) || is.null(N)) {
    stop("Either provide 'object' (a fitted aov or lm) or all of 'F_value', 'df_effect', 'df_error', and 'N'.")
  }
  .omega_squared_one(F_value, df_effect, df_error, N, effect_label = "overall")
}


# Single-effect worker. Returns a 1-row data.frame with the
# Hays (1994) point estimate and the inputs echoed back so the row
# is self-describing and joins cleanly with ci_omega_squared() output.
.omega_squared_one <- function(F_value, df_effect, df_error, N, effect_label) {
  if (!is.numeric(F_value) || length(F_value) != 1L || F_value <= 0)
    stop("'F_value' must be a single positive number.")
  if (!is.numeric(df_effect) || length(df_effect) != 1L || df_effect <= 0)
    stop("'df_effect' must be a single positive number.")
  if (!is.numeric(df_error) || length(df_error) != 1L || df_error <= 0)
    stop("'df_error' must be a single positive number.")
  if (!is.numeric(N) || length(N) != 1L || N <= df_effect + df_error)
    stop("'N' must be greater than df_effect + df_error.")

  num <- df_effect * (F_value - 1)
  omega_sq <- max(0, num / (num + N))

  data.frame(
    effect        = effect_label,
    omega_squared = omega_sq,
    F_value       = F_value,
    df_effect     = df_effect,
    df_error      = df_error,
    N             = N,
    stringsAsFactors = FALSE,
    row.names     = NULL
  )
}


# Model-based worker: loops over the non-Residuals rows of anova(object),
# matching the structure used by ci_omega_squared() so the two functions
# return aligned columns.
.omega_squared_from_model <- function(object) {
  if (!inherits(object, c("aov", "lm"))) {
    stop("'object' must be an aov or lm fit.")
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
    .omega_squared_one(F_val, df_e, df_error, N, effect_label = effect)
  })
  results <- Filter(Negate(is.null), results)
  if (length(results) == 0L) {
    stop("No effects with a computable F statistic were found in the model.")
  }

  out <- do.call(rbind, results)
  rownames(out) <- NULL
  out
}
