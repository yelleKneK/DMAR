# Partial eta squared (effect size for ANOVA).
#' Partial Eta Squared (Effect Size for ANOVA)
#'
#' Computes the sample \emph{partial} eta squared (\eqn{\eta^2_p}), the
#' proportion of variance accounted for by a fixed effect after the variance
#' attributable to the other effects in the model has been removed:
#' \deqn{\hat{\eta}^2_p = \frac{\mathit{SS}_{\text{effect}}}{\mathit{SS}_{\text{effect}} +
#'   \mathit{SS}_{\text{error}}} = \frac{df_{\text{effect}} \cdot F}{df_{\text{effect}} \cdot F + df_{\text{error}}}.}
#' Accepts either the raw ANOVA summary (\emph{F}, effect df, error df) or a
#' fitted \code{aov}/\code{lm}/\code{aovlist} object, in which case the
#' function returns one row per effect (with stratum identification for
#' within-subjects fits).
#'
#' This function is the explicitly-named counterpart of \code{\link{eta_squared}}.
#' The two share the same point-estimate formula, in a one-way ANOVA they
#' coincide with \emph{total} \eqn{\eta^2}; in a factorial or within-subjects
#' ANOVA both functions return the per-effect \emph{partial} value computed
#' against that effect's own error stratum. \code{eta_squared_partial} is
#' provided so that user code that explicitly intends partial \eqn{\eta^2}
#' carries that meaning in its name.
#'
#' @param object Optional. A fitted model object of class
#'   \code{\link[stats]{aov}}, \code{\link[stats]{lm}}, or
#'   \code{aovlist} (multi-stratum aov fit, e.g.\
#'   \code{aov(y ~ A + Error(subject/A), data = d)}).
#' @param F_value Observed \emph{F}-value (ignored if \code{object} is supplied).
#' @param df_effect Numerator degrees of freedom for the effect
#'   (ignored if \code{object} is supplied).
#' @param df_error Error (residual) degrees of freedom
#'   (ignored if \code{object} is supplied).
#'
#' @return A \code{data.frame} with one row per effect. Single-stratum
#'   fits and the raw interface return columns \code{effect},
#'   \code{eta_squared_partial}, \code{F_value}, \code{df_effect},
#'   \code{df_error}. \code{aovlist} (within-subjects / mixed) fits
#'   additionally include a \code{stratum} column identifying which error
#'   term each effect's \emph{F} test came from. With the raw-argument
#'   interface \code{effect} is \code{"overall"}.
#'
#' @details
#' \strong{Designs supported.} Single-stratum \code{aov}/\code{lm} fits and
#' multi-stratum \code{aovlist} fits (within-subjects and mixed designs)
#' are both handled by the model interface. For multi-stratum fits, each
#' effect uses its own stratum's residual \emph{df}, so a within-subjects
#' factor's partial \eqn{\eta^2} is computed against the within-subjects
#' error and a between-subjects factor's is computed against the
#' between-subjects error.
#'
#' @references
#' Cohen, J. (1973). Eta-squared and partial eta-squared in fixed factor ANOVA
#'   designs. \emph{Educational and Psychological Measurement, 33}(1), 107--112.
#'
#' Kelley, K. (2007). Confidence intervals for standardized effect
#'   sizes: Theory, application, and implementation. \emph{Journal of Statistical
#'   Software, 20}(8), 1--24. \doi{10.18637/jss.v020.i08}
#'
#' Kelley, K., & Preacher, K. J. (2012). On effect size.
#'   \emph{Psychological Methods, 17}, 137--152. \doi{10.1037/a0028086}
#'
#' Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). \emph{Designing
#'   experiments and analyzing data: A model comparison perspective}
#'   (4th ed.). Routledge. (See Chapter 3 on \eqn{\eta^2}, Chapter 7 on
#'   factorial designs, and Chapter 11 on generalized \eqn{\eta^2} for
#'   within-subjects designs.)
#'
#' Steiger, J. H. (2004). Beyond the \emph{F} test: Effect size confidence
#'   intervals and tests of close fit in the analysis of variance and contrast
#'   analysis. \emph{Psychological Methods, 9}(2), 164--182.
#'   \doi{10.1037/1082-989X.9.2.164}
#'
#' @seealso \code{\link{ci_eta_squared_partial}}, \code{\link{eta_squared}}
#'
#' @examples
#' # Raw-argument interface.
#' eta_squared_partial(F_value = 11.221, df_effect = 4, df_error = 50)
#'
#' # Factorial ANOVA: partial eta squared per effect (pygmalion data:
#' # expectancy treatment x grade, 2 x 6 with unequal cell sizes,
#' # N = 310). The treatment is manipulated; grade is a measured
#' # classification of the pupils.
#' fit <- aov(iq_8 ~ treatment * factor(grade), data = pygmalion)
#' eta_squared_partial(fit)
#'
#' # Within-subjects ANOVA: per-effect partial eta squared with stratum.
#' set.seed(113)
#' n <- 20
#' rm_data <- data.frame(
#'   subject = factor(rep(seq_len(n), each = 3)),
#'   time    = factor(rep(c("Pre", "Mid", "Post"), n),
#'                    levels = c("Pre", "Mid", "Post")),
#'   y       = rnorm(n, sd = 1.5)[rep(seq_len(n), each = 3)] +
#'             0.7 * rep(1:3, n) + rnorm(n * 3, sd = 1.2)
#' )
#' fit_rm <- aov(y ~ time + Error(subject/time), data = rm_data)
#' eta_squared_partial(fit_rm)
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @keywords htest design
#'
#' @family effect size estimates
#'
#' @export

eta_squared_partial <- function(
  object    = NULL,
  F_value   = NULL,
  df_effect = NULL,
  df_error  = NULL
) {
  if (!is.null(object)) {
    return(.eta_squared_from_model(object, value_name = "eta_squared_partial"))
  }
  if (is.null(F_value) || is.null(df_effect) || is.null(df_error)) {
    stop("Either provide 'object' (a fitted aov or lm) or all of 'F_value', 'df_effect', and 'df_error'.")
  }
  .eta_squared_one(F_value, df_effect, df_error,
                      effect_label = "overall",
                      value_name   = "eta_squared_partial")
}
