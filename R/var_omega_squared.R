# Asymptotic variance of omega squared.
#' Asymptotic Variance of Omega Squared (ANOVA Effect Size)
#'
#' Computes the large-sample (delta method) variance of the sample
#' \eqn{\hat\omega^2} (Hays' 1994 bias-corrected estimator) in a
#' fixed-effects ANOVA. Fleishman (1980, Eq. 22, p. 669) gives the exact
#' variance of the unbiased estimator of the signal-to-noise ratio
#' \eqn{f^2 = \sigma^2_a / \sigma^2_e} under the noncentral \emph{F}
#' sampling distribution of the observed \emph{F} statistic; because
#' \eqn{\omega^2 = f^2 / (1 + f^2)} (his Eq. 8), the delta method carries
#' that variance to the \eqn{\omega^2} scale with the Jacobian
#' \eqn{\mathrm{d}\omega^2/\mathrm{d}f^2 = (1 - \omega^2)^2}. Fleishman
#' gives no variance on the \eqn{\omega^2} scale himself, so the transfer
#' is this package's step rather than his. The result is the natural
#' companion to
#' \code{\link{ci_omega_squared}} (CI) and \code{\link{omega_squared}}
#' (point estimate).
#'
#' @param population_omega_squared Population \eqn{\omega^2}. Numeric
#'   scalar in \eqn{[0, 1)}. Ignored when \code{object} is supplied.
#' @param df_effect Numerator degrees of freedom for the effect.
#'   Ignored when \code{object} is supplied.
#' @param df_error Error (residual) degrees of freedom.
#'   Ignored when \code{object} is supplied.
#' @param N Total sample size. Ignored when \code{object} is supplied.
#' @param object Optional fitted \code{\link[stats]{aov}} or
#'   \code{\link[stats]{lm}} object. When supplied, the function loops
#'   over the non-\code{Residuals} rows of \code{\link[stats]{anova}(object)}
#'   and returns one row per effect, plugging in the sample
#'   \eqn{\hat\omega^2_p} for each as the working population value.
#'
#' @return A 1-row \code{data.frame} with columns \code{term}
#'   (\code{"var_omega_squared"}) and \code{value} (the asymptotic
#'   variance).
#'
#' @details
#' \strong{Derivation.} In a fixed-effects ANOVA with numerator df
#' \eqn{df_1} and denominator df \eqn{df_2}, the observed \emph{F}
#' statistic follows a noncentral \emph{F} distribution with
#' noncentrality \eqn{\lambda = df_1 (F - 1)} when
#' \eqn{\hat\omega^2 = df_1 (F - 1) / [df_1 (F - 1) + N]} is the
#' population value (Hays, 1994). The asymptotic variance of
#' \eqn{\hat\omega^2} is obtained by the delta method on this
#' relationship (Fleishman, 1980), yielding:
#' \deqn{\mathrm{Var}(\hat\omega^2) \;\approx\;
#'   \frac{2 \cdot df_1 \cdot (df_2 - 2) (1 - \omega^2)^2 (1 + \lambda^*/df_1)^2}
#'        {N^2 (df_2 - 4)},}
#' with \eqn{\lambda^* = \omega^2 N / (1 - \omega^2)} the noncentrality
#' implied by the population value. This is the form used by
#' \code{\link{ci_omega_squared}} when constructing a Wald-style
#' interval; the noncentral \emph{F} CI is generally preferred.
#'
#' \strong{Caveats.} The variance is a delta method approximation; it
#' becomes inaccurate when \eqn{df_2} is small (\eqn{< 10}), when
#' \eqn{\omega^2} is near the boundaries 0 or 1, or when the residual
#' distribution is heavy-tailed. For small-sample inference, the
#' noncentral \emph{F} CI (\code{\link{ci_omega_squared}}) is preferred
#' over a Wald-style interval built on this variance.
#'
#' @references
#' Fleishman, A. I. (1980). Confidence intervals for correlation ratios.
#'   \emph{Educational and Psychological Measurement, 40}(3), 659--670.
#'
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
#' @seealso \code{\link{omega_squared}}, \code{\link{omega_squared_partial}},
#'   \code{\link{ci_omega_squared}}, \code{\link{ss_aipe_omega_squared}}
#'
#' @examples
#' # 1. One way ANOVA: 3 groups (df_effect = 2), 60 total (df_error = 57),
#' #        population omega^2 = 0.10.
#' var_omega_squared(population_omega_squared = 0.10,
#'                    df_effect = 2,
#'                    df_error = 57,
#'                    N = 60)
#'
#' # 2. Per-effect variance from a fitted lm() / aov() (pygmalion data:
#' #        expectancy treatment x grade, N = 310):
#' fit_factorial <- aov(iq_8 ~ treatment * factor(grade), data = pygmalion)
#' var_omega_squared(object = fit_factorial)
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @keywords htest
#'
#' @family variance utilities
#'
#' @export

var_omega_squared <- function(population_omega_squared = NULL,
                              df_effect = NULL, df_error = NULL, N = NULL,
                              object = NULL) {
  if (!is.null(object)) return(.as_dmar_tbl(.var_omega_squared_from_model(object)))

  if (!is.numeric(population_omega_squared) ||
      length(population_omega_squared) != 1L ||
      population_omega_squared < 0 || population_omega_squared >= 1)
    stop("'population_omega_squared' must be a single value in [0, 1).")
  if (!is.numeric(df_effect) || df_effect < 1)
    stop("'df_effect' must be >= 1.")
  if (!is.numeric(df_error) || df_error < 5)
    stop("'df_error' must be >= 5 for the delta method variance to be defined.")
  if (!is.numeric(N) || N < df_effect + df_error + 1)
    stop("'N' must be >= df_effect + df_error + 1.")

  v <- .var_omega_squared_one(population_omega_squared, df_effect, df_error, N)
  out <- data.frame(term  = "var_omega_squared",
                    value = v,
                    stringsAsFactors = FALSE,
                    row.names = NULL)
  .as_dmar_tbl(out)
}

# Scalar worker. Fleishman (1980, Eq. 22, p. 669) gives the exact variance
# of his unbiased estimator of the signal-to-noise ratio f^2 = sigma^2_a /
# sigma^2_e. His Eq. 8 maps that scale to omega^2 by omega^2 = f^2/(1 + f^2),
# whose derivative is d(omega^2)/d(f^2) = (1 - omega^2)^2, so the delta method
# carries his variance across with the square of that Jacobian, hence the
# fourth power below. Fleishman himself gives no variance on the omega^2
# scale (p. 669: for the correlation ratio one can obtain the interval and
# the median "but not its mean or variance"), so the transfer is DMAR's step
# and is documented as such on the help page.
.var_omega_squared_one <- function(omega2, df_effect, df_error, N) {
  lambda <- omega2 * N / (1 - omega2)
  var_f2 <- 2 * ((df_effect + lambda)^2 +
                 (df_effect + 2 * lambda) * (df_error - 2)) /
            (N^2 * (df_error - 4))
  (1 - omega2)^4 * var_f2
}

# Per-effect variance from a fitted lm/aov.
.var_omega_squared_from_model <- function(object) {
  if (!inherits(object, c("aov", "lm")))
    stop("'object' must be an aov or lm fit.")
  tbl <- stats::anova(object)
  rn <- trimws(rownames(tbl))
  resid_row <- which(rn == "Residuals")
  if (length(resid_row) != 1L)
    stop("Could not find a 'Residuals' row in anova(object).")
  df_error <- tbl[resid_row, "Df"]
  N <- stats::nobs(object)

  rows <- list()
  for (i in seq_len(nrow(tbl))) {
    if (rn[i] == "Residuals") next
    F_val <- tbl[i, "F value"]; df_e <- tbl[i, "Df"]
    if (is.na(F_val) || F_val <= 0) next
    om <- max(0, df_e * (F_val - 1) / (df_e * (F_val - 1) + N))
    if (om >= 1) om <- 1 - 1e-8
    if (df_error <= 4) next
    v <- .var_omega_squared_one(om, df_e, df_error, N)
    rows[[length(rows) + 1L]] <- data.frame(
      effect = rn[i],
      omega_squared_partial = om,
      var_omega_squared = v,
      df_effect = df_e, df_error = df_error, N = N,
      stringsAsFactors = FALSE, row.names = NULL
    )
  }
  if (length(rows) == 0L)
    stop("No testable effects found in the model.")
  do.call(rbind, rows)
}
