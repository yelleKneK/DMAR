# Exact expected value of Cohen's d under noncentral t sampling.
#' Exact Expected Value of Cohen's \emph{d} (and Hedges' \emph{g} Bias Correction)
#'
#' Computes \eqn{\mathrm{E}[\hat d \mid \delta, n_1, n_2]}, the exact
#' expected value of the sample standardized mean difference (Cohen's
#' \emph{d}, with pooled variance) under bivariate normality and the
#' noncentral \emph{t} sampling distribution (Hedges, 1981). The sample
#' \emph{d} is upward-biased as an estimator of the population
#' \eqn{\delta}: \eqn{\mathrm{E}[\hat d] = \delta / J(\mathit{df})}, where
#' \deqn{J(\mathit{df}) \;=\; \frac{\Gamma(\mathit{df}/2)}
#'                                  {\sqrt{\mathit{df}/2}\,
#'                                   \Gamma((\mathit{df}-1)/2)}}
#' is Hedges' (1981) bias-correction factor (\eqn{< 1} for finite
#' \eqn{\mathit{df}}, tending to 1 as \eqn{n \to \infty}). Hedges' \emph{g},
#' the unbiased estimator of \eqn{\delta}, is then \eqn{g = J \cdot \hat d}.
#' The same \eqn{J(\mathit{df})} is the workhorse of \code{\link{smd}}
#' when \code{unbiased = TRUE}.
#'
#' \code{expected_smd()} is especially useful at the \emph{design stage} of
#' a study, where sample size planning typically proceeds from an assumed
#' population standardized mean difference \eqn{\delta}. Because the value
#' of \eqn{\hat d} a researcher should expect to observe on average is
#' \emph{larger in absolute value} than \eqn{\delta}, plugging \eqn{\delta}
#' directly into sampling-distribution machinery (precision of \eqn{\hat d},
#' width of a confidence interval, power of a test of
#' \eqn{H_0\!: \delta = 0}) over-promises on the realized precision when
#' the precision is expressed on the \eqn{\hat d} scale. Substituting
#' \code{expected_smd(delta, n_1, n_2)} for the bare \eqn{\delta} corrects
#' this leading-order bias.
#'
#' @param delta Population standardized mean difference. A numeric scalar
#'   or vector.
#' @param n_1 Sample size in the first group. Scalar or vector.
#' @param n_2 Sample size in the second group. Scalar or vector. If
#'   omitted, defaults to \code{n_1} (balanced design).
#'
#' @return A \code{data.frame} with one row per (\code{delta},
#'   \code{n_1}, \code{n_2}) input and the columns
#'   \itemize{
#'     \item \code{delta}, the population SMD,
#'     \item \code{n_1}, \code{n_2}, group sample sizes,
#'     \item \code{expected_smd}, \eqn{\mathrm{E}[\hat d \mid \delta]},
#'     \item \code{bias}, \eqn{\mathrm{E}[\hat d] - \delta} (the amount
#'       by which \eqn{\hat d} overestimates \eqn{\delta} on average),
#'     \item \code{j_correction}, the Hedges (1981) correction factor
#'       \eqn{J(\mathit{df})} used to compute the unbiased \emph{g}, where
#'       \eqn{\mathit{df} = n_1 + n_2 - 2}.
#'   }
#'
#' @details
#' \strong{Derivation.} Under bivariate normality with equal variances,
#' the observed \emph{t}-statistic
#' \eqn{t = \hat d \sqrt{n_1 n_2 / (n_1 + n_2)}} follows a noncentral
#' \emph{t} distribution with \eqn{\mathit{df} = n_1 + n_2 - 2} degrees of
#' freedom and noncentrality parameter
#' \eqn{\lambda = \delta \sqrt{n_1 n_2 / (n_1 + n_2)}}. The expected value
#' of a noncentral \emph{t} variate equals \eqn{\lambda / J(\mathit{df})}
#' (Johnson, Kotz, & Balakrishnan, 1995, Section 31.3), so
#' \eqn{\mathrm{E}[\hat d] = \mathrm{E}[t] / \sqrt{n_1 n_2 / (n_1 + n_2)} =
#'   \delta / J(\mathit{df})}. The bias \eqn{\mathrm{E}[\hat d] - \delta =
#'   \delta\,(1 - J)/J} is positive when \eqn{\delta > 0}.
#'
#' \strong{Magnitude of the correction.} \eqn{J(\mathit{df}) \approx 1 -
#' 3/(4\,\mathit{df} - 1)} to leading order (Hedges & Olkin, 1985, p.\ 81).
#' For \eqn{\delta = 0.5}, \eqn{n_1 = n_2 = 10} (\eqn{\mathit{df} = 18}),
#' \eqn{J \approx 0.957}, so \eqn{\mathrm{E}[\hat d] \approx 0.522} and
#' the upward bias is about 4\%. For \eqn{n_1 = n_2 = 50} the bias is
#' under 1\%; for \eqn{n_1 = n_2 = 5} (very small samples) it exceeds 10\%.
#'
#' \strong{Connection to Hedges' \emph{g}.} The natural inverse of this
#' function is Hedges' \emph{g}: given an observed \eqn{\hat d}, the
#' unbiased estimator of \eqn{\delta} is \eqn{g = J(\mathit{df}) \hat d},
#' which satisfies \eqn{\mathrm{E}[g \mid \delta] = \delta} exactly under
#' the same noncentral \emph{t} model. \code{\link{smd}} with
#' \code{unbiased = TRUE} returns \emph{g}.
#'
#' @references
#' Hedges, L. V. (1981). Distribution theory for Glass's estimator of
#'   effect size and related estimators. \emph{Journal of Educational
#'   Statistics, 6}(2), 107--128.
#'
#' Hedges, L. V., & Olkin, I. (1985). \emph{Statistical methods for
#'   meta-analysis}. Academic Press. (See Section 5, equations 6 and 9.)
#'
#' Johnson, N. L., Kotz, S., & Balakrishnan, N. (1995).
#'   \emph{Continuous univariate distributions, volume 2} (2nd ed.),
#'   Section 31.3. Wiley.
#'
#' Kelley, K. (2007). Confidence intervals for standardized
#'   effect sizes: Theory, application, and implementation. \emph{Journal
#'   of Statistical Software, 20}(8), 1--24. \doi{10.18637/jss.v020.i08}
#'
#' @seealso \code{\link{smd}}, \code{\link{ci_smd}}, \code{\link{ss_aipe_smd}},
#'   \code{\link{expected_r}}, \code{\link{expected_R2}}
#'
#' @examples
#' # 1. Balanced design: delta = 0.5, n = 10 per group.
#' expected_smd(delta = 0.5, n_1 = 10)
#'
#' # 2. Bias as a function of n at fixed delta.
#' expected_smd(delta = 0.5, n_1 = c(5, 10, 20, 50, 100, 500))
#'
#' # 3. Unbalanced design.
#' expected_smd(delta = 0.5, n_1 = 30, n_2 = 60)
#'
#' # 4. Design-stage use: an a-priori delta of 0.4, planned n of 40/group.
#' #        The d we should expect to observe on average is slightly larger.
#' expected_smd(delta = 0.4, n_1 = 40)
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @keywords htest
#'
#' @family effect size estimates
#'
#' @export

expected_smd <- function(delta, n_1, n_2 = NULL) {
  if (!is.numeric(delta)) stop("'delta' must be numeric.")
  if (!is.numeric(n_1))   stop("'n_1' must be numeric.")
  if (is.null(n_2)) n_2 <- n_1
  if (!is.numeric(n_2))   stop("'n_2' must be numeric.")
  if (any(n_1 < 2, na.rm = TRUE) || any(n_2 < 2, na.rm = TRUE))
    stop("'n_1' and 'n_2' must each be >= 2.")

  args <- data.frame(delta = delta, n_1 = n_1, n_2 = n_2)
  df_v <- args$n_1 + args$n_2 - 2L

  # Hedges (1981) J on log scale for numerical stability.
  log_J <- lgamma(df_v / 2) - 0.5 * log(df_v / 2) - lgamma((df_v - 1) / 2)
  J     <- exp(log_J)

  exp_d <- args$delta / J

  .as_dmar_tbl(data.frame(
    delta        = args$delta,
    n_1          = args$n_1,
    n_2          = args$n_2,
    expected_smd = exp_d,
    bias         = exp_d - args$delta,
    j_correction = J,
    stringsAsFactors = FALSE,
    row.names    = NULL
  ))
}
