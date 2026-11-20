# Exact expected value of Cohen's d under noncentral t sampling.
#' Exact Expected Value of Cohen's \emph{d} (and Hedges' \emph{g} Bias Correction)
#'
#' Computes \eqn{\mathrm{E}[\hat d \mid \delta, n_1, n_2]}, the exact
#' expected value of the sample standardized mean difference (Cohen's
#' \emph{d}, with pooled variance) under normality within each group with
#' a common variance, from the noncentral \emph{t} sampling distribution
#' (Hedges, 1981). The sample
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
#' \code{expected_smd(true_smd, n_1, n_2)} for the bare \eqn{\delta} corrects
#' this leading-order bias.
#'
#' @param true_smd Population standardized mean difference, \eqn{\delta}
#'   in the formulas. A numeric scalar or vector.
#' @param n_1 Sample size in the first group. Scalar or vector.
#' @param n_2 Sample size in the second group. Scalar or vector. If
#'   omitted, defaults to \code{n_1} (balanced design).
#'
#' @return A \code{data.frame} with one row per (\code{true_smd},
#'   \code{n_1}, \code{n_2}) input and the columns
#'   \itemize{
#'     \item \code{true_smd}, the population SMD,
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
#' \strong{Derivation.} When the scores are normal within each group with
#' a common variance,
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
#' # 1. Balanced design: true_smd = 0.5, n = 10 per group.
#' expected_smd(true_smd = 0.5, n_1 = 10)
#'
#' # 2. Bias as a function of n at a fixed true_smd.
#' expected_smd(true_smd = 0.5, n_1 = c(5, 10, 20, 50, 100, 500))
#'
#' # 3. Unbalanced design.
#' expected_smd(true_smd = 0.5, n_1 = 30, n_2 = 60)
#'
#' # 4. Design-stage use: an a-priori true_smd of 0.4, planned n of 40/group.
#' #        The d we should expect to observe on average is slightly larger.
#' expected_smd(true_smd = 0.4, n_1 = 40)
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @keywords htest
#'
#' @family effect size estimates
#'
#' @export

expected_smd <- function(true_smd, n_1, n_2 = NULL) {
  if (!is.numeric(true_smd)) .dmar_refuse_argument(
    "'true_smd' must be numeric.", arg = "true_smd")
  if (!is.numeric(n_1))   .dmar_refuse_argument(
    "'n_1' must be numeric.", arg = "n_1")
  if (is.null(n_2)) n_2 <- n_1
  if (!is.numeric(n_2))   .dmar_refuse_argument(
    "'n_2' must be numeric.", arg = "n_2")
  if (any(n_1 < 2, na.rm = TRUE) || any(n_2 < 2, na.rm = TRUE))
    .dmar_refuse_argument(
      "'n_1' and 'n_2' must each be >= 2.", arg = c("n_1", "n_2"))

  args <- data.frame(true_smd = true_smd, n_1 = n_1, n_2 = n_2)
  df_v <- args$n_1 + args$n_2 - 2L

  # Hedges (1981) J from its logarithm (R/c4_internals.R). The bias,
  # true_smd (1 / J - 1), is taken from the logarithm with expm1(): 1 / J
  # is within 1e-8 of 1 from about df = 1e8 on, and subtracting true_smd
  # from true_smd / J there keeps few correct digits or none.
  log_J <- .hedges_log_J(df_v)
  J     <- exp(log_J)

  exp_d <- args$true_smd * exp(-log_J)

  # true_smd is a measured value, so a whole one prints to the decimals of
  # the rest, as it did under its DMAR 1.0.x name, delta (a measured name
  # of the display).
  .as_dmar_tbl(data.frame(
    true_smd     = args$true_smd,
    n_1          = args$n_1,
    n_2          = args$n_2,
    expected_smd = exp_d,
    bias         = args$true_smd * expm1(-log_J),
    j_correction = J,
    stringsAsFactors = FALSE,
    row.names    = NULL
  ), measured_terms = "true_smd")
}
