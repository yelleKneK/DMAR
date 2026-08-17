# Variance of Cohen's d (exact noncentral t and Hedges-Olkin approximation).
#' Variance of Cohen's \emph{d} and Hedges' \emph{g}
#'
#' Computes the variance of the sample standardized mean difference
#' (Cohen's \emph{d}) under bivariate normality and homogeneous variances,
#' using the \emph{exact} noncentral \emph{t} sampling distribution
#' (Hedges, 1981) as the default and reporting the large-sample
#' Hedges-Olkin (1985) approximation alongside for comparison. Optionally
#' returns the variance of Hedges' \emph{g}, the bias-corrected
#' counterpart of \emph{d}.
#'
#' \code{var_smd()} is a stand-alone variance utility: most R packages
#' return only the Hedges-Olkin large-sample approximation and conflate
#' the variances of \emph{d} and \emph{g} (Goulet-Pelletier & Cousineau,
#' 2018). The drift between the exact and approximate forms becomes
#' non-trivial below about \eqn{n = 30} per group and matters whenever
#' \code{var_smd()} feeds into meta-analytic weighting, AIPE planning, or
#' a Wald-style standard-error report.
#'
#' @param delta Population standardized mean difference. Numeric scalar
#'   or vector.
#' @param n_1 Sample size in group 1. Scalar or vector.
#' @param n_2 Sample size in group 2. Scalar or vector. Defaults to
#'   \code{n_1} (balanced design).
#' @param unbiased Logical. If \code{TRUE}, returns the variance of
#'   Hedges' \emph{g} (the bias-corrected estimator); if \code{FALSE}
#'   (the default), returns the variance of Cohen's \emph{d}.
#'
#' @return A \code{data.frame} with rows for the exact (noncentral-
#'   \emph{t}) variance and the Hedges-Olkin large-sample approximation.
#'   Columns are \code{term} (\code{"var_smd_exact"} or
#'   \code{"var_smd_approx"}) and \code{value}.
#'
#' @details
#' \strong{Exact noncentral \emph{t} form.} For
#' \eqn{\hat d = (\bar Y_1 - \bar Y_2)/s_p} with pooled \eqn{s_p}, the
#' rescaled statistic \eqn{t = \hat d \sqrt{n_1 n_2 / (n_1 + n_2)}}
#' follows a noncentral \emph{t} with \eqn{\mathit{df} = n_1 + n_2 - 2}
#' degrees of freedom and noncentrality parameter
#' \eqn{\lambda = \delta \sqrt{n_1 n_2 / (n_1 + n_2)}}. The variance of a
#' noncentral \emph{t} is (Johnson, Kotz, & Balakrishnan, 1995, Sec.\ 31.3)
#' \deqn{\mathrm{Var}(t) \;=\;
#'   \frac{\mathit{df}\,(1 + \lambda^2)}{\mathit{df} - 2} \,-\,
#'   \lambda^2 \, c(\mathit{df})^{2},}
#' where \eqn{c(\mathit{df}) = \sqrt{\mathit{df}/2}\,
#' \Gamma((\mathit{df}-1)/2)\,/\,\Gamma(\mathit{df}/2)}; dividing by the
#' design factor \eqn{n_1 n_2 / (n_1 + n_2)} returns \eqn{\mathrm{Var}(\hat d)}.
#' For Hedges' \emph{g}, multiply the result by \eqn{J(\mathit{df})^2}
#' where \eqn{J(\mathit{df}) = 1/c(\mathit{df})} is the Hedges-Olkin
#' (1985) bias-correction factor (see \code{\link{expected_smd}}).
#'
#' \strong{Hedges-Olkin large-sample approximation.} The frequently
#' quoted approximation (Hedges & Olkin, 1985, equation 8) is
#' \deqn{\mathrm{Var}(\hat d) \;\approx\;
#'   \frac{n_1 + n_2}{n_1 n_2} \;+\;
#'   \frac{\delta^2}{2(n_1 + n_2 - 2)}.}
#' This approaches the exact form only as the degrees of freedom grow:
#' even at \eqn{\delta = 0} it returns \eqn{1/(n_1 n_2 / (n_1 + n_2))}
#' while the exact noncentral \emph{t} variance is
#' \eqn{[\mathit{df}/(\mathit{df} - 2)]/(n_1 n_2 / (n_1 + n_2))}, so the
#' approximation is biased downward by a factor of
#' \eqn{(\mathit{df} - 2)/\mathit{df}}, and the downward bias grows with
#' \eqn{\delta} and small \eqn{n}. Goulet-Pelletier & Cousineau
#' (2018) document the drift and recommend the exact form for
#' \eqn{n < 30} per group.
#'
#' \strong{Companions.} \code{var_smd()} is the variance partner of
#' \code{\link{expected_smd}} (mean) and \code{\link{ci_smd}} (CI). For
#' design-stage AIPE planning that solves for \eqn{n} given a target CI
#' width on \emph{d}, see \code{\link{ss_aipe_smd}}.
#'
#' @references
#' Goulet-Pelletier, J.-C., & Cousineau, D. (2018). A review of effect
#'   sizes and their confidence intervals, Part I: The Cohen's \emph{d}
#'   family. \emph{The Quantitative Methods for Psychology, 14}(4), 242--265.
#'   \doi{10.20982/tqmp.14.4.p242}
#'
#' Hedges, L. V. (1981). Distribution theory for Glass's estimator of
#'   effect size and related estimators. \emph{Journal of Educational
#'   Statistics, 6}(2), 107--128.
#'
#' Hedges, L. V., & Olkin, I. (1985). \emph{Statistical methods for
#'   meta-analysis}. Academic Press.
#'
#' Johnson, N. L., Kotz, S., & Balakrishnan, N. (1995).
#'   \emph{Continuous univariate distributions, volume 2} (2nd ed.).
#'   Wiley.
#'
#' Kelley, K., & Rausch, J. R. (2006). Sample size planning for the
#'   standardized mean difference: Accuracy in parameter estimation via
#'   narrow confidence intervals. \emph{Psychological Methods, 11}(4),
#'   363--385. \doi{10.1037/1082-989X.11.4.363}
#'
#' Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). \emph{Designing
#'   experiments and analyzing data: A model comparison perspective}
#'   (4th ed.). Routledge. (See Chapter 4 on individual comparisons and
#'   Chapter 3 on one-way ANOVA.)
#'
#' @seealso \code{\link{smd}}, \code{\link{ci_smd}},
#'   \code{\link{expected_smd}}, \code{\link{ss_aipe_smd}}
#'
#' @examples
#' # 1. Balanced design, delta = 0.5, n = 20 per group.
#' var_smd(delta = 0.5, n_1 = 20)
#'
#' # 2. Hedges-Olkin approximation drifts from exact when n is small.
#' var_smd(delta = 0.5, n_1 = 5)
#' var_smd(delta = 0.5, n_1 = 50)
#'
#' # 3. Variance of Hedges' g (bias-corrected).
#' var_smd(delta = 0.5, n_1 = 20, unbiased = TRUE)
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @keywords htest
#'
#' @family variance utilities
#'
#' @export

var_smd <- function(delta, n_1, n_2 = NULL, unbiased = FALSE) {
  if (!is.numeric(delta)) stop("'delta' must be numeric.")
  if (!is.numeric(n_1))   stop("'n_1' must be numeric.")
  if (is.null(n_2)) n_2 <- n_1
  if (!is.numeric(n_2))   stop("'n_2' must be numeric.")

  df_v     <- n_1 + n_2 - 2L
  if (any(df_v < 3, na.rm = TRUE))
    stop("'n_1 + n_2 - 2' must be >= 3 for var_smd().")

  design   <- n_1 * n_2 / (n_1 + n_2)
  lambda   <- delta * sqrt(design)

  # Hedges J(df) = 1 / c(df); via log-gamma for stability.
  log_J    <- lgamma(df_v / 2) - 0.5 * log(df_v / 2) - lgamma((df_v - 1) / 2)
  J        <- exp(log_J)
  c_df     <- 1 / J

  # Exact noncentral t variance:
  var_nct  <- df_v * (1 + lambda^2) / (df_v - 2) - lambda^2 * c_df^2
  var_d    <- var_nct / design

  if (unbiased) var_d <- var_d * J^2

  # Hedges-Olkin large-sample approximation:
  var_approx <- (n_1 + n_2) / (n_1 * n_2) + delta^2 / (2 * (df_v))
  if (unbiased) var_approx <- var_approx * J^2

  out <- data.frame(
    term  = c("var_smd_exact", "var_smd_approx"),
    value = c(var_d, var_approx),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  .as_dmar_tbl(out)
}
