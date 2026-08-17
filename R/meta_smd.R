#' Random Effects Meta-Analysis of Standardized Mean Differences
#'
#' Pools two-group standardized mean differences across independent studies.
#' Each study contributes its standardized mean difference and per-group
#' sample sizes; the function computes the within-study sampling variances,
#' applies the Hedges small-sample bias correction by default (the same
#' \eqn{J} factor as \code{\link{expected_smd}} and \code{\link{smd}}), and
#' fits the random effects model of \code{\link{meta_es}}, returning the
#' pooled effect with its confidence interval, tau and tau-squared with
#' intervals, I-squared, Cochran's Q, and a prediction interval for the
#' effect in a new study.
#'
#' @param smd Numeric vector of standardized mean differences (Cohen's
#'   \emph{d}), one per study, positive in the direction of the common
#'   hypothesis.
#' @param n_1,n_2 Per-group sample sizes for each study.
#' @param unbiased Logical: convert each \emph{d} to Hedges \emph{g} (the
#'   small-sample unbiased estimator) before pooling? Default \code{TRUE}.
#'   Set \code{FALSE} to pool the raw \emph{d} values, for example when
#'   reproducing a historical analysis such as Raudenbush (1984) that
#'   predates routine use of the correction.
#' @param method,hartung_knapp,conf_level Passed to \code{\link{meta_es}}:
#'   the tau-squared estimator (\code{"reml"} default), the Hartung-Knapp
#'   small-sample adjustment (default \code{TRUE}), and the confidence
#'   level.
#'
#' @details
#' The within-study variance is the standard large-sample form
#' \deqn{v_i = \frac{n_{1i} + n_{2i}}{n_{1i} n_{2i}} +
#'       \frac{g_i^2}{2 (n_{1i} + n_{2i})},}
#' computed from the bias-corrected \eqn{g_i} when \code{unbiased = TRUE}
#' (Hedges, 1981; Borenstein, Hedges, Higgins, & Rothstein, 2009). All
#' reported quantities are in the standardized mean difference metric.
#'
#' @return A \code{data.frame} (class \code{dmar_tbl}) with the same
#'   rows as \code{\link{meta_es}}.
#'
#' @references
#' Borenstein, M., Hedges, L. V., Higgins, J. P. T., & Rothstein, H. R.
#'   (2009). \emph{Introduction to meta-analysis}. Wiley.
#'
#' Hedges, L. V. (1981). Distribution theory for Glass's estimator of
#'   effect size and related estimators. \emph{Journal of Educational
#'   Statistics, 6}(2), 107--128.
#'
#' Raudenbush, S. W. (1984). Magnitude of teacher expectancy effects on
#'   pupil IQ as a function of the credibility of expectancy induction: A
#'   synthesis of findings from 18 experiments. \emph{Journal of
#'   Educational Psychology, 76}(1), 85--97.
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso \code{\link{meta_es}} for the engine and the reported rows;
#'   \code{\link{smd}} and \code{\link{ci_smd}} for the single-study
#'   quantities; \code{\link{plot_forest}} for the picture;
#'   \code{\link{teacher_expectancy}} for the example data.
#'
#' @family meta-analysis
#'
#' @keywords models
#'
#' @examples
#' # Pool the teacher expectancy studies (Raudenbush, 1984). Hedges g and
#' # the Hartung-Knapp adjustment are on by default; the prediction
#' # interval shows where a new expectancy study would be expected to land.
#' data(teacher_expectancy)
#' meta_smd(smd = teacher_expectancy$d,
#'          n_1 = teacher_expectancy$n_experimental,
#'          n_2 = teacher_expectancy$n_control)
#'
#' @export
meta_smd <- function(smd, n_1, n_2, unbiased = TRUE,
                     method = c("reml", "pm", "dl", "fe"),
                     hartung_knapp = TRUE, conf_level = 0.95) {
  k <- length(smd)
  if (!is.numeric(smd) || k < 2L || anyNA(smd)) {
    stop("'smd' must be two or more standardized mean differences.",
         call. = FALSE)
  }
  for (nm in c("n_1", "n_2")) {
    val <- get(nm)
    if (!is.numeric(val) || length(val) != k || anyNA(val) || any(val < 2) ||
        any(val != round(val))) {
      stop(sprintf("'%s' must give an integer sample size (>= 2) for each ",
                   nm), "study.", call. = FALSE)
    }
  }
  if (!is.logical(unbiased) || length(unbiased) != 1L || is.na(unbiased)) {
    stop("'unbiased' must be TRUE or FALSE.", call. = FALSE)
  }

  yi <- smd
  if (unbiased) {
    df <- n_1 + n_2 - 2
    J  <- exp(lgamma(df / 2) - log(sqrt(df / 2)) - lgamma((df - 1) / 2))
    yi <- J * smd
  }
  N  <- n_1 + n_2
  vi <- N / (n_1 * n_2) + yi^2 / (2 * N)

  meta_es(yi, vi, method = method, hartung_knapp = hartung_knapp,
          conf_level = conf_level)
}
