#' Random Effects Meta-Analysis of Generic Effect Sizes
#'
#' Pools independent effect sizes given their sampling variances: the general
#' engine behind \code{\link{meta_smd}} and \code{\link{meta_r}}, exposed for
#' any effect metric whose estimates are approximately normal with known
#' variances. The random effects model is the default and the fit reports the
#' full uncertainty picture in one table: the pooled estimate with its
#' confidence interval, the between-study standard deviation tau, the
#' between-study variance tau-squared with its Q-profile confidence
#' interval, I-squared with an interval mapped from the tau-squared limits,
#' H-squared, Cochran's Q test, and, always, a prediction interval for the
#' effect in a new study. Reporting the prediction interval by default is
#' deliberate: when heterogeneity is real, the confidence interval for the
#' average effect understates what the next study will show, and the package
#' treats \dQuote{where will the next study land} as part of the answer, not
#' an option.
#'
#' @param yi Numeric vector of effect sizes, one per independent study.
#' @param vi Sampling variances of \code{yi}, one per study.
#' @param method Between-study variance estimator: \code{"reml"} (restricted
#'   maximum likelihood, the default), \code{"pm"} (Paule-Mandel),
#'   \code{"dl"} (DerSimonian-Laird), or \code{"fe"} (a fixed effect /
#'   common effect analysis, which assumes tau-squared is zero and reports
#'   no prediction interval).
#' @param hartung_knapp Logical: apply the Hartung-Knapp-Sidik-Jonkman
#'   small-sample adjustment (the pooled standard error rescaled from the
#'   weighted residuals, with a \emph{t} reference on \eqn{k - 1} degrees of
#'   freedom)? Default \code{TRUE}: with the small numbers of studies
#'   typical in psychology and education it keeps the confidence interval
#'   near its nominal coverage, where the conventional normal interval is
#'   anticonservative. Ignored for \code{method = "fe"}.
#' @param conf_level Confidence level for all intervals. Defaults to 0.95.
#'
#' @details
#' The model is \eqn{y_i = \mu + u_i + e_i} with \eqn{u_i \sim N(0, \tau^2)}
#' and \eqn{e_i \sim N(0, v_i)}, \eqn{v_i} treated as known. The
#' \eqn{\tau^2} confidence interval inverts the generalized Q statistic
#' (Viechtbauer, 2007); the I-squared interval maps the \eqn{\tau^2}
#' interval through the typical within-study variance of Higgins and
#' Thompson (2002). The prediction interval follows Higgins,
#' Thompson, and Spiegelhalter (2009), using \emph{t} with \eqn{k - 2}
#' degrees of freedom, and requires at least three studies.
#'
#' I-squared is reported because readers expect it, but note its
#' well-known limitation: it is a \emph{proportion} of variability, not an
#' amount, so the same tau matched with larger studies yields a larger
#' I-squared. The quantity with direct scientific meaning is tau (the
#' between-study standard deviation, in the metric of \code{yi}) together
#' with the prediction interval.
#'
#' @return A \code{data.frame} (class \code{dmar_tbl}) with rows
#'   \code{estimate}, \code{se}, the test statistic (\code{t} under
#'   Hartung-Knapp, \code{z} otherwise), \code{p_value},
#'   \code{lower_limit} / \code{upper_limit}, \code{prediction_lower} /
#'   \code{prediction_upper}, \code{tau2} with \code{tau2_lower} /
#'   \code{tau2_upper}, \code{tau}, \code{I2} with limits, \code{H2},
#'   \code{Q} / \code{Q_df} / \code{Q_p}, and \code{k}. The estimator and
#'   adjustment are recorded in the \code{"method"} and
#'   \code{"hartung_knapp"} attributes.
#'
#' @references
#' DerSimonian, R., & Laird, N. (1986). Meta-analysis in clinical trials.
#'   \emph{Controlled Clinical Trials, 7}(3), 177--188.
#'
#' Hartung, J., & Knapp, G. (2001). On tests of the overall treatment effect
#'   in meta-analysis with normally distributed responses.
#'   \emph{Statistics in Medicine, 20}(12), 1771--1782. \doi{10.1002/sim.791}
#'
#' Higgins, J. P. T., & Thompson, S. G. (2002). Quantifying heterogeneity in
#'   a meta-analysis. \emph{Statistics in Medicine, 21}(11), 1539--1558.
#'   \doi{10.1002/sim.1186}
#'
#' Higgins, J. P. T., Thompson, S. G., & Spiegelhalter, D. J. (2009). A
#'   re-evaluation of random-effects meta-analysis. \emph{Journal of the
#'   Royal Statistical Society: Series A, 172}(1), 137--159.
#'   \doi{10.1111/j.1467-985X.2008.00552.x}
#'
#' Viechtbauer, W. (2007). Confidence intervals for the amount of
#'   heterogeneity in meta-analysis. \emph{Statistics in Medicine, 26}(1),
#'   37--52. \doi{10.1002/sim.2514}
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso \code{\link{meta_smd}} and \code{\link{meta_r}} for the metric-
#'   specific front ends; \code{\link{meta_contrast}} for focused moderator
#'   contrasts; \code{\link{combine_p}} for combined significance tests;
#'   \code{\link{plot_forest}} to see the studies and the pool together.
#'
#' @family meta-analysis
#'
#' @keywords models
#'
#' @examples
#' # The teacher expectancy studies (Raudenbush, 1984), pooled in the d
#' # metric with variances from the standard large-sample formula.
#' data(teacher_expectancy)
#' d <- teacher_expectancy$d
#' n_e <- teacher_expectancy$n_experimental
#' n_c <- teacher_expectancy$n_control
#' v <- (n_e + n_c) / (n_e * n_c) + d^2 / (2 * (n_e + n_c))
#' meta_es(d, v)
#'
#' # A fixed effect (common effect) analysis of the same studies.
#' meta_es(d, v, method = "fe")
#'
#' @export
meta_es <- function(yi, vi, method = c("reml", "pm", "dl", "fe"),
                    hartung_knapp = TRUE, conf_level = 0.95) {
  .meta_check_yivi(yi, vi)
  method <- match.arg(method)
  if (!is.logical(hartung_knapp) || length(hartung_knapp) != 1L ||
      is.na(hartung_knapp)) {
    stop("'hartung_knapp' must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.numeric(conf_level) || length(conf_level) != 1L ||
      is.na(conf_level) || conf_level <= 0 || conf_level >= 1) {
    stop("'conf_level' must be a single number in (0, 1).", call. = FALSE)
  }
  fit <- .meta_fit(yi, vi, method = method, hartung_knapp = hartung_knapp,
                   conf_level = conf_level)
  .meta_table(fit, conf_level = conf_level)
}
