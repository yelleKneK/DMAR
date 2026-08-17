# Maximal reliability coefficient H (Hancock & Mueller, 2001).
#' Maximal Reliability Coefficient \emph{H} (Hancock & Mueller, 2001)
#'
#' Computes the Hancock-Mueller (2001) maximal-reliability coefficient
#' \emph{H} from a vector of standardized factor loadings; \emph{H} is
#' the reliability of the optimally-weighted composite of a set of
#' indicators of a single latent construct. It is uniformly greater than
#' or equal to coefficient alpha and McDonald's omega for the same data,
#' so it sets a useful upper bound on what reliability can plausibly be
#' for that indicator set. A delta method confidence interval is reported
#' when standard errors of the standardized loadings are supplied.
#'
#' @param loadings Numeric vector of standardized factor loadings,
#'   each in \eqn{(-1, 1)}. At least 2 loadings are required.
#' @param se_loadings Optional vector of standard errors of the
#'   standardized loadings (same length as \code{loadings}). When
#'   supplied, a delta method CI on \emph{H} is reported.
#' @param conf_level Confidence level for the CI. Default \code{0.95}.
#'
#' @return A \code{data.frame} with rows for the point estimate
#'   \code{reliability_H} and (when SEs are supplied) the lower / upper
#'   CI bounds and the delta method variance.
#'
#' @details
#' \strong{Definition.} For \eqn{p} indicators of a single latent factor
#' with standardized loadings \eqn{\lambda_1, \ldots, \lambda_p},
#' Hancock & Mueller (2001) showed that the maximum reliability
#' achievable by any linear composite of the indicators is
#' \deqn{H \;=\;
#'   \frac{\sum_{i=1}^{p} \lambda_i^2 / (1 - \lambda_i^2)}
#'        {1 + \sum_{i=1}^{p} \lambda_i^2 / (1 - \lambda_i^2)}.}
#' Equivalently, defining \eqn{\theta_i = \lambda_i^2 / (1 - \lambda_i^2)}
#' (the signal-to-noise ratio for indicator \eqn{i}), \eqn{H =
#' \sum \theta_i / (1 + \sum \theta_i)}. As \eqn{p} grows or as the
#' individual loadings grow toward 1, \eqn{H \to 1}.
#'
#' \strong{Relationship to coefficient alpha and omega.} Coefficient
#' alpha (\code{\link{reliability_alpha}}) is the reliability of the
#' \emph{equally-weighted} sum of indicators; \emph{H} is the reliability
#' of the \emph{optimally-weighted} composite. Hancock & Mueller (2001)
#' prove \eqn{H \ge \omega \ge \alpha} for a unidimensional indicator
#' set, with equality only when all loadings are equal. \emph{H} is
#' therefore most useful for diagnostics: if \emph{H} is much higher
#' than alpha, the standard composite is leaving reliability on the
#' table.
#'
#' \strong{Confidence interval via the delta method.} Conditional on
#' standard errors \eqn{\mathrm{SE}(\hat \lambda_i)}, the delta method
#' variance of \emph{H} is
#' \deqn{\mathrm{Var}(\hat H) \;\approx\;
#'   \sum_{i=1}^{p} \left(\frac{\partial H}{\partial \lambda_i}\right)^2
#'     \mathrm{SE}(\hat \lambda_i)^2,}
#' with
#' \eqn{\partial H / \partial \lambda_i = 2 \lambda_i /
#'      [(1 - \lambda_i^2)^2 (1 + \sum_j \theta_j)^2]}.
#' The CI is built on the \eqn{\mathrm{logit}(H)} scale (mapping
#' \eqn{[0, 1]} to the real line) and back-transformed, as recommended
#' by Browne (1968) for bounded reliability coefficients.
#'
#' @references
#' Browne, M. W. (1968). A comparison of factor analytic techniques.
#'   \emph{Psychometrika, 33}(3), 267--334.
#'
#' Hancock, G. R., & Mueller, R. O. (2001). Rethinking construct
#'   reliability within latent variable systems. In R. Cudeck, S. du
#'   Toit, & D. Sörbom (Eds.), \emph{Structural equation modeling:
#'   Present and future} (pp. 195--216). Scientific Software
#'   International.
#'
#' Kelley, K., & Cheng, Y. (2012). Estimation of and confidence interval
#'   formation for reliability coefficients of homogeneous measurement
#'   instruments. \emph{Methodology, 8}, 39--50.
#'   \doi{10.1027/1614-2241/a000036}
#'
#' Kelley, K., & Pornprasertmanit, S. (2016). Confidence intervals for
#'   population reliability coefficients: Evaluation of methods,
#'   recommendations, and software for composite measures.
#'   \emph{Psychological Methods, 21}, 69--92. \doi{10.1037/a0040086}
#'
#' Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). \emph{Designing
#'   experiments and analyzing data: A model comparison perspective}
#'   (4th ed.). Routledge.
#'
#' Raykov, T. (1997). Estimation of composite reliability for congeneric
#'   measures. \emph{Applied Psychological Measurement, 21}(2), 173--184.
#'   \doi{10.1177/01466216970212006}
#'
#' Terry, L. J., & Kelley, K. (2012). Sample size planning for composite
#'   reliability coefficients: Accuracy in parameter estimation via narrow
#'   confidence intervals. \emph{British Journal of Mathematical and
#'   Statistical Psychology, 65}, 371--401.
#'   \doi{10.1111/j.2044-8317.2011.02030.x}
#'
#' @seealso \code{\link{reliability_alpha}}, \code{\link{reliability_omega}},
#'   \code{\link{reliability}}
#'
#' @examples
#' # 1. Five indicators with standardized loadings 0.6, 0.7, ..., 0.8:
#' reliability_H(loadings = c(0.6, 0.65, 0.70, 0.75, 0.80))
#'
#' # 2. With per-loading standard errors from a CFA output:
#' reliability_H(loadings    = c(0.6, 0.7, 0.8),
#'                se_loadings = c(0.05, 0.04, 0.03))
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @keywords htest
#'
#' @family reliability
#'
#' @export

reliability_H <- function(loadings, se_loadings = NULL, conf_level = 0.95) {
  if (!is.numeric(loadings) || length(loadings) < 2L)
    stop("'loadings' must be a numeric vector with length >= 2.")
  if (any(abs(loadings) >= 1))
    stop("Each entry of 'loadings' must be in (-1, 1).")

  theta <- loadings^2 / (1 - loadings^2)
  sum_t <- sum(theta)
  H_hat <- sum_t / (1 + sum_t)

  out <- data.frame(
    term  = "reliability_H",
    value = H_hat,
    stringsAsFactors = FALSE,
    row.names = NULL
  )

  if (!is.null(se_loadings)) {
    if (length(se_loadings) != length(loadings))
      stop("'se_loadings' must be the same length as 'loadings'.")
    if (any(se_loadings < 0))
      stop("'se_loadings' must be non-negative.")

    # dH/dlambda_i:
    denom <- (1 + sum_t)^2
    dH_dl <- 2 * loadings / ((1 - loadings^2)^2 * denom)
    var_H <- sum(dH_dl^2 * se_loadings^2)

    # logit-transform CI:
    z_alpha <- stats::qnorm(1 - (1 - conf_level) / 2)
    logit_H <- log(H_hat / (1 - H_hat))
    se_logit <- sqrt(var_H) / (H_hat * (1 - H_hat))
    H_lo <- 1 / (1 + exp(-(logit_H - z_alpha * se_logit)))
    H_hi <- 1 / (1 + exp(-(logit_H + z_alpha * se_logit)))

    out <- rbind(
      out,
      data.frame(
        term  = c("lower_limit", "upper_limit", "var_H"),
        value = c(H_lo, H_hi, var_H),
        stringsAsFactors = FALSE,
        row.names = NULL
      )
    )
  }

  .as_dmar_tbl(out, conf_level = conf_level)
}
