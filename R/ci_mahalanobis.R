#' Confidence Interval for the Squared Mahalanobis Distance
#'
#' Computes the squared Mahalanobis distance \eqn{D^2} together with an exact
#' confidence interval for the population squared distance \eqn{\Delta^2},
#' obtained by inverting Hotelling's \eqn{T^2} statistic through its
#' (noncentral) \emph{F}-distribution as in Reiser (2001). Both the one-sample
#' setting (mean vector against a hypothesized population mean) and the
#' two-sample setting (between-groups distance from discriminant analysis or
#' multivariate group comparison) are supported, and either raw data or a
#' pre-computed \eqn{D^2} with sample sizes can be supplied.
#'
#' @param D2 Optional pre-computed squared Mahalanobis distance. Ignored if
#'   \code{group_1} is supplied.
#' @param group_1 Optional numeric matrix or data frame for the first sample
#'   (\eqn{n_1 \times p}). Rows are observations and columns are variables.
#' @param group_2 Optional numeric matrix or data frame for the second sample
#'   (\eqn{n_2 \times p}). When \code{NULL}, the function operates in
#'   one-sample mode against \code{mu_0}.
#' @param mu_0 Optional hypothesized population mean for the one-sample case
#'   (length-\eqn{p} numeric vector). Defaults to a vector of zeros.
#' @param n_1 Sample size for group 1 (required when supplying \code{D2}
#'   directly).
#' @param n_2 Sample size for group 2 (required for two-sample mode when
#'   supplying \code{D2} directly; leave \code{NULL} for one-sample mode).
#' @param p Dimensionality (number of variables) when supplying \code{D2}
#'   directly.
#' @param conf_level Confidence coverage for a symmetric interval (default
#'   \code{0.95}). Set it to \code{NULL} to specify the tails directly
#'   through \code{alpha_lower} and \code{alpha_upper}.
#' @param alpha_lower,alpha_upper Optional Type I error rates for the lower and
#'   upper tail. To use them, set \code{conf_level = NULL} and supply both
#'   (an asymmetric or one-sided interval with coverage
#'   \code{1 - alpha_lower - alpha_upper}); supplying either alongside a
#'   non-\code{NULL} \code{conf_level} is an error, as in
#'   \code{\link{conf_limits_ncf}}, to which they are passed.
#' @param \dots Additional arguments passed to \code{\link{conf_limits_ncf}}
#'   (for example \code{tol}).
#'
#' @return A one-row \code{data.frame} with columns \code{sample_type}
#'   (\code{"one-sample"} or \code{"two-sample"}), \code{D2} (point estimate of
#'   the squared distance), \code{lower_limit} and \code{upper_limit} (the
#'   confidence limits on the population squared distance \eqn{\Delta^2}),
#'   \code{F_value}, \code{df_1}, \code{df_2}, \code{n_1}, \code{n_2}
#'   (\code{NA} in one-sample mode), and \code{p}.
#'
#' @details
#' \strong{Definition.} For a \eqn{p}-vector \eqn{\mathbf{x}} drawn from a
#' multivariate normal with mean \eqn{\boldsymbol{\mu}} and covariance
#' \eqn{\boldsymbol{\Sigma}}, Mahalanobis's (1936) squared distance from a
#' reference vector \eqn{\boldsymbol{\mu}_0} is
#' \deqn{\Delta^2 = (\boldsymbol{\mu} - \boldsymbol{\mu}_0)^\top
#'                   \boldsymbol{\Sigma}^{-1}
#'                   (\boldsymbol{\mu} - \boldsymbol{\mu}_0).}
#' In the two-sample case the population distance between groups is
#' \eqn{\Delta^2 = (\boldsymbol{\mu}_1 - \boldsymbol{\mu}_2)^\top
#'                  \boldsymbol{\Sigma}^{-1}
#'                  (\boldsymbol{\mu}_1 - \boldsymbol{\mu}_2)}, assuming a
#' common covariance. The corresponding sample estimates plug the sample means
#' and the sample (or pooled) covariance into the same quadratic form.
#'
#' \strong{Link to Hotelling's \eqn{T^2}.} Hotelling's (1931)
#' \eqn{T^2} statistic is
#' \eqn{T^2 = n D^2} (one sample) or
#' \eqn{T^2 = \{n_1 n_2 / (n_1 + n_2)\} D^2} (two samples). Under
#' multivariate normality
#' \deqn{\frac{n_1 + n_2 - p - 1}{(n_1 + n_2 - 2)\,p}\,T^2 \sim
#'         F'\!\left(p,\, n_1 + n_2 - p - 1,\;
#'                   \lambda = \frac{n_1 n_2}{n_1 + n_2}\,\Delta^2\right)}
#' in the two-sample case, and analogously
#' \eqn{\{(n - p)/[(n-1)p]\}\,T^2 \sim F'(p, n-p, n\Delta^2)} in the
#' one-sample case (see Anderson, 2003, Section 5.2).
#'
#' \strong{Confidence interval.} The CI on \eqn{\Delta^2} is obtained by
#' inverting these distributional results (Reiser, 2001): a CI on the
#' noncentrality parameter \eqn{\lambda} is constructed via
#' \code{\link{conf_limits_ncf}} and then mapped back to \eqn{\Delta^2} by
#' \eqn{\Delta^2 = \lambda\,(n_1 + n_2)/(n_1 n_2)} (two sample) or
#' \eqn{\Delta^2 = \lambda / n} (one sample). When the observed \eqn{F} is
#' below the lower-tail critical value of the central \emph{F}-distribution at
#' the requested confidence level, the lower CI on \eqn{\lambda} (and hence
#' on \eqn{\Delta^2}) is clamped to zero, in keeping with the
#' \code{\link{conf_limits_ncf}} convention.
#'
#' \strong{Bias.} The plug-in estimator \eqn{D^2} is upward biased for
#' \eqn{\Delta^2}; the CI from this function is exact for \eqn{\Delta^2}
#' under multivariate normality and reflects the bias structure correctly,
#' but the point estimate reported is the standard plug-in \eqn{D^2}.
#'
#' @references
#' Anderson, T. W. (2003). \emph{An Introduction to Multivariate Statistical
#'   Analysis} (4th ed.). Wiley.
#'
#' Hotelling, H. (1931). The generalization of Student's ratio.
#'   \emph{The Annals of Mathematical Statistics, 2}(3), 360--378.
#'
#' Mahalanobis, P. C. (1936). On the generalized distance in statistics.
#'   \emph{Proceedings of the National Institute of Sciences of India, 2}(1),
#'   49--55.
#'
#' Reiser, B. (2001). Confidence intervals for the Mahalanobis distance.
#'   \emph{Communications in Statistics--Simulation and Computation, 30}(1),
#'   37--45. \doi{10.1081/SAC-100001856}
#'
#' Steiger, J. H., & Fouladi, R. T. (1997). Noncentrality interval
#'   estimation and the evaluation of statistical methods. In L. L. Harlow,
#'   S. A. Mulaik, & J. H. Steiger (Eds.), \emph{What if there were no
#'   significance tests?} (pp. 221--257). Mahwah, NJ: Lawrence Erlbaum.
#'
#' Kelley, K. (2007). Confidence intervals for standardized
#'   effect sizes: Theory, application, and implementation.
#'   \emph{Journal of Statistical Software, 20}(8), 1--24.
#'   \doi{10.18637/jss.v020.i08}
#'
#' @examples
#' # Two-sample distance between the two schools of the Holzinger and
#' #     Swineford (1939) study on a four-test cognitive battery.
#' battery <- c("t1_visual_perception", "t2_cubes", "t4_lozenges",
#'              "t6_paragraph_comprehension")
#' g1 <- as.matrix(holzinger_swineford[
#'   holzinger_swineford$school == "Grant-White", battery])
#' g2 <- as.matrix(holzinger_swineford[
#'   holzinger_swineford$school == "Pasteur", battery])
#' ci_mahalanobis(group_1 = g1, group_2 = g2)
#'
#' # One-sample distance: how far is the Grant-White centroid from a
#' #     reference vector of (29, 24, 18, 9)?
#' ci_mahalanobis(group_1 = g1, mu_0 = c(29, 24, 18, 9))
#'
#' # Pre-computed D^2 (no raw data needed): the two-school distance,
#' #     reproduced from reported summaries.
#' ci_mahalanobis(D2 = 0.608, n_1 = 145, n_2 = 156, p = 4)
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso \code{\link{conf_limits_ncf}}, \code{\link{ci_smd}},
#'   \code{\link{ci_R2}}
#'
#' @keywords multivariate htest
#'
#' @family confidence intervals for effect sizes
#'
#' @export
ci_mahalanobis <- function(D2          = NULL,
                           group_1     = NULL,
                           group_2     = NULL,
                           mu_0        = NULL,
                           n_1         = NULL,
                           n_2         = NULL,
                           p           = NULL,
                           conf_level  = 0.95,
                           alpha_lower = NULL,
                           alpha_upper = NULL,
                           ...) {

  raw_data <- !is.null(group_1)

  if (raw_data) {
    g1 <- as.matrix(group_1)
    if (!is.numeric(g1)) stop("'group_1' must be numeric.", call. = FALSE)
    if (anyNA(g1)) stop("Missing values are not supported in 'group_1'.", call. = FALSE)
    n_1 <- nrow(g1)
    p   <- ncol(g1)
    xbar_1 <- colMeans(g1)
    S_1    <- stats::cov(g1)

    if (!is.null(group_2)) {
      g2 <- as.matrix(group_2)
      if (!is.numeric(g2)) stop("'group_2' must be numeric.", call. = FALSE)
      if (anyNA(g2)) stop("Missing values are not supported in 'group_2'.", call. = FALSE)
      if (ncol(g2) != p) {
        stop("'group_1' and 'group_2' must have the same number of variables.",
             call. = FALSE)
      }
      n_2 <- nrow(g2)
      xbar_2 <- colMeans(g2)
      S_2    <- stats::cov(g2)
      S_p    <- ((n_1 - 1) * S_1 + (n_2 - 1) * S_2) / (n_1 + n_2 - 2)
      delta  <- xbar_1 - xbar_2
      D2     <- as.numeric(crossprod(delta, solve(S_p, delta)))
      two_sample <- TRUE
    } else {
      if (is.null(mu_0)) {
        mu_0 <- rep(0, p)
      } else if (length(mu_0) != p) {
        stop("'mu_0' must have length p = ", p, ".", call. = FALSE)
      }
      delta <- xbar_1 - mu_0
      D2    <- as.numeric(crossprod(delta, solve(S_1, delta)))
      two_sample <- FALSE
    }
  } else {
    if (is.null(D2) || is.null(p) || is.null(n_1)) {
      stop("Either supply raw data via 'group_1' (and optionally 'group_2'), ",
           "or supply 'D2', 'n_1', and 'p' (plus 'n_2' for a two-sample ",
           "distance).", call. = FALSE)
    }
    two_sample <- !is.null(n_2)
  }

  if (!is.numeric(D2) || length(D2) != 1L || !is.finite(D2) || D2 < 0) {
    stop("'D2' must be a single non-negative finite number.", call. = FALSE)
  }
  if (p < 1L) stop("'p' must be at least 1.", call. = FALSE)
  if (n_1 < 2L) stop("'n_1' must be at least 2.", call. = FALSE)
  if (two_sample && n_2 < 2L) stop("'n_2' must be at least 2.", call. = FALSE)

  if (two_sample) {
    df_1 <- p
    df_2 <- n_1 + n_2 - p - 1L
    if (df_2 <= 0L) {
      stop("Insufficient degrees of freedom: need n_1 + n_2 > p + 1.",
           call. = FALSE)
    }
    T2    <- (n_1 * n_2 / (n_1 + n_2)) * D2
    F_obs <- T2 * df_2 / (p * (n_1 + n_2 - 2L))
    lambda_to_D2 <- function(lam) lam * (n_1 + n_2) / (n_1 * n_2)
  } else {
    df_1 <- p
    df_2 <- n_1 - p
    if (df_2 <= 0L) {
      stop("Insufficient degrees of freedom: need n_1 > p.", call. = FALSE)
    }
    T2    <- n_1 * D2
    F_obs <- T2 * df_2 / (p * (n_1 - 1L))
    lambda_to_D2 <- function(lam) lam / n_1
  }

  ncp_ci <- .conf_limits_ncf_for(
    caller      = "ci_mahalanobis",
    quantity    = "the squared Mahalanobis distance",
    F_value     = F_obs,
    df_1        = df_1,
    df_2        = df_2,
    conf_level  = conf_level,
    alpha_lower = alpha_lower,
    alpha_upper = alpha_upper,
    verbose     = FALSE,
    ...
  )
  lam_lower <- ncp_ci$value[ncp_ci$term == "lower_limit"]
  lam_upper <- ncp_ci$value[ncp_ci$term == "upper_limit"]

  D2_lower <- lambda_to_D2(lam_lower)
  D2_upper <- lambda_to_D2(lam_upper)

  out <- data.frame(
    sample_type = if (two_sample) "two-sample" else "one-sample",
    D2          = D2,
    lower_limit = D2_lower,
    upper_limit = D2_upper,
    F_value     = F_obs,
    df_1        = df_1,
    df_2        = df_2,
    n_1         = n_1,
    n_2         = if (two_sample) n_2 else NA_integer_,
    p           = p,
    stringsAsFactors = FALSE,
    row.names   = NULL
  )
  # Wide-format dmar_tbl: a leading label column (sample_type) beside several
  # typed numeric columns. The format method works column-by-column, so the
  # degrees of freedom and sample sizes print as whole numbers while D2 and the
  # interval limits print to significant figures. There is no p-value column
  # here (the `p` column is the dimensionality, an integer count of variables).
  .as_dmar_tbl(out, conf_level = conf_level)
}
