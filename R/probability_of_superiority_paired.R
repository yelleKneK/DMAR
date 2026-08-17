# Probability of superiority (paired-samples design).
#' Probability of Superiority for a Paired-Samples Design
#'
#' Computes the probability-of-superiority effect size for paired
#' observations (Grissom & Kim, 2005, 2012), \eqn{P_S = \Pr(Y_1 > Y_2)},
#' along with an analytic confidence interval based on the
#' Brunner-Munzel (2000) U-statistic standard error and a Fisher-
#' \eqn{\mathrm{arctanh}} transformation to keep the bounds inside \eqn{[0, 1]}.
#' The paired counterpart of the Vargha-Delaney (2000) \eqn{A} statistic
#' / \code{\link{cliff_delta}} for two independent groups.
#'
#' @param x,y Paired numeric vectors of equal length. \eqn{x} and
#'   \eqn{y} are interpreted as repeated measurements on the same
#'   units (e.g., pre/post, condition 1 / condition 2, sibling-1 /
#'   sibling-2).
#' @param conf_level Confidence level. Default \code{0.95}.
#'
#' @return A \code{data.frame} with rows for the point estimate of
#'   \eqn{P_S}, the lower / upper CI bounds, the variance, and the
#'   counts of within-pair wins / ties / losses for \eqn{y_1}.
#'
#' @details
#' \strong{Definition.} For paired observations \eqn{(x_i, y_i)},
#' \deqn{P_S \;=\; \Pr(Y > X) + 0.5 \cdot \Pr(Y = X),}
#' where ties are split. The sample estimator is the proportion of
#' pairs with \eqn{y_i > x_i}, plus half the proportion of ties. This is
#' the natural paired-data analog of Vargha-Delaney's \eqn{A} statistic
#' and is unbiased under exchangeability of paired observations.
#'
#' \strong{Why paired-specific.} The independent-groups \code{cles} and
#' \code{cliff_delta} estimators are biased when the two samples are
#' paired, because their variance formulas assume independence of the
#' two groups. For paired data the within-pair correlation reduces the
#' effective sampling variance, which is captured by the Brunner-Munzel
#' (2000) variance used here.
#'
#' \strong{Confidence interval.} The standard error is built from the
#' within-pair sign indicators (Brunner-Munzel, 2000):
#' \deqn{\mathrm{Var}(\hat P_S) \;=\;
#'   \frac{1}{n^2}\sum_{i=1}^{n} (s_i - \bar s)^2,}
#' where \eqn{s_i = \mathrm{I}(y_i > x_i) + 0.5 \cdot \mathrm{I}(y_i = x_i)}.
#' The CI is built on the \eqn{\mathrm{arctanh}(2 P_S - 1)} scale (mapping
#' \eqn{P_S \in [0, 1]} to the real line) and back-transformed to keep
#' the limits inside the unit interval, exactly mirroring
#' \code{\link{cliff_delta}}.
#'
#' @references
#' Brunner, E., & Munzel, U. (2000). The nonparametric Behrens-Fisher
#'   problem: Asymptotic theory and a small-sample approximation.
#'   \emph{Biometrical Journal, 42}(1), 17--25.
#'   \doi{10.1002/(SICI)1521-4036(200001)42:1<17::AID-BIMJ17>3.0.CO;2-U}
#'
#' Grissom, R. J., & Kim, J. J. (2005). \emph{Effect sizes for research:
#'   A broad practical approach}. Lawrence Erlbaum.
#'
#' Grissom, R. J., & Kim, J. J. (2012). \emph{Effect sizes for research:
#'   Univariate and multivariate applications} (2nd ed.). Routledge.
#'
#' Vargha, A., & Delaney, H. D. (2000). A critique and improvement of
#'   the CL common language effect size statistics of McGraw and Wong.
#'   \emph{Journal of Educational and Behavioral Statistics, 25}(2),
#'   101--132. \doi{10.3102/10769986025002101}
#'
#' @seealso \code{\link{cliff_delta}}, \code{\link{cles}},
#'   \code{\link{proportion_of_superiority}}
#'
#' @examples
#' # 1. Paired pre/post data:
#' set.seed(113)
#' pre  <- rnorm(30, mean = 100, sd = 15)
#' post <- pre + rnorm(30, mean =  5, sd = 10)
#' probability_of_superiority_paired(x = pre, y = post)
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @keywords htest
#'
#' @family effect size estimates
#'
#' @export

probability_of_superiority_paired <- function(x, y, conf_level = 0.95) {
  if (!is.numeric(x) || !is.numeric(y))
    stop("'x' and 'y' must both be numeric vectors.")
  if (length(x) != length(y))
    stop("'x' and 'y' must be the same length (paired observations).")
  if (length(x) < 3L)
    stop("Need at least 3 pairs to estimate the variance.")
  if (conf_level <= 0 || conf_level >= 1)
    stop("'conf_level' must be in (0, 1).")

  n <- length(x)
  d <- y - x
  s <- ifelse(d > 0, 1, ifelse(d < 0, 0, 0.5))   # tie-split indicator
  ps_hat <- mean(s)

  wins   <- sum(d > 0)
  losses <- sum(d < 0)
  ties   <- sum(d == 0)

  # Brunner-Munzel (2000) within-pair variance.
  var_ps <- sum((s - ps_hat)^2) / n^2

  # Fisher-arctanh CI on the (2 PS - 1) ∈ [-1, 1] scale, back-transformed.
  z_alpha <- stats::qnorm(1 - (1 - conf_level) / 2)
  d_hat   <- 2 * ps_hat - 1
  z_d     <- atanh(d_hat)
  se_z    <- 2 * sqrt(var_ps) / (1 - d_hat^2)
  d_lo    <- tanh(z_d - z_alpha * se_z)
  d_hi    <- tanh(z_d + z_alpha * se_z)
  ps_lo   <- (d_lo + 1) / 2
  ps_hi   <- (d_hi + 1) / 2

  out <- data.frame(
    term  = c("probability_of_superiority", "lower_limit", "upper_limit",
              "var_ps", "wins_y_over_x", "losses_y_under_x", "ties"),
    value = c(ps_hat, ps_lo, ps_hi, var_ps, wins, losses, ties),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  .as_dmar_tbl(out, conf_level = conf_level)
}
