# Cliff's delta (ordinal dominance effect size).
#' Cliff's \eqn{\delta} Ordinal Effect Size
#'
#' Computes Cliff's (1993) \eqn{\delta} statistic for two independent
#' groups, the difference between the probability that a randomly
#' drawn observation from group 1 exceeds one from group 2 and the
#' reverse probability, together with an analytic confidence
#' interval built from the U-statistic variance (Cliff, 1996). Most
#' R implementations of Cliff's \eqn{\delta} fall back to a bootstrap
#' CI; the analytic CI here is faster, deterministic, and exact in the
#' large-sample limit.
#'
#' @param group_1,group_2 Numeric vectors of observations in the two
#'   groups. Ordinal data are fine; the statistic uses only ranks.
#' @param conf_level Confidence level for the CI. Default \code{0.95}.
#'
#' @return A \code{data.frame} with rows for the point estimate
#'   \code{cliff_delta} and the lower/upper CI bounds. The output also
#'   reports the proportion of pairs with \eqn{y_1 > y_2}, the proportion
#'   with \eqn{y_1 < y_2}, and the proportion of ties.
#'
#' @details
#' \strong{Definition.} Cliff's \eqn{\delta} is
#' \deqn{\delta \;=\; \Pr(Y_1 > Y_2) - \Pr(Y_1 < Y_2)
#'   \;=\; 2 \cdot A - 1,}
#' where \eqn{A} is the Vargha-Delaney (2000) statistic. The sample
#' estimator is
#' \eqn{\hat\delta = (\#\{(i,j): y_{1i} > y_{2j}\}
#'   - \#\{(i,j): y_{1i} < y_{2j}\}) / (n_1 n_2)}.
#' Ties contribute zero to both counts. \eqn{\delta} ranges over
#' \eqn{[-1, 1]}, with 0 indicating no stochastic dominance.
#'
#' \strong{Analytic CI.} The asymptotic variance of \eqn{\hat\delta} is
#' (Cliff, 1993; restated as Feng & Cliff, 2004, Equation 2, p. 323)
#' \deqn{\mathrm{Var}(\hat\delta) \;=\;
#'   \frac{(n_2 - 1) \sigma^2_{d_1} + (n_1 - 1) \sigma^2_{d_2}
#'         + \sigma^2_d}{n_1 n_2},}
#' where \eqn{\sigma^2_{d_i}} is the variance of the per-observation
#' dominance scores within each group. (Feng & Cliff's printed equation
#' transposes the \eqn{(n_1 - 1)} and \eqn{(n_2 - 1)} coefficients, which
#' matters only for unequal group sizes; the pairing above is the correct
#' one, checked by simulation against the empirical variance of
#' \eqn{\hat\delta}.) The CI is constructed on the Fisher-style
#' \eqn{\mathrm{arctanh}}-transformed scale and back-transformed to respect
#' the bounded range of \eqn{\delta} (analogous to Fisher's \emph{Z} CI
#' for Pearson \eqn{r}). Feng & Cliff (2004, Equation 5, p. 324) recommend
#' an alternative asymmetric interval that models the dependence of the
#' variance on \eqn{\delta}; the two constructions agree to first order.
#'
#' \strong{Connection to other measures.} Cliff's \eqn{\delta} is a
#' linear transformation of the Vargha-Delaney (2000) \eqn{A} statistic
#' (\eqn{\delta = 2A - 1}) and of the Mann-Whitney \eqn{U} statistic
#' (\eqn{U / (n_1 n_2) = A}). It is the ordinal analog of the
#' common-language effect size \code{\link{cles}} and is preferable when
#' bivariate normality is implausible (skewed outcomes, ordinal scales).
#'
#' @references
#' Cliff, N. (1993). Dominance statistics: Ordinal analyses to answer
#'   ordinal questions. \emph{Psychological Bulletin, 114}(3), 494--509.
#'   \doi{10.1037/0033-2909.114.3.494}
#'
#' Cliff, N. (1996). \emph{Ordinal methods for behavioral data analysis}.
#'   Lawrence Erlbaum.
#'
#' Feng, D., & Cliff, N. (2004). Monte Carlo evaluation of ordinal d with
#'   improved confidence interval. \emph{Journal of Modern Applied
#'   Statistical Methods, 3}(2), 322--332. \doi{10.22237/jmasm/1099267560}
#'
#' Long, J. D., Feng, D., & Cliff, N. (2003). Ordinal analysis of
#'   behavioral data. In I. B. Weiner (Ed.), \emph{Handbook of
#'   psychology, Vol. 2: Research methods} (pp.\ 635--661). Wiley.
#'
#' Vargha, A., & Delaney, H. D. (2000). A critique and improvement of
#'   the CL common language effect size statistics of McGraw and Wong.
#'   \emph{Journal of Educational and Behavioral Statistics, 25}(2),
#'   101--132. \doi{10.3102/10769986025002101}
#'
#' @seealso \code{\link{cles}}, \code{\link{nnt_from_smd}},
#'   \code{\link{ss_aipe_cliff_delta}}
#'
#' @examples
#' # 1. Two groups of different sizes, no ties:
#' set.seed(113)
#' a <- rnorm(30, mean = 0, sd = 1)
#' b <- rnorm(40, mean = 0.5, sd = 1)
#' cliff_delta(a, b)
#'
#' # 2. With ties (ordinal data):
#' o1 <- c(1, 2, 2, 3, 3, 3, 4, 4, 5)
#' o2 <- c(2, 3, 3, 4, 4, 5, 5, 5)
#' cliff_delta(o1, o2)
#'
#' # 3. Robust to right skew. Cliff's delta on the raw, untransformed
#' #    drinking outcome from the Smith, Meyers, and Delaney (1998)
#' #    trial, comparing the Community Reinforcement Approach (CRA)
#' #    against standard care. Because the statistic uses only ranks it
#' #    needs no normalizing transformation of the heavily skewed
#' #    outcome, unlike the standardized mean difference.
#' data(drinks_trial)
#' cra <- drinks_trial$drinks_per_week[drinks_trial$treatment == "CRA"]
#' std <- drinks_trial$drinks_per_week[drinks_trial$treatment == "Standard"]
#' cliff_delta(cra, std)
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @keywords htest
#'
#' @family effect size estimates
#'
#' @export

cliff_delta <- function(group_1, group_2, conf_level = 0.95) {
  if (!is.numeric(group_1) || !is.numeric(group_2))
    stop("'group_1' and 'group_2' must both be numeric vectors.")
  if (length(group_1) < 2L || length(group_2) < 2L)
    stop("Each group must have at least 2 observations.")
  if (conf_level <= 0 || conf_level >= 1)
    stop("'conf_level' must be in (0, 1).")

  n_1 <- length(group_1)
  n_2 <- length(group_2)

  # Per-pair sign matrix: +1 if y_1 > y_2, -1 if y_1 < y_2, 0 if tie.
  signs <- sign(outer(group_1, group_2, "-"))

  delta_hat <- mean(signs)
  p_gt <- mean(signs == 1)
  p_lt <- mean(signs == -1)
  p_eq <- mean(signs == 0)

  # Cliff (1996) variance:
  d_i <- rowMeans(signs)              # per-group-1 dominance scores
  d_j <- colMeans(signs)              # per-group-2 dominance scores
  sigma_di <- var(d_i)
  sigma_dj <- var(d_j)
  sigma_d  <- var(as.vector(signs))

  var_delta <- ((n_2 - 1) * sigma_di + (n_1 - 1) * sigma_dj + sigma_d) /
               (n_1 * n_2)

  # CI on Fisher-arctanh scale, back-transformed to handle [-1, 1]:
  z_alpha <- stats::qnorm(1 - (1 - conf_level) / 2)
  if (abs(delta_hat) == 1) {
    # Complete separation: the plug-in variance is 0 and arctanh(+/-1) diverges,
    # so se_z = sqrt(var_delta) / (1 - delta_hat^2) is 0 / 0 = NaN. Report the
    # degenerate boundary interval and flag it, rather than returning NaN limits.
    warning("cliff_delta(): the groups are completely separated (|delta| = 1); the estimated variance is 0, so the confidence interval is degenerate at the point estimate.", call. = FALSE)
    lo <- hi <- delta_hat
  } else {
    z_delta <- atanh(delta_hat)
    se_z    <- sqrt(var_delta) / (1 - delta_hat^2)
    lo <- tanh(z_delta - z_alpha * se_z)
    hi <- tanh(z_delta + z_alpha * se_z)
  }

  out <- data.frame(
    term  = c("cliff_delta", "lower_limit", "upper_limit",
              "var_cliff_delta", "p_y1_greater", "p_y1_less", "p_tied"),
    value = c(delta_hat, lo, hi, var_delta, p_gt, p_lt, p_eq),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  .as_dmar_tbl(out, conf_level = conf_level)
}
