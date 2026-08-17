#' Power of the TOST or Noninferiority Test for a Linear Contrast
#'
#' Computes the exact power of the Schuirmann (1987) two one-sided
#' tests procedure, or of the one-sided noninferiority test, for a
#' linear contrast of group means \eqn{\psi = \sum_j c_j \mu_j} with
#' one pooled error term. For equivalence, the power is the
#' probability that the \eqn{(1 - 2\alpha)} confidence interval for
#' \eqn{\psi} lies entirely inside \eqn{(-\delta_L, \delta_U)},
#' computed by numerical integration over the chi distribution of the
#' estimated error standard deviation; for noninferiority, the power
#' is a noncentral \emph{t} probability in closed form. This is the
#' contrast generalization of \code{\link{power_equivalence_md}}.
#'
#' @param c_weights The contrast weights. The weights must sum to zero
#'   with the positive weights summing to 1 and the negative weights
#'   to -1, so that the bounds are on the raw scale of the response.
#' @param n Sample sizes per group (if length 1, equal group sizes are
#'   assumed). Together with \code{c_weights}, \code{n} determines the
#'   standard error factor \eqn{\sqrt{\sum_j c_j^2 / n_j}}.
#' @param sigma The error standard deviation (the square root of the
#'   mean square error).
#' @param delta_lower,delta_upper Equivalence bounds on the raw scale
#'   of the response. Both must be positive; the equivalence region is
#'   \eqn{(-\delta_L, +\delta_U)}. If only \code{delta_upper} is
#'   supplied, the bounds are symmetric. Noninferiority uses
#'   \eqn{-\delta_L} alone.
#' @param true_psi The population value of the contrast at which the
#'   power is evaluated. Default \code{0}, the most favorable point
#'   for an equivalence declaration.
#' @param alpha_level One-sided significance level for each test. Default
#'   \code{0.05}.
#' @param side \code{"equivalence"} (default) for the TOST power, or
#'   \code{"noninferiority"} for the one-sided test against
#'   \eqn{-\delta_L}.
#' @param df_error The error degrees of freedom. Defaults to
#'   \eqn{N - J}; supply it directly when the error term comes from a
#'   model with more groups or additional predictors than the contrast
#'   involves (for example, a five-group model supplying the pooled
#'   error for a two-group contrast).
#'
#' @return A one-row \code{data.frame} with columns \code{term}
#'   (\code{"power"}) and \code{value} (the computed power, in
#'   \eqn{[0, 1]}).
#'
#' @details
#' \strong{Equivalence power.} Conditional on the estimated error
#' standard deviation \eqn{S}, the \eqn{(1 - 2\alpha)} CI fits inside
#' the bounds on a computable event, and the unconditional power
#' integrates that event over the scaled chi distribution of \eqn{S}
#' on \code{df_error} degrees of freedom. With \code{c_weights =
#' c(1, -1)} and equal \code{n}, the result reproduces
#' \code{\link{power_equivalence_md}} exactly.
#'
#' \strong{Noninferiority power.} The one-sided test rejects when
#' \eqn{t = (\hat\psi + \delta_L)/\mathrm{SE}(\hat\psi)} exceeds
#' \eqn{t_{1-\alpha,\nu}}, so the power is
#' \eqn{\Pr(T'_{\nu}(\lambda) > t_{1-\alpha,\nu})} with noncentrality
#' \eqn{\lambda = (\psi + \delta_L)/(\sigma \sqrt{\sum_j c_j^2/n_j})}.
#'
#' \strong{The feasibility condition.} If the expected half-width of
#' the CI is not smaller than the bounds allow, the equivalence power
#' is zero or near zero regardless of \code{true_psi}: an imprecise
#' design cannot declare equivalence even when the arms are truly
#' identical. Planning should target a half-width of about half the
#' bound; see \code{\link{ss_power_equivalence_c}} and
#' \code{\link{ss_aipe_c}}.
#'
#' @references
#' Chattopadhyay, B., Bandyopadhyay, T., Kelley, K., & Padalunkal,
#'   J. J. (2025). A sequential approach for noninferiority or
#'   equivalence of a linear contrast under cost constraints.
#'   \emph{Psychological Methods, 30}(2), 425--439. \doi{10.1037/met0000570}
#'
#' Phillips, K. F. (1990). Power of the two one-sided tests procedure
#'   in bioequivalence. \emph{Journal of Pharmacokinetics and
#'   Biopharmaceutics, 18}(2), 139--144. \doi{10.1007/BF01063556}
#'
#' Schuirmann, D. J. (1987). A comparison of the two one-sided tests
#'   procedure and the power approach for assessing the equivalence of
#'   average bioavailability. \emph{Journal of Pharmacokinetics and
#'   Biopharmaceutics, 15}(6), 657--680.
#'
#' @seealso \code{\link{power_equivalence_md}},
#'   \code{\link{ss_power_equivalence_c}}, \code{\link{equivalence_c}},
#'   \code{\link{ss_aipe_c}}
#'
#' @examples
#' # 1. Two groups of 61 and 113 sharing a five-group pooled error term
#' #    (so df_error = 404 - 5 = 399), bounds of 5 raw-scale points:
#' #    the design's probability of declaring equivalence when the
#' #    groups are truly identical.
#' power_equivalence_c(c_weights = c(1, -1), n = c(61, 113),
#'                     sigma = 15.67, delta_upper = 5,
#'                     true_psi = 0, df_error = 399)
#'
#' # 2. The same design's noninferiority power at the same point.
#' power_equivalence_c(c_weights = c(1, -1), n = c(61, 113),
#'                     sigma = 15.67, delta_upper = 5,
#'                     true_psi = 0, df_error = 399,
#'                     side = "noninferiority")
#'
#' # 3. Agreement with power_equivalence_md() in the two-group case
#' #    (Phillips, 1990, Table 1: expected 0.8029678).
#' power_equivalence_c(c_weights = c(1, -1), n = 24, sigma = 0.20,
#'                     delta_lower = 0.2, delta_upper = 0.2,
#'                     true_psi = 0.05, df_error = 22)
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @keywords design
#'
#' @family equivalence testing
#'
#' @export

power_equivalence_c <- function(c_weights, n, sigma,
                                delta_lower = NULL, delta_upper = NULL,
                                true_psi = 0, alpha_level = 0.05,
                                side = c("equivalence", "noninferiority"),
                                df_error = NULL) {
  side <- match.arg(side)

  if (is.null(delta_upper))
    stop("'delta_upper' must be specified (the upper equivalence bound).")
  if (!is.numeric(delta_upper) || length(delta_upper) != 1L || delta_upper <= 0)
    stop("'delta_upper' must be a single positive number.")
  if (is.null(delta_lower)) delta_lower <- delta_upper
  if (!is.numeric(delta_lower) || length(delta_lower) != 1L || delta_lower <= 0)
    stop("'delta_lower' must be a single positive number.")
  if (!is.numeric(alpha_level) || length(alpha_level) != 1L || alpha_level <= 0 || alpha_level >= 0.5)
    stop("'alpha_level' must be a single number in (0, 0.5).")
  if (!is.numeric(sigma) || length(sigma) != 1L || sigma <= 0)
    stop("'sigma' must be a single positive number.")
  if (!is.numeric(true_psi) || length(true_psi) != 1L)
    stop("'true_psi' must be a single numeric value.")

  if (length(n) == 1L) n <- rep(n, length(c_weights))
  if (length(n) != length(c_weights))
    stop("The lengths of 'n' and 'c_weights' differ, which should not be the case.")
  if (any(n < 2)) stop("Each group needs at least 2 observations.")
  if (!identical(round(sum(c_weights), 5), 0))
    stop("The sum of the contrast weights ('c_weights') should equal zero.")
  pos <- sum(c_weights[c_weights > 0])
  neg <- sum(c_weights[c_weights < 0])
  if (!isTRUE(all.equal(pos, 1)) || !isTRUE(all.equal(neg, -1)))
    stop("The positive weights must sum to 1 and the negative weights to ",
         "-1, so that the bounds are on the raw scale of the response.")

  k  <- sqrt(sum(c_weights^2 / n))
  nu <- if (is.null(df_error)) sum(n) - length(c_weights) else df_error
  if (!is.numeric(nu) || length(nu) != 1L || nu <= 0)
    stop("The error degrees of freedom must be a single positive number.")

  power <- .power_equivalence_c_fast(k = k, nu = nu, sigma = sigma,
                                     delta_lower = delta_lower,
                                     delta_upper = delta_upper,
                                     true_psi = true_psi, alpha_level = alpha_level,
                                     side = side)

  .as_dmar_tbl(data.frame(term = "power", value = power,
                          stringsAsFactors = FALSE))
}


# --- Fast path: bare-numeric power for the search loops -------------------
# Skips the argument validation and data.frame construction of the public
# power_equivalence_c(); consumed by ss_power_equivalence_c()'s sample size
# search, which evaluates the power at many candidate n per invocation.
# 'k' is the standard error factor sqrt(sum(c_weights^2 / n)).
.power_equivalence_c_fast <- function(k, nu, sigma, delta_lower, delta_upper,
                                      true_psi, alpha_level, side) {
  if (side == "noninferiority") {
    ncp <- (true_psi + delta_lower) / (sigma * k)
    return(stats::pt(stats::qt(1 - alpha_level, nu), df = nu, ncp = ncp,
                     lower.tail = FALSE))
  }

  # Equivalence: integrate the conditional probability that the interval fits
  # inside the bounds over the sampling distribution of the estimated error
  # standard deviation. The bounds enter only through their distance from
  # 'true_psi' in standard errors, which is the unit-free scale the shared
  # quadrature works on. That quadrature (.power_equivalence_tost_quad(), in
  # power_equivalence_md.R) documents the integrand and explains why its
  # limits are quantiles of the sampling distribution of the error standard
  # deviation rather than the raw interval that deviation must fall in.
  one_over_se <- 1 / (sigma * k)

  .power_equivalence_tost_quad(
    power_t = stats::qt(1 - alpha_level, df = nu),
    z_lower = (-delta_lower - true_psi) * one_over_se,
    z_upper = ( delta_upper - true_psi) * one_over_se,
    nu      = nu
  )
}
