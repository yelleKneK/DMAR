# Sample size for the TOST or noninferiority test on a linear contrast.
#' Sample Size for Equivalence or Noninferiority of a Linear Contrast
#'
#' Computes the smallest per-group sample size at which the two
#' one-sided tests procedure (Schuirmann, 1987) for a linear contrast
#' \eqn{\psi = \sum_j c_j \mu_j}, or the companion one-sided
#' noninferiority test, attains a desired power. Power is computed
#' exactly through \code{\link{power_equivalence_c}}. This is the
#' declaration-probability route to planning; the accuracy in
#' parameter estimation (AIPE) route, which targets the confidence
#' interval width directly, is \code{\link{ss_aipe_c}}, and the two
#' answer the same question whenever \code{true_psi = 0} and the width
#' target is calibrated to the bounds.
#'
#' @param c_weights The contrast weights. The weights must sum to zero
#'   with the positive weights summing to 1 and the negative weights
#'   to -1, so that the bounds are on the raw scale of the response.
#' @param sigma The anticipated error standard deviation (the square
#'   root of the mean square error).
#' @param delta_lower,delta_upper Equivalence bounds on the raw scale
#'   of the response. Both must be positive; the equivalence region is
#'   \eqn{(-\delta_L, +\delta_U)}. If only \code{delta_upper} is
#'   supplied, the bounds are symmetric. Noninferiority uses
#'   \eqn{-\delta_L} alone.
#' @param true_psi The population contrast the design should be able
#'   to detect as equivalent (or noninferior). Default \code{0}. For
#'   equivalence it must lie strictly inside the bounds; for
#'   noninferiority, strictly above \eqn{-\delta_L}. The farther
#'   \code{true_psi} sits from the center of the region, the larger
#'   the required sample size.
#' @param desired_power The target probability of declaring
#'   equivalence (or noninferiority) at \code{true_psi}. Default
#'   \code{0.85}, matching \code{\link{ss_power_contrast}}.
#' @param alpha_level One-sided significance level for each test. Default
#'   \code{0.05}.
#' @param side \code{"equivalence"} (default) or
#'   \code{"noninferiority"}.
#'
#' @return A \code{data.frame} with rows \code{necessary_n_per_group}
#'   (the recommended sample size for each of the \eqn{J} groups named
#'   by \code{c_weights}), \code{total_N} (the implied total,
#'   \eqn{J \times n}), and \code{actual_power} (the exact power
#'   achieved at the recommendation). The result carries the
#'   \code{dmar_ss_power} class, so \code{\link[generics]{tidy}} and
#'   \code{\link[generics]{glance}} summarize it in broom convention.
#'
#' @details
#' \strong{Design.} Planning assumes equal allocation across the
#' \eqn{J} groups named by \code{c_weights} and a pooled error term on
#' \eqn{N - J} degrees of freedom. Groups with zero weight still
#' contribute error degrees of freedom, which is why they belong in
#' \code{c_weights} when the fitted model will include them.
#'
#' \strong{The search.} The function starts from the normal-theory
#' approximation and moves to the smallest integer \eqn{n} whose exact
#' power reaches \code{desired_power}. Power is monotone in \eqn{n}
#' once the design is feasible, so the search is a short walk.
#'
#' \strong{Relation to the half-width rule.} With symmetric bounds
#' \eqn{\pm\delta}, \code{true_psi = 0}, and \eqn{\alpha = .05},
#' targeting a 90\% CI half-width of \eqn{\delta/2} yields a
#' declaration probability of about .90; this function makes the
#' probability the target directly rather than through the width.
#'
#' @references
#' Chattopadhyay, B., Bandyopadhyay, T., Kelley, K., & Padalunkal,
#'   J. J. (2025). A sequential approach for noninferiority or
#'   equivalence of a linear contrast under cost constraints.
#'   \emph{Psychological Methods, 30}(2), 425--439. \doi{10.1037/met0000570}
#'
#' Maxwell, S. E., Kelley, K., & Rausch, J. R. (2008). Sample size
#'   planning for statistical power and accuracy in parameter
#'   estimation. \emph{Annual Review of Psychology, 59}, 537--563.
#'   \doi{10.1146/annurev.psych.59.103006.093735}
#'
#' Schuirmann, D. J. (1987). A comparison of the two one-sided tests
#'   procedure and the power approach for assessing the equivalence of
#'   average bioavailability. \emph{Journal of Pharmacokinetics and
#'   Biopharmaceutics, 15}(6), 657--680.
#'
#' @seealso \code{\link{power_equivalence_c}}, \code{\link{ss_aipe_c}},
#'   \code{\link{ss_power_contrast}}, \code{\link{equivalence_c}}
#'
#' @examples
#' # 1. Two-group equivalence with bounds of 5 raw-scale points and an
#' #    anticipated error SD of 15.67: n per group for 90% power at
#' #    true equivalence.
#' ss_power_equivalence_c(c_weights = c(1, -1), sigma = 15.67,
#'                        delta_upper = 5, desired_power = 0.90)
#'
#' # 2. Noninferiority is cheaper than equivalence at the same bound.
#' ss_power_equivalence_c(c_weights = c(1, -1), sigma = 15.67,
#'                        delta_upper = 5, desired_power = 0.90,
#'                        side = "noninferiority")
#'
#' # 3. A true contrast off center raises the requirement.
#' ss_power_equivalence_c(c_weights = c(1, -1), sigma = 15.67,
#'                        delta_upper = 5, true_psi = 2,
#'                        desired_power = 0.90)
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @keywords design
#'
#' @family equivalence testing
#'
#' @family sample size for power
#'
#' @export

ss_power_equivalence_c <- function(c_weights, sigma,
                                   delta_lower = NULL, delta_upper = NULL,
                                   true_psi = 0, desired_power = 0.85,
                                   alpha_level = 0.05,
                                   side = c("equivalence", "noninferiority")) {
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
  if (!is.numeric(desired_power) || length(desired_power) != 1L ||
      desired_power <= 0 || desired_power >= 1)
    stop("'desired_power' must be a single number in (0, 1).")
  if (!is.numeric(true_psi) || length(true_psi) != 1L)
    stop("'true_psi' must be a single numeric value.")
  if (!identical(round(sum(c_weights), 5), 0))
    stop("The sum of the contrast weights ('c_weights') should equal zero.")
  pos <- sum(c_weights[c_weights > 0])
  neg <- sum(c_weights[c_weights < 0])
  if (!isTRUE(all.equal(pos, 1)) || !isTRUE(all.equal(neg, -1)))
    stop("The positive weights must sum to 1 and the negative weights to ",
         "-1, so that the bounds are on the raw scale of the response.")

  if (side == "equivalence" &&
      (true_psi <= -delta_lower || true_psi >= delta_upper))
    stop("For equivalence, 'true_psi' must lie strictly inside ",
         "(-delta_lower, delta_upper); no sample size can reach the ",
         "desired power at a population contrast outside the bounds.")
  if (side == "noninferiority" && true_psi <= -delta_lower)
    stop("For noninferiority, 'true_psi' must exceed -delta_lower; no ",
         "sample size can reach the desired power at a population ",
         "contrast at or below the bound.")

  J   <- length(c_weights)
  ssq <- sum(c_weights^2)

  # Normal-theory start: half-width h = z * sigma * sqrt(ssq / n) must fit
  # between true_psi and the nearer bound with z_power to spare. The exact
  # search below corrects the approximation.
  z_a <- stats::qnorm(1 - alpha_level)
  z_p <- stats::qnorm(desired_power)
  gap <- if (side == "noninferiority") true_psi + delta_lower else
    min(true_psi + delta_lower, delta_upper - true_psi)
  n <- max(2, ceiling(((z_a + z_p) * sigma)^2 * ssq / gap^2))

  power_at <- function(n) {
    .power_equivalence_c_fast(k = sqrt(ssq / n), nu = n * J - J,
                              sigma = sigma,
                              delta_lower = delta_lower,
                              delta_upper = delta_upper,
                              true_psi = true_psi, alpha_level = alpha_level,
                              side = side)
  }

  # Walk down while the target still holds, then up until it is reached.
  # Power is monotone in n once feasible, so both walks terminate.
  while (n > 2 && power_at(n - 1) >= desired_power) n <- n - 1
  n_max <- n * 1000
  while (power_at(n) < desired_power) {
    n <- n + 1
    if (n > n_max)
      stop("The search did not reach the desired power; check that the ",
           "bounds, 'sigma', and 'true_psi' describe a feasible design.")
  }

  out <- data.frame(
    term  = c("necessary_n_per_group", "total_N", "actual_power"),
    value = c(n, n * J, power_at(n)),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  # The size and power rows use the shared ss_power_* schema names, so the
  # dmar_ss_power class gives tidy()/glance() the right rows with no change to
  # the shared tidy-verb lookup vectors (see R/dmar_tidiers.R).
  class(out) <- c("dmar_ss_power", class(out))
  .as_dmar_tbl(out)
}
