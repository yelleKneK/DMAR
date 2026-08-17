#' Sample Size Planning for Power for the Indirect (Mediation) Effect
#'
#' Power, or the sample size required for a desired power, for the test of
#' the indirect effect \eqn{a b} in the simple mediation model, with paths
#' specified in standardized metric (unit-variance \eqn{X}, \eqn{M}, and
#' \eqn{Y}). The default test is joint significance (the indirect effect is
#' declared when both \eqn{\hat a} and \eqn{\hat b} are individually
#' significant), which tracks the resampling tests' power closely and far
#' exceeds the Sobel test in small samples (Fritz & MacKinnon, 2007); the
#' Sobel normal-theory test is available for comparison. This is the power
#' counterpart of the accuracy in parameter estimation (AIPE) planner
#' \code{\link{ss_aipe_indirect_effect}}, and the planning complement of
#' the analysis function \code{\link{mediate}}.
#'
#' @param a Standardized \eqn{X \to M} path.
#' @param b Standardized \eqn{M \to Y} path, holding \eqn{X}.
#' @param c_prime Standardized direct effect of \eqn{X} on \eqn{Y} holding
#'   \eqn{M}. Defaults to 0; it enters only through the residual variance
#'   of \eqn{Y}.
#' @param desired_power Desired power; supply this to solve for \eqn{N}.
#' @param N Total sample size; supply this to evaluate the realized power.
#'   Specify exactly one of \code{desired_power} and \code{N}. This is the
#'   total sample size.
#' @param alpha_level Two-sided Type I error rate for each component test.
#'   Defaults to 0.05.
#' @param method \code{"joint_significance"} (default) or \code{"sobel"}.
#'
#' @details
#' With unit-variance variables, the large-sample standard errors are
#' \eqn{\mathrm{se}_a = \sqrt{(1 - a^2)/N}} and
#' \eqn{\mathrm{se}_b = \sqrt{\sigma^2_{e_Y} / [N (1 - a^2)]}} with
#' \eqn{\sigma^2_{e_Y} = 1 - (b^2 + c'^2 + 2 a b c')}. Because \eqn{\hat a}
#' and \eqn{\hat b} are asymptotically independent in this model, the joint
#' significance power is the product of the two component powers; the Sobel
#' power refers \eqn{ab / \mathrm{se}_{ab}} (first-order delta method,
#' via the same variance as \code{\link{var_indirect_effect}}) to the
#' normal. The joint-significance component powers use the exact noncentral
#' \emph{t} (with \eqn{n - 2} and \eqn{n - 3} degrees of freedom), so its
#' only approximations are the population standard errors and component
#' independence; the tests validate the result against raw-data
#' simulation. A specified parameter combination must be
#' admissible (positive residual variances), or the function stops.
#'
#' @return A tidy \code{data.frame} with
#'   \code{necessary_N} (or \code{specified_N}), \code{actual_power} (the
#'   power to detect the indirect effect, the quantity the sample size is
#'   planned against), the component powers (\code{power_a}, \code{power_b};
#'   \code{NA} for the Sobel method), the paths (\code{a}, \code{b},
#'   \code{c_prime}), the implied \code{indirect_effect}, and
#'   \code{alpha_level}. The method is recorded in the \code{"method"}
#'   attribute. The result carries the \code{dmar_ss_power} class, so
#'   \code{\link[generics]{tidy}} reports the sample size and the power to
#'   detect the indirect effect, and \code{\link[generics]{glance}} adds the
#'   component powers and the planning inputs.
#'
#' @references
#' Fritz, M. S., & MacKinnon, D. P. (2007). Required sample size to detect
#'   the mediated effect. \emph{Psychological Science, 18}(3), 233--239.
#'   \doi{10.1111/j.1467-9280.2007.01882.x}
#'
#' MacKinnon, D. P., Lockwood, C. M., Hoffman, J. M., West, S. G., &
#'   Sheets, V. (2002). A comparison of methods to test mediation and
#'   other intervening variable effects. \emph{Psychological Methods,
#'   7}(1), 83--104. \doi{10.1037/1082-989X.7.1.83}
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso \code{\link{mediate}} to analyze the study this plans;
#'   \code{\link{ss_aipe_indirect_effect}} to plan for confidence interval
#'   width instead of detection; \code{\link{design_consequences}} for what
#'   the chosen design delivers.
#'
#' @family mediation
#' @family sample size for power
#'
#' @keywords design
#'
#' @examples
#' # Fritz and MacKinnon's (2007) running scenario (a = b = .39). The
#' # joint significance approximation returns necessary_N = 65; raw-data
#' # simulation puts the power at N = 65 nearer .77 and reaches .80 near
#' # N = 70, which is why Fritz and MacKinnon's simulation-based table
#' # reports a somewhat larger requirement.
#' ss_power_indirect_effect(a = .39, b = .39, desired_power = .80)
#'
#' # A near-zero a path against a larger b: the weak link drives the requirement.
#' ss_power_indirect_effect(a = .14, b = .39, desired_power = .80)
#'
#' # Realized power at a given N, and the Sobel comparison (always lower).
#' ss_power_indirect_effect(a = .39, b = .39, N = 75)
#' ss_power_indirect_effect(a = .39, b = .39, N = 75, method = "sobel")
#'
#' @export
#' @importFrom stats pnorm qnorm
ss_power_indirect_effect <- function(a, b, c_prime = 0,
                                     desired_power = NULL, N = NULL,
                                     alpha_level = 0.05,
                                     method = c("joint_significance",
                                                "sobel")) {
  method <- match.arg(method)
  for (nm in c("a", "b", "c_prime")) {
    val <- get(nm)
    if (!is.numeric(val) || length(val) != 1L || is.na(val) ||
        abs(val) >= 1) {
      stop(sprintf("'%s' must be a single standardized path in (-1, 1).",
                   nm), call. = FALSE)
    }
  }
  if (a == 0 || b == 0) {
    stop("'a' and 'b' must be nonzero: a zero path has no indirect effect ",
         "to detect.", call. = FALSE)
  }
  if (is.null(desired_power) == is.null(N)) {
    stop("Specify exactly one of 'desired_power' or 'N'.", call. = FALSE)
  }
  if (!is.null(desired_power) &&
      (!is.numeric(desired_power) || length(desired_power) != 1L ||
       is.na(desired_power) || desired_power <= 0 || desired_power >= 1)) {
    stop("'desired_power' must be a single number in (0, 1).",
         call. = FALSE)
  }
  if (!is.null(N) && (!is.numeric(N) || length(N) != 1L || is.na(N) ||
                      N < 10 || N != round(N))) {
    stop("'N' must be a single integer of at least 10.", call. = FALSE)
  }
  if (!is.numeric(alpha_level) || length(alpha_level) != 1L ||
      is.na(alpha_level) || alpha_level <= 0 || alpha_level >= 1) {
    stop("'alpha_level' must be a single number in (0, 1).", call. = FALSE)
  }

  resid_M <- 1 - a^2
  resid_Y <- 1 - (b^2 + c_prime^2 + 2 * a * b * c_prime)
  if (resid_Y <= 0) {
    stop("The supplied standardized paths imply a non-positive residual ",
         "variance for the outcome; the combination is not admissible.",
         call. = FALSE)
  }
  power_at <- function(n) {
    se_a <- sqrt(resid_M / n)
    se_b <- sqrt(resid_Y / (n * resid_M))
    if (method == "joint_significance") {
      # Each component is a regression t test; the noncentral t gives the
      # exact component power (df = n - 2 for the M model, n - 3 for the Y
      # model with X and M), which a normal approximation overstates in
      # small samples.
      pw_t <- function(ncp, df) {
        crit <- stats::qt(1 - alpha_level / 2, df)
        stats::pt(crit, df, ncp = ncp, lower.tail = FALSE) +
          stats::pt(-crit, df, ncp = ncp)
      }
      pa <- pw_t(abs(a) / se_a, n - 2)
      pb <- pw_t(abs(b) / se_b, n - 3)
      c(pa * pb, pa, pb)
    } else {
      z_crit <- qnorm(1 - alpha_level / 2)
      se_ab <- sqrt(b^2 * se_a^2 + a^2 * se_b^2)
      z <- abs(a * b) / se_ab
      c(pnorm(z - z_crit) + pnorm(-z - z_crit), NA_real_, NA_real_)
    }
  }

  if (is.null(N)) {
    # Fail fast when the target cannot be reached at any realistic size (for
    # example a near-zero indirect effect) rather than incrementing the sample
    # size without bound. Power is monotone in N, so if the bound does not
    # reach the target no smaller size does either.
    n_max <- 1e7
    if (power_at(n_max)[1] < desired_power)
      stop("Could not reach 'desired_power' at a sample size up to ",
           format(n_max, scientific = FALSE), "; the indirect effect may be ",
           "too small to detect at this power.", call. = FALSE)
    n_i <- 10L
    while (power_at(n_i)[1] < desired_power) {
      n_i <- n_i + 1L
      if (n_i > n_max)
        stop("The search for 'desired_power' did not converge below ",
             format(n_max, scientific = FALSE), ".", call. = FALSE)
    }
    label <- "necessary_N"
  } else {
    n_i <- as.integer(N)
    label <- "specified_N"
  }
  pw <- power_at(n_i)

  out <- data.frame(
    term  = c(label, "actual_power", "power_a", "power_b",
              "a", "b", "c_prime", "indirect_effect", "alpha_level"),
    value = c(n_i, pw[1], pw[2], pw[3], a, b, c_prime, a * b, alpha_level),
    stringsAsFactors = FALSE
  )
  # The design's power is the joint power to detect the indirect effect
  # (actual_power); tidy() reports it with the sample size, and glance() surfaces
  # the diagnostic component powers power_a / power_b as extra columns.
  class(out) <- c("dmar_ss_power", class(out))
  out <- .as_dmar_tbl(out)
  attr(out, "method") <- method
  out
}
