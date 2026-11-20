# --- Internal integrand for the TOST power calculation -------------------
# Numerical integration via stats::integrate() requires a function that takes
# a numeric vector and returns a numeric vector of the same length. The
# tidy/data-frame-returning public helper (power_density_equivalence_md())
# would not satisfy that contract, so we expose this private numeric version
# for the integrator and re-use it inside the public density function.
.power_equivalence_integrand <- function(sd_values, alpha_level, lower_bound,
                                         upper_bound, true_effect, true_sd,
                                         n_per_treatment, df_error) {
  power_t      <- stats::qt(1 - alpha_level, df = df_error)
  a            <- true_sd / sqrt(df_error)  # Chi-distribution scale constant
  # lgamma, not log(gamma(.)): gamma(df_error / 2) overflows to Inf for
  # df_error beyond about 344 (ordinary residual df), which would drive
  # power_const to -Inf and the integrand, and hence the power, silently to 0.
  power_const  <- -(df_error / 2 - 1) * log(2) - lgamma(df_error / 2) - log(a)
  half_width   <- power_t * sqrt(2 / n_per_treatment)
  one_over_se  <- 1 / (true_sd * sqrt(2 / n_per_treatment))

  d1   <- ( sd_values * half_width + (lower_bound - true_effect)) * one_over_se
  d2   <- (-sd_values * half_width + (upper_bound - true_effect)) * one_over_se
  phi  <- stats::pnorm(d2) - stats::pnorm(d1)

  phi * exp(power_const + (df_error - 1) * log(sd_values / a) -
              0.5 * (sd_values / a)^2)
}


# --- Shared scale-free quadrature for the TOST power ----------------------
# Consumed by power_equivalence_md() and, for the contrast generalization, by
# .power_equivalence_c_fast() in power_equivalence_c.R. Both reduce to the
# same one-dimensional problem once the equivalence limits are expressed in
# standard error units.
#
# Write r = S / true_sd for the realized error standard deviation on a unit-free
# scale, t for the one-sided critical value, and z_lower, z_upper for the
# equivalence limits measured in standard errors from the true difference.
# Conditional on r, the (1 - 2 * alpha_level) confidence interval falls inside the
# equivalence interval with probability
#
#   phi(r) = Phi(z_upper - r * t) - Phi(z_lower + r * t),
#
# which reaches zero at r_max = (z_upper - z_lower) / (2 * t), the point where
# the interval is exactly as wide as the equivalence interval, and is zero
# beyond it. The power is the expectation of phi(r) over
# r = sqrt(X / df_error) with X a chi square variate on df_error degrees of
# freedom, so the density of r is
# 2 * df_error * r * dchisq(df_error * r^2, df_error).
#
# Two properties of this form matter. First, every quantity is a ratio, so the
# result cannot depend on the units of the response. Second, the integration
# runs over the extreme quantiles of the distribution of r rather than over
# the raw interval (0, r_max]. That second point is what makes the quadrature
# trustworthy: for a precise design r_max sits far out in the right tail, and
# (0, r_max] can be hundreds of times wider than the region where r has any
# mass at all. The initial Gauss-Kronrod nodes then all land where the density
# has underflowed, and QUADPACK returns message "OK" with a negligible error
# estimate on an essentially zero integral, which reported a power of about
# 1e-13 for designs whose true power is 1. Cutting the domain at the quantiles
# of r instead bounds its width by the spread of r itself (about 18 standard
# deviations of r once df_error is large), so the nodes always sample the mass, and
# the discarded tails carry at most 2e-18 of probability.
.power_equivalence_tost_quad <- function(power_t, z_lower, z_upper, df_error) {
  r_max <- (z_upper - z_lower) / (2 * power_t)
  if (!is.finite(r_max) || r_max <= 0) return(0)

  tail_p <- 1e-18
  lower  <- sqrt(stats::qchisq(tail_p, df = df_error) / df_error)
  upper  <- sqrt(stats::qchisq(tail_p, df = df_error, lower.tail = FALSE) /
                   df_error)
  if (r_max <= lower) {
    # The interval can fit only for an error standard deviation smaller than
    # essentially any that will be observed: the power is below 1e-18, and the
    # whole of it lives in the left tail.
    lower <- 0
    upper <- r_max
  } else {
    upper <- min(upper, r_max)
  }
  if (!(upper > lower)) return(0)

  integrand <- function(r) {
    # pmax() at zero: past r_max the two normal tails cross and their
    # difference turns negative, whereas the conditional probability of the
    # interval fitting is exactly zero there.
    pmax(stats::pnorm(z_upper - r * power_t) -
           stats::pnorm(z_lower + r * power_t), 0) *
      2 * df_error * r * stats::dchisq(df_error * r^2, df = df_error)
  }

  power <- stats::integrate(
    integrand,
    lower         = lower, upper = upper,
    subdivisions  = 10000L,
    rel.tol       = 1e-12,
    abs.tol       = 1e-14,
    stop.on.error = TRUE
  )$value

  # A quadrature result is not automatically a probability. The documented
  # return is a power in [0, 1] and the raw value can land a few units in the
  # last place outside it (values such as 1.000000000013 were reachable), so
  # the interval is enforced rather than assumed.
  min(max(power, 0), 1)
}


#' Power of the Two One-Sided Tests Procedure (TOST) for Equivalence
#'
#' Computes the power of the Schuirmann (1987) two one-sided tests procedure
#', the probability that a \eqn{(1 - 2\alpha)} confidence interval for the
#' mean difference (or ratio, on the log scale) lies entirely within the
#' equivalence interval \eqn{[\theta_1, \theta_2]} (\code{lower_bound} and
#' \code{upper_bound}), by numerical integration over the chi distribution
#' of the sample standard deviation.
#'
#' @param alpha_level Type I error rate for each of the two one-sided tests
#'   (typically \code{0.05}). The full equivalence test uses a
#'   \eqn{(1 - 2\alpha)} confidence interval.
#' @param logscale Logical. If \code{TRUE}, treatment means are compared on
#'   the logarithmic scale; \code{lower_bound}, \code{upper_bound}, and
#'   \code{true_effect} are then given as ratios of test to reference on the
#'   original scale (untransformed) and are log-transformed internally. If
#'   \code{FALSE}, they are differences of the means.
#' @param lower_bound Lower limit of the equivalence interval on the
#'   original scale: a ratio when \code{logscale = TRUE} (logged
#'   internally), a difference otherwise.
#' @param upper_bound Upper limit of the equivalence interval on the
#'   original scale: a ratio when \code{logscale = TRUE} (logged
#'   internally), a difference otherwise.
#' @param true_effect The true (population) effect: the ratio of the
#'   treatment means on the original scale when \code{logscale = TRUE}
#'   (logged internally), the difference in treatment means otherwise.
#' @param true_sd \eqn{\sqrt{\mathrm{error\ variance}}}; root-MSE from an
#'   ANOVA. On the log scale, this is the coefficient of variation.
#' @param n_per_treatment Number of subjects who receive each treatment:
#'   the size of each group in a parallel design, and the number of
#'   subjects in a crossover design, where every subject receives both
#'   treatments. It is the same quantity in both designs, the \eqn{n} of
#'   the standard error \eqn{\sigma \sqrt{2 / n}} of the estimated
#'   difference, so in a parallel design it is the size of one group and
#'   never the total. Must be greater than 1.
#' @param df_error Error degrees of freedom associated with \code{true_sd}.
#'
#' @return A one-row \code{data.frame} with columns \code{term} (\code{"power"})
#'   and \code{value} (the computed power, in \eqn{[0, 1]}).
#'
#' @details
#' The computation conditions on the error standard deviation the study will
#' actually observe. Given that value, whether the confidence interval fits
#' inside the equivalence interval is an ordinary normal probability, and the
#' power is that probability averaged over the chi distribution the error
#' standard deviation follows on \code{df_error} degrees of freedom. The
#' averaging is carried out on a unit-free scale, with the equivalence limits
#' expressed in standard errors and the error standard deviation as a
#' multiple of \code{true_sd}, so the power depends on the design rather than
#' on the units of the response: with \code{logscale = FALSE}, multiplying
#' \code{lower_bound}, \code{upper_bound}, \code{true_effect}, and
#' \code{true_sd} by a common factor leaves the answer unchanged.
#'
#' \strong{The scale of the bounds and the true effect.} This function and
#' \code{\link{power_equivalence_md_plot}} take \code{lower_bound},
#' \code{upper_bound}, and \code{true_effect} on the original scale: ratios
#' of test to reference with \code{logscale = TRUE}, which are logged
#' before the computation, and differences with \code{logscale = FALSE},
#' which are used as given. \code{\link{power_density_equivalence_md}} takes
#' them on the analysis scale and transforms nothing, so
#' \code{power_equivalence_md(logscale = TRUE, lower_bound = .8, ...)}
#' corresponds to \code{power_density_equivalence_md(lower_bound = log(.8),
#' ...)}.
#'
#' For Phillips's (1990) original example (regular-scale two-period crossover
#' with \eqn{\theta_1 = -0.2}, \eqn{\theta_2 = 0.2}, CV \eqn{= 0.20},
#' \eqn{\delta = 0.05}, \eqn{n = 24}, \eqn{\nu = 22}, which are
#' \code{lower_bound}, \code{upper_bound}, \code{true_sd},
#' \code{true_effect}, \code{n_per_treatment}, and \code{df_error}), this
#' function reproduces the published value of \eqn{0.8029678} (Phillips,
#' 1990, Table 1, 5th row, 5th column).
#'
#' @references
#' Diletti, E., Hauschke, D., & Steinijans, V. W. (1991). Sample size
#'   determination of bioequivalence assessment by means of confidence
#'   intervals. \emph{International Journal of Clinical Pharmacology, Therapy
#'   and Toxicology, 29}(1), 1--8.
#'
#' Kelley, K. (2007a). Confidence intervals for standardized effect
#'   sizes: Theory, application, and implementation. \emph{Journal of
#'   Statistical Software, 20}(8), 1--24. \doi{10.18637/jss.v020.i08}
#'
#' Kelley, K. (2007b). Methods for the behavioral, educational, and
#'   social sciences: An R package. \emph{Behavior Research Methods,
#'   39}(4), 979--984. \doi{10.3758/BF03192993}
#'
#' Phillips, K. F. (1990). Power of the two one-sided tests procedure in
#'   bioequivalence. \emph{Journal of Pharmacokinetics and Biopharmaceutics,
#'   18}(2), 137--144. \doi{10.1007/BF01063556}
#'
#' Schuirmann, D. J. (1987). A comparison of the two one-sided tests
#'   procedure and the power approach for assessing the equivalence of
#'   average bioavailability. \emph{Journal of Pharmacokinetics and
#'   Biopharmaceutics, 15}(6), 657--680.
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @note See the legacy \code{MBESS} package (Kelley, 2007a, 2007b) for
#'   additional details and discussion.
#'
#' @seealso \code{\link{power_equivalence_md_plot}},
#'   \code{\link{power_density_equivalence_md}}
#'
#' @examples
#' # Table 1 of Phillips, 1990, fifth row and fifth column, where the
#' # published power is 0.8029678.
#' power_equivalence_md(alpha_level = .05, logscale = FALSE,
#'                      lower_bound = -.2, upper_bound = .2,
#'                      true_effect = .05, true_sd = .20,
#'                      n_per_treatment = 24, df_error = 22)
#'
#' # Table 1 of Diletti et al., 1991, on the log scale, so the limits and
#' # the true difference are ratios of test to reference. The published
#' # power is 0.7922796.
#' power_equivalence_md(alpha_level = .05, logscale = TRUE,
#'                      lower_bound = .8, upper_bound = 1.25,
#'                      true_effect = 1.05, true_sd = .20,
#'                      n_per_treatment = 18, df_error = 16)
#'
#' @keywords design
#'
#' @family equivalence testing
#'
#' @export
#' @import stats
power_equivalence_md <- function(alpha_level, logscale, lower_bound,
                                 upper_bound, true_effect, true_sd,
                                 n_per_treatment, df_error) {
  if (!missing(alpha_level)) .dmar_check_alpha_level(alpha_level)
  .dmar_check_numbers("power_equivalence_md", c("lower_bound", "upper_bound", "true_sd", "n_per_treatment", "df_error"), required = c("lower_bound", "upper_bound", "true_sd", "n_per_treatment", "df_error"))

  if (!is.numeric(alpha_level) || length(alpha_level) != 1L || alpha_level <= 0 || alpha_level >= 0.5) {
    .dmar_refuse_argument(
      "'alpha_level' must be a single number in (0, 0.5).", arg = "alpha_level")
  }
  if (!is.logical(logscale) || length(logscale) != 1L) {
    .dmar_refuse_argument(
      "'logscale' must be a single logical value.", arg = "logscale")
  }
  if (true_sd <= 0) .dmar_refuse_argument(
    "'true_sd' must be positive.", arg = "true_sd")
  if (n_per_treatment <= 1) .dmar_refuse_argument(
    "'n_per_treatment' must be greater than 1.", arg = "n_per_treatment")
  if (df_error <= 0) .dmar_refuse_argument(
    "'df_error' must be positive.", arg = "df_error")
  if (upper_bound <= lower_bound) {
    .dmar_refuse_argument(
      "'upper_bound' must be greater than 'lower_bound'.",
      arg = c("upper_bound", "lower_bound"))
  }

  # theta1, theta2, and diff are the bounds and the true effect on the
  # analysis scale: the logarithms of the ratios given when logscale = TRUE,
  # the differences given otherwise.
  if (logscale) {
    if (lower_bound <= 0 || upper_bound <= 0 || true_effect <= 0) {
      .dmar_refuse_argument(
        "On the log scale, 'lower_bound', 'upper_bound', and 'true_effect' must all be positive.",
        arg = c("lower_bound", "upper_bound", "true_effect"))
    }
    theta1 <- log(lower_bound)
    theta2 <- log(upper_bound)
    diff   <- log(true_effect)
  } else {
    theta1 <- lower_bound
    theta2 <- upper_bound
    diff   <- true_effect
  }

  # The power depends on the response scale only through the two equivalence
  # limits measured in standard errors from the true difference, so the
  # integration is carried out on that unit-free scale. See
  # .power_equivalence_tost_quad() above for the integrand and for why the
  # limits of integration are quantiles of the sampling distribution of the
  # error standard deviation rather than the raw interval it must fall in.
  power_t     <- stats::qt(1 - alpha_level, df = df_error)
  one_over_se <- 1 / (true_sd * sqrt(2 / n_per_treatment))

  power <- .power_equivalence_tost_quad(
    power_t = power_t,
    z_lower = (theta1 - diff) * one_over_se,
    z_upper = (theta2 - diff) * one_over_se,
    df_error = df_error
  )

  .as_dmar_tbl(data.frame(term = "power", value = power, stringsAsFactors = FALSE))
}


#' Density Underlying the TOST Power Calculation
#'
#' Evaluates the integrand whose integral over \eqn{(0, \mathrm{upper})}
#' yields the power of the Schuirmann (1987) two one-sided tests procedure;
#' see \code{\link{power_equivalence_md}}. Useful for plotting the integrand
#' and for diagnostic work.
#'
#' \strong{The analysis scale.} Unlike \code{\link{power_equivalence_md}},
#' which takes the bounds and the true effect on the original scale and logs
#' ratios when \code{logscale = TRUE}, this function has no
#' \code{logscale} argument and transforms nothing: \code{lower_bound},
#' \code{upper_bound}, and \code{true_effect} are on the scale of the
#' analysis, differences of the means, or, for means compared on the log
#' scale, logarithms of the ratios (\code{log(0.8)} for a ratio of 0.8).
#'
#' @param sd_values Numeric vector of \eqn{\sigma} values at which to
#'   evaluate the integrand.
#' @param alpha_level Type I error rate for each of the two one-sided tests.
#' @param lower_bound Lower limit of the equivalence interval on the
#'   analysis scale: a difference, or the logarithm of a ratio.
#' @param upper_bound Upper limit of the equivalence interval on the
#'   analysis scale: a difference, or the logarithm of a ratio.
#' @param true_effect The true (population) effect on the analysis scale:
#'   the difference in treatment means, or the logarithm of their ratio.
#' @param true_sd \eqn{\sqrt{\mathrm{error\ variance}}}.
#' @param n_per_treatment Number of subjects who receive each treatment
#'   (the size of each group in a parallel design; every subject in a
#'   crossover design), as in \code{\link{power_equivalence_md}}.
#' @param df_error Error degrees of freedom for \code{true_sd}.
#'
#' @return A \code{data.frame} with one row per supplied
#'   \code{sd_values}, and columns \code{sd_values} and
#'   \code{power_density}.
#'
#' @references
#' Kelley, K. (2007a). Confidence intervals for standardized effect
#'   sizes: Theory, application, and implementation. \emph{Journal of
#'   Statistical Software, 20}(8), 1--24. \doi{10.18637/jss.v020.i08}
#'
#' Kelley, K. (2007b). Methods for the behavioral, educational, and
#'   social sciences: An R package. \emph{Behavior Research Methods,
#'   39}(4), 979--984. \doi{10.3758/BF03192993}
#'
#' Phillips, K. F. (1990). Power of the two one-sided tests procedure in
#'   bioequivalence. \emph{Journal of Pharmacokinetics and Biopharmaceutics,
#'   18}(2), 137--144. \doi{10.1007/BF01063556}
#'
#' Schuirmann, D. J. (1987). A comparison of the two one-sided tests
#'   procedure and the power approach for assessing the equivalence of
#'   average bioavailability. \emph{Journal of Pharmacokinetics and
#'   Biopharmaceutics, 15}(6), 657--680.
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @note See the legacy \code{MBESS} package (Kelley, 2007a, 2007b) for
#'   additional details and discussion.
#'
#' @seealso \code{\link{power_equivalence_md}},
#'   \code{\link{power_equivalence_md_plot}}
#'
#' @examples
#' # Density at a single value of true_sd:
#' power_density_equivalence_md(sd_values = 0.10, alpha_level = .05,
#'                              lower_bound = -.2, upper_bound = .2,
#'                              true_effect = .05, true_sd = .20,
#'                              n_per_treatment = 24, df_error = 22)
#'
#' # Vectorized over a grid:
#' grid <- power_density_equivalence_md(
#'   sd_values = seq(0.01, 0.40, length.out = 50),
#'   alpha_level = .05, lower_bound = -.2, upper_bound = .2,
#'   true_effect = .05, true_sd = .20, n_per_treatment = 24, df_error = 22
#' )
#' head(grid)
#'
#' # On the log scale the bounds and the true effect are logarithms of the
#' # ratios: the ratios .8, 1.25, and 1.05 of power_equivalence_md(logscale =
#' # TRUE) are given here as their logarithms.
#' power_density_equivalence_md(sd_values = 0.10, alpha_level = .05,
#'                              lower_bound = log(.8), upper_bound = log(1.25),
#'                              true_effect = log(1.05), true_sd = .20,
#'                              n_per_treatment = 18, df_error = 16)
#'
#' @keywords design
#'
#' @family equivalence testing
#'
#' @export
#' @import stats
power_density_equivalence_md <- function(sd_values, alpha_level, lower_bound,
                                         upper_bound, true_effect, true_sd,
                                         n_per_treatment, df_error) {
  if (!missing(alpha_level)) .dmar_check_alpha_level(alpha_level)
  density <- .power_equivalence_integrand(
    sd_values = sd_values,
    alpha_level = alpha_level, lower_bound = lower_bound,
    upper_bound = upper_bound, true_effect = true_effect, true_sd = true_sd,
    n_per_treatment = n_per_treatment, df_error = df_error
  )
  data.frame(sd_values = sd_values, power_density = density,
             stringsAsFactors = FALSE)
}


#' Plot TOST Equivalence-Test Power Curves Over a Range of True Differences
#'
#' For each sample size in \code{n_per_treatment}, draws power as a function of the true
#' mean difference (or ratio, on the log scale), evaluated at 201 equally
#' spaced points across the equivalence interval. Returns a \pkg{ggplot2}
#' object; the underlying numerical grid is attached as
#' \code{attr(<plot>, "power_grid")}.
#'
#' @param alpha_level Type I error rate for each of the two one-sided tests.
#' @param logscale Logical. If \code{TRUE}, the means are compared on the
#'   logarithmic scale.
#' @param lower_bound Lower limit of the equivalence interval on the
#'   original scale, as in \code{\link{power_equivalence_md}}: a ratio when
#'   \code{logscale = TRUE}, a difference otherwise.
#' @param upper_bound Upper limit of the equivalence interval on the
#'   original scale: a ratio when \code{logscale = TRUE}, a difference
#'   otherwise.
#' @param true_sd \eqn{\sqrt{\mathrm{error\ variance}}}.
#' @param n_per_treatment Vector of sample sizes, each the number of
#'   subjects who receive each treatment (the size of each group in a
#'   parallel design; every subject in a crossover design), as in
#'   \code{\link{power_equivalence_md}}. One curve is drawn per element.
#' @param df_error Vector of error degrees of freedom for \code{true_sd},
#'   the same length as \code{n_per_treatment}.
#' @param title Optional plot title (default \code{"Power of TOST"}).
#' @param subtitle Optional subtitle (typically a reference like
#'   \code{"Phillips, Figure 3"}).
#'
#' @return A \code{ggplot} object. The 201-row power grid (column 1,
#'   \code{true_effect}: the true effect on the original scale, a ratio when
#'   \code{logscale = TRUE} and a difference otherwise, from
#'   \code{lower_bound} to \code{upper_bound}; remaining columns: power for
#'   each element of \code{n_per_treatment}, labeled \code{n=<size>}) is
#'   attached as \code{attr(<plot>, "power_grid")}.
#'
#' @references
#' Diletti, E., Hauschke, D., & Steinijans, V. W. (1991). Sample size
#'   determination of bioequivalence assessment by means of confidence
#'   intervals. \emph{International Journal of Clinical Pharmacology, Therapy
#'   and Toxicology, 29}(1), 1--8.
#'
#' Kelley, K. (2007a). Confidence intervals for standardized effect
#'   sizes: Theory, application, and implementation. \emph{Journal of
#'   Statistical Software, 20}(8), 1--24. \doi{10.18637/jss.v020.i08}
#'
#' Kelley, K. (2007b). Methods for the behavioral, educational, and
#'   social sciences: An R package. \emph{Behavior Research Methods,
#'   39}(4), 979--984. \doi{10.3758/BF03192993}
#'
#' Phillips, K. F. (1990). Power of the two one-sided tests procedure in
#'   bioequivalence. \emph{Journal of Pharmacokinetics and Biopharmaceutics,
#'   18}(2), 137--144. \doi{10.1007/BF01063556}
#'
#' Schuirmann, D. J. (1987). A comparison of the two one-sided tests
#'   procedure and the power approach for assessing the equivalence of
#'   average bioavailability. \emph{Journal of Pharmacokinetics and
#'   Biopharmaceutics, 15}(6), 657--680.
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @note See the legacy \code{MBESS} package (Kelley, 2007a, 2007b) for
#'   additional details and discussion.
#'
#' @seealso \code{\link{power_equivalence_md}},
#'   \code{\link{power_density_equivalence_md}}
#'
#' @examples
#' # One curve per sample size, showing power against the true mean
#' # difference. The seven sample sizes are the ones behind Figure 3 of
#' # Phillips, 1990, so the figure reproduces that one; every curve
#' # evaluates the power integral at 201 true differences.
#' n      <- c(9, 12, 18, 24, 30, 40, 60)
#' df_err <- c(7, 10, 16, 22, 28, 38, 58)
#' fig <- power_equivalence_md_plot(
#'   alpha_level = .05, logscale = FALSE,
#'   lower_bound = -.2, upper_bound = .2, true_sd = .20,
#'   n_per_treatment = n, df_error = df_err,
#'   subtitle = "Phillips Figure 3"
#' )
#' fig
#'
#' # The numbers behind the curves travel with the figure, so a particular
#' # power value can be read off rather than eyeballed. The first column,
#' # true_effect, is the true difference and the remaining columns give power,
#' # one column per sample size. Power is highest where the true difference
#' # is zero.
#' power_grid <- attr(fig, "power_grid")
#' power_grid[which.min(abs(power_grid[, "true_effect"])), ]
#'
#' # Figure 1c of Diletti et al., 1991, is the same idea on the log scale,
#' # where the equivalence limits are the 0.80 to 1.25 ratio bounds used
#' # in bioequivalence work, given as ratios.
#' n_d      <- c(8, 12, 18, 24, 30, 40, 60)
#' df_err_d <- c(6, 10, 16, 22, 28, 38, 58)
#' power_equivalence_md_plot(
#'   alpha_level = .05, logscale = TRUE,
#'   lower_bound = .8, upper_bound = 1.25, true_sd = .20,
#'   n_per_treatment = n_d, df_error = df_err_d,
#'   subtitle = "Diletti, Figure 1c"
#' )
#'
#' @keywords hplot design
#'
#' @family equivalence testing
#' @family plotting
#'
#' @export
#' @import stats
power_equivalence_md_plot <- function(alpha_level, logscale, lower_bound,
                                      upper_bound, true_sd, n_per_treatment,
                                      df_error, title = NULL,
                                      subtitle = NULL) {
  if (!missing(alpha_level)) .dmar_check_alpha_level(alpha_level)
  .dmar_check_numbers("power_equivalence_md_plot", c("lower_bound", "upper_bound", "true_sd"), required = c("lower_bound", "upper_bound", "true_sd"))
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    .dmar_refuse(c(
      "Package 'ggplot2' is required for power_equivalence_md_plot(). ",
      "Install with install.packages(\"ggplot2\")."),
      "package_required", package = "ggplot2")
  }
  if (length(n_per_treatment) != length(df_error)) {
    .dmar_refuse_argument(
      "'n_per_treatment' and 'df_error' must have the same length.",
      arg = c("n_per_treatment", "df_error"))
  }
  # Inside the function the sizes are n, the symbol of the figure's legend
  # and of the labels of the attached grid ("n=24").
  n        <- n_per_treatment

  nn       <- length(n)
  n_grid   <- 201L
  delta    <- (upper_bound - lower_bound) / (n_grid - 1L)
  diffs    <- lower_bound + delta * (0:(n_grid - 1L))

  # Compute power[n_i, diff_j], one row per sample size.
  power_array <- matrix(NA_real_, nrow = nn, ncol = n_grid,
                        dimnames = list(paste0("n=", n), NULL))
  for (i in seq_len(nn)) {
    for (j in seq_len(n_grid)) {
      power_array[i, j] <- power_equivalence_md(
        alpha_level = alpha_level, logscale = logscale,
        lower_bound = lower_bound, upper_bound = upper_bound,
        true_effect = diffs[j], true_sd = true_sd,
        n_per_treatment = n[i], df_error = df_error[i]
      )$value
    }
  }

  long <- data.frame(
    true_effect = rep(diffs, each = nn),
    n      = factor(rep(n, times = n_grid), levels = sort(unique(n))),
    power  = as.numeric(power_array),
    stringsAsFactors = FALSE
  )

  x_label <- if (logscale) "Ratio of test to reference" else
                            "Difference of test and reference"

  p <- ggplot2::ggplot(long, ggplot2::aes(x = .data$true_effect,
                                          y = .data$power,
                                          group = .data$n, color = .data$n)) +
    ggplot2::geom_hline(
      yintercept = c(.05, .1, .2, .3, .4, .5, .6, .7, .8, .9),
      color = "grey85", linewidth = 0.3
    ) +
    ggplot2::geom_line(linewidth = 0.7) +
    ggplot2::scale_y_continuous(limits = c(0, 1),
                                breaks = seq(0, 1, by = 0.1)) +
    ggplot2::scale_x_continuous(limits = c(lower_bound, upper_bound)) +
    ggplot2::labs(
      title    = if (is.null(title)) "Power of TOST" else title,
      subtitle = subtitle,
      x        = x_label,
      y        = "Power",
      color    = "n",
      caption  = sprintf("Equivalence: (%g, %g)   true_sd = %g",
                         lower_bound, upper_bound, true_sd)
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank()
    )

  # Attach numeric grid as documented (column 1 = the true effects; rest =
  # power per n).
  out_grid <- cbind(true_effect = diffs, t(power_array))
  colnames(out_grid) <- c("true_effect", paste0("n=", n))
  attr(p, "power_grid") <- out_grid

  p
}
