# --- Internal integrand for the TOST power calculation -------------------
# Numerical integration via stats::integrate() requires a function that takes
# a numeric vector and returns a numeric vector of the same length. The
# tidy/data-frame-returning public helper (power_density_equivalence_md())
# would not satisfy that contract, so we expose this private numeric version
# for the integrator and re-use it inside the public density function.
.power_equivalence_integrand <- function(power_sigma, alpha_level, theta1, theta2,
                                         diff, sigma, n, nu) {
  power_t      <- stats::qt(1 - alpha_level, df = nu)
  a            <- sigma / sqrt(nu)        # Chi-distribution scale constant
  # lgamma, not log(gamma(.)): gamma(nu / 2) overflows to Inf for nu beyond about
  # 344 (ordinary residual df), which would drive power_const to -Inf and the
  # integrand, and hence the power, silently to 0.
  power_const  <- -(nu / 2 - 1) * log(2) - lgamma(nu / 2) - log(a)
  half_width   <- power_t * sqrt(2 / n)
  one_over_se  <- 1 / (sigma * sqrt(2 / n))

  d1   <- ( power_sigma * half_width + (theta1 - diff)) * one_over_se
  d2   <- (-power_sigma * half_width + (theta2 - diff)) * one_over_se
  phi  <- stats::pnorm(d2) - stats::pnorm(d1)

  phi * exp(power_const + (nu - 1) * log(power_sigma / a) -
              0.5 * (power_sigma / a)^2)
}


# --- Shared scale-free quadrature for the TOST power ----------------------
# Consumed by power_equivalence_md() and, for the contrast generalization, by
# .power_equivalence_c_fast() in power_equivalence_c.R. Both reduce to the
# same one-dimensional problem once the equivalence limits are expressed in
# standard error units.
#
# Write r = S / sigma for the realized error standard deviation on a unit-free
# scale, t for the one-sided critical value, and z_lower, z_upper for the
# equivalence limits measured in standard errors from the true difference.
# Conditional on r, the (1 - 2 * alpha_level) confidence interval falls inside the
# equivalence interval with probability
#
#   phi(r) = Phi(z_upper - r * t) - Phi(z_lower + r * t),
#
# which reaches zero at r_max = (z_upper - z_lower) / (2 * t), the point where
# the interval is exactly as wide as the equivalence interval, and is zero
# beyond it. The power is the expectation of phi(r) over r = sqrt(X / nu) with
# X a chi square variate on nu degrees of freedom, so the density of r is
# 2 * nu * r * dchisq(nu * r^2, nu).
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
# deviations of r once nu is large), so the nodes always sample the mass, and
# the discarded tails carry at most 2e-18 of probability.
.power_equivalence_tost_quad <- function(power_t, z_lower, z_upper, nu) {
  r_max <- (z_upper - z_lower) / (2 * power_t)
  if (!is.finite(r_max) || r_max <= 0) return(0)

  tail_p <- 1e-18
  lower  <- sqrt(stats::qchisq(tail_p, df = nu) / nu)
  upper  <- sqrt(stats::qchisq(tail_p, df = nu, lower.tail = FALSE) / nu)
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
      2 * nu * r * stats::dchisq(nu * r^2, df = nu)
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
#' equivalence interval \eqn{[\theta_1, \theta_2]}, by numerical integration
#' over the chi distribution of the sample standard deviation.
#'
#' @param alpha_level Type I error rate for each of the two one-sided tests
#'   (typically \code{0.05}). The full equivalence test uses a
#'   \eqn{(1 - 2\alpha)} confidence interval.
#' @param logscale Logical. If \code{TRUE}, treatment means are compared on
#'   the logarithmic scale; \code{ltheta1}, \code{ltheta2}, and \code{ldiff}
#'   are expected as ratios (untransformed) and are log-transformed
#'   internally.
#' @param ltheta1 Lower limit of the equivalence interval (on the original scale; logged internally if \code{logscale = TRUE}).
#' @param ltheta2 Upper limit of the equivalence interval (on the original scale; logged internally if \code{logscale = TRUE}).
#' @param ldiff True difference in treatment means (or ratio on the log
#'   scale).
#' @param sigma \eqn{\sqrt{\mathrm{error\ variance}}}; root-MSE from an
#'   ANOVA. On the log scale, this is the coefficient of variation.
#' @param n Number of subjects per treatment (or total subjects in a
#'   crossover design).
#' @param nu Degrees of freedom associated with \code{sigma}.
#'
#' @return A one-row \code{data.frame} with columns \code{term} (\code{"power"})
#'   and \code{value} (the computed power, in \eqn{[0, 1]}).
#'
#' @details
#' The computation conditions on the error standard deviation the study will
#' actually observe. Given that value, whether the confidence interval fits
#' inside the equivalence interval is an ordinary normal probability, and the
#' power is that probability averaged over the chi distribution the error
#' standard deviation follows on \code{nu} degrees of freedom. The averaging
#' is carried out on a unit-free scale, with the equivalence limits expressed
#' in standard errors and the error standard deviation as a multiple of
#' \code{sigma}, so the power depends on the design rather than on the units
#' of the response: multiplying \code{ltheta1}, \code{ltheta2}, \code{ldiff},
#' and \code{sigma} by a common factor leaves the answer unchanged.
#'
#' For Phillips's (1990) original example (regular-scale two-period crossover
#' with \eqn{\theta_1 = -0.2}, \eqn{\theta_2 = 0.2}, CV \eqn{= 0.20},
#' \eqn{\delta = 0.05}, \eqn{n = 24}, \eqn{\nu = 22}), this function reproduces
#' the published value of \eqn{0.8029678} (Phillips, 1990, Table 1, 5th row,
#' 5th column).
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
#'   18}(2), 139--144. \doi{10.1007/BF01063556}
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
#' # Phillips (1990) Table 1, 5th row, 5th column. Expected: 0.8029678.
#' power_equivalence_md(alpha_level = .05, logscale = FALSE,
#'                      ltheta1 = -.2, ltheta2 = .2, ldiff = .05,
#'                      sigma = .20, n = 24, nu = 22)
#'
#' # Diletti (1991) Table 1, on the log scale (ratio of test to reference).
#' # Expected: 0.7922796.
#' power_equivalence_md(alpha_level = .05, logscale = TRUE,
#'                      ltheta1 = .8, ltheta2 = 1.25, ldiff = 1.05,
#'                      sigma = .20, n = 18, nu = 16)
#'
#' @keywords design
#'
#' @family equivalence testing
#'
#' @export
#' @import stats
power_equivalence_md <- function(alpha_level, logscale, ltheta1, ltheta2, ldiff,
                                 sigma, n, nu) {

  if (!is.numeric(alpha_level) || length(alpha_level) != 1L || alpha_level <= 0 || alpha_level >= 0.5) {
    stop("'alpha_level' must be a single number in (0, 0.5).")
  }
  if (!is.logical(logscale) || length(logscale) != 1L) {
    stop("'logscale' must be a single logical value.")
  }
  if (sigma <= 0) stop("'sigma' must be positive.")
  if (n <= 1)    stop("'n' must be greater than 1.")
  if (nu <= 0)   stop("'nu' must be positive.")
  if (ltheta2 <= ltheta1) {
    stop("'ltheta2' must be greater than 'ltheta1'.")
  }

  if (logscale) {
    if (ltheta1 <= 0 || ltheta2 <= 0 || ldiff <= 0) {
      stop("On the log scale, 'ltheta1', 'ltheta2', and 'ldiff' must all be positive.")
    }
    theta1 <- log(ltheta1)
    theta2 <- log(ltheta2)
    diff   <- log(ldiff)
  } else {
    theta1 <- ltheta1
    theta2 <- ltheta2
    diff   <- ldiff
  }

  # The power depends on the response scale only through the two equivalence
  # limits measured in standard errors from the true difference, so the
  # integration is carried out on that unit-free scale. See
  # .power_equivalence_tost_quad() above for the integrand and for why the
  # limits of integration are quantiles of the sampling distribution of the
  # error standard deviation rather than the raw interval it must fall in.
  power_t     <- stats::qt(1 - alpha_level, df = nu)
  one_over_se <- 1 / (sigma * sqrt(2 / n))

  power <- .power_equivalence_tost_quad(
    power_t = power_t,
    z_lower = (theta1 - diff) * one_over_se,
    z_upper = (theta2 - diff) * one_over_se,
    nu      = nu
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
#' @param power_sigma Numeric vector of \eqn{\sigma} values at which to
#'   evaluate the integrand.
#' @param alpha_level Type I error rate for each of the two one-sided tests.
#' @param theta1 Lower limit of the equivalence interval on the appropriate scale (regular or log).
#' @param theta2 Upper limit of the equivalence interval on the appropriate scale (regular or log).
#' @param diff True difference in treatment means (ratio on the log scale)
#'   on the appropriate scale.
#' @param sigma \eqn{\sqrt{\mathrm{error\ variance}}}.
#' @param n Number of subjects per treatment.
#' @param nu Degrees of freedom for \code{sigma}.
#'
#' @return A \code{data.frame} with one row per supplied
#'   \code{power_sigma}, and columns \code{power_sigma} and
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
#'   18}(2), 139--144. \doi{10.1007/BF01063556}
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
#' # Density at a single value of sigma:
#' power_density_equivalence_md(power_sigma = 0.10, alpha_level = .05,
#'                              theta1 = -.2, theta2 = .2, diff = .05,
#'                              sigma = .20, n = 24, nu = 22)
#'
#' # Vectorized over a grid:
#' grid <- power_density_equivalence_md(
#'   power_sigma = seq(0.01, 0.40, length.out = 50),
#'   alpha_level = .05, theta1 = -.2, theta2 = .2, diff = .05,
#'   sigma = .20, n = 24, nu = 22
#' )
#' head(grid)
#'
#' @keywords design
#'
#' @family equivalence testing
#'
#' @export
#' @import stats
power_density_equivalence_md <- function(power_sigma, alpha_level, theta1, theta2,
                                         diff, sigma, n, nu) {
  density <- .power_equivalence_integrand(
    power_sigma = power_sigma,
    alpha_level = alpha_level, theta1 = theta1, theta2 = theta2,
    diff = diff, sigma = sigma, n = n, nu = nu
  )
  data.frame(power_sigma = power_sigma, power_density = density,
             stringsAsFactors = FALSE)
}


#' Plot TOST Equivalence-Test Power Curves Over a Range of True Differences
#'
#' For each sample size in \code{n}, draws power as a function of the true
#' mean difference (or ratio, on the log scale), evaluated at 201 equally
#' spaced points across the equivalence interval. Returns a \pkg{ggplot2}
#' object; the underlying numerical grid is attached as
#' \code{attr(<plot>, "power_grid")}.
#'
#' @param alpha_level Type I error rate for each of the two one-sided tests.
#' @param logscale Logical. If \code{TRUE}, the means are compared on the
#'   logarithmic scale.
#' @param theta1 Lower limit of the equivalence interval.
#' @param theta2 Upper limit of the equivalence interval.
#' @param sigma \eqn{\sqrt{\mathrm{error\ variance}}}.
#' @param n Vector of sample sizes (one curve per element).
#' @param nu Vector of degrees of freedom for \code{sigma}, the same length
#'   as \code{n}.
#' @param title Optional plot title (default \code{"Power of TOST"}).
#' @param subtitle Optional subtitle (typically a reference like
#'   \code{"Phillips, Figure 3"}).
#'
#' @return A \code{ggplot} object. The 201-row power grid (column 1: true
#'   difference; remaining columns: power for each \code{n}) is attached as
#'   \code{attr(<plot>, "power_grid")}.
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
#'   18}(2), 139--144. \doi{10.1007/BF01063556}
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
#' # difference. Two of the seven sample sizes behind Phillips (1990)
#' # Figure 3 are drawn here so the example stays quick; the full
#' # reproduction is given below.
#' fig <- power_equivalence_md_plot(
#'   alpha_level = .05, logscale = FALSE,
#'   theta1 = -.2, theta2 = .2, sigma = .20,
#'   n = c(24, 60), nu = c(22, 58)
#' )
#' fig
#'
#' # The numbers behind the curves travel with the figure, so a particular
#' # power value can be read off rather than eyeballed. The first column is
#' # the true difference and the remaining columns give power, one column
#' # per sample size. Power is highest where the true difference is zero.
#' power_grid <- attr(fig, "power_grid")
#' power_grid[which.min(abs(power_grid[, 1])), ]
#'
#' # The two published figures are not run here because every curve
#' # evaluates the power integral at 201 true differences, so a
#' # seven-curve figure costs a few tenths of a second. Phillips (1990)
#' # Figure 3 is:
#' # n  <- c(9, 12, 18, 24, 30, 40, 60)
#' # nu <- c(7, 10, 16, 22, 28, 38, 58)
#' # power_equivalence_md_plot(
#' #   alpha_level = .05, logscale = FALSE,
#' #   theta1 = -.2, theta2 = .2, sigma = .20,
#' #   n = n, nu = nu,
#' #   subtitle = "Phillips Figure 3"
#' # )
#'
#' # Diletti (1991) Figure 1c is the same idea on the log scale, where the
#' # equivalence limits are the 0.80 to 1.25 ratio bounds used in
#' # bioequivalence work:
#' # n_d  <- c(8, 12, 18, 24, 30, 40, 60)
#' # nu_d <- c(6, 10, 16, 22, 28, 38, 58)
#' # power_equivalence_md_plot(
#' #   alpha_level = .05, logscale = TRUE,
#' #   theta1 = .8, theta2 = 1.25, sigma = .20,
#' #   n = n_d, nu = nu_d,
#' #   subtitle = "Diletti, Figure 1c"
#' # )
#'
#' @keywords hplot design
#'
#' @family equivalence testing
#' @family plotting
#'
#' @export
#' @import stats
power_equivalence_md_plot <- function(alpha_level, logscale, theta1, theta2, sigma,
                                      n, nu, title = NULL, subtitle = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for power_equivalence_md_plot(). ",
         "Install with install.packages(\"ggplot2\").", call. = FALSE)
  }
  if (length(n) != length(nu)) {
    stop("'n' and 'nu' must have the same length.", call. = FALSE)
  }

  nn       <- length(n)
  n_grid   <- 201L
  delta    <- (theta2 - theta1) / (n_grid - 1L)
  diffs    <- theta1 + delta * (0:(n_grid - 1L))

  # Compute power[n_i, diff_j] -- one row per sample size.
  power_array <- matrix(NA_real_, nrow = nn, ncol = n_grid,
                        dimnames = list(paste0("n=", n), NULL))
  for (i in seq_len(nn)) {
    for (j in seq_len(n_grid)) {
      power_array[i, j] <- power_equivalence_md(
        alpha_level = alpha_level, logscale = logscale,
        ltheta1 = theta1, ltheta2 = theta2, ldiff = diffs[j],
        sigma = sigma, n = n[i], nu = nu[i]
      )$value
    }
  }

  long <- data.frame(
    diff   = rep(diffs, each = nn),
    n      = factor(rep(n, times = n_grid), levels = sort(unique(n))),
    power  = as.numeric(power_array),
    stringsAsFactors = FALSE
  )

  x_label <- if (logscale) "Ratio of test to reference" else
                            "Difference of test and reference"

  p <- ggplot2::ggplot(long, ggplot2::aes(x = .data$diff, y = .data$power,
                                          group = .data$n, color = .data$n)) +
    ggplot2::geom_hline(
      yintercept = c(.05, .1, .2, .3, .4, .5, .6, .7, .8, .9),
      color = "grey85", linewidth = 0.3
    ) +
    ggplot2::geom_line(linewidth = 0.7) +
    ggplot2::scale_y_continuous(limits = c(0, 1),
                                breaks = seq(0, 1, by = 0.1)) +
    ggplot2::scale_x_continuous(limits = c(theta1, theta2)) +
    ggplot2::labs(
      title    = if (is.null(title)) "Power of TOST" else title,
      subtitle = subtitle,
      x        = x_label,
      y        = "Power",
      color    = "n",
      caption  = sprintf("Equivalence: (%g, %g)   sigma = %g",
                         theta1, theta2, sigma)
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank()
    )

  # Attach numeric grid as documented (column 1 = diffs; rest = power per n).
  out_grid <- cbind(diff = diffs, t(power_array))
  colnames(out_grid) <- c("diff", paste0("n=", n))
  attr(p, "power_grid") <- out_grid

  p
}
