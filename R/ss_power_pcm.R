#' Sample Size Planning for Power for Polynomial Change Models
#'
#' Returns power given the sample size, or sample size given the desired power,
#' for the group difference in a polynomial change coefficient (a flat-line
#' intercept, a linear slope, a quadratic acceleration, or any higher-order
#' trend) in a two-group longitudinal design, following Raudenbush and Liu
#' (2001). The trend whose group difference is tested is selected with
#' \code{trend}; \code{trend = "linear"} (the default) reproduces the
#' straight-line case.
#'
#' @param beta The level two regression coefficient for the group by time
#'   interaction in the polynomial change coefficient selected by \code{trend}
#'   (linear by default), where the grouping variable is coded -.5 and .5 for
#'   the two groups. When \code{standardized = TRUE} (the default) this is the
#'   standardized change difference (see \code{standardized}).
#' @param tau The true between-subject variance of the individuals' change
#'   coefficient for the trend selected by \code{trend} (the variance of the
#'   slopes for \code{trend = "linear"}).
#' @param level_1_variance Level one (within-subject) error variance
#' @param frequency Number of measurements per unit of time, where the unit is
#'   the one in which \code{duration} is expressed. It need not be a whole
#'   number; for example, \code{frequency = 0.5} means one measurement every
#'   two time units. Together with \code{duration} it fixes the number of
#'   equally spaced measurement occasions, \eqn{M = f \times D + 1}, with
#'   \eqn{f} the frequency and \eqn{D} the duration.
#' @param duration Length of the study in the chosen time unit (for example,
#'   years, grades, or hours). Measurements are taken at
#'   \eqn{M = f \times D + 1} equally spaced occasions spanning time 0 to
#'   \code{duration}.
#' @param desired_power Desired power
#' @param N Total sample size (one-half in each of the two groups)
#' @param alpha_level Type I error rate
#' @param standardized The standardized change difference is the unstandardized
#'   change difference divided by the square root of \code{tau}, the
#'   between-subject variance of the change coefficient. \code{TRUE} (the
#'   default) treats \code{beta} as already standardized; \code{FALSE} treats it
#'   as the raw (unstandardized) change difference.
#' @param directional Should a one (\code{TRUE}) or two (\code{FALSE}) tailed test be performed.
#' @param trend The polynomial change coefficient whose group difference is
#'   tested, given either as a name (\code{"intercept"}, \code{"linear"},
#'   \code{"quadratic"}, \code{"cubic"}, \code{"quartic"}, ...) or as a
#'   non-negative integer order (\code{0} = intercept / flat line, \code{1} =
#'   linear, \code{2} = quadratic, ...). Defaults to \code{"linear"}. The design
#'   must supply at least \eqn{p + 1} measurement occasions to estimate a
#'   degree-\eqn{p} trend.
#'
#' @references
#' Kelley, K., & Maxwell, S. E. (2008). Sample size planning with
#'   applications to multiple regression: Power and accuracy for omnibus
#'   and targeted effects. In P. Alasuutari, L. Bickman, & J. Brannen
#'   (Eds.), \emph{The Sage handbook of social research methods}
#'   (pp. 166--192). Sage.
#'
#' Kelley, K., & Rausch, J. R. (2011). Sample size planning for
#'   longitudinal models: Accuracy in parameter estimation for polynomial
#'   change parameters. \emph{Psychological Methods, 16}(4), 391--405.
#'   \doi{10.1037/a0023352}
#'
#' Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). \emph{Designing
#'   experiments and analyzing data: A model comparison perspective}
#'   (4th ed.). Routledge. (See Chapters 11, 15.)
#'
#' Maxwell, S. E., Kelley, K., & Rausch, J. R. (2008). Sample size planning
#'   for statistical power and accuracy in parameter estimation.
#'   \emph{Annual Review of Psychology, 59}, 537--563.
#'   \doi{10.1146/annurev.psych.59.103006.093735}
#'
#' Raudenbush, S. W., & Liu, X.-F. (2001). Effects of study duration,
#'   frequency of observation, and sample size on power in studies of
#'   group differences in polynomial change.
#'   \emph{Psychological Methods, 6}(4), 387--401.
#'   \doi{10.1037/1082-989X.6.4.387}
#'
#' @return A \code{data.frame} (class \code{dmar_tbl}) with one row per
#'   reported quantity in a \code{term} / \code{value} layout: the per-group
#'   size per group (\code{necessary_n_per_group}, or \code{specified_n_per_group}
#'   when \code{N} is supplied) and the total (\code{total_N});
#'   the achieved power (\code{actual_power}); the measurement schedule
#'   (\code{freq}, \code{duration}, \code{measurement_occasions}); the
#'   polynomial order of the tested change coefficient (\code{polynomial_order},
#'   0 = intercept, 1 = linear, 2 = quadratic, ...); the
#'   unstandardized and standardized change difference
#'   (\code{unstd_coefficient}, \code{std_coefficient}); the level one error
#'   variance (\code{l1_error_var}); the true and error variance of the change
#'   coefficient (\code{true_var_of_slopes}, \code{error_var_of_slopes}, whose
#'   names retain "slopes" from the linear case); the change-coefficient
#'   reliability (\code{reliability}); and the noncentrality parameter of the
#'   \emph{t} test (\code{noncentral_t_parm}).
#'
#' @details
#' The two groups each contain \eqn{N / 2} subjects measured on
#' \eqn{M = f \times D + 1} equally spaced occasions. Each subject's degree-\eqn{p}
#' polynomial change coefficient is estimated within subject; the test compares
#' the two group means of that coefficient. The change coefficient is taken in
#' the derivative-scaled metric (\eqn{p!} times the leading coefficient of
#' \eqn{t^p}), the metric in which the Raudenbush and Liu (2001) constants
#' apply.
#'
#' The within-subject sampling variance of the estimated coefficient is
#' (Raudenbush & Liu, 2001, p. 392)
#' \deqn{V = \sigma^2_e\, f^{2p}\,\frac{(M - p - 1)!}{K_p\,(M + p)!}, \qquad
#'       \frac{1}{K_p} = \frac{(2p)!\,(2p+1)!}{(p!)^2},}
#' so that \eqn{1/K_p = 1, 12, 720, 100800, \ldots} for \eqn{p = 0, 1, 2, 3,
#' \ldots}. This \eqn{V} equals \eqn{(p!)^2} times the variance of the ordinary
#' least squares estimate of the coefficient of \eqn{t^p}, reduces to
#' \eqn{\sigma^2_e / M} at \eqn{p = 0} and to
#' \eqn{12\,\sigma^2_e f^2 / [M(M^2 - 1)]} at \eqn{p = 1}. The slope reliability
#' is \eqn{\tau / (\tau + V)} (their Equation 15), the variance of the
#' between-group difference is \eqn{4(\tau + V)/N} (their Equation 10), and the
#' \emph{t} test has \eqn{N - 2} degrees of freedom and noncentrality
#' \eqn{\sqrt{N\,\beta^2\,[\tau/(\tau + V)]/4}} (their Equations 12 and 14). The
#' linear case reproduces the National Youth Survey benchmark in their Tables 1
#' and 2; the general-\eqn{p} formula has been checked against the exact
#' \eqn{(X'X)^{-1}} variance and against an end-to-end Monte Carlo power study
#' for the quadratic trend.
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @examples
#' # The examples reproduce the National Youth Survey illustration of
#' # Raudenbush and Liu (2001, p. 393). One observation per year
#' # (frequency = 1) over a four-year study (duration = 4) gives
#' # M = frequency * duration + 1 = 5 equally spaced occasions. The
#' # standardized slope difference is -0.40, the slope variance is
#' # tau = 0.003, and the level-one error variance is 0.0262, so the slope
#' # reliability is 0.53.
#'
#' # (1) Power at a given sample size. With N = 238 (119 per group) the
#' #     design has power 0.61 (Raudenbush and Liu, 2001, Table 1, D = 4,
#' #     f = 1).
#' ss_power_pcm(beta = -.4, tau = .003, level_1_variance = .0262,
#'              frequency = 1, duration = 4, N = 238)
#'
#' # (2) Sample size for a target power. Solving the same design for 0.80
#' #     power returns N = 370 (185 per group), between their Table 2 cells
#' #     N = 300 (power 0.71) and N = 400 (power 0.83).
#' ss_power_pcm(beta = -.4, tau = .003, level_1_variance = .0262,
#'              frequency = 1, duration = 4, desired_power = .80)
#'
#' # (3) Unstandardized slope. The unstandardized slope difference is
#' #     beta * sqrt(tau) = -0.40 * sqrt(0.003) = -0.0219. Passing it with
#' #     standardized = FALSE reproduces the power of 0.61 from example (1).
#' ss_power_pcm(beta = -.0219, tau = .003, level_1_variance = .0262,
#'              frequency = 1, duration = 4, N = 238, standardized = FALSE)
#'
#' # (4) Longer study, same number of occasions. Doubling the duration to
#' #     D = 8 while keeping M = 5 (so frequency = 0.5, one observation every
#' #     two years) raises the slope reliability to 0.82 and power to about
#' #     0.80. Spreading the same five occasions over a longer span sharply
#' #     increases power (Raudenbush & Liu, 2001, p. 393).
#' ss_power_pcm(beta = -.4, tau = .003, level_1_variance = .0262,
#'              frequency = .5, duration = 8, N = 238)
#'
#' # (5) More frequent sampling over a shorter span, same occasions. Halving
#' #     the duration to D = 2 while keeping M = 5 (so frequency = 2) drops
#' #     the slope reliability to 0.22 and power to about 0.31. Sampling more
#' #     often over a shorter study does little for power (Raudenbush & Liu,
#' #     2001, p. 393).
#' ss_power_pcm(beta = -.4, tau = .003, level_1_variance = .0262,
#'              frequency = 2, duration = 2, N = 238)
#'
#' # (6) One-sided test. A directional test of the base design places the
#' #     whole Type I error rate in the predicted tail, raising power from
#' #     0.61 to about 0.73.
#' ss_power_pcm(beta = -.4, tau = .003, level_1_variance = .0262,
#'              frequency = 1, duration = 4, N = 238, directional = TRUE)
#'
#' # (7) A higher-order trend. The same machinery plans power for the group
#' #     difference in any polynomial change coefficient. Here the target is
#' #     the quadratic trend (curvature / acceleration): with eight occasions
#' #     (frequency = 1, duration = 7), a between-subject quadratic-coefficient
#' #     variance tau = 0.002, level-one error variance 0.05, and a
#' #     standardized quadratic difference of 0.45, the design is planned for
#' #     0.80 power. A quadratic trend needs at least three occasions; a cubic
#' #     at least four (trend = "cubic" or trend = 3), and so on.
#' ss_power_pcm(beta = 0.45, tau = 0.002, level_1_variance = 0.05,
#'              frequency = 1, duration = 7, desired_power = .80,
#'              trend = "quadratic")
#'
#' @keywords design
#'
#' @seealso \code{\link{design_consequences}} for what a chosen design delivers:
#'   power, the Type S (sign) and Type M (exaggeration) errors of the
#'   significance filter, and the expected confidence interval width.
#'
#' @family sample size for power
#'
#' @export

ss_power_pcm <- function(beta, tau, level_1_variance, frequency, duration, desired_power = NULL,
                         N = NULL, alpha_level = 0.05, standardized = TRUE, directional = FALSE,
                         trend = "linear") {

  # Exactly one of `N` and `desired_power` is the unknown to solve for: supply
  # `desired_power` to solve for the required sample size, or supply `N` to
  # solve for the power achieved at that sample size. Supplying neither leaves
  # nothing to solve for (the routine would otherwise fail with an obscure
  # length-zero error in the search loop); supplying both over-determines the
  # problem.
  if (is.null(N) && is.null(desired_power))
    stop("Specify one of 'N' or 'desired_power': supply 'desired_power' to ",
         "solve for the required sample size, or 'N' to solve for the ",
         "achieved power.", call. = FALSE)
  if (!is.null(N) && !is.null(desired_power))
    stop("Specify only one of 'N' or 'desired_power', not both: the one left ",
         "unspecified is the quantity this function solves for.", call. = FALSE)

  # Polynomial order of the change coefficient being tested: p = 0 the
  # intercept (a flat line), 1 the linear slope, 2 the quadratic acceleration,
  # and so on. The within-subject sampling variance V below is general in p
  # through the Raudenbush and Liu (2001) constant K_p, with
  # 1/K_p = (2p)! (2p+1)! / (p!)^2 (= 1, 12, 720, 100800, ... for p = 0, 1, 2,
  # 3, ...); .pcm_sampling_variance() carries that closed form.
  p <- .pcm_trend_order(trend)

  # Number of equally spaced measurement occasions (Raudenbush & Liu, 2001,
  # p. 392, lower left): M = f * D + 1.
  M <- frequency * duration + 1

  if (standardized == FALSE) {
    beta_standardized <- beta / sqrt(tau) # Equation 13
    beta_unstandardized <- beta
  }

  if (standardized == TRUE) {
    beta_standardized <- beta
    beta_unstandardized <- beta * sqrt(tau)
  }

  if (is.null(N)) {
    # Power of the two-sided slope test at a candidate total sample size N_i.
    power_at_N <- function(N_i) {
      V <- .pcm_sampling_variance(level_1_variance, frequency, M, p)
      reliability <- tau / (tau + V)
      lambda <- if (standardized) (N_i * beta^2 * reliability) / 4 else
        (beta^2) / (4 * (tau + V) / N_i)
      .power_noncentral_t(ncp = sqrt(lambda), df = N_i - 2,
                          alpha_level = alpha_level, directional = directional)
    }
    # Fail fast when the target cannot be reached at any realistic size. A null
    # slope (beta = 0) has power equal to alpha_level at every N, so its power
    # never reaches a higher target; a tiny nonzero slope reaches it only at an
    # unusable size. Searching without a bound would not terminate.
    N_max <- 1e7
    if (power_at_N(N_max) < desired_power)
      stop("Could not reach 'desired_power' at a total sample size up to ",
           format(N_max, scientific = FALSE), "; the slope may be too small ",
           "to detect at this power, or the target may be unreachable (a null ",
           "slope has power equal to 'alpha_level' at every sample size).",
           call. = FALSE)
    N_i <- 4
    repeat {
      N_i <- N_i + 2
      if (power_at_N(N_i) >= desired_power) break
    }
    # Recompute the design quantities the output reports at the resolved N_i.
    V <- .pcm_sampling_variance(level_1_variance, frequency, M, p) # R&L (2001) p. 392; general in p
    var_beta <- 4 * (tau + V) / N_i                               # Equation 10
    reliability <- tau / (tau + V)                                # Equation 15
    lambda <- if (standardized) (N_i * beta^2 * reliability) / 4 else # Eq 14 / 12
      (beta^2) / var_beta
    actual_power <- power_at_N(N_i)
  }


  if (is.null(desired_power)) {
    V <- .pcm_sampling_variance(level_1_variance, frequency, M, p) # R&L (2001) p. 392; general in p

    var_beta <- 4 * (tau + V) / N
    reliability <- tau / (tau + V) # Equation 15
    if (standardized == FALSE) {
      lambda <- (beta^2) / var_beta # Equation 12
    }
    if (standardized == TRUE) {
      lambda <- (N * beta^2 * reliability) / 4
    }
    actual_power <- .power_noncentral_t(ncp = sqrt(lambda), df = N - 2,
                                        alpha_level = alpha_level,
                                        directional = directional)
    N_i <- N
  }

  # One per-group size row named for its role replaces the MBESS-era ss_c
  # and ss_t pair: the design is balanced by construction, so the control
  # and treatment groups are the same size and two rows said one thing
  # twice. The implied total stays as total_N.
  size_term <- if (is.null(N)) "necessary_n_per_group" else "specified_n_per_group"
  return(
    .as_dmar_tbl(data.frame(
      term = c(size_term, 'total_N', 'actual_power', 'freq', 'duration', 'measurement_occasions',
               'polynomial_order', 'unstd_coefficient', 'std_coefficient', 'l1_error_var', 'true_var_of_slopes',
               'error_var_of_slopes', 'reliability', 'noncentral_t_parm'),
      value = c(N_i / 2, N_i, actual_power, frequency, duration, M,
                p, beta_unstandardized, beta_standardized, level_1_variance, tau,
                V, reliability, sqrt(lambda))),
      subclass = "dmar_ss_power")
  )
}
