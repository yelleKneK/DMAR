# Limits of agreement (Bland-Altman) with CIs on the limits.
#' Limits of Agreement (Bland-Altman) With Confidence Intervals on the Limits
#'
#' Computes the limits of agreement (LoA) of Bland and Altman (1986, 1999)
#' between two methods of measurement applied to the same units, along
#' with the Carkeet (2015) exact confidence intervals on the LoA
#' themselves. The CIs treat the two limits either as a pair (the
#' default), so that the confidence statement holds for both limits
#' jointly, or individually, one limit at a time. Both constructions
#' replace the approximate normal CIs originally given by Bland and
#' Altman (1999), which are too narrow at small \emph{n}. The short
#' alias \code{loa()} calls the same function.
#'
#' @param x,y Paired numeric vectors of equal length (e.g., method A
#'   and method B applied to the same units).
#' @param coverage Probability content of the limits of agreement.
#'   Default \code{0.95} (the conventional 95\% LoA). The mean difference
#'   is bracketed by \eqn{\pm z_{(1 + \text{coverage})/2}} standard
#'   deviations of the differences.
#' @param conf_level Confidence level for the CIs on the LoA themselves.
#'   Default \code{0.95}.
#' @param method How the CIs on the LoA are constructed, following
#'   Carkeet (2015): \code{"pair"} (the default) treats the two limits
#'   as a pair, so that the confidence statement holds for both limits
#'   jointly; \code{"individual"} treats each limit separately. See
#'   Details for when each is appropriate.
#'
#' @return A \code{data.frame} with rows for the mean difference,
#'   the SD of differences, the lower and upper LoA (\code{loa_lower},
#'   \code{loa_upper}), and the lower / upper CI bounds on each LoA.
#'   The rows are the same under both methods; the construction that
#'   produced the CI bounds is recorded in the \code{method} attribute.
#'
#' @details
#' \strong{Definition.} For paired observations \eqn{(x_i, y_i)}, the
#' Bland-Altman limits of agreement are
#' \deqn{\mathrm{LoA}_\pm \;=\; \bar d \pm k \cdot s_d,}
#' where \eqn{d_i = y_i - x_i}, \eqn{\bar d} is the mean of the
#' differences, \eqn{s_d} is their SD, and
#' \eqn{k = z_{(1 + \mathrm{coverage})/2}} (for 95\% coverage,
#' \eqn{k = 1.96}). The LoA are population intervals: they describe the
#' range within which approximately \emph{coverage}\% of \emph{individual}
#' differences are expected to lie if the differences are normally
#' distributed.
#'
#' \strong{CIs on the LoA themselves.} The sample LoA are random
#' variables, and Carkeet (2015) derived exact CIs for them in two
#' forms, selected by \code{method}. The choice turns on what the
#' agreement claim is about.
#'
#' \strong{The pair method (the default).} A Bland-Altman analysis is
#' usually read as a statement about the range of agreement as a whole,
#' the span from the lower to the upper LoA within which about
#' \emph{coverage}\% of individual differences lie. That claim involves
#' both limits at once, so the confidence statement should hold for the
#' two limits jointly; this is the treatment Carkeet (2015) recommends
#' for most situations. Writing \eqn{k_t(F)} for the exact two-sided
#' normal tolerance factor with confidence \eqn{F} and content equal to
#' \code{coverage} (Odeh, 1978), the CI on the upper LoA is
#' \deqn{\left[\, \bar d + k_t(\alpha/2)\, s_d, \;\;
#'   \bar d + k_t(1 - \alpha/2)\, s_d \,\right],}
#' with \eqn{\alpha} equal to one minus \code{conf_level}, and the CI
#' on the lower LoA is its mirror image about \eqn{\bar d}. The joint confidence
#' statement runs through the probability content of the two symmetric
#' intervals: with confidence \code{conf_level}, the interval between
#' the inner pair of bounds, \eqn{\bar d \pm k_t(\alpha/2)\, s_d},
#' captures less than \emph{coverage}\% of the population of
#' differences, while the interval between the outer pair,
#' \eqn{\bar d \pm k_t(1 - \alpha/2)\, s_d}, captures more, so the pair
#' of population limits is bracketed simultaneously. In the Bland and
#' Altman (1986) example that Carkeet reanalyzes (\eqn{n = 17},
#' \eqn{\bar d = -2.1}, \eqn{s_d = 38.8}), the pair bounds are
#' \eqn{-2.1 \pm 57.81} (inner) and \eqn{-2.1 \pm 119.60} (outer).
#'
#' \strong{The individual method.} When a single limit carries the
#' substantive question (for example, only the upper limit matters
#' because only differences in one direction are clinically
#' consequential), each limit can be treated on its own. Writing
#' \eqn{t_{p,\, n - 1}(\delta)} for the \eqn{p} quantile of the
#' noncentral \emph{t} distribution with \eqn{n - 1} degrees of freedom
#' and noncentrality parameter \eqn{\delta = k \sqrt{n}}, the exact CI
#' on the upper LoA is
#' \deqn{\left[\, \bar d + \frac{s_d}{\sqrt{n}}\,
#'     t_{\alpha/2,\, n - 1}(\delta), \;\;
#'   \bar d + \frac{s_d}{\sqrt{n}}\,
#'     t_{1 - \alpha/2,\, n - 1}(\delta) \,\right],}
#' and the CI on the lower LoA uses \eqn{-\delta} in place of
#' \eqn{\delta}. These intervals are asymmetric about the sample LoA,
#' wider on the side away from the mean difference. In the worked
#' example above, the individual CI on the upper LoA is
#' \eqn{[48.9,\, 120.0]}. The confidence statement is per limit: each
#' limit is covered with \code{conf_level} confidence separately, not
#' both at once.
#'
#' \strong{Numerical accuracy.} The pair tolerance factors are computed
#' by numerical integration of the Odeh (1978) chi square by normal
#' integral, which reproduces Carkeet's Table 2 to all four printed
#' decimals. The individual quantiles come from \code{stats::qt} with a
#' noncentrality parameter, so at very large \emph{n} their accuracy is
#' bounded by R's noncentral \emph{t} algorithm: near \eqn{n = 1000}
#' the tolerance coefficient carries an error of about
#' \eqn{3 \times 10^{-4}}, far past the sample sizes at which the
#' exact-versus-approximate distinction matters.
#'
#' \strong{Caveats.} The LoA construction assumes (i) the differences
#' \eqn{d_i} are approximately normally distributed, and (ii) the
#' difference does not systematically depend on the magnitude of the
#' measurement (proportional bias). Both should be checked, the second
#' by plotting \eqn{d_i} against \eqn{(x_i + y_i)/2}; a non-flat
#' relationship indicates that a single set of LoA is inappropriate.
#'
#' @references
#' Bland, J. M., & Altman, D. G. (1986). Statistical methods for
#'   assessing agreement between two methods of clinical measurement.
#'   \emph{Lancet, 327}(8476), 307--310.
#'
#' Bland, J. M., & Altman, D. G. (1999). Measuring agreement in method
#'   comparison studies. \emph{Statistical Methods in Medical Research,
#'   8}(2), 135--160. \doi{10.1191/096228099673819272}
#'
#' Carkeet, A. (2015). Exact parametric confidence intervals for
#'   Bland-Altman limits of agreement. \emph{Optometry and Vision
#'   Science, 92}(3), e71--e80. \doi{10.1097/OPX.0000000000000513}
#'
#' Odeh, R. E. (1978). Tables of two-sided tolerance factors for a
#'   normal distribution. \emph{Communications in Statistics -
#'   Simulation and Computation, 7}(2), 183--201.
#'
#' @seealso \code{\link{lin_ccc}}
#'
#' @examples
#' # 1. Two methods that agree well; the CIs treat the limits as a
#' #    pair (the default):
#' set.seed(113)
#' method_a <- rnorm(40, mean = 100, sd = 15)
#' method_b <- method_a + rnorm(40, mean = 0, sd = 3)
#' limits_of_agreement(method_a, method_b)
#'
#' # 2. Each limit treated individually, for when a single limit
#' #    carries the substantive question:
#' limits_of_agreement(method_a, method_b, method = "individual")
#'
#' # 3. 90% LoA with 95% CIs on the limits:
#' limits_of_agreement(method_a, method_b, coverage = 0.90, conf_level = 0.95)
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @keywords htest
#'
#' @family agreement and measurement
#'
#' @export

limits_of_agreement <- function(x, y, coverage = 0.95, conf_level = 0.95,
                method = c("pair", "individual")) {
  method <- match.arg(method)
  if (!is.numeric(x) || !is.numeric(y))
    stop("'x' and 'y' must both be numeric vectors.")
  if (length(x) != length(y))
    stop("'x' and 'y' must be the same length (paired observations).")
  n <- length(x)
  if (n < 4L)
    stop("Need at least 4 pairs.")
  if (coverage <= 0 || coverage >= 1)
    stop("'coverage' must be in (0, 1).")
  if (conf_level <= 0 || conf_level >= 1)
    stop("'conf_level' must be in (0, 1).")

  d   <- y - x
  d_m <- mean(d)
  d_s <- sd(d)
  k   <- stats::qnorm((1 + coverage) / 2)

  loa_lower <- d_m - k * d_s
  loa_upper <- d_m + k * d_s

  nu    <- n - 1L
  alpha <- 1 - conf_level

  if (method == "pair") {
    # Carkeet's (2015) pair method: exact two-sided tolerance factors
    # (Odeh, 1978), so the confidence statement covers both limits
    # jointly. The four CI bounds are symmetric about the mean
    # difference: d_m +/- kt_inner * s_d and d_m +/- kt_outer * s_d.
    kt_inner <- .loa_kt_pair(alpha / 2, nu, coverage)
    kt_outer <- .loa_kt_pair(1 - alpha / 2, nu, coverage)
    loa_upper_lo <- d_m + kt_inner * d_s
    loa_upper_hi <- d_m + kt_outer * d_s
    loa_lower_lo <- d_m - kt_outer * d_s
    loa_lower_hi <- d_m - kt_inner * d_s
  } else {
    # Carkeet's (2015) individual method: each limit on its own through
    # the noncentral t distribution with noncentrality parameter
    # ncp = k * sqrt(n). (Bland and Altman's 1999 approximate CI,
    # LoA +/- t * s_d * sqrt(1/n + k^2 / (2 * (n - 1))), is what this
    # exact construction replaces.)
    ncp <- k * sqrt(n)

    # Upper LoA: CI = d_m + (s_d / sqrt(n)) * (quantile of noncentral t).
    t_lo <- stats::qt(alpha / 2,     df = nu, ncp = ncp)
    t_hi <- stats::qt(1 - alpha / 2, df = nu, ncp = ncp)
    loa_upper_lo <- d_m + (d_s / sqrt(n)) * t_lo
    loa_upper_hi <- d_m + (d_s / sqrt(n)) * t_hi

    # Lower LoA: the same construction with -ncp. R's pnt/qt emits a
    # benign "full precision may not have been achieved in 'pnt{final}'"
    # warning for the negative noncentrality parameter; the returned
    # quantiles are accurate to the reported precision, so muffle just
    # that one message rather than surface it on every call (including
    # the examples).
    ncp_neg <- -k * sqrt(n)
    muffle_pnt <- function(expr) {
      withCallingHandlers(
        expr,
        warning = function(w) {
          if (grepl("pnt\\{final\\}", conditionMessage(w)))
            invokeRestart("muffleWarning")
        }
      )
    }
    t_lo_n <- muffle_pnt(stats::qt(alpha / 2,     df = nu, ncp = ncp_neg))
    t_hi_n <- muffle_pnt(stats::qt(1 - alpha / 2, df = nu, ncp = ncp_neg))
    loa_lower_lo <- d_m + (d_s / sqrt(n)) * t_lo_n
    loa_lower_hi <- d_m + (d_s / sqrt(n)) * t_hi_n
  }

  out <- data.frame(
    term  = c("mean_difference", "sd_difference",
              "loa_lower", "loa_lower_lower_limit", "loa_lower_upper_limit",
              "loa_upper", "loa_upper_lower_limit", "loa_upper_upper_limit"),
    value = c(d_m, d_s,
              loa_lower, loa_lower_lo, loa_lower_hi,
              loa_upper, loa_upper_lo, loa_upper_hi),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  attr(out, "method") <- method
  .as_dmar_tbl(out, conf_level = conf_level)
}

#' @rdname limits_of_agreement
#' @export
loa <- limits_of_agreement

# Exact two-sided normal tolerance factor k_t with confidence f and
# probability content `coverage`, at nu = n - 1 degrees of freedom: the
# value satisfying
#   Pr{ d_bar +/- k_t * s_d contains at least `coverage` of the
#       population of differences } = f.
# Computed from the Odeh (1978) chi square by normal integral, the same
# formulation Carkeet (2015) used to build his Table 2 (this
# implementation reproduces that table to all four printed decimals):
#   (2 sqrt(n) / sqrt(2 pi)) *
#     Int_0^inf Pr{ chisq_nu > nu r(u)^2 / k_t^2 } exp(-n u^2 / 2) du = f,
# where r(u) solves Phi(u + r) - Phi(u - r) = coverage. The inner
# uniroot finds r(u); the outer uniroot solves for k_t. Not exported;
# consumed only by limits_of_agreement(method = "pair").
.loa_kt_pair <- function(f, nu, coverage) {
  n <- nu + 1
  r_of <- function(u) {
    stats::uniroot(function(r)
      stats::pnorm(u + r) - stats::pnorm(u - r) - coverage,
      c(1e-8, u + 50), tol = 1e-13)$root
  }
  u_max <- sqrt(220 / n)  # exp(-n u^2 / 2) < 1e-47 beyond this
  conf_of <- function(kt) {
    integrand <- function(u) {
      vapply(u, function(ui)
        stats::pchisq(nu * r_of(ui)^2 / kt^2, nu, lower.tail = FALSE) *
          exp(-n * ui^2 / 2), numeric(1L))
    }
    2 * sqrt(n) / sqrt(2 * pi) *
      stats::integrate(integrand, 0, u_max, rel.tol = 1e-10)$value
  }
  stats::uniroot(function(kt) conf_of(kt) - f, c(1e-6, 1000),
                 tol = 1e-10)$root
}
