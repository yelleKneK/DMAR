#' Confidence Interval for the Population Root Mean Square Error of Approximation
#'
#' Constructs a confidence interval for the population root mean square
#' error of approximation (RMSEA), a population badness-of-fit index for
#' structural equation models, from a reported RMSEA, its degrees of
#' freedom, and the sample size. The interval is obtained by inverting the
#' noncentral chi square distribution of the test statistic \eqn{T =
#' a\hat{F}}, where \eqn{\hat{F}} is the minimized fit function and \eqn{a}
#' the sample size multiplier, and mapping the resulting noncentrality
#' limits to the RMSEA metric (Steiger & Lind, 1980; Browne & Cudeck, 1993).
#'
#' @param rmsea Observed root mean square error of approximation
#' @param df Degrees of freedom of the model, greater than 0 (a model with
#'   zero degrees of freedom is saturated, and its RMSEA is undefined)
#' @param N Sample size: greater than 1 under the default
#'   \code{sample_size_convention}, and greater than 0 under \code{"N"}
#' @param conf_level Desired confidence level, between 0 and 1. The
#'   default, 0.90, gives the conventional interval for the RMSEA, whose
#'   two limits decide the tests of close and not close fit at .05 (see
#'   Details).
#' @param alpha_lower The Type I error rate for the lower tail
#' @param alpha_upper The Type I error rate for the upper tail.
#'   The table records the coverage, \code{1 - alpha_lower - alpha_upper},
#'   as its level and, when the two tail areas differ, the areas
#'   themselves, which the line beneath the printed table,
#'   \code{\link{results_sentence}}, and the \code{alpha_lower} and
#'   \code{alpha_upper} columns of \code{tidy()} state; with a tail area
#'   of 0 the interval is a one-sided confidence bound.
#' @param sample_size_convention The multiplier \eqn{a} with which the
#'   reported RMSEA was computed: \code{"N-1"} (the default), so \eqn{a = N
#'   - 1}, the convention of Browne and Cudeck (1993), of Kelley and Lai
#'   (2011), and of \code{\link{ss_aipe_rmsea}}; or \code{"N"}, so \eqn{a =
#'   N}, the convention of lavaan's default normal likelihood and of Mplus.
#'
#' @details
#' The RMSEA expresses the badness of model fit per degree of freedom on
#' the noncentrality scale. Under the noncentral chi square model for the
#' sample fit statistic, \eqn{T = a\hat{F}} has approximate noncentral chi
#' square distribution with \eqn{df} degrees of freedom and noncentrality
#' parameter \eqn{\lambda = a \cdot df \cdot \mathrm{RMSEA}^2}. The
#' function rebuilds \eqn{T = a \cdot df \cdot \mathrm{RMSEA}^2 + df} from
#' the reported RMSEA, obtains the interval for \eqn{\lambda} by inverting
#' the noncentral chi square distribution at the requested confidence level
#' (\code{\link{ci_nc_chisq}} does the inversion), and maps its limits back
#' to the RMSEA scale as \eqn{\sqrt{\lambda / (a \cdot df)}}. When the
#' lower noncentrality limit hits zero (\emph{i.e.}, the data are
#' compatible with a well-fitting model), the lower RMSEA limit is
#' truncated at zero because RMSEA is nonnegative by construction.
#'
#' The multiplier must be the one with which the RMSEA was computed. An
#' RMSEA from lavaan under its default likelihood was computed with
#' \eqn{N}; given to the default \code{"N-1"}, it gives an interval that
#' matches neither convention, off by an amount of order \eqn{1/N}. Two
#' further cases need the model's chi square statistic rather than its
#' RMSEA, and \code{\link{sem_fit_measures}} takes it: an RMSEA truncated at
#' zero, which cannot distinguish the statistics below \eqn{df} that
#' produce it, and the RMSEA of a multiple group model, which carries the
#' factor \eqn{\sqrt{G}} for \eqn{G} groups (Steiger, 1998); this function
#' assumes a single group.
#'
#' The 90 percent CI (rather than the usual 95 percent) is the
#' conventional reporting choice for the RMSEA, from the .05 and .95
#' percentile points of Steiger and Lind (1980) through Browne and Cudeck
#' (1993), because a \eqn{100(1 - 2\alpha)} percent interval decides the
#' two one-sided tests at level \eqn{\alpha}, so that the limits of the 90
#' percent interval settle them at level .05 (MacCallum, Browne, &
#' Sugawara, 1996, p. 136 and Table 1): the hypothesis of close fit,
#' \eqn{H_0: \mathrm{RMSEA} \le \epsilon_0}, is rejected when the
#' lower limit exceeds \eqn{\epsilon_0}, and the hypothesis of not close
#' fit, \eqn{H_0: \mathrm{RMSEA} \ge \epsilon_0}, is rejected, so
#' that fit within \eqn{\epsilon_0} is supported, when the upper limit
#' falls below \eqn{\epsilon_0}. A lower limit of zero likewise says that
#' the test of exact fit does not reject at .05. \code{ci_rmsea} therefore
#' defaults to \code{conf_level = 0.90}, as every RMSEA interval in DMAR
#' does. \code{conf_level = 0.95}, the level of the package's other
#' intervals, gives the 95 percent interval, whose limits decide the
#' one-sided tests at .025. The
#' printed table states the level beneath it. See
#' \code{\link{equivalence_rmsea}} for the test of not close fit with its
#' \emph{p}-value and verdict.
#'
#' @return
#' A 3-row \code{data.frame} with columns \code{term} and \code{value}.
#' The \code{term} values are \code{"ci_lower"} (the lower bound of
#' the confidence interval on the population RMSEA, truncated at zero by
#' definition), \code{"rmsea"} (the observed point estimate), and
#' \code{"ci_upper"} (the upper bound).
#'
#' @examples
#' # 1. The conventional 90 percent CI on RMSEA, the default. Its upper
#' #    limit, 0.052, lies just above 0.05, so the test of not close fit at
#' #    .05 does not reject an RMSEA of 0.05 or more, and fit within 0.05 is
#' #    not established even though the point estimate sits well below it.
#' ci_rmsea(rmsea = .035, df = 40, N = 425)
#'
#' # 2. The 95 percent CI is wider.
#' ci_rmsea(rmsea = .035, df = 40, N = 425, conf_level = .95)
#'
#' # 3. An RMSEA reported by lavaan under its default likelihood, which
#' #    forms the statistic with N.
#' ci_rmsea(rmsea = .0921, df = 24, N = 301, sample_size_convention = "N")
#'
#' # 4. Fewer degrees of freedom and a smaller N: a wider CI.
#' ci_rmsea(rmsea = .055, df = 10, N = 100)
#'
#' @references
#' Browne, M. W., & Cudeck, R. (1993). Alternative ways of assessing
#'   model fit. In K. A. Bollen & J. S. Long (Eds.), \emph{Testing
#'   structural equation models} (pp. 136--162). Sage.
#'
#' Kelley, K., & Lai, K. (2011). Accuracy in parameter estimation for the
#'   root mean square error of approximation: Sample size planning for
#'   narrow confidence intervals.
#'   \emph{Multivariate Behavioral Research, 46}, 1--32.
#'   \doi{10.1080/00273171.2011.543027}
#'
#' MacCallum, R. C., Browne, M. W., & Sugawara, H. M. (1996). Power
#'   analysis and determination of sample size for covariance structure
#'   modeling. \emph{Psychological Methods, 1}(2), 130--149.
#'   \doi{10.1037/1082-989X.1.2.130}
#'
#' Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). \emph{Designing
#'   experiments and analyzing data: A model comparison perspective}
#'   (4th ed.). Routledge.
#'
#' Steiger, J. H. (1998). A note on multiple sample extensions of the RMSEA
#'   fit index. \emph{Structural Equation Modeling, 5}(4), 411--419.
#'   \doi{10.1080/10705519809540115}
#'
#' Steiger, J. H., & Lind, J. C. (1980). \emph{Statistically-based tests for the number of common
#' factors}. Paper presented at the annual Spring meeting of the Psychometric Society, Iowa City, IA.
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso \code{\link{sem_fit_measures}}, \code{\link{equivalence_rmsea}},
#'   \code{\link{ss_aipe_rmsea}}
#'
#' @keywords design
#'
#' @family confidence intervals for effect sizes
#'
#' @export


ci_rmsea <- function(rmsea, df, N, conf_level = 0.90, alpha_lower = NULL, alpha_upper = NULL,
                     sample_size_convention = c("N-1", "N")) {
  .dmar_check_numbers("ci_rmsea", "rmsea", required = "rmsea")
  sample_size_convention <- match.arg(sample_size_convention)
  if (!is.null(alpha_lower) || !is.null(alpha_upper)) {
    if (is.null(alpha_upper) || is.null(alpha_lower)) .dmar_refuse_argument(
      "Both 'alpha_lower' and 'alpha_upper' must be specified",
      arg = c("alpha_lower", "alpha_upper"))
    if (sum(alpha_lower, alpha_upper) >= 1) .dmar_refuse_argument(
      "'alpha_lower' and 'alpha_upper' must sum to less than 1.",
      arg = c("alpha_lower", "alpha_upper"))
    if (alpha_lower >= 1 || alpha_lower < 0) .dmar_refuse_argument(
      "'alpha_lower' must be at least 0 and below 1.", arg = "alpha_lower")
    if (alpha_upper >= 1 || alpha_upper < 0) .dmar_refuse_argument(
      "'alpha_upper' must be at least 0 and below 1.", arg = "alpha_upper")
    conf_level <- NULL
  }

  if (!is.null(conf_level)) {
    # Reached only when the alphas were left NULL: the alpha branch above
    # sets conf_level to NULL whenever the alphas are supplied, so the
    # explicit alphas take precedence over the (default) conf_level.
    if (!is.numeric(conf_level) || length(conf_level) != 1L ||
        !is.finite(conf_level) || conf_level <= 0 || conf_level >= 1) {
      .dmar_refuse_argument(
        "'conf_level' must be a single number between 0 and 1.",
        arg = "conf_level")
    }
    alpha_lower <- (1 - conf_level) / 2
    alpha_upper <- (1 - conf_level) / 2
  }

  if (rmsea < 0) .dmar_refuse_argument(
    "Your RMSEA cannot be less than zero; RMSEA is defined to be nonnegative.",
    arg = NULL)
  # A model with zero degrees of freedom is saturated: its RMSEA is 0/0 and
  # has no interval (the limits were once Inf, the lower one above the
  # estimate). The df need not be whole, since a mean and variance adjusted
  # statistic has fractional degrees of freedom.
  if (!is.numeric(df) || length(df) != 1L || !is.finite(df) || df <= 0) {
    .dmar_refuse_argument(c(
      "'df' must be a single number greater than 0: a model with zero degrees ",
      "of freedom is saturated, and its RMSEA is undefined."), arg = "df")
  }

  # The multiplier a in T = a F: N - 1 by default, N for a statistic
  # formed as lavaan's default normal likelihood forms it. It must be
  # positive, so N must exceed 1 under "N-1" and 0 under "N".
  a <- if (sample_size_convention == "N") N else N - 1
  if (!is.numeric(N) || length(N) != 1L || !is.finite(N) || a <= 0) {
    .dmar_refuse_argument(c(
      "'N' must be a single number greater than ",
      if (sample_size_convention == "N") "0" else "1 (under the default sample_size_convention = \"N-1\")",
      ", so that the multiplier a is positive."), arg = "N")
  }
  chi_sq_statistic <- rmsea^2 * df * a + df

  # The clamp-to-zero warning from ci_nc_chisq fires routinely for
  # well-fitting models where the noncentrality lower limit hits the floor;
  # surfacing it on every RMSEA call would be noise. Suppress it here only.
  chi_sq_ncp_limits <- suppressWarnings(ci_nc_chisq(
    chisq = chi_sq_statistic, conf_level = NULL, df = df,
    alpha_lower = alpha_lower, alpha_upper = alpha_upper
  ))

  Low_Limit <- chi_sq_ncp_limits$value[chi_sq_ncp_limits$term == "ci_lower"]
  Up_Limit  <- chi_sq_ncp_limits$value[chi_sq_ncp_limits$term == "ci_upper"]
  if (is.na(Low_Limit)) Low_Limit <- 0
  Low_Limit <- max(0, Low_Limit)

  if (Low_Limit == 0) Low_RMSEA <- 0
  if (Low_Limit > 0) Low_RMSEA <- (Low_Limit / (df * a))^.5
  if (Low_RMSEA <= 0) message("Note: The lower confidence limit of the noncentrality parameter is at its lower bound, so the lower RMSEA limit is set to 0 based on RMSEA's definition.")
  if (is.na(Low_RMSEA)) {
    Low_RMSEA <- 0
  }
  Low_RMSEA <- max(0, Low_RMSEA)

  Up_RMSEA <- if (is.infinite(Up_Limit)) Inf else (Up_Limit / (df * a))^.5

  output <- data.frame(
    term = c("ci_lower", "rmsea", "ci_upper"),
    value = c(Low_RMSEA, rmsea, Up_RMSEA)
  )
  return(.as_dmar_tbl(output, conf_level = conf_level %||% (1 - alpha_lower - alpha_upper),
                     tail_areas = .tail_areas(alpha_lower, alpha_upper)))
}
