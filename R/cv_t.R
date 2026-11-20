# Critical value(s) for a \emph{t}-distribution
#' Provides the Critical Value(s) for a \emph{t}-distribution
#'
#' @param alpha_level Type I error rate (i.e., the false positive rate).
#' @param df The number of degrees of freedom (a positive number)
#' @param alternative The type of alternative hypothesis of interest.
#' @param alpha_lower The error rate on the lower (negative) side of the distribution.
#' @param alpha_upper The error rate on the upper (positive) side of the distribution.
#' @param ncp The noncentral parameter (if zero, the default, it is the central \emph{t}-distribution).
#' @param verbose Provides extra information about areas under the curve.
#'
#' @return Returns the critical value(s), based on the input specifications, in a output style.
#'
#' @details Though a noncentral parameter can be included, that would not be done for a standard null hypothesis significance test.
#'
#' With a nonzero \code{ncp} the critical values are quantiles of the
#' noncentral \emph{t} distribution, computed exactly. R's \code{qt()} and
#' \code{pt()} with \code{ncp} switch to a normal approximation above a
#' noncentrality of 37.62 (and above 33 with more than 2,000 degrees of
#' freedom), where they can be far off: the .975 quantile with
#' \code{ncp = 50} and \code{df = 3} is 186.49, and \code{qt()} gives
#' 429.18. There the quantiles and the tail areas are computed by
#' integrating over the sampling distribution of the standard deviation;
#' inside that range they are R's. With fewer than one degree of freedom
#' (\code{df} may be any positive number) \code{pt()} is not accurate inside
#' that range either. At the exact quantiles of probability .0005 to .9995,
#' with noncentralities up to 37.5 in absolute value, the tail probability
#' it gave was off by as much as .14 with \code{df = 0.05}, .07 with 0.1,
#' .01 with 0.2, .0013 with 0.3, \eqn{3.8 \times 10^{-5}} with 0.5, and
#' \eqn{2 \times 10^{-11}} with 0.8, so with \code{df < 1} every noncentral
#' quantile and tail area is integrated. The critical values then differ
#' from \code{qt()}'s by orders of magnitude with \code{df} of 0.3 or less
#' (and by as much as 15 percent with 0.5): with \code{df = 0.1},
#' \code{ncp = 0.5}, and \code{alpha_level = .05} they are
#' \eqn{-1.0428 \times 10^{10}} and \eqn{5.7095 \times 10^{13}}, where
#' \code{qt()} gives -33554432 and 305284064. With so few degrees of
#' freedom the tails
#' are heavy enough that a critical value can lie beyond the largest number
#' R represents (about \eqn{1.8 \times 10^{308}}); it is then returned as
#' \code{-Inf} or \code{Inf}, as \code{qt()} returns it with
#' \code{ncp = 0} (with \code{df = 0.001} one or both of the critical
#' values are infinite at every \code{alpha_level} from .000001 to .4 and
#' noncentrality from -100 to 100 tried).
#'
#' @references
#' Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). \emph{Designing
#'   experiments and analyzing data: A model comparison perspective}
#'   (4th ed.). Routledge.
#'
#' @examples
#' # A basic call for finding critical values with equal area in the two tails.
#' cv_t(alpha_level = .05, df = 13)
#'
#' # A basic call for a single-sided confidence interval (for "a greater than" alternative hypothesis)
#' cv_t(alpha_level = .05, df = 13, alternative = "greater")
#'
#' # A single-sided confidence interval (for "a greater than" alternative hypothesis); simple output.
#' cv_t(alpha_lower = 0, alpha_upper = .05, df = 13, verbose = FALSE)
#'
#' # For a nonsymmetric 95% confidence interval.
#' cv_t(alpha_lower = .01, alpha_upper = .04, df = 13)
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @family critical values
#'
#' @export
#' @import stats

cv_t <- function(alpha_level, df, alternative = "not_equal", alpha_lower, alpha_upper, ncp = 0, verbose = TRUE) {
  if (!missing(alpha_level)) .dmar_check_alpha_level(alpha_level)
  .dmar_check_numbers("cv_t", "df", required = "df")
  if (missing(df)) .dmar_refuse_argument(
    "You must specify the degrees of freedom (i.e., 'df'), which must be a positive number.",
    arg = "df")
  if (!df > 0) .dmar_refuse_argument(
    "You must specify a positive value for the degrees of freedom.", arg = NULL)

  if (missing(alpha_level) && missing(alpha_lower) && missing(alpha_upper)) .dmar_refuse_argument(
    "You must specify either 'alpha_level' or 'alpha_lower' and 'alpha_upper'.",
    arg = c("alpha_level", "alpha_lower", "alpha_upper"))

  if (missing(alpha_level)) alpha_level <- NULL
  if (missing(alpha_lower)) alpha_lower <- NULL
  if (missing(alpha_upper)) alpha_upper <- NULL

  if (!is.null(alpha_level)) {
    if (!is.null(alpha_lower)) .dmar_refuse_argument(
      "You have specified both 'alpha_level' and 'alpha_lower'; only use one approach.",
      arg = c("alpha_level", "alpha_lower"))
    if (!is.null(alpha_upper)) .dmar_refuse_argument(
      "You have specified both 'alpha_level' and 'alpha_upper'; only use one approach.",
      arg = c("alpha_level", "alpha_upper"))
    if (alpha_level <= 0 || alpha_level >= 1) .dmar_refuse_argument(
      "Specify 'alpha_level' to be greater than zero and less than 1.",
      arg = "alpha_level")

    if (tolower(alternative) %in% c("greater than", "greater-than", "greater", "greater.than", "gt", "g", ">", ">=")) {
      alpha_lower <- 0
      alpha_upper <- alpha_level
    }

    if (tolower(alternative) %in% c("less than", "less-than", "less", "lesser", "less.than", "lt", "l", "<", "<=")) {
      alpha_lower <- alpha_level
      alpha_upper <- 0
    }

    if (tolower(alternative) %in% c("ne", "2s", "not equal", "two.sided", "two sided", "two-sided", "!=", "not_equal")) {
      alpha_lower <- alpha_upper <- alpha_level / 2
    }
  }

  if (is.null(alpha_level)) {
    if (is.null(alpha_lower)) .dmar_refuse_argument(
      "With 'alpha_level=NULL' you need to specify 'alpha_lower' (and alpha_upper)",
      arg = "alpha_lower")
    if (is.null(alpha_upper)) .dmar_refuse_argument(
      "With 'alpha_level=NULL' you need to specify 'alpha_upper' (and alpha_lower)",
      arg = "alpha_upper")
  }

  if (alpha_lower < 0 || alpha_lower >= .5) .dmar_refuse_argument(
    "Specify 'alpha_lower' to be greater than or equal to zero but less than .50",
    arg = "alpha_lower")
  if (alpha_upper < 0 || alpha_upper >= .5) .dmar_refuse_argument(
    "Specify 'alpha_upper' to be greater than or equal to zero but less than .50.",
    arg = "alpha_upper")

  # Outputs. The quantiles and tail areas are exact (.nc_t_quantile() and
  # .nc_t_prob(), R/nc_t_internals.R): R's qt() and pt() with ncp are a
  # normal approximation above a noncentrality of 37.62, which through
  # 3338679 gave 429.18 for the .975 quantile with ncp = 50 and df = 3
  # (exactly 186.49). With ncp = 0 they are qt() and pt().
  term <- c("lower_cv", "upper_cv")
  value <- c(.nc_t_quantile(alpha_lower, df, ncp),
             .nc_t_quantile(1 - alpha_upper, df, ncp))
  area_less <- c(.nc_t_prob(value[1], df, ncp), .nc_t_prob(value[2], df, ncp))
  area_greater <- c(.nc_t_prob(value[1], df, ncp, lower_tail = FALSE),
                    .nc_t_prob(value[2], df, ncp, lower_tail = FALSE))

  if (verbose == TRUE) {
    return(.as_dmar_tbl(data.frame(term, value, area_less, area_greater)))
  }
  if (verbose == FALSE) {
    return(.as_dmar_tbl(data.frame(term, value)))
  }
}

