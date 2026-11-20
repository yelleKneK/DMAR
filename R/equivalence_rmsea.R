#' Equivalence Test for the Population RMSEA (Test of Not Close Fit)
#'
#' Tests whether the population root mean square error of approximation
#' (RMSEA) of a structural equation model is below a bound the researcher
#' chooses in advance, so that a model can be endorsed as fitting within a
#' stated tolerance rather than merely not rejected. This is the test of not
#' close fit of MacCallum, Browne, and Sugawara (1996), which Yuan, Chan,
#' Marcoulides, and Bentler (2016) present as an equivalence test: the null
#' hypothesis \eqn{H_0: \epsilon \ge \epsilon_0} that the
#' population RMSEA \eqn{\epsilon} is at least the bound
#' \eqn{\epsilon_0} is tested against \eqn{H_1: \epsilon <
#' \epsilon_0}, and rejecting it supports fit within the bound. It
#' reverses the test of close fit (Browne & Cudeck, 1993), whose null
#' hypothesis \eqn{\epsilon \le \epsilon_0} a well-fitting model
#' merely fails to reject. The function works from the summary statistics a
#' paper reports.
#'
#' @param chisq The model's chi square test statistic \eqn{T}, the
#'   likelihood ratio statistic of a normal theory maximum likelihood fit.
#' @param df The model's degrees of freedom, a positive whole number.
#' @param N The total sample size, summed over groups: a whole number
#'   greater than \code{n_groups}.
#' @param rmsea_bound The bound \eqn{\epsilon_0}: the largest population
#'   RMSEA the researcher is prepared to call fit within tolerance, chosen
#'   before seeing the data. No universal value is assumed.
#' @param n_groups The number of independent groups of a multiple group model,
#'   whose \code{chisq}, \code{df}, and \code{N} are totals over the groups.
#'   Default 1.
#' @param alpha_level The one-sided Type I error rate, in \eqn{(0, 0.5)};
#'   default 0.05. The interval reported has confidence level \eqn{1 -
#'   2\alpha}, 90 percent at the default.
#' @param sample_size_convention The multiplier \eqn{a} in \eqn{T =
#'   a\hat{F}} with which \code{chisq} was formed: \code{"N-1"} (the default),
#'   each group's sample size less one, so \eqn{a = N - G} for \eqn{G}
#'   groups, or \code{"N"}, so \eqn{a = N}. lavaan's default normal
#'   likelihood and Mplus form \eqn{T} with \eqn{N}, so their statistics
#'   take \code{"N"}.
#'
#' @details
#' Write \eqn{T} for the chi square statistic on \eqn{df} degrees of
#' freedom and \eqn{a} for the multiplier with which it was formed. Under a
#' model whose misfit is small relative to sampling error, \eqn{T} is
#' approximately noncentral chi square with noncentrality \eqn{\lambda = a
#' \cdot df \cdot \epsilon^2 / G} (Steiger & Lind, 1980; Browne &
#' Cudeck, 1993; Steiger, 1998, for \eqn{G} groups). At the boundary of the
#' null hypothesis the noncentrality is \eqn{\lambda_0 = a \cdot df \cdot
#' \epsilon_0^2 / G} (row \code{null_noncentrality}), and a small
#' \eqn{T} is evidence against \eqn{H_0}, so the \emph{p}-value is the lower
#' tail, \eqn{P\{\chi^2_{df}(\lambda_0) \le T\}} (MacCallum et al., 1996,
#' Figure 2). The model is declared within tolerance (row
#' \code{equivalent} equal to 1) when the \emph{p}-value is below
#' \code{alpha_level}. The \emph{p}-value and the interval stay accurate for
#' a statistic in the millions, where R's \code{pchisq()} fails (see
#' \code{\link{ci_nc_chisq}}).
#'
#' The decision is the same as comparing an interval with the bound. The
#' interval reported is the two-sided \eqn{100(1 - 2\alpha)} percent
#' interval for the population RMSEA, 90 percent at the default, whose
#' upper limit is the one-sided \eqn{100(1 - \alpha)} percent upper bound;
#' the \emph{p}-value is below \eqn{\alpha} exactly when the upper limit is
#' below \eqn{\epsilon_0} (MacCallum et al., 1996, Table 1). This is why
#' the RMSEA is conventionally reported with a 90 percent interval (Browne &
#' Cudeck, 1993): its upper limit settles the test of not close fit at
#' level .05 and its lower limit the test of close fit. The upper limit is
#' the T-size RMSEA of Yuan et al. (2016, Eq. 4): the data endorse the
#' model at level \eqn{\alpha} for every bound above it, so it is the
#' quantity to report when no bound was fixed in advance. The point
#' estimate is truncated at zero; Steiger (2000) shows that the truncation
#' follows from a basic property of a good estimator and costs little in
#' practice, and the interval carries the information about precision.
#' Yuan et al. use \eqn{a = N - 1}, the default here. The function does not compute their
#' T-size CFI (their Eqs. 6 to 10, which divide \eqn{\alpha} between the
#' model and baseline statistics) or the adjusted cutoff values for judging
#' a T-size index, which their Tables 1 and 2 give as regression
#' coefficients in \eqn{N} and \eqn{df}).
#'
#' When \eqn{T} lies at or below the lower \eqn{\alpha} quantile of the
#' central chi square distribution, no noncentrality is consistent with it
#' at this level; the upper limit is then 0, as Yuan et al. (2016, Eq. 4)
#' define the T-size RMSEA and as lavaan reports it, and the model is
#' within tolerance for every positive bound. That happens with probability
#' \eqn{\alpha} for a model that fits exactly and more often when the
#' normal theory statistic is deflated, for example by data with light
#' tails, so a statistic far below its degrees of freedom deserves a look
#' at the data before it is read as evidence of fit.
#'
#' Failing to reject \eqn{H_0} does not establish poor fit; it says only
#' that fit within the bound has not been shown at this sample size. The
#' test rests on the noncentral chi square approximation for a normal
#' theory statistic. A scaled or robust statistic does not follow that
#' distribution with this noncentrality, so it must not be supplied here;
#' for a fitted model, including a robust (MLR, MLM), categorical (WLSMV),
#' or FIML fit, \code{\link{sem_fit_measures}} reports this test as
#' \code{p_not_close_fit} from lavaan's corrected measures. For a multiple
#' group model the convention is Steiger's (1998): \code{chisq},
#' \code{df}, and \code{N} are totals over the groups, and the RMSEA and
#' its limits carry the factor \eqn{\sqrt{G}}.
#'
#' @return A \code{dmar_tbl} (class \code{dmar_sem_fit}) with columns
#'   \code{term} and \code{value} and the rows \code{rmsea} (the point
#'   estimate \eqn{\sqrt{G \max(T - df, 0)/(a \cdot df)}}),
#'   \code{ci_lower} and \code{ci_upper} (the \eqn{100(1 - 2\alpha)}
#'   percent interval for the population RMSEA), \code{rmsea_bound},
#'   \code{null_noncentrality}, \code{p_value}, \code{equivalent} (1 when
#'   the test rejects \eqn{H_0} at \code{alpha_level}, 0 otherwise),
#'   \code{chi_square}, \code{df}, \code{N}, and \code{n_groups}. The
#'   \code{conf_level} attribute is \eqn{1 - 2\alpha}. The \code{verdict}
#'   attribute states the inference in a sentence, and the \code{sem_fit}
#'   attribute records the sample size convention and multiplier, the
#'   boundary convention that applied, and the lines of the printed
#'   footer.
#'
#' @references
#' Browne, M. W., & Cudeck, R. (1993). Alternative ways of assessing
#'   model fit. In K. A. Bollen & J. S. Long (Eds.), \emph{Testing
#'   structural equation models} (pp. 136--162). Sage.
#'
#' MacCallum, R. C., Browne, M. W., & Sugawara, H. M. (1996). Power
#'   analysis and determination of sample size for covariance structure
#'   modeling. \emph{Psychological Methods, 1}(2), 130--149.
#'   \doi{10.1037/1082-989X.1.2.130}
#'
#' Steiger, J. H. (1998). A note on multiple sample extensions of the RMSEA
#'   fit index. \emph{Structural Equation Modeling, 5}(4), 411--419.
#'   \doi{10.1080/10705519809540115}
#'
#' Steiger, J. H. (2000). Point estimation, hypothesis testing, and
#'   interval estimation using the RMSEA: Some comments and a reply to
#'   Hayduk and Glaser. \emph{Structural Equation Modeling, 7}(2), 149--162.
#'
#' Steiger, J. H., & Lind, J. C. (1980). \emph{Statistically-based tests
#'   for the number of common factors}. Paper presented at the annual Spring
#'   meeting of the Psychometric Society, Iowa City, IA.
#'
#' Yuan, K.-H., Chan, W., Marcoulides, G. A., & Bentler, P. M. (2016).
#'   Assessing structural equation models by equivalence testing with
#'   adjusted fit indexes. \emph{Structural Equation Modeling, 23}(3),
#'   319--330. \doi{10.1080/10705511.2015.1065414}
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso \code{\link{sem_fit_measures}} for the test of not close fit
#'   from a fitted model (robust, categorical, and FIML fits included) beside
#'   the other fit measures, \code{\link{ci_rmsea}} for an interval from a
#'   reported RMSEA, and \code{\link{equivalence_r}} and
#'   \code{\link{equivalence_smd}} for equivalence tests of other effect
#'   sizes.
#'
#' @keywords htest
#'
#' @family hypothesis tests
#'
#' @family equivalence testing
#'
#' @family multivariate and latent variable methods
#'
#' @examples
#' # Yuan et al. (2016, Example 1): nine holzinger_swineford
#' # tests in the Grant-White school (N = 145), the three-factor model with
#' # one cross-loading, T = 28.098 on 23 df. The 90 percent upper limit,
#' # their T-size RMSEA, is .083, so fit within .05 is not established.
#' equivalence_rmsea(chisq = 28.098, df = 23, N = 145, rmsea_bound = 0.05)
#'
#' # A bound of .09 is endorsed at alpha = .05.
#' equivalence_rmsea(chisq = 28.098, df = 23, N = 145, rmsea_bound = 0.09)
#'
#' # The statistic lavaan reports for the two-school model of the same
#' # tests (formed with N, its default), with Steiger's (1998) multiple
#' # group convention.
#' equivalence_rmsea(chisq = 115.8513, df = 48, N = 301, rmsea_bound = 0.10,
#'                   n_groups = 2, sample_size_convention = "N")
#'
#' @export
equivalence_rmsea <- function(chisq, df, N, rmsea_bound, n_groups = 1,
                              alpha_level = 0.05,
                              sample_size_convention = c("N-1", "N")) {
  convention <- match.arg(sample_size_convention)
  .sem_scalar(chisq, "chisq", 0)
  .sem_scalar(df, "df", 0, strict = TRUE, integer = TRUE)
  .sem_scalar(n_groups, "n_groups", 1, integer = TRUE)
  .sem_scalar(N, "N", 1, strict = TRUE, integer = TRUE)
  .sem_scalar(rmsea_bound, "rmsea_bound", 0, strict = TRUE)
  .sem_probability(alpha_level, "alpha_level")
  if (alpha_level >= 0.5) {
    .dmar_refuse_argument(
      "'alpha_level' must be in (0, 0.5).", arg = "alpha_level")
  }
  chisq <- .sem_number(chisq)
  df <- .sem_number(df)
  n_groups <- .sem_number(n_groups)
  N <- .sem_number(N)
  rmsea_bound <- .sem_number(rmsea_bound)
  if (N <= n_groups) {
    .dmar_refuse_argument(
      "'N' is the total sample size over the groups and must exceed 'n_groups'.",
      arg = c("N", "n_groups"))
  }
  a <- .sem_multiplier(N, n_groups, convention)
  lambda0 <- a * df * rmsea_bound^2 / n_groups
  if (!is.finite(lambda0)) {
    .dmar_refuse_argument(
      "The noncentrality at the bound is not finite; check 'rmsea_bound'.",
      arg = "rmsea_bound")
  }
  p_value <- .pchisq_nc(chisq, df = df, ncp = lambda0)
  conf_level <- 1 - 2 * alpha_level
  ncp <- .sem_ncp_bounds(chisq, df, conf_level)
  rmsea <- .sem_rmsea_from_ncp(chisq - df, df, a, n_groups)
  lower <- .sem_rmsea_from_ncp(ncp[["lower"]], df, a, n_groups)
  upper <- .sem_rmsea_from_ncp(ncp[["upper"]], df, a, n_groups)
  equivalent <- as.numeric(p_value < alpha_level)
  vals <- c(rmsea = rmsea, ci_lower = lower, ci_upper = upper,
            rmsea_bound = rmsea_bound, null_noncentrality = lambda0,
            p_value = p_value, equivalent = equivalent,
            chi_square = chisq, df = df, N = N, n_groups = n_groups)
  verdict <- if (equivalent == 1) {
    sprintf(paste0("Not close fit (RMSEA >= %s) is rejected at alpha = %s: ",
                   "the population RMSEA is below the bound (upper limit %s)."),
            format(rmsea_bound), format(alpha_level), format(signif(upper, 3)))
  } else {
    sprintf(paste0("Not close fit (RMSEA >= %s) is not rejected at alpha = %s: ",
                   "fit within the bound is not established (upper limit %s), ",
                   "which is not evidence of poor fit."),
            format(rmsea_bound), format(alpha_level), format(signif(upper, 3)))
  }
  out <- .as_dmar_tbl(data.frame(term = names(vals), value = unname(vals),
                                 stringsAsFactors = FALSE),
                      conf_level = conf_level, p_terms = "p_value",
                      subclass = "dmar_sem_fit")
  attr(out, "verdict") <- verdict
  attr(out, "sem_fit") <- list(
    source = "summary statistics", test = "not close fit",
    null = "population RMSEA >= rmsea_bound",
    alpha_level = alpha_level,
    interval = paste0("two-sided ", .format_conf_pct(conf_level),
                      " percent interval; its upper limit is the one-sided ",
                      .format_conf_pct(1 - alpha_level),
                      " percent upper bound"),
    sample_size_convention = convention, multiplier = a,
    N = N, groups = n_groups,
    boundary = if (stats::pchisq(chisq, df) <= alpha_level) {
      paste0("the statistic lies at or below the lower alpha quantile of the ",
             "central chi square, so the upper limit is 0")
    } else character(0),
    footer = c(verdict,
               sprintf("Sample size: N = %s in %s; multiplier a = %s (sample_size_convention = \"%s\")",
                       format(N),
                       if (n_groups == 1) "1 group" else paste(n_groups, "groups"),
                       format(a), convention)))
  out
}
