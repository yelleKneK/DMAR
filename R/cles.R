# Common-language effect size (McGraw & Wong, 1992).
#' Common-Language Effect Size (McGraw & Wong, 1992)
#'
#' Computes the common-language (CL) effect size for two independent
#' groups, defined as the probability that a randomly drawn observation
#' from group 1 exceeds a randomly drawn observation from group 2 under
#' bivariate normality with equal variances:
#' \deqn{\mathrm{CL} \;=\; \Pr(Y_1 > Y_2) \;=\;
#'   \Phi\!\bigl(\delta / \sqrt{2}\bigr),}
#' where \eqn{\delta} is the population standardized mean difference and
#' \eqn{\Phi} is the standard normal cumulative distribution function.
#' When sample sizes are supplied, the confidence interval on CL is
#' constructed by transforming the noncentral \emph{t}-based CI on
#' Cohen's \emph{d} (Steiger & Fouladi, 1997; Kelley, 2007) through
#' \eqn{\Phi(\cdot / \sqrt{2})}, which is monotone-increasing so the
#' coverage probability is preserved exactly. This is preferred over the
#' normal-approximation CI on CL commonly seen in applied work
#' (Brooks, Dalal, & Nolan, 2014).
#'
#' The common-language idea extends to other effect sizes; the common
#' language effect size for correlations is developed by Liu, Carlson, and
#' Kelley (2019).
#'
#' @param smd Sample standardized mean difference (Cohen's \emph{d}); a
#'   numeric scalar. Positive means group 1 exceeds group 2.
#' @param n_1,n_2 Sample sizes in the two groups; both required when a
#'   confidence interval on CL is desired.
#' @param conf_level Confidence level for the CI. Default \code{0.95}.
#' @param smd_lower,smd_upper Optional pre-computed confidence limits on
#'   \emph{d}. If supplied, these are used directly and the noncentral
#'   computation is skipped.
#'
#' @return A \code{data.frame} with rows for the point estimate
#'   (\code{cl}) and, when sample sizes are supplied, the lower and upper
#'   CI limits. The \emph{d}-equivalent of each row is also reported for
#'   transparency.
#'
#' @details
#' \strong{Background.} McGraw & Wong (1992) introduced the CL effect size
#' to make Cohen's \emph{d} more interpretable: instead of "the means
#' differ by 0.5 SD," one can say "in 64% of randomly drawn pairs, the
#' treated person scores higher than the control person." Under
#' bivariate normality with equal variances, the population probability
#' \eqn{\Pr(Y_1 > Y_2)} equals \eqn{\Phi(\delta/\sqrt{2})}, where
#' \eqn{\delta = (\mu_1 - \mu_2)/\sigma} (McGraw & Wong, 1992).
#'
#' \strong{Connection to other measures.} CL is identical to the AUC
#' (Area Under the Curve) interpretation of \emph{d} in receiver-operating
#' analysis. Vargha & Delaney (2000) generalized CL to the nonparametric
#' setting (their \emph{A} measure) by replacing the population \emph{p}
#' with its empirical Mann-Whitney estimate; under bivariate normality
#' the two coincide. The success-rate-difference and number-needed-to-
#' treat scales (Kraemer & Kupfer, 2006; see \code{\link{nnt_from_smd}})
#' are linear transformations of CL: \eqn{\mathrm{SRD} = 2 \mathrm{CL} - 1}.
#'
#' \strong{Confidence interval construction.} Because \eqn{\Phi(\cdot/\sqrt{2})}
#' is monotone-increasing, the CI on CL is obtained by transforming the
#' CI on \emph{d}: \eqn{[\Phi(d_L/\sqrt 2),\, \Phi(d_U/\sqrt 2)]}. This is
#' an exact-coverage interval (under the noncentral \emph{t} sampling model)
#' and is more accurate than the normal-approximation CI on CL that uses
#' a Wald-style variance for \eqn{\hat p} (Brooks, Dalal, & Nolan, 2014).
#'
#' @references
#' Brooks, M. E., Dalal, D. K., & Nolan, K. P. (2014). Are common language
#'   effect sizes easier to understand than traditional effect sizes?
#'   \emph{Journal of Applied Psychology, 99}(2), 332--340.
#'   \doi{10.1037/a0034745}
#'
#' Kelley, K. (2007). Confidence intervals for standardized
#'   effect sizes: Theory, application, and implementation.
#'   \emph{Journal of Statistical Software, 20}(8), 1--24.
#'   \doi{10.18637/jss.v020.i08}
#'
#' Kraemer, H. C., & Kupfer, D. J. (2006). Size of treatment effects and
#'   their importance to clinical research and practice. \emph{Biological
#'   Psychiatry, 59}(11), 990--996. \doi{10.1016/j.biopsych.2005.09.014}
#'
#' Liu, X. S., Carlson, R., & Kelley, K. (2019). Common language effect
#'   size for correlations. \emph{The Journal of General Psychology,
#'   146}(3), 325--338. \doi{10.1080/00221309.2019.1585321}
#'
#' McGraw, K. O., & Wong, S. P. (1992). A common language effect size
#'   statistic. \emph{Psychological Bulletin, 111}(2), 361--365.
#'   \doi{10.1037/0033-2909.111.2.361}
#'
#' Steiger, J. H., & Fouladi, R. T. (1997). Noncentrality interval
#'   estimation and the evaluation of statistical methods. In L. L. Harlow,
#'   S. A. Mulaik, & J. H. Steiger (Eds.), \emph{What if there were no
#'   significance tests?} (pp. 221--257). Mahwah, NJ: Lawrence Erlbaum.
#'
#' Vargha, A., & Delaney, H. D. (2000). A critique and improvement of the
#'   CL common language effect size statistics of McGraw and Wong.
#'   \emph{Journal of Educational and Behavioral Statistics, 25}(2),
#'   101--132. \doi{10.3102/10769986025002101}
#'
#' @seealso \code{\link{smd}}, \code{\link{ci_smd}},
#'   \code{\link{nnt_from_smd}}
#'
#' @examples
#' # 1. Point estimate only.
#' cles(smd = 0.5)
#'
#' # 2. With a noncentral t CI from sample sizes (preferred):
#' cles(smd = 0.5, n_1 = 50, n_2 = 50, conf_level = 0.95)
#'
#' # 3. With a pre-computed CI on d:
#' cles(smd = 0.5, smd_lower = 0.20, smd_upper = 0.80)
#'
#' # 4. CL at three reference d values:
#' cles(smd = 0.2)
#' cles(smd = 0.5)
#' cles(smd = 0.8)
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @concept probability of superiority
#'
#' @keywords htest
#'
#' @family effect size estimates
#'
#' @export

cles <- function(smd,
                 n_1 = NULL, n_2 = NULL,
                 conf_level = 0.95,
                 smd_lower = NULL, smd_upper = NULL) {
  if (!is.numeric(smd) || length(smd) != 1L)
    stop("'smd' must be a single numeric value.")

  cl <- stats::pnorm(smd / sqrt(2))
  out <- data.frame(
    term  = c("smd", "cl"),
    value = c(smd, cl),
    stringsAsFactors = FALSE,
    row.names = NULL
  )

  # If both CI limits supplied directly, propagate; otherwise compute via
  # ci_smd() when sample sizes are available.
  if (is.null(smd_lower) && is.null(smd_upper) &&
      !is.null(n_1) && !is.null(n_2)) {
    nct <- smd * sqrt(n_1 * n_2 / (n_1 + n_2))
    ci  <- ci_smd(ncp = nct, n_1 = n_1, n_2 = n_2, conf_level = conf_level)
    smd_lower <- ci$value[ci$term == "lower_limit"]
    smd_upper <- ci$value[ci$term == "upper_limit"]
  }

  has_ci <- !is.null(smd_lower) && !is.null(smd_upper)
  if (has_ci) {
    cl_lo <- stats::pnorm(smd_lower / sqrt(2))
    cl_hi <- stats::pnorm(smd_upper / sqrt(2))
    out <- rbind(
      out,
      data.frame(term  = c("smd_lower", "smd_upper",
                           "cl_lower",  "cl_upper"),
                 value = c(smd_lower, smd_upper, cl_lo, cl_hi),
                 stringsAsFactors = FALSE,
                 row.names = NULL)
    )
  }

  # Attach the confidence level only when limits were actually produced, so a
  # point-estimate-only call does not print a misleading coverage footer.
  .as_dmar_tbl(out, conf_level = if (has_ci) conf_level else NULL)
}
