# Number needed to treat (NNT) from Cohen's d.
#' Number Needed to Treat (NNT) From Cohen's \emph{d}
#'
#' Converts a standardized mean difference (Cohen's \emph{d}) into the
#' number needed to treat (NNT) using the Kraemer-Kupfer (2006) framework,
#' which connects \emph{d} to the success-rate difference (SRD) under the
#' assumption of a continuous, normally distributed outcome with equal
#' variances and a hypothetical median-cut criterion for treatment success.
#' When a confidence interval on \emph{d} is supplied (or noncentrality
#' parameters / sample sizes are provided so that one can be constructed),
#' the bounds are propagated through the same conversion to give a CI on
#' the NNT.
#'
#' @param smd Sample standardized mean difference (Cohen's \emph{d}); a
#'   numeric scalar. Positive values correspond to the treatment group
#'   exceeding the control group.
#' @param n_1,n_2 Sample sizes in the treatment and control groups; both
#'   required when a noncentral \emph{t}-based CI on the NNT is desired.
#' @param conf_level Confidence level for the CI on the NNT (when \code{n_1}
#'   and \code{n_2} are supplied). Default \code{0.95}.
#' @param smd_lower,smd_upper Optional pre-computed confidence limits on
#'   \emph{d}. If supplied, these are used directly to propagate the
#'   interval through the SRD-to-NNT map and the noncentral \emph{t}
#'   computation is skipped.
#'
#' @return A \code{data.frame} with rows for the success-rate
#'   difference (\code{srd}), the point estimate of NNT (\code{nnt}), and
#'   (when an interval is constructable) the lower and upper NNT limits.
#'   When the lower CI on \emph{d} is exactly zero, the corresponding NNT
#'   bound is reported as \code{Inf}; when it is negative, that bound is a
#'   finite negative value (the NNT to harm), so a CI on \emph{d} that
#'   spans zero yields an NNT interval passing through the infinite point
#'   that separates benefit from harm.
#'
#' @details
#' \strong{The conversion.} Under bivariate normality with equal variances,
#' Kraemer & Kupfer (2006) showed that the proportion of times a randomly
#' drawn treatment-group observation exceeds a randomly drawn control-group
#' observation is
#' \deqn{p \;=\; \Pr(Y_T > Y_C) \;=\; \Phi\!\bigl(d / \sqrt{2}\bigr),}
#' from which the success-rate difference (their effect size) is
#' \deqn{\mathrm{SRD} \;=\; 2 p - 1 \;=\; 2 \Phi\!\bigl(d / \sqrt{2}\bigr) - 1,}
#' and the number needed to treat is its reciprocal,
#' \deqn{\mathrm{NNT} \;=\; 1 / \mathrm{SRD}.}
#' Larger \emph{d} produces smaller NNT; \eqn{d = 0} produces
#' \eqn{\mathrm{NNT} = \infty} (no advantage). The conversion is monotone,
#' so the SRD/NNT confidence interval is obtained by applying the same
#' transformation to the endpoints of the CI on \emph{d}; the lower NNT
#' limit comes from the \emph{upper} \emph{d} limit and vice versa
#' (Furukawa & Leucht, 2011).
#'
#' \strong{When NNT becomes infinite or negative.} If the lower CI on
#' \emph{d} is exactly zero, the corresponding upper NNT bound is
#' \code{Inf}: the data do not exclude the possibility that the treatment
#' produces no advantage (or even harm). Negative values of \emph{d} are
#' allowed; the function returns negative NNT values which are conventionally
#' read as the NNT to \emph{harm}.
#'
#' \strong{Assumption check.} The Kraemer-Kupfer conversion assumes a
#' continuous, normally distributed outcome with equal variances across
#' groups. For skewed outcomes, ordinal outcomes, or unequal variances, the
#' empirical common-language effect size \code{\link{cles}} or the
#' Vargha-Delaney \emph{A} statistic is more defensible. Furukawa & Leucht
#' (2011) compare four methods and recommend the Kraemer-Kupfer formula as
#' the most accurate under normality.
#'
#' @references
#' Furukawa, T. A., & Leucht, S. (2011). How to obtain NNT from Cohen's
#'   \emph{d}: Comparison of four methods. \emph{PLoS ONE, 6}(4), e19070.
#'   \doi{10.1371/journal.pone.0019070}
#'
#' Kelley, K. (2007). Confidence intervals for standardized effect sizes:
#'   Theory, application, and implementation. \emph{Journal of Statistical
#'   Software, 20}(8), 1--24. (The noncentral \emph{t} interval on the
#'   standardized mean difference that is mapped to the NNT here.)
#'   \doi{10.18637/jss.v020.i08}
#'
#' Kraemer, H. C., & Kupfer, D. J. (2006). Size of treatment effects and
#'   their importance to clinical research and practice. \emph{Biological
#'   Psychiatry, 59}(11), 990--996. \doi{10.1016/j.biopsych.2005.09.014}
#'
#' @seealso \code{\link{smd}}, \code{\link{ci_smd}}, \code{\link{cles}}
#'
#' @examples
#' # 1. Point estimate only.
#' nnt_from_smd(smd = 0.5)
#'
#' # 2. With a noncentral t CI from sample sizes:
#' nnt_from_smd(smd = 0.5, n_1 = 50, n_2 = 50, conf_level = 0.95)
#'
#' # 3. With pre-computed CI on d:
#' nnt_from_smd(smd = 0.5, smd_lower = 0.20, smd_upper = 0.80)
#'
#' # 4. Lower d below zero: upper NNT bound is a finite negative value (NNT to harm).
#' nnt_from_smd(smd = 0.3, smd_lower = -0.10, smd_upper = 0.70)
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @keywords htest
#'
#' @family effect size estimates
#'
#' @export

nnt_from_smd <- function(smd,
                         n_1 = NULL, n_2 = NULL,
                         conf_level = 0.95,
                         smd_lower = NULL, smd_upper = NULL) {
  if (!is.numeric(smd) || length(smd) != 1L)
    stop("'smd' must be a single numeric value.")

  # Point estimate.
  srd <- 2 * stats::pnorm(smd / sqrt(2)) - 1
  nnt <- if (srd != 0) 1 / srd else Inf

  out <- data.frame(
    term  = c("smd", "srd", "nnt"),
    value = c(smd, srd, nnt),
    stringsAsFactors = FALSE,
    row.names = NULL
  )

  # If both CI limits supplied directly, propagate; otherwise, compute
  # via ci_smd() when sample sizes are available.
  if (is.null(smd_lower) && is.null(smd_upper) &&
      !is.null(n_1) && !is.null(n_2)) {
    nct <- smd * sqrt(n_1 * n_2 / (n_1 + n_2))
    ci  <- ci_smd(ncp = nct, n_1 = n_1, n_2 = n_2, conf_level = conf_level)
    smd_lower <- ci$value[ci$term == "lower_limit"]
    smd_upper <- ci$value[ci$term == "upper_limit"]
  }

  has_ci <- !is.null(smd_lower) && !is.null(smd_upper)
  if (has_ci) {
    # NNT is a monotone-decreasing function of d, so flip the bounds:
    srd_lo <- 2 * stats::pnorm(smd_lower / sqrt(2)) - 1
    srd_hi <- 2 * stats::pnorm(smd_upper / sqrt(2)) - 1
    nnt_hi <- if (srd_lo != 0) 1 / srd_lo else Inf
    nnt_lo <- if (srd_hi != 0) 1 / srd_hi else Inf

    out <- rbind(
      out,
      data.frame(term  = c("smd_lower", "smd_upper",
                           "nnt_lower", "nnt_upper"),
                 value = c(smd_lower, smd_upper, nnt_lo, nnt_hi),
                 stringsAsFactors = FALSE,
                 row.names = NULL)
    )
  }

  # Only label a coverage level when an interval was actually produced.
  .as_dmar_tbl(out, conf_level = if (has_ci) conf_level else NULL)
}
