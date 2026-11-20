# Proportion of superiority (sometimes called Cohen's U3).
#' Proportion of Superiority (Sometimes Called Cohen's \eqn{U_3})
#'
#' Computes the proportion of the treatment-group population that
#' exceeds the \emph{control-group mean} under bivariate normality
#' with equal variances. This quantity is sometimes called Cohen's
#' \eqn{U_3} (Cohen, 1988). Under those assumptions it equals
#' \eqn{\Phi(\delta)}, where \eqn{\delta} is the population
#' standardized mean difference. When sample sizes are supplied, the
#' CI on the proportion of superiority is constructed by transforming
#' the noncentral \emph{t} CI on Cohen's \emph{d} via \eqn{\Phi},
#' which is monotone and therefore preserves coverage exactly.
#'
#' The proportion of superiority is one of three "U" indices Cohen
#' (1988) defined; the other two (\eqn{U_1}, the proportion of
#' non-overlap, and \eqn{U_2}, the proportion of either population
#' that exceeds the same percentile in the other) can be derived from
#' it directly: with \eqn{U_3 = \Phi(\delta)} for the proportion of
#' superiority, Cohen's \eqn{U_2 = \Phi(\delta/2)} and
#' \eqn{U_1 = (2 \cdot U_2 - 1) / U_2} (Cohen, 1988, Table 2.2.1).
#'
#' @param smd Sample standardized mean difference (Cohen's \emph{d}).
#'   Numeric scalar.
#' @param n_1,n_2 Group sample sizes (required if a CI is wanted).
#' @param conf_level Confidence level for the CI. Default \code{0.95}.
#' @param smd_ci_lower,smd_ci_upper Optional pre-computed CI limits on
#'   \emph{d}; when supplied directly, the function skips the
#'   noncentral \emph{t} step and just transforms these limits.
#'
#' @return A \code{data.frame} with rows for \emph{d} (\code{smd}), the
#'   proportion of superiority (\code{proportion_of_superiority}), and
#'   (when a CI is constructable) the confidence limits on \emph{d}
#'   (\code{smd_ci_lower}, \code{smd_ci_upper}) and on the proportion of
#'   superiority (\code{ci_lower}, \code{ci_upper}). With an interval the
#'   table carries the attributes \code{conf_level}; \code{ci_estimate},
#'   which names \code{proportion_of_superiority} as the row that
#'   \code{ci_lower} and \code{ci_upper} belong to, so that
#'   \code{\link{results_sentence}} reports the proportion of superiority
#'   with them; and
#'   \code{ci_source}, which names \code{smd} as the estimate whose
#'   interval that interval is computed from, so that \code{tidy()} gives
#'   the confidence level on the rows of the proportion of superiority
#'   and its limits alone.
#'
#' @details
#' \strong{Why this rather than \code{cles}.} The proportion of
#' superiority answers the question "what fraction of the treatment
#' population exceeds the \emph{control-group mean},'' whereas
#' \code{\link{cles}} answers "what fraction of randomly drawn pairs
#' favor the treatment over the control.'' Both are unitless
#' probability-scale summaries of a Cohen's-\emph{d} difference, but
#' the proportion of superiority is marginal while CLES is paired.
#' Specifically, \eqn{\Phi(d)} versus \eqn{\Phi(d/\sqrt{2})}; for
#' \eqn{d = 0.5}, the proportion of superiority is 0.69 and CLES is
#' 0.64.
#'
#' \strong{CI construction.} Because \eqn{\Phi(\cdot)} is monotone,
#' the CI on the proportion of superiority is just
#' \eqn{[\Phi(d_L),\, \Phi(d_U)]} where \eqn{[d_L,\, d_U]} is the
#' noncentral \emph{t} CI on \emph{d} from \code{\link{ci_smd}}.
#'
#' @references
#' Cohen, J. (1988). \emph{Statistical power analysis for the behavioral
#'   sciences} (2nd ed.). Hillsdale, NJ: Lawrence Erlbaum.
#'   (See Section 2.2 for the \eqn{U_1}, \eqn{U_2}, and \eqn{U_3} indices.)
#'
#' Hedges, L. V., & Olkin, I. (1985). \emph{Statistical methods for
#'   meta-analysis}. Academic Press.
#'
#' Kelley, K. (2007). Confidence intervals for standardized effect sizes:
#'   Theory, application, and implementation. \emph{Journal of Statistical
#'   Software, 20}(8), 1--24. (The noncentral \emph{t} interval on the
#'   standardized mean difference that is transformed here.)
#'   \doi{10.18637/jss.v020.i08}
#'
#' Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). \emph{Designing
#'   experiments and analyzing data: A model comparison perspective}
#'   (4th ed.). Routledge.
#'
#' @seealso \code{\link{cles}}, \code{\link{nnt_from_smd}},
#'   \code{\link{smd}}, \code{\link{ci_smd}},
#'   \code{\link{probability_of_superiority_paired}}
#'
#' @examples
#' # 1. Proportion of superiority at three reference d values:
#' proportion_of_superiority(smd = 0.2)
#' proportion_of_superiority(smd = 0.5)
#' proportion_of_superiority(smd = 0.8)
#'
#' # 2. With a noncentral t CI from sample sizes:
#' proportion_of_superiority(smd = 0.5, n_1 = 50, n_2 = 50, conf_level = 0.95)
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @concept Cohen's U3
#'
#' @keywords htest
#'
#' @family effect size estimates
#'
#' @export

proportion_of_superiority <- function(smd,
                                      n_1 = NULL, n_2 = NULL,
                                      conf_level = 0.95,
                                      smd_ci_lower = NULL, smd_ci_upper = NULL) {
  .dmar_check_conf_level(conf_level, null_ok = TRUE)
  if (!is.numeric(smd) || length(smd) != 1L)
    .dmar_refuse_argument("'smd' must be a single numeric value.", arg = "smd")

  ps_point <- stats::pnorm(smd)
  out <- data.frame(
    term  = c("smd", "proportion_of_superiority"),
    value = c(smd, ps_point),
    stringsAsFactors = FALSE,
    row.names = NULL
  )

  if (is.null(smd_ci_lower) && is.null(smd_ci_upper) &&
      !is.null(n_1) && !is.null(n_2)) {
    ncp_t <- smd * sqrt(as.numeric(n_1) * n_2 / (n_1 + n_2))
    ci  <- ci_smd(ncp = ncp_t, n_1 = n_1, n_2 = n_2, conf_level = conf_level)
    smd_ci_lower <- ci$value[ci$term == "ci_lower"]
    smd_ci_upper <- ci$value[ci$term == "ci_upper"]
  }

  has_ci <- !is.null(smd_ci_lower) && !is.null(smd_ci_upper)
  if (has_ci) {
    out <- rbind(
      out,
      data.frame(term  = c("smd_ci_lower", "smd_ci_upper",
                           "ci_lower", "ci_upper"),
                 value = c(smd_ci_lower, smd_ci_upper,
                           stats::pnorm(smd_ci_lower),
                           stats::pnorm(smd_ci_upper)),
                 stringsAsFactors = FALSE,
                 row.names = NULL)
    )
  }

  # Only label a coverage level when an interval was actually produced.
  out <- .as_dmar_tbl(out, conf_level = if (has_ci) conf_level else NULL)
  # The ci_lower and ci_upper rows are the interval for the
  # proportion of superiority, listed after the standardized mean
  # difference it is computed from, so the table names that row as the
  # one its interval belongs to; results_sentence() reads the name.
  # Through 1.0.1 the sentence paired the standardized mean difference,
  # the first row, with this interval.
  if (has_ci) attr(out, "ci_estimate") <- "proportion_of_superiority"
  # The interval for d it is computed from is echoed in rows named for
  # d (smd_ci_lower, smd_ci_upper); the table names d as that source, so
  # that tidy() marks the level on the rows of the proportion of
  # superiority and its limits alone, as when the limits on d were named
  # smd_lower and smd_upper.
  if (has_ci) attr(out, "ci_source") <- "smd"
  out
}
