#' Expected Cross-Validation Index (ECVI) for a Covariance-Structure Model
#'
#' The ECVI of Browne and Cudeck (1989) estimates how well a fitted model's
#' implied covariance matrix would fit an independent sample of the same
#' size from the same population. It is the single-sample estimate of the
#' cross-validation discrepancy, so a smaller ECVI indicates a model
#' expected to generalize better; ECVI is most useful for comparing
#' competing models fit to the same data. A confidence interval, derived
#' from the noncentral chi square distribution, accompanies the point
#' estimate.
#'
#' @param fit A fitted \pkg{lavaan} model. Supply this, or the summary
#'   statistics below.
#' @param chisq,df,npar The model chi square, its degrees of freedom, and
#'   the number of free parameters. Used when \code{fit} is not supplied,
#'   so an ECVI can be obtained from a published fit table.
#' @param multiplier The multiplier with which the chi square was formed
#'   from the minimized fit function, which is the divisor of the ECVI
#'   (\eqn{n} in Details): the total sample size \eqn{N} for a statistic
#'   formed as \eqn{N\hat{F}}, and \eqn{N - G} for \eqn{G} groups for a
#'   statistic formed as \eqn{(N - G)\hat{F}}. Given with \code{chisq},
#'   \code{df}, and \code{npar} when \code{fit} is not supplied; with a
#'   fitted model it is taken from the fit.
#' @param conf_level Confidence level for the interval. Defaults to 0.95.
#'
#' @details
#' With \eqn{q} free parameters, \eqn{\mathrm{ECVI} = (\chi^2 + 2q)/n},
#' where \eqn{n} is the multiplier with which the chi square was formed
#' from the minimized fit function (the argument \code{multiplier}),
#' \eqn{\chi^2 = n\hat{F}}, so that
#' \eqn{\mathrm{ECVI} = \hat{F} + 2q/n} (Browne & Cudeck, 1989). With a
#' fitted model the divisor is the one lavaan uses, and the value is the one
#' \code{lavaan::fitMeasures()} reports: \eqn{n = N}, the total sample size,
#' for maximum likelihood under lavaan's default normal likelihood
#' (including MLR, whose unscaled statistic is used, and FIML), and \eqn{n =
#' N - G} for \eqn{G} groups under \code{likelihood = "wishart"} and for the
#' least squares estimators (GLS, WLS, ULS, DWLS). With summary statistics,
#' \code{multiplier} is used as given: pass \eqn{N - G} when the published chi square
#' was formed as \eqn{(N - G)\hat{F}}, as by programs that follow the
#' Wishart convention. Writing \eqn{\hat\lambda = \chi^2 - df} for the
#' estimated noncentrality, the ECVI is \eqn{(\hat\lambda + df + 2q)/n};
#' the confidence interval replaces \eqn{\hat\lambda} by the lower and
#' upper noncentrality limits from \code{\link{ci_nc_chisq}}, the same
#' inversion used for the RMSEA interval (see \code{\link{ci_rmsea}}).
#'
#' The noncentrality limits are at least 0, while the point estimate is the
#' untruncated \eqn{\hat{F} + 2q/n} of Browne and Cudeck (1989), so for a
#' model with \eqn{\chi^2 < df} the estimate lies below the lower limit,
#' \eqn{(df + 2q)/n}. When the chi square lies at or below the upper
#' \eqn{\alpha/2} quantile of the central chi square distribution, for
#' \eqn{\alpha} = 1 - \code{conf_level}, the lower noncentrality limit is 0
#' and the lower limit is that same \eqn{(df + 2q)/n}. When it lies at or
#' below the lower \eqn{\alpha/2} quantile, no noncentrality is consistent
#' with it at that level, and the upper limit is \code{NA}, with a warning.
#' \code{\link{sem_fit_measures}} reports the upper RMSEA limit as 0 in
#' that case, so that it agrees with the test of not close fit, but the
#' ECVI has no such test, and an upper limit of \eqn{(df + 2q)/n} would make
#' an interval of zero width that excludes the estimate.
#'
#' A fitted model must have converged to an admissible solution, as for
#' \code{\link{sem_fit_measures}}: a solution with a negative variance
#' estimate, or with a covariance matrix of the latent variables, of the
#' residuals of the observed variables, or of the structural residuals
#' that is not positive definite, is refused (with an error of class
#' \code{dmar_inadmissible}), since its fit describes estimates outside
#' the parameter space; with auxiliary variables the check is made on the
#' model's own parameters.
#'
#' When the chi square is the likelihood ratio statistic \eqn{2(\ell_s -
#' \ell_M)}, as under lavaan's default normal likelihood, \eqn{N \cdot
#' \mathrm{ECVI} = \mathrm{AIC} + 2\ell_s}, where \eqn{\ell_s} is the
#' maximized log likelihood of the saturated model and \eqn{\ell_M} that
#' of the fitted model. ECVI and AIC then rank models fitted to the same
#' data identically, and their values differ only by the saturated log
#' likelihood constant before division by \eqn{N}. Under the Wishart
#' convention the statistic is \eqn{(N - G)/N} times the likelihood ratio
#' statistic, and the identity holds only approximately. ECVI is reported
#' because its metric (a discrepancy per observation) and its confidence
#' interval are interpretable on their own.
#'
#' @return A \code{data.frame} (class \code{dmar_tbl}) with rows
#'   \code{ecvi}, \code{ci_lower}, and \code{ci_upper} in the
#'   \code{value} column.
#'
#' @references
#' Browne, M. W., & Cudeck, R. (1989). Single sample cross-validation
#'   indices for covariance structures. \emph{Multivariate Behavioral
#'   Research, 24}(4), 445--455.
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso \code{\link{ci_rmsea}}, \code{\link{ci_nc_chisq}}.
#'
#' @family multivariate and latent variable methods
#'
#' @keywords multivariate
#'
#' @examples
#' # From a published fit table (no model object needed).
#' ecvi(chisq = 24.361, df = 8, npar = 13, multiplier = 301)
#'
#' fit <- lavaan::cfa(
#'   "visual =~ t1_visual_perception + t2_cubes + t4_lozenges
#'    verbal =~ t6_paragraph_comprehension + t7_sentence + t9_word_meaning",
#'   data = holzinger_swineford, std.lv = TRUE)
#' ecvi(fit)
#'
#' @export
ecvi <- function(fit = NULL, chisq = NULL, df = NULL, npar = NULL,
                 multiplier = NULL,
                 conf_level = 0.95) {
  if (!is.null(fit)) {
    # lavaan is installed, the fit converged, and its solution is admissible,
    # as for every fit measure of the SEM family.
    .sem_fit_check(fit)
    fm <- lavaan::fitMeasures(fit, c("chisq", "df", "npar"))
    chisq <- unname(fm["chisq"]); df <- unname(fm["df"])
    npar <- unname(fm["npar"])
    # The divisor lavaan uses: N for ML under the normal likelihood, and
    # N less the number of groups under the Wishart likelihood and for the
    # least squares estimators (see .sem_lavaan_convention()).
    multiplier <- .sem_lavaan_convention(fit)$multiplier
  }
  vals <- list(chisq = chisq, df = df, npar = npar, multiplier = multiplier)
  for (nm in names(vals)) {
    v <- vals[[nm]]
    if (!is.numeric(v) || length(v) != 1L || is.na(v) || v < 0) {
      .dmar_refuse_argument(
        sprintf("'%s' must be a single non-negative number (or supply 'fit').", nm),
        arg = nm)
    }
  }
  if (multiplier <= 0) {
    .dmar_refuse_argument(
      "'multiplier' must be positive: it is the divisor of the ECVI.",
      arg = "multiplier")
  }
  if (df < 1) {
    .dmar_refuse_argument(
      "'df' must be at least 1 (the model must be overidentified).", arg = "df")
  }
  if (!is.numeric(conf_level) || length(conf_level) != 1L || is.na(conf_level) ||
      conf_level <= 0 || conf_level >= 1) {
    .dmar_refuse_argument(
      "'conf_level' must be a single number in (0, 1).", arg = "conf_level")
  }

  const <- (df + 2 * npar) / multiplier
  est <- (chisq + 2 * npar) / multiplier               # = (chisq - df)/n + const

  # ci_nc_chisq() warns when it clamps the lower noncentrality limit to 0,
  # which is the right limit for any statistic below the upper alpha/2
  # quantile and needs no warning here, and when no noncentrality reaches
  # the upper limit, which ecvi() reports in its own terms below.
  nc <- withCallingHandlers(
    ci_nc_chisq(chisq = chisq, df = df, conf_level = conf_level,
                verbose = FALSE),
    warning = function(cond) {
      if (grepl("clamped to 0|upper noncentrality limit is undefined",
                conditionMessage(cond)))
        invokeRestart("muffleWarning")
    })
  lambda_l <- nc$value[nc$term == "ci_lower"]
  lambda_u <- nc$value[nc$term == "ci_upper"]
  lower <- lambda_l / multiplier + const
  upper <- lambda_u / multiplier + const
  if (is.na(upper)) {
    .dmar_warn(c(
      "The chi square statistic lies at or below the lower ",
      "(1 - conf_level)/2 quantile of the central chi square ",
      "distribution, so no noncentrality is consistent with it at this ",
      "level and the upper limit of the ECVI is NA."), "result_qualified")
  }

  out <- data.frame(
    term  = c("ecvi", "ci_lower", "ci_upper"),
    value = c(est, lower, upper),
    stringsAsFactors = FALSE, row.names = NULL
  )
  .as_dmar_tbl(out, conf_level = conf_level)
}
