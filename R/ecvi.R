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
#' @param chisq,df,npar,n The model chi square, its degrees of freedom, the
#'   number of free parameters, and the total sample size. Used when
#'   \code{fit} is not supplied, so an ECVI can be obtained from a published
#'   fit table.
#' @param conf_level Confidence level for the interval. Defaults to 0.95.
#'
#' @details
#' With \eqn{q} free parameters and total sample size \eqn{N},
#' \eqn{\mathrm{ECVI} = (\chi^2 + 2q)/N}, the value the
#' \emph{Journal of Statistical Software} reference implementation in
#' \pkg{lavaan} reports. Writing \eqn{\hat\lambda = \chi^2 - df} for the
#' estimated noncentrality, this is \eqn{(\hat\lambda + df + 2q)/N}; the
#' confidence interval replaces \eqn{\hat\lambda} by the lower and upper
#' noncentrality limits from \code{\link{conf_limits_nc_chisq}}, the same
#' inversion used for the RMSEA interval (see \code{\link{ci_rmsea}}). ECVI
#' differs from the AIC only by the constant factor \eqn{N}, so the two rank
#' models identically; ECVI is reported because its metric (a discrepancy
#' per observation) and its confidence interval are interpretable on their
#' own.
#'
#' @return A \code{data.frame} (class \code{dmar_tbl}) with rows
#'   \code{ecvi}, \code{lower_limit}, and \code{upper_limit} in the
#'   \code{value} column.
#'
#' @references
#' Browne, M. W., & Cudeck, R. (1989). Single sample cross-validation
#'   indices for covariance structures. \emph{Multivariate Behavioral
#'   Research, 24}(4), 445--455.
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso \code{\link{ci_rmsea}}, \code{\link{conf_limits_nc_chisq}}.
#'
#' @family multivariate and latent variable methods
#'
#' @keywords multivariate
#'
#' @examples
#' # From a published fit table (no model object needed).
#' ecvi(chisq = 24.361, df = 8, npar = 13, n = 301)
#'
#' fit <- lavaan::cfa(
#'   "visual =~ t1_visual_perception + t2_cubes + t4_lozenges
#'    verbal =~ t6_paragraph_comprehension + t7_sentence + t9_word_meaning",
#'   data = holzinger_swineford, std.lv = TRUE)
#' ecvi(fit)
#'
#' @export
ecvi <- function(fit = NULL, chisq = NULL, df = NULL, npar = NULL, n = NULL,
                 conf_level = 0.95) {
  if (!is.null(fit)) {
    if (!requireNamespace("lavaan", quietly = TRUE)) {
      stop("Package 'lavaan' is required when supplying a fit. Install it ",
           "with install.packages(\"lavaan\").", call. = FALSE)
    }
    if (!inherits(fit, "lavaan")) {
      stop("'fit' must be a lavaan fit object.", call. = FALSE)
    }
    fm <- lavaan::fitMeasures(fit, c("chisq", "df", "npar", "ntotal"))
    chisq <- unname(fm["chisq"]); df <- unname(fm["df"])
    npar <- unname(fm["npar"]); n <- unname(fm["ntotal"])
  }
  vals <- list(chisq = chisq, df = df, npar = npar, n = n)
  for (nm in names(vals)) {
    v <- vals[[nm]]
    if (!is.numeric(v) || length(v) != 1L || is.na(v) || v < 0) {
      stop(sprintf("'%s' must be a single non-negative number (or supply 'fit').", nm),
           call. = FALSE)
    }
  }
  if (df < 1) {
    stop("'df' must be at least 1 (the model must be overidentified).",
         call. = FALSE)
  }
  if (!is.numeric(conf_level) || length(conf_level) != 1L || is.na(conf_level) ||
      conf_level <= 0 || conf_level >= 1) {
    stop("'conf_level' must be a single number in (0, 1).", call. = FALSE)
  }

  const <- (df + 2 * npar) / n
  est <- (chisq + 2 * npar) / n               # = (chisq - df)/n + const

  nc <- conf_limits_nc_chisq(chi_square = chisq, df = df,
                             conf_level = conf_level, verbose = FALSE)
  lambda_l <- nc$value[nc$term == "lower_limit"]
  lambda_u <- nc$value[nc$term == "upper_limit"]
  lower <- lambda_l / n + const
  upper <- lambda_u / n + const

  out <- data.frame(
    term  = c("ecvi", "lower_limit", "upper_limit"),
    value = c(est, lower, upper),
    stringsAsFactors = FALSE, row.names = NULL
  )
  .as_dmar_tbl(out, conf_level = conf_level)
}
