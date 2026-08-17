#' Random Effects Meta-Analysis of Correlations
#'
#' Pools correlations across independent studies on the Fisher's \emph{Z}
#' scale and reports the results back in the correlation metric. Optionally,
#' each study's correlation is first corrected for attenuation due to
#' measurement error in either or both variables (the Spearman correction of
#' \code{\link{correction_for_attenuation}}, the basic artifact correction of
#' Hunter and Schmidt's psychometric meta-analysis), using reliabilities you
#' supply, for example from the \code{\link{reliability}} family. That
#' combination, synthesis connected to an actual reliability toolkit, is the
#' measurement-aware path: the pooled quantity is then the construct-level
#' correlation rather than the attenuated observed one.
#'
#' @param r Numeric vector of observed correlations, one per study, each in
#'   (-1, 1).
#' @param n Per-study sample sizes (integer, at least 4).
#' @param reliability_x,reliability_y Optional per-study reliabilities in
#'   (0, 1] for the two measured variables; a single value is recycled
#'   across studies. When either is supplied, each correlation is
#'   disattenuated by \eqn{r_i / \sqrt{\rho_{xx,i}\, \rho_{yy,i}}} before
#'   pooling (a reliability left \code{NULL} is treated as 1). The
#'   reliabilities are treated as known.
#' @param method,hartung_knapp,conf_level Passed to \code{\link{meta_es}}.
#'
#' @details
#' Pooling uses \eqn{z_i = \mathrm{atanh}(r_i)} with sampling variance
#' \eqn{1 / (n_i - 3)}; the pooled estimate, its confidence limits, and the
#' prediction interval are transformed back through \eqn{\tanh}. The
#' heterogeneity quantities (tau, tau-squared, I-squared, H-squared, Q)
#' remain on the Fisher's \emph{Z} scale, where the model lives; tau is
#' therefore the between-study standard deviation of the \emph{z}-scale
#' correlations.
#'
#' When corrections are applied, the corrected correlation's variance is
#' computed from its own \eqn{n_i} on the \emph{z} scale, the conventional
#' simple treatment when reliabilities are taken as known constants; the
#' more elaborate artifact-distribution machinery of Hunter and Schmidt
#' (2004) is deliberately out of scope here. A corrected correlation that
#' exceeds 1 in magnitude (possible when an observed \eqn{r} outruns the
#' supplied reliabilities) is an error at the pooling stage, unlike the
#' single-study \code{\link{correction_for_attenuation}}, which reports it
#' with a warning: \eqn{\mathrm{atanh}} is undefined there.
#'
#' @return A \code{data.frame} (class \code{dmar_tbl}) with the same
#'   rows as \code{\link{meta_es}}: the \code{estimate},
#'   \code{lower_limit} / \code{upper_limit}, and prediction interval rows
#'   in the correlation metric; the \code{se}, test statistic, and
#'   heterogeneity rows on the Fisher's \emph{Z} scale where the model lives.
#'
#' @references
#' Hunter, J. E., & Schmidt, F. L. (2004). \emph{Methods of meta-analysis:
#'   Correcting error and bias in research findings} (2nd ed.). Sage.
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso \code{\link{meta_es}} for the engine;
#'   \code{\link{correction_for_attenuation}} for the single-study
#'   correction and its connection to latent variable modeling;
#'   \code{\link{reliability}} for estimating the reliabilities;
#'   \code{\link{convert_r_Z}} / \code{\link{convert_Z_r}} for the
#'   transformation used.
#'
#' @family meta-analysis
#'
#' @keywords models
#'
#' @examples
#' # Five validity studies of the same selection instrument.
#' r <- c(.28, .35, .22, .40, .31)
#' n <- c(120, 85, 200, 60, 150)
#' meta_r(r, n)
#'
#' # The same studies corrected for criterion unreliability (reliability
#' # 0.80 in every study): the construct-level validity.
#' meta_r(r, n, reliability_y = 0.80)
#'
#' @export
meta_r <- function(r, n, reliability_x = NULL, reliability_y = NULL,
                   method = c("reml", "pm", "dl", "fe"),
                   hartung_knapp = TRUE, conf_level = 0.95) {
  k <- length(r)
  if (!is.numeric(r) || k < 2L || anyNA(r) || any(abs(r) >= 1)) {
    stop("'r' must be two or more correlations, each in (-1, 1).",
         call. = FALSE)
  }
  if (!is.numeric(n) || length(n) != k || anyNA(n) || any(n < 4) ||
      any(n != round(n))) {
    stop("'n' must give an integer sample size (>= 4) for each study.",
         call. = FALSE)
  }
  expand_rel <- function(rel, nm) {
    if (is.null(rel)) return(rep(1, k))
    if (!is.numeric(rel) || anyNA(rel) || any(rel <= 0) || any(rel > 1) ||
        !(length(rel) %in% c(1L, k))) {
      stop(sprintf("'%s' must be in (0, 1], length 1 or one per study.", nm),
           call. = FALSE)
    }
    rep(rel, length.out = k)
  }
  rxx <- expand_rel(reliability_x, "reliability_x")
  ryy <- expand_rel(reliability_y, "reliability_y")

  r_use <- r / sqrt(rxx * ryy)
  if (any(abs(r_use) >= 1)) {
    stop("After the attenuation correction at least one correlation reaches ",
         "or exceeds 1 in magnitude; the supplied reliabilities cannot ",
         "support these observed correlations.", call. = FALSE)
  }

  yi <- atanh(r_use)
  vi <- 1 / (n - 3)
  fit <- .meta_fit(yi, vi, method = match.arg(method),
                   hartung_knapp = hartung_knapp, conf_level = conf_level)
  .meta_table(fit, conf_level = conf_level, transform = tanh)
}
