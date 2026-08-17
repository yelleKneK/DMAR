#' Correct a Correlation for Attenuation Due to Measurement Error
#'
#' Applies the Spearman (1904) correction for attenuation: the observed
#' correlation between two fallible measures understates the correlation
#' between the constructs they measure, and dividing by the square root of
#' the product of the two reliabilities recovers it,
#' \deqn{r_c \;=\; \frac{r_{XY}}{\sqrt{\rho_{XX'}\,\rho_{YY'}}}.}
#' Within classical test theory (Lord & Novick, 1968), the disattenuated
#' correlation estimates the correlation between the true scores, that is,
#' how strongly the two constructs would correlate if each were measured
#' without error. When \code{N} is supplied, a confidence interval for the
#' corrected correlation is formed by disattenuating the endpoints of the
#' Fisher's \emph{Z} interval for the observed correlation, the standard
#' practice when the reliabilities are treated as known.
#'
#' @param r The observed correlation between the two measures, in
#'   \eqn{[-1, 1]}.
#' @param reliability_x Reliability of the first measure, in \eqn{(0, 1]}.
#'   Any reliability estimate appropriate to the use may be supplied (for
#'   example, coefficient alpha or omega from \code{\link{reliability}}).
#' @param reliability_y Reliability of the second measure, in \eqn{(0, 1]}.
#'   For a correlation between a measure and an error-free criterion,
#'   supply 1 for that side (correcting for criterion unreliability only
#'   yields what the validity generalization literature calls the
#'   operational validity).
#' @param N Optional sample size on which \code{r} is based. When supplied,
#'   a \code{conf_level} confidence interval for the corrected correlation
#'   is reported by correcting the endpoints of the Fisher's \emph{Z}
#'   interval for \code{r}.
#' @param conf_level Confidence level for the interval when \code{N} is
#'   supplied. Defaults to 0.95.
#'
#' @details
#' The correction treats the two reliabilities as known constants, which is
#' the conventional assumption; uncertainty in the reliabilities themselves
#' would widen the interval further. Because the observed correlation can
#' exceed what the supplied reliabilities allow (sampling error, or
#' reliabilities that understate the truth), the corrected value can exceed
#' 1 in magnitude; when that happens the value is reported as computed,
#' with a warning, rather than silently truncated, since a corrected
#' correlation beyond 1 is itself diagnostic information about the inputs.
#'
#' \strong{Prefer the factor model when you have the items.} The Spearman
#' formula is the summary-statistics route: it is exactly right when all
#' you have are the observed correlation and reliability estimates. When
#' the item-level data are available, the better practice is to estimate
#' the construct-level correlation directly as the factor correlation in a
#' two-factor model (each scale loading on its own factor, factors free to
#' correlate): the latent correlation is then estimated jointly with the
#' measurement model rather than assembled from plug-in reliabilities, and
#' it comes with a standard error that propagates the sampling
#' variability of all the moving parts. The example below shows both routes on the same data, the
#' formula route using \code{\link{reliability_omega}}, the model route
#' using \pkg{lavaan}; with congeneric items the two agree closely, and
#' when they disagree the factor model is the one to trust.
#'
#' @return A \code{data.frame} (class \code{dmar_tbl}) in
#'   \code{term} / \code{value} layout with the observed correlation
#'   (\code{correlation_observed}), the corrected correlation
#'   (\code{correlation_corrected}), the corrected interval
#'   (\code{lower_limit}, \code{upper_limit}; present only when \code{N} is
#'   supplied), and the \code{reliability_x}, \code{reliability_y}, and
#'   \code{N} inputs.
#'
#' @references
#' Lord, F. M., & Novick, M. R. (1968). \emph{Statistical theories of
#'   mental test scores}. Addison-Wesley.
#'
#' Spearman, C. (1904). The proof and measurement of association between
#'   two things. \emph{The American Journal of Psychology, 15}(1), 72--101.
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso \code{\link{reliability}} and its family for estimating the
#'   reliabilities supplied here; \code{\link{cfa_1}} and \pkg{lavaan} for
#'   the latent variable route the Details recommend when items are
#'   available; \code{\link{ci_r}} for inference on the
#'   observed correlation itself; \code{\link{convert_r_Z}} and
#'   \code{\link{convert_Z_r}} for the Fisher transformation the interval
#'   uses.
#'
#' @family effect size estimates
#'
#' @keywords design
#'
#' @examples
#' # An observed correlation of .30 between measures with reliabilities .80
#' # and .70 corresponds to a construct-level correlation of about .40.
#' correction_for_attenuation(r = 0.30, reliability_x = 0.80, reliability_y = 0.70)
#'
#' # With the sample size, the corrected interval comes along.
#' correction_for_attenuation(r = 0.30, reliability_x = 0.80, reliability_y = 0.70,
#'                     N = 120)
#'
#' # Correct one side only (error-free criterion).
#' correction_for_attenuation(r = 0.30, reliability_x = 0.80, reliability_y = 1)
#'
#' # The two routes to the construct-level correlation, on the same data
#' # (requires lavaan). Two congeneric scales of three items each whose
#' # latent variables correlate .50:
#' set.seed(113)
#' n <- 400
#' fx <- rnorm(n); fy <- 0.5 * fx + sqrt(1 - 0.25) * rnorm(n)
#' lam <- c(.8, .7, .6)
#' items <- data.frame(
#'   x1 = lam[1] * fx + rnorm(n, 0, sqrt(1 - lam[1]^2)),
#'   x2 = lam[2] * fx + rnorm(n, 0, sqrt(1 - lam[2]^2)),
#'   x3 = lam[3] * fx + rnorm(n, 0, sqrt(1 - lam[3]^2)),
#'   y1 = lam[1] * fy + rnorm(n, 0, sqrt(1 - lam[1]^2)),
#'   y2 = lam[2] * fy + rnorm(n, 0, sqrt(1 - lam[2]^2)),
#'   y3 = lam[3] * fy + rnorm(n, 0, sqrt(1 - lam[3]^2)))
#'
#' # Route 1, summary statistics: omega reliabilities into the formula.
#' x_score <- rowMeans(items[, 1:3]); y_score <- rowMeans(items[, 4:6])
#' om_x <- reliability_omega(data = items[, 1:3])$value[1]
#' om_y <- reliability_omega(data = items[, 4:6])$value[1]
#' correction_for_attenuation(r = cor(x_score, y_score),
#'                            reliability_x = om_x, reliability_y = om_y,
#'                            N = n)
#'
#' # Route 2, the factor model: the latent correlation estimated directly.
#' fit <- lavaan::cfa("X =~ x1 + x2 + x3\nY =~ y1 + y2 + y3",
#'                    data = items, std.lv = TRUE)
#' lavaan::parameterEstimates(fit)[
#'   lavaan::parameterEstimates(fit)$op == "~~" &
#'   lavaan::parameterEstimates(fit)$lhs == "X" &
#'   lavaan::parameterEstimates(fit)$rhs == "Y", ]
#'
#' @export
#' @importFrom stats qnorm
correction_for_attenuation <- function(r, reliability_x, reliability_y,
                                N = NULL, conf_level = 0.95) {
  if (!is.numeric(r) || length(r) != 1L || is.na(r) || abs(r) > 1) {
    stop("'r' must be a single correlation in [-1, 1].", call. = FALSE)
  }
  for (nm in c("reliability_x", "reliability_y")) {
    val <- get(nm)
    if (!is.numeric(val) || length(val) != 1L || is.na(val) ||
        val <= 0 || val > 1) {
      stop(sprintf("'%s' must be a single number in (0, 1].", nm),
           call. = FALSE)
    }
  }
  if (!is.null(N) && (!is.numeric(N) || length(N) != 1L || is.na(N) ||
                      N != round(N) || N < 4)) {
    stop("'N' must be a single integer of at least 4.", call. = FALSE)
  }
  if (!is.numeric(conf_level) || length(conf_level) != 1L ||
      is.na(conf_level) || conf_level <= 0 || conf_level >= 1) {
    stop("'conf_level' must be a single number in (0, 1).", call. = FALSE)
  }

  atten <- sqrt(reliability_x * reliability_y)
  r_c   <- r / atten
  if (abs(r_c) > 1) {
    warning("The corrected correlation exceeds 1 in magnitude (", round(r_c, 3),
            "): the observed correlation is larger than the supplied ",
            "reliabilities allow. Reported as computed, not truncated.",
            call. = FALSE)
  }

  if (is.null(N)) {
    out <- data.frame(
      term  = c("correlation_observed", "correlation_corrected",
                "reliability_x", "reliability_y"),
      value = c(r, r_c, reliability_x, reliability_y),
      stringsAsFactors = FALSE
    )
    return(.as_dmar_tbl(out))
  }

  # Fisher's Z interval for the observed r, endpoints disattenuated. The
  # reliabilities are treated as known, so the correction is a monotone
  # rescaling of the interval for r.
  z    <- atanh(r)
  half <- qnorm(1 - (1 - conf_level) / 2) / sqrt(N - 3)
  lims <- tanh(c(z - half, z + half)) / atten

  out <- data.frame(
    term  = c("correlation_observed", "correlation_corrected",
              "lower_limit", "upper_limit",
              "reliability_x", "reliability_y", "N"),
    value = c(r, r_c, lims[1L], lims[2L], reliability_x, reliability_y, N),
    stringsAsFactors = FALSE
  )
  .as_dmar_tbl(out, conf_level = conf_level)
}
