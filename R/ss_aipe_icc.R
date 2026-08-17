# Sample size for AIPE on an intraclass correlation coefficient.
#' Sample Size for AIPE on an Intraclass Correlation Coefficient
#'
#' Determines the sample size needed for a confidence interval on a
#' population intraclass correlation coefficient (ICC) to have a desired
#' width, using Bonett's (2002) Fisher-style variance-stabilizing
#' transformation. The function inverts the asymptotic variance on the
#' transformed scale (where the CI is symmetric and approximately normal),
#' solves for the smallest \eqn{n} that achieves the target half-width on
#' the back-transformed (raw-ICC) scale, and optionally inflates the
#' result by a chi squared assurance correction (Kelley & Maxwell, 2003).
#'
#' @param rho Anticipated population ICC at the level matching
#'   \code{type}, in \eqn{[0, 1)}. This is the value the researcher
#'   expects the truth to be near, motivated by prior literature, a
#'   pilot, or substantive theory. Because the planning formula
#'   inverts the asymptotic variance \emph{at} this value, a wrong
#'   guess inflates or deflates the realized confidence interval
#'   width relative to \code{width}; the function
#'   \code{\link{ss_aipe_icc_sensitivity}} quantifies the impact of
#'   misspecification by Monte Carlo.
#' @param k Number of raters (or measurements per subject); must be at
#'   least 2.
#' @param width Desired full width of the back-transformed CI on the ICC.
#' @param which_width Whether \code{width} is the \code{"Full"} width of
#'   the interval (default) or a half-width: \code{"Lower"} and
#'   \code{"Upper"} both interpret \code{width} as half the full width,
#'   so they plan for a full width of twice \code{width} and return the
#'   same sample size. Because the interval is not generally symmetric
#'   about the estimate, its realized lower and upper half-widths can
#'   differ from each other and from half the full width; the planner
#'   does not target them separately. A genuinely one-sided width target
#'   is not currently offered.
#' @param conf_level Desired confidence level (default \code{0.95}).
#' @param type Which Shrout-Fleiss (1979) ICC form is being planned. One
#'   of the single-rater forms \code{"ICC(1,1)"}, \code{"ICC(2,1)"},
#'   \code{"ICC(3,1)"} (which share the planning variance) or the
#'   average-of-\eqn{k} forms \code{"ICC(1,k)"}, \code{"ICC(2,k)"},
#'   \code{"ICC(3,k)"}; default \code{"ICC(1,1)"}. Both \code{rho} and
#'   \code{width} are interpreted on the scale of the requested form: for
#'   an average-of-\eqn{k} form the planning value is mapped to the
#'   single-rater scale through the inverse Spearman-Brown relation and
#'   each candidate confidence limit is mapped back, so the returned
#'   sample size targets the width of the interval on the
#'   average-of-\eqn{k} ICC itself (see Details).
#' @param assurance Optional. Probability that the realized CI is no
#'   wider than \code{width}; when supplied, the sample size is inflated
#'   by the standard chi squared correction.
#'
#' @return A \code{data.frame} with rows for the recommended sample
#'   size (\emph{number of subjects}), the expected back-transformed
#'   CI width, and the inputs echoed back. The Shrout-Fleiss form the
#'   plan targets is stored as the \code{"icc_type"} attribute so the
#'   \code{value} column stays numeric.
#'
#' @details
#' \strong{Bonett's (2002) Fisher-style transform.} Bonett (2002) showed
#' that the transformation
#' \deqn{L(\rho) \;=\; \frac{1}{2} \log\!\left(
#'                       \frac{1 + (k - 1)\rho}{1 - \rho}\right)}
#' approximately variance-stabilizes the single-rater ICC, with
#' \deqn{\mathrm{Var}(L(\hat\rho)) \;\approx\;
#'         \frac{k}{2\,(k - 1)\,(n - 2)}.}
#' A confidence interval is constructed by adding \eqn{\pm z_{1-\alpha/2}}
#' standard errors on the \eqn{L} scale and back-transforming to the
#' raw-ICC scale via \eqn{\rho = (e^{2L} - 1) / (e^{2L} - 1 + k)}. The
#' minimum sample size is found by searching for the smallest \eqn{n}
#' whose back-transformed CI width is below the target.
#'
#' \strong{Single-rater vs.\ average-of-\eqn{k} ICC.} The Bonett (2002)
#' variance applies directly to the single-rater forms (\code{ICC(1,1)},
#' \code{ICC(2,1)}, \code{ICC(3,1)}). For the average-of-\eqn{k} forms
#' the planning value \eqn{\rho_k} is first mapped to the single-rater
#' scale through the inverse Spearman-Brown relation
#' \eqn{\rho = \rho_k / [k - (k - 1)\rho_k]}, the interval is formed on
#' the \eqn{L} scale as above, and each candidate limit is mapped back to
#' the average-of-\eqn{k} scale (composing the inverse \eqn{L} transform
#' with the Spearman-Brown formula reduces to
#' \eqn{\rho_k = 1 - e^{-2L}}). The width criterion therefore applies to
#' the confidence interval on the average-of-\eqn{k} ICC itself,
#' following the convention used by \code{\link{var_icc}}. Because the
#' two scales differ, an average-of-\eqn{k} plan generally recommends a
#' different sample size than a single-rater plan at the same numeric
#' \code{rho} and \code{width}; at \code{rho = 0.7}, \code{k = 3}, and
#' \code{width = 0.20}, planning for \code{ICC(1,1)} recommends
#' \eqn{n = 69} subjects while planning for \code{ICC(1,k)} recommends
#' \eqn{n = 110}.
#'
#' \strong{The assurance correction can fall slightly short.} Monte Carlo
#' evaluation with \code{\link{ss_aipe_icc_sensitivity}} shows that the
#' chi squared inflation tends to deliver a little less assurance than
#' requested. At the condition of the second example below
#' (\code{rho = 0.7}, \code{k = 3}, \code{width = 0.20},
#' \code{assurance = 0.80}), the recommended \eqn{n = 79} yields an
#' empirical assurance of about .77 against the requested .80 (10,000
#' replications of the \emph{F}-based interval computed by
#' \code{\link{icc}}), and the smallest sample size whose empirical
#' assurance reaches .80 is \eqn{n = 81}. The mechanism is that the
#' realized interval widths on the raw-ICC scale have a heavier upper
#' tail than the chi squared inflation on the transformed scale accounts
#' for, so the buffer the correction adds is slightly too small. When
#' meeting the assurance target matters, check the recommended sample
#' size with \code{\link{ss_aipe_icc_sensitivity}} and increase \eqn{n}
#' until the empirical assurance reaches the target.
#'
#' @references
#' Bonett, D. G. (2002). Sample size requirements for estimating
#'   intraclass correlations with desired precision. \emph{Statistics in
#'   Medicine, 21}(9), 1331--1335. \doi{10.1002/sim.1108}
#'
#' Donner, A. (1986). A review of inference procedures for the intraclass
#'   correlation coefficient in the one-way random effects model.
#'   \emph{International Statistical Review, 54}(1), 67--82.
#'
#' Kelley, K., & Maxwell, S. E. (2003). Sample size for multiple
#'   regression: Obtaining regression coefficients that are accurate,
#'   not simply significant. \emph{Psychological Methods, 8}(3),
#'   305--321. \doi{10.1037/1082-989X.8.3.305}
#'
#' Shrout, P. E., & Fleiss, J. L. (1979). Intraclass correlations: Uses in
#'   assessing rater reliability. \emph{Psychological Bulletin, 86}(2),
#'   420--428.
#'
#' Smith, C. A. B. (1956). On the estimation of intraclass correlation.
#'   \emph{Annals of Human Genetics, 21}(4), 363--373.
#'
#' @seealso \code{\link{icc}}, \code{\link{var_icc}}
#'
#' @examples
#' # 1. Plan n so the 95% CI on a single-rater ICC has full width <= 0.20
#' #        with k = 3 raters and an anticipated ICC of 0.7.
#' ss_aipe_icc(rho = 0.7, k = 3, width = 0.20)
#'
#' # 2. With 80% assurance:
#' ss_aipe_icc(rho = 0.7, k = 3, width = 0.20, assurance = 0.80)
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @keywords design
#'
#' @seealso \code{\link{design_consequences}} for what a chosen design delivers:
#'   power, the Type S (sign) and Type M (exaggeration) errors of the
#'   significance filter, and the expected confidence interval width.
#'
#' @family AIPE sample size planning
#'
#' @export

ss_aipe_icc <- function(rho, k, width,
                        which_width = c("Full", "Lower", "Upper"),
                        conf_level = 0.95,
                        type = c("ICC(1,1)", "ICC(2,1)", "ICC(3,1)",
                                 "ICC(1,k)", "ICC(2,k)", "ICC(3,k)"),
                        assurance = NULL) {
  which_width <- match.arg(which_width)
  type        <- match.arg(type)
  if (!is.numeric(rho) || length(rho) != 1L || rho < 0 || rho >= 1)
    stop("'rho' must be a single value in [0, 1).")
  if (!is.numeric(k) || length(k) != 1L || k < 2)
    stop("'k' (raters per subject) must be a single integer >= 2.")
  if (!is.numeric(width) || length(width) != 1L || width <= 0)
    stop("'width' must be a single positive number.")
  if (conf_level <= 0 || conf_level >= 1)
    stop("'conf_level' must be in (0, 1).")

  average_k  <- type %in% c("ICC(1,k)", "ICC(2,k)", "ICC(3,k)")
  half_width <- if (which_width == "Full") width / 2 else width
  z_alpha    <- stats::qnorm(1 - (1 - conf_level) / 2)

  # 'rho' and 'width' arrive on the scale matching 'type'. Planning happens
  # on the single-rater scale, where the Bonett (2002) variance applies; an
  # average-of-k planning value is first mapped down through the inverse
  # Spearman-Brown relation (the convention var_icc() uses).
  rho_single <- if (average_k) rho / (k - (k - 1) * rho) else rho

  # Bonett (2002) L-transform of the single-rater rho:
  L_rho <- 0.5 * log((1 + (k - 1) * rho_single) / (1 - rho_single))

  # Search for smallest n that gives back-transformed half-width <= target.
  # The back-transform lands on the scale the target width lives on: the
  # single-rater scale directly, or the average-of-k scale by composing the
  # inverse L-transform with Spearman-Brown, which reduces to
  # rho_k = 1 - exp(-2L).
  back_transform <- if (average_k) {
    function(L) 1 - exp(-2 * L)
  } else {
    function(L) (exp(2 * L) - 1) / (exp(2 * L) - 1 + k)
  }

  n_seq <- 4:5000
  half_back <- vapply(n_seq, function(nn) {
    se_L <- sqrt(k / (2 * (k - 1) * (nn - 2)))
    rho_hi <- back_transform(L_rho + z_alpha * se_L)
    rho_lo <- back_transform(L_rho - z_alpha * se_L)
    (rho_hi - rho_lo) / 2
  }, numeric(1))

  hit <- which(half_back <= half_width)
  if (length(hit) == 0L) {
    stop("No sample size up to 5000 achieves the requested precision; ",
         "the target width may be unattainable for the given rho and k.")
  }
  n <- n_seq[hit[1]]

  if (!is.null(assurance)) {
    if (assurance <= 0.5 || assurance >= 1)
      stop("'assurance' must be in (0.5, 1).")
    df <- max(1, n - 1)
    inflate <- stats::qchisq(assurance, df = df) / df
    n <- ceiling(n * inflate)
  }

  # Final expected width at the recommended n.
  se_L <- sqrt(k / (2 * (k - 1) * (n - 2)))
  rho_hi <- back_transform(L_rho + z_alpha * se_L)
  rho_lo <- back_transform(L_rho - z_alpha * se_L)
  expected_w <- rho_hi - rho_lo

  out <- data.frame(
    term  = c("necessary_N", "expected_width", "rho", "k",
              "width_target", "conf_level"),
    value = c(n, expected_w, rho, k, width, conf_level),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  attr(out, "icc_type") <- type
  .as_dmar_tbl(out, conf_level = conf_level, subclass = "dmar_ss_aipe")
}
