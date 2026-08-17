# Sample size for AIPE on a partial correlation.
#' Sample Size for AIPE on a Partial Correlation
#'
#' Determines the sample size needed for a confidence interval on a
#' population partial correlation \eqn{\rho_{XY \cdot Z_1 \cdots Z_J}} to
#' have a desired width (Accuracy in Parameter Estimation; Kelley, 2008).
#' The function inverts the asymptotic variance of the partial Pearson
#' correlation, either on the raw scale (Olkin & Finn, 1995) or on the
#' Fisher's \emph{Z}-transformed scale (Fisher, 1921, 1924; Bonett, 2008), and
#' solves for the smallest \eqn{n} that achieves the target half-width
#' (or full width).
#'
#' @param rho Anticipated population partial correlation, in \eqn{(-1, 1)}.
#' @param J Number of variables partialled out (count of
#'   \eqn{Z_1, \ldots, Z_J}); must be at least 1.
#' @param width Desired full width of the confidence interval on the
#'   partial correlation.
#' @param which_width Whether \code{width} is the \code{"Full"} width of
#'   the interval (default) or a half-width: \code{"Lower"} and
#'   \code{"Upper"} both interpret \code{width} as half the full width,
#'   so they plan for a full width of twice \code{width} and return the
#'   same sample size. Because the interval is not generally symmetric
#'   about the estimate (with \code{fisher_z = TRUE} in particular), its
#'   realized lower and upper half-widths can differ from each other and
#'   from half the full width; the planner does not target them
#'   separately. A genuinely one-sided width target is not currently
#'   offered.
#' @param conf_level Desired confidence level (default \code{0.95}).
#' @param fisher_z Logical. If \code{TRUE}, the half-width target is
#'   applied on the Fisher-\emph{z} scale (variance \eqn{1/(n - J - 3)};
#'   Bonett, 2008) and the resulting CI is back-transformed via
#'   \eqn{\tanh}. If \code{FALSE} (default), uses the raw-scale
#'   Olkin-Finn (1995) asymptotic variance.
#' @param assurance Optional. Probability that the realized CI is no
#'   wider than \code{width} (\eqn{1 - \gamma}). When supplied, the
#'   sample size is inflated using the standard chi squared correction
#'   (Kelley, 2008); when \code{NULL}, the assurance is fixed at 0.5.
#'
#' @return A \code{data.frame} with the rows \code{necessary_N}
#'   (the recommended total sample size, rounded up), \code{expected_width}
#'   at that sample size, and the inputs echoed back.
#'
#' @details
#' \strong{Raw-scale Olkin-Finn asymptotic variance.} The half-width of a
#' \eqn{100(1 - \alpha)\%} CI on the partial Pearson correlation is
#' approximately
#' \deqn{w_{1/2} \;\approx\; z_{1 - \alpha/2} \cdot
#'   \sqrt{\,\frac{(1 - \rho_{XY \cdot Z}^{\,2})^2}{n - J - 1}\,}.}
#' Solving for \eqn{n}:
#' \deqn{n \;=\; J + 1 + \Big\lceil
#'         (z_{1 - \alpha/2})^2 \cdot (1 - \rho_{XY \cdot Z}^{\,2})^2
#'         / w_{1/2}^{2}
#'       \Big\rceil.}
#' This is the planning analog of the half-width of \code{\link{ci_r}}
#' applied to a partial correlation.
#'
#' \strong{Fisher-\emph{z} scale (recommended for small \eqn{\rho}, near
#' boundary, or small \eqn{n - J}).} Bonett (2008) advocates planning on
#' the variance-stabilized Fisher-\emph{z} scale and back-transforming
#' the bounds. On the Fisher's \eqn{Z} scale, the asymptotic half-width is
#' \deqn{w^{(z)}_{1/2} \;\approx\; z_{1 - \alpha/2} / \sqrt{n - J - 3}.}
#' Solving for the \eqn{n} that achieves a given back-transformed
#' \eqn{w_{1/2}} is done by a 1-D search; this is generally the more
#' accurate route when \eqn{n} is small or \eqn{|\rho|} is large.
#'
#' \strong{When to use partial vs. simple correlation planning.} Use this
#' function when the inferential target is the population correlation
#' between \eqn{X} and \eqn{Y} \emph{after} statistically controlling for
#' \eqn{Z_1, \ldots, Z_J}. For the simple Pearson correlation, see
#' \code{\link{ss_aipe_r}}.
#'
#' \strong{Note on conservatism of the assurance plan.} The empirical
#' simulation study of the AIPE planner family finds that
#' \code{ss_aipe_partial_r()} is tight
#' (zero overshoot) at 80\% assurance but modestly conservative at 99\%
#' assurance, with an empirical ideal sample size of about 5 to 10
#' subjects smaller than the recommended sample size. The mechanism is
#' the usual one for AIPE assurance plans: the Olkin-Finn (1995) Wald-
#' style upper bound on \eqn{\Pr(\widehat W > \omega)} that the planner
#' inverts is not tight at the recommended sample size, especially at
#' the 99\% level where the inversion has to push further into the
#' upper tail of \eqn{\widehat W}. The recommended sample size is a
#' sufficient sample size rather than the smallest possible sample
#' size. \code{\link{ss_aipe_partial_r_sensitivity}} quantifies the
#' overshoot for any one condition.
#'
#' @references
#' Algina, J., & Olejnik, S. (2003). Sample size tables for correlation
#'   analysis with applications in partial correlation and multiple
#'   regression analysis. \emph{Multivariate Behavioral Research, 38}(3),
#'   309--323. \doi{10.1207/s15327906mbr3803_02}
#'
#' Bonett, D. G. (2008). Confidence intervals for standardized linear
#'   contrasts of means. \emph{Psychological Methods, 13}(2), 99--109.
#'   \doi{10.1037/1082-989X.13.2.99}
#'
#' Fisher, R. A. (1921). On the "probable error" of a coefficient of
#'   correlation deduced from a small sample. \emph{Metron, 1}, 3--32.
#'
#' Fisher, R. A. (1924). The distribution of the partial correlation
#'   coefficient. \emph{Metron, 3}, 329--332.
#'
#' Kelley, K. (2008). Sample size planning for the squared multiple
#'   correlation coefficient: Accuracy in parameter estimation via narrow
#'   confidence intervals. \emph{Multivariate Behavioral Research, 43}(4),
#'   524--555. \doi{10.1080/00273170802490632}
#'
#' Kelley, K., & Maxwell, S. E. (2003). Sample size for multiple
#'   regression: Obtaining regression coefficients that are accurate,
#'   not simply significant. \emph{Psychological Methods, 8}(3),
#'   305--321. \doi{10.1037/1082-989X.8.3.305}
#'
#' Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027).
#'   \emph{Designing experiments and analyzing data: A model comparison
#'   perspective} (4th ed.). Routledge. (See Chapter 3 on the one-way
#'   ANOVA and Chapter 4 on contrasts.)
#'
#' Olkin, I., & Finn, J. D. (1995). Correlations redux.
#'   \emph{Psychological Bulletin, 118}(1), 155--164.
#'   \doi{10.1037/0033-2909.118.1.155}
#'
#' @seealso \code{\link{var_partial_r}}, \code{\link{expected_partial_r}},
#'   \code{\link{ss_aipe_semipartial_r}}, \code{\link{ss_aipe_r}}
#'
#' @examples
#' # 1. Plan n so the 95% CI on rho_XY.Z (J = 2 controls) has
#' #        full width <= 0.20, when the anticipated partial r is 0.30.
#' ss_aipe_partial_r(rho = 0.30, J = 2, width = 0.20)
#'
#' # 2. Same problem on the Fisher's Z scale (Bonett 2008):
#' ss_aipe_partial_r(rho = 0.30, J = 2, width = 0.20, fisher_z = TRUE)
#'
#' # 3. With 80% assurance (Kelley 2008):
#' ss_aipe_partial_r(rho = 0.30, J = 2, width = 0.20, assurance = 0.80)
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

ss_aipe_partial_r <- function(rho, J, width,
                              which_width = c("Full", "Lower", "Upper"),
                              conf_level = 0.95,
                              fisher_z = FALSE,
                              assurance = NULL) {
  which_width <- match.arg(which_width)
  if (!is.numeric(rho) || length(rho) != 1L || abs(rho) >= 1)
    stop("'rho' must be a single value in (-1, 1).")
  if (!is.numeric(J) || length(J) != 1L || J < 1)
    stop("'J' must be a single integer >= 1.")
  if (!is.numeric(width) || length(width) != 1L || width <= 0)
    stop("'width' must be a single positive number.")
  if (conf_level <= 0 || conf_level >= 1)
    stop("'conf_level' must be in (0, 1).")

  half_width <- if (which_width == "Full") width / 2 else width
  z_alpha    <- stats::qnorm(1 - (1 - conf_level) / 2)

  if (!fisher_z) {
    # Raw-scale Olkin-Finn asymptotic variance: n - J - 1 = z^2 (1 - rho^2)^2 / w^2
    n <- J + 1 + ceiling((z_alpha^2 * (1 - rho^2)^2) / (half_width^2))
  } else {
    # Fisher's Z scale: the smallest n whose back-transformed interval half-width
    # is <= the target on the raw scale. The half-width decreases monotonically
    # in n toward 0, so a finite n always exists; solve for it directly. The old
    # fixed grid (J + 5 to J + 5000) silently returned its minimum, J + 5,
    # whenever the target fell past the grid's right end, since which.max() of an
    # all-FALSE vector is 1.
    z_rho <- atanh(rho)
    half_back <- function(nn) {
      se_z <- 1 / sqrt(nn - J - 3)
      (tanh(z_rho + z_alpha * se_z) - tanh(z_rho - z_alpha * se_z)) / 2
    }
    n_min <- J + 4                         # smallest n with a positive Fisher's Z df
    if (half_back(n_min) <= half_width) {
      n <- n_min
    } else {
      n_hi <- n_min + 1
      while (half_back(n_hi) > half_width) n_hi <- n_min + 2 * (n_hi - n_min)
      root <- stats::uniroot(function(nn) half_back(nn) - half_width,
                             lower = n_min, upper = n_hi)$root
      n <- max(n_min, ceiling(root))
      if (half_back(n) > half_width) n <- n + 1L   # guard against rounding down
    }
  }

  # Assurance correction (Kelley, 2008): conservative chi squared inflation.
  if (!is.null(assurance)) {
    if (assurance <= 0.5 || assurance >= 1)
      stop("'assurance' must be in (0.5, 1).")
    df <- max(1, n - J - 2)
    inflate <- stats::qchisq(assurance, df = df) / df
    n <- ceiling(n * inflate)
  }

  # Expected width at the recommended n. On the Fisher's Z path this is the
  # back-transformed full width; otherwise the raw-scale asymptotic width.
  expected_w <- if (fisher_z) {
    2 * half_back(n)
  } else {
    2 * z_alpha * sqrt((1 - rho^2)^2 / (n - J - 1))
  }

  out <- data.frame(
    term  = c("necessary_N", "expected_width", "rho", "J",
              "width_target", "conf_level"),
    value = c(n, expected_w, rho, J, width, conf_level),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  .as_dmar_tbl(out, conf_level = conf_level, subclass = "dmar_ss_aipe")
}
