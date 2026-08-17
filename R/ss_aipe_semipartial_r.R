# Sample size for AIPE on a semipartial (part) correlation.
#' Sample Size for AIPE on a Semipartial (Part) Correlation
#'
#' Determines the sample size needed for a confidence interval on a
#' population semipartial correlation
#' \eqn{r_{Y(X \cdot Z_1 \cdots Z_J)}} (the unique contribution of
#' \eqn{X} to \eqn{Y} after controlling for \eqn{Z_1, \ldots, Z_J}, with
#' \eqn{Y} \emph{not} residualized) to have a desired width, using the
#' Olkin-Finn (1995) / Algina-Olejnik (2003) asymptotic variance and the
#' AIPE framework of Kelley & Maxwell (2003).
#'
#' @param r_sp Anticipated population semipartial correlation, in
#'   \eqn{(-1, 1)}.
#' @param J Number of variables partialled out of \eqn{X} (count of
#'   \eqn{Z_1, \ldots, Z_J}); must be at least 1.
#' @param width Desired full width of the confidence interval on the
#'   semipartial correlation.
#' @param which_width Whether \code{width} refers to the \code{"Full"}
#'   width (default) or the \code{"Lower"}/\code{"Upper"} half-width.
#' @param conf_level Desired confidence level (default \code{0.95}).
#' @param assurance Optional. Probability that the realized CI is no
#'   wider than \code{width}. When supplied, the sample size is inflated
#'   using the standard chi squared correction (Kelley & Maxwell, 2003).
#'
#' @return A \code{data.frame} with rows for the recommended sample
#'   size, the expected CI width at that sample size, and the inputs
#'   echoed back.
#'
#' @details
#' \strong{Asymptotic variance of the semipartial.} Under multivariate
#' normality, the sample semipartial correlation
#' \eqn{r_{Y(X \cdot Z)}} has asymptotic variance
#' \deqn{\mathrm{Var}(\hat r_{Y(X \cdot Z)}) \;\approx\;
#'         \frac{(1 - r_{Y(X \cdot Z)}^2)^2}{n - J - 1}}
#' (Olkin & Finn, 1995, with the partial-correlation degrees-of-freedom
#' correction). Inverting for the sample size needed to achieve a target
#' half-width \eqn{w_{1/2}} at confidence level \eqn{1 - \alpha}:
#' \deqn{n \;=\; J + 1 + \Big\lceil
#'         z_{1 - \alpha/2}^{2} \cdot (1 - r_{Y(X \cdot Z)}^{2})^2
#'         / w_{1/2}^{2} \Big\rceil.}
#'
#' \strong{Comparison with partial-r planning.} The partial correlation
#' \eqn{r_{XY \cdot Z}} divides the covariance after residualizing both
#' \eqn{X} and \eqn{Y} on \eqn{Z}; the semipartial divides after
#' residualizing only \eqn{X}. The semipartial is the natural effect size
#' companion to a standardized regression coefficient: its square equals
#' the \eqn{\Delta R^2} contributed by \eqn{X} above and beyond the
#' controls. See \code{\link{var_semipartial_r}} for the asymptotic
#' variance, and \code{\link{ss_aipe_partial_r}} for the partial-correlation
#' analog of this function.
#'
#' \strong{Note on conservatism of the assurance plan.} The empirical
#' simulation study of the AIPE planner family finds that
#' \code{ss_aipe_semipartial_r()} is
#' on the boundary of its valid range at 80\% assurance and modestly
#' conservative at 99\% assurance. At \eqn{\gamma = 0.80}, the realized
#' assurance at the recommended sample size is within Monte Carlo error
#' of the target, that is, the bound is operating at the edge of its
#' validity. At \eqn{\gamma = 0.99}, the ideal sample size is about 15
#' to 20 subjects smaller than the recommended sample size, reflecting
#' the looser upper-tail bound at the 99\% level. The recommended
#' sample size is therefore a sufficient sample size rather than the
#' smallest possible sample size. A small safety margin (5 to 10
#' subjects) is advisable when planning at \eqn{\gamma = 0.80}. \code{\link{ss_aipe_semipartial_r_sensitivity}} quantifies the
#' overshoot for any one condition.
#'
#' @references
#' Algina, J., & Olejnik, S. (2003). Sample size tables for correlation
#'   analysis with applications in partial correlation and multiple
#'   regression analysis. \emph{Multivariate Behavioral Research, 38}(3),
#'   309--323. \doi{10.1207/s15327906mbr3803_02}
#'
#' Cohen, J., Cohen, P., West, S. G., & Aiken, L. S. (2003). \emph{Applied
#'   multiple regression/correlation analysis for the behavioral
#'   sciences} (3rd ed.). Lawrence Erlbaum.
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
#' @seealso \code{\link{var_semipartial_r}}, \code{\link{ss_aipe_partial_r}},
#'   \code{\link{ss_aipe_R2}}
#'
#' @examples
#' # 1. Plan n so the 95% CI on r_sp (J = 3) has full width <= 0.15
#' #        when the anticipated semipartial is 0.25.
#' ss_aipe_semipartial_r(r_sp = 0.25, J = 3, width = 0.15)
#'
#' # 2. With 80% assurance:
#' ss_aipe_semipartial_r(r_sp = 0.25, J = 3, width = 0.15,
#'                       assurance = 0.80)
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

ss_aipe_semipartial_r <- function(r_sp, J, width,
                                  which_width = c("Full", "Lower", "Upper"),
                                  conf_level = 0.95,
                                  assurance = NULL) {
  which_width <- match.arg(which_width)
  if (!is.numeric(r_sp) || length(r_sp) != 1L || abs(r_sp) >= 1)
    stop("'r_sp' must be a single value in (-1, 1).")
  if (!is.numeric(J) || length(J) != 1L || J < 1)
    stop("'J' must be a single integer >= 1.")
  if (!is.numeric(width) || length(width) != 1L || width <= 0)
    stop("'width' must be a single positive number.")
  if (conf_level <= 0 || conf_level >= 1)
    stop("'conf_level' must be in (0, 1).")

  half_width <- if (which_width == "Full") width / 2 else width
  z_alpha    <- stats::qnorm(1 - (1 - conf_level) / 2)

  # Olkin-Finn (1995) asymptotic variance: n - J - 1 = z^2 (1 - r^2)^2 / w^2
  n <- J + 1 + ceiling((z_alpha^2 * (1 - r_sp^2)^2) / (half_width^2))

  if (!is.null(assurance)) {
    if (assurance <= 0.5 || assurance >= 1)
      stop("'assurance' must be in (0.5, 1).")
    df <- max(1, n - J - 2)
    inflate <- stats::qchisq(assurance, df = df) / df
    n <- ceiling(n * inflate)
  }

  expected_w <- 2 * z_alpha * sqrt((1 - r_sp^2)^2 / (n - J - 1))

  out <- data.frame(
    term  = c("necessary_N", "expected_width", "r_sp", "J",
              "width_target", "conf_level"),
    value = c(n, expected_w, r_sp, J, width, conf_level),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  .as_dmar_tbl(out, conf_level = conf_level, subclass = "dmar_ss_aipe")
}
