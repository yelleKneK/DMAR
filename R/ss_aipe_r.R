# Sample size for AIPE on a Pearson correlation.
#' Sample Size for AIPE on a Pearson Correlation
#'
#' Determines the sample size needed for a confidence interval on a
#' population Pearson correlation \eqn{\rho} to have a desired width
#' (accuracy in parameter estimation; Kelley & Maxwell, 2003). The
#' interval planned for is the Fisher's \emph{Z} interval that
#' \code{\link{correlations_test}} reports for
#' \code{method = "pearson"}: the correlation is transformed as
#' \eqn{z(r) = \mathrm{atanh}(r)}, an interval with standard error
#' \eqn{1/\sqrt{n - 3}} is formed on the \emph{z} scale, and the limits
#' are back-transformed through \eqn{\tanh(\cdot)} (Fisher, 1921;
#' Bonett & Wright, 2000, Equation 2). Because the plan targets the same
#' interval the analysis will report, the planned width and the analyzed
#' width agree.
#'
#' @param rho Anticipated population Pearson correlation, in
#'   \eqn{(-1, 1)}.
#' @param width Desired full width of the confidence interval on the
#'   correlation.
#' @param conf_level Desired confidence level (default \code{0.95}).
#' @param assurance Optional. Probability that the realized CI is no
#'   wider than \code{width} (\eqn{1 - \gamma}). When supplied, the
#'   sample size is inflated using the standard chi squared correction
#'   (Kelley, 2008); when \code{NULL}, the assurance is fixed at 0.5.
#'
#' @return A \code{data.frame} with the rows \code{necessary_N}
#'   (the recommended total sample size, rounded up),
#'   \code{expected_width} at that sample size, and the inputs echoed
#'   back.
#'
#' @details
#' \strong{Closed-form first pass.} On the Fisher's \emph{Z} scale the
#' interval has half-width \eqn{z_{1 - \alpha/2} / \sqrt{n - 3}}, and the
#' delta method maps it back to the correlation scale as approximately
#' \deqn{w \;\approx\; 2\, z_{1 - \alpha/2} \,
#'   \frac{1 - \rho^2}{\sqrt{n - 3}},}
#' which solves to the first-stage approximation of Bonett and Wright
#' (2000),
#' \deqn{n_0 \;=\; 3 + \Big\lceil 4\, (z_{1 - \alpha/2})^2
#'   (1 - \rho^2)^2 / w^2 \Big\rceil.}
#'
#' \strong{Exact iteration.} The back-transformed width depends on
#' \eqn{\rho} through \eqn{\tanh(\cdot)}, so the delta method
#' approximation can land a few observations off in either direction.
#' Starting from \eqn{n_0}, the function evaluates the exact
#' back-transformed width
#' \eqn{\tanh(z_\rho + z_{1 - \alpha/2}/\sqrt{n - 3}) -
#'   \tanh(z_\rho - z_{1 - \alpha/2}/\sqrt{n - 3})}
#' and steps the integer \eqn{n} until it is the smallest sample size
#' whose width is at or below \code{width}. Where Bonett and Wright
#' (2000) stop after a single second-stage adjustment, this search is
#' exact.
#'
#' \strong{The planning value matters least near zero.} At a fixed
#' sample size the back-transformed width is largest at \eqn{\rho = 0}
#' and shrinks as \eqn{|\rho|} grows, so a planning value closer to zero
#' yields a larger, more conservative sample size. When little is known
#' about the population correlation, \code{rho = 0} gives the sample
#' size that suffices for any population value.
#'
#' \strong{When to use simple vs. partial correlation planning.} Use
#' this function when the inferential target is the correlation between
#' two variables with nothing partialed out. When the target is the
#' correlation after statistically controlling for other variables, see
#' \code{\link{ss_aipe_partial_r}}.
#'
#' The Monte Carlo companion \code{\link{ss_aipe_r_sensitivity}}
#' evaluates how the plan behaves when the population correlation
#' differs from the planning value.
#'
#' @references
#' Bonett, D. G., & Wright, T. A. (2000). Sample size requirements for
#'   estimating Pearson, Kendall and Spearman correlations.
#'   \emph{Psychometrika, 65}(1), 23--28. \doi{10.1007/BF02294183}
#'
#' Fisher, R. A. (1921). On the "probable error" of a coefficient of
#'   correlation deduced from a small sample. \emph{Metron, 1}, 3--32.
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
#' Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). \emph{Designing
#'   experiments and analyzing data: A model comparison perspective}
#'   (4th ed.). Routledge. (See Chapter 9 on correlations.)
#'
#' @seealso \code{\link{ss_aipe_r_sensitivity}},
#'   \code{\link{correlations_test}}, \code{\link{ss_aipe_partial_r}},
#'   \code{\link{ss_power_r}}, \code{\link{convert_r_Z}}
#'
#' @examples
#' # Plan n so the 95% CI on the Pearson correlation has full width
#' # at most 0.20, when the anticipated correlation is 0.30.
#' ss_aipe_r(rho = 0.30, width = 0.20)
#'
#' # A narrower target width requires a larger sample size.
#' ss_aipe_r(rho = 0.30, width = 0.10)
#'
#' # With 80% assurance that the realized interval is no wider than
#' # the target (Kelley, 2008):
#' ss_aipe_r(rho = 0.30, width = 0.20, assurance = 0.80)
#'
#' # Planning at rho = 0 gives the sample size that suffices for any
#' # population correlation, since the interval is widest there.
#' ss_aipe_r(rho = 0, width = 0.20)
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

ss_aipe_r <- function(rho, width,
                      conf_level = 0.95,
                      assurance = NULL) {
  if (!is.numeric(rho) || length(rho) != 1L || abs(rho) >= 1)
    stop("'rho' must be a single value in (-1, 1).")
  if (!is.numeric(width) || length(width) != 1L || !is.finite(width) ||
      width <= 0)
    stop("'width' must be a single positive number.")
  if (!is.numeric(conf_level) || length(conf_level) != 1L ||
      conf_level <= 0 || conf_level >= 1)
    stop("'conf_level' must be in (0, 1).")

  z_alpha <- stats::qnorm(1 - (1 - conf_level) / 2)
  z_rho   <- atanh(rho)

  # Full width of the back-transformed Fisher's Z interval at sample size nn.
  # This is the interval correlations_test() reports for method = "pearson":
  # SE(z) = 1 / sqrt(n - 3), limits back-transformed through tanh.
  full_back <- function(nn) {
    se_z <- 1 / sqrt(nn - 3)
    tanh(z_rho + z_alpha * se_z) - tanh(z_rho - z_alpha * se_z)
  }

  n_min <- 4    # smallest n with a positive Fisher's Z variance (n - 3 > 0)

  # Closed-form first pass: the delta method width 2 z (1 - rho^2) / sqrt(n - 3)
  # solved for n, the first-stage approximation of Bonett and Wright (2000).
  # The candidate stays a double (not an integer) so that a very narrow
  # width target cannot overflow R's integer range; whole-number doubles
  # count exactly.
  n0 <- 3 + ceiling(4 * z_alpha^2 * (1 - rho^2)^2 / width^2)
  n  <- max(n_min, n0)

  # Exact refinement: the back-transformed width is strictly decreasing in n,
  # so walk the integer n in whichever direction the closed form missed and
  # stop at the smallest n whose exact width is at or below the target.
  while (full_back(n) > width) n <- n + 1
  while (n > n_min && full_back(n - 1) <= width) n <- n - 1

  # Assurance correction (Kelley, 2008): conservative chi squared inflation.
  if (!is.null(assurance)) {
    if (!is.numeric(assurance) || length(assurance) != 1L ||
        assurance <= 0.5 || assurance >= 1)
      stop("'assurance' must be in (0.5, 1).")
    df <- max(1, n - 2)
    inflate <- stats::qchisq(assurance, df = df) / df
    n <- ceiling(n * inflate)
  }

  out <- data.frame(
    term  = c("necessary_N", "expected_width", "rho",
              "width_target", "conf_level"),
    value = c(n, full_back(n), rho, width, conf_level),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  .as_dmar_tbl(out, conf_level = conf_level, subclass = "dmar_ss_aipe")
}
