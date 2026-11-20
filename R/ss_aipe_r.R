# Sample size for AIPE on a Pearson correlation.
#' Sample Size for AIPE on a Pearson Correlation
#'
#' Determines the sample size needed for a confidence interval on a
#' population Pearson correlation \eqn{\rho} to have a desired width
#' (accuracy in parameter estimation; Kelley & Maxwell, 2003). The
#' interval planned for is the Fisher's \emph{Z} interval that
#' \code{\link{correlations_test}} reports for
#' \code{method = "pearson"}: the correlation is transformed as
#' \eqn{Z(r) = \mathrm{atanh}(r)}, an interval with standard error
#' \eqn{1/\sqrt{N - 3}} is formed on the \emph{Z} scale, and the limits
#' are back-transformed through \eqn{\tanh(\cdot)} (Fisher, 1921;
#' Bonett & Wright, 2000, Equation 2). Because the plan targets the same
#' interval the analysis will report, the planned width and the analyzed
#' width agree.
#'
#' @param true_cor Population Pearson correlation \eqn{\rho} the sample size is
#'   planned for, in \eqn{(-1, 1)}. The planner treats this value as the
#'   truth; \code{\link{ss_aipe_r_sensitivity}} shows what the plan delivers
#'   when the truth (\code{true_cor} there) differs from the value planned
#'   with (\code{supposed_cor} there).
#' @param width Desired full width of the confidence interval on the
#'   correlation.
#' @param conf_level Desired confidence level (default \code{0.95}). A level above .999999998, whose tail areas are below 1e-9, is refused (see "Small tail areas" in \code{\link{ci_nc_t}}).
#' @param assurance Optional. The probability that the realized interval
#'   is no wider than \code{width}. The function returns the smallest
#'   sample size, at or above the sample size planned without it, at which
#'   that probability is at least \code{assurance}. A single number above
#'   0.5 and below 1; when \code{NULL} (the default) the plan is for the
#'   width at \eqn{r = \rho} only.
#'
#' @return A \code{data.frame} with the columns \code{term} and
#'   \code{value} and one row per term:
#'   \describe{
#'     \item{\code{necessary_N}}{the planned total sample size.}
#'     \item{\code{expected_width}}{the full width of the interval at
#'       \code{necessary_N} when the sample correlation equals
#'       \code{true_cor}.}
#'     \item{\code{true_cor}, \code{width}, \code{conf_level}}{the
#'       inputs.}
#'     \item{\code{assurance}}{the requested probability; present only
#'       when \code{assurance} was given.}
#'     \item{\code{assurance_attained}}{the probability, at
#'       \code{necessary_N} and \code{true_cor}, that the realized interval
#'       is no wider than \code{width}; present only when \code{assurance}
#'       was given. It is at least \code{assurance}, and above it because
#'       the sample size is a whole number or because the plan without
#'       \code{assurance} already exceeds the request.}
#'   }
#'
#' @details
#' \strong{Closed-form first pass.} On the Fisher's \emph{Z} scale the
#' interval has half-width \eqn{z_{1 - \alpha/2} / \sqrt{N - 3}}, and the
#' delta method maps it back to the correlation scale as approximately
#' \deqn{w \;\approx\; 2\, z_{1 - \alpha/2} \,
#'   \frac{1 - \rho^2}{\sqrt{N - 3}},}
#' which solves to the first-stage approximation of Bonett and Wright
#' (2000),
#' \deqn{n_0 \;=\; 3 + \Big\lceil 4\, (z_{1 - \alpha/2})^2
#'   (1 - \rho^2)^2 / w^2 \Big\rceil.}
#' Here and below \eqn{\rho} is \code{true_cor}, \eqn{w} is \code{width},
#' and \eqn{1 - \alpha} is \code{conf_level}.
#'
#' \strong{Exact iteration.} The back-transformed width depends on
#' \eqn{\rho} through \eqn{\tanh(\cdot)}, so the delta method
#' approximation can land a few observations off in either direction.
#' Starting from \eqn{n_0}, the function evaluates the exact
#' back-transformed width
#' \eqn{\tanh(z_\rho + z_{1 - \alpha/2}/\sqrt{N - 3}) -
#'   \tanh(z_\rho - z_{1 - \alpha/2}/\sqrt{N - 3})}
#' and steps the integer \eqn{N} until it is the smallest sample size
#' whose width is at or below \code{width}. Where Bonett and Wright
#' (2000) stop after a single second-stage adjustment, this search is
#' exact.
#'
#' \strong{The assurance step.} A plan without \code{assurance} makes the
#' interval as narrow as the target when the sample correlation equals
#' \code{true_cor}; in many samples it is wider. With \code{assurance}
#' the function finds the sample size at which the realized interval is
#' no wider than the target with at least that probability, from the
#' sampling distribution of the sample correlation \emph{r} of \emph{N}
#' independent observations from a bivariate normal population:
#' \enumerate{
#'   \item At a given \emph{N} the width of the interval depends on the
#'     sample only through \emph{r}, and it decreases as \eqn{|r|} grows.
#'     The realized interval is therefore no wider than the target exactly
#'     when \eqn{|r| \ge r^*(N)}, where \eqn{r^*(N)} is the root of
#'     \eqn{\tanh(\mathrm{atanh}(r) + z_{1 - \alpha/2}/\sqrt{N - 3}) -
#'     \tanh(\mathrm{atanh}(r) - z_{1 - \alpha/2}/\sqrt{N - 3}) = w}.
#'     When the width at \eqn{r = 0}, the widest the interval can be,
#'     already meets the target, every sample does and the probability is
#'     1.
#'   \item The distribution function of \emph{r} is computed as
#'     \deqn{F(c) \;=\; \int_0^1 T_{2N - 3}\!\Big(\sqrt{2N - 3}\,
#'       \Big\{\frac{c\,\sqrt{1 - b}}{\sqrt{1 - c^2}} -
#'       \frac{\rho\,\sqrt{b}}{\sqrt{1 - \rho^2}}\Big\}\Big)\,
#'       g(b)\, db,}
#'     where \eqn{T_{2N - 3}} is the distribution function of Student's
#'     \emph{t} on \eqn{2N - 3} degrees of freedom and \eqn{g} is the
#'     density of the beta distribution with parameters \eqn{(N - 1)/2}
#'     and \eqn{(N - 2)/2}. The form follows from writing
#'     \eqn{r / \sqrt{1 - r^2}} as \eqn{(Z + \theta \sqrt{U}) / \sqrt{V}},
#'     with \eqn{Z} standard normal, \eqn{U} and \eqn{V} chi square on
#'     \eqn{N - 1} and \eqn{N - 2} degrees of freedom, all independent,
#'     and \eqn{\theta = \rho / \sqrt{1 - \rho^2}}, and then conditioning
#'     on \eqn{U / (U + V)}. It is the computation
#'     \code{\link{ss_aipe_partial_r}} uses with no variables partialed
#'     out.
#'   \item The probability at \emph{N} is
#'     \eqn{P_N = 1 - F(r^*(N)) + F(-r^*(N))}, and the function returns
#'     the smallest \emph{N}, at or above the sample size planned without
#'     \code{assurance}, with \eqn{P_N \ge} \code{assurance}, found by
#'     doubling a step upward and then bisecting. \code{assurance_attained}
#'     is \eqn{P_N} at that \emph{N}.
#' }
#' The probability is exact, up to numerical integration, under bivariate
#' normality. It was checked against the density of \emph{r} computed
#' from its hypergeometric series, and against the proportion of
#' intervals no wider than the target among samples of bivariate normal
#' data. The search does not go below the plan without \code{assurance},
#' so a plan with \code{assurance} is never smaller than the plan without
#' it. Near \eqn{\rho = 0} a slightly smaller \emph{N} can attain the
#' probability, because there the interval at \eqn{r = 0} is barely wider
#' than the target (at \eqn{\rho = 0} and a width of .10 the probability is
#' .917 at 1537, one below the plan without \code{assurance}). At and above
#' the plan without \code{assurance}, \eqn{P_N} did not decrease from one
#' \emph{N} to the next in any of the 15 planning sets scanned (\eqn{|\rho|}
#' of 0, .1, .5, .9, and .99 with widths of .05, .1, and .2). The result
#' depends on \code{true_cor}
#' only through \eqn{|\rho|}. At \code{true_cor = 0} the plan without
#' \code{assurance} already makes every interval as narrow as the target,
#' since the interval is widest at \eqn{r = 0}; \code{assurance} then
#' adds nothing and \code{assurance_attained} is 1. As \eqn{|\rho|}
#' grows the width becomes more variable from sample to sample and the
#' step adds more: for a width of .10 and an assurance of .80 the plan
#' moves from 62 to 83 at \eqn{\rho = .90}.
#'
#' \strong{Why the step is not a chi square factor.} Multiplying the
#' sample size planned without \code{assurance} by
#' \eqn{\chi^2_{\gamma, N - 2} / (N - 2)}, the \eqn{\gamma} quantile of a
#' chi square on \eqn{N - 2} degrees of freedom over its degrees of
#' freedom, with \eqn{\gamma} the \code{assurance}, does not serve here.
#' That factor is the \eqn{\gamma} quantile of the ratio of a sample
#' variance to its population value, so it applies to an interval whose
#' width is proportional to a sample standard deviation; the width of the
#' interval for a correlation depends on the sample through \emph{r}
#' instead. At \code{true_cor = .90}, a width of .10, and an assurance of
#' .80 the factor gives 72, where the probability is .672 (83 are needed);
#' with an assurance of .99 it gives 92, where the probability is .887
#' (122 are needed). At \code{true_cor = 0} it gives more than the plan
#' without \code{assurance}, which already attains a probability of 1.
#'
#' \strong{The planning value matters least near zero.} At a fixed
#' sample size the back-transformed width is largest at \eqn{\rho = 0}
#' and shrinks as \eqn{|\rho|} grows, so a planning value closer to zero
#' yields a larger, more conservative sample size. When little is known
#' about the population correlation, \code{true_cor = 0} gives the sample
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
#' # Plan N so the 95% CI on the Pearson correlation has full width
#' # at most 0.20, when the anticipated correlation is 0.30.
#' ss_aipe_r(true_cor = 0.30, width = 0.20)
#'
#' # A narrower target width requires a larger sample size.
#' ss_aipe_r(true_cor = 0.30, width = 0.10)
#'
#' # The smallest N at which the realized interval is no wider than
#' # 0.20 with probability at least .80. The last row is the probability
#' # at the returned N.
#' ss_aipe_r(true_cor = 0.30, width = 0.20, assurance = 0.80)
#'
#' # At a large correlation the width varies more from sample to sample,
#' # and the assurance step adds more cases.
#' ss_aipe_r(true_cor = 0.90, width = 0.10)
#' ss_aipe_r(true_cor = 0.90, width = 0.10, assurance = 0.80)
#'
#' # Planning at true_cor = 0 gives the sample size that suffices for any
#' # population correlation, since the interval is widest there.
#' ss_aipe_r(true_cor = 0, width = 0.20)
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

ss_aipe_r <- function(true_cor, width,
                      conf_level = 0.95,
                      assurance = NULL) {
  .dmar_check_conf_level(conf_level)
  .dmar_check_assurance_single(assurance)
  # The arithmetic below is written in rho, the planning value.
  rho <- true_cor
  if (!is.numeric(rho) || length(rho) != 1L || !is.finite(rho) ||
      abs(rho) >= 1)
    .dmar_refuse_argument(
      "'true_cor' must be a single value in (-1, 1).", arg = "true_cor")
  if (!is.numeric(width) || length(width) != 1L || !is.finite(width) ||
      width <= 0)
    .dmar_refuse_argument(
      "'width' must be a single positive number.", arg = "width")
  if (!is.numeric(conf_level) || length(conf_level) != 1L ||
      conf_level <= 0 || conf_level >= 1)
    .dmar_refuse_argument("'conf_level' must be in (0, 1).", arg = "conf_level")
  # The floor on a tail area (R/dmar_conditions.R), which ci_r(), the
  # interval planned for, keeps.
  .dmar_check_tail_floor((1 - conf_level) / 2, (1 - conf_level) / 2,
                         conf_level = conf_level)
  if (!is.null(assurance) &&
      (!is.finite(assurance) || assurance <= 0.5 || assurance >= 1))
    .dmar_refuse_argument(
      "'assurance' must be in (0.5, 1).", arg = "assurance")

  z_alpha <- stats::qnorm(1 - (1 - conf_level) / 2)
  # The width and the probabilities below depend on rho only through its
  # magnitude; computing with |rho| makes the plans for rho and -rho
  # identical to the last bit.
  rho_abs <- abs(rho)
  z_rho   <- atanh(rho_abs)

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
  n0 <- 3 + ceiling(4 * z_alpha^2 * (1 - rho_abs^2)^2 / width^2)
  n  <- max(n_min, n0)

  # Exact refinement: the back-transformed width is strictly decreasing in n,
  # so walk the integer n in whichever direction the closed form missed and
  # stop at the smallest n whose exact width is at or below the target.
  while (full_back(n) > width) n <- n + 1
  while (n > n_min && full_back(n - 1) <= width) n <- n - 1

  # The assurance step: the smallest n, at or above the plan just made, at
  # which the probability of an interval no wider than the target is at
  # least `assurance`. The interval is the partial correlation's Fisher's Z
  # interval with no variables partialed out (J = 0: standard error
  # 1 / sqrt(n - 3), and the sample correlation of n observations), so
  # .partial_r_assurance_n() in R/partial_r_internals.R computes it: the
  # width event is |r| >= r*, and its probability comes from the exact
  # distribution function of r.
  attained <- NULL
  if (!is.null(assurance)) {
    step <- .partial_r_assurance_n(n, rho_abs, J = 0,
                                   half_width = width / 2, z = z_alpha,
                                   fisher_z = TRUE, assurance = assurance)
    n <- step$n
    attained <- step$attained
  }

  out <- data.frame(
    term  = c("necessary_N", "expected_width", "true_cor",
              "width", "conf_level",
              if (!is.null(assurance)) c("assurance", "assurance_attained")),
    value = c(n, full_back(n), rho, width, conf_level,
              if (!is.null(assurance)) c(assurance, attained)),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  .as_dmar_tbl(out, conf_level = conf_level, subclass = "dmar_ss_aipe")
}
