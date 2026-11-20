# Sample size for AIPE on a partial correlation.
#' Sample Size for AIPE on a Partial Correlation
#'
#' Takes a planning value of the population partial correlation
#' \eqn{\rho_{XY \cdot Z_1 \cdots Z_{p-1}}}, the number of predictors
#' (\emph{X} and the variables partialed out), and a desired confidence
#' interval width, and returns the smallest
#' sample size at which the interval is expected to be that narrow
#' (accuracy in parameter estimation, AIPE; Kelley, 2008). With
#' \code{assurance} it
#' returns the smallest sample size at which the realized interval is that
#' narrow with at least the stated probability.
#'
#' @param true_cor Population partial correlation
#'   \eqn{\rho_{XY \cdot Z_1 \cdots Z_{p-1}}} (written \eqn{\rho} in Details)
#'   the sample size is planned for, in \eqn{(-1, 1)}. The planner treats this value as the
#'   truth; \code{\link{ss_aipe_partial_r_sensitivity}} shows what the plan
#'   delivers when the truth (\code{true_cor} there) differs from the value
#'   planned with (\code{supposed_cor} there).
#' @param p Number of predictors in the regression of \emph{Y} on \emph{X}
#'   and the variables partialed out, so \code{p - 1} variables are
#'   partialed out and \code{p} is at least 2. A categorical variable with
#'   \emph{a} levels counts as \emph{a} - 1 predictors. It is the \code{p}
#'   of \code{\link{ss_power_R2}}, \code{\link{ci_R2}}, and
#'   \code{\link{ci_reg_coef}}: the correlation concerns one predictor of a
#'   regression with \code{p} predictors. In the notation of a model
#'   comparison, \eqn{q} is the number of predictors removed from the full
#'   model to give the restricted model (\eqn{q = p_{full} - p_{restricted}},
#'   the predictors being tested), and the number of variables partialed out
#'   is \eqn{p - q}; this correlation concerns one predictor, so \eqn{q = 1}
#'   and \eqn{p - 1} variables are partialed out.
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
#' @param fisher_z Logical. If \code{TRUE}, the plan is for the interval
#'   formed on Fisher's \emph{Z} scale (variance \eqn{1/(N - p - 2)})
#'   and back-transformed with \eqn{\tanh}. If \code{FALSE} (default), the
#'   plan is for the interval formed on the raw scale from the first-order
#'   variance.
#' @param assurance Optional. The probability that the realized interval
#'   is no wider than \code{width}. The function returns the smallest
#'   sample size, at or above the sample size planned without it, at which
#'   that probability is at least \code{assurance}. A single number above
#'   0.5 and below 1; when \code{NULL} (the default) the plan is for the
#'   expected width only.
#'
#' @return A \code{data.frame} with the columns \code{term} and
#'   \code{value} and one row per term:
#'   \describe{
#'     \item{\code{necessary_N}}{the planned total sample size.}
#'     \item{\code{expected_width}}{the full width of the interval at
#'       \code{necessary_N} when the sample partial correlation equals
#'       \code{true_cor}.}
#'     \item{\code{true_cor}, \code{p}, \code{width},
#'       \code{conf_level}}{the inputs, with \code{width} the
#'       \code{width} as given.}
#'     \item{\code{assurance}}{the requested probability; present only
#'       when \code{assurance} was given.}
#'     \item{\code{assurance_attained}}{the probability, at
#'       \code{necessary_N} and \code{true_cor}, that the realized interval is
#'       no wider than the target; present only when \code{assurance} was
#'       given. It is at least \code{assurance} and exceeds it because the
#'       sample size is a whole number.}
#'   }
#'
#' @details
#' \strong{Estimand and interval.} The quantity planned for is the
#' population partial correlation between \eqn{X} and \eqn{Y} given
#' \eqn{Z_1, \ldots, Z_{p-1}}, for independent observations from a
#' multivariate normal population. The interval is two-sided at
#' \code{conf_level} and is one of two Wald intervals, chosen by
#' \code{fisher_z}; the plan is for the one chosen, and a study that
#' reports the other is not the study that was planned. Below,
#' \eqn{\rho} is that partial correlation (\code{true_cor}) and
#' \eqn{1 - \alpha} is \code{conf_level}.
#'
#' \strong{Raw scale, first-order variance (\code{fisher_z = FALSE}).}
#' The first-order variance of a simple correlation is
#' \eqn{(1 - \rho^2)^2 / N} (Olkin & Finn, 1995, Equation 3). The sample
#' partial correlation given \eqn{p - 1} variables is distributed as a simple
#' correlation from \eqn{N - (p - 1)} observations (Fisher, 1924, p. 330,
#' for normally distributed variables), and the functions here use
#' \eqn{(1 - \rho^2)^2 / (N - p)}, which differs from that at order
#' \eqn{1/N^2} (see \code{\link{var_partial_r}}). The interval is
#' \eqn{r \pm z_{1 - \alpha/2} (1 - r^2) / \sqrt{N - p}}, with half-width
#' \deqn{w_{1/2}(r, N) \;=\; z_{1 - \alpha/2}\,
#'   \frac{1 - r^2}{\sqrt{N - p}}.}
#' Setting \eqn{r = \rho} and solving for the smallest whole \eqn{N} with
#' \eqn{w_{1/2} \le} the target half-width gives
#' \deqn{N \;=\; p + \Big\lceil
#'         z_{1 - \alpha/2}^2 \, (1 - \rho^2)^2 / w_{1/2}^{2}
#'       \Big\rceil.}
#'
#' \strong{Fisher's \emph{Z} scale (\code{fisher_z = TRUE}).} Bonett
#' and Wright (2000, Section 3, pp. 24--25) plan the sample size for a simple
#' correlation so that the interval formed on Fisher's \emph{Z} scale
#' (Fisher, 1921, 1924) and back-transformed has the desired width. The same
#' plan is applied here to the partial correlation, with the standard error
#' \eqn{1/\sqrt{N - p - 2}} in place of \eqn{1/\sqrt{N - 3}}. With
#' \eqn{Z = \mathrm{arctanh}(r)} and standard error
#' \eqn{1/\sqrt{N - p - 2}}, the interval is
#' \eqn{\tanh(Z \pm z_{1 - \alpha/2} / \sqrt{N - p - 2})}, and half its
#' width is
#' \deqn{w_{1/2}(r, N) \;=\; \tfrac{1}{2}\Big\{
#'   \tanh\!\Big(\mathrm{arctanh}(r) + \frac{z_{1 - \alpha/2}}{\sqrt{N - p - 2}}\Big)
#'   - \tanh\!\Big(\mathrm{arctanh}(r) - \frac{z_{1 - \alpha/2}}{\sqrt{N - p - 2}}\Big)
#'   \Big\}.}
#' The planner finds the smallest \eqn{N} of at least \eqn{p + 3} at which
#' this is no larger than the target at \eqn{r = \rho}, by a one-dimensional
#' search (the half-width decreases in \eqn{N}). This interval stays inside
#' \eqn{(-1, 1)} and is generally the more accurate of the two when
#' \eqn{N - p} is small or \eqn{|\rho|} is large.
#'
#' \strong{The assurance step.} A plan without \code{assurance} makes the
#' interval as narrow as the target when the sample partial correlation
#' equals \code{true_cor}; in about half of all samples it is wider. With
#' \code{assurance} the function finds the sample size at which the
#' realized interval is no wider than the target with at least that
#' probability, from the sampling distribution of the sample partial
#' correlation:
#' \enumerate{
#'   \item At a given \eqn{N} the width of either interval depends on the
#'     sample only through \eqn{r}, and both half-widths above decrease as
#'     \eqn{|r|} grows. The realized interval is therefore no wider than
#'     the target exactly when \eqn{|r| \ge r^*(N)}, where \eqn{r^*(N)} is
#'     the \eqn{|r|} at which the half-width equals the target half-width
#'     \eqn{h}: on the raw scale
#'     \eqn{r^*(N) = \sqrt{1 - h \sqrt{N - p} / z_{1 - \alpha/2}}}, and
#'     on Fisher's scale the root of \eqn{w_{1/2}(r, N) = h}. When the
#'     half-width at \eqn{r = 0}, the widest the interval can be, already
#'     meets the target, every sample does and the probability is 1.
#'   \item Under multivariate normality the sample partial correlation
#'     given \eqn{p - 1} variables from \eqn{N} observations is distributed as
#'     a simple correlation from \eqn{m = N - (p - 1)} observations. Its
#'     distribution function is computed as
#'     \deqn{F(c) \;=\; \int_0^1 T_{2m - 3}\!\Big(\sqrt{2m - 3}\,
#'       \Big\{\frac{c\,\sqrt{1 - b}}{\sqrt{1 - c^2}} -
#'       \frac{\rho\,\sqrt{b}}{\sqrt{1 - \rho^2}}\Big\}\Big)\,
#'       g(b)\, db,}
#'     where \eqn{T_{2m - 3}} is the distribution function of Student's
#'     \emph{t} on \eqn{2m - 3} degrees of freedom and \eqn{g} is the
#'     density of the beta distribution with parameters \eqn{(m - 1)/2} and
#'     \eqn{(m - 2)/2}. The form follows from writing
#'     \eqn{r / \sqrt{1 - r^2}} as
#'     \eqn{(Z + \theta \sqrt{U}) / \sqrt{V}}, with \eqn{Z} standard
#'     normal, \eqn{U} and \eqn{V} chi square on \eqn{m - 1} and
#'     \eqn{m - 2} degrees of freedom, all independent, and
#'     \eqn{\theta = \rho / \sqrt{1 - \rho^2}}, and then conditioning on
#'     \eqn{U / (U + V)}. It uses the central \emph{t} only, because R's
#'     noncentral \emph{t} is an approximation at the noncentrality values
#'     that plans for narrow intervals reach.
#'   \item The probability at \eqn{N} is
#'     \eqn{P_N = 1 - F(r^*(N)) + F(-r^*(N))}.
#'   \item The function returns the smallest \eqn{N}, at or above the
#'     sample size planned without \code{assurance}, with
#'     \eqn{P_N \ge} \code{assurance}, found by doubling a step upward and
#'     then bisecting. \code{assurance_attained} is \eqn{P_N} at that
#'     \eqn{N}.
#' }
#' The probability is exact, up to numerical integration, under
#' multivariate normality for the interval that is planned. It was checked
#' by simulation in two ways: in 112 planning sets of 200,000 simulated
#' correlations each, the proportion of intervals no wider than the target
#' at the planned sample size was within .0023 of \eqn{P_N}; and in 16
#' further planning sets simulated from regression data (20,000 data sets
#' each, with the variables partialed out by least squares) it was within
#' .004 of \eqn{P_N}, about two Monte Carlo standard errors, and between
#' .002 and .012 above the request.
#'
#' Four properties of the step are worth knowing.
#' \itemize{
#'   \item The answer is never below the sample size planned without
#'     \code{assurance}. Below that sample size \eqn{P_N} is not monotone
#'     in \eqn{N} (at a very small \eqn{N} the sample correlation is
#'     variable enough to exceed \eqn{r^*(N)} by chance, and does so less
#'     often as \eqn{N} grows, before the narrowing of the interval takes
#'     over), so a smaller \eqn{N} that
#'     happens to meet the request is not a plan. At or above it
#'     \eqn{P_N} did not decrease in any of 72 planning sets scanned over
#'     every whole \eqn{N}.
#'   \item When \code{true_cor} is 0 the plan without \code{assurance} already
#'     makes the interval as narrow as the target in every sample, since
#'     the interval is widest at \eqn{r = 0}; \code{assurance} then adds
#'     nothing and \code{assurance_attained} is 1. Near 0 the same holds
#'     approximately and the step adds little.
#'   \item As \eqn{|\rho|} grows the width becomes more variable from
#'     sample to sample and the step adds more: for a width of .10 with
#'     one variable partialed out and an assurance of .80, the plan moves
#'     from 867 to 915 at \eqn{\rho = .50} and from 121 to 151 at
#'     \eqn{\rho = .85}.
#'   \item \code{assurance} must be above 0.5 and below 1; other values
#'     are refused.
#' }
#'
#' \strong{Why the step is not a chi square factor.} Multiplying the
#' sample size planned without \code{assurance} by a chi square quantile
#' over its degrees of freedom is the assurance step for an interval whose
#' width is proportional to a sample standard deviation (Kelley & Maxwell,
#' 2003); the width of an interval for a correlation depends on the sample
#' through \eqn{r} instead, and the factor gives too few cases when
#' \eqn{|\rho|} is large and too many when it is small. For one variable
#' partialed out (\eqn{p = 2}), a width of .10, and a requested .80 it
#' gives 426 at \eqn{\rho = .70}, where the probability is .675 (448 are
#' needed), and 1318 at \eqn{\rho = .30}, where 1310 suffice; in the 16
#' planning sets simulated from regression data its sample sizes delivered
#' from .14 below the request to .15 above it.
#'
#' \strong{When to use partial or simple correlation planning.} Use this
#' function when the target is the population correlation between \eqn{X}
#' and \eqn{Y} after statistically controlling for
#' \eqn{Z_1, \ldots, Z_{p-1}}. For the simple Pearson correlation, see
#' \code{\link{ss_aipe_r}}. \code{\link{ss_aipe_partial_r_sensitivity}}
#' simulates the plan, including under a population value other than the
#' planning value.
#'
#' @references
#' Algina, J., & Olejnik, S. (2003). Sample size tables for correlation
#'   analysis with applications in partial correlation and multiple
#'   regression analysis. \emph{Multivariate Behavioral Research, 38}(3),
#'   309--323. \doi{10.1207/s15327906mbr3803_02}
#'
#' Bonett, D. G., & Wright, T. A. (2000). Sample size requirements for
#'   estimating Pearson, Kendall and Spearman correlations.
#'   \emph{Psychometrika, 65}(1), 23--28. \doi{10.1007/BF02294183}
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
#'   \code{\link{ss_aipe_partial_r_sensitivity}},
#'   \code{\link{ss_aipe_semipartial_r}}, \code{\link{ss_aipe_r}}
#'
#' @examples
#' # 1. Plan N so the 95% CI on rho_XY.Z (two variables partialed out, so
#' #        p = 3 predictors) has full width <= 0.20, when the planning
#' #        value of the partial correlation is 0.30.
#' ss_aipe_partial_r(true_cor = 0.30, p = 3, width = 0.20)
#'
#' # 2. The same problem for the interval formed on Fisher's Z scale
#' #        (Bonett & Wright, 2000).
#' ss_aipe_partial_r(true_cor = 0.30, p = 3, width = 0.20, fisher_z = TRUE)
#'
#' # 3. With an assurance of .80: the smallest N at which the realized
#' #        interval is no wider than 0.20 with probability at least .80.
#' #        The last row is the probability at the returned N.
#' ss_aipe_partial_r(true_cor = 0.30, p = 3, width = 0.20, assurance = 0.80)
#'
#' # 4. A larger correlation needs fewer cases for the same width, and
#' #        the assurance step adds proportionally more of them.
#' ss_aipe_partial_r(true_cor = 0.70, p = 3, width = 0.20)
#' ss_aipe_partial_r(true_cor = 0.70, p = 3, width = 0.20, assurance = 0.80)
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

ss_aipe_partial_r <- function(true_cor, p, width,
                              which_width = c("Full", "Lower", "Upper"),
                              conf_level = 0.95,
                              fisher_z = FALSE,
                              assurance = NULL) {
  .dmar_check_conf_level(conf_level)
  if (!missing(width)) .dmar_check_width(width)
  if (missing(p))
    .dmar_refuse_argument(
      "'p' (the number of predictors) is required.", arg = "p")
  p <- .check_partial_p(p)
  # The arithmetic below is written in J, the number of variables partialed
  # out, and in rho: N - J - 1 = N - p, N - J - 3 = N - p - 2, and the
  # smallest N on Fisher's Z scale is J + 4 = p + 3.
  J <- p - 1
  rho <- true_cor
  which_width <- match.arg(which_width)
  if (!is.numeric(rho) || length(rho) != 1L || abs(rho) >= 1)
    .dmar_refuse_argument(
      "'true_cor' must be a single value in (-1, 1).", arg = "true_cor")
  if (!is.numeric(width) || length(width) != 1L || width <= 0)
    .dmar_refuse_argument(
      "'width' must be a single positive number.", arg = "width")
  if (conf_level <= 0 || conf_level >= 1)
    .dmar_refuse_argument("'conf_level' must be in (0, 1).", arg = "conf_level")
  if (!is.null(assurance) &&
      (!is.numeric(assurance) || length(assurance) != 1L ||
       !is.finite(assurance) || assurance <= 0.5 || assurance >= 1))
    .dmar_refuse_argument("'assurance' must be in (0.5, 1).", arg = "assurance")

  half_width <- if (which_width == "Full") width / 2 else width
  z_alpha    <- stats::qnorm(1 - (1 - conf_level) / 2)

  if (!fisher_z) {
    # Raw scale, first-order variance: n - J - 1 = z^2 (1 - rho^2)^2 / w^2
    n <- J + 1 + ceiling((z_alpha^2 * (1 - rho^2)^2) / (half_width^2))
  } else {
    # Fisher's Z scale: the smallest n whose back-transformed interval half-width
    # is <= the target on the raw scale. The half-width decreases monotonically
    # in n toward 0, so a finite n always exists; solve for it directly. The old
    # fixed grid (J + 5 to J + 5000) silently returned its minimum, J + 5,
    # whenever the target fell past the grid's right end, since which.max() of an
    # all-FALSE vector is 1.
    half_back <- function(nn) .partial_r_half_width(rho, nn, J, z_alpha, TRUE)
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

  # The assurance step: the smallest n, at or above the plan just made, at
  # which the probability of an interval no wider than the target is at
  # least `assurance` (.partial_r_assurance_n() in R/partial_r_internals.R,
  # where the distribution function and the search are described).
  attained <- NULL
  if (!is.null(assurance)) {
    step <- .partial_r_assurance_n(n, rho, J, half_width, z_alpha, fisher_z,
                                   assurance)
    n <- step$n
    attained <- step$attained
  }

  # Expected width at the recommended n: the full width of the planned
  # interval at a sample partial correlation equal to rho.
  expected_w <- 2 * .partial_r_half_width(rho, n, J, z_alpha, fisher_z)

  out <- data.frame(
    term  = c("necessary_N", "expected_width", "true_cor", "p",
              "width", "conf_level",
              if (!is.null(assurance)) c("assurance", "assurance_attained")),
    value = c(n, expected_w, rho, p, width, conf_level,
              if (!is.null(assurance)) c(assurance, attained)),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  .as_dmar_tbl(out, conf_level = conf_level, subclass = "dmar_ss_aipe")
}
