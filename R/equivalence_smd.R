# TOST on the standardized mean difference.
#' Equivalence Test for the Standardized Mean Difference via Two One-Sided Tests (TOST)
#'
#' Performs a two one-sided tests procedure (Schuirmann, 1987) for
#' equivalence between two independent groups on the standardized
#' mean difference (Cohen's \emph{d}) scale. The null hypothesis is
#' that the true \eqn{\delta} lies outside the user-specified
#' equivalence bounds \eqn{[-\delta_L, \delta_U]}; the alternative is
#' that \eqn{\delta} lies inside them. The test is the joint pair of
#' one-sided \emph{t} tests: \eqn{H_{0,L}: \delta \le -\delta_L} vs
#' \eqn{H_{1,L}: \delta > -\delta_L}, and \eqn{H_{0,U}: \delta \ge
#' \delta_U} vs \eqn{H_{1,U}: \delta < \delta_U}. Equivalence is
#' declared when \emph{both} null hypotheses are rejected at level
#' \eqn{\alpha}.
#'
#' @param x,y Numeric vectors of observations from the two groups.
#'   Alternatively, supply \code{smd}, \code{n_1}, \code{n_2} via the
#'   summary-statistic interface.
#' @param smd Observed standardized mean difference (Cohen's \emph{d}).
#'   Required if \code{x} and \code{y} are not supplied.
#' @param n_1,n_2 Group sample sizes. Required if \code{x} and \code{y}
#'   are not supplied.
#' @param delta_lower,delta_upper Equivalence bounds on the \emph{d}
#'   scale. Both must be positive; the equivalence region is
#'   \eqn{[-\delta_L, +\delta_U]}. If only \code{delta_upper} is
#'   supplied, the bounds are symmetric:
#'   \eqn{[-\delta_U, +\delta_U]}.
#' @param alpha_level One-sided significance level for each of the two
#'   tests. Default \code{0.05}. The TOST is then equivalent to a
#'   100(1 - 2\eqn{\alpha})\% CI on \eqn{\delta} lying inside the
#'   equivalence bounds.
#'
#' @return A \code{data.frame} with rows for the observed \emph{d},
#'   the two one-sided test statistics (\code{t_lower}, \code{t_upper})
#'   and their degrees of freedom (\code{df}), the two one-sided
#'   \emph{p}-values (\code{p_lower}, \code{p_upper}), the joint TOST
#'   \emph{p}-value (the larger of the two), the
#'   100(1 - 2\eqn{\alpha})\% CI on \eqn{\delta}, the equivalence
#'   bounds, and a binary decision flag (\code{equivalent}:
#'   1 = equivalent, 0 = not).
#'
#' @details
#' \strong{Schuirmann's TOST.} The TOST procedure tests
#' \eqn{H_0: \delta \le -\delta_L \cup \delta \ge \delta_U} against
#' \eqn{H_1: -\delta_L < \delta < \delta_U}. Both component tests are
#' rejected (and equivalence is declared) when the 100(1 -
#' 2\eqn{\alpha})\% CI on \eqn{\delta} lies entirely inside
#' \eqn{[-\delta_L, \delta_U]}.
#'
#' \strong{Critical insight.} A non-significant conventional NHST
#' (\eqn{t} test of \eqn{\delta = 0}) is \emph{not} evidence of
#' equivalence; it only means we cannot reject \eqn{\delta = 0}. TOST
#' inverts the testing logic so that "no meaningful effect" is the
#' alternative, not the null.
#'
#' \strong{Choosing \eqn{\delta_L} and \eqn{\delta_U}.} The equivalence
#' bounds must be set \emph{before} data collection, based on what
#' constitutes the smallest effect size of practical interest (Lakens,
#' Scheel, & Isager, 2018). Common choices in the literature are 0.2 or
#' 0.3, but the bound should reflect domain-specific meaningfulness.
#'
#' \strong{Connection to the CI.} The TOST rejection at level
#' \eqn{\alpha} is equivalent to the 100(1 - 2\eqn{\alpha})\%
#' Cohen's-\emph{d} CI (i.e., 90\% for the default \eqn{\alpha = 0.05})
#' lying entirely inside the equivalence region. This is the
#' "two-one-sided" equivalence and matches the Westlake (1972)
#' rationale for symmetric bioequivalence CIs.
#'
#' \strong{Standard error and approximation.} Each one-sided component is
#' a Wald \emph{t} test, \eqn{t = (\hat d \mp \delta) / \mathrm{SE}(\hat d)},
#' referred to a central \emph{t} distribution on \eqn{n_1 + n_2 - 2} degrees
#' of freedom, with the Hedges and Olkin (1985) large-sample standard error
#' \eqn{\mathrm{SE}(\hat d) = \sqrt{(n_1 + n_2) / (n_1 n_2) +
#' \hat d^{\,2} / (2 (n_1 + n_2))}}. This is the asymptotic form of
#' Schuirmann's (1987) two one-sided tests applied on the standardized scale;
#' it is accurate in moderate-to-large samples but is an approximation to the
#' exact noncentral \emph{t} inversion used by \code{\link{ci_smd}}, and the
#' two can differ in small samples.
#'
#' @references
#' Hedges, L. V., & Olkin, I. (1985). \emph{Statistical methods for
#'   meta-analysis}. Academic Press.
#'
#' Lakens, D. (2017). Equivalence tests: A practical primer for
#'   \emph{t} tests, correlations, and meta-analyses. \emph{Social
#'   Psychological and Personality Science, 8}(4), 355--362.
#'   \doi{10.1177/1948550617697177}
#'
#' Lakens, D., Scheel, A. M., & Isager, P. M. (2018). Equivalence
#'   testing for psychological research: A tutorial. \emph{Advances in
#'   Methods and Practices in Psychological Science, 1}(2), 259--269.
#'   \doi{10.1177/2515245918770963}
#'
#' Schuirmann, D. J. (1987). A comparison of the two one-sided tests
#'   procedure and the power approach for assessing the equivalence of
#'   average bioavailability. \emph{Journal of Pharmacokinetics and
#'   Biopharmaceutics, 15}(6), 657--680.
#'
#' Westlake, W. J. (1972). Use of confidence intervals in analysis of
#'   comparative bioavailability trials. \emph{Journal of Pharmaceutical
#'   Sciences, 61}(8), 1340--1341.
#'
#' @seealso \code{\link{equivalence_r}}, \code{\link{ss_aipe_equivalence_smd}},
#'   \code{\link{ci_smd}}, \code{\link{smd}}
#'
#' @examples
#' # 1. Two groups, raw data, equivalence bound delta = 0.4:
#' set.seed(113)
#' x <- rnorm(50, 100, 15); y <- rnorm(50, 101, 15)
#' equivalence_smd(x = x, y = y, delta_upper = 0.4)
#'
#' # 2. Summary-statistic interface: published d = 0.05, n = 40 per group.
#' equivalence_smd(smd = 0.05, n_1 = 40, n_2 = 40, delta_upper = 0.3)
#'
#' # 3. Asymmetric bounds (acceptable -0.2 to +0.5):
#' equivalence_smd(smd = 0.10, n_1 = 80, n_2 = 80,
#'          delta_lower = 0.2, delta_upper = 0.5)
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @keywords htest
#'
#' @family hypothesis tests
#'
#' @family equivalence testing
#'
#' @export

equivalence_smd <- function(x = NULL, y = NULL,
                     smd = NULL, n_1 = NULL, n_2 = NULL,
                     delta_lower = NULL, delta_upper = NULL,
                     alpha_level = 0.05) {
  if (is.null(delta_upper))
    stop("'delta_upper' must be specified (the upper equivalence bound).")
  if (!is.numeric(delta_upper) || delta_upper <= 0)
    stop("'delta_upper' must be a positive number.")
  if (is.null(delta_lower)) delta_lower <- delta_upper
  if (!is.numeric(delta_lower) || delta_lower <= 0)
    stop("'delta_lower' must be a positive number.")
  if (alpha_level <= 0 || alpha_level >= 0.5)
    stop("'alpha_level' must be in (0, 0.5).")

  if (!is.null(x) && !is.null(y)) {
    if (!is.numeric(x) || !is.numeric(y))
      stop("'x' and 'y' must be numeric vectors.")
    x <- x[!is.na(x)]; y <- y[!is.na(y)]
    n_1 <- length(x); n_2 <- length(y)
    if (n_1 < 2L || n_2 < 2L)
      stop("Each group needs at least 2 non-missing observations.")
    m_1 <- mean(x); m_2 <- mean(y)
    s_p <- sqrt(((n_1 - 1) * stats::var(x) +
                 (n_2 - 1) * stats::var(y)) / (n_1 + n_2 - 2))
    if (s_p == 0) stop("Pooled SD is zero; d is undefined.")
    smd <- (m_1 - m_2) / s_p
  } else if (is.null(smd) || is.null(n_1) || is.null(n_2)) {
    stop("Supply either ('x', 'y') or ('smd', 'n_1', 'n_2').")
  } else {
    if (!is.numeric(smd) || length(smd) != 1L)
      stop("'smd' must be a single numeric value.")
    if (n_1 < 2 || n_2 < 2)
      stop("'n_1' and 'n_2' must each be at least 2.")
  }

  df <- n_1 + n_2 - 2L
  se_d <- sqrt((n_1 + n_2) / (n_1 * n_2) + smd^2 / (2 * (n_1 + n_2)))

  t_lower <- (smd - (-delta_lower)) / se_d
  t_upper <- (smd -   delta_upper)  / se_d
  p_lower <- stats::pt(t_lower, df, lower.tail = FALSE)
  p_upper <- stats::pt(t_upper, df, lower.tail = TRUE)
  p_tost  <- max(p_lower, p_upper)

  # (1 - 2 alpha_level) CI on delta:
  t_crit  <- stats::qt(1 - alpha_level, df)
  ci_lo   <- smd - t_crit * se_d
  ci_hi   <- smd + t_crit * se_d

  equivalent <- as.integer((ci_lo > -delta_lower) && (ci_hi < delta_upper))

  out <- data.frame(
    term  = c("smd", "t_lower", "t_upper", "df",
              "p_lower", "p_upper", "p_tost",
              "lower_limit", "upper_limit",
              "delta_lower", "delta_upper",
              "equivalent"),
    value = c(smd, t_lower, t_upper, df,
              p_lower, p_upper, p_tost,
              ci_lo, ci_hi,
              -delta_lower, delta_upper,
              equivalent),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  .as_dmar_tbl(out, p_terms = c("p_lower", "p_upper", "p_tost"),
               conf_level = 1 - 2 * alpha_level)
}
