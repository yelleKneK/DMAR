# TOST and noninferiority testing on a linear contrast.
#' Equivalence and Noninferiority Tests for a Linear Contrast via Two One-Sided Tests (TOST)
#'
#' Performs the two one-sided tests procedure (Schuirmann, 1987) for
#' equivalence, and the companion one-sided noninferiority test, for a
#' linear contrast of group means \eqn{\psi = \sum_j c_j \mu_j} in a
#' fixed effects design with one pooled error term. The equivalence
#' null hypothesis is that the population contrast lies outside the
#' user-specified bounds, \eqn{H_0: \psi \le -\delta_L \cup \psi \ge
#' \delta_U}, against the alternative \eqn{H_1: -\delta_L < \psi <
#' \delta_U}. The noninferiority null is \eqn{H_0: \psi \le -\delta_L}
#' against \eqn{H_1: \psi > -\delta_L}. Following Chattopadhyay,
#' Bandyopadhyay, Kelley, and Padalunkal (2025), the equivalence
#' decision is read from the 100(1 - 2\eqn{\alpha})\% confidence
#' interval: equivalence is declared when the whole interval lies
#' inside \eqn{(-\delta_L, \delta_U)}, and noninferiority when the
#' interval's lower limit exceeds \eqn{-\delta_L}.
#'
#' @param means A vector of group means. Supply together with
#'   \code{s_anova}, \code{c_weights}, and \code{n}. Alternatively,
#'   use the direct interface via \code{psi_hat}, \code{se}, and
#'   \code{df_error}.
#' @param s_anova The standard deviation of the errors from the ANOVA
#'   model (the square root of the mean square error), as in
#'   \code{\link{ci_c}}.
#' @param c_weights The contrast weights. For a mean comparison the
#'   weights must sum to zero, and, so that the bounds are on the raw
#'   scale of the response, the positive weights must sum to 1 and the
#'   negative weights to -1 (use fractional values, not integers).
#'   When \code{benchmark} is supplied, the weights must instead be
#'   nonnegative and sum to 1 (typically a single 1 selecting one
#'   group).
#' @param n Sample sizes per group (if length 1, equal group sizes are
#'   assumed).
#' @param psi_hat The estimated contrast, for the direct interface.
#'   Supply together with \code{se} and \code{df_error} when the
#'   contrast and its standard error have already been computed (for
#'   example, from a fitted model with covariates).
#' @param se The standard error of \code{psi_hat}, for the direct
#'   interface.
#' @param df_error The error degrees of freedom. On the
#'   summary-statistic interface the default is \eqn{N - J}, with
#'   \eqn{J} the number of groups; it must be supplied for designs
#'   with additional factors. Required on the direct interface.
#' @param delta_lower,delta_upper Equivalence bounds on the raw scale
#'   of the response. Both must be positive; the equivalence region is
#'   \eqn{(-\delta_L, +\delta_U)}. If only \code{delta_upper} is
#'   supplied, the bounds are symmetric. Noninferiority uses
#'   \eqn{-\delta_L} alone.
#' @param benchmark An optional known constant to compare against
#'   (for example, a normative or regulatory cutoff). When supplied,
#'   the contrast is \eqn{\sum_j c_j \mu_j - b} with nonnegative
#'   weights summing to 1, and the constant contributes no sampling
#'   variability.
#' @param alpha_level One-sided significance level for each of the two
#'   tests. Default \code{0.05}, so the interval the decisions are
#'   read from is the 90\% CI.
#'
#' @return A \code{data.frame} with rows for the estimated contrast
#'   (\code{psi_hat}), its standard error (\code{se}), the error
#'   degrees of freedom (\code{df}), the two one-sided test statistics
#'   (\code{t_lower}, \code{t_upper}) and their \emph{p}-values
#'   (\code{p_lower}, \code{p_upper}), the joint TOST \emph{p}-value
#'   (\code{p_tost}, the larger of the two), the noninferiority
#'   \emph{p}-value (\code{p_noninferiority}, equal to \code{p_lower}
#'   by construction), the 100(1 - 2\eqn{\alpha})\% confidence limits
#'   (\code{lower_limit}, \code{upper_limit}), the bounds
#'   (\code{delta_lower}, stored as the signed lower bound, and
#'   \code{delta_upper}), and four binary decision flags
#'   (\code{equivalent}, \code{noninferior}, \code{superior},
#'   \code{inferior}; 1 = declared, 0 = not). When all four flags are
#'   0, the interval straddles a bound and the result is inconclusive.
#'   The five-way classification is also attached as the
#'   \code{"verdict"} attribute, one of \code{"Equivalent"},
#'   \code{"Superior"}, \code{"Inferior"}, \code{"Noninferior only"},
#'   or \code{"Inconclusive"}.
#'
#' @details
#' \strong{One pooled error term.} On the summary-statistic interface
#' the standard error is
#' \eqn{\mathrm{SE}(\hat\psi) = s_{\mathrm{anova}} \sqrt{\sum_j c_j^2 / n_j}},
#' the model comparison position of Maxwell, Delaney, and Kelley
#' (2027): every one-degree-of-freedom contrast is judged against the
#' same yardstick, the root mean square error of one model fit to all
#' groups.
#'
#' \strong{The verdict logic.} Reading the 100(1 - 2\eqn{\alpha})\% CI
#' against the bounds: an interval entirely inside
#' \eqn{(-\delta_L, \delta_U)} is \emph{equivalent}; entirely above
#' \eqn{\delta_U} is \emph{superior} (which implies noninferior);
#' entirely below \eqn{-\delta_L} is \emph{inferior}; a lower limit
#' above \eqn{-\delta_L} with an upper limit past \eqn{\delta_U} is
#' \emph{noninferior only}; and an interval straddling a bound is
#' \emph{inconclusive}. An inconclusive result is a statement about
#' precision, not evidence of a difference: only an interval clearing
#' a bound entirely licenses a directional claim.
#'
#' \strong{Why the weights must sum to \eqn{\pm 1}.} A bound stated in
#' raw units of the response is only meaningful if the contrast is
#' itself a simple difference of (weighted) means on that scale, which
#' requires the positive weights to sum to 1 and the negative weights
#' to -1. A weight vector such as \code{c(2, -2)} would silently
#' double the effective bounds, so it is rejected rather than
#' rescaled.
#'
#' \strong{Choosing the bounds.} The bounds must be fixed before the
#' data are examined, on substantive grounds: the smallest difference
#' that would matter (Serlin & Lapsley, 1985; Lakens, Scheel, &
#' Isager, 2018). They are never derived from a standard error, which
#' would make the definition of "close enough" a function of the
#' sample size.
#'
#' \strong{Agreement with emmeans.} The \emph{p}-values reproduce
#' \code{emmeans::test(..., side = "equivalence")} and
#' \code{emmeans::test(..., side = "noninferiority")} with
#' \code{adjust = "none"} to machine precision. Note two emmeans
#' pitfalls the interface here avoids: \code{side = "left"} tests
#' non-superiority, not noninferiority, and \code{trt.vs.ctrl}
#' families silently apply a Dunnett-type adjustment unless
#' \code{adjust = "none"} is passed.
#'
#' @references
#' Chattopadhyay, B., Bandyopadhyay, T., Kelley, K., & Padalunkal,
#'   J. J. (2025). A sequential approach for noninferiority or
#'   equivalence of a linear contrast under cost constraints.
#'   \emph{Psychological Methods, 30}(2), 425--439. \doi{10.1037/met0000570}
#'
#' Lakens, D., Scheel, A. M., & Isager, P. M. (2018). Equivalence
#'   testing for psychological research: A tutorial. \emph{Advances in
#'   Methods and Practices in Psychological Science, 1}(2), 259--269.
#'   \doi{10.1177/2515245918770963}
#'
#' Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). \emph{Designing
#'   experiments and analyzing data: A model comparison perspective}
#'   (4th ed.). Routledge. (See Chapter 4 on individual comparisons of
#'   means.)
#'
#' Schuirmann, D. J. (1987). A comparison of the two one-sided tests
#'   procedure and the power approach for assessing the equivalence of
#'   average bioavailability. \emph{Journal of Pharmacokinetics and
#'   Biopharmaceutics, 15}(6), 657--680.
#'
#' Serlin, R. C., & Lapsley, D. K. (1985). Rationality in
#'   psychological research: The good-enough principle. \emph{American
#'   Psychologist, 40}(1), 73--83.
#'
#' Wellek, S. (2010). \emph{Testing statistical hypotheses of
#'   equivalence and noninferiority} (2nd ed.). Chapman & Hall/CRC.
#'
#' @seealso \code{\link{equivalence_smd}}, \code{\link{equivalence_r}},
#'   \code{\link{ci_c}}, \code{\link{contrast_test}},
#'   \code{\link{power_equivalence_c}},
#'   \code{\link{ss_power_equivalence_c}},
#'   \code{\link{plot_equivalence}}
#'
#' @examples
#' # 1. Two of five groups compared against the reference group, with a
#' #    pooled error term from one model across all five groups.
#' #    Bounds of 5 raw-scale points; alpha_level = .05, so decisions read
#' #    from the 90% CI.
#' equivalence_c(means = c(70.40, 55.61, 51.91, 65.66, 65.12),
#'        s_anova = 15.67,
#'        c_weights = c(-1, 0, 0, 0, 1),
#'        n = c(113, 74, 76, 80, 61),
#'        delta_upper = 5)
#'
#' # 2. The same contrast through the direct interface, as when the
#' #    estimate and standard error come from a model with covariates.
#' res <- equivalence_c(psi_hat = -5.28, se = 2.49, df_error = 399,
#'               delta_upper = 5)
#' res
#' attr(res, "verdict")
#'
#' # 3. A group mean against a fixed benchmark of 68: the constant
#' #    contributes no sampling variability.
#' equivalence_c(means = c(70.40, 55.61, 51.91, 65.66, 65.12),
#'        s_anova = 15.67,
#'        c_weights = c(1, 0, 0, 0, 0),
#'        n = c(113, 74, 76, 80, 61),
#'        benchmark = 68, delta_upper = 5)
#'
#' # 4. Asymmetric bounds: a shortfall of 3 matters, an excess of 8 does.
#' equivalence_c(psi_hat = 1.2, se = 1.1, df_error = 120,
#'        delta_lower = 3, delta_upper = 8)
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @keywords htest
#'
#' @family equivalence testing
#'
#' @export

equivalence_c <- function(means = NULL, s_anova = NULL, c_weights = NULL, n = NULL,
                   psi_hat = NULL, se = NULL, df_error = NULL,
                   delta_lower = NULL, delta_upper = NULL,
                   benchmark = NULL, alpha_level = 0.05) {
  if (is.null(delta_upper))
    stop("'delta_upper' must be specified (the upper equivalence bound).")
  if (!is.numeric(delta_upper) || length(delta_upper) != 1L ||
      !is.finite(delta_upper) || delta_upper <= 0)
    stop("'delta_upper' must be a single positive, finite number.")
  if (is.null(delta_lower)) delta_lower <- delta_upper
  if (!is.numeric(delta_lower) || length(delta_lower) != 1L ||
      !is.finite(delta_lower) || delta_lower <= 0)
    stop("'delta_lower' must be a single positive, finite number.")
  if (!is.numeric(alpha_level) || length(alpha_level) != 1L || alpha_level <= 0 || alpha_level >= 0.5)
    stop("'alpha_level' must be a single number in (0, 0.5).")

  direct  <- !is.null(psi_hat) || !is.null(se)
  summary <- !is.null(means)

  if (direct && summary)
    stop("Supply either ('means', 's_anova', 'c_weights', 'n') or ",
         "('psi_hat', 'se', 'df_error'), not both.")

  if (direct) {
    if (is.null(psi_hat) || is.null(se) || is.null(df_error))
      stop("The direct interface requires 'psi_hat', 'se', and 'df_error'.")
    if (!is.numeric(psi_hat) || length(psi_hat) != 1L || !is.finite(psi_hat))
      stop("'psi_hat' must be a single finite numeric value.")
    if (!is.numeric(se) || length(se) != 1L || !is.finite(se) || se <= 0)
      stop("'se' must be a single positive, finite number.")
    if (!is.numeric(df_error) || length(df_error) != 1L || !is.finite(df_error) ||
        df_error <= 0)
      stop("'df_error' must be a single positive, finite number.")
    if (!is.null(benchmark))
      stop("On the direct interface, subtract the benchmark from 'psi_hat' ",
           "yourself; 'benchmark' applies to the summary-statistic interface.")
    df <- df_error
  } else {
    if (is.null(means) || is.null(s_anova) || is.null(c_weights) || is.null(n))
      stop("The summary-statistic interface requires 'means', 's_anova', ",
           "'c_weights', and 'n'.")
    if (!is.numeric(means) || any(!is.finite(means)))
      stop("'means' must be a finite numeric vector.")
    if (!is.null(df_error) &&
        (!is.numeric(df_error) || length(df_error) != 1L || !is.finite(df_error) ||
         df_error <= 0))
      stop("'df_error' must be a single positive, finite number.")
    if (!is.numeric(s_anova) || length(s_anova) != 1L || !is.finite(s_anova) ||
        s_anova <= 0)
      stop("'s_anova' must be a single positive, finite number (the root mean ",
           "square error, not the error variance).")
    if (length(n) == 1L) n <- rep(n, length(means))
    if (length(n) != length(c_weights) || length(means) != length(c_weights))
      stop("'means', 'c_weights', and 'n' must have the same length.")
    if (!is.numeric(n) || any(!is.finite(n)) || any(n < 2))
      stop("Each group needs at least 2 observations (finite).")

    if (is.null(benchmark)) {
      if (!identical(round(sum(c_weights), 5), 0))
        stop("The sum of the contrast weights ('c_weights') should equal zero.")
      pos <- sum(c_weights[c_weights > 0])
      neg <- sum(c_weights[c_weights < 0])
      if (!isTRUE(all.equal(pos, 1)) || !isTRUE(all.equal(neg, -1)))
        stop("The positive weights must sum to 1 and the negative weights ",
             "to -1, so that the bounds are on the raw scale of the ",
             "response. Use fractional weights, not integers.")
      psi_hat <- sum(c_weights * means)
    } else {
      if (!is.numeric(benchmark) || length(benchmark) != 1L || !is.finite(benchmark))
        stop("'benchmark' must be a single finite numeric value.")
      if (any(c_weights < 0) || !isTRUE(all.equal(sum(c_weights), 1)))
        stop("With a 'benchmark', the weights must be nonnegative and sum ",
             "to 1 (typically a single 1 selecting one group).")
      psi_hat <- sum(c_weights * means) - benchmark
    }

    se <- s_anova * sqrt(sum(c_weights^2 / n))
    df <- if (is.null(df_error)) sum(n) - length(c_weights) else df_error
    if (df <= 0) stop("The error degrees of freedom must be positive.")
  }

  t_lower <- (psi_hat - (-delta_lower)) / se
  t_upper <- (psi_hat -   delta_upper)  / se
  p_lower <- stats::pt(t_lower, df, lower.tail = FALSE)
  p_upper <- stats::pt(t_upper, df, lower.tail = TRUE)
  p_tost  <- max(p_lower, p_upper)

  # The noninferiority test is exactly the lower TOST component, so it is
  # reported rather than recomputed.
  p_noninferiority <- p_lower

  t_crit <- stats::qt(1 - alpha_level, df)
  ci_lo  <- psi_hat - t_crit * se
  ci_hi  <- psi_hat + t_crit * se

  equivalent  <- as.integer(ci_lo > -delta_lower && ci_hi < delta_upper)
  noninferior <- as.integer(ci_lo > -delta_lower)
  superior    <- as.integer(ci_lo >  delta_upper)
  inferior    <- as.integer(ci_hi < -delta_lower)

  # The order of the checks matters: superiority implies noninferiority,
  # so the more informative label wins.
  verdict <-
    if (equivalent == 1L)  "Equivalent"
    else if (superior == 1L) "Superior"
    else if (inferior == 1L) "Inferior"
    else if (noninferior == 1L) "Noninferior only"
    else "Inconclusive"

  out <- data.frame(
    term  = c("psi_hat", "se", "df",
              "t_lower", "t_upper",
              "p_lower", "p_upper", "p_tost", "p_noninferiority",
              "lower_limit", "upper_limit",
              "delta_lower", "delta_upper",
              "equivalent", "noninferior", "superior", "inferior"),
    value = c(psi_hat, se, df,
              t_lower, t_upper,
              p_lower, p_upper, p_tost, p_noninferiority,
              ci_lo, ci_hi,
              -delta_lower, delta_upper,
              equivalent, noninferior, superior, inferior),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  out <- .as_dmar_tbl(out, conf_level = 1 - 2 * alpha_level,
                      p_terms = c("p_lower", "p_upper", "p_tost",
                                  "p_noninferiority"))
  attr(out, "verdict") <- verdict
  out
}
