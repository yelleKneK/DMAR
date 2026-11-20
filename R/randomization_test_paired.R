# Paired-samples randomization (sign-flip) test.
#' Paired-Samples Randomization (Sign-Flip) Test
#'
#' Computes an exact (or Monte Carlo) sign-flip randomization test for
#' paired observations \eqn{(x_i, y_i)}, treating the within-pair sign
#' of \eqn{d_i = y_i - x_i} as the randomization mechanism (Fisher,
#' 1971; Edgington & Onghena, 2007). Under the null hypothesis of
#' exchangeability (\eqn{H_0}: the labeling of \eqn{x} and \eqn{y}
#' within each pair is arbitrary), each of the \eqn{2^n} sign patterns
#' is equally likely.
#'
#' For small \eqn{n} (default \eqn{n \le 20}) the test enumerates all
#' \eqn{2^n} sign patterns exactly; for larger \eqn{n} a Monte Carlo
#' approximation is used (default \code{n_resamples = 10000L}).
#'
#' @param x,y Paired numeric vectors of equal length. \code{NA}s are
#'   removed pairwise; the remaining values must be finite.
#' @param statistic One of \code{"mean"} (default) or \code{"t"}.
#'   \code{"mean"} reports the mean difference and \code{"t"} the paired
#'   \emph{t} statistic. Over the sign patterns the paired \emph{t} is a
#'   strictly increasing function of the mean difference, so the two give
#'   the same \emph{p}-value and differ only in the statistic reported.
#'   \code{"t"} is refused when every difference has the same nonzero
#'   value (exactly, or to within the tolerance of the comparison), because
#'   the statistic then has a zero standard error and no value (see
#'   \emph{When the studentized statistic is undefined}); when every
#'   difference is zero it is reported as 0 with a \emph{p}-value of 1.
#' @param alternative One of \code{"two_sided"} (default; the base-R
#'   spelling \code{"two.sided"} is accepted as an alias),
#'   \code{"less"}, or \code{"greater"}.
#' @param exact Logical. If \code{NULL} (default), uses exact
#'   enumeration when \eqn{n \le 20} and Monte Carlo otherwise. If
#'   \code{TRUE}, forces exact enumeration (caps at \eqn{n \le 25};
#'   above that the \eqn{2^n} space is too large). If \code{FALSE},
#'   forces Monte Carlo.
#' @param n_resamples A single whole number of at least 1: the number of
#'   Monte Carlo resamples when exact enumeration is not used. Default
#'   \code{10000L}.
#' @param seed Optional integer seed for reproducibility of the Monte
#'   Carlo branch. Default \code{NULL}, under which the randomizations
#'   are drawn from, and advance, the current random number stream;
#'   supply an integer for reproducibility, which is set for the call
#'   only, the caller's random number state being restored when it
#'   returns.
#'
#' @return A \code{data.frame} (class \code{dmar_tbl}) with columns
#'   \code{term} and \code{value} and five rows: \code{statistic} (the
#'   mean difference \eqn{\bar d}, or the paired \emph{t} under
#'   \code{statistic = "t"}, always finite), \code{p_value},
#'   \code{n_pairs} (the pairs with both scores present),
#'   \code{n_evaluated} (\eqn{2^n} under exact enumeration,
#'   \code{n_resamples} otherwise), and \code{exact} (1 when every sign
#'   pattern was enumerated, 0 when patterns were sampled). No table is
#'   returned when \code{statistic = "t"} is refused; the condition, of
#'   class \code{dmar_studentized_statistic_undefined}, carries
#'   \code{reason} (\code{"constant_differences"}), \code{alternative},
#'   \code{n} (the number of pairs), and \code{p_value_mean}, the exact
#'   \emph{p}-value \code{statistic = "mean"} gives for the same call.
#'
#' @details
#' \strong{Why randomization.} The randomization test makes no
#' distributional assumption on \eqn{d_i}; it only assumes that under
#' the null, the sign of each \eqn{d_i} is arbitrary. This is exactly
#' the inference that pre-experimental random assignment licenses.
#' Because the reference distribution is built from the observed
#' differences, the Type I error rate under that null does not exceed the
#' nominal level whatever the distribution of the differences, and for
#' any continuous distribution, heavy-tailed, a mixture, or with
#' outliers, it is the same attainable level: under exact enumeration,
#' the largest multiple of \eqn{1/2^{n-1}} (two-sided) or \eqn{1/2^n}
#' (one-sided) that does not exceed the nominal level (6/128 = .047
#' two-sided at .05 with 8 pairs). When different sign patterns can give
#' the same statistic, as with discrete scores, it can be lower.
#'
#' \strong{Exact enumeration.} All \eqn{2^n} sign patterns are
#' enumerated for \eqn{n} of 20 or fewer by default, and for \eqn{n} up
#' to 25 when \code{exact = TRUE}. The observed test statistic is compared
#' with the full reference distribution. The exact two-sided
#' \emph{p}-value is the proportion of patterns yielding a test
#' statistic at least as extreme (in absolute value) as the observed. A
#' pattern whose statistic equals the observed one counts as at least as
#' extreme; the comparison allows for the rounding of the arithmetic, with
#' a bound proportional to the size of the scores, so that patterns tied in
#' exact arithmetic are counted as tied.
#'
#' \strong{The two statistics order the patterns alike.} Write \eqn{S}
#' for the sum of the signed differences and \eqn{Q = \sum d_i^2}, which
#' no sign pattern changes. The paired \emph{t} is
#' \eqn{S \sqrt{n - 1} / \sqrt{n Q - S^2}}, which increases with \eqn{S}
#' because \eqn{Q} is fixed. Both statistics therefore rank the patterns
#' as \eqn{S} does, and the function counts on \eqn{S} for both.
#'
#' \strong{When the studentized statistic is undefined.} The studentized
#' statistic is a ratio, the mean difference over its standard error, and
#' it has a value only when that standard error is positive. If every
#' paired difference is the same nonzero number (exactly, or to within the
#' tolerance of the comparison), the standard error of the observed data
#' is 0 and \code{statistic = "t"} is refused with a condition of class
#' \code{dmar_studentized_statistic_undefined}. The observed data are
#' then at the boundary of the reference set: among the \eqn{2^n} sign
#' patterns only the observed one and its mirror image have a zero
#' standard error. The refusal concerns the statistic and not the
#' evidence. \code{statistic = "mean"} is defined for every arrangement
#' and gives the smallest \emph{p}-value the design allows in the
#' direction of the data, \eqn{2/2^n} two-sided and \eqn{1/2^n}
#' one-sided; the message states the value, and the condition carries it
#' as \code{p_value_mean}. The refusal is made before any sign pattern is
#' enumerated or drawn, so a refused Monte Carlo call leaves the random
#' number stream where it was.
#'
#' When the observed statistic is defined, some arrangements in the
#' reference set can still have a zero standard error: with paired
#' differences of equal size and mixed signs, the two patterns that give
#' them all one sign. Those patterns are the two most extreme, which is
#' where the mean difference also places them, and since the function
#' counts on the signed sum \eqn{S} for both statistics they are counted
#' without a standard error being formed. Refusing them would refuse
#' ordinary coarsely measured data. For differences of 2, 2, 2, \eqn{-2},
#' and 2 the paired \emph{t} is 1.5 and the \emph{p}-value, 12/32 =
#' 0.375, is the one the mean difference gives. When every difference is
#' zero, every sign pattern is the same data set; the statistic is
#' reported as 0 and the \emph{p}-value is 1.
#'
#' \strong{Monte Carlo branch.} For larger \eqn{n}, \code{n_resamples}
#' random sign patterns are drawn uniformly from \eqn{\{-1, +1\}^n};
#' the Monte Carlo \emph{p}-value uses the standard
#' (\eqn{1 + \mathrm{count}}) / (\eqn{1 + B}) plug-in to avoid
#' \emph{p} = 0.
#'
#' @references
#' Edgington, E. S., & Onghena, P. (2007). \emph{Randomization tests}
#'   (4th ed.). Chapman & Hall/CRC.
#'
#' Fisher, R. A. (1971). \emph{The design of experiments} (9th ed.,
#'   reprint). Hafner.
#'
#' Pitman, E. J. G. (1937). Significance tests which may be applied to
#'   samples from any populations. \emph{Supplement to the Journal of
#'   the Royal Statistical Society, 4}(1), 119--130.
#'
#' @seealso \code{\link[stats]{t.test}} (parametric paired test),
#'   \code{\link{probability_of_superiority_paired}}
#'
#' @examples
#' # 1. Small-n exact: Bayley scores on twin pairs.
#' control <- c(95, 102,  98, 107, 105)
#' treat   <- c(102, 108, 100, 112, 109)
#' randomization_test_paired(control, treat)
#'
#' # 2. Larger n: Monte Carlo branch.
#' set.seed(113)
#' x <- rnorm(50, 100, 15)
#' y <- x + rnorm(50,   5, 12)
#' randomization_test_paired(x, y, n_resamples = 10000L)
#'
#' # 3. Every pair differs by the same amount. The mean difference is
#' #    defined, and gives the smallest p-value 5 pairs allow. With
#' #    statistic = "t" the same call is refused, with a message that says
#' #    why: the observed standard error is 0, so the studentized statistic
#' #    is undefined.
#' before <- c(10, 12, 9, 15, 11)
#' randomization_test_paired(before, before + 2)
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @keywords htest
#'
#' @family hypothesis tests
#'
#' @export

randomization_test_paired <- function(x, y,
                                      statistic = c("mean", "t"),
                                      alternative = c("two_sided", "less", "greater"),
                                      exact = NULL,
                                      n_resamples = 10000L,
                                      seed = NULL) {
  statistic   <- match.arg(statistic)
  alternative <- .match_alternative(alternative)

  if (!is.numeric(x) || !is.numeric(y))
    .dmar_refuse_argument(
      "'x' and 'y' must be numeric vectors.", arg = c("x", "y"))
  if (length(x) != length(y))
    .dmar_refuse_argument(
      "'x' and 'y' must be the same length (paired observations).",
      arg = c("x", "y"))
  # Scores stored as integers are read as doubles: the rounding bounds
  # below sum absolute scores, and integer arithmetic overflows to NA once
  # a sum passes .Machine$integer.max.
  x <- as.numeric(x)
  y <- as.numeric(y)

  n_resamples <- .dmar_check_count(n_resamples, "n_resamples")

  ok <- !is.na(x) & !is.na(y)
  x <- x[ok]; y <- y[ok]
  if (!all(is.finite(x)) || !all(is.finite(y)))
    .dmar_refuse_argument("'x' and 'y' must be finite.", arg = c("x", "y"))
  d <- y - x
  n <- length(d)
  if (n < 3L)
    .dmar_refuse_argument("Need at least 3 paired observations.", arg = NULL)

  use_exact <- if (is.null(exact)) n <= 20L else isTRUE(exact)
  if (use_exact && n > 25L)
    .dmar_refuse_argument(
      "Exact enumeration requires n <= 25 (2^n grows too quickly above that).",
      arg = NULL)

  # Both statistics order the sign patterns as the signed sum S of the
  # differences does (the paired t is a strictly increasing function of
  # S), so the patterns are counted on S whichever statistic is reported.
  # No standard error enters the count, and a tie between two patterns is
  # judged against the rounding bound of a sum.
  tol <- .rt_paired_tol(x, y, d)

  # The paired t of differences that are all one nonzero value has a zero
  # standard error and no value. The call is refused here, before a sign
  # pattern is enumerated or drawn.
  tol_d <- 2 * .Machine$double.eps * max(abs(x) + abs(y))
  if (statistic == "t" && .rt_paired_constant(d, tol_d))
    .rt_refuse_constant_differences(d, alternative, use_exact,
                                    forced = identical(exact, FALSE) &&
                                      n <= 20L,
                                    n_resamples = n_resamples)

  if (use_exact) {
    # Every signed sum, by doubling: after pair i the vector holds the
    # 2^i sums of +/- d[1], ..., +/- d[i]. The last element is the pattern
    # with every sign positive, the observed one.
    S_null <- 0
    for (i in seq_len(n)) S_null <- c(S_null - d[i], S_null + d[i])
    n_eval  <- length(S_null)
    S_obs   <- S_null[n_eval]
    extreme <- .rt_extreme(S_null, S_obs, alternative, tol)
    p_value <- extreme / n_eval
  } else {
    .dmar_local_seed(seed)
    S_null <- numeric(n_resamples)
    for (j in seq_len(n_resamples)) {
      signs <- sample(c(-1L, 1L), n, replace = TRUE)
      S_null[j] <- sum(signs * d)
    }
    S_obs   <- sum(d)
    extreme <- .rt_extreme(S_null, S_obs, alternative, tol)
    p_value <- (1 + extreme) / (1 + n_resamples)
    n_eval  <- n_resamples
  }

  T_obs <- .rt_statistic(d, statistic, tol_d = tol_d)

  out <- data.frame(
    term  = c("statistic", "p_value", "n_pairs",
              "n_evaluated", "exact"),
    value = c(T_obs, p_value, n, n_eval, as.integer(use_exact)),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  .as_dmar_tbl(out, p_terms = "p_value")
}

# ---- internals ----
#
# The counting rule that turns a reference distribution into a p-value is
# shared with the independent-groups test and lives in
# R/randomization_test_internals.R as .rt_extreme(). Only the statistic
# reported for the observed differences is specific to the paired design.
#
# 'tol_d' is the rounding bound of one difference, 2 eps max(|x| + |y|).
# The paired t has two boundary cases that the bound decides. When every
# difference is within the bound of zero the statistic is 0. When the
# differences share a sign and their absolute values agree to within the
# bound, the standard deviation of the differences is zero in exact
# arithmetic and the statistic has no value: .rt_paired_constant() detects
# the case and randomization_test_paired() refuses statistic = "t" before
# .rt_statistic() is reached. Without the bound, differences such as
# 1.2 - 1.0 and 0.5 - 0.3, equal in exact arithmetic, would give a standard
# deviation of rounding error and a statistic of the order of 1e15.

.rt_paired_constant <- function(d, tol_d) {
  abs_d <- abs(d)
  max(abs_d) > tol_d && max(abs_d) - min(abs_d) <= tol_d &&
    (all(d > 0) || all(d < 0))
}

.rt_statistic <- function(d, statistic, tol_d) {
  if (statistic == "mean") return(mean(d))
  if (max(abs(d)) <= tol_d) return(0)
  mean(d) * sqrt(length(d)) / stats::sd(d)
}

# The refusal of statistic = "t" for constant nonzero differences. The
# message states why the statistic has no value, where the observed data
# lie in the reference set, and the exact p-value of the mean difference
# for the alternative of the call: of the 2^n sign patterns only the
# observed one reaches the observed sum in its own direction (1 / 2^n),
# only it and its mirror image reach it in absolute value (2 / 2^n), and
# every pattern lies at or beyond it in the other direction (1).

.rt_refuse_constant_differences <- function(d, alternative, use_exact,
                                            forced, n_resamples) {
  n <- length(d)
  n_patterns <- 2^n
  positive <- sum(d) > 0
  toward <- if (positive) "greater" else "less"
  numerator <- if (alternative == "two_sided") 2 else
    if (alternative == toward) 1 else n_patterns
  p_value_mean <- numerator / n_patterns
  value <- .rt_number_text(mean(d))
  reading <- if (numerator < n_patterns)
    sprintf("the smallest p-value this alternative can give with %d pairs",
            n)
  else
    sprintf(paste0("because every sign pattern has a mean difference at ",
                   "or %s the observed one"),
            if (positive) "below" else "above")
  monte_carlo <- if (!use_exact && numerator < n_patterns)
    .rt_monte_carlo_text(sprintf("%d pairs", n), forced, n_resamples)
  else ""
  # Differences that are equal only to within the tolerance of the
  # comparison (scores far larger than the spread of their differences)
  # are described as such, and the p-value quoted is then that of
  # differences that are all equal, which statistic = "mean" need not
  # return for the stored values. Differences that agree to the seven
  # digits a message prints (0.2 stored with a last bit of noise) are
  # equal for this purpose.
  within_rounding <- !identical(.rt_number_text(min(d)),
                                .rt_number_text(max(d)))
  equal_text <- if (within_rounding) {
    paste0(
      "All ", n, " paired differences equal ", value, " to within the ",
      "rounding of scores of this size (they span ",
      .rt_number_text(min(d)), " to ", .rt_number_text(max(d)), "), so ",
      "their standard deviation cannot be told from 0 and the paired t ",
      "statistic, the mean difference divided by its standard error, has ",
      "a denominator that cannot be told from 0. ")
  } else {
    paste0(
      "All ", n, " paired differences equal ", value, ", so ",
      "their standard deviation is 0 and the paired t statistic, the mean ",
      "difference divided by its standard error, has a zero denominator. ")
  }
  p_text <- paste0(.rt_p_text(numerator, n_patterns), " (", alternative,
                   "), ", reading)
  mean_text <- if (within_rounding) {
    paste0(
      "The mean difference is defined for every sign pattern. The p-value ",
      "of differences that are all equal is ", p_text, "; statistic = ",
      "\"mean\" uses the differences as stored and can return a different ",
      "p-value.")
  } else {
    paste0(
      "The mean difference is defined for every sign pattern, and ",
      "statistic = \"mean\" gives the mean difference ", value,
      " with the exact p-value ", p_text, ".")
  }
  .dmar_refuse(c(
    "randomization_test_paired(): statistic = \"t\" is not defined for ",
    "these data. ", equal_text,
    "The observed data lie at the boundary of the reference set: of the ",
    .rt_count_text(n_patterns), " sign patterns, the observed one and its ",
    "mirror image (every sign reversed) are the only two with a standard ",
    "error of 0, and every other pattern has a finite t. The statistic ",
    "therefore has no value for the data that were observed, and no ",
    "limiting value (Inf) is reported in its place. This is a statement ",
    "about the statistic and not about the evidence. ", mean_text,
    monte_carlo, " Use statistic = \"mean\"."),
    "studentized_statistic_undefined",
    reason = "constant_differences", alternative = alternative,
    p_value_mean = p_value_mean, n = n)
}
