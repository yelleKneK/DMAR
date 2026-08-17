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
#'   removed pairwise.
#' @param statistic One of \code{"mean"} (default) or \code{"t"}.
#'   \code{"mean"} uses the mean difference \eqn{\bar d} as the test
#'   statistic. \code{"t"} uses the paired \emph{t}-statistic
#'   \eqn{\bar d \sqrt{n} / s_d}, which is more robust to pair-to-pair
#'   variability in \eqn{|d_i|}.
#' @param alternative One of \code{"two_sided"} (default; the base-R
#'   spelling \code{"two.sided"} is accepted as an alias),
#'   \code{"less"}, or \code{"greater"}.
#' @param exact Logical. If \code{NULL} (default), uses exact
#'   enumeration when \eqn{n \le 20} and Monte Carlo otherwise. If
#'   \code{TRUE}, forces exact enumeration (caps at \eqn{n \le 25};
#'   above that the \eqn{2^n} space is too large). If \code{FALSE},
#'   forces Monte Carlo.
#' @param n_resamples Number of Monte Carlo resamples when exact
#'   enumeration is not used. Default \code{10000L}.
#' @param seed Optional integer seed for reproducibility of the Monte-
#'   Carlo branch. Default \code{NULL}, which leaves the user's current RNG state intact; supply an integer for reproducibility.
#'
#' @return A \code{data.frame} with rows for the observed test
#'   statistic, the \emph{p}-value, the number of pairs, the number of
#'   randomizations evaluated, and a flag indicating whether the test
#'   was exact or Monte Carlo.
#'
#' @details
#' \strong{Why randomization.} The randomization test makes no
#' distributional assumption on \eqn{d_i}; it only assumes that under
#' the null, the sign of each \eqn{d_i} is arbitrary. This is exactly
#' the inference that pre-experimental random assignment licenses, and
#' it is robust to heavy-tailed differences, mixtures, and outliers.
#'
#' \strong{Exact enumeration.} For \eqn{n \le 25}, all \eqn{2^n} sign
#' patterns are enumerated. The observed test statistic is compared
#' with the full reference distribution. The exact two-sided
#' \emph{p}-value is the proportion of patterns yielding a test
#' statistic at least as extreme (in absolute value) as the observed.
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
    stop("'x' and 'y' must be numeric vectors.")
  if (length(x) != length(y))
    stop("'x' and 'y' must be the same length (paired observations).")

  ok <- !is.na(x) & !is.na(y)
  x <- x[ok]; y <- y[ok]
  d <- y - x
  n <- length(d)
  if (n < 3L)
    stop("Need at least 3 paired observations.")

  use_exact <- if (is.null(exact)) n <= 20L else isTRUE(exact)
  if (use_exact && n > 25L)
    stop("Exact enumeration requires n <= 25 (2^n grows too quickly above that).")

  T_obs <- .rt_statistic(d, statistic)

  if (use_exact) {
    # Enumerate all 2^n sign vectors as binary expansions of 0 ... 2^n - 1.
    n_perm <- 2L^n
    T_null <- numeric(n_perm)
    bits   <- 2^seq.int(0L, n - 1L)
    for (j in seq_len(n_perm)) {
      signs <- (bitwAnd(j - 1L, bits) > 0L) * 2L - 1L  # -1 / +1
      T_null[j] <- .rt_statistic(signs * d, statistic)
    }
    extreme <- .rt_extreme(T_null, T_obs, alternative)
    p_value <- extreme / n_perm
    n_eval  <- n_perm
  } else {
    if (!is.null(seed)) {
      if (exists(".Random.seed", envir = .GlobalEnv)) {
        .old_seed <- get(".Random.seed", envir = .GlobalEnv)
        on.exit(assign(".Random.seed", .old_seed, envir = .GlobalEnv), add = TRUE)
      } else {
        on.exit(if (exists(".Random.seed", envir = .GlobalEnv)) rm(list = ".Random.seed", envir = .GlobalEnv), add = TRUE)
      }
      set.seed(seed)
    }
    T_null <- numeric(n_resamples)
    for (j in seq_len(n_resamples)) {
      signs <- sample(c(-1L, 1L), n, replace = TRUE)
      T_null[j] <- .rt_statistic(signs * d, statistic)
    }
    extreme <- .rt_extreme(T_null, T_obs, alternative)
    p_value <- (1 + extreme) / (1 + n_resamples)
    n_eval  <- n_resamples
  }

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
# R/randomization_test_internals.R as .rt_extreme(). Only the sign-flip
# statistic is specific to the paired design.

.rt_statistic <- function(d, statistic) {
  switch(statistic,
    mean = mean(d),
    t    = {
      m <- mean(d); s <- stats::sd(d); n <- length(d)
      if (s == 0) 0 else m * sqrt(n) / s
    }
  )
}
