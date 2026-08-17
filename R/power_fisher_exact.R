# Power of Fisher's exact test under Fisher's noncentral hypergeometric.
#' Power of Fisher's Exact Test (Noncentral Hypergeometric)
#'
#' Computes the power of Fisher's exact test (Fisher, 1934) for the
#' 2 \eqn{\times} 2 table under Fisher's noncentral hypergeometric
#' distribution, where the alternative is parameterized by the true odds
#' ratio \eqn{\psi}. The power is the probability that the
#' (conditional) exact test rejects \eqn{H_0: \psi = 1} when in fact
#' \eqn{\psi = \psi_1 \ne 1}.
#'
#' @param n_1,n_2 Group sample sizes for the two columns of the
#'   2 \eqn{\times} 2 table (e.g., treatment vs control).
#' @param p_1,p_2 Success probabilities in the two groups under the
#'   alternative. The odds ratio under the alternative is
#'   \eqn{\psi = (p_1 / (1 - p_1)) / (p_2 / (1 - p_2))}.
#' @param alpha_level Two-sided significance level. Default \code{0.05}.
#' @param alternative One of \code{"two_sided"} (default; the base-R
#'   spelling \code{"two.sided"} is accepted as an alias),
#'   \code{"less"}, or \code{"greater"}.
#'
#' @return A \code{data.frame} with rows for the power, the
#'   alternative-odds-ratio \eqn{\psi_1}, the alternative success
#'   probabilities \eqn{p_1} and \eqn{p_2}, the expected column-1
#'   total \eqn{E[S]} across both groups, and the row / column totals.
#'
#' @details
#' \strong{Setup.} Fisher's exact test conditions on the marginal
#' totals of the 2 \eqn{\times} 2 table:
#' \tabular{lccc}{
#'   \tab Success \tab Failure \tab Total \cr
#'   Group 1 \tab \eqn{X_1} \tab \eqn{n_1 - X_1} \tab \eqn{n_1} \cr
#'   Group 2 \tab \eqn{S - X_1} \tab \eqn{(n_1 + n_2) - n_1 - S + X_1} \tab \eqn{n_2} \cr
#'   Total   \tab \eqn{S} \tab \eqn{n_1 + n_2 - S} \tab \eqn{n_1 + n_2}
#' }
#'
#' Under \eqn{H_0: \psi = 1}, \eqn{X_1} given the marginals follows the
#' \emph{central} hypergeometric. Under the alternative \eqn{\psi_1},
#' \eqn{X_1} follows \emph{Fisher's noncentral} hypergeometric with odds
#' ratio parameter \eqn{\psi_1} (Fisher, 1935; Fog, 2008), the
#' conditional distribution of one binomial count given the total of two
#' independent binomials. (Wallenius' noncentral hypergeometric, which
#' arises from sequential biased urn sampling, is a different
#' distribution and is not the relevant one here.)
#'
#' \strong{Power calculation.} For each possible value of the column-1
#' total \eqn{S = 0, 1, \ldots, n_1 + n_2}:
#' \enumerate{
#'   \item Determine the rejection region under \eqn{H_0: \psi = 1}
#'         using the central hypergeometric.
#'   \item Compute \eqn{\Pr(X_1 \in \mathrm{reject} \mid \psi = \psi_1, S)}
#'         under the noncentral hypergeometric.
#'   \item Weight by \eqn{\Pr(S \mid \psi_1)}, the marginal probability
#'         of total column-1 successes under the alternative.
#' }
#' The power is the resulting weighted sum.
#'
#' @references
#' Fisher, R. A. (1934). \emph{Statistical methods for research
#'   workers} (5th ed.). Oliver & Boyd.
#'
#' Fisher, R. A. (1935). The logic of inductive inference.
#'   \emph{Journal of the Royal Statistical Society, 98}(1), 39--82.
#'
#' Fog, A. (2008). Sampling methods for Wallenius' and Fisher's
#'   noncentral hypergeometric distributions. \emph{Communications in
#'   Statistics -- Simulation and Computation, 37}(2), 241--257.
#'   \doi{10.1080/03610910701790236}
#'
#' Good, P. I. (2000). \emph{Permutation tests: A practical guide to
#'   resampling methods for testing hypotheses} (2nd ed.). Springer.
#'
#' O'Brien, R. G. (1998). A tour of UnifyPow: A SAS module/macro for
#'   sample size analysis. \emph{Proceedings of the 23rd SAS Users
#'   Group International Conference}, 1346--1355.
#'
#' @seealso \code{\link[stats]{fisher.test}}
#'
#' @examples
#' # 1. Power for n_1 = n_2 = 30 when the two population proportions are
#' #    0.6 and 0.3. The value comes from enumerating the conditional
#' #    reference set the test itself uses, not from a normal
#' #    approximation, so it is the power of the test as conducted.
#' power_fisher_exact(n_1 = 30, n_2 = 30, p_1 = 0.6, p_2 = 0.3)
#'
#' # 2. A difference of 0.10 between proportions is much harder to detect:
#' #    100 per group is not close to enough. Sample size requirements grow
#' #    quickly as the difference between the two proportions shrinks.
#' power_fisher_exact(n_1 = 100, n_2 = 100, p_1 = 0.45, p_2 = 0.35)
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @keywords htest
#'
#' @family sample size for power
#'
#' @export

power_fisher_exact <- function(n_1, n_2, p_1, p_2,
                               alpha_level = 0.05,
                               alternative = c("two_sided", "less", "greater")) {
  alternative <- .match_alternative(alternative)

  for (nm in c("n_1", "n_2", "p_1", "p_2", "alpha_level")) {
    v <- get(nm)
    if (!is.numeric(v) || length(v) != 1L)
      stop(sprintf("'%s' must be a single numeric value.", nm))
  }
  if (n_1 < 2 || n_2 < 2)
    stop("'n_1' and 'n_2' must each be at least 2.")
  if (p_1 <= 0 || p_1 >= 1 || p_2 <= 0 || p_2 >= 1)
    stop("'p_1' and 'p_2' must be in (0, 1).")
  if (alpha_level <= 0 || alpha_level >= 1)
    stop("'alpha_level' must be in (0, 1).")
  # Odds ratio under the alternative:
  psi_1 <- (p_1 / (1 - p_1)) / (p_2 / (1 - p_2))
  N <- n_1 + n_2

  # Expected column-1 total under the alternative:
  e_s  <- n_1 * p_1 + n_2 * p_2
  power_sum <- 0

  for (s in 0:N) {
    # Skip empty conditional tables:
    x_min <- max(0L, s - n_2)
    x_max <- min(n_1, s)
    if (x_max < x_min) next

    # Pr(column-1 total = s) under the alternative:
    p_s <- stats::dbinom(0:n_1, n_1, p_1) %*%
           stats::dbinom(s - 0:n_1, n_2, p_2)
    p_s <- as.numeric(p_s)
    if (p_s < .Machine$double.eps) next

    # Reject region under H_0: central hypergeom on x in [x_min, x_max].
    x_seq <- x_min:x_max
    p0    <- stats::dhyper(x_seq, n_1, n_2, s)

    reject <- switch(alternative,
      two_sided = {
        # Two-sided: order by p0 ascending; cumulate until <= alpha_level.
        ord <- order(p0)
        cum <- cumsum(p0[ord])
        in_R <- logical(length(x_seq))
        in_R[ord[cum <= alpha_level]] <- TRUE
        in_R
      },
      less    = stats::phyper(x_seq, n_1, n_2, s)              <= alpha_level,
      greater = stats::phyper(x_seq - 1L, n_1, n_2, s, lower.tail = FALSE) <= alpha_level
    )
    if (!any(reject)) next

    # Pr(X_1 in reject | psi = psi_1, s) under Fisher's noncentral
    # hypergeometric.
    p1_seq <- .dfnc_hypergeo(x_seq, n_1, n_2, s, psi_1)
    power_sum <- power_sum + p_s * sum(p1_seq[reject])
  }

  out <- data.frame(
    term  = c("power", "odds_ratio_alt", "p_1", "p_2",
              "expected_s", "n_1", "n_2", "alpha_level"),
    value = c(power_sum, psi_1, p_1, p_2, e_s, n_1, n_2, alpha_level),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  .as_dmar_tbl(out)
}

# Fisher's noncentral hypergeometric density, computed directly from its
# definition: with x successes among n_1 draws when s successes are split
# between margins of size n_1 and n_2 at odds ratio psi,
#   Pr(X = x) is proportional to choose(n_1, x) choose(n_2, s - x) psi^x
# over the support max(0, s - n_2) <= x <= min(n_1, s). The finite sum is
# normalized on the log scale, so no external package and no asymptotic
# approximation is involved. Replaces BiasedUrn::dFNCHypergeo(), which
# computes the same quantity; equivalence is pinned in the tests and
# checked live in tools/oracle_checks.R.
.dfnc_hypergeo <- function(x, m1, m2, s, odds) {
  lo <- max(0L, s - m2)
  hi <- min(m1, s)
  support <- lo:hi
  lw <- lchoose(m1, support) + lchoose(m2, s - support) +
    support * log(odds)
  lw <- lw - max(lw)
  w  <- exp(lw)
  w  <- w / sum(w)
  out <- numeric(length(x))
  ok  <- x >= lo & x <= hi
  out[ok] <- w[x[ok] - lo + 1L]
  out
}
