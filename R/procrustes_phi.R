# Tucker's congruence coefficient with a permutation test.
#' Tucker's Congruence Coefficient \eqn{\phi} (Factor Similarity)
#'
#' Computes Tucker's (1951) congruence coefficient \eqn{\phi}, a measure
#' of similarity between two factor-loading patterns (typically the
#' standardized loadings of the same factor estimated on two different
#' samples or with different methods), together with a permutation-based
#' \emph{p}-value testing the null hypothesis of unrelated loading
#' patterns. \eqn{\phi} is the standard tool for factor-replication
#' studies (Lorenzo-Seva & ten Berge, 2006).
#'
#' @param loadings_1,loadings_2 Numeric vectors of factor loadings on
#'   the same indicator set, of equal length. Either standardized or
#'   raw loadings work; the coefficient is scale-invariant.
#' @param n_resamples A single whole number of at least 0: the number of
#'   permutations for the significance test. Default \code{10000}. Set to
#'   0 to skip the test.
#' @param seed \code{NULL} (the default) to draw from the current random
#'   number stream, or a single number that makes the call reproducible; a
#'   supplied seed is set for this call only, and the caller's random
#'   number state is restored when the function returns.
#'
#' @return A \code{data.frame} (class \code{dmar_tbl}) in \code{term} /
#'   \code{value} layout with the row \code{tucker_phi}, the point
#'   estimate of \eqn{\phi}. When \code{n_resamples > 0} the table also
#'   carries \code{p_value_perm}, the two-sided permutation
#'   \emph{p}-value, and \code{n_resamples}, the number of permutations
#'   requested.
#'
#' @details
#' \strong{Definition.} For two vectors of loadings
#' \eqn{\bm\lambda_1, \bm\lambda_2} on a shared set of \eqn{p} indicators,
#' Tucker's congruence coefficient is
#' \deqn{\phi(\bm\lambda_1, \bm\lambda_2) \;=\;
#'   \frac{\sum_{i=1}^{p} \lambda_{1i} \lambda_{2i}}
#'        {\sqrt{\sum_{i=1}^{p} \lambda_{1i}^2 \cdot
#'               \sum_{i=1}^{p} \lambda_{2i}^2}}.}
#' \eqn{\phi} is the cosine of the angle between the two loading
#' vectors and ranges over \eqn{[-1, 1]}; values near \eqn{\pm 1}
#' indicate high (anti-)congruence, values near 0 indicate orthogonality.
#'
#' \strong{Permutation test.} Under the null hypothesis that the two
#' loading patterns are unrelated, in the sense that the pairing of the
#' items' loadings across the two vectors is arbitrary, permuting one
#' vector and recomputing \eqn{\phi} gives a permutation distribution
#' against which the observed \eqn{\phi} can be evaluated. The test takes
#' the two sets of loadings as given: it does not reflect the sampling
#' variability of loadings estimated from data, and the function reports
#' no standard error or interval for \eqn{\phi}.
#' The two-sided \emph{p}-value is \eqn{(r + 1) / (m + 1)}, where
#' \emph{r} counts the permuted \eqn{|\phi|} values at least as large as
#' the observed \eqn{|\phi|} and \emph{m} is \code{n_resamples}. A permutation
#' whose coefficient equals the observed one counts as at least as
#' extreme, as when it exchanges two items with the same loading; the
#' comparison allows for the rounding of the arithmetic. Adding one
#' to each part counts the observed arrangement, which is itself a
#' legitimate permutation; without it a \emph{p}-value of exactly zero
#' could be reported, a value a sampled permutation test cannot support
#' (Phipson & Smyth, 2010). The smallest reportable \emph{p}-value is
#' therefore \eqn{1 / (m + 1)}.
#'
#' \strong{Interpretation.} Benchmark values for \eqn{\phi} have been
#' proposed in the literature (Lorenzo-Seva & ten Berge, 2006), but
#' context always matters; this package reports the coefficient and the
#' permutation \emph{p}-value and leaves interpretation to the context of
#' the application.
#'
#' @references
#' Lorenzo-Seva, U., & ten Berge, J. M. F. (2006). Tucker's congruence
#'   coefficient as a meaningful index of factor similarity.
#'   \emph{Methodology, 2}(2), 57--64. \doi{10.1027/1614-2241.2.2.57}
#'
#' Phipson, B., & Smyth, G. K. (2010). Permutation \emph{p}-values should
#'   never be zero: Calculating exact \emph{p}-values when permutations are
#'   randomly drawn. \emph{Statistical Applications in Genetics and
#'   Molecular Biology, 9}(1), Article 39. \doi{10.2202/1544-6115.1585}
#'
#' Tucker, L. R. (1951). \emph{A method for synthesis of factor analysis
#'   studies} (Personnel Research Section Report No. 984). Department
#'   of the Army.
#'
#' @seealso \code{\link{factor_congruence}} for the coefficients between
#'   the factors of two whole solutions, matched one to one and reflected,
#'   \code{\link[stats]{cor}}, \code{\link{reliability_H}}
#'
#' @examples
#' set.seed(113)
#' # 1. Two highly similar loading patterns:
#' l1 <- c(0.72, 0.65, 0.81, 0.55, 0.69)
#' l2 <- c(0.70, 0.62, 0.83, 0.58, 0.66)
#' procrustes_phi(l1, l2)
#'
#' # 2. Loadings on different factors should show low congruence:
#' l3 <- c(0.10, 0.05, 0.20, 0.85, 0.78)
#' procrustes_phi(l1, l3, n_resamples = 5000)
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @keywords htest
#'
#' @family multivariate and latent variable methods
#'
#' @export

procrustes_phi <- function(loadings_1, loadings_2, n_resamples = 10000L,
                           seed = NULL) {
  if (!is.numeric(loadings_1) || !is.numeric(loadings_2))
    .dmar_refuse_argument(
      "'loadings_1' and 'loadings_2' must both be numeric vectors.",
      arg = c("loadings_1", "loadings_2"))
  if (length(loadings_1) != length(loadings_2))
    .dmar_refuse_argument(
      "'loadings_1' and 'loadings_2' must be the same length.",
      arg = c("loadings_1", "loadings_2"))
  if (length(loadings_1) < 2L)
    .dmar_refuse_argument(
      "Need at least 2 loadings to compute Tucker's phi.", arg = NULL)

  n_resamples <- .dmar_check_count(n_resamples, "n_resamples", min = 0L)
  .dmar_local_seed(seed)

  n   <- length(loadings_1)
  den <- sqrt(sum(loadings_1^2) * sum(loadings_2^2))
  phi <- sum(loadings_1 * loadings_2) / den

  out <- data.frame(
    term  = "tucker_phi",
    value = phi,
    stringsAsFactors = FALSE,
    row.names = NULL
  )

  if (n_resamples > 0L) {
    # A permutation leaves the denominator unchanged, so it is computed
    # once and every coefficient shares it.
    perm_phis <- replicate(n_resamples, {
      l2_perm <- sample(loadings_2)
      sum(loadings_1 * l2_perm) / den
    })
    # A permutation that exchanges tied loadings reproduces the observed
    # coefficient in exact arithmetic but can fall a unit in the last place
    # below it as computed, and a tie counts as at least as extreme. With
    # u = eps / 2, each product carries at most 3u |l_1 l_2| and the sum of
    # n products adds at most (n - 1) u sum|l_1 l_2|; since sum|l_1 l_2| is
    # at most the denominator, one coefficient is in error by at most
    # (n + 2) u, and the difference of two coefficients by at most
    # (n + 2) eps. The tolerance is four times that bound.
    tol <- 4 * (n + 2) * .Machine$double.eps
    n_extreme <- sum(abs(perm_phis) >= abs(phi) - tol)
    p_two_sided <- (1 + n_extreme) / (1 + n_resamples)
    out <- rbind(out, data.frame(
      term  = c("p_value_perm", "n_resamples"),
      value = c(p_two_sided, n_resamples),
      stringsAsFactors = FALSE,
      row.names = NULL
    ))
  }

  .as_dmar_tbl(out, p_terms = "p_value_perm")
}
