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
#' @param n_perm Number of permutations for the significance test.
#'   Default \code{10000}. Set to 0 to skip the test.
#'
#' @return A \code{data.frame} (class \code{dmar_tbl}) in \code{term} /
#'   \code{value} layout with the row \code{tucker_phi}, the point
#'   estimate of \eqn{\phi}. When \code{n_perm > 0} the table also
#'   carries \code{p_value_perm}, the two-sided permutation
#'   \emph{p}-value, and \code{n_perm}, the number of permutations
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
#' loading patterns are unrelated, randomly permuting one of the
#' loading vectors and recomputing \eqn{\phi} produces a sampling
#' distribution against which the observed \eqn{\phi} can be evaluated.
#' The two-sided \emph{p}-value is \eqn{(r + 1) / (m + 1)}, where
#' \emph{r} counts the permuted \eqn{|\phi|} values at least as large as
#' the observed \eqn{|\phi|} and \emph{m} is \code{n_perm}. Adding one
#' to each part counts the observed arrangement, which is itself a
#' legitimate permutation; without it a \emph{p}-value of exactly zero
#' could be reported, a value a sampled permutation test cannot support
#' (Phipson & Smyth, 2010). The smallest reportable \emph{p}-value is
#' therefore \eqn{1 / (m + 1)}.
#'
#' \strong{Interpretation.} Benchmark values for \eqn{\phi} have been
#' proposed in the literature (Lorenzo-Seva & ten Berge, 2006), but
#' context always matters; this package reports the coefficient with its
#' uncertainty and leaves interpretation to the context of the
#' application.
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
#' @seealso \code{\link[stats]{cor}}, \code{\link{reliability_H}}
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
#' procrustes_phi(l1, l3, n_perm = 5000)
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @keywords htest
#'
#' @family multivariate and latent variable methods
#'
#' @export

procrustes_phi <- function(loadings_1, loadings_2, n_perm = 10000L) {
  if (!is.numeric(loadings_1) || !is.numeric(loadings_2))
    stop("'loadings_1' and 'loadings_2' must both be numeric vectors.")
  if (length(loadings_1) != length(loadings_2))
    stop("'loadings_1' and 'loadings_2' must be the same length.")
  if (length(loadings_1) < 2L)
    stop("Need at least 2 loadings to compute Tucker's phi.")

  phi <- sum(loadings_1 * loadings_2) /
         sqrt(sum(loadings_1^2) * sum(loadings_2^2))

  out <- data.frame(
    term  = "tucker_phi",
    value = phi,
    stringsAsFactors = FALSE,
    row.names = NULL
  )

  if (n_perm > 0L) {
    perm_phis <- replicate(n_perm, {
      l2_perm <- sample(loadings_2)
      sum(loadings_1 * l2_perm) /
        sqrt(sum(loadings_1^2) * sum(l2_perm^2))
    })
    n_extreme <- sum(abs(perm_phis) >= abs(phi))
    p_two_sided <- (1 + n_extreme) / (1 + n_perm)
    out <- rbind(out, data.frame(
      term  = c("p_value_perm", "n_perm"),
      value = c(p_two_sided, n_perm),
      stringsAsFactors = FALSE,
      row.names = NULL
    ))
  }

  .as_dmar_tbl(out, p_terms = "p_value_perm")
}
