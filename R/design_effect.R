# Kish's design effect and its square root (DEFT).
#' Kish's Design Effect (DEFF), DEFT, and the Effective Sample Size
#'
#' Computes the design effect (DEFF) and its square root (DEFT) for a
#' clustered or multistage sample, given a vector of per-cluster sample
#' sizes and a value of the intraclass correlation. Returns Kish's
#' (1965) classic formula together with the effective sample size and a
#' description of the clustering, including the number of empty
#' clusters (no observations) and the number of singleton clusters (one
#' observation), which carry different amounts of within-cluster
#' information in a mixed-effects context.
#'
#' \strong{Definition (Kish, 1965).} For a clustered sample with
#' per-cluster sizes \eqn{m_1, m_2, \ldots, m_K} and intraclass correlation
#' \eqn{\rho}, the design effect on the variance of the mean is
#' \deqn{\mathrm{DEFF} \;=\; 1 + (\bar m^{*} - 1) \rho,}
#' where \eqn{\bar m^{*} = \sum_k m_k^2 / \sum_k m_k} is the
#' design-weighted average cluster size (Kish, 1965, eq. 5.4; sometimes called
#' the "Kish weighted average" or "effective cluster size"). For equal
#' cluster sizes \eqn{m_k = m}, this reduces to the classroom form
#' \eqn{1 + (m - 1)\rho}. The DEFT is the square root of DEFF and is
#' the inflation factor on the \emph{standard error} of the mean
#' (whereas DEFF inflates the variance).
#'
#' \strong{Effective sample size.} The number of observations from a
#' simple random sample that would yield the same standard error as the
#' clustered sample is
#' \deqn{N_{\mathrm{eff}} \;=\; N / \mathrm{DEFF} \;=\; N / \mathrm{DEFT}^2,}
#' where \eqn{N = \sum_k m_k} is the total observations.
#'
#' \strong{Why empty and singleton clusters are reported separately.}
#' In mixed-effects / multilevel modeling, clusters with zero
#' observations carry no information (they should be dropped before
#' fitting), and clusters with one observation contribute to \eqn{N}
#' and to the fixed-effect estimate but contribute nothing to the
#' estimation of the random-effect variance or to the within-cluster
#' residual. Hox et al. (2017) note that singleton-heavy designs have a
#' design effect close to 1 even at moderate \eqn{\rho} because the
#' weighted cluster size is small. The output reports the counts of
#' empty and singleton clusters so the user can see at a glance how
#' much of the nominal sample size carries clustering information.
#'
#' @param cluster_sizes Numeric vector of per-cluster sample sizes. Each
#'   element is the number of observations in one cluster (so the
#'   length of the vector is the number of clusters). Zero values are
#'   allowed and counted as empty clusters; they do not affect the
#'   computation of DEFF.
#' @param icc Intraclass correlation coefficient, in \eqn{[0, 1)}. Use
#'   the population value when planning prospectively, or the sample
#'   estimate (e.g., from \code{\link{icc}} or \code{\link{icc_lmer}})
#'   when describing an observed clustered sample.
#'
#' @return A \code{data.frame} with rows for the design effect,
#'   its square root, the effective sample size, the total observation
#'   and cluster counts (including separate counts of empty and
#'   singleton clusters), the mean cluster size, Kish's design-weighted
#'   mean cluster size, and the input \code{icc}.
#'
#' @details
#' If a fitted \code{lmerMod} object is available, the typical workflow
#' is to extract \code{icc} via \code{\link{icc_lmer}} and the
#' per-cluster sample sizes via \code{table(cluster_id)}, then pass both
#' to \code{design_effect()}; see the second example.
#'
#' @references
#' Hox, J. J., Moerbeek, M., & van de Schoot, R. (2017). \emph{Multilevel
#'   analysis: Techniques and applications} (3rd ed.). Routledge.
#'
#' Kish, L. (1965). \emph{Survey sampling}. Wiley.
#'
#' Kish, L. (1992). Weighting for unequal Pi. \emph{Journal of Official
#'   Statistics, 8}(2), 183--200.
#'
#' Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). \emph{Designing
#'   experiments and analyzing data: A model comparison perspective}
#'   (4th ed.). Routledge. (See Chapters 15 and 16 on mixed-effects
#'   models and nested designs.)
#'
#' Snijders, T. A. B., & Bosker, R. J. (2012). \emph{Multilevel analysis:
#'   An introduction to basic and advanced multilevel modeling} (2nd
#'   ed.). Sage.
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso \code{\link{icc}}, \code{\link{icc_lmer}}, \code{\link{var_icc}},
#'   \code{\link{ss_aipe_icc}}
#'
#' @examples
#' # 1. Balanced design: K = 30 clusters of size 20, ICC = 0.10.
#' design_effect(cluster_sizes = rep(20, 30), icc = 0.10)
#' # DEFF = 1 + (20 - 1) * 0.10 = 2.9; DEFT = 1.7; effective N = 600 / 2.9 = 207.
#'
#' # 2. Unbalanced design with some empty and some singleton clusters.
#' # Suppose K = 25 schools with attendance ranging from 0 (closed)
#' # through 1 (single student showed up) to 30 (full class).
#' sizes <- c(0, 0, 1, 1, 1, 2, 3, 5, 8, 10, 12, 15, 18, 20, 22,
#'            22, 25, 26, 28, 28, 30, 30, 30, 30, 30)
#' design_effect(cluster_sizes = sizes, icc = 0.10)
#'
#' # 3. From a cluster-id vector: tabulate, then call design_effect().
#' cluster_id <- rep(1:8, times = c(20, 15, 22, 1, 0, 30, 18, 25))
#' design_effect(cluster_sizes = as.numeric(table(factor(cluster_id, levels = 1:8))),
#'      icc = 0.15)
#'
#' @keywords design multivariate
#'
#' @family design utilities
#'
#' @export
design_effect <- function(cluster_sizes, icc) {
  if (!is.numeric(cluster_sizes) || length(cluster_sizes) < 1L) {
    stop("'cluster_sizes' must be a non-empty numeric vector.", call. = FALSE)
  }
  if (any(is.na(cluster_sizes))) {
    stop("'cluster_sizes' must not contain NA values.", call. = FALSE)
  }
  if (any(cluster_sizes < 0) || any(cluster_sizes != as.integer(cluster_sizes))) {
    stop("'cluster_sizes' values must be non-negative integers.", call. = FALSE)
  }
  if (!is.numeric(icc) || length(icc) != 1L || is.na(icc) || icc < 0 || icc >= 1) {
    stop("'icc' must be a single number in [0, 1).", call. = FALSE)
  }

  m <- as.numeric(cluster_sizes)

  n_clusters_total      <- length(m)
  n_clusters_empty      <- sum(m == 0)
  n_clusters_singletons <- sum(m == 1)
  n_clusters_with_data  <- sum(m >  0)
  n_clusters_informative <- sum(m >= 2)
  n_total               <- sum(m)

  if (n_total == 0) {
    stop("All clusters are empty: total sample size is 0.", call. = FALSE)
  }

  # Mean cluster size among non-empty clusters (the natural denominator
  # for descriptive interpretation; empty clusters carry no information).
  m_bar <- n_total / n_clusters_with_data

  # Kish's design-weighted average cluster size:
  #   m_kish = sum(m^2) / sum(m)
  # Empty clusters contribute 0 / 0 to neither numerator nor denominator
  # and so drop out cleanly. Singletons contribute 1 to both, pulling
  # m_kish toward 1 when they are numerous (and so the design effect
  # toward 1, as expected: singletons carry only fixed-effect information).
  m_kish <- sum(m^2) / n_total

  deff_val      <- 1 + (m_kish - 1) * icc
  deft_val      <- sqrt(deff_val)
  effective_n   <- n_total / deff_val

  out <- data.frame(
    term  = c("design_effect", "deft", "effective_n",
              "n_total",
              "n_clusters_total", "n_clusters_with_data",
              "n_clusters_empty", "n_clusters_singletons",
              "n_clusters_informative",
              "m_bar", "m_kish",
              "icc"),
    value = c(deff_val, deft_val, effective_n,
              n_total,
              n_clusters_total, n_clusters_with_data,
              n_clusters_empty, n_clusters_singletons,
              n_clusters_informative,
              m_bar, m_kish,
              icc),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  .as_dmar_tbl(out)
}

