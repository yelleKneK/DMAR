# Krippendorff's alpha inter-rater agreement coefficient.
#' Krippendorff's \eqn{\alpha} Inter-Rater Agreement
#'
#' Computes Krippendorff's (1980, 2004, 2011) \eqn{\alpha}, the most
#' general chance-corrected inter-rater agreement coefficient. Unlike
#' \code{\link{cohen_kappa}} (two raters, nominal data) or
#' \code{\link{fleiss_kappa}} (multiple raters, nominal data),
#' Krippendorff's \eqn{\alpha} supports any number of raters,
#' missing values, and any of four levels of measurement (nominal,
#' ordinal, interval, ratio) via a user-specified distance metric.
#'
#' @param ratings A units \eqn{\times} raters matrix (or
#'   \code{data.frame}). Rows = units of analysis; columns = raters.
#'   \code{NA} entries are allowed.
#' @param level One of \code{"nominal"} (default), \code{"ordinal"},
#'   \code{"interval"}, or \code{"ratio"}; controls the distance
#'   metric used to compute disagreement.
#' @param conf_level Confidence level for the bootstrap CI. Default
#'   \code{0.95}.
#' @param boot Logical. If \code{TRUE}, returns a bootstrap
#'   percentile CI. Set to \code{FALSE} to return only the point
#'   estimate (much faster).
#' @param B Number of bootstrap resamples when \code{boot = TRUE}.
#'   Default \code{1000L}.
#' @param seed Optional integer seed for reproducibility of the
#'   bootstrap. Default \code{NULL}, which leaves the user's current RNG state intact; supply an integer for reproducibility.
#'
#' @return A \code{data.frame} with rows for the point estimate
#'   \eqn{\hat\alpha}, the observed disagreement \eqn{D_o}, the
#'   expected disagreement \eqn{D_e}, the number of pairable values,
#'   and, when a bootstrap was run, the lower and upper bootstrap CI
#'   limits and \code{B_used}, the number of resamples that
#'   returned a finite value and so entered the interval.
#'
#' @details
#' \strong{Coefficient.} Krippendorff's \eqn{\alpha} is
#' \deqn{\alpha \;=\; 1 - \frac{D_o}{D_e},}
#' where \eqn{D_o} is the observed disagreement (average squared
#' distance over all within-unit pairs of ratings, scaled by the
#' number of pairable values), and \eqn{D_e} is the expected
#' disagreement (average squared distance over all between-unit
#' pairs). The metric used in the squared distance depends on
#' \code{level}:
#' \itemize{
#'   \item \emph{nominal}: \eqn{d(a, b) = \mathrm{I}(a \ne b)}
#'   \item \emph{ordinal}: distance based on cumulative rank counts
#'   \item \emph{interval}: \eqn{d(a, b) = (a - b)^2}
#'   \item \emph{ratio}: \eqn{d(a, b) = ((a - b) / (a + b))^2}
#' }
#'
#' \strong{CI.} The CI is by case-resampling bootstrap over units
#' (rows): the rows of \code{ratings} are resampled with replacement
#' \code{B} times and \eqn{\alpha} is recomputed on each resample,
#' so units are the sampling unit and the rater panel is treated as
#' fixed. Only the percentile interval is offered: the limits are the
#' empirical quantiles of the bootstrap estimates (Efron & Tibshirani,
#' 1993); there is no bias-corrected and accelerated (BCa) variant.
#' Resamples on which the coefficient cannot be computed (for example,
#' a resample without enough pairable values) are dropped; the interval
#' is computed from the ones that return a finite value, and how many
#' did is reported as the \code{B_used} row of the result. No
#' closed-form sampling variance is in general use for Krippendorff's
#' alpha across its measurement levels and missing data patterns, so the
#' bootstrap is the interval Krippendorff recommends (Krippendorff,
#' 2011; Hayes & Krippendorff, 2007). \code{B = 1000L} typically gives a
#' stable CI to two decimal places. The bootstrap is opt-in
#' (\code{boot = FALSE} by default, which returns the point estimate
#' alone and is much faster); ask for it whenever the coefficient is
#' being reported rather than explored, since a point estimate on its
#' own says nothing about how precisely \eqn{\alpha} is determined.
#' Bootstrap results vary from run to run; supply \code{seed} for
#' reproducibility.
#'
#' \strong{Interpretation.} \eqn{\alpha} ranges from \eqn{-D_e / D_o}
#' (perfect disagreement) through \eqn{0} (chance level) to \eqn{1}
#' (perfect agreement). Report the coefficient with its confidence
#' interval and judge it against the reliability the application
#' requires; Krippendorff (2004) discusses how that judgment depends on
#' the cost of acting on unreliable data.
#'
#' @references
#' Efron, B., & Tibshirani, R. J. (1993). \emph{An introduction to the
#'   bootstrap}. New York, NY: Chapman & Hall/CRC.
#'
#' Hayes, A. F., & Krippendorff, K. (2007). Answering the call for a
#'   standard reliability measure for coding data. \emph{Communication
#'   Methods and Measures, 1}(1), 77--89. \doi{10.1080/19312450709336664}
#'
#' Krippendorff, K. (1980). \emph{Content analysis: An introduction to
#'   its methodology}. Sage.
#'
#' Krippendorff, K. (2004). \emph{Content analysis: An introduction to
#'   its methodology} (2nd ed.). Sage.
#'
#' Krippendorff, K. (2011). Computing Krippendorff's
#'   alpha-reliability. \emph{Departmental Papers (ASC)}, Annenberg
#'   School for Communication, University of Pennsylvania.
#'
#' @seealso \code{\link{cohen_kappa}}, \code{\link{fleiss_kappa}},
#'   \code{\link{icc}}
#'
#' @examples
#' # 1. Nominal ratings, 4 raters, 12 units (from Krippendorff 2011 Tab. 1):
#' ratings <- matrix(c(
#'   1, 2, 3, 3, 2, 1, 4, 1, 2, NA, NA, NA,
#'   1, 2, 3, 3, 2, 2, 4, 1, 2, 5,  NA, 3,
#'   NA, 3, 3, 3, 2, 3, 4, 2, 2, 5,  1,  NA,
#'   1, 2, 3, 3, 2, 4, 4, 1, 2, 5,  1,  NA
#' ), nrow = 12, ncol = 4)
#' krippendorff_alpha(ratings, level = "nominal")
#'
#' # 2. Interval ratings:
#' set.seed(113)
#' r1 <- rnorm(30, 0, 1)
#' r2 <- r1 + rnorm(30, 0, 0.3)
#' krippendorff_alpha(cbind(r1, r2), level = "interval")
#'
#' # The percentile bootstrap interval for the same ratings, which
#' # recomputes alpha on each of B resamples of the units. Not run
#' # here, because 500 refits of alpha is more than a help page should
#' # do; the call is:
#' # krippendorff_alpha(cbind(r1, r2), level = "interval",
#' #                    boot = TRUE, B = 500L, seed = 113)
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @keywords htest
#'
#' @family agreement and measurement
#'
#' @export

krippendorff_alpha <- function(ratings,
                               level = c("nominal", "ordinal",
                                         "interval", "ratio"),
                               conf_level = 0.95,
                               boot = FALSE, B = 1000L,
                               seed = NULL) {
  level <- match.arg(level)
  if (is.data.frame(ratings)) ratings <- as.matrix(ratings)
  if (!is.matrix(ratings))
    stop("'ratings' must be a matrix or data.frame.")
  if (conf_level <= 0 || conf_level >= 1)
    stop("'conf_level' must be in (0, 1).")

  one_alpha <- function(R) {
    .kripp_alpha_core(R, level)
  }

  obs <- one_alpha(ratings)

  out <- data.frame(
    term  = c("krippendorff_alpha", "D_observed", "D_expected", "n_pairable"),
    value = c(obs$alpha, obs$D_o, obs$D_e, obs$n_pairable),
    stringsAsFactors = FALSE,
    row.names = NULL
  )

  if (isTRUE(boot)) {
    if (!is.null(seed)) {
      if (exists(".Random.seed", envir = .GlobalEnv)) {
        .old_seed <- get(".Random.seed", envir = .GlobalEnv)
        on.exit(assign(".Random.seed", .old_seed, envir = .GlobalEnv), add = TRUE)
      } else {
        on.exit(if (exists(".Random.seed", envir = .GlobalEnv)) rm(list = ".Random.seed", envir = .GlobalEnv), add = TRUE)
      }
      set.seed(seed)
    }
    n <- nrow(ratings)
    boot_a <- numeric(B)
    for (i in seq_len(B)) {
      idx <- sample(seq_len(n), n, replace = TRUE)
      bA  <- tryCatch(one_alpha(ratings[idx, , drop = FALSE])$alpha,
                      error = function(e) NA_real_)
      boot_a[i] <- bA
    }
    boot_a <- boot_a[is.finite(boot_a)]
    alpha_lo <- as.numeric(stats::quantile(boot_a, (1 - conf_level) / 2,
                                           na.rm = TRUE))
    alpha_hi <- as.numeric(stats::quantile(boot_a, 1 - (1 - conf_level) / 2,
                                           na.rm = TRUE))
    out <- rbind(out, data.frame(
      term  = c("lower_limit", "upper_limit", "B_used"),
      value = c(alpha_lo, alpha_hi, length(boot_a)),
      stringsAsFactors = FALSE, row.names = NULL
    ))
  }

  .as_dmar_tbl(out, conf_level = conf_level)
}

# ---- core ----

.kripp_distance <- function(a, b, level, all_vals) {
  switch(level,
    nominal  = as.integer(a != b),
    ordinal  = .kripp_ord_distance(a, b, all_vals),
    interval = (a - b)^2,
    # Identical ratio-scale values have zero distance. Guarding a == b before
    # the division avoids 0 / 0 = NaN for an identical zero pair (a == b == 0),
    # which otherwise propagates NaN into D_o / D_e and the alpha comparison.
    ratio    = if (a == b) 0 else ((a - b) / (a + b))^2
  )
}

.kripp_ord_distance <- function(a, b, all_vals) {
  # Krippendorff (2011) ordinal distance: based on rank frequencies.
  ord <- sort(unique(all_vals))
  freq <- table(factor(all_vals, levels = ord))
  if (a == b) return(0)
  lo <- min(a, b); hi <- max(a, b)
  idx_lo <- which(ord == lo); idx_hi <- which(ord == hi)
  if (length(idx_lo) == 0 || length(idx_hi) == 0)
    stop("Ordinal level requires values to appear in 'all_vals'.")
  inner <- if (idx_hi - idx_lo >= 2L)
    sum(freq[(idx_lo + 1L):(idx_hi - 1L)]) else 0
  (freq[idx_lo] / 2 + inner + freq[idx_hi] / 2)^2
}

.kripp_alpha_core <- function(R, level) {
  units <- nrow(R)

  # A unit rated only once forms no within-unit pair, so it contributes nothing
  # to Krippendorff's coincidence matrix and is excluded from the marginals and
  # from the normalizing count n. Only pairable units (>= 2 ratings) enter.
  unit_vals <- lapply(seq_len(units), function(u) {
    v <- R[u, ]
    v[!is.na(v)]
  })
  m_u <- vapply(unit_vals, length, integer(1))
  pairable <- m_u >= 2L
  if (!any(pairable))
    stop("Need at least one unit with 2 or more ratings.")

  # Values entering the coincidence matrix. Its marginals (and, for the ordinal
  # metric, the rank frequencies) come from pairable units only.
  all_vals <- unlist(unit_vals[pairable], use.names = FALSE)
  N <- length(all_vals)                       # total pairable values (Krippendorff's n)
  if (N < 2L)
    stop("Need at least 2 pairable values total.")

  # Observed disagreement (Krippendorff, 2011): each unit's pairwise distance
  # sum is weighted by 1 / (m_u - 1) -- the coincidence-matrix normalization --
  # and the total is divided by n. Pooling the pair-sums and the (m_u - 1) terms
  # into a single ratio, as before, drops the per-unit weighting whenever the
  # unit sizes m_u differ, which is the general case with missing ratings.
  num_obs <- 0
  for (u in which(pairable)) {
    vals <- unit_vals[[u]]
    mu   <- length(vals)
    pairsum <- 0
    for (j in seq_len(mu - 1L)) for (k in (j + 1L):mu)
      pairsum <- pairsum + .kripp_distance(vals[j], vals[k], level, all_vals)
    num_obs <- num_obs + 2 * pairsum / (mu - 1L)
  }
  D_o <- num_obs / N                          # within-unit average

  # Expected disagreement: the average distance over all pairs of pairable
  # values, treated as if drawn from the pooled marginals (Krippendorff, 2011).
  num_exp <- 0
  for (j in seq_len(N - 1L)) for (k in (j + 1L):N) {
    num_exp <- num_exp +
      .kripp_distance(all_vals[j], all_vals[k], level, all_vals)
  }
  D_e <- 2 * num_exp / (N * (N - 1L))         # between-unit average

  alpha <- if (D_e == 0) NA_real_ else 1 - D_o / D_e
  list(alpha = alpha, D_o = D_o, D_e = D_e, n_pairable = N)
}
