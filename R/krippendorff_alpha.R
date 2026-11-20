# Krippendorff's alpha inter-rater agreement coefficient.
#' Krippendorff's \eqn{\alpha} Inter-Rater Agreement
#'
#' Computes Krippendorff's (1980, 2004, 2011) \eqn{\alpha}, a
#' chance-corrected inter-rater agreement coefficient. Unlike
#' \code{\link{cohen_kappa}} (two raters; nominal categories, or ordinal
#' ones through weights) or \code{\link{fleiss_kappa}} (multiple raters,
#' nominal categories, no missing ratings),
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
#' @param boot Logical. If \code{TRUE}, fills the interval rows with a
#'   bootstrap percentile CI. The default \code{FALSE} runs no bootstrap
#'   and leaves those rows \code{NA} (much faster).
#' @param B Number of bootstrap resamples when \code{boot = TRUE}, a
#'   single whole number. Default \code{1000L}, because every resample
#'   recomputes the coincidences over all pairable values. The fewest
#'   accepted is 400 at the default \code{conf_level} (200 at 0.90, and
#'   2000 at 0.99, which the default does not meet); a smaller \code{B} is
#'   refused before anything is computed (class \code{dmar_count_invalid}).
#'   See \emph{Bootstrap replications} under \emph{Details} and
#'   \code{\link{dmar_bootstrap}}.
#' @param seed Optional integer seed for reproducibility of the
#'   bootstrap. Default \code{NULL}, under which the resamples are drawn
#'   from, and advance, the current random number stream; supply an
#'   integer for reproducibility, which is set for the call only, the
#'   caller's random number state being restored when it returns.
#'
#' @return A \code{data.frame} (class \code{dmar_tbl}) with the same
#'   seven rows whatever \code{boot}: the point estimate
#'   \eqn{\hat\alpha} (\code{krippendorff_alpha}), the observed
#'   disagreement \eqn{D_o} (\code{D_observed}), the expected
#'   disagreement \eqn{D_e} (\code{D_expected}), the number of pairable
#'   values (\code{n_pairable}), the lower and upper bootstrap CI limits
#'   (\code{ci_lower}, \code{ci_upper}), and \code{B_used}, the
#'   number of resamples that returned a finite value and so entered the
#'   interval. The last three are \code{NA} when \code{boot = FALSE}, since
#'   no bootstrap was run. With \code{boot = TRUE} the table also carries
#'   the \code{conf_level} attribute (when an interval was formed) and the
#'   attributes \code{B}, \code{B_used}, \code{n_dropped},
#'   \code{B_required}, and \code{bootstrap_status} (see \emph{Bootstrap
#'   replications} under Details).
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
#' When the coefficient is undefined for the
#' data themselves (every value alike), no bootstrap is run and the
#' interval rows stay \code{NA}. No
#' closed-form sampling variance is in general use for Krippendorff's
#' alpha across its measurement levels and missing data patterns, so the
#' bootstrap is the interval Krippendorff recommends (Krippendorff,
#' 2011; Hayes & Krippendorff, 2007). \code{B = 1000L} typically gives a
#' stable CI to two decimal places. The bootstrap is opt-in
#' (\code{boot = FALSE} by default, which returns the point estimate
#' alone and is much faster); ask for it whenever the coefficient is
#' being reported rather than explored, since a point estimate on its
#' own says nothing about how precisely \eqn{\alpha} is determined.
#'
#' \strong{Bootstrap replications.} The two limits are quantiles of the
#' \code{B} replicated values of \eqn{\alpha}, so the interval is run only
#' when \code{B} puts at least 10 resamples beyond each limit (400 at the
#' 95 percent level). A replication in which \eqn{\alpha} is undefined (a
#' resample without enough pairable values, or one whose values are all
#' alike, so that the expected disagreement is zero) is dropped and
#' counted, never replaced, and the interval is computed from the rest. The
#' result records \code{B} (requested), \code{B_used} (the number the
#' interval rests on), \code{n_dropped}, and \code{B_required} (the fewest
#' replications a request may ask for; that count is held to \code{B}, not
#' to \code{B_used}) as attributes, and
#' \code{B_used} also as a row, with \code{bootstrap_status}
#' (\code{"complete"}, \code{"dropped"}, \code{"exploratory"},
#' \code{"insufficient_count"}, or \code{"insufficient_share"}). One
#' warning of class \code{dmar_bootstrap_replications_dropped} states the
#' counts when any replication was dropped, with a caution when the dropped
#' replications are at least the 2.5 percent of \code{B} that lies beyond
#' one limit at the 95 percent level, because they could then decide by
#' themselves where that limit falls. When more than half were dropped, the
#' two limit rows are \code{NA},
#' the estimate is returned, and the warning has class
#' \code{dmar_bootstrap_insufficient}.
#' \code{options(dmar.bootstrap_exploratory = TRUE)} lets a smaller
#' \code{B} (2 or more) run for a first look or a demonstration; every such
#' call gives a warning of class \code{dmar_bootstrap_exploratory} and the
#' result is not one to report. See \code{\link{dmar_bootstrap}} for the
#' rule and its reasons.
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
#' # recomputes alpha on each of B resamples of the units and fills the
#' # ci_lower, ci_upper, and B_used rows. B = 400 is the fewest
#' # resamples a 95 percent interval is reported from and keeps the
#' # example quick; a reported interval deserves the default B = 1000.
#' krippendorff_alpha(cbind(r1, r2), level = "interval",
#'                    boot = TRUE, B = 400L, seed = 113)
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
  .dmar_check_conf_level(conf_level)
  level <- match.arg(level)
  if (is.data.frame(ratings)) ratings <- as.matrix(ratings)
  if (!is.matrix(ratings))
    .dmar_refuse_argument(
      "'ratings' must be a matrix or data.frame.", arg = "ratings")
  if (conf_level <= 0 || conf_level >= 1)
    .dmar_refuse_argument("'conf_level' must be in (0, 1).", arg = "conf_level")
  # B is checked only when it is used, by the rule of ?dmar_bootstrap: 10
  # replications beyond each limit of the percentile interval (400 at the
  # 95 percent level). Through 1.0.1 it was not checked at all, and B = 10
  # gave a percentile interval from ten resamples.
  if (isTRUE(boot)) {
    B <- .dmar_check_B(B, tail = (1 - conf_level) / 2,
                       conf_level = conf_level)
  }

  one_alpha <- function(R) {
    .kripp_alpha_core(R, level)
  }

  obs <- one_alpha(ratings)

  # One shape whatever boot: the interval rows and the count of resamples
  # behind them are NA unless the bootstrap runs.
  out <- data.frame(
    term  = c("krippendorff_alpha", "D_observed", "D_expected", "n_pairable",
              "ci_lower", "ci_upper", "B_used"),
    value = c(obs$alpha, obs$D_o, obs$D_e, obs$n_pairable,
              NA_real_, NA_real_, NA_real_),
    stringsAsFactors = FALSE,
    row.names = NULL
  )

  # No interval is formed for a coefficient the data leave undefined
  # (every value alike, so the expected disagreement is zero), as
  # cohen_kappa() forms none for an undefined kappa.
  formed <- isTRUE(boot) && is.finite(obs$alpha)
  if (formed) {
    .dmar_local_seed(seed)
    n <- nrow(ratings)
    boot_a <- numeric(B)
    for (i in seq_len(B)) {
      idx <- sample(seq_len(n), n, replace = TRUE)
      bA  <- tryCatch(one_alpha(ratings[idx, , drop = FALSE])$alpha,
                      error = function(e) NA_real_)
      boot_a[i] <- bA
    }
    # A resample on which alpha cannot be computed is dropped and counted,
    # with one warning, as the package's other bootstraps do; through
    # 1.0.1 such resamples were dropped without a word.
    # The verdict of the rule (R/dmar_bootstrap.R): more than half of
    # the resamples dropped leave the limits NA under a classed warning
    # (through 1.0.x the call stopped when fewer than 100 remained).
    boot_a <- boot_a[is.finite(boot_a)]
    verdict <- .dmar_boot_settle(
      B, length(boot_a), tail = (1 - conf_level) / 2,
      conf_level = conf_level, fun = "krippendorff_alpha",
      quantity = "Krippendorff's alpha",
      dropped_because = paste0("left alpha undefined (too few pairable ",
                               "values, or values all alike)"))
    alpha_lo <- alpha_hi <- NA_real_
    if (verdict$formed) {
      alpha_lo <- as.numeric(stats::quantile(boot_a, (1 - conf_level) / 2))
      alpha_hi <- as.numeric(stats::quantile(boot_a,
                                             1 - (1 - conf_level) / 2))
    }
    out$value[match(c("ci_lower", "ci_upper", "B_used"), out$term)] <-
      c(alpha_lo, alpha_hi, length(boot_a))
  }

  # The confidence level is recorded only beside an interval that was
  # formed, as in htmt() and average_variance_extracted(), so no
  # "Confidence level" line is printed under NA limits.
  out <- .as_dmar_tbl(out, conf_level = if (formed && verdict$formed)
    conf_level)
  if (formed) out <- .dmar_boot_attrs(out, verdict)
  out
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
    .dmar_refuse_argument(
      "Ordinal level requires values to appear in 'all_vals'.",
      arg = "all_vals")
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
    .dmar_refuse_argument(
      "Need at least one unit with 2 or more ratings.", arg = NULL)

  # Values entering the coincidence matrix. Its marginals (and, for the ordinal
  # metric, the rank frequencies) come from pairable units only.
  all_vals <- unlist(unit_vals[pairable], use.names = FALSE)
  N <- length(all_vals)                       # total pairable values (Krippendorff's n)
  if (N < 2L)
    .dmar_refuse_argument("Need at least 2 pairable values total.", arg = NULL)

  # Observed disagreement (Krippendorff, 2011): each unit's pairwise distance
  # sum is weighted by 1 / (m_u - 1), the coincidence matrix normalization,
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
  # N (N - 1) is formed as a double: N is a length, an integer, and the
  # integer product overflowed to NA from N = 46,342 (23,171 units rated
  # twice), which stopped the call with "missing value where TRUE/FALSE
  # needed" after the pairs had been summed.
  num_exp <- .kripp_pair_sum(all_vals, level)
  D_e <- 2 * num_exp / (N * (N - 1))          # between-unit average

  alpha <- if (D_e == 0) NA_real_ else 1 - D_o / D_e
  list(alpha = alpha, D_o = D_o, D_e = D_e, n_pairable = N)
}

# The sum of the distances over all N (N - 1) / 2 pairs of the pairable
# values, the numerator of the expected disagreement, in the order of the
# values.
.kripp_pair_sum <- function(all_vals, level) {
  N <- length(all_vals)
  num_exp <- 0
  for (j in seq_len(N - 1L)) for (k in (j + 1L):N) {
    num_exp <- num_exp +
      .kripp_distance(all_vals[j], all_vals[k], level, all_vals)
  }
  num_exp
}
