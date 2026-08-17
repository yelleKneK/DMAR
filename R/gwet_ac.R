# Gwet's AC1 / AC2 chance-corrected agreement coefficient.
#' Gwet's AC1 and AC2 Chance-Corrected Agreement Coefficients
#'
#' Computes Gwet's AC1 (nominal data; Gwet, 2008) and AC2 (ordinal
#' data with user-supplied weights; Gwet, 2014) chance-corrected
#' agreement coefficients for two or more raters. AC1/AC2 are more robust than
#' Cohen's \eqn{\kappa} to extreme marginal-prevalence imbalance and
#' the trait-distribution paradox.
#'
#' @param ratings A units \eqn{\times} raters matrix or
#'   \code{data.frame}. Rows = units; columns = raters. \code{NA}
#'   entries are allowed.
#' @param weights One of \code{"unweighted"} (default; computes AC1)
#'   or \code{"linear"} / \code{"quadratic"} (compute AC2 with
#'   weight matrices used by weighted-\eqn{\kappa} conventions).
#' @param conf_level Confidence level. Default \code{0.95}.
#'
#' @return A \code{data.frame} with rows for the point estimate
#'   \eqn{\widehat{\mathrm{AC}}}, the standard error, the CI lower
#'   and upper limits, the percent agreement \eqn{p_a}, and the
#'   chance-agreement term \eqn{p_e}.
#'
#' @details
#' \strong{Coefficient.}
#' \deqn{\widehat{\mathrm{AC}} \;=\; \frac{p_a - p_e}{1 - p_e},}
#' identical to Cohen's \eqn{\kappa} in structure but with a different
#' chance-correction \eqn{p_e}:
#' \deqn{p_e \;=\; \frac{T_w}{Q (Q - 1)} \sum_{k = 1}^{Q} \pi_k (1 - \pi_k),}
#' where \eqn{\pi_k} is the mean within-unit proportion of
#' category \eqn{k}, \eqn{Q} is the number of categories, and
#' \eqn{T_w} is the sum of all entries of the weight matrix. For
#' nominal data with unit weights (AC1), \eqn{T_w = Q} and \eqn{p_e}
#' reduces to \eqn{(1 / (Q - 1)) \sum_k \pi_k (1 - \pi_k)}.
#'
#' \strong{Why AC over \eqn{\kappa}.} \eqn{\kappa} can be near zero
#' even when raters agree on almost every unit if the trait is rare
#' or very common (the "kappa paradox"; Feinstein & Cicchetti, 1990).
#' Gwet's AC keeps the same chance-correction logic but uses a less
#' extreme reference distribution.
#'
#' \strong{Variance.} The SE is Gwet's (2008) linearization variance,
#' \deqn{\mathrm{Var}(\widehat{\mathrm{AC}}) \;=\;
#'   \frac{1 - f}{n (n - 1)} \sum_i (\widehat{\mathrm{AC}}_i^{*}
#'   - \widehat{\mathrm{AC}})^2,}
#' where \eqn{\widehat{\mathrm{AC}}_i^{*}} is the \eqn{i}th unit's
#' influence value, combining its agreement and chance-term
#' contributions, \eqn{n} is the number of units, and \eqn{f} is the
#' sampling fraction (\eqn{0} for an infinite target population). The
#' interval is \eqn{\widehat{\mathrm{AC}} \pm t_{1 - \alpha / 2,\, n - 1}
#' \mathit{SE}}, with the upper limit truncated at \eqn{1}. These
#' quantities match Gwet's (2014) reference software.
#'
#' @references
#' Feinstein, A. R., & Cicchetti, D. V. (1990). High agreement but
#'   low kappa: I. The problems of two paradoxes. \emph{Journal of
#'   Clinical Epidemiology, 43}(6), 543--549.
#'   \doi{10.1016/0895-4356(90)90158-L}
#'
#' Gwet, K. L. (2008). Computing inter-rater reliability and its
#'   variance in the presence of high agreement. \emph{British Journal
#'   of Mathematical and Statistical Psychology, 61}(1), 29--48.
#'   \doi{10.1348/000711006X126600}
#'
#' Gwet, K. L. (2014). \emph{Handbook of inter-rater reliability}
#'   (4th ed.). Advanced Analytics, LLC.
#'
#' @seealso \code{\link{cohen_kappa}}, \code{\link{fleiss_kappa}},
#'   \code{\link{krippendorff_alpha}}
#'
#' @examples
#' # 1. Unweighted AC1, two raters, nominal:
#' set.seed(113)
#' r1 <- sample(c("A", "B", "C"), 50, replace = TRUE)
#' r2 <- ifelse(runif(50) < 0.8, r1, sample(c("A", "B", "C"), 50, TRUE))
#' gwet_ac(cbind(r1, r2))
#'
#' # 2. AC2 with linear weights, ordinal scale 1-5:
#' set.seed(113)
#' r1 <- sample(1:5, 60, replace = TRUE)
#' r2 <- pmin(5, pmax(1, r1 + sample(-1:1, 60, replace = TRUE)))
#' gwet_ac(cbind(r1, r2), weights = "linear")
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @keywords htest
#'
#' @family agreement and measurement
#'
#' @export

gwet_ac <- function(ratings,
                    weights = c("unweighted", "linear", "quadratic"),
                    conf_level = 0.95) {
  weights <- match.arg(weights)
  if (is.data.frame(ratings)) ratings <- as.matrix(ratings)
  if (!is.matrix(ratings))
    stop("'ratings' must be a matrix or data.frame.")
  if (ncol(ratings) < 2L)
    stop("Need at least 2 raters (columns).")
  if (conf_level <= 0 || conf_level >= 1)
    stop("'conf_level' must be in (0, 1).")

  # A unit with no ratings (an all-missing row) contributes nothing to any
  # agreement or chance term, but leaving it in inflates n and divides the
  # category proportions by zero raters, turning AC and its SE into NaN. Drop
  # such units up front so N, the marginals, the influence values, the SE, and
  # the degrees of freedom are all formed from units that were actually rated.
  ratings <- ratings[rowSums(!is.na(ratings)) > 0L, , drop = FALSE]
  if (nrow(ratings) < 1L)
    stop("Need at least one unit with a rating.", call. = FALSE)

  cats <- sort(unique(as.vector(ratings[!is.na(ratings)])))
  Q    <- length(cats)
  if (Q < 2L) stop("Need at least 2 categories.")

  W <- switch(weights,
    unweighted = diag(Q),
    linear     = .gwet_weight_matrix(cats, "linear"),
    quadratic  = .gwet_weight_matrix(cats, "quadratic")
  )
  # Gwet's T_w is the sum of every entry of the weight matrix (Gwet,
  # 2014). Using the trace divided by Q, which equals 1 for any matrix with a
  # unit diagonal, collapses AC2's chance term back to the unweighted case and
  # leaves the weights inert.
  T_w <- sum(W)

  # Matrix formulation following Gwet (2008, 2014), valid for any number of
  # raters and for partially missing rows. agree_mat[i, k] is the number of
  # raters who placed unit i in category k.
  n <- nrow(ratings)
  N_pop <- Inf                                   # infinite target population; no fpc
  f     <- n / N_pop
  agree_mat <- matrix(0, nrow = n, ncol = Q)
  for (k in seq_len(Q)) {
    is_k <- ratings == cats[k]
    is_k[is.na(is_k)] <- FALSE
    agree_mat[, k] <- rowSums(is_k)
  }
  ri_vec  <- rowSums(agree_mat)                   # raters per unit
  agree_w <- agree_mat %*% W                      # weighted counts (W symmetric)
  sum_q   <- rowSums(agree_mat * (agree_w - 1))   # weighted concordant pairs per unit
  n2more  <- sum(ri_vec >= 2L)
  if (n2more < 1L)
    stop("Need at least one unit rated by 2 or more raters.", call. = FALSE)

  # Percent agreement, averaged over units rated at least twice.
  den_pa <- ri_vec * (ri_vec - 1)
  pa <- sum((sum_q / den_pa)[ri_vec >= 2L]) / n2more

  # Chance agreement from the mean within-unit category proportions.
  pi_vec <- colSums(agree_mat / ri_vec) / n
  pe <- T_w * sum(pi_vec * (1 - pi_vec)) / (Q * (Q - 1))
  AC <- (pa - pe) / (1 - pe)

  # Gwet's (2008) linearization variance: each unit's influence combines its
  # agreement contribution with its marginal-probability (chance) contribution.
  den_iv  <- den_pa - (den_pa == 0)               # guard 0 / 0 for one-rating units
  pa_iv   <- sum_q / den_iv
  pe_r2   <- pe * (ri_vec >= 2L)
  ac_iv   <- (n / n2more) * (pa_iv - pe_r2) / (1 - pe)
  pe_iv   <- (T_w / (Q * (Q - 1))) * drop(agree_mat %*% (1 - pi_vec)) / ri_vec
  ac_iv_x <- ac_iv - 2 * (1 - AC) * (pe_iv - pe) / (1 - pe)

  se_AC <- NA_real_; lo <- NA_real_; hi <- NA_real_
  if (n >= 2L) {
    var_AC <- ((1 - f) / (n * (n - 1))) * sum((ac_iv_x - AC)^2)
    se_AC  <- sqrt(var_AC)
    tcrit  <- stats::qt(1 - (1 - conf_level) / 2, n - 1)
    lo     <- AC - tcrit * se_AC
    hi     <- min(1, AC + tcrit * se_AC)
  }

  out <- data.frame(
    term  = c("gwet_ac", "se", "lower_limit", "upper_limit",
              "percent_agreement", "chance_agreement", "n_units"),
    value = c(AC, se_AC, lo, hi, pa, pe, n2more),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  .as_dmar_tbl(out, conf_level = conf_level)
}

.gwet_weight_matrix <- function(cats, kind) {
  Q <- length(cats)
  M <- max(cats) - min(cats)
  W <- matrix(0, Q, Q)
  for (i in seq_len(Q)) for (j in seq_len(Q)) {
    diff_ij <- abs(cats[i] - cats[j])
    W[i, j] <- switch(kind,
      linear    = 1 - diff_ij / M,
      quadratic = 1 - (diff_ij / M)^2
    )
  }
  W
}
