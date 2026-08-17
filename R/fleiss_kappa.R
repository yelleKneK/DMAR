#' Fleiss's Kappa for Inter-Rater Agreement Among Multiple Raters
#'
#' Computes Fleiss's (1971) kappa coefficient of agreement among
#' \eqn{m \ge 2} raters who classify each of \eqn{N} subjects into one of
#' \eqn{k} nominal categories. Returns the point estimate together with the
#' asymptotic standard error and the Wald confidence interval, plus a test
#' of \eqn{H_0\!: \kappa = 0}.
#'
#' @param ratings A numeric \eqn{N \times k} matrix or \code{data.frame}:
#'   row \eqn{i} gives the counts of raters who assigned each of the
#'   \eqn{k} categories to subject \eqn{i}. Each row must sum to the same
#'   value \eqn{m} (the common number of raters per subject).
#' @param conf_level Confidence level for the interval (default
#'   \code{0.95}).
#' @param ci_method Interval method: \code{"wald"} (the default, the
#'   asymptotic interval from the Gwet (2008) linearization variance),
#'   \code{"percentile"} (bootstrap percentile), or \code{"bca"}
#'   (bootstrap bias-corrected and accelerated). The multirater kappa
#'   variance literature is unsettled, and Zapf, Castell, Morawietz,
#'   and Karch (2016) recommend bootstrap intervals in this setting:
#'   the subjects (rows) are resampled with replacement \code{B}
#'   times, kappa is recomputed on each resample, and the interval is
#'   read off the bootstrap distribution; the BCa variant additionally
#'   adjusts the quantile positions for median bias and acceleration
#'   (Efron & Tibshirani, 1993). The \code{se}, \code{z_value}, and
#'   \code{p_value} columns keep their asymptotic definitions under
#'   every \code{ci_method}; only the interval changes.
#' @param B Number of bootstrap replications when \code{ci_method} is
#'   \code{"percentile"} or \code{"bca"} (default \code{10000};
#'   ignored for \code{"wald"}).
#' @param seed Optional integer seed for the bootstrap. The default
#'   \code{NULL} uses the current state of the random number generator;
#'   a supplied seed is set internally and the prior state restored on
#'   exit.
#'
#' @return A one-row \code{data.frame} (class \code{dmar_tbl}) with columns
#'   \code{kappa}, \code{se} (asymptotic standard error of \eqn{\hat\kappa_F}
#'   used for the interval), \code{lower_limit}, \code{upper_limit},
#'   \code{z_value}, \code{p_value} (Wald test of \eqn{H_0\!: \kappa = 0}),
#'   \code{n_subjects}, \code{n_raters} (\eqn{m}), and
#'   \code{n_categories} (\eqn{k}).
#'
#' @details
#' For \eqn{n_{ij}} = the number of raters who assigned subject \eqn{i} to
#' category \eqn{j}, with \eqn{\sum_j n_{ij} = m} for every \eqn{i}, define
#' the marginal proportion of category \eqn{j} as
#' \eqn{p_j = \sum_i n_{ij} / (Nm)}, and the per-subject agreement
#' \deqn{P_i = \frac{1}{m(m-1)}\Bigl(\sum_j n_{ij}^2 - m\Bigr).}
#' Then Fleiss's kappa is
#' \deqn{\hat\kappa_F = \frac{\bar P - P_e}{1 - P_e}, \qquad
#'                     \bar P = \frac{1}{N}\sum_i P_i, \qquad
#'                     P_e = \sum_j p_j^2.}
#'
#' \strong{Standard error.} Two variances are involved, because the variance
#' of \eqn{\hat\kappa_F} under \eqn{H_0\!: \kappa = 0} is not its variance at
#' a nonzero value. The test of no agreement uses the null variance of
#' Fleiss, Nee, and Landis (1979, Equation 12), who corrected the standard
#' errors given in Fleiss (1971),
#' \deqn{\mathrm{Var}_0(\hat\kappa_F) = \frac{2\,\bigl(P_e + P_e^2 - 2\sum_j p_j^3\bigr)}{N\,m\,(m-1)\,(1-P_e)^2},}
#' and the reported \eqn{z} statistic and \emph{p}-value come from it. On
#' the Fleiss (1971) Table 1 example below this gives \eqn{z = 17.65},
#' matching \code{irr::kappam.fleiss}. The
#' confidence interval instead uses the linearization variance of Gwet
#' (2008, Section 6), which is consistent at the estimated
#' \eqn{\hat\kappa_F}: each subject \eqn{i} contributes an influence value
#' \eqn{\kappa_i^\ast} (Gwet's Equations 34 and 35), and
#' \eqn{\mathrm{Var}(\hat\kappa_F) = \sum_i (\kappa_i^\ast - \hat\kappa_F)^2 / \{N(N-1)\}},
#' Gwet's Equation 33 with the sampling fraction set to zero. (Gwet derives
#' the variance for the multiple-rater pi statistic, which is the same
#' estimator as Fleiss's kappa.)
#' The Wald confidence interval is
#' \eqn{\hat\kappa_F \pm z_{1-\alpha/2}\,\widehat{\mathrm{SE}}}, with the upper
#' limit truncated at 1. Using the null variance for the interval would
#' understate the standard error and give a spuriously narrow interval.
#'
#' \strong{Bootstrap interval.} The variance of multirater kappa is
#' unsettled in the literature, and Zapf, Castell, Morawietz, and Karch
#' (2016) recommend a bootstrap interval in this setting. With
#' \code{ci_method = "percentile"} or \code{"bca"} the subjects (the rows
#' of \code{ratings}) are resampled with replacement \code{B} times,
#' kappa is recomputed on each resample, and the interval is read off the
#' bootstrap distribution: the percentile interval takes the empirical
#' quantiles, and the BCa interval adjusts the quantile positions for
#' median bias and for acceleration (Efron & Tibshirani, 1993). Ask for
#' it when \eqn{N} is small or \eqn{\hat\kappa_F} is near a boundary,
#' where the Wald interval's coverage is least dependable. The \code{se},
#' \code{z_value}, and \code{p_value} columns keep their asymptotic
#' definitions under every \code{ci_method}; only the interval changes.
#' Bootstrap results vary from run to run; supply \code{seed} for
#' reproducibility.
#'
#' Fleiss's kappa is purely nominal (no weighting). For ordinal categories
#' with two raters, use \code{\link{cohen_kappa}} with quadratic weights;
#' for ordinal categories with three or more raters, an extension based on
#' the intraclass correlation (\code{\link{icc}}) is more appropriate.
#'
#' @references
#' Fleiss, J. L. (1971). Measuring nominal scale agreement among many
#'   raters. \emph{Psychological Bulletin, 76}(5), 378--382.
#'
#' Fleiss, J. L., Nee, J. C. M., & Landis, J. R. (1979). Large sample
#'   variance of kappa in the case of different sets of raters.
#'   \emph{Psychological Bulletin, 86}(5), 974--977.
#'
#' Efron, B., & Tibshirani, R. J. (1993). \emph{An introduction to the
#'   bootstrap}. New York, NY: Chapman & Hall/CRC.
#'
#' Gwet, K. L. (2008). Computing inter-rater reliability and its
#'   variance in the presence of high agreement. \emph{British Journal
#'   of Mathematical and Statistical Psychology, 61}(1), 29--48.
#'   \doi{10.1348/000711006X126600}
#'
#' Zapf, A., Castell, S., Morawietz, L., & Karch, A. (2016). Measuring
#'   inter-rater reliability for nominal data: Which coefficients and
#'   confidence intervals are appropriate? \emph{BMC Medical Research
#'   Methodology, 16}, 93. \doi{10.1186/s12874-016-0200-9}
#'
#' @examples
#' # Fleiss (1971) Table 1 example: 30 subjects rated by 6 raters into
#' # 5 diagnostic categories (Depression, Personality Disorder,
#' # Schizophrenia, Neurosis, Other). Each row of `ratings` gives, for one
#' # subject, the count of raters who chose each category (rows sum to 6).
#' # kappa = 0.430, matching Fleiss (1971).
#' fleiss_1971 <- matrix(c(
#'   0, 0, 0, 6, 0,
#'   0, 3, 0, 0, 3,
#'   0, 1, 4, 0, 1,
#'   0, 0, 0, 0, 6,
#'   0, 3, 0, 3, 0,
#'   2, 0, 4, 0, 0,
#'   0, 0, 4, 0, 2,
#'   2, 0, 3, 1, 0,
#'   2, 0, 0, 4, 0,
#'   0, 0, 0, 0, 6,
#'   1, 0, 0, 5, 0,
#'   1, 1, 0, 4, 0,
#'   0, 3, 3, 0, 0,
#'   1, 0, 0, 5, 0,
#'   0, 2, 0, 3, 1,
#'   0, 0, 5, 0, 1,
#'   3, 0, 0, 1, 2,
#'   5, 1, 0, 0, 0,
#'   0, 2, 0, 4, 0,
#'   1, 0, 2, 0, 3,
#'   0, 0, 0, 0, 6,
#'   0, 1, 0, 5, 0,
#'   0, 2, 0, 1, 3,
#'   2, 0, 0, 4, 0,
#'   1, 0, 0, 4, 1,
#'   0, 5, 0, 1, 0,
#'   4, 0, 0, 0, 2,
#'   0, 2, 0, 4, 0,
#'   1, 0, 5, 0, 0,
#'   0, 0, 0, 0, 6
#' ), nrow = 30, byrow = TRUE)
#' fleiss_kappa(fleiss_1971)
#'
#' # A bootstrap interval, which resamples the subjects (rows) with
#' # replacement and recomputes kappa on each resample. Not run here,
#' # because 2000 refits of kappa is more than a help page should do;
#' # the call is:
#' # fleiss_kappa(fleiss_1971, ci_method = "percentile", B = 2000,
#' #              seed = 113)
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso \code{\link{cohen_kappa}}, \code{\link{icc}}
#'
#' @family reliability
#'
#' @keywords htest
#'
#' @export
#' @import stats
fleiss_kappa <- function(ratings, conf_level = 0.95,
                         ci_method = c("wald", "percentile", "bca"),
                         B = 10000L, seed = NULL) {
  ci_method <- match.arg(ci_method)
  if (ci_method != "wald" &&
      (!is.numeric(B) || length(B) != 1L || B < 100)) {
    stop("'B' must be a single integer of at least 100 when ",
         "bootstrapping.", call. = FALSE)
  }
  if (is.data.frame(ratings)) ratings <- as.matrix(ratings)
  if (!is.matrix(ratings) || !is.numeric(ratings)) {
    stop("'ratings' must be a numeric N x k matrix.", call. = FALSE)
  }
  if (anyNA(ratings)) {
    stop("Missing values are not supported by fleiss_kappa(); ",
         "remove or impute first.", call. = FALSE)
  }
  if (!is.numeric(conf_level) || length(conf_level) != 1L ||
      conf_level <= 0 || conf_level >= 1) {
    stop("'conf_level' must be a single number in (0, 1).", call. = FALSE)
  }

  N <- nrow(ratings); k <- ncol(ratings)
  if (N < 2L) stop("At least 2 subjects are required.", call. = FALSE)
  if (k < 2L) stop("At least 2 categories are required.", call. = FALSE)

  row_sums <- rowSums(ratings)
  m <- row_sums[1L]
  if (any(row_sums != m)) {
    stop("Every row of 'ratings' must sum to the same value (the common ",
         "number of raters per subject).", call. = FALSE)
  }
  if (m < 2L) stop("Need at least 2 raters per subject.", call. = FALSE)

  p_j   <- colSums(ratings) / (N * m)
  P_i   <- (rowSums(ratings^2) - m) / (m * (m - 1))
  P_bar <- mean(P_i)
  P_e   <- sum(p_j^2)

  if (1 - P_e <= .Machine$double.eps) {
    kappa <- NA_real_; se <- NA_real_; z <- NA_real_
    p    <- NA_real_; lo <- NA_real_; hi <- NA_real_
  } else {
    kappa <- (P_bar - P_e) / (1 - P_e)

    # Standard error for the confidence interval. The Fleiss, Nee, and Landis
    # (1979) closed-form variance below is derived under H_0: kappa = 0 and is
    # the correct null variance for the z-test of no agreement, but it is not
    # the variance of kappa-hat at a nonzero value and so must not be used to
    # build the interval. The interval instead uses the linearization variance
    # of Gwet (2008, Equations 33 to 35, with the sampling fraction zero),
    # which is consistent at the estimated kappa: each subject contributes an
    # influence value, and the variance is the sample variance of those
    # values. This reproduces the standard error and interval reported by
    # Gwet's irrCAC implementation.
    r_i     <- m
    sum_q   <- rowSums(ratings * (ratings - 1))
    pa_i    <- sum_q / (r_i * (r_i - 1))
    kappa_i <- (pa_i - P_e) / (1 - P_e)
    pe_i    <- as.vector(ratings %*% p_j) / r_i
    kappa_i_star <- kappa_i - 2 * (1 - kappa) * (pe_i - P_e) / (1 - P_e)
    var_ci  <- sum((kappa_i_star - kappa)^2) / (N * (N - 1))
    se <- sqrt(max(0, var_ci))
    z_crit <- stats::qnorm(1 - (1 - conf_level) / 2)
    lo <- kappa - z_crit * se
    hi <- min(1, kappa + z_crit * se)

    # Test of H_0: kappa = 0 uses the Fleiss, Nee, and Landis (1979) null
    # variance, which corrected the Fleiss (1971) standard errors and is
    # valid only under that hypothesis. With q_j = 1 - p_j their printed
    # form is 2 ((sum p_j q_j)^2 - sum p_j q_j (q_j - p_j)) /
    # (N m (m-1) (sum p_j q_j)^2); the expression below is its algebraic
    # simplification, and it reproduces the z of irr::kappam.fleiss.
    var_null <- (2 / (N * m * (m - 1) * (1 - P_e)^2)) *
                (P_e + P_e^2 - 2 * sum(p_j^3))
    se_null <- sqrt(max(0, var_null))
    z  <- if (se_null == 0) NA_real_ else kappa / se_null
    p  <- if (is.na(z)) NA_real_ else 2 * stats::pnorm(-abs(z))
  }

  if (ci_method != "wald" && is.finite(kappa)) {
    kappa_at <- function(idx) {
      rt <- ratings[idx, , drop = FALSE]
      pj <- colSums(rt) / (length(idx) * m)
      Pb <- mean((rowSums(rt^2) - m) / (m * (m - 1)))
      Pe <- sum(pj^2)
      if (1 - Pe <= .Machine$double.eps) return(NA_real_)
      (Pb - Pe) / (1 - Pe)
    }
    if (!is.null(seed)) {
      has_old <- exists(".Random.seed", envir = globalenv())
      old_seed <- if (has_old) get(".Random.seed", envir = globalenv())
      on.exit({
        if (has_old) assign(".Random.seed", old_seed, envir = globalenv())
        else if (exists(".Random.seed", envir = globalenv()))
          rm(".Random.seed", envir = globalenv())
      }, add = TRUE)
      set.seed(seed)
    }
    B <- as.integer(B)
    boots <- vapply(seq_len(B), function(b) {
      kappa_at(sample.int(N, N, replace = TRUE))
    }, numeric(1))
    n_bad <- sum(!is.finite(boots))
    if (n_bad > 0L) {
      boots <- boots[is.finite(boots)]
      if (length(boots) < 100L) {
        stop("Only ", length(boots), " of ", B, " bootstrap replications ",
             "returned a defined kappa; the interval would not be ",
             "trustworthy.", call. = FALSE)
      }
      warning(n_bad, " of ", B, " bootstrap replications left kappa ",
              "undefined (chance agreement 1) and were dropped; the ",
              "interval is computed from the ", length(boots),
              " that did.", call. = FALSE)
    }
    alpha_2 <- (1 - conf_level) / 2
    if (ci_method == "percentile") {
      lims <- stats::quantile(boots, c(alpha_2, 1 - alpha_2),
                              names = FALSE)
    } else {
      # BCa: median bias from the bootstrap distribution, acceleration
      # from the jackknife (Efron & Tibshirani, 1993).
      z0 <- stats::qnorm(mean(boots < kappa))
      jack <- vapply(seq_len(N), function(i) {
        kappa_at(seq_len(N)[-i])
      }, numeric(1))
      jack <- jack[is.finite(jack)]
      jm <- mean(jack)
      acc <- sum((jm - jack)^3) / (6 * (sum((jm - jack)^2))^1.5)
      zq <- stats::qnorm(c(alpha_2, 1 - alpha_2))
      adj <- stats::pnorm(z0 + (z0 + zq) / (1 - acc * (z0 + zq)))
      lims <- stats::quantile(boots, adj, names = FALSE)
    }
    lo <- lims[1L]
    hi <- min(1, lims[2L])
  }

  out <- data.frame(
    kappa        = kappa,
    se           = se,
    lower_limit  = lo,
    upper_limit  = hi,
    z_value      = z,
    p_value      = p,
    n_subjects   = N,
    n_raters     = m,
    n_categories = k,
    stringsAsFactors = FALSE,
    row.names    = NULL
  )
  out <- .as_dmar_tbl(out, conf_level = conf_level)
  attr(out, "ci_method") <- ci_method
  if (ci_method != "wald" && is.finite(kappa)) {
    attr(out, "B_used") <- length(boots)
  }
  out
}
