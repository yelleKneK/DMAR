# Tidy paired pairwise comparisons (within-subjects).
#' Paired Pairwise Comparisons With Multiple-Comparison Adjustment
#'
#' Computes all pairwise paired-\emph{t} comparisons among the levels of
#' a within-subjects factor and returns the mean difference, paired
#' SD, paired \emph{t}-statistic, degrees of freedom, raw and
#' adjusted \emph{p}-values, and a confidence interval on the mean
#' difference, one row per pair. \emph{p}-values are adjusted
#' across comparisons by the user-specified method (Bonferroni,
#' Holm, Hochberg, Hommel, BH, BY, or none).
#'
#' @param data Either an \eqn{n \times k} wide numeric matrix /
#'   data.frame (one row per subject, one column per condition), or a
#'   long-format \code{data.frame} together with \code{subject},
#'   \code{condition}, and \code{outcome} column names.
#' @param subject Long-format only: character name of the subject-id
#'   column.
#' @param condition Long-format only: character name of the within-
#'   subjects factor column.
#' @param outcome Long-format only: character name of the response
#'   column.
#' @param adjust Multiple-comparison adjustment method. One of
#'   \code{"holm"} (default), \code{"bonferroni"}, \code{"hochberg"},
#'   \code{"hommel"}, \code{"BH"}, \code{"BY"}, or \code{"none"}.
#'   Passed to \code{stats::p.adjust()}.
#' @param conf_level Family-wise confidence level for the per-pair
#'   CIs. Default \code{0.95}. CIs are computed at the per-pair
#'   nominal level (\eqn{1 - \alpha/m} under Bonferroni;
#'   \code{conf_level} otherwise).
#' @param bonferroni_ci Logical. If \code{TRUE}, the CIs use a
#'   per-pair confidence level of \eqn{1 - (1 - \mathrm{conf\_level}) / m}
#'   to give simultaneous \code{conf_level} coverage across the
#'   \eqn{m} comparisons (Bonferroni-corrected CIs). Default
#'   \code{FALSE}.
#'
#' @return A \code{data.frame} (class \code{dmar_tbl}) with one row per
#'   pair. Columns:
#'   \code{contrast} (the labeled difference, e.g.
#'   \code{"B - A"}), \code{mean_difference}, \code{sd_difference},
#'   \code{t_value}, \code{df}, \code{p_value} (raw),
#'   \code{p_adjusted}, \code{ci_lower}, \code{ci_upper},
#'   \code{n_pairs}. The table records \code{conf_level} in its
#'   \code{conf_level} attribute, which the console prints beneath the
#'   table and \code{\link{results_sentence}} reports. With
#'   \code{bonferroni_ci = TRUE} that is the simultaneous level of the
#'   \eqn{m} intervals, the convention of \code{\link{contrast_test}}
#'   with \code{adjust = "bonferroni"} and of the post hoc intervals
#'   (\code{\link{ci_tukey_kramer}} and the rest), and the line beneath
#'   the table names the Bonferroni adjustment and the level of each
#'   interval, \eqn{1 - (1 - \mathrm{conf\_level}) / m}.
#'
#' @details
#' \strong{Why a paired pairwise.} \code{stats::pairwise.t.test()}
#' returns a square matrix of \emph{p}-values, which doesn't compose
#' with the rest of the \pkg{DMAR} pipeline. This function returns
#' one row per comparison, matching the \code{data.frame(term, value)}
#' style used elsewhere.
#'
#' \strong{CI scale.} CIs are on the mean-difference scale (unstandardized).
#' When \code{bonferroni_ci = TRUE}, the per-pair confidence level is
#' \eqn{1 - (1 - \mathrm{conf\_level}) / m}, giving Bonferroni-style
#' simultaneous coverage. The CI is built from the paired-\emph{t}
#' distribution with \eqn{n - 1} degrees of freedom.
#'
#' @references
#' Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027).
#'   \emph{Designing experiments and analyzing data: A model comparison
#'   perspective} (4th ed.). Routledge. (See Chapter 11.)
#'
#' Holm, S. (1979). A simple sequentially rejective multiple test
#'   procedure. \emph{Scandinavian Journal of Statistics, 6}(2), 65--70.
#'
#' @seealso \code{\link[stats]{pairwise.t.test}},
#'   \code{\link{anova_within}}
#'
#' @examples
#' # 1. Wide-format input: 4 timepoints x 10 subjects.
#' set.seed(113)
#' n <- 10; k <- 4
#' Y <- matrix(rnorm(n * k, 0, 1), n, k) +
#'      matrix(rep(seq(0, 0.9, length.out = k), n), n, k, byrow = TRUE) +
#'      rnorm(n, 0, 1.5)
#' colnames(Y) <- paste0("T", 1:k)
#' pairwise_within(Y)
#'
#' # 2. Long-format input:
#' long <- data.frame(
#'   subject = factor(rep(1:n, times = k)),
#'   time    = factor(rep(paste0("T", 1:k), each = n)),
#'   y       = as.vector(Y)
#' )
#' pairwise_within(long, subject = "subject", condition = "time", outcome = "y")
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @keywords htest design
#'
#' @family within-subjects analysis
#' @family hypothesis tests
#'
#' @export

pairwise_within <- function(data,
                            subject = NULL, condition = NULL, outcome = NULL,
                            adjust = c("holm", "bonferroni", "hochberg",
                                       "hommel", "BH", "BY", "none"),
                            conf_level = 0.95,
                            bonferroni_ci = FALSE) {
  .dmar_check_conf_level(conf_level)
  adjust <- match.arg(adjust)
  if (conf_level <= 0 || conf_level >= 1)
    .dmar_refuse_argument("'conf_level' must be in (0, 1).", arg = "conf_level")

  if ((is.matrix(data) || is.data.frame(data)) &&
      is.null(subject) && is.null(condition) && is.null(outcome)) {
    Y <- as.matrix(data)
    if (!is.numeric(Y))
      .dmar_refuse_argument("Wide-format 'data' must be numeric.", arg = "data")
  } else {
    if (is.null(subject) || is.null(condition) || is.null(outcome))
      .dmar_refuse_argument(
        "For long-format input, supply 'subject', 'condition', and 'outcome'.",
        arg = c("subject", "condition", "outcome"))
    for (nm in c(subject, condition, outcome))
      if (!nm %in% names(data))
        .dmar_refuse_argument(
          sprintf("Column '%s' not in 'data'.", nm), arg = "data")
    Y <- stats::reshape(
      as.data.frame(data[, c(subject, condition, outcome)]),
      timevar   = condition,
      idvar     = subject,
      direction = "wide"
    )
    Y <- as.matrix(Y[, -1, drop = FALSE])
    storage.mode(Y) <- "double"
    colnames(Y) <- sub(paste0("^", outcome, "\\."), "", colnames(Y))
  }
  if (anyNA(Y))
    .dmar_refuse_argument(
      "Missing values are not supported in pairwise_within().", arg = NULL)
  if (ncol(Y) < 2L)
    .dmar_refuse_argument("Need at least 2 conditions.", arg = NULL)
  n <- nrow(Y); k <- ncol(Y); m <- k * (k - 1L) / 2L

  ci_conf <- if (bonferroni_ci) 1 - (1 - conf_level) / m else conf_level
  alpha   <- 1 - ci_conf
  t_crit  <- stats::qt(1 - alpha / 2, df = n - 1L)

  rows <- vector("list", m)
  idx <- 0L
  for (i in seq_len(k - 1L)) {
    for (j in (i + 1L):k) {
      idx <- idx + 1L
      d <- Y[, j] - Y[, i]
      m_d <- mean(d); s_d <- stats::sd(d)
      t_v <- if (s_d == 0) NA_real_ else m_d * sqrt(n) / s_d
      p_v <- if (is.na(t_v)) NA_real_
             else 2 * stats::pt(-abs(t_v), df = n - 1L)
      lo  <- m_d - t_crit * s_d / sqrt(n)
      hi  <- m_d + t_crit * s_d / sqrt(n)
      rows[[idx]] <- data.frame(
        contrast        = paste(colnames(Y)[j], "-", colnames(Y)[i]),
        mean_difference = m_d, sd_difference = s_d,
        t_value     = t_v, df = n - 1L, p_value = p_v,
        ci_lower     = lo, ci_upper = hi,
        n_pairs         = n,
        stringsAsFactors = FALSE
      )
    }
  }
  out <- do.call(rbind, rows)
  out$p_adjusted <- stats::p.adjust(out$p_value, method = adjust)
  out <- out[, c("contrast", "mean_difference", "sd_difference",
                 "t_value", "df", "p_value", "p_adjusted",
                 "ci_lower", "ci_upper", "n_pairs")]
  rownames(out) <- NULL
  # The table records the family level, as contrast_test() and the post
  # hoc intervals do for their simultaneous intervals, so that the console
  # prints it and results_sentence() on a row reports it; under Bonferroni
  # the line beneath the table also names the level of each interval.
  out <- .as_dmar_tbl(out, conf_level = conf_level)
  attr(out, "p_adjustment") <- adjust
  attr(out, "ci_adjustment") <- if (bonferroni_ci && m > 1L) "bonferroni" else "none"
  if (bonferroni_ci && m > 1L) {
    attr(out, "simultaneous") <- list(method = "Bonferroni",
                                      n_comparisons = m, each = ci_conf)
  }
  out
}
