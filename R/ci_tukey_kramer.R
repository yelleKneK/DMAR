# Tukey-Kramer simultaneous CIs for all pairwise contrasts.
#' Tukey-Kramer Simultaneous Confidence Intervals for Pairwise Contrasts
#'
#' Computes the Tukey-Kramer simultaneous confidence intervals for all
#' \eqn{a (a - 1) / 2} pairwise contrasts among \eqn{a} group means in
#' a one-way design with possibly unequal group sample sizes (Tukey,
#' 1953; Kramer, 1956; Hayter, 1984), and returns the result in tidy
#' long form. Each interval has individual coverage at least the
#' specified \code{conf_level} and the family-wise coverage is at
#' least \code{conf_level}.
#'
#' @param x Either (a) a fitted \code{\link[stats]{lm}} or
#'   \code{\link[stats]{aov}} object with a single one-way factor
#'   predictor, or (b) a numeric vector of observations, in which
#'   case \code{group} must also be supplied.
#' @param group When \code{x} is a vector, a factor (or coercible to
#'   factor) of group labels, same length as \code{x}.
#' @param conf_level Family-wise confidence level. Default \code{0.95}.
#'
#' @return A \code{data.frame} with one row per pairwise
#'   contrast. Columns: \code{contrast} (e.g., \code{"B - A"}),
#'   \code{mean_difference}, \code{se}, \code{q_statistic}
#'   (the studentized-range \eqn{q}), \code{lower_limit},
#'   \code{upper_limit}, \code{p_adjusted}.
#'
#' @details
#' \strong{Formula.} For groups \eqn{i, j} with means \eqn{\bar y_i,
#' \bar y_j} and sample sizes \eqn{n_i, n_j}, the Tukey-Kramer
#' simultaneous CI is
#' \deqn{\bar y_i - \bar y_j \;\pm\; q_{\alpha, a, \nu}
#'   \sqrt{\frac{\mathit{MS}_E}{2} \left(\frac{1}{n_i} + \frac{1}{n_j}\right)},}
#' where \eqn{q_{\alpha, a, \nu}} is the upper \eqn{\alpha} quantile of
#' the studentized-range distribution with \eqn{a} groups and \eqn{\nu}
#' error degrees of freedom (\code{stats::qtukey()}).
#'
#' \strong{Why Tukey-Kramer.} Hayter (1984) proved that the
#' Tukey-Kramer procedure is conservative for unbalanced designs (the
#' coverage probability is at least \code{conf_level}). For balanced
#' designs it reduces to Tukey's HSD and the coverage is exactly
#' \code{conf_level}.
#'
#' \strong{Adjusted \emph{p}-values.} Each pairwise \eqn{p}-value is
#' computed from the studentized-range distribution:
#' \eqn{p = 1 - \mathrm{ptukey}(|q|, a, \nu)}.
#'
#' @references
#' Hayter, A. J. (1984). A proof of the conjecture that the
#'   Tukey-Kramer multiple comparisons procedure is conservative.
#'   \emph{Annals of Statistics, 12}(1), 61--75.
#'
#' Kramer, C. Y. (1956). Extension of multiple range tests to group
#'   means with unequal numbers of replications. \emph{Biometrics,
#'   12}(3), 307--310.
#'
#' Tukey, J. W. (1953). \emph{The problem of multiple comparisons}.
#'   Unpublished manuscript, Princeton University.
#'
#' @seealso \code{\link{cv_tukey_hsd}}, \code{\link{ci_dunnett}},
#'   \code{\link{ci_scheffe}}, \code{\link[stats]{TukeyHSD}}
#'
#' @examples
#' # 1. Balanced one-way: the six marketing panels of the test_market
#' #    data, four outlets per panel, so the procedure is exactly Tukey's
#' #    HSD. Panels 5 and 6 separate from panel 1.
#' fit <- lm(brand_movement ~ panel, data = test_market)
#' ci_tukey_kramer(fit)
#'
#' # 2. Same data via vector / group interface:
#' ci_tukey_kramer(test_market$brand_movement, group = test_market$panel)
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @keywords htest
#'
#' @family hypothesis tests
#'
#' @export

ci_tukey_kramer <- function(x, group = NULL, conf_level = 0.95) {
  if (conf_level <= 0 || conf_level >= 1)
    stop("'conf_level' must be in (0, 1).")

  if (inherits(x, c("lm", "aov"))) {
    mf <- stats::model.frame(x)
    if (ncol(mf) != 2L)
      stop("ci_tukey_kramer() requires a fitted lm/aov with one factor predictor.")
    y <- mf[[1L]]; g <- mf[[2L]]
    if (!is.factor(g)) g <- factor(g)
    tab <- stats::anova(x); rn <- trimws(rownames(tab))
    ms_e <- tab[rn == "Residuals", "Mean Sq"]
    df_e <- tab[rn == "Residuals", "Df"]
  } else if (is.numeric(x) && !is.null(group)) {
    if (length(x) != length(group))
      stop("'x' and 'group' must be the same length.")
    g <- factor(group)
    y <- x
    ok <- !is.na(y) & !is.na(g)
    y <- y[ok]; g <- g[ok]
    # Compute MS_E from a one-way ANOVA:
    fit <- stats::aov(y ~ g)
    tab <- stats::anova(fit); rn <- trimws(rownames(tab))
    ms_e <- tab[rn == "Residuals", "Mean Sq"]
    df_e <- tab[rn == "Residuals", "Df"]
  } else {
    stop("Supply either a fitted lm/aov, or (x, group) vectors.")
  }

  lev <- levels(g)
  a <- length(lev); m <- a * (a - 1L) / 2L
  if (a < 2L) stop("Need at least 2 groups.")
  means <- tapply(y, g, mean)
  ns    <- tapply(y, g, length)

  q_crit <- stats::qtukey(conf_level, a, df_e)

  rows <- vector("list", m); idx <- 0L
  for (i in seq_len(a - 1L)) {
    for (j in (i + 1L):a) {
      idx <- idx + 1L
      diff_ij <- as.numeric(means[j] - means[i])
      se      <- sqrt(ms_e / 2 * (1 / ns[i] + 1 / ns[j]))
      hw      <- q_crit * se
      q_obs   <- diff_ij / se
      p_adj   <- 1 - stats::ptukey(abs(q_obs), a, df_e)
      rows[[idx]] <- data.frame(
        contrast        = paste(lev[j], "-", lev[i]),
        mean_difference = diff_ij,
        se              = se,
        q_statistic     = q_obs,
        lower_limit     = diff_ij - hw,
        upper_limit     = diff_ij + hw,
        p_adjusted      = p_adj,
        stringsAsFactors = FALSE
      )
    }
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out <- .as_dmar_tbl(out, conf_level = conf_level)
  class(out) <- c("dmar_post_hoc_ci", class(out))
  out
}
