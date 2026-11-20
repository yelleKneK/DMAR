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
#'   (the studentized-range \eqn{q}), \code{ci_lower},
#'   \code{ci_upper}, \code{p_adjusted}.
#'
#' @details
#' \strong{Formula.} For groups \eqn{i, j} with means \eqn{\bar y_i,
#' \bar y_j} and sample sizes \eqn{n_i, n_j}, the Tukey-Kramer
#' simultaneous CI is
#' \deqn{\bar y_i - \bar y_j \;\pm\; q_{\alpha, a, \nu}
#'   \sqrt{\frac{\mathit{MS}_E}{2} \left(\frac{1}{n_i} + \frac{1}{n_j}\right)},}
#' where \eqn{q_{\alpha, a, \nu}} is the upper \eqn{\alpha} quantile of
#' the studentized-range distribution with \eqn{a} groups and \eqn{\nu}
#' error degrees of freedom. The quantile is base R's
#' \code{\link[stats]{qtukey}} where that is accurate and the exact quantile
#' where it is not: \code{qtukey} and \code{ptukey} lose accuracy at small
#' error df (fewer than 9 for up to 6 groups, rising by steps to fewer than 111
#' with more than 3000 groups), above 10000 and up to \eqn{10^8} df, at confidence
#' levels above .999, and, with more than 50 groups, below the .75 point of
#' the distribution (.80 with more than 1000 groups; with 20 groups and 2
#' error df, \code{qtukey(.99, 20, 2)}
#' is 44.526 where the exact value is 37.943), and \code{qtukey} can stop
#' short of its root anywhere. There the distribution function is computed by
#' integrating the probability integral of the range against the chi square
#' distribution of the error variance, and the quantile is its root.
#' \code{\link{cv_tukey_hsd}} describes the computation, the regions, and how
#' they were measured.
#'
#' \strong{Why Tukey-Kramer.} Hayter (1984) proved that the
#' Tukey-Kramer procedure is conservative for unbalanced designs (the
#' coverage probability is at least \code{conf_level}). For balanced
#' designs it reduces to Tukey's HSD and the coverage is exactly
#' \code{conf_level}.
#'
#' \strong{Adjusted \emph{p}-values.} Each pairwise \eqn{p}-value is
#' computed from the studentized-range distribution,
#' \eqn{p = P(Q_{a, \nu} > |q|)}, by the same computation as the quantile
#' (\code{\link[stats]{ptukey}} where it is accurate, the exact distribution
#' function where it is not).
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
  .dmar_check_conf_level(conf_level)
  if (conf_level <= 0 || conf_level >= 1)
    .dmar_refuse_argument("'conf_level' must be in (0, 1).", arg = "conf_level")

  if (inherits(x, c("lm", "aov"))) {
    mf <- stats::model.frame(x)
    if (ncol(mf) != 2L)
      .dmar_refuse_argument(
        "ci_tukey_kramer() requires a fitted lm/aov with one factor predictor.",
        arg = NULL)
    y <- mf[[1L]]; g <- mf[[2L]]
    if (!is.factor(g)) g <- factor(g)
    tab <- stats::anova(x); rn <- trimws(rownames(tab))
    ms_e <- tab[rn == "Residuals", "Mean Sq"]
    df_e <- tab[rn == "Residuals", "Df"]
  } else if (is.numeric(x) && !is.null(group)) {
    if (length(x) != length(group))
      .dmar_refuse_argument(
        "'x' and 'group' must be the same length.", arg = c("x", "group"))
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
    .dmar_refuse_argument(
      "Supply either a fitted lm/aov, or (x, group) vectors.", arg = NULL)
  }

  lev <- levels(g)
  a <- length(lev); m <- a * (a - 1L) / 2L
  if (a < 2L) .dmar_refuse_argument("Need at least 2 groups.", arg = NULL)
  means <- tapply(y, g, mean)
  ns    <- tapply(y, g, length)

  # stats::qtukey() where it is accurate, the exact quantile elsewhere
  # (R/studentized_range_internals.R); the p-values below likewise.
  q_crit <- .sr_quantile(conf_level, a, df_e)

  rows <- vector("list", m); idx <- 0L
  for (i in seq_len(a - 1L)) {
    for (j in (i + 1L):a) {
      idx <- idx + 1L
      diff_ij <- as.numeric(means[j] - means[i])
      se      <- sqrt(ms_e / 2 * (1 / ns[i] + 1 / ns[j]))
      hw      <- q_crit * se
      q_obs   <- diff_ij / se
      p_adj   <- 1 - .sr_cdf(abs(q_obs), a, df_e)
      rows[[idx]] <- data.frame(
        contrast        = paste(lev[j], "-", lev[i]),
        mean_difference = diff_ij,
        se              = se,
        q_statistic     = q_obs,
        ci_lower     = diff_ij - hw,
        ci_upper     = diff_ij + hw,
        p_adjusted      = p_adj,
        stringsAsFactors = FALSE
      )
    }
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out <- .as_dmar_tbl(out, conf_level = conf_level)
  attr(out, "ci_adjustment") <- "tukey"
  attr(out, "p_adjustment") <- "tukey"
  class(out) <- c("dmar_post_hoc_ci", class(out))
  out
}
