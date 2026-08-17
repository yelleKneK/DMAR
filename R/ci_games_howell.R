# Games-Howell simultaneous confidence intervals for all pairwise comparisons
# when homogeneity of variance is not assumed.
#' Provides Games--Howell Simultaneous Confidence Intervals for All Pairwise Comparisons Without Assuming Homogeneity of Variance
#'
#' @param x Either (a) a fitted \code{\link[stats]{lm}} or
#'   \code{\link[stats]{aov}} object with a single factor predictor, or
#'   (b) a numeric vector of the outcome, in which case \code{group} must
#'   also be supplied.
#' @param group When \code{x} is a vector, a factor (or coercible to
#'   factor) giving the group membership of each observation.
#' @param conf_level Family-wise confidence level. Default \code{0.95}.
#'
#' @return A \code{data.frame} with one row per pairwise comparison and
#'   columns \code{contrast}, \code{mean_difference}, \code{se},
#'   \code{df}, \code{q_statistic}, \code{lower_limit}, \code{upper_limit},
#'   and \code{p_adjusted}. The \code{df} column is the Welch--Satterthwaite
#'   degrees of freedom for that pair, which is why it varies from row to
#'   row. The table prints through the \code{\link{dmar_tbl}} display layer
#'   and works with \code{\link[generics]{tidy}} and
#'   \code{\link[generics]{glance}} (see \code{\link{dmar_tidiers}}).
#'
#' @details
#' Tukey's HSD and the Kramer modification for unequal \emph{n}
#' (\code{\link{ci_tukey_kramer}}) both pool the within-group variances into
#' \eqn{\mathit{MS}_W}, so both assume homogeneity of variance. Neither is
#' robust when that assumption fails. The Games--Howell procedure drops the
#' assumption: it uses a separate error term for each pair and a
#' Welch--Satterthwaite degrees of freedom for each pair, then takes its
#' critical value from the studentized range.
#'
#' For groups \eqn{g} and \eqn{h}, the standard error of the difference uses
#' only those two groups' variances, and the degrees of freedom are
#' \deqn{\mathit{df} = \frac{(s_g^2/n_g + s_h^2/n_h)^2}{s_g^4/[n_g^2(n_g-1)] + s_h^4/[n_h^2(n_h-1)]},}
#' the same Welch--Satterthwaite expression that base R's
#' \code{\link[stats]{t.test}} uses by default for two groups. A pair is
#' declared different when the observed \emph{t} exceeds \eqn{q/\sqrt{2}},
#' with \eqn{q} the studentized range critical value
#' (\code{\link{cv_tukey_hsd}}) at the pair's degrees of freedom, so the
#' interval for the difference of means is
#' \eqn{(\bar Y_g - \bar Y_h) \pm q_{\alpha;a,\mathit{df}}\sqrt{(s_g^2/n_g + s_h^2/n_h)/2}}.
#' Maxwell, Delaney, and Kelley (2027, Chapter 5) develop this as one of the
#' two modifications of Tukey's HSD for heterogeneous variances (their
#' Equations 5.13 and 5.14).
#'
#' \strong{When to use it.} Reach for Games--Howell when the group variances
#' are not interchangeable and the design is between subjects. It is the
#' heterogeneity-robust counterpart of \code{\link{ci_tukey_kramer}} and, like
#' it, controls the family-wise error rate across all \eqn{a(a-1)/2} pairs.
#' It handles unequal \emph{n} as a matter of course, so it does not need a
#' separate unequal-\emph{n} variant.
#'
#' \strong{When something else is better.} Dunnett (1980) found that
#' Games--Howell becomes slightly liberal (the family-wise error rate runs
#' somewhat above the nominal level) when the samples are small. Maxwell,
#' Delaney, and Kelley (2027, Chapter 5) therefore recommend Games--Howell for
#' larger samples and Dunnett's T3, which takes its critical value from the
#' studentized maximum modulus (\code{\link{cv_smm}}) rather than the
#' studentized range, when the groups have fewer than roughly 50 observations
#' each. When the variances are in fact homogeneous, use
#' \code{\link{ci_tukey_kramer}} instead: it pools the variances, so it has
#' more error degrees of freedom and more power. When only treatments are
#' compared to a single control, use \code{\link{ci_dunnett}}.
#'
#' With \eqn{a = 2} groups the procedure is exactly Welch's \emph{t} test:
#' the interval and the \emph{p}-value equal those from
#' \code{t.test(..., var.equal = FALSE)}, because
#' \eqn{q_{\alpha;2,\mathit{df}} = \sqrt{2}\,t_{1-\alpha/2,\mathit{df}}}.
#'
#' @references
#' Games, P. A., & Howell, J. F. (1976). Pairwise multiple comparison
#'   procedures with unequal \emph{n}'s and/or variances: A Monte Carlo
#'   study. \emph{Journal of Educational Statistics, 1}(2), 113--125.
#'   \doi{10.2307/1164979}
#'
#' Dunnett, C. W. (1980). Pairwise multiple comparisons in the unequal
#'   variance case. \emph{Journal of the American Statistical Association,
#'   75}(372), 796--800. \doi{10.2307/2287161}
#'
#' Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). \emph{Designing
#'   experiments and analyzing data: A model comparison perspective}
#'   (4th ed.). Routledge. (See Chapter 5 on the multiple-comparisons
#'   problem, where the modifications of Tukey's HSD for unequal \emph{n}
#'   and unequal variances are developed.)
#'
#' @examples
#' # The raw drinks_per_week outcome of the drinks_trial data: the
#' # Standard arm's variance is more than three times that of either CRA
#' # arm, so the pooled error term Tukey's HSD relies on is questionable.
#' ci_games_howell(drinks_trial$drinks_per_week, drinks_trial$treatment)
#'
#' # A fitted one-way model may be passed instead of the two vectors.
#' ci_games_howell(aov(drinks_per_week ~ treatment, data = drinks_trial))
#'
#' # With two groups the procedure is Welch's t test, so the limits agree;
#' # the trial's two enrollment cohorts give a two-group comparison.
#' ci_games_howell(drinks_trial$drinks_per_week, drinks_trial$cohort)$lower_limit
#' -t.test(drinks_per_week ~ cohort, data = drinks_trial)$conf.int[2]
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso \code{\link{ci_tukey_kramer}} for the homogeneity-assuming
#'   counterpart, \code{\link{ci_dunnett}} for many-to-one comparisons,
#'   \code{\link{ci_scheffe}} for arbitrary contrasts,
#'   \code{\link{cv_tukey_hsd}} for the critical value it uses, and
#'   \code{\link{dmar_tidiers}} for the tidy methods.
#'
#' @keywords htest design
#'
#' @export
#' @import stats
ci_games_howell <- function(x, group = NULL, conf_level = 0.95) {
  if (conf_level <= 0 || conf_level >= 1)
    stop("'conf_level' must be in (0, 1).")

  if (inherits(x, c("lm", "aov"))) {
    mf <- stats::model.frame(x)
    if (ncol(mf) != 2L)
      stop("ci_games_howell() requires a fitted lm/aov with one factor predictor.")
    y <- mf[[1L]]; g <- mf[[2L]]
    if (!is.factor(g)) g <- factor(g)
  } else if (is.numeric(x) && !is.null(group)) {
    if (length(x) != length(group))
      stop("'x' and 'group' must be the same length.")
    y <- x; g <- factor(group)
  } else {
    stop("Supply either a fitted lm/aov, or a numeric vector 'x' with a 'group' factor.")
  }

  ok <- !is.na(y) & !is.na(g)
  y <- y[ok]; g <- droplevels(g[ok])

  lev   <- levels(g)
  a     <- length(lev)
  if (a < 2) stop("At least two groups are required.")
  means <- tapply(y, g, mean)
  vars  <- tapply(y, g, stats::var)
  ns    <- tapply(y, g, length)
  if (any(ns < 2))
    stop("Each group needs at least two observations to estimate a variance.")

  # The standard error of a pair uses only that pair's two variances, so it is
  # zero (and the Welch-Satterthwaite df an undefined 0/0) exactly when both
  # groups in the pair are essentially constant. Welch's t handles a single
  # constant group, so only a pair of two constant groups is degenerate. Refuse
  # those pairs by name rather than return NaN limits, in the same spirit as
  # base R's t.test(), which errors "data are essentially constant".
  zero_var <- !is.finite(vars) | vars <= 0
  if (any(zero_var)) {
    degenerate <- character(0)
    for (i in seq_len(a - 1L)) {
      for (j in (i + 1L):a) {
        if (zero_var[i] && zero_var[j])
          degenerate <- c(degenerate, paste(lev[j], "-", lev[i]))
      }
    }
    if (length(degenerate) > 0L)
      stop("The pairwise sampling variance is zero because both groups are essentially constant, so the standard error and degrees of freedom are undefined for the contrast(s): ",
           paste(degenerate, collapse = ", "), ".", call. = FALSE)
  }

  rows <- list(); idx <- 0L
  for (i in seq_len(a - 1L)) {
    for (j in (i + 1L):a) {
      idx <- idx + 1L
      # A separate error term per pair (Maxwell, Delaney, & Kelley, 2027,
      # Eq. 5.13) and Welch-Satterthwaite degrees of freedom (their Eq. 5.14).
      v_i <- vars[i] / ns[i]
      v_j <- vars[j] / ns[j]
      df_ij <- (v_i + v_j)^2 /
        (v_i^2 / (ns[i] - 1) + v_j^2 / (ns[j] - 1))
      # se is on the studentized-range scale, so the critical value multiplies
      # it directly (equivalently, the observed t is compared with q/sqrt(2)).
      se      <- sqrt((v_i + v_j) / 2)
      diff_ij <- as.numeric(means[j] - means[i])
      q_crit  <- stats::qtukey(conf_level, nmeans = a, df = df_ij)
      hw      <- q_crit * se
      q_obs   <- diff_ij / se
      p_adj   <- stats::ptukey(abs(q_obs), nmeans = a, df = df_ij,
                               lower.tail = FALSE)
      rows[[idx]] <- data.frame(
        contrast        = paste(lev[j], "-", lev[i]),
        mean_difference = diff_ij,
        se              = as.numeric(se),
        df              = as.numeric(df_ij),
        q_statistic     = as.numeric(q_obs),
        lower_limit     = diff_ij - as.numeric(hw),
        upper_limit     = diff_ij + as.numeric(hw),
        p_adjusted      = as.numeric(p_adj),
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
