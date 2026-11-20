# Dunnett's T3 simultaneous confidence intervals for all pairwise comparisons
# when homogeneity of variance is not assumed.
#' Provides Dunnett's T3 Simultaneous Confidence Intervals for All Pairwise Comparisons Without Assuming Homogeneity of Variance
#'
#' @description
#' Forms simultaneous confidence intervals and multiplicity-adjusted
#' \emph{p}-values for all \eqn{a(a-1)/2} pairwise differences among
#' \eqn{a} group means without assuming that the groups share a variance.
#' Each pair has its own standard error and Welch--Satterthwaite degrees of
#' freedom, and its critical value is that of the studentized maximum
#' modulus for the \eqn{a(a-1)/2} comparisons at those degrees of freedom
#' (Dunnett, 1980). This is the procedure SPSS reports as Dunnett's T3;
#' \code{\link{cv_smm}} gives its critical value.
#'
#' @param x Either (a) a fitted \code{\link[stats]{lm}} or
#'   \code{\link[stats]{aov}} object with a single factor predictor, or
#'   (b) a numeric vector of the outcome, in which case \code{group} must
#'   also be supplied. An infinite outcome value is an error, since it
#'   leaves its group's mean and variance undefined.
#' @param group When \code{x} is a vector, a factor (or coercible to
#'   factor) giving the group membership of each observation. Observations
#'   with a missing outcome or a missing group are dropped.
#' @param conf_level Family-wise confidence level. Default \code{0.95}.
#'
#' @return A \code{data.frame} with one row per pairwise comparison, in the
#'   order \code{\link{ci_games_howell}} uses: with the factor levels in
#'   their order, each later level is compared with each earlier one, and
#'   the \code{contrast} label (\code{"B - A"}) names the difference taken.
#'   The columns are \code{contrast}, \code{mean_difference}, \code{se} (the
#'   separate-variance standard error of the difference), \code{df} (the
#'   pair's Welch--Satterthwaite degrees of freedom, which is why it varies
#'   from row to row), \code{t_value} (\code{mean_difference} divided by
#'   \code{se}), \code{ci_lower}, \code{ci_upper}, and
#'   \code{p_adjusted}. The \code{conf_level} attribute records the
#'   family-wise level, and the \code{ci_adjustment} and
#'   \code{p_adjustment} attributes are \code{"dunnett_t3"}. The table
#'   prints through the \code{\link{dmar_tbl}} display layer and works with
#'   \code{\link[generics]{tidy}} and \code{\link[generics]{glance}} (see
#'   \code{\link{dmar_tidiers}}).
#'
#' @details
#' \strong{The procedure.} For groups \eqn{g} and \eqn{h}, with means
#' \eqn{\bar Y_g} and \eqn{\bar Y_h}, variances \eqn{s_g^2} and
#' \eqn{s_h^2}, and sizes \eqn{n_g} and \eqn{n_h}, the standard error of the
#' difference uses only those two groups' variances,
#' \deqn{\mathit{SE}_{gh} = \sqrt{s_g^2/n_g + s_h^2/n_h},}
#' and its degrees of freedom are the Welch--Satterthwaite value
#' \deqn{\mathit{df}_{gh} = \frac{(s_g^2/n_g + s_h^2/n_h)^2}{s_g^4/[n_g^2(n_g-1)] + s_h^4/[n_h^2(n_h-1)]},}
#' as in Maxwell, Delaney, and Kelley (2027, Equations 5.13 and 5.14) and in
#' base R's \code{\link[stats]{t.test}} for two groups. With
#' \eqn{C = a(a-1)/2} comparisons in the family and \eqn{\alpha} equal to one
#' minus \code{conf_level}, the interval for \eqn{\mu_h - \mu_g} is
#' \deqn{(\bar Y_h - \bar Y_g) \pm m_{\alpha;\,C,\,\mathit{df}_{gh}}\,\mathit{SE}_{gh},}
#' where \eqn{m_{\alpha;C,\mathit{df}}} is the upper \eqn{\alpha} critical
#' value of the studentized maximum modulus for \eqn{C} comparisons on
#' \eqn{\mathit{df}} degrees of freedom (Stoline & Ury, 1979), the value
#' \code{cv_smm(alpha_level = 1 - conf_level, df, n_comparisons = C)}
#' returns. The adjusted \emph{p}-value is the probability that a
#' studentized maximum modulus variate for \eqn{C} comparisons on the pair's
#' degrees of freedom exceeds the observed
#' \eqn{|t| = |\bar Y_h - \bar Y_g| / \mathit{SE}_{gh}}, so a pair's
#' \code{p_adjusted} falls below \eqn{\alpha} when its interval excludes
#' zero and not otherwise, to numerical precision: the critical value's
#' tail area matches \eqn{\alpha} to a relative \eqn{10^{-6}} or better.
#' Both are computed by numerical integration rather than simulation, and
#' the degrees of freedom are not rounded. As with
#' Games--Howell, the family-wise level holds approximately rather than
#' exactly, since each pair is referred to a distribution indexed by its own
#' estimated degrees of freedom.
#'
#' \strong{Dunnett's T3 and Tamhane's T2.} Tamhane (1979) compared
#' procedures for pairwise comparisons of means with unequal variances. The
#' one known as Tamhane's T2 uses the same standard errors and degrees of
#' freedom but takes its critical value from the inequality of Šidák
#' (1967), treating the \eqn{C} \emph{t} statistics as if they were
#' independent: each pair is tested at the level
#' \eqn{1 - (1 - \alpha)^{1/C}}, so its critical value is the
#' \eqn{(1 + (1 - \alpha)^{1/C})/2} quantile of \emph{t} on
#' \eqn{\mathit{df}_{gh}} and its adjusted \emph{p}-value is
#' \eqn{1 - (1 - p_{gh})^C}, with \eqn{p_{gh}} the pair's two-sided Welch
#' \emph{p}-value. Dunnett (1980) replaced that bound with the studentized
#' maximum modulus at the same degrees of freedom. The two distribution
#' functions differ only in where the power is taken. With
#' \eqn{S = \sqrt{\chi^2_{\mathit{df}}/\mathit{df}}} and \eqn{\Phi} the
#' standard normal distribution function, the probability that all \eqn{C}
#' statistics fall within \eqn{\pm c} is
#' \eqn{\{\mathrm{E}[2\Phi(cS) - 1]\}^C} as T2 computes it and
#' \eqn{\mathrm{E}\{[2\Phi(cS) - 1]^C\}} under the studentized maximum
#' modulus, which by Jensen's inequality is never the smaller of the two.
#' T3's critical value therefore never exceeds T2's at the same degrees of
#' freedom: its intervals are never wider, and its adjusted
#' \emph{p}-values never larger. On the reactivated-story scores of
#' Maxwell, Delaney, and Kelley (2027, Chapter 12, Exercise 25), three
#' groups of 13, the T2 \emph{p}-values are .045, .001, and .027 and the T3
#' values .044, .001, and .026, which are the values SPSS reports for
#' Dunnett's T3. SPSS offers the two as separate procedures, Tamhane's T2
#' and Dunnett's T3; the T in T3 comes from Tamhane, who developed its
#' predecessor (Maxwell, Delaney, & Kelley, 2027, Chapter 5). A formula
#' that raises the single \emph{t} probability \eqn{2F_t(|t|) - 1}, with
#' \eqn{F_t} the distribution function of \emph{t} on
#' \eqn{\mathit{df}_{gh}}, to the power \eqn{C} computes T2, whatever it is
#' called.
#'
#' \strong{Dunnett's T3 or Games--Howell.} \code{\link{ci_games_howell}}
#' uses the same standard errors and degrees of freedom but refers
#' \eqn{|t|\sqrt{2}} to the studentized range for \eqn{a} means. Dunnett
#' (1980) found that Games--Howell becomes slightly liberal when the samples
#' are small, and Maxwell, Delaney, and Kelley (2027, Chapter 5) therefore
#' recommend Dunnett's T3 when the groups have fewer than roughly 50
#' observations each and Games--Howell for larger samples. When the
#' variances are in fact homogeneous, \code{\link{ci_tukey_kramer}} pools
#' them and so has more error degrees of freedom and more power. With two
#' groups the T3 interval and \emph{p}-value are exactly those of Welch's
#' \emph{t} test, \code{t.test(..., var.equal = FALSE)}, because the
#' studentized maximum modulus for one comparison is the two-sided \emph{t}.
#'
#' @references
#' Dunnett, C. W. (1980). Pairwise multiple comparisons in the unequal
#'   variance case. \emph{Journal of the American Statistical Association,
#'   75}(372), 796--800. \doi{10.2307/2287161}
#'
#' Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). \emph{Designing
#'   experiments and analyzing data: A model comparison perspective}
#'   (4th ed.). Routledge. (See Chapter 5 on the multiple-comparisons
#'   problem, where the modifications of Tukey's HSD for unequal variances,
#'   Dunnett's T3 and Games--Howell, are developed.)
#'
#' Šidák, Z. (1967). Rectangular confidence regions for the means of
#'   multivariate normal distributions. \emph{Journal of the American
#'   Statistical Association, 62}(318), 626--633.
#'   \doi{10.1080/01621459.1967.10482935}
#'
#' Stoline, M. R., & Ury, H. K. (1979). Tables of the Studentized maximum
#'   modulus distribution and an application to multiple comparisons among
#'   means. \emph{Technometrics, 21}(1), 87--93.
#'   \doi{10.1080/00401706.1979.10489726}
#'
#' Tamhane, A. C. (1979). A comparison of procedures for multiple
#'   comparisons of means with unequal variances. \emph{Journal of the
#'   American Statistical Association, 74}(366), 471--480.
#'   \doi{10.2307/2286358}
#'
#' @examples
#' # The raw drinks_per_week outcome of the drinks_trial data: the arms hold
#' # 37, 32, and 19 patients, fewer than the roughly 50 per group below
#' # which Games-Howell runs slightly liberal, and the Standard arm's
#' # variance is more than three times that of either CRA arm.
#' ci_dunnett_t3(drinks_trial$drinks_per_week, drinks_trial$treatment)
#'
#' # A fitted one-way model may be passed instead of the two vectors.
#' ci_dunnett_t3(aov(drinks_per_week ~ treatment, data = drinks_trial))
#'
#' # Each half-width is the studentized maximum modulus critical value for
#' # the three comparisons, at that pair's degrees of freedom, times the
#' # pair's standard error.
#' t3 <- ci_dunnett_t3(drinks_trial$drinks_per_week, drinks_trial$treatment)
#' t3$ci_upper[1] - t3$mean_difference[1]
#' cv_smm(alpha_level = .05, df = t3$df[1], n_comparisons = 3)$value * t3$se[1]
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso \code{\link{cv_smm}} for the critical value it uses,
#'   \code{\link{ci_games_howell}} for the studentized range counterpart,
#'   \code{\link{ci_tukey_kramer}} for the homogeneity-assuming counterpart,
#'   \code{\link{ci_dunnett}} for many-to-one comparisons, and
#'   \code{\link{dmar_tidiers}} for the tidy methods.
#'
#' @keywords htest design
#'
#' @family hypothesis tests
#'
#' @export
#' @import stats
ci_dunnett_t3 <- function(x, group = NULL, conf_level = 0.95) {
  if (!is.numeric(conf_level) || length(conf_level) != 1L ||
      is.na(conf_level) || conf_level <= 0 || conf_level >= 1)
    .dmar_refuse_argument("'conf_level' must be in (0, 1).", arg = "conf_level")

  if (inherits(x, c("lm", "aov"))) {
    mf <- stats::model.frame(x)
    if (ncol(mf) != 2L)
      .dmar_refuse_argument(
        "ci_dunnett_t3() requires a fitted lm/aov with one factor predictor.",
        arg = NULL)
    y <- mf[[1L]]; g <- mf[[2L]]
    if (!is.factor(g)) g <- factor(g)
  } else if (is.numeric(x) && !is.null(group)) {
    if (length(x) != length(group))
      .dmar_refuse_argument(
        "'x' and 'group' must be the same length.", arg = c("x", "group"))
    y <- x; g <- factor(group)
  } else {
    .dmar_refuse_argument(
      "Supply either a fitted lm/aov, or a numeric vector 'x' with a 'group' factor.",
      arg = c("x", "group"))
  }

  ok <- !is.na(y) & !is.na(g)
  y <- y[ok]; g <- droplevels(g[ok])
  # A missing value (NA or NaN) was dropped above. An infinite one leaves its
  # group's mean and variance undefined, and its NaN variance would otherwise
  # pass for a constant group below, so it is refused.
  if (any(!is.finite(y)))
    .dmar_refuse_argument(
      "The outcome contains infinite values, which leave a group's mean and variance undefined; remove or recode them first.",
      arg = NULL)

  lev <- levels(g)
  a   <- length(lev)
  if (a < 2) .dmar_refuse_argument(
    "At least two groups are required.", arg = NULL)
  means <- tapply(y, g, mean)
  vars  <- tapply(y, g, stats::var)
  ns    <- tapply(y, g, length)
  if (any(ns < 2))
    .dmar_refuse_argument(
      "Each group needs at least two observations to estimate a variance.",
      arg = NULL)

  # The standard error of a pair uses only that pair's two variances, so it is
  # zero (and the Welch-Satterthwaite df an undefined 0/0) exactly when both
  # groups in the pair are essentially constant. Welch's t handles a single
  # constant group, so only a pair of two constant groups is degenerate. Refuse
  # those pairs by name rather than return NaN limits, as ci_games_howell()
  # does and in the spirit of base R's t.test(), which errors "data are
  # essentially constant".
  zero_var <- vars <= 0
  if (any(zero_var)) {
    degenerate <- character(0)
    for (i in seq_len(a - 1L)) {
      for (j in (i + 1L):a) {
        if (zero_var[i] && zero_var[j])
          degenerate <- c(degenerate, paste(lev[j], "-", lev[i]))
      }
    }
    if (length(degenerate) > 0L)
      .dmar_refuse_argument(c(
        "The pairwise sampling variance is zero because both groups are essentially constant, so the standard error and degrees of freedom are undefined for the contrast(s): ",
        paste(degenerate, collapse = ", "), "."), arg = NULL)
  }

  # Every pair belongs to the one family of C comparisons, whatever its df.
  n_comparisons <- a * (a - 1L) / 2L
  alpha_level   <- 1 - conf_level

  rows <- vector("list", n_comparisons); idx <- 0L
  for (i in seq_len(a - 1L)) {
    for (j in (i + 1L):a) {
      idx <- idx + 1L
      # A separate error term per pair (Maxwell, Delaney, & Kelley, 2027,
      # Eq. 5.13) and Welch-Satterthwaite degrees of freedom (their Eq. 5.14).
      v_i <- as.numeric(vars[i] / ns[i])
      v_j <- as.numeric(vars[j] / ns[j])
      df_ij <- (v_i + v_j)^2 /
        (v_i^2 / (ns[[i]] - 1) + v_j^2 / (ns[[j]] - 1))
      se      <- sqrt(v_i + v_j)
      diff_ij <- as.numeric(means[j] - means[i])
      t_obs   <- diff_ij / se
      # The studentized maximum modulus for the C comparisons at this pair's
      # df (the internals cv_smm() uses), for both the limits and p-value.
      m_crit  <- .smm_upper_cv(alpha_level, df_ij, n_comparisons)
      p_adj   <- .smm_cdf(abs(t_obs), df_ij, n_comparisons, lower_tail = FALSE)
      rows[[idx]] <- data.frame(
        contrast        = paste(lev[j], "-", lev[i]),
        mean_difference = diff_ij,
        se              = se,
        df              = df_ij,
        t_value     = t_obs,
        ci_lower     = diff_ij - m_crit * se,
        ci_upper     = diff_ij + m_crit * se,
        p_adjusted      = p_adj,
        stringsAsFactors = FALSE
      )
    }
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out <- .as_dmar_tbl(out, conf_level = conf_level)
  attr(out, "ci_adjustment") <- "dunnett_t3"
  attr(out, "p_adjustment") <- "dunnett_t3"
  class(out) <- c("dmar_post_hoc_ci", class(out))
  out
}
