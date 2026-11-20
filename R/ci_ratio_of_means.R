#' Confidence Interval for the Ratio of Two Independent Means
#'
#' Computes a confidence set for the ratio of two population means
#' estimated from independent groups, \eqn{\rho = \mu_1 / \mu_2}, or, with
#' a third, baseline group, for the ratio of the two groups' advantages
#' over the baseline, \eqn{\rho = (\mu_1 - \mu_b) / (\mu_2 - \mu_b)}. A
#' ratio states an effect in relative terms (a new procedure's mean is 90
#' percent of the standard procedure's), which is the natural scale when
#' the outcome has a meaningful zero. The second form is the fraction of a
#' reference treatment's advantage over a placebo or a trivial baseline
#' that a new treatment retains, the retention of effect of the three-arm
#' noninferiority trial (Pigeot et al., 2003). The default interval is
#' Fieller's, which inverts the \emph{t} test of
#' \eqn{\mu_1 - \rho \mu_2 = 0} and is reported as what it is: a bounded
#' interval, the complement of an interval, or the whole line. The
#' log-scale interval of the response ratio and percentile and BCa
#' bootstrap intervals are offered beside it. For two means measured on
#' the same units, see \code{\link{ci_ratio_of_means_paired}}.
#'
#' @param x,y Numeric vectors of observations from the numerator group and
#'   the denominator group. The groups are independent and need not be the
#'   same size; \code{NA}s are removed from each vector separately.
#' @param baseline An optional numeric vector of observations from a third,
#'   independent baseline group (a placebo arm, a trivial comparison
#'   condition). When it is supplied the ratio is that of the advantages,
#'   \eqn{(\mu_1 - \mu_b) / (\mu_2 - \mu_b)}.
#' @param mean_1,mean_2 Group sample means, for analysis from summary
#'   statistics in place of \code{x} and \code{y}.
#' @param sd_1,sd_2 Group sample standard deviations (with \eqn{n - 1} in
#'   the denominator of the variance).
#' @param n_1,n_2 Group sample sizes.
#' @param mean_baseline,sd_baseline,n_baseline The baseline group's
#'   sample mean, standard deviation, and size, for the ratio of the two
#'   groups' advantages over it from summary statistics (the three arms a
#'   noninferiority trial reports). Supply all three or none.
#' @param var_equal Logical. \code{FALSE} (the default) lets the group
#'   variances differ and inverts the Welch-type test, whose degrees of
#'   freedom change with the hypothesized ratio; \code{TRUE} pools the
#'   variances, which gives Fieller's exact interval under normality and
#'   equal variances. Used by \code{"fieller"} and \code{"log"}; the
#'   bootstrap resamples each group from itself and needs neither.
#' @param ci_method The interval: \code{"fieller"} (the default),
#'   \code{"log"} (the response ratio interval, formed on the log scale
#'   and defined only when both means, or both advantages, are positive),
#'   \code{"percentile"}, or \code{"bca"} (bootstrap intervals, which need
#'   the raw data).
#' @param conf_level Confidence level (default \code{0.95}).
#' @param alpha_lower The tail area at the lower limit. An omitted
#'   confidence level yields to supplied tail areas; an explicit level must
#'   agree with them. Each tail area is at most .5.
#' @param alpha_upper The tail area at the upper limit, supplied together
#'   with \code{alpha_lower}. The table records the coverage,
#'   \code{1 - alpha_lower - alpha_upper}, as its level and, when the two
#'   tail areas differ, the areas themselves, which the line beneath the
#'   printed table, \code{\link{results_sentence}}, and the
#'   \code{alpha_lower} and \code{alpha_upper} columns of \code{tidy()}
#'   state; with a tail area of 0 the interval is a one-sided confidence
#'   bound, unbounded on that side.
#' @param B Number of bootstrap replications for \code{"percentile"} and
#'   \code{"bca"}, a single whole number (ignored otherwise). Default
#'   \code{10000}: a replication of a ratio of means costs little, so a
#'   reported interval deserves the default, and the BCa adjustment moves
#'   the working quantiles into the tails, where a small \code{B} is least
#'   stable. The fewest accepted with 0.025 in each tail is 400 for the
#'   percentile interval and 2000 for BCa (200 and 1000 at a 0.90 level,
#'   2000 and 10000 at 0.99); with unequal tail areas the smaller one
#'   decides. A smaller \code{B} is refused before anything is computed
#'   (class \code{dmar_count_invalid}). See \emph{Bootstrap replications}
#'   under \emph{Details} and \code{\link{dmar_bootstrap}}.
#' @param seed Optional integer seed for the bootstrap. The default
#'   \code{NULL} uses the current state of the random number generator; a
#'   supplied seed is set for the call and the prior state restored on
#'   exit.
#'
#' @details
#' \strong{Fieller's confidence set.} For a candidate ratio \eqn{\rho} the
#' contrast \eqn{\bar x_1 - \rho \bar x_2} (with a baseline,
#' \eqn{(\bar x_1 - \bar x_b) - \rho (\bar x_2 - \bar x_b)}) estimates
#' zero exactly when \eqn{\rho} is the population ratio, and its \emph{t}
#' statistic,
#' \deqn{T(\rho) = \frac{\hat\theta_1 - \rho \hat\theta_2}
#'   {\sqrt{v_{11} - 2 \rho v_{12} + \rho^2 v_{22}}},}{T(rho) =
#'   (theta1 - rho theta2) / sqrt(v11 - 2 rho v12 + rho^2 v22),}
#' where \eqn{\hat\theta_1} and \eqn{\hat\theta_2} are the numerator and
#' denominator and the \eqn{v} terms their estimated variances and
#' covariance, is a pivot: it involves the data and \eqn{\rho} alone
#' (Fieller, 1941, 1954). The 1941 reference is Fieller's study of the
#' biological standardization of insulin, which Fieller (1954) cites with
#' the year 1940. The confidence set is every \eqn{\rho} that the
#' test does not reject, the solution of a quadratic inequality in
#' \eqn{\rho}. Its form is decided by the denominator's \emph{t},
#' \eqn{\hat\theta_2 / \sqrt{v_{22}}}, reported in the \code{t_denominator}
#' row (Fieller, 1954, p. 177). The row keeps the sign of the denominator,
#' and only its absolute value matters: when the absolute value of the
#' denominator's \emph{t} exceeds the critical value (with unequal tail
#' areas, each side of the estimate is bounded when that absolute value
#' exceeds the side's own critical value) the set is a bounded interval.
#' When it does not, the set contains arbitrarily large values of the
#' ratio, as it should when the data cannot distinguish the denominator
#' from zero: it is then either the complement of an interval, every value
#' outside a gap, or, when the data are consistent with every ratio at the
#' chosen level, the whole line. The table reports all three as they are.
#' The \code{ci_lower} and \code{ci_upper} rows are the outer ends of
#' the set (\code{-Inf} and \code{Inf} when it is unbounded), and the
#' \code{excluded_lower} and \code{excluded_upper} rows are the ends of the
#' gap of a complement (\code{NA} otherwise): the set is every value from
#' \code{ci_lower} to \code{ci_upper} except those strictly between
#' \code{excluded_lower} and \code{excluded_upper}. The
#' \code{confidence_set} attribute names the form, and so does the line
#' printed under the table; \code{\link{results_sentence}} writes the set in
#' that form, for example as every value except those between the ends of
#' the gap. With unequal tail areas one side of the set can be bounded while
#' the other is not, which can also give an interval unbounded on one side
#' or the union of two intervals. von Luxburg and Franz (2009) derive the
#' same sets geometrically, from confidence intervals for the means of
#' linear combinations of the two variables.
#'
#' \strong{Equal or unequal variances.} With \code{var_equal = TRUE} the
#' variances are pooled over the groups, on \eqn{n_1 + n_2 - 2} degrees of
#' freedom (\eqn{n_1 + n_2 + n_b - 3} with a baseline), and the set is
#' Fieller's exact one under normality and equal variances. With the
#' default \code{var_equal = FALSE} each group keeps its own variance and
#' the contrast is tested as Welch's \emph{t} tests a difference, with
#' Satterthwaite degrees of freedom; those degrees of freedom depend on
#' \eqn{\rho}, so the set has no closed form and is found by solving for
#' the values of \eqn{\rho} at which the test's \emph{p}-value equals its
#' tail area. The interval therefore agrees with the test: without a
#' baseline and with equal tail areas, each finite end \eqn{L} of the set
#' (a limit or an end of the gap) is a value at which
#' \code{t.test(x, L * y)} has \emph{p}-value \code{1 - conf_level}; with a
#' baseline the test is that of the contrast
#' \eqn{\bar x_1 - L \bar x_2 - (1 - L) \bar x_b}, with the Satterthwaite
#' degrees of freedom of its three variances. The degrees of freedom of the
#' test at each finite end of the set are in the \code{df_limits}
#' attribute, and the \code{df} row holds those of the denominator's
#' \emph{t}. This inversion is the one the \pkg{mratios} function
#' \code{ttestratio()} performs with \code{iterativeCI = TRUE}, which its
#' documentation describes as searching iteratively for the limits by
#' inverting the test; its default instead puts degrees of freedom
#' computed at the estimate into Fieller's formula, and its documentation
#' notes that those limits need not agree with the test. Hasler, Vonk, and
#' Hothorn (2008) treat the three-arm design when the group variances
#' differ.
#'
#' \strong{The log ratio.} \code{ci_method = "log"} forms the interval for
#' \eqn{\log \rho} from the delta method standard error,
#' \eqn{\sqrt{v_{11}/\hat\theta_1^2 + v_{22}/\hat\theta_2^2 -
#' 2 v_{12}/(\hat\theta_1 \hat\theta_2)}}, and maps it back by the
#' exponential. This is the response ratio of meta-analysis (Hedges,
#' Gurevitch, & Curtis, 1999), offered so a ratio reported here can be
#' compared with one entered in a meta-analysis. It is defined only when
#' both terms of the ratio are positive, and is refused otherwise with a
#' condition of class \code{dmar_log_ratio_undefined} (which inherits
#' \code{dmar_refusal}). It is always bounded; when the denominator is well
#' determined it is close to Fieller's interval, and the two part company
#' as the denominator's \emph{t} falls. With unequal variances its degrees
#' of freedom are the Satterthwaite degrees of freedom of the linearized
#' log ratio.
#'
#' \strong{The bootstrap.} \code{"percentile"} and \code{"bca"} resample
#' each group from itself (the baseline group too), as \pkg{boot}, the
#' software of Davison and Hinkley (1997), does for several samples given
#' as strata, and compute the ratio in every replication. The BCa
#' interval reads its limits at quantile positions moved from those of the
#' percentile interval by a bias correction, the median bias of the
#' replications in normal units, and by an acceleration, the rate at which
#' the standard error of the estimate changes with the parameter on a
#' normalized scale (Efron & Tibshirani, 1993); its acceleration is
#' estimated by the jackknife, deleting one observation at a time within
#' its group, with each group's influence values scaled by its size. A
#' bootstrap interval for a ratio is always bounded, whatever the
#' denominator's \emph{t}, so it can undercover badly when the denominator
#' is near zero, exactly where Fieller's set becomes unbounded. When the
#' Fieller set at the same level is unbounded, the line under the table
#' says so. Pigeot et al. (2003)
#' compare the percentile interval with Fieller's in the three-arm design.
#' Summary statistics cannot be resampled, so a bootstrap from them is
#' refused with a condition of class \code{dmar_bootstrap_needs_data}.
#'
#' \strong{Bootstrap replications.} The limits are quantiles of the
#' \code{B} replicated ratios. A limit that cuts off a tail area
#' \eqn{\alpha} is located by the replications beyond it, so the percentile
#' interval is run only when \code{B} puts at least 10 replications beyond
#' each limit asked for (400 with 0.025 in each tail; with unequal tail
#' areas the smaller one decides, and a one-sided bound needs 10 beyond its
#' one limit), and the BCa interval only with 50 beyond each nominal limit
#' (2000 with 0.025 in each tail). A replication in which the ratio is
#' undefined (a denominator of exactly zero, which discrete data can
#' produce) is dropped and counted, never replaced, and the interval is
#' computed from the rest. The result records \code{B} (requested),
#' \code{B_used} (the number the interval rests on), \code{n_dropped}, and
#' \code{B_required} (the fewest replications a request may ask for; that
#' count is held to \code{B}, not to \code{B_used}) as attributes, and the
#' line under the table states them,
#' with \code{bootstrap_status} (\code{"complete"}, \code{"dropped"},
#' \code{"exploratory"}, \code{"insufficient_count"}, or
#' \code{"insufficient_share"}). One warning of class
#' \code{dmar_bootstrap_replications_dropped} states the counts when any
#' replication was dropped, with a caution when the dropped replications
#' are at least the 2.5 percent of \code{B} that lies beyond one limit at
#' the 95 percent level, because they could then decide by themselves where
#' that limit falls. When more than half were dropped, the limits asked for
#' are \code{NA}, the
#' estimate is returned, and the warning has class
#' \code{dmar_bootstrap_insufficient}. For the BCa interval the count of 10
#' is applied again where the bias correction and the acceleration place
#' the limits, and when fewer than 10 replications lie beyond one of them
#' no BCa interval is reported (limits \code{NA}, the same class, with the
#' two positions as \code{positions} and a suggested \code{B} as
#' \code{B_suggest}). A BCa interval whose bias correction is infinite
#' (every replicate on one side of the estimate) or whose jackknife is
#' undefined (deleting one observation leaves a denominator of zero) stops
#' with an error that names \code{ci_method = "percentile"}.
#' \code{options(dmar.bootstrap_exploratory = TRUE)} lets a smaller
#' \code{B} (2 or more) run for a first look or a demonstration; every such
#' call gives a warning of class \code{dmar_bootstrap_exploratory} and the
#' result is not one to report. See \code{\link{dmar_bootstrap}} for the
#' rule and its reasons.
#'
#' \strong{Reading the interval for noninferiority.} With a baseline, a
#' new treatment is shown to retain more than a fraction \eqn{f} of the
#' reference treatment's advantage, at one-sided level \code{alpha_lower},
#' exactly when \eqn{f} is below the lower limit of a set that is a
#' bounded interval (or, for a one-sided bound, an interval unbounded only
#' above) and whose denominator is positive. The one-sided bound is
#' \code{alpha_lower = 0.05, alpha_upper = 0}, and a two-sided 95 percent
#' interval gives the bound at one-sided level .025. Reporting the limit
#' lets each reader apply his or her own fraction; the function itself
#' tests no fraction and states no verdict. When the denominator is
#' negative, or the set is unbounded below, the lower limit does not carry
#' that reading.
#'
#' @return A long \code{data.frame} (classes \code{dmar_ci_ratio_of_means},
#'   \code{dmar_ci_long}, \code{dmar_tbl}) with columns \code{term} and
#'   \code{value}, whose rows are \code{ci_lower}, \code{ratio}, and
#'   \code{ci_upper} (the set's outer ends, \code{-Inf} or \code{Inf}
#'   when it is unbounded on that side); \code{excluded_lower} and
#'   \code{excluded_upper} (the ends of the excluded gap of a complement or
#'   a union of two intervals, \code{NA} otherwise); \code{mean_1},
#'   \code{mean_2}, and \code{mean_baseline} (the group means;
#'   \code{mean_baseline} is \code{NA} without a baseline);
#'   \code{se_1}, \code{se_2}, and \code{se_baseline} (the standard errors of
#'   the means, from the pooled variance when \code{var_equal = TRUE});
#'   \code{t_denominator} (the denominator, \eqn{\bar x_2} or
#'   \eqn{\bar x_2 - \bar x_b}, over its standard error); \code{df} (the
#'   degrees of freedom of that \emph{t}: pooled, or Satterthwaite's for
#'   the denominator when \code{var_equal = FALSE}); and \code{n_1},
#'   \code{n_2}, and \code{n_baseline}. The rows are the same for every
#'   call. Attributes: \code{confidence_set} (\code{"interval"},
#'   \code{"complement of an interval"}, \code{"whole line"}, or, with
#'   unequal tail areas only, \code{"union of two intervals"};
#'   \code{NA} when a bootstrap interval was withheld and both limits are
#'   \code{NA}),
#'   \code{ci_method}, \code{var_equal}, \code{conf_level}, and the tail
#'   areas when they differ; \code{df_limits} when \code{var_equal = FALSE}
#'   and the interval is Fieller's or the log ratio's, the degrees of
#'   freedom of the test at each end of the set, a vector with elements
#'   \code{lower}, \code{upper}, \code{excluded_lower}, and
#'   \code{excluded_upper} that is \code{NA} at an end no test decides (an
#'   infinite limit, the open end of a one-sided bound, or a gap the set
#'   does not have); \code{B}, \code{B_used}, \code{n_dropped},
#'   \code{B_required}, and \code{bootstrap_status} for a bootstrap
#'   interval (see \emph{Bootstrap replications} under Details); and \code{fieller_set} (the form of the Fieller
#'   set) when the interval is a log or bootstrap interval and the Fieller
#'   set at the same level is unbounded. \code{tidy()} returns the ratio,
#'   its limits, the ends of any excluded gap, and the form of the set in
#'   one row (see \code{\link{dmar_tidiers}}).
#'
#' @references
#' Davison, A. C., & Hinkley, D. V. (1997). \emph{Bootstrap methods and
#'   their application}. Cambridge University Press.
#'   \doi{10.1017/CBO9780511802843}
#'
#' Efron, B., & Tibshirani, R. J. (1993). \emph{An introduction to the
#'   bootstrap}. New York, NY: Chapman & Hall/CRC.
#'
#' Fieller, E. C. (1941). The biological standardization of insulin.
#'   \emph{Supplement to the Journal of the Royal Statistical Society,
#'   7}(1), 1--54. \doi{10.2307/2983630}
#'
#' Fieller, E. C. (1954). Some problems in interval estimation.
#'   \emph{Journal of the Royal Statistical Society, Series B, 16}(2),
#'   175--185. \doi{10.1111/j.2517-6161.1954.tb00159.x}
#'
#' Hasler, M., Vonk, R., & Hothorn, L. A. (2008). Assessing
#'   non-inferiority of a new treatment in a three-arm trial in the
#'   presence of heteroscedasticity. \emph{Statistics in Medicine, 27}(4),
#'   490--503. \doi{10.1002/sim.3052}
#'
#' Hedges, L. V., Gurevitch, J., & Curtis, P. S. (1999). The meta-analysis
#'   of response ratios in experimental ecology. \emph{Ecology, 80}(4),
#'   1150--1156. \doi{10.1890/0012-9658(1999)080[1150:TMAORR]2.0.CO;2}
#'
#' Pigeot, I., Schäfer, J., Röhmel, J., & Hauschke, D. (2003). Assessing
#'   non-inferiority of a new treatment in a three-arm clinical trial
#'   including a placebo. \emph{Statistics in Medicine, 22}(6), 883--899.
#'   \doi{10.1002/sim.1450}
#'
#' von Luxburg, U., & Franz, V. H. (2009). A geometric approach to
#'   confidence sets for ratios: Fieller's theorem, generalizations and
#'   bootstrap. \emph{Statistica Sinica, 19}(3), 1095--1117.
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso \code{\link{ci_ratio_of_means_paired}} for two means measured
#'   on the same units; \code{\link{summary_t_test}} and
#'   \code{\link{welch_t}} for the difference between two means.
#'
#' @examples
#' # Drinks per week among the homeless alcohol-dependent adults of
#' # drinks_trial: the mean under the community reinforcement approach (CRA)
#' # relative to the mean under standard treatment. The group variances are
#' # not assumed equal.
#' drinks <- drinks_trial$drinks_per_week
#' cra <- drinks[drinks_trial$treatment == "CRA"]
#' standard <- drinks[drinks_trial$treatment == "Standard"]
#' ci_ratio_of_means(cra, standard)
#'
#' # The drinking counts are strongly skewed, so the response ratio interval
#' # and a BCa bootstrap interval are worth setting beside Fieller's;
#' # B = 2000 keeps the example quick, and a reported interval deserves the
#' # default B = 10000.
#' ci_ratio_of_means(cra, standard, ci_method = "log")
#' ci_ratio_of_means(cra, standard, ci_method = "bca", B = 2000, seed = 113)
#'
#' # From the summary statistics a paper reports, with pooled variances.
#' ci_ratio_of_means(mean_1 = 104, sd_1 = 14, n_1 = 40,
#'                   mean_2 = 96, sd_2 = 15, n_2 = 42, var_equal = TRUE)
#'
#' # A three-arm trial reported as summary statistics: the fraction of the
#' # reference arm's advantage over placebo that the new arm retains.
#' ci_ratio_of_means(mean_1 = 60.1, sd_1 = 12.2, n_1 = 60,
#'                   mean_2 = 62.4, sd_2 = 10.1, n_2 = 60,
#'                   mean_baseline = 50.3, sd_baseline = 9.8,
#'                   n_baseline = 60)
#'
#' # A three-arm trial: the fraction of the reference treatment's advantage
#' # over placebo that a new treatment retains, with a one-sided 95 percent
#' # lower bound. Noninferiority at a retained fraction f holds for every f
#' # below the lower limit.
#' set.seed(113)
#' placebo <- rnorm(60, mean = 50, sd = 10)
#' reference <- rnorm(60, mean = 62, sd = 10)
#' new_drug <- rnorm(60, mean = 60, sd = 12)
#' ci_ratio_of_means(new_drug, reference, baseline = placebo,
#'                   alpha_lower = 0.05, alpha_upper = 0)
#'
#' @keywords htest
#' @family confidence intervals for effect sizes
#'
#' @export
ci_ratio_of_means <- function(x = NULL, y = NULL, baseline = NULL,
                              mean_1 = NULL, sd_1 = NULL, n_1 = NULL,
                              mean_2 = NULL, sd_2 = NULL, n_2 = NULL,
                              mean_baseline = NULL, sd_baseline = NULL,
                              n_baseline = NULL, var_equal = FALSE,
                              ci_method = c("fieller", "log", "percentile",
                                            "bca"),
                              conf_level = 0.95, alpha_lower = NULL,
                              alpha_upper = NULL, B = 10000L, seed = NULL) {
  ci_method <- match.arg(ci_method)
  tails <- .ci_resolve_tails(conf_level, alpha_lower, alpha_upper,
                             level_missing = missing(conf_level),
                             max_tail = 0.5)
  conf_level <- tails$conf_level
  alpha_lower <- tails$alpha_lower
  alpha_upper <- tails$alpha_upper
  if (!is.logical(var_equal) || length(var_equal) != 1L || is.na(var_equal)) {
    .dmar_refuse_argument(
      "'var_equal' must be TRUE or FALSE.", arg = "var_equal")
  }
  bootstrap <- ci_method %in% c("percentile", "bca")

  stats_given <- !vapply(list(mean_1, sd_1, n_1, mean_2, sd_2, n_2),
                         is.null, logical(1L))
  base_given <- !vapply(list(mean_baseline, sd_baseline, n_baseline),
                        is.null, logical(1L))
  raw <- !is.null(x) || !is.null(y) || !is.null(baseline)
  if (raw && (any(stats_given) || any(base_given))) {
    .dmar_refuse_argument(c(
      "Supply either the raw data ('x', 'y', and optionally ",
      "'baseline') or the summary statistics ('mean_1', 'sd_1', 'n_1', ",
      "'mean_2', 'sd_2', 'n_2', and optionally 'mean_baseline', ",
      "'sd_baseline', 'n_baseline'), not both."), arg = c("x", "y", "baseline"))
  }
  if (raw) {
    if (is.null(x) || is.null(y)) {
      .dmar_refuse_argument(c(
        "Supply both 'x' (the numerator group) and 'y' (the ",
        "denominator group)."), arg = c("x", "y"))
    }
    groups <- list(x = x, y = y, baseline = baseline)
    groups <- groups[!vapply(groups, is.null, logical(1L))]
    for (g in names(groups)) {
      .ratio_check_data(groups[[g]], g)
      groups[[g]] <- groups[[g]][!is.na(groups[[g]])]
      if (length(groups[[g]]) < 2L) {
        .dmar_refuse_argument(c(
          "'", g, "' needs at least two observations that are not ",
          "missing."), arg = NULL)
      }
    }
    m <- unname(vapply(groups, mean, numeric(1L)))
    s <- unname(vapply(groups, stats::sd, numeric(1L)))
    n <- unname(vapply(groups, length, integer(1L)))
    if (any(s == 0)) {
      .dmar_refuse_argument(c(
        "'", names(groups)[s == 0][1L], "' has no variability; a ",
        "confidence interval needs a positive standard deviation in ",
        "every group."), arg = NULL)
    }
  } else {
    if (!any(stats_given) && !any(base_given)) {
      .dmar_refuse_argument(c(
        "Supply the raw data in 'x' and 'y', or the summary statistics ",
        "'mean_1', 'sd_1', 'n_1', 'mean_2', 'sd_2', and 'n_2'."),
        arg = c("x", "y", "mean_1"))
    }
    if (!all(stats_given)) {
      .dmar_refuse_argument(c(
        "From summary statistics, supply all six of 'mean_1', 'sd_1', ",
        "'n_1', 'mean_2', 'sd_2', and 'n_2'."),
        arg = c("mean_1", "sd_1", "n_1"))
    }
    if (any(base_given) && !all(base_given)) {
      .dmar_refuse_argument(c(
        "For a baseline group from summary statistics, supply all three ",
        "of 'mean_baseline', 'sd_baseline', and 'n_baseline'."),
        arg = c("mean_baseline", "sd_baseline", "n_baseline"))
    }
    if (bootstrap) .ratio_bootstrap_refuse(ci_method)
    .ratio_check_scalar(mean_1, "mean_1")
    .ratio_check_scalar(mean_2, "mean_2")
    .ratio_check_scalar(sd_1, "sd_1", positive = TRUE)
    .ratio_check_scalar(sd_2, "sd_2", positive = TRUE)
    .ratio_check_n(n_1, "n_1")
    .ratio_check_n(n_2, "n_2")
    m <- c(mean_1, mean_2)
    s <- c(sd_1, sd_2)
    n <- c(n_1, n_2)
    if (all(base_given)) {
      .ratio_check_scalar(mean_baseline, "mean_baseline")
      .ratio_check_scalar(sd_baseline, "sd_baseline", positive = TRUE)
      .ratio_check_n(n_baseline, "n_baseline")
      m <- c(m, mean_baseline)
      s <- c(s, sd_baseline)
      n <- c(n, n_baseline)
    }
  }

  # The means, the variances of the means, and the covariance matrix of
  # the numerator and denominator. A baseline mean enters both, so its
  # variance is their covariance.
  k <- length(m)
  nu <- n - 1
  df_pooled <- sum(nu)
  vm <- if (var_equal) sum(nu * s^2) / df_pooled / n else s^2 / n
  m_b <- if (k == 3L) m[3L] else 0
  v_b <- if (k == 3L) vm[3L] else 0
  num <- m[1L] - m_b
  den <- m[2L] - m_b
  if (den == 0) {
    .dmar_refuse_argument(c(
      "The denominator of the ratio is exactly zero, so the ratio is ",
      "undefined."), arg = NULL)
  }
  V <- matrix(c(vm[1L] + v_b, v_b, v_b, vm[2L] + v_b), 2L, 2L)
  df_den <- if (var_equal) {
    df_pooled
  } else if (k == 3L) {
    .satterthwaite_df(c(1, -1), vm[2:3], nu[2:3])
  } else {
    nu[2L]
  }
  ratio <- num / den
  attrs <- list(var_equal = var_equal)

  # Fieller's set at the requested tails, the reported interval for
  # "fieller" and the reference whose form a log or bootstrap interval is
  # read against.
  fieller <- if (var_equal) {
    .fieller_set(num, den, V, .ratio_crit(alpha_lower, df_pooled),
                 .ratio_crit(alpha_upper, df_pooled))
  } else {
    .fieller_set_df(m, vm, nu, alpha_lower, alpha_upper)
  }
  set <- switch(ci_method,
    fieller = {
      if (!var_equal) attrs$df_limits <- fieller$df_limits
      fieller
    },
    log = {
      if (num <= 0 || den <= 0) {
        .log_ratio_refuse(num, den, if (k == 3L) {
          "advantages over the baseline"
        } else {
          "means"
        })
      }
      df_log <- if (var_equal) {
        df_pooled
      } else {
        w <- c(1 / num, -1 / den)
        if (k == 3L) w <- c(w, 1 / den - 1 / num)
        .satterthwaite_df(w, vm, nu)
      }
      lim <- .log_ratio_limits(num, den, V, .ratio_crit(alpha_lower, df_log),
                               .ratio_crit(alpha_upper, df_log))
      if (!var_equal) {
        # A tail area of 0 leaves its end open (0 or Inf on the ratio
        # scale), and no test decides it.
        attrs$df_limits <- c(
          lower = if (alpha_lower > 0) df_log else NA_real_,
          upper = if (alpha_upper > 0) df_log else NA_real_,
          excluded_lower = NA_real_, excluded_upper = NA_real_)
      }
      list(lower = lim[["lower"]], upper = lim[["upper"]],
           excluded_lower = NA_real_, excluded_upper = NA_real_,
           shape = "interval")
    },
    {
      # The number of replications follows the rule of ?dmar_bootstrap,
      # at the smaller of the tail areas asked for.
      rule <- .ratio_boot_rule(alpha_lower, alpha_upper, conf_level,
                               ci_method)
      B <- .ratio_check_B(B, rule)
      .dmar_local_seed(seed)
      boot_x <- .ratio_boot_means(groups$x, B)
      boot_y <- .ratio_boot_means(groups$y, B)
      boot_b <- if (k == 3L) .ratio_boot_means(groups$baseline, B) else 0
      boot_den <- boot_y - boot_b
      kept <- .ratio_boot_keep((boot_x - boot_b) / boot_den, boot_den, B,
                               rule, "ci_ratio_of_means")
      settled <- .ratio_boot_limits(
        kept, ratio, alpha_lower, alpha_upper, ci_method,
        acc = function() .ratio_acceleration_groups(groups$x, groups$y,
                                                    groups$baseline))
      lim <- settled$limits
      attrs$B <- settled$verdict$B
      attrs$B_used <- settled$verdict$B_used
      attrs$n_dropped <- settled$verdict$n_dropped
      attrs$B_required <- settled$verdict$B_required
      attrs$bootstrap_status <- settled$verdict$status
      list(lower = lim[1L], upper = lim[2L], excluded_lower = NA_real_,
           excluded_upper = NA_real_, shape = "interval")
    })
  if (ci_method != "fieller") {
    attrs$fieller_set <- .ratio_fieller_note(fieller, alpha_lower,
                                             alpha_upper)
  }

  rows <- c(mean_1 = m[[1L]], mean_2 = m[[2L]],
            mean_baseline = if (k == 3L) m[[3L]] else NA_real_,
            se_1 = sqrt(vm[[1L]]), se_2 = sqrt(vm[[2L]]),
            se_baseline = if (k == 3L) sqrt(vm[[3L]]) else NA_real_,
            t_denominator = den / sqrt(V[2L, 2L]), df = df_den,
            n_1 = n[[1L]], n_2 = n[[2L]],
            n_baseline = if (k == 3L) n[[3L]] else NA_real_)
  .ratio_table(set, ratio, rows, conf_level, alpha_lower, alpha_upper,
               ci_method, attrs)
}
