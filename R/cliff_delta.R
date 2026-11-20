# Cliff's delta (ordinal dominance effect size).
#' Cliff's \eqn{\delta} Ordinal Effect Size
#'
#' Takes two independent samples of scores that can be ordered and returns
#' Cliff's \eqn{\delta}, the probability that a randomly drawn observation
#' from the first population exceeds one from the second minus the
#' probability of the reverse, with a confidence interval and the
#' proportions of pairs won, lost, and tied. The statistic uses only the
#' order of the observations, so it needs no assumption of normality or of
#' equal variances and no transformation of a skewed outcome, and its
#' interval comes from a variance formula, not from resampling.
#'
#' @param group_1,group_2 Numeric or logical vectors, or ordered factors
#'   with the same levels in the same order, holding the observations in the
#'   two groups. The statistic uses only their order. Missing values are
#'   removed group by group (see the \code{n_missing} attribute). An
#'   unordered factor and a character vector are refused; see \emph{What
#'   is refused, and why}.
#' @param conf_level Confidence level for the CI. Default \code{0.95}.
#'
#' @return A \code{data.frame} (classes \code{dmar_cliff_delta} and
#'   \code{dmar_tbl}) with columns \code{term} and \code{value} and these
#'   rows:
#'   \describe{
#'   \item{\code{cliff_delta}}{The estimate \eqn{\hat\delta}, in
#'     \eqn{[-1, 1]}.}
#'   \item{\code{ci_lower}, \code{ci_upper}}{The confidence limits
#'     at \code{conf_level}. Both equal the estimate when the estimated
#'     variance is 0, and are then not a confidence interval.}
#'   \item{\code{var_cliff_delta}}{The estimated variance of
#'     \eqn{\hat\delta} (see \emph{Variance}).}
#'   \item{\code{p_y1_greater}}{The proportion of the \eqn{n_1 n_2} pairs,
#'     one observation from each group, in which the observation from
#'     group 1 is the larger.}
#'   \item{\code{p_y1_less}}{The proportion in which it is the smaller.}
#'   \item{\code{p_tied}}{The proportion of tied pairs. The three
#'     proportions sum to 1, and \code{cliff_delta} is the first minus the
#'     second.}
#'   }
#'   The table carries the attributes \code{conf_level}; \code{n_missing},
#'   the numbers of missing values removed from the two groups, as
#'   \code{c(group_1 = , group_2 = )}; and \code{interval_degenerate},
#'   present only when the estimated variance is 0, with the value
#'   \code{"complete_separation"} or \code{"no_variation"}. The
#'   \code{tidy()} and \code{glance()} verbs (see
#'   \code{\link{dmar_tidiers}}) add a logical column
#'   \code{interval_degenerate}.
#'
#' @details
#' \strong{Estimand and interval.} For independent observations
#' \eqn{Y_1} from the first population and \eqn{Y_2} from the second,
#' \deqn{\delta \;=\; \Pr(Y_1 > Y_2) - \Pr(Y_1 < Y_2)
#'   \;=\; 2 \cdot A - 1,}
#' where \eqn{A} is the Vargha-Delaney (2000) statistic. The sampling
#' units are the \eqn{n_1 + n_2} individuals, each observed once and each
#' in one group. The interval is a large-sample confidence interval for
#' the population \eqn{\delta}, formed on the inverse hyperbolic tangent
#' scale. \eqn{\delta} ranges over \eqn{[-1, 1]}; 0 says that neither
#' group tends to score higher, and 1 or \eqn{-1} that the two
#' distributions do not overlap. Feng and Cliff (2004, p. 322) credit the
#' statistic to Cliff (1993).
#'
#' \strong{Sample estimate.}
#' \deqn{\hat\delta = \frac{\#\{(i,j): y_{1i} > y_{2j}\}
#'   - \#\{(i,j): y_{1i} < y_{2j}\}}{n_1 n_2},}
#' the mean of the \eqn{n_1 n_2} signs
#' \eqn{d_{ij} = \mathrm{sign}(y_{1i} - y_{2j})}. A tied pair has sign 0
#' and counts for neither side. The rows \code{p_y1_greater},
#' \code{p_y1_less}, and \code{p_tied} give the three proportions of
#' pairs.
#'
#' \strong{Assumptions.} The two samples are independent of each other and
#' the observations within each are independent and identically
#' distributed; for paired observations use
#' \code{\link{probability_of_superiority_paired}}. The response is at
#' least ordinal. The estimate needs nothing more: no distributional
#' shape, no equality of variances, and ties are allowed. The interval is
#' a large-sample result; \emph{When the estimated variance is zero} says
#' how it behaves in small samples with a large effect.
#'
#' \strong{Variance.} The variance of \eqn{\hat\delta} is
#' \deqn{\mathrm{Var}(\hat\delta) \;=\;
#'   \frac{(n_2 - 1) \sigma^2_{d_1} + (n_1 - 1) \sigma^2_{d_2}
#'         + \sigma^2_d}{n_1 n_2}}
#' (Feng & Cliff, 2004, Equation 2, p. 323, who credit Cliff, 1993), where
#' \eqn{\sigma^2_{d_1}} is the variance of the row means
#' \eqn{d_{i\cdot}}, the mean sign of one observation of the first group
#' against the whole second population, \eqn{\sigma^2_{d_2}} is the
#' variance of the column means \eqn{d_{\cdot j}}, and \eqn{\sigma^2_d} is
#' the variance of the signs \eqn{d_{ij}}. (The equation as printed there
#' attaches \eqn{(n_1 - 1)} to the variance of the first group's row means
#' and \eqn{(n_2 - 1)} to the second's. The coefficients are the other way
#' around, as the variance of a two-sample U statistic shows: as
#' \eqn{n_2} grows, \eqn{\hat\delta} becomes the mean of \eqn{n_1}
#' independent row means, whose variance is \eqn{\sigma^2_{d_1}/n_1}. The
#' difference matters only for unequal group sizes.) The function
#' estimates the three components by their sample values: the sample
#' variances (divisors \eqn{n_1 - 1} and \eqn{n_2 - 1}) of the observed
#' row means and column means, and the sample variance (divisor
#' \eqn{n_1 n_2 - 1}) of the \eqn{n_1 n_2} signs. Feng and Cliff's
#' Equation 4, an unbiased estimate with a minimum allowed value, is not
#' used. The result is the row \code{var_cliff_delta}.
#'
#' \strong{Confidence interval.} The interval is formed on the scale of
#' \eqn{\mathrm{arctanh}(\delta)}, as Fisher's \emph{Z} interval for a
#' correlation is, and transformed back, so that the limits stay inside
#' \eqn{[-1, 1]}:
#' \deqn{\tanh\!\left\{\mathrm{arctanh}(\hat\delta) \pm z_{1 - \alpha/2}\,
#'   \frac{\sqrt{\widehat{\mathrm{Var}}(\hat\delta)}}{1 - \hat\delta^2}
#'   \right\},}
#' the standard error on the transformed scale coming from the delta
#' method. Feng and Cliff (2004, Equation 5, p. 324) give an asymmetric
#' interval of another form, which this function does not compute.
#'
#' \strong{Relation to Vargha and Delaney's \emph{A}.}
#' \eqn{A = (\delta + 1)/2}, and because
#' \eqn{\mathrm{arctanh}(2A - 1) = \mathrm{logit}(A)/2}, the interval here
#' and the logit interval of \code{\link{vargha_delaney_A}} are the same
#' construction on the two scales. The two functions estimate the variance
#' differently, however (the three-term estimate above here; the two
#' placement variances of DeLong, DeLong, and Clarke-Pearson, 1988,
#' there), so their limits are not the same after the change of scale. For
#' \code{group_1 = 1:5} and \code{group_2 = 3:7}, \eqn{\hat\delta = -0.64}
#' with limits \eqn{-0.935} and 0.177; \code{vargha_delaney_A()} gives
#' \eqn{\hat A = 0.18} with limits 0.034 and 0.576, which are
#' \eqn{-0.931} and 0.151 on the scale of \eqn{\delta}.
#'
#' \strong{Scale and missing values.} The observations enter only through
#' their order, so any increasing rescoring of the response gives the same
#' result, and an ordinal response may be supplied as an ordered factor.
#' Missing values (\code{NA} and \code{NaN}) are removed group by group;
#' the numbers removed are in the \code{n_missing} attribute. Infinite
#' values are ordered as usual. Dates and date-times (classes \code{Date}
#' and \code{POSIXct}) are ordered by their numeric values, as are two
#' time differences (class \code{difftime}) in the same units. A
#' one-column matrix and a one-dimensional array are read as the vectors
#' they are. Each group needs at least two non-missing observations; with
#' fewer the function stops, since no variance can be estimated.
#'
#' \strong{What is refused, and why.} The statistic uses the order of the
#' observations and nothing else, so the function must be given objects
#' whose order is the order meant. It refuses, with a condition of class
#' \code{dmar_scores_not_ordered} (see \code{\link{dmar_metadata}}), an
#' unordered factor (its codes follow the order of \code{levels()},
#' alphabetical unless set, so the levels low, mid, high are stored as 2,
#' 3, 1, and \eqn{\hat\delta} computed from the codes of the example in
#' \emph{Examples} is \eqn{-0.19} where the ordered levels give 0.31); a
#' character vector (strings sort alphabetically, \code{"10"} before
#' \code{"9"}); an ordered factor beside anything that is not an ordered
#' factor with the same levels in the same order (each factor numbers its
#' own levels from 1); time differences in different units, or one beside
#' plain numbers (each such object is stored in its own units); a
#' \code{Date} beside a \code{POSIXct}; and a list, data frame, or matrix
#' of several columns. Nothing is computed in any of these cases. The
#' message says which case applies and what to supply: an ordered factor
#' with its levels from lowest to highest,
#' \code{factor(x, levels = c("low", "mid", "high"), ordered = TRUE)}, or
#' numeric scores.
#'
#' \strong{When the estimated variance is zero.} Two kinds of data give a
#' variance estimate of exactly 0, and they are the only two.
#' \emph{Complete separation}: every observation in one group exceeds
#' every observation in the other, so the estimate is at the end of its
#' range (\eqn{\hat\delta} \eqn{-1} or 1). \emph{No variation}: every
#' observation has the same value, so the estimate is 0. In both, the two
#' limits are reported equal to the estimate, \code{var_cliff_delta} is 0,
#' the table's \code{interval_degenerate} attribute names the case, and a
#' warning of class \code{dmar_interval_degenerate} is given whose
#' \code{reason} element is \code{"complete_separation"} or
#' \code{"no_variation"}; the condition also carries \code{estimate},
#' \code{n} (the two group sizes), and \code{limit_feng_cliff}. The
#' reported limits are not a confidence interval with the stated coverage.
#' The zero describes the sample: the variance estimate measures how much
#' the comparisons differ from one observation to the next, and without
#' overlap they do not differ. Complete separation is common when the
#' effect is large and the samples are small. With five observations per
#' group from normal populations, it occurred in 19 percent of samples
#' when the population value of \eqn{\delta} was 0.71 and in 59 percent
#' when it was 0.92, and across all samples the reported 95 percent limits
#' then covered the population value 79 and about 40 percent of the time,
#' against 98 and 96 percent among the samples that were not separated
#' (20,000 simulated samples each). Feng and Cliff (2004, p. 324,
#' Equations 6 to 8) treat this case by assuming the largest variance
#' \eqn{\hat\delta} can have, which gives the limit
#' \eqn{(n_b - z^2)/(n_b + z^2)} when \eqn{\hat\delta = 1} (and its
#' negative when \eqn{\hat\delta = -1}), with \eqn{n_b} the smaller group
#' size; they describe the result as conservative. At the 95 percent level
#' it is \eqn{-0.123}, 0.020, 0.131, 0.445, and 0.678 for \eqn{n_b} of 3,
#' 4, 5, 10, and 20. The warning states that limit. It is not the limit
#' this function reports.
#'
#' Particular cases. Infinite values are ordered like any others, so a
#' group made up of \code{Inf} beside finite scores is completely
#' separated, and two groups made up of the same infinite value have no
#' variation. Two constant groups at different values are completely
#' separated. Groups that only touch, the largest score of one equal to
#' the smallest of the other, are not separated: the tied pair gives a
#' positive variance and an ordinary interval, which is wide
#' (\code{c(3, 4, 5)} against \code{c(1, 2, 3)} gives
#' \eqn{\hat\delta = 0.889} with limits \eqn{-0.167} and 0.995). The mark
#' travels with the result: \code{\link{results_sentence}} writes the
#' estimate and says that there is no confidence interval, the printed
#' table says so beneath itself, \code{tidy()} and \code{glance()} carry
#' it in the column \code{interval_degenerate},
#' \code{\link{randomization_test}} passes this warning on, and
#' \code{\link{ss_aipe_cliff_delta_sensitivity}} counts a separated
#' replication as one without an interval. The randomization test remains
#' valid in both cases.
#'
#' \strong{Relation to other measures.} Cliff's \eqn{\delta} is a linear
#' function of the Vargha-Delaney (2000) \eqn{A} statistic
#' (\eqn{\delta = 2A - 1}) and of the Mann-Whitney \eqn{U} statistic
#' (\eqn{U / (n_1 n_2) = A}). The common language effect size computed by
#' \code{\link{cles}} is \eqn{\Pr(Y_1 > Y_2)} for normal distributions
#' with a common variance; \eqn{(\delta + 1)/2} estimates that probability
#' (with ties split) without those assumptions, which suits skewed
#' outcomes and ordinal scales.
#'
#' \strong{Reading the output.} In the third example below the 32
#' participants given the Community Reinforcement Approach are group 1 and
#' the 37 given standard care are group 2, so there are 1,184 pairs. The
#' CRA participant reports more drinks per week in 29.3 percent of the
#' pairs (\code{p_y1_greater}), fewer in 59.3 percent (\code{p_y1_less}),
#' and the same number in 11.4 percent (\code{p_tied}), which gives
#' \eqn{\hat\delta = 0.293 - 0.593 = -0.300}, with 95 percent limits
#' \eqn{-0.535} and \eqn{-0.021}.
#'
#' @references
#' Cliff, N. (1993). Dominance statistics: Ordinal analyses to answer
#'   ordinal questions. \emph{Psychological Bulletin, 114}(3), 494--509.
#'   \doi{10.1037/0033-2909.114.3.494}
#'
#' DeLong, E. R., DeLong, D. M., & Clarke-Pearson, D. L. (1988). Comparing
#'   the areas under two or more correlated receiver operating
#'   characteristic curves: A nonparametric approach. \emph{Biometrics,
#'   44}(3), 837--845. \doi{10.2307/2531595}
#'
#' Feng, D., & Cliff, N. (2004). Monte Carlo evaluation of ordinal d with
#'   improved confidence interval. \emph{Journal of Modern Applied
#'   Statistical Methods, 3}(2), 322--332. \doi{10.22237/jmasm/1099267560}
#'
#' Smith, J. E., Meyers, R. J., & Delaney, H. D. (1998). The community
#'   reinforcement approach with homeless alcohol-dependent individuals.
#'   \emph{Journal of Consulting and Clinical Psychology, 66}(3), 541--548.
#'
#' Vargha, A., & Delaney, H. D. (2000). A critique and improvement of
#'   the CL common language effect size statistics of McGraw and Wong.
#'   \emph{Journal of Educational and Behavioral Statistics, 25}(2),
#'   101--132. \doi{10.3102/10769986025002101}
#'
#' @seealso \code{\link{vargha_delaney_A}},
#'   \code{\link{probability_of_superiority_paired}}, \code{\link{cles}},
#'   \code{\link{randomization_test}}, \code{\link{nnt_from_smd}},
#'   \code{\link{ss_aipe_cliff_delta}}
#'
#' @examples
#' # 1. Two groups of different sizes, no ties:
#' set.seed(113)
#' a <- rnorm(30, mean = 0, sd = 1)
#' b <- rnorm(40, mean = 0.5, sd = 1)
#' cliff_delta(a, b)
#'
#' # 2. With ties (ordinal data):
#' o1 <- c(1, 2, 2, 3, 3, 3, 4, 4, 5)
#' o2 <- c(2, 3, 3, 4, 4, 5, 5, 5)
#' cliff_delta(o1, o2)
#'
#' # 3. Robust to right skew. Cliff's delta on the raw, untransformed
#' #    drinking outcome from the Smith, Meyers, and Delaney (1998)
#' #    trial, comparing the Community Reinforcement Approach (CRA)
#' #    against standard care. Because the statistic uses only ranks it
#' #    needs no normalizing transformation of the heavily skewed
#' #    outcome, unlike the standardized mean difference.
#' data(drinks_trial)
#' cra <- drinks_trial$drinks_per_week[drinks_trial$treatment == "CRA"]
#' std <- drinks_trial$drinks_per_week[drinks_trial$treatment == "Standard"]
#' cliff_delta(cra, std)
#'
#' # 4. An ordinal response is supplied as an ordered factor with its
#' #    levels stated from lowest to highest. The unordered factors
#' #    rating_1 and rating_2 themselves are refused, with a message that
#' #    says why: their codes follow the alphabetical order of the labels
#' #    (high, low, mid).
#' rating_1 <- factor(c("low", "mid", "high", "high"))
#' rating_2 <- factor(c("low", "low", "mid", "high"))
#' lev <- c("low", "mid", "high")
#' cliff_delta(factor(rating_1, levels = lev, ordered = TRUE),
#'             factor(rating_2, levels = lev, ordered = TRUE))
#'
#' # 5. Completely separated groups: the estimated variance is 0, the two
#' #    limits equal the estimate, and a warning says that they are not a
#' #    confidence interval and gives the limit of Feng and Cliff (2004).
#' separated <- cliff_delta(6:10, 1:5)
#' separated
#' results_sentence(separated)
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @keywords htest
#'
#' @family effect size estimates
#'
#' @export

cliff_delta <- function(group_1, group_2, conf_level = 0.95) {
  .dominance_check_conf_level(conf_level)
  scores  <- .dominance_scores(group_1, group_2, "group_1", "group_2",
                               "cliff_delta")
  group_1 <- scores$x
  group_2 <- scores$y
  if (length(group_1) < 2L || length(group_2) < 2L)
    .dmar_refuse_argument(
      "Each group must have at least 2 non-missing observations.", arg = NULL)

  n_1 <- length(group_1)
  n_2 <- length(group_2)
  M   <- as.numeric(n_1) * n_2   # a double: n_1 * n_2 can pass the largest integer

  # The pairwise comparisons come from ranks (.dominance_placements()), not
  # from the n_1 by n_2 matrix of signs: +1 if y_1 > y_2, -1 if y_1 < y_2,
  # 0 if tied.
  placements <- .dominance_placements(group_1, group_2)
  p_gt <- placements$p_greater
  p_lt <- placements$p_less
  p_eq <- placements$p_tied
  delta_hat <- p_gt - p_lt

  # The plug-in estimate of the variance in Feng and Cliff (2004,
  # Equation 2), coefficients corrected: (n_2 - 1) goes with the variance
  # of the first group's row means and (n_1 - 1) with the second's.
  d_i <- 2 * placements$V_10 - 1      # per-group-1 dominance scores
  d_j <- 2 * placements$V_01 - 1      # per-group-2 dominance scores
  sigma_di <- var(d_i)
  sigma_dj <- var(d_j)
  # The variance of the n_1 n_2 signs: their mean square is p_gt + p_lt.
  sigma_d  <- (M / (M - 1)) * (p_gt + p_lt - delta_hat^2)

  var_delta <- ((n_2 - 1) * sigma_di + (n_1 - 1) * sigma_dj + sigma_d) / M

  # CI on Fisher-arctanh scale, back-transformed to handle [-1, 1]:
  z_alpha <- stats::qnorm(1 - (1 - conf_level) / 2)
  degenerate <- NULL
  if (abs(delta_hat) == 1) {
    # Complete separation: the plug-in variance is 0 and arctanh(+/-1) diverges,
    # so se_z = sqrt(var_delta) / (1 - delta_hat^2) is 0 / 0 = NaN. Report the
    # degenerate boundary interval and flag it, rather than returning NaN limits.
    degenerate <- .dominance_degenerate_interval(
      "cliff_delta", "complete_separation",
      "the groups are completely separated (|delta| = 1)",
      estimate = delta_hat, n = c(n_1, n_2), conf_level = conf_level)
    lo <- hi <- delta_hat
  } else {
    if (p_eq == 1) {
      # Every observation has the same value: the variance is 0 here too,
      # and the limits below are both 0.
      degenerate <- .dominance_degenerate_interval(
        "cliff_delta", "no_variation",
        "all observations have the same value",
        estimate = delta_hat, n = c(n_1, n_2), conf_level = conf_level)
    }
    z_delta <- atanh(delta_hat)
    se_z    <- sqrt(var_delta) / (1 - delta_hat^2)
    lo <- tanh(z_delta - z_alpha * se_z)
    hi <- tanh(z_delta + z_alpha * se_z)
  }

  out <- data.frame(
    term  = c("cliff_delta", "ci_lower", "ci_upper",
              "var_cliff_delta", "p_y1_greater", "p_y1_less", "p_tied"),
    value = c(delta_hat, lo, hi, var_delta, p_gt, p_lt, p_eq),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  out <- .as_dmar_tbl(out, conf_level = conf_level,
                      subclass = "dmar_cliff_delta")
  attr(out, "n_missing") <- scores$n_missing
  # Present only when the estimated variance is 0.
  attr(out, "interval_degenerate") <- degenerate
  out
}

# A table of cliff_delta() says beneath itself when its limits are not a
# confidence interval.
#' @exportS3Method .dmar_tbl_footer dmar_cliff_delta
#' @noRd
.dmar_tbl_footer.dmar_cliff_delta <- function(x, ...) {
  .dmar_footer_add(NextMethod(), .dominance_degenerate_footer(x),
                   "interval_degenerate")
}

# The broom verbs of cliff_delta() are the long table's own (tidy.dmar_tbl
# and glance.dmar_tbl) with one more column, the logical
# interval_degenerate, which says in the rows themselves that the limits
# are not a confidence interval (an estimated variance of 0), so that the
# mark stays with the rows in a stack.
#' @rdname dmar_tidiers
#' @importFrom generics tidy
#' @exportS3Method generics::tidy
tidy.dmar_cliff_delta <- function(x, ...) {
  .dominance_mark_degenerate(NextMethod(), x)
}

#' @rdname dmar_tidiers
#' @importFrom generics glance
#' @exportS3Method generics::glance
glance.dmar_cliff_delta <- function(x, ...) {
  .dominance_mark_degenerate(NextMethod(), x)
}
