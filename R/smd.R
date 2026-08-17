#' Standardized Mean Difference
#'
#' Estimates the standardized mean difference (Cohen's \emph{d}), the
#' difference between two group means divided by the pooled standard
#' deviation, from either raw data or summary statistics. Expressing the
#' difference in standard deviation units frees the comparison from the raw
#' measurement units, so effects can be compared across measures and studies;
#' either the ordinary or the unbiased (Hedges, 1981) estimate can be
#' returned.
#'
#' @param group_1 Raw data for group 1
#' @param group_2 Raw data for group 2
#' @param mean_1 The mean of group 1
#' @param mean_2 The mean of group 2
#' @param s_1 The standard deviation of group 1 (i.e., the square root of the unbiased estimator of the population variance)
#' @param s_2 The standard deviation of group 2 (i.e., the square root of the unbiased estimator of the population variance)
#' @param s The pooled group standard deviation (i.e., the square root of the unbiased estimator of the population variance)
#' @param n_1 The sample size within group 1
#' @param n_2 The sample size within group 2
#' @param unbiased Returns the unbiased estimate of the standardized mean difference
#'
#' @details
#' When \code{unbiased=TRUE}, the unbiased estimate of the standardized mean difference is returned (Hedges, 1981).
#'
#' @return
#' A 1-row \code{data.frame} with columns \code{term} (\code{"smd"}) and
#' \code{value} (the estimated standardized mean difference).
#'
#' @references
#' Cohen, J. (1988). \emph{Statistical power analysis for the behavioral
#'   sciences} (2nd ed.). Hillsdale, NJ: Lawrence Erlbaum.
#'
#' Cumming, G., & Finch, S. (2001). A primer on the understanding, use, and
#'   calculation of confidence intervals that are based on central and
#'   noncentral distributions. \emph{Educational and Psychological
#'   Measurement, 61}(4), 532--574. \doi{10.1177/0013164401614002}
#'
#' Hedges, L. V. (1981). Distribution theory for Glass's Estimator of effect size and related estimators.
#' \emph{Journal of Educational Statistics, 6}(2), 107--128.
#'
#' Kelley, K. (2005) The effects of nonnormal distributions on confidence intervals around the standardized
#' mean difference: Bootstrap and parametric confidence intervals, \emph{Educational and Psychological
#' Measurement, 65}, 51--69. \doi{10.1177/0013164404264850}
#'
#' Kelley, K. (2007). Confidence intervals for standardized effect sizes:
#' Theory, application, and implementation. \emph{Journal of Statistical Software, 20}(8), 1--24.
#'   \doi{10.18637/jss.v020.i08}
#'
#' Kelley, K., & Rausch, J. R. (2006). Sample size planning for the
#'   standardized mean difference: Accuracy in parameter estimation via
#'   narrow confidence intervals. \emph{Psychological Methods, 11}(4),
#'   363--385. \doi{10.1037/1082-989X.11.4.363}
#'
#' Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). \emph{Designing
#'   experiments and analyzing data: A model comparison perspective}
#'   (4th ed.). Routledge. (See Chapter 4 on individual comparisons and
#'   Chapter 3 on one-way ANOVA.)
#'
#' Steiger, J. H., & Fouladi, R. T. (1997). Noncentrality interval
#'   estimation and the evaluation of statistical methods. In L. L. Harlow,
#'   S. A. Mulaik, & J. H. Steiger (Eds.), \emph{What if there were no
#'   significance tests?} (pp. 221--257). Mahwah, NJ: Lawrence Erlbaum.
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso
#' \code{\link{smd_c}}, \code{\link{ci_smd}}, \code{\link{ci_smd_c}},
#' \code{\link{ss_aipe_smd}}, \code{\link{ss_power_smd}},
#' \code{\link{plot_smd}}, \code{\link{conf_limits_nct}}
#'
#' @examples
#' # Generate sample data.
#' set.seed(113)
#' g.1 <- rnorm(n = 25, mean = .5, sd = 1)
#' g.2 <- rnorm(n = 25, mean = 0, sd = 1)
#' smd(group_1 = g.1, group_2 = g.2)
#'
#' M.x <- .66745
#' M.y <- .24878
#' sd <- 1.048
#' smd(mean_1 = M.x, mean_2 = M.y, s = sd)
#'
#' M.x <- .66745
#' M.y <- .24878
#' n1 <- 25
#' n2 <- 25
#' sd.1 <- .95817
#' sd.2 <- 1.1311
#' smd(mean_1 = M.x, mean_2 = M.y, s_1 = sd.1, s_2 = sd.2, n_1 = n1, n_2 = n2)
#'
#' smd(mean_1 = M.x, mean_2 = M.y, s_1 = sd.1, s_2 = sd.2, n_1 = n1, n_2 = n2,
#'     unbiased = TRUE)
#'
#' @concept Cohen's d
#' @concept Hedges' g
#'
#' @keywords models htest
#'
#' @export


smd <- function(group_1 = NULL, group_2 = NULL, mean_1 = NULL, mean_2 = NULL, s_1 = NULL, s_2 = NULL, s = NULL, n_1 = NULL, n_2 = NULL, unbiased = FALSE) {
  # Internal helper: Hedges (1981) unbiased correction factor given total df.
  # Evaluated on the log scale so it does not overflow: gamma() returns Inf for
  # arguments above about 171.6, which would make the raw ratio NaN for total df
  # of roughly 344 or more (ordinary sample sizes), whereas lgamma is finite for
  # all df and agrees with the direct ratio to machine precision at small df.
  .hedges_correction <- function(df) {
    exp(lgamma(df / 2) - lgamma((df - 1) / 2)) / sqrt(df / 2)
  }

  # Function to calculate the standardized mean difference.
  if (length(group_1) >= 1 && length(group_2) >= 1 && (!is.null(mean_1) || !is.null(mean_2) || !is.null(s_1) || !is.null(s_2) || !is.null(s))) stop("Since you've specified raw data, you do not need to specify summary values.", call. = FALSE)
  if (length(group_1) == 1 || length(group_2) == 1) stop("You only have 1 individual in one or both groups, you cannot calculate a variance in such a situation.", call. = FALSE)
  if (length(mean_1) == 1 && length(mean_2) == 1 && !is.null(s) && (is.null(n_1) || is.null(n_2)) && unbiased == TRUE) stop("You need to specify the group sample sizes in order to obtain the unbiased estimate of the standardized mean difference.", call. = FALSE)

  # Here determine 'smd' on summary data.
  if (length(mean_1) == 1 && length(mean_2) == 1) {
    if (!is.null(group_1) || !is.null(group_2)) stop("group_1 and group_2 should be NULL, since you've specified the group means directly. Alternatively, you can specify the groups directly and make mean_1 and mean_2 NULL.", call. = FALSE)
    if (is.null(s_1) && is.null(s_2) && length(s) == 1) {
      d <- (mean_1 - mean_2) / s
      if (unbiased == TRUE) d <- d * .hedges_correction(n_1 + n_2 - 2)
      return(.as_dmar_tbl(data.frame(term = 'smd', value = d)))
    }

    if (length(s_1) == 1 && length(s_2) == 1 && is.null(s)) {
      if (is.null(n_1) || is.null(n_2) || length(n_1) == 0 || length(n_2) == 0) stop("You did not specify the per group sample sizes (i.e., 'n_1' and/or 'n_2')", call. = FALSE)
      # Calculate the pooled variance.
      sd <- sqrt((s_1^2 * (n_1 - 1) + s_2^2 * (n_2 - 1)) / (n_1 + n_2 - 2))
      d <- (mean_1 - mean_2) / sd
      if (unbiased == TRUE) d <- d * .hedges_correction(n_1 + n_2 - 2)
      return(.as_dmar_tbl(data.frame(term = 'smd', value = d)))
    }
  }

  # Here determine 'smd' for raw data.
  if (length(group_1) > 1 && length(group_2) > 1) {
    if (!is.null(s_1) || !is.null(s_2) || !is.null(s)) stop("Since you've specified raw group data, you should not specify any standard deviations.", call. = FALSE)
    x_bar_1 <- mean(group_1)
    x_bar_2 <- mean(group_2)
    SS1 <- sum((group_1 - x_bar_1)^2)
    SS2 <- sum((group_2 - x_bar_2)^2)
    n_1 <- length(group_1)
    n_2 <- length(group_2)
    s <- sqrt((SS1 + SS2) / (n_1 + n_2 - 2))
    d <- (x_bar_1 - x_bar_2) / s
    if (unbiased == TRUE) d <- d * .hedges_correction(n_1 + n_2 - 2)
    return(.as_dmar_tbl(data.frame(term = 'smd', value = d)))
  }

  stop("Specify either both raw data vectors ('group_1' and 'group_2'), or both group means ('mean_1' and 'mean_2') together with a standard deviation (either 's', or 's_1' and 's_2' with 'n_1' and 'n_2').", call. = FALSE)
}
