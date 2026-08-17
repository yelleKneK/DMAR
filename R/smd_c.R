#' Standardized Mean Difference Using the Control Group as the Basis of Standardization
#'
#' @description
#' Estimates the standardized mean difference using the control group standard
#' deviation as the basis of standardization (Glass's \emph{g}), from either
#' raw data or summary statistics, in ordinary or unbiased form. Standardizing
#' by the control group alone keeps the scale of the effect free of any
#' treatment effect on variability.
#'
#' @param group_T Raw data for the treatment group
#' @param group_C Raw data for the control group
#' @param mean_T The mean of the treatment group
#' @param mean_C The mean of the control group
#' @param s_C The standard deviation of the control group (i.e., the square root of the unbiased estimator of the population variance)
#' @param n_C The sample size of the control group
#' @param unbiased Returns the unbiased estimate of the standardized mean difference using the standard deviation of the control group
#'
#' @details
#' When \code{unbiased=TRUE}, the unbiased estimate of the standardized mean difference (using the control
#' group as the basis of standardization) is returned (Hedges, 1981). Although the unbiased estimate of the
#' standardized mean difference is not often reported, at least at the present time, it is nevertheless made
#' available to those who are interested in calculating this quantity.
#'
#' @return
#' A 1-row \code{data.frame} with columns \code{term} (\code{"smd_c"}) and
#' \code{value} (the estimated standardized mean difference using the
#' control group standard deviation as the basis of standardization).
#'
#' @references
#' Glass, G. V. (1976). Primary, secondary, and meta-analysis of research. \emph{Educational Researcher, 5}, 3--8.
#'
#' Hedges, L. V. (1981). Distribution theory for Glass's Estimator of effect size and related estimators.
#' \emph{Journal of Educational Statistics, 6}(2), 107--128.
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
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso
#' \code{\link{smd}}, \code{\link{conf_limits_nct}}
#'
#' @examples
#' # Generate sample data.
#' set.seed(113)
#' g.T <- rnorm(n = 25, mean = .5, sd = 1)
#' g.C <- rnorm(n = 25, mean = 0, sd = 1)
#' smd_c(group_T = g.T, group_C = g.C)
#'
#' M.T <- .66745
#' M.C <- .24878
#' sd.c <- 1.1311
#' n.c <- 25
#' smd_c(mean_T = M.T, mean_C = M.C, s_C = sd.c)
#' smd_c(mean_T = M.T, mean_C = M.C, s_C = sd.c, n_C = n.c, unbiased = TRUE)
#'
#' @keywords models htest
#'
#' @export


smd_c <- function(group_T = NULL, group_C = NULL, mean_T = NULL, mean_C = NULL, s_C = NULL, n_C = NULL, unbiased = FALSE) {
  # Internal helper: Hedges (1981) unbiased correction factor given df. Evaluated
  # on the log scale so it does not overflow: gamma() returns Inf for arguments
  # above about 171.6, which makes the raw ratio NaN for a control-group size
  # n_C beyond roughly 345, whereas lgamma is finite for all df and agrees with
  # the direct ratio to machine precision at small df.
  .hedges_correction <- function(df) {
    exp(lgamma(df / 2) - lgamma((df - 1) / 2)) / sqrt(df / 2)
  }

  # Function to calculate the standardized mean difference.
  if (length(group_T) >= 1 && length(group_C) >= 1 && (!is.null(mean_T) || !is.null(mean_C) || !is.null(s_C))) stop("Since you've specified raw data, you do not need to specify summary values.", call. = FALSE)
  if (length(group_T) == 1 || length(group_C) == 1) stop("You only have 1 individual in one or both groups, you cannot calculate a variance in such a situation.", call. = FALSE)
  if (length(mean_T) == 1 && length(mean_C) == 1 && !is.null(s_C) && is.null(n_C) && unbiased == TRUE) stop("You need to specify the control group sample size in order to obtain the unbiased estimate of the standardized mean difference.", call. = FALSE)

  # Here determine 'smd_c' on summary data.
  if (length(mean_T) == 1 && length(mean_C) == 1) {
    if (!is.null(group_T) || !is.null(group_C)) stop("group_C and group_T should be NULL, since you've specified the group means directly. Alternatively, you can specify the groups directly and make mean_C and mean_T NULL.", call. = FALSE)
    d <- (mean_T - mean_C) / s_C
    if (unbiased == TRUE) d <- d * .hedges_correction(n_C - 1)
    return(.as_dmar_tbl(data.frame(term = 'smd_c', value = d)))
  }

  # Here determine 'smd' for raw data.
  if (length(group_T) > 1 && length(group_C) > 1) {
    if (!is.null(s_C)) stop("Since you've specified raw group data, you should not specify any standard deviations.", call. = FALSE)
    x_bar_T <- mean(group_T)
    x_bar_C <- mean(group_C)
    n_C <- length(group_C)
    s_Con <- sqrt(var(group_C))
    d <- (x_bar_T - x_bar_C) / s_Con
    if (unbiased == TRUE) d <- d * .hedges_correction(n_C - 1)
    return(.as_dmar_tbl(data.frame(term = 'smd_c', value = d)))
  }

  stop("Specify either both raw data vectors ('group_T' and 'group_C'), or both group means ('mean_T' and 'mean_C') together with the control group standard deviation ('s_C').", call. = FALSE)
}
