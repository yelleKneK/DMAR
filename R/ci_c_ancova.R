#' Confidence Interval for an (Unstandardized) Contrast in ANCOVA With One Covariate
#'
#' @description
#' Calculates the confidence interval for an unstandardized contrast in the
#' one-covariate ANCOVA. Two procedures are available through the
#' \code{procedure} argument. The default (\code{"t"}) returns a single
#' \emph{per-comparison} interval based on the \emph{t} distribution: it gives
#' the correct \eqn{(1 -} \code{conf_level}\eqn{)} coverage for one contrast
#' chosen in advance, and its standard error includes the
#' \eqn{(\sum c_i \bar X_i)^2 / SS_{\mathrm{within}(x)}} term that accounts for
#' the covariate separation between the groups in that one contrast. The
#' \code{"bryant_paulson"} procedure instead returns Bryant--Paulson
#' \emph{simultaneous} (familywise) intervals over a whole family of contrasts
#' of adjusted means; when selected, \code{ci_c_ancova} simply forwards its
#' arguments to \code{\link{ci_c_ancova_bp}} and returns that result. See
#' Details for which to use when.
#'
#' @details
#' \strong{Per-comparison versus simultaneous.} The two procedures answer
#' different questions and are not interchangeable. Use the default
#' \code{procedure = "t"} when a single contrast was planned in advance: the
#' interval has exact per-comparison coverage and its width reflects the
#' covariate adjustment for that specific contrast through the
#' \eqn{(\sum c_i \bar X_i)^2 / SS_{\mathrm{within}(x)}} term. Use
#' \code{procedure = "bryant_paulson"} when several contrasts (for example all
#' pairwise comparisons of adjusted means) are examined together and the
#' coverage statement must hold simultaneously across the family: the
#' Bryant--Paulson generalized studentized range supplies a larger critical
#' value that controls the familywise error rate and, because the covariates
#' are random, correctly absorbs the extra sampling uncertainty from estimating
#' the covariate adjustment (which holds on average over the covariate
#' distribution, so the per-contrast separation term is not added again). A
#' per-comparison interval used as if it were simultaneous understates the
#' family error rate; a simultaneous interval used for one planned contrast is
#' wider than necessary. See \code{\link{ci_c_ancova_bp}} for the full
#' description of the simultaneous procedure and its arguments.
#'
#' @param contrast The unstandardized contrast of adjusted means
#' @param adj_means The vector that contains the adjusted mean of each group on the dependent variable
#' @param s_ancova The standard deviation of the errors from the ANCOVA model (i.e., the square root of the mean square error from ANCOVA); required, as a single positive number, with either procedure, since the standard error of the contrast is formed from it and the adjusted means and sample sizes cannot supply it
#' @param c_weights The contrast weights
#' @param n_per_group Either a single number that indicates the sample size \emph{per group} or a vector that contains the sample size of each group
#' @param cov_means A vector that contains the group means of the covariate
#' @param SSwithin_x The sum of squares within groups obtained from the summary table for ANOVA on the covariate
#' @param conf_level The desired confidence interval coverage, (i.e., 1 - Type I error rate)
#' @param procedure The interval procedure, one of \code{"t"} (the default, a
#'   single per-comparison interval based on the \emph{t} distribution) or
#'   \code{"bryant_paulson"} (Bryant--Paulson simultaneous intervals for a
#'   family of contrasts of adjusted means). When \code{"bryant_paulson"} is
#'   chosen the call is forwarded to \code{\link{ci_c_ancova_bp}}; see Details.
#' @param \dots Allows one to potentially include parameter values for inner functions.
#'   When \code{procedure = "bryant_paulson"}, these are passed on to
#'   \code{\link{ci_c_ancova_bp}} (for example \code{n_covariates},
#'   \code{df}, or \code{contrast_type}).
#'
#' @return
#' A 3-row \code{data.frame} with columns \code{term} and \code{value}
#' (numeric). The \code{term} values are \code{"ci_lower"} (the lower
#' confidence limit on the unstandardized ANCOVA contrast), \code{"contrast"}
#' (the unstandardized contrast point estimate), and \code{"ci_upper"}
#' (the upper limit).
#'
#' When \code{procedure = "bryant_paulson"}, the return value is whatever
#' \code{\link{ci_c_ancova_bp}} returns (a table with one row per
#' contrast and columns \code{contrast}, \code{estimate}, \code{ci_lower},
#' and \code{ci_upper}).
#'
#' @references
#' Kelley, K. (2007). Confidence intervals for standardized
#'   effect sizes: Theory, application, and implementation.
#'   \emph{Journal of Statistical Software, 20}(8), 1--24.
#'   \doi{10.18637/jss.v020.i08}
#'
#' Lai, K., & Kelley, K. (2012). Accuracy in parameter estimation for
#'   ANCOVA and ANOVA contrasts: Sample size planning via narrow
#'   confidence intervals.
#'   \emph{British Journal of Mathematical and Statistical Psychology, 65},
#'   350--370. \doi{10.1111/j.2044-8317.2011.02029.x}
#'
#' Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). \emph{Designing
#'   experiments and analyzing data: A model comparison perspective}
#'   (4th ed.). Routledge. (See Chapter 9.)
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @note
#' Be sure to use the standard deviation and not the error variance for \code{s_ancova},
#' not the square of this value which would come from the source table
#' (i.e., do not use the variance of the error but rather use the square root).
#'
#' If \code{n_per_group} receives a single number, that number is considered as the sample size \emph{per group}.
#' If \code{n_per_group} receives a vector, the vector is considered as the sample size of each group.
#'
#' Be sure to use fractions not the integers to specify \code{c_weights}. For example, in an ANCOVA of four groups,
#' if the user wants to compare the mean of group 1 and 2 with the mean of group 3 and 4, \code{c_weights} should
#' be specified as c(0.5, 0.5, -0.5, -0.5) rather than c(1, 1, -1, -1). Make sure the sum of the contrast weights are zero.
#'
#' @seealso
#' \code{\link{ci_c}},\code{\link{ci_sc_ancova}},
#' \code{\link{ci_c_ancova_bp}} for the Bryant--Paulson simultaneous procedure
#'
#' @examples
#' # Maxwell, Delaney, & Kelley (2027) offer an example that 30 depressive
#' # individuals are randomly assigned to three groups, 10 in each, and ANCOVA
#' # is performed on the posttest scores using the participants' pretest
#' # scores as the covariate. The means of pretest scores of group 1 to 3 are
#' # 17, 17.7, and 17.4, respectively, and the adjusted means of groups 1 to 3
#' # are 7.5, 12, and 14, respectively. The error variance in ANCOVA is 29,
#' # and the sum of squares within groups from ANOVA on the covariate is 752.5.
#'
#' # To obtain the confidence interval for adjusted mean of group 1 versus group 2:
#' ci_c_ancova(adj_means = c(7.5, 12, 14), s_ancova = sqrt(29),
#'             c_weights = c(1, -1, 0), n_per_group = 10,
#'             cov_means = c(17, 17.7, 17.4), SSwithin_x = 752.5)
#'
#' # That interval is the right one for a single contrast planned in advance.
#' # For the family of all three pairwise comparisons of the adjusted means,
#' # with coverage that holds simultaneously across the family, select
#' # procedure = "bryant_paulson"; the call is forwarded to ci_c_ancova_bp().
#' # The simultaneous limits are wider, which is the price of the family
#' # statement: group 1 against group 2 runs from -10.61 to 1.61, where the
#' # single planned contrast above ran from -9.46 to 0.46.
#' ci_c_ancova(adj_means = c(7.5, 12, 14), s_ancova = sqrt(29), n_per_group = 10,
#'             procedure = "bryant_paulson")
#'
#' @keywords design
#'
#' @family confidence intervals for effect sizes
#'
#' @export


ci_c_ancova <- function(contrast = NULL, adj_means = NULL, s_ancova = NULL, c_weights, n_per_group, cov_means, SSwithin_x, conf_level = .95, procedure = c("t", "bryant_paulson"), ...) {
  .dmar_check_conf_level(conf_level, null_ok = TRUE)
  # The contrast was 'psi' through DMAR 1.0.x. The function takes `...`, so
  # that name would not stop with "unused argument": under procedure = "t"
  # it would be dropped without a word. It is refused here instead.
  # The name is read from the call, so the arguments in `...` are not
  # evaluated here.
  if ("psi" %in% names(match.call(expand.dots = FALSE)$...))
    .dmar_refuse_argument(c(
      "'psi' is not an argument of ci_c_ancova(): give the unstandardized ",
      "contrast of adjusted means as 'contrast'."), arg = "psi")
  procedure <- match.arg(procedure)

  # The Bryant-Paulson simultaneous procedure is implemented in ci_c_ancova_bp();
  # forward the relevant arguments and return its result unchanged. The single
  # per-comparison ('t') interval is computed below.
  if (procedure == "bryant_paulson") {
    if (is.null(adj_means))
      .dmar_refuse_argument(
        "'procedure = \"bryant_paulson\"' requires 'adj_means' (the simultaneous procedure works on the adjusted means, not on a single 'contrast').",
        arg = c("adj_means", "contrast"))
    bp_args <- list(adj_means = adj_means, s_ancova = s_ancova, n_per_group = n_per_group,
                    conf_level = conf_level)
    if (!missing(c_weights)) bp_args$c_weights <- c_weights
    return(do.call(ci_c_ancova_bp, c(bp_args, list(...))))
  }

  # All-zero weights define no contrast; they gave an interval of zero
  # width around the contrast ([0, 0] from the adjusted means, [.5, .5]
  # with contrast = .5).
  # The Bryant-Paulson path above is checked in ci_c_ancova_bp().
  .dmar_check_contrast_nonzero(c_weights)
  if (length(cov_means) != length(c_weights)) .dmar_refuse_argument(
    "The input 'cov_means' and 'c_weights' imply different number of groups",
    arg = c("cov_means", "c_weights"))

  if (is.null(contrast) && is.null(adj_means)) .dmar_refuse_argument(
    "Supply one of 'contrast' and 'adj_means'; neither was given.", arg = c("contrast", "adj_means"))
  if (!is.null(contrast) && !is.null(adj_means)) .dmar_refuse_argument(
    "Supply one of 'contrast' and 'adj_means', not both.",
    arg = c("contrast", "adj_means"))
  # The standard error of the contrast is s_ancova times
  # sqrt(sum(c_weights^2 / n_per_group) + (sum(c_weights * cov_means))^2 / SSwithin_x)
  # (Maxwell, Delaney, & Kelley, 2027, ch. 9), and s_ancova, the square
  # root of the ANCOVA mean square error, cannot be recovered from the
  # adjusted means and the sample sizes, so it is required on both the contrast
  # and the adj_means paths and must be a positive standard deviation
  # (.check_error_sd(), shared with ci_c()). Through 1.0.1 it was not
  # checked: omitted, each limit added contrast to a product with NULL, which is
  # numeric(0), and data.frame() recycled contrast into all three rows, a
  # zero-width interval at the estimate labeled 95%; 0 gave the same, and
  # a negative value swapped the limits.
  .check_error_sd(s_ancova, "s_ancova", "ANCOVA")

  if (is.null(contrast)) contrast <- sum(adj_means * c_weights)
  J <- length(c_weights)
  if (length(n_per_group) == 1) n_per_group <- rep(n_per_group, J)
  if (length(n_per_group) > 1 && length(n_per_group) != length(c_weights)) .dmar_refuse_argument(
    "The input 'n_per_group' and 'c_weights' imply different number of groups ",
    arg = c("n_per_group", "c_weights"))
  ########################################################################
  f_x_numerater <- (sum(c_weights * cov_means))^2
  f_x_denominator <- SSwithin_x
  sample_size_weighted <- sum(c_weights^2 / n_per_group)

  se_Psi <- s_ancova * sqrt(sample_size_weighted + f_x_numerater / f_x_denominator)

  alpha <- 1 - conf_level
  nu <- sum(n_per_group) - J - 1
  t_value <- qt(1 - alpha / 2, df = nu)

  term <- c("ci_lower", "contrast", "ci_upper")
  value <- c(contrast - t_value * se_Psi, contrast, contrast + t_value * se_Psi)

  out <- data.frame(term, value)
  .as_dmar_tbl(out, conf_level = conf_level)
}
