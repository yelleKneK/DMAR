# Confidence interval for the mean of a fully crossed units by raters design.
#' Confidence Interval for the Mean of a Crossed Design of Units and Raters
#'
#' @description
#' The confidence interval for the population mean score of a fully
#' crossed design, in which every rater scores every unit once, when the
#' raters, like the units, are a sample from a population. The units are
#' the objects of measurement (persons, essays, recordings, or the topics
#' of an information retrieval evaluation) and the raters are whoever or
#' whatever produces the scores (judges, coders, assessors, or automated
#' scoring models). Because a different panel of raters would have given
#' a different mean, the standard error carries the variability of the
#' raters as well as of the units; treating the raters as fixed, as an
#' interval computed from the unit means alone does, understates it
#' whenever raters differ. This is the interval
#' \code{\link{ss_aipe_crossed}} plans for, and its variance component
#' rows are the planning values that function takes.
#'
#' @param ratings A numeric matrix or \code{data.frame} with one row per
#'   unit and one column per rater, holding each rater's score of each
#'   unit. At least two units and two raters are needed, and every cell
#'   must hold a finite score: the design is complete, with no missing
#'   values. A cell may also hold a difference of two scores, as when
#'   every rater scores every unit under two conditions; the interval is
#'   then for the mean difference.
#' @param conf_level Confidence level of the interval (default
#'   \code{0.95}).
#'
#' @return A \code{data.frame} (class \code{dmar_tbl}) with rows
#'   \code{ci_lower}, \code{grand_mean}, and \code{ci_upper} (the
#'   estimate between its limits), \code{se} (the estimated standard
#'   error of the grand mean), \code{df} (its Satterthwaite degrees of
#'   freedom), the estimated variance components \code{var_unit},
#'   \code{var_rater}, and \code{var_residual} (the first two truncated
#'   at zero), and the counts \code{n_units} and \code{n_raters}.
#'
#' @details
#' \strong{Model.} The score of unit \eqn{u} by rater \eqn{r} is
#' \eqn{y_{ur} = \mu + a_u + b_r + e_{ur}}, with independent unit
#' effects, rater effects, and residuals of mean zero and variances
#' \eqn{\sigma^2_u}, \eqn{\sigma^2_r}, and \eqn{\sigma^2_e}; the residual
#' holds the unit by rater interaction, which one score per cell cannot
#' separate from error. This is the design Shrout and Fleiss (1979) call
#' Case 2, a random sample of judges each of whom rates every target,
#' and it is the crossed design with one random facet of generalizability
#' theory (Brennan, 2001). The grand mean \eqn{\bar y} of the
#' \eqn{n_u \times n_r} array estimates \eqn{\mu} with variance
#' \deqn{V = \frac{\sigma^2_u}{n_u} + \frac{\sigma^2_r}{n_r} +
#'   \frac{\sigma^2_e}{n_u n_r}.}
#'
#' \strong{Estimate and interval.} With the mean squares \eqn{MS_u},
#' \eqn{MS_r}, and \eqn{MS_e} of the two-way analysis of variance for
#' units, raters, and residual, the ANOVA estimates of the components are
#' \eqn{(MS_u - MS_e) / n_r}, \eqn{(MS_r - MS_e) / n_u}, and \eqn{MS_e},
#' each of the first two set to zero when it is negative. Plugged into
#' \eqn{V} they give
#' \deqn{\hat V = \frac{k_u MS_u + k_r MS_r + (1 - k_u - k_r) MS_e}{n_u n_r},}
#' with \eqn{k_u} equal to 1 when \eqn{MS_u > MS_e} and 0 otherwise, and
#' \eqn{k_r} likewise for raters. The interval is
#' \eqn{\bar y \pm t_{\nu} \sqrt{\hat V}}, where the degrees of
#' freedom \eqn{\nu} are those Satterthwaite (1946) gives for a linear
#' combination of independent mean squares,
#' \deqn{\nu = \frac{\left(\sum_k a_k MS_k\right)^2}{\sum_k (a_k MS_k)^2 / \nu_k},}
#' with the coefficients \eqn{a_k} of the numerator above and
#' \eqn{\nu_k} the degrees of freedom of each mean square:
#' \eqn{n_u - 1}, \eqn{n_r - 1}, and \eqn{(n_u - 1)(n_r - 1)}. The
#' truncation keeps \eqn{\hat V} at or above \eqn{MS_e / (n_u n_r)}, so it
#' is never negative. When few raters carry a large rater component,
#' \eqn{\nu} is near \eqn{n_r - 1} and the interval is wide: the data
#' then hold little information about how far the mean of another panel
#' of raters would lie from this one's.
#'
#' \strong{An information retrieval evaluation.} To compare two
#' retrieval systems over topics judged by several relevance assessors
#' (people, or automated judges such as language models), let each row
#' be a topic, each column an assessor, and each cell the difference
#' between the two systems' effectiveness scores on that topic computed
#' from that assessor's judgments. The interval is then for the mean
#' difference between the systems over the populations of topics and
#' assessors.
#'
#' \strong{Designs chosen from the data.} The interval assumes the
#' numbers of units and raters were fixed in advance, as
#' \code{\link{ss_aipe_crossed}} fixes them, and not chosen from the data,
#' as a two-stage sequential procedure chooses them from a pilot array.
#'
#' @references
#' Brennan, R. L. (2001). \emph{Generalizability theory}. Springer.
#'   \doi{10.1007/978-1-4757-3456-0}
#'
#' Satterthwaite, F. E. (1946). An approximate distribution of
#'   estimates of variance components. \emph{Biometrics Bulletin, 2}(6),
#'   110--114. \doi{10.2307/3002019}
#'
#' Shrout, P. E., & Fleiss, J. L. (1979). Intraclass correlations: Uses
#'   in assessing rater reliability. \emph{Psychological Bulletin,
#'   86}(2), 420--428. \doi{10.1037/0033-2909.86.2.420}
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso \code{\link{ss_aipe_crossed}} to plan the numbers of units
#'   and raters, \code{\link{icc}} for the intraclass correlations of the
#'   same array, and \code{\link{variance_components_mls}} for intervals
#'   on variance components.
#'
#' @family confidence intervals
#'
#' @keywords htest
#'
#' @examples
#' # Shrout and Fleiss (1979, Table 2): four judges each rated the same
#' # six targets.
#' ratings <- matrix(c( 9, 2, 5, 8,
#'                      6, 1, 3, 2,
#'                      8, 4, 6, 8,
#'                      7, 1, 2, 6,
#'                     10, 5, 6, 9,
#'                      6, 2, 4, 7),
#'                   nrow = 6, byrow = TRUE)
#' fit <- ci_crossed(ratings)
#' fit
#'
#' # The judges differ a great deal in leniency, and with only four of
#' # them the interval for the mean over the population of judges is wide.
#' # Treating these four judges as the only ones of interest ignores that
#' # variability; the interval from the six target means is much narrower.
#' t.test(rowMeans(ratings))$conf.int
#'
#' # The variance components are planning values for a larger study, here
#' # one whose 95% interval for the mean is no wider than 2 points.
#' v <- setNames(fit$value, fit$term)
#' ss_aipe_crossed(true_var_unit = v[["var_unit"]], true_var_rater = v[["var_rater"]],
#'                 true_var_residual = v[["var_residual"]], width = 2)
#'
#' @export
ci_crossed <- function(ratings, conf_level = 0.95) {
  .crossed_check_number(conf_level, "conf_level", lower = 0, upper = 1,
                        open = TRUE)
  # SPSS user-defined missing values are missing values (as icc() and
  # fleiss_kappa() read them), which the completeness check below then
  # refuses, rather than scores.
  ratings <- .relia_declared_missing(ratings)
  if (is.data.frame(ratings)) {
    if (!all(vapply(ratings, is.numeric, logical(1))))
      .dmar_refuse_argument(c(
        "Every column of 'ratings' must be numeric: one column per ",
        "rater, one row per unit."), arg = "ratings")
    ratings <- as.matrix(ratings)
  }
  if (!is.matrix(ratings) || !is.numeric(ratings))
    .dmar_refuse_argument(c(
      "'ratings' must be a numeric matrix or data.frame with one row ",
      "per unit and one column per rater."), arg = "ratings")
  if (nrow(ratings) < 2L || ncol(ratings) < 2L)
    .dmar_refuse_argument(c(
      "'ratings' must have at least 2 rows (units) and 2 columns ",
      "(raters)."), arg = "ratings")
  if (anyNA(ratings))
    .dmar_refuse_argument(c(
      "'ratings' has missing values. The interval needs a complete ",
      "design, every rater scoring every unit."), arg = "ratings")
  if (!all(is.finite(ratings)))
    .dmar_refuse_argument(
      "Every score in 'ratings' must be finite.", arg = "ratings")

  ms <- .crossed_ms(ratings)
  if (!(ms[["ms_residual"]] > 0) && !(ms[["ms_unit"]] > 0) &&
      !(ms[["ms_rater"]] > 0))
    .dmar_refuse_argument(c(
      "The scores in 'ratings' do not vary, so the variance of their ",
      "mean cannot be estimated."), arg = "ratings")
  ci <- .crossed_full_ci(ms[["grand_mean"]], ms[["ms_unit"]],
                         ms[["ms_rater"]], ms[["ms_residual"]],
                         ms[["n_units"]], ms[["n_raters"]], conf_level)
  .as_dmar_tbl(data.frame(
    term  = c("ci_lower", "grand_mean", "ci_upper", "se", "df",
              "var_unit", "var_rater", "var_residual", "n_units",
              "n_raters"),
    value = c(ci$ci_lower, ci$grand_mean, ci$ci_upper, ci$se, ci$df,
              ci$var_unit, ci$var_rater, ci$var_residual,
              ms[["n_units"]], ms[["n_raters"]]),
    stringsAsFactors = FALSE
  ), conf_level = conf_level)
}
