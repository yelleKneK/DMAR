#' Bifactor Model Dimensionality and Reliability Indices
#'
#' Computes the indices that describe how much of the common variance of
#' a scale's items, and of the variance of its total and subscale scores,
#' a general factor accounts for beside the group factors of its
#' subscales (Rodriguez, Reise, and Haviland, 2016): the explained common
#' variance (ECV), coefficient omega and omega hierarchical (omega_H) for
#' the general factor, omega hierarchical subscale (omega_HS) for each
#' group factor, the percentage of uncontaminated correlations (PUC), and
#' the index those authors call coefficient H, the closed form of Hancock and
#' Mueller (2001) applied to one factor's loadings, which counts the other
#' factors' common variance as uncorrelated error and so is not the
#' maximal reliability of the factor under the fitted model; that maximal
#' reliability is returned beside the table. The input is a fitted
#' bifactor model in which one general factor loads on every item and
#' each item loads on at most one of two or more orthogonal group
#' factors. \code{\link{fa_bifactor}} fits such a model and reports these
#' indices with delta method intervals, with omega, omega_H, and omega_HS
#' for the sum of the items as scored rather than of the standardized
#' items (see \emph{Details}).
#'
#' @param fit A fitted bifactor \pkg{lavaan} model: one general factor on
#'   all items plus two or more orthogonal group factors, each item on one
#'   group factor at most. An item on no group factor is on the general
#'   factor alone: it counts in the general factor's indices and in PUC,
#'   and in no subscale. An item listed under two or more group factors is
#'   refused
#'   unless all but one of those loadings are held at zero, fixed
#'   (\code{0*item}) or constrained (\code{a*item} with \code{a == 0});
#'   such a loading describes the same model as leaving the item off that
#'   factor. The factors must be orthogonal (the bifactor specification),
#'   and a fit that correlates any two of them is refused. Residual
#'   covariances among the items are allowed and enter the model implied
#'   denominators. A fit with fewer than two group factors, a
#'   multiple-group fit, a multilevel fit, and a fit whose optimizer did
#'   not converge are refused (see \emph{Details}).
#' @param general Optional name of the general factor. When \code{NULL}
#'   (default) the general factor is detected as the one that loads on every
#'   item.
#' @param denominator How the total-score variance in the denominators of
#'   omega, omega_H, and omega_HS is estimated: \code{"model_implied"}
#'   (the default; the variance the fitted bifactor model reproduces, as
#'   in Rodriguez, Reise, and Haviland, 2016, with any residual
#'   covariances the model specifies) or \code{"observed"} (the variance
#'   of the composite estimated from the data, the denominator of robust
#'   omega in \code{\link{reliability_omega}}). For a fit to ordered
#'   items both describe the sum of the latent responses, not of the
#'   summed categories. See \emph{Details}.
#'
#' @details
#' Let \eqn{\lambda^g_i} be item \eqn{i}'s standardized loading on the
#' general factor, \eqn{\lambda^s_i} its loading on its group factor (0
#' for an item on the general factor alone), and
#' \eqn{\theta_i} its standardized residual variance. With orthogonal
#' factors the overall ECV is
#' \eqn{\sum_i (\lambda^g_i)^2 / \sum_i [(\lambda^g_i)^2 + (\lambda^s_i)^2]};
#' omega and omega_H share the total-score variance
#' \eqn{(\sum_i \lambda^g_i)^2 + \sum_g (\sum_{i \in g} \lambda^s_i)^2 +
#' \sum_i \theta_i} as denominator, with the general part
#' \eqn{(\sum_i \lambda^g_i)^2} in the numerator of omega_H. Each group
#' factor's omega_HS uses the analogous numerator and that subscale's own
#' total variance. PUC is one minus the share of item pairs that fall within
#' the same group factor. Omega_H against omega says how much of the
#' total score's reliable variance is due to the general factor; omega_HS
#' against a subscale's omega says how much of the subscale score's
#' reliable variance is its own, beyond the total; ECV and PUC describe
#' how close the items come to measuring one thing. The table reports the
#' numbers; whether a total or subscale score is worth reporting is the
#' researcher's judgment, and the package draws no line, as on the page
#' of \code{\link{fa_bifactor}}.
#'
#' \strong{The metric of omega.} The loadings and residual variances above
#' are standardized, so omega, omega_H, and omega_HS describe the sum of
#' the standardized items, as Rodriguez, Reise, and Haviland (2016) write
#' them. \code{\link{fa_bifactor}} reports them for the sum of the items as
#' scored, the coefficients \code{\link{fa_confirmatory}} and
#' \code{semTools::compRelSEM()} compute. The two agree when the items have
#' equal variances and differ otherwise: in the first example below, whose
#' items have nearly equal variances, omega_H is 0.67348 here and 0.67305
#' for the sum of the items as scored, while on the 15 tests of the example
#' of \code{\link{fa_bifactor}}, whose variances range from 8.0 to 153.3,
#' it is 0.577 here and 0.552 there. ECV, \emph{H}, the maximal
#' reliability, and PUC do not differ.
#'
#' \strong{A group factor's row.} The row of group factor \eqn{s}
#' describes its items, the subscale \eqn{s}, and its ECV and omega are
#' not the general row's quantities restricted to a group factor. Its ECV
#' is the general factor's share of the subscale's common variance,
#' \deqn{\frac{\sum_{i \in s} (\lambda^g_i)^2}
#'   {\sum_{i \in s} [(\lambda^g_i)^2 + (\lambda^s_i)^2]},}
#' so one minus it is the group factor's own share of that variance, and
#' a high value says the subscale's common variance is mostly general
#' (it goes with a low omega_HS). It is not the group factor's share of
#' the common variance of all the items, which is
#' \eqn{\sum_{i \in s} (\lambda^s_i)^2} over the denominator of the
#' overall ECV. Its omega is the omega of the subscale's total score,
#' \eqn{(\sum_{i \in s} \lambda^g_i)^2 + (\sum_{i \in s} \lambda^s_i)^2}
#' over that score's variance, the reliability of the subscale score with
#' both factors counted as true score; omega_HS keeps only the group
#' factor's part of the numerator. Its H is the index described next,
#' computed from the group factor's loadings. The printed table states
#' what ECV and omega mean on these rows.
#'
#' \strong{Coefficient H.} The column \code{H} is the index Rodriguez,
#' Reise, and Haviland (2016) report: the closed form of Hancock and
#' Mueller (2001),
#' \deqn{H = \frac{\sum_i \lambda_i^2 / (1 - \lambda_i^2)}
#'   {1 + \sum_i \lambda_i^2 / (1 - \lambda_i^2)},}
#' over the general loadings for the general factor and over a group
#' factor's loadings for that factor. It counts \eqn{1 - \lambda_i^2} as
#' item \eqn{i}'s error variance, and so treats the variance an item
#' shares with the other factors as error that is uncorrelated across
#' items. In a bifactor model that variance is correlated (the items of
#' a subscale share their group factor, and every item shares the general
#' factor), so the index is not the maximal reliability of the factor
#' under the fitted model, which is what \code{\link{reliability_H}}
#' computes for a block of items whose errors are correlated (Gabler &
#' Raykov, 2017).
#' That maximal reliability is \eqn{\lambda_f' R^{-1} \lambda_f}, where
#' \eqn{\lambda_f} holds the factor's standardized loadings (zero for an
#' item off the factor) and \eqn{R} is the item correlation matrix the
#' model reproduces: the squared multiple correlation of the factor with
#' all of the items, the largest reliability any weighted sum of them
#' reaches as a measure of that factor, and the square of the determinacy
#' of its regression factor score in \code{\link{factor_scores}}. It is
#' returned as the \code{"maximal_reliability"} attribute and printed
#' under the table, and \code{\link{ci_bifactor_indices}} gives it a
#' confidence interval, as it does every index of the table but PUC,
#' which depends on the model's structure alone.
#'
#' The two differ in either direction, and what sets the direction is
#' not the same on the two kinds of row. On the general factor's row the
#' index counts the group factors' variance as uncorrelated error, but
#' the items of a subscale share it. When the loadings on each group
#' factor are positive and similar, no weighting of the items removes
#' that shared variance without removing general variance with it, and
#' the index is the larger: in the first example below it is 0.769 for
#' the general factor, whose maximal reliability is 0.678, barely above
#' its omega_H of 0.673, so no weighting of the items reaches 0.769. When
#' an item is nearly a pure indicator of its group factor, or a group
#' factor's loadings differ in sign, a weighting can cancel group
#' variance that the index counts as error, and the maximal reliability
#' is the larger. On a group factor's row the index counts the general
#' factor's variance, with the unique variance, as uncorrelated error.
#' When the subscale's own items have proportional loadings on the two
#' factors (equal loadings, for instance), no weighting of them alone
#' removes the general factor without removing the group factor with it,
#' but the optimal weights use every item, and the items of the other
#' subscales measure the general factor and nothing of this group factor,
#' so the weights use them to remove it. The better those items measure
#' the general factor, the closer the maximal reliability comes to the
#' closed form above with each item's unique variance \eqn{\theta_i} in
#' place of \eqn{1 - (\lambda^s_i)^2}, which is above the index. Which of
#' the two prevails depends on the general loadings, the group loadings,
#' and the number of items outside the subscale together, not on the
#' sign and similarity of the group loadings alone. In a population of
#' nine items in three subscales of three, with every general loading .7
#' and every group loading .3 (an ECV of .84), the index is 0.229 for
#' each group factor and its maximal reliability 0.279, while for the
#' general factor they are 0.896 and 0.865; with general and group
#' loadings of .5, as in the first example below, the index is the larger
#' on every row (0.500 against 0.471 for a group factor in the
#' population); and with two subscales of six items, general loadings
#' .8, and group loadings .4, the index is the larger for a group factor
#' (0.533 against 0.468), although the general factor is stronger than in
#' the first population.
#'
#' The composite that reaches a group factor's maximal reliability is the
#' factor's regression score from all of the items, and it weights the
#' items of the other subscales negatively: in the first population the
#' optimal weights are proportional to 1 for each of the subscale's three
#' items and -0.405 for each of the other six. No score of the subscale's
#' own items reaches that reliability. The best weighting of them alone
#' reaches \eqn{\lambda_s' R_{ss}^{-1} \lambda_s}, with \eqn{\lambda_s}
#' their group loadings and \eqn{R_{ss}} their correlations, which is
#' omega_HS when the items have equal loadings on both factors, and then
#' below the index as well: 0.125 in the first population, against an
#' index of 0.229 and a maximal reliability of 0.279. So a maximal
#' reliability above the index on a group factor's row, which is common
#' when the items outside the subscale measure a strong general factor
#' well (a profile of high ECV, in which the subscale rows are commonly
#' consulted), describes a factor score computed from all of the items,
#' not the subscale score; how much of a subscale score's variance is due
#' to its group factor is omega_HS. Read the column as the published
#' index, comparable with values reported beside it, and the attribute as
#' the reliability that optimal weights of all of the items could reach
#' in the population. The sample values of both overstate their
#' population values on average, the attribute most of all for a weakly
#' determined group factor, and a composite weighted with sample estimates
#' reaches less, as Aguirre-Urreta, Rönkkö, and McIntosh (2019) found for
#' H of a single factor; \code{\link{ci_bifactor_indices}} gives the
#' bootstrap estimate of that bias with the intervals.
#'
#' The H of a factor any of whose standardized loadings is at or beyond
#' one in absolute value is \code{NA}, since the index then gives that
#' item an error variance that is not positive.
#'
#' Each denominator is the sum of the corresponding block of an item
#' correlation matrix on the standardized metric of the loadings: all
#' items for omega and omega_H, a group factor's items for its omega_HS.
#' The default (\code{denominator = "model_implied"}) takes the
#' correlation matrix the fitted model reproduces. For an orthogonal
#' bifactor model without residual covariances that block sum is the
#' total-score variance written above, the form given by Rodriguez,
#' Reise, and Haviland (2016); when the model also specifies residual
#' covariances among the items, the block sum carries each of them as
#' well, since they are part of the variance of the composite the model
#' reproduces, although not of its common variance, so the numerators
#' leave them out. With \code{denominator = "observed"} the matrix is
#' instead the item correlation matrix estimated from the data
#' (\pkg{lavaan}'s sample statistics, which under \code{missing = "ml"}
#' are the full information estimates). This is the denominator of
#' robust omega (\code{\link{reliability_omega}}) carried to the
#' bifactor coefficients: it estimates the variance of the composite
#' consistently whether or not the bifactor model reproduces it, while
#' the numerators remain model based. When the model reproduces the item
#' correlations the two settings agree, and a gap between them shows
#' that it does not reproduce the composite variance. The printed table
#' names the denominator in a line under the table.
#'
#' The indices assume uncorrelated factors: a correlation between the
#' general factor and a group factor, or between two group factors, is
#' common variance that belongs to neither, and the partition behind
#' ECV, omega_H, and omega_HS has no place for it. A fit that estimates
#' such a correlation is refused with an error naming the pair.
#'
#' \strong{A fit with one group factor.} A general factor beside a single
#' group factor is not a bifactor model with subscales: its indices would
#' divide the common variance between those two factors alone, every item
#' would belong to the one subscale, and PUC would be 0. Such a fit is
#' refused, as \code{\link{fa_bifactor}} and
#' \code{\link{ci_bifactor_indices}} refuse it, with a message that names
#' the group factor. A general factor with one
#' further factor can be fitted with \code{lavaan::cfa()} and read from
#' its own parameter estimates.
#'
#' \strong{Ordered-categorical items.} For a fit with ordered
#' indicators (WLSMV or DWLS on polychoric correlations), \pkg{lavaan}'s
#' sample statistics are the polychoric correlations, and the
#' standardized loadings are on the metric of the continuous latent
#' responses behind the categories. Both settings then give the
#' reliability of the sum of those latent responses, not of the summed
#' category scores a user computes, and categorization lowers
#' reliability, often substantially (Green & Yang, 2009, Table 1). The
#' printed line says which metric the denominator is on, and the
#' \code{"metric"} attribute records it. For a single factor,
#' \code{\link{reliability_omega_categorical}} gives the reliability of
#' the summed categories; \pkg{DMAR} does not compute the bifactor
#' indices on that metric.
#'
#' A bifactor model that is over-parameterized for the data frequently
#' yields an improper (Heywood) solution, a standardized loading outside
#' \eqn{[-1, 1]} or a negative residual variance, whose indices are not
#' trustworthy. When that happens the indices are still returned but a
#' warning is issued and the \code{"improper"} attribute is set to
#' \code{TRUE}; H is \code{NA} for a factor with a loading at or beyond
#' one, as described above, and the maximal reliability is \code{NA} for
#' every factor, since it is not defined for an improper solution.
#'
#' A fit whose optimizer did not converge (\code{lavaan::lavInspect(fit,
#' "converged")} is \code{FALSE}, as when it reached its iteration limit)
#' is not a solution, proper or improper: its estimates mark where the
#' optimizer stopped. It is refused with an error of class
#' \code{dmar_not_converged}, as \code{\link{ci_bifactor_indices}}
#' refuses it.
#'
#' \strong{One group and one level.} The standardized loadings and
#' residual variances of a fit are those of one group at one level: each
#' group of a multiple-group fit has its own, standardized by that group's
#' variances, so the factors have indices in each group, and neither an
#' unweighted mean nor a mean weighted by group size is the index of any
#' group. A multiple-group fit is refused, as
#' \code{\link{average_variance_extracted}} refuses it, with an error of
#' class \code{dmar_argument_invalid} that names the groups (as
#' \code{groups}); the indices of one group are those of the model fitted
#' to that group's data alone (which drops any equality constraints
#' across groups). A multilevel fit, whose
#' loadings and residual variances differ by level, is refused with an
#' error of the same class; the standardized solution of each level is in
#' \code{lavaan::lavInspect(fit, "std")}.
#'
#' @return A \code{data.frame} (class \code{dmar_tbl}) with one row per
#'   factor (the general factor first, then each group factor) and columns
#'   \code{factor}, \code{ECV}, \code{omega}, \code{omega_H},
#'   \code{omega_HS}, \code{PUC}, and \code{H}. Quantities that do not apply
#'   to a row are \code{NA} (for example omega_H and PUC on a group factor).
#'   On a group factor's row, ECV is the general factor's share of the
#'   common variance of that factor's items and omega is the reliability
#'   of their total score (see \emph{Details}); H is the index of
#'   Rodriguez, Reise, and Haviland (2016) on every row.
#'   The \code{"improper"} attribute flags a Heywood solution, the
#'   \code{"maximal_reliability"} attribute holds the maximal reliability
#'   the fitted model implies for each factor (a vector named by factor,
#'   in the order of the rows; \code{NA} for an improper solution), the
#'   \code{"denominator"} attribute records the \code{denominator}
#'   setting, and the \code{"metric"} attribute records the metric of
#'   the items the indices describe: \code{"observed"} for continuous
#'   items, \code{"latent_response"} for ordered items, whose
#'   coefficients describe the sum of the latent responses.
#'
#' @references
#' Aguirre-Urreta, M. I., Rönkkö, M., & McIntosh, C. N. (2019). A
#'   cautionary note on the finite sample behavior of maximal reliability.
#'   \emph{Psychological Methods, 24}(2), 236--252.
#'   \doi{10.1037/met0000176}
#'
#' Gabler, S., & Raykov, T. (2017). Evaluation of maximal reliability for
#'   unidimensional measuring instruments with correlated errors.
#'   \emph{Structural Equation Modeling, 24}(1), 104--111.
#'   \doi{10.1080/10705511.2016.1159916}
#'
#' Green, S. B., & Yang, Y. (2009). Reliability of summed item scores
#'   using structural equation modeling: An alternative to coefficient
#'   alpha. \emph{Psychometrika, 74}, 155--167.
#'   \doi{10.1007/s11336-008-9099-3}
#'
#' Hancock, G. R., & Mueller, R. O. (2001). Rethinking construct
#'   reliability within latent variable systems. In R. Cudeck, S. du Toit,
#'   & D. Sörbom (Eds.), \emph{Structural equation modeling: Present and
#'   future} (pp. 195--216). Scientific Software International.
#'
#' Reise, S. P. (2012). The rediscovery of bifactor measurement models.
#'   \emph{Multivariate Behavioral Research, 47}(5), 667--696.
#'   \doi{10.1080/00273171.2012.715555}
#'
#' Rodriguez, A., Reise, S. P., & Haviland, M. G. (2016). Evaluating
#'   bifactor models: Calculating and interpreting statistical indices.
#'   \emph{Psychological Methods, 21}(2), 137--150.
#'   \doi{10.1037/met0000045}
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso \code{\link{fa_bifactor}}, which fits the bifactor model and
#'   gives these indices delta method intervals;
#'   \code{\link{ci_bifactor_indices}} for bootstrap confidence
#'   intervals for these indices and for the maximal reliability of each
#'   factor; \code{\link{reliability_omega}},
#'   \code{\link{reliability_H}}, and \code{\link{factor_scores}}, whose
#'   regression score determinacy is the square root of the maximal
#'   reliability.
#'
#' @family multivariate and latent variable methods
#'
#' @keywords multivariate
#'
#' @examples
#' # Nine items: one general factor and three orthogonal group factors.
#' set.seed(113)
#' n <- 600
#' g <- rnorm(n); grp <- list(rnorm(n), rnorm(n), rnorm(n))
#' X <- vapply(1:9, function(i)
#'   0.5 * g + 0.5 * grp[[ceiling(i / 3)]] + sqrt(0.5) * rnorm(n), numeric(n))
#' colnames(X) <- paste0("x", 1:9)
#' model <- "g  =~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9
#'           f1 =~ x1 + x2 + x3
#'           f2 =~ x4 + x5 + x6
#'           f3 =~ x7 + x8 + x9"
#' fit <- lavaan::cfa(model, data = as.data.frame(X),
#'                    orthogonal = TRUE, std.lv = TRUE)
#' bifactor_indices(fit)
#'
#' # The same indices with the observed composite variances in the
#' # denominators; when the model reproduces the item correlations the
#' # two settings agree.
#' bifactor_indices(fit, denominator = "observed")
#'
#' @export
bifactor_indices <- function(fit, general = NULL,
                             denominator = c("model_implied", "observed")) {
  denominator <- match.arg(denominator)
  .bifactor_indices(fit, general, denominator)
}

# .bifactor_indices(fit, general, denominator)
#
# The body of bifactor_indices(), shared with ci_bifactor_indices().
# Returns the table.
.bifactor_indices <- function(fit, general, denominator) {
  if (!requireNamespace("lavaan", quietly = TRUE)) {
    .dmar_refuse(c(
      "Package 'lavaan' is required. Install it with ",
      "install.packages(\"lavaan\")."), "package_required", package = "lavaan")
  }
  if (!inherits(fit, "lavaan")) {
    .dmar_refuse_argument("'fit' must be a lavaan fit object.", arg = "fit")
  }
  # A fit whose optimizer did not converge holds the point where it
  # stopped, not a solution (.lavaan_check_converged(),
  # R/no_maximum_internals.R), and no index is read from it. Through DMAR
  # 1.0.1 and the first builds of 1.1.0 its table was returned as any
  # other: a nine-item fit stopped after one iteration gave the general
  # factor an ECV of .558, an omega of .851, and an H of .764, where the
  # converged fit gives .610, .810, and .743. ci_bifactor_indices() makes
  # the same check with the same helper before its bootstrap, and the
  # refits of its bootstrap are checked where they are made
  # (.bifactor_point_fn()).
  .lavaan_check_converged(fit, "bifactor_indices()")
  # A fit of one group, refused otherwise with the class and the wording
  # of average_variance_extracted() (.ave_check_one_block(),
  # R/average_variance_extracted.R): each group has its own standardized
  # loadings and residual variances, so each has its own indices, and no
  # pooled value is the index of any group. Through DMAR 1.0.1 the table
  # read the first group's loadings and the last group's residual
  # variances (for nine items in groups a and b of 200, a general factor
  # omega_H of .623, where group a alone gives .621 and group b alone
  # .674); the first builds of 1.1.0 described the first group with a
  # message.
  n_groups <- lavaan::lavInspect(fit, "ngroups")
  if (n_groups > 1L) {
    labels <- as.character(lavaan::lavInspect(fit, "group.label"))
    .dmar_refuse_argument(c(
      "'fit' is a multiple-group fit (", n_groups, " groups: ",
      paste(labels, collapse = ", "), "). Each group has its own ",
      "standardized loadings and residual variances, so the factors have ",
      "bifactor indices in each group, and neither their unweighted mean ",
      "nor a mean weighted by group size is the index of any group; ",
      "bifactor_indices() describes a fit of one group. For one group of ",
      "this model, fit the model to that group's data alone (which drops ",
      "any equality constraints across groups)."),
      arg = "fit", groups = labels)
  }
  # A two-level fit has loadings, residual variances, and implied
  # correlations at each level, and the indices of one level are not
  # those of the other. Through DMAR 1.0.1 the table of such a fit read
  # the within level's loadings and the between level's residual
  # variances (a general factor omega of .942 where the within level's
  # own solution gives .885); the first builds of 1.1.0 stopped with an
  # error of R's indexing ("incorrect number of dimensions") or with the
  # refusal of a fit that has no unique general factor.
  n_levels <- lavaan::lavInspect(fit, "nlevels")
  if (n_levels > 1L) {
    .dmar_refuse_argument(c(
      "'fit' is a multilevel fit (", n_levels, " levels: ",
      paste(lavaan::lavInspect(fit, "level.label"), collapse = ", "),
      "). Its loadings, residual variances, and implied correlations ",
      "differ by level, and bifactor_indices() describes a fit of one ",
      "level. The standardized solution of each level is in ",
      "lavaan::lavInspect(fit, \"std\")."), arg = "fit")
  }
  .fiml_require_h1(fit, "These bifactor indices",
                    observed = denominator == "observed", groups = 1L)
  # Only est.std is read, so the delta method standard errors, tests, and
  # intervals of standardizedSolution() are not computed. Through DMAR
  # 1.0.1 and the first builds of 1.1.0 they were, and for a fit whose
  # information matrix lavaan cannot invert (common with an improper
  # bifactor solution) the call passed on lavaan's warning that it "Could
  # not compute standard errors", about standard errors this function does
  # not report.
  std <- lavaan::standardizedSolution(fit, se = FALSE, zstat = FALSE,
                                      pvalue = FALSE, ci = FALSE)
  lam <- std[std$op == "=~", c("lhs", "rhs", "est.std"), drop = FALSE]
  if (nrow(lam) == 0L) {
    .dmar_refuse_argument(
      "The fit contains no measurement (=~) loadings.", arg = NULL)
  }
  items <- unique(lam$rhs)
  factors <- unique(lam$lhs)
  pt <- lavaan::parTable(fit)
  pt_g1 <- if (is.null(pt$group)) TRUE else pt$group == 1L
  zero <- which(pt$op == "=~" & .pt_zero_rows(pt) & pt_g1)
  active <- lam[!paste(lam$lhs, lam$rhs) %in%
                  paste(pt$lhs[zero], pt$rhs[zero]), , drop = FALSE]

  # Identify the general factor: the one loading on every item.
  loads_all <- vapply(factors, function(f) all(items %in% active$rhs[active$lhs == f]),
                      logical(1))
  if (is.null(general)) {
    if (sum(loads_all) != 1L) {
      .dmar_refuse_argument(c(
        "Could not identify a unique general factor (one factor loading ",
        "on every item, with no loading fixed or constrained to zero). ",
        "Name it with 'general' if the model has one. A model outside ",
        "this bifactor definition can still be fitted with lavaan::cfa() ",
        "or lavaan::sem(); inspect its ordinary parameters and fit rather ",
        "than these bifactor indices."), arg = "general")
    }
    general <- factors[loads_all]
  } else if (!general %in% factors || !all(items %in% active$rhs[active$lhs == general])) {
    .dmar_refuse_argument(c(
      "'general' must name a factor that loads on every item; a loading ",
      "fixed or constrained to zero does not count. The model can still ",
      "be fitted with lavaan::cfa() or lavaan::sem(), but these bifactor ",
      "indices do not apply."), arg = "general")
  }
  # Two or more group factors. A general factor beside a single group
  # factor is not a bifactor model with subscales, and a fit with none or
  # with one is refused. (DMAR 1.0.0 refused only a fit with none, and
  # described a fit with one when every item was listed under it.)
  groups <- setdiff(factors, general)
  if (!length(groups)) {
    .dmar_refuse_argument(c(
      "A bifactor model has a general factor plus two or more group ",
      "factors, and this fit has none. The indices partition the ",
      "common variance among group factors, which takes two or more."),
      arg = NULL)
  }
  if (length(groups) == 1L) {
    off_group <- setdiff(items, lam$rhs[lam$lhs == groups])
    .dmar_refuse_argument(c(
      "A bifactor model has a general factor plus two or more group ",
      "factors, and this fit has one (", groups, ")",
      if (length(off_group)) {
        paste0(", with items on no group factor (",
               paste(off_group, collapse = ", "), ")")
      },
      ". A general factor (", general, ") beside a single group ",
      "factor is not a bifactor model with subscales: the indices ",
      "would divide the common variance between those two factors ",
      "alone. Fit a model with two or more group factors ",
      "(fa_bifactor() writes one from a list of the subscales), or ",
      "read this fit from its own parameter estimates."), arg = NULL)
  }

  # The indices partition the common variance into general and group
  # parts, which holds only when the factors are uncorrelated (the
  # bifactor specification). A correlation between two of them is
  # common variance that belongs to neither, so refuse the fit rather
  # than report indices that leave it out.
  cor_lv <- lavaan::lavInspect(fit, "cor.lv")
  cor_lv <- unclass(cor_lv)[factors, factors, drop = FALSE]
  off <- which(upper.tri(cor_lv) & abs(cor_lv) > sqrt(.Machine$double.eps),
               arr.ind = TRUE)
  if (nrow(off)) {
    pairs <- paste0(factors[off[, 1]], " and ", factors[off[, 2]], " (",
                    formatC(cor_lv[off], format = "f", digits = 3), ")")
    .dmar_refuse_argument(c(
      "bifactor_indices() needs orthogonal general and group factors ",
      "(the bifactor specification), but the fit correlates ",
      paste(pairs, collapse = ", "), ". Refit with ",
      "orthogonal = TRUE in lavaan::cfa(), or fix those covariances ",
      "to zero."), arg = NULL)
  }

  # Per-item general loading, group loading, group membership. Each item
  # belongs to at most one group factor, which defines the subscales
  # behind omega_HS and PUC and the partition behind ECV. An item listed
  # under two or more group factors keeps only those of its loadings that
  # the model does not hold at zero (.pt_zero_rows(): 0*x, or a*x with
  # a == 0, either of which describes the same model as leaving the item
  # off that factor); an item with more than one such loading is refused.
  # Through DMAR 1.0.1 the last group factor listed silently took the
  # item, and its other loading dropped out of every index. An item with a
  # single group loading keeps it, whatever its value. An item with none
  # is on the general factor alone: its group loading is 0, it belongs to
  # no subscale, and it counts in the general factor's indices and in PUC
  # (its pairs with every other item are uncontaminated).
  lg <- stats::setNames(lam$est.std[lam$lhs == general][match(items, lam$rhs[lam$lhs == general])], items)
  glam <- lam[lam$lhs %in% groups, , drop = FALSE]
  n_listed <- table(factor(glam$rhs, levels = items))
  if (any(n_listed > 1L)) {
    pt <- lavaan::parTable(fit)
    pt_g1 <- if (is.null(pt$group)) TRUE else pt$group == 1L
    zero <- which(pt$op == "=~" & .pt_zero_rows(pt) & pt_g1)
    held_zero <- paste(glam$lhs, glam$rhs) %in%
      paste(pt$lhs[zero], pt$rhs[zero])
    glam <- glam[!(held_zero & glam$rhs %in% names(n_listed)[n_listed > 1L]),
                 , drop = FALSE]
  }
  on_groups <- split(glam$lhs, factor(glam$rhs, levels = items))
  several <- items[lengths(on_groups) > 1L]
  if (length(several)) {
    .dmar_refuse_argument(c(
      "Every item must load on one group factor at most, but ",
      paste0(several, " loads on ",
             vapply(on_groups[several], paste, character(1),
                    collapse = " and "),
             collapse = "; "),
      ". The bifactor indices assign each item to one subscale (the ",
      "basis of omega_HS, PUC, and the ECV partition), and an item on ",
      "two group factors belongs to none of them alone. Respecify the ",
      "model with each item on one group factor (a loading fixed or ",
      "constrained to zero does not count)."), arg = NULL)
  }
  ls <- stats::setNames(rep(0, length(items)), items)
  grp <- stats::setNames(rep(NA_character_, length(items)), items)
  ls[glam$rhs] <- glam$est.std
  grp[glam$rhs] <- glam$lhs

  # Standardized residual variances (fall back to 1 - communality).
  rv <- std[std$op == "~~" & std$lhs == std$rhs & std$lhs %in% items, c("lhs", "est.std")]
  theta <- stats::setNames(1 - lg^2 - ls^2, items)
  if (nrow(rv)) theta[rv$lhs] <- rv$est.std

  improper <- any(abs(lg) > 1) || any(abs(ls) > 1) || any(theta < 0)
  if (improper) {
    .dmar_warn(c(
      "Improper (Heywood) solution: a standardized loading exceeds one ",
      "or a residual variance is negative. The indices are not ",
      "trustworthy until the model is respecified."), "heywood_warning")
  }

  # Coefficient H in the form of Rodriguez, Reise, and Haviland (2016):
  # Hancock and Mueller's closed form on one factor's standardized
  # loadings, which counts 1 - l^2 as each item's error variance (the
  # other factors' common variance included). A loading at or beyond one
  # in absolute value leaves an item an error variance that is not
  # positive, and H is then not defined, as reliability_H() treats an
  # improper solution; it is NA. Through DMAR 1.0.1 such a loading was
  # left out, so H described the factor's other items, and fell as the
  # loading crossed one where it should have risen toward one.
  H_of <- function(l) {
    if (!length(l) || any(abs(l) >= 1)) return(NA_real_)
    t <- sum(l^2 / (1 - l^2)); t / (1 + t)
  }
  choose2 <- function(k) k * (k - 1) / 2

  # The denominators are sums of blocks of an item correlation matrix on
  # the standardized metric of the loadings: under the model implied
  # setting the matrix the fitted model reproduces (which for an
  # orthogonal bifactor model without residual covariances is the sum
  # of the squared summed loadings and the residual variances, and
  # otherwise also carries every residual covariance the model
  # specifies), and under the observed setting the matrix estimated
  # from the data (lavaan's sample statistics: the full information
  # estimates under missing = "ml", the polychoric correlations for
  # ordered items).
  R_den <- if (denominator == "observed") {
    lavaan::lavInspect(fit, "sampstat")$cov
  } else {
    lavaan::lavInspect(fit, "cor.ov")
  }
  missing_items <- setdiff(items, rownames(R_den))
  if (length(missing_items)) {
    .dmar_refuse_argument(c(
      if (denominator == "observed") "The sample statistics" else
        "The model implied correlations",
      " of 'fit' do not include ",
      paste(missing_items, collapse = ", "),
      "; every factor must load on observed items."), arg = "fit")
  }
  R_den <- stats::cov2cor(unclass(R_den)[items, items, drop = FALSE])
  block_var <- function(sel) sum(R_den[sel, sel])

  # Overall (general-factor) row.
  total_common <- sum(lg)^2 + sum(vapply(groups, function(g) sum(ls[grp %in% g])^2, numeric(1)))
  total_var <- block_var(items)
  ecv_overall <- sum(lg^2) / (sum(lg^2) + sum(ls^2))
  k <- length(items)
  within_pairs <- sum(vapply(groups, function(g) choose2(sum(grp %in% g)), numeric(1)))
  puc <- if (choose2(k) > 0) 1 - within_pairs / choose2(k) else NA_real_

  rows <- list(data.frame(
    factor = general, ECV = ecv_overall,
    omega = if (total_var > 0) total_common / total_var else NA_real_,
    omega_H = if (total_var > 0) sum(lg)^2 / total_var else NA_real_,
    omega_HS = NA_real_, PUC = puc, H = H_of(lg),
    stringsAsFactors = FALSE
  ))
  for (g in groups) {
    sel <- grp %in% g
    lg_s <- lg[sel]; ls_s <- ls[sel]
    denom <- sum(lg_s^2) + sum(ls_s^2)
    sub_common <- sum(lg_s)^2 + sum(ls_s)^2
    sub_total <- block_var(items[sel])
    rows[[length(rows) + 1L]] <- data.frame(
      factor = g,
      ECV = if (denom > 0) sum(lg_s^2) / denom else NA_real_,
      omega = if (sub_total > 0) sub_common / sub_total else NA_real_,
      omega_H = NA_real_,
      omega_HS = if (sub_total > 0) sum(ls_s)^2 / sub_total else NA_real_,
      PUC = NA_real_, H = H_of(ls_s),
      stringsAsFactors = FALSE
    )
  }
  out <- do.call(rbind, rows)
  row.names(out) <- NULL

  # The maximal reliability the fitted model implies for each factor: the
  # squared multiple correlation of the factor with all of the items,
  # l_f' R^{-1} l_f, where l_f holds the factor's standardized loadings
  # (zero for an item off it) and R is the item correlation matrix the
  # model reproduces. With orthogonal factors of unit variance the
  # covariance of the items with f is l_f, so this is the largest
  # reliability any weighted sum of the items reaches as a measure of f
  # (Hancock and Mueller's H with every other source of item variance,
  # the other factors' included, in the error covariance matrix, which is
  # not diagonal), and the square of the determinacy of f's regression
  # factor score. It does not depend on `denominator`, and it is not
  # defined for an improper solution.
  r_imp <- lavaan::lavInspect(fit, "cor.ov")
  R_imp <- stats::cov2cor(unclass(r_imp)[items, items, drop = FALSE])
  L <- matrix(0, length(items), length(factors),
              dimnames = list(items, c(general, groups)))
  L[, general] <- lg
  on_group <- !is.na(grp)
  L[cbind(items[on_group], grp[on_group])] <- ls[on_group]
  maximal <- stats::setNames(rep(NA_real_, ncol(L)), colnames(L))
  if (!improper) {
    RinvL <- tryCatch(solve(R_imp, L), error = function(e) NULL)
    if (!is.null(RinvL)) maximal[] <- colSums(L * RinvL)
  }

  res <- .as_dmar_tbl(out)
  attr(res, "improper") <- improper
  attr(res, "maximal_reliability") <- maximal
  attr(res, "denominator") <- denominator
  attr(res, "metric") <- if (isTRUE(lavaan::lavInspect(fit, "categorical"))) {
    "latent_response"
  } else {
    "observed"
  }
  res
}
