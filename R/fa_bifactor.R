#' Bifactor Confirmatory Factor Analysis Model
#'
#' Fits the bifactor model, in which a general factor accounts for what
#' every item of a scale has in common and two or more group factors
#' account for what the items of each subscale share beyond it, with all
#' of the factors uncorrelated (Holzinger & Swineford, 1937; Reise, 2012).
#' The model is described as in \code{\link{fa_confirmatory}}, by naming
#' the items of each group factor, and is fit by maximum likelihood to raw
#' data (by full information maximum likelihood, the default, or after
#' listwise deletion) or to a covariance matrix. The table reports the
#' loadings, error variances, and intercepts with their confidence
#' intervals; the indices that describe how much of the total score and of
#' each subscale score the general factor accounts for, each with a delta
#' method standard error and a confidence interval (Rodriguez, Reise, &
#' Haviland, 2016): omega and omega hierarchical (omega_H) for the total
#' score, omega and omega hierarchical subscale (omega_HS) for each
#' subscale score, the explained common variance (ECV), coefficient
#' \emph{H}, and the maximal reliability of each factor, and the percentage
#' of uncontaminated correlations (PUC); and the fit of the model.
#'
#' @inheritParams fa_confirmatory
#' @param factors Named list of the group factors. Each element names a
#'   group factor and gives the character vector of its items; each item
#'   is listed under one group factor, and the general factor loads on all
#'   of them. Two or more group factors are needed, each with three or more
#'   items (two when their loadings are equated with
#'   \code{equal_loading}); see \emph{Details} for the structures that are
#'   refused. An item, here or in \code{general_only}, that is not a column
#'   of \code{data} or an item of \code{S} is refused with an error of
#'   class \code{dmar_items_invalid} (see \code{\link{dmar_metadata}}).
#' @param general The name of the general factor, \code{"general"} by
#'   default. It names the general factor's rows (\code{lambda_general_1},
#'   \code{omega_H_general}).
#' @param general_only \code{NULL} (the default), or a character vector of
#'   items on the general factor alone, with no group factor: an item
#'   written for the construct as a whole, say, or the items of a subscale
#'   too small to carry a group factor of its own. They count in the
#'   general factor's rows (its loadings, omega, omega_H, ECV, \emph{H},
#'   maximal reliability, and PUC) and in no subscale.
#' @param equal_loading Logical, or a named logical vector with one element
#'   for the general factor and one for each group factor. \code{TRUE}
#'   constrains the loadings within a factor to a single value. Defaults to
#'   \code{FALSE}.
#' @param equal_intercept Logical, or a named logical vector with one
#'   element per group factor. \code{TRUE} constrains the intercepts of a
#'   group factor's items to a single value; this requires the mean
#'   structure. The items on the general factor alone keep their own.
#'   Defaults to \code{FALSE}.
#' @param equal_error Logical, or a named logical vector with one element
#'   per group factor. \code{TRUE} constrains the error variances of a
#'   group factor's items to a single value; the items on the general
#'   factor alone keep their own. Defaults to \code{FALSE}.
#' @param meanstructure Logical or \code{NULL}. \code{NULL} (default)
#'   models the mean structure exactly when it is needed: under FIML (the
#'   default with raw data, so the item intercepts, the \code{nu} rows,
#'   join the table), when any \code{equal_intercept} is \code{TRUE}, or
#'   when \code{M} is supplied. \code{TRUE} models the intercepts
#'   regardless (raw data or \code{M} required); \code{FALSE} suppresses
#'   them, which FIML does not allow.
#' @param estimator \code{"ML"} (default), maximum likelihood, or
#'   \code{"MLR"}, maximum likelihood with standard errors and a test
#'   statistic robust to nonnormality (raw data required). Ordered items
#'   are not available; see \emph{Details}.
#' @param missing The missing data treatment for raw data: \code{"fiml"}
#'   (default; full information maximum likelihood, which uses every case
#'   with an observed item; \pkg{lavaan}'s name \code{"ml"} is accepted) or
#'   \code{"listwise"}. The default is resolved as in
#'   \code{\link{fa_confirmatory}}: with \code{se = "robust.sem"} or a test
#'   statistic FIML cannot compute, the items are fit as they stand when
#'   no item value is missing, and otherwise the call stops and names the
#'   option to choose. Ignored with \code{S}.
#' @param se Standard error method passed to \pkg{lavaan}. The default,
#'   \code{NULL}, follows the estimator: \code{"standard"} for ML and
#'   \code{"robust.huber.white"} for MLR. \code{"robust.sem"} requires
#'   complete data and is refused with FIML, and \code{"none"} requests
#'   point estimates only.
#' @param conf_level Confidence level of the intervals for the parameters
#'   and the indices. Defaults to 0.95. It does not reach the RMSEA
#'   interval, whose level is \code{rmsea_conf_level}.
#' @param rmsea_conf_level Confidence level of the RMSEA interval of the
#'   fit rows, which the \code{rmsea_ci_level} row records. Defaults to
#'   0.90, the conventional level, whose limits decide the tests of close
#'   and not close fit at .05 (see \code{\link{sem_fit_measures}});
#'   \code{rmsea_conf_level = 0.95} gives a 95 percent interval.
#'   It sets the level in the table DMAR returns and in its
#'   \code{glance()}; \code{output = "summary"} returns the \pkg{lavaan}
#'   summary, which shows lavaan's 90 percent interval whatever this
#'   argument is.
#' @param output Format of the returned object:
#'   \describe{
#'     \item{\code{"verbose"}}{(default) The loadings, factor variances
#'       (fixed at 1), error variances, and intercepts, then the indices,
#'       then the fit rows.}
#'     \item{\code{"measurement"}}{The indices alone.}
#'     \item{\code{"summary"}}{The \code{summary()} of the \pkg{lavaan}
#'       fit (not a data frame).}
#'     \item{\code{"standardized"}}{\code{lavaan::standardizedSolution()}
#'       of the fit; see \emph{Details}. Under \code{se = "none"} its
#'       \code{se}, \code{z}, \code{pvalue}, \code{ci.lower}, and
#'       \code{ci.upper} columns are \code{NA}.}
#'     \item{\code{"fit"}}{The \pkg{lavaan} fit, which
#'       \code{\link{bifactor_indices}}, \code{\link{ci_bifactor_indices}},
#'       \code{\link{factor_scores}}, and \code{lavaan::lavTestLRT()}
#'       take.}
#'   }
#' @param \dots For \code{fa_bifactor()}, additional arguments forwarded
#'   to \code{\link[lavaan]{lavaan}} (for example \code{cluster}); a
#'   \code{group} argument is an error, as in \code{fa_confirmatory()},
#'   and so is \code{ordered}. Unused by \code{tidy()} and
#'   \code{glance()}.
#' @param wald_df The reference distribution of the Wald intervals and
#'   \emph{p}-values: \code{Inf} (the default), the normal distribution, or
#'   a single positive number, the degrees of freedom of a \emph{t}
#'   distribution for every parameter and index row (for the indices, on
#'   the logit scale), as in \code{\link{fa_confirmatory}}. The keyword
#'   \code{"missing_information"} is refused: the calibration study behind
#'   it simulated one factor models, and nothing calibrates it for the
#'   loadings of a bifactor model or for its indices. It follows
#'   \code{\dots}, so it must be named in full.
#' @param x A \code{dmar_fa_bifactor} table returned by
#'   \code{fa_bifactor()}.
#'
#' @details
#' \strong{The model.} An item \eqn{i} of group factor \eqn{s} is
#' \deqn{x_i = \nu_i + \lambda_{gi} g + \lambda_{si} s + e_i,}{x_i = nu_i +
#'   lambda_gi g + lambda_si s + e_i,}
#' and an item on the general factor alone lacks the group term. Every
#' factor's variance is fixed at 1 (and its mean at 0 when the mean
#' structure is modeled), every pair of factors is uncorrelated, and the
#' first loading of each factor is freed, so each loading is the
#' covariance of its item with the factor. The general factor carries
#' what all of the items share; a group factor carries only what its
#' subscale's items share beyond that, so its loadings are smaller than
#' those of a factor of the subscale alone, and they can be near zero or
#' change sign when the general factor takes up what the subscale
#' measures. Every loading, error variance, and intercept carries a label,
#' as in \code{\link{fa_confirmatory}}, and the indices are \pkg{lavaan}
#' defined parameters of those labels, so the fit itself gives each one a
#' delta method standard error. The standardized solution is derived from
#' this fit of the covariances (or of the raw data), with the factor
#' variances fixed at 1, not from a fit to correlations.
#'
#' \strong{Scope and identification.} A structure outside the model, or
#' one the data cannot identify, is refused before anything is fit, with
#' an error of class \code{dmar_bifactor_unidentified} (which inherits
#' \code{dmar_refusal}; see \code{\link{dmar_metadata}}) that names the
#' factors or items at issue and gives the kind of refusal as its
#' \code{reason} element. Fewer than two group factors, and an item listed
#' under two group factors, are outside the model the function fits
#' (\code{reason} \code{"scope"}): a general factor beside a single group
#' factor, with items on the general factor alone, can be identified, but
#' omega_HS, ECV, and PUC describe how the common variance divides among
#' subscales, which takes two or more. Names placed twice or that coincide
#' are refused with \code{reason} \code{"names"}: an item listed twice
#' under one group factor, or both under a group factor and in
#' \code{general_only}; the general factor and a group factor of one name,
#' a factor and an item of one name, or labels made from the names that
#' coincide. A group factor of one item, or of two whose loadings are free
#' (they enter the model only through their product), and two group
#' factors whose general loadings are fixed only up to a scale are not
#' identified (\code{reason} \code{"unidentified"}). With two group
#' factors and no item on the general factor alone, the covariances
#' between the items of the two groups fix the general loadings of one
#' group only up to a factor \eqn{c} and those of the other up to
#' \eqn{1/c}, and \eqn{c} must come from inside a group, whose item
#' variances and covariances must outnumber the group loadings and error
#' variances it estimates. Two
#' groups each of three items (or of two items with equal loadings)
#' leave \eqn{c} undetermined, and so the model is not identified; an
#' item on the general factor alone, a group of four or more items, equal
#' loadings or equal error variances within a group of three, or equal
#' general loadings each identify it. Every structure with exactly two
#' group factors that is fit gives a warning of class
#' \code{dmar_bifactor_empirical}: such a model can be empirically
#' underidentified, since when the item covariances are nearly equal
#' within each group and between the groups, different general and group
#' loadings reproduce them equally well (Green & Yang, 2018). Items on the
#' general factor alone, or items whose loadings differ in strength,
#' protect against it. The identification checks run once the other
#' arguments and the data have passed theirs, so the warning comes only
#' with a fit. A count of parameters cannot see empirical
#' underidentification, and when \pkg{lavaan}'s information matrix cannot
#' be inverted, or is singular (scaled to a unit diagonal, its smallest
#' eigenvalue is below \eqn{10^{-6}}), every standard error, test, and
#' interval is \code{NA}, the maximal reliability's included, with a
#' warning. The singularity is checked, as in
#' \code{\link{fa_confirmatory}}, for the observed information of a FIML
#' fit and for the expected information of a listwise fit or of a
#' covariance matrix alike, since \pkg{lavaan} can invert either when it is
#' singular.
#'
#' \strong{The indices.} Write \eqn{\lambda_{gi}}{lambda_gi} for the
#' general loading of item \eqn{i}, \eqn{\lambda_{si}}{lambda_si} for its
#' loading on its group factor \eqn{s} (zero for an item on the general
#' factor alone), \eqn{\psi_i}{psi_i} for its error variance, and
#' \eqn{v_i = \lambda_{gi}^2 + \lambda_{si}^2 + \psi_i}{v_i = lambda_gi^2 +
#' lambda_si^2 + psi_i} for its model implied variance. The variance of the
#' total score, the unweighted sum of the items, is
#' \deqn{\sigma^2_X = (\textstyle\sum_i \lambda_{gi})^2 + \sum_s
#'   (\sum_{i \in s} \lambda_{si})^2 + \sum_i \psi_i,}{sigma^2_X =
#'   (sum_i lambda_gi)^2 + sum_s (sum_{i in s} lambda_si)^2 + sum_i psi_i,}
#' and the rows are
#' \describe{
#'   \item{\code{omega_<general>}}{Omega for the total score,
#'     \eqn{[(\sum_i \lambda_{gi})^2 + \sum_s (\sum_{i \in s}
#'     \lambda_{si})^2] / \sigma^2_X}{[(sum_i lambda_gi)^2 + sum_s (sum_{i
#'     in s} lambda_si)^2] / sigma^2_X}: its reliability with every factor
#'     counted as true score (McDonald, 1999).}
#'   \item{\code{omega_H_<general>}}{Omega hierarchical,
#'     \eqn{(\sum_i \lambda_{gi})^2 / \sigma^2_X}{(sum_i lambda_gi)^2 /
#'     sigma^2_X}: the share of the total score's variance that is due to
#'     the general factor alone, the group factors' variance counted as
#'     error. Its square root is the correlation of the total score with
#'     the general factor.}
#'   \item{\code{omega_<s>}}{Omega for the subscale score of group factor
#'     \eqn{s}, the sum of its items, \eqn{(\sum_{i \in s} \lambda_{gi})^2 +
#'     (\sum_{i \in s} \lambda_{si})^2}{(sum_{i in s} lambda_gi)^2 +
#'     (sum_{i in s} lambda_si)^2} over that score's variance, which adds
#'     \eqn{\sum_{i \in s} \psi_i}{sum_{i in s} psi_i}: its reliability
#'     with both factors counted as true score.}
#'   \item{\code{omega_HS_<s>}}{Omega hierarchical subscale,
#'     \eqn{(\sum_{i \in s} \lambda_{si})^2}{(sum_{i in s} lambda_si)^2}
#'     over the same variance: the share of the subscale score's variance
#'     that is due to its group factor alone, what the subscale score says
#'     beyond the total.}
#'   \item{\code{ECV_<general>}, \code{ECV_<s>}}{The explained common
#'     variance, formed from the standardized loadings
#'     \eqn{\lambda_{gi} / \sqrt{v_i}}{lambda_gi / sqrt(v_i)} and
#'     \eqn{\lambda_{si} / \sqrt{v_i}}{lambda_si / sqrt(v_i)}, as
#'     \code{\link{bifactor_indices}} forms it: for the general factor, its
#'     share of the common variance of all of the items,
#'     \eqn{\sum_i \lambda_{gi}^2 / v_i}{sum_i lambda_gi^2 / v_i} over
#'     \eqn{\sum_i (\lambda_{gi}^2 + \lambda_{si}^2) / v_i}{sum_i
#'     (lambda_gi^2 + lambda_si^2) / v_i}; for a group factor's row, the
#'     general factor's share of the common variance of that factor's
#'     items, the same sums over them, so one minus it is the group
#'     factor's own share.}
#'   \item{\code{H_<general>}, \code{H_<s>}}{Coefficient \emph{H} in the
#'     form of Rodriguez, Reise, and Haviland (2016): the closed form of
#'     Hancock and Mueller (2001) over the factor's loadings
#'     \eqn{\lambda_{fi}}{lambda_fi}, with every other source of an item's
#'     variance (the other factors' included) counted as its error,
#'     \eqn{1 / (1 + 1 / \sum_i \lambda_{fi}^2 / (v_i -
#'     \lambda_{fi}^2))}{1 / (1 + 1 / sum_i lambda_fi^2 / (v_i -
#'     lambda_fi^2))}.}
#'   \item{\code{maximal_reliability_<general>},
#'     \code{maximal_reliability_<s>}}{The maximal reliability of the
#'     factor under the fitted model, \eqn{\lambda_f' \Sigma^{-1}
#'     \lambda_f}{lambda_f' Sigma^-1 lambda_f}, with \eqn{\lambda_f}{lambda_f}
#'     the factor's loadings (zero for an item off it) and
#'     \eqn{\Sigma}{Sigma} the model implied covariance matrix of the items:
#'     the squared multiple correlation of the factor with all of the
#'     items, the largest reliability any weighted sum of them reaches as a
#'     measure of that factor, and the square of the determinacy of the
#'     factor's regression score. \emph{H} treats the variance an item
#'     shares with the other factors as uncorrelated error, which in a
#'     bifactor model it is not, so the two differ, in either direction;
#'     \code{\link{bifactor_indices}} says when.}
#'   \item{\code{PUC_<general>}}{The percentage of uncontaminated
#'     correlations, as a proportion: the share of the pairs of items that
#'     do not share a group factor, whose correlations reflect the general
#'     factor alone. It depends on the structure, not on the data, so it
#'     has no standard error or interval.}
#' }
#' Rodriguez, Reise, and Haviland (2016) describe how each is read. In
#' brief, omega_H against omega for the total score says how much of the
#' total score's reliable variance is general; omega_HS against omega for
#' a subscale says how much of the subscale score's reliable variance is
#' its own; ECV and PUC describe how close the items are to measuring one
#' thing; and \emph{H} and the maximal reliability describe how well a
#' factor could be measured by an optimally weighted composite. The table
#' reports the numbers with their uncertainty; whether a total or
#' subscale score is worth reporting is the researcher's judgment, and the
#' package draws no line.
#'
#' \strong{The metric of omega.} Omega, omega_H, and omega_HS here are the
#' reliabilities of unit-weighted sums of the items as they are scored,
#' the coefficients that \code{\link{fa_confirmatory}} and
#' \code{semTools::compRelSEM()} compute.
#' \code{\link{bifactor_indices}} computes them from the standardized
#' solution, for sums of the standardized items, as Rodriguez, Reise, and
#' Haviland (2016) write them. The two agree when the items have equal
#' variances and differ otherwise: on the example below, whose tests'
#' variances range from 8.0 to 153.3, omega_H is 0.552 for the sum of the
#' scores and 0.577 for the sum of the standardized scores, and on the
#' simulated example of \code{\link{bifactor_indices}}, whose items have
#' nearly equal variances, 0.67305 against 0.67348. ECV is formed from the
#' standardized loadings, and \emph{H}, the maximal reliability, and PUC do
#' not depend on the metric of the items, so those rows equal
#' \code{\link{bifactor_indices}}'s with its model implied denominator.
#' With \code{output = "standardized"}, the defined indices are evaluated
#' at the standardized parameters, so the omega rows there are those of
#' the standardized items. \pkg{lavaan} evaluates a defined parameter at the
#' standardized values of the labels it names, and a label that
#' \code{equal_loading} or \code{equal_error} gives several items has one
#' unstandardized value but a standardized value for each item, which
#' differ when the items' variances differ; \pkg{lavaan} uses the first
#' item's. Each index that uses such a label is computed here from every
#' item's own standardized loadings and error variance, as
#' \code{\link{fa_confirmatory}} computes its rows, with a delta method
#' standard error from their derivatives with respect to the free
#' parameters and \pkg{lavaan}'s covariance matrix of the estimates (under
#' \code{se = "bootstrap"}, from \pkg{lavaan}'s bootstrap draws); the
#' other indices are \pkg{lavaan}'s. ECV and \emph{H} then equal the rows
#' of the other tables.
#'
#' \strong{Intervals.} The standard errors of omega, omega_H, omega_HS,
#' ECV, and \emph{H} are \pkg{lavaan}'s delta method standard errors of the
#' defined parameters. The maximal reliability needs a matrix inverse,
#' which \pkg{lavaan}'s defined parameters cannot express, so its standard
#' error is the delta method's with the derivatives taken numerically
#' over the free parameters and \pkg{lavaan}'s covariance matrix of their
#' estimates, as \code{\link{fa_exploratory}} forms those of its
#' determinacies. For every index the limits are formed on the logit
#' scale from that standard error and transformed back (Browne, 1982; the
#' confidence interval transformation principle of Kelley, 2007), as
#' \code{\link{fa_confirmatory}} forms the limits of omega and \emph{H}, so
#' they lie within (0, 1); the \code{se} column is the standard error of
#' the index itself. The sample values of \emph{H} and of the maximal
#' reliability overstate their population values on average, the more so
#' for a weakly determined factor (Aguirre-Urreta, Rönkkö, & McIntosh,
#' 2019); \code{\link{ci_bifactor_indices}} on the fit of
#' \code{output = "fit"} gives bootstrap intervals for the indices of
#' \code{\link{bifactor_indices}} and the bootstrap estimate of that bias.
#'
#' \strong{Improper solutions.} A bifactor model often gives an improper
#' solution, an error variance estimated below zero, and the estimate is
#' one of two kinds. It can be a finite maximum of the likelihood: on
#' \code{holzinger_swineford}, nine tests in three groups of three (visual
#' perception, cubes, and lozenges; paragraph comprehension, sentence
#' completion, and word meaning; addition, counting groups of dots, and
#' straight and curved capitals) give visual perception an error variance
#' of -41, and the likelihood falls on either side of that value. Such a
#' solution is reported. The maximal reliability is not defined for it,
#' and neither is \emph{H} for a factor on which an item's standardized
#' loading reaches one, nor an omega whose items' error variances sum
#' below zero; those rows are \code{NA}, with a warning that names the
#' items (beside the warning of class \code{dmar_heywood_warning} that
#' \code{\link{fa_confirmatory}} gives), and the \code{improper}
#' attribute holds the negative error variances. The other rows are
#' reported as estimated, and an improper solution's indices are better
#' read as a sign that the model asks more of the data than they hold. Or
#' the likelihood can have no maximum at all. With a fourth group factor
#' on the Woody-McCall mixed fundamentals, revised paper form board, and
#' flags tests beside those three groups, the residual covariances of the
#' fourth group's tests, beyond the general factor, are all negative,
#' which no group factor can reproduce (a single factor implies three
#' covariances whose product is not negative), and the likelihood keeps
#' rising, toward a limit, as the error variance of the paper form board
#' falls without bound. \pkg{lavaan}'s optimizer stops at -1163 under
#' FIML, where its steps stop gaining, a point set by its starting values
#' and tolerance rather than by the data, and every index read there
#' would describe the optimizer. So when an error variance is negative,
#' the function fits the model again with it held further out (at ten
#' times its estimate, or one model implied variance below it, whichever
#' is further), and when the log likelihood there is higher, or lower by
#' no more than \eqn{10^{-6}} per case (a likelihood that flat leaves the
#' estimate as much to where the optimizer stopped as a ridge does), it
#' stops, for every \code{output}
#' (\code{"fit"} included), with an error of class \code{dmar_no_maximum}
#' (see \code{\link{dmar_metadata}}) that names the test and prints the
#' \pkg{lavaan} code of the fit; there, \code{bounds = "pos.var"} (which
#' holds every error variance at zero or above, and which \code{\dots}
#' passes on from here as well) or a different model, such as one with
#' the covariances of the group's items freed in place of the group
#' factor, are the choices.
#'
#' \strong{Ordered items.} \code{fa_bifactor()} fits continuous items, and
#' \code{ordered} is refused. Fit the model to ordered items with
#' \code{lavaan::cfa()} (\code{ordered = TRUE}, \code{orthogonal = TRUE},
#' \code{std.lv = TRUE}) and pass the fit to \code{\link{bifactor_indices}},
#' which reports the indices on the metric of the latent responses behind
#' the categories; the bifactor coefficients of a sum of the observed
#' categories are not implemented.
#'
#' \strong{The other bifactor functions.} \code{output = "fit"} returns
#' the \pkg{lavaan} fit, so the same model feeds
#' \code{\link{bifactor_indices}}, which reports the indices on the
#' standardized metric, without intervals, and can take the observed
#' composite variance as the denominator, and
#' \code{\link{ci_bifactor_indices}}, which gives them bootstrap
#' intervals, an alternative to the delta method intervals here that does
#' not rest on the normal approximation and that sets improper
#' replications aside.
#'
#' Row names under \code{term}: \code{lambda_<factor>_<j>} (loadings, the
#' general factor's numbered over every item, the items of the group
#' factors in their order and then those of \code{general_only};
#' \code{lambda_<factor>} when equated), \code{phi_<factor>} (the factor
#' variances, fixed at 1), \code{psi_<factor>_<j>} and
#' \code{nu_<factor>_<j>} (the error variances and intercepts of a group
#' factor's items; \code{psi_<general>_<j>} and \code{nu_<general>_<j>}
#' for an item on the general factor alone, numbered as its general
#' loading), the indices above, and \code{\link{fa_confirmatory}}'s fit
#' rows \code{chi_square} through \code{H1}, whose versions and the lines
#' under the table are as there.
#'
#' @return
#' For \code{output = "verbose"} (default) and \code{"measurement"}, a
#' \code{data.frame} (classes \code{dmar_fa_bifactor}, \code{dmar_tbl})
#' with columns \code{syntax}, \code{term}, \code{estimate}, \code{se},
#' \code{test_statistic}, \code{p_value}, \code{ci_lower}, and \code{ci_upper}.
#' Its attributes record the model: \code{general}, \code{factors}, and
#' \code{general_only}; \code{equal}, a list naming the factors whose
#' loadings, error variances, and intercepts are equated;
#' \code{has_means}; \code{omega_metric} (\code{"raw_item_sum"});
#' \code{improper}, the negative error variances named by item (empty for
#' a proper solution); and \code{wald_df}. The verbose table also carries
#' \code{fa_confirmatory()}'s \code{test}, \code{estimator},
#' \code{missing}, \code{fit_indices}, and \code{fit_rows}. The printed
#' header describes the model, and the lines under the table say what the
#' indices describe. For \code{"summary"}, the \pkg{lavaan} summary; for
#' \code{"standardized"}, the standardized solution with the
#' \code{"wald_df"} attribute; for \code{"fit"}, the \pkg{lavaan} fit.
#' \code{tidy()} gives the parameter and index rows and \code{glance()}
#' one row of the fit rows and what produced them, as for
#' \code{\link{fa_confirmatory}}.
#'
#' @references
#' Aguirre-Urreta, M. I., Rönkkö, M., & McIntosh, C. N. (2019). A
#'   cautionary note on the finite sample behavior of maximal reliability.
#'   \emph{Psychological Methods, 24}(2), 236--252.
#'   \doi{10.1037/met0000176}
#'
#' Browne, M. W. (1982). Covariance structures. In D. M. Hawkins (Ed.),
#'   \emph{Topics in applied multivariate analysis} (pp. 72--141).
#'   Cambridge, UK: Cambridge University Press.
#'
#' Green, S., & Yang, Y. (2018). Empirical underidentification with the
#'   bifactor model: A case study. \emph{Educational and Psychological
#'   Measurement, 78}(5), 717--736. \doi{10.1177/0013164417719947}
#'
#' Hancock, G. R., & Mueller, R. O. (2001). Rethinking construct
#'   reliability within latent variable systems. In R. Cudeck, S. du Toit,
#'   & D. Sörbom (Eds.), \emph{Structural equation modeling: Present and
#'   future} (pp. 195--216). Scientific Software International.
#'
#' Holzinger, K. J., & Swineford, F. (1937). The bi-factor method.
#'   \emph{Psychometrika, 2}(1), 41--54. \doi{10.1007/BF02287965}
#'
#' Kelley, K. (2007). Confidence intervals for standardized effect sizes:
#'   Theory, application, and implementation. \emph{Journal of Statistical
#'   Software, 20}(8), 1--24. \doi{10.18637/jss.v020.i08}
#'
#' McDonald, R. P. (1999). \emph{Test theory: A unified treatment}.
#'   Erlbaum.
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
#' @seealso \code{\link{fa_confirmatory}} for the correlated factors
#'   model of the same items; \code{\link{fa_second_order}} for the
#'   second-order model whose first-order factors are the group factors,
#'   which is nested in the bifactor model; \code{\link{bifactor_indices}} and
#'   \code{\link{ci_bifactor_indices}} for the indices of a fitted
#'   bifactor model on the standardized metric and their bootstrap
#'   intervals; \code{\link{factor_scores}} for the factor scores and their
#'   determinacy; \code{\link{reliability_omega}} and
#'   \code{\link{reliability_H}}.
#'
#' @family multivariate and latent variable methods
#'
#' @keywords multivariate
#'
#' @examples
#' data(holzinger_swineford)
#' hs_groups <- list(
#'   spatial = c("t1_visual_perception", "t2_cubes", "t3_paper_form_board",
#'               "t4_lozenges"),
#'   verbal  = c("t5_general_information", "t6_paragraph_comprehension",
#'               "t7_sentence", "t8_word_classification", "t9_word_meaning"),
#'   memory  = c("t14_word_recognition", "t15_number_recognition",
#'               "t16_figure_recognition", "t17_object_number",
#'               "t18_number_figure", "t19_figure_word"))
#'
#' # A general factor on the 15 tests and three group factors: the
#' # indices, each with a delta method interval.
#' fa_bifactor(holzinger_swineford, hs_groups, output = "measurement")
#'
#' # The fit statistics of the same model, one row.
#' generics::glance(fa_bifactor(holzinger_swineford, hs_groups))
#'
#' # The lavaan fit feeds bifactor_indices(), whose omega, omega_H, and
#' # omega_HS describe sums of the standardized tests (0.577 for omega_H,
#' # against 0.552 for the sum of the scores above).
#' fit <- fa_bifactor(holzinger_swineford, hs_groups, output = "fit")
#' bifactor_indices(fit)
#'
#' @export
fa_bifactor <- function(data = NULL, factors, S = NULL, N = NULL, M = NULL,
                        general = "general", general_only = NULL,
                        equal_loading = FALSE, equal_intercept = FALSE,
                        equal_error = FALSE, meanstructure = NULL,
                        estimator = "ML", missing = "fiml", se = NULL,
                        conf_level = 0.95, rmsea_conf_level = 0.90,
                        output = c("verbose", "measurement", "summary",
                                   "standardized", "fit"),
                        ...,
                        wald_df = Inf) {
  output <- match.arg(output)
  # Whether 'missing' was left at its default, read before anything else
  # (see fa_confirmatory()).
  missing_default <- base::missing(missing)
  # Ordered items would reach lavaan through '...' (as would a name that
  # partially matches its 'ordered' formal), and lavaan would fit them; the
  # indices here describe sums of continuous items, so they are refused by
  # name.
  dot_names <- names(list(...))
  if (!is.null(dot_names) &&
      any(nzchar(dot_names) & startsWith("ordered", dot_names))) {
    .dmar_refuse_argument(c(
      "fa_bifactor() fits continuous items by maximum likelihood, and ",
      "ordered items are not available: its omega, omega_H, and ",
      "omega_HS describe sums of the items as scored, and the bifactor ",
      "coefficients of a sum of ordered categories are not implemented. ",
      "bifactor_indices() gives the indices on the metric of the latent ",
      "responses behind the categories, from a lavaan fit to the ordered ",
      "items (lavaan::cfa(..., ordered = TRUE, orthogonal = TRUE, ",
      "std.lv = TRUE))."), arg = NULL)
  }
  .fa_cfa_engine(..., data = data, factors = factors, S = S, N = N, M = M,
                 equal_loading = equal_loading,
                 equal_intercept = equal_intercept,
                 equal_error = equal_error, correlated_factors = FALSE,
                 meanstructure = meanstructure, estimator = estimator,
                 missing = missing, ordered = NULL, se = se,
                 conf_level = conf_level, rmsea_conf_level = rmsea_conf_level,
                 output = output, wald_df = wald_df,
                 missing_default = missing_default, structure = "bifactor",
                 general = general, general_only = general_only)
}

# The model introduces the table, at the console and in a knitted report
# alike (the "header" block of .dmar_tbl_footer(), R/dmar_output_helpers.R),
# and the lines under it say what the omega rows describe, what a group
# factor's rows mean, and which rows an improper solution left NA; the
# lines that name the version of the fit rows follow, as for
# fa_confirmatory(). A column subset of the table, which loses the
# attributes, prints the default lines alone.
#' @exportS3Method .dmar_tbl_footer dmar_fa_bifactor
#' @noRd
.dmar_tbl_footer.dmar_fa_bifactor <- function(x, ...) {
  general <- attr(x, "general", exact = TRUE)
  factors <- attr(x, "factors", exact = TRUE)
  if (is.null(general) || !is.list(factors)) return(NextMethod())
  alone <- attr(x, "general_only", exact = TRUE)
  equal <- attr(x, "equal", exact = TRUE)
  improper <- attr(x, "improper", exact = TRUE)
  J <- lengths(factors)
  header <- c(
    sprintf(paste0("Bifactor model: the general factor (%s) loads on all %d ",
                   "items; group factors %s"),
            general, sum(J) + length(alone),
            paste0(names(factors), " (", J, " items)", collapse = ", ")),
    if (length(alone)) {
      sprintf("On the general factor alone: %s", paste(alone, collapse = ", "))
    },
    paste0("Every factor's variance is fixed at 1 and every pair of factors ",
           "is uncorrelated",
           if (length(unlist(equal))) {
             parts <- c(
               if (length(equal$loading)) {
                 paste0("loadings within ", paste(equal$loading, collapse = ", "))
               },
               if (length(equal$error)) {
                 paste0("error variances within ",
                        paste(equal$error, collapse = ", "))
               },
               if (length(equal$intercept)) {
                 paste0("intercepts within ",
                        paste(equal$intercept, collapse = ", "))
               })
             paste0("; equal ", paste(parts, collapse = "; equal "))
           }))
  metric <- paste0(
    "omega, omega_H, and omega_HS describe unit-weighted sums of the items ",
    "as scored (bifactor_indices() gives them for sums of the standardized ",
    "items); ECV is formed from the standardized loadings, and H and ",
    "maximal_reliability do not depend on the metric of the items; PUC ",
    "depends on the structure alone, so it has no interval")
  definitions <- paste0(
    "For a group factor, omega is the reliability of the total score of ",
    "its items, omega_HS the share of that score's variance due to the ",
    "group factor alone, and ECV the general factor's share of the items' ",
    "common variance; H is the index of Rodriguez, Reise, and Haviland ",
    "(2016), which counts the other factors' common variance as ",
    "uncorrelated error, and maximal_reliability is the reliability the ",
    "optimally weighted items reach as a measure of the factor")
  improper_line <- if (length(improper)) {
    paste0("Improper solution (a negative error variance for ",
           paste(names(improper), collapse = ", "), "): the rows the ",
           "warning named are NA")
  }
  c(.dmar_footer_add(character(0), header, "header"), NextMethod(),
    .dmar_footer_add(character(0), metric, "bifactor_metric"),
    .dmar_footer_add(character(0), definitions, "bifactor_definitions"),
    .dmar_footer_add(character(0), improper_line, "bifactor_improper"),
    .dmar_footer_add(character(0), attr(x, "fit_rows", exact = TRUE),
                     "fit_rows"))
}

#' @rdname fa_bifactor
#' @importFrom generics tidy
#' @exportS3Method generics::tidy
tidy.dmar_fa_bifactor <- function(x, ...) .fa_cfa_tidy(x)

#' @rdname fa_bifactor
#' @importFrom generics glance
#' @exportS3Method generics::glance
glance.dmar_fa_bifactor <- function(x, ...) .fa_cfa_glance(x)
