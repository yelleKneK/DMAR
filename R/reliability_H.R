#' Coefficient \emph{H} (Maximal Reliability) With a Confidence Interval
#'
#' @description
#' Estimates coefficient \emph{H} (Hancock & Mueller, 2001), the
#' reliability of the optimally weighted composite of the items that
#' measure one factor, and returns a confidence interval for the
#' population coefficient. \emph{H} is the largest reliability that any
#' weighted sum of the items can reach under the factor model, whereas
#' McDonald's \eqn{\omega} (\code{\link{reliability_omega}}) is the
#' reliability of the unit-weighted sum that people usually score. The
#' coefficient can be computed from a fitted \pkg{lavaan} model
#' (\code{fit}, the recommended input), from raw item scores
#' (\code{data}), from a covariance matrix and its sample size (\code{S}
#' and \code{N}), or from the standardized loadings a paper reports
#' (\code{loadings}). The interval comes from the delta method with the
#' full covariance matrix of the parameter estimates, built on the logit
#' scale by default, or from a nonparametric bootstrap when one is
#' requested.
#'
#' @details
#' \strong{The coefficient.} For a factor \eqn{f} measured by items
#' \eqn{b} that load on it alone, with unstandardized loadings
#' \eqn{\lambda}, residual covariance matrix \eqn{\Theta_{bb}}, and model
#' implied factor variance \eqn{\phi_f},
#' \deqn{q = \phi_f \, \lambda' \Theta_{bb}^{-1} \lambda, \qquad
#'   H = \frac{q}{1 + q}.}
#' With uncorrelated residuals and standardized loadings \eqn{\ell_j} this
#' is the form Hancock and Mueller (2001) gave,
#' \deqn{H = \frac{\sum_j \ell_j^2 / (1 - \ell_j^2)}
#'                {1 + \sum_j \ell_j^2 / (1 - \ell_j^2)},}
#' in which each item contributes its ratio of common to unique variance,
#' \eqn{\ell_j^2 / (1 - \ell_j^2)}. The weights that reach \emph{H} are
#' proportional to \eqn{\Theta_{bb}^{-1} \lambda}, so an item counts in
#' proportion to its loading divided by its residual variance, and adding
#' an item can never lower \emph{H}. Gabler and Raykov (2017, Eq. 9)
#' showed that the same expression with a nondiagonal \eqn{\Theta_{bb}}
#' is the maximal reliability when residuals inside the block are
#' correlated; a fitted model with residual covariances among the
#' factor's items gets that coefficient here, and the printed table says
#' how many entered the weights. Rescaling the factor rescales
#' \eqn{\lambda} and \eqn{\phi_f} inversely, so \emph{H} does not depend
#' on whether the factor was identified by a marker item or by fixing its
#' variance.
#'
#' Two further readings of the number are useful. \emph{H} is the squared
#' correlation between the factor and its best linear predictor from the
#' block's items, the regression factor score computed from those items
#' alone, so it is the reliability of that score: for a one-factor model
#' the determinacy that \code{\link{factor_scores}} reports is the square
#' root of \emph{H}. And because \emph{H} is the largest reliability any
#' weighting of the items can reach under the fitted model, it is at
#' least the model implied \eqn{\omega} of the unit-weighted sum from the
#' same fit (Raykov, 1997), with equality when the weights
#' \eqn{\Theta_{bb}^{-1} \lambda} are all equal, which with uncorrelated
#' residuals means equal standardized loadings. Robust omega, whose
#' denominator is the observed variance of the composite, need not lie
#' below it.
#'
#' \strong{Reading \emph{H} with care.} The sample \emph{H} is a positively
#' biased estimator of the population coefficient, and a composite formed
#' with weights estimated from the sample is on average less reliable
#' than the population \emph{H}. Aguirre-Urreta, Rönkkö, and McIntosh
#' (2019) found both effects pronounced in small samples (around
#' \emph{N} = 100), diminishing as the sample grows and barely noticeable
#' beyond about 1,000, although with three items of modest loadings even
#' 1,000 cases left more than 1\% bias; the reliability estimate of the
#' unit-weighted composite was nearly unbiased in the same conditions.
#' Report \eqn{\omega}, with its interval, for a scale that people will
#' score by summing or averaging its items, and \emph{H} as the ceiling
#' that optimal weights could reach, or as the reliability of scores
#' computed from the factor model. A bootstrap call estimates the bias of
#' the sample \emph{H} (the \code{bias} element of the \code{bootstrap}
#' attribute), and the printed table shows it.
#'
#' \strong{The four inputs.} Exactly one of \code{fit}, \code{data},
#' \code{S}, and \code{loadings} is supplied.
#' \describe{
#'   \item{\code{fit}}{A fitted \pkg{lavaan} model, from
#'   \code{lavaan::cfa()}, \code{lavaan::sem()}, or
#'   \code{fa_confirmatory(output = "fit")}; the recommended input, because it
#'   carries everything the coefficient and its interval need, and the
#'   model can be the one the analysis reports. Any factor measured by
#'   observed items can be named with \code{factor} (a factor of a
#'   confirmatory model, or an exogenous or endogenous factor of a
#'   structural model), and any group of a multiple group model with
#'   \code{group}. \emph{H} is read from the fitted matrices, and its
#'   standard error comes from the fit's own covariance matrix of the
#'   estimates, so it is the standard, robust, or FIML standard error the
#'   model was fitted with. With ordered items fitted by WLSMV (either
#'   parameterization) \emph{H} is the maximal reliability of the latent
#'   responses behind the categories, not of the observed categories; the
#'   \code{metric} attribute and the printed table say so. A loading or a
#'   regression the model holds at zero is not one, whether it is fixed
#'   (\code{0*item}) or constrained (\code{a*item} with \code{a == 0},
#'   directly or through a chain of equality constraints between labels,
#'   as in \code{a == b} with \code{b == 0}): the model
#'   written with that term is the model written without it and gets the
#'   same \emph{H}, so an item listed under the factor with its loading
#'   held at zero is not one of the factor's items, and a cross-loading
#'   held at zero refuses nothing. A fit is
#'   refused, with a message naming the remedy, when an item of the
#'   factor also loads on another factor or is regressed on another
#'   variable (for the factors of a bifactor model,
#'   \code{\link{bifactor_indices}} reports the index that Rodriguez,
#'   Reise, and Haviland, 2016, call \emph{H}, which counts the other
#'   factors' common variance as uncorrelated error and is not this
#'   maximal reliability, and, in its \code{"maximal_reliability"}
#'   attribute, the maximal reliability the fitted model implies for each
#'   factor, and \code{\link{ci_bifactor_indices}} gives both a
#'   confidence interval), when the factor is measured by other factors,
#'   when an observed covariate predicts the factor, when the solution is
#'   improper (an item's residual variance that is not positive, for
#'   which \emph{H} is not defined), and when a factor or a group must be
#'   named and is not.}
#'   \item{\code{data}}{Raw scores on the items of one scale. The
#'   single-factor congeneric model is fitted exactly as
#'   \code{\link{reliability_omega}} fits it, so \emph{H} and
#'   \eqn{\omega} rest on the same fit, the same cases, and the same kind
#'   of standard error, and \code{missing} and \code{aux} work as they do
#'   there. Three or more items are required, since the model is not
#'   identified with two.}
#'   \item{\code{S}}{A covariance matrix among the items, with its sample
#'   size \code{N}; the same single-factor model, fitted by maximum
#'   likelihood.}
#'   \item{\code{loadings}}{The standardized loadings of one factor, as a
#'   paper reports them, with the closed form above (which assumes
#'   uncorrelated residuals). For an interval, supply
#'   \code{vcov_loadings}, the covariance matrix of the standardized
#'   loading estimates (for example the matching block of
#'   \code{lavaan::lavInspect(fit, "vcov.std.all")}); \code{se_loadings},
#'   their marginal standard errors, is accepted, but it ignores the
#'   covariances and warns that it does. Loadings alone carry nothing an
#'   interval can be formed from, so without either the estimate is
#'   returned with \code{se} and the limits \code{NA}, and the printed
#'   table says what an interval would need.}
#' }
#'
#' \strong{The interval.} The delta method variance of \eqn{\hat H} is
#' \eqn{g' V g}, where \eqn{V} is the covariance matrix of all the free
#' parameter estimates of the fit and \eqn{g} is the gradient of
#' \emph{H} with respect to them, computed by central differences of the
#' closed form above. \eqn{V} is used only when it can be, by the check
#' the factor analysis tables (\code{\link{fa_confirmatory}},
#' \code{\link{fa_bifactor}}, \code{\link{fa_second_order}}) make before
#' they report a standard error: when \pkg{lavaan} could not compute it,
#' or computed it from an information matrix, observed or expected, that
#' is singular in the directions the fit's equality constraints leave free
#' (scaled to a unit diagonal, its smallest eigenvalue is below
#' \eqn{10^{-6}}), a parameter of the fit is not identified at the
#' estimates and \eqn{g' V g} means nothing. The estimate is then returned
#' with \code{se} and the limits \code{NA}, and with the warning, of class
#' \code{dmar_result_qualified}, that those tables give for the same fit:
#' the same words for a singular matrix, and for one \pkg{lavaan} could
#' not compute a sentence that ends by naming the standard error and
#' interval of \emph{H} rather than the tables' columns. The check is of
#' the whole fit, as the tables make it, so a factor whose own parameters
#' are identified loses its interval beside one whose parameters are not.
#' An interval
#' from \code{se_loadings} uses only the
#' marginal standard errors of the standardized loadings, the diagonal
#' of \eqn{V}. The loading estimates of a fitted factor model are
#' correlated, and the diagonal can understate or overstate the variance:
#' in seven blocks of three to six tests from
#' \code{\link{holzinger_swineford}}, each fitted as one factor, the
#' standard error from the full matrix ranged from about 0.7 to 1.3 times
#' the one from the diagonal. By default the Wald interval is built on
#' the logit scale and transformed back (Browne, 1982): the limits of a
#' monotone transformation of a parameter are the transformation of its
#' limits, the confidence interval transformation principle that Kelley
#' (2007) sets out, so the interval stays inside (0, 1) and is asymmetric
#' near the upper bound, as the sampling distribution of a high
#' reliability is; \code{"delta"}, \code{"ml"},
#' and \code{"mlr"} give the untransformed Wald interval. The default
#' standard error depends on the input: with \code{data}, the robust
#' maximum likelihood standard error (Satorra-Bentler, or Yuan-Bentler
#' under FIML), as for \eqn{\omega}, because Kelley and Pornprasertmanit
#' (2016) found the normal theory Wald interval wanting when items depart
#' from normality; with \code{S}, the maximum likelihood standard error;
#' and with \code{fit}, the fit's own. With \code{fit}, \code{"ml"} and
#' \code{"mlr"} (and their \code{_logistic} forms) are accepted when they
#' name the standard error the fit already carries, and otherwise stop
#' with a message saying which one it carries.
#'
#' The delta method interval is the fast default. For a reported
#' analysis, the percentile or BCa bootstrap with the model refitted on
#' every replication (by FIML when data are missing) is the safer choice,
#' as it is for \eqn{\omega} (Kelley & Pornprasertmanit, 2016), and it
#' also estimates the bias of the sample \emph{H}. With \code{data} the
#' rows are resampled and the single-factor model refitted; with
#' \code{fit} the rows the model was fitted to are resampled and the
#' model itself refitted through \code{lavaan::lavaan()}, with the fitted
#' estimates as starting values (for a single-group fit that carries its
#' raw data). Before any resample is drawn, the model is refitted this
#' way to the rows of \code{fit} itself, and the call stops unless that
#' refit reproduces the estimate, so the replications are refits of the
#' model the estimate describes. A fit with sampling weights is refused,
#' since a resample of its rows does not carry them and every replication
#' would be an unweighted fit, and so is a fit with a \code{cluster}
#' variable, whose cases are not independent and cannot be resampled one
#' at a time; the delta method interval, built from the fit's own
#' covariance matrix of the estimates, serves both. The bootstrap engine
#' is the one the whole reliability family uses (Efron & Tibshirani,
#' 1993): replications that fail to converge, give an improper solution
#' (or, with three items, a resample whose correlations admit no
#' single-factor solution), or raise an error are set aside and counted
#' (a resample in which an item takes a single value is refused before
#' the refit, by name, as the data would be; \pkg{lavaan} would stop on
#' the item after printing its table of the variables; so is a resample
#' with too few distinct cases, by the rule of \code{missing} with
#' \code{data}, and with \code{fit} by the rule of the fit's own
#' treatment of incomplete rows: under listwise deletion, or with
#' complete data, one with no more distinct rows than the model's
#' continuous variables, whose sample covariance matrix is singular, and
#' under FIML or two-stage estimation one with no more distinct cases
#' that observe a variable, while a fit by pairwise deletion is held to
#' neither),
#' and the counts travel in the \code{bootstrap} attribute, print under
#' the table, and are stated in a warning. A replication whose likelihood
#' has no maximum (an error variance running off below zero while the
#' likelihood is flat or still rising; see \code{\link{fa_bifactor}}) is
#' counted among those that raised an error, and a fit with no maximum
#' behind the estimate itself stops with an error of class
#' \code{dmar_no_maximum}, which with \code{data} prints the \pkg{lavaan}
#' code of the single-factor fit, where \code{bounds = "pos.var"} or a
#' different model are the choices. \emph{H} is not defined for an
#' improper solution, so the attribute records \code{improper_rule =
#' "drop"}. Each replication
#' is a model fit, a few hundredths of a second for a one-factor model on
#' complete data and longer under FIML, so the default \code{B = 10000}
#' (the family's default, chosen for the stability of the BCa quantiles)
#' takes minutes. The examples therefore stop at the delta method; the
#' bootstrap call for their five verbal tests is
#' \code{reliability_H(data = holzinger_swineford[, verbal], ci_method =
#' "bca", B = 2000, seed = 113)}.
#'
#' \strong{Missing data.} With \code{data}, incomplete items are fitted
#' by default by full information maximum likelihood (FIML), which keeps
#' every case with at least one observed item or observed auxiliary
#' variable, with auxiliary variables (\code{aux}) entered as saturated
#' correlates (Graham, 2003), and \code{missing = "listwise"} fits the
#' complete rows alone; data in which no value is missing, and \code{S},
#' are fitted as they stand. A row with every item missing does not count
#' as incomplete data, since without auxiliary variables both treatments
#' drop it. The default analysis is therefore FIML whenever an item value
#' is missing and the complete data analysis otherwise, unlike
#' \code{\link{fa_confirmatory}}, which fits complete data by FIML too.
#' That is the rule of \code{\link{reliability_omega}}, whose
#' complete-data intervals need the complete data analysis, and it is
#' followed here so that the two fit the same single-factor model the
#' same way; on incomplete data it agrees with
#' \code{\link{fa_confirmatory}}. With \code{fit},
#' the fit's own treatment applies: a model fitted with \code{missing =
#' "ml"} gives \emph{H} and its standard error by FIML. Either way
#' \code{N} counts the cases the fit used and \code{N_complete} those
#' complete on the factor's items, so the cost of listwise deletion is
#' visible. Under
#' \code{missing = "fiml"} a pair of items or auxiliary variables that no
#' case observes together is refused by name before any fit, since the
#' data carry no information about its covariance (lavaan only warns, and
#' its EM returns an arbitrary value for it), and so is a pair the cases
#' observe together only where one of them takes a single value, since
#' those cases carry none about how the two vary together; in a
#' bootstrap, a resample that leaves such a pair is set aside with the
#' same message and counted.
#'
#' \strong{From standardized loadings.} With \code{loadings} alone the
#' function returns the estimate, with \code{se} and the limits
#' \code{NA}, because standardized loadings carry nothing an interval
#' can be formed from (the default \code{ci_method} is then
#' \code{"none"}). With \code{vcov_loadings} the interval is the delta
#' method interval from the full covariance matrix of the loading
#' estimates. With \code{se_loadings} it is the delta method interval
#' from their marginal standard errors, and the call warns every time,
#' with class \code{dmar_H_marginal_se}, that the interval ignores the
#' covariances among the loading estimates. Every input returns the
#' reliability family table described under \emph{Value}; the sampling
#' variance of \emph{H} is the square of the row \code{se}.
#' \code{reliability(type = "H")} dispatches here.
#'
#' The coefficient that \code{\link{fa_confirmatory}} reports as
#' \code{H_<factor>} is the same quantity, with \pkg{lavaan}'s delta
#' method standard error and an untransformed interval;
#' \code{reliability_H(fit, factor = "speed", ci_method = "delta")}
#' reproduces its estimate and standard error for the factor
#' \code{speed}.
#'
#' @param loadings Numeric vector of the standardized loadings of one
#'   factor, each in \eqn{(-1, 1)}, at least two. A fitted \pkg{lavaan}
#'   model supplied in this position is taken as \code{fit}, so
#'   \code{reliability_H(fit)} works.
#' @param se_loadings Optional vector of the marginal standard errors of
#'   the standardized loadings, the same length as \code{loadings}. The
#'   interval it gives ignores the covariances among the loading
#'   estimates and warns that it does; prefer \code{vcov_loadings} or
#'   \code{fit}.
#' @param conf_level Confidence level for the interval. Defaults to
#'   \code{0.95}.
#' @param fit A fitted \pkg{lavaan} model (from \code{lavaan::cfa()},
#'   \code{lavaan::sem()}, or \code{fa_confirmatory(output = "fit")}), the
#'   recommended input; see \emph{Details}.
#' @param factor With \code{fit}, the name of the factor whose \emph{H}
#'   is wanted; required when the model has more than one factor measured
#'   by observed items.
#' @param group With a multiple group \code{fit}, the label of the group
#'   (as \code{lavaan::lavInspect(fit, "group.label")} gives it) whose
#'   \emph{H} is wanted; required for such a fit.
#' @param data A numeric matrix or data frame of item scores of one scale
#'   (rows are respondents, columns are items, three or more), and, when
#'   \code{aux} is supplied, the auxiliary columns. An item column that
#'   is not numeric (or logical), and an item that takes a single value
#'   in the rows analyzed, are refused by name before any fit, as is a
#'   zero variance on the diagonal of \code{S}: the single-factor model
#'   cannot be fit with such an item (lavaan would stop on it with an
#'   error that names no variable).
#'   SPSS user-defined missing values, which
#'   \code{haven::read_sav(user_na = TRUE)} keeps as values of a
#'   \code{haven_labelled_spss} column, are missing, as SPSS treats
#'   them.
#' @param S A symmetric covariance matrix among the items (three or
#'   more); \code{N} is then required, and only \code{"ml_logistic"},
#'   \code{"ml"}, and \code{"none"} apply.
#' @param N Sample size. Required with \code{S}; optional with
#'   \code{loadings}, where it is recorded in the table; not used with
#'   \code{fit} or \code{data}, which count their own cases.
#' @param vcov_loadings With \code{loadings}, the \eqn{J \times J}
#'   covariance matrix of the standardized loading estimates, for the
#'   delta method interval with the full covariance matrix.
#' @param ci_method Method for the confidence interval.
#'   \code{"delta_logistic"} and \code{"delta"} use the delta method with
#'   the covariance matrix that comes with the input (the fit's, or
#'   \code{vcov_loadings}); \code{"mlr_logistic"}, \code{"mlr"},
#'   \code{"ml_logistic"}, and \code{"ml"} name the standard error of the
#'   single-factor model fitted to \code{data} or \code{S} (robust maximum
#'   likelihood or maximum likelihood); \code{"percentile"},
#'   \code{"bca"}, \code{"bootstrap_se"}, and
#'   \code{"bootstrap_se_logistic"} are nonparametric bootstraps (with
#'   \code{data}, or with a single-group \code{fit} that carries its raw
#'   data); \code{"none"} returns the point estimate. The
#'   \code{_logistic} methods build the interval on the logit scale. When
#'   not supplied, the default is \code{"mlr_logistic"} with \code{data},
#'   \code{"ml_logistic"} with \code{S}, and \code{"delta_logistic"} with
#'   \code{fit} and with \code{loadings} that come with
#'   \code{vcov_loadings} or \code{se_loadings}; for \code{loadings}
#'   without either it is \code{"none"}, since nothing is then known about
#'   their sampling variability. A method that does not apply to the
#'   input stops with a message naming the ones that do.
#' @param missing With \code{data}, how incomplete rows are handled:
#'   \code{"fiml"} (the default; full information maximum likelihood) or
#'   \code{"listwise"}. The default analysis is FIML whenever an item value
#'   is missing and the complete data analysis otherwise (data in which no
#'   value is missing, or none but in rows with every item missing, which
#'   both treatments drop, are fitted as they stand), as in
#'   \code{\link{reliability_omega}} and unlike
#'   \code{\link{fa_confirmatory}}; see \emph{Details}. Under
#'   listwise deletion the single-factor model needs more complete rows
#'   than items (with no more, their covariance matrix is singular and the
#'   maximum likelihood fit is not defined), and the function stops
#'   before the fit when there are not. Ordinary FIML instead evaluates
#'   the constrained measurement model. A failure of the unrestricted H1
#'   diagnostic does not alone invalidate its proper identified estimate
#'   or H0-based uncertainty; the result carries an explicit qualification.
#'   Supplied two-stage fits require valid unrestricted first-stage
#'   moments even for a point estimate. See \code{\link{cov_fiml}} for the
#'   sufficient affine-collapse diagnostic and its limits. Bootstrap
#'   resamples retain the family's separate distinct-case screen.
#' @param aux With \code{data}, an optional character vector naming
#'   auxiliary variable columns of \code{data}, entered as saturated
#'   correlates under full information maximum likelihood; supplying it
#'   implies \code{missing = "fiml"}.
#' @param B Number of bootstrap replications when a bootstrap method is
#'   selected, a single whole number. Defaults to \code{10000}, the
#'   reliability family's default; every replication refits a model (see
#'   \emph{Details} for the cost). The fewest accepted depends on the
#'   interval method and the level: at the default \code{conf_level}, 400
#'   for \code{"percentile"} (10 replications beyond each limit), 2000 for
#'   \code{"bca"} (50 beyond each nominal limit, because the adjusted
#'   limits sit farther into the tails), and 100 for \code{"bootstrap_se"}
#'   and \code{"bootstrap_se_logistic"} (a standard deviation of the
#'   replications); at 0.99 the first two are 2000 and 10000. A smaller
#'   \code{B} is refused before anything is computed (class
#'   \code{dmar_count_invalid}). The count is held to the number requested.
#'   Replications set aside afterward are counted and reported in a warning
#'   of class \code{dmar_bootstrap_replications_dropped}, which names the
#'   function and says to request more when those that return a value are
#'   fewer than the number required; when fewer than half of the \code{B}
#'   replications return a value, the interval is withheld (limits and
#'   bootstrap standard error \code{NA}) under a warning of class
#'   \code{dmar_bootstrap_insufficient}. See \code{\link{dmar_bootstrap}},
#'   which also describes \code{options(dmar.bootstrap_exploratory = TRUE)}
#'   for a first look with a smaller \code{B}.
#' @param seed Random number seed for the bootstrap. Defaults to
#'   \code{NULL}, under which the resamples are drawn from, and advance,
#'   the current random number stream; a supplied seed is set for the
#'   call, and the caller's RNG state is restored when it returns.
#' @inheritParams reliability_omega
#' @param ncpus The number of processes over which the bootstrap
#'   replications are spread. The default \code{1} runs them in this R
#'   session; a larger number uses forked workers (\pkg{boot}'s
#'   \code{parallel = "multicore"}) or, on Windows, which cannot fork, a
#'   socket cluster (\code{"snow"}), whose workers load the installed
#'   \pkg{DMAR}. \pkg{boot} draws every resample before any replication
#'   runs, so a given \code{seed} returns the same interval, counts, and
#'   replicates whatever \code{ncpus} is. A session that forks workers and also
#'   starts background processes through \pkg{callr} or \pkg{processx}
#'   needs the environment variable \code{PROCESSX_NOTIFY_OLD_SIGCHLD} set
#'   before \pkg{processx} is loaded (see \code{\link{sem_exact_fit_bootstrap}}). A replication whose solution is
#'   improper is always set aside, since \emph{H} is not defined for it.
#'   Ignored unless a bootstrap method is selected.
#'
#' Bootstrap intervals from supplied two-stage fits currently refuse
#' separately stored auxiliary variables: each refit must resample and
#' carry those values with the items. Use an eligible delta interval
#' where available, or explicitly resample and refit them together.
#'
#' @return The reliability family's table (class
#'   \code{c("dmar_reliability", "dmar_tbl", "data.frame")}), with
#'   columns \code{term} and \code{value} and rows \code{"estimate"}
#'   (\emph{H}), \code{"se"} (its standard error: delta method, or the
#'   standard deviation of the bootstrap replicates; \code{NA} for
#'   \code{ci_method = "none"}), \code{"ci_lower"} and
#'   \code{"ci_upper"}, \code{"conf_level"}, \code{"N"} (the cases the
#'   fit used, in the named group of a multiple group fit; \code{NA} for
#'   \code{loadings} without \code{N}), \code{"N_complete"} (the cases
#'   complete on the factor's items), and \code{"n_items"} (the number of
#'   items). Attributes record the computation: \code{coefficient}
#'   (\code{"H"}), \code{ci_method} (the method used), \code{missing}
#'   (\code{"listwise"}, \code{"fiml"}, the fit's own label such as
#'   \code{"pairwise"}, or \code{"unknown"} for \code{loadings}),
#'   \code{vcov_basis} (\code{"full"} when the delta method used the full
#'   covariance matrix of the estimates, \code{"diagonal"} when only
#'   marginal standard errors were supplied, \code{"none"} when no delta
#'   method standard error was formed), \code{se_type} (the \pkg{lavaan}
#'   standard error type of the fit; \code{NA} for \code{loadings}),
#'   \code{factor} and \code{group} (with \code{fit}; \code{NA}
#'   otherwise), \code{metric} (\code{"observed"},
#'   \code{"latent_response"} for ordered items, or \code{"unknown"} for
#'   \code{loadings}), \code{residual_covariances} (the number of
#'   residual covariances inside the block that entered the weights), and
#'   \code{loadings} (the standardized loadings of the block, named by
#'   item). With \code{loadings}, \code{vcov_supplied} records what came
#'   with them, whether or not an interval was formed: \code{"full"} for
#'   \code{vcov_loadings}, \code{"diagonal"} for \code{se_loadings}, and
#'   \code{"none"} for neither, in which case the printed table says what
#'   an interval would need. With \code{data}, \code{aux} is recorded
#'   when supplied. A
#'   bootstrap call also records \code{B} and the family's
#'   \code{bootstrap} attribute (the replications attempted, returned,
#'   and set aside, with the reasons; see \code{\link{reliability_omega}}),
#'   which here also carries \code{bias}, the mean of the returned
#'   replicates minus the estimate (\code{NA} when no interval is
#'   formed). Its \code{no_solution} counts the
#'   resamples of three items whose correlations admit no single-factor
#'   solution, which are among the \code{improper} ones (on the data
#'   path; zero with \code{fit}, whose refits are not checked for it),
#'   and \code{improper_kept} is zero, since \emph{H} is not defined for
#'   an improper solution.
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
#' Efron, B., & Tibshirani, R. J. (1993). \emph{An introduction to the
#'   bootstrap}. New York, NY: Chapman & Hall/CRC.
#'
#' Gabler, S., & Raykov, T. (2017). Evaluation of maximal reliability for
#'   unidimensional measuring instruments with correlated errors.
#'   \emph{Structural Equation Modeling, 24}(1), 104--111.
#'   \doi{10.1080/10705511.2016.1159916}
#'
#' Graham, J. W. (2003). Adding missing-data-relevant variables to
#'   FIML-based structural equation models. \emph{Structural Equation
#'   Modeling, 10}(1), 80--100. \doi{10.1207/S15328007SEM1001_4}
#'
#' Hancock, G. R., & Mueller, R. O. (2001). Rethinking construct
#'   reliability within latent variable systems. In R. Cudeck, S. du
#'   Toit, & D. Sörbom (Eds.), \emph{Structural equation modeling:
#'   Present and future} (pp. 195--216). Scientific Software
#'   International.
#'
#' Kelley, K. (2007). Confidence intervals for standardized effect sizes:
#'   Theory, application, and implementation. \emph{Journal of Statistical
#'   Software, 20}(8), 1--24. \doi{10.18637/jss.v020.i08}
#'
#' Kelley, K., & Pornprasertmanit, S. (2016). Confidence intervals for
#'   population reliability coefficients: Evaluation of methods,
#'   recommendations, and software for composite measures.
#'   \emph{Psychological Methods, 21}, 69--92. \doi{10.1037/a0040086}
#'
#' Raykov, T. (1997). Estimation of composite reliability for congeneric
#'   measures. \emph{Applied Psychological Measurement, 21}(2), 173--184.
#'   \doi{10.1177/01466216970212006}
#'
#' Rodriguez, A., Reise, S. P., & Haviland, M. G. (2016). Evaluating
#'   bifactor models: Calculating and interpreting statistical indices.
#'   \emph{Psychological Methods, 21}(2), 137--150.
#'   \doi{10.1037/met0000045}
#'
#' @seealso \code{\link{reliability_omega}} for the reliability of the
#'   unit-weighted composite, \code{\link{reliability}} for the family's
#'   general entry point, \code{\link{factor_scores}} for the scores whose
#'   reliability \emph{H} is, \code{\link{fa_confirmatory}} for the fit, and
#'   \code{\link{bifactor_indices}} for the bifactor index of Rodriguez,
#'   Reise, and Haviland (2016) and the maximal reliability of each factor
#'   of a bifactor model.
#'
#' @examples
#' verbal <- c("t5_general_information", "t6_paragraph_comprehension",
#'             "t7_sentence", "t8_word_classification", "t9_word_meaning")
#'
#' # From raw item scores: the single-factor model is fitted here, and the
#' # interval uses its robust maximum likelihood standard error on the
#' # logit scale.
#' reliability_H(data = holzinger_swineford[, verbal])
#'
#' # The omega of the unit-weighted sum from the same model, for
#' # comparison; H is the ceiling that optimal weights could reach.
#' reliability_omega(data = holzinger_swineford[, verbal],
#'                   denominator = "model_implied")
#'
#' # From a fitted model with several factors, name the factor. The
#' # standard error is the fit's own, from the full covariance matrix of
#' # its estimates (requires lavaan).
#' fit <- lavaan::cfa(
#'   "visual  =~ t1_visual_perception + t2_cubes + t4_lozenges
#'    textual =~ t6_paragraph_comprehension + t7_sentence + t9_word_meaning
#'    speed   =~ t10_addition + t12_counting_groups_of_dots +
#'               t13_straight_and_curved_capitals",
#'   data = holzinger_swineford, std.lv = TRUE)
#' h <- reliability_H(fit, factor = "textual")
#' h
#' attr(h, "loadings")
#'
#' # From the standardized loadings a paper reports, without an interval.
#' reliability_H(loadings = c(0.6, 0.7, 0.8), ci_method = "none")
#'
#' # The broom verbs.
#' generics::tidy(h)
#' generics::glance(h)
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @keywords htest
#'
#' @family reliability
#'
#' @export
reliability_H <- function(loadings = NULL, se_loadings = NULL,
                          conf_level = 0.95,
                          fit = NULL, factor = NULL, group = NULL,
                          data = NULL, S = NULL, N = NULL,
                          vcov_loadings = NULL,
                          ci_method = c("mlr_logistic", "mlr",
                                        "ml_logistic", "ml",
                                        "delta_logistic", "delta",
                                        "bootstrap_se",
                                        "bootstrap_se_logistic",
                                        "percentile", "bca", "none"),
                          missing = c("fiml", "listwise"),
                          aux = NULL,
                          B = 10000,
                          seed = NULL,
                          ncpus = 1L,
                          progress = NULL) {
  # A fitted model passed in the first position lands in `loadings` and is
  # taken as `fit`.
  if (inherits(loadings, "lavaan")) {
    if (!is.null(fit)) {
      .dmar_refuse_argument(
        "Supply the fitted model once, as 'fit'.", arg = "fit")
    }
    fit <- loadings
    loadings <- NULL
  }

  # ---- which input --------------------------------------------------------
  inputs <- c(loadings = !is.null(loadings), fit = !is.null(fit),
              data = !is.null(data), S = !is.null(S))
  if (sum(inputs) != 1L) {
    .dmar_refuse_argument(c(
      "Supply exactly one of 'loadings', 'fit', 'data', or 'S' (with ",
      "'N')."), arg = c("loadings", "fit", "data"))
  }
  input <- names(inputs)[inputs]
  missing_supplied <- !missing(missing)
  ci_supplied <- !missing(ci_method)
  missing <- match.arg(missing)

  # Arguments that belong to another input are refused rather than
  # ignored, so a call never reads as doing something it does not.
  given <- c(se_loadings = !is.null(se_loadings),
             vcov_loadings = !is.null(vcov_loadings),
             factor = !is.null(factor), group = !is.null(group),
             N = !is.null(N), missing = missing_supplied,
             aux = !is.null(aux))
  allowed <- switch(input,
                    loadings = c("se_loadings", "vcov_loadings", "N"),
                    fit = c("factor", "group"),
                    data = c("missing", "aux"),
                    S = c("N", "missing"))
  stray <- setdiff(names(given)[given], allowed)
  if (length(stray)) {
    .dmar_refuse_argument(c(
      "With '", input, "', ",
      paste0("'", stray, "'", collapse = " and "),
      if (length(stray) == 1L) " does" else " do",
      " not apply: 'se_loadings', 'vcov_loadings', and 'N' go with ",
      "'loadings'; 'factor' and 'group' with 'fit'; 'missing' and ",
      "'aux' with 'data'; 'N' with 'S'."),
      arg = c("se_loadings", "vcov_loadings", "N"))
  }

  if (!is.numeric(conf_level) || length(conf_level) != 1L ||
      is.na(conf_level) || conf_level <= 0 || conf_level >= 1) {
    .dmar_refuse_argument(
      "'conf_level' must be a single number in (0, 1).", arg = "conf_level")
  }

  # ---- the interval method, resolved by input -----------------------------
  # Standardized loadings with neither their covariance matrix nor their
  # standard errors carry nothing an interval can be formed from, so such
  # a call returns the estimate alone.
  ci_method <- if (ci_supplied) {
    match.arg(ci_method)
  } else {
    switch(input, data = "mlr_logistic", S = "ml_logistic",
           fit = "delta_logistic",
           loadings = if (is.null(vcov_loadings) && is.null(se_loadings)) {
             "none"
           } else {
             "delta_logistic"
           })
  }
  bootstrap_methods <- c("bootstrap_se", "bootstrap_se_logistic",
                         "percentile", "bca")
  delta_methods <- c("delta_logistic", "delta", "ml_logistic", "ml",
                     "mlr_logistic", "mlr")
  is_boot <- ci_method %in% bootstrap_methods
  progress <- .check_progress(progress)
  allowed_ci <- switch(
    input,
    loadings = c("delta_logistic", "delta", "none"),
    fit = c(delta_methods, bootstrap_methods, "none"),
    data = c("mlr_logistic", "mlr", "ml_logistic", "ml",
             bootstrap_methods, "none"),
    S = c("ml_logistic", "ml", "none"))
  if (!ci_method %in% allowed_ci) {
    why <- switch(
      input,
      loadings = paste0("standardized loadings carry no model to refit or ",
                        "to name a standard error for"),
      data = paste0("the single-factor model is fitted here, so name its ",
                    "standard error (\"delta\" and \"delta_logistic\" use ",
                    "the covariance matrix that comes with a supplied ",
                    "'fit' or 'loadings')"),
      S = paste0("the robust standard error and the bootstrap need raw ",
                 "'data'; a covariance matrix is not sufficient"))
    .dmar_refuse_argument(c(
      "ci_method = \"", ci_method, "\" does not apply to '", input,
      "': ", why, ". Use ", paste0("\"", allowed_ci, "\"",
                                   collapse = ", "), "."), arg = NULL)
  }
  if (is_boot) .check_boot_B(B, ci_method, conf_level)

  # ---- standardized loadings ----------------------------------------------
  if (input == "loadings") {
    if (!is.numeric(loadings) || length(loadings) < 2L ||
        anyNA(loadings) || any(abs(loadings) >= 1)) {
      .dmar_refuse_argument(c(
        "'loadings' must be two or more standardized loadings, each ",
        "in (-1, 1)."), arg = "loadings")
    }
    J <- length(loadings)
    if (!is.null(se_loadings) && !is.null(vcov_loadings)) {
      .dmar_refuse_argument(
        "Supply 'vcov_loadings' or 'se_loadings', not both.",
        arg = c("vcov_loadings", "se_loadings"))
    }
    V <- NULL
    basis <- "none"
    if (!is.null(vcov_loadings)) {
      V <- as.matrix(vcov_loadings)
      if (!is.numeric(V) || !identical(dim(V), c(J, J)) || anyNA(V) ||
          !isSymmetric(unname(V), tol = 1e-8) ||
          min(eigen(V, symmetric = TRUE, only.values = TRUE)$values) <
            -1e-10 * max(1, abs(V))) {
        .dmar_refuse_argument(c(
          "'vcov_loadings' must be the symmetric, positive ",
          "semidefinite ", J, " x ", J, " covariance matrix of the ",
          "loading estimates."), arg = "vcov_loadings")
      }
      basis <- "full"
    } else if (!is.null(se_loadings)) {
      if (!is.numeric(se_loadings) || length(se_loadings) != J ||
          anyNA(se_loadings) || any(se_loadings < 0)) {
        .dmar_refuse_argument(c(
          "'se_loadings' must be ", J, " nonnegative standard errors, ",
          "one per loading."), arg = "se_loadings")
      }
      V <- diag(se_loadings^2, J)
      basis <- "diagonal"
    }
    if (!is.null(N) &&
        (!is.numeric(N) || length(N) != 1L || is.na(N) || N <= 0)) {
      .dmar_refuse_argument("'N' must be a single positive number.", arg = "N")
    }
    # What came with the loadings, kept apart from the basis of a standard
    # error, which ci_method = "none" leaves "none" either way: the
    # printed table says what an interval would need only when nothing
    # did.
    supplied <- basis
    if (ci_method == "none") {
      V <- NULL
      basis <- "none"
    } else if (is.null(V)) {
      .dmar_refuse_argument(c(
        "The delta method interval needs the sampling covariance ",
        "matrix of the loadings: supply 'vcov_loadings' (or ",
        "'se_loadings', which ignores their covariances), or set ",
        "ci_method = \"none\" for the point estimate."),
        arg = c("vcov_loadings", "se_loadings"))
    }
    h <- .H_from_loadings(loadings, V)
    se <- sqrt(h$var)
    lim <- .H_wald(h$H, se, conf_level, ci_method == "delta_logistic")
    if (basis == "diagonal") .H_warn_marginal()
    N_out <- if (is.null(N)) NA_real_ else N
    out <- .relia_result(
      estimate = h$H, se = se, lower = lim$lower, upper = lim$upper,
      conf_level = conf_level, N = N_out, J = J, coefficient = "H",
      ci_method = ci_method, N_complete = N_out, missing = "unknown")
    item_names <- names(loadings)
    if (is.null(item_names)) item_names <- paste0("item", seq_len(J))
    return(.H_attributes(
      out, vcov_basis = basis, se_type = NA_character_,
      factor = NA_character_, group = NA_character_, metric = "unknown",
      residual_covariances = 0L,
      loadings = stats::setNames(as.numeric(loadings), item_names),
      vcov_supplied = supplied))
  }

  if (!requireNamespace("lavaan", quietly = TRUE)) {
    .dmar_refuse(c(
      "Coefficient H from a model requires the 'lavaan' package. ",
      "Install with install.packages(\"lavaan\")."),
      "package_required", package = "lavaan")
  }

  # ---- a fitted model -----------------------------------------------------
  if (input == "fit") {
    if (!inherits(fit, "lavaan")) {
      .dmar_refuse_argument(c(
        "'fit' must be a fitted lavaan model (from lavaan::cfa(), ",
        "lavaan::sem(), or fa_confirmatory(output = \"fit\"))."), arg = "fit")
    }
    opts <- lavaan::lavInspect(fit, "options")
    h1_diagnostic <- .fiml_fit_h1_diagnostic(fit)
    if (!is.null(h1_diagnostic) &&
        opts$missing %in% c("two.stage", "robust.two.stage")) {
      .dmar_refuse_argument(c(
        "Two-stage coefficient H requires unrestricted first-stage ",
        "moments, even for a point estimate. ", h1_diagnostic), arg = NULL)
    }
    if (!is.null(h1_diagnostic) && ci_method %in% delta_methods &&
        !.fiml_h0_vcov_supported(fit)) {
      .dmar_refuse(c(
        "The unrestricted H1 diagnostic failed, and the supplied fit's ",
        "custom SE/information combination has not been verified to be ",
        "H0-based. Refit with the standard observed-Hessian or default ",
        "MLR Huber-White SEs, or request ci_method = \"none\" for the ",
        "point estimate alone. ", h1_diagnostic), "computation_failed")
    }
    if (ci_method %in% c("ml_logistic", "ml", "mlr_logistic", "mlr")) {
      fit_kind <- if (identical(opts$estimator, "ML") &&
                      identical(opts$se, "standard")) {
        "ml"
      } else if (identical(opts$estimator, "ML") &&
                 opts$se %in% c("robust.sem", "robust.huber.white")) {
        "mlr"
      } else {
        NA_character_
      }
      if (!identical(sub("_logistic$", "", ci_method), fit_kind)) {
        .dmar_refuse_argument(c(
          "ci_method = \"", ci_method, "\" names a ",
          if (startsWith(ci_method, "mlr")) "robust " else "",
          "maximum likelihood standard error, but the fit carries ",
          "se = \"", opts$se, "\" (estimator ", opts$estimator, "). ",
          "With 'fit' the interval uses the fit's own covariance ",
          "matrix of the estimates: use ci_method = \"delta_logistic\" ",
          "or \"delta\", or refit the model with the standard error ",
          "you want."), arg = "fit")
      }
    }
    need_se <- ci_method %in% delta_methods
    if (need_se && identical(opts$se, "none")) {
      .dmar_refuse_argument(c(
        "The fit carries no standard errors (se = \"none\"), so the ",
        "delta method has no covariance matrix to use. Refit with ",
        "standard errors, or use ci_method = \"none\" or a bootstrap ",
        "method."), arg = NULL)
    }
    h <- .H_from_fit(fit, factor = factor, group = group, se = need_se)
    boot_out <- NULL
    if (is_boot) {
      .fiml_check_boot_aux(fit, "Coefficient H bootstrap")
      if (lavaan::lavInspect(fit, "ngroups") > 1L) {
        .dmar_refuse_argument(c(
          "The bootstrap is implemented for single-group fits; use ",
          "ci_method = \"delta_logistic\" for a multiple group model, ",
          "or fit each group's model separately and bootstrap it."), arg = NULL)
      }
      # Cases are resampled one at a time and every replication refits
      # the model to its rows alone, so a fit whose cases are not
      # independent, or whose rows carry sampling weights, would be
      # bootstrapped as a different model: an interval that ignores the
      # clustering, or one around the unweighted H. The delta method uses
      # the fit's own covariance matrix of the estimates, which carries
      # either. ci_bifactor_indices() refuses both for the same reason.
      cluster <- tryCatch(lavaan::lavInspect(fit, "cluster"),
                          error = function(e) character(0))
      if (length(cluster)) {
        .dmar_refuse_argument(c(
          "'fit' names a cluster variable (",
          paste(cluster, collapse = ", "), "), so its cases are not ",
          "independent, and a bootstrap that resamples them one at a ",
          "time would ignore the clustering. Use ci_method = ",
          "\"delta_logistic\", whose standard error is the one the fit ",
          "was estimated with."), arg = "fit")
      }
      weights <- tryCatch(lavaan::lavInspect(fit, "sampling.weights"),
                          error = function(e) NULL)
      if (length(weights)) {
        .dmar_refuse_argument(c(
          "'fit' uses sampling weights, which a resample of its rows ",
          "does not carry, so every replication would be an unweighted ",
          "fit and the interval would describe the unweighted H. Use ",
          "ci_method = \"delta_logistic\", whose standard error is the ",
          "one the fit was estimated with."), arg = "fit")
      }
      boot_data <- .H_fit_data(fit)
      point_fn <- .H_refit_fun(fit, h$setup)
      # The replications must be refits of the model behind the estimate:
      # refit to the rows of the fit itself, the point function must
      # reproduce it, as in ci_bifactor_indices(). This catches anything
      # the refit through lavaan::lavaan() does not carry over.
      check <- tryCatch(point_fn(boot_data), error = function(e) e)
      # A fit whose likelihood has no maximum is refused as DMAR's own
      # fits are (.no_max_refuse()); 'fit' is the user's lavaan fit, so
      # there is no DMAR syntax to print.
      if (inherits(check, "dmar_no_maximum")) {
        .no_max_refuse("reliability_H", .no_max_fields(check))
      }
      if (inherits(check, "error")) {
        .dmar_refuse(c(
          "The model could not be refit to the data stored with 'fit': ",
          conditionMessage(check)), "computation_failed")
      }
      gap <- abs(as.numeric(check) - h$H)
      if (!isTRUE(gap <= 1e-6)) {
        .dmar_refuse_argument(c(
          "Refitting the model to the data stored with 'fit' did not ",
          "reproduce its H",
          if (is.finite(gap)) {
            sprintf(" (%.6f against %.6f)", as.numeric(check), h$H)
          } else {
            ""
          },
          ", so the bootstrap replications would not be refits of this ",
          "model. Use ci_method = \"delta_logistic\"."), arg = "fit")
      }
      boot_out <- .bootstrap_ci(
        fun = "reliability_H",
        data = boot_data, point_fn = point_fn,
        B = B, conf_level = conf_level, kind = ci_method, seed = seed,
        ncpus = ncpus, improper = "drop", improper_note = .H_improper_note(),
        progress = progress)
    }
    return(.H_finish(h, ci_method, boot_out, conf_level, B, aux = NULL,
                     factor = h$setup$factor, group = h$setup$group,
                     h1_diagnostic = h1_diagnostic))
  }

  # ---- raw data or a covariance matrix ------------------------------------
  if (input == "data" && !is.null(aux) && missing_supplied &&
      missing == "listwise") {
    .dmar_refuse_argument(c(
      "'aux' cannot be combined with missing = \"listwise\": an ",
      "auxiliary variable only carries information through full ",
      "information maximum likelihood. Drop 'aux' or set ",
      "missing = \"fiml\"."), arg = "aux")
  }
  # The default is FIML when a value is missing (supplying aux implies it
  # in any case); complete data and a covariance matrix are analyzed as
  # they stand, as reliability_omega() analyzes them, so that the two fit
  # the same single-factor model the same way (.relia_resolve_missing()).
  # Every interval H takes from data is available under FIML.
  missing <- .relia_resolve_missing(missing, missing_supplied,
                                    if (input == "data") data, aux)
  if (input == "S" && missing == "fiml") {
    .dmar_refuse_argument(c(
      "missing = \"fiml\" requires raw 'data': a covariance matrix has ",
      "no incomplete cases for full information maximum likelihood to ",
      "use."), arg = "data")
  }
  inp <- .relia_resolve_inputs(data = if (input == "data") data else NULL,
                               S = if (input == "S") S else NULL, N = N,
                               allow_S = TRUE, missing = missing, aux = aux,
                               single_factor = "H", needs_saturated = FALSE)
  if (inp$J < 3L) {
    .dmar_refuse_argument(c(
      "Coefficient H from ",
      if (input == "data") "raw data" else "a covariance matrix",
      " needs three or more items: the single-factor congeneric model ",
      "is not identified with two. With a fitted model in which the ",
      "factor is identified, supply 'fit'."), arg = "fit")
  }
  fiml <- missing == "fiml"
  spec <- switch(
    ci_method,
    ml = , ml_logistic = list(estimator = "ML", se = "standard"),
    mlr = , mlr_logistic = list(estimator = "MLR",
                                se = if (fiml) "robust.huber.white"
                                     else "robust.sem"),
    list(estimator = "ML", se = "none"))
  # The single-factor model behind H, fitted as reliability_omega() fits
  # it; d = NULL fits the covariance matrix.
  one_factor <- function(d, se) {
    .omega_fit_cfa(data = d, h1_diagnostic = inp$h1_diagnostic,
                   S = if (is.null(d)) inp$S else NULL,
                   N = if (is.null(d)) inp$N else NULL,
                   equal_loadings = FALSE, estimator = spec$estimator,
                   se = se, missing = if (fiml) "ml" else "listwise",
                   aux = inp$aux)
  }
  of <- one_factor(inp$data, spec$se)
  # Three items whose correlations have a product that is not positive
  # admit no single-factor solution at all; the fit is refused before it
  # is attempted (.triad_no_solution()), with the same message
  # reliability_omega() gives.
  if (isTRUE(of$no_solution)) {
    .dmar_refuse(.no_solution_message(of$triad), "computation_failed")
  }
  # A fit whose likelihood has no maximum is refused with the lavaan code
  # of the single-factor fit, as reliability_omega() refuses it
  # (.no_max_refuse()).
  if (!is.null(of$no_maximum)) .no_max_refuse("reliability_H", of$no_maximum)
  # A fit lavaan stopped with an error is reported with lavaan's message,
  # as reliability_omega() reports it; through version 1.0.1 any failure
  # was reported as a model that did not converge.
  if (!of$converged) {
    .dmar_refuse(
      .relia_fit_failure(of, model = "single-factor", missing = missing),
      "computation_failed")
  }
  if (isTRUE(of$improper)) {
    .dmar_refuse_argument(c(
      "The single-factor solution is improper (a Heywood case): an ",
      "item's residual variance is not positive, so H is not defined ",
      "for it. Reconsider whether one factor describes these items."),
      arg = NULL)
  }
  h <- .H_from_fit(of$fit, factor = "f1", se = ci_method %in% delta_methods)
  # The fit names the items by position when auxiliaries are present (and
  # when the input had no names); the block is the input's item columns,
  # in order.
  if (!is.null(inp$items)) names(h$loadings) <- inp$items
  h$N <- inp$N
  h$N_complete <- inp$N_complete
  h$missing <- missing
  boot_out <- NULL
  if (is_boot) {
    # A resample whose solution is improper has no H, and neither has one
    # of three items whose correlations admit no single-factor solution
    # (.omega_fit_none() records it as not converged); both are counted as
    # improper, as reliability_omega() counts them, the second also in
    # the attribute's `no_solution` (the engine's own reason for it), and
    # a fit that did not converge among the non-finite ones.
    # A replication in which an item (or auxiliary variable) takes a
    # single value is refused before the fit, as the data are
    # (.relia_checked_point()).
    point_fn <- .relia_checked_point(function(d) {
      bf <- one_factor(d, "none")
      # A replication whose likelihood has no maximum stops, and the
      # engine counts it among those that raised an error
      # (.no_max_replication()).
      if (!is.null(bf$no_maximum)) .no_max_replication(bf$no_maximum)
      verdict <- .omega_fit_verdict(bf)
      if (verdict == "not_converged") return(NA_real_)
      if (verdict == "no_solution") {
        return(structure(NA_real_, dmar_reason = "no_solution"))
      }
      if (verdict != "proper") {
        return(structure(NA_real_, dmar_reason = "improper"))
      }
      .H_replicate(h$setup, bf$fit)
    }, missing = missing,
    items = if (missing == "listwise") NULL else inp$items,
    aux = inp$aux, fitted = TRUE, needs_saturated = missing == "two_stage")
    boot_out <- .bootstrap_ci(fun = "reliability_H",
                              data = inp$data, point_fn = point_fn, B = B,
                              conf_level = conf_level, kind = ci_method,
                              seed = seed, ncpus = ncpus, improper = "drop",
                              improper_note = .H_improper_note(),
                              progress = progress)
  }
  .H_finish(h, ci_method, boot_out, conf_level, B, aux = inp$aux,
            factor = NA_character_, group = NA_character_,
            h1_diagnostic = inp$h1_diagnostic)
}
