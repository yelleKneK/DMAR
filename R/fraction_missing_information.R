#' Fraction of Missing Information for Each Parameter of a FIML Fit
#'
#' Estimates the fraction of missing information (FMI) for every free and
#' every defined parameter of a model fit by full information maximum
#' likelihood (FIML), and states each in the sample size and interval
#' width terms of accuracy in parameter estimation (AIPE). The FMI of a
#' parameter is the proportion of the sampling variance of its estimate
#' that is due to the missing values (Savalei & Rhemtulla, 2012): zero
#' when nothing is missing, and closer to one the less the observed data
#' say about the parameter. Unlike the percentage of missing values, it
#' depends on the parameter, the model, the correlations among the
#' variables, and the missing data mechanism, so it is the number that
#' says how much the missing values cost the estimate a reader cares
#' about, and how many more cases a study with the same missing data
#' would need to estimate it as precisely as complete data would.
#'
#' @param fit A converged full information maximum likelihood fit: a
#'   \pkg{lavaan} fit with \code{missing = "ml"} (or \code{"ml.x"}), for
#'   example from \code{fa_confirmatory(..., missing = "fiml", output = "fit")} or
#'   \code{lavaan::cfa(..., missing = "ml")}, or a fit returned by
#'   \code{\link{mlmr}} or \code{\link{mlmr_mv}}, whose \pkg{lavaan} fit is
#'   used. Requires the \pkg{lavaan} package.
#' @param information The estimate of the information matrix behind both
#'   standard errors: \code{"unstructured"} (the default),
#'   \code{"structured"}, or \code{"hessian"}, the three estimators of
#'   Chen and Savalei (2021). \code{"hessian"} is \pkg{lavaan}'s default:
#'   for a fit made with \pkg{lavaan}'s default information settings it
#'   agrees with \code{lavaan::parameterEstimates(fit, fmi = TRUE)} to
#'   about \eqn{10^{-7}} when the fit runs \pkg{lavaan}'s EM algorithm as
#'   far as the computation here does (\code{em_h1_args = list(tol =
#'   1e-10)}), and otherwise to the precision of \pkg{lavaan}'s EM, which
#'   stops earlier; at extreme scales \pkg{lavaan}'s own values are the
#'   ones that are off. See \emph{Details}.
#' @param ci_method \code{"none"} (the default) or \code{"percentile"}, a
#'   percentile bootstrap interval for each FMI. See \emph{Details}.
#' @param conf_level Confidence level of the percentile interval (default
#'   \code{0.95}); used when \code{ci_method = "percentile"}.
#' @param B Number of bootstrap replications when \code{ci_method =
#'   "percentile"}, a single whole number. Default \code{1000}, smaller
#'   than the package's usual \code{10000} because every replication refits
#'   the model and then makes the two refits behind each FMI. The fewest
#'   accepted is 400 at the default \code{conf_level} (200 at 0.90, and
#'   2000 at 0.99, which the default does not meet); a smaller \code{B} is
#'   refused before anything is computed (class \code{dmar_count_invalid}).
#'   See \emph{Details} and \code{\link{dmar_bootstrap}}.
#' @param seed Optional integer seed for the bootstrap. The default
#'   \code{NULL} uses the current state of the random number generator; a
#'   supplied seed is set internally and the prior state restored on exit.
#'
#' @details
#' \strong{What the fraction of missing information is.} The missing
#' information principle of Orchard and Woodbury (1972) splits the
#' information that complete data would carry about the parameters,
#' \eqn{J_X}, into the information the observed data carry, \eqn{J_Y},
#' and the information lost with the missing values, \eqn{J_{X|Y}}:
#' \eqn{J_X = J_Y + J_{X|Y}}. In large samples the inverse of each
#' information matrix is the covariance matrix of the corresponding
#' estimates, so for parameter \eqn{j}
#' \deqn{\mathrm{FMI}_j = 1 - \frac{\{J_X^{-1}\}_{jj}}{\{J_Y^{-1}\}_{jj}}
#'   = 1 - \frac{\mathrm{SE}_{C,j}^2}{\mathrm{SE}_{O,j}^2},}
#' where \eqn{\mathrm{SE}_{O,j}} is the standard error of the FIML
#' estimate and \eqn{\mathrm{SE}_{C,j}} is the standard error it would
#' have had if the same
#' cases had been observed completely (Savalei & Rhemtulla, 2012,
#' Eq. 4). One minus the FMI is the asymptotic efficiency of the FIML
#' estimate relative to the complete data estimate. The FMI is not the
#' proportion of missing values. When the values are missing completely
#' at random, the FMI of a parameter is bounded by the largest univariate
#' and pairwise rates of missing values, and it falls below them when the
#' observed variables carry information about the missing ones; when the
#' values are missing at random it can be well above the missing rate
#' (Savalei & Rhemtulla, 2012, pp. 489--490). Savalei and Rhemtulla
#' (2012, p. 489) recommend reading the FMI as an effect size for the
#' impact of the missing values on each parameter rather than against a
#' cutoff.
#'
#' \strong{The FMI in planning terms.} A Wald interval's width is
#' proportional to the standard error, and the standard error to
#' \eqn{1/\sqrt{N}}, so the precision the missing values cost a parameter
#' can be stated as sample sizes and widths (Savalei & Rhemtulla, 2012,
#' p. 488). With \emph{N} the number of cases behind the parameter, the
#' table reports
#' \itemize{
#'   \item \code{N_effective} \eqn{= N(1 - \mathrm{FMI})}, the number of
#'     complete cases that would estimate the parameter as precisely as
#'     the \emph{N} incomplete cases did;
#'   \item \code{N_target} \eqn{= N / (1 - \mathrm{FMI})}, the number of
#'     cases with the same pattern and mechanism of missing values that
#'     would match the precision of \emph{N} complete cases; and
#'   \item \code{width_inflation} \eqn{= 1 / \sqrt{1 - \mathrm{FMI}} =
#'     \mathrm{SE}_O / \mathrm{SE}_C}, how many times wider the confidence interval is than
#'     it would have been with complete data at the same \emph{N}. An FMI
#'     of .5 widens the interval by a factor of 1.41, and an FMI of .75
#'     doubles it.
#' }
#' These are the statements of accuracy in parameter estimation, which
#' plans a sample size for the width of a confidence interval rather than
#' for power (Maxwell, Kelley, & Rausch, 2008), applied to the cost of the
#' missing values. A sample size planned for a confidence interval width under complete
#' data (by \code{\link{ss_aipe_reliability}}, \code{\link{ss_aipe_sem_path}},
#' or another \code{ss_aipe_*} planner) is, to a first approximation,
#' divided by one minus the FMI when the study is expected to have the
#' same missing data as the pilot or comparable study whose FMI is used.
#' The approximation is asymptotic: it treats the FMI as known and the
#' interval as a Wald interval.
#'
#' \strong{How the two standard errors are computed.} Both come from
#' refits made at the estimates of \code{fit} with the optimizer switched
#' off, under normal theory maximum likelihood standard errors and one
#' type of information for both, which is how \pkg{lavaan} computes the
#' FMI since version 0.6-7 (Chen & Savalei, 2021, p. 3).
#' \eqn{\mathrm{SE}_O} comes from the fit's own data under its missing
#' data option. \eqn{\mathrm{SE}_C} comes from a fit to the unrestricted FIML estimates of the means and
#' covariances (the saturated model's estimates, found by the EM
#' algorithm), with the same \emph{N} and no rescaling of the covariance
#' matrix. Savalei and Rhemtulla (2012) evaluated the complete data
#' standard errors at the model implied moments, with expected
#' information for the complete data and observed information for the
#' incomplete data; Chen and Savalei (2021, p. 3) report that the
#' mismatch produced negative estimates, particularly in small samples,
#' and the published recipe is not zero for complete data when the model
#' misfits. The matched computation used here gives an FMI of zero, to
#' rounding error, for every parameter when no value is missing. The
#' unrestricted moments are the FIML refit's own EM estimates, computed
#' in its own order of the variables, so a fit made with \pkg{lavaan}'s
#' \code{ov_order = "data"}, which orders its variables as the data do,
#' gets the FMIs of the same fit made in the model's order (but see
#' \emph{Estimates that are not at the maximum} for such a fit that
#' \pkg{lavaan} rescales). Both refits are made with the variables in
#' units in which their standard deviations are near 1, and their
#' standard errors are carried back to the original units exactly (see
#' \emph{Units of measurement}).
#'
#' \strong{Units of measurement.} The FMI is a ratio of two variances of
#' the same estimate, so it does not depend on the units in which the
#' variables are measured, and neither should its estimate. Three steps
#' of the computation do depend on the units when they are carried out in
#' the variables' own metric: \pkg{lavaan}'s EM estimate of the
#' unrestricted moments, which stops when no mean or covariance changes by
#' more than a fixed amount, too coarse a rule for variables with small
#' variances; the inversion of the information matrix, which fails in
#' floating point when the variances span many orders of magnitude; and
#' the Hessian estimator, which differentiates the log likelihood
#' numerically with a step suited to standard deviations near 1. Both
#' refits are therefore made with every observed variable divided by the
#' power of two nearest its standard deviation (over all groups), the
#' factors \pkg{lavaan} itself divides by when it rescales a badly scaled
#' fit, and with the parameter table carried through the same map: every
#' estimate and fixed value is divided by the factor its units imply,
#' every defined parameter and constraint is rewritten in terms of the
#' parameters in their original units, and the standard errors are mapped
#' back. Division by a power of two is exact in floating point, so the
#' refits hold the fit's estimates exactly, and with no missing values
#' their standard errors equal \pkg{lavaan}'s own to about ten digits.
#' The EM is run until no mean or covariance changes by more than
#' \eqn{10^{-10}} times the smallest variance, about \eqn{10^{-10}} in that
#' metric, where \pkg{lavaan}'s default is an absolute \eqn{10^{-5}}: EM
#' converges slowly where much information is missing, and stopped at
#' \pkg{lavaan}'s tolerance the FMIs of a factor model of \pkg{lavaan}'s
#' nine Holzinger and Swineford tests, with values missing at random on
#' two of them, moved by up to \eqn{2 \times 10^{-5}} when the EM was
#' run further. It is run for at most 10,000 iterations (or the
#' \code{max_iter} of the fit's \code{em_h1_args}, if that is larger), and
#' the computation warns when it stops there short of its tolerance, which
#' happens rarely and only where nearly all the information about some
#' combination of the means and covariances is missing: in a one factor
#' model of 12 items at \emph{N} = 50 with 30 percent of the values of half
#' of them missing, one data set needed about 12,500 iterations, and its
#' FMIs at the 10,000th were within \eqn{3 \times 10^{-7}} of the converged
#' ones (and within .002 at \pkg{lavaan}'s tolerance). Where the estimate
#' does not exist the computation stops (see \emph{Fits that are
#' refused}). With the fit's estimates held fixed, the three
#' estimators' FMIs then agreed to within \eqn{2 \times 10^{-9}} when
#' every variable was multiplied by a constant from .001 to 1,000,000, in
#' that factor model and in a regression of one variable on three with
#' values missing at random. Computed from fits to the rescaled data, the
#' regression's FMIs agreed with those at a standard deviation of 1 to
#' within \eqn{1.3 \times 10^{-7}} when one variable or every variable
#' was multiplied by .001, 10,000, or 100,000, and when the outcome was
#' multiplied by 100,000,000. \pkg{lavaan}'s own
#' \code{parameterEstimates(fit, fmi = TRUE)}, which refits in the
#' variables' own metric, was off by up to .93 (by .04 to .07 for the
#' slopes) with every standard deviation near 100,000, and returned
#' \code{NA} with every standard deviation near 300,000, or with one near
#' 10,000,000.
#'
#' \strong{Estimates that are not at the maximum.} The FMI is computed at
#' the fit's estimates, so it is the FMI of the maximum likelihood
#' estimates only when they are at the maximum. \pkg{lavaan} rescales its
#' own fit only when a standard deviation is .001 or less or 1000 or
#' more, or two differ by a factor of 1000, and a fit to variables
#' measured in small or large units between those limits can stop short
#' of the maximum, although \pkg{lavaan} reports it as converged and its
#' check of the gradient, made in the variables' own metric, passes it.
#' With every standard deviation of that regression near .003 the slopes
#' were off by .14 percent and the FMIs computed at them by up to .001
#' (by .00006 near .01 and .03). The factor model of the nine tests,
#' identified by a marker indicator for each factor (the default of
#' \code{lavaan::cfa()}), with every variable multiplied by 100 had a log
#' likelihood .11 below the maximum and FMIs off by up to .011, and
#' multiplied by 500 a log likelihood 287 below it and FMIs off by .12;
#' identified by fixing each factor's variance (\code{std.lv = TRUE}) it
#' reached the maximum at every scale from .01 to 500. In \pkg{lavaan}
#' 0.7-3 a fit made with \code{ov_order = "data"} and rescaled by
#' \pkg{lavaan} is not the maximum likelihood fit when the columns of the
#' data are not in the model's order: in such a regression the log
#' likelihood was 50 below the maximum, and the slopes' FMIs were .28,
#' .22, and .70 in place of .06, .19, and .70. The computation therefore
#' takes one Newton step from the fit's estimates in the metric of the
#' refits. Where the step moves an estimate by more than \eqn{10^{-4}} of
#' its size (or of 1, for an estimate smaller than 1; fits at the maximum
#' took steps below \eqn{6 \times 10^{-6}} in the examples studied here),
#' it computes the FMIs again one step on, at the maximum to working
#' precision, and warns when they differ from those at the fit's
#' estimates by more than \eqn{10^{-4}}, naming the largest change of an
#' estimate, the gain in the log likelihood, and the largest change of an
#' FMI. The FMIs returned are those at the fit's estimates. For
#' variables in small or large units, or in an order other than the
#' model's, the remedy is to fit the model again with \pkg{lavaan}'s
#' \code{rescale_data = TRUE}, which rescales whatever the units, and with
#' the columns of the data in the model's order: the regression near .003
#' and the factor model multiplied by 100 then gave the FMIs of the
#' variables in their own units to within \eqn{10^{-7}} and
#' \eqn{2 \times 10^{-6}}. The warning names units as the cause, with
#' that remedy, only for a fit \pkg{lavaan} did not rescale with a
#' standard deviation of .01 or less or 100 or more (with every variable
#' multiplied by one constant, the factor model and the regression reached
#' the maximum with standard deviations from .01 to 80), and the order of
#' the variables only for a fit \pkg{lavaan} rescaled; a fit made with
#' \code{rescale_data = TRUE}, or one \pkg{lavaan} rescaled itself, is not
#' told to take it again. A fit can also stop short where the likelihood
#' is flat in some direction, as in a small sample with most of the
#' outcome missing (60 cases, 16 of them observing the outcome of a
#' regression on four predictors), where \pkg{lavaan}'s optimizer stopped
#' about \eqn{5 \times 10^{-4}} standard errors short of the maximum and
#' neither rescaling, a tighter tolerance, nor a restart from the
#' estimates moved it; the FMIs were then within \eqn{3 \times 10^{-5}} of
#' the maximum's, which passes, or in another such data set
#' \eqn{6 \times 10^{-4}} off, which the warning reports without naming
#' units as the cause. The check is not made for a model with an
#' inequality constraint or a parameter at a bound, where the gradient
#' need not be zero at the maximum.
#'
#' The refits stay in the variables' own metric for a model with
#' composites (\code{<~}) or product terms, which \pkg{lavaan} does not
#' rescale either, and for one that shares a label between parameters of
#' variables in different units under \code{ceq.simple = TRUE} (one
#' common factor for every variable is tried first). For such a model the
#' Hessian estimator warns when a standard deviation falls outside .01 to
#' 1000, an information matrix that cannot be inverted stops the
#' computation with a message that says so, and the check for a maximum
#' is not made; dividing the variables by constants that bring their
#' standard deviations nearer 1, which leaves the FMI unchanged, is the
#' remedy for all three.
#'
#' \strong{The three estimators and the default.} Chen and Savalei
#' (2021) studied three versions of the observed information, which
#' coincide for a saturated model such as a regression and estimate the
#' same quantity when the model is correct.
#' \itemize{
#'   \item \code{"hessian"} (their first estimator) uses the numeric
#'     Hessian of the log likelihood, \pkg{lavaan}'s default.
#'     \code{information = "hessian"} agrees with
#'     \code{lavaan::parameterEstimates(fit, fmi = TRUE)} to about
#'     \eqn{10^{-7}} for a fit made with \code{em_h1_args = list(tol =
#'     1e-10)}, so that \pkg{lavaan} runs its EM as far, and otherwise to
#'     the precision of \pkg{lavaan}'s EM, about \eqn{10^{-5}} in the
#'     examples studied here (see \emph{Units of measurement} above).
#'   \item \code{"structured"} (their second) keeps the first, dominant
#'     term of the observed information and evaluates it at the model
#'     implied means and covariances.
#'   \item \code{"unstructured"} (their third) has the same analytic form
#'     but uses the unrestricted estimates of the means and covariances
#'     wherever the computation needs them.
#' }
#' In their simulations the unstructured estimator was much less likely
#' than the others to break down and more precise than the structured
#' one, while the Hessian estimator was the most likely to break down
#' when the sample was small or the missing rate high (Chen & Savalei,
#' 2021, p. 15). It is therefore the default here, which differs from
#' \pkg{lavaan}'s. It also cannot be negative, apart from rounding
#' error: both of its informations are evaluated at the saturated
#' maximum likelihood estimates, where the observed data information of
#' the saturated model is bounded above by its complete data information,
#' and the bound carries through the model Jacobian to every parameter,
#' defined parameters included, so \eqn{0 \le \mathrm{FMI} < 1}. The
#' structured and Hessian estimators evaluate the information away from
#' that point and can fall below zero, far below when the model is
#' misspecified: with the three-factor model for \pkg{lavaan}'s nine
#' Holzinger and Swineford tests and values missing at random on six of
#' them, the Hessian estimates for the speed loadings are near \eqn{-2}
#' and the unstructured ones between .17 and .34. When the model is wrong
#' only the Hessian estimator is consistent for its population value, an
#' advantage that Chen and Savalei (2021, p. 15) note has yet to be shown
#' to matter in practice. Under misspecification the three estimate
#' different quantities; the unstructured one is the information about
#' the parameters carried by the unrestricted moments.
#'
#' \strong{MLR and other standard errors.} For a maximum likelihood fit
#' with other standard errors (estimator \code{"MLR"}, whose standard
#' errors are robust to nonnormality; \code{"MLF"}, whose standard
#' errors use the outer product of the casewise gradients; or a fit with
#' standard errors robust to clustering and no multilevel structure) the
#' estimates are the maximum likelihood estimates, and only the standard
#' errors differ. The FMI is defined through the information of the likelihood, and the information
#' of the working normal likelihood is the bread of the sandwich, so the
#' FMI reported for such a fit is the normal theory FMI at its estimates,
#' the same number the ML fit gives. The source fit's standard errors are
#' named in the \code{source_se} attribute and under the printed table. A
#' robust counterpart would need the complete data meat of the sandwich,
#' which involves fourth moments of the values that were not observed and
#' is not identified from incomplete data without a distributional
#' assumption. The ratio of a normal theory complete data variance to a
#' robust FIML variance is not an FMI, because it attributes the
#' inflation due to nonnormality and misfit to the missing values. The
#' robust and the normal theory fractions coincide when nonnormality
#' inflates the variances by the same factor with and without the
#' missing values, which is an assumption, not a result.
#'
#' \strong{Defined parameters.} A parameter defined with \code{:=}, such
#' as the \code{omega_}, \code{H_}, and \code{ave_} rows that
#' \code{\link{fa_confirmatory}} defines for each factor or an indirect effect,
#' has a delta method standard error under each refit, so its FMI is the
#' same ratio. For a reliability coefficient this is the FMI that a
#' sample size plan for its confidence interval needs.
#'
#' \strong{Multiple groups and the \emph{N} column.} Multiple group fits
#' are supported. The \code{N} of each row is the number of cases behind
#' that parameter: the group's \emph{N} for a parameter of one group; the
#' summed \emph{N} of the groups involved for a parameter constrained
#' equal across groups (by a shared label, an \code{==} constraint
#' written with labels or with defined parameters, or
#' \code{group.equal}), whose FMI is then the same in every group; and,
#' for a defined parameter, the summed \emph{N} of the groups its
#' expression reaches. \emph{N} counts the cases the fit used.
#' \pkg{lavaan} sets aside a case with no observed value on any analyzed
#' variable, as Savalei and Rhemtulla (2012, p. 487) did in their
#' example, and the number set aside is the \code{N_excluded} attribute.
#' In a regression fit directly with \pkg{lavaan}, \code{missing = "ml"}
#' with \code{fixed.x = TRUE} also drops every case missing a predictor
#' before anything is estimated; \code{missing = "ml.x"}, or
#' \code{\link{mlmr}}, which models the predictors jointly by default,
#' keeps those cases.
#'
#' \strong{Auxiliary variables.} A fit with auxiliary variables entered
#' as saturated correlates (Graham, 2003), as \pkg{lavaan}'s \code{aux}
#' argument enters them (and \code{fa_confirmatory(..., aux = )} passes it on), is
#' supported, and the table includes the rows of the auxiliaries'
#' variances, means, and covariances. An auxiliary that shares a factor
#' with the indicators covaries with them through their residuals only,
#' so the residual covariance matrix that includes it is not positive
#' definite, and the covariances can imply correlations above 1 with the
#' residuals, although every variance is positive and the model implied
#' covariance matrix is positive definite. \pkg{lavaan} flags both, and
#' a refit that took its checks at face value would replace the
#' estimates; the refits here keep them.
#'
#' \strong{Precision of the FMI itself.} A sample FMI is a noisy
#' estimate of the population value. In the simulations of Chen and
#' Savalei (2021, p. 14) it was nearly unbiased but varied greatly from
#' sample to sample at small \emph{N}; with 10 to 20 percent of all values
#' missing, samples of several hundred gave reasonably precise estimates,
#' and with 30 percent missing it took more than 1,000.
#' \code{ci_method = "percentile"} resamples the cases the fit used
#' (within groups for a multiple group fit), refits the model with its
#' own options, recomputes every FMI, and takes as the limits the
#' quantiles of the replications at \code{(1 - conf_level) / 2} and
#' \code{(1 + conf_level) / 2} (Efron & Tibshirani, 1993). The
#' properties of this interval have not been studied (Chen and Savalei, 2021, p. 15, list bootstrap intervals
#' for the sample FMI as future work); it shows the sampling variability
#' of the estimate, and its upper limit is the conservative value to
#' carry into a sample size plan.
#'
#' \strong{Bootstrap replications.} The two limits of each parameter's
#' interval are quantiles of the \code{B} replicated FMI values, so the
#' interval is run only when \code{B} puts at least 10 replications beyond
#' each limit (400 at the 95 percent level). A replication in which the FMI
#' of some parameter is undefined (the refit fails, does not converge, or
#' returns an FMI that is not a finite number) is dropped and counted,
#' never replaced, and the interval is computed from the rest. The result
#' records \code{B} (requested), \code{B_used} (the number the interval
#' rests on), \code{n_dropped}, and \code{B_required} (the fewest
#' replications a request may ask for; that count is held to \code{B}, not
#' to \code{B_used}) as attributes, with
#' \code{bootstrap_status} (\code{"complete"}, \code{"dropped"},
#' \code{"exploratory"}, \code{"insufficient_count"}, or
#' \code{"insufficient_share"}). One warning of class
#' \code{dmar_bootstrap_replications_dropped} states the counts when any
#' replication was dropped, with a caution when the dropped replications
#' are at least the 2.5 percent of \code{B} that lies beyond one limit at
#' the 95 percent level, because they could then decide by themselves where
#' that limit falls. When more than half were dropped, the limits of every
#' parameter are \code{NA},
#' the estimate is returned, and the warning has class
#' \code{dmar_bootstrap_insufficient}.
#' \code{options(dmar.bootstrap_exploratory = TRUE)} lets a smaller
#' \code{B} (2 or more) run for a first look or a demonstration; every such
#' call gives a warning of class \code{dmar_bootstrap_exploratory} and the
#' result is not one to report. See \code{\link{dmar_bootstrap}} for the
#' rule and its reasons.
#'
#' Each replication costs a model fit and two refits, so the
#' interval takes about a minute at the default \code{B} for a model the
#' size of the one in the examples; the call is
#' \code{fraction_missing_information(fit, ci_method = "percentile",
#' seed = 113)}.
#'
#' \strong{Fits that are refused.} The computation needs a converged and
#' unweighted maximum likelihood fit with one level, to continuous
#' variables, that used every observed value. A fit with ordered
#' indicators (and so WLSMV or another least squares estimator), with an
#' estimator outside the maximum likelihood family, with listwise
#' deletion (\pkg{lavaan}'s default, which leaves no missing information
#' to measure), with the pairwise, available cases, or two-stage missing
#' data options, with \code{conditional.x = TRUE}, with more than one
#' level, with sampling weights, or that did not converge is refused
#' with a message that says why, as are a \code{lavaanList} or pooled
#' multiple imputation object and a rotated exploratory (\code{efa})
#' solution, whose loadings \pkg{lavaan} rotates after estimation. As a
#' last guard, the computation stops if the refit at the estimates moves
#' them. Besides a rotated solution, this happens with an improper one: a
#' Heywood case whose negative residual variance exceeds the variable's
#' observed variance in size is a solution \pkg{lavaan} (0.7-3) reads as
#' having run away, so it will not take those estimates as the starting
#' point of a refit and starts again from its own simple starting values.
#' The standard errors of such a refit belong to those values, not to the
#' fit, and \code{lavaan::parameterEstimates(fit, fmi = TRUE)}, which makes
#' the same refit, then reports FMIs computed at them without saying so.
#' Every other solution \pkg{lavaan} returns is refit at its estimates,
#' whatever its checks of a solution say of them, and its FMIs are
#' computed there: a Heywood case whose residual variance is less far
#' below zero, two factors correlated above 1 (a covariance matrix of the
#' latent variables that is not positive definite), and a covariance
#' between residuals that implies a correlation above 1. Before it refits,
#' \code{lavaan::parameterEstimates(fit, fmi = TRUE)} shrinks such a
#' covariance to a correlation below 1 (its check of starting values), so
#' it reports the FMIs of other estimates: for two factors correlated
#' 1.08 in a sample of 60, it gave the covariance an FMI of -1.62, where
#' the FMI at the estimates is .25. These are the FMIs of an improper
#' solution, and they deserve the same caution as the solution itself,
#' of which \pkg{lavaan} warns when it fits the model. Finally, the
#' computation stops when the information matrix of a refit cannot be
#' inverted, with a message that names the likeliest cause: variables too
#' badly scaled for a model the refits keep in their own metric (see
#' \emph{Units of measurement}); a predictor (an observed variable
#' regressed on) that lies so far from zero, relative to its standard
#' deviation, that an intercept and its coefficients are nearly collinear,
#' which no division by a power of two changes (a predictor of a
#' regression with a standard deviation of 1 did so from about 8000 from
#' zero, and gave the FMIs of the unshifted data at 7000; an indicator's
#' or an outcome's location enters only its own intercept and does not);
#' a model that is not identified; or a structured or Hessian information
#' that is singular at the estimates. A predictor 1000 or more standard
#' deviations from zero is named only when the same model fit to the data
#' with it centered has an information that can be inverted, since
#' centering, which changes no FMI but those of intercepts, cannot repair
#' a model that is not identified; where centering would change a fixed
#' or constrained intercept, so that this cannot be checked, both causes
#' are given.
#'
#' The computation also stops when the data define no unrestricted
#' estimate of the means and covariances at which to evaluate the
#' fraction (all three estimators need it: the unstructured information
#' is evaluated there, and the complete data side of each is fitted to
#' it). A pair of variables that no case observes together, as in a
#' planned missing design, leaves their covariance undetermined:
#' \pkg{lavaan}'s EM leaves it at its starting value and only warns, and
#' in a two-form design on seven of the Holzinger and Swineford tests
#' the FMIs moved by up to .12 across values of it that the data could not
#' tell apart. The pairs are named, as \code{\link{cov_fiml}} and the
#' reliability family name them. The computation also checks for a
#' supported affine-collapse direction: all rows observing a subset must
#' satisfy the same affine relation, with a nonzero coefficient on every
#' member of the subset. Every partially overlapping row is included when
#' a candidate support is reduced. A subset with no more observed rows
#' than variables and a verified full-support null direction supplies a
#' sufficient unbounded-likelihood diagnosis, even when EM stopped at a
#' stable interior point. The check considers the full variable set,
#' observed-pattern sets, single variables, and reductions of their null
#' supports; it is not an exhaustive existence test.
#'
#' Affine differences are column scaled. A singular-value tolerance of
#' \eqn{100\epsilon\max(n,p)\max(1,d_{\max})}, with machine precision
#' \eqn{\epsilon}, locates candidate nullspaces. A deterministic normalized
#' direction must have every coefficient larger in absolute value than
#' \eqn{p\sqrt{\epsilon}} and a residual below that tolerance. A relation
#' inferred only numerically in an overdetermined subset is reported as
#' numerical instability, without claiming mathematical nonexistence.
#' A near-singular iterate at EM's iteration limit is a fitting failure.
#' No diagnosed obstacle means that no implemented diagnostic failed; it
#' does not certify finite global existence or uniqueness. These checks
#' concern the unrestricted likelihood needed for FMI. The constrained
#' substantive model's own estimates and standard errors do not
#' use the unrestricted estimate (under \pkg{lavaan}'s default
#' information); its chi-square test and the fit indices built on the
#' unrestricted model's likelihood do. These stops are errors of classes
#' \code{dmar_fmi_never_together} and \code{dmar_fmi_saturated_singular},
#' so that a caller can catch them.
#'
#' @return A \code{data.frame} (class \code{dmar_fmi}, \code{dmar_tbl})
#'   with one row per free parameter and per defined (\code{:=})
#'   parameter, in the order of the \pkg{lavaan} parameter table, and the
#'   same columns for every call: \code{term} (\pkg{lavaan}'s
#'   \code{lhs op rhs}, or the name of a defined parameter); \code{op}
#'   (the \pkg{lavaan} operator, for filtering); \code{group} (the group
#'   label; \code{NA} for a single group fit and for defined parameters);
#'   \code{fmi}; \code{fmi_ci_lower} and \code{fmi_ci_upper} (the
#'   percentile limits; \code{NA} when \code{ci_method = "none"});
#'   \code{width_inflation}; \code{N}; \code{N_effective};
#'   \code{N_target}; \code{estimate} (the parameter estimate of
#'   \code{fit}); and \code{se_observed} and \code{se_complete}
#'   (\eqn{\mathrm{SE}_O} and \eqn{\mathrm{SE}_C}, normal theory
#'   standard errors under
#'   \code{information}). The attributes record \code{information} and
#'   \code{information_description}; \code{source_estimator},
#'   \code{source_se}, and \code{source_missing}, the settings of
#'   \code{fit}; \code{basis}; \code{N_used} (per group);
#'   \code{N_excluded}; \code{n_patterns} (per group); \code{fmi_max};
#'   \code{ci_method}; and, with the percentile interval,
#'   \code{conf_level}, \code{B}, \code{B_used}, \code{n_dropped},
#'   \code{B_required}, and \code{bootstrap_status} (see \emph{Bootstrap
#'   replications} under Details). The print method
#'   leaves out the columns that carry nothing for the call at hand
#'   (\code{op}; \code{group} for a single group; the interval limits
#'   without an interval), shows an FMI within \eqn{10^{-10}} of zero (a
#'   fully observed variable's intercept, say) as 0, and states the
#'   information, the sample size, and the imprecision of the estimate
#'   under the table. The object keeps every column at full precision. \code{\link{tidy.dmar_fmi}} and
#'   \code{\link{glance.dmar_fmi}} give the one row per parameter and the
#'   one row summary views.
#'
#' @references
#' Chen, L., & Savalei, V. (2021). Three sample estimates of fraction of
#'   missing information from full information maximum likelihood.
#'   \emph{Frontiers in Psychology, 12}, 667802.
#'   \doi{10.3389/fpsyg.2021.667802}
#'
#' Efron, B., & Tibshirani, R. J. (1993). \emph{An introduction to the
#'   bootstrap}. New York, NY: Chapman & Hall/CRC.
#'
#' Graham, J. W. (2003). Adding missing-data-relevant variables to
#'   FIML-based structural equation models. \emph{Structural Equation
#'   Modeling, 10}(1), 80--100. \doi{10.1207/S15328007SEM1001_4}
#'
#' Maxwell, S. E., Kelley, K., & Rausch, J. R. (2008). Sample size planning
#'   for statistical power and accuracy in parameter estimation.
#'   \emph{Annual Review of Psychology, 59}, 537--563.
#'   \doi{10.1146/annurev.psych.59.103006.093735}
#'
#' Orchard, T., & Woodbury, M. A. (1972). A missing information
#'   principle: Theory and applications. In \emph{Proceedings of the Sixth
#'   Berkeley Symposium on Mathematical Statistics and Probability}
#'   (Vol. 1, pp. 697--715). Berkeley, CA: University of California
#'   Press.
#'
#' Savalei, V., & Rhemtulla, M. (2012). On obtaining estimates of the
#'   fraction of missing information from full information maximum
#'   likelihood. \emph{Structural Equation Modeling, 19}(3), 477--494.
#'   \doi{10.1080/10705511.2012.687669}
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso \code{\link{missing_patterns}} and
#'   \code{\link{missing_coverage}} for the description of the missing
#'   values that comes before any model; \code{\link{fa_confirmatory}} and
#'   \code{\link{mlmr}} for FIML fits; \code{\link{ss_aipe_reliability}}
#'   for planning the sample size of a reliability interval.
#'
#' @family missing data
#'
#' @keywords models multivariate
#'
#' @examples
#' # The six spatial tests of the Holzinger and Swineford data. The revised
#' # paper form board (t25) and the flags test (t26) were given only at
#' # Grant-White, so both are missing for the 156 Pasteur pupils.
#' spatial <- c("t1_visual_perception", "t2_cubes", "t3_paper_form_board",
#'              "t4_lozenges", "t25_paper_form_board_r", "t26_flags")
#' missing_patterns(holzinger_swineford, items = spatial)
#'
#' fit <- fa_confirmatory(data = holzinger_swineford,
#'              factors = list(spatial = spatial),
#'              missing = "fiml", output = "fit")
#' fmi <- fraction_missing_information(fit)
#' fmi
#'
#' # The loadings of t25 and t26 have an FMI near .51, close to their 52
#' # percent of values missing, so each is estimated about as precisely as
#' # 147 complete cases would estimate it. Omega for the six tests loses
#' # less, about a quarter of its information: its interval is 1.15 times
#' # as wide as complete data would give, and the same precision with the
#' # same missing values would take about 398 cases.
#' fmi[fmi$term == "omega_spatial",
#'     c("term", "fmi", "N_effective", "N_target", "width_inflation")]
#'
#' # The Hessian estimator, lavaan's default, for comparison.
#' fraction_missing_information(fit, information = "hessian")$fmi
#'
#' @export
fraction_missing_information <- function(fit,
                                         information = c("unstructured",
                                                         "structured",
                                                         "hessian"),
                                         ci_method = c("none", "percentile"),
                                         conf_level = 0.95,
                                         B = 1000L,
                                         seed = NULL) {
  information <- match.arg(information)
  ci_method <- match.arg(ci_method)
  if (!is.numeric(conf_level) || length(conf_level) != 1L ||
      is.na(conf_level) || conf_level <= 0 || conf_level >= 1) {
    .dmar_refuse_argument(
      "'conf_level' must be a single number in (0, 1).", arg = "conf_level")
  }
  # The rule of ?dmar_bootstrap: 10 replications beyond each limit of the
  # percentile interval (400 at the 95 percent level).
  if (ci_method == "percentile") {
    B <- .dmar_check_B(B, tail = (1 - conf_level) / 2,
                       conf_level = conf_level)
  }
  if (!requireNamespace("lavaan", quietly = TRUE)) {
    .dmar_refuse(c(
      "Package 'lavaan' is required. Install it with ",
      "install.packages(\"lavaan\")."), "package_required", package = "lavaan")
  }
  display_names <- .fmi_mlmr_names(fit)
  fit <- .fmi_extract_fit(fit)
  opt <- .fmi_check_fit(fit)

  pt <- .fmi_with_names(.fmi_core(fit, information), display_names)
  # The Hessian is differentiated numerically with a step of at least
  # 1e-6 (lavaan's), too coarse for the parameters of variables measured
  # in small units and swamped by rounding for large ones. The refits are
  # made with the variables' standard deviations near 1, where the step
  # suits them, unless the model cannot be mapped to that metric
  # (attr(pt, "scale_reason")); in the variables' own metric the Hessian
  # FMIs of one regression were off by .004 at a standard deviation of
  # 10,000 and by .5 at 100,000, and were partly missing at .001, while
  # from .01 to 1000 they agreed with the analytic estimators to 4e-5.
  sd_range <- attr(pt, "sd_range")
  if (identical(information, "hessian") && all(is.finite(sd_range)) &&
      (sd_range[1L] < 0.01 || sd_range[2L] > 1000)) {
    reason <- attr(pt, "scale_reason")
    .dmar_warn(c(
      "The variables' standard deviations run from ",
      format(sd_range[1L], digits = 3L), " to ",
      format(sd_range[2L], digits = 3L),
      if (!is.null(reason)) {
        paste0(", and the refits behind the FMI could not be moved to ",
               "a metric in which they are near 1 because the model ",
               "has ", reason)
      },
      ". The Hessian estimator differentiates the log likelihood ",
      "numerically with a step suited to standard deviations between ",
      "about .01 and 1000, and outside that range its FMIs can be far ",
      "off. Divide the variables by constants that bring their ",
      "standard deviations nearer 1, which leaves the FMI unchanged, ",
      "or use information = \"unstructured\" (the default) or ",
      "\"structured\", which are computed analytically."), "result_qualified")
  }
  keep <- which(pt$free > 0L | pt$op == ":=")
  nobs <- as.numeric(unlist(lavaan::lavInspect(fit, "nobs")))
  norig <- as.numeric(unlist(lavaan::lavInspect(fit, "norig")))
  ng <- length(nobs)
  group_labels <- if (ng > 1L) {
    as.character(lavaan::lavInspect(fit, "group.label"))
  } else {
    character(0)
  }
  row_N <- .fmi_row_N(pt, nobs)
  group <- rep(NA_character_, length(pt$lhs))
  if (ng > 1L) {
    in_group <- pt$group >= 1L & pt$op != ":="
    group[in_group] <- group_labels[pt$group[in_group]]
  }
  prec <- .fmi_precision_terms(pt$fmi, row_N)

  out <- data.frame(
    term            = .fmi_term(pt, display_names)[keep],
    op              = pt$op[keep],
    group           = group[keep],
    fmi             = pt$fmi[keep],
    fmi_ci_lower    = NA_real_,
    fmi_ci_upper    = NA_real_,
    width_inflation = prec$width_inflation[keep],
    N               = row_N[keep],
    N_effective     = prec$N_effective[keep],
    N_target        = prec$N_target[keep],
    estimate        = pt$est[keep],
    se_observed     = pt$se_observed[keep],
    se_complete     = pt$se_complete[keep],
    stringsAsFactors = FALSE,
    row.names = NULL
  )

  B_used <- NULL
  if (ci_method == "percentile") {
    .dmar_local_seed(seed)
    boot_fun <- function(f) .fmi_core(f, information)$fmi[keep]
    # Inside the bootstrap every replication can warn: lavaan notes the
    # cases with nothing observed, its post estimation checks flag an
    # improper solution in a resample, and lavBootstrap() reports the
    # replications that failed. Those are noise here; the replications
    # that did not return an FMI for every parameter are counted below
    # and reported in one warning.
    boots <- try(suppressWarnings(
      lavaan::lavBootstrap(fit, r = as.integer(B), fun = boot_fun)),
      silent = TRUE)
    if (inherits(boots, "try-error")) {
      .dmar_refuse(c(
        "The bootstrap failed: ",
        conditionMessage(attr(boots, "condition"))), "computation_failed")
    }
    boots <- matrix(as.numeric(boots), ncol = length(keep))
    ok <- rowSums(!is.finite(boots)) == 0L
    B_used <- sum(ok)
    # The verdict of the rule (R/dmar_bootstrap.R). Too few usable
    # replications leave the limits NA under a classed warning; through
    # the first versions of 1.1.0 the call stopped below two.
    verdict <- .dmar_boot_settle(
      B, B_used, tail = (1 - conf_level) / 2, conf_level = conf_level,
      fun = "fraction_missing_information",
      quantity = "the FMI of some parameter",
      dropped_because = paste0("failed to refit, did not converge, or ",
                               "returned an FMI that is not a finite ",
                               "number"))
    boots <- boots[ok, , drop = FALSE]
    alpha_2 <- (1 - conf_level) / 2
    if (verdict$formed) {
      out$fmi_ci_lower <- apply(boots, 2L, stats::quantile, probs = alpha_2,
                                names = FALSE)
      out$fmi_ci_upper <- apply(boots, 2L, stats::quantile,
                                probs = 1 - alpha_2, names = FALSE)
    } else {
      out$fmi_ci_lower <- out$fmi_ci_upper <- NA_real_
    }
  }

  patterns <- lavaan::lavInspect(fit, "patterns")
  if (ng == 1L) patterns <- list(patterns)
  per_group <- function(v) {
    v <- as.numeric(v)
    if (ng > 1L) names(v) <- group_labels
    v
  }

  out <- .as_dmar_tbl(out,
                      conf_level = if (ci_method == "percentile") conf_level,
                      subclass = "dmar_fmi")
  attr(out, "information") <- information
  attr(out, "information_description") <-
    .fmi_information_description(information)
  attr(out, "source_estimator") <- opt$estimator.orig %||% opt$estimator
  attr(out, "source_se") <- as.character(opt$se)
  attr(out, "source_missing") <- as.character(opt$missing)
  attr(out, "basis") <- paste0("normal theory information at the ",
                               "estimates of the source fit")
  attr(out, "N_used") <- per_group(nobs)
  attr(out, "N_excluded") <- sum(norig - nobs)
  attr(out, "n_patterns") <- per_group(vapply(patterns, NROW, integer(1L)))
  attr(out, "fmi_max") <- if (all(is.na(out$fmi))) NA_real_ else
    max(out$fmi, na.rm = TRUE)
  attr(out, "ci_method") <- ci_method
  if (ci_method == "percentile") out <- .dmar_boot_attrs(out, verdict)
  out
}


# Every display of the table (the console, knit_print(), as_kable(); see
# .dmar_tbl_display() and .dmar_tbl_footer() in R/dmar_output_helpers.R)
# leaves out the columns that carry nothing for the call at hand, shows an
# FMI within 1e-10 of zero (a fully observed variable's intercept, say,
# whose ratio differs from zero only by rounding error) as 0, and states
# beneath the table how the numbers were obtained. The object itself is
# not changed. A subset of the table, which no longer carries the
# attributes, displays as an ordinary DMAR table.
#' @exportS3Method .dmar_tbl_display dmar_fmi
#' @noRd
.dmar_tbl_display.dmar_fmi <- function(x, ...) {
  if (is.null(attr(x, "information"))) return(NextMethod())
  drop <- "op"
  if ("group" %in% names(x) && all(is.na(x$group))) drop <- c(drop, "group")
  if (all(c("fmi_ci_lower", "fmi_ci_upper") %in% names(x)) &&
      all(is.na(x$fmi_ci_lower)) && all(is.na(x$fmi_ci_upper))) {
    drop <- c(drop, "fmi_ci_lower", "fmi_ci_upper")
  }
  disp <- x
  class(disp) <- "data.frame"
  disp <- disp[setdiff(names(disp), drop)]
  for (nm in intersect(c("fmi", "fmi_ci_lower", "fmi_ci_upper"),
                       names(disp))) {
    v <- disp[[nm]]
    v[!is.na(v) & abs(v) < 1e-10] <- 0
    disp[[nm]] <- v
  }
  .as_dmar_tbl(disp)
}

#' @exportS3Method .dmar_tbl_footer dmar_fmi
#' @noRd
.dmar_tbl_footer.dmar_fmi <- function(x, ...) {
  .dmar_footer_add(NextMethod(), .fmi_footer(x), "fmi")
}

# The lines written beneath a fraction_missing_information() table: the
# information and the cases the FMI came from, cases the fit set aside,
# the source fit's estimator, negative estimates, and the bootstrap.
# NULL for a subset that no longer carries the attributes.
.fmi_footer <- function(x) {
  if (is.null(attr(x, "information"))) return(NULL)
  N_used <- attr(x, "N_used")
  n_text <- paste(format(sum(N_used), trim = TRUE), "cases")
  if (length(N_used) > 1L) {
    n_text <- sprintf("%s (%s)", n_text,
                      paste(sprintf("%s in %s", format(N_used, trim = TRUE),
                                    names(N_used)), collapse = ", "))
  }
  estimator_ordinal <- c(hessian = "first", structured = "second",
                         unstructured = "third")
  lines <- sprintf(paste0("FMI under the %s information (the %s ",
                          "estimator of Chen & Savalei, 2021), from %s."),
                   attr(x, "information"),
                   estimator_ordinal[[attr(x, "information")]], n_text)
  n_excluded <- attr(x, "N_excluded")
  if (isTRUE(n_excluded > 0)) {
    lines <- c(lines, sprintf(paste0("The fit set aside %d case%s with no ",
                                     "usable observed value."),
                              as.integer(n_excluded),
                              if (n_excluded == 1) "" else "s"))
  }
  source_se <- attr(x, "source_se")
  source_estimator <- attr(x, "source_estimator")
  if (!identical(source_se, "standard") ||
      !identical(source_estimator, "ML")) {
    lines <- c(lines, sprintf(paste0("The source fit is %s with %s standard ",
                                     "errors; the FMI is the normal theory ",
                                     "FMI at its estimates."),
                              source_estimator, source_se))
  }
  if (any(x$fmi < 0, na.rm = TRUE) &&
      !identical(attr(x, "information"), "unstructured")) {
    lines <- c(lines, paste0("Some estimates are negative, which the ",
                             "structured and Hessian estimators can give in ",
                             "small samples or under misspecification; the ",
                             "unstructured estimator cannot."))
  }
  if (identical(attr(x, "ci_method"), "percentile")) {
    status <- attr(x, "bootstrap_status")
    required <- attr(x, "B_required")
    lines <- c(lines, paste0(
      sprintf(paste0("Bootstrap: %d replications attempted, %d returned ",
                     "an FMI for every parameter"),
              attr(x, "B"), attr(x, "B_used")),
      if (identical(status, "insufficient_share")) {
        "; more than half were dropped, so no interval is reported."
      } else if (identical(status, "insufficient_count")) {
        sprintf(paste0("; fewer than the %d required at this level were ",
                       "requested, so no interval is reported."), required)
      } else if (identical(status, "exploratory")) {
        sprintf(paste0("; at least %d are required for a reported ",
                       "interval. Exploratory: below the number required ",
                       "for a reported interval."), required)
      } else {
        "."
      }))
    lines <- c(lines, paste0("The sample FMI is imprecise unless N is ",
                             "several hundred (Chen & Savalei, 2021)."))
  } else {
    lines <- c(lines, paste0("The sample FMI is imprecise unless N is ",
                             "several hundred (Chen & Savalei, 2021); ",
                             "ci_method = \"percentile\" gives an ",
                             "interval."))
  }
  lines
}


#' Tidy and Glance Methods for the Fraction of Missing Information
#'
#' \code{tidy()} returns one row per parameter of a
#' \code{\link{fraction_missing_information}} table, with the FMI as the
#' estimate, its percentile limits, and the effective sample size;
#' \code{glance()} returns a one-row summary of the table.
#'
#' @param x A \code{dmar_fmi} object returned by
#'   \code{\link{fraction_missing_information}}.
#' @param \dots Unused.
#'
#' @return \code{tidy()} returns a \code{data.frame} with columns
#'   \code{term} (the parameter), \code{group}, \code{estimate} (its
#'   FMI), \code{estimand} (\code{"fmi"}, which names the estimate, so
#'   that \code{\link{results_sentence}} on a row reports
#'   \code{"fmi = 0.24, 90\% CI [0.18, 0.32]"} and not the FMI as the
#'   parameter's value), \code{ci_lower}, \code{ci_upper},
#'   \code{conf_level} (all three \code{NA} without the percentile
#'   interval), and \code{N_effective}. \code{glance()} returns a one-row
#'   \code{data.frame} with \code{n_parameters}, \code{fmi_max},
#'   \code{fmi_median}, \code{N} (the cases the fit used, summed over
#'   groups), and \code{information}.
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @examples
#' fit <- fa_confirmatory(data = holzinger_swineford,
#'              factors = list(spatial = c("t1_visual_perception",
#'                "t2_cubes", "t3_paper_form_board", "t4_lozenges",
#'                "t25_paper_form_board_r", "t26_flags")),
#'              missing = "fiml", output = "fit")
#' fmi <- fraction_missing_information(fit)
#' generics::tidy(fmi)
#' generics::glance(fmi)
#'
#' @importFrom generics tidy
#' @exportS3Method generics::tidy
tidy.dmar_fmi <- function(x, ...) {
  # A table whose columns were subset keeps its class; a column it no
  # longer has is returned as NA so the shape stays the same.
  column <- function(nm, na) if (nm %in% names(x)) x[[nm]] else
    rep(na, nrow(x))
  data.frame(
    term        = as.character(column("term", NA_character_)),
    group       = as.character(column("group", NA_character_)),
    estimate    = as.numeric(column("fmi", NA_real_)),
    # The estimate is the parameter's FMI, not the parameter: the estimand
    # column says so, and results_sentence() labels the row "fmi" rather
    # than with the term, which names the parameter (through 1.1.0's
    # development a row read "omega_spatial = 0.24" for an omega of .76).
    estimand    = rep("fmi", nrow(x)),
    ci_lower    = as.numeric(column("fmi_ci_lower", NA_real_)),
    ci_upper    = as.numeric(column("fmi_ci_upper", NA_real_)),
    # The level of the percentile intervals (NA when none was formed), so
    # that a tidy row states the level of its limits.
    conf_level  = rep(attr(x, "conf_level", exact = TRUE) %||% NA_real_,
                      nrow(x)),
    N_effective = as.numeric(column("N_effective", NA_real_)),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
}

#' @rdname tidy.dmar_fmi
#' @importFrom generics glance
#' @exportS3Method generics::glance
glance.dmar_fmi <- function(x, ...) {
  fmi <- if ("fmi" %in% names(x)) as.numeric(x[["fmi"]]) else NA_real_
  has <- any(!is.na(fmi))
  N_used <- attr(x, "N_used")
  data.frame(
    n_parameters = nrow(x),
    fmi_max      = if (has) max(fmi, na.rm = TRUE) else NA_real_,
    fmi_median   = if (has) stats::median(fmi, na.rm = TRUE) else NA_real_,
    N            = if (is.null(N_used)) NA_real_ else sum(N_used),
    information  = attr(x, "information") %||% NA_character_,
    stringsAsFactors = FALSE
  )
}
