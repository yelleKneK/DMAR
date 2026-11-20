#' Bootstrap Confidence Intervals for the Bifactor Indices
#'
#' Forms nonparametric bootstrap confidence intervals for the indices
#' \code{\link{bifactor_indices}} reports for a fitted bifactor model: the
#' explained common variance (ECV), coefficient omega, omega hierarchical
#' (omega_H) for the general factor, omega hierarchical subscale
#' (omega_HS) for each group factor, the index Rodriguez, Reise, and
#' Haviland (2016) call coefficient H, the closed form of Hancock and
#' Mueller (2001) applied to one factor's loadings, which counts the other
#' factors' common variance as uncorrelated error, and the maximal
#' reliability of each factor under the fitted model, which the index is
#' not (see \code{\link{bifactor_indices}}). Each replication resamples
#' the cases the fit used, refits the same model with the same estimator,
#' and recomputes every index, so the intervals rest on the same
#' replications (those of the maximal reliability, under
#' \code{improper = "keep"}, on the proper ones among them). Replications
#' whose solution is improper are set aside by default and counted, and
#' the table says how many.
#'
#' @param fit A fitted bifactor \pkg{lavaan} model to raw data, as for
#'   \code{\link{bifactor_indices}}: one general factor on all items plus
#'   two or more orthogonal group factors, each item on one group factor at
#'   most (an item on none is on the general factor alone). A fit that
#'   \code{bifactor_indices()} refuses is refused here, and so are a
#'   multiple-group fit, a fit with a cluster
#'   variable or sampling weights, a fit to a covariance matrix, and a fit
#'   that did not converge (see \emph{Details}).
#' @param general Optional name of the general factor, as for
#'   \code{\link{bifactor_indices}}; by default the factor that loads on
#'   every item.
#' @param denominator How the total-score variance in the denominators of
#'   omega, omega_H, and omega_HS is estimated, as for
#'   \code{\link{bifactor_indices}}: \code{"model_implied"} (the default)
#'   or \code{"observed"}. Each replication uses the same setting.
#' @param ci_method The bootstrap interval: \code{"percentile"} (the
#'   default) or \code{"bca"}, the bias-corrected and accelerated
#'   interval. See \emph{Details}.
#' @param conf_level Confidence level of the intervals. Defaults to
#'   \code{0.95}.
#' @param B Number of bootstrap replications, a single whole number.
#'   Defaults to \code{1000}, the package's count for a bootstrap whose
#'   every replication refits a model; see \emph{Details} for the cost. The
#'   fewest accepted at the default \code{conf_level} is 400 for
#'   \code{"percentile"} (10 replications beyond each limit; 200 at 0.90,
#'   and 2000 at 0.99, which the default does not meet) and 2000 for
#'   \code{"bca"} (50 beyond each nominal limit), which the default does
#'   not meet either: pass \code{B = 2000} or more with \code{ci_method =
#'   "bca"}. A smaller \code{B} is refused before anything is computed
#'   (class \code{dmar_count_invalid}). The count is held to the number
#'   requested: replications set aside afterward are counted and reported in
#'   a warning of class \code{dmar_bootstrap_replications_dropped}, and an
#'   index for which fewer than half of the \code{B} replications return a
#'   value has its interval withheld (\code{NA}) under a warning of class
#'   \code{dmar_bootstrap_insufficient}. Where refits can fail, request more
#'   than the minimum. See \code{\link{dmar_bootstrap}},
#'   which also describes \code{options(dmar.bootstrap_exploratory = TRUE)}
#'   for a first look with a smaller \code{B}.
#' @param seed Random number seed for the bootstrap. The default
#'   \code{NULL} uses the current state of the random number generator;
#'   a supplied seed is set internally and the caller's state is restored
#'   on exit.
#' @inheritParams reliability_omega
#' @param ncpus The number of processes over which the replications are
#'   spread. The default \code{1} runs them in this R session; a larger
#'   number uses forked workers (\pkg{boot}'s
#'   \code{parallel = "multicore"}) or, on Windows, a socket cluster
#'   (\code{"snow"}), whose workers load the installed \pkg{DMAR}.
#'   \pkg{boot} draws every resample before any replication runs, so a
#'   given \code{seed} returns the same intervals, counts, and replicates
#'   whatever \code{ncpus} is. A session that forks workers and also
#'   starts background processes through \pkg{callr} or \pkg{processx}
#'   needs the environment variable \code{PROCESSX_NOTIFY_OLD_SIGCHLD} set
#'   before \pkg{processx} is loaded (see \code{\link{sem_exact_fit_bootstrap}}).
#' @param improper How a replication whose solution is improper (a
#'   standardized loading beyond one or a negative residual variance, a
#'   Heywood case) is treated: \code{"drop"} (the default) sets it aside,
#'   and \code{"keep"} uses its indices when every one of them but the
#'   maximal reliability, which is not defined for an improper solution,
#'   is a proportion. Either way such replications are counted. With
#'   \code{"drop"} an improper point estimate stops the call. See
#'   \emph{Details}.
#'
#' @details
#' \strong{Why intervals.} The bifactor indices are estimates, and the
#' decisions they are read for (whether a total score is close enough to
#' unidimensional to report alone, whether a subscale carries reliable
#' variance beyond the general factor) are decisions about the population
#' values. A confidence interval says how precisely the sample has
#' located each of them. In the bifactor model of 15 Holzinger and
#' Swineford tests described under \emph{Improper replications} below,
#' omega_H is 0.577 with a 95 percent interval from 0.468 to 0.672
#' (\code{B = 1000}, \code{seed = 113}): 301 cases locate the share of
#' total score variance due to the general factor only to within an
#' interval about 0.2 wide. \code{bifactor_indices()} reports the point
#' estimates alone. The intervals are a separate function, rather than
#' an argument of \code{bifactor_indices()}, because they need columns
#' that its one-row-per-factor table cannot gain without its shape
#' depending on an argument.
#'
#' \strong{What is resampled.} The cases are the rows of data
#' \pkg{lavaan} stored with \code{fit}, which are the rows the fit used:
#' under listwise deletion (the default, and \pkg{lavaan}'s default for
#' ordered items) the complete rows, and under \code{missing = "ml"}
#' (full information maximum likelihood) every row with at least one
#' observed item, each with its own pattern of observed and missing
#' values, so that a resample has the missing data the sample has and
#' each replication is again a full information fit. A row with nothing
#' observed, which \pkg{lavaan} keeps in its stored data but leaves out
#' of the sample, is left out here too, and the number of rows resampled
#' is checked against the number of cases \pkg{lavaan} reports. Rows are
#' drawn with replacement, \code{B} times (Efron & Tibshirani, 1993).
#'
#' \strong{How each replication is fit.} Each resample is fit with the
#' parameter table and the options of \code{fit}: the same model, the
#' same estimator, the same missing data treatment, the same ordered
#' items, and default starting values, as for \code{fit} itself. Standard
#' errors and test statistics are not computed in a replication, since no
#' index uses them; robust maximum likelihood (MLR) and WLSMV have the
#' estimates of maximum likelihood and diagonally weighted least squares
#' and differ from them only in those, so a replication is the cheaper
#' fit with the same estimates. Before any resample is drawn, the model is
#' refit this way to the rows of \code{fit} itself, and the call stops
#' unless that refit reproduces the point estimates, so the replications
#' are refits of the model the point estimates describe. The indices of
#' each replication are \code{bifactor_indices()} of its refit, with the
#' general factor named, so every index is computed as the point estimate
#' is.
#'
#' \strong{Improper replications.} A bifactor model is prone to improper
#' solutions, a standardized loading past one or a residual variance
#' below zero, which can come from sampling error around a small residual
#' variance, from a model that does not describe the items, or from a
#' group factor the data barely separate from the general factor, and a
#' resample can give such a solution where the sample does not. A
#' replication is improper when \code{bifactor_indices()} flags its
#' solution or \pkg{lavaan}'s own check after estimation fails, and the
#' point estimate is judged the same way. For a fit with auxiliary
#' variables (\pkg{lavaan}'s \code{aux} under full information maximum
#' likelihood) that check is made on the model's own parameters, as
#' \code{\link{sem_fit_measures}} explains: \pkg{lavaan}'s own check also
#' tests the residual covariances of the auxiliary variables, which are
#' not positive definite once one of them is related to the general
#' factor strongly enough, however proper the bifactor solution. A
#' bifactor model of 15 of the Holzinger and Swineford tests
#' (\code{holzinger_swineford}: a general factor on all of them, with
#' spatial, verbal, and memory group factors on tests 1--4, 5--9, and
#' 14--19) gave an improper solution in 212 of 1,000
#' resamples (\code{seed = 113}) and did not converge in 17 more. The
#' omega of the improper solutions ranged from 0.84 to 26.5, 132 of them
#' above one; kept, all of them put the upper limit of the 95 percent
#' interval for omega at 13.7 and the lower limit for ECV at 0.004, where
#' the proper replications give [0.855, 0.902] and [0.309, 0.509]. Such
#' values are not sampling variability around a proper structure: a
#' solution in which a residual variance has run away, or in which the
#' general factor has all but vanished (an ECV near zero), is not a
#' perturbation of the sample's solution but a different one, and its
#' indices do not estimate the quantities the point estimate does. The
#' default, \code{improper = "drop"}, therefore sets them aside. This
#' differs from the default of
#' the reliability family (\code{\link{reliability_omega}}), where an
#' improper single-factor replication that still yields a coefficient is
#' typically a loading that ran just past one, and keeping it preserves
#' the upper tail of the bootstrap distribution. Setting replications
#' aside conditions the intervals on admissible solutions: they describe
#' the sampling distribution of the indices among resamples that the
#' model fits properly. Yuan and Hayashi (2003, Section 6) argue that a
#' replication that fails because its resample lies far from the model
#' is informative and should be counted rather than discarded; for a test
#' statistic they count it as an extreme value, but an index has no side
#' on which such a replication belongs, so here it is counted and
#' reported rather than placed. Every replication set aside, for any
#' reason, is counted in the \code{bootstrap} attribute, named in a
#' warning, and stated under the printed table, and when the replications
#' set aside (or the improper ones kept) reach the tail area of an
#' interval, \eqn{(1 - C)/2} of \code{B} for a confidence level
#' \eqn{C}, the warning and the printed line add a caution that the
#' limits depend on how they were treated. When fewer than half of the
#' \code{B} replications, or fewer than the number the interval method
#' and the level require (400 for the percentile interval and 2000 for
#' BCa at the 95 percent level; see \code{\link{dmar_bootstrap}}),
#' return a value for an index, its interval is withheld, and the warning
#' (class \code{dmar_bootstrap_insufficient}) and the printed line say
#' so, as in \code{\link{reliability_omega}}. A run at exactly the
#' required \code{B} therefore loses its intervals to a single
#' replication set aside; the default \code{B = 1000} leaves room for
#' them at the 95 percent level. With
#' \code{improper = "keep"} an improper replication is used when every index it gives but the
#' maximal reliability, which is not defined for an improper solution
#' (see \emph{The maximal reliability} below), is a proportion; one
#' whose solution gives an index outside [0, 1] (an omega above one
#' belongs to a composite whose implied error variance is negative) or
#' an H that is not defined (a standardized loading at or past one leaves
#' its item no error variance in H's formula) yields no coefficient and
#' is set aside under either setting, as in the reliability family. So is
#' a replication whose likelihood has no maximum, a residual variance
#' running off below zero while the likelihood is flat or still rising
#' (see \code{\link{fa_bifactor}}): it is counted among the replications that
#' raised an error, under a message that names the item, and the run goes
#' on (of the 22 improper replications in 100 of the 15-test model above,
#' \code{seed = 113}, 8 were of this kind), while a refit of the data
#' stored with \code{fit} that has no maximum stops the call with an
#' error of class \code{dmar_no_maximum}.
#' Keeping the admissible indices of improper replications can change
#' their interval limits; the maximal reliability uses only proper
#' replications under either setting. With
#' \code{"drop"} an improper point estimate stops the call, since an
#' interval built from proper replications would not describe it; with
#' \code{"keep"} it is reported with a warning, unless one of its indices
#' other than the maximal reliability is not defined or lies outside
#' [0, 1]; its maximal reliability, which is not defined, is reported as
#' \code{NA}.
#'
#' \strong{The maximal reliability.} The rows \code{maximal_reliability}
#' give each factor's maximal reliability under the fitted model,
#' \eqn{\lambda_f' R^{-1} \lambda_f} with \eqn{R} the item correlations
#' the model reproduces (the \code{"maximal_reliability"} attribute of
#' \code{bifactor_indices()}): the reliability the optimally weighted
#' items reach as a measure of the factor, and the square of the
#' determinacy of its regression factor score. Like every other index, it
#' is read from \code{bifactor_indices()} on each refit, so its interval
#' costs nothing beyond the replications the other indices need. It is a
#' reliability estimate in its own right, as much in need of an interval
#' as omega (Raykov & Penev, 2006, give one for a composite of congeneric
#' measures), and it is the quantity to read beside the index H, which
#' can lie above or below it. It is not defined for an improper solution.
#' With \code{improper = "drop"} every replication used is proper, and
#' its intervals rest on the same replications as the others; with
#' \code{"keep"} the improper replications kept for the other indices
#' have no maximal reliability, so its intervals rest on the proper
#' replications among them, which are the replications \code{"drop"}
#' uses, and the printed table says how many. Its interval is therefore
#' the same under either rule. An improper point estimate, reported under
#' \code{"keep"}, has no maximal reliability, and these rows then carry
#' no estimate and no interval.
#'
#' \strong{The bias of H and the maximal reliability.} The sample H, the
#' closed form the index applies to a factor's loadings, is a positively
#' biased estimator of its population value, and a composite weighted with
#' sample estimates is on average less reliable than the population value:
#' Aguirre-Urreta, Rönkkö, and McIntosh (2019) found both effects
#' pronounced in samples of about 100 and slight beyond about 1,000, for
#' H of a single factor. The maximal reliability is the reliability of
#' weights that are optimal for the sample's estimates and is biased
#' upward in the same way, most of all for a group factor that the data
#' determine weakly. Run with \code{B = 200} in place of its 10, the
#' example below gives a mean of the replicates that exceeds the estimate
#' by 0.031 for the first group factor's maximal reliability (whose
#' replicates have a standard deviation of 0.078) and by 0.019 for its H,
#' against 0.005 for its omega and 0.003 for its omega_HS. The percentile
#' interval takes its
#' limits from the replicates as they fall, and so carries the bias; the
#' BCa interval (\code{ci_method = "bca"}) moves its limits to correct for
#' the median bias of the replicates. The bias of each index with an
#' interval, the mean of the replicates its interval uses minus the
#' estimate, is the \code{bias} element of the \code{bootstrap} attribute
#' (PUC, which has no interval, has none, and neither has the maximal
#' reliability of an improper point estimate under \code{improper =
#' "keep"}), and the printed table gives it for H and the maximal
#' reliability. Read the maximal reliability as
#' a ceiling that weights estimated from these data will not quite reach.
#'
#' \strong{Other replications set aside.} A refit that does not converge
#' is set aside and counted as such, and so is one that raises an error;
#' the warning quotes the most frequent error. For ordered items, a
#' resample that draws no response in some category of an item cannot
#' refit that item's thresholds, and is set aside as an error that says
#' so. The intervals are formed from the replications that remain, and
#' the effective count is the \code{B_used} attribute.
#'
#' \strong{The two intervals.} The percentile interval takes its limits
#' from the empirical quantiles of the bootstrap estimates; it respects
#' the [0, 1] range of the indices and is not forced to be symmetric.
#' Kelley and Pornprasertmanit (2016) found it the most consistently
#' accurate of the intervals they compared for the single-factor omega
#' with an observed denominator (their hierarchical omega, which is not
#' the omega_H here), which is the reason it is the default; the sources
#' cited here do not evaluate the coverage of bootstrap intervals for
#' the bifactor indices themselves. The bias-corrected and accelerated
#' (BCa) interval moves the two quantile positions to correct for median
#' bias and for the rate at which the variance of the estimator changes
#' with the parameter, which makes it second-order accurate where the
#' percentile interval is first-order accurate (DiCiccio & Efron, 1996).
#' Its acceleration is estimated as \pkg{boot} estimates it, by
#' regression on the resampling frequencies when at least as many
#' replications as cases were used and otherwise, or when that regression
#' fails, by the jackknife, one leave-one-out refit per case, which is run
#' once for all of the indices; an index whose BCa interval cannot be
#' formed gets the percentile interval, with a warning, and the
#' \code{ci_method} attribute and the printed line say which. PUC, the
#' proportion of item pairs whose items lie on different group factors,
#' depends only on the model's structure, so it is the same in every
#' replication and has no interval; its row carries the point value with
#' \code{NA} limits.
#'
#' \strong{Cost.} Each replication is a model fit, so the call costs
#' \code{B} fits and, for \code{"bca"} with fewer replications used than
#' cases, as many more as there are cases. A proper fit of the 15-test
#' Holzinger and Swineford model takes a few hundredths of a second, but
#' a solution that runs away took seconds before the optimizer stopped,
#' so a model that gives many improper replications is slow to bootstrap
#' as well as hard to interpret. \code{ncpus} spreads the replications
#' over several processes, and \code{seed} makes a call reproducible.
#' Requires the \pkg{lavaan} and \pkg{boot} packages.
#'
#' \strong{What is refused.} A fit that \code{bifactor_indices()} refuses
#' (correlated factors, no factor that loads on every item, fewer than
#' two group factors, an item on more than one) is refused, with the same
#' message. So is
#' a multiple-group fit, since resampling the pooled rows would let the
#' group sizes vary from replication to replication (fit each group
#' separately); a fit with a \code{cluster} variable, whose cases are not
#' independent and cannot be resampled one at a time; a fit with
#' sampling weights, which a resample of rows does not carry; a fit to a
#' covariance matrix, which has no cases to resample; and a fit that did
#' not converge.
#'
#' Bootstrap intervals from supplied two-stage fits currently refuse
#' separately stored auxiliary variables: each refit must resample and
#' carry those values with the items. Use an eligible delta interval
#' where available, or explicitly resample and refit them together.
#'
#' @return A \code{data.frame} (class \code{dmar_ci_bifactor},
#'   \code{dmar_tbl}) with one row for each index that applies to a
#'   factor, the general factor first (ECV, omega, omega_H, PUC, H,
#'   maximal_reliability) and then each group factor (ECV, omega,
#'   omega_HS, H, maximal_reliability; its ECV is the general factor's
#'   share of the common variance of that factor's items, and its omega
#'   the reliability of their total score, as defined under
#'   \emph{Details} of \code{\link{bifactor_indices}}, which the printed
#'   table repeats), and columns
#'   \code{factor}, \code{index}, \code{estimate} (the value
#'   \code{bifactor_indices()} reports, in its table or, for
#'   \code{maximal_reliability}, in its \code{"maximal_reliability"}
#'   attribute), \code{se} (the standard deviation of the replications
#'   used), \code{ci_lower}, and \code{ci_upper}. PUC has \code{NA}
#'   in \code{se} and the limits, and so does every column of the
#'   maximal_reliability rows when the point estimate is improper, and
#'   so do the \code{se} and limits of any index whose interval is
#'   withheld (see \code{require_half} under the \code{bootstrap}
#'   attribute below).
#'   The rows depend on the model alone, and the columns never change.
#'   Attributes:
#'   \describe{
#'     \item{\code{bootstrap}}{The counts and draws: \code{attempted}
#'       (\code{B}), \code{returned} (the replications used),
#'       \code{failed} (those set aside), \code{reasons},
#'       \code{error_messages}, \code{undetermined} (always zero here: it
#'       counts the replications that left a polychoric correlation
#'       undetermined, which only \code{\link{reliability_omega_categorical}}
#'       meets), \code{improper} (improper replications set
#'       aside), \code{improper_kept}, \code{no_solution} (always zero
#'       here: it counts the resamples of three items whose correlations
#'       admit no single-factor solution, which only the single-factor
#'       coefficients of \code{\link{reliability_omega}} and
#'       \code{\link{reliability_H}} meet),
#'       \code{not_converged}, \code{too_few_cases} (always zero here:
#'       it counts the replications the single-factor coefficients set
#'       aside for too few distinct cases),
#'       \code{improper_rule}, \code{caution} (\code{TRUE} when the
#'       replications set aside and the improper ones kept number at
#'       least \eqn{(1 - C)/2} of \code{B}, with \eqn{C} the confidence
#'       level, so that how they were treated can decide where a limit
#'       falls), \code{require_half} (\code{TRUE}: an index's interval
#'       requires at least half of all attempted replications to return a
#'       value for it; otherwise its \code{se},
#'       \code{ci_lower}, and \code{ci_upper} are \code{NA}, and
#'       the warning and the printed line say so; this is an
#'       estimability safeguard, not a coverage guarantee),
#'       \code{required} (the fewest replications a request
#'       may ask for, 400 for the percentile interval and 2000 for BCa at
#'       the 95 percent level; the count is held to \code{B}, not to the
#'       replications that return a value) and \code{bootstrap_status}
#'       (\code{"complete"}, \code{"dropped"}, \code{"exploratory"},
#'       \code{"insufficient_count"}, or \code{"insufficient_share"},
#'       the last two when some index has its interval withheld; see
#'       \code{\link{dmar_bootstrap}}), \code{target_returned} and
#'       \code{target_withheld} (for each bootstrapped index, that is,
#'       each column of \code{t} and named as those columns, the
#'       replications that returned a value for it and whether its
#'       interval was withheld),
#'       \code{acceleration} and
#'       \code{acceleration_method} (one per index), \code{jackknife},
#'       \code{ncpus}, \code{t} (a \code{B} by index matrix of the
#'       estimates of every replication, one column for each row with an
#'       interval, named \code{"factor:index"}, and \code{NA} where a
#'       replication gave none; the indices of an improper replication
#'       set aside under \code{"drop"} are kept here when they are
#'       proportions, and its maximal reliabilities are \code{NA}),
#'       \code{status} (the fate of every replication, a factor with
#'       levels \code{"proper"}, \code{"improper"}, \code{"error"},
#'       \code{"non_finite"}, \code{"undetermined"}, and
#'       \code{"not_converged"}), and \code{bias} (one per column of
#'       \code{t}, named as its columns: the mean of the replicates the
#'       index's interval uses minus its estimate, the bootstrap estimate
#'       of its bias; see \emph{Details}).}
#'     \item{\code{B}, \code{B_used}}{The replications requested and the
#'       smallest number any reported interval rests on (under
#'       \code{improper = "keep"}, the intervals of the maximal
#'       reliability rest on the proper replications, the count of finite
#'       values in a maximal reliability column of \code{t}, which can be
#'       fewer than \code{bootstrap$returned}).}
#'     \item{\code{n_dropped}, \code{B_required},
#'       \code{bootstrap_status}}{\code{B} minus \code{B_used}, the fewest
#'       replications a request may ask for, and the verdict, in the names
#'       every bootstrap result of the package
#'       carries.}
#'     \item{\code{ci_method}}{The interval of each index, named
#'       \code{"factor:index"}: the method requested, or the percentile
#'       label with the reason when a BCa interval could not be formed.}
#'     \item{\code{ci_method_requested}, \code{conf_level},
#'       \code{improper_rule}, \code{denominator}}{The settings of the
#'       call.}
#'     \item{\code{metric}}{As for \code{bifactor_indices()}:
#'       \code{"observed"} for continuous items and
#'       \code{"latent_response"} for ordered items.}
#'     \item{\code{general}}{The name of the general factor.}
#'     \item{\code{improper}}{Whether the point estimate is improper
#'       (possible only with \code{improper = "keep"}).}
#'     \item{\code{N}, \code{N_incomplete}}{The rows resampled and how
#'       many of them have a missing value.}
#'   }
#'   \code{generics::tidy()} returns one row per index with \code{term}
#'   (\code{"factor:index"}), \code{estimate}, \code{se},
#'   \code{ci_lower}, \code{ci_upper}, and \code{conf_level};
#'   \code{generics::glance()} returns one row with the counts and the
#'   settings.
#'
#' @references
#' Aguirre-Urreta, M. I., Rönkkö, M., & McIntosh, C. N. (2019). A
#'   cautionary note on the finite sample behavior of maximal reliability.
#'   \emph{Psychological Methods, 24}(2), 236--252.
#'   \doi{10.1037/met0000176}
#'
#' DiCiccio, T. J., & Efron, B. (1996). Bootstrap confidence intervals.
#'   \emph{Statistical Science, 11}(3), 189--228.
#'
#' Efron, B., & Tibshirani, R. J. (1993). \emph{An introduction to the
#'   bootstrap}. New York, NY: Chapman & Hall/CRC.
#'
#' Hancock, G. R., & Mueller, R. O. (2001). Rethinking construct
#'   reliability within latent variable systems. In R. Cudeck, S. du Toit,
#'   & D. Sörbom (Eds.), \emph{Structural equation modeling: Present and
#'   future} (pp. 195--216). Scientific Software International.
#'
#' Kelley, K., & Pornprasertmanit, S. (2016). Confidence intervals for
#'   population reliability coefficients: Evaluation of methods,
#'   recommendations, and software for composite measures.
#'   \emph{Psychological Methods, 21}(1), 69--92. \doi{10.1037/a0040086}
#'
#' Raykov, T., & Penev, S. (2006). A direct method for obtaining
#'   approximate standard error and confidence interval of maximal
#'   reliability for composites with congeneric measures.
#'   \emph{Multivariate Behavioral Research, 41}(1), 15--28.
#'   \doi{10.1207/s15327906mbr4101_2}
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
#' Yuan, K.-H., & Hayashi, K. (2003). Bootstrap approach to inference
#'   and power analysis based on three test statistics for covariance
#'   structure models. \emph{British Journal of Mathematical and
#'   Statistical Psychology, 56}, 93--110.
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso \code{\link{bifactor_indices}} for the point estimates;
#'   \code{\link{reliability_omega}} for a single-factor omega and its
#'   intervals; \code{\link{reliability_H}} for coefficient H.
#'
#' @family multivariate and latent variable methods
#'
#' @keywords multivariate
#'
#' @examples
#' # The nine-item bifactor model of ?bifactor_indices: one general factor
#' # and three orthogonal group factors, simulated so that resamples give
#' # proper solutions and the example runs quickly. (The 15-test
#' # Holzinger and Swineford model of Details gives an improper solution
#' # in a fifth of its resamples, each of which took seconds to fit.)
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
#'
#' # Every replication refits the bifactor model. B = 10 keeps the
#' # example quick and is below the 400 a reported interval needs; the
#' # option lets it run and the warning says so. A reported interval
#' # deserves the default B = 1000. A resample whose solution is improper
#' # is set aside, and the warning and the printed lines count it.
#' old <- options(dmar.bootstrap_exploratory = TRUE)
#' ci <- ci_bifactor_indices(fit, B = 10, seed = 113)
#' options(old)
#' ci
#'
#' # Every replication's indices, and what happened to it.
#' head(attr(ci, "bootstrap")$t[, c("g:omega", "g:omega_H", "g:H",
#'                                  "g:maximal_reliability")])
#' table(attr(ci, "bootstrap")$status)
#'
#' # The broom verbs: one row per index, and the counts in one row.
#' generics::tidy(ci)
#' generics::glance(ci)
#'
#' @export
ci_bifactor_indices <- function(fit, general = NULL,
                                denominator = c("model_implied", "observed"),
                                ci_method = c("percentile", "bca"),
                                conf_level = 0.95, B = 1000L, seed = NULL,
                                ncpus = 1L, improper = c("drop", "keep"),
                                progress = NULL) {
  denominator <- match.arg(denominator)
  ci_method <- match.arg(ci_method)
  improper <- match.arg(improper)
  progress <- .check_progress(progress)
  if (!requireNamespace("lavaan", quietly = TRUE)) {
    .dmar_refuse(c(
      "Package 'lavaan' is required. Install it with ",
      "install.packages(\"lavaan\")."), "package_required", package = "lavaan")
  }
  if (!requireNamespace("boot", quietly = TRUE)) {
    .dmar_refuse(c(
      "Bootstrap intervals require the 'boot' package. Install it ",
      "with install.packages(\"boot\")."), "package_required", package = "boot")
  }
  if (!inherits(fit, "lavaan") || inherits(fit, "blavaan")) {
    .dmar_refuse_argument(c(
      "'fit' must be a lavaan fit object (a fit from lavaan::cfa() or ",
      "lavaan::sem(), not a Bayesian fit)."), arg = "fit")
  }
  if (!is.numeric(conf_level) || length(conf_level) != 1L ||
      !is.finite(conf_level) || conf_level <= 0 || conf_level >= 1) {
    .dmar_refuse_argument(
      "'conf_level' must be a single number between 0 and 1.",
      arg = "conf_level")
  }
  # The number of replications follows the rule of ?dmar_bootstrap: 10
  # beyond each limit for the percentile interval (400 at 95 percent), 50
  # beyond each nominal limit for BCa (2000).
  B <- .check_boot_B(B, ci_method, conf_level)
  ncpus <- .check_ncpus(ncpus)

  # What a case-resampling bootstrap of this fit can and cannot describe.
  .fiml_check_boot_aux(fit, "Bifactor bootstrap")
  # A fit that did not converge has no point estimates for an interval to
  # describe: the check bifactor_indices() makes, with its message
  # (.lavaan_check_converged(), R/no_maximum_internals.R).
  .lavaan_check_converged(fit, "ci_bifactor_indices()")
  if (lavaan::lavInspect(fit, "ngroups") > 1L) {
    .dmar_refuse_argument(c(
      "ci_bifactor_indices() needs a single-group fit: resampling the ",
      "pooled rows of a multiple-group fit would let the group sizes ",
      "vary from one replication to the next. Fit each group ",
      "separately."), arg = NULL)
  }
  cluster <- tryCatch(lavaan::lavInspect(fit, "cluster"),
                      error = function(e) character(0))
  if (length(cluster)) {
    .dmar_refuse_argument(c(
      "'fit' names a cluster variable (", paste(cluster, collapse = ", "),
      "), so its cases are not independent and cannot be resampled ",
      "one at a time; ci_bifactor_indices() resamples independent ",
      "cases only."), arg = "fit")
  }
  weights <- tryCatch(lavaan::lavInspect(fit, "sampling.weights"),
                      error = function(e) NULL)
  if (length(weights)) {
    .dmar_refuse_argument(c(
      "'fit' uses sampling weights, which a resample of its rows ",
      "does not carry; ci_bifactor_indices() bootstraps unweighted ",
      "fits only."), arg = "fit")
  }

  # The point estimates, with bifactor_indices()'s refusals; its warning
  # about an improper solution is replaced by the rule's own message.
  point <- withCallingHandlers(
    .bifactor_indices(fit, general, denominator),
    warning = function(w) {
      if (startsWith(conditionMessage(w), "Improper (Heywood) solution")) {
        invokeRestart("muffleWarning")
      }
    })
  # lavaan's post-estimation check, made on the model's own parameters
  # when the fit has auxiliary variables (.sem_admissible(),
  # R/sem_fit_internals.R), as for every replication: lavaan's own check
  # also tests the residual covariances of an auxiliary variable, which
  # are not positive definite once it is related to the general factor
  # strongly enough, however proper the bifactor solution.
  post_ok <- withCallingHandlers(
    .sem_admissible(fit, explain = FALSE),
    warning = function(w) invokeRestart("muffleWarning"))
  point_improper <- isTRUE(attr(point, "improper")) || !isTRUE(post_ok)
  cells <- .bifactor_cells(point)
  estimate <- .bifactor_cell_values(point, cells)
  # Every index but PUC is bootstrapped, and so is the maximal reliability
  # of each factor, which is not defined for an improper solution, as a
  # replication's is not: under improper = "keep" an improper point
  # estimate is reported, and its maximal reliability rows carry no
  # estimate and no interval.
  no_maximal <- .bifactor_maximal_cells(cells) &
    (point_improper | !is.finite(estimate))
  estimate[no_maximal] <- NA_real_
  boot_cells <- cells[cells$index != "PUC" & !no_maximal, , drop = FALSE]
  boot_est <- .bifactor_cell_values(point, boot_cells)
  if (point_improper) {
    if (improper == "drop") {
      .dmar_refuse_argument(c(
        "The solution of 'fit' is improper (a standardized loading ",
        "beyond one or a negative residual variance), and with ",
        "improper = \"drop\" the intervals are built from proper ",
        "replications only, which would not describe it. Respecify ",
        "the model, or call with improper = \"keep\" to bootstrap it ",
        "as it is."), arg = "fit")
    }
    if (!all(is.finite(boot_est)) || any(boot_est < 0 | boot_est > 1)) {
      .dmar_refuse_argument(c(
        "The solution of 'fit' is improper and gives an index that ",
        "is not defined or not a proportion (H is not defined when a ",
        "standardized loading reaches one, and an omega above one ",
        "belongs to a composite whose implied error variance is ",
        "negative), so no interval can describe it. Respecify the ",
        "model."), arg = "fit")
    }
    .dmar_warn(c(
      "Improper (Heywood) solution: a standardized loading exceeds ",
      "one or a residual variance is negative. The indices and ",
      "their intervals are not trustworthy until the model is ",
      "respecified."), "heywood_warning")
  }

  pool <- .bifactor_boot_pool(fit)
  general_name <- as.character(point$factor[1L])
  ordered <- as.character(lavaan::lavNames(fit, "ov.ord"))
  categories <- lapply(stats::setNames(ordered, ordered), function(v) {
    sort(unique(pool$data[[v]][!is.na(pool$data[[v]])]))
  })
  pt <- lavaan::parTable(fit)
  pt <- pt[setdiff(names(pt), c("est", "se", "start"))]
  point_fn <- .bifactor_point_fn(
    model = pt, options = .bifactor_refit_options(fit, denominator),
    ordered = ordered, categories = categories, general = general_name,
    denominator = denominator, cells = boot_cells)

  # The replications must be refits of the model behind the point
  # estimates: refit to the rows of the fit itself, the point function
  # must reproduce them.
  check <- tryCatch(point_fn(pool$data), error = function(e) e)
  # A fit whose likelihood has no maximum is refused as DMAR's own fits
  # are (.no_max_refuse()); 'fit' is the user's lavaan fit, so there is no
  # DMAR syntax to print.
  if (inherits(check, "dmar_no_maximum")) {
    .no_max_refuse("ci_bifactor_indices", .no_max_fields(check))
  }
  if (inherits(check, "error")) {
    .dmar_refuse(c(
      "The model could not be refit to the data stored with 'fit': ",
      conditionMessage(check)), "computation_failed")
  }
  gap <- suppressWarnings(max(abs(as.numeric(check) - boot_est)))
  if (!is.finite(gap) || gap > 1e-6) {
    .dmar_refuse_argument(c(
      "Refitting the model to the data stored with 'fit' did not ",
      "reproduce its bifactor indices",
      if (is.finite(gap)) sprintf(" (largest difference %.3g)", gap),
      ", so the replications would not be refits of this model."), arg = "fit")
  }

  stat_names <- paste(boot_cells$factor, boot_cells$index, sep = ":")
  boot <- .bootstrap_ci(fun = "ci_bifactor_indices",
                        pool$data, point_fn, B = B,
                        conf_level = conf_level, kind = ci_method,
                        seed = seed, ncpus = ncpus, improper = improper,
                        progress = progress, stat_names = stat_names,
                        improper_note = .bifactor_improper_note(improper))

  out <- data.frame(factor = cells$factor, index = cells$index,
                    estimate = estimate, se = NA_real_,
                    ci_lower = NA_real_, ci_upper = NA_real_,
                    stringsAsFactors = FALSE)
  at <- match(stat_names, paste(cells$factor, cells$index, sep = ":"))
  out$se[at] <- boot$se
  out$ci_lower[at] <- boot$lower
  out$ci_upper[at] <- boot$upper
  # The bootstrap estimate of each bootstrapped index's bias, the mean of
  # the replicates its interval uses minus the estimate. The sample H and
  # maximal reliability are biased upward (see Details), and the printed
  # table gives theirs.
  boot$bootstrap$bias <- stats::setNames(boot$replicate_mean - boot_est,
                                         stat_names)
  res <- .as_dmar_tbl(out, conf_level = conf_level,
                      subclass = "dmar_ci_bifactor")
  attr(res, "bootstrap") <- boot$bootstrap
  attr(res, "B") <- B
  attr(res, "B_used") <- boot$bootstrap$B_used
  attr(res, "n_dropped") <- as.integer(B) - boot$bootstrap$B_used
  attr(res, "B_required") <- boot$bootstrap$required
  attr(res, "bootstrap_status") <- boot$bootstrap$bootstrap_status
  attr(res, "ci_method") <- stats::setNames(boot$ci_method, stat_names)
  attr(res, "ci_method_requested") <- ci_method
  attr(res, "improper_rule") <- improper
  attr(res, "denominator") <- denominator
  attr(res, "metric") <- attr(point, "metric")
  attr(res, "general") <- general_name
  attr(res, "improper") <- point_improper
  attr(res, "N") <- nrow(pool$data)
  attr(res, "N_incomplete") <- pool$incomplete
  res
}


#' @rdname ci_bifactor_indices
#' @param x A \code{dmar_ci_bifactor} table returned by
#'   \code{ci_bifactor_indices()}.
#' @param ... Unused.
#' @importFrom generics tidy
#' @exportS3Method generics::tidy
tidy.dmar_ci_bifactor <- function(x, ...) {
  data.frame(term = paste(x$factor, x$index, sep = ":"),
             estimate = as.numeric(x$estimate),
             se = as.numeric(x$se),
             ci_lower = as.numeric(x$ci_lower),
             ci_upper = as.numeric(x$ci_upper),
             conf_level = attr(x, "conf_level") %||% NA_real_,
             stringsAsFactors = FALSE)
}

#' @rdname ci_bifactor_indices
#' @importFrom generics glance
#' @exportS3Method generics::glance
glance.dmar_ci_bifactor <- function(x, ...) {
  bs <- attr(x, "bootstrap")
  count <- function(name) {
    v <- if (is.list(bs)) bs[[name]] else NULL
    if (is.null(v)) NA_integer_ else as.integer(v)
  }
  data.frame(n_factors = length(unique(x$factor)),
             N = attr(x, "N") %||% NA_integer_,
             B = attr(x, "B") %||% NA_integer_,
             B_used = attr(x, "B_used") %||% NA_integer_,
             set_aside = count("failed"),
             improper_set_aside = count("improper"),
             improper_kept = count("improper_kept"),
             not_converged = count("not_converged"),
             caution = if (is.list(bs)) isTRUE(bs$caution) else NA,
             conf_level = attr(x, "conf_level") %||% NA_real_,
             ci_method = attr(x, "ci_method_requested") %||% NA_character_,
             improper_rule = attr(x, "improper_rule") %||% NA_character_,
             denominator = attr(x, "denominator") %||% NA_character_,
             stringsAsFactors = FALSE)
}
