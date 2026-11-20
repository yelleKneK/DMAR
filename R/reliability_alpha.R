#' Coefficient Alpha With a Confidence Interval
#'
#' @description
#' Estimates coefficient \eqn{\alpha} for a homogeneous composite score and
#' returns a confidence interval for the population coefficient using one
#' of several documented methods.
#'
#' @details
#' Coefficient \eqn{\alpha} was first derived by Guttman (1945) as his
#' \eqn{\lambda_3} and was subsequently popularized by Cronbach (1951),
#' under whose name the coefficient is often informally cited. The
#' attribution to Cronbach is historically incomplete; the modern
#' literature increasingly refers to the coefficient simply as
#' \dQuote{coefficient alpha} (see, e.g., Sijtsma, 2009; Revelle &
#' Zinbarg, 2009). DMAR follows that convention.
#'
#' For a \eqn{J}-item composite \eqn{Y = \sum_j X_j} the population
#' coefficient is
#' \deqn{\alpha = \frac{J}{J-1}\left(1 - \frac{\sum_{j} \sigma_{j}^{2}}{\sigma_{Y}^{2}}\right),}
#' where \eqn{\sigma_{j}^{2}} is the variance of item \emph{j} and
#' \eqn{\sigma_{Y}^{2}} is the variance of the composite. The sample
#' estimate substitutes sample variances. Under classical test theory,
#' \eqn{\alpha} equals the population reliability of the composite when
#' the items are essentially \eqn{\tau}-equivalent (i.e., equal factor
#' loadings); when loadings differ, \eqn{\alpha} is a lower bound on the
#' population reliability. For well-behaved homogeneous measurement
#' instruments coefficient \eqn{\alpha} and McDonald's (1999) coefficient
#' \eqn{\omega} (\code{\link{reliability_omega}}) typically yield very
#' similar values; \eqn{\omega} extends to congeneric items
#' (heterogeneous loadings) without the lower-bound caveat. For
#' ordered categorical items see \code{\link{reliability_omega_categorical}}.
#'
#' Coefficient \eqn{\alpha} is at most one but, as a function of the
#' covariance matrix, has no lower bound: it falls below zero when the
#' item covariances are negative on balance, which can arise when an
#' item is scored in the opposite direction to the others and needs
#' reverse scoring. A reliability coefficient cannot be
#' negative, so a negative estimate is returned as computed, with a
#' warning that says so, and each confidence limit is clamped into
#' [0, 1]; the interval around an estimate far enough below zero is then
#' [0, 0], which says that the data support no positive reliability for
#' the items as keyed.
#'
#' \strong{Two estimators of the same coefficient.} The argument
#' \code{estimator} selects how \eqn{\alpha} is estimated from the data.
#' Both target the same population quantity, and they agree in the
#' population whenever the \eqn{\tau}-equivalent model holds; they differ
#' in a finite sample because they take different routes to it.
#' \describe{
#'   \item{\code{"analytic"}}{The default. The closed-form equation above,
#'   applied to the observed covariance matrix. This is the classical
#'   coefficient, the number a hand calculation produces, and it makes no
#'   assumption beyond those of classical test theory.}
#'   \item{\code{"model_implied"}}{The reliability implied by the
#'   \eqn{\tau}-equivalent (equal loadings) single-factor model fit by
#'   maximum likelihood. With a shared loading \eqn{\lambda} and error
#'   variances \eqn{\psi_j^{2}} the model implied reliability of the
#'   \eqn{J}-item composite is
#'   \deqn{\alpha = \frac{(J \lambda)^{2}}{(J \lambda)^{2} + \sum_{j} \psi_{j}^{2}},}
#'   with maximum likelihood estimates substituted. Estimating through the
#'   model brings inference the formula cannot provide: a delta method
#'   standard error, a robust (Satorra-Bentler) variant under
#'   nonnormality, the profile likelihood interval, and a fit assessment
#'   of the \eqn{\tau}-equivalence claim itself. When the equal-loadings
#'   claim is doubtful, the congeneric \code{\link{reliability_omega}} is
#'   the appropriate coefficient rather than either \eqn{\alpha}.}
#' }
#'
#' Users of \pkg{MBESS} will recognize these as
#' \code{ci.reliability(type = "alpha")} and
#' \code{ci.reliability(type = "alpha-cfa")} respectively.
#'
#' \strong{Missing data and auxiliary variables.} Listwise deletion is
#' unbiased only when the data are missing completely at random and is
#' inefficient always; the \code{mlmr} vignette develops the argument at
#' length. By default (\code{missing = "fiml"}) incomplete data are
#' therefore analyzed by full information maximum likelihood (FIML),
#' which keeps every case with at least one observed item and is
#' consistent and efficient under the weaker missing at random (MAR)
#' assumption; \code{missing = "listwise"} analyzes the complete rows
#' alone. Data in which no value is missing, and a covariance matrix, are
#' analyzed as they stand: FIML would give the same coefficient, and the
#' interval is the complete data one, Bonett's by default for the analytic
#' estimator. A row with every item
#' missing does not count as incomplete data: without auxiliary
#' variables, listwise deletion and FIML both drop it, so data that are
#' complete but for such rows (a blank line at the end of a spreadsheet,
#' say) are analyzed as they stand too. The default analysis is
#' therefore FIML whenever an item value is missing and the complete data
#' analysis otherwise. \code{\link{fa_confirmatory}} differs: it fits
#' complete data by FIML as well. The complete data analysis is used here
#' because the closed-form intervals (\code{"bonett"}, \code{"feldt"},
#' \code{"fisher"}, \code{"hakstian_whalen"}), \code{"adf"},
#' \code{"adf_logistic"}, and \code{"likelihood"} need it, FIML giving
#' none of them. When the default meets incomplete items and one of those intervals, the
#' call stops rather than deleting cases unasked, and names the two
#' choices: \code{missing = "listwise"}, which computes the interval from
#' the complete rows alone, or an interval FIML gives for the estimator
#' in use (\code{"ml"}, \code{"ml_logistic"}, with the model implied
#' estimator \code{"mlr"} and \code{"mlr_logistic"}, or a bootstrap
#' method); when \code{ci_method} is not supplied, the analytic
#' estimator's interval under FIML is \code{"ml"}. The \code{aux}
#' argument names auxiliary variables:
#' columns of \code{data} that are not part of the composite but are
#' correlated with the items or with the reasons values are missing.
#' They are entered as \emph{saturated correlates} (Graham, 2003):
#' correlated freely with each other and with every item's residual,
#' never loading on the factor and never entering the composite, so the
#' measurement model is undisturbed while FIML uses their information.
#' Beyond recovering information, a good auxiliary makes the MAR
#' assumption itself more plausible, since missingness that depends on
#' the auxiliary becomes MAR once the auxiliary is conditioned on
#' (Collins, Schafer, & Kam, 2001). A case whose items are all missing
#' is kept when an auxiliary variable is observed on it: the saturated
#' model estimates the items and the auxiliaries jointly, so such a case
#' informs the item moments through its auxiliaries, and setting it
#' aside would make the analysis conditional on having an observed item,
#' which biases the estimates when the items go missing together at
#' random given the auxiliaries. Supplying \code{aux} without
#' \code{missing} implies \code{missing = "fiml"}, whether or not a
#' value is missing; combining it with an explicit \code{missing =
#' "listwise"} is an error. How each estimator uses FIML: the model
#' implied estimator simply fits its model with
#' \code{missing = "ml"}; the analytic estimator applies the classical
#' formula, unchanged, to the FIML estimate of the item covariance
#' matrix (the EM estimate under a saturated model over the items and any
#' auxiliaries, the matrix \code{\link{cov_fiml}} returns), so the
#' estimand stays the classical coefficient and only the covariance
#' matrix it is computed from improves. Under
#' \code{missing = "fiml"} the available intervals are \code{"ml"} and
#' \code{"ml_logistic"} (both estimators; the standard error comes from
#' the FIML information matrix), \code{"mlr"} and \code{"mlr_logistic"}
#' (model implied estimator; the Yuan-Bentler robust standard error),
#' and the bootstrap methods, which resample rows (including the
#' partially observed ones) and refit by FIML on each replication. The
#' complete-data closed forms (\code{"feldt"}, \code{"fisher"},
#' \code{"bonett"}, \code{"hakstian_whalen"}), \code{"adf"}, and
#' \code{"likelihood"} are errors with \code{missing = "fiml"} rather
#' than silently reverting to listwise deletion.
#'
#' \code{missing = "two_stage"} estimates in two stages (Yuan & Bentler,
#' 2000; Savalei & Bentler, 2009): the means and covariances of the items
#' are first estimated by maximum likelihood under the saturated normal
#' model (the EM estimates, with any \code{aux} variables entering this
#' stage only), and the model is then fit to them by the complete data
#' maximum likelihood discrepancy function. The analytic coefficient is a
#' function of the first stage alone, which is the saturated maximum
#' likelihood estimate FIML uses too, so its value and its intervals are
#' those of \code{missing = "fiml"} computed from the same cases. Those
#' are the cases FIML uses: every case with an observed item or, with
#' \code{aux}, an observed auxiliary variable, since the first stage
#' estimates the items and the auxiliaries jointly and a case observed on
#' an auxiliary alone informs the item moments through it. Setting such a
#' case aside would condition the analysis on an observed item and bias
#' the estimates when the items went missing together at random given the
#' auxiliaries.
#' \pkg{lavaan}'s two-stage estimator reads such a case in its first stage
#' but leaves it out of its sample size, which would make its standard
#' errors too large; the model implied fits go through
#' \code{\link{fa_confirmatory}}, which computes them again with every case the
#' first stage read, and a message counts such cases (when a fit cannot be
#' computed again, a warning says that the standard error is too large and
#' that a bootstrap interval counts every case). The model
#' implied coefficient comes
#' from the tau-equivalent model fit in the second stage, with the
#' normal theory two-stage standard error (Savalei & Bentler, 2009, Eq.
#' 6) for \code{"ml"} and the one robust to nonnormality (Yuan & Bentler,
#' 2000; Savalei & Falk, 2014) for \code{"mlr"}. The intervals FIML
#' refuses are refused here too, \code{"likelihood"} also because no
#' likelihood of the observed data is maximized. Multiple imputation is a
#' different feature with a different interface and is out of scope
#' here.
#'
#' Both methods need the covariance of every pair of items and auxiliary
#' variables, and the data carry no information about the covariance of
#' a pair that no case observes together: lavaan only warns, and its EM
#' returns an arbitrary value for it. Nor do they carry any about how two
#' variables vary together when the cases that observe both hold a
#' single value of one of them (or are one case): the estimate of their
#' covariance then rests on the form of the normal model rather than on
#' the data, or is not determined at all. Either pair is refused by name
#' before any fit, with a pointer to \code{\link{missing_coverage}}, and
#' a refusal of listwise deletion offers \code{missing = "fiml"} only when
#' the data hold no such pair. In a bootstrap, a resample that leaves
#' such a pair (or, for the BCa acceleration, a jackknife deletion that
#' does) is set aside with the same message and counted, under either
#' method.
#'
#' Available confidence interval methods (set via \code{ci_method}). Some
#' belong to one estimator only, because a closed-form interval for the
#' sample coefficient and a model-based interval are not interchangeable;
#' requesting a method the chosen estimator cannot supply is an error that
#' names the estimator to use instead:
#' \describe{
#'   \item{\code{"feldt"}}{\emph{Analytic estimator only.} The
#'   \emph{F}-distribution interval of Feldt (1965),
#'   exact under multivariate normality and parallel items.}
#'   \item{\code{"fisher"}}{\emph{Analytic estimator only.} Fisher's
#'   \eqn{z'} transformation (Fisher, 1950).
#'   Tends to overcover (Padilla, Divers, & Newton, 2012) and is generally
#'   not recommended.}
#'   \item{\code{"bonett"}}{\emph{Analytic estimator only.} Bonett's (2002)
#'   log transformation. The default for that estimator; well-behaved under
#'   normality.}
#'   \item{\code{"hakstian_whalen"}}{\emph{Analytic estimator only.}
#'   Cube-root transformation of Hakstian and Whalen (1976).}
#'   \item{\code{"mlr"}, \code{"mlr_logistic"}}{\emph{Model implied estimator
#'   only.} Wald interval using the robust (Satorra & Bentler, 1994) standard error
#'   from the fitted model. The default for that estimator; recommended
#'   among the closed forms when item distributions deviate from normality
#'   (Kelley & Pornprasertmanit, 2016). Requires raw data.}
#'   \item{\code{"likelihood"}}{\emph{Model implied estimator only.} Profile
#'   likelihood interval: the set of population values not rejected by the
#'   likelihood ratio test under the \eqn{\tau}-equivalent model, located by
#'   refitting under the nonlinear constraint that the model implied
#'   reliability equals each candidate value. Respects [0, 1], is not forced
#'   to be symmetric, and works from raw data or covariance input. Requires
#'   \pkg{lavaan}. It is unavailable with the analytic estimator because the
#'   interval and the point estimate would then refer to different
#'   quantities: the interval profiles the model implied coefficient while
#'   the estimate is the sample coefficient, so under a misspecified model
#'   the interval can exclude the estimate it accompanies.}
#'   \item{\code{"ml"}, \code{"ml_logistic"}}{Available to both estimators,
#'   by the route each affords. With \code{"analytic"} it is the closed-form
#'   ML standard error of van Zyl, Neudecker, and Nel (2000), computed
#'   directly from the covariance matrix; with \code{"model_implied"} it is
#'   the delta method standard error from the fitted model. The
#'   \code{_logistic} variant applies Browne's (1982) logit transformation.}
#'   \item{\code{"adf"}, \code{"adf_logistic"}}{Available to both
#'   estimators. With \code{"analytic"} it is the asymptotic
#'   distribution-free standard error of Maydeu-Olivares, Coffman, and
#'   Hartmann (2007); with \code{"model_implied"} the model is fit by
#'   weighted least squares (Browne, 1984) and the delta method applied.
#'   Requires raw data and a relatively large sample size.}
#'   \item{\code{"bootstrap_se"}, \code{"bootstrap_se_logistic"},
#'   \code{"percentile"}, \code{"bca"}}{Available to both estimators.
#'   Nonparametric bootstrap intervals (Efron & Tibshirani, 1993): the
#'   rows of \code{data} are resampled with replacement \code{B} times
#'   and the chosen estimator is recomputed on each replication.
#'   \code{"percentile"} takes the interval limits from the empirical
#'   quantiles of the bootstrap estimates (for a 95 percent interval,
#'   the 2.5th and 97.5th percentiles); it respects [0, 1] and is not
#'   forced to be symmetric about the estimate, but its coverage
#'   degrades when the estimator is biased or its variance changes
#'   with the parameter. \code{"bca"} (bias-corrected and accelerated)
#'   adjusts the two quantile positions for exactly those features,
#'   estimating the median bias from the bootstrap distribution and
#'   the acceleration from the jackknife, and is second-order accurate
#'   where the percentile interval is first-order accurate (DiCiccio &
#'   Efron, 1996). The acceleration is estimated as \pkg{boot}
#'   estimates it, by regression on the resampling frequencies when at
#'   least as many replications as cases returned a value and by the
#'   jackknife otherwise, or when that regression fails; a leave-one-out
#'   refit that fails is counted and set aside rather than stopping the
#'   call, and when the acceleration
#'   or the bias correction cannot be estimated at all the percentile
#'   interval is returned in its place, with a warning and the
#'   substitution recorded in the \code{ci_method} attribute; see
#'   \code{\link{reliability_omega_categorical}} for the full account.
#'   \code{"bootstrap_se"} uses the standard deviation
#'   of the bootstrap estimates as a standard error in an ordinary
#'   normal-theory interval; the \code{_logistic} variant builds that
#'   interval on the logit scale so the endpoints respect [0, 1], with
#'   the standard deviation carried to that scale by the delta method,
#'   divided by \eqn{\hat\alpha(1 - \hat\alpha)}, as in Equation 32 of
#'   Kelley and Pornprasertmanit (2016) and in
#'   \code{MBESS::ci.reliability()}. When every replication returns the
#'   same value, as when every item is an exact linear function of the
#'   same variable (an item entered twice beside other items does not do
#'   it), the bootstrap distribution is a single point and every method
#'   returns that point at both limits, with a warning; so does the
#'   \code{_logistic} variant at a value of one (identical items), where
#'   the logit is infinite.
#'   Replications on which the estimator cannot be computed (for
#'   example, a model refit that does not converge) are set aside, and
#'   so is a replication in which an item takes a single value (a
#'   resample that has left out every case on which a rare item differs
#'   from the rest, or the jackknife deletion of such a case), which is
#'   refused as the data would be, by name.
#'   A replication follows the rule on the number of cases as the data do
#'   (see \code{missing}): a resample repeats rows, and one with no more
#'   distinct cases than the rule asks for is set aside, so near the
#'   rule's bound many are (at 12 cases and eight variables, about three
#'   resamples in four). When more than half are, no interval is reported
#'   (the standard error and the limits are \code{NA}): the replications
#'   left hold more distinct cases than a resample typically does, so
#'   their spread understates the sampling variability of the
#'   coefficient, and at the bound itself each holds every distinct case
#'   and, without repeated rows, gives the estimate. Without repeated rows
#'   that happens at fewer cases than about 1.65 times the count the rule compares them with (fewer than 14
#'   for eight variables).
#'   The interval is computed from the replications that return a
#'   value, whose quantiles the percentile and BCa limits are; the
#'   counts, and why the rest were set aside, travel in the
#'   \code{bootstrap} attribute (see \emph{Value}), print under the
#'   table, and are stated in a warning whenever a replication was set
#'   aside. Requires raw data and the \pkg{boot} package. The default
#'   \code{B = 10000} is an accuracy choice: the BCa adjustment pushes
#'   the working quantiles farther into the tails of the bootstrap
#'   distribution than the percentile interval uses, and stabilizing
#'   them takes more replications than the customary 2000; reduce
#'   \code{B} for exploration, not for a reported analysis. Bootstrap
#'   results vary from run to run; supply \code{seed} for
#'   reproducibility, and \code{ncpus} to spread the replications over
#'   several processes, which leaves the interval for a given seed
#'   unchanged. Under \code{estimator = "model_implied"} each
#'   replication refits the tau-equivalent model by maximum likelihood,
#'   whose estimate robust maximum likelihood shares, and a replication
#'   whose solution is improper is kept or set aside according to \code{improper}, as
#'   \code{\link{reliability_omega}} describes; an improper point
#'   estimate whose error variances sum below zero, so that the model
#'   implies a negative error variance for the composite, stops with a
#'   message rather than returning a coefficient above one. (The
#'   likelihood of the tau-equivalent model has a maximum: one item's
#'   error variance cannot run off below zero while its loading grows,
#'   since the loadings are equal and every covariance would grow with
#'   them, so the stop of class \code{dmar_no_maximum} that
#'   \code{\link{reliability_omega}} describes does not arise here.)}
#'   \item{\code{"none"}}{Return only the point estimate.}
#' }
#'
#' \strong{Two complete rows.} Coefficient \eqn{\alpha} needs only the
#' trace and the sum of the covariance matrix, so it is computed from any
#' two rows, but the covariance matrix of two rows has rank one: every
#' pair of items correlates +1 or -1, whatever the items, and the
#' estimate is a function of the difference between the two rows alone.
#' Feldt's interval is exact for any \emph{N} of at least 2 under
#' multivariate normality and parallel items (for 0/1 items, as in
#' \code{\link{reliability_kr20}}, it is an approximation; Feldt, 1965),
#' and it is the interval the refusals below name. The Hakstian-Whalen
#' interval is defined there too, but at \emph{N} = 2 its upper limit lies
#' above one for every data set whose alpha is below one at any confidence
#' level above about .90, and is clamped to one (at an alpha of one, which
#' two rows can give, the interval is \eqn{[1, 1]}). The rest are refused at \emph{N} = 2, each for
#' its own reason, as \code{"fisher"} is below \emph{N} = 4: Bonett's
#' variance, \eqn{2J/((J-1)(N-2))}, is infinite; the \code{"ml"} and
#' \code{"adf"} standard errors (and their \code{_logistic} variants),
#' which treat that matrix as the population's, are zero for any data
#' (the asymptotic variance of van Zyl, Neudecker, and Nel, 2000,
#' vanishes at a matrix of rank one, and
#' each case's centered cross products are then proportional to the
#' covariance matrix, which alpha does not change with); and every
#' bootstrap resample is one row twice or both rows, so the bootstrap
#' distribution has no spread. The model implied estimator needs
#' more complete rows than items (see \code{missing}).
#'
#' \strong{A \emph{t} reference distribution for the Wald intervals.}
#' The intervals \code{"ml"}, \code{"mlr"}, and \code{"adf"}, with their
#' \code{_logistic} variants, are Wald intervals: the estimate plus or
#' minus a normal quantile times its standard error. The argument
#' \code{wald_df} replaces the normal quantile with the \emph{t} quantile
#' on \code{wald_df} degrees of freedom and leaves the estimate and the
#' standard error as they are (von Hippel, 2016, Eq. 29); for the
#' \code{_logistic} variants the \emph{t} quantile is applied on the
#' logit scale. The default, \code{Inf}, is the normal reference, and a
#' finite value with any other interval is an error. A number is the analyst's
#' own choice of degrees of freedom, a sensitivity check;
#' \code{\link{reliability_omega}} discusses it further.
#'
#' The keyword \code{wald_df = "missing_information"}, which
#' \code{\link{reliability_omega}} and \code{\link{fa_confirmatory}} accept, is an
#' error here. The calibration study behind it (in the package's source
#' repository, \code{methodology/missing_data}; \code{\link{fa_confirmatory}}
#' describes its design) evaluated alpha's FIML interval (\code{missing =
#' "fiml"}, \code{ci_method = "ml"}) with normal data, one factor
#' models of six items (\emph{N} from 30 to 200) or 12 (\emph{N} = 50
#' and 100), and 15 or 30 percent of the values of half the items
#' missing. On the normal reference that interval was already close to
#' .95: over the cells with \eqn{N \le 100} and missing values its
#' coverage was .0043 from .95 on average (between .940 and .958 in
#' every cell). The \emph{t} reference on von Hippel's
#' \eqn{\max\{3, N(1 - \gamma) - 1\}}{max(3, N(1 - gamma) - 1)}
#' degrees of freedom came barely closer (.0042), and in one cell its
#' coverage reached .966, above the .96 that the study's decision rule,
#' fixed before the run, allowed, so no rule passed. That cell, with 12
#' items, \emph{N} = 50, and 30 percent missing completely at random, is
#' the one where the study's fit of the saturated model converged in
#' only 53 percent of the replications (in 100 replications of that
#' design this function's fit, which starts at the EM estimates,
#' converged in 91 and the study's in 49), so the coverage there
#' conditions on the samples the study's fit could handle; without it the
#' \emph{t} reference would have passed (.0035 against .0042, at most
#' .9595). The verdict is the rule's as fixed in advance, and alpha
#' keeps the normal reference unless a number is given; the \emph{t}
#' reference would change its interval little in any case. On the logit
#' scale none of von Hippel's four rules passed either. The printed table
#' names the reference distribution in a line under it, and the
#' \code{wald_df} attribute records it.
#'
#' \strong{Comparison with other packages.} The \pkg{psych} package
#' provides \code{\link[psych]{alpha}}, which reports coefficient
#' \eqn{\alpha} along with item-level diagnostics (item-total
#' correlations, alpha-if-item-deleted, and several alternative
#' coefficients). The emphasis in \pkg{psych} is broad exploratory
#' psychometric reporting. \code{reliability_alpha} in \pkg{DMAR}
#' differs in emphasis: it returns a single point estimate alongside a
#' confidence interval drawn from the methods compared in Kelley and
#' Pornprasertmanit (2016).
#'
#' @param data A numeric matrix or data frame of item scores (rows are
#'   respondents, columns are items). Either \code{data} or \code{S} must
#'   be supplied. How incomplete rows are handled is governed by
#'   \code{missing}: full information maximum likelihood by default,
#'   listwise deletion with \code{missing = "listwise"}, or two-stage
#'   maximum likelihood with \code{missing = "two_stage"}. When
#'   \code{aux} is supplied, \code{data} contains both the item columns
#'   and the auxiliary columns, and the columns not named in \code{aux}
#'   are the items. An item column that is not numeric (or logical), and
#'   an item that takes a single value in the rows analyzed (the
#'   complete rows under listwise deletion, the rows that observe it
#'   otherwise), are refused by name before anything is computed, as is
#'   a zero variance on the diagonal of \code{S}: such an item adds a
#'   constant to the composite, whose reliability it leaves as it is,
#'   while the coefficient would count it among the \eqn{J} items (alpha
#'   with it is \eqn{J(J-2)/(J-1)^2} times alpha without it).
#'   A composite that does not vary is refused too, as when the items are
#'   ipsative (their scores sum to the same total in every row): its
#'   reliability is not defined. Under \code{missing = "fiml"} and
#'   \code{"two_stage"} it is checked in the rows that observe every item
#'   when they outnumber the items, where a composite that does not vary
#'   leaves the likelihood of the saturated model without a maximum, and
#'   otherwise on the EM estimate of the covariance matrix.
#'   SPSS user-defined missing values, which
#'   \code{haven::read_sav(user_na = TRUE)} keeps as values of a
#'   \code{haven_labelled_spss} column, are missing, as SPSS treats
#'   them.
#' @param S A symmetric covariance matrix among the items. If \code{S} is
#'   supplied, \code{N} must also be supplied; raw-data-only confidence
#'   interval methods (\code{"adf"}, \code{"bootstrap_*"}, \code{"percentile"},
#'   \code{"bca"}) are then unavailable, as are \code{missing = "fiml"}
#'   and \code{missing = "two_stage"} named explicitly, and \code{aux},
#'   since a covariance matrix has no incomplete cases (the default of
#'   \code{missing} analyzes it as it stands).
#' @param N Total sample size; required when \code{S} (rather than
#'   \code{data}) is supplied.
#' @param estimator How \eqn{\alpha} is estimated: \code{"analytic"}
#'   (the default) applies the closed-form equation to the observed
#'   covariance matrix, and \code{"model_implied"} takes the reliability
#'   implied by the \eqn{\tau}-equivalent single-factor model fit by
#'   maximum likelihood. See \emph{Details}.
#' @param missing How incomplete rows of \code{data} are handled:
#'   \code{"fiml"} (the default; full information maximum likelihood,
#'   using every case with at least one observed item or observed
#'   auxiliary variable), \code{"listwise"} (complete-case analysis), or
#'   \code{"two_stage"} (two-stage maximum likelihood, using the same cases
#'   as FIML). The default analysis is FIML whenever an item value
#'   is missing and the complete data analysis otherwise, unlike
#'   \code{\link{fa_confirmatory}}, which fits complete data by FIML too:
#'   data in which no value is missing, data complete but for rows with
#'   every item missing (which both treatments drop), and \code{S} are
#'   analyzed as they stand, with the complete data default interval,
#'   since the complete-data intervals need that analysis; incomplete
#'   items are analyzed by FIML, or, when \code{ci_method} names an
#'   interval FIML cannot give, the call stops and names the choices. See
#'   \emph{Details}. With \code{missing = "listwise"}, when fewer than
#'   two rows are complete, listwise deletion leaves nothing to analyze,
#'   and the function stops before any fit, says so, and names the
#'   treatment that uses the incomplete rows. The analytic estimator needs
#'   only the trace and the sum of the covariance matrix, which any two rows give, so it is computed from
#'   two or more complete rows. The model implied estimator fits the
#'   single-factor model by maximum likelihood, which needs the
#'   covariance matrix positive definite, and \eqn{N} complete rows give
#'   a matrix of rank at most \eqn{N - 1}; with no more complete rows
#'   than items it stops before the fit, says why, and names the
#'   treatment that uses the incomplete rows. Ordinary FIML is evaluated
#'   by the observed values and missingness patterns, not a blanket count
#'   of complete cases. Analytic alpha and two-stage estimation need
#'   unrestricted moments. The likelihood diagnostic refuses a certified
#'   affine-collapse direction or a numerical affine instability, as
#'   described in \code{\link{cov_fiml}}; EM convergence alone does not
#'   establish existence or uniqueness. Model-implied alpha can retain a
#'   proper identified constrained FIML fit with H0-based errors when H1
#'   diagnostics are unavailable; the returned result states this limit.
#'   The complete rows behind the model-implied listwise estimator are
#'   counted as distinct rows. Bootstrap resamples retain a separate
#'   distinct-case estimability screen, which is a policy, not a proof
#'   that every smaller constrained FIML likelihood is unbounded. The
#'   refusal message for the
#'   analytic estimator names listwise deletion when two or more rows are
#'   complete and listwise deletion would accept them, and otherwise says
#'   why it cannot help either. A refusal of listwise
#'   deletion names \code{missing = "fiml"} only when that method could
#'   use the data, so the two refusals never point at each other, and
#'   when two columns share a name (which listwise deletion without
#'   \code{aux} reads by position, and that method refuses) it says to
#'   name them apart first.
#' @param aux Optional character vector naming auxiliary variable
#'   columns of \code{data}: entered as saturated correlates under full
#'   information maximum likelihood, and into the first stage under
#'   \code{missing = "two_stage"}. Supplying \code{aux} without
#'   \code{missing} implies \code{missing = "fiml"}. See \emph{Details}.
#' @param ci_method Method for constructing the confidence interval. See
#'   \emph{Details} for which methods each estimator supports. The default
#'   \code{NULL} resolves to \code{"bonett"} for the analytic estimator and
#'   to \code{"mlr"} for the model implied estimator (or \code{"ml"} when
#'   only a covariance matrix is supplied, since \code{"mlr"} needs raw
#'   data). Under FIML (the default of \code{missing} with incomplete
#'   items) and two-stage estimation the analytic default is
#'   \code{"ml"}, whose standard error comes from the information matrix
#'   of the saturated model; with complete data the default of
#'   \code{missing} keeps \code{"bonett"}.
#' @param conf_level Confidence level for the interval (1 - Type I error
#'   rate). Defaults to \code{0.95}.
#' @param B Number of bootstrap replications when a bootstrap method is
#'   selected, a single whole number. Defaults to \code{10000}. The fewest
#'   accepted depends on the interval method and the level: at the default
#'   \code{conf_level}, 400 for \code{"percentile"} (10 replications beyond
#'   each limit), 2000 for \code{"bca"} (50 beyond each nominal limit,
#'   because the adjusted limits sit farther into the tails), and 100 for
#'   \code{"bootstrap_se"} and \code{"bootstrap_se_logistic"} (a standard
#'   deviation of the replications); at 0.99 the first two are 2000 and
#'   10000. A smaller \code{B} is refused before anything is computed
#'   (class \code{dmar_count_invalid}). The count is held to the number
#'   requested.
#'   Replications set aside afterward are counted and reported in a warning
#'   of class \code{dmar_bootstrap_replications_dropped}, which names the
#'   function and says to request more when those that return a value are
#'   fewer than the number required; when fewer than half of the \code{B}
#'   replications return a value, the interval is withheld (limits and
#'   bootstrap standard error \code{NA}) under a warning of class
#'   \code{dmar_bootstrap_insufficient}. See \code{\link{dmar_bootstrap}},
#'   which also describes \code{options(dmar.bootstrap_exploratory = TRUE)}
#'   for a first look with a smaller \code{B}.
#' @param seed Random number seed used for bootstrap reproducibility.
#'   Defaults to \code{NULL}, under which the resamples are drawn from,
#'   and advance, the current random number stream; supply an integer
#'   for a reproducible interval. When set, the seed applies inside the
#'   function only, and the user's generator state is restored when the
#'   function returns.
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
#'   before \pkg{processx} is loaded (see \code{\link{sem_exact_fit_bootstrap}}). Ignored unless a bootstrap
#'   method is selected.
#' @param improper For \code{estimator = "model_implied"} with a
#'   bootstrap interval: how a replication whose tau-equivalent
#'   single-factor solution is improper (a negative error variance, a
#'   Heywood case) is treated when the solution still yields a
#'   coefficient. \code{"keep"} (the default) uses its estimate like any
#'   other, and \code{"drop"} sets it aside; either way such
#'   replications are counted, and one whose solution yields no
#'   coefficient is set aside under both. With \code{"drop"} an improper
#'   point estimate stops the call. See
#'   \code{\link{reliability_omega}} for the reasoning. The analytic
#'   estimator fits no model, so no replication of it can be improper,
#'   and the argument has no effect there beyond being recorded as the
#'   \code{improper_rule} of the \code{bootstrap} attribute.
#' @param wald_df The degrees of freedom of the \emph{t} distribution to
#'   which a Wald interval refers: \code{Inf} (the default) is the normal
#'   reference, and a single positive number replaces the normal quantile
#'   with the \emph{t} quantile on that many degrees of freedom. It
#'   applies to \code{"ml"}, \code{"mlr"}, \code{"adf"}, and their
#'   \code{_logistic} variants; a finite value with any other interval
#'   is an error. The keyword \code{"missing_information"} of
#'   \code{\link{reliability_omega}} is an error here: its calibration
#'   found that the \emph{t} reference gives alpha's FIML interval no
#'   gain and too much coverage. See \emph{Details}.
#'
#' @return A \code{data.frame} with columns \code{term} and \code{value}
#' and rows
#' \code{"estimate"} (sample coefficient \eqn{\alpha}),
#' \code{"se"} (standard error on the coefficient scale, \code{NA} for
#' methods that do not produce one; for the transformation-based
#' intervals \code{"fisher"}, \code{"bonett"}, and
#' \code{"hakstian_whalen"} it is the delta method back-transform of the
#' transformation-scale standard error, evaluated at the estimate),
#' \code{"se_transformed"} (only for those transformation-based
#' intervals: the standard error on the transformation scale, the
#' quantity the interval is built from, with the scale named in the
#' attribute \code{se_transform_scale}: \code{"fisher_z"},
#' \code{"log(1-alpha)"}, or \code{"cube_root"}),
#' \code{"ci_lower"} and \code{"ci_upper"} (each clamped into
#' [0, 1], so the limits stay in order),
#' \code{"conf_level"}, \code{"N"} (the cases the analysis used: the
#' complete cases under listwise deletion, and under \code{"fiml"} and
#' \code{"two_stage"} every case with at least one observed item or
#' observed auxiliary variable),
#' \code{"N_complete"}
#' (the complete cases, so the cost of listwise deletion is visible at
#' a glance; equal to \code{"N"} under listwise deletion and for
#' covariance input), and \code{"n_items"} (number of items). The selected
#' coefficient, CI method, and missing-data treatment travel as the
#' attributes \code{coefficient}, \code{ci_method}, \code{missing},
#' \code{wald_df} (the degrees of freedom of the Wald interval's \emph{t}
#' reference, \code{Inf} for the normal one and for intervals that are
#' not Wald intervals), and (when supplied) \code{aux}; bootstrap calls
#' also record \code{B}. A bootstrap call also carries the attribute
#' \code{bootstrap}, a list recording the replications the interval
#' conditions on: \code{attempted} (\code{B}), \code{returned} (the
#' replications whose value the interval uses; the standard error is
#' their standard deviation and the percentile and BCa limits are their
#' quantiles), \code{failed} (the difference), \code{reasons} (a named
#' count of why replications were set aside: \code{error} when the refit
#' raised an error, \code{non_finite} when it was set aside without one,
#' as a fit that did not converge is, and as an improper one is under
#' \code{improper = "drop"} even though its estimate is finite),
#' \code{error_messages} (the
#' distinct error messages with their counts, empty when none),
#' \code{undetermined} (the replications set aside because the
#' resample left a polychoric correlation undetermined, which only
#' \code{\link{reliability_omega_categorical}} counts; zero here),
#' \code{improper} and \code{improper_kept} (the replications whose
#' tau-equivalent single-factor solution, under
#' \code{estimator = "model_implied"}, was improper, with a negative
#' error variance: those set aside, which are among the
#' \code{non_finite} ones, and those whose estimate the interval uses,
#' which are among the \code{returned} ones), \code{no_solution} (zero
#' here: it counts the replications of three items with no congeneric
#' single-factor solution, which \code{\link{reliability_omega}} meets),
#' \code{improper_rule} (\code{"keep"} or \code{"drop"}; see the
#' argument \code{improper}), \code{not_converged} (zero here: a fit
#' that did not converge is counted among the \code{non_finite} ones),
#' \code{too_few_cases} (the replications set aside because the
#' resample held too few distinct cases for the rule on the number of
#' cases, which are among the \code{error} ones; when they number more
#' than half of \code{B}, no interval is formed),
#' \code{caution} (\code{TRUE} when the replications set aside and the
#' improper ones kept number at least
#' \eqn{(1 - C)/2} of \code{B}, with \eqn{C} the confidence level, so
#' that how they were treated can decide where a limit falls),
#' \code{require_half} (\code{TRUE}: an interval requires at least half of
#' all attempted replications to return defined estimates; this is an
#' estimability safeguard, not a coverage
#' guarantee), \code{required} (the fewest replications a request for an
#' interval of this method and level may ask for, a count held to \code{B}
#' and not to the replications that return a value: 400 for
#' \code{"percentile"}, 2000 for \code{"bca"}, 100 for
#' the two bootstrap standard error methods at the 95 percent level; see
#' \code{\link{dmar_bootstrap}}), \code{bootstrap_status}
#' (\code{"complete"}, \code{"dropped"}, \code{"exploratory"},
#' \code{"insufficient_count"}, or \code{"insufficient_share"}; under the
#' last two the limits are \code{NA}), \code{B_used} (the smallest count
#' any reported limit rests on), \code{target_returned} and \code{target_withheld} (counts
#' and withholding flags for each target in a vector bootstrap),
#' \code{acceleration} and \code{acceleration_method} (the BCa
#' acceleration used and whether it came from the regression or the
#' jackknife; \code{NA} unless a BCa interval was formed),
#' \code{jackknife} (the same \code{attempted}, \code{returned},
#' \code{failed}, \code{reasons}, \code{error_messages},
#' \code{undetermined}, \code{improper}, \code{improper_kept},
#' \code{no_solution}, \code{not_converged}, and
#' \code{too_few_cases} for the leave-one-out refits behind the
#' acceleration, all zero when no jackknife ran), \code{ncpus} (the
#' number of processes used), \code{t} (the estimate of every
#' replication, in order, \code{NA}
#' where there is none; an improper estimate set aside under
#' \code{improper = "drop"} is kept here, so the interval under the
#' other rule can be formed without refitting), and \code{status} (a
#' factor giving the fate of every replication, with levels
#' \code{"proper"}, \code{"improper"}, \code{"error"},
#' \code{"non_finite"}, \code{"undetermined"}, and
#' \code{"not_converged"}). The interval uses the replications whose
#' status is \code{"proper"} and, under
#' \code{improper = "keep"}, those whose status is \code{"improper"}
#' and whose \code{t} is finite. The attribute \code{B_used} is the
#' number of replications the interval uses (\code{returned}). The attributes \code{n_dropped}
#' (\code{B} minus \code{B_used}), \code{B_required}, and
#' \code{bootstrap_status} repeat the counts in the names every bootstrap
#' result of the package carries. The warning has class
#' \code{dmar_bootstrap_replications_dropped} when replications were set
#' aside and the interval was formed, and
#' \code{dmar_bootstrap_insufficient} when the interval was withheld.
#' A BCa interval is also withheld, with the second class, when fewer than
#' 10 replications lie beyond a limit where the bias correction and the
#' acceleration place it (see \code{\link{dmar_bootstrap}}).
#' The printed table shows these counts in one line, and a warning states
#' them whenever a replication was set aside or the caution applies.
#'
#' @references
#' Bonett, D. G. (2002). Sample size requirements for testing and estimating
#'   coefficient alpha. \emph{Journal of Educational and Behavioral
#'   Statistics, 27}(4), 335--340. \doi{10.3102/10769986027004335}
#'
#' Browne, M. W. (1982). Covariance structures. In D. M. Hawkins (Ed.),
#'   \emph{Topics in applied multivariate analysis} (pp. 72--141).
#'   Cambridge, UK: Cambridge University Press.
#'
#' Browne, M. W. (1984). Asymptotically distribution-free methods for the
#'   analysis of covariance structures. \emph{British Journal of
#'   Mathematical and Statistical Psychology, 37}, 62--83.
#'
#' Collins, L. M., Schafer, J. L., & Kam, C.-M. (2001). A comparison of
#'   inclusive and restrictive strategies in modern missing data
#'   procedures. \emph{Psychological Methods, 6}, 330--351.
#'   \doi{10.1037/1082-989X.6.4.330}
#'
#' Cronbach, L. J. (1951). Coefficient alpha and the internal structure of
#'   tests. \emph{Psychometrika, 16}(3), 297--334.
#'
#' DiCiccio, T. J., & Efron, B. (1996). Bootstrap confidence intervals.
#'   \emph{Statistical Science, 11}(3), 189--228.
#'
#' Efron, B., & Tibshirani, R. J. (1993). \emph{An introduction to the
#'   bootstrap}. New York, NY: Chapman & Hall/CRC.
#'
#' Feldt, L. S. (1965). The approximate sampling distribution of
#'   Kuder-Richardson reliability coefficient twenty.
#'   \emph{Psychometrika, 30}, 357--370.
#'
#' Fisher, R. A. (1950). \emph{Statistical methods for research workers}.
#'   Edinburgh, UK: Oliver & Boyd.
#'
#' Graham, J. W. (2003). Adding missing-data-relevant variables to
#'   FIML-based structural equation models. \emph{Structural Equation
#'   Modeling, 10}(1), 80--100. \doi{10.1207/S15328007SEM1001_4}
#'
#' Guttman, L. (1945). A basis for analyzing test-retest reliability.
#'   \emph{Psychometrika, 10}(4), 255--282.
#'
#' Hakstian, A. R., & Whalen, T. E. (1976). A \emph{k}-sample significance
#'   test for independent alpha coefficients. \emph{Psychometrika, 41},
#'   219--231.
#'
#' Kelley, K. (2007a). Confidence intervals for standardized effect
#'   sizes: Theory, application, and implementation. \emph{Journal of
#'   Statistical Software, 20}(8), 1--24. \doi{10.18637/jss.v020.i08}
#'
#' Kelley, K. (2007b). Methods for the behavioral, educational, and
#'   social sciences: An R package. \emph{Behavior Research Methods,
#'   39}(4), 979--984. \doi{10.3758/BF03192993}
#'
#' Kelley, K., & Cheng, Y. (2012). Estimation of and confidence interval
#'   formation for reliability coefficients of homogeneous measurement
#'   instruments. \emph{Methodology, 8}, 39--50.
#'   \doi{10.1027/1614-2241/a000036}
#'
#' Kelley, K., & Pornprasertmanit, S. (2016). Confidence intervals for
#'   population reliability coefficients: Evaluation of methods,
#'   recommendations, and software for composite measures.
#'   \emph{Psychological Methods, 21}, 69--92. \doi{10.1037/a0040086}
#'
#' Satorra, A., & Bentler, P. M. (1994). Corrections to test statistics
#'   and standard errors in covariance structure analysis. In A. von Eye &
#'   C. C. Clogg (Eds.), \emph{Latent variables analysis: Applications for
#'   developmental research} (pp. 399--419). Thousand Oaks, CA: Sage.
#'
#' Savalei, V., & Bentler, P. M. (2009). A two-stage approach to missing
#'   data: Theory and application to auxiliary variables. \emph{Structural
#'   Equation Modeling, 16}, 477--497. \doi{10.1080/10705510903008238}
#'
#' Savalei, V., & Falk, C. F. (2014). Robust two-stage approach outperforms
#'   robust full information maximum likelihood with incomplete nonnormal
#'   data. \emph{Structural Equation Modeling, 21}(2), 280--302.
#'   \doi{10.1080/10705511.2014.882692}
#'
#' Yuan, K.-H., & Bentler, P. M. (2000). Three likelihood-based methods
#'   for mean and covariance structure analysis with nonnormal missing
#'   data. \emph{Sociological Methodology, 30}, 165--200.
#'
#' Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). \emph{Designing
#'   experiments and analyzing data: A model comparison perspective}
#'   (4th ed.). Routledge.
#'
#' Maydeu-Olivares, A., Coffman, D. L., & Hartmann, W. M. (2007).
#'   Asymptotically distribution-free (ADF) interval estimation of
#'   coefficient alpha. \emph{Psychological Methods, 12}, 157--176.
#'   \doi{10.1037/1082-989X.12.2.157}
#'
#' McDonald, R. P. (1999). \emph{Test theory: A unified treatment}.
#'   Mahwah, NJ: Lawrence Erlbaum Associates.
#'
#' Padilla, M. A., Divers, J., & Newton, M. (2012). Coefficient alpha
#'   bootstrap confidence interval under nonnormality.
#'   \emph{Applied Psychological Measurement, 36}, 331--348.
#'   \doi{10.1177/0146621612445470}
#'
#' Revelle, W., & Zinbarg, R. E. (2009). Coefficients alpha, beta, omega,
#'   and the GLB: Comments on Sijtsma.
#'   \emph{Psychometrika, 74}, 145--154. \doi{10.1007/s11336-008-9102-z}
#'
#' Sijtsma, K. (2009). On the use, the misuse, and the very limited
#'   usefulness of Cronbach's alpha. \emph{Psychometrika, 74}, 107--120.
#'   \doi{10.1007/s11336-008-9101-0}
#'
#' Terry, L. J., & Kelley, K. (2012). Sample size planning for composite
#'   reliability coefficients: Accuracy in parameter estimation via narrow
#'   confidence intervals. \emph{British Journal of Mathematical and
#'   Statistical Psychology, 65}, 371--401.
#'   \doi{10.1111/j.2044-8317.2011.02030.x}
#'
#' van Zyl, J. M., Neudecker, H., & Nel, D. G. (2000). On the distribution
#'   of the maximum likelihood estimator of Cronbach's alpha.
#'   \emph{Psychometrika, 65}(3), 271--280. \doi{10.1007/BF02296146}
#'
#' von Hippel, P. T. (2016). New confidence intervals and bias comparisons
#'   show that maximum likelihood can beat multiple imputation in small
#'   samples. \emph{Structural Equation Modeling, 23}(3), 422--437.
#'   \doi{10.1080/10705511.2015.1047931}
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso
#' \code{\link{reliability}} (general wrapper that dispatches by
#' coefficient),
#' \code{\link{reliability_omega}},
#' \code{\link{reliability_omega_categorical}},
#' \code{\link{reliability_kr20}},
#' \code{\link{cfa_1}} (single-factor CFA used internally),
#' \code{\link{ss_aipe_reliability}},
#' \code{\link[psych]{alpha}}.
#'
#' @examples
#' set.seed(113)
#' # Simulate six tau-equivalent items with population reliability ~ .8.
#' J <- 6
#' loadings <- rep(0.6, J)
#' eta <- rnorm(200)
#' errors <- matrix(rnorm(200 * J, sd = sqrt(1 - 0.6^2)), 200, J)
#' items <- outer(eta, loadings) + errors
#' colnames(items) <- paste0("y", seq_len(J))
#'
#' # Default (Bonett's transformation) CI from raw data.
#' reliability_alpha(data = items)
#'
#' # Same point estimate from a covariance matrix; CI requires N.
#' S <- cov(items)
#' reliability_alpha(S = S, N = 200, ci_method = "feldt")
#'
#' # The bootstrap intervals resample the rows and recompute the
#' # coefficient on each replication. The percentile interval reads its
#' # limits off the empirical quantiles of the bootstrap estimates.
#' # B = 500 keeps the example quick (a 95 percent percentile interval
#' # is reported from 400 or more); a reported interval deserves the
#' # default B = 10000, and seed makes the interval reproducible.
#' reliability_alpha(data = items, ci_method = "percentile", B = 500,
#'                   seed = 113)
#'
#' # The bias-corrected and accelerated interval adjusts those two
#' # quantile positions for median bias and for acceleration, and is the
#' # better choice for a reported interval. The default B = 10000 is an
#' # accuracy choice, not a formality, since BCa works farther into the
#' # tails of the bootstrap distribution than the percentile interval
#' # does (see Details). For the same reason B = 2000 is the fewest
#' # replications a 95 percent BCa interval is run with.
#' reliability_alpha(data = items, ci_method = "bca", B = 2000, seed = 113)
#'
#' # An application that runs a long bootstrap (a Shiny app, for example)
#' # can follow it through a progress callback, which is called with the
#' # replications attempted so far and B, and can stop it by returning
#' # FALSE. The callback changes no result: this is the percentile
#' # interval above.
#' reliability_alpha(data = items, ci_method = "percentile", B = 500,
#'                   seed = 113, progress = function(done, total) {
#'                     if (done %% 250 == 0) message(done, " of ", total)
#'                   })
#'
#' # Full information maximum likelihood with an auxiliary variable.
#' # Missingness on y2 depends on an auxiliary z (missing at random given
#' # z), so listwise deletion is biased and FIML with z is not. Supplying
#' # aux implies missing = "fiml"; the N_complete row shows how many
#' # rows listwise deletion would have kept.
#' z <- eta + rnorm(200, sd = 0.5)
#' d <- data.frame(items, z = z)
#' d$y2[runif(200) < plogis(-1 + 1.5 * as.numeric(scale(z)))] <- NA
#' reliability_alpha(data = d, aux = "z")
#'
#' @concept Cronbach's alpha
#'
#' @keywords htest multivariate
#' @family reliability
#'
#' @export

reliability_alpha <- function(data = NULL, S = NULL, N = NULL,
                              estimator = c("analytic", "model_implied"),
                              missing = c("fiml", "listwise", "two_stage"),
                              aux = NULL,
                              ci_method = NULL,
                              conf_level = 0.95,
                              B = 10000,
                              seed = NULL,
                              ncpus = 1L,
                              improper = c("keep", "drop"),
                              wald_df = Inf,
                              progress = NULL) {
  .dmar_check_conf_level(conf_level)
  estimator <- match.arg(estimator)
  improper <- match.arg(improper)
  ncpus <- .check_ncpus(ncpus)
  progress <- .check_progress(progress)
  missing_supplied <- !missing(missing)
  missing <- match.arg(missing)
  if (!is.null(aux) && missing_supplied && missing == "listwise") {
    .dmar_refuse_argument(c(
      "'aux' cannot be combined with missing = \"listwise\": an ",
      "auxiliary variable carries information only through an ",
      "incomplete data method. Drop 'aux' or set missing = \"fiml\" ",
      "(or \"two_stage\")."), arg = "aux")
  }

  # The two estimators support different interval methods, because a
  # closed-form interval for the sample coefficient and an interval read
  # off a fitted model are not interchangeable. Naming the sets here keeps
  # the contract in one place and lets the error say which estimator to
  # use rather than only that the value was wrong. A supplied ci_method is
  # checked against the estimator before the missing data treatment is
  # resolved, so that the refusal of an interval FIML cannot give (below)
  # never names a remedy that meets this refusal next; an unsupplied one
  # follows the treatment, after it.
  shared  <- c("ml", "ml_logistic", "adf", "adf_logistic",
               "bootstrap_se", "bootstrap_se_logistic",
               "percentile", "bca", "none")
  ok_analytic <- c("bonett", "feldt", "fisher", "hakstian_whalen", shared)
  ok_model    <- c("mlr", "mlr_logistic", "likelihood", shared)
  ok <- if (estimator == "analytic") ok_analytic else ok_model

  if (!is.null(ci_method)) {
    ci_method <- ci_method[1L]
    if (!is.character(ci_method) ||
        !ci_method %in% c(ok_analytic, ok_model)) {
      .dmar_refuse_argument(c(
        "'ci_method' must be one of: ",
        paste(sQuote(unique(c(ok_analytic, ok_model))), collapse = ", "),
        "."), arg = "ci_method")
    }
    if (!ci_method %in% ok) {
      other <- if (estimator == "analytic") "model_implied" else "analytic"
      detail <- if (ci_method == "likelihood") {
        paste0("The profile likelihood interval is a property of the ",
               "tau-equivalent model: it profiles the model implied ",
               "coefficient, so pairing it with the analytic point estimate ",
               "would report an interval and an estimate for different ",
               "quantities, and under a misspecified model the interval can ",
               "exclude the estimate. ")
      } else if (ci_method %in% c("mlr", "mlr_logistic")) {
        paste0("The robust (Satorra-Bentler) standard error comes from a ",
               "fitted model, which the analytic estimator does not fit. ")
      } else {
        paste0("That interval is a closed form for the sample coefficient ",
               "computed from the covariance matrix, not a quantity the ",
               "fitted model supplies. ")
      }
      .dmar_refuse_argument(c(
        "ci_method = \"", ci_method, "\" is not available with ",
        "estimator = \"", estimator, "\". ", detail,
        "Use estimator = \"", other, "\", or choose one of: ",
        paste(sQuote(ok), collapse = ", "), "."), arg = NULL)
    }
  }

  # The default is FIML when a value is missing (supplying aux implies it
  # in any case, unless two-stage estimation was asked for, so the common
  # case needs one argument, not two); complete data and a covariance
  # matrix are analyzed as they stand, as DMAR 1.0.x analyzed them, so the
  # analytic estimator's default interval stays Bonett's there, and an
  # interval FIML cannot give stops the call rather than deleting cases
  # unasked (.relia_resolve_missing()). An unsupplied ci_method follows
  # the treatment, below.
  missing <- .relia_resolve_missing(
    missing, missing_supplied, data, aux,
    unavailable = if (!is.null(ci_method) &&
                      ci_method %in% .RELIA_COMPLETE_DATA_METHODS) {
      paste0("ci_method = \"", ci_method, "\"")
    },
    fiml_options = paste0(
      "keep FIML with ci_method = \"ml\", \"ml_logistic\", ",
      if (estimator == "model_implied") "\"mlr\", \"mlr_logistic\", ",
      "or a bootstrap method"))
  if (missing != "listwise" && is.null(data)) {
    .dmar_refuse_argument(c(
      "missing = \"", missing, "\"",
      if (missing == "fiml") " (and 'aux')", " requires raw 'data': a ",
      "covariance matrix has no incomplete cases for ",
      .relia_incomplete_label(missing), " to use."), arg = "data")
  }

  if (is.null(ci_method)) {
    # Each estimator keeps the default it had as a separate function; the
    # model implied default needs raw data, so fall back when given only S.
    # Under FIML (the default of missing with incomplete items, or by
    # request) the analytic default becomes "ml", whose
    # standard error comes from the FIML information matrix of the
    # saturated model; the complete-case closed forms (Bonett and its
    # siblings) are unavailable there.
    ci_method <- if (estimator == "analytic") {
      if (missing != "listwise") "ml" else "bonett"
    } else if (is.null(data)) {
      "ml"
    } else {
      "mlr"
    }
  }

  if (!is.numeric(conf_level) || length(conf_level) != 1L ||
      conf_level <= 0 || conf_level >= 1) {
    .dmar_refuse_argument(
      "'conf_level' must be a single number in (0, 1).", arg = "conf_level")
  }
  # A t reference distribution changes a Wald interval only.
  wald_df <- .wald_df_check(wald_df, "reliability_alpha")
  wald <- ci_method %in% .RELIA_WALD_METHODS
  .wald_df_require_wald(wald_df, wald, paste0(
    "and ci_method = \"", ci_method, "\" is not one; wald_df applies to ",
    "ci_method \"ml\", \"adf\", \"mlr\" (with estimator = ",
    "\"model_implied\"), and their _logistic variants."))

  bootstrap_methods <- c("bootstrap_se", "bootstrap_se_logistic",
                         "percentile", "bca")
  if (ci_method %in% bootstrap_methods) {
    .check_boot_B(B, ci_method, conf_level)
  }
  if (missing != "listwise") {
    # Only intervals whose sampling theory can respect the missingness
    # pattern are offered; nothing silently falls back to listwise
    # deletion. The delta method ("ml"/"mlr") standard errors come from
    # the FIML information matrix, or under two-stage estimation from the
    # two-stage sandwich, and the bootstrap resamples the partially
    # observed rows and refits by the same method. The remaining methods
    # are complete-data constructions: the Feldt, Fisher, Bonett, and
    # Hakstian-Whalen forms depend on a single complete-case N that
    # does not exist when rows are partially observed, ADF needs
    # complete raw-data fourth moments, and the profile likelihood
    # interval has not been extended to the FIML fit (and has no
    # likelihood of the observed data to profile under two-stage
    # estimation).
    if (ci_method %in% c("feldt", "fisher", "bonett", "hakstian_whalen")) {
      .dmar_refuse_argument(c(
        "ci_method = \"", ci_method, "\" is a complete-data closed ",
        "form: its standard error treats every case as fully ",
        "observed, which overstates the information in partially ",
        "observed rows. With missing = \"", missing, "\" use \"ml\" ",
        "(or \"ml_logistic\"), whose standard error comes from ",
        if (missing == "fiml") "the FIML information matrix"
        else "the two-stage estimates and so reflects the missingness",
        ", or a bootstrap method."), arg = NULL)
    }
    if (ci_method %in% c("adf", "adf_logistic")) {
      .dmar_refuse_argument(c(
        "ci_method = \"", ci_method, "\" is not available with ",
        "missing = \"", missing, "\": the asymptotic distribution-free ",
        "variance requires complete raw-data fourth moments, and DMAR ",
        "offers no ADF counterpart under ",
        .relia_incomplete_label(missing), ". Use \"ml\", \"mlr\", or ",
        "a bootstrap method."), arg = NULL)
    }
    if (ci_method == "likelihood") {
      .dmar_refuse_argument(c(
        if (missing == "fiml") {
     "ci_method = \"likelihood\" has not been implemented for missing = \"fiml\"."
   } else {
     paste0("ci_method = \"likelihood\" is not available with missing = ",
            "\"two_stage\": two-stage estimation maximizes no ",
            "likelihood of the observed data, so there is none to ",
            "profile.")
   }, " Use \"ml\", \"mlr\", or a bootstrap method."), arg = NULL)
    }
  }

  # The model implied alpha is fitted by the single-factor model, which
  # needs more complete rows than items under listwise deletion; the
  # analytic alpha needs only the sums of the covariance matrix.
  inp <- .relia_resolve_inputs(data = data, S = S, N = N, allow_S = TRUE,
                               missing = missing, aux = aux,
                               fun = "reliability_alpha",
                               needs_saturated = estimator == "analytic" || missing == "two_stage" ||
                                 identical(wald_df, "missing_information"),
                               single_factor = if (estimator ==
                                                   "model_implied") {
                                 "the model implied alpha"
                               })

  needs_raw <- ci_method %in% c("mlr", "mlr_logistic",
                                "adf", "adf_logistic", bootstrap_methods)
  if (needs_raw && is.null(inp$data)) {
    .dmar_refuse_argument(c(
      "ci_method = \"", ci_method, "\" requires raw 'data'; ",
      "a covariance matrix is not sufficient."), arg = "data")
  }

  if (estimator == "analytic") {
    if (missing != "listwise") {
      # The classical coefficient is still the classical coefficient,
      # computed by its own equation; only the covariance matrix it is
      # applied to changes, from the complete-case estimate to the FIML
      # estimate from the saturated model over the items and any
      # auxiliary variables. The "ml" interval reads the delta method
      # standard error of alpha off that same saturated fit, so it
      # reflects the actual missingness pattern. Under two-stage
      # estimation alpha is a function of the stage one moments alone,
      # which are that same saturated maximum likelihood estimate, and
      # their two-stage covariance is the inverse of the same observed
      # information, so the estimate and every interval are FIML's.
      fc <- .fiml_cov(inp$data, items = inp$items, aux = inp$aux,
                      se = if (ci_method %in% c("ml", "ml_logistic")) {
                        "standard"
                      } else {
                        "none"
                      }, missing = missing)
      if (!fc$converged) {
        .dmar_refuse(c(
          "The saturated FIML model for the item covariance matrix ",
          "did not converge. This happens more often with small ",
          "samples, many auxiliary variables, or sparse missingness ",
          "patterns; consider dropping auxiliary variables."), "not_converged")
      }
      estimate <- .alpha_from_S(fc$S, em = TRUE)
      df_wald <- if (wald) {
        .relia_wald_df(wald_df, "reliability_alpha", fit = fc, N = inp$N,
                       items = inp$items, aux = inp$aux,
                       J = inp$J)$df
      } else {
        Inf
      }
      crit <- .wald_crit(conf_level, df_wald)
      ci <- switch(
        ci_method,
        none = list(se = NA_real_, lower = NA_real_, upper = NA_real_),
        ml = list(se    = fc$alpha_se,
                  lower = estimate - crit * fc$alpha_se,
                  upper = estimate + crit * fc$alpha_se),
        ml_logistic = {
          lg <- .logistic_ci(estimate, fc$alpha_se, crit)
          list(se = fc$alpha_se, lower = lg$lower, upper = lg$upper)
        },
        bootstrap_se = ,
        bootstrap_se_logistic = ,
        percentile = ,
        bca = .bootstrap_ci(
          fun = "reliability_alpha",
          data = inp$data,
          # A replication in which an item or auxiliary variable takes a
          # single value is refused, as the data are
          # (.relia_checked_point()).
          point_fn = .relia_checked_point(function(d) {
            f <- .fiml_cov(d, items = inp$items, aux = inp$aux,
                           se = "none", missing = missing)
            if (!f$converged) NA_real_ else .alpha_from_S(f$S, em = TRUE)
          }, missing = missing, items = inp$items, aux = inp$aux,
          composite = TRUE, needs_saturated = TRUE),
          B = B, conf_level = conf_level, kind = ci_method, seed = seed,
          ncpus = ncpus, improper = improper, progress = progress
        )
      )
    } else {
      estimate <- .alpha_from_S(inp$S)
      # The closed forms are computed from the complete cases, which carry
      # no missing information.
      df_wald <- if (wald) {
        .relia_wald_df(wald_df, "reliability_alpha", fit = NULL, N = inp$N,
                       J = inp$J)$df
      } else {
        Inf
      }
      ci <- switch(
        ci_method,
        none = list(se = NA_real_, lower = NA_real_, upper = NA_real_),
        feldt = .ci_feldt(estimate, inp$N, inp$J, conf_level),
        fisher = .ci_fisher(estimate, inp$N, conf_level),
        bonett = .ci_bonett(estimate, inp$N, inp$J, conf_level),
        hakstian_whalen = .ci_hakstian_whalen(estimate, inp$N, inp$J,
                                              conf_level),
        ml = .ci_alpha_ml(estimate, inp$S, inp$N, inp$J,
                          conf_level, logistic = FALSE, wald_df = df_wald),
        ml_logistic = .ci_alpha_ml(estimate, inp$S, inp$N, inp$J,
                                   conf_level, logistic = TRUE,
                                   wald_df = df_wald),
        adf = .ci_alpha_adf(estimate, inp$data, conf_level,
                            logistic = FALSE, wald_df = df_wald),
        adf_logistic = .ci_alpha_adf(estimate, inp$data, conf_level,
                                     logistic = TRUE, wald_df = df_wald),
        bootstrap_se = ,
        bootstrap_se_logistic = ,
        percentile = ,
        bca = .bootstrap_ci(
          fun = "reliability_alpha",
          data = inp$data,
          # A replication in which an item takes a single value is
          # refused, as the data are (.relia_checked_point()): alpha with
          # it would count it among the J items.
          point_fn = .relia_checked_point(
            function(d) .alpha_from_S(stats::cov(d))),
          B = B, conf_level = conf_level, kind = ci_method, seed = seed,
          ncpus = ncpus, improper = improper, progress = progress,
          two_rows = .relia_two_rows_alternative()
        )
      )
    }
  } else {
    fit_spec <- switch(
      ci_method,
      none            = list(estimator = "ML",  se = "none"),
      likelihood      = list(estimator = "ML",  se = "none"),
      ml              = list(estimator = "ML",  se = "standard"),
      ml_logistic     = list(estimator = "ML",  se = "standard"),
      mlr             = list(estimator = "MLR", se = "robust.sem"),
      mlr_logistic    = list(estimator = "MLR", se = "robust.sem"),
      adf             = list(estimator = "WLS", se = "standard"),
      adf_logistic    = list(estimator = "WLS", se = "standard"),
      # The bootstrap families need only the estimate, which ML and MLR
      # share to the last bit; ML spares the robust test statistic MLR
      # computes and these never read (MLR through DMAR 1.0.1).
      list(estimator = "ML", se = "none")
    )
    if (missing == "fiml" && identical(fit_spec$se, "robust.sem")) {
      # With missing = "ml" the canonical robust standard error is the
      # Yuan-Bentler sandwich, lavaan's own default pairing for MLR
      # under FIML.
      fit_spec$se <- "robust.huber.white"
    }
    if (missing == "two_stage") {
      # Two-stage estimation fits by ML; the robust intervals take the
      # robust two-stage standard errors, chosen through `missing` below.
      fit_spec <- list(estimator = "ML",
                       se = if (fit_spec$se == "none") "none" else "standard")
    }

    fit <- .omega_fit_cfa(data = inp$data, S = inp$S, N = inp$N,
                          h1_diagnostic = inp$h1_diagnostic,
                          equal_loadings = TRUE,
                          estimator = fit_spec$estimator,
                          se = fit_spec$se,
                          missing = .relia_lavaan_missing(
                            missing,
                            robust = ci_method %in% c("mlr", "mlr_logistic")),
                          aux = inp$aux)
    # An improper solution that yields no coefficient is refused, one
    # that does is warned about (or refused under improper = "drop" with
    # a bootstrap), before the check for convergence
    # (.omega_point_gate()). The gate's first refusal, of a likelihood
    # with no maximum, does not arise for the tau-equivalent model, whose
    # equal loadings keep one item's error variance from running off
    # (every covariance would grow with its loading); `fn` names this
    # function for it all the same.
    .omega_point_gate(fit, improper = improper,
                      bootstrap = ci_method %in% bootstrap_methods,
                      model = "tau-equivalent single-factor",
                      fn = "reliability_alpha")
    if (!fit$converged) {
      .dmar_refuse(
        .relia_fit_failure(fit, model = "tau-equivalent single-factor",
                           estimator = fit_spec$estimator,
                           missing = missing), "computation_failed")
    }
    estimate <- fit$omega
    df_wald <- if (wald) {
      .relia_wald_df(wald_df, "reliability_alpha", fit = fit, N = inp$N,
                     items = inp$items, aux = inp$aux,
                     J = inp$J)$df
    } else {
      Inf
    }

    ci <- switch(
      ci_method,
      none = list(se = NA_real_, lower = NA_real_, upper = NA_real_),
      ml = .ci_omega_delta(fit, conf_level, logistic = FALSE,
                           wald_df = df_wald),
      ml_logistic = .ci_omega_delta(fit, conf_level, logistic = TRUE,
                                    wald_df = df_wald),
      mlr = .ci_omega_delta(fit, conf_level, logistic = FALSE,
                            wald_df = df_wald),
      mlr_logistic = .ci_omega_delta(fit, conf_level, logistic = TRUE,
                                     wald_df = df_wald),
      adf = .ci_omega_delta(fit, conf_level, logistic = FALSE,
                            wald_df = df_wald),
      adf_logistic = .ci_omega_delta(fit, conf_level, logistic = TRUE,
                                     wald_df = df_wald),
      likelihood = .ci_reliability_likelihood(
        data = inp$data, S = if (is.null(inp$data)) inp$S else NULL,
        N = inp$N, equal_loadings = TRUE, conf_level = conf_level
      ),
      bootstrap_se = ,
      bootstrap_se_logistic = ,
      percentile = ,
      bca = .bootstrap_ci(
        fun = "reliability_alpha",
        data = inp$data,
        # A replication in which an item (or auxiliary variable) takes a
        # single value is refused before the fit, as the data are
        # (.relia_checked_point()).
        point_fn = .relia_checked_point(function(d) {
          # Under two-stage estimation a replication fits the second stage
          # alone to the EM moments (.omega_fit_cfa_ts_estimates()), the
          # same estimate (to the optimizer's tolerance) without the
          # two-stage standard errors and test.
          .omega_point_value(if (missing == "two_stage") {
            .omega_fit_cfa_ts_estimates(d, items = inp$items,
                                        aux = inp$aux,
                                        equal_loadings = TRUE)
          } else {
            .omega_fit_cfa(
              data = d, equal_loadings = TRUE,
              estimator = "ML", se = "none",
              missing = .relia_lavaan_missing(missing),
              aux = inp$aux)
          })
        }, missing = missing,
        items = if (missing == "listwise") NULL else inp$items,
        aux = inp$aux, fitted = TRUE, composite = TRUE,
        needs_saturated = missing == "two_stage"),
        B = B, conf_level = conf_level, kind = ci_method, seed = seed,
        ncpus = ncpus, improper = improper, progress = progress
      )
    )
  }

  out <- .relia_result(
    estimate = estimate, se = ci$se,
    lower = ci$lower, upper = ci$upper,
    conf_level = conf_level, N = inp$N, J = inp$J,
    coefficient = "alpha",
    ci_method = ci$ci_method %||% ci_method,
    bootstrap = ci$bootstrap,
    B = if (ci_method %in% bootstrap_methods) B else NA_integer_,
    N_complete = inp$N_complete,
    missing = missing,
    aux = inp$aux,
    se_transformed = ci$se_transformed,
    se_transform_scale = ci$se_transform_scale
  )
  out <- .relia_h1_qualify(out, inp$h1_diagnostic)
  attr(out, "estimator") <- estimator
  attr(out, "wald_df") <- df_wald
  out
}
