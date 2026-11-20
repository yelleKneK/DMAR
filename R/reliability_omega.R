#' Coefficient Omega (McDonald) With a Confidence Interval
#'
#' @description
#' Estimates McDonald's (1999) coefficient \eqn{\omega} for a homogeneous
#' composite score from a single-factor confirmatory factor analysis model
#' and returns a confidence interval for the population coefficient. The
#' \code{denominator} argument selects whether the total variance in the
#' denominator of \eqn{\omega} is estimated directly from the data
#' (\code{"observed"}, the default: robust omega) or taken from the
#' fitted model (\code{"model_implied"}); see \emph{Details} for the
#' properties of each choice. A bootstrap confidence interval is never
#' run unless requested; see \code{ci_method}.
#'
#' @details
#' Coefficient \eqn{\omega} is a population coefficient of determination.
#' Write the composite \eqn{Y = \sum_j X_j} through its measurement
#' decomposition, \eqn{Y = E[Y \mid \eta] + e}, the regression of the
#' observed composite on the latent variable \eqn{\eta} it is intended to
#' measure. By the law of total variance, the proportion of composite
#' variance attributable to \eqn{\eta} is
#' \eqn{\mathrm{Var}(E[Y \mid \eta]) / \mathrm{Var}(Y)}
#' (McDonald, 1999, 2011). Under the
#' congeneric single-factor model, in which item \eqn{j} has factor
#' loading \eqn{\lambda_j} and error variance \eqn{\psi_j^{2}}, the
#' numerator equals \eqn{\left(\sum_j \lambda_j\right)^2} and the
#' population coefficient is
#' \deqn{\omega = \frac{\left(\sum_{j} \lambda_{j}\right)^{2}}{\sigma_{Y}^{2}},}
#' with \eqn{\sigma_{Y}^{2}} the variance of the composite. Coefficient
#' \eqn{\omega} relaxes the equal-loadings assumption underlying
#' coefficient \eqn{\alpha} (Guttman, 1945; Cronbach, 1951), so on
#' congeneric scales \eqn{\omega} estimates population reliability
#' whereas \eqn{\alpha} underestimates it.
#'
#' \strong{The denominator argument.} The sample estimate substitutes
#' \pkg{lavaan} estimates of the loadings into the numerator; the two
#' \code{denominator} settings differ in how \eqn{\sigma_{Y}^{2}} is
#' estimated.
#' \describe{
#'   \item{\code{"observed"} (default)}{The variance of the composite
#'   estimated directly from the data (the sum of all elements of the
#'   unrestricted covariance matrix, on the same maximum likelihood
#'   divisor as the fitted loadings). This is the coefficient that
#'   Kelley and Pornprasertmanit (2016) call hierarchical omega
#'   (\code{MBESS::ci.reliability(type = "hierarchical")}); the DMAR
#'   documentation refers to it as \emph{robust omega}, and it is the
#'   default because its denominator estimates the variance of the
#'   composite consistently whether or not the single-factor model is
#'   correctly specified.}
#'   \item{\code{"model_implied"}}{The total variance reproduced by the
#'   fitted single-factor model,
#'   \eqn{\left(\sum_j \hat\lambda_j\right)^2 + \sum_j \hat\psi_j^{2}}.
#'   This is the textbook form of \eqn{\hat\omega}, correct exactly when
#'   the one-factor model reproduces the composite variance.}
#' }
#' The choice matters only insofar as the single-factor model is
#' misspecified, and the properties of each setting can be stated
#' exactly.
#' \itemize{
#'   \item If the single-factor model is correctly specified, the two
#'   definitions coincide in the population, and both estimators
#'   converge to the same value. Nothing is given up, in that sense, by
#'   either choice.
#'   \item The observed total variance is a consistent estimator of
#'   \eqn{\mathrm{Var}(Y)} whether or not the single-factor model is
#'   correctly specified; the model implied total variance is consistent
#'   for \eqn{\mathrm{Var}(Y)} only when the model is correct.
#'   \item Under misspecification (for example, a minor unmodeled factor
#'   or correlated errors), only \code{"observed"} retains the
#'   coefficient of determination interpretation: the proportion of the
#'   variance of the composite users actually compute that is
#'   attributable to the fitted common factor. With
#'   \code{"model_implied"} the estimate becomes a ratio of two model
#'   derived quantities whose denominator is no longer the variance of
#'   any composite a user scores.
#'   \item Under a correctly specified model with normal items, the
#'   model implied denominator uses the model structure and can be a
#'   slightly more efficient estimator of \eqn{\sigma_{Y}^{2}} in finite
#'   samples. In the simulations of Kelley and Pornprasertmanit (2016)
#'   the practical differences between the two coefficients were
#'   negligible when the model held, while interval coverage under
#'   modest model error favored the observed denominator paired with a
#'   bootstrap interval.
#' }
#' Those simulation results are why Kelley and Pornprasertmanit (2016)
#' recommend the observed denominator with a bootstrap confidence
#' interval when unidimensionality is only approximate, which is common
#' with real items.
#'
#' Stated as one identity, the two denominators differ by
#' \eqn{\mathbf{1}'(\mathbf{S} - \hat{\Sigma})\mathbf{1}}{1'(S - Sigma-hat)1},
#' the sum of every residual covariance of the single-factor fit. The
#' model implied setting treats covariance the model omits as absent
#' from the composite; the observed setting counts it as variance of
#' the composite that the factor does not explain. The first describes
#' the composite as the fitted model describes it, the second the
#' composite users compute. Estimating the factor model from its own
#' implied moments and reporting the reliability of a scored composite
#' with the observed variance are therefore not in tension: the model is
#' needed for the numerator, which partitions the variance of the
#' composite into the part due to the factor and the rest, and is not
#' needed for the denominator, which is the variance of the composite
#' itself. The omega that \code{\link{cfa_1}} reports beside the
#' loadings, and that \code{\link{fa_confirmatory}} reports for continuous items,
#' is the model implied coefficient, a property of the fitted model.
#'
#' Two cautions frame the choice. First, no denominator repairs a
#' misspecified measurement model: the fitted loadings absorb part of
#' whatever structure the single-factor model omits, so the numerator is
#' affected under either setting. Assess the single-factor model (for
#' example with \code{\link{cfa_1}}) before interpreting any
#' \eqn{\omega} variant, and model real multidimensionality directly
#' rather than patching over it. A gap between the two denominators
#' shows that the model does not reproduce the composite variance, but
#' a small gap does not show that the loadings are free of such
#' structure: a doublet or a method factor shared by a few same-signed
#' items inflates their loadings while leaving the two denominators
#' close. The residual correlations of the single-factor fit are where
#' that structure usually shows, so the result carries them (attribute
#' \code{residual_cor}; see \emph{Value}) and the printed table names
#' the largest. One limit applies: when the shared structure sits on the
#' items that dominate the factor, the single factor absorbs it and the
#' misfit appears as residual correlations among the other items, so the
#' largest residual signals local dependence more reliably than it
#' locates it. Second, the naming history is worth
#' knowing. The observed denominator coefficient was introduced as
#' hierarchical omega (Kelley & Pornprasertmanit, 2016; \pkg{MBESS}
#' type \code{"hierarchical"}), a name motivated by a hierarchical
#' factor logic: model misfit is viewed as a set of minor common
#' factors (visible as residual correlations), the single factor is
#' retained as an approximation, and the coefficient isolates the
#' variance attributable to the general factor alone, expressed
#' relative to the observed variance of the unweighted composite. It is
#' not the bifactor coefficient \eqn{\omega_H} of Zinbarg, Revelle,
#' Yovel, and Li (2005), whose numerator comes from the general factor
#' loadings of an explicitly multidimensional model. Upon reflection,
#' the authors would have named the coefficient for its behavior
#' rather than for the hierarchical motivation, as observed omega or
#' robust omega; \pkg{DMAR} uses \emph{robust omega}, with the
#' qualifications that word requires. The robustness is to
#' misspecification of the \emph{total variance} only, since the
#' numerator remains model based under either denominator; it is not
#' the outlier robustness of Zhang and Yuan (2016), and it is separate
#' from the robust maximum likelihood standard errors available through
#' \code{ci_method}. Robust omega also shares a design principle with
#' categorical omega as Kelley and Pornprasertmanit (2016) computed it:
#' in both, the total variance in the denominator is not taken from the
#' fitted factor model. Green and Yang (2009) built the denominator of
#' their categorical coefficient from the model implied polychoric
#' correlations instead, and \code{\link{reliability_omega_categorical}}
#' offers both through its own \code{denominator} argument. The
#' reliability vignette develops this framing.
#'
#' Available confidence interval methods (set via \code{ci_method}):
#' \describe{
#'   \item{\code{"ml"}, \code{"ml_logistic"}}{Wald interval using the
#'   maximum likelihood standard error from \pkg{lavaan} (Raykov, 2002).
#'   The \code{_logistic} variant applies Browne's (1982) logit
#'   transformation (see below).}
#'   \item{\code{"mlr"}, \code{"mlr_logistic"}}{Wald interval using the
#'   robust (Satorra & Bentler, 1994) standard error, recommended among
#'   closed-form methods when item distributions deviate from normality
#'   (Kelley & Pornprasertmanit, 2016). \code{"mlr_logistic"} is the
#'   default for \code{denominator = "model_implied"}.}
#'   \item{\code{"likelihood"}}{Profile likelihood interval: the set of
#'   population values not rejected by the likelihood ratio test,
#'   located by refitting the model under the nonlinear constraint that
#'   the model implied reliability equals each candidate value. It
#'   respects the [0, 1] range and is not forced to be symmetric about
#'   the estimate. Maximum likelihood; available with raw data or
#'   covariance input, for \code{denominator = "model_implied"} only.}
#'   \item{\code{"adf"}, \code{"adf_logistic"}}{Wald interval using
#'   weighted least squares (\dQuote{ADF}) estimation (Browne, 1984).
#'   Requires raw data and a relatively large sample size.}
#'   \item{\code{"feldt"}, \code{"fisher"}, \code{"bonett"},
#'   \code{"hakstian_whalen"}}{Closed-form intervals derived for
#'   coefficient \eqn{\alpha}. They apply mechanically to \eqn{\omega} but
#'   their coverage performance for \eqn{\omega} is generally inferior to
#'   the maximum likelihood and bootstrap intervals (Kelley &
#'   Pornprasertmanit, 2016).}
#'   \item{\code{"bootstrap_se"}, \code{"bootstrap_se_logistic"},
#'   \code{"percentile"}, \code{"bca"}}{Nonparametric bootstrap intervals
#'   (Efron & Tibshirani, 1993): the rows of \code{data} are resampled
#'   with replacement \code{B} times and \eqn{\omega} is recomputed,
#'   with the factor model refit, on each replication.
#'   \code{"percentile"} takes the interval limits from the empirical
#'   quantiles of the bootstrap estimates; it respects [0, 1] and is
#'   not forced to be symmetric about the estimate. \code{"bca"}
#'   (bias-corrected and accelerated) adjusts the two quantile
#'   positions for median bias, estimated from the bootstrap
#'   distribution, and for the rate at which the estimator's variance
#'   changes with the parameter, the acceleration, estimated by the
#'   jackknife; those two adjustments make it second-order accurate
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
#'   of the bootstrap estimates as a standard error in a normal-theory
#'   interval, built on the logit scale for the \code{_logistic}
#'   variant so the endpoints respect [0, 1] (the standard deviation is
#'   carried to the logit scale by the delta method, divided by
#'   \eqn{\hat\omega(1 - \hat\omega)}, as in Equation 32 of Kelley and
#'   Pornprasertmanit, 2016, and in \code{MBESS::ci.reliability()}).
#'   When every replication returns the same value, the bootstrap
#'   distribution is a single point and every method returns that point
#'   at both limits, with a warning (items that are exact linear
#'   functions of one variable, which give alpha such replications, do
#'   not reach the bootstrap here: their covariance matrix is singular,
#'   and the single-factor fit behind the estimate stops on it).
#'   Replications whose model
#'   refit does not converge are set aside, and so is a replication in
#'   which an item takes a single value (a resample that has left out
#'   every case on which a rare item differs from the rest), which is
#'   refused before the refit, by name, as the data would be, and so is
#'   one with no more distinct cases
#'   than the rule on the number of cases for the data asks for (see
#'   \code{missing}; a resample repeats rows, so near the rule's bound
#'   many are, and when more than half are no interval is reported, since
#'   those left hold more distinct cases than a resample typically does
#'   and their spread understates the sampling variability, and at the
#'   bound each holds every distinct case; without repeated rows that
#'   happens at fewer cases than about 1.65 times the count the rule
#'   compares them with); the interval is computed
#'   from the replications that return a value, whose quantiles the
#'   percentile and BCa limits are; the counts, and why the rest were
#'   set aside, travel in the \code{bootstrap} attribute (see
#'   \emph{Value}), print under the table, and are stated in a warning
#'   whenever a replication was set aside. The default
#'   \code{B = 10000} is an accuracy choice: the BCa adjustment pushes
#'   the working quantiles farther into the tails of the bootstrap
#'   distribution than the percentile interval uses, and stabilizing
#'   them takes more replications than the customary 2000; reduce
#'   \code{B} for exploration, not for a reported analysis.
#'   Recommended when assumptions of parametric methods are
#'   questionable. Each replication refits the single-factor model by
#'   maximum likelihood: robust maximum likelihood has the same
#'   estimates and differs only in the standard errors and the test
#'   statistic, which a replication does not use, so the robust fit is
#'   not repeated \code{B} times. Requires
#'   raw data and the \pkg{boot} package; supply \code{seed} for
#'   run-to-run reproducibility, and \code{ncpus} to spread the
#'   replications over several processes, which leaves the interval for
#'   a given seed unchanged. How replications with an improper solution
#'   are treated is described next.}
#'   \item{\code{"none"}}{Return only the point estimate.}
#' }
#' \strong{Which intervals apply to which denominator.}
#' \describe{
#'   \item{\code{denominator = "model_implied"}}{Every method on the
#'   menu: the Wald intervals \code{"ml"}, \code{"mlr"}, and
#'   \code{"adf"} with their \code{_logistic} variants, the profile
#'   \code{"likelihood"} interval, the alpha-derived closed forms
#'   \code{"feldt"}, \code{"fisher"}, \code{"bonett"}, and
#'   \code{"hakstian_whalen"}, the four bootstrap methods, and
#'   \code{"none"}.}
#'   \item{\code{denominator = "observed"} (robust omega, the
#'   default)}{The four bootstrap methods (\code{"percentile"},
#'   \code{"bca"}, \code{"bootstrap_se"},
#'   \code{"bootstrap_se_logistic"}) and \code{"none"} only. The delta
#'   method standard errors, the profile likelihood, and the
#'   alpha-derived closed forms are derived under the model implied
#'   ratio and do not account for sampling of the observed denominator,
#'   so \code{"mlr_logistic"}, the default interval of the model implied
#'   coefficient, is not an interval for robust omega.}
#' }
#' A request for a combination that does not apply, such as
#' \code{ci_method = "mlr"} with the default denominator, stops with a
#' message naming the denominator that does; nothing falls back
#' silently. Covariance input (\code{S}) and the incomplete data
#' methods (\code{missing = "fiml"} and \code{"two_stage"}) narrow the
#' menu further, as the \code{S} argument and the paragraphs on missing
#' data below describe, and each stops the same way. The
#' pairing of robust omega with the bootstrap is not a limitation in
#' practice, since the bootstrap is the interval Kelley and
#' Pornprasertmanit (2016) recommend for the observed denominator
#' coefficient in any case. Because no analysis in \pkg{DMAR} runs a
#' bootstrap unless the user requests one, the default for robust omega
#' is the point estimate with no interval, accompanied by a message
#' naming the call that produces the recommended interval; request
#' \code{ci_method = "percentile"} or \code{"bca"} to obtain it.
#'
#' \strong{Limits on the scale of omega or on the logit scale.} The Wald
#' intervals come in two forms, both built from the same standard error.
#' \code{"ml"}, \code{"mlr"}, and \code{"adf"} take the estimate plus or
#' minus a normal quantile times the standard error, an interval
#' symmetric about the estimate whose upper limit can pass one (the
#' limits are then clamped into [0, 1]; see \emph{Value}). Their
#' \code{_logistic} variants form the limits on the logit scale,
#' \eqn{\log\{\hat\omega / (1 - \hat\omega)\}}{log(omega-hat / (1 -
#' omega-hat))}, whose standard error the delta method gives as
#' \eqn{\mathrm{SE} / \{\hat\omega (1 - \hat\omega)\}}{SE / (omega-hat (1
#' - omega-hat))}, and transform them back (Browne, 1982): the limits of a
#' monotone transformation of a parameter are the transformation of its
#' limits, the confidence interval transformation principle that Kelley
#' (2007a) sets out. The interval then lies within (0, 1) and is
#' asymmetric near one, as the sampling distribution of a high
#' reliability is; the \code{se} row is the standard error of
#' \eqn{\hat\omega}{omega-hat} itself under either form. The default is
#' \code{"mlr_logistic"} on the evidence of a Monte
#' Carlo study (in the package's source repository,
#' \code{methodology/reliability_limits}) of factors of three to six items
#' with \emph{N} from 50 to 500. The study formed both kinds of limits
#' from the same standard error, the normal theory maximum likelihood one
#' for normal data and the Huber-White one of MLR for nonnormal data, so
#' what it compared is the scale on which the limits are formed. With
#' normal data, under listwise deletion and FIML alike, the logit limits
#' covered from .939 to .956 in every cell (a mean absolute departure from
#' .95 of .003), where the limits on the scale of omega covered from .910
#' to .952 (.009 to .010), their misses falling mostly on one side, with
#' the population value below the lower limit. With nonnormal data fitted
#' by MLR both forms covered below .95, the logit form by less (a mean
#' absolute departure of .013 against .022). No cell used the robust.sem
#' standard error that \code{"mlr_logistic"} takes with complete data; it
#' estimates the same asymptotic variance as the Huber-White one, which
#' the default takes under FIML. \code{\link{fa_confirmatory}} forms the
#' limits of its omega row the same way, and \code{\link{reliability_H}}
#' those of \emph{H} by default, so the three give limits of one form.
#' With a covariance matrix, which the robust standard error cannot use,
#' the default is \code{"ml_logistic"}, as in \code{\link{reliability_H}}.
#' Under
#' \code{wald_df = "missing_information"} the default is
#' \code{"mlr"} (\code{"ml"} with a covariance matrix), since that rule is
#' calibrated for the intervals on the scale of omega only (see below).
#'
#' \strong{A \emph{t} reference distribution for the Wald intervals.}
#' The Wald intervals (\code{"ml"}, \code{"mlr"}, \code{"adf"}, and their
#' \code{_logistic} variants) are the estimate plus or minus a normal
#' quantile times its standard error. The normal quantile treats the
#' standard error as known, which it is only in the limit; in small
#' samples, and more so when values are missing, such intervals are too
#' short. Von Hippel (2016, pp. 431--433) showed this for maximum
#' likelihood estimates from incomplete data and proposed keeping the
#' estimate and its standard error while replacing the normal quantile
#' with a \emph{t} quantile (his Eq. 29). The argument \code{wald_df}
#' does that: the interval becomes \eqn{\hat{\omega} \pm
#' t_{1 - \alpha/2, \nu} \mathrm{SE}}{omega-hat +/- t[1 - alpha/2, nu] SE}
#' with \eqn{\nu} the value of \code{wald_df}, and for the
#' \code{_logistic} variants the \emph{t} quantile is applied on the logit
#' scale before the limits are transformed back. The default, \code{Inf},
#' is the normal reference. A number is the
#' analyst's own choice of degrees of
#' freedom, best used to see how much a plausible range of them would
#' lengthen the interval; with a few hundred cases the answer is very
#' little (on 200 degrees of freedom the 95\% quantile is 1.972 rather
#' than 1.960).
#'
#' The calibration study described next selected its simulated data sets
#' by an eligibility rule that keeps some draws the affine-collapse
#' diagnostic sets aside, so its results are conditional on that rule, and
#' coverage under the diagnostic's selection has not been established by
#' that study. A raw-data eligibility replay alone cannot validate the FMI
#' refits and information diagnostics. The keyword is an approximate,
#' opt-in rule, not a coverage guarantee.
#'
#' The keyword \code{wald_df = "missing_information"} takes the degrees
#' of freedom from omega's own fraction of missing information
#' \eqn{\gamma} by von Hippel's (2016, Eqs. 31 and 34) effective sample
#' size rule with \eqn{k = 1},
#' \deqn{\nu = \max\{3, N(1 - \gamma) - 1\},}{nu = max(3, N(1 - gamma) -
#' 1),}
#' where \eqn{N} is the number of cases the fit used and \eqn{\gamma} is
#' estimated from the fit as \code{\link{fraction_missing_information}}
#' estimates it; with complete data (listwise deletion, or a covariance
#' matrix) \eqn{\gamma = 0} and \eqn{\nu = N - 1}. The rule comes from a
#' calibration study of these intervals under full information maximum
#' likelihood (in the package's source repository,
#' \code{methodology/missing_data}; \code{\link{fa_confirmatory}} describes its
#' design), which simulated the model implied omega of one factor models
#' with six items (\emph{N} from 30 to 200) or 12 (\emph{N} = 50 and
#' 100), normal data fit by ML and nonnormal data by MLR. Over the cells
#' with \eqn{N \le 100} and missing values the coverage of the \emph{t}
#' interval was .0180 from .95 on average, against .0242 for the normal
#' one, and no higher than .9550
#' in any cell, with normal theory standard errors (\code{"ml"}) and
#' Huber-White ones (\code{"mlr"}) alike; with nonnormal data both
#' intervals still covered below .95, the \emph{t} one less so. On the
#' logit scale the rule did not pass: with normal data the
#' \code{_logistic} interval on the normal reference was already within
#' .0044 of .95 on average, and the \emph{t} reference raised its coverage
#' as high as .9645, above the .96 the study allowed. The keyword
#' therefore applies to \code{ci_method = "ml"} and \code{"mlr"} only
#' (with \code{ci_method} not supplied, the interval is \code{"mlr"}, or
#' \code{"ml"} with a covariance matrix), under listwise deletion or
#' \code{missing = "fiml"}, and is an error with the \code{_logistic}
#' variants, with \code{"adf"} and
#' \code{"adf_logistic"} (the study evaluated maximum likelihood standard
#' errors only), with two-stage estimation (for which no fraction of
#' missing information is computed), and with every interval that is not
#' a Wald interval. The study fit no auxiliary variables: with \code{aux}
#' the rule takes omega's fraction of missing information from the fit
#' that carries them, which is smaller when they predict the missing
#' values, and the study supports that case by analogy rather than
#' directly. Under listwise deletion the standard errors here come from
#' the expected information (\code{"ml"}, from a covariance matrix too)
#' and the robust.sem sandwich (\code{"mlr"}, which needs raw data), where the
#' study's cells without missing values used the
#' observed information and the Huber-White sandwich; for omega at
#' \emph{N} = 30 and 50 they are shorter by less than 1 and by 1 to 2
#' percent, so the complete data intervals are, if anything, slightly
#' shorter than those the study evaluated. The printed table names the
#' reference distribution in a line under it, under the keyword with
#' \emph{N} and the FMI; the \code{wald_df} attribute records the degrees
#' of freedom, and under the keyword the \code{wald_df_table} attribute
#' holds the \code{\link{df_missing}} table they came from.
#'
#' \strong{Missing data and auxiliary variables.} By default
#' (\code{missing = "fiml"}) incomplete data are analyzed by full
#' information maximum likelihood (FIML), which keeps every case with at
#' least one observed item or observed auxiliary variable and fits the
#' single-factor model to all of them, consistent and efficient under the
#' missing at random (MAR) assumption; \code{missing = "listwise"} fits
#' the complete rows alone. Under the default, the model implied
#' omega here and the omega \code{\link{fa_confirmatory}} reports for the
#' same items agree on incomplete data, the estimate and, with
#' \code{ci_method = "ml_logistic"}, the interval; the default interval,
#' \code{"mlr_logistic"}, forms its limits the same way from the robust
#' standard error. Data in which no value
#' is missing, and a covariance matrix, are analyzed as they stand: FIML
#' would give the same estimate, and the standard errors and intervals
#' come from the expected information (and for
#' \code{"mlr"} the robust.sem sandwich) where FIML takes the observed
#' information (and the Huber-White sandwich). A row with every item
#' missing does not count as incomplete data: without auxiliary
#' variables, listwise deletion and FIML both drop it, so data that are
#' complete but for such rows (a blank line at the end of a spreadsheet,
#' say) are analyzed as they stand too. The default analysis is
#' therefore FIML whenever an item value is missing and the complete data
#' analysis otherwise. \code{\link{fa_confirmatory}} differs: it fits
#' complete data by FIML as well. The complete data analysis is used here
#' because the closed-form intervals (\code{"feldt"}, \code{"fisher"},
#' \code{"bonett"}, \code{"hakstian_whalen"}), \code{"adf"},
#' \code{"adf_logistic"}, and \code{"likelihood"} need it, FIML giving
#' none of them. On complete data the standard error of
#' omega from \code{\link{fa_confirmatory}} therefore differs slightly
#' from the one here: by 0.4 percent for \code{t1_visual_perception},
#' \code{t2_cubes}, \code{t3_paper_form_board}, and \code{t4_lozenges}
#' in \code{holzinger_swineford}. When the default meets incomplete
#' items and one of those intervals, the call stops rather than deleting
#' cases unasked, and names the two choices: \code{missing = "listwise"},
#' which computes the interval from the complete rows alone, or an
#' interval FIML gives (\code{"ml"}, \code{"mlr"}, their
#' \code{_logistic} variants, or a bootstrap method). The
#' \code{aux} argument names auxiliary variables: columns of
#' \code{data} that are not part of the composite but are correlated
#' with the items or with the reasons values are missing. They enter as
#' \emph{saturated correlates} (Graham, 2003): correlated freely with
#' each other and with every item's residual, never loading on the
#' factor, so the measurement model is undisturbed while FIML uses
#' their information; a good auxiliary also makes MAR itself more
#' plausible (Collins, Schafer, & Kam, 2001). A case whose items are all
#' missing is kept when an auxiliary variable is observed on it: the
#' model estimates the items and the auxiliaries jointly, so such a case
#' informs the item moments through its auxiliaries, and setting it
#' aside would make the analysis conditional on having an observed item,
#' which biases the estimates when the items go missing together at
#' random given the auxiliaries. Supplying \code{aux}
#' without \code{missing} implies \code{missing = "fiml"}, whether or not
#' a value is missing; combining it with an explicit \code{missing =
#' "listwise"} is an error. Under
#' \code{missing = "fiml"} the denominators are estimated as follows:
#' with \code{"model_implied"} the fitted FIML model supplies the total
#' variance directly, and with \code{"observed"} the total variance is
#' the sum of the FIML estimate of the item covariance matrix (the EM
#' estimate under a saturated model over the items and any auxiliaries,
#' the matrix \code{\link{cov_fiml}} returns), which is the
#' estimate of \eqn{\mathrm{Var}(Y)} that uses the partially observed
#' rows; it is already on the maximum likelihood divisor, matching the
#' fitted loadings. The available intervals under \code{"fiml"} are
#' \code{"ml"}, \code{"mlr"} (Yuan & Bentler, 2000, robust), their
#' \code{_logistic} variants, and the bootstrap methods, which resample
#' rows including the partially observed ones and refit by FIML; the
#' complete-data closed forms, \code{"adf"}, and \code{"likelihood"}
#' are errors rather than silent fallbacks to listwise deletion.
#'
#' \code{missing = "two_stage"} estimates in two stages instead (Yuan &
#' Bentler, 2000; Savalei & Bentler, 2009): the means and covariances of
#' the items are first estimated by maximum likelihood under the saturated
#' normal model (the EM estimates, the same saturated estimate FIML's
#' observed denominator uses), and the single-factor model is then fit to
#' them by the complete data maximum likelihood discrepancy function. The
#' auxiliary variables named in \code{aux} enter the first stage only,
#' where they inform the estimated moments, and the measurement model is
#' fit as if they were absent, which is simpler than the saturated
#' correlates FIML needs. Two-stage estimation uses the cases FIML uses,
#' every case with an observed item or an observed auxiliary variable:
#' the first stage estimates the items and the auxiliaries jointly, so a
#' case observed on an auxiliary alone informs the item moments through
#' it, and setting such a case aside would condition the analysis on an
#' observed item and bias the estimates when the items went missing
#' together at random given the auxiliaries. \pkg{lavaan}'s two-stage
#' estimator reads such a case in its first stage but leaves it out of
#' its sample size, which would make its standard errors too large; the
#' fits go through \code{\link{fa_confirmatory}}, which computes them again with
#' every case the first stage read, and a message counts such cases (when
#' a fit cannot be computed again, a warning says that the standard error
#' is too large and that a bootstrap interval counts every case). The
#' estimates are consistent under MAR and only marginally less efficient
#' than FIML's (Savalei & Bentler, 2009), and
#' with incomplete nonnormal data Savalei and Falk (2014) found the
#' coverage of robust two-stage intervals better than that of robust FIML
#' intervals. With \code{"model_implied"} the \code{"ml"} interval
#' uses the normal theory two-stage standard error (Savalei & Bentler,
#' 2009, Eq. 6) and \code{"mlr"} the standard error robust to
#' nonnormality (Yuan & Bentler, 2000; lavaan's
#' \code{missing = "robust.two.stage"}); with \code{"observed"} the
#' numerator comes from the two-stage loadings and the denominator from the
#' EM covariance matrix, so the coefficient differs from its FIML
#' counterpart only through the loadings, and its interval is a bootstrap
#' that refits both stages on every replication. A delta method interval
#' for that coefficient can be derived under two-stage estimation, but its
#' behavior has not been studied, so it is not offered. With three items
#' the single-factor model is saturated and two-stage estimation gives the
#' FIML estimate and standard error. The intervals that FIML refuses are
#' refused here too, and \code{"likelihood"} for a further reason: no
#' likelihood of the observed data is maximized. Multiple imputation is a
#' different feature with a different interface and is out of scope here.
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
#' \strong{Improper solutions.} A single-factor fit can converge to an
#' improper solution, one with a negative error variance (a Heywood
#' case). McDonald (2004), following van Driel (1978), separates
#' improper estimates that reflect sampling error around a proper
#' structure, typically a negative error variance near zero in a small
#' sample or for an item that nearly duplicates the factor, from those
#' that reflect a misspecified or empirically underidentified model:
#' the first kind needs no change of model, while the second calls for
#' a different model rather than an estimate. One sample cannot always
#' tell the two apart, so the function draws the line at what the
#' coefficient can mean. When the error variances still sum to a
#' nonnegative number, the model implies a nonnegative error variance
#' for the composite and a model implied \eqn{\omega} between zero and
#' one, so the coefficient is returned with a warning. When they sum to less than zero, the model implies a
#' negative error variance for the composite and a model implied
#' \eqn{\omega} above one, which no reliability coefficient can take;
#' sampling error near a zero error variance does not produce this
#' unless every item is nearly error free, and it is the signature of a
#' solution that has run away, one error variance falling without bound
#' as its loading grows. The function then stops with a message naming
#' the item. So it does for three items whose correlations have a
#' negative product: a single factor implies that the product of the
#' three correlations equals the square of the product of the three
#' loadings, so no real loadings reproduce them and the maximum
#' likelihood estimate does not exist.
#' Whatever the sum, a fit whose likelihood has no maximum, an error
#' variance running off below zero while the likelihood is flat or still
#' rising (checked by fitting the model again with the variance held
#' further out, as \code{\link{fa_bifactor}} describes), stops the
#' function with an error of class \code{dmar_no_maximum} that names the
#' item and prints the \pkg{lavaan} code of the single-factor fit, where
#' \code{bounds = "pos.var"} or a different model are the choices.
#'
#' A bootstrap meets improper solutions more often than the sample does,
#' since a resample can push an error variance that is near zero in the
#' sample below it, which is a standardized loading past one (the fits
#' are on the items' own metric, where a loading above one is
#' routine). The argument \code{improper} governs the replications whose
#' solution is improper but still yields a coefficient: \code{"keep"},
#' the default, uses their estimates like any other, and \code{"drop"}
#' sets them aside. Keeping is the default because setting them aside
#' conditions the bootstrap distribution on admissibility and so trims
#' its upper tail, the replications in which an error variance fell
#' below zero; Yuan and Hayashi (2003) make the parallel argument for
#' test statistics, counting a replication that fails because its
#' resample lies far from the model rather than discarding it. A
#' replication whose solution yields no coefficient is set aside under
#' either setting. When its error variances sum below zero, the solution
#' has run away and implies an \eqn{\omega} above one, beyond the top of
#' the distribution like the improper replications kept. When three
#' items have no single-factor solution, the side depends on the scale
#' the resample is read on: on the estimator's own scale it continues
#' the top of the distribution, since as a resample approaches that
#' boundary one loading, and \eqn{\omega} with it, grows without bound,
#' while its correlations, one of which is negative or zero, point to
#' low reliability. The argument of Yuan and Hayashi does not settle it,
#' since their extreme value belongs to a statistic that measures misfit
#' and so has a direction, so such a replication is counted, on its own
#' as well as among the improper ones, rather than placed on either
#' side. A replication whose likelihood has no maximum is set aside as
#' well, and counted among those that raised an error, under a message
#' that names the item whose error variance ran off; the run goes on.
#' Every such replication is counted in the \code{bootstrap}
#' attribute and the printed line, and when the replications set aside
#' and the improper ones kept together reach the tail area of the
#' interval, \eqn{(1 - C)/2} of \code{B}, a caution says that the limits
#' depend on how they were treated. When some of them had no
#' single-factor solution, the caution adds that the lower limit is
#' likely too high: with three weakly related items the estimator itself
#' runs high near that boundary, and the resamples beyond it are those
#' whose correlations point to low reliability. When about half of the
#' replications are improper, neither setting gives an interval worth
#' reporting, and the counts are what the result can state.
#'
#' \strong{Comparison with other packages.} The \pkg{psych} package
#' provides \code{\link[psych]{omega}}, which fits a Schmid-Leiman
#' hierarchical factor model and reports several variants of \eqn{\omega}
#' (\eqn{\omega_t}, \eqn{\omega_h}) alongside extensive psychometric
#' diagnostics. \code{reliability_omega} in \pkg{DMAR} differs in
#' emphasis: it implements McDonald's \eqn{\omega} from a single-factor
#' (congeneric) model and accompanies the point estimate with a
#' confidence interval drawn from the methods compared in Kelley and
#' Pornprasertmanit (2016). The same denominator distinction appears in
#' semTools' \code{compRelSEM()} as its \code{obs.var} argument.
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
#'   otherwise on the EM estimate of the covariance matrix. In
#'   \code{holzinger_swineford}, \code{t1_visual_perception},
#'   \code{t2_cubes}, and \code{t4_lozenges}, each centered on the
#'   pupil's mean of the three, are ipsative items of this kind.
#'   SPSS user-defined missing values, which
#'   \code{haven::read_sav(user_na = TRUE)} keeps as values of a
#'   \code{haven_labelled_spss} column, are missing, as SPSS treats
#'   them.
#' @param S A symmetric covariance matrix among the items. If supplied,
#'   \code{N} must also be supplied; methods that require raw data
#'   (\code{"mlr*"}, \code{"adf*"}, \code{"bootstrap_*"},
#'   \code{"percentile"}, \code{"bca"}) are then unavailable, as are
#'   \code{missing = "fiml"} and \code{missing = "two_stage"} named
#'   explicitly, and \code{aux}, since a covariance matrix has no
#'   incomplete cases (the default of \code{missing} analyzes it as it
#'   stands).
#' @param N Total sample size; required when \code{S} is supplied.
#' @param ci_method Method for constructing the confidence interval. See
#'   \emph{Details}. When not supplied (or passed as \code{NULL}): for
#'   \code{denominator = "model_implied"} the default is
#'   \code{"mlr_logistic"} with raw data, the robust maximum likelihood
#'   standard error with the limits formed on the logit scale and
#'   transformed back (Browne, 1982; Kelley, 2007a), whatever the
#'   treatment of missing values, except under \code{wald_df =
#'   "missing_information"}, where it is \code{"mlr"}, the same standard
#'   error with the limits formed on the scale of omega, the interval that
#'   rule is calibrated for; with
#'   covariance input, which the robust standard error cannot use, it is
#'   \code{"ml_logistic"}, the maximum likelihood standard error with the
#'   limits formed the same way (\code{"ml"} under \code{wald_df =
#'   "missing_information"});
#'   for
#'   \code{denominator = "observed"} (robust omega) no interval is
#'   computed by default, because its interval is bootstrap based and a
#'   bootstrap is never run unless requested. Ask for
#'   \code{"percentile"} or \code{"bca"} to obtain the recommended
#'   interval. The closed forms, the Wald intervals, and the profile
#'   likelihood apply to \code{denominator = "model_implied"} only; a
#'   combination that does not apply stops with a message naming the
#'   denominator that does (see \emph{Details}).
#' @param denominator How the total variance in the denominator of
#'   \eqn{\omega} is estimated: \code{"observed"} (default; robust
#'   omega) uses the variance of the composite estimated directly from
#'   the data; \code{"model_implied"} uses the total variance reproduced
#'   by the fitted single-factor model. See \emph{Details}. With
#'   \code{"observed"}, \code{ci_method} must be a bootstrap method or
#'   \code{"none"}.
#' @param missing How incomplete rows of \code{data} are handled:
#'   \code{"fiml"} (the default; full information
#'   maximum likelihood, using every case with at least one observed item
#'   or observed auxiliary variable), \code{"listwise"} (complete-case
#'   analysis), or \code{"two_stage"}
#'   (two-stage maximum likelihood, using the same cases as FIML). The default analysis is FIML whenever an item value
#'   is missing and the complete data analysis otherwise, unlike
#'   \code{\link{fa_confirmatory}}, which fits complete data by FIML too:
#'   data in which no value is missing, data complete but for rows with
#'   every item missing (which both treatments drop), and \code{S} are
#'   analyzed as they stand, since the
#'   complete-data intervals need that analysis; incomplete items are
#'   analyzed by FIML, or, when \code{ci_method} names an interval FIML cannot give, the call
#'   stops and names the choices. See \emph{Details}. With
#'   \code{missing = "listwise"}, when fewer than
#'   two rows are complete, listwise deletion leaves nothing to analyze,
#'   and the function stops before any fit, says so, and names the
#'   treatment that uses the incomplete rows. The single-factor model is fit by
#'   maximum likelihood, which needs the covariance matrix positive
#'   definite, and \eqn{N} complete rows give a matrix of rank at most
#'   \eqn{N - 1}, so with no more complete rows than items the function
#'   stops the same way;
#'   \code{\link{reliability_alpha}}, which needs only the sums of the
#'   matrix, is computed from two or more complete rows. Full information
#'   maximum likelihood does not impose an ordinary complete-case count
#'   rule. Observed-denominator omega and two-stage estimation require
#'   unrestricted moments, checked for affine collapse or numerical
#'   instability as described in \code{\link{cov_fiml}}. EM convergence
#'   alone does not establish existence or uniqueness. Model-implied omega
#'   can retain a proper identified constrained FIML fit and H0-based
#'   uncertainty while withholding unavailable H1 diagnostics and residual
#'   correlations. The result states this qualification. Bootstrap
#'   resamples retain a distinct-case estimability screen as a separate
#'   policy, not a universal likelihood-existence theorem. Listwise
#'   complete rows are counted as distinct rows. Two columns of one name
#'   are refused before any of
#'   this, under every treatment, since the items are read by name.
#' @param aux Optional character vector naming auxiliary variable
#'   columns of \code{data}: entered as saturated correlates under full
#'   information maximum likelihood, and into the first stage under
#'   \code{missing = "two_stage"}. Supplying \code{aux} without
#'   \code{missing} implies \code{missing = "fiml"}. See \emph{Details}.
#' @param conf_level Confidence level for the interval. Defaults to
#'   \code{0.95}.
#' @param B Number of bootstrap replications when a bootstrap method is
#'   selected, a single whole number. Defaults to \code{10000}. The fewest
#'   accepted depends on the interval method and the level: at the default
#'   \code{conf_level}, 400 for \code{"percentile"} (10 replications beyond
#'   each limit), 2000 for \code{"bca"} (50 beyond each nominal limit,
#'   because the adjusted limits sit farther into the tails), and 100 for
#'   \code{"bootstrap_se"} and \code{"bootstrap_se_logistic"} (a standard
#'   deviation of the replications); at 0.99 the first two are 2000 and
#'   10000. A smaller \code{B} is refused before anything is computed
#'   (class \code{dmar_count_invalid}). The count is held to the number requested.
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
#'   for a reproducible interval, which is set for the call only, the
#'   caller's random number state being restored when it returns.
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
#' @param improper How a bootstrap replication whose single-factor
#'   solution is improper (a negative error variance, a Heywood case) is
#'   treated when the solution still yields a coefficient:
#'   \code{"keep"} (the default) uses its estimate like any other, and
#'   \code{"drop"} sets it aside. Either way such replications are
#'   counted in the \code{bootstrap} attribute and the printed line, and
#'   a replication whose solution yields no coefficient is set aside
#'   under both (see \emph{Details}). With \code{"drop"} an improper
#'   point estimate stops the call, since an interval built from proper
#'   replications only would not describe it. Ignored unless a
#'   bootstrap method is selected.
#' @param wald_df The reference distribution of a Wald interval:
#'   \code{Inf} (the default), the normal distribution; a single positive number,
#'   the degrees of freedom of a \emph{t} distribution; or
#'   \code{"missing_information"}, von Hippel's (2016) degrees of freedom
#'   \eqn{\max\{3, N(1 - \gamma) - 1\}}{max(3, N(1 - gamma) - 1)} from
#'   omega's fraction of missing information \eqn{\gamma} (\eqn{N - 1}
#'   with complete data). A number
#'   applies to the Wald intervals of the model implied coefficient
#'   (\code{"ml"}, \code{"mlr"}, \code{"adf"}, and their \code{_logistic}
#'   variants), the keyword to \code{"ml"} and \code{"mlr"} only, without
#'   two-stage estimation, the intervals its calibration covers; anything
#'   else is an error that says why. See \emph{Details}.
#' @param progress \code{NULL} (the default), or a function of two
#'   arguments that the bootstrap calls as it runs, so that an application
#'   can show how far a long bootstrap has gone and let its user stop it (a
#'   Shiny app, for example). It is called as \code{progress(done, total)},
#'   with the number of replications attempted so far and \code{B}: after
#'   each replication when \code{ncpus = 1}, and with more processes after
#'   each block of replications they complete together, about a fiftieth
#'   of \code{B} (\code{ncpus * ceiling(B / (50 * ncpus))} replications,
#'   the last block possibly fewer), so that \code{done} advances by the
#'   block. A BCa interval takes its acceleration from a regression on
#'   the resampling frequencies when at least as many replications return
#'   a value as there are cases, and from the jackknife, one leave-one-out
#'   refit per case, when fewer do or when that regression fails or gives
#'   influence values that are not all finite; the jackknife runs at most
#'   once in a call, and not at all when no interval is formed or every
#'   returned replicate has the same value. While it runs, the callback is
#'   called again after each refit (with more processes, after each block
#'   of \code{ncpus * ceiling(n / (50 * ncpus))} refits, \code{n} the
#'   number of cases resampled), with \code{done} and \code{total} both
#'   \code{B}, so it can be called more than \code{B} times, and
#'   \code{done} equal to \code{total} does not mean that the bootstrap
#'   has finished. If the callback returns \code{FALSE}, the bootstrap
#'   stops before the next replication or block with an error of class
#'   \code{dmar_bootstrap_cancelled}, whose elements \code{done} and
#'   \code{total} say how far it got (see \code{\link{dmar_metadata}});
#'   any other value is ignored. An error inside the callback stops the
#'   call with an error of class \code{dmar_progress_error}, whose message
#'   names the callback and which carries \code{done}, \code{total}, and,
#'   as \code{parent}, the callback's own condition with its class intact,
#'   so that an application can tell a failing callback from a failing
#'   analysis. The callback runs outside the random number stream, so a
#'   given \code{seed} returns the same result with or without it, and a
#'   callback that draws random numbers leaves the caller's stream as it
#'   found it. With \code{NULL} nothing is called. Ignored unless a
#'   bootstrap method is selected.
#'
#' @return A \code{data.frame} with columns \code{term} and \code{value}
#' and rows
#' \code{"estimate"} (sample coefficient \eqn{\omega}),
#' \code{"se"} (standard error on the coefficient scale, \code{NA} for
#' methods that do not produce one; for the transformation-based
#' intervals \code{"fisher"}, \code{"bonett"}, and
#' \code{"hakstian_whalen"} it is the delta method back-transform of the
#' transformation-scale standard error),
#' \code{"se_transformed"} (only for those transformation-based
#' intervals: the standard error on the transformation scale, with the
#' scale named in the attribute \code{se_transform_scale}:
#' \code{"fisher_z"}, \code{"log(1-alpha)"}, or \code{"cube_root"}),
#' \code{"ci_lower"} and \code{"ci_upper"} (each clamped into
#' [0, 1], so the limits stay in order),
#' \code{"conf_level"}, \code{"N"} (the cases the analysis used: the
#' complete cases under listwise deletion, and under \code{"fiml"} and
#' \code{"two_stage"} every case with at least one observed item or
#' observed auxiliary variable),
#' \code{"N_complete"}
#' (the complete cases; equal to \code{"N"} under listwise deletion and
#' for covariance input), and \code{"n_items"} (number of items). Attributes
#' \code{coefficient} (\code{"omega"}), \code{ci_method},
#' \code{denominator}, \code{missing}, \code{wald_df} (the degrees of
#' freedom of the Wald interval's \emph{t} reference, \code{Inf} for the
#' normal one and for intervals that are not Wald intervals; under
#' \code{wald_df = "missing_information"} the attribute
#' \code{wald_df_table} holds the \code{\link{df_missing}} table they
#' came from), and (when supplied) \code{aux} record the computation; bootstrap calls also
#' record \code{B}. The
#' attribute \code{residual_cor} carries the local dependence summary of
#' a maximum likelihood fit of the single-factor model to the same cases,
#' whatever the \code{ci_method} (the weighted least squares fit behind
#' the ADF intervals has residuals of its own, so it is not used for the
#' summary): a list with \code{matrix} (the observed minus the model
#' implied item correlations, named by item), \code{largest} (the
#' residual correlation largest in absolute value, with its sign),
#' \code{items} (the pair it belongs to), and \code{saturated}. With
#' three items the single-factor model is saturated (zero degrees of
#' freedom) and a proper solution reproduces every correlation, so its
#' residuals are zero and carry no information about local dependence:
#' \code{saturated} is then \code{TRUE}, \code{largest} is \code{NA},
#' and \code{items} is \code{NULL}, rather than naming whichever pair
#' carries the largest rounding error. The printed table names the
#' denominator and the largest residual correlation (or says that the
#' model is saturated) in two lines under the table. A bootstrap call also carries the attribute
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
#' single-factor solution was improper, with a negative error
#' variance, or with three items had no single-factor solution at all:
#' those set aside, which are among the \code{non_finite} ones, and
#' those whose estimate the interval uses, which are among the
#' \code{returned} ones), \code{no_solution} (the replications of three
#' items set aside because their correlations admit no single-factor
#' solution, a subset of \code{improper}),
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
#'   analysis of covariance structures.
#'   \emph{British Journal of Mathematical and Statistical Psychology, 37},
#'   62--83.
#'
#' Satorra, A., & Bentler, P. M. (1994). Corrections to test statistics
#' and standard errors in covariance structure analysis. In A. von Eye &
#' C. C. Clogg (Eds.), \emph{Latent variables analysis: Applications for
#' developmental research} (pp. 399--419). Thousand Oaks, CA: Sage.
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
#' van Driel, O. P. (1978). On various causes of improper solutions in
#'   maximum likelihood factor analysis. \emph{Psychometrika, 43},
#'   225--243.
#'
#' von Hippel, P. T. (2016). New confidence intervals and bias comparisons
#'   show that maximum likelihood can beat multiple imputation in small
#'   samples. \emph{Structural Equation Modeling, 23}(3), 422--437.
#'   \doi{10.1080/10705511.2015.1047931}
#'
#' Yuan, K.-H., & Bentler, P. M. (2000). Three likelihood-based methods
#' for mean and covariance structure analysis with nonnormal missing
#' data. \emph{Sociological Methodology, 30}, 165--200.
#'
#' Collins, L. M., Schafer, J. L., & Kam, C.-M. (2001). A comparison of
#'   inclusive and restrictive strategies in modern missing data
#'   procedures. \emph{Psychological Methods, 6}, 330--351.
#'   \doi{10.1037/1082-989X.6.4.330}
#'
#' Cronbach, L. J. (1951). Coefficient alpha and the internal structure
#'   of tests. \emph{Psychometrika, 16}(3), 297--334.
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
#' Green, S. B., & Yang, Y. (2009). Reliability of summed item scores
#'   using structural equation modeling: An alternative to coefficient
#'   alpha. \emph{Psychometrika, 74}, 155--167.
#'   \doi{10.1007/s11336-008-9099-3}
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
#' Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). \emph{Designing
#'   experiments and analyzing data: A model comparison perspective}
#'   (4th ed.). Routledge.
#'
#' McDonald, R. P. (1999). \emph{Test theory: A unified treatment}.
#'   Mahwah, NJ: Lawrence Erlbaum Associates.
#'
#' McDonald, R. P. (2004). Respecifying improper structures.
#'   \emph{Structural Equation Modeling, 11}(2), 194--209.
#'
#' McDonald, R. P. (2011). Measuring latent quantities.
#'   \emph{Psychometrika, 76}, 511--536. \doi{10.1007/s11336-011-9223-7}
#'
#' Raykov, T. (2002). Analytic estimation of standard error and confidence
#'   interval for scale reliability.
#'   \emph{Multivariate Behavioral Research, 37}, 89--103.
#'   \doi{10.1207/S15327906MBR3701_04}
#'
#' Satorra, A., & Bentler, P. M. (2001). A scaled difference chi-square test
#'   statistic for moment structure analysis. \emph{Psychometrika, 66}(4),
#'   507--514. \doi{10.1007/BF02296192}
#'
#' Terry, L. J., & Kelley, K. (2012). Sample size planning for composite
#'   reliability coefficients: Accuracy in parameter estimation via narrow
#'   confidence intervals. \emph{British Journal of Mathematical and
#'   Statistical Psychology, 65}, 371--401.
#'   \doi{10.1111/j.2044-8317.2011.02030.x}
#'
#' Yuan, K.-H., & Hayashi, K. (2003). Bootstrap approach to inference
#'   and power analysis based on three test statistics for covariance
#'   structure models. \emph{British Journal of Mathematical and
#'   Statistical Psychology, 56}, 93--110.
#'
#' Zhang, Z., & Yuan, K.-H. (2016). Robust coefficients alpha and omega
#'   and confidence intervals with outlying observations and missing
#'   data: Methods and software.
#'   \emph{Educational and Psychological Measurement, 76}, 387--411.
#'   \doi{10.1177/0013164415594658}
#'
#' Zinbarg, R. E., Revelle, W., Yovel, I., & Li, W. (2005). Cronbach's
#'   \eqn{\alpha}, Revelle's \eqn{\beta}, and McDonald's \eqn{\omega_H}:
#'   Their relations with each other and two alternative
#'   conceptualizations of reliability. \emph{Psychometrika, 70},
#'   123--133. \doi{10.1007/s11336-003-0974-7}
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso
#' \code{\link{reliability}} (general wrapper),
#' \code{\link{reliability_omega_categorical}} (categorical omega for ordered
#' items),
#' \code{\link{reliability_alpha}},
#' \code{\link{cfa_1}} (single-factor CFA used internally),
#' \code{\link[psych]{omega}}.
#'
#' @examples
#' set.seed(113)
#' J <- 6
#' loadings <- seq(0.5, 0.8, length.out = J)
#' eta <- rnorm(200)
#' errors <- matrix(rnorm(200 * J), 200, J) %*% diag(sqrt(1 - loadings^2))
#' items <- sweep(matrix(rep(eta, J), 200, J), 2, loadings, `*`) + errors
#' colnames(items) <- paste0("y", seq_len(J))
#'
#' # Default: robust omega, point estimate only (no bootstrap is run
#' # unless requested; a message names the call that produces the
#' # recommended interval).
#' reliability_omega(data = items)
#'
#' # The closed-form standard errors are derived under the model implied
#' # ratio, so with the observed denominator the interval comes from a
#' # bootstrap, which refits the single factor model once per
#' # replication. B = 10 keeps the example quick and is below the 400 a
#' # reported interval needs; the option lets it run and the warning says
#' # so. A reported interval deserves the default B = 10000, an accuracy
#' # choice rather than a formality (see Details). The percentile
#' # interval is what Kelley and Pornprasertmanit (2016) recommend for
#' # this coefficient, with ci_method = "bca" as the alternative, and
#' # seed makes the interval reproducible.
#' old <- options(dmar.bootstrap_exploratory = TRUE)
#' reliability_omega(data = items, ci_method = "percentile", B = 10,
#'                   seed = 113)
#' options(old)
#'
#' # Model implied denominator with its closed-form robust ML interval,
#' # whose limits are formed on the logit scale; ci_method = "mlr" gives
#' # limits formed on the scale of omega.
#' reliability_omega(data = items, denominator = "model_implied")
#'
#' # Two further routes into the same coefficient. The first works from
#' # the summary statistics a paper reports, a covariance matrix and its
#' # sample size, with the maximum likelihood interval that summary
#' # input supports.
#' reliability_omega(S = cov(items), N = 200,
#'                   denominator = "model_implied", ci_method = "ml")
#'
#' # How much would a t reference lengthen that interval? On 30 degrees of
#' # freedom, an arbitrary small value chosen to show the effect, the
#' # limits move out a little; the line under the table says so.
#' reliability_omega(S = cov(items), N = 200,
#'                   denominator = "model_implied", ci_method = "ml",
#'                   wald_df = 30)
#'
#' # The second is full information maximum likelihood with an auxiliary
#' # variable, where missingness on y2 depends on an auxiliary z (missing
#' # at random given z). Supplying aux implies missing = "fiml"; the
#' # N_complete row shows how many rows listwise deletion would have
#' # kept. With values missing, the interval can also refer to a t
#' # distribution on degrees of freedom taken from omega's own fraction of
#' # missing information, here in the fit that carries z (the calibration
#' # study behind the rule had no auxiliary variable, so this case rests
#' # on analogy; see Details); the line under the table gives N, the FMI,
#' # and the degrees of freedom that result. The rule is calibrated for
#' # limits on the scale of omega, so the interval here is the "mlr" one.
#' z <- eta + rnorm(200, sd = 0.5)
#' d <- data.frame(items, z = z)
#' d$y2[runif(200) < plogis(-1 + 1.5 * as.numeric(scale(z)))] <- NA
#' reliability_omega(data = d, aux = "z", denominator = "model_implied",
#'                   wald_df = "missing_information")
#'
#' # Two-stage estimation instead: z enters only the first stage, where the
#' # means and covariances are estimated, and the interval uses the robust
#' # two-stage standard error, its limits again formed on the logit scale.
#' reliability_omega(data = d, aux = "z", missing = "two_stage",
#'                   denominator = "model_implied")
#'
#' @keywords htest multivariate
#' @family reliability
#'
#' @export

reliability_omega <- function(data = NULL, S = NULL, N = NULL,
                              ci_method = c("mlr_logistic", "ml_logistic",
                                            "mlr", "ml",
                                            "likelihood",
                                            "adf", "adf_logistic",
                                            "feldt", "fisher", "bonett",
                                            "hakstian_whalen",
                                            "bootstrap_se",
                                            "bootstrap_se_logistic",
                                            "percentile", "bca", "none"),
                              denominator = c("observed", "model_implied"),
                              missing = c("fiml", "listwise", "two_stage"),
                              aux = NULL,
                              conf_level = 0.95,
                              B = 10000,
                              seed = NULL,
                              ncpus = 1L,
                              improper = c("keep", "drop"),
                              wald_df = Inf,
                              progress = NULL) {
  .dmar_check_conf_level(conf_level)
  # An explicit ci_method = NULL means the default, as it does in
  # reliability() and as it did through DMAR 1.0.x, where match.arg(NULL)
  # gave the first choice, the default of the time. It therefore counts as
  # not supplied: it takes the default path below, with its choice under
  # wald_df = "missing_information" and for a covariance matrix, and
  # robust omega's point estimate with advice.
  ci_supplied <- !missing(ci_method) && !is.null(ci_method)
  ci_method <- match.arg(ci_method)
  denominator <- match.arg(denominator)
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
  # Robust omega has no closed-form interval. This is checked before the
  # missing data treatment is resolved, so that the refusal of an interval
  # FIML cannot give (below) never names a remedy that meets this refusal
  # next: with the observed denominator such an interval is refused here
  # whatever the treatment. An unsupplied ci_method becomes "none" with
  # this denominator (advise_interval, below).
  bootstrap_methods <- c("bootstrap_se", "bootstrap_se_logistic",
                         "percentile", "bca")
  if (denominator == "observed" && ci_supplied &&
      !ci_method %in% c(bootstrap_methods, "none")) {
    .dmar_refuse_argument(c(
      "ci_method = \"", ci_method, "\" is available only for ",
      "denominator = \"model_implied\": the closed forms and the ",
      "profile likelihood are defined under the model implied total ",
      "variance. With denominator = \"observed\", use a bootstrap ",
      "method (\"percentile\", \"bca\", \"bootstrap_se\", ",
      "\"bootstrap_se_logistic\") or ci_method = \"none\"."), arg = NULL)
  }
  # The default is FIML when a value is missing (supplying aux implies it
  # in any case, unless two-stage estimation was asked for, so the common
  # case needs one argument, not two); complete data and a covariance
  # matrix are analyzed as they stand, as DMAR 1.0.x analyzed them, and an
  # interval FIML cannot give stops the call rather than deleting cases
  # unasked (.relia_resolve_missing()). Such an interval reaches the
  # resolution only with the model implied denominator (above), so the
  # choices named under FIML are that denominator's.
  missing <- .relia_resolve_missing(
    missing, missing_supplied, data, aux,
    unavailable = if (ci_supplied &&
                      ci_method %in% .RELIA_COMPLETE_DATA_METHODS) {
      paste0("ci_method = \"", ci_method, "\"")
    },
    fiml_options = paste0("keep FIML with ci_method = \"ml\", \"mlr\", ",
                          "their _logistic variants, or a bootstrap ",
                          "method"))
  if (missing != "listwise" && is.null(data)) {
    .dmar_refuse_argument(c(
      "missing = \"", missing, "\"",
      if (missing == "fiml") " (and 'aux')", " requires raw 'data': a ",
      "covariance matrix has no incomplete cases for ",
      .relia_incomplete_label(missing), " to use."), arg = "data")
  }
  if (!is.numeric(conf_level) || length(conf_level) != 1L ||
      conf_level <= 0 || conf_level >= 1) {
    .dmar_refuse_argument(
      "'conf_level' must be a single number in (0, 1).", arg = "conf_level")
  }
  # A t reference distribution changes a Wald interval only: the delta
  # method intervals of the model implied coefficient. Robust omega's
  # interval is a bootstrap (and its default is no interval at all). The
  # keyword "missing_information" is narrower: the calibration study
  # behind it covers the intervals on the scale of omega ("ml", "mlr"),
  # with complete data or by FIML, and anything else is refused here with
  # its reason (.wald_df_scope(), R/missing_information_internals.R).
  wald_df <- .wald_df_check(wald_df, "reliability_omega")
  # The default interval of the model implied coefficient is "mlr_logistic"
  # since DMAR 1.1.0 ("mlr" through 1.0.x): the same robust standard error,
  # with the limits formed on the logit scale and transformed back, as
  # fa_confirmatory() forms its omega limits. The keyword
  # "missing_information" is calibrated for the intervals on the scale of
  # omega only (.wald_df_scope_omega()), so under it the default stays
  # "mlr", and the limits of that path are those of earlier versions. A
  # covariance matrix, which the robust standard error cannot use, takes
  # the maximum likelihood standard error of the same form ("ml_logistic",
  # or "ml" under the keyword), as reliability_H() and reliability() do;
  # through 1.0.x such a call stopped. Robust omega's default is no
  # interval (advise_interval, below).
  if (!ci_supplied) {
    ci_method <- if (identical(wald_df, "missing_information")) {
      "mlr"
    } else {
      "mlr_logistic"
    }
    if (is.null(data) && !is.null(S)) {
      ci_method <- sub("^mlr", "ml", ci_method)
    }
  }
  .wald_df_scope(wald_df, "reliability_omega",
                 list(denominator = denominator, ci_method = ci_method,
                      missing = missing))
  wald <- denominator == "model_implied" && ci_method %in% .RELIA_WALD_METHODS
  .wald_df_require_wald(wald_df, wald, if (denominator == "observed") {
    paste0("and robust omega (denominator = \"observed\") has a ",
           "bootstrap interval only; wald_df applies with denominator = ",
           "\"model_implied\" and ci_method \"ml\", \"mlr\", \"adf\", ",
           "or their _logistic variants.")
  } else {
    paste0("and ci_method = \"", ci_method, "\" is not one; wald_df ",
           "applies to ci_method \"ml\", \"mlr\", \"adf\", and their ",
           "_logistic variants.")
  })

  if (missing != "listwise") {
    # Only intervals whose sampling theory can respect the missingness
    # pattern are offered; nothing silently falls back to listwise
    # deletion. The delta method ("ml"/"mlr") standard errors come from
    # the FIML information matrix, or under two-stage estimation from the
    # two-stage sandwich (normal theory for "ml", robust to nonnormality
    # for "mlr"), and the bootstrap resamples the partially observed rows
    # and refits by the same method. The remaining methods are
    # complete-data constructions.
    if (ci_method %in% c("feldt", "fisher", "bonett", "hakstian_whalen")) {
      .dmar_refuse_argument(c(
        "ci_method = \"", ci_method, "\" is a complete-data closed ",
        "form: its standard error treats every case as fully ",
        "observed, which overstates the information in partially ",
        "observed rows. With missing = \"", missing, "\" use \"ml\" ",
        "(or \"mlr\"), whose standard error comes from ",
        if (missing == "fiml") "the FIML information matrix"
        else "the two-stage sandwich",
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

  inp <- .relia_resolve_inputs(data = data, S = S, N = N, allow_S = TRUE,
                               missing = missing, aux = aux,
                               fun = "reliability_omega",
                               single_factor = "omega",
                               needs_saturated = denominator == "observed" || missing == "two_stage" ||
                                 identical(wald_df, "missing_information"))

  # The single-factor congeneric model behind omega is not identified
  # with two items (four free parameters against three observed
  # moments), so a two-item fit returns whichever value the optimizer
  # wanders to. Refuse it rather than report solver noise.
  if (inp$J < 3L) {
    .dmar_refuse_argument(c(
      "Coefficient omega needs three or more items: the single-factor ",
      "congeneric model is not identified with two. For a two-item ",
      "scale use reliability_alpha(), whose tau-equivalent ",
      "model is just identified."), arg = NULL)
  }

  # No closed-form standard error exists for the observed denominator,
  # and the bootstrap is never run unless the user asks for it. Report
  # the point estimate and say exactly how to request the recommended
  # interval. The message waits until the point estimate is in hand, so a
  # call that is about to stop (on its input or on its fit) does not
  # first receive advice.
  advise_interval <- denominator == "observed" && !ci_supplied
  if (advise_interval) ci_method <- "none"

  if (ci_method %in% bootstrap_methods) {
    .check_boot_B(B, ci_method, conf_level)
  }

  needs_raw <- ci_method %in% c("mlr", "mlr_logistic",
                                "adf", "adf_logistic",
                                bootstrap_methods)
  if (needs_raw && is.null(inp$data)) {
    .dmar_refuse_argument(c(
      "ci_method = \"", ci_method, "\" requires raw 'data'; ",
      "a covariance matrix is not sufficient."), arg = "data")
  }

  if (denominator == "observed") {
    # Under missing = "fiml" both pieces of robust omega come from full
    # information maximum likelihood (.omega_observed_fiml_parts()); with
    # complete data the loadings come from the complete cases and the
    # total variance is their observed composite variance on the same
    # maximum likelihood divisor (.omega_observed_parts()). Either way
    # the single-factor model is fit once, by ML, and that fit also
    # supplies the residual correlations below (through DMAR 1.0.1 the
    # estimate came from an MLR fit with the same loadings and the
    # residual correlations from a second, ML, fit).
    point_fn <- if (missing == "fiml") {
      function(d) {
        p <- .omega_observed_fiml_parts(d, items = inp$items, aux = inp$aux)
        .omega_point_value(p$fit, p$estimate)
      }
    } else if (missing == "two_stage") {
      # A replication needs the estimates alone, so it fits the second
      # stage to the EM moments without the two-stage standard errors and
      # test (.omega_fit_cfa_ts_estimates()).
      function(d) {
        p <- .omega_observed_ts_parts(d, items = inp$items, aux = inp$aux,
                                      estimates_only = TRUE)
        .omega_point_value(p$fit, p$estimate)
      }
    } else {
      .omega_observed_from_data
    }
    # A replication in which an item (or auxiliary variable) takes a
    # single value is refused before the fit, as the data are
    # (.relia_checked_point()).
    point_fn <- .relia_checked_point(
      point_fn, missing = missing,
      items = if (missing == "listwise") NULL else inp$items,
      aux = inp$aux, fitted = TRUE, composite = TRUE,
      needs_saturated = missing != "listwise")

    parts <- if (missing == "fiml") {
      .omega_observed_fiml_parts(inp$data, items = inp$items, aux = inp$aux)
    } else if (missing == "two_stage") {
      # Under two-stage estimation the loadings come from the model fit to
      # the stage one (EM) moments and the total variance from those
      # moments themselves, both read off the one fit
      # (.omega_observed_ts_parts()).
      .omega_observed_ts_parts(inp$data, items = inp$items, aux = inp$aux)
    } else if (!is.null(inp$data)) {
      .omega_observed_parts(inp$data)
    } else {
      fit <- .omega_fit_cfa(S = inp$S, N = inp$N,
                            equal_loadings = FALSE,
                            estimator = "ML", se = "none",
                            missing = "listwise")
      # Rescale the supplied covariance total to the fit's N-divisor
      # (maximum likelihood) metric so the ratio is internally consistent.
      list(estimate = if (fit$converged) {
        (sum(fit$loadings))^2 / (sum(inp$S) * (inp$N - 1) / inp$N)
      } else {
        NA_real_
      }, fit = fit)
    }
    # An improper solution is refused, warned about, or (under
    # improper = "drop" with a bootstrap) refused, before anything else
    # (.omega_point_gate()); a fit that did not converge stops below.
    .omega_point_gate(parts$fit, improper = improper,
                      bootstrap = ci_method %in% bootstrap_methods)
    estimate <- parts$estimate
    if (is.na(estimate)) {
      .dmar_refuse(
        .relia_fit_failure(parts$fit, missing = missing), "computation_failed")
    }
    if (advise_interval) {
      if (is.null(data)) {
        message("Robust omega is reported without a confidence interval: ",
                "its interval is bootstrap based, which requires raw data ",
                "rather than a covariance matrix.")
      } else {
        message("Robust omega is reported without a confidence interval ",
                "by default because its interval is bootstrap based. ",
                "Request it with ci_method = \"percentile\" (or \"bca\"); ",
                "B = 10000 replications is the default when you do.")
      }
    }

    ci <- if (ci_method == "none") {
      list(se = NA_real_, lower = NA_real_, upper = NA_real_)
    } else {
      .bootstrap_ci(
        fun = "reliability_omega",
        data = inp$data, point_fn = point_fn,
        B = B, conf_level = conf_level, kind = ci_method, seed = seed,
        ncpus = ncpus, improper = improper, progress = progress
      )
    }

    out <- .relia_result(
      estimate = estimate, se = ci$se,
      lower = ci$lower, upper = ci$upper,
      conf_level = conf_level, N = inp$N, J = inp$J,
      coefficient = "omega",
      ci_method = ci$ci_method %||% ci_method,
      bootstrap = ci$bootstrap,
      B = if (ci_method %in% bootstrap_methods) B else NA_integer_,
      N_complete = inp$N_complete,
      missing = missing,
      aux = inp$aux
    )
    attr(out, "denominator") <- denominator
    attr(out, "wald_df") <- wald_df
    # The local dependence summary comes from the maximum likelihood fit
    # behind the estimate; it runs once, never inside the bootstrap.
    residual_cor <- .omega_ml_residual_cor(inp, missing, fit = parts$fit,
                                           estimator = "ML")
    if (!is.null(residual_cor)) attr(out, "residual_cor") <- residual_cor
    return(out)
  }

  # Map ci_method to the (estimator, se) lavaan options for the fit that
  # produces the point estimate and (for ML/MLR/ADF) the standard error.
  fit_spec <- switch(
    ci_method,
    none            = list(estimator = "ML",  se = "none"),
    feldt           = list(estimator = "ML",  se = "none"),
    fisher          = list(estimator = "ML",  se = "none"),
    bonett          = list(estimator = "ML",  se = "none"),
    hakstian_whalen = list(estimator = "ML",  se = "none"),
    ml              = list(estimator = "ML",  se = "standard"),
    ml_logistic     = list(estimator = "ML",  se = "standard"),
    likelihood      = list(estimator = "ML",  se = "none"),
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
    # Yuan-Bentler sandwich, lavaan's own default pairing for MLR under
    # FIML.
    fit_spec$se <- "robust.huber.white"
  }
  if (missing == "two_stage") {
    # Two-stage estimation fits by ML; the robust intervals take the
    # robust two-stage standard errors (lavaan's "robust.two.stage"),
    # chosen through `missing` below rather than through the estimator.
    fit_spec <- list(estimator = "ML",
                     se = if (fit_spec$se == "none") "none" else "standard")
  }

  fit <- .omega_fit_cfa(data = inp$data, S = inp$S, N = inp$N,
                          h1_diagnostic = inp$h1_diagnostic,
                       equal_loadings = FALSE,
                       estimator = fit_spec$estimator,
                       se = fit_spec$se,
                       missing = .relia_lavaan_missing(
                         missing,
                         robust = ci_method %in% c("mlr", "mlr_logistic")),
                       aux = inp$aux)
  # A solution that yields no coefficient is refused, an improper one is
  # warned about (or refused under improper = "drop" with a bootstrap),
  # before the check for convergence (.omega_point_gate()).
  .omega_point_gate(fit, improper = improper,
                    bootstrap = ci_method %in% bootstrap_methods)
  if (!fit$converged) {
    .dmar_refuse(
      .relia_fit_failure(fit, estimator = fit_spec$estimator,
                         missing = missing), "computation_failed")
  }
  estimate <- fit$omega
  # The degrees of freedom of the Wald interval's t reference: wald_df as
  # given (Inf, the normal reference, by default), or under the keyword
  # max(3, N(1 - FMI) - 1) from omega's own fraction of missing
  # information, with its df_missing() table kept for the line under the
  # result.
  wald_res <- if (wald) {
    .relia_wald_df(wald_df, "reliability_omega", fit = fit, N = inp$N,
                   items = inp$items, aux = inp$aux,
                   J = inp$J)
  } else {
    list(df = Inf, table = NULL)
  }
  df_wald <- wald_res$df

  ci <- switch(
    ci_method,
    none = list(se = NA_real_, lower = NA_real_, upper = NA_real_),
    feldt = .ci_feldt(estimate, inp$N, inp$J, conf_level),
    fisher = .ci_fisher(estimate, inp$N, conf_level),
    bonett = .ci_bonett(estimate, inp$N, inp$J, conf_level),
    hakstian_whalen = .ci_hakstian_whalen(estimate, inp$N, inp$J, conf_level),
    ml = .ci_omega_delta(fit, conf_level, logistic = FALSE,
                         wald_df = df_wald),
    ml_logistic = .ci_omega_delta(fit, conf_level, logistic = TRUE,
                                  wald_df = df_wald),
    likelihood = .ci_reliability_likelihood(
      data = inp$data, S = if (is.null(inp$data)) inp$S else NULL,
      N = inp$N, equal_loadings = FALSE, conf_level = conf_level
    ),
    mlr = .ci_omega_delta(fit, conf_level, logistic = FALSE,
                          wald_df = df_wald),
    mlr_logistic = .ci_omega_delta(fit, conf_level, logistic = TRUE,
                                   wald_df = df_wald),
    adf = .ci_omega_delta(fit, conf_level, logistic = FALSE,
                          wald_df = df_wald),
    adf_logistic = .ci_omega_delta(fit, conf_level, logistic = TRUE,
                                   wald_df = df_wald),
    bootstrap_se = ,
    bootstrap_se_logistic = ,
    percentile = ,
    bca = .bootstrap_ci(
      fun = "reliability_omega",
      data = inp$data,
      # A replication in which an item (or auxiliary variable) takes a
      # single value is refused before the fit, as the data are
      # (.relia_checked_point()).
      point_fn = .relia_checked_point(function(d) {
        # Under two-stage estimation a replication fits the second stage
        # alone to the EM moments (.omega_fit_cfa_ts_estimates()), the
        # same estimates (to the optimizer's tolerance) without the
        # two-stage standard errors and test.
        .omega_point_value(if (missing == "two_stage") {
          .omega_fit_cfa_ts_estimates(d, items = inp$items, aux = inp$aux)
        } else {
          .omega_fit_cfa(
            data = d, equal_loadings = FALSE,
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

  out <- .relia_result(
    estimate = estimate, se = ci$se,
    lower = ci$lower, upper = ci$upper,
    conf_level = conf_level, N = inp$N, J = inp$J,
    coefficient = "omega",
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
  attr(out, "denominator") <- denominator
  attr(out, "wald_df") <- df_wald
  if (!is.null(wald_res$table)) attr(out, "wald_df_table") <- wald_res$table
  attr(out, "wald_calibration_note") <- .wald_eligibility_note(wald_res$table)
  # The local dependence summary is the maximum likelihood one whatever
  # the interval method: the ML and MLR fits share their point estimates
  # and are reused, and the WLS fit behind the ADF intervals is replaced
  # by one ML fit, since its residuals differ and can fall on another
  # pair.
  residual_cor <- .omega_ml_residual_cor(inp, missing, fit = fit,
                                         estimator = fit_spec$estimator)
  if (!is.null(residual_cor)) attr(out, "residual_cor") <- residual_cor
  out
}
