#' Maximum Likelihood Multiple Regression
#'
#' Fits a multiple regression model by maximum likelihood, with full
#' information likelihood handling of missing values by default. The
#' formula interface and S3 methods mirror \code{\link[stats]{lm}} so
#' that calls such as \code{coef()}, \code{vcov()}, \code{confint()},
#' \code{summary()}, \code{fitted()}, \code{residuals()}, and
#' \code{predict()} continue to work. Confidence intervals default to
#' the likelihood ratio (profile) form; Wald and bootstrap variants
#' are also available.
#'
#' \strong{Why a separate function from \code{lm}.} \code{lm} uses
#' ordinary least squares and listwise deletes any row with a missing
#' value on the outcome or on any predictor. Two situations motivate a
#' maximum likelihood alternative.
#'
#' First, when a predictor is missing on some rows, listwise deletion
#' can be biased if the missingness mechanism depends on other
#' observed variables (the missing at random or MAR pattern). Full
#' information maximum likelihood (FIML) jointly models the
#' distribution of \eqn{(Y, X_1, \ldots, X_p)}, the outcome and the
#' \eqn{p} predictors, and yields consistent
#' regression estimates under MAR, while listwise estimates can be
#' biased away from the population values (Enders, 2010; Schafer &
#' Graham, 2002).
#'
#' Second, the joint likelihood estimates the predictor distribution
#' as well, so quantities that depend on the predictor moments (the
#' standardized coefficients, the model implied \eqn{R^2}, and the
#' predictor variances and covariances) draw on every row with an
#' observed predictor, not only the rows that are complete on the
#' outcome. When the missing values are confined to the outcome,
#' however, the unstandardized slopes and their standard errors match
#' listwise deletion up to the maximum likelihood \eqn{N} versus
#' \eqn{N - p - 1} variance divisor. The rows with an observed
#' predictor but a missing outcome inform the marginal distribution of
#' \emph{X}, not the conditional distribution of \emph{Y} given
#' \emph{X} that identifies the slopes, so they leave the slope
#' estimates and their conditional-model standard errors unchanged.
#'
#' Full information estimation differs from listwise deletion in the
#' slope estimates themselves when any predictors are missing on some
#' rows, or when auxiliary variables that correlate with the outcome or
#' with the missingness mechanism are supplied through
#' \code{auxiliary} (see below). The bootstrap gives inference that
#' does not depend on the multivariate normality assumption.
#'
#' \strong{Auxiliary variables.} A variable that is not part of the
#' regression but is correlated with the outcome or with the
#' missingness can be supplied through \code{auxiliary}. Auxiliaries
#' are added to the model as saturated correlates (Graham, 2003): each
#' one is correlated with the residual of the outcome, with every
#' predictor, and with every other auxiliary, but is never entered as
#' a predictor, and the predictors keep their free means, variances,
#' and covariances. The model over the outcome, the predictors, and the
#' auxiliaries is then saturated (its chi square is zero on zero degrees
#' of freedom), so the auxiliaries impose no restriction and the focal
#' regression coefficients keep their meaning: on complete data they
#' equal those of \code{lm()} to working precision, and with missing
#' values they are the regression implied by the full information
#' maximum likelihood estimates of the means and covariances of every
#' variable, auxiliaries included, which \code{\link{cov_fiml}} computes.
#' The full information maximum likelihood uses the auxiliaries' observed
#' values to make the MAR assumption hold conditional on more of the
#' observed data and to recover information that listwise deletion
#' discards. This is the inclusive analysis strategy of Collins,
#' Schafer, and Kam (2001): a variable that predicts the missingness
#' or the incomplete outcome belongs in the analysis even when it is
#' of no substantive interest. Auxiliary variables must be numeric and
#' require \code{fixed_x = FALSE} (the default).
#'
#' \strong{Why likelihood ratio confidence intervals by default.}
#' Wald intervals (point estimate \eqn{\pm} \eqn{z_{1 - \alpha/2}}
#' standard error) are symmetric by construction and assume the
#' sampling distribution of the estimator is approximately normal
#' over the relevant range. Likelihood ratio intervals invert the
#' likelihood ratio test directly: an interval contains every value
#' of the parameter that would not be rejected at level \eqn{\alpha}.
#' Likelihood ratio intervals are invariant under monotone
#' reparameterizations, often have better coverage in small samples,
#' and respect parameter boundaries (Pawitan, 2001). The cost is
#' computational: each parameter requires a sequence of refits with
#' that parameter constrained. The examples below ask for the Wald and
#' bootstrap intervals so the help page stays quick; a reported
#' analysis leaves \code{ci_method} at its default.
#'
#' \strong{A \emph{t} reference distribution for the Wald interval.}
#' The Wald interval and every \emph{p}-value refer the ratio of a
#' coefficient to its standard error to the normal distribution, which
#' treats the standard error as known. That holds only in the limit: in
#' small samples, and more so when values are missing, the normal
#' intervals are too short. In von Hippel's (2016, pp. 432--433)
#' simulations of bivariate normal data with half of one variable
#' missing, the normal FIML intervals for the intercept and slope of a
#' regression covered 5 to 6 points below the nominal 95\% at \emph{n}
#' = 25, and 0 to 3 points below it at \emph{n} = 100. He kept the
#' maximum likelihood estimate and its standard error and replaced the
#' normal quantile with a \emph{t} quantile (his Eq. 29), on degrees of
#' freedom that take the missing information into account. With
#' \code{ci_method = "wald"}, \code{wald_df} does this. A number refers
#' every coefficient's interval and \emph{p}-value to the \emph{t}
#' distribution on that many degrees of freedom. The keyword
#' \code{"missing_information"} applies von Hippel's rule for regression
#' coefficients (his Eqs. 31 and 34) to each coefficient separately,
#' \deqn{\nu_j = \max\{3, N(1 - \gamma_j) - k\},}{nu_j = max(3, N(1 -
#' gamma_j) - k),}
#' the complete data degrees of freedom \eqn{N - k} of a regression with
#' \eqn{k} parameters (the slopes and the intercept) with the sample size
#' replaced by the coefficient's effective sample size, where \eqn{N} is
#' the number of cases the fit used and \eqn{\gamma_j} the coefficient's
#' fraction of missing information. The fraction is estimated from the
#' fit as \code{\link{fraction_missing_information}} estimates it, by
#' the unstructured information of Chen and Savalei (2021), which cannot
#' be negative. The model leaves the means and covariances of the
#' outcome, the predictors, and any auxiliary variables unrestricted (it
#' is saturated), and the three estimators of Chen and Savalei then
#' coincide up to numerical error (Chen & Savalei, 2021, pp. 3, 9), so the
#' fraction agrees with \pkg{lavaan}'s
#' \code{parameterEstimates(fit, fmi = TRUE)} to the precision of
#' \pkg{lavaan}'s own computation. Unlike \pkg{lavaan}'s, its computation
#' does not depend on the units of the variables: an outcome measured in
#' dollars, with a standard deviation of 100,000,000, gets the degrees of
#' freedom it gets in millions of dollars (see
#' \code{\link{fraction_missing_information}}). It is computed at the
#' fit's estimates, which \pkg{lavaan} can leave short of the maximum of
#' the likelihood for variables measured in small units (with every
#' standard deviation near .003 the intercept's degrees of freedom moved
#' by .15); the call then warns, and \code{rescale_data = TRUE}, passed
#' through \code{\dots} to \pkg{lavaan}, is the remedy. A fraction that
#' cannot be computed stops the call with the reason, rather than leave a
#' coefficient with no degrees of freedom. A
#' coefficient estimated mostly from fully observed variables keeps
#' nearly all of the complete data degrees of freedom, and with no
#' missing values every coefficient has \eqn{N - k}. The floor of 3
#' keeps a \emph{t} distribution with finite variance when the estimated
#' fraction is near one. \code{\link{df_missing}} computes the same
#' degrees of freedom from any fractions of missing information, and its
#' help page describes both of von Hippel's rules and the evidence for
#' them; the table it returns for the fit is the \code{wald_df_table}
#' attribute. Two cautions. The \emph{t} interval is an approximation,
#' since the small-sample distribution of a maximum likelihood estimate
#' from incomplete data is not exactly \emph{t} (von Hippel, 2016,
#' p. 431). And the standard error is the maximum likelihood one, whose
#' variance estimate divides by \eqn{N} rather than \eqn{N - k}, so with
#' complete data the interval on \eqn{N - k} degrees of freedom is
#' shorter than the interval \code{lm()} reports by the factor
#' \eqn{\sqrt{(N - k)/N}}{sqrt((N - k)/N)}: von Hippel's construction
#' pairs the maximum likelihood standard error with the \emph{t}
#' quantile, and it is reproduced here as he published it, without a
#' variance correction.
#' The default, \code{Inf}, is the normal reference. The
#' \code{test_statistic} column holds the ratio of the estimate to its
#' standard error, the same number under either reference.
#' \code{summary()} heads the printed coefficient table
#' \code{t value} and \code{Pr(>|t|)} under a \emph{t} reference, as
#' \code{summary.lm()} does, and names the degrees of freedom directly
#' beneath it.
#'
#' \strong{The bootstrap interval.} With \code{ci_method = "boot"} the
#' rows of \code{data} are resampled with replacement \code{B}
#' times (1000 by default) and the model is refit on each resample;
#' with \code{boot_type = "bollen.stine"} the resamples are instead
#' drawn from data transformed to satisfy the fitted model (Bollen &
#' Stine, 1992), a model-based bootstrap. Only the percentile interval
#' is offered: each coefficient's limits are the empirical quantiles
#' of its resampled estimates (Efron & Tibshirani, 1993); there is no
#' bias-corrected and accelerated (BCa) or bootstrap standard error
#' variant. Raising \code{B} tightens the Monte Carlo error of the
#' reported limits. Bootstrap results vary from run to run; supply
#' \code{boot_seed} for reproducibility.
#'
#' \strong{Bootstrap replications.} Each coefficient's limits are quantiles
#' of its \code{B} resampled estimates, so the interval is run only when
#' \code{B} puts at least 10 resamples beyond each limit (400 at the 95
#' percent level). A replication in which some coefficient is undefined
#' (the refit fails or does not converge, or returns a coefficient that is
#' not a finite number) is dropped and counted, never replaced, and the
#' interval is computed from the rest. The result records \code{B}
#' (requested), \code{B_used} (the number the interval rests on),
#' \code{n_dropped}, and \code{B_required} (the fewest replications a
#' request may ask for; that count is held to \code{B}, not to
#' \code{B_used}) as elements of the fitted object, and
#' \code{print()} and \code{summary()} state them in one line, with
#' \code{bootstrap_status} (\code{"complete"}, \code{"dropped"},
#' \code{"exploratory"}, \code{"insufficient_count"}, or
#' \code{"insufficient_share"}). One warning of class
#' \code{dmar_bootstrap_replications_dropped} states the counts when any
#' replication was dropped, with a caution when the dropped replications
#' are at least the 2.5 percent of \code{B} that lies beyond one limit at
#' the 95 percent level, because they could then decide by themselves where
#' that limit falls. When more than half were dropped, the limits of every
#' coefficient are \code{NA},
#' the estimate is returned, and the warning has class
#' \code{dmar_bootstrap_insufficient}. Such a resample is dropped for every
#' coefficient, so all intervals rest on the same resamples.
#' \code{options(dmar.bootstrap_exploratory = TRUE)} lets a smaller
#' \code{B} (2 or more) run for a first look or a demonstration; every such
#' call gives a warning of class \code{dmar_bootstrap_exploratory} and the
#' result is not one to report. See \code{\link{dmar_bootstrap}} for the
#' rule and its reasons.
#'
#' \strong{Model representation.} Internally the model is fit through
#' \pkg{lavaan} as a structural equation model in which \emph{Y} is
#' regressed on the predictors and (when \code{fixed_x = FALSE}) the
#' predictor distribution is also estimated. With
#' \code{fixed_x = FALSE} and complete data, point estimates of the
#' slopes are identical to \code{lm} and the residual variance
#' differs only by the usual \eqn{N} versus \eqn{N - p - 1}
#' divisor.
#'
#' \strong{Caveats.} The function assumes that, conditional on the
#' modeled predictors, the dependent variable is normally distributed
#' with constant variance. Missingness is assumed to be at most MAR;
#' missing not at random patterns require selection or pattern mixture
#' models outside the scope of this function. Factor predictors and
#' interactions are
#' expanded through \code{\link[stats]{model.matrix}} and entered as
#' numeric covariates, so the same caveats about dummy variable
#' encoding that apply to \code{lm} apply here as well.
#'
#' @param formula A two-sided \code{\link[stats]{formula}} of the form
#'   \code{y ~ x1 + x2 + ...}. Factor predictors, interactions
#'   (\code{x1 * x2}), polynomial terms (\code{poly(x, 2)}), and
#'   transformations (\code{I(x^2)}) are supported through the usual
#'   \code{\link[stats]{model.matrix}} expansion. An intercept-only
#'   formula, \code{y ~ 1}, fits the null model (the mean and the
#'   residual variance only); it is the natural restricted model in a
#'   model comparison and pairs with \code{anova()} for a likelihood
#'   ratio test against a fuller model.
#' @param data A \code{data.frame} containing the variables in
#'   \code{formula}. Rows with missing values on \emph{any} of the
#'   modeled variables are retained (under \code{missing = "fiml"})
#'   and contribute to the likelihood through whichever components
#'   are observed.
#' @param missing Character; how missing values are handled. The
#'   default \code{"fiml"} (equivalently \code{"ml"} in \pkg{lavaan})
#'   uses the full information maximum likelihood over all rows that have at
#'   least one observed value on the modeled variables. Other choices
#'   are \code{"listwise"} (drop rows with any missing modeled
#'   variable), \code{"pairwise"} (sample moments computed pairwise;
#'   not recommended with this model), and
#'   \code{"available.cases"}.
#' @param ci_method Character; the method for confidence intervals on
#'   the regression coefficients. \code{"profile"} (default) inverts
#'   the likelihood ratio test on each parameter through a sequence
#'   of constrained refits. \code{"wald"} returns the symmetric
#'   estimate \eqn{\pm} \eqn{z_{1 - \alpha/2}} standard error
#'   interval. \code{"boot"} resamples the rows of \code{data}
#'   \code{B} times, refits the model on each resample, and
#'   reports the percentile interval of the resampled slopes.
#' @param conf_level Desired level of confidence (the complement of
#'   the Type I error rate). Defaults to \code{0.95}.
#' @param B Number of bootstrap resamples when \code{ci_method = "boot"},
#'   a single whole number. Default \code{1000}, because every resample
#'   refits the model. The fewest accepted is 400 at the default
#'   \code{conf_level} (200 at 0.90, and 2000 at 0.99, which the default
#'   does not meet); a smaller \code{B} is refused before anything is
#'   computed (class \code{dmar_count_invalid}). See \emph{Bootstrap
#'   replications} under \emph{Details} and \code{\link{dmar_bootstrap}}.
#' @param boot_type Character; \code{"ordinary"} (default) for the
#'   nonparametric resampling of rows, or \code{"bollen.stine"} for
#'   the Bollen and Stine (1992) model-based bootstrap implemented in
#'   \pkg{lavaan}.
#' @param boot_seed Integer or \code{NULL}; optional seed for the
#'   bootstrap resampling RNG. The default \code{NULL} draws the
#'   resamples from, and advances, the current random number stream,
#'   so successive bootstrap calls draw fresh resamples; supply an
#'   integer to make the bootstrap reproducible. When supplied, the
#'   function seeds the RNG locally and restores the prior state on exit
#'   so the user's global RNG is not polluted.
#' @param estimator Character; the \pkg{lavaan} estimator. One of
#'   \code{"ML"} (default), \code{"MLR"}, \code{"MLM"}, or
#'   \code{"GLS"}. \code{ML} is standard maximum likelihood under
#'   conditional normality of \emph{Y}. \code{MLR} is robust maximum
#'   likelihood with Yuan-Bentler scaled standard errors and a
#'   Yuan-Bentler scaled test statistic, recommended when the
#'   conditional distribution of \emph{Y} departs from normality and
#'   FIML is in use (Yuan & Bentler, 2000). \code{MLM} is Satorra-
#'   Bentler scaled \eqn{\chi^2} statistics under complete data
#'   (Satorra & Bentler, 1994). \code{GLS} is normal theory generalized
#'   least squares. The ordinal-data estimators
#'   (\code{"DWLS"}, \code{"WLS"}, \code{"ULS"} and their robust
#'   variants) are intentionally not exposed; for ordinal outcomes,
#'   fit a different model class.
#' @param se Character or \code{NULL}; the standard error type passed
#'   to \pkg{lavaan}. When \code{NULL} (default), the standard error
#'   type is chosen automatically from \code{estimator}: \code{"ML"}
#'   and \code{"GLS"} use \code{"standard"}, \code{"MLR"} uses
#'   \code{"robust.huber.white"} (Huber-White heteroskedasticity
#'   consistent), and \code{"MLM"} uses \code{"robust.sem"}
#'   (Satorra-Bentler). Override only when the default does not match
#'   the desired analysis. Other values include \code{"robust"},
#'   \code{"first.order"}, and \code{"none"}.
#' @param fixed_x Logical; whether to treat predictors as fixed (not
#'   modeled jointly) or as random (jointly modeled). Defaults to
#'   \code{FALSE}, which is required for the full information
#'   likelihood to use rows with missing predictors. Set to
#'   \code{TRUE} only when the predictors are fully observed and the
#'   user wants the \code{lm}-style conditional model.
#' @param auxiliary Character vector of variable names in \code{data}
#'   to include as auxiliary variables, or \code{NULL} (default) for
#'   none. Auxiliaries are entered as saturated correlates (Graham,
#'   2003): correlated with the outcome residual, every predictor, and
#'   each other, but not as predictors, while the predictors keep their
#'   free covariances, so the model is saturated, the regression
#'   coefficients keep their meaning, and the full information maximum
#'   likelihood draws on the auxiliaries' observed values (the inclusive
#'   analysis strategy; Collins, Schafer, & Kam, 2001). A name in
#'   \code{auxiliary} must be numeric, must be present in \code{data},
#'   must not appear in \code{formula}, and requires
#'   \code{fixed_x = FALSE}. Auxiliaries are refused with
#'   \code{missing = "listwise"}, which would drop every row on which an
#'   auxiliary is missing and so change the sample the coefficients
#'   describe, and with \code{missing = "pairwise"} or
#'   \code{"available.cases"}, which compute the moments from each pair of
#'   variables on its own, so that an auxiliary would leave every estimate
#'   unchanged.
#' @param effect_sizes Logical; whether to compute regression effect
#'   sizes (standardized betas, semi-partial \eqn{R^2}, Cohen's
#'   \eqn{f^2} per predictor, and the overall LR omnibus test).
#'   Defaults to \code{TRUE}. Disabling saves \eqn{p + 1} additional
#'   lavaan refits.
#' @param enforce_es_bounds Logical; whether to clamp
#'   semi-partial \eqn{R^2} and Cohen's \eqn{f^2} estimates to their
#'   theoretical lower bound of zero. Defaults to \code{FALSE}: the
#'   raw maximum likelihood estimates of \eqn{R^2_{\text{reduced}}}
#'   and \eqn{R^2_{\text{full}}} are reported as-is, and the
#'   difference can be slightly negative as a finite-sample artifact
#'   when the two are nearly equal. Setting to \code{TRUE} replaces
#'   any negative value with zero, which yields an estimate that
#'   respects the parameter space but is no longer the maximum
#'   likelihood estimate. When the clamp fires, the affected rows of
#'   the returned \code{effect_sizes} table carry the attribute
#'   \code{"clamped"} for diagnostics.
#' @param \dots Additional arguments forwarded to
#'   \code{\link[lavaan]{lavaan}}.
#' @param wald_df The reference distribution of the Wald interval and of
#'   the \emph{p}-values: \code{Inf} (the default), the normal
#'   distribution; a single positive number, the degrees of freedom of a
#'   \emph{t} distribution for every coefficient; or
#'   \code{"missing_information"}, von Hippel's (2016) degrees of freedom
#'   \eqn{\max\{3, N(1 - \gamma_j) - k\}} for each coefficient, from its
#'   fraction of missing information \eqn{\gamma_j}. A finite value or
#'   the keyword requires \code{ci_method = "wald"}. The keyword is an
#'   error where, with values missing, \pkg{lavaan} fit the complete
#'   cases in place of the FIML fit that \code{missing = "fiml"} asks
#'   for, as it does for \code{se = "robust.sem"} or a test statistic of
#'   the Satorra-Bentler kind passed through \code{\dots}; with no value
#'   missing that fit is the FIML fit, and the keyword gives \eqn{N - k}.
#'   It follows \code{\dots}, so it must be named in full. See
#'   \emph{Details}.
#'
#' @return An object of class \code{"mlmr"}, a list with components
#'   modeled on the structure of an \code{lm} fit:
#'   \describe{
#'     \item{\code{call}}{The matched call.}
#'     \item{\code{formula}}{The model formula.}
#'     \item{\code{terms}}{The terms object.}
#'     \item{\code{model}}{The model frame (with missing values
#'       preserved when \code{missing = "fiml"}).}
#'     \item{\code{coefficients}}{Named numeric vector of regression
#'       coefficients, with \code{(Intercept)} first when an intercept
#'       is in the formula.}
#'     \item{\code{vcov}}{The variance-covariance matrix of the
#'       regression coefficients, returned by \code{vcov()}.}
#'     \item{\code{ci}}{A two-column matrix (\code{lower}, \code{upper})
#'       of confidence limits in the order of \code{coefficients}.}
#'     \item{\code{ci_method}}{Which method was used to compute
#'       \code{ci}.}
#'     \item{\code{B}, \code{B_used}, \code{n_dropped},
#'       \code{B_required}, \code{bootstrap_status}}{Present only with
#'       \code{ci_method = "boot"}: the resamples requested, the number
#'       the intervals rest on, their difference, the fewest resamples a
#'       request may ask for, and the verdict
#'       (\code{"complete"}, \code{"dropped"}, \code{"exploratory"},
#'       \code{"insufficient_count"}, or \code{"insufficient_share"};
#'       under the last two the limits are \code{NA}). See
#'       \emph{Bootstrap replications} under Details.}
#'     \item{\code{conf_level}}{The confidence level used.}
#'     \item{\code{coef_table}}{A \code{data.frame} with columns
#'       \code{term}, \code{estimate}, \code{se}, \code{test_statistic},
#'       \code{p_value}, \code{ci_lower}, \code{ci_upper}, carrying the
#'       confidence level in its \code{conf_level} attribute, so that
#'       \code{\link{results_sentence}} on a row reports it.}
#'     \item{\code{sigma2}}{Residual variance of \emph{Y}, on the
#'       maximum likelihood scale (divisor \emph{N}, not
#'       \eqn{N - p - 1}).}
#'     \item{\code{R2}}{Model implied squared multiple correlation,
#'       \eqn{1 - \hat{\sigma}^2_e / \hat{\sigma}^2_Y}, where both
#'       variances come from the FIML estimated model implied
#'       covariance matrix.}
#'     \item{\code{adj_R2}}{Adjusted \eqn{R^2} using the number of
#'       complete cases (\code{N_complete}, the rows complete on the
#'       outcome and every predictor, which are the rows that identify
#'       the regression) and the number of slopes; the lavaan reported
#'       \code{N} can be larger under FIML because it counts rows that
#'       inform only the predictor distribution.}
#'     \item{\code{logLik}}{The log likelihood at the maximum, with
#'       attributes \code{df} and \code{nobs} for compatibility with
#'       \code{stats::AIC} and \code{stats::BIC}.}
#'     \item{\code{N}}{Sample size used by lavaan (rows with at least
#'       one observed value when \code{missing = "fiml"}; rows with
#'       no missing values when \code{missing = "listwise"}).}
#'     \item{\code{N_complete}}{Number of rows that are complete on
#'       all modeled variables.}
#'     \item{\code{fitted.values}}{Vector of fitted values, length
#'       \code{nrow(data)}, with \code{NA} for rows missing any
#'       predictor.}
#'     \item{\code{residuals}}{Vector of residuals
#'       (\code{y - fitted}), with \code{NA} where \code{y} or any
#'       predictor was missing.}
#'     \item{\code{lavaan_fit}}{The underlying \pkg{lavaan} fit
#'       object, returned for advanced users who want to apply
#'       \pkg{lavaan} accessors directly.}
#'   }
#'   The object carries the attribute \code{wald_df}, the degrees of
#'   freedom of the \emph{t} reference of the Wald interval and the
#'   \emph{p}-values: \code{Inf} for the normal reference, the number
#'   supplied, or under \code{wald_df = "missing_information"} a named
#'   vector with one value per coefficient. In the last case the
#'   attribute \code{wald_df_table} is the \code{\link{df_missing}} table
#'   the values came from, with each coefficient's fraction of missing
#'   information, \eqn{N}, \eqn{k}, and whether the floor applied.
#'
#' @references
#' Bollen, K. A., & Stine, R. A. (1992). Bootstrapping goodness of fit
#' measures in structural equation models. \emph{Sociological Methods
#' & Research, 21}, 205--229. \doi{10.1177/0049124192021002004}
#'
#' Chen, L., & Savalei, V. (2021). Three sample estimates of fraction of
#'   missing information from full information maximum likelihood.
#'   \emph{Frontiers in Psychology, 12}, 667802.
#'   \doi{10.3389/fpsyg.2021.667802}
#'
#' Collins, L. M., Schafer, J. L., & Kam, C.-M. (2001). A comparison
#' of inclusive and restrictive strategies in modern missing data
#' procedures. \emph{Psychological Methods, 6}(4), 330--351.
#'   \doi{10.1037/1082-989X.6.4.330}
#'
#' Efron, B., & Tibshirani, R. J. (1993). \emph{An introduction to the
#' bootstrap}. New York, NY: Chapman & Hall/CRC.
#'
#' Enders, C. K. (2010). \emph{Applied missing data analysis}. New
#' York, NY: Guilford Press.
#'
#' Graham, J. W. (2003). Adding missing-data-relevant variables to
#' FIML-based structural equation models. \emph{Structural Equation
#' Modeling, 10}(1), 80--100. \doi{10.1207/S15328007SEM1001_4}
#'
#' Pawitan, Y. (2001). \emph{In all likelihood: Statistical modelling
#' and inference using likelihood}. Oxford, UK: Oxford University
#' Press.
#'
#' Rosseel, Y. (2012). lavaan: An R package for structural equation
#' modeling. \emph{Journal of Statistical Software, 48}(2), 1--36.
#'   \doi{10.18637/jss.v048.i02}
#'
#' Satorra, A., & Bentler, P. M. (1994). Corrections to test statistics
#' and standard errors in covariance structure analysis. In A. von Eye &
#' C. C. Clogg (Eds.), \emph{Latent variables analysis: Applications for
#' developmental research} (pp. 399--419). Sage.
#'
#' Schafer, J. L., & Graham, J. W. (2002). Missing data: Our view of
#' the state of the art. \emph{Psychological Methods, 7}, 147--177.
#'   \doi{10.1037/1082-989X.7.2.147}
#'
#' von Hippel, P. T. (2016). New confidence intervals and bias comparisons
#'   show that maximum likelihood can beat multiple imputation in small
#'   samples. \emph{Structural Equation Modeling, 23}(3), 422--437.
#'   \doi{10.1080/10705511.2015.1047931}
#'
#' Yuan, K.-H., & Bentler, P. M. (2000). Three likelihood-based methods
#' for mean and covariance structure analysis with nonnormal missing
#' data. \emph{Sociological Methodology, 30}(1), 165--200.
#'   \doi{10.1111/0081-1750.00078}
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso \code{\link[stats]{lm}}, \code{\link[lavaan]{sem}},
#'   \code{\link[lavaan]{lavaan}}, \code{\link{ci_reg_coef}},
#'   \code{\link{ci_rc}}, \code{\link{ci_src}}.
#'
#' @examples
#' # Complete data: the maximum likelihood estimates agree with lm() to
#' # working precision. The fit asks for the Wald interval; the choice
#' # among the intervals is taken up below.
#' fit_mlmr <- mlmr(t6_paragraph_comprehension ~ t5_general_information +
#'                    t9_word_meaning,
#'                  data = holzinger_swineford, ci_method = "wald")
#' fit_lm   <- lm(t6_paragraph_comprehension ~ t5_general_information +
#'                  t9_word_meaning,
#'                data = holzinger_swineford)
#' cbind(mlmr = coef(fit_mlmr), lm = coef(fit_lm))
#'
#' fit_mlmr
#' confint(fit_mlmr)
#'
#' # summary() adds the intervals, the per-predictor semi-partial R^2
#' # and Cohen's f^2, and the omnibus likelihood ratio test of all
#' # slopes equal to zero.
#' summary(fit_mlmr)
#'
#' # The interval menu is profile, Wald, and bootstrap. The default,
#' # ci_method = "profile", inverts the likelihood ratio test one
#' # parameter at a time through a sequence of constrained refits; it
#' # is what a reported interval deserves, and leaving ci_method at its
#' # default asks for it. The fits on this page ask for the Wald or the
#' # bootstrap interval because those refits take longer than a help
#' # page should. The bootstrap resamples rows and takes percentile
#' # limits; it is what to ask for when the normality the likelihood
#' # assumes is doubtful. Every resample refits the model. B = 10 keeps
#' # the example quick and is below the 400 a reported interval needs;
#' # the option lets it run and the warning says so. A reported interval
#' # deserves the default B = 1000. boot_seed fixes the resamples, so
#' # the limits are reproducible rather than moving from run to run, and
#' # effect_sizes = FALSE skips the effect size refits, which
#' # summary() above already showed.
#' old <- options(dmar.bootstrap_exploratory = TRUE)
#' fit_boot <- mlmr(t6_paragraph_comprehension ~ t5_general_information +
#'                    t9_word_meaning, data = holzinger_swineford,
#'                  ci_method = "boot", B = 10, boot_seed = 113,
#'                  effect_sizes = FALSE)
#' options(old)
#' confint(fit_boot)
#'
#' # Missing values on a predictor are where maximum likelihood and
#' # least squares part company. The full information likelihood keeps
#' # every row that carries information; listwise deletion keeps only
#' # the rows that are complete. The Holzinger and Swineford battery
#' # carries real missingness for this: the revised second-form test
#' # t26_flags was administered to only 145 of the 301 students, so a
#' # model using it loses more than half the sample under listwise
#' # deletion while the full information fit keeps all 301 rows. Notice
#' # the two sample sizes, and that the two fits do not agree on the
#' # intercept.
#' fit_fiml <- mlmr(t6_paragraph_comprehension ~ t7_sentence +
#'                    t26_flags, data = holzinger_swineford,
#'                  ci_method = "wald", effect_sizes = FALSE)
#' fit_lwd  <- mlmr(t6_paragraph_comprehension ~ t7_sentence +
#'                    t26_flags, data = holzinger_swineford,
#'                  missing = "listwise",
#'                  ci_method = "wald", effect_sizes = FALSE)
#' rbind(FIML = coef(fit_fiml), listwise = coef(fit_lwd))
#' c(N_fiml = nobs(fit_fiml), N_listwise = nobs(fit_lwd))
#'
#' # An auxiliary variable is not a predictor. The complete speed test
#' # t13_straight_and_curved_capitals enters as a saturated correlate,
#' # correlated with the outcome residual and with the predictors, so
#' # the likelihood can draw on it for the rows where t26_flags is
#' # missing while the coefficients keep their meaning. Continuing from
#' # the model above, the coefficients move a little and the slope
#' # standard errors shrink a little: among the 145 students who took
#' # t26_flags its correlation with t13 is .34, so the auxiliary
#' # restores only part of the information the missing scores would
#' # have carried.
#' fit_aux <- mlmr(t6_paragraph_comprehension ~ t7_sentence +
#'                   t26_flags, data = holzinger_swineford,
#'                 ci_method = "wald",
#'                 auxiliary = "t13_straight_and_curved_capitals",
#'                 effect_sizes = FALSE)
#' cbind(no_aux = coef(fit_fiml), aux = coef(fit_aux))
#' cbind(se_no_aux = sqrt(diag(vcov(fit_fiml))),
#'       se_aux = sqrt(diag(vcov(fit_aux))))
#'
#' # With most of a variable missing, the normal reference understates the
#' # uncertainty in a small sample. Here the outcome t26_flags is observed
#' # for 145 of the 301 students, and wald_df = "missing_information"
#' # refers each coefficient to t on von Hippel's degrees of freedom,
#' # N(1 - FMI) - k, fewer than half of the N - k = 298 that complete
#' # data would give. The summary lists each coefficient's fraction of
#' # missing information and degrees of freedom under the intervals.
#' fit_t <- mlmr(t26_flags ~ t1_visual_perception + t4_lozenges,
#'               data = holzinger_swineford, ci_method = "wald",
#'               wald_df = "missing_information", effect_sizes = FALSE)
#' summary(fit_t)
#'
#' @keywords regression models
#'
#' @export
mlmr <- function(formula,
                 data,
                 missing = c("fiml", "ml", "listwise", "pairwise",
                             "available.cases"),
                 ci_method = c("profile", "wald", "boot"),
                 conf_level = 0.95,
                 B = 1000L,
                 boot_type = c("ordinary", "bollen.stine"),
                 boot_seed = NULL,
                 estimator = c("ML", "MLR", "MLM", "GLS"),
                 se = NULL,
                 fixed_x = FALSE,
                 auxiliary = NULL,
                 effect_sizes = TRUE,
                 enforce_es_bounds = FALSE,
                 ...,
                 wald_df = Inf) {
  .dmar_check_conf_level(conf_level)

  call <- match.call()

  .mlmr_require_lavaan("mlmr")

  if (missing(data) || !is.data.frame(data)) {
    .dmar_refuse_argument("'data' must be a data.frame.", arg = "data")
  }
  if (!inherits(formula, "formula")) {
    .dmar_refuse_argument("'formula' must be a formula.", arg = "formula")
  }
  if (length(formula) != 3L) {
    .dmar_refuse_argument(
      "'formula' must be two-sided (response ~ predictors).", arg = "formula")
  }
  if (!is.numeric(conf_level) || length(conf_level) != 1L ||
      conf_level <= 0 || conf_level >= 1) {
    .dmar_refuse_argument(
      "'conf_level' must be a single number in (0, 1).", arg = "conf_level")
  }
  if (!is.logical(fixed_x) || length(fixed_x) != 1L || is.na(fixed_x)) {
    .dmar_refuse_argument(
      "'fixed_x' must be a single TRUE or FALSE.", arg = "fixed_x")
  }

  missing <- match.arg(missing)
  ci_method <- match.arg(ci_method)
  # The number of replications follows the rule of ?dmar_bootstrap when
  # the bootstrap interval is requested: 10 beyond each limit (400 at the
  # 95 percent level). Through DMAR 1.0.x B was not checked.
  if (ci_method == "boot") {
    B <- .dmar_check_B(B, tail = (1 - conf_level) / 2,
                       conf_level = conf_level)
  }
  boot_type <- match.arg(boot_type)
  estimator <- match.arg(estimator)
  # A t reference distribution belongs with the Wald interval, so that the
  # interval and the p-value refer to the same distribution.
  wald_df <- .wald_df_check(wald_df, "mlmr")
  .wald_df_require_wald(wald_df, ci_method == "wald", paste0(
    "and ci_method = \"", ci_method, "\" is not one; wald_df applies ",
    "with ci_method = \"wald\", whose interval and p-value it refers to ",
    "the t distribution together."))
  if (!is.logical(effect_sizes) || length(effect_sizes) != 1L ||
      is.na(effect_sizes)) {
    .dmar_refuse_argument(
      "'effect_sizes' must be a single TRUE or FALSE.", arg = "effect_sizes")
  }
  if (!is.logical(enforce_es_bounds) ||
      length(enforce_es_bounds) != 1L ||
      is.na(enforce_es_bounds)) {
    .dmar_refuse_argument(
      "'enforce_es_bounds' must be a single TRUE or FALSE.",
      arg = "enforce_es_bounds")
  }
  if (!is.null(boot_seed) &&
      (!is.numeric(boot_seed) || length(boot_seed) != 1L ||
       !is.finite(boot_seed))) {
    .dmar_refuse_argument(
      "'boot_seed' must be NULL or a single finite numeric.", arg = "boot_seed")
  }
  if (is.null(se)) {
    se <- switch(estimator,
                 ML  = "standard",
                 GLS = "standard",
                 MLR = "robust.huber.white",
                 MLM = "robust.sem")
  }
  # lavaan uses "ml" internally; "fiml" is the more common name in the
  # missing data literature, so accept both and pass the lavaan name.
  missing_lavaan <- if (missing == "fiml") "ml" else missing

  design <- mlmr_build_design(formula, data)
  y_name      <- design$y_name
  X_names     <- design$X_names
  X_display   <- design$X_display
  has_intercept <- design$has_intercept
  model_data  <- design$model_data
  model_frame <- design$model_frame
  trms        <- design$terms

  if (length(X_names) == 0L && !has_intercept) {
    .dmar_refuse_argument(c(
      "The right hand side of 'formula' has no predictors and no ",
      "intercept after model.matrix expansion; there is nothing for ",
      "mlmr() to estimate. Use 'y ~ 1' for an intercept-only (null) ",
      "model."), arg = "formula")
  }

  syntax <- mlmr_build_syntax(y_name, X_names, has_intercept)

  # Auxiliary variables (saturated correlates; Graham, 2003). When
  # supplied, they are added to the data the likelihood sees and
  # correlated with the outcome residual, the predictors, and each
  # other, but never entered as predictors, and the predictors keep their
  # covariances (mlmr_auxiliary_syntax()), so the model stays saturated,
  # the focal coefficients keep their meaning, and FIML draws on the
  # auxiliaries' rows.
  aux <- mlmr_prepare_auxiliary(auxiliary, data,
                                model_vars = all.vars(formula),
                                fixed_x = fixed_x, missing = missing)
  fit_data <- model_data
  if (!is.null(aux)) {
    syntax <- paste(syntax,
                    mlmr_auxiliary_syntax(y_name, X_names, aux$aux_names),
                    sep = "\n")
    fit_data <- cbind(model_data, aux$aux_data)
  }

  fit_args <- list(model = syntax,
                   data = fit_data,
                   estimator = estimator,
                   se = se,
                   missing = missing_lavaan,
                   fixed.x = fixed_x,
                   meanstructure = TRUE)
  fit_args <- c(fit_args, list(...))

  # Wrap the initial fit so the "observed variances differ by a factor
  # of 1000" warning from lavaan, which is informational and not
  # actionable when predictors are on naturally different scales, is
  # muffled. Other warnings still
  # surface so users can see genuine issues.
  fit <- try(withCallingHandlers(
    do.call(lavaan::sem, fit_args),
    warning = function(w) {
      if (grepl("factor 1000 times larger", w$message, fixed = TRUE)) {
        invokeRestart("muffleWarning")
      }
    }
  ), silent = TRUE)
  if (inherits(fit, "try-error")) {
    .dmar_refuse(c(
      "lavaan failed to fit the model. Error: ",
      attr(fit, "condition")$message), "computation_failed")
  }
  if (!lavaan::lavInspect(fit, "converged")) {
    .dmar_refuse(c(
      "lavaan reports the model did not converge. Check for ",
      "constant predictors, severe collinearity, or excessive ",
      "missingness."), "not_converged")
  }

  pe <- lavaan::parameterEstimates(fit, ci = FALSE)
  beta_labels <- c(if (has_intercept) "b0",
                   if (length(X_names) > 0L) paste0("b_", seq_along(X_names)))
  display_names <- c(if (has_intercept) "(Intercept)", X_display)

  est <- vapply(beta_labels, function(lab) {
    row <- pe[pe$label == lab, , drop = FALSE]
    if (nrow(row) != 1L) NA_real_ else row$est[[1]]
  }, numeric(1))
  se_w <- vapply(beta_labels, function(lab) {
    row <- pe[pe$label == lab, , drop = FALSE]
    if (nrow(row) != 1L) NA_real_ else row$se[[1]]
  }, numeric(1))
  z_val <- est / se_w
  # The reference distribution of the Wald ratio: the normal one by
  # default, or t on wald_df degrees of freedom, a number for every
  # coefficient or, under "missing_information", von Hippel's (2016)
  # max(3, N(1 - FMI) - k) for each coefficient with k the number of
  # regression parameters (.wald_df_resolve(),
  # R/missing_information_internals.R).
  # A message of the computation of the missing information names the
  # variables as the formula does (.fmi_with_names()).
  wald <- .fmi_with_names(
    .wald_df_resolve(
      wald_df, "mlmr", fit = fit,
      rows = pe[match(beta_labels, pe$label), c("lhs", "op", "rhs")],
      context = list(n_coef = length(beta_labels), terms = display_names,
                     requested_missing = if (missing %in% c("fiml", "ml")) {
                       missing
                     })),
    stats::setNames(c(paste(deparse(formula[[2L]]), collapse = ""),
                      X_display, aux$aux_display),
                    c(y_name, X_names, aux$aux_names)))
  if (!is.null(wald$table)) .df_require_finite(wald$df, "mlmr")
  df_coef <- unname(wald$df)
  p_val <- .wald_p(z_val, df_coef)

  ci_mat <- switch(ci_method,
    profile = mlmr_profile_ci(fit, syntax, fit_args, beta_labels,
                              est, se_w, conf_level),
    wald    = mlmr_wald_ci(est, se_w, conf_level, wald_df = df_coef),
    boot    = mlmr_boot_ci(fit, fit_args, beta_labels, conf_level,
                           B, boot_type, boot_seed)
  )

  boot_counts <- mlmr_boot_counts(attr(ci_mat, "bootstrap"))
  attr(ci_mat, "bootstrap") <- NULL

  coefs <- setNames(est, display_names)

  full_vcov <- lavaan::vcov(fit)
  vcov_idx <- match(beta_labels, colnames(full_vcov))
  vcov_beta <- full_vcov[vcov_idx, vcov_idx, drop = FALSE]
  dimnames(vcov_beta) <- list(display_names, display_names)

  rownames(ci_mat) <- display_names

  implied_cov <- lavaan::lavInspect(fit, "cov.ov")
  if (y_name %in% rownames(implied_cov)) {
    sigma2_y <- implied_cov[y_name, y_name]
  } else {
    sigma2_y <- NA_real_
  }
  sigma2_e <- pe$est[pe$label == "sigma2_e"][1]
  R2 <- 1 - sigma2_e / sigma2_y
  N_used <- as.integer(lavaan::lavInspect(fit, "ntotal"))
  p <- length(X_names)
  # Adjusted R^2 uses the number of rows complete on Y AND all
  # predictors (the rows that identify the regression of Y on X),
  # not the total FIML cases (which can include rows that contribute
  # only to the X distribution).
  N_for_R2 <- sum(stats::complete.cases(model_data))
  if (N_for_R2 > p + 1L) {
    adj_R2 <- 1 - (1 - R2) * (N_for_R2 - 1L) / (N_for_R2 - p - 1L)
  } else {
    adj_R2 <- NA_real_
  }
  f2_model <- if (R2 < 1) R2 / (1 - R2) else Inf

  ll_val <- as.numeric(lavaan::logLik(fit))
  npar <- length(lavaan::coef(fit))
  attr(ll_val, "df") <- npar
  attr(ll_val, "nobs") <- N_used
  class(ll_val) <- "logLik"

  # Effect sizes: standardized betas (no extra refit), then
  # per-predictor semi-partial R-squared and Cohen's f-squared via
  # constrained refits (one per slope), plus the overall LR omnibus
  # test (one more refit with all slopes fixed to zero). The full
  # block costs p + 1 extra lavaan fits.
  std_estimate <- rep(NA_real_, length(display_names))
  effect_size_table <- NULL
  omnibus_test <- NULL
  if (effect_sizes && length(X_names) > 0L) {
    SD_Y <- sqrt(sigma2_y)
    SD_X <- sqrt(diag(implied_cov)[X_names])
    slope_idx <- if (has_intercept) seq.int(2L, length(display_names))
                 else seq_along(display_names)
    std_estimate[slope_idx] <- est[slope_idx] * SD_X / SD_Y

    es <- mlmr_compute_effect_sizes(
      fit_args = fit_args,
      syntax = syntax,
      slope_labels = paste0("b_", seq_along(X_names)),
      sigma2_y = sigma2_y,
      full_R2 = R2,
      ll_full = as.numeric(ll_val)
    )
    sr2_vec <- es$sr2
    f2_vec  <- es$f2
    clamped <- logical(length(sr2_vec))
    if (enforce_es_bounds) {
      neg_sr2 <- !is.na(sr2_vec) & sr2_vec < 0
      neg_f2  <- !is.na(f2_vec)  & f2_vec  < 0
      clamped <- neg_sr2 | neg_f2
      sr2_vec[neg_sr2] <- 0
      f2_vec[neg_f2]   <- 0
      if (!is.na(R2) && R2 < 0) {
        R2 <- 0
        f2_model <- 0
      }
    }
    effect_size_table <- data.frame(
      term  = X_display,
      sr2   = sr2_vec,
      f2    = f2_vec,
      stringsAsFactors = FALSE,
      row.names = NULL
    )
    attr(effect_size_table, "clamped") <- clamped
    attr(effect_size_table, "enforce_es_bounds") <- enforce_es_bounds
    omnibus_test <- es$omnibus
  }

  coef_table <- data.frame(
    term         = display_names,
    estimate     = unname(est),
    std_estimate = unname(std_estimate),
    se           = unname(se_w),
    test_statistic      = unname(z_val),
    p_value      = unname(p_val),
    ci_lower     = unname(ci_mat[, 1]),
    ci_upper     = unname(ci_mat[, 2]),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  # The table records the confidence level of its intervals, which a row
  # keeps, so that results_sentence() on a coefficient reports it.
  attr(coef_table, "conf_level") <- conf_level

  N_complete <- sum(stats::complete.cases(model_data))

  intercept_val <- if (has_intercept) est[["b0"]] else 0
  if (length(X_names) > 0L) {
    X_for_fit <- model_data[, X_names, drop = FALSE]
    any_X_missing <- !stats::complete.cases(X_for_fit)
    fitted_vals <- rep(NA_real_, nrow(model_data))
    X_mat <- as.matrix(X_for_fit)
    beta_slopes <- est[!names(est) %in% "b0"]
    fitted_vals[!any_X_missing] <- intercept_val +
      X_mat[!any_X_missing, , drop = FALSE] %*% beta_slopes
  } else {
    # Intercept-only (null) model: every case's fitted value is the
    # estimated mean, and the residual is the deviation from that mean.
    fitted_vals <- rep(intercept_val, nrow(model_data))
  }
  resids <- model_data[[y_name]] - fitted_vals

  xlevels <- stats::.getXlevels(trms, model_frame)

  structure(c(list(
    call          = call,
    formula       = formula,
    terms         = trms,
    model         = model_frame,
    xlevels       = xlevels,
    coefficients  = coefs,
    vcov          = vcov_beta,
    ci            = ci_mat,
    ci_method     = ci_method,
    conf_level    = conf_level,
    coef_table    = coef_table,
    sigma2        = unname(sigma2_e),
    R2            = unname(R2),
    adj_R2        = unname(adj_R2),
    f2            = unname(f2_model),
    effect_sizes  = effect_size_table,
    omnibus_test  = omnibus_test,
    logLik        = ll_val,
    N             = N_used,
    N_complete    = N_complete,
    fitted.values = fitted_vals,
    residuals     = resids,
    missing       = missing,
    estimator     = estimator,
    se_type       = se,
    fixed_x       = fixed_x,
    auxiliary     = if (!is.null(aux)) aux$aux_display else NULL,
    y_name        = y_name,
    X_names       = X_names,
    X_display     = X_display,
    has_intercept = has_intercept,
    syntax        = syntax,
    model_data    = model_data,
    fit_args      = fit_args,
    lavaan_fit    = fit
  ), boot_counts), class = "mlmr",
  wald_df = if (is.null(wald$table)) wald_df else wald$df,
  wald_df_table = wald$table)
}


#' @export
print.mlmr <- function(x, digits = max(3L, getOption("digits") - 3L),
                       ...) {
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
      "\n\n", sep = "")
  cat("Coefficients:\n")
  print.default(format(x$coefficients, digits = digits),
                print.gap = 2L, quote = FALSE)
  cat("\n")
  mlmr_boot_line(x)
  invisible(x)
}


#' @export
summary.mlmr <- function(object, ...) {
  ans <- list(
    call         = object$call,
    formula      = object$formula,
    coef_table   = object$coef_table,
    effect_sizes = object$effect_sizes,
    omnibus_test = object$omnibus_test,
    f2_model     = object$f2,
    sigma2       = object$sigma2,
    sigma        = sqrt(object$sigma2),
    R2           = object$R2,
    adj_R2       = object$adj_R2,
    N            = object$N,
    N_complete   = object$N_complete,
    df_residual  = object$N_complete - length(object$coefficients),
    logLik       = object$logLik,
    conf_level   = object$conf_level,
    ci_method    = object$ci_method,
    missing      = object$missing,
    estimator    = object$estimator,
    se_type      = object$se_type,
    fixed_x      = object$fixed_x,
    auxiliary    = object$auxiliary,
    wald_df      = attr(object, "wald_df"),
    wald_df_table = attr(object, "wald_df_table")
  )
  ans <- c(ans, object[intersect(c("B", "B_used", "n_dropped", "B_required",
                                   "bootstrap_status"), names(object))])
  class(ans) <- "summary.mlmr"
  ans
}


#' @export
print.summary.mlmr <- function(x, digits = max(3L, getOption("digits") - 3L),
                               signif.stars = getOption("show.signif.stars"),
                               ...) {
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
      "\n\n", sep = "")
  cat("Missing data: ", x$missing, " | Estimator: ", x$estimator,
      " | SE: ", x$se_type, "\n", sep = "")
  if (!is.null(x$auxiliary)) {
    cat("Auxiliary variable(s), saturated correlates: ",
        paste(x$auxiliary, collapse = ", "), "\n", sep = "")
  }
  cat("Sample size (used by lavaan): ", x$N,
      "   Complete cases: ", x$N_complete, "\n\n", sep = "")

  coef_mat <- as.matrix(x$coef_table[, c("estimate", "se", "test_statistic",
                                         "p_value")])
  rownames(coef_mat) <- x$coef_table$term
  # The headers name the reference distribution of the p-values, as
  # summary.lm() and summary.glm() do: t when wald_df is finite (the
  # stored column is test_statistic under either reference), with the line
  # naming the degrees of freedom directly beneath the table.
  colnames(coef_mat) <- .wald_df_coef_headers(x$wald_df)
  cat("Coefficients:\n")
  stats::printCoefmat(coef_mat, digits = digits,
                      signif.stars = signif.stars,
                      has.Pvalue = TRUE, ...)
  wd_lines <- .wald_df_summary_lines(x$wald_df, x$wald_df_table)
  if (length(wd_lines)) cat(paste(wd_lines, collapse = "\n"), "\n", sep = "")

  ci_mat <- as.matrix(x$coef_table[, c("ci_lower", "ci_upper")])
  rownames(ci_mat) <- x$coef_table$term
  # The level and the column labels are written whatever
  # getOption("digits") is (.format_conf_prop(), .ci_quantile_labels()):
  # with format() under options(digits = 2) a .975 interval read
  # "level = 0.97" over columns "1.3 %" and "99 %".
  colnames(ci_mat) <- .ci_quantile_labels(x$conf_level)
  cat("\nConfidence intervals (method: ", x$ci_method, ", level = ",
      .format_conf_prop(x$conf_level), "):\n", sep = "")
  print(ci_mat, digits = digits)
  mlmr_boot_line(x)

  if (!is.null(x$effect_sizes)) {
    es_mat <- as.matrix(x$effect_sizes[, c("sr2", "f2")])
    rownames(es_mat) <- x$effect_sizes$term
    colnames(es_mat) <- c("sr^2", "Cohen's f^2")
    cat("\nPer-predictor effect sizes (semi-partial R^2 and f^2):\n")
    print(es_mat, digits = digits)
  }

  cat("\nResidual std. error (ML): ", format(x$sigma, digits = digits),
      " on ", x$df_residual, " residual degrees of freedom\n", sep = "")
  cat("Model implied R-squared: ", format(x$R2, digits = digits),
      ",   Adjusted R-squared: ",
      format(x$adj_R2, digits = digits),
      ",   Cohen's f^2: ", format(x$f2_model, digits = digits),
      "\n", sep = "")
  if (!is.null(x$omnibus_test)) {
    cat("Omnibus likelihood ratio test (all slopes = 0): ",
        "chi square = ", format(x$omnibus_test$statistic, digits = digits),
        " on ", x$omnibus_test$df, " df, p = ",
        format.pval(x$omnibus_test$p_value, digits = digits), "\n", sep = "")
  }
  cat("Log likelihood: ", format(as.numeric(x$logLik), digits = digits),
      "   AIC: ", format(stats::AIC(x$logLik), digits = digits),
      "   BIC: ", format(stats::BIC(x$logLik), digits = digits), "\n",
      sep = "")
  invisible(x)
}


#' @export
coef.mlmr <- function(object, ...) object$coefficients


#' @export
vcov.mlmr <- function(object, ...) object$vcov


#' @export
confint.mlmr <- function(object, parm, level = 0.95, ...) {
  # Only warn when the user explicitly supplied a level that differs
  # from what mlmr was fit with; the bare call confint(fit) should
  # never warn, even if fit used a non-default conf_level.
  if (!missing(level) &&
      !isTRUE(all.equal(level, object$conf_level))) {
    .dmar_warn(c(
      "confint() ignores 'level' for mlmr fits; mlmr returns ",
      "the interval at the conf_level used at fit time (",
      object$conf_level,
      "). Refit with the desired conf_level."), "argument_ignored")
  }
  ci <- object$ci
  if (!missing(parm)) {
    ci <- ci[parm, , drop = FALSE]
  }
  # Labeled as stats::confint.default() labels its columns, whatever
  # getOption("digits") is (.ci_quantile_labels()).
  colnames(ci) <- .ci_quantile_labels(object$conf_level)
  ci
}


#' @export
nobs.mlmr <- function(object, ...) object$N


#' @export
logLik.mlmr <- function(object, ...) object$logLik


#' @export
fitted.mlmr <- function(object, ...) object$fitted.values


#' @export
residuals.mlmr <- function(object, ...) object$residuals


#' @export
formula.mlmr <- function(x, ...) x$formula


#' @export
model.matrix.mlmr <- function(object, ...) {
  stats::model.matrix(object$terms,
                      stats::model.frame(object$terms, object$model,
                                         na.action = stats::na.pass))
}


#' @export
predict.mlmr <- function(object, newdata, ...) {
  if (missing(newdata) || is.null(newdata)) {
    return(object$fitted.values)
  }
  trms <- stats::delete.response(object$terms)
  mf <- stats::model.frame(trms, newdata, na.action = stats::na.pass,
                           xlev = object$xlevels)
  X <- stats::model.matrix(trms, mf)
  drop(X %*% object$coefficients[colnames(X)])
}


#' @export
as.data.frame.mlmr <- function(x, ...) {
  x$coef_table
}


#' An Mlmr Fit
#'
#' Returns a one-row-per-coefficient \code{data.frame} in the column
#' convention used by the \pkg{broom} ecosystem (\code{term},
#' \code{estimate}, \code{se}, \code{statistic},
#' \code{p_value}, \code{wald_df}, and optionally \code{ci_lower},
#' \code{ci_upper}). \code{wald_df} is the degrees of freedom of the
#' \emph{t} reference of the row's statistic and \emph{p}-value (the
#' \code{wald_df} argument of \code{\link{mlmr}}), \code{Inf} for the
#' normal reference, so that rows on the two references can be told apart
#' after tables are stacked.
#' Use \code{as.data.frame()} on the fit for the DMAR-style table
#' (snake_case columns) stored at \code{fit$coef_table}.
#'
#' @param x An object of class \code{"mlmr"}.
#' @param conf.int Logical; if \code{TRUE}, append \code{ci_lower} and
#'   \code{ci_upper} columns using the confidence intervals already
#'   computed at fit time, and a \code{conf_level} column holding the
#'   level they were computed at. Defaults to \code{FALSE}.
#' @param conf_level Ignored; the confidence interval comes from the
#'   fit object at \code{x$conf_level}. Present for compatibility with
#'   the broom generic.
#' @param standardized Logical; if \code{TRUE} and effect sizes were
#'   computed at fit time, append a \code{std_estimate} column.
#'   Defaults to \code{FALSE}.
#' @param \dots Unused.
#'
#' @return A \code{data.frame}.
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @examples
#' fit <- mlmr(t6_paragraph_comprehension ~ t5_general_information +
#'               t9_word_meaning,
#'             data = holzinger_swineford, ci_method = "wald")
#' generics::tidy(fit)
#' generics::tidy(fit, conf.int = TRUE)
#' generics::tidy(fit, conf.int = TRUE, standardized = TRUE)
#'
#' @importFrom generics tidy
#' @exportS3Method generics::tidy
tidy.mlmr <- function(x, conf.int = FALSE, conf_level = NULL,
                      standardized = FALSE, ...) {
  tab <- x$coef_table
  out <- data.frame(
    term      = tab$term,
    estimate  = tab$estimate,
    se = tab$se,
    statistic = tab$test_statistic,
    p_value   = tab$p_value,
    # The reference of the statistic and p-value of each row, which a stack
    # of tables would otherwise lose (.wald_df_column()).
    wald_df   = .wald_df_column(attr(x, "wald_df", exact = TRUE), nrow(tab),
                                tab$term),
    stringsAsFactors = FALSE
  )
  if (isTRUE(conf.int)) {
    out$ci_lower  <- tab$ci_lower
    out$ci_upper <- tab$ci_upper
    # The level the intervals were computed at, so that a row of the tidy
    # table says what level its limits belong to.
    out$conf_level <- rep(x$conf_level %||% NA_real_, nrow(out))
  }
  if (isTRUE(standardized) && "std_estimate" %in% colnames(tab)) {
    out$std_estimate <- tab$std_estimate
  }
  out
}


#' Glance at an Mlmr Fit
#'
#' Returns a one-row \code{data.frame} of model-level summaries in the
#' column convention used by the \pkg{broom} ecosystem
#' (\code{R2}, \code{adj_R2}, \code{sigma},
#' \code{statistic}, \code{p_value}, \code{df}, \code{logLik},
#' \code{AIC}, \code{BIC}, \code{deviance}, \code{df_residual},
#' \code{nobs}). The \code{statistic} and \code{p_value} columns
#' report the omnibus likelihood ratio test of all slopes equal to
#' zero (the FIML analog of the \code{lm} omnibus \emph{F}-test).
#'
#' @param x An object of class \code{"mlmr"}.
#' @param \dots Unused.
#'
#' @return A one-row \code{data.frame}.
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @examples
#' fit <- mlmr(t6_paragraph_comprehension ~ t5_general_information +
#'               t9_word_meaning,
#'             data = holzinger_swineford, ci_method = "wald")
#' generics::glance(fit)
#'
#' @importFrom generics glance
#' @exportS3Method generics::glance
glance.mlmr <- function(x, ...) {
  ll <- as.numeric(x$logLik)
  npar <- attr(x$logLik, "df")
  omni <- x$omnibus_test
  data.frame(
    R2     = x$R2,
    adj_R2 = x$adj_R2,
    sigma         = sqrt(x$sigma2),
    f2            = x$f2,
    statistic     = if (!is.null(omni)) omni$statistic else NA_real_,
    p_value       = if (!is.null(omni)) omni$p_value else NA_real_,
    df            = npar,
    logLik        = ll,
    AIC           = stats::AIC(x$logLik),
    BIC           = stats::BIC(x$logLik),
    deviance      = -2 * ll,
    df_residual   = x$N - length(x$coefficients),
    nobs          = x$N,
    stringsAsFactors = FALSE
  )
}


#' Likelihood Ratio Test for Nested Mlmr Fits
#'
#' Compares two or more nested \code{\link{mlmr}} fits with the
#' likelihood ratio test, delegating the chi square computation to
#' \code{\link[lavaan]{lavTestLRT}}. The models must be nested
#' (every parameter in the more restricted model is also in the more
#' general one) and must be fit to the same data with the same
#' missing data handling.
#'
#' @details
#' "The same data" means the same observations, with the variables the
#' fits share holding the same values. It does not mean the same number
#' of complete cases. Under \code{missing = "fiml"} a legitimate nested
#' comparison routinely has different complete-case counts: adding a
#' predictor that is itself incompletely observed lowers the number of
#' rows complete on every modeled variable, yet the observations, and
#' the information each of them contributes to the likelihood, are
#' unchanged. Comparing \code{y ~ x1} with \code{y ~ x1 + x2} when
#' \code{x2} has missing values is exactly the comparison full
#' information maximum likelihood exists to support, and it is
#' accepted. Fits run on genuinely different data are refused.
#'
#' Each smaller model is refit as a constrained version of the largest
#' one on the largest fit's data, with the slopes of its absent
#' predictors fixed at zero. That keeps the comparison on a single
#' joint observed variable set, which is what makes the chi square
#' difference interpretable when the predictors are modeled as random
#' (\code{fixed_x = FALSE}, the \code{\link{mlmr}} default).
#'
#' @param object An \code{mlmr} fit.
#' @param \dots Additional \code{mlmr} fits, nested with respect to
#'   \code{object}.
#'
#' @return A \code{data.frame} of class \code{anova} reporting the
#'   degrees of freedom, AIC, BIC, log-likelihood, chi square test
#'   statistic, and \emph{p}-value for each consecutive pairwise
#'   comparison.
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @examples
#' # Both fits ask for the Wald interval and skip the effect size block,
#' # since the likelihood ratio test needs neither and each costs refits.
#' fit1 <- mlmr(t6_paragraph_comprehension ~ t5_general_information,
#'              data = holzinger_swineford, ci_method = "wald",
#'              effect_sizes = FALSE)
#' fit2 <- mlmr(t6_paragraph_comprehension ~ t5_general_information +
#'                t9_word_meaning,
#'              data = holzinger_swineford, ci_method = "wald",
#'              effect_sizes = FALSE)
#' anova(fit1, fit2)
#'
#' @export
anova.mlmr <- function(object, ...) {
  fits <- c(list(object), list(...))
  if (length(fits) < 2L) {
    .dmar_refuse_argument(
      "anova.mlmr requires at least two mlmr fits to compare.", arg = NULL)
  }
  if (!all(vapply(fits, inherits, logical(1L), "mlmr"))) {
    .dmar_refuse_argument(
      "All arguments to anova.mlmr must be of class 'mlmr'.", arg = NULL)
  }

  # Cross-fit sanity checks. Likelihoods are only comparable when
  # the fits are over the same data with the same missing data
  # treatment; otherwise the chi square statistic is meaningless.
  # The same-data requirement is enforced further down, on the
  # analysis data itself, once the fits have been ordered and the
  # largest one identified.
  miss <- vapply(fits, function(f) f$missing, character(1L))
  if (length(unique(miss)) != 1L) {
    .dmar_refuse_argument(c(
      "anova.mlmr requires fits to use the same 'missing' ",
      "setting; got (", paste(miss, collapse = ", "),
      "). Likelihoods under different missing-data treatments ",
      "are not comparable."), arg = "missing")
  }
  est <- vapply(fits, function(f) f$estimator, character(1L))
  if (length(unique(est)) != 1L) {
    .dmar_refuse_argument(c(
      "anova.mlmr requires fits to use the same 'estimator' ",
      "(got ", paste(est, collapse = ", "), ")."), arg = "estimator")
  }
  ynm <- vapply(fits, function(f) deparse(f$formula[[2L]]),
                character(1L))
  if (length(unique(ynm)) != 1L) {
    .dmar_refuse_argument(c(
      "anova.mlmr requires fits to share the same response: ",
      "got (", paste(ynm, collapse = ", "), ")."), arg = NULL)
  }

  # Order fits from fewest predictors to most. The largest serves as
  # the common variable set; smaller fits are re-expressed as
  # constrained versions of the largest so the LR test compares
  # likelihoods over a single joint observed variable set, which is
  # required for lavTestLRT to be valid when fixed_x = FALSE.
  K_each <- vapply(fits, function(f) length(f$X_display), integer(1L))
  fits <- fits[order(K_each)]
  largest <- fits[[length(fits)]]

  for (f in fits) {
    if (!all(f$X_display %in% largest$X_display)) {
      .dmar_refuse_argument(c(
        "anova.mlmr requires nested fits: the predictors in every ",
        "smaller fit must be a subset of the largest fit's ",
        "predictors. Offending predictor(s): ",
        paste(setdiff(f$X_display, largest$X_display),
              collapse = ", ")), arg = NULL)
    }
  }

  # Same-data check. The likelihood ratio test below refits every smaller
  # model as a constrained version of the largest fit on the largest
  # fit's data, so it is valid only when each fit saw the same
  # observations. What has to match is the analysis data itself, not the
  # number of rows that happen to be complete on every modeled variable.
  # Those are different requirements, and each direction matters:
  #
  #   1. Equal complete-case counts do not imply the same data. Two fits
  #      of the same row count can be run on entirely different data
  #      frames, and their likelihoods are not comparable.
  #   2. The same data does not imply equal complete-case counts. Under
  #      FIML, adding a predictor that is itself incompletely observed
  #      lowers the number of rows complete on all modeled variables
  #      while leaving the observations, and the likelihood they
  #      contribute, untouched. Comparing y ~ x1 with y ~ x1 + x2 when
  #      x2 has missing values is precisely the comparison FIML exists
  #      to support, so refusing it on a complete-case count would be
  #      wrong.
  #
  # So compare the analysis data directly. Each fit stores it in
  # $model_data with internal, positional column names (y, x_1, ...) and
  # with missing values preserved, one row per row of 'data'; relabeling
  # by the display names lets fits be compared on the variables they
  # share regardless of predictor order.
  relabel <- function(f) {
    md <- f$model_data
    names(md) <- c(paste(deparse(f$formula[[2L]]), collapse = ""),
                   f$X_display)
    md
  }
  largest_data <- relabel(largest)
  for (f in fits[-length(fits)]) {
    fdat <- relabel(f)
    if (nrow(fdat) != nrow(largest_data)) {
      .dmar_refuse_argument(c(
        "anova.mlmr requires all fits to be on the same data. The fit ",
        "'", paste(deparse(f$formula), collapse = " "), "' was run on ",
        nrow(fdat), " observations and the largest fit '",
        paste(deparse(largest$formula), collapse = " "), "' on ",
        nrow(largest_data),
        ". Their likelihoods are not comparable; refit both models on ",
        "the same 'data' before calling anova()."), arg = "data")
    }
    shared <- intersect(names(fdat), names(largest_data))
    mismatched <- shared[!vapply(shared, function(nm) {
      identical(fdat[[nm]], largest_data[[nm]])
    }, logical(1L))]
    if (length(mismatched) > 0L) {
      .dmar_refuse_argument(c(
        "anova.mlmr requires all fits to be on the same data. The fit ",
        "'", paste(deparse(f$formula), collapse = " "), "' and the ",
        "largest fit '", paste(deparse(largest$formula), collapse = " "),
        "' hold different values for: ",
        paste(mismatched, collapse = ", "),
        ". Their likelihoods are not comparable; refit both models on ",
        "the same 'data' before calling anova()."), arg = "data")
    }
  }

  aligned_fits <- lapply(fits, function(f) {
    if (identical(f$X_display, largest$X_display)) {
      return(f$lavaan_fit)
    }
    missing_in_f <- setdiff(largest$X_display, f$X_display)
    slope_positions <- which(largest$X_display %in% missing_in_f)
    extra_constraints <- paste0("b_", slope_positions, " == 0")
    constrained_syntax <- paste(largest$syntax,
                                paste(extra_constraints, collapse = "\n"),
                                sep = "\n")
    new_args <- largest$fit_args
    new_args$model <- constrained_syntax
    refit <- try(suppressWarnings(do.call(lavaan::sem, new_args)),
                 silent = TRUE)
    if (inherits(refit, "try-error") ||
        !lavaan::lavInspect(refit, "converged")) {
      .dmar_refuse(c(
        "Could not refit the smaller model on the common variable ",
        "set with constraints (formula: ",
        paste(deparse(f$formula), collapse = " "),
        "). Try fitting both models on the same data first."),
        "computation_failed")
    }
    refit
  })

  # model.names must be supplied explicitly: without it lavTestLRT
  # derives display names by deparsing its call arguments, and under
  # do.call those arguments are the inlined lavaan objects, whose
  # deparse runs to hundreds of lines. On R and lavaan builds where
  # that deparse is not collapsed, assigning those lines as the names
  # of the model list stops the call ("'names' attribute [N] must be
  # the same length as the vector"); the platform dependence made the
  # failure look intermittent. Short explicit names sidestep the
  # deparse entirely.
  out <- do.call(lavaan::lavTestLRT,
                 c(aligned_fits,
                   list(model.names = paste0("fit", seq_along(aligned_fits)))))
  # lavTestLRT sorts its rows by model degrees of freedom, which need not match
  # the order the fits were supplied in, and it writes deparsed model names
  # into the row names. Map each
  # output row back to its originating fit by model df so the compact "Model i"
  # labels and the heading formulas describe the statistics actually on that
  # row. The row order itself is preserved so the adjacent chi square difference
  # tests remain valid.
  model_df <- vapply(aligned_fits, function(ff) {
    as.numeric(lavaan::fitMeasures(ff, "df"))
  }, numeric(1L))
  row_fit <- match(out[["Df"]], model_df)   # row_fit[j] = index of the fit in row j
  rownames(out) <- paste("Model", row_fit)
  attr(out, "heading") <- c(
    "Likelihood ratio test for nested mlmr fits",
    vapply(seq_len(nrow(out)), function(j) {
      i <- row_fit[j]
      sprintf("Model %d: %s", i,
              paste(deparse(fits[[i]]$formula), collapse = " "))
    }, character(1L))
  )
  out
}
