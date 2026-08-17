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
#' Two cautions frame the choice. First, no denominator repairs a
#' misspecified measurement model: the fitted loadings absorb part of
#' whatever structure the single-factor model omits, so the numerator is
#' affected under either setting. Assess the single-factor model (for
#' example with \code{\link{cfa_1}}) before interpreting any
#' \eqn{\omega} variant, and model real multidimensionality directly
#' rather than patching over it. Second, the naming history is worth
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
#' categorical omega: in both, the total variance in the denominator is
#' not taken from the fitted factor model. The reliability vignette
#' develops this framing.
#'
#' Available confidence interval methods (set via \code{ci_method}):
#' \describe{
#'   \item{\code{"ml"}, \code{"ml_logistic"}}{Wald interval using the
#'   maximum likelihood standard error from \pkg{lavaan} (Raykov, 2002).
#'   The \code{_logistic} variant applies Browne's (1982) logit
#'   transformation.}
#'   \item{\code{"mlr"}, \code{"mlr_logistic"}}{Wald interval using the
#'   robust (Satorra & Bentler, 1994) standard error. Default; recommended
#'   among closed-form methods when item distributions deviate from
#'   normality (Kelley & Pornprasertmanit, 2016).}
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
#'   Efron, 1996). \code{"bootstrap_se"} uses the standard deviation
#'   of the bootstrap estimates as a standard error in a normal-theory
#'   interval, built on the logit scale for the \code{_logistic}
#'   variant so the endpoints respect [0, 1]. Replications whose model
#'   refit does not converge are dropped, and the interval is computed
#'   from the replications that return a value. The default
#'   \code{B = 10000} is an accuracy choice: the BCa adjustment pushes
#'   the working quantiles farther into the tails of the bootstrap
#'   distribution than the percentile interval uses, and stabilizing
#'   them takes more replications than the customary 2000; reduce
#'   \code{B} for exploration, not for a reported analysis.
#'   Recommended when assumptions of parametric methods are
#'   questionable. Requires raw data and the \pkg{boot} package;
#'   supply \code{seed} for run-to-run reproducibility.}
#'   \item{\code{"none"}}{Return only the point estimate.}
#' }
#' With \code{denominator = "observed"}, only the bootstrap methods and
#' \code{"none"} are available: the delta method standard errors and the
#' alpha-derived closed forms are derived under the model implied ratio
#' and do not account for sampling of the observed denominator. This
#' pairing is not a limitation in practice, since the bootstrap is the
#' interval Kelley and Pornprasertmanit (2016) recommend for the
#' observed denominator coefficient in any case. Because no analysis in
#' \pkg{DMAR} runs a bootstrap unless the user requests one, the default
#' for robust omega is the point estimate with no interval, accompanied
#' by a message naming the call that produces the recommended interval;
#' request \code{ci_method = "percentile"} or \code{"bca"} to obtain it.
#'
#' \strong{Missing data and auxiliary variables.} By default incomplete
#' rows are listwise-deleted; \code{missing = "fiml"} keeps every case
#' with at least one observed item and fits the single-factor model by
#' full information maximum likelihood, which is consistent and
#' efficient under the missing at random (MAR) assumption. The
#' \code{aux} argument names auxiliary variables: columns of
#' \code{data} that are not part of the composite but are correlated
#' with the items or with the reasons values are missing. They enter as
#' \emph{saturated correlates} (Graham, 2003): correlated freely with
#' each other and with every item's residual, never loading on the
#' factor, so the measurement model is undisturbed while FIML uses
#' their information; a good auxiliary also makes MAR itself more
#' plausible (Collins, Schafer, & Kam, 2001). Supplying \code{aux}
#' implies \code{missing = "fiml"}; combining it with an explicit
#' \code{missing = "listwise"} is an error, and listwise deletion
#' remains the default so no existing result changes. Under
#' \code{missing = "fiml"} the denominators are estimated as follows:
#' with \code{"model_implied"} the fitted FIML model supplies the total
#' variance directly, and with \code{"observed"} the total variance is
#' the sum of the FIML estimate of the item covariance matrix (from a
#' saturated model over the items and any auxiliaries), which is the
#' estimate of \eqn{\mathrm{Var}(Y)} that uses the partially observed
#' rows; it is already on the maximum likelihood divisor, matching the
#' fitted loadings. The available intervals under \code{"fiml"} are
#' \code{"ml"}, \code{"mlr"} (Yuan & Bentler, 2000, robust), their
#' \code{_logistic} variants, and the bootstrap methods, which resample
#' rows including the partially observed ones and refit by FIML; the
#' complete-data closed forms, \code{"adf"}, and \code{"likelihood"}
#' are errors rather than silent fallbacks to listwise deletion.
#' Multiple imputation is a different feature with a different
#' interface and is out of scope here.
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
#'   \code{missing}: listwise deletion by default, or full information
#'   maximum likelihood with \code{missing = "fiml"}. When \code{aux} is
#'   supplied, \code{data} contains both the item columns and the
#'   auxiliary columns, and the columns not named in \code{aux} are the
#'   items.
#' @param S A symmetric covariance matrix among the items. If supplied,
#'   \code{N} must also be supplied; methods that require raw data
#'   (\code{"mlr*"}, \code{"adf*"}, \code{"bootstrap_*"},
#'   \code{"percentile"}, \code{"bca"}) are then unavailable, as are
#'   \code{missing = "fiml"} and \code{aux}, since a covariance matrix
#'   has no incomplete cases.
#' @param N Total sample size; required when \code{S} is supplied.
#' @param ci_method Method for constructing the confidence interval. See
#'   \emph{Details}. When not supplied: for
#'   \code{denominator = "model_implied"} the default is \code{"mlr"}
#'   with raw data (use \code{"ml"} with covariance input); for
#'   \code{denominator = "observed"} (robust omega) no interval is
#'   computed by default, because its interval is bootstrap based and a
#'   bootstrap is never run unless requested. Ask for
#'   \code{"percentile"} or \code{"bca"} to obtain the recommended
#'   interval.
#' @param denominator How the total variance in the denominator of
#'   \eqn{\omega} is estimated: \code{"observed"} (default; robust
#'   omega) uses the variance of the composite estimated directly from
#'   the data; \code{"model_implied"} uses the total variance reproduced
#'   by the fitted single-factor model. See \emph{Details}. With
#'   \code{"observed"}, \code{ci_method} must be a bootstrap method or
#'   \code{"none"}.
#' @param missing How incomplete rows of \code{data} are handled:
#'   \code{"listwise"} (the default; complete-case analysis, the
#'   historical behavior) or \code{"fiml"} (full information maximum
#'   likelihood, using every case with at least one observed item).
#'   See \emph{Details}.
#' @param aux Optional character vector naming auxiliary variable
#'   columns of \code{data}, entered as saturated correlates under full
#'   information maximum likelihood. Supplying \code{aux} implies
#'   \code{missing = "fiml"}. See \emph{Details}.
#' @param conf_level Confidence level for the interval. Defaults to
#'   \code{0.95}.
#' @param B Number of bootstrap replications when a bootstrap method is
#'   selected. Defaults to \code{10000}.
#' @param seed Random number seed used for bootstrap reproducibility.
#'   Defaults to \code{NULL}, which leaves the user's current RNG state intact; supply an integer for reproducibility.
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
#' \code{"lower_limit"} and \code{"upper_limit"} (clamped to [0, 1]),
#' \code{"conf_level"}, \code{"N"} (the cases the analysis used: the
#' complete cases under listwise deletion, every case with at least one
#' observed item under \code{missing = "fiml"}), \code{"N_complete"}
#' (the complete cases; equal to \code{"N"} under listwise deletion and
#' for covariance input), and \code{"J"}. Attributes
#' \code{coefficient} (\code{"omega"}), \code{ci_method},
#' \code{denominator}, \code{missing}, and (when supplied) \code{aux}
#' record the computation; bootstrap calls also record \code{B}.
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
#' # replication. That refitting is why it is not run here; the call is
#' #   reliability_omega(data = items, ci_method = "percentile", B = 10000,
#' #                     seed = 113)
#' # with ci_method = "bca" as the alternative. The percentile interval is
#' # what Kelley and Pornprasertmanit (2016) recommend for hierarchical
#' # omega, and the default B = 10000 is an accuracy choice rather than a
#' # formality: a reported interval deserves the full count (see Details).
#'
#' # Model implied denominator with its closed-form robust ML interval.
#' reliability_omega(data = items, denominator = "model_implied")
#'
#' # Two further routes into the same coefficient are shown rather than
#' # run, since each one fits the single factor model again. The first
#' # works from the summary statistics a paper reports, a covariance
#' # matrix and its sample size:
#' #   reliability_omega(S = cov(items), N = 200,
#' #                     denominator = "model_implied", ci_method = "ml")
#' # The second is full information maximum likelihood with an auxiliary
#' # variable, where missingness on y2 depends on an auxiliary z (missing
#' # at random given z). Supplying aux implies missing = "fiml":
#' #   z <- eta + rnorm(200, sd = 0.5)
#' #   d <- data.frame(items, z = z)
#' #   d$y2[runif(200) < plogis(-1 + 1.5 * as.numeric(scale(z)))] <- NA
#' #   reliability_omega(data = d, aux = "z",
#' #                     denominator = "model_implied")
#'
#' @keywords htest multivariate
#' @family reliability
#'
#' @export

reliability_omega <- function(data = NULL, S = NULL, N = NULL,
                              ci_method = c("mlr", "ml",
                                            "mlr_logistic", "ml_logistic",
                                            "likelihood",
                                            "adf", "adf_logistic",
                                            "feldt", "fisher", "bonett",
                                            "hakstian_whalen",
                                            "bootstrap_se",
                                            "bootstrap_se_logistic",
                                            "percentile", "bca", "none"),
                              denominator = c("observed", "model_implied"),
                              missing = c("listwise", "fiml"),
                              aux = NULL,
                              conf_level = 0.95,
                              B = 10000,
                              seed = NULL) {
  ci_supplied <- !missing(ci_method)
  ci_method <- match.arg(ci_method)
  denominator <- match.arg(denominator)
  missing_supplied <- !missing(missing)
  missing <- match.arg(missing)
  if (!is.null(aux)) {
    if (missing_supplied && missing == "listwise") {
      stop("'aux' cannot be combined with missing = \"listwise\": an ",
           "auxiliary variable only carries information through full ",
           "information maximum likelihood. Drop 'aux' or set ",
           "missing = \"fiml\".", call. = FALSE)
    }
    # Supplying aux implies FIML, so the common case needs one
    # argument, not two.
    missing <- "fiml"
  }
  if (missing == "fiml" && is.null(data)) {
    stop("missing = \"fiml\" (and 'aux') requires raw 'data': a ",
         "covariance matrix has no incomplete cases for full ",
         "information maximum likelihood to use.", call. = FALSE)
  }
  if (!is.numeric(conf_level) || length(conf_level) != 1L ||
      conf_level <= 0 || conf_level >= 1) {
    stop("'conf_level' must be a single number in (0, 1).", call. = FALSE)
  }

  if (missing == "fiml") {
    # Only intervals whose sampling theory can respect the missingness
    # pattern are offered; nothing silently falls back to listwise
    # deletion. The delta method ("ml"/"mlr") standard errors come from
    # the FIML information matrix, and the bootstrap resamples the
    # partially observed rows and refits by FIML. The remaining methods
    # are complete-data constructions.
    if (ci_method %in% c("feldt", "fisher", "bonett", "hakstian_whalen")) {
      stop("ci_method = \"", ci_method, "\" is a complete-data closed ",
           "form: its standard error treats every case as fully ",
           "observed, which overstates the information in partially ",
           "observed rows. With missing = \"fiml\" use \"ml\" (or ",
           "\"mlr\"), whose standard error comes from the FIML ",
           "information matrix, or a bootstrap method.", call. = FALSE)
    }
    if (ci_method %in% c("adf", "adf_logistic")) {
      stop("ci_method = \"", ci_method, "\" is not available with ",
           "missing = \"fiml\": the asymptotic distribution-free ",
           "variance requires complete raw-data fourth moments, and no ",
           "clean ADF counterpart exists under full information ",
           "maximum likelihood. Use \"ml\", \"mlr\", or a bootstrap ",
           "method.", call. = FALSE)
    }
    if (ci_method == "likelihood") {
      stop("ci_method = \"likelihood\" has not been implemented for ",
           "missing = \"fiml\". Use \"ml\", \"mlr\", or a bootstrap ",
           "method.", call. = FALSE)
    }
  }

  inp <- .relia_resolve_inputs(data = data, S = S, N = N, allow_S = TRUE,
                               missing = missing, aux = aux)

  # The single-factor congeneric model behind omega is not identified
  # with two items (four free parameters against three observed
  # moments), so a two-item fit returns whichever value the optimizer
  # wanders to. Refuse it rather than report solver noise.
  if (inp$J < 3L) {
    stop("Coefficient omega needs three or more items: the single-factor ",
         "congeneric model is not identified with two. For a two-item ",
         "scale use reliability_alpha(), whose tau-equivalent ",
         "model is just identified.", call. = FALSE)
  }

  if (denominator == "observed" && !ci_supplied) {
    # No closed-form standard error exists for the observed denominator,
    # and the bootstrap is never run unless the user asks for it. Report
    # the point estimate and say exactly how to request the recommended
    # interval. The message waits until after input validation so a call
    # that is about to error does not first receive advice.
    ci_method <- "none"
    if (is.null(data)) {
      message("Robust omega is reported without a confidence interval: ",
              "its interval is bootstrap based, which requires raw data ",
              "rather than a covariance matrix.")
    } else {
      message("Robust omega is reported without a confidence interval by ",
              "default because its interval is bootstrap based. Request it ",
              "with ci_method = \"percentile\" (or \"bca\"); B = 10000 ",
              "replications is the default when you do.")
    }
  }

  bootstrap_methods <- c("bootstrap_se", "bootstrap_se_logistic",
                         "percentile", "bca")
  if (denominator == "observed" &&
      !ci_method %in% c(bootstrap_methods, "none")) {
    stop("ci_method = \"", ci_method, "\" is available only for ",
         "denominator = \"model_implied\": the closed forms and the ",
         "profile likelihood are defined under the model implied total ",
         "variance. With denominator = \"observed\", use a bootstrap ",
         "method (\"percentile\", \"bca\", \"bootstrap_se\", ",
         "\"bootstrap_se_logistic\") or ci_method = \"none\".",
         call. = FALSE)
  }

  needs_raw <- ci_method %in% c("mlr", "mlr_logistic",
                                "adf", "adf_logistic",
                                bootstrap_methods)
  if (needs_raw && is.null(inp$data)) {
    stop("ci_method = \"", ci_method, "\" requires raw 'data'; ",
         "a covariance matrix is not sufficient.", call. = FALSE)
  }

  if (denominator == "observed") {
    # Under missing = "fiml" both pieces of robust omega come from full
    # information maximum likelihood: the loadings in the numerator from
    # the single-factor fit with missing = "ml" (with any auxiliary
    # variables as saturated correlates), and the total variance in the
    # denominator from the FIML estimate of the item covariance matrix
    # (the saturated model over items plus auxiliaries). The FIML
    # covariance matrix is already on the maximum likelihood (N divisor)
    # metric, so numerator and denominator are internally consistent
    # without the complete-case rescaling.
    omega_observed_fiml <- function(d) {
      mf <- .omega_fit_cfa(data = d, equal_loadings = FALSE,
                           estimator = "MLR", se = "none",
                           missing = "ml", aux = inp$aux)
      if (!mf$converged) return(NA_real_)
      fc <- .fiml_cov(d, items = inp$items, aux = inp$aux, se = "none")
      if (!fc$converged) return(NA_real_)
      (sum(mf$loadings))^2 / sum(fc$S)
    }

    estimate <- if (missing == "fiml") {
      omega_observed_fiml(inp$data)
    } else if (!is.null(inp$data)) {
      .omega_observed_from_data(inp$data)
    } else {
      fit <- .omega_fit_cfa(S = inp$S, N = inp$N,
                            equal_loadings = FALSE,
                            estimator = "ML", se = "none",
                            missing = "listwise")
      if (!fit$converged) {
        stop("The single-factor CFA model did not converge.", call. = FALSE)
      }
      # Rescale the supplied covariance total to the fit's N-divisor
      # (maximum likelihood) metric so the ratio is internally consistent.
      (sum(fit$loadings))^2 / (sum(inp$S) * (inp$N - 1) / inp$N)
    }
    if (is.na(estimate)) {
      stop("The single-factor CFA model did not converge.",
           if (missing == "fiml") {
             paste0(" Full information maximum likelihood fails to ",
                    "converge more often than complete-case estimation, ",
                    "especially with small samples, many auxiliary ",
                    "variables, or sparse missingness patterns.")
           } else {
             ""
           },
           call. = FALSE)
    }

    ci <- if (ci_method == "none") {
      list(se = NA_real_, lower = NA_real_, upper = NA_real_)
    } else {
      .bootstrap_ci(
        data = inp$data,
        point_fn = if (missing == "fiml") {
          omega_observed_fiml
        } else {
          .omega_observed_from_data
        },
        B = B, conf_level = conf_level, kind = ci_method, seed = seed
      )
    }

    out <- .relia_result(
      estimate = estimate, se = ci$se,
      lower = ci$lower, upper = ci$upper,
      conf_level = conf_level, N = inp$N, J = inp$J,
      coefficient = "omega",
      ci_method = ci_method,
      B = if (ci_method %in% bootstrap_methods) B else NA_integer_,
      N_complete = inp$N_complete,
      missing = missing,
      aux = inp$aux
    )
    attr(out, "denominator") <- denominator
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
    list(estimator = "MLR", se = "none")   # for bootstrap families
  )
  if (missing == "fiml" && identical(fit_spec$se, "robust.sem")) {
    # With missing = "ml" the canonical robust standard error is the
    # Yuan-Bentler sandwich, lavaan's own default pairing for MLR under
    # FIML.
    fit_spec$se <- "robust.huber.white"
  }

  fit <- .omega_fit_cfa(data = inp$data, S = inp$S, N = inp$N,
                       equal_loadings = FALSE,
                       estimator = fit_spec$estimator,
                       se = fit_spec$se,
                       missing = if (missing == "fiml") "ml"
                                 else "listwise",
                       aux = inp$aux)
  if (!fit$converged) {
    stop("The single-factor CFA model did not converge with estimator = '",
         fit_spec$estimator, "'.",
         if (missing == "fiml") {
           paste0(" Full information maximum likelihood fails to ",
                  "converge more often than complete-case estimation, ",
                  "especially with small samples, many auxiliary ",
                  "variables, or sparse missingness patterns.")
         } else {
           ""
         },
         call. = FALSE)
  }
  if (isTRUE(fit$improper)) {
    warning("Improper solution (Heywood case): at least one error ",
            "variance estimate is negative, so the model reproduces an ",
            "impossible variance for that item. The fit converged, so ",
            "the coefficient is returned, but interpret it with care ",
            "and consider whether the single-factor model is ",
            "appropriate for these items.", call. = FALSE)
  }
  estimate <- fit$omega

  ci <- switch(
    ci_method,
    none = list(se = NA_real_, lower = NA_real_, upper = NA_real_),
    feldt = .ci_feldt(estimate, inp$N, inp$J, conf_level),
    fisher = .ci_fisher(estimate, inp$N, conf_level),
    bonett = .ci_bonett(estimate, inp$N, inp$J, conf_level),
    hakstian_whalen = .ci_hakstian_whalen(estimate, inp$N, inp$J, conf_level),
    ml = .ci_omega_delta(fit, conf_level, logistic = FALSE),
    ml_logistic = .ci_omega_delta(fit, conf_level, logistic = TRUE),
    likelihood = .ci_reliability_likelihood(
      data = inp$data, S = if (is.null(inp$data)) inp$S else NULL,
      N = inp$N, equal_loadings = FALSE, conf_level = conf_level
    ),
    mlr = .ci_omega_delta(fit, conf_level, logistic = FALSE),
    mlr_logistic = .ci_omega_delta(fit, conf_level, logistic = TRUE),
    adf = .ci_omega_delta(fit, conf_level, logistic = FALSE),
    adf_logistic = .ci_omega_delta(fit, conf_level, logistic = TRUE),
    bootstrap_se = ,
    bootstrap_se_logistic = ,
    percentile = ,
    bca = .bootstrap_ci(
      data = inp$data,
      point_fn = function(d) {
        bf <- .omega_fit_cfa(data = d, equal_loadings = FALSE,
                             estimator = "MLR", se = "none",
                             missing = if (missing == "fiml") "ml"
                                       else "listwise",
                             aux = inp$aux)
        if (!bf$converged) NA_real_ else bf$omega
      },
      B = B, conf_level = conf_level, kind = ci_method, seed = seed
    )
  )

  out <- .relia_result(
    estimate = estimate, se = ci$se,
    lower = ci$lower, upper = ci$upper,
    conf_level = conf_level, N = inp$N, J = inp$J,
    coefficient = "omega",
    ci_method = ci_method,
    B = if (ci_method %in% bootstrap_methods) B else NA_integer_,
    N_complete = inp$N_complete,
    missing = missing,
    aux = inp$aux,
    se_transformed = ci$se_transformed,
    se_transform_scale = ci$se_transform_scale
  )
  attr(out, "denominator") <- denominator
  out
}
