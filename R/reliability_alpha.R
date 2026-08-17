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
#' ordered-categorical items see \code{\link{reliability_omega_categorical}}.
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
#' \strong{Missing data and auxiliary variables.} By default incomplete
#' rows are listwise-deleted, which is unbiased only when the data are
#' missing completely at random and is inefficient always; the
#' \code{mlmr} vignette develops the argument at length. Setting
#' \code{missing = "fiml"} keeps every case with at least one observed
#' item and estimates by full information maximum likelihood, which is
#' consistent and efficient under the weaker missing at random (MAR)
#' assumption. The \code{aux} argument names auxiliary variables:
#' columns of \code{data} that are not part of the composite but are
#' correlated with the items or with the reasons values are missing.
#' They are entered as \emph{saturated correlates} (Graham, 2003):
#' correlated freely with each other and with every item's residual,
#' never loading on the factor and never entering the composite, so the
#' measurement model is undisturbed while FIML uses their information.
#' Beyond recovering information, a good auxiliary makes the MAR
#' assumption itself more plausible, since missingness that depends on
#' the auxiliary becomes MAR once the auxiliary is conditioned on
#' (Collins, Schafer, & Kam, 2001). Supplying \code{aux} implies
#' \code{missing = "fiml"}; combining it with an explicit
#' \code{missing = "listwise"} is an error. Listwise deletion remains
#' the default so no existing result changes and the missing-data
#' treatment is always a visible, deliberate choice. How each estimator
#' uses FIML: the model implied estimator simply fits its model with
#' \code{missing = "ml"}; the analytic estimator applies the classical
#' formula, unchanged, to the FIML estimate of the item covariance
#' matrix (from a saturated model over the items and any auxiliaries),
#' so the estimand stays the classical coefficient and only the
#' covariance matrix it is computed from improves. Under
#' \code{missing = "fiml"} the available intervals are \code{"ml"} and
#' \code{"ml_logistic"} (both estimators; the standard error comes from
#' the FIML information matrix), \code{"mlr"} and \code{"mlr_logistic"}
#' (model implied estimator; the Yuan-Bentler robust standard error),
#' and the bootstrap methods, which resample rows (including the
#' partially observed ones) and refit by FIML on each replication. The
#' complete-data closed forms (\code{"feldt"}, \code{"fisher"},
#' \code{"bonett"}, \code{"hakstian_whalen"}), \code{"adf"}, and
#' \code{"likelihood"} are errors with \code{missing = "fiml"} rather
#' than silently reverting to listwise deletion. Multiple imputation is
#' a different feature with a different interface and is out of scope
#' here.
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
#'   Efron, 1996). \code{"bootstrap_se"} uses the standard deviation
#'   of the bootstrap estimates as a standard error in an ordinary
#'   normal-theory interval; the \code{_logistic} variant builds that
#'   interval on the logit scale so the endpoints respect [0, 1].
#'   Replications on which the estimator cannot be computed (for
#'   example, a model refit that does not converge) are dropped, and
#'   the interval is computed from the replications that return a
#'   value. Requires raw data and the \pkg{boot} package. The default
#'   \code{B = 10000} is an accuracy choice: the BCa adjustment pushes
#'   the working quantiles farther into the tails of the bootstrap
#'   distribution than the percentile interval uses, and stabilizing
#'   them takes more replications than the customary 2000; reduce
#'   \code{B} for exploration, not for a reported analysis. Bootstrap
#'   results vary from run to run; supply \code{seed} for
#'   reproducibility.}
#'   \item{\code{"none"}}{Return only the point estimate.}
#' }
#'
#' \strong{Comparison with other packages.} The \pkg{psych} package
#' provides \code{\link[psych]{alpha}}, which reports coefficient
#' \eqn{\alpha} along with item-level diagnostics (item-total
#' correlations, alpha-if-item-deleted, and several alternative
#' coefficients). The emphasis in \pkg{psych} is broad exploratory
#' psychometric reporting. \code{reliability_alpha} in \pkg{DMAR}
#' differs in emphasis: it returns a single point estimate alongside a
#' principled confidence interval drawn from the methods compared in
#' Kelley and Pornprasertmanit (2016).
#'
#' @param data A numeric matrix or data frame of item scores (rows are
#'   respondents, columns are items). Either \code{data} or \code{S} must
#'   be supplied. How incomplete rows are handled is governed by
#'   \code{missing}: listwise deletion by default, or full information
#'   maximum likelihood with \code{missing = "fiml"}. When \code{aux} is
#'   supplied, \code{data} contains both the item columns and the
#'   auxiliary columns, and the columns not named in \code{aux} are the
#'   items.
#' @param S A symmetric covariance matrix among the items. If \code{S} is
#'   supplied, \code{N} must also be supplied; raw-data-only confidence
#'   interval methods (\code{"adf"}, \code{"bootstrap_*"}, \code{"percentile"},
#'   \code{"bca"}) are then unavailable, as are \code{missing = "fiml"}
#'   and \code{aux}, since a covariance matrix has no incomplete cases.
#' @param N Total sample size; required when \code{S} (rather than
#'   \code{data}) is supplied.
#' @param estimator How \eqn{\alpha} is estimated: \code{"analytic"}
#'   (the default) applies the closed-form equation to the observed
#'   covariance matrix, and \code{"model_implied"} takes the reliability
#'   implied by the \eqn{\tau}-equivalent single-factor model fit by
#'   maximum likelihood. See \emph{Details}.
#' @param missing How incomplete rows of \code{data} are handled:
#'   \code{"listwise"} (the default; complete-case analysis, the
#'   historical behavior) or \code{"fiml"} (full information maximum
#'   likelihood, using every case with at least one observed item).
#'   See \emph{Details}.
#' @param aux Optional character vector naming auxiliary variable
#'   columns of \code{data}, entered as saturated correlates under full
#'   information maximum likelihood. Supplying \code{aux} implies
#'   \code{missing = "fiml"}. See \emph{Details}.
#' @param ci_method Method for constructing the confidence interval. See
#'   \emph{Details} for which methods each estimator supports. The default
#'   \code{NULL} resolves to \code{"bonett"} for the analytic estimator and
#'   to \code{"mlr"} for the model implied estimator (or \code{"ml"} when
#'   only a covariance matrix is supplied, since \code{"mlr"} needs raw
#'   data). With \code{missing = "fiml"} the analytic default is
#'   \code{"ml"}, whose standard error comes from the FIML information
#'   matrix.
#' @param conf_level Confidence level for the interval (1 - Type I error
#'   rate). Defaults to \code{0.95}.
#' @param B Number of bootstrap replications when a bootstrap method is
#'   selected. Defaults to \code{10000}.
#' @param seed Random number seed used for bootstrap reproducibility.
#'   Defaults to \code{NULL}, which leaves the user's current RNG
#'   state intact; supply an integer for reproducibility. When set,
#'   the function saves and restores \code{.Random.seed} so the
#'   user's global RNG state is not polluted.
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
#' \code{"lower_limit"} and \code{"upper_limit"} (clamped to [0, 1]),
#' \code{"conf_level"}, \code{"N"} (the cases the analysis used: the
#' complete cases under listwise deletion, every case with at least one
#' observed item under \code{missing = "fiml"}), \code{"N_complete"}
#' (the complete cases, so the cost of listwise deletion is visible at
#' a glance; equal to \code{"N"} under listwise deletion and for
#' covariance input), and \code{"J"} (number of items). The selected
#' coefficient, CI method, and missing-data treatment travel as the
#' attributes \code{coefficient}, \code{ci_method}, \code{missing}, and
#' (when supplied) \code{aux}; bootstrap calls also record \code{B}.
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
#' # coefficient B times, so they are shown rather than run. The
#' # percentile interval reads its limits off the empirical quantiles of
#' # the bootstrap estimates:
#' #   reliability_alpha(data = items, ci_method = "percentile",
#' #                     B = 10000, seed = 113)
#' # The bias-corrected and accelerated interval adjusts those two
#' # quantile positions for median bias and for acceleration, and is the
#' # better choice for a reported interval:
#' #   reliability_alpha(data = items, ci_method = "bca", B = 10000,
#' #                     seed = 113)
#' # The default B = 10000 is an accuracy choice, not a formality, since
#' # BCa works farther into the tails of the bootstrap distribution than
#' # the percentile interval does (see Details).
#'
#' # Full information maximum likelihood with an auxiliary variable,
#' # shown rather than run because it fits a model over the items and
#' # the auxiliary. Missingness on y2 depends on an auxiliary z (missing
#' # at random given z), so listwise deletion is biased and FIML with z
#' # is not. Supplying aux implies missing = "fiml":
#' #   z <- eta + rnorm(200, sd = 0.5)
#' #   d <- data.frame(items, z = z)
#' #   d$y2[runif(200) < plogis(-1 + 1.5 * as.numeric(scale(z)))] <- NA
#' #   reliability_alpha(data = d, aux = "z")
#'
#' @concept Cronbach's alpha
#'
#' @keywords htest multivariate
#' @family reliability
#'
#' @export

reliability_alpha <- function(data = NULL, S = NULL, N = NULL,
                              estimator = c("analytic", "model_implied"),
                              missing = c("listwise", "fiml"),
                              aux = NULL,
                              ci_method = NULL,
                              conf_level = 0.95,
                              B = 10000,
                              seed = NULL) {
  estimator <- match.arg(estimator)
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

  # The two estimators support different interval methods, because a
  # closed-form interval for the sample coefficient and an interval read
  # off a fitted model are not interchangeable. Naming the sets here keeps
  # the contract in one place and lets the error say which estimator to
  # use rather than only that the value was wrong.
  shared  <- c("ml", "ml_logistic", "adf", "adf_logistic",
               "bootstrap_se", "bootstrap_se_logistic",
               "percentile", "bca", "none")
  ok_analytic <- c("bonett", "feldt", "fisher", "hakstian_whalen", shared)
  ok_model    <- c("mlr", "mlr_logistic", "likelihood", shared)
  ok <- if (estimator == "analytic") ok_analytic else ok_model

  if (is.null(ci_method)) {
    # Each estimator keeps the default it had as a separate function; the
    # model implied default needs raw data, so fall back when given only S.
    # Under missing = "fiml" the analytic default becomes "ml", whose
    # standard error comes from the FIML information matrix of the
    # saturated model; the complete-case closed forms (Bonett and its
    # siblings) are unavailable there.
    ci_method <- if (estimator == "analytic") {
      if (missing == "fiml") "ml" else "bonett"
    } else if (is.null(data)) {
      "ml"
    } else {
      "mlr"
    }
  }
  ci_method <- ci_method[1L]
  if (!is.character(ci_method) || !ci_method %in% c(ok_analytic, ok_model)) {
    stop("'ci_method' must be one of: ",
         paste(sQuote(unique(c(ok_analytic, ok_model))), collapse = ", "),
         ".", call. = FALSE)
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
    stop("ci_method = \"", ci_method, "\" is not available with ",
         "estimator = \"", estimator, "\". ", detail,
         "Use estimator = \"", other, "\", or choose one of: ",
         paste(sQuote(ok), collapse = ", "), ".", call. = FALSE)
  }

  if (!is.numeric(conf_level) || length(conf_level) != 1L ||
      conf_level <= 0 || conf_level >= 1) {
    stop("'conf_level' must be a single number in (0, 1).", call. = FALSE)
  }

  bootstrap_methods <- c("bootstrap_se", "bootstrap_se_logistic",
                         "percentile", "bca")
  if (missing == "fiml") {
    # Only intervals whose sampling theory can respect the missingness
    # pattern are offered; nothing silently falls back to listwise
    # deletion. The delta method ("ml"/"mlr") standard errors come from
    # the FIML information matrix, and the bootstrap resamples the
    # partially observed rows and refits by FIML. The remaining methods
    # are complete-data constructions: the Feldt, Fisher, Bonett, and
    # Hakstian-Whalen forms depend on a single complete-case N that
    # does not exist when rows are partially observed, ADF needs
    # complete raw-data fourth moments, and the profile likelihood
    # interval has not been extended to the FIML fit.
    if (ci_method %in% c("feldt", "fisher", "bonett", "hakstian_whalen")) {
      stop("ci_method = \"", ci_method, "\" is a complete-data closed ",
           "form: its standard error treats every case as fully ",
           "observed, which overstates the information in partially ",
           "observed rows. With missing = \"fiml\" use \"ml\" (or ",
           "\"ml_logistic\"), whose standard error comes from the FIML ",
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

  needs_raw <- ci_method %in% c("mlr", "mlr_logistic",
                                "adf", "adf_logistic", bootstrap_methods)
  if (needs_raw && is.null(inp$data)) {
    stop("ci_method = \"", ci_method, "\" requires raw 'data'; ",
         "a covariance matrix is not sufficient.", call. = FALSE)
  }

  if (estimator == "analytic") {
    if (missing == "fiml") {
      # The classical coefficient is still the classical coefficient,
      # computed by its own equation; only the covariance matrix it is
      # applied to changes, from the complete-case estimate to the FIML
      # estimate from the saturated model over the items and any
      # auxiliary variables. The "ml" interval reads the delta method
      # standard error of alpha off that same saturated fit, so it
      # reflects the actual missingness pattern.
      fc <- .fiml_cov(inp$data, items = inp$items, aux = inp$aux,
                      se = if (ci_method %in% c("ml", "ml_logistic")) {
                        "standard"
                      } else {
                        "none"
                      })
      if (!fc$converged) {
        stop("The saturated FIML model for the item covariance matrix ",
             "did not converge. This happens more often with small ",
             "samples, many auxiliary variables, or sparse missingness ",
             "patterns; consider dropping auxiliary variables.",
             call. = FALSE)
      }
      estimate <- .alpha_from_S(fc$S)
      crit <- stats::qnorm(1 - (1 - conf_level) / 2)
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
          data = inp$data,
          point_fn = function(d) {
            f <- .fiml_cov(d, items = inp$items, aux = inp$aux,
                           se = "none")
            if (!f$converged) NA_real_ else .alpha_from_S(f$S)
          },
          B = B, conf_level = conf_level, kind = ci_method, seed = seed
        )
      )
    } else {
      estimate <- .alpha_from_S(inp$S)
      ci <- switch(
        ci_method,
        none = list(se = NA_real_, lower = NA_real_, upper = NA_real_),
        feldt = .ci_feldt(estimate, inp$N, inp$J, conf_level),
        fisher = .ci_fisher(estimate, inp$N, conf_level),
        bonett = .ci_bonett(estimate, inp$N, inp$J, conf_level),
        hakstian_whalen = .ci_hakstian_whalen(estimate, inp$N, inp$J,
                                              conf_level),
        ml = .ci_alpha_ml(estimate, inp$S, inp$N, inp$J,
                          conf_level, logistic = FALSE),
        ml_logistic = .ci_alpha_ml(estimate, inp$S, inp$N, inp$J,
                                   conf_level, logistic = TRUE),
        adf = .ci_alpha_adf(estimate, inp$data, conf_level,
                            logistic = FALSE),
        adf_logistic = .ci_alpha_adf(estimate, inp$data, conf_level,
                                     logistic = TRUE),
        bootstrap_se = ,
        bootstrap_se_logistic = ,
        percentile = ,
        bca = .bootstrap_ci(
          data = inp$data,
          point_fn = function(d) .alpha_from_S(stats::cov(d)),
          B = B, conf_level = conf_level, kind = ci_method, seed = seed
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
      list(estimator = "MLR", se = "none")   # for the bootstrap families
    )
    if (missing == "fiml" && identical(fit_spec$se, "robust.sem")) {
      # With missing = "ml" the canonical robust standard error is the
      # Yuan-Bentler sandwich, lavaan's own default pairing for MLR
      # under FIML.
      fit_spec$se <- "robust.huber.white"
    }

    fit <- .omega_fit_cfa(data = inp$data, S = inp$S, N = inp$N,
                          equal_loadings = TRUE,
                          estimator = fit_spec$estimator,
                          se = fit_spec$se,
                          missing = if (missing == "fiml") "ml"
                                    else "listwise",
                          aux = inp$aux)
    if (!fit$converged) {
      stop("The tau-equivalent single-factor model did not converge with ",
           "estimator = '", fit_spec$estimator, "'.",
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
      ml = .ci_omega_delta(fit, conf_level, logistic = FALSE),
      ml_logistic = .ci_omega_delta(fit, conf_level, logistic = TRUE),
      mlr = .ci_omega_delta(fit, conf_level, logistic = FALSE),
      mlr_logistic = .ci_omega_delta(fit, conf_level, logistic = TRUE),
      adf = .ci_omega_delta(fit, conf_level, logistic = FALSE),
      adf_logistic = .ci_omega_delta(fit, conf_level, logistic = TRUE),
      likelihood = .ci_reliability_likelihood(
        data = inp$data, S = if (is.null(inp$data)) inp$S else NULL,
        N = inp$N, equal_loadings = TRUE, conf_level = conf_level
      ),
      bootstrap_se = ,
      bootstrap_se_logistic = ,
      percentile = ,
      bca = .bootstrap_ci(
        data = inp$data,
        point_fn = function(d) {
          bf <- .omega_fit_cfa(data = d, equal_loadings = TRUE,
                               estimator = "MLR", se = "none",
                               missing = if (missing == "fiml") "ml"
                                         else "listwise",
                               aux = inp$aux)
          if (!bf$converged) NA_real_ else bf$omega
        },
        B = B, conf_level = conf_level, kind = ci_method, seed = seed
      )
    )
  }

  out <- .relia_result(
    estimate = estimate, se = ci$se,
    lower = ci$lower, upper = ci$upper,
    conf_level = conf_level, N = inp$N, J = inp$J,
    coefficient = "alpha",
    ci_method = ci_method,
    B = if (ci_method %in% bootstrap_methods) B else NA_integer_,
    N_complete = inp$N_complete,
    missing = missing,
    aux = inp$aux,
    se_transformed = ci$se_transformed,
    se_transform_scale = ci$se_transform_scale
  )
  attr(out, "estimator") <- estimator
  out
}
