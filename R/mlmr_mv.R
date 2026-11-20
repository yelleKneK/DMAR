#' Multivariate Maximum Likelihood Regression With Full Information
#' Missing Data Handling
#'
#' Fits a multivariate multiple regression model by maximum
#' likelihood, with full information maximum likelihood handling of missing
#' values by default. Multiple outcomes are regressed on a shared
#' predictor set simultaneously, with the residual covariance among
#' outcomes estimated as part of the model. The formula interface
#' mirrors \code{\link[stats]{lm}}'s multivariate syntax,
#' \code{cbind(y1, y2, y3) ~ x1 + x2}, and the returned object
#' supports the same family of S3 methods as a univariate
#' \code{\link{mlmr}} fit.
#'
#' \strong{Why a separate function from \code{mlmr}.} A univariate
#' \code{mlmr} fit handles the case of one outcome regressed on one
#' or more predictors. \code{mlmr_mv} extends to the case of two or
#' more outcomes regressed on the same predictor set, modeling the
#' residual covariance among outcomes explicitly. Under FIML a row
#' that is missing on one outcome still contributes information about
#' the other outcomes (through the modeled residual covariance) and
#' about the joint distribution of the predictors, where listwise
#' deletion drops the row. A user who fits separate univariate regressions
#' for each outcome under listwise deletion can discard a great deal
#' of information when the outcomes are correlated and missingness
#' patterns differ.
#'
#' \strong{Same predictor set across outcomes.} The formula
#' \code{cbind(y1, y2) ~ x1 + x2} regresses both \code{y1} and
#' \code{y2} on \code{x1} and \code{x2}. Per-outcome predictor sets
#' (sometimes called seemingly unrelated regression with
#' heterogeneous predictors) are not supported; users who need that
#' can fit separate \code{\link{mlmr}} models or call
#' \code{\link[lavaan]{sem}} directly with a custom model string.
#'
#' \strong{Auxiliary variables.} As in \code{\link{mlmr}}, variables
#' that are not part of the regression but are correlated with an
#' outcome or with the missingness can be supplied through
#' \code{auxiliary} and are entered as saturated correlates (Graham,
#' 2003): each is correlated with every outcome's residual, every
#' predictor, and each other auxiliary, but never as a predictor, and
#' the predictors keep their free means, variances, and covariances. The
#' model over the outcomes, the predictors, and the auxiliaries is then
#' saturated, so the per-outcome coefficients keep their meaning (on
#' complete data they equal those of \code{lm()}) while the likelihood
#' draws on the auxiliaries' observed values (the inclusive analysis
#' strategy of Collins, Schafer, & Kam, 2001).
#'
#' \strong{The bootstrap interval.} \code{ci_method = "boot"}
#' resamples the rows of \code{data} with replacement \code{B} times
#' (1000 by default), refits the model on each resample, and reports
#' each coefficient's percentile interval. It is the interval to ask
#' for when the multivariate normality the likelihood assumes is
#' doubtful, since its coverage does not rest on that assumption. The
#' price is \code{B} refits of a model that already carries \emph{J}
#' outcomes, so a bootstrap interval is a deliberate request rather
#' than a default. Bootstrap results vary from run to run; supply
#' \code{boot_seed} for reproducibility. The mechanics, including the
#' Bollen-Stine variant, are given in the \code{ci_method} argument
#' description and in \code{\link{mlmr}}.
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
#' \strong{A \emph{t} reference distribution for the Wald interval.}
#' With \code{ci_method = "wald"}, \code{wald_df} refers every
#' coefficient's interval and \emph{p}-value to a \emph{t} distribution
#' rather than the normal one, keeping the estimate and its standard
#' error (von Hippel, 2016, Eq. 29), exactly as in \code{\link{mlmr}},
#' whose help page gives the reasoning. A number applies to every
#' coefficient of every outcome. Under \code{"missing_information"}
#' each coefficient gets von Hippel's \eqn{\max\{3, N(1 - \gamma) -
#' k\}} from its own fraction of missing information \eqn{\gamma}, with
#' \eqn{k} the number of regression parameters of its outcome (the
#' slopes and the intercept), so the coefficients of a fully observed
#' outcome keep nearly \eqn{N - k} degrees of freedom while those of an
#' outcome observed for few cases get fewer.
#'
#' \strong{Caveats.} Same as \code{\link{mlmr}}: the function
#' assumes that, conditional on the predictors, the joint
#' distribution of the outcomes is multivariate normal with constant
#' covariance, and that missingness is at most MAR. Factor
#' predictors and interactions are expanded through
#' \code{\link[stats]{model.matrix}} once and reused for every
#' outcome.
#'
#' @param formula A two-sided \code{\link[stats]{formula}} whose
#'   left-hand side is \code{cbind(y1, y2, ...)} for two or more
#'   numeric outcomes and whose right-hand side names the shared
#'   predictor set. Factor predictors, interactions, polynomial
#'   terms, and transformations are supported as in \code{mlmr}.
#' @param data A \code{data.frame} containing the variables in
#'   \code{formula}.
#' @param missing Character; passed to \pkg{lavaan}. Defaults to
#'   \code{"fiml"} (equivalently \code{"ml"}). See \code{\link{mlmr}}
#'   for the full list of accepted values.
#' @param ci_method Character; confidence interval method for the
#'   regression coefficients. \code{"profile"} (default), \code{"wald"},
#'   or \code{"boot"}. The same trade-offs apply as in
#'   \code{\link{mlmr}}; with \emph{J} outcomes and \emph{K}
#'   predictors, profile likelihood requires \eqn{O(JK)} constrained
#'   refits. \code{"boot"} resamples the rows of \code{data}
#'   \code{B} times (or draws Bollen-Stine model-based resamples),
#'   refits on each, drops resamples whose refit does not converge, and
#'   reports each coefficient's percentile interval, the empirical
#'   quantiles of its resampled estimates; as in \code{\link{mlmr}},
#'   the percentile interval is the only bootstrap interval offered.
#'   The examples ask for the Wald and bootstrap intervals so the help
#'   page stays quick; a reported analysis leaves \code{ci_method} at
#'   its default.
#' @param conf_level Desired level of confidence. Defaults to
#'   \code{0.95}.
#' @param B Number of bootstrap resamples when \code{ci_method = "boot"},
#'   a single whole number. Default \code{1000}, because every resample
#'   refits the model for every outcome. The fewest accepted is 400 at the
#'   default \code{conf_level} (200 at 0.90, and 2000 at 0.99, which the
#'   default does not meet); a smaller \code{B} is refused before anything
#'   is computed (class \code{dmar_count_invalid}). See \emph{Bootstrap
#'   replications} under \emph{Details} and \code{\link{dmar_bootstrap}}.
#' @param boot_type Character; \code{"ordinary"} (default) or
#'   \code{"bollen.stine"}.
#' @param boot_seed Integer or \code{NULL}. Defaults to \code{NULL},
#'   under which the resamples are drawn from, and advance, the current
#'   random number stream; supply an integer for reproducible
#'   bootstraps. When set, the seed applies for the duration of the call
#'   and the caller's generator state is restored on exit.
#' @param estimator Character; one of \code{"ML"} (default),
#'   \code{"MLR"}, \code{"MLM"}, \code{"GLS"}.
#' @param se Character or \code{NULL}; defaults to \code{NULL},
#'   meaning "choose automatically from \code{estimator}" (same
#'   logic as \code{\link{mlmr}}).
#' @param fixed_x Logical; defaults to \code{FALSE} (jointly model
#'   the predictor distribution, required for FIML to use rows with
#'   missing predictors).
#' @param auxiliary Character vector of variable names in \code{data}
#'   to include as auxiliary variables (saturated correlates; Graham,
#'   2003), or \code{NULL} (default) for none. Each auxiliary is
#'   correlated with every outcome's residual, every predictor, and
#'   each other auxiliary, but is never entered as a predictor, while the
#'   predictors keep their free covariances, so the model is saturated,
#'   the per-outcome regression coefficients keep their meaning, and the
#'   full information maximum likelihood draws on the auxiliaries' observed
#'   values (the inclusive analysis strategy; Collins, Schafer, & Kam,
#'   2001). A name in \code{auxiliary} must be numeric, present in
#'   \code{data}, absent from \code{formula}, and requires
#'   \code{fixed_x = FALSE}. Auxiliaries are refused with
#'   \code{missing = "listwise"}, which would drop every row on which an
#'   auxiliary is missing and so change the sample the coefficients
#'   describe, and with \code{missing = "pairwise"} or
#'   \code{"available.cases"}, which compute the moments from each pair of
#'   variables on its own, so that an auxiliary would leave every estimate
#'   unchanged.
#' @param effect_sizes Logical; whether to compute per-outcome
#'   standardized betas, semi-partial \eqn{R^2}, and Cohen's
#'   \eqn{f^2}. Defaults to \code{TRUE}.
#' @param enforce_es_bounds Logical; if \code{TRUE}, negative
#'   per-outcome \eqn{sr^2} or \eqn{f^2} estimates (finite-sample
#'   artifacts) are clamped to zero. Defaults to \code{FALSE}.
#' @param \dots Additional arguments forwarded to
#'   \code{\link[lavaan]{lavaan}}.
#' @param wald_df The reference distribution of the Wald interval and of
#'   the \emph{p}-values: \code{Inf} (the default), the normal
#'   distribution; a single positive number, the degrees of freedom of a
#'   \emph{t} distribution for every coefficient; or
#'   \code{"missing_information"}, von Hippel's (2016) degrees of freedom
#'   for each coefficient from its fraction of missing information (see
#'   \code{\link{mlmr}}). A finite value or the keyword requires
#'   \code{ci_method = "wald"}; the keyword is an error where, with
#'   values missing, \pkg{lavaan} fit the complete cases in place of the
#'   FIML fit asked for, as in \code{\link{mlmr}}. It follows
#'   \code{\dots}, so it must be named in full.
#'
#' @return An object of class \code{"mlmr_mv"}, a list with
#'   components similar to a univariate \code{\link{mlmr}} fit but
#'   extended for multiple outcomes:
#'   \describe{
#'     \item{\code{call}, \code{formula}, \code{terms},
#'       \code{model}, \code{xlevels}}{As in \code{mlmr}.}
#'     \item{\code{coefficients}}{A matrix with predictors (and an
#'       intercept row, when present) as rows and outcomes as
#'       columns, matching \code{coef.mlm}.}
#'     \item{\code{B}, \code{B_used}, \code{n_dropped},
#'       \code{B_required}, \code{bootstrap_status}}{Present only with
#'       \code{ci_method = "boot"}: the resamples requested, the number
#'       the intervals rest on, their difference, the fewest resamples a
#'       request may ask for, and the verdict
#'       (\code{"complete"}, \code{"dropped"}, \code{"exploratory"},
#'       \code{"insufficient_count"}, or \code{"insufficient_share"};
#'       under the last two the limits are \code{NA}). See
#'       \emph{Bootstrap replications} under Details.}
#'     \item{\code{coef_table}}{A long \code{data.frame} with one
#'       row per (outcome, term) combination; columns include
#'       \code{outcome}, \code{term}, \code{estimate}, \code{se},
#'       \code{test_statistic}, \code{p_value}, \code{ci_lower},
#'       \code{ci_upper}, \code{std_estimate}, carrying the confidence
#'       level in its \code{conf_level} attribute, so that
#'       \code{\link{results_sentence}} on a row reports it.}
#'     \item{\code{vcov}}{The variance-covariance matrix of the
#'       regression coefficients across all outcomes, returned by
#'       \code{vcov()}. Rows and columns follow the outcome-major
#'       order of \code{coef_table} (the column-major flattening of
#'       \code{coefficients}) and are named \code{"outcome:term"},
#'       for example \code{"mpg:wt"}, the naming
#'       \code{\link[stats]{vcov}} uses for an \code{"mlm"} fit. The
#'       cross-outcome blocks carry the sampling covariance between
#'       coefficients of different outcomes, so joint Wald tests
#'       across outcomes compose with \code{coef()}.}
#'     \item{\code{residual_cov}}{The estimated residual covariance
#'       matrix among outcomes (\emph{J} by \emph{J}).}
#'     \item{\code{R2}}{Named vector of model implied \eqn{R^2} per
#'       outcome.}
#'     \item{\code{adj_R2}}{Named vector of adjusted \eqn{R^2} per
#'       outcome.}
#'     \item{\code{effect_sizes}}{When \code{effect_sizes = TRUE},
#'       a long \code{data.frame} with one row per
#'       (outcome, predictor) combination giving \eqn{sr^2} and
#'       Cohen's \eqn{f^2}.}
#'     \item{\code{fitted.values}, \code{residuals}}{Matrices with
#'       rows = observations and columns = outcomes; \code{NA} in
#'       rows where any predictor is missing.}
#'     \item{\code{logLik}, \code{N}, \code{N_complete}}{As in
#'       \code{mlmr}.}
#'     \item{\code{lavaan_fit}}{The underlying \pkg{lavaan} fit.}
#'   }
#'   The object carries the attributes \code{wald_df} and, under
#'   \code{wald_df = "missing_information"}, \code{wald_df_table}, as an
#'   \code{\link{mlmr}} fit does, with the coefficients named
#'   \code{"outcome:term"} in the order of \code{coef_table}.
#'
#' @references
#' Collins, L. M., Schafer, J. L., & Kam, C.-M. (2001). A comparison of
#'   inclusive and restrictive strategies in modern missing data
#'   procedures. \emph{Psychological Methods, 6}(4), 330--351.
#'   \doi{10.1037/1082-989X.6.4.330}
#'
#' Graham, J. W. (2003). Adding missing-data-relevant variables to
#'   FIML-based structural equation models. \emph{Structural Equation
#'   Modeling, 10}(1), 80--100. \doi{10.1207/S15328007SEM1001_4}
#'
#' von Hippel, P. T. (2016). New confidence intervals and bias comparisons
#'   show that maximum likelihood can beat multiple imputation in small
#'   samples. \emph{Structural Equation Modeling, 23}(3), 422--437.
#'   \doi{10.1080/10705511.2015.1047931}
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso \code{\link{mlmr}} for the univariate sibling;
#'   \code{\link[stats]{lm}} (and the \code{"mlm"} object class) for
#'   the OLS multivariate analog; \code{\link[lavaan]{sem}} for the
#'   underlying engine.
#'
#' @examples
#' # Two outcomes on a shared predictor set. The residual covariance
#' # between the outcomes is estimated as part of the model, which is
#' # what separates this from two separate regressions. The fit asks
#' # for the Wald interval, a choice taken up below, and keeps the
#' # default effect_sizes = TRUE: the per-outcome effect sizes come
#' # back on the fit rather than in summary(), one row per outcome and
#' # predictor, giving the semi-partial R^2 and Cohen's f^2, and they
#' # are what fills the standardized coefficients in coef_table.
#' fit <- mlmr_mv(cbind(t6_paragraph_comprehension, t9_word_meaning) ~
#'                  t5_general_information + t7_sentence,
#'                data = holzinger_swineford,
#'                ci_method = "wald")
#' coef(fit)              # matrix: rows = predictors, cols = outcomes
#' summary(fit)
#' fit$R2                 # per-outcome R^2
#' fit$residual_cov       # residual covariance among outcomes
#' print(fit$effect_sizes, row.names = FALSE)
#'
#' # The interval menu is profile, Wald, and bootstrap. The default,
#' # ci_method = "profile", inverts the likelihood ratio test one
#' # coefficient at a time, and with two outcomes there are twice as
#' # many coefficients to profile; it is what a reported interval
#' # deserves, and leaving ci_method at its default asks for it. The
#' # fits on this page ask for the Wald or the bootstrap interval
#' # because those refits take longer than a help page should. The
#' # bootstrap resamples rows and takes percentile limits; it is what
#' # to ask for when the multivariate normality the likelihood assumes
#' # is doubtful. Every resample refits the two-outcome model. B = 10
#' # keeps the example quick and is below the 400 a reported interval
#' # needs; the option lets it run and the warning says so. A reported
#' # interval deserves the default B = 1000. boot_seed fixes the
#' # resamples, so the limits are reproducible rather than moving from
#' # run to run. The effect sizes cost one constrained refit per outcome
#' # and predictor, so this fit and the ones after it leave them off.
#' old <- options(dmar.bootstrap_exploratory = TRUE)
#' fit_boot <- mlmr_mv(cbind(t6_paragraph_comprehension, t9_word_meaning) ~
#'                       t5_general_information + t7_sentence,
#'                     data = holzinger_swineford,
#'                     ci_method = "boot", B = 10, boot_seed = 113,
#'                     effect_sizes = FALSE)
#' options(old)
#' confint(fit_boot)
#'
#' # FIML versus listwise when one outcome has missing values. The
#' # revised second-form test t26_flags was administered to only 145
#' # of the 301 students, so it carries real missingness. A row with
#' # t26_flags missing still informs the likelihood about the other
#' # outcome, about the predictors, and, through the residual
#' # covariance, about t26_flags itself, so no row is discarded. Notice
#' # the two sample sizes, and that the coefficients of the complete
#' # outcome differ between the fits: listwise deletion drops 156 of
#' # its observed rows along with the missing t26_flags values.
#' fit_fiml <- mlmr_mv(cbind(t6_paragraph_comprehension,
#'                           t26_flags) ~
#'                       t7_sentence + t9_word_meaning,
#'                     data = holzinger_swineford,
#'                     ci_method = "wald", effect_sizes = FALSE)
#' fit_lwd  <- mlmr_mv(cbind(t6_paragraph_comprehension,
#'                           t26_flags) ~
#'                       t7_sentence + t9_word_meaning,
#'                     data = holzinger_swineford,
#'                     missing = "listwise", ci_method = "wald",
#'                     effect_sizes = FALSE)
#' c(N_fiml = nobs(fit_fiml), N_listwise = nobs(fit_lwd))
#' cbind(FIML = coef(fit_fiml)[, "t6_paragraph_comprehension"],
#'       listwise = coef(fit_lwd)[, "t6_paragraph_comprehension"])
#'
#' # Auxiliary variable (saturated correlates): the complete speed test
#' # t13_straight_and_curved_capitals informs the likelihood without
#' # entering either regression. Continuing from the model above, the
#' # coefficients of the complete outcome are unchanged to working
#' # precision, while those of t26_flags move, since the auxiliary
#' # carries information about the rows where t26_flags is missing.
#' fit_aux <- mlmr_mv(cbind(t6_paragraph_comprehension,
#'                          t26_flags) ~
#'                      t7_sentence + t9_word_meaning,
#'                    data = holzinger_swineford,
#'                    ci_method = "wald",
#'                    auxiliary = "t13_straight_and_curved_capitals",
#'                    effect_sizes = FALSE)
#' coef(fit_aux)
#'
#' @keywords regression models multivariate
#'
#' @export
mlmr_mv <- function(formula,
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

  .mlmr_require_lavaan("mlmr_mv")

  if (missing(data) || !is.data.frame(data)) {
    .dmar_refuse_argument("'data' must be a data.frame.", arg = "data")
  }
  if (!inherits(formula, "formula") || length(formula) != 3L) {
    .dmar_refuse_argument(
      "'formula' must be a two-sided formula.", arg = "formula")
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
  wald_df <- .wald_df_check(wald_df, "mlmr_mv")
  .wald_df_require_wald(wald_df, ci_method == "wald", paste0(
    "and ci_method = \"", ci_method, "\" is not one; wald_df applies ",
    "with ci_method = \"wald\", whose interval and p-value it refers to ",
    "the t distribution together."))
  missing_lavaan <- if (missing == "fiml") "ml" else missing
  if (is.null(se)) {
    se <- switch(estimator,
                 ML  = "standard", GLS = "standard",
                 MLR = "robust.huber.white", MLM = "robust.sem")
  }

  design <- mlmr_mv_build_design(formula, data)
  Y_names      <- design$Y_names
  Y_display    <- design$Y_display
  X_names      <- design$X_names
  X_display    <- design$X_display
  has_intercept <- design$has_intercept
  model_data   <- design$model_data
  model_frame  <- design$model_frame
  trms         <- design$terms

  J <- length(Y_names)
  K <- length(X_names)
  if (J < 2L) {
    .dmar_refuse_argument(c(
      "mlmr_mv() requires two or more outcomes. ",
      "For a single outcome, use mlmr() instead."), arg = NULL)
  }
  if (K == 0L) {
    .dmar_refuse_argument(c(
      "The right-hand side of 'formula' has no predictors after ",
      "model.matrix expansion. mlmr_mv requires at least one ",
      "predictor."), arg = "formula")
  }

  syntax <- mlmr_mv_build_syntax(Y_names, X_names, has_intercept)

  # Auxiliary variables (saturated correlates; Graham, 2003), shared
  # with mlmr(). Each auxiliary correlates with every outcome's
  # residual, every predictor, and every other auxiliary, but enters
  # no regression, and the predictors keep their covariances, so the
  # model stays saturated, the focal coefficients keep their meaning, and
  # FIML draws on the auxiliaries' rows.
  aux <- mlmr_prepare_auxiliary(auxiliary, data,
                                model_vars = all.vars(formula),
                                fixed_x = fixed_x, missing = missing)
  fit_data <- model_data
  if (!is.null(aux)) {
    syntax <- paste(syntax,
                    mlmr_auxiliary_syntax(Y_names, X_names, aux$aux_names),
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
    .dmar_refuse("lavaan reports the model did not converge.", "not_converged")
  }

  pe <- lavaan::parameterEstimates(fit, ci = FALSE)
  full_vcov <- lavaan::vcov(fit)

  # Build per-outcome label sets and pull the estimates.
  beta_labels_all <- character(0)
  est_mat <- matrix(NA_real_, nrow = K + as.integer(has_intercept),
                    ncol = J,
                    dimnames = list(
                      c(if (has_intercept) "(Intercept)", X_display),
                      Y_display))
  se_mat <- est_mat
  for (j in seq_len(J)) {
    y_j <- Y_names[j]
    labs <- mlmr_mv_beta_labels(y_j, K, has_intercept)
    beta_labels_all <- c(beta_labels_all, labs)
    for (i in seq_along(labs)) {
      row <- pe[pe$label == labs[i], , drop = FALSE]
      if (nrow(row) == 1L) {
        est_mat[i, j] <- row$est[[1L]]
        se_mat[i, j]  <- row$se[[1L]]
      }
    }
  }
  z_mat <- est_mat / se_mat
  # The reference distribution of the Wald ratios: the normal one by
  # default, or t on wald_df degrees of freedom, a number for every
  # coefficient or, under "missing_information", von Hippel's (2016)
  # max(3, N(1 - FMI) - k) for each coefficient, with k the number of
  # regression parameters of its outcome (.wald_df_resolve(),
  # R/missing_information_internals.R). The coefficients are taken in
  # outcome-major order, the order of beta_labels_all and coef_table.
  # A message of the computation of the missing information names the
  # variables as the formula does (.fmi_with_names()).
  wald <- .fmi_with_names(
    .wald_df_resolve(
      wald_df, "mlmr_mv", fit = fit,
      rows = pe[match(beta_labels_all, pe$label), c("lhs", "op", "rhs")],
      context = list(n_coef = nrow(est_mat),
                     terms = paste0(rep(Y_display, each = nrow(est_mat)),
                                    ":", rownames(est_mat)),
                     requested_missing = if (missing %in% c("fiml", "ml")) {
                       missing
                     })),
    stats::setNames(c(Y_display, X_display, aux$aux_display),
                    c(Y_names, X_names, aux$aux_names)))
  if (!is.null(wald$table)) .df_require_finite(wald$df, "mlmr_mv")
  df_mat <- matrix(unname(wald$df), nrow = nrow(est_mat),
                   ncol = ncol(est_mat), dimnames = dimnames(est_mat))
  p_mat <- .wald_p(z_mat, if (.wald_df_is_normal(wald_df)) Inf else df_mat)

  # Coefficient-only covariance. lavaan's vcov covers every model
  # parameter (residual variances and covariances, and the predictor
  # moments when fixed_x = FALSE); subset to the regression
  # coefficients in the outcome-major order of beta_labels_all, which
  # is the order of coef_table and of the column-major flattening of
  # the coefficient matrix. Name rows and columns "outcome:term", the
  # convention stats::vcov uses for an "mlm" fit.
  vcov_idx <- match(beta_labels_all, colnames(full_vcov))
  vcov_beta <- full_vcov[vcov_idx, vcov_idx, drop = FALSE]
  vcov_names <- paste0(rep(Y_display, each = nrow(est_mat)), ":",
                       rownames(est_mat))
  dimnames(vcov_beta) <- list(vcov_names, vcov_names)

  ci_array <- mlmr_mv_compute_ci(
    fit = fit, syntax = syntax, fit_args = fit_args,
    beta_labels_all = beta_labels_all, est_mat = est_mat,
    se_mat = se_mat, conf_level = conf_level,
    ci_method = ci_method, B = B,
    boot_type = boot_type, boot_seed = boot_seed,
    wald_df = if (.wald_df_is_normal(wald_df)) Inf else df_mat
  )
  boot_counts <- mlmr_boot_counts(attr(ci_array, "bootstrap"))
  attr(ci_array, "bootstrap") <- NULL

  # Implied covariance and per-outcome residual / variance.
  implied_cov <- lavaan::lavInspect(fit, "cov.ov")
  sigma2_y <- vapply(Y_names, function(y) implied_cov[y, y],
                     numeric(1L))
  sigma2_e <- vapply(Y_names, function(y) {
    lab <- paste0("sigma2_e_", y)
    val <- pe$est[pe$label == lab]
    if (length(val) == 0L) NA_real_ else val[1L]
  }, numeric(1L))
  R2 <- 1 - sigma2_e / sigma2_y
  N_used <- as.integer(lavaan::lavInspect(fit, "ntotal"))
  N_complete <- sum(stats::complete.cases(model_data))
  if (N_complete > K + 1L) {
    adj_R2 <- 1 - (1 - R2) * (N_complete - 1L) / (N_complete - K - 1L)
  } else {
    adj_R2 <- rep(NA_real_, J)
  }
  names(R2) <- names(adj_R2) <- Y_display

  # Residual covariance matrix among outcomes.
  resid_cov <- matrix(NA_real_, J, J,
                      dimnames = list(Y_display, Y_display))
  diag(resid_cov) <- sigma2_e
  if (J >= 2L) {
    for (a in seq_len(J - 1L)) {
      for (b in seq.int(a + 1L, J)) {
        lab <- paste0("sigma_e_", Y_names[a], "_", Y_names[b])
        val <- pe$est[pe$label == lab]
        if (length(val) > 0L) {
          resid_cov[a, b] <- val[1L]
          resid_cov[b, a] <- val[1L]
        }
      }
    }
  }

  ll_val <- as.numeric(lavaan::logLik(fit))
  npar <- length(lavaan::coef(fit))
  attr(ll_val, "df") <- npar
  attr(ll_val, "nobs") <- N_used
  class(ll_val) <- "logLik"

  # Effect sizes per outcome.
  std_mat <- matrix(NA_real_, nrow = nrow(est_mat), ncol = J,
                    dimnames = dimnames(est_mat))
  effect_size_table <- NULL
  if (effect_sizes) {
    SD_Y <- sqrt(sigma2_y)
    SD_X <- sqrt(diag(implied_cov)[X_names])
    slope_rows <- if (has_intercept)
                    seq.int(2L, nrow(est_mat))
                  else seq_len(nrow(est_mat))
    for (j in seq_len(J)) {
      std_mat[slope_rows, j] <- est_mat[slope_rows, j] * SD_X / SD_Y[j]
    }

    es <- mlmr_mv_compute_effect_sizes(
      fit_args = fit_args, syntax = syntax,
      Y_names = Y_names, Y_display = Y_display,
      X_display = X_display, K = K,
      sigma2_y = sigma2_y, full_R2 = R2
    )
    if (enforce_es_bounds) {
      neg_sr2 <- !is.na(es$sr2) & es$sr2 < 0
      neg_f2  <- !is.na(es$f2)  & es$f2  < 0
      es$sr2[neg_sr2] <- 0
      es$f2[neg_f2]   <- 0
      attr(es, "clamped") <- neg_sr2 | neg_f2
    } else {
      attr(es, "clamped") <- rep(FALSE, nrow(es))
    }
    attr(es, "enforce_es_bounds") <- enforce_es_bounds
    effect_size_table <- es
  }

  # Long coefficient table.
  coef_table <- do.call(rbind, lapply(seq_len(J), function(j) {
    data.frame(
      outcome      = Y_display[j],
      term         = rownames(est_mat),
      estimate     = unname(est_mat[, j]),
      std_estimate = unname(std_mat[, j]),
      se           = unname(se_mat[, j]),
      test_statistic      = unname(z_mat[, j]),
      p_value      = unname(p_mat[, j]),
      ci_lower     = unname(ci_array[, 1L, j]),
      ci_upper     = unname(ci_array[, 2L, j]),
      stringsAsFactors = FALSE,
      row.names = NULL
    )
  }))
  # The table records the confidence level of its intervals, which a row
  # keeps, so that results_sentence() on a coefficient reports it.
  attr(coef_table, "conf_level") <- conf_level

  # Fitted values and residuals as matrices.
  X_mat_internal <- as.matrix(model_data[, X_names, drop = FALSE])
  any_X_missing <- !stats::complete.cases(X_mat_internal)
  fitted_mat <- matrix(NA_real_, nrow = nrow(model_data), ncol = J,
                       dimnames = list(rownames(model_data), Y_display))
  for (j in seq_len(J)) {
    slope_vals <- est_mat[if (has_intercept)
                           -1L else seq_len(nrow(est_mat)), j]
    intercept_val <- if (has_intercept) est_mat[1L, j] else 0
    fitted_mat[!any_X_missing, j] <- intercept_val +
      X_mat_internal[!any_X_missing, , drop = FALSE] %*% slope_vals
  }
  Y_mat <- as.matrix(model_data[, Y_names, drop = FALSE])
  colnames(Y_mat) <- Y_display
  resid_mat <- Y_mat - fitted_mat

  xlevels <- stats::.getXlevels(trms, model_frame)

  structure(c(list(
    call          = call,
    formula       = formula,
    terms         = trms,
    model         = model_frame,
    xlevels       = xlevels,
    coefficients  = est_mat,
    vcov          = vcov_beta,
    ci            = ci_array,
    ci_method     = ci_method,
    conf_level    = conf_level,
    coef_table    = coef_table,
    residual_cov  = resid_cov,
    sigma2        = setNames(sigma2_e, Y_display),
    R2            = R2,
    adj_R2        = adj_R2,
    effect_sizes  = effect_size_table,
    logLik        = ll_val,
    N             = N_used,
    N_complete    = N_complete,
    fitted.values = fitted_mat,
    residuals     = resid_mat,
    missing       = missing,
    estimator     = estimator,
    se_type       = se,
    fixed_x       = fixed_x,
    auxiliary     = if (!is.null(aux)) aux$aux_display else NULL,
    Y_names       = Y_names,
    Y_display     = Y_display,
    X_names       = X_names,
    X_display     = X_display,
    has_intercept = has_intercept,
    syntax        = syntax,
    model_data    = model_data,
    fit_args      = fit_args,
    lavaan_fit    = fit
  ), boot_counts), class = "mlmr_mv",
  wald_df = if (is.null(wald$table)) wald_df else wald$df,
  wald_df_table = wald$table)
}


#' @export
print.mlmr_mv <- function(x, digits = max(3L, getOption("digits") - 3L),
                          ...) {
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
      "\n\n", sep = "")
  cat("Coefficients (rows = predictors, columns = outcomes):\n")
  print.default(format(x$coefficients, digits = digits),
                print.gap = 2L, quote = FALSE)
  cat("\nResidual covariance among outcomes:\n")
  print(round(x$residual_cov, digits))
  cat("\nR^2 per outcome:\n")
  print(round(x$R2, digits))
  if (!is.null(x$bootstrap_status)) cat("\n")
  mlmr_boot_line(x)
  invisible(x)
}


#' @export
summary.mlmr_mv <- function(object, ...) {
  ans <- list(
    call         = object$call,
    formula      = object$formula,
    coef_table   = object$coef_table,
    effect_sizes = object$effect_sizes,
    residual_cov = object$residual_cov,
    R2           = object$R2,
    adj_R2       = object$adj_R2,
    sigma2       = object$sigma2,
    N            = object$N,
    N_complete   = object$N_complete,
    logLik       = object$logLik,
    conf_level   = object$conf_level,
    ci_method    = object$ci_method,
    missing      = object$missing,
    estimator    = object$estimator,
    se_type      = object$se_type,
    auxiliary    = object$auxiliary,
    Y_display    = object$Y_display,
    wald_df      = attr(object, "wald_df"),
    wald_df_table = attr(object, "wald_df_table")
  )
  ans <- c(ans, object[intersect(c("B", "B_used", "n_dropped", "B_required",
                                   "bootstrap_status"), names(object))])
  class(ans) <- "summary.mlmr_mv"
  ans
}


#' @export
print.summary.mlmr_mv <- function(x,
                                  digits = max(3L,
                                               getOption("digits") - 3L),
                                  signif.stars = getOption("show.signif.stars"),
                                  ...) {
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
      "\n\n", sep = "")
  cat("Missing: ", x$missing, " | Estimator: ", x$estimator,
      " | SE: ", x$se_type, "\n", sep = "")
  if (!is.null(x$auxiliary)) {
    cat("Auxiliary variable(s), saturated correlates: ",
        paste(x$auxiliary, collapse = ", "), "\n", sep = "")
  }
  cat("Sample size used: ", x$N, "   Complete cases: ",
      x$N_complete, "\n\n", sep = "")
  for (yi in x$Y_display) {
    cat("--- Outcome: ", yi, " ---\n", sep = "")
    rows <- x$coef_table[x$coef_table$outcome == yi, , drop = FALSE]
    coef_mat <- as.matrix(rows[, c("estimate", "se", "test_statistic",
                                   "p_value")])
    rownames(coef_mat) <- rows$term
    # t headers when wald_df is finite, as in print.summary.mlmr().
    colnames(coef_mat) <- .wald_df_coef_headers(x$wald_df)
    stats::printCoefmat(coef_mat, digits = digits,
                        signif.stars = signif.stars,
                        has.Pvalue = TRUE, ...)
    cat("R^2:        ", format(x$R2[yi], digits = digits), "\n",
        sep = "")
    cat("Adj. R^2:   ", format(x$adj_R2[yi], digits = digits),
        "\n\n", sep = "")
  }
  wd_lines <- .wald_df_summary_lines(x$wald_df, x$wald_df_table)
  if (length(wd_lines)) {
    cat(paste(wd_lines, collapse = "\n"), "\n\n", sep = "")
  }
  cat("Residual covariance among outcomes:\n")
  print(round(x$residual_cov, digits))
  cat("\nLog likelihood: ", format(as.numeric(x$logLik), digits = digits),
      "   AIC: ", format(stats::AIC(x$logLik), digits = digits),
      "   BIC: ", format(stats::BIC(x$logLik), digits = digits), "\n",
      sep = "")
  mlmr_boot_line(x)
  invisible(x)
}


#' @export
coef.mlmr_mv <- function(object, ...) object$coefficients


#' @export
vcov.mlmr_mv <- function(object, ...) object$vcov


#' @export
nobs.mlmr_mv <- function(object, ...) object$N


#' @export
logLik.mlmr_mv <- function(object, ...) object$logLik


#' @export
fitted.mlmr_mv <- function(object, ...) object$fitted.values


#' @export
residuals.mlmr_mv <- function(object, ...) object$residuals


#' @export
formula.mlmr_mv <- function(x, ...) x$formula


#' @export
confint.mlmr_mv <- function(object, parm, level = 0.95, ...) {
  if (!missing(level) &&
      !isTRUE(all.equal(level, object$conf_level))) {
    .dmar_warn(c(
      "confint() ignores 'level' for mlmr_mv fits; mlmr_mv ",
      "returns the interval at the conf_level used at fit ",
      "time (", object$conf_level, ")."), "argument_ignored")
  }
  # Return a long-format data.frame: outcome, term, lower, upper.
  out <- object$coef_table[, c("outcome", "term", "ci_lower",
                               "ci_upper"), drop = FALSE]
  # Labeled as stats::confint.default() labels its columns, whatever
  # getOption("digits") is (.ci_quantile_labels()).
  colnames(out) <- c("outcome", "term",
                     .ci_quantile_labels(object$conf_level))
  if (!missing(parm)) {
    out <- out[out$term %in% parm | out$outcome %in% parm, ,
               drop = FALSE]
  }
  out
}


#' @export
predict.mlmr_mv <- function(object, newdata, ...) {
  if (missing(newdata) || is.null(newdata)) {
    return(object$fitted.values)
  }
  trms <- stats::delete.response(object$terms)
  mf <- stats::model.frame(trms, newdata, na.action = stats::na.pass,
                           xlev = object$xlevels)
  X <- stats::model.matrix(trms, mf)
  # X column order matches the names in object$coefficients rows.
  X[, rownames(object$coefficients), drop = FALSE] %*%
    object$coefficients
}


#' @export
as.data.frame.mlmr_mv <- function(x, ...) {
  x$coef_table
}


#' Compare Nested Multivariate FIML Regression Fits
#'
#' Compares two or more nested \code{\link{mlmr_mv}} fits with the
#' likelihood ratio test, delegating the chi square computation to
#' \code{\link[lavaan]{lavTestLRT}}. The models must be nested (every
#' predictor in the more restricted model is also in the more general
#' one) and must be fit to the same outcomes and the same data with
#' the same missing data handling. This is the multivariate
#' counterpart of \code{\link{anova.mlmr}}: because every outcome is
#' regressed on the shared predictor set, dropping a predictor drops
#' its slope on every outcome, so the nesting constraint sets that
#' slope to zero across all outcomes at once.
#'
#' @details
#' "The same data" means the same observations, with the variables the
#' fits share holding the same values. It does not mean the same number
#' of complete cases. Under \code{missing = "fiml"} a legitimate nested
#' comparison routinely has different complete-case counts: adding a
#' predictor that is itself incompletely observed lowers the number of
#' rows complete on every modeled variable, yet the observations, and
#' the information each of them contributes to the likelihood, are
#' unchanged. That comparison is accepted. Fits run on genuinely
#' different data are refused.
#'
#' @param object An \code{mlmr_mv} fit.
#' @param \dots Additional \code{mlmr_mv} fits, nested with respect to
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
#' fit1 <- mlmr_mv(cbind(t6_paragraph_comprehension, t9_word_meaning) ~
#'                   t5_general_information,
#'                 data = holzinger_swineford,
#'                 ci_method = "wald", effect_sizes = FALSE)
#' fit2 <- mlmr_mv(cbind(t6_paragraph_comprehension, t9_word_meaning) ~
#'                   t5_general_information + t7_sentence,
#'                 data = holzinger_swineford,
#'                 ci_method = "wald", effect_sizes = FALSE)
#' anova(fit1, fit2)
#'
#' @export
anova.mlmr_mv <- function(object, ...) {
  fits <- c(list(object), list(...))
  if (length(fits) < 2L) {
    .dmar_refuse_argument(
      "anova.mlmr_mv requires at least two mlmr_mv fits to compare.",
      arg = NULL)
  }
  if (!all(vapply(fits, inherits, logical(1L), "mlmr_mv"))) {
    .dmar_refuse_argument(
      "All arguments to anova.mlmr_mv must be of class 'mlmr_mv'.", arg = NULL)
  }

  # Cross-fit sanity checks. Likelihoods are only comparable when the
  # fits are over the same outcomes and data with the same missing
  # data treatment; otherwise the chi square statistic is meaningless.
  # The same-data requirement is enforced further down, on the analysis
  # data itself, once the fits have been ordered and the largest one
  # identified.
  miss <- vapply(fits, function(f) f$missing, character(1L))
  if (length(unique(miss)) != 1L) {
    .dmar_refuse_argument(c(
      "anova.mlmr_mv requires fits to use the same 'missing' ",
      "setting; got (", paste(miss, collapse = ", "),
      "). Likelihoods under different missing-data treatments ",
      "are not comparable."), arg = "missing")
  }
  est <- vapply(fits, function(f) f$estimator, character(1L))
  if (length(unique(est)) != 1L) {
    .dmar_refuse_argument(c(
      "anova.mlmr_mv requires fits to use the same 'estimator' ",
      "(got ", paste(est, collapse = ", "), ")."), arg = "estimator")
  }
  ynm <- vapply(fits, function(f) paste(f$Y_display, collapse = ","),
                character(1L))
  if (length(unique(ynm)) != 1L) {
    .dmar_refuse_argument(c(
      "anova.mlmr_mv requires fits to share the same outcomes: ",
      "got (", paste(ynm, collapse = " | "), ")."), arg = NULL)
  }

  # Order fits from fewest predictors to most. The largest serves as
  # the common variable set; smaller fits are re-expressed as
  # constrained versions of the largest so the LR test compares
  # likelihoods over a single joint observed variable set, which is
  # required for lavTestLRT to be valid when fixed_x = FALSE.
  K_each <- vapply(fits, function(f) length(f$X_display), integer(1L))
  fits <- fits[order(K_each)]
  largest <- fits[[length(fits)]]
  n_out <- length(largest$Y_display)

  for (f in fits) {
    if (!all(f$X_display %in% largest$X_display)) {
      .dmar_refuse_argument(c(
        "anova.mlmr_mv requires nested fits: the predictors in every ",
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
  #      contribute, untouched. That comparison is exactly what FIML
  #      exists to support, so refusing it on a complete-case count
  #      would be wrong.
  #
  # So compare the analysis data directly. Each fit stores it in
  # $model_data with internal, positional column names (y_1, x_1, ...)
  # and with missing values preserved, one row per row of 'data';
  # relabeling by the display names lets fits be compared on the
  # outcomes and predictors they share regardless of predictor order.
  relabel <- function(f) {
    md <- f$model_data
    names(md) <- c(f$Y_display, f$X_display)
    md
  }
  largest_data <- relabel(largest)
  for (f in fits[-length(fits)]) {
    fdat <- relabel(f)
    if (nrow(fdat) != nrow(largest_data)) {
      .dmar_refuse_argument(c(
        "anova.mlmr_mv requires all fits to be on the same data. The ",
        "fit '", paste(deparse(f$formula), collapse = " "), "' was run ",
        "on ", nrow(fdat), " observations and the largest fit '",
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
        "anova.mlmr_mv requires all fits to be on the same data. The ",
        "fit '", paste(deparse(f$formula), collapse = " "), "' and the ",
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
    # Each dropped predictor loses its slope on every outcome, so the
    # nesting constraint zeroes b_y_<j>_<k> for outcome j = 1..n_out
    # and predictor position k. The slope labels follow the mlmr_mv
    # lavaan syntax (see the model string in mlmr_mv()).
    extra_constraints <- as.vector(outer(
      seq_len(n_out), slope_positions,
      function(j, k) paste0("b_y_", j, "_", k, " == 0")
    ))
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

  # model.names must be supplied explicitly; see the matching comment
  # in anova.mlmr(): without it lavTestLRT deparses the inlined lavaan
  # objects for display names, which stops the call on platforms where
  # that deparse spans multiple lines.
  out <- do.call(lavaan::lavTestLRT,
                 c(aligned_fits,
                   list(model.names = paste0("fit", seq_along(aligned_fits)))))
  # lavTestLRT sorts its rows by model degrees of freedom, which need not match
  # the order the fits were supplied in, and it writes deparsed lavaan fit
  # objects into the row names (hundreds of kilobytes at print time). Map each
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
    "Likelihood ratio test for nested mlmr_mv fits",
    vapply(seq_len(nrow(out)), function(j) {
      i <- row_fit[j]
      sprintf("Model %d: %s", i,
              paste(deparse(fits[[i]]$formula), collapse = " "))
    }, character(1L))
  )
  out
}


#' A Multivariate FIML Regression Fit
#'
#' Broom-style \code{tidy()} and \code{glance()} for \code{\link{mlmr_mv}}
#' fits. \code{tidy()} returns one row per coefficient per outcome with the
#' broom dotted columns plus a leading \code{response} column identifying the
#' outcome; \code{glance()} returns a one-row model-level summary with the
#' per-outcome \eqn{R^2} averaged and the number of responses reported.
#'
#' @param x An \code{mlmr_mv} fit.
#' @param conf.int Logical: include \code{ci_lower} / \code{ci_upper},
#'   with a \code{conf_level} column holding the level of the intervals?
#' @param conf_level Ignored (the interval level is fixed at fit time and
#'   stored on the object); present for broom signature compatibility.
#' @param standardized Logical: include \code{std_estimate}?
#' @param ... Unused.
#'
#' @return For \code{tidy.mlmr_mv}, a \code{data.frame} with columns
#'   \code{response}, \code{term}, \code{estimate}, \code{se},
#'   \code{statistic}, \code{p_value}, \code{wald_df} (the degrees of
#'   freedom of the \emph{t} reference of the row's statistic and
#'   \emph{p}-value, \code{Inf} for the normal reference), and optionally
#'   \code{ci_lower}, \code{ci_upper}, \code{std_estimate}. For
#'   \code{glance.mlmr_mv}, a
#'   one-row \code{data.frame} with \code{R2} (mean across
#'   outcomes), \code{df}, \code{logLik}, \code{AIC}, \code{BIC},
#'   \code{deviance}, \code{nobs}, and \code{n.responses}.
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso \code{\link{mlmr_mv}}; \code{\link{mlmr}} for the univariate
#'   methods these mirror.
#'
#' @importFrom generics tidy
#' @exportS3Method generics::tidy
#' @name tidy.mlmr_mv
tidy.mlmr_mv <- function(x, conf.int = FALSE, conf_level = NULL,
                         standardized = FALSE, ...) {
  tab <- x$coef_table
  out <- data.frame(
    response  = tab$outcome,
    term      = tab$term,
    estimate  = tab$estimate,
    se = tab$se,
    statistic = tab$test_statistic,
    p_value   = tab$p_value,
    # The reference of the statistic and p-value of each row, which a stack
    # of tables would otherwise lose (.wald_df_column()); under the keyword
    # the record is named "outcome:term" in the order of the table.
    wald_df   = .wald_df_column(attr(x, "wald_df", exact = TRUE), nrow(tab),
                                paste0(tab$outcome, ":", tab$term)),
    stringsAsFactors = FALSE
  )
  if (isTRUE(conf.int)) {
    out$ci_lower  <- tab$ci_lower
    out$ci_upper <- tab$ci_upper
    # The level the intervals were computed at, so that a row of the tidy
    # table says what level its limits belong to.
    out$conf_level <- rep(x$conf_level %||% NA_real_, nrow(out))
  }
  if (isTRUE(standardized)) out$std_estimate <- tab$std_estimate
  out
}

#' @rdname tidy.mlmr_mv
#' @importFrom generics glance
#' @exportS3Method generics::glance
glance.mlmr_mv <- function(x, ...) {
  ll   <- as.numeric(x$logLik)
  npar <- attr(x$logLik, "df")
  data.frame(
    R2   = mean(x$R2),
    df          = npar,
    logLik      = ll,
    AIC         = stats::AIC(x$logLik),
    BIC         = stats::BIC(x$logLik),
    deviance    = -2 * ll,
    nobs        = x$N,
    n.responses = length(x$R2),
    stringsAsFactors = FALSE
  )
}
