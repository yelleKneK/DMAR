# Confidence interval for generalized eta squared (effect size for ANOVA).
#' Confidence Interval for Generalized Eta Squared (Approximate)
#'
#' Returns the point estimate of generalized eta squared (\eqn{\eta^2_G};
#' Olejnik & Algina, 2003) along with an \emph{optional} confidence interval
#' computed by one of two approximate methods. The default is to return only
#' the point estimate (\code{method = "none"}), because both available CI
#' methods are approximations whose coverage properties have not been broadly
#' validated for this estimand and warrant independent evaluation before being
#' used in substantive inference.
#'
#' @param object Optional. A fitted model object of class
#'   \code{\link[stats]{aov}}, \code{\link[stats]{lm}}, or \code{aovlist}
#'   (multi-stratum aov fit for within-subjects / mixed designs).
#' @param observed Character vector of factor names treated as measured.
#' @param SS_effect,SS_observed,SS_error Sums of squares (option 2 in
#'   \code{\link{eta_squared_generalized}}): \code{SS_effect} for the focal
#'   effect, \code{SS_error} for its error term, and \code{SS_observed} (a
#'   scalar or a vector) for every other source in the denominator of
#'   \eqn{\eta^2_G}. The estimate depends on \code{SS_observed} and
#'   \code{SS_error} only through their sum. The parametric interval is
#'   built on the \emph{F} they imply,
#'   \eqn{(\mathit{SS}_{\text{effect}}/df_{\text{effect}})/(\mathit{SS}_{\text{error}}/df_{\text{error}})},
#'   which is the \emph{F} of the effect's test only when \code{SS_error} is
#'   the sum of squares that test divides by: the residual, adjusted for any
#'   covariate, in a design without an \code{Error()} term, and in a design
#'   with one the residual of the stratum in which the effect is tested
#'   (subjects within groups for an effect tested between subjects, and in a
#'   one-way within-subjects design the interaction of subjects and the
#'   within-subjects factor). \code{SS_observed} then holds the other
#'   sources that Olejnik and Algina (2003, Eq. 5, p. 437) place in the
#'   denominator: the measured factors, the interactions that contain one,
#'   the covariates, and the residuals of the other strata, such as that of
#'   the subjects stratum of a within-subjects design. The interval also
#'   needs \code{df_effect}, \code{df_error}, and \code{N}.
#' @param F_effect,df_effect,F_observed,df_observed,df_error \emph{F}-values
#'   and degrees of freedom (option 3 in
#'   \code{\link{eta_squared_generalized}}). Each \emph{F} divides a mean
#'   square by the focal effect's error mean square: \code{F_effect} is the
#'   effect's test, and \code{F_observed} covers the sources that
#'   \code{SS_observed} holds, so for a source in another stratum it is not
#'   the \emph{F} that \code{summary()} of the fit reports.
#'   \code{df_effect} and \code{df_error} go with the sums of squares as
#'   well. Given with the sums of squares, \code{F_effect},
#'   \code{F_observed}, and \code{df_observed} must describe the same test:
#'   the \emph{F} the sums of squares imply for the effect, and the ratio of
#'   \code{SS_observed} to \code{SS_error}, must agree with the \emph{F}
#'   interface's to a relative tolerance of 1e-6, the tolerance
#'   \code{\link{eta_squared_generalized}} applies to the two estimates, and
#'   the interval is then formed from the \emph{F}-values. The two can agree on the estimate and still imply
#'   different \emph{F} statistics for the effect, and so different
#'   intervals; such a call is refused.
#' @param N Total sample size, with which the confidence limits for the
#'   noncentrality parameter are converted (see Details). Required when
#'   \code{method = "parametric"} and no fitted model is supplied; ignored
#'   when \code{method = "none"} or when a fitted model is supplied, in
#'   which case a single-stratum fit is converted with
#'   \code{\link[stats]{nobs}(object)} and each effect of an \code{aovlist}
#'   fit with the size of its own stratum, as in
#'   \code{\link{ci_eta_squared}}. For an effect tested in an
#'   \code{Error()} stratum, give the size of that stratum instead: the
#'   number of subjects for an effect tested between subjects, and for an
#'   effect tested within subjects the stratum's total degrees of freedom,
#'   which is \code{df_effect + df_error} when the effect is the only one
#'   tested in its stratum. \code{N} must be finite and at least
#'   \code{df_effect + df_error}.
#' @param method One of \code{"none"} (default), \code{"parametric"}, or
#'   \code{"bootstrap"}. See Details.
#' @param B Number of bootstrap replications when \code{method =
#'   "bootstrap"}, a single whole number. Default \code{10000}. The fewest
#'   accepted is 1000 through the 98 percent level (this function's own
#'   floor) and 2000 at the 99 percent level, where
#'   10 replications beyond each limit take that many; a smaller \code{B}
#'   is refused before anything is computed (class
#'   \code{dmar_count_invalid}). A reported interval deserves \code{10000}
#'   or more. See \emph{Bootstrap replications} under \emph{Details} and
#'   \code{\link{dmar_bootstrap}}.
#' @param conf_level Desired confidence coverage; default \code{0.95}.
#' @param alpha_lower,alpha_upper Optional Type I error on the lower and upper
#'   side.
#'   The table records the coverage, \code{1 - alpha_lower - alpha_upper},
#'   as its level and, when the two tail areas differ, the areas
#'   themselves, which the line beneath the printed table,
#'   \code{\link{results_sentence}}, and the \code{alpha_lower} and
#'   \code{alpha_upper} columns of \code{tidy()} state; with a tail area
#'   of 0 the interval is a one-sided confidence bound.
#' @param seed Optional integer seed for the bootstrap, for
#'   reproducibility. Used locally: the caller's random number generator
#'   state is restored on exit. Default \code{NULL} draws the resamples
#'   from the current random number stream, which advances that stream.
#'
#' @return A \code{data.frame} with one row per focal effect and the
#'   columns \code{effect}, \code{eta_squared_generalized},
#'   \code{ci_lower}, \code{ci_upper}, and \code{method}. For
#'   \code{aovlist} fits a \code{stratum} column is also present. When
#'   \code{method = "none"}, the limit columns contain \code{NA}.
#'   Otherwise the table records the confidence level of its intervals in
#'   its \code{conf_level} attribute: \code{conf_level}, or
#'   \code{1 - alpha_lower - alpha_upper} when the tail areas are
#'   supplied. With
#'   \code{method = "bootstrap"} the attributes \code{B}, \code{B_used},
#'   \code{n_dropped}, \code{B_required}, and \code{bootstrap_status}
#'   record the replications (see \emph{Bootstrap replications} under
#'   Details).
#'
#' @details
#' \strong{Why CI = "none" is the default.} Confidence interval construction
#' for \eqn{\eta^2_G} is not as settled as for partial \eqn{\eta^2} or
#' \eqn{\omega^2}, because the denominator mixes sums of squares from
#' heterogeneous sources (the focal effect, one or more measured factors, and
#' the error term). No noncentral \emph{F} transformation maps the population
#' noncentrality parameter directly to \eqn{\eta^2_G}. Both methods below are
#' approximations and are exposed for exploration rather than as defaults.
#'
#' \strong{\code{method = "parametric"}.} The function first obtains a
#' confidence interval for the population noncentrality parameter
#' \eqn{\lambda} of the focal effect's \emph{F}-test via
#' \code{\link{ci_nc_F}} and maps its limits through the partial
#' \eqn{\eta^2} transformation
#' \eqn{\eta^2_{p,\text{bound}} = \lambda_{\text{bound}}/(\lambda_{\text{bound}} + N)},
#' as \code{\link{ci_eta_squared}} converts them: \emph{N} is the number of
#' observations for a fit without an \code{Error()} term (the convention of
#' \code{\link{ci_pvaf}} and \code{\link{ci_omega_squared}} as well), the
#' size of the effect's stratum for an effect of an \code{aovlist} fit (see
#' below), and \code{N} as supplied for the raw interfaces. The limits are
#' then re-expressed as \eqn{\eta^2_G} bounds via
#' \deqn{\eta^2_{G,\text{bound}} = \frac{r_{\text{bound}}}{r_{\text{bound}} + r_{\text{obs}} + 1},}
#' where \eqn{r_{\text{bound}} = \eta^2_{p,\text{bound}}/(1-\eta^2_{p,\text{bound}}) = \lambda_{\text{bound}}/N}
#' and \eqn{r_{\text{obs}}} is the sum of the other sums of squares in the
#' denominator of \eqn{\eta^2_G} (those of the measured sources and, for an
#' \code{aovlist} fit, the residuals of the other strata) divided by the
#' focal effect's error sum of squares. The point estimate is the same map
#' applied to
#' \eqn{\mathit{SS}_{\text{effect}}/\mathit{SS}_{\text{error}} = F \cdot df_{\text{effect}}/df_{\text{error}}},
#' which for a fitted model converges to \eqn{\lambda/N} as the sample
#' grows. When \eqn{r_{\text{obs}}} is 0, as in a fit without an
#' \code{Error()} term in which nothing measured enters the denominator,
#' \eqn{\eta^2_G} is the partial \eqn{\eta^2}, and the interval is the
#' exact one \code{ci_eta_squared(estimand = "partial")} gives. Otherwise the
#' interval is not exact, and no conversion of the limits for
#' \eqn{\lambda} can make it so: \eqn{\eta^2_G} then depends on the value
#' that \eqn{r_{\text{obs}}} converges to as well as on \eqn{\lambda},
#' while the distribution of the focal \emph{F} depends on \eqn{\lambda}
#' alone, so the transformation principle behind the exact intervals
#' (Steiger, 2004, Proposition 1, p. 167), which maps the limits through a
#' function of that one parameter, does not apply. The interval takes
#' \eqn{r_{\text{obs}}} as known and equal to its sample value, which
#' leaves the sampling variability of \eqn{r_{\text{obs}}} out of it. Its
#' coverage has not been validated, and it should be treated as
#' preliminary.
#'
#' \strong{\code{method = "bootstrap"}.} A residual bootstrap from the
#' fitted model: the resampling unit is a residual, drawn
#' nonparametrically from the model's own residuals rather than from a
#' fitted distribution. For each of the \code{B} replications, the
#' response is regenerated as \eqn{\hat{y}_i + \varepsilon^*_i} where
#' \eqn{\varepsilon^*} is sampled with replacement from the model's
#' residuals; the model is refit; \eqn{\eta^2_G} is recomputed; and the
#' percentile interval, the empirical quantiles of the \code{B}
#' bootstrap values (Efron & Tibshirani, 1993), is reported. The
#' percentile interval is the only bootstrap interval offered; there is
#' no bias-corrected and accelerated (BCa) variant. Bootstrap results
#' vary from run to run; supply \code{seed} for reproducibility. Requires a
#' single-stratum \code{aov}/\code{lm} fit. \strong{Not yet supported
#' for \code{aovlist} (multi-stratum / within-subjects) fits}, since
#' the bootstrap needs to respect the subject-level correlation
#' structure, which naive residual resampling does not. Use
#' \code{method = "parametric"} for \code{aovlist} fits.
#' Coverage has not been broadly validated for this estimand.
#'
#' \strong{Bootstrap replications.} The limits are quantiles of the
#' \code{B} replicated values of each effect, read at \code{alpha_lower}
#' and \code{1 - alpha_upper}. The function requires at least 1000
#' replications, more than the 400 that 10 replications beyond each limit
#' take at the 95 percent level; at the 99 percent level that count is 2000
#' and governs. A replication in which \eqn{\eta^2_G} for some effect is
#' undefined (the refit fails, or the effect size cannot be computed from
#' it) is dropped and counted, never replaced, and the interval is computed
#' from the rest. The result records \code{B} (requested), \code{B_used}
#' (the number the interval rests on), \code{n_dropped}, and
#' \code{B_required} (the fewest replications a request may ask for; that
#' count is held to \code{B}, not to \code{B_used}) as attributes, with
#' \code{bootstrap_status}
#' (\code{"complete"}, \code{"dropped"}, \code{"exploratory"},
#' \code{"insufficient_count"}, or \code{"insufficient_share"}). One
#' warning of class \code{dmar_bootstrap_replications_dropped} states the
#' counts when any replication was dropped, with a caution when the dropped
#' replications are at least the 2.5 percent of \code{B} that lies beyond
#' one limit at the 95 percent level, because they could then decide by
#' themselves where that limit falls. When more than half were dropped, the
#' limits of every effect are
#' \code{NA}, the estimate is returned, and the warning has class
#' \code{dmar_bootstrap_insufficient}. With several effects, \code{B_used}
#' is the fewest replications the limits of any effect rest on.
#' \code{options(dmar.bootstrap_exploratory = TRUE)} lets a smaller
#' \code{B} (2 or more) run for a first look or a demonstration; every such
#' call gives a warning of class \code{dmar_bootstrap_exploratory} and the
#' result is not one to report. See \code{\link{dmar_bootstrap}} for the
#' rule and its reasons.
#'
#' \strong{Within-subjects designs (\code{aovlist}).} Each focal effect's
#' parametric interval is built from its own stratum's \emph{F} test, and
#' its limits for \eqn{\lambda} are converted with the size of that
#' stratum, as \code{\link{ci_eta_squared}} converts them: the stratum's
#' total degrees of freedom for an effect tested within subjects, and the
#' number of subjects for an effect in the between-subjects stratum of a
#' mixed design. \eqn{r_{\text{obs}}} includes the residual sums of squares
#' of the other strata, since every source of variance that involves
#' subjects belongs to the denominator of \eqn{\eta^2_G} (Olejnik & Algina,
#' 2003, Eq. 5, p. 437; Bakeman, 2005, p. 381). With \eqn{r_{\text{obs}}}
#' at the value it converges to, the conversion maps each limit for
#' \eqn{\lambda} to the corresponding value of the population
#' \eqn{\eta^2_G} of Olejnik and Algina (2003, Eq. 2, p. 436), the value the
#' estimate converges to as subjects are added. In a one-way
#' within-subjects design with \emph{n} subjects, \emph{p} occasions, and
#' compound symmetry, for which Steiger (2004, Eq. 58, p. 177) gives
#' \eqn{\lambda = n \sum_j (\mu_j - \mu)^2 / [\sigma^2 (1 - \rho)]}, that
#' value is \eqn{\lambda (1 - \rho)/[\lambda (1 - \rho) + np]}: the
#' variance of the occasion means, \eqn{\sum_j (\mu_j - \mu)^2 / p}, over
#' the sum of that variance and \eqn{\sigma^2}, which holds the subject,
#' subject by occasion, and error variances that Olejnik and Algina (2003,
#' p. 443, note 1) place in the denominator. Given \eqn{\lambda}, that
#' value still depends on \eqn{\rho}, and the distribution of the
#' occasions' \emph{F}, which depends on \eqn{\lambda} alone, carries no
#' further information about \eqn{\rho}. Steiger (2004, pp. 176--177)
#' notes that in repeated measures designs \eqn{\lambda} confounds the
#' effect with the correlation among the observations, and that the
#' root-mean-square standardized effect of such a design cannot be
#' estimated directly with the exact noncentral methods unless \eqn{\rho}
#' is known. Converted with the number of rows, \eqn{np}, the limits
#' would be for a smaller quantity than the one the estimate estimates:
#' for the three IQ occasions of \code{\link{pygmalion}} analyzed in
#' \code{vignette("DMAR", package = "DMAR")}, 930 rows against
#' \emph{F}(2, 618) give the interval [.019, .037], below the estimate of
#' .041; converted with 2 + 618 = 620, the interval is [.029, .054]. As
#' \code{\link{eta_squared_generalized}} does, the function stops when an
#' effect is tested in more than one stratum, as happens when observations
#' are missing within subjects or a covariate varies within them, since no
#' one stratum's sum of squares is then the effect's.
#'
#' @references
#' Algina, J., Keselman, H. J., & Penfield, R. D. (2005). An alternative to
#'   Cohen's standardized mean difference effect size: A robust parameter and
#'   confidence interval in the two independent groups case.
#'   \emph{Psychological Methods, 10}(3), 317--328.
#'   \doi{10.1037/1082-989X.10.3.317}
#'
#' Bakeman, R. (2005). Recommended effect size statistics for repeated
#'   measures designs. \emph{Behavior Research Methods, 37}(3), 379--384.
#'   \doi{10.3758/BF03192707}
#'
#' Efron, B., & Tibshirani, R. J. (1993). \emph{An introduction to the
#'   bootstrap}. New York, NY: Chapman & Hall/CRC.
#'
#' Kelley, K. (2007). Confidence intervals for standardized effect
#'   sizes: Theory, application, and implementation. \emph{Journal of Statistical
#'   Software, 20}(8), 1--24. \doi{10.18637/jss.v020.i08}
#'
#' Kelley, K., & Preacher, K. J. (2012). On effect size.
#'   \emph{Psychological Methods, 17}, 137--152. \doi{10.1037/a0028086}
#'
#' Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). \emph{Designing
#'   experiments and analyzing data: A model comparison perspective}
#'   (4th ed.). Routledge. (See Chapter 3 on \eqn{\eta^2}, Chapter 7 on
#'   factorial designs, and Chapter 11 on generalized \eqn{\eta^2} for
#'   within-subjects designs.)
#'
#' Olejnik, S., & Algina, J. (2003). Generalized eta and omega squared
#'   statistics: Measures of effect size for some common research designs.
#'   \emph{Psychological Methods, 8}(4), 434--447.
#'   \doi{10.1037/1082-989X.8.4.434}
#'
#' Smithson, M. (2001). Correct confidence intervals for various regression
#'   effect sizes and parameters: The importance of noncentral distributions in
#'   computing intervals. \emph{Educational and Psychological Measurement, 61},
#'   605--632. \doi{10.1177/00131640121971392}
#'
#' Steiger, J. H. (2004). Beyond the \emph{F} test: Effect size confidence
#'   intervals and tests of close fit in the analysis of variance and contrast
#'   analysis. \emph{Psychological Methods, 9}(2), 164--182.
#'   \doi{10.1037/1082-989X.9.2.164}
#'
#' @seealso \code{\link{eta_squared_generalized}}, \code{\link{ci_eta_squared}}
#'
#' @examples
#' # The pygmalion expectancy experiment: treatment is manipulated, while
#' # grade is a measured classification, so grade belongs in the denominator.
#' pyg <- pygmalion
#' pyg$grade <- factor(pyg$grade)
#' fit <- aov(iq_8 ~ treatment * grade, data = pyg)
#'
#' # The default returns the point estimate and leaves the limits NA, because
#' # both interval methods are approximations.
#' ci_eta_squared_generalized(fit, observed = "grade")
#'
#' # The parametric approximation maps a noncentral F interval for the focal
#' # effect through the observed sums of squares. Warnings mark the method as
#' # preliminary and report that the interaction's lower limit is clamped at
#' # 0; treat the limits accordingly.
#' ci_eta_squared_generalized(fit, observed = "grade", method = "parametric")
#'
#' # The third option is a residual bootstrap, which refits the model once per
#' # replication. B = 100 keeps the example quick and is below the 1000 the
#' # function requires; the option lets it run and the warning says so. A
#' # reported interval deserves B = 10000 or more. The seed makes the
#' # interval reproducible and leaves the random number stream of the
#' # surrounding session as it was.
#' old <- options(dmar.bootstrap_exploratory = TRUE)
#' ci_eta_squared_generalized(fit, observed = "grade",
#'                            method = "bootstrap", B = 100, seed = 113)
#' options(old)
#'
#' # The parametric treatment row from summary statistics, as an ANOVA table
#' # reports them: the sums of squares of treatment, of its error term (the
#' # residual), and of the measured sources in its denominator (grade and
#' # the treatment by grade interaction), with the degrees of freedom and
#' # the number of children for the interval.
#' tab <- anova(fit)
#' ci_eta_squared_generalized(
#'   SS_effect   = tab["treatment", "Sum Sq"],
#'   SS_observed = tab[c("grade", "treatment:grade"), "Sum Sq"],
#'   SS_error    = tab["Residuals", "Sum Sq"],
#'   df_effect   = tab["treatment", "Df"], df_error = tab["Residuals", "Df"],
#'   N = nobs(fit), method = "parametric")
#'
#' # Within-subjects ANOVA. The parametric CI uses each effect's own stratum
#' # and converts its limits with that stratum's size.
#' set.seed(113)
#' n <- 20
#' rm_data <- data.frame(
#'   subject = factor(rep(seq_len(n), each = 3)),
#'   time    = factor(rep(c("Pre", "Mid", "Post"), n),
#'                    levels = c("Pre", "Mid", "Post")),
#'   y       = rnorm(n, sd = 1.5)[rep(seq_len(n), each = 3)] +
#'             0.7 * rep(1:3, n) + rnorm(n * 3, sd = 1.2)
#' )
#' fit_rm <- aov(y ~ time + Error(subject/time), data = rm_data)
#' ci_eta_squared_generalized(fit_rm, method = "parametric")
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @keywords htest design
#'
#' @family confidence intervals for effect sizes
#'
#' @export

ci_eta_squared_generalized <- function(
  object       = NULL,
  observed     = NULL,
  SS_effect    = NULL,
  SS_observed  = NULL,
  SS_error     = NULL,
  F_effect     = NULL,
  df_effect    = NULL,
  F_observed   = NULL,
  df_observed  = NULL,
  df_error     = NULL,
  N            = NULL,
  method       = c("none", "parametric", "bootstrap"),
  B            = 10000L,
  conf_level   = 0.95,
  alpha_lower  = NULL,
  alpha_upper  = NULL,
  seed         = NULL
) {
  .dmar_check_conf_level(conf_level, null_ok = TRUE)
  method <- match.arg(method)

  # df_effect and df_error belong to both summary-statistics interfaces:
  # the F interface needs them for the estimate, and the sums of squares
  # need them, with N, for the parametric interval only. F_effect,
  # F_observed, and df_observed belong to the F interface alone. Given both
  # interfaces, the estimate could agree while the F of the effect, and so
  # the interval, did not. Such a call is answered, as DMAR 1.0.x answered
  # it, when the two describe the same test: the F that the sums of squares
  # imply for the effect, and the ratio of the measured sources' sums of
  # squares to the error sum of squares, agree with the F interface's to
  # a relative tolerance of 1e-6, the tolerance eta_squared_generalized()
  # applies to the two estimates; the interval is then formed from the
  # F-values, as below. When they disagree, which interval is meant cannot be told, and
  # the call is refused.
  ss_given <- !is.null(SS_effect) || !is.null(SS_observed) || !is.null(SS_error)
  F_given  <- !is.null(F_effect)  || !is.null(F_observed)  || !is.null(df_observed)
  if (is.null(object) && ss_given && F_given) {
    .eta_g_both_interfaces(SS_effect, SS_observed, SS_error, F_effect,
                           df_effect, F_observed, df_observed, df_error)
  }

  # eta_squared_generalized() takes df_effect and df_error beside the sums
  # of squares and leaves the estimate unchanged; DMAR 1.0.x read them as
  # the F interface there and stopped for want of F_effect.
  point <- eta_squared_generalized(
    object      = object,      observed    = observed,
    SS_effect   = SS_effect,   SS_observed = SS_observed, SS_error = SS_error,
    F_effect    = F_effect,    df_effect   = df_effect,
    F_observed  = F_observed,  df_observed = df_observed, df_error = df_error
  )
  point$ci_lower <- NA_real_
  point$ci_upper <- NA_real_
  point$method      <- method

  if (method == "none") {
    message(
      "No CI computed (method = 'none'). Set method = 'parametric' for an ",
      "approximate noncentral F transformation, or method = 'bootstrap' for a ",
      "residual bootstrap (requires a fitted model). Both methods are ",
      "approximate; see ?ci_eta_squared_generalized."
    )
    # Wide-format dmar_tbl: a leading `effect` label column beside the
    # eta_squared_generalized estimate and its (here NA) interval limits. No
    # conf_level is attached because method = "none" computes no interval, so a
    # confidence-level footer would be misleading.
    return(.as_dmar_tbl(point, subclass = "dmar_ci_anova"))
  }

  alphas <- .ci_eta_squared_alphas(conf_level, alpha_lower, alpha_upper)

  if (method == "parametric") {
    .dmar_warn(c(
      "Parametric CI uses an approximate transformation that maps the partial ",
      "eta squared CI through the observed sums of squares for measured factors. ",
      "Coverage properties have not been broadly validated; results should be ",
      "treated as preliminary and independently evaluated."),
      "result_qualified")
    # Wide-format dmar_tbl with the CI present; the conf_level footer applies.
    # One clamp warning per call, not one per clamped effect.
    return(.as_dmar_tbl(
      .warn_ncf_clamp_once(.ci_eta_g_parametric(
        point, object, observed,
        SS_effect, SS_observed, SS_error,
        F_effect, df_effect, F_observed, df_observed, df_error, N,
        alphas$alpha_lower, alphas$alpha_upper
      )),
      conf_level = alphas$conf_level,
      tail_areas = .tail_areas(alphas$alpha_lower, alphas$alpha_upper)
    , subclass = "dmar_ci_anova"))
  }

  if (method == "bootstrap") {
    if (is.null(object)) {
      .dmar_refuse_argument(
        "method = 'bootstrap' requires a fitted model ('object'); the raw SS / F-df interfaces cannot bootstrap.",
        arg = "object")
    }
    if (inherits(object, "aovlist")) {
      .dmar_refuse_argument(
        "method = 'bootstrap' is not yet supported for aovlist (multi-stratum / within-subjects) fits. Use method = 'parametric' for an approximate CI on aovlist models.",
        arg = NULL)
    }
    # The rule of ?dmar_bootstrap at the smaller tail area asked for, with
    # this function's own floor of 1000 replications, which it has
    # documented since 1.0.0 and which the rule does not lower.
    boot_rule <- .ci_eta_g_boot_rule(alphas$alpha_lower, alphas$alpha_upper,
                                     alphas$conf_level)
    B <- .dmar_check_B(B, tail = boot_rule$tail,
                       conf_level = boot_rule$conf_level, floor = 1000L,
                       what = boot_rule$what, position = boot_rule$position)
    .dmar_warn(
      sprintf(
"Bootstrap CI resamples the fitted model's residuals nonparametrically (B = %d). Coverage has not been broadly validated for generalized eta squared; results should be treated as preliminary.",
as.integer(B)
), "result_qualified")
    # Wide-format dmar_tbl with the bootstrap CI present; conf_level footer applies.
    boot <- .ci_eta_g_bootstrap(
      object, observed, as.integer(B),
      alphas$alpha_lower, alphas$alpha_upper, seed, boot_rule
    )
    out <- .as_dmar_tbl(
      boot$table,
      conf_level = alphas$conf_level,
      tail_areas = .tail_areas(alphas$alpha_lower, alphas$alpha_upper)
    , subclass = "dmar_ci_anova")
    return(.dmar_boot_attrs(out, boot$verdict))
  }
}


# Parametric approximate CI: map NCP bound -> partial eta_sq -> generalized.
.ci_eta_g_parametric <- function(point, object, observed,
                                 SS_effect, SS_observed, SS_error,
                                 F_effect, df_effect,
                                 F_observed, df_observed, df_error, N,
                                 alpha_lower, alpha_upper) {

  per_effect <- function(eff_name, F_v, df_e, df_err, ss_obs_ratio, N_use) {
    ncp <- .ci_nc_F_for(
      caller      = "ci_eta_squared_generalized",
      quantity    = "generalized eta squared",
      F_value     = F_v,     df_1 = df_e, df_2 = df_err,
      alpha_lower = alpha_lower, alpha_upper = alpha_upper,
      conf_level  = NULL
    )
    lo_ncp <- ncp$value[ncp$term == "ci_lower"]
    up_ncp <- ncp$value[ncp$term == "ci_upper"]

    lo_etap <- lo_ncp / (lo_ncp + N_use)
    # A NA upper NCP means the observed F is so small that even at lambda = 0 the
    # lower-tail probability is already at or below alpha_upper, leaving the upper
    # noncentrality limit undefined (ci_nc_F documents this). Propagate NA
    # to the upper generalized eta squared limit, matching ci_eta_squared.
    up_etap <- if (is.na(up_ncp)) NA_real_
               else if (is.infinite(up_ncp)) 1
               else up_ncp / (up_ncp + N_use)

    lo_r <- if (lo_etap <= 0) 0 else lo_etap / (1 - lo_etap)
    up_r <- if (is.na(up_etap)) NA_real_
            else if (up_etap >= 1) Inf
            else up_etap / (1 - up_etap)

    lo_eta_g <- lo_r / (lo_r + ss_obs_ratio + 1)
    up_eta_g <- if (is.na(up_r)) NA_real_
                else if (is.infinite(up_r)) 1
                else up_r / (up_r + ss_obs_ratio + 1)

    list(lower = lo_eta_g, upper = up_eta_g)
  }

  if (!is.null(object)) {
    if (inherits(object, "aovlist")) {
      effects_tbl <- .aovlist_effects_table(object)
      total_ss_err <- .aovlist_total_residual_ss(object)
      effect_names <- effects_tbl$effect
      ss_lookup    <- setNames(effects_tbl$SS_effect, effect_names)
      observed_in  <- .expand_observed(observed, effect_names)
      # Each effect's noncentrality limits are converted with the size of
      # its own stratum, as ci_eta_squared() converts them, so that
      # lambda / size is the value the effect's SS / SS_error converges to,
      # the ratio the point estimate is built on. DMAR 1.0.x converted them
      # with the number of rows, which within subjects put the interval on
      # a smaller quantity than the one the estimate estimates.
      sizes        <- .aovlist_stratum_sizes(object)

      for (i in seq_len(nrow(point))) {
        eff   <- point$effect[i]
        row   <- effects_tbl[effects_tbl$effect == eff &
                               effects_tbl$stratum == point$stratum[i], ,
                             drop = FALSE]
        F_v   <- row$F_value
        df_e  <- row$df_effect
        df_err <- row$df_error
        ss_err_focal <- row$SS_error

        obs_others <- setdiff(observed_in, eff)
        ss_obs <- if (length(obs_others) == 0L) 0 else sum(ss_lookup[obs_others])
        # Per-effect denominator ratio relative to the focal stratum's
        # error, not to the error SS of all strata. The per_effect helper
        # appends `+ 1` for the focal stratum's residual, so subtract that
        # here.
        ss_obs_ratio <- (ss_obs + total_ss_err - ss_err_focal) / ss_err_focal

        lims <- per_effect(eff, F_v, df_e, df_err, ss_obs_ratio,
                           sizes[[row$stratum]])
        point$ci_lower[i] <- lims$lower
        point$ci_upper[i] <- lims$upper
      }
      return(point)
    }

    tbl    <- stats::anova(object)
    ss_err <- tbl["Residuals", "Sum Sq"]
    df_err <- tbl["Residuals", "Df"]
    N_use  <- stats::nobs(object)
    observed_in <- .expand_observed(observed, setdiff(rownames(tbl), "Residuals"))

    for (i in seq_len(nrow(point))) {
      eff <- point$effect[i]
      F_v <- tbl[eff, "F value"]
      df_e <- tbl[eff, "Df"]
      obs_others <- setdiff(observed_in, eff)
      ss_obs <- if (length(obs_others) == 0L) 0 else sum(tbl[obs_others, "Sum Sq"])
      lims <- per_effect(eff, F_v, df_e, df_err, ss_obs / ss_err, N_use)
      point$ci_lower[i] <- lims$lower
      point$ci_upper[i] <- lims$upper
    }
    return(point)
  }

  # The summary-statistics interfaces describe one effect: the F of its
  # test, given or implied by the sums of squares, and the ratio of the
  # other sums of squares in the denominator to the effect's error sum of
  # squares, which the F interface gives as F_observed * df_observed /
  # df_error. eta_squared_generalized() has validated the arguments of
  # the estimate, and the F interface's degrees of freedom among them.
  if (!is.null(F_effect)) {
    if (is.null(N)) {
      .dmar_refuse_argument(
        "Parametric CI from the F/df interface requires 'N' (total sample size).",
        arg = "N")
    }
    F_v <- F_effect
    ss_obs_ratio <- sum(F_observed * df_observed / df_error)
  } else {
    if (is.null(df_effect) || is.null(df_error) || is.null(N)) {
      .dmar_refuse_argument(c(
        "Parametric CI from the SS interface additionally requires 'df_effect', ",
        "'df_error', and 'N' so that the noncentral F transformation can be applied."),
        arg = c("df_effect", "df_error", "N"))
    }
    .eta_g_check_df(df_effect, df_error)
    # SS_error is the sum of squares of the effect's error term, as the
    # help defines it, so this is the F of the effect's test, as anova()
    # and summary() of a fit compute it.
    F_v <- (SS_effect / df_effect) / (SS_error / df_error)
    ss_obs_ratio <- sum(SS_observed) / SS_error
  }
  # A within-subjects stratum that tests one effect has exactly
  # df_effect + df_error degrees of freedom, its size, so N may equal
  # that sum, as in ci_eta_squared(); a smaller N describes no design, and
  # an infinite one converts every limit to 0.
  if (!is.numeric(N) || length(N) != 1L || !is.finite(N) ||
      N < df_effect + df_error)
    .dmar_refuse_argument(
      "'N' must be a single finite number of at least df_effect + df_error.",
      arg = "N")
  lims <- per_effect("overall", F_v, df_effect, df_error, ss_obs_ratio, N)
  point$ci_lower[1] <- lims$lower
  point$ci_upper[1] <- lims$upper
  point
}


# Residual bootstrap CI from a fitted model.
# What the rule of R/dmar_bootstrap.R needs to know about a percentile
# interval with the tail areas alpha_lower and alpha_upper: the smallest
# tail area a requested limit cuts off (min(a, 1 - a) over the positive
# areas) and how the result is named in a message.
.ci_eta_g_boot_rule <- function(alpha_lower, alpha_upper, conf_level) {
  areas <- c(alpha_lower, alpha_upper)
  asked <- areas > 0
  what <- if (!all(asked)) "bound" else
    if (isTRUE(all.equal(alpha_lower, alpha_upper))) "interval" else "limits"
  list(tail = min(pmin(areas[asked], 1 - areas[asked])), what = what,
       position = if (what != "bound") NULL else
         if (alpha_lower > 0) alpha_lower else 1 - alpha_upper,
       conf_level = conf_level)
}

# The residual bootstrap. Returns list(table, verdict): the table of
# point estimates with the percentile limits, and the verdict of the rule
# on the replications (.dmar_boot_settle()). A replication whose refit or
# whose effect size computation fails has no value and is dropped and
# counted; B_used is the fewest replications any effect's limits rest on.
# Through DMAR 1.0.x such replications were removed without a count or a
# warning.
.ci_eta_g_bootstrap <- function(object, observed, B,
                                alpha_lower, alpha_upper, seed, rule) {
  .dmar_local_seed(seed)
  fitted_vals <- stats::fitted(object)
  resid_vals  <- stats::residuals(object)
  if (length(fitted_vals) != length(resid_vals)) {
    .dmar_refuse_argument(
      "Fitted values and residuals have inconsistent lengths.", arg = NULL)
  }
  model_data <- stats::model.frame(object)
  resp_name  <- names(model_data)[1L]

  base_point <- eta_squared_generalized(object = object, observed = observed)
  effect_names <- base_point$effect
  K <- length(effect_names)

  boot_mat <- matrix(NA_real_, nrow = B, ncol = K)
  colnames(boot_mat) <- effect_names

  for (b in seq_len(B)) {
    eps <- sample(resid_vals, replace = TRUE)
    y_boot <- fitted_vals + eps
    boot_data <- model_data
    boot_data[[resp_name]] <- y_boot
    fit_b <- try(stats::update(object, data = boot_data), silent = TRUE)
    if (inherits(fit_b, "try-error")) next
    eg_b <- try(eta_squared_generalized(object = fit_b, observed = observed),
                silent = TRUE)
    if (inherits(eg_b, "try-error")) next
    # Align by effect name
    idx <- match(effect_names, eg_b$effect)
    boot_mat[b, ] <- eg_b$eta_squared_generalized[idx]
  }

  verdict <- .dmar_boot_settle(
    B, min(colSums(is.finite(boot_mat))), tail = rule$tail,
    conf_level = rule$conf_level, fun = "ci_eta_squared_generalized",
    quantity = "generalized eta squared for some effect",
    dropped_because = paste0("gave a refit that failed or an effect size ",
                             "that could not be computed"),
    floor = 1000L, what = rule$what, position = rule$position)
  point <- base_point
  # The columns keep their order: the lower limit, then the upper.
  point$ci_lower <- NA_real_
  point$ci_upper <- NA_real_
  if (verdict$formed) {
    point$ci_lower <- apply(boot_mat, 2L, function(x) {
      stats::quantile(x[is.finite(x)], probs = alpha_lower, names = FALSE)
    })
    point$ci_upper <- apply(boot_mat, 2L, function(x) {
      stats::quantile(x[is.finite(x)], probs = 1 - alpha_upper,
                      names = FALSE)
    })
  }
  point$method <- "bootstrap"
  rownames(point) <- NULL
  list(table = point, verdict = verdict)
}


# .eta_g_both_interfaces(SS_effect, SS_observed, SS_error, F_effect,
#                        df_effect, F_observed, df_observed, df_error)
#
# Given the sums of squares and the F interface together (and no fitted
# model), checks that they describe the same test of the effect, and stops
# otherwise: the F of the effect that the sums of squares imply,
# (SS_effect / df_effect) / (SS_error / df_error), against F_effect, and
# the ratio sum(SS_observed) / SS_error against
# sum(F_observed * df_observed) / df_error, each to a relative tolerance of
# 1e-6, as eta_squared_generalized() compares the two estimates. Returns
# NULL invisibly when they agree.
.eta_g_both_interfaces <- function(SS_effect, SS_observed, SS_error, F_effect,
                                   df_effect, F_observed, df_observed,
                                   df_error) {
  if (is.null(SS_effect) || is.null(SS_error) || is.null(F_effect) ||
      is.null(df_effect) || is.null(df_error)) {
    .dmar_refuse_argument(c(
      "Supply the sums of squares (SS_effect, SS_observed, SS_error) or ",
      "the F interface (F_effect, F_observed, df_observed) in full, with ",
      "df_effect and df_error; given parts of both, the test of the effect ",
      "cannot be told. df_effect, df_error, and N go with either one."),
      arg = NULL)
  }
  .eta_g_check_df(df_effect, df_error)
  F_ss <- (SS_effect / df_effect) / (SS_error / df_error)
  r_ss <- if (is.null(SS_observed)) 0 else sum(SS_observed) / SS_error
  r_F <- if (is.null(F_observed)) 0 else sum(F_observed * df_observed) / df_error
  agree <- function(a, b) isTRUE(all.equal(a, b, tolerance = 1e-6))
  if (!agree(F_ss, F_effect) || !agree(r_ss, r_F)) {
    .dmar_refuse_argument(c(
      "The sums of squares and the F interface describe different tests of ",
      "the effect: the sums of squares imply F = ", signif(F_ss, 6),
      " for the effect and a ratio of ", signif(r_ss, 6), " for SS_observed ",
      "to SS_error, against F_effect = ", signif(F_effect, 6), " and ",
      signif(r_F, 6), " from F_observed and df_observed. They can agree on ",
      "the estimate and still imply different intervals; supply one of them."),
      arg = NULL)
  }
  invisible(NULL)
}
