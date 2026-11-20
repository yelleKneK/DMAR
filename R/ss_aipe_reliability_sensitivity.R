# Sensitivity analysis for AIPE on a reliability coefficient (parallel tests).
#' Sensitivity Analysis for Sample Size Planning From the AIPE Perspective for a Reliability Coefficient
#'
#' @description
#' Quantifies how much misspecification of the population reliability
#' coefficient distorts an AIPE-based sample size plan for the
#' reliability of a composite score. On each replication the function
#' simulates an \emph{n} \eqn{\times} \emph{J} data matrix of subjects
#' by items from a single-factor model of parallel tests whose
#' population reliability of the sum score equals
#' \code{true_reliability}, fits
#' the requested estimator (alpha or omega) via the corresponding
#' \code{reliability_*} function with the supplied \code{ci_method},
#' and records the realized reliability estimate and its confidence
#' interval.
#'
#' \strong{Population model.} Each item has a single common-factor
#' loading and uncorrelated unique error. With per-item variance
#' normalized to 1, the loading and unique variance are chosen so the
#' reliability of the sum score equals \code{true_reliability}:
#' \deqn{\lambda^2 \;=\; \frac{\rho}{J(1 - \rho) + \rho}, \qquad
#'        \psi^2 \;=\; 1 - \lambda^2,}
#' where \eqn{\rho = }\code{true_reliability} and \eqn{J} is the number of
#' items, \code{n_items}. Item scores are
#' \eqn{y_{ij} = \lambda T_i + e_{ij}}, with \eqn{T_i \sim N(0, 1)} and
#' \eqn{e_{ij} \sim N(0, \psi^2)}.
#'
#' \strong{Missing data.} Two arguments carry a plan into incomplete
#' data, and they stand to each other as \code{supposed_reliability}
#' stands to \code{true_reliability}. \code{fmi} is the planning value: the
#' fraction of missing information the plan assumes for the coefficient,
#' passed to \code{\link{ss_aipe_reliability}}, which enlarges the sample
#' size accordingly. \code{missing_prop} is the data generating mechanism:
#' each simulated item score is deleted independently with that probability
#' (missing completely at random), and every replication is analyzed by
#' full information maximum likelihood. A plan made with \code{fmi} is
#' then judged by the widths the incomplete data actually deliver, and a
#' plan that ignored the missing data (\code{fmi = 0}) is seen to fall
#' short. Planning with \code{fmi} while simulating complete data shows
#' the opposite misspecification, widths narrower than planned.
#'
#' @details
#' \strong{The interval evaluated.} An AIPE plan is made for a particular
#' confidence interval procedure, so a check of the plan evaluates that
#' procedure. The plan here is the one \code{\link{ss_aipe_reliability}}
#' makes for the parallel model with the normal theory interval: its
#' closed form takes the sampling variance of alpha to be the compound
#' symmetry form \eqn{2 (1 - \rho)^2 J / ((J - 1)(N - 1))}{2 (1 - rho)^2 J
#' / ((J - 1)(N - 1))} of van Zyl, Neudecker, and Nel (2000), and its
#' assurance search forms, on each simulated sample, the Wald interval
#' with that variance at the sample alpha. The default \code{ci_method =
#' "ml"} evaluates the interval of \code{\link{reliability_alpha}(ci_method
#' = "ml")}, the Wald interval with the general variance of van Zyl et al.
#' at the sample covariance matrix, which at a population of parallel items
#' is the same variance; for omega it evaluates the interval of
#' \code{\link{reliability_omega}(denominator = "model_implied", ci_method
#' = "ml")}, the delta method Wald interval of the model implied omega,
#' whose variance at a population of parallel items is that same compound
#' symmetry form (it is the variance \code{ss_aipe_reliability()} plans
#' with for the congeneric model, and the two agree to machine precision
#' there). The plan does not describe other intervals, among them
#' Bonett's (2002) for alpha and \code{"mlr"} for omega. With four items, a true
#' and planned reliability of .80, a width of .15, and an assurance of
#' .80, the plan of 99 cases gave alpha intervals no wider than .15 in
#' 81.3 percent of 2,000 replications with \code{"ml"} and in 76.75
#' percent with Bonett's interval (\code{set.seed(113)} before each call);
#' for omega the two intervals came closer (83.3 percent of 1,000 with
#' \code{"ml"} and 81.9 with \code{"mlr"}).
#' Any other interval the estimator offers can still be named in
#' \code{ci_method}, which asks how the planned \emph{N} serves an
#' interval other than the one it was planned for.
#'
#' \strong{The analysis under missing data.} With \code{missing_prop}
#' above 0, coefficient alpha is \code{\link{reliability_alpha}(missing =
#' "fiml")}, the classical coefficient computed from the full information
#' maximum likelihood estimate of the item covariance matrix, whose
#' \code{"ml"} interval, the default here as with complete data, takes its
#' standard error from the FIML information matrix of the saturated model.
#' Coefficient omega is \code{\link{reliability_omega}(missing = "fiml")},
#' whose \code{"ml"} interval takes its standard error from the FIML
#' information matrix of the single-factor model. These are the Wald
#' intervals the planner's \code{fmi} adjustment plans for, with the
#' standard error the missing data inflate. The complete data closed forms (Bonett, Feldt,
#' Fisher, Hakstian-Whalen) and the asymptotic distribution-free intervals
#' do not apply to partially observed rows, so, with \code{missing_prop}
#' above 0, \code{ci_method} must be \code{"ml"} or \code{"ml_logistic"} for
#' alpha and one of \code{"ml"}, \code{"mlr"}, \code{"ml_logistic"}, or
#' \code{"mlr_logistic"} for omega. A row that loses every item score
#' carries no information and is dropped, as lavaan drops it.
#'
#' \strong{Choosing \code{fmi} for a given \code{missing_prop}.} The two
#' are not the same number. The FMI of the coefficient depends on the
#' number of items and how strongly they are correlated, because
#' information lost on one item is partly recovered from the others; under
#' data missing completely at random it is bounded by the largest rates of
#' univariate and pairwise missingness (Savalei & Rhemtulla, 2012). For the
#' population of the second example (four parallel items, reliability .80,
#' width .15), \code{missing_prop = .2} gives coefficient alpha and
#' coefficient omega an FMI of .184 relative to every case, the footing of
#' the planner's \code{fmi} (from one FIML fit to 200,000 simulated cases;
#' it is .182 relative to the cases the fit keeps, since a share
#' \eqn{.2^4 = .0016} of them loses every item). In 1,000 replications
#' with the default \code{ci_method = "ml"} (\code{set.seed(113)} before
#' each call), the complete data plan (\eqn{N = 74}) gave alpha intervals
#' that averaged .153 wide with complete data and .170 with
#' \code{missing_prop = .2}, wider by close to the width inflation factor
#' \eqn{1 / \sqrt{1 - .184} = 1.107}{1 / sqrt(1 - .184) = 1.107}, while the
#' plan made with \code{fmi = .184} (\eqn{N = 91}) averaged .152 with
#' \code{missing_prop = .2}. The help
#' page of \code{\link{ss_aipe_reliability}} describes where a planning FMI
#' comes from in practice.
#'
#' \strong{The fraction of missing information the data realized.} With
#' \code{missing_prop} above 0, every replication that returns an interval
#' also estimates the FMI of its coefficient from the same simulated data,
#' as \code{\link{fraction_missing_information}} does with its default
#' unstructured estimator: the FIML fit behind the coefficient is refit at
#' its estimates with the information of the observed data and with the
#' information complete data would have carried, and the FMI is one minus
#' the ratio of the coefficient's complete data sampling variance to its
#' sampling variance with the data as observed. For alpha
#' the fit is the saturated model of the items whose covariance matrix the
#' classical coefficient is computed from, and the three estimators of
#' Chen and Savalei (2021) coincide for a saturated model (p. 9); for omega
#' it is the single-factor model, and the FMI is that of its model implied
#' coefficient. The FMI is taken relative to all \emph{n} simulated cases,
#' the \emph{N} of the plan. A case that loses every item score is set
#' aside by the fit, and the FMI relative to the \eqn{n_{used}}{n_used}
#' cases the fit used, \eqn{\gamma_{used}}{fmi_used} (the FMI that
#' \code{\link{fraction_missing_information}} reports), leaves out the
#' information those cases would have carried; it is folded in as
#' \eqn{1 - (n_{used} / n)(1 - \gamma_{used})}{1 - (n_used / n)(1 -
#' fmi_used)}. The difference matters with few items and much missing
#' data: a share \eqn{p^J}{p^J} of the cases loses every item, with
#' \emph{p} the value of \code{missing_prop}, which is .125 for three
#' items at \code{missing_prop = .5}. The mean over the replications is the
#' \code{mean_realized_fmi} attribute and is printed under the table, so the
#' \code{fmi} a plan assumed can be read beside the FMI that the simulated
#' missing data cost the coefficient. A single sample FMI is a poor guide:
#' in the simulations of Chen and Savalei (2021, p. 9) the sample FMI of a
#' regression coefficient was largely unbiased but varied from sample to
#' sample, with root mean squared errors of about .10 or more at \emph{N} =
#' 50 and 100.
#' The mean over many replications is the stable figure, and it is the
#' value to carry into \code{fmi} for a study whose missing data the
#' simulation mimics. A replication whose FMI cannot be computed (a fit that
#' fails, or an improper solution that \pkg{lavaan} will not refit at its
#' estimates; see \code{\link{fraction_missing_information}}) is left out of
#' the mean, and the \code{n_realized_fmi} attribute counts the replications
#' that entered it. The computation adds a model fit and two refits to every
#' replication, about doubling the time a replication takes.
#'
#' \strong{Bootstrap intervals.} When \code{ci_method} is a bootstrap
#' method (\code{"bootstrap_se"}, \code{"bootstrap_se_logistic"},
#' \code{"percentile"}, or \code{"bca"}), every one of the \code{G}
#' replications runs its own bootstrap of \code{B} resamples, so the
#' simulation computes the coefficient \eqn{G \times B}{G x B} times. The
#' default \code{B} is the larger of 200 and one and a quarter times the
#' number of replications the interval method and the level require (see
#' \code{\link{dmar_bootstrap}}): at the 95 percent level, 500 for
#' \code{"percentile"}, 2500 for \code{"bca"}, and 200 for the two
#' bootstrap standard error methods. The quarter leaves room for
#' replications that return no value: the number required is held to the
#' request, so a simulated study that loses a few replications still reports
#' its interval, and with the quarter its limits still rest on about the
#' number required. That default is far
#' below the
#' \code{B = 10000} of \code{\link{reliability_alpha}} and
#' \code{\link{reliability_omega}} and keeps the product manageable. It
#' also makes each replication's limits less precise than the reported
#' interval's would be: a percentile limit at \code{B = 500} rests on
#' about 12 resamples in each tail, and the Monte Carlo error of the limits
#' adds to the spread of the realized widths. A check of a bootstrap
#' interval that will be reported with a large \code{B} is closer to that
#' interval with a larger \code{B} here, at a cost that grows in
#' proportion. A \code{B} below the required number is refused unless
#' \code{options(dmar.bootstrap_exploratory = TRUE)} is set; the run is
#' then marked by one warning of class \code{dmar_bootstrap_exploratory}
#' for the call, not one for each simulated study.
#'
#' \strong{The omega coefficient.} Omega is computed with
#' \code{denominator = "model_implied"}, the coefficient whose delta method
#' interval the planner evaluates and the one the Wald intervals are
#' defined for. For the parallel population simulated here the model
#' implied and observed denominators estimate the same value.
#'
#' \strong{Coverage at a true reliability of zero.} Coverage and the misses
#' in each tail are counted on the interval the estimation function
#' reports, whose limits are each clamped into [0, 1] (see
#' \code{\link{reliability_alpha}}), since that is the interval a plan is
#' made for and a study reports. Clamping can only raise the coverage of a
#' value in [0, 1], and at the boundary it does: with
#' \code{true_reliability = 0}, an alpha interval lying wholly below zero is
#' reported as [0, 0], which contains zero, so \code{pct_ci_miss_high} is
#' zero by construction and \code{total_type_I_error} is the rate of
#' misses low. A miss high there would place the reliability below zero,
#' outside the range a reliability coefficient can take. With four items,
#' \eqn{N = 40}, and Bonett's interval (400 replications,
#' \code{set.seed(1)}) the intervals miss .0075 of the time, all of it low.
#' For any true reliability above zero, an interval lying wholly below zero
#' misses high.
#'
#' \strong{Replications that fail.} A replication that returns no interval
#' (a fit that did not converge, which is more common under full
#' information maximum likelihood with small samples and much missing
#' data, a fit whose likelihood has no maximum, which the estimation
#' functions refuse with an error of class \code{dmar_no_maximum} when an
#' error variance runs off below zero while the likelihood is flat or
#' still rising, or an interval the method could not form) is set aside:
#' it counts in no summary, and one warning states how many were set
#' aside and the first reason. When no replication returns an interval,
#' the function stops with that reason rather than return a table of
#' missing values.
#'
#' @param true_reliability Population reliability coefficient (in
#'   \eqn{[0, 1)}).
#' @param supposed_reliability Reliability supposed when the study is
#'   planned, a planning assumption rather than an estimate from data; the
#'   function passes the factor loading and residual variance it implies to
#'   \code{\link{ss_aipe_reliability}} as \code{true_factor_loading} and
#'   \code{true_residual_var}.
#' @param n_items Number of items in the composite (\eqn{J} in the formulas):
#'   a whole number of at least 2, and at least 3 for
#'   \code{estimator = "omega"}.
#' @param width Desired full width of the CI on reliability.
#' @param specified_N Sample size to evaluate (incompatible with
#'   \code{supposed_reliability}).
#' @param estimator One of \code{"alpha"} (default; coefficient alpha
#'   via \code{\link{reliability_alpha}}) or \code{"omega"} (composite
#'   reliability via \code{\link{reliability_omega}}). For a population
#'   of parallel tests the two coincide; differences in
#'   sample estimates reflect estimator-specific finite-sample bias
#'   and CI behavior.
#' @param ci_method CI method passed to the estimator. The default
#'   \code{NULL} is \code{"ml"} for alpha and omega, with complete or
#'   incomplete data: the normal theory Wald interval the plan is made for
#'   (see \emph{Details}). Any other interval the estimator
#'   offers can be named.
#' @param conf_level Confidence level (default \code{0.95}).
#' @param assurance Optional assurance probability passed to the
#'   planner, whose Monte Carlo search then runs at its defaults
#'   (\code{G_initial = 500}, \code{G = 5000}) and draws from
#'   the current random number stream, so \code{set.seed()} before the call
#'   makes the plan and the simulation reproducible together. The plan is
#'   the one the planner reports for the parallel model with the normal
#'   theory interval, which the default \code{ci_method} evaluates;
#'   \code{specified_N} evaluates any other size.
#' @param G Number of Monte Carlo replications.
#' @param print_iter Logical.
#' @param filename Optional path for a comma separated file recording
#'   every replication (the sample reliability estimate, the two confidence
#'   limits, the interval width, and two indicators of whether the
#'   interval missed \code{true_reliability} below or above, followed, when
#'   \code{missing_prop} is above 0, by the fraction of missing information
#'   the replication realized for the coefficient): nothing is
#'   written when \code{filename} is \code{NULL} (the default), a new
#'   file with a header row is created otherwise, an existing file at
#'   that path is appended to, and a throwaway run should point it at
#'   \code{tempfile(fileext = ".csv")}.
#' @param fmi The fraction of missing information the plan assumes for
#'   the coefficient, passed to \code{\link{ss_aipe_reliability}} when
#'   \code{supposed_reliability} resolves the sample size: a single number
#'   in \eqn{[0, 1)}, with the default \code{0} for a complete data plan.
#'   With \code{specified_N} there is no plan for it to act on, and a value
#'   above 0 is an error.
#' @param missing_prop The probability with which each simulated item
#'   score is missing, completely at random: a single number in
#'   \eqn{[0, 1)}. The default \code{0} simulates complete data;
#'   above 0, every replication is analyzed by
#'   full information maximum likelihood (see \emph{Details}).
#' @param B Number of bootstrap replications inside each simulated study
#'   when \code{ci_method} is a bootstrap method. The default,
#'   \code{NULL}, is the larger of 200 and one and a quarter times the
#'   number the interval method and the level require (see
#'   \code{\link{dmar_bootstrap}}): 500 for \code{"percentile"}, 2500 for
#'   \code{"bca"}, and 200 for the two bootstrap standard error methods at
#'   the 95 percent level. The quarter is room for replications that return
#'   no value (a resample without a one-factor solution): the required
#'   number (400, 2000, 100) is held to the request, so a simulated study
#'   that loses a few replications still reports its interval, and with the
#'   quarter its limits still rest on about the required number. A simulated
#'   study in which more than half of the replications return no value
#'   reports no interval.
#'   Every simulated study runs \code{B} replications, so the run time is
#'   proportional to \code{G} times \code{B}. Not used by the other
#'   interval methods.
#' @param seed \code{NULL} (the default) to draw from the current random
#'   number stream, or a single number that makes the call reproducible; a
#'   supplied seed is set for this call only, and the caller's random
#'   number state is restored when the function returns.
#'
#' @return A \code{data.frame} with rows for mean / median / SD of
#'   the realized reliability and CI width, the proportion of intervals
#'   at or below \code{width}, tail-specific and overall non-coverage
#'   of \code{true_reliability}, and the input echoes \code{total_N},
#'   \code{n_items}, \code{true_reliability}, \code{supposed_reliability},
#'   \code{width}, \code{conf_level}, and \code{assurance} (present only when an
#'   assurance was supplied). The proportion of intervals at or below
#'   \code{width}, the widths, and the noncoverage proportions are over
#'   the replications that returned an interval. The \code{simulation}
#'   attribute records the numbers of replications attempted, with an
#'   interval, and without one (see the section on the record of the
#'   simulation), and, here, the message of the first failure
#'   (\code{first_failure}). When \code{fmi} or \code{missing_prop} is
#'   above 0, the rows are unchanged and the missing data inputs travel as
#'   attributes: \code{fmi} (\code{NA} with \code{specified_N}),
#'   \code{missing_prop}, and, when \code{fmi} is above 0,
#'   \code{necessary_N_complete}, the planner's complete data sample size.
#'   When \code{missing_prop} is above 0, two more attributes record the
#'   fraction of missing information the simulated data realized for the
#'   coefficient: \code{mean_realized_fmi}, its mean over the replications
#'   in which it was computed, and \code{n_realized_fmi}, the number of
#'   those replications. With \code{fmi} and \code{missing_prop} both at 0
#'   no attribute is added. With a bootstrap \code{ci_method} the result
#'   also carries the attributes \code{B}, \code{B_used},
#'   \code{n_dropped}, \code{B_required}, and \code{bootstrap_status} of
#'   \code{\link{dmar_bootstrap}}; they record the request each simulated
#'   study ran (\code{B_used} equals \code{B} and \code{n_dropped} is 0,
#'   because replications dropped inside a simulated study are not
#'   totaled), and \code{bootstrap_status} is \code{"exploratory"} when
#'   \code{B} is below \code{B_required}.
#'
#' @references
#' Bonett, D. G. (2002). Sample size requirements for testing and estimating
#'   coefficient alpha. \emph{Journal of Educational and Behavioral
#'   Statistics, 27}(4), 335--340. \doi{10.3102/10769986027004335}
#'
#' Chen, L., & Savalei, V. (2021). Three sample estimates of fraction of
#'   missing information from full information maximum likelihood.
#'   \emph{Frontiers in Psychology, 12}, 667802.
#'   \doi{10.3389/fpsyg.2021.667802}
#'
#' Kelley, K., & Pornprasertmanit, S. (2016). Confidence intervals for
#'   population reliability coefficients: Evaluation of methods,
#'   recommendations, and software for composite measures.
#'   \emph{Psychological Methods, 21}, 69--92. \doi{10.1037/a0040086}
#'
#' Savalei, V., & Rhemtulla, M. (2012). On obtaining estimates of the
#'   fraction of missing information from full information maximum
#'   likelihood. \emph{Structural Equation Modeling, 19}(3), 477--494.
#'   \doi{10.1080/10705511.2012.687669}
#'
#' Terry, L. J., & Kelley, K. (2012). Sample size planning for composite
#'   reliability coefficients: Accuracy in parameter estimation via
#'   narrow confidence intervals. \emph{British Journal of Mathematical
#'   and Statistical Psychology, 65}, 371--401.
#'   \doi{10.1111/j.2044-8317.2011.02030.x}
#'
#' Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). \emph{Designing
#'   experiments and analyzing data: A model comparison perspective}
#'   (4th ed.). Routledge.
#'
#' van Zyl, J. M., Neudecker, H., & Nel, D. G. (2000). On the distribution of
#'   the maximum likelihood estimator of Cronbach's alpha.
#'   \emph{Psychometrika, 65}(3), 271--280. \doi{10.1007/BF02296146}
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso \code{\link{ss_aipe_reliability}}, \code{\link{reliability_alpha}}, \code{\link{reliability_omega}}
#'
#' @examples
#' # Reduced Monte Carlo sweep (small G) so the example runs quickly;
#' # raise G for a production sensitivity analysis.
#' set.seed(113)
#' ss_aipe_reliability_sensitivity(
#'   true_reliability     = 0.80,
#'   supposed_reliability = 0.80,
#'   n_items = 4, width = 0.15,
#'   estimator = "alpha",
#'   G = 20, print_iter = FALSE
#' )
#'
#' # The same plan made for an FMI of .18 and checked against data with each
#' # item score missing completely at random with probability .2, analyzed
#' # by full information maximum likelihood. The second line under the table
#' # gives the FMI those data realized for the coefficient, to read beside
#' # the .18 the plan assumed. Every replication fits a model, so G = 3 only
#' # shows the call; a reported check deserves G = 1000.
#' set.seed(113)
#' ss_aipe_reliability_sensitivity(
#'   true_reliability     = 0.80,
#'   supposed_reliability = 0.80,
#'   n_items = 4, width = 0.15,
#'   estimator = "alpha",
#'   G = 3, fmi = 0.18, missing_prop = 0.2
#' )
#'
#' @keywords design htest
#'
#' @seealso \code{\link{design_consequences}} for what a chosen design delivers:
#'   power, the Type S (sign) and Type M (exaggeration) errors of the
#'   significance filter, and the expected confidence interval width.
#'
#' @family AIPE sample size planning
#'
#' @inheritSection ss_aipe_cv_sensitivity The record of the simulation
#' @export
ss_aipe_reliability_sensitivity <- function(true_reliability = NULL,
                                            supposed_reliability = NULL,
                                            n_items, width,
                                            specified_N = NULL,
                                            estimator = c("alpha", "omega"),
                                            ci_method = NULL,
                                            conf_level = 0.95,
                                            assurance = NULL,
                                            G = 1000, print_iter = FALSE,
                                            filename = NULL,
                                            fmi = 0, missing_prop = 0,
                                            B = NULL,
                                            seed = NULL) {
  .dmar_check_assurance(assurance)
  estimator <- match.arg(estimator)
  scalar <- function(x) is.numeric(x) && length(x) == 1L && is.finite(x)
  G <- .dmar_check_count(G, "G")
  .dmar_local_seed(seed)
  if (!scalar(conf_level) || conf_level <= 0 || conf_level >= 1) .dmar_refuse_argument(
    "'conf_level' must be strictly between zero and one.", arg = "conf_level")
  if (!scalar(width) || width <= 0) .dmar_refuse_argument(
    "'width' must be a positive finite number.", arg = "width")
  if (!scalar(n_items) || n_items < 2 || n_items != trunc(n_items)) .dmar_refuse_argument(
    "'n_items' must be a whole number of at least two items.", arg = "n_items")
  if (estimator == "omega" && n_items < 3) {
    .dmar_refuse_argument(c(
      "Omega sensitivity requires at least three items: a two-item ",
      "single-factor model is not identified by this estimator. ",
      "Use at least three items or plan alpha explicitly."), arg = NULL)
  }
  if (!is.null(specified_N) && (!scalar(specified_N) || specified_N < 4 || specified_N != trunc(specified_N))) {
    .dmar_refuse_argument(
      "'specified_N' must be a whole number of at least four respondents.",
      arg = "specified_N")
  }
  .aipe_check_fmi(fmi)
  .aipe_check_missing_prop(missing_prop)
  if (is.null(ci_method)) {
    # The replications evaluate the interval the plan was made for. The
    # plan is ss_aipe_reliability()'s for the parallel model with the
    # normal theory interval: a Wald interval whose standard error is van
    # Zyl, Neudecker, and Nel's (2000), the interval of
    # reliability_alpha(ci_method = "ml"), and the variance it takes at the
    # population equals that of the model implied omega, so for omega it
    # is reliability_omega(denominator = "model_implied", ci_method =
    # "ml"). With missing data both are the same Wald intervals with the
    # FIML standard error, the interval the fmi mapping plans for. Before
    # 1.1.0 the defaults were Bonett's interval for alpha (ml_logistic
    # under missing data) and mlr for omega, wider than the interval
    # planned for, so a correctly specified assurance plan fell short of
    # its assurance.
    ci_method <- "ml"
  }
  allowed <- if (estimator == "alpha") {
    c("bonett", "feldt", "fisher", "hakstian_whalen", "ml", "ml_logistic", "adf", "adf_logistic",
      "bootstrap_se", "bootstrap_se_logistic", "percentile", "bca")
  } else {
    setdiff(eval(formals(reliability_omega)$ci_method), "none")
  }
  if (!is.character(ci_method) || length(ci_method) != 1L || is.na(ci_method) || !ci_method %in% allowed) {
    .dmar_refuse_argument(c(
      "'ci_method' must name an interval available for estimator = '", estimator,
   "': ", paste(allowed, collapse = ", "), "."), arg = "ci_method")
  }
  if (missing_prop > 0) {
    fiml_ok <- if (estimator == "alpha") {
      c("ml", "ml_logistic")
    } else {
      c("ml", "mlr", "ml_logistic", "mlr_logistic")
    }
    if (!(is.character(ci_method) && length(ci_method) == 1L &&
          ci_method %in% fiml_ok)) {
      .dmar_refuse_argument(c(
        "With missing_prop > 0 the simulated data are analyzed by full ",
        "information maximum likelihood, so 'ci_method' must be an ",
        "interval whose standard error comes from the FIML information ",
        "matrix. For estimator = \"", estimator, "\" use one of: ",
        paste(sQuote(fiml_ok, FALSE), collapse = ", "), "."), arg = "ci_method")
    }
  }
  bootstrap_methods <- c("bootstrap_se", "bootstrap_se_logistic",
                         "percentile", "bca")
  # The number of replications inside each simulated study follows the one
  # rule for bootstrap replications (R/dmar_bootstrap.R): the default is the
  # larger of 200 and the number the method and the level require.
  boot_inner <- is.character(ci_method) && length(ci_method) == 1L &&
    ci_method %in% bootstrap_methods
  if (boot_inner) {
    boot_kind <- .boot_method_kind(ci_method)
    boot_tail <- (1 - conf_level) / 2
    if (is.null(B)) {
      # A quarter more than the requirement. The count is held to the
      # request, so a simulated study in which a few replications fail (no
      # one-factor solution in a resample) reports its interval either
      # way; the quarter keeps its limits on about the required number.
      B <- max(200L, as.integer(ceiling(
        1.25 * .dmar_boot_required(boot_kind, boot_tail))))
    }
    B <- .dmar_check_B(B, kind = boot_kind, tail = boot_tail,
                       conf_level = conf_level)
  }

  if (is.null(supposed_reliability) && is.null(specified_N))
    .dmar_refuse_argument(
      "You must specify either 'supposed_reliability' or 'specified_N'.",
      arg = c("supposed_reliability", "specified_N"))
  if (!is.null(supposed_reliability) && !is.null(specified_N))
    .dmar_refuse_argument(c(
      "You must specify 'supposed_reliability' or 'specified_N', but not ",
      "both."), arg = c("supposed_reliability", "specified_N"))
  if (is.null(true_reliability) || !is.numeric(true_reliability) ||
      true_reliability < 0 || true_reliability >= 1)
    .dmar_refuse_argument(
      "'true_reliability' must be a single value in [0, 1).",
      arg = "true_reliability")
  if (!is.numeric(n_items) || length(n_items) != 1L || n_items < 2)
    .dmar_refuse_argument(
      "'n_items' (number of items) must be a single integer >= 2.",
      arg = "n_items")
  if (!is.null(specified_N) && fmi > 0)
    .dmar_refuse_argument(c(
      "'fmi' is a planning value: it enlarges the N that ",
      "'supposed_reliability' resolves through ss_aipe_reliability(). ",
      "With 'specified_N' there is no plan for it to act on; to simulate ",
      "missing data at that N, use 'missing_prop'."),
      arg = c("fmi", "supposed_reliability", "specified_N"))
  .check_filename(filename)

  # Solve for the lambda^2 and psi^2 of parallel tests that produce the
  # requested rho at unit per-item variance.
  rho_solve <- function(rho) {
    u <- rho / (n_items * (1 - rho) + rho)  # lambda^2
    list(lambda = sqrt(u), psi_square = 1 - u)
  }

  if (!is.null(supposed_reliability)) {
    pars <- rho_solve(supposed_reliability)
    # An assurance plan is searched at the planner's default numbers of
    # Monte Carlo iterations, the plan a user would report. Before 1.1.0 the
    # search ran at 50 and 200 iterations, whose noise let
    # a size below the one the assurance needs pass by chance, so the size
    # evaluated was a few cases short of the reported plan.
    plan <- suppressMessages(suppressWarnings(
      ss_aipe_reliability(model = "parallel", type = "Normal Theory",
                          width = width, n_items = n_items,
                          true_factor_loading = pars$lambda,
                          true_residual_var = pars$psi_square,
                          conf_level = conf_level,
                          assurance = assurance,
                          fmi = fmi)
    ))
    n <- plan$value[plan$term == "necessary_N"]
  } else {
    n <- specified_N
  }
  n <- as.integer(n)
  if (n < 4) .dmar_refuse_argument("Resolved sample size is < 4.", arg = NULL)

  pars_true <- rho_solve(true_reliability)
  lambda_t  <- pars_true$lambda
  psi_t     <- sqrt(pars_true$psi_square)

  rel_hat   <- numeric(G)
  ci_lo     <- numeric(G)
  ci_hi     <- numeric(G)
  ci_w      <- numeric(G)
  tI_lower  <- logical(G)
  tI_upper  <- logical(G)
  # The FMI each replication realized for the coefficient, recorded only
  # under missing data (with complete data it is zero by construction).
  fmi_hat   <- if (missing_prop > 0) rep(NA_real_, G) else NULL

  # The width, width-target, and coverage summaries are over the
  # replications that returned an interval; the simulation attribute
  # carries the numbers attempted, with an interval, and without one.
  first_failure <- NULL
  # Whether the first failure was an interval withheld (limits NA) and not
  # an error.
  withheld_first <- FALSE
  for (g in seq_len(G)) {
    if (isTRUE(print_iter)) cat(g, "\n")
    T_i <- stats::rnorm(n, mean = 0, sd = 1)
    E   <- matrix(stats::rnorm(n * n_items, mean = 0, sd = psi_t), nrow = n, ncol = n_items)
    Y   <- matrix(lambda_t * T_i, nrow = n, ncol = n_items, byrow = FALSE) + E
    if (missing_prop > 0) Y <- .mcar_delete(Y, missing_prop)
    Y_df <- as.data.frame(Y)

    res <- tryCatch({
      if (estimator == "alpha") {
        if (missing_prop > 0) {
          # lavaan's warnings about low pairwise coverage and about the
          # convergence of a saturated FIML fit would repeat once per
          # replication; a fit that fails is counted below instead.
          suppressWarnings(
            reliability_alpha(data = Y_df, missing = "fiml",
                              ci_method = ci_method,
                              conf_level = conf_level, B = B)
          )
        } else {
          # A sample alpha below zero, possible with few cases and low
          # reliability, is warned about by reliability_alpha() on every
          # call; here it is one replication among G, and its limits,
          # clamped into [0, 1], enter the summaries like any other. The
          # simulated items are complete, and missing = "listwise" says so
          # rather than leaving it to the default's resolution.
          withCallingHandlers(
            reliability_alpha(data = Y_df, missing = "listwise",
                              ci_method = ci_method,
                              conf_level = conf_level, B = B),
            dmar_negative_reliability = function(w) {
              invokeRestart("muffleWarning")
            },
            # An exploratory run is marked once for the call, after the
            # loop, not once for each simulated study.
            dmar_bootstrap_exploratory = function(w) {
              invokeRestart("muffleWarning")
            })
        }
      } else {
        # The model implied denominator is the coefficient whose delta
        # method interval the planner evaluates, and the one the Wald
        # intervals (mlr, ml, and their logistic forms) are defined for.
        suppressMessages(suppressWarnings(
          reliability_omega(data = Y_df, ci_method = ci_method,
                            denominator = "model_implied",
                            missing = if (missing_prop > 0) "fiml"
                                      else "listwise",
                            conf_level = conf_level, B = B)
        ))
      }
    }, error = function(e) e)

    r <- NA_real_; lo <- NA_real_; hi <- NA_real_
    if (inherits(res, "error")) {
      # A replication whose likelihood has no maximum is recorded by its
      # short reason: the refusal's lavaan code is that of one simulated
      # sample (.no_max_short()).
      if (is.null(first_failure)) {
        first_failure <- if (inherits(res, "dmar_no_maximum")) {
          .no_max_short(.no_max_fields(res))
        } else {
          conditionMessage(res)
        }
      }
    } else {
      # reliability_alpha() / reliability_omega() return the point estimate
      # under the term "estimate" (not "alpha" or "omega"), with
      # "ci_lower" / "ci_upper" rows.
      pick <- function(needle) {
        hit <- res$value[res$term == needle]
        if (length(hit) == 0) NA_real_ else as.numeric(hit[1])
      }
      r  <- pick("estimate")
      lo <- pick("ci_lower")
      hi <- pick("ci_upper")
      if (!(is.finite(lo) && is.finite(hi)) && is.null(first_failure)) {
        first_failure <- paste0("ci_method = \"", ci_method, "\" returned ",
                                "no confidence limits.")
        withheld_first <- TRUE
      }
    }
    ok <- is.finite(lo) && is.finite(hi)
    # Coverage is counted on the reported limits, each clamped into [0, 1]
    # by the estimation function. At true_reliability = 0 an alpha interval
    # wholly below zero is [0, 0] and covers zero, so it cannot miss high;
    # through 1.0.1 its inverted limits counted as a miss high. A miss high
    # there would put the reliability below zero, outside its range, and
    # the interval evaluated is the one a study would report (see the help).
    rel_hat[g]  <- r
    ci_lo[g]    <- lo
    ci_hi[g]    <- hi
    ci_w[g]     <- hi - lo
    tI_lower[g] <- if (ok) true_reliability < lo else NA
    tI_upper[g] <- if (ok) true_reliability > hi else NA
    # The FMI is taken for the replications that enter the summaries.
    if (ok && missing_prop > 0) {
      fmi_hat[g] <- .aipe_realized_fmi(Y_df, estimator)
    }
  }

  # An exploratory run (options(dmar.bootstrap_exploratory = TRUE) with a
  # B below the requirement) is marked once for the call.
  if (boot_inner && .dmar_boot_exploratory() &&
      B < .dmar_boot_required(boot_kind, boot_tail)) {
    .dmar_boot_signal_text(.dmar_boot_settle(
      B, B, kind = boot_kind, tail = boot_tail, conf_level = conf_level,
      fun = "ss_aipe_reliability_sensitivity",
      quantity = "the coefficient", dropped_because = "returned no value",
      signal = FALSE))
  }

  n_failed <- sum(is.na(tI_lower))
  if (n_failed > 0L) {
    # A bootstrap interval is withheld in a simulated study in which
    # fewer than half of the replications returned a value (the share
    # rule; the count rule is held to the request, which was checked),
    # and a BCa interval also where its adjusted limits take more
    # replications than were requested.
    raise_B <- if (boot_inner && isTRUE(withheld_first)) {
      paste0(" A bootstrap interval is withheld in a simulated study in ",
             "which fewer than half of its ", B, " replications returned ",
             "a value, which a larger B does not change",
             if (identical(boot_kind, "bca"))
               paste0(", and a BCa interval also where the bias correction ",
                      "and the acceleration place a limit so far into the ",
                      "tail that ", B, " replications put fewer than 10 ",
                      "beyond it, which a larger B does"),
             " (see ?dmar_bootstrap).")
    }
    # The class and the fields are those the coefficient of variation and
    # Cliff's delta siblings give for the same event.
    .dmar_warn(c(n_failed, " of ", G, " replications returned no confidence ",
                 "interval. The width, width-target, and coverage rows describe ",
                 "the ", G - n_failed, " replications that returned one. The ",
                 "first failure: ", first_failure, raise_B),
               "replications_without_interval",
               n_without_interval = n_failed, G = G)
  }

  if (!is.null(filename)) {
    per_rep <- data.frame(rel_hat = rel_hat, ci_lower = ci_lo,
                          ci_upper = ci_hi, ci_width = ci_w,
                          type_I_lower = tI_lower, type_I_upper = tI_upper)
    if (missing_prop > 0) per_rep$realized_fmi <- fmi_hat
    .write_sensitivity_csv(per_rep, filename)
  }

  out <- data.frame(
    term  = c("mean_reliability", "median_reliability", "sd_reliability",
              "mean_ci_width", "median_ci_width", "sd_ci_width",
              "pct_ci_less_width",
              "pct_ci_miss_low", "pct_ci_miss_high", "total_type_I_error",
              "total_N", "n_items",
              "true_reliability", "supposed_reliability",
              "width", "conf_level",
              if (!is.null(assurance)) "assurance"),
    value = c(mean(rel_hat, na.rm = TRUE),
              stats::median(rel_hat, na.rm = TRUE),
              stats::sd(rel_hat, na.rm = TRUE),
              mean(ci_w, na.rm = TRUE),
              stats::median(ci_w, na.rm = TRUE),
              stats::sd(ci_w, na.rm = TRUE),
              sum(is.finite(ci_w) & ci_w <= width) / sum(is.finite(ci_w)),
              mean(tI_lower, na.rm = TRUE),
              mean(tI_upper, na.rm = TRUE),
              mean(tI_lower | tI_upper, na.rm = TRUE),
              n, n_items,
              true_reliability,
              if (is.null(supposed_reliability)) NA_real_ else
                supposed_reliability,
              width, conf_level,
              if (!is.null(assurance)) assurance),
    stringsAsFactors = FALSE, row.names = NULL
  )
  out$value[is.nan(out$value)] <- NA_real_
  out <- .as_dmar_tbl(out, conf_level = conf_level)
  # A replication has an interval exactly when its width is finite (both
  # limits finite), which is also when its two miss indicators are not NA.
  simulation <- .ss_sens_simulation(G, is.finite(ci_w),
                                    is.finite(ci_w) & ci_w <= width,
                                    conf_level)
  attr(out, "simulation") <- append(
    simulation, list(first_failure = first_failure),
    after = match("width_target_probability", names(simulation)))

  # The missing data inputs travel as attributes, not rows, so the table
  # keeps the family's schema; with both at 0 nothing is added and the
  # result is the one the function returned before they existed.
  if (fmi > 0 || missing_prop > 0) {
    attr(out, "fmi") <- if (is.null(supposed_reliability)) NA_real_ else fmi
    attr(out, "missing_prop") <- missing_prop
    if (fmi > 0) {
      attr(out, "necessary_N_complete") <- attr(plan, "necessary_N_complete")
    }
    if (missing_prop > 0) {
      n_fmi <- sum(is.finite(fmi_hat))
      attr(out, "mean_realized_fmi") <- if (n_fmi > 0L) {
        mean(fmi_hat, na.rm = TRUE)
      } else {
        NA_real_
      }
      attr(out, "n_realized_fmi") <- n_fmi
    }
  }
  # The request each simulated study ran (the exploratory warning was given
  # above, before the count of missing intervals).
  if (boot_inner) {
    out <- .dmar_boot_mark_planner(
      out, B, kind = boot_kind, tail = boot_tail, conf_level = conf_level,
      fun = "ss_aipe_reliability_sensitivity", quantity = "the coefficient",
      signal = FALSE)
  }
  out
}
