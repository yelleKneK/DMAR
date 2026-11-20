# Sample size for AIPE on a semipartial (part) correlation.
#' Sample Size for AIPE on a Semipartial (Part) Correlation
#'
#' Takes a planning value of the population semipartial correlation
#' \eqn{\rho_{Y(X \cdot Z_1 \cdots Z_{p-1})}} (the correlation of \eqn{Y} with
#' the part of \eqn{X} that \eqn{Z_1, \ldots, Z_{p-1}} do not explain, with
#' \eqn{Y} not residualized), the number of predictors (\emph{X} and the variables
#' partialed out), and a desired confidence interval width, and returns the smallest sample size
#' at which the interval is expected to be that narrow (accuracy in
#' parameter estimation, AIPE; Kelley & Maxwell, 2003). With
#' \code{assurance} it returns, by simulation, the smallest sample size at
#' which the realized interval is that narrow with about the stated
#' probability (the search and the probability reported are both Monte
#' Carlo estimates; see Details).
#'
#' @param true_r_sp Population semipartial correlation the sample size is
#'   planned for, in \eqn{(-1, 1)}. The planner treats this value as the
#'   truth; \code{\link{ss_aipe_semipartial_r_sensitivity}} shows what the
#'   plan delivers when the truth (\code{true_r_sp} there) differs from the
#'   value planned with (\code{supposed_r_sp} there).
#' @param p Number of predictors in the regression of \emph{Y} on \emph{X}
#'   and the variables partialed out, so \code{p - 1} variables are
#'   partialed out and \code{p} is at least 2. A categorical variable with
#'   \emph{a} levels counts as \emph{a} - 1 predictors. It is the \code{p}
#'   of \code{\link{ss_power_R2}}, \code{\link{ci_R2}}, and
#'   \code{\link{ci_reg_coef}}: the correlation concerns one predictor of a
#'   regression with \code{p} predictors. In the notation of a model
#'   comparison, \eqn{q} is the number of predictors removed from the full
#'   model to give the restricted model (\eqn{q = p_{full} - p_{restricted}},
#'   the predictors being tested), and the number of variables partialed out
#'   is \eqn{p - q}; this correlation concerns one predictor, so \eqn{q = 1}
#'   and \eqn{p - 1} variables are partialed out.
#' @param width Desired full width of the confidence interval on the
#'   semipartial correlation.
#' @param which_width Whether \code{width} refers to the \code{"Full"}
#'   width (default) or the \code{"Lower"}/\code{"Upper"} half-width.
#' @param conf_level Desired confidence level (default \code{0.95}).
#' @param assurance Optional. The probability that the realized interval
#'   is no wider than \code{width}. The function returns the smallest
#'   sample size, at or above the sample size planned without it, at which
#'   the simulated probability is at least \code{assurance} (see "The
#'   assurance step" in Details). A single number above 0.5 and below 1;
#'   when \code{NULL} (the default) the plan is for the expected width
#'   only and nothing is simulated.
#' @param true_R2_full Population squared multiple correlation of the full
#'   model the sample size is planned for, the model with \eqn{X} and every
#'   variable partialed out (\code{true_R2_Y_X} in
#'   \code{\link{ss_aipe_reg_coef}}), treated as the truth as
#'   \code{true_r_sp} is, in
#'   \eqn{[0, 1]} and not smaller than \code{true_r_sp} squared (a smaller
#'   value is refused with an error of class
#'   \code{dmar_semipartial_r_incoherent}, unless the pair could have come
#'   from rounding to two decimals, when the plan proceeds as if
#'   \code{true_R2_full} were \code{true_r_sp} squared and gives a warning of class
#'   \code{dmar_semipartial_r_rounded}; see
#'   \code{\link{var_semipartial_r}}). When \code{NULL} (the default),
#'   the plan uses the largest variance over every full-model \eqn{R^2}
#'   that is possible with \code{true_r_sp}, so the sample size is sufficient
#'   whatever that \eqn{R^2} is.
#'
#' @param G Number of simulated samples at each sample size the assurance
#'   step tries, and in the final estimate of the attained probability; a
#'   whole number of at least 100 (default 20000). A smaller value is
#'   refused, because the step estimates a probability from these samples
#'   and with fewer than 100 its standard error can exceed .05. Used and
#'   checked only when \code{assurance} is given.
#' @param seed Optional single number. When given, the assurance step is
#'   reproducible, and the random number stream of the calling session is
#'   left as it was. With \code{NULL} (the default) the draws come from
#'   the session's current stream. Used only when \code{assurance} is
#'   given.
#'
#' @return A \code{data.frame} with the columns \code{term} and
#'   \code{value} and one row per term:
#'   \describe{
#'     \item{\code{necessary_N}}{the planned total sample size.}
#'     \item{\code{expected_width}}{the full width of the interval at
#'       \code{necessary_N} when the two sample values equal the planning
#'       values.}
#'     \item{\code{true_r_sp}, \code{p}, \code{width},
#'       \code{conf_level}}{the inputs, with \code{width} the
#'       \code{width} as given.}
#'     \item{\code{true_R2_full}}{the planning value of the full-model
#'       \eqn{R^2}; \code{NA} when none was supplied.}
#'     \item{\code{assurance}}{the requested probability; present only
#'       when \code{assurance} was given, as are the next two rows.}
#'     \item{\code{assurance_attained}}{the proportion of \code{G} fresh
#'       simulated samples of size \code{necessary_N} whose interval was no
#'       wider than the target. Its Monte Carlo standard error is the
#'       attribute \code{assurance_mc_se}.}
#'     \item{\code{G}}{the number of simulated samples.}
#'   }
#'
#' @details
#' \strong{Estimand and interval.} The quantity planned for is the
#' population semipartial correlation of \eqn{Y} with \eqn{X} given
#' \eqn{Z_1, \ldots, Z_{p-1}}, for independent observations from a
#' multivariate normal population. The interval is the two-sided Wald
#' interval at \code{conf_level}, the estimate plus and minus
#' \eqn{z_{1 - \alpha/2} \sqrt{\hat V / (N - p)}}, where \eqn{\hat V}
#' is the first-order variance below evaluated at the sample semipartial
#' correlation and the sample full-model \eqn{R^2}
#' (\code{\link{var_semipartial_r}} with \code{under_null = FALSE}). Its
#' limits are not truncated at \eqn{-1} and \eqn{1}.
#'
#' \strong{Asymptotic variance of the semipartial.} Under multivariate
#' normality, to first order, the sample semipartial correlation has
#' variance \eqn{V / N}, with
#' \deqn{V \;=\; (1 - \rho_{sp}^2)^2 \;-\; P_r^2\,(1 - 2 \rho_{sp}^2),}
#' where \eqn{\rho_{sp}} is the population semipartial correlation
#' (\code{true_r_sp}) and
#' \eqn{P_r^2} is the full-model \eqn{R^2} (\code{true_R2_full}) minus \eqn{\rho_{sp}^2}, that is,
#' the squared multiple correlation of \eqn{Y} with the variables
#' partialed out. This is the first-order variance of Aloe and Becker
#' (2012, Equation 4), which they obtain from Alf and Graf (1999), with
#' \eqn{N - p} in place of their \eqn{n}: in their notation the numerator is
#' \eqn{P^4 - 2 P^2 + P_{(f)}^2 + 1 - P_{(f)}^4}, with \eqn{P^2} the
#' full-model and \eqn{P_{(f)}^2} the restricted-model squared multiple
#' correlation, which equals \eqn{V} with \eqn{P_r^2 = P_{(f)}^2}. It follows
#' from the delta method variance of the difference between two nested
#' squared multiple correlations (Model
#' A of Olkin & Finn, 1995, for one variable partialed out; Alf & Graf,
#' 1999, for any number; see \code{\link{var_semipartial_r}}).
#' The planner uses the divisor \eqn{N - p}, as
#' \code{\link{var_semipartial_r}} does. Inverting for the sample size
#' needed to achieve a target half-width \eqn{w_{1/2}} at confidence level
#' \eqn{1 - \alpha}:
#' \deqn{N \;=\; p + \Big\lceil
#'         z_{1 - \alpha/2}^{2} \cdot V / w_{1/2}^{2} \Big\rceil.}
#'
#' \strong{Planning without \code{true_R2_full}.} \eqn{V} decreases in \eqn{P_r^2}
#' when \eqn{\rho_{sp}^2 < 1/2} and increases in \eqn{P_r^2} when
#' \eqn{\rho_{sp}^2 > 1/2}, so its largest value over every possible
#' full-model \eqn{R^2} is
#' \deqn{\max\{(1 - \rho_{sp}^2)^2,\; \rho_{sp}^2 (1 - \rho_{sp}^2)\},}
#' and the planner uses that value when \code{true_R2_full} is \code{NULL}. For
#' \eqn{\rho_{sp}^2 \le 1/2} it is \eqn{(1 - \rho_{sp}^2)^2}, the variance
#' when the variables partialed out explain none of the outcome; Algina and
#' Olejnik (2003, p. 317) note, for the squared semipartial correlation,
#' that the required sample size is largest when the squared multiple
#' correlation of the restricted model is zero. A planning value of
#' \code{true_R2_full} then gives a smaller sample size. For
#' \eqn{\rho_{sp}^2 > 1/2} the largest value is
#' \eqn{\rho_{sp}^2 (1 - \rho_{sp}^2)}, approached as the full-model
#' \eqn{R^2} approaches 1.
#'
#' \strong{Comparison with partial-r planning.} The partial correlation
#' \eqn{r_{XY \cdot Z}} divides the covariance after residualizing both
#' \eqn{X} and \eqn{Y} on \eqn{Z}; the semipartial divides after
#' residualizing only \eqn{X}. The semipartial is the natural effect size
#' companion to a standardized regression coefficient: its square equals
#' the \eqn{\Delta R^2} contributed by \eqn{X} above and beyond the
#' variables partialed out. See \code{\link{var_semipartial_r}} for the asymptotic
#' variance, and \code{\link{ss_aipe_partial_r}} for the partial-correlation
#' analog of this function.
#'
#' \strong{Limits of the plan.} The plan is a first-order one. When the
#' returned sample size is small, which happens when \eqn{|\rho_{sp}|} is
#' large and \code{true_R2_full} is near \code{true_r_sp} squared, the first-order
#' variance understates the variance of the estimate and the realized
#' intervals are wider than planned; check such a plan with
#' \code{\link{ss_aipe_semipartial_r_sensitivity}}. The planned sample
#' size is also small when \code{true_R2_full} is near 1, for any \code{true_r_sp},
#' and the same check applies: the plan for \code{true_r_sp = .1},
#' \code{p = 4}, a width of .20, and \code{true_R2_full = .99} is 12, at which
#' the intervals missed the population value in .13 of 200,000 simulated
#' samples. When \eqn{|\rho_{sp}|} is above about .9 and the variables
#' partialed out explain some of the outcome, the interval misses more
#' often on one side than the other even with several hundred observations
#' (in one simulated condition, \eqn{\rho_{sp} = .95}, a full-model
#' \eqn{R^2} of .97, \eqn{p = 2}, and \eqn{N = 395}, it missed in .064 of
#' samples, .051 of them on one side). A plan with \code{true_r_sp = 0} and
#' \code{true_R2_full = 1} is refused, because the first-order variance is zero
#' there and no sample size follows from it.
#'
#' \strong{The assurance step.} A plan without \code{assurance} makes
#' the interval as narrow as the target when the two sample values equal
#' their planning values; in roughly half of all samples it is wider. With
#' \code{assurance} the function looks for the sample size at which the
#' realized interval is no wider than the target with that
#' probability. The width of the interval depends on two statistics, the
#' sample semipartial correlation and the sample full-model \eqn{R^2},
#' and no closed form for their joint distribution is used here, so the
#' probability is simulated:
#' \enumerate{
#'   \item \emph{What is simulated.} Samples of size \eqn{N} from a
#'     multivariate normal population with semipartial correlation
#'     \code{true_r_sp} and full-model \eqn{R^2} \code{true_R2_full}: \eqn{X},
#'     \eqn{Z_1, \ldots, Z_{p-1}}, and an error \eqn{E} are independent
#'     standard normal, and
#'     \eqn{Y = \rho_{sp} X + P_r\, Z_1 + \sqrt{1 - \rho_{sp}^2 - P_r^2}\, E}.
#'     Under normality the distribution of the two statistics depends on
#'     the population only through \eqn{\rho_{sp}}, the full-model
#'     \eqn{R^2}, \eqn{p}, and \eqn{N} (see
#'     \code{\link{ss_aipe_semipartial_r_sensitivity}}), so this population
#'     stands for all of them. The sums of squares and products of the
#'     \eqn{p + 1} independent variables are Wishart with \eqn{N - 1}
#'     degrees of freedom, and the function draws the few entries of that
#'     matrix's triangular factor that the two statistics depend on (the
#'     Bartlett decomposition) in place of the data.
#'   \item \emph{The interval.} From each draw, the sample semipartial
#'     correlation and sample full-model \eqn{R^2} give the half-width
#'     \eqn{z_{1 - \alpha/2} \sqrt{\hat V / (N - p)}}, the interval
#'     that is planned and that
#'     \code{\link{ss_aipe_semipartial_r_sensitivity}} forms from data.
#'   \item \emph{The search.} The simulated probability at \eqn{N} is the
#'     share of \code{G} half-widths no larger than the target half-width.
#'     Starting from the sample size planned without \code{assurance}, the
#'     step upward doubles until that share reaches \code{assurance}, and
#'     a bisection finds the first sample size at which it does, with
#'     fresh draws at every sample size tried. The answer is never below
#'     the plan without \code{assurance}.
#'   \item \emph{The attained probability.} The share that chose
#'     \eqn{N} is biased upward by the choice, so
#'     \code{assurance_attained} is estimated again at the returned
#'     \eqn{N} from \code{G} draws that had no part in choosing it. Its
#'     standard error, \eqn{\sqrt{A (1 - A) / G}} with \eqn{A} that proportion, is the attribute
#'     \code{assurance_mc_se}.
#' }
#' Plans made this way were checked from data with
#' \code{\link{ss_aipe_semipartial_r_sensitivity}} in 11 conditions (5,000
#' simulated data sets each; a Monte Carlo standard error of .003 to
#' .006): the proportion of intervals no wider than the target was .798 to
#' .819 for a requested .80, .854 to .858 for .85, .906 to .908 for .90,
#' and .953 for .95.
#'
#' The step has properties a closed-form plan does not.
#' \itemize{
#'   \item The answer has Monte Carlo error. Without \code{seed}, two
#'     calls can return sample sizes that differ by a few cases (eight
#'     unseeded plans at the default \code{G} for one condition gave 252
#'     to 254, and 40 seeded plans of each of two other conditions gave
#'     470 to 473 and 306 to 310). The probability at the smallest of
#'     those sizes was below the request (.897 and .895 for a requested
#'     .90), so a plan is the size at which the probability is about, and
#'     not at least, the request. Give \code{seed} for a plan that can be
#'     reproduced, and report it with \code{G}.
#'   \item \code{assurance_attained} can fall below \code{assurance} by
#'     a few of its standard errors, since it is estimated from draws that
#'     did not choose \eqn{N}. Its standard error is .003 at the default
#'     \code{G} and an assurance of .80; in 200 seeded plans the lowest
#'     value was .009 below a request of .90. It is an unbiased estimate
#'     of the probability at the sample size returned. A larger \code{G}
#'     narrows both this and the spread of the sample sizes.
#'   \item Without \code{true_R2_full} the simulated population is the one the
#'     plan without \code{assurance} is made for: variables partialed out
#'     that explain none of the outcome (\eqn{P_r^2 = 0}) when \code{true_r_sp}
#'     squared is at most 1/2, and a full-model \eqn{R^2} of 1
#'     (\eqn{P_r^2 = 1 - \rho_{sp}^2}, no error term) when it is above 1/2.
#'     The second is the limit that gives the largest first-order
#'     variance; \code{\link{ss_aipe_semipartial_r_sensitivity}} refuses
#'     a population \eqn{R^2} of exactly 1, so check such a plan there
#'     with a value just below 1 or give \code{true_R2_full}.
#'   \item \code{which_width = "Lower"} or \code{"Upper"} treats
#'     \code{width} as a half-width, in the assurance step as in the plan.
#'   \item \code{G} and \code{seed} are checked and used only when
#'     \code{assurance} is given; without it nothing is drawn and the
#'     three rows are absent.
#' }
#'
#' @references
#' Alf, E. F., Jr., & Graf, R. G. (1999). Asymptotic confidence limits for
#'   the difference between two squared multiple correlations: A simplified
#'   approach. \emph{Psychological Methods, 4}(1), 70--75.
#'   \doi{10.1037/1082-989X.4.1.70}
#'
#' Algina, J., & Olejnik, S. (2003). Sample size tables for correlation
#'   analysis with applications in partial correlation and multiple
#'   regression analysis. \emph{Multivariate Behavioral Research, 38}(3),
#'   309--323. \doi{10.1207/s15327906mbr3803_02}
#'
#' Aloe, A. M., & Becker, B. J. (2012). An effect size for regression
#'   predictors in meta-analysis. \emph{Journal of Educational and Behavioral
#'   Statistics, 37}(2), 278--297. \doi{10.3102/1076998610396901}
#'
#' Cohen, J., Cohen, P., West, S. G., & Aiken, L. S. (2003). \emph{Applied
#'   multiple regression/correlation analysis for the behavioral
#'   sciences} (3rd ed.). Lawrence Erlbaum.
#'
#' Kelley, K., & Maxwell, S. E. (2003). Sample size for multiple
#'   regression: Obtaining regression coefficients that are accurate,
#'   not simply significant. \emph{Psychological Methods, 8}(3),
#'   305--321. \doi{10.1037/1082-989X.8.3.305}
#'
#' Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027).
#'   \emph{Designing experiments and analyzing data: A model comparison
#'   perspective} (4th ed.). Routledge. (See Chapter 3 on the one-way
#'   ANOVA and Chapter 4 on contrasts.)
#'
#' Olkin, I., & Finn, J. D. (1995). Correlations redux.
#'   \emph{Psychological Bulletin, 118}(1), 155--164.
#'   \doi{10.1037/0033-2909.118.1.155}
#'
#' @seealso \code{\link{var_semipartial_r}}, \code{\link{ss_aipe_partial_r}},
#'   \code{\link{ss_aipe_R2}}, \code{\link{ss_power_R2_change}}
#'
#' @examples
#' # 1. Plan N so the 95% CI on r_sp (three variables partialed out, so
#' #        p = 4 predictors) has full width <= 0.15 when the planning
#' #        value of the semipartial is 0.25, whatever the full-model R^2.
#' ss_aipe_semipartial_r(true_r_sp = 0.25, p = 4, width = 0.15)
#'
#' # 2. With a planning value of 0.42 for the full-model R^2, fewer are
#' #        needed.
#' ss_aipe_semipartial_r(true_r_sp = 0.25, p = 4, width = 0.15,
#'                       true_R2_full = 0.42)
#'
#' # 3. With an assurance of .80, found by simulation. G = 2000 keeps the
#' #        example quick; a reported plan deserves the default G = 20000.
#' #        The last three rows are the request, the proportion of fresh
#' #        simulated samples that met the width at the returned N, and G.
#' plan <- ss_aipe_semipartial_r(true_r_sp = 0.25, p = 4, width = 0.15,
#'                               assurance = 0.80, G = 2000, seed = 113)
#' plan
#' attr(plan, "assurance_mc_se")
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @keywords design
#'
#' @seealso \code{\link{design_consequences}} for what a chosen design delivers:
#'   power, the Type S (sign) and Type M (exaggeration) errors of the
#'   significance filter, and the expected confidence interval width.
#'
#' @family AIPE sample size planning
#'
#' @export

ss_aipe_semipartial_r <- function(true_r_sp, p, width,
                                  which_width = c("Full", "Lower", "Upper"),
                                  conf_level = 0.95,
                                  assurance = NULL,
                                  true_R2_full = NULL,
                                  G = 20000, seed = NULL) {
  if (missing(p))
    .dmar_refuse_argument(
      "'p' (the number of predictors) is required.", arg = "p")
  p <- .check_partial_p(p)
  # The arithmetic below is written in J, the number of variables partialed
  # out (N - J - 1 = N - p), and in r_sp and R2_full, the planning values.
  J <- p - 1
  r_sp <- true_r_sp
  R2_full <- true_R2_full
  which_width <- match.arg(which_width)
  if (!is.numeric(r_sp) || length(r_sp) != 1L || is.na(r_sp) ||
      abs(r_sp) >= 1)
    .dmar_refuse_argument(
      "'true_r_sp' must be a single value in (-1, 1).", arg = "true_r_sp")
  if (!is.numeric(width) || length(width) != 1L || !is.finite(width) ||
      width <= 0)
    .dmar_refuse_argument(
      "'width' must be a single positive number.", arg = "width")
  if (!is.numeric(conf_level) || length(conf_level) != 1L ||
      !is.finite(conf_level) || conf_level <= 0 || conf_level >= 1)
    .dmar_refuse_argument(
      "'conf_level' must be a single number in (0, 1).", arg = "conf_level")
  if (!is.null(assurance) &&
      (!is.numeric(assurance) || length(assurance) != 1L ||
       !is.finite(assurance) || assurance <= 0.5 || assurance >= 1))
    .dmar_refuse_argument(
      "'assurance' must be a single number in (0.5, 1).", arg = "assurance")
  # The step estimates a probability from G simulated samples, so G has a
  # floor: with fewer than 100 the standard error of that estimate can
  # exceed .05 (G = 1 returned an attained assurance of 0 or 1).
  if (!is.null(assurance)) {
    G <- .dmar_check_count(G, "G", min = 100L, reason = paste0(
      "The assurance step estimates a probability from 'G' simulated ",
      "samples at each sample size it tries; with 100 the standard error ",
      "of that estimate is at most .05, and the default of 20000 brings ",
      "it to .004 or less."))
  }

  # The first-order variance times n: at the planning value of R2_full when
  # one is given, otherwise its largest value over every R2_full that is
  # possible with r_sp. The second of the two worst-case forms differs from
  # the 1.0.x plan; the two are equal at r_sp^2 = 1/2, and sqrt(.5)^2 > 1/2
  # in floating point, so the comparison is of the numerators with a margin
  # and not of r_sp^2 with 1/2.
  worst_case_above_half <- is.null(R2_full) &&
    r_sp^2 * (1 - r_sp^2) > (1 - r_sp^2)^2 + 1e-12
  if (is.null(R2_full)) {
    numerator <- .semipartial_r_var_numerator_max(r_sp)
    # The population the worst case stands for, for the assurance step.
    q <- if (worst_case_above_half) 1 - r_sp^2 else 0
  } else {
    .check_semipartial_R2_full(r_sp, R2_full, arg = "true_R2_full",
                               r_arg = "true_r_sp")
    numerator <- .semipartial_r_var_numerator(r_sp, R2_full)
    # The first-order variance is 0 only at r_sp = 0 with R2_full = 1, where
    # the plan would be n = J + 1 with an expected width of 0 / 0.
    if (numerator <= 0)
      .dmar_refuse_argument(c(
        "'true_R2_full' = 1 with 'true_r_sp' = 0 leaves no variance to ",
        "plan for: the first-order variance of the semipartial ",
        "correlation is zero there. Use a planning value of ",
        "'true_R2_full' below 1."), arg = c("true_R2_full", "true_r_sp"))
    q <- max(R2_full - r_sp^2, 0)
  }

  half_width <- if (which_width == "Full") width / 2 else width
  z_alpha    <- stats::qnorm(1 - (1 - conf_level) / 2)

  # n - J - 1 = z^2 V / w^2
  n <- J + 1 + ceiling((z_alpha^2 * numerator) / (half_width^2))

  # The assurance step, a simulation (.semipartial_r_assurance_n() in
  # R/partial_r_internals.R): the smallest n found, at or above the plan
  # just made, at which the simulated probability of an interval no wider
  # than the target is at least `assurance`, and a fresh estimate of that
  # probability at the n returned.
  step <- NULL
  if (!is.null(assurance)) {
    .dmar_local_seed(seed)
    step <- .semipartial_r_assurance_n(n, r_sp, q, J, half_width, z_alpha,
                                       assurance, G)
    n <- step$n
  }

  expected_w <- 2 * z_alpha * sqrt(numerator / (n - J - 1))

  out <- data.frame(
    term  = c("necessary_N", "expected_width", "true_r_sp", "p",
              "width", "conf_level", "true_R2_full",
              if (!is.null(assurance))
                c("assurance", "assurance_attained", "G")),
    value = c(n, expected_w, r_sp, p, width, conf_level,
              if (is.null(R2_full)) NA_real_ else R2_full,
              if (!is.null(assurance)) c(assurance, step$attained, G)),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  out <- .as_dmar_tbl(out, conf_level = conf_level, subclass = "dmar_ss_aipe")
  if (!is.null(assurance)) attr(out, "assurance_mc_se") <- step$mc_se
  out
}
