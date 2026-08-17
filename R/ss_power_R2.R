#' Plan Sample Size to Make the Test of the Squared Multiple Correlation Coefficient Sufficiently Powerful
#'
#' @description
#' Determine the necessary sample size for the omnibus test of the squared
#' multiple correlation coefficient (\eqn{R^2}), or the realized statistical
#' power given a specified sample size, under either fixed or random
#' predictors. The fixed-predictors path uses Cohen's (1988) noncentral
#' \emph{F} formulation; the random-predictors path uses the Lee (1971)
#' two-moment approximation to the sampling distribution of the sample
#' \eqn{R^2} under joint multivariate normality.
#'
#' @param population_R2 Population squared multiple correlation coefficient
#' @param alpha_level Type I error rate
#' @param desired_power Desired degree of statistical power
#' @param p The number of predictor variables
#' @param specified_N The sample size used to calculate power (rather than determine necessary sample size). This is the \emph{total} sample size across all groups or observations.
#' @param cohen_f2 Cohen's (1988) effect size for multiple regression: \code{population_R2}/(1-\code{population_R2})
#' @param null_R2 Value of the null hypothesis that the squared multiple correlation will be evaluated against (this will typically be zero)
#' @param random_predictors Whether the predictor variables are treated as random (\code{TRUE}, the default) or fixed (\code{FALSE}). See Details.
#' @param print_progress If the progress of the iterative procedure is printed to the screen as the iterations are occurring
#' @param \dots Possible additional parameters for internal functions
#'
#' @details
#' Determine the necessary sample size given a particular
#' \code{population_R2}, \code{alpha_level}, \code{p}, and
#' \code{desired_power}. Alternatively, given \code{population_R2},
#' \code{alpha_level}, \code{p}, and \code{specified_N}, the function can
#' be used to determine the statistical power.
#'
#' \strong{Fixed vs.\ random predictors.} The two regression models give
#' \emph{different} sampling distributions for the omnibus \eqn{F}-statistic,
#' and so different power. Under fixed predictors the design matrix is
#' treated as constant in hypothetical replications of the study, and
#' \eqn{F} follows a noncentral \emph{F} with \eqn{p} and \eqn{N-p-1}
#' degrees of freedom and noncentrality \eqn{\lambda = N \cdot f^2},
#' where \eqn{f^2 = \rho^2 / (1 - \rho^2)} (Cohen, 1988). Under random
#' predictors the design matrix is itself a draw from a joint multivariate
#' normal distribution, and the unconditional distribution of the sample
#' \eqn{R^2} is given by Lee (1971); \code{ss_power_R2()} uses Lee's
#' two-moment Patnaik (1949) approximation to that distribution, the same
#' approximation \code{ci_R2()} uses for random-predictor confidence
#' intervals. Gatsonis and Sampson (1989) document the comparison and
#' show that Cohen's fixed-predictor formula tends to over-state power
#' (and so under-state required \eqn{N}) relative to the random model;
#' the discrepancy is modest at moderate to large \eqn{N} but non-trivial
#' for small \eqn{N} with moderate-to-large effects. In the behavioral,
#' educational, and social sciences predictor variables are almost always
#' random, so the default is \code{random_predictors = TRUE}; pass
#' \code{random_predictors = FALSE} for designs in which the predictor
#' variables are fixed by design (for example, planned dosing levels).
#'
#' @return
#' A \code{data.frame} with columns \code{term} and \code{value}.
#' For an \eqn{N} search the rows are \code{necessary_N},
#' \code{actual_power}, \code{noncentral_f_parm} (only meaningful for
#' \code{random_predictors = FALSE}; \code{NA} otherwise), and
#' \code{effect_size} (Cohen's \eqn{f^2}). For a power-at-specified-\code{N}
#' computation the first row is \code{specified_N} instead of
#' \code{necessary_N}.
#'
#' @references
#' Cohen, J. (1988). \emph{Statistical power analysis for the behavioral
#'   sciences} (2nd ed.). Hillsdale, NJ: Lawrence Erlbaum.
#'
#' Gatsonis, C., & Sampson, A. R. (1989). Multiple correlation: Exact
#'   power and sample size calculations. \emph{Psychological Bulletin,
#'   106}(3), 516--524.
#'
#' Kelley, K., & Maxwell, S. E. (2003). Sample size for multiple
#'   regression: Obtaining regression coefficients that are accurate,
#'   not simply significant. \emph{Psychological Methods, 8}(3),
#'   305--321. \doi{10.1037/1082-989X.8.3.305}
#'
#' Kelley, K. (2008). Sample size planning for the squared multiple
#'   correlation coefficient: Accuracy in parameter estimation via narrow
#'   confidence intervals. \emph{Multivariate Behavioral Research, 43}(4),
#'   524--555. \doi{10.1080/00273170802490632}
#'
#' Kelley, K., & Maxwell, S. E. (2008). Sample size planning with
#'   applications to multiple regression: Power and accuracy for omnibus
#'   and targeted effects. In P. Alasuutari, L. Bickman, & J. Brannen
#'   (Eds.), \emph{The Sage handbook of social research methods}
#'   (pp. 166--192). Sage.
#'
#' Lee, Y. S. (1971). Some results on the sampling distribution of the
#'   multiple correlation coefficient. \emph{Journal of the Royal
#'   Statistical Society, Series B, 33}(1), 117--130.
#'
#' Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). \emph{Designing
#'   experiments and analyzing data: A model comparison perspective}
#'   (4th ed.). Routledge. (See Chapter 3 on \eqn{R^2} as a model comparison
#'   effect size.)
#'
#' Maxwell, S. E., Kelley, K., & Rausch, J. R. (2008). Sample size planning
#'   for statistical power and accuracy in parameter estimation.
#'   \emph{Annual Review of Psychology, 59}, 537--563.
#'   \doi{10.1146/annurev.psych.59.103006.093735}
#'
#' Patnaik, P. B. (1949). The non-central \eqn{\chi^2}- and
#'   \emph{F}-distributions and their applications. \emph{Biometrika,
#'   36}(1--2), 202--232. \doi{10.1093/biomet/36.1-2.202}
#'
#' Anderson, S. F., Kelley, K., & Maxwell, S. E. (2017). Sample-size
#'   planning for more accurate statistical power: A method adjusting
#'   sample effect sizes for publication bias and uncertainty.
#'   \emph{Psychological Science, 28}(11), 1547--1562.
#'   \doi{10.1177/0956797617723724}
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @note
#' When determining sample size for a desired degree of power, there will
#' always be a slightly larger degree of actual power. This is the case
#' because the algorithm employed determines sample size until the actual
#' power is no less than the desired power (given sample size is a whole
#' number power will almost certainly not be exactly the specified value).
#' This is the same as other statistical power procedures that return
#' whole numbers for necessary sample size.
#'
#' @seealso
#' \code{\link{ss_aipe_R2}}, \code{\link{ss_power_R2_sensitivity}},
#' \code{\link{ss_power_reg_coef}}, \code{\link{conf_limits_ncf}},
#' \code{\link{ci_R2}}
#'
#' @examples
#' # Random predictors (default; appropriate for most behavioral / social
#' # science applications).
#' ss_power_R2(population_R2 = .5, alpha_level = .05, desired_power = .85, p = 5)
#'
#' # Fixed predictors (Cohen 1988): predictor variables fixed by design.
#' ss_power_R2(population_R2 = .5, alpha_level = .05, desired_power = .85,
#'             p = 5, random_predictors = FALSE)
#'
#' # Effect size input (Cohen's f^2).
#' ss_power_R2(cohen_f2 = 1, alpha_level = .05, desired_power = .85, p = 5)
#'
#' # Realized power at a specified N.
#' ss_power_R2(population_R2 = .5, specified_N = 15, alpha_level = .05,
#'             desired_power = .85, p = 5)
#'
#' @keywords design
#'
#' @seealso \code{\link{design_consequences}} for what a chosen design delivers:
#'   power, the Type S (sign) and Type M (exaggeration) errors of the
#'   significance filter, and the expected confidence interval width.
#'
#' @family sample size for power
#'
#' @export

ss_power_R2 <- function(population_R2 = NULL, alpha_level = .05,
                        desired_power = .85, p, specified_N = NULL,
                        cohen_f2 = NULL, null_R2 = 0,
                        random_predictors = TRUE,
                        print_progress = FALSE, ...) {
  if (alpha_level > 1 || alpha_level < 0) stop("\'alpha_level\' has been specified incorrectly.")
  if (null_R2 > 1 || null_R2 < 0) stop("\'null_R2\' has been specified incorrectly.")
  # Validate desired_power outside the population_R2 branch so that the
  # cohen_f2-only entry point also catches bad inputs (otherwise an invalid
  # desired_power can drive the iterative search into an infinite loop).
  if (is.null(specified_N) && (desired_power > 1 || desired_power < 0)) {
    stop("\'desired_power\' has been specified incorrectly.")
  }

  if (!is.null(population_R2)) {
    if (!is.null(cohen_f2)) stop("Since you specified \'population_R2\', do not specify \'cohen_f2\'.")
    if (population_R2 > 1 || population_R2 < 0) stop("\'population_R2\' has been specified incorrectly.")
    f_Squared <- (population_R2 - null_R2) / (1 - population_R2)
  }

  if (!is.null(cohen_f2)) {
    f_Squared <- cohen_f2
    # The random-predictors path needs an explicit population_R2; back it
    # out from cohen_f2 if necessary. The inversion of
    #   cohen_f2 = (R^2 - null_R2) / (1 - R^2)
    # is
    #   R^2 = (cohen_f2 + null_R2) / (1 + cohen_f2)
    # which reduces to cohen_f2 / (1 + cohen_f2) when null_R2 = 0 (the default).
    if (is.null(population_R2)) {
      population_R2 <- (cohen_f2 + null_R2) / (1 + cohen_f2)
    }
  }

  # Power at a given N under the chosen sampling model. If the alternative
  # is at or below the null (no real effect in the hypothesized direction),
  # the right-tail F-test rejects with probability alpha_level regardless
  # of the predictor model.
  .power_at_N <- function(N_now) {
    df_2 <- N_now - p - 1
    if (df_2 < 1) return(0)
    if (!is.null(population_R2) && population_R2 <= null_R2) return(alpha_level)
    f_crit <- qf(1 - alpha_level, df1 = p, df2 = df_2,
                 lower.tail = TRUE, log.p = FALSE)
    if (random_predictors) {
      # R^2 critical corresponding to f_crit under H_0.
      R2_crit <- (p * f_crit) / (p * f_crit + df_2)
      1 - .lee_random_R2_cdf(R2_obs = R2_crit, rho2 = population_R2,
                             N = N_now, p = p)
    } else {
      1 - pf(f_crit, df1 = p, df2 = df_2,
             ncp = N_now * f_Squared,
             lower.tail = TRUE, log.p = FALSE)
    }
  }

  # In the sample-size search, power at every N is capped at alpha_level when
  # the alternative is at or below the null, so desired_power is unreachable and
  # the 'while (Dif > 0)' loop below would never terminate. Stop cleanly, as the
  # sibling ss_power_r() does for rho at or below its null.
  if (is.null(specified_N) && !is.null(population_R2) && population_R2 <= null_R2) {
    stop("With 'population_R2' at or below 'null_R2' there is no effect in the hypothesized direction, so power cannot exceed 'alpha_level' and no finite sample size attains 'desired_power'. Increase 'population_R2' or decrease 'null_R2'.", call. = FALSE)
  }

  if (is.null(specified_N)) {
    # Evaluate the smallest admissible size (df_2 >= 1) before incrementing, so a
    # target the minimum already attains returns that minimum, not one above it.
    N_i <- p + 1 + 1
    repeat {
      Actual_Power <- .power_at_N(N_i)

      if (print_progress) {
        cat(c(Current.Power = Actual_Power,
              Current.NC.F.Parm = ifelse(random_predictors, NA, N_i * f_Squared),
              Current.N = N_i), "\n")
      }
      if (Actual_Power >= desired_power) break
      N_i <- N_i + 1
      if (N_i > 1e7) stop("Failed to reach 'desired_power' within a reasonable sample size.", call. = FALSE)
    }
    VALUE <- data.frame(
      term  = c("necessary_N", "actual_power", "noncentral_f_parm", "effect_size"),
      value = c(N_i, Actual_Power,
                ifelse(random_predictors, NA_real_, N_i * f_Squared),
                f_Squared)
    )
  } else {
    # A total N below p + 2 leaves no residual degrees of freedom; a fractional
    # N is not a sample size. Reject both rather than returning 0, NaN, or a
    # power table for a size the design cannot have.
    specified_N <- .check_whole_n(specified_N, "specified_N", p + 2L)
    actual_power_specified_N <- .power_at_N(specified_N)
    VALUE <- data.frame(
      term  = c("specified_N", "actual_power", "noncentral_f_parm", "effect_size"),
      value = c(specified_N, actual_power_specified_N,
                ifelse(random_predictors, NA_real_, specified_N * f_Squared),
                f_Squared)
    )
  }
  class(VALUE) <- c("dmar_ss_power", "dmar_tbl", "data.frame")
  return(VALUE)
}
