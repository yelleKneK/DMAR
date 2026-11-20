#' Sample Size Planning for a Single Regression Coefficient (AIPE)
#'
#' @description
#' Computes the necessary sample size for the confidence interval on a
#' targeted regression coefficient \eqn{\beta_j}, the quantity
#' \code{true_reg_coef_j} stands for, in a multiple regression with
#' \emph{p} predictors to be no wider than a user-specified value. The
#' coefficient is on the scale \code{true_sd_Y} and \code{true_sd_X}
#' describe: with both 1, the default, its population value is the
#' standardized coefficient. This is the accuracy in parameter estimation
#' (AIPE) framework of Kelley and Maxwell (2003), targeted at a specific
#' coefficient rather than at the omnibus \eqn{R^2}. The default
#' (\code{noncentral = FALSE}) plans for the central \emph{t} interval for
#' the coefficient in the units of the data, the interval
#' \code{\link{ci_rc}} reports; \code{noncentral = TRUE} plans for the
#' noncentral \emph{t} interval for the standardized coefficient (Kelley,
#' 2007a), the interval \code{\link{ci_src}} reports, with unit standard
#' deviations. The central \emph{t} interval for the standardized
#' coefficient, the interval of Kelley and Maxwell (2003, Equation 1), is a
#' third interval, and \code{\link{ss_aipe_src}} with
#' \code{noncentral = FALSE} plans for it; see "The interval planned for"
#' in Details. Supplying
#' \code{assurance} returns the smallest \emph{N}, not below the plan
#' without it (the plan for the width at the expected values of the sample
#' squared multiple correlations), at which the realized interval is no
#' wider than \code{width} with at least that probability, computed exactly
#' for the central interval; for the noncentral interval it returns by
#' default the published approximation that \pkg{MBESS} uses, which does
#' not compute the probability, and with
#' \code{assurance_method = "simulation"} the smallest \emph{N} at which a
#' simulated estimate of the probability reaches the request.
#'
#' @param true_R2_Y_X Population value of \eqn{\rho^2_{Y \cdot X_1, \ldots, X_p}},
#'   the squared multiple correlation of the outcome \emph{Y} with all
#'   \emph{p} predictors.
#' @param true_R2_j_X_without_j Population value of
#'   \eqn{\rho^2_{X_j \cdot X_{-j}}}, the squared multiple correlation when
#'   the \emph{j}th predictor is regressed on the remaining \eqn{p - 1}
#'   predictors. Quantifies the multicollinearity faced by the targeted
#'   coefficient.
#' @param p The number of predictor variables.
#' @param true_reg_coef_j Population value \eqn{\beta_j} of the
#'   regression coefficient for the \emph{j}th predictor, the predictor of
#'   interest, on the scale \code{true_sd_Y} and \code{true_sd_X} describe:
#'   with both 1 it is the standardized coefficient, and otherwise the
#'   unstandardized one, the standardized coefficient times
#'   \code{true_sd_Y / true_sd_X} (Kelley & Maxwell, 2003, p. 316);
#'   \eqn{b_j} is its sample estimate. On the covariance structure path it
#'   may be left out; when it is given it must equal the coefficient the
#'   correlations imply on that scale.
#' @param width Desired (full) width of the two-sided confidence interval
#'   on \eqn{\beta_j}.
#' @param which_width Which portion of the confidence interval
#'   \code{width} refers to. Only \code{"Full"} (or an abbreviation of
#'   it, such as \code{"F"}) is currently implemented; \code{"Lower"},
#'   \code{"Upper"}, and their abbreviations are refused.
#' @param true_sd_Y Population standard deviation of \emph{Y}, which with
#'   \code{true_sd_X} sets the scale of the coefficient and of \code{width}.
#' @param true_sd_X Population standard deviation of the \emph{j}th predictor.
#' @param true_cor_XX Population correlation matrix for the \emph{p} predictor
#'   variables. If supplied with \code{true_cor_YX}, \code{true_R2_Y_X} and
#'   \code{true_R2_j_X_without_j} are derived from the covariance structure.
#'   With one predictor it is \code{matrix(1)} (or the number 1), and
#'   \code{p} and \code{which_predictor} are then 1 when not given.
#' @param true_cor_YX Length-\emph{p} vector of population correlations between
#'   \emph{Y} and the \emph{p} predictors.
#' @param which_predictor Which of the \emph{p} predictors is the targeted
#'   coefficient.
#' @param noncentral If \code{TRUE}, plans for the noncentral \emph{t}
#'   interval for the standardized coefficient (Kelley, 2007a), the interval
#'   of \code{\link{ci_src}}. If \code{FALSE} (the default), plans for the
#'   central \emph{t} interval for the coefficient in the units of the
#'   data, the interval of \code{\link{ci_rc}}, even when both standard
#'   deviations are 1; the central interval for the standardized
#'   coefficient is planned by \code{\link{ss_aipe_src}} with
#'   \code{noncentral = FALSE}. The noncentral path requires
#'   \code{true_sd_Y = true_sd_X = 1} (\emph{i.e.}, a standardized solution).
#' @param alpha_lower Type I error rate for the lower confidence limit,
#'   at least 1e-9 and at most .5. A tail area below 1e-9 is refused: the
#'   noncentral \emph{t} limits lose accuracy there, and the central plan
#'   keeps the same floor.
#' @param alpha_upper Type I error rate for the upper confidence limit,
#'   at least 1e-9 and at most .5.
#' @param conf_level Confidence level (i.e., \eqn{1 - \alpha}, where
#'   \eqn{\alpha} is the Type I error rate). Default \code{0.95}. Mutually
#'   exclusive with \code{alpha_lower} and \code{alpha_upper}. A level
#'   above .999999998, whose tail areas are below 1e-9, is refused.
#' @param assurance Optional probability with which the realized confidence
#'   interval is to be no wider than \code{width}, a single number above .5
#'   and below 1. When \code{NULL} (the default), the plan is the smallest
#'   \emph{N} at which the interval is no wider than \code{width} when the
#'   two sample squared multiple correlations take their expected values
#'   (see "The interval planned for" in Details): the width at those
#'   expected values, not the mean of the width over samples.
#' @param assurance_method How the probability that the realized interval
#'   is no wider than \code{width} is found, used only with
#'   \code{assurance} (see "The assurance step" in Details).
#'   \code{"analytic"} (the default) computes without simulation: for the
#'   central interval (\code{noncentral = FALSE}) the probability is exact,
#'   and the function returns the smallest sample size, at or above the plan
#'   without \code{assurance}, at which it is at least \code{assurance},
#'   reported as \code{assurance_attained}; for the noncentral interval it
#'   returns the published approximation that \pkg{MBESS} uses, which
#'   does not compute the probability (\code{assurance_attained} is
#'   \code{NA}) and can attain less than \code{assurance}.
#'   \code{"simulation"}, offered for the noncentral interval only, returns
#'   the smallest sample size, at or above the plan without
#'   \code{assurance}, at which the probability estimated from \code{G}
#'   simulated samples is at least \code{assurance}, and reports that
#'   estimate. With \code{noncentral = FALSE} \code{"simulation"} is refused
#'   with a condition of class \code{dmar_argument_invalid}, since a
#'   simulated estimate of an exact probability would only add Monte Carlo
#'   error.
#' @param G Number of simulated samples at each sample size the noncentral
#'   assurance step tries, and in its final estimate of the attained
#'   probability; a whole number of at least 100 (default 20000). Used and
#'   checked only with \code{noncentral = TRUE}, \code{assurance}, and
#'   \code{assurance_method = "simulation"}.
#' @param seed Optional single number. When given, the noncentral assurance
#'   step is reproducible, and the random number stream of the calling
#'   session is left as it was. With \code{NULL} (the default) the draws
#'   come from the session's current stream. Used only where \code{G} is.
#'
#' @details
#' \strong{Calling conventions.} The function offers several mutually
#' exclusive ways to supply the population information needed to plan
#' \eqn{N}; the user picks the one most aligned with their available
#' planning values. Specify exactly one of:
#' \itemize{
#'   \item \emph{Covariance structure path.} Supply \code{true_cor_XX} and
#'     \code{true_cor_YX} along with \code{which_predictor}. The function
#'     derives \eqn{\rho^2_{Y \cdot X}} and
#'     \eqn{\rho^2_{X_j \cdot X_{-j}}} (the values of \code{true_R2_Y_X}
#'     and \code{true_R2_j_X_without_j}) from the covariance structure, and
#'     the population coefficient \eqn{\beta_j}, the \emph{j}th element of
#'     \eqn{R_{XX}^{-1} r_{XY}} (\eqn{R_{XX}} is \code{true_cor_XX} and
#'     \eqn{r_{XY}} is \code{true_cor_YX}) times
#'     \code{true_sd_Y / true_sd_X}, against which a supplied
#'     \code{true_reg_coef_j} is checked.
#'   \item \emph{Squared multiple correlations path.} Supply
#'     \code{true_R2_Y_X}, \code{true_R2_j_X_without_j}, \code{p}, and
#'     \code{true_reg_coef_j} directly when the user already has these planning values
#'     from prior research and does not need to specify the full
#'     covariance structure.
#'   \item \emph{Standardized solution.} For the \code{noncentral = TRUE}
#'     path, set \code{true_sd_Y = true_sd_X = 1}; the result returns the
#'     standardized regression coefficient sample size.
#' }
#'
#' \strong{The interval planned for.} With \code{noncentral = FALSE} the
#' interval is the central \emph{t} interval for the coefficient in the
#' units of the data, \eqn{b_j} plus the \emph{t} quantiles (on
#' \eqn{N - p - 1} degrees of freedom) at \code{alpha_lower} and
#' \code{1 - alpha_upper} times its standard error, the interval
#' \code{\link{ci_rc}} and \code{\link{ci_reg_coef}} report from
#' unstandardized data, and this is so whatever \code{true_sd_Y} and
#' \code{true_sd_X} are: with both 1 the coefficient is on the
#' standardized scale, but the interval planned for is still the one
#' computed in the units of the data, whose standard error comes from the
#' residual sums of squares of \emph{Y} and of \eqn{X_j}. The plan without
#' \code{assurance} is the smallest \emph{N} at which that interval,
#' computed at the expected values of the two sample squared multiple
#' correlations, is no wider than \code{width} (Kelley & Maxwell, 2003,
#' Equations 1 and 2, with the \emph{t} quantile); the two tails exchanged
#' give the same plan. The plan is the smallest such \emph{N} of at least
#' \eqn{p + 2}, the smallest sample with a residual degree of freedom, and
#' the expected squared multiple correlation of \eqn{X_j} with the other
#' predictors is that of a regression on those \eqn{p - 1} predictors (0
#' with one predictor; see "Comparison with MBESS"). Kelley and Maxwell's
#' Equation 1 is the central interval for the standardized coefficient,
#' whose standard error is formed from the two sample squared multiple
#' correlations; its plan without \code{assurance} is this one at unit
#' standard deviations, but the probability that it is no wider than
#' \code{width} is not the probability the assurance step below computes,
#' and \code{\link{ss_aipe_src}} with \code{noncentral = FALSE} plans its
#' assurance.
#' With \code{noncentral = TRUE} the interval is the noncentral \emph{t}
#' interval for the standardized coefficient of \code{\link{ci_src}}: the
#' limits of the noncentrality parameter at the observed \emph{t} statistic,
#' times the standard error
#' \eqn{\sqrt{(1 - R^2)/[(1 - R^2_{X_j \cdot X_{-j}})(N - p - 1)]}}, and the
#' plan without \code{assurance} is the smallest \emph{N}, of at least
#' \eqn{p + 2}, at which it is no wider than \code{width} at the expected
#' squared multiple correlations and \eqn{\beta_j}. That standard error
#' treats the sample standard deviations as constants, which for random
#' predictors is an approximation (Kelley & Maxwell, 2003, pp. 315--316).
#' The central path is fast and
#' adequate at moderate to large \emph{N}; the noncentral path is the one
#' to use when the analysis will report the noncentral interval for the
#' standardized coefficient. A tail of 0, or two of .5, is refused, and so
#' is a tail area below 1e-9 or a \code{conf_level} above .999999998.
#'
#' \strong{The assurance step.} A plan without \code{assurance} makes the
#' interval as narrow as the target when the sample values equal their
#' expectations; in roughly half of all samples it is wider. With
#' \code{assurance} the function plans so that the realized interval is no
#' wider than \code{width} with at least that probability, for random
#' predictors from a multivariate normal population.
#' \itemize{
#'   \item \emph{The central interval} (\code{noncentral = FALSE}). The
#'     only method is \code{"analytic"}: the function computes \eqn{P_N},
#'     the probability that the realized interval is no wider than
#'     \code{width} at sample size \emph{N}, exactly, and returns the
#'     smallest \emph{N}, at or above the plan without \code{assurance}, at
#'     which \eqn{P_N} is at least \code{assurance}, found by doubling a
#'     step upward and then bisecting. \code{assurance_attained} is
#'     \eqn{P_N} at the returned \emph{N}. The squared standard
#'     error of \eqn{b_j} is the residual sum of squares of \emph{Y} on all
#'     \emph{p} predictors over \eqn{N - p - 1}, divided by the residual sum
#'     of squares of \eqn{X_j} on the other \eqn{p - 1}. Under multivariate
#'     normality the first, over \eqn{\sigma^2_Y (1 - \rho^2_{Y \cdot X})},
#'     is chi square on \eqn{N - p - 1} degrees of freedom whatever the
#'     predictors are, so it is independent of the second, which over
#'     \eqn{\sigma^2_{X_j} (1 - \rho^2_{X_j \cdot X_{-j}})} is chi square
#'     on \eqn{N - p} degrees of freedom. The squared standard error is then
#'     \eqn{K F / (N - p)}, with
#'     \eqn{K = \sigma^2_Y (1 - \rho^2_{Y \cdot X}) /
#'     [\sigma^2_{X_j} (1 - \rho^2_{X_j \cdot X_{-j}})]} and \emph{F} an
#'     \emph{F} variate on \eqn{N - p - 1} and \eqn{N - p} degrees of
#'     freedom, and with \eqn{q} the sum of the two \emph{t} quantiles of
#'     the interval,
#'     \deqn{P_N = \Pr(F \le w^2 (N - p) / (q^2 K)),}{P_N = Pr(F <= w^2 (N - p) / (q^2 K)),}
#'     where \eqn{w} is \code{width}. (With fixed predictors the second sum
#'     of squares is a constant and the same event is a chi square
#'     probability; the plan is for random predictors, the model
#'     \code{\link{ss_aipe_reg_coef_sensitivity}} simulates.) The units of
#'     \emph{Y} and \eqn{X_j} cancel, so changing them, with \code{width}
#'     changed to match, leaves the plan as it was, and the sign of
#'     \eqn{\beta_j} does not enter.
#'   \item \emph{The noncentral interval} (\code{noncentral = TRUE}). The
#'     width depends on the observed \emph{t} statistic and on both sample
#'     squared multiple correlations, and no closed form for its distribution
#'     is used. With \code{assurance_method = "analytic"}, the default, the
#'     plan is the approximation that \pkg{MBESS} uses (Kelley, 2007a,
#'     2007b): the plan without \code{assurance} made again with a
#'     lower confidence limit (at level \code{assurance}) of the squared
#'     multiple correlation in place of its population value, and not below
#'     the first. It does not compute the probability it attains, so
#'     \code{assurance_attained} is \code{NA}, and the printed footer says
#'     so. Kelley and Maxwell (2008, pp. 186--187) wrote that no satisfactory
#'     method then existed for the standardized coefficient and described
#'     finding, by simulation, the smallest sample size that attains the
#'     requested degree of certainty, which is what
#'     \code{assurance_method = "simulation"} does: \eqn{P_N} is estimated as
#'     the proportion of \code{G} samples, drawn from the exact sampling
#'     distribution of the sums of
#'     squares and products (Wishart, through the entries of its triangular
#'     factor that the statistics depend on), whose interval is no wider than
#'     \code{width}. The same \code{G} samples serve every sample size the
#'     search tries (common random numbers: each chi square entry is the
#'     quantile, at the degrees of freedom of that \emph{N}, of one uniform
#'     draw), so the estimate changes smoothly with \emph{N}. The search
#'     assumes that it does not decrease as \emph{N} grows, which held in
#'     every set checked although a single simulated interval can be wider
#'     at \eqn{N + 1} than at \emph{N} (see below), and under that
#'     assumption the step returns the smallest \emph{N} at which the
#'     estimate is at least \code{assurance}. The
#'     distribution of the width depends on the population through
#'     \emph{p}, \eqn{\beta_j}, the two squared multiple correlations, and
#'     the correlation of \emph{Y} with \eqn{X_j}. On the covariance
#'     structure path that correlation is given; on the squared multiple
#'     correlations path it is not, and the step uses the smallest estimated
#'     probability over 21 equally spaced values of it, from the smallest to
#'     the largest that the three values allow. That is an approximate
#'     minimum over the range: between two of the values the probability
#'     can be a little lower than at either, so for a correlation there the
#'     plan can fall a little short of the request. (With \eqn{p = 2} only
#'     the two ends of the range are possible, and the step uses both.)
#'     \code{assurance_attained} is the estimate at the returned
#'     \emph{N}, so it is at least \code{assurance}, and its Monte Carlo
#'     standard error is the attribute \code{assurance_mc_se}. The exact
#'     \eqn{P_N} differs from the estimate by that error, so the returned
#'     \emph{N} can be a case or two from the one the exact probability
#'     would give, and \eqn{P_N} there can be a little below
#'     \code{assurance}; give \code{seed} for a plan that can be
#'     reproduced, and a larger \code{G} for a more precise one. Reversing
#'     the sign of \eqn{\beta_j} (with \code{alpha_lower} and
#'     \code{alpha_upper} exchanged when they differ) gives the same plan.
#' }
#' For both intervals \eqn{P_N} did not decrease from one \emph{N} to the next
#' at or above the plan without \code{assurance}: for the central interval
#' in 1,530 planning sets (\emph{p} from 1 to 15, squared multiple
#' correlations from 0 to .9, widths from .05 to 2, and tails of .025 and
#' .025, .005 and .005, and .01 and .04) over every \emph{N} from that plan
#' to three times it, and for the noncentral interval in five sets by
#' simulation. With the
#' samples held fixed across \emph{N}, the estimate did not decrease from
#' one \emph{N} to the next in seven sets (\emph{p} from 1 to 10), and of
#' 95 million comparisons of a simulated interval at \emph{N} and at
#' \eqn{N + 1}, 3208, all at \emph{N} below 60, were wider at
#' \eqn{N + 1}. Below the plan without \code{assurance} \eqn{P_N} can
#' decrease, which is why the search starts there.
#'
#' The approximations can attain more or less than the request. For two
#' predictors correlated .3 with correlations of .4 and .3 with \emph{Y},
#' a width of .3, and an assurance of .90, the chi square approximation of
#' Kelley and Maxwell (2003, Equation 3; 2008, Equation 37) for the central
#' interval, the plan without \code{assurance} with the variance of
#' \eqn{b_j} multiplied by the \code{assurance} quantile of a chi square on
#' \eqn{N - 1} degrees of freedom over \eqn{N - p - 1}, is 181, at which
#' \eqn{P_N} is .852, and the exact plan is 187 (\eqn{P_N} = .902); the
#' noncentral approximation, the default plan for that interval, is 179, at
#' which the simulated probability is .948, and
#' \code{assurance_method = "simulation"} returns 176 (.908 with
#' \code{seed = 1}). For the five exchangeable predictors of Kelley and
#' Maxwell (2003, pp. 312--313), correlated .4 with one another and .3 with
#' \emph{Y}, at a width of .2 and an assurance of .80, the chi square
#' approximation is 486 (\eqn{P_N} = .765) and the exact plan 491 (.800).
#' Each was checked against
#' the proportion of intervals no wider than the target among samples of
#' raw multivariate normal data. These probabilities, and those of the
#' first example, are for the central interval in the units of the data,
#' the interval this function plans for. Kelley and Maxwell's example is of
#' the standardized coefficient, with the central interval for it (their
#' Equation 1), which is computed from standardized data with the sample
#' standard deviations treated as constants; its width varies less, and
#' it is no wider than .2 with probability .877 at 486 and .911 at 491
#' (their own simulation, p. 317, reported .816 at their plan of 485, where
#' raw multivariate normal data give .871 with the \emph{t} quantile of
#' their Equation 1, or .887 with the normal quantile of their plan). For
#' that interval \code{\link{ss_aipe_src}} with
#' \code{noncentral = FALSE}, \code{assurance_method = "simulation"}, and
#' \code{seed = 1} plans 478 (.807, and .797
#' at 477), and in the first example 171 (.914, and .895 at 170), where the
#' plans above reach .986 at 181 and .997 at 187; its plan is simulated, and
#' another seed can give one a case away (170 for one of seeds 1 to 30 in
#' the first example).
#'
#' \strong{Comparison with MBESS.} The plans without \code{assurance} of
#' \code{ss.aipe.reg.coef()} in \pkg{MBESS} differ from these in two ways.
#' They take the expected squared multiple correlation of \eqn{X_j} with
#' the other predictors to be that of a regression on \emph{p} predictors,
#' where that regression has \eqn{p - 1}, which overstates it and so the
#' width. And their search goes up one case at a time from a start, the
#' larger of \eqn{2p + 1} and ten below the normal approximation for the
#' central interval and six below the central plan for the noncentral one,
#' so it returns the start when a smaller sample size already meets the
#' width: for \code{true_R2_Y_X = .8}, \code{true_R2_j_X_without_j = 0},
#' \code{p = 10}, and \code{width = 2}, 21 where 14 does. In 1,530
#' planning sets the two differences make the central plan here one case
#' smaller in 1,329 and 2 to 13 cases smaller in 104, and in 810 sets the
#' noncentral plan one smaller in 596 and 2 to 10 smaller in 177. The
#' approximations start from these plans (Kelley and Maxwell's example
#' above is 487 in \pkg{MBESS} and 486 here). With unequal tails
#' \code{ss.aipe.reg.coef()} also takes both \emph{t} quantiles of its
#' central plans at the average of the two tails, so it plans for an
#' interval narrower than the one \code{\link{ci_rc}} reports; the central
#' plans here take the quantile at each tail.
#'
#' @return
#' A \code{data.frame} with columns \code{term} and \code{value}:
#' \describe{
#'   \item{\code{necessary_N}}{the necessary total sample size.}
#'   \item{\code{assurance}}{the requested probability; present only when
#'     \code{assurance} was given.}
#'   \item{\code{assurance_attained}}{the probability, at
#'     \code{necessary_N}, that the realized interval is no wider than
#'     \code{width}; present only when \code{assurance} was given. Exact for
#'     the central interval; for the noncentral interval \code{NA} with
#'     \code{assurance_method = "analytic"}, whose approximation does not
#'     compute it, and with \code{"simulation"} a Monte Carlo estimate,
#'     whose standard error is the attribute \code{assurance_mc_se}.}
#'   \item{\code{G}}{the number of simulated samples; present only when
#'     the noncentral assurance step simulated.}
#' }
#' With \code{assurance} the printed footer names the assurance method and
#' says how \code{assurance_attained} was found, or that the approximation
#' does not compute it; the attributes \code{assurance_method}
#' (\code{"analytic"} or \code{"simulation"}) and
#' \code{assurance_probability} (\code{"exact"}, \code{"simulation"}, or
#' \code{"not computed"}) record the same.
#'
#' @references
#' Kelley, K. (2007a). Confidence intervals for standardized effect sizes:
#'   Theory, application, and implementation. \emph{Journal of Statistical
#'   Software, 20}(8), 1--24. \doi{10.18637/jss.v020.i08}
#'
#' Kelley, K. (2007b). Methods for the behavioral, educational, and
#'   social sciences: An R package. \emph{Behavior Research Methods,
#'   39}(4), 979--984. \doi{10.3758/BF03192993}
#'
#' Kelley, K., & Maxwell, S. E. (2003). Sample size for multiple
#'   regression: Obtaining regression coefficients that are accurate,
#'   not simply significant. \emph{Psychological Methods, 8}(3),
#'   305--321. \doi{10.1037/1082-989X.8.3.305}
#'
#' Kelley, K., & Maxwell, S. E. (2008). Sample size planning with
#'   applications to multiple regression: Power and accuracy for omnibus
#'   and targeted effects. In P. Alasuutari, L. Bickman, & J. Brannen
#'   (Eds.), \emph{The Sage handbook of social research methods}
#'   (pp. 166--192). Sage.
#'
#' Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027).
#'   \emph{Designing experiments and analyzing data: A model comparison
#'   perspective} (4th ed.). Routledge. (See Chapter 4 on individual
#'   comparisons of means and Chapter 6 on trend analysis.)
#'
#' Maxwell, S. E., Kelley, K., & Rausch, J. R. (2008). Sample size planning
#'   for statistical power and accuracy in parameter estimation.
#'   \emph{Annual Review of Psychology, 59}, 537--563.
#'   \doi{10.1146/annurev.psych.59.103006.093735}
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso
#' \code{\link{ss_aipe_reg_coef_sensitivity}}, \code{\link{ci_nc_t}}
#'
#' @examples
#' # 1. Covariance structure path: supply the population correlation
#' #    matrix and the population YX cross-correlations. Five predictors
#' #    in an exchangeable structure (all pairwise correlations 0.5,
#' #    all Y-X correlations 0.3).
#' true_cor_YX <- c(.3, .3, .3, .3, .3)
#' true_cor_XX <- rbind(c(1, .5, .5, .5, .5), c(.5, 1, .5, .5, .5),
#'                 c(.5, .5, 1, .5, .5), c(.5, .5, .5, 1, .5),
#'                 c(.5, .5, .5, .5, 1))
#'
#' # Closed-form (central t) planning, targeting the first predictor's
#' # standardized coefficient.
#' ss_aipe_reg_coef(width = .10, which_width = "Full",
#'                  true_sd_Y = 1, true_sd_X = 1,
#'                  true_cor_XX = true_cor_XX, true_cor_YX = true_cor_YX,
#'                  which_predictor = 1, noncentral = FALSE,
#'                  conf_level = .95)
#'
#' # With an assurance of .85: the smallest N at which the realized
#' # interval is no wider than 0.10 with probability at least .85. The
#' # last row is that probability at the returned N, computed exactly.
#' ss_aipe_reg_coef(width = .10, which_width = "Full",
#'                  true_sd_Y = 1, true_sd_X = 1,
#'                  true_cor_XX = true_cor_XX, true_cor_YX = true_cor_YX,
#'                  which_predictor = 1, noncentral = FALSE,
#'                  conf_level = .95, assurance = .85)
#'
#' # The noncentral t interval for the standardized coefficient. Required
#' # N differs at small to moderate samples.
#' ss_aipe_reg_coef(width = .10, which_width = "Full",
#'                  true_sd_Y = 1, true_sd_X = 1,
#'                  true_cor_XX = true_cor_XX, true_cor_YX = true_cor_YX,
#'                  which_predictor = 1, noncentral = TRUE,
#'                  conf_level = .95)
#'
#' # With assurance its default plan is the approximation that MBESS uses,
#' # which does not compute the probability it attains; by simulation, the
#' # smallest N at which the estimated probability is at least .85. G = 2000
#' # keeps the example quick; a reported plan deserves the default G = 20000.
#' ss_aipe_reg_coef(width = .10, true_cor_XX = true_cor_XX,
#'                  true_cor_YX = true_cor_YX, which_predictor = 1,
#'                  noncentral = TRUE, assurance = .85)
#' plan <- ss_aipe_reg_coef(width = .10, true_cor_XX = true_cor_XX,
#'                          true_cor_YX = true_cor_YX, which_predictor = 1,
#'                          noncentral = TRUE, assurance = .85,
#'                          assurance_method = "simulation",
#'                          G = 2000, seed = 113)
#' plan
#' attr(plan, "assurance_mc_se")
#'
#' # 2. Squared multiple correlations path: when the user has planning
#' #    values for true_R2_Y_X and true_R2_j_X_without_j directly (e.g.,
#' #    from a prior power analysis), without specifying the full
#' #    covariance structure. true_reg_coef_j is required on this path.
#' ss_aipe_reg_coef(true_R2_Y_X = 0.30, true_R2_j_X_without_j = 0.20,
#'                  p = 5, true_reg_coef_j = 0.25,
#'                  width = .15, which_width = "Full",
#'                  true_sd_Y = 1, true_sd_X = 1,
#'                  noncentral = FALSE, conf_level = .95)
#'
#' # 3. One predictor: its correlation matrix is matrix(1), and p and
#' #    which_predictor follow from it.
#' ss_aipe_reg_coef(width = .20, true_cor_XX = matrix(1), true_cor_YX = .4)
#'
#' @keywords design
#'
#' @seealso \code{\link{design_consequences}} for what a chosen design delivers:
#'   power, the Type S (sign) and Type M (exaggeration) errors of the
#'   significance filter, and the expected confidence interval width.
#'
#' @export


ss_aipe_reg_coef <- function(true_R2_Y_X = NULL, true_R2_j_X_without_j = NULL, p = NULL,
                             true_reg_coef_j = NULL, width, which_width = "Full", true_sd_Y = 1, true_sd_X = 1,
                             true_cor_XX = NULL, true_cor_YX = NULL, which_predictor = NULL, noncentral = FALSE,
                             alpha_lower = NULL, alpha_upper = NULL, conf_level = 0.95,
                             assurance = NULL,
                             assurance_method = c("analytic", "simulation"),
                             G = 20000, seed = NULL) {
  .dmar_check_conf_level(conf_level, null_ok = TRUE)
  if (!missing(width)) .dmar_check_width(width)
  .dmar_check_assurance_single(assurance)
  .dmar_check_numbers("ss_aipe_reg_coef", c("true_sd_Y", "true_sd_X", "which_predictor"))
  # "simulation" is offered for the noncentral interval only; the
  # probability for the central one is exact (R/ss_aipe_assurance_internals.R).
  nc <- isTRUE(noncentral == TRUE)
  assurance_method <- .ss_aipe_assurance_method(
    assurance_method, "ss_aipe_reg_coef", simulation = nc,
    where = " with noncentral = FALSE")
  # The assurance step: for the central interval the exact one, and for the
  # noncentral interval the approximation of DMAR 1.0.x ("analytic") or the
  # simulated one.
  step <- if (!nc) "exact" else if (assurance_method == "simulation")
    "simulation" else "approximation"
  # The simulated step estimates a probability from G samples at each
  # sample size it tries. G is checked, and the seed set, only for that
  # step.
  if (!is.null(assurance) && step == "simulation") {
    G <- .dmar_check_count(G, "G", min = 100L, reason = paste0(
      "The assurance step of the noncentral plan estimates a probability ",
      "from 'G' simulated samples at each sample size it tries; with 100 ",
      "the standard error of that estimate is at most .05, and the default ",
      "of 20000 brings it to .004 or less."))
    .dmar_local_seed(seed)
  }
  # The noncentral t limits come from ci_nc_t(), which is exact at every
  # noncentrality (R/nc_t_internals.R), so no accuracy warning is raised or
  # counted here.
  .ss_aipe_reg_coef_plan(
    true_R2_Y_X = true_R2_Y_X, true_R2_j_X_without_j = true_R2_j_X_without_j,
    p = p, true_reg_coef_j = true_reg_coef_j, width = width,
    which_width = which_width, true_sd_Y = true_sd_Y, true_sd_X = true_sd_X,
    true_cor_XX = true_cor_XX, true_cor_YX = true_cor_YX,
    which_predictor = which_predictor, noncentral = noncentral,
    alpha_lower = alpha_lower, alpha_upper = alpha_upper,
    conf_level = conf_level, assurance = assurance, step = step, G = G)
}

# .ss_aipe_reg_coef_plan(): the plan of ss_aipe_reg_coef() after the checks
# above, with the same arguments less the seed (set by the caller) and with
# `step`, the assurance step, in place of assurance_method: "exact" (the
# exact probability of the central interval), "approximation" (the
# published approximation, which does not compute the probability: the chi
# square step of Kelley and Maxwell, 2003, for the central interval, and
# that of DMAR 1.0.x for the noncentral one), or "simulation" (the
# noncentral interval only). ss_aipe_reg_coef() gives "exact" for the
# central interval; ss_aipe_src() calls this function directly, for the
# central interval with "approximation". Not exported.
.ss_aipe_reg_coef_plan <- function(true_R2_Y_X, true_R2_j_X_without_j, p,
                                   true_reg_coef_j, width, which_width,
                                   true_sd_Y, true_sd_X, true_cor_XX,
                                   true_cor_YX, which_predictor, noncentral,
                                   alpha_lower, alpha_upper, conf_level,
                                   assurance, step, G) {

  # A single predictor: a 1 x 1 'true_cor_XX' (or a single number) is
  # enough to set p = 1 on the covariance structure path.
  if (is.null(p) && !is.null(true_cor_XX) && length(true_cor_XX) == 1L) p <- 1
  if (!is.null(p)) {
    if (p == 1) {
      if (!is.null(true_cor_XX)) true_cor_XX <- as.matrix(true_cor_XX)
      if (!is.null(true_cor_YX)) true_cor_YX <- as.matrix(true_cor_YX)
      if (!is.null(true_R2_j_X_without_j)) {
        if (true_R2_j_X_without_j != 0) .dmar_refuse_argument(
          "With p = 1 there is no other predictor, so 'true_R2_j_X_without_j' must be 0.",
          arg = "true_R2_j_X_without_j")
      }
      if (!is.null(true_cor_XX)) {
        # dim() is a vector of two, so the comparison is of the whole
        # dimension, not elementwise inside if().
        if (!identical(dim(true_cor_XX), c(1L, 1L))) .dmar_refuse_argument(
          "With p = 1, 'true_cor_XX' must be a 1 x 1 matrix (the correlation of the one predictor with itself).",
          arg = "true_cor_XX")
        if (is.null(true_cor_YX) || !identical(dim(true_cor_YX), c(1L, 1L))) .dmar_refuse_argument(
          "With p = 1, 'true_cor_YX' must be a single correlation, that of Y with the one predictor.",
          arg = "true_cor_YX")
      }
      if (is.null(which_predictor)) which_predictor <- 1
    }
  }
  # Thanks to Jan Herman for the modification just above that lets p = 1 work.

  # The expected sample squared multiple correlation of a regression on p
  # predictors, as expected_R2() computes it: with p predictors for that of
  # Y, and with p - 1 for that of X_j on the other predictors, the
  # regression whose residual enters the standard error of b_j (with p = 1
  # there is no other predictor, and the value is 0). DMAR 1.0.x and MBESS
  # took p for both.
  Expected_R2 <- function(true_R2, N, p) {
    # gsl::hyperg_2F1(1, 1, 0.5 * (N + 1), true_R2) is another way to
    # obtain this 2F1 value; computed in base R via .hyperg_2F1() (no GSL
    # system dependency; see R/R2_internals.R).
    Value <- 1 - ((N - p - 1) / (N - 1)) * (1 - true_R2) *
      .hyperg_2F1(1, 1, 0.5 * (N + 1), true_R2)
    Value <- max(0, Value)
    return(Value)
  }
  if (noncentral == TRUE && is.null(true_sd_X)) {
    true_sd_X <- 1
  }
  if (noncentral == TRUE && is.null(true_sd_Y)) {
    true_sd_Y <- 1
  }
  if (noncentral == TRUE && true_sd_Y != 1) {
    .dmar_refuse_argument(
      "Since you've specified 'noncentral=TRUE', all variances should be one (your 'true_sd_Y' (i.e., the Y standard deviation) is not one), for a standardized solution.",
      arg = "true_sd_Y")
  }
  if (noncentral == TRUE && true_sd_X != 1) {
    .dmar_refuse_argument(
      "Since you've specified 'noncentral=TRUE', all variances should be one (your 'true_sd_X' (i.e., the X standard deviation) is not one), for a standardized solution.",
      arg = "true_sd_X")
  }
  if (is.null(p) && is.null(true_cor_XX)) {
    .dmar_refuse_argument(
      "Since true_cor_XX is not specified, you must specify 'p'.", arg = "p")
  }
  if (!is.null(true_cor_XX)) {
    if (!(sum(round(true_cor_XX, 5) == round(t(true_cor_XX), 5)) ==
      dim(as.matrix(true_cor_XX))[1] * dim(as.matrix(true_cor_XX))[2])) {
      .dmar_refuse_argument(
        "The correlation matrix, 'true_cor_XX' should be symmetric.",
        arg = "true_cor_XX")
    }
  }
  if (is.null(p) && !is.null(true_cor_XX)) {
    p <- dim(true_cor_XX)[1]
  }
  # An abbreviation is planned as the value it stands for (DMAR 1.0.x
  # planned "L" and "Up" for the full width; R/which_width_internals.R).
  given_width <- which_width
  which_width <- .match_which_width(which_width)
  if (which_width != "Full") {
    .dmar_refuse_argument(c(
      "'which_width' must be \"Full\" (or an abbreviation of it, such as ",
      "\"F\"): only the plan for the full width of the confidence interval ",
      "is implemented. \"", given_width, "\" was given",
      if (given_width != which_width)
        paste0(", which stands for \"", which_width, "\""), "."),
      arg = "which_width", value = given_width)
  }
  if (is.null(conf_level)) {
    if (!is.numeric(alpha_lower) || !is.numeric(alpha_upper)) {
      .dmar_refuse_argument(
        "Since 'conf_level' is not specified, you need to correctly specify 'alpha_lower' and 'alpha_upper'.",
        arg = c("conf_level", "alpha_lower", "alpha_upper"))
    }
    # A tail of 0 puts a limit at infinity (DMAR 1.0.x planned the average
    # tail); ci_reg_coef() takes tails up to .5.
    if (alpha_lower <= 0 || alpha_upper <= 0 || alpha_lower > .5 ||
        alpha_upper > .5 || alpha_lower + alpha_upper >= 1) .dmar_refuse_argument(
      "'alpha_lower' and 'alpha_upper' must each be above 0 (a tail of 0 puts a limit at infinity) and at most .5, as in ci_reg_coef(), with a sum below 1.",
      arg = c("alpha_lower", "alpha_upper"))
  }
  if (!is.null(conf_level)) {
    if (!is.null(alpha_lower) || !is.null(alpha_upper)) {
      .dmar_refuse_argument(c(
        "'alpha_lower' and 'alpha_upper' are used only when 'conf_level' is ",
        "NULL, and 'conf_level' is ", format(conf_level), " (its default is ",
        ".95). To plan with the two tail areas, give them with ",
        "conf_level = NULL; otherwise leave both out."),
        arg = c("conf_level", "alpha_lower", "alpha_upper"))
    }
    alpha_lower <- alpha_upper <- (1 - conf_level) / 2
  }
  # The floor on a tail area (R/dmar_conditions.R), checked here for the
  # central plans as well as the noncentral ones, whose limits come from
  # ci_nc_t().
  .dmar_check_tail_floor(alpha_lower, alpha_upper, conf_level = conf_level,
                         zero_ok = FALSE)
  if (!is.null(assurance)) {
    if ((assurance <= 0) || (assurance >=
      1)) {
      .dmar_refuse_argument(
        "The 'assurance' must either be NULL or some value greater than .50 and less than 1.",
        arg = "assurance")
    }
    if (assurance <= 0.5) {
      .dmar_refuse_argument(
        "The 'assurance' should be > .5 (but less than 1).", arg = "assurance")
    }
  }
  if ((!is.null(true_R2_j_X_without_j) && !is.null(true_R2_Y_X)) &
    (!is.null(true_cor_XX) & !is.null(true_cor_YX))) {
    .dmar_refuse_argument(
      "Since 'true_R2_j_X_without_j' and 'true_R2_Y_X' are specified, do not specify 'true_cor_XX' or 'true_cor_YX' (or vice versa).",
      arg = c("true_R2_j_X_without_j", "true_R2_Y_X", "true_cor_XX"))
  }
  if (!is.null(true_cor_XX) && !is.null(true_cor_YX)) {
    if (!is.null(true_R2_j_X_without_j) || !is.null(true_R2_Y_X)) {
      .dmar_refuse_argument(
        "Since 'true_cor_XX' and 'true_cor_YX' are specified, do not specify 'true_R2_Y_X' or 'true_R2_j_X_without_j'.",
        arg = c("true_cor_XX", "true_cor_YX", "true_R2_Y_X"))
    }
    if (is.null(which_predictor) || !is.numeric(which_predictor) ||
        length(which_predictor) != 1L || which_predictor < 1 ||
        which_predictor > p || which_predictor != round(which_predictor)) {
      .dmar_refuse_argument(
        sprintf("'which_predictor' must be a whole number between 1 and p = %d: on the covariance structure path it says which coefficient is planned for.", p),
        arg = "which_predictor")
    }
    true_R2_Y_X <- (true_cor_YX %*% solve(true_cor_XX) %*% true_cor_YX)
    true_R2_j_X_without_j <- 1 - 1 / solve(true_cor_XX)[
      which_predictor,
      which_predictor
    ]
    # if-statement added by Jan Herman to cope with the p=1 case
    if (p > 1) {
      R2_Y_X_without_j <- (true_cor_YX[-which_predictor] %*% solve(true_cor_XX[
        -which_predictor,
        -which_predictor
      ]) %*% true_cor_YX[-which_predictor])
    } else {
      R2_Y_X_without_j <- as.matrix(0) # If p=1 then R2_Y_X_without_j will be zero
    }
    # The correlations imply the standardized coefficient; on the scale that
    # true_sd_Y and true_sd_X describe, the coefficient is that times
    # true_sd_Y / true_sd_X (Kelley & Maxwell, 2003, p. 316), and it is on
    # that scale that a supplied true_reg_coef_j is checked.
    b_j_tmp <- (solve(true_cor_XX) %*% true_cor_YX)[which_predictor] *
      (true_sd_Y / true_sd_X)
    if (!is.null(true_reg_coef_j)) {
      if (!isTRUE(all.equal(as.numeric(true_reg_coef_j), as.numeric(b_j_tmp)))) {
        .dmar_refuse_argument(
          sprintf("'true_reg_coef_j' (%s) is not the coefficient the correlations 'true_cor_XX' and 'true_cor_YX' imply on the scale 'true_sd_Y' and 'true_sd_X' describe (%s); give the implied value or leave 'true_reg_coef_j' out.",
                  format(true_reg_coef_j, digits = 7), format(b_j_tmp, digits = 7)),
          arg = "true_reg_coef_j", value = true_reg_coef_j, implied = b_j_tmp)
      }
    }
  }
  if (!is.null(true_R2_j_X_without_j) && !is.null(true_R2_Y_X) && (is.null(true_cor_XX) |
    is.null(true_cor_YX))) {
    if (is.null(true_reg_coef_j)) {
      .dmar_refuse_argument(
        "Since 'true_cor_XX' and 'true_cor_YX' are not specified, implying 'true_R2_j_X_without_j' and 'true_R2_Y_X' are specified, 'true_reg_coef_j' must also be specified.",
        arg = "true_reg_coef_j")
    }
  }
  # The width of the interval ci_reg_coef() reports is the sum of its two
  # quantiles times the standard error; n0, the first sample size the
  # search tries, has the normal quantiles.
  q_sum <- function(quantile, ...)
    quantile(alpha_lower, ..., lower.tail = FALSE) +
      quantile(alpha_upper, ..., lower.tail = FALSE)
  # The squared standard error of b_j at sample size n, at the expected
  # squared multiple correlations (Kelley & Maxwell, 2003, Equation 2).
  central_var <- function(n) {
    E_rho2_Y_X <- Expected_R2(true_R2 = true_R2_Y_X, N = n, p = p)
    E_rho2_j_X_without_j <- Expected_R2(
      true_R2 = true_R2_j_X_without_j, N = n, p = p - 1)
    ((1 - E_rho2_Y_X) / ((1 - E_rho2_j_X_without_j) * (n - p - 1))) *
      (true_sd_Y^2 / true_sd_X^2)
  }
  # The plan without assurance: the smallest N of at least p + 2 (one
  # residual degree of freedom) at which that width is no wider than
  # `width`. The width does not increase with N (in every one of 1,866
  # planning sets checked from p + 2), so the doubling and bisection of
  # .reg_coef_first_n() (R/reg_coef_assurance_internals.R) find it.
  central_plan <- function(q) {
    n0 <- (q(qnorm) / width)^2 * (1 - true_R2_Y_X) /
      (1 - true_R2_j_X_without_j) * (true_sd_Y^2 / true_sd_X^2) + p + 1
    .reg_coef_first_n(function(n) q(qt, n - p - 1) * sqrt(central_var(n)) <=
                        width, lower = p + 2, guess = n0)
  }
  # Kelley and Maxwell (2003, Equation 3) at the plan n without assurance
  # (the width at the expected squared multiple correlations), with the
  # quantile sum q: the plan for the variance of b_j there multiplied by the
  # `assurance` quantile of a chi square on n - 1 degrees of freedom over
  # n - p - 1.
  chisq_plan <- function(q, n) {
    ceiling((q(qt, n - p - 1) / width)^2 * central_var(n) *
      qchisq(assurance, n - 1) + p + 1)
  }
  N <- central_plan(q_sum)
  # The returned table: the sample size, then, with assurance, the request
  # and the probability attained (NA for the published approximation,
  # which does not compute it; and, when it was simulated, G and the Monte
  # Carlo standard error as an attribute), with the record of how that
  # probability was found (R/ss_aipe_assurance_internals.R).
  plan_table <- function(n, attained = NA_real_, G_used = NULL, mc_se = NULL) {
    terms <- "necessary_N"
    values <- n
    if (!is.null(assurance)) {
      terms <- c(terms, "assurance", "assurance_attained")
      values <- c(values, assurance, attained)
    }
    if (!is.null(G_used)) {
      terms <- c(terms, "G")
      values <- c(values, G_used)
    }
    out <- .as_dmar_tbl(data.frame(term = terms, value = values),
                        conf_level = conf_level, subclass = "dmar_ss_aipe")
    if (!is.null(mc_se)) attr(out, "assurance_mc_se") <- mc_se
    if (!is.null(assurance)) {
      out <- .ss_aipe_assurance_record(
        out, if (step == "simulation") "simulation" else "analytic",
        switch(step, exact = "exact", simulation = "simulation",
               approximation = "not computed"))
    }
    out
  }
  if (noncentral == FALSE && is.null(assurance)) {
    return(plan_table(N))
  }
  if (noncentral == FALSE && !is.null(assurance)) {
    # The chi square step, for ss_aipe_src(noncentral = FALSE) with
    # assurance_method = "analytic".
    if (step == "approximation") return(plan_table(chisq_plan(q_sum, N)))
    # The analytic assurance step for the interval in the units of the data:
    # the smallest sample size at or above N, the plan without assurance, at
    # which the probability of an interval no wider than `width` (an F
    # probability; .reg_coef_raw_width_probability() in
    # R/reg_coef_assurance_internals.R) is at least `assurance`.
    R2 <- as.numeric(true_R2_Y_X)
    R2_j <- as.numeric(true_R2_j_X_without_j)
    assurance_step <- .reg_coef_assurance_search(N, function(n)
      .reg_coef_raw_width_probability(n, p, R2, R2_j, true_sd_Y / true_sd_X,
                                      width, alpha_lower, alpha_upper),
      assurance)
    return(plan_table(assurance_step$n, attained = assurance_step$attained))
  }
  if (noncentral == TRUE) {
    if (is.null(true_reg_coef_j)) {
      if (is.null(true_cor_XX) || is.null(true_cor_YX)) .dmar_refuse_argument(
        "'true_reg_coef_j' must be specified directly or implied by 'true_cor_XX' and 'true_cor_YX'.",
        arg = "true_reg_coef_j")
      true_reg_coef_j <- (solve(true_cor_XX) %*% true_cor_YX)[which_predictor]
    }
    # The width of the noncentral interval at sample size n, at the
    # expected squared multiple correlations, with R2_Y the population value
    # of that of Y (true_R2_Y_X, or a lower confidence limit of it below).
    nc_width <- function(n, R2_Y) {
      ci <- ci_reg_coef(
        reg_coef_j = true_reg_coef_j, se_reg_coef_j = NULL,
        sd_Y = true_sd_Y, sd_X = true_sd_X, N = n, p = p,
        R2_Y_X = Expected_R2(true_R2 = R2_Y, N = n, p = p),
        R2_j_X_without_j = Expected_R2(
          true_R2 = true_R2_j_X_without_j, N = n, p = p - 1),
        conf_level = NULL, R2_Y_X_without_j = NULL, t_value = NULL,
        alpha_lower = alpha_lower, alpha_upper = alpha_upper,
        noncentral = TRUE
      )
      ci$value[ci$term == "ci_upper"] - ci$value[ci$term == "ci_lower"]
    }
    # The plan without assurance: the smallest N of at least p + 2 at which
    # that width is no wider than `width`, found from the central plan N,
    # which is near it. The width does not increase with N (in every one of
    # 492 planning sets checked), and at N = p + 2, one degree of freedom,
    # ci_nc_t() gives its limits without complaint.
    N_NC <- .reg_coef_first_n(function(n) nc_width(n, true_R2_Y_X) <= width,
                              lower = p + 2, guess = N)
    if (noncentral == TRUE && is.null(assurance)) {
      return(plan_table(N_NC))
    }
    # The approximation of assurance_method = "analytic" for the noncentral
    # plan (Kelley, 2007a, 2007b, and MBESS): the plan without assurance again,
    # not below N_NC, with a lower confidence limit (at level `assurance`)
    # of the squared multiple correlation in place of its population value.
    nc_assurance_approx <- function() {
      R2_lower <- ci_R2(
        R2 = Expected_R2(true_R2 = true_R2_Y_X, N = N_NC, p = p), df_1 = p,
        df_2 = N_NC - p - 1, conf_level = assurance, F_value = NULL,
        N = NULL, p = NULL, alpha_lower = NULL, alpha_upper = NULL
      )[1, 2]
      .reg_coef_first_n(function(n) nc_width(n, R2_lower) <= width,
                        lower = N_NC)
    }
    if (step == "approximation") return(plan_table(nc_assurance_approx()))
    # The simulated assurance step for the noncentral interval of the
    # standardized coefficient, from the exact sampling
    # distribution (.reg_coef_std_assurance_n() in
    # R/reg_coef_assurance_internals.R). The distribution of the width also
    # depends on the correlation of Y with the targeted predictor: on the
    # covariance structure path it is given, and otherwise the step takes
    # the smallest probability over 21 equally spaced values of it across
    # the range the inputs allow, an approximate minimum over that range
    # (.reg_coef_std_ad_values()). The width's distribution for -beta with
    # the two tail probabilities exchanged is that for beta (reverse the
    # sign of Y), so a negative coefficient is planned as its magnitude with
    # the tails exchanged, and reversing the sign gives the same plan.
    R2 <- as.numeric(true_R2_Y_X)
    R2_j <- as.numeric(true_R2_j_X_without_j)
    beta <- as.numeric(true_reg_coef_j)
    if (beta^2 * (1 - R2_j) > R2 + 1e-8)
      .dmar_refuse_argument(c(
        "'true_reg_coef_j' squared times 1 - 'true_R2_j_X_without_j' is the ",
        "part of 'true_R2_Y_X' that the targeted predictor adds, and here it ",
        "is larger than 'true_R2_Y_X' (", format(beta^2 * (1 - R2_j)),
        " against ", format(R2), "), which no population allows. Check the ",
        "three values, or use assurance_method = \"analytic\" for the ",
        "approximation with a lower confidence limit of the squared ",
        "multiple correlation."),
        arg = c("true_reg_coef_j", "true_R2_Y_X", "true_R2_j_X_without_j"))
    cor_YXj <- if (!is.null(true_cor_XX) && !is.null(true_cor_YX))
      as.numeric(true_cor_YX)[which_predictor] else NULL
    ad_values <- .reg_coef_std_ad_values(R2, R2_j, beta, p, cor_YXj)
    a_lower <- alpha_lower
    a_upper <- alpha_upper
    if (beta < 0) {
      beta <- -beta
      ad_values <- -ad_values
      a_lower <- alpha_upper
      a_upper <- alpha_lower
    }
    assurance_step <- .reg_coef_std_assurance_n(N_NC, p, beta, R2, R2_j,
                                                ad_values, width, a_lower,
                                                a_upper, assurance, G)
    return(plan_table(assurance_step$n, attained = assurance_step$attained,
                      G_used = G, mc_se = assurance_step$mc_se))
  }
}
