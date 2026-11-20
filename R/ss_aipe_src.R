#' Sample Size Necessary for the Accuracy in Parameter Estimation Approach for a Standardized Regression Coefficient of Interest
#'
#' @description
#' Plans the sample size from the accuracy in parameter estimation (AIPE)
#' approach for a standardized regression coefficient of interest: the
#' smallest \emph{N} at which the confidence interval for the coefficient,
#' computed at the expected values of the sample squared multiple
#' correlations, is no wider than \code{width}, or, with \code{assurance},
#' the published approximation for an interval no wider than \code{width}
#' with at least that probability or, by simulation, the smallest \emph{N}
#' at which it is. Two intervals
#' can be planned for. The default (\code{noncentral = TRUE}) is the
#' noncentral \emph{t} interval for the standardized coefficient, the
#' interval \code{\link{ci_src}} reports; \code{noncentral = FALSE}
#' is the central \emph{t} interval of Kelley and Maxwell (2003, Equation 1),
#' the interval of \code{\link{ci_reg_coef}} with \code{noncentral = FALSE}
#' for standardized data. For the coefficient in the units of the data use
#' \code{\link{ss_aipe_rc}}.
#'
#' @param true_R2_Y_X Population value of the squared multiple correlation coefficient
#' @param true_R2_j_X_without_j Population value of the squared multiple correlation coefficient predicting the \emph{j}th predictor variable from the remaining \emph{p}-1 predictor variables
#' @param p The number of predictor variables
#' @param true_std_reg_coef_j Population value of the standardized regression coefficient for the \emph{j}th predictor variable (the predictor of interest)
#' @param width The desired width of the confidence interval
#' @param which_width Which width (\code{"Full"}, \code{"Lower"}, or \code{"Upper"}) the width refers to (at present, only \code{"Full"} can be specified)
#' @param true_sd_Y The population standard deviation of \emph{Y} (i.e., the dependent variable). It must be 1 (the default): the plan is for the standardized coefficient. Use \code{\link{ss_aipe_rc}} for the coefficient in the units of the data.
#' @param true_sd_X The population standard deviation of the \emph{j}th \emph{X} variable (i.e., the predictor variable of interest). It must be 1 (the default), as \code{true_sd_Y} must.
#' @param true_cor_XX Population correlation matrix for the \emph{p} predictor variables
#' @param true_cor_YX Population \emph{p} length vector of correlation between the dependent variable (\emph{Y}) and the \emph{p} independent variables
#' @param which_predictor Identifies which of the \emph{p} predictors is of interest
#' @param alpha_lower Type I error rate for the lower confidence interval limit, at least 1e-9 and at most .5; a tail area
#'   below 1e-9 is refused, as in \code{\link{ss_aipe_reg_coef}}
#' @param alpha_upper Type I error rate for the upper confidence interval limit, at least 1e-9 and at most .5
#' @param conf_level Desired level of confidence for the computed interval (i.e., 1 - the Type I error rate); a level above
#'   .999999998, whose tail areas are below 1e-9, is refused
#' @param assurance Optional probability with which the realized confidence interval is to be no wider than \code{width}; when
#'   \code{NULL} (the default) the plan is for the width of the interval at the expected values of the sample squared
#'   multiple correlations (see \code{\link{ss_aipe_reg_coef}}), not for the mean width over samples
#' @param noncentral \code{TRUE} (the default) plans for the noncentral
#'   \emph{t} interval for the standardized coefficient, the interval of
#'   \code{\link{ci_src}}; \code{FALSE} plans for the central \emph{t}
#'   interval of Kelley and Maxwell (2003, Equation 1), the estimate plus
#'   and minus the \emph{t} quantiles of its two tails times the standard error formed from the
#'   two sample squared multiple correlations.
#' @param assurance_method How the plan with \code{assurance} is made, used only with \code{assurance}. Neither interval
#'   has an exact form for the probability that it is no wider than \code{width}. \code{"analytic"} (the default),
#'   computed without simulation, is the published approximation: that of \pkg{MBESS} for the noncentral
#'   interval, and the chi square approximation of Kelley and Maxwell (2003, Equation 3) for the central one; neither
#'   computes that probability (\code{assurance_attained} is \code{NA}), and either can attain more or less than
#'   \code{assurance}. \code{"simulation"} is the smallest sample size, at or above the plan without \code{assurance},
#'   at which the probability estimated from \code{G} simulated samples is at least \code{assurance}, and reports that
#'   estimate as \code{assurance_attained} (see Details).
#' @param G Number of simulated samples at each sample size the assurance step tries, and in its final estimate of the
#'   attained probability; a whole number of at least 100 (default 20000). Used and checked only with \code{assurance} and
#'   \code{assurance_method = "simulation"}
#' @param seed Optional single number. When given, the assurance step is reproducible, and the random number stream of the
#'   calling session is left as it was. With \code{NULL} (the default) the draws come from the session's current stream.
#'   Used only where \code{G} is
#'
#' @details
#' Not all of the arguments need to be specified, only those that provide all of the necessary information
#' so that the sample size can be determined for the conditions specified.
#'
#' The coefficient planned for is the population standardized regression
#' coefficient of the \emph{j}th predictor, \eqn{\beta_j}, the quantity
#' \code{true_std_reg_coef_j} stands for; on the covariance structure path
#' it is the \emph{j}th element of \eqn{R_{XX}^{-1} r_{XY}}, with
#' \eqn{R_{XX}} given by \code{true_cor_XX} and \eqn{r_{XY}} by
#' \code{true_cor_YX}, and a supplied \code{true_std_reg_coef_j} must equal
#' it. A condition raised by \code{\link{ss_aipe_reg_coef}}, which makes
#' the plans without \code{assurance}, names the arguments of this function
#' (\code{true_std_reg_coef_j}, not \code{true_reg_coef_j}).
#'
#' \strong{The interval planned for.} Both intervals use the standard error
#' \eqn{\sqrt{(1 - R^2)/[(1 - R^2_{X_j \cdot X_{-j}})(N - p - 1)]}}, formed
#' from the sample squared multiple correlations of \emph{Y} with the
#' \emph{p} predictors and of \eqn{X_j} with the other \eqn{p - 1}, which
#' treats the sample standard deviations as constants (see
#' \code{\link{ci_src}}). With \code{noncentral = TRUE} the limits are the
#' limits for the noncentrality parameter at the observed \emph{t}
#' statistic times that standard error, and the plan without
#' \code{assurance} is the noncentral plan of \code{\link{ss_aipe_reg_coef}}
#' (\code{noncentral = TRUE}). With \code{noncentral = FALSE} the interval
#' is the estimate plus and minus \emph{t} quantiles on \eqn{N - p - 1}
#' degrees of freedom times that standard error (Kelley & Maxwell, 2003,
#' Equation 1), and the plan without \code{assurance} is theirs (Equation 2
#' with the \emph{t} quantile, at the expected squared multiple
#' correlations), the central plan of \code{\link{ss_aipe_reg_coef}} at unit
#' standard deviations, with the \emph{t} quantiles of its two tails.
#'
#' \strong{The assurance step.} With \code{assurance} the plan is by default (\code{assurance_method = "analytic"}) the
#' published approximation, computed without simulation: for the noncentral interval that of
#' \pkg{MBESS}, for a lower confidence limit of the squared multiple correlation, and for the central interval the chi
#' square approximation of Kelley and Maxwell (2003, Equation 3), the plan of \code{\link{ss_aipe_reg_coef}} without
#' \code{assurance} with the variance of the coefficient multiplied by the \code{assurance} quantile of a chi square on
#' \eqn{N - 1} degrees of freedom over \eqn{N - p - 1}. Neither computes the probability it attains, so
#' \code{assurance_attained} is \code{NA}, and the printed footer says so; both start from the plan without
#' \code{assurance}, which can be a case or more smaller than that of \pkg{MBESS} (see
#' \code{\link{ss_aipe_reg_coef}}).
#'
#' With \code{assurance_method = "simulation"} the plan is the smallest sample size, not below
#' the plan without \code{assurance}, at which the interval for the standardized coefficient is no wider
#' than \code{width} with at least that probability, for random predictors from a multivariate normal population. The
#' width of the noncentral interval depends on the observed \emph{t} statistic and on both sample squared multiple
#' correlations, and that of the central interval on the two squared multiple correlations, so for either interval the
#' probability is estimated from \code{G} samples drawn from the exact sampling distribution of the sums of squares
#' and products, the same samples at every sample size the search tries, so that the estimate changes smoothly with the
#' sample size. The search assumes that the estimate does not decrease as the sample size grows, which held in every set
#' checked (see "The assurance step" in \code{\link{ss_aipe_reg_coef}}); \code{assurance_attained} is the estimate at
#' the returned sample size, at least \code{assurance},
#' with its Monte Carlo standard error in the attribute \code{assurance_mc_se}. The distribution of the
#' width also depends on the correlation of \emph{Y} with the targeted predictor: on the covariance structure path it is
#' given, and on the squared multiple correlations path the step uses the smallest probability over 21 equally spaced
#' values of it, from the smallest to the largest the three planning values allow (the two ends with \eqn{p = 2}, the
#' only values possible there), an approximate minimum over that range. Give \code{seed} for a plan that can be
#' reproduced. This is the simulation
#' Kelley and Maxwell (2008, pp. 186--187) described for the standardized coefficient, for which they wrote that no
#' satisfactory method then existed.
#'
#' The central interval for the standardized coefficient is not the central
#' interval for the coefficient in the units of the data that
#' \code{\link{ss_aipe_reg_coef}} plans for with \code{noncentral = FALSE},
#' and the exact probability that planner computes is that of the interval
#' in the units of the data. The standard error of that interval is formed
#' from the residual sums of squares of \emph{Y} and of \eqn{X_j}; that of
#' the standardized interval divides each by the total sum of squares of its
#' variable, which gives \eqn{1 - R^2} and
#' \eqn{1 - R^2_{X_j \cdot X_{-j}}}, and these vary less from sample to
#' sample, so the standardized interval needs fewer cases to be as narrow
#' with the same probability. For two predictors correlated .3 with
#' correlations of .4 and .3 with \emph{Y}, a width of .3, and an assurance
#' of .90, this function returns 171 with \code{assurance_method = "simulation"}, \code{seed = 1}, and the default
#' \code{G}, at which the central interval for the standardized coefficient
#' is no wider than .3 with probability .914 (.895 at 170). The plan is
#' simulated, so another seed can give one a case away: over seeds 1 to 30,
#' 29 gave 171 and one gave 170. The central plan of
#' \code{\link{ss_aipe_reg_coef}} there is 187 (.997 for this interval), and
#' the chi square approximation of Kelley and Maxwell (2003, Equation 3,
#' with the \emph{t} quantile in place of the normal one), the default plan
#' here, is 181 (.986). For
#' their five exchangeable predictors (pp. 312--313), correlated .4 with one
#' another and .3 with \emph{Y}, at a width of .2 and an assurance of .80,
#' the simulated plan is 478 (with each of seeds 1 to 30; .807, and .797 at
#' 477), where the approximation gives 486 (.877) and
#' \code{\link{ss_aipe_reg_coef}} 491 (.911). Each probability is the
#' proportion of intervals no wider than the target among samples of raw
#' multivariate normal data (50,000 or 100,000 samples at each \emph{N},
#' Monte Carlo standard error at most .0014), standardized and fit by least
#' squares without DMAR. Kelley and Maxwell's own simulation of that example
#' (p. 317) reported 81.64 percent of the intervals no wider than .2 at
#' their plan of 485 (made with the normal quantile) and a mean half-width
#' of .098; raw multivariate normal data give .871 there with the \emph{t}
#' quantile of their Equation 1 (.887 with the normal one) and a mean
#' half-width of .097 (50,000 samples, Monte Carlo standard error at most
#' .0015).
#'
#' @return
#' A \code{data.frame} with columns \code{term} and \code{value}: \code{necessary_N}, the necessary total sample
#' size; with \code{assurance}, the rows \code{assurance} and \code{assurance_attained}, which is \code{NA} for the
#' published approximation of \code{assurance_method = "analytic"} and, with \code{"simulation"}, the simulated
#' probability at \code{necessary_N} that the realized interval is no wider than \code{width}, with its Monte Carlo
#' standard error in the attribute \code{assurance_mc_se} and the row \code{G}. With \code{assurance} the printed
#' footer names the assurance method and says how \code{assurance_attained} was found, or that the approximation does
#' not compute it; the attributes \code{assurance_method} and \code{assurance_probability} (\code{"simulation"} or
#' \code{"not computed"}) record the same.
#'
#' @references
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
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @note
#' This function calls upon \code{\link{ss_aipe_reg_coef}} in DMAR, which
#' plans the coefficient in the units of the data by default, but has a
#' different naming scheme. See \code{\link{ss_aipe_reg_coef}} for more
#' details.
#'
#' @section Warning: As discussed in Kelley and Maxwell (2008), the sample size planning approach from the AIPE perspective used in this function is only an approximation.
#'
#' @seealso
#' \code{\link{ss_aipe_src_sensitivity}}, \code{\link{ci_src}}, \code{\link{ci_reg_coef}}, \code{\link{ci_nc_t}}, \code{\link{ss_aipe_reg_coef}}, \code{\link{ss_aipe_rc}}
#'
#' @examples
#' # Exchangeable correlation structure
#' true_cor_YX <- c(.3, .3, .3, .3, .3)
#' true_cor_XX <- rbind(c(1, .5, .5, .5, .5), c(.5, 1, .5, .5, .5), c(.5, .5, 1, .5, .5),
#'                 c(.5, .5, .5, 1, .5), c(.5, .5, .5, .5, 1))
#'
#' ss_aipe_src(width = .1, which_width = "Full", true_sd_Y = 1, true_sd_X = 1,
#'             true_cor_XX = true_cor_XX,
#'             true_cor_YX = true_cor_YX, which_predictor = 1, conf_level = 1 - .05)
#'
#' # With an assurance of .85: the approximation with a lower confidence
#' # limit of the squared multiple correlation, and the plan found by
#' # simulation. G = 2000 keeps the example quick; a reported plan deserves
#' # the default G = 20000.
#' ss_aipe_src(width = .1, true_cor_XX = true_cor_XX, true_cor_YX = true_cor_YX,
#'             which_predictor = 1, assurance = .85)
#' plan <- ss_aipe_src(width = .1, which_width = "Full", true_sd_Y = 1, true_sd_X = 1,
#'                     true_cor_XX = true_cor_XX,
#'                     true_cor_YX = true_cor_YX, which_predictor = 1, conf_level = 1 - .05,
#'                     assurance = .85, assurance_method = "simulation",
#'                     G = 2000, seed = 113)
#' plan
#' attr(plan, "assurance_mc_se")
#'
#' # The central t interval of Kelley and Maxwell (2003, Equation 1): the
#' # plan without assurance, and with an assurance of .85 by their chi
#' # square approximation and by simulation.
#' ss_aipe_src(width = .1, true_cor_XX = true_cor_XX, true_cor_YX = true_cor_YX,
#'             which_predictor = 1, noncentral = FALSE)
#' ss_aipe_src(width = .1, true_cor_XX = true_cor_XX, true_cor_YX = true_cor_YX,
#'             which_predictor = 1, noncentral = FALSE, assurance = .85)
#' ss_aipe_src(width = .1, true_cor_XX = true_cor_XX, true_cor_YX = true_cor_YX,
#'             which_predictor = 1, noncentral = FALSE, assurance = .85,
#'             assurance_method = "simulation", G = 2000, seed = 113)
#'
#' @keywords design
#'
#' @seealso \code{\link{design_consequences}} for what a chosen design delivers:
#'   power, the Type S (sign) and Type M (exaggeration) errors of the
#'   significance filter, and the expected confidence interval width.
#'
#' @export


ss_aipe_src <- function(true_R2_Y_X = NULL, true_R2_j_X_without_j = NULL, p = NULL,
                        true_std_reg_coef_j = NULL, width, which_width = "Full", true_sd_Y = 1, true_sd_X = 1, true_cor_XX = NULL,
                        true_cor_YX = NULL, which_predictor = NULL, alpha_lower = NULL, alpha_upper = NULL,
                        conf_level = .95, assurance = NULL, noncentral = TRUE,
                        assurance_method = c("analytic", "simulation"),
                        G = 20000, seed = NULL) {
  .dmar_check_conf_level(conf_level, null_ok = TRUE)
  if (!missing(width)) .dmar_check_width(width)
  .dmar_check_assurance_single(assurance)
  .dmar_check_numbers("ss_aipe_src", c("true_sd_Y", "true_sd_X", "which_predictor"))
  if (!is.logical(noncentral) || length(noncentral) != 1L || is.na(noncentral))
    .dmar_refuse_argument(c(
      "'noncentral' must be TRUE (the noncentral t interval for the ",
      "standardized coefficient) or FALSE (the central t interval, Kelley ",
      "and Maxwell, 2003, Equation 1)."), arg = "noncentral", value = noncentral)
  # Neither interval for the standardized coefficient has an exact form for
  # the probability that it is no wider than 'width', so both offer
  # "simulation" (R/ss_aipe_assurance_internals.R).
  assurance_method <- .ss_aipe_assurance_method(
    assurance_method, "ss_aipe_src", simulation = TRUE)
  # G sizes the simulated assurance step, and is checked, and the seed set,
  # only when that step runs.
  simulate <- !is.null(assurance) && assurance_method == "simulation"
  if (simulate) {
    G <- .dmar_check_count(G, "G", min = 100L, reason = paste0(
      "The assurance step estimates a probability from 'G' simulated ",
      "samples at each sample size it tries; with 100 the standard error ",
      "of that estimate is at most .05."))
    .dmar_local_seed(seed)
  }
  # The plan is for the standardized coefficient, so both standard
  # deviations are 1.
  for (nm in c("true_sd_Y", "true_sd_X")) {
    v <- get(nm)
    if (!is.null(v) && !isTRUE(all.equal(as.numeric(v), 1))) .dmar_refuse_argument(
      c("'", nm, "' must be 1: ss_aipe_src() plans for the standardized ",
        "coefficient, whose scale has unit standard deviations. For the ",
        "coefficient in the units of the data use ss_aipe_rc()."),
      arg = nm, value = v)
  }
  # A plan of the engine of ss_aipe_reg_coef() (R/ss_aipe_reg_coef.R), whose
  # true_reg_coef_j is this function's true_std_reg_coef_j, with the
  # assurance step `step`; its conditions are stated in this function's
  # names. The engine is called directly because the analytic method of
  # the central interval here is the chi square step of Kelley and Maxwell
  # (2003), which ss_aipe_reg_coef(), whose central interval has an exact
  # probability, does not offer.
  engine <- function(nc, a, step) .reg_coef_in_wrapper_names(
    .ss_aipe_reg_coef_plan(
      true_R2_Y_X = true_R2_Y_X, true_R2_j_X_without_j = true_R2_j_X_without_j, p = p, true_reg_coef_j = true_std_reg_coef_j,
      width = width, which_width = which_width, true_sd_Y = true_sd_Y, true_sd_X = true_sd_X,
      true_cor_XX = true_cor_XX, true_cor_YX = true_cor_YX, which_predictor = which_predictor, noncentral = nc,
      alpha_lower = alpha_lower, alpha_upper = alpha_upper, conf_level = conf_level,
      assurance = a, step = step, G = G
    ),
    c(true_reg_coef_j = "true_std_reg_coef_j", ss_aipe_reg_coef = "ss_aipe_src"))
  # The noncentral interval is the noncentral plan of the engine: with
  # assurance_method = "analytic" the approximation of DMAR 1.0.x and
  # MBESS, and with "simulation" the simulated step. The central interval
  # without assurance is the engine's central plan, the plan of Kelley and
  # Maxwell (2003, Equation 2, with the t quantile at the expected squared
  # multiple correlations), which at unit standard deviations is the same
  # for the interval in the units of the data and for the standardized one;
  # with assurance_method = "analytic" it is their Equation 3, the chi
  # square step.
  if (noncentral || !simulate)
    return(engine(noncentral, assurance,
                  if (simulate) "simulation" else "approximation"))
  .ss_aipe_src_central_assurance(
    engine = engine, true_R2_Y_X = true_R2_Y_X,
    true_R2_j_X_without_j = true_R2_j_X_without_j, p = p,
    true_std_reg_coef_j = true_std_reg_coef_j, width = width,
    true_cor_XX = true_cor_XX, true_cor_YX = true_cor_YX,
    which_predictor = which_predictor, alpha_lower = alpha_lower,
    alpha_upper = alpha_upper, conf_level = conf_level, assurance = assurance,
    G = G)
}

# .ss_aipe_src_central_assurance(engine, ...): the simulated assurance step
# of ss_aipe_src(noncentral = FALSE), assurance_method = "simulation", with
# the arguments of ss_aipe_src() (checked there, where the seed is set) and
# `engine`, its call of the engine of ss_aipe_reg_coef(). The probability
# that the central interval for the standardized coefficient
# (ci_reg_coef(noncentral = FALSE) on standardized data, Kelley and
# Maxwell, 2003, Equation 1) is no wider than `width` is not that of the
# interval in the units of the data, whose F law ss_aipe_reg_coef() uses:
# the standard error is formed from the two sample squared multiple
# correlations, whose sampling distribution is not the F law. It is
# simulated from the exact sampling distribution, from the draws the
# noncentral step uses (.reg_coef_std_assurance_n(interval = "central") in
# R/reg_coef_assurance_internals.R), and reported with G and its Monte
# Carlo standard error as that step reports them. Not exported.
.ss_aipe_src_central_assurance <- function(engine, true_R2_Y_X,
                                           true_R2_j_X_without_j, p,
                                           true_std_reg_coef_j, width,
                                           true_cor_XX, true_cor_YX,
                                           which_predictor, alpha_lower,
                                           alpha_upper, conf_level,
                                           assurance, G) {
  # The engine checks every argument (the tails among them) in the plan
  # without assurance, where the search starts.
  n_0 <- engine(FALSE, NULL, "approximation")$value[1]
  if (is.null(conf_level)) {
    a_lower <- alpha_lower
    a_upper <- alpha_upper
  } else {
    a_lower <- a_upper <- (1 - conf_level) / 2
  }
  pop <- .reg_coef_std_population(true_R2_Y_X, true_R2_j_X_without_j, p,
                                  true_std_reg_coef_j, true_cor_XX,
                                  true_cor_YX, which_predictor)
  beta <- pop$beta
  if (beta^2 * (1 - pop$R2_j) > pop$R2 + 1e-8)
    .dmar_refuse_argument(c(
      "'true_std_reg_coef_j' squared times 1 - 'true_R2_j_X_without_j' is ",
      "the part of 'true_R2_Y_X' that the targeted predictor adds, and here ",
      "it is larger than 'true_R2_Y_X' (", format(beta^2 * (1 - pop$R2_j)),
      " against ", format(pop$R2), "), which no population allows. Check ",
      "the three values, or use assurance_method = \"analytic\" for the ",
      "approximation of Kelley and Maxwell (2003, Equation 3)."),
      arg = c("true_std_reg_coef_j", "true_R2_Y_X", "true_R2_j_X_without_j"))
  ad_values <- .reg_coef_std_ad_values(pop$R2, pop$R2_j, beta, pop$p,
                                       pop$cor_YXj)
  # The distribution of the width for -beta is that for beta (reverse the
  # sign of Y), and the central width does not depend on which tail is
  # which, so a negative coefficient is planned as its magnitude.
  if (beta < 0) {
    beta <- -beta
    ad_values <- -ad_values
  }
  step <- .reg_coef_std_assurance_n(n_0, pop$p, beta, pop$R2, pop$R2_j,
                                    ad_values, width, a_lower, a_upper,
                                    assurance, G, interval = "central")
  out <- .as_dmar_tbl(
    data.frame(term = c("necessary_N", "assurance", "assurance_attained", "G"),
               value = c(step$n, assurance, step$attained, G)),
    conf_level = conf_level, subclass = "dmar_ss_aipe")
  attr(out, "assurance_mc_se") <- step$mc_se
  .ss_aipe_assurance_record(out, "simulation", "simulation")
}
