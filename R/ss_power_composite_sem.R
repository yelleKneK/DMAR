#' Sample Size or Composite Power for a Set of SEM Parameters
#'
#' Determine the necessary sample size for a structural equation model study
#' so that every parameter of interest is statistically significant in the
#' same study with a desired probability, or, given a sample size, return
#' that probability. Composite power is the probability that all of the named
#' parameters are significant jointly, the quantity a design must be planned
#' against when its conclusion requires more than one result to hold at once:
#' a study can have adequate power for each hypothesis on its own and still
#' be underpowered for the conclusion that rests on all of them together
#' (Maxwell, 2004). The parameters of interest are any labeled parameters of
#' a lavaan analysis model, structural paths, loadings, covariances, or
#' quantities defined with \code{:=} such as an indirect effect, and any
#' subset of them can make up the composite.
#'
#' @param model A single character string giving the free analysis model in
#'   lavaan model syntax (see \code{\link[lavaan]{model.syntax}}), the model
#'   that would be fit to the data. Each parameter of interest must carry a
#'   parameter label so it can be referred to by name, for example
#'   \code{"f2 ~ b*f1"} labels the structural path \code{b}, and
#'   \code{"ab := a*b"} defines an indirect effect from the labeled paths
#'   \code{a} and \code{b}.
#' @param Sigma Population covariance matrix of the observed variables, with
#'   row and column names matching the observed variables in \code{model}. It
#'   is typically obtained from a fully fixed population model via
#'   \code{\link{cov_sem}}. Supply exactly one of \code{Sigma} or
#'   \code{pop_model}.
#' @param pop_model A single character string giving the population model in
#'   lavaan model syntax with every parameter fixed to its population value,
#'   from which \code{\link{cov_sem}} derives \code{Sigma} (and the
#'   population means, when the model has a mean structure). This is where
#'   the purported population values of the parameters of interest are
#'   chosen; they are values the researcher posits (from theory, prior
#'   studies, or pilot data), never sample estimates. Supply exactly one of
#'   \code{Sigma} or \code{pop_model}.
#' @param mu Optional population means of the observed variables, used with
#'   \code{Sigma}: a named numeric vector with one entry per observed
#'   variable, or an unnamed vector in the row order of \code{Sigma}. The
#'   default \code{NULL} is zero means. Means matter only when the analysis
#'   model has a mean structure (an intercept term such as \code{s ~ 1}, as
#'   in a latent growth curve model); when \code{pop_model} is supplied its
#'   mean structure provides the means and \code{mu} must not also be given.
#' @param parameters Character vector of the parameter labels that make up
#'   the composite. The default \code{NULL} uses every user-labeled parameter
#'   in \code{model}, in order of appearance, so labeling exactly the
#'   parameters of interest is the simplest way to state the set.
#' @param desired_power Desired composite statistical power (default 0.85).
#'   Used only when \code{N} is \code{NULL}.
#' @param alpha_level Type I error rate for each individual two-sided Wald
#'   \emph{z} test (default 0.05), the per-test rate, not a rate for the
#'   composite event.
#' @param N Sample size; if supplied, the realized composite power at that
#'   \emph{N} is returned rather than a sample size planned.
#' @param G Number of converged Monte Carlo replications per evaluated sample
#'   size (default 1000). The simulation error of each estimated power is
#'   about \eqn{\sqrt{p(1 - p)/G}}; raise \code{G} for a sharper answer.
#' @param seed Optional integer seed for reproducibility. The default
#'   \code{NULL} uses the current state of the random number generator; a
#'   supplied seed is set internally and the prior state restored on exit.
#' @param \dots Additional arguments passed to \code{\link[lavaan]{sem}},
#'   both when the population values are resolved and for every Monte Carlo
#'   fit (for example \code{std.lv = TRUE} or \code{missing = "listwise"}).
#'   A robust estimator such as \code{"MLM"} cannot be used here: the same
#'   arguments reach the setup fit, which is always from summary statistics,
#'   and lavaan refuses a robust estimator there. The Monte Carlo data are
#'   multivariate normal by construction, so a robust estimator would buy
#'   nothing.
#'
#' @details
#' Analytic sample size planning methods in SEM exist for a single targeted
#' parameter (Satorra & Saris, 1985; Lai & Kelley, 2011) or for overall model
#' fit (MacCallum, Browne, & Sugawara, 1996; \code{\link{ss_power_sem}}), but
#' most studies state several hypotheses and support their conclusion only
#' when all of them hold. This function plans for that case by a priori Monte
#' Carlo simulation (Muthén & Muthén, 2002; Maxwell,
#' Kelley, & Rausch, 2008): for a candidate \emph{N}, \code{G} data sets are
#' drawn from the multivariate normal population with covariance matrix
#' \code{Sigma}, the analysis model is fit to each, and each parameter of
#' interest is tested with its two-sided Wald \emph{z} test at
#' \code{alpha_level}. The proportion of replications in which every
#' parameter is significant estimates the composite power, and the
#' per-parameter proportions estimate the marginal powers. Because the
#' estimates share one fitted model, the tests are dependent; the simulation
#' reflects that dependence exactly, at the stated \emph{N}, with no
#' asymptotic shortcut.
#'
#' The composite event is contained in each marginal event, so composite
#' power is at most the smallest marginal power: the weakest parameter
#' governs the design, and the marginal \code{power_<label>} rows show which
#' parameter that is.
#'
#' When \code{N} is \code{NULL} the necessary sample size is searched for.
#' The search starts where the product of the marginal Wald powers (an
#' independence approximation computed from the asymptotic variances, spent
#' before any simulation) reaches \code{desired_power}, brackets the crossing
#' geometrically, and bisects to adjacent integers, each candidate evaluated
#' with its own \code{G} replications. A planning call therefore fits the
#' analysis model several thousand times, which is why the examples on this
#' page are shown but not run.
#'
#' @section Monte Carlo Precision:
#' Each reported power is a proportion of \code{G} replications, with
#' simulation standard error about \eqn{\sqrt{p(1 - p)/G}}; the
#' \code{composite_power_mc_se} row reports it for the composite. The
#' necessary sample size inherits that uncertainty: near the target the power
#' curve is flat enough that neighboring \emph{N} are separated by less than
#' the simulation error, so repeated calls with different seeds return
#' slightly different sizes. Raising \code{G} narrows the spread; reporting
#' the seed makes a plan reproducible. A proportion of \code{G}
#' replications takes only the values \eqn{0, 1/G, \ldots, 1}, so a
#' \code{desired_power} above \eqn{1 - 1/G} is refused with a message
#' saying how large \code{G} must be for that target; the same
#' resolution guard applies to \code{\link{ss_aipe_composite_sem}}'s
#' \code{assurance}.
#'
#' @note
#' A replication whose fit does not converge, or converges without a usable
#' standard error for some parameter of interest, is discarded and fresh data
#' are drawn, up to \code{20 * G} attempts per evaluated sample size; the
#' reported powers condition on convergence. When fewer than \code{G}
#' replications converge within the cap, a single warning is issued and the
#' summary is based on the converged replications (their count is the
#' \code{converged_replications} row). Frequent nonconvergence at small
#' \emph{N} is itself design information: a sample size at which the model
#' rarely converges is too small in a sense that precedes power.
#'
#' @return A \code{data.frame} with \code{term} and \code{value} columns: the
#'   \code{necessary_N} (or supplied \code{specified_N}), the
#'   \code{composite_power} and its simulation standard error
#'   \code{composite_power_mc_se}, then for each parameter its marginal
#'   \code{power_<label>} and purported \code{population_<label>} value under
#'   the analysis model, followed by \code{alpha_level}, the requested
#'   \code{replications}, the \code{converged_replications} the summary is
#'   based on, and, when a size was planned, the \code{desired_power}. The
#'   result carries the \code{dmar_ss_power} class, so
#'   \code{\link[generics]{tidy}} and \code{\link[generics]{glance}}
#'   summarize the sample size and the composite power in broom convention.
#'
#' @references
#' Lai, K., & Kelley, K. (2011). Accuracy in parameter estimation for
#'   targeted effects in structural equation modeling: Sample size
#'   planning for narrow confidence intervals.
#'   \emph{Psychological Methods, 16}(2), 127--148. \doi{10.1037/a0021764}
#'
#' MacCallum, R. C., Browne, M. W., & Sugawara, H. M. (1996). Power
#'   analysis and determination of sample size for covariance structure
#'   modeling. \emph{Psychological Methods, 1}(2), 130--149.
#'   \doi{10.1037/1082-989X.1.2.130}
#'
#' Maxwell, S. E. (2004). The persistence of underpowered studies in
#'   psychological research: Causes, consequences, and remedies.
#'   \emph{Psychological Methods, 9}(2), 147--163.
#'   \doi{10.1037/1082-989X.9.2.147}
#'
#' Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). \emph{Designing
#'   experiments and analyzing data: A model comparison perspective}
#'   (4th ed.). Routledge. (See Chapter 3 on statistical power.)
#'
#' Maxwell, S. E., Kelley, K., & Rausch, J. R. (2008). Sample size planning
#'   for statistical power and accuracy in parameter estimation.
#'   \emph{Annual Review of Psychology, 59}, 537--563.
#'   \doi{10.1146/annurev.psych.59.103006.093735}
#'
#' Muthén, L. K., & Muthén, B. O. (2002). How to use a
#'   Monte Carlo study to decide on sample size and determine power.
#'   \emph{Structural Equation Modeling, 9}(4), 599--620.
#'   \doi{10.1207/S15328007SEM0904_8}
#'
#' Rosseel, Y. (2012). lavaan: An R package for structural equation modeling.
#'   \emph{Journal of Statistical Software, 48}(2), 1--36.
#'   \doi{10.18637/jss.v048.i02}
#'
#' Satorra, A., & Saris, W. E. (1985). Power of the likelihood ratio test in
#'   covariance structure analysis. \emph{Psychometrika, 50}(1), 83--90.
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso \code{\link{ss_aipe_composite_sem}} for the same set of
#'   parameters planned for accuracy in parameter estimation (AIPE) instead
#'   of significance; \code{\link{cov_sem}} for deriving \code{Sigma} from a
#'   fully fixed population model; \code{\link{ss_power_sem}} for overall
#'   model fit; \code{\link{ss_aipe_sem_path}} for a single targeted path;
#'   \code{\link{ss_power_composite_anova}} and its siblings for composite
#'   power in ANOVA and ANCOVA designs, where the composite is evaluated by
#'   quadrature rather than simulation.
#'
#' @examples
#' # A three-factor model whose conclusion rests on three structural paths
#' # at once: f1 predicting f2, f2 predicting f3, and f1 predicting f3
#' # directly. Composite power here is a simulated quantity, so every call
#' # refits the analysis model G times and a planning search refits it several
#' # thousand times. The worked example that follows is therefore shown rather
#' # than run.
#' #
#' # The population model fixes every parameter to its purported population
#' # value.
#' #   pop_model <- "
#' #     f1 =~ 1*y1 + 0.8*y2 + 0.8*y3
#' #     f2 =~ 1*y4 + 0.8*y5 + 0.8*y6
#' #     f3 =~ 1*y7 + 0.8*y8 + 0.8*y9
#' #     f2 ~ 0.4*f1
#' #     f3 ~ 0.3*f2 + 0.25*f1
#' #     f1 ~~ 1*f1
#' #     f2 ~~ 0.84*f2
#' #     f3 ~~ 0.8*f3
#' #     y1 ~~ 0.5*y1; y2 ~~ 0.5*y2; y3 ~~ 0.5*y3
#' #     y4 ~~ 0.5*y4; y5 ~~ 0.5*y5; y6 ~~ 0.5*y6
#' #     y7 ~~ 0.5*y7; y8 ~~ 0.5*y8; y9 ~~ 0.5*y9
#' #   "
#' #
#' # The analysis model is free; the labels name the parameters of interest.
#' #   analysis_model <- "
#' #     f1 =~ y1 + y2 + y3
#' #     f2 =~ y4 + y5 + y6
#' #     f3 =~ y7 + y8 + y9
#' #     f2 ~ a*f1
#' #     f3 ~ b*f2 + c*f1
#' #   "
#' #
#' # Realized composite power at N = 200. The probability that all three
#' # paths come out significant in the same study is lower than the marginal
#' # power of any one of them: the composite event sits inside each marginal
#' # event, so the weakest parameter governs the design.
#' #   set.seed(113)
#' #   ss_power_composite_sem(model = analysis_model, pop_model = pop_model,
#' #                          N = 200, G = 1000)
#' #
#' # Leaving N out plans the necessary sample size for a desired composite
#' # power instead, here over the two structural paths a and b with c left
#' # out of the composite. That search evaluates a sequence of candidate
#' # sample sizes, each with its own G replications, so it costs several
#' # thousand model fits:
#' #   set.seed(113)
#' #   ss_power_composite_sem(model = analysis_model, pop_model = pop_model,
#' #                          parameters = c("a", "b"),
#' #                          desired_power = 0.80, G = 1000)
#'
#' @keywords design multivariate htest
#'
#' @family sample size for power
#'
#' @family composite power
#'
#' @export
ss_power_composite_sem <- function(model, Sigma = NULL, pop_model = NULL,
                                   mu = NULL, parameters = NULL,
                                   desired_power = 0.85, alpha_level = 0.05,
                                   N = NULL, G = 1000, seed = NULL, ...) {
  if (!requireNamespace("lavaan", quietly = TRUE)) {
    stop("The package 'lavaan' is needed; please install it and try again.",
         call. = FALSE)
  }
  if (!is.numeric(alpha_level) || length(alpha_level) != 1L ||
      is.na(alpha_level) || alpha_level <= 0 || alpha_level >= 1) {
    stop("'alpha_level' must be a single number in (0, 1).", call. = FALSE)
  }
  if (!is.numeric(G) || length(G) != 1L || is.na(G) || G < 10 ||
      G != round(G)) {
    stop("'G' must be a single whole number of at least 10.", call. = FALSE)
  }
  G <- as.integer(G)

  if (!is.null(seed)) {
    if (!is.numeric(seed) || length(seed) != 1L || is.na(seed)) {
      stop("'seed' must be NULL or a single number.", call. = FALSE)
    }
    if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
      old_seed <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
      on.exit(assign(".Random.seed", old_seed, envir = globalenv()),
              add = TRUE)
    } else {
      on.exit(if (exists(".Random.seed", envir = globalenv(),
                         inherits = FALSE))
        rm(".Random.seed", envir = globalenv()), add = TRUE)
    }
    set.seed(seed)
  }

  setup <- .composite_sem_setup(model, Sigma, pop_model, parameters,
                                mu = mu, ...)
  labels <- setup$labels
  k <- length(labels)
  N_min <- max(10L, setup$n_obs_vars + 2L)
  z_crit <- stats::qnorm(1 - alpha_level / 2)

  # Evaluations where fewer than G replications converged within the attempts
  # cap; a single summary warning is issued at the end rather than one per
  # evaluation.
  short_evals <- 0L
  evaluate <- function(n) {
    mc <- .composite_sem_mc(model, setup$Sigma, setup$mu, labels, n, G, ...)
    if (mc$converged == 0L) {
      stop("No replication converged at N = ", n, " within ", mc$attempts,
           " attempts; the model cannot be fit at so small a sample size. ",
           "Reconsider the model or the candidate sample size.",
           call. = FALSE)
    }
    if (mc$converged < G) short_evals <<- short_evals + 1L
    significant <- abs(mc$est / mc$se) > z_crit
    list(marginal = colMeans(significant),
         composite = mean(rowSums(significant) == k),
         converged = mc$converged)
  }

  if (!is.null(N)) {
    if (!is.numeric(N) || length(N) != 1L || is.na(N) || N != round(N) ||
        N < N_min) {
      stop("'N' must be a single whole number of at least ", N_min,
           " for this model.", call. = FALSE)
    }
    size_term <- "specified_N"
    size <- as.integer(N)
    res <- evaluate(size)
    desired_power <- NULL
  } else {
    if (!is.numeric(desired_power) || length(desired_power) != 1L ||
        is.na(desired_power) || desired_power <= 0 || desired_power >= 1) {
      stop("'desired_power' must be a single number in (0, 1).",
           call. = FALSE)
    }
    # The composite power is a proportion of G replications, so it can only
    # take the values 0, 1/G, ..., 1. A target within 1/G of 1 is met the
    # moment every replication happens to be significant, which the search
    # reaches at a sample size well below the one the target requires, and
    # the result then reports a power of 1 with a simulation standard error
    # of 0. Refuse the request rather than return a size the Monte Carlo
    # cannot justify.
    .composite_sem_check_resolution(desired_power, G, "desired_power")
    N_start <- .composite_sem_start_power(setup$theta, setup$h, alpha_level,
                                          desired_power, N_min)
    if (N_start >= 1e6) {
      stop("Even N = 1e6 does not reach 'desired_power' under the analytic ",
           "Wald approximation; at least one parameter of interest is zero ",
           "or nearly zero in the population implied by 'Sigma'. Drop it ",
           "from 'parameters' or revisit the population values.",
           call. = FALSE)
    }
    searched <- .composite_sem_search(
      evaluate, function(e) e$composite >= desired_power, N_start, N_min)
    size_term <- "necessary_N"
    size <- searched$N
    res <- searched$eval
  }

  if (short_evals > 0L) {
    warning("In ", short_evals, " Monte Carlo evaluation",
            if (short_evals > 1L) "s" else "", " fewer than G = ", G,
            " replications converged within the attempts cap; those powers ",
            "are based on the converged replications.", call. = FALSE)
  }

  mc_se <- sqrt(res$composite * (1 - res$composite) / res$converged)
  out <- data.frame(
    term = c(size_term, "composite_power", "composite_power_mc_se",
             paste0("power_", labels),
             paste0("population_", labels),
             "alpha_level", "replications", "converged_replications"),
    value = c(size, res$composite, mc_se,
              unname(res$marginal),
              setup$theta,
              alpha_level, G, res$converged),
    stringsAsFactors = FALSE)
  if (!is.null(desired_power)) {
    out <- rbind(out,
                 data.frame(term = "desired_power", value = desired_power))
  }
  class(out) <- c("dmar_ss_power", "data.frame")
  attr(out, "composite_terms") <- labels
  .as_dmar_tbl(out)
}
