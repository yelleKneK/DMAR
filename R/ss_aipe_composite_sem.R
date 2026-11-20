#' Sample Size for Accurate Estimation of a Set of SEM Parameters
#'
#' Determine the necessary sample size for a structural equation model study
#' so that the confidence interval for every parameter of interest is
#' sufficiently narrow in the same study, or, given a sample size, return how
#' narrow the set of intervals can be expected to be. This is the accuracy in
#' parameter estimation (AIPE) counterpart of
#' \code{\link{ss_power_composite_sem}}: where that function plans for every
#' parameter to be statistically significant jointly, this one plans for
#' every parameter to be estimated with a confidence interval no wider than
#' desired, the goal when the research questions concern the magnitudes of
#' the effects rather than their existence. The parameters of interest are
#' any labeled parameters of a lavaan analysis model, structural paths,
#' loadings, covariances, or quantities defined with \code{:=} such as an
#' indirect effect, and any subset of them can make up the set.
#'
#' @param model A single character string giving the free analysis model in
#'   lavaan model syntax (see \code{\link[lavaan]{model.syntax}}), the model
#'   that would be fit to the data. Each parameter of interest must carry a
#'   parameter label so it can be referred to by name, for example
#'   \code{"f2 ~ b*f1"} labels the structural path \code{b}, and
#'   \code{"ab := a*b"} defines an indirect effect from the labeled paths
#'   \code{a} and \code{b}.
#' @param true_cov_matrix Population covariance matrix of the observed
#'   variables, with row and column names matching the observed variables in
#'   \code{model}. It is typically obtained from a fully fixed population
#'   model via \code{\link{cov_sem}}. Supply exactly one of
#'   \code{true_cov_matrix} or \code{true_model}.
#' @param true_model A single character string giving the population model in
#'   lavaan model syntax with every parameter fixed to its population value,
#'   from which \code{\link{cov_sem}} derives \code{true_cov_matrix} (and the
#'   population means, when the model has a mean structure). The fixed
#'   values are values the researcher posits (from theory, prior studies, or
#'   pilot data), never sample estimates. Supply exactly one of
#'   \code{true_cov_matrix} or \code{true_model}.
#' @param true_mean Optional population means of the observed variables,
#'   used with \code{true_cov_matrix}: a named numeric vector with one entry
#'   per observed variable, or an unnamed vector in the row order of
#'   \code{true_cov_matrix}. The default \code{NULL} is zero means. Means
#'   matter only when the analysis model has a mean structure (an intercept
#'   term such as \code{s ~ 1}, as in a latent growth curve model); when
#'   \code{true_model} is supplied its mean structure provides the means and
#'   \code{true_mean} must not also be given.
#' @param parameters Character vector of the parameter labels that make up
#'   the set. The default \code{NULL} uses every user-labeled parameter in
#'   \code{model}, in order of appearance, so labeling exactly the parameters
#'   of interest is the simplest way to state the set.
#' @param width The desired full confidence interval width for each
#'   parameter of interest: a single value applied to every parameter, or a
#'   named numeric vector with one entry per parameter label. An unnamed
#'   vector of several widths is not accepted, so a width can never silently
#'   attach to the wrong parameter.
#' @param conf_level Confidence level of each interval (default 0.95).
#' @param assurance The desired probability that a single study yields
#'   confidence intervals no wider than desired for every parameter of
#'   interest simultaneously (a value in [0.5, 1)), or \code{NULL} (the
#'   default) to plan against the expected widths instead; see Details.
#' @param specified_N Total sample size \emph{N}; if supplied, the realized
#'   interval widths at that \emph{N} are summarized rather than a sample
#'   size planned.
#' @param G Number of converged Monte Carlo replications per evaluated sample
#'   size (default 1000). The simulation error of each estimated proportion
#'   is about \eqn{\sqrt{p(1 - p)/G}}; raise \code{G} for a sharper answer.
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
#'   Every name given here must be an argument or an option of \pkg{lavaan}
#'   (see \code{\link[lavaan]{lavOptions}}); any other name, and any
#'   argument without a name, is refused before a model is fitted, with a
#'   condition of class \code{dmar_argument_invalid} that gives the name
#'   as it was typed
#'
#' @details
#' AIPE planning for a single targeted SEM parameter is available in closed
#' form (Lai & Kelley, 2011; \code{\link{ss_aipe_sem_path}}), but most
#' studies estimate several effects and report all of them; a design is only
#' as informative as its widest interval of interest. This function plans for
#' the set by a priori Monte Carlo simulation (Muthén &
#' Muthén, 2002; Maxwell, Kelley, & Rausch, 2008): for a candidate
#' \emph{N}, \code{G} data sets are drawn from the multivariate normal
#' population with covariance matrix \code{true_cov_matrix}, the analysis
#' model is fit to each, and each parameter's Wald confidence interval width,
#' twice
#' \eqn{z_{1 - \alpha/2}} times its standard error, is recorded. Because the
#' estimates share one fitted model, the widths are dependent; the simulation
#' reflects that dependence exactly, at the stated \emph{N}, with no
#' asymptotic shortcut.
#'
#' Two planning criteria are available. With \code{assurance = NULL} the
#' necessary sample size is the smallest \emph{N} at which the mean simulated
#' width of every parameter's interval is at or below its desired width, the
#' expected-width criterion of the AIPE framework applied to each member of
#' the set. Widths vary from sample to sample around their means, so each
#' interval separately lands at or below its desired width in roughly half
#' of the realizations. That is a statement about one interval at a time,
#' not about the set: the probability that \emph{every} interval is narrow
#' enough at once falls well below one half as soon as more than one
#' parameter binds, and falls further the more parameters are targeted and
#' the more weakly their widths move together. Planning the whole set to a
#' stated probability is exactly what \code{assurance} is for. Supplying
#' \code{assurance} plans against the
#' joint event instead: the smallest \emph{N} at which the proportion of
#' replications where every interval is simultaneously within its desired
#' width reaches the assurance. The joint event is contained in each marginal
#' event, so its probability is at most the smallest per-parameter
#' proportion, and the \code{pct_ci_less_width_<label>} rows show which
#' parameter binds the design.
#'
#' When \code{specified_N} is \code{NULL} the search starts at the largest of the
#' per-parameter closed-form sample sizes (the no-assurance approximation
#' \code{\link{ss_aipe_sem_path}} uses, computed from the asymptotic
#' variances before any simulation), brackets the crossing geometrically, and
#' bisects to adjacent integers, each candidate evaluated with its own
#' \code{G} replications. A planning call therefore fits the analysis model
#' several thousand times at the default \code{G}, and even at the smallest
#' admissible \code{G} the search runs for several seconds, so the example
#' below evaluates a stated \code{specified_N}, which is the cheap half of the method.
#' A planning call is the same call with \code{specified_N} left out: supplying
#' \code{assurance} plans against the joint event, and leaving
#' \code{assurance} out as well plans against the expected widths; either way
#' the first row of the result is \code{necessary_N} rather than
#' \code{specified_N}. The vignette
#' \code{vignette("composite_sem_planning", package = "DMAR")} works through
#' both planning calls for a mediation model and a latent growth curve
#' model, with reference values computed at \code{G = 10000}.
#'
#' Each reported proportion carries a simulation standard error of about
#' \eqn{\sqrt{p(1 - p)/G}}, and the necessary sample size inherits that
#' uncertainty; raising \code{G} narrows it, and reporting the seed makes a
#' plan reproducible.
#'
#' @note
#' A replication whose fit does not converge, or converges without a usable
#' standard error for some parameter of interest, is discarded and fresh data
#' are drawn, up to \code{20 * G} attempts per evaluated sample size; the
#' reported summaries condition on convergence. A replication whose
#' likelihood has no maximum (an error variance running off below zero
#' while the likelihood is flat or still rising; see
#' \code{\link{fa_bifactor}}), which \pkg{lavaan} reports as converged where
#' its optimizer stopped, is discarded the same way, and an analysis model
#' with no maximum at the population moments stops the function with an
#' error of class \code{dmar_no_maximum} that prints the \pkg{lavaan} code
#' of that fit. When fewer than \code{G} replications converge within the
#' cap, a single warning is issued and the summary is based on the
#' converged replications (their count is the \code{converged_replications}
#' row).
#'
#' Because the planner itself is a Monte Carlo study, it has no separate
#' \code{_sensitivity} sibling; to study misspecification of the population
#' values, rerun the planner with the alternative \code{true_cov_matrix} or
#' \code{true_model} values under consideration and compare the plans.
#'
#' @return A \code{data.frame} (a \code{dmar_tbl}) with \code{term} and
#'   \code{value} columns: the \code{necessary_N} (or supplied
#'   \code{specified_N}), the \code{composite_assurance} (the proportion of
#'   replications in which every interval was simultaneously within its
#'   desired width, reported under both criteria), then for each parameter
#'   its \code{mean_width_<label>}, its marginal
#'   \code{pct_ci_less_width_<label>} proportion, its
#'   \code{width_<label>}, and its purported
#'   \code{true_<label>} value under the analysis model, followed by
#'   \code{conf_level}, the requested \code{replications}, the
#'   \code{converged_replications} the summary is based on, and, when
#'   supplied, the \code{assurance}.
#'
#' @references
#' Lai, K., & Kelley, K. (2011). Accuracy in parameter estimation for
#'   targeted effects in structural equation modeling: Sample size
#'   planning for narrow confidence intervals.
#'   \emph{Psychological Methods, 16}(2), 127--148. \doi{10.1037/a0021764}
#'
#' Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). \emph{Designing
#'   experiments and analyzing data: A model comparison perspective}
#'   (4th ed.). Routledge.
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
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso \code{\link{ss_power_composite_sem}} for the same set of
#'   parameters planned for joint statistical significance;
#'   \code{\link{cov_sem}} for deriving \code{true_cov_matrix} from a fully
#'   fixed
#'   population model; \code{\link{ss_aipe_sem_path}} and
#'   \code{\link{ss_aipe_sem_path_sensitivity}} for a single targeted path.
#'
#' @examples
#' # A mediation model whose research questions concern the magnitudes of
#' # both individual paths and the indirect effect. The population model
#' # fixes every parameter to its purported population value.
#' true_model <- "
#'   f1 =~ 1*y1 + 0.8*y2 + 0.8*y3
#'   f2 =~ 1*y4 + 0.8*y5 + 0.8*y6
#'   f3 =~ 1*y7 + 0.8*y8 + 0.8*y9
#'   f2 ~ 0.4*f1
#'   f3 ~ 0.5*f2 + 0.2*f1
#'   f1 ~~ 1*f1
#'   f2 ~~ 0.84*f2
#'   f3 ~~ 0.7*f3
#'   y1 ~~ 0.5*y1; y2 ~~ 0.5*y2; y3 ~~ 0.5*y3
#'   y4 ~~ 0.5*y4; y5 ~~ 0.5*y5; y6 ~~ 0.5*y6
#'   y7 ~~ 0.5*y7; y8 ~~ 0.5*y8; y9 ~~ 0.5*y9
#' "
#'
#' # The analysis model labels the two paths and defines the indirect
#' # effect; all three make up the set of interest.
#' analysis_model <- "
#'   f1 =~ y1 + y2 + y3
#'   f2 =~ y4 + y5 + y6
#'   f3 =~ y7 + y8 + y9
#'   f2 ~ a*f1
#'   f3 ~ b*f2 + cp*f1
#'   ab := a*b
#' "
#'
#' # Realized interval widths at N = 200, with the indirect effect held to a
#' # narrower interval than the paths through a named vector of widths. Each
#' # interval lands within its desired width in most of the replications, yet
#' # all three do so together in far fewer of them: that joint proportion,
#' # reported as composite_assurance, is what a design of this kind has to be
#' # planned against. G = 10 keeps the example quick; a reported plan
#' # deserves the default G = 1000 or more.
#' ss_aipe_composite_sem(model = analysis_model, true_model = true_model,
#'                       parameters = c("a", "b", "ab"),
#'                       width = c(a = 0.35, b = 0.40, ab = 0.25),
#'                       specified_N = 200, G = 10, seed = 113)
#'
#' @keywords design multivariate
#'
#' @family AIPE sample size planning
#'
#' @export
ss_aipe_composite_sem <- function(model, true_cov_matrix = NULL,
                                  true_model = NULL,
                                  true_mean = NULL, parameters = NULL,
                                  width, conf_level = 0.95,
                                  assurance = NULL, specified_N = NULL, G = 1000,
                                  seed = NULL, ...) {
  .dmar_check_lavaan_dots(..., fn = "ss_aipe_composite_sem")
  if (!requireNamespace("lavaan", quietly = TRUE)) {
    .dmar_refuse(
      "The package 'lavaan' is needed; please install it and try again.",
      "package_required", package = "lavaan")
  }
  if (!is.numeric(conf_level) || length(conf_level) != 1L ||
      is.na(conf_level) || conf_level <= 0 || conf_level >= 1) {
    .dmar_refuse_argument(
      "'conf_level' must be a single number in (0, 1).", arg = "conf_level")
  }
  if (!is.null(assurance) &&
      (!is.numeric(assurance) || length(assurance) != 1L ||
       is.na(assurance) || assurance < 0.5 || assurance >= 1)) {
    .dmar_refuse_argument(
      "'assurance' must be NULL or a single value in [0.5, 1).",
      arg = "assurance")
  }
  G <- .dmar_check_count(G, "G", min = 10L)
  # The joint proportion has the same 1/G granularity as the composite power
  # in ss_power_composite_sem(); see .composite_sem_check_resolution().
  if (!is.null(assurance)) {
    .composite_sem_check_resolution(assurance, G, "assurance")
  }

  # A supplied seed governs the Monte Carlo draws of this call only; the
  # caller's generator state is restored when the function exits.
  .dmar_local_seed(seed)

  setup <- .composite_sem_setup(model, true_cov_matrix, true_model, parameters,
                                true_mean = true_mean, ...,
                                fn = "ss_aipe_composite_sem")
  labels <- setup$labels
  k <- length(labels)

  # One desired width per parameter: a scalar recycled to all of them, or a
  # named vector matched by label so a width can never silently attach to
  # the wrong parameter.
  if (!is.numeric(width) || anyNA(width) ||
      any(width <= 0)) {
    .dmar_refuse_argument(c(
      "'width' must be a positive number or a named numeric ",
      "vector of positive widths."), arg = "width")
  }
  if (length(width) == 1L && is.null(names(width))) {
    omega <- stats::setNames(rep(width, k), labels)
  } else {
    if (is.null(names(width)) ||
        !setequal(names(width), labels) ||
        anyDuplicated(names(width))) {
      .dmar_refuse_argument(c(
        "A 'width' vector must be named, with exactly one entry ",
        "per parameter of interest (", paste(labels, collapse = ", "),
        ")."), arg = "width")
    }
    omega <- width[labels]
  }

  N_min <- max(10L, setup$n_obs_vars + 2L)
  z <- stats::qnorm(1 - (1 - conf_level) / 2)

  # Evaluations where fewer than G replications converged within the attempts
  # cap are counted in a state environment that evaluate() updates, so a
  # single summary warning is issued at the end rather than one per
  # evaluation.
  state <- new.env(parent = emptyenv())
  state$short_evals <- 0L
  evaluate <- function(n) {
    mc <- .composite_sem_mc(model, setup$true_cov_matrix, setup$true_mean,
                            labels, n, G, ...)
    if (mc$converged == 0L) {
      .dmar_refuse_argument(c(
        "No replication converged at N = ", n, " within ", mc$attempts,
        " attempts; the model cannot be fit at so small a sample size. ",
        "Reconsider the model or the candidate sample size."), arg = NULL)
    }
    if (mc$converged < G) state$short_evals <- state$short_evals + 1L
    widths <- 2 * z * mc$se
    within <- sweep(widths, 2L, omega, "<=")
    list(mean_width = colMeans(widths),
         marginal = colMeans(within),
         joint = mean(rowSums(within) == k),
         converged = mc$converged)
  }
  met <- if (is.null(assurance)) {
    function(e) all(e$mean_width <= omega)
  } else {
    function(e) e$joint >= assurance
  }

  if (!is.null(specified_N)) {
    if (!is.numeric(specified_N) || length(specified_N) != 1L ||
        is.na(specified_N) || specified_N != round(specified_N) ||
        specified_N < N_min) {
      .dmar_refuse_argument(c(
        "'specified_N' must be a single whole number of at least ", N_min,
        " for this model."), arg = "specified_N")
    }
    size_term <- "specified_N"
    size <- as.integer(specified_N)
    res <- evaluate(size)
  } else {
    N_start <- .composite_sem_start_width(setup$h, omega, conf_level, N_min)
    searched <- .composite_sem_search(evaluate, met, N_start, N_min)
    size_term <- "necessary_N"
    size <- searched$N
    res <- searched$eval
  }

  if (state$short_evals > 0L) {
    .dmar_warn(c(
      "In ", state$short_evals, " Monte Carlo evaluation",
      if (state$short_evals > 1L) "s" else "", " fewer than G = ", G,
      " replications converged within the attempts cap; those ",
      "summaries are based on the converged replications."), "result_qualified")
  }

  out <- data.frame(
    term = c(size_term, "composite_assurance",
             paste0("mean_width_", labels),
             paste0("pct_ci_less_width_", labels),
             paste0("width_", labels),
             paste0("true_", labels),
             "conf_level", "replications", "converged_replications"),
    value = c(size, res$joint,
              unname(res$mean_width),
              unname(res$marginal),
              unname(omega),
              setup$theta,
              conf_level, G, res$converged),
    stringsAsFactors = FALSE)
  if (!is.null(assurance)) {
    out <- rbind(out, data.frame(term = "assurance", value = assurance))
  }
  attr(out, "composite_terms") <- labels
  .as_dmar_tbl(out, conf_level = conf_level, subclass = "dmar_ss_aipe")
}
