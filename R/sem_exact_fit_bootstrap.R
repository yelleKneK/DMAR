#' Bootstrap Test of Exact Fit for a Structural Equation Model With Complete or Incomplete Data
#'
#' Gives the chi square test of exact fit of a structural equation model a
#' reference distribution that does not assume normal data: the model
#' based bootstrap of Bollen and Stine (1992), in which the data are first
#' transformed so that the model holds in them exactly, as estimated, and
#' the statistic is then recomputed on \code{B} samples drawn with
#' replacement from the transformed data. The share of bootstrap statistics
#' at least as large as the observed one is the bootstrap \emph{p}-value.
#' For incomplete data fitted by full information maximum likelihood (FIML)
#' the transformation is extended so that the FIML estimates of the
#' saturated model on the transformed data are the model implied moments
#' (Enders, 2002; Savalei & Yuan, 2009), and each bootstrap sample keeps
#' the missing values of the rows it draws, so that the model is fitted to
#' it by FIML as it was to the data. Requires \pkg{lavaan} (in Suggests).
#'
#' @param fit A fitted \pkg{lavaan} model: a maximum likelihood fit
#'   (estimator ML, MLM, MLMV, or MLR) of raw continuous data, complete or
#'   incomplete under \code{missing = "ml"} (with auxiliary variables
#'   through lavaan's \code{aux} or without), in one or several groups,
#'   that converged to an admissible solution. See Details for what is
#'   refused.
#' @param B Number of bootstrap samples, a single whole number of at
#'   least 100. Defaults to \code{1000}: every sample refits the model. A
#'   Monte Carlo \emph{p}-value is run only from a request of at least 100
#'   samples, at which its standard error is at most 0.05; a smaller
#'   \code{B} is
#'   refused before anything is computed (class \code{dmar_count_invalid}).
#'   See \code{\link{dmar_bootstrap}}.
#' @param transformation For incomplete data, how the transformation
#'   makes the saturated FIML estimates equal the model implied moments:
#'   \code{"common"} (the default) finds one matrix that transforms every
#'   row on its observed variables, symmetric on the scale of the
#'   variables divided by their saturated standard deviations, so that the
#'   FIML estimating equations hold summed over the missing data patterns;
#'   \code{"pattern"} transforms the rows of each missing data pattern on
#'   their own, with the pattern's own covariance matrix and means, so that
#'   each pattern's covariance matrix and means equal the model implied
#'   ones exactly, which needs every pattern to have more rows than
#'   observed variables. Both are solved on that standardized scale, so the
#'   result does not depend on the units of the variables. See Details.
#'   Ignored with complete data.
#' @param improper Whether a bootstrap sample whose fit converged to an
#'   inadmissible solution counts in the \emph{p}-value: \code{"keep"} (the
#'   default) or \code{"drop"}. A solution is inadmissible when lavaan's
#'   post-estimation check finds a negative variance estimate or a
#'   covariance matrix that is not positive definite, of the latent
#'   variables, of the residuals of the observed variables, or of the
#'   residuals of the structural part of the model; with auxiliary
#'   variables the check is made on the model's own parameters (see What Is
#'   Counted). Either way the row \code{B_improper} reports how many there
#'   were.
#' @param conf_level Confidence level of the interval for the Monte Carlo
#'   error of the bootstrap \emph{p}-value. Defaults to 0.95.
#' @param seed Random number seed for the bootstrap samples. Defaults to
#'   \code{NULL}, under which the samples are drawn from, and advance, the
#'   current random number stream; supply an integer for a reproducible
#'   result, which is set for the call only, the caller's random number
#'   state being restored when it returns.
#' @param ncpus The number of processes over which the bootstrap refits are
#'   spread. The default \code{1} runs them in this R session; a larger
#'   number uses forked workers or, on Windows, a socket cluster. Every
#'   bootstrap sample is drawn before any refit runs, so a given
#'   \code{seed} gives the same result whatever \code{ncpus} is. A session
#'   that forks workers and also starts background processes through
#'   \pkg{callr} or \pkg{processx} (a Shiny app, a job runner) needs the
#'   environment variable \code{PROCESSX_NOTIFY_OLD_SIGCHLD} set to a
#'   nonempty value before \pkg{processx} is loaded: \pkg{processx}
#'   otherwise replaces the signal handler through which \pkg{parallel}
#'   learns that a forked worker has finished, and \pkg{parallel} then
#'   reports at exit that it could not terminate its workers
#'   (\pkg{processx} documents the interaction and the switch). The
#'   results are not affected.
#'
#' @details
#' \subsection{The Model Based Bootstrap}{
#' The chi square statistic \eqn{T} of a maximum likelihood fit follows its
#' chi square distribution under the model only for normal data in large
#' samples. Bollen and Stine (1992) proposed to find its distribution under
#' the model by bootstrap instead: resampling the data as they are would
#' reproduce the misfit of the model as well as the sampling variation, so
#' the data are first transformed to
#' \deqn{Z = (Y - \bar{Y}) S^{-1/2} \hat{\Sigma}^{1/2}}{Z = (Y - Ybar) S^(-1/2) Sigma_hat^(1/2)}
#' (Enders, 2002, Eq. 1), where \eqn{S} is the sample covariance matrix,
#' \eqn{\hat{\Sigma}} the covariance matrix the fitted model implies, and
#' the square roots are symmetric. The transformed data have covariance
#' matrix exactly \eqn{\hat{\Sigma}}, so the model fits them exactly, and
#' their distribution otherwise keeps the shape of the data. Each bootstrap
#' sample draws rows of \eqn{Z} with replacement, within each group as
#' many as the group has (each group is transformed with its own
#' moments), the model is fitted to it, and its statistic \eqn{T^*} is
#' recorded. The bootstrap
#' \emph{p}-value is the share of the \eqn{T^*} at least as large as the
#' observed \eqn{T} (row \code{p_bootstrap}); a small value says that the
#' observed statistic is unusual for data in which the model holds and
#' have this distribution. Under a mean structure the model implied means
#' are added back, so the means hold too. The statistic is the standard
#' (uncorrected) likelihood ratio statistic in every case: the bootstrap
#' replaces the correction a scaled statistic applies, and a fit by MLM or
#' MLR contributes its estimates and its standard statistic. This is the
#' bootstrap \code{lavaan::bootstrapLavaan(type = "bollen.stine")} runs
#' for complete data, with the same transformation and resamples: the
#' \code{seed} given here, passed to lavaan as its \code{iseed} (with
#' lavaan run serially), gives the same \eqn{T^*}. lavaan called without
#' \code{iseed} draws its own seed from the caller's random number stream,
#' so a \code{set.seed()} before the two calls does not make them agree.
#' }
#'
#' \subsection{Incomplete Data}{
#' With missing values the transformation cannot be applied row by row,
#' since a row observes only some of the variables. Enders (2002, Eq. 4)
#' first extended it by transforming each row with the submatrices, for
#' the variables the row observes, of the saturated and the model implied
#' covariance matrices estimated by FIML from the whole sample, and adding
#' back the model implied means; the model then does not hold exactly in
#' the transformed data, since the rows of a pattern keep their own
#' covariance matrix, and that transformation is not offered here. Savalei
#' and Yuan (2009) gave two transformations under which the FIML estimates
#' of the saturated model on the transformed data are exactly
#' \eqn{\hat{\mu}} and \eqn{\hat{\Sigma}}, so that the model fits the
#' transformed data exactly under FIML, as the complete data transformation
#' makes it fit exactly under complete data. The \code{"pattern"}
#' transformation (their first) transforms the rows of each missing data
#' pattern with the pattern's own covariance matrix and means, so each
#' pattern's rows take exactly the model implied moments of the variables
#' they observe; it needs every pattern to have more rows than observed
#' variables, since the pattern's covariance matrix is inverted, and a
#' pattern that falls short is refused with its index and counts. The
#' \code{"common"} transformation (their second, and the default) centers
#' each pattern at its own means, transforms every row with one matrix,
#' symmetric on the standardized scale described below, chosen so that
#' the covariance estimating equations of the saturated FIML model, summed
#' over the patterns, hold at \eqn{\hat{\Sigma}}, and shifts each pattern
#' to \eqn{\hat{\mu}}; it does not need large patterns. DMAR finds the
#' matrix by Newton's method on the standardized variables, starting from
#' the complete data transformation applied to the EM moments, and
#' records the number of steps and the largest estimating equation at the
#' solution in the \code{transformation} attribute, where the pattern
#' transformation records the same check computed afresh. Both
#' transformations are solved on the variables divided by their saturated
#' standard deviations, and the transformed values are multiplied back:
#' the matrix of the common transformation and the square roots of the
#' pattern transformation are symmetric on that standardized scale. The
#' bootstrap therefore does not depend on the units of the variables
#' (rescaling an item rescales its transformed values and changes nothing
#' else). Symmetry belongs to the scale on which it is imposed, so the
#' same transformations solved in raw units, as
#' \code{semTools::bsBootMiss()} solves the common one (its transformation
#' 2), give transformed data that differ element by element from these and
#' depend on the units; both versions make the saturated FIML estimates
#' exactly \eqn{\hat{\mu}} and \eqn{\hat{\Sigma}}, so the model holds
#' exactly in either. Each bootstrap sample draws
#' rows of the transformed data with their missing values, so it has the
#' missing data patterns of the data in bootstrap proportions, and the
#' model is fitted to it by FIML. In Enders's (2002) simulation of his
#' transformation, with 10 to 20 percent of the values missing completely
#' at random from nine nonnormal indicators, the bootstrap held the Type I
#' error rate near .05 at \eqn{N} = 250 and 500, where the uncorrected
#' FIML statistic rejected far too often, and was conservative at \eqn{N}
#' = 100, the more so the more nonnormal the data and the more missing
#' values; he suggests a sample of about 200 or more with incomplete data.
#' The calibration of the two transformations offered here on the current
#' code is planned, not done.
#' }
#'
#' \subsection{What Is Counted}{
#' A bootstrap sample whose fit did not converge, raised an error, or
#' returned no finite statistic is set aside, and one that converged to
#' an inadmissible solution is counted or set aside as \code{improper}
#' says. That includes a sample whose likelihood has no maximum, an error
#' variance running off below zero while the likelihood is flat or still
#' rising (see \code{\link{fa_bifactor}}): its estimates are wherever
#' \pkg{lavaan}'s optimizer stopped on that ridge, but its chi square is
#' not, since the likelihood levels off along the ridge and the optimizer
#' stops only where its steps no longer gain, which puts the statistic
#' within the optimizer's tolerance of its limit. A sample's solution is
#' judged as the fit's is (see What Is Refused); with auxiliary variables
#' that is the check of the model's own parameters, since \pkg{lavaan}'s
#' own check would also count as improper the samples in which only the
#' residual covariances of the auxiliary variables are indefinite, which
#' with an auxiliary variable related to a factor can be many samples
#' whose solution is proper. The rows \code{B},
#' \code{B_used}, \code{B_improper}, and \code{B_failed} count them, and a
#' single warning (class \code{dmar_bootstrap_replications_dropped})
#' reports the samples set aside. The \emph{p}-value is formed from the
#' \code{B_used} statistics; when fewer than half of the attempted samples
#' were used, it is \code{NA} with a warning of class
#' \code{dmar_bootstrap_insufficient}, since a \emph{p}-value from the
#' samples that happened to converge would describe a selected distribution.
#' The 100 required are held to the number requested, \code{B}, and not to
#' \code{B_used}: a request that met the count and then had samples set
#' aside still reports its \emph{p}-value, the warning says how many were
#' used, and the Wilson interval shows the Monte Carlo error from that
#' number. Where refits can fail, request more than the minimum.
#' \code{options(dmar.bootstrap_exploratory =
#' TRUE)} lets a smaller \code{B} (2 or more) run for a first look or a
#' demonstration; every such call gives a warning of class
#' \code{dmar_bootstrap_exploratory}, and the line under the table says
#' so (see \code{\link{dmar_bootstrap}}). The rows \code{p_bootstrap_ci_lower}
#' and \code{p_bootstrap_ci_upper} are the Wilson interval for the Monte
#' Carlo error of the \emph{p}-value at \code{conf_level}, the uncertainty
#' that comes from using \code{B_used} samples rather than infinitely many
#' (\code{\link{ci_proportion}}); it is not an interval for anything about
#' the model. The row \code{bootstrap_mean} is the mean of the used
#' statistics, the bootstrap estimate of the expected statistic under the
#' model, which the Bollen-Stine adjustments of fit indices use in place
#' of the degrees of freedom.
#' }
#'
#' \subsection{What Is Refused}{
#' A fit is refused when it did not converge (an error of class
#' \code{dmar_not_converged}) or has an inadmissible solution (a refusal
#' of class \code{dmar_inadmissible} that names the negative variance or
#' the matrix at fault), as for every fit measure of the SEM family (see
#' \code{\link{sem_fit_measures}}, which explains the check and how it
#' treats auxiliary variables). It is refused with an error of class
#' \code{dmar_exact_fit_bootstrap_scope}, whose \code{reason} element
#' names the case (see \code{\link{dmar_metadata}}), when it was made from
#' summary statistics rather than raw data (there are no rows to
#' transform); when it has ordered indicators or was fitted by an estimator
#' other than maximum likelihood (the transformation is derived for the
#' normal theory fit function); when it is a two-level model, has a
#' cluster variable (the resamples draw rows independently, which the
#' refit could not undo) or sampling weights, uses two-stage estimation
#' (whose statistic is not the likelihood ratio statistic the bootstrap
#' refits), or conditions on observed exogenous variables (\code{fixed.x =
#' TRUE} with such variables, or \code{conditional.x = TRUE}: the
#' transformation moves every variable, and the conditional likelihood
#' holds some fixed; refit with \code{fixed.x = FALSE}, and with
#' \code{conditional.x = FALSE} too, since \pkg{lavaan} does not take the
#' first beside \code{conditional.x = TRUE}); and when its analysis data
#' keep missing values without FIML (as under pairwise deletion) or it has
#' no standard test statistic. A missing
#' data pattern with no more rows than observed variables is refused for
#' the pattern transformation with an error of class
#' \code{dmar_pattern_too_small}, which carries the pattern, its counts,
#' and its group, and a transformation that is undefined for the data
#' (a moment that is not positive definite, a Newton iteration that does
#' not converge) with one of class \code{dmar_transformation_undefined}.
#' These three classes, and \code{dmar_inadmissible}, inherit
#' \code{dmar_refusal}; \code{dmar_not_converged} is caught by its own
#' name. The checks of \code{B}, \code{conf_level}, and \code{ncpus} are
#' argument checks and give plain errors. A fit whose cases
#' were removed by listwise deletion is bootstrapped on its complete cases,
#' as it was fitted, with a warning. A fit by MLM, MLMV, or MLR is accepted:
#' its estimates are the ML estimates, and the bootstrap refits and
#' resamples its standard statistic.
#' }
#'
#' @return A \code{dmar_tbl} (class \code{dmar_sem_fit}) with columns
#'   \code{term} and \code{value} and, always in this order, the rows
#'   \code{chi_square}, \code{df}, \code{p_chi_square} (the asymptotic
#'   \emph{p}-value of the standard statistic), \code{p_bootstrap},
#'   \code{p_bootstrap_ci_lower}, \code{p_bootstrap_ci_upper},
#'   \code{bootstrap_mean}, \code{B}, \code{B_used}, \code{B_improper},
#'   \code{B_failed}, \code{N}, and \code{n_groups}. The attributes
#'   \code{B}, \code{B_used}, \code{n_dropped}, \code{B_required} (100),
#'   and \code{bootstrap_status} record the samples in the names every
#'   bootstrap result of the package carries. The attribute
#'   \code{bootstrap} is a data frame with one row per bootstrap sample,
#'   its \code{chi_square} (\code{NA} when it was set aside) and its
#'   \code{status} (\code{proper}, \code{improper}, \code{not_converged},
#'   or \code{error}); the attribute \code{transformation} records, per
#'   group, the transformation applied, the number of missing data
#'   patterns, the Newton steps and the largest estimating equation at the
#'   solution (on the standardized scale), and the largest difference
#'   between the transformed means and the model implied ones; and the
#'   \code{sem_fit} attribute records the statistic and its degrees of
#'   freedom, the transformation requested and the ones applied, the
#'   counts, and the lines of the printed footer.
#'
#' @references
#' Bollen, K. A., & Stine, R. A. (1992). Bootstrapping goodness-of-fit
#'   measures in structural equation models. \emph{Sociological Methods &
#'   Research, 21}(2), 205--229. \doi{10.1177/0049124192021002004}
#'
#' Enders, C. K. (2002). Applying the Bollen-Stine bootstrap for
#'   goodness-of-fit measures to structural equation models with missing
#'   data. \emph{Multivariate Behavioral Research, 37}(3), 359--377.
#'   \doi{10.1207/S15327906MBR3703_3}
#'
#' Savalei, V., & Yuan, K.-H. (2009). On the model-based bootstrap with
#'   missing data: Obtaining a \emph{p}-value for a test of exact fit.
#'   \emph{Multivariate Behavioral Research, 44}(6), 741--763.
#'   \doi{10.1080/00273170903333590}
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso \code{\link{sem_fit_measures}} for the asymptotic test of
#'   exact fit beside the other fit measures, and
#'   \code{\link{equivalence_rmsea}} for the test of not close fit.
#'
#' @family multivariate and latent variable methods
#'
#' @keywords multivariate htest
#'
#' @examples
#' # Two correlated factors for six of the holzinger_swineford tests, with
#' # the cube and sentence scores of about a third of the students missing
#' # at random given the general information score, which the FIML fit
#' # takes as an auxiliary variable. lavaan notes that it adds the
#' # auxiliary variable and warns that the residual covariance matrix,
#' # which includes it, is not positive definite; the model's own solution
#' # is proper (see What Is Refused). B = 10 keeps the example quick and
#' # is below the 100 a reported p-value needs; the option lets it run
#' # and the warning says so. Use the default of 1000 or more in an
#' # analysis.
#' d <- holzinger_swineford[c("t1_visual_perception", "t2_cubes", "t4_lozenges",
#'                            "t6_paragraph_comprehension", "t7_sentence",
#'                            "t9_word_meaning", "t5_general_information")]
#' set.seed(113)
#' drop <- runif(nrow(d)) < plogis((d$t5_general_information - 30) / 10) - 0.35
#' d[drop, c("t2_cubes", "t7_sentence")] <- NA
#' fit <- lavaan::cfa(
#'   "visual  =~ t1_visual_perception + t2_cubes + t4_lozenges
#'    textual =~ t6_paragraph_comprehension + t7_sentence + t9_word_meaning",
#'   data = d, missing = "ml", std.lv = TRUE, aux = "t5_general_information")
#' old <- options(dmar.bootstrap_exploratory = TRUE)
#' sem_exact_fit_bootstrap(fit, B = 10, seed = 113)
#' options(old)
#'
#' @export
sem_exact_fit_bootstrap <- function(fit, B = 1000,
                                    transformation = c("common", "pattern"),
                                    improper = c("keep", "drop"),
                                    conf_level = 0.95, seed = NULL, ncpus = 1L) {
  transformation <- match.arg(transformation)
  improper <- match.arg(improper)
  # The fit is the user's and is not refit to its own data, but the
  # statistic and the transformation are read from it. Whether its
  # likelihood has a maximum (.lavaan_no_maximum()) is not asked: a fit
  # can have none only with an error variance below zero, and such a fit
  # fails the check of admissibility (.sem_admissible(), with or without
  # auxiliary variables), which .sem_fit_check() refuses.
  # Nor is it asked of the bootstrap samples, whose chi square does not
  # depend on where the optimizer stopped on a ridge (.bs_refit_one()).
  options <- .sem_fit_check(fit)
  # A Monte Carlo p-value is reported only from at least 100 samples
  # (?dmar_bootstrap), at which its standard error is at most 0.05.
  # Through the first versions of 1.1.0 the minimum was 2.
  B <- .dmar_check_B(B, kind = "p")
  ncpus <- .check_ncpus(ncpus)
  .sem_probability(conf_level, "conf_level")
  .dmar_local_seed(seed)

  # ---- scope ----
  # A fit outside the scope of the transformation is refused with an error
  # of class dmar_exact_fit_bootstrap_scope (a dmar_refusal) whose `reason`
  # names the check that failed. The checks of B, conf_level, and ncpus
  # above are argument checks and stop plainly.
  refuse <- function(message, reason, ...) {
    .dmar_refuse(message, "exact_fit_bootstrap_scope", reason = reason, ...)
  }
  if (is.null(lavaan::lavInspect(fit, "case.idx"))) {
    refuse(c("sem_exact_fit_bootstrap() needs a model fitted to raw data: a ",
             "fit to summary statistics has no rows to transform and ",
             "resample."), "summary_data")
  }
  if (length(lavaan::lavNames(fit, "ov.ord"))) {
    refuse(c("sem_exact_fit_bootstrap() does not take ordered indicators: the ",
             "model based transformation acts on continuous variables."),
           "ordered", variables = lavaan::lavNames(fit, "ov.ord"))
  }
  if (!identical(options$estimator, "ML")) {
    refuse(c("sem_exact_fit_bootstrap() requires a maximum likelihood fit ",
             "(estimator ML, MLM, MLMV, or MLR): the model based ",
             "transformation is derived for the normal theory fit function; ",
             "this fit used ", options$estimator, "."),
           "estimator", estimator = options$estimator)
  }
  if (lavaan::lavInspect(fit, "nlevels") > 1L) {
    refuse("sem_exact_fit_bootstrap() does not take a two-level model.",
           "two_level")
  }
  if (isTRUE(options$.clustered)) {
    cluster <- lavaan::lavInspect(fit, "cluster")
    refuse(c("sem_exact_fit_bootstrap() does not take a fit with a cluster ",
             "variable (", paste(cluster, collapse = ", "), "): the ",
             "resamples draw rows independently, and the refit cannot carry ",
             "the clustering."), "cluster", cluster = cluster)
  }
  if (isTRUE(options$.sampling.weights)) {
    refuse(c("sem_exact_fit_bootstrap() does not take a fit with sampling ",
             "weights: the resamples would not carry the weights' design."),
           "weights")
  }
  if (options$missing %in% c("two.stage", "robust.two.stage")) {
    refuse(c("sem_exact_fit_bootstrap() does not take a two-stage fit: its ",
             "statistic is not the likelihood ratio statistic the bootstrap ",
             "refits. Refit with missing = \"ml\"."), "two_stage")
  }
  if (isTRUE(options$conditional.x) ||
      (isTRUE(options$fixed.x) && length(lavaan::lavNames(fit, "ov.x")))) {
    # lavaan does not take fixed.x = FALSE beside conditional.x = TRUE.
    refuse(c("sem_exact_fit_bootstrap() does not take a fit that conditions ",
             "on observed exogenous variables (fixed.x = TRUE with such ",
             "variables, or conditional.x = TRUE): the transformation moves ",
             "every variable while the conditional likelihood holds some ",
             "fixed. Refit with ",
             if (isTRUE(options$conditional.x)) "conditional.x = FALSE and ",
             "fixed.x = FALSE."),
           "fixed_x", variables = lavaan::lavNames(fit, "ov.x"))
  }
  fiml <- options$missing %in% c("ml", "ml.x")

  # ---- the data, moments, and transformation of every group ----
  G <- lavaan::lavInspect(fit, "ngroups")
  as_list <- function(x) if (G == 1L) list(x) else x
  X <- as_list(lavaan::lavInspect(fit, "data"))
  Sigma <- as_list(lavaan::lavInspect(fit, "cov.ov"))
  Mu <- if (isTRUE(options$meanstructure)) {
    as_list(lavaan::lavInspect(fit, "mean.ov"))
  } else vector("list", G)
  sampstat <- lavaan::lavInspect(fit, "sampstat")
  S <- if (G == 1L) list(sampstat$cov) else lapply(sampstat, `[[`, "cov")
  group_labels <- if (G == 1L) NULL else lavaan::lavInspect(fit, "group.label")
  ov <- colnames(unclass(Sigma[[1L]]))
  transformed <- vector("list", G)
  record <- vector("list", G)
  for (g in seq_len(G)) {
    Xg <- as.matrix(X[[g]])
    if (is.null(colnames(Xg))) colnames(Xg) <- lavaan::lavNames(fit, "ov")
    Xg <- Xg[, ov, drop = FALSE]
    Sg <- unclass(S[[g]])[ov, ov, drop = FALSE]
    Sigg <- unclass(Sigma[[g]])[ov, ov, drop = FALSE]
    Mug <- if (is.null(Mu[[g]])) NULL else as.numeric(Mu[[g]])[match(ov, names(Mu[[g]]))]
    if (anyNA(Xg) && !fiml) {
      refuse(c("The analysis data of group ", g, " have missing values ",
               "although the fit did not use FIML; refit with missing = ",
               "\"ml\"."), "missing_without_fiml", group = g)
    }
    tr <- .bs_transform_group(Xg, Sigg, Mug, Sg, method = transformation,
                              group = g, group_label = group_labels[g])
    transformed[[g]] <- tr$Y
    record[[g]] <- tr[c("method", "patterns", "iterations", "score",
                        "mean_residual")]
  }
  names(record) <- if (G == 1L) "group_1" else lavaan::lavInspect(fit, "group.label")
  n_g <- vapply(transformed, nrow, integer(1))
  N <- sum(n_g)

  # ---- the observed statistic ----
  test <- lavaan::lavTech(fit, "test")[["standard"]]
  T_obs <- as.numeric(test$stat)
  if (is.null(test) || length(T_obs) != 1L || !is.finite(T_obs)) {
    refuse(c("The fit has no standard test statistic (it was fitted with ",
             "test = \"none\"); refit it with a test."), "no_test")
  }
  df <- as.numeric(test$df)
  p_asym <- if (!is.null(test$pvalue)) as.numeric(test$pvalue) else NA_real_

  # ---- the resamples, drawn before any refit ----
  index_sets <- lapply(seq_len(B), function(b) {
    lapply(n_g, function(n) sample.int(n, replace = TRUE))
  })
  refit_options <- options
  refit_options$se <- "none"
  refit_options$test <- "standard"
  refit_options$baseline <- FALSE
  refit_options$h1 <- FALSE
  refit_options$verbose <- FALSE
  refit_options$warn <- FALSE
  group_var <- if (G == 1L) NULL else lavaan::lavInspect(fit, "group")
  # lavaan's slot interface takes the parameter table as a plain list.
  pt <- as.list(lavaan::parTable(fit))
  model <- fit@Model
  refit <- function(rows) {
    .bs_refit_one(rows, transformed, group_var, group_labels, pt, model,
                  refit_options)
  }
  draws <- .bs_run_replications(index_sets, refit, ncpus)
  status <- factor(.BS_STATUS[draws[, "status"]], levels = .BS_STATUS)
  stat <- draws[, "statistic"]

  # ---- the tally ----
  used <- status == "proper" | (improper == "keep" & status == "improper")
  used <- used & is.finite(stat)
  B_used <- sum(used)
  B_improper <- sum(status == "improper")
  B_failed <- sum(status %in% c("not_converged", "error"))
  set_aside <- B - B_used
  parts <- c(
    if (B_failed) paste0(B_failed, " did not converge or failed"),
    if (improper == "drop" && B_improper) paste0(B_improper, " gave an improper solution"))
  # The verdict of the rule (R/dmar_bootstrap.R): the p-value is withheld
  # when more than half of the samples were set aside (the count of 100
  # is held to the request, checked above). One warning, of class dmar_bootstrap_replications_dropped
  # when samples were set aside and the p-value is reported, and of class
  # dmar_bootstrap_insufficient when it is withheld.
  verdict <- .dmar_boot_settle(
    B, B_used, kind = "p", conf_level = conf_level,
    fun = "sem_exact_fit_bootstrap", quantity = "the bootstrap statistic",
    dropped_because = "returned no usable statistic", signal = FALSE)
  enough <- verdict$formed
  .dmar_boot_signal_text(verdict, if (enough && set_aside > 0L) {
    paste0(set_aside, " of ", B, " bootstrap samples were set aside (",
           paste(parts, collapse = "; "), "); the p-value uses the ",
           B_used, " that remain.",
           .dmar_boot_request_note(B, B_used, verdict$B_required,
                                   "bootstrap samples"))
  } else if (!enough) {
    paste0(
      if (set_aside > 0L) paste0(
        set_aside, " of ", B, " bootstrap samples were set aside (",
        paste(parts, collapse = "; "), "). "),
      if (identical(verdict$status, "insufficient_share")) paste0(
        "Fewer than half of the attempted bootstrap samples (", B_used,
        " of ", B, ") returned a usable statistic, so the bootstrap ",
        "p-value is NA: a p-value from the samples that happened to ",
        "converge would describe a selected distribution.")
      else if (!is.na(verdict$B_required) && B < verdict$B_required) paste0(
        "B = ", B, " bootstrap samples were requested, fewer than the ",
        verdict$B_required, " a Monte Carlo p-value is run from, so the ",
        "bootstrap p-value is NA.")
      else paste0(
        "Only ", B_used, " of ", B, " bootstrap samples returned a usable ",
        "statistic, and a p-value needs at least two, so the bootstrap ",
        "p-value is NA."))
  })
  if (enough) {
    exceed <- sum(stat[used] >= T_obs)
    p_boot <- exceed / B_used
    mc <- ci_proportion(successes = exceed, N = B_used, conf_level = conf_level)
    p_lower <- mc$value[mc$term == "ci_lower"]
    p_upper <- mc$value[mc$term == "ci_upper"]
    boot_mean <- mean(stat[used])
  } else {
    p_boot <- p_lower <- p_upper <- boot_mean <- NA_real_
  }
  listwise <- !fiml && sum(lavaan::lavInspect(fit, "norig")) - N > 0.5
  if (listwise) {
    .dmar_warn(
      sprintf(paste0("lavaan removed %d of %d cases by listwise ",
                     "deletion; the bootstrap resamples the %d complete ",
                     "cases the model was fitted to."),
              round(sum(lavaan::lavInspect(fit, "norig")) - N),
              round(sum(lavaan::lavInspect(fit, "norig"))), N),
      "result_qualified")
  }

  vals <- c(chi_square = T_obs, df = df, p_chi_square = p_asym,
            p_bootstrap = p_boot, p_bootstrap_ci_lower = p_lower,
            p_bootstrap_ci_upper = p_upper, bootstrap_mean = boot_mean,
            B = B, B_used = B_used, B_improper = B_improper, B_failed = B_failed,
            N = N, n_groups = G)
  methods <- unique(vapply(record, `[[`, character(1), "method"))
  method_line <- paste0(
    "Transformation: ",
    paste(vapply(methods, function(m) switch(m,
      bollen_stine = "complete data, Bollen and Stine (1992)",
      common = "incomplete data, one matrix for every pattern (Savalei & Yuan, 2009)",
      pattern = "incomplete data, pattern by pattern (Savalei & Yuan, 2009)"),
      character(1)), collapse = "; "),
    "; the statistic is the standard likelihood ratio statistic",
    if (fiml) " under FIML")
  count_line <- paste0(
    "Bootstrap: B = ", B, " samples, ", B_used, " used",
    if (B_improper) paste0(" (", B_improper, " improper, ",
                           if (improper == "keep") "kept" else "set aside", ")"),
    if (B_failed) paste0(", ", B_failed, " failed"),
    "; p_bootstrap is the share of used statistics at least the observed one, ",
    "with the Wilson interval for its Monte Carlo error at ",
    .format_conf_pct(conf_level), " percent",
    if (identical(verdict$status, "exploratory")) paste0(
      ". Exploratory: below the ", verdict$B_required, " samples a ",
      "reported p-value needs"))
  footer <- c(method_line, count_line,
              sprintf("Sample size: N = %s in %s", format(N),
                      if (G == 1L) "1 group" else paste(G, "groups")))
  out <- .as_dmar_tbl(data.frame(term = names(vals), value = unname(vals),
                                 stringsAsFactors = FALSE),
                      conf_level = conf_level,
                      p_terms = c("p_chi_square", "p_bootstrap"),
                      subclass = "dmar_sem_fit")
  attr(out, "bootstrap") <- data.frame(replication = seq_len(B),
                                       chi_square = unname(stat),
                                       status = status,
                                       stringsAsFactors = FALSE)
  attr(out, "transformation") <- record
  attr(out, "sem_fit") <- list(
    source = "fitted lavaan model", test = "standard likelihood ratio statistic",
    chi_square = T_obs, df = df,
    missing = options$missing, estimator = options$estimator.orig,
    transformation_requested = transformation,
    transformation = methods, improper = improper,
    B = B, B_used = B_used, B_improper = B_improper, B_failed = B_failed,
    B_required = verdict$B_required, bootstrap_status = verdict$status,
    N = N, groups = G, footer = footer)
  # The counts in the names every result of the package that ran
  # replications carries (?dmar_bootstrap); the `bootstrap` attribute here
  # is the data frame of the samples.
  for (nm in c("B", "B_used", "n_dropped", "B_required")) {
    attr(out, nm) <- verdict[[nm]]
  }
  attr(out, "bootstrap_status") <- verdict$status
  out
}
