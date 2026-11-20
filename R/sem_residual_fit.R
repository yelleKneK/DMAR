#' The SRMR With Its Unbiased Estimate, Confidence Interval, and Tests of Fit
#'
#' Summarizes how far a fitted model's implied covariances are from the
#' observed ones with the standardized root mean square residual (SRMR) and
#' gives it the inferential treatment the RMSEA has long had: the
#' asymptotically unbiased estimate of the population SRMR of
#' Maydeu-Olivares (2017), its standard error and confidence interval, a
#' test of exact fit, and an optional test of close fit, together with the
#' largest standardized residual. The SRMR is an effect size of misfit in
#' the metric of a correlation, and Maydeu-Olivares, Shi, and Rosseel
#' (2018) found its interval and test of close fit more accurate than
#' those of the RMSEA in models with many variables.
#' The computations are lavaan's
#' (\code{lavaan::lavResiduals()}), and the function requires \pkg{lavaan}
#' 0.7-2 or newer (in Suggests).
#'
#' @param fit A fitted \pkg{lavaan} model: a maximum likelihood fit
#'   (estimator ML, MLM, MLMV, or MLR) under lavaan's default normal
#'   likelihood, of one group of continuous variables at a single level,
#'   from complete raw data (fitted under listwise deletion, FIML, or
#'   pairwise deletion) or from a covariance matrix. It must have
#'   converged to an admissible solution, as
#'   \code{\link{sem_fit_measures}} judges it. See Details for what is
#'   refused and why.
#' @param conf_level Confidence level of the interval for the population
#'   SRMR. Defaults to 0.90, as lavaan's does and as every RMSEA interval in
#'   DMAR does: a 90 percent interval pairs with a test of close fit at .05
#'   (see Details). \code{conf_level = 0.95} gives a 95 percent interval.
#' @param close_fit Optional bound \eqn{c} for the test of close fit,
#'   \eqn{H_0: \mathrm{SRMR} \le c}. Without it, no bound is imposed
#'   and the \code{close_fit_bound} and \code{p_close_fit} rows are
#'   \code{NA}. There is no universal value: Maydeu-Olivares, Shi, and
#'   Rosseel (2018) treat the bound as the researcher's choice. Shi,
#'   Maydeu-Olivares, and DiStefano (2018) proposed population reference
#'   values that scale with the average \eqn{R^2} (communality) of the
#'   observed variables, in a criterion with two parts that must both
#'   hold: a model fits closely when its SRMR is at most .05 times the
#'   average \eqn{R^2} and its largest absolute standardized residual
#'   covariance is at most .10 (and acceptably at .10 times the average
#'   \eqn{R^2} and .15).
#'   They recommend applying these population values through an interval
#'   or a test of close fit rather than to the sample values, so
#'   \code{close_fit = 0.05} times the average \eqn{R^2} tests the first
#'   part (\code{lavaan::lavInspect(fit, "rsquare")} gives the
#'   \eqn{R^2} values). The second part is not tested here; the row
#'   \code{largest_absolute_residual} is its sample value, which they
#'   caution against comparing with .10 unless the sample is very large.
#'   Ximénez, Maydeu-Olivares, Shi, and Revuelta (2022) found the bound of
#'   the first part, applied to the unbiased SRMR, to work reasonably well
#'   in finite samples, except that with small loadings (an average
#'   \eqn{R^2} of .09) close fit called for samples of 500 or more, and
#'   recommend it.
#'
#' @details
#' \subsection{The SRMR and Its Unbiased Estimate}{
#' With \eqn{p} observed variables, sample covariances \eqn{s_{ij}}, and
#' model implied covariances \eqn{\hat{\sigma}_{ij}}, the standardized
#' residuals are \eqn{e_{ij} = (s_{ij} - \hat{\sigma}_{ij}) / \sqrt{s_{ii}
#' s_{jj}}}, and the SRMR is \eqn{\sqrt{\sum_{i \le j} e_{ij}^2 / t}} over
#' the \eqn{t = p(p + 1)/2} variances and covariances (row \code{srmr};
#' lavaan's \code{srmr_bentler_nomean}). Its population counterpart is the
#' same quantity computed from the population covariances, a measure of the
#' size of misfit that does not depend on the sample size. The sample SRMR
#' overestimates it, since sampling error adds to every residual.
#' Maydeu-Olivares (2017) removes the bias on the squared scale: with
#' \eqn{\Xi} the asymptotic covariance matrix of the standardized residuals
#' and \eqn{e} their vector, the unbiased estimate (row
#' \code{unbiased_srmr}) is \eqn{k^{-1}\sqrt{\max(e'e - \mathrm{tr}\,\Xi,
#' 0)/t}}, with \eqn{k = 1 - [\mathrm{tr}(\Xi^2) + 2e'\Xi e]/[4(e'e)^2]}
#' (Maydeu-Olivares et al., 2018, Eq. 16). It is truncated at zero when the
#' estimated sampling contribution exceeds \eqn{e'e}, and the square root
#' and the truncation leave some bias in small samples.
#' }
#'
#' \subsection{The Interval and the Test of Close Fit}{
#' The interval for the population SRMR is the large sample normal
#' reference interval \eqn{\hat{\rho}_u \pm z_{1 - \alpha/2}
#' \widehat{SE}(\hat{\rho}_u)} for the unbiased estimate \eqn{\hat{\rho}_u}
#' (Maydeu-Olivares et al., 2018, Eqs. 20 and 22), where \eqn{\alpha} is 1
#' minus \code{conf_level}. Its limits are reported as lavaan returns them,
#' without truncation, so the lower limit can fall below zero although the
#' population SRMR cannot, and when the estimate is truncated to zero the
#' interval is centered at zero. In the simulations of Maydeu-Olivares et
#' al. (2018) the interval computed without the normality assumption
#' covered at close to the nominal rate in most conditions, including
#' models of 60 variables, and did less well than the RMSEA interval only
#' for very small misfit in small models, and then mostly in small
#' samples; computed under normality, it covered poorly when the data were
#' not normal. The test of
#' close fit, \eqn{H_0: \mathrm{SRMR} \le c}, refers \eqn{z =
#' (\hat{\rho}_u - c)/\widehat{SE}(\hat{\rho}_u)} to the standard normal
#' (Maydeu-Olivares et al., 2018, Eq. 21; row \code{p_close_fit}); it
#' rejects at level \eqn{\alpha} exactly when the lower limit of the
#' \eqn{100(1 - 2\alpha)} percent interval exceeds \eqn{c}, so the
#' default \code{conf_level = 0.90} pairs with a test at .05, as the RMSEA
#' interval of \code{\link{sem_fit_measures}} does with its tests, and
#' \code{conf_level = 0.95} with a test at .025. The printed table states
#' the level beneath it. Failing to reject close fit does not establish
#' it.
#' }
#'
#' \subsection{The Test of Exact Fit}{
#' The test of exact fit (row \code{p_exact_fit}) asks whether the
#' population SRMR is zero. It uses the sample (biased) SRMR, whose
#' expected value and standard error under exact fit are known
#' (Maydeu-Olivares, 2017): \eqn{z = \max(\mathrm{SRMR} - E_0, 0) / SE_0},
#' referred to the upper tail of the standard normal. The row
#' \code{srmr_se_exact_fit} is that standard error, \eqn{SE_0 =
#' \sqrt{\mathrm{tr}(\Xi^2)/(2t\,\mathrm{tr}\,\Xi)}}, which holds only
#' under exact fit: it is not the standard error of the sample SRMR of a
#' misspecified model and must not be used to form an interval, which is
#' what the unbiased estimate and its standard error are for. Because of
#' the truncation at zero the \emph{p}-value is at most 0.5, and it equals
#' 0.5 whenever the sample SRMR is below its expectation under exact fit,
#' as it is more than half the time when the model is correct.
#' }
#'
#' \subsection{Which Fits, and the Residual Covariance}{
#' The asymptotic covariance \eqn{\Xi} of the residuals, and so the
#' unbiased estimate, its standard error, the interval, and both tests,
#' depend on the fit, as lavaan computes them. For an ML fit with standard
#' errors under normality it is normal theory. When the fit computed the
#' asymptotic covariance matrix of the sample moments from their fourth
#' order moments (MLM, MLMV, or a Satorra-Bentler test), or uses sandwich
#' standard errors (MLR), it is asymptotically distribution free, which is
#' the version to use with nonnormal item data, for which normal theory
#' intervals perform poorly (Maydeu-Olivares et al., 2018). The
#' projection onto the residual space uses the fit's information setting,
#' so observed information (the default under MLR) gives slightly different
#' values from expected information; an MLR fit with \code{information =
#' "expected"} reproduces MLM. The standardizing standard deviations are
#' treated as fixed. The \code{sem_fit} attribute records which residual
#' covariance and information were used.
#'
#' The function refuses what these computations do not cover, each with a
#' message that says why: estimators other than ML (least squares and
#' categorical estimators), \code{likelihood = "wishart"}, multiple groups,
#' multilevel models, ordered variables, sampling weights,
#' \code{conditional.x = TRUE} (whose
#' residuals are conditional on the exogenous covariates rather than the
#' marginal covariances the SRMR summarizes), fits with incomplete data
#' under FIML or pairwise deletion, and two-stage fits. A FIML or pairwise
#' fit of data with no missing values is accepted: its residuals, their
#' standard errors, and the tests are those of a listwise fit with the same
#' mean structure and information setting (FIML's default is the observed
#' information, which the footer names). A fit with auxiliary variables
#' (\pkg{lavaan}'s \code{aux}) is refused: \pkg{lavaan}'s residual summary
#' averages over them, and their residuals of about zero lower the SRMR,
#' while with complete data they change nothing in the model's fit, so the
#' fit without them is the one to summarize. A fit whose cases were removed by
#' listwise deletion is evaluated on the complete cases with a warning,
#' since that assumes the data are missing completely at random. A model
#' fitted from a covariance matrix (\code{sample.cov}) gives the same values
#' as one fitted to the raw data under normal theory, since these need only
#' the moments. For a saturated model (zero degrees of freedom) the
#' residuals are zero up to the convergence tolerance, and the standard
#' errors, limits, and \emph{p}-values are \code{NA}, with a warning that
#' says so.
#'
#' The SRMR here summarizes the covariance residuals only. For a model with
#' a mean structure, row \code{srmr_mean} is the root mean square of the
#' standardized mean residuals, zero when the means are saturated; when the
#' model restricts the means, a warning points to it, since the covariance
#' SRMR does not reflect misfit in the means. Always look at the residuals
#' themselves as well as their average (the attribute
#' \code{standardized_residuals}): a good average can hide a large local
#' misfit.
#' }
#'
#' @return A \code{dmar_tbl} (class \code{dmar_sem_fit}) with columns
#'   \code{term} and \code{value} and, always in this order, the rows
#'   \code{srmr} (the sample SRMR), \code{srmr_se_exact_fit} (its standard
#'   error under exact fit only), \code{p_exact_fit},
#'   \code{unbiased_srmr}, \code{unbiased_srmr_se},
#'   \code{unbiased_srmr_ci_lower}, \code{unbiased_srmr_ci_upper},
#'   \code{close_fit_bound} and \code{p_close_fit} (\code{NA} unless
#'   \code{close_fit} is given), \code{srmr_mean} (\code{NA} without a mean
#'   structure), \code{largest_absolute_residual}, and \code{N}. Attributes:
#'   \code{standardized_residuals}, the matrix of standardized covariance
#'   residuals; \code{largest_residual}, the pair of variables with the
#'   largest absolute residual and its signed value; and \code{sem_fit},
#'   which records the residual covariance, the information setting, the
#'   cases used, and the lines of the printed footer.
#'
#' @references
#' Maydeu-Olivares, A. (2017). Assessing the size of model misfit in
#'   structural equation models. \emph{Psychometrika, 82}(3), 533--558.
#'   \doi{10.1007/s11336-016-9552-7}
#'
#' Maydeu-Olivares, A., Shi, D., & Rosseel, Y. (2018). Assessing fit in
#'   structural equation models: A Monte-Carlo evaluation of RMSEA versus
#'   SRMR confidence intervals and tests of close fit.
#'   \emph{Structural Equation Modeling, 25}(3), 389--402.
#'   \doi{10.1080/10705511.2017.1389611}
#'
#' Shi, D., Maydeu-Olivares, A., & DiStefano, C. (2018). The relationship
#'   between the standardized root mean square residual and model
#'   misspecification in factor analysis models. \emph{Multivariate
#'   Behavioral Research, 53}(5), 676--694.
#'   \doi{10.1080/00273171.2018.1476221}
#'
#' Ximénez, C., Maydeu-Olivares, A., Shi, D., & Revuelta, J. (2022).
#'   Assessing cutoff values of SEM fit indices: Advantages of the unbiased
#'   SRMR index and its cutoff criterion based on communality.
#'   \emph{Structural Equation Modeling, 29}(3), 368--380.
#'   \doi{10.1080/10705511.2021.1992596}
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso \code{\link{sem_fit_measures}} for the chi square, RMSEA, and
#'   incremental indices beside the sample SRMR, and
#'   \code{\link{equivalence_rmsea}} for the test of not close fit of the
#'   RMSEA.
#'
#' @family multivariate and latent variable methods
#'
#' @keywords multivariate
#'
#' @examples
#' # The three-factor model of nine holzinger_swineford tests.
#' fit <- lavaan::cfa(
#'   "visual  =~ t1_visual_perception + t2_cubes + t4_lozenges
#'    textual =~ t6_paragraph_comprehension + t7_sentence + t9_word_meaning
#'    speed   =~ t10_addition + t12_counting_groups_of_dots +
#'               t13_straight_and_curved_capitals",
#'   data = holzinger_swineford, std.lv = TRUE)
#' sem_residual_fit(fit)
#'
#' # The 90 percent interval, the default, pairs with a test of close fit
#' # at .05. The bound is the researcher's choice; here it is .05 times the
#' # average R squared of the nine tests, the population value for close
#' # fit of Shi, Maydeu-Olivares, and DiStefano (2018).
#' r2 <- lavaan::lavInspect(fit, "rsquare")
#' sem_residual_fit(fit, close_fit = 0.05 * mean(r2))
#'
#' # A 95 percent interval instead.
#' sem_residual_fit(fit, conf_level = 0.95)
#'
#' @export
sem_residual_fit <- function(fit, conf_level = 0.90, close_fit = NULL) {
  options <- .sem_fit_check(fit)
  .sem_require_lavaan_version("sem_residual_fit")
  .sem_probability(conf_level, "conf_level")
  if (!is.null(close_fit)) {
    .sem_scalar(close_fit, "close_fit", 0, strict = TRUE)
    close_fit <- .sem_number(close_fit)
  }
  if (isTRUE(options$conditional.x)) {
    .dmar_refuse_argument(c(
      "Residual inference is not available for a conditional.x = TRUE fit: ",
      "lavaan then reports residuals conditional on the exogenous covariates ",
      "(res.cov, res.int, res.slopes), not the marginal covariances the SRMR ",
      "summarizes. Refit with conditional.x = FALSE."), arg = NULL)
  }
  # A FIML or pairwise fit of complete data computes nothing a listwise fit
  # does not: lavResiduals() gives the values of a listwise fit with the same
  # mean structure and information setting (to 1e-14). Such a fit is
  # accepted; with missing values, or under two-stage estimation, it is not.
  incomplete <- !is.null(lavaan::lavInspect(fit, "case.idx")) &&
    anyNA(unlist(lavaan::lavInspect(fit, "data")))
  complete_ok <- options$missing %in% c("listwise", "ml", "ml.x", "pairwise") &&
    !incomplete
  scope <- c(
    ordered = length(lavaan::lavNames(fit, "ov.ord")) > 0L,
    estimator = !identical(options$estimator, "ML"),
    likelihood = !identical(options$likelihood, "normal"),
    groups = lavaan::lavInspect(fit, "ngroups") != 1L,
    levels = lavaan::lavInspect(fit, "nlevels") != 1L,
    weights = isTRUE(options$.sampling.weights),
    missing = !complete_ok)
  if (any(scope)) {
    why <- c(
      estimator = paste0("a maximum likelihood fit (estimator ML, MLM, MLMV, or ",
                         "MLR); this fit used ", options$estimator),
      likelihood = "lavaan's default normal likelihood, not likelihood = \"wishart\"",
      groups = "a single group",
      levels = "a single level",
      ordered = "continuous variables, not ordered ones",
      weights = paste0("an unweighted fit: the residual covariance theory ",
                       "assumes independent, equally weighted cases"),
      missing = paste0("complete data: raw data with no missing values (under ",
                       "listwise deletion, FIML, or pairwise deletion) or a ",
                       "covariance matrix; this fit used missing = \"",
                       options$missing, "\"",
                       if (incomplete) " with missing values in its analysis data"
                       else ", which this function does not cover"))
    .dmar_refuse_argument(c(
      "sem_residual_fit() requires ", why[names(scope)[scope][1L]], "."),
      arg = NULL)
  }
  # lavaan's residual inference averages over every observed variable, the
  # auxiliary variables among them, whose residuals of about zero lower the
  # SRMR. With the complete data this function takes, an auxiliary variable
  # changes nothing in the model's fit, so the fit without it is the one to
  # summarize.
  aux <- .sem_aux_names(fit)
  if (length(aux)) {
    .dmar_refuse_argument(c(
      "sem_residual_fit() does not take a fit with auxiliary variables (",
      paste(aux, collapse = ", "), "): lavaan's residual summary averages ",
      "over them, and their residuals of about zero lower the SRMR. With ",
      "complete data an auxiliary variable changes nothing in the model's ",
      "fit; refit without aux."), arg = NULL)
  }
  residuals <- lavaan::lavResiduals(fit, type = "cor.bentler", se = TRUE,
                                    zstat = TRUE, summary = TRUE,
                                    h1_acov = "unstructured",
                                    usrmr_ci_level = conf_level,
                                    usrmr_close_h0 = if (is.null(close_fit)) 0.05 else close_fit)
  s <- residuals$summary
  if (is.null(s) || !"cov" %in% colnames(s)) {
    .dmar_refuse_argument(c(
      "lavaan returned no covariance residual summary for this fit (its ",
      "columns are: ", paste(colnames(s), collapse = ", "), ")."), arg = NULL)
  }
  keys <- c(srmr = "srmr", srmr_se_exact_fit = "srmr.se",
            p_exact_fit = "srmr.exactfit.pvalue",
            unbiased_srmr = "usrmr", unbiased_srmr_se = "usrmr.se",
            unbiased_srmr_ci_lower = "usrmr.ci.lower",
            unbiased_srmr_ci_upper = "usrmr.ci.upper",
            close_fit_bound = "usrmr.closefit.h0.value",
            p_close_fit = "usrmr.closefit.pvalue")
  vals <- stats::setNames(as.numeric(s[keys, "cov"]), names(keys))
  if (is.null(close_fit)) vals[c("close_fit_bound", "p_close_fit")] <- NA_real_

  # The mean residuals, for a model with a mean structure.
  vals[["srmr_mean"]] <- if ("mean" %in% colnames(s)) {
    as.numeric(s["srmr", "mean"])
  } else NA_real_
  if (isTRUE(vals[["srmr_mean"]] > 1e-6)) {
    .dmar_warn(c(
      "The model restricts the means, and srmr summarizes the covariance ",
      "residuals only; srmr_mean (", format(signif(vals[["srmr_mean"]], 3)),
      ") summarizes the standardized mean residuals."), "result_qualified")
  }

  # The largest standardized covariance residual.
  mat <- residuals$cov
  pairs <- which(lower.tri(mat), arr.ind = TRUE)
  if (nrow(pairs)) {
    ij <- pairs[which.max(abs(mat[pairs])), ]
    largest <- data.frame(item_1 = colnames(mat)[ij[2]],
                          item_2 = rownames(mat)[ij[1]],
                          residual = mat[ij[1], ij[2]],
                          stringsAsFactors = FALSE)
    vals[["largest_absolute_residual"]] <- abs(largest$residual)
  } else {
    largest <- data.frame(item_1 = character(), item_2 = character(),
                          residual = numeric(), stringsAsFactors = FALSE)
    vals[["largest_absolute_residual"]] <- NA_real_
  }
  n_used <- lavaan::lavInspect(fit, "ntotal")
  n_orig <- sum(lavaan::lavInspect(fit, "norig"))
  vals[["N"]] <- n_used

  # Warnings: a saturated model, otherwise any quantity lavaan left NA, and
  # cases removed by listwise deletion.
  inferential <- setdiff(names(keys), c("srmr", "unbiased_srmr",
                                        if (is.null(close_fit))
                                          c("close_fit_bound", "p_close_fit")))
  df <- .sem_number(lavaan::fitMeasures(fit, "df"))
  if (isTRUE(df == 0)) {
    .dmar_warn(c(
      "The model is saturated (zero degrees of freedom): it reproduces ",
      "the sample covariances, so the SRMR is zero up to the convergence ",
      "tolerance, and no standard error, test, or interval is defined; ",
      "they are NA."), "result_qualified")
  } else if (anyNA(vals[inferential])) {
    .dmar_warn(c(
      "lavaan could not compute ",
      paste(inferential[is.na(vals[inferential])], collapse = ", "),
      " for this fit; NA is returned."), "result_qualified")
  }
  if (n_orig - n_used > 0.5) {
    .dmar_warn(
      sprintf(paste0("lavaan removed %d of %d cases by listwise deletion; ",
                     "the residuals use the %d complete cases, which ",
                     "assumes the data are missing completely at random."),
              round(n_orig - n_used), round(n_orig), round(n_used)),
      "result_qualified")
  }

  # Which residual covariance lavaan used. lavResiduals() takes the
  # asymptotic covariance of the sample moments from the fit's stored
  # fourth moment matrix (NACOV) when the fit computed one, and otherwise
  # from the h1 information, which is a sandwich for robust standard
  # errors; the slot is read here as lavResiduals() reads it.
  nacov <- tryCatch(!is.null(fit@SampleStats@NACOV[[1L]]), error = function(e) NA)
  residual_acov <- if (isTRUE(nacov)) {
    "asymptotically distribution free (fourth order moments of the data)"
  } else if (options$se %in% c("robust.huber.white", "robust.cluster",
                               "robust.cluster.sem")) {
    "asymptotically distribution free (sandwich, as for MLR)"
  } else if (isFALSE(nacov)) {
    "normal theory"
  } else "not determined"
  information <- options$information[1L]
  out <- .as_dmar_tbl(data.frame(term = names(vals), value = unname(vals),
                                 stringsAsFactors = FALSE),
                      conf_level = conf_level,
                      p_terms = c("p_exact_fit", "p_close_fit"),
                      subclass = "dmar_sem_fit")
  attr(out, "standardized_residuals") <- mat
  attr(out, "largest_residual") <- largest
  attr(out, "sem_fit") <- list(
    source = "fitted lavaan model", residual_type = "cor.bentler, covariances",
    residual_acov = residual_acov, information = information,
    h1_acov = "unstructured", estimator = options$estimator,
    N = n_used, N_original = n_orig,
    close_fit_null = if (is.null(close_fit)) "none" else "population SRMR <= close_fit",
    footer = c(sprintf("Residual covariance: %s; %s information", residual_acov,
                       information),
               if (nrow(largest)) {
                 sprintf("Largest standardized residual: %s with %s, %s",
                         largest$item_1, largest$item_2,
                         formatC(largest$residual, format = "f", digits = 3))
               }))
  out
}
