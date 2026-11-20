#' A Priori Monte Carlo Simulation for Sample Size Planning for RMSEA in SEM
#'
#' @description
#' Conduct a priori Monte Carlo simulation to empirically study the effects of
#' (mis)specifications of input information on the calculated sample size. The
#' sample size is planned so that the expected width of a confidence interval
#' for the population RMSEA is no larger than desired. Random data are generated
#' from the true covariance matrix but fit to the proposed model, whereas the
#' sample size is calculated based on the input covariance matrix and proposed
#' model.
#'
#' @param width desired confidence interval width for the population RMSEA.
#' @param model the model the researcher proposes, which may or may not be the
#'   true model, written in \pkg{lavaan} model syntax (see
#'   \code{\link[lavaan]{model.syntax}}). The observed variable names in the
#'   model must match the row and column names of \code{true_cov_matrix}.
#' @param true_cov_matrix the true population covariance matrix, which is
#'   used to generate random data for the simulation study. The row and column
#'   names of
#'   \code{true_cov_matrix} must match the observed variables in \code{model}.
#' @param specified_N the total sample size \emph{N}. If \code{specified_N}
#'   is given, random samples of that size are generated. Otherwise the sample size is calculated with the sample size
#'   planning method so that the expected width of a confidence interval for the
#'   population RMSEA is no larger than \code{width}, for the interval formed
#'   as lavaan forms it (see Details).
#' @param conf_level confidence level of the planned and the simulated
#'   intervals (i.e., 1 - the Type I error rate). The default, 0.90, is the
#'   conventional level of the RMSEA interval and the default of
#'   \code{\link{ss_aipe_rmsea}} and \code{\link{ci_rmsea}}.
#' @param G number of replications in the Monte Carlo simulation.
#' @param filename an optional path for a comma separated file recording
#'   every converged replication (its index, the RMSEA estimate, the two
#'   confidence limits, and the interval width): nothing is written when
#'   \code{filename} is \code{NULL} (the default), a new file with a header
#'   row is created otherwise, an existing file at that path is appended to,
#'   and a throwaway run should point it at
#'   \code{tempfile(fileext = ".csv")}.
#' @param \dots additional arguments passed to \code{\link[lavaan]{sem}} when
#'   fitting the model (for example \code{estimator} or \code{missing}).
#'   Every name given here must be an argument or an option of \pkg{lavaan}
#'   (see \code{\link[lavaan]{lavOptions}}); any other name, and any
#'   argument without a name, is refused before a model is fitted, with a
#'   condition of class \code{dmar_argument_invalid} that gives the name
#'   as it was typed
#' @param seed \code{NULL} (the default) to draw from the current random
#'   number stream, or a single number that makes the call reproducible; a
#'   supplied seed is set for this call only, and the caller's random
#'   number state is restored when the function returns.
#'
#' @details
#' This function implements the sample size planning method proposed in Kelley
#' and Lai (2011). It uses \code{\link[lavaan]{sem}} to fit the proposed model
#' to the population covariance matrix, which recovers the population RMSEA (the
#' model misspecification) and the model degrees of freedom, and to fit the
#' model to each simulated sample, and it forms the confidence interval for
#' the population RMSEA in each replication from that sample's chi square
#' statistic, inverted with \code{\link{ci_nc_chisq}} as
#' \code{\link{sem_fit_measures}} inverts it. The model is specified in
#' \pkg{lavaan} syntax, so \pkg{lavaan} must be installed.
#'
#' The population RMSEA is read from \code{lavaan::fitMeasures()}, which is
#' computed reliably for the large-sample population fit.
#'
#' Each replication's interval is formed from the statistic \eqn{T} of that
#' replication's fit and the sample size multiplier \eqn{a} lavaan formed
#' \eqn{T} with (\code{sample_size_convention} in \code{\link{ci_rmsea}}):
#' \eqn{a = N} for maximum likelihood under lavaan's default normal
#' likelihood, and \eqn{a = N - 1} under \code{likelihood = "wishart"} and
#' for the least squares estimators, either of which \code{\dots} can
#' request. That is the interval lavaan itself reports for the sample,
#' whatever the size of \eqn{T}. When \eqn{T} is at or below its degrees of
#' freedom, lavaan reports an RMSEA of 0, from which \code{\link{ci_rmsea}}
#' could rebuild only \eqn{T = df}; the interval formed from \eqn{T} itself
#' is narrower the smaller \eqn{T} is, and its upper limit is 0 when
#' \eqn{T} lies at or below the \code{(1 - conf_level) / 2} quantile of the
#' central chi square, the convention lavaan and
#' \code{\link{sem_fit_measures}} share. The replications of a model that
#' fits well often have such statistics, so the distinction matters most
#' when the proposed model is close to the true one. The RMSEA summaries are
#' those of the RMSEA lavaan reports.
#'
#' The width of the interval depends on the sample size only through
#' \eqn{a}, and \code{\link{ss_aipe_rmsea}} plans with \eqn{a = N - 1}, the
#' convention of Kelley and Lai (2011). The simulated samples have the
#' planner's \code{necessary_N} cases, the size a study would collect, so
#' under lavaan's default \eqn{a = N} each realized interval is the one
#' lavaan reports for a sample of that size, whose multiplier is one larger
#' than the plan's, and the realized widths sit below the planned width by
#' an amount of order \eqn{1/N}.
#'
#' The RMSEA is conventionally reported with a 90 percent interval, because
#' a \eqn{100(1 - 2\alpha)} percent interval decides the one-sided tests of
#' close and not close fit at level \eqn{\alpha} (MacCallum, Browne, &
#' Sugawara, 1996, p. 136), and that is the default here, as in
#' \code{\link{ss_aipe_rmsea}}.
#'
#' @return
#' A \code{data.frame} with columns \code{term} and \code{value}
#' summarizing the a priori Monte Carlo study. The \code{term} entries
#' are: \code{"mean_rmsea"}, \code{"median_rmsea"}, \code{"sd_rmsea"}
#' (summaries of the realized RMSEA estimates across the converged
#' replications); \code{"mean_ci_width"}, \code{"median_ci_width"},
#' \code{"sd_ci_width"} (summaries of the realized interval widths);
#' \code{"pct_ci_less_width"} (proportion of intervals narrower than the
#' target width); \code{"pct_ci_miss_low"} and \code{"pct_ci_miss_high"}
#' (tail-specific empirical non-coverage of the population RMSEA);
#' \code{"total_type_I_error"} (overall empirical non-coverage, the sum
#' of the two tails); and the echoes \code{"total_N"} (the \emph{N} evaluated),
#' \code{"df"} (model degrees of freedom), \code{"true_rmsea"} (the
#' population RMSEA recovered from fitting \code{model} to
#' \code{true_cov_matrix}), \code{"width"}, and \code{"conf_level"}. The
#' proportion rows are on the 0 to 1 scale, not percentages. The number
#' of converged replications is \code{attr(x, "simulation")$returned}
#' (see the section on the record of the simulation).
#'
#' @references
#' Cudeck, R., & Browne, M. W. (1992). Constructing a covariance matrix
#'   that yields a specified minimizer and a specified minimum
#'   discrepancy function value.
#'   \emph{Psychometrika, 57}, 357--369. \doi{10.1007/BF02295424}
#'
#' Kelley, K., & Lai, K. (2011). Accuracy in parameter estimation for the
#'   root mean square error of approximation: Sample size planning for
#'   narrow confidence intervals.
#'   \emph{Multivariate Behavioral Research, 46}, 1--32.
#'   \doi{10.1080/00273171.2011.543027}
#'
#' MacCallum, R. C., Browne, M. W., & Sugawara, H. M. (1996). Power
#'   analysis and determination of sample size for covariance structure
#'   modeling. \emph{Psychological Methods, 1}(2), 130--149.
#'   \doi{10.1037/1082-989X.1.2.130}
#'
#' Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). \emph{Designing
#'   experiments and analyzing data: A model comparison perspective}
#'   (4th ed.). Routledge.
#'
#' Rosseel, Y. (2012). lavaan: An R package for structural equation
#'   modeling. \emph{Journal of Statistical Software, 48}(2), 1--36.
#'   \doi{10.18637/jss.v048.i02}
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @note
#' Replications in which \pkg{lavaan} fails to converge, or for which the RMSEA
#' is undefined, are skipped; the number of converged replications is reported
#' as \code{attr(x, "simulation")$returned}. A fit whose likelihood has no maximum, an error
#' variance running off below zero while the likelihood is flat or still
#' rising (see \code{\link{fa_bifactor}}), is used as it is, at
#' \code{true_cov_matrix}
#' and in a replication: its estimates are wherever \pkg{lavaan}'s optimizer
#' stopped on that ridge, but the RMSEA is formed from the chi square alone,
#' which is not, since the likelihood levels off along the ridge and the
#' optimizer stops only where its steps no longer gain. Increase \code{G} if
#' many replications fail to converge.
#'
#' @seealso \code{\link[lavaan]{sem}}, \code{\link{ss_aipe_rmsea}},
#'   \code{\link{ci_rmsea}}, \code{\link{sem_fit_measures}}
#'
#' @examples
#' # True data generating model: two correlated factors, each measured by
#' # three standardized indicators. The factor correlation is 0.5 and every
#' # loading is 0.7. The implied population covariance matrix is assembled
#' # from the loading matrix, the factor correlation matrix, and the
#' # residual variances.
#' Lambda <- matrix(0, 6, 2)
#' Lambda[1:3, 1] <- 0.7
#' Lambda[4:6, 2] <- 0.7
#' Phi   <- matrix(c(1, 0.5, 0.5, 1), 2, 2)
#' true_cov_matrix <- Lambda %*% Phi %*% t(Lambda) + diag(1 - 0.7^2, 6)
#' dimnames(true_cov_matrix) <- list(paste0("x", 1:6), paste0("x", 1:6))
#'
#' # Proposed (misspecified) model: a single common factor.
#' proposed <- "g =~ x1 + x2 + x3 + x4 + x5 + x6"
#'
#' # The proposed model is fit once at a very large N to recover the
#' # population RMSEA, the sample size is planned so that the expected width
#' # of the 90 percent interval is 0.05, and a fresh sample of that size is
#' # drawn and fit on every replication. Notice that true_rmsea is about
#' # 0.20, since a single factor is a poor description of two-factor data,
#' # that the realized widths sit close to the target, and that
#' # pct_ci_less_width is well away from 0 and 1: planning for the expected
#' # width aims the average width at the target, not the share of intervals
#' # narrower than it, which is what assurance controls. G = 10 keeps
#' # the example quick and estimates that share only roughly; a reported
#' # sensitivity study deserves the default G = 200 or more.
#' set.seed(113)
#' ss_aipe_rmsea_sensitivity(width = 0.05, model = proposed,
#'                           true_cov_matrix = true_cov_matrix, G = 10)
#'
#' @inheritSection ss_aipe_cv_sensitivity The record of the simulation
#' @export
ss_aipe_rmsea_sensitivity <- function(width, model, true_cov_matrix,
                                      specified_N = NULL, conf_level = 0.90, G = 200,
                                      filename = NULL, ...,
                                      seed = NULL) {
  .dmar_check_lavaan_dots(..., fn = "ss_aipe_rmsea_sensitivity")
  if (!requireNamespace("MASS", quietly = TRUE)) .dmar_refuse(
    "The package 'MASS' is needed; please install the package and try again.",
    "package_required", package = "MASS")
  if (!requireNamespace("lavaan", quietly = TRUE)) .dmar_refuse(
    "The package 'lavaan' is needed; please install the package and try again.",
    "package_required", package = "lavaan")

  if (is.null(rownames(true_cov_matrix)) || is.null(colnames(true_cov_matrix))) {
    .dmar_refuse_argument(
      "'true_cov_matrix' must have row and column names that match the observed variables in 'model'.",
      arg = c("true_cov_matrix", "model"))
  }
  .check_filename(filename)
  G <- .dmar_check_count(G, "G")
  .dmar_local_seed(seed)

  # A borderline sample makes lavaan warn through the per-replication fit:
  # that the optimizer has not found a solution, that an estimated observed
  # variable variance is negative, that the gradient at the reported
  # solution is not near zero, or that a poor marker item was swapped for
  # another. Nonconvergence is handled by the convergence check in the loop,
  # which skips the replication; the others describe one sample's fit, and
  # the summary statistics absorb them as sampling variability. Each is
  # noise at the level of the Monte Carlo study, so every warning lavaan
  # raises inside this one fit (all carry lavaan's prefix) is muffled here.
  # Warnings from the population fit, the sample size planning, and the
  # inversion of each statistic reach the caller.
  fit_one <- function(Data) {
    withCallingHandlers(
      lavaan::sem(model, data = as.data.frame(Data), ...),
      warning = function(w) {
        if (grepl("^lavaan", conditionMessage(w))) {
          invokeRestart("muffleWarning")
        }
      })
  }

  # Population fit: fit the proposed (possibly misspecified) model to the true
  # covariance matrix at a very large N to recover the population RMSEA and the
  # model degrees of freedom. lavaan::fitMeasures() returns the population RMSEA
  # reliably for this large-sample fit.
  M_fit <- lavaan::sem(model, sample.cov = true_cov_matrix, sample.nobs = 1e6, ...)
  # Whether this fit, or a replication's below, has a maximum
  # (.lavaan_no_maximum(), R/no_maximum_internals.R) is not asked: the
  # RMSEA is formed from the chi square alone, and on a ridge, where an
  # error variance runs off below zero with the likelihood flat or still
  # rising, the likelihood levels off and lavaan's optimizer stops only
  # where its steps no longer gain, so the RMSEA, unlike the estimates,
  # does not depend on where on the ridge it stops.
  rmsea <- unname(lavaan::fitMeasures(M_fit, "rmsea"))
  df    <- unname(lavaan::fitMeasures(M_fit, "df"))
  # The multiplier with which lavaan forms the statistic, and so the RMSEA,
  # of every replication's fit (the arguments in ... are the same for each):
  # N for ML under its default normal likelihood, N - 1 under likelihood =
  # "wishart" and for the least squares estimators (.sem_lavaan_convention(),
  # R/sem_fit_internals.R). Each interval is formed with it. Through DMAR
  # 1.0.x every interval was formed with N - 1 from the RMSEA, so a default
  # fit's RMSEA, computed with N, gave an interval that matched neither
  # convention, and a statistic at or below df was given the interval of a
  # statistic equal to df.
  convention <- .sem_lavaan_convention(M_fit)$convention

  # The sample is the planner's necessary_N, the size the user would
  # collect. The planner works with a = N - 1 (Kelley & Lai, 2011); under
  # lavaan's default the replications' statistics are formed with a = N, so
  # the realized intervals are those lavaan reports for a sample of that
  # size, narrower than the planned width by an amount of order 1/N.
  N <- if (is.null(specified_N)) {
    ss_aipe_rmsea(true_rmsea = rmsea, df = df, width = width,
                  conf_level = conf_level)[1, 2]
  } else {
    specified_N
  }
  # The multiplier of every replication's statistic, which has N cases.
  a <- .sem_multiplier(N, 1L, convention)
  p <- dim(true_cov_matrix)[1]
  rmsea_hat <- rep(NA_real_, G)
  CI_upper  <- rep(NA_real_, G)
  CI_lower  <- rep(NA_real_, G)

  if (!is.null(filename)) {
    if (file.exists(filename) && file.size(filename) > 0) {
      message("The file '", filename, "' already exists; the simulation results are appended to it.")
    } else {
      res_col_names <- data.frame("iteration", "rmsea_hat", "ci_lower", "ci_upper", "width")
      utils::write.table(res_col_names, filename, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
    }
  }

  for (g in seq_len(G)) {
    Data <- MASS::mvrnorm(n = N, mu = rep(0, p), Sigma = true_cov_matrix)
    colnames(Data) <- rownames(true_cov_matrix)

    m_fit <- try(fit_one(Data), silent = TRUE)
    if (inherits(m_fit, "try-error") || !isTRUE(lavaan::lavInspect(m_fit, "converged"))) {
      next
    }
    fm <- unname(lavaan::fitMeasures(m_fit, c("chisq", "rmsea")))
    if (length(fm) != 2L || anyNA(fm)) {
      next
    }
    rmsea_hat[g] <- fm[2L]
    # The interval is formed from the replication's own statistic, as
    # lavaan and sem_fit_measures() form it, not from its RMSEA: lavaan
    # truncates the RMSEA at 0 when the statistic is at or below df, and
    # from an RMSEA of 0 the statistic could be rebuilt only as df.
    ncp <- .sem_ncp_bounds(fm[1L], df, conf_level)
    CI_lower[g] <- .sem_rmsea_from_ncp(ncp[["lower"]], df, a, 1L)
    CI_upper[g] <- .sem_rmsea_from_ncp(ncp[["upper"]], df, a, 1L)
    if (!is.null(filename)) {
      sim_result <- cbind(g, rmsea_hat[g], CI_lower[g], CI_upper[g], CI_upper[g] - CI_lower[g])
      utils::write.table(sim_result, filename, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
    }
  }

  w <- CI_upper - CI_lower
  suc_rep <- sum(!is.na(w))

  term <- c(
    "mean_rmsea", "median_rmsea", "sd_rmsea",
    "mean_ci_width", "median_ci_width", "sd_ci_width",
    "pct_ci_less_width", "pct_ci_miss_low", "pct_ci_miss_high", "total_type_I_error",
    "total_N", "df", "true_rmsea", "width", "conf_level"
  )
  value <- c(
    mean(rmsea_hat, na.rm = TRUE), median(rmsea_hat, na.rm = TRUE), sd(rmsea_hat, na.rm = TRUE),
    mean(w, na.rm = TRUE), median(w, na.rm = TRUE), sd(w, na.rm = TRUE),
    sum(w < width, na.rm = TRUE) / suc_rep,
    sum(CI_lower > rmsea, na.rm = TRUE) / suc_rep,
    sum(CI_upper < rmsea, na.rm = TRUE) / suc_rep,
    sum(CI_lower > rmsea, na.rm = TRUE) / suc_rep + sum(CI_upper < rmsea, na.rm = TRUE) / suc_rep,
    N, df, rmsea, width, conf_level
  )

  out <- .as_dmar_tbl(data.frame(term, value), conf_level = conf_level)
  .ss_sens_record(out, G, !is.na(w), conf_level)
}
