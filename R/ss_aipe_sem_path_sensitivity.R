#' A Priori Monte Carlo Simulation for Sample Size Planning for SEM Targeted Effects
#'
#' @description
#' Conduct a priori Monte Carlo simulation to empirically study the effects of (mis)specifications of input
#' information on the calculated sample size. Random data are generated from the true covariance matrix but
#' fit to the proposed model, whereas sample size is calculated based on the input covariance matrix and
#' proposed model.
#'
#' @param model A single character string giving the free analysis model in
#'   lavaan model syntax (see \code{\link[lavaan]{model.syntax}}), the model
#'   that would be fit to the data. The target path must carry a parameter
#'   label so it can be referred to by name, for example \code{"f2 ~ b*f1"}
#'   labels the structural path \code{b}. The model may or may not be the true
#'   data generating model
#' @param supposed_cov_matrix the population covariance matrix supposed when
#'   the sample size is planned (passed to \code{\link{ss_aipe_sem_path}} as
#'   its \code{true_cov_matrix}). It is a planning assumption, not a matrix
#'   estimated from data, and it may or may not be the true covariance
#'   matrix. The row names and column names of \code{supposed_cov_matrix}
#'   should be the same as the observed variables in \code{model}
#' @param true_cov_matrix the true population covariance matrix, which will
#'   be used to generate random data for the simulation study. The row names
#'   and column names of \code{true_cov_matrix} should be the same as the
#'   observed variables in \code{model}. Defaults to
#'   \code{supposed_cov_matrix}
#' @param which_path the parameter label of the targeted path, given as a
#'   character string, for example \code{"b"} for the path labeled
#'   \code{f2 ~ b*f1} in \code{model}
#' @param width desired confidence interval width for the model parameter of interest
#' @param specified_N the total sample size \emph{N} of each set of random
#'   data. If it is \code{NULL}, it is the sample size
#'   \code{\link{ss_aipe_sem_path}} plans from \code{supposed_cov_matrix}
#' @param conf_level confidence level (i.e., 1- Type I error rate)
#' @param assurance the assurance that the confidence interval obtained in a particular study will be no wider than desired (must be \code{NULL} or a value between 0.50 and 1)
#' @param G number of replications in the Monte Carlo simulation
#' @param filename an optional path for a comma separated file recording
#'   every replication (the estimate of the targeted path, its standard
#'   error, the two confidence limits, and the interval width): nothing is
#'   written when \code{filename} is \code{NULL} (the default), a new file
#'   with a header row is created otherwise, an existing file at that path
#'   is appended to, and a throwaway run should point it at
#'   \code{tempfile(fileext = ".csv")}
#' @param \dots further arguments passed to \code{\link[lavaan]{sem}} for
#'   every fit: the planning fit to \code{supposed_cov_matrix} (through
#'   \code{\link{ss_aipe_sem_path}}), the fit to \code{true_cov_matrix}, and the
#'   fit to each simulated sample (for example \code{std.lv = TRUE}).
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
#' This function implements the sample size planning methods proposed in Lai
#' and Kelley (2011). It calls \code{\link{ss_aipe_sem_path}} to plan the
#' sample size and to identify the targeted path, then fits the analysis model
#' to data simulated from \code{true_cov_matrix} with
#' \code{\link[lavaan]{sem}}. The
#' analysis model is written in lavaan model syntax with the targeted path
#' given a parameter label; see \code{\link[lavaan]{model.syntax}} for the
#' syntax and \code{\link[lavaan]{sem}} for the fitting machinery. The
#' population covariance matrices are most naturally produced by
#' \code{\link{cov_sem}} from a fully fixed population model. This function
#' requires \pkg{lavaan} and \pkg{MASS} to be installed.
#'
#' @return
#' A \code{data.frame} with columns \code{term} and \code{value}
#' summarizing the a priori Monte Carlo study. The \code{term} entries
#' are: \code{"mean_path"}, \code{"median_path"}, \code{"sd_path"}
#' (summaries of the realized estimates of the targeted path across the
#' converged replications); \code{"mean_ci_width"},
#' \code{"median_ci_width"}, \code{"sd_ci_width"} (summaries of the
#' realized interval widths); \code{"pct_ci_less_width"} (proportion of
#' realized widths at or below \code{width});
#' \code{"pct_ci_miss_low"} and \code{"pct_ci_miss_high"}
#' (tail-specific empirical non-coverage of the population path);
#' \code{"total_type_I_error"} (overall empirical non-coverage, the sum
#' of the two tails); and the echoes \code{"total_N"} (the \emph{N} evaluated),
#' \code{"true_path"} (the population value of the targeted path under
#' \code{true_cov_matrix}), \code{"width"}, \code{"conf_level"}, and
#' \code{"assurance"} (present only when an assurance was supplied). The
#' proportion rows are on the 0 to 1 scale, not percentages. The number
#' of replications that gave an interval is
#' \code{attr(x, "simulation")$returned}, and the number of data sets
#' drawn again is \code{attr(x, "simulation")$n_redrawn} (see the
#' sections below).
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
#' Rosseel, Y. (2012). lavaan: An R package for structural equation modeling.
#'   \emph{Journal of Statistical Software, 48}(2), 1--36.
#'   \doi{10.18637/jss.v048.i02}
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @note
#' Occasionally a replication fails to converge when the analysis model is fit
#' to a simulated data set. Such replications are not counted toward the
#' \code{G} converged replications; the simulation draws fresh data and
#' continues. A replication whose likelihood has no maximum, an error
#' variance running off below zero while the likelihood is flat or still
#' rising (see \code{\link{fa_bifactor}}), is treated the same way, and an
#' analysis model with no maximum at \code{supposed_cov_matrix} or
#' \code{true_cov_matrix} stops the function, before any replication is drawn,
#' with an error of
#' class \code{dmar_no_maximum} that prints the \pkg{lavaan} code of that
#' fit. A safety cap stops the loop after \code{20 * G} attempts, and a
#' single warning is issued if fewer than \code{G} replications converged.
#'
#' @seealso
#' \code{\link[lavaan]{sem}}, \code{\link{cov_sem}},
#' \code{\link{ss_aipe_sem_path}}
#'
#' @examples
#' # The planning values a researcher would bring to ss_aipe_sem_path(): two
#' # factors, each measured by three indicators, joined by a structural path
#' # of 0.5, with every residual variance at 0.5.
#' planning_model <- "
#'   f1 =~ 1*y1 + 0.8*y2 + 0.8*y3
#'   f2 =~ 1*y4 + 0.8*y5 + 0.8*y6
#'   f2 ~ 0.5*f1
#'   f1 ~~ 1*f1
#'   f2 ~~ 0.75*f2
#'   y1 ~~ 0.5*y1; y2 ~~ 0.5*y2; y3 ~~ 0.5*y3
#'   y4 ~~ 0.5*y4; y5 ~~ 0.5*y5; y6 ~~ 0.5*y6
#' "
#'
#' # The population the study will actually sample from: the same structural
#' # path, but noisier indicators than the planning values assumed, with
#' # every residual variance at 0.8.
#' true_model <- "
#'   f1 =~ 1*y1 + 0.8*y2 + 0.8*y3
#'   f2 =~ 1*y4 + 0.8*y5 + 0.8*y6
#'   f2 ~ 0.5*f1
#'   f1 ~~ 1*f1
#'   f2 ~~ 0.75*f2
#'   y1 ~~ 0.8*y1; y2 ~~ 0.8*y2; y3 ~~ 0.8*y3
#'   y4 ~~ 0.8*y4; y5 ~~ 0.8*y5; y6 ~~ 0.8*y6
#' "
#'
#' # The free analysis model, with the targeted path labeled "b".
#' analysis_model <- "
#'   f1 =~ y1 + y2 + y3
#'   f2 =~ y4 + y5 + y6
#'   f2 ~ b*f1
#' "
#'
#' supposed_cov_matrix <- cov_sem(planning_model)$sigma_theta
#' true_cov_matrix <- cov_sem(true_model)$sigma_theta
#'
#' # The sample size planned from the optimistic measurement quality is
#' # evaluated against the population that actually holds. Notice that the
#' # mean realized width (mean_ci_width) exceeds the desired 0.30, and that
#' # pct_ci_less_width, the proportion of intervals no wider than desired, falls
#' # well below the one half that the planned sample size would deliver if
#' # the planning values were right. G = 10 keeps the example quick; a
#' # reported sensitivity study deserves the default G = 100 or more.
#' set.seed(113)
#' ss_aipe_sem_path_sensitivity(model = analysis_model,
#'                              supposed_cov_matrix = supposed_cov_matrix,
#'                              true_cov_matrix = true_cov_matrix,
#'                              which_path = "b", width = 0.30, G = 10)
#'
#' @keywords design multivariate
#'
#' @seealso \code{\link{design_consequences}} for what a chosen design delivers:
#'   power, the Type S (sign) and Type M (exaggeration) errors of the
#'   significance filter, and the expected confidence interval width.
#'
#' @section Replications drawn again:
#' A simulated data set whose model fit does not converge, or whose
#' estimate is not a maximum of the likelihood, is replaced by a new one, so such a data set is not among the
#' replications without an interval and is in no summary row.
#' At most \code{20 * G} data sets are drawn in all; if that many do not
#' give \code{G} converged fits, the function warns and summarizes the
#' fits it has. The number of data sets drawn and set aside is the
#' element \code{n_redrawn} of the \code{simulation} attribute, an
#' integer read with \code{attr(x, "simulation")$n_redrawn}: 0 when
#' every fit converged, and the total number of data sets drawn is
#' \code{n_redrawn} plus the number of converged fits. The other
#' elements of the record keep their meaning: \code{attempted} is
#' \code{G}, \code{returned} the replications with an interval, and
#' \code{failed} those without one (a converged fit whose variance for
#' the path was negative or missing, or a replication never obtained).
#' A large \code{n_redrawn} beside \code{G} says the model is hard to
#' fit at this sample size, and the summaries then describe only the
#' samples in which it could be fitted.
#'
#' @inheritSection ss_aipe_cv_sensitivity The record of the simulation
#' @export

ss_aipe_sem_path_sensitivity <- function(model, supposed_cov_matrix,
                                         true_cov_matrix = supposed_cov_matrix,
                                         which_path,
                                         width, specified_N = NULL, conf_level = 0.95, assurance = NULL,
                                         G = 100, filename = NULL, ...,
                                         seed = NULL) {
  .dmar_check_lavaan_dots(..., fn = "ss_aipe_sem_path_sensitivity")
  if (!requireNamespace("MASS", quietly = TRUE)) .dmar_refuse(
    "The package 'MASS' is needed; please install the package and try again.",
    "package_required", package = "MASS")
  if (!requireNamespace("lavaan", quietly = TRUE)) .dmar_refuse(
    "The package 'lavaan' is needed; please install the package and try again.",
    "package_required", package = "lavaan")
  .check_filename(filename)
  G <- .dmar_check_count(G, "G")
  .dmar_local_seed(seed)

  # A borderline sample makes lavaan warn through the per-replication fit:
  # that the optimizer has not found a solution, that an estimated observed
  # variable variance is negative, that the gradient at the reported
  # solution is not near zero, or that a poor marker item was swapped for
  # another. Nonconvergence is handled by the convergence check in the loop,
  # which discards the replication and draws a fresh sample; the others
  # describe one sample's fit, and the summary statistics absorb them as
  # sampling variability. Each is noise at the level of the Monte Carlo
  # study, so every warning lavaan raises inside this one fit (all carry
  # lavaan's prefix) is muffled here. Warnings from the planning step, the
  # population fit, and the convergence tally below reach the caller. The
  # arguments given in '...' go to this fit, as to the planning fit and to
  # the fit to true_cov_matrix, so that all three estimate the path in one
  # parameterization.
  fit_one <- function(Data) {
    withCallingHandlers(
      lavaan::sem(model, data = as.data.frame(Data), ...),
      warning = function(w) {
        if (grepl("^lavaan", conditionMessage(w))) {
          invokeRestart("muffleWarning")
        }
      })
  }

  # An analysis model whose likelihood has no maximum at supposed_cov_matrix is
  # refused by the planner (.no_max_check()), and the refusal is raised
  # again under this function's name, with the matrix named as this
  # function's argument.
  result_plan <- tryCatch(ss_aipe_sem_path(
    model = model, true_cov_matrix = supposed_cov_matrix, width = width,
    which_path = which_path, conf_level = conf_level, assurance = assurance, internal = TRUE, ...
  ), dmar_no_maximum = function(e) {
    .no_max_refuse("ss_aipe_sem_path_sensitivity", .no_max_fields(e),
                   where = "the analysis model fit to supposed_cov_matrix",
                   placeholder = paste("the matrix given as supposed_cov_matrix in",
                                       "place of my_S"),
                   defined = "keep")
  })

  obs_vars <- result_plan$obs_vars
  alpha <- 1 - conf_level
  p <- dim(supposed_cov_matrix)[1]
  N <- if (is.null(specified_N)) result_plan$necessary_N else specified_N

  # The population value of the target, from the analysis model fit to
  # true_cov_matrix, before any replication is drawn: a fit whose likelihood has
  # no maximum there, an error variance running off below zero along a
  # ridge with the likelihood flat or still rising (.lavaan_no_maximum(),
  # R/no_maximum_internals.R), has no population value to cover, and the
  # study is refused with the lavaan code of the fit.
  true_fit <- lavaan::sem(model, sample.cov = true_cov_matrix, sample.nobs = 1e5,
                          ...)
  .no_max_check(true_fit, "ss_aipe_sem_path_sensitivity", syntax = model,
                fun = "lavaan::sem",
                args = c(list(model = model, sample.cov = true_cov_matrix,
                              sample.nobs = 1e5), list(...)),
                where = "the analysis model fit to true_cov_matrix",
                placeholder = paste("the matrix given as true_cov_matrix in",
                                    "place of my_S"))
  theta_j <- unname(lavaan::coef(true_fit)[which_path])

  theta_hat_j <- rep(NA_real_, G)
  SE_theta_hat_j <- rep(NA_real_, G)

  G_attempted <- G
  g <- 0
  attempts <- 0
  max_attempts <- 20 * G
  while (g < G && attempts < max_attempts) {
    attempts <- attempts + 1
    Data <- MASS::mvrnorm(n = N, mu = rep(0, p), Sigma = true_cov_matrix)
    colnames(Data) <- obs_vars
    S_fit <- try(fit_one(Data), silent = TRUE)
    # A fit whose likelihood has no maximum (.lavaan_no_maximum()) reported
    # convergence where lavaan's optimizer stopped on a ridge, and its
    # estimate and standard error are that point's: it is discarded as a
    # fit that did not converge, and fresh data are drawn.
    if (inherits(S_fit, "try-error") ||
        !isTRUE(lavaan::lavInspect(S_fit, "converged")) ||
        !is.null(.lavaan_no_maximum(S_fit))) {
      next
    }
    g <- g + 1
    est <- lavaan::coef(S_fit)[which_path]
    v <- lavaan::vcov(S_fit)[which_path, which_path]
    theta_hat_j[g] <- unname(est)
    SE_theta_hat_j[g] <- if (is.na(v) || v < 0) NA_real_ else sqrt(v)
  } # end of while(g < G)

  if (g < G) {
    .dmar_warn(c(
      "Only ", g, " of ", G, " replications converged within ", max_attempts,
" attempts; the summary is based on the ", g, " converged replications."),
      "result_qualified")
    theta_hat_j <- theta_hat_j[seq_len(g)]
    SE_theta_hat_j <- SE_theta_hat_j[seq_len(g)]
    G <- g
  }

  CI_upper <- theta_hat_j + qnorm(1 - alpha / 2) * SE_theta_hat_j
  CI_lower <- theta_hat_j - qnorm(1 - alpha / 2) * SE_theta_hat_j
  w <- CI_upper - CI_lower

  # Per-replication detail at full length (rows with NA are replications whose
  # standard error could not be computed); built before trimming so every
  # column has the same length.
  result <- data.frame(
    theta_hat = theta_hat_j, se_theta_hat = SE_theta_hat_j,
    ci_lower = CI_lower, ci_upper = CI_upper, width = w
  )

  if (!is.null(filename)) {
    if (file.exists(filename) && file.size(filename) > 0) {
      message("The file '", filename, "' already exists; the simulation results are appended to it.")
      utils::write.table(result, filename, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
    } else {
      utils::write.table(result, filename, sep = ",", row.names = FALSE, append = FALSE)
    }
  }

  CI_upper <- na.omit(CI_upper)
  CI_lower <- na.omit(CI_lower)
  w <- na.omit(w)
  G <- length(w)

  percent_narrower <- sum(w <= width) / G
  alpha_emp_upper <- sum(theta_j > CI_upper) / G
  alpha_emp_lower <- sum(theta_j < CI_lower) / G

  # Summary

  out <- data.frame(
    term = c(
      "mean_path", "median_path", "sd_path",
      "mean_ci_width", "median_ci_width", "sd_ci_width",
      "pct_ci_less_width",
      "pct_ci_miss_low", "pct_ci_miss_high", "total_type_I_error",
      "total_N", "true_path",
      "width", "conf_level",
      if (!is.null(assurance)) "assurance"
    ),
    value = c(
      mean(theta_hat_j, na.rm = TRUE), median(theta_hat_j, na.rm = TRUE),
      sd(theta_hat_j, na.rm = TRUE),
      mean(w), median(w), sd(w),
      percent_narrower,
      alpha_emp_lower, alpha_emp_upper, alpha_emp_lower + alpha_emp_upper,
      N, theta_j,
      width, conf_level,
      if (!is.null(assurance)) assurance
    )
  )
  # The record counts the G replications asked for: those with an interval,
  # and those without one (a negative or missing variance for the path, or
  # a replication never obtained within max_attempts). A data set whose fit
  # did not converge was drawn again; n_redrawn counts those.
  out <- .ss_sens_record(.as_dmar_tbl(out, conf_level = conf_level),
                         G_attempted, seq_len(G_attempted) <= G, conf_level)
  attr(out, "simulation")$n_redrawn <- as.integer(attempts - g)
  out
}
