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
#' @param est_Sigma the covariance matrix used to calculate sample size, may or may not be the true covariance matrix. The row names and column names of \code{est_Sigma} should be the same as the observed variables in \code{model}
#' @param true_Sigma the true population covariance matrix, which will be used to generate random data for the simulation study. The row names and column names of \code{true_Sigma} should be the same as the observed variables in \code{model}
#' @param which_path the parameter label of the targeted path, given as a
#'   character string, for example \code{"b"} for the path labeled
#'   \code{f2 ~ b*f1} in \code{model}
#' @param desired_width desired confidence interval width for the model parameter of interest
#' @param N the sample size of random data. If it is \code{NULL}, it will be determined by the sample size planning method
#' @param conf_level confidence level (i.e., 1- Type I error rate)
#' @param assurance the assurance that the confidence interval obtained in a particular study will be no wider than desired (must be \code{NULL} or a value between 0.50 and 1)
#' @param G number of replications in the Monte Carlo simulation
#' @param save option to save simulation results. It can be saved with \code{save = TRUE} outside of the printed results
#' @param filename the name of the file that simulation results will be saved to
#' @param \dots allows one to potentially include parameter values for inner functions
#'
#' @details
#' This function implements the sample size planning methods proposed in Lai
#' and Kelley (2011). It calls \code{\link{ss_aipe_sem_path}} to plan the
#' sample size and to identify the targeted path, then fits the analysis model
#' to data simulated from \code{true_Sigma} with \code{\link[lavaan]{sem}}. The
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
#' realized interval widths); \code{"pct_ci_less_w"} (proportion of
#' realized widths at or below \code{desired_width});
#' \code{"pct_ci_miss_low"} and \code{"pct_ci_miss_high"}
#' (tail-specific empirical non-coverage of the population path);
#' \code{"total_type_I_error"} (overall empirical non-coverage, the sum
#' of the two tails); and the echoes \code{"suc_rep"} (number of
#' converged replications), \code{"total_N"} (the \emph{N} evaluated),
#' \code{"true_path"} (the population value of the targeted path under
#' \code{true_Sigma}), \code{"width"}, \code{"conf_level"}, and
#' \code{"assurance"} (present only when an assurance was supplied). The
#' proportion rows are on the 0 to 1 scale, not percentages.
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
#' continues. A safety cap stops the loop after \code{20 * G} attempts, and a
#' single warning is issued if fewer than \code{G} replications converged.
#'
#' @seealso
#' \code{\link[lavaan]{sem}}, \code{\link{cov_sem}},
#' \code{\link{ss_aipe_sem_path}}
#'
#' @examples
#' # This function is itself a Monte Carlo study: it plans a sample size and
#' # then fits the analysis model to G freshly simulated data sets, so even a
#' # modest G takes long enough that the worked example below is shown here
#' # rather than run.
#' #
#' # The planning values a researcher would bring to ss_aipe_sem_path():
#' # each factor is measured by three indicators, each with a residual
#' # variance of 0.5.
#' #   planning_model <- "
#' #     f1 =~ 1*y1 + 0.8*y2 + 0.8*y3
#' #     f2 =~ 1*y4 + 0.8*y5 + 0.8*y6
#' #     f2 ~ 0.5*f1
#' #     f1 ~~ 1*f1
#' #     f2 ~~ 0.75*f2
#' #     y1 ~~ 0.5*y1; y2 ~~ 0.5*y2; y3 ~~ 0.5*y3
#' #     y4 ~~ 0.5*y4; y5 ~~ 0.5*y5; y6 ~~ 0.5*y6
#' #   "
#' #
#' # The population the study will actually sample from: the same structural
#' # path of 0.5, but noisier indicators than the planning values assumed,
#' # with residual variances of 0.8.
#' #   true_model <- "
#' #     f1 =~ 1*y1 + 0.8*y2 + 0.8*y3
#' #     f2 =~ 1*y4 + 0.8*y5 + 0.8*y6
#' #     f2 ~ 0.5*f1
#' #     f1 ~~ 1*f1
#' #     f2 ~~ 0.75*f2
#' #     y1 ~~ 0.8*y1; y2 ~~ 0.8*y2; y3 ~~ 0.8*y3
#' #     y4 ~~ 0.8*y4; y5 ~~ 0.8*y5; y6 ~~ 0.8*y6
#' #   "
#' #
#' #   analysis_model <- "
#' #     f1 =~ y1 + y2 + y3
#' #     f2 =~ y4 + y5 + y6
#' #     f2 ~ b*f1
#' #   "
#' #
#' #   est_Sigma <- cov_sem(planning_model)$sigma_theta
#' #   true_Sigma <- cov_sem(true_model)$sigma_theta
#' #
#' # The sample size planned from the optimistic measurement quality is
#' # evaluated against the population that actually holds: the realized
#' # intervals are wider than desired, and few of them meet the target.
#' #   set.seed(113)
#' #   ss_aipe_sem_path_sensitivity(model = analysis_model,
#' #                                est_Sigma = est_Sigma,
#' #                                true_Sigma = true_Sigma, which_path = "b",
#' #                                desired_width = 0.30, G = 1000)
#'
#' @keywords design multivariate
#'
#' @seealso \code{\link{design_consequences}} for what a chosen design delivers:
#'   power, the Type S (sign) and Type M (exaggeration) errors of the
#'   significance filter, and the expected confidence interval width.
#'
#' @export

ss_aipe_sem_path_sensitivity <- function(model, est_Sigma, true_Sigma = est_Sigma, which_path,
                                         desired_width, N = NULL, conf_level = 0.95, assurance = NULL,
                                         G = 100, save = FALSE, filename = "ss_aipe_sem_path_sensitivity_result.csv", ...) {
  if (!requireNamespace("MASS", quietly = TRUE)) stop("The package 'MASS' is needed; please install the package and try again.")
  if (!requireNamespace("lavaan", quietly = TRUE)) stop("The package 'lavaan' is needed; please install the package and try again.")

  # The per-replication lavaan fits emit convergence and standard-error
  # warnings on borderline samples; muffle them for the duration of the
  # Monte Carlo and restore the user's setting on exit (as in
  # ss_aipe_rmsea_sensitivity).
  prev_warn <- getOption("warn")
  on.exit(options(warn = prev_warn), add = TRUE)
  options(warn = -1)

  result_plan <- ss_aipe_sem_path(
    model = model, Sigma = est_Sigma, desired_width = desired_width,
    which_path = which_path, conf_level = conf_level, assurance = assurance, internal = TRUE, ...
  )

  obs_vars <- result_plan$obs_vars
  alpha <- 1 - conf_level
  p <- dim(est_Sigma)[1]
  if (is.null(N)) N <- result_plan$necessary_N

  theta_hat_j <- rep(NA_real_, G)
  SE_theta_hat_j <- rep(NA_real_, G)

  g <- 0
  attempts <- 0
  max_attempts <- 20 * G
  while (g < G && attempts < max_attempts) {
    attempts <- attempts + 1
    Data <- MASS::mvrnorm(n = N, mu = rep(0, p), Sigma = true_Sigma)
    colnames(Data) <- obs_vars
    S_fit <- try(lavaan::sem(model, data = as.data.frame(Data)), silent = TRUE)
    if (inherits(S_fit, "try-error") || !isTRUE(lavaan::lavInspect(S_fit, "converged"))) {
      next
    }
    g <- g + 1
    est <- lavaan::coef(S_fit)[which_path]
    v <- lavaan::vcov(S_fit)[which_path, which_path]
    theta_hat_j[g] <- unname(est)
    SE_theta_hat_j[g] <- if (is.na(v) || v < 0) NA_real_ else sqrt(v)
  } # end of while(g < G)

  if (g < G) {
    warning("Only ", g, " of ", G, " replications converged within ", max_attempts,
      " attempts; the summary is based on the ", g, " converged replications.",
      call. = FALSE
    )
    theta_hat_j <- theta_hat_j[seq_len(g)]
    SE_theta_hat_j <- SE_theta_hat_j[seq_len(g)]
    G <- g
  }

  CI_upper <- theta_hat_j + qnorm(1 - alpha / 2) * SE_theta_hat_j
  CI_lower <- theta_hat_j - qnorm(1 - alpha / 2) * SE_theta_hat_j
  w <- CI_upper - CI_lower

  true_fit <- lavaan::sem(model, sample.cov = true_Sigma, sample.nobs = 1e5)
  theta_j <- unname(lavaan::coef(true_fit)[which_path])

  # Per-replication detail at full length (rows with NA are replications whose
  # standard error could not be computed); built before trimming so every
  # column has the same length.
  result <- data.frame(
    theta_hat = theta_hat_j, se_theta_hat = SE_theta_hat_j,
    ci_low = CI_lower, ci_up = CI_upper, width = w
  )

  if (save) {
    result_file <- filename
    suppressWarnings(file_exist <- try(utils::read.csv(result_file), silent = TRUE))
    if (!is.null(dim(file_exist))) {
      utils::write.table(result, result_file, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
      cat("A file in the local directory has the same name as the file where simulation", "\n", "results will be saved to. Simulation results will be appended to this file.", "\n", sep = "")
    } else {
      utils::write.table(result, result_file, sep = ",", row.names = FALSE, append = FALSE)
    }
  }

  CI_upper <- na.omit(CI_upper)
  CI_lower <- na.omit(CI_lower)
  w <- na.omit(w)
  G <- length(w)

  percent_narrower <- sum(w <= desired_width) / G
  alpha_emp_upper <- sum(theta_j > CI_upper) / G
  alpha_emp_lower <- sum(theta_j < CI_lower) / G

  # Summary

  out <- data.frame(
    term = c(
      "mean_path", "median_path", "sd_path",
      "mean_ci_width", "median_ci_width", "sd_ci_width",
      "pct_ci_less_w",
      "pct_ci_miss_low", "pct_ci_miss_high", "total_type_I_error",
      "suc_rep", "total_N", "true_path",
      "width", "conf_level",
      if (!is.null(assurance)) "assurance"
    ),
    value = c(
      mean(theta_hat_j, na.rm = TRUE), median(theta_hat_j, na.rm = TRUE),
      sd(theta_hat_j, na.rm = TRUE),
      mean(w), median(w), sd(w),
      percent_narrower,
      alpha_emp_lower, alpha_emp_upper, alpha_emp_lower + alpha_emp_upper,
      G, N, theta_j,
      desired_width, conf_level,
      if (!is.null(assurance)) assurance
    )
  )
  .as_dmar_tbl(out, conf_level = conf_level)
}
