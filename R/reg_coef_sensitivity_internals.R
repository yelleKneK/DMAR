# The Monte Carlo engine shared by ss_aipe_reg_coef_sensitivity() and
# ss_aipe_rc_sensitivity(), which report the unstandardized coefficient, and
# ss_aipe_src_sensitivity(), which reports the standardized one. The name of
# the public function fixes the scale (Ken's ruling of 2026-11-11); `scale`
# here is the internal switch and is never a formal of an exported function
# (DMAR 1.0.x took a public `standardize` argument).
#
# scale = "raw": the coefficient is b_j, the jth element of
#   solve(cov_XX, cov_YX), from the true covariances for the truth the
#   intervals are judged against and from the supposed ones for the echo;
#   the replications are fit in the units of the data, and a sample size
#   planned from the supposed values is planned by ss_aipe_reg_coef() with
#   the supposed standard deviations, so that `width` is in the units of the
#   data. The rows and the per-replication columns are named reg_coef_j.
# scale = "std": every replication is standardized with its own sample
#   standard deviations before the fit, so each estimate and interval is for
#   the standardized coefficient; the truth is the population standardized
#   coefficient b_j sd_Xj / sd_Y with the population standard deviations
#   (Kelley & Maxwell, 2003, p. 316), the echo puts the supposed coefficient
#   on the same scale with the supposed standard deviations, and the plan is
#   made by ss_aipe_src(noncentral = noncentral) (unit standard deviations),
#   so that `width` is in standardized units and the plan is for the
#   interval simulated, the noncentral or the central interval for the
#   standardized coefficient. The rows and the per-replication columns are
#   named std_reg_coef_j.
#
# `fn` names the calling function in its conditions. The caller checks its
# own arguments that .dmar_check_numbers() reads from its frame and sets the
# seed; everything else is checked here.
.ss_aipe_reg_coef_sens <- function(fn, scale = c("raw", "std"),
                                   true_var_Y, true_cov_YX, true_cov_XX,
                                   supposed_var_Y, supposed_cov_YX,
                                   supposed_cov_XX, specified_N,
                                   which_predictor, width, noncentral,
                                   conf_level, assurance, G, print_iter,
                                   filename,
                                   assurance_method = c("analytic", "simulation")) {
  scale <- match.arg(scale)
  std <- scale == "std"
  # The planner's method, checked in the calling function's name: the
  # standardized coefficient offers "simulation" for both intervals
  # (ss_aipe_src()), and the coefficient in the units of the data for the
  # noncentral one only (ss_aipe_reg_coef()); see
  # R/ss_aipe_assurance_internals.R.
  assurance_method <- .ss_aipe_assurance_method(
    assurance_method, fn, simulation = std || isTRUE(noncentral == TRUE),
    where = if (!std) " with noncentral = FALSE" else "")
  if (!requireNamespace("MASS", quietly = TRUE)) .dmar_refuse(
    "The package 'MASS' is needed; please install the package and try again.",
    "package_required", package = "MASS")

  if (noncentral == TRUE && is.null(true_var_Y)) true_var_Y <- 1
  if (is.null(true_var_Y)) .dmar_refuse_argument(
    "You must specify 'true_var_Y' (i.e., the population variance of Y).",
    arg = "true_var_Y")

  if (is.null(width)) .dmar_refuse_argument(
    "You must specify 'width' (i.e., a confidence interval width).",
    arg = "width")
  if (is.null(conf_level)) .dmar_refuse_argument(
    "You must specify a confidence level (i.e., 1 - Type I error rate).",
    arg = NULL)
  if (is.null(G)) .dmar_refuse_argument(
    "You must specify 'G' (i.e., the number of generations of the simulation).",
    arg = "G")
  .check_filename(filename)
  G <- .dmar_check_count(G, "G")

  if (is.null(true_cov_XX)) .dmar_refuse_argument(
    "You must specify 'true_cov_XX' (i.e., the covariance matrix of the predictors).",
    arg = "true_cov_XX")
  if (is.null(true_cov_YX)) .dmar_refuse_argument(
    "You must specify 'true_cov_YX' (i.e., the covariance vector of the predictors with the dependent variable).",
    arg = "true_cov_YX")

  if ((sum(round(true_cov_XX, 5) == round(t(true_cov_XX), 5))) != (dim(true_cov_XX)[1] * dim(true_cov_XX)[2])) .dmar_refuse_argument(
    "The correlation matrix, \'true_cov_XX\' should be symmetric.",
    arg = "true_cov_XX")

  p <- dim(true_cov_XX)[1]
  if (!is.numeric(which_predictor) || length(which_predictor) != 1L ||
      which_predictor < 1 || which_predictor > p ||
      which_predictor != round(which_predictor)) .dmar_refuse_argument(
    sprintf("'which_predictor' must be a whole number between 1 and p = %d, the number of predictors in 'true_cov_XX'.", p),
    arg = "which_predictor")

  if (!is.null(supposed_cov_XX)) {
    if ((sum(round(supposed_cov_XX, 5) == round(t(supposed_cov_XX), 5))) != (dim(supposed_cov_XX)[1] * dim(supposed_cov_XX)[2])) .dmar_refuse_argument(
      "The covariance matrix, \'supposed_cov_XX\' should be symmetric.",
      arg = "supposed_cov_XX")
  }

  if (is.null(supposed_var_Y)) supposed_var_Y <- true_var_Y
  if (is.null(supposed_cov_XX)) supposed_cov_XX <- true_cov_XX
  if (is.null(supposed_cov_YX)) supposed_cov_YX <- true_cov_YX

  Supposed_Sigma <- cbind(c(supposed_var_Y, supposed_cov_YX), rbind(supposed_cov_YX, supposed_cov_XX))
  sigma_Y <- sqrt(unname(Supposed_Sigma[1, 1]))
  sigma_X <- sqrt(unname(Supposed_Sigma[(1 + which_predictor), (1 + which_predictor)]))
  Supposed_R2_Y_X <- (supposed_cov_YX %*% solve(supposed_cov_XX) %*% supposed_cov_YX) / (sigma_Y^2)
  Supposed_b_j <- (solve(supposed_cov_XX) %*% supposed_cov_YX)[which_predictor]

  # Covariance structure.
  True_Sigma <- cbind(c(true_var_Y, true_cov_YX), rbind(true_cov_YX, true_cov_XX))

  True_R2_Y_X <- (True_Sigma[1, -1] %*% solve(True_Sigma[-1, -1]) %*% True_Sigma[-1, 1]) / (True_Sigma[1, 1])
  True_b_j <- (solve(True_Sigma[-1, -1]) %*% True_Sigma[-1, 1])[which_predictor]

  if (True_R2_Y_X > 1) .dmar_refuse_argument(
    "You have specified an impossible correlational structure of \'true_cov_XX\' and/or \'true_cov_YX\' (the multiple R square is above 1).",
    arg = c("true_cov_XX", "true_cov_YX"))
  if (Supposed_R2_Y_X > 1) .dmar_refuse_argument(
    "You have specified an impossible correlational structure of \'supposed_cov_XX\' and/or \'supposed_cov_YX\' (the multiple R square is above 1).",
    arg = c("supposed_cov_XX", "supposed_cov_YX"))

  # The scale of the coefficient (see the head of this file). On the
  # standardized scale the true and supposed coefficients are multiplied by
  # the ratio of the standard deviations, population for the truth and
  # supposed for the echo, and the plan is made with unit standard
  # deviations; on the raw scale the plan uses the supposed standard
  # deviations.
  if (std) {
    std_true <- True_b_j * sqrt(true_cov_XX[which_predictor, which_predictor] / true_var_Y)
    std_supposed <- Supposed_b_j * sigma_X / sigma_Y
    True_b_j <- std_true
    Supposed_b_j <- std_supposed
    plan_sd_Y <- 1
    plan_sd_X <- 1
  } else {
    plan_sd_Y <- sigma_Y
    plan_sd_X <- sigma_X
  }

  if (is.null(specified_N)) {
    # The noncentral plan of ss_aipe_reg_coef() is for a standardized
    # coefficient; on the raw scale its width is in the units of the data
    # only when the supposed standard deviations are 1.
    if (noncentral == TRUE && !std &&
        !isTRUE(all.equal(c(plan_sd_Y, plan_sd_X), c(1, 1)))) .dmar_refuse_argument(
      paste0("With 'noncentral' = TRUE the sample size is planned for a standardized coefficient, but ", fn,
             "() reports the unstandardized coefficient and the supposed variances of Y and of the targeted predictor are not both 1, so the planned width would be on a different scale from the simulated intervals. Give unit supposed variances or 'specified_N', or use ss_aipe_src_sensitivity() for the standardized coefficient."),
      arg = c("noncentral", "supposed_var_Y", "supposed_cov_XX"))
    Supposed_Sigma_as_Cor <- cov2cor(Supposed_Sigma)
    supposed_cor_XX <- Supposed_Sigma_as_Cor[2:(p + 1), 2:(p + 1), drop = FALSE]
    supposed_cor_YX <- Supposed_Sigma_as_Cor[1, 2:(p + 1)]
    # The standardized coefficient is planned by ss_aipe_src() for the
    # interval simulated here: with noncentral = FALSE the central interval
    # for the standardized coefficient (Kelley & Maxwell, 2003, Equation 1),
    # whose assurance step is not the F law of the interval in the units of
    # the data that ss_aipe_reg_coef(noncentral = FALSE) uses.
    N <- if (std) {
      ss_aipe_src(
        assurance_method = assurance_method, width = width,
        true_cor_XX = supposed_cor_XX, true_cor_YX = supposed_cor_YX,
        which_predictor = which_predictor, conf_level = conf_level,
        noncentral = noncentral, assurance = assurance)[1, 2]
    } else {
      ss_aipe_reg_coef(
        assurance_method = assurance_method,
        width = width, true_cor_XX = supposed_cor_XX, true_cor_YX = supposed_cor_YX, which_predictor = which_predictor,
        conf_level = conf_level, noncentral = noncentral, assurance = assurance, true_sd_Y = plan_sd_Y, true_sd_X = plan_sd_X
      )[1, 2]
    }
  } else {
    # A fractional N, or one leaving no residual degrees of freedom, would
    # give a degenerate fit in every replication; refuse it.
    N <- .check_whole_n(specified_N, "specified_N", p + 2L)
  }

  # Means (arbitrary)
  MU <- rep(0, p + 1)

  # Begin simulation. The per-replication columns are named for the
  # coefficient computed, as the rows are.
  coef_name <- if (std) "std_reg_coef_j" else "reg_coef_j"
  Results <- matrix(NA, G, 6)
  colnames(Results) <- c(coef_name, "ci_lower", "ci_upper", "R2",
                         paste0("se_", coef_name), "t_value")
  for (i in 1:G)
  {
    if (print_iter == TRUE) cat(c(i), "\n")
    DATA <- MASS::mvrnorm(N, mu = MU, Sigma = True_Sigma)

    if (std) DATA <- scale(DATA)

    Regression_Results <- lm(DATA[, 1] ~ DATA[, -1])
    Summary_Results <- summary(Regression_Results)

    b_j <- coef(Summary_Results)[(which_predictor + 1), 1]
    SE_b_j <- coef(Summary_Results)[(which_predictor + 1), 2]

    CI_Lims <- ci_reg_coef(reg_coef_j = b_j, se_reg_coef_j = SE_b_j, sd_Y = (var(DATA[, 1]))^.5, sd_X = (var(DATA[, 1 + which_predictor]))^.5, N = dim(DATA)[1], p = (dim(DATA)[2] - 1), R2_Y_X = NULL, R2_j_X_without_j = NULL, conf_level = conf_level, R2_Y_X_without_j = NULL, t_value = NULL, alpha_lower = NULL, alpha_upper = NULL, noncentral = noncentral)

    Results[i, 1] <- b_j
    Results[i, 2] <- CI_Lims[which(CI_Lims$term == "ci_lower"), 2]
    Results[i, 3] <- CI_Lims[which(CI_Lims$term == "ci_upper"), 2]
    Results[i, 4] <- Summary_Results$r.squared
    Results[i, 5] <- SE_b_j
    Results[i, 6] <- b_j / SE_b_j
  }

  # End Simulation.

  Results <- as.data.frame(Results)

  if (!is.null(filename)) {
    suppressWarnings(file_exist <- try(utils::read.csv(filename), silent = TRUE))
    if (!is.null(dim(file_exist))) {
      utils::write.table(Results, filename, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
      cat("The file '", filename, "' already exists; the simulation results were appended to it.\n", sep = "")
    } else {
      utils::write.table(Results, filename, sep = ",", row.names = FALSE, append = FALSE)
    }
  }

  Summary <- data.frame(
    term = c(
      paste0(c("mean_", "median_", "sd_"), coef_name), "mean_ci_width", "median_ci_width", "sd_ci_width",
      "pct_ci_less_width", "pct_ci_miss_low", "pct_ci_miss_high", "total_type_I_error",
      "mean_R2", "median_R2", "sd_R2",
      "total_N", "p", "which_predictor",
      paste0(c("true_", "supposed_"), coef_name),
      "width", "conf_level",
      if (!is.null(assurance)) "assurance"
    ),
    value = c(
      mean(Results[, 1]), median(Results[, 1]), (var(Results[, 1]))^.5,
      mean(Results[, 3] - Results[, 2]), median(Results[, 3] - Results[, 2]),
      (var(Results[, 3] - Results[, 2]))^.5, mean((Results[, 3] - Results[, 2]) <= width),
      mean(True_b_j < Results[, 2]), mean(True_b_j > Results[, 3]),
      (mean((True_b_j < Results[, 2]) | (True_b_j > Results[, 3]))),
      mean(Results[, 4]), median(Results[, 4]), (var(Results[, 4]))^.5,
      N, p, which_predictor,
      True_b_j, Supposed_b_j,
      width, conf_level,
      if (!is.null(assurance)) assurance
    )
  )

  .ss_sens_record(.as_dmar_tbl(Summary, conf_level = conf_level), G,
                  !is.na(Results[, 3] - Results[, 2]), conf_level)
}
