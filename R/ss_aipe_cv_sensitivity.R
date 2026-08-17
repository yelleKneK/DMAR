#' Sensitivity Analysis for Sample Size Planning From the Accuracy in Parameter Estimation Perspective for the Coefficient of Variation
#'
#' @description
#' Quantifies how much misspecification of the population coefficient of variation can distort an AIPE-based
#' sample size plan. Given a true (population) value of the coefficient of variation and the value that was
#' used in planning, the function simulates draws of size \emph{N} from a normal population, computes the
#' confidence interval for the coefficient of variation on each replication, and summarizes how often the
#' realized interval width is below the desired target and how often the interval covers the population
#' value. This is the standard sensitivity-analysis workflow described in Kelley (2007) and Maxwell,
#' Delaney, and Kelley (2027, Section 3.11 on sample size planning).
#'
#' @param true_cv Population coefficient of variation (the data generating value)
#' @param estimated_cv Coefficient of variation used to plan the study (the value the researcher
#'   guessed when invoking \code{\link{ss_aipe_cv}}); must be positive. Supply this or \code{specified_N}
#'   but not both.
#' @param width Desired (full) width of the two-sided confidence interval for the population coefficient
#'   of variation
#' @param assurance Probability with which the realized interval should be no wider than \code{width}
#'   (must be \code{NULL} or strictly between 0 and 1; \code{NULL} means plan to the \emph{expected} width
#'   without an assurance constraint)
#' @param mean Population mean used by the simulator to generate data (the standard deviation is
#'   determined by \code{mean} and \code{true_cv}, since CV = sigma/mu). Default 100.
#' @param specified_N Pre-specified sample size to evaluate (use this when you want the sensitivity
#'   results at a fixed \emph{N} rather than at the \emph{N} that \code{\link{ss_aipe_cv}} would
#'   recommend); incompatible with \code{estimated_cv}.
#' @param conf_level Desired confidence level (i.e., 1 - Type I error rate); default 0.95
#' @param G Number of Monte Carlo replications; defaults to 1000. Increase (e.g., 5000 or 10000) for
#'   stable Type I error estimates.
#' @param print_iter Logical. If \code{TRUE} the simulation prints the iteration index after each
#'   replication (helpful for long runs); default \code{FALSE}.
#' @param save Logical. If \code{TRUE} the per-replication results are appended to a CSV file at
#'   \code{filename}; default \code{FALSE}.
#' @param filename Path used when \code{save = TRUE}; default
#'   \code{"ss_aipe_cv_sensitivity_result.csv"} in the current working directory.
#'
#' @details
#' Sample size planning for the coefficient of variation under the Accuracy in Parameter Estimation
#' framework chooses \emph{N} so that the expected (or, with assurance, the high-probability) confidence
#' interval width is no larger than \code{width} (Kelley, 2007). Because the procedure assumes the
#' planning value \code{estimated_cv} matches the population value \code{true_cv}, in practice
#' the realized width will deviate from the planned width whenever the planning value is wrong. This
#' sensitivity analysis quantifies the deviation by Monte Carlo simulation: the planned \emph{N} is
#' obtained from \code{\link{ss_aipe_cv}} with \code{estimated_cv}, then samples are drawn from the
#' \emph{true} population (with coefficient of variation \code{true_cv}) and the realized confidence
#' interval widths are summarized.
#'
#' For a discussion of AIPE-based sample size planning more generally and how sensitivity analyses guard
#' against misspecification, see Maxwell, Delaney, & Kelley (2027, Section 3.5).
#'
#' @return
#' A \code{data.frame} with columns \code{term} and \code{value}
#' summarizing the Monte Carlo results across the \code{G} replications.
#' The \code{term} entries are: \code{"mean_cv"}, \code{"median_cv"},
#' \code{"sd_cv"} (mean / median / SD of the \emph{G} observed sample
#' coefficients of variation); \code{"mean_ci_width"},
#' \code{"median_ci_width"}, \code{"sd_ci_width"} (corresponding
#' summaries of the realized interval widths); \code{"pct_ci_less_w"}
#' (proportion of intervals at or below the planning width \code{width});
#' \code{"pct_ci_miss_low"} and \code{"pct_ci_miss_high"} (tail-specific
#' non-coverage); \code{"total_type_I_error"} (overall empirical
#' non-coverage of \code{true_cv}); plus the input echoes
#' \code{"total_N"} (the sample size evaluated), \code{"true_cv"},
#' \code{"estimated_cv"} (NA when \code{specified_N} was supplied
#' instead), \code{"width"}, \code{"conf_level"}, and
#' \code{"assurance"} (present only when an assurance was supplied). The
#' proportion rows are on the 0 to 1 scale, not percentages, so
#' \code{total_type_I_error} is the sum of \code{pct_ci_miss_low} and
#' \code{pct_ci_miss_high}.
#'
#' @references
#' Chattopadhyay, B., & Kelley, K. (2016). Estimation of the coefficient of variation with minimum
#' risk: A sequential method for minimizing sampling error and study cost. \emph{Multivariate
#' Behavioral Research, 51}(5), 627--648. \doi{10.1080/00273171.2016.1203279}
#'
#' Kelley, K. (2007). Sample size planning for the coefficient of variation from the accuracy in
#' parameter estimation approach. \emph{Behavior Research Methods, 39}(4), 755--766.
#'   \doi{10.3758/BF03192966}
#'
#' Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). \emph{Designing
#'   experiments and analyzing data: A model comparison perspective}
#'   (4th ed.). Routledge. (See Chapter 3.)
#'
#' @examples
#' # Textbook scenario (Kelley, 2007). A researcher plans a study to estimate the
#' # coefficient of variation of reaction times with a 95% confidence interval no
#' # wider than .10. They guess from prior work that the population coefficient of
#' # variation is around .25 and apply ss_aipe_cv() to obtain a planned N.
#'
#' # Question 1: how does the realized interval width behave when the planning
#' # value is correct (well-specified case)?
#' set.seed(113)
#' ss_aipe_cv_sensitivity(
#'   true_cv      = .25,
#'   estimated_cv = .25,
#'   width            = .10,
#'   assurance        = NULL,
#'   conf_level       = .95,
#'   G                = 200,
#'   print_iter       = FALSE
#' )
#'
#' # Question 2: what happens if the planning value is materially smaller than
#' # the true coefficient of variation (a common direction of misspecification,
#' # since planning values are often optimistic)? The intervals will be wider on
#' # average than the target and the pct_ci_less_w will fall.
#' set.seed(113)
#' ss_aipe_cv_sensitivity(
#'   true_cv      = .35,
#'   estimated_cv = .25,
#'   width            = .10,
#'   assurance        = NULL,
#'   conf_level       = .95,
#'   G                = 200,
#'   print_iter       = FALSE
#' )
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso \code{\link{ss_aipe_cv}}, \code{\link{ci_cv}}, \code{\link{cv}}
#'
#' @keywords design htest
#'
#' @seealso \code{\link{design_consequences}} for what a chosen design delivers:
#'   power, the Type S (sign) and Type M (exaggeration) errors of the
#'   significance filter, and the expected confidence interval width.
#'
#' @export
ss_aipe_cv_sensitivity <- function(true_cv = NULL, estimated_cv = NULL, width = NULL, assurance = NULL, mean = 100,
                                   specified_N = NULL, conf_level = .95, G = 1000, print_iter = FALSE,
                                   save = FALSE, filename = "ss_aipe_cv_sensitivity_result.csv") {
  if (is.null(estimated_cv) && is.null(specified_N)) {
    stop("You must specify either 'estimated_cv' or 'specified_N'.", call. = FALSE)
  }
  if (!is.null(estimated_cv) && !is.null(specified_N)) {
    stop("You must specify 'estimated_cv' or 'specified_N', but not both.", call. = FALSE)
  }

  if (!is.null(estimated_cv)) {
    if (estimated_cv <= 0) {
      stop("'estimated_cv' must be positive (see Chattopadhyay & Kelley, 2016).", call. = FALSE)
    }
    N <- ss_aipe_cv(C_of_V = estimated_cv, mu = NULL, sigma = NULL, width = width, conf_level = conf_level, alpha_lower = NULL, alpha_upper = NULL, assurance = assurance)[1, 2]
  } else {
    N <- specified_N
  }

  if (is.null(true_cv) || true_cv <= 0) {
    stop("'true_cv' must be specified and positive.", call. = FALSE)
  }

  CN <- c("lower_limit", "upper_limit", "cv", "int_ok", "Width")
  Results <- matrix(NA, G, length(CN))
  colnames(Results) <- CN

  for (i in 1:G)
  {
    if (print_iter == TRUE) cat(c(i), "\n")
    X <- rnorm(N, mean = mean, sd = true_cv * mean)
    CI_for_CV <- ci_cv(data = X, conf_level = conf_level)

    Results[i, 1] <- CI_for_CV$value[CI_for_CV$term == "lower_limit"]
    Results[i, 2] <- CI_for_CV$value[CI_for_CV$term == "upper_limit"]
    Results[i, 3] <- CI_for_CV$value[CI_for_CV$term == "c_of_v"]

    Results[i, 4] <- sum((Results[i, 1] <= true_cv) & (true_cv <= Results[i, 2]))
    Results[i, 5] <- Results[i, 2] - Results[i, 1]
  }

  # Observed coefficients of variation.
  Obs_CV <- Results[, 3]

  Results <- as.data.frame(Results)

  if (save) {
    result_file <- filename
    # print("Simulation results will be saved to a .csv file")
    suppressWarnings(file_exist <- try(utils::read.csv(result_file), silent = TRUE))
    if (!is.null(dim(file_exist))) {
      utils::write.table(Results, result_file, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
      cat("A file in the local directory has the same name as the file where simulation", "\n", "results will be saved to. Simulation results will be appended to this file.", "\n", sep = "")
    } else {
      utils::write.table(Results, result_file, sep = ",", row.names = FALSE, append = FALSE)
    }
  }

  Summary_of_Results <- data.frame(
    term = c(
      "mean_cv", "median_cv", "sd_cv", "mean_ci_width", "median_ci_width", "sd_ci_width",
      "pct_ci_less_w", "pct_ci_miss_low", "pct_ci_miss_high", "total_type_I_error",
      "total_N", "true_cv", "estimated_cv", "width", "conf_level",
      if (!is.null(assurance)) "assurance"
    ),
    value = c(
      mean(Obs_CV), median(Obs_CV), (var(Obs_CV))^.5, mean(Results[, 2] - Results[, 1]),
      median(Results[, 2] - Results[, 1]), (var(Results[, 2] - Results[, 1]))^.5,
      mean((Results[, 2] - Results[, 1]) <= width), mean(c(true_cv <= Results[, 1])),
      mean(c(true_cv >= Results[, 2])), (mean((true_cv <= Results[, 1]) | (true_cv >= Results[, 2]))),
      N, true_cv,
      if (is.null(estimated_cv)) NA_real_ else estimated_cv,
      if (is.null(width)) NA_real_ else width,
      conf_level,
      if (!is.null(assurance)) assurance
    )
  )

  return(.as_dmar_tbl(Summary_of_Results, conf_level = conf_level))
}
