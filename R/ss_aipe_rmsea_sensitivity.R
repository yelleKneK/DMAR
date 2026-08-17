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
#'   model must match the row and column names of \code{Sigma}.
#' @param Sigma the true population covariance matrix, which is used to generate
#'   random data for the simulation study. The row and column names of
#'   \code{Sigma} must match the observed variables in \code{model}.
#' @param N if \code{N} is specified, random samples of the specified size are
#'   generated. Otherwise the sample size is calculated with the sample size
#'   planning method so that the expected width of a confidence interval for the
#'   population RMSEA is no larger than \code{width}.
#' @param conf_level confidence level (i.e., 1 - the Type I error rate).
#' @param G number of replications in the Monte Carlo simulation.
#' @param save option to save simulation results. With \code{save = TRUE} the
#'   per-replication results are written to \code{filename}.
#' @param filename the name of the file that simulation results are saved to.
#' @param \dots additional arguments passed to \code{\link[lavaan]{sem}} when
#'   fitting the model (for example \code{estimator} or \code{missing}).
#'
#' @details
#' This function implements the sample size planning method proposed in Kelley
#' and Lai (2011). It uses \code{\link[lavaan]{sem}} to fit the proposed model
#' to the population covariance matrix, which recovers the population RMSEA (the
#' model misspecification) and the model degrees of freedom, and to fit the
#' model to each simulated sample, and it uses \code{\link{ci_rmsea}} to
#' construct the confidence interval for the population RMSEA in each
#' replication. The model is specified in \pkg{lavaan} syntax, so \pkg{lavaan}
#' must be installed.
#'
#' Earlier versions of this function used the \pkg{sem} package to fit the
#' model. The fit is now carried out with \pkg{lavaan}, the structural equation
#' modeling backend used throughout DMAR. The population RMSEA is read from
#' \code{lavaan::fitMeasures()}, which is computed reliably for the large-sample
#' population fit.
#'
#' @return
#' A \code{data.frame} with columns \code{term} and \code{value}
#' summarizing the a priori Monte Carlo study. The \code{term} entries
#' are: \code{"mean_rmsea"}, \code{"median_rmsea"}, \code{"sd_rmsea"}
#' (summaries of the realized RMSEA estimates across the converged
#' replications); \code{"mean_ci_width"}, \code{"median_ci_width"},
#' \code{"sd_ci_width"} (summaries of the realized interval widths);
#' \code{"pct_ci_less_w"} (proportion of intervals narrower than the
#' target width); \code{"pct_ci_miss_low"} and \code{"pct_ci_miss_high"}
#' (tail-specific empirical non-coverage of the population RMSEA);
#' \code{"total_type_I_error"} (overall empirical non-coverage, the sum
#' of the two tails); and the echoes \code{"suc_rep"} (number of
#' converged replications), \code{"total_N"} (the \emph{N} evaluated),
#' \code{"df"} (model degrees of freedom), \code{"true_rmsea"} (the
#' population RMSEA recovered from fitting \code{model} to
#' \code{Sigma}), \code{"width"}, and \code{"conf_level"}. The
#' proportion rows are on the 0 to 1 scale, not percentages.
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
#' as \code{suc_rep}. Increase \code{G} if many replications fail to converge.
#'
#' @seealso \code{\link[lavaan]{sem}}, \code{\link{ss_aipe_rmsea}},
#'   \code{\link{ci_rmsea}}
#'
#' @examples
#' set.seed(113)
#'
#' # True data generating model: two correlated factors (r = 0.5), three
#' # standardized indicators each (loadings 0.7). Build the implied population
#' # covariance matrix Sigma = Lambda Phi Lambda' + Psi.
#' Lambda <- matrix(0, 6, 2)
#' Lambda[1:3, 1] <- 0.7
#' Lambda[4:6, 2] <- 0.7
#' Phi   <- matrix(c(1, 0.5, 0.5, 1), 2, 2)
#' Sigma <- Lambda %*% Phi %*% t(Lambda) + diag(1 - 0.7^2, 6)
#' dimnames(Sigma) <- list(paste0("x", 1:6), paste0("x", 1:6))
#'
#' # Proposed (misspecified) model: a single common factor.
#' proposed <- "g =~ x1 + x2 + x3 + x4 + x5 + x6"
#'
#' # The simulation itself is not run at example time: it fits the proposed
#' # model once at a very large N to recover the population RMSEA, then
#' # generates and fits a fresh sample on every replication. The G below is
#' # already far smaller than a study one would report; the default of 200,
#' # or more, is the realistic setting. The call is:
#' # ss_aipe_rmsea_sensitivity(width = 0.05, model = proposed, Sigma = Sigma,
#' #                           G = 25)
#'
#' @export
ss_aipe_rmsea_sensitivity <- function(width, model, Sigma, N = NULL, conf_level = 0.95, G = 200, save = FALSE,
                                      filename = "ss_aipe_rmsea_sensitivity_result.csv", ...) {
  if (!requireNamespace("MASS", quietly = TRUE)) stop("The package 'MASS' is needed; please install the package and try again.")
  if (!requireNamespace("lavaan", quietly = TRUE)) stop("The package 'lavaan' is needed; please install the package and try again.")

  if (is.null(rownames(Sigma)) || is.null(colnames(Sigma))) {
    stop("'Sigma' must have row and column names that match the observed variables in 'model'.", call. = FALSE)
  }

  prev_warn <- getOption("warn")
  on.exit(options(warn = prev_warn), add = TRUE)
  options(warn = -1)

  # Population fit: fit the proposed (possibly misspecified) model to the true
  # covariance matrix at a very large N to recover the population RMSEA and the
  # model degrees of freedom. lavaan::fitMeasures() returns the population RMSEA
  # reliably for this large-sample fit.
  M_fit <- lavaan::sem(model, sample.cov = Sigma, sample.nobs = 1e6, ...)
  rmsea <- unname(lavaan::fitMeasures(M_fit, "rmsea"))
  df    <- unname(lavaan::fitMeasures(M_fit, "df"))

  if (is.null(N)) N <- ss_aipe_rmsea(RMSEA = rmsea, df = df, width = width, conf_level = conf_level)[1, 2]
  p <- dim(Sigma)[1]
  rmsea_hat <- rep(NA_real_, G)
  CI_upper  <- rep(NA_real_, G)
  CI_lower  <- rep(NA_real_, G)

  res_col_names <- data.frame("iteration", "rmsea_hat", "ci_low", "ci_up", "width")
  result_file <- filename
  if (save) {
    suppressWarnings(file_exist <- try(utils::read.csv(result_file), silent = TRUE))
    if (!is.null(dim(file_exist))) {
      cat("A file in the local directory has the same name as the file where simulation", "\n", "results will be saved to. Simulation results will be appended to this file.", "\n", sep = "")
    } else {
      utils::write.table(res_col_names, result_file, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
    }
  }

  for (g in seq_len(G)) {
    Data <- MASS::mvrnorm(n = N, mu = rep(0, p), Sigma = Sigma)
    colnames(Data) <- rownames(Sigma)

    m_fit <- try(lavaan::sem(model, data = as.data.frame(Data), ...), silent = TRUE)
    if (inherits(m_fit, "try-error") || !isTRUE(lavaan::lavInspect(m_fit, "converged"))) {
      next
    }
    rh <- unname(lavaan::fitMeasures(m_fit, "rmsea"))
    if (length(rh) != 1L || is.na(rh)) {
      next
    }
    rmsea_hat[g] <- rh
    CI <- ci_rmsea(rmsea_hat[g], df = df, N = N, conf_level = conf_level)
    CI_lower[g] <- CI[1, 2]
    CI_upper[g] <- CI[3, 2]
    if (save) {
      sim_result <- cbind(g, rmsea_hat[g], CI_lower[g], CI_upper[g], CI_upper[g] - CI_lower[g])
      utils::write.table(sim_result, result_file, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
    }
  }

  w <- CI_upper - CI_lower
  suc_rep <- sum(!is.na(w))

  term <- c(
    "mean_rmsea", "median_rmsea", "sd_rmsea",
    "mean_ci_width", "median_ci_width", "sd_ci_width",
    "pct_ci_less_w", "pct_ci_miss_low", "pct_ci_miss_high", "total_type_I_error",
    "suc_rep", "total_N", "df", "true_rmsea", "width", "conf_level"
  )
  value <- c(
    mean(rmsea_hat, na.rm = TRUE), median(rmsea_hat, na.rm = TRUE), sd(rmsea_hat, na.rm = TRUE),
    mean(w, na.rm = TRUE), median(w, na.rm = TRUE), sd(w, na.rm = TRUE),
    sum(w < width, na.rm = TRUE) / suc_rep,
    sum(CI_lower > rmsea, na.rm = TRUE) / suc_rep,
    sum(CI_upper < rmsea, na.rm = TRUE) / suc_rep,
    sum(CI_lower > rmsea, na.rm = TRUE) / suc_rep + sum(CI_upper < rmsea, na.rm = TRUE) / suc_rep,
    suc_rep, N, df, rmsea, width, conf_level
  )

  return(.as_dmar_tbl(data.frame(term, value), conf_level = conf_level))
}
