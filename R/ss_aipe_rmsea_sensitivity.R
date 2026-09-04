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
#' @param filename an optional path for a comma separated file recording
#'   every converged replication (its index, the RMSEA estimate, the two
#'   confidence limits, and the interval width): nothing is written when
#'   \code{filename} is \code{NULL} (the default), a new file with a header
#'   row is created otherwise, an existing file at that path is appended to,
#'   and a throwaway run should point it at
#'   \code{tempfile(fileext = ".csv")}.
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
#' # True data generating model: two correlated factors, each measured by
#' # three standardized indicators. The factor correlation is 0.5 and every
#' # loading is 0.7. The implied population covariance matrix is assembled
#' # from the loading matrix, the factor correlation matrix, and the
#' # residual variances.
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
#' # The proposed model is fit once at a very large N to recover the
#' # population RMSEA, the sample size is planned so that the expected width
#' # of the 95 percent interval is 0.05, and a fresh sample of that size is
#' # drawn and fit on every replication. Notice that true_rmsea is about
#' # 0.20, since a single factor is a poor description of two-factor data,
#' # that the realized widths sit close to the target, and that
#' # pct_ci_less_w is near one half, which is what planning for the expected
#' # width delivers. G = 20 keeps the example quick; a reported sensitivity
#' # study deserves the default G = 200 or more.
#' set.seed(113)
#' ss_aipe_rmsea_sensitivity(width = 0.05, model = proposed, Sigma = Sigma,
#'                           G = 20)
#'
#' @export
ss_aipe_rmsea_sensitivity <- function(width, model, Sigma, N = NULL, conf_level = 0.95, G = 200,
                                      filename = NULL, ...) {
  if (!requireNamespace("MASS", quietly = TRUE)) stop("The package 'MASS' is needed; please install the package and try again.")
  if (!requireNamespace("lavaan", quietly = TRUE)) stop("The package 'lavaan' is needed; please install the package and try again.")

  if (is.null(rownames(Sigma)) || is.null(colnames(Sigma))) {
    stop("'Sigma' must have row and column names that match the observed variables in 'model'.", call. = FALSE)
  }
  .check_filename(filename)

  # A borderline sample makes lavaan warn through the per-replication fit:
  # that the optimizer has not found a solution, that an estimated observed
  # variable variance is negative, that the gradient at the reported
  # solution is not near zero, or that a poor marker item was swapped for
  # another. Nonconvergence is handled by the convergence check in the loop,
  # which skips the replication; the others describe one sample's fit, and
  # the summary statistics absorb them as sampling variability. Each is
  # noise at the level of the Monte Carlo study, so every warning lavaan
  # raises inside this one fit (all carry lavaan's prefix) is muffled here.
  # Warnings from the population fit, the sample size planning, and
  # ci_rmsea() reach the caller.
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
  M_fit <- lavaan::sem(model, sample.cov = Sigma, sample.nobs = 1e6, ...)
  rmsea <- unname(lavaan::fitMeasures(M_fit, "rmsea"))
  df    <- unname(lavaan::fitMeasures(M_fit, "df"))

  if (is.null(N)) N <- ss_aipe_rmsea(RMSEA = rmsea, df = df, width = width, conf_level = conf_level)[1, 2]
  p <- dim(Sigma)[1]
  rmsea_hat <- rep(NA_real_, G)
  CI_upper  <- rep(NA_real_, G)
  CI_lower  <- rep(NA_real_, G)

  if (!is.null(filename)) {
    if (file.exists(filename) && file.size(filename) > 0) {
      message("The file '", filename, "' already exists; the simulation results are appended to it.")
    } else {
      res_col_names <- data.frame("iteration", "rmsea_hat", "ci_low", "ci_up", "width")
      utils::write.table(res_col_names, filename, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
    }
  }

  for (g in seq_len(G)) {
    Data <- MASS::mvrnorm(n = N, mu = rep(0, p), Sigma = Sigma)
    colnames(Data) <- rownames(Sigma)

    m_fit <- try(fit_one(Data), silent = TRUE)
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
