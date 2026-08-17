# Sensitivity analysis for AIPE on a reliability coefficient (parallel-tests).
#' Sensitivity Analysis for Sample Size Planning From the AIPE Perspective for a Reliability Coefficient
#'
#' @description
#' Quantifies how much misspecification of the population reliability
#' coefficient distorts an AIPE-based sample size plan for the
#' composite-score reliability. On each replication the function
#' simulates an \emph{n} \eqn{\times} \emph{i} item-by-subject data
#' matrix from a single-factor parallel-tests model whose population
#' reliability of the sum score equals \code{true_reliability}, fits
#' the requested estimator (alpha or omega) via the corresponding
#' \code{reliability_*} function with the supplied \code{ci_method},
#' and records the realized reliability estimate and its confidence
#' interval.
#'
#' \strong{Population model.} Each item has a single common-factor
#' loading and uncorrelated unique error. With per-item variance
#' normalized to 1, the loading and unique variance are chosen so the
#' Cronbach-style sum-score reliability equals \code{true_reliability}:
#' \deqn{\lambda^2 \;=\; \frac{\rho}{i(1 - \rho) + \rho}, \qquad
#'        \psi^2 \;=\; 1 - \lambda^2,}
#' where \eqn{\rho = }\code{true_reliability} and \eqn{i} is the item
#' count. Item scores are
#' \eqn{y_{ij} = \lambda T_i + e_{ij}}, with \eqn{T_i \sim N(0, 1)} and
#' \eqn{e_{ij} \sim N(0, \psi^2)}.
#'
#' @param true_reliability Population reliability coefficient (in
#'   \eqn{[0, 1)}).
#' @param estimated_reliability Reliability used to plan the study;
#'   the function passes the implied lambda / psi^2 to
#'   \code{\link{ss_aipe_reliability}}.
#' @param i Number of items in the composite.
#' @param width Desired full width of the CI on reliability.
#' @param specified_N Sample size to evaluate (incompatible with
#'   \code{estimated_reliability}).
#' @param estimator One of \code{"alpha"} (default; coefficient alpha
#'   via \code{\link{reliability_alpha}}) or \code{"omega"} (composite
#'   reliability via \code{\link{reliability_omega}}). For a
#'   parallel-tests population the two coincide; differences in
#'   sample estimates reflect estimator-specific finite-sample bias
#'   and CI behavior.
#' @param ci_method CI method passed to the estimator. Default
#'   \code{"bonett"} for alpha and \code{"mlr"} for omega.
#' @param conf_level Confidence level (default \code{0.95}).
#' @param assurance Optional assurance probability passed to the
#'   planner.
#' @param G Number of Monte Carlo replications.
#' @param print_iter Logical.
#' @param save Logical. Save per-replication CSV.
#' @param filename Path used when \code{save = TRUE}.
#'
#' @return A \code{data.frame} with rows for mean / median / SD of
#'   the realized reliability and CI width, the proportion of intervals
#'   at or below \code{width}, tail-specific and overall non-coverage
#'   of \code{true_reliability}, and the input echoes, including \code{assurance} (present only when an
#'   assurance was supplied).
#'
#' @references
#' Kelley, K., & Pornprasertmanit, S. (2016). Confidence intervals for
#'   population reliability coefficients: Evaluation of methods,
#'   recommendations, and software for composite measures.
#'   \emph{Psychological Methods, 21}, 69--92. \doi{10.1037/a0040086}
#'
#' Terry, L. J., & Kelley, K. (2012). Sample size planning for composite
#'   reliability coefficients: Accuracy in parameter estimation via
#'   narrow confidence intervals. \emph{British Journal of Mathematical
#'   and Statistical Psychology, 65}, 371--401.
#'   \doi{10.1111/j.2044-8317.2011.02030.x}
#'
#' Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). \emph{Designing
#'   experiments and analyzing data: A model comparison perspective}
#'   (4th ed.). Routledge.
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso \code{\link{ss_aipe_reliability}}, \code{\link{reliability_alpha}}, \code{\link{reliability_omega}}
#'
#' @examples
#' # Reduced Monte Carlo sweep (small G) so the example runs quickly;
#' # raise G for a production sensitivity analysis.
#' set.seed(113)
#' ss_aipe_reliability_sensitivity(
#'   true_reliability      = 0.80,
#'   estimated_reliability = 0.80,
#'   i = 4, width = 0.15,
#'   estimator = "alpha",
#'   G = 20, print_iter = FALSE
#' )
#'
#' @keywords design htest
#'
#' @seealso \code{\link{design_consequences}} for what a chosen design delivers:
#'   power, the Type S (sign) and Type M (exaggeration) errors of the
#'   significance filter, and the expected confidence interval width.
#'
#' @family AIPE sample size planning
#'
#' @export
ss_aipe_reliability_sensitivity <- function(true_reliability = NULL,
                                            estimated_reliability = NULL,
                                            i, width,
                                            specified_N = NULL,
                                            estimator = c("alpha", "omega"),
                                            ci_method = NULL,
                                            conf_level = 0.95,
                                            assurance = NULL,
                                            G = 1000, print_iter = FALSE,
                                            save = FALSE,
                                            filename = "ss_aipe_reliability_sensitivity_result.csv") {
  estimator <- match.arg(estimator)
  if (is.null(ci_method))
    ci_method <- if (estimator == "alpha") "bonett" else "mlr"

  if (is.null(estimated_reliability) && is.null(specified_N))
    stop("You must specify either 'estimated_reliability' or 'specified_N'.", call. = FALSE)
  if (!is.null(estimated_reliability) && !is.null(specified_N))
    stop("You must specify 'estimated_reliability' or 'specified_N', but not both.", call. = FALSE)
  if (is.null(true_reliability) || !is.numeric(true_reliability) ||
      true_reliability < 0 || true_reliability >= 1)
    stop("'true_reliability' must be a single value in [0, 1).", call. = FALSE)
  if (!is.numeric(i) || length(i) != 1L || i < 2)
    stop("'i' (number of items) must be a single integer >= 2.", call. = FALSE)

  # Solve for parallel-tests lambda^2 and psi^2 that produce the requested rho
  # at unit per-item variance.
  rho_solve <- function(rho) {
    u <- rho / (i * (1 - rho) + rho)  # lambda^2
    list(lambda = sqrt(u), psi_square = 1 - u)
  }

  if (!is.null(estimated_reliability)) {
    pars <- rho_solve(estimated_reliability)
    plan <- suppressMessages(suppressWarnings(
      ss_aipe_reliability(model = "parallel", type = "Normal Theory",
                          width = width, i = i,
                          lambda = pars$lambda,
                          psi_square = pars$psi_square,
                          conf_level = conf_level,
                          assurance = assurance,
                          initial_iter = 50, final_iter = 200)
    ))
    n <- plan$value[plan$term == "necessary_N"]
  } else {
    n <- specified_N
  }
  n <- as.integer(n)
  if (n < 4) stop("Resolved sample size is < 4.", call. = FALSE)

  pars_true <- rho_solve(true_reliability)
  lambda_t  <- pars_true$lambda
  psi_t     <- sqrt(pars_true$psi_square)

  rel_hat   <- numeric(G)
  ci_lo     <- numeric(G)
  ci_hi     <- numeric(G)
  ci_w      <- numeric(G)
  tI_lower  <- logical(G)
  tI_upper  <- logical(G)

  for (g in seq_len(G)) {
    if (isTRUE(print_iter)) cat(g, "\n")
    T_i <- stats::rnorm(n, mean = 0, sd = 1)
    E   <- matrix(stats::rnorm(n * i, mean = 0, sd = psi_t), nrow = n, ncol = i)
    Y   <- matrix(lambda_t * T_i, nrow = n, ncol = i, byrow = FALSE) + E
    Y_df <- as.data.frame(Y)

    res <- tryCatch({
      if (estimator == "alpha") {
        reliability_alpha(data = Y_df, ci_method = ci_method,
                          conf_level = conf_level, B = 200)
      } else {
        suppressMessages(suppressWarnings(
          reliability_omega(data = Y_df, ci_method = ci_method,
                            conf_level = conf_level, B = 200)
        ))
      }
    }, error = function(e) NULL)

    if (is.null(res)) {
      r <- NA_real_; lo <- NA_real_; hi <- NA_real_
    } else {
      # reliability_alpha() / reliability_omega() return the point estimate
      # under the term "estimate" (not "alpha" or "omega"), with
      # "lower_limit" / "upper_limit" rows.
      pick <- function(needle) {
        hit <- res$value[res$term == needle]
        if (length(hit) == 0) NA_real_ else as.numeric(hit[1])
      }
      r  <- pick("estimate")
      lo <- pick("lower_limit")
      hi <- pick("upper_limit")
    }
    rel_hat[g]  <- r
    ci_lo[g]    <- lo
    ci_hi[g]    <- hi
    ci_w[g]     <- hi - lo
    tI_lower[g] <- isTRUE(true_reliability < lo)
    tI_upper[g] <- isTRUE(true_reliability > hi)
  }

  if (isTRUE(save)) {
    per_rep <- data.frame(rel_hat = rel_hat, ci_lower = ci_lo,
                          ci_upper = ci_hi, ci_width = ci_w,
                          type_I_lower = tI_lower, type_I_upper = tI_upper)
    .write_sensitivity_csv(per_rep, filename)
  }

  out <- data.frame(
    term  = c("mean_reliability", "median_reliability", "sd_reliability",
              "mean_ci_width", "median_ci_width", "sd_ci_width",
              "pct_ci_less_w",
              "pct_ci_miss_low", "pct_ci_miss_high", "total_type_I_error",
              "total_N", "items",
              "true_reliability", "estimated_reliability",
              "width", "conf_level",
              if (!is.null(assurance)) "assurance"),
    value = c(mean(rel_hat, na.rm = TRUE),
              stats::median(rel_hat, na.rm = TRUE),
              stats::sd(rel_hat, na.rm = TRUE),
              mean(ci_w, na.rm = TRUE),
              stats::median(ci_w, na.rm = TRUE),
              stats::sd(ci_w, na.rm = TRUE),
              mean(ci_w <= width, na.rm = TRUE),
              mean(tI_lower, na.rm = TRUE),
              mean(tI_upper, na.rm = TRUE),
              mean(tI_lower | tI_upper, na.rm = TRUE),
              n, i,
              true_reliability,
              if (is.null(estimated_reliability)) NA_real_ else estimated_reliability,
              width, conf_level,
              if (!is.null(assurance)) assurance),
    stringsAsFactors = FALSE, row.names = NULL
  )
  .as_dmar_tbl(out, conf_level = conf_level)
}
