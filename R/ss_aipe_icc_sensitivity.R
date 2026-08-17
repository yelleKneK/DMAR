# Sensitivity analysis for AIPE on an intraclass correlation coefficient.
#' Sensitivity Analysis for Sample Size Planning From the Accuracy in Parameter Estimation Perspective for an Intraclass Correlation Coefficient
#'
#' @description
#' Quantifies how much misspecification of the population ICC can distort an
#' AIPE-based sample size plan. Given a true (population) ICC and the value
#' used in planning, the function simulates draws of size \eqn{n \times k}
#' from the relevant variance-components model, computes the ICC and its
#' \emph{F}-distribution confidence interval on each replication via
#' \code{\link{icc}}, and summarizes how often the realized interval width
#' is below the desired target and how often the interval covers the
#' population value. This is the standard sensitivity-analysis workflow
#' described in Kelley (2007) and Maxwell, Delaney, and Kelley (2027,
#' Section 3.11 on sample size planning).
#'
#' @param true_rho Population intraclass correlation coefficient (the
#'   data generating value), at the level matching \code{type}. Must lie in
#'   \eqn{[0, 1)}.
#' @param estimated_rho ICC used to plan the study (the value the
#'   researcher guessed when invoking \code{\link{ss_aipe_icc}}). Supply
#'   this or \code{specified_N} but not both.
#' @param k Number of raters (or repeated measurements) per subject; must
#'   be at least 2.
#' @param width Desired full width of the (back-transformed) confidence
#'   interval on the ICC.
#' @param assurance Probability with which the realized interval should be
#'   no wider than \code{width} (must be \code{NULL} or strictly between
#'   0 and 1; \code{NULL} means plan to the \emph{expected} width without
#'   an assurance constraint).
#' @param specified_N Pre-specified number of subjects to evaluate (use this
#'   when you want the sensitivity results at a fixed \emph{n} rather than
#'   at the \emph{n} that \code{\link{ss_aipe_icc}} would recommend);
#'   incompatible with \code{estimated_rho}.
#' @param conf_level Desired confidence level (i.e., 1 minus the Type I
#'   error rate); default 0.95.
#' @param type Which Shrout-Fleiss (1979) ICC form is being planned. One
#'   of \code{"ICC(1,1)"} (default; one-way random model),
#'   \code{"ICC(2,1)"} (two-way random model, absolute agreement),
#'   \code{"ICC(3,1)"} (two-way mixed model, consistency), or the
#'   average-of-\eqn{k} versions \code{"ICC(1,k)"}, \code{"ICC(2,k)"},
#'   \code{"ICC(3,k)"}. The simulation generates data from the variance-
#'   components model matching the requested type, so the realized ICC
#'   on each replication is comparable to \code{true_rho}.
#' @param G Number of Monte Carlo replications; defaults to 1000.
#'   Increase (e.g., 5000 or 10000) for stable empirical-coverage
#'   estimates.
#' @param print_iter Logical. If \code{TRUE} the simulation prints the
#'   iteration index after each replication (helpful for long runs);
#'   default \code{FALSE}.
#' @param save Logical. If \code{TRUE} the per-replication results are
#'   appended to a CSV file at \code{filename}; default \code{FALSE}.
#' @param filename Path used when \code{save = TRUE}; default
#'   \code{"ss_aipe_icc_sensitivity_result.csv"} in the current working
#'   directory.
#'
#' @details
#' Sample size planning for the intraclass correlation coefficient under
#' the Accuracy in Parameter Estimation framework chooses \emph{n} so that
#' the expected (or, with assurance, the high-probability) confidence
#' interval width is no larger than \code{width} (Bonett, 2002;
#' Kelley & Maxwell, 2003). Because the procedure assumes the planning
#' value \code{estimated_rho} matches the population value
#' \code{true_rho}, in practice the realized width will deviate from the
#' planned width whenever the planning value is wrong. This sensitivity
#' analysis quantifies the deviation by Monte Carlo simulation: the
#' planned \emph{n} is obtained from \code{\link{ss_aipe_icc}} with
#' \code{estimated_rho}, then samples are drawn from the \emph{true}
#' population (with population ICC \code{true_rho}) and the realized
#' confidence interval widths are summarized.
#'
#' \strong{Data generating model.} For \code{type} starting with
#' \code{"ICC(1,"} the simulation uses the one-way random model: each row
#' (subject) gets a subject random effect, every cell adds independent
#' Gaussian noise, and the population ICC equals
#' \eqn{\sigma^2_{\mathrm{subj}} / (\sigma^2_{\mathrm{subj}} +
#' \sigma^2_{\mathrm{err}})}. For \code{type} starting with \code{"ICC(2,"}
#' the simulation also adds a rater random effect, so the population
#' \code{ICC(2,1)} equals \eqn{\sigma^2_{\mathrm{subj}} /
#' (\sigma^2_{\mathrm{subj}} + \sigma^2_{\mathrm{rater}} +
#' \sigma^2_{\mathrm{err}})}. For \code{type} starting with \code{"ICC(3,"}
#' the rater effect is fixed (centered constants), so the population
#' \code{ICC(3,1)} equals \eqn{\sigma^2_{\mathrm{subj}} /
#' (\sigma^2_{\mathrm{subj}} + \sigma^2_{\mathrm{err}})} but the
#' two-way decomposition is used in the estimator. Within each cell, the
#' total variance is unity by construction. For the single-rater forms
#' \code{true_rho} therefore maps directly to the subject-variance share;
#' for the average-of-\eqn{k} forms \code{true_rho} is first mapped to the
#' single-rater scale through the inverse Spearman-Brown relation
#' \eqn{\rho = \rho_k / [k - (k - 1)\rho_k]}, so that the population ICC
#' at the average-of-\eqn{k} level equals \code{true_rho}. Coverage is
#' checked against \code{true_rho} on its own scale, matching the scale on
#' which \code{\link{icc}} reports each form.
#'
#' @return
#' A \code{data.frame} with columns \code{term} and \code{value}
#' summarizing the Monte Carlo results across the \code{G} replications.
#' The \code{term} entries are: \code{"mean_icc"}, \code{"median_icc"},
#' \code{"sd_icc"} (mean / median / SD of the \emph{G} observed ICC
#' estimates); \code{"mean_ci_width"}, \code{"median_ci_width"},
#' \code{"sd_ci_width"} (corresponding summaries of the realized interval
#' widths); \code{"pct_ci_less_w"} (proportion of intervals at or below
#' the planning width \code{width}); \code{"pct_ci_miss_low"} and
#' \code{"pct_ci_miss_high"} (tail-specific non-coverage of
#' \code{true_rho}); \code{"total_type_I_error"} (overall empirical
#' non-coverage of \code{true_rho}); plus the input echoes
#' \code{"total_N"}, \code{"k"}, \code{"true_rho"},
#' \code{"estimated_rho"}, \code{"width"}, \code{"conf_level"}, and
#' \code{"assurance"} (present only when an assurance was supplied).
#' The ICC type is not a row; it is stored as the \code{"icc_type"}
#' attribute on the returned object so the \code{value} column stays
#' numeric.
#'
#' @references
#' Bonett, D. G. (2002). Sample size requirements for estimating intraclass
#'   correlations with desired precision. \emph{Statistics in Medicine, 21}(9),
#'   1331--1335. \doi{10.1002/sim.1108}
#'
#' Kelley, K. (2007). Confidence intervals for standardized effect sizes:
#'   Theory, application, and implementation. \emph{Journal of Statistical
#'   Software, 20}(8), 1--24. \doi{10.18637/jss.v020.i08}
#'
#' Kelley, K., & Maxwell, S. E. (2003). Sample size for multiple
#'   regression: Obtaining regression coefficients that are accurate,
#'   not simply significant. \emph{Psychological Methods, 8}(3),
#'   305--321. \doi{10.1037/1082-989X.8.3.305}
#'
#' Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). \emph{Designing
#'   experiments and analyzing data: A model comparison perspective}
#'   (4th ed.). Routledge. (See Chapters 10 and 11 on intraclass
#'   correlation and reliability.)
#'
#' Shrout, P. E., & Fleiss, J. L. (1979). Intraclass correlations: Uses in
#'   assessing rater reliability. \emph{Psychological Bulletin, 86}(2),
#'   420--428.
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso \code{\link{ss_aipe_icc}}, \code{\link{icc}}, \code{\link{var_icc}}
#'
#' @examples
#' # Reduced G and a wide target width keep this Monte Carlo example fast;
#' # raise G (e.g., 1000 or more) for stable empirical-coverage estimates.
#'
#' # Well-specified case: plan with the true population ICC.
#' set.seed(113)
#' ss_aipe_icc_sensitivity(
#'   true_rho      = 0.70,
#'   estimated_rho = 0.70,
#'   k             = 3,
#'   width         = 0.40,
#'   conf_level    = 0.95,
#'   G             = 25,
#'   print_iter    = FALSE
#' )
#'
#' # Misspecified case: the planner used .70 but the truth is .50.
#' # The realized interval widths will tend to be wider than the target.
#' set.seed(113)
#' ss_aipe_icc_sensitivity(
#'   true_rho      = 0.50,
#'   estimated_rho = 0.70,
#'   k             = 3,
#'   width         = 0.40,
#'   conf_level    = 0.95,
#'   G             = 25,
#'   print_iter    = FALSE
#' )
#'
#' # Fixed-n mode: skip the planner and evaluate at a chosen sample size.
#' set.seed(113)
#' ss_aipe_icc_sensitivity(
#'   true_rho    = 0.50,
#'   specified_N = 40,
#'   k           = 3,
#'   width       = 0.40,
#'   G           = 25,
#'   print_iter  = FALSE
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
ss_aipe_icc_sensitivity <- function(true_rho = NULL, estimated_rho = NULL, k, width,
                                    assurance = NULL, specified_N = NULL,
                                    conf_level = 0.95,
                                    type = c("ICC(1,1)", "ICC(2,1)", "ICC(3,1)",
                                             "ICC(1,k)", "ICC(2,k)", "ICC(3,k)"),
                                    G = 1000, print_iter = FALSE,
                                    save = FALSE,
                                    filename = "ss_aipe_icc_sensitivity_result.csv") {
  type <- match.arg(type)
  if (is.null(estimated_rho) && is.null(specified_N)) {
    stop("You must specify either 'estimated_rho' or 'specified_N'.", call. = FALSE)
  }
  if (!is.null(estimated_rho) && !is.null(specified_N)) {
    stop("You must specify 'estimated_rho' or 'specified_N', but not both.", call. = FALSE)
  }
  if (is.null(true_rho) || !is.numeric(true_rho) || true_rho < 0 || true_rho >= 1) {
    stop("'true_rho' must be a single value in [0, 1).", call. = FALSE)
  }
  if (!is.numeric(k) || length(k) != 1L || k < 2) {
    stop("'k' must be a single integer >= 2.", call. = FALSE)
  }
  if (!is.numeric(width) || length(width) != 1L || width <= 0) {
    stop("'width' must be a single positive number.", call. = FALSE)
  }
  if (!is.numeric(G) || length(G) != 1L || G < 1) {
    stop("'G' must be a positive integer.", call. = FALSE)
  }

  # Resolve sample size.
  if (!is.null(estimated_rho)) {
    if (!is.numeric(estimated_rho) || estimated_rho < 0 || estimated_rho >= 1) {
      stop("'estimated_rho' must be a single value in [0, 1).", call. = FALSE)
    }
    plan <- ss_aipe_icc(rho = estimated_rho, k = k, width = width,
                        conf_level = conf_level, type = type,
                        assurance = assurance)
    n <- plan$value[plan$term == "necessary_N"]
  } else {
    if (!is.numeric(specified_N) || length(specified_N) != 1L || specified_N < 2) {
      stop("'specified_N' must be a single integer >= 2.", call. = FALSE)
    }
    n <- specified_N
  }
  n <- as.integer(n)

  # Choose the data generating model from the planning type.
  type_canon <- type
  family <- if (grepl("^ICC\\(1,", type_canon)) "oneway"
            else if (grepl("^ICC\\(2,", type_canon)) "tworandom"
            else "twomixed"

  # 'true_rho' arrives on the scale matching 'type'. The variance components
  # live on the single-rater scale, so an average-of-k value is first mapped
  # down through the inverse Spearman-Brown relation; the population ICC at
  # the requested level then equals true_rho (for every family the
  # Spearman-Brown map carries the single-rater share to the average-of-k
  # ICC). Coverage below is still checked against true_rho itself, because
  # icc() reports the average-of-k forms on the average-of-k scale.
  average_k <- type_canon %in% c("ICC(1,k)", "ICC(2,k)", "ICC(3,k)")
  rho_single <- if (average_k) true_rho / (k - (k - 1) * true_rho)
                else true_rho

  # Variance components such that the population ICC at the requested level
  # equals true_rho and total variance equals 1.
  sigma2_subj <- rho_single
  sigma2_err  <- 1 - rho_single
  sigma2_rater <- 0
  if (family == "tworandom") {
    # Distribute the non-subject share equally between rater and residual
    # so total variance stays at 1 and ICC(2,1) = sigma2_subj.
    sigma2_rater <- (1 - rho_single) / 2
    sigma2_err   <- (1 - rho_single) / 2
  }

  # Storage.
  icc_hat       <- numeric(G)
  ci_lo         <- numeric(G)
  ci_hi         <- numeric(G)
  ci_width      <- numeric(G)
  type_I_lower  <- logical(G)
  type_I_upper  <- logical(G)

  # Rater fixed effects (only used by family == "twomixed").
  if (family == "twomixed") {
    rater_fx <- seq(-1, 1, length.out = k) * sqrt(3 * (1 - rho_single) / 2)
    rater_fx <- rater_fx - mean(rater_fx)
  } else {
    rater_fx <- rep(0, k)
  }

  for (g in seq_len(G)) {
    if (isTRUE(print_iter)) cat(g, "\n")

    # Subject random effect, replicated across the k columns.
    alpha_i <- stats::rnorm(n, mean = 0, sd = sqrt(sigma2_subj))
    M <- matrix(alpha_i, nrow = n, ncol = k, byrow = FALSE)

    # Optional rater random effect, replicated across the n rows.
    if (family == "tworandom") {
      beta_j <- stats::rnorm(k, mean = 0, sd = sqrt(sigma2_rater))
      M <- M + matrix(beta_j, nrow = n, ncol = k, byrow = TRUE)
    } else if (family == "twomixed") {
      M <- M + matrix(rater_fx, nrow = n, ncol = k, byrow = TRUE)
    }

    # Residual.
    M <- M + matrix(stats::rnorm(n * k, mean = 0, sd = sqrt(sigma2_err)),
                    nrow = n, ncol = k)

    fit <- icc(M, type = type_canon, conf_level = conf_level)
    icc_hat[g]      <- fit$value
    ci_lo[g]        <- fit$lower_limit
    ci_hi[g]        <- fit$upper_limit
    ci_width[g]     <- fit$upper_limit - fit$lower_limit
    type_I_lower[g] <- true_rho < fit$lower_limit
    type_I_upper[g] <- true_rho > fit$upper_limit
  }

  if (isTRUE(save)) {
    out_per_rep <- data.frame(
      icc_hat = icc_hat, ci_lower = ci_lo, ci_upper = ci_hi,
      ci_width = ci_width, type_I_lower = type_I_lower,
      type_I_upper = type_I_upper
    )
    suppressWarnings(file_exist <- try(utils::read.csv(filename), silent = TRUE))
    if (!is.null(dim(file_exist))) {
      utils::write.table(out_per_rep, filename, sep = ",", row.names = FALSE,
                         col.names = FALSE, append = TRUE)
    } else {
      utils::write.table(out_per_rep, filename, sep = ",", row.names = FALSE,
                         append = FALSE)
    }
  }

  out <- data.frame(
    term  = c("mean_icc", "median_icc", "sd_icc",
              "mean_ci_width", "median_ci_width", "sd_ci_width",
              "pct_ci_less_w", "pct_ci_miss_low", "pct_ci_miss_high",
              "total_type_I_error",
              "total_N", "k", "true_rho",
              "estimated_rho", "width", "conf_level",
              if (!is.null(assurance)) "assurance"),
    value = c(mean(icc_hat, na.rm = TRUE),
              stats::median(icc_hat, na.rm = TRUE),
              stats::sd(icc_hat, na.rm = TRUE),
              mean(ci_width, na.rm = TRUE),
              stats::median(ci_width, na.rm = TRUE),
              stats::sd(ci_width, na.rm = TRUE),
              mean(ci_width <= width, na.rm = TRUE),
              mean(type_I_lower, na.rm = TRUE),
              mean(type_I_upper, na.rm = TRUE),
              mean(type_I_lower | type_I_upper, na.rm = TRUE),
              n, k, true_rho,
              if (is.null(estimated_rho)) NA_real_ else estimated_rho,
              width, conf_level,
              if (!is.null(assurance)) assurance),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  attr(out, "icc_type") <- type_canon
  .as_dmar_tbl(out, conf_level = conf_level)
}
