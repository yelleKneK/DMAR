# Dunnett-adjusted CIs for comparisons against a control group.
#' Dunnett's Simultaneous Confidence Intervals Against a Control
#'
#' Computes Dunnett's (1955, 1964) simultaneous confidence intervals
#' for the \eqn{a - 1} comparisons of \eqn{a - 1} treatment means
#' against a single control mean, with family-wise coverage at the
#' specified \code{conf_level}. Returns the result in tidy long form.
#'
#' @param x Either (a) a fitted \code{\link[stats]{lm}} or
#'   \code{\link[stats]{aov}} object with a one-way factor predictor,
#'   or (b) a numeric vector of observations, in which case
#'   \code{group} must also be supplied.
#' @param group When \code{x} is a vector, a factor of group labels
#'   of the same length.
#' @param control Character name of the control level (must be one of
#'   the factor levels). If \code{NULL} (default), the first level is
#'   used.
#' @param alternative One of \code{"two_sided"} (default; the base-R
#'   spelling \code{"two.sided"} is accepted as an alias),
#'   \code{"less"}, or \code{"greater"}.
#' @param conf_level Family-wise confidence level. Default
#'   \code{0.95}.
#'
#' @return A \code{data.frame} with one row per non-control
#'   level. Columns: \code{contrast}, \code{mean_difference},
#'   \code{se}, \code{t_statistic}, \code{lower_limit},
#'   \code{upper_limit}, \code{p_adjusted}.
#'
#' @details
#' \strong{Critical value.} The two-sided Dunnett critical value
#' \eqn{d_{\alpha, a - 1, \nu}^{(2)}} is obtained from the multivariate
#' \emph{t} distribution with \eqn{a - 1} dimensions, common correlation
#' \eqn{0.5} (the Dunnett correlation under balanced \emph{n}; the
#' function does not adjust for unequal \emph{n}), and \eqn{\nu} error
#' degrees of freedom. The function uses the existing
#' \code{\link{cv_dunnett}()} critical value.
#'
#' \strong{Adjusted \emph{p}-values.} Computed exactly from the same
#' equicorrelated multivariate \emph{t} distribution. The one common
#' correlation \eqn{1/2} admits a one-factor representation, so the
#' probability that all comparisons fall inside (or below) the observed
#' statistic collapses to two nested one-dimensional integrals, evaluated
#' by quadrature. The adjusted \emph{p}-value is one minus that
#' probability. The computation is deterministic (no Monte Carlo) and
#' needs no additional package.
#'
#' @references
#' Dunnett, C. W. (1955). A multiple comparison procedure for
#'   comparing several treatments with a control. \emph{Journal of the
#'   American Statistical Association, 50}(272), 1096--1121.
#'
#' Dunnett, C. W. (1964). New tables for multiple comparisons with a
#'   control. \emph{Biometrics, 20}(3), 482--491.
#'
#' Hsu, J. C. (1996). \emph{Multiple comparisons: Theory and methods}.
#'   Chapman & Hall.
#'
#' @seealso \code{\link{cv_dunnett}}, \code{\link{ci_tukey_kramer}},
#'   \code{\link{ci_scheffe}}
#'
#' @examples
#' # 1. Compare the SSRI and placebo arms of the depression_bdi study
#' #    against the wait list control:
#' fit <- lm(bdi_post ~ condition, data = depression_bdi)
#' ci_dunnett(fit, control = "wait_list")
#'
#' # 2. One-sided: a treatment that works pulls the posttest BDI down,
#' #    so the directional alternative is "less":
#' ci_dunnett(fit, control = "wait_list", alternative = "less")
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @keywords htest
#'
#' @family hypothesis tests
#'
#' @export

ci_dunnett <- function(x, group = NULL, control = NULL,
                       alternative = c("two_sided", "less", "greater"),
                       conf_level = 0.95) {
  alternative <- .match_alternative(alternative)
  if (conf_level <= 0 || conf_level >= 1)
    stop("'conf_level' must be in (0, 1).")

  if (inherits(x, c("lm", "aov"))) {
    mf <- stats::model.frame(x)
    if (ncol(mf) != 2L)
      stop("ci_dunnett() requires a fitted lm/aov with one factor predictor.")
    y <- mf[[1L]]; g <- mf[[2L]]
    if (!is.factor(g)) g <- factor(g)
    tab <- stats::anova(x); rn <- trimws(rownames(tab))
    ms_e <- tab[rn == "Residuals", "Mean Sq"]
    df_e <- tab[rn == "Residuals", "Df"]
  } else if (is.numeric(x) && !is.null(group)) {
    if (length(x) != length(group))
      stop("'x' and 'group' must be the same length.")
    g <- factor(group); y <- x
    ok <- !is.na(y) & !is.na(g)
    y <- y[ok]; g <- g[ok]
    fit <- stats::aov(y ~ g)
    tab <- stats::anova(fit); rn <- trimws(rownames(tab))
    ms_e <- tab[rn == "Residuals", "Mean Sq"]
    df_e <- tab[rn == "Residuals", "Df"]
  } else {
    stop("Supply either a fitted lm/aov, or (x, group) vectors.")
  }

  lev <- levels(g)
  if (is.null(control)) control <- lev[1L]
  if (!control %in% lev)
    stop(sprintf("'control' = '%s' is not one of the factor levels.", control))
  others <- setdiff(lev, control)
  a <- length(lev); k <- a - 1L
  means <- tapply(y, g, mean)
  ns    <- tapply(y, g, length)
  m_ctrl <- as.numeric(means[control]); n_ctrl <- as.numeric(ns[control])

  alpha_fw <- 1 - conf_level
  alt_for_cv <- switch(alternative,
    two_sided = "not_equal",
    less      = "less",
    greater   = "greater"
  )
  d_crit <- as.numeric(cv_dunnett(alpha_level = alpha_fw, df = df_e,
                                  n_comparisons = k,
                                  alternative = alt_for_cv,
                                  verbose = FALSE)$value[1])
  # cv_dunnett() reports the "less" critical value on the signed scale of its
  # rejection region (a negative lower-tail quantile). The interval arithmetic
  # below is written as diff_i - d_crit * se and diff_i + d_crit * se with a
  # positive multiplier for every alternative, so keep only the magnitude here;
  # the alternative determines which side is finite, not the sign of d_crit.
  d_crit <- abs(d_crit)

  rows <- vector("list", k)
  for (i in seq_along(others)) {
    g_i <- others[i]
    m_i <- as.numeric(means[g_i]); n_i <- as.numeric(ns[g_i])
    diff_i <- m_i - m_ctrl
    se     <- sqrt(ms_e * (1 / n_i + 1 / n_ctrl))
    t_obs  <- diff_i / se
    # Exact adjusted p-value from the equicorrelated (rho = 1/2) multivariate t:
    # 1 - P(all comparisons inside/below the observed statistic). For "less",
    # P(all T_i > t_obs) = P(all -T_i <= -t_obs) by symmetry of the reference.
    p_adj <- switch(alternative,
      two_sided = 1 - .dunnett_cdf(abs(t_obs), df_e, k, two_sided = TRUE),
      greater   = 1 - .dunnett_cdf(t_obs,      df_e, k, two_sided = FALSE),
      less      = 1 - .dunnett_cdf(-t_obs,     df_e, k, two_sided = FALSE)
    )
    p_adj <- min(1, max(0, p_adj))
    lo <- switch(alternative,
      two_sided = diff_i - d_crit * se,
      less      = -Inf,
      greater   = diff_i - d_crit * se)
    hi <- switch(alternative,
      two_sided = diff_i + d_crit * se,
      less      = diff_i + d_crit * se,
      greater   = Inf)
    rows[[i]] <- data.frame(
      contrast        = paste(g_i, "-", control),
      mean_difference = diff_i,
      se              = se,
      t_statistic     = t_obs,
      lower_limit     = lo,
      upper_limit     = hi,
      p_adjusted      = p_adj,
      stringsAsFactors = FALSE)
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out <- .as_dmar_tbl(out, conf_level = conf_level)
  class(out) <- c("dmar_post_hoc_ci", class(out))
  out
}
