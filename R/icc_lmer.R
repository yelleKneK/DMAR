# Tidy ICC + CI from a fitted lme4::lmer object.
#' Intraclass Correlation From a Fitted \code{lme4} Mixed-Effects Model
#'
#' Reads the variance-component decomposition off a fitted
#' \code{\link[lme4]{lmer}} model and returns the implied intraclass
#' correlation (ICC) for the specified grouping factor along with a
#' Bonett (2002) Fisher-\eqn{L}-transform confidence interval, in
#' tidy long form. The bridge function between DMAR's classical-
#' ANOVA path (\code{\link{variance_components_mls}},
#' \code{\link{var_icc}}) and the modern mixed-effects path.
#'
#' @param fit A fitted \code{\link[lme4]{lmer}} model (class
#'   \code{lmerMod}) with at least one random-effects grouping factor.
#' @param group Character name of the grouping factor whose ICC is
#'   wanted. If \code{NULL} (default), the first random-effects
#'   grouping factor is used. For three-level models, supply the
#'   target level explicitly.
#' @param conf_level Confidence level. Default \code{0.95}.
#'
#' @return A \code{data.frame} with rows for the point estimate
#'   of the ICC, the variance components (between-group and residual),
#'   the implied total variance, the Bonett (2002) CI lower/upper
#'   limits, and the cluster-level effective sample size used in the
#'   CI.
#'
#' @details
#' \strong{Definition.} For a two-level model with random intercept by
#' group \eqn{j} and within-group residual,
#' \deqn{\rho \;=\; \frac{\sigma^2_b}{\sigma^2_b + \sigma^2_w},}
#' read directly from \code{VarCorr(fit)}. For three-level models, the
#' user specifies which level (\code{group}) provides the variance
#' contribution; the denominator is the total variance summed across
#' all variance components.
#'
#' \strong{Bonett (2002) CI.} The CI is built on the Fisher-style
#' \eqn{L}-transformation
#' \deqn{L \;=\; \tfrac{1}{2} \log\left(\frac{1 + (k - 1) \rho}{1 - \rho}\right),}
#' with variance \eqn{k / (2 (k - 1) (n - 2))}, where \eqn{k} is the
#' average cluster size and \eqn{n} is the number of clusters. Back-
#' transformation keeps the bounds in \eqn{[0, 1]}.
#'
#' \strong{Limitations.} The Bonett CI is built on a balanced /
#' approximately balanced approximation; for severely unbalanced
#' designs the profile likelihood CI from \code{confint(fit, ...)} is
#' preferable.
#'
#' @references
#' Bonett, D. G. (2002). Sample size requirements for estimating
#'   intraclass correlations with desired precision. \emph{Statistics
#'   in Medicine, 21}(9), 1331--1335. \doi{10.1002/sim.1108}
#'
#' Donner, A. (1986). A review of inference procedures for the
#'   intraclass correlation coefficient in the one-way random effects
#'   model. \emph{International Statistical Review, 54}(1), 67--82.
#'
#' Snijders, T. A. B., & Bosker, R. J. (2012). \emph{Multilevel
#'   analysis: An introduction to basic and advanced multilevel
#'   modeling} (2nd ed.). Sage.
#'
#' @seealso \code{\link{icc}}, \code{\link{var_icc}},
#'   \code{\link{variance_components_mls}}, \code{\link{ss_aipe_icc}},
#'   \code{\link[lme4]{lmer}}
#'
#' @examples
#' # Twenty groups of six, with a group-level standard deviation of 0.7 on
#' # top of within-group noise with a standard deviation of 1, so the data
#' # generating ICC is 0.7^2 / (0.7^2 + 1) = 0.329. The estimate below comes
#' # with a wide interval: 20 clusters is not many, and the interval is what
#' # keeps that fact visible.
#' set.seed(113)
#' n_grp <- 20; n_per <- 6
#' grp <- factor(rep(1:n_grp, each = n_per))
#' y   <- rnorm(n_grp * n_per, 0, 1) + rep(rnorm(n_grp, 0, 0.7), each = n_per)
#' d <- data.frame(y, grp)
#' fit <- lme4::lmer(y ~ 1 + (1 | grp), data = d)
#' icc_lmer(fit)
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @keywords htest
#'
#' @family agreement and measurement
#' @family mixed models
#'
#' @export

icc_lmer <- function(fit, group = NULL, conf_level = 0.95) {
  if (!requireNamespace("lme4", quietly = TRUE))
    stop("The 'lme4' package is required. Install it with: ",
         "install.packages('lme4')")
  if (!inherits(fit, "lmerMod"))
    stop("'fit' must be a fitted lme4::lmer model.")
  if (conf_level <= 0 || conf_level >= 1)
    stop("'conf_level' must be in (0, 1).")

  vc <- as.data.frame(lme4::VarCorr(fit))
  random_vc <- vc[vc$grp != "Residual" & is.na(vc$var2) & vc$var1 == "(Intercept)", , drop = FALSE]
  if (nrow(random_vc) == 0L)
    stop("No random-intercept variance components found in 'fit'.")

  resid_var <- vc$vcov[vc$grp == "Residual"]
  if (length(resid_var) == 0L)
    stop("Could not find the residual variance in VarCorr(fit).")

  if (is.null(group)) {
    target <- random_vc$grp[1L]
  } else {
    if (!group %in% random_vc$grp)
      stop(sprintf("'group' = '%s' is not a random-intercept grouping in 'fit'.",
                   group))
    target <- group
  }
  sigma2_b <- random_vc$vcov[random_vc$grp == target]
  total_var <- sum(random_vc$vcov) + resid_var
  rho <- sigma2_b / total_var

  # Approximate effective n and k from lme4's grouping-factor list, which is
  # keyed by the same names VarCorr() reports and so also covers interaction
  # groupings such as (1 | A:B), where model.frame(fit) has no literal
  # "A:B" column. Fall back to model.frame() only if the name is somehow
  # absent from flist.
  flist <- lme4::getME(fit, "flist")
  if (target %in% names(flist)) {
    g_vec <- flist[[target]]
  } else {
    mf <- stats::model.frame(fit)
    if (!target %in% names(mf))
      stop(sprintf("Grouping factor '%s' not found in the fitted model.",
                   target))
    g_vec <- mf[[target]]
  }
  cl_tbl <- table(g_vec)
  n_cl   <- length(cl_tbl)
  k_avg  <- mean(cl_tbl)

  # Bonett (2002) CI:
  if (n_cl < 4L) {
    lo <- hi <- NA_real_
  } else {
    L <- 0.5 * log((1 + (k_avg - 1) * rho) / (1 - rho))
    # Bonett (2002) variance of L, var(L) = k / (2 (k - 1) (n - 2)). The k/(k-1)
    # factor was previously omitted, which made the SE too small by sqrt(2) at
    # k = 2 and produced anti-conservative intervals.
    se_L <- sqrt(k_avg / (2 * (k_avg - 1) * (n_cl - 2)))
    z <- stats::qnorm(1 - (1 - conf_level) / 2)
    L_lo <- L - z * se_L; L_hi <- L + z * se_L
    # Back-transform from L = 0.5 log((1 + (k-1) rho) / (1 - rho))
    inv_L <- function(LL) {
      ee <- exp(2 * LL)
      (ee - 1) / (ee + k_avg - 1)
    }
    lo <- max(0, inv_L(L_lo)); hi <- min(1, inv_L(L_hi))
  }

  out <- data.frame(
    term  = c("icc", "sigma2_between", "sigma2_within", "total_variance",
              "lower_limit", "upper_limit", "n_clusters", "average_cluster_size"),
    value = c(rho, sigma2_b, resid_var, total_var, lo, hi, n_cl, k_avg),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  .as_dmar_tbl(out, conf_level = conf_level)
}
