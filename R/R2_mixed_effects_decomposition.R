# Rights and Sterba (2019) integrative framework of R-squared measures for
# mixed-effects (multilevel) models. Native reimplementation of the method whose
# reference implementation is the r2mlm package (Rights & Sterba); validated to
# match it.

#' R-Squared Measures for Mixed-Effects Models
#'
#' Computes the Rights and Sterba (2019) integrative framework of R-squared
#' measures for a fitted two-level linear mixed-effects (multilevel) model. The
#' model implied outcome variance is fully decomposed into five sources:
#' variance due to level-1 predictors via fixed slopes (\eqn{f_1}), level-2
#' predictors via fixed slopes (\eqn{f_2}), predictors via random slope
#' (co)variation (\eqn{v}), cluster-specific outcome means via random intercept
#' variation (\eqn{m}), and level-1 residuals (\eqn{\sigma^2}). Proportions of
#' the total, within-cluster, and between-cluster outcome variance attributable
#' to combinations of these sources give the family of R-squared measures.
#'
#' The companion \code{\link{R2_mixed_effects}} returns the Nakagawa and
#' Schielzeth marginal and conditional R-squared; this function contains those
#' two as special cases (\code{total_f} and \code{total_fvm}) within the fuller
#' source decomposition.
#'
#' @param model A fitted two-level model of class \code{merMod} (from
#'   \pkg{lme4}, e.g.\ \code{\link[lme4]{lmer}}) or \code{lme} (from
#'   \pkg{nlme}). The model must use numeric predictors (only the cluster
#'   variable may be a factor) and must not contain \code{I()} terms; create
#'   any transformed predictors as their own columns first.
#'
#' @return A \code{data.frame} (\code{dmar_tbl}) with columns \code{term} and
#'   \code{value}. When the level-1 predictors are cluster-mean-centered, the
#'   full set of 12 measures is returned, named \code{total_f1}, \code{total_f2},
#'   \code{total_v}, \code{total_m}, \code{total_f}, \code{total_fv},
#'   \code{total_fvm}, \code{within_f1}, \code{within_v}, \code{within_fv},
#'   \code{between_f2}, and \code{between_m}; otherwise the five total-variance
#'   measures \code{total_f}, \code{total_v}, \code{total_m}, \code{total_fv},
#'   and \code{total_fvm} are returned (the within/between split requires
#'   cluster-mean-centering). The returned object carries the source-by-target
#'   variance decomposition in \code{attr(x, "decomposition")}.
#'
#' @details
#' The measure superscripts index the variance sources in the numerator and the
#' subscripts index the outcome variance in the denominator: \code{total_*} use
#' the total outcome variance, \code{within_*} the within-cluster variance
#' (\eqn{f_1 + v + \sigma^2}), and \code{between_*} the between-cluster variance
#' (\eqn{f_2 + m}). \code{total_fvm} is the omnibus measure (all explained
#' sources over the total variance) and, for a random-intercept model, coincides
#' with the Nakagawa and Schielzeth conditional R-squared computed by
#' \code{\link{R2_mixed_effects}}; \code{total_f} coincides with their marginal
#' R-squared. See Rights and Sterba (2019, Table 1) for the definitions.
#'
#' The measures were derived under the assumption that the fitted model uses
#' cluster-mean-centering of the level-1 predictors (with the cluster means
#' entered as level-2 predictors). When that centering is not detected, only the
#' total-variance measures are returned, matching the reference implementation.
#'
#' Fitting the model requires \pkg{lme4} (for a \code{merMod} fit) or \pkg{nlme}
#' (for an \code{lme} fit) to be installed.
#'
#' @references
#' Rights, J. D., & Sterba, S. K. (2019). Quantifying explained variance in
#'   multilevel models: An integrative framework for defining R-squared
#'   measures. \emph{Psychological Methods, 24}(3), 309--338.
#'   \doi{10.1037/met0000184}
#'
#' Nakagawa, S., & Schielzeth, H. (2013). A general and simple method for
#'   obtaining \eqn{R^2} from generalized linear mixed-effects models.
#'   \emph{Methods in Ecology and Evolution, 4}(2), 133--142.
#'   \doi{10.1111/j.2041-210x.2012.00261.x}
#'
#' @seealso \code{\link{R2_mixed_effects}}, \code{\link{icc_lmer}}
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @examples
#' fit <- lme4::lmer(Reaction ~ Days + (Days | Subject),
#'                   data = lme4::sleepstudy)
#' R2_mixed_effects_decomposition(fit)
#'
#' @concept mixed-effects R-squared
#' @concept multilevel R-squared
#'
#' @keywords models
#'
#' @family mixed-effects R-squared
#' @family mixed models
#'
#' @export
R2_mixed_effects_decomposition <- function(model) {
  if (inherits(model, "merMod")) {
    parts <- .R2_mixed_effects_decomposition_extract_lme4(model)
  } else if (inherits(model, "lme")) {
    parts <- .R2_mixed_effects_decomposition_extract_nlme(model)
  } else {
    stop("'model' must be a fitted lme4 (merMod) or nlme (lme) two-level model.",
         call. = FALSE)
  }
  comp <- do.call(.R2_mixed_effects_decomposition_compute, parts)

  if (parts$cluster_mean_centered) {
    term <- c("total_f1", "total_f2", "total_v", "total_m", "total_f",
              "total_fv", "total_fvm", "within_f1", "within_v", "within_fv",
              "between_f2", "between_m")
    value <- with(comp, c(R2_t_f1, R2_t_f2, R2_t_v, R2_t_m, R2_t_f, R2_t_fv,
                          R2_t_fvm, R2_w_f1, R2_w_v, R2_w_fv, R2_b_f2, R2_b_m))
  } else {
    term  <- c("total_f", "total_v", "total_m", "total_fv", "total_fvm")
    value <- with(comp, c(R2_t_f, R2_t_v, R2_t_m, R2_t_fv, R2_t_fvm))
  }

  out <- .as_dmar_tbl(data.frame(term = term, value = value,
                                 stringsAsFactors = FALSE))
  attr(out, "decomposition") <- comp$decomposition
  class(out) <- c("dmar_R2_mixed_effects_decomposition", class(out))
  out
}


# Extract the Rights-Sterba quantities from a fitted lme4 model. Follows the
# reference r2mlm algorithm: split predictors into level-1 (varying within
# cluster) and level-2 (cluster constant), identify random-slope predictors,
# detect cluster-mean-centering, and read the fixed effects, random-effect
# covariance, and residual variance.
.R2_mixed_effects_decomposition_extract_lme4 <- function(model) {
  if (!requireNamespace("lme4", quietly = TRUE))
    stop("Package 'lme4' is required for R2_mixed_effects_decomposition() on a merMod fit.", call. = FALSE)
  tt <- stats::terms(model)
  if (any(grepl("I(", labels(tt), fixed = TRUE)))
    stop("R2_mixed_effects_decomposition() does not support I() terms; add transformed predictors as ",
         "their own columns.", call. = FALSE)
  has_intercept <- attr(tt, "intercept") == 1
  flist <- lme4::getME(model, "flist")
  if (length(flist) != 1L)
    stop("R2_mixed_effects_decomposition() supports two-level models with a single clustering factor.",
         call. = FALSE)
  cluster_variable <- names(flist)
  data <- as.data.frame(model.frame(model))
  fe <- lme4::fixef(model)
  predictors <- setdiff(names(fe), "(Intercept)")
  bad <- predictors[!vapply(predictors, function(v)
    is.numeric(data[[v]]) || is.integer(data[[v]]), logical(1))]
  if (length(bad))
    stop("All predictors must be numeric; offending: ",
         paste(bad, collapse = ", "), ". Only the cluster variable may be a ",
         "factor.", call. = FALSE)

  sorted <- .R2_mixed_effects_decomposition_sort_vars(data, predictors, cluster_variable)
  random_slope_vars <- setdiff(names(lme4::ranef(model)[[1]]), "(Intercept)")
  cmc <- if (length(sorted$l1) == 0L) TRUE else
    .R2_mixed_effects_decomposition_is_cwc(sorted$l1, cluster_variable, data)

  gamma_w <- unname(fe[sorted$l1])
  gamma_b <- if (has_intercept) unname(fe[c("(Intercept)", sorted$l2)])
             else unname(fe[sorted$l2])
  # Single clustering factor (enforced above), so every VarCorr component
  # belongs to it. Assemble the full random-effect covariance over
  # c("(Intercept)", random_slope_vars), the same column set and order as the
  # downstream var(cbind(1, random_slope_vars)). A single combined term such as
  # (x | g) supplies the whole block; separate terms such as (1 | g) +
  # (0 + x | g), or the (x || g) shorthand that lme4 expands to them, supply
  # disjoint blocks with zero cross-covariance. Reading a single VarCorr
  # component would use the first bar term only, mismatching the ranef columns
  # and producing a non-conformable variance decomposition.
  Tau <- .R2_mixed_effects_decomposition_assemble_tau(lme4::VarCorr(model), random_slope_vars)
  sigma2 <- lme4::getME(model, "sigma")^2

  list(data = data, l1 = sorted$l1, l2 = sorted$l2,
       random_slope_vars = random_slope_vars, gamma_w = gamma_w,
       gamma_b = gamma_b, Tau = Tau, sigma2 = sigma2,
       has_intercept = has_intercept, cluster_mean_centered = cmc)
}

.R2_mixed_effects_decomposition_extract_nlme <- function(model) {
  if (!requireNamespace("nlme", quietly = TRUE))
    stop("Package 'nlme' is required for R2_mixed_effects_decomposition() on an lme fit.", call. = FALSE)
  tt <- stats::terms(model)
  if (any(grepl("I(", labels(tt), fixed = TRUE)))
    stop("R2_mixed_effects_decomposition() does not support I() terms; add transformed predictors as ",
         "their own columns.", call. = FALSE)
  has_intercept <- attr(tt, "intercept") == 1
  grp <- names(model$groups)
  if (length(grp) != 1L)
    stop("R2_mixed_effects_decomposition() supports two-level models with a single clustering factor.",
         call. = FALSE)
  cluster_variable <- grp
  data <- as.data.frame(model$data)
  fe <- nlme::fixef(model)
  predictors <- setdiff(names(fe), "(Intercept)")
  sorted <- .R2_mixed_effects_decomposition_sort_vars(data, predictors, cluster_variable)
  random_slope_vars <- setdiff(colnames(nlme::ranef(model)), "(Intercept)")
  cmc <- if (length(sorted$l1) == 0L) TRUE else
    .R2_mixed_effects_decomposition_is_cwc(sorted$l1, cluster_variable, data)
  gamma_w <- unname(fe[sorted$l1])
  gamma_b <- if (has_intercept) unname(fe[c("(Intercept)", sorted$l2)])
             else unname(fe[sorted$l2])
  Tau <- as.matrix(nlme::getVarCov(model))
  sigma2 <- model$sigma^2
  list(data = data, l1 = sorted$l1, l2 = sorted$l2,
       random_slope_vars = random_slope_vars, gamma_w = gamma_w,
       gamma_b = gamma_b, Tau = Tau, sigma2 = sigma2,
       has_intercept = has_intercept, cluster_mean_centered = cmc)
}

# Assemble the random-effect covariance over c("(Intercept)",
# random_slope_vars) from the VarCorr component list of a single-grouping-factor
# model. Each component's named block is dropped into its matching positions;
# any position not covered by a component, such as the cross-covariance between
# two separate bar terms, stays zero. The intercept row and column are always
# present because the downstream var(cbind(1, random_slope_vars)) always carries
# a leading intercept column; a model with no random intercept simply leaves
# that entry at zero.
.R2_mixed_effects_decomposition_assemble_tau <- function(vc, random_slope_vars) {
  tau_names <- c("(Intercept)", random_slope_vars)
  Tau <- matrix(0, length(tau_names), length(tau_names),
                dimnames = list(tau_names, tau_names))
  for (comp in vc) {
    cn <- colnames(comp)
    Tau[cn, cn] <- as.matrix(comp)
  }
  Tau
}

# A predictor with zero total within-cluster variance is a level-2 (cluster
# constant) predictor; otherwise it is a level-1 predictor.
.R2_mixed_effects_decomposition_sort_vars <- function(data, predictors, cluster_variable) {
  l1 <- character(0); l2 <- character(0)
  for (v in predictors) {
    within_var <- sum(tapply(data[[v]], data[[cluster_variable]],
                             function(z) stats::var(z)), na.rm = TRUE)
    if (isTRUE(all.equal(within_var, 0))) l2 <- c(l2, v) else l1 <- c(l1, v)
  }
  list(l1 = l1, l2 = l2)
}

# Cluster-mean-centered means each cluster's sum (hence mean) of each level-1
# predictor is 0.
.R2_mixed_effects_decomposition_is_cwc <- function(l1_vars, cluster_variable, data) {
  for (v in l1_vars) {
    cs <- tapply(data[[v]], data[[cluster_variable]], sum)
    if (any(abs(cs) >= 1e-7)) return(FALSE)
  }
  TRUE
}

# The Rights-Sterba (2019) decomposition and R-squared measures from the
# extracted components (matrix algebra of Rights & Sterba, Equations 10-27).
.R2_mixed_effects_decomposition_compute <- function(data, l1, l2, random_slope_vars, gamma_w, gamma_b,
                            Tau, sigma2, has_intercept, cluster_mean_centered) {
  gw <- if (length(l1) == 0L) 0 else matrix(gamma_w, ncol = 1)
  gb <- if (length(gamma_b) == 0L) 0 else matrix(gamma_b, ncol = 1)

  phi_w <- if (length(l1) == 0L) 0 else stats::var(data[, l1, drop = FALSE], na.rm = TRUE)
  phi_b <- if (has_intercept)
             stats::var(cbind(1, data[, l2, drop = FALSE]), na.rm = TRUE)
           else if (length(l2) == 0L) 0
           else stats::var(data[, l2, drop = FALSE], na.rm = TRUE)

  var_rc <- stats::var(cbind(1, data[, random_slope_vars, drop = FALSE]), na.rm = TRUE)
  psi   <- matrix(diag(as.matrix(Tau)), ncol = 1)
  kappa <- if (length(Tau) == 1L) 0 else matrix(Tau[lower.tri(Tau)], ncol = 1)
  vv    <- matrix(diag(as.matrix(var_rc)), ncol = 1)
  rr    <- if (length(var_rc) == 1L) 0 else matrix(var_rc[lower.tri(var_rc)], ncol = 1)
  if (length(random_slope_vars) == 0L)
    mm <- matrix(1, ncol = 1)
  else
    mm <- matrix(colMeans(cbind(1, data[, random_slope_vars, drop = FALSE]),
                          na.rm = TRUE), ncol = 1)

  var_slopes <- as.numeric(t(vv) %*% psi + 2 * (t(rr) %*% kappa))
  f1 <- if (length(l1) == 0L) 0 else as.numeric(t(gw) %*% phi_w %*% gw)
  f2 <- if (identical(gb, 0)) 0 else as.numeric(t(gb) %*% phi_b %*% gb)
  m_total <- as.numeric(t(mm) %*% as.matrix(Tau) %*% mm)   # non-decomposed mean variance
  tau00   <- as.numeric(as.matrix(Tau)[1, 1])

  # Non-decomposed fixed variance (used when the model is NOT cluster-mean-
  # centered, where level-1 and level-2 predictors may covary): the full
  # quadratic form t(gamma) phi gamma over the joint predictor covariance,
  # rather than the sum f1 + f2 (which holds only under cluster-mean-centering).
  cols <- c(l1, l2)
  gamma_full <- if (has_intercept) c(1, gamma_w, if (length(l2)) gamma_b[-1] else NULL)
                else c(gamma_w, gamma_b)
  gamma_full <- matrix(gamma_full, ncol = 1)
  phi_full <- if (has_intercept)
                stats::var(cbind(1, data[, cols, drop = FALSE]), na.rm = TRUE)
              else stats::var(data[, cols, drop = FALSE], na.rm = TRUE)
  fixed_nd <- as.numeric(t(gamma_full) %*% phi_full %*% gamma_full)

  within_var  <- f1 + var_slopes + sigma2
  between_var <- f2 + tau00
  total_var   <- within_var + between_var
  total_var_nd <- fixed_nd + var_slopes + m_total + sigma2   # non-CMC total

  decomposition <- rbind(
    total   = c(fixed_within = f1, fixed_between = f2, slope_variation = var_slopes,
                mean_variation = if (cluster_mean_centered) tau00 else m_total,
                residual = sigma2),
    NULL
  )
  denom_total <- if (cluster_mean_centered) total_var else total_var_nd
  decomposition[1, ] <- decomposition[1, ] / denom_total

  out <- list(decomposition = decomposition)
  if (cluster_mean_centered) {
    out$R2_t_f1  <- f1 / total_var
    out$R2_t_f2  <- f2 / total_var
    out$R2_t_v   <- var_slopes / total_var
    out$R2_t_m   <- tau00 / total_var
    out$R2_t_f   <- (f1 + f2) / total_var
    out$R2_t_fv  <- (f1 + f2 + var_slopes) / total_var
    out$R2_t_fvm <- (f1 + f2 + var_slopes + tau00) / total_var
    out$R2_w_f1  <- f1 / within_var
    out$R2_w_v   <- var_slopes / within_var
    out$R2_w_fv  <- (f1 + var_slopes) / within_var
    out$R2_b_f2  <- f2 / between_var
    out$R2_b_m   <- tau00 / between_var
  } else {
    out$R2_t_f   <- fixed_nd / total_var_nd
    out$R2_t_v   <- var_slopes / total_var_nd
    out$R2_t_m   <- m_total / total_var_nd
    out$R2_t_fv  <- (fixed_nd + var_slopes) / total_var_nd
    out$R2_t_fvm <- (fixed_nd + var_slopes + m_total) / total_var_nd
  }
  out
}
