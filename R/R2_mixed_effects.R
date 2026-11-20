# Nakagawa & Schielzeth marginal and conditional R-squared for a fitted
# mixed-effects model, with the Johnson (2014) random-slope extension.
#' Marginal and Conditional \eqn{R^2} for a Mixed-Effects Model
#'
#' Computes the Nakagawa and Schielzeth (2013) marginal and conditional
#' coefficients of determination from a fitted mixed-effects model,
#' returned in tidy long form. The marginal \eqn{R^2} is the proportion
#' of total variance explained by the fixed effects alone; the
#' conditional \eqn{R^2} is the proportion explained by the fixed and
#' random effects together. The two quantities are the mixed-effects
#' companions to the intraclass correlation returned by
#' \code{\link{icc_lmer}}: where the ICC isolates the share of variance
#' attributable to a single grouping factor, the marginal and
#' conditional \eqn{R^2} summarize how much of the outcome variance the
#' fixed and random parts of the model account for.
#'
#' @param model A fitted mixed-effects model. A \code{\link[lme4]{lmer}}
#'   fit (class \code{lmerMod}) is the primary target; an
#'   \code{\link[nlme]{lme}} fit is also accepted.
#' @param conf_level Confidence level. Default \code{0.95}. Used only
#'   when a bootstrap interval is requested (see \code{ci_method}).
#' @param ci_method Interval method for the two \eqn{R^2} values. The
#'   default \code{"none"} returns point estimates only. \code{"boot"}
#'   adds a parametric bootstrap percentile interval (\code{lmerMod}
#'   models only, via \code{\link[lme4]{bootMer}}); see \emph{Details}.
#' @param B Number of bootstrap replications when \code{ci_method =
#'   "boot"}, a single whole number. Default \code{1000}, because every
#'   replication simulates a response and refits the model. The fewest
#'   accepted is 400 at the default \code{conf_level} (200 at 0.90, and
#'   2000 at 0.99, which the default does not meet); a smaller \code{B} is
#'   refused before anything is computed (class \code{dmar_count_invalid}).
#'   See \emph{Bootstrap replications} under \emph{Details} and
#'   \code{\link{dmar_bootstrap}}.
#' @param seed Optional integer seed for the bootstrap. Default
#'   \code{NULL}, under which the bootstrap draws from, and advances,
#'   the caller's current random number stream. When supplied, the seed
#'   is set inside the function and the caller's RNG state is restored
#'   on exit.
#'
#' @return A \code{data.frame} with rows \code{"R2_marginal"} and
#'   \code{"R2_conditional"} in the \code{value} column. When
#'   \code{ci_method = "boot"}, lower- and upper-limit rows for each
#'   quantity are appended (\code{R2_marginal_ci_lower},
#'   \code{R2_marginal_ci_upper}, \code{R2_conditional_ci_lower},
#'   \code{R2_conditional_ci_upper}), the confidence level is carried on the
#'   object, and the attributes \code{B}, \code{B_used},
#'   \code{n_dropped}, \code{B_required}, and \code{bootstrap_status}
#'   record the replications (see \emph{Bootstrap replications} under
#'   Details).
#'
#' @details
#' \strong{Variance decomposition.} Writing \eqn{\sigma^2_f} for the
#' variance of the fixed-effect linear predictor, \eqn{\sigma^2_r} for
#' the variance attributable to the random effects, and
#' \eqn{\sigma^2_\varepsilon} for the residual variance,
#' \deqn{R^2_{\mathrm{marginal}} \;=\;
#'   \frac{\sigma^2_f}{\sigma^2_f + \sigma^2_r + \sigma^2_\varepsilon},
#'   \qquad
#'   R^2_{\mathrm{conditional}} \;=\;
#'   \frac{\sigma^2_f + \sigma^2_r}{\sigma^2_f + \sigma^2_r + \sigma^2_\varepsilon}.}
#' The fixed-effect variance is
#' \eqn{\sigma^2_f = \mathrm{var}(\mathbf{X}\boldsymbol{\beta})}, the
#' variance of the fitted fixed-effect linear predictor across the
#' observations. The residual variance is
#' \eqn{\sigma^2_\varepsilon = \mathrm{sigma}(\mathrm{model})^2}.
#'
#' \strong{Random-effect variance.} For a random-intercept model the
#' random-effect variance is the sum of the variance components read off
#' \code{\link[lme4]{VarCorr}}. For a model with random slopes the
#' variance contributed by a random-effects term depends on the values
#' of the associated covariates, so the sum of the diagonal variance
#' components is not correct on its own. This function uses the Johnson
#' (2014) extension: for each random-effects term with design matrix
#' \eqn{\mathbf{Z}} and estimated covariance matrix
#' \eqn{\boldsymbol{\Sigma}}, its contribution is the mean over the
#' observations of the quadratic form
#' \eqn{\mathbf{z}_i^\top \boldsymbol{\Sigma}\, \mathbf{z}_i}, that is,
#' \eqn{\tfrac{1}{n}\,\mathrm{tr}(\mathbf{Z}\boldsymbol{\Sigma}\mathbf{Z}^\top)},
#' and \eqn{\sigma^2_r} is the sum of these contributions across all
#' random-effects terms. For a random-intercept term this reduces to the
#' intercept variance component, so the two paths agree.
#'
#' \strong{Scope.} The decomposition here is the one appropriate for a
#' Gaussian (identity-link) linear mixed model, which is what
#' \code{lmer} and \code{lme} fit. Generalized linear mixed models
#' introduce a distribution-specific variance term and are not handled
#' by this function.
#'
#' \strong{The bootstrap interval.} The default
#' \code{ci_method = "none"} reports the two point estimates alone, so
#' the bootstrap is what to ask for when the marginal and conditional
#' \eqn{R^2} are to be reported with an interval and the refits it
#' costs are affordable. With \code{ci_method = "boot"} the
#' interval comes from a parametric bootstrap
#' (\code{\link[lme4]{bootMer}}): each of the \code{B} replicates
#' (1000 by default) simulates a new response vector from the fitted
#' model, refits the model, and recomputes the two \eqn{R^2} values.
#' The unit of resampling is therefore a whole simulated data set drawn
#' from the estimated model, not a resampled set of cases. Only the
#' percentile interval is offered: the limits are the empirical
#' quantiles of the \code{B} bootstrap values (Efron & Tibshirani,
#' 1993); there is no BCa or bootstrap standard error variant.
#' Raising \code{B} tightens the Monte Carlo error of the reported
#' limits. Bootstrap results vary from run to run; supply \code{seed}
#' for reproducibility.
#'
#' \strong{Bootstrap replications.} The limits of each interval are
#' quantiles of the \code{B} replicated values, so the interval is run only
#' when \code{B} puts at least 10 replicates beyond each limit (400 at the
#' 95 percent level). A replication in which either \eqn{R^2} is undefined
#' (a refit that fails, or one whose variance components cannot be
#' computed) is dropped and counted, never replaced, and the interval is
#' computed from the rest. The result records \code{B} (requested),
#' \code{B_used} (the number the interval rests on), \code{n_dropped}, and
#' \code{B_required} (the fewest replications a request may ask for; that
#' count is held to \code{B}, not to \code{B_used}) as attributes, with
#' \code{bootstrap_status}
#' (\code{"complete"}, \code{"dropped"}, \code{"exploratory"},
#' \code{"insufficient_count"}, or \code{"insufficient_share"}). One
#' warning of class \code{dmar_bootstrap_replications_dropped} states the
#' counts when any replication was dropped, with a caution when the dropped
#' replications are at least the 2.5 percent of \code{B} that lies beyond
#' one limit at the 95 percent level, because they could then decide by
#' themselves where that limit falls. When more than half were dropped, the
#' four limit rows are \code{NA},
#' the estimate is returned, and the warning has class
#' \code{dmar_bootstrap_insufficient}. Such a replicate is dropped for both
#' quantities, so the two intervals rest on the same replicates.
#' \code{options(dmar.bootstrap_exploratory = TRUE)} lets a smaller
#' \code{B} (2 or more) run for a first look or a demonstration; every such
#' call gives a warning of class \code{dmar_bootstrap_exploratory} and the
#' result is not one to report. See \code{\link{dmar_bootstrap}} for the
#' rule and its reasons.
#'
#' @references
#' Efron, B., & Tibshirani, R. J. (1993). \emph{An introduction to the
#'   bootstrap}. New York, NY: Chapman & Hall/CRC.
#'
#' Johnson, P. C. D. (2014). Extension of Nakagawa & Schielzeth's
#'   \eqn{R^2_{GLMM}} to random slopes models. \emph{Methods in Ecology
#'   and Evolution, 5}(9), 944--946. \doi{10.1111/2041-210X.12225}
#'
#' Nakagawa, S., & Schielzeth, H. (2013). A general and simple method
#'   for obtaining \eqn{R^2} from generalized linear mixed-effects
#'   models. \emph{Methods in Ecology and Evolution, 4}(2), 133--142.
#'   \doi{10.1111/j.2041-210x.2012.00261.x}
#'
#' @seealso \code{\link{icc_lmer}}, \code{\link{ss_power_mixed_effects}},
#'   \code{\link[lme4]{lmer}}, \code{\link[lme4]{VarCorr}}
#'
#' @examples
#' fit <- lme4::lmer(Reaction ~ Days + (Days | Subject),
#'                   data = lme4::sleepstudy)
#'
#' # Marginal R2 is the proportion of variance the fixed effects account for;
#' # conditional R2 adds what the random effects account for, so the gap
#' # between the two is what the subject-level terms buy.
#' R2_mixed_effects(fit)
#'
#' # A parametric bootstrap percentile interval for both quantities. Each
#' # replication refits the model. B = 20 keeps the example quick and is
#' # below the 400 a reported interval needs; the option lets it run and
#' # the warning says so. A reported interval deserves the default
#' # B = 1000. The seed makes the limits reproducible and leaves the
#' # caller's generator state as it was.
#' old <- options(dmar.bootstrap_exploratory = TRUE)
#' R2_mixed_effects(fit, ci_method = "boot", B = 20, seed = 113)
#' options(old)
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @keywords htest
#'
#' @family agreement and measurement
#' @family mixed models
#'
#' @export

R2_mixed_effects <- function(model, conf_level = 0.95, ci_method = c("none", "boot"),
                     B = 1000, seed = NULL) {
  ci_method <- match.arg(ci_method)
  if (!is.numeric(conf_level) || length(conf_level) != 1L ||
      is.na(conf_level) || conf_level <= 0 || conf_level >= 1)
    .dmar_refuse_argument("'conf_level' must be in (0, 1).", arg = "conf_level")
  # The number of replications follows the rule of ?dmar_bootstrap when
  # the interval is requested: 10 beyond each limit (400 at the 95 percent
  # level). Through DMAR 1.0.x B was not checked.
  if (ci_method == "boot") {
    B <- .dmar_check_B(B, tail = (1 - conf_level) / 2,
                       conf_level = conf_level)
  }

  parts <- .R2_mixed_effects_parts(model)
  total <- parts$var_fixed + parts$var_random + parts$var_resid
  R2_marginal    <- parts$var_fixed / total
  R2_conditional <- (parts$var_fixed + parts$var_random) / total

  if (ci_method == "none") {
    out <- data.frame(
      term  = c("R2_marginal", "R2_conditional"),
      value = c(R2_marginal, R2_conditional),
      stringsAsFactors = FALSE,
      row.names = NULL
    )
    out <- .as_dmar_tbl(out)
    attr(out, "conf_level") <- conf_level
    class(out) <- c("dmar_R2_mixed_effects", class(out))
    return(out)
  }

  # ci_method == "boot": parametric bootstrap of the two R2 values.
  if (!inherits(model, "lmerMod"))
    .dmar_refuse_argument(
      "ci_method = 'boot' is supported only for lme4::lmer ('lmerMod') models.",
      arg = NULL)
  if (!requireNamespace("lme4", quietly = TRUE))
    .dmar_refuse(c(
      "The 'lme4' package is required. Install it with: ",
      "install.packages('lme4')"), "package_required", package = "lme4")

  # A supplied seed is set for this call and the caller's generator state
  # is restored on exit; NULL leaves the generator alone.
  .dmar_local_seed(seed)

  stat <- function(m) {
    p <- .R2_mixed_effects_parts(m)
    tot <- p$var_fixed + p$var_random + p$var_resid
    c(p$var_fixed / tot, (p$var_fixed + p$var_random) / tot)
  }
  bb <- lme4::bootMer(model, stat, nsim = B, seed = NULL)
  alpha <- 1 - conf_level
  probs <- c(alpha / 2, 1 - alpha / 2)
  # A replicate with a value that is not a finite number for either
  # quantity (a refit that failed, or variance components that could not
  # be computed) is dropped for both, counted, and reported, so the two
  # intervals rest on the same replicates. Through DMAR 1.0.x such values
  # were removed quantity by quantity, without a count or a warning.
  draws <- matrix(as.numeric(bb$t), ncol = 2L)
  keep <- rowSums(!is.finite(draws)) == 0L
  verdict <- .dmar_boot_settle(
    B, sum(keep), tail = alpha / 2, conf_level = conf_level,
    fun = "R2_mixed_effects",
    quantity = "the marginal or the conditional R-squared",
    dropped_because = paste0("gave a refit whose variance components ",
                             "could not be computed"))
  q_marg <- q_cond <- c(NA_real_, NA_real_)
  if (verdict$formed) {
    q_marg <- stats::quantile(draws[keep, 1L], probs = probs, names = FALSE)
    q_cond <- stats::quantile(draws[keep, 2L], probs = probs, names = FALSE)
  }

  out <- data.frame(
    term  = c("R2_marginal", "R2_marginal_ci_lower", "R2_marginal_ci_upper",
              "R2_conditional", "R2_conditional_ci_lower",
              "R2_conditional_ci_upper"),
    value = c(R2_marginal, q_marg[1], q_marg[2],
              R2_conditional, q_cond[1], q_cond[2]),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  out <- .as_dmar_tbl(out, conf_level = conf_level)
  attr(out, "conf_level") <- conf_level
  out <- .dmar_boot_attrs(out, verdict)
  class(out) <- c("dmar_R2_mixed_effects", class(out))
  out
}


# Variance decomposition shared by the point estimate and the bootstrap
# statistic. Returns a list with var_fixed, var_random, var_resid.
# var_random uses the Johnson (2014) quadratic-form contribution per
# random-effects term, which reduces to the intercept variance component
# for a random-intercept term and handles random slopes correctly.
.R2_mixed_effects_parts <- function(model) {
  if (inherits(model, "lmerMod")) {
    if (!requireNamespace("lme4", quietly = TRUE))
      .dmar_refuse(c(
        "The 'lme4' package is required. Install it with: ",
        "install.packages('lme4')"), "package_required", package = "lme4")
    beta <- lme4::fixef(model)
    var_fixed <- stats::var(as.vector(stats::model.matrix(model) %*% beta))
    var_resid <- stats::sigma(model)^2

    vc <- lme4::VarCorr(model)
    # VarCorr(model) and getME(model, "mmList") are both ordered by
    # random-effects term (bar), so the i-th covariance component aligns with
    # the i-th design matrix positionally. Indexing VarCorr by grouping-factor
    # name is wrong when a factor appears in more than one term, for example
    # (1 | g) + (0 + x | g) or the (x || g) double-bar form that lme4 expands
    # to it: VarCorr names the components 'g' and 'g.1', but cnms repeats 'g',
    # so a name lookup selects the intercept covariance for both terms and
    # never touches the slope component.
    mm <- lme4::getME(model, "mmList")
    var_random <- 0
    for (i in seq_along(mm)) {
      Z <- mm[[i]]
      Sigma <- vc[[i]]
      var_random <- var_random + mean(rowSums((Z %*% Sigma) * Z))
    }
    return(list(var_fixed = var_fixed, var_random = var_random,
                var_resid = var_resid))
  }

  if (inherits(model, "lme")) {
    if (!requireNamespace("nlme", quietly = TRUE))
      .dmar_refuse(c(
        "The 'nlme' package is required. Install it with: ",
        "install.packages('nlme')"), "package_required", package = "nlme")
    # Variance of the fixed-effect linear predictor: predict at level 0 gives
    # the fitted values from the fixed effects only.
    var_fixed <- stats::var(as.vector(stats::predict(model, level = 0)))
    var_resid <- model[["sigma"]]^2

    # nlme random-effects covariance on the response scale, as a matrix (single
    # grouping level) or a list of grouping-level matrices (nested). Build the
    # matching random-effects design from the model's data and apply the
    # Johnson (2014) quadratic-form contribution per level.
    vlist <- nlme::getVarCov(model, type = "random.effects")
    if (!is.list(vlist)) vlist <- list(vlist)
    dat <- nlme::getData(model)
    var_random <- 0
    for (lv in seq_along(vlist)) {
      Sigma <- vlist[[lv]]
      cn <- colnames(Sigma)
      slopes <- setdiff(cn, "(Intercept)")
      if (length(slopes) == 0L) {
        Z <- matrix(1, nrow = nrow(dat), ncol = 1L,
                    dimnames = list(NULL, "(Intercept)"))
      } else {
        Z <- stats::model.matrix(stats::reformulate(slopes), data = dat)
      }
      Z <- Z[, cn, drop = FALSE]
      var_random <- var_random + mean(rowSums((Z %*% Sigma) * Z))
    }
    return(list(var_fixed = var_fixed, var_random = var_random,
                var_resid = var_resid))
  }

  .dmar_refuse_argument(
    "'model' must be a fitted lme4::lmer ('lmerMod') or nlme::lme model.",
    arg = "model")
}


#' Tidy / Glance Methods for R2_mixed_effects Output
#'
#' Returns the broom-style summary of the marginal and conditional
#' \eqn{R^2}: \code{term} (\code{"R2_marginal"} or
#' \code{"R2_conditional"}), \code{estimate}, and, when a bootstrap
#' interval was requested, \code{ci_lower}, \code{ci_upper}, and
#' \code{conf_level}. \code{glance} returns the same information in one
#' row per quantity.
#'
#' @param x A \code{dmar_R2_mixed_effects} object returned by
#'   \code{\link{R2_mixed_effects}}.
#' @param \dots Unused.
#' @return A \code{data.frame} in broom convention.
#' @author Ken Kelley \email{kkelley@@nd.edu}
#' @examples
#' fit <- lme4::lmer(Reaction ~ Days + (Days | Subject),
#'                   data = lme4::sleepstudy)
#' res <- R2_mixed_effects(fit)
#' generics::tidy(res)
#' generics::glance(res)
#' @importFrom generics tidy
#' @exportS3Method generics::tidy
tidy.dmar_R2_mixed_effects <- function(x, ...) {
  est_marg <- x$value[x$term == "R2_marginal"]
  est_cond <- x$value[x$term == "R2_conditional"]
  has_ci <- any(x$term == "R2_marginal_ci_lower")
  if (!has_ci) {
    return(data.frame(
      term     = c("R2_marginal", "R2_conditional"),
      estimate = c(est_marg, est_cond),
      stringsAsFactors = FALSE,
      row.names = NULL
    ))
  }
  data.frame(
    term       = c("R2_marginal", "R2_conditional"),
    estimate   = c(est_marg, est_cond),
    ci_lower   = c(x$value[x$term == "R2_marginal_ci_lower"],
                   x$value[x$term == "R2_conditional_ci_lower"]),
    ci_upper  = c(x$value[x$term == "R2_marginal_ci_upper"],
                   x$value[x$term == "R2_conditional_ci_upper"]),
    conf_level = attr(x, "conf_level"),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
}

#' @rdname tidy.dmar_R2_mixed_effects
#' @importFrom generics glance
#' @exportS3Method generics::glance
glance.dmar_R2_mixed_effects <- function(x, ...) {
  tidy.dmar_R2_mixed_effects(x, ...)
}
