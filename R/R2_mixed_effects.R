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
#' @param B Number of bootstrap replications when
#'   \code{ci_method = "boot"}. Default \code{1000}.
#' @param seed Optional integer seed for the bootstrap. Default
#'   \code{NULL}, meaning the caller's current RNG state is used and left
#'   unchanged. When supplied, the seed is set inside the function and
#'   the caller's RNG state is restored on exit.
#' @param \dots Currently unused.
#'
#' @return A \code{data.frame} with rows \code{"R2_marginal"} and
#'   \code{"R2_conditional"} in the \code{value} column. When
#'   \code{ci_method = "boot"}, lower- and upper-limit rows for each
#'   quantity are appended and the confidence level is carried on the
#'   object.
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
#' Replicates whose refit fails are dropped, and the interval is
#' computed from the replications that return a value. The default
#' \code{B = 1000} is adequate for the central quantiles a percentile
#' interval uses; raising it tightens the Monte Carlo error of the
#' reported limits. Bootstrap results vary from run to run; supply
#' \code{seed} for reproducibility.
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
#' # A bootstrap interval is available through ci_method = "boot". It refits
#' # the model once per replication, so it is not run here; the call is
#' #   R2_mixed_effects(fit, ci_method = "boot", B = 1000, seed = 113)
#' # where B = 1000 is the default and seed is supplied because the limits
#' # otherwise move from run to run. Raise B when the Monte Carlo error of
#' # the reported limits needs to be smaller.
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
                     B = 1000, seed = NULL, ...) {
  ci_method <- match.arg(ci_method)
  if (conf_level <= 0 || conf_level >= 1)
    stop("'conf_level' must be in (0, 1).")

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
    stop("ci_method = 'boot' is supported only for lme4::lmer ('lmerMod') models.")
  if (!requireNamespace("lme4", quietly = TRUE))
    stop("The 'lme4' package is required. Install it with: ",
         "install.packages('lme4')")

  if (!is.null(seed)) {
    if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
      on.exit(assign(".Random.seed", old_seed, envir = .GlobalEnv), add = TRUE)
    } else {
      on.exit(if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
        rm(".Random.seed", envir = .GlobalEnv), add = TRUE)
    }
    set.seed(seed)
  }

  stat <- function(m) {
    p <- .R2_mixed_effects_parts(m)
    tot <- p$var_fixed + p$var_random + p$var_resid
    c(p$var_fixed / tot, (p$var_fixed + p$var_random) / tot)
  }
  bb <- lme4::bootMer(model, stat, nsim = B, seed = NULL)
  alpha <- 1 - conf_level
  probs <- c(alpha / 2, 1 - alpha / 2)
  q_marg <- stats::quantile(bb$t[, 1], probs = probs, na.rm = TRUE, names = FALSE)
  q_cond <- stats::quantile(bb$t[, 2], probs = probs, na.rm = TRUE, names = FALSE)

  out <- data.frame(
    term  = c("R2_marginal", "R2_marginal_lower", "R2_marginal_upper",
              "R2_conditional", "R2_conditional_lower", "R2_conditional_upper"),
    value = c(R2_marginal, q_marg[1], q_marg[2],
              R2_conditional, q_cond[1], q_cond[2]),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  out <- .as_dmar_tbl(out, conf_level = conf_level)
  attr(out, "conf_level") <- conf_level
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
      stop("The 'lme4' package is required. Install it with: ",
           "install.packages('lme4')")
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
      stop("The 'nlme' package is required. Install it with: ",
           "install.packages('nlme')")
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

  stop("'model' must be a fitted lme4::lmer ('lmerMod') or nlme::lme model.")
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
  has_ci <- any(x$term == "R2_marginal_lower")
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
    ci_lower   = c(x$value[x$term == "R2_marginal_lower"],
                   x$value[x$term == "R2_conditional_lower"]),
    ci_upper  = c(x$value[x$term == "R2_marginal_upper"],
                   x$value[x$term == "R2_conditional_upper"]),
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
