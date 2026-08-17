# Confidence interval for generalized eta squared (effect size for ANOVA).
#' Confidence Interval for Generalized Eta Squared (Approximate)
#'
#' Returns the point estimate of generalized eta squared (\eqn{\eta^2_G};
#' Olejnik & Algina, 2003) along with an \emph{optional} confidence interval
#' computed by one of two approximate methods. The default is to return only
#' the point estimate (\code{method = "none"}), because both available CI
#' methods are approximations whose coverage properties have not been broadly
#' validated for this estimand and warrant independent evaluation before being
#' used in substantive inference.
#'
#' @param object Optional. A fitted model object of class
#'   \code{\link[stats]{aov}}, \code{\link[stats]{lm}}, or \code{aovlist}
#'   (multi-stratum aov fit for within-subjects / mixed designs).
#' @param observed Character vector of factor names treated as measured.
#' @param SS_effect,SS_observed,SS_error Sums of squares (option 2 in
#'   \code{\link{eta_squared_generalized}}).
#' @param F_effect,df_effect,F_observed,df_observed,df_error \emph{F}-values
#'   and degrees of freedom (option 3 in \code{\link{eta_squared_generalized}}).
#' @param N Total sample size. Required when \code{method = "parametric"} and
#'   no fitted model is supplied; ignored when \code{method = "none"} or
#'   when a fitted model is supplied (derived automatically, via
#'   \code{\link[stats]{nobs}(object)} for single-stratum fits, or by
#'   summing degrees of freedom across error strata for \code{aovlist}).
#' @param method One of \code{"none"} (default), \code{"parametric"}, or
#'   \code{"bootstrap"}. See Details.
#' @param B Integer. Number of bootstrap replications when
#'   \code{method = "bootstrap"}. Minimum \code{1000} (enforced); default
#'   \code{10000}. We recommend \code{10000} or more for publication-quality
#'   intervals.
#' @param conf_level Desired confidence coverage; default \code{0.95}.
#' @param alpha_lower,alpha_upper Optional Type I error on the lower and upper
#'   side.
#' @param seed Optional integer seed for the bootstrap, for
#'   reproducibility. Used locally: the caller's random number generator
#'   state is restored on exit. Default \code{NULL} leaves the random
#'   number generator state alone.
#'
#' @return A \code{data.frame} with one row per focal effect and the
#'   columns \code{effect}, \code{eta_squared_generalized},
#'   \code{lower_limit}, \code{upper_limit}, and \code{method}. For
#'   \code{aovlist} fits a \code{stratum} column is also present. When
#'   \code{method = "none"}, the limit columns contain \code{NA}.
#'
#' @details
#' \strong{Why CI = "none" is the default.} Confidence interval construction
#' for \eqn{\eta^2_G} is not as settled as for partial \eqn{\eta^2} or
#' \eqn{\omega^2}, because the denominator mixes sums of squares from
#' heterogeneous sources (the focal effect, one or more measured factors, and
#' the error term). No noncentral \emph{F} transformation maps the population
#' noncentrality parameter directly to \eqn{\eta^2_G}. Both methods below are
#' approximations and are exposed for exploration rather than as defaults.
#'
#' \strong{\code{method = "parametric"}.} The function first obtains a
#' confidence interval for the population noncentrality parameter
#' \eqn{\lambda} of the focal effect's \emph{F}-test via
#' \code{\link{conf_limits_ncf}}. The NCP bounds are mapped through the
#' partial-\eqn{\eta^2} transformation
#' \eqn{\eta^2_{p,\text{bound}} = \lambda_{\text{bound}}/(\lambda_{\text{bound}} + N)}
#' (matching the convention used by \code{\link{ci_pvaf}} and
#' \code{\link{ci_omega_squared}}), and then re-expressed as \eqn{\eta^2_G}
#' bounds via
#' \deqn{\eta^2_{G,\text{bound}} = \frac{r_{\text{bound}}}{r_{\text{bound}} + r_{\text{obs}} + 1},}
#' where \eqn{r_{\text{bound}} = \eta^2_{p,\text{bound}}/(1-\eta^2_{p,\text{bound}})}
#' and \eqn{r_{\text{obs}} = \sum \mathit{SS}_{\text{measured}}/\mathit{SS}_{\text{error}}}.
#' This treats the observed-factor sums of squares as fixed at their sample
#' values, so the interval inherits whatever sampling variability those
#' contribute. It has not been validated for coverage and should be treated as
#' preliminary.
#'
#' \strong{\code{method = "bootstrap"}.} A residual bootstrap from the
#' fitted model: the resampling unit is a residual, drawn
#' nonparametrically from the model's own residuals rather than from a
#' fitted distribution. For each of the \code{B} replications, the
#' response is regenerated as \eqn{\hat{y}_i + \varepsilon^*_i} where
#' \eqn{\varepsilon^*} is sampled with replacement from the model's
#' residuals; the model is refit; \eqn{\eta^2_G} is recomputed; and the
#' percentile interval, the empirical quantiles of the \code{B}
#' bootstrap values (Efron & Tibshirani, 1993), is reported. The
#' percentile interval is the only bootstrap interval offered; there is
#' no bias-corrected and accelerated (BCa) variant. Replicates whose
#' refit fails are dropped, and the interval is computed from the
#' replications that return a value. Bootstrap results vary from run to
#' run; supply \code{seed} for reproducibility. Requires a
#' single-stratum \code{aov}/\code{lm} fit. \strong{Not yet supported
#' for \code{aovlist} (multi-stratum / within-subjects) fits}, since
#' the bootstrap needs to respect the subject-level correlation
#' structure, which naive residual resampling does not. Use
#' \code{method = "parametric"} for \code{aovlist} fits.
#' Coverage has not been broadly validated for this estimand.
#'
#' \strong{Within-subjects designs (\code{aovlist}).} For multi-stratum
#' fits the parametric CI uses each focal effect's stratum-specific
#' \emph{F} test and degrees of freedom. The denominator ratio is
#' adjusted to reflect the full set of error strata: the implied
#' \eqn{\mathit{SS}_{\text{error}}} for the focal effect is its own
#' stratum's residual SS, while the \eqn{r_{\text{obs}}} term includes
#' both measured-factor SS and the residual SS of all \emph{other}
#' strata. This generalizes the partial-\eqn{\eta^2} CI machinery to the
#' Bakeman (2005) denominator.
#'
#' @references
#' Algina, J., Keselman, H. J., & Penfield, R. D. (2005). An alternative to
#'   Cohen's standardized mean difference effect size: A robust parameter and
#'   confidence interval in the two independent groups case.
#'   \emph{Psychological Methods, 10}(3), 317--328.
#'   \doi{10.1037/1082-989X.10.3.317}
#'
#' Bakeman, R. (2005). Recommended effect size statistics for repeated
#'   measures designs. \emph{Behavior Research Methods, 37}(3), 379--384.
#'   \doi{10.3758/BF03192707}
#'
#' Efron, B., & Tibshirani, R. J. (1993). \emph{An introduction to the
#'   bootstrap}. New York, NY: Chapman & Hall/CRC.
#'
#' Kelley, K. (2007). Confidence intervals for standardized effect
#'   sizes: Theory, application, and implementation. \emph{Journal of Statistical
#'   Software, 20}(8), 1--24. \doi{10.18637/jss.v020.i08}
#'
#' Kelley, K., & Preacher, K. J. (2012). On effect size.
#'   \emph{Psychological Methods, 17}, 137--152. \doi{10.1037/a0028086}
#'
#' Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). \emph{Designing
#'   experiments and analyzing data: A model comparison perspective}
#'   (4th ed.). Routledge. (See Chapter 3 on \eqn{\eta^2}, Chapter 7 on
#'   factorial designs, and Chapter 11 on generalized \eqn{\eta^2} for
#'   within-subjects designs.)
#'
#' Olejnik, S., & Algina, J. (2003). Generalized eta and omega squared
#'   statistics: Measures of effect size for some common research designs.
#'   \emph{Psychological Methods, 8}(4), 434--447.
#'   \doi{10.1037/1082-989X.8.4.434}
#'
#' Smithson, M. (2001). Correct confidence intervals for various regression
#'   effect sizes and parameters: The importance of noncentral distributions in
#'   computing intervals. \emph{Educational and Psychological Measurement, 61},
#'   605--632. \doi{10.1177/00131640121971392}
#'
#' Steiger, J. H. (2004). Beyond the \emph{F} test: Effect size confidence
#'   intervals and tests of close fit in the analysis of variance and contrast
#'   analysis. \emph{Psychological Methods, 9}(2), 164--182.
#'   \doi{10.1037/1082-989X.9.2.164}
#'
#' @seealso \code{\link{eta_squared_generalized}}, \code{\link{ci_eta_squared}}
#'
#' @examples
#' # The pygmalion expectancy experiment: treatment is manipulated, while
#' # grade is a measured classification, so grade belongs in the denominator.
#' pyg <- pygmalion
#' pyg$grade <- factor(pyg$grade)
#' fit <- aov(iq_8 ~ treatment * grade, data = pyg)
#'
#' # The default returns the point estimate and leaves the limits NA, because
#' # both interval methods are approximations.
#' ci_eta_squared_generalized(fit, observed = "grade")
#'
#' # The parametric approximation maps a noncentral F interval for the focal
#' # effect through the observed sums of squares. Warnings mark the method as
#' # preliminary and report that the interaction's lower limit is clamped at
#' # 0; treat the limits accordingly.
#' ci_eta_squared_generalized(fit, observed = "grade", method = "parametric")
#'
#' # The third option is a residual bootstrap, which refits the model once per
#' # replication. It is not run here: B is required to be at least 1000, and
#' # even that is slower than an example should be, while a reported interval
#' # deserves B = 10000 or more. The call is
#' #   ci_eta_squared_generalized(fit, observed = "grade",
#' #                              method = "bootstrap", B = 10000, seed = 113)
#'
#' # Within-subjects ANOVA. The parametric CI uses each effect's own stratum.
#' set.seed(113)
#' n <- 20
#' rm_data <- data.frame(
#'   subject = factor(rep(seq_len(n), each = 3)),
#'   time    = factor(rep(c("Pre", "Mid", "Post"), n),
#'                    levels = c("Pre", "Mid", "Post")),
#'   y       = rnorm(n, sd = 1.5)[rep(seq_len(n), each = 3)] +
#'             0.7 * rep(1:3, n) + rnorm(n * 3, sd = 1.2)
#' )
#' fit_rm <- aov(y ~ time + Error(subject/time), data = rm_data)
#' ci_eta_squared_generalized(fit_rm, method = "parametric")
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @keywords htest design
#'
#' @family confidence intervals for effect sizes
#'
#' @export

ci_eta_squared_generalized <- function(
  object       = NULL,
  observed     = NULL,
  SS_effect    = NULL,
  SS_observed  = NULL,
  SS_error     = NULL,
  F_effect     = NULL,
  df_effect    = NULL,
  F_observed   = NULL,
  df_observed  = NULL,
  df_error     = NULL,
  N            = NULL,
  method       = c("none", "parametric", "bootstrap"),
  B            = 10000L,
  conf_level   = 0.95,
  alpha_lower  = NULL,
  alpha_upper  = NULL,
  seed         = NULL
) {
  method <- match.arg(method)

  point <- eta_squared_generalized(
    object      = object,      observed    = observed,
    SS_effect   = SS_effect,   SS_observed = SS_observed, SS_error = SS_error,
    F_effect    = F_effect,    df_effect   = df_effect,
    F_observed  = F_observed,  df_observed = df_observed, df_error = df_error
  )
  point$lower_limit <- NA_real_
  point$upper_limit <- NA_real_
  point$method      <- method

  if (method == "none") {
    message(
      "No CI computed (method = 'none'). Set method = 'parametric' for an ",
      "approximate noncentral F transformation, or method = 'bootstrap' for a ",
      "residual bootstrap (requires a fitted model). Both methods are ",
      "approximate; see ?ci_eta_squared_generalized."
    )
    # Wide-format dmar_tbl: a leading `effect` label column beside the
    # eta_squared_generalized estimate and its (here NA) interval limits. No
    # conf_level is attached because method = "none" computes no interval, so a
    # confidence-level footer would be misleading.
    return(.as_dmar_tbl(point, subclass = "dmar_ci_anova"))
  }

  alphas <- .ci_eta_squared_alphas(conf_level, alpha_lower, alpha_upper)

  if (method == "parametric") {
    warning(
      "Parametric CI uses an approximate transformation that maps the partial-",
      "eta squared CI through the observed sums of squares for measured factors. ",
      "Coverage properties have not been broadly validated; results should be ",
      "treated as preliminary and independently evaluated."
    )
    # Wide-format dmar_tbl with the CI present; the conf_level footer applies.
    # One clamp warning per call, not one per clamped effect.
    return(.as_dmar_tbl(
      .warn_ncf_clamp_once(.ci_eta_g_parametric(
        point, object, observed,
        SS_effect, SS_observed, SS_error,
        F_effect, df_effect, F_observed, df_observed, df_error, N,
        alphas$alpha_lower, alphas$alpha_upper
      )),
      conf_level = conf_level
    , subclass = "dmar_ci_anova"))
  }

  if (method == "bootstrap") {
    if (is.null(object)) {
      stop("method = 'bootstrap' requires a fitted model ('object'); the raw SS / F-df interfaces cannot bootstrap.")
    }
    if (inherits(object, "aovlist")) {
      stop("method = 'bootstrap' is not yet supported for aovlist (multi-stratum / within-subjects) fits. Use method = 'parametric' for an approximate CI on aovlist models.")
    }
    if (!is.numeric(B) || length(B) != 1L || B < 1000)
      stop("'B' must be a single integer >= 1000. We recommend 10000 or more.")
    warning(sprintf(
      "Bootstrap CI resamples the fitted model's residuals nonparametrically (B = %d). Coverage has not been broadly validated for generalized eta squared; results should be treated as preliminary.",
      as.integer(B)
    ))
    # Wide-format dmar_tbl with the bootstrap CI present; conf_level footer applies.
    return(.as_dmar_tbl(
      .ci_eta_g_bootstrap(
        object, observed, as.integer(B),
        alphas$alpha_lower, alphas$alpha_upper, seed
      ),
      conf_level = conf_level
    , subclass = "dmar_ci_anova"))
  }
}


# Parametric approximate CI: map NCP bound -> partial eta_sq -> generalized.
.ci_eta_g_parametric <- function(point, object, observed,
                                 SS_effect, SS_observed, SS_error,
                                 F_effect, df_effect,
                                 F_observed, df_observed, df_error, N,
                                 alpha_lower, alpha_upper) {

  per_effect <- function(eff_name, F_v, df_e, df_err, ss_obs_ratio, N_use) {
    ncp <- .conf_limits_ncf_for(
      caller      = "ci_eta_squared_generalized",
      quantity    = "generalized eta squared",
      F_value     = F_v,     df_1 = df_e, df_2 = df_err,
      alpha_lower = alpha_lower, alpha_upper = alpha_upper,
      conf_level  = NULL
    )
    lo_ncp <- ncp$value[ncp$term == "lower_limit"]
    up_ncp <- ncp$value[ncp$term == "upper_limit"]

    lo_etap <- lo_ncp / (lo_ncp + N_use)
    # A NA upper NCP means the observed F is so small that even at lambda = 0 the
    # lower-tail probability is already at or below alpha_upper, leaving the upper
    # noncentrality limit undefined (conf_limits_ncf documents this). Propagate NA
    # to the upper generalized eta squared limit, matching ci_eta_squared.
    up_etap <- if (is.na(up_ncp)) NA_real_
               else if (is.infinite(up_ncp)) 1
               else up_ncp / (up_ncp + N_use)

    lo_r <- if (lo_etap <= 0) 0 else lo_etap / (1 - lo_etap)
    up_r <- if (is.na(up_etap)) NA_real_
            else if (up_etap >= 1) Inf
            else up_etap / (1 - up_etap)

    lo_eta_g <- lo_r / (lo_r + ss_obs_ratio + 1)
    up_eta_g <- if (is.na(up_r)) NA_real_
                else if (is.infinite(up_r)) 1
                else up_r / (up_r + ss_obs_ratio + 1)

    list(lower = lo_eta_g, upper = up_eta_g)
  }

  if (!is.null(object)) {
    if (inherits(object, "aovlist")) {
      effects_tbl <- .aovlist_effects_table(object)
      total_ss_err <- .aovlist_total_residual_ss(object)
      effect_names <- effects_tbl$effect
      ss_lookup    <- setNames(effects_tbl$SS_effect, effect_names)
      observed_in  <- .expand_observed(observed, effect_names)
      N_use        <- .aovlist_nobs(object)

      for (i in seq_len(nrow(point))) {
        eff   <- point$effect[i]
        row   <- effects_tbl[effects_tbl$effect == eff, , drop = FALSE]
        F_v   <- row$F_value
        df_e  <- row$df_effect
        df_err <- row$df_error
        ss_err_focal <- row$SS_error

        obs_others <- setdiff(observed_in, eff)
        ss_obs <- if (length(obs_others) == 0L) 0 else sum(ss_lookup[obs_others])
        # Per-effect denominator ratio relative to the FOCAL stratum's error,
        # NOT relative to total SS error. The per_effect helper appends
        # `+ 1` for the focal stratum's residual, so subtract that here.
        ss_obs_ratio <- (ss_obs + total_ss_err - ss_err_focal) / ss_err_focal

        lims <- per_effect(eff, F_v, df_e, df_err, ss_obs_ratio, N_use)
        point$lower_limit[i] <- lims$lower
        point$upper_limit[i] <- lims$upper
      }
      return(point)
    }

    tbl    <- stats::anova(object)
    ss_err <- tbl["Residuals", "Sum Sq"]
    df_err <- tbl["Residuals", "Df"]
    N_use  <- stats::nobs(object)
    observed_in <- .expand_observed(observed, setdiff(rownames(tbl), "Residuals"))

    for (i in seq_len(nrow(point))) {
      eff <- point$effect[i]
      F_v <- tbl[eff, "F value"]
      df_e <- tbl[eff, "Df"]
      obs_others <- setdiff(observed_in, eff)
      ss_obs <- if (length(obs_others) == 0L) 0 else sum(tbl[obs_others, "Sum Sq"])
      lims <- per_effect(eff, F_v, df_e, df_err, ss_obs / ss_err, N_use)
      point$lower_limit[i] <- lims$lower
      point$upper_limit[i] <- lims$upper
    }
    return(point)
  }

  if (!is.null(F_effect)) {
    if (is.null(N)) {
      stop("Parametric CI from the F/df interface requires 'N' (total sample size).")
    }
    if (!is.numeric(N) || N <= df_effect + df_error)
      stop("'N' must exceed df_effect + df_error.")
    ss_obs_ratio <- if (is.null(F_observed)) 0 else sum(F_observed * df_observed / df_error)
    lims <- per_effect("overall", F_effect, df_effect, df_error, ss_obs_ratio, N)
    point$lower_limit[1] <- lims$lower
    point$upper_limit[1] <- lims$upper
    return(point)
  }

  # SS-only path: need df_effect, df_error, N too. Reconstruct F.
  if (is.null(df_effect) || is.null(df_error) || is.null(N)) {
    stop(
      "Parametric CI from the SS interface additionally requires 'df_effect', ",
      "'df_error', and 'N' so that the noncentral F transformation can be applied."
    )
  }
  F_v <- (SS_effect / df_effect) / (SS_error / df_error)
  ss_obs_ratio <- sum(SS_observed) / SS_error
  lims <- per_effect("overall", F_v, df_effect, df_error, ss_obs_ratio, N)
  point$lower_limit[1] <- lims$lower
  point$upper_limit[1] <- lims$upper
  point
}


# Residual bootstrap CI from a fitted model.
.ci_eta_g_bootstrap <- function(object, observed, B,
                                alpha_lower, alpha_upper, seed) {
  if (!is.null(seed)) {
    if (exists(".Random.seed", envir = .GlobalEnv)) {
      .old_seed <- get(".Random.seed", envir = .GlobalEnv)
      on.exit(assign(".Random.seed", .old_seed, envir = .GlobalEnv), add = TRUE)
    } else {
      on.exit(if (exists(".Random.seed", envir = .GlobalEnv)) rm(list = ".Random.seed", envir = .GlobalEnv), add = TRUE)
    }
    set.seed(seed)
  }
  fitted_vals <- stats::fitted(object)
  resid_vals  <- stats::residuals(object)
  if (length(fitted_vals) != length(resid_vals)) {
    stop("Fitted values and residuals have inconsistent lengths.")
  }
  model_data <- stats::model.frame(object)
  resp_name  <- names(model_data)[1L]

  base_point <- eta_squared_generalized(object = object, observed = observed)
  effect_names <- base_point$effect
  K <- length(effect_names)

  boot_mat <- matrix(NA_real_, nrow = B, ncol = K)
  colnames(boot_mat) <- effect_names

  for (b in seq_len(B)) {
    eps <- sample(resid_vals, replace = TRUE)
    y_boot <- fitted_vals + eps
    boot_data <- model_data
    boot_data[[resp_name]] <- y_boot
    fit_b <- try(stats::update(object, data = boot_data), silent = TRUE)
    if (inherits(fit_b, "try-error")) next
    eg_b <- try(eta_squared_generalized(object = fit_b, observed = observed),
                silent = TRUE)
    if (inherits(eg_b, "try-error")) next
    # Align by effect name
    idx <- match(effect_names, eg_b$effect)
    boot_mat[b, ] <- eg_b$eta_squared_generalized[idx]
  }

  point <- base_point
  point$lower_limit <- apply(boot_mat, 2L,
                             function(x) stats::quantile(x, probs = alpha_lower, na.rm = TRUE))
  point$upper_limit <- apply(boot_mat, 2L,
                             function(x) stats::quantile(x, probs = 1 - alpha_upper, na.rm = TRUE))
  point$method <- "bootstrap"
  rownames(point) <- NULL
  point
}
