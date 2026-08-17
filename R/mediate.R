#' Mediation Analysis With Bootstrap Confidence Intervals
#'
#' Estimates the simple mediation model, a predictor \eqn{X} affecting an
#' outcome \eqn{Y} directly and through a mediator \eqn{M}, and reports the
#' indirect, direct, and total effects with confidence intervals in one tidy
#' table. The indirect effect \eqn{a b} is the product of the \eqn{X \to M}
#' path and the \eqn{M \to Y} path (holding \eqn{X}), and its sampling
#' distribution is skewed, which is why the default interval is the
#' percentile bootstrap rather than a normal approximation; the
#' bias-corrected and accelerated (BCa) bootstrap, the Monte Carlo
#' (parametric simulation) interval, and the Sobel normal-theory interval
#' are available for comparison. This is the analysis counterpart of the
#' planning functions \code{\link{ss_aipe_indirect_effect}} and
#' \code{\link{ss_power_indirect_effect}}.
#'
#' @param data A \code{data.frame} containing the variables.
#' @param x,m,y Names (single character strings) of the predictor, the
#'   mediator, and the outcome columns in \code{data}.
#' @param covariates Optional character vector of covariate column names,
#'   entered in both the mediator and the outcome models.
#' @param ci_method Confidence interval method for the indirect effect:
#'   \code{"boot_percentile"} (default), \code{"boot_bca"},
#'   \code{"monte_carlo"} (simulate \eqn{a} and \eqn{b} from their joint
#'   normal approximation; MacKinnon, Lockwood, & Williams, 2004), or
#'   \code{"sobel"} (the first-order normal-theory interval; reported for
#'   comparison, not recommended for inference). Direct and total effects
#'   always carry their ordinary \emph{t}-based intervals.
#' @param B Number of bootstrap or Monte Carlo replications. Defaults to
#'   2000; published analyses often use 5000 or more, and the BCa
#'   interval in particular rewards a large \code{B} (see
#'   \emph{Details}).
#' @param conf_level Confidence level. Defaults to 0.95.
#' @param seed Optional integer seed for the resampling, used locally (the
#'   caller's random number generator state is restored on exit). Default
#'   \code{NULL} leaves the random number generator state alone.
#'
#' @details
#' The model is the standard pair of regressions
#' \deqn{M = i_M + a X + \mathbf{g}'\mathbf{C} + e_M, \qquad
#'       Y = i_Y + c' X + b M + \mathbf{h}'\mathbf{C} + e_Y,}
#' with \eqn{\mathbf{C}} the optional covariates. The indirect effect is
#' \eqn{a b}, the direct effect \eqn{c'}, and the total effect
#' \eqn{c = c' + a b} (an identity in linear models with the same cases,
#' which the implementation exploits as an internal consistency check).
#'
#' Bootstrap intervals resample cases (rows) with replacement \code{B}
#' times, refitting both regressions in each resample (Efron &
#' Tibshirani, 1993). \code{"boot_percentile"} takes the interval limits
#' from the empirical quantiles of the bootstrapped \eqn{a b} estimates;
#' it is not forced to be symmetric about the estimate, which is the
#' point for a skewed sampling distribution. \code{"boot_bca"} (the
#' bias-corrected and accelerated interval) additionally adjusts the two
#' quantile positions for median bias, estimated from the bootstrap
#' distribution, and for the rate at which the variance of the estimator
#' changes with the parameter, the acceleration, estimated by the
#' jackknife (which adds \emph{N} extra pairs of fits); the adjustments
#' make it second-order accurate where the percentile interval is
#' first-order accurate (DiCiccio & Efron, 1996). Because the adjusted
#' quantile positions sit farther into the tails of the bootstrap
#' distribution, the BCa interval benefits more than the percentile
#' interval does from a \code{B} well above the default. Each resample
#' refits the two regressions by least squares, a closed-form fit with
#' no iterative estimation, so in ordinary data all \code{B}
#' replications enter the interval. A degenerate resample (one whose
#' refit is rank deficient, possible with a near-constant predictor)
#' returns no indirect effect; such replications are dropped with a
#' warning stating how many, and the interval is computed from the
#' replications that returned a value. The Sobel
#' standard error is computed by
#' \code{\link{var_indirect_effect}}. The proportion mediated is one of
#' the effect size measures for mediation models surveyed by Preacher
#' and Kelley (2011); its instability when the total effect is small is
#' why it is reported as \code{NA} near a zero total effect. Listwise
#' deletion is applied to the analysis variables; for full information
#' maximum likelihood under missingness, fit the model in \pkg{lavaan}
#' (see \code{\link{mlmr}} for the package's FIML front end philosophy).
#'
#' Mediation language implies causal structure: with observational data
#' the estimates are conditional associations, and the causal reading
#' requires the usual no-unmeasured-confounding assumptions for both the
#' \eqn{X \to M} and \eqn{M \to Y} links (MacKinnon, 2008). The function
#' computes; the design earns the interpretation.
#'
#' @return A \code{data.frame} (class \code{dmar_tbl}) with rows
#'   \code{indirect_effect} (with its \code{ci_method} interval),
#'   \code{direct_effect}, \code{total_effect}, the paths \code{a} and
#'   \code{b}, their standard errors (\code{se_indirect} per Sobel,
#'   \code{se_a}, \code{se_b}), the interval limits
#'   (\code{indirect_lower} / \code{indirect_upper}, \code{direct_lower} /
#'   \code{direct_upper}, \code{total_lower} / \code{total_upper}),
#'   \code{proportion_mediated} (\eqn{ab/c}; \code{NA} when the total
#'   effect is near zero, where the ratio is unstable), \code{N}, and
#'   \code{B}. The interval method is recorded in the \code{"ci_method"}
#'   attribute and the confidence level in \code{"conf_level"}.
#'
#' @references
#' DiCiccio, T. J., & Efron, B. (1996). Bootstrap confidence intervals.
#'   \emph{Statistical Science, 11}(3), 189--228.
#'
#' Efron, B., & Tibshirani, R. J. (1993). \emph{An introduction to the
#'   bootstrap}. New York, NY: Chapman & Hall/CRC.
#'
#' MacKinnon, D. P. (2008). \emph{Introduction to statistical mediation
#'   analysis}. Erlbaum.
#'
#' MacKinnon, D. P., Lockwood, C. M., & Williams, J. (2004). Confidence
#'   limits for the indirect effect: Distribution of the product and
#'   resampling methods. \emph{Multivariate Behavioral Research, 39}(1),
#'   99--128. \doi{10.1207/s15327906mbr3901_4}
#'
#' Preacher, K. J., & Hayes, A. F. (2008). Asymptotic and resampling
#'   strategies for assessing and comparing indirect effects in multiple
#'   mediator models. \emph{Behavior Research Methods, 40}(3), 879--891.
#'   \doi{10.3758/BRM.40.3.879}
#'
#' Preacher, K. J., & Kelley, K. (2011). Effect size measures for mediation
#'   models: Quantitative strategies for communicating indirect effects.
#'   \emph{Psychological Methods, 16}(2), 93--115. \doi{10.1037/a0022658}
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso \code{\link{ss_aipe_indirect_effect}} and
#'   \code{\link{ss_power_indirect_effect}} for planning the study this
#'   function analyzes; \code{\link{var_indirect_effect}} for the Sobel
#'   variance.
#'
#' @family mediation
#'
#' @keywords models
#'
#' @examples
#' # Simulated mediation: X raises M (a = .5), M raises Y (b = .4), and a
#' # little direct effect remains (c' = .2).
#' set.seed(113)
#' n <- 200
#' x <- rnorm(n)
#' m <- 0.5 * x + rnorm(n, 0, sqrt(1 - 0.25))
#' y <- 0.2 * x + 0.4 * m + rnorm(n, 0, 0.8)
#' d <- data.frame(x = x, m = m, y = y)
#'
#' # The Sobel interval is closed form, so it is the call that runs here.
#' # It assumes the product ab is normally distributed, which is why it is
#' # reported for comparison rather than used for inference. The B row of
#' # the result is NA because no replications are drawn.
#' mediate(d, x = "x", m = "m", y = "y", ci_method = "sobel")
#'
#' # The percentile bootstrap is the default and is what a reported
#' # analysis would use. It is not run here because it refits both
#' # regressions in each of the B resamples; the call is:
#' # mediate(d, x = "x", m = "m", y = "y", seed = 113)
#'
#' # The Monte Carlo interval needs no refitting, so it stays quick even
#' # for a large B, but it still draws B replications. Also not run here:
#' # mediate(d, x = "x", m = "m", y = "y", ci_method = "monte_carlo",
#' #         B = 10000, seed = 113)
#'
#' @export
#' @importFrom stats coef complete.cases lm qnorm quantile rnorm vcov
mediate <- function(data, x, m, y, covariates = NULL,
                    ci_method = c("boot_percentile", "boot_bca",
                                  "monte_carlo", "sobel"),
                    B = 2000, conf_level = 0.95, seed = NULL) {
  ci_method <- match.arg(ci_method)
  if (!is.data.frame(data)) {
    stop("'data' must be a data.frame.", call. = FALSE)
  }
  for (nm in c("x", "m", "y")) {
    val <- get(nm)
    if (!is.character(val) || length(val) != 1L || !(val %in% names(data))) {
      stop(sprintf("'%s' must name one column of 'data'.", nm),
           call. = FALSE)
    }
  }
  if (!is.null(covariates) &&
      (!is.character(covariates) || !all(covariates %in% names(data)))) {
    stop("'covariates' must name columns of 'data'.", call. = FALSE)
  }
  if (!is.numeric(B) || length(B) != 1L || is.na(B) || B < 100 ||
      B != round(B)) {
    stop("'B' must be a single integer of at least 100.", call. = FALSE)
  }
  if (!is.numeric(conf_level) || length(conf_level) != 1L ||
      is.na(conf_level) || conf_level <= 0 || conf_level >= 1) {
    stop("'conf_level' must be a single number in (0, 1).", call. = FALSE)
  }

  vars <- c(x, m, y, covariates)
  dat  <- data[stats::complete.cases(data[, vars, drop = FALSE]),
               vars, drop = FALSE]
  N <- nrow(dat)
  if (N < length(vars) + 5L) {
    stop("Too few complete cases to fit the mediation model.",
         call. = FALSE)
  }

  rhs_m <- paste(c(x, covariates), collapse = " + ")
  rhs_y <- paste(c(x, m, covariates), collapse = " + ")
  f_m <- stats::as.formula(paste(m, "~", rhs_m))
  f_y <- stats::as.formula(paste(y, "~", rhs_y))

  fit_paths <- function(d) {
    fm <- lm(f_m, data = d)
    fy <- lm(f_y, data = d)
    c(a = unname(coef(fm)[x]), b = unname(coef(fy)[m]),
      cp = unname(coef(fy)[x]),
      se_a = sqrt(vcov(fm)[x, x]), se_b = sqrt(vcov(fy)[m, m]),
      se_cp = sqrt(vcov(fy)[x, x]), df_y = fy$df.residual)
  }
  est <- fit_paths(dat)
  a <- est["a"]; b <- est["b"]; cp <- est["cp"]
  indirect <- a * b
  total    <- cp + indirect

  # Sobel (first-order delta method) SE, via the shared variance utility.
  se_ind <- sqrt(var_indirect_effect(a = unname(a), b = unname(b),
                                     var_a = unname(est["se_a"])^2,
                                     var_b = unname(est["se_b"])^2)$value[1])

  alpha <- 1 - conf_level
  probs <- c(alpha / 2, 1 - alpha / 2)

  # Local RNG: seed if asked, and always restore the caller's state.
  if (!is.null(seed)) {
    has_old <- exists(".Random.seed", envir = globalenv())
    old <- if (has_old) get(".Random.seed", envir = globalenv()) else NULL
    on.exit({
      if (has_old) assign(".Random.seed", old, envir = globalenv())
      else if (exists(".Random.seed", envir = globalenv()))
        rm(".Random.seed", envir = globalenv())
    }, add = TRUE)
    set.seed(seed)
  }

  if (ci_method %in% c("boot_percentile", "boot_bca")) {
    boots <- vapply(seq_len(B), function(i) {
      idx <- sample.int(N, N, replace = TRUE)
      p <- fit_paths(dat[idx, , drop = FALSE])
      p["a"] * p["b"]
    }, numeric(1))
    # A degenerate resample (for example, a rank-deficient refit) returns
    # NA and would otherwise error the quantile call or bias the BCa
    # z0; drop those replications and compute on the ones that returned
    # a value, the package-wide failed-replicate policy.
    n_bad <- sum(!is.finite(boots))
    if (n_bad > 0L) {
      boots <- boots[is.finite(boots)]
      if (length(boots) < 100L) {
        stop("Only ", length(boots), " of ", B, " bootstrap replications ",
             "returned a finite indirect effect; the interval would not ",
             "be trustworthy. Check the model and the data.",
             call. = FALSE)
      }
      warning(n_bad, " of ", B, " bootstrap replications did not return ",
              "a finite indirect effect and were dropped; the interval ",
              "is computed from the ", length(boots), " that did.",
              call. = FALSE)
    }
    if (ci_method == "boot_percentile") {
      lims_ind <- stats::quantile(boots, probs, names = FALSE)
    } else {
      # BCa: bias correction from the bootstrap distribution, acceleration
      # from the jackknife (Efron, 1987).
      z0 <- qnorm(mean(boots < indirect))
      jack <- vapply(seq_len(N), function(i) {
        p <- fit_paths(dat[-i, , drop = FALSE])
        p["a"] * p["b"]
      }, numeric(1))
      jack <- jack[is.finite(jack)]
      jm <- mean(jack)
      acc <- sum((jm - jack)^3) / (6 * (sum((jm - jack)^2))^1.5)
      zq <- qnorm(probs)
      adj <- stats::pnorm(z0 + (z0 + zq) / (1 - acc * (z0 + zq)))
      lims_ind <- stats::quantile(boots, adj, names = FALSE)
    }
  } else if (ci_method == "monte_carlo") {
    a_sim <- rnorm(B, a, est["se_a"])
    b_sim <- rnorm(B, b, est["se_b"])
    lims_ind <- stats::quantile(a_sim * b_sim, probs, names = FALSE)
  } else {
    lims_ind <- indirect + qnorm(probs) * se_ind
  }

  t_crit  <- stats::qt(1 - alpha / 2, df = est["df_y"])
  lims_cp <- cp + c(-1, 1) * t_crit * est["se_cp"]
  # Total effect from its own regression for the exact t interval.
  f_t   <- stats::as.formula(paste(y, "~", rhs_m))
  fit_t <- lm(f_t, data = dat)
  se_c  <- sqrt(vcov(fit_t)[x, x])
  lims_c <- unname(coef(fit_t)[x]) +
    c(-1, 1) * stats::qt(1 - alpha / 2, fit_t$df.residual) * se_c

  prop_med <- if (abs(total) < 1e-10) NA_real_ else unname(indirect / total)

  out <- data.frame(
    term  = c("indirect_effect", "indirect_lower", "indirect_upper",
              "direct_effect", "direct_lower", "direct_upper",
              "total_effect", "total_lower", "total_upper",
              "a", "b", "se_a", "se_b", "se_indirect",
              "proportion_mediated", "N", "B"),
    value = unname(c(indirect, lims_ind[1], lims_ind[2],
                     cp, lims_cp[1], lims_cp[2],
                     total, lims_c[1], lims_c[2],
                     a, b, est["se_a"], est["se_b"], se_ind,
                     prop_med, N,
                     if (ci_method == "sobel") NA_real_ else B)),
    stringsAsFactors = FALSE
  )
  out <- .as_dmar_tbl(out, conf_level = conf_level)
  attr(out, "ci_method") <- ci_method
  out
}
