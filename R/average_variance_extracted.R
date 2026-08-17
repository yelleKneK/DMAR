#' Average Variance Extracted (AVE)
#'
#' The average variance extracted (AVE) is the mean proportion of indicator
#' variance a factor accounts for, a standard convergent validity summary
#' for a reflective measurement block in confirmatory factor analysis and
#' structural equation modeling. With standardized loadings \eqn{\ell_j},
#' \deqn{\mathrm{AVE} = \frac{1}{J} \sum_j \ell_j^2.}
#' The quantity itself is elementary. In a model with cross-loadings, an
#' item contributes its loading to the AVE of every factor it loads on;
#' the same shared variance then counts toward each factor's summary, so
#' compare AVE values across factors of such a model with that overlap in
#' mind. Fornell and Larcker (1981) are
#' credited for establishing it as a validity criterion: a construct shows
#' convergent validity when its AVE reaches the conventional 0.50 (the
#' construct explains at least half its indicators' variance), and the
#' Fornell-Larcker discriminant criterion compares each construct's AVE with
#' its squared correlations with the other constructs. AVE is closely
#' related to composite reliability (omega); the modern complement on the
#' discriminant side is \code{\link{htmt}}.
#'
#' @param fit Optional \pkg{lavaan} fit (for example from
#'   \code{\link{cfa_1}} or \code{lavaan::cfa}); the standardized loadings
#'   of every factor are extracted and one AVE is reported per factor.
#' @param loadings Optional numeric vector of standardized loadings for a
#'   single block, as an alternative to \code{fit}. Supply exactly one of
#'   \code{fit} and \code{loadings}. No interval can be constructed
#'   from loadings alone (their sampling variability is not carried by
#'   the numbers), so \code{ci_method = "percentile"} requires
#'   \code{fit}.
#' @param conf_level Confidence level for the bootstrap interval
#'   (default \code{0.95}); used when \code{ci_method = "percentile"}.
#' @param ci_method Interval method: \code{"none"} (the default) or
#'   \code{"percentile"}. No closed-form interval for the AVE is in
#'   common use; the percentile bootstrap is the standard route in the
#'   validity literature. The cases in the fitted data are resampled
#'   with replacement \code{B} times, the model is refit to each
#'   resample, and each factor's interval is the pair of empirical
#'   quantiles of its \code{B} AVE values (Efron & Tibshirani, 1993).
#'   Replications whose refit fails or does not converge are dropped,
#'   and the interval is computed from those that return a value; a
#'   single warning reports how many were dropped.
#' @param B Number of bootstrap replications when
#'   \code{ci_method = "percentile"} (default \code{1000}). The
#'   default is smaller than the package's usual \code{10000} because
#'   every replication refits the model; raise it for a reported
#'   analysis when time allows.
#' @param seed Optional integer seed for the bootstrap. The default
#'   \code{NULL} uses the current state of the random number generator;
#'   a supplied seed is set internally and the prior state restored on
#'   exit.
#'
#' @return A \code{data.frame} (class \code{dmar_tbl}) with one row
#'   per factor: \code{factor} (label), \code{ave}, and
#'   \code{ci_lower} / \code{ci_upper} (the percentile bootstrap
#'   limits; \code{NA} when \code{ci_method = "none"}).
#'
#' @references
#' Efron, B., & Tibshirani, R. J. (1993). \emph{An introduction to the
#'   bootstrap}. New York, NY: Chapman & Hall/CRC.
#'
#' Fornell, C., & Larcker, D. F. (1981). Evaluating structural equation
#'   models with unobservable variables and measurement error.
#'   \emph{Journal of Marketing Research, 18}(1), 39--50.
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso \code{\link{htmt}} for discriminant validity;
#'   \code{\link{reliability_omega}} for the composite reliability of the
#'   same block (coefficient omega is what the composite-reliability
#'   literature computes); \code{\link{cfa_1}} to obtain the fit.
#'
#' @family multivariate and latent variable methods
#'
#' @keywords multivariate
#'
#' @examples
#' # Directly from the standardized loadings a paper reports:
#' average_variance_extracted(loadings = c(.8, .7, .6))
#'
#' # From a fitted model, one AVE per factor (requires lavaan).
#' data(holzinger_swineford)
#' fit <- lavaan::cfa(
#'   "verbal    =~ t6_paragraph_comprehension + t7_sentence +
#'                 t9_word_meaning
#'    deduction =~ t20_deduction + t22_problem_reasoning +
#'                 t23_series_completion",
#'   data = holzinger_swineford)
#' ave_tbl <- average_variance_extracted(fit)
#' ave_tbl
#'
#' # The Fornell and Larcker (1981) discriminant criterion compares each
#' # factor's AVE with the squared correlation between the factors: a
#' # factor should account for more of its own indicators' variance than
#' # it shares with the other factor. Here the comparison favors verbal
#' # and goes against deduction, whose AVE falls below the shared
#' # variance. Fitting with cfa_k(..., output = "measurement") puts the
#' # AVE values and the latent correlations in one table.
#' lavaan::lavInspect(fit, "cor.lv")["verbal", "deduction"]^2
#'
#' # An interval comes from ci_method = "percentile", which resamples the
#' # cases and refits the model once per replication. That refitting is
#' # why it is not run here; the call is
#' #   average_variance_extracted(fit, ci_method = "percentile",
#' #                              B = 1000, seed = 113)
#' # and a reported interval deserves the default B = 1000 or more.
#'
#' # The broom verbs: one row per factor.
#' generics::tidy(ave_tbl)
#' generics::glance(ave_tbl)
#'
#' @export
average_variance_extracted <- function(fit = NULL, loadings = NULL,
                                       conf_level = 0.95,
                                       ci_method = c("none", "percentile"),
                                       B = 1000L, seed = NULL) {
  ci_method <- match.arg(ci_method)
  if (is.null(fit) == is.null(loadings)) {
    stop("Supply exactly one of 'fit' or 'loadings'.", call. = FALSE)
  }
  if (!is.numeric(conf_level) || length(conf_level) != 1L ||
      conf_level <= 0 || conf_level >= 1) {
    stop("'conf_level' must be a single number in (0, 1).", call. = FALSE)
  }
  if (!is.null(loadings)) {
    if (ci_method != "none") {
      stop("No interval can be constructed from 'loadings' alone; their ",
           "sampling variability is not carried by the numbers. Supply ",
           "the fitted model as 'fit' to bootstrap.", call. = FALSE)
    }
    if (!is.numeric(loadings) || length(loadings) < 2L ||
        anyNA(loadings) || any(abs(loadings) > 1)) {
      stop("'loadings' must be two or more standardized loadings in ",
           "[-1, 1].", call. = FALSE)
    }
    out <- data.frame(factor = "f", ave = mean(loadings^2),
                      ci_lower = NA_real_, ci_upper = NA_real_,
                      stringsAsFactors = FALSE)
    return(.as_dmar_tbl(out))
  }
  if (!requireNamespace("lavaan", quietly = TRUE)) {
    stop("Package 'lavaan' is required when supplying a fit. Install it ",
         "with install.packages(\"lavaan\").", call. = FALSE)
  }
  if (!inherits(fit, "lavaan")) {
    stop("'fit' must be a lavaan fit object.", call. = FALSE)
  }
  ave_from <- function(f) {
    std <- lavaan::standardizedSolution(f)
    lam <- std[std$op == "=~", c("lhs", "est.std")]
    if (nrow(lam) == 0L) return(NULL)
    tapply(lam$est.std^2, lam$lhs, mean)
  }
  ave <- ave_from(fit)
  if (is.null(ave)) {
    stop("The fit contains no measurement (=~) loadings.", call. = FALSE)
  }
  out <- data.frame(factor = names(ave), ave = as.numeric(ave),
                    ci_lower = NA_real_, ci_upper = NA_real_,
                    stringsAsFactors = FALSE, row.names = NULL)

  if (ci_method == "percentile") {
    if (!is.numeric(B) || length(B) != 1L || B < 100) {
      stop("'B' must be a single integer of at least 100 when ",
           "bootstrapping.", call. = FALSE)
    }
    if (lavaan::lavInspect(fit, "ngroups") > 1L) {
      stop("The bootstrap interval is implemented for single-group ",
           "fits; refit per group or bootstrap by hand for a ",
           "multiple-group model.", call. = FALSE)
    }
    if (!is.null(seed)) {
      has_old <- exists(".Random.seed", envir = globalenv())
      old_seed <- if (has_old) get(".Random.seed", envir = globalenv())
      on.exit({
        if (has_old) assign(".Random.seed", old_seed, envir = globalenv())
        else if (exists(".Random.seed", envir = globalenv()))
          rm(".Random.seed", envir = globalenv())
      }, add = TRUE)
      set.seed(seed)
    }
    fac <- names(ave)
    boot_fun <- function(f) {
      a <- ave_from(f)
      if (is.null(a)) rep(NA_real_, length(fac)) else as.numeric(a[fac])
    }
    boots <- try(suppressWarnings(
      lavaan::bootstrapLavaan(fit, R = as.integer(B), FUN = boot_fun)),
      silent = TRUE)
    if (inherits(boots, "try-error")) {
      stop("The bootstrap failed: ",
           conditionMessage(attr(boots, "condition")),
           " (the fit must carry its data; refit with the raw data ",
           "available).", call. = FALSE)
    }
    boots <- matrix(as.numeric(boots), ncol = length(fac))
    ok <- rowSums(!is.finite(boots)) == 0L
    n_bad <- as.integer(B) - sum(ok)
    if (n_bad > 0L) {
      if (sum(ok) < 100L) {
        stop("Only ", sum(ok), " of ", B, " bootstrap replications ",
             "refit successfully; the interval would not be ",
             "trustworthy.", call. = FALSE)
      }
      warning(n_bad, " of ", B, " bootstrap replications failed to ",
              "refit (or did not converge) and were dropped; the ",
              "interval is computed from the ", sum(ok), " that did.",
              call. = FALSE)
      boots <- boots[ok, , drop = FALSE]
    }
    alpha_2 <- (1 - conf_level) / 2
    out$ci_lower <- apply(boots, 2, stats::quantile, probs = alpha_2,
                          names = FALSE)
    out$ci_upper <- apply(boots, 2, stats::quantile, probs = 1 - alpha_2,
                          names = FALSE)
  }
  out <- .as_dmar_tbl(out)
  if (ci_method == "percentile") {
    attr(out, "conf_level") <- conf_level
    attr(out, "B_used") <- nrow(boots)
  }
  out
}
