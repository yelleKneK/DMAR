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
#'   mediator, and the outcome columns in \code{data}. Each column must be
#'   numeric (a Date, date-time, or time difference column is read as its
#'   number, and a one-column matrix, such as \code{scale()} returns, as
#'   its column); a predictor with two categories is supplied coded 0
#'   and 1.
#'   Any column name is accepted as it is.
#' @param covariates Optional character vector of covariate column names,
#'   entered in both the mediator and the outcome models. Numeric, logical,
#'   factor, and character columns are accepted; a factor or character
#'   covariate is coded once, on the complete cases, with R's contrasts for
#'   the session. A covariate with a single value in the complete cases
#'   adds nothing to either model and is left out with a message. A column
#'   may fill one role only.
#' @param ci_method Confidence interval method for the indirect effect:
#'   \code{"boot_percentile"} (default), \code{"boot_bca"},
#'   \code{"monte_carlo"} (simulate \eqn{a} and \eqn{b} from their joint
#'   normal approximation; MacKinnon, Lockwood, & Williams, 2004), or
#'   \code{"sobel"} (the first-order normal-theory interval; reported for
#'   comparison, not recommended for inference). Direct and total effects
#'   always carry their ordinary \emph{t}-based intervals.
#' @param B Number of bootstrap or Monte Carlo replications, a single
#'   whole number. Defaults to 2000, because each bootstrap replication
#'   refits both regressions. The fewest accepted depends on the interval
#'   and the level: 400 for the percentile and Monte Carlo intervals at
#'   the default \code{conf_level} (200 at 0.90, 2000 at 0.99), and 2000
#'   for the BCa interval (1000 at 0.90, 10000 at 0.99, which the default
#'   does not meet). A smaller \code{B} is refused before anything is
#'   computed (class \code{dmar_count_invalid}); see
#'   \emph{Bootstrap replications} under \emph{Details} and
#'   \code{\link{dmar_bootstrap}}. Not used, and not checked, for
#'   \code{ci_method = "sobel"}.
#' @param conf_level Confidence level. Defaults to 0.95.
#' @param seed Optional integer seed for the resampling, used locally (the
#'   caller's random number generator state is restored on exit). Default
#'   \code{NULL} draws the resamples from the current random number
#'   stream, which advances that stream.
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
#' interval does from a \code{B} well above the default.
#'
#' \strong{Bootstrap replications.} Each resample refits the two
#' regressions by least squares on the rows drawn, with the covariates
#' coded as they are in the data, and the interval is read from the
#' quantiles of the \code{B} values of \eqn{a b}. A limit read from a
#' quantile is located by the replications that lie beyond it, so the
#' percentile and Monte Carlo intervals are run only when \code{B} puts
#' at least 10 replications beyond each limit (400 at the 95 percent
#' level), and the BCa interval, whose limits move farther into the
#' tails, only with 50 beyond each nominal limit (2000 at the 95 percent
#' level). A resample that leaves a covariate column redundant (a
#' category with no case in the resample) still gives unique values of
#' \eqn{a} and \eqn{b} and is kept; the \code{"n_covariate_redundant"}
#' attribute counts these. A resample in which \eqn{a} or \eqn{b} is not
#' estimable (the predictor constant, or the predictor or the mediator a
#' linear combination of the other predictors there) has no indirect
#' effect. Such a replication is dropped and counted, never replaced, and
#' the interval is computed from the rest. The result records \code{B},
#' \code{B_used}, \code{n_dropped}, and \code{B_required} as attributes
#' and \code{B} and \code{B_used} as rows. One warning of class
#' \code{dmar_bootstrap_replications_dropped} states the counts when any
#' replication was dropped, with a caution when the dropped replications
#' are at least the share of \code{B} beyond one limit (2.5 percent at
#' the 95 percent level), because they could then decide by themselves
#' where that limit falls. When more than half were dropped, no interval is
#' reported: the two limits
#' of the indirect effect are \code{NA}, every other row is as it would
#' have been, and the warning has class
#' \code{dmar_bootstrap_insufficient}. An interval from the
#' replications that remain describes resamples in which the effect is
#' defined. Dropping these replications is the rule DMAR applies to a
#' replication that returns no value; it is a convention of the package,
#' not a result that gives the interval its nominal coverage. Many
#' dropped replications mean the data carry little information about a
#' path (a predictor that varies in few cases), a design in which the
#' case bootstrap is itself not reliable.
#'
#' The BCa interval needs the indirect effect with each case left out;
#' when one of those is not estimable, or the bias correction or the
#' acceleration is not finite, the BCa interval is refused (class
#' \code{dmar_bca_undefined}) and the percentile interval remains
#' available. After the replications are run, the count of 10 is applied
#' again at the two positions where the bias correction and the
#' acceleration place the limits. When fewer than 10 replications lie
#' beyond one of them, no BCa interval is reported (limits \code{NA},
#' warning of class \code{dmar_bootstrap_insufficient}, which carries the
#' two positions as \code{positions} and a suggested \code{B} as
#' \code{B_suggest}).
#'
#' \code{options(dmar.bootstrap_exploratory = TRUE)} lets a smaller
#' \code{B} (2 or more) run for a first look or a demonstration. Every
#' such call gives a warning of class \code{dmar_bootstrap_exploratory},
#' and the attribute \code{"bootstrap_status"} is \code{"exploratory"};
#' the interval is not one to report.
#'
#' The Sobel
#' standard error is computed by
#' \code{\link{var_indirect_effect}}. The proportion mediated is one of
#' the effect size measures for mediation models surveyed by Preacher
#' and Kelley (2011); its instability when the total effect is small is
#' why it is reported as \code{NA} near a zero total effect. Listwise
#' deletion is applied to the analysis variables; for full information
#' maximum likelihood under missingness, fit the model in \pkg{lavaan}
#' (see \code{\link{mlmr}} for the package's FIML front end philosophy).
#' The model has one mediator and no interaction, the cases are treated
#' as independent, and incomplete cases are removed; for several
#' mediators or missing data see \code{\link{mediation_mbco}}.
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
#'   (\code{indirect_effect_ci_lower} / \code{indirect_effect_ci_upper},
#'   \code{direct_effect_ci_lower} / \code{direct_effect_ci_upper},
#'   \code{total_effect_ci_lower} / \code{total_effect_ci_upper}),
#'   \code{proportion_mediated} (\eqn{ab/c}; \code{NA} when the total
#'   effect is near zero, where the ratio is unstable), \code{N},
#'   \code{B}, and \code{B_used}, the number of replications the interval
#'   of the indirect effect rests on (equal to \code{B} for the Monte
#'   Carlo interval and \code{NA}, as \code{B} is, for the Sobel
#'   interval). The interval method is recorded in the \code{"ci_method"}
#'   attribute and the confidence level in \code{"conf_level"}. For the
#'   bootstrap and Monte Carlo intervals the attributes \code{"B"}
#'   (requested), \code{"B_used"} (the number the interval rests on),
#'   \code{"n_dropped"} (their difference), and \code{"B_required"} (the
#'   fewest replications a request may ask for; that count is held to
#'   \code{B}, not to \code{B_used}) hold the
#'   counts as integers, which is where \code{\link{dmar_metadata}} reads
#'   them, and \code{"bootstrap_status"} is one of \code{"complete"},
#'   \code{"dropped"}, \code{"exploratory"}, \code{"insufficient_count"},
#'   or \code{"insufficient_share"} (the last two mean that the limits of
#'   the indirect effect are \code{NA}). For the two bootstrap intervals
#'   \code{"n_covariate_redundant"} is the number of resamples kept in
#'   which a covariate column was redundant (see \emph{Details}).
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
#' # The percentile bootstrap is the default and is what a reported
#' # analysis would use. Each resample refits both regressions.
#' # B = 400 is the fewest replications a 95 percent percentile interval
#' # is reported from and keeps the example quick; a reported interval
#' # deserves the default B = 2000 or more.
#' mediate(d, x = "x", m = "m", y = "y", B = 400, seed = 113)
#'
#' # The Monte Carlo interval draws a and b from their joint normal
#' # approximation instead of refitting, so it stays quick even at a
#' # large B. Its limits sit close to the bootstrap limits here.
#' mediate(d, x = "x", m = "m", y = "y", ci_method = "monte_carlo",
#'         B = 10000, seed = 113)
#'
#' # The Sobel interval is closed form. It assumes the product ab is
#' # normally distributed, which is why it is reported for comparison
#' # rather than used for inference; its B and B_used rows are NA
#' # because no replications are drawn.
#' mediate(d, x = "x", m = "m", y = "y", ci_method = "sobel")
#'
#' @export
#' @importFrom stats complete.cases qnorm quantile rnorm
mediate <- function(data, x, m, y, covariates = NULL,
                    ci_method = c("boot_percentile", "boot_bca",
                                  "monte_carlo", "sobel"),
                    B = 2000, conf_level = 0.95, seed = NULL) {
  ci_method <- match.arg(ci_method)
  if (!is.data.frame(data)) {
    .dmar_refuse_argument("'data' must be a data.frame.", arg = "data")
  }
  for (nm in c("x", "m", "y")) {
    val <- get(nm)
    if (!is.character(val) || length(val) != 1L || is.na(val) ||
        !(val %in% names(data))) {
      .dmar_refuse_argument(
        sprintf("'%s' must name one column of 'data'.", nm), arg = nm)
    }
  }
  if (!is.null(covariates) &&
      (!is.character(covariates) || anyNA(covariates) ||
       !all(covariates %in% names(data)))) {
    .dmar_refuse_argument(
      "'covariates' must name columns of 'data'.",
      arg = c("covariates", "data"))
  }
  vars <- c(x, m, y, covariates)
  if (anyDuplicated(vars)) {
    .dmar_refuse_argument(c(
      "'x', 'm', 'y', and 'covariates' must name different columns; '",
      vars[duplicated(vars)][1L], "' is named more than once."),
      arg = c("x", "m", "y"))
  }
  # The predictor, the mediator, and the outcome enter as numbers. Dates,
  # date-times, and time differences are numbers with a class and are read
  # with as.numeric(); a factor, a logical, or a character column is not
  # recoded here, because which category is the 1 is the analyst's choice.
  # A one-column matrix or a one-dimensional array (what scale() and table
  # arithmetic leave in a data frame) is one variable and is accepted; a
  # matrix of two or more columns is not.
  for (nm in c("x", "m", "y")) {
    col <- data[[get(nm)]]
    if (!(is.numeric(col) ||
          inherits(col, c("Date", "POSIXct", "difftime"))) ||
        is.factor(col) || NCOL(col) != 1L || length(dim(col)) > 2L) {
      .dmar_refuse_argument(
        sprintf(paste0("'%s' must name a numeric column. A predictor ",
                       "with two categories can be coded 0 and 1; a ",
                       "categorical mediator or outcome is outside the ",
                       "linear model this function fits."), nm), arg = nm)
    }
  }
  if (!is.numeric(conf_level) || length(conf_level) != 1L ||
      is.na(conf_level) || conf_level <= 0 || conf_level >= 1) {
    .dmar_refuse_argument(
      "'conf_level' must be a single number in (0, 1).", arg = "conf_level")
  }
  # The number of replications is checked against the rule of
  # ?dmar_bootstrap, and only when replications are requested: 10 beyond
  # each limit for the percentile and Monte Carlo intervals, 50 beyond
  # each nominal limit for BCa. The Sobel interval draws nothing.
  boot_kind <- if (ci_method == "boot_bca") "bca" else "quantile"
  if (ci_method != "sobel") {
    B <- .dmar_check_B(B, kind = boot_kind, tail = (1 - conf_level) / 2,
                       conf_level = conf_level)
  }

  dat  <- data[stats::complete.cases(data[, vars, drop = FALSE]),
               vars, drop = FALSE]
  N <- nrow(dat)
  if (N < length(vars) + 5L) {
    .dmar_refuse_argument(
      "Too few complete cases to fit the mediation model.", arg = NULL)
  }
  x_v <- as.numeric(dat[[x]])
  m_v <- as.numeric(dat[[m]])
  y_v <- as.numeric(dat[[y]])

  # The covariates are coded once, on the complete cases, under internal
  # names. No column name is then read as formula syntax, and every
  # resample is fit to the same columns: a category with no case in a
  # resample is a constant column there, not a different model.
  C <- matrix(0, N, 0L)
  if (length(covariates)) {
    cd <- droplevels(as.data.frame(dat[, covariates, drop = FALSE]))
    names(cd) <- paste0("covariate_", seq_along(covariates))
    one_value <- vapply(cd, function(col) length(unique(col)) < 2L,
                        logical(1))
    if (any(one_value)) {
      several <- sum(one_value) > 1L
      message(if (several) "Covariates " else "Covariate ",
              paste0("'", covariates[one_value], "'", collapse = ", "),
              if (several) " have" else " has",
              " a single value in the complete cases, ",
              if (several) "add" else "adds",
              " nothing to the model, and ",
              if (several) "are" else "is", " left out.")
      cd <- cd[, !one_value, drop = FALSE]
    }
    if (ncol(cd)) {
      C <- stats::model.matrix(~ ., data = cd)[, -1L, drop = FALSE]
    }
  }
  if (!all(is.finite(x_v)) || !all(is.finite(m_v)) ||
      !all(is.finite(y_v)) || !all(is.finite(C))) {
    .dmar_refuse_argument(
      "The analysis variables must be finite (no Inf or -Inf).", arg = NULL)
  }
  # Centering changes no slope (every design has the intercept) and makes
  # the estimability decision below independent of the origin of the
  # variables, as it already is of their units.
  x_c <- x_v - mean(x_v)
  m_c <- m_v - mean(m_v)
  C <- sweep(C, 2L, colMeans(C), "-")
  one <- rep(1, N)
  # Each design puts the coefficient it is read for in its last column.
  X_a  <- cbind(one, C, x_c)            # m on x: a
  X_b  <- cbind(one, C, x_c, m_c)       # y on x and m: b
  X_cp <- cbind(one, C, m_c, x_c)       # y on x and m: c'
  fit_a  <- .mediate_focal_fit(X_a, m_v)
  fit_b  <- .mediate_focal_fit(X_b, y_v)
  fit_cp <- .mediate_focal_fit(X_cp, y_v)
  fit_c  <- .mediate_focal_fit(X_a, y_v)  # total effect, its own regression
  if (anyNA(c(fit_a[["est"]], fit_b[["est"]], fit_cp[["est"]]))) {
    .dmar_refuse_argument(c(
      "The paths of the mediation model are not estimable in these ",
      "data: ",
      if (is.na(fit_a[["est"]])) "'x' is constant or a linear " else
        "'m' or 'x' is a linear ",
      "combination of the other predictors."), arg = NULL)
  }
  if (fit_b[["df"]] < 1) {
    .dmar_refuse_argument(c(
      "Too few complete cases to fit the mediation model: the outcome ",
      "regression has no residual degrees of freedom."), arg = NULL)
  }
  a <- fit_a[["est"]]; b <- fit_b[["est"]]; cp <- fit_cp[["est"]]
  se_a <- fit_a[["se"]]; se_b <- fit_b[["se"]]; se_cp <- fit_cp[["se"]]
  indirect <- a * b
  total    <- cp + indirect

  # Sobel (first-order delta method) SE, via the shared variance utility.
  se_ind <- sqrt(var_indirect_effect(reg_coef_a = a, reg_coef_b = b,
                                     var_reg_coef_a = se_a^2,
                                     var_reg_coef_b = se_b^2)$value[1])

  alpha <- 1 - conf_level
  probs <- c(alpha / 2, 1 - alpha / 2)

  # A supplied seed is local to this call; the caller's random number
  # generator state is restored on exit.
  .dmar_local_seed(seed)

  B_used <- NA_integer_
  n_bad <- NA_integer_
  n_redundant <- NA_integer_
  if (ci_method %in% c("boot_percentile", "boot_bca")) {
    # The indirect effect on the rows 'idx', NA when a or b is not
    # estimable there, with the rank of the outcome design on those rows.
    rank_b <- N - fit_b[["df"]]
    ab_of <- function(idx) {
      p_a <- .mediate_focal_coef(X_a[idx, , drop = FALSE], m_v[idx])
      p_b <- .mediate_focal_coef(X_b[idx, , drop = FALSE], y_v[idx])
      c(p_a[[1L]] * p_b[[1L]], p_b[[2L]])
    }
    draws <- vapply(seq_len(B), function(i) {
      ab_of(sample.int(N, N, replace = TRUE))
    }, numeric(2))
    # A resample in which a or b is not estimable (x constant, or x or m a
    # linear combination of the other predictors there) has no indirect
    # effect. Those replications are dropped and counted, in the warning,
    # in the B_used row, and in the attributes: the package-wide rule for
    # a replication that returns no value. A resample that loses rank
    # among the covariate columns only (the rank of its outcome design is
    # below that of the data although a and b are estimable) is kept, and
    # counted on its own.
    kept <- is.finite(draws[1L, ])
    n_bad <- sum(!kept)
    boots <- draws[1L, kept]
    B_used <- length(boots)
    n_redundant <- sum(draws[2L, kept] < rank_b)
    # The verdict of the rule (R/dmar_bootstrap.R): no interval from a
    # run in which more than half of the replications were dropped (the
    # count rule is held to the request, checked above). The conditions are signaled
    # once the interval is settled, so that a BCa interval that is refused
    # or withheld below is not preceded by a warning about an interval.
    verdict <- .dmar_boot_settle(
      B, B_used, kind = boot_kind, tail = alpha / 2,
      conf_level = conf_level, fun = "mediate",
      quantity = "the indirect effect",
      dropped_because = paste0(
        "left a path of the indirect effect not estimable (the predictor ",
        "or the mediator constant, or a linear combination of the other ",
        "predictors, in the resample)"),
      signal = FALSE)
    lims_ind <- c(NA_real_, NA_real_)
    if (verdict$formed && ci_method == "boot_percentile") {
      lims_ind <- stats::quantile(boots, probs, names = FALSE)
    } else if (verdict$formed) {
      # BCa: bias correction from the bootstrap distribution, acceleration
      # from the jackknife with the leave-one-out values centered at their
      # mean (Efron & Tibshirani, 1993, Equations 14.10, 14.14, and 14.15).
      # A leave-one-out value that is not estimable is not removed: the
      # acceleration is then undefined and the interval is refused.
      z0 <- qnorm(mean(boots < indirect))
      jack <- vapply(seq_len(N), function(i) ab_of(seq_len(N)[-i])[1L],
                     numeric(1))
      jm <- mean(jack)
      acc <- sum((jm - jack)^3) / (6 * (sum((jm - jack)^2))^1.5)
      zq <- qnorm(probs)
      adj <- stats::pnorm(z0 + (z0 + zq) / (1 - acc * (z0 + zq)))
      if (!all(is.finite(jack)) || !is.finite(z0) || !is.finite(acc) ||
          !all(is.finite(adj))) {
        .dmar_refuse(c(
          "The BCa interval is not available for these data: ",
          if (!all(is.finite(jack))) {
            paste0("its acceleration needs the indirect effect with each ",
                   "case left out in turn, and the effect is not estimable ",
                   "with case ",
                   paste(utils::head(which(!is.finite(jack)), 5L),
                         collapse = ", "),
                   " left out (row numbers among the complete cases).")
          } else if (!is.finite(z0)) {
            paste0("every bootstrap estimate falls on one side of the ",
                   "sample estimate, so the bias correction is infinite.")
          } else {
            "its acceleration is not finite."
          },
          " ci_method = \"boot_percentile\" does not need these quantities."),
          "bca_undefined")
      }
      # The count of 10 beyond each limit is applied where the bias
      # correction and the acceleration place the limits.
      verdict <- .dmar_boot_bca_positions(verdict, adj)
      if (verdict$formed) {
        lims_ind <- stats::quantile(boots, adj, names = FALSE)
      }
    }
    .dmar_boot_signal(verdict)
  } else if (ci_method == "monte_carlo") {
    a_sim <- rnorm(B, a, se_a)
    b_sim <- rnorm(B, b, se_b)
    lims_ind <- stats::quantile(a_sim * b_sim, probs, names = FALSE)
    # A product of two normal draws is always a number, so every draw is
    # used; the verdict records the counts and, under the exploratory
    # option, marks a run below the requirement.
    B_used <- as.integer(B)
    verdict <- .dmar_boot_settle(
      B, B_used, kind = "quantile", tail = alpha / 2,
      conf_level = conf_level, fun = "mediate",
      quantity = "the indirect effect",
      dropped_because = "gave a product that is not a finite number",
      replications = "Monte Carlo draws")
  } else {
    lims_ind <- indirect + qnorm(probs) * se_ind
  }

  t_crit  <- stats::qt(1 - alpha / 2, df = fit_b[["df"]])
  lims_cp <- cp + c(-1, 1) * t_crit * se_cp
  # Total effect from its own regression for the exact t interval.
  lims_c <- fit_c[["est"]] +
    c(-1, 1) * stats::qt(1 - alpha / 2, fit_c[["df"]]) * fit_c[["se"]]

  prop_med <- if (abs(total) < 1e-10) NA_real_ else indirect / total

  out <- data.frame(
    term  = c("indirect_effect", "indirect_effect_ci_lower",
              "indirect_effect_ci_upper",
              "direct_effect", "direct_effect_ci_lower",
              "direct_effect_ci_upper",
              "total_effect", "total_effect_ci_lower", "total_effect_ci_upper",
              "a", "b", "se_a", "se_b", "se_indirect",
              "proportion_mediated", "N", "B", "B_used"),
    value = c(indirect, lims_ind[1], lims_ind[2],
              cp, lims_cp[1], lims_cp[2],
              total, lims_c[1], lims_c[2],
              a, b, se_a, se_b, se_ind,
              prop_med, N,
              if (ci_method == "sobel") NA_real_ else B, B_used),
    stringsAsFactors = FALSE
  )
  out <- .as_dmar_tbl(out, conf_level = conf_level)
  attr(out, "ci_method") <- ci_method
  if (!is.na(B_used)) out <- .dmar_boot_attrs(out, verdict)
  if (!is.na(n_bad)) attr(out, "n_covariate_redundant") <- n_redundant
  out
}

# Least squares on a fixed design matrix whose last column is the focal
# one. Returns the focal coefficient, NA when it is not estimable, and the
# rank of the design. With the focal column last, the limited pivoting of
# LINPACK's dqrdc2 (the routine under lm()) moves it out of the first
# 'rank' columns exactly when it is, to the tolerance, a linear combination
# of the columns before it, which is the condition for its coefficient not
# to be estimable. A redundant column among the others (a factor level
# absent from a resample) is pivoted out and leaves the focal coefficient
# unique.
.mediate_focal_coef <- function(X, response) {
  fit <- stats::.lm.fit(X, response)
  pos <- match(ncol(X), fit$pivot)
  c(if (pos > fit$rank) NA_real_ else fit$coefficients[pos], fit$rank)
}

# The same fit with the standard error of the focal coefficient and the
# residual degrees of freedom (the number of rows less the rank, not less
# the number of columns), for the data as observed.
.mediate_focal_fit <- function(X, response) {
  fit <- stats::.lm.fit(X, response)
  r <- fit$rank
  pos <- match(ncol(X), fit$pivot)
  if (pos > r) return(c(est = NA_real_, se = NA_real_, df = NA_real_))
  df <- nrow(X) - r
  R_inv <- backsolve(fit$qr[seq_len(r), seq_len(r), drop = FALSE], diag(r))
  s2 <- sum(fit$residuals^2) / df
  c(est = fit$coefficients[pos], se = sqrt(s2 * sum(R_inv[pos, ]^2)),
    df = df)
}
