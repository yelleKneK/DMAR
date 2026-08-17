# Generalized eta squared (effect size for ANOVA).
#' Generalized Eta Squared (Effect Size for ANOVA, Comparable Across Designs)
#'
#' Computes the sample generalized eta squared (\eqn{\eta^2_G}; Olejnik &
#' Algina, 2003; Bakeman, 2005), the proportion of variance in the dependent
#' variable accounted for by a fixed effect after the variance attributable to
#' \emph{measured} (observed) factors, but not other \emph{manipulated}
#' factors, has been left in the denominator. This makes \eqn{\eta^2_G}
#' comparable across designs that differ in which factors are present, which
#' regular \eqn{\eta^2} and partial \eqn{\eta^2} are not.
#'
#' The function accepts \emph{one} of three input interfaces:
#' \enumerate{
#'   \item a fitted model object (\code{\link[stats]{aov}},
#'     \code{\link[stats]{lm}}, or \code{aovlist} for within-subjects /
#'     mixed designs) together with an \code{observed} vector listing
#'     which factors are measured (rather than manipulated);
#'   \item raw sums of squares: \code{SS_effect}, \code{SS_observed} (one value
#'     per measured factor, or a scalar), and \code{SS_error}; or
#'   \item raw \emph{F}-values and degrees of freedom: \code{F_effect},
#'     \code{df_effect}, \code{F_observed} (vector aligned with
#'     \code{df_observed}), and \code{df_error}.
#' }
#'
#' If the user supplies both the SS interface (option 2) and the F/df interface
#' (option 3), the function computes \eqn{\eta^2_G} from each and compares the
#' results to within a 1e-6 tolerance. When the two interfaces agree, the SS
#' value is returned. When they disagree, the function stops with a detailed
#' message reporting both values.
#'
#' @param object Optional. A fitted model object of class
#'   \code{\link[stats]{aov}}, \code{\link[stats]{lm}}, or \code{aovlist}
#'   (multi-stratum aov fit, e.g.\
#'   \code{aov(y ~ A + Error(subject/A), data = d)}).
#' @param observed Character vector naming the \emph{measured} (rather
#'   than \emph{manipulated}) factors. Their sums of squares, and the SS
#'   of every interaction containing a listed factor, are kept in the
#'   denominator of \eqn{\eta^2_G} per Olejnik and Algina (2003,
#'   Eq. 5). Manipulated effects are excluded from the denominator when
#'   they are not the focal effect. For \code{aovlist} fits the names must match the
#'   effect labels visible in \code{summary(object)}.
#' @param SS_effect Sum of squares for the focal effect (option 2).
#' @param SS_observed Sums of squares for the measured factors. Scalar or
#'   numeric vector (option 2).
#' @param SS_error Error (residual) sum of squares (option 2).
#' @param F_effect Observed \emph{F}-value for the focal effect (option 3).
#' @param df_effect Numerator degrees of freedom for the focal effect (option 3).
#' @param F_observed Vector of \emph{F}-values for the measured factors (option 3).
#' @param df_observed Numerator degrees of freedom for the measured factors,
#'   aligned with \code{F_observed} (option 3).
#' @param df_error Error degrees of freedom (option 3).
#'
#' @return A \code{data.frame} with one row per focal effect. With a
#'   single-stratum fit and the raw interfaces the columns are
#'   \code{effect} and \code{eta_squared_generalized}; \code{effect} is
#'   \code{"overall"} for the raw interfaces. With an \code{aovlist} fit a
#'   \code{stratum} column is added, identifying which error stratum each
#'   effect came from.
#'
#' @details
#' \strong{Formula.}
#' \deqn{\hat{\eta}^2_G = \frac{\mathit{SS}_{\text{effect}}}{\mathit{SS}_{\text{effect}} +
#'   \sum_\text{obs} \mathit{SS}_{\text{measured}} + \mathit{SS}_{\text{error}}}.}
#' For the F/df interface the equivalent ratio form is used, dividing through
#' by \eqn{\mathit{SS}_{\text{error}}} so that no total-N argument is required:
#' \eqn{\mathit{SS}_i / \mathit{SS}_{\text{error}} = F_i \cdot df_i / df_{\text{error}}}.
#'
#' \strong{Designs supported.}
#' \itemize{
#'   \item \emph{Between-subjects ANOVA (single stratum).} The function
#'     reads \code{anova(object)}. For each focal effect, the denominator
#'     is
#'     \eqn{\mathit{SS}_{\text{focal}} + \sum \mathit{SS}_{\text{measured (others)}} + \mathit{SS}_{\text{error}}}.
#'     Manipulated factors that are not the focal effect contribute
#'     nothing to the denominator.
#'   \item \emph{Within-subjects and mixed ANOVA (\code{aovlist}, multi-stratum).}
#'     The function reads \code{summary(object)} and walks every error
#'     stratum. For each focal effect, the denominator is
#'     \eqn{\mathit{SS}_{\text{focal}} + \sum \mathit{SS}_{\text{measured (others)}} + \sum_{s} \mathit{SS}_{\text{error}(s)}},
#'     where the last sum runs over every error stratum (both the
#'     between-subjects "subjects" stratum and any within-subjects error
#'     strata). This is the Olejnik & Algina (2003) / Bakeman (2005) rule
#'     that makes the subject-level variance act as an "always-measured"
#'     contributor in repeated measures designs.
#' }
#'
#' \strong{Focal-effect self-exclusion.} If the focal effect itself is
#' listed in \code{observed}, the function excludes it from the
#' observed-sum component (the focal effect's \emph{own} SS already
#' appears in the numerator and the leading term of the denominator).
#' Practically this means listing every effect as \code{observed} reduces
#' to total \eqn{\eta^2} for between-subjects designs.
#'
#' \strong{Higher-order interactions.} An effect is a measured source of
#' variance when \emph{any} factor in its term is measured (Olejnik &
#' Algina, 2003, Eq. 5; Bakeman, 2005), so listing a factor in
#' \code{observed} also places every interaction containing that factor
#' in the denominator automatically. In a design with manipulated
#' \code{A} and measured \code{c}, \code{observed = "c"} therefore puts
#' \code{c} and \code{A:c} in the denominator, which is what the cited
#' papers' worked examples do. An explicit interaction label in
#' \code{observed} is honored as given, declaring that one term measured
#' without marking its constituent factors.
#'
#' \strong{Covariates.} Under Olejnik and Algina's Eq. 5, a covariate is
#' a measured source whose SS always belongs in the denominator, so in an
#' ANCOVA list the covariate in \code{observed}. Because
#' \code{anova()} on an \code{lm}/\code{aov} fit uses sequential sums
#' of squares, enter the covariate before the treatment factors in the
#' model formula so its SS is adjusted the way the ANCOVA decomposition
#' intends.
#'
#' \strong{Confidence intervals.} See \code{\link{ci_eta_squared_generalized}}
#' for the corresponding CI function; both available CI methods are approximate
#' and require independent evaluation.
#'
#' @references
#' Bakeman, R. (2005). Recommended effect size statistics for repeated measures
#'   designs. \emph{Behavior Research Methods, 37}(3), 379--384.
#'   \doi{10.3758/BF03192707}
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
#' @seealso \code{\link{ci_eta_squared_generalized}}, \code{\link{eta_squared}},
#'   \code{\link{eta_squared_partial}}, \code{\link{ci_omega_squared}}
#'
#' @examples
#' # 1. Fitted model with `observed`. In the pygmalion expectancy
#' #        experiment, treatment is manipulated while grade is a measured
#' #        classification, so grade's variance stays in the denominator.
#' pyg <- pygmalion
#' pyg$grade <- factor(pyg$grade)
#' fit <- aov(iq_8 ~ treatment * grade, data = pyg)
#' eta_squared_generalized(fit, observed = "grade")
#'
#' # 2. Raw sums of squares.
#' eta_squared_generalized(SS_effect = 100, SS_observed = c(40, 30),
#'                         SS_error = 200)
#'
#' # 3. Raw F-values and degrees of freedom.
#' eta_squared_generalized(F_effect = 6.0, df_effect = 2,
#'                         F_observed = c(2.5, 1.8),
#'                         df_observed = c(1, 2), df_error = 50)
#'
#' # 4. Within-subjects (repeated measures) ANOVA. The denominator
#' #        automatically includes the subject-level variance plus the
#' #        within-subjects error, following Bakeman (2005).
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
#' eta_squared_generalized(fit_rm)
#'
#' # 5. Mixed design with a measured between-subjects factor. Treat
#' #        'group' as observed; its SS stays in the denominator for
#' #        'time' and the 'group:time' interaction.
#' set.seed(113)
#' n_per_group <- 10
#' n <- n_per_group * 2
#' mixed_data <- data.frame(
#'   subject = factor(rep(seq_len(n), each = 3)),
#'   group   = factor(rep(c("Treatment", "Control"), each = 3 * n_per_group)),
#'   time    = factor(rep(c("Pre", "Mid", "Post"), n),
#'                    levels = c("Pre", "Mid", "Post")),
#'   y       = rnorm(n, sd = 1)[rep(seq_len(n), each = 3)] +
#'             0.5 * rep(1:3, n) + rnorm(n * 3, sd = 1)
#' )
#' fit_mixed <- aov(y ~ group * time + Error(subject/time), data = mixed_data)
#' eta_squared_generalized(fit_mixed, observed = "group")
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @keywords htest design
#'
#' @family effect size estimates
#'
#' @export
#' @import stats

eta_squared_generalized <- function(
  object       = NULL,
  observed     = NULL,
  SS_effect    = NULL,
  SS_observed  = NULL,
  SS_error     = NULL,
  F_effect     = NULL,
  df_effect    = NULL,
  F_observed   = NULL,
  df_observed  = NULL,
  df_error     = NULL
) {
  has_a <- !is.null(SS_effect)  || !is.null(SS_observed) || !is.null(SS_error)
  has_b <- !is.null(F_effect)   || !is.null(df_effect)   ||
           !is.null(F_observed) || !is.null(df_observed) ||
           !is.null(df_error)

  if (!is.null(object)) {
    if (has_a || has_b) {
      warning("'object' supplied; raw SS or F/df arguments are ignored.")
    }
    return(.eta_squared_generalized_from_model(object, observed))
  }

  if (has_a && has_b) {
    val_a <- .eta_g_from_ss(SS_effect, SS_observed, SS_error)
    val_b <- .eta_g_from_fdf(F_effect, df_effect, F_observed, df_observed, df_error)
    if (!isTRUE(all.equal(val_a, val_b, tolerance = 1e-6))) {
      stop(sprintf(
        paste0(
          "Inconsistent inputs: the SS-based interface yields ",
          "eta_squared_generalized = %.6f, but the F/df-based interface ",
          "yields %.6f. Provide one interface or supply mutually consistent values."
        ),
        val_a, val_b
      ))
    }
    val <- val_a   # equivalent within tolerance; SS preferred as the direct form
  } else if (has_a) {
    if (is.null(SS_effect) || is.null(SS_observed) || is.null(SS_error)) {
      stop("The SS interface requires all of 'SS_effect', 'SS_observed', and 'SS_error'.")
    }
    val <- .eta_g_from_ss(SS_effect, SS_observed, SS_error)
  } else if (has_b) {
    if (is.null(F_effect) || is.null(df_effect) || is.null(F_observed) ||
        is.null(df_observed) || is.null(df_error)) {
      stop("The F/df interface requires all of 'F_effect', 'df_effect', 'F_observed', 'df_observed', and 'df_error'.")
    }
    val <- .eta_g_from_fdf(F_effect, df_effect, F_observed, df_observed, df_error)
  } else {
    stop("Provide a fitted model ('object'), or the SS interface (SS_effect, SS_observed, SS_error), or the F/df interface (F_effect, df_effect, F_observed, df_observed, df_error).")
  }

  data.frame(
    effect                 = "overall",
    eta_squared_generalized = val,
    stringsAsFactors       = FALSE,
    row.names              = NULL
  )
}


# Compute eta_g from sums of squares.
.eta_g_from_ss <- function(ss_eff, ss_obs, ss_err) {
  if (!is.numeric(ss_eff) || length(ss_eff) != 1L || ss_eff < 0)
    stop("'SS_effect' must be a single non-negative number.")
  if (!is.numeric(ss_err) || length(ss_err) != 1L || ss_err <= 0)
    stop("'SS_error' must be a single positive number.")
  if (!is.numeric(ss_obs) || any(ss_obs < 0))
    stop("'SS_observed' must be a non-negative numeric scalar or vector.")
  ss_eff / (ss_eff + sum(ss_obs) + ss_err)
}


# Compute eta_g from F-values and degrees of freedom, using the SS-ratio form.
.eta_g_from_fdf <- function(F_eff, df_eff, F_obs, df_obs, df_err) {
  if (!is.numeric(F_eff) || length(F_eff) != 1L || F_eff < 0)
    stop("'F_effect' must be a single non-negative number.")
  if (!is.numeric(df_eff) || length(df_eff) != 1L || df_eff <= 0)
    stop("'df_effect' must be a single positive number.")
  if (!is.numeric(df_err) || length(df_err) != 1L || df_err <= 0)
    stop("'df_error' must be a single positive number.")
  if (!is.numeric(F_obs) || any(F_obs < 0))
    stop("'F_observed' must be a non-negative numeric scalar or vector.")
  if (!is.numeric(df_obs) || any(df_obs <= 0))
    stop("'df_observed' must be a positive numeric scalar or vector.")
  if (length(F_obs) != length(df_obs))
    stop("'F_observed' and 'df_observed' must have the same length.")

  r_eff <- F_eff * df_eff / df_err
  r_obs <- sum(F_obs * df_obs / df_err)
  r_eff / (r_eff + r_obs + 1)
}


# Per-effect generalized eta squared from a fitted model. Dispatches on
# aovlist (multi-stratum) vs single-stratum aov/lm.
# Expand the user's `observed` declaration per Olejnik & Algina (2003,
# Eq. 5): an effect is a measured source if ANY factor in its term is
# measured, so each bare factor name in `observed` pulls in every
# interaction containing that factor. Explicit interaction labels are
# honored as given (they declare that one term measured without marking
# its constituent factors).
.expand_observed <- function(observed, effect_names) {
  if (is.null(observed) || length(observed) == 0L) return(character(0))
  obs_factors <- observed[!grepl(":", observed, fixed = TRUE)]
  auto <- if (length(obs_factors) == 0L) character(0) else {
    Filter(function(e) {
      any(strsplit(e, ":", fixed = TRUE)[[1]] %in% obs_factors)
    }, effect_names)
  }
  union(intersect(observed, effect_names), auto)
}

.eta_squared_generalized_from_model <- function(object, observed) {

  if (inherits(object, "aovlist")) {
    return(.eta_squared_generalized_from_aovlist(object, observed))
  }
  if (!inherits(object, c("aov", "lm"))) {
    stop("'object' must be an aov, lm, or aovlist fit.")
  }
  tbl <- stats::anova(object)
  if (!("Sum Sq" %in% colnames(tbl))) {
    stop("anova(object) does not contain a 'Sum Sq' column; check the model class.")
  }
  resid_row <- which(rownames(tbl) == "Residuals")
  if (length(resid_row) == 0L) {
    stop("Could not find a 'Residuals' row in anova(object).")
  }
  ss_err <- tbl[resid_row, "Sum Sq"]
  effect_names <- setdiff(rownames(tbl), "Residuals")

  if (!is.null(observed)) {
    if (!is.character(observed)) {
      stop("'observed' must be a character vector of factor names.")
    }
    missing_obs <- setdiff(observed, effect_names)
    if (length(missing_obs) > 0L) {
      stop(sprintf("'observed' refers to factor(s) not in anova(object): %s",
                   paste(missing_obs, collapse = ", ")))
    }
  }
  observed <- .expand_observed(observed, effect_names)

  results <- lapply(effect_names, function(effect) {
    ss_eff <- tbl[effect, "Sum Sq"]
    # Focal effect is excluded from the observed-SS sum even if listed.
    obs_others <- setdiff(observed, effect)
    ss_obs <- if (length(obs_others) == 0L) 0 else sum(tbl[obs_others, "Sum Sq"])
    val <- ss_eff / (ss_eff + ss_obs + ss_err)
    data.frame(
      effect                 = effect,
      eta_squared_generalized = val,
      stringsAsFactors       = FALSE,
      row.names              = NULL
    )
  })
  out <- do.call(rbind, results)
  rownames(out) <- NULL
  out
}


# Within-subjects (multi-stratum) worker for generalized eta squared. Follows
# Olejnik & Algina (2003) / Bakeman (2005): the denominator includes
# SS_focal, the SS for each listed measured factor (excluding the focal
# effect itself), and the SUM of residual SS across all error strata. The
# top-stratum residual (between-subjects "subjects" variance) is included
# automatically; users specify any additional measured fixed factors via
# `observed`.
.eta_squared_generalized_from_aovlist <- function(object, observed) {
  effects_tbl <- .aovlist_effects_table(object)   # cross-stratum effect SS
  total_ss_err <- .aovlist_total_residual_ss(object)

  effect_names <- effects_tbl$effect
  ss_lookup    <- setNames(effects_tbl$SS_effect, effect_names)

  if (!is.null(observed)) {
    if (!is.character(observed)) {
      stop("'observed' must be a character vector of factor names.")
    }
    missing_obs <- setdiff(observed, effect_names)
    if (length(missing_obs) > 0L) {
      stop(sprintf("'observed' refers to effect(s) not in the aovlist: %s",
                   paste(missing_obs, collapse = ", ")))
    }
  }
  observed <- .expand_observed(observed, effect_names)

  results <- lapply(effect_names, function(effect) {
    ss_eff <- ss_lookup[[effect]]
    obs_others <- setdiff(observed, effect)
    ss_obs <- if (length(obs_others) == 0L) 0 else sum(ss_lookup[obs_others])
    val <- ss_eff / (ss_eff + ss_obs + total_ss_err)
    data.frame(
      effect                 = effect,
      eta_squared_generalized = val,
      stratum                = effects_tbl$stratum[effects_tbl$effect == effect],
      stringsAsFactors       = FALSE,
      row.names              = NULL
    )
  })
  out <- do.call(rbind, results)
  rownames(out) <- NULL
  out
}
