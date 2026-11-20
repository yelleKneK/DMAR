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
#'   \item raw sums of squares: \code{SS_effect}, \code{SS_observed} (a
#'     scalar, or one value per other source in the denominator; see
#'     \code{SS_observed}), and \code{SS_error}; or
#'   \item raw \emph{F}-values and degrees of freedom: \code{F_effect},
#'     \code{df_effect}, \code{F_observed} (vector aligned with
#'     \code{df_observed}), and \code{df_error}.
#' }
#'
#' The estimate from sums of squares needs no degrees of freedom, while the
#' parametric interval that \code{\link{ci_eta_squared_generalized}} builds
#' from the same sums of squares needs \code{df_effect} and
#' \code{df_error}, so the two are accepted with the sums of squares here
#' and leave the estimate unchanged. With sums of squares, a call is read
#' as the \emph{F}/df interface only when it also gives \code{F_effect},
#' \code{F_observed}, or \code{df_observed}.
#'
#' If the user supplies both the SS interface (option 2) and the \emph{F}/df
#' interface (option 3), the function computes \eqn{\eta^2_G} from each and
#' compares the results to within a 1e-6 tolerance. When the two interfaces
#' agree, the SS value is returned. When they disagree, the function stops
#' with a detailed message reporting both values.
#' \code{\link{ci_eta_squared_generalized}} refuses the two together, since
#' they can agree on the estimate and still imply different \emph{F}
#' statistics for the effect, and so different intervals.
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
#' @param SS_observed Sums of squares for the other sources in the
#'   denominator: the measured factors, the interactions that contain one,
#'   and covariates, and in a design with an \code{Error()} term the
#'   residuals of the strata other than the effect's. Scalar or numeric
#'   vector (option 2).
#' @param SS_error Error (residual) sum of squares of the focal effect: the
#'   sum of squares its \emph{F} test divides by, which in a design with an
#'   \code{Error()} term is the residual of the effect's stratum (option
#'   2). The estimate depends on \code{SS_observed} and \code{SS_error} only
#'   through their sum, but the parametric interval of
#'   \code{\link{ci_eta_squared_generalized}} is built on the \emph{F} they
#'   imply.
#' @param F_effect Observed \emph{F}-value for the focal effect (option 3).
#' @param df_effect Numerator degrees of freedom for the focal effect
#'   (option 3; accepted with the sums of squares of option 2, which do not
#'   use it, as Details explains).
#' @param F_observed Vector of \emph{F}-values for the other sources in the
#'   denominator, each the source's mean square divided by the focal
#'   effect's error mean square (option 3).
#' @param df_observed Numerator degrees of freedom for those sources,
#'   aligned with \code{F_observed} (option 3).
#' @param df_error Degrees of freedom of the focal effect's error term
#'   (option 3; accepted with the sums of squares of option 2, which do not
#'   use it, as Details explains).
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
#' For the \emph{F}/df interface the equivalent ratio form is used, dividing
#' through by \eqn{\mathit{SS}_{\text{error}}} so that no total-N argument is
#' required:
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
#'     contributor in repeated measures designs. Each effect must be
#'     tested in a single stratum, as it is when every subject is
#'     observed equally often in every within-subjects condition,
#'     whatever the sizes of the between-subjects groups, and every
#'     covariate is constant within subjects. When observations are
#'     missing within subjects, or a covariate varies within them,
#'     \code{aov()} divides the effect's sum of squares among strata and
#'     tests the effect in each of them; no one stratum's share is the
#'     effect's sum of squares, and Olejnik and Algina (2003, p. 437)
#'     define the estimate for balanced designs, so the function stops and
#'     names the effects and strata involved.
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
#'   \code{\link{ci_omega_squared}}
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
  # df_effect and df_error belong to both raw interfaces: the F/df
  # interface needs them for the estimate, and ci_eta_squared_generalized()
  # needs them beside the sums of squares for its parametric interval. Only
  # F_effect, F_observed, and df_observed mark the F/df interface when sums
  # of squares are given; DMAR 1.0.x read the degrees of freedom as that
  # interface too and stopped for want of F_effect.
  has_a <- !is.null(SS_effect)  || !is.null(SS_observed) || !is.null(SS_error)
  has_b <- !is.null(F_effect)   || !is.null(F_observed)  ||
           !is.null(df_observed) ||
           (!has_a && (!is.null(df_effect) || !is.null(df_error)))

  if (!is.null(object)) {
    if (has_a || has_b) {
      .dmar_warn(
        "'object' supplied; raw SS or F/df arguments are ignored.",
        "argument_ignored")
    }
    return(.eta_squared_generalized_from_model(object, observed))
  }

  if (has_a && has_b) {
    val_a <- .eta_g_from_ss(SS_effect, SS_observed, SS_error)
    val_b <- .eta_g_from_fdf(F_effect, df_effect, F_observed, df_observed, df_error)
    if (!isTRUE(all.equal(val_a, val_b, tolerance = 1e-6))) {
      .dmar_refuse_argument(
        sprintf(
     paste0(
       "Inconsistent inputs: the SS-based interface yields ",
       "eta_squared_generalized = %.6f, but the F/df-based interface ",
       "yields %.6f. Provide one interface or supply mutually consistent values."
     ),
     val_a, val_b
   ), arg = NULL)
    }
    val <- val_a   # equivalent within tolerance; SS preferred as the direct form
  } else if (has_a) {
    if (is.null(SS_effect) || is.null(SS_observed) || is.null(SS_error)) {
      .dmar_refuse_argument(
        "The SS interface requires all of 'SS_effect', 'SS_observed', and 'SS_error'.",
        arg = c("SS_effect", "SS_observed", "SS_error"))
    }
    val <- .eta_g_from_ss(SS_effect, SS_observed, SS_error)
  } else if (has_b) {
    if (is.null(F_effect) || is.null(df_effect) || is.null(F_observed) ||
        is.null(df_observed) || is.null(df_error)) {
      .dmar_refuse_argument(
        "The F/df interface requires all of 'F_effect', 'df_effect', 'F_observed', 'df_observed', and 'df_error'.",
        arg = c("F_effect", "df_effect", "F_observed"))
    }
    val <- .eta_g_from_fdf(F_effect, df_effect, F_observed, df_observed, df_error)
  } else {
    .dmar_refuse_argument(
      "Provide a fitted model ('object'), or the SS interface (SS_effect, SS_observed, SS_error), or the F/df interface (F_effect, df_effect, F_observed, df_observed, df_error).",
      arg = "object")
  }

  data.frame(
    effect                 = "overall",
    eta_squared_generalized = val,
    stringsAsFactors       = FALSE,
    row.names              = NULL
  )
}


# Compute eta_g from sums of squares. Every value must be finite: DMAR 1.0.x
# stopped on most missing values with R's "missing value where TRUE/FALSE
# needed" and returned 0 or NaN for an infinite one.
.eta_g_from_ss <- function(ss_eff, ss_obs, ss_err) {
  if (!is.numeric(ss_eff) || length(ss_eff) != 1L || !is.finite(ss_eff) ||
      ss_eff < 0)
    .dmar_refuse_argument(
      "'SS_effect' must be a single finite non-negative number.",
      arg = "SS_effect")
  if (!is.numeric(ss_err) || length(ss_err) != 1L || !is.finite(ss_err) ||
      ss_err <= 0)
    .dmar_refuse_argument(
      "'SS_error' must be a single finite positive number.", arg = "SS_error")
  if (!is.numeric(ss_obs) || !all(is.finite(ss_obs)) || any(ss_obs < 0))
    .dmar_refuse_argument(
      "'SS_observed' must be a finite non-negative numeric scalar or vector.",
      arg = "SS_observed")
  ss_eff / (ss_eff + sum(ss_obs) + ss_err)
}


# The degrees of freedom of the focal effect and of its error term, which
# the F/df interface needs for the estimate and the sums of squares need
# for the parametric interval of ci_eta_squared_generalized().
.eta_g_check_df <- function(df_eff, df_err) {
  if (!is.numeric(df_eff) || length(df_eff) != 1L || !is.finite(df_eff) ||
      df_eff <= 0)
    .dmar_refuse_argument(
      "'df_effect' must be a single finite positive number.", arg = "df_effect")
  if (!is.numeric(df_err) || length(df_err) != 1L || !is.finite(df_err) ||
      df_err <= 0)
    .dmar_refuse_argument(
      "'df_error' must be a single finite positive number.", arg = "df_error")
  invisible(NULL)
}


# Compute eta_g from F-values and degrees of freedom, using the SS-ratio
# form. Every value must be finite, as for the sums of squares.
.eta_g_from_fdf <- function(F_eff, df_eff, F_obs, df_obs, df_err) {
  if (!is.numeric(F_eff) || length(F_eff) != 1L || !is.finite(F_eff) ||
      F_eff < 0)
    .dmar_refuse_argument(
      "'F_effect' must be a single finite non-negative number.",
      arg = "F_effect")
  .eta_g_check_df(df_eff, df_err)
  if (!is.numeric(F_obs) || !all(is.finite(F_obs)) || any(F_obs < 0))
    .dmar_refuse_argument(
      "'F_observed' must be a finite non-negative numeric scalar or vector.",
      arg = "F_observed")
  if (!is.numeric(df_obs) || !all(is.finite(df_obs)) || any(df_obs <= 0))
    .dmar_refuse_argument(
      "'df_observed' must be a finite positive numeric scalar or vector.",
      arg = "df_observed")
  if (length(F_obs) != length(df_obs))
    .dmar_refuse_argument(
      "'F_observed' and 'df_observed' must have the same length.",
      arg = c("F_observed", "df_observed"))

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
    .dmar_refuse_argument(
      "'object' must be an aov, lm, or aovlist fit.", arg = "object")
  }
  tbl <- stats::anova(object)
  if (!("Sum Sq" %in% colnames(tbl))) {
    .dmar_refuse_argument(
      "anova(object) does not contain a 'Sum Sq' column; check the model class.",
      arg = NULL)
  }
  resid_row <- which(rownames(tbl) == "Residuals")
  if (length(resid_row) == 0L) {
    .dmar_refuse_argument(
      "Could not find a 'Residuals' row in anova(object).", arg = NULL)
  }
  ss_err <- tbl[resid_row, "Sum Sq"]
  effect_names <- setdiff(rownames(tbl), "Residuals")

  if (!is.null(observed)) {
    if (!is.character(observed)) {
      .dmar_refuse_argument(
        "'observed' must be a character vector of factor names.",
        arg = "observed")
    }
    missing_obs <- setdiff(observed, effect_names)
    if (length(missing_obs) > 0L) {
      .dmar_refuse_argument(
        sprintf("'observed' refers to factor(s) not in anova(object): %s",
                paste(missing_obs, collapse = ", ")), arg = "observed")
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

  # When observations are missing within subjects, or a covariate varies
  # within them, aov() divides the effect's sum of squares among strata
  # and tests the effect in each of them. No one stratum's share is the
  # effect's sum of squares, and Olejnik and Algina (2003, p. 437) define
  # the estimate for balanced designs, so the function stops rather than
  # choose a share. DMAR 1.0.x took the first stratum's share and
  # returned a row per stratum each time the effect's name came up, so
  # duplicated rows, and its parametric interval stopped inside the
  # noncentral F computation.
  split <- unique(effects_tbl$effect[duplicated(effects_tbl$effect)])
  if (length(split) > 0L) {
    strata <- unique(effects_tbl$stratum[effects_tbl$effect %in% split])
    .dmar_refuse_argument(
      sprintf(paste0(
   "Generalized eta squared needs each effect of an aovlist fit to be ",
   "tested in a single error stratum, as in a balanced design. These ",
   "effects are tested in more than one, as happens when observations ",
   "are missing within subjects or a covariate varies within them: %s ",
   "(strata %s). See the Details of ?eta_squared_generalized."),
   paste(split, collapse = ", "), paste(strata, collapse = ", ")), arg = NULL)
  }

  effect_names <- effects_tbl$effect
  ss_lookup    <- setNames(effects_tbl$SS_effect, effect_names)

  if (!is.null(observed)) {
    if (!is.character(observed)) {
      .dmar_refuse_argument(
        "'observed' must be a character vector of factor names.",
        arg = "observed")
    }
    missing_obs <- setdiff(observed, effect_names)
    if (length(missing_obs) > 0L) {
      .dmar_refuse_argument(
        sprintf("'observed' refers to effect(s) not in the aovlist: %s",
                paste(missing_obs, collapse = ", ")), arg = "observed")
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
