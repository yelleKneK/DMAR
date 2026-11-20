# Internal helpers shared by the ANOVA effect sizes eta_squared(),
# omega_squared(), ci_eta_squared(), and ci_omega_squared(), which report
# a total and a partial estimand (estimand = c("all", "total", "partial")),
# and by the generalized eta squared functions, which read the strata of an
# aovlist fit through the same helpers.
#
# The four public functions share one path:
#   1. .anova_es_estimand() checks the estimand asked for.
#   2. .anova_es_spec() reads the inputs, a fitted model or summary
#      statistics, into one row per effect (the F statistic, its degrees of
#      freedom, the sums of squares when they are known, and the size with
#      which a partial value is computed) plus what the inputs say about the
#      total: the corrected total sum of squares, or the reason it is not
#      known.
#   3. .anova_es_eta() or .anova_es_omega() computes each estimand, or NA
#      with the reason the inputs do not determine it.
#   4. .anova_es_table() lays out one row per effect with one column per
#      reported estimand (and, for an interval, its two limits beside it),
#      and refuses an estimand asked for by name that the inputs do not
#      determine.
#
# Population values and sources (read for this implementation; see the
# help pages for the citations): total eta squared SS_effect / SS_total and
# partial eta squared SS_effect / (SS_effect + SS_error) (Olejnik & Algina,
# 2000, Table 10, p. 268; Maxwell, Delaney, & Kelley, 2027, Eqs. 7.40 and
# 7.41); total omega squared (SS_effect - df_effect MS_error) /
# (SS_total + MS_error) and partial omega squared
# (SS_effect - df_effect MS_error) / (SS_effect + (N - df_effect) MS_error)
# for fixed effects in designs with between-subjects factors only (Olejnik
# & Algina, 2000, p. 266 and Table 10; Maxwell, Delaney, & Kelley, 2027,
# Eqs. 7.44 to 7.46).

# The estimands, in the order their columns appear.
.ANOVA_ES_ESTIMANDS <- c("total", "partial")

# The estimand asked for: "all" (the default, every estimand the inputs
# determine), "total", or "partial".
.anova_es_estimand <- function(estimand) {
  if (identical(estimand, c("all", "total", "partial"))) return("all")
  if (is.character(estimand) && length(estimand) == 1L && !is.na(estimand) &&
      estimand %in% c("all", "total", "partial")) {
    return(estimand)
  }
  .dmar_refuse_argument(c(
    "'estimand' must be one of \"all\" (every estimand the inputs ",
    "determine), \"total\", or \"partial\"; ", .dmar_given(estimand), "."),
    arg = "estimand", value = estimand)
}

# A single finite number, or a vector of them, checked against a bound.
.anova_es_check_number <- function(x, arg, scalar = TRUE, lower = 0,
                                   strict = TRUE) {
  ok <- is.numeric(x) && length(x) >= 1L && (!scalar || length(x) == 1L) &&
    all(is.finite(x)) && all(if (strict) x > lower else x >= lower)
  if (!ok) {
    what <- if (scalar) "a single" else "finite"
    sign <- if (strict) {
      if (lower == 0) "positive" else paste("greater than", lower)
    } else {
      if (lower == 0) "non-negative" else paste("at least", lower)
    }
    .dmar_refuse_argument(
      if (scalar) {
        sprintf("'%s' must be a single %s number; %s.", arg, sign,
                .dmar_given(x))
      } else {
        sprintf("'%s' must be one or more finite %s numbers; %s.", arg,
                sign, .dmar_given(x))
      },
      arg = arg, value = x)
  }
  invisible(x)
}

# The labels of the effects of the summary-statistic interfaces: the names
# of the effect-level vector when it has them, "overall" for a single
# unnamed effect (as in DMAR 1.0.x), and effect_1, effect_2, ... otherwise.
.anova_es_labels <- function(x) {
  k <- length(x)
  nm <- names(x)
  if (!is.null(nm) && all(!is.na(nm) & nzchar(nm))) return(nm)
  if (k == 1L) "overall" else paste0("effect_", seq_len(k))
}

# The reasons an estimand is not determined, for the NA columns of a table
# under estimand = "all" (its estimand_notes attribute and the line beneath
# the printed table) and for the refusal when it is asked for by name.
.anova_es_reason <- function(code, metric) {
  label <- if (metric == "eta_squared") "eta squared" else "omega squared"
  switch(code,
    F_total = paste0(
      "F and its degrees of freedom do not determine the total sum of ",
      "squares, so they do not determine the total ", label, " of an ",
      "effect in a design with other effects. Give SS_effect and SS_total",
      if (metric == "omega_squared") " with df_effect and MS_error",
      ", or the fitted model. In a one-way design the effect and the ",
      "error account for the total, and the partial value is the total."),
    no_intercept = paste0(
      "The model has no intercept, so its sums of squares do not ",
      "partition the corrected total sum of squares, and the total ",
      label, " is not determined."),
    aovlist_omega = paste0(
      "The total omega squared of Hays' formula is defined for designs ",
      "whose factors are all between subjects (Olejnik & Algina, 2000, ",
      "p. 266); a fit with an Error() term has more than one error term. ",
      "The partial value of each effect is computed within its own ",
      "stratum."),
    no_SS_total = paste0(
      "SS_total was not given, so the total ", label, " is not ",
      "determined."),
    no_SS_error = paste0(
      "SS_error was not given, so the partial eta squared is not ",
      "determined."),
    no_N = paste0(
      "N was not given, so the partial omega squared is not determined."),
    total_interval = paste0(
      "No confidence interval is formed for the total ", label, " of an ",
      "effect in a design with other effects or with an Error() term: ",
      "the noncentral F interval from the effect's F is an interval for ",
      "the partial value, and the total depends as well on the other ",
      "effects of the design (see the Details of ?ci_", metric, ")."),
    F_total_interval = paste0(
      "F and its degrees of freedom do not determine the total ", label,
      " of an effect in a design with other effects, so no interval is ",
      "formed for it. In a one-way design, where the two coincide, use ",
      "estimand = \"partial\", or give the fitted model."))
}

# The refusal of an estimand asked for by name that the inputs do not
# determine. Class dmar_estimand_undetermined; the condition carries the
# estimand and the reason.
.anova_es_refuse_estimand <- function(fn, estimand, reason) {
  .dmar_refuse(c(
    fn, "(estimand = \"", estimand, "\"): ", reason,
    " With estimand = \"all\" the table reports every estimand the inputs ",
    "determine and NA for the others."),
    "estimand_undetermined", estimand = estimand, reason = reason)
}


# Reading the inputs ---------------------------------------------------------

# One row per effect, from a fitted model or from summary statistics:
#   rows: effect, stratum (NA for a single-stratum fit or summaries),
#         F_value, df_effect, df_error, SS_effect, SS_error, MS_error,
#         size (the size with which a partial omega squared is computed and
#         an interval's limits are converted), N (the size the table
#         reports);
#   SS_total: the corrected total sum of squares, or NA;
#   total_reason: why the total is not determined (NULL when it is);
#   one_way: TRUE when the effect and the error account for the total (a
#         single-stratum fit with an intercept and one effect), where the
#         total and the partial coincide;
#   aovlist, route ("fit", "ss", or "F").
.anova_es_spec <- function(metric, object = NULL, SS_effect = NULL,
                           SS_error = NULL, SS_total = NULL, MS_error = NULL,
                           F_value = NULL, df_effect = NULL, df_error = NULL,
                           N = NULL, interval = FALSE) {
  summaries <- list(SS_effect = SS_effect, SS_error = SS_error,
                    SS_total = SS_total, MS_error = MS_error,
                    F_value = F_value, df_effect = df_effect,
                    df_error = df_error, N = N)
  given <- names(summaries)[!vapply(summaries, is.null, logical(1L))]
  if (!is.null(object)) {
    if (length(given)) {
      .dmar_refuse_argument(c(
        "'object' and summary statistics (", paste0("'", given, "'",
                                                    collapse = ", "),
        ") were both given; give the fitted model or the summary ",
        "statistics, not both."), arg = c("object", given))
    }
    return(.anova_es_fit(object, metric))
  }
  ss_args <- intersect(given, c("SS_effect", "SS_error", "SS_total",
                                "MS_error"))
  if (!is.null(F_value) && length(ss_args)) {
    .dmar_refuse_argument(c(
      "'F_value' and sums of squares (", paste0("'", ss_args, "'",
                                                collapse = ", "),
      ") were both given. They are two ways of giving the same effect; ",
      "give either F_value with its degrees of freedom or the sums of ",
      "squares."), arg = c("F_value", ss_args))
  }
  if (!is.null(F_value)) {
    return(.anova_es_from_F(metric, F_value, df_effect, df_error, N,
                            need_N = interval || metric == "omega_squared"))
  }
  if (length(ss_args)) {
    return(.anova_es_from_ss(metric, SS_effect, SS_error, SS_total,
                             MS_error, df_effect, df_error, N))
  }
  .dmar_refuse_argument(
    switch(if (interval) "interval" else metric,
      eta_squared = c(
        "Give a fitted aov, lm, or aovlist model as 'object'; or ",
        "'SS_effect' with 'SS_total' (total) and/or 'SS_error' (partial); ",
        "or 'F_value', 'df_effect', and 'df_error' (partial)."),
      omega_squared = c(
        "Give a fitted aov, lm, or aovlist model as 'object'; or ",
        "'SS_effect', 'df_effect', and 'MS_error' (or 'SS_error' with ",
        "'df_error'), with 'SS_total' (total) and/or 'N' (partial); or ",
        "'F_value', 'df_effect', 'df_error', and 'N' (partial)."),
      c("Give a fitted aov, lm, or aovlist model as 'object', or all of ",
        "'F_value', 'df_effect', 'df_error', and 'N'.")),
    arg = c("object", "F_value"))
}

# The F interface: F_value (one or more effects), df_effect (one per
# effect, or one for all), df_error, and, for omega squared and the
# intervals, N. The partial value only: F does not determine the total.
.anova_es_from_F <- function(metric, F_value, df_effect, df_error, N,
                             need_N) {
  missing <- c(if (is.null(df_effect)) "df_effect",
               if (is.null(df_error)) "df_error",
               if (need_N && is.null(N)) "N")
  if (length(missing)) {
    .dmar_refuse_argument(c(
      "'F_value' was given without ",
      paste0("'", missing, "'", collapse = " and "), ". The F interface ",
      "takes 'F_value', 'df_effect', 'df_error'",
      if (need_N) ", and 'N'", "."), arg = missing)
  }
  .anova_es_check_number(F_value, "F_value", scalar = FALSE, strict = FALSE)
  k <- length(F_value)
  .anova_es_check_number(df_effect, "df_effect", scalar = FALSE)
  if (!length(df_effect) %in% c(1L, k)) {
    .dmar_refuse_argument(sprintf(paste0(
      "'df_effect' must have one value for all effects or one for each of ",
      "the %d F values; it has %d."), k, length(df_effect)),
      arg = "df_effect")
  }
  df_effect <- rep_len(df_effect, k)
  .anova_es_check_number(df_error, "df_error")
  if (need_N) {
    # A within-subjects stratum that tests one effect has exactly
    # df_effect + df_error degrees of freedom, its size, so N may equal it.
    if (!is.numeric(N) || length(N) != 1L || !is.finite(N) ||
        N < max(df_effect) + df_error) {
      .dmar_refuse_argument(
        "'N' must be at least df_effect + df_error.", arg = "N", value = N)
    }
  } else {
    N <- NA_real_
  }
  rows <- data.frame(
    effect = .anova_es_labels(F_value), stratum = NA_character_,
    F_value = unname(as.numeric(F_value)), df_effect = df_effect,
    df_error = df_error, SS_effect = NA_real_, SS_error = NA_real_,
    MS_error = NA_real_, size = N, N = N, stringsAsFactors = FALSE)
  list(rows = rows, SS_total = NA_real_,
       total_reason = .anova_es_reason("F_total", metric),
       one_way = FALSE, aovlist = FALSE, route = "F")
}

# The sums-of-squares interface. Eta squared: SS_effect with SS_total for
# the total and SS_error for the partial (df_effect and df_error only fill
# the F_value and degrees of freedom columns). Omega squared: SS_effect,
# df_effect, and the error mean square (MS_error, or SS_error with
# df_error), with SS_total for the total and N for the partial.
.anova_es_from_ss <- function(metric, SS_effect, SS_error, SS_total,
                              MS_error, df_effect, df_error, N) {
  omega <- metric == "omega_squared"
  if (is.null(SS_effect)) {
    .dmar_refuse_argument(c(
      "The sums-of-squares interface needs 'SS_effect', the sum of ",
      "squares of the effect."), arg = "SS_effect")
  }
  .anova_es_check_number(SS_effect, "SS_effect", scalar = FALSE,
                         strict = FALSE)
  k <- length(SS_effect)
  if (!is.null(SS_error)) .anova_es_check_number(SS_error, "SS_error")
  if (!is.null(SS_total)) .anova_es_check_number(SS_total, "SS_total")
  if (!is.null(MS_error)) .anova_es_check_number(MS_error, "MS_error")
  if (!is.null(df_error)) .anova_es_check_number(df_error, "df_error")
  if (!is.null(df_effect)) {
    .anova_es_check_number(df_effect, "df_effect", scalar = FALSE)
    if (!length(df_effect) %in% c(1L, k)) {
      .dmar_refuse_argument(sprintf(paste0(
        "'df_effect' must have one value for all effects or one for each ",
        "of the %d values of SS_effect; it has %d."), k, length(df_effect)),
        arg = "df_effect")
    }
    df_effect <- rep_len(df_effect, k)
  }
  if (omega) {
    if (is.null(df_effect)) {
      .dmar_refuse_argument(c(
        "Omega squared from sums of squares needs 'df_effect', the ",
        "degrees of freedom of the effect."), arg = "df_effect")
    }
    if (!is.null(MS_error) && !is.null(SS_error)) {
      .dmar_refuse_argument(c(
        "'MS_error' and 'SS_error' were both given. They are two ways of ",
        "giving the error mean square; give 'MS_error', or 'SS_error' ",
        "with 'df_error'."), arg = c("MS_error", "SS_error"))
    }
    if (is.null(MS_error)) {
      if (is.null(SS_error) || is.null(df_error)) {
        .dmar_refuse_argument(c(
          "Omega squared needs the error mean square: give 'MS_error', or ",
          "'SS_error' with 'df_error'."),
          arg = c("MS_error", "SS_error", "df_error"))
      }
      MS_error <- SS_error / df_error
    } else if (!is.null(df_error)) {
      SS_error <- MS_error * df_error
    }
    if (is.null(SS_total) && is.null(N)) {
      .dmar_refuse_argument(c(
        "Give 'SS_total' (for the total omega squared), 'N' (for the ",
        "partial), or both."), arg = c("SS_total", "N"))
    }
    if (!is.null(N)) {
      floor_N <- max(df_effect) + if (is.null(df_error)) 0 else df_error
      ok <- is.numeric(N) && length(N) == 1L && is.finite(N) &&
        (if (is.null(df_error)) N > floor_N else N >= floor_N)
      if (!ok) {
        .dmar_refuse_argument(
          if (is.null(df_error)) "'N' must be greater than df_effect." else
            "'N' must be at least df_effect + df_error.",
          arg = "N", value = N)
      }
    }
  } else {
    if (!is.null(MS_error)) {
      .dmar_refuse_argument(
        "'MS_error' is not an input of eta squared.", arg = "MS_error")
    }
    if (is.null(SS_error) && is.null(SS_total)) {
      .dmar_refuse_argument(c(
        "Give 'SS_total' (for the total eta squared), 'SS_error' (for the ",
        "partial), or both, with 'SS_effect'."),
        arg = c("SS_total", "SS_error"))
    }
    if (!is.null(N)) {
      .dmar_refuse_argument("'N' is not an input of eta squared.",
                            arg = "N")
    }
  }
  # The corrected total holds the effect and the error (for every type of
  # sums of squares, since a reduced model that keeps the intercept leaves
  # no more than the corrected total), so a total below their sum
  # contradicts them.
  if (!is.null(SS_total)) {
    part <- SS_effect + if (is.null(SS_error)) 0 else SS_error
    if (any(part > SS_total * (1 + 1e-8))) {
      .dmar_refuse_argument(c(
        "'SS_total' (", format(SS_total), ") is smaller than ",
        if (is.null(SS_error)) "SS_effect" else "SS_effect + SS_error",
        " (", format(max(part)), "). The corrected total sum of squares ",
        "contains the effect and the error; check the values (and that ",
        "SS_total is the corrected total)."),
        arg = c("SS_total", "SS_effect"))
    }
  }
  F_value <- if (!is.null(df_effect) && !is.null(SS_error) &&
                 !is.null(df_error)) {
    (SS_effect / df_effect) / (SS_error / df_error)
  } else if (omega) {
    (SS_effect / df_effect) / MS_error
  } else {
    NA_real_
  }
  rows <- data.frame(
    effect = .anova_es_labels(SS_effect), stratum = NA_character_,
    F_value = unname(F_value),
    df_effect = if (is.null(df_effect)) NA_real_ else df_effect,
    df_error = if (is.null(df_error)) NA_real_ else df_error,
    SS_effect = unname(as.numeric(SS_effect)),
    SS_error = if (is.null(SS_error)) NA_real_ else SS_error,
    MS_error = if (is.null(MS_error)) NA_real_ else MS_error,
    size = if (is.null(N)) NA_real_ else N,
    N = if (is.null(N)) NA_real_ else N, stringsAsFactors = FALSE)
  list(rows = rows, SS_total = if (is.null(SS_total)) NA_real_ else SS_total,
       total_reason = if (is.null(SS_total))
         .anova_es_reason("no_SS_total", metric),
       one_way = FALSE, aovlist = FALSE, route = "ss")
}

# A fitted aov or lm (one error stratum) or aovlist (aov() with an Error()
# term). The F statistics and sums of squares are those of anova() and
# summary(), which are sequential (Type I).
.anova_es_fit <- function(object, metric) {
  if (inherits(object, "aovlist")) {
    eff <- .aovlist_effects_table(object)
    eff <- eff[!is.na(eff$F_value) & eff$F_value >= 0, , drop = FALSE]
    if (nrow(eff) == 0L) {
      .dmar_refuse_argument(
        "No effects with a computable F statistic were found in the aovlist.",
        arg = NULL)
    }
    sizes <- .aovlist_stratum_sizes(object)
    N <- .aovlist_nobs(object)
    rows <- data.frame(
      effect = eff$effect, stratum = eff$stratum, F_value = eff$F_value,
      df_effect = eff$df_effect, df_error = eff$df_error,
      SS_effect = eff$SS_effect, SS_error = eff$SS_error,
      MS_error = eff$SS_error / eff$df_error,
      size = unname(sizes[eff$stratum]), N = N, stringsAsFactors = FALSE)
    return(list(rows = rows, SS_total = .aovlist_total_ss(object),
                total_reason = if (metric == "omega_squared")
                  .anova_es_reason("aovlist_omega", metric),
                one_way = FALSE, aovlist = TRUE, route = "fit"))
  }
  if (!inherits(object, c("aov", "lm"))) {
    .dmar_refuse_argument(
      "'object' must be an aov, lm, or aovlist fit.", arg = "object")
  }
  tbl <- stats::anova(object)
  rn <- trimws(rownames(tbl))
  resid_row <- which(rn == "Residuals")
  if (length(resid_row) == 0L) {
    .dmar_refuse_argument(c(
      "Could not find a 'Residuals' row in anova(object); is this a ",
      "standard fixed-effects ANOVA?"), arg = NULL)
  }
  df_error <- tbl[resid_row, "Df"]
  SS_error <- tbl[resid_row, "Sum Sq"]
  effect_rows <- which(rn != "Residuals")
  keep <- effect_rows[!is.na(tbl[effect_rows, "F value"]) &
                        tbl[effect_rows, "F value"] >= 0]
  if (length(keep) == 0L) {
    .dmar_refuse_argument(
      "No effects with a computable F statistic were found in the model.",
      arg = NULL)
  }
  N <- stats::nobs(object)
  intercept <- isTRUE(attr(stats::terms(object), "intercept") == 1L)
  rows <- data.frame(
    effect = rn[keep], stratum = NA_character_,
    F_value = tbl[keep, "F value"], df_effect = tbl[keep, "Df"],
    df_error = df_error, SS_effect = tbl[keep, "Sum Sq"],
    SS_error = SS_error, MS_error = SS_error / df_error,
    size = N, N = N, stringsAsFactors = FALSE)
  list(rows = rows,
       SS_total = if (intercept) sum(tbl[, "Sum Sq"]) else NA_real_,
       total_reason = if (!intercept) .anova_es_reason("no_intercept", metric),
       one_way = intercept && length(effect_rows) == 1L,
       aovlist = FALSE, route = "fit")
}


# The estimands ---------------------------------------------------------------

# Each estimand of eta squared: a list with `value` (one per effect, NA
# when not determined) and `reason` (NULL when determined).
.anova_es_eta <- function(spec) {
  r <- spec$rows
  partial <- if (spec$route == "ss") {
    if (all(is.na(r$SS_error))) {
      list(value = rep(NA_real_, nrow(r)),
           reason = .anova_es_reason("no_SS_error", "eta_squared"))
    } else {
      list(value = r$SS_effect / (r$SS_effect + r$SS_error), reason = NULL)
    }
  } else {
    # The F form (Maxwell, Delaney, & Kelley, 2027, Eq. 7.42), which is
    # SS_effect / (SS_effect + SS_error), as DMAR 1.0.x computed it.
    num <- r$df_effect * r$F_value
    list(value = num / (num + r$df_error), reason = NULL)
  }
  total <- if (is.null(spec$total_reason)) {
    list(value = r$SS_effect / spec$SS_total, reason = NULL)
  } else {
    list(value = rep(NA_real_, nrow(r)), reason = spec$total_reason)
  }
  list(total = total, partial = partial)
}

# Each estimand of omega squared, truncated at 0 (a negative estimate,
# which occurs when F < 1, is reported as 0, as in DMAR 1.0.x).
.anova_es_omega <- function(spec) {
  r <- spec$rows
  partial <- if (spec$route == "ss") {
    if (all(is.na(r$size))) {
      list(value = rep(NA_real_, nrow(r)),
           reason = .anova_es_reason("no_N", "omega_squared"))
    } else {
      num <- r$SS_effect - r$df_effect * r$MS_error
      list(value = pmax(0, num / (r$SS_effect + (r$size - r$df_effect) *
                                    r$MS_error)), reason = NULL)
    }
  } else {
    # df_effect (F - 1) / (df_effect (F - 1) + size), the same value from
    # the F statistic (Maxwell, Delaney, & Kelley, 2027, Eq. 7.46), with
    # the size of the effect's stratum for an aovlist fit.
    num <- r$df_effect * (r$F_value - 1)
    list(value = pmax(0, num / (num + r$size)), reason = NULL)
  }
  total <- if (is.null(spec$total_reason)) {
    num <- r$SS_effect - r$df_effect * r$MS_error
    list(value = pmax(0, num / (spec$SS_total + r$MS_error)), reason = NULL)
  } else {
    list(value = rep(NA_real_, nrow(r)), reason = spec$total_reason)
  }
  list(total = total, partial = partial)
}

# The noncentral F limits of the partial value of each effect, converted
# with the size of its stratum, lambda / (lambda + size) (Steiger, 2004,
# Eqs. 16 and 37). One clamp warning per call, not one per effect.
.anova_es_partial_limits <- function(spec, alpha_lower, alpha_upper,
                                     caller, quantity) {
  r <- spec$rows
  lims <- .warn_ncf_clamp_once(lapply(seq_len(nrow(r)), function(i) {
    ncp <- .ci_nc_F_for(caller = caller, quantity = quantity,
                        F_value = r$F_value[i], df_1 = r$df_effect[i],
                        df_2 = r$df_error[i], alpha_lower = alpha_lower,
                        alpha_upper = alpha_upper, conf_level = NULL)
    lo <- ncp$value[ncp$term == "ci_lower"]
    up <- ncp$value[ncp$term == "ci_upper"]
    c(lo / (lo + r$size[i]),
      if (is.infinite(up)) 1 else up / (up + r$size[i]))
  }))
  list(lower = vapply(lims, `[`, numeric(1L), 1L),
       upper = vapply(lims, `[`, numeric(1L), 2L))
}

# The interval of each estimand: the partial's from the noncentral F, and
# the total's only where the two coincide (spec$one_way), else NA with the
# reason.
.anova_es_limits <- function(spec, estimand, metric, alpha_lower,
                             alpha_upper) {
  k <- nrow(spec$rows)
  label <- if (metric == "eta_squared") "eta squared" else "omega squared"
  need_partial <- estimand %in% c("all", "partial") ||
    (estimand == "total" && spec$one_way)
  partial <- if (need_partial) {
    .anova_es_partial_limits(
      spec, alpha_lower, alpha_upper, caller = paste0("ci_", metric),
      quantity = if (estimand == "total") paste("total", label) else
        paste("partial", label))
  }
  total <- if (spec$one_way) {
    c(partial, list(reason = NULL))
  } else {
    list(lower = rep(NA_real_, k), upper = rep(NA_real_, k),
         reason = .anova_es_reason(
           if (spec$route == "F") "F_total_interval" else "total_interval",
           metric))
  }
  list(total = total, partial = c(partial, list(reason = NULL)))
}


# The table -------------------------------------------------------------------

# One row per effect: effect, then for each reported estimand its column
# (<metric>_total, <metric>_partial) and, for an interval, its limits
# (<metric>_total_ci_lower, ...), then the stratum of an aovlist fit, the F
# statistic and its degrees of freedom, and N (omega squared and the
# intervals). An estimand asked for by name that the inputs do not
# determine is refused; under "all" its columns are NA and the reason is
# recorded in the estimand_notes attribute, which the line beneath the
# printed table states.
.anova_es_table <- function(fn, metric, spec, est, estimand, limits = NULL,
                            conf_level = NULL, tail_areas = NULL) {
  reported <- if (estimand == "all") .ANOVA_ES_ESTIMANDS else estimand
  if (estimand != "all") {
    reason <- est[[estimand]]$reason
    if (is.null(reason) && !is.null(limits)) reason <- limits[[estimand]]$reason
    if (!is.null(reason)) .anova_es_refuse_estimand(fn, estimand, reason)
  }
  r <- spec$rows
  out <- data.frame(effect = r$effect, stringsAsFactors = FALSE)
  notes <- character(0)
  for (e in reported) {
    col <- paste0(metric, "_", e)
    out[[col]] <- est[[e]]$value
    if (!is.null(est[[e]]$reason)) notes[col] <- est[[e]]$reason
    if (!is.null(limits)) {
      out[[paste0(col, "_ci_lower")]] <- limits[[e]]$lower
      out[[paste0(col, "_ci_upper")]] <- limits[[e]]$upper
      if (is.null(est[[e]]$reason) && !is.null(limits[[e]]$reason))
        notes[paste0(col, "_ci")] <- limits[[e]]$reason
    }
  }
  if (spec$aovlist) out$stratum <- r$stratum
  out$F_value <- r$F_value
  out$df_effect <- r$df_effect
  out$df_error <- r$df_error
  if (metric == "omega_squared" || !is.null(limits)) out$N <- r$N
  rownames(out) <- NULL
  attr(out, "estimand") <- estimand
  if (length(notes)) attr(out, "estimand_notes") <- notes
  .as_dmar_tbl(out, conf_level = conf_level, tail_areas = tail_areas,
               subclass = if (is.null(limits)) "dmar_anova_es" else
                 "dmar_ci_anova")
}

# The line beneath a printed table for each estimand that is NA, naming the
# column and the reason (the estimand_notes attribute).
.anova_es_footer_lines <- function(x) {
  notes <- attr(x, "estimand_notes", exact = TRUE)
  if (!is.character(notes) || !length(notes)) return(character(0))
  keep <- vapply(names(notes), function(nm) {
    col <- sub("_ci$", "", nm)
    any(c(col, paste0(col, "_ci_lower")) %in% names(x))
  }, logical(1L))
  notes <- notes[keep]
  if (!length(notes)) return(character(0))
  what <- ifelse(grepl("_ci$", names(notes)),
                 paste0(sub("_ci$", "", names(notes)), " limits"),
                 names(notes))
  paste0("NA ", what, ": ", unname(notes))
}

#' @exportS3Method .dmar_tbl_footer dmar_anova_es
#' @noRd
.dmar_tbl_footer.dmar_anova_es <- function(x, ...) {
  .dmar_footer_add(NextMethod(), .anova_es_footer_lines(x), "estimand_note")
}

#' @exportS3Method .dmar_tbl_footer dmar_ci_anova
#' @noRd
.dmar_tbl_footer.dmar_ci_anova <- function(x, ...) {
  .dmar_footer_add(NextMethod(), .anova_es_footer_lines(x), "estimand_note")
}


# Alpha / confidence-level plumbing ------------------------------------------

# Shared by ci_eta_squared(), ci_omega_squared(), and
# ci_eta_squared_generalized(). Returns the two tail areas and the
# confidence level of the interval they define, conf_level itself or, when
# the tail areas are supplied, 1 - alpha_lower - alpha_upper, which is the
# level the table records.
.ci_eta_squared_alphas <- function(conf_level, alpha_lower, alpha_upper) {
  if (is.null(alpha_lower) && is.null(alpha_upper)) {
    alpha_lower <- (1 - conf_level) / 2
    alpha_upper <- (1 - conf_level) / 2
  } else if (!is.null(alpha_lower) && !is.null(alpha_upper)) {
    conf_level <- 1 - alpha_lower - alpha_upper
  } else if (!is.null(alpha_lower) && is.null(alpha_upper)) {
    .dmar_refuse_argument(
      "'alpha_lower' specified but not 'alpha_upper'.",
      arg = c("alpha_lower", "alpha_upper"))
  } else {
    .dmar_refuse_argument(
      "'alpha_upper' specified but not 'alpha_lower'.",
      arg = c("alpha_upper", "alpha_lower"))
  }
  if (alpha_lower < 0 || alpha_lower > 0.5) .dmar_refuse_argument(
    "'alpha_lower' must be in [0, 0.5].", arg = "alpha_lower")
  if (alpha_upper < 0 || alpha_upper > 0.5) .dmar_refuse_argument(
    "'alpha_upper' must be in [0, 0.5].", arg = "alpha_upper")
  if (conf_level <= 0 || conf_level >= 1)   .dmar_refuse_argument(
    "'conf_level' must be strictly between 0 and 1.", arg = "conf_level")
  list(alpha_lower = alpha_lower, alpha_upper = alpha_upper,
       conf_level = conf_level)
}


# Strata of an aovlist fit -----------------------------------------------------

# Walk an aovlist and return a data.frame with one row per effect across all
# strata. Uses summary(object), which is the canonical access path for
# stratum-wise anova tables in an aovlist, rather than anova(object[[i]]),
# since the latter fails on the (Intercept) and pure-residual strata. Skips
# strata without a Residuals row and skips any effect named "(Intercept)".
.aovlist_effects_table <- function(object) {
  s <- summary(object)
  rows <- list()
  stratum_names <- names(s)
  for (i in seq_along(s)) {
    stratum_name <- stratum_names[i]
    if (is.null(stratum_name) || !nzchar(stratum_name)) {
      stratum_name <- as.character(i)
    }
    # Strip the "Error: " prefix that summary.aovlist prepends for readability.
    stratum_name <- sub("^Error:\\s*", "", stratum_name)

    inner <- s[[i]]
    if (is.list(inner) && !is.data.frame(inner)) inner <- inner[[1L]]
    if (is.null(inner) || !"Sum Sq" %in% colnames(inner)) next
    # summary.aov rownames are padded with trailing spaces ("time      ",
    # "Residuals "); strip before any comparison.
    rownames(inner) <- trimws(rownames(inner))
    resid_row <- which(rownames(inner) == "Residuals")
    if (length(resid_row) == 0L) next
    df_err <- inner[resid_row, "Df"]
    ss_err <- inner[resid_row, "Sum Sq"]
    effect_rows <- setdiff(rownames(inner), c("Residuals", "(Intercept)"))
    effect_rows <- effect_rows[nzchar(effect_rows)]
    for (eff in effect_rows) {
      rows[[length(rows) + 1L]] <- data.frame(
        effect    = eff,
        stratum   = stratum_name,
        F_value   = inner[eff, "F value"],
        df_effect = inner[eff, "Df"],
        SS_effect = inner[eff, "Sum Sq"],
        df_error  = df_err,
        SS_error  = ss_err,
        stringsAsFactors = FALSE
      )
    }
  }
  if (length(rows) == 0L) {
    .dmar_refuse_argument(
      "No effects with a computable F statistic found in the aovlist.",
      arg = NULL)
  }
  do.call(rbind, rows)
}

# Recover the total number of observations from an aovlist. stats::nobs has
# no method for aovlist, so we sum the degrees of freedom across all strata
# returned by summary() and add 1 for the intercept (which is absorbed into
# the (Intercept) stratum that summary() suppresses). For balanced designs
# this gives the textbook total N; for unbalanced designs it is exact too
# because the strata Df partition is complete by construction.
.aovlist_nobs <- function(object) {
  s <- summary(object)
  total_df <- 0
  for (i in seq_along(s)) {
    inner <- s[[i]]
    if (is.list(inner) && !is.data.frame(inner)) inner <- inner[[1L]]
    if (is.null(inner) || !"Df" %in% colnames(inner)) next
    total_df <- total_df + sum(inner[, "Df"], na.rm = TRUE)
  }
  total_df + 1L
}

# The corrected total sum of squares of an aovlist: the sum of every sum of
# squares in every stratum summary() reports, effects and residuals, since
# the strata are orthogonal projections that together span the deviations
# from the grand mean (the (Intercept) stratum, which summary() leaves out,
# holds the grand mean). It equals sum((y - mean(y))^2) for balanced and
# unbalanced designs alike (checked in tests/testthat/test-eta_squared.R).
.aovlist_total_ss <- function(object) {
  s <- summary(object)
  total <- 0
  for (i in seq_along(s)) {
    inner <- s[[i]]
    if (is.list(inner) && !is.data.frame(inner)) inner <- inner[[1L]]
    if (is.null(inner) || !"Sum Sq" %in% colnames(inner)) next
    total <- total + sum(inner[, "Sum Sq"], na.rm = TRUE)
  }
  total
}

# Sum the residual sums of squares across all error strata of an aovlist.
# Used by the generalized eta squared denominator for within-subjects designs.
.aovlist_total_residual_ss <- function(object) {
  s <- summary(object)
  total <- 0
  for (i in seq_along(s)) {
    inner <- s[[i]]
    if (is.list(inner) && !is.data.frame(inner)) inner <- inner[[1L]]
    if (is.null(inner) || !"Sum Sq" %in% colnames(inner)) next
    rownames(inner) <- trimws(rownames(inner))
    resid_row <- which(rownames(inner) == "Residuals")
    if (length(resid_row) == 0L) next
    total <- total + inner[resid_row, "Sum Sq"]
  }
  total
}

# The size with which the noncentrality limits of the effects in each
# stratum of an aovlist are converted, lambda / (lambda + size), named by
# stratum as .aovlist_effects_table() names them. A stratum's total
# degrees of freedom, those of its effects and of its residual, count the
# independent contrasts it holds. A stratum whose Error() term contains no
# other stratum's term, such as id in Error(id / occasion), lies between
# the units of that term: it analyzes their means, centered only at the
# grand mean, and is converted with the number of units, one more than its
# degrees of freedom, as a single-stratum fit is converted with N. Every
# other stratum, such as id:occasion or Within, holds contrasts within
# those units and is converted with its degrees of freedom. Under
# sphericity each makes lambda / (lambda + size) the value the stratum's
# sample partial eta squared converges to as units are added (see the
# Details of ?ci_eta_squared). The partial omega squared of an aovlist
# effect is computed with the same size.
.aovlist_stratum_sizes <- function(object) {
  s <- summary(object)
  nm <- names(s)
  strata <- vapply(seq_along(s), function(i) {
    name_i <- if (is.null(nm) || !nzchar(nm[i])) as.character(i) else nm[i]
    sub("^Error:\\s*", "", name_i)
  }, character(1))
  df_total <- vapply(s, function(inner) {
    if (is.list(inner) && !is.data.frame(inner)) inner <- inner[[1L]]
    if (is.null(inner) || !"Df" %in% colnames(inner)) return(NA_real_)
    sum(inner[, "Df"], na.rm = TRUE)
  }, numeric(1), USE.NAMES = FALSE)
  terms <- strsplit(strata, ":", fixed = TRUE)
  between_units <- vapply(seq_along(strata), function(i) {
    if (strata[i] == "Within") return(FALSE)
    nested <- vapply(seq_along(strata)[-i], function(j)
      strata[j] != "Within" && length(terms[[j]]) < length(terms[[i]]) &&
        all(terms[[j]] %in% terms[[i]]), logical(1))
    !any(nested)
  }, logical(1))
  stats::setNames(df_total + between_units, strata)
}
