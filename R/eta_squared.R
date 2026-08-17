# Eta squared (effect size for ANOVA).
#' Eta Squared (Effect Size for ANOVA)
#'
#' Computes the sample eta squared (\eqn{\eta^2}), the proportion of variance in
#' the dependent variable accounted for by a fixed effect. Accepts either the
#' raw ANOVA summary (\emph{F} and the effect and error degrees of freedom) or
#' a fitted model object. Supports both between-subjects designs (single-stratum
#' \code{\link[stats]{aov}} or \code{\link[stats]{lm}} fits) and
#' within-subjects / mixed designs (\code{aovlist} fits produced by
#' \code{\link[stats]{aov}} with an \code{Error()} term in the formula). For
#' factorial and within-subjects designs the function returns \emph{partial}
#' \eqn{\eta^2} per effect (one row per non-\code{Residuals} effect across all
#' strata), each computed against its own stratum's error term.
#'
#' The confidence interval is provided by the separate
#' \code{\link{ci_eta_squared}}, paralleling the existing
#' \code{\link{smd}}/\code{\link{ci_smd}} pairing.
#'
#' @param object Optional. A fitted model object of class
#'   \code{\link[stats]{aov}}, \code{\link[stats]{lm}}, or
#'   \code{aovlist} (a multi-stratum aov fit such as
#'   \code{aov(y ~ A + Error(subject/A), data = d)}). When supplied, the
#'   function loops over the non-\code{Residuals} effects across all error
#'   strata and returns one row per effect, with the stratum reported.
#' @param F_value Observed \emph{F}-value from the fixed-effects ANOVA
#'   (ignored if \code{object} is supplied).
#' @param df_effect Numerator degrees of freedom for the effect
#'   (ignored if \code{object} is supplied).
#' @param df_error Error (residual) degrees of freedom
#'   (ignored if \code{object} is supplied).
#'
#' @return A \code{data.frame} with one row per effect. For
#'   single-stratum (\code{aov}/\code{lm}) fits and the raw-argument
#'   interface the columns are \code{effect}, \code{eta_squared},
#'   \code{F_value}, \code{df_effect}, \code{df_error}. For multi-stratum
#'   (\code{aovlist}) fits an additional \code{stratum} column reports
#'   which error stratum each effect's \emph{F} test came from. When the
#'   raw-argument interface is used, \code{effect} is \code{"overall"}.
#'
#' @details
#' \strong{Point estimate.} The function uses the algebraically equivalent
#' \emph{F}-and-df form
#' \deqn{\hat{\eta}^2 = \frac{df_{\text{effect}} \cdot F}{df_{\text{effect}} \cdot F + df_{\text{error}}},}
#' which equals \eqn{\mathit{SS}_{\text{effect}} / (\mathit{SS}_{\text{effect}} +
#' \mathit{SS}_{\text{error}})}. In a one-way ANOVA this is also
#' \eqn{\mathit{SS}_{\text{effect}} / \mathit{SS}_{\text{total}}}, the
#' conventional \emph{total} \eqn{\eta^2}. In a factorial design the same
#' expression yields \emph{partial} \eqn{\eta^2} for each effect, because
#' \eqn{\mathit{SS}_{\text{error}}} appears in the denominator instead of
#' \eqn{\mathit{SS}_{\text{total}}}; this matches the convention used by
#' \code{\link{ci_omega_squared}}.
#'
#' \strong{Designs supported.}
#' \itemize{
#'   \item \emph{Between-subjects ANOVA} (one-way or factorial): supply
#'     either a fitted \code{aov}/\code{lm} model or the raw \emph{F} and
#'     degrees of freedom for a single effect.
#'   \item \emph{Within-subjects or mixed ANOVA}: supply a fitted
#'     \code{aovlist} model produced with an \code{Error()} term, e.g.\
#'     \code{aov(y ~ time + Error(subject/time), data = d)}. The function
#'     walks every error stratum returned by \code{summary(object)} and
#'     reports each effect with the stratum's residual \emph{df}, so the
#'     \eqn{\eta^2} value uses the stratum's specific error term. The
#'     reported \code{stratum} column tells you which one.
#' }
#' For more advanced model classes (\code{lmerMod}, \code{lme}, etc.) the
#' fitted-model interface is not yet supported; supply the relevant
#' \emph{F} and degrees of freedom via the raw interface.
#'
#' \strong{Sums of squares in factorial designs.} \code{anova()} on an
#' \code{aov}/\code{lm} uses Type I (sequential) sums of squares. For
#' balanced designs all three types agree; for unbalanced designs they
#' differ. If Type II or III \emph{F}-values are required, compute them
#' with e.g.\ \code{car::Anova(object, type = 3)} and pass the relevant
#' \emph{F} and degrees of freedom into the raw-argument interface.
#'
#' \strong{Generalized eta squared, comparable across designs.} The basic
#' \eqn{\eta^2} (and partial \eqn{\eta^2}) returned here are not comparable
#' across studies that differ in factor structure. See
#' \code{\link{eta_squared_generalized}} for a comparable alternative
#' (Olejnik & Algina, 2003; Bakeman, 2005).
#'
#' @references
#' Bakeman, R. (2005). Recommended effect size statistics for repeated
#'   measures designs. \emph{Behavior Research Methods, 37}(3), 379--384.
#'   \doi{10.3758/BF03192707}
#'
#' Cohen, J. (1988). \emph{Statistical power analysis for the behavioral
#'   sciences} (2nd ed.). Hillsdale, NJ: Lawrence Erlbaum.
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
#' Steiger, J. H. (2004). Beyond the \emph{F} test: Effect size confidence
#'   intervals and tests of close fit in the analysis of variance and contrast
#'   analysis. \emph{Psychological Methods, 9}(2), 164--182.
#'   \doi{10.1037/1082-989X.9.2.164}
#'
#' @seealso \code{\link{ci_eta_squared}}, \code{\link{eta_squared_partial}},
#'   \code{\link{ci_omega_squared}}, \code{\link{ci_pvaf}}
#'
#' @examples
#' # 1. Raw-argument interface. Bargman's (1970) 5-group one-way ANOVA,
#' #        also used in Venables (1975), Fleishman (1980), and Steiger (2004):
#' #        11 subjects per group, observed F = 11.221.
#' eta_squared(F_value = 11.221, df_effect = 4, df_error = 50)
#'
#' # 2. One way ANOVA from a fitted model (depression_bdi: three
#' #        treatment arms, 10 per arm, N = 30).
#' fit_one <- aov(bdi_post ~ condition, data = depression_bdi)
#' eta_squared(fit_one)
#'
#' # 3. Factorial ANOVA: partial eta squared per effect (pygmalion
#' #        data: expectancy treatment x grade, 2 x 6 with unequal cell
#' #        sizes, N = 310). The treatment is manipulated; grade is a
#' #        measured classification of the pupils.
#' fit_factorial <- aov(iq_8 ~ treatment * factor(grade), data = pygmalion)
#' eta_squared(fit_factorial)
#'
#' # 4. Within-subjects (repeated measures) ANOVA. Simulated 20-subject
#' #        x 3-time design; each subject is measured at Pre, Mid, Post.
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
#' eta_squared(fit_rm)   # 'stratum' column identifies the within-subjects error
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @keywords htest design
#'
#' @family effect size estimates
#'
#' @export
#' @import stats

eta_squared <- function(
  object    = NULL,
  F_value   = NULL,
  df_effect = NULL,
  df_error  = NULL
) {
  if (!is.null(object)) {
    return(.eta_squared_from_model(object, value_name = "eta_squared"))
  }
  if (is.null(F_value) || is.null(df_effect) || is.null(df_error)) {
    stop("Either provide 'object' (a fitted aov or lm) or all of 'F_value', 'df_effect', and 'df_error'.")
  }
  .eta_squared_one(F_value, df_effect, df_error,
                      effect_label = "overall", value_name = "eta_squared")
}


# Single-effect worker. Returns a 1-row data.frame in which the second column
# is named via `value_name` so that eta_squared() and eta_squared_partial()
# can share the same machinery while producing self-documenting output.
.eta_squared_one <- function(F_value, df_effect, df_error,
                                effect_label, value_name) {
  if (!is.numeric(F_value) || length(F_value) != 1L || F_value < 0)
    stop("'F_value' must be a single non-negative number.")
  if (!is.numeric(df_effect) || length(df_effect) != 1L || df_effect <= 0)
    stop("'df_effect' must be a single positive number.")
  if (!is.numeric(df_error) || length(df_error) != 1L || df_error <= 0)
    stop("'df_error' must be a single positive number.")

  num <- df_effect * F_value
  val <- num / (num + df_error)

  out <- data.frame(
    effect     = effect_label,
    placeholder = val,
    F_value    = F_value,
    df_effect  = df_effect,
    df_error   = df_error,
    stringsAsFactors = FALSE,
    row.names  = NULL
  )
  names(out)[2] <- value_name
  .as_dmar_tbl(out)
}


# Model-based worker: dispatches on aovlist (multi-stratum) vs aov/lm.
.eta_squared_from_model <- function(object, value_name) {
  if (inherits(object, "aovlist")) {
    return(.eta_squared_from_aovlist(object, value_name))
  }
  if (!inherits(object, c("aov", "lm"))) {
    stop("'object' must be an aov, lm, or aovlist fit.")
  }
  tbl <- stats::anova(object)

  resid_row <- which(rownames(tbl) == "Residuals")
  if (length(resid_row) == 0L) {
    stop("Could not find a 'Residuals' row in anova(object); is this a standard fixed-effects ANOVA?")
  }
  df_error <- tbl[resid_row, "Df"]
  effect_rows <- setdiff(rownames(tbl), "Residuals")

  results <- lapply(effect_rows, function(effect) {
    F_val <- tbl[effect, "F value"]
    df_e  <- tbl[effect, "Df"]
    if (is.na(F_val) || F_val < 0) return(NULL)
    .eta_squared_one(F_val, df_e, df_error,
                        effect_label = effect, value_name = value_name)
  })
  results <- Filter(Negate(is.null), results)
  if (length(results) == 0L) {
    stop("No effects with a computable F statistic were found in the model.")
  }

  out <- do.call(rbind, results)
  rownames(out) <- NULL
  .as_dmar_tbl(out)
}


# Within-subjects (multi-stratum) worker: each effect uses its own stratum's
# error term, so partial eta squared is computed per-stratum.
.eta_squared_from_aovlist <- function(object, value_name) {
  effects <- .aovlist_effects_table(object)
  results <- lapply(seq_len(nrow(effects)), function(i) {
    row <- effects[i, , drop = FALSE]
    if (is.na(row$F_value) || row$F_value < 0) return(NULL)
    out <- .eta_squared_one(row$F_value, row$df_effect, row$df_error,
                            effect_label = row$effect, value_name = value_name)
    out$stratum <- row$stratum
    out
  })
  results <- Filter(Negate(is.null), results)
  if (length(results) == 0L) {
    stop("No computable effects found in the aovlist.")
  }
  out <- do.call(rbind, results)
  # Reorder columns: effect, value_name, stratum, F_value, df_effect, df_error
  out <- out[, c("effect", value_name, "stratum", "F_value", "df_effect", "df_error")]
  rownames(out) <- NULL
  .as_dmar_tbl(out)
}


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
    stop("No effects with a computable F statistic found in the aovlist.")
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
