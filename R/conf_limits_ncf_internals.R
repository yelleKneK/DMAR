# Internal helpers shared by the ci_* functions that build their confidence
# intervals by inverting a noncentral F-distribution through
# conf_limits_ncf(): ci_snr(), ci_srsnr(), ci_pvaf(), ci_R2() with fixed
# predictors, ci_eta_squared(), ci_eta_squared_partial(),
# ci_omega_squared(), ci_eta_squared_generalized() with the parametric
# method, and ci_mahalanobis().
#
# Two jobs (QC decision 12, 2026-08-12):
#
# * When conf_limits_ncf() sets the lower noncentrality limit to 0 because
#   the observed F is below the alpha_lower critical value of the central
#   F-distribution (a normal consequence of a small observed effect, not a
#   failure), the warning the user sees must state the consequence for the
#   interval the user asked for, not conf_limits_ncf()'s internals.
#   .conf_limits_ncf_for() muffles the clamp warning (condition class
#   "dmar_ncf_clamp") and re-signals it with the calling function's own
#   quantity in the message. The re-signaled warning carries the same
#   condition class, so the iterative callers that deduplicate the clamp
#   (ss_aipe_R2(), ss_aipe_omega_squared(), factorial_anova(),
#   simple_effects_AB()) can keep matching it by class.
#
# * When the inner root finding genuinely fails, the error must name the
#   function the user called, the argument values that caused it, and what
#   to try, rather than surfacing a bare uniroot() message.

.conf_limits_ncf_for <- function(caller, quantity, F_value, df_1, df_2,
                                 conf_level = NULL, alpha_lower = NULL,
                                 alpha_upper = NULL, ...) {
  fmt <- function(x) if (is.null(x)) "NULL" else format(x)
  withCallingHandlers(
    tryCatch(
      conf_limits_ncf(F_value = F_value, conf_level = conf_level,
                      df_1 = df_1, df_2 = df_2,
                      alpha_lower = alpha_lower, alpha_upper = alpha_upper,
                      ...),
      error = function(e) {
        stop(sprintf(
          paste0("In %s(), the confidence limits for the noncentrality ",
                 "parameter could not be computed for the observed ",
                 "F-statistic %s with %s numerator and %s denominator ",
                 "degrees of freedom (alpha_lower = %s, alpha_upper = %s). ",
                 "The inner computation reported: %s"),
          caller, fmt(F_value), fmt(df_1), fmt(df_2),
          fmt(alpha_lower), fmt(alpha_upper), conditionMessage(e)
        ), call. = FALSE)
      }
    ),
    dmar_ncf_clamp = function(w) {
      warning(warningCondition(sprintf(
        paste0("The observed F_value is below the alpha_lower critical ",
               "value of the central F-distribution, so the lower ",
               "confidence limit on %s is 0."),
        quantity), class = "dmar_ncf_clamp"))
      invokeRestart("muffleWarning")
    }
  )
}


# Restrict the clamp warning to a single report per user-level call. The
# model-based paths of ci_eta_squared(), ci_eta_squared_partial(),
# ci_omega_squared(), and ci_eta_squared_generalized() compute one interval
# per effect; each clamped effect signals the classed clamp warning, and
# without deduplication a small factorial fit would print several identical
# messages. The first clamp warning propagates untouched; later ones are
# muffled.
.warn_ncf_clamp_once <- function(expr) {
  clamp_seen <- FALSE
  withCallingHandlers(expr, dmar_ncf_clamp = function(w) {
    if (clamp_seen) invokeRestart("muffleWarning")
    clamp_seen <<- TRUE
  })
}
