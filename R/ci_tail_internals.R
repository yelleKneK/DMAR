# Resolve a confidence specification before a wrapper forwards its defaults.
# Explicitly consistent level and tail areas are accepted; an omitted level
# yields to supplied tail areas. A zero tail denotes a one-sided bound.
.ci_resolve_tails <- function(conf_level, alpha_lower, alpha_upper,
                              level_missing = FALSE, max_tail = 1) {
  scalar <- function(x) is.numeric(x) && length(x) == 1L && is.finite(x)
  supplied <- !is.null(alpha_lower) || !is.null(alpha_upper)
  if (supplied) {
    if (!scalar(alpha_lower) || !scalar(alpha_upper) ||
        alpha_lower < 0 || alpha_upper < 0 ||
        alpha_lower > max_tail || alpha_upper > max_tail ||
        alpha_lower + alpha_upper <= 0 || alpha_lower + alpha_upper >= 1) {
      .dmar_refuse_argument(
        "Supply both 'alpha_lower' and 'alpha_upper' as finite tail areas with a sum strictly between zero and one.",
        arg = c("alpha_lower", "alpha_upper"))
    }
    coverage <- 1 - alpha_lower - alpha_upper
    if (!level_missing && !is.null(conf_level) &&
        (!scalar(conf_level) || conf_level <= 0 || conf_level >= 1 ||
         abs(conf_level - coverage) > sqrt(.Machine$double.eps))) {
      .dmar_refuse_argument(
        "The explicitly supplied 'conf_level' is inconsistent with 'alpha_lower' and 'alpha_upper'.",
        arg = c("conf_level", "alpha_lower", "alpha_upper"))
    }
    return(list(conf_level = coverage, alpha_lower = alpha_lower, alpha_upper = alpha_upper))
  }
  if (!scalar(conf_level) || conf_level <= 0 || conf_level >= 1) {
    .dmar_refuse_argument(
      "Supply a finite 'conf_level' strictly between zero and one, or both tail areas.",
      arg = "conf_level")
  }
  list(conf_level = conf_level, alpha_lower = (1 - conf_level) / 2,
    alpha_upper = (1 - conf_level) / 2)
}
