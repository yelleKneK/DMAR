# Shared normalization of the 'alternative' argument for the directional
# test family (ci_dunnett, power_fisher_exact, randomization_test,
# randomization_test_paired, summary_t_test, welch_t).
#
# The package vocabulary is snake_case, so the canonical values are
# "two_sided", "less", and "greater". The base-R spelling "two.sided" is
# accepted as an alias, so a call written against t.test() muscle memory
# keeps working, but the normalized value (and anything stored on the
# returned object) is always the underscore form. The cv_* critical value
# family has its own, deliberately wider synonym vocabulary ("ne", "2s",
# "!=", ...) with "not_equal" canonical; it does not route through this
# helper.
.match_alternative <- function(alternative) {
  if (!is.character(alternative) || length(alternative) < 1L) {
    stop("'alternative' must be a character string.", call. = FALSE)
  }
  alternative <- alternative[1L]
  if (identical(alternative, "two.sided")) alternative <- "two_sided"
  match.arg(alternative, c("two_sided", "less", "greater"))
}
