# Shared helpers for the regression coefficient family: ci_reg_coef() and
# its wrappers ci_rc() and ci_src(), and ss_aipe_reg_coef() and its
# wrappers ss_aipe_rc() and ss_aipe_src().


# .reg_coef_in_wrapper_names(expr, map)
#
# Evaluates `expr`, a call of the general engine (ci_reg_coef() or
# ss_aipe_reg_coef()) that a wrapper makes with its own arguments passed
# under the engine's names, and states the engine's conditions in the
# wrapper's names. `map` is a named character vector, engine name =
# wrapper name: ci_src() passes c(reg_coef_j = "std_reg_coef_j",
# se_reg_coef_j = "se_std_reg_coef_j"), because the engine's reg_coef_j is
# the wrapper's std_reg_coef_j. A refusal (class dmar_refusal) is signaled
# again with every whole-word occurrence of an engine name in its message,
# and every element of its `arg` field, replaced by the wrapper's name; its
# class and other fields are kept. A classed warning (one made by
# .dmar_warn(), which carries a `code`) is signaled again with its message
# rewritten the same way, and the original is muffled. Other conditions,
# R's own included, pass through untouched. Through the first 1.1.0 builds
# ss_aipe_src() and ci_src() reported, for example, "'true_reg_coef_j' must
# also be specified" to a user whose argument is true_std_reg_coef_j.
.reg_coef_in_wrapper_names <- function(expr, map) {
  rename <- function(x) {
    if (is.null(x) || !is.character(x)) return(x)
    for (inner in names(map)) {
      x <- gsub(paste0("(?<![A-Za-z0-9._])", inner, "(?![A-Za-z0-9._])"),
                map[[inner]], x, perl = TRUE)
    }
    x
  }
  withCallingHandlers(
    tryCatch(expr, dmar_refusal = function(e) {
      e$message <- rename(e$message)
      if (!is.null(e$arg)) e$arg <- rename(e$arg)
      stop(e)
    }),
    warning = function(w) {
      if (is.null(w$code)) return(invisible(NULL))
      w$message <- rename(w$message)
      warning(w)
      invokeRestart("muffleWarning")
    })
}
