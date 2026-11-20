# Classed refusals. A refusal DMAR signals on purpose (an analysis it will
# not run as asked, with the reason and what to do instead) carries a class
# a caller can catch by name, so that a program built on DMAR (SIEVER, a
# simulation runner) reads the kind of refusal from the condition rather
# than from its wording, which may change. Every refusal made through this
# helper inherits "dmar_refusal"; its own class names the kind, and the
# same name is stored in the condition's `code` element for a caller that
# prefers a string. See ?dmar_metadata for the classes in use.
#
# .dmar_refuse("Fit indices ... are unavailable.", "h1_unavailable")
# signals an error of class c("dmar_h1_unavailable", "dmar_refusal",
# "error", "condition") with code "h1_unavailable". Extra named arguments
# are stored in the condition (a diagnostic, the rows at issue).
# `extra_class` names classes the refusal also carries, after its own and
# before dmar_refusal, for a refusal that is also one of an older kind a
# caller already reads by class (fa_exploratory()'s refusal of a fit whose
# likelihood has no maximum is also a dmar_efa_no_solution).
.dmar_refuse <- function(message, code, ..., extra_class = NULL) {
  stop(structure(list(message = paste0(message, collapse = ""), call = NULL,
                      code = code, ...),
                 class = c(paste0("dmar_", code), extra_class,
                           "dmar_refusal", "error", "condition")))
}

# A classed warning, for a qualification a caller may want to muffle or
# read by class rather than by its wording: class c("dmar_<code>",
# "warning", "condition") with the code stored in the condition.
.dmar_warn <- function(message, code, ...) {
  warning(structure(list(message = paste0(message, collapse = ""), call = NULL,
                         code = code, ...),
                    class = c(paste0("dmar_", code), "warning", "condition")))
}

# .dmar_check_count(x, arg, min = 1L, max = NULL)
#
# The one check of a count of draws: resamples, permutations, bootstrap or
# Monte Carlo replications. `x` is the value given, `arg` the argument's
# name for the message, `min` the smallest admissible count, and `max` an
# optional largest one. The count must be a single finite whole number of
# at least `min` and at most `max`; with or without `max` it may not exceed
# .Machine$integer.max, the largest count R holds as an integer. A whole
# double (1000) and an integer (1000L) are both accepted. NULL, NA, NaN,
# Inf, a logical, a character string, a vector of any length but 1, a
# fractional value (2.5), and a value outside the range are refused. The
# comparison with round(x) is exact: .dmar_is_whole() (R/dmar_tbl.R) is a
# display helper with a relative tolerance and is not used here.
#
# Returns as.integer(x), visibly, so that a caller writes
# B <- .dmar_check_count(B, "B") and loops or allocates on an integer.
# Call it before any allocation or loop that the count sizes.
#
# The refusal is made through .dmar_refuse() with code "count_invalid": an
# error of class c("dmar_count_invalid", "dmar_refusal", "error",
# "condition") that carries `arg`, `min`, and `max` (NULL when no maximum
# was stated). Its message has the wording of the count checks DMAR already
# has: "'B' must be a single whole number of at least 2." and, when a
# maximum was stated or the value is above the largest integer,
# "'B' must be a single whole number of at least 2 and at most 100000."
#
# `reason`, an optional string, is appended to the message after one space
# and stored in the condition as `reason`: why the minimum is what it is
# (the bootstrap rule of .dmar_check_B() in R/dmar_bootstrap.R states the
# count behind its minimum this way). The class and the fields arg, min,
# and max are the same with or without it.
#
# Not for degrees of freedom, sample sizes that may be fractional, or any
# other quantity that is legitimately not whole (wald_df, df_error, the n
# of ci_cv()).
.dmar_check_count <- function(x, arg, min = 1L, max = NULL, reason = NULL) {
  upper <- if (is.null(max)) .Machine$integer.max else
    base::min(max, .Machine$integer.max)
  shape_ok <- is.numeric(x) && length(x) == 1L && is.finite(x)
  if (shape_ok && x == round(x) && x >= min && x <= upper)
    return(as.integer(x))
  state_upper <- !is.null(max) || (shape_ok && x > upper)
  shown <- function(v) format(v, scientific = FALSE, trim = TRUE)
  .dmar_refuse(c("'", arg, "' must be a single whole number of at least ",
                 shown(min),
                 if (state_upper) c(" and at most ", shown(upper)), ".",
                 if (!is.null(reason)) c(" ", reason)),
               "count_invalid", arg = arg, min = min, max = max,
               reason = reason)
}

# .dmar_check_desired_power(desired_power)
#
# The one check of `desired_power` in the power planners (ss_power_*):
# a single number strictly between 0 and 1. Anything else is refused with
# a condition of class dmar_desired_power_invalid that carries the value
# given. Both ends are excluded in every planner: no finite sample size
# has power 1, and a power of 0 does not define a sample size. Through
# 1.0.x ss_power_R2() accepted 0 and 1 (desired_power = 1 returned 346
# for true_R2 = .3 and p = 5).
.dmar_check_desired_power <- function(desired_power) {
  if (is.numeric(desired_power) && length(desired_power) == 1L &&
      is.finite(desired_power) && desired_power > 0 && desired_power < 1) {
    return(invisible(desired_power))
  }
  given <- if (is.numeric(desired_power) && length(desired_power) == 1L) {
    paste0(format(desired_power), " was given")
  } else if (is.null(desired_power)) {
    "NULL was given"
  } else {
    paste0("an object of class ", class(desired_power)[1L], " and length ",
           length(desired_power), " was given")
  }
  .dmar_refuse(c(
    "'desired_power' must be a single number in (0, 1), strictly between ",
    "0 and 1; ", given, ". It is the probability with which the planned ",
    "study should reject the null hypothesis when the planning values ",
    "are the population values. No finite sample size has power 1, and a ",
    "power of 0 does not define a sample size, so neither end is ",
    "accepted. Give a value such as desired_power = .80 or .90; a value ",
    "as close to 1 as wanted, such as .999, is accepted."),
    "desired_power_invalid", arg = "desired_power", value = desired_power)
}

# .dmar_refuse_argument(message, arg, ...)
#
# A refusal of the value given for an argument (or of a combination of
# arguments) that has no class of its own: an error of class
# c("dmar_argument_invalid", "dmar_refusal", "error", "condition") that
# carries the name or names at issue as `arg`. Further named values (the
# value given, the bound it met) are stored in the condition. The message
# says what is wrong, why, and what to give instead. A check with a class
# of its own (.dmar_check_count(), .dmar_check_desired_power(),
# .dmar_check_assurance()) is used where one applies.
#
# When a call site passes arg = NULL and the message opens with the name
# or names of the arguments in single quotes ("'x' must be ...",
# "'a' and 'b' cannot both be given"), those names are stored as `arg`:
# many checks build their message from the name with sprintf() and do not
# repeat it. A message that opens with anything else leaves `arg` NULL (a
# property of the data, a search that ran out).
.dmar_refuse_argument <- function(message, arg, ...) {
  if (is.null(arg)) arg <- .dmar_leading_arg_names(message)
  .dmar_refuse(message, "argument_invalid", arg = arg, ...)
}

# The argument names a message opens with: one or more identifiers in
# single quotes, joined by commas, "and", or "or". NULL when the message
# does not open with one.
.dmar_leading_arg_names <- function(message) {
  msg <- paste0(message, collapse = "")
  id <- "'[A-Za-z.][A-Za-z0-9._]*'"
  lead <- regmatches(msg, regexpr(
    paste0("^", id, "((, and |, or |, | and | or )", id, ")*"), msg))
  if (!length(lead)) return(NULL)
  gsub("'", "", regmatches(lead, gregexpr(id, lead))[[1L]], fixed = TRUE)
}

# .dmar_check_assurance(assurance)
#
# The check of `assurance` in a planner whose procedure can be run for any
# probability strictly inside (0, 1): NULL (no assurance), or a single
# number strictly between 0 and 1. Anything else is refused with a
# condition of class dmar_assurance_invalid that carries the value given.
# Both ends are excluded: no finite sample size makes the interval no
# wider than the width with probability 1, and a probability of 0 asks
# for nothing.
.dmar_check_assurance <- function(assurance) {
  if (is.null(assurance)) return(invisible(NULL))
  if (is.numeric(assurance) && length(assurance) == 1L &&
      is.finite(assurance) && assurance > 0 && assurance < 1) {
    return(invisible(assurance))
  }
  given <- if (is.numeric(assurance) && length(assurance) == 1L) {
    paste0(format(assurance), " was given")
  } else {
    paste0("an object of class ", class(assurance)[1L], " and length ",
           length(assurance), " was given")
  }
  .dmar_refuse(c(
    "'assurance' must be NULL or a single number in (0, 1), strictly ",
    "between 0 and 1; ", given, ". It is the probability that the ",
    "confidence interval of the planned study is no wider than 'width'. ",
    "No finite sample size makes that probability 1, and a probability of ",
    "0 asks for nothing, so neither end is accepted. Give a value such as ",
    "assurance = .80 or .99, or leave it NULL to plan for the expected ",
    "width."),
    "assurance_invalid", arg = "assurance", value = assurance)
}

# .dmar_check_width(width, what)
#
# The check of `width`, the desired full width of the confidence interval,
# where a sample size is planned from it: a single positive finite number.
# NULL or anything else is refused with a condition of class
# dmar_argument_invalid (arg = "width"). `what` names the quantity whose
# interval is planned, for the message.
.dmar_check_width <- function(width, what = "the parameter") {
  if (is.numeric(width) && length(width) == 1L && is.finite(width) &&
      width > 0) {
    return(invisible(width))
  }
  given <- if (is.null(width)) {
    "none was given"
  } else if (is.numeric(width) && length(width) == 1L) {
    paste0(format(width), " was given")
  } else {
    paste0("an object of class ", class(width)[1L], " and length ",
           length(width), " was given")
  }
  .dmar_refuse_argument(c(
    "'width' must be a single positive number, the desired full width of ",
    "the confidence interval for ", what, "; ", given, ". The sample size ",
    "is planned from it, so it cannot be left out here. Give the width in ",
    "the units of the estimate, for example width = .30."),
    arg = "width", value = width)
}

# .dmar_check_lavaan_dots(..., fn)
#
# The check of the `...` of a function that passes it to lavaan::sem().
# Every name given there must be an argument of lavaan::lavaan() or a
# lavaan option (names(lavaan::lavOptions())); a name that is neither, or
# an argument with no name, is refused before any model is fitted, with a
# condition of class dmar_argument_invalid whose `arg` holds the names as
# they were typed. Without it a stray name (a DMAR 1.0.x name such as
# `Sigma`, `desired_width`, or `N`, or a misspelling) reached lavaan,
# which stopped with its own message and printed `N` as 'n'. It is the
# first line of the function, so that a 1.0.x name is reported before the
# argument it stood for is found missing. `fn` is the calling function's
# name, for the message.
.dmar_check_lavaan_dots <- function(..., fn) {
  given <- names(list(...))
  if (...length() == 0L) return(invisible(NULL))
  # Without lavaan the function's own refusal (package_required) follows.
  if (!requireNamespace("lavaan", quietly = TRUE)) return(invisible(NULL))
  if (is.null(given)) given <- rep("", ...length())
  known <- union(names(formals(lavaan::lavaan)), names(lavaan::lavOptions()))
  if (any(!nzchar(given))) {
    .dmar_refuse_argument(c(
      "An argument without a name was given to ", fn, "() beyond the ",
      "arguments it has. Further arguments are passed to lavaan::sem() by ",
      "name, so each must be named; give every argument of ", fn, "() by ",
      "its name (see ?", fn, ")."),
      arg = NULL, fn = fn)
  }
  bad <- setdiff(given, known)
  if (length(bad)) {
    shown <- paste0("'", bad, "'", collapse = ", ")
    .dmar_refuse_argument(c(
      shown, if (length(bad) == 1L) " is not an argument" else
        " are not arguments", " of ", fn, "() or of lavaan::sem(), to ",
      "which ", fn, "() passes any further argument. Nothing was fitted. ",
      "Check the spelling against ?", fn, " and ?lavaan::lavOptions."),
      arg = bad, fn = fn)
  }
  invisible(NULL)
}

# How a value that failed a check is described in the message: "NULL was
# given", "NA was given", "0.5 was given", or the class and length of
# anything that is not one number.
.dmar_given <- function(x) {
  if (is.null(x)) return("NULL was given")
  if (is.atomic(x) && length(x) == 1L && is.na(x)) {
    return(paste0(if (is.nan(x)) "NaN" else "NA", " was given"))
  }
  if (is.numeric(x) && length(x) == 1L) {
    # Seven significant digits, unless they show a different number (a
    # value of 1 + 1e-12 refused for not being whole would read "1").
    txt <- format(x)
    if (is.finite(x) && as.numeric(txt) != x) txt <- format(x, digits = 15)
    return(paste0(txt, " was given"))
  }
  if (is.character(x) && length(x) == 1L) {
    return(paste0("the character string \"", x, "\" was given"))
  }
  paste0("an object of class ", class(x)[1L], " and length ", length(x),
         " was given")
}

# .dmar_check_conf_level(conf_level, null_ok = FALSE)
#
# The one check of `conf_level` on entry: a single number strictly between
# 0 and 1. With null_ok = TRUE, NULL passes (a function that takes
# `alpha_lower` and `alpha_upper` in its place, or that reads the level
# from elsewhere when none is given). Anything else is refused with a
# condition of class dmar_argument_invalid (arg = "conf_level") that
# carries the value given. Both ends are excluded: an interval with
# confidence 1 has no finite limits and one with confidence 0 is a point.
# A logical, a character string, NA, and a vector of two levels are
# refused here, before the level reaches a quantile function that would
# stop with R's own message or return limits of NaN or Inf.
.dmar_check_conf_level <- function(conf_level, null_ok = FALSE) {
  if (is.null(conf_level) && null_ok) return(invisible(NULL))
  if (is.numeric(conf_level) && length(conf_level) == 1L &&
      is.finite(conf_level) && conf_level > 0 && conf_level < 1) {
    return(invisible(conf_level))
  }
  .dmar_refuse_argument(c(
    "'conf_level' must be a single number in (0, 1), strictly between 0 ",
    "and 1; ", .dmar_given(conf_level), ". It is the confidence level of ",
    "the interval, as a proportion and not a percentage. An interval with ",
    "confidence 1 has no finite limits and one with confidence 0 is a ",
    "single point, so neither end is accepted. Give a value such as ",
    "conf_level = .95 or .90."),
    arg = "conf_level", value = conf_level)
}

# .dmar_check_alpha_level(alpha_level)
#
# The one check of `alpha_level` on entry: a single number strictly
# between 0 and 1, the Type I error rate of the test. Anything else is
# refused with a condition of class dmar_argument_invalid
# (arg = "alpha_level") that carries the value given. A function with a
# narrower range of its own (an equivalence test that needs a level below
# .5) checks that after this.
.dmar_check_alpha_level <- function(alpha_level) {
  if (is.numeric(alpha_level) && length(alpha_level) == 1L &&
      is.finite(alpha_level) && alpha_level > 0 && alpha_level < 1) {
    return(invisible(alpha_level))
  }
  .dmar_refuse_argument(c(
    "'alpha_level' must be a single number in (0, 1), strictly between 0 ",
    "and 1; ", .dmar_given(alpha_level), ". It is the Type I error rate ",
    "of the test, as a proportion. A test at level 0 never rejects and ",
    "one at level 1 always does, so no power or sample size is defined ",
    "at either end. Give a value such as alpha_level = .05."),
    arg = "alpha_level", value = alpha_level)
}

# .dmar_min_tail_area
#
# The smallest nonzero tail area (alpha_lower or alpha_upper) that
# ci_nc_t(), and so every interval and plan whose limits come from it,
# accepts. The central intervals of the contrast, the correlation, and the
# regression coefficient (ci_c(), ci_r(), and ci_reg_coef() with
# noncentral = FALSE, so ci_rc()) and their planners (ss_aipe_c(),
# ss_aipe_r(), and ss_aipe_reg_coef()) keep it as well, so that an interval
# and the plan for it accept the same tail areas. ci_c_ancova() and
# ss_aipe_c_ancova() do not; they accept any level, as a pair.
#
# Below the floor the noncentral t limits lose accuracy. A limit is the
# root of a tail probability, and R's pt() sums its series to an absolute
# tolerance of 1e-12 and forms one of the two tails as 1 minus the other
# (.nc_t_pt() in R/nc_t_internals.R), so the absolute error of the
# probability is a growing share of a smaller tail. Against the oracle of
# tests/testthat/helper-nc_t_oracle.R, which conditions on the normal
# variate, the limits of ci_nc_t() with equal tails at t = 5 on 47 df, .5
# on 10, 2 on 200, 10 on 20, -3 on 30, and 1 on 1,000, at 66 tail areas
# from 1e-6 down to 1e-9, were off by at most 1.7e-5 at tails of 1e-8 and
# above, 7.3e-5 at tails of 2e-9 and above, and 1.4e-4 at 1e-9 (the upper
# limit at t = 10 on 20 df, 6.0e-6 of the width of the interval). The
# error does not fall steadily as the tail grows: at t = 5 on 47 df the
# lower limit was off by 2.3e-5 at 2.5e-9 and 3.6e-6 at 2e-9. At 1e-10 the
# largest error of the six cases was 1.1e-3. Against an oracle that
# integrates over the chi variate, the upper limit of
# ci_reg_coef(noncentral = TRUE) at t = 5 on 47 df was off by 6.0e-7 of
# itself at 1e-9 and 7.6 percent at 1e-20; without the floor, the plan of
# ss_aipe_sc_ancova() for width .5 (delta .5) at tails of 1e-20 was 18,143
# per group where the oracle puts it between 2,800 and 3,000 (the oracle
# width is .503 at 2,800 and .486 at 3,000).
#
# A test of the display of the confidence level asks ci_smd() for
# conf_level = .999999995, tails of 2.5e-9, which is why the floor is not
# 1e-8. At 1e-9 the largest level accepted, .999999998, leaves tails of
# 1.0000000272e-9 in double precision, just above the floor; at a floor of
# 2e-9 the level .999999996 would leave 1.99999999895e-9, just below it. A
# central quantile taken in the upper tail (lower.tail = FALSE) is accurate
# at any tail; the central intervals keep the floor for the reason above.
.dmar_min_tail_area <- 1e-9

# A tail area as a refusal states it: four significant digits, and an
# exponent without the leading zero format() gives it (1e-9, 2.5e-10).
.dmar_tail_text <- function(x)
  sub("e([-+])0*([0-9])", "e\\1\\2", format(x, digits = 4))

# .dmar_check_tail_floor(alpha_lower, alpha_upper, conf_level = NULL,
#                        zero_ok = TRUE)
#
# The check of the floor, called once a function has resolved its two tail
# areas (from conf_level or as given) and checked them by its own rules. A
# tail area that is a single finite number above 0 and below
# .dmar_min_tail_area is refused with a condition of class
# dmar_argument_invalid whose `reason` is "tail_below_floor" and which
# carries the tail areas below the floor as `value` (named alpha_lower and
# alpha_upper) and the floor as `floor`. Anything else passes: a tail of 0
# (a one-sided interval, wherever the function allows one) and values the
# function's own checks deal with.
#
# `conf_level` is the level the tails were formed from, (1 - conf_level) / 2
# each, or NULL when the tails were given; with a level the refusal names
# 'conf_level' as its `arg` and says the largest level accepted, and
# without one it names the tail or tails below the floor. `zero_ok` says
# whether the function accepts a tail of 0, for the advice.
.dmar_check_tail_floor <- function(alpha_lower, alpha_upper,
                                   conf_level = NULL, zero_ok = TRUE) {
  floor <- .dmar_min_tail_area
  below <- function(a) is.numeric(a) && length(a) == 1L && is.finite(a) &&
    a > 0 && a < floor
  low <- c(alpha_lower = below(alpha_lower), alpha_upper = below(alpha_upper))
  if (!any(low)) return(invisible(NULL))
  shown <- .dmar_tail_text
  # A list keeps a NULL tail in its place, so `low` picks the right one
  # when the other was not given (c() would drop the NULL).
  value <- unlist(list(alpha_lower = alpha_lower, alpha_upper = alpha_upper)[low])
  why <- c(
    "Below ", shown(floor), " the noncentral t limits lose accuracy (the ",
    "tail probability a limit is found from carries an absolute error of ",
    "up to about 1e-12, a growing share of a smaller tail), so a limit or ",
    "a sample size could be off without notice. The central intervals of ",
    "a contrast, a correlation, and a regression coefficient, and their ",
    "plans, keep the same floor, so that an interval and the plan for it ",
    "accept the same tail areas.")
  if (!is.null(conf_level)) {
    .dmar_refuse_argument(c(
      "'conf_level' of ", format(conf_level, digits = 15), " leaves a tail ",
      "area of ", shown(alpha_lower), " on each side, below ", shown(floor),
      ", the smallest nonzero tail area accepted. ", why, " Give a ",
      "conf_level of at most ", format(1 - 2 * floor, digits = 15), "."),
      arg = "conf_level", reason = "tail_below_floor", value = value,
      floor = floor, conf_level = conf_level)
  }
  names_low <- names(low)[low]
  .dmar_refuse_argument(c(
    if (length(names_low) == 1L) {
      c("'", names_low, "' is ", shown(value), ", a tail area")
    } else if (value[1L] == value[2L]) {
      c("'alpha_lower' and 'alpha_upper' are both ", shown(value[1L]),
        ", tail areas")
    } else {
      c("'alpha_lower' and 'alpha_upper' are ", shown(value[1L]), " and ",
        shown(value[2L]), ", tail areas")
    },
    " above 0 and below ", shown(floor), ", the smallest nonzero tail area ",
    "accepted. ", why, " Give a tail area of at least ", shown(floor),
    if (zero_ok) ", or 0 for a one-sided interval", "."),
    arg = names_low, reason = "tail_below_floor", value = value,
    floor = floor)
}

# .dmar_check_assurance_tail_floor(assurance)
#
# The published assurance step of ss_aipe_sc_ancova() and ss_aipe_sm()
# (Kelley & Rausch, 2006) finds noncentral t limits with tail areas of
# (1 - assurance) / 2 and 1 - assurance, which ci_nc_t() holds to the floor
# of .dmar_check_tail_floor(). An assurance above
# 1 - 2 * .dmar_min_tail_area would be refused there in the names of tail
# areas the caller did not give, so it is refused here, in its own name:
# class dmar_argument_invalid, arg = "assurance", reason
# "tail_below_floor". ss_aipe_sc() and ss_aipe_smd() apply the same bound,
# so that the planners of a standardized contrast and of a standardized
# mean difference accept the same assurances. NULL, and anything the
# planner's own check of `assurance` deals with, passes.
.dmar_check_assurance_tail_floor <- function(assurance) {
  floor <- .dmar_min_tail_area
  if (!(is.numeric(assurance) && length(assurance) == 1L &&
        is.finite(assurance) && assurance < 1 &&
        (1 - assurance) / 2 < floor)) return(invisible(NULL))
  .dmar_refuse_argument(c(
    "'assurance' of ", format(assurance, digits = 15), " is too close to 1: ",
    "(1 - assurance) / 2 is below ", .dmar_tail_text(floor), ", the ",
    "smallest nonzero tail area accepted. The assurance step of Kelley and ",
    "Rausch (2006) finds noncentral t limits with tail areas of ",
    "(1 - assurance) / 2 and 1 - assurance, and below the floor those ",
    "limits lose accuracy. Give an assurance of at most ",
    format(1 - 2 * floor, digits = 15), "."),
    arg = "assurance", reason = "tail_below_floor", value = assurance,
    floor = floor)
}

# .dmar_check_numbers(fn, args, required = character(0))
#
# The check, on entry, of arguments that each hold one number: `args`
# names them, and the values are read from the calling function's frame.
# An argument the caller left out is not looked at (its default, or the
# function's own message for a missing argument, stands). An argument
# given as NULL passes unless it is named in `required`. Otherwise the
# value must not be a single NA or NaN, a character string, a list, or a
# vector of length 0: each of those reached arithmetic or an `if` and
# stopped with R's own message ("missing value where TRUE/FALSE needed",
# "non-numeric argument to binary operator", "argument is of length
# zero"), which names neither the argument nor the function. The refusal
# has class dmar_argument_invalid, names the argument as `arg`, and
# carries the value given. Range, wholeness, and length above 1 are left
# to the function's own checks, which follow this one. `fn` is the
# calling function's name, for the message.
.dmar_check_numbers <- function(fn, args, required = character(0)) {
  frame <- parent.frame()
  for (nm in args) {
    if (eval(call("missing", as.name(nm)), frame)) next
    x <- get(nm, envir = frame, inherits = FALSE)
    bad <- if (is.null(x)) {
      nm %in% required
    } else {
      is.character(x) || is.list(x) || length(x) == 0L ||
        (is.atomic(x) && length(x) == 1L && is.na(x))
    }
    if (bad) {
      .dmar_refuse_argument(c(
        "'", nm, "' must be a number that is not missing; ",
        .dmar_given(x), ". ", fn, "() computes with this value, and a ",
        "missing value, an empty vector, a character string, or a list ",
        "cannot stand in for it. Give the number itself (see ?", fn,
        " for what '", nm, "' is)."),
        arg = nm, value = x, fn = fn)
    }
  }
  invisible(NULL)
}

# .dmar_refuse_no_sample_size(fn, desired_power, size = "sample size")
#
# The refusal of a power planner whose search reached its limit of
# 1,000,000 without the power reaching `desired_power`: class
# dmar_computation_failed, with the function and the power asked for.
# The message says what that means (the effect at the planning values is
# too small for any size the search covers, or is 0) and what to check.
.dmar_refuse_no_sample_size <- function(fn, desired_power, size = "sample size") {
  .dmar_refuse(c(
    fn, "() found no ", size, " up to 1,000,000 at which the power ",
    "reaches desired_power = ", format(desired_power), ". The effect at ",
    "the planning values is too small for the test to reach that power ",
    "at any size the search covers; an effect of exactly 0 has power ",
    "equal to 'alpha_level' at every sample size. Check the planning ",
    "values (the effect and the standard deviation or variance it is ",
    "set against), or ask for less power."),
    "computation_failed", fn = fn, desired_power = desired_power)
}

# .dmar_refuse_zero_effect(fn, arg, what)
#
# The refusal of a planning value of exactly 0 for the effect a power
# planner plans to detect, when a sample size is to be solved for: the
# null hypothesis is then true at the planning value, the power is
# `alpha_level` at every sample size, and no sample size is defined.
# Class dmar_argument_invalid, arg = the planning value's name.
.dmar_refuse_zero_effect <- function(fn, arg, what) {
  .dmar_refuse_argument(c(
    "'", arg, "' is 0, so ", fn, "() has no sample size to plan: ", what,
    " of 0 is the null hypothesis itself, and the power of the test is ",
    "then 'alpha_level' at every sample size and never reaches ",
    "'desired_power'. Give the smallest value that would matter if it ",
    "were the population value. (With a sample size given, the function ",
    "returns the power at 0, which is the Type I error rate.)"),
    arg = arg, value = 0, fn = fn)
}

# .dmar_check_assurance_single(assurance)
#
# The entry check of `assurance` in a planner whose own rule for the
# range is narrower than (0, 1) (a procedure defined only for a
# probability above one half): NULL, or a single number that is not
# missing. The range is the function's own check, which follows and
# states it. Refused here, with class dmar_assurance_invalid, are NA, a
# character string, a logical, a list, and a vector of any length but 1,
# which otherwise reached a comparison and stopped with R's own message.
.dmar_check_assurance_single <- function(assurance) {
  if (is.null(assurance)) return(invisible(NULL))
  if (is.numeric(assurance) && length(assurance) == 1L && !is.na(assurance)) {
    return(invisible(assurance))
  }
  .dmar_refuse(c(
    "'assurance' must be NULL or a single number; ", .dmar_given(assurance),
    ". It is the probability that the confidence interval of the planned ",
    "study is no wider than 'width'. Give a value such as assurance = .80 ",
    "or .99, or leave it NULL to plan for the expected width."),
    "assurance_invalid", arg = "assurance", value = assurance)
}

# .dmar_check_contrast_nonzero(c_weights, arg = "c_weights", margin = NULL)
#
# The check that contrast weights define a contrast. Weights with no
# nonzero value (c(0, 0), numeric(0), or c(FALSE, FALSE), whose values are
# zeros when R computes with them) pass the rules that the weights sum to
# zero and that the positive ones sum to at most 1, and they define no
# contrast: DMAR 1.0.0 gave an interval of zero width for them ([0, 0]
# from ci_c() with means), a plan of 2 per group (ss_aipe_c()), and a
# power of 1 with an infinite noncentrality parameter (ss_power_c()), or
# stopped with a message from R (ss_aipe_sc_ancova_sensitivity()). Every
# function that takes contrast weights calls this check, and such weights
# are refused with a condition of class dmar_argument_invalid that carries
# `arg`, the name of the argument as the calling function spells it
# ("c_weights", "contrast", "weights", "contrasts"). The message is
# "'c_weights' has no nonzero weight, so it defines no contrast; give
# weights such as c(1, -1) or c(1/2, 1/2, -1).", with that name in place
# of 'c_weights'.
#
# `margin` is for a matrix of several contrasts: 1 when each row is a
# contrast (ci_c_ancova_bp(), contrast_test()), 2 when each column is
# (ci_scheffe(), is_orthogonal_set()). Each contrast is checked, and the
# message names those with no nonzero weight, by their row or column names
# when the matrix has them and by position otherwise. A matrix that holds
# one contrast is checked as a vector is.
#
# Logical and complex weights are checked as numbers are, as the checks
# this replaced in ss_aipe_sc_ancova() and ci_sc_ancova() checked them. A
# value of another type, and weights with a missing value and no nonzero
# one, are passed unchanged to the calling function's own checks. Call it
# before the weights enter any computation. Returns the weights
# invisibly.
.dmar_check_contrast_nonzero <- function(c_weights, arg = "c_weights",
                                         margin = NULL) {
  if (is.data.frame(c_weights)) c_weights <- as.matrix(c_weights)
  if (!(is.numeric(c_weights) || is.logical(c_weights) ||
        is.complex(c_weights))) return(invisible(c_weights))
  none <- function(w) isTRUE(!any(w != 0))
  advice <- "give weights such as c(1, -1) or c(1/2, 1/2, -1)."
  d <- dim(c_weights)
  if (is.null(margin) || is.null(d) || any(d == 0L) || d[margin] == 1L) {
    if (none(c_weights)) .dmar_refuse_argument(c(
      "'", arg, "' has no nonzero weight, so it defines no contrast; ",
      advice), arg = arg)
    return(invisible(c_weights))
  }
  bad <- which(apply(c_weights, margin, none))
  if (!length(bad)) return(invisible(c_weights))
  labels <- dimnames(c_weights)[[margin]]
  named <- !is.null(labels) && all(nzchar(labels[bad]))
  ids <- if (named) paste0("'", labels[bad], "'") else as.character(bad)
  one <- length(ids) == 1L
  if (!one) ids <- paste0(paste(ids[-length(ids)], collapse = ", "),
                          if (length(ids) > 2L) ",", " and ",
                          ids[length(ids)])
  what <- if (named) {
    c(if (one) "The contrast " else "The contrasts ", ids, " in '", arg, "'")
  } else {
    c(if (margin == 1L) "Row" else "Column", if (!one) "s", " ", ids,
      " of '", arg, "'")
  }
  .dmar_refuse_argument(c(
    what, if (one) " has" else " have", " no nonzero weight, so ",
    if (one) "it defines" else "they define", " no contrast; ", advice),
    arg = arg)
}
