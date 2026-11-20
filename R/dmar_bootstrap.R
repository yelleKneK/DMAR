# One rule for the number of bootstrap and Monte Carlo replications.
#
# Internal helpers, shared by every function that reads an interval, a
# bound, a reference quantile, a standard error, or a Monte Carlo p-value
# from B replications, and the help topic ?dmar_bootstrap that states the
# rule to the user.
#
#   .dmar_boot_required()     the fewest replications a request may ask
#                             for
#   .dmar_boot_count_met()    the count rule itself: whether the number
#                             requested meets the requirement
#   .dmar_boot_exploratory()  whether options(dmar.bootstrap_exploratory =
#                             TRUE) is set
#   .dmar_check_B()           the check of the request, before anything is
#                             computed
#   .dmar_boot_settle()       the verdict on a finished run and the
#                             conditions it signals
#   .dmar_boot_bca_positions() the same count applied where a BCa interval
#                             finally reads its limits
#   .dmar_boot_signal()       signals the conditions of a verdict
#   .dmar_boot_attrs()        writes the counts on a result
#   .dmar_boot_footer()       the line a print method shows
#   .boot_caution()           whether the replications in question reach
#                             the tail area beyond a limit
#   .check_boot_B()           the reliability family's entry to
#                             .dmar_check_B()
#
# No published source is cited for the counts. The count of 10 beyond a
# limit, the floor of 100, the factor of 50 for BCa, and the 100 of a
# standard error and of a p-value were computed for DMAR (the records are
# in the development notes of 2026-11-10): the tail area cut off by the
# k-th order statistic of B draws is Beta(k, B - k + 1) for any continuous
# distribution, and the standard deviations quoted in ?dmar_bootstrap are
# from simulations of 20,000 runs per cell.


# .dmar_boot_required(kind, tail, floor)
#
# The fewest replications a request may ask for, an integer. The count is
# held to the number requested (.dmar_boot_count_met()); what becomes of
# a run in which some replications returned no value is decided afterward
# by the share rules of .dmar_boot_settle().
#   kind = "quantile": 10 beyond each limit and never fewer than `floor`;
#     `tail` is the smallest tail area any reported limit cuts off, and
#     for a limit at quantile position p it is min(p, 1 - p).
#   kind = "bca": 50 beyond each nominal limit, never fewer than `floor`.
#   kind = "se", "p": `floor`, 100 (a standard error taken from the
#     replications; a Monte Carlo p-value).
# NA_integer_ when the count exceeds the largest integer R holds.
# The small allowance keeps 10 / 0.025 at 400 when the tail was computed
# as (1 - 0.95) / 2, which is slightly above 0.025 in floating point.
#
# Tests that exercise the mechanics of a bootstrap with a small B replace
# this binding (local_small_B() in tests/testthat/helper-bootstrap.R).
.dmar_boot_required <- function(kind = c("quantile", "bca", "se", "p"),
                                tail = NULL, floor = 100L) {
  kind <- match.arg(kind)
  floor <- as.integer(floor)
  if (kind %in% c("se", "p")) return(floor)
  if (!is.numeric(tail) || length(tail) != 1L || !is.finite(tail) ||
      tail <= 0 || tail > 0.5) {
    .dmar_refuse(
      "internal: 'tail' must be a single number in (0, 0.5].", "internal_error")
  }
  need <- max(floor, ceiling((if (kind == "bca") 50 else 10) / tail - 1e-8))
  if (need > .Machine$integer.max) return(NA_integer_)
  as.integer(need)
}


# .dmar_boot_count_met(B, required)
#
# The count rule, in the one place that states it: TRUE when the number of
# replications REQUESTED, B, is at least the number required (or the
# exploratory option is set, or the requirement is not a number R can
# hold, which .dmar_check_B() refuses on its own). The rule is not applied
# to the replications that returned a value: a request that meets the
# count and then loses a few replications still reports its interval, the
# warning states how many were dropped, and the share rules decide the
# rest (more than half dropped: no interval; dropped replications that
# reach the tail area beyond a limit: the caution). Ruled 2026-11-10;
# through the first versions of 1.1.0 the count was applied to the usable
# replications, so one dropped replication at exactly the required B
# withheld every interval.
.dmar_boot_count_met <- function(B, required) {
  .dmar_boot_exploratory() || is.na(required) || B >= required
}


# .dmar_boot_beyond()
#
# The number of replications that must lie beyond each limit where
# a BCa interval finally reads it (.dmar_boot_bca_positions()). A binding
# of its own so that local_small_B() can set it to 0 for tests of the
# mechanics of a BCa interval at a small B.
.dmar_boot_beyond <- function() 10


# .dmar_boot_exploratory()
#
# TRUE when the user has set options(dmar.bootstrap_exploratory = TRUE),
# the one documented way to run fewer replications than a reported result
# needs (a first look, a help page example, a simulation over many data
# sets). Under it a B from 2 up is accepted and every result is marked.
.dmar_boot_exploratory <- function() {
  isTRUE(getOption("dmar.bootstrap_exploratory", FALSE))
}


# An ordinal for a percentile ("95th", "97.5th", "1st").
.dmar_boot_ordinal <- function(p) {
  v <- 100 * p
  shown <- format(v, digits = 6, trim = TRUE)
  if (abs(v - round(v)) > 1e-8) return(paste0(shown, "th"))
  v <- as.integer(round(v))
  suffix <- if (v %% 100L %in% 11:13) "th" else
    switch(as.character(v %% 10L), "1" = "st", "2" = "nd", "3" = "rd", "th")
  paste0(v, suffix)
}

# A tail area as a proportion for the arithmetic in a message ("0.025").
.dmar_boot_tail <- function(tail) format(tail, digits = 3, trim = TRUE)

# What a result of each kind is called in a message.
.dmar_boot_noun <- function(kind, what, conf_level) {
  level <- if (is.null(conf_level)) NULL else
    paste0(.format_conf_pct(conf_level), "% ")
  switch(kind,
    quantile = if (what == "bound") "bound" else
      if (what == "reference") "reference quantile" else
        paste0(level, "interval"),
    bca = paste0(level, "BCa interval"),
    se = "standard error",
    p = "Monte Carlo p-value")
}


# .dmar_boot_reason(kind, tail, conf_level, floor, what, position)
#
# The sentence a refused B carries after "'B' must be a single whole
# number of at least N.": why that number, and where to read how to run
# fewer.
.dmar_boot_reason <- function(kind, tail, conf_level, floor = 100L,
                              what = "interval", position = NULL) {
  see <- paste0(" To run fewer for exploration or a demonstration, see ",
                "?dmar_bootstrap.")
  if (kind == "se") {
    return(paste0("A standard error taken from the replications is ",
                  "reported only from a request of at least ", floor,
                  " of them.", see))
  }
  if (kind == "p") {
    return(paste0("A Monte Carlo p-value is reported only from at least ",
                  floor, " replications, at which its standard error is ",
                  "at most ", format(1 / (2 * sqrt(floor)), digits = 2),
                  ".", see))
  }
  factor <- if (kind == "bca") 50 else 10
  raw <- ceiling(factor / tail - 1e-8)
  count <- paste0("(", factor, " / ", .dmar_boot_tail(tail), " = ",
                  format(raw, scientific = FALSE, trim = TRUE), ")")
  never <- if (raw < floor)
    paste0(", and never from fewer than ", floor, " replications")
  # Two limits with tail areas of their own, or a caller that passes no
  # level: the sentence names the smaller tail area and no level.
  if (is.null(conf_level) && what == "interval") what <- "limits"
  if (kind == "bca" && what %in% c("limits", "bound")) {
    return(paste0(
      "A BCa interval is read from quantiles that the bias correction and ",
      "the acceleration move away from the nominal ones (the smaller ",
      "nominal tail area is ", .pct(tail), "), usually to a point with ",
      "fewer replications beyond it. At least 10 replications must lie ",
      "beyond each limit where it finally sits, so a BCa interval is run ",
      "only with 50 beyond each nominal limit ", count, never, ".", see))
  }
  if (kind == "bca") {
    return(paste0(
      "A ", .format_conf_pct(conf_level), "% BCa interval is read from ",
      "quantiles that the bias correction and the acceleration move away ",
      "from ", .pct(tail), " and ", .pct(1 - tail), ", usually to a point ",
      "with fewer replications beyond it. At least 10 replications must ",
      "lie beyond each limit where it finally sits, so a BCa interval is ",
      "run only with 50 beyond each nominal limit ", count, never, ".",
      see))
  }
  if (what == "reference") {
    if (is.null(position)) position <- 1 - tail
    return(paste0(
      "A reference quantile at the ", .dmar_boot_ordinal(position),
      " percentile of the reference data sets has ", .pct(tail), " of them ",
      "beyond it, and it is reported only when at least 10 lie beyond it ",
      count, never, ".", see))
  }
  if (what == "bound") {
    if (is.null(position)) position <- 1 - tail
    return(paste0(
      "A bound at the ", .dmar_boot_ordinal(position), " percentile of ",
      "the replications has ", .pct(tail), " of them beyond it, and it is ",
      "reported only when at least 10 replications lie beyond it ", count,
      never, ".", see))
  }
  if (what == "limits") {
    return(paste0(
      "The limit that cuts off the smaller tail area has ", .pct(tail),
      " of the replications beyond it, and an interval read from the ",
      "quantiles of the replications is reported only when at least 10 ",
      "replications lie beyond each limit ", count, never, ".", see))
  }
  paste0(
    "A ", .format_conf_pct(conf_level), "% interval read from the ",
    "quantiles of the replications has ", .pct(tail), " of them beyond ",
    "each limit, and an interval is reported only when at least 10 ",
    "replications lie beyond each limit ", count, never, ".", see)
}


# .dmar_check_B(B, kind, tail, conf_level, arg, zero_ok, floor, what,
#               position)
#
# The check of a requested number of replications, made before anything is
# computed and only when replications are requested. Returns as.integer(B).
#   kind, tail, floor : as .dmar_boot_required().
#   conf_level : the level, for the message (and the tail of kind "se").
#   zero_ok    : TRUE where B = 0 means "no bootstrap" (htmt(),
#                factor_count()); 0L is then returned.
#   what       : "interval" (two-sided, equal tails), "limits" (two limits
#                with tail areas of their own: the ratio functions and
#                ci_eta_squared_generalized()), or "bound" (a one-sided
#                bound), or "reference" (the reference quantile of a
#                parallel analysis); `position` is the quantile a bound
#                or a reference quantile sits at.
# A B below the requirement is the classed refusal dmar_count_invalid of
# .dmar_check_count(), with the reason appended. Under
# options(dmar.bootstrap_exploratory = TRUE) a B from 2 up is accepted;
# the run is then marked by .dmar_boot_settle().
.dmar_check_B <- function(B, kind = "quantile", tail = NULL,
                          conf_level = NULL, arg = "B", zero_ok = FALSE,
                          floor = 100L, what = "interval", position = NULL) {
  if (zero_ok && is.numeric(B) && length(B) == 1L && !is.na(B) && B == 0)
    return(0L)
  required <- .dmar_boot_required(kind, tail, floor)
  if (.dmar_boot_exploratory()) return(.dmar_check_count(B, arg, min = 2L))
  if (is.na(required)) {
    .dmar_refuse(c("'", arg, "' cannot serve this level: 10 replications ",
                   "beyond each limit take more replications than R can ",
                   "count. Choose a level further from 1."),
                 "count_invalid", arg = arg, min = NA_integer_, max = NULL)
  }
  reason <- .dmar_boot_reason(kind, tail, conf_level, floor, what, position)
  if (zero_ok) reason <- paste0(reason, " ", arg, " = 0 runs none.")
  .dmar_check_count(B, arg, min = required, reason = reason)
}


# .boot_caution(in_question, B, conf_level)
#
# Whether `in_question` replications of B (those set aside and the
# improper ones kept) reach the tail area beyond a limit of a
# conf_level interval, (1 - conf_level)/2 of B. The small allowance
# keeps an exact tie (25 of 1000 at 95 percent) from failing on the
# floating point value of (1 - 0.95)/2, which is slightly above 0.025.
.boot_caution <- function(in_question, B, conf_level) {
  in_question > 0L && in_question >= B * (1 - conf_level) / 2 - 1e-8
}


# .check_boot_B(B, ci_method, conf_level)
#
# The `B` argument of the reliability family when a bootstrap interval is
# asked for: .dmar_check_B() with the kind the interval method implies
# ("percentile" reads quantiles, "bca" reads adjusted quantiles, and the
# two standard error methods take a standard deviation). Returns
# as.integer(B). A caller that checks conf_level after B still gives its
# own message for a bad level (see the first lines of the body).
.check_boot_B <- function(B, ci_method = "percentile", conf_level = 0.95) {
  kind <- .boot_method_kind(ci_method)
  # A level that is not a single number in (0, 1) is refused by the
  # caller's own check of conf_level, with its message; here B is then
  # held only to being a count.
  if (!is.numeric(conf_level) || length(conf_level) != 1L ||
      is.na(conf_level) || conf_level <= 0 || conf_level >= 1) {
    return(.dmar_check_count(B, "B", min = 2L))
  }
  .dmar_check_B(B, kind = kind, tail = (1 - conf_level) / 2,
                conf_level = conf_level)
}

# The kind of .dmar_boot_required() a reliability interval method is.
.boot_method_kind <- function(ci_method) {
  if (identical(ci_method, "bca")) "bca" else
    if (ci_method %in% c("bootstrap_se", "bootstrap_se_logistic")) "se" else
      "quantile"
}


# .dmar_boot_settle(B, B_used, kind, tail, conf_level, fun, quantity,
#                   dropped_because, floor, signal, what, position,
#                   replications)
#
# The verdict on a finished run. Returns
#   list(formed, status, B, B_used, n_dropped, B_required, caution,
#        share_dropped, context)
# with status one of "complete", "dropped", "exploratory",
# "insufficient_count", "insufficient_share". In order:
#   1. more than half dropped: not formed, "insufficient_share";
#   2. fewer than the required number REQUESTED (and the exploratory
#      option not set; .dmar_boot_count_met()): not formed,
#      "insufficient_count". A caller that checked its request with
#      .dmar_check_B() never reaches this;
#   3. fewer than 2 usable: not formed, "insufficient_count";
#   4. otherwise formed: "exploratory" when the request is below the
#      requirement, else "dropped" when any was dropped, else
#      "complete". `caution` is TRUE when the dropped replications number
#      at least the tail area beyond a limit (for kinds "se" and "p", the
#      tail (1 - conf_level)/2, as .boot_caution() has it).
# With signal = TRUE the conditions are given through .dmar_boot_signal().
# A caller that still has to apply the BCa position rule, or that composes
# its own message, passes signal = FALSE.
#   fun             : the function's name, without parentheses.
#   quantity        : its noun phrase ("the indirect effect", "kappa").
#   dropped_because : its clause for a dropped replication ("left kappa
#                     undefined (chance agreement of 1 in the resample)").
#   replications    : "bootstrap replications" or "Monte Carlo draws".
#   recorded        : FALSE where the result cannot carry the counts (a
#                     plot); the messages then do not say that it does.
#   defined_as      : the clause for the resamples kept when it is not
#                     "<quantity> is defined" (mlmr(): quantity "some
#                     coefficient", kept when "every coefficient is
#                     defined").
.dmar_boot_settle <- function(B, B_used, kind = "quantile", tail = NULL,
                              conf_level = NULL, fun, quantity,
                              dropped_because, floor = 100L, signal = TRUE,
                              what = "interval", position = NULL,
                              replications = "bootstrap replications",
                              recorded = TRUE, defined_as = NULL) {
  B <- as.integer(B)
  B_used <- as.integer(B_used)
  n_dropped <- B - B_used
  required <- .dmar_boot_required(kind, tail, floor)
  if (kind %in% c("se", "p") && is.null(tail))
    tail <- if (is.null(conf_level)) 0.025 else (1 - conf_level) / 2
  caution <- FALSE
  if (n_dropped > B / 2) {
    formed <- FALSE
    status <- "insufficient_share"
  } else if (!.dmar_boot_count_met(B, required)) {
    formed <- FALSE
    status <- "insufficient_count"
  } else if (B_used < 2L) {
    formed <- FALSE
    status <- "insufficient_count"
  } else {
    formed <- TRUE
    status <- if (is.na(required) || B < required) "exploratory" else
      if (n_dropped > 0L) "dropped" else "complete"
    caution <- n_dropped > 0L && n_dropped >= tail * B - 1e-8
  }
  verdict <- list(
    formed = formed, status = status, B = B, B_used = B_used,
    n_dropped = n_dropped, B_required = required, caution = caution,
    share_dropped = n_dropped / B,
    context = list(kind = kind, tail = tail, conf_level = conf_level,
                   fun = fun, quantity = quantity,
                   dropped_because = dropped_because, floor = floor,
                   what = what, position = position,
                   replications = replications, recorded = recorded,
                   defined_as = defined_as))
  if (signal) .dmar_boot_signal(verdict)
  verdict
}


# .dmar_boot_bca_positions(verdict, positions)
#
# The count of 10 applied where a BCa interval finally reads its limits.
# `positions` are the two adjusted quantile positions (proportions).
# required_at_positions = ceiling(10 / min(p_1, 1 - p_2)); like the count
# of the request it is held to the number requested, and when that is
# fewer the interval is withheld (status
# "insufficient_count", with `positions`, `B_suggest`, and B_required set
# to the count at those positions). Under the exploratory option the
# interval is not withheld. A verdict that is not formed is returned as it
# came. Signals nothing; the caller then calls .dmar_boot_signal().
.dmar_boot_bca_positions <- function(verdict, positions) {
  if (!isTRUE(verdict$formed) || .dmar_boot_exploratory()) return(verdict)
  smaller <- min(positions[1L], 1 - positions[2L])
  beyond <- .dmar_boot_beyond()
  need <- if (beyond <= 0) 0 else
    if (smaller <= 0) Inf else ceiling(beyond / smaller - 1e-8)
  if (verdict$B >= need) return(verdict)
  verdict$formed <- FALSE
  verdict$status <- "insufficient_count"
  verdict$caution <- FALSE
  verdict$positions <- as.numeric(positions)
  verdict$bca_withheld <- TRUE
  if (is.finite(need) && need <= .Machine$integer.max) {
    verdict$B_required <- as.integer(need)
    verdict$B_suggest <- as.integer(100 * ceiling(1.5 * need / 100))
  } else {
    verdict$B_required <- NA_integer_
    verdict$B_suggest <- NA_integer_
  }
  verdict
}


# .dmar_boot_request_note(B, B_used, required, replications)
#
# The sentence a warning about dropped replications carries on the count
# rule: the number required is held to the number requested, and a
# request can ask for more where replications can fail. When the
# replications that remained are fewer than the number required, a second
# sentence says so and says what to do. Used by .dmar_boot_message() and
# by the functions that compose their own account of a run. "" when the
# requirement is not a number.
.dmar_boot_request_note <- function(B, B_used, required,
                                    replications = "replications") {
  if (is.null(required) || is.na(required)) return("")
  need <- format(required, scientific = FALSE, trim = TRUE)
  paste0(
    " The ", need, " required are counted on the number requested (", B,
    "), not on the number that remained",
    if (B_used < required) {
      paste0("; the ", B_used, " that remained are fewer, so ",
             "request more ", replications, " where some can fail, enough that ",
             "at least ", need, " remain.")
    } else {
      paste0("; request more than ", need, " where ", replications, " can fail.")
    })
}


# The message of a verdict, by status. NULL for "complete".
.dmar_boot_message <- function(verdict, which = verdict$status) {
  ctx <- verdict$context
  B <- verdict$B
  B_used <- verdict$B_used
  n_dropped <- verdict$n_dropped
  required <- verdict$B_required
  tail <- ctx$tail
  quantile_kind <- ctx$kind %in% c("quantile", "bca")
  noun <- .dmar_boot_noun(ctx$kind, ctx$what, ctx$conf_level)
  share <- .pct(n_dropped / B)
  head <- paste0(ctx$fun, "(): ")
  whole <- function(v) format(v, scientific = FALSE, trim = TRUE)
  one_sided <- ctx$what %in% c("bound", "reference")
  result <- if (ctx$kind == "se") "standard error" else
    if (ctx$kind == "p") "p-value" else
      if (ctx$what == "bound") "bound" else
        if (ctx$what == "reference") "reference quantile" else "interval"
  recorded <- !identical(ctx$recorded, FALSE)
  # Monte Carlo draws are draws of the estimates from their estimated
  # sampling distribution, not resamples of the data, and the messages
  # say so: what was lost is a draw, and sparse data are not the reason.
  monte_carlo <- identical(ctx$replications, "Monte Carlo draws")
  unit <- if (monte_carlo) "draws" else "resamples"
  reps <- if (monte_carlo) "draws" else "replications"
  distribution <- if (monte_carlo) "distribution of the draws" else
    "bootstrap distribution"

  if (isTRUE(verdict$bca_withheld)) {
    pos <- verdict$positions
    if (is.na(required)) {
      return(paste0(
        head, "no BCa interval is reported. The bias correction and the ",
        "acceleration place a limit beyond every replication, so no ",
        "number of replications gives it; use the percentile interval."))
    }
    return(paste0(
      head, "no BCa interval is reported. The bias correction and the ",
      "acceleration place the limits at the ", .pct(pos[1L]), " and ",
      .pct(pos[2L]), " points of the bootstrap distribution, not at ",
      .pct(tail), " and ", .pct(1 - tail), ", and an interval is ",
      "reported only when at least 10 replications lie beyond each limit, ",
      "which takes a B of ", whole(required), " at those points (", B,
      " were requested). The positions are estimated from ",
      "the same replications and change a little from run to run, so ",
      "raise B to about ", whole(verdict$B_suggest), " (one and a half ",
      "times that count), or use the percentile interval."))
  }

  if (which == "insufficient_share") {
    return(paste0(
      head, "no ", result, " is reported. ", n_dropped, " of ", B, " ",
      ctx$replications, " (", share, ") ", ctx$dropped_because,
      " and were dropped. More than half is too many for ",
      if (result == "interval") "an interval" else paste("a", result),
      " at any B: the ", B_used, " that remained are a minority of the ",
      unit, ", selected by whether ", ctx$quantity, " could be ",
      "computed, and their spread describes those ", unit, " and not the ",
      "sampling distribution of the estimate. The point estimate is not ",
      "affected. A larger B does not help; ",
      if (monte_carlo)
        paste0("the draws are of the estimates, from their estimated ",
               "sampling distribution, and not resamples of the data, so ",
               "the quantity as it is defined is not a finite number over ",
               "most of that distribution (a ratio, a logarithm, or a ",
               "square root of estimates can be undefined).")
      else
        paste0("the data are too sparse for this bootstrap (a small ",
               "sample, a category or pattern that many resamples lose).")))
  }

  if (which == "insufficient_count") {
    # After the count rule moved to the request, two things give this
    # status besides a BCa interval withheld at its positions (above): a
    # request below the requirement that no check refused (an internal
    # caller), and fewer than two usable replications.
    dropped <- if (n_dropped > 0L) {
      paste0(n_dropped, " of ", B, " ", ctx$replications, " ",
             ctx$dropped_because, " and were dropped, and ", B_used,
             " remained. ")
    } else {
      ""
    }
    if (!is.na(required) && B < required) {
      return(paste0(
        head, "no ", result, " is reported. ", dropped, "B = ", B,
        " is below the ", whole(required), " that must be requested. ",
        .dmar_boot_reason(ctx$kind, tail, ctx$conf_level, ctx$floor,
                          ctx$what, ctx$position),
        " The point estimate is not affected."))
    }
    return(paste0(
      head, "no ", result, " is reported. ", dropped,
      if (result == "interval") "An interval" else paste("A", result),
      " needs at least two usable ", ctx$replications, ". The point ",
      "estimate is not affected.",
      if (n_dropped > 0L)
        paste0(" Consider why ", ctx$quantity, " is undefined in ", share,
               " of the ", unit, ".")))
  }

  if (which == "dropped") {
    return(paste0(
      head, n_dropped, " of ", B, " ", ctx$replications, " (", share, ") ",
      ctx$dropped_because, " and were dropped. The ", result,
      " is computed from the ", B_used, " that remained and describes ",
      unit, " in which ",
      if (is.null(ctx$defined_as)) paste(ctx$quantity, "is defined") else
        ctx$defined_as, ".",
      .dmar_boot_request_note(B, B_used, required, ctx$replications),
      if (recorded) " B, B_used, and n_dropped are recorded in the result.",
      if (isTRUE(verdict$caution) && quantile_kind)
        paste0(" Caution: the dropped ", reps, " are at least the ",
               .pct(tail), " of the ", distribution, " that lies ",
               "beyond ",
               if (one_sided) paste("the", result)
               else paste0("each limit of a ", noun),
               ", so where a limit falls depends on where those ",
               reps, " would have fallen. Report the ", result,
               " with these counts.")
      else if (isTRUE(verdict$caution))
        paste0(" Caution: the dropped ", reps, " are at least ",
               .pct(tail), " of those requested, enough to change the ",
               result, " by themselves. Report it with these counts.")))
  }

  if (which == "exploratory") {
    basis <- switch(ctx$kind,
      quantile = if (one_sided) paste("10 beyond the", result) else
        "10 beyond each limit",
      bca = "50 beyond each nominal limit",
      se = paste(whole(required), "for a standard error"),
      p = paste(whole(required), "for a Monte Carlo p-value"))
    return(paste0(
      head, "B = ", B, " is below the ", whole(required), " ", reps, " a",
      if (grepl("^8|^11|^18", noun)) "n " else " ", noun,
      " is reported from (", basis, "). The ",
      if (result == "interval") "limits are" else paste(result, "is"),
      " shown because options(dmar.bootstrap_exploratory = TRUE) is set.",
      if (quantile_kind)
        paste0(" With ", B_used, " ", reps, " ",
               if (one_sided) paste("the", result) else "each limit",
               " rests on about ",
               format(round(B_used * tail, 1), trim = TRUE),
               " of them and will move from seed to seed;")
      else
        paste0(" With ", B_used, " ", reps, " it will move from seed to ",
               "seed;"),
      " use the result to explore or to demonstrate, and report ",
      if (result == "interval") "an interval" else paste("a", result),
      " only from B of at least ", whole(required), ".",
      if (recorded) " B, B_used, and B_required are recorded in the result."))
  }
  NULL
}


# .dmar_boot_signal(verdict)
#
# The conditions of a verdict, each through .dmar_warn() and each carrying
# B, B_used, n_dropped, B_required, share_dropped, caution, status:
#   not formed                 dmar_bootstrap_insufficient
#   formed, some dropped       dmar_bootstrap_replications_dropped
#   status "exploratory"       dmar_bootstrap_exploratory (after the
#                              dropped warning when both apply)
# A withheld BCa interval also carries `positions` and `B_suggest`.
.dmar_boot_signal <- function(verdict) {
  give <- function(message, code, ...) {
    .dmar_warn(message, code, B = verdict$B, B_used = verdict$B_used,
               n_dropped = verdict$n_dropped,
               B_required = verdict$B_required,
               share_dropped = verdict$share_dropped,
               caution = verdict$caution, status = verdict$status, ...)
  }
  if (!isTRUE(verdict$formed)) {
    if (isTRUE(verdict$bca_withheld)) {
      give(.dmar_boot_message(verdict), "bootstrap_insufficient",
           positions = verdict$positions, B_suggest = verdict$B_suggest)
    } else {
      give(.dmar_boot_message(verdict), "bootstrap_insufficient")
    }
    return(invisible(verdict))
  }
  if (verdict$n_dropped > 0L) {
    give(.dmar_boot_message(verdict, "dropped"),
         "bootstrap_replications_dropped")
  }
  if (identical(verdict$status, "exploratory")) {
    give(.dmar_boot_message(verdict, "exploratory"),
         "bootstrap_exploratory")
  }
  invisible(verdict)
}


# .dmar_boot_signal_text(verdict, text)
#
# For a function that composes its own account of a run (the factors, the
# improper replications kept, the reasons a refit failed): signals `text`
# with the class the verdict implies, dmar_bootstrap_insufficient when no
# interval was formed, dmar_bootstrap_replications_dropped when one was
# and replications were dropped, and dmar_bootstrap_note otherwise, each
# with the fields of .dmar_boot_signal(); then the exploratory warning of
# the rule when the status is "exploratory". `text` may be NULL. The text
# is given the head every message of the rule has, the name of the
# function the verdict was settled for ("average_variance_extracted(): "),
# so that a warning met in a long script says where it came from; a text
# that already begins with it is left as it is.
.dmar_boot_signal_text <- function(verdict, text = NULL) {
  head <- paste0(verdict$context$fun, "(): ")
  if (!is.null(text) && length(verdict$context$fun) == 1L &&
      !startsWith(text, head)) {
    text <- paste0(head, text)
  }
  give <- function(message, code) {
    .dmar_warn(message, code, B = verdict$B, B_used = verdict$B_used,
               n_dropped = verdict$n_dropped,
               B_required = verdict$B_required,
               share_dropped = verdict$share_dropped,
               caution = verdict$caution, status = verdict$status)
  }
  if (!is.null(text)) {
    if (!isTRUE(verdict$formed)) {
      give(text, "bootstrap_insufficient")
    } else if (verdict$n_dropped > 0L) {
      give(text, "bootstrap_replications_dropped")
    } else {
      give(text, "bootstrap_note")
    }
  }
  if (isTRUE(verdict$formed) && identical(verdict$status, "exploratory")) {
    give(.dmar_boot_message(verdict, "exploratory"), "bootstrap_exploratory")
  }
  invisible(verdict)
}


# .dmar_boot_attrs(out, verdict)
#
# Writes the counts of a run on its result: attributes B, B_used,
# n_dropped, B_required (integers) and bootstrap_status (character).
# A result that ran no replications carries none of them.
.dmar_boot_attrs <- function(out, verdict) {
  attr(out, "B") <- as.integer(verdict$B)
  attr(out, "B_used") <- as.integer(verdict$B_used)
  attr(out, "n_dropped") <- as.integer(verdict$n_dropped)
  attr(out, "B_required") <- as.integer(verdict$B_required)
  attr(out, "bootstrap_status") <- verdict$status
  out
}


# .dmar_boot_mark_planner(out, B, kind, tail, conf_level, fun, quantity,
#                         replications)
#
# For a planner that runs B replications inside each of many simulated
# studies (ss_aipe_indirect_effect(), ss_aipe_indirect_effect_sensitivity(),
# ss_aipe_reliability_sensitivity()). The request was checked by
# .dmar_check_B(); here the run is marked once for the call: the
# exploratory warning when B is below the requirement under
# options(dmar.bootstrap_exploratory = TRUE), and the five attributes of
# .dmar_boot_attrs() with B_used = B and n_dropped = 0, since what is
# recorded is the request each simulated study ran and replications dropped
# inside a study are not totaled. The attribute bootstrap_scope =
# "per_study" tells the print method to show the line only for an
# exploratory run (.dmar_boot_footer_per_study()).
.dmar_boot_mark_planner <- function(out, B, kind = "quantile", tail,
                                    conf_level, fun, quantity,
                                    replications = "bootstrap replications",
                                    signal = TRUE) {
  verdict <- .dmar_boot_settle(
    B, B, kind = kind, tail = tail, conf_level = conf_level, fun = fun,
    quantity = quantity, dropped_because = "returned no value",
    signal = FALSE, replications = replications)
  if (signal) .dmar_boot_signal_text(verdict)
  out <- .dmar_boot_attrs(out, verdict)
  attr(out, "bootstrap_scope") <- "per_study"
  out
}

# The printed line of a planner marked by .dmar_boot_mark_planner(): NULL
# unless the run was exploratory.
.dmar_boot_footer_per_study <- function(x) {
  if (!identical(attr(x, "bootstrap_status", exact = TRUE), "exploratory"))
    return(NULL)
  paste0("Each simulated study ran B = ", attr(x, "B", exact = TRUE),
         " replications, below the ",
         format(attr(x, "B_required", exact = TRUE), scientific = FALSE,
                trim = TRUE),
         " required. Exploratory: below the number required for a ",
         "reported interval.")
}


# .dmar_boot_footer(x, label)
#
# The line a print method shows for a result that carries the counts:
#   Bootstrap: 2000 replications, 1988 used (12 dropped; at least 400
#   must be requested)
# with "No interval: ..." in place of the parenthesis when the interval
# was withheld, and a closing sentence in exploratory mode. NULL for a
# result without the attributes. `x` is a result or a list with the five
# elements B, B_used, n_dropped, B_required, bootstrap_status.
.dmar_boot_footer <- function(x, label = "Bootstrap") {
  get <- function(name) {
    if (is.list(x) && !is.data.frame(x) && !is.null(x[[name]])) x[[name]]
    else attr(x, name, exact = TRUE)
  }
  status <- get("bootstrap_status")
  if (is.null(status) && is.list(x) && !is.data.frame(x)) status <- x$status
  B <- get("B")
  B_used <- get("B_used")
  if (is.null(status) || is.null(B) || is.null(B_used)) return(NULL)
  n_dropped <- get("n_dropped")
  if (is.null(n_dropped)) n_dropped <- B - B_used
  required <- get("B_required")
  need <- if (is.null(required) || is.na(required)) NULL else
    format(required, scientific = FALSE, trim = TRUE)
  line <- paste0(label, ": ", B, " replications, ", B_used, " used")
  line <- paste0(line, switch(status,
    insufficient_share = paste0(". No interval: more than half of the ",
                                "replications were dropped."),
    insufficient_count = if (B_used < 2L) {
      ". No interval: fewer than two replications were usable."
    } else {
      paste0(". No interval: fewer than the ",
             if (is.null(need)) "number" else need,
             " required were requested.")
    },
    paste0(" (", n_dropped, " dropped",
           if (!is.null(need))
             paste0("; at least ", need, " must be requested"),
           ")")))
  if (identical(status, "exploratory")) {
    line <- paste0(line, ". Exploratory: below the number required for a ",
                   "reported interval.")
  }
  line
}


#' How DMAR Counts and Requires Bootstrap Replications
#'
#' One rule governs every DMAR function that reads an interval, a bound, a
#' reference quantile, a standard error, or a Monte Carlo \emph{p}-value
#' from \code{B} bootstrap or Monte Carlo replications: how many
#' replications the request must ask for, what happens to a replication
#' that returns no value, what the result records, and when no interval is
#' reported. This page states the rule and its reasons; each function's
#' page says how the rule applies to it.
#'
#' @section How many replications an interval needs:
#' A percentile limit is a quantile of the \code{B} replications, so it is
#' located by the replications that lie beyond it. With \code{B} = 100 and
#' a 95 percent interval there are two or three of them on each side, and
#' the limit moves with every new seed. DMAR runs an interval read from
#' quantiles only when the number of replications requested puts at least
#' 10 beyond each limit, and never from a request of fewer than 100: 200,
#' 400, and 2,000 replications at the 90, 95, and 99 percent levels. The
#' count is held to the number requested, \code{B}; what happens when some
#' of them return no value is decided afterward and is described under
#' \emph{Replications that return no value}. In symbols, with \eqn{t} the
#' smallest tail area a reported limit cuts off (\eqn{(1 - C)/2} for a
#' two-sided interval at level \eqn{C}),
#' \deqn{B_{\mathrm{required}} = \max\{100, \lceil 10 / t \rceil\}.}
#' With 10, the share of the bootstrap distribution a limit actually cuts
#' off has a relative standard deviation of about 30 percent whatever the
#' shape of that distribution (a nominal 2.5 percent is between about 1.2
#' and 4.5 percent in 95 of 100 runs, depending slightly on how the
#' quantile is interpolated). How far the limit itself moves depends on
#' the shape: with 400 replications at the 95 percent level its standard
#' deviation from run to run is about 0.13 bootstrap standard errors for a
#' bell-shaped distribution, about 0.20 for one with tails as heavy as a
#' \emph{t} distribution on 5 degrees of freedom, and about 0.29 on the
#' long side of a skewed one. Ten is a floor and not a recommendation. The
#' defaults (1,000 where each replication refits a model, 2,000 to 10,000
#' elsewhere) put 25 to 250 replications beyond each limit at the 95
#' percent level, and at 1,000 the standard deviation of a limit of a
#' bell-shaped distribution is about 0.08 standard errors. A BCa interval
#' needs more, 1,000, 2,000, and 10,000 at the three levels, for the
#' reason given under \emph{BCa intervals}. An interval built from the
#' standard deviation of the replications (\code{"bootstrap_se"}) is
#' reported from 100, at which that standard deviation has a coefficient
#' of variation of about 7 percent for a bell-shaped distribution (12 to 14
#' percent for the heavy-tailed and the skewed shapes above), and a Monte
#' Carlo \emph{p}-value from
#' 100, at which its standard error is at most 0.05. A request below the
#' number required is refused before anything is computed, with a
#' condition of class \code{dmar_count_invalid} that states the number.
#' None of these numbers concerns the coverage of the bootstrap for a
#' given statistic, which no number of replications repairs; they limit
#' only the error added by drawing a finite number of resamples. The
#' numbers in this paragraph were computed for DMAR by simulation (20,000
#' runs each); the statement about the share cut off follows from the Beta
#' distribution of an order statistic (the tail area beyond the
#' \eqn{k}-th of \eqn{B} ordered draws from a continuous distribution is
#' distributed as Beta(\eqn{k}, \eqn{B - k + 1})). No published source is
#' cited for the counts.
#'
#' A limit at quantile position \eqn{p} has the tail area
#' \eqn{\min(p, 1 - p)}. For an interval with unequal tail areas
#' (\code{alpha_lower} and \code{alpha_upper} of the ratio functions and
#' of \code{\link{ci_eta_squared_generalized}}) the smaller area decides;
#' for a one-sided bound (\code{\link{htmt}}) and for a reference quantile
#' (\code{\link{parallel_analysis}}, \code{\link{factor_count}}) it is the
#' area beyond that one quantile, so a bound at the 95th percentile needs
#' 200 replications and one at the median needs the floor of 100.
#'
#' The functions that compute their limits directly read them with
#' \code{\link[stats]{quantile}} (its default, type 7); the reliability
#' family and \code{\link{ci_bifactor_indices}} read them with
#' \code{boot::boot.ci()}, which places a limit at the
#' \eqn{(B + 1)\alpha} position among the ordered replications. The two
#' differ by a fraction of the distance between adjacent replications.
#'
#' @section Replications that return no value:
#' A resample can leave a statistic undefined (a constant variable, a
#' model that does not converge, a denominator of zero). Such a
#' replication is dropped and counted, never replaced. Every result
#' records \code{B} (requested), \code{B_used} (the number the interval
#' rests on), \code{n_dropped}, and \code{B_required} (the fewest
#' replications a request may ask for) as attributes, which
#' \code{\link{dmar_metadata}} returns, with \code{bootstrap_status}, and
#' one warning of class \code{dmar_bootstrap_replications_dropped} states
#' the counts. An interval from the replications that remain describes the
#' resamples in which the statistic is defined. When the dropped
#' replications number at least the share of the distribution beyond one
#' limit (2.5 percent of \code{B} for a 95 percent interval), they could
#' by themselves decide where that limit falls, and the warning says so
#' (\code{caution} is \code{TRUE} in the condition). When more than half
#' were dropped no interval is reported: the limits are \code{NA}, the
#' estimate is returned, and the warning has class
#' \code{dmar_bootstrap_insufficient}. The call does not
#' stop for this, so an estimate that took minutes to compute is kept.
#' With more than half dropped a larger \code{B} does not help: the
#' replications that remain are a minority of the resamples, selected by
#' whether the statistic could be computed, and their spread describes
#' those resamples and not the sampling distribution of the estimate.
#'
#' The number required is counted on the request and not on the
#' replications that remain. A request of exactly 400 at the 95 percent
#' level from which 3 replications are dropped reports its interval from
#' the 397 that remain, with the warning. The two rules about shares, in
#' the paragraph
#' above, then decide: more than half dropped withholds the interval, and
#' dropped replications that reach the share beyond one limit bring the
#' caution. The warning states how many were dropped, and when those that
#' remain are fewer than the number required it says so and says to
#' request more. Where replications can fail (a model refit in every
#' resample, sparse categories, a small sample), ask for more than the
#' minimum, so that the limits rest on at least the number the rule
#' intends; the defaults do.
#'
#' \code{bootstrap_status} is one of \code{"complete"} (nothing dropped),
#' \code{"dropped"} (some dropped, the interval formed),
#' \code{"exploratory"} (formed from a request below the required number,
#' under the option described below), \code{"insufficient_count"} (no
#' interval: a BCa interval whose adjusted positions take more
#' replications than were requested, or fewer than two usable
#' replications), and
#' \code{"insufficient_share"} (more than half dropped, no interval).
#' With several statistics in one call (the pairs of \code{\link{htmt}},
#' the effects of \code{\link{mediation_mbco}}, the coefficients of
#' \code{\link{mlmr}}), \code{B_used} is the smallest count any reported
#' limit rests on and \code{n_dropped} is \code{B - B_used}. Where a table
#' already has rows for the counts it keeps them (\code{\link{mediate}}
#' has rows \code{B} and \code{B_used}, \code{\link{krippendorff_alpha}}
#' the row \code{B_used}); the reliability family keeps its fuller
#' \code{bootstrap} attribute, which gains \code{required} and
#' \code{bootstrap_status}.
#'
#' @section BCa intervals:
#' The BCa interval reads its limits from quantiles that the bias
#' correction and the acceleration move away from the nominal ones. For a
#' 95 percent interval the positions might be 0.9 and 99.1 percent in
#' place of 2.5 and 97.5, and the limit at 0.9 percent has few
#' replications beyond it. The count of 10 is applied where each limit
#' finally sits, again to the number of replications requested, so a BCa
#' interval can be withheld after the replications are run (warning of
#' class \code{dmar_bootstrap_insufficient}, with the positions and the
#' count that would have served). To make that uncommon,
#' a BCa interval is run only with 50 replications beyond each nominal
#' limit: 1,000, 2,000, and 10,000 at the 90, 95, and 99 percent levels.
#' In simulations of an indirect effect and of a correlation at the 95
#' percent level, the adjusted positions asked for more than 2,000
#' replications in at most 14 percent of data sets (an indirect effect
#' with 50 cases and small paths) and in 0 to 3 percent of the others, and
#' for more than 10,000 in about 1 percent of the data sets of that hardest
#' condition and in none of the others (400 data sets per condition).
#'
#' When the bias correction or the acceleration is itself undefined (a
#' leave-one-out estimate that cannot be computed, every replicate on one
#' side of the estimate, influence values without spread), the functions
#' differ by family. \code{\link{mediate}},
#' \code{\link{cohen_kappa}}, and \code{\link{fleiss_kappa}} refuse the
#' BCa interval with an error of class \code{dmar_bca_undefined} whose
#' message names the cause; \code{\link{ci_ratio_of_means}} and
#' \code{\link{ci_ratio_of_means_paired}} stop with an error that names
#' \code{ci_method = "percentile"}; the reliability family and
#' \code{\link{ci_bifactor_indices}} report the percentile interval in its
#' place, with a warning, and record the substitution in the
#' \code{ci_method} attribute.
#'
#' @section Fewer replications, for exploration:
#' A first look at an analysis, a help page example, and a simulation that
#' averages over many data sets do not need a reportable interval from
#' each run. \code{options(dmar.bootstrap_exploratory = TRUE)} lets a
#' function run with any \code{B} of 2 or more. Nothing else changes: the
#' result still records \code{B_required}, its footer says that the
#' interval is below the requirement, and every call gives a warning of
#' class \code{dmar_bootstrap_exploratory}. A planner that runs \code{B}
#' replications inside each of many simulated studies
#' (\code{\link{ss_aipe_indirect_effect}},
#' \code{\link{ss_aipe_indirect_effect_sensitivity}},
#' \code{\link{ss_aipe_reliability_sensitivity}}) gives the warning once
#' for the call and records the request (\code{B_used} is then \code{B},
#' the number each simulated study ran), and
#' \code{\link{plot_mediation_mbco}}, whose result is a plot, gives the
#' warning and records nothing. The share rules still apply.
#' An interval obtained this way is not one to report. Set the option for
#' the calls that need it and restore it afterward:
#' \preformatted{old <- options(dmar.bootstrap_exploratory = TRUE)
#' mediate(d, x = "x", m = "m", y = "y", B = 50, seed = 113)
#' options(old)}
#'
#' @section Monte Carlo intervals and reference quantiles:
#' The same rule holds for limits read from the quantiles of Monte Carlo
#' draws (\code{ci_method = "monte_carlo"} of \code{\link{mediate}} and of
#' \code{\link{mediation_mbco}}, the band of
#' \code{\link{plot_mediation_mbco}}) and for the reference quantile of
#' \code{\link{parallel_analysis}} and \code{\link{factor_count}}, with
#' "draw" or "reference data set" in place of "replication". The
#' Bollen-Stine \emph{p}-value of \code{\link{sem_exact_fit_bootstrap}}
#' is a proportion and not a quantile: it needs a request of 100 samples,
#' is withheld when more than half were set aside, and is
#' reported with the Wilson interval of its own Monte Carlo error.
#'
#' @section What is not a bootstrap:
#' The rule does not govern \code{n_resamples} of
#' \code{\link{randomization_test}},
#' \code{\link{randomization_test_paired}}, and
#' \code{\link{procrustes_phi}}, or the number of simulated studies
#' \code{G} of the sample size planning and sensitivity
#' functions. Those counts are checked as counts and their functions
#' report their own Monte Carlo error.
#'
#' @section Reproducibility:
#' Every function that draws replications has a \code{seed} argument
#' (\code{boot_seed} in \code{\link{mlmr}} and \code{\link{mlmr_mv}}) that
#' defaults to \code{NULL}: the replications are then drawn from, and
#' advance, the current random number stream. A supplied seed is set for
#' the call only and the caller's generator state is restored when the
#' function returns, so a reproducible call does not disturb the random
#' numbers of the script around it. Two runs with the same seed and the
#' same \code{B} give the same limits; the same seed with a different
#' \code{B} does not.
#'
#' @section Conditions:
#' \describe{
#' \item{\code{dmar_count_invalid}}{Error (a \code{dmar_refusal}). \code{B}
#'   is not a single whole number of at least the number required. Fields
#'   \code{arg}, \code{min} (the number required), \code{max}, and
#'   \code{reason} (the sentence that explains the number).}
#' \item{\code{dmar_bootstrap_replications_dropped}}{Warning. Replications
#'   were dropped and the interval was formed from the rest. The message
#'   states how many, and says to request more when those that remain are
#'   fewer than \code{B_required}. Fields
#'   \code{B}, \code{B_used}, \code{n_dropped}, \code{B_required},
#'   \code{share_dropped}, \code{caution}, and \code{status}.}
#' \item{\code{dmar_bootstrap_insufficient}}{Warning. No interval is
#'   reported: more than half of the replications were dropped
#'   (\code{status} \code{"insufficient_share"}), or the adjusted
#'   positions of a BCa interval take more replications than were
#'   requested, or fewer than two replications were usable
#'   (\code{"insufficient_count"}). The same
#'   fields; for a BCa interval withheld at its adjusted positions also
#'   \code{positions} (the two positions, as proportions) and
#'   \code{B_suggest} (one and a half times the count those positions
#'   take, rounded up to a hundred).}
#' \item{\code{dmar_bootstrap_exploratory}}{Warning. The result was
#'   computed from fewer replications than a reported result needs,
#'   because \code{options(dmar.bootstrap_exploratory = TRUE)} is set. The
#'   same fields as the dropped warning.}
#' \item{\code{dmar_bootstrap_note}}{Warning. Something else about the
#'   run that the reader of the interval should know, with every
#'   replication used and the interval formed: in the reliability family
#'   and \code{\link{ci_bifactor_indices}}, a BCa interval replaced by the
#'   percentile interval (the \code{ci_method} attribute records the
#'   substitution), or replicates that all have the same value, so that
#'   the interval has no width; in
#'   \code{\link{average_variance_extracted}} and
#'   \code{\link{parallel_analysis}}, the function's account of a run in
#'   which nothing was dropped (replications from an improper solution
#'   that were kept, reference data sets whose fits gave warnings). The
#'   same fields as the dropped warning where the function records them.}
#' \item{\code{dmar_bca_undefined}}{Error (a \code{dmar_refusal}). The
#'   bias correction or the acceleration of a BCa interval is undefined
#'   (\code{\link{mediate}}, \code{\link{cohen_kappa}},
#'   \code{\link{fleiss_kappa}}).}
#' }
#' The warnings of the reliability family and of
#' \code{\link{ci_bifactor_indices}}, which share one bootstrap engine,
#' begin with the name of the function that was called
#' (\code{"reliability_omega(): "}) and carry that name in the field
#' \code{fun}; \code{\link{reliability}} passes its call on, so its
#' warning names the function that ran the bootstrap.
#' A program can muffle or read any of the warnings by class:
#' \preformatted{withCallingHandlers(
#'   mediate(d, x = "x", m = "m", y = "y", seed = 113),
#'   dmar_bootstrap_replications_dropped = function(w) {
#'     message(w$n_dropped, " of ", w$B, " dropped")
#'     invokeRestart("muffleWarning")
#'   })}
#'
#' @section The functions the rule governs:
#' "Required" is the fewest replications a request may ask for at the 95 percent level
#' (for the reference quantiles, at the 95th percentile). Kind Q reads
#' quantiles, S takes a standard deviation of the replications, and P is a
#' Monte Carlo \emph{p}-value.
#' \tabular{llrrl}{
#' \strong{Function} \tab \strong{Kind} \tab \strong{Default \code{B}} \tab \strong{Required} \tab \strong{Why that default} \cr
#' \code{\link{mediate}}, percentile and Monte Carlo \tab Q \tab 2000 \tab 400 \tab each replication refits two regressions \cr
#' \code{\link{mediate}}, BCa \tab Q, BCa \tab 2000 \tab 2000 \tab the same \cr
#' \code{\link{mediation_mbco}}, \code{\link{plot_mediation_mbco}} (Monte Carlo) \tab Q \tab 10000 \tab 400 \tab a draw evaluates the effects, no refit \cr
#' \code{\link{reliability_alpha}}, \code{\link{reliability_kr20}}, \code{\link{reliability_omega}}, \code{\link{reliability}}, \code{\link{reliability_H}}: percentile \tab Q \tab 10000 \tab 400 \tab a closed form or a small model per replication \cr
#' the same functions: BCa \tab Q, BCa \tab 10000 \tab 2000 \tab the same \cr
#' the same functions: \code{"bootstrap_se"}, \code{"bootstrap_se_logistic"} \tab S \tab 10000 \tab 100 \tab the same \cr
#' \code{\link{reliability_omega_categorical}} \tab Q, BCa, or S \tab 10000 (1000 under \code{missing = "mi"}) \tab 400, 2000, or 100 \tab under multiple imputation each replication is imputed several times \cr
#' \code{\link{ci_bifactor_indices}} \tab Q or BCa \tab 1000 \tab 400 or 2000 \tab each replication refits the bifactor model \cr
#' \code{\link{average_variance_extracted}} \tab Q \tab 1000 \tab 400 \tab each replication refits the model \cr
#' \code{\link{fraction_missing_information}} \tab Q \tab 1000 \tab 400 \tab each replication refits the model three times \cr
#' \code{\link{cohen_kappa}}, \code{\link{fleiss_kappa}} \tab Q or BCa \tab 10000 \tab 400 or 2000 \tab kappa is cheap to recompute \cr
#' \code{\link{krippendorff_alpha}} \tab Q \tab 1000 \tab 400 \tab each resample recomputes all coincidences \cr
#' \code{\link{htmt}} \tab Q, one-sided \tab 0 (no bootstrap) \tab 200 \tab the bound is opt-in \cr
#' \code{\link{ci_ratio_of_means}}, \code{\link{ci_ratio_of_means_paired}} \tab Q or BCa \tab 10000 \tab 400 or 2000 \tab a ratio of means is cheap to recompute \cr
#' \code{\link{R2_mixed_effects}} \tab Q \tab 1000 \tab 400 \tab each replication refits the mixed model \cr
#' \code{\link{ci_eta_squared_generalized}} \tab Q \tab 10000 \tab 1000 \tab its own floor of 1000 \cr
#' \code{\link{mlmr}}, \code{\link{mlmr_mv}} \tab Q \tab 1000 \tab 400 \tab each resample refits the model \cr
#' \code{\link{parallel_analysis}}, \code{\link{factor_count}} \tab Q, reference quantile \tab 1000 \tab 200 \tab each reference data set is generated and fit \cr
#' \code{\link{sem_exact_fit_bootstrap}} \tab P \tab 1000 \tab 100 \tab each sample refits the model
#' }
#' The planning functions \code{\link{ss_aipe_indirect_effect}},
#' \code{\link{ss_aipe_indirect_effect_sensitivity}}, and
#' \code{\link{ss_aipe_reliability_sensitivity}} run an interval inside
#' each simulated study; their pages state the inner \code{B} they use.
#'
#' @references
#' Efron, B., & Tibshirani, R. J. (1993). \emph{An introduction to the
#'   bootstrap}. New York, NY: Chapman & Hall/CRC. (For the percentile and
#'   BCa intervals themselves, as the function pages cite it. The counts
#'   on this page are not taken from it.)
#'
#' @seealso \code{\link{dmar_metadata}} for reading the counts from a
#'   result, and the pages of the functions in the table.
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @keywords utilities
#'
#' @name dmar_bootstrap
#' @aliases dmar.bootstrap_exploratory
NULL
