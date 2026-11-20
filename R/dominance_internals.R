# Internal helpers shared by the order-based (dominance) effect sizes:
# vargha_delaney_A(), cliff_delta(), and probability_of_superiority_paired().
# Each of the three depends on its observations only through their order, so
# the three read their input through one contract (.dominance_scores()),
# compute the pairwise comparisons from ranks (.dominance_placements()), and
# qualify a result whose estimated variance is 0 with one classed warning
# (.dominance_degenerate_interval()).

# The check of 'conf_level' the three functions share.
.dominance_check_conf_level <- function(conf_level) {
  if (!is.numeric(conf_level) || length(conf_level) != 1L ||
      is.na(conf_level) || conf_level <= 0 || conf_level >= 1) {
    .dmar_refuse_argument(
      "'conf_level' must be a single number in (0, 1).", arg = "conf_level")
  }
  invisible(conf_level)
}

# .dominance_vector(x, name)
#
# One set of observations as a vector, before any coercion. A one-column
# matrix (the result of scale()) and a one-dimensional array (the result of
# tapply()) are vectors with a dim attribute, which is dropped. A list (a
# data frame included), a matrix of two or more columns, and an array of
# three or more dimensions hold more than one set of observations and are
# refused: as.numeric() would fail on the first and string the columns of
# the others end to end.
.dominance_vector <- function(x, name) {
  if (inherits(x, "POSIXlt")) x <- as.POSIXct(x)
  shape <- NULL
  d <- dim(x)
  if (is.data.frame(x)) {
    shape <- "a data frame"
  } else if (is.list(x)) {
    shape <- "a list"
  } else if (length(d) == 2L && d[2L] != 1L) {
    shape <- paste0("a matrix with ", d[2L], " columns")
  } else if (length(d) > 2L) {
    shape <- paste0("an array with ", length(d), " dimensions")
  }
  if (!is.null(shape)) {
    .dmar_refuse(c("'", name, "' must be a numeric or logical vector or an ",
                   "ordered factor; it is ", shape, ". Nothing was ",
                   "computed."),
                 "scores_not_ordered", arg = name)
  }
  if (!is.null(d)) dim(x) <- NULL
  x
}

# The kind of scale a vector of observations is on, read from the object
# the user gave and not from what as.numeric() would make of it.
.dominance_kind <- function(x) {
  if (is.ordered(x)) "ordered"
  else if (is.factor(x)) "unordered"
  else if (inherits(x, "difftime")) "difftime"
  else if (inherits(x, "Date")) "date"
  else if (inherits(x, "POSIXct")) "time"
  else if (is.numeric(x) || is.logical(x)) "numeric"
  else if (is.character(x)) "character"
  else "other"
}

# What a vector of observations is, in words, for the message that says
# why an ordered factor cannot be compared with it.
.dominance_kind_phrase <- function(x) {
  switch(.dominance_kind(x),
         ordered   = "an ordered factor",
         unordered = "an unordered factor",
         difftime  = "a time difference (class difftime)",
         date      = "a vector of dates (class Date)",
         time      = "a vector of date-times (class POSIXct)",
         numeric   = if (is.logical(x)) "a logical vector" else
                       "a numeric vector",
         character = "a character vector",
         paste0("of class ", paste(class(x), collapse = "/")))
}

# Up to six elements of `x` joined by `sep`, then ", and so on" when there
# are more, so that a message shows the levels of a factor without
# printing a long list (and without dots that would run into the period
# of the sentence).
.dominance_show <- function(x, sep = ", ") {
  x <- as.character(x)
  if (length(x) > 6L) {
    return(paste0(paste(x[1:6], collapse = sep), ", and so on"))
  }
  paste(x, collapse = sep)
}

# .dominance_scores(x, y, x_name, y_name, fun, paired = FALSE)
#
# The input contract of the three functions. `x` and `y` are the two sets
# of observations as the user gave them, `x_name` and `y_name` the names to
# use for them in a message, `fun` the name of the calling function.
# Returns list(x, y, n_missing) with `x` and `y` numeric and complete.
#
# Accepted: numeric and logical vectors; two ordered factors with the same
# levels in the same order, through their integer codes; Date and POSIXct
# vectors, through as.numeric(). Refused, with a condition of class
# dmar_scores_not_ordered: an unordered factor (its codes follow the
# alphabetical order of its labels, so factor(c("10", "9", "8")) has codes
# 1, 3, 2); an ordered factor beside anything else, or beside an ordered
# factor with other levels (after droplevels() the same label has another
# code); two difftime objects in different units, or a difftime beside
# anything else (each object is stored in its own units, so 1, 2, 3 hours
# beside 30, 90, 150 minutes would be compared as 1, 2, 3 against 30, 90,
# 150; two in the same units, as the two groups of one column always are,
# are on one scale and are accepted); a Date beside a POSIXct (days
# against seconds); a character vector (strings have only alphabetical
# order); and complex, list, and matrix input. Each message says which
# case applies, what would have been computed, and what to supply.
#
# Missing values (NA and NaN) are removed group by group, or pair by pair
# when `paired` is TRUE, and counted in `n_missing`: c(group_1 = , group_2 = )
# for two groups, one integer for pairs. Infinite values are kept; they have
# a place in the order.
.dominance_scores <- function(x, y, x_name, y_name, fun, paired = FALSE) {
  x <- .dominance_vector(x, x_name)
  y <- .dominance_vector(y, y_name)
  kind_x <- .dominance_kind(x)
  kind_y <- .dominance_kind(y)
  kinds  <- c(kind_x, kind_y)
  # What the two sets of observations are called: two groups, or for the
  # paired function the two members of each pair.
  both  <- if (paired) "both members of the pairs" else "both groups"
  two   <- if (paired) "the two members" else "the two groups"
  other <- if (paired) "the other member" else "the other group"

  # Each refusal says what the function needs, why the object given does
  # not supply it, what would have been computed from it, and what to
  # supply instead.
  if ("unordered" %in% kinds) {
    i    <- match("unordered", kinds)
    name <- c(x_name, y_name)[i]
    lev  <- levels(list(x, y)[[i]])
    ranked <- if (length(lev) >= 2L) {
      paste0("the statistic would rank ", .dominance_show(lev, " below "),
             ", and it would change")
    } else {
      "the statistic would change"
    }
    .dmar_refuse(c(
      fun, "() cannot use an unordered factor ('", name, "'). The ",
      "statistic is computed from the order of the observations and from ",
      "nothing else, and the levels of an unordered factor are labels ",
      "with no order. R stores them as whole number codes in the order of ",
      "levels(), which is alphabetical unless it was set",
      # A factor without levels (every value missing) has none to show.
      if (length(lev)) paste0(
        ": here the levels are ", .dominance_show(lev), ", stored as ",
        .dominance_show(seq_along(lev))),
      ". Computed from those codes, ",
      ranked, " if a label were renamed or the levels reordered, although ",
      "the data would be the same. Nothing was computed. State the order: ",
      "supply ", both, " as ordered factors with the same levels from ",
      "lowest to highest, for example factor(", name, ", levels = c(...), ",
      "ordered = TRUE), or supply numeric scores (as.numeric(as.character(",
      name, ")) when the labels are numbers)."),
      "scores_not_ordered", arg = name)
  }
  if ("character" %in% kinds) {
    name <- c(x_name, y_name)[match("character", kinds)]
    .dmar_refuse(c(
      fun, "() cannot use a character vector ('", name, "'). The statistic ",
      "is computed from the order of the observations, and character ",
      "strings have only alphabetical order, which is rarely the order ",
      "meant (\"10\" sorts before \"9\", and \"high\" before \"low\"). ",
      "Nothing was computed. Supply numeric scores (as.numeric(", name,
      ") when the strings are numbers), or an ordered factor with the ",
      "levels from lowest to highest, for example factor(", name,
      ", levels = c(...), ordered = TRUE)."),
      "scores_not_ordered", arg = name)
  }
  if ("ordered" %in% kinds) {
    if (kind_x != kind_y) {
      i <- match("ordered", kinds)
      .dmar_refuse(c(
        fun, "() needs ", both, " on one scale. '", c(x_name, y_name)[i],
        "' is an ordered factor and '", c(x_name, y_name)[3L - i], "' is ",
        .dominance_kind_phrase(list(x, y)[[3L - i]]), ". The codes of an ",
        "ordered factor (1 for its lowest level, 2 for the next) are not ",
        "the numbers in ", other, ", so comparing them would compare ",
        "positions in a list of levels with scores. Nothing was computed. ",
        "Supply ", both, " as ordered factors with the same levels in ",
        "the same order, or both as numeric scores."),
        "scores_not_ordered")
    }
    if (!identical(levels(x), levels(y))) {
      lev_x <- levels(x)
      lev_y <- levels(y)
      # A label the two share under different codes shows the problem in
      # the user's own data; when the shared labels have the same codes,
      # or no label is shared, the sentence says what does hold.
      shared <- intersect(lev_x, lev_y)
      moved  <- shared[match(shared, lev_x) != match(shared, lev_y)]
      why <- if (length(moved)) {
        paste0("so the same label can carry a different code in ", two,
               " (", moved[1L], " is ", match(moved[1L], lev_x),
               " in one and ", match(moved[1L], lev_y), " in the other)")
      } else if (!length(shared) && length(lev_x) && length(lev_y)) {
        paste0("so the same code stands for different labels in ", two,
               " (1 is ", lev_x[1L], " in one and ", lev_y[1L],
               " in the other)")
      } else {
        paste0("and with two lists of levels the function cannot take ",
               "the two sets of codes to be one scale (the labels the ",
               "two share have the same codes here, but a level present ",
               "in one list alone has no stated place in the other)")
      }
      .dmar_refuse(c(
        fun, "() needs ", both, " on one scale. Both are ordered ",
        "factors, but their levels differ: '", x_name, "' has ",
        .dominance_show(lev_x, " < "), " and '", y_name, "' has ",
        .dominance_show(lev_y, " < "), ". Each factor numbers its own ",
        "levels from 1, ", why, ". Nothing was computed. Give both the ",
        "same levels in the same order, for example lev <- c(...); ",
        "factor(", x_name, ", levels = lev, ordered = TRUE) for each."),
        "scores_not_ordered")
    }
  }
  # Two time differences in the same units are on one scale.
  if (kind_x == "difftime" && kind_y == "difftime" &&
      identical(units(x), units(y))) {
    kinds <- c("numeric", "numeric")
  }
  for (i in 1:2) {
    if (kinds[i] == "difftime") {
      name <- c(x_name, y_name)[i]
      .dmar_refuse(c("'", name, "' is of class difftime, and the other ",
                     "group is not a time difference in the same units. A ",
                     "time difference is stored in its own units (hours in ",
                     "one object, minutes in another), so the two are not ",
                     "on one scale. Convert each to numbers in one unit, ",
                     "for example as.numeric(", name, ", units = \"hours\"). ",
                     "Nothing was computed."),
                   "scores_not_ordered", arg = name)
    }
    if (kinds[i] == "other") {
      name <- c(x_name, y_name)[i]
      .dmar_refuse(c("'", name, "' must be a numeric or logical vector or ",
                     "an ordered factor; it is of class ",
                     paste(class(list(x, y)[[i]]), collapse = "/"), ". ",
                     fun, "() uses only the order of the observations, ",
                     "and an object of this class has no order the ",
                     "function can read. Nothing was computed."),
                   "scores_not_ordered", arg = name)
    }
  }
  if (all(c("date", "time") %in% kinds)) {
    .dmar_refuse(c("One group is of class Date and the other of class ",
                   "POSIXct, so the two are not on one scale (days against ",
                   "seconds). Supply both in one of the two classes. ",
                   "Nothing was computed."),
                 "scores_not_ordered")
  }

  x <- as.numeric(x)
  y <- as.numeric(y)
  if (isTRUE(paired)) {
    if (length(x) != length(y)) {
      .dmar_refuse_argument(c(
        "'", x_name, "' and '", y_name,
        "' must be the same length (paired observations)."), arg = NULL)
    }
    keep <- !is.na(x) & !is.na(y)
    n_missing <- sum(!keep)
    x <- x[keep]
    y <- y[keep]
  } else {
    miss_x <- is.na(x)
    miss_y <- is.na(y)
    n_missing <- c(group_1 = sum(miss_x), group_2 = sum(miss_y))
    x <- x[!miss_x]
    y <- y[!miss_y]
  }
  list(x = x, y = y, n_missing = n_missing)
}

# .dominance_placements(x, y)
#
# The pairwise comparisons of two complete numeric samples, from ranks, in
# place of the n_1 by n_2 matrix that outer() would build. With R the
# mid-rank of an observation in the pooled sample and R' its mid-rank in
# its own sample, R - R' counts the observations of the other sample below
# it, with those tied with it at one half. The same difference of minimum
# ranks counts those strictly below, and of maximum ranks those below or
# tied. Returned:
#   V_10       for each x_i, the mean over j of psi(x_i, y_j), where psi is
#              1, 1/2, 0 as x_i >, =, < y_j (the placement of DeLong,
#              DeLong, & Clarke-Pearson, 1988);
#   V_01       for each y_j, the mean over i of psi(x_i, y_j);
#   p_greater, p_tied, p_less
#              the proportions of the n_1 n_2 pairs with x_i >, =, < y_j.
# The counts are whole numbers held as doubles (a sum of integer ranks
# would overflow once n_1 n_2 passed the largest integer), so the three
# proportions are exact ratios, and are exactly 0 or 1 at complete
# separation and when every observation is tied. Infinite values are
# ranked as the largest or smallest values, tied with each other.
.dominance_placements <- function(x, y) {
  n_1 <- length(x)
  n_2 <- length(y)
  i_1 <- seq_len(n_1)
  i_2 <- n_1 + seq_len(n_2)
  pooled <- c(x, y)
  r_avg <- rank(pooled)
  r_min <- rank(pooled, ties.method = "min")
  r_max <- rank(pooled, ties.method = "max")
  n_less <- as.numeric(r_min[i_1] - rank(x, ties.method = "min"))
  n_leq  <- as.numeric(r_max[i_1] - rank(x, ties.method = "max"))
  M <- as.numeric(n_1) * n_2
  list(V_10      = (r_avg[i_1] - rank(x)) / n_2,
       V_01      = 1 - (r_avg[i_2] - rank(y)) / n_1,
       p_greater = sum(n_less) / M,
       p_tied    = sum(n_leq - n_less) / M,
       p_less    = (M - sum(n_leq)) / M)
}

# .dominance_degenerate_interval(fun, reason, situation, estimate, n,
#                                conf_level, test_reported = FALSE)
#
# The warning the three functions give when the estimated variance is 0,
# which happens exactly when the groups are completely separated (`reason`
# "complete_separation") or every observation has the same value
# ("no_variation"). The limits reported are then both equal to the
# estimate. `situation` is the opening clause that says which case holds,
# in the calling function's terms; `estimate` is the estimate on the
# function's own scale, and `n` is c(n_1, n_2) for two groups or the
# number of pairs. `test_reported` is TRUE when vargha_delaney_A() ran its
# permutation test, the one test defined in this case.
#
# The warning has class dmar_interval_degenerate and carries `reason`,
# `estimate`, `n`, and `limit_feng_cliff`, so that a simulation can muffle
# or count it by class. `limit_feng_cliff` is the limit Feng and Cliff
# (2004, p. 324, Equations 6 to 8) give for complete separation, on the
# function's own scale: with the largest variance delta-hat can have,
# (1 - delta^2) / n_b, the limit is (n_b - z^2) / (n_b + z^2) when
# delta-hat is 1 and its negative when delta-hat is -1, the other limit
# being delta-hat, with n_b the smaller group size and z the normal
# quantile of the level. On the scale of A = (delta + 1) / 2 that is
# n_b / (n_b + z^2) when A is 1 and z^2 / (n_b + z^2) when A is 0. It is
# stated in the message and is not the limit any of the functions
# reports; it is NA for the paired function and when there is no
# variation. Returns `reason`, invisibly, for the caller to record in the
# interval_degenerate attribute of its result.
.dominance_degenerate_interval <- function(fun, reason, situation, estimate,
                                           n, conf_level,
                                           test_reported = FALSE) {
  paired <- identical(fun, "probability_of_superiority_paired")
  on_A   <- identical(fun, "vargha_delaney_A")
  num    <- function(v) formatC(v, format = "f", digits = 3)
  est    <- format(estimate)
  test_sentence <- if (!on_A) {
    ""
  } else if (isTRUE(test_reported)) {
    if (identical(reason, "no_variation")) {
      paste0("The permutation test is reported (test_statistic is 0 and p_value ",
             "is 1); the t and z tests are not defined here. ")
    } else {
      paste0("The permutation test is reported (test_statistic is infinite); ",
             "the t and z tests are not defined here. ")
    }
  } else {
    "No test is reported (test_statistic and p_value are NA). "
  }
  limit <- NA_real_
  if (identical(reason, "complete_separation")) {
    limit_sentence <- ""
    if (!paired) {
      n_b <- min(n)
      z2  <- stats::qnorm(1 - (1 - conf_level) / 2)^2
      delta_hat   <- if (on_A) 2 * estimate - 1 else estimate
      limit_delta <- sign(delta_hat) * (n_b - z2) / (n_b + z2)
      limit <- if (on_A) (limit_delta + 1) / 2 else limit_delta
      limit_sentence <- paste0(
        "For this case Feng and Cliff (2004, Equation 8) give the limit ",
        if (delta_hat < 0) "-", "(n_b - z^2) / (n_b + z^2) for delta",
        if (on_A) " = 2A - 1",
        ", with n_b the smaller group size and z the normal quantile of ",
        "the level: ", num(limit_delta), " here at the ",
        .format_conf_pct(conf_level), "% level, with the other limit ",
        format(delta_hat),
        if (on_A) paste0("; on the scale of A, ", num(limit),
                         ", with the other limit ", est),
        ". ")
    }
    msg <- c(
      fun, "(): ", situation, ", so the estimated variance is 0. Both ",
      "confidence limits are reported equal to the estimate (", est, ", ",
      est, "). They are not a confidence interval with the stated ",
      "coverage: the variance estimate is built from how much the ",
      "comparisons differ from one ",
      if (paired) "pair" else "observation", " to the next, and in a ",
      "sample ",
      if (paired) "in which every pair favors the same member" else
        "with no overlap",
      " they do not differ at all. That describes this sample, not the ",
      "population. A population ",
      if (paired) {
        paste0("in which some pairs favor the other member often ",
               "produces a sample with none when there are few pairs ",
               "(here ", n, " pairs)")
      } else {
        paste0("with some overlap often produces a sample with none when ",
               "the samples are small (here ", n[1L], " and ", n[2L],
               " observations)")
      },
      ", so the limits must not be read as showing that the population ",
      "value is ", est, ". ", limit_sentence, test_sentence,
      "Report the estimate with the sample ",
      if (paired) {
        paste0("size and say that every pair favored the same member; ",
               "the randomization test (randomization_test_paired()) ")
      } else {
        paste0("sizes and say that the groups did not overlap; the ",
               "randomization test (randomization_test()) ")
      },
      "gives a p-value that remains valid here.")
  } else {
    msg <- c(
      fun, "(): ", situation, ", so the estimated variance is 0. The ",
      "estimate is ", est, " (no difference in either direction) and ",
      "both confidence limits are reported equal to it (", est, ", ", est,
      "). They are not a confidence interval with the stated coverage: ",
      "with no variation in the data nothing can be learned about how ",
      "far the population value is from ", est, ", and an interval of no ",
      "width states the opposite. ", test_sentence,
      "Check that the scores were read as intended (a constant column, a ",
      "scale with one category in use).")
  }
  .dmar_warn(msg, "interval_degenerate", reason = reason,
             estimate = estimate, n = n, limit_feng_cliff = limit)
  invisible(reason)
}

# The line a degenerate table of the three functions prints beneath
# itself, or NULL. The class methods of .dmar_tbl_footer() in the three
# function files add it.
.dominance_degenerate_footer <- function(x) {
  if (is.null(attr(x, "interval_degenerate", exact = TRUE))) return(NULL)
  "Limits are not a confidence interval: the estimated variance is 0."
}

# Adds the logical column interval_degenerate to the tidy() or glance()
# rows `out` of the table `x`: TRUE on every row when the table's limits
# are not a confidence interval. The attributes tidy.dmar_tbl() sets on
# its result (p_terms, p_sentence) are kept.
.dominance_mark_degenerate <- function(out, x) {
  out$interval_degenerate <- rep(
    !is.null(attr(x, "interval_degenerate", exact = TRUE)), nrow(out))
  out
}

# .dominance_permutation_W(block_start, block_end, G)
#
# The statistic W = (A - 1/2) / se of vargha_delaney_A() for each column
# of `G`, an N by m matrix of 0 and 1 whose rows follow the pooled scores
# in increasing order and whose columns mark the scores assigned to the
# first group (n_1 ones in every column). `block_start` and `block_end`
# give, for each row, the first and last row of the block of scores tied
# with it.
#
# Everything is counted in whole numbers. With c the doubled placement
# count of a first-group score (twice the number of second-group scores
# below it, plus the number tied with it), A = sum(c) / (2 n_1 n_2) and
# the placement variance S2_10 = {n_1 sum(c^2) - sum(c)^2} /
# {n_1 (n_1 - 1) (2 n_2)^2}, and likewise for the second group. The sums
# are exact (they stay far below 2^53 for any design that can be run),
# so a variance of 0 is found exactly, and two reassignments that give
# the same counts give the same W to the last bit. A reassignment whose
# estimated variance is 0 has W = 0 when A = 1/2 (every score equal) and
# W = Inf with the sign of A - 1/2 otherwise (complete separation).
.dominance_permutation_W <- function(block_start, block_end, G) {
  N   <- nrow(G)
  n_1 <- sum(G[, 1L])
  n_2 <- N - n_1
  # Running counts of first-group scores at or before each row, with a
  # row of zeros in front, so that row k + 1 holds the count through k.
  cum_1 <- matrix(0, nrow = N + 1L, ncol = ncol(G))
  for (k in seq_len(N)) cum_1[k + 1L, ] <- cum_1[k, ] + G[k, ]
  below_1 <- cum_1[block_start, , drop = FALSE]        # strictly below
  upto_1  <- cum_1[block_end + 1L, , drop = FALSE]     # below or tied
  below_2 <- (block_start - 1) - below_1
  upto_2  <- block_end - upto_1
  # Doubled placements: of a first-group score among the second group,
  # and of a second-group score among the first (first-group scores
  # above it, plus half of those tied).
  c_10 <- (below_2 + upto_2) * G
  c_01 <- (2 * n_1 - below_1 - upto_1) * (1 - G)
  sum_10 <- colSums(c_10)
  sum_01 <- colSums(c_01)
  ss_10  <- n_1 * colSums(c_10 * c_10) - sum_10^2
  ss_01  <- n_2 * colSums(c_01 * c_01) - sum_01^2
  M      <- as.numeric(n_1) * n_2
  num    <- (sum_10 - M) / (2 * M)                      # A - 1/2
  var_A  <- ss_10 / (n_1 * (n_1 - 1) * (2 * n_2)^2) / n_1 +
            ss_01 / (n_2 * (n_2 - 1) * (2 * n_1)^2) / n_2
  W <- num / sqrt(var_A)
  flat <- ss_10 == 0 & ss_01 == 0
  W[flat] <- sign(num[flat]) * Inf
  W[flat & num == 0] <- 0
  W
}

# .dominance_permutation_test(x, y, exact, n_resamples)
#
# The studentized permutation test of A = 1/2 for vargha_delaney_A(): the
# statistic W over reassignments of the pooled scores to groups of sizes
# n_1 and n_2. The choice between enumeration and random draws, the
# thresholds, the form of the p-value, and the counting rule are those of
# randomization_test() (R/randomization_test_internals.R): every
# reassignment when choose(N, n_1) is at most .rt_exact_threshold (or
# when `exact` is TRUE, up to .rt_exact_cap), `n_resamples` random
# reassignments otherwise; the share of reassignments with |W_b| >= |W|
# under enumeration and (1 + count) / (1 + n_resamples) under random
# draws, counted by .rt_extreme(). The caller has set the seed. The
# observed W goes through the same arithmetic as the reassignments, and
# the comparison allows a rounding bound of a few units in the last
# place of the two values, so that reassignments tied with the observed
# one are counted. Returns list(W, p_value, exact, n_arrangements,
# n_resamples_used).
.dominance_permutation_test <- function(x, y, exact, n_resamples) {
  n_1 <- length(x)
  n_2 <- length(y)
  N   <- n_1 + n_2
  pooled <- c(x, y)
  ord    <- order(pooled)
  sorted <- pooled[ord]
  # The block of tied scores each row belongs to.
  block_start <- match(sorted, sorted)
  block_end   <- N + 1L - match(sorted, rev(sorted))
  position    <- integer(N)
  position[ord] <- seq_len(N)          # row of each pooled score

  membership <- function(perm) {
    m <- ncol(perm)
    G <- matrix(0, nrow = N, ncol = m)
    G[cbind(position[as.vector(perm)], rep(seq_len(m), each = n_1))] <- 1
    G
  }
  W_obs <- .dominance_permutation_W(
    block_start, block_end, membership(matrix(seq_len(n_1), ncol = 1L)))

  n_possible <- choose(N, n_1)
  use_exact  <- if (is.null(exact)) n_possible <= .rt_exact_threshold else
    isTRUE(exact)
  if (use_exact && n_possible > .rt_exact_cap) {
    .dmar_refuse_argument(c(
      "Exact enumeration would require ",
      format(n_possible, big.mark = ",", scientific = FALSE),
      " reassignments, above the ",
      format(.rt_exact_cap, big.mark = ",", scientific = FALSE),
      " cap. Use exact = FALSE for the Monte Carlo test."), arg = NULL)
  }
  extreme <- function(perm) {
    W_null <- .dominance_permutation_W(block_start, block_end,
                                       membership(perm))
    tol <- 64 * .Machine$double.eps * (abs(W_null) + abs(W_obs))
    .rt_extreme(W_null, W_obs, "two_sided", tol)
  }
  # The reassignments are taken in slices, so that the N by m matrices of
  # a slice stay near two million cells.
  slice <- max(1L, floor(2e6 / N))
  if (use_exact) {
    perm   <- utils::combn(N, n_1)
    n_eval <- ncol(perm)
    count  <- 0
    for (from in seq(1L, n_eval, by = slice)) {
      to <- min(n_eval, from + slice - 1L)
      count <- count + extreme(perm[, from:to, drop = FALSE])
    }
    p_value <- count / n_eval
  } else {
    n_eval <- n_resamples
    count  <- 0
    done   <- 0L
    while (done < n_eval) {
      m <- min(slice, n_eval - done)
      perm  <- matrix(replicate(m, sample.int(N, n_1)), nrow = n_1)
      count <- count + extreme(perm)
      done  <- done + m
    }
    p_value <- (1 + count) / (1 + n_eval)
  }
  list(W = W_obs, p_value = p_value, exact = use_exact,
       n_arrangements = n_possible, n_resamples_used = n_eval)
}
