# =============================================================================
#  Internal helpers for the minimum risk functions, mr_cv() and mr_smd()
# =============================================================================
#
# Nothing in this file is exported.
#
# Layout:
#   (A) The checks of the arguments the two functions share: .mr_check_flag,
#       .mr_check_positive, .mr_resolve_A, .mr_check_gamma
#   (B) The variance estimate and the stopping rule of mr_cv():
#       .mr_cv_check_data, .mr_cv_estimate, .mr_cv_rule
#   (C) The decision line under a table of either function:
#       .dmar_tbl_footer.dmar_mr, .mr_decision_line

# --- (A) Shared argument checks ---------------------------------------------

# `pilot` and `verbose`: TRUE or FALSE, and nothing that R would coerce to
# one of them (1, 0, "TRUE").
.mr_check_flag <- function(x, arg) {
  if (is.logical(x) && length(x) == 1L && !is.na(x)) return(invisible(x))
  .dmar_refuse_argument(c("'", arg, "' must be TRUE or FALSE."), arg = arg,
                        value = x)
}

# A cost, a price, or a maximum probable error: one positive finite number.
# `what` says what the quantity is, for the message.
.mr_check_positive <- function(x, arg, what) {
  if (is.numeric(x) && length(x) == 1L && is.finite(x) && x > 0)
    return(invisible(x))
  .dmar_refuse_argument(c(
    "'", arg, "' must be a single positive finite number; ", .dmar_given(x),
    ". It is ", what, "."), arg = arg, value = x)
}

# A, given directly or as structural_cost / epsilon^2. `given` is a named
# logical vector, from missing() in the calling function, saying which of
# the three arguments the call supplied; an argument that was not supplied
# is not evaluated here.
.mr_resolve_A <- function(A, structural_cost, epsilon, given) {
  if (given[["A"]]) {
    if (given[["structural_cost"]]) .dmar_refuse_argument(
      "Because you specified 'A' directly, you should not also specify 'structural_cost'.",
      arg = c("A", "structural_cost"))
    if (given[["epsilon"]]) .dmar_refuse_argument(
      "Because you specified 'A' directly, you should not also specify 'epsilon'.",
      arg = c("A", "epsilon"))
    .mr_check_positive(A, "A", paste0(
      "the structural cost divided by the square of 'epsilon', the price ",
      "of a squared unit of maximum probable error"))
    return(A)
  }
  if (!given[["structural_cost"]]) .dmar_refuse_argument(
    "Because you did not specify 'A' directly, you must specify 'structural_cost'.",
    arg = c("A", "structural_cost"))
  if (!given[["epsilon"]]) .dmar_refuse_argument(
    "Because you did not specify 'A' directly, you must specify 'epsilon'.",
    arg = c("A", "epsilon"))
  .mr_check_positive(structural_cost, "structural_cost",
                     "the cost of the study apart from the cost of sampling")
  .mr_check_positive(epsilon, "epsilon", paste0(
    "the largest difference between the estimate and the population value ",
    "that is acceptable"))
  structural_cost / epsilon^2
}

# A / sampling_cost enters the pilot sample size and the stopping rule; it
# must be a finite number (A = 1e300 with sampling_cost = 1e-10 is not).
.mr_check_ratio <- function(A, sampling_cost) {
  if (is.finite(A / sampling_cost)) return(invisible(TRUE))
  .dmar_refuse_argument(c(
    "'A' divided by 'sampling_cost' must be a finite number; ",
    format(A), " / ", format(sampling_cost), " is not. The pilot sample ",
    "size and the stopping rule both use that ratio."),
    arg = c("A", "sampling_cost"))
}

# The pilot sample size, max(m0, ceiling(x)), where x is the formula of the
# procedure ((A / c)^(1 / (2 + 2 gamma)) for mr_cv(), with A / (2c) for
# mr_smd()). Sample sizes are whole numbers, so the formula is rounded up,
# and a value that is a whole number in exact arithmetic must stay that
# number. In floating point it can come out a little above it: (1e5)^0.4 is
# 100.00000000000003, whose ceiling is 101. x is therefore lowered by a
# relative 1e-10 before the ceiling, far less than any formula value that
# is not a whole number can lie above one.
.mr_pilot_size <- function(x, m0) {
  max(m0, ceiling(x * (1 - 1e-10)))
}

# The exponent of the correction term: strictly between 0 and 1/2. `why` is
# the sentence that says where the interval comes from.
.mr_check_gamma <- function(gamma, why) {
  if (is.numeric(gamma) && length(gamma) == 1L && is.finite(gamma) &&
      gamma > 0 && gamma < 0.5) return(invisible(gamma))
  .dmar_refuse_argument(c(
    "'gamma' must be a single number strictly between 0 and 1/2; ",
    .dmar_given(gamma), ". ", why, " The default is .49."),
    arg = "gamma", value = gamma)
}

# --- (B) The variance estimate and the stopping rule of mr_cv() --------------

# The scores of mr_cv(): a numeric vector of at least 4 finite observations
# whose mean is positive. The mean is compared with a multiple of the
# machine precision at the scale of the scores, so that centered scores,
# whose mean is zero to rounding, are refused with scores whose mean is
# negative.
.mr_cv_check_data <- function(data) {
  shape_ok <- is.numeric(data) && is.null(dim(data))
  if (!shape_ok || length(data) < 4L || any(!is.finite(data))) {
    given <- if (!shape_ok) {
      sprintf("an object of class %s was given", class(data)[1L])
    } else if (length(data) < 4L) {
      sprintf("%d %s given", length(data),
              if (length(data) == 1L) "was" else "were")
    } else {
      n_bad <- sum(!is.finite(data))
      sprintf("%d of the %d values %s missing or infinite", n_bad,
              length(data), if (n_bad == 1L) "is" else "are")
    }
    .dmar_refuse_argument(c(
      "'data' must be a numeric vector of at least 4 finite observations; ",
      given, ". The variance estimate uses the fourth central moment, ",
      "whose unbiased estimate needs four observations."), arg = "data")
  }
  m <- mean(data)
  if (!(m > 10 * .Machine$double.eps * max(abs(data)))) {
    .dmar_refuse_argument(c(
      "'data' must have a positive mean; the mean of the ", length(data),
      " observations is ", format(m, digits = 4),
      if (m > 0) ", which is zero to rounding at the scale of the scores",
      ". The coefficient of ",
      "variation divides by the mean, and the procedure is for a quantity ",
      "that cannot be negative (Chattopadhyay & Kelley, 2016, p. 627)."),
      arg = "data", value = m)
  }
  invisible(data)
}

# The estimate of xi^2, the variance of the limiting distribution of
# sqrt(n) (k_n - kappa): Equation 32 of Kelley, Darku, and Chattopadhyay
# (2018),
#   V_n^2 = mu4 / (4 m^2 s^2) - s^2 / (4 m^2) - mu3 / m^3 + s^4 / m^4,
# with m the sample mean, s^2 the sample variance, and mu3 and mu4 the
# unbiased estimates of the third and fourth central moments (their
# Equations 78 and 79). Every sum is taken over the centered scores: the
# fourth moment written in raw power sums, as Equation 79 prints it, loses
# its digits to cancellation when the mean is large beside the standard
# deviation. For scores that do not vary every term but the first is 0 and
# the first is 0 / 0; the estimate is then taken as 0, the limit of
# Equation 32 as the standard deviation goes to 0 for scores of any shape.
# Scores whose range is at most 1e-10 of their largest absolute value
# differ only by rounding error (equal gains computed as posttest minus
# pretest, for example) and are taken as not varying, and an estimate
# within 1e-10 of the sum of the absolute values of its four terms is zero
# to rounding; either is taken as 0.
#
# Returns the pieces mr_cv() reports: n, the mean, s^2, mu3, mu4, and the
# estimate itself, v2, before any bound.
.mr_cv_estimate <- function(data) {
  n <- length(data)
  m <- mean(data)
  d <- data - m
  ss <- sum(d^2)
  s2 <- ss / (n - 1)
  mu3 <- n * sum(d^3) / ((n - 1) * (n - 2))
  mu4 <- ((n^2 - 2 * n + 3) * sum(d^4) - 3 * (2 * n - 3) / n * ss^2) /
    ((n - 1) * (n - 2) * (n - 3))
  varies <- s2 > 0 && (max(data) - min(data)) > 1e-10 * max(abs(data))
  v2 <- if (varies) {
    mu4 / (4 * m^2 * s2) - s2 / (4 * m^2) - mu3 / m^3 + s2^2 / m^4
  } else {
    0
  }
  if (varies) {
    size <- abs(mu4 / (4 * m^2 * s2)) + s2 / (4 * m^2) + abs(mu3 / m^3) +
      s2^2 / m^4
    if (abs(v2) <= 1e-10 * size) v2 <- 0
  }
  list(n = n, mean = m, s2 = s2, mu3 = mu3, mu4 = mu4, v2 = v2)
}

# The stopping rule of Chattopadhyay and Kelley (2016, Equation 23) and the
# estimated risk (their Equation 14 with the estimate for xi^2), given the
# estimate v2 at sample size n.
#
# A usable estimate is positive and finite. It enters the rule and the risk
# as max(v2, n^-3), the estimator of Equation 31 of Kelley, Darku, and
# Chattopadhyay (2018). An estimate that is zero, negative, or not finite
# is unusable: the rule is reported as not met and the risk and the
# criterion are NA. Equation 31 would put n^-3 in its place, and the rule
# can then be met on no evidence about the variance (for the scores
# 1, 1, 1, 1, 1, 1, 3, 3 with A = 400 and c = 1 the criterion from n^-3 is
# 52.9, below n^2 = 64).
.mr_cv_rule <- function(v2, n, A, sampling_cost, gamma) {
  if (!is.finite(v2) || v2 <= 0) {
    return(list(unusable = TRUE, risk = NA_real_, criterion = NA_real_,
                is_satisfied = FALSE))
  }
  used <- max(v2, n^-3)
  criterion <- (A / sampling_cost) * (used + n^(-2 * gamma))
  list(unusable = FALSE, risk = A * used / n + sampling_cost * n,
       criterion = criterion, is_satisfied = n^2 >= criterion)
}

# --- (C) The decision line ---------------------------------------------------

# A table of mr_cv() or mr_smd() says beneath itself whether the stopping
# rule is met, at the console and in a knitted report alike (see
# .dmar_tbl_footer() in R/dmar_output_helpers.R). The table of a pilot
# sample size carries no decision and prints as an ordinary DMAR table. So
# does a selection that drops the attributes (columns chosen with `[`, or
# subset()). A selection of rows alone (x[i, ], head()) keeps the
# attributes and the line, which stays true of those rows; without the row
# the line compares (N, or n_1) it prints no line, unless the variance
# estimate was unusable, whose line names no row.
#' @exportS3Method .dmar_tbl_footer dmar_mr
#' @noRd
.dmar_tbl_footer.dmar_mr <- function(x, digits = getOption("dmar.digits", 3L), ...) {
  .dmar_footer_add(NextMethod(), .mr_decision_line(x, digits = digits),
                   "stopping_rule")
}

# The line itself, or NULL when the table does not carry a decision. The
# two functions are told apart by their rows: mr_cv() reports N and
# compares N^2 with its criterion; mr_smd() reports n_1 and n_2 and
# compares the sample size per group with its criterion. The criterion
# prints to `digits` decimal places.
.mr_decision_line <- function(x, digits = getOption("dmar.digits", 3L)) {
  met <- attr(x, "is_satisfied", exact = TRUE)
  term <- x[["term"]]
  value <- x[["value"]]
  if (!is.logical(met) || length(met) != 1L || is.na(met) ||
      is.null(term) || !is.numeric(value)) return(NULL)
  per_group <- "n_1" %in% term
  more <- if (per_group) "collect another observation in each group." else
    "collect another observation."
  if (isTRUE(attr(x, "variance_unusable", exact = TRUE))) {
    return(paste(
      "The variance estimate is zero or negative at this sample size, so",
      "the stopping rule is reported as not met;", more))
  }
  criterion <- attr(x, "criterion", exact = TRUE)
  size <- value[match(if (per_group) "n_1" else "N", term)]
  if (!is.numeric(criterion) || length(criterion) != 1L ||
      !is.finite(criterion) || !is.finite(size)) return(NULL)
  whole <- function(v) format(v, scientific = FALSE, trim = TRUE,
                              big.mark = "")
  compared <- if (per_group) {
    sprintf("the sample size per group, %s,", whole(size))
  } else {
    sprintf("N^2 = %s", whole(size^2))
  }
  sprintf("Stopping rule: %s (%s is %s the criterion, %s); %s",
          if (met) "met" else "not met", compared,
          if (met) "at least" else "below",
          .dmar_fixed_text(criterion, digits),
          if (met) "sampling can stop." else more)
}
