# Internal helpers shared by the polynomial-change-model planners
# (ss_power_pcm) and the longitudinal polynomial simulator
# (simulate_longitudinal_polynomial). Not exported.
#
# Everything here rests on Raudenbush and Liu (2001, Psychological Methods, 6,
# 387-401). A "polynomial change model" describes each subject's repeated
# measurements with a degree-p polynomial in time; the degree-p coefficient is
# the change parameter of interest (p = 0 the intercept / flat line, p = 1 the
# linear slope, p = 2 the quadratic acceleration, and so on). The change
# coefficient is taken in the DERIVATIVE-scaled metric (p! times the leading
# coefficient of t^p), which is the metric in which the Raudenbush-Liu K_p
# constants below apply.

# Named polynomial trends indexed by order; order 0 (intercept) is a flat line.
.pcm_trend_names <- c("intercept", "linear", "quadratic", "cubic", "quartic",
                      "quintic", "sextic", "septic", "octic", "nonic", "decic")

# Resolve a user-supplied `trend` -- either a name ("quadratic") or a
# non-negative integer order (2) -- to its integer polynomial order p.
.pcm_trend_order <- function(trend) {
  if (is.character(trend)) {
    p <- match(tolower(trimws(trend)), .pcm_trend_names) - 1L
    if (is.na(p))
      stop("Unrecognized 'trend' name \"", trend, "\". Use one of: ",
           paste(.pcm_trend_names, collapse = ", "),
           ", or a non-negative integer order (0 = intercept / flat line, ",
           "1 = linear, 2 = quadratic, ...).", call. = FALSE)
    return(p)
  }
  if (is.numeric(trend) && length(trend) == 1L && !is.na(trend) &&
      trend >= 0 && trend == round(trend))
    return(as.integer(trend))
  stop("'trend' must be a single non-negative integer order ",
       "(0 = intercept / flat line, 1 = linear, 2 = quadratic, ...) or one of ",
       "the names: ", paste(.pcm_trend_names, collapse = ", "), ".",
       call. = FALSE)
}

# Human-readable label for an integer order.
.pcm_trend_label <- function(p) {
  if (p + 1L <= length(.pcm_trend_names)) .pcm_trend_names[p + 1L] else
    paste0("order_", p)
}

# Within-subject sampling variance V of the ordinary-least-squares estimate of
# the degree-p (derivative-scaled) polynomial change coefficient, given M
# equally spaced measurement occasions, level-one error variance sigma2_e, and
# `frequency` measurements per time unit (Raudenbush & Liu, 2001, p. 392):
#
#     V = sigma2_e * f^(2p) * (M - p - 1)! / ( K_p * (M + p)! ),
#
# where the Raudenbush-Liu constant satisfies
#
#     1 / K_p = (2p)! (2p+1)! / (p!)^2   (= 1, 12, 720, 100800, ... for
#                                          p = 0, 1, 2, 3, ...),
#
# so equivalently
#
#     V = sigma2_e * f^(2p) * (M - p - 1)! (2p)! (2p+1)! / ( (p!)^2 (M + p)! ).
#
# This V equals (p!)^2 times the variance of the OLS estimate of the
# coefficient of t^p (verified against (X'X)^{-1} for every order), reduces to
# sigma2_e / M at p = 0 (the variance of the mean) and to
# 12 sigma2_e f^2 / [M(M^2 - 1)] at p = 1 (the variance of the OLS slope). The
# combinatorial part is evaluated on the log scale so it stays finite for high
# orders and many occasions.
.pcm_sampling_variance <- function(sigma2_e, frequency, M, p) {
  if (M < p + 1L)
    stop(sprintf(paste0("A degree-%d (%s) change coefficient needs at least ",
                        "%d measurement occasions, but the design has M = %g. ",
                        "Increase 'frequency' or 'duration'."),
                 p, .pcm_trend_label(p), p + 1L, M), call. = FALSE)
  log_comb <- lfactorial(M - p - 1) + lfactorial(2 * p) + lfactorial(2 * p + 1) -
              2 * lfactorial(p) - lfactorial(M + p)
  sigma2_e * frequency^(2 * p) * exp(log_comb)
}
