# Suppress R CMD CHECK notes for ggplot2 aes() column references.
if (getRversion() >= "2.15.1") {
  utils::globalVariables(c("statistic_value", "region"))
}

# The randomization distribution behind a randomization_test() result.
#' Plot the Randomization Distribution Behind a Randomization Test
#'
#' Displays the reference distribution that \code{\link{randomization_test}}
#' built by reassigning the observed scores to the two groups, with the
#' observed statistic marked and every reassignment at least as extreme as
#' the observed one shaded. The shaded reassignments are the ones the
#' test counted, so under exact enumeration the shaded proportion is the
#' \emph{p}-value, and the figure shows where that number came from
#' instead of only reporting it.
#'
#' Reading the figure is the point of it. The spread of the distribution is
#' what the reassignments alone can produce when the grouping is irrelevant,
#' which is the null hypothesis of the test. If the observed statistic sits
#' inside that spread, reassignment alone explains it. If it sits out in a
#' tail, few reassignments reproduce it, and that scarcity is the evidence.
#' No normal or \emph{t} distribution appears anywhere in the construction.
#' This is the display Chapter 1 of Maxwell, Delaney, and Kelley (2027) uses
#' to introduce the logic of the randomization test.
#'
#' @param object A result of \code{\link{randomization_test}}. A
#'   reassignment whose statistic is infinite has no place on the axis;
#'   the histogram shows the finite values and the subtitle states how
#'   many reassignments are not drawn. Infinite values occur in the
#'   reference distribution only, under \code{statistic = "t"}, for a
#'   reassignment that leaves both groups constant at different values
#'   (with a two-valued outcome, the one that puts every high score in one
#'   group). The observed statistic is always finite, because
#'   \code{\link{randomization_test}} refuses \code{statistic = "t"} when
#'   the observed groups are themselves constant at different values, so
#'   the line that marks it is always drawn.
#' @param bins Number of histogram bins used to display the reference
#'   distribution. Defaults to \code{40}.
#' @param palette Character; the color palette. Defaults to \code{"okabe_ito"},
#'   base R's colorblind-safe Okabe-Ito palette; \code{"tableau"} is also
#'   available.
#'
#' @return A \code{ggplot} object, which can be printed or further modified
#'   with the usual \pkg{ggplot2} verbs.
#'
#' @author Ken Kelley
#'
#' @references
#' Fisher, R. A. (1935). \emph{The design of experiments}. Oliver & Boyd.
#'
#' Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). \emph{Designing
#'   experiments and analyzing data: A model comparison perspective}
#'   (4th ed.). Routledge. (See Chapter 1 on the logic of the randomization
#'   test.)
#'
#' @seealso \code{\link{randomization_test}} for the test itself.
#'
#' @family plotting
#'
#' @keywords hplot
#'
#' @examples
#' treatment <- c(80, 84, 79, 88, 83)
#' control   <- c(72, 75, 68, 81, 74)
#' rt <- randomization_test(group_1 = treatment, group_2 = control)
#' plot_randomization_test(rt)
#'
#' @export
plot_randomization_test <- function(object, bins = 40L,
                                    palette = "okabe_ito") {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    .dmar_refuse(c(
      "The 'ggplot2' package is required for plot_randomization_test(). ",
      "Install it with install.packages(\"ggplot2\")."),
      "package_required", package = "ggplot2")
  }

  T_null <- attr(object, "reference_distribution")
  T_obs  <- attr(object, "observed_statistic")
  if (is.null(T_null) || is.null(T_obs)) {
    .dmar_refuse_argument(c(
      "'object' does not carry a reference distribution. Supply a result ",
      "of randomization_test()."), arg = "object")
  }
  if (!is.numeric(bins) || length(bins) != 1L || !is.finite(bins) ||
      bins < 1) {
    .dmar_refuse_argument(
      "'bins' must be a single positive number.", arg = "bins")
  }

  alternative <- attr(object, "alternative")
  stat_name   <- attr(object, "statistic_name")
  method      <- attr(object, "method")
  if (is.null(alternative)) alternative <- "two_sided"
  if (is.null(stat_name))   stat_name   <- "test statistic"

  # The shaded reassignments are the ones the p-value counted, which the
  # result carries, so the shaded area and the reported p-value are the
  # same quantity. A result saved by a version that did not store them is
  # shaded by the counting rule with a tolerance relative to the largest
  # finite statistic.
  # Objects saved before the vocabulary moved to snake_case may carry the
  # dotted spelling in their attribute; normalize before dispatching.
  if (identical(alternative, "two.sided")) alternative <- "two_sided"
  extreme <- attr(object, "reference_extreme")
  if (is.null(extreme) || length(extreme) != length(T_null)) {
    tol <- sqrt(.Machine$double.eps) *
      max(abs(T_obs[is.finite(T_obs)]), abs(T_null[is.finite(T_null)]), 0)
    extreme <- .rt_is_extreme(T_null, T_obs, alternative, tol)
  }
  drawn      <- is.finite(T_null)
  n_infinite <- sum(!drawn)

  p_value <- object$value[object$term == "p_value"]
  p_value <- if (length(p_value) == 1L) p_value else NA_real_

  plot_data <- data.frame(
    statistic_value = as.numeric(T_null[drawn]),
    region = factor(ifelse(extreme[drawn], "At least as extreme",
                           "Less extreme"),
                    levels = c("Less extreme", "At least as extreme"))
  )

  cols <- .dmar_palette(2, palette = palette)
  fills <- stats::setNames(c("grey75", cols[2L]),
                           c("Less extreme", "At least as extreme"))

  n_ref <- length(T_null)
  # When every reassignment the p-value counted has an infinite statistic,
  # nothing that is drawn is shaded, and the subtitle does not call p a
  # shaded area. Under exact enumeration the observed assignment, whose
  # statistic is finite, is always among those counted, so this arises
  # only in a Monte Carlo sample that holds a counted infinite value and
  # no counted finite one.
  all_counted_infinite <- n_infinite > 0L && sum(extreme) > 0L &&
    !any(extreme[drawn])
  subtitle <- sprintf(
    "%s over %s reassignments%s",
    if (identical(method, "exact enumeration")) "Exact enumeration"
    else "Monte Carlo sampling",
    format(n_ref, big.mark = ","),
    if (is.na(p_value)) "" else sprintf(
      if (all_counted_infinite) "; p = %s" else "; shaded area is p = %s",
      format_p(p_value))
  )
  if (n_infinite > 0L) {
    # On a line of its own, so that the subtitle fits a figure of ordinary
    # width.
    subtitle <- paste0(
      subtitle, "\n", format(n_infinite, big.mark = ","),
      if (n_infinite == 1L) " reassignment" else " reassignments",
      " with an infinite statistic ",
      if (n_infinite == 1L) "is" else "are", " not drawn",
      if (all_counted_infinite && !is.na(p_value))
        "; p counts only these, so no area is shaded")
  }

  out <- ggplot2::ggplot(plot_data,
                         ggplot2::aes(x = statistic_value, fill = region)) +
    ggplot2::geom_histogram(bins = bins, color = "white", linewidth = 0.15)
  # The observed statistic is finite (randomization_test() refuses the one
  # case in which it is not), so its line is always drawn.
  out +
    ggplot2::geom_vline(xintercept = T_obs, linewidth = 0.9,
                        color = cols[1L]) +
    ggplot2::annotate("text", x = T_obs, y = Inf, label = " observed",
                      hjust = 0, vjust = 1.6, size = 3.3,
                      color = cols[1L]) +
    ggplot2::scale_fill_manual(values = fills, name = NULL) +
    ggplot2::labs(
      title = "Randomization distribution of the test statistic",
      subtitle = subtitle,
      x = paste0("Reassignment value of the ", stat_name),
      y = "Number of reassignments"
    ) +
    ggplot2::theme(legend.position = "top")
}
