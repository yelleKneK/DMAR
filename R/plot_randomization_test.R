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
#' the observed one shaded. The shaded proportion is the \emph{p}-value, so
#' the figure shows where that number came from instead of only reporting
#' it.
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
#' @param object A result of \code{\link{randomization_test}}.
#' @param bins Number of histogram bins used to display the reference
#'   distribution. Defaults to \code{40}.
#' @param palette Character; the color palette. Defaults to \code{"okabe_ito"},
#'   base R's colorblind-safe Okabe-Ito palette; \code{"tableau"} is also
#'   available.
#' @param ... Currently unused; present so the signature can grow without
#'   breaking existing calls.
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
                                    palette = "okabe_ito", ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("The 'ggplot2' package is required for plot_randomization_test(). ",
         "Install it with install.packages(\"ggplot2\").", call. = FALSE)
  }

  T_null <- attr(object, "reference_distribution")
  T_obs  <- attr(object, "observed_statistic")
  if (is.null(T_null) || is.null(T_obs)) {
    stop("'object' does not carry a reference distribution. Supply a result ",
         "of randomization_test().", call. = FALSE)
  }
  if (!is.numeric(bins) || length(bins) != 1L || !is.finite(bins) ||
      bins < 1) {
    stop("'bins' must be a single positive number.", call. = FALSE)
  }

  alternative <- attr(object, "alternative")
  stat_name   <- attr(object, "statistic_name")
  method      <- attr(object, "method")
  if (is.null(alternative)) alternative <- "two_sided"
  if (is.null(stat_name))   stat_name   <- "test statistic"

  # The same "at least as extreme" rule the p-value used, so the shaded
  # area and the reported p-value are the same quantity.
  tol <- 1e-12
  # Objects saved before the vocabulary moved to snake_case may carry the
  # dotted spelling in their attribute; normalize before dispatching.
  if (identical(alternative, "two.sided")) alternative <- "two_sided"
  extreme <- switch(
    alternative,
    two_sided = abs(T_null) >= abs(T_obs) - tol,
    greater   = T_null >= T_obs - tol,
    less      = T_null <= T_obs + tol
  )

  p_value <- object$value[object$term == "p_value"]
  p_value <- if (length(p_value) == 1L) p_value else NA_real_

  plot_data <- data.frame(
    statistic_value = as.numeric(T_null),
    region = factor(ifelse(extreme, "At least as extreme", "Less extreme"),
                    levels = c("Less extreme", "At least as extreme"))
  )

  cols <- .dmar_palette(2, palette = palette)
  fills <- stats::setNames(c("grey75", cols[2L]),
                           c("Less extreme", "At least as extreme"))

  n_ref <- length(T_null)
  subtitle <- sprintf(
    "%s over %s reassignments%s",
    if (identical(method, "exact enumeration")) "Exact enumeration"
    else "Monte Carlo sampling",
    format(n_ref, big.mark = ","),
    if (is.na(p_value)) "" else sprintf("; shaded area is p = %s",
                                        format_p(p_value))
  )

  ggplot2::ggplot(plot_data,
                  ggplot2::aes(x = statistic_value, fill = region)) +
    ggplot2::geom_histogram(bins = bins, color = "white", linewidth = 0.15) +
    ggplot2::geom_vline(xintercept = T_obs, linewidth = 0.9,
                        color = cols[1L]) +
    ggplot2::annotate("text", x = T_obs, y = Inf, label = " observed",
                      hjust = 0, vjust = 1.6, size = 3.3, color = cols[1L]) +
    ggplot2::scale_fill_manual(values = fills, name = NULL) +
    ggplot2::labs(
      title = "Randomization distribution of the test statistic",
      subtitle = subtitle,
      x = paste0("Reassignment value of the ", stat_name),
      y = "Number of reassignments"
    ) +
    ggplot2::theme(legend.position = "top")
}
