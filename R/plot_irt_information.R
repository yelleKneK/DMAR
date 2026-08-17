# Suppress R CMD check notes for the ggplot2 aes() column references below.
if (getRversion() >= "2.15.1") {
  utils::globalVariables(c("theta", "information", "item", "se_rescaled"))
}

#' Plot an Item Response Theory Information Curve
#'
#' Draws the information function computed by
#' \code{\link{irt_information}}: either the test information curve, with
#' the standard error of the latent trait estimate on a secondary axis, or
#' one curve per item. The test view answers "where on the latent
#' continuum does this scale measure precisely?", and because the standard
#' error is \eqn{1 / \sqrt{I(\theta)}} the same picture shows the precision
#' directly. The item view decomposes that curve, since information is
#' additive across items, and so shows which items cover which part of the
#' continuum.
#'
#' @param x The result of \code{\link{irt_information}}.
#' @param what Which curves to draw: \code{"test"} (default) for the test
#'   information function, or \code{"item"} for one curve per item.
#' @param show_se Logical. When \code{TRUE} (the default) and
#'   \code{what = "test"}, the standard error of the latent trait estimate
#'   is drawn as a dashed curve against a secondary axis. The layer is
#'   omitted when the standard error is not finite and varying over the
#'   grid (for example when test information is zero somewhere).
#' @param show_peak Logical. When \code{TRUE} (the default) and
#'   \code{what = "test"}, a vertical dotted line marks the value of
#'   \code{theta} at which test information peaks on the supplied grid.
#' @param palette Character string naming the color palette. Defaults to
#'   \code{"okabe_ito"}, base R's colorblind-safe Okabe-Ito palette;
#'   \code{"tableau"} is also available.
#' @param title Optional plot title.
#' @param xlab Label for the horizontal axis. Defaults to a description of
#'   the latent trait metric.
#' @param ylab Label for the vertical axis. Defaults to a description of
#'   the information plotted.
#'
#' @return A \code{ggplot2} object.
#'
#' @note Requires \pkg{ggplot2} (listed in \code{Suggests}).
#'
#' @details
#' The secondary axis is a linear rescaling of the primary axis, so the
#' dashed standard error curve shares the panel with the information curve
#' without either being distorted relative to its own axis. The standard
#' error is largest where information is smallest, which is why the two
#' curves run in opposite directions.
#'
#' @references
#' Embretson, S. E., & Reise, S. P. (2000). \emph{Item response theory for
#'   psychologists}. Lawrence Erlbaum.
#'
#' Samejima, F. (1969). Estimation of latent ability using a response
#'   pattern of graded scores. \emph{Psychometrika Monograph Supplement,
#'   34}(4, Pt. 2), 1--97.
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso \code{\link{irt_information}}
#'
#' @family plotting
#'
#' @keywords hplot
#'
#' @examples
#' info <- irt_information(
#'   a = c(mood_1 = 1.4, mood_2 = 0.9, mood_3 = 1.1),
#'   b = c(-1.5, -0.5, 0.5, 1.5, 0.0, 0.8),
#'   item = c(rep("mood_1", 4), "mood_2", "mood_3")
#' )
#'
#' # Test information with the standard error on the secondary axis.
#' plot_irt_information(info)
#'
#' # One curve per item.
#' plot_irt_information(info, what = "item")
#'
#' @export
plot_irt_information <- function(x, what = c("test", "item"),
                                 show_se = TRUE, show_peak = TRUE,
                                 palette = "okabe_ito",
                                 title = NULL, xlab = NULL, ylab = NULL) {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plot_irt_information(). ",
         "Install it with install.packages(\"ggplot2\").", call. = FALSE)
  }
  if (!is.data.frame(x) ||
      !all(c("theta", "test_information", "se") %in% names(x)) ||
      is.null(attr(x, "item_information"))) {
    stop("'x' must be the result of irt_information().", call. = FALSE)
  }
  what <- match.arg(what)
  for (flag in c("show_se", "show_peak")) {
    value <- get(flag)
    if (!is.logical(value) || length(value) != 1L || is.na(value)) {
      stop("'", flag, "' must be a single logical value (TRUE or FALSE).",
           call. = FALSE)
    }
  }

  if (is.null(xlab)) {
    xlab <- "Latent trait score (theta), in standard normal units"
  }

  if (what == "item") {
    item_information <- attr(x, "item_information")
    items <- colnames(item_information)
    plot_data <- data.frame(
      theta = rep(x$theta, times = ncol(item_information)),
      information = as.numeric(item_information),
      item = factor(rep(items, each = nrow(item_information)),
                    levels = items),
      stringsAsFactors = FALSE
    )
    if (is.null(ylab)) {
      ylab <- "Item information (the item's contribution to test information)"
    }
    plt <- ggplot2::ggplot(
      plot_data,
      ggplot2::aes(x = theta, y = information, color = item)
    ) +
      ggplot2::geom_line(linewidth = 0.8) +
      .dmar_discrete_scale(
        "colour",
        function(k) .dmar_palette(k, palette = palette),
        name = "Item"
      ) +
      ggplot2::labs(title = title, x = xlab, y = ylab) +
      ggplot2::theme_minimal(base_size = 12)
    return(plt)
  }

  # ---------- Test information ----------
  plot_data <- data.frame(
    theta = x$theta,
    information = x$test_information,
    se = x$se,
    stringsAsFactors = FALSE
  )
  if (is.null(ylab)) {
    ylab <- "Test information (the reciprocal of the squared standard error)"
  }
  line_colors <- .dmar_palette(2L, palette = palette)

  plt <- ggplot2::ggplot(
    plot_data, ggplot2::aes(x = theta, y = information)
  )

  if (isTRUE(show_peak)) {
    peak <- attr(x, "theta_max_information")
    if (is.numeric(peak) && length(peak) == 1L && is.finite(peak)) {
      plt <- plt + ggplot2::geom_vline(
        xintercept = peak, linetype = "dotted", color = "grey45",
        linewidth = 0.4
      )
    }
  }

  # The standard error shares the panel through a linear rescaling onto
  # the information axis, with the secondary axis inverting the map.
  information_max <- max(plot_data$information, na.rm = TRUE)
  finite_se <- plot_data[is.finite(plot_data$se), , drop = FALSE]
  se_low  <- if (nrow(finite_se)) min(finite_se$se) else NA_real_
  se_high <- if (nrow(finite_se)) max(finite_se$se) else NA_real_
  se_range <- se_high - se_low
  draw_se <- isTRUE(show_se) && nrow(finite_se) > 1L &&
    is.finite(se_range) && se_range > 0 && information_max > 0

  if (draw_se) {
    finite_se$se_rescaled <-
      (finite_se$se - se_low) / se_range * information_max
    plt <- plt +
      ggplot2::geom_line(
        data = finite_se,
        mapping = ggplot2::aes(x = theta, y = se_rescaled),
        color = line_colors[2L], linetype = "dashed", linewidth = 0.7
      ) +
      ggplot2::scale_y_continuous(
        sec.axis = ggplot2::sec_axis(
          ~ . / information_max * se_range + se_low,
          name = "Standard error of the estimated latent trait score (dashed)"
        )
      )
  }

  plt +
    ggplot2::geom_line(color = line_colors[1L], linewidth = 0.9) +
    ggplot2::labs(title = title, x = xlab, y = ylab) +
    ggplot2::theme_minimal(base_size = 12)
}
