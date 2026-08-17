# Suppress R CMD CHECK notes for ggplot2 aes() column references.
if (getRversion() >= "2.15.1") {
  utils::globalVariables(c("x", "density", "group"))
}

# Visualize a standardized mean difference with overlapping distributions.
#' Visualize a Standardized Mean Difference With Overlapping Distributions
#'
#' Creates a publication-quality plot showing two normal distributions separated
#' by the standardized mean difference (\emph{d}). The plot includes a
#' confidence interval for the population effect size and sample size
#' annotations, both shown by default.
#'
#' @param smd The standardized mean difference (Cohen's \emph{d}).
#' @param n_1 Sample size for Group 1.
#' @param n_2 Sample size for Group 2.
#' @param group_1 Raw data for Group 1. When provided, \code{smd}, \code{n_1},
#'   and \code{n_2} are computed from the data.
#' @param group_2 Raw data for Group 2.
#' @param conf_level Confidence level for the confidence interval (default
#'   \code{0.95}).
#' @param show_ci Logical. If \code{TRUE} (the default), a confidence interval
#'   for the population standardized mean difference is displayed beneath the
#'   distributions. Requires both \code{n_1} and \code{n_2}.
#' @param show_n Logical. If \code{TRUE} (the default), per-group sample sizes
#'   are annotated on the plot.
#' @param title Optional character string for the plot title. Defaults to
#'   \code{"Standardized Mean Difference"}.
#' @param group_labels Character vector of length 2 giving labels for the two
#'   groups. Defaults to \code{c("Group 1", "Group 2")}.
#' @param palette Character string naming the color palette used when
#'   \code{colors} is \code{NULL}. Defaults to \code{"okabe_ito"}, base R's
#'   colorblind-safe Okabe-Ito palette; \code{"tableau"} is also available.
#' @param colors Optional character vector of length 2 giving fill colors for
#'   the two groups. When \code{NULL} (the default), the first two colors of
#'   \code{palette} are used.
#'
#' @details
#' Two unit-variance normal distributions are drawn, centered at 0
#' (Group 2 / reference) and \emph{d} (Group 1 / focal). The semi-transparent
#' fills make the overlap visible, giving a direct visual impression of how much
#' the distributions differ.
#'
#' When \code{show_ci = TRUE} and both \code{n_1} and \code{n_2} are available,
#' the function calls \code{\link{ci_smd}} to compute the noncentral \emph{t}
#' based confidence interval and displays it as a horizontal bar beneath the
#' curves. A filled dot marks the point estimate and vertical caps mark the
#' confidence bounds.
#'
#' @return A \code{ggplot2} object that can be further customized with standard
#'   \pkg{ggplot2} layers, scales, and themes.
#'
#' @note Requires \pkg{ggplot2} (listed in \code{Suggests}). Install it with
#'   \code{install.packages("ggplot2")} if needed.
#'
#' @references
#' Cohen, J. (1988). \emph{Statistical power analysis for the behavioral
#'   sciences} (2nd ed.). Hillsdale, NJ: Lawrence Erlbaum.
#'
#' Kelley, K. (2007). Confidence intervals for standardized effect
#'   sizes: Theory, application, and implementation. \emph{Journal of Statistical
#'   Software, 20}(8), 1--24. \doi{10.18637/jss.v020.i08}
#'
#' @seealso \code{\link{smd}}, \code{\link{ci_smd}}, \code{\link{plot_ci}},
#'   \code{\link{plot_R2}}
#'
#' @examples
#' # From a known effect size and sample sizes.
#' plot_smd(smd = 0.50, n_1 = 50, n_2 = 50)
#'
#' # The variations below are not run, since the call above already shows
#' # the default display and each additional figure has to be drawn. From
#' # raw data, where the standardized mean difference and both sample
#' # sizes are taken from the data:
#' # set.seed(113)
#' # g1 <- rnorm(40, mean = 0.6, sd = 1)
#' # g2 <- rnorm(40, mean = 0.0, sd = 1)
#' # plot_smd(group_1 = g1, group_2 = g2)
#'
#' # Without the confidence interval or the sample size annotations:
#' # plot_smd(smd = 0.80, show_ci = FALSE, show_n = FALSE)
#'
#' # Custom group labels and title:
#' # plot_smd(smd = 0.45, n_1 = 75, n_2 = 75,
#' #          group_labels = c("Treatment", "Control"),
#' #          title = "Treatment Effect on Reading Scores")
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @keywords hplot
#'
#' @family plotting
#' @family confidence intervals for effect sizes
#'
#' @export
plot_smd <- function(smd = NULL, n_1 = NULL, n_2 = NULL,
                     group_1 = NULL, group_2 = NULL,
                     conf_level = 0.95,
                     show_ci = TRUE, show_n = TRUE,
                     title = NULL,
                     group_labels = c("Group 1", "Group 2"),
                     palette = "okabe_ito",
                     colors = NULL) {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plot_smd(). ",
         "Install it with install.packages(\"ggplot2\").", call. = FALSE)
  }

  if (is.null(colors)) colors <- .dmar_palette(2, palette = palette)

  # ---------- Compute d and sample sizes ----------
  if (!is.null(group_1) && !is.null(group_2)) {
    if (!is.null(smd)) {
      warning("Both raw data and 'smd' supplied; computing from the data.",
              call. = FALSE)
    }
    n_1 <- length(group_1)
    n_2 <- length(group_2)
    d_result <- DMAR::smd(group_1 = group_1, group_2 = group_2)
    d <- d_result$value[1]
  } else if (!is.null(smd)) {
    d <- smd
  } else {
    stop("Provide either raw data (group_1, group_2) or a standardized ",
         "mean difference (smd).", call. = FALSE)
  }

  # ---------- Compute CI ----------
  ci_lower <- ci_upper <- NULL
  if (show_ci) {
    if (!is.null(n_1) && !is.null(n_2)) {
      ci_result <- ci_smd(smd = d, n_1 = n_1, n_2 = n_2,
                          conf_level = conf_level)
      ci_lower <- ci_result$value[ci_result$term == "lower_limit"]
      ci_upper <- ci_result$value[ci_result$term == "upper_limit"]
    } else {
      warning("Cannot compute CI without n_1 and n_2; suppressing CI.",
              call. = FALSE)
      show_ci <- FALSE
    }
  }

  # ---------- Distribution data ----------
  x_lo <- min(-4, d - 4)
  x_hi <- max( 4, d + 4)
  x_seq <- seq(x_lo, x_hi, length.out = 500)

  curve_df <- data.frame(
    x       = rep(x_seq, 2),
    density = c(dnorm(x_seq, mean = 0, sd = 1),
                dnorm(x_seq, mean = d, sd = 1)),
    group   = factor(
      rep(c(group_labels[2], group_labels[1]), each = length(x_seq)),
      levels = group_labels
    )
  )

  # ---------- Layout constants ----------
  y_peak    <- dnorm(0)
  y_d_lab   <- -0.030
  y_ci      <- -0.070
  y_ci_lab  <- -0.105
  y_bottom  <- if (show_ci) -0.135 else -0.060

  # ---------- Core plot ----------
  p <- ggplot2::ggplot(curve_df,
         ggplot2::aes(x = x, y = density, fill = group)) +
    ggplot2::geom_area(alpha = 0.40, position = "identity", color = NA) +
    ggplot2::geom_line(
      ggplot2::aes(color = group),
      linewidth = 0.6, show.legend = FALSE
    ) +
    ggplot2::scale_fill_manual(
      values = stats::setNames(colors, group_labels), name = NULL
    ) +
    ggplot2::scale_color_manual(
      values = stats::setNames(colors, group_labels), guide = "none"
    )

  # Dashed mean lines
  p <- p +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed",
                        color = colors[2], alpha = 0.6, linewidth = 0.4) +
    ggplot2::geom_vline(xintercept = d, linetype = "dashed",
                        color = colors[1], alpha = 0.6, linewidth = 0.4)

  # ---------- d annotation ----------
  d_label <- paste0("italic(d) == ", format(round(d, 3), nsmall = 3))
  p <- p +
    ggplot2::annotate("text", x = d / 2, y = y_d_lab,
                      label = d_label, parse = TRUE,
                      size = 4, color = "grey20")

  # ---------- CI bar ----------
  if (show_ci) {
    ci_text <- paste0(round(conf_level * 100), "% CI [",
                      format(round(ci_lower, 3), nsmall = 3), ", ",
                      format(round(ci_upper, 3), nsmall = 3), "]")
    p <- p +
      # Horizontal line
      ggplot2::annotate("segment",
        x = ci_lower, xend = ci_upper, y = y_ci, yend = y_ci,
        linewidth = 0.7, color = "grey30"
      ) +
      # Point estimate dot
      ggplot2::annotate("point",
        x = d, y = y_ci, size = 2.5, shape = 16, color = "grey10"
      ) +
      # Left cap
      ggplot2::annotate("segment",
        x = ci_lower, xend = ci_lower,
        y = y_ci - 0.008, yend = y_ci + 0.008,
        linewidth = 0.5, color = "grey30"
      ) +
      # Right cap
      ggplot2::annotate("segment",
        x = ci_upper, xend = ci_upper,
        y = y_ci - 0.008, yend = y_ci + 0.008,
        linewidth = 0.5, color = "grey30"
      ) +
      # Label
      ggplot2::annotate("text",
        x = (ci_lower + ci_upper) / 2, y = y_ci_lab,
        label = ci_text, size = 3.3, color = "grey30"
      )
  }

  # ---------- Sample-size annotation ----------
  if (show_n && !is.null(n_1) && !is.null(n_2)) {
    n_label <- paste0("italic(n)[1] == ", n_1, "~~~italic(n)[2] == ", n_2)
    p <- p +
      ggplot2::annotate("text",
        x = x_hi - 0.3, y = y_peak * 0.95,
        label = n_label, parse = TRUE, size = 3.5,
        hjust = 1, color = "grey40"
      )
  }

  # ---------- Theme ----------
  if (is.null(title)) title <- "Standardized Mean Difference"

  p <- p +
    ggplot2::labs(title = title, x = "Standardized Scale", y = NULL) +
    ggplot2::coord_cartesian(
      ylim = c(y_bottom, y_peak * 1.08), clip = "off"
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      legend.position   = "bottom",
      panel.grid.minor  = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank(),
      axis.text.y       = ggplot2::element_blank(),
      plot.title        = ggplot2::element_text(hjust = 0.5, face = "bold")
    )

  p
}
