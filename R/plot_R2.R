# Visualize the proportion of variance explained (R-squared).
#' Visualize the Proportion of Variance Explained (\eqn{R^2})
#'
#' Creates a horizontal bar chart showing the observed \eqn{R^2} as a proportion
#' of total variance, with an optional confidence interval displayed beneath the
#' bar and sample size / predictor-count annotations.
#'
#' @param R2 The observed squared multiple correlation coefficient
#'   (\eqn{0 \le R^2 \le 1}).
#' @param N Total sample size.
#' @param p Number of predictors.
#' @param conf_level Confidence level for the confidence interval
#'   (default \code{0.95}).
#' @param show_ci Logical. If \code{TRUE} (the default), a confidence interval
#'   is shown beneath the proportion bar. Requires both \code{N} and \code{p}.
#' @param show_n Logical. If \code{TRUE} (the default), \code{N} and \code{p}
#'   are annotated on the plot.
#' @param random_predictors Logical. Whether the predictors are random
#'   (\code{TRUE}, the default) or fixed. Passed to \code{\link{ci_R2}}.
#' @param title Optional plot title.
#' @param palette Character string naming the color palette used for the
#'   \dQuote{Explained} portion of the bar when \code{colors} is \code{NULL}.
#'   Defaults to \code{"okabe_ito"}, base R's colorblind-safe Okabe-Ito
#'   palette; \code{"tableau"} is also available.
#' @param colors Optional character vector of length 2: the first color fills
#'   the \dQuote{Explained} portion of the bar, the second the
#'   \dQuote{Unexplained} portion. When \code{NULL} (the default), the
#'   \dQuote{Explained} portion uses the first color of \code{palette} and
#'   the \dQuote{Unexplained} portion a neutral light gray.
#'
#' @return A \code{ggplot2} object.
#'
#' @note Requires \pkg{ggplot2} (listed in \code{Suggests}).
#'
#' @seealso \code{\link{ci_R2}}, \code{\link{ci_R}}, \code{\link{plot_ci}},
#'   \code{\link{plot_smd}}
#'
#' @examples
#' # Basic call.
#' plot_R2(R2 = 0.25, N = 100, p = 5)
#'
#' # The variations below are not run, since the call above already shows
#' # the default display and each additional figure has to be drawn. With
#' # fixed predictors and a 90% confidence interval:
#' # plot_R2(R2 = 0.35, N = 200, p = 3, conf_level = 0.90,
#' #         random_predictors = FALSE)
#'
#' # Without the annotations:
#' # plot_R2(R2 = 0.10, show_ci = FALSE, show_n = FALSE)
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @references
#' Kelley, K. (2008). Sample size planning for the squared multiple
#'   correlation coefficient: Accuracy in parameter estimation via narrow
#'   confidence intervals. \emph{Multivariate Behavioral Research, 43},
#'   524--555. \doi{10.1080/00273170802490632}
#'
#' Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). \emph{Designing
#'   experiments and analyzing data: A model comparison perspective}
#'   (4th ed.). Routledge. (See Chapter 3 on \eqn{R^2} as a model comparison
#'   effect size.)
#'
#' @keywords hplot
#'
#' @family plotting
#'
#' @export
plot_R2 <- function(R2, N = NULL, p = NULL,
                    conf_level = 0.95,
                    show_ci = TRUE, show_n = TRUE,
                    random_predictors = TRUE,
                    title = NULL,
                    palette = "okabe_ito",
                    colors = NULL) {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plot_R2(). ",
         "Install it with install.packages(\"ggplot2\").", call. = FALSE)
  }
  if (R2 < 0 || R2 > 1) stop("'R2' must be between 0 and 1.", call. = FALSE)

  # Explained portion takes the palette's primary; the unexplained
  # portion keeps a neutral light gray as a visual reference.
  if (is.null(colors)) colors <- c(.dmar_palette(1, palette = palette), "#D4D4D4")

  # ---------- Compute CI ----------
  ci_lower <- ci_upper <- NULL
  if (show_ci) {
    if (!is.null(N) && !is.null(p)) {
      ci_result <- ci_R2(R2 = R2, N = N, p = p, conf_level = conf_level,
                         random_predictors = random_predictors)
      ci_lower <- ci_result$value[ci_result$term == "lower_limit"]
      ci_upper <- ci_result$value[ci_result$term == "upper_limit"]
    } else {
      warning("Cannot compute CI without N and p; suppressing CI.",
              call. = FALSE)
      show_ci <- FALSE
    }
  }

  # ---------- Layout constants ----------
  bar_lo  <- 0.3
  bar_hi  <- 0.7
  y_ci    <- 0.10
  y_ci_lab <- -0.04
  y_bottom <- if (show_ci) -0.12 else 0.15

  # ---------- Build the plot ----------
  plt <- ggplot2::ggplot() +
    # Full bar (unexplained)
    ggplot2::annotate("rect",
      xmin = 0, xmax = 1, ymin = bar_lo, ymax = bar_hi,
      fill = colors[2], color = "grey50", linewidth = 0.3
    ) +
    # Explained portion
    ggplot2::annotate("rect",
      xmin = 0, xmax = R2, ymin = bar_lo, ymax = bar_hi,
      fill = colors[1], alpha = 0.85
    ) +
    # Full outline
    ggplot2::annotate("rect",
      xmin = 0, xmax = 1, ymin = bar_lo, ymax = bar_hi,
      fill = NA, color = "grey40", linewidth = 0.5
    )

  # ---------- Labels inside the bar ----------
  bar_mid <- (bar_lo + bar_hi) / 2

  # R-squared label: place inside the explained portion when wide enough,
  # otherwise to its right.
  r2_label <- paste0("italic(R)^2 == ",
                     format(round(R2, 3), nsmall = 3))
  if (R2 >= 0.18) {
    plt <- plt + ggplot2::annotate("text",
      x = R2 / 2, y = bar_mid, label = r2_label, parse = TRUE,
      size = 4.5, fontface = "bold", color = "white"
    )
  } else {
    plt <- plt + ggplot2::annotate("text",
      x = R2 + 0.02, y = bar_mid, label = r2_label, parse = TRUE,
      size = 4, fontface = "bold", color = "grey20", hjust = 0
    )
  }

  # Unexplained label (only if space permits)
  if (1 - R2 >= 0.30) {
    plt <- plt + ggplot2::annotate("text",
      x = R2 + (1 - R2) / 2, y = bar_mid,
      label = "Unexplained", size = 3.5, color = "grey50"
    )
  }

  # ---------- CI bar ----------
  if (show_ci) {
    ci_text <- paste0(round(conf_level * 100), "% CI [",
                      format(round(ci_lower, 3), nsmall = 3), ", ",
                      format(round(ci_upper, 3), nsmall = 3), "]")
    cap <- 0.03
    plt <- plt +
      ggplot2::annotate("segment",
        x = ci_lower, xend = ci_upper, y = y_ci, yend = y_ci,
        linewidth = 0.7, color = "grey30"
      ) +
      ggplot2::annotate("point",
        x = R2, y = y_ci, size = 2.5, shape = 16, color = "grey10"
      ) +
      ggplot2::annotate("segment",
        x = ci_lower, xend = ci_lower,
        y = y_ci - cap, yend = y_ci + cap,
        linewidth = 0.5, color = "grey30"
      ) +
      ggplot2::annotate("segment",
        x = ci_upper, xend = ci_upper,
        y = y_ci - cap, yend = y_ci + cap,
        linewidth = 0.5, color = "grey30"
      ) +
      ggplot2::annotate("text",
        x = (ci_lower + ci_upper) / 2, y = y_ci_lab,
        label = ci_text, size = 3.3, color = "grey30"
      )
  }

  # ---------- Sample-size / predictor annotation ----------
  if (show_n && !is.null(N) && !is.null(p)) {
    n_label <- paste0("italic(N) == ", N, "~~~italic(p) == ", p)
    plt <- plt + ggplot2::annotate("text",
      x = 0.98, y = bar_hi + 0.08, label = n_label,
      parse = TRUE, size = 3.5, hjust = 1, color = "grey40"
    )
  }

  # ---------- Scales and theme ----------
  if (is.null(title)) title <- expression(bold("Proportion of Variance Explained"))

  pct_breaks <- seq(0, 1, 0.25)
  pct_labels <- paste0(pct_breaks * 100, "%")

  plt <- plt +
    ggplot2::scale_x_continuous(
      limits = c(-0.02, 1.02),
      breaks = pct_breaks,
      labels = pct_labels,
      expand = c(0, 0)
    ) +
    ggplot2::coord_cartesian(
      ylim = c(y_bottom, bar_hi + 0.15), clip = "off"
    ) +
    ggplot2::labs(title = title, x = NULL, y = NULL) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      axis.text.y        = ggplot2::element_blank(),
      axis.ticks.y       = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor   = ggplot2::element_blank(),
      plot.title         = ggplot2::element_text(hjust = 0.5)
    )

  plt
}
