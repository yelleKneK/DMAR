# Suppress R CMD CHECK notes for ggplot2 aes() column references.
if (getRversion() >= "2.15.1") {
  utils::globalVariables(c("estimate", "name", "lower", "upper"))
}

# Forest-plot-style confidence interval display.
#' Forest-Plot-Style Confidence Interval Display
#'
#' Creates a clean visualization of one or more effect size estimates with their
#' confidence intervals.
#'
#' The function accepts either (a) a \code{data.frame} produced by an
#' DMAR \code{ci_*} function (e.g., \code{\link{ci_smd}}, \code{\link{ci_R2}},
#' \code{\link{ci_omega_squared}}), or (b) explicit numeric vectors for the
#' estimate(s), lower bound(s), and upper bound(s).
#'
#' @param ci A \code{data.frame} from a DMAR \code{ci_*} function. When
#'   supplied, the function auto-detects the format and extracts the point
#'   estimate(s), lower limit(s), and upper limit(s). Explicit \code{estimate},
#'   \code{lower}, and \code{upper} arguments override values parsed from
#'   \code{ci}.
#' @param estimate Numeric vector of point estimates.
#' @param lower Numeric vector of lower confidence limits.
#' @param upper Numeric vector of upper confidence limits.
#' @param names Optional character vector of labels for each effect.
#' @param n Optional numeric vector (or scalar) of sample sizes. Recycled to
#'   match the number of effects.
#' @param conf_level Confidence level; used only for the axis label
#'   (default \code{0.95}).
#' @param show_n Logical. If \code{TRUE} (the default), the sample size is
#'   annotated above each estimate, with the estimate and its interval
#'   printed below.
#' @param reference_line Optional numeric value at which to draw a vertical
#'   reference line (e.g., \code{0} for mean differences, \code{1} for ratios).
#' @param xlab Label for the horizontal (effect size) axis. Defaults to
#'   \code{"Effect Size"}.
#' @param title Optional plot title.
#' @param palette Character string naming the color palette. The point
#'   estimates and interval bars are drawn in the palette's primary color.
#'   Defaults to \code{"okabe_ito"}, base R's colorblind-safe Okabe-Ito
#'   palette; \code{"tableau"} is also available.
#'
#' @details
#' The function recognizes three DMAR output formats:
#' \describe{
#'   \item{\bold{Long term/value with estimate row}}{Output from
#'     \code{\link{ci_smd}}, which includes a row for the point estimate (e.g.,
#'     \code{term = "smd"}) in addition to \code{"lower_limit"} and
#'     \code{"upper_limit"}.}
#'   \item{\bold{Long term/value without estimate}}{Output from
#'     \code{\link{ci_R}} or \code{\link{ci_R2}}, which contains only
#'     \code{"lower_limit"} and \code{"upper_limit"}. Supply the point estimate
#'     via the \code{estimate} argument.}
#'   \item{\bold{Wide per-effect format}}{Output from
#'     \code{\link{ci_omega_squared}}, which has one row per effect with columns
#'     for the point estimate, \code{lower_limit}, \code{upper_limit}, and
#'     \code{N}.}
#' }
#'
#' @return A \code{ggplot2} object.
#'
#' @note Requires \pkg{ggplot2} (listed in \code{Suggests}).
#'
#' @seealso \code{\link{ci_smd}}, \code{\link{ci_R}}, \code{\link{ci_R2}},
#'   \code{\link{ci_omega_squared}}, \code{\link{plot_smd}},
#'   \code{\link{plot_R2}}
#'
#' @examples
#' # From explicit values.
#' plot_ci(estimate = 0.45, lower = 0.15, upper = 0.75,
#'         names = "Cohen's d", n = 60, reference_line = 0)
#'
#' # From ci_smd() output.
#' ci_result <- ci_smd(smd = 0.5, n_1 = 50, n_2 = 50)
#' plot_ci(ci_result, n = 100, reference_line = 0)
#'
#' # Multiple effects from ci_omega_squared(): the expectancy treatment
#' # and the grade classification in the pygmalion data.
#' pyg <- pygmalion
#' pyg$grade <- factor(pyg$grade)
#' fit <- aov(iq_8 ~ treatment + grade, data = pyg)
#' omega_result <- ci_omega_squared(fit)
#' plot_ci(omega_result, reference_line = 0,
#'         xlab = expression(omega^2))
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @keywords hplot
#'
#' @family plotting
#'
#' @export
plot_ci <- function(ci = NULL, estimate = NULL, lower = NULL, upper = NULL,
                    names = NULL, n = NULL,
                    conf_level = 0.95,
                    show_n = TRUE,
                    reference_line = NULL,
                    xlab = "Effect Size",
                    title = NULL,
                    palette = "okabe_ito") {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plot_ci(). ",
         "Install it with install.packages(\"ggplot2\").", call. = FALSE)
  }

  point_color <- .dmar_palette(1, palette = palette)

  # ---------- Parse DMAR CI data.frame ----------
  if (is.data.frame(ci)) {
    parsed <- .parse_dmar_ci(ci)
    if (is.null(estimate)) estimate <- parsed$estimate
    if (is.null(lower))    lower    <- parsed$lower
    if (is.null(upper))    upper    <- parsed$upper
    if (is.null(names))    names    <- parsed$names
    if (is.null(n))        n        <- parsed$n
  }

  # ---------- Validate ----------
  if (is.null(lower) || is.null(upper)) {
    stop("Confidence limits are required. Supply 'lower' and 'upper', ",
         "or a 'ci' data.frame from a DMAR ci_* function.", call. = FALSE)
  }
  k <- length(lower)
  if (!is.null(estimate) && length(estimate) != k) {
    stop("'estimate' must be the same length as 'lower' and 'upper'.",
         call. = FALSE)
  }
  if (is.null(estimate)) {
    estimate <- (lower + upper) / 2
    message("No point estimate found; plotting the midpoint of each CI.")
  }

  if (is.null(names)) names <- if (k == 1) "Effect" else paste("Effect", seq_len(k))

  # ---------- Build plot data ----------
  plot_df <- data.frame(
    name     = factor(names, levels = rev(names)),
    estimate = estimate,
    lower    = lower,
    upper    = upper,
    stringsAsFactors = FALSE
  )

  # Recycle n to length k
  if (!is.null(n)) {
    n <- rep_len(n, k)
    plot_df$n <- n
  }

  # ---------- Build plot ----------
  p <- ggplot2::ggplot(plot_df,
         ggplot2::aes(x = estimate, y = name)) +
    ggplot2::geom_pointrange(
      ggplot2::aes(xmin = lower, xmax = upper),
      size = 0.7, linewidth = 0.6, color = point_color
    )

  # Reference line
  if (!is.null(reference_line)) {
    p <- p +
      ggplot2::geom_vline(
        xintercept = reference_line,
        linetype = "dashed", color = "grey50", linewidth = 0.4
      )
  }

  # Sample-size annotation, centered above each interval. Placing it beyond
  # the upper limit put it outside the panel whenever an interval ran wide,
  # where it was clipped; above the interval the horizontal extent of the
  # interval cannot reach it.
  if (show_n && !is.null(n)) {
    p <- p +
      ggplot2::annotate("text",
        x    = plot_df$estimate,
        y    = as.numeric(plot_df$name) + 0.28,
        label = paste0("italic(n) == ", plot_df$n),
        parse = TRUE, size = 3.4, color = "grey40"
      )
  }

  # Value labels beneath each bar
  p <- p +
    ggplot2::annotate("text",
      x     = plot_df$estimate,
      y     = as.numeric(plot_df$name) - 0.28,
      label = paste0(
        format(round(plot_df$estimate, 3), nsmall = 3),
        "  [",
        format(round(plot_df$lower, 3), nsmall = 3), ", ",
        format(round(plot_df$upper, 3), nsmall = 3), "]"
      ),
      size  = 3.2, color = "grey30"
    )

  # ---------- Theme ----------
  ci_label <- paste0(round(conf_level * 100), "% Confidence Interval")
  if (is.null(title)) title <- ci_label

  p <- p +
    ggplot2::labs(
      title = title,
      x     = xlab,
      y     = NULL
    ) +
    # Both annotations are centered on the estimate, so an estimate near
    # either end of the range would otherwise push its label into the panel
    # edge. The extra expansion buys the label room.
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = 0.12)) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      panel.grid.minor   = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank(),
      plot.title         = ggplot2::element_text(hjust = 0.5, face = "bold")
    ) +
    ggplot2::coord_cartesian(clip = "off")

  p
}


# Internal: parse the output of a DMAR ci_* function into a standard
# list with components estimate, lower, upper, names, and n.
.parse_dmar_ci <- function(ci) {
  # --- Wide per-effect format (ci_omega_squared) ---
  if ("lower_limit" %in% colnames(ci) && "upper_limit" %in% colnames(ci) &&
      !("term" %in% colnames(ci))) {
    # Find the point-estimate column by name.
    est_candidates <- c("omega_squared", "smd", "R2", "estimate", "value")
    est_col <- NULL
    for (col in est_candidates) {
      if (col %in% colnames(ci)) { est_col <- col; break }
    }
    names_val <- if ("effect" %in% colnames(ci)) ci$effect else
                   paste("Effect", seq_len(nrow(ci)))
    n_val     <- if ("N" %in% colnames(ci)) ci$N else NULL

    return(list(
      estimate = if (!is.null(est_col)) ci[[est_col]] else NULL,
      lower    = ci$lower_limit,
      upper    = ci$upper_limit,
      names    = names_val,
      n        = n_val
    ))
  }

  # --- Long term/value format (ci_smd, ci_R, ci_R2) ---
  if ("term" %in% colnames(ci) && "value" %in% colnames(ci)) {
    lower_val  <- ci$value[ci$term == "lower_limit"]
    upper_val  <- ci$value[ci$term == "upper_limit"]
    est_mask   <- !ci$term %in% c("lower_limit", "upper_limit")
    est_val    <- if (any(est_mask)) ci$value[est_mask][1]  else NULL
    est_name   <- if (any(est_mask)) ci$term[est_mask][1]   else NULL

    return(list(
      estimate = est_val,
      lower    = lower_val,
      upper    = upper_val,
      names    = est_name,
      n        = NULL
    ))
  }

  stop("Unrecognized data.frame format. Use explicit 'estimate', 'lower', ",
       "and 'upper' arguments.", call. = FALSE)
}
