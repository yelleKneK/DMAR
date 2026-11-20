# Suppress R CMD CHECK notes for ggplot2 aes() column references.
if (getRversion() >= "2.15.1") {
  utils::globalVariables(c("estimate", "name", "lower", "upper"))
}

# Forest-plot-style confidence interval display.
#' Forest-Plot-Style Confidence Interval Display
#'
#' Plots one or more effect size estimates with their confidence intervals.
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
#' @param conf_level Confidence level of the limits, used only for the
#'   title (default \code{0.95}). When \code{ci} records the level of its
#'   interval, as every DMAR table does, that level is used, and a
#'   \code{conf_level} that differs from it is refused, since the plot
#'   draws the limits it is given and cannot change their level.
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
#'     \code{term = "smd"}) in addition to \code{"ci_lower"} and
#'     \code{"ci_upper"}.}
#'   \item{\bold{Long term/value without estimate}}{Output from
#'     \code{\link{ci_R}} or \code{\link{ci_R2}}, which contains only
#'     \code{"ci_lower"} and \code{"ci_upper"}. Supply the point estimate
#'     via the \code{estimate} argument.}
#'   \item{\bold{Wide per-effect format}}{Output from
#'     \code{\link{ci_omega_squared}} or \code{\link{ci_eta_squared}},
#'     which has one row per effect with columns for each estimate and its
#'     limits (\code{omega_squared_partial},
#'     \code{omega_squared_partial_ci_lower},
#'     \code{omega_squared_partial_ci_upper}, and the same for the total)
#'     and \code{N}, or from \code{\link{ci_eta_squared_generalized}},
#'     with \code{ci_lower} and \code{ci_upper}. Each estimate that has
#'     limits is drawn (the total of an effect of a factorial design has
#'     none), and when a table has limits for two, each row is labeled
#'     with the effect and the estimate.}
#' }
#'
#' The default title names the level the limits were formed at. An
#' interval formed from unequal tail areas (\code{alpha_lower} and
#' \code{alpha_upper} of a \code{ci_*} function) is drawn with a subtitle
#' that gives the two areas, or, when one of them is 0, says which limit is
#' the one-sided confidence bound, as \code{\link{results_sentence}} does.
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
#' omega_result <- ci_omega_squared(fit, estimand = "partial")
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
  .dmar_check_conf_level(conf_level, null_ok = TRUE)

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    .dmar_refuse(c(
      "Package 'ggplot2' is required for plot_ci(). ",
      "Install it with install.packages(\"ggplot2\")."),
      "package_required", package = "ggplot2")
  }

  point_color <- .dmar_palette(1, palette = palette)

  # ---------- Parse DMAR CI data.frame ----------
  tails <- NULL
  if (is.data.frame(ci)) {
    # The plot draws the limits it is given, so its title must name the
    # level they were formed at: the level the table records is used when
    # conf_level is not supplied, and a supplied level that differs is
    # refused (with results_sentence()'s tolerance), where the argument's
    # default of .95 titled a 97.5% or 90% table "95% Confidence
    # Interval". Unequal tail areas the table records are named in a
    # subtitle, as results_sentence() names them.
    recorded <- .plot_ci_recorded_level(ci)
    if (!is.null(recorded)) {
      if (missing(conf_level) || is.null(conf_level)) {
        conf_level <- recorded
      } else if (!is.numeric(conf_level) || length(conf_level) != 1L ||
                 !is.finite(conf_level) ||
                 abs(conf_level - recorded) > 1e-8) {
        .dmar_refuse_argument(c(
          "conf_level = ", format(conf_level), " is not the level the ",
          "'ci' table records for its limits, ",
          .format_conf_pct(recorded), "%. The plot draws the limits it ",
          "is given, so leave conf_level out to title it with their ",
          "level, or form the interval again at the level wanted."), arg = "ci")
      }
    }
    tails <- .recorded_tail_areas(ci)
    parsed <- .parse_dmar_ci(ci)
    if (is.null(estimate)) estimate <- parsed$estimate
    if (is.null(lower))    lower    <- parsed$lower
    if (is.null(upper))    upper    <- parsed$upper
    if (is.null(names))    names    <- parsed$names
    if (is.null(n))        n        <- parsed$n
  }

  # ---------- Validate ----------
  if (is.null(lower) || is.null(upper)) {
    .dmar_refuse_argument(c(
      "Confidence limits are required. Supply 'lower' and 'upper', ",
      "or a 'ci' data.frame from a DMAR ci_* function."),
      arg = c("lower", "upper", "ci"))
  }
  k <- length(lower)
  if (!is.null(estimate) && length(estimate) != k) {
    .dmar_refuse_argument(
      "'estimate' must be the same length as 'lower' and 'upper'.",
      arg = c("estimate", "lower", "upper"))
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
  # The level is written as it is (97.5, not 98; 99.9, not 100) by the
  # formatter of every displayed level, .format_conf_pct().
  ci_label <- paste0(.format_conf_pct(conf_level), "% Confidence Interval")
  if (is.null(title)) title <- ci_label
  # Unequal tail areas in the words of results_sentence(), short enough
  # for a subtitle: the two areas, or, when one is 0, which limit is the
  # one-sided confidence bound at the level; no subtitle for an
  # equal-tailed interval.
  subtitle <- if (is.null(tails)) NULL else {
    lo <- unname(tails[1L])
    hi <- unname(tails[2L])
    pct <- .format_conf_pct(conf_level)
    if (lo == 0) {
      sprintf("One-sided: the upper limit is the %s%% upper confidence bound",
              pct)
    } else if (hi == 0) {
      sprintf("One-sided: the lower limit is the %s%% lower confidence bound",
              pct)
    } else {
      sprintf("Unequal tail areas: %s below, %s above", .format_tail(lo),
              .format_tail(hi))
    }
  }

  p <- p +
    ggplot2::labs(
      title    = title,
      subtitle = subtitle,
      x        = xlab,
      y        = NULL
    ) +
    # Both annotations are centered on the estimate, so an estimate near
    # either end of the range would otherwise push its label into the panel
    # edge. The extra expansion buys the label room.
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = 0.12)) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      panel.grid.minor   = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank(),
      plot.title         = ggplot2::element_text(hjust = 0.5, face = "bold"),
      plot.subtitle      = ggplot2::element_text(hjust = 0.5)
    ) +
    ggplot2::coord_cartesian(clip = "off")

  p
}


# The level a table passed to plot_ci() records for its limits: its
# conf_level attribute, else a conf_level column holding one value (as the
# rows of a tidy() result do), each taken only when it is a finite number;
# NULL when it records none. A conf_level column holding several levels
# (a stack of intervals formed at different levels) is refused, since one
# title cannot name them. Not exported.
.plot_ci_recorded_level <- function(ci) {
  a <- attr(ci, "conf_level", exact = TRUE)
  if (is.numeric(a) && length(a) == 1L && is.finite(a)) return(a)
  if ("conf_level" %in% names(ci) && is.numeric(ci$conf_level)) {
    lv <- unique(ci$conf_level[is.finite(ci$conf_level)])
    if (length(lv) == 1L) return(lv)
    if (length(lv) > 1L)
      .dmar_refuse_argument(c(
        "The rows of 'ci' record different confidence levels (",
        paste0(.format_conf_pct(lv), "%", collapse = ", "), "); plot ",
        "the intervals of each level separately."), arg = "ci")
  }
  NULL
}

# Internal: parse the output of a DMAR ci_* function into a standard
# list with components estimate, lower, upper, names, and n.
.parse_dmar_ci <- function(ci) {
  # --- Wide per-effect format with limits named for each estimate
  # (ci_eta_squared() and ci_omega_squared(): eta_squared_partial,
  # eta_squared_partial_ci_lower, ...). Every estimate with limits is
  # drawn, one row per effect; the total of a factorial effect, whose
  # limits are NA, is left out. With two estimates the rows are labeled
  # "effect: estimate". ---
  if (!("ci_lower" %in% colnames(ci)) && !("term" %in% colnames(ci))) {
    nms <- colnames(ci)
    stems <- sub("_ci_lower$", "", nms[grepl("_ci_lower$", nms)])
    stems <- stems[stems %in% nms & paste0(stems, "_ci_upper") %in% nms]
    stems <- stems[vapply(stems, function(e)
      any(!is.na(ci[[paste0(e, "_ci_lower")]])), logical(1L))]
    if (length(stems)) {
      effect <- if ("effect" %in% nms) as.character(ci$effect) else
                  paste("Effect", seq_len(nrow(ci)))
      n_val <- if ("N" %in% nms) ci$N else NULL
      return(list(
        estimate = unlist(lapply(stems, function(e) ci[[e]])),
        lower    = unlist(lapply(stems, function(e)
                     ci[[paste0(e, "_ci_lower")]])),
        upper    = unlist(lapply(stems, function(e)
                     ci[[paste0(e, "_ci_upper")]])),
        names    = if (length(stems) == 1L) effect else
                     unlist(lapply(stems, function(e)
                       paste0(effect, ": ", e))),
        n        = if (is.null(n_val)) NULL else rep(n_val, length(stems))
      ))
    }
  }

  # --- Wide per-effect format with ci_lower and ci_upper
  # (ci_eta_squared_generalized()) ---
  if ("ci_lower" %in% colnames(ci) && "ci_upper" %in% colnames(ci) &&
      !("term" %in% colnames(ci))) {
    # Find the point-estimate column by name.
    est_candidates <- c("eta_squared_generalized", "smd", "R2", "estimate",
                        "value")
    est_col <- NULL
    for (col in est_candidates) {
      if (col %in% colnames(ci)) { est_col <- col; break }
    }
    names_val <- if ("effect" %in% colnames(ci)) ci$effect else
                   paste("Effect", seq_len(nrow(ci)))
    n_val     <- if ("N" %in% colnames(ci)) ci$N else NULL

    return(list(
      estimate = if (!is.null(est_col)) ci[[est_col]] else NULL,
      lower    = ci$ci_lower,
      upper    = ci$ci_upper,
      names    = names_val,
      n        = n_val
    ))
  }

  # --- Long term/value format (ci_smd, ci_R, ci_R2) ---
  if ("term" %in% colnames(ci) && "value" %in% colnames(ci)) {
    lower_val  <- ci$value[ci$term == "ci_lower"]
    upper_val  <- ci$value[ci$term == "ci_upper"]
    est_mask   <- !ci$term %in% c("ci_lower", "ci_upper")
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

  .dmar_refuse_argument(c(
    "Unrecognized data.frame format. Use explicit 'estimate', 'lower', ",
    "and 'upper' arguments."), arg = c("estimate", "lower", "upper"))
}
