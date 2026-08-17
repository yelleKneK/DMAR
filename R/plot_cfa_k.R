# Suppress R CMD check notes for ggplot2 aes() column references.
if (getRversion() >= "2.15.1") {
  utils::globalVariables(c("item", "ci_lower", "ci_upper", "ref"))
}

#' Plot the Estimates of a Multiple-Factor CFA
#'
#' Displays the item-level estimates of a \code{\link{cfa_k}} fit, one
#' panel per factor, with each estimate's confidence interval. The
#' display is built to make the equality questions behind the classical
#' measurement structures visible: a dashed vertical line marks, per
#' factor, either the common (equated) estimate when the plotted
#' parameter was constrained equal, or the mean of the free estimates as
#' an informal anchor for the question "could these plausibly be one
#' value?". Confidence intervals that all cover the anchor are what
#' equal loadings (or equal error variances, or equal intercepts) would
#' look like; an interval far from it shows which item resists the
#' constraint, and the likelihood ratio test of the two nested
#' \code{cfa_k()} fits is the formal companion (see the examples in
#' \code{\link{cfa_k}}).
#'
#' @param x A \code{dmar_cfa_k} object from \code{\link{cfa_k}} with the
#'   default \code{output = "verbose"}.
#' @param what Which parameter to display: \code{"loadings"} (default,
#'   the \code{lambda} terms), \code{"errors"} (the \code{psi} terms), or
#'   \code{"intercepts"} (the \code{nu} terms; requires a fit with the
#'   mean structure).
#' @param show_equal_reference Logical. If \code{TRUE} (default), draw
#'   the dashed per-factor reference line described above. When the
#'   parameter was constrained equal the line is the common estimate and
#'   is always drawn.
#' @param xlab Label for the horizontal axis. Defaults to a description
#'   of the plotted parameter.
#' @param title Optional plot title.
#' @param palette Character string naming the color palette. Defaults to
#'   \code{"okabe_ito"}, base R's colorblind-safe Okabe-Ito palette;
#'   \code{"tableau"} is also available.
#'
#' @return A \code{ggplot2} object.
#'
#' @note Requires \pkg{ggplot2} (listed in \code{Suggests}).
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso \code{\link{cfa_k}} for the fit; \code{\link{plot_ci}} for
#'   the general forest-style confidence interval display.
#'
#' @family plotting
#'
#' @keywords hplot
#'
#' @examples
#' data(holzinger_swineford)
#' hs_factors <- list(
#'   verbal = c("t6_paragraph_comprehension",
#'              "t7_sentence", "t9_word_meaning"),
#'   deduction = c("t20_deduction", "t22_problem_reasoning",
#'              "t23_series_completion"))
#' res <- cfa_k(holzinger_swineford, hs_factors)
#'
#' # Are equal loadings plausible? Compare each interval with the anchor.
#' plot_cfa_k(res)
#'
#' # Two further displays are shown but not run here, since each draws
#' # another figure and the second refits the model as well. The same
#' # question for the error variances, the additional constraint that
#' # separates essentially parallel from essentially tau-equivalent:
#' # plot_cfa_k(res, what = "errors")
#' #
#' # After imposing the constraint, every item in a factor sits at the
#' # common estimate and the dashed line is that estimate rather than
#' # the mean of the free ones:
#' # res_equal <- cfa_k(holzinger_swineford, hs_factors,
#' #                    equal_loading = TRUE)
#' # plot_cfa_k(res_equal)
#'
#' @export
plot_cfa_k <- function(x, what = c("loadings", "errors", "intercepts"),
                       show_equal_reference = TRUE,
                       xlab = NULL, title = NULL,
                       palette = "okabe_ito") {
  what <- match.arg(what)
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plot_cfa_k(). ",
         "Install it with install.packages(\"ggplot2\").", call. = FALSE)
  }
  factors <- attr(x, "factors")
  if (!inherits(x, "dmar_cfa_k") || is.null(factors)) {
    stop("'x' must be a cfa_k() result with the default ",
         "output = \"verbose\".", call. = FALSE)
  }
  prefix <- switch(what, loadings = "lambda", errors = "psi",
                   intercepts = "nu")
  if (what == "intercepts" && !isTRUE(attr(x, "has_means"))) {
    stop("The fit has no mean structure, so there are no intercepts to ",
         "plot; refit with meanstructure = TRUE (or equal_intercept).",
         call. = FALSE)
  }

  panels <- lapply(names(factors), function(f) {
    v <- factors[[f]]
    term_free <- paste0(prefix, "_", f, "_", seq_along(v))
    term_equal <- paste0(prefix, "_", f)
    rows <- x[x$term %in% c(term_free, term_equal), , drop = FALSE]
    if (nrow(rows) == 0L) return(NULL)
    item <- if (what == "loadings") {
      sub(paste0("^", f, " =~ "), "", rows$syntax)
    } else {
      sub(" ~.*$", "", rows$syntax)
    }
    data.frame(factor = f, item = item, term = rows$term,
               estimate = rows$estimate,
               ci_lower = rows$ci_lower, ci_upper = rows$ci_upper,
               equal = rows$term %in% term_equal,
               stringsAsFactors = FALSE)
  })
  df <- do.call(rbind, panels)
  if (is.null(df) || nrow(df) == 0L) {
    stop("No ", what, " terms found in 'x'.", call. = FALSE)
  }
  df$factor <- base::factor(df$factor, levels = names(factors))
  df$item <- base::factor(df$item, levels = rev(unique(df$item)))

  refs <- do.call(rbind, lapply(split(df, df$factor), function(d) {
    data.frame(factor = d$factor[1L],
               ref = if (all(d$equal)) d$estimate[1L] else mean(d$estimate),
               equal = all(d$equal), stringsAsFactors = FALSE)
  }))
  if (!show_equal_reference) refs <- refs[refs$equal, , drop = FALSE]

  k <- length(names(factors))
  colors <- stats::setNames(.dmar_palette(k, palette = palette),
                            names(factors))
  if (is.null(xlab)) {
    xlab <- switch(what,
                   loadings = "Factor loading (unstandardized)",
                   errors = "Error variance",
                   intercepts = "Intercept")
  }

  p <- ggplot2::ggplot(df, ggplot2::aes(x = estimate, y = item,
                                        color = factor))
  if (nrow(refs) > 0L) {
    p <- p + ggplot2::geom_vline(data = refs,
                                 ggplot2::aes(xintercept = ref),
                                 linetype = "dashed", color = "gray40")
  }
  if (all(is.na(df$ci_lower))) {
    p <- p + ggplot2::geom_point(size = 2)
  } else {
    p <- p + ggplot2::geom_pointrange(
      ggplot2::aes(xmin = ci_lower, xmax = ci_upper))
  }
  p +
    ggplot2::facet_wrap(~ factor, scales = "free_y") +
    ggplot2::scale_color_manual(values = colors, guide = "none") +
    ggplot2::labs(x = xlab, y = "Item", title = title)
}
