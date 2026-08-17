# Colors for the DMAR plotting functions.
#
# The plot_* family colors itself through the internal .dmar_palette()
# helper so every DMAR figure shares one accessible default: base R's
# Okabe-Ito colorblind-safe qualitative palette (grDevices::palette.colors).
# The base R Tableau 10 set is offered as an alternative. DMAR defines no
# color of its own and ships no institutional palette; the helper only
# forwards to base R, so it adds no dependency and no branding. A user who
# wants other colors adds an ordinary ggplot2 scale to the returned plot.
# Neither function is exported: color choice is a property of the plots,
# not a public utility DMAR needs to expose.

# Normalize a user-facing palette name to its base R palette.colors() key.
.dmar_resolve_palette_name <- function(palette) {
  if (!is.character(palette) || length(palette) != 1L || is.na(palette)) {
    stop("'palette' must be a single character string.", call. = FALSE)
  }
  switch(palette,
    "okabe-ito" = ,
    "okabe_ito" = "Okabe-Ito",
    "tableau"   = "Tableau 10",
    stop("Unknown 'palette': \"", palette, "\". Choose \"okabe_ito\" ",
         "or \"tableau\".", call. = FALSE)
  )
}

# n colors from a base R qualitative palette, interpolating with
# colorRampPalette only when n exceeds the anchors; the common path
# (n <= number of anchors) returns the anchors verbatim, unnamed. NULL
# returns the full anchor set. This is the single color source the
# plot_* functions draw from.
.dmar_palette <- function(n = NULL, palette = "okabe_ito", reverse = FALSE) {
  pal <- .dmar_resolve_palette_name(palette)

  if (!is.logical(reverse) || length(reverse) != 1L || is.na(reverse)) {
    stop("'reverse' must be a single logical value (TRUE or FALSE).",
         call. = FALSE)
  }
  if (!is.null(n)) {
    if (!is.numeric(n) || length(n) != 1L || is.na(n) ||
        n < 1 || n != round(n)) {
      stop("'n' must be a single positive integer, or NULL.", call. = FALSE)
    }
    n <- as.integer(n)
  }

  anchors <- unname(grDevices::palette.colors(palette = pal))
  out <- if (is.null(n)) {
    anchors
  } else if (n <= length(anchors)) {
    anchors[seq_len(n)]
  } else {
    grDevices::colorRampPalette(anchors)(n)
  }

  if (reverse) rev(out) else out
}

# Build a ggplot2 discrete scale from a palette function, across the version
# boundary at which the required `scale_name` positional argument was retired.
# ggplot2 (>= 3.5.0) defaults `scale_name` to deprecated() and warns if it is
# supplied; ggplot2 (< 3.5.0) requires it. Used by the plot_* functions that
# color a discrete grouping variable through .dmar_palette().
.dmar_discrete_scale <- function(aesthetics, palette_fn, ...) {
  if (utils::packageVersion("ggplot2") >= "3.5.0") {
    ggplot2::discrete_scale(aesthetics, palette = palette_fn, ...)
  } else {
    ggplot2::discrete_scale(aesthetics, scale_name = "dmar",
                            palette = palette_fn, ...)
  }
}
