# Picture of a regions-of-significance analysis: the estimated group
# difference as a function of the covariate, with its confidence band, the
# zero line it has to clear, and the covariate values where it clears it.

# Suppress R CMD check notes for ggplot2 aes() column references.
if (getRversion() >= "2.15.1") {
  utils::globalVariables(c("difference", "band_lower", "band_upper",
                           "pair", "covariate", "boundary"))
}

#' Plot Regions of Significance for a Covariate by Group Interaction
#'
#' Draws the estimated group difference \eqn{\hat D(x)} across the
#' observed range of the covariate, with the confidence band that the
#' region of significance is read from, a reference line at zero, and
#' vertical lines at the boundaries of the region. Wherever the band
#' clears zero the groups differ significantly, so the boundaries are
#' exactly the covariate values at which the band touches the zero line:
#' the plot \emph{is} the decision rule, which is what makes it the
#' natural report of the analysis.
#'
#' @param x A result of \code{\link{regions_of_significance}}.
#'   Alternatively a fitted \code{lm} or \code{aov}, or a formula with
#'   \code{data} supplied, in which case
#'   \code{\link{regions_of_significance}} is called first with
#'   \code{conf_level} and \code{method}.
#' @param data,conf_level,method Passed to
#'   \code{\link{regions_of_significance}}, and used only when \code{x}
#'   is a model or a formula rather than an already computed result.
#' @param xlab,ylab Axis labels. The defaults name the covariate and the
#'   group difference.
#' @param title Optional plot title.
#' @param palette Character string naming the color palette. Defaults to
#'   \code{"okabe_ito"}, base R's colorblind-safe Okabe-Ito palette;
#'   \code{"tableau"} is also available.
#' @param facet Logical. Draw one panel per pair of groups. Defaults to
#'   \code{TRUE} when there is more than one pair, which keeps the
#'   panels from overplotting each other; set it to \code{FALSE} to lay
#'   the pairs over one another in a single panel.
#' @param n_points Number of covariate values at which the difference
#'   and its band are evaluated. Default 200.
#'
#' @return A \code{ggplot} object. Requires \pkg{ggplot2} to be
#'   installed.
#'
#' @details
#' The band is \eqn{\hat D(x) \pm t_{crit} \sqrt{\mathrm{Var}[\hat
#' D(x)]}} with the same critical value used to find the boundaries, so
#' the picture and the table can never disagree. With the default
#' simultaneous critical value (Potthoff, 1964) the band is a
#' simultaneous band: it holds over the whole covariate range at once,
#' which is what licenses scanning it for the covariate values where the
#' groups differ.
#'
#' The band is drawn over the covariate values actually observed in the
#' two groups. A boundary that falls outside that range is therefore not
#' drawn, deliberately: it is an extrapolation of two fitted lines into
#' a region with no data, and drawing it would invite reading it as a
#' place where something was observed.
#'
#' @references
#' Johnson, P. O., & Neyman, J. (1936). Tests of certain linear
#'   hypotheses and their application to some educational problems.
#'   \emph{Statistical Research Memoirs, 1}, 57--93.
#'
#' Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027).
#'   \emph{Designing experiments and analyzing data: A model comparison
#'   perspective} (4th ed.). Routledge. (See Chapter 9 and its extension
#'   on heterogeneity of regression.)
#'
#' Potthoff, R. F. (1964). On the Johnson-Neyman technique and some
#'   extensions thereof. \emph{Psychometrika, 29}(3), 241--256.
#'   \doi{10.1007/BF02289721}
#'
#' @seealso \code{\link{regions_of_significance}},
#'   \code{\link{plot_ci}}
#'
#' @examples
#' # The Pygmalion teacher-expectancy data: post-test IQ (averaged over
#' # the two follow-ups) on pretest IQ, by condition. The expectancy
#' # effect is significant only in a band of pretest IQ values.
#' data(pygmalion)
#' pygmalion$iq_post <- (pygmalion$iq_4 + pygmalion$iq_8) / 2
#' fit <- lm(iq_post ~ iq_pre * treatment, data = pygmalion)
#' plot_regions_of_significance(fit)
#'
#' # Three groups: one panel per pair.
#' set.seed(113)
#' n <- 150
#' g <- factor(rep(c("control", "low", "high"), each = n / 3))
#' x <- rnorm(n, 50, 10)
#' y <- 2 + 0.5 * x + (g == "high") * (0.4 * x - 15) + rnorm(n, 0, 5)
#' plot_regions_of_significance(y ~ x * g, data = data.frame(y, x, g))
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @keywords hplot
#'
#' @family plotting
#'
#' @export

plot_regions_of_significance <- function(x, data = NULL, conf_level = 0.95,
                                         method = c("simultaneous",
                                                    "pointwise"),
                                         xlab = NULL, ylab = NULL,
                                         title = NULL,
                                         palette = "okabe_ito",
                                         facet = NULL, n_points = 200L) {

  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("Package 'ggplot2' is required for plot_regions_of_significance(). ",
         "Install it with install.packages(\"ggplot2\").", call. = FALSE)

  method <- match.arg(method)

  if (!inherits(x, "dmar_regions_of_significance"))
    x <- regions_of_significance(x, data = data, conf_level = conf_level,
                                 method = method)

  geometry  <- attr(x, "geometry")
  pair_info <- attr(x, "pairs")
  cov_name  <- attr(x, "covariate")
  if (is.null(geometry))
    stop("'x' does not carry the geometry of a regions-of-significance ",
         "analysis. Pass a result of regions_of_significance(), a fitted ",
         "model, or a formula with 'data'.", call. = FALSE)

  if (!is.numeric(n_points) || length(n_points) != 1L || is.na(n_points) ||
      n_points < 10)
    stop("'n_points' must be a single number of at least 10.", call. = FALSE)
  n_points <- as.integer(n_points)

  pair_levels <- geometry$pair
  if (is.null(facet)) facet <- length(pair_levels) > 1L
  if (!is.logical(facet) || length(facet) != 1L || is.na(facet))
    stop("'facet' must be TRUE or FALSE.", call. = FALSE)

  # ---------- The difference line and its band, pair by pair ----------
  curves <- lapply(seq_len(nrow(geometry)), function(k) {
    g <- geometry[k, ]
    xs <- seq(g$covariate_min, g$covariate_max, length.out = n_points)
    fitted_difference <- g$difference_intercept + g$difference_slope * xs
    half_width <- g$critical_value *
      sqrt(g$var_intercept + 2 * xs * g$cov_intercept_slope +
             xs^2 * g$var_slope)
    data.frame(pair = g$pair, covariate = xs,
               difference = fitted_difference,
               band_lower = fitted_difference - half_width,
               band_upper = fitted_difference + half_width,
               stringsAsFactors = FALSE)
  })
  curve_data <- do.call(rbind, curves)
  curve_data$pair <- factor(curve_data$pair, levels = pair_levels)

  # Only boundaries that fall inside the covariate values observed for the
  # pair are drawn; the rest are extrapolations, reported in the table but
  # deliberately not given a line on the picture.
  boundary_rows <- do.call(rbind, lapply(seq_len(nrow(geometry)), function(k) {
    g <- geometry[k, ]
    b <- c(g$lower_bound, g$upper_bound)
    b <- b[!is.na(b) & b >= g$covariate_min & b <= g$covariate_max]
    if (!length(b)) return(NULL)
    data.frame(pair = g$pair, boundary = b, stringsAsFactors = FALSE)
  }))
  if (!is.null(boundary_rows))
    boundary_rows$pair <- factor(boundary_rows$pair, levels = pair_levels)

  colors <- stats::setNames(.dmar_palette(length(pair_levels),
                                         palette = palette),
                            pair_levels)

  if (is.null(xlab)) xlab <- cov_name
  if (is.null(ylab)) {
    ylab <- if (nrow(pair_info) == 1L)
      paste0("Estimated difference (", pair_info$group_1, " - ",
             pair_info$group_2, ")")
    else "Estimated group difference"
  }

  p <- ggplot2::ggplot(
      curve_data,
      ggplot2::aes(x = covariate, y = difference, color = pair, fill = pair)
    ) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = band_lower, ymax = band_upper),
      alpha = 0.20, color = NA
    ) +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.5) +
    ggplot2::geom_line(linewidth = 0.8)

  if (!is.null(boundary_rows))
    p <- p + ggplot2::geom_vline(
      data = boundary_rows,
      ggplot2::aes(xintercept = boundary),
      linetype = "dashed", linewidth = 0.4, color = "grey30",
      inherit.aes = FALSE
    )

  p <- p +
    ggplot2::scale_color_manual(values = colors, name = NULL) +
    ggplot2::scale_fill_manual(values = colors, name = NULL) +
    ggplot2::labs(x = xlab, y = ylab, title = title) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank())

  if (facet) {
    p <- p + ggplot2::facet_wrap(~ pair) +
      ggplot2::theme(legend.position = "none")
  } else if (length(pair_levels) == 1L) {
    p <- p + ggplot2::theme(legend.position = "none")
  } else {
    p <- p + ggplot2::theme(legend.position = "bottom")
  }

  p
}
