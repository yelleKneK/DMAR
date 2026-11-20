# Suppress R CMD CHECK notes for ggplot2 aes() column references.
if (getRversion() >= "2.15.1") {
  utils::globalVariables(c("label", "verdict"))
}

# Forest plot of contrasts against an equivalence region.
#' Plot Contrasts Against an Equivalence Region
#'
#' Draws a forest-style plot of one or more contrast estimates with
#' their 100(1 - 2\eqn{\alpha})\% confidence intervals against the
#' equivalence region \eqn{(-\delta_L, \delta_U)} and the
#' noninferiority bound \eqn{-\delta_L}, colored by the six-way
#' verdict of \code{\link{equivalence_c}}: an interval entirely inside the
#' region is equivalent; entirely above \eqn{\delta_U}, superior;
#' entirely below \eqn{-\delta_L}, inferior; a lower limit above
#' \eqn{-\delta_L} with an upper limit at or past \eqn{\delta_U},
#' noninferior only; an upper limit below \eqn{\delta_U} with a lower
#' limit at or below \eqn{-\delta_L}, not relevantly superior only; and
#' an interval containing both bounds, inconclusive. The geometry
#' \emph{is} the decision rule, so the verdict for each contrast can be
#' read from where its interval falls.
#'
#' @param x Either a single result from \code{\link{equivalence_c}} or a list
#'   of them (a named list supplies the row labels). Alternatively,
#'   supply \code{estimate}, \code{lower}, and \code{upper} directly.
#' @param estimate,lower,upper Numeric vectors of contrast estimates
#'   and their confidence limits, used when \code{x} is not supplied.
#' @param names Optional character vector of row labels.
#' @param effect_lower,effect_upper Equivalence bounds on the scale of
#'   the estimates, as positive magnitudes (the region drawn is
#'   \eqn{(-\delta_L, +\delta_U)}, with \code{effect_lower} as
#'   \eqn{\delta_L} and \code{effect_upper} as \eqn{\delta_U}).
#'   Taken from \code{x} when it carries \code{equivalence_c} results;
#'   required otherwise. If only \code{effect_upper} is supplied, the
#'   bounds are symmetric.
#' @param xlab The horizontal axis label. Default
#'   \code{"Contrast"}.
#' @param title Optional plot title.
#' @param palette Character string naming the color palette for the
#'   verdict colors. Defaults to \code{"okabe_ito"}, base R's
#'   colorblind-safe Okabe-Ito palette; \code{"tableau"} is also
#'   available.
#'
#' @return A \code{ggplot} object. Requires \pkg{ggplot2} to be
#'   installed.
#'
#' @details
#' The shaded band is the equivalence region and the dashed vertical
#' lines are its bounds; the solid line at zero marks exact equality,
#' which is the null value of ordinary significance testing and is
#' deliberately \emph{not} a decision boundary here. Verdicts are
#' recomputed from the supplied limits and bounds by the rule
#' \code{\link{equivalence_c}} applies, so the plot cannot disagree with
#' it. The legend sets each one-sided verdict beside the stronger one
#' that implies it, with inconclusive between them (equivalent;
#' superior, noninferior only; inconclusive; not relevantly superior
#' only, inferior).
#'
#' @references
#' Chattopadhyay, B., Bandyopadhyay, T., Kelley, K., & Padalunkal,
#'   J. J. (2025). A sequential approach for noninferiority or
#'   equivalence of a linear contrast under cost constraints.
#'   \emph{Psychological Methods, 30}(2), 425--439. \doi{10.1037/met0000570}
#'
#' Schuirmann, D. J. (1987). A comparison of the two one-sided tests
#'   procedure and the power approach for assessing the equivalence of
#'   average bioavailability. \emph{Journal of Pharmacokinetics and
#'   Biopharmaceutics, 15}(6), 657--680.
#'
#' @seealso \code{\link{equivalence_c}}, \code{\link{plot_ci}}
#'
#' @examples
#' # Six constructed intervals, one per verdict, against bounds of 5
#' # (A equivalent, B superior, C noninferior only, D inconclusive,
#' #  E not relevantly superior only, F inferior):
#' plot_equivalence(estimate = c(-1.0, 7.0, 3.5, 0.5, -1.5, -8.0),
#'                  lower    = c(-3.2, 5.5, -1.4, -6.0, -6.6, -10.5),
#'                  upper    = c( 1.2, 8.5,  8.4,  7.0,  3.6,  -5.5),
#'                  names    = c("A", "B", "C", "D", "E", "F"),
#'                  effect_upper = 5)
#'
#' # From equivalence_c() results; a named list supplies the labels.
#' res <- list(
#'   "Focal vs. reference" = equivalence_c(contrast = -5.28, se = 2.49,
#'                                  df_error = 399, contrast_upper = 5),
#'   "Within pipeline"     = equivalence_c(contrast = -0.53, se = 2.66,
#'                                  df_error = 399, contrast_upper = 5)
#' )
#' plot_equivalence(res)
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @keywords hplot
#'
#' @family equivalence testing
#'
#' @family plotting
#'
#' @export

plot_equivalence <- function(x = NULL, estimate = NULL, lower = NULL,
                             upper = NULL, names = NULL,
                             effect_lower = NULL, effect_upper = NULL,
                             xlab = "Contrast", title = NULL,
                             palette = "okabe_ito") {
  .dmar_check_numbers("plot_equivalence", "effect_upper")
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    .dmar_refuse(c(
      "Package 'ggplot2' is required for plot_equivalence(). ",
      "Install it with install.packages(\"ggplot2\")."),
      "package_required", package = "ggplot2")
  }

  # ---------- Parse equivalence_c() input ----------
  is_equivalence_c <- function(obj) {
    is.data.frame(obj) && all(c("term", "value") %in% names(obj)) &&
      all(c("contrast", "ci_lower", "ci_upper",
            "contrast_lower", "contrast_upper") %in% obj$term)
  }
  pull <- function(obj, what) obj$value[obj$term == what]

  if (!is.null(x)) {
    if (is_equivalence_c(x)) x <- list(x)
    if (!is.list(x) || !all(vapply(x, is_equivalence_c, logical(1))))
      .dmar_refuse_argument(c(
        "'x' must be an equivalence_c() result or a list of them; otherwise ",
        "supply 'estimate', 'lower', and 'upper' directly."),
        arg = c("x", "estimate", "lower"))
    estimate <- vapply(x, pull, numeric(1), what = "contrast")
    lower    <- vapply(x, pull, numeric(1), what = "ci_lower")
    upper    <- vapply(x, pull, numeric(1), what = "ci_upper")
    if (is.null(names))
      names <- if (!is.null(base::names(x))) base::names(x) else
        paste("Contrast", seq_along(x))
    dl <- vapply(x, pull, numeric(1), what = "contrast_lower")
    du <- vapply(x, pull, numeric(1), what = "contrast_upper")
    # A stored bound only matters when it is actually used, that is, when the
    # caller did not override it. Differing stored bounds block the plot only
    # for the bound whose explicit override is absent; a supplied effect_lower /
    # effect_upper resolves the ambiguity and the plot proceeds.
    if ((is.null(effect_lower) && length(unique(dl)) > 1L) ||
        (is.null(effect_upper) && length(unique(du)) > 1L))
      .dmar_refuse_argument(c(
        "The equivalence_c() results carry different equivalence bounds; ",
        "plot contrasts that share one region, or pass the bounds ",
        "explicitly."), arg = NULL)
    if (is.null(effect_lower)) effect_lower <- -dl[1]   # stored signed
    if (is.null(effect_upper)) effect_upper <-  du[1]
  }

  # ---------- Validate ----------
  if (is.null(estimate) || is.null(lower) || is.null(upper))
    .dmar_refuse_argument(c(
      "Supply 'x' (equivalence_c results) or all of 'estimate', 'lower', ",
      "and 'upper'."), arg = c("x", "estimate", "lower"))
  if (is.null(effect_upper))
    .dmar_refuse_argument(
      "'effect_upper' must be specified (the upper equivalence bound).",
      arg = "effect_upper")
  if (is.null(effect_lower)) effect_lower <- effect_upper
  if (effect_lower <= 0 || effect_upper <= 0)
    .dmar_refuse_argument(
      "The equivalence bounds must be positive magnitudes.", arg = NULL)
  k <- length(estimate)
  if (length(lower) != k || length(upper) != k)
    .dmar_refuse_argument(
      "'estimate', 'lower', and 'upper' must have the same length.",
      arg = c("estimate", "lower", "upper"))
  if (is.null(names)) names <- paste("Contrast", seq_len(k))

  # ---------- Verdicts, by the rule equivalence_c() applies ----------
  verdict <- .equivalence_verdict(lower, upper, effect_lower, effect_upper)

  # Legend order: each one-sided verdict beside the stronger verdict that
  # implies it, with Inconclusive between the two pairs.
  verdict_levels <- c("Equivalent", "Superior", "Noninferior only",
                      "Inconclusive", "Not relevantly superior only",
                      "Inferior")

  dat <- data.frame(
    label    = factor(names, levels = rev(names)),
    estimate = estimate, lower = lower, upper = upper,
    verdict  = factor(verdict, levels = verdict_levels)
  )

  # Colors are assigned by name, in the order of the DMAR 1.0.x legend with
  # the sixth verdict (added in 1.1.0) after it, so that a verdict keeps the
  # color it had and the new one takes the palette's next color.
  color_order <- c("Equivalent", "Superior", "Noninferior only",
                   "Inconclusive", "Inferior",
                   "Not relevantly superior only")
  colors <- stats::setNames(.dmar_palette(length(color_order),
                                         palette = palette),
                            color_order)

  ggplot2::ggplot(dat, ggplot2::aes(x = estimate, y = label,
                                    color = verdict)) +
    ggplot2::annotate("rect", xmin = -effect_lower, xmax = effect_upper,
                      ymin = -Inf, ymax = Inf, alpha = 0.10) +
    ggplot2::geom_vline(xintercept = 0, linewidth = 0.4) +
    ggplot2::geom_vline(xintercept = c(-effect_lower, effect_upper),
                        linetype = "dashed", linewidth = 0.4) +
    # size 0.375 draws the point at 0.375 * 4 = 1.5 (the geom multiplies
    # size by its default fatten of 4 on every supported ggplot2), matching
    # the pre-4.0 fatten = 3 appearance without the deprecated argument.
    ggplot2::geom_pointrange(ggplot2::aes(xmin = lower, xmax = upper),
                             linewidth = 0.8, size = 0.375) +
    ggplot2::scale_color_manual(values = colors, drop = FALSE,
                                name = "Verdict") +
    ggplot2::labs(x = xlab, y = NULL, title = title) +
    ggplot2::theme_minimal()
}
