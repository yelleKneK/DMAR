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
#' noninferiority bound \eqn{-\delta_L}, colored by the five-way
#' verdict of \code{\link{equivalence_c}}: an interval entirely inside the
#' region is equivalent; entirely above \eqn{\delta_U}, superior;
#' entirely below \eqn{-\delta_L}, inferior; a lower limit above
#' \eqn{-\delta_L} with an upper limit past \eqn{\delta_U},
#' noninferior only; and an interval straddling a bound,
#' inconclusive. The geometry \emph{is} the decision rule, which is
#' what makes the plot the natural report of an equivalence analysis.
#'
#' @param x Either a single result from \code{\link{equivalence_c}} or a list
#'   of them (a named list supplies the row labels). Alternatively,
#'   supply \code{estimate}, \code{lower}, and \code{upper} directly.
#' @param estimate,lower,upper Numeric vectors of contrast estimates
#'   and their confidence limits, used when \code{x} is not supplied.
#' @param names Optional character vector of row labels.
#' @param delta_lower,delta_upper Equivalence bounds, as positive
#'   magnitudes (the region drawn is \eqn{(-\delta_L, +\delta_U)}).
#'   Taken from \code{x} when it carries \code{equivalence_c} results;
#'   required otherwise. If only \code{delta_upper} is supplied, the
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
#' recomputed from the supplied limits and bounds, so the plot cannot
#' disagree with \code{\link{equivalence_c}}.
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
#' # Five constructed intervals, one per verdict, against bounds of 5
#' # (A equivalent, B noninferior only, C superior, D inconclusive,
#' #  E inferior):
#' plot_equivalence(estimate = c(-1.0, 3.5, 7.0, -1.5, -8.0),
#'                  lower    = c(-3.2, -1.4, 5.5, -6.6, -10.5),
#'                  upper    = c( 1.2,  8.4, 8.5,  3.6,  -5.5),
#'                  names    = c("A", "B", "C", "D", "E"),
#'                  delta_upper = 5)
#'
#' # From equivalence_c() results; a named list supplies the labels.
#' res <- list(
#'   "Focal vs. reference" = equivalence_c(psi_hat = -5.28, se = 2.49,
#'                                  df_error = 399, delta_upper = 5),
#'   "Within pipeline"     = equivalence_c(psi_hat = -0.53, se = 2.66,
#'                                  df_error = 399, delta_upper = 5)
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
                             delta_lower = NULL, delta_upper = NULL,
                             xlab = "Contrast", title = NULL,
                             palette = "okabe_ito") {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plot_equivalence(). ",
         "Install it with install.packages(\"ggplot2\").", call. = FALSE)
  }

  # ---------- Parse equivalence_c() input ----------
  is_equivalence_c <- function(obj) {
    is.data.frame(obj) && all(c("term", "value") %in% names(obj)) &&
      all(c("psi_hat", "lower_limit", "upper_limit",
            "delta_lower", "delta_upper") %in% obj$term)
  }
  pull <- function(obj, what) obj$value[obj$term == what]

  if (!is.null(x)) {
    if (is_equivalence_c(x)) x <- list(x)
    if (!is.list(x) || !all(vapply(x, is_equivalence_c, logical(1))))
      stop("'x' must be a equivalence_c() result or a list of them; otherwise ",
           "supply 'estimate', 'lower', and 'upper' directly.",
           call. = FALSE)
    estimate <- vapply(x, pull, numeric(1), what = "psi_hat")
    lower    <- vapply(x, pull, numeric(1), what = "lower_limit")
    upper    <- vapply(x, pull, numeric(1), what = "upper_limit")
    if (is.null(names))
      names <- if (!is.null(base::names(x))) base::names(x) else
        paste("Contrast", seq_along(x))
    dl <- vapply(x, pull, numeric(1), what = "delta_lower")
    du <- vapply(x, pull, numeric(1), what = "delta_upper")
    # A stored bound only matters when it is actually used, that is, when the
    # caller did not override it. Differing stored bounds block the plot only
    # for the bound whose explicit override is absent; a supplied delta_lower /
    # delta_upper resolves the ambiguity and the plot proceeds.
    if ((is.null(delta_lower) && length(unique(dl)) > 1L) ||
        (is.null(delta_upper) && length(unique(du)) > 1L))
      stop("The equivalence_c() results carry different equivalence bounds; ",
           "plot contrasts that share one region, or pass the bounds ",
           "explicitly.", call. = FALSE)
    if (is.null(delta_lower)) delta_lower <- -dl[1]   # stored signed
    if (is.null(delta_upper)) delta_upper <-  du[1]
  }

  # ---------- Validate ----------
  if (is.null(estimate) || is.null(lower) || is.null(upper))
    stop("Supply 'x' (equivalence_c results) or all of 'estimate', 'lower', ",
         "and 'upper'.", call. = FALSE)
  if (is.null(delta_upper))
    stop("'delta_upper' must be specified (the upper equivalence bound).",
         call. = FALSE)
  if (is.null(delta_lower)) delta_lower <- delta_upper
  if (delta_lower <= 0 || delta_upper <= 0)
    stop("The equivalence bounds must be positive magnitudes.",
         call. = FALSE)
  k <- length(estimate)
  if (length(lower) != k || length(upper) != k)
    stop("'estimate', 'lower', and 'upper' must have the same length.",
         call. = FALSE)
  if (is.null(names)) names <- paste("Contrast", seq_len(k))

  # ---------- Verdicts, recomputed so the plot matches equivalence_c() ----------
  verdict <- ifelse(lower > -delta_lower & upper < delta_upper, "Equivalent",
             ifelse(lower >  delta_upper, "Superior",
             ifelse(upper < -delta_lower, "Inferior",
             ifelse(lower > -delta_lower, "Noninferior only",
                    "Inconclusive"))))
  verdict_levels <- c("Equivalent", "Superior", "Noninferior only",
                      "Inconclusive", "Inferior")

  dat <- data.frame(
    label    = factor(names, levels = rev(names)),
    estimate = estimate, lower = lower, upper = upper,
    verdict  = factor(verdict, levels = verdict_levels)
  )

  colors <- stats::setNames(.dmar_palette(length(verdict_levels),
                                         palette = palette),
                            verdict_levels)

  ggplot2::ggplot(dat, ggplot2::aes(x = estimate, y = label,
                                    color = verdict)) +
    ggplot2::annotate("rect", xmin = -delta_lower, xmax = delta_upper,
                      ymin = -Inf, ymax = Inf, alpha = 0.10) +
    ggplot2::geom_vline(xintercept = 0, linewidth = 0.4) +
    ggplot2::geom_vline(xintercept = c(-delta_lower, delta_upper),
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
