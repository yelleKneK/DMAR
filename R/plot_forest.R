if (getRversion() >= "2.15.1")
  utils::globalVariables(c("lower", "upper", "pooled", "label", "size"))

#' Forest Plot of Study Effect Sizes With the Pooled Estimate
#'
#' Draws the meta-analyst's central picture: every study's effect size with
#' its confidence interval, the random effects pooled estimate beneath them,
#' and, by default, the prediction interval showing where the effect of a
#' \emph{new} study is expected to land. Point sizes are proportional to
#' precision (inverse variance), so the eye weighs the studies the way the
#' model does. Requires \pkg{ggplot2}.
#'
#' @param yi Numeric vector of study effect sizes.
#' @param vi Sampling variances of \code{yi}.
#' @param labels Optional study labels, one per study; defaults to
#'   \code{"Study 1"}, \code{"Study 2"}, and so on, in the supplied order.
#' @param method,hartung_knapp Passed to \code{\link{meta_es}} for the
#'   pooled row (\code{"reml"} and \code{TRUE} by default).
#' @param conf_level Confidence level for the per-study and pooled
#'   intervals. Defaults to 0.95.
#' @param show_prediction Logical: draw the prediction interval band on the
#'   pooled row? Default \code{TRUE} (ignored for \code{method = "fe"},
#'   which has none).
#' @param xlab Label for the effect size axis. Defaults to
#'   \code{"Effect size"}.
#' @param title Optional plot title.
#' @param palette Palette name. Defaults to \code{"okabe_ito"}, base R's
#'   colorblind-safe Okabe-Ito palette; \code{"tableau"} is also available.
#' @param colors Optional length-2 vector overriding the palette: the study
#'   color and the pooled-estimate color.
#'
#' @return A \pkg{ggplot} object; print it, or add further layers.
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso \code{\link{meta_es}}, \code{\link{meta_smd}}, and
#'   \code{\link{meta_r}} for the numbers behind the picture;
#'   \code{\link{teacher_expectancy}} for the example data.
#'
#' @family meta-analysis
#'
#' @family plotting
#'
#' @keywords hplot
#'
#' @examples
#' # Twelve simulated studies whose true effects vary from study to study
#' # (between-study standard deviation 0.35), so the prediction interval
#' # for the effect of a new study is visibly wider than the confidence
#' # interval for the mean effect.
#' set.seed(113)
#' k <- 12
#' n <- sample(20:100, k)                    # per-group sample sizes
#' theta <- rnorm(k, mean = 0.4, sd = 0.35)  # true study effects
#' d <- rnorm(k, mean = theta, sd = sqrt(2 / n))
#' v <- 2 / n + d^2 / (4 * n)
#' plot_forest(d, v, xlab = "Standardized mean difference (d)")
#'
#' # The teacher expectancy literature (Raudenbush, 1984): most studies
#' # cluster near zero, the estimated between-study variance is zero, and
#' # the prediction interval nearly coincides with the confidence interval.
#' data(teacher_expectancy)
#' d <- teacher_expectancy$d
#' n_e <- teacher_expectancy$n_experimental
#' n_c <- teacher_expectancy$n_control
#' v <- (n_e + n_c) / (n_e * n_c) + d^2 / (2 * (n_e + n_c))
#' plot_forest(d, v, labels = teacher_expectancy$author,
#'             xlab = "Standardized mean difference (d)")
#'
#' @export
plot_forest <- function(yi, vi, labels = NULL,
                        method = c("reml", "pm", "dl", "fe"),
                        hartung_knapp = TRUE, conf_level = 0.95,
                        show_prediction = TRUE,
                        xlab = "Effect size", title = NULL,
                        palette = "okabe_ito", colors = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plot_forest(). ",
         "Install it with install.packages(\"ggplot2\").", call. = FALSE)
  }
  .meta_check_yivi(yi, vi)
  method <- match.arg(method)
  k <- length(yi)
  if (is.null(labels)) labels <- paste("Study", seq_len(k))
  if (length(labels) != k) {
    stop("'labels' must supply one label per study.", call. = FALSE)
  }
  if (is.null(colors)) colors <- .dmar_palette(2, palette = palette)

  fit  <- .meta_fit(yi, vi, method = method, hartung_knapp = hartung_knapp,
                    conf_level = conf_level)
  crit <- stats::qnorm(1 - (1 - conf_level) / 2)

  pooled_label <- if (method == "fe") "Common effect" else "Random effects"
  dat <- data.frame(
    label  = factor(c(labels, pooled_label),
                    levels = rev(c(labels, pooled_label))),
    yi     = c(yi, fit$estimate),
    lower  = c(yi - crit * sqrt(vi), fit$lower),
    upper  = c(yi + crit * sqrt(vi), fit$upper),
    size   = c(1 / vi, max(1 / vi)),
    pooled = c(rep(FALSE, k), TRUE)
  )

  p <- ggplot2::ggplot(dat, ggplot2::aes(x = yi, y = label))

  if (show_prediction && method != "fe" &&
      is.finite(fit$prediction_lower)) {
    p <- p + ggplot2::annotate(
      "segment", x = fit$prediction_lower, xend = fit$prediction_upper,
      y = pooled_label, yend = pooled_label,
      linewidth = 4, color = colors[2], alpha = 0.25)
  }

  p +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed",
                        color = "gray55", linewidth = 0.4) +
    ggplot2::geom_errorbar(
      ggplot2::aes(xmin = lower, xmax = upper, color = pooled),
      width = 0.22, linewidth = 0.55, orientation = "y",
      show.legend = FALSE) +
    ggplot2::geom_point(
      ggplot2::aes(size = size, color = pooled, shape = pooled),
      show.legend = FALSE) +
    ggplot2::scale_color_manual(values = c(`FALSE` = colors[1],
                                            `TRUE`  = colors[2])) +
    ggplot2::scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 18)) +
    ggplot2::scale_size_continuous(range = c(1.6, 4.2)) +
    ggplot2::labs(x = xlab, y = NULL, title = title) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank())
}
