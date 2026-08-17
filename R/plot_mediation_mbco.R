# Suppress R CMD check notes for ggplot2 aes() column references.
if (getRversion() >= "2.15.1") {
  utils::globalVariables(c("w_value", "effect_label", "band_lower",
                           "band_upper"))
}

#' Plot Conditional Effects From a Moderated Mediation Analysis
#'
#' Draws the conditional effects from a \code{\link{mediation_mbco}}
#' analysis that declared a \code{moderator}: for each moderated
#' pathway effect, the curve tracing how the effect changes over the
#' moderator's range, with a pointwise confidence band, the probed
#' values marked, a dashed reference line at zero, and a rug showing
#' where the moderator was actually observed. The picture answers, at
#' a glance, the questions the table answers row by row: how large is
#' the effect at any given moderator value, where (if anywhere) does
#' its interval exclude zero, and over what part of the moderator's
#' range the data can support either statement.
#'
#' @param x A \code{dmar_mediation_mbco} object returned by
#'   \code{\link{mediation_mbco}} with a \code{moderator}. An object
#'   fit without a moderator has no conditional effects to draw, and
#'   the function says so.
#' @param effects Character vector naming which moderated effects to
#'   draw, using the base effect names from the result table (e.g.,
#'   \code{"indirect_via_m"}, \code{"total_effect"}). Defaults to all
#'   moderated effects. Unmoderated effects are flat lines and are not
#'   drawn.
#' @param conf_level Confidence level for the band. Defaults to the
#'   level used when the object was fit.
#' @param B Number of Monte Carlo draws behind the band. Defaults to
#'   10000.
#' @param from,to Range of moderator values to draw. Defaults to the
#'   observed range of the moderator. Values outside the observed
#'   range are extrapolation; the rug makes that visible.
#' @param n_grid Number of grid points along the moderator at which
#'   the curve and band are evaluated. Defaults to 200.
#' @param show_probe_values Logical. If \code{TRUE} (default), mark
#'   the probed moderator values (the \code{_at_} rows of the result
#'   table) as points on each curve.
#' @param show_rug Logical. If \code{TRUE} (default), draw a rug of
#'   the observed moderator values along the horizontal axis.
#' @param palette Character string naming the color palette. Defaults
#'   to \code{"okabe_ito"}, base R's colorblind-safe Okabe-Ito
#'   palette; \code{"tableau"} is also available.
#' @param xlab,ylab,title Optional axis labels and title. The defaults
#'   name the moderator on the horizontal axis and describe the
#'   vertical axis as the conditional effect of \code{x} on \code{y}.
#' @param seed Optional integer seed for the Monte Carlo band, used
#'   locally (the caller's random number generator state is restored
#'   on exit). Default \code{NULL} leaves the random number generator
#'   state alone.
#'
#' @details
#' \strong{What is drawn, and where it comes from.} A pathway effect
#' in a model with interactions is a polynomial in the moderator: a
#' straight line when the pathway is moderated in one place (its slope
#' is the index of moderated mediation), a curve when it is moderated
#' in more than one. \code{\link{mediation_mbco}} derives each
#' polynomial symbolically and stores it with the result, so this
#' function evaluates the same quantity the table probes, just
#' everywhere in the moderator's range instead of at two or three
#' values. The marked points are exactly the table's \code{_at_} rows.
#'
#' \strong{The band is pointwise.} At each grid value of the
#' moderator, the band is a \code{conf_level} Monte Carlo confidence
#' interval for the conditional effect at that one value: the path
#' coefficients are drawn from their joint normal approximation
#' (MacKinnon, Lockwood, & Williams, 2004), each draw's polynomial is
#' evaluated along the grid, and the band connects the pointwise
#' quantiles. Read vertically at a single moderator value of interest,
#' it is an ordinary confidence interval. Read horizontally, the
#' moderator values where the band crosses zero estimate the
#' Johnson-Neyman boundaries (Johnson & Neyman, 1936; Preacher,
#' Rucker, & Hayes, 2007), the values separating "interval excludes
#' zero" from "interval includes zero". That horizontal reading scans
#' many intervals at once, so the pointwise band understates the
#' uncertainty of the boundary locations themselves; treat the
#' crossing points as estimates, not as sharp cutoffs, and lean on the
#' table's moderation and constancy tests for the formal question of
#' whether the effect depends on the moderator at all.
#'
#' \strong{The rug guards against extrapolation.} The curve can be
#' evaluated at any moderator value, but the data only inform it where
#' the moderator was observed. The rug shows that support directly; a
#' confident-looking band in a region with no rug beneath it is
#' arithmetic, not evidence.
#'
#' \strong{The band and the table may differ slightly.} The band is
#' always Monte Carlo, whichever \code{ci_method} the table used. At a
#' probed value, a Monte Carlo band and a profile likelihood or Wald
#' interval agree closely in large samples but are not the same
#' construction; small discrepancies between the band and an
#' \code{_at_} row's interval are expected, not a defect.
#'
#' The plot is an ordinary \pkg{ggplot2} object, so any further
#' customization (themes, additional layers, institutional color
#' scales) can be added to the returned value with \code{+}.
#'
#' @return A \code{ggplot2} object. Its data contains one row per
#'   effect and grid value with columns \code{effect_label},
#'   \code{w_value}, \code{estimate}, \code{band_lower}, and
#'   \code{band_upper}, so the numbers behind the picture are
#'   recoverable from the object itself.
#'
#' @note Requires \pkg{ggplot2} (listed in \code{Suggests}).
#'
#' @references
#' Johnson, P. O., & Neyman, J. (1936). Tests of certain linear
#'   hypotheses and their application to some educational problems.
#'   \emph{Statistical Research Memoirs, 1}, 57--93.
#'
#' MacKinnon, D. P., Lockwood, C. M., & Williams, J. (2004). Confidence
#'   limits for the indirect effect: Distribution of the product and
#'   resampling methods. \emph{Multivariate Behavioral Research, 39}(1),
#'   99--128. \doi{10.1207/s15327906mbr3901_4}
#'
#' Preacher, K. J., Rucker, D. D., & Hayes, A. F. (2007). Addressing
#'   moderated mediation hypotheses: Theory, methods, and
#'   prescriptions. \emph{Multivariate Behavioral Research, 42}(1),
#'   185--227. \doi{10.1080/00273170701341316}
#'
#' Tofighi, D., & Kelley, K. (2020). Improved inference in mediation
#'   analysis: Introducing the model-based constrained optimization
#'   procedure. \emph{Psychological Methods, 25}(4), 496--515.
#'   \doi{10.1037/met0000259}
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso \code{\link{mediation_mbco}} for the analysis this
#'   function displays; \code{\link{regions_of_significance}} for the
#'   analogous display for mixed-effects model interactions.
#'
#' @family plotting
#'
#' @keywords hplot
#'
#' @examples
#' # First-stage moderated mediation: the effect of x on m depends on
#' # w, so the indirect effect of x on y is a line in w; the direct
#' # effect is unmoderated and is not drawn.
#' set.seed(113)
#' n <- 300
#' x <- rnorm(n)
#' w <- rnorm(n)
#' m <- 0.5 * x + 0.3 * w + 0.4 * x * w + rnorm(n)
#' y <- 0.5 * m + 0.2 * x + 0.1 * w + rnorm(n)
#' d_mod <- data.frame(x = x, w = w, m = m, y = y)
#'
#' # Neither the fit nor the plot is run here: every probed effect costs
#' # its own constrained null model fit in OpenMx, and the band draws B
#' # coefficient vectors from their joint normal approximation. The Wald
#' # interval and the two probe values keep a hand run of these lines
#' # quick; the default probe values are the moderator's mean and one
#' # standard deviation either side, and the curve and its band cover
#' # the whole range of w either way. The calls are:
#' # res <- mediation_mbco("m ~ x + w + x:w \n y ~ m + x + w",
#' #                       data = d_mod, x = "x", y = "y",
#' #                       moderator = "w", ci_method = "wald",
#' #                       probe_values = c(low = -1, high = 1))
#' # plot_mediation_mbco(res, seed = 113)
#'
#' # Only the indirect pathway, over a chosen moderator range:
#' # plot_mediation_mbco(res, effects = "indirect_via_m",
#' #                     from = -2, to = 2, seed = 113)
#'
#' @export
#' @importFrom stats coef vcov quantile
plot_mediation_mbco <- function(x, effects = NULL, conf_level = NULL,
                                B = 10000, from = NULL, to = NULL,
                                n_grid = 200,
                                show_probe_values = TRUE,
                                show_rug = TRUE,
                                palette = c("okabe_ito", "tableau"),
                                xlab = NULL, ylab = NULL, title = NULL,
                                seed = NULL) {
  palette <- match.arg(palette)
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("The package 'ggplot2' is needed; please install the ",
         "package and try again.", call. = FALSE)
  }
  if (!inherits(x, "dmar_mediation_mbco")) {
    stop("'x' must be a result from mediation_mbco().", call. = FALSE)
  }
  modn <- attr(x, "moderation")
  if (is.null(modn) || length(modn$effects) == 0L) {
    stop("This result was not fit with a 'moderator', so there are ",
         "no conditional effects to draw; refit mediation_mbco() ",
         "with the 'moderator' argument.", call. = FALSE)
  }
  fit <- attr(x, "mx_model")
  if (is.null(fit)) {
    stop("The fitted model (the \"mx_model\" attribute) is missing; ",
         "plot the object as returned by mediation_mbco().",
         call. = FALSE)
  }
  if (is.null(conf_level)) conf_level <- attr(x, "conf_level")
  if (is.null(conf_level)) conf_level <- 0.95
  if (!is.numeric(conf_level) || length(conf_level) != 1L ||
      is.na(conf_level) || conf_level <= 0 || conf_level >= 1) {
    stop("'conf_level' must be a single number in (0, 1).",
         call. = FALSE)
  }
  if (!is.numeric(B) || length(B) != 1L || is.na(B) || B < 100 ||
      B != round(B)) {
    stop("'B' must be a single integer of at least 100.", call. = FALSE)
  }
  if (!is.numeric(n_grid) || length(n_grid) != 1L || is.na(n_grid) ||
      n_grid < 2 || n_grid != round(n_grid)) {
    stop("'n_grid' must be a single integer of at least 2.",
         call. = FALSE)
  }
  available <- modn$effects
  if (is.null(effects)) effects <- names(available)
  unknown <- setdiff(effects, names(available))
  if (length(unknown) > 0L) {
    stop("Unknown effect(s): ", paste(sQuote(unknown), collapse = ", "),
         ". Moderated effects available: ",
         paste(sQuote(names(available)), collapse = ", "), ".",
         call. = FALSE)
  }
  sel <- available[effects]

  est <- coef(fit)
  V <- try(suppressWarnings(vcov(fit)), silent = TRUE)
  if (inherits(V, "try-error")) {
    stop("The parameter covariance matrix of the fitted model is ",
         "unavailable, so the confidence band cannot be computed.",
         call. = FALSE)
  }
  w_obs <- fit@data@observed[[modn$moderator]]
  w_obs <- w_obs[is.finite(w_obs)]
  if (is.null(from)) from <- min(w_obs)
  if (is.null(to)) to <- max(w_obs)
  if (!is.numeric(from) || !is.numeric(to) || length(from) != 1L ||
      length(to) != 1L || is.na(from) || is.na(to) || from >= to) {
    stop("'from' must be smaller than 'to'.", call. = FALSE)
  }
  grid <- seq(from, to, length.out = n_grid)

  # Local RNG: seed if asked, and always restore the caller's state.
  if (!is.null(seed)) {
    has_old <- exists(".Random.seed", envir = globalenv())
    old <- if (has_old) get(".Random.seed", envir = globalenv())
    on.exit({
      if (has_old) {
        assign(".Random.seed", old, envir = globalenv())
      } else if (exists(".Random.seed", envir = globalenv())) {
        rm(".Random.seed", envir = globalenv())
      }
    }, add = TRUE)
    set.seed(seed)
  }
  draws <- as.data.frame(MASS::mvrnorm(B, mu = est, Sigma = V))
  est_env <- as.list(est)
  alpha <- 1 - conf_level
  probs <- c(alpha / 2, 1 - alpha / 2)

  band <- NULL
  probe_pts <- NULL
  labels <- unname(vapply(sel, `[[`, character(1L), "pathway"))
  if (anyDuplicated(labels)) {
    labels <- paste0(labels, " [", names(sel), "]")
  }
  for (i in seq_along(sel)) {
    e <- sel[[i]]
    d_pow <- seq_along(e$coefs) - 1L
    W_pow <- outer(grid, d_pow, `^`)
    coef_hat <- vapply(e$coefs, function(cx) {
      eval(parse(text = cx)[[1L]], envir = est_env)
    }, numeric(1L))
    curve <- as.numeric(W_pow %*% coef_hat)
    C_draw <- vapply(e$coefs, function(cx) {
      eval(parse(text = cx)[[1L]], envir = draws)
    }, numeric(B))
    E_mat <- C_draw %*% t(W_pow)
    lims <- apply(E_mat, 2L, quantile, probs = probs, names = FALSE)
    band <- rbind(band, data.frame(
      effect_label = labels[i], w_value = grid, estimate = curve,
      band_lower = lims[1L, ], band_upper = lims[2L, ],
      stringsAsFactors = FALSE))
    if (show_probe_values) {
      P_pow <- outer(unname(modn$values), d_pow, `^`)
      probe_pts <- rbind(probe_pts, data.frame(
        effect_label = labels[i], w_value = unname(modn$values),
        estimate = as.numeric(P_pow %*% coef_hat),
        stringsAsFactors = FALSE))
    }
  }
  band$effect_label <- factor(band$effect_label, levels = labels)
  if (!is.null(probe_pts)) {
    probe_pts$effect_label <- factor(probe_pts$effect_label,
                                     levels = labels)
  }

  if (is.null(xlab)) xlab <- modn$moderator
  if (is.null(ylab)) {
    ylab <- paste0("Conditional effect of ", attr(x, "x"), " on ",
                   attr(x, "y"))
  }
  cols <- .dmar_palette(length(sel), palette = palette)
  p <- ggplot2::ggplot(band, ggplot2::aes(
    x = w_value, y = estimate, color = effect_label,
    fill = effect_label)) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed",
                        color = "grey40") +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = band_lower,
                                      ymax = band_upper),
                         alpha = 0.2, linewidth = 0) +
    ggplot2::geom_line(linewidth = 0.9) +
    ggplot2::scale_color_manual(values = cols, name = "Effect") +
    ggplot2::scale_fill_manual(values = cols, name = "Effect") +
    ggplot2::labs(
      x = xlab, y = ylab, title = title,
      caption = paste0(round(conf_level * 100), "% pointwise Monte ",
                       "Carlo confidence band"))
  if (show_probe_values && !is.null(probe_pts)) {
    p <- p + ggplot2::geom_point(data = probe_pts, size = 2)
  }
  if (show_rug && length(w_obs) > 0L) {
    p <- p + ggplot2::geom_rug(
      data = data.frame(w_value = w_obs),
      ggplot2::aes(x = w_value),
      inherit.aes = FALSE, sides = "b", alpha = 0.25)
  }
  p
}
