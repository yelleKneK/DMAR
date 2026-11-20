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
#' @param B Number of Monte Carlo draws behind the band, a single whole
#'   number. Defaults to 10000, because a draw costs one evaluation of each
#'   effect and no refit. The limits of the band are quantiles of the
#'   draws, so the fewest accepted is 400 at a 95 percent level (200 at
#'   0.90, 2000 at 0.99), 10 draws beyond each limit; a smaller \code{B} is
#'   refused (class \code{dmar_count_invalid}). A draw whose conditional
#'   effect is not a finite number at a value of the moderator is dropped
#'   there and counted; one warning (class
#'   \code{dmar_bootstrap_replications_dropped}) says so, and where more
#'   than half were dropped the band is not drawn (class
#'   \code{dmar_bootstrap_insufficient}).
#'   \code{options(dmar.bootstrap_exploratory = TRUE)} lets a smaller
#'   \code{B} run, with a warning of class
#'   \code{dmar_bootstrap_exploratory}. See \code{\link{dmar_bootstrap}}.
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
#'   on exit). Default \code{NULL} draws the band from the current random
#'   number stream, which advances that stream.
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
#'   recoverable from the object itself. The band limits are \code{NA}
#'   at a value of the moderator where no interval is formed (see
#'   \code{B}).
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
#' # w, so the indirect effect of x on y through m is a line in w. The
#' # simulated x reaches y only through m, and the model fit below
#' # carries no direct path, so the indirect effect is the whole effect
#' # of x and the picture has one curve.
#' set.seed(113)
#' n <- 300
#' x <- rnorm(n)
#' w <- rnorm(n)
#' m <- 0.5 * x + 0.3 * w + 0.4 * x * w + rnorm(n)
#' y <- 0.5 * m + 0.1 * w + rnorm(n)
#' d_mod <- data.frame(x = x, w = w, m = m, y = y)
#'
#' # Fit with the moderator declared. Every reported row costs its own
#' # constrained null model fit in OpenMx, so the Wald interval and two
#' # probe values keep the fit quick; the default probe values are the
#' # moderator's mean and one standard deviation either side, and the
#' # curve and its band cover the whole range of w either way.
#' model_mod <- "
#'   m ~ x + w + x:w
#'   y ~ m + w
#' "
#' res <- mediation_mbco(model_mod, data = d_mod, x = "x", y = "y",
#'                       moderator = "w", ci_method = "wald",
#'                       probe_values = c(low = -1, high = 1))
#'
#' # The band draws B coefficient vectors from their joint normal
#' # approximation. B = 2000 keeps the example quick; a reported figure
#' # deserves the default B = 10000. The probed values are marked on
#' # the curve, and the rug shows where w was observed.
#' plot_mediation_mbco(res, B = 2000, seed = 113)
#'
#' # The same curve over a chosen range of w, with a 90% band. The
#' # 'effects' argument names the pathway to draw; in a model that also
#' # carries a direct path the total effect is moderated too, and it is
#' # drawn as a second curve unless 'effects' selects one of them.
#' plot_mediation_mbco(res, effects = "indirect_via_m", from = -2,
#'                     to = 2, conf_level = 0.90, B = 2000, seed = 113)
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
    .dmar_refuse(c(
      "The package 'ggplot2' is needed; please install the ",
      "package and try again."), "package_required", package = "ggplot2")
  }
  if (!inherits(x, "dmar_mediation_mbco")) {
    .dmar_refuse_argument(
      "'x' must be a result from mediation_mbco().", arg = "x")
  }
  modn <- attr(x, "moderation")
  if (is.null(modn) || length(modn$effects) == 0L) {
    .dmar_refuse_argument(c(
      "This result was not fit with a 'moderator', so there are ",
      "no conditional effects to draw; refit mediation_mbco() ",
      "with the 'moderator' argument."), arg = NULL)
  }
  fit <- attr(x, "mx_model")
  if (is.null(fit)) {
    .dmar_refuse_argument(c(
      "The fitted model (the \"mx_model\" attribute) is missing; ",
      "plot the object as returned by mediation_mbco()."), arg = NULL)
  }
  if (is.null(conf_level)) conf_level <- attr(x, "conf_level")
  if (is.null(conf_level)) conf_level <- 0.95
  if (!is.numeric(conf_level) || length(conf_level) != 1L ||
      is.na(conf_level) || conf_level <= 0 || conf_level >= 1) {
    .dmar_refuse_argument(
      "'conf_level' must be a single number in (0, 1).", arg = "conf_level")
  }
  # The band is read from quantiles of the B draws, so B follows the rule
  # of ?dmar_bootstrap: 10 draws beyond each limit (400 at the 95 percent
  # level).
  B <- .dmar_check_B(B, tail = (1 - conf_level) / 2,
                     conf_level = conf_level)
  if (!is.numeric(n_grid) || length(n_grid) != 1L || is.na(n_grid) ||
      n_grid < 2 || n_grid != round(n_grid)) {
    .dmar_refuse_argument(
      "'n_grid' must be a single integer of at least 2.", arg = "n_grid")
  }
  available <- modn$effects
  if (is.null(effects)) effects <- names(available)
  unknown <- setdiff(effects, names(available))
  if (length(unknown) > 0L) {
    .dmar_refuse_argument(c(
      "Unknown effect(s): ", paste(sQuote(unknown), collapse = ", "),
      ". Moderated effects available: ",
      paste(sQuote(names(available)), collapse = ", "), "."), arg = NULL)
  }
  sel <- available[effects]

  est <- coef(fit)
  V <- try(suppressWarnings(vcov(fit)), silent = TRUE)
  if (inherits(V, "try-error")) {
    .dmar_refuse_argument(c(
      "The parameter covariance matrix of the fitted model is ",
      "unavailable, so the confidence band cannot be computed."), arg = NULL)
  }
  w_obs <- fit@data@observed[[modn$moderator]]
  # When mediation_mbco() fit badly scaled variables in a rescaled metric
  # (its "rescale" attribute), the estimates and the moderator are returned
  # to the original metric, in which the conditional effects are written.
  # The draws below are taken in the metric of the fit, whose covariance
  # matrix is well conditioned, and then multiplied, as mediation_mbco()
  # does: drawn from the covariance matrix of the original metric, whose
  # variances can span 25 orders of magnitude, the draws of the smallest
  # coefficients are lost in rounding.
  mult_est <- rep(1, length(est))
  resc <- attr(x, "rescale", exact = TRUE)
  if (!is.null(resc)) {
    hit <- names(est) %in% names(resc$multipliers)
    mult_est[hit] <- resc$multipliers[names(est)[hit]]
    if (modn$moderator %in% names(resc$divisors)) {
      w_obs <- w_obs * resc$divisors[[modn$moderator]]
    }
  }
  est_fit <- est
  est <- est * mult_est
  w_obs <- w_obs[is.finite(w_obs)]
  if (is.null(from)) from <- min(w_obs)
  if (is.null(to)) to <- max(w_obs)
  if (!is.numeric(from) || !is.numeric(to) || length(from) != 1L ||
      length(to) != 1L || is.na(from) || is.na(to) || from >= to) {
    .dmar_refuse_argument(
      "'from' must be smaller than 'to'.", arg = c("from", "to"))
  }
  grid <- seq(from, to, length.out = n_grid)

  # A supplied seed is local to this call; the caller's random number
  # generator state is restored on exit.
  .dmar_local_seed(seed)
  draws <- as.data.frame(MASS::mvrnorm(B, mu = est_fit, Sigma = V))
  if (!is.null(resc)) {
    draws[] <- Map(`*`, draws, mult_est[match(names(draws), names(est))])
  }
  est_env <- as.list(est)
  alpha <- 1 - conf_level
  probs <- c(alpha / 2, 1 - alpha / 2)

  band <- NULL
  probe_pts <- NULL
  # The fewest finite draws at any point of any band, for the rule.
  fewest <- B
  usable <- function(count) {
    .dmar_boot_settle(
      B, count, tail = alpha / 2, conf_level = conf_level,
      fun = "plot_mediation_mbco",
      quantity = "the conditional effect at some value of the moderator",
      dropped_because = paste0("gave a conditional effect that is not a ",
                               "finite number"),
      replications = "Monte Carlo draws", signal = FALSE, recorded = FALSE)
  }
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
    # A draw whose conditional effect is not a finite number at a value
    # of the moderator is dropped there and counted; where the rule of
    # R/dmar_bootstrap.R leaves no interval the band is not drawn.
    counts <- colSums(is.finite(E_mat))
    fewest <- min(fewest, counts)
    lims <- vapply(seq_len(ncol(E_mat)), function(j) {
      if (!usable(counts[j])$formed) return(c(NA_real_, NA_real_))
      v <- E_mat[, j]
      quantile(v[is.finite(v)], probs = probs, names = FALSE)
    }, numeric(2L))
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
  # One warning for the plot, from the point with the fewest finite draws:
  # dropped draws, a band withheld there, or a run below the required
  # number under options(dmar.bootstrap_exploratory = TRUE).
  .dmar_boot_signal(usable(fewest))
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
      caption = paste0(.format_conf_pct(conf_level), "% pointwise Monte ",
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
