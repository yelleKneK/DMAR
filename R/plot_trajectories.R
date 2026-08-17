# Internal: resolve which IDs to plot from one-of {ids, n_random, pct_random}.
.resolve_ids <- function(unique_ids, ids, n_random, pct_random) {
  spec_count <- sum(!is.null(ids), !is.null(n_random), !is.null(pct_random))
  if (spec_count > 1L) {
    stop("Specify at most one of 'ids', 'n_random', or 'pct_random'.", call. = FALSE)
  }

  if (!is.null(ids)) {
    missing_ids <- setdiff(ids, unique_ids)
    if (length(missing_ids) > 0L) {
      stop("ID(s) not present in the data: ",
           paste(missing_ids, collapse = ", "), ".", call. = FALSE)
    }
    return(ids)
  }

  if (!is.null(n_random)) {
    if (!is.numeric(n_random) || length(n_random) != 1L ||
        n_random <= 0 || n_random > length(unique_ids)) {
      stop("'n_random' must be a single integer between 1 and the number of unique IDs (",
           length(unique_ids), ").", call. = FALSE)
    }
    return(sample(unique_ids, n_random))
  }

  if (!is.null(pct_random)) {
    if (!is.numeric(pct_random) || length(pct_random) != 1L ||
        pct_random <= 0 || pct_random > 100) {
      stop("'pct_random' must be a single number in (0, 100].", call. = FALSE)
    }
    pct <- if (pct_random > 1) pct_random / 100 else pct_random
    n <- ceiling(length(unique_ids) * pct)
    return(sample(unique_ids, n))
  }

  unique_ids
}


# Visualize observed individual trajectories in a longitudinal context
#' Visualize Observed Individual Trajectories in a Longitudinal Data Set
#'
#' Plots one trajectory per subject from a long-format data frame, optionally
#' colored by a grouping variable, and optionally faceted into one panel per
#' subject. Returns a \pkg{ggplot2} object that can be further customized.
#'
#' @param data A long-format \code{data.frame} (one row per subject-occasion).
#' @param id Character. Column name in \code{data} identifying the subject.
#' @param time Character. Column name for the time / occasion variable
#'   (the \emph{x} axis).
#' @param outcome Character. Column name for the outcome / score variable
#'   (the \emph{y} axis).
#' @param group Optional character. Column name for a grouping variable used
#'   to color the trajectories (and panels, if faceted).
#' @param ids Optional vector of subject IDs to plot.
#' @param n_random Optional integer; randomly sample this many subjects.
#' @param pct_random Optional numeric; sample this percentage of subjects.
#'   Values \eqn{\leq 1} are interpreted as proportions; values
#'   \eqn{> 1} as percentages. At most one of \code{ids}, \code{n_random},
#'   or \code{pct_random} may be supplied.
#' @param facet Logical. If \code{TRUE}, draw one panel per subject via
#'   \code{\link[ggplot2]{facet_wrap}}. If \code{FALSE} (default), overlay all
#'   trajectories in one panel.
#' @param nrow,ncol Optional integers passed to \code{facet_wrap()} when
#'   \code{facet = TRUE}.
#' @param show_points Logical. If \code{TRUE} (default), draw the observed
#'   points as well as the connecting lines.
#' @param point_size Size of the observed points (default \code{1.5}).
#' @param linewidth Line width for the connecting segments (default \code{0.5}).
#' @param alpha Transparency for points and lines (default \code{0.7}).
#' @param palette Character string naming the color palette used to color the
#'   trajectories when \code{group} is a discrete (factor, character, or
#'   logical) variable. Defaults to \code{"okabe_ito"}, base R's
#'   colorblind-safe Okabe-Ito palette; \code{"tableau"} is also available.
#'   Ignored when \code{group} is \code{NULL} or numeric.
#' @param title,xlab,ylab Optional plot labels. Sensible defaults are taken
#'   from \code{outcome} and \code{time} when these are \code{NULL}.
#' @param seed Optional integer random seed used when
#'   \code{n_random} or \code{pct_random} is supplied. Defaults to
#'   \code{NULL}, which leaves the user's current RNG state intact;
#'   supply an integer for reproducible subject sampling.
#'
#' @return A \code{ggplot} object.
#'
#' @details
#' The function modernizes the original \code{vit()} (visualize individual
#' trajectories) function by returning a single \pkg{ggplot2} object instead
#' of producing graphical side effects. Saving is handled by the user via
#' \code{\link[ggplot2]{ggsave}}; multi-page output via faceting and
#' \code{\link[ggplot2]{facet_wrap}}'s \code{nrow}/\code{ncol}.
#'
#' @note Requires \pkg{ggplot2} (a \code{Suggests} dependency).
#'
#' @seealso \code{\link{plot_trajectories_fitted}} for plotting observed
#'   trajectories together with a fitted multilevel model's predictions.
#'
#' @examples
#' # Built-in Orthodont data: 27 children, 4 measurements each.
#' d <- nlme::Orthodont
#'
#' # Overlay all trajectories, colored by sex.
#' plot_trajectories(d, id = "Subject", time = "age",
#'                   outcome = "distance", group = "Sex")
#'
#' # One panel per child, for twelve children drawn at random. Not run
#' # here because faceting draws twelve small plots instead of one, which
#' # costs about twice what the overlay above does. The call is:
#' # plot_trajectories(d, id = "Subject", time = "age",
#' #                   outcome = "distance",
#' #                   n_random = 12, facet = TRUE, ncol = 4,
#' #                   seed = 113)
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @keywords hplot
#'
#' @family plotting
#'
#' @export
plot_trajectories <- function(data, id, time, outcome,
                              group       = NULL,
                              ids         = NULL,
                              n_random    = NULL,
                              pct_random  = NULL,
                              facet       = FALSE,
                              nrow        = NULL,
                              ncol        = NULL,
                              show_points = TRUE,
                              point_size  = 1.5,
                              linewidth   = 0.5,
                              alpha       = 0.7,
                              palette     = "okabe_ito",
                              title       = NULL,
                              xlab        = NULL,
                              ylab        = NULL,
                              seed        = NULL) {
  if (!is.null(seed) &&
      (!is.null(n_random) || !is.null(pct_random))) {
    if (exists(".Random.seed", envir = .GlobalEnv)) {
      .old_seed <- get(".Random.seed", envir = .GlobalEnv)
      on.exit(assign(".Random.seed", .old_seed, envir = .GlobalEnv), add = TRUE)
    } else {
      on.exit(if (exists(".Random.seed", envir = .GlobalEnv)) rm(list = ".Random.seed", envir = .GlobalEnv), add = TRUE)
    }
    set.seed(seed)
  }
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required. Install with install.packages(\"ggplot2\").",
         call. = FALSE)
  }
  if (!is.data.frame(data)) stop("'data' must be a data.frame.", call. = FALSE)

  # Column-name validation.
  for (col in c(id, time, outcome)) {
    if (!is.character(col) || length(col) != 1L) {
      stop("'id', 'time', and 'outcome' must each be a single column name (character).",
           call. = FALSE)
    }
    if (!col %in% names(data)) {
      stop("Column '", col, "' not found in 'data'.", call. = FALSE)
    }
  }
  if (!is.null(group)) {
    if (!is.character(group) || length(group) != 1L || !group %in% names(data)) {
      stop("'group' must be a single column name found in 'data'.", call. = FALSE)
    }
  }

  unique_ids   <- unique(data[[id]])
  selected_ids <- .resolve_ids(unique_ids, ids, n_random, pct_random)
  data_sub     <- data[data[[id]] %in% selected_ids, , drop = FALSE]

  p <- ggplot2::ggplot(
    data_sub,
    ggplot2::aes(x = .data[[time]], y = .data[[outcome]],
                 group = .data[[id]])
  )

  if (!is.null(group)) {
    p <- p + ggplot2::aes(color = .data[[group]])
    # Apply the DMAR palette only to a discrete grouping variable; a
    # numeric group keeps ggplot2's default continuous color scale.
    group_vals <- data_sub[[group]]
    if (is.factor(group_vals) || is.character(group_vals) ||
        is.logical(group_vals)) {
      p <- p + .dmar_discrete_scale(
        "colour",
        function(k) .dmar_palette(k, palette = palette)
      )
    }
  }

  if (isTRUE(show_points)) {
    p <- p + ggplot2::geom_point(alpha = alpha, size = point_size)
  }
  p <- p + ggplot2::geom_line(alpha = alpha, linewidth = linewidth)

  if (isTRUE(facet)) {
    p <- p + ggplot2::facet_wrap(
      stats::reformulate(id), nrow = nrow, ncol = ncol
    )
  }

  p <- p + ggplot2::labs(
    title = title,
    x     = if (is.null(xlab)) time    else xlab,
    y     = if (is.null(ylab)) outcome else ylab
  ) +
    ggplot2::theme_minimal(base_size = 12)

  p
}


# ---------------------------------------------------------------------------
# Internal: extract observed and fitted trajectories from a multilevel model.
# Returns a list with:
#   observed:  data.frame with columns {id, time, outcome, .fitted}
#   fitted:    data.frame with columns {id, time, .fitted} on a smooth grid
#   id_col, time_col, outcome_col: detected/specified column names
.extract_traj_info <- function(model, id_col, time_col, outcome_col, n_grid) {
  if (inherits(model, c("lme", "nlme"))) {
    if (!requireNamespace("nlme", quietly = TRUE)) {
      stop("Package 'nlme' is required for nlme/lme model objects.", call. = FALSE)
    }
    mf <- nlme::getData(model)
    if (is.null(mf)) {
      stop("Could not retrieve the data frame from the nlme/lme fit. ",
           "Refit with 'data = your_data' explicitly.", call. = FALSE)
    }
    if (is.null(outcome_col)) outcome_col <- as.character(stats::formula(model)[[2L]])
    if (is.null(id_col))      id_col      <- names(model$groups)[1L]
    # If only one fixed predictor, default time to it; otherwise require user.
    if (is.null(time_col)) {
      fixed_terms <- attr(stats::terms(stats::formula(model), fixed.only = FALSE),
                          "term.labels")
      time_col <- fixed_terms[1L]
    }
  } else if (inherits(model, c("lmerMod", "lmerModLmerTest"))) {
    if (!requireNamespace("lme4", quietly = TRUE)) {
      stop("Package 'lme4' is required for lmer model objects.", call. = FALSE)
    }
    mf <- model@frame
    if (is.null(outcome_col)) outcome_col <- names(mf)[1L]
    if (is.null(id_col)) {
      # findbars() moved from lme4 to reformulas (which lme4 imports, so it
      # is present whenever this branch runs); calling the lme4 shim warns.
      find_bars <- if (requireNamespace("reformulas", quietly = TRUE))
        reformulas::findbars else lme4::findbars
      bars <- find_bars(stats::formula(model))
      if (length(bars) == 0L) stop("Model has no random-effect grouping factor.",
                                   call. = FALSE)
      id_col <- as.character(bars[[1L]][[3L]])
    }
    if (is.null(time_col)) {
      fixed_terms <- attr(stats::terms(stats::formula(model), fixed.only = TRUE),
                          "term.labels")
      time_col <- fixed_terms[1L]
    }
  } else {
    stop("'model' must be an aov/lme/nlme fit (package nlme) or an lmer fit ",
         "(package lme4). Got class: ", paste(class(model), collapse = "/"), ".",
         call. = FALSE)
  }

  for (col in c(id_col, time_col, outcome_col)) {
    if (!col %in% names(mf)) {
      stop("Could not find column '", col, "' in the model frame; ",
           "supply it explicitly via id/time/outcome.", call. = FALSE)
    }
  }

  obs <- data.frame(
    .id      = mf[[id_col]],
    .time    = mf[[time_col]],
    .outcome = mf[[outcome_col]],
    .fitted  = unname(stats::fitted(model)),
    stringsAsFactors = FALSE
  )

  # Build a smooth per-id grid for fitted curves.
  ids_unique <- unique(obs$.id)
  grid_list <- lapply(ids_unique, function(i) {
    rng <- range(obs$.time[obs$.id == i])
    data.frame(
      .id   = rep(i, n_grid),
      .time = seq(rng[1L], rng[2L], length.out = n_grid),
      stringsAsFactors = FALSE
    )
  })
  grid <- do.call(rbind, grid_list)
  names(grid)[names(grid) == ".id"]   <- id_col
  names(grid)[names(grid) == ".time"] <- time_col

  pred <- tryCatch({
    if (inherits(model, c("lme", "nlme"))) {
      stats::predict(model, newdata = grid, level = 1L)
    } else {
      stats::predict(model, newdata = grid, re.form = NULL,
                     allow.new.levels = FALSE)
    }
  }, error = function(e) NULL)

  if (is.null(pred)) {
    # Fallback: just connect the fitted values at observed times.
    fit_df <- obs[, c(".id", ".time", ".fitted")]
    names(fit_df) <- c(id_col, time_col, ".fitted")
    fit_df <- fit_df[order(fit_df[[id_col]], fit_df[[time_col]]), ]
  } else {
    fit_df <- grid
    fit_df$.fitted <- as.numeric(pred)
  }

  obs_out <- obs
  names(obs_out) <- c(id_col, time_col, outcome_col, ".fitted")

  list(
    observed     = obs_out,
    fitted       = fit_df,
    id_col       = id_col,
    time_col     = time_col,
    outcome_col  = outcome_col
  )
}


# Internal: per-id R^2 (sample correlation squared between observed and fitted)
# and root-mean-square error (sqrt of mean squared residuals).
.per_id_quality_of_fit <- function(observed, id_col, outcome_col) {
  ids <- unique(observed[[id_col]])
  do.call(rbind, lapply(ids, function(i) {
    rows <- observed[[id_col]] == i
    y    <- observed[[outcome_col]][rows]
    yhat <- observed$.fitted[rows]
    if (length(y) < 2L || stats::sd(y) == 0 || stats::sd(yhat) == 0) {
      r2 <- NA_real_
    } else {
      r2 <- stats::cor(y, yhat)^2
    }
    rmse <- sqrt(mean((y - yhat)^2))
    df_row <- data.frame(.id = i, r_squared = r2, rmse = rmse,
                         stringsAsFactors = FALSE)
    names(df_row)[1L] <- id_col
    df_row
  }))
}


# Visualize individual trajectories with the fitted curve from a multilevel model
#' Plot Observed and Fitted Individual Trajectories From a Multilevel Model
#'
#' Given a fitted \code{lme}/\code{nlme} (\pkg{nlme}) or \code{lmer}
#' (\pkg{lme4}) model, plots each subject's observed values and fitted curve
#' on a smooth time grid, faceted one panel per subject. Per-subject \eqn{R^2}
#' (squared correlation between observed and fitted values) and root-mean-square
#' error are computed and attached to the returned \pkg{ggplot2} object as the
#' \code{quality_of_fit} attribute.
#'
#' @param model A fitted model object of class \code{lme}, \code{nlme}, or
#'   \code{lmerMod}.
#' @param id,time,outcome Optional character names of the ID, time, and
#'   outcome columns. When \code{NULL} (the default), each is auto-detected
#'   from the model object: the outcome is the response variable in the
#'   formula, the ID is the first random-effect grouping factor, and time is
#'   the first fixed-effect predictor. Override these when the auto-detection
#'   is wrong.
#' @param ids,n_random,pct_random Subject-subsetting options identical to
#'   those of \code{\link{plot_trajectories}}; at most one may be supplied.
#' @param n_grid Integer. Number of points used to draw each subject's smooth
#'   fitted curve (default \code{100}).
#' @param show_points Logical. Whether to draw the observed values
#'   (default \code{TRUE}).
#' @param point_size,linewidth,alpha,nrow,ncol Visual / layout controls.
#' @param palette Character string naming the color palette; the fitted curve
#'   is drawn in the palette's primary color. Defaults to \code{"okabe_ito"},
#'   base R's colorblind-safe Okabe-Ito palette; \code{"tableau"} is also
#'   available.
#' @param show_quality Logical. If \code{TRUE} (the default), each panel
#'   strip text includes the subject's \eqn{R^2} and RMSE.
#' @param title,xlab,ylab Optional plot labels.
#' @param seed Optional integer random seed used when
#'   \code{n_random} or \code{pct_random} is supplied. Defaults to
#'   \code{NULL}, which leaves the user's current RNG state intact;
#'   supply an integer for reproducible subject sampling.
#'
#' @return A \code{ggplot} object. The per-subject quality-of-fit
#'   \code{data.frame} (columns: id column, \code{r_squared}, \code{rmse}) is
#'   attached as \code{attr(<plot>, "quality_of_fit")}.
#'
#' @details Modernizes the original \code{vit_fitted()} function by:
#' \itemize{
#'   \item returning a \pkg{ggplot2} object instead of writing to graphics
#'         devices,
#'   \item attaching per-subject quality-of-fit as an attribute rather than
#'         assigning it to the global environment via \code{<<-} (a serious
#'         side effect of the original),
#'   \item drawing a smooth fitted curve from a per-subject time grid via
#'         \code{predict(..., re.form = NULL)} for \pkg{lme4} fits and
#'         \code{predict(..., level = 1)} for \pkg{nlme} fits,
#'   \item correctly identifying \pkg{lme4} fits (which use class
#'         \code{lmerMod}, not \code{lmer}).
#' }
#'
#' @note Requires \pkg{ggplot2} plus, depending on the model class,
#'   \pkg{nlme} or \pkg{lme4} (\code{Suggests} dependencies).
#'
#' @seealso \code{\link{plot_trajectories}}
#'
#' @examples
#' # nlme: linear growth in tooth distance over age (Orthodont, 27 children).
#' # Four of the children are paneled here so the figure is quick to draw;
#' # drop n_random to get a panel for every child.
#' fm_nlme <- nlme::lme(distance ~ age, random = ~ age | Subject,
#'                      data = nlme::Orthodont)
#' p <- plot_trajectories_fitted(fm_nlme, n_random = 4, seed = 113)
#' p
#' attr(p, "quality_of_fit")  # per-subject R^2 and RMSE
#'
#' # An lme4 fit is handled the same way. Not run here because the call
#' # loads the lme4 namespace and then draws a panel for each of the
#' # eighteen subjects, which is where the time goes; fitting the model
#' # is quick by comparison. The calls are:
#' # fm_lme4 <- lme4::lmer(Reaction ~ Days + (Days | Subject),
#' #                       data = lme4::sleepstudy)
#' # plot_trajectories_fitted(fm_lme4)
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @keywords hplot
#'
#' @family plotting
#' @family within-subjects analysis
#'
#' @export
plot_trajectories_fitted <- function(model,
                                     id           = NULL,
                                     time         = NULL,
                                     outcome      = NULL,
                                     ids          = NULL,
                                     n_random     = NULL,
                                     pct_random   = NULL,
                                     n_grid       = 100,
                                     show_points  = TRUE,
                                     point_size   = 1.5,
                                     linewidth    = 0.6,
                                     alpha        = 0.8,
                                     palette      = "okabe_ito",
                                     nrow         = NULL,
                                     ncol         = NULL,
                                     show_quality = TRUE,
                                     title        = NULL,
                                     xlab         = NULL,
                                     ylab         = NULL,
                                     seed         = NULL) {
  if (!is.null(seed) &&
      (!is.null(n_random) || !is.null(pct_random))) {
    if (exists(".Random.seed", envir = .GlobalEnv)) {
      .old_seed <- get(".Random.seed", envir = .GlobalEnv)
      on.exit(assign(".Random.seed", .old_seed, envir = .GlobalEnv), add = TRUE)
    } else {
      on.exit(if (exists(".Random.seed", envir = .GlobalEnv)) rm(list = ".Random.seed", envir = .GlobalEnv), add = TRUE)
    }
    set.seed(seed)
  }
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required. Install with install.packages(\"ggplot2\").",
         call. = FALSE)
  }
  if (!is.numeric(n_grid) || n_grid < 10L) stop("'n_grid' must be at least 10.")

  line_color <- .dmar_palette(1, palette = palette)

  info <- .extract_traj_info(model, id, time, outcome, n_grid)
  obs  <- info$observed
  fit  <- info$fitted
  id_c <- info$id_col
  t_c  <- info$time_col
  y_c  <- info$outcome_col

  unique_ids   <- unique(obs[[id_c]])
  selected_ids <- .resolve_ids(unique_ids, ids, n_random, pct_random)
  obs <- obs[obs[[id_c]] %in% selected_ids, , drop = FALSE]
  fit <- fit[fit[[id_c]] %in% selected_ids, , drop = FALSE]

  qof <- .per_id_quality_of_fit(obs, id_c, y_c)

  # Strip labels: optionally include R^2 and RMSE.
  if (isTRUE(show_quality)) {
    label_lookup <- stats::setNames(
      sprintf("%s\nR^2 = %.2f, RMSE = %.2f",
              as.character(qof[[id_c]]),
              qof$r_squared, qof$rmse),
      as.character(qof[[id_c]])
    )
    obs$.facet <- factor(label_lookup[as.character(obs[[id_c]])],
                         levels = unname(label_lookup))
    fit$.facet <- factor(label_lookup[as.character(fit[[id_c]])],
                         levels = unname(label_lookup))
    facet_var <- ".facet"
  } else {
    obs$.facet <- as.character(obs[[id_c]])
    fit$.facet <- as.character(fit[[id_c]])
    facet_var <- ".facet"
  }

  p <- ggplot2::ggplot()

  if (isTRUE(show_points)) {
    p <- p + ggplot2::geom_point(
      data    = obs,
      mapping = ggplot2::aes(x = .data[[t_c]], y = .data[[y_c]]),
      size    = point_size, alpha = alpha
    )
  }
  p <- p + ggplot2::geom_line(
    data    = fit,
    mapping = ggplot2::aes(x = .data[[t_c]], y = .data[[".fitted"]]),
    linewidth = linewidth, color = line_color
  )

  p <- p +
    ggplot2::facet_wrap(stats::reformulate(facet_var),
                        nrow = nrow, ncol = ncol) +
    ggplot2::labs(
      title = title,
      x     = if (is.null(xlab)) t_c else xlab,
      y     = if (is.null(ylab)) y_c else ylab
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(strip.text = ggplot2::element_text(size = 8))

  attr(p, "quality_of_fit") <- qof
  p
}
