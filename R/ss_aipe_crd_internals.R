# Internal helpers shared by R/ss_aipe_crd.R and R/ss_aipe_crd_es.R.
# Not exported. The "*_diff" helpers handle the unstandardized
# conditions difference; the "*_es" helpers handle the standardized
# effect size variants and use a priori Monte Carlo simulation via
# OpenMx for likelihood-based confidence intervals.

# Get a nice report for a public function
.report_crd <- function(n_clusters, cluster_size, width = NULL, cost = NULL, es = FALSE, es_type = 0, assurance = NULL, diff_size = NULL, solved = c("both", "n_clusters", "cluster_size"), requested_width = NULL) {
    # `solved` says which size the planner solved for; only a size that was
    # solved for is reported, in a row named necessary_<size>. A size the
    # user gave is an argument of the call and is not echoed. The cluster
    # size planned is necessary_cluster_size, the name that the size
    # registry of tidy() and glance() reads (R/dmar_tidiers.R).
    solved <- match.arg(solved)
    term <- character(0)
    value <- numeric(0)
    if (solved != "cluster_size") {
        term <- c(term, "necessary_n_clusters")
        value <- c(value, n_clusters)
    }
    if (is.null(diff_size)) {
        if (solved != "n_clusters") {
            term <- c(term, "necessary_cluster_size")
            value <- c(value, cluster_size)
        }
    } else {
        out <- data.frame(table(.find_n_indiv_vec(as.numeric(cluster_size), diff_size, n_clusters)))
        colnames(out) <- c("cluster_size", "freq")
        print(out)
    }
    # The width asked for, where the plan was made from one, is the row
    # `width`, as in the other AIPE planners.
    if (!is.null(requested_width)) {
        term <- c(term, "width")
        value <- c(value, requested_width)
    }
    # One row name for the width whatever the effect: expected_width when no
    # assurance is asked for; with an assurance, assured_width (the width that
    # is not exceeded with that probability) and the assurance in its own row.
    if (!is.null(width)) {
        if (is.null(assurance)) {
            term <- c(term, "expected_width")
            value <- c(value, width)
        } else {
            term <- c(term, "assured_width", "assurance")
            value <- c(value, width, assurance)
        }
    }
    if (!is.null(cost)) {
        term <- c(term, "budget")
        value <- c(value, cost)
    }
    return(data.frame(term, value))
}

.find_n_indiv_vec <- function(cluster_size, diff_size, n_clusters) {
    is_integer <- all(round(diff_size) == diff_size)
    result <- NULL
    if (is_integer) {
        result <- rep(cluster_size + diff_size, length.out = n_clusters)
    } else {
        result <- rep(round(cluster_size * diff_size), length.out = n_clusters)
    }
    result[result < 1] <- 1
    result
}

# Find the number of clusters given the specified width of ES and the cluster size.
# A width of Inf (a design at which no simulated data set gave a usable
# interval, or an assurance that the usable ones cannot reach; see
# .find_width_crd_es()) does not meet the target, so the downward search
# stops there and the upward search goes on. The downward search steps one
# cluster at a time, as before. The upward search brackets the target and
# bisects (.crd_es_search_up()); it has no bound on the number of clusters,
# since a finite width falls toward 0 as clusters are added, so it ends
# once the width is finite; designs with no finite width end it after
# .CRD_ES_MAX_INF_RUN in a row, and the call's budget of designs
# (.CRD_ES_MAX_DESIGNS) bounds the computation. Through the first builds of
# 1.1.0 it stepped one cluster at a time; the first build of this repair
# also stopped it after 100 steps, which refused designs that are reached,
# because the starting value (the number of clusters the unstandardized
# difference needs) leaves out the uncertainty of the standardizer, which
# grows with the square of the effect.
.find_n_clus_crd_es <- function(width, cluster_size, es, es_type = 1, true_icc, pr_treat, true_R2_between = 0, true_R2_within = 0, n_covariates = 0, assurance = NULL, conf_level = 0.95, G = 1000, true_icc_covariate = NULL, seed = NULL, ncpus = 1, diff_size = NULL, tracker = NULL) {
    if (n_covariates > 0 && is.null(true_icc_covariate)) true_icc_covariate <- true_icc
    if (n_covariates > 1) .dmar_refuse_argument(
      "'n_covariates' must be 0 or 1: the planners for the standardized effect size allow at most one covariate beside the treatment indicator.",
      arg = "n_covariates")
    true_total_var <- 1
    if (es_type == 0) {
        true_total_var <- 1
    } else if (es_type == 1) {
        true_total_var <- 1 / (1 - true_icc)
    } else if (es_type == 2) {
        true_total_var <- 1 / true_icc
    } else {
        .dmar_refuse_argument(
          "'es_type' can be 0 (total variance), 1 (level-1 variance), or 2 (level-2 variance) only.",
          arg = "es_type")
    }
    startval <- .find_n_clus_crd_diff(width = width, cluster_size = cluster_size, pr_treat = pr_treat, true_total_var = true_total_var, true_icc = true_icc, true_R2_between = true_R2_between, true_R2_within = true_R2_within, n_covariates = n_covariates, assurance = assurance, conf_level = conf_level, diff_size = diff_size, fn = tracker$fun)
    startval <- as.numeric(startval)
    startwidth <- .find_width_crd_es(G, assurance = assurance, n_clusters = startval, n_treat_clus = round(startval * pr_treat), cluster_size = cluster_size, true_icc = true_icc, es = es, es_type = es_type, true_total_var = true_total_var, covariate = as.logical(n_covariates), true_icc_covariate = true_icc_covariate, true_R2_within = true_R2_within, true_R2_between = true_R2_between, total_var_Z = 1, seed = seed, ncpus = ncpus, conf_level = conf_level, diff_size = diff_size, tracker = tracker)
    # A width meets the target when it is no larger than the target (see
    # .crd_es_search_up()); a start whose width meets it steps down.
    if (startwidth <= width) {
        repeat {
            startval <- startval - 1
            if (round(startval * pr_treat) == 1 || (startval - round(startval * pr_treat)) == 1) return(c(startval + 1, startwidth))
            savedwidth <- startwidth
            startwidth <- .find_width_crd_es(G, assurance = assurance, n_clusters = startval, n_treat_clus = round(startval * pr_treat), cluster_size = cluster_size, true_icc = true_icc, es = es, es_type = es_type, true_total_var = true_total_var, covariate = as.logical(n_covariates), true_icc_covariate = true_icc_covariate, true_R2_within = true_R2_within, true_R2_between = true_R2_between, total_var_Z = 1, seed = seed, ncpus = ncpus, conf_level = conf_level, diff_size = diff_size, tracker = tracker)
            if (startwidth > width) return(c(startval + 1, savedwidth))
        }
    } else {
        evaluate <- function(n) {
            .find_width_crd_es(G, assurance = assurance, n_clusters = n, n_treat_clus = round(n * pr_treat), cluster_size = cluster_size, true_icc = true_icc, es = es, es_type = es_type, true_total_var = true_total_var, covariate = as.logical(n_covariates), true_icc_covariate = true_icc_covariate, true_R2_within = true_R2_within, true_R2_between = true_R2_between, total_var_Z = 1, seed = seed, ncpus = ncpus, conf_level = conf_level, diff_size = diff_size, tracker = tracker)
        }
        return(.crd_es_search_up(evaluate, start = startval, start_width = startwidth, target = width, tracker = tracker, searched = "n_clusters", fixed = cluster_size, assurance = assurance))
    }
}

# The Monte Carlo evaluation of one design, and what is done with a
# simulated data set that gives no usable interval.
#
# .find_width_crd_es() draws G seeds, simulates one data set from each
# (.runrep_width_es_crd()), fits the likelihood-based interval for the
# standardized effect size with OpenMx, and summarizes the G widths. A data
# set whose interval cannot be used is a failed draw, and
# .crd_interval_status() records why (the reasons are listed there). The
# rule for failed draws:
#   * Without an assurance the width is the mean width over the draws with a
#     usable interval: the expected width among studies in which the
#     interval can be formed. Its denominator is the number of usable draws.
#   * With an assurance a failed draw counts as an interval wider than any
#     width: its width is Inf, and the width reported is the `assurance`
#     quantile of all G widths (type 7, as before), which is the width w for
#     which a study gives an interval no wider than w with probability
#     `assurance`. The denominator is G. When too many draws failed for the
#     quantile to fall among the usable widths, it is Inf: no width is
#     reached with that probability at this design.
#   * When no draw gave a usable interval the width is Inf.
# In the searches an Inf width does not meet the target. A public planner
# never reports an Inf width; .crd_es_finish() refuses it with the counts.
# DMAR 1.0.x and the first builds of 1.1.0 took the mean and the quantile
# over the widths that were not NA (na.rm = TRUE) without a word, returned
# NaN when every draw failed, and counted a fit whose optimizer or limit
# search had not converged as a usable interval.
#
# The value is the width, a number, with the attribute "simulation": the
# record of the evaluation (.crd_mc_record()). When `tracker` (an
# environment from .crd_tracker()) is given, the record and the warnings
# the fits gave are added to it, and the public planner reports them once
# at its end; a call that has already simulated .CRD_ES_MAX_DESIGNS designs
# is refused here (.crd_es_refuse_budget()). Without a tracker the failed
# draws and the warnings are signaled here.
.find_width_crd_es <- function(G, n_clusters, n_treat_clus, cluster_size, true_icc, es, es_type = 1, true_total_var = 1, covariate = FALSE, true_icc_covariate = NULL, true_R2_within = NULL, true_R2_between = NULL, total_var_Z = 1, assurance = NULL, seed = NULL, ncpus = 1, conf_level = 0.95, diff_size = NULL, tracker = NULL) {
    if (!is.null(tracker) && length(tracker$records) >= .CRD_ES_MAX_DESIGNS) {
        .crd_es_refuse_budget(tracker, assurance = assurance)
    }
    .dmar_local_seed(seed)

    seed_list <- as.list(sample(1:999999, G))

    # Each replication sets its own seed (seed_list), so the result is the
    # same whatever `ncpus` is.
    draws <- .crd_apply_replications(seed_list, .runrep_width_es_crd, ncpus = ncpus, n_clusters = n_clusters, n_treat_clus = n_treat_clus, cluster_size = cluster_size, true_icc = true_icc, es = es, es_type = es_type, true_total_var = true_total_var, covariate = covariate, true_icc_covariate = true_icc_covariate, true_R2_within = true_R2_within, true_R2_between = true_R2_between, total_var_Z = total_var_Z, conf_level = conf_level, diff_size = diff_size)
    record <- .crd_mc_record(draws, assurance = assurance, conf_level = conf_level, n_clusters = n_clusters, cluster_size = cluster_size)
    warned <- .crd_draw_warnings(draws)
    if (is.null(tracker)) {
        .crd_signal_warnings(warned, total = record$attempted, fun = NULL)
        if (record$failed > 0L) .crd_warn_failed(record, assurance = assurance, fun = NULL)
    } else {
        tracker$records[[length(tracker$records) + 1L]] <- record
        tracker$warnings <- .crd_merge_warnings(tracker$warnings, warned)
    }
    structure(record$width, simulation = record)
}

# .crd_mc_record(draws, assurance, conf_level, n_clusters, cluster_size)
#
# The record of one Monte Carlo evaluation, from the list the replications
# returned (each a list with `width`, `reason`, and `warnings`; see
# .runrep_width_es_crd()). A list:
#   width          the width the evaluation gives (see .find_width_crd_es()):
#                  the mean usable width, or the `assurance` quantile with
#                  failed draws at Inf; Inf when no draw was usable.
#   attempted      G.
#   returned       draws with a usable interval.
#   failed         draws without one (attempted - returned).
#   denominator    the number of draws the width is taken over: `returned`
#                  for the mean, `attempted` for the quantile.
#   failure_reasons the failed draws by reason, a named integer vector
#                  (integer(0) when none failed).
#   failure_messages for a reason that comes with a message (an error from
#                  the fit, "fit_error"), the first such message, a named
#                  character vector (character(0) when there is none), so
#                  that a fit that fails in every data set for one cause,
#                  such as an error in building the model, can be read.
#   width_mc_ci    the Monte Carlo interval of the width at conf_level, which
#                  describes the error of the simulation only: for the
#                  mean, the mean plus and minus the normal quantile times its
#                  standard error (NA with fewer than two usable draws); for
#                  the quantile, the distribution-free interval of two order
#                  statistics of the G widths, positions qbinom((1 - c)/2, G,
#                  a) and qbinom((1 + c)/2, G, a) + 1 (kept within 1 and G).
#   n_clusters, cluster_size  the design.
.crd_mc_record <- function(draws, assurance = NULL, conf_level = 0.95, n_clusters = NA, cluster_size = NA) {
    G <- length(draws)
    widths <- vapply(draws, function(d) as.numeric(d$width), numeric(1))
    reasons <- vapply(draws, function(d) as.character(d$reason), character(1))
    usable <- is.na(reasons) & is.finite(widths)
    reasons[!usable & is.na(reasons)] <- "limit_not_finite"
    n_usable <- sum(usable)
    failure_reasons <- if (any(!usable)) {
        tab <- table(reasons[!usable])
        stats::setNames(as.integer(tab), names(tab))
    } else {
        integer(0)
    }
    messages <- vapply(draws, function(d) {
        m <- d$message
        if (is.null(m) || !length(m)) NA_character_ else as.character(m)[1L]
    }, character(1))
    failure_messages <- character(0)
    for (r in names(failure_reasons)) {
        m <- messages[!usable & reasons == r & !is.na(messages)]
        if (length(m)) failure_messages[[r]] <- m[1L]
    }
    alpha <- 1 - conf_level
    z <- stats::qnorm(1 - alpha / 2)
    if (is.null(assurance)) {
        width <- if (n_usable > 0L) mean(widths[usable]) else Inf
        mc_ci <- if (n_usable >= 2L) {
            width + c(-1, 1) * z * stats::sd(widths[usable]) / sqrt(n_usable)
        } else {
            c(NA_real_, NA_real_)
        }
        denominator <- n_usable
    } else {
        all_widths <- widths
        all_widths[!usable] <- Inf
        width <- if (n_usable > 0L) unname(stats::quantile(all_widths, assurance, type = 7)) else Inf
        sorted <- sort(all_widths)
        lo <- min(max(stats::qbinom(alpha / 2, G, assurance), 1L), G)
        hi <- min(max(stats::qbinom(1 - alpha / 2, G, assurance) + 1L, 1L), G)
        mc_ci <- c(sorted[lo], sorted[hi])
        denominator <- G
    }
    list(width = width, attempted = G, returned = n_usable,
         failed = G - n_usable, denominator = denominator,
         failure_reasons = failure_reasons, failure_messages = failure_messages,
         width_mc_ci = mc_ci, n_clusters = n_clusters, cluster_size = cluster_size)
}

# .crd_apply_replications(X, FUN, ncpus, ..., cluster_type)
#
# lapply(X, FUN, ...) in this session when ncpus is 1, over forked workers
# (parallel::mclapply()) otherwise, and over a socket cluster on Windows,
# which cannot fork (`cluster_type` says which; a test gives "PSOCK" on
# another system). The socket cluster is stopped by on.exit() registered as
# soon as parallel::makeCluster() returns, so that an error in a worker or
# an interrupt does not leave the worker processes and their sockets open
# (through the first builds of 1.1.0 stopCluster() followed the work and
# was skipped on any error). A forked worker that stopped with an error
# returns a "try-error"; the first such error is signaled again here,
# rather than passed on as a result, and a worker that returned nothing is
# refused, since the replications cannot be summarized without it.
.crd_apply_replications <- function(X, FUN, ncpus = 1L, ..., cluster_type = if (.Platform$OS.type == "windows") "PSOCK" else "fork") {
    if (ncpus <= 1L) return(lapply(X, FUN, ...))
    if (identical(cluster_type, "PSOCK")) {
        cl <- parallel::makeCluster(ncpus, type = "PSOCK")
        on.exit(parallel::stopCluster(cl), add = TRUE)
        return(parallel::clusterApplyLB(cl, X, FUN, ...))
    }
    out <- parallel::mclapply(X, FUN, ..., mc.cores = ncpus)
    errored <- vapply(out, inherits, logical(1), what = "try-error")
    if (any(errored)) {
        first <- out[[which(errored)[1L]]]
        condition <- attr(first, "condition")
        if (inherits(condition, "condition")) stop(condition)
        .dmar_refuse(c("A worker process stopped with an error: ",
                       trimws(as.character(first))), "computation_failed")
    }
    if (any(vapply(out, is.null, logical(1)))) {
        .dmar_refuse(c(
            "A worker process ended without returning the results of its ",
            "simulated data sets (it may have run out of memory). Run the ",
            "call again, with a smaller 'ncpus' or with ncpus = 1."),
            "computation_failed")
    }
    out
}

# Create data and find the width of the likelihood-based CI of ES. Returns a
# list: `width` (the width of a usable interval, NA otherwise), `reason` (NA,
# or why the interval cannot be used; see .crd_interval_status()), the
# limits and the estimate, `message` (for reason "fit_error", the error's
# message), and `warnings`, the warnings the
# simulation and the fit gave, as a character vector of messages named by
# the first class of each. The warnings are recorded and muffled here so
# that they reach the user however the replications are run: a forked or
# socket worker does not pass a warning back to this session. The caller
# signals each distinct one once with the number of data sets that gave it
# (.crd_signal_warnings()).
.runrep_width_es_crd <- function(seed, n_clusters, n_treat_clus, cluster_size, true_icc, es, es_type = 1, true_total_var = 1, covariate = FALSE, true_icc_covariate = NULL, true_R2_within = NULL, true_R2_between = NULL, total_var_Z = 1, conf_level = 0.95, diff_size = NULL) {
    warned <- new.env(parent = emptyenv())
    warned$messages <- character(0)
    out <- withCallingHandlers(
        .runrep_width_es_crd_fit(seed = seed, n_clusters = n_clusters, n_treat_clus = n_treat_clus, cluster_size = cluster_size, true_icc = true_icc, es = es, es_type = es_type, true_total_var = true_total_var, covariate = covariate, true_icc_covariate = true_icc_covariate, true_R2_within = true_R2_within, true_R2_between = true_R2_between, total_var_Z = total_var_Z, conf_level = conf_level, diff_size = diff_size),
        warning = function(w) {
            warned$messages <- c(warned$messages, stats::setNames(conditionMessage(w), class(w)[1L]))
            invokeRestart("muffleWarning")
        })
    out$warnings <- warned$messages
    out
}

.runrep_width_es_crd_fit <- function(seed, n_clusters, n_treat_clus, cluster_size, true_icc, es, es_type = 1, true_total_var = 1, covariate = FALSE, true_icc_covariate = NULL, true_R2_within = NULL, true_R2_between = NULL, total_var_Z = 1, conf_level = 0.95, diff_size = NULL) {
    .dmar_local_seed(seed)
    datawide <- .create_data_crd_wide(n_clusters = n_clusters, n_treat_clus = n_treat_clus, cluster_size = cluster_size, true_icc = true_icc, es = es, es_type = es_type, true_total_var = true_total_var, covariate = covariate, true_icc_covariate = true_icc_covariate, true_R2_within = true_R2_within, true_R2_between = true_R2_between, total_var_Z = total_var_Z, diff_size = diff_size)
    y_label <- NULL
    if (!is.null(diff_size)) {
        size <- as.numeric(names(datawide))
        y_label <- lapply(size, function(x) paste("y", 1:x, sep = ""))
    } else {
        y_label <- paste("y", 1:cluster_size, sep = "")
    }

    x_label <- "x"
    z_within_label <- NULL
    if (covariate && true_icc_covariate != 1) {
        if (!is.null(diff_size)) {
            size <- as.numeric(names(datawide))
            z_within_label <- lapply(size, function(x) paste("zw", 1:x, sep = ""))
        } else {
            z_within_label <- paste("zw", 1:cluster_size, sep = "")
        }
    }
    z_between_label <- NULL
    if (covariate && true_icc_covariate != 0) z_between_label <- "zb"
    if (!is.null(diff_size)) {
        screencapture <- utils::capture.output(
            result <- .lik_ci_es_crd_unequal(datawide = datawide, y_label = y_label, x_label = x_label, z_within_label = z_within_label, z_between_label = z_between_label, es_type = es_type, true_icc = true_icc, es = es, true_total_var = true_total_var, covariate = covariate, true_icc_covariate = true_icc_covariate, true_R2_within = true_R2_within, true_R2_between = true_R2_between, total_var_Z = total_var_Z, conf_level = conf_level))
    } else {
        screencapture <- utils::capture.output(
            result <- .lik_ci_es_crd(datawide = datawide, y_label = y_label, x_label = x_label, z_within_label = z_within_label, z_between_label = z_between_label, es_type = es_type, true_icc = true_icc, es = es, true_total_var = true_total_var, covariate = covariate, true_icc_covariate = true_icc_covariate, true_R2_within = true_R2_within, true_R2_between = true_R2_between, total_var_Z = total_var_Z, conf_level = conf_level)
        )
    }
    return(result)
}

# .crd_run_openmx(model)
#
# mxRun() with the likelihood-based intervals, returning the fitted model or,
# when the fit stops with an error, the error condition (a failed draw of
# reason "fit_error"). OpenMx reports an optimizer status other than 0 as a
# warning ("Optimizer returned a non-zero status code"); that status is read
# from the fit and classified by .crd_interval_status(), so this one warning
# is muffled here: a status of 1 is accepted as a usable fit, as OpenMx's
# help describes it, and any other status is a failed draw, counted by
# reason in the record of the simulation. Every other warning passes on to
# .runrep_width_es_crd(), which records it for the user.
.crd_run_openmx <- function(model) {
    tryCatch(
        withCallingHandlers(
            OpenMx::mxRun(model, intervals = TRUE, silent = TRUE),
            warning = function(w) {
                if (grepl("returned a non-zero status code", conditionMessage(w), fixed = TRUE)) {
                    invokeRestart("muffleWarning")
                }
            }),
        error = function(e) e)
}

# .crd_fit_interval(fit)
#
# The interval for the standardized effect size (the algebra "es") of one
# fitted model, classified: reads the optimizer status of the fit
# (fit@output$status$code), the interval (fit@output$confidenceIntervals),
# the optimizer status of each limit's search
# (fit@output$confidenceIntervalCodes), and OpenMx's diagnostic for each
# limit (the `diagnostic` column of the CI details that
# summary(fit, verbose = TRUE) prints, read where summary.MxModel() reads it,
# fit$compute$steps$CI$output$detail), and passes them to
# .crd_interval_status(). `fit` may be the error .crd_run_openmx() returned.
.crd_fit_interval <- function(fit) {
    if (inherits(fit, "condition")) {
        # The error's message is kept, so that a cause common to every data
        # set (an error in building the model, say) can be read from the
        # record and the refusal.
        return(list(width = NA_real_, reason = "fit_error",
                    lower = NA_real_, estimate = NA_real_, upper = NA_real_,
                    message = conditionMessage(fit)))
    }
    output <- fit@output
    ci <- .crd_interval_row(output$confidenceIntervals)
    detail <- tryCatch(fit$compute$steps$CI$output$detail, error = function(e) NULL)
    status <- .crd_interval_status(
        status = output$status$code,
        ci = ci,
        codes = output$confidenceIntervalCodes,
        diagnostic = if (is.data.frame(detail) && "diagnostic" %in% names(detail)) as.character(detail$diagnostic) else NULL)
    status$lower <- if (length(ci) == 3L) ci[[1L]] else NA_real_
    status$estimate <- if (length(ci) == 3L) ci[[2L]] else NA_real_
    status$upper <- if (length(ci) == 3L) ci[[3L]] else NA_real_
    status
}

# The first row of OpenMx's table of intervals as c(lower, estimate,
# upper), read by column name: the table has the columns lbound, estimate,
# and ubound, in that order. Through the first builds of 1.1.0 the width
# was taken as the second entry minus the first, which with this table is
# the estimate minus the lower limit, not the width of the interval. NULL
# when the table is absent or lacks the limits; the estimate is NA when
# the table has no estimate column.
.crd_interval_row <- function(ci) {
    if (!is.matrix(ci) || nrow(ci) < 1L || !all(c("lbound", "ubound") %in% colnames(ci))) {
        return(NULL)
    }
    estimate <- if ("estimate" %in% colnames(ci)) ci[1L, "estimate"] else NA_real_
    as.numeric(c(ci[1L, "lbound"], estimate, ci[1L, "ubound"]))
}

# .crd_interval_status(status, ci, codes, diagnostic = NULL)
#
# Whether the likelihood-based interval of one simulated data set can be
# used, and if not why. A list with `width` (upper - lower, or NA) and
# `reason` (NA for a usable interval). The checks, in this order, and the
# reason recorded when one fails:
#   status      the optimizer status code of the fit must be 0 or 1. OpenMx
#               documents (?as.statusCode) 0 as "OK" and 1 as "OK/green":
#               optimization succeeded, with 1 (detected by NPSOL only) a
#               solution that "is likely okay", so a fit of status 1 is
#               used. Every other code (2 and 3 infeasible constraints, 4
#               the iteration limit, 5 a Hessian that does not appear to be
#               convex, 6 the first-order optimality conditions not met, 7
#               bad derivatives, 9 an internal error, 10 infeasible starting
#               values) means the fit may not be at the maximum the interval
#               is measured from; for 5 and 6 OpenMx's help says it would
#               not trust the solution. Reason "optimizer_status_<code>".
#   ci          the row of the interval, c(lower, estimate, upper); absent
#               gives "limit_not_found".
#   diagnostic  OpenMx's diagnostic for each limit, when available, must be
#               "success". OpenMx's help pages do not list these values; the
#               others the installed version contains include "alpha level
#               not reached" and "active box constraint", and none of them
#               says that the limit was found. Reason "limit_not_found".
#   codes       the optimizer status code of each limit's search must be 0
#               or 1, as for the fit; NA (no code) gives "limit_not_found",
#               another code "limit_status_<code>".
#   limits      both must be finite ("limit_not_finite") and bracket the
#               estimate, lower <= estimate <= upper with lower < upper
#               ("limits_not_ordered"; an NA estimate is not compared).
# The width is the upper limit minus the lower.
.crd_interval_status <- function(status, ci, codes, diagnostic = NULL) {
    fail <- function(reason) list(width = NA_real_, reason = reason)
    status <- if (length(status) == 1L && is.numeric(status)) status else NA_real_
    if (is.na(status) || !(status %in% c(0, 1))) {
        return(fail(paste0("optimizer_status_", if (is.na(status)) "NA" else status)))
    }
    if (length(ci) != 3L) return(fail("limit_not_found"))
    if (!is.null(diagnostic) &&
        (length(diagnostic) != 2L || anyNA(diagnostic) || any(diagnostic != "success"))) {
        return(fail("limit_not_found"))
    }
    codes <- if (is.numeric(codes)) as.vector(codes) else numeric(0)
    if (length(codes) != 2L || anyNA(codes)) return(fail("limit_not_found"))
    bad <- codes[!(codes %in% c(0, 1))]
    if (length(bad)) return(fail(paste0("limit_status_", bad[1L])))
    if (!all(is.finite(ci[c(1L, 3L)]))) return(fail("limit_not_finite"))
    brackets <- is.na(ci[2L]) || (ci[1L] <= ci[2L] && ci[2L] <= ci[3L])
    if (!(ci[1L] < ci[3L] && brackets)) return(fail("limits_not_ordered"))
    list(width = ci[3L] - ci[1L], reason = NA_character_)
}

# The reason of a failed draw (.crd_interval_status()) in words, with the
# label ?as.statusCode gives a status code.
.crd_reason_text <- function(reason) {
    status_label <- function(code) {
        label <- if (requireNamespace("OpenMx", quietly = TRUE)) {
            tryCatch(as.character(OpenMx::as.statusCode(code)), error = function(e) NA_character_)
        } else {
            NA_character_
        }
        if (length(label) == 1L && !is.na(label)) paste0(" (", label, ")") else ""
    }
    vapply(reason, function(r) {
        if (r == "fit_error") return("the fit stopped with an error")
        if (r == "limit_not_found") return("a confidence limit was not found")
        if (r == "limit_not_finite") return("a confidence limit was not finite")
        if (r == "limits_not_ordered") return("the limits did not bracket the estimate")
        code <- sub("^(optimizer|limit)_status_", "", r)
        if (startsWith(r, "optimizer_status_")) {
            return(paste0("the fit ended with optimizer status ", code, status_label(code)))
        }
        if (startsWith(r, "limit_status_")) {
            return(paste0("the search for a confidence limit ended with optimizer status ",
                          code, status_label(code)))
        }
        r
    }, character(1), USE.NAMES = FALSE)
}

# The failed draws by reason, as "3 (a confidence limit was not found)" joined
# by commas. With `messages` (the failure_messages of a record), the first
# error message of a reason follows its count.
.crd_reasons_phrase <- function(failure_reasons, messages = character(0)) {
    if (!length(failure_reasons)) return("")
    reason <- names(failure_reasons)
    text <- .crd_reason_text(reason)
    has <- reason %in% names(messages)
    text[has] <- paste0(text[has], "; the first error was: ",
                        gsub("\\s+", " ", trimws(messages[reason[has]])))
    paste0(failure_reasons, " (", text, ")", collapse = ", ")
}

# The state of one call of a public standardized planner: the record of
# every design it simulated and the warnings the fits gave, kept in an
# environment so that the evaluations deep in a search can add to it.
.crd_tracker <- function(fun) {
    st <- new.env(parent = emptyenv())
    st$fun <- fun
    st$records <- list()
    st$warnings <- list()
    st
}

# The totals over a list of records (.crd_mc_record()): `n_designs`, the
# data sets `attempted`, `returned` (with a usable interval), and `failed`,
# `failure_reasons` (the failed data sets by reason, summed over the
# records, by name), and `failure_messages` (the first error message of
# each reason that has one).
.crd_search_totals <- function(records) {
    records <- Filter(is.list, records)
    count <- function(what) sum(vapply(records, function(r) {
        v <- r[[what]]
        if (is.null(v)) 0L else as.integer(v)
    }, integer(1)))
    reasons <- unlist(lapply(records, function(r) r$failure_reasons))
    failure_reasons <- if (length(reasons)) {
        tab <- tapply(reasons, names(reasons), sum)
        stats::setNames(as.integer(tab), names(tab))
    } else {
        integer(0)
    }
    messages <- character(0)
    for (r in records) {
        for (k in names(r$failure_messages)) {
            if (!k %in% names(messages)) messages[[k]] <- r$failure_messages[[k]]
        }
    }
    list(n_designs = length(records), attempted = count("attempted"),
         returned = count("returned"), failed = count("failed"),
         failure_reasons = failure_reasons, failure_messages = messages)
}

# A warning's message with its numbers replaced by "#". Warnings of one
# class whose messages differ only in their numbers (a value or a count
# that changes from data set to data set) are one kind of warning, and are
# given together.
.crd_warning_pattern <- function(message) {
    gsub("[-+]?([0-9]+[.]?[0-9]*|[.][0-9]+)([eE][-+]?[0-9]+)?", "#", message)
}

# The warnings the replications of one evaluation recorded, grouped by
# class and by message with its numbers set aside (.crd_warning_pattern()):
# a list of list(class, pattern, message, messages, n), with `message` the
# first message of the group, `messages` its distinct messages, and `n` the
# number of data sets that gave a warning of the group (a data set counts
# once in a group however many such warnings it gave).
.crd_draw_warnings <- function(draws) {
    groups <- list()
    for (d in draws) {
        w <- d$warnings
        if (!length(w)) next
        cls <- names(w)
        if (is.null(cls)) cls <- rep("simpleWarning", length(w))
        msg <- unname(as.character(w))
        pattern <- .crd_warning_pattern(msg)
        key <- paste(cls, pattern, sep = "\r")
        for (i in seq_along(key)) {
            k <- key[i]
            if (is.null(groups[[k]])) {
                groups[[k]] <- list(class = cls[i], pattern = pattern[i],
                                    message = msg[i], messages = character(0), n = 0L)
            }
            if (!msg[i] %in% groups[[k]]$messages) {
                groups[[k]]$messages <- c(groups[[k]]$messages, msg[i])
            }
        }
        for (k in unique(key)) groups[[k]]$n <- groups[[k]]$n + 1L
    }
    unname(groups)
}

.crd_merge_warnings <- function(known, new) {
    for (w in new) {
        hit <- which(vapply(known, function(k) identical(k$class, w$class) &&
                                identical(k$pattern, w$pattern), logical(1)))
        if (length(hit)) {
            h <- hit[1L]
            known[[h]]$n <- known[[h]]$n + w$n
            known[[h]]$messages <- union(known[[h]]$messages, w$messages)
        } else {
            known[[length(known) + 1L]] <- w
        }
    }
    known
}

# Each kind of warning the fits gave (.crd_draw_warnings()), once, with the
# number of simulated data sets that gave it out of `total`: class
# dmar_result_qualified, with the original class, that count, and the
# distinct messages as `warning_class`, `n_replications`, and `messages`.
.crd_signal_warnings <- function(warned, total, fun = NULL) {
    lead <- if (is.null(fun)) "" else paste0(fun, "(): ")
    for (w in warned) {
        k <- length(w$messages)
        more <- if (k > 1L) {
            c(" (", k, " messages of this form, which differ in their numbers; ",
              "the first is shown, and all are in the field 'messages' of this ",
              "warning)")
        } else {
            character(0)
        }
        .dmar_warn(c(
            lead, w$n, " of the ", total, " simulated data sets gave this ",
            "warning (", w$class, ") when the data were simulated or the model ",
            "was fitted: ", w$message, more),
            "result_qualified", warning_class = w$class, n_replications = w$n,
            messages = w$messages)
    }
    invisible(NULL)
}

# Signals the warnings a tracker holds, once; called on exit from a public
# planner, so that they are given whether it returns or refuses.
.crd_forward_warnings <- function(tracker) {
    total <- sum(vapply(tracker$records, function(r) as.numeric(r$attempted), numeric(1)))
    .crd_signal_warnings(tracker$warnings, total = total, fun = tracker$fun)
    tracker$warnings <- list()
    invisible(NULL)
}

# The warning that some simulated data sets at the reported design gave no
# usable interval: class dmar_replications_without_interval, with the
# counts of the record and `scope` "design".
.crd_warn_failed <- function(record, assurance = NULL, fun = NULL, search = NULL) {
    lead <- if (is.null(fun)) "" else paste0(fun, "(): ")
    how <- if (is.null(assurance)) {
        c("The expected width is the mean over the ", record$returned,
          " that gave one, so it describes studies in which the interval ",
          "could be formed.")
    } else {
        c("The assured width counts a data set without an interval as one ",
          "whose interval is wider than any width, so the assurance of ",
          format(assurance), " is the probability, over all ", record$attempted,
          ", that a study gives an interval no wider than the assured width.")
    }
    more <- if (!is.null(search) && search$n_designs > 1L && search$failed > record$failed) {
        c(" Over the search, ", search$failed, " of the ", search$attempted,
          " data sets simulated at ", search$n_designs, " designs gave no ",
          "usable interval.")
    } else {
        character(0)
    }
    .dmar_warn(c(
        lead, record$failed, " of the ", record$attempted, " data sets simulated ",
        "at the design reported (", record$n_clusters, " clusters of ",
        record$cluster_size, ") gave no usable confidence interval: ",
        .crd_reasons_phrase(record$failure_reasons, record$failure_messages),
        ". ", how, more, " The counts are in attr(x, \"simulation\")."),
        "replications_without_interval", scope = "design",
        n_without_interval = record$failed, G = record$attempted,
        attempted = record$attempted, returned = record$returned,
        failed = record$failed, denominator = record$denominator,
        failure_reasons = record$failure_reasons,
        failure_messages = record$failure_messages)
}

# The warning that every data set simulated at the reported design gave a
# usable interval but some at other designs the search simulated did not:
# those data sets bore on where the search stopped. Class
# dmar_replications_without_interval, with the totals of the search and
# `scope` "search".
.crd_warn_failed_search <- function(record, search, assurance = NULL, fun = NULL) {
    lead <- if (is.null(fun)) "" else paste0(fun, "(): ")
    how <- if (is.null(assurance)) {
        c("A design at which no data set gave an interval counted as one whose ",
          "width does not meet the target, and the expected width of any other ",
          "design is the mean over the data sets that gave one, so these data ",
          "sets bear on where the search stopped.")
    } else {
        c("Each counted as a data set whose interval is wider than any width, ",
          "so these data sets bear on where the search stopped.")
    }
    .dmar_warn(c(
        lead, "each of the ", record$attempted, " data sets simulated at the ",
        "design reported (", record$n_clusters, " clusters of ", record$cluster_size,
        ") gave a usable confidence interval, but ", search$failed, " of the ",
        search$attempted, " data sets simulated at the ", search$n_designs,
        " designs the search simulated did not: ",
        .crd_reasons_phrase(search$failure_reasons, search$failure_messages),
        ". ", how, " The counts are in attr(x, \"simulation\")$search."),
        "replications_without_interval", scope = "search",
        n_without_interval = search$failed, G = record$attempted,
        attempted = search$attempted, returned = search$returned,
        failed = search$failed, failure_reasons = search$failure_reasons,
        failure_messages = search$failure_messages)
}

# .crd_es_finish(out, tracker, n_clusters, cluster_size, assurance)
#
# The end of a public standardized planner that simulated: finds the record
# of the design whose width the table reports (the last evaluation of that
# number of clusters and cluster size), refuses a width that is not finite
# (no usable interval at that design, or an assurance the usable ones
# cannot reach) with an error of class dmar_computation_failed that gives
# the counts, attaches the record as attr(out, "simulation") with the
# totals of the whole call and the warnings the fits gave in its element
# `search`, and warns when draws failed (at the design reported, or else
# at other designs of the search). A table from a call that simulated
# nothing is returned as it is.
.crd_es_finish <- function(out, tracker, n_clusters, cluster_size, assurance = NULL) {
    records <- tracker$records
    if (!length(records)) return(out)
    hit <- which(vapply(records, function(r) isTRUE(r$n_clusters == n_clusters) &&
                            isTRUE(r$cluster_size == cluster_size), logical(1)))
    if (!length(hit)) {
        .dmar_refuse(c(tracker$fun, "(): the simulation of the design reported ",
                       "was not found among the designs simulated."), "internal_error")
    }
    record <- records[[hit[length(hit)]]]
    search <- .crd_search_totals(records)
    search$warnings <- tracker$warnings
    if (record$returned == 0L || !is.finite(record$width)) {
        kind <- if (is.null(assurance)) "expected" else "assured"
        why <- if (record$returned == 0L) {
            c("none of the ", record$attempted, " data sets simulated at the ",
              "design gave a usable confidence interval, so no ", kind,
              " width can be computed")
        } else {
            c(record$failed, " of the ", record$attempted, " data sets simulated ",
              "at the design gave no usable confidence interval, more than the ",
              "share 1 - assurance = ", format(1 - assurance), " allows, so ",
              "no width is reached with probability ", format(assurance),
              " (a data set without an interval counts as one whose interval is ",
              "wider than any width)")
        }
        .dmar_refuse(c(
            tracker$fun, "(): at ", record$n_clusters, " clusters of ",
            record$cluster_size, ", ", why, ". The failed data sets: ",
            .crd_reasons_phrase(record$failure_reasons, record$failure_messages),
            ". A design with more clusters, larger clusters, or (for es_type = 2) ",
            "a larger intraclass correlation gives an interval more often."),
            "computation_failed", attempted = record$attempted,
            returned = record$returned, failed = record$failed,
            failure_reasons = record$failure_reasons,
            failure_messages = record$failure_messages, search = search)
    }
    record$search <- search
    attr(out, "simulation") <- record
    if (record$failed > 0L) {
        .crd_warn_failed(record, assurance = assurance, fun = tracker$fun, search = search)
    } else if (search$failed > 0L) {
        .crd_warn_failed_search(record, search, assurance = assurance, fun = tracker$fun)
    }
    out
}

# Bounds on the computation of the standardized planners. Each design a
# search simulates is G data sets, each fitted with OpenMx, so these bound
# the time a call takes; none of them is a statistical rule.
#   .CRD_ES_MAX_DESIGNS        the most designs one call of a public planner
#                              simulates, counted in its tracker; the two
#                              planners that search over both sizes, which
#                              run a search for the number of clusters at
#                              each cluster size they try, share it among
#                              those searches.
#   .CRD_ES_CLUSTER_SIZE_REACH how far the cluster size search goes above
#                              the cluster size the unstandardized
#                              difference needs, or above its own start
#                              when that size is beyond 100,000
#                              (.find_n_indiv_crd_es()). The OpenMx model
#                              has one observed variable per individual of
#                              a cluster, and a fit takes far longer as the
#                              clusters grow: with 20 clusters, one fit took
#                              about 0.4 s with 100 individuals per cluster,
#                              1.7 s with 200, and 11 s with 400 (OpenMx
#                              2.22.11, one core), against 0.03 s to 0.06 s
#                              with 5 individuals in 200 to 4000 clusters.
#   .CRD_ES_MAX_INF_RUN        the most designs in a row with no finite width
#                              an upward search simulates before it stops.
#   .CRD_ES_MAX_CLUSTER_SIZE   the largest cluster size the cluster size
#                              search simulates, its start included. At
#                              the fit times above, a design of G = 1000
#                              data sets at 500 individuals per cluster
#                              takes hours, and the cost grows faster than
#                              the cluster size; a search that would start
#                              above it is refused before anything is
#                              simulated (reason "cluster_size_beyond_reach").
#   .CRD_ES_MAX_G              the largest G: each simulated data set takes
#                              its own seed, drawn without replacement from
#                              1 to 999,999 (.find_width_crd_es()), so no
#                              more than 999,999 can be drawn. The public
#                              planners refuse a larger G on entry, with
#                              class dmar_count_invalid; through the first
#                              builds of 1.1.0 it passed the check and the
#                              draw stopped with "cannot take a sample
#                              larger than the population when 'replace =
#                              FALSE'". The draw is unchanged, so every G
#                              that was accepted gives the seeds it gave.
.CRD_ES_MAX_DESIGNS <- 500L
.CRD_ES_CLUSTER_SIZE_REACH <- 100L
.CRD_ES_MAX_INF_RUN <- 3L
.CRD_ES_MAX_CLUSTER_SIZE <- 500L
.CRD_ES_MAX_G <- 999999L
.CRD_ES_G_REASON <- paste0(
    "Each simulated data set takes its own seed, drawn without replacement ",
    "from the whole numbers 1 to 999,999, so at most 999,999 data sets can ",
    "be simulated.")

.crd_es_fun <- function(tracker) {
    if (is.null(tracker)) "The planner" else paste0(tracker$fun, "()")
}

# .crd_es_search_up(evaluate, start, start_width, target, tracker, searched,
#                   fixed, max_size, assurance, width_limit, reach_from)
#
# The upward search of the standardized planners: from `start`, whose width
# `start_width` does not meet `target`, the smallest larger size whose
# simulated width meets it. A width meets the target when it is finite and
# no larger than the target, here and in the downward searches of
# .find_n_clus_crd_es() and .find_n_indiv_crd_es() alike (through the first
# builds of 1.1.0 this search alone asked for a width below the target, so
# a width equal to it was met going down and not going up; with simulated
# widths a tie has probability 0). evaluate(size) is the simulated
# width at a size, with the record of its simulation as the attribute
# "simulation" (.find_width_crd_es()). The search first brackets the
# target: the step from the last size doubles after each size whose width
# does not meet the target (start + 1, start + 3, start + 7, ...), and after
# a size with no finite width the next size is at least twice that size.
# It then bisects the bracket, the last size that did not meet the target
# and the first that did, to adjacent sizes, and returns c(size, width).
# When the simulated width falls as the size grows, this is the smallest
# size above the start that meets the target, the design that stepping up
# one size at a time (as the searches did through the first builds of
# 1.1.0) gives, with far fewer designs simulated: about twice the base 2
# logarithm of the distance from the start. When the error of the
# simulation makes the width rise and fall near the target, the two can
# differ by a size or two. A size with no finite width counts as one that
# does not meet the target, so when such sizes lie among sizes with a
# finite width (failures near the share an assurance allows, or an
# artificial pattern such as review1/r01's, where every odd number of
# clusters failed: one-step 38, this search 62), the size returned meets
# the target and the one below it does not, but a smaller size may also
# meet it. A bracket of adjacent sizes is not bisected, so a target one
# size above the start takes one design, as before.
#
# The search ends with an error of class dmar_computation_failed:
#   * after .CRD_ES_MAX_INF_RUN designs in a row with no finite width
#     (.crd_es_refuse_no_finite(), reason "no_finite_width");
#   * when it reaches `max_size` without meeting the target (the cluster
#     size search; .crd_es_refuse_unreached(), reason "width_not_reached").
# The tracker's budget of designs (.CRD_ES_MAX_DESIGNS) is checked by
# .find_width_crd_es(). `width_limit` is the width the cluster size search
# approaches as the clusters grow, where it is known
# (.crd_es_width_limit()), and `reach_from` the size `max_size` was set
# from (.CRD_ES_CLUSTER_SIZE_REACH above it), both for the message.
.crd_es_search_up <- function(evaluate, start, start_width, target, tracker = NULL, searched = c("n_clusters", "cluster_size"), fixed = NULL, max_size = Inf, assurance = NULL, width_limit = NA_real_, reach_from = start) {
    searched <- match.arg(searched)
    lo <- start
    lo_width <- start_width
    best <- list(size = start, width = start_width)
    inf_run <- if (is.finite(start_width)) list() else list(list(size = start, width = start_width))
    step <- 1
    repeat {
        if (length(inf_run) >= .CRD_ES_MAX_INF_RUN) {
            .crd_es_refuse_no_finite(tracker, searched = searched, run = inf_run,
                                     assurance = assurance, fixed = fixed)
        }
        if (lo >= max_size) {
            .crd_es_refuse_unreached(tracker, from = start, to = lo, best = best,
                                     target = target, assurance = assurance,
                                     fixed = fixed, width_limit = width_limit,
                                     reach_from = reach_from)
        }
        nxt <- if (is.finite(lo_width)) lo + step else max(lo + step, 2 * lo)
        nxt <- min(nxt, max_size)
        w <- evaluate(nxt)
        if (is.finite(w) && w <= target) {
            hi <- nxt
            hi_width <- w
            break
        }
        if (is.finite(w)) {
            inf_run <- list()
            if (!is.finite(best$width) || w < best$width) best <- list(size = nxt, width = w)
        } else {
            inf_run[[length(inf_run) + 1L]] <- list(size = nxt, width = w)
        }
        step <- 2 * (nxt - lo)
        lo <- nxt
        lo_width <- w
    }
    while (hi - lo > 1) {
        mid <- floor((lo + hi) / 2)
        w <- evaluate(mid)
        if (is.finite(w) && w <= target) {
            hi <- mid
            hi_width <- w
        } else {
            lo <- mid
        }
    }
    c(hi, hi_width)
}

# The refusal of an upward search after .CRD_ES_MAX_INF_RUN designs in a
# row with no finite width (`run`, a list of list(size, width), the width
# with its record): class dmar_computation_failed, reason
# "no_finite_width", with the sizes, the counts and reasons over those
# designs, and the totals of the call (`search`).
.crd_es_refuse_no_finite <- function(tracker, searched, run, assurance = NULL, fixed = NULL) {
    fun <- .crd_es_fun(tracker)
    kind <- if (is.null(assurance)) "expected" else "assured"
    sizes <- vapply(run, function(r) as.numeric(r$size), numeric(1))
    sims <- lapply(run, function(r) attr(r$width, "simulation"))
    known <- vapply(sims, is.list, logical(1))
    totals <- .crd_search_totals(sims[known])
    none_usable <- all(known) && totals$returned == 0L
    where <- if (searched == "cluster_size") {
        c("cluster sizes of ", paste(sizes, collapse = ", "), " with ", fixed, " clusters")
    } else {
        c(paste(sizes, collapse = ", "), " clusters of ", fixed)
    }
    why <- if (none_usable) {
        "no simulated data set gave a usable confidence interval"
    } else if (!is.null(assurance)) {
        c("more of the simulated data sets gave no usable confidence interval ",
          "than the share 1 - assurance = ", format(1 - assurance), " allows, so ",
          "no width is reached with probability ", format(assurance))
    } else {
        "no simulated data set gave a usable confidence interval"
    }
    counts <- if (all(known)) {
        c(" Data sets without a usable interval: ",
          paste0(vapply(sims, function(s) as.integer(s$failed), integer(1)), " of ",
                 vapply(sims, function(s) as.integer(s$attempted), integer(1)),
                 collapse = ", "),
          "; by reason, ",
          .crd_reasons_phrase(totals$failure_reasons, totals$failure_messages), ".")
    } else {
        character(0)
    }
    advice <- if (!none_usable && !is.null(assurance)) {
        c(" A smaller assurance, or assurance = NULL (the expected width over ",
          "the data sets that give an interval), can be planned.")
    } else {
        character(0)
    }
    .dmar_refuse(c(
        fun, ": the search found no finite ", kind, " width at ", length(run),
        " designs in a row (", where, "): at each, ", why, ".", counts,
        " The search stops here rather than simulate larger designs.", advice),
        "computation_failed", reason = "no_finite_width", searched = searched,
        sizes = sizes, attempted = totals$attempted, returned = totals$returned,
        failed = totals$failed, failure_reasons = totals$failure_reasons,
        failure_messages = totals$failure_messages,
        search = if (is.null(tracker)) list() else .crd_search_totals(tracker$records))
}

# The refusal of the cluster size search when it reaches the largest
# cluster size it simulates (`to`, .CRD_ES_CLUSTER_SIZE_REACH above
# `reach_from`, which is the start `from` or the larger cluster size the
# unstandardized difference needs, and at most .CRD_ES_MAX_CLUSTER_SIZE)
# with the width still above the target: class
# dmar_computation_failed, reason "width_not_reached", with the narrowest
# width found, its Monte Carlo interval, its counts, and the totals of the
# call. The message says what the reach of the search, not the target,
# decided, and what the limit of the width says when it is known: a target
# above it is reached by a larger cluster size. A larger G is suggested
# only when the target lies within the Monte Carlo interval of the
# narrowest width, where it can tell whether that design meets the target.
.crd_es_refuse_unreached <- function(tracker, from, to, best, target, assurance = NULL, fixed = NULL, width_limit = NA_real_, reach_from = from) {
    fun <- .crd_es_fun(tracker)
    kind <- if (is.null(assurance)) "expected" else "assured"
    sim <- attr(best$width, "simulation")
    mc <- if (is.list(sim)) as.numeric(sim$width_mc_ci) else c(NA_real_, NA_real_)
    narrowest <- as.numeric(best$width)
    narrow_text <- if (is.finite(narrowest)) {
        c(" The narrowest was ", format(signif(narrowest, 4)),
          if (all(is.finite(mc))) c(" (Monte Carlo interval ", format(signif(mc[1L], 4)),
                                    " to ", format(signif(mc[2L], 4)), ")") else "",
          " at a cluster size of ", best$size, ".")
    } else {
        " No cluster size searched gave a finite width."
    }
    limit_text <- if (is.finite(width_limit)) {
        c(" As the clusters grow, the ", kind, " width approaches ",
          format(signif(width_limit, 4)), ", its value with every cluster mean ",
          "known exactly; the target is above that, so a cluster size larger than ",
          "this search simulates reaches it.")
    } else {
        c(" With a fixed number of clusters the width levels off as the clusters ",
          "grow, because the clusters themselves vary; whether a larger cluster ",
          "size reaches the target was not determined.")
    }
    g_text <- if (all(is.finite(mc)) && mc[1L] <= target && target <= mc[2L]) {
        c(" The target lies within the Monte Carlo interval of the narrowest ",
          "width, so a larger G can also tell whether that design meets it.")
    } else {
        character(0)
    }
    .dmar_refuse(c(
        fun, ": the ", kind, " width stayed above the target of ", format(target),
        " at cluster sizes from ", from, " to ", to, " with ", fixed, " clusters, ",
        "the largest the search simulates (", .CRD_ES_CLUSTER_SIZE_REACH, " above ",
        if (isTRUE(reach_from > from)) {
            c(format(reach_from, big.mark = ",", scientific = FALSE, trim = TRUE),
              ", the cluster size the unstandardized difference ",
              "between the condition means needs,")
        } else {
            "its starting value"
        },
        " and at most ", .CRD_ES_MAX_CLUSTER_SIZE,
        ", since a fit takes far longer as the clusters grow).",
        narrow_text, limit_text, g_text,
        " ss_aipe_crd_es_n_clusters_fixed_width() plans the number of clusters ",
        "instead."),
        "computation_failed", reason = "width_not_reached", searched = "cluster_size",
        from = from, to = to, target = target, narrowest = narrowest,
        narrowest_at = best$size, narrowest_mc_ci = mc, width_limit = width_limit,
        attempted = if (is.list(sim)) sim$attempted else NA_integer_,
        returned = if (is.list(sim)) sim$returned else NA_integer_,
        failed = if (is.list(sim)) sim$failed else NA_integer_,
        failure_reasons = if (is.list(sim)) sim$failure_reasons else integer(0),
        search = if (is.null(tracker)) list() else .crd_search_totals(tracker$records))
}

# The refusal of a call that has simulated .CRD_ES_MAX_DESIGNS designs:
# class dmar_computation_failed, reason "search_budget", with the totals of
# the call. The message says that the bound on the computation, not a
# finding about the target, ended the call.
.crd_es_refuse_budget <- function(tracker, assurance = NULL) {
    fun <- .crd_es_fun(tracker)
    kind <- if (is.null(assurance)) "expected" else "assured"
    totals <- .crd_search_totals(tracker$records)
    last <- tracker$records[[length(tracker$records)]]
    failed_text <- if (totals$failed > 0L) {
        c(": ", .crd_reasons_phrase(totals$failure_reasons, totals$failure_messages))
    } else {
        character(0)
    }
    .dmar_refuse(c(
        fun, ": the call simulated ", totals$n_designs, " designs, the most one ",
        "call simulates, without settling on a design; that bound on the ",
        "computation ended the call, not a finding about the target. The last ",
        "design simulated was ", last$n_clusters, " clusters of ",
        last$cluster_size, ", with an ", kind, " width of ",
        format(signif(as.numeric(last$width), 4)), ". Of the ", totals$attempted,
        " data sets simulated, ", totals$failed, " gave no usable interval",
        failed_text, "."),
        "computation_failed", reason = "search_budget", n_designs = totals$n_designs,
        attempted = totals$attempted, returned = totals$returned,
        failed = totals$failed, failure_reasons = totals$failure_reasons,
        failure_messages = totals$failure_messages)
}

# .crd_es_cluster_size_start(width, width_limit, true_icc, true_R2_within)
#
# A starting cluster size for the search of es_type 1 where the limit is
# known: the smallest m at which width_limit * sqrt(1 + v / (tau m)) falls
# below `width`, with v / tau = (1 - true_R2_within) (1 - icc) / icc, the
# within-cluster variance left by a within-cluster covariate over the
# between-cluster variance (the caller passes true_R2_within = 0 when there
# is no covariate). That is the limit's form with the variance of a
# cluster mean at m, tau + v / m, in place of tau, and the exact expected
# width (or quantile) with the within-cluster standard deviation known; it
# leaves out the uncertainty of the standardizer, so the width simulated
# at it is larger and the search goes on upward from there. With no
# covariate, equal cluster sizes, and the expected width, it is a lower
# bound on the cluster size the target needs (.find_n_indiv_crd_es() says
# why). At least 2.
.crd_es_cluster_size_start <- function(width, width_limit, true_icc, true_R2_within = 0) {
    ratio <- (1 - true_R2_within) * (1 - true_icc) / true_icc
    excess <- (width / width_limit)^2 - 1
    max(2, floor(ratio / excess) + 1)
}

# .crd_es_width_limit(n_clusters, n_treat_clus, true_icc, conf_level,
#                     assurance)
#
# The width that the likelihood-based interval for the effect standardized
# by the within-cluster standard deviation (es_type = 1) approaches as every
# cluster grows without bound at a fixed number of clusters K, with no
# covariate at the cluster level: the expected width, or with an assurance
# its `assurance` quantile. NA when K < 3 or a condition has no cluster.
#
# In that limit the within-cluster standard deviation sigma_w is known and
# the cluster means are exact, so the data are K means with a common
# unknown variance tau, of which K1 are treated (S_xx = K1 K0 / K). The
# profile likelihood of the difference d, with the intercept and tau at
# their maximum likelihood values, falls by K log(1 + (d - dhat)^2 S_xx /
# RSS), RSS the residual sum of squares, so the interval at confidence c is
# |d - dhat| <= sqrt(RSS (exp(q / K) - 1) / S_xx), q = qchisq(c, 1). Its
# width on the standardized scale, with RSS = tau chi^2_{K - 2}, is
#   2 sqrt((exp(q / K) - 1) / S_xx) sqrt(tau / sigma_w^2) sqrt(chi^2_{K - 2}),
# where tau / sigma_w^2 = icc / (1 - icc). Its mean takes
# E sqrt(chi^2_nu) = sqrt(2) Gamma((nu + 1) / 2) / Gamma(nu / 2), and its
# `assurance` quantile takes sqrt(qchisq(assurance, nu)), nu = K - 2. The
# widths simulated at a finite cluster size lie above this value and fall
# toward it as the clusters grow (evidence/crd/verify_after_round2.log of
# the November 2026 review).
.crd_es_width_limit <- function(n_clusters, n_treat_clus, true_icc, conf_level = 0.95, assurance = NULL) {
    K <- n_clusters
    K1 <- n_treat_clus
    K0 <- K - K1
    if (!(K >= 3 && K1 >= 1 && K0 >= 1 && true_icc > 0 && true_icc < 1)) return(NA_real_)
    nu <- K - 2
    q <- stats::qchisq(conf_level, 1)
    s_xx <- K1 * K0 / K
    root_chisq <- if (is.null(assurance)) {
        sqrt(2) * exp(lgamma((nu + 1) / 2) - lgamma(nu / 2))
    } else {
        sqrt(stats::qchisq(assurance, nu))
    }
    2 * sqrt((exp(q / K) - 1) / s_xx) * sqrt(true_icc / (1 - true_icc)) * root_chisq
}

# Calculate a likelihood-based CI of ES
.lik_ci_es_crd <- function(datawide, y_label, x_label, z_within_label = NULL, z_between_label = NULL, es_type = 1, true_icc = 0.25, es = 0.5, true_total_var = 1, covariate = FALSE, true_icc_covariate = 0.25, true_R2_within = 0.5, true_R2_between = 0.5, total_var_Z = 1, conf_level = 0.95) {
    tau <- true_icc * true_total_var
    sigma <- (1 - true_icc) * true_total_var
    gamma1 <- NULL
    if (es_type == 0) {
        gamma1 <- es * sqrt(true_total_var)
    } else if (es_type == 1) {
        gamma1 <- es * sqrt(sigma)
    } else if (es_type == 2) {
        gamma1 <- es * sqrt(tau)
    } else {
        .dmar_refuse_argument(
          "'es_type' can be 0 (total variance), 1 (level-1 variance), or 2 (level-2 variance) only.",
          arg = "es_type")
    }
    if (covariate) {
        if (true_icc_covariate == 0 && true_R2_between != 0) .dmar_refuse_argument(
          "Because the covariate varies at the level 1 only, the r-square at level 2 must be 0.",
          arg = NULL)
        if (true_icc_covariate == 1 && true_R2_within != 0) .dmar_refuse_argument(
          "Because the covariate varies at the level 2 only, the r-square at level 1 must be 0.",
          arg = NULL)
        tauz <- total_var_Z * true_icc_covariate
        sigmaz <- total_var_Z * (1 - true_icc_covariate)
        gammazb <- sqrt(true_R2_between * tau / tauz)
        gammazw <- sqrt(true_R2_within * sigma / sigmaz)
        tau <- (1 - true_R2_between) * tau
        sigma <- (1 - true_R2_within) * sigma
    }
    probx <- sum(datawide[, x_label]) / nrow(datawide)

    if (!requireNamespace("OpenMx", quietly = TRUE)) .dmar_refuse(
      "The package 'OpenMx' is needed; please install the package and try again.",
      "package_required", package = "OpenMx")

    latent_label <- c("intcept", "slope")
    if (is.null(z_within_label)) latent_label <- "intcept"

    frowlab <- c(y_label, x_label, z_between_label)
    fcollab <- c(frowlab, latent_label)
    lenrow <- length(frowlab)
    lencol <- length(fcollab)

    Alab <- matrix(NA, lencol, lencol)
    Aval <- matrix(0, lencol, lencol)
    Afree <- matrix(FALSE, lencol, lencol)
    colnames(Alab) <- colnames(Aval) <- colnames(Afree) <- rownames(Alab) <- rownames(Aval) <- rownames(Afree) <- fcollab

    Alab["intcept", x_label] <- "groupdiff"
    if (!is.null(z_between_label)) Alab["intcept", z_between_label] <- "zbeffect"
    if (!is.null(z_within_label)) Alab[y_label, "slope"] <- paste("data.", z_within_label, sep = "")

    Aval["intcept", x_label] <- gamma1
    if (!is.null(z_between_label)) Aval["intcept", z_between_label] <- gammazb
    Aval[y_label, latent_label] <- 1

    Afree["intcept", x_label] <- TRUE
    if (!is.null(z_between_label)) Afree["intcept", z_between_label] <- TRUE

    Slab <- matrix(NA, lencol, lencol)
    Sval <- matrix(0, lencol, lencol)
    Sfree <- matrix(FALSE, lencol, lencol)
    colnames(Slab) <- colnames(Sval) <- colnames(Sfree) <- rownames(Slab) <- rownames(Sval) <- rownames(Sfree) <- fcollab
    diag(Slab)[1:length(y_label)] <- "l1error"
    diag(Sval)[1:length(y_label)] <- sigma
    diag(Sfree)[1:length(y_label)] <- TRUE
    Slab["intcept", "intcept"] <- "l2error"
    Sval["intcept", "intcept"] <- tau
    Sfree["intcept", "intcept"] <- TRUE
    Slab[x_label, x_label] <- "varx"
    Sval[x_label, x_label] <- probx * (1 - probx)
    Sfree[x_label, x_label] <- TRUE
    if (!is.null(z_between_label)) {
        Slab[c(x_label, z_between_label), c(x_label, z_between_label)] <- "covxzb"
        Slab[x_label, x_label] <- "varx"
        Slab[z_between_label, z_between_label] <- "varzb"
        Sval[c(x_label, z_between_label), c(x_label, z_between_label)] <- 0
        Sval[x_label, x_label] <- probx * (1 - probx)
        Sval[z_between_label, z_between_label] <- tauz
        Sfree[c(x_label, z_between_label), c(x_label, z_between_label)] <- TRUE
    }

    Fval <- cbind(diag(lenrow), matrix(0, lenrow, length(latent_label)))
    Flab <- matrix(NA, lenrow, lencol)
    Ffree <- matrix(FALSE, lenrow, lencol)
    colnames(Flab) <- colnames(Fval) <- colnames(Ffree) <- fcollab
    rownames(Flab) <- rownames(Fval) <- rownames(Ffree) <- frowlab

    Mlab <- c(rep(NA, length(y_label)), "meanX")
    Mval <- c(rep(0, length(y_label)), probx)
    Mfree <- c(rep(FALSE, length(y_label)), TRUE)

    if (!is.null(z_between_label)) {
        Mlab <- c(Mlab, "meanzb")
        Mval <- c(Mval, 0)
        Mfree <- c(Mfree, TRUE)
    }
    Mlab <- c(Mlab, "meanctrl")
    Mval <- c(Mval, 0)
    Mfree <- c(Mfree, TRUE)
    if (!is.null(z_within_label)) {
        Mlab <- c(Mlab, "zweffect")
        Mval <- c(Mval, gammazw)
        Mfree <- c(Mfree, TRUE)
    }

    Mlab <- matrix(Mlab, 1, lencol)
    Mval <- matrix(Mval, 1, lencol)
    Mfree <- matrix(Mfree, 1, lencol)
    colnames(Mlab) <- colnames(Mval) <- colnames(Mfree) <- fcollab

    varzw <- 0
    if (!is.null(z_within_label)) varzw <- var(as.vector(as.matrix(datawide[, z_within_label])))

    # These lines prevent "no visible binding" notes when OpenMx evaluates the mxAlgebra expressions.
    groupdiff <- NULL
    l1error <- NULL
    l2error <- NULL
    zbeffect <- NULL
    varzb <- NULL
    zweffect <- NULL

    constraint <- NULL
    if (es_type == 0) {
        if (is.null(z_within_label)) {
            if (is.null(z_between_label)) {
                constraint <- OpenMx::mxAlgebra(expression = groupdiff / sqrt(l1error + l2error), name = "es")
            } else {
                constraint <- OpenMx::mxAlgebra(expression = groupdiff / sqrt(l1error + l2error + (zbeffect^2 * varzb)), name = "es")
            }
        } else {
            if (is.null(z_between_label)) {
                constraint <- OpenMx::mxAlgebra(expression = groupdiff / sqrt(l1error + (zweffect^2 * varzw) + l2error), name = "es")
            } else {
                constraint <- OpenMx::mxAlgebra(expression = groupdiff / sqrt(l1error + (zweffect^2 * varzw) + l2error + (zbeffect^2 * varzb)), name = "es")
            }
        }
    } else if (es_type == 1) {
        if (is.null(z_within_label)) {
            constraint <- OpenMx::mxAlgebra(expression = groupdiff / sqrt(l1error), name = "es")
        } else {
            constraint <- OpenMx::mxAlgebra(expression = groupdiff / sqrt(l1error + (zweffect^2 * varzw)), name = "es")
        }
    } else if (es_type == 2) {
        if (is.null(z_between_label)) {
            constraint <- OpenMx::mxAlgebra(expression = groupdiff / sqrt(l2error), name = "es")
        } else {
            constraint <- OpenMx::mxAlgebra(expression = groupdiff / sqrt(l2error + (zbeffect^2 * varzb)), name = "es")
        }
    } else {
        .dmar_refuse_argument(
          "'es_type' can be 0 (total variance), 1 (level-1 variance), or 2 (level-2 variance) only.",
          arg = "es_type")
    }
    # The expectation and fit function are the two objects the deprecated
    # mxRAMObjective() returned (it built exactly these and warned), so the
    # model and its fit are the ones the earlier call gave.
    onecov <- OpenMx::mxModel("effect size CRD", type = "RAM",
                              OpenMx::mxData(datawide, type = "raw"),
                              OpenMx::mxMatrix(type = "Full", nrow = lencol, ncol = lencol, values = Aval, free = Afree, labels = Alab, name = "A"),
                              OpenMx::mxMatrix(type = "Symm", nrow = lencol, ncol = lencol, values = Sval, free = Sfree, labels = Slab, name = "S"),
                              OpenMx::mxMatrix(type = "Full", nrow = lenrow, ncol = lencol, values = Fval, free = Ffree, labels = Flab, name = "F"),
                              OpenMx::mxMatrix(type = "Full", nrow = 1, ncol = lencol, values = Mval, free = Mfree, labels = Mlab, name = "M"),
                              OpenMx::mxMatrix(type = "Full", nrow = 1, ncol = 1, values = varzw, free = FALSE, labels = "varzw", name = "J"),
                              OpenMx::mxExpectationRAM("A", "S", "F", "M", dimnames = fcollab), OpenMx::mxFitFunctionML(),
                              constraint, OpenMx::mxCI(c("es"), interval = conf_level)
    )

    return(.crd_fit_interval(.crd_run_openmx(onecov)))
}

# Calculate a likelihood-based CI of ES (unequal cluster sizes)
.lik_ci_es_crd_unequal <- function(datawide, y_label, x_label, z_within_label = NULL, z_between_label = NULL, es_type = 1, true_icc = 0.25, es = 0.5, true_total_var = 1, covariate = FALSE, true_icc_covariate = 0.25, true_R2_within = 0.5, true_R2_between = 0.5, total_var_Z = 1, conf_level = 0.95) {
    if (!requireNamespace("OpenMx", quietly = TRUE)) .dmar_refuse(
      "The package 'OpenMx' is needed; please install the package and try again.",
      "package_required", package = "OpenMx")

    tau <- true_icc * true_total_var
    sigma <- (1 - true_icc) * true_total_var
    gamma1 <- NULL
    if (es_type == 0) {
        gamma1 <- es * sqrt(true_total_var)
    } else if (es_type == 1) {
        gamma1 <- es * sqrt(sigma)
    } else if (es_type == 2) {
        gamma1 <- es * sqrt(tau)
    } else {
        .dmar_refuse_argument(
          "'es_type' can be 0 (total variance), 1 (level-1 variance), or 2 (level-2 variance) only.",
          arg = "es_type")
    }

    if (covariate) {
        if (true_icc_covariate == 0 && true_R2_between != 0) .dmar_refuse_argument(
          "Because the covariate varies at the level 1 only, the r-square at level 2 must be 0.",
          arg = NULL)
        if (true_icc_covariate == 1 && true_R2_within != 0) .dmar_refuse_argument(
          "Because the covariate varies at the level 2 only, the r-square at level 1 must be 0.",
          arg = NULL)
        tauz <- total_var_Z * true_icc_covariate
        sigmaz <- total_var_Z * (1 - true_icc_covariate)
        gammazb <- sqrt(true_R2_between * tau / tauz)
        gammazw <- sqrt(true_R2_within * sigma / sigmaz)
        tau <- (1 - true_R2_between) * tau
        sigma <- (1 - true_R2_within) * sigma
    }

    ntreat <- sum(sapply(datawide, function(x) sum(x[, x_label])))
    totaln <- sum(sapply(datawide, nrow))
    probx <- ntreat / totaln

    FUNgroupsize <- function(dat, y, zw = NULL) {
        latent_label <- c("intcept", "slope")
        if (is.null(zw)) latent_label <- "intcept"

        frowlab <- c(y, x_label, z_between_label)
        fcollab <- c(frowlab, latent_label)
        lenrow <- length(frowlab)
        lencol <- length(fcollab)

        Alab <- matrix(NA, lencol, lencol)
        Aval <- matrix(0, lencol, lencol)
        Afree <- matrix(FALSE, lencol, lencol)
        colnames(Alab) <- colnames(Aval) <- colnames(Afree) <- rownames(Alab) <- rownames(Aval) <- rownames(Afree) <- fcollab

        Alab["intcept", x_label] <- "groupdiff"
        if (!is.null(z_between_label)) Alab["intcept", z_between_label] <- "zbeffect"
        if (!is.null(zw)) Alab[y, "slope"] <- paste("data.", zw, sep = "")

        Aval["intcept", x_label] <- gamma1
        if (!is.null(z_between_label)) Aval["intcept", z_between_label] <- gammazb
        Aval[y, latent_label] <- 1

        Afree["intcept", x_label] <- TRUE
        if (!is.null(z_between_label)) Afree["intcept", z_between_label] <- TRUE

        Slab <- matrix(NA, lencol, lencol)
        Sval <- matrix(0, lencol, lencol)
        Sfree <- matrix(FALSE, lencol, lencol)
        colnames(Slab) <- colnames(Sval) <- colnames(Sfree) <- rownames(Slab) <- rownames(Sval) <- rownames(Sfree) <- fcollab
        diag(Slab)[1:length(y)] <- "l1error"
        diag(Sval)[1:length(y)] <- sigma
        diag(Sfree)[1:length(y)] <- TRUE
        Slab["intcept", "intcept"] <- "l2error"
        Sval["intcept", "intcept"] <- tau
        Sfree["intcept", "intcept"] <- TRUE
        Slab[x_label, x_label] <- "varx"
        Sval[x_label, x_label] <- probx * (1 - probx)
        Sfree[x_label, x_label] <- TRUE
        if (!is.null(z_between_label)) {
            Slab[c(x_label, z_between_label), c(x_label, z_between_label)] <- "covxzb"
            Slab[x_label, x_label] <- "varx"
            Slab[z_between_label, z_between_label] <- "varzb"
            Sval[c(x_label, z_between_label), c(x_label, z_between_label)] <- 0
            Sval[x_label, x_label] <- probx * (1 - probx)
            Sval[z_between_label, z_between_label] <- tauz
            Sfree[c(x_label, z_between_label), c(x_label, z_between_label)] <- TRUE
        }

        Fval <- cbind(diag(lenrow), matrix(0, lenrow, length(latent_label)))
        Flab <- matrix(NA, lenrow, lencol)
        Ffree <- matrix(FALSE, lenrow, lencol)
        colnames(Flab) <- colnames(Fval) <- colnames(Ffree) <- fcollab
        rownames(Flab) <- rownames(Fval) <- rownames(Ffree) <- frowlab

        Mlab <- c(rep(NA, length(y)), "meanX")
        Mval <- c(rep(0, length(y)), probx)
        Mfree <- c(rep(FALSE, length(y)), TRUE)

        if (!is.null(z_between_label)) {
            Mlab <- c(Mlab, "meanzb")
            Mval <- c(Mval, 0)
            Mfree <- c(Mfree, TRUE)
        }
        Mlab <- c(Mlab, "meanctrl")
        Mval <- c(Mval, 0)
        Mfree <- c(Mfree, TRUE)
        if (!is.null(zw)) {
            Mlab <- c(Mlab, "zweffect")
            Mval <- c(Mval, gammazw)
            Mfree <- c(Mfree, TRUE)
        }

        Mlab <- matrix(Mlab, 1, lencol)
        Mval <- matrix(Mval, 1, lencol)
        Mfree <- matrix(Mfree, 1, lencol)
        colnames(Mlab) <- colnames(Mval) <- colnames(Mfree) <- fcollab

        onecov <- OpenMx::mxModel(paste0("group", length(y)), type = "RAM",
                                  OpenMx::mxData(dat, type = "raw"),
                                  OpenMx::mxMatrix(type = "Full", nrow = lencol, ncol = lencol, values = Aval, free = Afree, labels = Alab, name = "A"),
                                  OpenMx::mxMatrix(type = "Symm", nrow = lencol, ncol = lencol, values = Sval, free = Sfree, labels = Slab, name = "S"),
                                  OpenMx::mxMatrix(type = "Full", nrow = lenrow, ncol = lencol, values = Fval, free = Ffree, labels = Flab, name = "F"),
                                  OpenMx::mxMatrix(type = "Full", nrow = 1, ncol = lencol, values = Mval, free = Mfree, labels = Mlab, name = "M"),
                                  OpenMx::mxMatrix(type = "Full", nrow = 1, ncol = 1, values = varzw, free = FALSE, labels = "varzw", name = "J"),
                                  OpenMx::mxExpectationRAM("A", "S", "F", "M", dimnames = fcollab), OpenMx::mxFitFunctionML()
        )
        return(onecov)
    }
    varzw <- 0
    if (!is.null(z_within_label)) varzw <- weighted.mean(do.call(c, mapply(function(x, y) var(as.vector(as.matrix(x[, y]))), x = datawide, y = z_within_label, SIMPLIFY = FALSE)), as.numeric(names(datawide)))

    # These lines prevent "no visible binding" notes when OpenMx evaluates the mxAlgebra expressions.
    groupdiff <- NULL
    l1error <- NULL
    l2error <- NULL
    zbeffect <- NULL
    varzb <- NULL
    zweffect <- NULL

    constraint <- NULL
    if (es_type == 0) {
        if (is.null(z_within_label)) {
            if (is.null(z_between_label)) {
                constraint <- OpenMx::mxAlgebra(expression = groupdiff / sqrt(l1error + l2error), name = "es")
            } else {
                constraint <- OpenMx::mxAlgebra(expression = groupdiff / sqrt(l1error + l2error + (zbeffect^2 * varzb)), name = "es")
            }
        } else {
            if (is.null(z_between_label)) {
                constraint <- OpenMx::mxAlgebra(expression = groupdiff / sqrt(l1error + (zweffect^2 * varzw) + l2error), name = "es")
            } else {
                constraint <- OpenMx::mxAlgebra(expression = groupdiff / sqrt(l1error + (zweffect^2 * varzw) + l2error + (zbeffect^2 * varzb)), name = "es")
            }
        }
    } else if (es_type == 1) {
        if (is.null(z_within_label)) {
            constraint <- OpenMx::mxAlgebra(expression = groupdiff / sqrt(l1error), name = "es")
        } else {
            constraint <- OpenMx::mxAlgebra(expression = groupdiff / sqrt(l1error + (zweffect^2 * varzw)), name = "es")
        }
    } else if (es_type == 2) {
        if (is.null(z_between_label)) {
            constraint <- OpenMx::mxAlgebra(expression = groupdiff / sqrt(l2error), name = "es")
        } else {
            constraint <- OpenMx::mxAlgebra(expression = groupdiff / sqrt(l2error + (zbeffect^2 * varzb)), name = "es")
        }
    } else {
        .dmar_refuse_argument(
          "'es_type' can be 0 (total variance), 1 (level-1 variance), or 2 (level-2 variance) only.",
          arg = "es_type")
    }
    list_model <- NULL
    if (!is.null(z_within_label)) {
        list_model <- mapply(FUNgroupsize, dat = datawide, y = y_label, zw = z_within_label)
    } else {
        list_model <- mapply(FUNgroupsize, dat = datawide, y = y_label)
    }
    title <- "Effect Size CRD"
    # Build the summed objective from the per-group objectives as one
    # parsed expression. (The old idiom of constructing an empty
    # mxAlgebra("") and patching its @formula slot fails on current
    # OpenMx, which validates the expression eagerly.) Each group's fit
    # function is named "fitfunction"; mxFitFunctionAlgebra() is what the
    # deprecated mxAlgebraObjective() built, without its warning.
    groupnames <- paste0("group", names(datawide), ".fitfunction")
    algebra <- OpenMx::mxAlgebraFromString(
        paste0("sum(", paste(groupnames, collapse = ", "), ")"),
        name = "allobjective")
    objective <- OpenMx::mxFitFunctionAlgebra("allobjective")
    finalmodel <- OpenMx::mxModel(title, OpenMx::mxMatrix(type = "Full", nrow = 1, ncol = 1, values = varzw, free = FALSE, labels = "varzw", name = "J"), unlist(list_model), constraint, algebra, objective, OpenMx::mxCI(c("es"), interval = conf_level))
    return(.crd_fit_interval(.crd_run_openmx(finalmodel)))
}

# Create dataset from the CRD model and transform it to the wide format
.create_data_crd_wide <- function(n_clusters, n_treat_clus, cluster_size, true_icc, es, es_type = 1, true_total_var = 1, covariate = FALSE, true_icc_covariate = NULL, true_R2_within = NULL, true_R2_between = NULL, total_var_Z = 1, diff_size = NULL) {
    dat <- .create_data_crd(n_clusters = n_clusters, n_treat_clus = n_treat_clus, cluster_size = cluster_size, true_icc = true_icc, es = es, es_type = es_type, true_total_var = true_total_var, covariate = covariate, true_icc_covariate = true_icc_covariate, true_R2_within = true_R2_within, true_R2_between = true_R2_between, total_var_Z = total_var_Z, diff_size = diff_size)
    datawide <- NULL
    if (!is.null(diff_size)) {
        if (covariate) {
            if (true_icc_covariate == 0) {
                datawide <- .wide_format_unequal(dat, 3, c(2, 4), 1)
            } else if (true_icc_covariate == 1) {
                datawide <- .wide_format_unequal(dat, c(3, 5), 2, 1)
            } else {
                datawide <- .wide_format_unequal(dat, c(3, 5), c(2, 4), 1)
            }
        } else {
            datawide <- .wide_format_unequal(dat, 3, 2, 1)
        }
    } else {
        if (covariate) {
            if (true_icc_covariate == 0) {
                datawide <- .wide_format(dat, 3, c(2, 4), 1)
            } else if (true_icc_covariate == 1) {
                datawide <- .wide_format(dat, c(3, 5), 2, 1)
            } else {
                datawide <- .wide_format(dat, c(3, 5), c(2, 4), 1)
            }
        } else {
            datawide <- .wide_format(dat, 3, 2, 1)
        }
    }
    return(datawide)
}

# Reshape long-format data to wide format
.wide_format <- function(data, between_col, within_col, id_col) {
    temp <- split(data[, within_col], data[, id_col])
    temp <- lapply(temp, function(x) as.vector(as.matrix(x)))
    dataw <- do.call(rbind, temp)
    datab <- as.matrix(data[match(unique(data[, id_col]), data[, id_col]), between_col])
    colnames(datab) <- colnames(data)[between_col]
    cluster_size <- nrow(data) / nrow(datab)
    colnames(dataw) <- paste(rep(colnames(data)[within_col], each = cluster_size), rep(1:cluster_size, length(within_col)), sep = "")
    return(data.frame(datab, dataw))
}

# Reshape long-format data to wide format with unequal cluster sizes
.wide_format_unequal <- function(data, between_col, within_col, id_col) {
    temp <- split(data[, within_col], data[, id_col])
    temp <- lapply(temp, function(x) as.vector(as.matrix(x)))

    size <- sapply(temp, length) / length(within_col)
    dataw <- lapply(split(temp, size), function(x) do.call(rbind, x))
    datab <- split(data[match(unique(data[, id_col]), data[, id_col]), between_col], size)
    resultdat <- mapply(data.frame, datab, dataw)
    varnamesw <- lapply(sapply(dataw, ncol) / length(within_col), function(x) paste(rep(colnames(data)[within_col], each = x), rep(1:x, length(within_col)), sep = ""))
    varnames <- lapply(varnamesw, function(x) c(colnames(data)[between_col], x))
    resultdat <- mapply(function(x, y) { colnames(x) <- y; x }, x = resultdat, y = varnames)
    return(resultdat)
}

# Create dataset from the CRD model
.create_data_crd <- function(n_clusters, n_treat_clus, cluster_size, true_icc, es, es_type = 1, true_total_var = 1, covariate = FALSE, true_icc_covariate = NULL, true_R2_within = NULL, true_R2_between = NULL, total_var_Z = 1, diff_size = NULL) {
    if (!requireNamespace("MASS", quietly = TRUE)) .dmar_refuse(
      "The package 'MASS' is needed; please install the package and try again.",
      "package_required", package = "MASS")

    if (!is.null(diff_size)) {
        cluster_size <- .find_n_indiv_vec(cluster_size, diff_size, n_clusters)
    } else {
        cluster_size <- rep(cluster_size, each = n_clusters)
    }
    id <- do.call(c, mapply(rep, 1:n_clusters, cluster_size, SIMPLIFY = FALSE))
    x <- c(rep(1, n_treat_clus), rep(0, n_clusters - n_treat_clus))
    tau <- true_icc * true_total_var
    sigma <- (1 - true_icc) * true_total_var
    gamma1 <- NULL
    if (es_type == 0) {
        gamma1 <- es * sqrt(true_total_var)
    } else if (es_type == 1) {
        gamma1 <- es * sqrt(sigma)
    } else if (es_type == 2) {
        gamma1 <- es * sqrt(tau)
    } else {
        .dmar_refuse_argument(
          "'es_type' can be 0 (total variance), 1 (level-1 variance), or 2 (level-2 variance) only.",
          arg = "es_type")
    }
    if (covariate) {
        if (true_icc_covariate == 0 && true_R2_between != 0) .dmar_refuse_argument(
          "Because the covariate varies at the level 1 only, the r-square at level 2 must be 0.",
          arg = NULL)
        if (true_icc_covariate == 1 && true_R2_within != 0) .dmar_refuse_argument(
          "Because the covariate varies at the level 2 only, the r-square at level 1 must be 0.",
          arg = NULL)
        tauz <- total_var_Z * true_icc_covariate
        sigmaz <- total_var_Z * (1 - true_icc_covariate)
        gammazb <- sqrt(true_R2_between * tau / tauz)
        gammazw <- sqrt(true_R2_within * sigma / sigmaz)
        tau <- (1 - true_R2_between) * tau
        sigma <- (1 - true_R2_within) * sigma
    }
    ybetween <- (gamma1 * x) + rnorm(n_clusters, 0, sqrt(tau))
    if (covariate && true_icc_covariate != 0) {
        zbetween <- rnorm(n_clusters, 0, sqrt(tauz))
        ybetween <- ybetween + gammazb * zbetween
    }
    ybetween <- do.call(c, mapply(rep, ybetween, cluster_size, SIMPLIFY = FALSE))
    y <- ybetween + rnorm(sum(cluster_size), 0, sqrt(sigma))
    if (covariate && true_icc_covariate != 1) {
        zwithin <- rnorm(sum(cluster_size), 0, sqrt(sigmaz))
        y <- y + gammazw * zwithin
    }
    x <- do.call(c, mapply(rep, x, cluster_size, SIMPLIFY = FALSE))
    dat <- data.frame(id, y, x)
    if (covariate) {
        z <- NULL
        if (true_icc_covariate == 0) {
            z <- zwithin
        } else if (true_icc_covariate == 1) {
            z <- do.call(c, mapply(rep, zbetween, cluster_size, SIMPLIFY = FALSE))
        } else {
            z <- do.call(c, mapply(rep, zbetween, cluster_size, SIMPLIFY = FALSE)) + zwithin
        }
        if (true_icc_covariate == 0) {
            zb <- 0
            zw <- z
        } else {
            zlist <- split(z, id)
            zb <- do.call(c, mapply(rep, sapply(zlist, mean), cluster_size, SIMPLIFY = FALSE))
            zw <- do.call(c, lapply(zlist, scale, scale = FALSE))
        }
        z <- data.frame(zw = zw, zb = zb)
        dat <- data.frame(dat, z)
    }
    rownames(dat) <- NULL
    return(dat)
}

# Find the number of clusters given the specified width and the cluster
# size: the smallest number of clusters, from the floor
# .crd_min_n_clusters() (at least two clusters in each condition), whose
# width is no larger than the target. The widths at 100, 200, ..., 1000
# clusters bracket the target (in blocks of 1000 above that), and the
# bracket is searched one cluster at a time. When the width falls as
# clusters are added, the search returns the smallest number that meets
# the target; with unequal cluster sizes (diff_size) or an assurance below
# .5 the width need not fall at every step, and it returns a number that
# meets the target, which need not be the smallest. A width meets the
# target when it is no larger than it (through the first builds of 1.1.0 a
# width had to be below the target, so a target equal to the width of a
# design was met by the next number of clusters). A number of clusters
# without a degree of freedom (below 3 + n_covariates, which only the
# grid can try) does not meet the target; any other width that is not a
# number is refused (.crd_refuse_width_nan()), since the search would
# otherwise go on without end. When no number up to 1000 meets the target
# the width at .CRD_MAX_N_CLUSTERS is computed first, and a target it does
# not meet is refused (.crd_refuse_n_clusters_beyond()); otherwise the
# search ends by that number. Through the first builds of 1.1.0 a width
# equal to the target at the end of a bracket stopped the search with
# "'from' must be a finite number"; when every width from 100 to 1000 met
# the target and so did every one from the floor to 100, the search
# returned 3 + n_covariates clusters, a number below the floor whose width
# it had not computed (with pr_treat = .1 a design with no treated cluster
# and no finite width); with the floor above 100 (pr_treat below .02) it
# stopped with "wrong sign in 'by' argument" or returned a number below
# the floor; and a target too narrow for any study (width = 1e-300) was
# searched for without end. `fn` names the public planner in a refusal.
.find_n_clus_crd_diff <- function(width, cluster_size, pr_treat, true_var_between = NULL, true_var_within = NULL, true_total_var = NULL, true_icc = NULL, true_R2_between = 0, true_R2_within = 0, n_covariates = 0, assurance = NULL, conf_level = 0.95, diff_size = NULL, fn = NULL) {
    lo <- .crd_min_n_clusters(pr_treat, n_covariates)
    # TRUE where the width at a number of clusters meets the target.
    meets <- function(n_clusters) {
        result <- vapply(n_clusters, .find_width_crd_diff, numeric(1), cluster_size = cluster_size, pr_treat = pr_treat, true_var_between = true_var_between, true_var_within = true_var_within, true_total_var = true_total_var, true_icc = true_icc, true_R2_between = true_R2_between, true_R2_within = true_R2_within, n_covariates = n_covariates, assurance = assurance, conf_level = conf_level, diff_size = diff_size)
        has_df <- n_clusters - 2 - n_covariates >= 1
        nan <- which(is.na(result) & has_df)
        if (length(nan)) .crd_refuse_width_nan(fn, n_clusters[nan[1]], cluster_size, diff_size)
        has_df & !is.na(result) & width >= result
    }
    n_clusters <- seq(100, 1000, 100)
    ok <- meets(n_clusters)
    if (all(ok)) {
        # Every number from 100 to 1000 meets the target: the first from
        # the floor that does.
        n_clusters <- seq.int(lo, max(lo, 100))
        ok <- meets(n_clusters)
        if (any(ok)) return(n_clusters[which(ok)[1]])
        return(.crd_first_meeting(lo + 1, meets))
    }
    if (!any(ok)) {
        if (!meets(.CRD_MAX_N_CLUSTERS)) {
            at_max <- .find_width_crd_diff(n_clusters = .CRD_MAX_N_CLUSTERS, cluster_size = cluster_size, pr_treat = pr_treat, true_var_between = true_var_between, true_var_within = true_var_within, true_total_var = true_total_var, true_icc = true_icc, true_R2_between = true_R2_between, true_R2_within = true_R2_within, n_covariates = n_covariates, assurance = assurance, conf_level = conf_level, diff_size = diff_size)
            .crd_refuse_n_clusters_beyond(fn, width, at_max, cluster_size = cluster_size, assurance = assurance)
        }
        start <- 1000
        repeat {
            n_clusters <- seq(start, start + 1000, 100)
            ok <- meets(n_clusters)
            if (any(ok)) break
            start <- start + 1000
        }
    }
    first <- which(ok)[1]
    # The bracket: the last number of the grid that does not meet the
    # target and the first that does. The first of the grid meets it only
    # when the width does not fall at every step (100 meets and a larger
    # number does not); the search then starts from the floor.
    n_clusters <- if (first == 1L) seq.int(lo, max(lo, n_clusters[1])) else seq(n_clusters[first - 1], n_clusters[first], 1)
    ok <- meets(n_clusters)
    n <- if (any(ok)) n_clusters[which(ok)[1]] else NA
    # A number below the floor (possible when the floor is above 100) gives
    # way to the first number from the floor that meets the target.
    if (is.na(n) || n < lo) n <- .crd_first_meeting(lo, meets)
    n
}

# Find width of the CI of the unstandardized conditions difference
.find_width_crd_diff <- function(n_clusters, cluster_size, pr_treat, true_var_between = NULL, true_var_within = NULL, true_total_var = NULL, true_icc = NULL, true_R2_between = 0, true_R2_within = 0, n_covariates = 0, assurance = NULL, conf_level = 0.95, diff_size = NULL) {
    n_clusters <- as.numeric(n_clusters)
    pr_treat <- round(n_clusters * pr_treat) / n_clusters
    cluster_size <- as.numeric(cluster_size)
    if (!is.null(diff_size)) cluster_size <- .find_n_indiv_vec(cluster_size, diff_size, n_clusters)
    v <- .var_crd_diff(n_clusters = n_clusters, cluster_size = cluster_size, pr_treat = pr_treat, true_var_between = true_var_between, true_var_within = true_var_within, true_total_var = true_total_var, true_icc = true_icc, true_R2_between = true_R2_between, true_R2_within = true_R2_within, n_covariates = n_covariates, assurance = assurance)
    alpha <- 1 - conf_level
    width <- 2 * sqrt(v) * qt(1 - alpha / 2, n_clusters - 2 - n_covariates)
    return(width)
}

# Find the variance of the conditions difference with or without a
# covariate. The public planners resolve the four variance inputs once
# (.crd_resolve_variances()) and pass the variance between clusters and
# the variance within them; the standardized planners pass the total
# variance and the intraclass correlation, from which the two are formed
# here. A design that leaves a condition without a cluster has an
# infinite variance, so its width does not meet any target (with equal
# cluster sizes the denominator is then 0; the equal-size branch takes a
# vector of numbers of clusters). Through the first builds of 1.1.0 the
# unequal-size branch summed the treated clusters' terms over 1:ntreat,
# which with no treated cluster is c(1, 0) and took the first cluster's
# term, and gave NaN there.
.var_crd_diff <- function(n_clusters, cluster_size, pr_treat, true_var_between = NULL, true_var_within = NULL, true_total_var = NULL, true_icc = NULL, true_R2_between = 0, true_R2_within = 0, n_covariates = 0, assurance = NULL) {
    if (is.null(true_var_between)) {
        if (is.null(true_icc)) {
            true_var_between <- true_total_var - true_var_within
        } else {
            true_var_between <- true_icc * true_total_var
        }
    }
    if (is.null(true_var_within)) {
        if (is.null(true_icc)) {
            true_var_within <- true_total_var - true_var_between
        } else {
            true_var_within <- (1 - true_icc) * true_total_var
        }
    }
    ntreat <- round(pr_treat * n_clusters)
    if (length(cluster_size) > 1) {
        if (ntreat == 0 || ntreat == n_clusters) return(Inf)
        true_var_within <- true_var_within * (1 - true_R2_within)
        true_var_between <- true_var_between * (1 - true_R2_between)
        D <- 1 / true_var_within
        F <- 0 - (true_var_between) / (true_var_within * (true_var_within + cluster_size * true_var_between))
        pr_treat <- ntreat / n_clusters
        terms <- cluster_size * (D + cluster_size * F)
        termstotal <- sum(terms)
        termstreatment <- sum(terms[seq_len(ntreat)])
        vargamma1 <- termstotal / (termstotal * termstreatment - termstreatment^2)
        varclusmean <- vargamma1 * (n_clusters * pr_treat * (1 - pr_treat))
    } else {
        pr_treat <- ntreat / n_clusters
        varclusmean <- ((true_var_within * (1 - true_R2_within)) / cluster_size) + (true_var_between * (1 - true_R2_between))
    }
    if (!is.null(assurance)) {
        df <- n_clusters - 2 - n_covariates
        varclusmean <- varclusmean * qchisq(assurance, df) / df
    }
    denominator <- n_clusters * pr_treat * (1 - pr_treat)
    return(varclusmean / denominator)
}

# Find the total cost of the whole cluster randomized design
.cost_crd <- function(n_clusters, cluster_size, cluster_cost, individual_cost, diff_size = NULL) {
    cluster_size <- as.numeric(cluster_size)
    if (!is.null(diff_size)) cluster_size <- .find_n_indiv_vec(cluster_size, diff_size, n_clusters)
    result <- NULL
    if (length(cluster_size) > 1) {
        result <- (n_clusters * cluster_cost) + sum(cluster_size * individual_cost)
    } else {
        result <- n_clusters * (cluster_cost + (cluster_size * individual_cost))
    }
    return(result)
}

# Find the cluster size given the specified width of ES and the number of
# clusters. With a fixed number of clusters the width of the interval for
# the standardized effect levels off as the clusters grow, so a target can
# be out of reach at any cluster size. For the effect standardized by the
# within-cluster standard deviation (es_type = 1) with no covariate at the
# cluster level, the value it levels off at is known exactly
# (.crd_es_width_limit()), and a target at or below it is refused before
# anything is simulated (class dmar_argument_invalid, reason
# "width_not_reachable"). There the search starts from the smaller of the
# cluster size the unstandardized conditions difference needs
# (.find_n_indiv_crd_diff()) and the closed-form start
# (.crd_es_cluster_size_start()); elsewhere it starts from the former,
# which is refused when it is beyond 100,000. A start above
# .CRD_ES_MAX_CLUSTER_SIZE is refused before anything is simulated (reason
# "cluster_size_beyond_reach"), and the refusal calls the start a lower
# bound on the size the target needs only where it is one (see below).
# From the start the search steps down one individual at a time, or
# brackets and bisects upward (.crd_es_search_up()) through cluster sizes
# at most .CRD_ES_CLUSTER_SIZE_REACH above the size the unstandardized
# difference needs (above the start when that size is beyond 100,000) and
# at most .CRD_ES_MAX_CLUSTER_SIZE, since a fit takes far longer as the
# clusters grow (the OpenMx model has one variable per individual of a
# cluster). Through the first builds of 1.1.0 the upward search stepped
# one individual at a time and had no end when the width leveled off
# above the target, and the start was the closed-form one only when the
# unstandardized difference needed more than 100,000. An Inf width (see
# .find_width_crd_es()) does not meet the target.
.find_n_indiv_crd_es <- function(width, n_clusters, es, es_type = 1, true_icc, pr_treat, true_R2_between = 0, true_R2_within = 0, n_covariates = 0, assurance = NULL, conf_level = 0.95, G = 1000, true_icc_covariate = NULL, seed = NULL, ncpus = 1, diff_size = NULL, tracker = NULL) {
    if (n_covariates > 0 && is.null(true_icc_covariate)) true_icc_covariate <- true_icc
    if (n_covariates > 1) .dmar_refuse_argument(
      "'n_covariates' must be 0 or 1: the planners for the standardized effect size allow at most one covariate beside the treatment indicator.",
      arg = "n_covariates")
    true_total_var <- 1
    if (es_type == 0) {
        true_total_var <- 1
    } else if (es_type == 1) {
        true_total_var <- 1 / (1 - true_icc)
    } else if (es_type == 2) {
        true_total_var <- 1 / true_icc
    } else {
        .dmar_refuse_argument(
          "'es_type' can be 0 (total variance), 1 (level-1 variance), or 2 (level-2 variance) only.",
          arg = "es_type")
    }
    # The width the interval approaches as the clusters grow, where it is
    # known exactly: es_type 1 with no covariate at the cluster level (a
    # covariate with an intraclass correlation of 0 varies within clusters
    # only).
    n_treat_clus <- round(n_clusters * pr_treat)
    width_limit <- if (es_type == 1 && (n_covariates == 0 || isTRUE(true_icc_covariate == 0))) {
        .crd_es_width_limit(n_clusters, n_treat_clus, true_icc, conf_level = conf_level, assurance = assurance)
    } else {
        NA_real_
    }
    kind <- if (is.null(assurance)) "expected" else "assured"
    if (is.finite(width_limit) && width <= width_limit) {
        .dmar_refuse_argument(c(
            .crd_es_fun(tracker), ": with ", n_clusters, " clusters (", n_treat_clus,
            " of them treated), the ", kind, " width of the interval for the effect ",
            "standardized by the within-cluster standard deviation falls toward ",
            format(signif(width_limit, 4)), " as the clusters grow, its value with ",
            "every cluster mean known exactly, and does not reach it. The target ",
            "'width' of ", format(width), " is not above that, so no cluster size ",
            "reaches it with this number of clusters. ",
            "ss_aipe_crd_es_n_clusters_fixed_width() plans the number of clusters ",
            "instead."),
            arg = "width", reason = "width_not_reachable", width_limit = width_limit,
            target = width, n_clusters = n_clusters)
    }
    # The cluster size the unstandardized conditions difference needs, from
    # its width with the t quantile on K - 2 - n_covariates degrees of
    # freedom; Inf when it is beyond 100,000.
    diff_start <- .find_n_indiv_crd_diff(width = width, n_clusters = n_clusters, pr_treat = pr_treat, true_total_var = true_total_var, true_icc = true_icc, true_R2_between = true_R2_between, true_R2_within = true_R2_within, n_covariates = n_covariates, assurance = assurance, conf_level = conf_level, diff_size = diff_size, fn = tracker$fun)
    diff_start <- as.numeric(diff_start)
    if (is.finite(width_limit)) {
        # With the limit known and the target above it (checked above), the
        # search starts from the smaller of diff_start and the closed-form
        # start (.crd_es_cluster_size_start()). With few clusters the t
        # quantile makes diff_start far too large: with 10 clusters and an
        # intraclass correlation of 0.2 it is 1,875 for a target of 1.46,
        # where about 10 individuals per cluster reach it with a
        # standardized effect of 0.5 (diff_start does not depend on the
        # effect, but the size that reaches the target grows with it: 10 or
        # 11 with an effect of 1, 13 or 14 with 2; closed-form start 8), and
        # it is beyond 100,000 for every target from the limit, 1.187, to
        # about 1.458. With equal cluster sizes and no covariate (the
        # proportions of variance explained at 0) the ratio of the
        # unstandardized width to the closed-form width is the same at
        # every cluster size and above 1 (qt(1 - alpha / 2, K - 2)
        # over sqrt(exp(q / K) - 1) E sqrt(chi^2_{K - 2}), or over
        # sqrt((K - 2) (exp(q / K) - 1)) with an assurance), so the start is
        # then the closed-form one. true_R2_within enters the closed-form
        # start only with a covariate: without one the simulated data do
        # not use it (.create_data_crd()), and a negative value would raise
        # the start above the size the target needs.
        r2w_start <- if (n_covariates > 0) true_R2_within else 0
        startval <- min(.crd_es_cluster_size_start(width, width_limit, true_icc, r2w_start), diff_start)
    } else {
        if (!is.finite(diff_start)) .dmar_refuse_argument(c(
            .crd_es_fun(tracker), ": the search for the cluster size starts from ",
            "the cluster size the unstandardized difference between the condition ",
            "means needs, and with ", n_clusters, " clusters (", n_treat_clus,
            " of them treated) no cluster size up to 100,000 brings the width of ",
            "that interval to the target 'width' of ", format(width), ", so the ",
            "search has no start. With a fixed number of clusters the width levels ",
            "off as the clusters grow, because the clusters themselves vary. ",
            "ss_aipe_crd_es_n_clusters_fixed_width() plans the number of clusters ",
            "instead."),
            arg = "width", reason = "cluster_size_beyond_reach", start = Inf,
            lower_bound = FALSE, max_cluster_size = .CRD_ES_MAX_CLUSTER_SIZE,
            width_limit = width_limit, target = width, n_clusters = n_clusters)
        startval <- diff_start
    }
    # The start is a lower bound on the cluster size the target needs for
    # the expected width with no covariate and equal cluster sizes: at
    # every smaller cluster size the expected width of the interval is
    # above the target. With the within-cluster standard deviation sigma
    # known, the expected width at a cluster size m is exactly
    # width_limit * sqrt(1 + v / (tau m)), v / tau = (1 - icc) / icc (see
    # .crd_es_width_limit()). With sigma estimated, the log likelihood
    # separates into the within-cluster sum of squares, which carries
    # sigma alone, and the K cluster means, whose variance
    # theta = tau + sigma^2 / m is free (tau has no lower bound in the
    # model), so the deviance of the effect at es is no larger than its
    # deviance with sigma fixed at its maximum likelihood estimate
    # sigma_hat. The interval therefore contains the one with sigma fixed
    # at sigma_hat, whose width is the known-sigma width times
    # sigma / sigma_hat; sigma_hat is independent of the cluster means, and
    # E(sigma / sigma_hat) >= sigma / E(sigma_hat) >= sigma /
    # sqrt(E(sigma_hat^2)) = 1 by Jensen's inequality (1 / x is convex and
    # sqrt(x) concave; sigma_hat^2 = SSW / (K (m - 1)) has mean sigma^2), with
    # equality only in the limit. The closed-form start is the smallest m at
    # which width_limit * sqrt(1 + v / (tau m)) falls below the target, so
    # at every smaller m the expected width is above the target (it does
    # not meet it), and the start, which is no larger, is a lower bound
    # too. The argument says nothing of a quantile of the width (an
    # assurance), of a covariate, whose slope the cluster means also
    # carry, or of unequal cluster sizes, where the variances of the
    # cluster means differ, so the start is not called a lower bound there.
    lower_bound <- is.finite(width_limit) && n_covariates == 0 && is.null(diff_size) &&
        is.null(assurance)
    if (startval > .CRD_ES_MAX_CLUSTER_SIZE) {
        .dmar_refuse_argument(c(
            .crd_es_fun(tracker), ": the search for the cluster size would start at ",
            format(startval, big.mark = ",", scientific = FALSE, trim = TRUE),
            " individuals per cluster",
            if (lower_bound) {
                c(", a lower bound on the size the target needs with ", n_clusters,
                  " clusters: the expected width is above the target at every smaller ",
                  "cluster size, and falls toward ", format(signif(width_limit, 4)),
                  " as the clusters grow")
            } else if (is.finite(width_limit)) {
                c(" (the ", kind, " width falls toward ", format(signif(width_limit, 4)),
                  " as the clusters grow)")
            },
            ". This search simulates cluster sizes up to ", .CRD_ES_MAX_CLUSTER_SIZE,
            ": each simulated data set is fitted with a model that has one ",
            "variable per individual of a cluster, and the time a fit takes grows ",
            "faster than the cluster size. ",
            "ss_aipe_crd_es_n_clusters_fixed_width() plans the number of clusters ",
            "instead."),
            arg = "width", reason = "cluster_size_beyond_reach", start = startval,
            lower_bound = lower_bound, max_cluster_size = .CRD_ES_MAX_CLUSTER_SIZE,
            width_limit = width_limit, target = width, n_clusters = n_clusters)
    }
    # The upward search reaches .CRD_ES_CLUSTER_SIZE_REACH above diff_start
    # where that is finite, and above the start otherwise, and never above
    # .CRD_ES_MAX_CLUSTER_SIZE: every cluster size a search from diff_start
    # (the start through the first builds of 1.1.0 wherever it was at most
    # 100,000) could reach stays within reach of a search from a smaller
    # start. The size a target needs can lie well above the closed-form
    # start when the effect is large, since the standardizer is then
    # estimated less precisely relative to the effect.
    reach_from <- if (is.finite(diff_start)) diff_start else startval
    startwidth <- .find_width_crd_es(G, assurance = assurance, n_clusters = n_clusters, n_treat_clus = round(n_clusters * pr_treat), cluster_size = startval, true_icc = true_icc, es = es, es_type = es_type, true_total_var = true_total_var, covariate = as.logical(n_covariates), true_icc_covariate = true_icc_covariate, true_R2_within = true_R2_within, true_R2_between = true_R2_between, total_var_Z = 1, seed = seed, ncpus = ncpus, conf_level = conf_level, diff_size = diff_size, tracker = tracker)
    # A width meets the target when it is no larger than the target (see
    # .crd_es_search_up()); a start whose width meets it steps down.
    if (startwidth <= width) {
        repeat {
            startval <- startval - 1
            if (startval == 1) return(c(startval + 1, startwidth))
            savedwidth <- startwidth
            startwidth <- .find_width_crd_es(G, assurance = assurance, n_clusters = n_clusters, n_treat_clus = round(n_clusters * pr_treat), cluster_size = startval, true_icc = true_icc, es = es, es_type = es_type, true_total_var = true_total_var, covariate = as.logical(n_covariates), true_icc_covariate = true_icc_covariate, true_R2_within = true_R2_within, true_R2_between = true_R2_between, total_var_Z = 1, seed = seed, ncpus = ncpus, conf_level = conf_level, diff_size = diff_size, tracker = tracker)
            if (startwidth > width) return(c(startval + 1, savedwidth))
        }
    } else {
        evaluate <- function(m) {
            .find_width_crd_es(G, assurance = assurance, n_clusters = n_clusters, n_treat_clus = round(n_clusters * pr_treat), cluster_size = m, true_icc = true_icc, es = es, es_type = es_type, true_total_var = true_total_var, covariate = as.logical(n_covariates), true_icc_covariate = true_icc_covariate, true_R2_within = true_R2_within, true_R2_between = true_R2_between, total_var_Z = 1, seed = seed, ncpus = ncpus, conf_level = conf_level, diff_size = diff_size, tracker = tracker)
        }
        return(.crd_es_search_up(evaluate, start = startval, start_width = startwidth, target = width, tracker = tracker, searched = "cluster_size", fixed = n_clusters, max_size = min(reach_from + .CRD_ES_CLUSTER_SIZE_REACH, .CRD_ES_MAX_CLUSTER_SIZE), assurance = assurance, width_limit = width_limit, reach_from = reach_from))
    }
}

# Find the cluster size given the specified width and the number of
# clusters: the smallest cluster size, from 2, whose width is no larger
# than the target, or Inf when the width at 100,000 individuals per
# cluster, the largest size the search considers, is still above the
# target (a number, so that a vector of cluster sizes stays numeric; the
# public planners refuse it with reason "width_not_reachable"). A width
# meets the target when it is no larger than it. The callers pass a number
# of clusters with a cluster in each condition and a degree of freedom, so
# a width at 100,000 that is not a number is refused
# (.crd_refuse_width_nan()); `fn` names the public planner in that
# refusal. Through the first builds of 1.1.0 the case beyond 100,000
# returned the string "> 100000", which made the value column of
# ss_aipe_crd_both_fixed_width() character whenever one of the numbers of
# clusters it tried could not reach the target; a width had to be below the
# target, so a target equal to the width of a design was met by the next
# cluster size; and a width equal to the target at the end of a bracket
# stopped the search with "'from' must be a finite number".
.find_n_indiv_crd_diff <- function(width, n_clusters, pr_treat, true_var_between = NULL, true_var_within = NULL, true_total_var = NULL, true_icc = NULL, true_R2_between = 0, true_R2_within = 0, n_covariates = 0, assurance = NULL, conf_level = 0.95, diff_size = NULL, fn = NULL) {
    baremaximum <- .find_width_crd_diff(n_clusters = n_clusters, cluster_size = 100000, pr_treat = pr_treat, true_var_between = true_var_between, true_var_within = true_var_within, true_total_var = true_total_var, true_icc = true_icc, true_R2_between = true_R2_between, true_R2_within = true_R2_within, n_covariates = n_covariates, assurance = assurance, conf_level = conf_level, diff_size = diff_size)
    if (is.na(baremaximum)) .crd_refuse_width_nan(fn, n_clusters, 100000, diff_size)
    if (width < baremaximum) return(Inf)
    cluster_size <- seq(100, 1000, 100)
    result <- sapply(cluster_size, .find_width_crd_diff, n_clusters = n_clusters, pr_treat = pr_treat, true_var_between = true_var_between, true_var_within = true_var_within, true_total_var = true_total_var, true_icc = true_icc, true_R2_between = true_R2_between, true_R2_within = true_R2_within, n_covariates = n_covariates, assurance = assurance, conf_level = conf_level, diff_size = diff_size)
    if (all(width >= result)) {
        cluster_size <- seq(2, 100, 1)
        result <- sapply(cluster_size, .find_width_crd_diff, n_clusters = n_clusters, pr_treat = pr_treat, true_var_between = true_var_between, true_var_within = true_var_within, true_total_var = true_total_var, true_icc = true_icc, true_R2_between = true_R2_between, true_R2_within = true_R2_within, n_covariates = n_covariates, assurance = assurance, conf_level = conf_level, diff_size = diff_size)
        if (all(width >= result)) {
            return(2)
        } else if (all(width < result)) {
            return(100)
        } else {
            return(cluster_size[which(width >= result)[1]])
        }
    } else if (all(width < result)) {
        start <- 1000
        repeat {
            cluster_size <- seq(start, start + 1000, 100)
            result <- sapply(cluster_size, .find_width_crd_diff, n_clusters = n_clusters, pr_treat = pr_treat, true_var_between = true_var_between, true_var_within = true_var_within, true_total_var = true_total_var, true_icc = true_icc, true_R2_between = true_R2_between, true_R2_within = true_R2_within, n_covariates = n_covariates, assurance = assurance, conf_level = conf_level, diff_size = diff_size)
            if (all(width >= result)) {
                return(start)
            } else if (all(width < result)) {
                start <- start + 1000
            } else {
                minval <- cluster_size[which(width >= result)[1] - 1]
                maxval <- cluster_size[which(width >= result)[1]]
                cluster_size <- seq(minval, maxval, 1)
                result <- sapply(cluster_size, .find_width_crd_diff, n_clusters = n_clusters, pr_treat = pr_treat, true_var_between = true_var_between, true_var_within = true_var_within, true_total_var = true_total_var, true_icc = true_icc, true_R2_between = true_R2_between, true_R2_within = true_R2_within, n_covariates = n_covariates, assurance = assurance, conf_level = conf_level, diff_size = diff_size)
                if (all(width >= result)) {
                    return(minval)
                } else if (all(width < result)) {
                    return(maxval)
                } else {
                    return(cluster_size[which(width >= result)[1]])
                }
            }
        }
    } else {
        minval <- cluster_size[which(width >= result)[1] - 1]
        maxval <- cluster_size[which(width >= result)[1]]
        cluster_size <- seq(minval, maxval, 1)
        result <- sapply(cluster_size, .find_width_crd_diff, n_clusters = n_clusters, pr_treat = pr_treat, true_var_between = true_var_between, true_var_within = true_var_within, true_total_var = true_total_var, true_icc = true_icc, true_R2_between = true_R2_between, true_R2_within = true_R2_within, n_covariates = n_covariates, assurance = assurance, conf_level = conf_level, diff_size = diff_size)
        if (all(width >= result)) {
            return(minval)
        } else if (all(width < result)) {
            return(maxval)
        } else {
            return(cluster_size[which(width >= result)[1]])
        }
    }
}

# Find the number of clusters given the budget and the cluster size: the
# largest number of clusters, from 4, whose cost is within the budget, or
# NA when 4 clusters cost more than the budget (the public planners refuse
# that with reason "budget_infeasible"). A cost equal to the budget is
# within it. The costs must not both be 0, or every number of clusters is
# within the budget and the search has no end; the public planners refuse
# that before searching. Through the first builds of 1.1.0 a budget equal
# to the cost at the end of a bracket (such as 20000 for 1000 clusters of
# 20 at a cost of 1 each) stopped the search with "'from' must be a
# finite number".
.find_n_clus_crd_budget <- function(budget, cluster_size, cluster_cost, individual_cost, diff_size = NULL) {
    n_clusters <- seq(100, 1000, 100)
    result <- sapply(n_clusters, .cost_crd, cluster_size = cluster_size, cluster_cost = cluster_cost, individual_cost = individual_cost, diff_size = diff_size)
    if (all(budget < result)) {
        n_clusters <- seq(4, 100, 1)
        result <- sapply(n_clusters, .cost_crd, cluster_size = cluster_size, cluster_cost = cluster_cost, individual_cost = individual_cost, diff_size = diff_size)
        if (all(budget < result)) {
            return(NA)
        } else if (all(budget >= result)) {
            return(100)
        } else {
            index <- which(budget >= result)
            return(n_clusters[index[length(index)]])
        }
    } else if (all(budget >= result)) {
        start <- 1000
        repeat {
            n_clusters <- seq(start, start + 1000, 100)
            result <- sapply(n_clusters, .cost_crd, cluster_size = cluster_size, cluster_cost = cluster_cost, individual_cost = individual_cost, diff_size = diff_size)
            if (all(budget < result)) {
                return(start)
            } else if (all(budget >= result)) {
                start <- start + 1000
            } else {
                minval <- n_clusters[which(budget < result)[1] - 1]
                maxval <- n_clusters[which(budget < result)[1]]
                n_clusters <- seq(minval, maxval, 1)
                result <- sapply(n_clusters, .cost_crd, cluster_size = cluster_size, cluster_cost = cluster_cost, individual_cost = individual_cost, diff_size = diff_size)
                if (all(budget < result)) {
                    return(minval)
                } else if (all(budget >= result)) {
                    return(maxval)
                } else {
                    index <- which(budget >= result)
                    return(n_clusters[index[length(index)]])
                }
            }
        }
    } else {
        minval <- n_clusters[which(budget < result)[1] - 1]
        maxval <- n_clusters[which(budget < result)[1]]
        n_clusters <- seq(minval, maxval, 1)
        result <- sapply(n_clusters, .cost_crd, cluster_size = cluster_size, cluster_cost = cluster_cost, individual_cost = individual_cost, diff_size = diff_size)
        if (all(budget < result)) {
            return(minval)
        } else if (all(budget >= result)) {
            return(maxval)
        } else {
            index <- which(budget >= result)
            return(n_clusters[index[length(index)]])
        }
    }
}

# Find the cluster size given the budget and the number of clusters: the
# largest cluster size, from 2, whose cost is within the budget, or NA
# when clusters of 2 cost more than the budget. A cost equal to the budget
# is within it. 'individual_cost' must be positive, or every cluster size
# that one cluster costs allows is within the budget and the search has
# no end; the public planners refuse a cost of 0 before searching.
.find_n_indiv_crd_budget <- function(budget, n_clusters, cluster_cost, individual_cost, diff_size = NULL) {
    cluster_size <- seq(100, 1000, 100)
    result <- sapply(cluster_size, .cost_crd, n_clusters = n_clusters, cluster_cost = cluster_cost, individual_cost = individual_cost, diff_size = diff_size)
    if (all(budget < result)) {
        cluster_size <- seq(2, 100, 1)
        result <- sapply(cluster_size, .cost_crd, n_clusters = n_clusters, cluster_cost = cluster_cost, individual_cost = individual_cost, diff_size = diff_size)
        if (all(budget < result)) {
            return(NA)
        } else if (all(budget >= result)) {
            return(100)
        } else {
            index <- which(budget >= result)
            return(cluster_size[index[length(index)]])
        }
    } else if (all(budget >= result)) {
        start <- 1000
        repeat {
            cluster_size <- seq(start, start + 1000, 100)
            result <- sapply(cluster_size, .cost_crd, n_clusters = n_clusters, cluster_cost = cluster_cost, individual_cost = individual_cost, diff_size = diff_size)
            if (all(budget < result)) {
                return(start)
            } else if (all(budget >= result)) {
                start <- start + 1000
            } else {
                minval <- cluster_size[which(budget < result)[1] - 1]
                maxval <- cluster_size[which(budget < result)[1]]
                cluster_size <- seq(minval, maxval, 1)
                result <- sapply(cluster_size, .cost_crd, n_clusters = n_clusters, cluster_cost = cluster_cost, individual_cost = individual_cost, diff_size = diff_size)
                if (all(budget < result)) {
                    return(minval)
                } else if (all(budget >= result)) {
                    return(maxval)
                } else {
                    index <- which(budget >= result)
                    return(cluster_size[index[length(index)]])
                }
            }
        }
    } else {
        minval <- cluster_size[which(budget < result)[1] - 1]
        maxval <- cluster_size[which(budget < result)[1]]
        cluster_size <- seq(minval, maxval, 1)
        result <- sapply(cluster_size, .cost_crd, n_clusters = n_clusters, cluster_cost = cluster_cost, individual_cost = individual_cost, diff_size = diff_size)
        if (all(budget < result)) {
            return(minval)
        } else if (all(budget >= result)) {
            return(maxval)
        } else {
            index <- which(budget >= result)
            return(cluster_size[index[length(index)]])
        }
    }
}

# Find the least width of ES combination of the number of clusters and
# cluster size given the budget. The three neighboring numbers of clusters
# around the start are compared, and the search moves toward the narrowest
# until it is in the middle. which.min() picks the first of tied widths
# (through the first builds of 1.1.0 a tie, such as two designs at which no
# draw gave an interval, stopped with an error from if()), three designs
# with no finite width are refused with their counts, and the call's budget
# of designs (.CRD_ES_MAX_DESIGNS, checked by .find_width_crd_es()) bounds
# the search.
.find_min_width_crd_es <- function(budget, cluster_cost = 0, individual_cost = 1, es, es_type = 1, true_icc, pr_treat, true_R2_between = 0, true_R2_within = 0, n_covariates = 0, assurance = NULL, conf_level = 0.95, G = 1000, true_icc_covariate = NULL, seed = NULL, ncpus = 1, diff_size = NULL, tracker = NULL) {
    if (n_covariates > 0 && is.null(true_icc_covariate)) true_icc_covariate <- true_icc
    if (n_covariates > 1) .dmar_refuse_argument(
      "'n_covariates' must be 0 or 1: the planners for the standardized effect size allow at most one covariate beside the treatment indicator.",
      arg = "n_covariates")
    true_total_var <- 1
    if (es_type == 0) {
        true_total_var <- 1
    } else if (es_type == 1) {
        true_total_var <- 1 / (1 - true_icc)
    } else if (es_type == 2) {
        true_total_var <- 1 / true_icc
    } else {
        .dmar_refuse_argument(
          "'es_type' can be 0 (total variance), 1 (level-1 variance), or 2 (level-2 variance) only.",
          arg = "es_type")
    }
    # A number of clusters at which the budget does not buy clusters of 2
    # (the neighbor above the largest number it buys) is not a design: its
    # width is Inf and nothing is simulated. Through the first builds of
    # 1.1.0 it was simulated with a cluster size of NA, which stopped with
    # "invalid 'times' argument".
    FUN <- function(n_clusters, cluster_size) {
        if (is.na(cluster_size)) return(Inf)
        .find_width_crd_es(G = G, assurance = assurance, n_clusters = n_clusters, n_treat_clus = round(n_clusters * pr_treat), cluster_size = cluster_size, true_icc = true_icc, es = es, es_type = es_type, true_total_var = true_total_var, covariate = as.logical(n_covariates), true_icc_covariate = true_icc_covariate, true_R2_within = true_R2_within, true_R2_between = true_R2_between, total_var_Z = 1, seed = seed, ncpus = ncpus, conf_level = conf_level, diff_size = diff_size, tracker = tracker)
    }
    none_finite <- function(widths, sizes) {
        if (any(is.finite(vapply(widths, as.numeric, numeric(1))))) return(invisible(NULL))
        sims <- lapply(widths, function(w) attr(w, "simulation"))
        simulated <- vapply(sims, is.list, logical(1))
        totals <- .crd_search_totals(sims[simulated])
        counts <- if (any(simulated)) {
            c(" Data sets without a usable interval: ",
              paste0(vapply(sims[simulated], function(s) as.integer(s$failed), integer(1)), " of ",
                     vapply(sims[simulated], function(s) as.integer(s$attempted), integer(1)),
                     collapse = ", "),
              "; by reason, ",
              .crd_reasons_phrase(totals$failure_reasons, totals$failure_messages), ".")
        } else {
            character(0)
        }
        .dmar_refuse(c(
            .crd_es_fun(tracker),
            ": no usable width at ", paste(sizes[simulated], collapse = ", "),
            " clusters, the designs the budget buys around the starting value",
            if (any(!simulated)) c(" (it does not buy ", paste(sizes[!simulated], collapse = ", "),
                                   " clusters of 2)"),
            ": at each, no simulated data set gave a usable confidence interval",
            if (!is.null(assurance)) " often enough to reach the assurance" else "",
            ". The narrowest design for this budget cannot be found.", counts),
            "computation_failed", reason = "no_finite_width", n_clusters = sizes,
            attempted = totals$attempted, returned = totals$returned,
            failed = totals$failed, failure_reasons = totals$failure_reasons,
            failure_messages = totals$failure_messages)
    }
    startval <- .find_min_width_crd_diff(budget = budget, cluster_cost = cluster_cost, individual_cost = individual_cost, pr_treat = pr_treat, true_total_var = true_total_var, true_icc = true_icc, true_R2_between = true_R2_between, true_R2_within = true_R2_within, n_covariates = n_covariates, assurance = assurance, conf_level = conf_level, diff_size = diff_size, fn = tracker$fun)

    start_n_clusters <- c(startval[1] - 1, startval[1], startval[1] + 1)
    result_n_individuals <- sapply(start_n_clusters, .find_n_indiv_crd_budget, budget = budget, cluster_cost = cluster_cost, individual_cost = individual_cost, diff_size = diff_size)
    # The widths are kept in a list so that each keeps the record of its
    # simulation for a refusal; num() reads the numbers.
    startwidth <- mapply(FUN, n_clusters = start_n_clusters, cluster_size = result_n_individuals, SIMPLIFY = FALSE)
    num <- function(w) vapply(w, as.numeric, numeric(1))
    none_finite(startwidth, start_n_clusters)

    if (which.min(num(startwidth)) == 1) {
        repeat {
            start_n_clusters <- start_n_clusters - 1
            result_n_individuals[2:3] <- result_n_individuals[1:2]
            startwidth[2:3] <- startwidth[1:2]
            result_n_individuals[1] <- .find_n_indiv_crd_budget(start_n_clusters[1], budget = budget, cluster_cost = cluster_cost, individual_cost = individual_cost, diff_size = diff_size)
            startwidth[[1]] <- FUN(n_clusters = start_n_clusters[1], cluster_size = result_n_individuals[1])
            none_finite(startwidth, start_n_clusters)
            if (which.min(num(startwidth)) != 1) return(c(start_n_clusters[2], result_n_individuals[2], num(startwidth)[2]))
        }
    } else if (which.min(num(startwidth)) == 3) {
        repeat {
            start_n_clusters <- start_n_clusters + 1
            result_n_individuals[1:2] <- result_n_individuals[2:3]
            startwidth[1:2] <- startwidth[2:3]
            result_n_individuals[3] <- .find_n_indiv_crd_budget(start_n_clusters[3], budget = budget, cluster_cost = cluster_cost, individual_cost = individual_cost, diff_size = diff_size)
            startwidth[[3]] <- FUN(n_clusters = start_n_clusters[3], cluster_size = result_n_individuals[3])
            none_finite(startwidth, start_n_clusters)
            if (which.min(num(startwidth)) != 3) return(c(start_n_clusters[2], result_n_individuals[2], num(startwidth)[2]))
        }
    } else {
        return(c(start_n_clusters[2], result_n_individuals[2], num(startwidth)[2]))
    }
}

# Find the least width combination of the number of clusters and cluster
# size given the budget. At each number of clusters the cluster size is
# the largest the budget buys (.find_n_indiv_crd_budget()), and the widths
# of those designs are compared: at 100, 200, ..., 1100 clusters (in
# blocks above that while the narrowest is the largest), and then one
# cluster at a time around the narrowest, or from the floor
# .crd_min_n_clusters() to 200 when the narrowest is at 100 or no design of
# 100 or more clusters is within the budget. A number of clusters below the
# floor, or one at which the budget does not buy clusters of 2, is not a
# candidate. The public planners refuse a budget that does not buy the
# floor's number of clusters of 2 (reason "budget_infeasible") before
# searching, so a search without a candidate (class
# dmar_computation_failed, reason "no_candidate") is a safeguard, as is
# the refusal of a width that is not a number at a candidate
# (.crd_refuse_width_nan()); `fn` names the public planner in either.
# Returns c(n_clusters, cluster_size, width). Through the first builds of
# 1.1.0 a search whose numbers of clusters had no candidate stopped with
# "invalid 'type' (list) of argument" (budget = 20000, cluster_cost =
# 1000, individual_cost = 1, pr_treat = .1, where 20 clusters of 2 cost
# 20,040), a floor above 200 (pr_treat below .01) stopped with "wrong sign
# in 'by' argument", and a number of clusters below a floor above 100 was
# a candidate.
.find_min_width_crd_diff <- function(budget, cluster_cost = 0, individual_cost = 1, pr_treat, true_var_between = NULL, true_var_within = NULL, true_total_var = NULL, true_icc = NULL, true_R2_between = 0, true_R2_within = 0, n_covariates = 0, assurance = NULL, conf_level = 0.95, diff_size = NULL, fn = NULL) {
    lo <- .crd_min_n_clusters(pr_treat, n_covariates)
    # The cluster size the budget buys at each number of clusters (NA where
    # it buys none) and the width of that design (NA where the number of
    # clusters is not a candidate).
    evaluate <- function(n_clusters) {
        cluster_size <- vapply(n_clusters, .find_n_indiv_crd_budget, numeric(1), budget = budget, cluster_cost = cluster_cost, individual_cost = individual_cost, diff_size = diff_size)
        width <- rep(NA_real_, length(n_clusters))
        use <- which(!is.na(cluster_size) & n_clusters >= lo)
        width[use] <- vapply(use, function(i) .find_width_crd_diff(n_clusters = n_clusters[i], cluster_size = cluster_size[i], pr_treat = pr_treat, true_var_between = true_var_between, true_var_within = true_var_within, true_total_var = true_total_var, true_icc = true_icc, true_R2_between = true_R2_between, true_R2_within = true_R2_within, n_covariates = n_covariates, assurance = assurance, conf_level = conf_level, diff_size = diff_size), numeric(1))
        nan <- use[is.na(width[use])]
        if (length(nan)) .crd_refuse_width_nan(fn, n_clusters[nan[1]], cluster_size[nan[1]], diff_size)
        list(n_clusters = n_clusters, cluster_size = cluster_size, width = width)
    }
    # The narrowest candidate, the first of tied widths.
    narrowest <- function(r) {
        if (!any(!is.na(r$width))) .crd_refuse_no_candidate(fn, lo, pr_treat, n_covariates, range(r$n_clusters), budget = budget)
        i <- which.min(r$width)
        c(r$n_clusters[i], r$cluster_size[i], r$width[i])
    }
    r <- evaluate(seq(100, 1100, 100))
    posmin <- if (any(!is.na(r$width))) which.min(r$width) else 1L
    if (posmin == 1) {
        return(narrowest(evaluate(seq.int(lo, max(lo, 200)))))
    } else if (posmin == length(r$width)) {
        start <- r$n_clusters[length(r$n_clusters) - 1]
        repeat {
            r <- evaluate(seq(start, start + 1100, 100))
            if (!any(!is.na(r$width))) .crd_refuse_no_candidate(fn, lo, pr_treat, n_covariates, range(r$n_clusters), budget = budget)
            posmin <- which.min(r$width)
            if (posmin == 1) {
                return(narrowest(r))
            } else if (posmin == length(r$width)) {
                start <- start + 1000
            } else {
                return(narrowest(evaluate(seq(r$n_clusters[posmin - 1], r$n_clusters[posmin + 1], 1))))
            }
        }
    } else {
        return(narrowest(evaluate(seq(r$n_clusters[posmin - 1], r$n_clusters[posmin + 1], 1))))
    }
}

# Find the least expensive combination of the number of clusters and
# cluster size given the specified width of ES. For each of three
# neighboring cluster sizes the number of clusters that reaches the width is
# found (.find_n_clus_crd_es()), and the search moves toward the least
# costly until it is in the middle. The searches it runs share the call's
# budget of designs (.CRD_ES_MAX_DESIGNS, checked by .find_width_crd_es()),
# which bounds the whole call; through the first build of this repair each
# of up to 100 steps ran a search with a budget of its own.
.find_min_cost_crd_es <- function(width, cluster_cost = 0, individual_cost = 1, es, es_type = 1, true_icc, pr_treat, true_R2_between = 0, true_R2_within = 0, n_covariates = 0, assurance = NULL, conf_level = 0.95, G = 1000, true_icc_covariate = NULL, seed = NULL, ncpus = 1, diff_size = NULL, tracker = NULL) {
    if (n_covariates > 0 && is.null(true_icc_covariate)) true_icc_covariate <- true_icc
    if (n_covariates > 1) .dmar_refuse_argument(
      "'n_covariates' must be 0 or 1: the planners for the standardized effect size allow at most one covariate beside the treatment indicator.",
      arg = "n_covariates")
    true_total_var <- 1
    if (es_type == 0) {
        true_total_var <- 1
    } else if (es_type == 1) {
        true_total_var <- 1 / (1 - true_icc)
    } else if (es_type == 2) {
        true_total_var <- 1 / true_icc
    } else {
        .dmar_refuse_argument(
          "'es_type' can be 0 (total variance), 1 (level-1 variance), or 2 (level-2 variance) only.",
          arg = "es_type")
    }
    startval <- .find_min_cost_crd_diff(width = width, cluster_cost = cluster_cost, individual_cost = individual_cost, pr_treat = pr_treat, true_total_var = true_total_var, true_icc = true_icc, true_R2_between = true_R2_between, true_R2_within = true_R2_within, n_covariates = n_covariates, assurance = assurance, conf_level = conf_level, diff_size = diff_size, fn = tracker$fun)
    startval <- as.numeric(startval)
    start_n_individuals <- c(startval[2] - 1, startval[2], startval[2] + 1)
    result <- sapply(start_n_individuals, .find_n_clus_crd_es, width = width, es = es, es_type = es_type, true_icc = true_icc, pr_treat = pr_treat, true_R2_between = true_R2_between, true_R2_within = true_R2_within, n_covariates = n_covariates, assurance = assurance, conf_level = conf_level, G = G, true_icc_covariate = true_icc_covariate, seed = seed, ncpus = ncpus, diff_size = diff_size, tracker = tracker)
    result_n_clusters <- result[1, ]
    resultwidth <- result[2, ]
    startbudget <- mapply(.cost_crd, n_clusters = result_n_clusters, cluster_size = start_n_individuals, MoreArgs = list(cluster_cost = cluster_cost, individual_cost = individual_cost, diff_size = diff_size))
    if (which.min(startbudget) == 1) {
        repeat {
            start_n_individuals <- start_n_individuals - 1
            result_n_clusters[2:3] <- result_n_clusters[1:2]
            startbudget[2:3] <- startbudget[1:2]
            resultwidth[2:3] <- resultwidth[1:2]
            result <- .find_n_clus_crd_es(width = width, cluster_size = start_n_individuals[1], es = es, es_type = es_type, true_icc = true_icc, pr_treat = pr_treat, true_R2_between = true_R2_between, true_R2_within = true_R2_within, n_covariates = n_covariates, assurance = assurance, conf_level = conf_level, G = G, true_icc_covariate = true_icc_covariate, seed = seed, ncpus = ncpus, diff_size = diff_size, tracker = tracker)
            result_n_clusters[1] <- result[1]
            resultwidth[1] <- result[2]
            startbudget[1] <- .cost_crd(n_clusters = result_n_clusters[1], cluster_size = start_n_individuals[1], cluster_cost = cluster_cost, individual_cost = individual_cost, diff_size = diff_size)
            if (which.min(startbudget) != 1) return(c(result_n_clusters[2], start_n_individuals[2], startbudget[2], resultwidth[2]))
        }
    } else if (which.min(startbudget) == 3) {
        repeat {
            start_n_individuals <- start_n_individuals + 1
            result_n_clusters[1:2] <- result_n_clusters[2:3]
            startbudget[1:2] <- startbudget[2:3]
            resultwidth[1:2] <- resultwidth[2:3]
            result <- .find_n_clus_crd_es(width = width, cluster_size = start_n_individuals[3], es = es, es_type = es_type, true_icc = true_icc, pr_treat = pr_treat, true_R2_between = true_R2_between, true_R2_within = true_R2_within, n_covariates = n_covariates, assurance = assurance, conf_level = conf_level, G = G, true_icc_covariate = true_icc_covariate, seed = seed, ncpus = ncpus, diff_size = diff_size, tracker = tracker)
            result_n_clusters[3] <- result[1]
            resultwidth[3] <- result[2]
            startbudget[3] <- .cost_crd(n_clusters = result_n_clusters[3], cluster_size = start_n_individuals[3], cluster_cost = cluster_cost, individual_cost = individual_cost, diff_size = diff_size)
            if (which.min(startbudget) != 3) return(c(result_n_clusters[2], start_n_individuals[2], startbudget[2], resultwidth[2]))
        }
    } else {
        return(c(result_n_clusters[2], start_n_individuals[2], startbudget[2], resultwidth[2]))
    }
}

# Find the least expensive combination of the number of clusters and
# cluster size given the specified width. At each number of clusters the
# cluster size is the smallest that reaches the width
# (.find_n_indiv_crd_diff(), Inf where no cluster size up to 100,000 does),
# and the costs of those designs are compared: at 100, 200, ..., 1100
# clusters (moved up by 1000 until one reaches the width, and in blocks
# above that while the least costly is the largest), and then one cluster
# at a time around the least costly, or from the floor .crd_min_n_clusters()
# to 200 when the least costly is at 100. A number of clusters below the
# floor is not a candidate, and its cluster size is not searched for.
# 'individual_cost' is positive (the public planners refuse 0), so a
# design that does not reach the width costs Inf. When no number of
# clusters from 100 to 1100 reaches the width with clusters of up to
# 100,000, the cluster size at .CRD_MAX_N_CLUSTERS is found first, and a
# target no design within both limits reaches is refused
# (.crd_refuse_n_clusters_beyond()); otherwise the blocks of 1000 end by
# that number. `fn` names the public planner in a refusal. Returns
# c(n_clusters, cluster_size, cost), all numbers. Through the first builds
# of 1.1.0 the cluster size of a number of clusters that could not reach
# the width was the string "> 100000", which made the result, and the
# value column of ss_aipe_crd_both_fixed_width(), character; a floor above
# 200 (pr_treat below .01) stopped with "wrong sign in 'by' argument"; a
# number of clusters below a floor above 100 was a candidate; and the
# blocks of 1000 had no end.
.find_min_cost_crd_diff <- function(width, cluster_cost = 0, individual_cost = 1, pr_treat, true_var_between = NULL, true_var_within = NULL, true_total_var = NULL, true_icc = NULL, true_R2_between = 0, true_R2_within = 0, n_covariates = 0, assurance = NULL, conf_level = 0.95, diff_size = NULL, fn = NULL) {
    lo <- .crd_min_n_clusters(pr_treat, n_covariates)
    size_at <- function(n_clusters) .find_n_indiv_crd_diff(width = width, n_clusters = n_clusters, pr_treat = pr_treat, true_var_between = true_var_between, true_var_within = true_var_within, true_total_var = true_total_var, true_icc = true_icc, true_R2_between = true_R2_between, true_R2_within = true_R2_within, n_covariates = n_covariates, assurance = assurance, conf_level = conf_level, diff_size = diff_size, fn = fn)
    # The cluster size that reaches the width at each number of clusters
    # and the cost of that design (NA where the number of clusters is not a
    # candidate).
    evaluate <- function(n_clusters) {
        cluster_size <- rep(NA_real_, length(n_clusters))
        cost <- rep(NA_real_, length(n_clusters))
        use <- which(n_clusters >= lo)
        cluster_size[use] <- vapply(n_clusters[use], size_at, numeric(1))
        cost[use] <- vapply(use, function(i) .cost_crd(n_clusters = n_clusters[i], cluster_size = cluster_size[i], cluster_cost = cluster_cost, individual_cost = individual_cost, diff_size = diff_size), numeric(1))
        list(n_clusters = n_clusters, cluster_size = cluster_size, cost = cost)
    }
    # The least costly candidate that reaches the width, the first of tied
    # costs.
    cheapest <- function(r) {
        if (!any(is.finite(r$cost))) .crd_refuse_no_candidate(fn, lo, pr_treat, n_covariates, range(r$n_clusters), width = width)
        i <- which.min(r$cost)
        c(r$n_clusters[i], r$cluster_size[i], r$cost[i])
    }
    r <- evaluate(seq(100, 1100, 100))
    if (!any(is.finite(r$cost)) && !is.finite(size_at(max(lo, .CRD_MAX_N_CLUSTERS)))) {
        at_max <- .find_width_crd_diff(n_clusters = max(lo, .CRD_MAX_N_CLUSTERS), cluster_size = 100000, pr_treat = pr_treat, true_var_between = true_var_between, true_var_within = true_var_within, true_total_var = true_total_var, true_icc = true_icc, true_R2_between = true_R2_between, true_R2_within = true_R2_within, n_covariates = n_covariates, assurance = assurance, conf_level = conf_level, diff_size = diff_size)
        .crd_refuse_n_clusters_beyond(fn, width, at_max, cluster_size = NULL, assurance = assurance)
    }
    while (!any(is.finite(r$cost))) r <- evaluate(r$n_clusters + 1000)
    posmin <- which.min(r$cost)
    if (posmin == 1) {
        return(cheapest(evaluate(seq.int(lo, max(lo, 200)))))
    } else if (posmin == length(r$cost)) {
        start <- r$n_clusters[length(r$n_clusters) - 1]
        repeat {
            r <- evaluate(seq(start, start + 1100, 100))
            posmin <- which.min(r$cost)
            if (posmin == 1) {
                return(cheapest(r))
            } else if (posmin == length(r$cost)) {
                start <- start + 1000
            } else {
                return(cheapest(evaluate(seq(r$n_clusters[posmin - 1], r$n_clusters[posmin + 1], 1))))
            }
        }
    } else {
        return(cheapest(evaluate(seq(r$n_clusters[posmin - 1], r$n_clusters[posmin + 1], 1))))
    }
}

# .crd_es_check_true_es(true_es, fun)
#
# The check of `true_es` in a standardized cluster randomized planner that
# computes a width: a single finite number. `true_es` is NULL when the
# caller's argument was missing or NULL. Anything else is refused with a
# condition of class dmar_argument_invalid that names the calling function
# `fun`, says why the effect size is needed, and says what to give. Returns
# the value.
.crd_es_check_true_es <- function(true_es, fun) {
  if (is.numeric(true_es) && length(true_es) == 1L && is.finite(true_es)) {
    return(true_es)
  }
  given <- if (is.null(true_es)) {
    "none was given"
  } else if (is.numeric(true_es) && length(true_es) == 1L) {
    paste0(format(true_es), " was given")
  } else {
    paste0("an object of class ", class(true_es)[1L], " and length ",
           length(true_es), " was given")
  }
  .dmar_refuse_argument(c(
    fun, "() needs 'true_es', the population standardized difference ",
    "between the condition means that the interval is planned for, as a ",
    "single finite number; ", given, ". The width of the interval is read ",
    "from data simulated at that effect size, so no width can be computed ",
    "without it. Give a planning value on the scale 'es_type' names, such ",
    "as true_es = 0.5."), arg = "true_es", value = true_es)
}

# .crd_es_check_width_request(G, pr_treat, true_icc, true_es, fun)
#
# The two standardized planners that fix a size and a budget return the
# size the budget buys, and add the expected width of its interval when
# `G` is supplied. The width is simulated, which takes the proportion of
# treatment clusters, the intraclass correlation, and the effect size. A
# call with `G` and without any of the three asked for a width that cannot
# be computed, and is refused with a condition of class
# dmar_argument_invalid that names the ones missing (through DMAR 1.0.x
# and the first builds of 1.1.0 a missing `pr_treat` or `true_icc` returned
# the table without the width and without a word, and a missing `true_es`
# stopped inside the simulation). With G = NULL nothing is checked.
.crd_es_check_width_request <- function(G, pr_treat, true_icc, true_es, fun) {
  if (is.null(G)) return(invisible(NULL))
  absent <- c(pr_treat = is.null(pr_treat), true_icc = is.null(true_icc),
              true_es = is.null(true_es))
  if (any(absent)) {
    lacking <- names(absent)[absent]
    what <- c(
      pr_treat = "'pr_treat' (the proportion of treatment clusters)",
      true_icc = paste0("'true_icc' (the population intraclass correlation ",
                        "of the outcome)"),
      true_es = paste0("'true_es' (the population standardized difference ",
                       "between the condition means)"))[lacking]
    listed <- if (length(what) == 1L) what else
      paste0(paste(what[-length(what)], collapse = ", "),
             if (length(what) > 2L) "," else "", " and ", what[length(what)])
    .dmar_refuse_argument(c(
      fun, "() was asked for the expected width of the interval (G = ",
      format(G), " replications), but ", listed,
      if (length(what) == 1L) " was" else " were", " not given. The width ",
      "is read from data simulated at those population values, so it ",
      "cannot be computed without them. Supply ",
      if (length(what) == 1L) "it" else "them",
      " (for example pr_treat = 0.5, true_icc = 0.25, true_es = 0.5), or ",
      "leave G = NULL to get the size the budget buys without an expected ",
      "width."), arg = lacking)
  }
  .crd_es_check_true_es(true_es, fun)
  invisible(NULL)
}

# .crd_es_check_covariate(n_covariates, true_icc_covariate, true_R2_within,
#                         true_R2_between)
#
# The check of the covariate's planning values in the two standardized
# planners that fix a size and a budget, made when the expected width is
# asked for (`G` supplied): `n_covariates` is 0 or 1, and with one
# covariate its intraclass correlation (when given) and the two
# proportions of variance it explains are single numbers in [0, 1], the
# proportions below 1. Through DMAR 1.0.x the two proportions had no
# default here and a call with a covariate and without them stopped inside
# the simulation with an error from R ("missing value where TRUE/FALSE
# needed", "arguments imply differing number of rows").
.crd_es_check_covariate <- function(n_covariates, true_icc_covariate, true_R2_within, true_R2_between) {
  if (!is.numeric(n_covariates) || length(n_covariates) != 1L ||
      is.na(n_covariates) || !(n_covariates %in% c(0, 1))) {
    .dmar_refuse_argument(c(
      "'n_covariates' must be 0 or 1: the planners for the standardized ",
      "effect size allow at most one covariate beside the treatment ",
      "indicator."), arg = "n_covariates")
  }
  if (n_covariates == 0) return(invisible(NULL))
  one_in <- function(v, upper_open) {
    is.numeric(v) && length(v) == 1L && is.finite(v) && v >= 0 &&
      (if (upper_open) v < 1 else v <= 1)
  }
  if (!is.null(true_icc_covariate) && !one_in(true_icc_covariate, FALSE)) {
    .dmar_refuse_argument(c(
      "'true_icc_covariate', the population intraclass correlation of the ",
      "covariate, must be a single number in [0, 1], or NULL to use the ",
      "value of 'true_icc'. 0 is a covariate that varies within clusters ",
      "only and 1 one that varies between clusters only."),
      arg = "true_icc_covariate", value = true_icc_covariate)
  }
  for (nm in c("true_R2_within", "true_R2_between")) {
    v <- get(nm)
    if (!one_in(v, TRUE)) {
      .dmar_refuse_argument(c(
        "'", nm, "' must be a single number in [0, 1): the proportion of ",
        "the ", if (nm == "true_R2_within") "within" else "between",
        "-level variance of the outcome that the covariate explains in ",
        "the population. With n_covariates = 1 the expected width is ",
        "simulated at this value, so it cannot be NULL or missing; the ",
        "default 0 is a covariate that explains nothing at that level. ",
        "Give, for example, ", nm, " = 0.3."),
        arg = nm, value = v)
    }
  }
  invisible(NULL)
}

# .crd_min_n_clusters(pr_treat, n_covariates)
#
# The fewest clusters the searches over the number of clusters consider:
# ceiling(2 / min(pr_treat, 1 - pr_treat)), at least two clusters in each
# condition, plus one for each covariate. The rule is MBESS's ("at least
# two groups in each condition", widthsscrd.R), and the searches of DMAR
# 1.0.x and of the first builds of 1.1.0 started there. With that many
# clusters each condition has at least two after rounding, and the
# interval has at least two degrees of freedom. It is written as MBESS
# writes it, 2 * (1 / p), so that the floor is the same number in floating
# point: 1 - .9 is a little below .1, and the floor is 21 at pr_treat = .9
# where it is 20 at .1.
.crd_min_n_clusters <- function(pr_treat, n_covariates = 0) {
    ceiling(2 * (1 / min(pr_treat, 1 - pr_treat))) + n_covariates
}

# .crd_first_meeting(from, meets)
#
# The first number of clusters from `from` upward at which meets(), which
# takes a vector of numbers of clusters, is TRUE, looked for in blocks of
# 100. The width falls toward 0 as clusters are added, so the search ends.
.crd_first_meeting <- function(from, meets) {
    repeat {
        n_clusters <- seq.int(from, from + 99)
        ok <- meets(n_clusters)
        if (any(ok)) return(n_clusters[which(ok)[1]])
        from <- from + 100
    }
}

# A number as the messages of the cluster randomized planners show it:
# with a thousands separator and without scientific notation.
.crd_num <- function(x) {
    format(x, big.mark = ",", scientific = FALSE, trim = TRUE, digits = 15)
}

# The planner a refusal from a search names: "fn()", or "The planner" when
# a search is called on its own.
.crd_fn_text <- function(fn) {
    if (is.null(fn)) "The planner" else paste0(fn, "()")
}

# .CRD_MAX_N_CLUSTERS
#
# The most clusters the searches over the number of clusters of the
# unstandardized difference consider (.find_n_clus_crd_diff(),
# .find_min_cost_crd_diff()), and so the most from which a standardized
# search starts. The width falls toward 0 as clusters are added, so every
# positive target is met by some number, but the searches step through
# blocks of 1000, and a target too narrow for any study (width = 1e-300)
# had them search without end; with unequal cluster sizes each width at a
# million clusters also builds a vector of a million sizes. The bound is
# the one the power planners use (.dmar_refuse_no_sample_size()). A target
# not met at it is refused before the blocks are searched
# (.crd_refuse_n_clusters_beyond()).
.CRD_MAX_N_CLUSTERS <- 1000000

# .crd_refuse_width_nan(fn, n_clusters, cluster_size, diff_size)
#
# The refusal of a search that computed a width that is not a number at a
# design with a cluster in each condition and a degree of freedom: class
# dmar_computation_failed, reason "width_not_a_number". The public
# planners check the variances, the intraclass correlation, and
# 'diff_size' on entry, so this is a safeguard. Through the first build of
# this repair such a width counted as one that does not meet the target,
# and the searches went on without end (true_icc = 1 or 1.5 in
# ss_aipe_crd_es_n_clusters_fixed_width(), diff_size = numeric(0)); before
# that it stopped them with "missing value where TRUE/FALSE needed".
.crd_refuse_width_nan <- function(fn, n_clusters, cluster_size, diff_size = NULL) {
    .dmar_refuse(c(
        .crd_fn_text(fn), ": the width of the interval for the difference between ",
        "the condition means at ", .crd_num(n_clusters), " clusters of ",
        .crd_num(cluster_size), " individuals",
        if (!is.null(diff_size)) " (before 'diff_size' is applied)",
        " is not a number, so the search cannot compare it with the target. The ",
        "variances it is computed from are not finite numbers of 0 or more at that ",
        "design. Check the planning values: the variances or the intraclass ",
        "correlation, and 'diff_size'."),
        "computation_failed", reason = "width_not_a_number", n_clusters = n_clusters,
        cluster_size = cluster_size)
}

# .crd_refuse_n_clusters_beyond(fn, width, width_at_max, cluster_size,
#                               assurance)
#
# The refusal of a target that the unstandardized difference does not meet
# at .CRD_MAX_N_CLUSTERS clusters (of `cluster_size`, or, when it is NULL,
# of 100,000, the largest cluster size the searches over both sizes
# consider): class dmar_argument_invalid, arg "width", reason
# "n_clusters_beyond_reach", with `width_at_max` and the bound as
# `max_n_clusters`.
.crd_refuse_n_clusters_beyond <- function(fn, width, width_at_max, cluster_size = NULL, assurance = NULL) {
    kind <- if (is.null(assurance)) "expected" else "assured"
    design <- if (is.null(cluster_size)) {
        "clusters of 100,000 individuals, the largest cluster size the search considers,"
    } else {
        paste0("clusters of ", .crd_num(cluster_size), " individuals")
    }
    .dmar_refuse_argument(c(
        .crd_fn_text(fn), ": no design with up to ", .crd_num(.CRD_MAX_N_CLUSTERS),
        " clusters, the most the search considers, meets the target 'width' of ",
        format(width), ". With ", .crd_num(.CRD_MAX_N_CLUSTERS), " ", design, " the ",
        kind, " width of the interval for the difference between the condition means ",
        "(unstandardized; the standardized planners start from it) is ",
        format(signif(width_at_max, 4)), ". The width falls with the square root of ",
        "the number of clusters, so halving it takes about four times as many. Give ",
        "a wider target."),
        arg = "width", reason = "n_clusters_beyond_reach",
        max_n_clusters = .CRD_MAX_N_CLUSTERS, width_at_max = width_at_max, target = width)
}

# .crd_check_diff_size(diff_size)
#
# The entry check of 'diff_size' in the twelve planners: NULL, or a
# numeric vector of one or more finite numbers. Through the first builds of
# 1.1.0 it was not checked: an empty vector or one holding NA or Inf
# stopped the searches with "missing value where TRUE/FALSE needed" (the
# first build of this repair searched without end there), and a character
# string with "non-numeric argument to mathematical function".
.crd_check_diff_size <- function(diff_size) {
    if (is.null(diff_size)) return(invisible(NULL))
    if (is.numeric(diff_size) && length(diff_size) >= 1L && all(is.finite(diff_size))) {
        return(invisible(diff_size))
    }
    given <- if (is.numeric(diff_size) && length(diff_size) > 1L) {
        paste0("a numeric vector holding ", sum(!is.finite(diff_size)), " value",
               if (sum(!is.finite(diff_size)) > 1) "s", " that ",
               if (sum(!is.finite(diff_size)) > 1) "are" else "is", " not finite was given")
    } else {
        .dmar_given(diff_size)
    }
    .dmar_refuse_argument(c(
        "'diff_size' must be NULL or a vector of one or more finite numbers; ",
        given, ". It sets the size of each cluster from 'cluster_size', recycled ",
        "over the clusters: whole numbers are added to it, and values that are ",
        "not all whole multiply it. Give, for example, diff_size = c(-2, 0, 2), or ",
        "NULL for clusters of equal size."),
        arg = "diff_size", value = diff_size)
}

# .crd_es_check_icc(true_icc, es_type)
#
# The entry check, in the standardized planners wherever a width is
# simulated, of `es_type` (0, 1, or 2, as the internals require) and
# `true_icc`: a single number in [0, 1), and above 0 with es_type = 2. The
# effect of es_type = 1 is divided by the within-cluster standard deviation,
# sqrt(1 - true_icc) of the total, and that of es_type = 2 by the
# between-cluster one, sqrt(true_icc); at an intraclass correlation of 1 no
# variance lies within clusters, every simulated cluster is one value
# repeated, and the model fitted to each data set has no within-cluster
# variance to estimate (with es_type = 0 every fit failed). Class
# dmar_argument_invalid. Through the first builds of 1.1.0 none of this was
# checked in ss_aipe_crd_es_n_clusters_fixed_width() and
# ss_aipe_crd_es_both_fixed_width(), and the range nowhere: true_icc = 1
# (es_type 1), 0 (es_type 2), 1.5, or -0.2 stopped with "missing value where
# TRUE/FALSE needed" (the first build of this repair searched without end),
# and true_icc = 1 with es_type = 0 simulated three designs at which every
# fit failed.
.crd_es_check_icc <- function(true_icc, es_type) {
    if (!(is.numeric(es_type) && length(es_type) == 1L && !is.na(es_type) &&
          es_type %in% c(0, 1, 2))) {
        .dmar_refuse_argument(
          "'es_type' can be 0 (total variance), 1 (level-1 variance), or 2 (level-2 variance) only.",
          arg = "es_type")
    }
    ok <- is.numeric(true_icc) && length(true_icc) == 1L && is.finite(true_icc) &&
        true_icc >= 0 && true_icc < 1
    if (!ok) {
        .dmar_refuse_argument(c(
            "'true_icc' must be a single number in [0, 1); ",
            if (is.null(true_icc)) "none was given" else .dmar_given(true_icc),
            ". It is the population intraclass correlation of the outcome, the share ",
            "of its variance that lies between clusters. At 1 no variance lies within ",
            "clusters: every simulated cluster is one value repeated, and the model ",
            "fitted to each data set has no within-cluster variance to estimate. Give, ",
            "for example, true_icc = 0.2."),
            arg = "true_icc", value = true_icc)
    }
    if (es_type == 2 && true_icc == 0) {
        .dmar_refuse_argument(c(
            "With es_type = 2 the effect is in units of the between-cluster standard ",
            "deviation, and true_icc = 0 makes that deviation 0, so the effect has no ",
            "scale. Give a positive true_icc, such as 0.2, or es_type = 0 (total ",
            "standard deviation) or 1 (within-cluster standard deviation)."),
            arg = c("true_icc", "es_type"), value = true_icc)
    }
    invisible(NULL)
}

# .crd_fixed_size_fewest(pr_treat, n_covariates)
#
# The fewest clusters of a fixed size whose interval has a finite width,
# for the planners that buy clusters of a fixed size with a budget and
# report the width: at least 4 (where the search starts), at least
# 3 + n_covariates (a degree of freedom), and a cluster in each condition
# at pr_treat. Returns list(n, why), `why` the phrase a refusal gives for
# a number above 4 (character(0) at 4).
.crd_fixed_size_fewest <- function(pr_treat, n_covariates) {
    from <- max(4, 3 + n_covariates)
    n <- .crd_fewest_both_arms(pr_treat, from)
    why <- c(if (n > from) paste0("the fewest that put a cluster in each condition at pr_treat = ", .crd_num(pr_treat)),
             if (from > 4) paste0("at least 3 + n_covariates = ", from, " for a degree of freedom"))
    list(n = n, why = if (length(why)) paste(why, collapse = "; ") else character(0))
}

# TRUE where a design of n_clusters clusters, round(n_clusters * pr_treat)
# of them treated, has a cluster in each condition.
.crd_both_arms <- function(n_clusters, pr_treat) {
    n_treat <- round(n_clusters * pr_treat)
    n_treat >= 1 & n_clusters - n_treat >= 1
}

# The fewest clusters, from `from`, that put a cluster in each condition
# at pr_treat. Once both conditions have a cluster they keep one as
# clusters are added, since each count grows by at most one at a time.
.crd_fewest_both_arms <- function(pr_treat, from) {
    n_clusters <- max(from, floor(0.5 / min(pr_treat, 1 - pr_treat)))
    while (!.crd_both_arms(n_clusters, pr_treat)) n_clusters <- n_clusters + 1
    n_clusters
}

# .crd_check_design(n_covariates, true_R2_between, true_R2_within,
#                   pr_treat = NULL, cluster_size = NULL, n_clusters = NULL)
#
# The entry check of the design arguments of the twelve cluster randomized
# planners, made once per call before any search. `n_covariates` is a
# single whole number of 0 or more; `true_R2_between` and `true_R2_within`
# are single numbers in [0, 1), as in .crd_es_check_covariate(), and
# `true_R2_between` is 0 when `n_covariates` is 0 (reason
# "R2_between_without_covariates": the share of the between-cluster
# variance the cluster-level covariates explain has no meaning without
# them, while `true_R2_within` keeps one); `pr_treat`, when given, is a
# single number strictly between 0 and 1;
# `cluster_size`, when given, is a single whole number of at least 2 (the
# cluster size searches start at 2); and `n_clusters`, when given, is a
# single whole number of at least 3 + n_covariates, so that the interval
# for the unstandardized difference has a degree of freedom. With
# `n_clusters` and `pr_treat` both given and `check_arms` TRUE (the
# planners that compute a width; the standardized budget planner with a
# fixed number of clusters passes it only when `G` asks for the width), the
# design must put a cluster in each condition (reason
# "arm_without_clusters"). Each refusal has class
# dmar_argument_invalid and names the argument. Through the first builds
# of 1.1.0 none of these was checked: n_covariates = .5 planned 17.5
# clusters, true_R2_between = 1.1 planned 5, a cluster size of 20.5 or 60.5
# clusters was planned for, 2 clusters stopped with "missing value where
# TRUE/FALSE needed" or returned a width of NaN, pr_treat = 0 and 1 searched
# without end, and 5 clusters at pr_treat = .1 (none treated) returned an
# infinite width. Through DMAR 1.0.x and the first builds of 1.1.0 a
# nonzero true_R2_between with n_covariates = 0 shrank the unstandardized
# plan (true_R2_between = .9 planned 7 clusters where 0 planned 18) and the
# start of a standardized search, whose simulation ignored it.
.crd_check_design <- function(n_covariates = 0, true_R2_between = 0, true_R2_within = 0, pr_treat = NULL, cluster_size = NULL, n_clusters = NULL, check_arms = TRUE) {
    whole <- function(x) is.numeric(x) && length(x) == 1L && is.finite(x) && x == round(x)
    if (!(whole(n_covariates) && n_covariates >= 0)) {
        .dmar_refuse_argument(c(
            "'n_covariates' must be a single whole number of 0 or more; ",
            .dmar_given(n_covariates), ". It counts the covariates at the cluster ",
            "level beside the treatment indicator, and the interval for the ",
            "difference between the condition means has n_clusters - 2 - ",
            "n_covariates degrees of freedom. Give 0 (the default) for the ",
            "treatment indicator alone."),
            arg = "n_covariates", value = n_covariates)
    }
    for (nm in c("true_R2_between", "true_R2_within")) {
        v <- get(nm)
        if (!(is.numeric(v) && length(v) == 1L && is.finite(v) && v >= 0 && v < 1)) {
            level <- if (nm == "true_R2_between") "between" else "within"
            .dmar_refuse_argument(c(
                "'", nm, "' must be a single number in [0, 1); ", .dmar_given(v),
                ". It is the proportion of the ", level, "-level variance of the ",
                "outcome that the covariates explain in the population. A value ",
                "below 0 or above 1 is not a proportion, and 1 would leave no ",
                level, "-level variance. Give, for example, ", nm, " = 0.3, or 0 ",
                "(the default) for covariates that explain nothing at that level."),
                arg = nm, value = v)
        }
    }
    if (n_covariates == 0 && true_R2_between != 0) {
        .dmar_refuse_argument(c(
            "'true_R2_between' is ", .crd_num(true_R2_between), " with n_covariates ",
            "= 0. It is the share of the between-cluster variance of the outcome ",
            "that the cluster-level covariates explain, and n_covariates = 0 says ",
            "there are none, so no share of it is explained. Give n_covariates, the ",
            "number of cluster-level covariates, or set true_R2_between = 0."),
            arg = "true_R2_between", reason = "R2_between_without_covariates",
            value = true_R2_between)
    }
    if (!is.null(pr_treat) && !(is.numeric(pr_treat) && length(pr_treat) == 1L &&
                                is.finite(pr_treat) && pr_treat > 0 && pr_treat < 1)) {
        .dmar_refuse_argument(c(
            "'pr_treat' must be a single number strictly between 0 and 1; ",
            .dmar_given(pr_treat), ". It is the proportion of the clusters ",
            "assigned to the treatment condition. At 0 or 1 every cluster is in ",
            "one condition and there is no difference between the conditions to ",
            "estimate. Give, for example, pr_treat = 0.5."),
            arg = "pr_treat", value = pr_treat)
    }
    if (!is.null(cluster_size) && !(whole(cluster_size) && cluster_size >= 2)) {
        .dmar_refuse_argument(c(
            "'cluster_size' must be a single whole number of at least 2; ",
            .dmar_given(cluster_size), ". It is the number of individuals in each ",
            "cluster (with 'diff_size', the size the differences are added to or ",
            "multiplied by). The searches over the cluster size start at 2, and ",
            "with one individual in a cluster the variance between clusters cannot ",
            "be told from the variance within them. Give the number of individuals ",
            "per cluster, such as cluster_size = 20."),
            arg = "cluster_size", value = cluster_size)
    }
    if (!is.null(n_clusters)) {
        fewest <- 3 + n_covariates
        if (!(whole(n_clusters) && n_clusters >= fewest)) {
            .dmar_refuse_argument(c(
                "'n_clusters' must be a single whole number of at least ", fewest,
                " (3 + n_covariates); ", .dmar_given(n_clusters), ". The interval ",
                "for the unstandardized difference between the condition means, ",
                "which the standardized planners also start from, has n_clusters - 2 ",
                "- n_covariates degrees of freedom, and ", fewest, " clusters give it ",
                "one. Give the number of clusters of the design."),
                arg = "n_clusters", value = n_clusters)
        }
        if (check_arms && !is.null(pr_treat) && !.crd_both_arms(n_clusters, pr_treat)) {
            n_treat <- round(n_clusters * pr_treat)
            need <- .crd_fewest_both_arms(pr_treat, n_clusters)
            .dmar_refuse_argument(c(
                "With n_clusters = ", n_clusters, " and pr_treat = ", .crd_num(pr_treat),
                ", round(", n_clusters, " * ", .crd_num(pr_treat), ") = ", n_treat,
                " of the clusters are treated, which leaves the ",
                if (n_treat == 0) "treatment" else "control", " condition without a ",
                "cluster, so the difference between the condition means cannot be ",
                "estimated and its interval has no finite width. ", need, " clusters ",
                "are the fewest that put a cluster in each condition at this ",
                "pr_treat; give at least that many, or a pr_treat nearer 0.5."),
                arg = c("n_clusters", "pr_treat"), reason = "arm_without_clusters",
                n_clusters = n_clusters, n_treated = n_treat, fewest = need)
        }
    }
    invisible(NULL)
}

# .crd_check_costs(cluster_cost, individual_cost, solves)
#
# The entry check of the costs. Each cost that is given must be a single
# finite number of 0 or more. `solves` says what the costs decide:
#   "none"          the planners that take a width and a fixed size, where
#                   the costs only price the design reported (its budget
#                   row); either cost may be NULL, and costs of 0 give a
#                   budget of 0.
#   "n_clusters"    the planners that buy clusters of a fixed size with a
#                   budget: the costs are required and must not both be 0
#                   (reason "costs_all_zero"), or every number of clusters
#                   is within the budget and the search has no end.
#   "cluster_size"  the planners that buy the cluster size with a budget or
#                   choose it by cost (the n_individuals budget planners and
#                   the both planners): in addition 'individual_cost' must
#                   be positive (reason "individual_cost_zero"), or clusters
#                   of every size cost the same and the cluster size has no
#                   bound.
# Through the first builds of 1.1.0 costs of 0 searched without end in six
# planners, ss_aipe_crd_both_fixed_width() returned a table of NA, and
# ss_aipe_crd_es_both_fixed_width() stopped with "missing value where
# TRUE/FALSE needed"; a negative cost or a vector of costs was not refused.
.crd_check_costs <- function(cluster_cost, individual_cost, solves = c("none", "n_clusters", "cluster_size")) {
    solves <- match.arg(solves)
    role <- c(cluster_cost = "the cost of adding one cluster, whatever the number of individuals in it",
              individual_cost = "the cost of adding one individual to a cluster")
    vals <- list(cluster_cost = cluster_cost, individual_cost = individual_cost)
    for (nm in names(vals)) {
        x <- vals[[nm]]
        if (is.null(x) && solves == "none") next
        if (!(is.numeric(x) && length(x) == 1L && is.finite(x) && x >= 0)) {
            .dmar_refuse_argument(c(
                "'", nm, "' must be a single finite number of 0 or more; ",
                .dmar_given(x), ". It is ", role[[nm]], ", in the units of the ",
                "budget. Give one cost of 0 or more, such as ", nm, " = ",
                if (nm == "cluster_cost") "30" else "1", "."),
                arg = nm, value = x)
        }
    }
    if (solves == "none") return(invisible(NULL))
    if (cluster_cost == 0 && individual_cost == 0) {
        .dmar_refuse_argument(c(
            "'cluster_cost' and 'individual_cost' are both 0, so every design ",
            "costs 0 and the costs cannot choose among designs: a budget buys ",
            "any number of clusters of any size, and every design that reaches ",
            "a width costs the same. Give the costs in the units of the budget, ",
            "for example cluster_cost = 30 and individual_cost = 1."),
            arg = c("cluster_cost", "individual_cost"), reason = "costs_all_zero")
    }
    if (solves == "cluster_size" && individual_cost == 0) {
        .dmar_refuse_argument(c(
            "'individual_cost' is 0, so an individual added to a cluster costs ",
            "nothing and the cluster size the costs choose has no bound: a budget ",
            "buys clusters of any size, and at a fixed number of clusters every ",
            "cluster size costs the same. Give the cost of one individual in the ",
            "units of 'cluster_cost', for example individual_cost = 1."),
            arg = "individual_cost", reason = "individual_cost_zero")
    }
    invisible(NULL)
}

# .crd_check_budget(budget)
#
# The entry check of `budget` in the planners that buy a design with it: a
# single finite number of 0 or more. A budget too small for the smallest
# design is refused after the search (.crd_refuse_budget()).
.crd_check_budget <- function(budget) {
    if (is.numeric(budget) && length(budget) == 1L && is.finite(budget) && budget >= 0) {
        return(invisible(budget))
    }
    .dmar_refuse_argument(c(
        "'budget' must be a single finite number of 0 or more; ",
        .dmar_given(budget), ". It is the amount available for the study, in ",
        "the units of 'cluster_cost' and 'individual_cost'. Give one amount of 0 ",
        "or more, such as budget = 10000."),
        arg = "budget", value = budget)
}

# .crd_refuse_budget(fn, budget, n_clusters, cluster_size, cluster_cost,
#                    individual_cost, diff_size, floor_text)
#
# The refusal of a budget that does not buy the smallest design a planner
# can return (`n_clusters` clusters of `cluster_size`): class
# dmar_argument_invalid, arg "budget", reason "budget_infeasible", with
# the budget and the cost of that design. `floor_text` says why the
# number of clusters is what it is, where that needs saying (a design the
# planner returns with its width has a cluster in each condition and a
# degree of freedom; .crd_fixed_size_fewest()), so the budget the message
# names buys a design that is returned. Through the
# first builds of 1.1.0 the planners that fix a size returned a table
# whose size and budget were NA, and the planners that search over both
# sizes stopped with "invalid 'type' (list) of argument".
.crd_refuse_budget <- function(fn, budget, n_clusters, cluster_size, cluster_cost, individual_cost, diff_size = NULL, floor_text = NULL) {
    cost <- .cost_crd(n_clusters, cluster_size, cluster_cost = cluster_cost, individual_cost = individual_cost, diff_size = diff_size)
    notes <- c(if (!is.null(diff_size)) "before 'diff_size' is applied",
               if (length(floor_text)) paste0(floor_text, collapse = ""))
    .dmar_refuse_argument(c(
        fn, "(): a budget of ", .crd_num(budget), " does not buy the smallest ",
        "design this planner can return, ", n_clusters, " clusters of ",
        cluster_size, " individuals",
        if (length(notes)) c(" (", paste(notes, collapse = "; "), ")"),
        ", which costs ", .crd_num(cost), " at cluster_cost = ",
        .crd_num(cluster_cost), " and individual_cost = ", .crd_num(individual_cost),
        ". Give a budget of at least ", .crd_num(cost), ", or lower costs."),
        arg = "budget", reason = "budget_infeasible", budget = budget,
        min_cost = cost, n_clusters = n_clusters, cluster_size = cluster_size)
}

# The floor of the planners that search over both sizes, as their
# refusals state it.
.crd_floor_text <- function(pr_treat, n_covariates) {
    c("at least two clusters in each condition at pr_treat = ", .crd_num(pr_treat),
      if (n_covariates == 1) {
          ", plus one for the covariate"
      } else if (n_covariates > 1) {
          c(", plus one for each of the ", n_covariates, " covariates")
      })
}

# .crd_check_bought_design(fn, budget, n_clusters, cluster_size, pr_treat,
#                          n_covariates, cluster_cost, individual_cost,
#                          diff_size)
#
# The check, in the planners that buy clusters of a fixed size with a
# budget and report the width when 'pr_treat' is given, that the design
# the budget buys has a finite width: a cluster in each condition (reason
# "arm_without_clusters") and at least 3 + n_covariates clusters (reason
# "too_few_clusters"). Class dmar_argument_invalid, with the budget, the
# design, and the budget that buys the fewest clusters that would do.
# Through the first builds of 1.1.0 such a design was returned with an
# expected width of Inf or NaN (budget = 100, cluster_size = 20, pr_treat
# = .1 bought 5 clusters, none treated).
.crd_check_bought_design <- function(fn, budget, n_clusters, cluster_size, pr_treat, n_covariates, cluster_cost, individual_cost, diff_size = NULL) {
    both <- .crd_both_arms(n_clusters, pr_treat)
    enough <- n_clusters >= 3 + n_covariates
    if (both && enough) return(invisible(NULL))
    need <- .crd_fewest_both_arms(pr_treat, max(n_clusters, 3 + n_covariates))
    cost <- .cost_crd(n_clusters, cluster_size, cluster_cost = cluster_cost, individual_cost = individual_cost, diff_size = diff_size)
    need_cost <- .cost_crd(need, cluster_size, cluster_cost = cluster_cost, individual_cost = individual_cost, diff_size = diff_size)
    bought <- c(fn, "(): a budget of ", .crd_num(budget), " buys ", n_clusters,
                " clusters of ", cluster_size, " individuals (cost ", .crd_num(cost), ")")
    if (!both) {
        n_treat <- round(n_clusters * pr_treat)
        .dmar_refuse_argument(c(
            bought, ", and with pr_treat = ", .crd_num(pr_treat), " round(", n_clusters,
            " * ", .crd_num(pr_treat), ") = ", n_treat, " of them are treated, so the ",
            "design leaves the ", if (n_treat == 0) "treatment" else "control",
            " condition without clusters and the interval for the difference ",
            "between the condition means has no finite width. ", need, " clusters, ",
            "the fewest that put a cluster in each condition here, cost ",
            .crd_num(need_cost), "; give a budget of at least that, or a pr_treat ",
            "nearer 0.5."),
            arg = c("budget", "pr_treat"), reason = "arm_without_clusters",
            budget = budget, n_clusters = n_clusters, cluster_size = cluster_size,
            n_treated = n_treat, fewest = need, fewest_cost = need_cost)
    }
    .dmar_refuse_argument(c(
        bought, ", and with n_covariates = ", n_covariates, " the interval for the ",
        "difference between the condition means has ", n_clusters, " - 2 - ",
        n_covariates, " = ", n_clusters - 2 - n_covariates, " degrees of freedom, ",
        "so it has no finite width. ", need, " clusters cost ", .crd_num(need_cost),
        "; give a budget of at least that, or fewer covariates."),
        arg = c("budget", "n_covariates"), reason = "too_few_clusters",
        budget = budget, n_clusters = n_clusters, cluster_size = cluster_size,
        fewest = need, fewest_cost = need_cost)
}

# .crd_refuse_width_unreachable(fn, width, n_clusters, pr_treat,
#                               true_var_between, true_var_within, ...)
#
# The refusal of a target width that no cluster size up to 100,000, the
# largest .find_n_indiv_crd_diff() considers, reaches with the given number
# of clusters: class dmar_argument_invalid, arg "width", reason
# "width_not_reachable", as in ss_aipe_crd_es_n_individuals_fixed_width().
# The message gives the width at 100,000 and the value the width falls
# toward as the clusters grow, with every cluster mean known exactly
# (`width_limit`, set by the variance between clusters), and says whether a
# larger cluster size would reach the target. Through the first builds of
# 1.1.0 the refusal said only that the target was "impossible to achieve"
# and had no reason.
.crd_refuse_width_unreachable <- function(fn, width, n_clusters, pr_treat, true_var_between, true_var_within, true_R2_between = 0, true_R2_within = 0, n_covariates = 0, assurance = NULL, conf_level = 0.95, diff_size = NULL) {
    at <- function(cluster_size, diff_size) .find_width_crd_diff(n_clusters = n_clusters, cluster_size = cluster_size, pr_treat = pr_treat, true_var_between = true_var_between, true_var_within = true_var_within, true_R2_between = true_R2_between, true_R2_within = true_R2_within, n_covariates = n_covariates, assurance = assurance, conf_level = conf_level, diff_size = diff_size)
    at_max <- at(100000, diff_size)
    limit <- at(Inf, NULL)
    kind <- if (is.null(assurance)) "expected" else "assured"
    .dmar_refuse_argument(c(
        fn, "(): with ", n_clusters, " clusters (", round(n_clusters * pr_treat),
        " of them treated), the ", kind, " width at a cluster size of 100,000, ",
        "the largest the search considers, is ", format(signif(at_max, 4)),
        ", above the target 'width' of ", format(width), ". As the clusters grow ",
        "the width falls toward ", format(signif(limit, 4)), ", its value with ",
        "every cluster mean known exactly, which the variance between clusters ",
        "sets; ",
        if (width > limit) {
            "a cluster size above 100,000 would reach the target"
        } else {
            "it does not reach it, so no cluster size reaches the target with this number of clusters"
        },
        ". Give more clusters, or plan the number of clusters with ",
        "ss_aipe_crd_n_clusters_fixed_width()."),
        arg = "width", reason = "width_not_reachable", width_limit = limit,
        width_at_max = at_max, max_cluster_size = 100000, target = width,
        n_clusters = n_clusters)
}

# The refusal of a search over both sizes that found no candidate design
# (.find_min_width_crd_diff(), .find_min_cost_crd_diff()): class
# dmar_computation_failed, reason "no_candidate". The public planners check
# what makes a candidate before searching (the budget buys the floor's
# number of clusters of 2; a target beyond .CRD_MAX_N_CLUSTERS is refused),
# so this is a safeguard. `fn` names the public planner.
.crd_refuse_no_candidate <- function(fn, lo, pr_treat, n_covariates, searched, budget = NULL, width = NULL) {
    .dmar_refuse(c(
        .crd_fn_text(fn), ": the search over the number of clusters found no ",
        "candidate design from ", .crd_num(searched[1]), " to ", .crd_num(searched[2]),
        " clusters. A candidate has at least ", .crd_num(lo), " clusters (",
        .crd_floor_text(pr_treat, n_covariates), ")",
        if (!is.null(budget)) {
            c(" and is bought by the budget of ", .crd_num(budget), " with clusters of at least 2")
        } else {
            c(" and reaches the target width of ", format(width), " with clusters of at most 100,000")
        },
        ". Give ", if (is.null(budget)) "a wider target" else "a larger budget, lower costs,",
        " or a pr_treat nearer 0.5, which lowers the number of clusters a candidate needs."),
        "computation_failed", reason = "no_candidate", floor = lo, searched = searched)
}

# .crd_resolve_variances(true_var_between, true_var_within, true_total_var,
#                        true_icc, diff_size = NULL, rel_tol)
#
# The variance model of the six unstandardized planners, resolved once per
# call from the four inputs, any two of which determine the other two:
# true_total_var = true_var_between + true_var_within and true_icc =
# true_var_between / true_total_var. Returns a list with true_var_between
# and true_var_within, the two the width is computed from, and
# true_total_var and true_icc, with the attribute `identified_by` (the
# pair used). Every refusal has class dmar_argument_invalid, names the
# arguments at issue as `arg`, and carries `reason` and `supplied` (the
# four inputs as given).
#
# Each input that is given must be a single finite number (reason
# "variance_input_invalid"): a variance 0 or more, the total variance
# positive, and the intraclass correlation in [0, 1] (reason
# "variance_out_of_range"; so is a pair whose total, between plus within,
# is 0 or not finite, as two variances of 0 are). The pair used is the first given of
# between + within, total + ICC, total + between, total + within, between +
# ICC, and within + ICC; fewer than two inputs (or between 0 with ICC 0, or
# within 0 with ICC 1, which leave the scale open) are refused with reason
# "underidentified_variances". The two variances are formed as the first
# builds of 1.1.0 formed them (between = icc * total and within = (1 - icc)
# * total from total and ICC; within = total - between and between = total -
# within), and a given variance is used as given, so every call those
# builds planned for with consistent inputs gives the same plan. A pair
# that implies a negative variance, or an
# input given beside the pair that disagrees with the value the pair
# implies by more than `rel_tol` (64 machine epsilons) relative to the
# total variance, or absolutely for the intraclass correlation, is refused
# with reason "inconsistent_variances"; the message gives the implied
# value, so a rounded input (an ICC of .333 beside variances of 1 and 2) is
# seen for what it is. A difference of a total and a variance that is
# negative by no more than that tolerance is taken as 0. With unequal
# cluster sizes (`diff_size`) the within-cluster variance must be positive
# (reason "within_variance_zero"), since the weight of each cluster mean
# divides by it.
#
# Through the first builds of 1.1.0 the inputs were resolved inside the
# search at every number of clusters tried: the documented pairs (between,
# ICC) and (within, ICC) stopped with "missing value where TRUE/FALSE
# needed", as did a negative variance after ten warnings, and an input
# given beside a pair was used without a check (true_total_var = 1,
# true_icc = .2, true_var_between = .5 planned with a between-cluster
# variance of .5 and a within-cluster variance of .8, 36 clusters where the
# pair gives 18).
.crd_resolve_variances <- function(true_var_between = NULL, true_var_within = NULL, true_total_var = NULL, true_icc = NULL, diff_size = NULL, rel_tol = 64 * .Machine$double.eps) {
    supplied <- list(true_var_between = true_var_between, true_var_within = true_var_within,
                     true_total_var = true_total_var, true_icc = true_icc)
    given <- names(supplied)[!vapply(supplied, is.null, logical(1))]
    refuse <- function(message, arg, reason, ...) {
        .dmar_refuse_argument(message, arg = arg, reason = reason, supplied = supplied, ...)
    }
    role <- c(true_var_between = "the population variance between clusters, of the cluster means about their condition means",
              true_var_within = "the population variance within clusters",
              true_total_var = "the population total variance, between clusters plus within them",
              true_icc = "the population intraclass correlation, the share of the total variance that lies between clusters")
    for (nm in given) {
        x <- supplied[[nm]]
        if (!(is.numeric(x) && length(x) == 1L && is.finite(x))) {
            refuse(c("'", nm, "' must be NULL or a single finite number; ", .dmar_given(x),
                     ". It is ", role[[nm]], " (a planning value). Give one number, ",
                     "or leave it NULL and give two of the other three."),
                   nm, "variance_input_invalid", value = x)
        }
        bad <- switch(nm,
                      true_icc = if (x < 0 || x > 1) "must be in [0, 1]",
                      true_total_var = if (x <= 0) "must be positive",
                      if (x < 0) "must not be negative")
        if (!is.null(bad)) {
            refuse(c("'", nm, "' ", bad, "; ", .dmar_given(x), ". It is ", role[[nm]],
                     " (a planning value). Give ",
                     switch(nm, true_icc = "a proportion, such as true_icc = 0.2",
                            true_total_var = "a positive variance, such as true_total_var = 1",
                            paste0("a variance of 0 or more, such as ", nm, " = ",
                                   if (nm == "true_var_between") "0.2" else "0.8")),
                     "."),
                   nm, "variance_out_of_range", value = x)
        }
    }
    # Values in the messages: 15 significant digits, so that a rounded
    # input is seen beside the value implied, and in scientific notation
    # where that is shorter.
    num <- function(x) format(x, digits = 15)
    shown <- function(nm) paste0(nm, " = ", num(supplied[[nm]]))
    has <- function(nm) nm %in% given
    between <- true_var_between
    within <- true_var_within
    total <- true_total_var
    icc <- true_icc
    # A difference of the total and a variance, taken as 0 when it is
    # negative by no more than the tolerance.
    rest <- function(x, of, from, implied) {
        if (x >= 0) return(x)
        if (-x <= rel_tol * of) return(0)
        refuse(c("The variances given do not agree: ", shown(from[1]), " and ",
                 shown(from[2]), " imply ", implied, " = ", num(x), ", and a ",
                 "variance cannot be negative. The total variance is the variance ",
                 "between clusters plus the variance within them; give two of the ",
                 "four and leave the others NULL."),
               from, "inconsistent_variances", implied = stats::setNames(x, implied))
    }
    # The two variances the pair determines, formed as the first builds of
    # 1.1.0 formed them.
    if (has("true_var_between") && has("true_var_within")) {
        pair <- c("true_var_between", "true_var_within")
    } else if (has("true_total_var") && has("true_icc")) {
        pair <- c("true_total_var", "true_icc")
        between <- icc * total
        within <- (1 - icc) * total
    } else if (has("true_total_var") && has("true_var_between")) {
        pair <- c("true_total_var", "true_var_between")
        within <- rest(total - between, total, pair, "true_var_within")
    } else if (has("true_total_var") && has("true_var_within")) {
        pair <- c("true_total_var", "true_var_within")
        between <- rest(total - within, total, pair, "true_var_between")
    } else if (has("true_var_between") && has("true_icc")) {
        pair <- c("true_var_between", "true_icc")
        if (icc == 0) {
            if (between == 0) {
                refuse(c("true_var_between = 0 with true_icc = 0 does not determine the ",
                         "variances: with no variance between clusters the intraclass ",
                         "correlation is 0 whatever the variance within them. Give ",
                         "'true_var_within' or 'true_total_var' as well."),
                       pair, "underidentified_variances")
            }
            refuse(c("The variances given do not agree: ", shown("true_var_between"),
                     " is positive, and true_icc = 0 says that no variance lies between ",
                     "clusters. Give two of the four and leave the others NULL."),
                   pair, "inconsistent_variances")
        }
        total <- between / icc
        within <- (1 - icc) * total
    } else if (has("true_var_within") && has("true_icc")) {
        pair <- c("true_var_within", "true_icc")
        if (icc == 1) {
            if (within == 0) {
                refuse(c("true_var_within = 0 with true_icc = 1 does not determine the ",
                         "variances: with no variance within clusters the intraclass ",
                         "correlation is 1 whatever the variance between them. Give ",
                         "'true_var_between' or 'true_total_var' as well."),
                       pair, "underidentified_variances")
            }
            refuse(c("The variances given do not agree: ", shown("true_var_within"),
                     " is positive, and true_icc = 1 says that all of the variance lies ",
                     "between clusters. Give two of the four and leave the others NULL."),
                   pair, "inconsistent_variances")
        }
        total <- within / (1 - icc)
        between <- icc * total
    } else {
        refuse(c("The planners take two of 'true_var_between', 'true_var_within', ",
                 "'true_total_var', and 'true_icc', which determine the other two; ",
                 if (length(given)) paste0("only '", given, "' was given") else "none was given",
                 ". Give, for example, true_var_between = 0.2 and true_var_within = 0.8, ",
                 "or true_total_var = 1 and true_icc = 0.2."),
               c("true_var_between", "true_var_within", "true_total_var", "true_icc"),
               "underidentified_variances")
    }
    implied_total <- between + within
    if (!is.finite(implied_total) || implied_total <= 0) {
        refuse(c("The variances given imply a total variance of ", num(implied_total),
                 " (", num(between), " between clusters plus ", num(within),
                 " within them), and the total variance must be a positive finite ",
                 "number: with no variance there is no difference for an interval to ",
                 "be uncertain about. Give variances whose sum is positive and finite, ",
                 "such as true_var_between = 0.2 and true_var_within = 0.8."),
               pair, "variance_out_of_range")
    }
    implied <- list(true_var_between = between, true_var_within = within,
                    true_total_var = implied_total, true_icc = between / implied_total)
    # An input given beside the pair agrees with the value the pair implies.
    for (nm in setdiff(given, pair)) {
        scale <- if (nm == "true_icc") 1 else implied_total
        if (abs(supplied[[nm]] - implied[[nm]]) > rel_tol * scale) {
            refuse(c("The variances given do not agree: ", shown(pair[1]), " and ",
                     shown(pair[2]), " imply ", nm, " = ", num(implied[[nm]]),
                     ", but ", num(supplied[[nm]]), " was given. An input given ",
                     "beside two that determine the model must agree with them to within ",
                     "about 1e-14 of its scale, so a rounded value is refused; give two of ",
                     "the four and leave the others NULL."),
                   c(pair, nm), "inconsistent_variances",
                   implied = stats::setNames(implied[[nm]], nm))
        }
    }
    # A variance that was given is used as given (it agrees with the pair to
    # within the tolerance), as the first builds of 1.1.0 used it.
    if (has("true_var_between")) between <- true_var_between
    if (has("true_var_within")) within <- true_var_within
    if (!is.null(diff_size) && within == 0) {
        refuse(c("With unequal cluster sizes ('diff_size') the variance within clusters ",
                 "must be positive, and the variances given make it 0. The planners ",
                 "weight each cluster mean by its precision, which divides by the ",
                 "variance within clusters. Give a positive variance within clusters, ",
                 "or plan with equal cluster sizes (diff_size = NULL)."),
               c(given, "diff_size"), "within_variance_zero")
    }
    implied_total <- between + within
    implied <- list(true_var_between = between, true_var_within = within,
                    true_total_var = implied_total, true_icc = between / implied_total)
    structure(implied, identified_by = pair)
}
