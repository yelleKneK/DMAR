# Shared helpers for the missing data family. Sections (A) to (D) serve
# the missing data descriptives and the fraction of missing information:
# fraction_missing_information(), missing_patterns(), and
# missing_coverage(). The public files validate their arguments and
# assemble the returned table; the numerics live here. Section (E)
# holds the small-sample degrees of freedom of df_missing() and
# df_barnard_rubin(), and Section (F) the t reference distribution that
# the Wald intervals of reliability_omega(), reliability_alpha(),
# fa_confirmatory(), mlmr(), and mlmr_mv() take through their wald_df argument.


# ---------------------------------------------------------------------------
# (A) The item frame shared by missing_patterns() and missing_coverage()
# ---------------------------------------------------------------------------

# .missing_item_frame(data, items)
#
# Validates `data` and `items` and returns the logical observed matrix
# (TRUE = observed, one column per item, in the order of `items`). Only
# is.na() matters, so any column type is allowed.
.missing_item_frame <- function(data, items = NULL) {
  if (!is.data.frame(data)) {
    .dmar_refuse_argument("'data' must be a data.frame.", arg = "data")
  }
  if (is.null(items)) items <- names(data)
  if (!is.character(items) || length(items) == 0L || anyNA(items)) {
    .dmar_refuse_argument(c(
      "'items' must be NULL or a character vector naming columns of ",
      "'data'."), arg = c("items", "data"))
  }
  if (anyDuplicated(items)) {
    .dmar_refuse_argument(c(
      "'items' names a column more than once: ",
      paste(unique(items[duplicated(items)]), collapse = ", "), "."),
      arg = "items")
  }
  absent <- setdiff(items, names(data))
  if (length(absent)) {
    .dmar_refuse_argument(c(
      "'items' names columns that are not in 'data': ",
      paste(absent, collapse = ", "), "."), arg = c("items", "data"))
  }
  if (nrow(data) == 0L) {
    .dmar_refuse_argument("'data' has no rows.", arg = "data")
  }
  # One logical per row and item, TRUE where a value was observed. A
  # vapply() over the columns keeps the matrix n by p even with one row.
  matrix(vapply(data[items], function(x) !is.na(x), logical(nrow(data))),
         nrow = nrow(data), dimnames = list(NULL, items))
}


# ---------------------------------------------------------------------------
# (B) Checks on the fit
# ---------------------------------------------------------------------------

# .fmi_extract_fit(fit)
#
# Accepts a lavaan fit, or an mlmr / mlmr_mv object whose lavaan fit is
# used, and refuses everything else with a sentence that says what to pass.
.fmi_extract_fit <- function(fit) {
  if (inherits(fit, c("mlmr", "mlmr_mv"))) fit <- fit$lavaan_fit
  if (inherits(fit, "dmar_fa_confirmatory")) {
    .dmar_refuse_argument(c(
      "'fit' is a fa_confirmatory() table, not a fitted model; refit with ",
      "fa_confirmatory(..., missing = \"fiml\", output = \"fit\") and pass that ",
      "object."), arg = "fit")
  }
  if (inherits(fit, "lavaanList") || inherits(fit, "lavaan.mi")) {
    .dmar_refuse_argument(c(
      "'fit' is a lavaanList (or pooled multiple imputation) object. ",
      "The fraction of missing information computed here is defined ",
      "for a single full information maximum likelihood fit; with ",
      "multiple imputation it comes from the pooling itself (the ",
      "between imputation variance over the total variance)."), arg = "fit")
  }
  if (!inherits(fit, "lavaan")) {
    .dmar_refuse_argument(c(
      "'fit' must be a fitted lavaan model (for example from ",
      "fa_confirmatory(..., missing = \"fiml\", output = \"fit\") or ",
      "lavaan::cfa(..., missing = \"ml\")), or an mlmr() or mlmr_mv() ",
      "fit."), arg = "fit")
  }
  fit
}

# .fmi_mlmr_names(x)
#
# mlmr() and mlmr_mv() fit their lavaan models to renamed variables (y or
# y_1, ..., x_1, ..., aux_1, ...). This maps those names back to the ones
# the user wrote, so the terms of the FMI table read as the formula does.
# NULL for a fit passed directly.
.fmi_mlmr_names <- function(x) {
  if (inherits(x, "mlmr")) {
    internal <- c(x$y_name, x$X_names)
    display <- c(paste(deparse(x$formula[[2L]]), collapse = ""),
                 x$X_display)
  } else if (inherits(x, "mlmr_mv")) {
    internal <- c(x$Y_names, x$X_names)
    display <- c(x$Y_display, x$X_display)
  } else {
    return(NULL)
  }
  if (length(x$auxiliary)) {
    internal <- c(internal, paste0("aux_", seq_along(x$auxiliary)))
    display <- c(display, x$auxiliary)
  }
  if (length(internal) != length(display)) return(NULL)
  stats::setNames(as.character(display), internal)
}

# .fmi_with_names(expr, names)
#
# Evaluates expr and gives the warnings and errors it signals again in the
# names the user wrote for the variables of an mlmr() or mlmr_mv() fit.
# names maps the names the fit gave them (y or y_1, ..., x_1, ..., aux_1,
# ...) to the user's, as .fmi_mlmr_names() returns it. The messages of
# the computation name a variable (one whose mean is far from zero) or a
# parameter (an estimate a Newton step moves), and lavaan's warnings name
# variables too, all in the fit's names. Words are replaced whole and all
# at once (.relia_lavaan_names_back() in R/reliability_internals.R), so
# x_1 is never read inside x_10 and two names are never exchanged twice.
# With names NULL, expr is evaluated as it is.
.fmi_with_names <- function(expr, names) {
  if (!length(names)) return(expr)
  rename <- function(text) {
    .relia_lavaan_names_back(text, names(names), unname(names))
  }
  withCallingHandlers(
    tryCatch(expr, error = function(e) {
      # The condition keeps its class (dmar_fmi_saturated_singular, say).
      e$message <- rename(conditionMessage(e))
      e$call <- NULL
      stop(e)
    }),
    warning = function(w) {
      .dmar_warn(rename(conditionMessage(w)), "result_qualified")
      invokeRestart("muffleWarning")
    })
}

# .fmi_check_fit(fit)
#
# The refusals. FMI as computed here is a ratio of two normal theory
# information based variances at the estimates of a converged, single
# level, unweighted full information maximum likelihood fit to continuous
# variables. Each case outside that is stopped with the reason and the
# remedy. Returns the fit's options invisibly.
.fmi_check_fit <- function(fit) {
  opt <- lavaan::lavInspect(fit, "options")
  est_name <- opt$estimator.orig %||% opt$estimator
  if (isTRUE(lavaan::lavInspect(fit, "categorical"))) {
    .dmar_refuse_argument(c(
      "The fit treats some indicators as ordered (estimator ", est_name,
      "). The fraction of missing information is defined through the ",
      "likelihood information of a full information maximum likelihood ",
      "fit, which a least squares fit to polychoric correlations does ",
      "not have. Fit the items as continuous with missing = \"fiml\" ",
      "(lavaan: missing = \"ml\") and estimator = \"MLR\", or use ",
      "multiple imputation, which DMAR does not yet provide."), arg = NULL)
  }
  if (!identical(opt$estimator, "ML")) {
    .dmar_refuse_argument(c(
      "The fit used estimator ", est_name, ". The fraction of missing ",
      "information is defined through the likelihood information, so ",
      "it is available for the maximum likelihood family (ML, and MLR ",
      "or MLF, whose estimates are the ML estimates) under full ",
      "information maximum likelihood."), arg = NULL)
  }
  if (isTRUE(opt$conditional.x)) {
    .dmar_refuse_argument(c(
      "The fit conditions on the exogenous covariates ",
      "(conditional.x = TRUE), and its unrestricted moments are then ",
      "regression moments that the complete data comparison cannot ",
      "use. Refit with conditional.x = FALSE (lavaan's default for ",
      "continuous variables)."), arg = NULL)
  }
  missing <- as.character(opt$missing)
  if (identical(missing, "listwise")) {
    .dmar_refuse_argument(c(
      "The fit used listwise deletion (lavaan's default), which drops ",
      "every incomplete case, so it carries no missing information to ",
      "measure. Refit with missing = \"fiml\" (lavaan: missing = ",
      "\"ml\"), which uses every observed value; with complete data ",
      "the fraction is then zero for every parameter."), arg = NULL)
  }
  if (!missing %in% c("ml", "ml.x")) {
    .dmar_refuse_argument(c(
      "The fit used missing = \"", missing, "\". The fraction of ",
      "missing information is available for full information maximum ",
      "likelihood fits (missing = \"ml\" or \"ml.x\"); the pairwise, ",
      "available cases, and two-stage options are not supported."), arg = NULL)
  }
  if (lavaan::lavInspect(fit, "nlevels") > 1L) {
    .dmar_refuse_argument(c(
      "The fit is a multilevel model. The fraction of missing ",
      "information for multilevel fits is not available."), arg = NULL)
  }
  if (.fmi_has_sampling_weights(fit, opt)) {
    .dmar_refuse_argument(c(
      "The fit used sampling weights. Its estimates maximize a ",
      "weighted pseudo-likelihood whose information is not the one ",
      "the fraction of missing information is defined through; refit ",
      "without weights to compute it."), arg = NULL)
  }
  if (!isTRUE(lavaan::lavInspect(fit, "converged"))) {
    .dmar_refuse(c(
      "The fit did not converge, so there are no estimates at which ",
      "to evaluate the information."), "not_converged")
  }
  pt <- lavaan::parTable(fit)
  if (!is.null(pt$efa) && any(nzchar(pt$efa))) {
    .dmar_refuse_argument(c(
      "The fit contains a rotated exploratory (efa) block. Its ",
      "loadings are rotated after estimation, so their standard errors ",
      "are not the inverse of the likelihood information at the ",
      "reported estimates, and a refit at those estimates rotates them ",
      "again. Fit the retained structure as a confirmatory model (for ",
      "example with fa_confirmatory()) and compute the fraction there."),
      arg = NULL)
  }
  invisible(opt)
}

# Whether the fit used sampling weights. lavaan records the fact in a hidden
# option; its documented inspector, which errors when there are none, is the
# backstop should the option ever be renamed.
.fmi_has_sampling_weights <- function(fit, opt) {
  if (isTRUE(opt[[".sampling.weights"]])) return(TRUE)
  tryCatch({
    lavaan::lavInspect(fit, "sampling.weights")
    TRUE
  }, error = function(e) FALSE)
}


# ---------------------------------------------------------------------------
# (C) The two refits and the ratio
# ---------------------------------------------------------------------------

# lavaan options for each information keyword (Chen & Savalei, 2021,
# Appendix A; lavaan 0.6-7 and later). All three are observed information;
# they differ in how the h1 (unrestricted) part is evaluated.
.fmi_information_options <- function(information) {
  switch(information,
    unstructured = list(observed.information = "h1",
                        h1.information = "unstructured"),
    structured   = list(observed.information = "h1",
                        h1.information = "structured"),
    hessian      = list(observed.information = "hessian",
                        h1.information = "structured"))
}

.fmi_information_description <- function(information) {
  switch(information,
    unstructured = paste0("observed information of the unrestricted (h1) ",
                          "model at the saturated FIML moments, carried ",
                          "to the parameters through the model Jacobian ",
                          "(the third estimator of Chen and Savalei, 2021)"),
    structured   = paste0("first term of the observed information at the ",
                          "model implied moments (the second estimator of ",
                          "Chen and Savalei, 2021)"),
    hessian      = paste0("numeric Hessian of the log likelihood, lavaan's ",
                          "default (the first estimator of Chen and ",
                          "Savalei, 2021)"))
}

# .fmi_refit_scale(pt, dat)
#
# The metric of the two refits. Each observed variable is divided by
# c_j = 2^round(log2(SD_j)), the power of two nearest its available case
# standard deviation over all groups, the factor lavaan itself divides a
# variable by when it fits a badly scaled model (lav_rescale_prep() in
# lavaan 0.7-3), and each row of the parameter table by the factor that
# map implies for it (.fmi_rescale_mono()). A variable whose standard
# deviation cannot be computed keeps the factor 1. Division by a power of
# two is exact in floating point, so the refits hold the fit's estimates
# exactly in the new metric, and their standard errors map back exactly;
# the choice of the factors governs only how well conditioned the
# computation is, never what it computes.
#
# The refits stay in the raw metric when every factor is 1 (every SD
# between about .71 and 1.41, where there is nothing to gain) and when the
# parameter table cannot be mapped: composites (<~) and product terms
# (x1:x2), which lavaan does not rescale either, and rows that share one
# free parameter while the map gives them different factors (a label
# shared under ceq.simple = TRUE by parameters of variables in different
# units). For the last, one common factor for every variable is tried
# first, since it gives every loading, every regression among observed
# variables, and every variance a factor of its own kind. Constraints and
# defined parameters do not block the map, because .fmi_rescale_pt()
# rewrites their expressions.
#
# Inputs: pt, the fit's parameter table; dat, a list with one data matrix
# per group, as lavInspect(fit, "data") gives them (named columns).
# Returns list(cj, mono, reason): cj the named factors, or NULL for the
# raw metric; mono one factor per row of pt (1 for the rows that are not
# parameters); and reason, the words for why the refits stay in the raw
# metric when the factors are not all 1, or NULL.
.fmi_refit_scale <- function(pt, dat) {
  n <- length(pt$lhs)
  raw <- function(reason = NULL) {
    list(cj = NULL, mono = rep(1, n), reason = reason)
  }
  ov <- unique(unlist(lapply(dat, colnames)))
  sds <- vapply(ov, function(v) {
    vals <- unlist(lapply(dat, function(x) {
      if (v %in% colnames(x)) x[, v] else NULL
    }))
    vals <- vals[is.finite(vals)]
    if (length(vals) > 1L) stats::sd(vals) else NA_real_
  }, numeric(1L))
  ok <- is.finite(sds) & sds > 0
  cj <- stats::setNames(rep(1, length(ov)), ov)
  cj[ok] <- 2^round(log2(sds[ok]))
  if (all(cj == 1)) return(raw())
  is_par <- !pt$op %in% c(":=", "==", "<", ">")
  if (any(pt$op == "<~")) return(raw("composites"))
  if (any(grepl(":", c(pt$lhs[is_par], pt$rhs[is_par]), fixed = TRUE))) {
    return(raw("product terms"))
  }
  mono <- .fmi_rescale_mono(pt, cj)
  if (is.null(mono) && any(ok)) {
    cj[] <- 2^round(mean(log2(sds[ok])))
    mono <- if (all(cj == 1)) NULL else .fmi_rescale_mono(pt, cj)
  }
  if (is.null(mono)) {
    return(raw(paste0("parameters that share one label across variables ",
                      "in different units")))
  }
  list(cj = cj, mono = mono, reason = NULL)
}

# .fmi_rescale_mono(pt, cj)
#
# The factor by which each row of the parameter table is divided when the
# observed variables are divided by cj, as lavaan's lav_rescale_mono()
# computes it: the scale of the indicator over that of the factor for a
# loading (=~), the scale of the outcome over that of the predictor for a
# regression (~), the product of the two scales for a variance or
# covariance (~~), and the variable's scale for an intercept or mean (~1).
# A latent variable takes the scale of its first indicator whose loading
# is fixed at a value other than zero (its marker, so that a loading
# fixed at 1 stays 1 and the factor's variance carries the marker's
# units), resolved through higher order factors, and 1 when no loading is
# fixed (std.lv = TRUE, whose fixed variance of 1 then stays 1). Rows that
# are not parameters (:=, ==, <, >) get 1. Returns NULL when two rows that
# share a free parameter would get different factors.
.fmi_rescale_mono <- function(pt, cj) {
  n <- length(pt$lhs)
  mono <- rep(1, n)
  group <- if (is.null(pt$group)) rep(1L, n) else pt$group
  is_par <- !pt$op %in% c(":=", "==", "<", ">")
  for (g in unique(group[is_par])) {
    rows <- which(is_par & group == g)
    load <- rows[pt$op[rows] == "=~"]
    lv <- unique(pt$lhs[load])
    marker <- vapply(lv, function(l) {
      fixed <- load[pt$lhs[load] == l & pt$free[load] == 0L &
                      !is.na(pt$est[load]) & pt$est[load] != 0]
      if (length(fixed)) pt$rhs[fixed[1L]] else NA_character_
    }, character(1L))
    s_lv <- stats::setNames(rep(NA_real_, length(lv)), lv)
    s_lv[is.na(marker)] <- 1
    # A factor whose marker is itself a factor takes that factor's scale
    # once it is known; length(lv) passes resolve any chain.
    for (pass in seq_along(lv)) {
      for (l in names(s_lv)[is.na(s_lv)]) {
        mk <- marker[[l]]
        s_lv[[l]] <- if (mk %in% names(cj)) {
          cj[[mk]]
        } else if (mk %in% names(s_lv)) {
          s_lv[[mk]]
        } else {
          1
        }
      }
      if (!anyNA(s_lv)) break
    }
    s_lv[is.na(s_lv)] <- 1
    scale_of <- function(v) {
      if (v %in% names(cj)) {
        cj[[v]]
      } else if (v %in% names(s_lv)) {
        s_lv[[v]]
      } else {
        1
      }
    }
    for (i in rows) {
      mono[i] <- switch(pt$op[i],
        "=~" = scale_of(pt$rhs[i]) / scale_of(pt$lhs[i]),
        "~"  = scale_of(pt$lhs[i]) / scale_of(pt$rhs[i]),
        "~~" = scale_of(pt$lhs[i]) * scale_of(pt$rhs[i]),
        "~1" = scale_of(pt$lhs[i]),
        1)
    }
  }
  free <- pt$free
  shared <- unique(free[free > 0L][duplicated(free[free > 0L])])
  for (id in shared) {
    if (length(unique(mono[free == id])) > 1L) return(NULL)
  }
  mono
}

# .fmi_rescale_pt(start_pt, mono)
#
# The parameter table of the refits in the metric of .fmi_refit_scale().
# Each parameter row's starting value (the fit's estimate, or the fixed
# value) and any bound is divided by its factor. Each expression of a
# defined parameter (:=) and of a constraint (==, <, >) is rewritten so
# that it reads every parameter in the original metric: a label whose row
# has the factor m = 2^k becomes (label * 2^k), which is exact. A defined
# parameter then takes its original value, and lavaan's delta method
# standard error for it, from the refit's covariance matrix and the
# Jacobian of the rewritten expression, is its standard error in the
# original metric (by the chain rule); a constraint holds in the new
# metric exactly when it held in the original, whatever the factors of
# the parameters it ties. Symbols are bound as lavaan binds them
# (lav_pt_con_def() and lav_pt_con_ceq() in lavaan 0.7-3): a label by the
# first row that carries it, a plabel in a constraint by its own row, and
# a name defined by := not at all, since its value is already in the
# original metric. An expression that does not parse is left as it is.
.fmi_rescale_pt <- function(start_pt, mono) {
  out <- start_pt
  n <- length(out$lhs)
  is_par <- !out$op %in% c(":=", "==", "<", ">")
  for (nm in intersect(c("ustart", "lower", "upper"), names(out))) {
    out[[nm]][is_par] <- out[[nm]][is_par] / mono[is_par]
  }
  defined <- out$lhs[out$op == ":="]
  label <- if (is.null(out$label)) rep("", n) else out$label
  plabel <- if (is.null(out$plabel)) rep("", n) else out$plabel
  factor_of <- function(s, constraint) {
    if (s %in% defined) return(1)
    row <- if (constraint && s %in% plabel) {
      match(s, plabel)
    } else {
      match(s, label)
    }
    if (is.na(row) || !is_par[row]) 1 else mono[row]
  }
  rewrite <- function(text, constraint) {
    expr <- tryCatch(parse(text = text, keep.source = FALSE),
                     error = function(e) NULL)
    if (length(expr) != 1L) return(text)
    syms <- all.vars(expr)
    m <- vapply(syms, factor_of, numeric(1L), constraint = constraint)
    if (!length(syms) || all(m == 1)) return(text)
    swap <- lapply(seq_along(syms), function(i) {
      if (m[i] == 1) return(as.name(syms[i]))
      call("(", call("*", as.name(syms[i]), call("^", 2, log2(m[i]))))
    })
    names(swap) <- syms
    paste(deparse(do.call(substitute, list(expr[[1L]], swap)),
                  width.cutoff = 500L), collapse = " ")
  }
  for (i in which(out$op == ":=")) {
    out$rhs[i] <- rewrite(out$rhs[i], constraint = FALSE)
  }
  for (i in which(out$op %in% c("==", "<", ">"))) {
    out$lhs[i] <- rewrite(out$lhs[i], constraint = TRUE)
    out$rhs[i] <- rewrite(out$rhs[i], constraint = TRUE)
  }
  out
}

# .fmi_em_tol(x)
#
# The absolute tolerance to which the FIML refit runs lavaan's EM estimate
# of the unrestricted moments, for the data x of the refit (one column
# per observed variable). lavaan's EM stops when no mean or covariance
# changes by more than an absolute tol, 1e-5 by default, and EM converges
# linearly, at a rate set by the largest fraction of missing information
# of the unrestricted model, so where that fraction is high the moments
# it stops at are still short of the maximum, and the FMIs computed from
# them move with the tolerance: with x5 missing at random on x4 and a
# quarter of x9 missing in lavaan's nine Holzinger and Swineford tests,
# the Hessian FMIs of the three factor model moved by 1.9e-5 and the
# unstructured ones by 5.9e-6 when the EM was run to 1e-14 in place of
# 1e-5, and in a two group fit of the model (whose solution has a
# negative residual variance) by up to .0032. The tolerance here is 1e-10
# times the smallest variance among the variables, a relative 1e-10 for
# every mean and covariance, and in the metric of the refits, where every
# variance is near 1, the FMIs of those fits then moved by at most 9.3e-8
# between it and 1e-14 (the Hessian and structured FMIs of the two group
# fits; 2e-9 for the others). It is held above the rounding error of the
# moments, which EM cannot get below and which grows with the squares of
# the largest mean and variance, by a floor of 1e-16 times the larger of
# the two: lavaan's EM did not reach 1e-10 in 10,000 iterations with a
# variable's mean 10,000 times its standard deviation, and reached it at
# 1000. The floor binds only for a variable far from zero relative to its
# spread, or for a model kept in the variables' own metric whose
# variances differ by orders of magnitude.
.fmi_em_tol <- function(x) {
  x <- x[vapply(x, is.numeric, logical(1L))]
  v <- vapply(x, stats::var, numeric(1L), na.rm = TRUE)
  m <- abs(colMeans(x, na.rm = TRUE))
  v <- v[is.finite(v) & v > 0]
  m <- m[is.finite(m)]
  if (!length(v)) return(1e-10)
  max(1e-10 * min(v), 1e-16 * max(c(v, m^2)))
}

# .fmi_and(x)
#
# Names joined as English lists them: "a", "a and b", "a, b, and c".
.fmi_and <- function(x) {
  x <- as.character(x)
  n <- length(x)
  if (n <= 1L) return(paste(x, collapse = ""))
  if (n == 2L) return(paste(x, collapse = " and "))
  paste0(paste(x[-n], collapse = ", "), ", and ", x[n])
}

# The smallest eigenvalue of the correlation matrix of an unrestricted
# covariance estimate below which the estimate counts as singular: lavaan's
# own tolerance for a (near) singular EM covariance matrix (the non_pd_tol
# element of em_h1_args, 1e-5), taken on the correlation scale so that it
# does not depend on the units. The estimates that drift toward a singular
# matrix (.fmi_unbounded_subset()) had smallest eigenvalues from 2e-8 to
# 2e-6 when lavaan's EM stopped at its default tolerance, and the ones that
# settled at a maximum .016 or more, in one factor models of 6 and 12 items
# at N = 30 and 50 with 30 percent of the values of half the items missing
# (130 data sets); one data set whose seven complete cases nearly lay on a
# hyperplane had a maximum at 8e-6, which the EM reached.
.FMI_SINGULAR_TOL <- 1e-5

# .fmi_min_corr_eigen(S)
#
# The smallest eigenvalue of the correlation matrix of the covariance
# matrix S (0 when S has a variance that is not positive or a value that
# is not finite).
.fmi_min_corr_eigen <- function(S) {
  S <- as.matrix(S)
  d <- diag(S)
  if (!length(d) || !all(is.finite(S)) || any(!is.finite(d) | d <= 0)) {
    return(0)
  }
  R <- S / sqrt(outer(d, d))
  min(eigen(R, symmetric = TRUE, only.values = TRUE)$values)
}

# .fmi_unbounded_subset(x)
#
# A sufficient affine-collapse certificate for an unbounded unrestricted
# normal likelihood, not a count rule or a necessary-and-sufficient test.
# For a subset T, ALL rows observing T must lie in a hyperplane whose
# normal has support exactly T. Shrinking covariance along that normal
# makes their densities unbounded while every margin missing a T member
# remains positive definite. A rank deficit alone is insufficient: its
# nullspace may have forced zero coordinates, and rows observing the
# smaller support can disprove the relation. Reduce to the union of the
# nullspace's coordinate support and recollect rows before trying again.
# Candidates are the full set, each observed pattern, and singletons;
# absence of a certificate does not prove existence or uniqueness.
.fmi_unbounded_subset <- function(x) {
  x <- as.matrix(x)
  x <- x[rowSums(!is.na(x)) > 0L, , drop = FALSE]
  obs <- !is.na(x)
  if (!nrow(x) || !ncol(x)) return(NULL)
  scales <- vapply(seq_len(ncol(x)), function(j) {
    z <- x[obs[, j], j]
    if (!length(z)) return(1)
    z <- max(abs(z - mean(z)))
    if (is.finite(z) && z > 0) z else 1
  }, numeric(1))
  cand <- unique(rbind(rep(TRUE, ncol(obs)), obs, diag(TRUE, ncol(obs))))
  for (k in seq_len(nrow(cand))) {
    set <- cand[k, ]
    repeat {
      size <- sum(set)
      if (!size) break
      rows <- rowSums(obs[, set, drop = FALSE]) == size
      n_set <- sum(rows)
      if (!n_set) break
      values <- x[rows, set, drop = FALSE]
      z <- sweep(values, 2L, values[1L, ])
      z <- sweep(z, 2L, scales[set], "/")
      if (n_set > 1L) z <- z[-1L, , drop = FALSE]
      decomp <- svd(z, nu = 0L, nv = size)
      tol <- 100 * .Machine$double.eps * max(dim(z)) *
        max(c(1, decomp$d))
      # Only the null dimensions guaranteed by the number of affine
      # differences can support the dimension-based certificate. Extra
      # small singular values do not establish exact rank deficiency.
      rank <- if (n_set <= size) n_set - 1L else sum(decomp$d > tol)
      if (rank == size) break
      null <- decomp$v[, seq.int(rank + 1L, size), drop = FALSE]
      support <- sqrt(rowSums(null^2)) > sqrt(.Machine$double.eps) * size
      if (all(support)) {
        # Verify an explicit well-separated full-support direction.
        # Deterministic combinations avoid consuming the caller's RNG.
        coefficients <- cbind(diag(ncol(null)),
          cos(outer(seq_len(ncol(null)), seq_len(64L))))
        directions <- null %*% coefficients
        directions <- sweep(directions, 2L, sqrt(colSums(directions^2)), "/")
        clearance <- apply(abs(directions), 2L, min)
        direction <- directions[, which.max(clearance)]
        residual <- max(abs(z %*% direction))
        if (max(clearance) <= sqrt(.Machine$double.eps) * size ||
            residual > tol * max(1, sqrt(sum(direction^2)))) break
        return(list(vars = colnames(x)[set], n = n_set,
                    nullspace = null, scales = scales[set],
                    direction = direction / scales[set],
                    relative_residual = residual / max(1, max(abs(z))),
                    certified_by_dimension = n_set <= size,
                    certificate = if (n_set <= size) "full-support affine collapse"
                      else "affine relation to numerical precision"))
      }
      set[which(set)[!support]] <- FALSE
    }
  }
  NULL
}

.fiml_likelihood_obstacle <- function(x) {
  certificate <- .fmi_unbounded_subset(x)
  if (is.null(certificate)) return(NULL)
  if (!isTRUE(certificate$certified_by_dimension)) {
    return(paste0("The unrestricted covariance is numerically unstable: the ",
                  certificate$n, " rows observing ", .fmi_and(certificate$vars),
                  " have an affine relation to numerical precision. This ",
                  "rank diagnostic does not prove an unbounded likelihood; ",
                  "a finite nonsingular maximum has not been established."))
  }
  paste0("The unrestricted normal likelihood has no finite maximum: the ",
         certificate$n, " rows observing ", .fmi_and(certificate$vars),
         " have an affine-collapse direction supported on every variable ",
         "of that set, to numerical precision. Shrinking covariance in ",
         "that direction makes their likelihood unbounded while all other ",
         "observed margins remain nonsingular. More complete rows alone ",
         "are not the criterion; the observed values and patterns determine ",
         "this diagnosis.")
}

.fiml_data_h1_diagnostic <- function(x) {
  pairs <- crossprod(!is.na(as.matrix(x)))
  absent <- which(upper.tri(pairs) & pairs == 0L, arr.ind = TRUE)
  if (nrow(absent)) {
    nms <- colnames(x) %||% paste0("V", seq_len(ncol(x)))
    return(paste0("The unrestricted H1 covariance is not identified: no ",
      "case observes both members of ", paste(apply(absent, 1L, function(ij) {
        paste(nms[ij], collapse = " and ")
      }), collapse = "; "), ". A constrained model can identify these ",
      "relations through its measurement structure; ordinary H1 fit ",
      "diagnostics remain unavailable."))
  }
  weak <- tryCatch({
    .stop_single_valued_together(x, "The unrestricted-moment estimability screen")
    NULL
  }, error = conditionMessage)
  if (!is.null(weak)) return(paste0(weak, " This is a conservative sample ",
    "information screen, not proof of likelihood nonexistence."))
  .fiml_likelihood_obstacle(x)
}

.stop_fiml_unbounded <- function(x, method = "Full information maximum likelihood") {
  certificate <- .fmi_unbounded_subset(x)
  if (!is.null(certificate)) {
    why <- .fiml_likelihood_obstacle(x)
    .dmar_refuse(
      paste0(method, ": ", why),
      if (isTRUE(certificate$certified_by_dimension)) "fiml_unbounded" else
        "fiml_unstable",
      diagnostic = certificate)
  }
  invisible(NULL)
}

# `uncertainty = FALSE` is for a result that computed none
# (fa_exploratory(point_estimates_only = TRUE)), whose note then speaks of
# the estimates alone.
.relia_h1_qualify <- function(out, diagnostic, uncertainty = TRUE) {
  if (is.null(diagnostic)) return(out)
  note <- paste0(if (uncertainty) {
                   "The constrained model's estimate and H0-based uncertainty are "
                 } else {
                   paste0("The constrained model's estimates (no uncertainty was ",
                          "computed: point_estimates_only = TRUE) are ")
                 },
                 "retained, subject to its convergence, identification, and ",
                 "proper-solution checks. Unrestricted-H1 diagnostics and ",
                 "residual correlations are unavailable. ", diagnostic)
  attr(out, "h1_diagnostic") <- note
  .dmar_warn(
    note,
    "h1_unavailable")
  out
}

# The unrestricted normal H1 likelihood of the actual raw-data fit.
# Ordered-data and summary-statistic models have different H1 contracts.
.fiml_fit_h1_diagnostic <- function(fit, groups = NULL) {
  options <- lavaan::lavInspect(fit, "options")
  if (!options$missing %in% c("ml", "ml.x", "two.stage", "robust.two.stage") ||
      length(lavaan::lavNames(fit, "ov.ord"))) return(NULL)
  data <- tryCatch(lavaan::lavInspect(fit, "data"), error = function(e) NULL)
  if (is.null(data)) return(NULL)
  if (is.matrix(data) || is.data.frame(data)) data <- list(data)
  # X moments are not free in a fixed/conditional-X H1 model. A collapse
  # supported only on X does not prove that its conditional likelihood is
  # unbounded. Screen responses only: a response-only certificate remains
  # sufficient by setting regression slopes to zero and collapsing the
  # residual covariance. Exact Y-on-X relations and design rank are outside
  # this deliberately incomplete screen.
  fixed_x <- (isTRUE(options$fixed.x) || isTRUE(options$conditional.x)) &&
    length(lavaan::lavNames(fit, "ov.x")) > 0L
  predictors <- if (fixed_x) lavaan::lavNames(fit, "ov.x") else character()
  if (options$missing %in% c("two.stage", "robust.two.stage")) {
    aux <- .ts_fit_aux(fit)
    if (!is.null(aux)) {
      if (length(aux$data) != length(data) ||
          any(vapply(seq_along(data), function(g) {
            nrow(data[[g]]) != nrow(aux$data[[g]])
          }, logical(1L)))) {
        .dmar_refuse_argument(c(
          "The two-stage item and auxiliary rows cannot be aligned for ",
          "the unrestricted first-stage diagnostic. Refit from the ",
          "original raw data before reporting this result."), arg = NULL)
      }
      data <- Map(cbind, data, aux$data)
    }
  }
  labels <- lavaan::lavInspect(fit, "group.label")
  selected <- if (is.null(groups)) seq_along(data) else
    intersect(groups, seq_along(data))
  issues <- lapply(selected, function(g) {
    x <- data[[g]]
    if (length(predictors)) {
      x <- x[, !colnames(x) %in% predictors, drop = FALSE]
    }
    if (!ncol(x)) return(NULL)
    why <- .fiml_data_h1_diagnostic(x)
    if (is.null(why)) return(NULL)
    if (fixed_x) why <- paste0(why, " This sufficient screen used response ",
      "variables only; fixed-predictor design rank and exact response-on-",
      "predictor relations were not assessed.")
    if (length(data) > 1L) paste0("Group ", labels[g], ": ", why) else why
  })
  issues <- unlist(issues, use.names = FALSE)
  if (length(issues)) paste(issues, collapse = " ") else NULL
}

.fiml_require_h1 <- function(fit, context, observed = FALSE, groups = NULL) {
  options <- lavaan::lavInspect(fit, "options")
  if (!observed && !options$missing %in%
      c("two.stage", "robust.two.stage")) return(invisible(NULL))
  diagnostic <- .fiml_fit_h1_diagnostic(fit, groups = groups)
  if (!is.null(diagnostic)) {
    .dmar_refuse_argument(c(
      context, " requires eligible unrestricted moments. ", diagnostic,
      " Use a supported constrained FIML analysis when its H0 estimates ",
      "remain eligible, or revise the data/model before using these ",
      "moments."), arg = NULL)
  }
  invisible(NULL)
}

.fiml_check_boot_aux <- function(fit, context) {
  if (.ts_is_two_stage_fit(fit) && !is.null(.ts_fit_aux(fit))) {
    .dmar_refuse_argument(c(
      context, " cannot yet carry separately stored two-stage auxiliary ",
      "variables into every resampled refit. An interval that omits them ",
      "would describe a different analysis. Use an eligible delta ",
      "interval where available, or explicitly resample and refit the ",
      "original items and auxiliaries together."), arg = NULL)
  }
  invisible(NULL)
}

# Recognize the supported H0-based covariance routes. Custom information
# settings are not silently interpreted as the default Hessian route.
.fiml_h0_vcov_supported <- function(fit) {
  options <- lavaan::lavInspect(fit, "options")
  first <- function(x) if (length(x)) x[1L] else NA_character_
  se <- first(options$se)
  if (identical(se, "none")) return(TRUE)
  if (!first(options$information.bread) %in% c("default", NA_character_)) {
    return(FALSE)
  }
  bread <- identical(first(options$h1.information), "structured") ||
    (identical(first(options$information), "observed") &&
       identical(first(options$observed.information), "hessian"))
  meat <- identical(se, "standard") ||
    (identical(se, "robust.huber.white") &&
       identical(first(options$information.meat), "first.order") &&
       identical(first(options$h1.information.meat), "structured"))
  isTRUE(bread && meat)
}

# .fmi_check_saturated(covs, dat, group_labels, capped)
#
# Every FMI needs unrestricted moments for its evaluation. A sufficient
# full-support affine-collapse certificate or a numerical affine
# instability is refused even when EM settles at a well-conditioned local
# stationary point. An iteration-limited near-singular iterate without a
# certificate is a fitting failure, not proof of nonexistence. covs and dat
# hold each group's current moments and actual observed rows. The legacy
# condition class remains for compatibility; the text distinguishes the
# diagnosis. Absence of either obstacle does not prove existence or a
# unique global maximum.
.fmi_check_saturated <- function(covs, dat, group_labels = NULL,
                                 capped = FALSE) {
  for (g in seq_along(covs)) {
    ev <- .fmi_min_corr_eigen(covs[[g]])
    unbounded <- .fmi_unbounded_subset(dat[[g]])
    near_singular <- !is.finite(ev) || ev < .FMI_SINGULAR_TOL
    if (is.null(unbounded) && !(near_singular && capped)) next
    n_cases <- sum(rowSums(!is.na(as.matrix(dat[[g]]))) > 0L)
    p <- NCOL(dat[[g]])
    where <- if (!is.null(group_labels)) {
      paste0(" in group ", group_labels[g])
    } else {
      ""
    }
    why <- if (!is.null(unbounded)) {
      k <- length(unbounded$vars)
      paste0(" Only ", unbounded$n, " of the ", n_cases, " cases", where,
             " observe ", if (k == p) paste0("all ", p, " variables") else
               paste0("all of ", .fmi_and(unbounded$vars)), ". ",
             .fiml_likelihood_obstacle(dat[[g]]))
    } else {
      paste0(" lavaan's EM algorithm reached its iteration limit at a near-singular ",
             "iterate. A finite maximum has not been established; this is a ",
             "fitting failure, not a proof that the likelihood is unbounded.")
    }

    .dmar_refuse(
      paste0(
      "The unrestricted estimate of the means and covariances, at which ",
      "the fraction of missing information is evaluated, is unavailable",
      where, " (the smallest eigenvalue of the current correlation estimate is ",
      format(ev, digits = 2L), ").", why, " A fraction computed where ",
      "the algorithm happens to stop would depend on where that is, so ",
      "none is computed. The model's own estimates and standard errors ",
      "do not use this estimate (under lavaan's default information); its ",
      "chi-square test and the fit indices built on the unrestricted ",
      "model's likelihood do."),
      "fmi_saturated_singular")
  }
  invisible(TRUE)
}

# .fmi_h1_covs(fit, ng)
#
# The fit's own unrestricted covariance estimates, one per group (lavaan's
# EM at the fit's settings), or NULL when the fit has none (h1 = FALSE).
.fmi_h1_covs <- function(fit, ng) {
  h1 <- tryCatch(lavaan::lavInspect(fit, "sampstat.h1"),
                 error = function(e) NULL)
  if (is.null(h1)) return(NULL)
  if (ng == 1L) h1 <- list(h1)
  covs <- lapply(h1, function(h) if (is.list(h)) h$cov else NULL)
  if (any(vapply(covs, is.null, logical(1L)))) return(NULL)
  covs
}

# .fmi_core(fit, information)
#
# The fraction of missing information for every row of the fit's parameter
# table. Two refits are made at the fit's estimates with the optimizer
# switched off, both under normal theory ML standard errors with the same
# information type:
#
#   the FIML side, to the fit's own data with its missing data option,
#     whose standard errors are SE_O; and
#   the complete data side, to the EM (unrestricted h1) means and
#     covariances of the data, with the same N and no rescaling of the
#     covariance matrix, whose standard errors are SE_C.
#
# FMI = 1 - SE_C^2 / SE_O^2 (Savalei & Rhemtulla, 2012, Eq. 4). This is
# what lavaan's parameterEstimates(fmi = TRUE) computes for an ML fit
# under its own information settings; doing it here lets the information
# be chosen and lets an MLR fit be evaluated at its (identical) ML
# estimates. A defined (:=) parameter gets its standard errors from the
# delta method under each fit, so its FMI uses the same ratio.
#
# Both refits are made in a metric in which every observed variable's
# standard deviation is near 1 (.fmi_refit_scale(), .fmi_rescale_pt()),
# with their standard errors mapped back exactly. Three things depend on
# that metric. (1) The unrestricted moments, at which the unstructured
# information is evaluated and to which the complete data side is fitted,
# are the FIML refit's own EM estimates. lavaan's EM stops when no mean or
# covariance changes by more than an absolute 1e-5, so in the raw metric
# of variables measured in small units it stops short of the maximum: in a
# two group regression whose variables had SDs near .003, the slopes' FMIs
# were off by up to .022, and lavaan does not rescale a refit with the
# optimizer off (it rescales its own fits only for SDs of .001 or less or
# 1000 or more, or a ratio of 1000). (2) The information matrix is
# inverted, which fails in floating point once the variances differ by
# enough orders of magnitude: with one variable's SD near 7e7 lavaan could
# not invert it, and every FMI was NA. (3) The Hessian estimator
# differentiates numerically with lavaan's fixed step. In the new metric,
# with the EM run to a relative 1e-10 (.fmi_em_tol()), none of the three
# depends on the units: with the estimates held, the FMIs of the three
# estimators agreed to within 2e-9 when every variable was multiplied by a
# constant from .001 to 1e6. What remains is the precision of the fit's own estimates, which
# lavaan computes in the variables' metric; .fmi_check_maximum() warns
# when they are not at the maximum. An earlier version handed the
# refit the source fit's own moments instead (its h1 slot, through
# lavaan's sloth1 argument). lavaan holds them by position, in the fit's
# order of the variables, and a fit made with ov.order = "data" orders
# them as the data do while the refit orders them as the model does, so
# the refit read a permuted covariance matrix (FMIs below -1, with no
# warning). The refit now estimates the moments itself, in its own order.
#
# Returns the fit's parameter table with columns se_observed, se_complete,
# and fmi added, and the attributes sd_range (the range of the variables'
# SDs in the metric of the refits, which says whether the Hessian
# estimator's step suits them) and scale_reason (why the refits stayed in
# the raw metric, or NULL). Stops when a pair of variables is never
# observed together and when the unrestricted estimate of the means and
# covariances does not exist (.fmi_check_saturated()), both before the
# refits where the data or the source fit show it; when the refit moved
# the estimates (the guard that catches a solution lavaan alters after
# estimation, such as a rotation); and when a refit's information matrix
# cannot be inverted. Warns when the refit's EM stopped at its iteration
# limit short of an estimate that exists, and when the fit's estimates are
# not at the maximum of the likelihood.
.fmi_core <- function(fit, information, confirm_location = TRUE) {
  opt <- lavaan::lavInspect(fit, "options")
  pt <- lavaan::parTable(fit)
  start_pt <- pt
  start_pt$ustart <- start_pt$est
  start_pt$start <- NULL
  start_pt$est <- NULL
  start_pt$se <- NULL

  ng <- lavaan::lavInspect(fit, "ngroups")
  dat <- lavaan::lavInspect(fit, "data")
  if (ng == 1L) dat <- list(dat)
  group_labels <- if (ng > 1L) lavaan::lavInspect(fit, "group.label")

  # Every pair of variables must be observed together somewhere. The
  # unrestricted estimate behind every FMI has no information about the
  # covariance of a pair that no case observes together (a two-form
  # planned missing design, say), and lavaan's EM leaves it at its
  # starting value, 0, while lavaan only warns that it is not identified:
  # with t1 and t6 of a two-form design on seven Holzinger and Swineford
  # tests, replacing that 0 by other values the data cannot tell apart (the
  # observed data log likelihood was the same to every digit) moved the
  # unstructured FMI of a loading by up to .12. The FIML fit of a model
  # that identifies the covariance through its structure is unaffected,
  # so only the fraction is refused, with the pairs named, as the
  # reliability family and cov_fiml() refuse such data.
  for (g in seq_len(ng)) {
    .stop_never_observed_together(
      dat[[g]],
      method = paste0("the fraction of missing information (through the ",
                      "unrestricted estimate of the means and covariances ",
                      "at which it is evaluated)",
                      if (ng > 1L) paste0(" in group ", group_labels[g])),
      remedy = paste0(
        "The model's own FIML estimates do not need that covariance where ",
        "the model identifies it, as in a planned missing design, and they ",
        "stand; the fraction of missing information cannot be computed ",
        "from these data. missing_coverage() reports the cases behind ",
        "every pair."),
      class = "dmar_fmi_never_together")
  }

  # The unrestricted estimate of the means and covariances must exist. The
  # source fit's own estimate (lavaan's EM at the fit's tolerance) is
  # checked first, at no cost: an estimate drifting toward a singular
  # matrix is already near singular there, and the refit's EM, run to a
  # far tighter tolerance, would otherwise spend its 10,000 iterations
  # following the drift (13 seconds, against .2 for the whole computation
  # where the estimate exists, on the data of the comment above
  # .fmi_check_saturated()). The refit's own estimate is checked again
  # below.
  src_covs <- .fmi_h1_covs(fit, ng)
  if (!is.null(src_covs)) .fmi_check_saturated(src_covs, dat, group_labels)

  scaling <- .fmi_refit_scale(pt, dat)
  cj <- scaling$cj
  mono <- scaling$mono
  if (!is.null(cj)) start_pt <- .fmi_rescale_pt(start_pt, mono)
  to_refit_metric <- function(x) {
    x <- as.data.frame(x)
    for (v in intersect(names(x), names(cj))) x[[v]] <- x[[v]] / cj[[v]]
    x
  }

  # Every refit holds the estimates, whatever lavaan's checks of a
  # solution would say of them, so those checks are left off; a solution
  # lavaan will not take as a starting point is caught by the guard below.
  # check.start = FALSE keeps lavaan from shrinking a covariance whose
  # starting value implies a correlation above 1, which it otherwise does
  # for any fit: a factor correlation above 1 (a latent covariance matrix
  # that is not positive definite), a Heywood case on an indicator with a
  # correlated residual (lavaan resets both variances to 1 and the
  # covariance to 0), and the auxiliary variables of lavaan's aux argument
  # (fa_confirmatory(..., aux = )), whose covariances with an indicator run entirely
  # through its residual and so can exceed what the residual variance
  # allows of a correlation while every variance is positive and the
  # implied covariance matrix of the observed variables is positive
  # definite. lavaan's post-estimation check flags the auxiliary fit,
  # whose augmented residual covariance matrix is then not positive
  # definite (check.post = FALSE; the reliability family's auxiliary fits
  # leave it off for the same reason, .omega_fit_cfa_aux() in
  # R/reliability_internals.R). check.vcov = FALSE: whether the model is
  # identified is the source fit's to report, and lavaan's check compares
  # the smallest eigenvalue of the covariance matrix of the estimates with
  # an absolute 1.8e-12, which variables in small or large units cross.
  common <- c(list(model = start_pt, estimator = "ML", se = "standard",
                   information = "observed", optim.method = "none",
                   check.gradient = FALSE, check.post = FALSE,
                   check.start = FALSE, check.vcov = FALSE,
                   meanstructure = TRUE, fixed.x = opt$fixed.x,
                   conditional.x = FALSE, test = "none",
                   baseline = FALSE),
              .fmi_information_options(information))

  if (ng == 1L) {
    fiml_data <- to_refit_metric(dat[[1L]])
    fiml_args <- list(data = fiml_data, missing = opt$missing, h1 = TRUE)
  } else {
    group_var <- ".dmar_fmi_group."
    pieces <- lapply(seq_len(ng), function(g) {
      x <- to_refit_metric(dat[[g]])
      x[[group_var]] <- group_labels[g]
      x
    })
    cols <- lapply(pieces, names)
    if (!all(vapply(cols, identical, logical(1L), cols[[1L]]))) {
      .dmar_refuse_argument(c(
        "The groups of the fit do not share one set of variables; the ",
        "fraction of missing information is computed for multiple ",
        "group fits whose groups analyze the same variables."), arg = NULL)
    }
    fiml_data <- do.call(rbind, pieces)
    fiml_args <- list(data = fiml_data, missing = opt$missing, h1 = TRUE,
                      group = group_var, group.label = group_labels)
  }

  # The EM estimation of the unrestricted moments follows the source fit's
  # settings, except for when it stops (.fmi_em_tol()), so the complete
  # data side uses the moments it would have at the maximum.
  for (nm in c("em.h1.args", "h1.missing.method")) {
    if (!is.null(opt[[nm]])) fiml_args[[nm]] <- opt[[nm]]
  }
  em_args <- if (is.list(fiml_args$em.h1.args)) fiml_args$em.h1.args else
    list()
  em_args$tol <- min(c(em_args$tol,
                       .fmi_em_tol(fiml_data[setdiff(names(fiml_data),
                                                     ".dmar_fmi_group.")])))
  em_args$max_iter <- as.integer(max(c(em_args$max_iter, 10000L)))
  fiml_args$em.h1.args <- em_args

  # The source fit has already reported the rows it set aside (cases with
  # every variable missing, and under missing = "ml" with fixed.x = TRUE the
  # cases missing a fixed covariate), the pairs of variables that fewer
  # than a tenth of the cases observe together (a pair no case observes
  # together is refused above), and the sizes of the variances (above a
  # million, or differing by a factor of 1000); the refits to the same
  # data would say so again, so those notices are muffled here. lavaan's
  # warning that it could not compute standard errors is muffled too and
  # recorded: it says the model may not be identified, which for a refit at
  # a fit's estimates is seldom the reason, and the stop below says what
  # is. lavaan's notices about its EM estimate of the unrestricted moments
  # in the refit, that it stopped at the iteration limit (with the advice
  # to raise max_iter, which cannot help where the estimate does not
  # exist) and that its covariance matrix is (near) singular, are muffled
  # and recorded, and .fmi_check_saturated() and the warning below say
  # what they mean for the fraction.
  state <- new.env(parent = emptyenv())
  state$not_inverted <- FALSE
  state$em_capped <- FALSE
  refit_notice <- function(w) {
    msg <- gsub("[[:space:]]+", " ", conditionMessage(w))
    if (grepl("Could not compute standard errors", msg, fixed = TRUE)) {
      state$not_inverted <- TRUE
      invokeRestart("muffleWarning")
    }
    if (grepl("Maximum number of iterations reached", msg, fixed = TRUE) &&
        grepl("EM", msg, fixed = TRUE)) {
      state$em_capped <- TRUE
      invokeRestart("muffleWarning")
    }
    if (grepl("EM estimated variance-covariance matrix", msg,
              fixed = TRUE) &&
        grepl("singular", msg, fixed = TRUE)) {
      invokeRestart("muffleWarning")
    }
    if (grepl("some cases are empty", msg, fixed = TRUE) ||
        grepl("pairwise combinations have less than 10% coverage", msg,
              fixed = TRUE) ||
        grepl("deleted due to missing values in exogenous", msg,
              fixed = TRUE) ||
        grepl("observed variances are larger than 1000000", msg,
              fixed = TRUE) ||
        grepl("observed variances in the sample covariance matrix are larger",
              msg, fixed = TRUE) ||
        grepl("factor 1000 times larger", msg, fixed = TRUE)) {
      invokeRestart("muffleWarning")
    }
  }
  fit_observed <- withCallingHandlers(
    do.call(lavaan::lavaan, c(common, fiml_args)),
    warning = refit_notice)

  # The EM moments are the FIML refit's own estimates of the unrestricted
  # model (h1 = TRUE), named in its order of the variables.
  h1 <- lavaan::lavInspect(fit_observed, "sampstat.h1")
  if (ng == 1L) {
    em_cov <- h1$cov
    em_mean <- h1$mean
  } else {
    em_cov <- lapply(h1, `[[`, "cov")
    em_mean <- lapply(h1, `[[`, "mean")
  }
  # The refit's estimate, run to the tighter tolerance, can show a drift the
  # source fit's did not, or stop at the iteration limit on the way.
  .fmi_check_saturated(if (ng == 1L) list(em_cov) else em_cov, dat,
                       group_labels, capped = state$em_capped)
  if (state$em_capped) {
    # An estimate that exists but that EM approaches slowly, as it does
    # when much of the information about the unrestricted moments is
    # missing (the rate of EM is the largest fraction of missing
    # information of the unrestricted model): in one data set of the
    # calibration design (12 items, N = 50, 30 percent of the values of
    # half of them missing) lavaan's EM needed about 12,500 iterations,
    # and the FMIs at the 10,000th were within 3e-7 of the converged ones
    # (and within .002 at lavaan's default tolerance). A fit made with a
    # larger max_iter in em_h1_args passes it to the refit.
    .dmar_warn(c(
      "lavaan's EM estimate of the unrestricted means and ",
      "covariances, at which the fraction of missing information is ",
      "evaluated, stopped at its limit of ", fiml_args$em.h1.args$max_iter,
      " iterations before it met its tolerance; EM converges slowly ",
      "when much of the information about those moments is missing. ",
      "The fractions are computed at its last iterate, and can differ ",
      "slightly from those at the converged estimate. To let it run ",
      "further, fit the model with lavaan's em_h1_args = ",
      "list(max_iter = ", 5L * fiml_args$em.h1.args$max_iter,
      "), which the computation then uses."), "result_qualified")
  }
  fit_complete <- withCallingHandlers(
    do.call(lavaan::lavaan, c(common, list(
      sample.cov = em_cov, sample.mean = em_mean,
      sample.nobs = lavaan::lavInspect(fit, "nobs"),
      sample.cov.rescale = FALSE, h1 = FALSE))),
    warning = refit_notice)

  pt_o <- lavaan::parTable(fit_observed)
  pt_c <- lavaan::parTable(fit_complete)
  same_rows <- function(a) {
    length(a$lhs) == length(start_pt$lhs) &&
      identical(paste(a$lhs, a$op, a$rhs), paste(start_pt$lhs, start_pt$op,
                                                 start_pt$rhs))
  }
  if (!same_rows(pt_o) || !same_rows(pt_c)) {
    .dmar_refuse_argument(c(
      "The refits of the model at its estimates did not keep the rows ",
      "of its parameter table, so their standard errors cannot be ",
      "matched to the fit's parameters."), arg = NULL)
  }

  # Estimate invariance guard: with the optimizer off the refit must
  # return the estimates it was started at, or its standard errors would
  # belong to other estimates. The comparison is made in the metric of the
  # refits, where every parameter is of the order of the variables' unit
  # variances, so the relative tolerance means the same for each. Two
  # kinds of fit move here. A solution that lavaan alters after
  # estimation, such as a rotated exploratory one, is altered again. And
  # an improper solution can be refused as a starting point: lavaan 0.7-3
  # checks the refit's result for a residual variance more negative than
  # the variable's observed variance (a solution that has "run away"),
  # finds one in a Heywood case of that size, and retries from its simple
  # starting values (start = "simple": free loadings of .7, unit
  # variances, zero covariances and intercepts), which the refit then
  # returns. lavaan's own parameterEstimates(fmi = TRUE) makes the same
  # refit and, for such a fit, reports FMIs evaluated at those defaults
  # without saying so.
  free_rows <- which(start_pt$free > 0L)
  held <- start_pt$ustart[free_rows]
  moved <- abs(pt_o$est[free_rows] - held) > 1e-8 * (1 + abs(held))
  if (any(is.na(moved)) || any(moved)) {
    # Improper means a negative variance estimate, which is what the
    # message says. lavaan's post-estimation check is not the test: it also
    # fails for a proper fit with saturated correlates, whose residual
    # covariance matrix is not positive definite by construction.
    improper <- any(pt$op == "~~" & pt$lhs == pt$rhs & pt$free > 0L &
                      !is.na(pt$est) & pt$est < 0)
    .dmar_refuse_argument(c(
      "Refitting the model at its estimates moved them, so the ",
      "standard errors of the refit would not belong to the fit's ",
      "estimates; the fraction of missing information cannot be ",
      "computed for this fit. lavaan moves the estimates of a solution ",
      "it alters after estimation, such as a rotated exploratory one, ",
      "and of an improper solution, such as a Heywood case with a ",
      "residual variance far below zero, which it will not take as the ",
      "starting point of a refit and replaces with its own simple ",
      "starting values (lavaan's own parameterEstimates(fmi = TRUE) is ",
      "affected the same way).",
      if (improper) {
        paste0(" This fit's solution is improper (a negative variance ",
               "estimate); respecify the model so that the solution is ",
               "proper, and compute the fraction for that fit.")
      }), arg = NULL)
  }

  # The range of the variables' standard deviations in the metric of the
  # refits, under the unrestricted moments.
  em_var <- unlist(lapply(if (ng == 1L) list(em_cov) else em_cov, diag))
  em_var <- em_var[is.finite(em_var) & em_var > 0]
  sd_range <- if (length(em_var)) {
    range(sqrt(em_var))
  } else {
    c(NA_real_, NA_real_)
  }
  if (state$not_inverted) {
    # The location of each variable relative to its spread, under the
    # unrestricted moments (the largest over groups). The groups' moments
    # are taken unnamed: lavaan names them by group, and unlist() would
    # prefix every variable with its group ("a.x2"), which no regressor
    # name matches.
    loc <- unlist(lapply(if (ng == 1L) list(list(em_mean, em_cov)) else
      Map(list, unname(em_mean), unname(em_cov)), function(m) {
        mu <- as.numeric(m[[1L]])
        s <- sqrt(pmax(diag(as.matrix(m[[2L]])), 0))
        if (length(mu) != length(s) || is.null(rownames(m[[2L]]))) {
          return(numeric(0))
        }
        stats::setNames(abs(mu) / s, rownames(m[[2L]]))
      }))
    loc <- loc[is.finite(loc)]
    location <- if (length(loc)) {
      vapply(split(loc, names(loc)), max, numeric(1L))
    } else {
      numeric(0)
    }
    # Only a regressor's location can make the information singular: the
    # mean of an observed variable that predicts another enters that
    # variable's mean through its slope, so far from zero the slope and
    # the intercept are nearly collinear, while an indicator's or an
    # outcome's location enters only its own intercept (an indicator of
    # a three factor model 100,000 standard deviations from zero left the
    # information invertible, with FMIs within 2e-5 of the unshifted
    # ones). Whether location is the cause is then checked
    # (.fmi_location_confirms()), since a model that is not identified
    # fails whatever its variables' means.
    regressors <- unique(pt$rhs[pt$op == "~"])
    far <- location[names(location) %in% regressors & location >= 1000]
    cause <- NULL
    if (length(far) && confirm_location) {
      cause <- .fmi_location_confirms(fit, pt, opt, dat, names(far),
                                      information, group_labels)
    }
    .fmi_stop_not_inverted(pt, opt, information, sd_range,
                           rescaled = !is.null(cj),
                           reason = scaling$reason, location = far,
                           location_cause = cause)
  }

  # A maximum of the likelihood is checked for in the metric of the
  # refits, where every parameter is of the order of a unit variance;
  # in the variables' own metric of an unmapped model it would not be.
  # The Newton step is a screen (.fmi_check_maximum()); where it moves an
  # estimate, the FMIs are computed again one step on, at the maximum to
  # working precision, and the warning is given when they differ from the
  # FMIs at the fit's estimates by more than 1e-4 (.fmi_warn_maximum()).
  # The FMIs returned stay those of the fit's estimates.
  if (is.null(scaling$reason)) {
    screen <- .fmi_check_maximum(fit_observed, pt_o, mono, group_labels)
    if (!is.null(screen)) {
      fmi_here <- 1 - pt_c$se^2 / pt_o$se^2
      fmi_step <- tryCatch({
        stepped <- start_pt
        stepped$ustart[screen$rows] <- pt_o$est[screen$rows] + screen$step
        args_step <- common
        args_step$model <- stepped
        fo <- suppressWarnings(do.call(lavaan::lavaan,
                                       c(args_step, fiml_args)))
        fc <- suppressWarnings(do.call(lavaan::lavaan, c(args_step, list(
          sample.cov = em_cov, sample.mean = em_mean,
          sample.nobs = lavaan::lavInspect(fit, "nobs"),
          sample.cov.rescale = FALSE, h1 = FALSE))))
        po <- lavaan::parTable(fo)
        pc <- lavaan::parTable(fc)
        if (length(po$se) != length(pt_o$se) ||
            length(pc$se) != length(pt_c$se)) {
          .dmar_refuse_argument("rows", arg = NULL)
        }
        1 - pc$se^2 / po$se^2
      }, error = function(e) NULL)
      .fmi_warn_maximum(screen, fit, pt_o, fmi_here, fmi_step, em_cov, cj,
                        group_labels)
    }
  }

  # Back to the original metric: a parameter row's standard error times
  # its factor. A defined parameter's is already in that metric.
  se_o <- pt_o$se * mono
  se_c <- pt_c$se * mono
  se_c[!is.na(se_c) & se_c == 0] <- NA_real_
  pt$se_observed <- se_o
  pt$se_complete <- se_c
  pt$fmi <- 1 - se_c^2 / se_o^2
  attr(pt, "sd_range") <- sd_range
  attr(pt, "scale_reason") <- scaling$reason
  pt
}

# .fmi_stop_not_inverted(pt, opt, information, sd_range, rescaled, reason,
#                        location, location_cause)
#
# The stop for a refit whose information matrix lavaan could not invert,
# in place of lavaan's warning that the model may not be identified,
# which for a refit at a fit's estimates is seldom the reason. The error
# has class dmar_fmi_not_inverted. The words name the likeliest cause.
# (1) The refits stayed in the variables' own metric (reason says why)
# and the variables are badly scaled by lavaan's own rule (an SD of .001
# or less or 1000 or more, or two that differ by a factor of 1000): the
# failure is numerical, whatever happened to the fit's own standard
# errors, which lavaan computed in the same metric. (2) A regressor, an
# observed variable that predicts another, lies 1000 or more standard
# deviations from zero (location, the ratio of its absolute mean to its
# standard deviation under the unrestricted moments, for the regressors
# only), and the same model fit to the data with it centered has an
# information that can be inverted (location_cause$confirmed TRUE,
# .fmi_location_confirms()). Dividing by a power of two leaves that ratio
# as it is, and at such a location the intercepts and the coefficients of
# the variable are nearly collinear in the information matrix, which is
# then singular to working precision. That happens only well beyond the
# ratio of 1000: a predictor of a regression with standard deviation 1
# gave the unshifted FMIs at 5000 and 7000 standard deviations from zero
# and failed from about 8000 (two predictors so shifted, from about 6000),
# where lavaan's fit had no standard errors either. Centering changes no
# estimate and no FMI but those of intercepts. Where the check could not
# be made (confirmed NA: an intercept that centering would change is
# fixed or constrained, or the centered fit failed), both causes are
# named, since centering cannot repair a model that is not identified.
# Where the centered fit failed in the same way (confirmed FALSE),
# location is not the cause and the stop goes on as if the variable were
# near zero. (3) lavaan could not compute the fit's own standard errors
# either: the model may not be identified, or its information is singular
# to working precision for another reason. (4) Under the structured or
# Hessian information, which can be singular at a fit's estimates in
# small samples (Chen & Savalei, 2021), that information is named. (5)
# Otherwise the unstructured information, evaluated at the saturated
# estimates, is singular, which is the mark of a model that is not
# identified, empirically or by construction; lavaan can return standard
# errors for such a fit, very large ones (580 for a loading of .7 in a
# factor with two indicators).
.fmi_stop_not_inverted <- function(pt, opt, information, sd_range, rescaled,
                                   reason, location = numeric(0),
                                   location_cause = NULL) {
  src_se <- pt$se[pt$free > 0L]
  has_se <- !identical(opt$se, "none") && length(src_se) > 0L
  fit_no_se <- has_se && all(!is.finite(src_se))
  either <- if (fit_no_se) {
    " (lavaan could not compute the fit's own standard errors either)"
  }
  badly_scaled <- !rescaled && !is.null(reason) &&
    all(is.finite(sd_range)) &&
    (sd_range[1L] <= 0.001 || sd_range[2L] >= 1000 ||
       sd_range[2L] / sd_range[1L] >= 1000)
  far <- location[is.finite(location) & location >= 1000]
  confirmed <- if (is.null(location_cause)) NA else
    location_cause$confirmed
  if (identical(confirmed, FALSE)) far <- far[0L]
  where <- NULL
  if (length(far)) {
    far <- sort(far, decreasing = TRUE)
    one <- length(far) == 1L
    where <- paste0(
      " ", .fmi_and(names(far)), if (one) " lies " else " lie ",
      if (one) {
        format(far[[1L]], digits = 3L)
      } else {
        paste0("between ", format(min(far), digits = 3L), " and ",
               format(max(far), digits = 3L))
      },
      " standard deviations from zero, and so far from zero an intercept ",
      "and the coefficients of ", if (one) "that variable" else
        "those variables")
    center <- paste0(
      "subtract a constant near its mean from ",
      if (one) "that variable (center it)" else "each of them (center them)")
  }
  why <- if (badly_scaled) {
    paste0(" The variables' standard deviations run from ",
           format(sd_range[1L], digits = 3L), " to ",
           format(sd_range[2L], digits = 3L), ", and the refits could not ",
           "be moved to a metric in which they are near 1 because the ",
           "model has ", reason, ", so the failure is numerical", either,
           ". Divide the variables by constants that bring their standard ",
           "deviations nearer 1, which leaves the fraction unchanged, and ",
           "fit the model again.")
  } else if (length(far) && isTRUE(confirmed)) {
    paste0(where, " are nearly collinear, which makes the information ",
           "matrix singular to working precision", either, "; the same ",
           "model fit with ", if (length(far) == 1L) "it" else "them",
           " centered has an information that can be inverted. ",
           toupper(substr(center, 1L, 1L)), substring(center, 2L),
           ", which changes no estimate and no fraction of missing ",
           "information but those of intercepts, and fit the model again.")
  } else if (length(far)) {
    paste0(where, " can be nearly collinear, which makes the information ",
           "matrix singular to working precision", either, ", but whether ",
           "that is the cause here could not be checked, because ",
           location_cause$why %||% "the check was not made",
           ". If the model is identified, ", center, " and fit the model ",
           "again; centering cannot repair a model that is not identified.")
  } else if (fit_no_se) {
    paste0(" lavaan could not compute the fit's own standard errors ",
           "either, which happens when the model is not identified or its ",
           "information matrix is singular to working precision for ",
           "another reason.")
  } else if (!identical(information, "unstructured")) {
    paste0(" The ", information, " information can be singular at a fit's ",
           "estimates in small samples (Chen & Savalei, 2021); ",
           "information = \"unstructured\", the default, is evaluated at ",
           "the saturated estimates, where the information of an identified ",
           "model is not singular.")
  } else {
    paste0(" The unstructured information is evaluated at the saturated ",
           "estimates, where it is singular only for a model that is not ",
           "identified, empirically or by construction; lavaan can return ",
           "standard errors for such a fit, and they are then very large.")
  }
  .dmar_refuse(
    paste0(
      "The information matrix of a refit at the fit's estimates could ",
      "not be inverted, so the standard errors behind the fraction of ",
      "missing information, and the fraction itself, cannot be ",
      "computed.", why),
    "fmi_not_inverted")
}

# .fmi_location_confirms(fit, pt, opt, dat, vars, information,
#                        group_labels = NULL)
#
# Whether centering the regressors `vars`, which lie far from zero, makes
# the information invertible, the check behind naming their location as
# the cause of a refit whose information could not be inverted
# (.fmi_stop_not_inverted()). The model is fit again by FIML, with the
# optimizer on, to the data with each of them centered at its mean over
# all groups, starting from the fit's estimates carried to the centered
# data (an outcome's intercept plus the slope times the constant, the
# regressor's own mean minus it), and the fraction is computed for that
# fit (.fmi_core(), which then does not check again). Centering is a
# reparameterization only when every intercept it changes is free and
# constrained to nothing (no shared label or free parameter, no ==, <, >,
# or := expression reaching it), and the regressor's own mean is free or
# fixed by lavaan at the sample mean (fixed.x); a latent outcome's
# intercept, fixed at 0, or intercepts held equal across groups would
# make the centered model another model, so the check is not made there.
# Returns list(confirmed, why): confirmed TRUE when the centered fit's
# information can be inverted, FALSE when it cannot either, NA when the
# check could not be made, and why, the reason for NA, completing "could
# not be checked, because ...".
.fmi_location_confirms <- function(fit, pt, opt, dat, vars, information,
                                   group_labels = NULL) {
  n <- length(pt$lhs)
  grp <- if (is.null(pt$group)) rep(1L, n) else pt$group
  comp <- .fmi_equality_classes(pt)
  expr_rows <- which(pt$op %in% c(":=", "==", "<", ">"))
  in_expr <- .fmi_rows_named(pt, unique(unlist(lapply(expr_rows, function(k) {
    c(.fmi_expression_names(pt$lhs[k]), .fmi_expression_names(pt$rhs[k]))
  }))))
  ov_x <- tryCatch(lavaan::lavNames(fit, "ov.x"), error = function(e) {
    character(0)
  })
  free_alone <- function(j) {
    pt$free[j] > 0L && sum(comp == comp[j]) == 1L && !j %in% in_expr
  }
  not_reparameterization <- list(
    confirmed = NA,
    why = paste0("centering would also change an intercept that the model ",
                 "fixes or constrains (a latent outcome's, or one held ",
                 "equal to another), so the centered model would be ",
                 "another model"))
  slope_rows <- which(pt$op == "~" & pt$rhs %in% vars)
  for (i in slope_rows) {
    j <- which(pt$op == "~1" & pt$lhs == pt$lhs[i] & grp == grp[i])
    if (length(j) != 1L || !free_alone(j)) return(not_reparameterization)
  }
  mean_rows <- which(pt$op == "~1" & pt$lhs %in% vars)
  for (m in mean_rows) {
    fixed_x <- pt$free[m] == 0L && isTRUE(opt$fixed.x) &&
      pt$lhs[m] %in% ov_x
    if (!fixed_x && !free_alone(m)) return(not_reparameterization)
  }

  centre <- vapply(vars, function(v) {
    mean(unlist(lapply(dat, function(x) x[, v])), na.rm = TRUE)
  }, numeric(1L))
  spt <- pt
  spt$ustart <- spt$est
  spt$start <- NULL
  spt$est <- NULL
  spt$se <- NULL
  for (i in slope_rows) {
    j <- which(pt$op == "~1" & pt$lhs == pt$lhs[i] & grp == grp[i])
    spt$ustart[j] <- spt$ustart[j] + pt$est[i] * centre[[pt$rhs[i]]]
  }
  for (m in mean_rows) {
    spt$ustart[m] <- spt$ustart[m] - centre[[pt$lhs[m]]]
  }
  cdat <- lapply(dat, function(x) {
    x <- as.data.frame(x)
    for (v in vars) x[[v]] <- x[[v]] - centre[[v]]
    x
  })
  args <- list(model = spt, estimator = "ML", se = "standard",
               meanstructure = TRUE, fixed.x = opt$fixed.x,
               conditional.x = FALSE, test = "none", baseline = FALSE,
               missing = opt$missing, h1 = TRUE)
  if (length(dat) == 1L) {
    args$data <- cdat[[1L]]
  } else {
    group_var <- ".dmar_fmi_group."
    args$data <- do.call(rbind, lapply(seq_along(cdat), function(g) {
      x <- cdat[[g]]
      x[[group_var]] <- group_labels[g]
      x
    }))
    args$group <- group_var
    args$group.label <- group_labels
  }
  for (nm in c("em.h1.args", "h1.missing.method")) {
    if (!is.null(opt[[nm]])) args[[nm]] <- opt[[nm]]
  }
  fit_c <- tryCatch(suppressMessages(suppressWarnings(
    do.call(lavaan::lavaan, args))), error = function(e) NULL)
  if (is.null(fit_c) || !isTRUE(lavaan::lavInspect(fit_c, "converged"))) {
    return(list(confirmed = NA,
                why = "the model fit to the centered data did not converge"))
  }
  confirmed <- tryCatch({
    suppressWarnings(.fmi_core(fit_c, information, confirm_location = FALSE))
    TRUE
  }, dmar_fmi_not_inverted = function(e) FALSE, error = function(e) NA)
  list(confirmed = confirmed,
       why = if (is.na(confirmed)) {
         "the computation for the model fit to the centered data failed"
       })
}

# .fmi_check_maximum(fit_observed, pt_o, mono, group_labels = NULL)
#
# Whether the fit's estimates may not be at a maximum of the likelihood.
# Every FMI is computed at the fit's estimates, and the refits hold them
# whatever they are, so a fit that stopped short of the maximum, or that
# is not the maximum likelihood fit at all, gets FMIs that are not those
# of its maximum likelihood estimates, with nothing in either refit to
# say so. lavaan's own check of the gradient after it fits a model
# compares each element with an absolute 1e-3 in the variables' own
# metric, where the gradient of a parameter of variables measured in
# large units is small whether or not the fit is at the maximum. The
# check here is made in the metric of the refits, where every variable's
# standard deviation is near 1 (the caller makes it only there): one
# Newton step from the estimates, step = -N V g, with g lavaan's gradient
# of its FIML objective (the minus log likelihood over N) and V the
# refit's covariance matrix of the estimates, whose inverse is the
# information times N projected onto the equality constraints, so a
# constrained fit at its constrained maximum takes no step. The screen
# trips when a step moves an estimate by more than 1e-4 times
# max(1, |estimate|): the fits at the maximum studied for it (one and two
# groups, equality constraints, MLR, growth, small samples) took steps
# below 6e-6. It is a screen and not the verdict, because a fit can stop
# that far short of the maximum with no consequence that matters for the
# FMIs: in small samples with most of the outcome missing (60 cases, 16 of
# them observing it, FMIs of .97 to .99) lavaan's optimizer stopped about
# 5e-4 standard errors short of the maximum (a log likelihood 5e-7 below
# it, which rescale_data, a tighter rel.tol, and a restart from the
# estimates did not move), with FMIs within 3e-5 of the maximum's; the
# caller computes the FMIs one step on and judges the change
# (.fmi_warn_maximum()). The
# step would mislead where an inequality constraint or a bound may be
# active, since a fit at a constrained maximum has a gradient there that
# V does not project away, so the check is not made for a model with an
# inequality constraint or a free parameter at its bound.
# Returns NULL when the screen passes, otherwise list(rows, step, rel,
# gain, term, est, mono): the free rows of pt_o, the step of each (refit
# metric), the relative step, the predicted gain in the log likelihood,
# and the term, estimate, and factor of the row with the largest relative
# step.
.fmi_check_maximum <- function(fit_observed, pt_o, mono, group_labels = NULL) {
  rows <- which(pt_o$free > 0L)
  if (!length(rows) || any(pt_o$op %in% c("<", ">"))) return(NULL)
  for (b in intersect(c("lower", "upper"), names(pt_o))) {
    bound <- pt_o[[b]][rows]
    at <- is.finite(bound) &
      abs(pt_o$est[rows] - bound) <= 1e-8 * (1 + abs(bound))
    if (any(at)) return(NULL)
  }
  g <- tryCatch(as.numeric(lavaan::lavInspect(fit_observed, "gradient")),
                error = function(e) NULL)
  V <- tryCatch(lavaan::lavInspect(fit_observed, "vcov"),
                error = function(e) NULL)
  if (is.null(g) || is.null(V) || length(g) != length(rows) ||
      !identical(dim(V), c(length(rows), length(rows))) ||
      !all(is.finite(g)) || !all(is.finite(V))) {
    return(NULL)
  }
  N <- sum(as.numeric(unlist(lavaan::lavInspect(fit_observed, "nobs"))))
  step <- -N * as.numeric(V %*% g)
  est <- pt_o$est[rows]
  rel <- abs(step) / pmax(1, abs(est))
  if (!any(rel > 1e-4)) return(NULL)
  j <- which.max(rel)
  i <- rows[j]
  list(rows = rows, step = step, rel = rel,
       gain = -0.5 * N * sum(g * step),
       term = .fmi_group_term(pt_o, i, group_labels),
       est = est[j] * mono[i], est_step = (est[j] + step[j]) * mono[i])
}

# .fmi_group_term(pt, i, group_labels = NULL)
#
# lavaan's "lhs op rhs" for row i of pt, with its group for a multiple
# group fit.
.fmi_group_term <- function(pt, i, group_labels = NULL) {
  term <- trimws(paste(pt$lhs[i], pt$op[i], pt$rhs[i]))
  if (pt$op[i] == ":=") term <- pt$lhs[i]
  if (!is.null(group_labels) && !is.null(pt$group) && pt$group[i] >= 1L) {
    term <- paste0(term, " (group ", group_labels[pt$group[i]], ")")
  }
  term
}

# .fmi_warn_maximum(screen, fit, pt_o, fmi_here, fmi_step, em_cov, cj,
#                   group_labels = NULL)
#
# The warning for a fit whose estimates are not at the maximum, given the
# screen of .fmi_check_maximum() and the FMIs at the fit's estimates
# (fmi_here) and one Newton step on (fmi_step, NULL when that computation
# failed), over the rows of the parameter table. The warning asserts that
# the FMIs are not those of the maximum likelihood estimates, so it is
# given when they differ from those one step on by more than 1e-4, the
# size the calibration of the check separates: a regression with every
# standard deviation near .01, which lavaan does not rescale, left FMIs
# 3e-5 from the maximum's and passes, while near .003 they were 7e-4
# off, a factor model of variables measured in units 100 times too small
# .011, and a fit lavaan rescaled with ov_order = "data" .22, all of which
# warn; a small sample with most of the outcome missing gave 3e-5 in one
# data set (passes) and 6e-4 in another (warns). When the change cannot be
# measured the screen alone decides. The words name a cause only where it
# applies: units far from 1 in a fit lavaan did not rescale (a raw
# standard deviation of .01 or less or 100 or more, where lavaan's
# tolerances, set in the variables' own metric, do not suit them and
# lavaan does not rescale them unless one reaches .001 or 1000), with
# rescale_data = TRUE as the remedy; and the variables in an order other
# than the model's (lavaan's ov_order = "data") in a fit lavaan did
# rescale, which lavaan 0.7-3 does not fit to the maximum. A fit made with
# rescale_data = TRUE, or one lavaan rescaled itself, was made in a metric
# in which the standard deviations are near 1, so units are not named for
# it and that remedy is not prescribed again. Otherwise the likelihood is
# flat in some direction and lavaan's optimizer stopped within its
# tolerance, which none of those remedies moved (the 1.1.0 review's small
# sample, multiplied by 20 and fit with rescale_data = TRUE, had been told
# to fit it again with rescale_data = TRUE).
.fmi_warn_maximum <- function(screen, fit, pt_o, fmi_here, fmi_step, em_cov,
                              cj, group_labels = NULL) {
  keep <- which(pt_o$free > 0L | pt_o$op == ":=")
  change <- NULL
  if (!is.null(fmi_step)) {
    d <- abs(fmi_step[keep] - fmi_here[keep])
    d[!is.finite(d)] <- NA_real_
    if (all(is.na(d)) || max(d, na.rm = TRUE) <= 1e-4) return(invisible())
    k <- keep[which.max(d)]
    change <- list(term = .fmi_group_term(pt_o, k, group_labels),
                   from = fmi_here[k], to = fmi_step[k])
  }
  covs <- if (is.list(em_cov)) em_cov else list(em_cov)
  raw_sd <- unlist(lapply(covs, function(S) {
    v <- sqrt(pmax(diag(as.matrix(S)), 0))
    f <- if (is.null(cj)) rep(1, length(v)) else cj[rownames(as.matrix(S))]
    f[is.na(f)] <- 1
    v * f
  }))
  raw_sd <- raw_sd[is.finite(raw_sd) & raw_sd > 0]
  # Whether lavaan fit the model in its internally rescaled metric
  # (rescale_data = TRUE, or "auto" for variables it judges badly
  # scaled, when it records its factors in the options as .rescale.cj).
  # Units are then not the cause, since the fit was made in a metric in
  # which the variables' standard deviations are near 1, and rescale_data
  # = TRUE is the remedy already taken; the order of the variables is the
  # cause only then.
  opt <- lavaan::lavInspect(fit, "options")
  lav_cj <- suppressWarnings(as.numeric(opt$.rescale.cj))
  rescaled <- isTRUE(opt$rescale.data) ||
    (length(lav_cj) > 0L && any(is.finite(lav_cj) & lav_cj != 1))
  # Units far from 1 are named only for a standard deviation of .01 or
  # less or 100 or more. Fit without rescaling, with every variable
  # multiplied by one constant, the factor model of the nine Holzinger and
  # Swineford tests (marker indicators, two of them partly missing) and a
  # regression of 200 cases on three partly missing predictors reached the
  # maximum (no FMI moved by more than 1e-4) with standard deviations from
  # .01 to 80 (the regression up to 200 times its own), and left it where
  # the factor model's largest standard deviation was 111 and the
  # regression's were all near .003 to .005; the small sample with most of
  # the outcome missing above, multiplied by 20 and fit with rescale_data =
  # TRUE, stopped short on its flat likelihood as it did in its own
  # units.
  units <- !rescaled && length(raw_sd) &&
    any(raw_sd <= 0.01 | raw_sd >= 100)
  ov_fit <- tryCatch(lavaan::lavNames(fit, "ov"), error = function(e) NULL)
  ov_model <- tryCatch(lavaan::lavNames(lavaan::parTable(fit), "ov"),
                       error = function(e) NULL)
  data_order <- rescaled && !is.null(ov_fit) && !is.null(ov_model) &&
    !identical(as.character(ov_fit), as.character(ov_model))
  fmt <- function(v) format(v, digits = 4L)
  cause <- if (units || data_order) {
    paste0(
      if (units) {
        paste0(" lavaan can stop short of the maximum for variables ",
               "measured in small or large units that it does not rescale ",
               "(here the standard deviations run from ",
               format(min(raw_sd), digits = 3L), " to ",
               format(max(raw_sd), digits = 3L), "; lavaan rescales a fit ",
               "only when a standard deviation is .001 or less or 1000 or ",
               "more, or two differ by a factor of 1000).")
      },
      if (data_order) {
        paste0(" lavaan 0.7-3 does not reach it when it rescales a fit made ",
               "with ov_order = \"data\" whose data columns are not in the ",
               "model's order.")
      },
      " Fit the model again with rescale_data = TRUE",
      if (data_order) " and the columns of the data in the model's order",
      ", and compute the fraction for that fit.")
  } else {
    paste0(" lavaan's optimizer stopped within its convergence tolerance ",
           "on a likelihood that is flat in some direction, as it is where ",
           "few values of a variable are observed, and refitting with ",
           "rescale_data = TRUE, a tighter rel.tol, or from these estimates ",
           "did not move it in the examples studied; the change above is ",
           "the error this leaves in the fractions.")
  }
  msg <- paste0("The fit's estimates are not at the maximum of the likelihood: ",
          "one Newton step from them moves the estimate of ", screen$term,
          " from ", fmt(screen$est), " to ", fmt(screen$est_step),
          ", raises the log likelihood by ", fmt(screen$gain),
          if (!is.null(change)) {
            paste0(", and changes the fraction of missing information of ",
                   change$term, " from ", fmt(change$from), " to ",
                   fmt(change$to))
          },
          ". The fractions of missing information are computed at the ",
          "fit's estimates, so they are not those of its maximum likelihood ",
          "estimates either.", cause)
  # Carry the same diagnostic to callers that build a fit internally,
  # where its parameter names and available remedies differ.
  i <- screen$rows[which.max(screen$rel)]
  .dmar_warn(
    msg,
    "fmi_not_at_maximum", lhs = pt_o$lhs[i], op = pt_o$op[i],
    rhs = pt_o$rhs[i], from = screen$est, to = screen$est_step,
    gain = screen$gain, change = change, units = isTRUE(units),
    data_order = isTRUE(data_order))
  invisible()
}

# ---------------------------------------------------------------------------
# (D) The sample size behind each row
# ---------------------------------------------------------------------------

# .fmi_equality_classes(pt)
#
# Connected components of the parameter rows under the equality
# constraints of the parameter table: rows that share a user label, rows
# that share a free parameter index, and rows whose labels appear together
# in an "==" constraint, a defined (:=) name in the constraint standing
# for the labels its definition reaches (.fmi_rows_named()), so that
# ind1 == ind2 with ind1 := a1*b1 and ind2 := a2*b2 links the same rows
# as a1*b1 == a2*b2. Returns an integer component id per row of `pt`
# (rows that are not parameters, such as := and == rows, keep their own
# id). Components are found by propagating the smallest row index along
# the links until nothing changes, which needs no recursion and no
# assignment outside the function.
.fmi_equality_classes <- function(pt) {
  n <- length(pt$lhs)
  comp <- seq_len(n)
  is_par <- !pt$op %in% c(":=", "==", "<", ">")
  label <- if (is.null(pt$label)) rep("", n) else pt$label
  plabel <- if (is.null(pt$plabel)) rep("", n) else pt$plabel
  edges <- list()

  link_runs <- function(key, keep) {
    idx <- which(keep)
    if (length(idx) < 2L) return(NULL)
    runs <- split(idx, key[idx])
    runs <- runs[lengths(runs) > 1L]
    do.call(rbind, lapply(runs, function(r) cbind(r[-length(r)], r[-1L])))
  }
  edges[[1L]] <- link_runs(label, is_par & nzchar(label))
  edges[[2L]] <- link_runs(as.character(pt$free), is_par & pt$free > 0L)

  eq_rows <- which(pt$op == "==")
  for (k in eq_rows) {
    nm <- unique(c(.fmi_expression_names(pt$lhs[k]),
                   .fmi_expression_names(pt$rhs[k])))
    rows <- .fmi_rows_named(pt, nm)
    if (length(rows) > 1L) {
      edges[[length(edges) + 1L]] <-
        cbind(rows[-length(rows)], rows[-1L])
    }
  }
  edges <- do.call(rbind, edges)
  if (is.null(edges) || nrow(edges) == 0L) return(comp)
  repeat {
    m <- pmin(comp[edges[, 1L]], comp[edges[, 2L]])
    changed <- comp[edges[, 1L]] != m | comp[edges[, 2L]] != m
    if (!any(changed)) break
    for (e in which(changed)) {
      comp[edges[e, 1L]] <- min(comp[edges[e, 1L]], m[e])
      comp[edges[e, 2L]] <- min(comp[edges[e, 2L]], m[e])
    }
  }
  comp
}

# The names an expression in a := or == row refers to (labels, plabels,
# or other defined names). An expression that does not parse refers to
# nothing we can resolve.
.fmi_expression_names <- function(expr) {
  out <- tryCatch(all.vars(parse(text = expr)), error = function(e) NULL)
  if (is.null(out)) character(0) else out
}

# .fmi_rows_named(pt, names)
#
# The parameter rows that `names` refer to: the rows whose label or plabel
# is one of them and, for a defined (:=) name, the rows its definition
# reaches, resolving definitions nested inside other definitions breadth
# first with a visited set (no recursion). A defined name is not the label
# of any parameter row, so without this a constraint or definition
# written through one reaches no row.
.fmi_rows_named <- function(pt, names) {
  n <- length(pt$lhs)
  is_par <- !pt$op %in% c(":=", "==", "<", ">")
  label <- if (is.null(pt$label)) rep("", n) else pt$label
  plabel <- if (is.null(pt$plabel)) rep("", n) else pt$plabel
  defined <- which(pt$op == ":=")
  queue <- as.character(names)
  seen <- character(0)
  rows <- integer(0)
  while (length(queue)) {
    nm <- queue[1L]
    queue <- queue[-1L]
    if (nm %in% seen) next
    seen <- c(seen, nm)
    d <- defined[pt$lhs[defined] == nm]
    if (length(d)) {
      queue <- c(queue, .fmi_expression_names(pt$rhs[d[1L]]))
    } else {
      rows <- c(rows, which(is_par & ((nzchar(label) & label == nm) |
                                        (nzchar(plabel) & plabel == nm))))
    }
  }
  unique(rows)
}

# .fmi_row_N(pt, nobs)
#
# The N basis of each row: the cases whose data inform the parameter. A
# group specific parameter gets its group's N; a parameter constrained
# equal across groups gets the summed N of the groups in its equality
# class; a defined (:=) parameter gets the summed N of the groups its
# expression reaches, resolving defined names inside other definitions.
# With one group every row gets that group's N.
.fmi_row_N <- function(pt, nobs) {
  n <- length(pt$lhs)
  if (length(nobs) == 1L) return(rep(as.numeric(nobs), n))
  comp <- .fmi_equality_classes(pt)
  group <- pt$group
  is_par <- !pt$op %in% c(":=", "==", "<", ">")
  groups_of_rows <- function(rows) {
    unique(group[is_par & comp %in% comp[rows]])
  }
  out <- numeric(n)
  for (i in seq_len(n)) {
    g <- if (is_par[i]) {
      groups_of_rows(i)
    } else if (pt$op[i] == ":=") {
      r <- .fmi_rows_named(pt, pt$lhs[i])
      if (length(r)) groups_of_rows(r) else seq_along(nobs)
    } else {
      seq_along(nobs)
    }
    g <- g[g >= 1L & g <= length(nobs)]
    out[i] <- sum(as.numeric(nobs[g]))
  }
  out
}

# .fmi_precision_terms(fmi, N)
#
# The fraction of missing information in planning terms (Savalei &
# Rhemtulla, 2012, p. 488): the effective sample size N(1 - FMI), the
# complete data N that would estimate the parameter as precisely; the
# target sample size N / (1 - FMI), the N with the same missingness that
# would match the precision of N complete cases; and the width inflation
# factor 1 / sqrt(1 - FMI), the ratio of the confidence interval width
# with the missing values to the width the same N would give complete.
# Their oracle: N = 494 and FMI = .31506 give 338 and 721.
.fmi_precision_terms <- function(fmi, N) {
  list(N_effective = N * (1 - fmi),
       N_target = N / (1 - fmi),
       width_inflation = 1 / sqrt(1 - fmi))
}

# A readable term for each row: "visual =~ x2", "x1 ~1", and for a defined
# parameter its name. `names`, when supplied, maps variable names in the
# lavaan model to the names to show (see .fmi_mlmr_names()).
.fmi_term <- function(pt, names = NULL) {
  lhs <- pt$lhs
  rhs <- pt$rhs
  if (length(names)) {
    hit <- lhs %in% names(names)
    lhs[hit] <- names[lhs[hit]]
    hit <- rhs %in% names(names)
    rhs[hit] <- names[rhs[hit]]
  }
  term <- ifelse(pt$op == "~1", paste(lhs, "~1"), paste(lhs, pt$op, rhs))
  defined <- pt$op == ":="
  term[defined] <- pt$lhs[defined]
  term
}


# ---------------------------------------------------------------------------
# (E) Small-sample degrees of freedom: df_missing() and df_barnard_rubin()
# ---------------------------------------------------------------------------

# Internal helpers for the missing data family: the small-sample degrees
# of freedom that turn a fraction of missing information into a t
# reference distribution. df_missing() and df_barnard_rubin() call them,
# and so will the Wald intervals that accept a missing information
# degrees of freedom rule. Nothing here is exported.
#
# Notation follows von Hippel (2016): gamma is the fraction of missing
# information of one parameter, df_complete the degrees of freedom the
# parameter's interval would use with complete data, N the number of
# cases, and k the number of parameters subtracted from the effective
# sample size.

# .df_complete(gamma, df_complete)
#
# von Hippel's (2016, Eq. 32) general rule, which is also the observed
# data degrees of freedom of Barnard and Rubin (1999; von Hippel's
# Eq. 24):
#   nu = df_complete * (1 - gamma) * (df_complete + 1) / (df_complete + 3).
# The ratio is written as (1 + 1/df_complete) / (1 + 3/df_complete) so that an
# infinite df_complete gives an infinite result (the normal reference) rather
# than Inf/Inf. At gamma = 1 the data carry no information about the
# parameter and the result is 0, also when df_complete is infinite (where the
# product would otherwise be Inf * 0). Vectorized; NA passes through.
.df_complete <- function(gamma, df_complete) {
  out <- df_complete * (1 - gamma) * (1 + 1 / df_complete) / (1 + 3 / df_complete)
  out[!is.na(gamma) & gamma >= 1] <- 0
  out
}

# .df_effective_n(gamma, N, k)
#
# von Hippel's (2016, Eq. 31) effective sample size rule,
#   nu* = N * (1 - gamma) - k,
# the complete data regression degrees of freedom N - k with N replaced by
# the effective sample size N(1 - gamma) of Savalei and Rhemtulla (2012).
# With gamma = 0 it is N - k exactly. It can be negative when N(1 - gamma)
# falls below k, which is what the floor is for. Vectorized; NA passes
# through.
.df_effective_n <- function(gamma, N, k) {
  N * (1 - gamma) - k
}

# .df_bound(df, floor)
#
# The lower bound of von Hippel's (2016) Eqs. 27 and 34: max(floor, df),
# elementwise, with NA kept as NA. Returns the bounded values and a
# logical marking where the bound changed the value.
.df_bound <- function(df, floor) {
  at_floor <- df < floor
  list(df = ifelse(at_floor, floor, df), at_floor = at_floor)
}

# .df_clean_fmi(fmi, tol = 1e-8)
#
# Validates a vector of fractions of missing information for use in a
# degrees of freedom rule. A value at or above 1 is refused (no
# information about the parameter, so no reference distribution). A value
# below -tol is refused: the structured and Hessian estimators can return
# clearly negative FMIs in small samples (Chen & Savalei, 2021), and a
# negative FMI has no degrees of freedom interpretation. A value in
# [-tol, 0) is rounding error around a complete data parameter and is set
# to 0. NA passes through.
.df_clean_fmi <- function(fmi, tol = 1e-8) {
  if (!is.numeric(fmi) && !(is.logical(fmi) && all(is.na(fmi)))) {
    .dmar_refuse_argument(c(
      "'fmi' must be numeric: the fraction of missing information of ",
      "each parameter, or the table fraction_missing_information() ",
      "returns."), arg = "fmi")
  }
  fmi <- as.numeric(fmi)
  if (!length(fmi)) {
    .dmar_refuse_argument("'fmi' has no values.", arg = "fmi")
  }
  if (any(is.infinite(fmi))) {
    .dmar_refuse_argument("'fmi' must be finite.", arg = "fmi")
  }
  ok <- !is.na(fmi)
  if (any(fmi[ok] >= 1)) {
    .dmar_refuse_argument(c(
      "A fraction of missing information must be below 1; at 1 the ",
      "data carry no information about the parameter, and no ",
      "reference distribution applies. Got ",
      paste(format(fmi[ok][fmi[ok] >= 1]), collapse = ", "), "."), arg = NULL)
  }
  if (any(fmi[ok] < -tol)) {
    .dmar_refuse_argument(c(
      "A negative fraction of missing information has no degrees of ",
      "freedom interpretation; such values come from the structured ",
      "or Hessian estimators in small samples. Use information = ",
      "\"unstructured\" in fraction_missing_information(), whose ",
      "estimate cannot be negative. Got ",
      paste(format(fmi[ok][fmi[ok] < -tol]), collapse = ", "), "."), arg = NULL)
  }
  fmi[ok & fmi < 0] <- 0
  fmi
}

# .df_recycle(x, n, name)
#
# Recycles an argument of length 1 to length n; any other length that is
# not n is an error naming the argument.
.df_recycle <- function(x, n, name) {
  if (length(x) == 1L) return(rep(x, n))
  if (length(x) == n) return(x)
  .dmar_refuse_argument(c(
    "'", name, "' must have length 1 or the length of the fractions ",
    "of missing information (", n, "); it has length ", length(x), "."),
    arg = NULL)
}

# .df_terms(x, n)
#
# Row labels for a degrees of freedom table: the names of x when it has
# them, otherwise parameter_1, parameter_2, and so on.
.df_terms <- function(x, n) {
  nm <- names(x)
  if (is.null(nm) || all(!nzchar(nm))) return(paste0("parameter_", seq_len(n)))
  nm[!nzchar(nm) | is.na(nm)] <- paste0("parameter_",
                                        which(!nzchar(nm) | is.na(nm)))
  nm
}

# .df_check_floor(floor)
#
# The floor is a single finite nonnegative number.
.df_check_floor <- function(floor) {
  if (!is.numeric(floor) || length(floor) != 1L || is.na(floor) ||
      !is.finite(floor) || floor < 0) {
    .dmar_refuse_argument(c(
      "'floor' must be a single finite number of 0 or more; 3 is ",
      "von Hippel's (2016) bound and 0 removes it."), arg = "floor")
  }
  invisible(floor)
}

# .df_tidy(x)
#
# The tidy() view shared by the degrees of freedom tables: one row per
# parameter, term, the bounded degrees of freedom as the estimate, and
# every other column passed through in table order (at_floor last).
.df_tidy <- function(x) {
  raw <- x
  class(raw) <- "data.frame"
  out <- data.frame(term = as.character(raw$term),
                    estimate = as.numeric(raw$df),
                    stringsAsFactors = FALSE)
  for (nm in setdiff(names(raw), c("term", "df"))) out[[nm]] <- raw[[nm]]
  out
}

# .df_require_finite(df, fn)
#
# Stops when a degrees of freedom rule left a parameter without degrees of
# freedom, which happens only when its fraction of missing information
# could not be computed (df_missing() passes an NA FMI through as NA). A
# Wald interval and p-value referred to NA degrees of freedom are NA, and
# a result that carried them would say nothing of why. mlmr() and
# mlmr_mv() call it on the degrees of freedom of wald_df =
# "missing_information"; `df` is named by the terms, and fn names the
# function for the message. .fmi_core() stops with the reason when a
# refit cannot be computed, so this is the backstop behind it.
.df_require_finite <- function(df, fn) {
  bad <- is.na(df)
  if (!any(bad)) return(invisible(df))
  terms <- names(df)
  if (is.null(terms)) terms <- paste0("parameter_", seq_along(df))
  .dmar_refuse(c(
    "wald_df = \"missing_information\" could not compute the fraction ",
    "of missing information of ", paste(terms[bad], collapse = ", "),
    " from the fit, so ", fn, "() has no degrees of freedom for ",
    if (sum(bad) == 1L) "its interval" else "their intervals",
    ". Leave wald_df at Inf for the normal reference, or pass a number ",
    "of degrees of freedom."), "computation_failed")
}


# ---------------------------------------------------------------------------
# (F) The t reference distribution of a Wald interval: wald_df
# ---------------------------------------------------------------------------

# reliability_omega(), reliability_alpha() (and reliability(), which
# dispatches to them), fa_confirmatory() (and cfa_1() and cfa_2(), which call it),
# mlmr(), and mlmr_mv() take wald_df = Inf: the degrees of freedom of the
# t distribution to which a Wald interval refers the ratio of an estimate
# to its standard error, and in fa_confirmatory() and the mlmr pair the p-value as
# well. The estimate and its standard error are unchanged; only the
# quantile changes, from the normal one to the t one (von Hippel, 2016,
# Eqs. 28 and 29). R evaluates qt() and pt() at infinite degrees of
# freedom by calling qnorm() and pnorm(), so the default reproduces every
# interval and p-value of DMAR 1.0.x to the last bit, which the tests
# check with identical().
#
# A number gives that df to every Wald interval of the call. The keyword
# "missing_information" computes one df per parameter by von Hippel's
# (2016, Eqs. 31 and 34) effective sample size rule, max(3, N(1 - FMI) -
# k), with each parameter's FMI taken from the fit by the machinery of
# fraction_missing_information() (the unstructured estimator, which cannot
# be negative) and k as .wald_df_rule() holds it for the function:
#
#   mlmr(), mlmr_mv(): k is the number of regression parameters of the
#     coefficient's outcome, the case von Hippel's own simulations
#     evaluated, so that complete data give lm()'s N - k.
#   fa_confirmatory(): k = 1, for the loadings only (its omega rows take
#     logit-scale limits, where the t reference overcovered).
#   reliability_omega(): k = 1, for omega's "ml" and "mlr" intervals.
#   reliability_alpha(): no rule; the keyword is refused.
#
# The k = 1 rule and its scope come from the calibration study in
# methodology/missing_data (its README.md has the design and the decision
# rule, committed before the run, and RESULTS.md the tables). It simulated
# one factor models fit by FIML, as fa_confirmatory() writes them (factor variance
# 1, free loadings, intercepts, and error variances; omega defined from
# the loadings and error variances): 6 items at N = 30, 50, 100, and 200
# and 12 items at N = 50 and 100, normal data fit by ML with normal theory
# standard errors and nonnormal data (standardized chi square on 4 df) by
# MLR with Huber-White standard errors, no missing values or 15 or 30
# percent of the values of half the items missing completely at random or
# at random, 2000 replications in each of 60 cells (and 4 cells of von
# Hippel's own design as a positive control). The z interval was
# compared with t intervals on four of von Hippel's rules (Eq. 31 with
# k = 1 and k = the number of items, Eq. 32 with df_complete = N - 1 and
# N - p) under a decision rule fixed before the run: a variant passes for
# a target when its mean |coverage - .95| over the cells with N <= 100
# and missing values is below z's, its coverage exceeds .96 in no cell,
# and it is within .005 of the best variant's, with a cell's coverage the
# mean over the target's parameters; the simplest passing variant wins.
# max(3, N(1 - FMI) - 1) passed for the loadings (.0171 against z's
# .0247, largest cell coverage .9502, a mean over the cell's loadings;
# single loadings reached .9595 with missing values, none above .96) and
# for omega on its own scale (.0180 against .0242, largest .9550), and it
# passed again within
# the normal theory cells (loadings .0055 against .0122, omega .0094
# against .0146) and within the robust ones (.0288 against .0373, .0267
# against .0337) taken alone (methodology/missing_data/scope.R). In the
# cells without missing values, where it is t on N - 1, it moved the
# coverage of the loadings with normal data from between .936 and .946 to
# between .945 and .949 (cell means; single loadings .934 to .963, where z
# reached .962).
# It failed on the logit scale (omega's coverage reached .9645 with normal
# data, where z's was within .0044 of .95 on average and t's .0047) and
# for alpha, so the _logistic intervals and reliability_alpha() keep the
# normal reference. For alpha the t reference met the first criterion
# barely (.0042 against .0043) and failed the second in one cell, 12 items
# at N = 50 with 30 percent MCAR (.966), the cell decision.txt flags
# because the study's fit of the saturated model, which unlike
# .fiml_cov() starts from lavaan's defaults and iterates, converged in 53
# percent of the replications there. Without that cell it passes (.0035
# against .0042, at most .9595; scope.R, section 5), but the verdict is
# the rule's over every cell, as fixed before the run. With
# nonnormal data both intervals still undercover (the robust standard
# errors are too small at these sample sizes); the t reference narrows the
# shortfall without removing it.
#
# What the keyword covers follows from that design, and everything else is
# refused with its reason (.wald_df_scope() before the fit,
# .wald_df_fit_scope() on the fit): the loading rows of fa_confirmatory()
# and omega's Wald intervals in reliability_omega() (the other rows of
# fa_confirmatory() keep the normal reference, since the study evaluated
# neither factor correlations, intercepts, error variances, nor the sums,
# AVE, and H, and its omega rows take logit-scale limits, where the t
# reference overcovered); maximum likelihood fits of continuous items
# (estimator "ML" or "MLR") with normal theory or Huber-White standard
# errors, by FIML or with complete data (listwise deletion or a
# covariance matrix, where every FMI is 0 and the df is N - 1). Under FIML
# the standard errors must come from the observed information computed
# from the Hessian, as the study's did; robust.sem is accepted with
# complete data only, since lavaan does not compute it under FIML but
# falls back to listwise deletion (with a warning), which would not be
# the FIML fit the call asked for. With
# complete data lavaan's normal theory standard errors use the expected
# information and reliability_omega(ci_method = "mlr") the robust.sem
# sandwich, where the study's cells without missing values (fit with
# missing = "ml") used the observed information and the Huber-White
# sandwich; methodology/missing_data/complete_data_se.R measures the
# difference (DMAR's are 1 to 2 percent shorter for the loadings and less
# than 1 percent for omega under normal theory, and the robust.sem ones of
# the loadings 4 to 7 percent shorter and of omega 1 to 2 percent, at
# N = 30 and 50), so that DMAR's complete data t intervals are, if
# anything, a little shorter than the ones the study evaluated, whose
# cell coverage there was at most .949 for the loadings and .953 for
# omega.
#
# The study fit one factor models with free loadings; several factors,
# equal loadings, and auxiliary variables are accepted, and the rule is
# applied to each loading as it stands, which the study supports
# by analogy rather than directly.

# .wald_crit(conf_level, wald_df = Inf)
#
# The two-sided critical value of a Wald interval at the confidence level,
# from the t distribution on wald_df degrees of freedom: qnorm() itself
# when wald_df is infinite. Vectorized over wald_df.
.wald_crit <- function(conf_level, wald_df = Inf) {
  stats::qt(1 - (1 - conf_level) / 2, wald_df)
}

# .wald_p(z, wald_df = Inf)
#
# The two-sided p-value of a Wald ratio referred to the t distribution on
# wald_df degrees of freedom: 2 * pnorm(-abs(z)) itself when wald_df is
# infinite. Vectorized over z and wald_df; NA passes through.
.wald_p <- function(z, wald_df = Inf) {
  2 * stats::pt(-abs(z), wald_df)
}

# .wald_df_is_normal(wald_df)
#
# Whether wald_df leaves the normal reference in place (the default, Inf).
.wald_df_is_normal <- function(wald_df) {
  is.numeric(wald_df) && length(wald_df) > 0L && all(is.infinite(wald_df))
}

# .wald_df_rule(fn)
#
# The rule wald_df = "missing_information" applies for the function named
# fn, or NULL where the keyword is not accepted. This is the one place the
# keyword is granted. A rule names a df_missing() method and supplies its
# argument as a function of the parameters (rows, a data frame with lhs,
# op, and rhs) and of what the calling function knows (context: N,
# n_coef, n_items, factors, terms, omega):
#
#   list(method = "effective_n", k = function(rows, context) ...)
#   list(method = "complete_df", df_complete = function(rows, context) ...)
#
# and may add `rows`, a function of the same arguments returning which
# parameters the rule covers (the others keep the normal reference; all
# are covered when it is absent), and `calibrated = TRUE`, which has
# .wald_df_resolve() check the fit against the calibration study's scope
# (.wald_df_fit_scope()).
#
# For mlmr() and mlmr_mv() it is von Hippel's (2016, Eq. 31) rule with k
# the number of regression parameters of the coefficient's outcome, slopes
# and intercept, which with complete data gives N - k, the degrees of
# freedom of lm(). For fa_confirmatory() and reliability_omega() it is the same rule
# with k = 1, for the loadings (op "=~") and the omega rows (the defined
# parameters named in context$omega), the variant the calibration study
# chose (see the head of this section). reliability_alpha() has no entry:
# no variant passed for alpha.
.wald_df_rule <- function(fn) {
  regression <- list(
    method = "effective_n",
    k = function(rows, context) rep(context$n_coef, nrow(rows)))
  calibrated <- list(
    method = "effective_n",
    k = function(rows, context) rep(1, nrow(rows)),
    rows = function(rows, context) {
      rows$op == "=~" | (rows$op == ":=" & rows$lhs %in% context$omega)
    },
    calibrated = TRUE)
  switch(fn,
         mlmr = regression,
         mlmr_mv = regression,
         fa_confirmatory = calibrated,
         reliability_omega = calibrated,
         NULL)
}

# .wald_df_check(wald_df, fn)
#
# Validates the wald_df argument of the function named fn: a single
# positive number (Inf, the default, is the normal reference), or the
# keyword "missing_information" where .wald_df_rule() has a rule for fn.
# Returns the number (as a double) or the keyword.
.wald_df_check <- function(wald_df, fn) {
  has_rule <- !is.null(.wald_df_rule(fn))
  if (identical(wald_df, "missing_information")) {
    if (has_rule) return(wald_df)
    .dmar_refuse_argument(.wald_df_no_rule(fn), arg = NULL)
  }
  if (!is.numeric(wald_df) || length(wald_df) != 1L || is.na(wald_df) ||
      wald_df <= 0) {
    .dmar_refuse_argument(c(
      "'wald_df' must be a single positive number of degrees of ",
      "freedom", if (has_rule) " or \"missing_information\"",
      "; Inf, the default, gives the normal reference."), arg = "wald_df")
  }
  as.numeric(wald_df)
}

# .wald_df_no_rule(fn)
#
# Why the keyword is refused for a function without a rule. For
# reliability_alpha() the calibration study's verdict, by the decision
# rule fixed before the run: with normal data alpha's FIML interval on the
# normal reference was already close to .95 (mean |coverage - .95| .0043
# over the cells with N <= 100 and missing values, coverage .940 to .958),
# the t reference came barely closer (.0042), and its coverage reached
# .966 in one cell (12 items, N = 50, and 30 percent missing completely at
# random, where the study's saturated fit converged in 53 percent of the
# replications; at most .9595 in every other cell, and without that cell
# the rule would have passed, .0035 against .0042). For fa_bifactor() and
# fa_second_order(), that the study's one factor models do not cover their
# loadings or their indices.
.wald_df_no_rule <- function(fn) {
  if (identical(fn, "reliability_alpha")) {
    return(paste0(
      "wald_df = \"missing_information\" is not available for ",
      "reliability_alpha(): in the calibration study behind the keyword no ",
      "rule passed for alpha. On max(3, N(1 - FMI) - 1) degrees of freedom ",
      "the t reference brought alpha's FIML interval barely closer to .95 ",
      "than the normal reference, which was already close with normal data ",
      "(mean |coverage - .95| .0042 against .0043), and in one cell (12 ",
      "items, N = 50, 30 percent missing completely at random, where the ",
      "study's fit of the saturated model converged in 53 percent of the ",
      "replications) its coverage reached .966, above the .96 the study's ",
      "decision rule allowed. Pass a number of degrees of freedom, or ",
      "leave wald_df at Inf."))
  }
  if (identical(fn, "fa_bifactor")) {
    return(paste0(
      "wald_df = \"missing_information\" is not available for ",
      "fa_bifactor(): the calibration study behind the keyword simulated ",
      "one factor models and evaluated their loadings and omega, so ",
      "nothing calibrates it for the loadings of a general factor and ",
      "group factors or for the bifactor indices. Pass a number of degrees ",
      "of freedom, or leave wald_df at Inf for the normal reference."))
  }
  if (identical(fn, "fa_second_order")) {
    return(paste0(
      "wald_df = \"missing_information\" is not available for ",
      "fa_second_order(): the calibration study behind the keyword ",
      "simulated one factor models and evaluated their loadings and omega, ",
      "so nothing calibrates it for the first-order and second-order ",
      "loadings of a second-order model or for its indices. Pass a number ",
      "of degrees of freedom, or leave wald_df at Inf for the normal ",
      "reference."))
  }
  paste0("wald_df = \"missing_information\" is not available for ", fn,
         "(): no rule for its degrees of freedom has been calibrated. ",
         "Pass a number of degrees of freedom, or leave wald_df at Inf ",
         "for the normal reference.")
}

# .wald_df_scope(wald_df, fn, spec)
#
# Stops when wald_df is the keyword and the call lies outside what the
# calibration study evaluated for fn, with the reason and what to do
# instead; returns invisibly otherwise. `spec` holds the call's settings:
# for fa_confirmatory(), output, se, estimator, missing, two_stage, ordered,
# complete (listwise deletion or a covariance matrix, so that no value is
# missing in the fit), and test (the test statistic passed on in ..., or
# NULL); for reliability_omega(), denominator, ci_method, and missing. fa_confirmatory() and reliability_omega() call it before anything is
# fitted; .wald_df_fit_scope() checks the fit itself.
.wald_df_scope <- function(wald_df, fn, spec) {
  if (!identical(wald_df, "missing_information")) return(invisible(TRUE))
  reason <- switch(fn,
                   fa_confirmatory = .wald_df_scope_cfa_k(spec),
                   reliability_omega = .wald_df_scope_omega(spec),
                   NULL)
  if (is.null(reason)) return(invisible(TRUE))
  .dmar_refuse_argument(c(
    "wald_df = \"missing_information\" ", reason), arg = NULL)
}

# The sentence ends shared by the refusals below.
.WALD_DF_FIML_ONLY <- paste0(
  "the rule takes each parameter's fraction of missing information from ",
  "a full information maximum likelihood fit, which DMAR does not compute ",
  "for %s, and the calibration study fit by full information maximum ",
  "likelihood.")
.WALD_DF_SE_SCOPE <- paste0(
  "the calibration study used normal theory standard errors (se = ",
  "\"standard\") and the Huber-White sandwich (se = ",
  "\"robust.huber.white\", MLR's under full information maximum ",
  "likelihood)")
.WALD_DF_LISTWISE_REMEDY <- paste0(
  "Use estimator = \"MLR\" with missing = \"fiml\" (its scaled test ",
  "statistic and Huber-White standard errors are computed under full ",
  "information maximum likelihood, and lie within the calibration), or ",
  "the test with missing = \"listwise\", or pass a number of degrees of ",
  "freedom.")

# The test statistics for which lavaan 0.7-3 sets missing = "ml" to
# "listwise" (lav_options_set()), used only to refuse the keyword before
# the fit; the fit itself is checked afterwards (.wald_df_fit_scope()),
# which catches any other spelling or route.
.LAVAAN_LISTWISE_TESTS <- c("satorra.bentler", "mean.var.adjusted",
                            "scaled.shifted", "mean.var.adjusted.corrected",
                            "scaled.shifted.corrected")

# .wald_df_scope_omega(spec)
#
# reliability_omega(): the Wald intervals of the model implied coefficient
# on the scale of omega ("ml", "mlr"), with complete data or by FIML.
.wald_df_scope_omega <- function(spec) {
  m <- spec$ci_method
  if (identical(spec$denominator, "observed")) {
    return(paste0(
      "applies to the Wald interval of the model implied omega ",
      "(denominator = \"model_implied\" with ci_method = \"ml\" or ",
      "\"mlr\"), and robust omega (denominator = \"observed\") has a ",
      "bootstrap interval only. Leave wald_df at Inf."))
  }
  if (m %in% c("ml_logistic", "mlr_logistic")) {
    return(paste0(
      "is not available with ci_method = \"", m, "\": in the calibration ",
      "study the t reference on the logit scale covered up to .9645 with ",
      "normal data, above the .96 the study allowed, where the normal ",
      "reference was already within .0044 of .95 on average, so the rule ",
      "is given for the intervals on the scale of omega only. Use ",
      "ci_method = \"", sub("_logistic$", "", m), "\", or pass a number of ",
      "degrees of freedom."))
  }
  if (m %in% c("adf", "adf_logistic")) {
    return(paste0(
      "is not available with ci_method = \"", m, "\": the calibration ",
      "study evaluated the maximum likelihood standard errors, normal ",
      "theory and robust, not the asymptotically distribution free one. ",
      "Use ci_method = \"ml\" or \"mlr\", or pass a number of degrees of ",
      "freedom."))
  }
  if (!m %in% c("ml", "mlr")) {
    return(paste0(
      "sets the reference distribution of a Wald interval, and ",
      "ci_method = \"", m, "\" is not one; the keyword applies to ",
      "ci_method = \"ml\" and \"mlr\" with denominator = ",
      "\"model_implied\". Leave wald_df at Inf."))
  }
  if (identical(spec$missing, "two_stage")) {
    return(paste0(
      "is not available with missing = \"two_stage\": ",
      sprintf(.WALD_DF_FIML_ONLY, "a two-stage fit"), " Use missing = ",
      "\"fiml\", or pass a number of degrees of freedom."))
  }
  NULL
}

# .wald_df_scope_cfa_k(spec)
#
# fa_confirmatory(): the table of output "verbose" from a maximum
# likelihood fit of continuous items, by FIML or with complete data, with
# normal theory or Huber-White standard errors, and robust.sem only under
# listwise deletion or from a covariance matrix whose asymptotic
# covariance matrix of the moments is passed as lavaan's NACOV through ...
# (lavaan computes robust.sem from raw data or a supplied NACOV only,
# refusing S alone itself, and falls back to listwise deletion when FIML
# is asked for; its "robust" is the Huber-White sandwich under FIML and
# robust.sem with complete data, so it is accepted either way). The scope
# does not refuse S with robust.sem: lavaan does, unless NACOV is given,
# and with it the intervals equal those of the raw data listwise fit.
.wald_df_scope_cfa_k <- function(spec) {
  # fa_confirmatory() hands 'se' to lavaan unchecked; anything but one string is
  # lavaan's to refuse, and here it is simply outside the calibration.
  if (!is.character(spec$se) || length(spec$se) != 1L) spec$se <- "?"
  if (spec$output %in% c("fit", "summary")) {
    return(paste0(
      "sets the reference distribution of the Wald intervals in fa_confirmatory()'s ",
      "tables, and output = \"", spec$output, "\" returns lavaan's own ",
      "object, whose intervals keep lavaan's normal reference; the keyword ",
      "applies to output = \"verbose\". Leave wald_df at Inf."))
  }
  if (identical(spec$output, "measurement")) {
    return(paste0(
      "is not available with output = \"measurement\": the keyword applies ",
      "to the loadings, which that table does not show; its omega, AVE, and ",
      "H rows take limits formed on the logit scale on the normal ",
      "reference, since the calibration study found the t reference ",
      "overcovered there. Use output = \"verbose\", or pass a number of ",
      "degrees of freedom."))
  }
  if (identical(spec$output, "standardized")) {
    return(paste0(
      "is not available with output = \"standardized\": the calibration ",
      "study evaluated the unstandardized loadings, with each factor's ",
      "variance fixed at 1, and omega, not the standardized solution. Use ",
      "output = \"verbose\", or pass a number of degrees of freedom."))
  }
  if (spec$se %in% c("none", "bootstrap")) {
    return(paste0(
      "sets the reference distribution of a Wald interval, and ",
      if (spec$se == "none") {
        "se = \"none\" gives no standard errors, so the tables have none."
      } else {
        paste0("with se = \"bootstrap\" lavaan reports percentile ",
               "bootstrap intervals, not Wald intervals.")
      }, " Leave wald_df at Inf."))
  }
  if (isTRUE(spec$ordered) || is.character(spec$ordered)) {
    return(paste0(
      "is not available with ordered items: they are fit by WLSMV, which ",
      "has no fraction of missing information of the kind the rule takes, ",
      "and the calibration study fit continuous items by maximum ",
      "likelihood. Pass a number of degrees of freedom."))
  }
  if (!spec$estimator %in% c("ML", "MLR")) {
    return(paste0(
      "is not available with estimator = \"", spec$estimator, "\": the ",
      "rule takes each parameter's fraction of missing information from a ",
      "maximum likelihood fit, and the calibration study fit by maximum ",
      "likelihood (estimator = \"ML\" and \"MLR\"). Use estimator = ",
      "\"ML\" or \"MLR\", or pass a number of degrees of freedom."))
  }
  if (isTRUE(spec$two_stage) || identical(spec$missing, "pairwise")) {
    return(paste0(
      "is not available with missing = \"", spec$missing, "\": ",
      sprintf(.WALD_DF_FIML_ONLY, if (isTRUE(spec$two_stage)) {
        "a two-stage fit"
      } else {
        "a fit to pairwise present data"
      }), " Use missing = \"fiml\" (or \"listwise\"), or pass a number of ",
      "degrees of freedom."))
  }
  if (spec$se == "robust.sem" && !isTRUE(spec$complete)) {
    # lavaan has no robust.sem sandwich under FIML: lav_options_set()
    # switches missing to "listwise", with a warning, and fits the
    # complete cases. The keyword's degrees of freedom would then describe
    # a fit other than the FIML fit the call asks for, so the combination
    # is refused whether or not a value is missing; with complete data
    # missing = "listwise" makes the same fit and is accepted.
    return(paste0(
      "is not available with se = \"robust.sem\" and missing = \"",
      spec$missing, "\": lavaan does not compute the robust.sem sandwich ",
      "under full information maximum likelihood but falls back to ",
      "listwise deletion, so the fit would not be the FIML fit whose ",
      "fractions of missing information the rule takes. Use se = ",
      "\"robust.huber.white\" (MLR's sandwich under full information ",
      "maximum likelihood) with missing = \"fiml\", or se = \"robust.sem\" ",
      "with missing = \"listwise\", or pass a number of degrees of ",
      "freedom."))
  }
  sb_test <- intersect(tolower(as.character(spec$test)),
                       .LAVAAN_LISTWISE_TESTS)
  if (length(sb_test) && !isTRUE(spec$complete)) {
    # The same fallback through a test statistic passed on in fa_confirmatory()'s
    # ...: lav_options_set() sets missing to "listwise" for any test of
    # the Satorra-Bentler kind. Refused before the fit, whether or not a
    # value is missing, as robust.sem is; .wald_df_fit_scope() checks the
    # fit itself for any other route to the same fallback.
    return(paste0(
      "is not available with test = \"", sb_test[1L], "\" and missing = \"",
      spec$missing, "\": lavaan does not compute that test under full ",
      "information maximum likelihood but falls back to listwise deletion, ",
      "so the fit would not be the FIML fit whose fractions of missing ",
      "information the rule takes. ", .WALD_DF_LISTWISE_REMEDY))
  }
  if (!spec$se %in% c("standard", "robust.huber.white", "robust",
                      "robust.sem")) {
    return(paste0(
      "is not available with se = \"", spec$se, "\": ", .WALD_DF_SE_SCOPE,
      ". Use one of those, or pass a number of degrees of freedom."))
  }
  NULL
}

# .wald_df_fit_complete(fit)
#
# Whether no value is missing among the cases the fit used (complete data,
# a listwise deletion fit, or a fit to a covariance matrix, which carries
# no raw data). lavTech() returns the data as a list with one matrix per
# group, or NULL in place of a group's matrix for a fit to moments, where
# lavInspect(fit, "data") stops with an error.
.wald_df_fit_complete <- function(fit) {
  dat <- lavaan::lavTech(fit, "data")
  if (!is.list(dat)) dat <- list(dat)
  !any(vapply(dat, anyNA, logical(1L)))
}

# .wald_df_fit_fiml(fit, requested_missing, calibrated = FALSE)
#
# requested_missing is the caller's missing data option when it asked for
# full information maximum likelihood with raw data ("fiml" or "ml"), or
# NULL. lavaan replaces such a request with listwise deletion, with a
# warning, for the robust.sem sandwich and for any test statistic of the
# Satorra-Bentler kind; the fit is then complete, so every FMI would be 0
# and every df that of the complete cases (N - 1 for a loading, N - k for
# a regression coefficient), a fit other than the one the call asked for.
# The keyword is refused, for every function that takes it, whenever the
# fit's own missing option is not FIML although the call asked for it,
# whatever the route (the rule fa_confirmatory() also applies to robust.sem and to
# such a test before the fit), and read from the fit, so a route lavaan
# adds is caught as well. A covariance matrix, which lavaan always fits as
# complete data, is not such a request. `calibrated` chooses the remedy:
# for fa_confirmatory() and reliability_omega() the fits within the calibration.
#
# Where lavaan's listwise fit used every case read and no value is
# missing among them, it is the FIML fit: with complete data the FIML
# likelihood is the complete data likelihood, every FMI is 0 by
# definition, and the rule gives the complete data degrees of freedom
# (N - k, those of lm(), for mlmr() and mlmr_mv()). The uncalibrated rules
# accept such a fit. (A case with every variable missing, which FIML sets
# aside as well, also counts as a case lavaan did not use, so such data
# are refused.) The calibrated rule of fa_confirmatory() refuses the combination
# whether or not a value is missing, as its check before the fit does
# (.wald_df_scope_cfa_k()), and then says that no case was set aside
# rather than that the complete cases replaced the FIML fit.
.wald_df_fit_fiml <- function(fit, requested_missing, calibrated = FALSE) {
  if (is.null(requested_missing)) return(invisible(TRUE))
  opt <- lavaan::lavInspect(fit, "options")
  if (as.character(opt$missing) %in% c("ml", "ml.x")) return(invisible(TRUE))
  n_used <- sum(as.numeric(unlist(lavaan::lavInspect(fit, "nobs"))))
  n_read <- sum(as.numeric(unlist(lavaan::lavInspect(fit, "norig"))))
  every_case <- isTRUE(n_used == n_read) && .wald_df_fit_complete(fit)
  if (every_case && !calibrated) return(invisible(TRUE))
  test <- intersect(as.character(opt$test), .LAVAAN_LISTWISE_TESTS)
  route <- if (length(test)) {
    paste0("test = \"", test[1L], "\"")
  } else {
    paste0("se = \"", as.character(opt$se), "\"")
  }
  .dmar_refuse_argument(c(
    "wald_df = \"missing_information\" is not available for this ",
    "fit: the call asked for full information maximum likelihood ",
    "(missing = \"", requested_missing, "\"), ",
    if (every_case) {
      paste0("which lavaan does not compute with ", route, " but ",
             "replaces with listwise deletion. No value is missing ",
             "here, so that fit used every one of the ", format(n_read),
             " cases read, but the keyword is refused for the ",
             "combination whether or not a value is missing, as it is ",
             "before the fit for se = \"robust.sem\" and the test ",
             "statistics of the Satorra-Bentler kind. ")
    } else {
      paste0("and lavaan fit the complete cases instead (",
             format(n_used), " of the ", format(n_read), " cases ",
             "read), as it does with ", route, ", so the fit is not ",
             "the FIML fit whose fractions of missing information the ",
             "rule takes. ")
    },
    if (calibrated) {
      .WALD_DF_LISTWISE_REMEDY
    } else {
      paste0("Use missing = \"fiml\" with lavaan's default test ",
             "statistic and standard errors (or estimator = \"MLR\"), ",
             "or the test with missing = \"listwise\", or pass a ",
             "number of degrees of freedom.")
    }), arg = NULL)
}

# .wald_df_fit_scope(fit)
#
# The calibration study's scope checked on the fit itself, whatever
# reached lavaan through the caller's arguments (a cluster or sampling
# weights passed on through fa_confirmatory()'s ..., say): normal theory or
# Huber-White standard errors (robust.sem only where no value is missing),
# an information matrix the study used or measured, and no sampling
# weights. The study's fits (lavaan with missing = "ml", the cells without
# missing values included) took the observed information from the Hessian
# of the fitted model (information = "observed", observed.information =
# "hessian"), lavaan's default under FIML, and that is the only one
# accepted when values are missing. With complete data lavaan's default
# for ML is the expected information, which complete_data_se.R compares
# with the study's, so it is accepted there as well. The first order
# information and the observed information of the unrestricted model
# (observed.information = "h1"; its loading standard errors differ from
# the Hessian's by up to 1.3 percent on holzinger_swineford's five
# spatial tests by FIML) were not evaluated. The FMI machinery's own refusals
# (.fmi_check_fit()) apply on top when values are missing.

.wald_df_fit_scope <- function(fit) {
  opt <- lavaan::lavInspect(fit, "options")
  complete <- .wald_df_fit_complete(fit)
  se <- as.character(opt$se)
  if (!se %in% c("standard", "robust.huber.white") &&
      !(complete && identical(se, "robust.sem"))) {
    .dmar_refuse_argument(c(
      "wald_df = \"missing_information\" is not available for this ",
      "fit: its standard errors are lavaan's se = \"", se, "\", and ",
      .WALD_DF_SE_SCOPE, ". Pass a number of degrees of freedom."), arg = NULL)
  }
  info <- as.character(opt$information[1L])
  obs_info <- as.character(opt$observed.information[1L])
  hessian <- identical(info, "observed") && identical(obs_info, "hessian")
  if (!(hessian || (complete && identical(info, "expected")))) {
    used <- if (identical(info, "observed")) {
      paste0("information = \"observed\" with observed.information = \"",
             obs_info, "\"")
    } else {
      paste0("information = \"", info, "\"")
    }
    .dmar_refuse_argument(c(
      "wald_df = \"missing_information\" is not available for this ",
      "fit: its standard errors come from lavaan's ", used, ", and ",
      if (complete) {
        paste0("with complete data the keyword takes the expected ",
               "information (lavaan's default there) or the observed ",
               "information from the Hessian of the fitted model ",
               "(information = \"observed\" with observed.information ",
               "= \"hessian\"), the one the calibration study used")
      } else {
        paste0("under full information maximum likelihood the ",
               "calibration study used the observed information from ",
               "the Hessian of the fitted model (information = ",
               "\"observed\" with observed.information = \"hessian\", ",
               "lavaan's default there)")
      }, ". Pass a number of degrees of freedom."), arg = NULL)
  }
  if (.fmi_has_sampling_weights(fit, opt)) {
    .dmar_refuse_argument(c(
      "wald_df = \"missing_information\" is not available for this ",
      "fit: it used sampling weights, whose pseudo-likelihood the ",
      "calibration study did not evaluate. Pass a number of degrees of ",
      "freedom."), arg = NULL)
  }
  invisible(TRUE)
}

# .wald_df_require_wald(wald_df, applies, reason)
#
# Stops when wald_df asks for a t reference (a finite number or the
# keyword) and the call builds no Wald interval for it to change. `reason`
# completes the sentence, saying what the call does instead and where
# wald_df applies.
.wald_df_require_wald <- function(wald_df, applies, reason) {
  if (applies || .wald_df_is_normal(wald_df)) return(invisible(TRUE))
  .dmar_refuse_argument(c(
    "'wald_df' sets the reference distribution of a Wald interval, ",
    reason, " Leave wald_df at Inf, its default."), arg = "wald_df")
}

# .wald_df_resolve(wald_df, fn, fit = NULL, rows = NULL, context = list())
#
# The degrees of freedom of each Wald interval the call builds. `rows`
# names the parameters (a data frame with columns lhs, op, rhs, and
# optionally group, one row per interval, as parameterEstimates() gives
# them; a defined parameter is matched by its lhs alone), `fit` is the
# lavaan fit they belong to (NULL for a closed form computed from complete
# data), and `context` carries what a rule may need (N, the counts,
# `terms`, the names the intervals are reported under, `omega`, the
# names of the omega rows, and `requested_missing`, the caller's FIML
# request for raw data, which .wald_df_fit_fiml() holds the fit to). A
# number is returned once per row. The
# keyword is computed by the function's rule for the rows the rule covers
# (all of them unless the rule says otherwise): each parameter's fraction
# of missing information from the fit (.wald_df_fit_fmi()), then
# df_missing(), whose table, over the covered rows, is returned alongside
# so the result can say how its degrees of freedom were reached. A row the
# rule does not cover gets Inf, the normal reference. A calibrated rule
# first checks the fit against the calibration study's scope
# (.wald_df_fit_scope()), and every rule checks that the fit is the FIML
# fit the call asked for (.wald_df_fit_fiml()).
# Returns list(df = named numeric, table = NULL or the df_missing() table).
.wald_df_resolve <- function(wald_df, fn, fit = NULL, rows = NULL,
                             context = list()) {
  n <- if (is.null(rows)) 1L else nrow(rows)
  terms <- context$terms
  if (length(terms) != n) terms <- paste0("parameter_", seq_len(n))
  if (is.numeric(wald_df)) {
    return(list(df = stats::setNames(rep(wald_df, n), terms),
                table = NULL))
  }
  rule <- .wald_df_rule(fn)
  if (is.null(rule)) {
    .dmar_refuse_argument(c(
      "wald_df = \"missing_information\" has no rule for ", fn, "()."),
      arg = NULL)
  }
  if (!is.null(fit)) {
    .wald_df_fit_fiml(fit, context$requested_missing,
                      calibrated = isTRUE(rule$calibrated))
  }
  if (isTRUE(rule$calibrated) && !is.null(fit)) .wald_df_fit_scope(fit)
  covered <- if (is.null(rule$rows)) rep(TRUE, n) else
    rule$rows(rows, context) %in% TRUE
  df <- stats::setNames(rep(Inf, n), terms)
  if (!any(covered)) return(list(df = df, table = NULL))
  rows_c <- rows[covered, , drop = FALSE]
  info <- .wald_df_fit_fmi(fit, rows_c, context)
  fmi <- stats::setNames(info$fmi, terms[covered])
  table <- if (identical(rule$method, "complete_df")) {
    df_missing(fmi, df_complete = rule$df_complete(rows_c, context),
               method = "complete_df")
  } else {
    df_missing(fmi, N = info$N, k = rule$k(rows_c, context))
  }
  df[covered] <- table$df
  list(df = df, table = table)
}

# .wald_df_attach(x, wald_df, table = NULL, rows = NULL)
#
# Records on a table the reference distribution of its Wald intervals, in
# the attributes the printed line under it reads (.wald_df_footer_line()):
# for a number, the "wald_df" attribute is that number; for the keyword,
# `table` is the df_missing() table .wald_df_resolve() returned, and
# "wald_df" becomes the degrees of freedom of each of its parameters that
# x reports (named by term, once per term, since an equality constraint
# repeats a label), "wald_df_table" the df_missing() table of those
# parameters, and "wald_df_rows" `rows`, what they are ("loadings",
# "omega"), for the printed line. No parameter of the table in x leaves
# x on the normal reference.
.wald_df_attach <- function(x, wald_df, table = NULL, rows = NULL) {
  if (is.null(table)) {
    attr(x, "wald_df") <- wald_df
    return(x)
  }
  keep <- table$term %in% x$term & !duplicated(table$term)
  if (!any(keep)) {
    attr(x, "wald_df") <- Inf
    return(x)
  }
  tab <- df_missing(stats::setNames(table$fmi[keep], table$term[keep]),
                    N = table$N[keep], k = table$k[keep])
  attr(x, "wald_df") <- stats::setNames(tab$df, tab$term)
  attr(x, "wald_df_table") <- tab
  if (inherits(x, "dmar_fa_confirmatory")) {
    attr(x, "wald_calibration_note") <- .wald_eligibility_note(tab)
  }
  if (!is.null(rows)) attr(x, "wald_df_rows") <- rows
  x
}

.wald_eligibility_note <- function(table) {
  if (is.null(table) || !any(table$fmi > 1e-8, na.rm = TRUE)) return(NULL)
  paste0("The FMI degrees-of-freedom rule is approximate. Its calibration ",
    "study admitted data sets that the affine-collapse screening sets ",
    "aside, so coverage under that screening has not been established. ",
    "See the method's help for its scope and conditional evidence.")
}

# .wald_df_column(wald_df, n, keys = NULL)
#
# The wald_df column of a tidy() table: the degrees of freedom of the t
# reference of each of its `n` rows, from the "wald_df" record of the
# result (.wald_df_attach()), so that a row on the t reference can be told
# from one on the normal reference after rows of several results are
# stacked (the record is an attribute, which rbind() keeps from the first
# frame only). NULL, a non-number, or Inf: Inf for every row, the normal
# reference (or an interval that is not a Wald interval). One number: that
# number for every row. A named vector (the keyword "missing_information"):
# matched by `keys`, the rows' terms, and Inf for a row it does not name,
# which keeps the normal reference; an unnamed vector of length n is taken
# in row order. Returns a numeric vector of length n.
.wald_df_column <- function(wald_df, n, keys = NULL) {
  if (!is.numeric(wald_df) || !length(wald_df)) return(rep(Inf, n))
  wald_df <- stats::setNames(as.numeric(wald_df), names(wald_df))
  if (!is.null(names(wald_df)) && !is.null(keys) &&
      any(keys %in% names(wald_df))) {
    out <- unname(wald_df[match(keys, names(wald_df))])
    out[is.na(out)] <- Inf
    return(out)
  }
  if (length(wald_df) == 1L) return(rep(unname(wald_df), n))
  if (length(wald_df) == n) return(unname(wald_df))
  rep(NA_real_, n)
}

# .wald_df_apply(tbl, wald_df, conf_level, est = "estimate", se = "se",
#                z = "test_statistic", p = "p_value", lower = "ci_lower",
#                upper = "ci_upper")
#
# Recomputes the Wald intervals and p-values of a parameter table on the
# t reference, from its own estimates, standard errors, and z ratios:
# est -/+ qt(1 - (1 - C)/2, df) se and 2 pt(-|z|, df), with wald_df one
# value per row (or one for all). Rows without a standard error (fit
# rows, categorical omega, HTMT) or without a z ratio (fixed parameters)
# keep what they had, and so do rows whose wald_df is infinite (the rows
# the keyword's rule does not cover), whose intervals and p-values stay
# lavaan's own to the last bit. Used by fa_confirmatory() for its own table (DMAR's
# column names) and for lavaan's standardized solution (lavaan's names).
.wald_df_apply <- function(tbl, wald_df, conf_level, est = "estimate",
                           se = "se", z = "test_statistic", p = "p_value",
                           lower = "ci_lower", upper = "ci_upper") {
  n <- nrow(tbl)
  wald_df <- rep_len(as.numeric(wald_df), n)
  t_ref <- is.finite(wald_df)
  half <- .wald_crit(conf_level, wald_df) * tbl[[se]]
  has_se <- !is.na(tbl[[se]]) & t_ref
  tbl[[lower]][has_se] <- (tbl[[est]] - half)[has_se]
  tbl[[upper]][has_se] <- (tbl[[est]] + half)[has_se]
  has_z <- !is.na(tbl[[z]]) & t_ref
  tbl[[p]][has_z] <- .wald_p(tbl[[z]], wald_df)[has_z]
  tbl
}

# .wald_df_fit_fmi(fit, rows, context = list())
#
# The fraction of missing information and the N basis of each parameter
# named in `rows`. With no missing values among the cases the fit used
# (complete data, or a listwise deletion fit) the complete and observed
# data information coincide, so every FMI is zero by definition and is set
# to 0 exactly rather than estimated, which makes the effective sample
# size rule return the complete data degrees of freedom N - k to the last
# bit. Otherwise the FMI is the unstructured estimate of
# fraction_missing_information() (.fmi_check_fit() and .fmi_core()), whose
# refusals (a pairwise or two-stage fit, ordered indicators, an estimator
# outside the maximum likelihood family) stop the call with their reason.
# For a robust (MLR) fit it is the normal theory FMI at its estimates,
# which are the ML estimates. A NULL fit is a closed form computed from
# complete data, whose N comes from context$N.
.wald_df_fit_fmi <- function(fit, rows, context = list()) {
  n <- nrow(rows)
  if (is.null(fit)) {
    return(list(fmi = rep(0, n), N = rep(as.numeric(context$N), n)))
  }
  nobs <- as.numeric(unlist(lavaan::lavInspect(fit, "nobs")))
  complete <- .wald_df_fit_complete(fit)
  pt <- if (complete) {
    lavaan::parTable(fit)
  } else {
    tryCatch({
      .fmi_check_fit(fit)
      .fmi_core(fit, "unstructured")
    }, error = function(e) {
      # The condition keeps its class. Where the data define no fraction
      # (no unrestricted estimate, or a pair of variables never observed
      # together), the remedy is the keyword's alternatives.
      e$message <- paste0(
        "wald_df = \"missing_information\" takes each coefficient's ",
        "fraction of missing information from the fit. ",
        conditionMessage(e),
        if (inherits(e, c("dmar_fmi_saturated_singular",
                          "dmar_fmi_never_together"))) {
          paste0(" Pass a number of degrees of freedom, or leave wald_df ",
                 "at Inf for the normal reference.")
        })
      e$call <- NULL
      stop(e)
    })
  }
  key <- function(lhs, op, rhs, group) {
    ifelse(op == ":=", lhs, paste(lhs, op, rhs, group))
  }
  row_group <- if (is.null(rows$group)) rep(1L, n) else rows$group
  idx <- match(key(rows$lhs, rows$op, rows$rhs, row_group),
               key(pt$lhs, pt$op, pt$rhs, pt$group))
  if (anyNA(idx)) {
    .dmar_refuse_argument(c(
      "The parameters of a Wald interval could not be found in the ",
      "fit's parameter table."), arg = NULL)
  }
  list(fmi = if (complete) rep(0, n) else pt$fmi[idx],
       N = .fmi_row_N(pt, nobs)[idx])
}

# .wald_df_phrase(wd)
#
# "25 degrees of freedom", or "118 to 136.5 degrees of freedom" when the
# parameters differ, for the printed lines below. Display only: the
# values stored on the result keep full precision.
.wald_df_phrase <- function(wd) {
  wd <- as.numeric(wd)
  fmt <- function(v) format(v, digits = 4L, trim = TRUE)
  if (length(unique(wd)) == 1L) {
    return(sprintf("%s degrees of freedom", fmt(wd[1L])))
  }
  sprintf("%s to %s degrees of freedom, by parameter",
          fmt(min(wd, na.rm = TRUE)), fmt(max(wd, na.rm = TRUE)))
}

# .wald_df_footer_line(x)
#
# The line a DMAR table writes under itself when its Wald intervals refer
# to a t distribution, that is when its wald_df attribute is finite: the
# reliability tables and the fa_confirmatory() tables. .dmar_tbl_footer() calls it,
# so the console and a knitted report show the same line. NULL for the
# normal reference and for a table without the attribute. Under the
# keyword the table carries the df_missing() table of its parameters
# (attribute "wald_df_table"), and the line names the rule, its N, and
# (when it is one value) the FMI, and in a fa_confirmatory() table which rows took
# the t reference (attribute "wald_df_rows") and that the others kept the
# normal one.
.wald_df_footer_line <- function(x) {
  wd <- attr(x, "wald_df", exact = TRUE)
  if (!is.numeric(wd) || !length(wd) || .wald_df_is_normal(wd)) {
    return(NULL)
  }
  wide <- "p_value" %in% names(x)
  tab <- attr(x, "wald_df_table", exact = TRUE)
  if (is.data.frame(tab) && nrow(tab) > 0L) {
    return(.wald_df_rule_line(tab, wide,
                              attr(x, "wald_df_rows", exact = TRUE)))
  }
  what <- if (wide) "Wald intervals and p-values" else "Wald interval"
  sprintf("%s: t reference on %s (wald_df), in place of the normal",
          what, .wald_df_phrase(wd))
}

# .wald_df_rule_line(tab, wide, rows = NULL)
#
# The line under a table whose degrees of freedom came from the keyword:
# "Wald interval: t reference on 245.3 degrees of freedom, from von
# Hippel's (2016) max(3, N(1 - FMI) - 1) with N = 301 and FMI = 0.1823
# (wald_df = "missing_information"), in place of the normal", and for a
# fa_confirmatory() table the rows it applies to, with the others on the normal
# reference. Display only: the attributes keep full precision.
.wald_df_rule_line <- function(tab, wide, rows = NULL) {
  fmt <- function(v) format(v, digits = 4L, trim = TRUE)
  k <- unique(tab$k)
  rule <- sprintf("max(3, N(1 - FMI) - %s)", paste(fmt(k), collapse = " or "))
  basis <- sprintf("N = %s", paste(fmt(unique(tab$N)), collapse = ", "))
  fmi <- tab$fmi
  # An FMI within 1e-10 of zero differs from zero only by rounding error,
  # as .wald_df_summary_lines() shows it.
  fmi[!is.na(fmi) & abs(fmi) < 1e-10] <- 0
  if (length(unique(fmi)) == 1L) {
    basis <- paste0(basis, " and FMI = ", fmt(fmi[1L]))
  }
  source <- sprintf(paste0("t reference on %s, from von Hippel's (2016) %s ",
                           "with %s (wald_df = \"missing_information\")"),
                    .wald_df_phrase(tab$df), rule, basis)
  if (!wide) {
    return(paste0("Wald interval: ", source, ", in place of the normal"))
  }
  of <- if (length(rows)) {
    paste0(" of ", paste(ifelse(rows == "loadings", "the loadings", rows),
                         collapse = " and "))
  } else {
    ""
  }
  paste0("Wald intervals and p-values", of, ": ", source,
         "; the other rows keep the normal reference")
}

# .wald_df_coef_headers(wd)
#
# The column headers of the coefficient table summary() of an mlmr() or
# mlmr_mv() fit prints. They name the reference distribution of the
# p-values, as summary.lm() and summary.glm() do: "z value" and
# "Pr(>|z|)" for the normal (wald_df = Inf, or no record), "t value" and
# "Pr(>|t|)" when the p-values are 2 P(T > |z|) on wald_df degrees of
# freedom. The ratio is the same estimate over standard error either way,
# and the stored column is test_statistic under either reference.
.wald_df_coef_headers <- function(wd) {
  if (is.null(wd) || .wald_df_is_normal(wd)) {
    return(c("Estimate", "Std. Error", "z value", "Pr(>|z|)"))
  }
  c("Estimate", "Std. Error", "t value", "Pr(>|t|)")
}

# .wald_df_summary_lines(wd, table = NULL, digits = 3L)
#
# The lines summary() of an mlmr() or mlmr_mv() fit prints under its
# coefficient table (under the last outcome's, for mlmr_mv()) when its
# confidence intervals and p-values refer to a t distribution. With the
# keyword, `table` is the df_missing() table the degrees of freedom came
# from, and the lines name the rule and give each coefficient's FMI and
# degrees of freedom. character(0) for the normal reference.
.wald_df_summary_lines <- function(wd, table = NULL, digits = 3L) {
  if (is.null(wd) || .wald_df_is_normal(wd)) return(character(0))
  if (is.null(table)) {
    return(sprintf(paste0("Wald intervals and p-values: t reference on %s ",
                          "(wald_df), in place of the normal"),
                   .wald_df_phrase(wd)))
  }
  fmt <- function(v) format(unique(v), digits = 4L, trim = TRUE)
  head <- sprintf(paste0(
    "Wald intervals and p-values: t reference on von Hippel's (2016) ",
    "missing information degrees of freedom, max(3, N(1 - FMI) - k), ",
    "with N = %s and k = %s:"),
    paste(fmt(table$N), collapse = ", "),
    paste(fmt(table$k), collapse = ", "))
  # An FMI within 1e-10 of zero (a coefficient whose outcome is fully
  # observed, say) differs from zero only by rounding error and is shown
  # as 0, as print.dmar_fmi() shows it.
  fmi <- table$fmi
  fmi[!is.na(fmi) & abs(fmi) < 1e-10] <- 0
  shown <- data.frame(term = table$term,
                      fmi = format(fmi, digits = digits),
                      df = format(table$df, digits = digits + 1L),
                      stringsAsFactors = FALSE)
  shown$df[table$at_floor %in% TRUE] <-
    paste0(shown$df[table$at_floor %in% TRUE], " (floor)")
  body <- utils::capture.output(print(shown, row.names = FALSE,
                                      right = FALSE))
  c(strwrap(head, width = max(40L, getOption("width", 80L) - 2L)),
    paste0("  ", body))
}
