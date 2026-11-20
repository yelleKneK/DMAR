# The engine behind fa_confirmatory(). The exported function validates
# nothing itself: it resolves `output`, records whether `missing` was left
# at its default, and hands every argument to .fa_cfa_engine(), which
# checks the call, writes the lavaan syntax for the requested structure,
# fits it, and assembles the table. The structure is named by `structure`:
# "correlated" is the model of fa_confirmatory(), simple structure with
# correlated (or, with correlated_factors = FALSE, orthogonal) factors;
# "bifactor" is fa_bifactor()'s, and "second_order" fa_second_order()'s.
# The syntax of a structure is written by its own builder
# (.fa_cfa_correlated_syntax()), and the fit rows and the attributes that
# describe them by .fa_cfa_fit_rows() and .fa_cfa_fit_attributes(), so
# that another structure fitted by the same engine reports them in the
# same way. The helpers below the engine (.cfa_k_*) carry the names they
# had when the function was cfa_k().
#
# Every argument of the engine follows `...`, so a name passed through to
# lavaan can never partially match one of them, and `...` is evaluated
# where fa_confirmatory() evaluated it.
.fa_cfa_engine <- function(..., data, factors, S, N, M, equal_loading,
                           equal_intercept, equal_error, correlated_factors,
                           meanstructure, estimator, missing, ordered, se,
                           conf_level, rmsea_conf_level, output, wald_df,
                           missing_default,
                           structure = "correlated", general = NULL,
                           general_only = NULL, second_order = NULL,
                           equal_second_order_loading = FALSE) {
  # The function whose name keys the wald_df rules and the two-stage
  # notices. Under "bifactor", `factors` holds the group factors, and
  # `general` names the general factor, which loads on their items and on
  # the items of `general_only`. Under "second_order", `factors` holds the
  # first-order factors, and `second_order` maps each second-order factor
  # to its first-order factors. The two are "layered" structures, fit by
  # maximum likelihood to continuous items with the indices of a general
  # factor and its group factors.
  fn <- switch(structure, correlated = "fa_confirmatory",
               bifactor = "fa_bifactor",
               second_order = "fa_second_order",
               .dmar_refuse_argument(c(
                 "Unknown structure '", structure, "'."), arg = NULL))
  bifactor <- identical(structure, "bifactor")
  second <- identical(structure, "second_order")
  layered <- bifactor || second
  se_auto <- is.null(se)
  if (se_auto) se <- if (is.character(estimator) && length(estimator) == 1L && !is.na(estimator)) {
    switch(estimator, MLR = "robust.huber.white", WLSMV = "robust.sem", "standard")
  } else "standard"
  if (!is.character(se) || length(se) != 1L || is.na(se)) {
    .dmar_refuse_argument(
      "'se' must be NULL or one standard-error method name.", arg = "se")
  }

  wald_df <- .wald_df_check(wald_df, fn)
  # ---- one group --------------------------------------------------------
  # Every parameter carries a label, and lavaan holds a label shared by
  # several groups equal across them, so a 'group' passed through ... would
  # fit one set of measurement parameters for every group, each factor's
  # variance 1 (and mean 0) in each, with lavaan's warning; through 1.0.1
  # that model was fit and its table repeated each labeled row once per
  # group. A name that partially matches lavaan's 'group' formal ("gr")
  # would reach it too, so any prefix of "group" is refused.
  dots <- list(...)
  dot_names <- names(dots)
  if (layered && !is.null(dot_names) &&
      any(nzchar(dot_names) & startsWith("group", dot_names))) {
    .dmar_refuse_argument(c(
      "'group' is not available in ", fn, "(): its parameters are ",
      "labeled, and lavaan holds a label shared by several groups equal ",
      "across them, so the groups would share every loading, error ",
      "variance, and intercept (one pooled model), and the indices would ",
      "belong to that pooled model rather than to any one group. Fit the ",
      "model to each group's data separately."), arg = "group")
  }
  if (!is.null(dot_names) &&
      any(nzchar(dot_names) & startsWith("group", dot_names))) {
    .dmar_refuse_argument(c(
      "'group' is not available in fa_confirmatory() (nor in cfa_1() and ",
      "cfa_2(), which call it): its parameters are labeled, and lavaan ",
      "holds a label shared by several groups equal across them, so ",
      "the groups would share every loading, error variance, ",
      "intercept, and factor correlation (one pooled model), and omega, ",
      "ave, and H would belong to that pooled model rather than to any ",
      "one group. For the model in several groups use ",
      "measurement_invariance(), whose configural model frees every ",
      "parameter in each group; reliability_H(fit = , factor = , ",
      "group = ) on that fit gives H for one group, and ",
      "reliability_omega() on one group's rows gives its omega."),
      arg = "group")
  }
  if (!requireNamespace("lavaan", quietly = TRUE)) {
    .dmar_refuse(c(
      "The package 'lavaan' is needed. Install it with ",
      "install.packages(\"lavaan\")."), "package_required", package = "lavaan")
  }

  # ---- validate the factor specification -------------------------------
  if (!is.list(factors) || length(factors) < 1L ||
      is.null(names(factors)) || any(names(factors) == "")) {
    .dmar_refuse_argument(c(
      "'factors' must be a named list: ",
      "list(factor_name = c(\"item_1\", \"item_2\", ...))."), arg = "factors")
  }
  factor_names <- names(factors)
  if (anyDuplicated(factor_names)) {
    .dmar_refuse_argument("Factor names must be unique.", arg = NULL)
  }
  if (!all(factor_names == make.names(factor_names))) {
    .dmar_refuse_argument(c(
      "Factor names must be syntactically valid R names (they become ",
      "lavaan labels); rename, for example, 'my factor' to 'my_factor'."),
      arg = NULL)
  }
  if (bifactor) {
    # The general factor's name and items, the items of the group factors,
    # and the refusals that need nothing but the names
    # (.fa_bifactor_items()).
    items <- .fa_bifactor_items(factors, general, general_only)
  } else if (second) {
    # The same for the first-order factors and the second-order factors
    # above them (.fa_second_order_items()).
    items <- .fa_second_order_items(factors, second_order)
  } else {
    for (f in factor_names) {
      if (!is.character(factors[[f]]) || length(factors[[f]]) < 2L) {
        .dmar_refuse_argument(
          sprintf("Factor '%s' must name two or more items.", f), arg = NULL)
      }
    }
    items <- unlist(factors, use.names = FALSE)
    if (anyDuplicated(items)) {
      .dmar_refuse_argument(c(
        "Each item may load on exactly one factor (simple structure); ",
        "duplicated: ",
        paste(unique(items[duplicated(items)]), collapse = ", "), "."),
        arg = NULL)
    }
  }


  # ---- expand the descriptor flags to one value per factor -------------
  # Under "bifactor" the loadings flag also takes the general factor
  # (equal general loadings); the error variances and intercepts belong to
  # the items, which the group factors partition, so their flags name the
  # group factors only. Under "second_order" the three name the
  # first-order factors, and equal_second_order_loading the second-order
  # factors.
  expand_flag <- function(flag, arg_name, flag_names = factor_names) {
    if (!is.logical(flag) || anyNA(flag)) {
      .dmar_refuse_argument(c(
        "'", arg_name, "' must be logical: a single value, or a named ",
        "vector with one element per factor."), arg = NULL)
    }
    if (length(flag) == 1L && is.null(names(flag))) {
      return(stats::setNames(rep(flag, length(flag_names)), flag_names))
    }
    if (is.null(names(flag)) || !setequal(names(flag), flag_names)) {
      .dmar_refuse_argument(c(
        "A vector '", arg_name, "' must be named with exactly the ",
        if (bifactor && arg_name == "equal_loading") {
          "names of the general factor and the group factors."
        } else if (bifactor) {
          "group factor names."
        } else if (second && arg_name == "equal_second_order_loading") {
          "second-order factor names."
        } else if (second) {
          "first-order factor names."
        } else {
          "factor names."
        }), arg = NULL)
    }
    flag[flag_names]
  }
  eq_load <- expand_flag(equal_loading, "equal_loading",
                         if (bifactor) c(general, factor_names) else
                           factor_names)
  eq_int  <- expand_flag(equal_intercept, "equal_intercept")
  eq_err  <- expand_flag(equal_error, "equal_error")
  eq_gamma <- if (second) {
    expand_flag(equal_second_order_loading, "equal_second_order_loading",
                names(second_order))
  }

  # (The identification of the bifactor and second-order structures is
  # checked after the arguments and the input, below.)
  if (!layered && length(factors) == 1L && length(factors[[1L]]) < 3L &&
      !(length(factors[[1L]]) == 2L && isTRUE(unname(eq_load[1L])))) {
    # A single congeneric factor needs three items; with equal loadings
    # the two-item factor is just identified (the tau-equivalent case the
    # reliability family fits for a two-item scale).
    .dmar_refuse_argument(c(
      "A single factor needs three or more items for identification ",
      "(two items are identified only with equal_loading = TRUE)."), arg = NULL)
  }

  # ---- whitelists, as in cfa_1() ---------------------------------------
  # The bifactor and second-order models are fit by maximum likelihood to
  # continuous items, by FIML or after listwise deletion; fa_bifactor()
  # and fa_second_order() refuse ordered items themselves.
  estimator_choices <- if (layered) c("ML", "MLR") else
    c("ML", "MLR", "WLS", "WLSMV", "GLS")
  missing_choices   <- if (layered) c("fiml", "ml", "listwise") else
    c("fiml", "ml", "listwise", "pairwise", .TWO_STAGE_VALUES)
  if (!is.character(estimator) || length(estimator) != 1L ||
      !(estimator %in% estimator_choices)) {
    .dmar_refuse_argument(c(
      "'estimator' must be one of ",
      paste(dQuote(estimator_choices), collapse = ", "), ".",
      if (layered) {
        paste0(" ", fn, "() fits continuous items by maximum ",
               "likelihood, with robust standard errors and test under ",
               "\"MLR\".")
      }), arg = "estimator")
  }
  if (!is.character(missing) || length(missing) != 1L ||
      !(missing %in% missing_choices)) {
    .dmar_refuse_argument(c(
      "'missing' must be one of ",
      paste(dQuote(missing_choices), collapse = ", "), ".",
      if (layered) {
        paste0(" ", fn, "() fits continuous items by full ",
               "information maximum likelihood or after listwise ",
               "deletion.")
      }), arg = "missing")
  }
  missing <- .cfa_k_resolve_missing(missing, missing_default, data, items,
                                    ordered, estimator, se, dots$test)
  if (!is.null(data) && isTRUE(missing %in% c("ml", "fiml")) &&
      identical(se, "robust.sem")) {
    .dmar_refuse_argument(
      'se = "robust.sem" requires complete data and is incompatible with FIML. Use se = "robust.huber.white" with FIML, or explicitly request listwise deletion.',
      arg = NULL)
  }
  # lavaan's value for a two-stage request ("two.stage" or
  # "robust.two.stage"), NULL otherwise.
  ts_missing <- .ts_lavaan_missing(missing)
  two_stage <- !is.null(ts_missing)
  if (!is.numeric(conf_level) || length(conf_level) != 1L ||
      is.na(conf_level) || conf_level <= 0 || conf_level >= 1) {
    .dmar_refuse_argument(
      "'conf_level' must be a single number in (0, 1).", arg = "conf_level")
  }
  # The level of the RMSEA interval of the fit rows, apart from the level
  # of the parameter intervals: 0.90 by default, the level at which DMAR
  # 1.0.0's cfa_k() reported it (lavaan's), whose limits decide the tests
  # of close and not close fit at .05.
  if (!is.numeric(rmsea_conf_level) || length(rmsea_conf_level) != 1L ||
      is.na(rmsea_conf_level) || rmsea_conf_level <= 0 ||
      rmsea_conf_level >= 1) {
    .dmar_refuse_argument(
      "'rmsea_conf_level' must be a single number in (0, 1).",
      arg = "rmsea_conf_level")
  }
  # A t reference distribution changes the Wald intervals and p-values of
  # the parameter tables (output "verbose", "measurement", and
  # "standardized"); the lavaan fit and its summary keep lavaan's normal
  # reference, and with bootstrap or no standard errors there is no Wald
  # interval to change. The keyword "missing_information" is narrower: the
  # calibration study behind it covers the loadings (and omega on its own
  # scale, which fa_confirmatory() no longer uses for omega's limits) of a
  # maximum likelihood fit of continuous items, by FIML or with complete data, with
  # normal theory or Huber-White standard errors, and anything else is
  # refused here with its reason (.wald_df_scope(),
  # R/missing_information_internals.R).
  .wald_df_scope(wald_df, fn, list(
    output = output, se = se, estimator = estimator, missing = missing,
    two_stage = two_stage, ordered = ordered,
    complete = is.null(data) || missing == "listwise", test = dots$test))
  # The FIML request the fit is held to under the keyword (lavaan fits the
  # complete cases instead with some standard errors and test statistics;
  # .wald_df_fit_scope()). A covariance matrix is fit as complete data.
  requested_missing <- if (!is.null(data) && missing %in% c("ml", "fiml")) {
    missing
  }
  .wald_df_require_wald(wald_df, output %in% c("verbose", "measurement",
                                               "standardized"),
    paste0("and output = \"", output, "\" returns lavaan's own object, ",
           "whose intervals keep lavaan's normal reference; wald_df ",
           "applies to the tables of output = \"verbose\", ",
           "\"measurement\", and \"standardized\"."))
  .wald_df_require_wald(wald_df, !se %in% c("none", "bootstrap"),
    if (se == "none") {
      paste0("and se = \"none\" gives no standard errors, so the tables ",
             "have no Wald interval.")
    } else {
      paste0("and with se = \"bootstrap\" lavaan reports percentile ",
             "bootstrap intervals, not Wald intervals.")
    })

  # ---- input: raw data ('data') or covariance matrix ('S') -------------
  if (is.null(data) == is.null(S)) {
    .dmar_refuse_argument(c(
      "Supply exactly one of 'data' (raw data) or 'S' (a symmetric ",
      "covariance matrix, with 'N')."), arg = c("data", "S", "N"))
  }
  input_is_raw_data <- !is.null(data)

  if (input_is_raw_data) {
    if (!is.data.frame(data) && !is.matrix(data)) {
      .dmar_refuse_argument(c(
        "'data' must be a data frame or matrix of raw data (rows are ",
        "respondents, columns include the items)."), arg = "data")
    }
    if (is.matrix(data) && nrow(data) == ncol(data) &&
        isSymmetric(unname(data), tol = 1e-05)) {
      .dmar_refuse_argument(c(
        "'data' looks like a covariance matrix (square and symmetric); ",
        "pass it through 'S', with 'N'."), arg = c("data", "S", "N"))
    }
    raw_data <- as.data.frame(data)
    # An item the factors name that the input lacks is the mistake an
    # `items` argument makes in fa_exploratory() and its siblings, and is
    # refused with the same class and reason (.fa_items_absent(),
    # R/factor_analysis_internals.R: "unknown", or "form" for an NA name).
    # The checks of the factors above (too few items, an item on two
    # factors) concern the model rather than the naming of its items and
    # do not take this class.
    absent <- setdiff(items, names(raw_data))
    if (length(absent) > 0L) {
      .fa_items_absent(c("Items not found in the data: ",
                         paste(absent, collapse = ", "), "."), absent)
    }
    lavaan_args <- list(data = raw_data)
  } else {
    if (!is.matrix(S) || !isSymmetric(unname(S), tol = 1e-05)) {
      .dmar_refuse_argument(c(
        "'S' must be a symmetric covariance matrix; raw data are ",
        "passed through 'data'."), arg = c("S", "data"))
    }
    if (is.null(N)) {
      .dmar_refuse_argument(c(
        "Sample size 'N' must be provided when input is a covariance ",
        "matrix."), arg = "N")
    }
    if (is.null(colnames(S))) {
      .dmar_refuse_argument(
        "A covariance matrix 'S' must have dimnames naming the items.",
        arg = "S")
    }
    rownames(S) <- colnames(S)
    absent <- setdiff(items, colnames(S))
    if (length(absent) > 0L) {
      .fa_items_absent(c("Items not found in the covariance matrix: ",
                         paste(absent, collapse = ", "), "."), absent)
    }
    # The names of S as given, for the lavaan code of a refusal (see
    # "a likelihood with no maximum" below).
    S_names <- colnames(S)
    S <- S[items, items, drop = FALSE]
    lavaan_args <- list(sample.cov = S, sample.nobs = N)
  }

  # ---- ordered-categorical items ---------------------------------------
  ordered_items <- NULL
  if (!is.null(ordered)) {
    if (isTRUE(ordered)) {
      ordered_items <- items
    } else if (is.character(ordered)) {
      absent <- setdiff(ordered, items)
      if (length(absent) > 0L) {
        .dmar_refuse_argument(c(
          "'ordered' names items that are not in 'factors': ",
          paste(absent, collapse = ", "), "."), arg = c("ordered", "factors"))
      }
      ordered_items <- unique(ordered)
    } else {
      .dmar_refuse_argument(c(
        "'ordered' must be TRUE (all items) or a character vector of ",
        "item names."), arg = "ordered")
    }
    if (!input_is_raw_data) {
      .dmar_refuse_argument(c(
        "Ordered items require raw data: thresholds and polychoric ",
        "correlations cannot be estimated from a covariance matrix."),
        arg = NULL)
    }
    # A factor must be all ordered or all continuous: the sum score of a
    # mixed factor has no single metric for its omega.
    for (f in factor_names) {
      n_ord <- sum(factors[[f]] %in% ordered_items)
      if (n_ord > 0L && n_ord < length(factors[[f]])) {
        .dmar_refuse_argument(c(
          "Factor '", f, "' mixes ordered and continuous items; ",
          "declare all of a factor's items ordered, or none."), arg = NULL)
      }
    }
    if (two_stage) {
      .dmar_refuse_argument(c(
        "Two-stage estimation is not available with ordered items: its ",
        "first stage estimates the means and covariances of continuous ",
        "variables under the normal model (lavaan refuses it in the ",
        "categorical setting). Use missing = \"pairwise\" or ",
        "\"listwise\"."), arg = NULL)
    }
    if (estimator %in% c("ML", "MLR", "GLS")) {
      message("Ordered items declared: switching to estimator = \"WLSMV\" ",
              "with robust standard errors.")
      estimator <- "WLSMV"
    }
    # Robust standard errors are the appropriate companion to WLSMV.
    if (se_auto) se <- if (estimator == "WLSMV") "robust.sem" else "standard"
    if (missing %in% c("ml", "fiml")) {
      .dmar_refuse_argument(c(
        "Full information maximum likelihood is not available with ",
        "ordered items; use missing = \"pairwise\" or \"listwise\"."),
        arg = NULL)
    }
  }

  # ---- two-stage estimation ---------------------------------------------
  # Two-stage maximum likelihood is a normal theory method for raw
  # continuous data: it needs incomplete cases to estimate the EM moments
  # from, it fits them by maximum likelihood, and its standard errors are
  # the two-stage sandwich, robust to nonnormality under
  # "robust_two_stage". lavaan would quietly override an MLR estimator or a
  # robust se with a warning; here each is an error that names the option
  # that does what was asked.
  if (two_stage) {
    if (!input_is_raw_data) {
      .dmar_refuse_argument(c(
        "missing = \"", missing, "\" requires raw 'data': two-stage ",
        "estimation starts from the EM estimates of the means and ",
        "covariances of incomplete data, and a covariance matrix has no ",
        "incomplete cases."), arg = "data")
    }
    if (estimator != "ML") {
      .dmar_refuse_argument(c(
        "Two-stage estimation fits the EM estimates of the means and ",
        "covariances by maximum likelihood (estimator = \"ML\"). ",
        if (estimator == "MLR") {
          paste0("For standard errors and a test robust to nonnormality ",
                 "use missing = \"robust_two_stage\" rather than ",
                 "estimator = \"MLR\".")
        } else {
          paste0("estimator = \"", estimator, "\" is not available ",
                 "with it.")
        }), arg = NULL)
    }
    if (!se %in% c("standard", "none", ts_missing)) {
      .dmar_refuse_argument(c(
        "Under two-stage estimation the standard errors are the ",
        "two-stage ones, so 'se' must be \"standard\" (or \"none\").",
        if (se %in% c("robust.sem", "robust.huber.white")) {
          paste0(" For standard errors robust to nonnormality use ",
                 "missing = \"robust_two_stage\".")
        }), arg = "se")
    }
  }
  ordered_factor <- vapply(factor_names, function(f)
    !is.null(ordered_items) && all(factors[[f]] %in% ordered_items),
    logical(1L))

  # ---- mean structure --------------------------------------------------
  if (is.null(meanstructure)) {
    # Full information maximum likelihood requires the mean structure to be
    # estimated with free item intercepts: through the bare lavaan()
    # interface, lavaan forces meanstructure = TRUE for missing = "ml" but
    # leaves int.ov.free = FALSE, so every intercept would be fixed at 0
    # and the loadings would absorb the item means (the same defect
    # cfa_1() carried until its FIML intercept fix; on uncentered data the
    # reported omega was 0.995 against a true 0.802). Turning the mean
    # structure on here writes free nu lines into the syntax, which
    # restores exact agreement with cfa_1() under FIML. Two-stage
    # estimation falls into the same trap through lavaan(), whose stage two
    # then fits the EM means with every intercept at 0 (on three spatial
    # tests with 60 values of one deleted, omega 0.981 and a chi square of
    # 1086 on 3 df, against 0.665 and 0 with free intercepts), so it turns
    # the mean structure on too.
    has_means <- any(eq_int) || !is.null(M) ||
      (input_is_raw_data && (missing %in% c("ml", "fiml") || two_stage))
  } else {
    if (!is.logical(meanstructure) || length(meanstructure) != 1L ||
        is.na(meanstructure)) {
      .dmar_refuse_argument(
        "'meanstructure' must be TRUE, FALSE, or NULL.", arg = "meanstructure")
    }
    has_means <- meanstructure
  }
  if (any(eq_int) && !has_means) {
    .dmar_refuse_argument(c(
      "'equal_intercept' asks about the mean structure; it cannot be ",
      "combined with meanstructure = FALSE."), arg = "equal_intercept")
  }
  if (!has_means && input_is_raw_data &&
      (missing %in% c("ml", "fiml") || two_stage)) {
    .dmar_refuse_argument(c(
      "missing = \"", missing, "\" estimates the mean structure with ",
      "free item intercepts; it cannot be combined with ",
      "meanstructure = FALSE. Drop meanstructure = FALSE, or use ",
      "missing = \"listwise\"."), arg = NULL)
  }
  if (!is.null(ordered_items) && has_means) {
    .dmar_refuse_argument(c(
      "Thresholds replace intercepts for ordered items; ",
      "'equal_intercept', 'M', and meanstructure = TRUE are not ",
      "available with 'ordered'."), arg = c("equal_intercept", "M", "ordered"))
  }
  if (has_means && !input_is_raw_data) {
    if (is.null(M)) {
      .dmar_refuse_argument(c(
        "The mean structure with a covariance matrix requires the item ",
        "means 'M'."), arg = "M")
    }
    if (!is.numeric(M) || is.null(names(M)) || !all(items %in% names(M))) {
      .dmar_refuse_argument(
        "'M' must be a named numeric vector covering every item.", arg = "M")
    }
    lavaan_args$sample.mean <- M[items]
  }

  # ---- identification of the layered structures ------------------------
  # It needs the equality flags (.fa_bifactor_identification(),
  # .fa_second_order_identification()), and it is checked once every
  # argument and the input have passed their checks, so that the bifactor
  # model's warning about empirical underidentification is given only for
  # a model that is then fit.
  if (bifactor) {
    .fa_bifactor_identification(factors, general, general_only, eq_load,
                                eq_err)
  } else if (second) {
    .fa_second_order_identification(factors, second_order, eq_gamma)
  }

  # ---- build the lavaan model syntax -----------------------------------
  model_syntax <- if (bifactor) {
    .fa_bifactor_syntax(factors, general, general_only, eq_load, eq_int,
                        eq_err, has_means)
  } else if (second) {
    .fa_second_order_syntax(factors, second_order, eq_load, eq_int, eq_err,
                            eq_gamma, has_means)
  } else {
    .fa_cfa_correlated_syntax(factors, eq_load, eq_int, eq_err, has_means,
                              correlated_factors)
  }

  # Two factors no case observes together through any of their items: under
  # the factor model (simple structure, so each item loads on one factor)
  # their correlation enters the implied covariances only through the
  # covariances of the items of one with the items of the other, which no
  # case observes, so it does not enter the likelihood and is not
  # identified. lavaan leaves it at its starting value of 0 under FIML,
  # cannot invert the information, and reports every defined quantity
  # (omega, AVE, H) with a standard error of 0; under two-stage estimation
  # its EM ignores the pairs with zero coverage and the second stage fits
  # the correlation to an arbitrary covariance, with a standard error of
  # 0. A third factor observed with both does not identify it. The fit is
  # refused by name, as cov_fiml() and the reliability family refuse a pair
  # of variables never observed together (.stop_never_observed_together(),
  # R/cov_fiml.R); a pair of factors with some of their items observed
  # together is identified and fitted. With correlated_factors = FALSE the
  # correlation of every pair is fixed at 0, not estimated, so nothing is
  # left unidentified and the fit is not refused: the implied covariance
  # matrix is then block diagonal by factor, the likelihood of every case
  # separates into the parts of the factors' items, and each factor's
  # estimates, standard errors, omega, AVE, and H are those of its own one
  # factor fit (the cross covariances no case observes enter neither).
  # Through the first version of this refusal such a fit was refused too,
  # with a message that called a correlation the call does not estimate
  # unidentified.
  if (length(factor_names) > 1L && correlated_factors && input_is_raw_data &&
      (missing %in% c("ml", "fiml") || two_stage)) {
    unlinked <- .cfa_k_unlinked_factors(raw_data, factors)
    if (length(unlinked)) {
      .dmar_refuse_argument(c(
        "No case observes an item of ",
        paste(unlinked, collapse = ", or an item of "),
        ", so the correlation of ",
        if (length(unlinked) > 1L) "each such pair" else "the two factors",
        " is not identified: under the factor model it enters the ",
        "likelihood only through the covariances of their items, which ",
        "the data never observe, and lavaan would leave it at its ",
        "starting value with standard errors of zero or none. Fit each ",
        "factor alone (its omega, AVE, and H are identified), or link the ",
        "forms with items given on both; missing_coverage() reports the ",
        "cases behind every pair of items."), arg = NULL)
    }
  }
  # fa_second_order() passes correlated_factors = FALSE, but its
  # first-order factors are correlated, through their second-order
  # loadings: there a pair no case observes together is refused when the
  # pairs that are observed do not identify those loadings
  # (.fa_second_order_unlinked()).
  if (second && input_is_raw_data && missing %in% c("ml", "fiml")) {
    .fa_second_order_unlinked(raw_data, factors, second_order, eq_gamma)
  }

  # ---- fit -------------------------------------------------------------
  fit_args <- c(list(model = model_syntax, estimator = estimator, se = se,
                     meanstructure = has_means),
                lavaan_args, dots)
  if (input_is_raw_data) fit_args$missing <- missing
  if (two_stage) {
    # lavaan's spelling, and its two-stage standard errors requested by
    # name, so that its notice of overriding 'se' never fires.
    fit_args$missing <- ts_missing
    fit_args$se <- ts_missing
  }
  if (!is.null(ordered_items)) {
    fit_args$ordered <- ordered_items
    # The theta parameterization keeps the residual variances of the
    # ordered items as free parameters, so the labeled psi terms and
    # every defined measurement quantity built from them remain valid;
    # under the delta parameterization those variances are not free
    # and lavaan rejects the defined parameters.
    fit_args$parameterization <- "theta"
  }
  # Under two-stage estimation with auxiliary variables lavaan warns that a
  # case observed on an auxiliary variable alone is empty and "will be
  # ignored". Its first stage reads such a case, and the fit is computed
  # again with it counted below (.ts_all_rows()), so the warning is wrong
  # about it; it is muffled when no case is empty of the auxiliary
  # variables as well, and .ts_aux_only_note() says what is done instead.
  ts_aux <- if (two_stage) intersect(.ts_dots_aux(dots), names(raw_data))
  ts_quiet_empty <- length(ts_aux) > 0L &&
    !any(rowSums(!is.na(raw_data[, c(items, ts_aux), drop = FALSE])) == 0L)
  # The warnings of a second-order model's first fit are held until it is
  # known whether that fit is the one reported (.fa_second_order_maximum()
  # below), and given then.
  first_fit <- new.env(parent = emptyenv())
  first_fit$warnings <- list()
  fit <- try(withCallingHandlers(
    do.call(lavaan::lavaan, fit_args),
    warning = function(cond) {
      if (ts_quiet_empty &&
          grepl("some cases are empty and will be ignored",
                conditionMessage(cond), fixed = TRUE))
        invokeRestart("muffleWarning")
      if (second) {
        first_fit$warnings <- c(first_fit$warnings, list(cond))
        invokeRestart("muffleWarning")
      }
    }), silent = TRUE)
  # The two failures carry classes, and lavaan's error its own message, so
  # a caller can tell them apart without reading the text: the reliability
  # family quotes lavaan's error and reports a fit that did not converge
  # as such, with its advice for an incomplete data method
  # (.omega_fit_cfa() in R/reliability_internals.R). The messages are the
  # ones these errors always had, but for the second-order model's, which
  # says what its check found.
  if (inherits(fit, "try-error")) {
    for (w in first_fit$warnings) warning(w)
    lavaan_message <- conditionMessage(attr(fit, "condition"))
    .dmar_refuse(
      paste0("Model fitting with lavaan failed: ",
        lavaan_message),
      "lavaan_error", lavaan_message = lavaan_message)
  }
  # The fit the table's parameter and index rows are read from, and the fit
  # its fit rows are read from: `fit` itself, but for a second-order
  # solution on its boundary (.fa_second_order_maximum()).
  table_fit <- rows_fit <- fit
  so_max <- NULL
  if (second) {
    # The maximum of the likelihood, found in a second parameterization
    # that reaches a standardized second-order loading of one, carried
    # back, and turned to the reported direction; or the plain statement
    # that it could not be found.
    so_max <- .fa_second_order_maximum(
      fit, fit_args,
      .fa_second_order_check_syntax(factors, second_order, eq_load, eq_int,
                                    eq_err, eq_gamma, has_means),
      factors, second_order, eq_gamma, first_fit$warnings)
    for (w in so_max$warnings) warning(w)
    if (!so_max$converged) {
      .dmar_refuse(
        .fa_second_order_not_converged(so_max),
        "not_converged")
    }
    fit <- so_max$fit
    table_fit <- so_max$table_fit
    rows_fit <- so_max$rows_fit
    .fa_second_order_check_notes(so_max, output, fit_args$se)
  } else if (!lavaan::lavInspect(fit, "converged")) {
    .dmar_refuse(
      "The model did not converge.",
      "not_converged")
  }
  actual_options <- lavaan::lavInspect(fit, "options")
  if (input_is_raw_data && missing %in% c("fiml", "ml", "ml.x") &&
      !actual_options$missing %in% c("ml", "ml.x")) {
    .dmar_refuse_argument(c(
      "FIML was requested, but lavaan changed missing-data handling to ",
      actual_options$missing, " for this estimator/SE/test combination. ",
      "No result is returned from the changed case set. Use a compatible ",
      "ML or MLR combination (robust.huber.white for robust FIML SEs) ",
      "or explicitly request listwise deletion."), arg = NULL)
  }
  h1_diagnostic <- .fiml_fit_h1_diagnostic(fit)
  if (two_stage && !is.null(h1_diagnostic)) {
    .dmar_refuse_argument(c(
      "Two-stage estimation requires unrestricted first-stage moments. ",
      h1_diagnostic), arg = NULL)
  }
  if (!is.null(h1_diagnostic) && !.fiml_h0_vcov_supported(fit)) {
    .dmar_refuse(c(
      "The unrestricted H1 diagnostic failed, and this custom SE/information ",
      "combination has not been verified to use only the constrained H0 ",
      "model. Refit with standard observed-Hessian SEs or the default ",
      "MLR Huber-White SEs before reporting an interval. ",
      h1_diagnostic), "computation_failed")
  }
  # ---- a likelihood with no maximum ------------------------------------
  # Before anything is read from the fit, and for every output (output =
  # "fit" too, the path the reliability family consumes): when an error
  # variance is negative, whether the likelihood has a maximum at all, or
  # the fit sits on a ridge where that variance runs off below zero, the
  # likelihood flat or still rising, and lavaan's optimizer stopped where
  # its steps stopped gaining (.lavaan_no_maximum(),
  # R/no_maximum_internals.R). Such a fit is refused, with the lavaan code
  # that reproduces it; a finite negative estimate is not, and goes on to
  # the Heywood warning below. lavaan reads the model's variables from a
  # data frame and, by name, from a covariance matrix, but it pairs the
  # means with the rows of the covariance matrix by position, so with
  # means named otherwise than S the code selects the items from both
  # (.no_max_call_text()). A second-order fit held at its boundary is
  # checked in the parameterization that reaches it (`rows_fit`, a fit at
  # its maximum; see .fa_second_order_maximum()), and the code printed is
  # that fit's syntax and arguments: the first fit's would stop elsewhere.
  # From lavaan's own starting values its optimizer can stop at another
  # point of the same ridge, and the message says that the code fits the
  # model as DMAR checked it. Least squares and WLSMV fits have no
  # likelihood, and the check does not apply to them.
  at_bound <- length(so_max$at_one) > 0L
  no_max <- .lavaan_no_maximum(if (at_bound) rows_fit else fit)
  if (!is.null(no_max)) {
    select_input <- if (!is.null(lavaan_args$sample.mean) &&
                        !identical(names(M), S_names)) {
      list(sample.cov = items, sample.mean = items)
    }
    if (at_bound) {
      .no_max_refuse(fn, no_max, syntax = so_max$rows_model,
                     lavaan_call = .no_max_call_text(
                       "lavaan::lavaan", so_max$rows_args,
                       select = select_input),
                     code_intro = paste0(
                       "This lavaan code fits the model as DMAR checked ",
                       "it, with the variance of each first-order ",
                       "factor's disturbance bounded at zero so that a ",
                       "standardized second-order loading of one can be ",
                       "reached (from lavaan's own starting values its ",
                       "optimizer can stop at another point of the same ",
                       "ridge)"))
    }
    .no_max_refuse(fn, no_max, syntax = model_syntax,
                   lavaan_call = .no_max_call_text("lavaan::lavaan",
                                                   fit_args,
                                                   select = select_input))
  }
  # A case observed on an auxiliary variable but on no item is read by the
  # first stage and left out of lavaan's N, so the fit is computed again
  # with its N set to every case the first stage read (.ts_all_rows()):
  # the estimates stay (the model is fit in one group, where they do not
  # depend on N), and the standard errors, the tests, the indices, and the
  # criteria of every output, the lavaan fit itself included, use all the
  # cases. The message says so once that is done; when it could not be,
  # the warning alone says what rests on lavaan's N. Through the first
  # versions of 1.1.0 such cases were set aside, which conditioned the
  # estimates on an observed item and made them inconsistent when the
  # items went missing together at random given the auxiliary variables.
  ts_rows <- NULL
  if (two_stage) {
    ts_rows <- .ts_all_rows(fit)
    fit <- ts_rows$fit
    # The table's rows are read from the fit computed again. table_fit and
    # rows_fit were set to the first fit above; the layered models, which
    # set them otherwise, take no two-stage estimation (missing_choices).
    table_fit <- rows_fit <- fit
    if (isTRUE(ts_rows$corrected)) {
      .ts_aux_only_note(fn, ts_rows$extra, sum(ts_rows$read),
                        sum(ts_rows$nobs))
    }
    .ts_all_rows_warning(fn, ts_rows, ts_missing)
  }

  # The Heywood check runs before any return so it also covers
  # output = "fit", the path the reliability family consumes; the classed
  # condition is counted across bootstrap replications there, so the
  # class, not just the text, is part of the contract (as in the
  # dmar_heywood_warning cfa_1() signaled before it became a wrapper).
  pe_hey <- lavaan::parameterEstimates(fit)
  hey <- pe_hey$op == "~~" & pe_hey$lhs == pe_hey$rhs &
    pe_hey$lhs %in% items & pe_hey$est < 0
  if (any(hey)) {
    .dmar_warn(
      paste0(
        "Negative error variance estimate (Heywood case) for: ",
        paste(unique(pe_hey$lhs[hey]), collapse = ", "),
        ". Interpret the solution with care."),
      "heywood_warning")
  }

  if (output == "fit") return(.relia_h1_qualify(fit, h1_diagnostic))
  if (output == "summary") {
    if (!is.null(h1_diagnostic)) {
      .dmar_refuse_argument(c(
        "The lavaan summary includes unrestricted-H1 fit diagnostics ",
        "that are unavailable for these data. Use output = \"verbose\" ",
        "for the eligible constrained-model estimates and intervals. ",
        h1_diagnostic), arg = NULL)
    }
    if (second && identical(lavaan::lavInspect(fit, "options")$optim.method,
                            "none")) {
      return(.fa_second_order_summary(fit, length(so_max$at_one) > 0L))
    }
    return(lavaan::summary(fit, fit.measures = TRUE))
  }
  # The rows whose limits are formed on the logit scale: omega, AVE, and H
  # of each factor, or the bifactor indices that lavaan defines
  # (.fa_bifactor_index_terms()), for the second-order model those of each
  # second-order factor and its first-order factors.
  logit_terms <- if (bifactor) {
    .fa_bifactor_index_terms(general, factor_names)
  } else if (second) {
    .fa_second_order_index_terms(second_order)
  } else {
    as.vector(outer(c("omega_", "ave_", "H_"), factor_names, paste0))
  }
  # Whether the covariance matrix of the estimates can be used, decided
  # once, with the information matrix tested for singularity whether it
  # is observed (FIML) or expected (listwise deletion, a covariance
  # matrix), and with the layered models' maximal reliability rows
  # following the same decision as every other row
  # (.cfa_k_vcov_status()). Through 1.1.0 development fa_confirmatory()
  # tested the observed information only, and when lavaan inverted a
  # singular expected information its tables reported finite but
  # meaningless standard errors and intervals (an H of .853 with a
  # standard error of .028 on the listwise fit of the fixture of the
  # test of this guard, with equal loadings on its other factor). (For
  # a second-order solution on its boundary the parameter rows are read
  # from `table_fit`, whose held loadings are constrained; see
  # .fa_second_order_maximum(). Otherwise `table_fit` is `fit`.)
  vcov_status <- .cfa_k_vcov_status(table_fit)
  no_vcov <- function(table, cols) {
    .cfa_k_no_vcov(table_fit, table, cols, status = vcov_status)
  }
  if (output == "standardized") {
    # lavaan evaluates a defined parameter at the standardized values of
    # its labels, and a label shared by several items (equal_loading,
    # equal_error) has one unstandardized value but as many standardized
    # values as the items have different variances: lavaan takes the
    # first item's and warns that the result is ambiguous. In the models
    # of fa_confirmatory() and fa_bifactor() the defined rows that use
    # such a label are computed again here from each item's own
    # standardized loadings and error variance (.fa_cfa_std_defined(),
    # .fa_bifactor_std_defined()), and the second-order model's defined
    # rows are left out (below), so lavaan's warning concerns values that
    # are not reported and is muffled (that warning alone). Through 1.1.0
    # development, as in cfa_k() of DMAR 1.0.0, the rows of
    # fa_confirmatory() and fa_bifactor() held lavaan's values, and the
    # warning reached the caller of all three.
    shared <- if (!layered) {
      factor_names[eq_load[factor_names] | eq_err[factor_names]]
    } else {
      character(0)
    }
    std <- withCallingHandlers(
      lavaan::standardizedSolution(table_fit, level = conf_level),
      warning = function(w) {
        if ((layered || length(shared)) && .fa_cfa_std_ambiguous(w)) {
          invokeRestart("muffleWarning")
        }
      })
    std <- if (bifactor) {
      .fa_bifactor_std_defined(table_fit, std, factors, general,
                               general_only, conf_level)
    } else if (!second) {
      .fa_cfa_std_defined(table_fit, std, factors[shared], conf_level)
    } else {
      std
    }
    # Under se = "none" lavaan's standardized solution has no standard
    # error, test, or interval columns; they are NA here, as in the other
    # tables, so that the steps below (and a caller's code) find them.
    # Through 1.1.0 development the logit-scale limits stopped such a
    # call with an R error.
    for (.col in c("se", "z", "pvalue", "ci.lower", "ci.upper")) {
      if (is.null(std[[.col]])) std[[.col]] <- NA_real_
    }
    # The second-order model's indices are written with each
    # first-order loading as the item's loading on the disturbance of its
    # factor, which holds while that disturbance's variance is 1, and the
    # standardized loading is standardized with the factor's whole
    # variance instead, so those values describe nothing; the rows are
    # left out (the standardized Schmid-Leiman loadings are the
    # schmid_leiman attribute of the default table).
    # The constraints that hold a loading at the boundary are not rows of
    # the model either.
    if (second) {
      std <- std[!std$op %in% c(":=", "=="), , drop = FALSE]
      rownames(std) <- NULL
    }
    std <- no_vcov(std, intersect(c("se", "z", "pvalue", "ci.lower",
                                    "ci.upper"), names(std)))
    df_std <- rep(Inf, nrow(std))
    if (!.wald_df_is_normal(wald_df)) {
      df_std <- unname(.wald_df_resolve(
        wald_df, fn, fit = table_fit, rows = std,
        context = list(N = sum(unlist(lavaan::lavInspect(table_fit, "nobs"))),
                       factors = factors,
                       requested_missing = requested_missing))$df)
      std <- .wald_df_apply(std, df_std, conf_level, est = "est.std",
                            se = "se", z = "z", p = "pvalue",
                            lower = "ci.lower", upper = "ci.upper")
    }
    # Omega, AVE, and H (or the bifactor indices) take logit-scale limits
    # here too (see the verbose table).
    std$term <- std$lhs
    std <- .cfa_k_logit_limits(std, logit_terms, conf_level, df_std,
                               est = "est.std", lower = "ci.lower",
                               upper = "ci.upper")
    std$term <- NULL
    attr(std, "wald_df") <- wald_df
    return(.relia_h1_qualify(std, h1_diagnostic))
  }

  # ---- name the implied classical structure, per factor ----------------
  # (The bifactor and second-order models have no classical name; their
  # tables describe the model in the header their print methods write,
  # .dmar_tbl_footer.dmar_fa_bifactor() in R/fa_bifactor.R and
  # .dmar_tbl_footer.dmar_fa_second_order() in R/fa_second_order.R.)
  model_name <- if (!layered) vapply(factor_names, function(f) {
    .cfa_k_model_name(eq_load[f], if (has_means) eq_int[f] else FALSE,
                      eq_err[f], has_means)
  }, character(1L))

  # ---- parameter table with confidence intervals -----------------------
  pe <- as.data.frame(lavaan::parameterEstimates(table_fit,
                                                 level = conf_level))
  # The constraints that hold a second-order loading at the boundary
  # (.fa_second_order_maximum()) are not rows of the model.
  if (second) pe <- pe[pe$op != "==", , drop = FALSE]
  for (.col in c("se", "z", "pvalue", "ci.lower", "ci.upper")) {
    if (is.null(pe[[.col]])) pe[[.col]] <- NA_real_
  }
  pe <- no_vcov(pe, c("se", "z", "pvalue", "ci.lower", "ci.upper"))
  pe$syntax <- paste(pe$lhs, pe$op, pe$rhs)
  defined_rows <- pe$op == ":="
  pe$label[defined_rows] <- pe$lhs[defined_rows]
  # The full := expressions (ave, H) are long enough to swamp the table;
  # the term names them, so the syntax column stays blank for those rows.
  pe$syntax[defined_rows] <- ""
  latent_names <- if (bifactor) {
    c(general, factor_names)
  } else if (second) {
    c(factor_names, names(second_order))
  } else {
    factor_names
  }
  fv <- pe$op == "~~" & pe$lhs == pe$rhs & pe$lhs %in% latent_names
  pe$label[fv] <- paste0("phi_", pe$lhs[fv])
  # In the second-order model the variance of a first-order factor's
  # disturbance, fixed at 1, is its zeta row; the phi rows are the
  # second-order factors' variances (fixed at 1) and correlations.
  if (second) {
    zeta <- fv & pe$lhs %in% factor_names
    pe$label[zeta] <- paste0("zeta_", pe$lhs[zeta])
  }
  # The factor means are fixed to 0 for identification and carry no
  # information; drop those rows.
  pe <- pe[!(pe$op == "~1" & pe$lhs %in% latent_names), , drop = FALSE]
  # In the bifactor model every factor covariance is fixed at 0, which is
  # the model itself (the header says so), so those rows are dropped too.
  if (bifactor) {
    pe <- pe[!(pe$op == "~~" & pe$lhs != pe$rhs & pe$lhs %in% latent_names &
                 pe$rhs %in% latent_names), , drop = FALSE]
  }

  params <- data.frame(
    syntax   = pe$syntax,
    term     = pe$label,
    estimate = pe$est,
    se       = pe$se,
    test_statistic  = pe$z,
    p_value  = pe$pvalue,
    ci_lower = pe$ci.lower,
    ci_upper = pe$ci.upper,
    stringsAsFactors = FALSE
  )
  # With a finite wald_df the intervals and p-values are recomputed from
  # the same estimates and standard errors on the t reference; at the
  # default they are lavaan's own, untouched. Under the keyword only the
  # loadings take the t reference, each on max(3, N(1 - FMI) - 1) degrees
  # of freedom from its own fraction of missing information; the other
  # rows get Inf, which .wald_df_apply() leaves as lavaan computed them.
  # Omega passed the calibration on its own scale, but its limits are now
  # formed on the logit scale, where the t reference overcovered (.9645),
  # so it keeps the normal reference with AVE and H.
  wald_table <- NULL
  row_df <- rep(Inf, nrow(params))
  omega_terms <- character(0)
  if (!.wald_df_is_normal(wald_df)) {
    wald <- .wald_df_resolve(
      wald_df, fn, fit = table_fit, rows = pe,
      context = list(N = sum(unlist(lavaan::lavInspect(table_fit, "nobs"))),
                     factors = factors, terms = pe$label,
                     omega = omega_terms,
                     requested_missing = requested_missing))
    params <- .wald_df_apply(params, unname(wald$df), conf_level)
    wald_table <- wald$table
    loading_terms <- pe$label[pe$op == "=~"]
    row_df <- unname(wald$df)
  }
  # The omega, AVE, and H rows: Wald limits on the logit scale from
  # lavaan's delta method standard error, transformed back (Browne, 1982;
  # Kelley, 2007), as reliability_H() and fa_exploratory() form H, so they
  # lie within (0, 1); with a finite wald_df the t quantile replaces the
  # normal one (under the missing-information keyword these rows keep the
  # normal reference, as every row but the loadings does). An estimate at
  # or beyond one (an improper solution) gets NA limits here and is set to
  # NA below; categorical omega replaces its row below.
  params <- .cfa_k_logit_limits(params, logit_terms, conf_level, row_df)
  # What the rows of the table x that took the t reference are, for the
  # line under it.
  wald_rows <- function(x) {
    if (is.null(wald_table)) return(NULL)
    c(if (any(x$term %in% loading_terms)) "loadings",
      if (any(x$term %in% omega_terms)) "omega")
  }

  # ---- the layered models' indices and table ---------------------------
  # The maximal reliabilities and PUC join the indices lavaan defined, an
  # improper solution leaves the rows it makes undefined NA
  # (.fa_bifactor_rows(), .fa_second_order_rows()), and the table is the
  # indices alone ("measurement") or the parameters, the indices, and the
  # fit rows.
  if (layered) {
    result_class <- if (bifactor) "dmar_fa_bifactor" else
      "dmar_fa_second_order"
    vcov_usable <- identical(vcov_status, "usable")
    bf <- if (bifactor) {
      .fa_bifactor_rows(params, fit, pe_hey, factors, general, items,
                        conf_level, wald_df, vcov_usable)
    } else {
      .fa_second_order_rows(params, table_fit, pe_hey, factors,
                            second_order, items, conf_level, wald_df,
                            vcov_usable, at_one = so_max$at_one)
    }
    params <- bf$params
    if (output == "measurement") {
      out <- params[params$term %in% bf$index_terms, , drop = FALSE]
      rownames(out) <- NULL
      class(out) <- c(result_class, "data.frame")
      out <- .as_dmar_tbl(out, conf_level = conf_level)
    } else {
      fr <- .fa_cfa_fit_rows(rows_fit, h1_diagnostic, rmsea_conf_level, has_means,
                             ts_missing, two_stage, ordered_items)
      out <- rbind(params, fr$rows)
      rownames(out) <- NULL
      class(out) <- c(result_class, "data.frame")
      out <- .as_dmar_tbl(out, conf_level = conf_level,
                          p_terms = "p_chi_square",
                          fixed_terms = c("AIC", "BIC", "H0", "H1"))
      out <- .fa_cfa_fit_attributes(
        out, fr, estimator = estimator, missing = missing,
        input_is_raw_data = input_is_raw_data, two_stage = two_stage,
        ts_missing = ts_missing, ts_rows = ts_rows, ts_aux = ts_aux,
        complete = !input_is_raw_data ||
          !anyNA(raw_data[, items, drop = FALSE]),
        h1_diagnostic = h1_diagnostic)
    }
    out <- if (bifactor) {
      .fa_bifactor_attributes(out, factors, general, general_only, eq_load,
                              eq_err, eq_int, has_means, bf$improper)
    } else {
      .fa_second_order_attributes(out, factors, second_order, eq_load,
                                  eq_err, eq_int, eq_gamma, has_means,
                                  bf$improper, bf$improper_second_order,
                                  bf$schmid_leiman)
    }
    out <- .wald_df_attach(out, wald_df, wald_table, wald_rows(out))
    return(.relia_h1_qualify(out, h1_diagnostic))
  }

  # ---- improper solutions ----------------------------------------------
  # The items of each factor whose error variance is negative, with their
  # communalities (.cfa_k_improper()); the rows they leave undefined are
  # set to NA below, after the categorical omegas.
  improper <- .cfa_k_improper(pe_hey, factors)

  # ---- categorical sum score omega for ordered factors -----------------
  # The defined-parameter omega is on the latent response metric, which
  # is not the metric of a sum of observed categories. For a factor
  # whose items are ordered, replace it with categorical omega (the
  # Green and Yang, 2009, assembly with the unrestricted polychoric
  # denominator of Kelley and Pornprasertmanit, 2016) computed from the
  # same fit; its delta method interval does not transfer, so the
  # interval columns are NA. A factor whose polychorics the data do not
  # all determine gets NA rather than a total variance built on an
  # undetermined correlation. An improper solution for one of the
  # factor's items leaves its categorical omega undefined, since the
  # assembly evaluates the bivariate normal distribution at each
  # communality: the factor gets NA, and the warning of the improper pass
  # below names the items, rather than the whole fit stopping (through
  # 1.0.1 it stopped with mvtnorm's "'corr' is not a correlation
  # matrix").
  if (any(ordered_factor)) {
    for (f in factor_names[ordered_factor]) {
      gy <- if (!is.null(improper[[f]])) {
        NA_real_
      } else {
        undetermined <- .undetermined_polychorics(fit, factors[[f]])
        if (nrow(undetermined) > 0L) {
          .dmar_warn(c(
            "omega_", f, " is NA: the data do not determine ",
            .undetermined_message(undetermined), " (the bivariate ",
            "table has empty cells), which the total variance of ",
            "categorical omega needs. ",
            "reliability_omega_categorical(denominator = ",
            "\"model_implied\") on the factor's items gives the ",
            "coefficient with the model implied total variance when ",
            "the correlations the data do determine identify the ",
            "single-factor loadings, and says when they do not."),
            "result_qualified")
          NA_real_
        } else {
          .omega_categorical_from_kfit(fit, factors[[f]])
        }
      }
      # The assembly's own check, a communality above one, agrees with
      # the negative residual variances found above; it stays as a
      # backstop.
      if (identical(attr(gy, "dmar_reason"), "improper")) {
        improper[[f]] <- attr(gy, "improper")
        gy <- NA_real_
      }
      row_f <- params$term == paste0("omega_", f)
      params$estimate[row_f] <- gy
      params$se[row_f]       <- NA_real_
      params$test_statistic[row_f]  <- NA_real_
      params$p_value[row_f]  <- NA_real_
      params$ci_lower[row_f] <- NA_real_
      params$ci_upper[row_f] <- NA_real_
    }
    message("omega_",
            paste(factor_names[ordered_factor], collapse = ", omega_"),
            " reported on the categorical sum score metric (categorical ",
            "omega; Green & Yang, 2009; Kelley & Pornprasertmanit, 2016).")
  }
  # An item with a negative error variance has a communality above one.
  # The average variance extracted averages the communalities, and
  # coefficient H, the largest reliability any weighted sum of the items
  # attains, is at least that communality (the item weighted alone
  # attains it), so for such a factor neither is a proportion of variance
  # and both rows are NA. Omega, the reliability of the unit-weighted sum,
  # depends on the error variances only through their sum: while the sum
  # is positive it is a proportion and is reported, as reliability_omega()
  # reports it, beside the Heywood warning above; when the sum is negative
  # the model implies a negative error variance for the sum score and an
  # omega above one, and the row is NA. An ordered factor's categorical
  # omega is NA already. Through 1.0.1 these rows reported, for a factor
  # of visual perception and lozenges beside a verbal factor of
  # holzinger_swineford, an H of 1.32, and on other two-item factors an
  # omega and an AVE above one.
  na_cols <- c("estimate", "se", "test_statistic", "p_value", "ci_lower",
               "ci_upper")
  for (f in factor_names) {
    if (is.null(improper[[f]])) next
    sum_negative <- !ordered_factor[[f]] && isTRUE(
      params$estimate[params$term == paste0("error_sum_", f)] < 0)
    omega_na <- ordered_factor[[f]] || sum_negative
    rows <- params$term %in% paste0(c(if (omega_na) "omega_", "ave_", "H_"),
                                    f)
    params[rows, na_cols] <- NA_real_
    .dmar_warn(
      .kfit_improper_message(f, improper[[f]],
                             ordered = ordered_factor[[f]],
                             omega = omega_na), "result_qualified")
  }
  omega_metric <- stats::setNames(
    ifelse(ordered_factor, "categorical_sum_score", "linear"),
    factor_names)

  if (output == "measurement") {
    property_terms <- as.vector(outer(c("omega_", "ave_", "H_"),
                                      factor_names, paste0))
    pair_terms <- character(0)
    if (length(factor_names) > 1L && correlated_factors) {
      pairs <- utils::combn(factor_names, 2L)
      pair_terms <- paste0("phi_", pairs[1L, ], "_", pairs[2L, ])
    }
    out <- params[params$term %in% c(property_terms, pair_terms), ,
                  drop = FALSE]
    if (input_is_raw_data && length(factor_names) > 1L) {
      # A pair whose HTMT the data leave undefined (items never observed
      # together, as in a planned missing data design under FIML) is NA,
      # with htmt()'s warning naming it; when no pair is defined, htmt()
      # stops, and here every htmt row is NA with that message as a
      # warning, so the model based rows above are kept. Through 1.0.1 an
      # undefined pair was NA without a word.
      ht <- tryCatch(
        htmt(raw_data, blocks = factors),
        dmar_htmt_undefined = function(e) {
          .dmar_warn(c(
            conditionMessage(e), " The htmt rows are NA."), "result_qualified")
          pr <- utils::combn(names(factors), 2L)
          data.frame(construct_1 = pr[1L, ], construct_2 = pr[2L, ],
                     htmt = NA_real_, stringsAsFactors = FALSE)
        })
      out <- rbind(out, data.frame(
        syntax   = paste(ht$construct_1, "~~", ht$construct_2),
        term     = paste0("htmt_", ht$construct_1, "_", ht$construct_2),
        estimate = ht$htmt,
        se = NA_real_, test_statistic = NA_real_, p_value = NA_real_,
        ci_lower = NA_real_, ci_upper = NA_real_,
        stringsAsFactors = FALSE))
    }
    rownames(out) <- NULL
    # dmar_cfa_k, the class of these results in DMAR 1.0.0 (when the
    # function was cfa_k()), trails so that code testing it keeps working.
    class(out) <- c("dmar_fa_confirmatory", "dmar_cfa_k", "data.frame")
    out <- .as_dmar_tbl(out, conf_level = conf_level)
    attr(out, "model") <- model_name
    attr(out, "omega_metric") <- omega_metric
    out <- .wald_df_attach(out, wald_df, wald_table, wald_rows(out))
    if (two_stage) {
      attr(out, "two_stage") <- list(
        missing = missing, lavaan_missing = ts_missing,
        rows = ts_rows[c("read", "nobs", "extra", "corrected")],
        footer = .ts_footer(ts_missing, "measurement", rows = ts_rows))
    }
    return(.relia_h1_qualify(out, h1_diagnostic))
  }

  # ---- fit information rows --------------------------------------------
  fr <- .fa_cfa_fit_rows(fit, h1_diagnostic, rmsea_conf_level, has_means,
                         ts_missing, two_stage, ordered_items)
  out <- rbind(params, fr$rows)
  rownames(out) <- NULL
  class(out) <- c("dmar_fa_confirmatory", "dmar_cfa_k", "data.frame")
  out <- .as_dmar_tbl(out, conf_level = conf_level,
                      p_terms = "p_chi_square",
                      fixed_terms = c("AIC", "BIC", "H0", "H1"))
  attr(out, "model") <- model_name
  attr(out, "factors") <- factors
  attr(out, "has_means") <- has_means
  attr(out, "omega_metric") <- omega_metric
  out <- .wald_df_attach(out, wald_df, wald_table, wald_rows(out))
  out <- .fa_cfa_fit_attributes(
    out, fr, estimator = estimator, missing = missing,
    input_is_raw_data = input_is_raw_data, two_stage = two_stage,
    ts_missing = ts_missing, ts_rows = ts_rows, ts_aux = ts_aux,
    complete = !input_is_raw_data ||
      !anyNA(raw_data[, items, drop = FALSE]),
    h1_diagnostic = h1_diagnostic)
  .relia_h1_qualify(out, h1_diagnostic)
}


# .fa_cfa_correlated_syntax(factors, eq_load, eq_int, eq_err, has_means,
#                           correlated_factors)
#
# The lavaan syntax of fa_confirmatory()'s model: each factor measured by
# its own items (simple structure), its variance fixed at 1 (and its mean
# at 0 with a mean structure) and its first loading freed, every loading,
# error variance, and intercept labeled (equated within a factor when its
# flag says so), the factor correlations labeled phi_<f1>_<f2> or fixed at
# 0, and the defined parameters behind the measurement rows. Returns one
# string.
.fa_cfa_correlated_syntax <- function(factors, eq_load, eq_int, eq_err,
                                      has_means, correlated_factors) {
  factor_names <- names(factors)
  lines <- character(0)
  lambda_labels <- list()
  psi_labels    <- list()
  for (f in factor_names) {
    v <- factors[[f]]
    J <- length(v)
    lam <- if (eq_load[f]) rep(paste0("lambda_", f), J) else
      paste0("lambda_", f, "_", seq_len(J))
    psi <- if (eq_err[f]) rep(paste0("psi_", f), J) else
      paste0("psi_", f, "_", seq_len(J))
    lambda_labels[[f]] <- lam
    psi_labels[[f]]    <- psi
    lines <- c(lines,
               paste0(f, " =~ NA*", v[1L], " + ",
                      paste0(lam, "*", v, collapse = " + ")),
               paste0(f, " ~~ 1*", f),
               paste0(v, " ~~ ", psi, "*", v))
    if (has_means) {
      nu <- if (eq_int[f]) rep(paste0("nu_", f), J) else
        paste0("nu_", f, "_", seq_len(J))
      lines <- c(lines, paste0(v, " ~ ", nu, "*1"), paste0(f, " ~ 0*1"))
    }
  }
  if (length(factor_names) > 1L) {
    pairs <- utils::combn(factor_names, 2L)
    for (p in seq_len(ncol(pairs))) {
      f1 <- pairs[1L, p]; f2 <- pairs[2L, p]
      lines <- c(lines, if (correlated_factors) {
        paste0(f1, " ~~ phi_", f1, "_", f2, "*", f2)
      } else {
        paste0(f1, " ~~ 0*", f2)
      })
    }
  }
  # Defined parameters: per-factor loading and error sums, omega, ave,
  # and H. lavaan's delta method then provides a standard error and
  # confidence interval for each at no extra cost to the fit.
  for (f in factor_names) {
    lam <- lambda_labels[[f]]
    psi <- psi_labels[[f]]
    J <- length(lam)
    ave_terms <- paste0("(", lam, "^2/(", lam, "^2 + ", psi, "))")
    h_terms   <- paste0("(", lam, "^2/", psi, ")")
    lines <- c(lines,
      paste0("loading_sum_", f, " := ", paste(lam, collapse = " + ")),
      paste0("error_sum_", f, " := ", paste(psi, collapse = " + ")),
      paste0("omega_", f, " := (loading_sum_", f, "^2) / ",
             "((loading_sum_", f, "^2) + error_sum_", f, ")"),
      paste0("ave_", f, " := (", paste(ave_terms, collapse = " + "),
             ") / ", J),
      paste0("H_", f, " := 1 / (1 + 1/(",
             paste(h_terms, collapse = " + "), "))"))
  }
  paste(lines, collapse = "\n")
}

# .fa_cfa_std_ambiguous(w)
#
# TRUE for lavaan's warning that the standardized value of a defined
# parameter is ambiguous because it uses a label shared by several
# coefficients (lav_partable_check_std_def()). It is matched on the name
# of the function lavaan prefixes to the message, which is not
# translated, or on the English text.
.fa_cfa_std_ambiguous <- function(w) {
  msg <- conditionMessage(w)
  grepl("lav_partable_check_std_def", msg, fixed = TRUE) ||
    grepl("standardized value of the defined parameter", msg, fixed = TRUE)
}

# .fa_cfa_std_indices(est, rows, grad = TRUE)
#
# The defined quantities of one factor of fa_confirmatory()'s model
# (.fa_cfa_correlated_syntax()) on the standardized scale, each item
# standardized by its own model implied variance: with loading a_j, error
# variance p_j, and factor variance phi (the parameter table values `est`
# at the rows `rows$lambda`, `rows$psi`, and `rows$phi`), the item's
# variance is v_j = a_j^2 phi + p_j (for an ordered item, that of its
# latent response under the theta parameterization), its standardized
# loading s_j = a_j sqrt(phi) / sqrt(v_j), and its standardized error
# variance e_j = p_j / v_j, the values lavaan's standardizedSolution()
# reports in each item's own rows. The quantities are those the syntax
# defines, evaluated at s_j and e_j:
#
#   loading_sum = sum s_j,      error_sum = sum e_j,
#   omega = loading_sum^2 / (loading_sum^2 + error_sum),
#   ave   = J^{-1} sum s_j^2 / (s_j^2 + e_j),
#   H     = 1 / (1 + 1 / sum s_j^2 / e_j).
#
# Omega here is the reliability of the unit-weighted sum of the
# standardized items; AVE and H do not depend on the metric of the items
# and equal the rows of the unstandardized table.
#
# Returns a list: `value`, the five quantities (named loading_sum,
# error_sum, omega, ave, H), and, with `grad`, `grad`, the 5 x (2J + 1)
# matrix of their derivatives with respect to (a_1, ..., a_J, p_1, ...,
# p_J, phi), in closed form by the chain rule:
#
#   ds_j/da_j = sqrt(phi) p_j / v_j^{3/2},
#   ds_j/dp_j = -a_j sqrt(phi) / (2 v_j^{3/2}),
#   ds_j/dphi = a_j p_j / (2 sqrt(phi) v_j^{3/2}),
#   de_j/da_j = -2 a_j phi p_j / v_j^2,
#   de_j/dp_j = a_j^2 phi / v_j^2,
#   de_j/dphi = -a_j^2 p_j / v_j^2.
.fa_cfa_std_indices <- function(est, rows, grad = TRUE) {
  a <- est[rows$lambda]
  p <- est[rows$psi]
  phi <- est[rows$phi]
  J <- length(a)
  v <- a^2 * phi + p
  s <- a * sqrt(phi) / sqrt(v)
  e <- p / v
  L <- sum(s)
  E <- sum(e)
  r <- s^2 / (s^2 + e)
  Q <- sum(s^2 / e)
  value <- c(loading_sum = L, error_sum = E, omega = L^2 / (L^2 + E),
             ave = mean(r), H = 1 / (1 + 1 / Q))
  if (!grad) return(list(value = value))
  # Derivatives of each quantity with respect to s_j (row 1) and e_j
  # (row 2), one column per item.
  d_se <- list(
    loading_sum = rbind(rep(1, J), rep(0, J)),
    error_sum = rbind(rep(0, J), rep(1, J)),
    omega = rbind(rep(2 * L * E / (L^2 + E)^2, J),
                  rep(-L^2 / (L^2 + E)^2, J)),
    ave = rbind(2 * s * e / (s^2 + e)^2, -s^2 / (s^2 + e)^2) / J,
    H = rbind(2 * s / e, -s^2 / e^2) / (1 + Q)^2)
  ds <- rbind(a = sqrt(phi) * p / v^1.5, p = -a * sqrt(phi) / (2 * v^1.5),
              phi = a * p / (2 * sqrt(phi) * v^1.5))
  de <- rbind(a = -2 * a * phi * p / v^2, p = a^2 * phi / v^2,
              phi = -a^2 * p / v^2)
  G <- t(vapply(d_se, function(d) {
    by_item <- ds * rep(d[1L, ], each = 3L) + de * rep(d[2L, ], each = 3L)
    c(by_item["a", ], by_item["p", ], sum(by_item["phi", ]))
  }, numeric(2L * J + 1L)))
  list(value = value, grad = G)
}

# .fa_cfa_std_defined(fit, std, factors, conf_level)
#
# lavaan's standardized solution `std` of a fa_confirmatory() fit with the
# rows loading_sum_<f>, error_sum_<f>, omega_<f>, ave_<f>, and H_<f> of
# each factor of `factors` (those whose loadings or error variances share
# a label) computed again from the items' own standardized values
# (.fa_cfa_std_indices()), with their standard errors and limits
# (.fa_std_redefine()). Returns `std`.
.fa_cfa_std_defined <- function(fit, std, factors, conf_level) {
  if (!length(factors)) return(std)
  pt <- lavaan::parTable(fit)
  key <- paste(pt$lhs, pt$op, pt$rhs)
  blocks <- lapply(names(factors), function(f) {
    items <- factors[[f]]
    rows <- list(lambda = match(paste(f, "=~", items), key),
                 psi = match(paste(items, "~~", items), key),
                 phi = match(paste(f, "~~", f), key))
    cols <- c(rows$lambda, rows$psi, rows$phi)
    list(terms = paste0(c("loading_sum_", "error_sum_", "omega_", "ave_",
                          "H_"), f),
         fun = function(est, grad = TRUE) {
           at <- .fa_cfa_std_indices(est, rows, grad)
           if (grad) {
             G <- matrix(0, nrow(at$grad), length(est))
             for (k in seq_along(cols)) {
               G[, cols[k]] <- G[, cols[k]] + at$grad[, k]
             }
             at$grad <- G
           }
           at
         })
  })
  .fa_std_redefine(fit, std, blocks, conf_level, pt)
}

# .fa_std_shared_terms(pt)
#
# The defined parameters of the parameter table `pt` whose expressions
# use a label that more than one parameter carries, the rows whose
# standardized values lavaan's standardizedSolution() takes from the
# first of those parameters (the rows its lav_partable_check_std_def()
# names). Returns their names.
.fa_std_shared_terms <- function(pt) {
  def <- which(pt$op == ":=")
  lab <- pt$label[!pt$op %in% c(":=", "==", "<", ">") & nzchar(pt$label)]
  shared <- unique(lab[duplicated(lab)])
  if (!length(def) || !length(shared)) return(character(0))
  uses <- vapply(def, function(k) {
    any(all.vars(str2lang(pt$rhs[k])) %in% shared)
  }, logical(1L))
  pt$lhs[def[uses]]
}

# .fa_bifactor_std_defined(fit, std, factors, general, general_only,
#                          conf_level)
#
# lavaan's standardized solution `std` of a fa_bifactor() fit with each
# defined index that uses a label several items share (equal_loading or
# equal_error; .fa_std_shared_terms()) evaluated at every item's own
# standardized values. With the factors' variances 1 and the factors
# uncorrelated, item i's variance is v_i = g_i^2 + s_i^2 + p_i (general
# loading g_i, group loading s_i, 0 for an item on the general factor
# alone, error variance p_i), and its standardized general and group
# loadings and error variance, the values standardizedSolution() reports
# in its own rows, are
#
#   a_i = g_i / v_i^{1/2},   b_i = s_i / v_i^{1/2},   e_i = p_i / v_i.
#
# The indices are the model's own definitions (.fa_layered_index_lines())
# written in a_i, b_i, and e_i, one symbol per item, so that omega,
# omega_H, and omega_HS are those of the sums of the standardized items
# and ECV and H, which do not depend on the metric of the items, equal
# the rows of the unstandardized tables. Their derivatives with respect
# to a_i, b_i, and e_i are R's symbolic ones (stats::deriv()), carried to
# the parameters by the chain rule,
#
#   da/dg = (s^2 + p) / v^{3/2},  da/ds = -g s / v^{3/2},
#   da/dp = -g / (2 v^{3/2}),     db/dg = -g s / v^{3/2},
#   db/ds = (g^2 + p) / v^{3/2},  db/dp = -s / (2 v^{3/2}),
#   de/dg = -2 g p / v^2,         de/ds = -2 s p / v^2,
#   de/dp = (g^2 + s^2) / v^2,
#
# and the standard errors and limits are .fa_std_redefine()'s. The rows
# that use no shared label are lavaan's, untouched. Returns `std`.
.fa_bifactor_std_defined <- function(fit, std, factors, general,
                                     general_only, conf_level) {
  pt <- lavaan::parTable(fit)
  terms <- .fa_std_shared_terms(pt)
  if (!length(terms)) return(std)
  key <- paste(pt$lhs, pt$op, pt$rhs)
  groups <- names(factors)
  items <- c(unlist(factors, use.names = FALSE), general_only)
  n <- length(items)
  group <- rep(NA_character_, n)
  for (f in groups) group[match(factors[[f]], items)] <- f
  alone <- is.na(group)
  row_g <- match(paste(general, "=~", items), key)
  row_s <- rep(NA_integer_, n)
  row_s[!alone] <- match(paste(group[!alone], "=~", items[!alone]), key)
  row_p <- match(paste(items, "~~", items), key)
  sym_a <- paste0("a_", seq_len(n))
  sym_b <- ifelse(alone, NA_character_, paste0("b_", seq_len(n)))
  sym_e <- paste0("e_", seq_len(n))
  lines <- .fa_layered_index_lines(general, groups, sym_a, sym_b, sym_e,
                                   group)
  def_lhs <- sub(" := .*$", "", lines)
  def_rhs <- sub("^[^ ]+ := ", "", lines)
  terms <- terms[terms %in% def_lhs]
  if (!length(terms)) return(std)
  exprs <- lapply(def_rhs[match(terms, def_lhs)], str2lang)
  vars <- c(sym_a, sym_b[!alone], sym_e)
  d_exprs <- lapply(exprs, stats::deriv, namevec = vars)
  fun <- function(est, grad = TRUE) {
    g <- est[row_g]
    s <- numeric(n)
    s[!alone] <- est[row_s[!alone]]
    p <- est[row_p]
    v <- g^2 + s^2 + p
    at <- stats::setNames(
      as.list(c(g / sqrt(v), (s / sqrt(v))[!alone], p / v)), vars)
    if (!grad) {
      return(list(value = stats::setNames(
        vapply(exprs, eval, numeric(1L), envir = at, enclos = baseenv()),
        terms)))
    }
    out <- lapply(d_exprs, eval, envir = at, enclos = baseenv())
    value <- stats::setNames(vapply(out, as.numeric, numeric(1L)), terms)
    D <- do.call(rbind, lapply(out, function(o) attr(o, "gradient")[1L, ]))
    D <- matrix(D, nrow = length(terms), dimnames = list(NULL, vars))
    d_a <- D[, sym_a, drop = FALSE]
    d_b <- matrix(0, nrow(D), n)
    d_b[, !alone] <- D[, sym_b[!alone]]
    d_e <- D[, sym_e, drop = FALSE]
    w <- v^1.5
    chain <- function(da, db, de) {
      d_a * rep(da, each = nrow(D)) + d_b * rep(db, each = nrow(D)) +
        d_e * rep(de, each = nrow(D))
    }
    G <- matrix(0, nrow(D), length(est))
    # The rows of the general loadings, the group loadings (items on a
    # group factor), and the error variances.
    G[, row_g] <- G[, row_g] +
      chain((s^2 + p) / w, -g * s / w, -2 * g * p / v^2)
    d_s <- chain(-g * s / w, (g^2 + p) / w, -2 * s * p / v^2)
    G[, row_s[!alone]] <- G[, row_s[!alone]] + d_s[, !alone, drop = FALSE]
    G[, row_p] <- G[, row_p] +
      chain(-g / (2 * w), -s / (2 * w), (g^2 + s^2) / v^2)
    list(value = value, grad = G)
  }
  .fa_std_redefine(fit, std, list(list(terms = terms, fun = fun)),
                   conf_level, pt)
}

# .fa_std_redefine(fit, std, blocks, conf_level, pt)
#
# lavaan's standardized solution `std` with the defined rows of each block
# of `blocks` computed again: a block is a list of `terms`, the names of
# its rows, and `fun(est, grad)`, which returns, at the parameter table
# values `est`, the list of their `value` and, with `grad`, the matrix
# `grad` of their derivatives with respect to the rows of the parameter
# table `pt`. The standard errors are those of the delta method over the
# full covariance matrix of the fit's free parameters,
#
#     Var = g' V g,   g = the derivatives with respect to the parameters,
#
# with g placed in the layout of lavaan's vcov(): one row for each row of
# the parameter table with free > 0, in table order, which every row that
# shares a label occupies under ceq.simple = TRUE, or one row per free
# number (see .H_from_fit()). This is the computation lavaan's
# standardizedSolution() makes for every other row. Under se =
# "bootstrap" they are, as lavaan's are, the standard deviation of the
# quantities over the bootstrap draws of the free parameters (divisor
# the number of draws), with lavaan's percentile limits
# (.fa_cfa_boot_limits()). The z, p-value, and Wald limits follow from
# the estimate and standard error as lavaan forms them (z is NA for a
# standard error of 0, and p is 2 (1 - Phi(|z|))); the caller then
# withholds them when the covariance matrix is unusable
# (.cfa_k_no_vcov()), refers them to t under wald_df, and forms the
# logit-scale limits of the indices. Without standard errors (se =
# "none") only the estimates are replaced. Returns `std`.
.fa_std_redefine <- function(fit, std, blocks, conf_level,
                             pt = lavaan::parTable(fit)) {
  if (!length(blocks)) return(std)
  free <- pt$free
  is_free <- free > 0L
  has_se <- "se" %in% names(std)
  boot <- has_se &&
    identical(lavaan::lavInspect(fit, "options")$se, "bootstrap")
  V <- NULL
  draws <- NULL
  if (has_se && !boot) {
    V <- tryCatch(unclass(lavaan::vcov(fit)), error = function(e) NULL)
  }
  if (boot) draws <- .fa_cfa_boot_est(fit, pt)
  crit <- stats::qnorm(1 - (1 - conf_level) / 2)
  def_lhs <- ifelse(std$op == ":=", std$lhs, NA_character_)
  for (b in blocks) {
    at <- b$fun(pt$est, grad = has_se && !boot)
    k <- length(b$terms)
    idx <- match(b$terms, def_lhs)
    std$est.std[idx] <- unname(at$value)
    if (!has_se) next
    lim <- NULL
    if (boot) {
      vals <- if (!is.null(draws)) {
        do.call(rbind, lapply(seq_len(nrow(draws)), function(r) {
          b$fun(draws[r, ], grad = FALSE)$value
        }))
      }
      se <- if (is.null(vals)) {
        rep(NA_real_, k)
      } else {
        apply(vals, 2L, function(x) {
          x <- x[is.finite(x)]
          stats::sd(x) * sqrt((length(x) - 1) / length(x))
        })
      }
      if (!is.null(vals)) {
        lim <- apply(vals, 2L, .fa_cfa_boot_limits, conf_level = conf_level)
      }
    } else {
      # The derivatives with respect to the parameter table rows, in the
      # layout of V.
      g <- if (is.null(V)) {
        NULL
      } else if (nrow(V) == sum(is_free)) {
        at$grad[, is_free, drop = FALSE]
      } else if (nrow(V) == max(c(0L, free))) {
        t(rowsum(t(at$grad[, is_free, drop = FALSE]), free[is_free],
                 reorder = TRUE))
      }
      se <- if (is.null(g)) {
        rep(NA_real_, k)
      } else {
        sqrt(pmax(rowSums((g %*% V) * g), 0))
      }
    }
    se <- unname(se)
    std$se[idx] <- se
    z <- ifelse(is.finite(se) & se > 0, unname(at$value) / se, NA_real_)
    if ("z" %in% names(std)) std$z[idx] <- z
    if ("pvalue" %in% names(std)) {
      std$pvalue[idx] <- 2 * (1 - stats::pnorm(abs(z)))
    }
    if (all(c("ci.lower", "ci.upper") %in% names(std))) {
      if (is.null(lim)) {
        lim <- rbind(at$value - crit * se, at$value + crit * se)
      }
      std$ci.lower[idx] <- unname(lim[1L, ])
      std$ci.upper[idx] <- unname(lim[2L, ])
    }
  }
  std
}

# .fa_cfa_boot_est(fit, pt)
#
# The bootstrap draws of a fit with se = "bootstrap" as parameter table
# `est` columns, one row per draw (lavaan's draws that failed left out),
# or NULL when lavaan kept none. The draws are laid out by free number, or
# by the rows with free > 0 under ceq.simple = TRUE (see
# .fa_std_redefine()).
.fa_cfa_boot_est <- function(fit, pt) {
  bc <- tryCatch(lavaan::lavInspect(fit, "coef.boot"),
                 error = function(e) NULL)
  if (is.null(bc) || !length(bc)) return(NULL)
  err <- attr(bc, "error.idx")
  bc <- matrix(unclass(bc), nrow = NROW(bc))
  if (length(err)) bc <- bc[-err, , drop = FALSE]
  bc <- bc[rowSums(!is.finite(bc)) == 0L, , drop = FALSE]
  if (!nrow(bc)) return(NULL)
  free <- pt$free
  is_free <- free > 0L
  pos <- if (ncol(bc) == sum(is_free)) {
    seq_len(sum(is_free))
  } else if (ncol(bc) == max(c(0L, free))) {
    free[is_free]
  } else {
    return(NULL)
  }
  t(apply(bc, 1L, function(x) {
    est <- pt$est
    est[is_free] <- x[pos]
    est
  }))
}

# .fa_cfa_boot_limits(x, conf_level)
#
# The percentile limits lavaan's standardizedSolution() gives a row under
# se = "bootstrap" (its default boot_ci_type = "perc"): the order
# statistics at (R + 1) alpha of the R finite draws, interpolated on the
# normal quantile scale between neighbors when (R + 1) alpha is not
# whole, and the extreme draw beyond them (the rule of the boot package's
# norm.inter(), which lavaan follows). Returns the lower and upper
# limits.
.fa_cfa_boot_limits <- function(x, conf_level) {
  x <- sort(x[is.finite(x)])
  R <- length(x)
  if (R < 2L) return(c(NA_real_, NA_real_))
  alpha <- (1 + c(-conf_level, conf_level)) / 2
  vapply(alpha, function(a) {
    rk <- (R + 1) * a
    k <- trunc(rk)
    if (k <= 0) return(x[1L])
    if (k >= R) return(x[R])
    if (k == rk) return(x[k])
    z <- stats::qnorm(c(a, k / (R + 1), (k + 1) / (R + 1)))
    x[k] + (z[1L] - z[2L]) / (z[3L] - z[2L]) * (x[k + 1L] - x[k])
  }, numeric(1))
}

# .fa_cfa_fit_rows(fit, h1_diagnostic, rmsea_conf_level, has_means,
#                  ts_missing, two_stage, ordered_items)
#
# The fit rows of a fa_confirmatory() table (chi_square through H1) and
# what the attributes and the lines under the table need to describe them.
# `rmsea_conf_level` is the level of the RMSEA interval, the only interval
# among the fit rows, which the rmsea_ci_level row records; it is
# rmsea_conf_level of the user's call, not the conf_level of the parameter
# intervals.
# Returns a list: `rmsea_conf_level`, the level of the RMSEA interval;
# `rows`, the rows in the table's columns; `test_label`,
# the statistic in the chi_square row; `fiml`, whether the fit used FIML;
# `fit_options`, lavaan's options of the fit; `index_method`,
# `index_candidates`, and `index_version`, the version of the index rows
# and how it was chosen; `saturated` and `zero_df_misfit`, for a model with
# no degrees of freedom; `mva`, the robust indices DMAR computed (empty
# when none); and `scaled`, whether the fit has a scaled statistic.
.fa_cfa_fit_rows <- function(fit, h1_diagnostic, rmsea_conf_level, has_means,
                             ts_missing, two_stage, ordered_items) {
  fm <- if (is.null(h1_diagnostic)) {
    lavaan::fitMeasures(fit, fm.args = list(rmsea.ci.level = rmsea_conf_level))
  } else {
    ll <- lavaan::logLik(fit)
    npar <- attr(ll, "df")
    nobs <- sum(lavaan::lavInspect(fit, "nobs"))
    c(df = lavaan::lavInspect(fit, "test")$standard$df,
      aic = -2 * as.numeric(ll) + 2 * npar,
      bic = -2 * as.numeric(ll) + log(nobs) * npar, logl = as.numeric(ll))
  }
  # The robust indices of a continuous ML fit with the older mean and
  # variance adjusted statistic (test = "mean.var.adjusted" through ...),
  # which lavaan does not compute, at the same requested confidence
  # level as the other RMSEA intervals.
  mva <- if (is.null(h1_diagnostic)) .sem_robust_mva(fit, conf_level = rmsea_conf_level) else numeric()
  if (length(mva)) fm <- c(fm, mva)
  # The first of the named measures lavaan computed for this fit, and its
  # value.
  pick_name <- function(...) {
    for (nm in c(...)) {
      if (nm %in% names(fm) && is.finite(fm[nm])) return(nm)
    }
    NA_character_
  }
  pick <- function(...) {
    nm <- pick_name(...)
    if (is.na(nm)) NA_real_ else unname(fm[nm])
  }
  fit_options <- lavaan::lavInspect(fit, "options")
  # Whether the fit used FIML, read from the fit rather than from the
  # argument: lavaan fits the complete cases instead, with a warning, when
  # missing = "ml" meets a test it cannot combine with FIML (test =
  # "satorra.bentler" or "mean.var.adjusted" passed through ...), and the
  # indices of such a fit are not the FIML corrected ones.
  fiml <- isTRUE(fit_options$missing %in% c("ml", "ml.x"))
  # One version holds every index row (cfi, tli, nnfi, rmsea, and its
  # interval), chosen as measurement_invariance() chooses the version of
  # its columns (.measurement_invariance_candidates() and
  # .measurement_invariance_choose()), and the lines under the table name
  # it. The robust indices are reported whenever lavaan (or, for the mean
  # and variance adjusted statistic, DMAR) computed them. Otherwise a
  # maximum likelihood fit of continuous items, complete or listwise,
  # reports the standard ones, which estimate the population values of the
  # ML indices with the small sample bias corrected as if the data were
  # normal, and not the scaled ones, which put the scaled statistic into
  # the standard formulas and estimate population values that change with
  # the kurtosis of the data (Savalei, 2018); a FIML or two-stage fit, or a
  # fit of ordered items or by least squares, reports the scaled versions
  # and then the standard ones. Through 1.0.1 each row took the first of
  # its robust, scaled, and standard versions that lavaan had computed, so
  # a fit with the mean and variance adjusted statistic reported the
  # scaled indices, and a row whose robust version was NA could hold
  # another version than the rows beside it.
  index_method <- if (fiml) {
    "fiml"
  } else if (two_stage) {
    "two_stage"
  } else if (is.null(ordered_items) &&
             identical(fit_options$estimator, "ML")) {
    "ml"
  } else "other"
  index_candidates <- .measurement_invariance_candidates(
    "indices", length(.lavaan_corrected_tests(fit_options$test)) > 0L,
    index_method)
  index_version <- .measurement_invariance_choose(list(fm), index_candidates)
  # A model with no degrees of freedom reproduces the moments exactly, so
  # its statistic is zero in every version and the versions do not differ:
  # its CFI is 1 in every version, the model term T - c df being zero
  # whatever the scaling constant c (0/0 there, which is how lavaan's
  # robust indices come out NA or not by rounding), and its TLI and RMSEA,
  # which divide by df, are 0/0 in every version, for which lavaan's
  # standard version reports 1 and 0 (and the limits 0) by convention.
  # Those rows take the standard values, as 1.0.1's row by row choice
  # gave them when the robust ones were NA, and no line under the table
  # names a version or a fallback from one, since the rows of such a model
  # do not depend on it. A model with no degrees of freedom whose statistic
  # is not zero, as when an inequality bound binds, does not fit exactly
  # (.sem_zero_df_misfit()): its cfi is formed from its statistic in the
  # version chosen, which its CFI alone decides
  # (.measurement_invariance_choose()), and its tli, nnfi, rmsea, and
  # RMSEA limits, which divide by the degrees of freedom, are undefined
  # and NA, as sem_fit_measures() reports them, with a line under the table
  # that says so; the lines that name the version name the cfi alone. On one
  # factor of three variables correlated .8, .8, and .5 with bounds =
  # "pos.var" (N = 200) the cfi is .925, as in 1.0.1 and lavaan, where the
  # first version of the rule above set it to 1; through 1.0.1 (and that
  # version) the other rows held lavaan's TLI and NNFI of 1 and RMSEA of 0
  # with limits of 0, the values of exact fit, beside a CFI that showed
  # misfit.
  saturated <- .sem_exact_fit(fm)
  zero_df_misfit <- .sem_zero_df_misfit(fm)
  pick_version <- function(nm) {
    if (saturated && nm == "cfi") return(1)
    if (zero_df_misfit && nm %in% c("tli", "nnfi", "rmsea", "rmsea.ci.lower",
                                    "rmsea.ci.upper")) {
      return(NA_real_)
    }
    key <- paste0(nm, if (saturated) "" else c(robust = ".robust",
                                               scaled = ".scaled",
                                               standard = "")[[index_version]])
    if (key %in% names(fm)) unname(fm[[key]]) else NA_real_
  }
  # The test of exact fit. A robust or categorical estimator reports its
  # scaled statistic. After normal theory two-stage estimation lavaan's
  # standard statistic is the naive one, which treats the EM moments as
  # complete data, and its scaled statistic is an approximation; the row
  # is the residual based statistic of Savalei and Bentler (2009), which is
  # asymptotically chi square. After robust two-stage estimation it is
  # lavaan's scaled statistic (Savalei & Falk, 2014). The indices below
  # (lavaan's robust ones: after two-stage estimation the corrected indices
  # of Zhang & Savalei, 2023, TS V2 after robust two-stage estimation and
  # its normal theory analog after normal theory two-stage estimation; see
  # .ts_indices_label()) are lavaan's whatever the test, since they are
  # computed from the naive statistic and its expectation, not from the
  # test.
  if (identical(ts_missing, "two.stage")) {
    rt <- .ts_residual_test(fit)
    chi <- c(rt$statistic, rt$df, rt$p_value)
  } else {
    chi <- c(pick("chisq.scaled", "chisq"), pick("df.scaled", "df"),
             pick("pvalue.scaled", "pvalue"))
  }
  test_label <- if (two_stage) {
    .ts_test_label(ts_missing)
  } else if ("chisq.scaled" %in% names(fm) &&
             is.finite(fm[["chisq.scaled"]])) {
    "scaled"
  } else {
    "standard"
  }
  fit_rows <- data.frame(
    syntax = "",
    term = c("chi_square", "df", "p_chi_square", "cfi", "tli", "nnfi",
             "rmsea", "rmsea_ci_lower", "rmsea_ci_upper", "rmsea_ci_level",
             "srmr", "AIC", "BIC", "H0", "H1"),
    estimate = c(
      chi,
      pick_version("cfi"), pick_version("tli"), pick_version("nnfi"),
      pick_version("rmsea"), pick_version("rmsea.ci.lower"),
      pick_version("rmsea.ci.upper"),
      { lv <- pick("rmsea.ci.level"); if (is.na(lv)) rmsea_conf_level else lv },
      # With auxiliary variables the SRMR is over the model's own items,
      # as in sem_fit_measures() (.sem_srmr_own()).
      { aux <- if (is.null(h1_diagnostic)) .sem_aux_names(fit) else character(0)
        sm <- if (length(aux)) {
          .sem_srmr_own(fit, aux)
        } else if (has_means) pick("srmr_bentler_nomean") else NA_real_
        if (is.na(sm)) pick("srmr") else sm },
      pick("aic"), pick("bic"),
      pick("logl"), pick("unrestricted.logl")),
    se = NA_real_, test_statistic = NA_real_, p_value = NA_real_,
    ci_lower = NA_real_, ci_upper = NA_real_,
    stringsAsFactors = FALSE
  )

  if (!is.null(h1_diagnostic)) {
    h1_terms <- c("chi_square", "p_chi_square", "cfi", "tli", "nnfi",
                  "rmsea", "rmsea_ci_lower", "rmsea_ci_upper", "srmr", "H1")
    fit_rows$estimate[fit_rows$term %in% h1_terms] <- NA_real_
  }
  list(rmsea_conf_level = rmsea_conf_level,
       rows = fit_rows, test_label = test_label, fiml = fiml,
       fit_options = fit_options, index_method = index_method,
       index_candidates = index_candidates, index_version = index_version,
       saturated = saturated, zero_df_misfit = zero_df_misfit, mva = mva,
       scaled = is.finite(pick("chisq.scaled")))
}

# .fa_cfa_fit_attributes(out, fr, estimator, missing, input_is_raw_data,
#                        two_stage, ts_missing, ts_rows, ts_aux, complete,
#                        h1_diagnostic)
#
# Records on the table `out` what produced its fit rows (`fr`, from
# .fa_cfa_fit_rows()): the "test", "estimator", "missing",
# "fit_indices", and "rmsea_conf_level" attributes that glance(),
# dmar_metadata(), and the line stating the table's confidence level read, the
# lines under the table that name the version of the rows ("fit_rows"),
# and after two-stage estimation the "two_stage" record. `complete` says
# that no item value is missing in the fit (a covariance matrix, or raw
# data without a missing item value). Returns `out`.
.fa_cfa_fit_attributes <- function(out, fr, estimator, missing,
                                   input_is_raw_data, two_stage, ts_missing,
                                   ts_rows, ts_aux, complete, h1_diagnostic) {
  test_label <- fr$test_label
  fiml <- fr$fiml
  fit_options <- fr$fit_options
  index_version <- fr$index_version
  saturated <- fr$saturated
  zero_df_misfit <- fr$zero_df_misfit
  attr(out, "test") <- test_label
  # What produced the fit rows, for glance(), whose row would otherwise
  # look alike for a listwise, FIML, and two-stage fit: the estimator, the
  # missing data method lavaan used (it falls back to listwise deletion
  # for robust.sem under FIML), and the version of the index rows (the
  # standard keys for a saturated model, whose rows are those of exact
  # fit in every version).
  attr(out, "estimator") <- estimator
  attr(out, "missing") <- if (!input_is_raw_data) {
    "none"
  } else if (identical(ts_missing, "robust.two.stage")) {
    "robust_two_stage"
  } else if (two_stage) {
    "two_stage"
  } else if (fiml) {
    "fiml"
  } else {
    fit_options$missing %||% NA_character_
  }
  attr(out, "fit_indices") <- if (saturated) "standard" else index_version
  # The level of the RMSEA interval, apart from the table's conf_level
  # (the level of the parameter intervals), so that the line stating the
  # table's level names it too (.rmsea_level_phrase(),
  # R/dmar_output_helpers.R).
  attr(out, "rmsea_conf_level") <- fr$rmsea_conf_level
  # The rows say which version they hold where the row name alone does
  # not: after FIML (and two-stage estimation, whose lines .ts_footer()
  # writes) the corrected indices of Zhang and Savalei (2023) are not
  # computed from the chi_square row, a fallback from the robust indices
  # is named, robust indices DMAR computed are named, and lavaan's SRMR of
  # a model with a mean structure is not the SRMR sem_fit_measures()
  # reports.
  # A saturated model's index rows hold the values of exact fit, which no
  # version changes, so neither the FIML line nor the two-stage line that
  # names the corrected indices of Zhang and Savalei (2023) is written for
  # it (its rows are not theirs: after two-stage estimation lavaan leaves
  # its robust CFI and TLI NA there); the lines that describe the
  # estimation, the test, the criteria, the number of cases, and the SRMR
  # stay. Through the first version of the rule for saturated models these
  # lines were written, though the help said no line names a version. A
  # model with no degrees of freedom whose statistic is not zero has its
  # cfi alone in a version, and the lines name that row alone: through the
  # first version of the rule for such a model the FIML and two-stage lines
  # named a correction for the tli, nnfi, and RMSEA interval rows, which
  # held lavaan's values of exact fit.
  robust_indices <- index_version == "robust"
  ts_footer <- if (two_stage) .ts_footer(ts_missing, "cfa", rows = ts_rows)
  if (two_stage && saturated) {
    ts_footer <- ts_footer[names(ts_footer) != "two_stage_indices"]
  } else if (two_stage && !robust_indices) {
    ts_footer[["two_stage_indices"]] <- .cfa_k_fallback_line(
      "two-stage", index_version, cfi_only = zero_df_misfit)
  } else if (two_stage && zero_df_misfit) {
    ts_footer[["two_stage_indices"]] <- paste0(
      "cfi: ", .ts_indices_label(ts_missing))
  }
  attr(out, "fit_rows") <- .cfa_k_fit_rows_lines(
    fiml = fiml, robust_indices = robust_indices, fallback = index_version,
    scaled = fr$scaled, complete = complete, method = fr$index_method,
    fell_back = !saturated &&
      index_version != names(fr$index_candidates)[1L],
    dmar_robust = !saturated && robust_indices && length(fr$mva) > 0L,
    saturated = saturated, zero_df_misfit = zero_df_misfit)
  if (!is.null(h1_diagnostic)) {
    attr(out, "fit_rows") <- paste0(
      "H1-based fit statistics, indices, and residual correlations are ",
      "unavailable; the constrained-model parameter rows and H0 criteria ",
      "remain, with the H1 qualification below.")
  }
  if (two_stage) {
    attr(out, "two_stage") <- list(missing = missing,
                                   lavaan_missing = ts_missing,
                                   test = test_label,
                                   rows = ts_rows[c("read", "nobs", "extra",
                                                    "corrected", "reestimated")],
                                   aux = if (length(ts_aux)) ts_aux,
                                   footer = ts_footer)
  }
  out
}

# ---------------------------------------------------------------------------
# The bifactor structure (fa_bifactor())
# ---------------------------------------------------------------------------
#
# A general factor on every item and two or more group factors, each on its
# own items, with every factor's variance fixed at 1 and every factor
# covariance at 0. `factors` names the group factors; `general` names the
# general factor, which also loads on the items of `general_only`, the
# items on the general factor alone. The items are ordered as the group
# factors list them, then `general_only`, and the general loadings, error
# variances, and intercepts of the items on the general factor alone are
# numbered in that order.

# .fa_bifactor_refuse(message, reason, factors, items)
#
# The classed refusal of a bifactor structure fa_bifactor() will not fit:
# an error of class c("dmar_bifactor_unidentified", "dmar_refusal",
# "error", "condition") with code "bifactor_unidentified", the factors and
# items at issue, and `reason`, which says what kind of refusal it is:
# "unidentified" (the data cannot identify the structure), "names" (a name
# placed twice, or names that clash), or "scope" (a structure outside the
# model fa_bifactor() fits, which may well be identified: fewer than two
# group factors, an item on two group factors), all in the condition
# (.dmar_refuse(), R/dmar_conditions.R).
.fa_bifactor_refuse <- function(message, reason, factors = character(0),
                                items = character(0)) {
  .dmar_refuse(message, "bifactor_unidentified", reason = reason,
               factors = factors, items = items)
}

# .fa_bifactor_items(factors, general, general_only)
#
# Checks the names of a bifactor structure, refusing (.fa_bifactor_refuse())
# a general factor name that is not one syntactically valid name or that a
# group factor also has, an item listed twice under one group factor, a
# factor named like an item, an item on two group factors, an item named
# twice in `general_only` or also on a group factor, and fewer than two
# group factors. Returns the items, as the group factors list them and then
# `general_only`.
.fa_bifactor_items <- function(factors, general, general_only) {
  if (!is.character(general) || length(general) != 1L || is.na(general) ||
      !nzchar(general) || general != make.names(general)) {
    .dmar_refuse_argument(c(
      "'general' must be the name of the general factor, one ",
      "syntactically valid R name (it becomes a lavaan label)."),
      arg = "general")
  }
  groups <- names(factors)
  for (f in groups) {
    if (!is.character(factors[[f]]) || !length(factors[[f]]) ||
        anyNA(factors[[f]])) {
      .dmar_refuse_argument(
        sprintf("Group factor '%s' must name its items.", f), arg = NULL)
    }
  }
  if (!is.null(general_only) &&
      (!is.character(general_only) || anyNA(general_only))) {
    .dmar_refuse_argument(c(
      "'general_only' must be NULL or a character vector naming the ",
      "items on the general factor alone."), arg = "general_only")
  }
  if (general %in% groups) {
    .fa_bifactor_refuse(c(
      "The general factor and a group factor are both named '", general,
      "'. Rename one of them ('general' names the general factor)."),
      "names", factors = general)
  }
  # An item listed twice under one group factor is a slip in the list, not
  # an item on two group factors, and is named as such.
  twice <- lapply(factors, function(v) unique(v[duplicated(v)]))
  twice <- twice[lengths(twice) > 0L]
  if (length(twice)) {
    .fa_bifactor_refuse(c(
      "Each item is listed once under its group factor, but ",
      paste0(vapply(twice, paste, character(1L), collapse = ", "),
             ifelse(lengths(twice) == 1L, " is", " are"),
             " listed more than once under ", names(twice),
             collapse = "; "),
      ". List each item once."),
      "names", factors = names(twice),
      items = unlist(twice, use.names = FALSE))
  }
  on_groups <- unlist(factors, use.names = FALSE)
  several <- unique(on_groups[duplicated(on_groups)])
  if (length(several)) {
    where <- vapply(several, function(i) {
      paste(groups[vapply(factors, function(v) i %in% v, logical(1L))],
            collapse = " and ")
    }, character(1L))
    .fa_bifactor_refuse(c(
      "Each item loads on one group factor in a bifactor model, but ",
      paste0(several, " is listed under ", where, collapse = "; "),
      ". An item on two group factors belongs to neither subscale alone, ",
      "and omega_HS, ECV, and PUC are defined by the subscales. List each ",
      "item under one group factor, or under 'general_only'."),
      "scope", factors = groups, items = several)
  }
  if (anyDuplicated(general_only)) {
    .fa_bifactor_refuse(c(
      "'general_only' names ",
      paste(unique(general_only[duplicated(general_only)]), collapse = ", "),
      " more than once."), "names", items = general_only)
  }
  both <- intersect(general_only, on_groups)
  if (length(both)) {
    .fa_bifactor_refuse(c(
      paste(both, collapse = ", "),
      if (length(both) == 1L) " is" else " are",
      " named both under a group factor and in 'general_only', the items ",
      "on the general factor alone. Keep each item in one place."),
      "names", items = both)
  }
  items <- c(on_groups, general_only)
  clash <- intersect(c(general, groups), items)
  if (length(clash)) {
    .fa_bifactor_refuse(c(
      "A factor and an item are both named ",
      paste0("'", clash, "'", collapse = ", "),
      ". Rename the factor, since lavaan reads the name as either one."),
      "names", factors = clash, items = clash)
  }
  # Fewer than two group factors is outside the model fa_bifactor() fits,
  # not a structure the data cannot identify: a general factor beside one
  # group factor and items on the general factor alone is identified (the
  # rank of the covariance Jacobian is full for three group items with two
  # or three items on the general factor alone), but omega_HS, ECV, and
  # PUC describe how the common variance divides among subscales.
  if (length(groups) < 2L) {
    .fa_bifactor_refuse(c(
      "fa_bifactor() fits a general factor beside two or more group ",
      "factors, and 'factors' names ",
      if (length(groups)) paste0("one ('", groups, "')") else "none",
      ". This is the scope of the function rather than a question of ",
      "identification: its omega_HS, ECV, and PUC describe how the items' ",
      "common variance divides among subscales, which takes two or more. ",
      "Name every group factor in 'factors'; a general factor beside a ",
      "single group factor can be fit with lavaan::cfa(..., ",
      "orthogonal = TRUE, std.lv = TRUE)."), "scope", factors = groups)
  }
  items
}

# .fa_bifactor_identification(factors, general, general_only, eq_load,
#                             eq_err)
#
# Refuses (.fa_bifactor_refuse()) a group factor that its items cannot
# identify, and the two group factor structures whose general loadings are
# fixed only up to a scale, and warns, with class dmar_bifactor_empirical,
# for every other structure with two group factors.
#
# A group factor needs three items, or two with equal loadings on it: a
# group factor of one item cannot be told apart from the item's error
# variance, and two free loadings on a factor of variance 1 enter the
# model only through their product. With two group factors and no item on
# the general factor alone, the covariances between the two groups' items,
# l_i l_j for i in one group and j in the other, fix the general loadings
# of one group only up to a factor c and those of the other up to 1/c, and
# c must come from inside a group: the group's item variances and
# covariances, J(J + 1)/2 for J items, must outnumber the parameters the
# group estimates beside the general loadings, its J group loadings (1 when
# equated) and J error variances (1 when equated). For a group of three
# items with free loadings and error variances they are 6 against 6, and
# for two items with equal loadings 3 against 3, so two such groups leave c
# undetermined: the model is not identified, and lavaan's information
# matrix is singular (on holzinger_swineford, two groups of three tests
# with nothing on the general factor alone could not be inverted). A group
# of four items (10 against 8), equal loadings or error variances in a
# group of three, equal general loadings (which tie the two groups' scales
# together), or an item on the general factor alone (whose covariances
# with the items of the two groups give c^2 directly) each fix c. The
# rule is a count, so it gives the generic answer; Green and Yang (2018)
# show that a structure it accepts can still be empirically
# underidentified with two group factors, when the item covariances are
# nearly equal within and between the groups, which is what the warning
# for two group factors says.
.fa_bifactor_identification <- function(factors, general, general_only,
                                        eq_load, eq_err) {
  groups <- names(factors)
  J <- lengths(factors)
  thin <- groups[J < 3L & !(J == 2L & eq_load[groups])]
  if (length(thin)) {
    one <- thin[J[thin] == 1L]
    two <- setdiff(thin, one)
    .fa_bifactor_refuse(c(
      "A group factor needs three or more items, or two with equal ",
      "loadings on it: ",
      paste(c(
        if (length(two)) paste0(
          paste0("'", two, "'", collapse = ", "), " ",
          if (length(two) == 1L) "has" else "have",
          " two items with free loadings, which enter the model only ",
          "through their product (equal_loading = c(", two[1L], " = TRUE, ",
          "...) equates them)"),
        if (length(one)) paste0(
          paste0("'", one, "'", collapse = ", "), " ",
          if (length(one) == 1L) "has" else "have",
          " one item, whose group loading cannot be told apart from its ",
          "error variance (put such an item in 'general_only')")),
        collapse = "; "), "."),
      "unidentified", factors = thin,
      items = unlist(factors[thin], use.names = FALSE))
  }
  if (length(groups) != 2L) return(invisible(TRUE))
  # The count of the comment above: does a group fix the scale its general
  # loadings share with the other group's?
  spare <- vapply(groups, function(f) {
    j <- J[[f]]
    j * (j + 1) / 2 - (if (eq_load[[f]]) 1 else j) -
      (if (eq_err[[f]]) 1 else j)
  }, numeric(1L))
  if (!length(general_only) && !eq_load[[general]] && all(spare <= 0)) {
    .fa_bifactor_refuse(c(
      "With two group factors and no item on the general factor alone, ",
      "the covariances between the items of '", groups[1L], "' and '",
      groups[2L], "' fix the general loadings of one group only up to a ",
      "factor c and those of the other up to 1/c, and c must come from ",
      "inside a group, whose item variances and covariances must outnumber ",
      "the group loadings and error variances it estimates. Here neither ",
      "group's do (",
      paste0("'", groups, "' with ", J, " items", collapse = " and "),
      "), so the model is not identified. Any of these identifies it: an ",
      "item on the general factor alone ('general_only'), a group of four ",
      "or more items, equal loadings or equal error variances within a ",
      "group of three, or equal general loadings."),
      "unidentified", factors = groups,
      items = unlist(factors, use.names = FALSE))
  }
  .dmar_warn(c(
    "A bifactor model with two group factors ('", groups[1L], "' and '",
    groups[2L], "') is identified by the count of its parameters, but it ",
    "can be empirically underidentified: when the item covariances are ",
    "nearly equal within each group and between the groups, different ",
    "general and group loadings reproduce them equally well (Green & ",
    "Yang, 2018). Items on the general factor alone, or items whose ",
    "loadings differ in strength, protect against it; an information ",
    "matrix lavaan cannot invert, or estimates that change with the ",
    "starting values, are its symptoms."), "bifactor_empirical",
    factors = groups)
  invisible(TRUE)
}

# .fa_bifactor_labels(factors, general, general_only, eq_load, eq_int,
#                     eq_err)
#
# The lavaan labels of the bifactor model, one row per item in the order of
# the items: the item, its group factor (NA for an item on the general
# factor alone), and the labels of its general loading, group loading (NA
# when none), error variance, and intercept. A loading, error variance, or
# intercept equated within a factor shares one label (lambda_<f>, psi_<f>,
# nu_<f>); otherwise each is numbered within its factor, and an item on the
# general factor alone is numbered as the general factor numbers its items.
.fa_bifactor_labels <- function(factors, general, general_only, eq_load,
                                eq_int, eq_err) {
  items <- c(unlist(factors, use.names = FALSE), general_only)
  lab <- function(prefix, f, j, equal) {
    if (equal) rep(paste0(prefix, "_", f), length(j)) else
      paste0(prefix, "_", f, "_", j)
  }
  out <- data.frame(item = items, group = NA_character_,
                    general = lab("lambda", general, seq_along(items),
                                  eq_load[[general]]),
                    loading = NA_character_, error = NA_character_,
                    intercept = NA_character_, stringsAsFactors = FALSE)
  for (f in names(factors)) {
    rows <- match(factors[[f]], items)
    j <- seq_along(rows)
    out$group[rows] <- f
    out$loading[rows] <- lab("lambda", f, j, eq_load[[f]])
    out$error[rows] <- lab("psi", f, j, eq_err[[f]])
    out$intercept[rows] <- lab("nu", f, j, eq_int[[f]])
  }
  alone <- is.na(out$group)
  out$error[alone] <- paste0("psi_", general, "_", which(alone))
  out$intercept[alone] <- paste0("nu_", general, "_", which(alone))
  out
}

# .fa_bifactor_index_terms(general, groups)
#
# The bifactor indices lavaan defines, in the order of the table: for the
# general factor omega (total), omega_H, ECV, and H; for each group factor
# its subscale omega, omega_HS, ECV, and H.
.fa_bifactor_index_terms <- function(general, groups) {
  c(paste0(c("omega_", "omega_H_", "ECV_", "H_"), general),
    as.vector(t(outer(groups, c("omega_", "omega_HS_", "ECV_", "H_"),
                      function(f, p) paste0(p, f)))))
}

# .fa_bifactor_syntax(factors, general, general_only, eq_load, eq_int,
#                     eq_err, has_means)
#
# The lavaan syntax of the bifactor model: the general factor on every item
# and each group factor on its own, every factor's variance fixed at 1
# (and its mean at 0 with a mean structure) with its first loading freed,
# every factor covariance fixed at 0, every loading, error variance, and
# intercept labeled (.fa_bifactor_labels()), and the indices defined from
# the labels, so that lavaan's delta method gives each a standard error.
# With unit factor variances and uncorrelated factors, the model implied
# variance of the unit-weighted sum of the items is
# (sum of the general loadings)^2 + the sum over the group factors of
# (sum of the group loadings)^2 + the sum of the error variances, and
#   omega_<g>    = [(sum l_g)^2 + sum_s (sum l_s)^2] / that variance
#                  (omega total; Rodriguez, Reise, and Haviland, 2016,
#                  Eq. 4, there for standardized items),
#   omega_H_<g>  = (sum l_g)^2 / that variance (their Eq. 6),
#   omega_<s>    = [(sum_{i in s} l_g)^2 + (sum_{i in s} l_s)^2] / the
#                  variance of the sum of group s's items (Eq. 5),
#   omega_HS_<s> = (sum_{i in s} l_s)^2 / that variance (Eq. 7),
# all four for the sum of the raw items, as fa_confirmatory()'s omega and
# semTools::compRelSEM() compute them. ECV (their Eq. 10) is formed from
# the standardized loadings, each loading over the square root of its
# item's model implied variance v_i = l_gi^2 + l_si^2 + psi_i, as
# bifactor_indices() forms it: for the general factor
# sum_i (l_gi^2 / v_i) / sum_i ((l_gi^2 + l_si^2) / v_i), and for group
# factor s the same sums over its items (the general factor's share of
# the subscale's common variance). H (their Eq. 9) is Hancock and
# Mueller's closed form with every other source of an item's variance
# counted as its error, 1 / (1 + 1 / sum_i l_fi^2 / (v_i - l_fi^2)), which
# does not depend on the metric of the items. Returns one string.
.fa_bifactor_syntax <- function(factors, general, general_only, eq_load,
                                eq_int, eq_err, has_means) {
  lab <- .fa_bifactor_labels(factors, general, general_only, eq_load, eq_int,
                             eq_err)
  groups <- names(factors)
  items <- lab$item
  # Labels of different parameters that coincide (a group factor named
  # "g_1" with equal loadings beside a general factor "g" would label both
  # lambda_g_1), or that coincide with an index or a variable, would be
  # read by lavaan as one: refused.
  alone <- is.na(lab$group)
  n <- nrow(lab)
  label <- c(lab$general, lab$loading, lab$error,
             if (has_means) lab$intercept)
  owner <- c(rep(paste("lambda", general), n), paste("lambda", lab$group),
             ifelse(alone, paste("psi alone", seq_len(n)),
                    paste("psi", lab$group)),
             if (has_means) {
               ifelse(alone, paste("nu alone", seq_len(n)),
                      paste("nu", lab$group))
             })
  keep <- !is.na(label)
  label <- label[keep]
  owner <- owner[keep]
  shared <- tapply(owner, label, function(o) length(unique(o)))
  index_terms <- .fa_bifactor_index_terms(general, groups)
  clash <- unique(c(names(shared)[shared > 1L],
                    intersect(label, c(index_terms, items, general, groups)),
                    index_terms[duplicated(index_terms)],
                    intersect(index_terms, c(items, general, groups))))
  if (length(clash)) {
    .fa_bifactor_refuse(c(
      "The names of the factors make lavaan labels that coincide (",
      paste(clash, collapse = ", "), "), which lavaan would read as one ",
      "parameter. Rename the factors."), "names",
      factors = c(general, groups))
  }
  first_free <- function(f, v, l) {
    paste0(f, " =~ NA*", v[1L], " + ", paste0(l, "*", v, collapse = " + "))
  }
  # (paste0() would turn an empty set of items into one empty line.)
  each <- function(v, ...) if (length(v)) paste0(v, ...) else character(0)
  lines <- c(first_free(general, items, lab$general),
             paste0(general, " ~~ 1*", general),
             each(items[alone], " ~~ ", lab$error[alone], "*", items[alone]))
  if (has_means) {
    lines <- c(lines,
               each(items[alone], " ~ ", lab$intercept[alone], "*1"),
               paste0(general, " ~ 0*1"))
  }
  for (f in groups) {
    r <- which(lab$group %in% f)
    v <- items[r]
    lines <- c(lines, first_free(f, v, lab$loading[r]),
               paste0(f, " ~~ 1*", f),
               paste0(v, " ~~ ", lab$error[r], "*", v))
    if (has_means) {
      lines <- c(lines, paste0(v, " ~ ", lab$intercept[r], "*1"),
                 paste0(f, " ~ 0*1"))
    }
  }
  latent <- c(general, groups)
  pairs <- utils::combn(latent, 2L)
  lines <- c(lines, paste0(pairs[1L, ], " ~~ 0*", pairs[2L, ]))
  # The indices, written out in the labels.
  paste(c(lines, .fa_layered_index_lines(general, groups, lab$general,
                                         lab$loading, lab$error, lab$group)),
        collapse = "\n")
}

# .fa_layered_index_lines(general, groups, lg, ls, ps, group)
#
# The lavaan lines that define the indices of a general factor and its
# group factors from the parameters of the items they describe, one
# element per item in each of `lg` (the item's general loading), `ls` (its
# group loading, NA for an item on the general factor alone), `ps` (its
# error variance), and `group` (its group factor, NA when it has none),
# each a lavaan label or a parenthesized expression in labels; `groups`
# names the group factors in the order of their rows. The general factor
# and the group factors are uncorrelated with variance 1, so the variance
# of the unit-weighted sum of the items is (sum of the general loadings)^2
# + the sum over the group factors of (sum of the group loadings)^2 + the
# sum of the error variances, and the indices are those of the comment of
# .fa_bifactor_syntax(). fa_bifactor() passes its labels; fa_second_order()
# passes, for each second-order factor and the items of its first-order
# factors, the loadings of the Schmid-Leiman solution its parameters imply
# (.fa_second_order_syntax()). Returns the lines, omega, omega_H, ECV, and
# H of the general factor, then omega, omega_HS, ECV, and H of each group
# factor.
.fa_layered_index_lines <- function(general, groups, lg, ls, ps, group) {
  alone <- is.na(group)
  sq_sum <- function(l) paste0("(", paste(l, collapse = " + "), ")^2")
  plus <- function(x) paste(x, collapse = " + ")
  lines <- character(0)
  common_g <- sq_sum(lg)
  common_s <- vapply(groups, function(f) sq_sum(ls[group %in% f]),
                     character(1L))
  # Each item's model implied variance, and the variance of its other
  # sources beside the general factor and beside its group factor.
  v_i <- ifelse(alone, paste0("(", lg, "^2 + ", ps, ")"),
                paste0("(", lg, "^2 + ", ls, "^2 + ", ps, ")"))
  other_g <- ifelse(alone, ps, paste0("(", ls, "^2 + ", ps, ")"))
  other_s <- paste0("(", lg, "^2 + ", ps, ")")
  ecv <- function(r) {
    paste0("(", plus(paste0(lg[r], "^2/", v_i[r])), ") / (",
           plus(paste0(ifelse(alone[r], paste0(lg[r], "^2"),
                              paste0("(", lg[r], "^2 + ", ls[r], "^2)")),
                       "/", v_i[r])), ")")
  }
  H <- function(l, other) {
    paste0("1 / (1 + 1/(", plus(paste0(l, "^2/", other)), "))")
  }
  total <- paste0("(", plus(c(common_g, common_s, ps)), ")")
  lines <- c(lines,
    paste0("omega_", general, " := (", plus(c(common_g, common_s)), ") / ",
           total),
    paste0("omega_H_", general, " := ", common_g, " / ", total),
    paste0("ECV_", general, " := ", ecv(seq_along(lg))),
    paste0("H_", general, " := ", H(lg, other_g)))
  for (f in groups) {
    r <- which(group %in% f)
    cg <- sq_sum(lg[r])
    sub_total <- paste0("(", plus(c(cg, common_s[[f]], ps[r])), ")")
    lines <- c(lines,
      paste0("omega_", f, " := (", cg, " + ", common_s[[f]], ") / ",
             sub_total),
      paste0("omega_HS_", f, " := ", common_s[[f]], " / ", sub_total),
      paste0("ECV_", f, " := ", ecv(r)),
      paste0("H_", f, " := ", H(ls[r], other_s[r])))
  }
  lines
}

# .fa_bifactor_maximal(fit, latent, items)
#
# The maximal reliability of each factor of a bifactor fit, the squared
# multiple correlation of the factor with all of the items,
# l_f' Sigma^{-1} l_f with l_f the factor's loadings (zero off it) and
# Sigma = L L' + Psi the model implied covariance matrix of the items
# (factor variances 1, factors uncorrelated), which does not depend on the
# metric of the items and equals bifactor_indices()'s l' R^{-1} l on the
# standardized loadings. The standard error is .fa_layered_maximal()'s
# (`vcov_usable` as there). Returns its list.
.fa_bifactor_maximal <- function(fit, latent, items, vcov_usable = TRUE) {
  pt <- lavaan::parTable(fit)
  load <- which(pt$op == "=~" & pt$lhs %in% latent & pt$rhs %in% items)
  errs <- which(pt$op == "~~" & pt$lhs == pt$rhs & pt$lhs %in% items)
  cells <- cbind(match(pt$rhs[load], items), match(pt$lhs[load], latent))
  ei <- match(pt$lhs[errs], items)
  .fa_layered_maximal(fit, latent, function(e) {
    L <- matrix(0, length(items), length(latent))
    L[cells] <- e[load]
    psi <- numeric(length(items))
    psi[ei] <- e[errs]
    Sigma <- tcrossprod(L)
    diag(Sigma) <- diag(Sigma) + psi
    SiL <- tryCatch(solve(Sigma, L), error = function(e) NULL)
    if (is.null(SiL)) return(rep(NA_real_, length(latent)))
    colSums(L * SiL)
  }, vcov_usable)
}

# .fa_layered_maximal(fit, latent, rho, vcov_usable)
#
# The maximal reliabilities of the factors `latent` of a fit and their
# delta method standard errors. `rho` maps the estimates of every row of
# lavaan's parTable() to the maximal reliabilities, in the order of
# `latent` (NA when the model implied covariance matrix cannot be
# inverted); a matrix inverse is what lavaan's defined parameters cannot
# express. The derivatives are taken numerically (central differences)
# over lavaan's free parameters, with lavaan's covariance matrix of their
# estimates, as fa_exploratory() forms the standard errors of the
# determinacy and H (.efa_factor_quality(), R/fa_exploratory_internals.R);
# equal labels are separate free parameters there, tied by the
# constraints lavaan's covariance matrix carries. The standard error is NA
# when lavaan's covariance matrix cannot be used (`vcov_usable` FALSE: no
# standard errors asked for, none computed, or a singular information
# matrix; see .cfa_k_vcov_status()), as every other row's is, and when the
# quadratic form J V J' is not positive, which a covariance matrix that is
# not positive definite can give. Returns a list with `estimate` and
# `se`, each named by factor.
.fa_layered_maximal <- function(fit, latent, rho, vcov_usable = TRUE) {
  pt <- lavaan::parTable(fit)
  free <- pt$free
  est <- pt$est
  idx <- which(free > 0L)
  x0 <- numeric(max(c(0L, free)))
  x0[free[idx]] <- est[idx]
  rho_x <- function(x) {
    e <- est
    e[idx] <- x[free[idx]]
    rho(e)
  }
  estimate <- rho_x(x0)
  se <- rep(NA_real_, length(latent))
  V <- if (vcov_usable) {
    tryCatch(unclass(lavaan::lavInspect(fit, "vcov")),
             error = function(e) NULL)
  }
  if (!is.null(V) && length(V) && all(dim(V) == length(x0)) &&
      all(is.finite(estimate))) {
    J <- matrix(0, length(latent), length(x0))
    for (i in seq_along(x0)) {
      h <- 1e-6 * max(1, abs(x0[i]))
      up <- dn <- x0
      up[i] <- up[i] + h
      dn[i] <- dn[i] - h
      J[, i] <- (rho_x(up) - rho_x(dn)) / (2 * h)
    }
    q <- diag(J %*% V %*% t(J))
    se <- ifelse(is.finite(q) & q > 0, sqrt(pmax(q, 0)), NA_real_)
  }
  list(estimate = stats::setNames(estimate, latent),
       se = stats::setNames(se, latent))
}

# .fa_bifactor_rows(params, fit, pe, factors, general, items, conf_level,
#                   wald_df, vcov_usable)
#
# The parameter table of a bifactor fit (`params`, in fa_confirmatory()'s
# columns, the defined indices already carrying their logit-scale limits)
# with the rows lavaan cannot define inserted and the rows an improper
# solution leaves undefined set to NA (.fa_layered_rows()): the maximal
# reliability of every factor (.fa_bifactor_maximal(), its standard error,
# test, and limits NA when `vcov_usable` is FALSE, as every other row's
# are), and PUC_<general>,
# the proportion of item pairs that do not share a group factor
# (Rodriguez, Reise, and Haviland, 2016). `pe` is lavaan's
# parameterEstimates(), whose loadings and error variances give the
# standardized loadings the improper check reads. Returns
# .fa_layered_rows()'s list.
.fa_bifactor_rows <- function(params, fit, pe, factors, general, items,
                              conf_level, wald_df, vcov_usable = TRUE) {
  groups <- names(factors)
  latent <- c(general, groups)
  pairs <- function(k) k * (k - 1) / 2
  psi_rows <- pe$op == "~~" & pe$lhs == pe$rhs & pe$lhs %in% items
  psi <- stats::setNames(pe$est[psi_rows], pe$lhs[psi_rows])[items]
  load_rows <- pe[pe$op == "=~" & pe$lhs %in% latent, , drop = FALSE]
  L <- matrix(0, length(items), length(latent),
              dimnames = list(items, latent))
  L[cbind(match(load_rows$rhs, items), match(load_rows$lhs, latent))] <-
    load_rows$est
  .fa_layered_rows(
    params, latent = latent,
    maximal = .fa_bifactor_maximal(fit, latent, items, vcov_usable),
    puc = stats::setNames(
      1 - sum(pairs(lengths(factors))) / pairs(length(items)), general),
    L = L, psi = psi, generals = stats::setNames(list(items), general),
    groups = factors, conf_level = conf_level, wald_df = wald_df,
    closing = paste0(" A negative error variance at a maximum of the ",
                     "likelihood can come from sampling error around a ",
                     "small error variance, from a model that does not ",
                     "describe these items, or from a group factor the ",
                     "data barely separate from the general factor."))
}

# .fa_layered_rows(params, latent, maximal, puc, L, psi, generals, groups,
#                  conf_level, wald_df, closing)
#
# The table of a layered model (fa_bifactor(), fa_second_order()) with the
# rows lavaan cannot define inserted: maximal_reliability_<f> after each
# factor's H row (`maximal`, from .fa_layered_maximal(), named by factor),
# with its delta method standard error, its test_statistic and p_value on the
# reference that wald_df names, and limits formed on the logit scale, as
# the other indices' are; and PUC_<g> after the maximal reliability of
# each general factor g (`puc`, named by the general factors), which
# depends on the structure alone and so has no standard error or interval.
# `latent` gives the factors in the order of their rows, each general
# factor's followed by those of its group factors. An item with a negative
# error variance (`psi`, named by item) makes the solution improper: the
# maximal reliabilities are then not defined, H is not defined for a
# factor on which an item's standardized loading is at or beyond one (the
# index gives that item an error variance that is not positive; `L`, the
# items by `latent` matrix of the loadings of the orthogonal layered
# solution, whose rows sum their squares with `psi` to the item
# variances), and omega and omega_H of a general factor (omega and
# omega_HS of a group factor) are not proportions when the error variances
# of its items (`generals` and `groups`, lists of the items of each) sum
# below zero; those rows are NA, with a warning that names the items and
# ends with `closing`. Returns a list with `params`, `index_terms` (the
# index rows, in order), and `improper` (the negative error variances,
# named by item; empty for a proper solution).
.fa_layered_rows <- function(params, latent, maximal, puc, L, psi, generals,
                             groups, conf_level, wald_df, closing = NULL) {
  mr <- maximal
  z <- mr$estimate / mr$se
  p <- 2 * stats::pt(-abs(z), df = wald_df)
  maximal <- data.frame(
    syntax = "", term = paste0("maximal_reliability_", latent),
    estimate = unname(mr$estimate[latent]), se = unname(mr$se[latent]),
    test_statistic = unname(z[latent]), p_value = unname(p[latent]),
    ci_lower = NA_real_, ci_upper = NA_real_, stringsAsFactors = FALSE)
  maximal <- .cfa_k_logit_limits(maximal, maximal$term, conf_level, wald_df)
  puc_rows <- data.frame(
    syntax = "", term = paste0("PUC_", names(puc)), estimate = unname(puc),
    se = NA_real_, test_statistic = NA_real_, p_value = NA_real_,
    ci_lower = NA_real_, ci_upper = NA_real_, stringsAsFactors = FALSE)
  # Insert each factor's maximal reliability after its H row, and PUC after
  # each general factor's.
  general <- latent %in% names(puc)
  for (k in seq_along(latent)) {
    at <- match(paste0("H_", latent[k]), params$term)
    add <- if (general[k]) {
      rbind(maximal[k, ], puc_rows[match(latent[k], names(puc)), ])
    } else {
      maximal[k, ]
    }
    params <- rbind(params[seq_len(at), , drop = FALSE], add,
                    params[-seq_len(at), , drop = FALSE])
  }
  rownames(params) <- NULL
  index_terms <- unlist(lapply(seq_along(latent), function(k) {
    f <- latent[k]
    if (general[k]) {
      paste0(c("omega_", "omega_H_", "ECV_", "H_", "maximal_reliability_",
               "PUC_"), f)
    } else {
      paste0(c("omega_", "omega_HS_", "ECV_", "H_", "maximal_reliability_"),
             f)
    }
  }))

  # ---- an improper solution ---------------------------------------------
  improper <- psi[!is.na(psi) & psi < 0]
  if (length(improper)) {
    v <- rowSums(L^2) + psi
    std <- L / sqrt(pmax(v, 0))
    on <- L != 0
    h_na <- latent[vapply(latent, function(f) {
      s <- std[on[, f], f]
      any(!is.finite(s) | abs(s) >= 1)
    }, logical(1L))]
    omega_na <- c(
      unlist(lapply(names(generals), function(g) {
        if (sum(psi[generals[[g]]]) < 0) paste0(c("omega_", "omega_H_"), g)
      })),
      unlist(lapply(names(groups), function(f) {
        if (sum(psi[groups[[f]]]) < 0) paste0(c("omega_", "omega_HS_"), f)
      })))
    na_terms <- c(paste0("maximal_reliability_", latent),
                  paste0("H_", h_na), omega_na)
    na_cols <- c("estimate", "se", "test_statistic", "p_value", "ci_lower",
                 "ci_upper")
    params[params$term %in% na_terms, na_cols] <- NA_real_
    .dmar_warn(c(
      "Improper solution: the error variance of ",
      paste0(names(improper), " (",
             formatC(improper, format = "g", digits = 3), ")",
             collapse = ", "),
      if (length(improper) == 1L) " is" else " are", " negative. ",
      "The maximal reliabilities are not defined for an improper ",
      "solution",
      if (length(h_na)) {
        paste0(", nor is H for ", paste(h_na, collapse = ", "),
               ", on which a standardized loading is at or beyond one")
      },
      if (length(omega_na)) {
        paste0(", nor ", paste(omega_na, collapse = ", "),
               ", whose items' error variances sum below zero")
      },
      "; those rows are NA.", closing), "result_qualified")
  }
  list(params = params, index_terms = index_terms, improper = improper)
}

# .fa_bifactor_attributes(out, factors, general, general_only, eq_load,
#                         eq_err, eq_int, has_means, improper)
#
# Records the bifactor structure on a fa_bifactor() table: the general
# factor's name ("general"), the group factors ("factors"), the items on
# the general factor alone ("general_only"), the factors whose loadings,
# error variances, and intercepts are equated ("equal", a list), whether
# the mean structure is modeled ("has_means"), the metric of the omega
# rows ("omega_metric", "raw_item_sum"), and the negative error variances
# of an improper solution ("improper", empty for a proper one). The header
# printed above the table is built from them
# (.dmar_tbl_footer.dmar_fa_bifactor(), R/fa_bifactor.R).
.fa_bifactor_attributes <- function(out, factors, general, general_only,
                                    eq_load, eq_err, eq_int, has_means,
                                    improper) {
  attr(out, "general") <- general
  attr(out, "factors") <- factors
  attr(out, "general_only") <- as.character(general_only)
  attr(out, "equal") <- list(loading = names(eq_load)[eq_load],
                             error = names(eq_err)[eq_err],
                             intercept = if (has_means) {
                               names(eq_int)[eq_int]
                             } else {
                               character(0)
                             })
  attr(out, "has_means") <- has_means
  attr(out, "omega_metric") <- "raw_item_sum"
  attr(out, "improper") <- improper
  out
}

# ---------------------------------------------------------------------------
# The second-order structure (fa_second_order())
# ---------------------------------------------------------------------------
#
# First-order factors, each measured by its own items (simple structure),
# and second-order factors, each measured by its own first-order factors,
# correlated with each other when there are several. `factors` names the
# first-order factors and their items, and `second_order` names the
# second-order factors and their first-order factors. The model is
# identified as lavaan's std.lv identifies it: the variance of every
# second-order factor is fixed at 1, and so is the variance of the
# disturbance of every first-order factor, the part of the factor its
# second-order factor does not account for, so every first-order and
# second-order loading is free; both constraints are linear. A first-order
# factor f under the second-order factor g is then f = gamma_f g + u_f,
# with Var(u_f) = 1, so its variance is 1 + gamma_f^2 and its correlation
# with g, the standardized second-order loading, gamma_f / sqrt(1 +
# gamma_f^2); and an item i of f is
#   x_i = nu_i + lambda_i gamma_f g + lambda_i u_f + e_i,
# the orthogonal layered form of Schmid and Leiman (1957), with general
# loading lambda_i gamma_f and group loading lambda_i (Yung, Thissen, and
# McLeod, 1999, Eq. 8 and Eq. 10 for two levels: B_0 = [P_1 P_2 | P_1
# U_1], with U_1 = diag(1 - P_2 P_2')^(1/2) by their Eq. 4). In this
# parameterization the Schmid-Leiman loadings are products of the free
# parameters, and over the square root of the item's model implied
# variance v_i = lambda_i^2 (1 + gamma_f^2) + psi_i they are the
# standardized loadings lambda_std gamma_std and lambda_std sqrt(1 -
# gamma_std^2). The indices of each second-order factor and its
# first-order factors are those of a bifactor model with these loadings
# over the items of its first-order factors (.fa_layered_index_lines()).

# .fa_second_order_refuse(message, reason, factors, items)
#
# The classed refusal of a second-order structure fa_second_order() will
# not fit: an error of class c("dmar_second_order_unidentified",
# "dmar_refusal", "error", "condition") with code
# "second_order_unidentified", the factors and items at issue, and
# `reason`, as for .fa_bifactor_refuse(): "unidentified", "names" (a name
# placed twice or not at all, or names that clash), or "scope" (an item on
# two first-order factors, or a first-order factor under two second-order
# factors, which the model of fa_second_order() does not have), all in
# the condition (.dmar_refuse(), R/dmar_conditions.R).
.fa_second_order_refuse <- function(message, reason, factors = character(0),
                                    items = character(0)) {
  .dmar_refuse(message, "second_order_unidentified", reason = reason,
               factors = factors, items = items)
}

# .fa_second_order_items(factors, second_order)
#
# Checks the names of a second-order structure, refusing
# (.fa_second_order_refuse()) a first-order and a second-order factor of
# one name, an item listed twice under one first-order factor, a factor
# named like an item, an item on two first-order factors, a second-order
# factor that names a first-order factor `factors` does not have, and a
# first-order factor under no second-order factor or under two. Returns
# the items, in the order of `factors`.
.fa_second_order_items <- function(factors, second_order) {
  first <- names(factors)
  for (f in first) {
    if (!is.character(factors[[f]]) || !length(factors[[f]]) ||
        anyNA(factors[[f]])) {
      .dmar_refuse_argument(
        sprintf("First-order factor '%s' must name its items.", f), arg = NULL)
    }
  }
  if (!is.list(second_order) || !length(second_order) ||
      is.null(names(second_order)) || anyNA(names(second_order)) ||
      any(names(second_order) == "")) {
    .dmar_refuse_argument(c(
      "'second_order' must be a named list: ",
      "list(second_order_factor = c(\"first_order_factor_1\", ",
      "\"first_order_factor_2\", ...))."), arg = "second_order")
  }
  second <- names(second_order)
  if (anyDuplicated(second)) {
    .dmar_refuse_argument(
      "Second-order factor names must be unique.", arg = NULL)
  }
  if (!all(second == make.names(second))) {
    .dmar_refuse_argument(c(
      "Second-order factor names must be syntactically valid R names ",
      "(they become lavaan labels); rename, for example, 'my factor' ",
      "to 'my_factor'."), arg = NULL)
  }
  for (g in second) {
    if (!is.character(second_order[[g]]) || !length(second_order[[g]]) ||
        anyNA(second_order[[g]])) {
      .dmar_refuse_argument(c(
        sprintf("Second-order factor '%s' must name its first-order ",
                g), "factors."), arg = NULL)
    }
  }
  both <- intersect(second, first)
  if (length(both)) {
    .fa_second_order_refuse(c(
      "A first-order and a second-order factor are both named ",
      paste0("'", both, "'", collapse = ", "), ". Rename one of them."),
      "names", factors = both)
  }
  # An item listed twice under one first-order factor is a slip in the
  # list, not an item on two factors, and is named as such.
  repeated <- lapply(factors, function(v) unique(v[duplicated(v)]))
  repeated <- repeated[lengths(repeated) > 0L]
  if (length(repeated)) {
    .fa_second_order_refuse(c(
      "Each item is listed once under its first-order factor, but ",
      paste0(vapply(repeated, paste, character(1L), collapse = ", "),
             ifelse(lengths(repeated) == 1L, " is", " are"),
             " listed more than once under ", names(repeated),
             collapse = "; "),
      ". List each item once."),
      "names", factors = names(repeated),
      items = unlist(repeated, use.names = FALSE))
  }
  items <- unlist(factors, use.names = FALSE)
  several <- unique(items[duplicated(items)])
  if (length(several)) {
    where <- vapply(several, function(i) {
      paste(first[vapply(factors, function(v) i %in% v, logical(1L))],
            collapse = " and ")
    }, character(1L))
    .fa_second_order_refuse(c(
      "Each item loads on one first-order factor in a second-order model, ",
      "but ", paste0(several, " is listed under ", where, collapse = "; "),
      ". List each item under one first-order factor."),
      "scope", factors = first, items = several)
  }
  clash <- intersect(c(first, second), items)
  if (length(clash)) {
    .fa_second_order_refuse(c(
      "A factor and an item are both named ",
      paste0("'", clash, "'", collapse = ", "),
      ". Rename the factor, since lavaan reads the name as either one."),
      "names", factors = clash, items = clash)
  }
  listed <- unlist(second_order, use.names = FALSE)
  unknown <- setdiff(listed, first)
  if (length(unknown)) {
    .fa_second_order_refuse(c(
      "'second_order' names ", paste0("'", unknown, "'", collapse = ", "),
      ", which ", if (length(unknown) == 1L) "is" else "are",
      " not a first-order factor of 'factors' (",
      paste0("'", first, "'", collapse = ", "), "). A second-order factor ",
      "is measured by first-order factors, named as 'factors' names them."),
      "names", factors = unknown)
  }
  twice <- unique(listed[duplicated(listed)])
  if (length(twice)) {
    owners <- lapply(twice, function(f) {
      second[vapply(second_order, function(v) f %in% v, logical(1L))]
    })
    where <- vapply(owners, function(o) {
      if (length(o) == 1L) paste("twice under", o) else
        paste("under", paste(o, collapse = " and "))
    }, character(1L))
    # Under one second-order factor twice is a slip in the list; under two
    # it is a structure the model does not have.
    .fa_second_order_refuse(c(
      "Each first-order factor is measured by one second-order factor, ",
      "but ", paste0("'", twice, "' is listed ", where, collapse = "; "),
      ". List each first-order factor under one second-order factor."),
      if (any(lengths(owners) > 1L)) "scope" else "names", factors = twice)
  }
  none <- setdiff(first, listed)
  if (length(none)) {
    .fa_second_order_refuse(c(
      "Every first-order factor is measured by a second-order factor, but ",
      paste0("'", none, "'", collapse = ", "),
      if (length(none) == 1L) " is" else " are",
      " listed under none in 'second_order'. List ",
      if (length(none) == 1L) "it" else "each",
      " under a second-order factor, or leave ",
      if (length(none) == 1L) "it" else "them",
      " out of 'factors'."), "names", factors = none)
  }
  items
}

# .fa_second_order_identification(factors, second_order, eq_gamma)
#
# Refuses (.fa_second_order_refuse()) a first-order factor of one item and
# the second-order factors the first-order factors cannot identify. A
# first-order factor of one item cannot be told apart from the item: the
# item's covariances with the other items fix the product of its loading
# and its factor's second-order loading, and its variance fixes only the
# sum of its error variance and the variance its loading carries. A
# second-order factor of one first-order factor cannot be told apart from
# that factor, whose variance it shares with the disturbance. The
# second-order loadings of g are fixed by the correlations of its
# first-order factors, each the product of two standardized loadings (and
# of the correlation of the second-order factors, for two first-order
# factors under different ones): three first-order factors give three
# correlations for three loadings, and two under one second-order factor
# give one correlation for two loadings, so two are identified only with
# equal loadings (one correlation for one loading) or with another
# second-order factor beside them, whose first-order factors' correlations
# with these fix each loading on its own (when the second-order factors
# correlate).
.fa_second_order_identification <- function(factors, second_order,
                                            eq_gamma) {
  J <- lengths(factors)
  one <- names(factors)[J < 2L]
  if (length(one)) {
    .fa_second_order_refuse(c(
      "A first-order factor needs two or more items: ",
      paste0("'", one, "'", collapse = ", "), " ",
      if (length(one) == 1L) "has" else "have",
      " one item, whose loading cannot be told apart from its error ",
      "variance. Give the factor more items, or leave it out."),
      "unidentified", factors = one,
      items = unlist(factors[one], use.names = FALSE))
  }
  K <- lengths(second_order)
  single <- names(second_order)[K < 2L]
  if (length(single)) {
    .fa_second_order_refuse(c(
      "A second-order factor needs two or more first-order factors: ",
      paste0("'", single, "' has one ('",
             vapply(second_order[single], `[`, character(1L), 1L), "')",
             collapse = "; "),
      ", and a second-order factor of one first-order factor cannot be ",
      "told apart from it, since the two share its variance with the ",
      "disturbance. Fit that factor as a first-order factor with ",
      "fa_confirmatory(), or give the second-order factor more first-order ",
      "factors."), "unidentified", factors = single)
  }
  if (length(second_order) == 1L && K[[1L]] == 2L && !eq_gamma[[1L]]) {
    g <- names(second_order)
    .fa_second_order_refuse(c(
      "With one second-order factor, two first-order factors (",
      paste0("'", second_order[[1L]], "'", collapse = " and "),
      ") do not identify its loadings: their one correlation is the ",
      "product of the two standardized loadings. ",
      "equal_second_order_loading = TRUE identifies them (the two factors ",
      "then correlate equally with '", g, "'), and so does a third ",
      "first-order factor; fa_confirmatory() fits the correlation of the ",
      "two factors, which is what the data identify."),
      "unidentified", factors = c(g, second_order[[1L]]))
  }
  invisible(TRUE)
}

# .fa_second_order_unlinked(data, factors, second_order, eq_gamma)
#
# Refuses (.fa_second_order_refuse(), reason "unidentified") incomplete
# data whose pairs of first-order factors observed together do not
# identify the second-order structure. A pair no case observes together
# (no case observes an item of one with an item of the other) gives the
# likelihood nothing about the correlation of the two factors, which the
# second-order model writes as the product of their standardized
# second-order loadings (and of the correlation of their second-order
# factors, for factors under different ones). The products of the pairs
# that are observed identify the loadings and the correlations when the
# matrix of their derivatives has full column rank, which is checked at a
# point chosen to lie on no special set (standardized loadings spread over
# .5 to .85, correlations of the second-order factors over .3 to .5), as
# the rank at almost every point is the rank the model has. Three
# first-order factors under one second-order factor, one pair of them
# never observed together, leave two products for three loadings and are
# refused; four, with one such pair, leave five products for four loadings
# and are fit. fa_confirmatory() refuses every such pair, whose correlation
# is a parameter of its model (.cfa_k_unlinked_factors()). Through the
# first version of fa_second_order() these data were fit, and the loadings
# lavaan left at arbitrary values were reported with the indices built on
# them.
.fa_second_order_unlinked <- function(data, factors, second_order,
                                      eq_gamma) {
  first <- names(factors)
  k <- length(first)
  observed <- !is.na(as.matrix(data[, unlist(factors, use.names = FALSE),
                                    drop = FALSE]))
  together <- crossprod(observed)
  pairs <- utils::combn(k, 2L)
  linked <- apply(pairs, 2L, function(p) {
    any(together[factors[[p[1L]]], factors[[p[2L]]]] > 0)
  })
  if (all(linked)) return(invisible(TRUE))
  second <- names(second_order)
  owner <- stats::setNames(rep(second, lengths(second_order)),
                           unlist(second_order, use.names = FALSE))[first]
  # One parameter per standardized second-order loading (one per
  # second-order factor whose loadings are equated) and per correlation of
  # two second-order factors.
  load_par <- ifelse(eq_gamma[owner], paste0("g:", owner), paste0("f:", first))
  par <- unique(load_par)
  phi_pairs <- if (length(second) > 1L) utils::combn(second, 2L)
  phi_par <- if (!is.null(phi_pairs)) {
    paste0("p:", phi_pairs[1L, ], ":", phi_pairs[2L, ])
  }
  value <- c(stats::setNames(0.5 + 0.35 * ((seq_along(par) * 0.618034) %% 1),
                             par),
             if (length(phi_par)) {
               stats::setNames(0.3 + 0.2 * ((seq_along(phi_par) * 0.414214) %%
                                              1), phi_par)
             })
  phi_of <- function(g, h) {
    if (g == h) return(NA_character_)
    o <- sort(c(g, h))
    nm <- paste0("p:", o[1L], ":", o[2L])
    if (nm %in% phi_par) nm else paste0("p:", o[2L], ":", o[1L])
  }
  obs <- pairs[, linked, drop = FALSE]
  J <- matrix(0, ncol(obs), length(value), dimnames = list(NULL, names(value)))
  for (r in seq_len(ncol(obs))) {
    i <- obs[1L, r]
    j <- obs[2L, r]
    pi <- load_par[i]
    pj <- load_par[j]
    ph <- phi_of(owner[[i]], owner[[j]])
    p_val <- if (is.na(ph)) 1 else value[[ph]]
    J[r, pi] <- J[r, pi] + value[[pj]] * p_val
    J[r, pj] <- J[r, pj] + value[[pi]] * p_val
    if (!is.na(ph)) J[r, ph] <- J[r, ph] + value[[pi]] * value[[pj]]
  }
  if (qr(J, tol = 1e-7)$rank == ncol(J)) return(invisible(TRUE))
  never <- pairs[, !linked, drop = FALSE]
  named <- paste0("'", first[never[1L, ]], "' with an item of '",
                  first[never[2L, ]], "'")
  .fa_second_order_refuse(c(
    "No case observes an item of ", paste(named, collapse = ", or an item of "),
    ", and the pairs of first-order factors the data do observe together ",
    "do not identify the second-order structure: under the model the ",
    "correlation of two first-order factors is the product of their ",
    "standardized second-order loadings (and of the correlation of their ",
    "second-order factors), and a pair never observed together enters the ",
    "likelihood only through covariances the data do not have, so lavaan ",
    "would leave the loadings at arbitrary values. Link the forms with ",
    "items given on both, or give each second-order factor more ",
    "first-order factors observed together; missing_coverage() reports ",
    "the cases behind every pair of items."),
    "unidentified", factors = unique(first[never]))
}

# .fa_second_order_structural_df(second_order, eq_gamma)
#
# The degrees of freedom of the structural part of a second-order model:
# the k(k - 1)/2 correlations of its k first-order factors less the
# second-order loadings (one for a second-order factor whose loadings are
# equated) and the correlations of the second-order factors. Zero means
# the structural part is just identified: the second-order loadings are
# then another parameterization of the correlations of the first-order
# factors, and the model has the fit of the correlated factors model of
# the same first-order factors whenever those correlations are ones the
# second-order loadings can reproduce (one second-order factor over three
# first-order factors, or over two with equal loadings).
.fa_second_order_structural_df <- function(second_order, eq_gamma) {
  K <- lengths(second_order)
  k <- sum(K)
  m <- length(second_order)
  k * (k - 1) / 2 - sum(ifelse(eq_gamma[names(second_order)], 1, K)) -
    m * (m - 1) / 2
}

# .fa_second_order_labels(factors, second_order, eq_load, eq_int, eq_err,
#                         eq_gamma)
#
# The lavaan labels of the second-order model. `items`: one row per item in
# the order of `factors`, with its first-order factor, that factor's
# second-order factor, and the labels of its loading, error variance, and
# intercept (lambda_<f>_<j>, psi_<f>_<j>, nu_<f>_<j>, numbered within the
# factor; lambda_<f>, psi_<f>, nu_<f> when equated). `gammas`: one row per
# first-order factor in the order of `second_order`, with its second-order
# factor and the label of its second-order loading (gamma_<g>_<f>, or
# gamma_<g> when the loadings of g are equated).
.fa_second_order_labels <- function(factors, second_order, eq_load, eq_int,
                                    eq_err, eq_gamma) {
  lab <- function(prefix, f, j, equal) {
    if (equal) rep(paste0(prefix, "_", f), length(j)) else
      paste0(prefix, "_", f, "_", j)
  }
  owner <- stats::setNames(rep(names(second_order), lengths(second_order)),
                           unlist(second_order, use.names = FALSE))
  items <- do.call(rbind, lapply(names(factors), function(f) {
    j <- seq_along(factors[[f]])
    data.frame(item = factors[[f]], factor = f, second = owner[[f]],
               loading = lab("lambda", f, j, eq_load[[f]]),
               error = lab("psi", f, j, eq_err[[f]]),
               intercept = lab("nu", f, j, eq_int[[f]]),
               stringsAsFactors = FALSE)
  }))
  gammas <- do.call(rbind, lapply(names(second_order), function(g) {
    fs <- second_order[[g]]
    data.frame(factor = fs, second = g,
               gamma = if (eq_gamma[[g]]) rep(paste0("gamma_", g), length(fs))
                       else paste0("gamma_", g, "_", fs),
               stringsAsFactors = FALSE)
  }))
  rownames(items) <- rownames(gammas) <- NULL
  list(items = items, gammas = gammas)
}

# .fa_second_order_index_terms(second_order)
#
# The indices lavaan defines for a second-order model, in the order of the
# table: for each second-order factor, its omega, omega_H, ECV, and H, then
# for each of its first-order factors omega, omega_HS, ECV, and H.
.fa_second_order_index_terms <- function(second_order) {
  unlist(lapply(names(second_order), function(g) {
    .fa_bifactor_index_terms(g, second_order[[g]])
  }), use.names = FALSE)
}

# .fa_second_order_syntax(factors, second_order, eq_load, eq_int, eq_err,
#                         eq_gamma, has_means)
#
# The lavaan syntax of the second-order model: each first-order factor on
# its own items and each second-order factor on its own first-order
# factors, the first loading of each freed; the variance of every
# second-order factor and of every first-order factor's disturbance fixed
# at 1 (and the means of both at 0 with a mean structure); the
# second-order factors' correlations labeled phi_<g1>_<g2>; every loading,
# error variance, and intercept labeled (.fa_second_order_labels()); and
# the indices of each second-order factor and its first-order factors
# defined from the labels (.fa_layered_index_lines()), with the general
# loading of an item of f written (lambda_i*gamma_f) and its group
# loading lambda_i, the Schmid-Leiman loadings the parameters imply. With
# the disturbance's variance at 1 the variance of the sum of the items of
# g is (sum_i lambda_i gamma_f(i))^2 + sum_f (sum_{i in f} lambda_i)^2 +
# sum_i psi_i, and omega_H_<g> is the first term over it, the share of
# that variance due to g. Returns one string.
.fa_second_order_syntax <- function(factors, second_order, eq_load, eq_int,
                                    eq_err, eq_gamma, has_means) {
  lab <- .fa_second_order_labels(factors, second_order, eq_load, eq_int,
                                 eq_err, eq_gamma)
  it <- lab$items
  ga <- lab$gammas
  first <- names(factors)
  second <- names(second_order)
  items <- it$item
  pairs <- if (length(second) > 1L) utils::combn(second, 2L)
  phi <- if (!is.null(pairs)) {
    paste0("phi_", pairs[1L, ], "_", pairs[2L, ])
  }
  # Labels of different parameters that coincide (a first-order factor
  # named "a_1" with equal loadings beside a factor "a" would label both
  # lambda_a_1), or that coincide with the names of other rows, an index,
  # a variable, or a factor, would be read by lavaan or the reader as one:
  # refused.
  label <- c(it$loading, it$error, if (has_means) it$intercept, ga$gamma,
             phi)
  owner <- c(paste("lambda", it$factor), paste("psi", it$factor),
             if (has_means) paste("nu", it$factor), paste("gamma", ga$second),
             phi)
  shared <- tapply(owner, label, function(o) length(unique(o)))
  rows <- c(paste0("zeta_", first), paste0("phi_", second))
  index_terms <- .fa_second_order_index_terms(second_order)
  names_all <- c(items, first, second)
  clash <- unique(c(names(shared)[shared > 1L],
                    intersect(c(label, rows), c(index_terms, names_all)),
                    intersect(rows, label), rows[duplicated(rows)],
                    index_terms[duplicated(index_terms)],
                    intersect(index_terms, names_all)))
  if (length(clash)) {
    .fa_second_order_refuse(c(
      "The names of the factors make labels that coincide (",
      paste(clash, collapse = ", "), "), which lavaan would read as one ",
      "parameter. Rename the factors."), "names", factors = c(first, second))
  }
  first_free <- function(f, v, l) {
    paste0(f, " =~ NA*", v[1L], " + ", paste0(l, "*", v, collapse = " + "))
  }
  lines <- character(0)
  for (f in first) {
    r <- which(it$factor == f)
    v <- items[r]
    lines <- c(lines, first_free(f, v, it$loading[r]),
               paste0(f, " ~~ 1*", f),
               paste0(v, " ~~ ", it$error[r], "*", v))
    if (has_means) {
      lines <- c(lines, paste0(v, " ~ ", it$intercept[r], "*1"),
                 paste0(f, " ~ 0*1"))
    }
  }
  for (g in second) {
    r <- which(ga$second == g)
    lines <- c(lines, first_free(g, ga$factor[r], ga$gamma[r]),
               paste0(g, " ~~ 1*", g))
    if (has_means) lines <- c(lines, paste0(g, " ~ 0*1"))
  }
  if (!is.null(pairs)) {
    lines <- c(lines, paste0(pairs[1L, ], " ~~ ", phi, "*", pairs[2L, ]))
  }
  # The indices of each second-order factor over the items of its
  # first-order factors, taken in the order 'second_order' lists them.
  for (g in second) {
    fs <- second_order[[g]]
    r <- unlist(lapply(fs, function(f) which(it$factor == f)),
                use.names = FALSE)
    gam <- ga$gamma[match(it$factor[r], ga$factor)]
    lines <- c(lines, .fa_layered_index_lines(
      g, fs, lg = paste0("(", it$loading[r], "*", gam, ")"),
      ls = it$loading[r], ps = it$error[r], group = it$factor[r]))
  }
  paste(lines, collapse = "\n")
}

# .fa_second_order_check_syntax(factors, second_order, eq_load, eq_int,
#                               eq_err, eq_gamma, has_means)
#
# The second-order model written so that its boundary is a point of the
# parameter space, for .fa_second_order_maximum(). fa_second_order() fixes
# the variance of each first-order factor's disturbance at 1, so its
# standardized second-order loading gamma_f / sqrt(1 + gamma_f^2) reaches
# one only as gamma_f runs off to infinity, the first-order loadings
# lambda_i running to zero beside it. Here each first-order factor is
# written on the scale of its second-order factor instead, f = g + u_f,
# with its loading on g fixed at 1, the variance of u_f, z_f = 1 /
# gamma_f^2, free and bounded below at 0, and each item loading
# a_i = lambda_i gamma_f on f (its loading on g in the Schmid-Leiman
# solution). A standardized loading of one is then z_f = 0, a point the
# bounded optimizer reaches and stops at. The two are the same model: a
# point of fa_second_order()'s parameterization with no gamma_f at 0 is a
# point here with the same likelihood, and this one adds the boundary. The
# error variances, intercepts, and correlations of the second-order factors
# keep their labels, a_i takes lambda_i's and z_f takes gamma_f's, so the
# equality constraints carry over: equal first-order loadings are equal
# a_i, and equal second-order loadings are equal z_f (the sign that
# gamma_f gives up is carried by the a_i of its items). Returns one string.
.fa_second_order_check_syntax <- function(factors, second_order, eq_load,
                                          eq_int, eq_err, eq_gamma,
                                          has_means) {
  lab <- .fa_second_order_labels(factors, second_order, eq_load, eq_int,
                                 eq_err, eq_gamma)
  it <- lab$items
  ga <- lab$gammas
  lines <- character(0)
  for (f in names(factors)) {
    r <- which(it$factor == f)
    v <- it$item[r]
    z <- ga$gamma[ga$factor == f]
    lines <- c(lines,
               paste0(f, " =~ NA*", v[1L], " + ",
                      paste0(it$loading[r], "*", v, collapse = " + ")),
               paste0(f, " ~~ ", z, "*", f, " + lower(0)*", f),
               paste0(v, " ~~ ", it$error[r], "*", v))
    if (has_means) {
      lines <- c(lines, paste0(v, " ~ ", it$intercept[r], "*1"),
                 paste0(f, " ~ 0*1"))
    }
  }
  for (g in names(second_order)) {
    lines <- c(lines,
               paste0(g, " =~ ", paste0("1*", second_order[[g]],
                                        collapse = " + ")),
               paste0(g, " ~~ 1*", g))
    if (has_means) lines <- c(lines, paste0(g, " ~ 0*1"))
  }
  if (length(second_order) > 1L) {
    pairs <- utils::combn(names(second_order), 2L)
    lines <- c(lines, paste0(pairs[1L, ], " ~~ phi_", pairs[1L, ], "_",
                             pairs[2L, ], "*", pairs[2L, ]))
  }
  paste(lines, collapse = "\n")
}

# .fa_second_order_point(pt, factors, second_order)
#
# A point of either parameterization of the second-order model, read from
# lavaan's parTable() `pt`: `load`, the first-order loading of each item
# (lambda_i, or a_i in the parameterization of
# .fa_second_order_check_syntax()), named by item; `second`, the
# second-order loading of each first-order factor (gamma_f, or z_f), named
# by factor; and `rest`, the other free rows (error variances, intercepts,
# and the correlations of the second-order factors), whose meaning the two
# share, as lhs, op, rhs, and est.
.fa_second_order_point <- function(pt, factors, second_order) {
  first <- names(factors)
  items <- unlist(factors, use.names = FALSE)
  ld <- pt$op == "=~" & pt$lhs %in% first & pt$rhs %in% items
  gm_std <- pt$op == "=~" & pt$lhs %in% names(second_order) &
    pt$rhs %in% first & pt$free > 0L
  gm_chk <- pt$op == "~~" & pt$lhs == pt$rhs & pt$lhs %in% first &
    pt$free > 0L
  second <- if (any(gm_std)) {
    stats::setNames(pt$est[gm_std], pt$rhs[gm_std])[first]
  } else {
    stats::setNames(pt$est[gm_chk], pt$lhs[gm_chk])[first]
  }
  rest <- pt$free > 0L & !ld & !gm_std & !gm_chk
  list(load = stats::setNames(pt$est[ld], pt$rhs[ld])[items],
       second = second,
       rest = data.frame(lhs = pt$lhs[rest], op = pt$op[rest],
                         rhs = pt$rhs[rest], est = pt$est[rest],
                         stringsAsFactors = FALSE))
}

# .fa_second_order_start(load, second, rest, factors, second_order,
#                        check = FALSE)
#
# The starting values that put lavaan at a point of the second-order model:
# a data frame of lhs, op, rhs, and est for the free rows of
# fa_second_order()'s syntax (`check` FALSE: `load` the lambda_i and
# `second` the gamma_f) or of .fa_second_order_check_syntax()'s (`check`
# TRUE: `load` the a_i and `second` the z_f), which lavaan matches to the
# rows of the model by lhs, op, and rhs.
.fa_second_order_start <- function(load, second, rest, factors,
                                   second_order, check = FALSE) {
  first <- names(factors)
  items <- unlist(factors, use.names = FALSE)
  f_of <- rep(first, lengths(factors))
  owner <- stats::setNames(rep(names(second_order), lengths(second_order)),
                           unlist(second_order, use.names = FALSE))
  rbind(
    data.frame(lhs = f_of, op = "=~", rhs = items,
               est = unname(load[items]), stringsAsFactors = FALSE),
    if (check) {
      data.frame(lhs = first, op = "~~", rhs = first,
                 est = unname(second[first]), stringsAsFactors = FALSE)
    } else {
      data.frame(lhs = unname(owner[first]), op = "=~", rhs = first,
                 est = unname(second[first]), stringsAsFactors = FALSE)
    },
    rest)
}

# .fa_second_order_orient_point(lambda, gamma, rest, factors, second_order,
#                               eq_gamma)
#
# A point of fa_second_order()'s parameterization turned to the reported
# direction, by reflections, each of which leaves the likelihood unchanged:
# first each first-order factor whose loadings sum below zero (its lambda_i
# and its gamma_f change sign together), then each second-order factor
# whose loadings sum below zero (its gamma_f and its correlations with the
# other second-order factors change sign). A first-order factor under a
# second-order factor with equal loadings is not reflected alone, since
# that would break the equality; it keeps the direction the common loading
# gives it. Returns the list of the three, with `reflected`, the factors
# reflected.
.fa_second_order_orient_point <- function(lambda, gamma, rest, factors,
                                          second_order, eq_gamma) {
  owner <- stats::setNames(rep(names(second_order), lengths(second_order)),
                           unlist(second_order, use.names = FALSE))
  reflected <- character(0)
  for (f in names(factors)) {
    if (eq_gamma[[owner[[f]]]]) next
    if (sum(lambda[factors[[f]]]) < 0) {
      lambda[factors[[f]]] <- -lambda[factors[[f]]]
      gamma[[f]] <- -gamma[[f]]
      reflected <- c(reflected, f)
    }
  }
  for (g in names(second_order)) {
    fs <- second_order[[g]]
    if (sum(gamma[fs]) < 0) {
      gamma[fs] <- -gamma[fs]
      r <- rest$op == "~~" & rest$lhs != rest$rhs &
        xor(rest$lhs == g, rest$rhs == g)
      rest$est[r] <- -rest$est[r]
      reflected <- c(reflected, g)
    }
  }
  list(lambda = lambda, gamma = gamma, rest = rest, reflected = reflected)
}

# .fa_second_order_refit(fit_args, start, how, model = fit_args$model)
#
# The model `model` (fa_second_order()'s syntax by default) fit with the
# arguments of the first fit (`fit_args`) from the point `start`
# (.fa_second_order_start()): `how` "polish" lets nlminb run from there,
# with no further attempt from lavaan's own starting values should it
# fail; "none" computes the fit at that point without moving
# (optim.method = "none": the estimates, standard errors, and test are
# those of the point). Under "none" a request for bootstrap standard
# errors, which would resample without fitting, is replaced by none. The
# warnings lavaan gives are collected rather than shown. Returns a list
# with `fit` (NULL when lavaan stopped with an error or, under "polish",
# did not converge) and `warnings`.
.fa_second_order_refit <- function(fit_args, start, how,
                                   model = fit_args$model) {
  args <- fit_args
  args$model <- model
  args$start <- start
  # lavaan would otherwise move starting values whose latent covariance
  # matrix is not positive definite (a correlation of the second-order
  # factors beyond one, which a maximum can have) to ones that are.
  args$check.start <- FALSE
  if (identical(how, "polish")) {
    args$optim.method <- "nlminb"
    args$optim.attempts <- 1L
  } else {
    args$optim.method <- "none"
    if (identical(args$se, "bootstrap")) args$se <- "none"
  }
  seen <- new.env(parent = emptyenv())
  seen$warnings <- list()
  fit <- tryCatch(withCallingHandlers(
    do.call(lavaan::lavaan, args),
    warning = function(w) {
      seen$warnings <- c(seen$warnings, list(w))
      invokeRestart("muffleWarning")
    }), error = function(e) NULL)
  if (!is.null(fit) && identical(how, "polish") &&
      !isTRUE(lavaan::lavInspect(fit, "converged"))) {
    fit <- NULL
  }
  list(fit = fit, warnings = seen$warnings)
}

# .fa_second_order_check_fit(args)
#
# One fit of .fa_second_order_maximum()'s check, lavaan::lavaan() called
# with `args`: the fit, or NULL when lavaan stopped with an error or the
# fit did not converge. Its warnings (a covariance matrix of the estimates
# that cannot be formed, for one) are about a fit that is not reported,
# and are muffled.
.fa_second_order_check_fit <- function(args) {
  chk <- tryCatch(suppressWarnings(do.call(lavaan::lavaan, args)),
                  error = function(e) NULL)
  if (is.null(chk) || !isTRUE(lavaan::lavInspect(chk, "converged"))) {
    return(NULL)
  }
  chk
}

# .fa_second_order_maximum(fit, fit_args, model_check, factors,
#                          second_order, eq_gamma)
#
# Finds the maximum of the likelihood of a second-order model, or says
# that it could not. `fit` is the fit of fa_second_order()'s syntax by
# nlminb (`fit_args`, its arguments), converged or not. In that
# parameterization a first-order factor whose standardized second-order
# loading is heading to one has gamma_f running off along a ridge on which
# the likelihood still rises, and nlminb can stop on it far from the
# maximum, with other second-order loadings drifted to values from 30 to
# 230 (on holzinger_swineford, one second-order factor over spatial,
# verbal, speed, and t24, t25, and t26 stopped at a chi square of 251.77
# whose maximum is 169.77), or not converge at all when the maximum is on
# the boundary. The model is therefore also fit in the parameterization of
# .fa_second_order_check_syntax(), whose boundary z_f = 0 the bounded
# optimizer reaches, from the point of the first fit and, when the first
# fit did not converge, has a factor at the boundary, or falls short, from
# lavaan's own starting values; the larger log likelihood is kept. Its
# equal labels are imposed with ceq.simple = TRUE, as one parameter each,
# since lavaan 0.7.3 left them unenforced beside the bound. lavaan
# accepts either fit as converged only when its gradient is near zero
# (every element within 0.001, a parameter at its bound aside). A
# first-order factor is at the boundary when z_f is at most 1e-6 there (a
# standardized loading of .9999995 or more; at an active bound nlminb
# returns z_f = 0). Two log likelihoods are taken to be the same maximum
# when they are within 1e-6 per case (`tol`, 1e-6 times the number of
# cases): the log likelihood is a sum over the cases, and so are the gaps
# between fits that lavaan's optimizer, which works per case, accepts.
#
# When the first fit reaches the maximum (its log likelihood within `tol`
# of the check's) and no factor is at the boundary, it is kept. Otherwise
# the maximum is carried back to fa_second_order()'s parameterization
# (gamma_f = 1 / sqrt(z_f), lambda_i = a_i sqrt(z_f)), turned to the
# reported direction (.fa_second_order_orient_point()). With no factor at
# the boundary that point is a maximum of the first parameterization too,
# and nlminb, started there, reports the fit there. A factor at the
# boundary has no finite gamma_f, and its point is held at z_f = 1e-8
# (gamma_f = 1e4, a standardized loading of 1 - 5e-9): `fit`, the fit at
# that point, has a log likelihood within 1e-8 per case of the maximum on
# the holzinger_swineford examples and on covariance matrices that imply
# standardized loadings from 1.16 to 2.85, at every N from 300 to 3
# million (it must be within `tol`), and the model's degrees of freedom,
# but its information matrix cannot be inverted along gamma_f; `table_fit`
# is the same point with each such gamma_f held there by an equality
# constraint, whose standard errors are those of the other parameters
# with the loading held at one; and `rows_fit` is the check fit at its
# maximum, with the arguments of the first fit, whose fit measures (the
# robust ones included, which lavaan cannot compute at the held point)
# are those of the model at its maximum. A fit computed at a point
# without the optimizer has no bootstrap standard errors
# (.fa_second_order_refit()).
#
# Returns a list: `fit`, `table_fit`, and `rows_fit` (the same fit unless
# a factor is at the boundary); `converged`, FALSE when no fit can be
# reported, either because neither parameterization converged or because
# (`short` TRUE) the maximum the check found could not be reproduced in
# fa_second_order()'s parameterization within `tol`; `verified`, FALSE
# when the check could not be made for a first fit that converged (no
# check fit converged); `replaced`, TRUE when the fit is not the first
# one, reflection aside; `at_one`, the standardized second-order loadings
# of the factors at the boundary (signed, as reported), named by
# first-order factor; `held`, the gamma_f at which those are held, named
# by label; `owner`, the second-order factor of each first-order factor;
# `oriented`, FALSE when a reflection could not be applied; `stalled`,
# the absolute standardized second-order loadings of the first fit, where
# it stopped if it did not converge; `warnings`, lavaan's warnings to be
# given: those of the fit kept as `fit`, or, for a fit held at the
# boundary, those of `table_fit` and `rows_fit` (the first fit's are
# dropped whenever it is not the one reported); and, for a fit held at the
# boundary, `rows_model` and `rows_args`, the syntax and the other
# arguments rows_fit was fit with, its start aside (NULL otherwise).
.fa_second_order_maximum <- function(fit, fit_args, model_check, factors,
                                     second_order, eq_gamma, warnings) {
  first <- names(factors)
  items <- unlist(factors, use.names = FALSE)
  f_of <- stats::setNames(rep(first, lengths(factors)), items)
  owner <- stats::setNames(rep(names(second_order), lengths(second_order)),
                           unlist(second_order, use.names = FALSE))
  z_flag <- 1e-6
  z_hold <- 1e-8
  # Two fits are compared per case (see above).
  tol <- 1e-6 * lavaan::lavInspect(fit, "ntotal")
  logl <- function(x) as.numeric(lavaan::logLik(x))
  std_of <- function(gamma) gamma / sqrt(1 + gamma^2)
  out <- list(fit = fit, table_fit = fit, rows_fit = fit, converged = TRUE,
              verified = TRUE, replaced = FALSE,
              at_one = stats::setNames(numeric(0), character(0)),
              held = stats::setNames(numeric(0), character(0)),
              oriented = TRUE, stalled = NULL, owner = owner,
              short = FALSE, warnings = warnings)
  conv <- isTRUE(lavaan::lavInspect(fit, "converged"))
  pt <- lavaan::parTable(fit)
  p_s <- .fa_second_order_point(pt, factors, second_order)
  usable <- all(is.finite(pt$est[pt$free > 0L]))
  ll_s <- if (conv) logl(fit) else -Inf
  if (usable) out$stalled <- abs(std_of(p_s$second))

  # ---- the check fits --------------------------------------------------
  # Point estimates and the likelihood are all a check needs: no standard
  # errors, test, or saturated and baseline models
  # (.fa_second_order_check_fit()). Beside the bound on z_f, lavaan 0.7.3
  # does not enforce equal labels written as constraints (a check of equal
  # error variances reached a larger log likelihood at a point where they
  # were 31, 14, and 52); with ceq.simple = TRUE each label is one
  # parameter, and the equalities hold exactly.
  check_args <- fit_args
  check_args[c("model", "start", "se", "test", "h1", "baseline",
               "optim.method", "ceq.simple")] <- list(
                 model_check, NULL, "none", "none", FALSE, FALSE, "nlminb",
                 TRUE)
  check_args <- check_args[!vapply(check_args, is.null, logical(1L))]
  run_check <- function(start) {
    args <- check_args
    if (!is.null(start)) {
      args$start <- start
      args$check.start <- FALSE
    }
    .fa_second_order_check_fit(args)
  }
  z_of <- function(chk) {
    .fa_second_order_point(lavaan::parTable(chk), factors,
                           second_order)$second
  }
  checks <- list()
  if (usable) {
    gamma <- p_s$second
    z <- ifelse(abs(gamma) > 1e-4, 1 / gamma^2, 1e8)
    a <- p_s$load * ifelse(abs(gamma) > 1e-4, gamma, 1e-4)[f_of[items]]
    checks <- c(checks, list(run_check(.fa_second_order_start(
      a, z, p_s$rest, factors, second_order, check = TRUE))))
  }
  first_check <- if (length(checks)) checks[[1L]]
  if (is.null(first_check) || !conv || any(z_of(first_check) <= z_flag) ||
      logl(first_check) > ll_s + tol) {
    checks <- c(checks, list(run_check(NULL)))
  }
  # The first fit kept, turned to the reported direction: reflected and fit
  # again from there, which nlminb accepts at once, or computed at the
  # reflected point should it not; when neither can be done the fit is kept
  # as it is and `oriented` says so.
  keep_first <- function(out) {
    o <- .fa_second_order_orient_point(p_s$load, p_s$second, p_s$rest,
                                       factors, second_order, eq_gamma)
    if (!length(o$reflected)) return(out)
    start <- .fa_second_order_start(o$lambda, o$gamma, o$rest, factors,
                                    second_order)
    r <- .fa_second_order_refit(fit_args, start, "polish")
    if (is.null(r$fit) || logl(r$fit) < ll_s - tol) {
      r <- .fa_second_order_refit(fit_args, start, "none")
    }
    if (is.null(r$fit)) {
      out$oriented <- FALSE
    } else {
      out$fit <- out$table_fit <- out$rows_fit <- r$fit
    }
    out
  }
  checks <- checks[!vapply(checks, is.null, logical(1L))]
  if (!length(checks)) {
    out$converged <- conv
    out$verified <- FALSE
    return(if (conv) keep_first(out) else out)
  }
  best <- checks[[which.max(vapply(checks, logl, numeric(1L)))]]
  ll_b <- logl(best)
  p_b <- .fa_second_order_point(lavaan::parTable(best), factors,
                                second_order)
  at <- names(p_b$second)[p_b$second <= z_flag]

  # ---- the first fit is the maximum ------------------------------------
  if (conv && !length(at) && ll_b <= ll_s + tol) return(keep_first(out))

  # ---- the maximum is elsewhere ----------------------------------------
  out$replaced <- TRUE
  z <- p_b$second
  z[at] <- z_hold
  gamma <- 1 / sqrt(z)
  o <- .fa_second_order_orient_point(p_b$load * sqrt(z)[f_of[items]],
                                     gamma, p_b$rest, factors,
                                     second_order, eq_gamma)
  start <- .fa_second_order_start(o$lambda, o$gamma, o$rest, factors,
                                  second_order)
  if (!length(at)) {
    r <- .fa_second_order_refit(fit_args, start, "polish")
    if (is.null(r$fit) || logl(r$fit) < ll_b - tol) {
      r <- .fa_second_order_refit(fit_args, start, "none")
    }
    if (is.null(r$fit) || logl(r$fit) < ll_b - tol) {
      out$converged <- FALSE
      out$short <- TRUE
      return(out)
    }
    out$fit <- out$table_fit <- out$rows_fit <- r$fit
    out$warnings <- r$warnings
    return(out)
  }
  # A factor at the boundary: the fit at the held point, the same point
  # with the held loadings constrained there, and the check fit with the
  # arguments of the first fit.
  # Each held loading by its label (.fa_second_order_labels()), once for a
  # label that equated loadings share.
  held <- stats::setNames(o$gamma[at], ifelse(
    eq_gamma[owner[at]], paste0("gamma_", owner[at]),
    paste0("gamma_", owner[at], "_", at)))
  held <- held[!duplicated(names(held))]
  free_fit <- .fa_second_order_refit(fit_args, start, "none")$fit
  table_r <- .fa_second_order_refit(
    fit_args, start, "none",
    model = paste0(fit_args$model, "\n",
                   paste0(names(held), " == ", sprintf("%.17g", held),
                          collapse = "\n")))
  # The check fit again, with the arguments of the first fit, from its
  # maximum: nlminb stops there at once, with the bound marked active, which
  # lavaan's robust test reads (computed at the point without the optimizer,
  # the scaling factor of an MLR fit on the holzinger_swineford example was
  # .790 against .985).
  check_start <- .fa_second_order_start(p_b$load, p_b$second, p_b$rest,
                                        factors, second_order, check = TRUE)
  # Only its fit rows are read, so it is not bootstrapped.
  rows_args <- fit_args
  rows_args$ceq.simple <- TRUE
  if (identical(rows_args$se, "bootstrap")) rows_args$se <- "none"
  rows_r <- .fa_second_order_refit(rows_args, check_start, "polish",
                                   model = model_check)
  if (is.null(rows_r$fit) || logl(rows_r$fit) < ll_b - tol) {
    rows_r <- .fa_second_order_refit(rows_args, check_start, "none",
                                     model = model_check)
  }
  if (is.null(free_fit) || is.null(table_r$fit) || is.null(rows_r$fit) ||
      logl(free_fit) < ll_b - tol) {
    out$converged <- FALSE
    out$short <- TRUE
    return(out)
  }
  out$fit <- free_fit
  out$table_fit <- table_r$fit
  out$rows_fit <- rows_r$fit
  # The syntax and arguments of rows_fit, without the start it was fit
  # from: the lavaan code of a refusal of rows_fit (.no_max_refuse() in
  # .fa_cfa_engine()).
  out$rows_model <- model_check
  out$rows_args <- rows_args
  out$rows_args$model <- model_check
  out$at_one <- stats::setNames(sign(o$gamma[at]) / sqrt(1 + p_b$second[at]),
                                at)
  out$held <- held
  # lavaan's warnings about the fits the table is read from, which repeat
  # any it gave about the data (cases that are empty, for one), each once.
  # The fit at the held point warns that the covariance matrix of its
  # estimates is not positive definite, which is the boundary itself, and
  # the caller's own warning describes it.
  w <- c(table_r$warnings, rows_r$warnings)
  out$warnings <- w[!duplicated(vapply(w, conditionMessage, character(1L)))]
  out
}

# .fa_second_order_maximal(fit, factors, second_order, items)
#
# The maximal reliability of each second-order factor and of the
# disturbance of each first-order factor (the group factor of the
# Schmid-Leiman solution, what the first-order factor measures beyond its
# second-order factor): c' Sigma^{-1} c / Var, with c the covariances of
# the items with the factor and Sigma the model implied covariance matrix
# of the items. With L the first-order loadings, G the second-order
# loadings, P the covariance matrix of the second-order factors, and Z the
# disturbance variances, Sigma = L (G P G' + Z) L' + Psi, the covariances
# with the second-order factors are L G P (an item of another second-order
# factor covaries with g through their correlation), and those with the
# disturbances L Z. The standard error is .fa_layered_maximal()'s
# (`vcov_usable` as there). Returns its list, named by the second-order
# factors and then the first-order factors.
.fa_second_order_maximal <- function(fit, factors, second_order, items,
                                     vcov_usable = TRUE) {
  first <- names(factors)
  second <- names(second_order)
  pt <- lavaan::parTable(fit)
  load <- which(pt$op == "=~" & pt$lhs %in% first & pt$rhs %in% items)
  gam <- which(pt$op == "=~" & pt$lhs %in% second & pt$rhs %in% first)
  errs <- which(pt$op == "~~" & pt$lhs == pt$rhs & pt$lhs %in% items)
  dist <- which(pt$op == "~~" & pt$lhs == pt$rhs & pt$lhs %in% first)
  phis <- which(pt$op == "~~" & pt$lhs %in% second & pt$rhs %in% second)
  cells_l <- cbind(match(pt$rhs[load], items), match(pt$lhs[load], first))
  cells_g <- cbind(match(pt$rhs[gam], first), match(pt$lhs[gam], second))
  cells_p <- cbind(match(pt$lhs[phis], second), match(pt$rhs[phis], second))
  ei <- match(pt$lhs[errs], items)
  di <- match(pt$lhs[dist], first)
  latent <- c(second, first)
  .fa_layered_maximal(fit, latent, function(e) {
    L <- matrix(0, length(items), length(first))
    L[cells_l] <- e[load]
    G <- matrix(0, length(first), length(second))
    G[cells_g] <- e[gam]
    P <- matrix(0, length(second), length(second))
    P[cells_p] <- e[phis]
    P[cells_p[, 2:1, drop = FALSE]] <- e[phis]
    zeta <- numeric(length(first))
    zeta[di] <- e[dist]
    psi <- numeric(length(items))
    psi[ei] <- e[errs]
    Z <- diag(zeta, length(first))
    Sigma <- L %*% (G %*% P %*% t(G) + Z) %*% t(L)
    diag(Sigma) <- diag(Sigma) + psi
    C <- cbind(L %*% G %*% P, L %*% Z)
    SiC <- tryCatch(solve(Sigma, C), error = function(e) NULL)
    if (is.null(SiC)) return(rep(NA_real_, length(latent)))
    colSums(C * SiC) / c(diag(P), zeta)
  }, vcov_usable)
}

# .fa_second_order_schmid_leiman(pe, factors, second_order, items,
#                                at_one = character(0))
#
# The Schmid-Leiman solution a second-order fit implies (`pe`, lavaan's
# parameterEstimates()): for each item, its first-order factor and that
# factor's second-order factor, its standardized loading on the
# second-order factor, lambda_std gamma_std, its standardized loading on
# the first-order factor's disturbance (the group factor of the
# solution), lambda_std sqrt(1 - gamma_std^2), and its uniqueness, one
# minus their squares, where lambda_std is the item's standardized
# first-order loading and gamma_std the factor's standardized second-order
# loading (Schmid & Leiman, 1957; Yung, Thissen, & McLeod, 1999, Eqs. 4
# and 10). The first-order factors of `at_one`, whose standardized
# second-order loading is one at the maximum and is held just short of it
# in the fit (.fa_second_order_maximum()), take their limit there: a group
# loading of 0. Returns a list: `table`, that data frame; `L`, the items
# by factors matrix of the unstandardized Schmid-Leiman loadings, columns
# the second-order factors and then the first-order factors; and `psi`,
# the error variances named by item.
.fa_second_order_schmid_leiman <- function(pe, factors, second_order,
                                           items, at_one = character(0)) {
  first <- names(factors)
  second <- names(second_order)
  ld <- pe[pe$op == "=~" & pe$lhs %in% first & pe$rhs %in% items, ,
           drop = FALSE]
  gm <- pe[pe$op == "=~" & pe$lhs %in% second & pe$rhs %in% first, ,
           drop = FALSE]
  var_of <- function(v) {
    r <- pe$op == "~~" & pe$lhs == pe$rhs & pe$lhs %in% v
    stats::setNames(pe$est[r], pe$lhs[r])[v]
  }
  psi <- var_of(items)
  zeta <- var_of(first)
  phi <- var_of(second)
  f_of <- stats::setNames(ld$lhs, ld$rhs)[items]
  g_of <- stats::setNames(gm$lhs, gm$rhs)[f_of]
  lambda <- stats::setNames(ld$est, ld$rhs)[items]
  gamma <- stats::setNames(gm$est, gm$rhs)[f_of]
  # The covariance of an item with its second-order factor over that
  # factor's standard deviation, and with the disturbance over the
  # disturbance's, the unstandardized Schmid-Leiman loadings; both
  # variances are 1 as the model is identified.
  general <- unname(lambda * gamma * sqrt(phi[g_of]))
  group <- unname(lambda * sqrt(zeta[f_of]))
  group[unname(f_of) %in% at_one] <- 0
  v <- general^2 + group^2 + unname(psi)
  sd_v <- ifelse(v > 0, sqrt(pmax(v, 0)), NA_real_)
  table <- data.frame(item = items, factor = unname(f_of),
                      second_order = unname(g_of),
                      general_loading = general / sd_v,
                      group_loading = group / sd_v,
                      uniqueness = unname(psi) / v,
                      stringsAsFactors = FALSE)
  L <- matrix(0, length(items), length(second) + length(first),
              dimnames = list(items, c(second, first)))
  L[cbind(seq_along(items), match(g_of, second))] <- general
  L[cbind(seq_along(items), length(second) + match(f_of, first))] <- group
  list(table = table, L = L, psi = psi)
}

# .fa_second_order_improper(pe, second_order, at_one)
#
# The second-order part of an improper solution: the standardized
# second-order loadings at one, `at_one` (named by first-order factor, from
# .fa_second_order_maximum(), which finds them at the maximum of the
# likelihood in a parameterization that reaches the boundary), and the
# correlations of the second-order factors beyond one in absolute value,
# from lavaan's parameterEstimates() `pe`. Returns a list: `loading`, the
# standardized loadings at one, named gamma_<g>_<f>; `factor`, the
# first-order factors they belong to; `correlation`, the correlations
# beyond one, named phi_<g1>_<g2>; and `pair`, their second-order factors
# as "<g1> and <g2>".
.fa_second_order_improper <- function(pe, second_order, at_one) {
  second <- names(second_order)
  owner <- stats::setNames(rep(second, lengths(second_order)),
                           unlist(second_order, use.names = FALSE))
  var_of <- function(v) {
    r <- pe$op == "~~" & pe$lhs == pe$rhs & pe$lhs %in% v
    stats::setNames(pe$est[r], pe$lhs[r])[v]
  }
  cr <- pe[pe$op == "~~" & pe$lhs != pe$rhs & pe$lhs %in% second &
             pe$rhs %in% second, , drop = FALSE]
  r <- cr$est / sqrt(unname(var_of(cr$lhs)) * unname(var_of(cr$rhs)))
  beyond <- !is.finite(r) | abs(r) > 1
  # (sprintf(), unlike paste0(), gives no name when there is none.)
  list(loading = stats::setNames(unname(at_one),
                                 sprintf("gamma_%s_%s",
                                         owner[names(at_one)],
                                         names(at_one))),
       factor = names(at_one),
       correlation = stats::setNames(r[beyond],
                                     sprintf("phi_%s_%s", cr$lhs[beyond],
                                             cr$rhs[beyond])),
       pair = sprintf("%s and %s", cr$lhs[beyond], cr$rhs[beyond]))
}

# .fa_second_order_rows(params, fit, pe, factors, second_order, items,
#                       conf_level, wald_df, vcov_usable, at_one)
#
# The parameter table of a second-order fit with the rows lavaan cannot
# define inserted and the rows an improper solution leaves undefined set to
# NA (.fa_layered_rows()): the maximal reliability of each second-order
# factor and of each first-order factor's disturbance
# (.fa_second_order_maximal(), its standard error, test, and limits NA
# when `vcov_usable` is FALSE), and PUC_<g> for each second-order factor,
# the proportion of the pairs of its items that do not share a
# first-order factor. An improper solution at the second order
# (.fa_second_order_improper()) is handled as a negative error variance
# is. A standardized second-order loading at one at the maximum of the
# likelihood (`at_one`, from .fa_second_order_maximum()) leaves the
# factor's disturbance without variance, which this parameterization
# reaches only as the second-order loading gamma_f runs off to infinity
# and the factor's loadings lambda_i on the disturbance run to zero: the
# fit holds the point just short of it, so gamma_f and the lambda_i have
# no finite estimate and are NA, and so are the rows that describe the
# disturbance (omega_HS, ECV, H, and the maximal reliability of that
# first-order factor). A correlation of the second-order factors beyond
# one leaves the model implied covariance matrix of the factors outside
# what a covariance matrix can be, so every maximal reliability, which
# reads it, is NA. A warning names them. Returns .fa_layered_rows()'s list
# with the Schmid-Leiman table (.fa_second_order_schmid_leiman()) as
# `schmid_leiman` and the second-order part of an improper solution as
# `improper_second_order` (the standardized loadings at one and the
# correlations beyond one, named by their terms; empty for a proper
# solution).
.fa_second_order_rows <- function(params, fit, pe, factors, second_order,
                                  items, conf_level, wald_df,
                                  vcov_usable = TRUE,
                                  at_one = stats::setNames(numeric(0),
                                                           character(0))) {
  latent <- unlist(lapply(names(second_order), function(g) {
    c(g, second_order[[g]])
  }), use.names = FALSE)
  owner <- stats::setNames(rep(names(second_order), lengths(second_order)),
                           unlist(second_order, use.names = FALSE))
  sl <- .fa_second_order_schmid_leiman(pe, factors, second_order, items,
                                       names(at_one))
  pairs <- function(k) k * (k - 1) / 2
  puc <- vapply(second_order, function(fs) {
    J <- lengths(factors[fs])
    1 - sum(pairs(J)) / pairs(sum(J))
  }, numeric(1L))
  out <- .fa_layered_rows(
    params, latent = latent,
    maximal = .fa_second_order_maximal(fit, factors, second_order, items,
                                       vcov_usable),
    puc = puc, L = sl$L[, latent, drop = FALSE], psi = sl$psi,
    generals = lapply(second_order, function(fs) {
      unlist(factors[fs], use.names = FALSE)
    }),
    groups = factors, conf_level = conf_level, wald_df = wald_df)
  bad <- .fa_second_order_improper(pe, second_order, at_one)
  if (length(bad$loading) || length(bad$correlation)) {
    p <- out$params
    # The rows of a factor at one: its loadings and its second-order
    # loading (by their syntax, which equated labels do not change), and
    # the rows of its disturbance.
    loading_rows <- lapply(stats::setNames(bad$factor, bad$factor),
                           function(f) {
      unique(p$term[p$syntax %in% c(paste(f, "=~", factors[[f]]),
                                    paste(owner[[f]], "=~", f))])
    })
    disturbance <- lapply(stats::setNames(bad$factor, bad$factor),
                          function(f) {
      paste0(c("omega_HS_", "ECV_", "H_", "maximal_reliability_"), f)
    })
    na_terms <- c(unlist(loading_rows, use.names = FALSE),
                  unlist(disturbance, use.names = FALSE),
                  if (length(bad$correlation)) {
                    paste0("maximal_reliability_", latent)
                  })
    na_cols <- c("estimate", "se", "test_statistic", "p_value", "ci_lower",
                 "ci_upper")
    out$params[out$params$term %in% na_terms, na_cols] <- NA_real_
    .dmar_warn(c(
      .fa_second_order_at_one_message(bad$factor, owner,
                                      loading_rows, disturbance,
                                      se = vcov_usable),
      if (length(bad$correlation)) {
        paste0(if (length(bad$factor)) " Also, the" else
                 "Improper solution at the second order: the",
               " correlation of the second-order factors ",
               paste0(bad$pair, " (",
                      formatC(unname(bad$correlation), format = "f",
                              digits = 3), ")", collapse = ", "),
               if (length(bad$correlation) == 1L) " is" else " are",
               " beyond one, so the maximal reliabilities, which read ",
               if (length(bad$correlation) == 1L) "it" else "them",
               ", are NA.")
      }), "result_qualified")
  }
  out$improper_second_order <- stats::setNames(
    c(unname(bad$loading), unname(bad$correlation)),
    c(names(bad$loading), names(bad$correlation)))
  out$schmid_leiman <- sl$table
  out
}

# .fa_second_order_at_one_message(at, owner, loading_rows = NULL,
#                                 disturbance = NULL, held = NULL,
#                                 se = TRUE)
#
# The warning for first-order factors `at` whose standardized second-order
# loading is one at the maximum of the likelihood (`owner`, the
# second-order factor of each first-order factor). For a table
# (`loading_rows` and `disturbance`, the terms of each factor set NA) it
# names the rows left NA, and, when the table has standard errors (`se`),
# says that they hold the loading at one; for the lavaan fit or its
# summary or standardized solution (`held`, the second-order loadings at
# which the fit holds them, named by label) it says where the point is
# held. "" when `at` is empty.
.fa_second_order_at_one_message <- function(at, owner, loading_rows = NULL,
                                            disturbance = NULL,
                                            held = NULL, se = TRUE) {
  n <- length(at)
  if (!n) return("")
  one <- n == 1L
  gammas <- paste0("gamma_", owner[at], "_", at)
  paste0(
    "Improper solution at the second order: at the maximum of the ",
    "likelihood the standardized second-order loading",
    if (!one) "s", " of ", paste(at, collapse = ", "),
    if (one) " is" else " are", " one, so the disturbance of ",
    if (one) "that factor has" else "each of those factors has",
    " no variance (a parameterization that let it be negative would ",
    "estimate it at or below zero), a point fa_second_order()'s ",
    "parameterization reaches only as ", paste(gammas, collapse = ", "),
    if (one) " grows" else " grow", " without bound. ",
    if (!is.null(loading_rows)) {
      paste0(
        "The second-order loading", if (!one) "s", " and the loadings of ",
        "the factor", if (!one) "s", " on ", if (one) "its" else "their",
        " disturbance", if (!one) "s", " (",
        paste(unlist(loading_rows, use.names = FALSE), collapse = ", "),
        ") have no finite estimate and are NA, as are the rows of the ",
        "disturbance", if (!one) "s", " (",
        paste(unlist(disturbance, use.names = FALSE), collapse = ", "),
        "); the fit rows are those of the maximum",
        if (se) {
          paste0(", and the other standard errors hold the loading",
                 if (!one) "s", " at one")
        }, ".")
    } else {
      paste0(
        "The fit holds ", paste0(names(held), " at ",
                                 formatC(unname(held), format = "f", digits = 0),
                                 collapse = ", "),
        ", a standardized loading within 5e-9 of one, where the log ",
        "likelihood is within 1e-6 per case of its maximum, so ",
        if (one) "that value and the factor's" else
          "those values and the factors'",
        " loadings on ", if (one) "its" else "their", " disturbance",
        if (!one) "s", " are not estimates.")
    })
}

# .fa_second_order_near_one(x)
#
# Standardized loadings near one as a message prints them: to four
# decimals, and "0.9999 or more" for those that four decimals would round
# to one, which they are not.
.fa_second_order_near_one <- function(x) {
  ifelse(abs(x) >= 0.99995, "0.9999 or more",
         formatC(unname(x), format = "f", digits = 4))
}

# .fa_second_order_not_converged(so_max)
#
# The message of the error fa_second_order() gives when no fit of the
# model can be reported (`so_max`, .fa_second_order_maximum()'s list):
# neither its own parameterization nor the one that reaches a standardized
# second-order loading of one converged, and, when the first fit stopped
# with a standardized second-order loading of .99 or more, which way it was
# heading; or the maximum the second found could not be reproduced in the
# first.
.fa_second_order_not_converged <- function(so_max) {
  if (isTRUE(so_max$short)) {
    return(paste0(
      "The model did not converge: fa_second_order() found the maximum of ",
      "the likelihood in a second parameterization, one that reaches a ",
      "standardized second-order loading of one, but could not reproduce ",
      "it in its own, so no estimate is reported."))
  }
  st <- so_max$stalled
  near <- st[is.finite(st) & abs(st) >= 0.99]
  paste0(
    "The model did not converge, in fa_second_order()'s parameterization ",
    "or in a second one that reaches a standardized second-order loading ",
    "of one.",
    if (length(near)) {
      paste0(
        " Where the first fit stopped, the standardized second-order ",
        "loading", if (length(near) > 1L) "s", " of ",
        paste0(names(near), " (", .fa_second_order_near_one(near), ")",
               collapse = ", "),
        if (length(near) > 1L) " were" else " was", " near one, which ",
        "suggests data that imply a loading at or beyond one (a ",
        "disturbance variance at or below zero); fa_confirmatory() ",
        "estimates the correlations of the first-order factors, which show ",
        "whether a second-order factor can account for them.")
    })
}

# .fa_second_order_summary(fit, at_one)
#
# lavaan's summary (with fit measures) of a second-order fit computed at a
# point without the optimizer: one that holds a standardized second-order
# loading at one (`at_one` TRUE) or a fallback of .fa_second_order_refit().
# Under FIML lavaan forms its robust CFI, TLI, RMSEA, and GFI from a
# corrected fit (lav_fit_fiml_corrected()) that it cannot compute for such a
# fit; it reports them as NA with a warning from each of
# lav_fit_cfi_lavobject(), lav_fit_rmsea_lavobject(), and
# lav_fit_gfi_lavobject(). Those warnings, matched by the name of the
# function in them (the message itself is translated), are replaced by
# one. Any other warning passes through.
.fa_second_order_summary <- function(fit, at_one) {
  seen <- new.env(parent = emptyenv())
  seen$failed <- FALSE
  out <- withCallingHandlers(
    lavaan::summary(fit, fit.measures = TRUE),
    warning = function(w) {
      if (grepl("lav_fit_(cfi|rmsea|gfi)_lavobject", conditionMessage(w))) {
        seen$failed <- TRUE
        invokeRestart("muffleWarning")
      }
    })
  if (seen$failed) {
    .dmar_warn(c(
      "lavaan could not compute its robust CFI, TLI, RMSEA, and GFI, ",
      "which under FIML it forms from a corrected fit, for this fit, ",
      "computed without the optimizer at ",
      if (at_one) {
        paste0("the point that holds the standardized second-order ",
               "loading at one")
      } else {
        "the point found"
      },
      ", so the summary shows them as NA."), "result_qualified")
  }
  out
}

# .fa_second_order_check_notes(so_max, output, se)
#
# The warnings fa_second_order() gives about its search for the maximum
# (`so_max`, .fa_second_order_maximum()'s list) before anything is built
# from the fit: that the search could not be made, naming any standardized
# second-order loading of .995 or more (which three significant digits
# show as 1.00, and which without the search cannot be told from a point
# where the optimizer stopped on its way to one); that a factor could not
# be turned to the reported direction; for the outputs that are not
# a DMAR table ("fit", "summary", "standardized"), which factors are at
# one and where the fit holds them (a table names them, and the rows it
# leaves NA, through .fa_second_order_rows()); and, when bootstrap
# standard errors were asked for (`se`) and the fit reported is one
# computed at a point without the optimizer (a loading held at one, or a
# fallback of .fa_second_order_refit()), that they are not available (NA
# in a table; absent from the lavaan object of "fit", "summary", and
# "standardized", whose se is "none").
.fa_second_order_check_notes <- function(so_max, output, se = NULL) {
  if (!so_max$verified) {
    near <- so_max$stalled
    near <- near[is.finite(near) & abs(near) >= 0.995]
    .dmar_warn(c(
      "fa_second_order() could not confirm that this fit is the ",
      "maximum of the likelihood: its check, a second fit in a ",
      "parameterization that reaches a standardized second-order ",
      "loading of one, did not converge.",
      if (length(near)) {
        paste0(" The standardized second-order loading",
               if (length(near) > 1L) "s", " of ",
               paste0(names(near), " (",
                      .fa_second_order_near_one(near), ")",
                      collapse = ", "),
               " may be where the optimizer stopped on its way to ",
               "one rather than ",
               if (length(near) > 1L) "estimates" else "an estimate",
               ".")
      }), "result_qualified")
  }
  if (!so_max$oriented) {
    .dmar_warn(c(
      "A first-order or second-order factor of this fit could not be ",
      "turned to the reported direction (loadings that sum to a ",
      "positive value), so it is reported as the optimizer left it, ",
      "which changes the signs of its loadings and nothing else."),
      "result_qualified")
  }
  if (length(so_max$at_one) &&
      output %in% c("fit", "summary", "standardized")) {
    .dmar_warn(
      .fa_second_order_at_one_message(names(so_max$at_one),
                                      so_max$owner,
                                      held = so_max$held), "result_qualified")
  }
  if (identical(se, "bootstrap") &&
      identical(lavaan::lavInspect(so_max$table_fit, "options")$se, "none")) {
    .dmar_warn(c(
      "Bootstrap standard errors and intervals are not available for ",
      "this fit: ",
      if (length(so_max$at_one)) {
        paste0("with a standardized second-order loading at one, the ",
               "fit is computed at the point that holds the loading ",
               "there")
      } else {
        "the fit is computed at the point found"
      },
      " without the optimizer, and a bootstrap of it would not refit ",
      "the samples, so ",
      if (output %in% c("fit", "summary", "standardized")) {
        paste0("the ", output, " output carries no standard errors, ",
               "tests, or intervals (lavaan's se = \"none\").")
      } else {
        "every standard error, test, and interval is NA."
      }), "result_qualified")
  }
  invisible(NULL)
}

# .fa_second_order_attributes(out, factors, second_order, eq_load, eq_err,
#                             eq_int, eq_gamma, has_means, improper,
#                             improper_second_order, schmid_leiman)
#
# Records the second-order structure on a fa_second_order() table: the
# first-order factors ("factors") and the second-order factors
# ("second_order"), the factors whose loadings, error variances,
# intercepts, and second-order loadings are equated ("equal", a list),
# whether the mean structure is modeled ("has_means"), the degrees of
# freedom of the structural part ("structural_df",
# .fa_second_order_structural_df()), the metric of the omega rows
# ("omega_metric", "raw_item_sum"), the negative error variances of an
# improper solution ("improper", empty for a proper one), its standardized
# second-order loadings at one and second-order correlations beyond one
# ("improper_second_order", empty for a proper one), and the
# Schmid-Leiman solution ("schmid_leiman"). The header printed above the
# table is built from them (.dmar_tbl_footer.dmar_fa_second_order(),
# R/fa_second_order.R).
.fa_second_order_attributes <- function(out, factors, second_order, eq_load,
                                        eq_err, eq_int, eq_gamma, has_means,
                                        improper, improper_second_order,
                                        schmid_leiman) {
  attr(out, "factors") <- factors
  attr(out, "second_order") <- second_order
  attr(out, "equal") <- list(loading = names(eq_load)[eq_load],
                             error = names(eq_err)[eq_err],
                             intercept = if (has_means) {
                               names(eq_int)[eq_int]
                             } else {
                               character(0)
                             },
                             second_order_loading = names(eq_gamma)[eq_gamma])
  attr(out, "has_means") <- has_means
  attr(out, "structural_df") <- .fa_second_order_structural_df(second_order,
                                                               eq_gamma)
  attr(out, "omega_metric") <- "raw_item_sum"
  attr(out, "improper") <- improper
  attr(out, "improper_second_order") <- improper_second_order
  attr(out, "schmid_leiman") <- schmid_leiman
  out
}

# .fa_cfa_tidy(x), .fa_cfa_glance(x)
#
# The tidy() and glance() of a fa_confirmatory() or fa_bifactor() table,
# whose columns and fit rows are the same: tidy() gives the parameter and
# index rows, glance() one row of the fit rows and what produced them.
.fa_cfa_tidy <- function(x) {
  fit_terms <- c("chi_square", "df", "p_chi_square", "cfi", "tli", "nnfi",
                 "rmsea", "rmsea_ci_lower", "rmsea_ci_upper",
                 "rmsea_ci_level", "srmr", "AIC", "BIC", "H0", "H1")
  rows <- x[!x$term %in% fit_terms, , drop = FALSE]
  data.frame(
    term      = rows$term,
    estimate  = rows$estimate,
    se = rows$se,
    statistic = rows$test_statistic,
    p_value   = rows$p_value,
    ci_lower  = rows$ci_lower,
    ci_upper = rows$ci_upper,
    # The level of the parameter intervals, so that a tidy row states the
    # level of its limits (the RMSEA level is in rmsea_ci_level, a
    # fit row, which tidy() leaves to glance()).
    conf_level = rep(attr(x, "conf_level", exact = TRUE) %||% NA_real_,
                     nrow(rows)),
    # The reference of each row's interval and p-value, which a stack of
    # tables would otherwise lose (the wald_df attribute is kept from the
    # first frame only), and which differs between rows of one table under
    # the keyword (.wald_df_column(), R/missing_information_internals.R).
    wald_df = .wald_df_column(attr(x, "wald_df", exact = TRUE), nrow(rows),
                              rows$term),
    h1_diagnostic = rep(attr(x, "h1_diagnostic", exact = TRUE) %||% NA_character_, nrow(rows)),
    wald_calibration_note = rep(attr(x, "wald_calibration_note", exact = TRUE) %||% NA_character_, nrow(rows)),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
}

.fa_cfa_glance <- function(x) {
  pick <- function(term) {
    val <- x$estimate[x$term == term]
    if (length(val) == 0L) NA_real_ else val[1L]
  }
  data.frame(
    chi_square = pick("chi_square"),
    df         = pick("df"),
    p_value    = pick("p_chi_square"),
    cfi        = pick("cfi"),
    tli        = pick("tli"),
    rmsea      = pick("rmsea"),
    rmsea_low  = pick("rmsea_ci_lower"),
    rmsea_high = pick("rmsea_ci_upper"),
    srmr       = pick("srmr"),
    AIC        = pick("AIC"),
    BIC        = pick("BIC"),
    logLik     = pick("H0"),
    estimator  = attr(x, "estimator", exact = TRUE) %||% NA_character_,
    missing    = attr(x, "missing", exact = TRUE) %||% NA_character_,
    test       = attr(x, "test", exact = TRUE) %||% NA_character_,
    fit_indices = attr(x, "fit_indices", exact = TRUE) %||% NA_character_,
    rmsea_ci_level = pick("rmsea_ci_level"),
    h1_diagnostic = attr(x, "h1_diagnostic", exact = TRUE) %||% NA_character_,
    wald_calibration_note = attr(x, "wald_calibration_note", exact = TRUE) %||% NA_character_,
    stringsAsFactors = FALSE
  )
}

# .cfa_k_unlinked_factors(data, factors)
#
# The pairs of factors no case of `data` observes together through any of
# their items (every item of one never observed with any item of the
# other), as "'A' with an item of 'B'" for the refusal in fa_confirmatory(). Such a
# pair's correlation does not enter the likelihood of a model with simple
# structure, under FIML or two-stage estimation. Returns character(0) when
# every pair of factors has at least one pair of items observed together.
.cfa_k_unlinked_factors <- function(data, factors) {
  k <- length(factors)
  if (k < 2L) return(character(0))
  observed <- !is.na(as.matrix(data[, unlist(factors, use.names = FALSE),
                                    drop = FALSE]))
  together <- crossprod(observed)
  out <- character(0)
  for (i in seq_len(k - 1L)) {
    for (j in (i + 1L):k) {
      if (all(together[factors[[i]], factors[[j]]] == 0)) {
        out <- c(out, paste0("'", names(factors)[i], "' with an item of '",
                             names(factors)[j], "'"))
      }
    }
  }
  out
}

# .cfa_k_logit_limits(table, terms, conf_level, df, est, se, lower, upper)
#
# The confidence limits of the rows of `table` whose `term` is in `terms`
# (the omega, AVE, and H rows of fa_confirmatory()) replaced by Wald
# limits formed on the logit scale from the row's estimate and standard
# error and transformed back (.logistic_ci()), with the critical value of
# the t distribution on `df` degrees of freedom (one per row; Inf gives
# the normal quantile). A row whose standard error is NA, or whose
# estimate is not inside (0, 1), gets NA limits. Returns `table`.
.cfa_k_logit_limits <- function(table, terms, conf_level, df,
                                  est = "estimate", se = "se",
                                  lower = "ci_lower", upper = "ci_upper") {
  rows <- which(table$term %in% terms)
  df <- rep_len(as.numeric(df), nrow(table))
  for (i in rows) {
    lim <- .logistic_ci(table[[est]][i], table[[se]][i],
                        .wald_crit(conf_level, df[i]))
    table[[lower]][i] <- lim$lower
    table[[upper]][i] <- lim$upper
  }
  table
}

# .cfa_k_vcov_status(fit)
#
# Whether lavaan's covariance matrix of the estimates of `fit` can be used:
# "none" when no standard errors were asked for, "failed" when lavaan
# could not compute it (it warns that the information matrix could not be
# inverted, which suggests a parameter not identified at the estimates),
# "singular" when it computed it from an information matrix that is
# singular (.cfa_k_information_singular()), and "usable" otherwise. lavaan
# inverts the observed information (FIML, MLR) without complaint even when
# it is singular, and then reports finite standard errors for parameters
# that are not identified at the estimates (an H standard error of 534 on
# the fixture of the test of this guard). The expected information of a
# listwise ML fit, or of a fit to a covariance matrix, is often refused by
# lavaan itself in that case ("failed"), but not always, so it is tested
# too: on a bifactor fixture of two groups of three items with equal
# error variances in one, lavaan inverted the singular expected
# information of the listwise and covariance matrix fits (the smallest
# eigenvalue of the result was about -5e-18) and gave every index a
# finite interval, and so it did for fa_confirmatory()'s listwise fit
# of the fixture of the test of this guard with equal loadings on its
# other factor (H .853, standard error .028). Through 1.1.0 development
# fa_confirmatory() tested the observed information only, and the
# layered models and reliability_H() both kinds. Every caller uses this
# one decision, so a fit draws the same diagnostic from each.
.cfa_k_vcov_status <- function(fit) {
  opts <- lavaan::lavInspect(fit, "options")
  if (identical(opts$se, "none")) return("none")
  vc <- tryCatch(lavaan::lavInspect(fit, "vcov"), error = function(e) NULL)
  if (is.null(vc) || !length(vc)) return("failed")
  if (.cfa_k_information_singular(fit)) "singular" else "usable"
}

# .cfa_k_no_vcov(fit, table, cols, status)
#
# When lavaan could not compute the covariance matrix of the estimates, it
# leaves the standard errors of the free parameters NA and gives every
# defined parameter (omega, AVE, H) a standard error of 0, so its interval
# has no width; when it computed it from a singular information matrix,
# every standard error is finite but meaningless. Through 1.0.1
# fa_confirmatory() reported both. Here every standard error, test, and
# interval column `cols` of `table` (a parameterEstimates() or
# standardizedSolution() frame) is set to NA, with a warning, when
# standard errors were asked for and `status` (.cfa_k_vcov_status()) says
# they cannot be used. Returns `table`.
.cfa_k_no_vcov <- function(fit, table, cols,
                           status = .cfa_k_vcov_status(fit)) {
  if (status %in% c("none", "usable")) return(table)
  table[cols] <- NA_real_
  .cfa_k_vcov_warn(status)
  table
}

# .cfa_k_vcov_warn(status, failed_rest)
#
# The warning, of class dmar_result_qualified, that the standard errors
# of a fit are withheld because its covariance matrix of the estimates
# cannot be used (`status`, "singular" or "failed", from
# .cfa_k_vcov_status()). The tables of the factor analysis family give it
# (.cfa_k_no_vcov()), and so does reliability_H() before its delta method
# (.H_from_fit()), so the same fit draws the same diagnostic from both.
# `failed_rest` ends the sentence for a fit lavaan could not invert; the
# default speaks of the tables, whose defined rows lavaan gives a standard
# error of 0.
.cfa_k_vcov_warn <- function(status, failed_rest = c(
                               "every standard error, test, and interval ",
                               "is NA rather than the 0 lavaan reports for ",
                               "the defined quantities (omega, AVE, H).")) {
  if (identical(status, "singular")) {
    .dmar_warn(c(
      "The information matrix of this fit is singular (scaled to a ",
      "unit diagonal, its smallest eigenvalue is below 1e-6), which ",
      "suggests a parameter that is not identified at the estimates, ",
      "so every standard error, test, and interval is NA rather than ",
      "the finite values lavaan computes from it, which do not estimate ",
      "the sampling variability of the estimates."),
      "result_qualified")
  } else {
    .dmar_warn(c(
      "lavaan could not compute the standard errors of this fit (its ",
      "information matrix could not be inverted, which suggests a ",
      "parameter that is not identified at the estimates), so ",
      failed_rest), "result_qualified")
  }
  invisible(NULL)
}

# TRUE when the information matrix lavaan used, restricted to the
# directions its equality constraints leave free and scaled to a unit
# diagonal, has a smallest eigenvalue below `tol`. The restriction matters
# for a model whose constraints identify it (two items of one factor with
# equal loadings): lavaan's information is for the unconstrained
# parameters and is singular there, while the constrained model is
# identified. The scaling removes the units of the parameters, so the bound
# does not depend on the metric of the items. It is applied before the
# restriction as well as after it: the basis of the free directions mixes
# parameters, and when their scales differ a direction can mix them into
# one of no information. A second-order loading held at the boundary
# (.fa_second_order_maximum()) has first-order loadings near 1e-4 beside
# others near 1, and with equal loadings the matrix restricted first had a
# smallest scaled eigenvalue of 1.4e-8 where lavaan's standard errors were
# those of the model with the items on the second-order factor (0.093
# restricted after scaling). The diagonal is scaled in absolute value,
# since at a maximum under constraints the information of a single
# unconstrained parameter can be negative. Across ordinary fits (one to
# three factors, N from 40 to 301, complete and incomplete data) the
# smallest scaled eigenvalue was 0.038 or more; for the unidentified
# fixture it was 6.5e-10.
#
# Under ceq.simple = TRUE (a fit passed to reliability_H(), or an option
# passed through ... of fa_confirmatory()) the parameters that share a
# label share one free number and carry no equality constraint, while
# the information is still laid out over the parameter table rows with
# free > 0; the directions the free numbers leave are then the columns of
# K, which places each free number in every row that carries it. Without
# that restriction a model that its labels identify was taken for one that
# is not (smallest eigenvalue about 1e-17): two groups with loadings equal
# by label and the second group's factor variance free, a group.equal
# fit, and, through 1.1.0 development, fa_confirmatory() with
# ceq.simple = TRUE and a two-item factor with equal loadings, whose
# every standard error was then withheld.
.cfa_k_information_singular <- function(fit, tol = 1e-6) {
  info <- tryCatch(lavaan::lavInspect(fit, "information"),
                   error = function(e) NULL)
  if (is.null(info) || !length(info)) return(FALSE)
  info <- (unclass(info) + t(unclass(info))) / 2
  jac <- tryCatch(lavaan::lavInspect(fit, "con.jac"), error = function(e) NULL)
  ceq <- attr(jac, "ceq.idx")
  free <- tryCatch(lavaan::parTable(fit)$free, error = function(e) NULL)
  row_free <- free[free > 0L]
  shared <- length(row_free) == nrow(info) && anyDuplicated(row_free) > 0L
  if (length(ceq) || shared) {
    d <- abs(diag(info))
    s <- ifelse(is.finite(d) & d > 0, 1 / sqrt(d), 1)
    info <- info * outer(s, s)
    J <- if (length(ceq)) unclass(jac)[ceq, , drop = FALSE]
    if (shared) {
      K <- outer(row_free, sort(unique(row_free)), "==") + 0
      basis <- K / s
      if (!is.null(J)) {
        JK <- if (ncol(J) == nrow(K)) J %*% K else J
        if (ncol(JK) != ncol(K)) return(FALSE)
        N <- MASS::Null(t(JK))
        if (!ncol(N)) return(FALSE)
        basis <- basis %*% N
      }
      basis <- qr.Q(qr(basis))
    } else {
      basis <- MASS::Null(t(J * rep(s, each = nrow(J))))
    }
    if (!ncol(basis)) return(FALSE)
    info <- crossprod(basis, info %*% basis)
    info <- (info + t(info)) / 2
  }
  d <- diag(info)
  if (any(!is.finite(d)) || any(d <= 0)) return(TRUE)
  s <- 1 / sqrt(d)
  ev <- eigen(info * outer(s, s), symmetric = TRUE, only.values = TRUE)$values
  min(ev) < tol
}

# The default of 'missing' is FIML. FIML needs a maximum likelihood
# estimator of continuous items with normal theory or Huber-White standard
# errors and a test lavaan can compute under FIML. When the default meets
# anything else (ordered items, WLS, WLSMV, GLS, se = "robust.sem", or a
# Satorra-Bentler kind of test passed through ...), the items are fitted
# as they stand if none of their values is missing, which is what listwise
# deletion does with nothing to delete, and otherwise the call stops and
# names the remedy rather than deleting cases unasked. An explicit request
# is returned unchanged (the checks that follow refuse an impossible one).
.cfa_k_resolve_missing <- function(missing, default, data, items, ordered,
                                   estimator, se, test) {
  if (!default || is.null(data) || !missing %in% c("ml", "fiml")) {
    return(missing)
  }
  present <- intersect(items, colnames(data))
  incomplete <- length(present) > 0L &&
    anyNA(as.data.frame(data)[present])
  why <- character()
  if (!is.null(ordered) && !isFALSE(ordered)) why <- c(why, "ordered items")
  if (is.character(estimator) && estimator %in% c("WLS", "WLSMV", "GLS")) {
    why <- c(why, sprintf("estimator = \"%s\"", estimator))
  }
  if (identical(se, "robust.sem")) why <- c(why, "se = \"robust.sem\"")
  bad_test <- intersect(as.character(test), c("satorra.bentler",
                                              "mean.var.adjusted",
                                              "scaled.shifted"))
  if (length(bad_test)) {
    why <- c(why, sprintf("test = \"%s\"", bad_test))
  }
  if (!length(why)) return(missing)
  if (!incomplete) return("listwise")
  remedy <- c(
    "missing = \"listwise\" to fit the complete cases",
    if ("ordered items" %in% why) "missing = \"pairwise\" for ordered items",
    if (identical(se, "robust.sem")) {
      "se = \"robust.huber.white\" in place of se = \"robust.sem\", which FIML accepts"
    },
    if (length(bad_test)) {
      "estimator = \"MLR\" in place of the test, whose scaled statistic FIML accepts"
    })
  .dmar_refuse_argument(c(
    "Full information maximum likelihood, the default of 'missing', is ",
    "not available with ", paste(why, collapse = " and "), ", and the ",
    "items have missing values, so no cases are deleted unasked. Use ",
    paste(remedy, collapse = ", or "), "."), arg = "missing")
}

# The lines under a fa_confirmatory() table that name the version of its fit rows
# where the row name alone does not: after FIML the cfi, tli, nnfi, and
# rmsea rows are lavaan's robust versions, the corrected indices of Zhang
# and Savalei (2023), which are not computed from the chi_square row
# (lavaan's version 6, or version 3 when the fit has a scaled statistic,
# as with estimator = "MLR"); with no item value missing (`complete`) and
# the standard statistic those versions equal the standard indices, so no
# line names them. The lines after two-stage estimation
# are .ts_footer()'s. When lavaan could not compute the corrected indices
# (an NA), the rows fall back to its scaled or standard versions
# (`fallback`), and the line says so (.cfa_k_fallback_line()). Outside
# FIML and two-stage estimation (`method` "ml" or "other"), a line names
# robust indices DMAR computed (`dmar_robust`, .sem_robust_mva_line()),
# and one names a fallback from the robust indices (`fell_back`). A
# saturated model (`saturated`, exact fit on zero degrees of freedom) gets
# no FIML line, since its rows hold the values of exact fit in every
# version. A model with no degrees of freedom whose statistic is not zero
# (`zero_df_misfit`) has its cfi alone in a version, so the lines name
# that row alone, and a line says why its tli, nnfi, rmsea, and RMSEA
# limits are NA. Returns NULL when there is nothing to name.
.cfa_k_fit_rows_lines <- function(fiml, robust_indices, fallback, scaled,
                                  complete = FALSE, method = "other",
                                  fell_back = FALSE, dmar_robust = FALSE,
                                  saturated = FALSE, zero_df_misfit = FALSE) {
  lines <- if (saturated) .sem_saturated_text() else character(0)
  rows <- if (zero_df_misfit) "cfi" else "cfi, tli, nnfi, rmsea"
  if (dmar_robust) {
    lines <- c(lines, .sem_robust_mva_line(rows))
  } else if (fell_back && method %in% c("ml", "other")) {
    lines <- c(lines, paste0(
      rows, ": lavaan did not compute its robust ",
      if (zero_df_misfit) "version" else "versions",
      " for this fit, so ",
      if (zero_df_misfit) "this row holds its " else "these rows hold its ",
      fallback, if (zero_df_misfit) " version, which " else " versions, which ",
      if (fallback == "scaled" && zero_df_misfit) {
        "puts the scaled statistic into the standard formula"
      } else if (fallback == "scaled") {
        "put the scaled statistic into the standard formulas"
      } else if (method == "ml" && zero_df_misfit) {
        paste0("estimates the same population value with the bias ",
               "corrected as if the data were normal (Savalei, 2018)")
      } else if (method == "ml") {
        paste0("estimate the same population values with the bias ",
               "corrected as if the data were normal (Savalei, 2018)")
      } else if (zero_df_misfit) {
        "is computed from the standard statistic"
      } else {
        "are computed from the standard statistic"
      }))
  }
  if (fiml && !saturated && !(complete && !scaled)) {
    lines <- c(lines, if (robust_indices && zero_df_misfit) {
      paste0("cfi: the FIML corrected index of Zhang and Savalei (2023), ",
             "version ",
             if (scaled) "3 (no normality assumption)" else "6 (normal theory)",
             ", which approximates the value complete data would give and ",
             "is not computed from the chi_square row")
    } else if (robust_indices) {
      paste0("cfi, rmsea: the FIML corrected indices of Zhang and Savalei ",
             "(2023), version ",
             if (scaled) "3 (no normality assumption)" else "6 (normal theory)",
             ", which approximate the values complete data would give and ",
             "are not computed from the chi_square row; tli, nnfi, and the ",
             "RMSEA interval: lavaan's extension of the same correction")
    } else {
      .cfa_k_fallback_line("FIML", fallback, cfi_only = zero_df_misfit)
    })
  }
  if (zero_df_misfit) {
    lines <- c(lines, paste0(
      "tli, nnfi, rmsea, and the RMSEA interval: NA, since the model has ",
      "zero degrees of freedom and a statistic that is not zero (an ",
      "inequality constraint that binds costs no degree of freedom in ",
      "lavaan's count), so it does not fit exactly, and these indices ",
      "divide by the degrees of freedom; the cfi is formed from the ",
      "statistic"))
  }
  if (length(lines)) lines else NULL
}

# The line that replaces the name of the corrected indices when lavaan
# could not compute them for a FIML or two-stage fit and the rows hold its
# scaled or standard versions instead (the cfi row alone, `cfi_only`, for a
# model with no degrees of freedom that does not fit exactly, whose other
# index rows are NA).
.cfa_k_fallback_line <- function(method, fallback, cfi_only = FALSE) {
  paste0(if (cfi_only) "cfi" else "cfi, tli, nnfi, rmsea",
         ": lavaan could not compute the ", method,
         " corrected ", if (cfi_only) "index" else "indices",
         " of Zhang and Savalei (2023) for this fit, so ",
         if (cfi_only) "this row holds its " else "these rows hold its ",
         fallback, if (cfi_only) " version, which " else " versions, which ",
         if (fallback == "scaled" && cfi_only) {
           "puts the scaled statistic into the standard formula"
         } else if (fallback == "scaled") {
           "put the scaled statistic into the standard formulas"
         } else if (method == "FIML" && cfi_only) {
           "changes with the amount of missing data"
         } else if (method == "FIML") {
           "change with the amount of missing data"
         } else if (cfi_only) {
           "treats the EM estimates as complete data"
         } else {
           "treat the EM estimates as complete data"
         })
}

# The items of each factor whose error variance is estimated below zero
# (for an ordered item under the theta parameterization, the residual
# variance of its latent response), with their communalities
# lambda^2 / (lambda^2 + psi), the share of the item's model implied
# variance the factor accounts for when the factor variance is fixed at
# 1, which a negative psi puts above one. `pe` is lavaan's
# parameterEstimates() of a single group fit. Returns a list named by the
# factors that have such items, each element a numeric vector of the
# communalities named by item; an empty list for a proper solution.
.cfa_k_improper <- function(pe, factors) {
  var_rows <- pe[pe$op == "~~" & pe$lhs == pe$rhs, , drop = FALSE]
  out <- list()
  for (f in names(factors)) {
    v <- factors[[f]]
    psi <- var_rows$est[match(v, var_rows$lhs)]
    load_rows <- pe[pe$op == "=~" & pe$lhs == f, , drop = FALSE]
    lam <- load_rows$est[match(v, load_rows$rhs)]
    bad <- !is.na(psi) & psi < 0
    if (any(bad)) {
      out[[f]] <- stats::setNames(lam[bad]^2 / (lam[bad]^2 + psi[bad]),
                                  v[bad])
    }
  }
  out
}

# Name the classical measurement structure a constraint pattern implies.
# Without a mean structure the intercepts are not modeled, so only the
# "essentially" forms can be claimed (Graham, 2006).
.cfa_k_model_name <- function(eq_load, eq_int, eq_err, has_means) {
  if (!eq_load) {
    extras <- c(if (eq_int) "equal intercepts",
                if (eq_err) "equal error variances")
    if (length(extras) == 0L) return("congeneric (no equality constraints)")
    return(paste0("no classical name (", paste(extras, collapse = " and "),
                  ", free loadings)"))
  }
  if (eq_int && eq_err) {
    "parallel (equal loadings, intercepts, and error variances)"
  } else if (eq_int) {
    "tau-equivalent (equal loadings and intercepts)"
  } else if (eq_err) {
    "essentially parallel (equal loadings and error variances)"
  } else {
    "essentially tau-equivalent (equal loadings)"
  }
}
