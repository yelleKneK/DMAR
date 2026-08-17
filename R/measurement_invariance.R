#' Measurement Invariance Across Groups
#'
#' Fits the standard ladder of multiple group invariance models for a
#' measurement model of any number of factors and reports the comparison
#' table researchers actually use: configural invariance (same pattern, all
#' parameters free by group), metric (equal loadings; required before
#' comparing relations involving the factor), scalar (equal loadings and
#' intercepts; required before comparing factor means), and strict (equal
#' residual variances as well; required before comparing observed-score
#' variances). With ordered indicators a thresholds rung comes first,
#' because thresholds rather than intercepts carry the location information
#' (Wu & Estabrook, 2016); with dichotomous indicators the two are not
#' separately identified, so that rung is folded into metric. Each rung is
#' tested against the previous with a likelihood ratio test, scaled when the
#' estimator in force carries a robust test, which is what declaring ordered
#' items arranges, and the practical-fit changes (delta CFI, delta RMSEA)
#' are reported alongside, since with large samples the chi square will flag
#' trivial differences (Cheung & Rensvold, 2002, suggest delta CFI of about
#' -.01 as a red flag). The fitted \pkg{lavaan} objects come back as an
#' attribute, so score tests and partial invariance refits do not require
#' refitting the ladder. Requires \pkg{lavaan}.
#'
#' @param data A \code{data.frame} with the items and the grouping
#'   variable.
#' @param model The measurement model, given either as \pkg{lavaan} model
#'   syntax (a single string, or a character vector of lines that is
#'   collapsed with newlines, such as
#'   \code{c("visual =~ x1 + x2 + x3", "verbal =~ x4 + x5 + x6")}) or as a
#'   named list mapping each factor name to a character vector of its item
#'   names, from which the syntax is built. Exactly one of \code{model} and
#'   \code{items} must be supplied.
#' @param group Single character string naming the grouping column (two or
#'   more groups).
#' @param items Character vector (three or more) naming the indicator
#'   columns of a single factor. A convenience for the one-factor
#'   (congeneric) case of \code{\link{cfa_1}}, equivalent to
#'   \code{model = "f =~ item_1 + item_2 + ..."}. Exactly one of
#'   \code{model} and \code{items} must be supplied.
#' @param levels Which rungs of the ladder to fit, in order: any leading
#'   subset of the ladder in force. With continuous indicators the ladder is
#'   \code{c("configural", "metric", "scalar", "strict")}; with ordered
#'   indicators it is \code{c("configural", "thresholds", "metric",
#'   "scalar", "strict")}; when every indicator is ordered and dichotomous it
#'   is \code{c("configural", "metric", "scalar")}, because neither the
#'   thresholds nor the strict rung is identified for a two-category item
#'   (see \code{ordered} below). \code{NULL} (the default) fits the whole
#'   ladder in force.
#' @param ordered Ordered-categorical (including binary) items: \code{TRUE}
#'   for every indicator in the model, or a character vector naming the
#'   ordered ones. \code{NULL} (default) treats all indicators as
#'   continuous. Declaring ordered items switches the ladder and, unless the
#'   estimator already carries the mean and variance adjusted test, switches
#'   the estimator (see \code{estimator} below). The number of observed
#'   categories per declared item is counted before the ladder is built,
#'   because a dichotomous item has a single threshold that is not
#'   separately identified from the intercept of its underlying response
#'   (Millsap & Yun-Tein, 2004; Wu & Estabrook, 2016), and its residual
#'   variance is not a free parameter under either parameterization. When
#'   every declared item is dichotomous the thresholds rung is folded into
#'   metric, and the strict rung is dropped as well unless a continuous
#'   indicator is left to carry a residual variance; when only some declared
#'   items are dichotomous, both rungs are kept and a message names those
#'   items, whose thresholds and residual variances the rungs leave untested.
#' @param estimator Estimator passed to \pkg{lavaan}. Defaults to
#'   \code{"ML"}; use \code{"MLR"} for robust corrections, in which case the
#'   likelihood ratio tests use the scaled difference. When \code{ordered}
#'   is supplied, a maximum likelihood or generalized least squares (the \code{"GLS"} label; Browne, 1974)
#'   estimator is replaced by \code{"WLSMV"}, and \code{"DWLS"} and
#'   \code{"ULS"} are replaced by \code{"WLSMV"} and \code{"ULSMV"}. That
#'   second substitution changes only the test: in \pkg{lavaan}
#'   \code{"WLSMV"} is \code{"DWLS"} and \code{"ULSMV"} is \code{"ULS"} with
#'   the mean and variance adjusted statistic, which is what makes the
#'   rung-to-rung difference test asymptotically valid for ordered data. The
#'   one categorical estimator left alone is \code{"WLS"}, the full weight
#'   matrix estimator of Browne (1984) and Muthén (1984), whose statistic is
#'   asymptotically chi square already and whose weight matrix a
#'   substitution would silently change.
#' @param missing Missing data handling passed to \pkg{lavaan}: one of
#'   \code{"listwise"} (the default), \code{"ml"} or \code{"fiml"} (full
#'   information maximum likelihood, the reason to reach for this argument
#'   with continuous items and incomplete cases), or \code{"pairwise"}. Full
#'   information maximum likelihood is not available with the categorical
#'   estimator; asking for both is an error.
#' @param group_partial Character vector of parameters to leave free across
#'   groups at every rung, passed to \pkg{lavaan}'s \code{group.partial}.
#'   Either user-supplied parameter labels or \pkg{lavaan} parameter
#'   specifications such as \code{"visual =~ x2"} or \code{"x2 ~1"}. This is
#'   the partial invariance case of Byrne, Shavelson, and Muthén (1989).
#'   \code{NULL} (default) constrains everything the rung calls for.
#' @param parameterization Identification of ordered indicators, passed to
#'   \pkg{lavaan}: \code{"delta"} (the default, lavaan's) or \code{"theta"}.
#'   Residual variances are free parameters only under \code{"theta"}, so a
#'   strict rung requested with ordered items switches to \code{"theta"}
#'   with a message. Ignored when no item is ordered.
#' @param ... Further arguments passed to \code{\link[lavaan]{cfa}} at every
#'   rung (for example \code{std.lv}, \code{orthogonal}, \code{cluster}).
#'
#' @details
#' The models are nested by construction, each adding equality constraints
#' across groups to the previous. With continuous indicators the constraint
#' sets are none (beyond the configuration), then \code{group.equal =
#' "loadings"}, then \code{c("loadings", "intercepts")}, then
#' \code{c("loadings", "intercepts", "residuals")}.
#'
#' The ordered ladder differs, and the difference is substantive rather than
#' cosmetic. For an ordered indicator the observed response is a coarsening
#' of an underlying continuous response at a set of thresholds, and it is
#' the thresholds, not an intercept, that locate the item on the latent
#' scale. Constraining loadings while leaving thresholds free across groups
#' therefore does not deliver what metric invariance is supposed to deliver,
#' and the accepted sequence (Millsap & Yun-Tein, 2004; Wu & Estabrook,
#' 2016) constrains thresholds first: \code{"thresholds"}, then
#' \code{c("thresholds", "loadings")} for metric, then
#' \code{c("thresholds", "loadings", "intercepts")} for scalar, then adding
#' \code{"residuals"} for strict. The intercept rung is not vacuous even
#' when every indicator is ordered: once thresholds and loadings are
#' constrained, \pkg{lavaan} frees the underlying-response intercepts in the
#' non-reference groups, and the scalar rung is what returns them to zero
#' and lets the latent means be estimated instead. Residual variances,
#' however, are free parameters only under the theta parameterization, so
#' \code{parameterization} switches to \code{"theta"} when a strict rung is
#' requested with ordered items.
#'
#' Dichotomous items are the exception the ordered ladder has to make room
#' for. A two-category item contributes one threshold, and that threshold,
#' the intercept of the underlying response, and its residual variance are
#' not separately identified: the data give one proportion per group per
#' item, which pins down a single standardized location and nothing else
#' (Millsap & Yun-Tein, 2004; Wu & Estabrook, 2016). Constraining thresholds
#' alone across groups is then a reparameterization rather than a
#' restriction, since \pkg{lavaan} frees the underlying-response intercepts
#' by exactly as many parameters as the constraint removes, and the residual
#' variances stay fixed at one in every group whichever parameterization is
#' in force. So the function counts the observed categories of each declared
#' ordered indicator before building the ladder, and a message explains what
#' the count implies. When all of them are dichotomous the thresholds rung is
#' folded into metric, which constrains thresholds and loadings together; the
#' strict rung goes too when every indicator in the model was declared
#' ordered, leaving configural, metric, scalar. A continuous indicator
#' alongside dichotomous ones keeps the strict rung, since its residual
#' variance is a free parameter. When only some declared items are
#' dichotomous the full ordered ladder is kept, since the polytomous items
#' still carry testable thresholds and residual variances, and the message
#' names the dichotomous items so their contribution to those two rungs is
#' not overread. Every rung that is dropped is one that would have cost zero
#' degrees of freedom.
#'
#' When the estimator carries a robust test, the difference between two chi
#' square statistics is not itself chi square distributed.
#' \code{\link[lavaan]{lavTestLRT}} then returns the scaled difference test
#' (the Satorra-Bentler or Satorra correction, chosen by \pkg{lavaan} to
#' match the test in force), and that is what \code{delta_chi_square} and
#' \code{p_value} report; the \code{"test"} attribute names the test used. In
#' that case the model-level \code{chi_square}, \code{p_chi_square},
#' \code{cfi}, and \code{rmsea} columns are lavaan's scaled and robust
#' versions, so \code{delta_chi_square} will not equal the difference of
#' consecutive \code{chi_square} entries. That is a property of scaled
#' difference testing, not an inconsistency.
#'
#' Which estimators carry such a test is worth being precise about, because
#' a naive chi square difference on ordered data is not asymptotically valid.
#' \code{"MLM"}, \code{"MLMV"}, \code{"MLR"}, \code{"WLSM"}, \code{"WLSMV"},
#' \code{"ULSM"}, and \code{"ULSMV"} carry one. \code{"DWLS"} and
#' \code{"ULS"} do not, which is why declaring ordered items promotes them to
#' \code{"WLSMV"} and \code{"ULSMV"}: the discrepancy function and hence the
#' parameter estimates are untouched, and only the statistic changes. The
#' remaining case is \code{"WLS"}, the full weight matrix estimator, whose
#' statistic is asymptotically chi square under the theory of Browne (1984)
#' and Muthén (1984) and so needs no correction; the \code{"test"} attribute
#' reports the standard difference test there. The sample size that theory
#' asks for is large, which is why full weighted least squares is rarely the
#' right choice in practice, but that is a separate matter from whether the
#' difference test is the right one.
#'
#' A rung whose constraints cost no degrees of freedom tests nothing, so
#' \code{delta_chi_square} and \code{p_value} are \code{NA} when
#' \code{delta_df} is zero rather than reporting a difference in chi square
#' that is numerical noise and can come out negative. The situation is
#' reported in a message. It arises with dichotomous items (the ladder above
#' avoids the two rungs where it is structural) and, for instance, with
#' \code{group_partial} specifications that free every parameter a rung would
#' constrain.
#'
#' Failure at a rung does not end the conversation: partial invariance
#' (freeing the offending parameter) is the usual next step, for which
#' \code{group_partial} refits the whole ladder with named parameters free,
#' and the \code{"fits"} attribute gives the fitted objects to
#' \code{\link[lavaan]{lavTestScore}} for a score test of which constraint
#' is doing the damage. This function deliberately reports the standard
#' ladder rather than automating modification searches.
#'
#' @return A tidy wide \code{data.frame} (class \code{dmar_tbl}) with one
#'   row per fitted level and columns \code{level} (label),
#'   \code{chi_square}, \code{df}, \code{p_chi_square} (exact-fit test),
#'   \code{cfi}, \code{rmsea}, and, from the second row on, the
#'   step-comparison columns \code{delta_chi_square}, \code{delta_df},
#'   \code{p_value} (the likelihood ratio test against the previous rung),
#'   \code{delta_cfi}, and \code{delta_rmsea}. \code{delta_chi_square} and
#'   \code{p_value} are \code{NA} wherever \code{delta_df} is zero, since a
#'   difference test on zero degrees of freedom tests nothing. Attributes
#'   carry the non-numeric information: \code{"fits"}, the named list of fitted
#'   \pkg{lavaan} objects, one per level; \code{"estimator"}, the estimator
#'   actually used; \code{"ordered"}, \code{TRUE} when any indicator was
#'   declared ordered; \code{"test"}, naming the chi square difference test
#'   used (\code{NA} when a single rung was fit and nothing was compared);
#'   \code{"fit_indices"}, \code{"standard"} or \code{"robust"}
#'   according to which version of the fit indices is tabled; and
#'   \code{"model"}, the \pkg{lavaan} syntax that was fitted.
#'
#' @references
#' Browne, M. W. (1974). Generalized least squares estimators in the
#'   analysis of covariance structures. \emph{South African Statistical
#'   Journal, 8}, 1--24.
#'
#' Browne, M. W. (1984). Asymptotically distribution-free methods for the
#'   analysis of covariance structures. \emph{British Journal of
#'   Mathematical and Statistical Psychology, 37}(1), 62--83.
#'
#' Byrne, B. M., Shavelson, R. J., & Muthén, B. (1989). Testing for the
#'   equivalence of factor covariance and mean structures: The issue of
#'   partial measurement invariance. \emph{Psychological Bulletin, 105}(3),
#'   456--466.
#'
#' Cheung, G. W., & Rensvold, R. B. (2002). Evaluating goodness-of-fit
#'   indexes for testing measurement invariance. \emph{Structural Equation
#'   Modeling, 9}(2), 233--255. \doi{10.1207/S15328007SEM0902_5}
#'
#' Meredith, W. (1993). Measurement invariance, factor analysis and
#'   factorial invariance. \emph{Psychometrika, 58}(4), 525--543.
#'   \doi{10.1007/BF02294825}
#'
#' Millsap, R. E. (2011). \emph{Statistical approaches to measurement
#'   invariance}. Routledge.
#'
#' Millsap, R. E., & Yun-Tein, J. (2004). Assessing factorial invariance in
#'   ordered-categorical measures. \emph{Multivariate Behavioral Research,
#'   39}(3), 479--515. \doi{10.1207/s15327906mbr3903_4}
#'
#' Muthén, B. (1984). A general structural equation model with dichotomous,
#'   ordered categorical, and continuous latent variable indicators.
#'   \emph{Psychometrika, 49}(1), 115--132.
#'
#' Satorra, A., & Bentler, P. M. (2001). A scaled difference chi-square test
#'   statistic for moment structure analysis. \emph{Psychometrika, 66}(4),
#'   507--514. \doi{10.1007/BF02296192}
#'
#' Wu, H., & Estabrook, R. (2016). Identification of confirmatory factor
#'   analysis models of different levels of invariance for ordered
#'   categorical outcomes. \emph{Psychometrika, 81}(4), 1014--1045.
#'   \doi{10.1007/s11336-016-9506-0}
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso \code{\link{cfa_1}} and \code{\link{cfa_k}} for the
#'   single-group measurement models; \code{\link{reliability_omega}} for
#'   the reliability of the composite the model justifies;
#'   \code{\link{compare_cov_structures}} for covariance-structure
#'   comparisons outside the factor model;
#'   \code{\link[lavaan]{lavTestScore}} for the score test that localizes a
#'   failed rung.
#'
#' @family multivariate and latent variable methods
#'
#' @keywords models multivariate
#'
#' @examples
#' # Do the two schools measure the verbal and the reasoning construct
#' # in the same way? (Holzinger & Swineford, bundled.) The measurement
#' # model is a named list of factors; for a single factor, name its
#' # indicators with items = instead.
#' data(holzinger_swineford)
#' hs_factors <- list(
#'   verbal    = c("t6_paragraph_comprehension", "t7_sentence",
#'                 "t9_word_meaning"),
#'   deduction = c("t20_deduction", "t22_problem_reasoning",
#'                 "t23_series_completion"))
#'
#' # The bottom two rungs. Configural invariance asks whether the same
#' # pattern of loadings holds at both schools; metric adds the
#' # constraint that the loadings are equal across schools, and it is
#' # the rung that has to hold before a relation involving one of these
#' # factors is compared across them. Naming those two in levels = stops
#' # the ladder there, which is what keeps this example quick.
#' mi <- measurement_invariance(holzinger_swineford, hs_factors,
#'                              group = "school",
#'                              levels = c("configural", "metric"))
#' mi
#'
#' # The fitted models travel with the table, so localizing a failed rung
#' # costs no refitting.
#' names(attr(mi, "fits"))
#'
#' # The broom verbs: one row per rung of the ladder, and the model-level
#' # summary (estimator, test flavor, fit index flavor).
#' generics::tidy(mi)
#' generics::glance(mi)
#'
#' # Each of the calls below refits the ladder from the bottom, so they
#' # are shown here rather than run. Leaving levels = at its default
#' # fits the whole ladder, configural through metric, scalar, and
#' # strict:
#' # measurement_invariance(holzinger_swineford, hs_factors,
#' #                        group = "school")
#' #
#' # Partial invariance frees one loading across the schools at every
#' # rung (Byrne, Shavelson, & Muthén, 1989), so each constrained rung
#' # loses one degree of freedom relative to full invariance:
#' # measurement_invariance(holzinger_swineford, hs_factors,
#' #                        group = "school",
#' #                        group_partial = "verbal =~ t7_sentence")
#' #
#' # With the indicators coded as ordered categories, ordered = TRUE puts
#' # the thresholds rung first, because thresholds rather than intercepts
#' # carry the location information there, and makes the rung-to-rung
#' # tests the scaled difference tests:
#' # hs_ordered <- holzinger_swineford
#' # for (item in unlist(hs_factors, use.names = FALSE)) {
#' #   hs_ordered[[item]] <- as.integer(cut(
#' #     holzinger_swineford[[item]],
#' #     breaks = quantile(holzinger_swineford[[item]],
#' #                       c(0, .25, .5, .75, 1)),
#' #     include.lowest = TRUE))
#' # }
#' # mi_ordered <- measurement_invariance(hs_ordered, hs_factors,
#' #                                      group = "school",
#' #                                      ordered = TRUE)
#' # mi_ordered
#' # attr(mi_ordered, "test")
#' #
#' # The measurement model can also be given as lavaan syntax, and
#' # missing = "fiml" fits the ladder by full information maximum
#' # likelihood when continuous items have incomplete cases:
#' # hs_missing <- holzinger_swineford
#' # set.seed(113)
#' # hs_missing$t7_sentence[sample(nrow(hs_missing), 20)] <- NA
#' # measurement_invariance(
#' #   hs_missing,
#' #   model = "verbal =~ t6_paragraph_comprehension + t7_sentence +
#' #            t9_word_meaning",
#' #   group = "school", missing = "fiml")
#'
#' @export
measurement_invariance <- function(data, model = NULL, group,
                                   items = NULL,
                                   levels = NULL,
                                   ordered = NULL,
                                   estimator = "ML",
                                   missing = "listwise",
                                   group_partial = NULL,
                                   parameterization = c("delta", "theta"),
                                   ...) {
  if (!requireNamespace("lavaan", quietly = TRUE)) {
    stop("Package 'lavaan' is required for measurement_invariance(). ",
         "Install it with install.packages(\"lavaan\").", call. = FALSE)
  }
  if (!is.data.frame(data)) {
    stop("'data' must be a data.frame.", call. = FALSE)
  }
  parameterization <- match.arg(parameterization)

  # ---- the measurement model: syntax, a named list, or one factor ------
  if (is.null(model) == is.null(items)) {
    stop("Supply exactly one of 'model' (lavaan syntax or a named list of ",
         "factors) and 'items' (the indicators of a single factor).",
         call. = FALSE)
  }
  if (!is.null(items)) {
    if (!is.character(items) || length(items) < 3L ||
        !all(items %in% names(data))) {
      stop("'items' must name three or more columns of 'data'.",
           call. = FALSE)
    }
    syntax <- paste("f =~", paste(items, collapse = " + "))
  } else if (is.list(model)) {
    syntax <- .measurement_invariance_syntax(model, data)
  } else if (is.character(model) && length(model) >= 1L &&
             !anyNA(model)) {
    syntax <- paste(model, collapse = "\n")
    if (!grepl("=~", syntax, fixed = TRUE)) {
      stop("'model' must contain at least one measurement relation (=~).",
           call. = FALSE)
    }
  } else {
    stop("'model' must be lavaan model syntax (a character string) or a ",
         "named list: list(factor_name = c(\"item_1\", \"item_2\", ...)).",
         call. = FALSE)
  }

  # Indicators of the fitted model, needed for ordered = TRUE and to check
  # that every indicator is actually in the data.
  parsed <- lavaan::lavaanify(syntax)
  latent <- unique(parsed$lhs[parsed$op == "=~"])
  indicators <- setdiff(unique(parsed$rhs[parsed$op == "=~"]), latent)
  absent <- setdiff(indicators, names(data))
  if (length(absent)) {
    stop("Indicators not found in 'data': ", paste(absent, collapse = ", "),
         ".", call. = FALSE)
  }

  # ---- the grouping variable -------------------------------------------
  if (!is.character(group) || length(group) != 1L ||
      !(group %in% names(data))) {
    stop("'group' must name one column of 'data'.", call. = FALSE)
  }
  if (length(unique(stats::na.omit(data[[group]]))) < 2L) {
    stop("'group' must have at least two groups.", call. = FALSE)
  }

  # ---- estimator, missing data, ordered indicators ---------------------
  estimator_choices <- c("ML", "MLM", "MLMV", "MLR", "GLS", "WLS", "DWLS",
                         "WLSM", "WLSMV", "ULS", "ULSM", "ULSMV")
  missing_choices <- c("listwise", "ml", "fiml", "pairwise")
  if (!is.character(estimator) || length(estimator) != 1L ||
      !(estimator %in% estimator_choices)) {
    stop("'estimator' must be one of ",
         paste(dQuote(estimator_choices), collapse = ", "), ".",
         call. = FALSE)
  }
  if (!is.character(missing) || length(missing) != 1L ||
      !(missing %in% missing_choices)) {
    stop("'missing' must be one of ",
         paste(dQuote(missing_choices), collapse = ", "), ".",
         call. = FALSE)
  }

  ordered_items <- NULL
  if (!is.null(ordered)) {
    if (isTRUE(ordered)) {
      ordered_items <- indicators
    } else if (is.character(ordered) && length(ordered)) {
      unknown <- setdiff(ordered, indicators)
      if (length(unknown)) {
        stop("'ordered' names variables that are not indicators of the ",
             "model: ", paste(unknown, collapse = ", "), ".", call. = FALSE)
      }
      ordered_items <- unique(ordered)
    } else {
      stop("'ordered' must be TRUE (all indicators) or a character vector ",
           "of item names.", call. = FALSE)
    }
    if (missing %in% c("ml", "fiml")) {
      stop("Full information maximum likelihood is not available with the ",
           "categorical estimator; with ordered items use ",
           "missing = \"pairwise\" or missing = \"listwise\".",
           call. = FALSE)
    }
    if (estimator %in% c("ML", "MLM", "MLMV", "MLR", "GLS")) {
      message("Ordered items declared: switching to estimator = \"WLSMV\".")
      estimator <- "WLSMV"
    } else if (estimator %in% c("DWLS", "ULS")) {
      # Neither carries a robust test, so the rung-to-rung difference test
      # would be the unscaled one, which is not asymptotically valid for
      # ordered data. The promotion changes only the test statistic: in
      # lavaan "WLSMV" is "DWLS" and "ULSMV" is "ULS" with the mean and
      # variance adjustment, so the discrepancy function, and with it every
      # parameter estimate, is untouched. "WLS" is deliberately left alone
      # (see the estimator documentation).
      promoted <- c(DWLS = "WLSMV", ULS = "ULSMV")[[estimator]]
      message("Ordered items declared: switching to estimator = \"",
              promoted, "\". That is the same ", estimator,
              " fit with the mean and variance adjusted test, which the ",
              "chi square difference test between rungs needs to be ",
              "asymptotically valid with ordered data; the parameter ",
              "estimates are unchanged.")
      estimator <- promoted
    }
  }
  is_ordered <- !is.null(ordered_items)

  # A dichotomous indicator carries one threshold, and that threshold is not
  # separately identified from the intercept of its underlying response, nor
  # is its residual variance a free parameter (Millsap & Yun-Tein, 2004; Wu &
  # Estabrook, 2016). Counting the observed categories here is what keeps the
  # ladder from including rungs that cannot be tested. The strict rung goes
  # only when no indicator is left to carry a free residual variance, since a
  # continuous indicator alongside dichotomous ones still has one.
  binary_items <- character(0)
  all_binary <- FALSE
  drop_strict <- FALSE
  if (is_ordered) {
    n_categories <- .measurement_invariance_n_categories(data, ordered_items)
    binary_items <- ordered_items[n_categories <= 2L]
    all_binary <- length(binary_items) == length(ordered_items)
    drop_strict <- all_binary && setequal(ordered_items, indicators)
  }

  if (!is.null(group_partial) &&
      (!is.character(group_partial) || !length(group_partial) ||
       anyNA(group_partial))) {
    stop("'group_partial' must be a character vector of parameter labels ",
         "or lavaan parameter specifications, such as ",
         "\"visual =~ item_2\".", call. = FALSE)
  }

  # ---- the ladder in force ---------------------------------------------
  # With every declared item dichotomous the thresholds rung is folded into
  # metric (which constrains thresholds and loadings together), and the
  # strict rung goes too when no indicator is left to carry a free residual
  # variance. Each dropped rung would have cost zero degrees of freedom.
  ladder <- if (!is_ordered) {
    c("configural", "metric", "scalar", "strict")
  } else {
    rungs <- c("configural", "thresholds", "metric", "scalar", "strict")
    if (all_binary) rungs <- setdiff(rungs, "thresholds")
    if (drop_strict) rungs <- setdiff(rungs, "strict")
    rungs
  }
  ladder_label <- if (!is_ordered) {
    "continuous"
  } else if (all_binary) {
    "dichotomous ordered"
  } else {
    "ordered"
  }
  if (all_binary) {
    message("Every declared ordered indicator is dichotomous. A threshold ",
            "and the intercept of the underlying response are then not ",
            "separately identified (Millsap & Yun-Tein, 2004; Wu & ",
            "Estabrook, 2016), so the thresholds rung is folded into metric, ",
            "which constrains thresholds and loadings together. ",
            if (drop_strict) {
              paste0("The strict rung is dropped as well, since a ",
                     "dichotomous indicator has no free residual variance to ",
                     "constrain. ")
            } else {
              paste0("The strict rung is kept, since it still constrains the ",
                     "residual variances of the indicators that were not ",
                     "declared ordered. ")
            },
            "The ladder is ", paste(ladder, collapse = ", "), ".")
  } else if (length(binary_items)) {
    message(length(binary_items), " of ", length(ordered_items),
            " declared ordered indicators are dichotomous (",
            paste(binary_items, collapse = ", "),
            "). The full ordered ladder is kept, because the polytomous ",
            "items carry testable thresholds and residual variances, but ",
            "the thresholds and strict rungs test nothing about the ",
            "dichotomous ones.")
  }
  constraints <- if (is_ordered) {
    list(configural = character(0),
         thresholds = "thresholds",
         metric     = c("thresholds", "loadings"),
         scalar     = c("thresholds", "loadings", "intercepts"),
         strict     = c("thresholds", "loadings", "intercepts",
                        "residuals"))
  } else {
    list(configural = character(0),
         metric     = "loadings",
         scalar     = c("loadings", "intercepts"),
         strict     = c("loadings", "intercepts", "residuals"))
  }
  if (is.null(levels)) {
    levels <- ladder
  } else {
    if (!is.character(levels) || !length(levels) || anyNA(levels)) {
      stop("'levels' must be a character vector of rung names.",
           call. = FALSE)
    }
    unknown <- setdiff(levels, ladder)
    if (length(unknown)) {
      stop("'levels' names rungs that are not on the ", ladder_label,
           " ladder: ", paste(unknown, collapse = ", "),
           ". The ladder in force is ", paste(ladder, collapse = ", "), ".",
           if (all_binary)
             paste0(" With every declared ordered indicator dichotomous, a",
                    " rung missing from that list would impose no",
                    " restriction; see ?measurement_invariance."),
           call. = FALSE)
    }
    levels <- ladder[ladder %in% levels]
    if (!identical(levels, ladder[seq_along(levels)])) {
      stop("'levels' must be a leading subset of the ladder in force: ",
           paste(ladder, collapse = ", "), ".", call. = FALSE)
    }
  }

  # Residual variances are estimated only under the theta parameterization,
  # so a strict rung on ordered items would otherwise cost zero degrees of
  # freedom and test nothing.
  if (is_ordered && "strict" %in% levels && parameterization == "delta") {
    message("A strict rung with ordered items requires the theta ",
            "parameterization (residual variances are not free parameters ",
            "under delta): switching to parameterization = \"theta\".")
    parameterization <- "theta"
  }

  # ---- fit the ladder ---------------------------------------------------
  dots <- list(...)
  fits <- lapply(levels, function(lv) {
    args <- list(model = syntax, data = data, group = group,
                 estimator = estimator, missing = missing,
                 meanstructure = TRUE)
    if (is_ordered) {
      args$ordered <- ordered_items
      args$parameterization <- parameterization
    }
    if (length(constraints[[lv]])) args$group.equal <- constraints[[lv]]
    if (!is.null(group_partial)) args$group.partial <- group_partial
    do.call(lavaan::cfa, c(args, dots))
  })
  names(fits) <- levels

  # ---- fit indices ------------------------------------------------------
  # With a robust or categorical estimator the standard statistics are not
  # the ones to report; lavaan's scaled and robust versions are.
  scaled <- !identical(
    as.character(lavaan::lavInspect(fits[[1L]], "options")$test), "standard")
  wanted <- if (scaled) {
    c(chi_square = "chisq.scaled", df = "df.scaled",
      p_chi_square = "pvalue.scaled", cfi = "cfi.robust",
      rmsea = "rmsea.robust")
  } else {
    c(chi_square = "chisq", df = "df", p_chi_square = "pvalue",
      cfi = "cfi", rmsea = "rmsea")
  }
  fallback <- c(chi_square = "chisq", df = "df", p_chi_square = "pvalue",
                cfi = "cfi", rmsea = "rmsea")
  available <- names(lavaan::fitMeasures(fits[[1L]]))
  wanted[!(wanted %in% available)] <- fallback[!(wanted %in% available)]

  M <- t(vapply(fits, function(fit)
    as.numeric(lavaan::fitMeasures(fit, wanted)), numeric(length(wanted))))

  # ---- the rung-to-rung comparisons -------------------------------------
  k <- length(levels)
  d_chi <- d_df <- p_lrt <- d_cfi <- d_rmsea <- rep(NA_real_, k)
  # With a single rung there is nothing to compare, so no difference test
  # was performed and the attribute says so rather than naming one.
  test_label <- NA_character_
  untestable <- character(0)
  if (k > 1L) {
    for (i in 2:k) {
      # A rung that adds no degree of freedom makes lavTestLRT warn that two
      # models share their degrees of freedom. The package says what that
      # means in its own words below, so the warning is muffled here rather
      # than reaching the user twice, once cryptically. Every other lavaan
      # warning passes through untouched.
      lrt <- withCallingHandlers(
        lavaan::lavTestLRT(fits[[i - 1L]], fits[[i]]),
        warning = function(cond) {
          if (grepl("same degrees of freedom", conditionMessage(cond),
                    fixed = TRUE))
            invokeRestart("muffleWarning")
        })
      d_df[i] <- lrt[2L, "Df diff"]
      # Zero degrees of freedom means the rung imposed no restriction: the
      # two models are reparameterizations of one another, and their chi
      # square difference is numerical noise that can come out negative. A
      # difference test is reported only where one was performed.
      if (!is.na(d_df[i]) && d_df[i] <= 0) {
        untestable <- c(untestable, levels[i])
      } else {
        d_chi[i] <- lrt[2L, "Chisq diff"]
        p_lrt[i] <- lrt[2L, "Pr(>Chisq)"]
      }
      d_cfi[i]   <- M[i, 4L] - M[i - 1L, 4L]
      d_rmsea[i] <- M[i, 5L] - M[i - 1L, 5L]
      if (i == 2L) test_label <- .measurement_invariance_test_label(lrt,
                                                                    scaled)
    }
    if (length(untestable)) {
      message("The constraints at the ", paste(untestable, collapse = ", "),
              if (length(untestable) > 1L) " rungs" else " rung",
              " cost no degrees of freedom, so nothing was tested there: ",
              "delta_chi_square and p_value are NA rather than a difference ",
              "test on zero degrees of freedom.")
    }
  }

  out <- data.frame(
    level            = levels,
    chi_square       = M[, 1L],
    df               = M[, 2L],
    p_chi_square     = M[, 3L],
    cfi              = M[, 4L],
    rmsea            = M[, 5L],
    delta_chi_square = d_chi,
    delta_df         = d_df,
    p_value          = p_lrt,
    delta_cfi        = d_cfi,
    delta_rmsea      = d_rmsea,
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  res <- .as_dmar_tbl(out, subclass = "dmar_measurement_invariance")
  attr(res, "fits") <- fits
  attr(res, "estimator") <- estimator
  attr(res, "ordered") <- is_ordered
  attr(res, "test") <- test_label
  attr(res, "fit_indices") <- if (scaled) "robust" else "standard"
  attr(res, "model") <- syntax
  res
}


# Observed categories per declared ordered indicator. The ladder depends on
# the count: a two-category item has one threshold, which is not separately
# identified from the intercept of its underlying response, and no free
# residual variance, so a thresholds rung and a strict rung would impose no
# restriction on it. Levels a factor carries but the data never realize are
# not categories, hence the count of distinct non-missing values rather than
# nlevels().
.measurement_invariance_n_categories <- function(data, items) {
  vapply(items, function(item) {
    x <- data[[item]]
    length(unique(x[!is.na(x)]))
  }, integer(1L))
}


# Builds lavaan measurement syntax from a named list of factors. Kept here
# rather than shared with cfa_k() because cfa_k() also has to encode its
# within-factor equality descriptors; this one only needs the pattern.
.measurement_invariance_syntax <- function(factors, data) {
  if (!length(factors) || is.null(names(factors)) ||
      any(names(factors) == "")) {
    stop("A list 'model' must be named: ",
         "list(factor_name = c(\"item_1\", \"item_2\", ...)).",
         call. = FALSE)
  }
  factor_names <- names(factors)
  if (anyDuplicated(factor_names)) {
    stop("Factor names must be unique.", call. = FALSE)
  }
  if (!all(factor_names == make.names(factor_names))) {
    stop("Factor names must be syntactically valid R names (they become ",
         "lavaan labels); rename, for example, 'my factor' to 'my_factor'.",
         call. = FALSE)
  }
  for (f in factor_names) {
    if (!is.character(factors[[f]]) || length(factors[[f]]) < 2L) {
      stop(sprintf("Factor '%s' must name two or more items.", f),
           call. = FALSE)
    }
  }
  items <- unlist(factors, use.names = FALSE)
  if (anyDuplicated(items)) {
    stop("In the list form each item loads on exactly one factor (simple ",
         "structure); duplicated: ",
         paste(unique(items[duplicated(items)]), collapse = ", "),
         ". Give 'model' as lavaan syntax for cross-loadings.",
         call. = FALSE)
  }
  if (length(factors) == 1L && length(factors[[1L]]) < 3L) {
    stop("A single factor needs three or more items for identification.",
         call. = FALSE)
  }
  absent <- setdiff(items, names(data))
  if (length(absent)) {
    stop("Items not found in 'data': ", paste(absent, collapse = ", "), ".",
         call. = FALSE)
  }
  paste(vapply(factor_names, function(f)
    paste(f, "=~", paste(factors[[f]], collapse = " + ")), character(1L)),
    collapse = "\n")
}


# Names the chi square difference test lavTestLRT() actually performed. The
# scaling method is chosen by lavaan to match the test in force, so it is
# read back from the returned object's heading rather than assumed.
.measurement_invariance_test_label <- function(lrt, scaled) {
  if (!scaled) return("standard chi square difference test")
  heading <- paste(attr(lrt, "heading"), collapse = " ")
  heading <- gsub("[\r\n]+", " ", heading)
  # lavaan's heading may quote the method with straight or curly quotes
  # depending on the locale and the fancy-quotes option, so both are matched.
  # The curly quotes are written as escapes because R code in a portable
  # package must be ASCII.
  pattern <- "method[[:space:]]*=[[:space:]]*[\"\u201c]([^\"\u201d]+)[\"\u201d]"
  method <- if (grepl(pattern, heading)) {
    sub(paste0(".*", pattern, ".*"), "\\1", heading)
  } else {
    NA_character_
  }
  if (is.na(method)) {
    "scaled chi square difference test"
  } else {
    paste0("scaled chi square difference test (method = \"", method, "\")")
  }
}
