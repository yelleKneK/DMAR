#' Multiple-Factor Confirmatory Factor Analysis Model
#'
#' Fits a confirmatory factor analysis model with one or more factors,
#' where each factor is specified by naming its indicator variables and
#' the measurement structure is specified by describing what is
#' constrained (\code{equal_loading}, \code{equal_intercept},
#' \code{equal_error}) rather than by more technical terms. The function then
#' reports which classical measurement structure the description implies
#' (congeneric, essentially tau-equivalent, tau-equivalent, essentially
#' parallel, or parallel), the parameter estimates with confidence
#' intervals, fit information, and, per factor, coefficient omega, the
#' average variance extracted (AVE), and coefficient \emph{H}, each with a
#' delta method standard error and confidence interval computed by
#' \pkg{lavaan} from defined parameters (no additional packages are
#' involved).
#'
#' With \code{ordered} items, the model is fit by WLSMV to polychoric
#' correlations with thresholds. A sum score of ordered items lives on
#' the metric of the observed categories, not on the latent response
#' metric of the polychoric loadings, so for a factor whose items are
#' ordered the reported omega is the Green and Yang (2009) categorical
#' sum score omega computed from the same fit; the substitution is
#' announced in a message and recorded in the \code{omega_metric}
#' attribute, and the delta method interval columns are \code{NA} for
#' those rows because the delta method interval describes the latent
#' response metric. The categorical sum score omega is the coefficient
#' Kelley and Pornprasertmanit (2016) call categorical omega. AVE and
#' coefficient \emph{H} concern the latent
#' response correlations themselves and are reported unchanged on that
#' metric. For a bootstrap confidence interval on a categorical omega,
#' use \code{\link{reliability_omega_categorical}} on the factor's
#' items.
#'
#' @param data A raw data matrix or data frame, rows are respondents and
#'   columns include the items named in \code{factors}. Supply exactly
#'   one of \code{data} or \code{S}.
#' @param S A symmetric covariance matrix of the items, with dimnames
#'   naming the items; \code{N} is then required. Supply exactly one of
#'   \code{data} or \code{S}.
#' @param factors Named list. Each element names a factor and gives the
#'   character vector of its indicator columns (two or more per factor;
#'   three or more when only one factor is specified). Each item loads on
#'   exactly one factor (simple structure).
#' @param N Total sample size. Required with \code{S}; ignored (inferred
#'   from the rows) with \code{data}.
#' @param M Optional named numeric vector of item means, used with a
#'   covariance matrix to model the mean structure (required there when
#'   \code{equal_intercept} is used). Ignored for raw data.
#' @param equal_loading Logical, or a named logical vector with one
#'   element per factor. \code{TRUE} constrains the loadings within a
#'   factor to a single value (across factors nothing is equated).
#'   Defaults to \code{FALSE}.
#' @param equal_intercept Logical, or a named logical vector with one
#'   element per factor. \code{TRUE} constrains the item intercepts
#'   within a factor to a single value; this requires the mean structure
#'   (raw data, or \code{M} with a covariance matrix). Defaults to
#'   \code{FALSE}.
#' @param equal_error Logical, or a named logical vector with one element
#'   per factor. \code{TRUE} constrains the error variances within a
#'   factor to a single value. Defaults to \code{FALSE}.
#' @param correlated_factors Logical. If \code{TRUE} (default) the
#'   factors covary freely; because each factor variance is fixed to 1,
#'   the \code{phi} terms for factor pairs are the latent correlations.
#'   If \code{FALSE} the factor covariances are fixed to zero.
#' @param meanstructure Logical or \code{NULL}. \code{NULL} (default)
#'   models the mean structure exactly when it is needed: when any
#'   \code{equal_intercept} is \code{TRUE} or when \code{M} is supplied.
#'   Set \code{TRUE} to model intercepts regardless (raw data or \code{M}
#'   required), or \code{FALSE} to suppress them (an error if
#'   \code{equal_intercept} is used).
#' @param estimator Character; estimator passed to \pkg{lavaan}. Must be
#'   one of \code{"ML"} (default; maximum likelihood, fully efficient
#'   under multivariate normality), \code{"MLR"} (robust maximum
#'   likelihood: maximum likelihood estimates with standard errors and
#'   test statistic corrected for nonnormality; Satorra & Bentler,
#'   1994), \code{"WLS"} (the asymptotic distribution free estimator of
#'   Browne, 1984; raw data and a large sample required),
#'   \code{"WLSMV"} (diagonally weighted least squares with mean- and
#'   variance-adjusted test statistic; Muthén, 1984; Muthén, du Toit, &
#'   Spisic, 1997; the standard choice for ordered categorical items,
#'   and what \code{ordered} switches to), or \code{"GLS"}
#'   (generalized least squares; Browne, 1974). With a robust estimator
#'   the reported fit indices are the robust versions.
#' @param missing Character; missing-data handling passed to \pkg{lavaan}
#'   when raw data are supplied. Common values are \code{"listwise"}
#'   (default, listwise deletion) and \code{"ml"} (full information
#'   maximum likelihood). Ignored with \code{S}.
#'   With \code{ordered} items, \code{"ml"}/\code{"fiml"} are not
#'   available; use \code{"pairwise"} or \code{"listwise"}.
#' @param ordered Ordered-categorical items: \code{NULL} (none, the
#'   default), \code{TRUE} (every item), or a character vector of item
#'   names. Requires raw data. Each factor must be all ordered or all
#'   continuous. Declaring ordered items switches the estimator to
#'   \code{"WLSMV"} (with a message, unless a categorical estimator was
#'   requested), fits thresholds in place of intercepts, and reports
#'   each ordered factor's omega on the categorical sum score metric via
#'   the Green and Yang (2009) computation; see \emph{Details}.
#' @param se Standard error type passed to \pkg{lavaan}; see
#'   \code{\link[lavaan]{cfa}}. Common values are \code{"standard"}
#'   (default), \code{"robust.sem"} (with \code{estimator = "MLR"}), and
#'   \code{"none"} (point estimates only; fastest).
#' @param conf_level Confidence level for the parameter confidence
#'   intervals, including the delta method intervals for omega, AVE, and
#'   \emph{H}. Defaults to 0.95. The RMSEA interval is a separate
#'   convention (see Details).
#' @param output Format of the returned object:
#'   \describe{
#'     \item{\code{"verbose"}}{(default) Parameter estimates with
#'       confidence intervals, the per-factor defined parameters
#'       (\code{loading_sum}, \code{error_sum}, \code{omega}, \code{ave},
#'       \code{H}), and fit information.}
#'     \item{\code{"measurement"}}{The measurement-property rows only:
#'       per factor \code{omega}, \code{ave}, and \code{H} (with delta
#'       method standard errors and confidence intervals), the latent
#'       correlation \code{phi} for every factor pair (with its
#'       confidence interval), and, for raw data, the
#'       heterotrait-monotrait ratio \code{htmt} for every factor pair
#'       via \code{\link{htmt}}.}
#'     \item{\code{"summary"}}{The raw \code{summary()} output from
#'       \pkg{lavaan} (not a data frame).}
#'     \item{\code{"standardized"}}{The standardized parameter estimates
#'       from \code{lavaan::standardizedSolution()}.}
#'     \item{\code{"fit"}}{The raw \pkg{lavaan} fit object. Escape hatch
#'       for direct lavaan access, including likelihood ratio tests
#'       between two \code{cfa_k()} fits via
#'       \code{lavaan::lavTestLRT()}.}
#'   }
#' @param \dots Additional arguments forwarded to
#'   \code{\link[lavaan]{lavaan}} (e.g., \code{group}, \code{cluster},
#'   \code{bootstrap}).
#'
#' @details
#' \code{\link{cfa_1}} and \code{\link{cfa_2}} are convenience wrappers
#' around this function for the one and two factor cases: \code{cfa_1()}
#' takes a vector of items and fits one factor over them, and
#' \code{cfa_2()} takes the items of each of two factors. Both forward
#' every argument here, so their results are this function's results,
#' with the factors named \code{f1} (and \code{f2}).
#'
#' \strong{Describing the model instead of naming it.} The classical
#' measurement structures are nested patterns of within-factor equality
#' constraints (Lord & Novick, 1968; Graham, 2006):
#' \describe{
#'   \item{congeneric}{loadings, intercepts, and error variances all
#'     free.}
#'   \item{essentially tau-equivalent}{equal loadings; intercepts and
#'     error variances free.}
#'   \item{tau-equivalent}{equal loadings and equal intercepts; error
#'     variances free.}
#'   \item{essentially parallel}{equal loadings and equal error
#'     variances; intercepts free.}
#'   \item{parallel}{equal loadings, equal intercepts, and equal error
#'     variances.}
#' }
#' The caller states the constraints; the function reports the implied
#' name, per factor, in the printed header and in the \code{"model"}
#' attribute of the returned table. The distinction between
#' tau-equivalent and essentially tau-equivalent (and between parallel
#' and essentially parallel) lives entirely in the mean structure: the
#' covariance structure of the two members of each pair is identical, so
#' without intercepts in the model only the "essentially" form can be
#' claimed. That is why \code{equal_intercept} requires raw data or
#' \code{M}: covariances alone cannot speak to it. A constraint pattern
#' outside the classical list (for example equal error variances with
#' free loadings) is fit as requested and labeled descriptively, since it
#' has no conventional name.
#'
#' Identification fixes each factor variance to 1 (and each factor mean
#' to 0 when the mean structure is modeled), so all loadings are
#' estimated and within-factor equality constraints are meaningful. Two
#' \code{cfa_k()} fits that differ only in descriptor settings are nested,
#' so \code{output = "fit"} feeds \code{lavaan::lavTestLRT()} directly
#' (with \code{estimator = "MLR"}, lavaan applies the scaled difference
#' test).
#'
#' \strong{Measurement properties.} For factor \eqn{f} with
#' unstandardized loadings \eqn{\lambda_j} and error variances
#' \eqn{\psi_j} (factor variance 1):
#' coefficient omega
#' \eqn{= (\sum_j \lambda_j)^2 / ((\sum_j \lambda_j)^2 + \sum_j \psi_j)}
#' (McDonald, 1999), the reliability of the unit-weighted composite;
#' the average variance extracted
#' \eqn{= J^{-1} \sum_j \lambda_j^2 / (\lambda_j^2 + \psi_j)}
#' (Fornell & Larcker, 1981), the mean proportion of item variance the
#' factor accounts for; and coefficient
#' \eqn{H = (1 + (\sum_j \lambda_j^2/\psi_j)^{-1})^{-1}}
#' (Hancock & Mueller, 2001), the reliability of the optimally weighted
#' composite, which no single item can drag below its value for any
#' subset. All three are computed as \pkg{lavaan} defined parameters, so
#' each carries a delta method standard error and a \code{conf_level}
#' confidence interval in the same table as the model parameters.
#'
#' \strong{Discriminant validity.} Three complementary readings come from
#' \code{output = "measurement"}: (a) the latent correlation \code{phi}
#' for a factor pair, with a confidence interval whose upper limit near 1
#' means the data cannot distinguish the two factors; (b) the Fornell and
#' Larcker (1981) comparison, which asks whether each factor's \code{ave}
#' exceeds the squared \code{phi} of its pairs (the factor should share
#' more variance with its own items than with the other factor); and (c)
#' for raw data, the model-free \code{htmt} ratio (Henseler, Ringle, &
#' Sarstedt, 2015). The rows report the numbers and their uncertainty;
#' the judgment is the researcher's.
#'
#' \strong{Confidence interval conventions.} Parameter rows (including
#' omega, \code{ave}, \emph{H}, and \code{phi}) use \code{conf_level}.
#' The RMSEA interval follows its own convention: the
#' \code{rmsea_ci_level} row records the level actually used (0.90, the
#' conventional level for RMSEA, as in \pkg{lavaan}), and
#' \code{rmsea_ci_lower} / \code{rmsea_ci_upper} are that interval.
#'
#' Common row names under \code{term}: \code{lambda_<factor>_<j>}
#' (loadings; \code{lambda_<factor>} when equated),
#' \code{psi_<factor>_<j>} (error variances; \code{psi_<factor>} when
#' equated), \code{nu_<factor>_<j>} (intercepts, when the mean structure
#' is modeled; \code{nu_<factor>} when equated), \code{phi_<factor>}
#' (factor variance, fixed to 1), \code{phi_<factor1>_<factor2>} (latent
#' correlation), the per-factor defined parameters
#' (\code{loading_sum_<factor>}, \code{error_sum_<factor>},
#' \code{omega_<factor>}, \code{ave_<factor>}, \code{H_<factor>}), and
#' the fit rows \code{chi_square}, \code{df}, \code{p_chi_square},
#' \code{cfi}, \code{tli}, \code{nnfi}, \code{rmsea},
#' \code{rmsea_ci_lower}, \code{rmsea_ci_upper}, \code{rmsea_ci_level},
#' \code{srmr}, \code{AIC}, \code{BIC}, \code{H0}, \code{H1}.
#'
#' @return
#' For \code{output = "verbose"} (default) and \code{output =
#' "measurement"}, a \code{data.frame} (classes \code{dmar_cfa_k},
#' \code{dmar_tbl}) with columns \code{syntax}, \code{term},
#' \code{estimate}, \code{se}, \code{z_value}, \code{p_value},
#' \code{ci_lower}, \code{ci_upper}. The \code{"model"} attribute is a
#' named character vector giving, per factor, the implied classical
#' structure; the printed header displays it. For \code{output =
#' "summary"}, the \pkg{lavaan} summary object; for
#' \code{"standardized"}, the standardized solution; for \code{"fit"},
#' the \pkg{lavaan} fit object.
#'
#' @references
#' Browne, M. W. (1974). Generalized least squares estimators in the
#'   analysis of covariance structures. \emph{South African Statistical
#'   Journal, 8}, 1--24.
#'
#' Browne, M. W. (1984). Asymptotically distribution-free methods for the
#'   analysis of covariance structures. \emph{British Journal of
#'   Mathematical and Statistical Psychology, 37}, 62--83.
#'
#' Muthén, B. (1984). A general structural equation model with
#'   dichotomous, ordered categorical, and continuous latent variable
#'   indicators. \emph{Psychometrika, 49}(1), 115--132.
#'
#' Muthén, B., du Toit, S. H. C., & Spisic, D. (1997). \emph{Robust
#'   inference using weighted least squares and quadratic estimating
#'   equations in latent variable modeling with categorical and continuous
#'   outcomes}. Unpublished technical report.
#'
#' Satorra, A., & Bentler, P. M. (1994). Corrections to test statistics
#'   and standard errors in covariance structure analysis. In A. von Eye &
#'   C. C. Clogg (Eds.), \emph{Latent variables analysis: Applications for
#'   developmental research} (pp. 399--419). Thousand Oaks, CA: Sage.
#'
#' Fornell, C., & Larcker, D. F. (1981). Evaluating structural equation
#'   models with unobservable variables and measurement error.
#'   \emph{Journal of Marketing Research, 18}(1), 39--50.
#'
#' Graham, J. M. (2006). Congeneric and (essentially) tau-equivalent
#'   estimates of score reliability: What they are and how to use them.
#'   \emph{Educational and Psychological Measurement, 66}(6), 930--944.
#'   \doi{10.1177/0013164406288165}
#'
#' Green, S. B., & Yang, Y. (2009). Reliability of summed item scores
#'   using structural equation modeling: An alternative to coefficient
#'   alpha. \emph{Psychometrika, 74}(1), 155--167.
#'   \doi{10.1007/s11336-008-9099-3}
#'
#' Hancock, G. R., & Mueller, R. O. (2001). Rethinking construct
#'   reliability within latent variable systems. In R. Cudeck, S. du
#'   Toit, & D. Sörbom (Eds.), \emph{Structural equation modeling:
#'   Present and future} (pp. 195--216). Scientific Software
#'   International.
#'
#' Henseler, J., Ringle, C. M., & Sarstedt, M. (2015). A new criterion
#'   for assessing discriminant validity in variance-based structural
#'   equation modeling. \emph{Journal of the Academy of Marketing
#'   Science, 43}(1), 115--135. \doi{10.1007/s11747-014-0403-8}
#'
#' Kelley, K., & Pornprasertmanit, S. (2016). Confidence intervals for
#'   population reliability coefficients: Evaluation of methods,
#'   recommendations, and software for composite measures.
#'   \emph{Psychological Methods, 21}, 69--92. \doi{10.1037/a0040086}
#'
#' Lord, F. M., & Novick, M. R. (1968). \emph{Statistical theories of
#'   mental test scores}. Addison-Wesley.
#'
#' McDonald, R. P. (1999). \emph{Test theory: A unified treatment}.
#'   Erlbaum.
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso
#' \code{\link{cfa_1}} and \code{\link{cfa_2}} for the one and two
#' factor convenience wrappers; \code{\link{plot_cfa_k}} to display the estimates and
#' the equality question visually; \code{\link{reliability_omega}},
#' \code{\link{reliability_H}}, \code{\link{average_variance_extracted}},
#' and \code{\link{htmt}} for the measurement properties as standalone
#' functions; \code{\link{measurement_invariance}} for the across-group
#' analog of these within-factor constraints;
#' \code{\link[lavaan]{lavaan}}, \code{\link[lavaan]{lavTestLRT}}.
#'
#' @family multivariate and latent variable methods
#'
#' @keywords multivariate
#'
#' @examples
#' data(holzinger_swineford)
#' hs_factors <- list(
#'   verbal    = c("t6_paragraph_comprehension", "t7_sentence",
#'                 "t9_word_meaning"),
#'   deduction = c("t20_deduction", "t22_problem_reasoning",
#'                 "t23_series_completion"))
#'
#' # Congeneric measurement model for both factors (the default:
#' # nothing is constrained, and the header names the structure).
#' cfa_k(holzinger_swineford, hs_factors)
#'
#' # Equal loadings within every factor. Because only the covariance
#' # structure identifies this constraint, the implied structure is
#' # essentially tau-equivalent, and the header says so.
#' cfa_k(holzinger_swineford, hs_factors, equal_loading = TRUE)
#'
#' # Measurement properties: omega, ave, and H per factor (each with a
#' # delta method standard error and confidence interval), the latent
#' # correlations, and the htmt ratios.
#' cfa_k(holzinger_swineford, hs_factors, output = "measurement")
#'
#' # The rest of the descriptor menu is shown but not run here, since
#' # each call refits the model. Descriptors can differ by factor:
#' # cfa_k(holzinger_swineford, hs_factors,
#' #       equal_loading = c(verbal = TRUE, deduction = FALSE))
#' #
#' # Equal loadings and intercepts (tau-equivalent), then also equal
#' # error variances (parallel). The mean structure is added because
#' # equal_intercept asks about it:
#' # cfa_k(holzinger_swineford, hs_factors, equal_loading = TRUE,
#' #       equal_intercept = TRUE)
#' # cfa_k(holzinger_swineford, hs_factors, equal_loading = TRUE,
#' #       equal_intercept = TRUE, equal_error = TRUE)
#' #
#' # Ordered-categorical items: the model is fit by WLSMV to polychoric
#' # correlations, and each ordered factor's omega is reported on the
#' # categorical sum score metric (Green & Yang, 2009):
#' # set.seed(113)
#' # eta <- rnorm(200)
#' # lat <- sweep(matrix(rep(eta, 6), 200, 6), 2,
#' #              seq(0.5, 0.8, length.out = 6), `*`) +
#' #   matrix(rnorm(200 * 6), 200, 6) %*%
#' #   diag(sqrt(1 - seq(0.5, 0.8, length.out = 6)^2))
#' # likert <- as.data.frame(apply(lat, 2, function(x)
#' #   as.integer(cut(x, breaks = c(-Inf, -1, 0, 1, Inf)))))
#' # names(likert) <- paste0("item_", 1:6)
#' # cfa_k(likert,
#' #       list(scale_a = paste0("item_", 1:3),
#' #            scale_b = paste0("item_", 4:6)),
#' #       ordered = TRUE, output = "measurement")
#' #
#' # Does the equal-loadings description hold? Two fits that differ only
#' # in a descriptor are nested, so output = "fit" hands them straight to
#' # lavaan's likelihood ratio test:
#' # fit_free  <- cfa_k(holzinger_swineford, hs_factors, output = "fit")
#' # fit_equal <- cfa_k(holzinger_swineford, hs_factors,
#' #                    equal_loading = TRUE, output = "fit")
#' # lavaan::lavTestLRT(fit_free, fit_equal)
#'
#' @export
cfa_k <- function(data = NULL, factors, S = NULL, N = NULL, M = NULL,
                  equal_loading = FALSE,
                  equal_intercept = FALSE,
                  equal_error = FALSE,
                  correlated_factors = TRUE,
                  meanstructure = NULL,
                  estimator = "ML",
                  missing = "listwise",
                  ordered = NULL,
                  se = "standard",
                  conf_level = 0.95,
                  output = c("verbose", "measurement", "summary",
                             "standardized", "fit"),
                  ...) {
  output <- match.arg(output)
  if (!requireNamespace("lavaan", quietly = TRUE)) {
    stop("The package 'lavaan' is needed. Install it with ",
         "install.packages(\"lavaan\").", call. = FALSE)
  }

  # ---- validate the factor specification -------------------------------
  if (!is.list(factors) || length(factors) < 1L ||
      is.null(names(factors)) || any(names(factors) == "")) {
    stop("'factors' must be a named list: ",
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
    stop("Each item may load on exactly one factor (simple structure); ",
         "duplicated: ",
         paste(unique(items[duplicated(items)]), collapse = ", "), ".",
         call. = FALSE)
  }


  # ---- expand the descriptor flags to one value per factor -------------
  expand_flag <- function(flag, arg_name) {
    if (!is.logical(flag) || anyNA(flag)) {
      stop("'", arg_name, "' must be logical: a single value, or a named ",
           "vector with one element per factor.", call. = FALSE)
    }
    if (length(flag) == 1L && is.null(names(flag))) {
      return(stats::setNames(rep(flag, length(factor_names)), factor_names))
    }
    if (is.null(names(flag)) || !setequal(names(flag), factor_names)) {
      stop("A vector '", arg_name, "' must be named with exactly the ",
           "factor names.", call. = FALSE)
    }
    flag[factor_names]
  }
  eq_load <- expand_flag(equal_loading, "equal_loading")
  eq_int  <- expand_flag(equal_intercept, "equal_intercept")
  eq_err  <- expand_flag(equal_error, "equal_error")

  # A single congeneric factor needs three items; with equal loadings
  # the two-item factor is just identified (the tau-equivalent case the
  # reliability family fits for a two-item scale).
  if (length(factors) == 1L && length(factors[[1L]]) < 3L &&
      !(length(factors[[1L]]) == 2L && isTRUE(unname(eq_load[1L])))) {
    stop("A single factor needs three or more items for identification ",
         "(two items are identified only with equal_loading = TRUE).",
         call. = FALSE)
  }

  # ---- whitelists, as in cfa_1() ---------------------------------------
  estimator_choices <- c("ML", "MLR", "WLS", "WLSMV", "GLS")
  missing_choices   <- c("listwise", "ml", "fiml", "pairwise")
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
  if (!is.numeric(conf_level) || length(conf_level) != 1L ||
      is.na(conf_level) || conf_level <= 0 || conf_level >= 1) {
    stop("'conf_level' must be a single number in (0, 1).", call. = FALSE)
  }

  # ---- input: raw data ('data') or covariance matrix ('S') -------------
  if (is.null(data) == is.null(S)) {
    stop("Supply exactly one of 'data' (raw data) or 'S' (a symmetric ",
         "covariance matrix, with 'N').", call. = FALSE)
  }
  input_is_raw_data <- !is.null(data)

  if (input_is_raw_data) {
    if (!is.data.frame(data) && !is.matrix(data)) {
      stop("'data' must be a data frame or matrix of raw data (rows are ",
           "respondents, columns include the items).", call. = FALSE)
    }
    if (is.matrix(data) && nrow(data) == ncol(data) &&
        isSymmetric(unname(data), tol = 1e-05)) {
      stop("'data' looks like a covariance matrix (square and symmetric); ",
           "pass it through 'S', with 'N'.", call. = FALSE)
    }
    raw_data <- as.data.frame(data)
    absent <- setdiff(items, names(raw_data))
    if (length(absent) > 0L) {
      stop("Items not found in the data: ",
           paste(absent, collapse = ", "), ".", call. = FALSE)
    }
    lavaan_args <- list(data = raw_data)
  } else {
    if (!is.matrix(S) || !isSymmetric(unname(S), tol = 1e-05)) {
      stop("'S' must be a symmetric covariance matrix; raw data are ",
           "passed through 'data'.", call. = FALSE)
    }
    if (is.null(N)) {
      stop("Sample size 'N' must be provided when input is a covariance ",
           "matrix.", call. = FALSE)
    }
    if (is.null(colnames(S))) {
      stop("A covariance matrix 'S' must have dimnames naming the items.",
           call. = FALSE)
    }
    rownames(S) <- colnames(S)
    absent <- setdiff(items, colnames(S))
    if (length(absent) > 0L) {
      stop("Items not found in the covariance matrix: ",
           paste(absent, collapse = ", "), ".", call. = FALSE)
    }
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
        stop("'ordered' names items that are not in 'factors': ",
             paste(absent, collapse = ", "), ".", call. = FALSE)
      }
      ordered_items <- unique(ordered)
    } else {
      stop("'ordered' must be TRUE (all items) or a character vector of ",
           "item names.", call. = FALSE)
    }
    if (!input_is_raw_data) {
      stop("Ordered items require raw data: thresholds and polychoric ",
           "correlations cannot be estimated from a covariance matrix.",
           call. = FALSE)
    }
    # A factor must be all ordered or all continuous: the sum score of a
    # mixed factor has no single metric for its omega.
    for (f in factor_names) {
      n_ord <- sum(factors[[f]] %in% ordered_items)
      if (n_ord > 0L && n_ord < length(factors[[f]])) {
        stop("Factor '", f, "' mixes ordered and continuous items; ",
             "declare all of a factor's items ordered, or none.",
             call. = FALSE)
      }
    }
    if (estimator %in% c("ML", "MLR", "GLS")) {
      message("Ordered items declared: switching to estimator = \"WLSMV\" ",
              "with robust standard errors.")
      estimator <- "WLSMV"
    }
    # Robust standard errors are the appropriate companion to WLSMV.
    if (se == "standard") se <- "robust.sem"
    if (missing %in% c("ml", "fiml")) {
      stop("Full information maximum likelihood is not available with ",
           "ordered items; use missing = \"pairwise\" or \"listwise\".",
           call. = FALSE)
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
    # restores exact agreement with cfa_1() under FIML.
    has_means <- any(eq_int) || !is.null(M) ||
      (input_is_raw_data && missing %in% c("ml", "fiml"))
  } else {
    if (!is.logical(meanstructure) || length(meanstructure) != 1L ||
        is.na(meanstructure)) {
      stop("'meanstructure' must be TRUE, FALSE, or NULL.", call. = FALSE)
    }
    has_means <- meanstructure
  }
  if (any(eq_int) && !has_means) {
    stop("'equal_intercept' asks about the mean structure; it cannot be ",
         "combined with meanstructure = FALSE.", call. = FALSE)
  }
  if (!has_means && input_is_raw_data && missing %in% c("ml", "fiml")) {
    stop("missing = \"", missing, "\" estimates the mean structure with ",
         "free item intercepts; it cannot be combined with ",
         "meanstructure = FALSE. Drop meanstructure = FALSE, or use ",
         "missing = \"listwise\".", call. = FALSE)
  }
  if (!is.null(ordered_items) && has_means) {
    stop("Thresholds replace intercepts for ordered items; ",
         "'equal_intercept', 'M', and meanstructure = TRUE are not ",
         "available with 'ordered'.", call. = FALSE)
  }
  if (has_means && !input_is_raw_data) {
    if (is.null(M)) {
      stop("The mean structure with a covariance matrix requires the item ",
           "means 'M'.", call. = FALSE)
    }
    if (!is.numeric(M) || is.null(names(M)) || !all(items %in% names(M))) {
      stop("'M' must be a named numeric vector covering every item.",
           call. = FALSE)
    }
    lavaan_args$sample.mean <- M[items]
  }

  # ---- build the lavaan model syntax -----------------------------------
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
  model_syntax <- paste(lines, collapse = "\n")

  # ---- fit -------------------------------------------------------------
  fit_args <- c(list(model = model_syntax, estimator = estimator, se = se,
                     meanstructure = has_means),
                lavaan_args, list(...))
  if (input_is_raw_data) fit_args$missing <- missing
  if (!is.null(ordered_items)) {
    fit_args$ordered <- ordered_items
    # The theta parameterization keeps the residual variances of the
    # ordered items as free parameters, so the labeled psi terms and
    # every defined measurement quantity built from them remain valid;
    # under the delta parameterization those variances are not free
    # and lavaan rejects the defined parameters.
    fit_args$parameterization <- "theta"
  }
  fit <- try(do.call(lavaan::lavaan, fit_args), silent = TRUE)
  if (inherits(fit, "try-error")) {
    stop("Model fitting with lavaan failed: ",
         conditionMessage(attr(fit, "condition")), call. = FALSE)
  }
  if (!lavaan::lavInspect(fit, "converged")) {
    stop("The model did not converge.", call. = FALSE)
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
    warning(structure(
      class = c("dmar_heywood_warning", "warning", "condition"),
      list(message = paste0(
        "Negative error variance estimate (Heywood case) for: ",
        paste(unique(pe_hey$lhs[hey]), collapse = ", "),
        ". Interpret the solution with care."),
        call = NULL)))
  }

  if (output == "fit") return(fit)
  if (output == "summary") {
    return(lavaan::summary(fit, fit.measures = TRUE))
  }
  if (output == "standardized") {
    return(lavaan::standardizedSolution(fit, level = conf_level))
  }

  # ---- name the implied classical structure, per factor ----------------
  model_name <- vapply(factor_names, function(f) {
    .cfa_k_model_name(eq_load[f], if (has_means) eq_int[f] else FALSE,
                      eq_err[f], has_means)
  }, character(1L))

  # ---- parameter table with confidence intervals -----------------------
  pe <- as.data.frame(lavaan::parameterEstimates(fit, level = conf_level))
  for (.col in c("se", "z", "pvalue", "ci.lower", "ci.upper")) {
    if (is.null(pe[[.col]])) pe[[.col]] <- NA_real_
  }
  pe$syntax <- paste(pe$lhs, pe$op, pe$rhs)
  defined_rows <- pe$op == ":="
  pe$label[defined_rows] <- pe$lhs[defined_rows]
  # The full := expressions (ave, H) are long enough to swamp the table;
  # the term names them, so the syntax column stays blank for those rows.
  pe$syntax[defined_rows] <- ""
  fv <- pe$op == "~~" & pe$lhs == pe$rhs & pe$lhs %in% factor_names
  pe$label[fv] <- paste0("phi_", pe$lhs[fv])
  # The factor means are fixed to 0 for identification and carry no
  # information; drop those rows.
  pe <- pe[!(pe$op == "~1" & pe$lhs %in% factor_names), , drop = FALSE]

  params <- data.frame(
    syntax   = pe$syntax,
    term     = pe$label,
    estimate = pe$est,
    se       = pe$se,
    z_value  = pe$z,
    p_value  = pe$pvalue,
    ci_lower = pe$ci.lower,
    ci_upper = pe$ci.upper,
    stringsAsFactors = FALSE
  )

  # ---- categorical sum score omega for ordered factors -----------------
  # The defined-parameter omega is on the latent response metric, which
  # is not the metric of a sum of observed categories. For a factor
  # whose items are ordered, replace it with the Green and Yang (2009)
  # categorical omega computed from the same fit; its delta method
  # interval does not transfer, so the interval columns are NA.
  if (any(ordered_factor)) {
    for (f in factor_names[ordered_factor]) {
      gy <- .omega_categorical_from_kfit(fit, factors[[f]])
      row_f <- params$term == paste0("omega_", f)
      params$estimate[row_f] <- gy
      params$se[row_f]       <- NA_real_
      params$z_value[row_f]  <- NA_real_
      params$p_value[row_f]  <- NA_real_
      params$ci_lower[row_f] <- NA_real_
      params$ci_upper[row_f] <- NA_real_
    }
    message("omega_",
            paste(factor_names[ordered_factor], collapse = ", omega_"),
            " reported on the categorical sum score metric (Green & ",
            "Yang, 2009).")
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
      ht <- htmt(raw_data, blocks = factors)
      out <- rbind(out, data.frame(
        syntax   = paste(ht$construct_1, "~~", ht$construct_2),
        term     = paste0("htmt_", ht$construct_1, "_", ht$construct_2),
        estimate = ht$htmt,
        se = NA_real_, z_value = NA_real_, p_value = NA_real_,
        ci_lower = NA_real_, ci_upper = NA_real_,
        stringsAsFactors = FALSE))
    }
    rownames(out) <- NULL
    class(out) <- c("dmar_cfa_k", "data.frame")
    out <- .as_dmar_tbl(out, conf_level = conf_level)
    attr(out, "model") <- model_name
    attr(out, "omega_metric") <- omega_metric
    return(out)
  }

  # ---- fit information rows --------------------------------------------
  fm <- lavaan::fitMeasures(fit)
  pick <- function(...) {
    for (nm in c(...)) {
      if (nm %in% names(fm) && is.finite(fm[nm])) return(unname(fm[nm]))
    }
    NA_real_
  }
  fit_rows <- data.frame(
    syntax = "",
    term = c("chi_square", "df", "p_chi_square", "cfi", "tli", "nnfi",
             "rmsea", "rmsea_ci_lower", "rmsea_ci_upper", "rmsea_ci_level",
             "srmr", "AIC", "BIC", "H0", "H1"),
    estimate = c(
      pick("chisq.scaled", "chisq"),
      pick("df.scaled", "df"),
      pick("pvalue.scaled", "pvalue"),
      pick("cfi.robust", "cfi.scaled", "cfi"),
      pick("tli.robust", "tli.scaled", "tli"),
      pick("nnfi.robust", "nnfi.scaled", "nnfi"),
      pick("rmsea.robust", "rmsea.scaled", "rmsea"),
      pick("rmsea.ci.lower.robust", "rmsea.ci.lower.scaled",
           "rmsea.ci.lower"),
      pick("rmsea.ci.upper.robust", "rmsea.ci.upper.scaled",
           "rmsea.ci.upper"),
      { lv <- pick("rmsea.ci.level"); if (is.na(lv)) 0.90 else lv },
      pick("srmr"),
      pick("aic"), pick("bic"),
      pick("logl"), pick("unrestricted.logl")),
    se = NA_real_, z_value = NA_real_, p_value = NA_real_,
    ci_lower = NA_real_, ci_upper = NA_real_,
    stringsAsFactors = FALSE
  )

  out <- rbind(params, fit_rows)
  rownames(out) <- NULL
  class(out) <- c("dmar_cfa_k", "data.frame")
  out <- .as_dmar_tbl(out, conf_level = conf_level,
                      p_terms = "p_chi_square",
                      fixed_terms = c("AIC", "BIC", "H0", "H1"))
  attr(out, "model") <- model_name
  attr(out, "factors") <- factors
  attr(out, "has_means") <- has_means
  attr(out, "omega_metric") <- omega_metric
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

#' @export
print.dmar_cfa_k <- function(x, ...) {
  model <- attr(x, "model")
  if (!is.null(model)) {
    cat("Measurement structure, per factor:\n")
    for (f in names(model)) {
      cat("  ", f, ": ", model[[f]], "\n", sep = "")
    }
    cat("\n")
  }
  NextMethod()
  invisible(x)
}

#' Tidy a Multiple-Factor CFA Fit
#'
#' Returns the parameter rows of a \code{\link{cfa_k}} table (loadings,
#' error variances, intercepts, latent correlations, and the defined
#' measurement properties) in the column convention used by the
#' \pkg{broom} ecosystem. Fit-information rows belong in
#' \code{\link{glance.dmar_cfa_k}}.
#'
#' @param x A \code{dmar_cfa_k} object returned by \code{\link{cfa_k}}.
#' @param \dots Unused.
#'
#' @return A \code{data.frame} with columns \code{term},
#'   \code{estimate}, \code{se}, \code{statistic},
#'   \code{p_value}, \code{ci_lower}, \code{ci_upper}.
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @examples
#' data(holzinger_swineford)
#' res <- cfa_k(holzinger_swineford,
#'              list(verbal = c("t6_paragraph_comprehension",
#'                              "t7_sentence", "t9_word_meaning"),
#'                   deduction = c("t20_deduction",
#'                                 "t22_problem_reasoning",
#'                                 "t23_series_completion")))
#' generics::tidy(res)
#' generics::glance(res)
#'
#' @importFrom generics tidy
#' @exportS3Method generics::tidy
tidy.dmar_cfa_k <- function(x, ...) {
  fit_terms <- c("chi_square", "df", "p_chi_square", "cfi", "tli", "nnfi",
                 "rmsea", "rmsea_ci_lower", "rmsea_ci_upper",
                 "rmsea_ci_level", "srmr", "AIC", "BIC", "H0", "H1")
  rows <- x[!x$term %in% fit_terms, , drop = FALSE]
  data.frame(
    term      = rows$term,
    estimate  = rows$estimate,
    se = rows$se,
    statistic = rows$z_value,
    p_value   = rows$p_value,
    ci_lower  = rows$ci_lower,
    ci_upper = rows$ci_upper,
    stringsAsFactors = FALSE,
    row.names = NULL
  )
}

#' Glance at a Multiple-Factor CFA Fit
#'
#' Returns a one-row \code{data.frame} of model-level summaries from a
#' \code{\link{cfa_k}} table, in the column convention used by the
#' \pkg{broom} ecosystem.
#'
#' @param x A \code{dmar_cfa_k} object returned by \code{\link{cfa_k}}.
#' @param \dots Unused.
#'
#' @return A one-row \code{data.frame} with columns \code{chi_square},
#'   \code{df}, \code{p_value}, \code{cfi}, \code{tli}, \code{rmsea},
#'   \code{rmsea_low}, \code{rmsea_high}, \code{srmr}, \code{AIC},
#'   \code{BIC}, \code{logLik}.
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @importFrom generics glance
#' @exportS3Method generics::glance
glance.dmar_cfa_k <- function(x, ...) {
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
    stringsAsFactors = FALSE
  )
}
