#' Graded Response Model for Ordered Categorical Items
#'
#' @description
#' Estimates Samejima's (1969) graded response model for a set of ordered
#' categorical items (Likert items, symptom severity ratings, rubric
#' scored performance items) and reports each item's discrimination and
#' its category boundary locations. The model is fitted as the single-factor
#' categorical factor analysis model it provably is (Takane and de Leeuw,
#' 1987), using \pkg{lavaan}'s categorical estimator on the polychoric
#' correlations, and the solution is then converted to the normal ogive
#' (or logistic) item response theory parameterization. A researcher who
#' already fits confirmatory factor analysis models therefore gets item
#' response theory item parameters without adopting a second estimation
#' engine, and the two analyses of the same items stay in one modeling
#' tradition.
#'
#' @param data A \code{data.frame} or matrix of ordered item responses,
#'   one column per item, coded with integer category values (for
#'   example 1, 2, 3, 4, 5). Rows with any missing value on the analyzed
#'   items are listwise-deleted. Each item must have at least 2 and at
#'   most 20 distinct observed categories; a column with more than 20
#'   distinct values, or with non-integer values, is treated as
#'   continuous and rejected.
#' @param items Optional character vector naming the columns of
#'   \code{data} to analyze. Defaults to \code{NULL}, which uses every
#'   column. At least 3 items are required.
#' @param estimator Character; the \pkg{lavaan} estimator for
#'   categorical data. The choices differ in the weight matrix applied
#'   to the polychoric correlations and in whether the test statistic is
#'   corrected. \code{"WLSMV"} (the default) is diagonally weighted
#'   least squares with a mean- and variance-adjusted test statistic
#'   (Muthén, 1984; Muthén, du Toit, & Spisic, 1997), the standard
#'   estimator for ordered categorical items; \code{"WLSM"} applies the
#'   mean adjustment only. \code{"DWLS"} is the same diagonal-weight
#'   estimator with no correction. \code{"WLS"} uses the full weight
#'   matrix (the asymptotic distribution free approach; Browne, 1984),
#'   which is unstable unless the sample is large relative to the number
#'   of thresholds. \code{"ULS"}, unweighted least squares, uses an
#'   identity weight matrix, an option worth considering in small
#'   samples where even the diagonal weights are noisy;
#'   \code{"ULSMV"} and \code{"ULSM"} add the corrected test
#'   statistics. When in doubt keep the default.
#' @param metric Which discrimination metric to report in the \code{a}
#'   column: \code{"normal_ogive"} (default) or \code{"logistic"}. Both
#'   are always computed and the one not reported in \code{a} is
#'   attached as an attribute, so the returned table has the same
#'   columns and the same number of rows either way. The boundary
#'   locations \code{b} are identical under the two metrics.
#'
#' @details
#' \strong{One model, two parameterizations.} Samejima's (1969) graded
#' response model and the single-factor categorical factor analysis model
#' of Muthén (1984) are the same model written in different
#' parameterizations; Takane and de Leeuw (1987) proved the equivalence,
#' and Kamata and Bauer (2008) give the algebra item by item. Each
#' observed response \eqn{X_i} is a categorization of a latent continuous
#' response variate \eqn{X_i^{*}} at thresholds \eqn{\tau_{ik}}, and
#' \eqn{X_i^{*} = \lambda_i \theta + \varepsilon_i} with \eqn{\theta}
#' standard normal and \eqn{X_i^{*}} standardized. Fitting that model on
#' the polychoric correlations and converting the solution gives the
#' normal ogive graded response model directly. For item \eqn{i} with
#' standardized loading \eqn{\lambda_i} and standardized thresholds
#' \eqn{\tau_{ik}},
#' \deqn{a_i = \frac{\lambda_i}{\sqrt{1 - \lambda_i^2}}, \qquad
#'       b_{ik} = \frac{\tau_{ik}}{\lambda_i}.}
#' The boundary response function is
#' \deqn{P^{*}_{ik}(\theta) = \Phi\!\left[a_i (\theta - b_{ik})\right],}
#' the probability of responding above boundary \eqn{k}, with
#' \eqn{P^{*}_{i0}(\theta) \equiv 1} and \eqn{P^{*}_{iK}(\theta) \equiv
#' 0}; the probability of the individual category is the difference of
#' adjacent boundary functions,
#' \eqn{P_{ik}(\theta) = P^{*}_{i,k-1}(\theta) - P^{*}_{ik}(\theta)}.
#' When \eqn{a_i > 0}, that is, for an item keyed in the same direction as
#' the rest of the scale, \eqn{P^{*}_{ik}} is monotone increasing in
#' \eqn{\theta}, the boundary locations of the item are ordered,
#' \eqn{b_{i1} < b_{i2} < \cdots}, and \eqn{b_{ik}} is the value of
#' \eqn{\theta} at which the probability of responding above boundary
#' \eqn{k} reaches 0.50. An item keyed in the opposite direction has
#' \eqn{\lambda_i < 0}, hence \eqn{a_i < 0} and boundary locations that run
#' from high to low; see the two paragraphs on direction below.
#'
#' \strong{The direction of the latent variable.} A single-factor model
#' fixes \eqn{\theta} only up to its direction. Relabeling \eqn{\theta} as
#' \eqn{-\theta} changes the sign of every loading and leaves the fitted
#' model, the thresholds, and every fit measure exactly as they were, so it
#' is a renaming of the latent direction rather than a different model. The
#' thresholds are untouched because \eqn{\tau_{ik}} cuts the item's own
#' latent response variate \eqn{X_i^{*}}, which the relabeling does not
#' move; the sign change therefore passes straight through to
#' \eqn{a_i = \lambda_i / \sqrt{1 - \lambda_i^2}} and to
#' \eqn{b_{ik} = \tau_{ik} / \lambda_i}, both of which change sign.
#' \pkg{lavaan} returns whichever direction its starting values point
#' toward, and for a scale that contains a reverse-keyed item that
#' direction can turn on something as incidental as the order of the
#' columns. The solution is therefore put in a fixed direction before it is
#' converted: if the standardized loadings sum to a negative number the
#' whole factor is flipped, so that \eqn{\theta} runs in the direction the
#' scale as a whole measures. The result is the same table no matter how
#' the columns are ordered. Whether the flip was applied is recorded on the
#' \code{"factor_sign_flipped"} attribute. The \pkg{lavaan} object on the
#' \code{"fit"} attribute is the fit as \pkg{lavaan} produced it, so when a
#' flip was applied its loadings carry the opposite sign to the
#' \code{lambda} column.
#'
#' \strong{Reverse-keyed items.} An item whose loading is still negative
#' after the direction is fixed is keyed opposite to the rest of the scale,
#' which is a property of the item rather than an artifact of the sign
#' indeterminacy. Its discrimination is negative and its boundary locations
#' run from high to low, so it does not satisfy the graded response model
#' as written above and its parameters do not belong on the same scale as
#' the others. Such items are named in a warning. Reverse score them (for
#' example \code{x <- (min(x) + max(x)) - x}) and refit; that puts the item
#' in the direction the rest of the scale measures and restores
#' \eqn{a_i > 0} and the ordering \eqn{b_{i1} < b_{i2} < \cdots}.
#'
#' \strong{The two discrimination metrics and the constant 1.702.} The
#' conversion above puts \eqn{a_i} in the normal ogive metric, where the
#' boundary function is a normal cumulative distribution function. The
#' item response theory literature more often writes the graded response
#' model with a logistic boundary function, and the two agree closely
#' once the logistic argument is stretched by a scaling constant:
#' \eqn{|\Phi(x) - \Psi(1.702 x)| < 0.01} for every \eqn{x}, where
#' \eqn{\Psi} is the standard logistic cumulative distribution function.
#' The value 1.702 is the constant that minimizes that maximum
#' discrepancy (Haley, 1952; see Camilli, 1994, for the history), so
#' \eqn{a_i(\mathrm{logistic}) = 1.702 \, a_i(\mathrm{normal\ ogive})}
#' and software that reports logistic slopes (for example \pkg{mirt} and
#' the classical two parameter logistic tradition) gives values about 1.7
#' times larger for the same items. The scaling multiplies the slope and
#' leaves the location alone, so \eqn{b_{ik}} does not depend on the
#' metric.
#'
#' \strong{Estimation and what to expect.} \pkg{lavaan} estimates the
#' thresholds and the polychoric correlations, then fits the single-factor
#' model to those correlations by (diagonally) weighted least squares.
#' This is limited information estimation: it uses the univariate and
#' bivariate margins of the response table, whereas marginal maximum
#' likelihood (the usual item response theory approach, as in \pkg{mirt})
#' uses the full response pattern likelihood. The two are consistent for
#' the same population parameters and agree closely in practice, but they
#' are different estimators and will not return identical numbers on a
#' finite sample. Limited information estimation scales well to many
#' items and brings the whole apparatus of factor analysis fit assessment
#' (CFI, TLI, RMSEA) along with it; the fit measures are returned on the
#' \code{"fit_measures"} attribute and the \pkg{lavaan} object itself on
#' \code{"fit"}, so any \pkg{lavaan} accessor can be applied to the
#' result.
#'
#' The model is unidimensional by construction. A standardized loading at
#' or beyond one is an improper (Heywood) solution: the implied
#' discrimination is infinite and the conversion is not interpretable.
#' That case is flagged with a warning rather than silently returned as a
#' number.
#'
#' This function requires \pkg{lavaan} to be installed.
#'
#' @return A \code{data.frame} (class \code{dmar_tbl}) with one row per
#'   item and category boundary and the columns
#'   \describe{
#'     \item{\code{item}}{Item name, taken from the column name.}
#'     \item{\code{factor}}{Name of the latent variable, the same for
#'       every row in this unidimensional model.}
#'     \item{\code{category}}{Boundary index \eqn{k}, running from 1 to
#'       one fewer than the item's number of categories.}
#'     \item{\code{lambda}}{Standardized factor loading \eqn{\lambda_i}
#'       of the item's latent response variate on the factor, repeated
#'       across the item's boundaries.}
#'     \item{\code{tau}}{Standardized threshold \eqn{\tau_{ik}}.}
#'     \item{\code{a}}{Discrimination in the metric named by
#'       \code{metric}, repeated across the item's boundaries.}
#'     \item{\code{b}}{Boundary location \eqn{b_{ik}} on the \eqn{\theta}
#'       scale.}
#'   }
#'   The attributes are \code{"fit"} (the fitted \pkg{lavaan} object),
#'   \code{"fit_measures"} (the full named numeric vector from
#'   \code{lavaan::fitMeasures}, unrounded), \code{"metric"} (the
#'   reported discrimination metric), \code{"estimator"},
#'   \code{"n_categories"} (named integer vector of the number of
#'   observed categories per item), \code{"N"} (the analyzed sample
#'   size), \code{"factor_sign_flipped"} (a single logical recording
#'   whether the direction of the latent variable was reversed to satisfy
#'   the sign convention described in Details), and whichever of
#'   \code{"a_logistic"} or
#'   \code{"a_normal_ogive"} was not reported in the \code{a} column (a
#'   named numeric vector, one element per item).
#'
#' @references
#' Camilli, G. (1994). Teacher's corner: Origin of the scaling constant
#'   \emph{d} = 1.7 in item response theory. \emph{Journal of Educational
#'   and Behavioral Statistics, 19}(3), 293--295.
#'   \doi{10.3102/10769986019003293}
#'
#' Haley, D. C. (1952). \emph{Estimation of the dosage mortality
#'   relationship when the dose is subject to error} (Technical Report
#'   No. 15). Applied Mathematics and Statistics Laboratory, Stanford
#'   University.
#'
#' Kamata, A., & Bauer, D. J. (2008). A note on the relation between
#'   factor analytic and item response theory models. \emph{Structural
#'   Equation Modeling, 15}(1), 136--153. \doi{10.1080/10705510701758406}
#'
#' Muthén, B. (1984). A general structural equation model with
#'   dichotomous, ordered categorical, and continuous latent variable
#'   indicators. \emph{Psychometrika, 49}(1), 115--132.
#'
#' Samejima, F. (1969). Estimation of latent ability using a response
#'   pattern of graded scores. \emph{Psychometrika Monograph Supplement,
#'   34}(4, Pt. 2), 1--97.
#'
#' Takane, Y., & de Leeuw, J. (1987). On the relationship between item
#'   response theory and factor analysis of discretized variables.
#'   \emph{Psychometrika, 52}(3), 393--408.
#'
#' Wirth, R. J., & Edwards, M. C. (2007). Item factor analysis: Current
#'   approaches and future directions. \emph{Psychological Methods,
#'   12}(1), 58--79. \doi{10.1037/1082-989X.12.1.58}
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso \code{\link{cfa_1}} (the same single-factor model reported in
#'   the factor analysis parameterization),
#'   \code{\link{reliability_omega_categorical}} (reliability for the
#'   same class of items), \code{\link[lavaan]{cfa}}.
#'
#' @family multivariate and latent variable methods
#'
#' @keywords multivariate models
#'
#' @examples
#' # Six five-category items generated from a known graded response model.
#' set.seed(113)
#' n <- 800
#' a_pop <- c(1.2, 0.9, 1.5, 1.0, 1.3, 1.1)
#' b_pop <- rbind(c(-1.6, -0.6, 0.3, 1.2), c(-1.4, -0.4, 0.5, 1.5),
#'                c(-1.8, -0.7, 0.2, 1.1), c(-1.2, -0.2, 0.7, 1.6),
#'                c(-1.5, -0.5, 0.4, 1.3), c(-1.3, -0.3, 0.6, 1.4))
#' theta <- rnorm(n)
#' responses <- vapply(seq_along(a_pop), function(i) {
#'   p_star <- outer(theta, b_pop[i, ], function(z, b) pnorm(a_pop[i] * (z - b)))
#'   as.integer(1 + rowSums(runif(n) < p_star))
#' }, integer(n))
#' colnames(responses) <- paste0("item", seq_along(a_pop))
#' responses <- as.data.frame(responses)
#'
#' # Item parameters in the normal ogive metric.
#' grm <- irt_grm(responses)
#' grm
#'
#' # The generating discriminations, for comparison.
#' a_pop
#'
#' # Model fit travels with the item parameters.
#' attr(grm, "fit_measures")[c("cfi", "tli", "rmsea", "srmr")]
#'
#' # Two further calls, each of which refits the model and so is not run
#' # here. The first reports the same fit with logistic slopes (about
#' # 1.702 times larger than the normal ogive slopes above), the second
#' # fits a subset of the items selected by name:
#' # irt_grm(responses, metric = "logistic")
#' # irt_grm(responses, items = c("item1", "item3", "item5"))
#'
#' # The boundary response function of the first item at theta = 0.
#' first <- grm[grm$item == "item1", ]
#' pnorm(first$a * (0 - first$b))
#'
#' @export
irt_grm <- function(data, items = NULL, estimator = "WLSMV",
                    metric = c("normal_ogive", "logistic")) {
  metric <- match.arg(metric)

  if (!requireNamespace("lavaan", quietly = TRUE)) {
    stop("Package 'lavaan' is required by irt_grm(). Install it with ",
         "install.packages(\"lavaan\").", call. = FALSE)
  }

  estimator_choices <- c("WLSMV", "WLSM", "DWLS", "WLS",
                         "ULSMV", "ULSM", "ULS")
  if (!is.character(estimator) || length(estimator) != 1L ||
      !(estimator %in% estimator_choices)) {
    stop("'estimator' must be one of ",
         paste(dQuote(estimator_choices), collapse = ", "),
         "; these are the lavaan estimators for categorical data.",
         call. = FALSE)
  }

  if (missing(data) || is.null(data)) {
    stop("'data' (a data frame or matrix of ordered item responses) is ",
         "required.", call. = FALSE)
  }
  if (!is.data.frame(data) && !is.matrix(data)) {
    stop("'data' must be a data.frame or a matrix of ordered item ",
         "responses.", call. = FALSE)
  }
  data <- as.data.frame(data)
  if (is.null(colnames(data)) || any(!nzchar(colnames(data)))) {
    stop("Every column of 'data' must be named; the names become the item ",
         "labels.", call. = FALSE)
  }

  if (!is.null(items)) {
    if (!is.character(items)) {
      stop("'items' must be a character vector naming columns of 'data'.",
           call. = FALSE)
    }
    missing_items <- setdiff(items, colnames(data))
    if (length(missing_items)) {
      stop("These 'items' are not columns of 'data': ",
           paste(missing_items, collapse = ", "), ".", call. = FALSE)
    }
    if (anyDuplicated(items)) {
      stop("'items' must not repeat an item.", call. = FALSE)
    }
    data <- data[, items, drop = FALSE]
  }

  item_names <- colnames(data)
  if (anyDuplicated(item_names)) {
    stop("The item names must be unique.", call. = FALSE)
  }
  # lavaan model syntax names the variables directly, so a non-syntactic
  # column name (a space, a leading digit) cannot be parsed. Say so rather
  # than letting lavaan fail on generated syntax the user never wrote.
  if (!identical(make.names(item_names), item_names)) {
    bad <- item_names[make.names(item_names) != item_names]
    stop("Item names must be syntactically valid R names because they are ",
         "used in the lavaan model syntax. Rename: ",
         paste(bad, collapse = ", "), ".", call. = FALSE)
  }

  data <- data[stats::complete.cases(data), , drop = FALSE]
  N <- nrow(data)
  n_items <- ncol(data)
  if (n_items < 3L) {
    stop("At least 3 items are required to identify a single-factor graded ",
         "response model; ", n_items,
         if (n_items == 1L) " was supplied." else " were supplied.",
         call. = FALSE)
  }
  if (N < 2L) {
    stop("Fewer than 2 complete cases remain after listwise deletion.",
         call. = FALSE)
  }

  # The classic user error is handing a graded response model a continuous
  # variable. Catch both faces of it: non-integer codes, and integer codes
  # with implausibly many distinct values.
  max_categories <- 20L
  n_categories <- integer(n_items)
  names(n_categories) <- item_names
  for (j in seq_len(n_items)) {
    x <- data[[j]]
    if (is.factor(x)) x <- as.integer(x)
    if (!is.numeric(x)) {
      stop("Item '", item_names[j], "' is not numeric. The graded response ",
           "model requires ordered categorical responses coded as integers.",
           call. = FALSE)
    }
    if (any(!is.finite(x))) {
      stop("Item '", item_names[j], "' contains non-finite values.",
           call. = FALSE)
    }
    if (!isTRUE(all.equal(x, round(x)))) {
      stop("Item '", item_names[j], "' contains non-integer values, so it ",
           "looks continuous. The graded response model is for ordered ",
           "categorical responses coded as integers (for example 1, 2, 3, ",
           "4, 5). Categorize the variable first, or model it with a ",
           "continuous factor analysis model such as cfa_1().",
           call. = FALSE)
    }
    k <- length(unique(x))
    if (k < 2L) {
      stop("Item '", item_names[j], "' has only ", k, " observed category. ",
           "Every item must have at least 2 categories; an item with no ",
           "variability carries no information about the latent variable.",
           call. = FALSE)
    }
    if (k > max_categories) {
      stop("Item '", item_names[j], "' has ", k, " distinct values, more ",
           "than the ", max_categories, " allowed, so it looks continuous ",
           "rather than ordered categorical. Collapse the categories, or ",
           "model the variable with a continuous factor analysis model ",
           "such as cfa_1().", call. = FALSE)
    }
    n_categories[j] <- k
    data[[j]] <- ordered(x)
  }

  # Name the latent variable so it cannot collide with an item name.
  factor_name <- "theta"
  while (factor_name %in% item_names) factor_name <- paste0(factor_name, "_")

  model <- paste0(factor_name, " =~ ", paste(item_names, collapse = " + "))
  fit <- try(lavaan::cfa(model, data = data, ordered = item_names,
                         estimator = estimator, std.lv = TRUE),
             silent = TRUE)
  if (inherits(fit, "try-error")) {
    stop("The single-factor categorical model could not be fitted by lavaan: ",
         conditionMessage(attr(fit, "condition")), call. = FALSE)
  }
  if (!isTRUE(lavaan::lavInspect(fit, "converged"))) {
    stop("The single-factor categorical model did not converge. Check for ",
         "sparse categories (collapsing rare categories often helps) or ",
         "items with near-zero association.", call. = FALSE)
  }

  # Standardized loadings and thresholds. standardizedSolution() is used
  # rather than the raw estimates so the conversion holds under either
  # lavaan parameterization: it divides by the standard deviation of the
  # latent response variate, which the conversion formulas assume is one.
  std <- lavaan::standardizedSolution(fit)
  load_rows <- std[std$op == "=~", , drop = FALSE]
  thr_rows <- std[std$op == "|", , drop = FALSE]
  lambda <- stats::setNames(load_rows$est.std[match(item_names, load_rows$rhs)],
                            item_names)
  if (anyNA(lambda)) {
    stop("Could not recover a standardized loading for every item from the ",
         "lavaan solution; this is an internal error.", call. = FALSE)
  }

  # The model is identified only up to the direction of the latent variable:
  # relabeling theta as -theta flips every loading, leaves every threshold
  # alone (a threshold cuts the item's own latent response variate, which the
  # relabeling does not move), and so flips every discrimination and every
  # boundary location. lavaan returns whichever direction its starting values
  # point toward, which for a scale containing a reverse-keyed item can turn
  # on the order of the columns. Fix the direction here so the reported item
  # parameters do not depend on it. See Details.
  lambda_sum <- sum(lambda)
  factor_sign_flipped <- if (lambda_sum < 0) {
    TRUE
  } else if (lambda_sum > 0) {
    FALSE
  } else {
    # Exact tie. Break it on the first nonzero loading so the convention
    # stays deterministic; if every loading is zero, leave the sign alone.
    first_nonzero <- which(lambda != 0)[1L]
    !is.na(first_nonzero) && lambda[[first_nonzero]] < 0
  }
  if (isTRUE(factor_sign_flipped)) lambda <- -lambda

  improper <- any(abs(lambda) >= 1)
  if (improper) {
    warning("Improper (Heywood) solution: a standardized loading is at or ",
            "beyond one, so the implied discrimination is not finite. The ",
            "item parameters are returned but should not be interpreted ",
            "until the model is respecified.", call. = FALSE)
  }

  # A loading that is still negative once the direction is fixed belongs to an
  # item keyed opposite to the rest of the scale. Its discrimination is
  # negative and its boundary locations descend, which the graded response
  # model does not allow, so say so rather than returning it silently.
  reverse_keyed <- item_names[lambda < 0]
  if (length(reverse_keyed)) {
    one <- length(reverse_keyed) == 1L
    warning(if (one) "Item " else "Items ",
            paste(sQuote(reverse_keyed), collapse = ", "),
            if (one) " loads" else " load",
            " negatively on the factor while the rest of the scale loads ",
            "positively, so ", if (one) "it is" else "they are",
            " keyed in the opposite direction. The reported discrimination ",
            "is negative and the boundary locations run from high to low, ",
            "which the graded response model does not allow. Reverse score ",
            if (one) "the item" else "these items",
            " (for example x <- (min(x) + max(x)) - x) and refit.",
            call. = FALSE)
  }

  a_normal_ogive <- lambda / sqrt(1 - lambda^2)
  # 1.702 rescales the normal ogive slope to the logistic metric; see
  # Details for why that constant and not another.
  scaling_constant <- 1.702
  a_logistic <- scaling_constant * a_normal_ogive
  a_reported <- if (metric == "logistic") a_logistic else a_normal_ogive

  rows <- vector("list", n_items)
  for (j in seq_len(n_items)) {
    nm <- item_names[j]
    tj <- thr_rows[thr_rows$lhs == nm, , drop = FALSE]
    # lavaan labels the thresholds t1, t2, ...; order by that index rather
    # than trusting the row order of the parameter table.
    k_index <- as.integer(sub("^t", "", tj$rhs))
    tj <- tj[order(k_index), , drop = FALSE]
    k_index <- sort(k_index)
    if (length(k_index) != n_categories[[j]] - 1L) {
      stop("lavaan returned ", length(k_index), " thresholds for item '", nm,
           "', which has ", n_categories[[j]],
           " categories; this is an internal error.", call. = FALSE)
    }
    tau <- tj$est.std
    rows[[j]] <- data.frame(
      item = nm,
      factor = factor_name,
      category = k_index,
      lambda = unname(lambda[[j]]),
      tau = tau,
      a = unname(a_reported[[j]]),
      b = tau / lambda[[j]],
      stringsAsFactors = FALSE
    )
  }
  out <- do.call(rbind, rows)
  row.names(out) <- NULL

  res <- .as_dmar_tbl(out)
  attr(res, "fit") <- fit
  attr(res, "fit_measures") <- lavaan::fitMeasures(fit)
  attr(res, "metric") <- metric
  attr(res, "estimator") <- estimator
  attr(res, "n_categories") <- n_categories
  attr(res, "N") <- N
  attr(res, "factor_sign_flipped") <- factor_sign_flipped
  if (metric == "logistic") {
    attr(res, "a_normal_ogive") <- a_normal_ogive
  } else {
    attr(res, "a_logistic") <- a_logistic
  }
  res
}
