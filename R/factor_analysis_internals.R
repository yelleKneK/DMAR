# Input and output shared by the exploratory and SAM interfaces. Keep raw
# observations distinct from covariance input, including declared SPSS NAs.
.fa_require_lavaan <- function() {
  if (!requireNamespace("lavaan", quietly = TRUE)) {
    .dmar_refuse(c(
      "This analysis requires lavaan. Install it with ",
      "install.packages(\"lavaan\")."), "package_required", package = "lavaan")
  }
}

.fa_level <- function(conf_level, name = "conf_level") {
  if (!is.numeric(conf_level) || length(conf_level) != 1L ||
      !is.finite(conf_level) || conf_level <= 0 || conf_level >= 1) {
    .dmar_refuse_argument(c(
      "'", name, "' must be one number strictly between 0 and 1."), arg = NULL)
  }
}

# Items that cannot be read from the input are refused with class
# dmar_items_invalid (a dmar_refusal; .dmar_refuse(), R/dmar_conditions.R),
# whichever argument named them: `items` of fa_exploratory(),
# factor_count(), and parallel_analysis() (here and in .fa_covariance()),
# the factors of fa_confirmatory() and the functions that call its engine
# (.fa_cfa_engine()), `items` of measurement_invariance() and
# measurement_alignment(), and the indicators of measurement_invariance()'s
# model (its list form and its lavaan syntax). sem_sam() names an
# indicator its data lack itself, before .fa_data(), which still refuses
# repeated column names in its data ("data_names"). A program then reads a
# misnamed item from the class rather than from the wording, which differs
# by path. The `reason` element names the check that failed and `items`
# the names at issue (character(0) when the failure is not about
# particular names); see ?dmar_metadata.
.fa_items_refuse <- function(message, reason, items = character(0)) {
  .dmar_refuse(message, "items_invalid", reason = reason, items = items)
}

# The first check that item names fail against the names of the input, in
# the order "form" (not a character vector, or a name is NA), "unknown"
# (a name the input does not have, the empty name included), "duplicate"
# (a name given twice), and "count" (fewer than three), as a list of the
# reason and the names at issue; NULL when the names pass every check.
.fa_items_check <- function(items, available) {
  if (!is.character(items) || anyNA(items)) {
    return(list(reason = "form", items = character(0)))
  }
  unknown <- unique(items[!items %in% available | !nzchar(items)])
  if (length(unknown)) return(list(reason = "unknown", items = unknown))
  if (anyDuplicated(items)) {
    return(list(reason = "duplicate", items = unique(items[duplicated(items)])))
  }
  if (length(items) < 3L) return(list(reason = "count", items = character(0)))
  NULL
}

# An item that a model names and the input lacks (the factors of
# .fa_cfa_engine()), refused with the reason .fa_items_check() gives the
# same name: "form" when a name is NA, "unknown" otherwise. The message is
# the caller's, which names the absent items as it did before the class.
.fa_items_absent <- function(message, absent) {
  if (anyNA(absent)) .fa_items_refuse(message, "form")
  .fa_items_refuse(message, "unknown", absent)
}

.fa_data <- function(data, items = NULL) {
  if (!is.data.frame(data) && !is.matrix(data)) {
    .dmar_refuse_argument(
      "'data' must be a data frame or numeric matrix.", arg = "data")
  }
  data <- .relia_declared_missing(as.data.frame(data))
  # Column names that cannot serve as item names (two alike, or, when
  # every column is an item, one empty or NA) are the data's fault
  # rather than the items', and are checked first ("data_names").
  columns <- names(data)
  every <- is.null(items)
  bad <- duplicated(columns) | duplicated(columns, fromLast = TRUE) |
    (every & (is.na(columns) | !nzchar(columns)))
  if (every) items <- columns
  failed <- if (any(bad)) {
    list(reason = "data_names", items = unique(columns[bad]))
  } else {
    .fa_items_check(items, columns)
  }
  if (!is.null(failed)) {
    .fa_items_refuse(c("'items' must name at least three distinct columns of data, ",
                       "whose column names must be unique."),
                     failed$reason, failed$items)
  }
  data <- data[items]
  if (any(!vapply(data, is.numeric, logical(1L)))) {
    .dmar_refuse_argument(c(
      "Every analysis item must be numeric. Check the declared item types ",
      "and use an analysis that supports them."), arg = NULL)
  }
  if (any(vapply(data, function(x) any(!is.finite(x) & !is.na(x)), logical(1L)))) {
    .dmar_refuse_argument(
      "Analysis items must contain finite values or NA.", arg = NULL)
  }
  constant <- items[vapply(data, function(x) length(unique(x[!is.na(x)])) < 2L,
                           logical(1L))]
  if (length(constant)) {
    .dmar_refuse_argument(c(
      "Items with fewer than two observed values: ",
      paste(constant, collapse = ", "), "."), arg = NULL)
  }
  data
}

# Covariance input, read as .fa_data() reads raw data: the items are
# selected first and only they are checked. The shape and the names are
# checked on the whole of S (.fa_covariance_names()), since the items are
# read from the names; `items` then reduces S to the named rows and
# columns, in the order given (NULL keeps every item); and the entries
# (finite, symmetric, positive variances, positive definite) are checked
# on the reduced matrix alone, before N. A problem confined to items not
# analyzed (a missing entry, or a pairwise estimate that leaves the whole
# matrix indefinite) therefore does not refuse the call, and S with items
# is analyzed exactly as S[items, items] is. Shared by fa_exploratory(),
# factor_count(), and parallel_analysis(), so that each reads the names
# the same way. Items that cannot be read from S are refused with class
# dmar_items_invalid and the reasons of .fa_data() (.fa_items_check();
# the names of S itself are checked first, by .fa_covariance_names(), so
# "data_names" has no counterpart here); a malformed S or N stops without
# a class. The message names any non-NA item that S lacks, whatever the
# reason, so with an NA among the items (reason "form" and `items` empty,
# as on the raw data path) its text can name an absent item that the
# condition's `items` does not carry.
.fa_covariance <- function(S, N, items = NULL) {
  S <- .fa_covariance_names(S)
  input <- "'S'"
  if (!is.null(items)) {
    failed <- .fa_items_check(items, colnames(S))
    if (!is.null(failed)) {
      absent <- if (is.character(items)) setdiff(items[!is.na(items)], colnames(S))
      .fa_items_refuse(c("'items' must name at least three distinct items of 'S'.",
                         if (length(absent)) paste0(" Not items of 'S': ",
                                                    paste(sQuote(absent, FALSE), collapse = ", "),
                                                    ".")),
                       failed$reason, failed$items)
    }
    S <- S[items, items, drop = FALSE]
    input <- "'S' reduced to 'items'"
  }
  if (any(!is.finite(S)) ||
      !isTRUE(all.equal(S, t(S), check.attributes = FALSE)) ||
      any(diag(S) <= 0)) {
    .dmar_refuse_argument(c(
      input, " must be a finite symmetric covariance matrix of at least ",
      "three items with positive diagonal entries."), arg = NULL)
  }
  R <- stats::cov2cor(S)
  if (min(eigen(R, symmetric = TRUE, only.values = TRUE)$values) <= 0) {
    .dmar_refuse_argument(c(
      input, " must be positive definite for this factor analysis."),
      arg = NULL)
  }
  if (!is.numeric(N) || length(N) != 1L || !is.finite(N) || N < 3 || N %% 1) {
    .dmar_refuse_argument(
      "'N' must be a whole number of at least three with covariance input.",
      arg = "N")
  }
  S
}

# The shape and the names of covariance input, which hold for the whole
# matrix, since items (or a target) are read from its names.
.fa_covariance_names <- function(S) {
  if (!is.matrix(S) || !is.numeric(S) || nrow(S) < 3L || nrow(S) != ncol(S)) {
    .dmar_refuse_argument(c(
      "'S' must be a finite symmetric covariance matrix of at least ",
      "three items with positive diagonal entries."), arg = "S")
  }
  if (is.null(colnames(S)) || !identical(rownames(S), colnames(S)) ||
      anyDuplicated(colnames(S)) || any(!nzchar(colnames(S)))) {
    .dmar_refuse_argument(
      "'S' needs matching, unique row and column names.", arg = "S")
  }
  S
}

.fa_parameter_table <- function(pe, conf_level) {
  # A fit made with se = "none" (ULS from a supplied matrix) carries no
  # se, z, p-value, or limit columns; the table keeps them as NA.
  column <- function(name) {
    if (name %in% names(pe)) pe[[name]] else rep(NA_real_, nrow(pe))
  }
  out <- data.frame(term = paste(pe$lhs, pe$op, pe$rhs),
    lhs = pe$lhs, op = pe$op, rhs = pe$rhs,
    estimate = pe$est, se = column("se"), statistic = column("z"),
    p_value = column("pvalue"),
    ci_lower = column("ci.lower"), ci_upper = column("ci.upper"),
    conf_level = conf_level, stringsAsFactors = FALSE)
  out
}

# A general factor-analysis footer also survives ordinary tidy() dispatch.
#' @exportS3Method .dmar_tbl_footer dmar_factor_analysis
.dmar_tbl_footer.dmar_factor_analysis <- function(x, ...) {
  .dmar_footer_add(NextMethod(), attr(x, "analysis_note", exact = TRUE), "analysis")
}

#' @rdname dmar_tidiers
#' @exportS3Method generics::tidy
tidy.dmar_factor_analysis <- function(x, ...) {
  out <- as.data.frame(x)
  attr(out, "class") <- "data.frame"
  method <- attr(x, "method", exact = TRUE)
  out$estimator <- rep(method$estimator, nrow(out))
  out$se_method <- rep(method$se, nrow(out))
  out$missing <- rep(method$missing, nrow(out))
  out$h1_diagnostic <- rep(attr(x, "h1_diagnostic", exact = TRUE) %||% NA_character_, nrow(out))
  # fa_exploratory() records whether its standard errors and intervals were
  # skipped by request (point_estimates_only), TRUE or FALSE.
  point_only <- attr(x, "point_estimates_only", exact = TRUE)
  if (!is.null(point_only)) out$point_estimates_only <- rep(point_only, nrow(out))
  out
}

#' @rdname dmar_tidiers
#' @exportS3Method generics::glance
glance.dmar_factor_analysis <- function(x, ...) {
  method <- attr(x, "method", exact = TRUE)
  out <- data.frame(N = method$N, estimator = method$estimator,
    se_method = method$se, missing = method$missing,
    conf_level = attr(x, "conf_level"), stringsAsFactors = FALSE)
  if (!is.null(method$rotation)) out$rotation <- method$rotation
  if (!is.null(method$structure)) out$structure <- method$structure
  if (!is.null(method$sam_method)) out$sam_method <- method$sam_method
  if (!is.null(method$gamma_unbiased)) out$gamma_unbiased <- method$gamma_unbiased
  out$h1_diagnostic <- attr(x, "h1_diagnostic", exact = TRUE) %||% NA_character_
  point_only <- attr(x, "point_estimates_only", exact = TRUE)
  if (!is.null(point_only)) out$point_estimates_only <- point_only
  out
}
