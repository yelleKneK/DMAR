# Internals of fa_exploratory(): the rotation criteria it offers, the
# quantities derived from a solution (Hofmann's complexity, the variance
# each factor accounts for, and each factor's determinacy and coefficient
# H with their intervals), and the display that every print,
# knit_print(), and as_kable() of its table shows (see .dmar_tbl_display()
# and .dmar_tbl_footer() in R/dmar_output_helpers.R). The shared tidy(),
# glance(), and footer methods of the class dmar_factor_analysis, which
# the SAM interface also carries, are in R/factor_analysis_internals.R.

# The rotation criteria, by the name fa_exploratory() takes: lavaan's name
# for the criterion and whether the rotated factors are kept uncorrelated.
# "none" is the unrotated solution lavaan identifies, whose factors are
# uncorrelated. The form is passed to lavaan explicitly and checked on the
# way back, so a criterion is never applied in the other form silently.
# "target" and "target_orthogonal" are the one criterion, lavaan's "pst"
# (rotation to a partially specified target), in its two forms, each under
# a name of its own so that the table keeps one name to one form; they take
# the target in fa_exploratory()'s `target` (see .efa_target_matrix()).
# lavaan's bifactor criteria (bi-geomin, bi-quartimin) are not offered: with
# two factors their rotated solution is not identified (a singular
# covariance matrix of the estimates), and they need identification checks
# of their own.
.efa_rotations <- function() {
  data.frame(
    name = c("oblimin", "quartimin", "geomin", "simplimax", "oblimax",
             "infomax", "bentler", "target", "varimax", "quartimax",
             "target_orthogonal", "none"),
    lavaan = c("oblimin", "quartimin", "geomin", "simplimax", "oblimax",
               "infomax", "bentler", "pst", "varimax", "quartimax",
               "pst", "none"),
    orthogonal = c(rep(FALSE, 8L), rep(TRUE, 4L)),
    stringsAsFactors = FALSE)
}

# Rotation to a partially specified target (Browne, 1972a, 1972b; Browne,
# 2001, pp. 123--125). The target is a matrix B, items by factors, whose
# specified entries hold target values for the standardized pattern
# loadings and whose unspecified entries are NA; the rotation minimizes
# the sum of squared differences between the specified loadings and their
# targets. lavaan rotates to it ("pst") but checks neither the form of the
# target nor whether it identifies the rotation: an underidentified target
# comes back with finite standard errors and no warning, and one with the
# same pattern in two columns stops with a message about degenerate random
# starts. The helpers below read the target, check the identification
# counts before fitting and the rank condition on the unrotated loadings,
# and refuse what fails, through .dmar_refuse(): class
# dmar_target_invalid for a target that cannot be read as one, and
# dmar_target_unidentified for one that does not identify the rotation,
# each with a `reason` element.
.efa_target_refuse <- function(message, reason, unidentified = FALSE, ...) {
  .dmar_refuse(message, if (unidentified) "target_unidentified" else "target_invalid",
               reason = reason, ...)
}

# The target as a numeric matrix, items by factors, NA where unspecified,
# with the rows in the order of `items` and the target's factor names as
# column names (f1, f2, ... for a matrix without them). Two forms are read.
# A named list of item keys, the form fa_confirmatory() takes for its
# factors: an item keyed to a factor is targeted at zero on every factor
# it is not keyed to and unspecified on the factors it is keyed to (an
# item keyed to two factors is unspecified on both), and an analyzed item
# keyed to no factor is unspecified on every factor. A numeric matrix with
# the items as row names, in any order, and NA for an unspecified loading;
# its specified values are standardized loadings, strictly between -1 and
# 1. With `items = NULL` the items are those the target names (the keyed
# items, in order of first appearance, or the row names), which
# fa_exploratory() uses as the default analysis items.
.efa_target_matrix <- function(target, items = NULL) {
  names_text <- function(x) paste(x, collapse = ", ")
  if (is.list(target) && !is.data.frame(target)) {
    f <- names(target)
    if (is.null(f) || anyNA(f) || any(!nzchar(f)) || anyDuplicated(f)) {
      .efa_target_refuse(paste0("A list 'target' must name each factor, with ",
                                "unique names: list(factor = c(\"item\", ...))."),
                         "names")
    }
    if (length(target) < 2L) {
      .efa_target_refuse(paste0("A target rotation needs at least two factors; ",
                                "with one there is nothing to rotate."), "factors")
    }
    for (j in seq_along(target)) {
      k <- target[[j]]
      if (!is.character(k) || !length(k) || anyNA(k) || any(!nzchar(k)) ||
          anyDuplicated(k)) {
        .efa_target_refuse(sprintf(paste0("Factor '%s' of the target must name one or ",
                                          "more distinct items."), f[j]), "names")
      }
    }
    keyed <- unique(unlist(target, use.names = FALSE))
    if (is.null(items)) items <- keyed
    unknown <- setdiff(keyed, items)
    if (length(unknown)) {
      .efa_target_refuse(c("The target keys items that are not analysis items: ",
                           names_text(unknown), "."), "names")
    }
    B <- matrix(0, length(items), length(f), dimnames = list(items, f))
    for (j in seq_along(f)) B[target[[j]], j] <- NA_real_
    B[!items %in% keyed, ] <- NA_real_
  } else if (is.matrix(target) &&
             (is.numeric(target) || (is.logical(target) && all(is.na(target))))) {
    if (ncol(target) < 2L) {
      .efa_target_refuse(paste0("A target rotation needs at least two factors; ",
                                "with one there is nothing to rotate."), "factors")
    }
    rn <- rownames(target)
    if (is.null(rn) || anyNA(rn) || any(!nzchar(rn)) || anyDuplicated(rn)) {
      .efa_target_refuse(paste0("A matrix 'target' needs unique row names naming ",
                                "the items, one row per item."), "names")
    }
    if (is.null(items)) items <- rn
    if (!setequal(rn, items)) {
      absent <- setdiff(items, rn)
      extra <- setdiff(rn, items)
      .efa_target_refuse(c("The row names of the target must be the analysis items, ",
                           "in any order.",
                           if (length(absent)) c(" Items without a row: ",
                                                 names_text(absent), "."),
                           if (length(extra)) c(" Rows that are not analysis items: ",
                                                names_text(extra), ".")), "names")
    }
    f <- colnames(target)
    if (is.null(f)) f <- paste0("f", seq_len(ncol(target)))
    if (anyNA(f) || any(!nzchar(f)) || anyDuplicated(f)) {
      .efa_target_refuse("The column names of the target must be unique factor names.",
                         "names")
    }
    storage.mode(target) <- "double"
    if (any(is.nan(target) | is.infinite(target))) {
      .efa_target_refuse(paste0("The target holds non-finite values. Use NA for an ",
                                "unspecified loading and a finite value otherwise."),
                         "values")
    }
    if (any(abs(target) >= 1, na.rm = TRUE)) {
      .efa_target_refuse(paste0("Specified target values are standardized pattern ",
                                "loadings and must lie strictly between -1 and 1."),
                         "values")
    }
    B <- target[items, , drop = FALSE]
    dimnames(B) <- list(items, f)
  } else {
    .efa_target_refuse(paste0("'target' must be a named list of item keys, as ",
                              "fa_confirmatory() takes its factors, or a numeric ",
                              "matrix of items by factors with NA for an ",
                              "unspecified loading."), "form")
  }
  if (all(is.na(B))) {
    .efa_target_refuse("The target specifies no loading.", "none_specified")
  }
  B
}

# Which columns of a target fall short of the restrictions that identify
# the rotation (Algina, 1980). `n` is what each column offers: its number
# of specified entries, or the rank of the unrotated loadings at them.
# In oblique rotation, with the factor variances at one, a column of zeros
# needs m - 1 restrictions (Theorem 2; Browne, 2001, p. 115) and a column
# with a nonzero target value needs m (Theorem 1). In orthogonal rotation
# the columns are identified one after another, each also orthogonal to
# those before it, so with the columns in increasing order of what they
# offer (less one for a nonzero value) the kth needs k - 1 zeros, or k
# restrictions with a nonzero value (Theorem 3 and the remark after it);
# the smallest columns are placed first, which no other order improves.
# Returns a logical vector, TRUE for a column that falls short, with the
# number each column needs as the attribute "need".
.efa_target_short <- function(n, nonzero, orthogonal) {
  m <- length(n)
  slack <- n - nonzero
  need <- if (orthogonal) {
    position <- integer(m)
    position[order(slack)] <- seq_len(m)
    position - 1L + nonzero
  } else {
    rep(m - 1L, m) + nonzero
  }
  structure(n < need, need = need)
}

# The identification counts, before fitting: each column of the target
# must carry the restrictions .efa_target_short() asks of it, and no two
# columns may specify the same entries with the same values (the rotation
# could not tell those factors apart: an oblique rotation collapses them,
# an orthogonal one leaves the pair free to turn in its plane).
.efa_target_counts <- function(B, orthogonal) {
  m <- ncol(B)
  n <- colSums(!is.na(B))
  nonzero <- colSums(!is.na(B) & B != 0) > 0
  short <- .efa_target_short(n, nonzero, orthogonal)
  if (any(short)) {
    need <- attr(short, "need")
    rule <- if (orthogonal) {
      sprintf(paste0("In orthogonal rotation the columns of the target, taken in ",
                     "increasing order of their specified entries, need at least ",
                     "%s of them (one more in a column with a nonzero value; ",
                     "Algina, 1980)."), paste(seq_len(m) - 1L, collapse = ", "))
    } else {
      sprintf(paste0("In oblique rotation each column of the target needs at least ",
                     "m - 1 = %d specified entries (m in a column with a nonzero ",
                     "value; Browne, 2001, p. 115; Algina, 1980)."), m - 1L)
    }
    .efa_target_refuse(c(rule, " The target does not identify the rotation: ",
                         paste(sprintf("%s has %d, needs %d", colnames(B)[short],
                                       n[short], need[short]), collapse = "; "), "."),
                       "counts", unidentified = TRUE, factors = colnames(B)[short])
  }
  for (j in seq_len(m - 1L)) for (k in (j + 1L):m) {
    a <- B[, j]
    b <- B[, k]
    if (identical(is.na(a), is.na(b)) && all(a[!is.na(a)] == b[!is.na(b)])) {
      .efa_target_refuse(sprintf(paste0("Columns %s and %s of the target specify the ",
                                        "same entries with the same values, so the ",
                                        "rotation cannot tell the two factors apart."),
                                 colnames(B)[j], colnames(B)[k]),
                         "duplicate", unidentified = TRUE,
                         factors = colnames(B)[c(j, k)])
    }
  }
  invisible(B)
}

# The rank condition, on the standardized unrotated loadings A (items by
# factors, the rows in the order of B's): the rows of A at a column's
# specified entries must have the rank that column's restrictions need
# (.efa_target_short() with the rank in place of the count; Algina, 1980,
# Theorems 1 to 3). The rank is unchanged by rotation, so the unrotated
# fit serves. In orthogonal rotation the full condition also involves the
# columns identified before, which are known only after rotating; the rank
# of the column's own rows is the part checkable in advance, and it is
# necessary. Only deficiency to within the precision of the fitted
# loadings is refused: a singular value below 1e-6 times the largest. The
# loadings come from an iterative fit, so an exact deficiency shows in them
# only to the precision the optimizer reaches, not to machine precision:
# on a covariance matrix built exactly from a population structure with two
# proportional rows, the ratio of the singular values was about 3e-9 under
# ML and ULS alike, and a threshold at sqrt(.Machine$double.eps), 1.5e-8,
# would leave a margin of less than five. Sampling error keeps a sample far
# from the threshold (in 50 samples of 300 from that structure the ratio
# never fell below 0.003). A nearly deficient set of rows shows itself in
# wide intervals instead. The unit test reads the population loadings
# themselves, where the rank is exact.
.efa_target_rank <- function(B, A, orthogonal) {
  rank_of <- function(M) {
    if (!nrow(M)) return(0L)
    d <- svd(M, nu = 0L, nv = 0L)$d
    sum(d > 1e-6 * max(d))
  }
  r <- vapply(seq_len(ncol(B)), function(j) {
    rank_of(A[!is.na(B[, j]), , drop = FALSE])
  }, integer(1))
  nonzero <- colSums(!is.na(B) & B != 0) > 0
  short <- .efa_target_short(r, nonzero, orthogonal)
  if (any(short)) {
    need <- attr(short, "need")
    .efa_target_refuse(c("The unrotated loadings of the items specified in ",
                         paste(sprintf("column %s have rank %d, below the %d",
                                       colnames(B)[short], r[short], need[short]),
                               collapse = "; "),
                         " that identification needs (Algina, 1980): those items do ",
                         "not separate the factors. Specify other items, or more of ",
                         "them."),
                       "rank", unidentified = TRUE, factors = colnames(B)[short])
  }
  invisible(B)
}

# Rows of a parameter table (lhs, op, rhs, se) that hold a specified
# loading of the target which the rotation attains exactly. Where the
# target has no restriction to spare (every column at its minimum in an
# oblique rotation, which then reproduces the confirmatory model with those
# loadings fixed; in an orthogonal rotation, a column at its minimum can be
# attained while others are not), the zeros are fixed by the target rather
# than estimated, and lavaan gives them a standard error of zero (below
# sqrt(.Machine$double.eps) here). Elsewhere a
# specified loading is estimated, with a standard error of the order of
# 1/sqrt(N). `factors` are lavaan's factor names, in the order of B's
# columns.
#
# With `by_value = TRUE` (fa_exploratory(point_estimates_only = TRUE), whose
# fit has no standard errors) a specified loading is taken as attained when
# its standardized estimate lies within `tol` of its target value. The
# rotation is iterated to a gradient tolerance of 1e-8, and an attained
# loading reaches its target to about that precision, while an estimated
# one is a sampling quantity with a standard error of the order of
# 1/sqrt(N), which falls within 1e-6 of a given value with a probability
# of the order of 1e-5. On every target of the tests (oblique and
# orthogonal, zero and nonzero values, attained in full, in part, or not
# at all) the two rules flag the same loadings.
.efa_target_attained <- function(tab, B, factors, by_value = FALSE, tol = 1e-6) {
  cell <- cbind(match(tab$rhs, rownames(B)), match(tab$lhs, factors))
  load <- tab$op == "=~" & !is.na(cell[, 1L]) & !is.na(cell[, 2L])
  b <- rep(NA_real_, nrow(tab))
  b[load] <- B[cell[load, , drop = FALSE]]
  out <- if (by_value) {
    load & !is.na(b) & is.finite(tab$estimate) & abs(tab$estimate - b) < tol
  } else {
    load & !is.na(b) & is.finite(tab$se) & tab$se < sqrt(.Machine$double.eps)
  }
  attr(out, "target") <- b
  out
}

# Shows the rows `rows` of a parameter table as fixed at `value`: the
# estimate at the value, both limits at it, and no standard error,
# statistic, or p-value.
.efa_fix_rows <- function(tab, rows, value) {
  tab$estimate[rows] <- value
  tab$ci_lower[rows] <- tab$ci_upper[rows] <- value
  tab$se[rows] <- tab$statistic[rows] <- tab$p_value[rows] <- NA_real_
  tab
}

# The error fa_exploratory() raises for a fit without a reportable solution,
# of class dmar_efa_no_solution with a `reason` element ("not converged",
# "improper", or "standard errors unavailable") that factor_count() reads
# to report the count as a row rather than stop. The refusal of a fit whose
# likelihood has no maximum carries the class too, with the reason "no
# maximum" (.efa_no_maximum() below).
#
# With `point_only = TRUE` (fa_exploratory(point_estimates_only = TRUE))
# the reason "standard errors unavailable" comes from
# .efa_information_singular() below, not from standard errors, and the
# message says so; the class and the reason are the same, so that
# factor_count() and a program reading the reason treat the two alike.
.efa_no_solution <- function(reason, point_only = FALSE) {
  text <- if (point_only && identical(reason, "standard errors unavailable")) {
    paste0("The expected information of the exploratory model could not be ",
           "inverted at its estimates, so the solution is not identified ",
           "there: its estimates are not unique, and the full analysis ",
           "would have no standard errors. No table is returned.")
  } else switch(reason,
    "not converged" = paste0("The exploratory model did not converge. Review ",
                             "the items and factor count before interpreting it."),
    improper = paste0("The exploratory model converged to an improper solution ",
                      "(lavaan's post-estimation check failed: a residual ",
                      "variance below zero, or a covariance matrix that is not ",
                      "positive definite). Review the items and factor count ",
                      "before interpreting it."),
    "standard errors unavailable" = paste0("Rotation-aware loading standard ",
                                           "errors were unavailable; no ordinary ",
                                           "interval table is returned."))
  .dmar_refuse(
    text,
    "efa_no_solution", reason = reason)
}

# .efa_information_singular(fit)
#
# fa_exploratory(point_estimates_only = TRUE) computes no standard errors,
# so the refusal of a solution whose rotation-aware errors are unavailable
# has nothing to read. That refusal is the full analysis's guard against a
# solution that is not identified at its estimates, whose information
# cannot be inverted and whose estimates are therefore not unique. In its
# place the information of the unrotated fit `fit` (identified by lavaan's
# fixed zeros, without constraints) is inverted as lavaan inverts the
# information behind its errors: TRUE when it cannot be, or when the
# inverse has a diagonal entry that is not finite and positive, where
# lavaan's errors would be missing. The expected information is tried
# first, being the cheap one: 0.14 s on 18 items, 1,000 cases, and 176
# missing data patterns under FIML, against 1.7 s for the observed
# information and 3.6 s for the fit itself. Only when it cannot be
# inverted is the information lavaan's errors use (its default for the
# fit: the observed under FIML and as the bread of MLR, the expected
# otherwise) inverted as well, and TRUE requires that to fail too, so that
# no fit is refused whose unrotated errors the full analysis would have.
# (The two differ at a stationary point where loadings vanish: one
# factor on two blocks of three items never observed together, which
# lavaan 0.7-3 starts at zero loadings for one block, has a singular
# expected and an invertible observed information, and the full analysis
# reports it.) A failure confined to the rotation's part of the
# rotation-aware errors, or to the delta method of the standardized
# solution, is not detected.
.efa_information_singular <- function(fit) {
  inverse <- function(what) {
    tryCatch(lavaan::lavInspect(fit, what), error = function(e) NULL)
  }
  .efa_inverse_unusable(inverse("inverted.information.expected")) &&
    .efa_inverse_unusable(inverse("inverted.information"))
}

# Whether `inv`, an inverted information matrix as lavaan returns it (a
# "try-error" string when solve() failed), is unusable.
.efa_inverse_unusable <- function(inv) {
  if (!is.matrix(inv) || !is.numeric(inv) || !length(inv)) return(TRUE)
  d <- diag(inv)
  any(!is.finite(d) | d <= 0)
}

# .efa_no_maximum(nm, efa_args, select = NULL)
#
# The refusal of an exploratory fit whose likelihood has no maximum
# (.lavaan_no_maximum()'s list `nm`): .no_max_refuse()'s condition under
# fa_exploratory()'s name, with the lavaan::efa() call that reproduces the
# unrotated fit (`efa_args`, as fa_exploratory() passes them to it). The
# placeholder of the input selects the items `select` from it, when they
# are not the whole of the data or of the covariance matrix the user
# passed: lavaan::efa() fits every column of a data frame and every row of
# a covariance matrix it is given (.no_max_call_text()). It carries the
# class dmar_efa_no_solution as well, with the reason "no maximum", so
# that a caller that reads fa_exploratory()'s refusals of a fit by that
# class (factor_count() among them) reads this one too.
.efa_no_maximum <- function(nm, efa_args, select = NULL) {
  pick <- if (!is.null(select)) list(data = select, sample_cov = select)
  .no_max_refuse("fa_exploratory", nm, syntax = NULL,
                 lavaan_call = .no_max_call_text("lavaan::efa", efa_args,
                                                 select = pick),
                 reason = "no maximum", extra_class = "dmar_efa_no_solution")
}

# .efa_resolve_missing(missing, omitted, data, ord, estimator)
#
# The missing data method of an exploratory fit of raw data, resolved and
# checked once, for fa_exploratory() and for factor_count(), which fits
# each count through it and forms its correlations and its parallel
# analysis under the method resolved here. `missing` is the matched
# value, `omitted` whether the caller left the argument out (read with
# base::missing() before match.arg()), `data` the analysis items, `ord`
# the names of the ordered items, and `estimator` the matched estimator.
#
# FIML, the default, is a likelihood for continuous items. With ordered
# items and no value missing, the default fits the items as they stand
# (listwise deletion with nothing to delete), as fa_confirmatory() does;
# with a value missing, the choice between pairwise and listwise deletion
# is the user's, and nothing is deleted unasked. ULS has no likelihood: it
# takes complete data, an explicit listwise choice, or a correlation
# matrix estimated elsewhere (through S). Continuous items take ML or MLR
# with FIML or listwise deletion. Returns the method, "listwise" in place
# of an omitted default for complete ordered items and the matched value
# otherwise; refuses (dmar_argument_invalid) what the fit does not take.
#
# Earlier builds of 1.1.0 (factor_count() is new in it) matched
# factor_count()'s own default, "fiml", and passed it on explicitly, so
# complete ordered items that fa_exploratory() fits with the argument left
# out were refused by factor_count().
.efa_resolve_missing <- function(missing, omitted, data, ord, estimator) {
  if (length(ord)) {
    if (omitted && missing %in% c("fiml", "ml") && !anyNA(data)) {
      missing <- "listwise"
    }
    if (estimator != "WLSMV" || !missing %in% c("pairwise", "listwise")) {
      .dmar_refuse_argument(c(
        "Ordered EFA requires estimator = \"WLSMV\" and, when item ",
        "values are missing, an explicit missing = \"pairwise\" or ",
        "\"listwise\" choice. Continuous normal-theory FIML is not an ",
        "ordered-data likelihood."), arg = NULL)
    }
  } else if (estimator == "ULS") {
    if (missing == "pairwise") {
      .dmar_refuse_argument(c(
        "estimator = \"ULS\" does not take pairwise deletion: lavaan ",
        "fits the correlation structure to complete cases. Use missing ",
        "= \"listwise\", or supply through 'S' a correlation matrix ",
        "estimated from the incomplete data (cov_fiml() gives the ",
        "maximum likelihood estimate)."), arg = "S")
    }
    if (missing %in% c("fiml", "ml") && anyNA(data)) {
      .dmar_refuse_argument(c(
        "estimator = \"ULS\" has no likelihood, so FIML does not apply ",
        "and the incomplete cases are not deleted silently. Use missing ",
        "= \"listwise\" to fit the complete cases, or supply through 'S' ",
        "a correlation matrix estimated from every case (cov_fiml() ",
        "gives the maximum likelihood estimate), or use estimator = ",
        "\"ML\" or \"MLR\" with FIML."), arg = "S")
    }
  } else if (estimator == "WLSMV" || missing == "pairwise") {
    .dmar_refuse_argument(
      "Continuous EFA uses ML or MLR with FIML or listwise deletion.",
      arg = NULL)
  }
  missing
}

# Hofmann's (1978) complexity of each item from its pattern loadings: the
# square of the sum of squared loadings over the sum of the loadings to
# the fourth power. 1 for an item loading on one factor alone; the number
# of factors for an item loading equally on all of them. NA for a row of
# zero loadings.
.efa_complexity <- function(L) {
  s2 <- rowSums(L^2)
  out <- s2^2 / rowSums(L^4)
  out[!is.finite(out)] <- NA_real_
  unname(out)
}

# The variance each factor accounts for: the diagonal of Phi Lambda'
# Lambda, with Lambda the standardized pattern loadings and Phi the factor
# correlation matrix (the column sums of squared loadings when the factors
# are uncorrelated). The entries sum to the sum of the communalities, the
# trace of Lambda Phi Lambda'. Each is also a proportion of the total
# variance of the items (their number, in the standardized metric) and of
# the common variance. With more than one factor a total row closes the
# table; no cumulative column, since the order of the factors is
# arbitrary.
.efa_variance_accounted <- function(L, Phi) {
  ss <- as.numeric(diag(Phi %*% crossprod(L)))
  p <- nrow(L)
  total <- sum(ss)
  out <- data.frame(factor = colnames(L), sum_sq_loadings = ss,
                    proportion_variance = ss / p,
                    proportion_common = ss / total,
                    stringsAsFactors = FALSE)
  if (ncol(L) > 1L) {
    out <- rbind(out, data.frame(factor = "total", sum_sq_loadings = total,
                                 proportion_variance = total / p,
                                 proportion_common = 1,
                                 stringsAsFactors = FALSE))
  }
  rownames(out) <- NULL
  out
}

# Each factor's determinacy and coefficient H, with delta method standard
# errors and logit-scale Wald limits.
#
# The determinacy of a factor is the correlation of its regression factor
# score with the factor, the multiple correlation of the factor with the
# items (Grice, 2001; the factor determinacy index of Ferrando &
# Lorenzo-Seva, 2018): the square root of the diagonal of
# Phi Lambda' R^{-1} Lambda Phi, with Lambda the standardized pattern
# loadings, Phi the factor correlations, and R the item correlation
# matrix the model implies (Lambda Phi Lambda' with unit diagonal). Its
# square is the maximal reliability of the factor under the model, the
# largest reliability any weighted sum of the items reaches as a measure
# of it (Ferrando and Lorenzo-Seva's generalized H; the
# maximal_reliability attribute of bifactor_indices()). Coefficient H is
# Hancock and Mueller's (2001) closed form on the factor's standardized
# loadings alone, t / (1 + t) with t the sum of l^2 / (1 - l^2), as
# Rodriguez, Reise, and Haviland (2016) apply it to one factor of a
# multidimensional model: it counts 1 - l^2 as each item's error
# variance, the other factors' common variance included. With one factor
# the two coincide (H is the squared determinacy); otherwise they can
# differ in either direction. H is NA for a factor with a standardized
# loading at or beyond one in absolute value, as reliability_H() treats
# an improper solution.
#
# Both are functions of the standardized loadings and factor
# correlations alone, so their delta method variances are J V J', with V
# lavaan's covariance matrix of the standardized free parameters
# (lavInspect(fit, "vcov.std.all"), the matrix behind the standardized
# loading and communality intervals of the table) and J the Jacobian by
# central differences with step 1e-6 * max(1, |z|), as .H_from_fit()
# takes it. Written as functions of the loadings and factor correlations
# alone, with the residual variances implied by the unit diagonal, they
# agree with the same quantities written over lavaan's free parameters
# and vcov(fit) to the precision of the differences (checked to 1e-9 on
# the Holzinger and Swineford items). The se column is the standard error
# of the index itself; the limits are the Wald limits on the logit scale
# transformed back (Browne, 1982), as reliability_H() forms them by
# default, so they lie within (0, 1). A fit without standard errors (ULS
# from a supplied matrix) gives the estimates with se and limits NA.
#
# `ss` is lavaan::standardizedSolution(fit) as fa_exploratory() has it.
# Returns a list: `rows`, in the layout of the parameter table (lhs and
# rhs the factor, op "determinacy" or "H", metric "standardized"), and
# `table`, one row per factor with both indices and their intervals (the
# factor_quality attribute).
.efa_factor_quality <- function(fit, ss, conf_level) {
  pt <- lavaan::parTable(fit)
  free <- pt[pt$free > 0L, , drop = FALSE]
  key <- function(t) paste(t$lhs, t$op, t$rhs)
  z0 <- ss$est.std[match(key(free), key(ss))]
  factors <- lavaan::lavNames(fit, "lv")
  items <- lavaan::lavNames(fit, "ov")
  k <- length(factors)
  load <- free$op == "=~"
  phi <- free$op == "~~" & free$lhs %in% factors & free$rhs %in% factors &
    free$lhs != free$rhs
  cells <- cbind(match(free$rhs[load], items), match(free$lhs[load], factors))
  pairs <- cbind(match(free$lhs[phi], factors), match(free$rhs[phi], factors))
  indices <- function(z) {
    L <- matrix(0, length(items), k)
    L[cells] <- z[load]
    Phi <- diag(1, k)
    if (nrow(pairs)) {
      Phi[pairs] <- z[phi]
      Phi[pairs[, 2:1, drop = FALSE]] <- z[phi]
    }
    R <- L %*% Phi %*% t(L)
    diag(R) <- 1
    RinvL <- tryCatch(solve(R, L), error = function(e) NULL)
    rho <- if (is.null(RinvL)) rep(NA_real_, k) else
      sqrt(pmax(diag(Phi %*% crossprod(L, RinvL) %*% Phi), 0))
    H <- apply(L, 2L, function(l) {
      if (any(abs(l) >= 1)) return(NA_real_)
      t <- sum(l^2 / (1 - l^2))
      t / (1 + t)
    })
    c(rho, H)
  }
  est <- indices(z0)
  se <- rep(NA_real_, 2L * k)
  V <- if (!identical(lavaan::lavInspect(fit, "options")$se, "none")) {
    tryCatch(lavaan::lavInspect(fit, "vcov.std.all"), error = function(e) NULL)
  }
  if (!is.null(V) && all(dim(V) == length(z0))) {
    J <- matrix(0, length(est), length(z0))
    for (i in which(load | phi)) {
      h <- 1e-6 * max(1, abs(z0[i]))
      up <- dn <- z0
      up[i] <- up[i] + h
      dn[i] <- dn[i] - h
      J[, i] <- (indices(up) - indices(dn)) / (2 * h)
    }
    se <- sqrt(pmax(diag(J %*% unclass(V) %*% t(J)), 0))
  }
  se[!is.finite(se)] <- NA_real_
  limits <- lapply(seq_along(est), function(i) {
    .H_wald(est[i], se[i], conf_level, logistic = TRUE)
  })
  lower <- vapply(limits, `[[`, numeric(1), "lower")
  upper <- vapply(limits, `[[`, numeric(1), "upper")
  index <- rep(c("determinacy", "H"), each = k)
  factor <- rep(factors, 2L)
  rows <- data.frame(term = paste(factor, index), lhs = factor, op = index,
    rhs = factor, estimate = est, se = se, statistic = NA_real_,
    p_value = NA_real_, ci_lower = lower, ci_upper = upper,
    conf_level = conf_level, metric = "standardized",
    stringsAsFactors = FALSE)
  d <- seq_len(k)
  h <- k + d
  table <- data.frame(factor = factors,
    determinacy = est[d], determinacy_se = se[d],
    determinacy_ci_lower = lower[d], determinacy_ci_upper = upper[d],
    min_score_correlation = 2 * est[d]^2 - 1,
    min_score_correlation_se = 4 * est[d] * se[d],
    min_score_correlation_ci_lower = 2 * lower[d]^2 - 1,
    min_score_correlation_ci_upper = 2 * upper[d]^2 - 1,
    H = est[h], H_se = se[h], H_ci_lower = lower[h], H_ci_upper = upper[h],
    stringsAsFactors = FALSE)
  list(rows = rows, table = table)
}

# The determinacy and H block of the display, one row per factor with
# "estimate [lower, upper]" cells (the estimate alone where the limits are
# NA), read from the table's rows; NULL when the table has no such rows.
.efa_quality_display <- function(x, digits = getOption("dmar.digits", 3L)) {
  det <- x[x$op == "determinacy", , drop = FALSE]
  H <- x[x$op == "H", , drop = FALSE]
  if (!nrow(det) || !nrow(H)) return(NULL)
  H <- H[match(det$lhs, H$lhs), , drop = FALSE]
  data.frame(factor = det$lhs,
             determinacy = .efa_cells(det$estimate, det$ci_lower, det$ci_upper, digits),
             min_r = .efa_cells(2 * det$estimate^2 - 1, 2 * det$ci_lower^2 - 1,
                                2 * det$ci_upper^2 - 1, digits),
             H = .efa_cells(H$estimate, H$ci_lower, H$ci_upper, digits),
             stringsAsFactors = FALSE)
}

# The standardized solution as matrices, read from the rows of the table
# (any data frame with lhs, op, rhs, metric, estimate, ci_lower, and
# ci_upper): the pattern loadings with their limits, items by factors, and
# the factor correlations with theirs, the diagonal 1 and the limits NA.
.efa_matrices <- function(x) {
  s <- x[x$metric == "standardized", , drop = FALSE]
  load <- s[s$op == "=~", , drop = FALSE]
  factors <- unique(load$lhs)
  items <- unique(load$rhs)
  cell <- cbind(match(load$rhs, items), match(load$lhs, factors))
  mat <- function(v) {
    m <- matrix(NA_real_, length(items), length(factors),
                dimnames = list(items, factors))
    m[cell] <- v
    m
  }
  Phi <- diag(1, length(factors))
  dimnames(Phi) <- list(factors, factors)
  Phi_lo <- Phi_hi <- Phi
  Phi_lo[] <- Phi_hi[] <- NA_real_
  off <- s[s$op == "~~" & s$lhs %in% factors & s$rhs %in% factors &
             s$lhs != s$rhs, , drop = FALSE]
  if (nrow(off)) {
    i <- cbind(match(off$lhs, factors), match(off$rhs, factors))
    j <- i[, 2:1, drop = FALSE]
    Phi[i] <- Phi[j] <- off$estimate
    Phi_lo[i] <- Phi_lo[j] <- off$ci_lower
    Phi_hi[i] <- Phi_hi[j] <- off$ci_upper
  }
  list(items = items, factors = factors, L = mat(load$estimate),
       lo = mat(load$ci_lower), hi = mat(load$ci_upper),
       Phi = Phi, Phi_lo = Phi_lo, Phi_hi = Phi_hi)
}

# "estimate [lower, upper]" for each element, the estimate alone where
# the limits are NA (a fit without standard errors). Shapes follow `est`.
# Every number is a measured value, an estimate or a limit, so a whole
# one keeps its decimals ("0.000 [0.000, 0.075]" for an RMSEA of 0, as
# .dmar_measured_cells() treats an interval in a table).
.efa_cells <- function(est, lo, hi, digits = getOption("dmar.digits", 3L)) {
  e <- .format_dmar_value(c(est), digits, whole = FALSE)
  e[is.na(e)] <- "NA"
  ok <- !is.na(c(lo)) & !is.na(c(hi))
  out <- e
  out[ok] <- sprintf("%s [%s, %s]", e[ok],
                     .format_dmar_value(c(lo)[ok], digits, whole = FALSE),
                     .format_dmar_value(c(hi)[ok], digits, whole = FALSE))
  if (is.matrix(est)) {
    out <- matrix(out, nrow(est), ncol(est), dimnames = dimnames(est))
  }
  out
}

# Whether a table still carries what the display reads. A subset of the
# table (x[rows, cols]) keeps the class but not the attributes, and it
# then prints as an ordinary DMAR table.
.efa_display_ok <- function(x) {
  m <- attr(x, "method", exact = TRUE)
  is.list(m) && !is.null(m$rotation) && !is.null(attr(x, "fit", exact = TRUE)) &&
    all(c("lhs", "op", "rhs", "metric", "estimate", "ci_lower", "ci_upper") %in%
          names(x)) &&
    any(x$op == "=~" & x$metric == "standardized") &&
    any(x$op == "communality")
}

.efa_estimator_label <- function(m) {
  switch(m$estimator,
         ULS = "ULS (minimum residual: unweighted least squares on the correlation matrix)",
         MLR = "MLR (ML estimates, Huber-White sandwich standard errors)",
         m$estimator)
}

.efa_missing_label <- function(m) {
  switch(m$missing,
         fiml = "FIML",
         listwise = "listwise deletion",
         pairwise = "pairwise deletion",
         "covariance input" = "none (covariance input)",
         m$missing)
}

.efa_rotation_label <- function(m) {
  if (is.matrix(m$target)) {
    spec <- m$target[!is.na(m$target)]
    return(sprintf("target (%s), %d of %d loadings specified (%s)",
                   if (isTRUE(m$orthogonal)) "orthogonal" else "oblique",
                   length(spec), length(m$target),
                   if (all(spec == 0)) "all 0" else
                     sprintf("%d at 0, %d nonzero", sum(spec == 0), sum(spec != 0))))
  }
  paste0(m$rotation,
         if (identical(m$rotation, "none")) " (unrotated)"
         else if (isTRUE(m$orthogonal)) " (orthogonal)" else " (oblique)",
         if (identical(m$row_weights, "kaiser")) ", Kaiser row weights")
}

.efa_se_label <- function(m) {
  switch(m$se,
         standard = "standard",
         robust.huber.white = "Huber-White sandwich",
         robust.sem.nt = "sandwich with the normal theory covariance matrix of the correlations (robust.sem.nt)",
         robust.sem = "robust.sem",
         none = "none",
         m$se)
}

# The fit lines beneath the display, read from the fit_measures attribute
# and its version lines (see sem_fit_measures()). Under ULS the statistic
# has no reference distribution, so the line gives it with the NFI and the
# SRMR; the note under the table says why the other indices are NA.
.efa_fit_lines <- function(x, digits = getOption("dmar.digits", 3L)) {
  fm <- attr(x, "fit_measures", exact = TRUE)
  if (is.null(fm)) {
    return(paste("Fit measures: not reported (the unrestricted H1 information",
                 "they need is unavailable; see the note)."))
  }
  # Every measure but the degrees of freedom is measured, so a CFI of 1
  # reads 1.000 beside the other indices, as it does in the table of
  # sem_fit_measures().
  value <- function(term) {
    v <- fm$value[match(term, fm$term)]
    if (length(v) != 1L || is.na(v)) "NA" else
      .format_dmar_value(v, digits, whole = identical(term, "df"))
  }
  p_value <- function(term) {
    v <- fm$value[match(term, fm$term)]
    if (length(v) != 1L || is.na(v)) return("p = NA")
    p <- .format_dmar_pvalue(v)
    if (startsWith(p, "<")) paste("p", p) else paste("p =", p)
  }
  uls <- identical(attr(x, "method", exact = TRUE)$estimator, "ULS")
  line <- if (uls) {
    sprintf(paste("Fit: ULS discrepancy statistic %s on %s df (no reference",
                  "distribution, so no p-value, RMSEA, CFI, or TLI); NFI %s; SRMR %s"),
            value("chi_square"), value("df"), value("nfi"), value("srmr"))
  } else {
    # The RMSEA interval is at its own level (rmsea_conf_level), which the
    # line names as factor_count() names it, so that it is not read at the
    # level of the intervals in the table above.
    level <- attr(fm, "conf_level", exact = TRUE)
    rmsea_label <- if (is.numeric(level) && length(level) == 1L &&
                       is.finite(level)) {
      sprintf("RMSEA with %s%% interval", .format_conf_pct(level))
    } else "RMSEA"
    sprintf("Fit: chi square %s, df %s, %s; %s %s [%s, %s]; CFI %s; TLI %s; SRMR %s",
            value("chi_square"), value("df"), p_value("p_chi_square"),
            rmsea_label, value("rmsea"), value("rmsea_ci_lower"),
            value("rmsea_ci_upper"), value("cfi"), value("tli"), value("srmr"))
  }
  footer <- attr(fm, "sem_fit", exact = TRUE)$footer
  c(line, footer[seq_len(min(2L, length(footer)))],
    "The fit_measures attribute holds every measure with its version lines.")
}

#' @exportS3Method .dmar_tbl_display dmar_fa_exploratory
#' @noRd
.dmar_tbl_display.dmar_fa_exploratory <- function(x, digits = getOption("dmar.digits", 3L), ...) {
  if (!.efa_display_ok(x)) return(NextMethod())
  q <- .efa_matrices(x)
  cells <- .efa_cells(q$L, q$lo, q$hi, digits)
  disp <- data.frame(item = q$items, stringsAsFactors = FALSE)
  for (f in q$factors) disp[[f]] <- cells[, f]
  h <- x[x$op == "communality", , drop = FALSE]
  h <- h[match(q$items, h$lhs), , drop = FALSE]
  disp$communality <- .efa_cells(h$estimate, h$ci_lower, h$ci_upper, digits)
  # Hofmann's complexity is a measured index, 1 for an item on one factor,
  # so it keeps its decimals when it is 1, as the loadings beside it do.
  complexity <- .format_dmar_value(.efa_complexity(q$L), digits, whole = FALSE)
  complexity[is.na(complexity)] <- "NA"
  disp$complexity <- complexity
  .as_dmar_tbl(disp)
}

#' @exportS3Method .dmar_tbl_footer dmar_fa_exploratory
#' @noRd
.dmar_tbl_footer.dmar_fa_exploratory <- function(x, digits = getOption("dmar.digits", 3L), ...) {
  if (!.efa_display_ok(x)) return(NextMethod())
  m <- attr(x, "method", exact = TRUE)
  q <- .efa_matrices(x)
  k <- length(q$factors)
  p <- length(q$items)
  lines <- .dmar_footer_add(character(0), c(
    sprintf("Exploratory factor analysis: %d %s, %d items, N = %d",
            k, if (k == 1L) "factor" else "factors", p, m$N),
    sprintf("Estimator: %s; missing data: %s; rotation: %s; SE: %s",
            .efa_estimator_label(m), .efa_missing_label(m),
            .efa_rotation_label(m),
            if (isTRUE(attr(x, "point_estimates_only", exact = TRUE))) {
              "not computed (point_estimates_only = TRUE)"
            } else {
              .efa_se_label(m)
            })), "header")
  level <- .format_conf_pct(attr(x, "conf_level", exact = TRUE))
  no_ci <- all(is.na(q$lo))
  lines <- .dmar_footer_add(lines, if (no_ci) {
    paste("(standardized pattern loadings, communalities, and Hofmann's",
          "complexity; no intervals, see the note; tidy(x) has every",
          "parameter in both metrics)")
  } else {
    sprintf(paste("(standardized pattern loadings and communalities with %s%%",
                  "intervals, and Hofmann's complexity; tidy(x) has every",
                  "parameter in both metrics)"), level)
  }, "table_note")
  if (is.matrix(m$target) && identical(ncol(m$target), k)) {
    # The target beneath the loadings it was rotated toward, in Browne's
    # (2001, p. 124) layout: the value where an entry is specified, "?"
    # where it is not. The columns carry lavaan's factor names, as the
    # table does, and a line maps them to the names the target gave.
    B <- m$target[match(q$items, rownames(m$target)), , drop = FALSE]
    shown <- matrix(ifelse(is.na(B), "?", .format_dmar_value(c(B), digits)),
                    nrow(B), k, dimnames = list(q$items, q$factors))
    lines <- .dmar_footer_add(lines,
      "Target (the value where a loading is specified, ? where it is not):", "target")
    lines <- .dmar_footer_add(
      lines, utils::capture.output(print(shown, quote = FALSE, right = TRUE)),
      "target_table")
    if (!identical(colnames(B), q$factors)) {
      lines <- .dmar_footer_add(lines, paste0(
        "Factors: ", paste(q$factors, "=", colnames(B), collapse = ", "),
        " (the table keeps lavaan's names)."), "target")
    }
  }
  va <- attr(x, "variance_accounted", exact = TRUE)
  if (is.null(va)) va <- .efa_variance_accounted(q$L, q$Phi)
  shown <- va
  shown[-1L] <- lapply(shown[-1L], .format_dmar_value, digits = digits,
                       whole = FALSE)
  lines <- .dmar_footer_add(lines, paste(
    "Variance accounted for (sums of squared standardized loadings, the",
    "diagonal of Phi Lambda' Lambda):"), "variance")
  lines <- .dmar_footer_add(
    lines, utils::capture.output(print(shown, row.names = FALSE, right = TRUE)),
    "variance_table")
  if (k > 1L) {
    if (isTRUE(m$orthogonal)) {
      lines <- .dmar_footer_add(lines, if (identical(m$rotation, "none")) {
        "Factor correlations: fixed at zero (unrotated solution)."
      } else {
        "Factor correlations: fixed at zero (orthogonal rotation)."
      }, "factor_correlations")
    } else {
      mat <- .efa_cells(q$Phi, q$Phi_lo, q$Phi_hi, digits)
      diag(mat) <- "1"
      lines <- .dmar_footer_add(lines, if (no_ci) "Factor correlations:" else
        sprintf("Factor correlations with %s%% intervals:", level),
        "factor_correlations")
      lines <- .dmar_footer_add(
        lines, utils::capture.output(print(mat, quote = FALSE, right = TRUE)),
        "factor_correlations_table")
    }
  }
  fq <- .efa_quality_display(x, digits)
  if (!is.null(fq)) {
    quality_ci <- !all(is.na(c(fq$determinacy, fq$H))) &&
      any(grepl("[", c(fq$determinacy, fq$H), fixed = TRUE))
    lines <- .dmar_footer_add(lines, paste0(
      "Factor determinacy and coefficient H",
      if (length(m$ordered)) " of the latent responses",
      if (quality_ci) sprintf(" with %s%% intervals", level), ":"), "factor_quality")
    lines <- .dmar_footer_add(
      lines, utils::capture.output(print(fq, row.names = FALSE, right = TRUE)),
      "factor_quality_table")
    lines <- .dmar_footer_add(lines, paste0(
      "(determinacy: the correlation of the regression factor score with its ",
      "factor, whose square is the factor's maximal reliability; min_r: the ",
      "lowest correlation two sets of scores equally consistent with the model ",
      "can have, 2 x determinacy^2 - 1 (Guttman, 1955); H: Hancock ",
      "and Mueller's reliability of the optimally weighted composite from the ",
      "factor's loadings alone",
      if (quality_ci) "; limits from the logit scale", ")"), "factor_quality")
  }
  lines <- .dmar_footer_add(lines, .efa_fit_lines(x, digits), "fit")
  rest <- NextMethod()
  if (isTRUE(attr(x, "point_estimates_only", exact = TRUE))) {
    # Without standard errors the table has no interval, so the level of
    # its intervals (conf_level) is not stated; the only interval left is
    # the RMSEA's, whose level (rmsea_conf_level) the Fit line names.
    rest <- rest[names(rest) != "conf_level"]
  }
  c(lines, rest)
}
