# Internals of factor_congruence(): reading a solution from each kind of
# input, aligning the items of two solutions, and the one-to-one
# assignment that matches their factors. The assignment solver is tested
# against brute force in tests/testthat/test-factor_congruence.R.

# One solution as a list: `L`, the loadings (items in rows, factors in
# columns, with names); `Phi`, the factor correlations, or NULL for a bare
# matrix, which carries none; `source`, the kind of input; and `named`,
# whether the input named its items. `arg` is the argument's name ("x" or
# "y") for the messages. Fits are read in the standardized metric:
# congruence is unchanged by rescaling a whole column of loadings but not
# by rescaling an item, and an item's original-metric loadings scale with
# its standard deviation, which differs from sample to sample.
.congruence_solution <- function(obj, arg, loadings) {
  accepted <- paste0(
    "a numeric loading matrix (items in rows, factors in columns), a ",
    "numeric vector (one factor), an fa_exploratory() table, or a ",
    "single-group, single-level lavaan fit")
  # cfa_1() and cfa_2() return fa_confirmatory()'s tables, so the class
  # does not say which of the three was called, and the message names all
  # three. Under output = "standardized" they return lavaan's own table
  # (class lavaan.data.frame, as from lavaan::standardizedSolution()).
  if (inherits(obj, "dmar_fa_confirmatory")) {
    .dmar_refuse_argument(
      sprintf(paste0(
   "'%s' is a table of results from fa_confirmatory(), cfa_1(), or ",
   "cfa_2(), not the fit. Call the same function again with ",
   "output = \"fit\" and pass the lavaan fit it returns; its ",
   "standardized loadings and factor correlations are then used."),
   arg), arg = arg)
  }
  if (inherits(obj, "lavaan.data.frame")) {
    .dmar_refuse_argument(
      sprintf(paste0(
   "'%s' is a lavaan table of parameter estimates, such as ",
   "output = \"standardized\" gives in fa_confirmatory(), cfa_1(), and ",
   "cfa_2(), not a fit. Pass the lavaan fit itself (output = \"fit\" in ",
   "those functions); its standardized loadings and factor correlations ",
   "are then used."), arg), arg = arg)
  }
  if (inherits(obj, "dmar_fa_exploratory")) {
    s <- .congruence_from_efa(obj, arg)
  } else if (inherits(obj, "efaList")) {
    .dmar_refuse_argument(
      sprintf(paste0(
   "'%s' holds several lavaan EFA solutions (an efaList). Pass one of ",
   "them, for example %s[[1]], or call lavaan::efa() with ",
   "output = \"lavaan\"."), arg, arg), arg = arg)
  } else if (inherits(obj, "lavaan")) {
    s <- .congruence_from_lavaan(obj, arg)
  } else if (is.data.frame(obj)) {
    .dmar_refuse_argument(
      sprintf("'%s' is a data frame. Supply %s.", arg, accepted), arg = arg)
  } else if (is.numeric(obj) && (is.null(dim(obj)) || is.matrix(obj))) {
    # A "loadings" object (psych::fa()$loadings, factanal()$loadings) is a
    # matrix with a class and a print method; its numbers are the matrix.
    if (is.null(dim(obj))) {
      L <- matrix(as.numeric(obj), ncol = 1L,
                  dimnames = list(names(obj), NULL))
    } else {
      L <- matrix(as.numeric(obj), nrow(obj), ncol(obj),
                  dimnames = dimnames(obj))
    }
    s <- list(L = L, Phi = NULL, source = "matrix")
  } else {
    .dmar_refuse_argument(sprintf("'%s' must be %s.", arg, accepted), arg = arg)
  }
  s <- .congruence_check(s, arg)
  if (loadings == "structure") {
    if (is.null(s$Phi)) {
      .dmar_refuse_argument(
        sprintf(paste0(
     "loadings = \"structure\" needs the factor correlations, which ",
     "the matrix '%s' does not carry. Pass the structure matrix itself ",
     "(the pattern times the factor correlations) under the default ",
     "loadings = \"pattern\"; a matrix is used as supplied."), arg), arg = NULL)
    }
    s$L <- s$L %*% s$Phi
    dimnames(s$L) <- list(rownames(s$L), colnames(s$Phi))
    s <- .congruence_check(s, arg)
  }
  s
}

# The standardized pattern loadings and factor correlations of an
# fa_exploratory() table, read from its rows by .efa_matrices() (see
# R/fa_exploratory_internals.R), as the display reads them.
.congruence_from_efa <- function(obj, arg) {
  need <- c("lhs", "op", "rhs", "metric", "estimate", "ci_lower", "ci_upper")
  if (!all(need %in% names(obj)) ||
      !any(obj$op == "=~" & obj$metric == "standardized")) {
    .dmar_refuse_argument(
      sprintf(paste0(
   "'%s' does not hold the standardized loadings of an fa_exploratory() ",
   "table (a subset of the table can drop them). Pass the whole table."),
   arg), arg = arg)
  }
  q <- .efa_matrices(obj)
  list(L = q$L, Phi = q$Phi, source = "fa_exploratory")
}

# The standardized loadings and the model implied factor correlations of a
# lavaan fit. The items are the observed indicators of the factors and the
# factors are lavaan's latent variables, so an observed variable that only
# takes part in the structural model, which lavaan represents by a phantom
# factor in its loading matrix, is neither. The correlations are
# lavInspect(fit, "cor.lv"); the standardized psi matrix equals them in a
# measurement model, but in a model with regressions among the factors it
# holds residual variances and covariances instead.
.congruence_from_lavaan <- function(fit, arg) {
  .fa_require_lavaan()
  if (lavaan::lavInspect(fit, "ngroups") > 1L) {
    .dmar_refuse_argument(
      sprintf(paste0(
   "'%s' is a multiple-group lavaan fit; factor_congruence() compares ",
   "two single-group solutions. Pass one group's standardized loadings, ",
   "for example lavaan::lavInspect(fit, \"std\")[[1]]$lambda."), arg),
      arg = arg)
  }
  if (lavaan::lavInspect(fit, "nlevels") > 1L) {
    .dmar_refuse_argument(
      sprintf(paste0(
   "'%s' is a multilevel lavaan fit, whose loadings differ by level; ",
   "factor_congruence() compares two single-level solutions. Pass the ",
   "loading matrix of the level to compare."), arg), arg = arg)
  }
  if (!isTRUE(lavaan::lavInspect(fit, "converged"))) {
    .dmar_refuse(
      sprintf(paste0(
   "The lavaan fit '%s' did not converge, so its estimates are not a ",
   "solution to compare. Resolve the convergence problem and refit."),
   arg), "not_converged")
  }
  lv <- lavaan::lavNames(fit, "lv")
  items <- lavaan::lavNames(fit, "ov.ind")
  if (!length(lv) || !length(items)) {
    .dmar_refuse_argument(
      sprintf("The lavaan fit '%s' has no factor measured by items.", arg),
      arg = NULL)
  }
  lambda <- lavaan::lavInspect(fit, "std")$lambda
  Phi <- lavaan::lavInspect(fit, "cor.lv")
  list(L = unclass(lambda)[items, lv, drop = FALSE],
       Phi = unclass(Phi)[lv, lv, drop = FALSE], source = "lavaan")
}

# The refusals every solution passes through, and the names it is given:
# finite loadings, at least two items, no factor whose loadings are all
# zero, and unique, nonempty item and factor names where they are given.
# Unnamed factors are called F1, F2, ...
.congruence_check <- function(s, arg) {
  L <- s$L
  storage.mode(L) <- "double"
  p <- nrow(L)
  k <- ncol(L)
  if (!k) {
    .dmar_refuse_argument(sprintf("'%s' has no factors.", arg), arg = arg)
  }
  if (p < 2L) {
    .dmar_refuse_argument(
      sprintf(paste0(
   "'%s' has %d item; the coefficient compares columns of loadings and ",
   "needs at least 2 items."), arg, p), arg = arg)
  }
  items <- rownames(L)
  if (!is.null(items) && (anyNA(items) || any(!nzchar(items)) ||
                          anyDuplicated(items))) {
    .dmar_refuse_argument(
      sprintf("The item (row) names of '%s' must be unique and nonempty.",
              arg), arg = NULL)
  }
  if (is.null(colnames(L))) colnames(L) <- paste0("F", seq_len(k))
  factors <- colnames(L)
  if (anyNA(factors) || any(!nzchar(factors)) || anyDuplicated(factors)) {
    .dmar_refuse_argument(
      sprintf("The factor (column) names of '%s' must be unique and nonempty.",
              arg), arg = NULL)
  }
  bad <- !is.finite(L)
  if (any(bad)) {
    rows <- which(rowSums(bad) > 0L)
    where <- if (is.null(items)) paste("rows", paste(rows, collapse = ", "))
             else paste("items", paste(items[rows], collapse = ", "))
    .dmar_refuse_argument(
      sprintf("'%s' has missing or nonfinite loadings (%s).", arg, where),
      arg = arg)
  }
  zero <- colSums(L^2) == 0
  if (any(zero)) {
    several <- sum(zero) > 1L
    .dmar_refuse_argument(
      sprintf(paste0(
   "%s %s of '%s' %s no nonzero loading, so %s congruence with any ",
   "factor is zero divided by zero.%s"),
   if (several) "Factors" else "Factor",
   paste(factors[zero], collapse = ", "), arg,
   if (several) "have" else "has", if (several) "their" else "its",
   if (identical(s$source, "lavaan")) paste0(
     " A higher-order factor, measured by other factors rather than by ",
     "items, is such a factor; to compare the first-order factors, pass ",
     "their columns of lavaan::lavInspect(fit, \"std\")$lambda.") else ""),
      arg = NULL)
  }
  s$L <- L
  if (!is.null(s$Phi)) dimnames(s$Phi) <- list(factors, factors)
  s
}

# The two loading matrices with their items in one order: by name when both
# solutions name their items (the two sets must be the same), otherwise by
# position (the counts must be the same). Returns the two matrices, the
# items (their names, or their positions when neither input names them),
# and how they were aligned.
.congruence_align <- function(sx, sy) {
  nx <- rownames(sx$L)
  ny <- rownames(sy$L)
  if (!is.null(nx) && !is.null(ny)) {
    only_x <- setdiff(nx, ny)
    only_y <- setdiff(ny, nx)
    if (length(only_x) || length(only_y)) {
      .dmar_refuse_argument(
        paste0(
     "'x' and 'y' must hold the same items.",
     if (length(only_x)) paste0(" In 'x' and not in 'y': ",
                                paste(only_x, collapse = ", "), "."),
     if (length(only_y)) paste0(" In 'y' and not in 'x': ",
                                paste(only_y, collapse = ", "), ".")),
        arg = NULL)
    }
    return(list(Lx = sx$L, Ly = sy$L[nx, , drop = FALSE], items = nx,
                by = "name"))
  }
  if (nrow(sx$L) != nrow(sy$L)) {
    .dmar_refuse_argument(
      sprintf(paste0(
   "'x' has %d items and 'y' has %d. Without item names on both, items ",
   "are matched by position, which needs the same number in each."),
   nrow(sx$L), nrow(sy$L)), arg = c("x", "y"))
  }
  items <- if (!is.null(nx)) nx else if (!is.null(ny)) ny else seq_len(nrow(sx$L))
  Lx <- sx$L
  Ly <- sy$L
  rownames(Lx) <- rownames(Ly) <- if (is.character(items)) items
  list(Lx = Lx, Ly = Ly, items = items, by = "position")
}

# Tucker's congruence coefficient between every column of Lx and every
# column of Ly (Lorenzo-Seva & ten Berge, 2006, eq. 1): the cross products
# divided by the square roots of the column sums of squares. Rows are Lx's
# columns, columns Ly's.
.congruence_matrix <- function(Lx, Ly) {
  crossprod(Lx, Ly) / sqrt(outer(colSums(Lx^2), colSums(Ly^2)))
}

# The one-to-one assignment of the columns of `score` to its rows that
# maximizes the sum of the assigned entries. Returns, for each row, the
# column assigned to it, or NA for a row left unassigned because there are
# more rows than columns (a column left over when there are more columns
# is simply not named). A rectangular matrix is padded with zeros to a
# square one, whose extra rows or columns stand for "unassigned" and
# contribute nothing to any total.
#
# Ties are broken in favor of the listed order: among the assignments
# whose totals reach the maximum (to within a relative sqrt(machine
# epsilon), so that rounding cannot decide), the first row takes the
# earliest column it can, then the second row, and so on, an unassigned
# row counting as later than every column. Each choice is tested by
# solving what remains, so the result is the lexicographically first
# optimal assignment, whatever order the solver itself visits.
.assignment_max <- function(score) {
  if (!is.matrix(score) || !is.numeric(score) || any(!is.finite(score))) {
    .dmar_refuse_argument(
      "'score' must be a finite numeric matrix.", arg = "score")
  }
  nr <- nrow(score)
  nc <- ncol(score)
  n <- max(nr, nc)
  if (!nr) return(integer(0))
  if (!nc) return(rep(NA_integer_, nr))
  S <- matrix(0, n, n)
  S[seq_len(nr), seq_len(nc)] <- score
  best <- .assignment_total(S)
  tol <- sqrt(.Machine$double.eps) * max(1, abs(best))
  assigned <- integer(n)
  rows <- seq_len(n)
  cols <- seq_len(n)
  fixed <- 0
  for (i in seq_len(nr)) {
    rest_rows <- rows[rows != i]
    chosen <- NA_integer_
    top <- -Inf
    top_j <- NA_integer_
    for (j in cols) {
      rest_cols <- cols[cols != j]
      total <- fixed + S[i, j] +
        .assignment_total(S[rest_rows, rest_cols, drop = FALSE])
      if (total >= best - tol) {
        chosen <- j
        break
      }
      if (total > top) {
        top <- total
        top_j <- j
      }
    }
    # Every row can reach the maximum through the optimal assignment
    # itself; the fallback only guards against rounding.
    if (is.na(chosen)) chosen <- top_j
    assigned[i] <- chosen
    fixed <- fixed + S[i, chosen]
    rows <- rest_rows
    cols <- cols[cols != chosen]
  }
  out <- assigned[seq_len(nr)]
  out[out > nc] <- NA_integer_
  out
}

# The largest total of a square matrix over one-to-one assignments.
.assignment_total <- function(S) {
  n <- nrow(S)
  if (!n) return(0)
  sum(S[cbind(seq_len(n), .assignment_hungarian(S))])
}

# The Hungarian method (Kuhn, 1955) for a square matrix, maximizing the
# total: the shortest augmenting path form with row and column potentials,
# which adds one row at a time to the assignment in O(n^2) steps, O(n^3)
# in all. The costs are max(S) - S, so that a minimum cost assignment is a
# maximum total one. Column 0 is a virtual column that holds the row being
# added; the vectors indexed by column carry it at position 1, so column j
# is at position j + 1. Returns, for each row, its column.
.assignment_hungarian <- function(S) {
  n <- nrow(S)
  cost <- max(S) - S
  u <- numeric(n)                # row potentials
  v <- numeric(n + 1L)           # column potentials, column j at j + 1
  p <- integer(n + 1L)           # the row assigned to column j, 0 if none
  way <- integer(n + 1L)         # the previous column on the shortest path
  for (i in seq_len(n)) {
    p[1L] <- i
    j0 <- 0L
    minv <- rep(Inf, n + 1L)
    used <- rep(FALSE, n + 1L)
    repeat {
      used[j0 + 1L] <- TRUE
      i0 <- p[j0 + 1L]
      free <- which(!used[-1L])
      reduced <- cost[i0, free] - u[i0] - v[free + 1L]
      lower <- reduced < minv[free + 1L]
      minv[free[lower] + 1L] <- reduced[lower]
      way[free[lower] + 1L] <- j0
      k <- which.min(minv[free + 1L])
      delta <- minv[free[k] + 1L]
      j1 <- free[k]
      on_path <- which(used)
      u[p[on_path]] <- u[p[on_path]] + delta
      v[on_path] <- v[on_path] - delta
      minv[!used] <- minv[!used] - delta
      j0 <- j1
      if (p[j0 + 1L] == 0L) break
    }
    repeat {
      j1 <- way[j0 + 1L]
      p[j0 + 1L] <- p[j1 + 1L]
      j0 <- j1
      if (j0 == 0L) break
    }
  }
  out <- integer(n)
  out[p[-1L]] <- seq_len(n)
  out
}
