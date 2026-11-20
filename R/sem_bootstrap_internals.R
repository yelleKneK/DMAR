# Internals of sem_exact_fit_bootstrap(): the model based (Bollen and
# Stine, 1992) transformation of complete data, its extensions to
# incomplete data (Savalei & Yuan, 2009), the refit of one bootstrap
# sample, and the tally of the replications. None is exported.

# ---- symmetric square roots ---------------------------------------------

# The symmetric square root of a positive definite matrix, and of its
# inverse, from the eigen decomposition (as lavaan forms them for its own
# Bollen-Stine bootstrap). A matrix that is not positive definite is
# refused with the name of the moment it came from (.bs_refuse_undefined(),
# reason "not_positive_definite"); `group` is the index of the group whose
# moment it is.
.bs_sym_sqrt <- function(M, what, inverse = FALSE, group = 1L) {
  e <- eigen((M + t(M)) / 2, symmetric = TRUE)
  if (any(!is.finite(e$values)) || min(e$values) <= 0) {
    .bs_refuse_undefined(c("The ", what, " is not positive definite, so the ",
                           "model based transformation is undefined."),
                         "not_positive_definite", group, matrix = what)
  }
  d <- if (inverse) 1 / sqrt(e$values) else sqrt(e$values)
  e$vectors %*% (d * t(e$vectors))
}

# The refusal of a transformation that is undefined for a group's data or
# moments: an error of class dmar_transformation_undefined (a
# dmar_refusal; .dmar_refuse(), R/dmar_conditions.R) carrying `reason`,
# the check that failed ("not_positive_definite", "no_mean_structure",
# "nonpositive_variance", "singular_pattern", "not_converged",
# "singular_jacobian", or "diverged"), `group`, the index of the group,
# and, through `...`, what the reason concerns (the `matrix` that is not
# positive definite, the `pattern` whose implied covariance matrix is
# singular, the Newton `step`).
.bs_refuse_undefined <- function(message, reason, group, ...) {
  .dmar_refuse(message, "transformation_undefined", reason = reason,
               group = group, ...)
}

# ---- the transformation of one group ------------------------------------

# Transforms the analysis data X of one group (cases by observed variables,
# NA for a missing value) so that the moments a maximum likelihood fit of
# the saturated model estimates from the transformed data are the moments
# the hypothesized model implies, Sigma (and Mu under a mean structure).
# Bootstrap samples drawn from the transformed data are then samples from a
# population in which the model holds exactly, as estimated, which is what
# the model based bootstrap of Bollen and Stine (1992) needs.
#
#   X      : numeric matrix of the group's analysis data, columns in the
#            order of Sigma's dimnames; rows with every value missing are
#            removed.
#   Sigma  : the model implied covariance matrix of the group.
#   Mu     : the model implied means, or NULL without a mean structure.
#   S      : the saturated (H1) covariance matrix lavaan estimated for the
#            group: the sample covariance with lavaan's divisor for complete
#            data, and the EM estimate under FIML. It is the S of the
#            complete data transformation; with incomplete data it supplies
#            the scale of the variables and the starting value of the
#            common transformation's iteration, and the pattern
#            transformation does not use it.
#   method : "common" or "pattern" (see sem_exact_fit_bootstrap()); with
#            complete data both reduce to the Bollen-Stine transformation.
#   group  : the index of the group, which a refusal carries.
#   group_label : the group's label, named in a refusal's message; NULL
#            for a fit of one group.
#
# With incomplete data both transformations are solved on the variables
# divided by their saturated standard deviations, and the transformed
# values are multiplied back. The solution therefore does not depend on
# the variables' units (rescaling an item rescales its transformed values
# and leaves the others as they are), and the estimating equations and the
# convergence tolerance mean the same whatever the units. It is not the
# solution of the same transformation in raw units: a symmetric square
# root, and the symmetric matrix of the common transformation, do not
# commute with a rescaling of the variables (for a diagonal rescaling L,
# L (L S L)^{-1/2} (L Sigma L)^{1/2} L^{-1} is not S^{-1/2} Sigma^{1/2}),
# so the raw unit version of either transformation (semTools' bsBootMiss()
# transformation 2 is the common one in raw units) gives other transformed
# values, element by element, and depends on the units. Both versions make
# the saturated FIML estimates the model implied moments exactly, so the
# model holds exactly in either. The complete data transformation is
# solved in raw units, as lavaan's bollen.stine bootstrap solves it, so
# that the two agree resample for resample; it depends on the units as
# lavaan's does.
#
# Returns a list: Y (the transformed data, missing values where X has
# them), method (the transformation applied: "bollen_stine", "common", or
# "pattern"), patterns (how many missing data patterns), iterations and
# score (the common transformation's Newton steps and the largest absolute
# element of the saturated estimating equations at the solution, on the
# standardized scale; the pattern transformation records the same check
# computed afresh, and complete data record 0 and NA), and mean_residual,
# the largest absolute difference between a pattern's transformed means and
# Mu on the original scale (NA without a mean structure).
#
# A transformation that is undefined for the data or the moments is refused
# with .bs_refuse_undefined(), and a pattern too small for the pattern-wise
# transformation with a refusal of class dmar_pattern_too_small (a
# dmar_refusal) carrying the pattern's index among the group's patterns
# (`pattern`, numbered in the order of the pattern's first row), its rows
# (`n`), its observed variables (`p`, and their names, `observed`), and the
# group's index (`group`).
.bs_transform_group <- function(X, Sigma, Mu, S, method = c("common", "pattern"),
                                group = 1L, group_label = NULL) {
  method <- match.arg(method)
  X <- as.matrix(X)
  X <- X[rowSums(!is.na(X)) > 0L, , drop = FALSE]
  p <- ncol(X)
  in_group <- if (is.null(group_label)) "" else paste0(" of group ", group_label)
  if (!anyNA(X)) {
    # Complete data: Bollen and Stine's (1992) transformation, as Enders
    # (2002, Eq. 1) writes it, Z = (Y - Ybar) S^{-1/2} Sigma^{1/2}, with
    # symmetric square roots and the mean structure added back (lavaan's
    # own bollen.stine bootstrap forms it the same way).
    Xc <- scale(X, center = TRUE, scale = FALSE)
    Y <- Xc %*% .bs_sym_sqrt(S, "saturated covariance matrix", inverse = TRUE,
                             group = group) %*%
      .bs_sym_sqrt(Sigma, "model implied covariance matrix", group = group)
    if (!is.null(Mu)) Y <- sweep(Y, 2L, as.numeric(Mu), `+`)
    dimnames(Y) <- dimnames(X)
    return(list(Y = Y, method = "bollen_stine", patterns = 1L, iterations = 0L,
                score = NA_real_,
                mean_residual = if (is.null(Mu)) NA_real_ else
                  max(abs(colMeans(Y) - as.numeric(Mu)))))
  }
  if (is.null(Mu)) {
    .bs_refuse_undefined(c("A fit of incomplete data needs a mean structure ",
                           "for the model based transformation (lavaan's ",
                           "FIML fits always have one)."),
                         "no_mean_structure", group)
  }
  # Standardize by the saturated standard deviations (see above).
  s <- sqrt(diag(S))
  if (any(!is.finite(s)) || any(s <= 0)) {
    .bs_refuse_undefined(c("The saturated covariance matrix has a variance ",
                           "that is not positive, so the model based ",
                           "transformation is undefined."),
                         "nonpositive_variance", group)
  }
  Mu0 <- as.numeric(Mu)
  X <- sweep(X, 2L, s, `/`)
  Sigma <- Sigma / tcrossprod(s)
  Mu <- Mu0 / s
  S <- S / tcrossprod(s)
  observed <- !is.na(X)
  pattern <- apply(observed, 1L, function(r) paste(as.integer(r), collapse = ""))
  patterns <- unique(pattern)
  J <- length(patterns)
  rows <- lapply(patterns, function(k) which(pattern == k))
  # For each pattern: its rows centered at their own means, the observed
  # variables, and the cross product matrix of the centered rows (padded
  # with zeros to p by p), H_j, and the inverse of the model implied
  # covariance matrix of its observed variables, padded the same way, M_j.
  Xc <- vector("list", J); H <- vector("list", J); M <- vector("list", J)
  obs <- vector("list", J); n_j <- integer(J)
  for (j in seq_len(J)) {
    r <- rows[[j]]
    o <- observed[r[1L], ]
    obs[[j]] <- o
    n_j[j] <- length(r)
    xj <- X[r, o, drop = FALSE]
    xc <- scale(xj, center = TRUE, scale = FALSE)
    Xc[[j]] <- xc
    Hj <- matrix(0, p, p); Hj[o, o] <- crossprod(xc) / n_j[j]
    H[[j]] <- Hj
    Mj <- matrix(0, p, p)
    Mj[o, o] <- tryCatch(solve(Sigma[o, o, drop = FALSE]), error = function(e) {
      .bs_refuse_undefined(c("The model implied covariance matrix of the ",
                             "variables observed in one missing data pattern ",
                             "is singular, so the transformation is ",
                             "undefined."),
                           "singular_pattern", group, pattern = j)
    })
    M[[j]] <- Mj
  }
  Y <- X
  finish <- function(Y, method, iterations, score) {
    Y <- sweep(Y, 2L, s, `*`)
    dimnames(Y) <- dimnames(X)
    list(Y = Y, method = method, patterns = J, iterations = iterations,
         score = score, mean_residual = .bs_mean_residual(Y, Mu0, pattern))
  }
  if (method == "pattern") {
    # Each pattern transformed on its own, with its own covariance matrix
    # and means and the submatrix of Sigma of its observed variables
    # (Savalei and Yuan's first transformation): its rows then have
    # covariance matrix exactly the model implied one and means exactly
    # Mu. The pattern's own covariance matrix must be positive definite,
    # which needs more rows than observed variables.
    for (j in seq_len(J)) {
      o <- obs[[j]]
      if (n_j[j] <= sum(o)) {
        .dmar_refuse(c("The pattern-wise transformation needs every missing ",
                       "data pattern to have more rows than observed ",
                       "variables; pattern ", j, in_group, " has ", n_j[j],
                       " rows on ", sum(o), " variables. Use ",
                       "transformation = \"common\"."),
                     "pattern_too_small", pattern = j, n = n_j[j],
                     p = sum(o), observed = colnames(X)[o], group = group)
      }
      Sj <- crossprod(Xc[[j]]) / n_j[j]
      Yj <- Xc[[j]] %*% .bs_sym_sqrt(Sj, paste0("covariance matrix of pattern ", j),
                                     inverse = TRUE, group = group) %*%
        .bs_sym_sqrt(Sigma[o, o, drop = FALSE], "model implied covariance matrix",
                     group = group)
      Y[rows[[j]], o] <- sweep(Yj, 2L, Mu[o], `+`)
    }
    return(finish(Y, "pattern", 0L, .bs_saturated_score(Y, Sigma, Mu, observed)))
  }
  # The common transformation (Savalei and Yuan's second): one symmetric
  # matrix A for every pattern, Y_j = Xc_j A[o, o] + Mu[o], chosen so that
  # the saturated estimating equations for the covariances of a FIML fit,
  # summed over the patterns, hold at Sigma:
  #   sum_j n_j D^+ vec(M_j A H_j A M_j - M_j) = 0,
  # with D the duplication matrix (vech(A) are the unknowns). Newton's
  # method from A = sym(Sigma^{1/2} S^{-1/2}), the complete data
  # transformation applied to the EM moments; the Jacobian of the left side
  # with respect to vech(A) is 2 sum_j n_j D^+ (M_j A H_j kronecker M_j) D.
  # Each pattern is centered at its own means and shifted to Mu, so the
  # equations for the means hold by construction.
  D <- lavaan::lav_matrix_duplication(p)
  Dplus <- solve(crossprod(D), t(D))
  score <- function(A) {
    g <- 0
    for (j in seq_len(J)) {
      Tj <- M[[j]] %*% A %*% H[[j]] %*% A %*% M[[j]] - M[[j]]
      g <- g + n_j[j] * (Dplus %*% as.vector(Tj))
    }
    as.vector(g)
  }
  jacobian <- function(A) {
    Jm <- 0
    for (j in seq_len(J)) {
      Jm <- Jm + 2 * n_j[j] * (Dplus %*% kronecker(M[[j]] %*% A %*% H[[j]], M[[j]]) %*% D)
    }
    Jm
  }
  A <- .bs_sym_sqrt(Sigma, "model implied covariance matrix", group = group) %*%
    .bs_sym_sqrt(S, "saturated covariance matrix", inverse = TRUE, group = group)
  A <- (A + t(A)) / 2
  g <- score(A)
  it <- 0L
  # On the standardized scale the equations are of order n_j, so the
  # tolerance is relative to the number of rows.
  tol <- 1e-10 * sum(n_j)
  while (max(abs(g)) > tol) {
    it <- it + 1L
    if (it > 200L) {
      .bs_refuse_undefined(c("The common transformation did not converge in ",
                             "200 Newton steps (largest estimating equation ",
                             format(max(abs(g)), digits = 3), "). Use ",
                             "transformation = \"pattern\" when every missing ",
                             "data pattern has many rows, or check the fit's ",
                             "moments."),
                           "not_converged", group, step = it - 1L)
    }
    step <- tryCatch(solve(jacobian(A), g), error = function(e) {
      .bs_refuse_undefined(c("The Jacobian of the common transformation is ",
                             "singular at step ", it, "; the transformation is ",
                             "undefined for these data."),
                           "singular_jacobian", group, step = it)
    })
    a <- as.vector(D %*% (Dplus %*% as.vector(A) - step))
    A <- matrix(a, p, p)
    A <- (A + t(A)) / 2
    g <- score(A)
    if (any(!is.finite(g))) {
      .bs_refuse_undefined(c("The common transformation diverged; the ",
                             "transformation is undefined for these data."),
                           "diverged", group, step = it)
    }
  }
  for (j in seq_len(J)) {
    o <- obs[[j]]
    Yj <- Xc[[j]] %*% A[o, o, drop = FALSE]
    Y[rows[[j]], o] <- sweep(Yj, 2L, Mu[o], `+`)
  }
  finish(Y, "common", it, max(abs(g)))
}

# The largest absolute element of the saturated FIML estimating equations
# for the covariances at (Mu, Sigma) on data Y (the quantity the common
# transformation sets to zero), computed afresh from Y: a check that holds
# for every transformation.
.bs_saturated_score <- function(Y, Sigma, Mu, observed) {
  p <- ncol(Y)
  pattern <- apply(observed, 1L, function(r) paste(as.integer(r), collapse = ""))
  D <- lavaan::lav_matrix_duplication(p)
  Dplus <- solve(crossprod(D), t(D))
  g <- 0
  for (k in unique(pattern)) {
    r <- which(pattern == k)
    o <- observed[r[1L], ]
    Yc <- sweep(Y[r, o, drop = FALSE], 2L, Mu[o], `-`)
    Sj <- matrix(0, p, p); Sj[o, o] <- crossprod(Yc) / length(r)
    Mj <- matrix(0, p, p); Mj[o, o] <- solve(Sigma[o, o, drop = FALSE])
    g <- g + length(r) * (Dplus %*% as.vector(Mj %*% Sj %*% Mj - Mj))
  }
  max(abs(g))
}

# The largest absolute difference between the means of any pattern's
# transformed rows and Mu on the variables it observes.
.bs_mean_residual <- function(Y, Mu, pattern) {
  worst <- 0
  for (k in unique(pattern)) {
    r <- which(pattern == k)
    o <- !is.na(Y[r[1L], ])
    worst <- max(worst, abs(colMeans(Y[r, o, drop = FALSE]) - Mu[o]))
  }
  worst
}

# ---- one bootstrap replication ------------------------------------------

# Fits the model to one bootstrap sample and returns its standard chi
# square statistic with the fate of the replication. `rows` is a list, one
# integer vector per group, of the rows of the transformed data to draw;
# `Y` the list of transformed group data; `group_var` and `group_labels`
# the group column and its labels (NULL and NULL for one group); `pt`,
# `model`, and `options` the slots of the original fit, with the standard
# errors and extra tests turned off. The statistic of a replication that
# did not converge, raised an error, or returned no finite statistic is
# NA, and its status says which; a converged fit whose solution is
# inadmissible is "improper". Admissibility is judged as for the fit
# itself (.sem_admissible(), R/sem_fit_internals.R): by lavaan's
# post-estimation check, and for a fit with auxiliary variables, which the
# refit inherits through the fit's options, by the same check of the
# model's own parameters. Judged by lavaan's check, the resamples of a fit
# whose auxiliary variable is related to a factor are flagged when only
# the auxiliary variable's residual covariances are indefinite: on six of
# the holzinger_swineford tests with age as the auxiliary variable, 125 of
# 300 resamples, of which 16 had an improper solution of the model.
#
# A replication whose likelihood has no maximum, an error variance running
# off below zero along a ridge with the likelihood flat or still rising
# (.lavaan_no_maximum(), R/no_maximum_internals.R), is improper like any
# other with a negative variance, and is not asked about. The bootstraps of
# the reliability family set such a replication aside because the
# coefficient they resample is read from the estimates, which are wherever
# lavaan's optimizer stopped on the ridge; the statistic here is the chi
# square, and the likelihood levels off along the ridge, where the
# optimizer stops only once its steps no longer gain, so the statistic is
# within the optimizer's tolerance of its limit wherever it stops. Setting
# the replication aside would condition the reference distribution on the
# samples without a ridge, and override `improper`.
.bs_refit_one <- function(rows, Y, group_var, group_labels, pt, model, options) {
  parts <- lapply(seq_along(Y), function(g) {
    d <- as.data.frame(Y[[g]][rows[[g]], , drop = FALSE])
    if (!is.null(group_var)) d[[group_var]] <- group_labels[g]
    d
  })
  d <- do.call(rbind, parts)
  rownames(d) <- NULL
  fit_b <- tryCatch(
    suppressWarnings(lavaan::lavaan(data = d, group = group_var,
                                    slotParTable = pt, slotModel = model,
                                    slotOptions = options)),
    error = function(e) NULL)
  if (is.null(fit_b)) return(c(statistic = NA_real_, status = 4))
  if (!isTRUE(lavaan::lavTech(fit_b, "converged"))) {
    return(c(statistic = NA_real_, status = 3))
  }
  stat <- tryCatch(lavaan::lavTech(fit_b, "test")[["standard"]]$stat,
                   error = function(e) NA_real_)
  if (!is.numeric(stat) || length(stat) != 1L || !is.finite(stat)) {
    return(c(statistic = NA_real_, status = 4))
  }
  proper <- isTRUE(tryCatch(suppressWarnings(.sem_admissible(fit_b, explain = FALSE)),
                            error = function(e) FALSE))
  c(statistic = as.numeric(stat), status = if (proper) 1 else 2)
}

# The status codes of .bs_refit_one(), as a factor.
.BS_STATUS <- c("proper", "improper", "not_converged", "error")

# Runs the replications in this session or over `ncpus` processes, with
# every index set drawn beforehand so that a seed gives the same result
# whatever ncpus is. Socket workers (Windows) load the installed DMAR and
# lavaan, since the refit closure refers to DMAR's internals by name.
# Returns the matrix of statistics and status codes.
.bs_run_replications <- function(index_sets, refit, ncpus) {
  par_type <- .boot_parallel_type(ncpus)
  results <- if (par_type == "no") {
    lapply(index_sets, refit)
  } else if (par_type == "multicore") {
    parallel::mclapply(index_sets, refit, mc.cores = ncpus)
  } else {
    cl <- parallel::makePSOCKcluster(ncpus)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterEvalQ(cl, requireNamespace("lavaan", quietly = TRUE))
    parallel::parLapply(cl, index_sets, refit)
  }
  out <- do.call(rbind, lapply(results, function(r) {
    if (!is.numeric(r) || length(r) != 2L) c(statistic = NA_real_, status = 4) else r
  }))
  colnames(out) <- c("statistic", "status")
  out
}
