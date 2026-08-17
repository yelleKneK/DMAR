# Shared engine for the nonlinear growth-curve simulators:
# simulate_longitudinal_negative_exponential(), _logistic(), _gompertz(),
# and _richards(). Each public function validates its own parameter
# vocabulary and calls .simulate_longitudinal_growth(), which draws the
# per-unit random coefficients, evaluates the curve, and applies the
# same level-one error and assessment-time machinery as
# simulate_longitudinal_polynomial() (whose .error_cor_matrix() it
# reuses). The parameterizations are those of Kelley (2005, dissertation;
# 2008, Methodology): each curve carries the intercept-shifting zeta so
# the lower asymptote is a modeled quantity rather than fixed at zero.

# Curve registry: mean function, gradient with respect to the parameter
# vector (used for the delta method true-score variance behind the
# reliability solve), and the parameter names in canonical order.
.growth_curves <- list(
  negative_exponential = list(
    par_names = c("alpha", "zeta", "gamma"),
    mean = function(t, p) p[["alpha"]] + p[["zeta"]] * exp(-p[["gamma"]] * t),
    grad = function(t, p) {
      e <- exp(-p[["gamma"]] * t)
      cbind(alpha = rep(1, length(t)), zeta = e,
            gamma = -p[["zeta"]] * t * e)
    }
  ),
  logistic = list(
    par_names = c("alpha", "beta", "gamma", "zeta"),
    mean = function(t, p) {
      E <- exp(-p[["gamma"]] * (t - p[["beta"]]))
      p[["alpha"]] / (1 + E) + p[["zeta"]]
    },
    grad = function(t, p) {
      E <- exp(-p[["gamma"]] * (t - p[["beta"]]))
      d <- (1 + E)^2
      cbind(alpha = 1 / (1 + E),
            beta  = -p[["alpha"]] * p[["gamma"]] * E / d,
            gamma = p[["alpha"]] * (t - p[["beta"]]) * E / d,
            zeta  = rep(1, length(t)))
    }
  ),
  gompertz = list(
    par_names = c("alpha", "beta", "gamma", "zeta"),
    mean = function(t, p) {
      W <- exp(-p[["gamma"]] * (t - p[["beta"]]))
      p[["alpha"]] * exp(-W) + p[["zeta"]]
    },
    grad = function(t, p) {
      W <- exp(-p[["gamma"]] * (t - p[["beta"]]))
      eW <- exp(-W)
      cbind(alpha = eW,
            beta  = -p[["alpha"]] * p[["gamma"]] * W * eW,
            gamma = p[["alpha"]] * (t - p[["beta"]]) * W * eW,
            zeta  = rep(1, length(t)))
    }
  ),
  richards = list(
    par_names = c("alpha", "beta", "gamma", "delta", "zeta"),
    mean = function(t, p) {
      E <- exp(-p[["gamma"]] * (t - p[["beta"]]))
      Q <- 1 + p[["delta"]] * E
      p[["alpha"]] * Q^(-1 / p[["delta"]]) + p[["zeta"]]
    },
    grad = function(t, p) {
      dl <- p[["delta"]]
      E <- exp(-p[["gamma"]] * (t - p[["beta"]]))
      Q <- 1 + dl * E
      core <- Q^(-1 / dl)
      cbind(alpha = core,
            beta  = -p[["alpha"]] * p[["gamma"]] * E * Q^(-1 / dl - 1),
            gamma = p[["alpha"]] * (t - p[["beta"]]) * E * Q^(-1 / dl - 1),
            delta = p[["alpha"]] * core *
              (log(Q) / dl^2 - E / (dl * Q)),
            zeta  = rep(1, length(t)))
    }
  )
)

# Validate one population's fixed-parameter vector against the curve's
# vocabulary: numeric, complete, and named exactly (order is normalized
# here so the public functions can accept any order).
.growth_check_params <- function(params, par_names, model) {
  if (!is.numeric(params) || anyNA(params) ||
      length(params) != length(par_names)) {
    stop(sprintf(paste0("Each parameter vector must be numeric, complete, ",
                        "and of length %d: %s."),
                 length(par_names),
                 paste(par_names, collapse = ", ")), call. = FALSE)
  }
  if (is.null(names(params)) || any(names(params) == "")) {
    names(params) <- par_names
  }
  if (!setequal(names(params), par_names)) {
    stop(sprintf("Parameter names must be %s (got %s).",
                 paste(par_names, collapse = ", "),
                 paste(names(params), collapse = ", ")), call. = FALSE)
  }
  params[par_names]
}

.simulate_longitudinal_growth <- function(model,
                                          n,
                                          target_times = NULL,
                                          fixed_parameters,
                                          time_range = NULL,
                                          occasions = NULL,
                                          time_distribution = "uniform",
                                          random_variances = 0,
                                          random_correlation = NULL,
                                          error_variance = NULL,
                                          reliability = NULL,
                                          error_structure = c("independent",
                                                              "ar1",
                                                              "compound_symmetry",
                                                              "toeplitz"),
                                          error_correlation = NULL,
                                          timing_sd = 0) {
  error_structure <- match.arg(error_structure)
  curve <- .growth_curves[[model]]
  par_names <- curve$par_names
  K <- length(par_names)

  # ---------- Populations (one fixed-parameter vector each) ----------
  if (is.list(fixed_parameters)) {
    if (length(fixed_parameters) < 1L) {
      stop("'fixed_parameters' must contain at least one parameter vector.",
           call. = FALSE)
    }
    coef_list <- lapply(fixed_parameters, .growth_check_params,
                        par_names = par_names, model = model)
  } else {
    coef_list <- list(.growth_check_params(fixed_parameters,
                                           par_names = par_names,
                                           model = model))
  }
  G <- length(coef_list)

  # ---------- Measurement schedule: shared or unit-specific ----------
  if (is.null(target_times) == is.null(time_range)) {
    stop("Specify exactly one of 'target_times' (one shared schedule for ",
         "every unit) or 'time_range' (unit-specific times drawn between ",
         "a lower and an upper bound).", call. = FALSE)
  }
  individual_times <- !is.null(time_range)
  if (individual_times) {
    time_distribution <- match.arg(time_distribution, c("uniform"))
    if (!is.numeric(time_range) || length(time_range) != 2L ||
        any(!is.finite(time_range)) || time_range[1L] >= time_range[2L]) {
      stop("'time_range' must be c(lower, upper) with lower < upper.",
           call. = FALSE)
    }
    if (is.null(occasions) || !is.numeric(occasions) || anyNA(occasions) ||
        !length(occasions) %in% c(1L, 2L) || any(occasions < 1) ||
        any(occasions != round(occasions)) ||
        (length(occasions) == 2L && occasions[1L] > occasions[2L])) {
      stop("With 'time_range', give 'occasions' as a single positive ",
           "integer (the same number of measurement times for every unit) ",
           "or c(min, max) for a unit-varying number.", call. = FALSE)
    }
    occasions <- as.integer(occasions)
    M <- max(occasions)
  } else {
    if (!is.numeric(target_times) || length(target_times) < 1L ||
        any(!is.finite(target_times))) {
      stop("'target_times' must be a numeric vector of finite measurement ",
           "times.", call. = FALSE)
    }
    target_times <- as.numeric(target_times)
    M <- length(target_times)
  }

  # ---------- Per-population unit counts ----------
  if (length(n) == 1L) {
    if (!is.numeric(n) || is.na(n) || n < 1 || n != round(n)) {
      stop("'n' must be a single positive integer (units per population).",
           call. = FALSE)
    }
    n_per_population <- rep(as.integer(n), G)
  } else if (length(n) == G) {
    if (!is.numeric(n) || anyNA(n) || any(n < 1) || any(n != round(n))) {
      stop("Each entry of 'n' must be a positive integer.", call. = FALSE)
    }
    n_per_population <- as.integer(n)
  } else {
    stop(sprintf(paste0("'n' must be a single positive integer or a numeric ",
                        "vector of length %d (one count of units per population)."),
                 G), call. = FALSE)
  }

  # ---------- Between-unit covariance of the curve parameters ----------
  if (!is.numeric(random_variances) || anyNA(random_variances) ||
      any(random_variances < 0)) {
    stop("'random_variances' must be non-negative.", call. = FALSE)
  }
  if (!is.null(names(random_variances)) && any(names(random_variances) != "")) {
    # A named vector may name any subset of the parameters; the rest get
    # variance zero. This makes single-parameter heterogeneity direct:
    # random_variances = c(beta = 0.8) varies only the inflection time.
    bad <- setdiff(names(random_variances), par_names)
    if (length(bad)) {
      stop(sprintf("Unknown parameter name(s) in 'random_variances': %s. ",
                   paste(bad, collapse = ", ")),
           sprintf("The parameters are %s.",
                   paste(par_names, collapse = ", ")), call. = FALSE)
    }
    full <- stats::setNames(rep(0, K), par_names)
    full[names(random_variances)] <- as.numeric(random_variances)
    random_variances <- full
  } else if (length(random_variances) == 1L) {
    random_variances <- rep(as.numeric(random_variances), K)
  } else if (length(random_variances) != K) {
    stop(sprintf(paste0("'random_variances' must be a single number, a ",
                        "vector of length %d (one variance per parameter, ",
                        "in the order %s), or a named vector over any ",
                        "subset of those parameters."),
                 K, paste(par_names, collapse = ", ")), call. = FALSE)
  }
  names(random_variances) <- par_names

  if (is.null(random_correlation)) {
    Rmat <- diag(K)
  } else {
    Rmat <- as.matrix(random_correlation)
    if (nrow(Rmat) != K || ncol(Rmat) != K) {
      stop(sprintf("'random_correlation' must be a %d-by-%d matrix.", K, K),
           call. = FALSE)
    }
    if (any(abs(Rmat - t(Rmat)) > 1e-8)) {
      stop("'random_correlation' must be symmetric.", call. = FALSE)
    }
    if (any(abs(diag(Rmat) - 1) > 1e-8)) {
      stop("'random_correlation' must have 1 on the diagonal (it is a ",
           "correlation matrix).", call. = FALSE)
    }
    ev <- eigen(Rmat, symmetric = TRUE, only.values = TRUE)$values
    if (min(ev) < -1e-8) {
      stop("'random_correlation' must be positive semidefinite.", call. = FALSE)
    }
  }
  std_dev <- sqrt(random_variances)
  Tcov <- outer(std_dev, std_dev) * Rmat
  dimnames(Tcov) <- list(par_names, par_names)

  # ---------- Assessment-time jitter ----------
  if (!is.numeric(timing_sd) || anyNA(timing_sd) || any(timing_sd < 0)) {
    stop("'timing_sd' must be non-negative.", call. = FALSE)
  }
  if (individual_times) {
    if (!(length(timing_sd) == 1L && timing_sd == 0)) {
      stop("'timing_sd' does not combine with 'time_range': unit-specific ",
           "times are already drawn per unit, so jitter around a nominal ",
           "schedule has no meaning there.", call. = FALSE)
    }
  } else if (length(timing_sd) == 1L) {
    timing_sd <- rep(as.numeric(timing_sd), M)
  } else if (length(timing_sd) != M) {
    stop(sprintf(paste0("'timing_sd' must have length 1 or %d (one standard ",
                        "deviation per occasion)."), M), call. = FALSE)
  }

  # ---------- Level-one error covariance ----------
  if (is.null(error_variance) == is.null(reliability)) {
    stop("Specify exactly one of 'error_variance' or 'reliability'.",
         call. = FALSE)
  }
  if (individual_times) {
    # With unit-specific times there is no shared occasion grid, so the
    # per-occasion machinery (reliability targets, heteroscedastic or
    # correlated errors) has nothing to attach to.
    if (!is.null(reliability)) {
      stop("'reliability' needs a shared measurement schedule: ",
           "per-occasion reliability is undefined when every unit has its ",
           "own times. Give 'error_variance' instead.", call. = FALSE)
    }
    if (!is.numeric(error_variance) || length(error_variance) != 1L ||
        is.na(error_variance) || error_variance < 0) {
      stop("With 'time_range', 'error_variance' must be a single ",
           "non-negative number: per-occasion vectors and covariance ",
           "matrices need a shared schedule.", call. = FALSE)
    }
    if (error_structure != "independent" || !is.null(error_correlation)) {
      stop("With 'time_range', 'error_structure' must be \"independent\": ",
           "across-occasion error correlation needs a shared schedule.",
           call. = FALSE)
    }
    err_scalar <- as.numeric(error_variance)
    Sigma_e <- NULL
    err_var_marginal <- err_scalar
    reliability_by_occ <- NA_real_
  } else {
  # Delta-method (first-order) true-score variance at each nominal
  # occasion, g(t)' T g(t) with g the gradient at the population's fixed
  # parameters, averaged over populations with weights n_per_population. Exact for
  # the polynomial simulator's linear case; here it is the first-order
  # approximation the reliability target is defined against.
  tv_by_population <- vapply(coef_list, function(p) {
    gmat <- curve$grad(target_times, as.list(p))
    rowSums((gmat %*% Tcov) * gmat)
  }, numeric(M))
  tv_by_population <- matrix(tv_by_population, nrow = M)
  true_var_by_occ <- as.numeric(tv_by_population %*% n_per_population) / sum(n_per_population)

  if (!is.null(error_variance) && is.matrix(error_variance)) {
    Sigma_e <- error_variance
    if (nrow(Sigma_e) != M || ncol(Sigma_e) != M ||
        !is.numeric(Sigma_e) || anyNA(Sigma_e) ||
        any(abs(Sigma_e - t(Sigma_e)) > 1e-8)) {
      stop(sprintf(paste0("When 'error_variance' is a covariance matrix it ",
                          "must be a symmetric numeric %d-by-%d matrix with ",
                          "no missing values."), M, M), call. = FALSE)
    }
    ev <- eigen(Sigma_e, symmetric = TRUE, only.values = TRUE)$values
    if (min(ev) < -1e-8) {
      stop("The 'error_variance' covariance matrix must be positive ",
           "semidefinite.", call. = FALSE)
    }
    err_var_marginal <- diag(Sigma_e)
  } else {
    if (!is.null(reliability)) {
      if (!is.numeric(reliability) || length(reliability) != 1L ||
          is.na(reliability) || reliability <= 0 || reliability >= 1) {
        stop("'reliability' must be a single number in (0, 1).", call. = FALSE)
      }
      if (all(true_var_by_occ <= 0)) {
        stop("'reliability' requires between-unit variation: give at ",
             "least one positive entry in 'random_variances'.", call. = FALSE)
      }
      # Homoscedastic error variance whose average per-occasion
      # reliability across the nominal times equals the target.
      solve_fun <- function(s2) {
        mean(true_var_by_occ / (true_var_by_occ + s2)) - reliability
      }
      upper <- max(true_var_by_occ) * (1 - reliability) / reliability * 10 + 1
      s2 <- stats::uniroot(solve_fun, c(1e-12, upper), tol = 1e-10)$root
      err_var_marginal <- rep(s2, M)
    } else {
      if (!is.numeric(error_variance) || anyNA(error_variance) ||
          any(error_variance < 0)) {
        stop("'error_variance' must be non-negative.", call. = FALSE)
      }
      if (length(error_variance) == 1L) {
        err_var_marginal <- rep(as.numeric(error_variance), M)
      } else if (length(error_variance) == M) {
        err_var_marginal <- as.numeric(error_variance)
      } else {
        stop(sprintf(paste0("'error_variance' must be a single number, a ",
                            "vector of length %d, or an %d-by-%d covariance ",
                            "matrix."), M, M, M), call. = FALSE)
      }
    }
    Rerr <- .error_cor_matrix(error_structure, error_correlation, M)
    sde <- sqrt(err_var_marginal)
    Sigma_e <- outer(sde, sde) * Rerr
  }

  reliability_by_occ <- ifelse(true_var_by_occ + err_var_marginal > 0,
                               true_var_by_occ /
                                 (true_var_by_occ + err_var_marginal),
                               NA_real_)
  }

  # ---------- Draw units ----------
  # Symmetric square roots rather than Cholesky factors: exact for the
  # positive semidefinite case (zero-variance parameters, a zero error
  # variance), where a ridge would leak noise into deterministic runs.
  psd_sqrt <- function(S) {
    e <- eigen(S, symmetric = TRUE)
    e$vectors %*% (sqrt(pmax(e$values, 0)) * t(e$vectors))
  }
  N <- sum(n_per_population)
  chol_T <- psd_sqrt(Tcov)
  chol_e <- if (individual_times) NULL else psd_sqrt(Sigma_e)

  out_list <- vector("list", G)
  subj <- 0L
  for (g in seq_len(G)) {
    ng <- n_per_population[g]
    theta_g <- matrix(rep(unlist(coef_list[[g]]), each = ng),
                      nrow = ng, dimnames = list(NULL, par_names))
    draws <- matrix(stats::rnorm(ng * K), ng, K) %*% chol_T
    theta <- theta_g + draws
    if (model == "richards" && any(theta[, "delta"] <= 0)) {
      stop(sum(theta[, "delta"] <= 0), " unit draw(s) produced a ",
           "Richards shape parameter delta <= 0, where the curve is ",
           "undefined. Shrink the delta entry of 'random_variances' (or ",
           "raise the fixed delta).", call. = FALSE)
    }

    rows <- vector("list", ng)
    for (i in seq_len(ng)) {
      subj <- subj + 1L
      if (individual_times) {
        m_i <- if (length(occasions) == 2L) {
          counts <- seq.int(occasions[1L], occasions[2L])
          counts[sample.int(length(counts), 1L)]
        } else {
          occasions
        }
        nominal <- sort(stats::runif(m_i, time_range[1L], time_range[2L]))
        actual <- nominal
        err <- stats::rnorm(m_i, 0, sqrt(err_scalar))
      } else {
        nominal <- target_times
        actual <- target_times + stats::rnorm(M, 0, timing_sd)
        err <- as.numeric(stats::rnorm(M) %*% chol_e)
      }
      p_i <- as.list(theta[i, ])
      true_score <- curve$mean(actual, p_i)
      rows[[i]] <- data.frame(
        id = sprintf("u%03d", subj),
        population = as.character(g),
        occasion = seq_along(actual),
        target_time = nominal,
        time = actual,
        true_score = true_score,
        y = true_score + err,
        stringsAsFactors = FALSE
      )
    }
    out_list[[g]] <- do.call(rbind, rows)
  }
  out <- do.call(rbind, out_list)
  out$id <- factor(out$id)
  out$population <- factor(out$population)
  rownames(out) <- NULL

  attr(out, "model") <- model
  attr(out, "fixed_parameters") <- coef_list
  attr(out, "random_covariance") <- Tcov
  attr(out, "schedule") <- if (individual_times) "unit_specific" else "shared"
  attr(out, "error_variance") <- if (individual_times ||
                                     (length(unique(err_var_marginal)) == 1L &&
                                      error_structure == "independent" &&
                                      !is.matrix(error_variance))) {
    err_var_marginal[1L]
  } else {
    err_var_marginal
  }
  attr(out, "error_covariance") <- if (individual_times) NA else Sigma_e
  attr(out, "reliability_by_occasion") <- reliability_by_occ
  out
}
