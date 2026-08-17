# Internal helpers shared by the SEM composite planners.
#
#   .composite_sem_setup()        resolve the population covariance matrix,
#                                 the parameter labels, the population value
#                                 each labeled parameter takes under the
#                                 analysis model, and its asymptotic variance
#   .composite_sem_mc()           the a priori Monte Carlo engine: G converged
#                                 fits of the analysis model to data drawn
#                                 from Sigma at sample size N
#   .composite_sem_search()       the smallest integer N whose Monte Carlo
#                                 criterion holds, by geometric bracketing
#                                 and integer bisection
#   .composite_sem_start_power()  analytic Wald start value for the power
#                                 search
#   .composite_sem_start_width()  analytic start value for the width search
#
# None are exported. ss_power_composite_sem() and ss_aipe_composite_sem() are
# the public faces; both require lavaan (Suggests) and gate it before calling
# in here.


# Resolve the planning inputs shared by both planners ------------------------
#
# The population is stated either as a covariance matrix Sigma (typically from
# cov_sem()), optionally with a mean vector mu, or as a fully fixed lavaan
# model pop_model from which cov_sem() derives both. The analysis model, fit
# to those population moments at a very large sample size, gives the
# population value each labeled parameter takes under the analysis model
# (equal to the generating values when the model is correctly specified, and
# the probability limits of the estimates when it is not) and, through N
# times the sampling variance at that N, the asymptotic variance of each
# estimate. The asymptotic variances seed the sample size searches; the
# Monte Carlo evaluations own the answer.
#
# The population means matter only when the analysis model has a mean
# structure (any ~1 term, as in a latent growth curve model); a covariance
# structure analysis is unmoved by them either way.
#
# parameters = NULL resolves to every user-labeled parameter in the analysis
# model, in order of appearance, including quantities defined with := (an
# indirect effect, a difference between two paths), whose standard errors
# lavaan obtains by the delta method.
.composite_sem_setup <- function(model, Sigma, pop_model, parameters,
                                 mu = NULL, ...) {
  if (!is.character(model) || length(model) != 1L || is.na(model)) {
    stop("'model' must be a single character string of lavaan model syntax.",
         call. = FALSE)
  }
  if (is.null(Sigma) == is.null(pop_model)) {
    stop("Supply exactly one of 'Sigma' (a population covariance matrix of ",
         "the observed variables, typically from cov_sem()) or 'pop_model' ",
         "(a fully fixed lavaan population model from which cov_sem() ",
         "derives it).", call. = FALSE)
  }
  if (!is.null(pop_model)) {
    if (!is.null(mu)) {
      stop("'pop_model' states the population means through its mean ",
           "structure (or their absence); do not also give 'mu'.",
           call. = FALSE)
    }
    population <- cov_sem(pop_model)
    Sigma <- population$sigma_theta
    mu <- population$mu_theta
  }
  Sigma <- as.matrix(Sigma)
  if (is.null(rownames(Sigma)) || is.null(colnames(Sigma))) {
    stop("'Sigma' must have row and column names matching the observed ",
         "variables in 'model'.", call. = FALSE)
  }
  if (is.null(mu)) {
    mu <- stats::setNames(rep(0, nrow(Sigma)), rownames(Sigma))
  } else {
    if (!is.numeric(mu) || anyNA(mu) || any(!is.finite(mu))) {
      stop("'mu' must be a finite numeric vector of population means.",
           call. = FALSE)
    }
    if (is.null(names(mu))) {
      if (length(mu) != nrow(Sigma)) {
        stop("An unnamed 'mu' must have one entry per row of 'Sigma', in ",
             "the same order; name the entries to match by variable ",
             "instead.", call. = FALSE)
      }
      names(mu) <- rownames(Sigma)
    } else {
      if (!setequal(names(mu), rownames(Sigma)) ||
          anyDuplicated(names(mu))) {
        stop("A named 'mu' must have exactly one entry per observed ",
             "variable in 'Sigma' (", paste(rownames(Sigma), collapse = ", "),
             ").", call. = FALSE)
      }
      mu <- mu[rownames(Sigma)]
    }
  }

  # A mean structure in the analysis model (any ~1 term) makes the population
  # means part of the fit; without one lavaan works from the covariance
  # structure alone and no sample.mean is passed.
  has_means <- any(lavaan::lavaanify(model)$op == "~1")
  N_big <- 1e6
  fit_args <- list(model, sample.cov = Sigma, sample.nobs = N_big, ...)
  if (has_means) fit_args$sample.mean <- mu
  fit0 <- do.call(lavaan::sem, fit_args)
  pe0 <- lavaan::parameterEstimates(fit0, ci = FALSE)

  available <- unique(pe0$label[pe0$label != ""])
  if (is.null(parameters)) {
    if (!length(available)) {
      stop("No labeled parameters were found in 'model'. Label each ",
           "parameter of interest in the lavaan syntax (for example ",
           "'f2 ~ b*f1', or 'ab := a*b' for a defined quantity) and, if you ",
           "want a subset of the labels, name it in 'parameters'.",
           call. = FALSE)
    }
    parameters <- available
  } else {
    if (!is.character(parameters) || !length(parameters) ||
        anyNA(parameters) || anyDuplicated(parameters)) {
      stop("'parameters' must be a character vector of distinct parameter ",
           "labels.", call. = FALSE)
    }
    missing_labels <- setdiff(parameters, available)
    if (length(missing_labels)) {
      stop("These entries of 'parameters' are not labeled parameters of ",
           "'model': ", paste(missing_labels, collapse = ", "),
           ". Label each parameter of interest in the lavaan syntax, for ",
           "example 'f2 ~ b*f1'.", call. = FALSE)
    }
  }

  # First row per label: an equality-constrained label appears once per
  # constrained parameter with identical estimate and standard error.
  idx <- match(parameters, pe0$label)
  theta <- pe0$est[idx]
  se_big <- pe0$se[idx]
  if (anyNA(se_big) || any(!is.finite(se_big)) || any(se_big <= 0)) {
    bad <- parameters[is.na(se_big) | !is.finite(se_big) | se_big <= 0]
    stop("The analysis model gives no sampling variability for ",
         paste(bad, collapse = ", "),
         "; a parameter of interest must be free (or defined from free ",
         "parameters), not fixed.", call. = FALSE)
  }

  obs_vars <- lavaan::lavNames(fit0, "ov")
  if (!all(obs_vars %in% rownames(Sigma))) {
    stop("'Sigma' must cover every observed variable in 'model'; missing: ",
         paste(setdiff(obs_vars, rownames(Sigma)), collapse = ", "), ".",
         call. = FALSE)
  }

  list(Sigma = Sigma, mu = mu, labels = parameters, theta = theta,
       h = se_big^2 * N_big, n_obs_vars = length(obs_vars))
}


# The a priori Monte Carlo engine ---------------------------------------------
#
# Draw a sample of size N from the multivariate normal population with mean
# vector mu and covariance matrix Sigma, fit the analysis model, and record
# each labeled parameter's estimate and standard error; repeat until G
# replications converge. A replication that fails to converge, or converges
# without a usable standard error for some target (missing, infinite, or
# zero), is redrawn, up to 20 * G attempts in total (as in
# ss_aipe_sem_path_sensitivity()); the caller decides what to do when fewer
# than G converge.
.composite_sem_mc <- function(model, Sigma, mu, labels, N, G, ...) {
  # Borderline samples make lavaan warn (nonconvergence, negative variances);
  # muffle for the duration of the replications and restore the user's
  # setting on exit, as in ss_aipe_sem_path_sensitivity().
  prev_warn <- getOption("warn")
  on.exit(options(warn = prev_warn), add = TRUE)
  options(warn = -1)

  vars <- rownames(Sigma)
  k <- length(labels)
  est <- matrix(NA_real_, nrow = G, ncol = k)
  se <- matrix(NA_real_, nrow = G, ncol = k)

  g <- 0L
  attempts <- 0L
  max_attempts <- 20L * G
  while (g < G && attempts < max_attempts) {
    attempts <- attempts + 1L
    Data <- MASS::mvrnorm(n = N, mu = mu, Sigma = Sigma)
    colnames(Data) <- vars
    # suppressMessages: lavaan narrates marker-item fallbacks through
    # message(), which try(silent = TRUE) does not muffle, so a hard
    # replication would otherwise print a note per failed fit inside the
    # Monte Carlo loop. Errors and nonconvergence are still counted below.
    fit <- try(suppressMessages(
      lavaan::sem(model, data = as.data.frame(Data), ...)),
      silent = TRUE)
    if (inherits(fit, "try-error") ||
        !isTRUE(lavaan::lavInspect(fit, "converged"))) {
      next
    }
    pe <- try(lavaan::parameterEstimates(fit, ci = FALSE), silent = TRUE)
    if (inherits(pe, "try-error")) next
    row_idx <- match(labels, pe$label)
    est_g <- pe$est[row_idx]
    se_g <- pe$se[row_idx]
    if (anyNA(est_g) || anyNA(se_g) || any(!is.finite(se_g)) ||
        any(se_g <= 0)) {
      next
    }
    g <- g + 1L
    est[g, ] <- est_g
    se[g, ] <- se_g
  }
  if (g < G) {
    est <- est[seq_len(g), , drop = FALSE]
    se <- se[seq_len(g), , drop = FALSE]
  }
  list(est = est, se = se, converged = g, requested = G, attempts = attempts)
}


# The smallest integer N whose Monte Carlo criterion holds --------------------
#
# evaluate(N) returns an evaluation summary and met(summary) says whether the
# criterion held there. The population criterion is monotone in N (power
# rises, widths shrink), so the search geometrically brackets a crossing from
# the start value and bisects the bracket to adjacent integers. The estimates
# are Monte Carlo, so the bracketing and the bisection act on estimated
# criteria: the returned N carries the simulation error of its G-replication
# evaluations, and a larger G resolves it more sharply.
.composite_sem_search <- function(evaluate, met, N_start, N_min,
                                  N_max = 1e6) {
  N <- max(N_min, N_start)
  ev <- evaluate(N)
  if (met(ev)) {
    hi <- N
    hi_ev <- ev
    repeat {
      if (hi <= N_min) {
        lo <- hi - 1L
        break
      }
      cand <- max(N_min, min(floor(hi / 1.25), hi - 1L))
      ev2 <- evaluate(cand)
      if (met(ev2)) {
        hi <- cand
        hi_ev <- ev2
      } else {
        lo <- cand
        break
      }
    }
  } else {
    lo <- N
    repeat {
      cand <- min(N_max, ceiling(lo * 1.25) + 1L)
      ev2 <- evaluate(cand)
      if (met(ev2)) {
        hi <- cand
        hi_ev <- ev2
        break
      }
      if (cand >= N_max) {
        stop("The goals were not met by N = ", N_max, ". Check that every ",
             "parameter of interest is nonzero in the population (for ",
             "power) and that the desired widths are attainable, or raise ",
             "'G' if the Monte Carlo estimates are too noisy near the ",
             "target.", call. = FALSE)
      }
      lo <- cand
    }
  }
  while (hi - lo > 1L) {
    mid <- floor((lo + hi) / 2)
    ev2 <- evaluate(mid)
    if (met(ev2)) {
      hi <- mid
      hi_ev <- ev2
    } else {
      lo <- mid
    }
  }
  list(N = as.integer(hi), eval = hi_ev)
}


# Analytic start value for the power search -----------------------------------
#
# The smallest N at which the product of the marginal Wald powers reaches the
# target. The product treats the tests as independent, which they are not;
# the Monte Carlo search owns the final answer and this only chooses where it
# begins. Returns N_max when even that sample size falls short, which the
# caller turns into an informative error before any simulation is spent.
.composite_sem_start_power <- function(theta, h, alpha_level, desired_power,
                                       N_min, N_max = 1e6) {
  z_crit <- stats::qnorm(1 - alpha_level / 2)
  product_power_at <- function(N) {
    delta <- abs(theta) / sqrt(h / N)
    prod(stats::pnorm(delta - z_crit) + stats::pnorm(-delta - z_crit))
  }
  N <- N_min
  while (product_power_at(N) < desired_power && N < N_max) {
    N <- min(N_max, ceiling(N * 1.25))
  }
  N
}


# Analytic start value for the width search -----------------------------------
#
# The Wald interval for parameter j has expected width about
# 2 z sqrt(h_j / N), so N_j = 4 z^2 h_j / omega_j^2 meets its desired width
# omega_j and the largest N_j is where the joint search begins (the same
# no-assurance approximation ss_aipe_sem_path() uses).
.composite_sem_start_width <- function(h, omega, conf_level, N_min) {
  z <- stats::qnorm(1 - (1 - conf_level) / 2)
  max(N_min, ceiling(max(4 * z^2 * h / omega^2)))
}


# Monte Carlo resolution guard ------------------------------------------------
#
# Both planners stop at the smallest N whose Monte Carlo estimate crosses the
# requested target. That estimate is a proportion of G replications, so it
# lives on the grid 0, 1/G, ..., 1 and cannot distinguish a target within 1/G
# of 1 from certainty: the search stops as soon as every replication happens
# to be significant (or narrow enough), which occurs at a sample size below
# the one the target actually requires, and the returned table then reports a
# realized value of 1 with a simulation standard error of 0. Rather than
# return a size the simulation cannot justify, refuse the combination and say
# how large G would have to be. The bound is deliberately stated in terms of
# the user's own arguments so the remedy is arithmetic, not guesswork.
.composite_sem_check_resolution <- function(target, G, arg_name) {
  if (target > 1 - 1 / G) {
    stop("'", arg_name, "' = ", format(target), " is finer than ", G,
         " Monte Carlo replications can resolve: a proportion of ", G,
         " replications takes only the values 0, 1/", G, ", ..., 1, so any ",
         "target above ", format(1 - 1 / G), " is met the moment every ",
         "replication succeeds, at a sample size below the one the target ",
         "requires. Raise 'G' to at least ", ceiling(1 / (1 - target)),
         " for this target, or lower '", arg_name, "'.", call. = FALSE)
  }
  invisible(TRUE)
}
