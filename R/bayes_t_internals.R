# =============================================================================
#  Internal numerics for the Bayesian t-test family
# =============================================================================
#
# Shared by bayes_one_sample_t(), bayes_paired_t(), and
# bayes_independent_t(). Nothing here is exported.
#
# The model is the JZS (Jeffreys-Zellner-Siow) setup of Rouder, Speckman,
# Sun, Morey, and Iverson (2009): a Jeffreys prior on the nuisance
# (mean-zero location and variance) and, by default, a Cauchy(0, r) prior
# on the standardized effect delta. Two informed variants are available:
# the Cauchy may be given a nonzero location (Gronau, Ly, and Wagenmakers,
# 2020), or the prior may instead be a normal with a stated mean and
# standard deviation, for the researcher who thinks in prior moments
# (a Cauchy has neither a mean nor a variance, so those beliefs cannot be
# expressed through the Cauchy at all). The two families are linked: a
# Cauchy(mu, r) is exactly a normal prior N(mu, g) whose variance g is
# itself uncertain, g ~ inverse-gamma(1/2, r^2/2), a fact the test suite
# exploits to cross-validate the two code paths against each other. Everything reduces to one-dimensional
# integrals against the noncentral t likelihood of the observed t
# statistic:
#
#   p(t | delta) = dt(t; df, ncp = delta * sqrt(n_eff))
#
# with n_eff the effective sample size (n for one sample / paired, the
# harmonic-mean-type n1 n2 / (n1 + n2) for two independent samples). The
# Bayes factor BF10 is the ratio of the Cauchy-marginal likelihood to the
# point-null likelihood, and the posterior of delta is proportional to
# likelihood times prior, normalized by quadrature.

# Build and validate the prior on delta from the public arguments.
# Exactly one family: the Cauchy (prior_location, prior_scale), which is
# the JZS default, or the normal (prior_mean, prior_sd). Mixing arguments
# from the two families is an error rather than a guess.
.bayes_t_prior <- function(prior_location = 0, prior_scale = sqrt(2) / 2,
                           prior_mean = NULL, prior_sd = NULL,
                           cauchy_args_supplied = FALSE) {
  normal <- !is.null(prior_mean) || !is.null(prior_sd)
  if (normal) {
    if (is.null(prior_mean) || is.null(prior_sd)) {
      stop("The normal prior needs both 'prior_mean' and 'prior_sd'.",
           call. = FALSE)
    }
    if (cauchy_args_supplied || prior_location != 0) {
      stop("Specify one prior family: either the Cauchy through ",
           "'prior_location' and 'prior_scale', or the normal through ",
           "'prior_mean' and 'prior_sd', not both.", call. = FALSE)
    }
    for (nm in c("prior_mean", "prior_sd")) {
      v <- get(nm)
      if (!is.numeric(v) || length(v) != 1L || is.na(v)) {
        stop("'", nm, "' must be a single non-missing number.",
             call. = FALSE)
      }
    }
    if (prior_sd <= 0) {
      stop("'prior_sd' must be positive.", call. = FALSE)
    }
    return(list(family = "normal", location = prior_mean, scale = prior_sd,
                dens = function(delta) {
                  stats::dnorm(delta, mean = prior_mean, sd = prior_sd)
                }))
  }
  for (nm in c("prior_location", "prior_scale")) {
    v <- get(nm)
    if (!is.numeric(v) || length(v) != 1L || is.na(v)) {
      stop("'", nm, "' must be a single non-missing number.", call. = FALSE)
    }
  }
  if (prior_scale <= 0) {
    stop("'prior_scale' must be a single positive number.", call. = FALSE)
  }
  list(family = "cauchy", location = prior_location, scale = prior_scale,
       dens = function(delta) {
         stats::dcauchy(delta, location = prior_location,
                        scale = prior_scale)
       })
}

# Marginal likelihood of t under the prior on delta.
.bayes_t_marginal <- function(t_obs, df, n_eff, prior) {
  f <- function(delta) {
    # suppressWarnings: R's noncentral t density (pnt) emits a benign
    # "full precision may not have been achieved" note across the grid; the
    # accuracy is far beyond what these posterior summaries require.
    suppressWarnings(stats::dt(t_obs, df = df, ncp = delta * sqrt(n_eff))) *
      prior$dens(delta)
  }
  stats::integrate(f, lower = -Inf, upper = Inf,
                   rel.tol = 1e-10)$value
}

# Posterior summaries of delta on an adaptive grid: normalizing constant,
# mean, median, credible interval, and P(delta > 0).
.bayes_t_posterior <- function(t_obs, df, n_eff, prior,
                               conf_level = 0.95) {
  unnorm <- function(delta) {
    # See .bayes_t_marginal: muffle the benign noncentral t precision note.
    suppressWarnings(stats::dt(t_obs, df = df, ncp = delta * sqrt(n_eff))) *
      prior$dens(delta)
  }
  # Span the grid over the sample standardized effect and the prior
  # location together, with generous width, so an informed prior far from
  # the data cannot fall off the grid.
  d_hat  <- t_obs / sqrt(n_eff)
  center <- c(d_hat, prior$location)
  # The prior's contribution to the span is capped: with a very diffuse
  # prior the posterior is confined by the likelihood, and letting the
  # grid stretch to the prior's tails would coarsen the quantile step for
  # no gain.
  half   <- max(6 / sqrt(n_eff), 4 * min(prior$scale, 2),
                2 * abs(d_hat), 1)
  lo <- min(center) - half
  hi <- max(center) + half
  grid <- seq(lo, hi, length.out = 8001L)
  dens <- unnorm(grid)
  step <- grid[2L] - grid[1L]
  Z    <- sum(dens) * step
  dens <- dens / Z
  cdf  <- cumsum(dens) * step

  q_at <- function(p) grid[which.min(abs(cdf - p))]
  alpha <- 1 - conf_level
  list(
    mean    = sum(grid * dens) * step,
    median  = q_at(0.5),
    lower   = q_at(alpha / 2),
    upper   = q_at(1 - alpha / 2),
    p_pos   = 1 - stats::approx(grid, cdf, xout = 0, rule = 2)$y,
    grid    = grid,
    density = dens
  )
}

# Assemble the full table common to the three public functions.
# Quantities reported on the standardized (delta) scale plus the raw mean
# (difference) scale via the plug-in s.
.bayes_t_table <- function(t_obs, df, n_eff, s_raw, prior,
                           conf_level, extra_terms, extra_values) {
  post <- .bayes_t_posterior(t_obs, df, n_eff, prior, conf_level)
  m1   <- .bayes_t_marginal(t_obs, df, n_eff, prior)
  m0   <- stats::dt(t_obs, df = df, ncp = 0)
  bf10 <- m1 / m0

  # The prior rows keep one shape across both families: location and
  # scale mean the Cauchy's parameters for the default family and the
  # normal's mean and standard deviation otherwise; the family itself is
  # recorded as an attribute, following the package's rule that
  # non-numeric metadata rides on attributes.
  out <- data.frame(
    term  = c("delta_posterior_median", "delta_posterior_mean",
              "delta_lower", "delta_upper", "p_delta_positive",
              "raw_posterior_median", "raw_lower", "raw_upper",
              "bf_10", "bf_01", "t", "df",
              "prior_location", "prior_scale",
              extra_terms),
    value = c(post$median, post$mean, post$lower, post$upper, post$p_pos,
              post$median * s_raw, post$lower * s_raw, post$upper * s_raw,
              bf10, 1 / bf10, t_obs, df,
              prior$location, prior$scale,
              extra_values),
    stringsAsFactors = FALSE
  )
  out <- .as_dmar_tbl(out, conf_level = conf_level)
  attr(out, "prior_family") <- prior$family
  # The full posterior of delta, for plotting or for probabilities the
  # summary rows do not report; density integrates to one over the grid.
  attr(out, "posterior") <- data.frame(delta = post$grid,
                                       density = post$density)
  out
}

# Validate a (mean, sd, n) summary triple for the summary-statistics
# entry path shared by the three bayes_*_t() functions. The names of the
# reported arguments follow the caller's vocabulary (mean_1/sd_1/n_1 for
# the independent groups form, mean_diff/sd_diff/n for the paired form),
# so the labels are passed in.
.bayes_t_check_summaries <- function(mean, sd, n,
                                     labels = c("mean", "sd", "n"),
                                     min_n = 2L) {
  vals <- list(mean, sd, n)
  for (i in seq_along(vals)) {
    v <- vals[[i]]
    if (is.null(v) || !is.numeric(v) || length(v) != 1L || is.na(v)) {
      stop("'", labels[i], "' must be a single non-missing number ",
           "(all of '", paste(labels, collapse = "', '"),
           "' are required for the summary-statistics form).",
           call. = FALSE)
    }
  }
  if (sd <= 0) {
    stop("'", labels[2], "' must be positive.", call. = FALSE)
  }
  if (n < min_n || n != round(n)) {
    stop("'", labels[3], "' must be a whole number of at least ",
         min_n, ".", call. = FALSE)
  }
  invisible(TRUE)
}

.bayes_t_check_x <- function(x, nm, min_n = 2L) {
  if (!is.numeric(x) || anyNA(x) || length(x) < min_n) {
    stop(sprintf("'%s' must be a numeric vector with at least %d ",
                 nm, min_n), "non-missing values.", call. = FALSE)
  }
  invisible(TRUE)
}
