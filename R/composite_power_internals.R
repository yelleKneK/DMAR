# Internal helpers shared by the composite power functions.
#
#   .composite_power_shared_sigma()    composite power for single-degree-of-
#                                      freedom (t) tests that share one error
#                                      estimate; two-group ANCOVA
#   .composite_power_shared_sigma_F()  composite power for F tests of any
#                                      numerator degrees of freedom that share
#                                      one error estimate; factorial designs
#   .ancova_design()                   the two-group ANCOVA parameterization
#   .factorial_composite_design()      the factorial ANCOVA/ANOVA effect
#                                      parameterization
#   .composite_approximate_terms()     rename the rows the unequal residual
#                                      variance approximation touches
#   .composite_exact_terms()           the inverse renaming, for the readers
#   .composite_as_exact()              the inverse renaming applied to a table
#   .check_ancova_args()               shared validation (two-group)
#   .plot_ancova_effects()             the two-group population-effects figure
#   .plot_factorial_composite()        the factorial effects-and-power figure
#
# None are exported. Each figure is reached through the plot method on the
# object the corresponding planner returns.


# Composite power for tests that share an error estimate --------------------
#
# Every test in a linear model divides by the same estimated error standard
# deviation, so the test statistics are dependent even when the effects they
# test are orthogonal: an error estimate that lands low inflates all of them
# together. That dependence is what makes composite power larger than the
# product of the marginal powers, and it is why the composite cannot be
# obtained by multiplying.
#
# Conditional on the error estimate the statistics are independent, which gives
# a one-dimensional integral. Writing s for the estimated error standard
# deviation, u = s / sigma, and q = df * u^2 ~ chi square(df),
#
#   t_j = Z_j / u,   Z_j ~ N(delta_j, 1)  independent across j,
#
# so
#
#   composite = E_q[ prod_j P(reject test j | q) ]
#             = integral over q of prod_j P(|Z_j| > crit * sqrt(q / df))
#               times the chi square(df) density.
#
# Adaptive quadrature evaluates it, so nothing is simulated. Passing a single
# delta returns that test's ordinary noncentral t power, which is the identity
# the tests check.
#
# delta        vector of noncentrality parameters, one per test
# df           residual degrees of freedom shared by every test
# alpha_level  Type I error rate for each individual test
# directional  TRUE for one-sided tests, each in the direction of its own delta
.composite_power_shared_sigma <- function(delta, df, alpha_level,
                                          directional = FALSE) {
  if (df < 1) {
    return(list(marginal = rep(NA_real_, length(delta)), composite = NA_real_))
  }
  crit <- if (directional) {
    stats::qt(1 - alpha_level, df = df)
  } else {
    stats::qt(1 - alpha_level / 2, df = df)
  }

  # P(reject test j | u), for u = s / sigma. Z_j ~ N(delta_j, 1), t_j = Z_j / u.
  reject_given_u <- function(u, d) {
    if (directional) {
      # Reject in the direction of the effect: t > crit when d >= 0.
      if (d >= 0) stats::pnorm(d - crit * u) else stats::pnorm(-crit * u - d)
    } else {
      stats::pnorm(d - crit * u) + stats::pnorm(-crit * u - d)
    }
  }

  # Marginal power. Integrating a single test over the error distribution
  # reproduces pt(); the tests assert that identity.
  marginal <- vapply(delta, function(d) {
    if (directional) {
      if (d >= 0) {
        1 - stats::pt(crit, df = df, ncp = d)
      } else {
        stats::pt(-crit, df = df, ncp = d)
      }
    } else {
      (1 - stats::pt(crit, df = df, ncp = d)) +
        stats::pt(-crit, df = df, ncp = d)
    }
  }, numeric(1))

  if (length(delta) == 1L) {
    return(list(marginal = marginal, composite = unname(marginal)))
  }

  integrand <- function(q) {
    u <- sqrt(q / df)
    pr <- rep(1, length(q))
    for (d in delta) pr <- pr * reject_given_u(u, d)
    pr * stats::dchisq(q, df = df)
  }

  # Finite limits spanning essentially all of the chi square mass integrate more
  # reliably than (0, Inf) when df is large.
  lo <- stats::qchisq(1e-12, df = df)
  hi <- stats::qchisq(1 - 1e-12, df = df)
  composite <- stats::integrate(integrand, lower = lo, upper = hi,
                                rel.tol = .Machine$double.eps^0.5)$value

  # Quadrature can drift a hair outside [0, 1] at extreme deltas.
  composite <- min(max(composite, 0), 1)
  list(marginal = marginal, composite = composite)
}


# Two-group ANCOVA design ----------------------------------------------------
#
# The model is
#
#   Y = b_0 + b_group * G + b_cov * X + b_group_by_cov * G * X + e,
#
# with G coded -1/2 and +1/2 and X centered at its mean with standard deviation
# sd_cov. That coding is what makes the coefficients read directly: b_group is
# the mean difference between the groups at the covariate mean, b_cov is the
# average of the two within-group slopes, and b_group_by_cov is the difference
# between them. With balanced groups and a covariate whose distribution does not
# differ across them, the four columns are mutually orthogonal, so the three
# tests are orthogonal and the coefficient estimates are uncorrelated.
#
# sigma is the within-group standard deviation of the outcome, the same sigma as
# in a one-way ANOVA on the outcome, and it is the same in both groups. The user
# states the covariate-outcome correlation in each group. Within group g the
# slope is then rho_g * sigma / sd_cov and the residual variance is
# sigma^2 * (1 - rho_g^2). Equal correlations mean equal slopes and no
# interaction; different correlations mean the slopes differ and the interaction
# carries the difference. The same model covers both, so nothing has to be
# settled in advance and rho decides.
#
# The error variance the ANCOVA pools and estimates is the average of the two,
#
#   sigma_adj^2 = sigma^2 * (1 - mean(rho_1^2, rho_2^2)),
#
# which is sigma^2 * (1 - rho^2) whenever the two correlations are equal. When
# they are not, that average is the right mean and the wrong distribution; see
# the section on the unequal residual variance approximation below for what the
# planners do about it.
#
# sd_cov cancels out of every noncentrality below, as it must: a correlation is
# scale free, so the spread of the covariate cannot change the power. It is
# carried because the coefficients and the figure are expressed in the
# covariate's own units.
.ancova_design <- function(smd, rho, sigma, sd_cov, n, include_interaction) {
  rho_1 <- rho[1L]
  rho_2 <- rho[2L]

  # Pooled error SD: the average of the two within-group residual variances.
  rho_bar_sq <- mean(c(rho_1^2, rho_2^2))
  sigma_adj  <- sigma * sqrt(1 - rho_bar_sq)

  slope_1 <- rho_1 * sigma / sd_cov
  slope_2 <- rho_2 * sigma / sd_cov

  b_group        <- smd * sigma               # difference at the covariate mean
  b_cov          <- (slope_1 + slope_2) / 2   # average within-group slope
  b_group_by_cov <- slope_2 - slope_1         # difference between the slopes

  N  <- 2 * n
  p  <- if (include_interaction) 3L else 2L   # G, X, and G * X
  df <- N - p - 1L

  # Noncentralities, in the sqrt(N) * f form the rest of the ss_power_* family
  # uses (see ss_power_reg_coef). Under the coding above E[G^2] = 1/4,
  # E[X^2] = sd_cov^2, and E[(G * X)^2] = sd_cov^2 / 4.
  delta_group <- b_group * sqrt(N) / (2 * sigma_adj)
  delta_cov   <- b_cov * sd_cov * sqrt(N) / sigma_adj
  delta_int   <- b_group_by_cov * sd_cov * sqrt(N) / (2 * sigma_adj)

  delta <- c(group = delta_group, covariate = delta_cov)
  if (include_interaction) delta <- c(delta, group_by_covariate = delta_int)

  list(delta = delta, df = df, N = N, sigma_adj = sigma_adj,
       slope_1 = slope_1, slope_2 = slope_2,
       b_group = b_group, b_cov = b_cov, b_group_by_cov = b_group_by_cov)
}


# The unequal residual variance approximation --------------------------------
#
# Within group g the residual variance is tau_g^2 = sigma^2 (1 - rho_g^2), so
# correlations that differ in absolute value leave the groups with different
# residual variances. .ancova_design() above, and .factorial_het_effect_f()
# below, pool them by averaging the squared correlations. That gets the expected
# pooled error variance exactly right and its distribution only approximately,
# in two separate ways.
#
# The residual sum of squares is sum_g tau_g^2 Q_g with the Q_g independent chi
# square variables, a mixture of scaled chi squares rather than the single
# scaled chi square .composite_power_shared_sigma() integrates over. The mixture
# has the same mean and a larger variance: for two balanced groups with df_g
# residual degrees of freedom each it exceeds the single chi square's variance
# by df_g (tau_1^2 - tau_2^2)^2, which is zero exactly when the two residual
# variances agree.
#
# The numerators stop being uncorrelated as well. Under the -1/2, +1/2 coding
# the covariate coefficient is the average of the two within-group slopes and
# the interaction coefficient is their difference, so
#
#   Cov(b_cov_hat, b_group_by_cov_hat) = (tau_2^2 - tau_1^2) / (2 n sd_cov^2),
#
# again zero exactly when the residual variances agree. The conditional
# independence the integral rests on therefore fails too.
#
# Neither departure has a guaranteed sign, and neither is repaired here. What
# the planners do instead is say so in the table they return: in the unequal
# variance case every row carrying a power, and any sample size resolved against
# one, is renamed with an approximate_ prefix, so a reader who never opens the
# help page still sees which numbers are approximations.
# .composite_approximate_terms() applies that renaming and .composite_exact_terms()
# undoes it, which lets the plot methods and the tidy verbs go on looking up one
# fixed set of names.
.composite_approximate_terms <- function(term) {
  named <- c(composite_power       = "approximate_composite_power",
             necessary_n_per_group = "approximate_n_per_group",
             necessary_n_per_cell  = "approximate_n_per_cell",
             necessary_N           = "approximate_N")
  hit <- match(term, names(named))
  term[!is.na(hit)] <- unname(named[hit[!is.na(hit)]])
  marginal <- grepl("^power_", term)
  term[marginal] <- paste0("approximate_", term[marginal])
  term
}


.composite_exact_terms <- function(term) {
  named <- c(approximate_composite_power = "composite_power",
             approximate_n_per_group     = "necessary_n_per_group",
             approximate_n_per_cell      = "necessary_n_per_cell",
             approximate_N               = "necessary_N")
  hit <- match(term, names(named))
  term[!is.na(hit)] <- unname(named[hit[!is.na(hit)]])
  marginal <- grepl("^approximate_power_", term)
  term[marginal] <- sub("^approximate_", "", term[marginal])
  term
}


# Restore the exact-case row names on a returned table, so a reader written
# against one set of names keeps working whichever set the table carries. The
# plot methods use it; the tidy verbs instead recognize both sets of names
# through .SS_POWER_SIZE_TERMS and .SS_POWER_POWER_TERMS in R/dmar_tidiers.R,
# which is where the row names those verbs read are listed.
.composite_as_exact <- function(x) {
  x$term <- .composite_exact_terms(x$term)
  x
}


# Shared validation ----------------------------------------------------------
#
# Returns rho as a length-two vector; the rest is checked for the side effect of
# stopping.
.check_ancova_args <- function(smd, rho, sigma, sd_cov, alpha_level) {
  if (!is.numeric(smd) || length(smd) != 1L || !is.finite(smd)) {
    stop("'smd' must be a single finite numeric value.", call. = FALSE)
  }
  if (!is.numeric(rho) || !length(rho) %in% c(1L, 2L) || any(!is.finite(rho))) {
    stop("'rho' must be a numeric value of length 1 (the same correlation in ",
         "both groups) or length 2 (one correlation per group).", call. = FALSE)
  }
  if (any(rho <= -1) || any(rho >= 1)) {
    stop("'rho' must lie in (-1, 1).", call. = FALSE)
  }
  if (!is.numeric(sigma) || length(sigma) != 1L || sigma <= 0) {
    stop("'sigma' must be a single positive numeric value.", call. = FALSE)
  }
  if (!is.numeric(sd_cov) || length(sd_cov) != 1L || sd_cov <= 0) {
    stop("'sd_cov' must be a single positive numeric value.", call. = FALSE)
  }
  if (!is.numeric(alpha_level) || length(alpha_level) != 1L ||
      alpha_level <= 0 || alpha_level >= 1) {
    stop("'alpha_level' must be a single numeric value in (0, 1).", call. = FALSE)
  }
  if (length(rho) == 1L) rep(rho, 2L) else rho
}


# The population-effects figure ----------------------------------------------
#
# Drawn from the planning values alone: the two population regression lines the
# parameters describe. Nothing is simulated and no model is fitted; the figure
# is a picture of the effect sizes the plan was built on.
#
# Each effect is a feature of the picture. The vertical distance between the
# lines at the covariate mean is the group effect. The average tilt is the
# covariate effect. The difference between the tilts is the interaction.
#
# design       the .ancova_design() list
# power_info   NULL for the effects alone, else a list with the marginal powers,
#              the composite, the sample size, alpha_level and tails
.plot_ancova_effects <- function(design, smd, sigma, sd_cov, cov_range,
                                 palette, group_labels, power_info) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Plotting requires the ggplot2 package; install it to use ",
         "the plot method.", call. = FALSE)
  }

  # The group factor is coded -1/2 and +1/2, so the group effect splits
  # symmetrically about the grand mean and the vertical gap at the covariate
  # mean is exactly smd * sigma.
  b_group <- smd * sigma
  mean_1  <- -b_group / 2
  mean_2  <-  b_group / 2

  # Two points per line: the model is linear in the covariate.
  x_lo <- -cov_range * sd_cov
  x_hi <-  cov_range * sd_cov
  line_data <- data.frame(
    covariate = c(x_lo, x_hi, x_lo, x_hi),
    outcome = c(mean_1 + design$slope_1 * x_lo, mean_1 + design$slope_1 * x_hi,
                mean_2 + design$slope_2 * x_lo, mean_2 + design$slope_2 * x_hi),
    group = factor(rep(group_labels, each = 2L), levels = group_labels)
  )
  point_data <- data.frame(
    covariate = c(0, 0),
    outcome = c(mean_1, mean_2),
    group = factor(group_labels, levels = group_labels)
  )

  cols <- .dmar_palette(n = 2L, palette = palette)

  # Naming what each feature of the picture is, so the figure teaches the model
  # rather than needing the model to be read first.
  subtitle <- sprintf(
    "Gap at the mean = group effect (d = %.2f)   |   slope difference = interaction (%.3f)",
    smd, design$b_group_by_cov)

  p <- ggplot2::ggplot(line_data,
                       ggplot2::aes(x = .data$covariate, y = .data$outcome,
                                    color = .data$group)) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed",
                        color = "grey60", linewidth = 0.4) +
    ggplot2::geom_segment(x = 0, xend = 0, y = mean_1, yend = mean_2,
                          color = "grey35", linewidth = 0.6,
                          inherit.aes = FALSE) +
    ggplot2::geom_line(linewidth = 1.1) +
    ggplot2::geom_point(data = point_data, size = 3) +
    ggplot2::scale_color_manual(values = stats::setNames(cols, group_labels)) +
    ggplot2::labs(x = "Covariate (centered, in its own units)", y = "Outcome",
                  color = NULL,
                  title = "Population Effects the Design Is Planned On",
                  subtitle = subtitle) +
    ggplot2::theme(legend.position = "top")

  if (!is.null(power_info)) {
    power_bits <- sprintf("%s: %.3f", names(power_info$marginal),
                          power_info$marginal)
    caption <- paste0(
      sprintf("n = %.0f per group (N = %.0f), alpha = %.2f, %.0f-tailed",
              power_info$n, 2 * power_info$n, power_info$alpha_level,
              power_info$tails),
      "\nMarginal power   ", paste(power_bits, collapse = "   "),
      sprintf("\nComposite power (all %d jointly) = %.3f",
              length(power_info$marginal), power_info$composite))
    p <- p + ggplot2::labs(caption = caption) +
      ggplot2::theme(plot.caption = ggplot2::element_text(hjust = 0))
  }
  p
}


# Composite power for F tests sharing an error estimate --------------------
#
# The factorial analogue of .composite_power_shared_sigma(). Each effect in a
# factorial design is an F test with its own numerator degrees of freedom, and
# all of them divide by the same error mean square, so the tests are dependent
# even when the effects are orthogonal. Conditional on the error estimate the
# effect sums of squares are independent noncentral chi squares, which gives a
# one-dimensional integral over the chi square distribution of that estimate.
#
# Writing q for the scaled error sum of squares (q ~ chi square on df_error), the
# jth F test rejects when its noncentral chi square SS_j / sigma^2 exceeds
# df_num_j * crit_j * q / df_error, so
#
#   composite = integral over q of
#     prod_j P(chi^2'(df_num_j, lambda_j) > df_num_j * crit_j * q / df_error)
#     times the chi square(df_error) density.
#
# A single effect reproduces the ordinary noncentral F power exactly, which is
# the identity the tests assert. df_num, lambda, and (through recycling) the
# critical values are vectors, one entry per effect.
#
# df_num       numerator degrees of freedom, one per effect
# lambda       F noncentralities, one per effect
# df_error     residual degrees of freedom shared by every test
# alpha_level  Type I error rate for each individual test
.composite_power_shared_sigma_F <- function(df_num, lambda, df_error,
                                            alpha_level) {
  if (df_error < 1) {
    return(list(marginal = rep(NA_real_, length(lambda)), composite = NA_real_))
  }
  crit     <- stats::qf(1 - alpha_level, df_num, df_error)
  marginal <- stats::pf(crit, df_num, df_error, ncp = lambda, lower.tail = FALSE)

  if (length(lambda) == 1L) {
    return(list(marginal = marginal, composite = unname(marginal)))
  }

  integrand <- function(q) {
    pr <- rep(1, length(q))
    for (j in seq_along(lambda)) {
      x  <- df_num[j] * crit[j] * q / df_error
      pr <- pr * stats::pchisq(x, df_num[j], ncp = lambda[j], lower.tail = FALSE)
    }
    pr * stats::dchisq(q, df_error)
  }
  lo <- stats::qchisq(1e-12, df = df_error)
  hi <- stats::qchisq(1 - 1e-12, df = df_error)
  composite <- stats::integrate(integrand, lower = lo, upper = hi,
                                rel.tol = .Machine$double.eps^0.5)$value
  list(marginal = marginal, composite = min(max(composite, 0), 1))
}


# Parse and validate the 'effects' list of a factorial composite --------------
#
# Each element names the factors in the effect and its effect size, as either
# Cohen's f or partial eta squared. Returns a list of normalized effects, each
# with integer 'factors', a positive 'f', and a unique character 'label'.
.parse_factorial_effects <- function(effects, n_factors, require_f = TRUE) {
  if (!is.list(effects) || !length(effects)) {
    stop("'effects' must be a non-empty list; each element names the factors ",
         "in the effect (a vector of factor indices) and its effect size ('f' ",
         "or 'partial_eta_squared').", call. = FALSE)
  }
  parsed <- lapply(seq_along(effects), function(i) {
    e <- effects[[i]]
    if (!is.list(e)) {
      stop("Effect ", i, " must be a list with 'factors' and one of 'f' or ",
           "'partial_eta_squared'.", call. = FALSE)
    }
    # Exact element extraction throughout: e$f would partial-match 'factors',
    # so a forgotten effect size would be silently taken from the factor list.
    fac <- e[["factors"]]
    if (is.null(fac) || !is.numeric(fac) || !length(fac) ||
        any(is.na(fac)) || any(fac != round(fac)) || any(fac < 1) ||
        any(fac > n_factors) || anyDuplicated(fac)) {
      stop("Effect ", i, ": 'factors' must be unique integer indices into ",
           "'factor_levels' (1 to ", n_factors, ").", call. = FALSE)
    }
    fac <- sort(as.integer(fac))
    has_f   <- !is.null(e[["f"]])
    has_pes <- !is.null(e[["partial_eta_squared"]])
    if (!require_f) {
      # Means were supplied, so each effect's size is read from the means; an
      # effect size given here as well would be ambiguous.
      if (has_f || has_pes) {
        stop("Effect ", i, ": 'means' was supplied, so the effect sizes are ",
             "read from the means; do not also give 'f' or ",
             "'partial_eta_squared'.", call. = FALSE)
      }
      f <- NA_real_
    } else {
      if (has_f == has_pes) {
        stop("Effect ", i, ": supply exactly one of 'f' or ",
             "'partial_eta_squared'.", call. = FALSE)
      }
      if (has_pes) {
        pes <- e[["partial_eta_squared"]]
        if (!is.numeric(pes) || length(pes) != 1L || is.na(pes) || pes <= 0 ||
            pes >= 1) {
          stop("Effect ", i, ": 'partial_eta_squared' must be a single number ",
               "in (0, 1).", call. = FALSE)
        }
        f <- sqrt(pes / (1 - pes))
      } else {
        f <- e[["f"]]
        if (!is.numeric(f) || length(f) != 1L || is.na(f) || f <= 0) {
          stop("Effect ", i, ": 'f' must be a single positive number.",
               call. = FALSE)
        }
      }
    }
    label <- if (!is.null(e[["label"]])) as.character(e[["label"]])[1L] else
      paste(fac, collapse = "x")
    list(factors = fac, f = f, label = label)
  })
  # The same factor set named twice is the more informative complaint, so check
  # it before the label collision it would also produce under default labels.
  fac_keys <- vapply(parsed, function(e) paste(e$factors, collapse = ","),
                     character(1))
  if (anyDuplicated(fac_keys)) {
    stop("The same effect is named more than once in 'effects'.", call. = FALSE)
  }
  labels <- vapply(parsed, function(e) e$label, character(1))
  if (anyDuplicated(labels)) {
    stop("Two effects share the label '",
         labels[duplicated(labels)][1L],
         "'. Give the effects distinct factor sets or explicit 'label's.",
         call. = FALSE)
  }
  parsed
}


# The factorial ANCOVA / ANOVA effect parameterization ------------------------
#
# For a balanced factorial design with 'cells' = prod(factor_levels) cells and
# 'n_per_cell' per cell, each named effect has numerator degrees of freedom
# prod(levels - 1) over the factors it spans, and noncentrality N * f_adj^2 with
# f_adj = f / sqrt(1 - covariate_R2), the same covariate adjustment
# ss_power_factorial_ancova() uses: the covariate removes covariate_R2 of the
# error variance and spends n_covariates residual degrees of freedom. With no
# covariate (covariate_R2 = 0, n_covariates = 0) this is the factorial ANOVA.
.factorial_composite_design <- function(factor_levels, effects, covariate_R2,
                                        n_covariates, n_per_cell) {
  cells    <- prod(factor_levels)
  N        <- n_per_cell * cells
  df_error <- N - cells - n_covariates
  df_num   <- vapply(effects, function(e) prod(factor_levels[e$factors] - 1),
                     numeric(1))
  f        <- vapply(effects, function(e) e$f, numeric(1))
  f_adj    <- f / sqrt(1 - covariate_R2)
  list(label    = vapply(effects, function(e) e$label, character(1)),
       df_num   = df_num,
       f        = f,
       lambda   = N * f_adj^2,
       df_error = df_error,
       N        = N,
       cells    = cells)
}


# The factorial effects-and-power figure --------------------------------------
#
# A factorial design's supposed population values are its effect sizes; unlike a
# two-group ANCOVA they do not pin the cell means (many mean patterns share one
# Cohen's f), so the figure draws the effect sizes themselves. One lollipop per
# named effect gives its partial eta squared, colored and labeled by the
# marginal power that effect's size and the sample size deliver, with the
# composite power (all effects jointly) in the subtitle. Nothing is simulated.
.plot_factorial_composite <- function(effects_df, composite, n_per_cell, N,
                                      alpha_level, palette, title) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Plotting requires the ggplot2 package; install it to use the plot ",
         "method.", call. = FALSE)
  }
  d <- effects_df
  d$partial_eta_sq <- d$f^2 / (1 + d$f^2)
  d$label <- factor(d$label, levels = rev(d$label))
  cols <- .dmar_palette(n = 2L, palette = palette)

  ggplot2::ggplot(d, ggplot2::aes(x = .data$partial_eta_sq, y = .data$label)) +
    ggplot2::geom_segment(ggplot2::aes(x = 0, xend = .data$partial_eta_sq,
                                       yend = .data$label),
                          color = "grey70", linewidth = 0.8) +
    ggplot2::geom_point(ggplot2::aes(color = .data$marginal_power), size = 4) +
    ggplot2::geom_text(ggplot2::aes(
      label = sprintf("df %s, power %.2f", .data$df_num, .data$marginal_power)),
      hjust = -0.15, size = 3) +
    ggplot2::scale_color_gradient(low = cols[1L], high = cols[2L],
                                  limits = c(0, 1), name = "Marginal power") +
    ggplot2::expand_limits(x = max(d$partial_eta_sq) * 1.35) +
    ggplot2::labs(
      x = "Partial eta squared (purported population value)", y = NULL,
      title = title,
      subtitle = sprintf(
        "Composite power (all %d effects jointly) = %.3f    n = %d per cell (N = %d), alpha = %.2f",
        nrow(d), composite, n_per_cell, N, alpha_level)) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "top")
}


# Cohen's f for a factorial effect from the cell means -------------------------
#
# When the researcher supplies population cell means and a common within-cell
# standard deviation instead of effect sizes, each effect's f is read off the
# means by the orthogonal analysis of variance decomposition. The component of
# the means carried by the effect on factor set S is, by inclusion-exclusion,
#
#   comp_S(cell) = sum over T subset of S of (-1)^{|S| - |T|} * groupmean_T(cell),
#
# where groupmean_T averages the means within the cell's levels of the factors
# in T (the empty set gives the grand mean). The effect sum of squares is the
# sum of that component squared over the cells, and
#
#   f^2 = SS_S / (cells * sigma^2),
#
# so the noncentrality N * f^2 = n_per_cell * SS_S / sigma^2 matches the balanced
# fixed-effects convention the rest of the family uses. A two-group main effect
# recovers Cohen's f = |mean difference| / (2 sigma), and a one-way design
# recovers the standard f.
.power_set <- function(x) {
  out <- list(integer(0))
  for (e in x) out <- c(out, lapply(out, function(s) c(s, e)))
  out
}

.factorial_effect_f <- function(means, factors, sigma) {
  d      <- dim(means)
  n_cell <- prod(d)
  grid   <- expand.grid(lapply(d, seq_len))          # factor 1 varies fastest
  m_vec  <- as.vector(means)                          # column-major, matches grid
  comp   <- numeric(n_cell)
  for (tset in .power_set(factors)) {
    sign_t <- (-1)^(length(factors) - length(tset))
    gm <- if (!length(tset)) {
      rep(mean(m_vec), n_cell)
    } else {
      stats::ave(m_vec, interaction(grid[, tset, drop = FALSE], drop = TRUE))
    }
    comp <- comp + sign_t * gm
  }
  sqrt(sum(comp^2) / (n_cell * sigma^2))
}


# The factorial cell-means figure ---------------------------------------------
#
# When cell means are supplied, the figure shows them directly: a profile plot
# of the population means over the first factor, one line per level of the
# second, faceted by any further factors, with an error bar of plus or minus one
# within-cell standard deviation at each mean so the effect sizes read against
# the noise. The composite power and the sample size are in the subtitle.
.plot_factorial_means <- function(means, factor_levels, composite, n_per_cell,
                                  N, sigma, alpha_level, palette, title) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Plotting requires the ggplot2 package; install it to use the plot ",
         "method.", call. = FALSE)
  }
  nf   <- length(factor_levels)
  grid <- expand.grid(lapply(factor_levels, seq_len))
  names(grid) <- paste0("factor_", seq_len(nf))
  grid$mean <- as.vector(means)
  grid$ymin <- grid$mean - sigma
  grid$ymax <- grid$mean + sigma
  grid$factor_1 <- factor(grid$factor_1,
                          labels = paste0("1.", seq_len(factor_levels[1L])))
  grp <- if (nf >= 2L) {
    grid$factor_2 <- factor(grid$factor_2,
                            labels = paste0("2.", seq_len(factor_levels[2L])))
    "factor_2"
  } else NULL

  n_col <- if (is.null(grp)) 1L else nlevels(grid$factor_2)
  cols  <- .dmar_palette(n = max(n_col, 1L), palette = palette)

  p <- ggplot2::ggplot(grid, ggplot2::aes(x = .data$factor_1, y = .data$mean,
    group = if (is.null(grp)) 1 else .data$factor_2,
    color = if (is.null(grp)) NULL else .data$factor_2)) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = .data$ymin, ymax = .data$ymax),
                           width = 0.12, linewidth = 0.5,
                           position = ggplot2::position_dodge(width = 0.15)) +
    ggplot2::geom_line(linewidth = 1, position = ggplot2::position_dodge(width = 0.15)) +
    ggplot2::geom_point(size = 2.5, position = ggplot2::position_dodge(width = 0.15))

  if (!is.null(grp)) {
    p <- p + ggplot2::scale_color_manual(values = cols, name = "Factor 2")
  }
  if (nf >= 3L) {
    facet_by <- paste0("factor_", 3:nf)
    for (fb in facet_by) grid[[fb]] <- factor(grid[[fb]])
    p <- p + ggplot2::facet_wrap(stats::as.formula(
      paste("~", paste(facet_by, collapse = " + "))),
      labeller = ggplot2::label_both)
  }

  p +
    ggplot2::labs(
      x = "Factor 1 level", y = "Population cell mean",
      title = title,
      subtitle = sprintf(
        "Error bars are plus/minus one within-cell SD (%.3g). Composite power = %.3f, n = %d per cell (N = %d), alpha = %.2f",
        sigma, composite, n_per_cell, N, alpha_level)) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "top")
}


# Coerce a 'means' argument to an array whose dimensions are factor_levels.
# Accepts an array or matrix already so shaped, or a numeric vector of length
# prod(factor_levels) in array order (the first factor varying fastest).
.as_factorial_means <- function(means, factor_levels) {
  if (is.array(means)) {
    if (!identical(as.integer(dim(means)), as.integer(factor_levels))) {
      stop("'means' has dimensions ", paste(dim(means), collapse = " by "),
           " but 'factor_levels' is ", paste(factor_levels, collapse = " by "),
           "; they must match.", call. = FALSE)
    }
    storage.mode(means) <- "double"
    return(means)
  }
  if (is.numeric(means) && length(means) == prod(factor_levels) &&
      !anyNA(means)) {
    return(array(as.double(means), dim = factor_levels))
  }
  stop("'means' must be an array whose dimensions equal 'factor_levels', or a ",
       "numeric vector of length prod(factor_levels) in array order (the first ",
       "factor varying fastest).", call. = FALSE)
}


# Heterogeneous-slope factorial ANCOVA -----------------------------------------
#
# The factorial analogue of the two-group ss_power_composite_ancova_2group(): the
# covariate's slope may differ across the cells, so the covariate main effect
# (the average slope) and the factor-by-covariate slope heterogeneity are
# testable effects that can join the composite alongside the factorial mean
# effects. The full model fits the factorial means, the covariate, and every
# factor-by-covariate slope term, so it has 2 * cells parameters and
# df_error = N - 2 * cells. Under balance and a covariate with a common
# distribution across cells, all of those terms are orthogonal, so the tests
# share only the pooled residual and the shared-error F integral applies.
#
# Each effect is one of three kinds. A "mean" effect on a factor set S has
# numerator df prod(levels[S] - 1) and its size is read from the cell means,
# tested against the pooled adjusted error. A "covariate" effect is the average
# slope, numerator df 1. A "slope" effect on S is the factor-by-covariate slope
# heterogeneity, numerator df prod(levels[S] - 1), read from the cell slopes. In
# every case the noncentrality is N * f^2, so the design reuses the F integrator.
#
# .parse_factorial_het_effects() validates the effect list and the type taxonomy.
.parse_factorial_het_effects <- function(effects, n_factors, require_f = TRUE) {
  if (!is.list(effects) || !length(effects)) {
    stop("'effects' must be a non-empty list; each element names an effect's ",
         "'type' (\"mean\", \"covariate\", or \"slope\"), its 'factors' (for ",
         "\"mean\" and \"slope\"), and its size.", call. = FALSE)
  }
  parsed <- lapply(seq_along(effects), function(i) {
    e <- effects[[i]]
    if (!is.list(e)) stop("Effect ", i, " must be a list.", call. = FALSE)
    type <- if (is.null(e[["type"]])) "mean" else as.character(e[["type"]])[1L]
    if (!type %in% c("mean", "covariate", "slope")) {
      stop("Effect ", i, ": 'type' must be one of \"mean\", \"covariate\", ",
           "or \"slope\".", call. = FALSE)
    }
    fac <- e[["factors"]]
    if (type == "covariate") {
      if (!is.null(fac)) {
        stop("Effect ", i, ": a \"covariate\" effect is the average slope and ",
             "spans no factors; omit 'factors'.", call. = FALSE)
      }
      fac <- integer(0)
    } else {
      if (is.null(fac) || !is.numeric(fac) || !length(fac) || any(is.na(fac)) ||
          any(fac != round(fac)) || any(fac < 1) || any(fac > n_factors) ||
          anyDuplicated(fac)) {
        stop("Effect ", i, ": 'factors' must be unique integer indices into ",
             "'factor_levels' (1 to ", n_factors, ").", call. = FALSE)
      }
      fac <- sort(as.integer(fac))
    }
    has_f <- !is.null(e[["f"]]); has_pes <- !is.null(e[["partial_eta_squared"]])
    if (!require_f) {
      if (has_f || has_pes) {
        stop("Effect ", i, ": population values were supplied, so the sizes ",
             "are read from them; do not also give 'f' or ",
             "'partial_eta_squared'.", call. = FALSE)
      }
      f <- NA_real_
    } else if (has_f == has_pes) {
      stop("Effect ", i, ": supply exactly one of 'f' or 'partial_eta_squared'.",
           call. = FALSE)
    } else if (has_pes) {
      pes <- e[["partial_eta_squared"]]
      if (!is.numeric(pes) || length(pes) != 1L || is.na(pes) || pes <= 0 ||
          pes >= 1) {
        stop("Effect ", i, ": 'partial_eta_squared' must be a single number ",
             "in (0, 1).", call. = FALSE)
      }
      f <- sqrt(pes / (1 - pes))
    } else {
      f <- e[["f"]]
      if (!is.numeric(f) || length(f) != 1L || is.na(f) || f <= 0) {
        stop("Effect ", i, ": 'f' must be a single positive number.",
             call. = FALSE)
      }
    }
    label <- if (!is.null(e[["label"]])) as.character(e[["label"]])[1L] else
      switch(type,
             mean      = paste(fac, collapse = "x"),
             covariate = "covariate",
             slope     = paste0("cov_x_", paste(fac, collapse = "x")))
    list(type = type, factors = fac, f = f, label = label)
  })
  keys <- vapply(parsed, function(e)
    paste(e$type, paste(e$factors, collapse = ","), sep = "|"), character(1))
  if (anyDuplicated(keys)) {
    stop("The same effect is named more than once in 'effects'.", call. = FALSE)
  }
  labs <- vapply(parsed, function(e) e$label, character(1))
  if (anyDuplicated(labs)) {
    stop("Two effects share the label '", labs[duplicated(labs)][1L],
         "'. Give the effects distinct 'label's.", call. = FALSE)
  }
  parsed
}


# Cohen's f for a heterogeneous-slope effect from the population values. Mean
# effects are read from the cell means against the pooled adjusted error; the
# covariate main and slope effects are read from the cell slopes, which the
# supplied cell correlations imply (beta = rho * sigma / sd_cov). A "covariate"
# effect is the grand slope, the factor set being empty.
#
# The pooled adjusted error averages the squared cell correlations, which is the
# same pooling .ancova_design() does for two groups and carries the same
# approximation when the cell correlations differ in absolute value: the cells'
# residual variances are then sigma^2 (1 - rho_c^2) and differ, so the pooled
# error is a mixture of scaled chi squares. See the section on the unequal
# residual variance approximation above.
.factorial_het_effect_f <- function(effect, means, correlations, sigma, sd_cov) {
  sig_adj <- sigma * sqrt(1 - mean(correlations^2))
  if (effect$type == "mean") {
    .factorial_effect_f(means, effect$factors, sig_adj)
  } else {
    beta <- correlations * sigma / sd_cov
    .factorial_effect_f(beta, effect$factors, sig_adj / sd_cov)
  }
}


# The per-effect degrees of freedom and noncentralities of the design.
.factorial_het_design <- function(factor_levels, effects, n_per_cell) {
  cells    <- prod(factor_levels)
  N        <- n_per_cell * cells
  df_error <- N - 2L * cells
  df_num   <- vapply(effects, function(e)
    if (e$type == "covariate") 1 else prod(factor_levels[e$factors] - 1),
    numeric(1))
  f <- vapply(effects, function(e) e$f, numeric(1))
  list(label    = vapply(effects, function(e) e$label, character(1)),
       type     = vapply(effects, function(e) e$type, character(1)),
       df_num   = df_num, f = f, lambda = N * f^2,
       df_error = df_error, N = N, cells = cells)
}


# The heterogeneous-slope figure: one population regression line per cell -------
#
# Drawn from the planning values alone. In cell c the line is
# mu_c + beta_c * (X - Xbar) over the covariate range, with beta_c the slope the
# cell's correlation implies. Heterogeneous slopes show as lines of different
# angle, mean effects as vertical separation. Colored by the first factor and
# faceted by any others. Nothing is simulated.
.plot_factorial_het <- function(means, correlations, sigma, sd_cov,
                                factor_levels, composite, n_per_cell, N,
                                alpha_level, palette, cov_range, title) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Plotting requires the ggplot2 package; install it to use the plot ",
         "method.", call. = FALSE)
  }
  nf   <- length(factor_levels)
  grid <- expand.grid(lapply(factor_levels, seq_len))
  names(grid) <- paste0("factor_", seq_len(nf))
  mu   <- as.vector(means)
  beta <- as.vector(correlations) * sigma / sd_cov
  xr   <- c(-cov_range * sd_cov, cov_range * sd_cov)

  d <- do.call(rbind, lapply(seq_len(nrow(grid)), function(k) {
    data.frame(grid[c(k, k), , drop = FALSE],
               covariate = xr, outcome = mu[k] + beta[k] * xr,
               cell = k, row.names = NULL)
  }))
  d$cell     <- factor(d$cell)
  d$factor_1 <- factor(d$factor_1,
                       labels = paste0("1.", seq_len(factor_levels[1L])))
  cols <- .dmar_palette(n = factor_levels[1L], palette = palette)

  p <- ggplot2::ggplot(d, ggplot2::aes(x = .data$covariate, y = .data$outcome,
                                       group = .data$cell,
                                       color = .data$factor_1)) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed",
                        color = "grey70", linewidth = 0.4) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::scale_color_manual(values = cols, name = "Factor 1")

  if (nf >= 2L) {
    for (j in 2:nf) grid[[paste0("factor_", j)]] <- factor(grid[[paste0("factor_", j)]])
    facet_by <- paste0("factor_", 2:nf)
    for (fb in facet_by) d[[fb]] <- factor(d[[fb]])
    p <- p + ggplot2::facet_wrap(stats::as.formula(
      paste("~", paste(facet_by, collapse = " + "))),
      labeller = ggplot2::label_both)
  }

  p +
    ggplot2::labs(
      x = "Covariate (centered, in its own units)", y = "Outcome",
      title = title,
      subtitle = sprintf(
        "Line angle is the covariate slope (heterogeneous when they differ); vertical gap is the mean effect. Composite power = %.3f, n = %d per cell (N = %d), alpha = %.2f",
        composite, n_per_cell, N, alpha_level)) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "top")
}
