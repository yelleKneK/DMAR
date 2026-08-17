## from tests/testthat/test-irt_grm.R
local({
  # irt_grm against mirt's marginal maximum likelihood fit of the same
  # graded response data. Exact agreement is not expected: lavaan fits by
  # limited information weighted least squares on the polychoric
  # correlations, mirt by marginal maximum likelihood on the full response
  # patterns, so the two are consistent but differ on a finite sample. mirt
  # reports slopes in the logistic metric, which is the normal ogive slope
  # multiplied by 1.702.
  grm_simulate <- function(n, a, b, seed) {
    set.seed(seed)
    theta <- stats::rnorm(n)
    responses <- vapply(seq_along(a), function(i) {
      p_star <- outer(theta, b[i, ],
                      function(z, bk) stats::pnorm(a[i] * (z - bk)))
      as.integer(1 + rowSums(stats::runif(n) < p_star))
    }, integer(n))
    colnames(responses) <- paste0("item", seq_along(a))
    as.data.frame(responses)
  }
  a_pop <- c(1.2, 0.8, 1.6, 1.0, 1.4, 1.1)
  b_pop <- rbind(c(-1.6, -0.6, 0.3, 1.2),
                 c(-1.4, -0.4, 0.5, 1.5),
                 c(-1.8, -0.7, 0.2, 1.1),
                 c(-1.2, -0.2, 0.7, 1.6),
                 c(-1.5, -0.5, 0.4, 1.3),
                 c(-1.3, -0.3, 0.6, 1.4))
  d <- grm_simulate(5000, a_pop, b_pop, seed = 113)
  res <- DMAR::irt_grm(d)
  a_lavaan <- res$a[!duplicated(res$item)]
  b_lavaan <- matrix(res$b, nrow = nrow(b_pop), byrow = TRUE)

  mirt_fit <- mirt::mirt(d, 1, itemtype = "graded", verbose = FALSE)
  mirt_par <- mirt::coef(mirt_fit, IRTpars = TRUE, simplify = TRUE)$items
  a_mirt <- mirt_par[, "a"] / 1.702                # logistic -> normal ogive
  b_mirt <- mirt_par[, paste0("b", seq_len(ncol(b_pop)))]

  stopifnot(
    stats::cor(a_lavaan, a_mirt) > 0.99,
    stats::cor(as.vector(b_lavaan), as.vector(as.matrix(b_mirt))) > 0.99,
    mean(abs(a_lavaan - a_mirt)) < 0.10,
    mean(abs(b_lavaan - as.matrix(b_mirt))) < 0.10
  )
})

## from tests/testthat/test-irt_information.R
local({
  # irt_information curve shape against mirt::testinfo on a fixed-parameter
  # graded response model. mirt works in the logistic metric, so the normal
  # ogive discrimination is scaled by 1.702 for the comparison. The logistic
  # and normal ogive information functions are proportional in shape, not
  # equal in value, so the check is on the shape of the curve and the
  # location of its peak.
  a <- c(1.2, 0.8, 1.5)
  b <- list(c(-1.2, -0.2, 0.9), c(-0.5, 0.7), 0.25)
  theta <- seq(-3, 3, length.out = 121)

  dmar <- DMAR::irt_information(
    a = rep(a, lengths(b)),
    b = unlist(b),
    item = rep(paste0("x", seq_along(a)), lengths(b)),
    theta = theta
  )

  pars <- mirt::mirt(
    data.frame(x1 = c(0:3, 0:3), x2 = c(0:2, 0:2, 0, 1),
               x3 = c(0, 1, 0, 1, 0, 1, 0, 1)),
    1, itemtype = "graded", pars = "values", verbose = FALSE
  )
  pars$est <- FALSE
  set_value <- function(pars, item_name, name, value) {
    pars$value[pars$item == item_name & pars$name == name] <- value
    pars
  }
  for (i in seq_along(a)) {
    item_name <- paste0("x", i)
    pars <- set_value(pars, item_name, "a1", 1.702 * a[i])
    for (k in seq_along(b[[i]])) {
      pars <- set_value(pars, item_name, paste0("d", k),
                        -1.702 * a[i] * b[[i]][k])
    }
  }
  fit <- mirt::mirt(
    data.frame(x1 = c(0:3, 0:3), x2 = c(0:2, 0:2, 0, 1),
               x3 = c(0, 1, 0, 1, 0, 1, 0, 1)),
    1, itemtype = "graded", pars = pars, verbose = FALSE, TOL = NaN
  )
  mirt_info <- mirt::testinfo(fit, Theta = matrix(theta))

  stopifnot(
    stats::cor(dmar$test_information, mirt_info) > 0.99,
    abs(theta[which.max(mirt_info)] -
          attr(dmar, "theta_max_information")) < 0.2
  )
})

## from tests/testthat/test-measurement_alignment.R
local({
  # measurement_alignment against sirt::invariance.alignment. sirt takes the
  # configural loading and intercept matrices directly, so it is fed exactly
  # the ones DMAR aligned. Its weights are normalized by their column sum
  # and it counts each group pair in both orders, so its criterion is DMAR's
  # divided by sum(sqrt(N_g))^2 / 2; align.pow = 0.5 is its current coding
  # of the fourth root loss of Asparouhov and Muthen.
  alignment_sim <- function(n_g, factor_mean, factor_sd, Lambda, Nu,
                            resid_sd = 0.6) {
    G <- length(n_g)
    I <- ncol(Lambda)
    do.call(rbind, lapply(seq_len(G), function(g) {
      eta <- stats::rnorm(n_g[g], factor_mean[g], factor_sd[g])
      x <- vapply(seq_len(I), function(i) {
        Nu[g, i] + Lambda[g, i] * eta + stats::rnorm(n_g[g], 0, resid_sd)
      }, numeric(n_g[g]))
      out <- as.data.frame(x)
      names(out) <- paste0("x", seq_len(I))
      out$cohort <- paste0("cohort_", g)
      out
    }))
  }
  set.seed(20)
  Lambda <- matrix(0.8, nrow = 4, ncol = 6)
  Nu <- matrix(1.0, nrow = 4, ncol = 6)
  Lambda[2, 3] <- 0.3
  Nu[3, 5] <- 1.8
  Nu[4, 1] <- 0.4
  d <- alignment_sim(c(200, 250, 180, 300), c(0, 0.3, -0.5, 0.8),
                     c(1, 1.2, 0.8, 1.1), Lambda, Nu)
  out <- DMAR::measurement_alignment(d, items = paste0("x", 1:6),
                                     group = "cohort", seed = 113)

  L0 <- attr(out, "configural_loadings")
  N0 <- attr(out, "configural_intercepts")
  n_g <- out$n
  oracle <- suppressWarnings(sirt::invariance.alignment(
    lambda = L0, nu = N0,
    wgt = matrix(sqrt(n_g), nrow = nrow(L0), ncol = ncol(L0)),
    eps = 0.01, align.pow = c(0.5, 0.5), fixed = TRUE, meth = 1
  ))

  stopifnot(
    max(abs(out$factor_mean - oracle$pars$alpha0)) < 5e-4,
    max(abs(sqrt(out$factor_variance) - oracle$pars$psi0)) < 5e-4,
    max(abs(attr(out, "aligned_loadings") - oracle$lambda.aligned)) < 5e-4,
    max(abs(attr(out, "aligned_intercepts") - oracle$nu.aligned)) < 5e-4
  )

  # The criterion itself, once sirt's weight normalization is undone. DMAR's
  # criterion evaluated at sirt's own solution equals sirt's criterion value
  # to machine precision, and DMAR's minimum is no larger than the value at
  # sirt's solution.
  rescaled <- oracle$fopt * sum(sqrt(n_g))^2 / 2
  machinery <- DMAR:::.alignment_machinery(L0, N0, n_g, 0.01, "fixed")
  at_oracle <- machinery$fn(c(oracle$pars$alpha0[-1],
                              log(oracle$pars$psi0[-1])))
  stopifnot(
    isTRUE(all.equal(at_oracle, rescaled, tolerance = 1e-12)),
    isTRUE(all.equal(attr(out, "simplicity_function"), rescaled,
                     tolerance = 1e-5)),
    attr(out, "simplicity_function") <= at_oracle + 1e-8,
    isTRUE(all.equal(unname(attr(out, "R2_total")),
                     unname(oracle$es.invariance["R2", ]),
                     tolerance = 1e-5))
  )
})
