# Live comparisons against mvtnorm (written against mvtnorm 1.4.2,
# 2026-11-01). Categorical omega called mvtnorm::pmvnorm() once per
# threshold cell through DMAR 1.0.1; DMAR now computes the bivariate
# normal probabilities itself (DMAR:::.pbinorm(), Genz's algorithm, the
# one mvtnorm applies in two dimensions) and assembles Green and Yang's
# coefficient in one vectorized pass. tests/testthat/test-pbinorm.R pins
# the mvtnorm values; these blocks re-run the comparisons live.

## from tests/testthat/test-pbinorm.R
local({
  # The bivariate normal on a grid: h and k from -3 to 3 by .5, and 17
  # correlations from -.999 to .999 on both sides of the points (.3,
  # .75, .925) where the algorithm changes its rule or its branch. The
  # largest difference was 4.0e-16.
  g <- expand.grid(h = seq(-3, 3, by = 0.5), k = seq(-3, 3, by = 0.5),
                   r = c(-0.999, -0.95, -0.9, -0.5, -0.1, 0, 0.1, 0.3, 0.5,
                         0.74, 0.76, 0.9, 0.924, 0.926, 0.95, 0.99, 0.999))
  oracle <- mapply(function(h, k, r) {
    mvtnorm::pmvnorm(upper = c(h, k), corr = matrix(c(1, r, r, 1), 2))[1]
  }, g$h, g$k, g$r)
  stopifnot(max(abs(DMAR:::.pbinorm(g$h, g$k, g$r) - oracle)) < 1e-15)
})

## from tests/testthat/test-pbinorm.R
local({
  # The 102 points whose mvtnorm values the test stores, drawn as the
  # test draws them, and the first of the values it stores (R's parser
  # can place a 17-digit decimal one unit in the last place away from
  # the double it was printed from, hence 1e-16 rather than zero).
  set.seed(113)
  sgn <- function(n) sample(c(-1, 1), n, replace = TRUE)
  r <- c(runif(24, -0.3, 0.3), sgn(24) * runif(24, 0.3, 0.75),
         sgn(24) * runif(24, 0.75, 0.925), sgn(24) * runif(24, 0.925, 0.9999),
         c(-0.925, -0.75, -0.3, 0.3, 0.75, 0.925))
  h <- c(runif(96, -3, 3), rep(0.4, 6))
  k <- c(runif(96, -3, 3), rep(-1.1, 6))
  oracle <- mapply(function(h, k, r) {
    mvtnorm::pmvnorm(upper = c(h, k), corr = matrix(c(1, r, r, 1), 2))[1]
  }, h, k, r)
  stopifnot(max(abs(DMAR:::.pbinorm(h, k, r) - oracle)) < 1e-15)
  stored <- c(0.00042464521235648843, 0.24019082896897503, 0.10058048057338051)
  stopifnot(max(abs(oracle[1:3] - stored)) < 1e-16)
})

## from tests/testthat/test-pbinorm.R
local({
  # Green and Yang's (2009, Table 1) ten population configurations, by
  # the loops of DMAR 1.0.1 with mvtnorm in the inner call, against
  # DMAR's vectorized assembly. The largest difference was 8.9e-16.
  gy_mvtnorm <- function(thresh, rho_num, rho_den) {
    cell <- function(a, b, r) {
      mvtnorm::pmvnorm(upper = c(a, b), corr = matrix(c(1, r, r, 1), 2))[1]
    }
    num <- 0
    den <- 0
    for (j in seq_along(thresh)) for (jp in seq_along(thresh)) {
      t1 <- thresh[[j]]
      t2 <- thresh[[jp]]
      p_num <- 0
      p_den <- 0
      for (a in t1) for (b in t2) {
        p_num <- p_num + cell(a, b, rho_num[j, jp])
        p_den <- p_den + cell(a, b, rho_den[j, jp])
      }
      product <- sum(stats::pnorm(t1)) * sum(stats::pnorm(t2))
      num <- num + (p_num - product)
      den <- den + (p_den - product)
    }
    num / den
  }
  thr <- function(P) stats::qnorm(cumsum(P)[-length(P)])
  both <- function(P_list, lambda) {
    rho_num <- outer(lambda, lambda)
    rho_den <- rho_num
    diag(rho_den) <- 1
    th <- lapply(P_list, thr)
    c(DMAR:::.gy_omega_from_components(th, rho_num, rho_den),
      gy_mvtnorm(th, rho_num, rho_den))
  }
  u2 <- c(.5, .5); neg2 <- c(.10, .90); pos2 <- rev(neg2)
  u5 <- rep(.2, 5); bell <- c(.09, .15, .51, .15, .09)
  matched <- c(.485, .01, .01, .01, .485)
  smooth_neg <- c(.011, .027, .047, .118, .797); smooth_pos <- rev(smooth_neg)
  l8 <- rep(.8, 4); l85 <- c(.8, .8, .5, .5)
  pairs <- rbind(
    both(rep(list(u2), 4), l8), both(rep(list(u5), 4), l8),
    both(rep(list(bell), 4), l8), both(rep(list(matched), 4), l8),
    both(rep(list(neg2), 4), l8), both(list(neg2, neg2, pos2, pos2), l8),
    both(rep(list(smooth_neg), 4), l8),
    both(list(smooth_neg, smooth_neg, smooth_pos, smooth_pos), l8),
    both(rep(list(u2), 4), l85), both(rep(list(u5), 4), l85))
  stopifnot(max(abs(pairs[, 1] - pairs[, 2])) < 2e-15)
})

## from tests/testthat/test-reliability_omega_categorical.R
local({
  # The model implied denominator of reliability_omega_categorical()
  # against Eq. 21 assembled with mvtnorm from a separate lavaan fit
  # (the comparison the test now makes by numerical integration), and
  # the observed denominator on six Holzinger and Swineford tests cut
  # into five categories against the mvtnorm loops on the same fit's
  # pieces.
  set.seed(113)
  N <- 400
  J <- 5
  loadings <- rep(0.7, J)
  eta <- rnorm(N)
  latent <- outer(eta, loadings) +
    matrix(rnorm(N * J), N, J) %*% diag(sqrt(1 - loadings^2))
  items <- apply(latent, 2, function(x)
    as.integer(cut(x, breaks = c(-Inf, -1.5, -0.5, 0.5, 1.5, Inf),
                   labels = FALSE)))
  colnames(items) <- paste0("y", seq_len(J))
  r <- DMAR::reliability_omega_categorical(data = items, ci_method = "none",
                                           denominator = "model_implied")
  est <- r$value[r$term == "estimate"]
  d <- as.data.frame(items)
  d[] <- lapply(d, ordered)
  fit <- lavaan::cfa(paste("f =~", paste(colnames(d), collapse = " + ")),
                     data = d, ordered = colnames(d), std.lv = TRUE,
                     se = "none")
  lam <- lavaan::lavInspect(fit, "est")$lambda[, 1]
  implied <- stats::cov2cor(lavaan::lavInspect(fit, "implied")$cov)
  tau <- lavaan::lavInspect(fit, "est")$tau
  th <- split(as.numeric(tau), sub("\\|.*$", "", rownames(tau)))[colnames(d)]
  cov_cat <- function(t1, t2, r) {
    p2 <- 0
    for (a in t1) for (b in t2) {
      p2 <- p2 + mvtnorm::pmvnorm(upper = c(a, b),
                                  corr = matrix(c(1, r, r, 1), 2))[1]
    }
    p2 - sum(stats::pnorm(t1)) * sum(stats::pnorm(t2))
  }
  num <- den <- 0
  for (j in seq_along(lam)) for (k in seq_along(lam)) {
    num <- num + cov_cat(th[[j]], th[[k]], lam[j] * lam[k])
    den <- den + cov_cat(th[[j]], th[[k]], implied[j, k])
  }
  stopifnot(isTRUE(all.equal(est, num / den, tolerance = 1e-6)))

  data(holzinger_swineford, package = "DMAR")
  tests <- c("t1_visual_perception", "t2_cubes", "t3_paper_form_board",
             "t4_lozenges", "t6_paragraph_comprehension", "t7_sentence")
  six <- as.data.frame(lapply(holzinger_swineford[, tests], function(v) {
    as.integer(cut(v, stats::quantile(v, 0:5 / 5), include.lowest = TRUE))
  }))
  parts <- DMAR:::.omega_c_parts(six)
  loops <- function(rho_den) {
    num <- den <- 0
    for (j in seq_along(parts$thresh)) for (k in seq_along(parts$thresh)) {
      num <- num + cov_cat(parts$thresh[[j]], parts$thresh[[k]],
                           parts$rho_num[j, k])
      den <- den + cov_cat(parts$thresh[[j]], parts$thresh[[k]], rho_den[j, k])
    }
    num / den
  }
  observed <- DMAR::reliability_omega_categorical(six, ci_method = "none")
  stopifnot(abs(observed$value[observed$term == "estimate"] -
                  loops(parts$rho_sat)) < 2e-15)
  stopifnot(abs(DMAR:::.gy_omega_from_components(parts$thresh, parts$rho_num,
                                                 parts$rho_model) -
                  loops(parts$rho_model)) < 2e-15)
})
