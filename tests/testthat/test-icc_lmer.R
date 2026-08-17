test_that("icc_lmer() returns documented rows", {
  skip_if_not_installed("lme4")
  set.seed(113)
  n_grp <- 20; n_per <- 6
  grp <- factor(rep(1:n_grp, each = n_per))
  y   <- rnorm(n_grp * n_per) + rep(rnorm(n_grp, 0, 0.7), each = n_per)
  fit <- lme4::lmer(y ~ 1 + (1 | grp), data = data.frame(y, grp))
  res <- icc_lmer(fit)
  expect_setequal(res$term,
                  c("icc", "sigma2_between", "sigma2_within", "total_variance",
                    "lower_limit", "upper_limit", "n_clusters",
                    "average_cluster_size"))
})

test_that("icc_lmer() ICC bounded in [0, 1]", {
  skip_if_not_installed("lme4")
  set.seed(113)
  grp <- factor(rep(1:20, each = 6))
  y   <- rnorm(120) + rep(rnorm(20, 0, 0.5), each = 6)
  fit <- lme4::lmer(y ~ 1 + (1 | grp), data = data.frame(y, grp))
  rho <- icc_lmer(fit)$value[1]
  expect_gte(rho, 0); expect_lte(rho, 1)
})

test_that("icc_lmer() rejects non-lmer input", {
  expect_error(icc_lmer(lm(mpg ~ 1, data = mtcars)), "lme4::lmer")
})

test_that("icc_lmer() handles interaction-term groupings like (1 | A:B)", {
  skip_if_not_installed("lme4")
  d <- transform(lme4::sleepstudy, g = factor(rep(1:2, 90)))
  f <- suppressMessages(lme4::lmer(Reaction ~ Days + (1 | g:Subject),
                                   data = d))
  res <- icc_lmer(f)
  # ICC equals the by-hand VarCorr ratio.
  vc <- as.data.frame(lme4::VarCorr(f))
  rho_hand <- vc$vcov[vc$grp == "g:Subject"] / sum(vc$vcov)
  expect_equal(res$value[res$term == "icc"], rho_hand)
  # The grouping has the right number of clusters (36 g x Subject cells).
  flist <- lme4::getME(f, "flist")
  expect_equal(res$value[res$term == "n_clusters"],
               nlevels(flist[["g:Subject"]]))
})

test_that("icc_lmer() CI follows the documented Bonett variance k / (2 (k - 1) (n - 2))", {
  skip_if_not_installed("lme4")
  set.seed(113)
  n_grp <- 20; n_per <- 6
  grp <- factor(rep(1:n_grp, each = n_per))
  y   <- rnorm(n_grp * n_per) + rep(rnorm(n_grp, 0, 0.7), each = n_per)
  fit <- lme4::lmer(y ~ 1 + (1 | grp), data = data.frame(y, grp))
  res <- icc_lmer(fit)
  v <- function(t) res$value[res$term == t]
  rho <- v("icc"); k <- v("average_cluster_size"); n_cl <- v("n_clusters")
  L  <- 0.5 * log((1 + (k - 1) * rho) / (1 - rho))
  # Documented Bonett (2002) variance of L. The page previously omitted
  # the k / (k - 1) factor the code applies; the omission understates
  # the SE (by sqrt(2) at k = 2) and moves these limits at the third
  # decimal, so this recomputation discriminates the two formulas.
  se <- sqrt(k / (2 * (k - 1) * (n_cl - 2)))
  z  <- stats::qnorm(0.975)
  inv_L <- function(LL) (exp(2 * LL) - 1) / (exp(2 * LL) + k - 1)
  expect_equal(v("lower_limit"), max(0, inv_L(L - z * se)), tolerance = 1e-10)
  expect_equal(v("upper_limit"), min(1, inv_L(L + z * se)), tolerance = 1e-10)
})
