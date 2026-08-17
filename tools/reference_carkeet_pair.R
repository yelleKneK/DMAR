# Verification of DMAR loa() against Carkeet (2015) OVS 92(3) e71-e80.
suppressPackageStartupMessages(library(DMAR))

opt <- options(warn = -1)  # silence pnt{final} precision warnings in output

## ---------------------------------------------------------------------------
## 1. Reproduce Carkeet Table 1 (individual LoA coefficients), p. e74.
##    Definition (display eq., p. e74):
##      Pr{ (noncentral t, df = nu, ncp = z_{0.975} * sqrt(n)) <= k*sqrt(n) } = F
##    so c_F = qt(F, df = nu, ncp = z*sqrt(n)) / sqrt(n),  n = nu + 1.
cF_qt <- function(nu, F, P = 0.975) {
  n <- nu + 1
  qt(F, df = nu, ncp = qnorm(P) * sqrt(n)) / sqrt(n)
}

## Independent implementation avoiding qt(ncp): Pr(T <= q) = E_S[ Phi(q*S - delta) ],
## S = chi_nu / sqrt(nu); solve for the F-quantile by uniroot.
cF_int <- function(nu, F, P = 0.975) {
  n <- nu + 1
  delta <- qnorm(P) * sqrt(n)
  pT <- function(q) {
    integrand <- function(x) pnorm(q * x / sqrt(nu) - delta) * dchisq(x^2, nu) * 2 * x
    # chi_nu mass lives near sqrt(nu) with sd ~ 1/sqrt(2); window it for large nu
    integrate(integrand, max(1e-10, sqrt(nu) - 14), sqrt(nu) + 14,
              rel.tol = 1e-12)$value
  }
  q <- uniroot(function(q) pT(q) - F, c(-10 * sqrt(n), 200 * sqrt(n)),
               tol = 1e-12)$root
  q / sqrt(n)
}

cat("== Table 1 reproduction (individual), paper values in brackets ==\n")
tab1 <- rbind(
  c(nu = 4,    F = 0.025, paper = 0.9232),
  c(nu = 4,    F = 0.975, paper = 5.9749),
  c(nu = 16,   F = 0.025, paper = 1.3150),
  c(nu = 16,   F = 0.500, paper = 1.9982),
  c(nu = 16,   F = 0.975, paper = 3.1483),
  c(nu = 50,   F = 0.050, paper = 1.6144),
  c(nu = 50,   F = 0.950, paper = 2.4268),
  c(nu = 1000, F = 0.025, paper = 1.8578),
  c(nu = 1000, F = 0.975, paper = 2.0700)
)
for (i in seq_len(nrow(tab1))) {
  nu <- tab1[i, "nu"]; F <- tab1[i, "F"]; pv <- tab1[i, "paper"]
  v1 <- cF_qt(nu, F); v2 <- cF_int(nu, F)
  cat(sprintf("nu=%4d F=%.3f  qt: %.6f  integral: %.6f  [paper %.4f]  diff(qt-paper)=%+.1e\n",
              nu, F, v1, v2, pv, v1 - pv))
}

## ---------------------------------------------------------------------------
## 2. Reproduce Carkeet Table 2 (pair / two-sided tolerance factors), p. e75.
##    Odeh (1978) integral:
##      (2*sqrt(n)/sqrt(2*pi)) * Int_0^inf Pr(chisq_nu > nu*r(xb)^2/kt^2) *
##                               exp(-n*xb^2/2) dxb = F
##    where r(xb) solves Phi(xb + r) - Phi(xb - r) = P.
kt_pair <- function(nu, F, P = 0.95) {
  n <- nu + 1
  r_of <- function(xb) {
    uniroot(function(r) pnorm(xb + r) - pnorm(xb - r) - P,
            c(1e-8, xb + 50), tol = 1e-13)$root
  }
  xb_max <- sqrt(220 / n)  # exp(-n*xb^2/2) < 1e-47 beyond this
  conf <- function(kt) {
    integrand <- function(xb) {
      sapply(xb, function(u)
        pchisq(nu * r_of(u)^2 / kt^2, nu, lower.tail = FALSE) *
          exp(-n * u^2 / 2))
    }
    2 * sqrt(n) / sqrt(2 * pi) *
      integrate(integrand, 0, xb_max, rel.tol = 1e-10)$value
  }
  uniroot(function(kt) conf(kt) - F, c(0.2, 400), tol = 1e-10)$root
}

cat("\n== Table 2 reproduction (pair), paper values in brackets ==\n")
tab2 <- rbind(
  c(nu = 4,  F = 0.025, paper = 1.2397),
  c(nu = 4,  F = 0.975, paper = 6.1569),
  c(nu = 16, F = 0.025, paper = 1.4900),
  c(nu = 16, F = 0.975, paper = 3.0824),
  c(nu = 19, F = 0.500, paper = 2.0410)
)
for (i in seq_len(nrow(tab2))) {
  nu <- tab2[i, "nu"]; F <- tab2[i, "F"]; pv <- tab2[i, "paper"]
  v <- kt_pair(nu, F)
  cat(sprintf("nu=%4d F=%.3f  computed: %.6f  [paper %.4f]  diff=%+.1e\n",
              nu, F, v, pv, v - pv))
}

## ---------------------------------------------------------------------------
## 3. Bland-Altman 1986 worked example (p. e74): n = 17, dbar = -2.1, s = 38.8.
##    Paper: upper LoA CI [48.9, 120.0]; lower LoA CI [-124.2, -53.1].
##    Construct differences with those exact moments and run DMAR loa().
n <- 17
z <- scale(rnorm(n))            # mean 0, sd 1 exactly (scale divides by sd)
d <- as.numeric(z) * 38.8 - 2.1 # mean -2.1, sd 38.8 exactly
stopifnot(abs(mean(d) + 2.1) < 1e-12, abs(sd(d) - 38.8) < 1e-12)
res <- loa(rep(0, n), d)
print(res)
v <- function(t) res$value[res$term == t]
cat(sprintf("\nloa() upper LoA %.4f, CI [%.4f, %.4f]  (paper: 73.9, [48.9, 120.0])\n",
            v("loa_upper"), v("loa_upper_lower_limit"), v("loa_upper_upper_limit")))
cat(sprintf("loa() lower LoA %.4f, CI [%.4f, %.4f]  (paper: -78.1, [-124.2, -53.1])\n",
            v("loa_lower"), v("loa_lower_lower_limit"), v("loa_lower_upper_limit")))
## Same from Table-1-style coefficients at full precision:
cat(sprintf("check via c_F: upper CI [%.4f, %.4f]\n",
            -2.1 + cF_qt(16, 0.025) * 38.8, -2.1 + cF_qt(16, 0.975) * 38.8))
## Pair-method CI (what DMAR does NOT compute), p. e76: -2.1 +/- 57.8 ; -2.1 +/- 119.6
cat(sprintf("pair method (not in DMAR): inner +/- %.4f (paper 57.8), outer +/- %.4f (paper 119.6)\n",
            kt_pair(16, 0.025) * 38.8, kt_pair(16, 0.975) * 38.8))

## ---------------------------------------------------------------------------
## 4. Equivalence of loa() to the exact construction over a grid.
cat("\n== loa() vs direct exact construction over a grid ==\n")
set.seed(7)
max_abs <- 0
for (nn in c(5, 17, 40)) for (cov in c(0.90, 0.95)) for (cl in c(0.90, 0.95)) {
  x <- rnorm(nn); y <- x + rnorm(nn, 0.3, 1.7)
  r <- loa(x, y, coverage = cov, conf_level = cl)
  vv <- function(t) r$value[r$term == t]
  d  <- y - x; dm <- mean(d); ds <- sd(d)
  k  <- qnorm((1 + cov) / 2); del <- k * sqrt(nn); al <- 1 - cl
  up_lo <- dm + ds / sqrt(nn) * qt(al / 2,      nn - 1, ncp = del)
  up_hi <- dm + ds / sqrt(nn) * qt(1 - al / 2,  nn - 1, ncp = del)
  lo_lo <- dm + ds / sqrt(nn) * qt(al / 2,      nn - 1, ncp = -del)
  lo_hi <- dm + ds / sqrt(nn) * qt(1 - al / 2,  nn - 1, ncp = -del)
  m <- max(abs(c(vv("loa_upper_lower_limit") - up_lo,
                 vv("loa_upper_upper_limit") - up_hi,
                 vv("loa_lower_lower_limit") - lo_lo,
                 vv("loa_lower_upper_limit") - lo_hi)))
  max_abs <- max(max_abs, m)
}
cat(sprintf("max |loa() - direct exact| over grid: %.2e\n", max_abs))

## ---------------------------------------------------------------------------
## 5. Monte Carlo: n = 20, d ~ N(0,1), 20000 reps.
##    True upper limit theta = 1.959964. Coverage of:
##    (a) loa()'s CI (exact individual method, per-limit),
##    (b) Bland-Altman 1999 approximate CI: LoA +/- t_{0.975,19}*sqrt(2.92)*s/sqrt(n),
##    (c) pair-method interval [dbar + ct_.025 s, dbar + ct_.975 s] for the same theta,
##        plus its intended simultaneous tolerance-content coverage.
cat("\n== Monte Carlo, n = 20, 20000 reps ==\n")
set.seed(113)
nmc <- 20000; n <- 20
theta <- qnorm(0.975)                      # true upper limit, mu=0 sigma=1
## constants (verified above to equal what loa() computes):
t_lo <- qt(0.025, n - 1, ncp = qnorm(0.975) * sqrt(n))
t_hi <- qt(0.975, n - 1, ncp = qnorm(0.975) * sqrt(n))
tBA  <- qt(0.975, n - 1) * sqrt(2.92)
ct_lo <- kt_pair(n - 1, 0.025); ct_hi <- kt_pair(n - 1, 0.975)

## spot-check that loa() reproduces the vectorized exact CI on 25 draws:
for (i in 1:25) {
  d <- rnorm(n)
  r <- loa(rep(0, n), d)
  vv <- function(t) r$value[r$term == t]
  stopifnot(abs(vv("loa_upper_lower_limit") - (mean(d) + sd(d)/sqrt(n)*t_lo)) < 1e-10,
            abs(vv("loa_upper_upper_limit") - (mean(d) + sd(d)/sqrt(n)*t_hi)) < 1e-10)
}
cat("spot-check: loa() CI == vectorized exact CI on 25 draws: OK\n")

dm <- ds <- numeric(nmc)
for (i in seq_len(nmc)) { d <- rnorm(n); dm[i] <- mean(d); ds[i] <- sd(d) }

exact_lo <- dm + ds / sqrt(n) * t_lo
exact_hi <- dm + ds / sqrt(n) * t_hi
cov_exact <- mean(exact_lo <= theta & theta <= exact_hi)
miss_lo_e <- mean(theta < exact_lo); miss_hi_e <- mean(theta > exact_hi)

ba_c <- dm + theta * ds                    # sample upper LoA
ba_lo <- ba_c - tBA * ds / sqrt(n)
ba_hi <- ba_c + tBA * ds / sqrt(n)
cov_ba <- mean(ba_lo <= theta & theta <= ba_hi)
miss_lo_b <- mean(theta < ba_lo); miss_hi_b <- mean(theta > ba_hi)

pr_lo <- dm + ct_lo * ds; pr_hi <- dm + ct_hi * ds
cov_pair_theta <- mean(pr_lo <= theta & theta <= pr_hi)
## intended pair guarantee: content of dbar +/- ct s brackets 0.95
content <- function(kk) pnorm(dm + kk * ds) - pnorm(dm - kk * ds)
cov_pair_content <- mean(content(ct_lo) <= 0.95 & content(ct_hi) >= 0.95)

se <- function(p) sqrt(p * (1 - p) / nmc)
cat(sprintf("loa() / exact individual CI coverage of mu+1.96sigma: %.4f (MC SE %.4f)\n",
            cov_exact, se(cov_exact)))
cat(sprintf("   miss low %.4f, miss high %.4f (nominal 0.025 each)\n", miss_lo_e, miss_hi_e))
cat(sprintf("Bland-Altman approximate CI coverage:                 %.4f (MC SE %.4f)\n",
            cov_ba, se(cov_ba)))
cat(sprintf("   miss low %.4f, miss high %.4f\n", miss_lo_b, miss_hi_b))
cat(sprintf("pair interval coverage of the same theta:             %.4f\n", cov_pair_theta))
cat(sprintf("pair intended simultaneous content guarantee:         %.4f (nominal 0.95)\n",
            cov_pair_content))

options(opt)
