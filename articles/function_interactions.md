# How DMAR's Functions Fit Together: A Deep Dive

## Why This Vignette Exists

DMAR has grown to encompass several large families of functions: effect
size point estimates, exact-noncentral distribution confidence
intervals, asymptotic-variance utilities, design-stage expected-value
calculators, accuracy in parameter estimation (AIPE) sample size
planners, reliability coefficients, and agreement / multivariate
measures. Each function is documented in isolation, but the package’s
value is in how the functions *compose*. This vignette walks through the
package’s architecture from the user’s perspective and shows the five
high-traffic composition patterns.

## The Four Building Blocks

For any single population quantity $`\theta`$ the package supports four
operations:

           point estimate      sample size                noncentral CI            asymptotic
           /  variance         planning                                            building blocks
           ┌──────────────┐    ┌──────────────────┐      ┌───────────────────┐    ┌─────────────────┐
    \(\theta\) → │ smd, omega², r│ →  │ ss_aipe_smd,     │ →    │ ci_smd, ci_r,     │ ←  │ var_smd, var_r, │
           │ R², icc, ...  │    │ ss_aipe_R2, ...  │      │ ci_omega_squared, │    │ var_omega_squared│
           └──────────────┘    └──────────────────┘      │ ci_R2, ...        │    │ ...             │
                  ↑                       ↑              └───────────────────┘            ↑
                  │                       │                       │                       │
                  │                       │                       │                       │
                  │                ┌──────┴─────────────┐         │                       │
                  │                │ expected_smd,      │         │                       │
                  │                │ expected_r,        │ ────────┘                       │
                  │                │ expected_R2 ...    │                                 │
                  │                └────────────────────┘                                 │
                  │                                                                       │
                  └──────────────── transforms of each other ─────────────────────────────┘
                  (smd ↔ cles ↔ proportion_of_superiority ↔ nnt_from_smd; r ↔ partial_r)

The four columns are tightly coupled: every AIPE planner is built from a
variance utility; every CI is built from a noncentral distribution;
every expected-value function is the design-stage analog of an
estimator. Once you internalize the layout, asking “what is the
recommended workflow for quantity X” reduces to asking “which row of X’s
variance / CI / expected value / AIPE planner do I need.”

## 1. Variance Utility → AIPE Planner

The most pervasive composition is `var_*` feeding `ss_aipe_*`. For any
quantity $`\theta`$ with a known asymptotic variance
$`\mathrm{Var}(\hat\theta \mid \theta, n)`$, the AIPE planner solves for
the smallest $`n`$ such that

``` math
 z_{1-\alpha/2} \cdot \sqrt{\mathrm{Var}(\hat\theta \mid \theta, n)} \le w/2, 
```

where $`w`$ is the target full width of the CI. The AIPE framework,
planning for a sufficiently narrow confidence interval rather than for
power, is developed in Kelley and Maxwell (2003) and Maxwell, Kelley,
and Rausch (2008); see also Maxwell, Delaney, and Kelley (2027) for the
model comparison treatment. Concretely:

``` r

# SMD (Cohen's d):     var_smd        → ss_aipe_smd
# Pearson r (J = 0):   var_r          → ss_aipe_partial_r (J = 1, simple case)
# Partial r:           var_partial_r  → ss_aipe_partial_r
# Semipartial r:       var_semipartial_r → ss_aipe_semipartial_r
# Squared multi R:     var_R2         → ss_aipe_R2
# Cliff's delta:       (DeLong bound) → ss_aipe_cliff_delta
# omega²:              var_omega_squared → ss_aipe_omega_squared (noncentral F)
# Cronbach alpha:      var_alpha      → ss_aipe_reliability
# ICC:                 var_icc        → ss_aipe_icc
# CV:                  var_cv         → ss_aipe_cv
# Indirect effect ab:  var_indirect_effect → ss_aipe_indirect_effect
```

This pairing is the package’s most common usage pattern. The variance
utility is the design-stage primitive; the AIPE planner is the inverse,
returning the $`n`$ that achieves a target precision.

### Worked Example: Pearson r

We plan a study where we expect $`\rho \approx 0.30`$ and want a 95% CI
of full width 0.20.

``` r

# Step 1: ask the variance utility for the sampling variance at a candidate n.
var_r(rho = 0.30, n = 100)
```

| term         | value   |
|:-------------|:--------|
| var_r_normal | 0.00836 |
| var_fisher_z | 0.0103  |

``` r


# Step 2: the AIPE planner inverts the same family of formulas to give
# the n needed for the target width. DMAR has no zero-order Pearson-r
# planner; the partial-r planner at J = 1 is the simple-correlation case.
ss_aipe_partial_r(rho = 0.30, J = 1, width = 0.20)
```

| term           | value |
|:---------------|:------|
| necessary_N    | 321   |
| expected_width | 0.2   |
| rho            | 0.3   |
| J              | 1     |
| width_target   | 0.2   |
| conf_level     | 0.95  |

Confidence level: 95%

Both sides use the same normal-theory $`(1 - \rho^2)^2`$ variance family
(Fisher, 1915; Olkin & Finn, 1995), differing only by an off-by-one in
the denominator degrees of freedom ($`n - 1`$ for
[`var_r()`](https://yelleknek.github.io/DMAR/reference/var_r.md),
$`n - J - 1`$ for the planner), so the planner’s output is consistent
with what
[`var_r()`](https://yelleknek.github.io/DMAR/reference/var_r.md) returns
at that recommended $`n`$.

## 2. Expected-Value Function → Design-Stage Substitution

For estimators with non-trivial bias, the package provides
`expected_*()` functions that report
$`\mathrm{E}[\hat\theta \mid \theta, n]`$ and the size of the bias.
These are useful in two situations:

1.  **Reporting:** at the analysis stage, the difference between
    $`\hat\theta`$ and $`\mathrm{E}[\hat\theta]`$ tells the analyst how
    much of the observed value is bias.
2.  **Design-stage planning:** when a sample size planner takes
    $`\theta`$ as input but the realized data will yield a biased sample
    $`\hat\theta`$, substituting $`\mathrm{E}[\hat\theta]`$ in place of
    $`\theta`$ corrects the planner’s prediction.

The package’s expected-value functions are:

| Quantity | Bias | Helper |
|----|----|----|
| Pearson $`r`$ | downward (small $`n`$) | [`expected_r()`](https://yelleknek.github.io/DMAR/reference/expected_r.md) |
| Partial $`r`$ | downward (grows with $`J`$) | [`expected_partial_r()`](https://yelleknek.github.io/DMAR/reference/expected_partial_r.md) |
| Squared multi $`R^2`$ | upward (always) | [`expected_R2()`](https://yelleknek.github.io/DMAR/reference/expected_R2.md) |
| Cohen’s $`d`$ | upward (small $`n`$) | [`expected_smd()`](https://yelleknek.github.io/DMAR/reference/expected_smd.md) |

### Worked Example: The Standardized Mean Difference

``` r

# Population delta = 0.4, planned n = 40 per group.
# Naive planning assumes the observed d will average 0.40; in fact
# the sample d will average slightly larger.
expected_smd(delta = 0.4, n_1 = 40)
```

| delta | n_1 | n_2 | expected_smd | bias   | j_correction |
|:------|:----|:----|:-------------|:-------|:-------------|
| 0.4   | 40  | 40  | 0.404        | 0.0039 | 0.99         |

The reported `j_correction` is the Hedges (1981) factor
$`J(\mathit{df})`$ that maps Cohen’s $`d`$ to Hedges’ $`g`$. The same
factor is applied internally by `smd(..., unbiased = TRUE)`, which
returns the bias-corrected estimate $`g = J(\mathit{df})\,d`$, so the
same constant runs throughout the $`d`$-family of functions for
consistency.

## 3. Effect Size Transforms: smd ↔︎ cles ↔︎ proportion_of_superiority ↔︎ nnt

A single standardized mean difference (Cohen’s $`d`$) can be reported on
five common scales without losing information:

| Scale | Formula | Function in DMAR |
|----|----|----|
| Standardized mean difference | $`d`$ | [`smd()`](https://yelleknek.github.io/DMAR/reference/smd.md), [`ci_smd()`](https://yelleknek.github.io/DMAR/reference/ci_smd.md) |
| Proportion of superiority (Cohen’s $`U_3`$) | $`\Phi(d)`$ | [`proportion_of_superiority()`](https://yelleknek.github.io/DMAR/reference/proportion_of_superiority.md) |
| Common language (independent groups) | $`\Phi(d / \sqrt 2)`$ | [`cles()`](https://yelleknek.github.io/DMAR/reference/cles.md) |
| Success-rate difference | $`2 \Phi(d / \sqrt 2) - 1`$ | (intermediate) |
| Number needed to treat | $`1 / [2 \Phi(d / \sqrt 2) - 1]`$ | [`nnt_from_smd()`](https://yelleknek.github.io/DMAR/reference/nnt_from_smd.md) |

Each downstream function accepts either a point estimate of $`d`$ or the
same plus group sizes for an exact noncentral $`t`$ CI, and propagates
the CI bounds through the monotone transform. Because the transforms are
monotone, the coverage probability is preserved exactly (no extra
approximation is introduced).

### Walked Example at Three Reference d Values

``` r

for (d in c(0.2, 0.5, 0.8)) {
  cat(sprintf(
    "d = %.1f  ->  PS = %.2f   CLES = %.2f   NNT = %.1f\n",
    d,
    proportion_of_superiority(d)$value[2],
    cles(d)$value[2],
    nnt_from_smd(d)$value[3]
  ))
}
#> d = 0.2  ->  PS = 0.58   CLES = 0.56   NNT = 8.9
#> d = 0.5  ->  PS = 0.69   CLES = 0.64   NNT = 3.6
#> d = 0.8  ->  PS = 0.79   CLES = 0.71   NNT = 2.3
```

For an applied report, the same effect can be communicated on whichever
scale is most natural for the audience without recomputing anything.

## 4. Ordinal Alternatives: cliff_delta, vargha_delaney_A, probability_of_superiority_paired

When the dependent variable is ordinal, skewed, or you want a
nonparametric effect size, the same role is played by U-statistic
estimators:

| Independent groups | Function |
|----|----|
| Probability of superiority | [`cles()`](https://yelleknek.github.io/DMAR/reference/cles.md) (under bivariate normality) |
| Stochastic dominance | [`cliff_delta()`](https://yelleknek.github.io/DMAR/reference/cliff_delta.md) (ordinal, U-statistic) |
|  | Vargha-Delaney $`A = (\delta + 1)/2`$, [`vargha_delaney_A()`](https://yelleknek.github.io/DMAR/reference/vargha_delaney_A.md) |
| Paired data | [`probability_of_superiority_paired()`](https://yelleknek.github.io/DMAR/reference/probability_of_superiority_paired.md) |

Each returns an analytic CI;
[`cliff_delta()`](https://yelleknek.github.io/DMAR/reference/cliff_delta.md)
uses the Cliff (1996) U-statistic variance, and
[`probability_of_superiority_paired()`](https://yelleknek.github.io/DMAR/reference/probability_of_superiority_paired.md)
uses the Brunner-Munzel (2000) within-pair variance, all with a
Fisher-style \$\arctanh\$ transformation so the bounds stay in
$`[-1, 1]`$ (for $`\delta`$) or $`[0, 1]`$ (for $`P_S`$).

## 5. Reliability Composition

The reliability family follows the same four-block pattern but with one
extra wrinkle: alpha is the workhorse, omega is the modern default, and
$`H`$ is the upper bound that sets the ceiling of what is achievable
with the same indicators:

``` r

# Point estimate:      reliability_alpha(), reliability_omega(), reliability_H()
# CI:                  reliability() (general wrapper)
# Asymptotic variance: var_alpha()
# Design-stage AIPE:   ss_aipe_reliability()
```

$`H`$ is uniformly $`\ge \omega \ge \alpha`$ for unidimensional
indicators; if
[`reliability_H()`](https://yelleknek.github.io/DMAR/reference/reliability_H.md)
is far larger than
[`reliability_alpha()`](https://yelleknek.github.io/DMAR/reference/reliability_alpha.md),
the equally weighted composite is wasting reliability that an optimally
weighted composite would capture.

## 6. The Mediation Pipeline

For a simple three-variable mediation model $`X \to M \to Y`$ with
indirect effect $`ab`$:

``` r

# Asymptotic variance: var_indirect_effect(a, b, var_a, var_b, cov_ab)
#   - returns Sobel, Aroian, Goodman, and second-order delta variances.
# AIPE planning:       ss_aipe_indirect_effect(a, b, width, method)
#   - method = "sobel"      → asymptotic SE
#   - method = "monte_carlo"→ Tofighi & MacKinnon (2011) MC CI
```

The two methods correspond to the two main approaches in the mediation
literature;
[`var_indirect_effect()`](https://yelleknek.github.io/DMAR/reference/var_indirect_effect.md)
makes the underlying SE formulas explicit so the user can see what the
planner is solving against.

## 7. Multivariate and Agreement Utilities

The package’s measurement family is small but covers the common needs:

- [`lin_ccc()`](https://yelleknek.github.io/DMAR/reference/lin_ccc.md):
  Lin’s (1989) concordance correlation, with the King-Chinchilli (2001)
  skewness-corrected CI as the default.
- [`loa()`](https://yelleknek.github.io/DMAR/reference/loa.md):
  Bland-Altman 95% limits of agreement with the Carkeet (2015) exact
  noncentral *t* CIs on the limits.
- [`variance_components_mls()`](https://yelleknek.github.io/DMAR/reference/variance_components_mls.md):
  Burdick-Graybill (1992) modified-large- sample CIs on variance
  components, the standard tool in generalizability theory (Brennan,
  2001).
- [`ci_eigenvalue()`](https://yelleknek.github.io/DMAR/reference/ci_eigenvalue.md):
  log-scale CI on a sample eigenvalue under the Anderson (2003)
  asymptotic distribution.
- [`procrustes_phi()`](https://yelleknek.github.io/DMAR/reference/procrustes_phi.md):
  Tucker’s congruence coefficient with permutation significance test
  (Lorenzo-Seva & ten Berge, 2006).

These are stand-alone utilities; they do not feed into AIPE planners
yet, although the variance components implied by
[`variance_components_mls()`](https://yelleknek.github.io/DMAR/reference/variance_components_mls.md)
do compose with
[`var_icc()`](https://yelleknek.github.io/DMAR/reference/var_icc.md) and
[`ss_aipe_icc()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_icc.md)
for generalizability- theory planning.

## Argument-Name Conventions

When using multiple DMAR functions in a single pipeline, the
argument-name conventions are worth knowing:

| Quantity                       | Argument name(s)                           |
|--------------------------------|--------------------------------------------|
| Population correlation         | `rho` (Greek; new code)                    |
| Sample correlation             | `r` (Roman; existing code)                 |
| Population SMD                 | `delta`                                    |
| Sample SMD                     | `smd`                                      |
| Group sample sizes             | `n_1`, `n_2`                               |
| Total sample size              | `n` (single-sample), `N` (some legacy)     |
| Number of controls (partial r) | `J`                                        |
| Number of items                | `p_items`                                  |
| Number of predictors           | `p`                                        |
| Number of raters (ICC)         | `k`                                        |
| Effect df / error df           | `df_effect`, `df_error`                    |
| Target CI width                | `width`                                    |
| Which-width specifier          | `which_width` ∈ {“Full”, “Lower”, “Upper”} |
| Confidence level               | `conf_level`                               |
| Per-tail alphas                | `alpha_lower`, `alpha_upper`               |
| Coverage assurance probability | `assurance` (in `(0.5, 1)`)                |

Pattern: lower-case Greek letters denote population values (`rho`,
`delta`, `alpha`); Roman letters denote sample estimators (`r`, `smd`).
When in doubt, the function’s `@param` line in the help file is explicit
about which is wanted.

## Two Consistent Return Shapes

DMAR functions return one of two tidy data-frame shapes:

1.  **Long form: `data.frame(term, value)`**: used by all variance
    utilities, AIPE planners, effect size point estimates, and CI
    helpers. The `term` column names the quantity; `value` is the
    numeric value.
2.  **Wide form: one row, many columns**: used by
    [`omega_squared()`](https://yelleknek.github.io/DMAR/reference/omega_squared.md),
    [`eta_squared_partial()`](https://yelleknek.github.io/DMAR/reference/eta_squared_partial.md),
    and similar effect size functions that return one row per *effect*
    in a factorial design (with columns `effect`, `omega_squared`,
    `F_value`, `df_effect`, `df_error`, `N`).

Pipeline code can use `dplyr::filter(term == "...")` for the long form
or column-access for the wide form. The two shapes compose: the wide
form can always be reshaped to long with
[`tidyr::pivot_longer()`](https://tidyr.tidyverse.org/reference/pivot_longer.html).

## Recommended Workflows

### A. Reporting a Standardized Mean Difference

``` r

set.seed(113)
g1 <- rnorm(50, mean = 0, sd = 1)
g2 <- rnorm(50, mean = 0.5, sd = 1)

# 1. Point estimate (biased d and Hedges' g):
smd(group_1 = g1, group_2 = g2)                                 # d
smd(group_1 = g1, group_2 = g2, unbiased = TRUE)                 # g

# 2. CI on the standardized mean difference (Cohen's d):
ci_smd(ncp = t.test(g1, g2, var.equal = TRUE)$statistic,
       n_1 = 50, n_2 = 50)

# 3. Alternative effect size scales:
cles(smd = 0.5, n_1 = 50, n_2 = 50)            # probability of superiority
proportion_of_superiority(smd = 0.5, n_1 = 50, n_2 = 50)  # proportion above ctrl mean
nnt_from_smd(smd = 0.5, n_1 = 50, n_2 = 50)    # number-needed-to-treat

# 4. Variance (for meta-analytic weighting):
var_smd(delta = 0.5, n_1 = 50, n_2 = 50)
```

### B. Planning a Partial Correlation Study

``` r

# Anticipated partial r = 0.30, 3 controls, target 95% CI width = 0.20.

# 1. AIPE planner (raw-scale Olkin-Finn variance):
ss_aipe_partial_r(rho = 0.30, J = 3, width = 0.20)

# 2. Same plan on the Fisher-z scale (Bonett 2008):
ss_aipe_partial_r(rho = 0.30, J = 3, width = 0.20, fisher_z = TRUE)

# 3. Variance utility: what SE would result at the recommended n?
var_partial_r(r = 0.30, n = 80, J = 3)

# 4. Design-stage bias estimate: what r should we expect to observe?
expected_partial_r(rho = 0.30, n = 80, J = 3)
```

### C. Planning an Omega Squared CI

``` r

# 3-group one-way ANOVA, anticipated omega² = 0.10, target full width 0.10.
ss_aipe_omega_squared(population_omega_squared = 0.10,
                      df_effect = 2,
                      width = 0.10)

# Confirm: at that n, what does ci_omega_squared() actually return?
N <- ss_aipe_omega_squared(0.10, df_effect = 2, width = 0.10
                          )$value[1]
F_val <- 1 + 0.10 * N / (2 * 0.90)
ci_omega_squared(F_value = F_val, df_effect = 2,
                  df_error = N - 3, N = N)
```

### D. Reporting an Indirect (Mediated) Effect

``` r

# After fitting the X -> M and Y ~ X + M regressions:
a_hat   <- 0.40; b_hat <- 0.40
var_a   <- 0.02; var_b <- 0.02

# All four variance formulas in one call:
var_indirect_effect(a = a_hat, b = b_hat,
                    var_a = var_a, var_b = var_b)

# Plan an n for the next study:
ss_aipe_indirect_effect(a = 0.40, b = 0.40, width = 0.20,
                         method = "sobel")
```

### E. Two-Method Agreement Study

``` r

set.seed(113)
method_a <- rnorm(40, mean = 100, sd = 15)
method_b <- method_a + rnorm(40, mean = 2, sd = 5)

# 1. Lin's concordance correlation with CI:
lin_ccc(method_a, method_b)

# 2. Bland-Altman 95% limits of agreement with CIs on the limits:
loa(method_a, method_b)

# 3. Cohen's U3 / CLES interpretation if a method comparison can be
# framed as a "difference":
d_estimate <- mean(method_b - method_a) / sd(method_b - method_a)
cles(smd = d_estimate)
```

## Summary

DMAR’s architecture rests on four building blocks for any quantity
$`\theta`$: a point estimator with variance, an exact-noncentral CI, an
expected-value function for design-stage planning, and an AIPE sample-
size planner. The functions in each block compose horizontally (the same
$`\theta`$ shows up in all four) and vertically (within a block, related
transforms like $`d \to \mathrm{CLES} \to U_3 \to \mathrm{NNT}`$ are
explicit). Knowing which block a given task lives in is half the battle;
the function name and help file complete the picture.
