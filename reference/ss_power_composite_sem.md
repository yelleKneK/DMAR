# Sample Size or Composite Power for a Set of SEM Parameters

Determine the necessary sample size for a structural equation model
study so that every parameter of interest is statistically significant
in the same study with a desired probability, or, given a sample size,
return that probability. Composite power is the probability that all of
the named parameters are significant jointly, the quantity a design must
be planned against when its conclusion requires more than one result to
hold at once: a study can have adequate power for each hypothesis on its
own and still be underpowered for the conclusion that rests on all of
them together (Maxwell, 2004). The parameters of interest are any
labeled parameters of a lavaan analysis model, structural paths,
loadings, covariances, or quantities defined with `:=` such as an
indirect effect, and any subset of them can make up the composite.

## Usage

``` r
ss_power_composite_sem(
  model,
  Sigma = NULL,
  pop_model = NULL,
  mu = NULL,
  parameters = NULL,
  desired_power = 0.85,
  alpha_level = 0.05,
  N = NULL,
  G = 1000,
  seed = NULL,
  ...
)
```

## Arguments

- model:

  A single character string giving the free analysis model in lavaan
  model syntax (see
  [`model.syntax`](https://rdrr.io/pkg/lavaan/man/model.syntax.html)),
  the model that would be fit to the data. Each parameter of interest
  must carry a parameter label so it can be referred to by name, for
  example `"f2 ~ b*f1"` labels the structural path `b`, and
  `"ab := a*b"` defines an indirect effect from the labeled paths `a`
  and `b`.

- Sigma:

  Population covariance matrix of the observed variables, with row and
  column names matching the observed variables in `model`. It is
  typically obtained from a fully fixed population model via
  [`cov_sem`](https://yelleknek.github.io/DMAR/reference/cov_sem.md).
  Supply exactly one of `Sigma` or `pop_model`.

- pop_model:

  A single character string giving the population model in lavaan model
  syntax with every parameter fixed to its population value, from which
  [`cov_sem`](https://yelleknek.github.io/DMAR/reference/cov_sem.md)
  derives `Sigma` (and the population means, when the model has a mean
  structure). This is where the purported population values of the
  parameters of interest are chosen; they are values the researcher
  posits (from theory, prior studies, or pilot data), never sample
  estimates. Supply exactly one of `Sigma` or `pop_model`.

- mu:

  Optional population means of the observed variables, used with
  `Sigma`: a named numeric vector with one entry per observed variable,
  or an unnamed vector in the row order of `Sigma`. The default `NULL`
  is zero means. Means matter only when the analysis model has a mean
  structure (an intercept term such as `s ~ 1`, as in a latent growth
  curve model); when `pop_model` is supplied its mean structure provides
  the means and `mu` must not also be given.

- parameters:

  Character vector of the parameter labels that make up the composite.
  The default `NULL` uses every user-labeled parameter in `model`, in
  order of appearance, so labeling exactly the parameters of interest is
  the simplest way to state the set.

- desired_power:

  Desired composite statistical power (default 0.85). Used only when `N`
  is `NULL`.

- alpha_level:

  Type I error rate for each individual two-sided Wald *z* test (default
  0.05), the per-test rate, not a rate for the composite event.

- N:

  Sample size; if supplied, the realized composite power at that *N* is
  returned rather than a sample size planned.

- G:

  Number of converged Monte Carlo replications per evaluated sample size
  (default 1000). The simulation error of each estimated power is about
  \\\sqrt{p(1 - p)/G}\\; raise `G` for a sharper answer.

- seed:

  Optional integer seed for reproducibility. The default `NULL` uses the
  current state of the random number generator; a supplied seed is set
  internally and the prior state restored on exit.

- ...:

  Additional arguments passed to
  [`sem`](https://rdrr.io/pkg/lavaan/man/sem.html), both when the
  population values are resolved and for every Monte Carlo fit (for
  example `std.lv = TRUE` or `missing = "listwise"`). A robust estimator
  such as `"MLM"` cannot be used here: the same arguments reach the
  setup fit, which is always from summary statistics, and lavaan refuses
  a robust estimator there. The Monte Carlo data are multivariate normal
  by construction, so a robust estimator would buy nothing.

## Value

A `data.frame` with `term` and `value` columns: the `necessary_N` (or
supplied `specified_N`), the `composite_power` and its simulation
standard error `composite_power_mc_se`, then for each parameter its
marginal `power_<label>` and purported `population_<label>` value under
the analysis model, followed by `alpha_level`, the requested
`replications`, the `converged_replications` the summary is based on,
and, when a size was planned, the `desired_power`. The result carries
the `dmar_ss_power` class, so
[`tidy`](https://generics.r-lib.org/reference/tidy.html) and
[`glance`](https://generics.r-lib.org/reference/glance.html) summarize
the sample size and the composite power in broom convention.

## Details

Analytic sample size planning methods in SEM exist for a single targeted
parameter (Satorra & Saris, 1985; Lai & Kelley, 2011) or for overall
model fit (MacCallum, Browne, & Sugawara, 1996;
[`ss_power_sem`](https://yelleknek.github.io/DMAR/reference/ss_power_sem.md)),
but most studies state several hypotheses and support their conclusion
only when all of them hold. This function plans for that case by a
priori Monte Carlo simulation (Muthén & Muthén, 2002; Maxwell, Kelley, &
Rausch, 2008): for a candidate *N*, `G` data sets are drawn from the
multivariate normal population with covariance matrix `Sigma`, the
analysis model is fit to each, and each parameter of interest is tested
with its two-sided Wald *z* test at `alpha_level`. The proportion of
replications in which every parameter is significant estimates the
composite power, and the per-parameter proportions estimate the marginal
powers. Because the estimates share one fitted model, the tests are
dependent; the simulation reflects that dependence exactly, at the
stated *N*, with no asymptotic shortcut.

The composite event is contained in each marginal event, so composite
power is at most the smallest marginal power: the weakest parameter
governs the design, and the marginal `power_<label>` rows show which
parameter that is.

When `N` is `NULL` the necessary sample size is searched for. The search
starts where the product of the marginal Wald powers (an independence
approximation computed from the asymptotic variances, spent before any
simulation) reaches `desired_power`, brackets the crossing
geometrically, and bisects to adjacent integers, each candidate
evaluated with its own `G` replications. A planning call therefore fits
the analysis model several thousand times; the examples use a small `G`
to run quickly, and a real plan is worth a larger one.

## Note

A replication whose fit does not converge, or converges without a usable
standard error for some parameter of interest, is discarded and fresh
data are drawn, up to `20 * G` attempts per evaluated sample size; the
reported powers condition on convergence. When fewer than `G`
replications converge within the cap, a single warning is issued and the
summary is based on the converged replications (their count is the
`converged_replications` row). Frequent nonconvergence at small *N* is
itself design information: a sample size at which the model rarely
converges is too small in a sense that precedes power.

## Monte Carlo Precision

Each reported power is a proportion of `G` replications, with simulation
standard error about \\\sqrt{p(1 - p)/G}\\; the `composite_power_mc_se`
row reports it for the composite. The necessary sample size inherits
that uncertainty: near the target the power curve is flat enough that
neighboring *N* are separated by less than the simulation error, so
repeated calls with different seeds return slightly different sizes.
Raising `G` narrows the spread; reporting the seed makes a plan
reproducible. A proportion of `G` replications takes only the values
\\0, 1/G, \ldots, 1\\, so a `desired_power` above \\1 - 1/G\\ is refused
with a message saying how large `G` must be for that target; the same
resolution guard applies to
[`ss_aipe_composite_sem`](https://yelleknek.github.io/DMAR/reference/ss_aipe_composite_sem.md)'s
`assurance`.

## References

Lai, K., & Kelley, K. (2011). Accuracy in parameter estimation for
targeted effects in structural equation modeling: Sample size planning
for narrow confidence intervals. *Psychological Methods, 16*(2),
127–148. [doi:10.1037/a0021764](https://doi.org/10.1037/a0021764)

MacCallum, R. C., Browne, M. W., & Sugawara, H. M. (1996). Power
analysis and determination of sample size for covariance structure
modeling. *Psychological Methods, 1*(2), 130–149.
[doi:10.1037/1082-989X.1.2.130](https://doi.org/10.1037/1082-989X.1.2.130)

Maxwell, S. E. (2004). The persistence of underpowered studies in
psychological research: Causes, consequences, and remedies.
*Psychological Methods, 9*(2), 147–163.
[doi:10.1037/1082-989X.9.2.147](https://doi.org/10.1037/1082-989X.9.2.147)

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge. (See Chapter 3 on statistical power.)

Maxwell, S. E., Kelley, K., & Rausch, J. R. (2008). Sample size planning
for statistical power and accuracy in parameter estimation. *Annual
Review of Psychology, 59*, 537–563.
[doi:10.1146/annurev.psych.59.103006.093735](https://doi.org/10.1146/annurev.psych.59.103006.093735)

Muthén, L. K., & Muthén, B. O. (2002). How to use a Monte Carlo study to
decide on sample size and determine power. *Structural Equation
Modeling, 9*(4), 599–620.
[doi:10.1207/S15328007SEM0904_8](https://doi.org/10.1207/S15328007SEM0904_8)

Rosseel, Y. (2012). lavaan: An R package for structural equation
modeling. *Journal of Statistical Software, 48*(2), 1–36.
[doi:10.18637/jss.v048.i02](https://doi.org/10.18637/jss.v048.i02)

Satorra, A., & Saris, W. E. (1985). Power of the likelihood ratio test
in covariance structure analysis. *Psychometrika, 50*(1), 83–90.

## See also

[`ss_aipe_composite_sem`](https://yelleknek.github.io/DMAR/reference/ss_aipe_composite_sem.md)
for the same set of parameters planned for accuracy in parameter
estimation (AIPE) instead of significance;
[`cov_sem`](https://yelleknek.github.io/DMAR/reference/cov_sem.md) for
deriving `Sigma` from a fully fixed population model;
[`ss_power_sem`](https://yelleknek.github.io/DMAR/reference/ss_power_sem.md)
for overall model fit;
[`ss_aipe_sem_path`](https://yelleknek.github.io/DMAR/reference/ss_aipe_sem_path.md)
for a single targeted path;
[`ss_power_composite_anova`](https://yelleknek.github.io/DMAR/reference/ss_power_composite_anova.md)
and its siblings for composite power in ANOVA and ANCOVA designs, where
the composite is evaluated by quadrature rather than simulation.

Other sample size for power:
[`power_fisher_exact()`](https://yelleknek.github.io/DMAR/reference/power_fisher_exact.md),
[`ss_aipe_mixed_effects()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_mixed_effects.md),
[`ss_aipe_tost_smd()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_tost_smd.md),
[`ss_power_R2()`](https://yelleknek.github.io/DMAR/reference/ss_power_R2.md),
[`ss_power_R2_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_power_R2_sensitivity.md),
[`ss_power_c()`](https://yelleknek.github.io/DMAR/reference/ss_power_c.md),
[`ss_power_c_ancova()`](https://yelleknek.github.io/DMAR/reference/ss_power_c_ancova.md),
[`ss_power_composite_ancova()`](https://yelleknek.github.io/DMAR/reference/ss_power_composite_ancova.md),
[`ss_power_composite_ancova_2group()`](https://yelleknek.github.io/DMAR/reference/ss_power_composite_ancova_2group.md),
[`ss_power_composite_anova()`](https://yelleknek.github.io/DMAR/reference/ss_power_composite_anova.md),
[`ss_power_composite_factorial_ancova()`](https://yelleknek.github.io/DMAR/reference/ss_power_composite_factorial_ancova.md),
[`ss_power_composite_factorial_ancova_het()`](https://yelleknek.github.io/DMAR/reference/ss_power_composite_factorial_ancova_het.md),
[`ss_power_composite_factorial_anova()`](https://yelleknek.github.io/DMAR/reference/ss_power_composite_factorial_anova.md),
[`ss_power_contrast()`](https://yelleknek.github.io/DMAR/reference/ss_power_contrast.md),
[`ss_power_equivalence_c()`](https://yelleknek.github.io/DMAR/reference/ss_power_equivalence_c.md),
[`ss_power_factorial_ancova()`](https://yelleknek.github.io/DMAR/reference/ss_power_factorial_ancova.md),
[`ss_power_factorial_anova()`](https://yelleknek.github.io/DMAR/reference/ss_power_factorial_anova.md),
[`ss_power_indirect_effect()`](https://yelleknek.github.io/DMAR/reference/ss_power_indirect_effect.md),
[`ss_power_mixed_effects()`](https://yelleknek.github.io/DMAR/reference/ss_power_mixed_effects.md),
[`ss_power_one_way_anova()`](https://yelleknek.github.io/DMAR/reference/ss_power_one_way_anova.md),
[`ss_power_pcm()`](https://yelleknek.github.io/DMAR/reference/ss_power_pcm.md),
[`ss_power_r()`](https://yelleknek.github.io/DMAR/reference/ss_power_r.md),
[`ss_power_rc()`](https://yelleknek.github.io/DMAR/reference/ss_power_rc.md),
[`ss_power_reg_coef()`](https://yelleknek.github.io/DMAR/reference/ss_power_reg_coef.md),
[`ss_power_reg_coef_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_power_reg_coef_sensitivity.md),
[`ss_power_rm_anova()`](https://yelleknek.github.io/DMAR/reference/ss_power_rm_anova.md),
[`ss_power_sc()`](https://yelleknek.github.io/DMAR/reference/ss_power_sc.md),
[`ss_power_sem()`](https://yelleknek.github.io/DMAR/reference/ss_power_sem.md),
[`ss_power_smd()`](https://yelleknek.github.io/DMAR/reference/ss_power_smd.md),
[`ss_power_split_plot_anova()`](https://yelleknek.github.io/DMAR/reference/ss_power_split_plot_anova.md)

Other composite power:
[`ss_power_composite_ancova()`](https://yelleknek.github.io/DMAR/reference/ss_power_composite_ancova.md),
[`ss_power_composite_ancova_2group()`](https://yelleknek.github.io/DMAR/reference/ss_power_composite_ancova_2group.md),
[`ss_power_composite_anova()`](https://yelleknek.github.io/DMAR/reference/ss_power_composite_anova.md),
[`ss_power_composite_factorial_ancova()`](https://yelleknek.github.io/DMAR/reference/ss_power_composite_factorial_ancova.md),
[`ss_power_composite_factorial_ancova_het()`](https://yelleknek.github.io/DMAR/reference/ss_power_composite_factorial_ancova_het.md),
[`ss_power_composite_factorial_anova()`](https://yelleknek.github.io/DMAR/reference/ss_power_composite_factorial_anova.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# \donttest{
# A three-factor model whose conclusion needs both structural paths: f1
# predicting f2, and f2 predicting f3 over and above f1. The population
# model fixes every parameter to its purported population value.
pop_model <- "
  f1 =~ 1*y1 + 0.8*y2 + 0.8*y3
  f2 =~ 1*y4 + 0.8*y5 + 0.8*y6
  f3 =~ 1*y7 + 0.8*y8 + 0.8*y9
  f2 ~ 0.4*f1
  f3 ~ 0.3*f2 + 0.25*f1
  f1 ~~ 1*f1
  f2 ~~ 0.84*f2
  f3 ~~ 0.8*f3
  y1 ~~ 0.5*y1; y2 ~~ 0.5*y2; y3 ~~ 0.5*y3
  y4 ~~ 0.5*y4; y5 ~~ 0.5*y5; y6 ~~ 0.5*y6
  y7 ~~ 0.5*y7; y8 ~~ 0.5*y8; y9 ~~ 0.5*y9
"

# The analysis model is free; the labels name the parameters of interest.
analysis_model <- "
  f1 =~ y1 + y2 + y3
  f2 =~ y4 + y5 + y6
  f3 =~ y7 + y8 + y9
  f2 ~ a*f1
  f3 ~ b*f2 + c*f1
"

# Realized composite power at N = 200. G is small so the example runs
# quickly; a real plan is worth G = 1000 or more.
set.seed(113)
ss_power_composite_sem(model = analysis_model, pop_model = pop_model,
                       N = 200, G = 100)
#>  term                   value
#>  specified_N            200  
#>  composite_power        0.67 
#>  composite_power_mc_se  0.047
#>  power_a                1    
#>  power_b                0.89 
#>  power_c                0.78 
#>  population_a           0.4  
#>  population_b           0.3  
#>  population_c           0.25 
#>  alpha_level            0.05 
#>  replications           100  
#>  converged_replications 100  

# The necessary N for composite power of 0.80 over the two structural
# paths a and b (c is left out of the composite).
set.seed(113)
ss_power_composite_sem(model = analysis_model, pop_model = pop_model,
                       parameters = c("a", "b"),
                       desired_power = 0.80, G = 100)
#>  term                   value 
#>  necessary_N            154   
#>  composite_power        0.88  
#>  composite_power_mc_se  0.0325
#>  power_a                1     
#>  power_b                0.88  
#>  population_a           0.4   
#>  population_b           0.3   
#>  alpha_level            0.05  
#>  replications           100   
#>  converged_replications 100   
#>  desired_power          0.8   
# }
```
