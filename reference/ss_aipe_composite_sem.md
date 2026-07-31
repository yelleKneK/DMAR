# Sample Size for Accurate Estimation of a Set of SEM Parameters

Determine the necessary sample size for a structural equation model
study so that the confidence interval for every parameter of interest is
sufficiently narrow in the same study, or, given a sample size, return
how narrow the set of intervals can be expected to be. This is the
accuracy in parameter estimation (AIPE) counterpart of
[`ss_power_composite_sem`](https://yelleknek.github.io/DMAR/reference/ss_power_composite_sem.md):
where that function plans for every parameter to be statistically
significant jointly, this one plans for every parameter to be estimated
with a confidence interval no wider than desired, the goal when the
research questions concern the magnitudes of the effects rather than
their existence. The parameters of interest are any labeled parameters
of a lavaan analysis model, structural paths, loadings, covariances, or
quantities defined with `:=` such as an indirect effect, and any subset
of them can make up the set.

## Usage

``` r
ss_aipe_composite_sem(
  model,
  Sigma = NULL,
  pop_model = NULL,
  mu = NULL,
  parameters = NULL,
  desired_width,
  conf_level = 0.95,
  assurance = NULL,
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
  structure). The fixed values are values the researcher posits (from
  theory, prior studies, or pilot data), never sample estimates. Supply
  exactly one of `Sigma` or `pop_model`.

- mu:

  Optional population means of the observed variables, used with
  `Sigma`: a named numeric vector with one entry per observed variable,
  or an unnamed vector in the row order of `Sigma`. The default `NULL`
  is zero means. Means matter only when the analysis model has a mean
  structure (an intercept term such as `s ~ 1`, as in a latent growth
  curve model); when `pop_model` is supplied its mean structure provides
  the means and `mu` must not also be given.

- parameters:

  Character vector of the parameter labels that make up the set. The
  default `NULL` uses every user-labeled parameter in `model`, in order
  of appearance, so labeling exactly the parameters of interest is the
  simplest way to state the set.

- desired_width:

  The desired full confidence interval width for each parameter of
  interest: a single value applied to every parameter, or a named
  numeric vector with one entry per parameter label. An unnamed vector
  of several widths is not accepted, so a width can never silently
  attach to the wrong parameter.

- conf_level:

  Confidence level of each interval (default 0.95).

- assurance:

  The desired probability that a single study yields confidence
  intervals no wider than desired for every parameter of interest
  simultaneously (a value in \[0.5, 1)), or `NULL` (the default) to plan
  against the expected widths instead; see Details.

- N:

  Sample size; if supplied, the realized interval widths at that *N* are
  summarized rather than a sample size planned.

- G:

  Number of converged Monte Carlo replications per evaluated sample size
  (default 1000). The simulation error of each estimated proportion is
  about \\\sqrt{p(1 - p)/G}\\; raise `G` for a sharper answer.

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

A `data.frame` (a `dmar_tbl`) with `term` and `value` columns: the
`necessary_N` (or supplied `specified_N`), the `composite_assurance`
(the proportion of replications in which every interval was
simultaneously within its desired width, reported under both criteria),
then for each parameter its `mean_width_<label>`, its marginal
`width_within_desired_<label>` proportion, its `desired_width_<label>`,
and its purported `population_<label>` value under the analysis model,
followed by `conf_level`, the requested `replications`, the
`converged_replications` the summary is based on, and, when supplied,
the `assurance`.

## Details

AIPE planning for a single targeted SEM parameter is available in closed
form (Lai & Kelley, 2011;
[`ss_aipe_sem_path`](https://yelleknek.github.io/DMAR/reference/ss_aipe_sem_path.md)),
but most studies estimate several effects and report all of them; a
design is only as informative as its widest interval of interest. This
function plans for the set by a priori Monte Carlo simulation (Muthén &
Muthén, 2002; Maxwell, Kelley, & Rausch, 2008): for a candidate *N*, `G`
data sets are drawn from the multivariate normal population with
covariance matrix `Sigma`, the analysis model is fit to each, and each
parameter's Wald confidence interval width, twice \\z\_{1 - \alpha/2}\\
times its standard error, is recorded. Because the estimates share one
fitted model, the widths are dependent; the simulation reflects that
dependence exactly, at the stated *N*, with no asymptotic shortcut.

Two planning criteria are available. With `assurance = NULL` the
necessary sample size is the smallest *N* at which the mean simulated
width of every parameter's interval is at or below its desired width,
the expected-width criterion of the AIPE framework applied to each
member of the set. Widths vary from sample to sample around their means,
so each interval separately lands at or below its desired width in
roughly half of the realizations. That is a statement about one interval
at a time, not about the set: the probability that *every* interval is
narrow enough at once falls well below one half as soon as more than one
parameter binds, and falls further the more parameters are targeted and
the more weakly their widths move together. Planning the whole set to a
stated probability is exactly what `assurance` is for. Supplying
`assurance` plans against the joint event instead: the smallest *N* at
which the proportion of replications where every interval is
simultaneously within its desired width reaches the assurance. The joint
event is contained in each marginal event, so its probability is at most
the smallest per-parameter proportion, and the
`width_within_desired_<label>` rows show which parameter binds the
design.

When `N` is `NULL` the search starts at the largest of the per-parameter
closed-form sample sizes (the no-assurance approximation
[`ss_aipe_sem_path`](https://yelleknek.github.io/DMAR/reference/ss_aipe_sem_path.md)
uses, computed from the asymptotic variances before any simulation),
brackets the crossing geometrically, and bisects to adjacent integers,
each candidate evaluated with its own `G` replications. A planning call
therefore fits the analysis model several thousand times; the examples
use a small `G` to run quickly, and a real plan is worth a larger one.

Each reported proportion carries a simulation standard error of about
\\\sqrt{p(1 - p)/G}\\, and the necessary sample size inherits that
uncertainty; raising `G` narrows it, and reporting the seed makes a plan
reproducible.

## Note

A replication whose fit does not converge, or converges without a usable
standard error for some parameter of interest, is discarded and fresh
data are drawn, up to `20 * G` attempts per evaluated sample size; the
reported summaries condition on convergence. When fewer than `G`
replications converge within the cap, a single warning is issued and the
summary is based on the converged replications (their count is the
`converged_replications` row).

Because the planner itself is a Monte Carlo study, it has no separate
`_sensitivity` sibling; to study misspecification of the population
values, rerun the planner with the alternative `Sigma` or `pop_model`
values under consideration and compare the plans.

## References

Lai, K., & Kelley, K. (2011). Accuracy in parameter estimation for
targeted effects in structural equation modeling: Sample size planning
for narrow confidence intervals. *Psychological Methods, 16*(2),
127–148. [doi:10.1037/a0021764](https://doi.org/10.1037/a0021764)

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge.

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

## See also

[`ss_power_composite_sem`](https://yelleknek.github.io/DMAR/reference/ss_power_composite_sem.md)
for the same set of parameters planned for joint statistical
significance;
[`cov_sem`](https://yelleknek.github.io/DMAR/reference/cov_sem.md) for
deriving `Sigma` from a fully fixed population model;
[`ss_aipe_sem_path`](https://yelleknek.github.io/DMAR/reference/ss_aipe_sem_path.md)
and
[`ss_aipe_sem_path_sensitivity`](https://yelleknek.github.io/DMAR/reference/ss_aipe_sem_path_sensitivity.md)
for a single targeted path.

Other AIPE sample size planning:
[`ss_aipe_c_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_c_sensitivity.md),
[`ss_aipe_cliff_delta()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_cliff_delta.md),
[`ss_aipe_cliff_delta_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_cliff_delta_sensitivity.md),
[`ss_aipe_icc()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_icc.md),
[`ss_aipe_icc_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_icc_sensitivity.md),
[`ss_aipe_indirect_effect()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_indirect_effect.md),
[`ss_aipe_indirect_effect_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_indirect_effect_sensitivity.md),
[`ss_aipe_mixed_effects_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_mixed_effects_sensitivity.md),
[`ss_aipe_omega_squared()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_omega_squared.md),
[`ss_aipe_omega_squared_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_omega_squared_sensitivity.md),
[`ss_aipe_partial_r()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_partial_r.md),
[`ss_aipe_partial_r_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_partial_r_sensitivity.md),
[`ss_aipe_pcm_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_pcm_sensitivity.md),
[`ss_aipe_reliability_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_reliability_sensitivity.md),
[`ss_aipe_semipartial_r()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_semipartial_r.md),
[`ss_aipe_semipartial_r_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_semipartial_r_sensitivity.md),
[`ss_aipe_tost_smd_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_tost_smd_sensitivity.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# \donttest{
# A mediation model whose research questions concern the magnitudes of
# both individual paths and the indirect effect. The population model
# fixes every parameter to its purported population value.
pop_model <- "
  f1 =~ 1*y1 + 0.8*y2 + 0.8*y3
  f2 =~ 1*y4 + 0.8*y5 + 0.8*y6
  f3 =~ 1*y7 + 0.8*y8 + 0.8*y9
  f2 ~ 0.4*f1
  f3 ~ 0.5*f2 + 0.2*f1
  f1 ~~ 1*f1
  f2 ~~ 0.84*f2
  f3 ~~ 0.7*f3
  y1 ~~ 0.5*y1; y2 ~~ 0.5*y2; y3 ~~ 0.5*y3
  y4 ~~ 0.5*y4; y5 ~~ 0.5*y5; y6 ~~ 0.5*y6
  y7 ~~ 0.5*y7; y8 ~~ 0.5*y8; y9 ~~ 0.5*y9
"

# The analysis model labels the two paths and defines the indirect
# effect; all three make up the set of interest.
analysis_model <- "
  f1 =~ y1 + y2 + y3
  f2 =~ y4 + y5 + y6
  f3 =~ y7 + y8 + y9
  f2 ~ a*f1
  f3 ~ b*f2 + cp*f1
  ab := a*b
"

# Realized interval widths at N = 200. G is small so the example runs
# quickly; a real plan is worth G = 1000 or more (G = 50 here).
set.seed(113)
ss_aipe_composite_sem(model = analysis_model, pop_model = pop_model,
                      parameters = c("a", "b", "ab"),
                      desired_width = 0.30, N = 200, G = 50)
#>  term                    value
#>  specified_N             200  
#>  composite_assurance     0    
#>  mean_width_a            0.343
#>  mean_width_b            0.369
#>  mean_width_ab           0.215
#>  width_within_desired_a  0.12 
#>  width_within_desired_b  0.02 
#>  width_within_desired_ab 0.98 
#>  desired_width_a         0.3  
#>  desired_width_b         0.3  
#>  desired_width_ab        0.3  
#>  population_a            0.4  
#>  population_b            0.5  
#>  population_ab           0.2  
#>  conf_level              0.95 
#>  replications            50   
#>  converged_replications  50   
#> 
#> Confidence level: 95%

# The necessary N for all three intervals to be simultaneously within
# their desired widths in 80 percent of studies, with the indirect
# effect held to a narrower interval than the paths.
set.seed(113)
ss_aipe_composite_sem(model = analysis_model, pop_model = pop_model,
                      parameters = c("a", "b", "ab"),
                      desired_width = c(a = 0.35, b = 0.35, ab = 0.25),
                      assurance = 0.80, G = 50)
#>  term                    value
#>  necessary_N             282  
#>  composite_assurance     0.82 
#>  mean_width_a            0.289
#>  mean_width_b            0.316
#>  mean_width_ab           0.185
#>  width_within_desired_a  0.98 
#>  width_within_desired_b  0.84 
#>  width_within_desired_ab 1    
#>  desired_width_a         0.35 
#>  desired_width_b         0.35 
#>  desired_width_ab        0.25 
#>  population_a            0.4  
#>  population_b            0.5  
#>  population_ab           0.2  
#>  conf_level              0.95 
#>  replications            50   
#>  converged_replications  50   
#>  assurance               0.8  
#> 
#> Confidence level: 95%
# }
```
