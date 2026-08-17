# Sequential Sample Size for a Fixed-Width Contrast Interval

Implements the purely sequential fixed-width confidence interval
procedure for a linear contrast \\\psi = \sum_j c_j \mu_j\\, following
Chattopadhyay, Bandyopadhyay, Kelley, and Padalunkal (2025). The goal is
a 100(1 - 2\\\alpha\\)% confidence interval \\\hat\psi \pm h\\ whose
half-width \\h\\ is fixed in advance, which is what makes a
noninferiority or equivalence verdict reachable by design: with bounds
\\\pm\delta\\, equivalence can never be declared unless \\h \< \delta\\,
and targeting \\h = \delta/2\\ gives a truly equivalent contrast about a
90% chance of being declared equivalent at \\\alpha = .05\\. Because no
fixed sample size can guarantee a bounded-width interval when the error
variance is unknown (Dantzig, 1940), the procedure is sequential: begin
with a pilot, then keep sampling, re-estimating the variance, until the
stopping criterion is met. Used with `pilot = TRUE` the function plans
the pilot and the cost-optimal allocation; with `pilot = FALSE` it
evaluates the stopping criterion at the data in hand, in the style of
[`mr_smd`](https://yelleknek.github.io/DMAR/reference/mr_smd.md).

## Usage

``` r
ss_seq_c(
  c_weights,
  half_width,
  s = NULL,
  n = NULL,
  cost = NULL,
  alpha_level = 0.05,
  quantile = c("t", "normal"),
  pilot = FALSE,
  m0 = 10
)
```

## Arguments

- c_weights:

  The contrast weights. The weights must sum to zero with the positive
  weights summing to 1 and the negative weights to -1, so that
  `half_width` is on the raw scale of the response.

- half_width:

  The target half-width \\h\\ of the 100(1 - 2\\\alpha\\)% confidence
  interval, in raw units of the response. For a noninferiority or
  equivalence decision at bounds \\\pm\delta\\, the recommended target
  is \\\delta/2\\.

- s:

  The current estimate(s) of the error standard deviation: either a
  single pooled value or one value per group (aligned with `c_weights`).
  Required when `pilot = FALSE`; optional planning values when
  `pilot = TRUE` (used only to shape the allocation).

- n:

  The current per-group sample sizes, aligned with `c_weights`. Required
  when `pilot = FALSE`.

- cost:

  Optional per-observation sampling costs, one per group, aligned with
  `c_weights`. When supplied, the reported allocation is the
  cost-optimal \\n_j \propto \|c_j\|\\ \sigma_j /
  \sqrt{\mathrm{cost}\_j}\\ (Chattopadhyay et al., 2025); with equal
  costs and standard deviations this reduces to allocation proportional
  to \\\|c_j\|\\.

- alpha_level:

  One-sided rate per bound; the interval is at confidence level 1 -
  2\\\alpha\\. Default `0.05` (a 90% interval).

- quantile:

  `"t"` (default) uses the *t* quantile on the current error degrees of
  freedom in the stopping criterion, a finite-sample refinement;
  `"normal"` uses the normal quantile of the classical statement of the
  rule (Chow & Robbins, 1965). The normal rule stops slightly early at
  wide targets, the finite-sample undershoot anticipated by Woodroofe
  (1977); the *t* rule corrects it at a small cost in sample size.

- pilot:

  `TRUE` to plan the pilot stage; `FALSE` (default) to evaluate the
  stopping criterion at the current data.

- m0:

  The minimum pilot sample size per group. Default `10`. A pilot that is
  too small makes the variance estimate that drives the early steps
  unstable.

## Value

With `pilot = TRUE`, a `data.frame` with row `pilot_n_per_group` (the
pilot size for each group with a nonzero weight) followed by one
`allocation_j` row per group giving the recommended sampling proportions
for the accrual stage. With `pilot = FALSE`, rows `stop` (1 = the
criterion is met, stop sampling; 0 = continue), `half_width_current`
(the half-width the interval would have now), `half_width_target`,
`N_current`, `N_projected` (the approximate total at which the criterion
would be met under the current allocation, from the normal
approximation), and the `allocation_j` rows for the next round of
sampling.

## Details

**The stopping criterion.** Sampling stops at the first \\N\\ for which
\$\$q^2 \sum_j c_j^2 s_j^2 / n_j \\\le\\ h^2,\$\$ where \\q\\ is the *t*
or normal quantile at \\1-\alpha\\. With a pooled \\s\\ and equal
allocation this is the Chow and Robbins (1965) rule specialized to a
contrast; the procedure is asymptotically consistent (coverage
approaches 1 - 2\\\alpha\\) and first-order efficient (the mean stopping
size approaches the oracle \\n^\*\\ an investigator with known variance
would use). See
[`ss_seq_c_sensitivity`](https://yelleknek.github.io/DMAR/reference/ss_seq_c_sensitivity.md)
for a Monte Carlo evaluation of both properties.

**Degrees of freedom.** With a pooled `s` the criterion uses \\\nu = N -
J\\. With per-group `s` it uses the Satterthwaite approximation, which
is the appropriate error law when the variances are not assumed
homogeneous.

**Batches.** Observations may be added in batches rather than one at a
time; the asymptotic properties survive batching. Re-call the function
after each batch.

## References

Chattopadhyay, B., Bandyopadhyay, T., Kelley, K., & Padalunkal, J. J.
(2025). A sequential approach for noninferiority or equivalence of a
linear contrast under cost constraints. *Psychological Methods, 30*(2),
425–439. [doi:10.1037/met0000570](https://doi.org/10.1037/met0000570)

Chow, Y. S., & Robbins, H. (1965). On the asymptotic theory of
fixed-width sequential confidence intervals for the mean. *The Annals of
Mathematical Statistics, 36*(2), 457–462.

Dantzig, G. B. (1940). On the non-existence of tests of "Student's"
hypothesis having power functions independent of \\\sigma\\. *The Annals
of Mathematical Statistics, 11*(2), 186–192.

Mukhopadhyay, N., & de Silva, B. M. (2009). *Sequential methods and
their applications*. CRC Press.

Woodroofe, M. (1977). Second order approximations for sequential point
and interval estimation. *The Annals of Statistics, 5*(5), 984–995.

## See also

[`ss_seq_c_sensitivity`](https://yelleknek.github.io/DMAR/reference/ss_seq_c_sensitivity.md),
[`ss_aipe_c`](https://yelleknek.github.io/DMAR/reference/ss_aipe_c.md),
[`ss_power_equivalence_c`](https://yelleknek.github.io/DMAR/reference/ss_power_equivalence_c.md),
[`equivalence_c`](https://yelleknek.github.io/DMAR/reference/equivalence_c.md),
[`mr_smd`](https://yelleknek.github.io/DMAR/reference/mr_smd.md)

Other sequential estimation:
[`ss_seq_c_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_seq_c_sensitivity.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# 1. Plan the pilot for a two-group contrast, target half-width 2.5
#    (bounds of 5 with the h = delta/2 rule):
ss_seq_c(c_weights = c(1, -1), half_width = 2.5, pilot = TRUE)
#>  term              value
#>  pilot_n_per_group 10   
#>  allocation_1      0.5  
#>  allocation_2      0.5  
#> 
#> Confidence level: 90%

# 2. Evaluate the stopping criterion mid-study: pooled s = 15.4 at
#    n = 60 per group. Too imprecise to stop; the projection says
#    roughly how much further to go.
ss_seq_c(c_weights = c(1, -1), half_width = 2.5,
         s = 15.4, n = c(60, 60))
#>  term               value
#>  stop               0    
#>  half_width_current 4.66 
#>  half_width_target  2.5  
#>  N_current          120  
#>  N_projected        411  
#>  allocation_1       0.5  
#>  allocation_2       0.5  
#> 
#> Confidence level: 90%

# 3. Costs differ: sampling the second group costs four times as
#    much per observation, so its share of new observations drops.
ss_seq_c(c_weights = c(1, -1), half_width = 2.5,
         s = c(15.4, 15.4), n = c(60, 60), cost = c(1, 4))
#>  term               value
#>  stop               0    
#>  half_width_current 4.66 
#>  half_width_target  2.5  
#>  N_current          120  
#>  N_projected        411  
#>  allocation_1       0.667
#>  allocation_2       0.333
#> 
#> Confidence level: 90%
```
