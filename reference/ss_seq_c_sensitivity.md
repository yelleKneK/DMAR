# Monte Carlo Sensitivity of the Sequential Fixed-Width Procedure

Simulates the purely sequential fixed-width confidence interval
procedure of
[`ss_seq_c`](https://yelleknek.github.io/DMAR/reference/ss_seq_c.md)
under a known data generating mechanism, reporting the distribution of
the stopping sample size and the empirical coverage of the fixed-width
interval. The two quantities to read are the ratio of the mean stopping
size to the oracle \\n^\*\\ (first-order efficiency: the ratio
approaches 1 as the target half-width shrinks) and the coverage
(asymptotic consistency: coverage approaches 1 - 2\\\alpha\\). The
normal quantile rule stops slightly early at wide targets, the
finite-sample undershoot anticipated by Woodroofe (1977); the *t*
quantile rule corrects it at a small cost in sample size.

## Usage

``` r
ss_seq_c_sensitivity(
  c_weights,
  half_width,
  true_sigma,
  true_means = NULL,
  alpha_level = 0.05,
  quantile = c("t", "normal"),
  m0 = 10,
  G = 1000,
  seed = NULL
)
```

## Arguments

- c_weights:

  The contrast weights. The weights must sum to zero with the positive
  weights summing to 1 and the negative weights to -1.

- half_width:

  The target half-width \\h\\ of the 100(1 - 2\\\alpha\\)% interval, in
  raw units of the response.

- true_sigma:

  The data generating error standard deviation: a single value applied
  to every group.

- true_means:

  Optional vector of data generating group means, aligned with
  `c_weights`. Default all zero, so the true contrast is 0.

- alpha_level:

  One-sided rate per bound; the interval is at confidence level 1 -
  2\\\alpha\\. Default `0.05`.

- quantile:

  `"t"` (default) or `"normal"`; see
  [`ss_seq_c`](https://yelleknek.github.io/DMAR/reference/ss_seq_c.md).

- m0:

  Pilot sample size per group. Default `10`.

- G:

  Number of Monte Carlo replications. Default `1000`.

- seed:

  Optional integer seed. Default `NULL` (the current RNG state is used).
  When supplied, the caller's RNG state is restored on exit.

## Value

A `data.frame` with rows `n_star` (the oracle total sample size an
investigator with known \\\sigma\\ would use), `mean_N`, `median_N`,
`sd_N` (the stopping total across replications), `ratio_mean_N_n_star`,
`coverage` (the proportion of replications whose \\\hat\psi_N \pm h\\
interval covered the true contrast), and `se_coverage` (its simulation
standard error).

## Details

**Simulation design.** Each replication samples the groups with nonzero
weights in balanced fashion (one observation per group per step) from
normal populations with common `true_sigma`, starting at `m0` per group,
and stops at the first step satisfying the
[`ss_seq_c`](https://yelleknek.github.io/DMAR/reference/ss_seq_c.md)
criterion with the pooled variance estimate. This matches the
equal-cost, equal-variance case of Chattopadhyay, Bandyopadhyay, Kelley,
and Padalunkal (2025); unequal costs change the optimal allocation but
not the logic.

**The oracle.** With known \\\sigma\\ and balanced allocation over the
\\J_0\\ groups with nonzero weights, the fixed-width requirement is
\\n^\* = z\_{1-\alpha}^2\\ \sigma^2 J_0 \sum_j c_j^2 / h^2\\ in total.
The sequential procedure spends about \\n^\*\\ without knowing
\\\sigma\\, which is its point.

## References

Chattopadhyay, B., Bandyopadhyay, T., Kelley, K., & Padalunkal, P. J.
(2025). A sequential approach for noninferiority or equivalence of a
linear contrast under cost constraints. *Psychological Methods, 30*(2),
425–439. [doi:10.1037/met0000570](https://doi.org/10.1037/met0000570)

Chow, Y. S., & Robbins, H. (1965). On the asymptotic theory of
fixed-width sequential confidence intervals for the mean. *The Annals of
Mathematical Statistics, 36*(2), 457–462.

Ghosh, M., Mukhopadhyay, N., & Sen, P. K. (1997). *Sequential
estimation*. Wiley.

Woodroofe, M. (1977). Second order approximations for sequential point
and interval estimation. *The Annals of Statistics, 5*(5), 984–995.

## See also

[`ss_seq_c`](https://yelleknek.github.io/DMAR/reference/ss_seq_c.md),
[`ss_aipe_c`](https://yelleknek.github.io/DMAR/reference/ss_aipe_c.md),
[`ss_power_equivalence_c`](https://yelleknek.github.io/DMAR/reference/ss_power_equivalence_c.md)

Other sequential estimation:
[`ss_seq_c()`](https://yelleknek.github.io/DMAR/reference/ss_seq_c.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# A two-group contrast, target half-width 2.5, error SD 15.67:
# the t-quantile rule stops near the oracle with near-nominal
# coverage. (G kept small here for speed; use G = 2000 or more in
# earnest.)
# \donttest{
ss_seq_c_sensitivity(c_weights = c(1, -1), half_width = 2.5,
                     true_sigma = 15.67, G = 200, seed = 113)
#>  term                value 
#>  n_star              425   
#>  mean_N              428   
#>  median_N            430   
#>  sd_N                31    
#>  ratio_mean_N_n_star 1.01  
#>  coverage            0.885 
#>  se_coverage         0.0226
#> 
#> Confidence level: 90%
# }
```
