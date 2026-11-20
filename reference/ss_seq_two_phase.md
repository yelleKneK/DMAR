# Sequential Sample Size for a Fixed-Width Interval for a Mean Under Two-Phase Sampling

Evaluates, at the data in hand, a sequential stopping rule for a
two-phase design: the proxy measurement is available on every unit of a
pool, the costly target measurement is added to units drawn at random
from the pool a few at a time, and sampling stops at the first number of
labeled units, at or after a minimum of `m0`, at which the confidence
interval from
[`ci_two_phase`](https://yelleknek.github.io/DMAR/reference/ci_two_phase.md)
has a half-width no larger than `half_width`. Unlike
[`ss_aipe_two_phase`](https://yelleknek.github.io/DMAR/reference/ss_aipe_two_phase.md),
the rule needs no planning values for the target's standard deviation or
its correlation with the proxy: the data supply both as they accrue. It
follows
[`ss_seq_c`](https://yelleknek.github.io/DMAR/reference/ss_seq_c.md) in
its stopping rule, in the `half_width` and `alpha_level` arguments, and
in the rows it returns (with the sizes named for the labeled units), and
like it keeps no state: call it again after each batch. It differs from
[`ss_seq_c()`](https://yelleknek.github.io/DMAR/reference/ss_seq_c.md)
in three ways: it takes the data rather than summary statistics; it has
no `quantile` argument, since its criterion is the half-width of the
interval
[`ci_two_phase`](https://yelleknek.github.io/DMAR/reference/ci_two_phase.md)
reports, which uses the *t* quantile (the default of
[`ss_seq_c()`](https://yelleknek.github.io/DMAR/reference/ss_seq_c.md));
and its minimum `m0` counts labeled units and defaults to 30 rather than
[`ss_seq_c()`](https://yelleknek.github.io/DMAR/reference/ss_seq_c.md)'s
10 per group, because the slope as well as the variance must be
estimated before the rule can be trusted.

## Usage

``` r
ss_seq_two_phase(y, x, x_unlabeled, half_width, alpha_level = 0.05, m0 = 30)
```

## Arguments

- y:

  The target measurement on the units labeled so far.

- x:

  The proxy measurement on the same units, aligned with `y`.

- x_unlabeled:

  The proxy measurement on the units of the pool not yet labeled.

- half_width:

  The target half-width \\h\\ of the 100(1 - 2\\\alpha\\)% confidence
  interval, in the units of the target. For a noninferiority or
  equivalence decision at bounds \\\pm\delta\\,
  [`ss_seq_c`](https://yelleknek.github.io/DMAR/reference/ss_seq_c.md)
  explains why \\\delta/2\\ is the usual choice.

- alpha_level:

  One-sided rate per bound; the interval is at confidence level 1 -
  2\\\alpha\\. Default `0.05` (a 90% interval), as in
  [`ss_seq_c`](https://yelleknek.github.io/DMAR/reference/ss_seq_c.md).

- m0:

  The minimum number of labeled units before the rule may stop, default
  `30`. The early variance and slope estimates are unstable, and a rule
  allowed to stop on them stops too early too often.

## Value

A `data.frame` (class `dmar_tbl`) with rows `stop` (1 when the rule is
met, so labeling stops; 0 to continue), `half_width_current` (the
half-width of the interval now), `half_width_target`,
`n_labeled_current` (the labeled units used, after any with a missing
value are set aside, as in
[`ci_two_phase`](https://yelleknek.github.io/DMAR/reference/ci_two_phase.md)),
`n_labeled_projected`, and `n_pool` (\\M\\).

## Details

**The rule.** With \\n\\ labeled units and a pool of \\M\\, sampling
stops when \\n \ge m_0\\ and \$\$t\_{1-\alpha,\\ n-2}
\left\\\frac{S_e^2}{n} + \frac{\hat b^2 S_X^2}{M}\right\\^{1/2} \le
h,\$\$ the left side being the half-width of the regression interval of
[`ci_two_phase`](https://yelleknek.github.io/DMAR/reference/ci_two_phase.md)
at confidence level 1 - 2\\\alpha\\ (\\S_e^2\\ the residual mean square
of the target on the proxy, \\\hat b\\ the slope, \\S_X^2\\ the proxy's
variance over the pool). This carries the fixed-width sequential
approach of Chow and Robbins (1965), which
[`ss_seq_c`](https://yelleknek.github.io/DMAR/reference/ss_seq_c.md)
applies to a contrast, over to the two-phase interval. Labeled units
must be added at random from the units not yet labeled, so that at every
step they remain a simple random subsample of the pool. In an evaluation
of a system whose outputs an automated judge labels, for example, a
human assessor labels outputs drawn at random from those the judge has
labeled until the rule says stop.

**The projection.** `n_labeled_projected` is the number of labeled units
at which the rule would be met if the current estimates held, from the
normal approximation \\n = S_e^2 / \\(h / z\_{1-\alpha})^2 - \hat b^2
S_X^2 / M\\\\, the same device
[`ss_seq_c`](https://yelleknek.github.io/DMAR/reference/ss_seq_c.md)
uses: a guide to how much further to go, not a commitment. When the
current estimates say that labeling every unit of the pool but one (the
interval needs at least one unit with the proxy alone) would not reach
the target, the projection is `NA` with a warning: the pool, not the
labeling, is then the constraint.

**A constant target, or one the proxy reproduces exactly.** While every
labeled unit has the same target value (for example, a binary target
that has been 1 on every unit so far), the residual variation and the
slope are both zero, and the interval would have zero width. While the
target is an exact linear function of the proxy on the labeled units
(for example, a binary judge and a binary human label that have agreed
on every unit so far), the residual variation is zero and the width of
the interval would come from the size of the pool alone. In both cases
the rule could stop on no evidence, so the function stops with an error
asking for more labeled units, whether or not the rule would have been
met, as it does while the proxy is constant on them. The residual
variance is taken to be zero when it is no larger than
`.Machine$double.eps` times the variance of the target on the labeled
units, which allows for rounding error.

**Coverage.** The stopping size is random and the interval at stopping
is computed as though it had been fixed, so its coverage is not
guaranteed in finite samples; it can fall short of nominal when \\n\\ at
stopping is small, or for a bounded target (a proportion, for example)
whose variance depends on its mean. The minimum `m0` is the guard
against the worst of this, and a simulation of the whole procedure in a
plausible population is the check for a given design.

## References

Chow, Y. S., & Robbins, H. (1965). On the asymptotic theory of
fixed-width sequential confidence intervals for the mean. *The Annals of
Mathematical Statistics, 36*(2), 457–462.
[doi:10.1214/aoms/1177700156](https://doi.org/10.1214/aoms/1177700156)

Neyman, J. (1938). Contribution to the theory of sampling human
populations. *Journal of the American Statistical Association, 33*(201),
101–116.
[doi:10.1080/01621459.1938.10503378](https://doi.org/10.1080/01621459.1938.10503378)

## See also

[`ss_seq_c`](https://yelleknek.github.io/DMAR/reference/ss_seq_c.md),
the sequential rule for a contrast;
[`ci_two_phase`](https://yelleknek.github.io/DMAR/reference/ci_two_phase.md);
[`ss_aipe_two_phase`](https://yelleknek.github.io/DMAR/reference/ss_aipe_two_phase.md)
for a design fixed in advance from planning values.

Other sequential estimation:
[`ss_seq_c()`](https://yelleknek.github.io/DMAR/reference/ss_seq_c.md),
[`ss_seq_c_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_seq_c_sensitivity.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# A proxy is available on a pool of 5,000 units; the target has been
# measured on a random 40 of them so far. Is a 90% interval with
# half-width 1 reached?
set.seed(113)
proxy  <- rnorm(5000, mean = 50, sd = 10)
target <- 20 + 0.6 * proxy + rnorm(5000, sd = 5)
order  <- sample(5000)
first  <- order[1:40]
ss_seq_two_phase(y = target[first], x = proxy[first],
                 x_unlabeled = proxy[-first], half_width = 1)
#>  term                value
#>  stop                0    
#>  half_width_current  1.124
#>  half_width_target   1    
#>  n_labeled_current   40   
#>  n_labeled_projected 49   
#>  n_pool              5000 
#> 
#> Confidence level: 90%

# Not yet; the projection suggests about 50. After labeling the next 30
# units at random:
first <- order[1:70]
ss_seq_two_phase(y = target[first], x = proxy[first],
                 x_unlabeled = proxy[-first], half_width = 1)
#>  term                value
#>  stop                1    
#>  half_width_current  0.896
#>  half_width_target   1    
#>  n_labeled_current   70   
#>  n_labeled_projected 70   
#>  n_pool              5000 
#> 
#> Confidence level: 90%
```
