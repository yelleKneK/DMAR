# Cohen's *h* Effect Size for a Difference Between Two Proportions

Computes Cohen's *h*, the effect size for the difference between two
proportions on the arcsine (variance-stabilizing) scale, \$\$h =
\varphi_1 - \varphi_2, \qquad \varphi_i = 2\\\arcsin\\\sqrt{p_i}.\$\$
The arcsine transform spaces proportions so that a given *h* carries the
same detectability wherever the proportions sit, which a raw difference
\\p_1 - p_2\\ does not: a shift from .01 to .05 is easier to detect than
one from .41 to .45, and *h* reflects that while the raw difference does
not. Cohen's *h* is the proportion analogue of the standardized mean
difference ([`smd`](https://yelleknek.github.io/DMAR/reference/smd.md)):
the effect size on which power and sample size planning for a difference
between two proportions is conventionally based.

## Usage

``` r
cohen_h(p1, p2)
```

## Arguments

- p1, p2:

  The two proportions, each in \\\[0, 1\]\\. `h` is \\\varphi(p_1) -
  \varphi(p_2)\\, so it is positive when `p1` is the larger.

## Value

A 1-row `data.frame` with columns `term` and `value`; `term` is
`"cohen_h"` and `value` is the signed effect size.

## Details

*Cohen's h is a population quantity*: supplied with population
proportions it returns the population value, supplied with sample
proportions it returns the corresponding sample value. It is signed,
positive when `p1` exceeds `p2`; its magnitude
[`abs()`](https://rdrr.io/r/base/MathFun.html) is the size of the effect
irrespective of direction.

## References

Cohen, J. (1988). *Statistical power analysis for the behavioral
sciences* (2nd ed.). Hillsdale, NJ: Lawrence Erlbaum. (See Chapter 6.)

## See also

[`smd`](https://yelleknek.github.io/DMAR/reference/smd.md) for the
standardized mean difference,
[`cohen_f`](https://yelleknek.github.io/DMAR/reference/cohen_f.md),
[`ci_proportion`](https://yelleknek.github.io/DMAR/reference/ci_proportion.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# A shift from .40 to .55.
cohen_h(p1 = 0.55, p2 = 0.40)
#>  term    value
#>  cohen_h 0.302

# The same raw difference near the floor is a larger h, since a difference is
# easier to detect where the proportions are small.
cohen_h(p1 = 0.20, p2 = 0.05)
#>  term    value
#>  cohen_h 0.476

# Its magnitude is the size irrespective of direction.
abs(cohen_h(p1 = 0.40, p2 = 0.55)$value)
#> [1] 0.3015253
```
