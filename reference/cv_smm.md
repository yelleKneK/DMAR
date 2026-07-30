# Provides the Critical Value of the Studentized Maximum Modulus Distribution

Provides the Critical Value of the Studentized Maximum Modulus
Distribution

## Usage

``` r
cv_smm(alpha_level, df, n_comparisons, verbose = TRUE)
```

## Arguments

- alpha_level:

  Type I error rate (i.e., the family-wise false-positive rate).

- df:

  Error degrees of freedom (typically \\N - k\\ in a one-way ANOVA). May
  be `Inf`, in which case the known-variance (normal) limit is used.

- n_comparisons:

  The number of simultaneous comparisons (\\m\\) for which a family-wise
  critical value is desired.

- verbose:

  Provides extra information about areas under the curve.

## Value

Returns the critical value as a `data.frame`, following the format used
by [`cv_t`](https://yelleknek.github.io/DMAR/reference/cv_t.md).

## Details

The Studentized maximum modulus (SMM) distribution is the distribution
of \$\$\max\_{i=1,\ldots,m}\\ \|Z_i\| \\/\\ S,\$\$ where \\Z_1, \ldots,
Z_m\\ are independent standard normal variates and \\S\\ is independent
of the \\Z_i\\ and equal to \\\sqrt{\chi^2\_{df} / df}\\.

**What the SMM is used for.** The SMM critical value is the multiplier
that turns a set of \\m\\ individual estimates into a family of
simultaneous confidence intervals, or equivalently a family of tests,
that jointly control the family-wise error rate at level `alpha_level`.
Because the modulus is the largest of \\m\\ standardized statistics in
absolute value, requiring that maximum to clear the critical value
bounds the chance of any one of the \\m\\ intervals failing to cover (or
any one of the \\m\\ tests producing a false positive). Maxwell,
Delaney, and Kelley (2027, Chapter 5) describe this use in the context
of the multiple-comparisons problem: when a researcher forms several
means or contrasts and wants the stated coverage to hold across the
whole set rather than one interval at a time, the SMM supplies the
simultaneous critical value. A pair-by-pair construction, applied to
\\m\\ comparisons, is \\\hat\psi_i \pm
c\_{\alpha;m,df}\\\mathit{SE}\_{\hat\psi_i}\\ with \\c\_{\alpha;m,df}\\
the SMM critical value returned here. When \\m = 1\\ it reduces to the
ordinary two-sided *t* critical value.

**How it is computed.** The \\m\\ statistics are independent given the
common scale estimate \\S\\, so the joint distribution factorizes after
conditioning on \\S\\ and the \\m\\-dimensional integral that defines
the maximum modulus collapses to a single one-dimensional integral,
\$\$P\\\left(\max_i \|Z_i\|/S \le c\right) = \int_0^\infty
\bigl\[\\2\\\Phi(c\\s) - 1\\\bigr\]^{m}\\ f_S(s)\\ ds,\$\$ where \\f_S\\
is the density of \\S = \sqrt{\chi^2\_{df}/df}\\. This package evaluates
that integral with [`integrate`](https://rdrr.io/r/stats/integrate.html)
and inverts it with [`uniroot`](https://rdrr.io/r/stats/uniroot.html),
so the returned value is deterministic and accurate to the solver
tolerance (there is no Monte Carlo simulation and no random seed). When
`df = Inf` the scale is degenerate at 1 and the closed form \\c =
\Phi^{-1}\\\bigl((1 + (1-\alpha)^{1/m})/2\bigr)\\ is returned; when \\m
= 1\\ the value is exactly \\t\_{1-\alpha/2,\\df}\\.

## References

Stoline, M. R., & Ury, H. K. (1979). Tables of the Studentized maximum
modulus distribution and an application to multiple comparisons among
means. *Technometrics, 21*(1), 87–93.

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge. (See Chapter 5 on the multiple-comparisons problem,
where simultaneous confidence intervals for several means or contrasts
are developed; Appendix Table A.2 reports SMM critical values.)

## See also

[`cv_t`](https://yelleknek.github.io/DMAR/reference/cv_t.md),
[`cv_tukey_hsd`](https://yelleknek.github.io/DMAR/reference/cv_tukey_hsd.md),
[`cv_dunnett`](https://yelleknek.github.io/DMAR/reference/cv_dunnett.md),
[`cv_scheffe`](https://yelleknek.github.io/DMAR/reference/cv_scheffe.md)

Other critical values:
[`cv_bonferroni_f()`](https://yelleknek.github.io/DMAR/reference/cv_bonferroni_f.md),
[`cv_bryant_paulson()`](https://yelleknek.github.io/DMAR/reference/cv_bryant_paulson.md),
[`cv_chisq()`](https://yelleknek.github.io/DMAR/reference/cv_chisq.md),
[`cv_dunnett()`](https://yelleknek.github.io/DMAR/reference/cv_dunnett.md),
[`cv_f()`](https://yelleknek.github.io/DMAR/reference/cv_f.md),
[`cv_scheffe()`](https://yelleknek.github.io/DMAR/reference/cv_scheffe.md),
[`cv_t()`](https://yelleknek.github.io/DMAR/reference/cv_t.md),
[`cv_tukey_hsd()`](https://yelleknek.github.io/DMAR/reference/cv_tukey_hsd.md),
[`cv_z()`](https://yelleknek.github.io/DMAR/reference/cv_z.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# Following the simultaneous-intervals setting of Maxwell, Delaney, and
# Kelley (2027, Chapter 5): a researcher forms m = 5 comparisons and wants
# all five confidence intervals to hold simultaneously at the .05 level,
# with 36 error degrees of freedom.
cv_smm(alpha_level = .05, df = 36, n_comparisons = 5)
#>  term     value area_less area_greater
#>  upper_cv 2.7   0.95      0.05        

# When m = 1, the SMM critical value reduces to the two-sided
# t critical value:
cv_smm(alpha_level = .05, df = 36, n_comparisons = 1)$value
#> [1] 2.028094
cv_t(alpha_level = .05, df = 36)$value[2]   # upper_cv from cv_t
#> [1] 2.028094
```
