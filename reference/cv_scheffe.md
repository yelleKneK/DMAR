# Provides the Critical Value for the Scheffé Procedure

Provides the Critical Value for the Scheffé Procedure

## Usage

``` r
cv_scheffe(alpha_level, df_numerator, df_denominator, verbose = TRUE)
```

## Arguments

- alpha_level:

  Type I error rate (i.e., the family-wise false-positive rate).

- df_numerator:

  The numerator degrees of freedom (typically the number of groups minus
  1, \\k - 1\\, in a one-way ANOVA).

- df_denominator:

  The denominator (error) degrees of freedom (typically \\N - k\\ in a
  one-way ANOVA).

- verbose:

  Provides extra information about areas under the curve.

## Value

Returns the critical value in a output style (a `data.frame` with one
row per critical value, following the format used by
[`cv_t`](https://yelleknek.github.io/DMAR/reference/cv_t.md) and
[`cv_tukey_hsd`](https://yelleknek.github.io/DMAR/reference/cv_tukey_hsd.md)).

## Details

The Scheffé critical value protects the family-wise error rate for the
simultaneous test of *any* contrast (or family of contrasts) in a
fixed-effects ANOVA, including data-driven contrasts selected after
looking at the data. It is therefore the most conservative of the
standard procedures.

The critical value, on the scale of a *t*-statistic, is
\$\$t\_{\mathrm{crit}}^{\mathrm{Scheffe}} = \sqrt{(k-1)\\
F\_{1-\alpha,\\k-1,\\df\_{\mathrm{denominator}}}},\$\$ so that a
contrast \\\hat\psi\\ with standard error \\\mathit{SE}\_{\hat\psi}\\ is
declared significant when \\\|\hat\psi/\mathit{SE}\_{\hat\psi}\| \>
t\_{\mathrm{crit}}^{\mathrm{Scheffe}}\\. The corresponding simultaneous
confidence interval is \\\hat\psi \pm
t\_{\mathrm{crit}}^{\mathrm{Scheffe}} \cdot \mathit{SE}\_{\hat\psi}\\.

Like the Studentized range distribution underlying Tukey HSD, the
Scheffé reference is one-sided (the underlying *F* statistic is
non-negative), so `alpha_level` is *not* split between two tails.

The Scheffé critical value is a function of a univariate *F* quantile,
which base R supplies through
[`qf`](https://rdrr.io/r/stats/Fdist.html), so unlike
[`cv_dunnett`](https://yelleknek.github.io/DMAR/reference/cv_dunnett.md)
and [`cv_smm`](https://yelleknek.github.io/DMAR/reference/cv_smm.md)
this function needs no multivariate distribution machinery. Scheffé's
procedure earns its simultaneous protection over the infinite family of
all possible contrasts by projecting onto the overall *F* test rather
than by integrating a multivariate *t* density; that is why a single
univariate quantile suffices and the mvtnorm package is not required
here. Maxwell, Delaney, and Kelley (2027, Chapter 5) develop the Scheffé
method as the procedure for arbitrary contrasts within the
multiple-comparisons problem.

## References

Scheffe, H. (1953). A method for judging all contrasts in the analysis
of variance. *Biometrika, 40*, 87–104.

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge. (See Chapter 5 on the multiple-comparisons problem,
where the Scheffé method for arbitrary contrasts is developed.)

## See also

[`cv_t`](https://yelleknek.github.io/DMAR/reference/cv_t.md),
[`cv_tukey_hsd`](https://yelleknek.github.io/DMAR/reference/cv_tukey_hsd.md),
[`contrast_test`](https://yelleknek.github.io/DMAR/reference/contrast_test.md)

Other critical values:
[`cv_bonferroni_f()`](https://yelleknek.github.io/DMAR/reference/cv_bonferroni_f.md),
[`cv_bryant_paulson()`](https://yelleknek.github.io/DMAR/reference/cv_bryant_paulson.md),
[`cv_chisq()`](https://yelleknek.github.io/DMAR/reference/cv_chisq.md),
[`cv_dunnett()`](https://yelleknek.github.io/DMAR/reference/cv_dunnett.md),
[`cv_f()`](https://yelleknek.github.io/DMAR/reference/cv_f.md),
[`cv_smm()`](https://yelleknek.github.io/DMAR/reference/cv_smm.md),
[`cv_t()`](https://yelleknek.github.io/DMAR/reference/cv_t.md),
[`cv_tukey_hsd()`](https://yelleknek.github.io/DMAR/reference/cv_tukey_hsd.md),
[`cv_z()`](https://yelleknek.github.io/DMAR/reference/cv_z.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# Following the arbitrary-contrasts setting of Maxwell, Delaney, and Kelley
# (2027, Chapter 5): a one-way ANOVA with k = 4 groups and 36 error degrees
# of freedom (e.g., n = 10 per group). The Scheffé critical value protects
# the family-wise error rate over any contrast, including contrasts chosen
# after looking at the data.
cv_scheffe(alpha_level = .05, df_numerator = 3, df_denominator = 36)
#>  term     value area_less area_greater
#>  upper_cv 2.93  0.95      0.05        

# The price of that protection is a larger multiplier than the unadjusted
# t critical value at the same error degrees of freedom:
cv_t(alpha_level = .05, df = 36)
#>  term     value area_less area_greater
#>  lower_cv -2.03 0.025     0.975       
#>  upper_cv 2.03  0.975     0.025       
```
