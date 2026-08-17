# Provides the Critical Value for Dunnett's Many-to-One Comparisons Procedure

Provides the Critical Value for Dunnett's Many-to-One Comparisons
Procedure

## Usage

``` r
cv_dunnett(
  alpha_level,
  df,
  n_comparisons,
  alternative = "not_equal",
  verbose = TRUE
)
```

## Arguments

- alpha_level:

  Type I error rate (i.e., the family-wise false-positive rate).

- df:

  Error degrees of freedom (typically \\N - k\\, where \\k\\ is the
  total number of groups including the control). May be `Inf`, in which
  case the known-variance (normal) limit is used.

- n_comparisons:

  The number of treatment-versus-control comparisons (i.e., \\k - 1\\
  where \\k\\ is the total number of groups).

- alternative:

  The form of the alternative hypothesis: one of `"not_equal"`
  (two-sided; default), `"greater"` (treatments greater than control),
  or `"less"` (treatments less than control).

- verbose:

  Provides extra information about areas under the curve.

## Value

Returns the critical value in a output style (a `data.frame` following
the format used by
[`cv_t`](https://yelleknek.github.io/DMAR/reference/cv_t.md)).

## Details

Dunnett's procedure controls the family-wise error rate for the special
case of comparing each of \\k - 1\\ treatments to a single control (the
many-to-one comparisons setting), using a multivariate *t* reference
with constant pairwise correlation \\1/2\\. That correlation is exact
for a balanced design: each comparison is \\(\bar Y_i - \bar Y_0)\\, all
sharing the one control mean \\\bar Y_0\\ of variance \\\sigma^2/n\\,
while each comparison has variance \\2\sigma^2/n\\, so any two
comparisons correlate \\(\sigma^2/n)/(2\sigma^2/n) = 1/2\\. Maxwell,
Delaney, and Kelley (2027, Chapter 5) develop the many-to-one
comparisons problem and tabulate these critical values.

**How it is computed.** The critical value is a quantile of a
\\(k-1)\\-dimensional multivariate *t* distribution whose correlation
matrix has a unit diagonal and off-diagonal entries of 1/2. Because that
correlation is a single common value, the comparisons have the
one-factor representation \\Z_i = \sqrt{1/2}\\W + \sqrt{1/2}\\U_i\\ with
a shared factor \\W\\ and independent \\U_i\\, all standard normal.
Conditioning on \\W\\ and on the common scale estimate \\S =
\sqrt{\chi^2\_{df}/df}\\ makes the comparisons independent, so the
\\(k-1)\\-dimensional integral collapses to two nested one-dimensional
integrals: \$\$P\\\left(\max_i T_i \le d\right) =
\int_0^\infty\\\\\int\_{-\infty}^{\infty} \Bigl\[\Phi\bigl((d\\s -
\sqrt{1/2}\\w)/\sqrt{1/2}\bigr)\Bigr\]^{k-1} \phi(w)\\dw\\
f_S(s)\\ds\$\$ for the one-sided value (the two-sided value replaces the
bracket with the probability that \\\|T_i\| \le d\\). This package
evaluates those integrals with
[`integrate`](https://rdrr.io/r/stats/integrate.html) and inverts with
[`uniroot`](https://rdrr.io/r/stats/uniroot.html), so the returned value
is deterministic and accurate to the solver tolerance, with no Monte
Carlo simulation and no random seed. When \\k - 1 = 1\\ the value is the
ordinary one- or two-sided *t* critical value.

## Note

The constant-correlation assumption holds for balanced designs (equal n
per group); for severely unbalanced designs use `glht` instead.

## References

Dunnett, C. W. (1955). A multiple comparison procedure for comparing
several treatments with a control. *Journal of the American Statistical
Association, 50*(272), 1096–1121.

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge. (See Chapter 5 on the multiple-comparisons problem;
Dunnett's many-to-one comparisons procedure is developed there, and
Appendix Tables A.6 and A.7 report its critical values for two- and
one-tailed tests.)

## See also

[`cv_t`](https://yelleknek.github.io/DMAR/reference/cv_t.md),
[`cv_tukey_hsd`](https://yelleknek.github.io/DMAR/reference/cv_tukey_hsd.md),
[`cv_smm`](https://yelleknek.github.io/DMAR/reference/cv_smm.md),
[`cv_scheffe`](https://yelleknek.github.io/DMAR/reference/cv_scheffe.md),
[`contrast_test`](https://yelleknek.github.io/DMAR/reference/contrast_test.md)

Other critical values:
[`cv_bonferroni_f()`](https://yelleknek.github.io/DMAR/reference/cv_bonferroni_f.md),
[`cv_bryant_paulson()`](https://yelleknek.github.io/DMAR/reference/cv_bryant_paulson.md),
[`cv_chisq()`](https://yelleknek.github.io/DMAR/reference/cv_chisq.md),
[`cv_f()`](https://yelleknek.github.io/DMAR/reference/cv_f.md),
[`cv_scheffe()`](https://yelleknek.github.io/DMAR/reference/cv_scheffe.md),
[`cv_smm()`](https://yelleknek.github.io/DMAR/reference/cv_smm.md),
[`cv_t()`](https://yelleknek.github.io/DMAR/reference/cv_t.md),
[`cv_tukey_hsd()`](https://yelleknek.github.io/DMAR/reference/cv_tukey_hsd.md),
[`cv_z()`](https://yelleknek.github.io/DMAR/reference/cv_z.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# Following the many-to-one comparisons setting of Maxwell, Delaney, and
# Kelley (2027, Chapter 5): three treatment groups each compared to a single
# control (k = 4 groups, so 3 comparisons) with 36 error degrees of freedom.
# The two-sided critical value controls the family-wise error rate at .05.
cv_dunnett(alpha_level = .05, df = 36, n_comparisons = 3)
#>  term     value area_less area_greater
#>  upper_cv 2.45  0.95      0.05        

# When the treatments are expected only to exceed the control, a one-sided
# ("treatments greater than control") critical value is smaller.
cv_dunnett(alpha_level = .05, df = 36, n_comparisons = 3, alternative = "greater")
#>  term     value area_less area_greater
#>  upper_cv 2.13  0.95      0.05        
```
